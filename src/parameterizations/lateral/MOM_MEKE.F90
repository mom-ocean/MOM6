!> Implements the Mesoscale Eddy Kinetic Energy framework
!! with topographic beta effect included in computing beta in Rhines scale

module MOM_MEKE

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum, uvchksum
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, MOM_mesg
use MOM_file_parser,   only : read_param, get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : vardesc, var_desc
use MOM_restart,       only : MOM_restart_CS, register_restart_field, query_initialized
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : vertvisc_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_MEKE_types,    only : MEKE_type

implicit none ; private

#include <MOM_memory.h>

public step_forward_MEKE, MEKE_init, MEKE_alloc_register_restart, MEKE_end

!> Control structure that contains MEKE parameters and diagnostics handles
type, public :: MEKE_CS ; private
  ! Parameters
  real, dimension(:,:), pointer :: equilibrium_value => NULL() !< The equilbrium value
                        !! of MEKE to be calculated at each time step [L2 T-2 ~> m2 s-2]
  real :: MEKE_FrCoeff  !< Efficiency of conversion of ME into MEKE [nondim]
  real :: MEKE_GMcoeff  !< Efficiency of conversion of PE into MEKE [nondim]
  real :: MEKE_GMECoeff !< Efficiency of conversion of MEKE into ME by GME [nondim]
  real :: MEKE_damping  !< Local depth-independent MEKE dissipation rate [T-1 ~> s-1].
  real :: MEKE_Cd_scale !< The ratio of the bottom eddy velocity to the column mean
                        !! eddy velocity, i.e. sqrt(2*MEKE). This should be less than 1
                        !! to account for the surface intensification of MEKE.
  real :: MEKE_Cb       !< Coefficient in the \f$\gamma_{bot}\f$ expression [nondim]
  real :: MEKE_min_gamma!< Minimum value of gamma_b^2 allowed [nondim]
  real :: MEKE_Ct       !< Coefficient in the \f$\gamma_{bt}\f$ expression [nondim]
  logical :: visc_drag  !< If true use the vertvisc_type to calculate bottom drag.
  logical :: MEKE_GEOMETRIC !< If true, uses the GM coefficient formulation from the GEOMETRIC
                        !! framework (Marshall et al., 2012)
  real    :: MEKE_GEOMETRIC_alpha !< The nondimensional coefficient governing the efficiency of the
                        !! GEOMETRIC thickness diffusion.
  logical :: MEKE_equilibrium_alt !< If true, use an alternative calculation for the
                        !! equilibrium value of MEKE.
  logical :: MEKE_equilibrium_restoring !< If true, restore MEKE back to its equilibrium value,
                        !!  which is calculated at each time step.
  logical :: GM_src_alt !< If true, use the GM energy conversion form S^2*N^2*kappa rather
                        !! than the streamfunction for the MEKE GM source term.
  logical :: Rd_as_max_scale !< If true the length scale can not exceed the
                        !! first baroclinic deformation radius.
  logical :: use_old_lscale !< Use the old formula for mixing length scale.
  logical :: use_min_lscale !< Use simple minimum for mixing length scale.
  real :: cdrag         !< The bottom drag coefficient for MEKE [nondim].
  real :: MEKE_BGsrc    !< Background energy source for MEKE [L2 T-3 ~> W kg-1] (= m2 s-3).
  real :: MEKE_dtScale  !< Scale factor to accelerate time-stepping [nondim]
  real :: MEKE_KhCoeff  !< Scaling factor to convert MEKE into Kh [nondim]
  real :: MEKE_Uscale   !< MEKE velocity scale for bottom drag [L T-1 ~> m s-1]
  real :: MEKE_KH       !< Background lateral diffusion of MEKE [L2 T-1 ~> m2 s-1]
  real :: MEKE_K4       !< Background bi-harmonic diffusivity (of MEKE) [L4 T-1 ~> m4 s-1]
  real :: KhMEKE_Fac    !< A factor relating MEKE%Kh to the diffusivity used for
                        !! MEKE itself [nondim].
  real :: viscosity_coeff_Ku !< The scaling coefficient in the expression for
                        !! viscosity used to parameterize lateral harmonic momentum mixing
                        !! by unresolved eddies represented by MEKE.
  real :: viscosity_coeff_Au !< The scaling coefficient in the expression for
                        !! viscosity used to parameterize lateral biharmonic momentum mixing
                        !! by unresolved eddies represented by MEKE.
  real :: Lfixed        !< Fixed mixing length scale [L ~> m].
  real :: aDeform       !< Weighting towards deformation scale of mixing length [nondim]
  real :: aRhines       !< Weighting towards Rhines scale of mixing length [nondim]
  real :: aFrict        !< Weighting towards frictional arrest scale of mixing length [nondim]
  real :: aEady         !< Weighting towards Eady scale of mixing length [nondim]
  real :: aGrid         !< Weighting towards grid scale of mixing length [nondim]
  real :: MEKE_advection_factor !< A scaling in front of the advection of MEKE [nondim]
  real :: MEKE_topographic_beta !< Weight for how much topographic beta is considered
                                !! when computing beta in Rhines scale [nondim]
  real :: MEKE_restoring_rate !< Inverse of the timescale used to nudge MEKE toward its equilibrium value [s-1].

  logical :: kh_flux_enabled !< If true, lateral diffusive MEKE flux is enabled.
  logical :: initialize !< If True, invokes a steady state solver to calculate MEKE.
  logical :: debug      !< If true, write out checksums of data for debugging

  type(diag_ctrl), pointer :: diag => NULL() !< A type that regulates diagnostics output
  !>@{ Diagnostic handles
  integer :: id_MEKE = -1, id_Ue = -1, id_Kh = -1, id_src = -1
  integer :: id_Ub = -1, id_Ut = -1
  integer :: id_GM_src = -1, id_mom_src = -1, id_GME_snk = -1, id_decay = -1
  integer :: id_KhMEKE_u = -1, id_KhMEKE_v = -1, id_Ku = -1, id_Au = -1
  integer :: id_Le = -1, id_gamma_b = -1, id_gamma_t = -1
  integer :: id_Lrhines = -1, id_Leady = -1
  integer :: id_MEKE_equilibrium = -1
  !>@}

  ! Infrastructure
  integer :: id_clock_pass !< Clock for group pass calls
  type(group_pass_type) :: pass_MEKE !< Group halo pass handle for MEKE%MEKE and maybe MEKE%Kh_diff
  type(group_pass_type) :: pass_Kh   !< Group halo pass handle for MEKE%Kh, MEKE%Ku, and/or MEKE%Au
end type MEKE_CS

contains

!> Integrates forward-in-time the MEKE eddy energy equation.
!! See \ref section_MEKE_equations.
subroutine step_forward_MEKE(MEKE, h, SN_u, SN_v, visc, dt, G, GV, US, CS, hu, hv)
  type(MEKE_type),                          pointer       :: MEKE !< MEKE data.
  type(ocean_grid_type),                    intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),                  intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),                    intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h    !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)),        intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)),        intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  type(vertvisc_type),                      intent(in)    :: visc !< The vertical viscosity type.
  real,                                     intent(in)    :: dt   !< Model(baroclinic) time-step [T ~> s].
  type(MEKE_CS),                            pointer       :: CS   !< MEKE control structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)   :: hu   !< Accumlated zonal mass flux [H L2 ~> m3 or kg].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)   :: hv   !< Accumlated meridional mass flux [H L2 ~> m3 or kg]

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    mass, &         ! The total mass of the water column [R Z ~> kg m-2].
    I_mass, &       ! The inverse of mass [R-1 Z-1 ~> m2 kg-1].
    src, &          ! The sum of all MEKE sources [L2 T-3 ~> W kg-1] (= m2 s-3).
    MEKE_decay, &   ! A diagnostic of the MEKE decay timescale [T-1 ~> s-1].
    drag_rate_visc, & ! Near-bottom velocity contribution to bottom dratg [L T-1 ~> m s-1]
    drag_rate, &    ! The MEKE spindown timescale due to bottom drag [T-1 ~> s-1].
    drag_rate_J15, &  ! The MEKE spindown timescale due to bottom drag with the Jansen 2015 scheme.
                    ! Unfortunately, as written the units seem inconsistent. [T-1 ~> s-1].
    del2MEKE, &     ! Laplacian of MEKE, used for bi-harmonic diffusion [T-2 ~> s-2].
    del4MEKE, &     ! Time-integrated MEKE tendency arising from the biharmonic of MEKE [L2 T-2 ~> m2 s-2].
    LmixScale, &    ! Eddy mixing length [L ~> m].
    barotrFac2, &   ! Ratio of EKE_barotropic / EKE [nondim]
    bottomFac2, &   ! Ratio of EKE_bottom / EKE [nondim]
    tmp             ! Temporary variable for diagnostic computation

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    MEKE_uflux, &   ! The zonal advective and diffusive flux of MEKE with units of [R Z L4 T-3 ~> kg m-2 s-3].
                    ! In one place, MEKE_uflux is used as temporary work space with units of [L2 T-2 ~> m2 s-2].
    Kh_u, &         ! The zonal diffusivity that is actually used [L2 T-1 ~> m2 s-1].
    baroHu, &       ! Depth integrated accumulated zonal mass flux [R Z L2 ~> kg].
    drag_vel_u      ! A (vertical) viscosity associated with bottom drag at
                    ! u-points [Z T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    MEKE_vflux, &   ! The meridional advective and diffusive flux of MEKE with units of [R Z L4 T-3 ~> kg m-2 s-3].
                    ! In one place, MEKE_vflux is used as temporary work space with units of [L2 T-2 ~> m2 s-2].
    Kh_v, &         ! The meridional diffusivity that is actually used [L2 T-1 ~> m2 s-1].
    baroHv, &       ! Depth integrated accumulated meridional mass flux [R Z L2 ~> kg].
    drag_vel_v      ! A (vertical) viscosity associated with bottom drag at
                    ! v-points [Z T-1 ~> m s-1].
  real :: Kh_here   ! The local horizontal viscosity [L2 T-1 ~> m2 s-1]
  real :: Inv_Kh_max ! The inverse of the local horizontal viscosity [T L-2 ~> s m-2]
  real :: K4_here   ! The local horizontal biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: Inv_K4_max ! The inverse of the local horizontal biharmonic viscosity [T L-4 ~> s m-4]
  real :: cdrag2
  real :: advFac    ! The product of the advection scaling factor and 1/dt [T-1 ~> s-1]
  real :: mass_neglect ! A negligible mass [R Z ~> kg m-2].
  real :: ldamping  ! The MEKE damping rate [T-1 ~> s-1].
  real :: Rho0      ! A density used to convert mass to distance [R ~> kg m-3].
  real :: sdt       ! dt to use locally [T ~> s] (could be scaled to accelerate)
  real :: sdt_damp  ! dt for damping [T ~> s] (sdt could be split).
  logical :: use_drag_rate ! Flag to indicate drag_rate is finite
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "MOM_MEKE: Module must be initialized before it is used.")
  if (.not.associated(MEKE)) call MOM_error(FATAL, &
         "MOM_MEKE: MEKE must be initialized before it is used.")

  if ((CS%MEKE_damping > 0.0) .or. (CS%MEKE_Cd_scale > 0.0) .or. (CS%MEKE_Cb>0.) &
      .or. CS%visc_drag) then
    use_drag_rate = .true.
  else
    use_drag_rate = .false.
  endif

  ! Only integrate the MEKE equations if MEKE is required.
  if (.not.associated(MEKE%MEKE)) then
!   call MOM_error(FATAL, "MOM_MEKE: MEKE%MEKE is not associated!")
    return
  endif

    if (CS%debug) then
      if (associated(MEKE%mom_src)) &
        call hchksum(MEKE%mom_src, 'MEKE mom_src', G%HI, scale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      if (associated(MEKE%GME_snk)) &
        call hchksum(MEKE%GME_snk, 'MEKE GME_snk', G%HI, scale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      if (associated(MEKE%GM_src)) &
        call hchksum(MEKE%GM_src, 'MEKE GM_src', G%HI, scale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      if (associated(MEKE%MEKE)) call hchksum(MEKE%MEKE, 'MEKE MEKE', G%HI, scale=US%L_T_to_m_s**2)
      call uvchksum("MEKE SN_[uv]", SN_u, SN_v, G%HI, scale=US%s_to_T, &
                    scalar_pair=.true.)
      call uvchksum("MEKE h[uv]", hu, hv, G%HI, haloshift=1, &
                    scale=GV%H_to_m*(US%L_to_m**2))
    endif

    sdt = dt*CS%MEKE_dtScale ! Scaled dt to use for time-stepping
    Rho0 = GV%Rho0
    mass_neglect = GV%H_to_RZ * GV%H_subroundoff
    cdrag2 = CS%cdrag**2

    ! With a depth-dependent (and possibly strong) damping, it seems
    ! advisable to use Strang splitting between the damping and diffusion.
    sdt_damp = sdt ; if (CS%MEKE_KH >= 0.0 .or. CS%MEKE_K4 >= 0.) sdt_damp = 0.5*sdt

    ! Calculate depth integrated mass exchange if doing advection [R Z L2 ~> kg]
    if (CS%MEKE_advection_factor>0.) then
      do j=js,je ; do I=is-1,ie
        baroHu(I,j) = 0.
      enddo ; enddo
      do k=1,nz
        do j=js,je ; do I=is-1,ie
          baroHu(I,j) = hu(I,j,k) * GV%H_to_RZ
        enddo ; enddo
      enddo
      do J=js-1,je ; do i=is,ie
        baroHv(i,J) = 0.
      enddo ; enddo
      do k=1,nz
        do J=js-1,je ; do i=is,ie
          baroHv(i,J) = hv(i,J,k) * GV%H_to_RZ
        enddo ; enddo
      enddo
    endif

    if (CS%MEKE_Cd_scale == 0.0 .and. .not. CS%visc_drag) then
      !$OMP parallel do default(shared) private(ldamping)
      do j=js,je ; do i=is,ie
        drag_rate(i,j) = 0. ; drag_rate_J15(i,j) = 0.
      enddo ; enddo
    endif

    ! Calculate drag_rate_visc(i,j) which accounts for the model bottom mean flow
    if (CS%visc_drag) then
      !$OMP parallel do default(shared)
      do j=js,je ; do I=is-1,ie
        drag_vel_u(I,j) = 0.0
        if ((G%mask2dCu(I,j) > 0.0) .and. (visc%bbl_thick_u(I,j) > 0.0)) &
          drag_vel_u(I,j) = visc%Kv_bbl_u(I,j) / visc%bbl_thick_u(I,j)
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-1,je ; do i=is,ie
        drag_vel_v(i,J) = 0.0
        if ((G%mask2dCv(i,J) > 0.0) .and. (visc%bbl_thick_v(i,J) > 0.0)) &
          drag_vel_v(i,J) = visc%Kv_bbl_v(i,J) / visc%bbl_thick_v(i,J)
      enddo ; enddo

      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        drag_rate_visc(i,j) = (0.25*G%IareaT(i,j) * US%Z_to_L * &
                ((G%areaCu(I-1,j)*drag_vel_u(I-1,j) + &
                  G%areaCu(I,j)*drag_vel_u(I,j)) + &
                 (G%areaCv(i,J-1)*drag_vel_v(i,J-1) + &
                  G%areaCv(i,J)*drag_vel_v(i,J)) ) )
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        drag_rate_visc(i,j) = 0.
      enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do j=js-1,je+1
      do i=is-1,ie+1 ; mass(i,j) = 0.0 ; enddo
      do k=1,nz ; do i=is-1,ie+1
        mass(i,j) = mass(i,j) + G%mask2dT(i,j) * (GV%H_to_RZ * h(i,j,k)) ! [R Z ~> kg m-2]
      enddo ; enddo
      do i=is-1,ie+1
        I_mass(i,j) = 0.0
        if (mass(i,j) > 0.0) I_mass(i,j) = 1.0 / mass(i,j) ! [R-1 Z-1 ~> m2 kg-1]
      enddo
    enddo

    if (CS%initialize) then
      call MEKE_equilibrium(CS, MEKE, G, GV, US, SN_u, SN_v, drag_rate_visc, I_mass)
      CS%initialize = .false.
    endif

    ! Calculates bottomFac2, barotrFac2 and LmixScale
    call MEKE_lengthScales(CS, MEKE, G, GV, US, SN_u, SN_v, MEKE%MEKE, bottomFac2, barotrFac2, LmixScale)
    if (CS%debug) then
      if (CS%visc_drag) &
        call uvchksum("MEKE drag_vel_[uv]", drag_vel_u, drag_vel_v, G%HI, &
                      scale=US%Z_to_m*US%s_to_T, scalar_pair=.true.)
      call hchksum(mass, 'MEKE mass',G%HI,haloshift=1, scale=US%RZ_to_kg_m2)
      call hchksum(drag_rate_visc, 'MEKE drag_rate_visc', G%HI, scale=US%L_T_to_m_s)
      call hchksum(bottomFac2, 'MEKE bottomFac2', G%HI)
      call hchksum(barotrFac2, 'MEKE barotrFac2', G%HI)
      call hchksum(LmixScale, 'MEKE LmixScale', G%HI,scale=US%L_to_m)
    endif

    ! Aggregate sources of MEKE (background, frictional and GM)
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      src(i,j) = CS%MEKE_BGsrc
    enddo ; enddo

    if (associated(MEKE%mom_src)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_FrCoeff*I_mass(i,j)*MEKE%mom_src(i,j)
      enddo ; enddo
    endif

    if (associated(MEKE%GME_snk)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_GMECoeff*I_mass(i,j)*MEKE%GME_snk(i,j)
      enddo ; enddo
    endif

    if (associated(MEKE%GM_src)) then
      if (CS%GM_src_alt) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src(i,j) = src(i,j) - CS%MEKE_GMcoeff*MEKE%GM_src(i,j) / &
                     (GV%Rho0 * MAX(1.0*US%m_to_Z, G%bathyT(i,j)))
        enddo ; enddo
      else
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src(i,j) = src(i,j) - CS%MEKE_GMcoeff*I_mass(i,j)*MEKE%GM_src(i,j)
        enddo ; enddo
      endif
    endif

    if (CS%MEKE_equilibrium_restoring) then
      call MEKE_equilibrium_restoring(CS, G, US, SN_u, SN_v)
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_restoring_rate*(MEKE%MEKE(i,j) - CS%equilibrium_value(i,j))
      enddo ; enddo
    endif

    ! Increase EKE by a full time-steps worth of source
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      MEKE%MEKE(i,j) = (MEKE%MEKE(i,j) + sdt*src(i,j))*G%mask2dT(i,j)
    enddo ; enddo

    if (use_drag_rate) then
      ! Calculate a viscous drag rate (includes BBL contributions from mean flow and eddies)
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        drag_rate(i,j) = (US%L_to_Z*Rho0 * I_mass(i,j)) * sqrt( drag_rate_visc(i,j)**2 + &
                 cdrag2 * ( max(0.0, 2.0*bottomFac2(i,j)*MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2 ) )
      enddo ; enddo
    endif

    ! First stage of Strang splitting
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      ldamping = CS%MEKE_damping + drag_rate(i,j) * bottomFac2(i,j)
      if (MEKE%MEKE(i,j) < 0.) ldamping = 0.
      ! notice that the above line ensures a damping only if MEKE is positive,
      ! while leaving MEKE unchanged if it is negative
      MEKE%MEKE(i,j) =  MEKE%MEKE(i,j) / (1.0 + sdt_damp*ldamping)
      MEKE_decay(i,j) = ldamping*G%mask2dT(i,j)
    enddo ; enddo

    if (CS%kh_flux_enabled .or. CS%MEKE_K4 >= 0.0) then
      ! Update MEKE in the halos for lateral or bi-harmonic diffusion
      call cpu_clock_begin(CS%id_clock_pass)
      call do_group_pass(CS%pass_MEKE, G%Domain)
      call cpu_clock_end(CS%id_clock_pass)
    endif

    if (CS%MEKE_K4 >= 0.0) then
      ! Calculate Laplacian of MEKE using MEKE_uflux and MEKE_vflux as temporary work space.
      !$OMP parallel do default(shared)
      do j=js-1,je+1 ; do I=is-2,ie+1
        ! MEKE_uflux is used here as workspace with units of [L2 T-2 ~> m2 s-2].
        MEKE_uflux(I,j) = ((G%dy_Cu(I,j)*G%IdxCu(I,j)) * G%mask2dCu(I,j)) * &
            (MEKE%MEKE(i+1,j) - MEKE%MEKE(i,j))
      ! This would have units of [R Z L2 T-2 ~> kg s-2]
      ! MEKE_uflux(I,j) = ((G%dy_Cu(I,j)*G%IdxCu(I,j)) * &
      !     ((2.0*mass(i,j)*mass(i+1,j)) / ((mass(i,j)+mass(i+1,j)) + mass_neglect)) ) * &
      !     (MEKE%MEKE(i+1,j) - MEKE%MEKE(i,j))
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-2,je+1 ; do i=is-1,ie+1
        ! MEKE_vflux is used here as workspace with units of [L2 T-2 ~> m2 s-2].
        MEKE_vflux(i,J) = ((G%dx_Cv(i,J)*G%IdyCv(i,J)) * G%mask2dCv(i,J)) * &
            (MEKE%MEKE(i,j+1) - MEKE%MEKE(i,j))
      ! This would have units of [R Z L2 T-2 ~> kg s-2]
      ! MEKE_vflux(i,J) = ((G%dx_Cv(i,J)*G%IdyCv(i,J)) * &
      !     ((2.0*mass(i,j)*mass(i,j+1)) / ((mass(i,j)+mass(i,j+1)) + mass_neglect)) ) * &
      !     (MEKE%MEKE(i,j+1) - MEKE%MEKE(i,j))
      enddo ; enddo

      !$OMP parallel do default(shared)
      do j=js-1,je+1 ; do i=is-1,ie+1 ! del2MEKE has units [T-2 ~> s-2].
        del2MEKE(i,j) = G%IareaT(i,j) * &
            ((MEKE_uflux(I,j) - MEKE_uflux(I-1,j)) + (MEKE_vflux(i,J) - MEKE_vflux(i,J-1)))
      enddo ; enddo

      ! Bi-harmonic diffusion of MEKE
      !$OMP parallel do default(shared) private(K4_here,Inv_K4_max)
      do j=js,je ; do I=is-1,ie
        K4_here = CS%MEKE_K4 ! [L4 T-1 ~> m4 s-1]
        ! Limit Kh to avoid CFL violations.
        Inv_K4_max = 64.0 * sdt * ((G%dy_Cu(I,j)*G%IdxCu(I,j)) * &
                     max(G%IareaT(i,j), G%IareaT(i+1,j)))**2
        if (K4_here*Inv_K4_max > 0.3) K4_here = 0.3 / Inv_K4_max

        ! Here the units of MEKE_uflux are [R Z L4 T-3 ~> kg m2 s-3].
        MEKE_uflux(I,j) = ((K4_here * (G%dy_Cu(I,j)*G%IdxCu(I,j))) * &
            ((2.0*mass(i,j)*mass(i+1,j)) / ((mass(i,j)+mass(i+1,j)) + mass_neglect)) ) * &
            (del2MEKE(i+1,j) - del2MEKE(i,j))
      enddo ; enddo
      !$OMP parallel do default(shared) private(K4_here,Inv_K4_max)
      do J=js-1,je ; do i=is,ie
        K4_here = CS%MEKE_K4 ! [L4 T-1 ~> m4 s-1]
        Inv_K4_max = 64.0 * sdt * ((G%dx_Cv(i,J)*G%IdyCv(i,J)) * max(G%IareaT(i,j), G%IareaT(i,j+1)))**2
        if (K4_here*Inv_K4_max > 0.3) K4_here = 0.3 / Inv_K4_max

        ! Here the units of MEKE_vflux are [R Z L4 T-3 ~> kg m2 s-3].
        MEKE_vflux(i,J) = ((K4_here * (G%dx_Cv(i,J)*G%IdyCv(i,J))) * &
            ((2.0*mass(i,j)*mass(i,j+1)) / ((mass(i,j)+mass(i,j+1)) + mass_neglect)) ) * &
            (del2MEKE(i,j+1) - del2MEKE(i,j))
      enddo ; enddo
      ! Store change in MEKE arising from the bi-harmonic in del4MEKE [L2 T-2 ~> m2 s-2].
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        del4MEKE(i,j) = (sdt*(G%IareaT(i,j)*I_mass(i,j))) * &
            ((MEKE_uflux(I-1,j) - MEKE_uflux(I,j)) + &
             (MEKE_vflux(i,J-1) - MEKE_vflux(i,J)))
      enddo ; enddo
    endif !

    if (CS%kh_flux_enabled) then
      ! Lateral diffusion of MEKE
      Kh_here = max(0., CS%MEKE_Kh)
      !$OMP parallel do default(shared) firstprivate(Kh_here) private(Inv_Kh_max)
      do j=js,je ; do I=is-1,ie
        ! Limit Kh to avoid CFL violations.
        if (associated(MEKE%Kh)) &
          Kh_here = max(0., CS%MEKE_Kh) + &
              CS%KhMEKE_Fac*0.5*(MEKE%Kh(i,j)+MEKE%Kh(i+1,j))
        if (associated(MEKE%Kh_diff)) &
          Kh_here = max(0.,CS%MEKE_Kh) + &
              CS%KhMEKE_Fac*0.5*(MEKE%Kh_diff(i,j)+MEKE%Kh_diff(i+1,j))
        Inv_Kh_max = 2.0*sdt * ((G%dy_Cu(I,j)*G%IdxCu(I,j)) * &
                     max(G%IareaT(i,j),G%IareaT(i+1,j)))
        if (Kh_here*Inv_Kh_max > 0.25) Kh_here = 0.25 / Inv_Kh_max
        Kh_u(I,j) = Kh_here

        ! Here the units of MEKE_uflux and MEKE_vflux are [R Z L4 T-3 ~> kg m2 s-3].
        MEKE_uflux(I,j) = ((Kh_here * (G%dy_Cu(I,j)*G%IdxCu(I,j))) * &
            ((2.0*mass(i,j)*mass(i+1,j)) / ((mass(i,j)+mass(i+1,j)) + mass_neglect)) ) * &
            (MEKE%MEKE(i,j) - MEKE%MEKE(i+1,j))
      enddo ; enddo
      !$OMP parallel do default(shared) firstprivate(Kh_here) private(Inv_Kh_max)
      do J=js-1,je ; do i=is,ie
        if (associated(MEKE%Kh)) &
          Kh_here = max(0.,CS%MEKE_Kh) + CS%KhMEKE_Fac * 0.5*(MEKE%Kh(i,j)+MEKE%Kh(i,j+1))
        if (associated(MEKE%Kh_diff)) &
          Kh_here = max(0.,CS%MEKE_Kh) + CS%KhMEKE_Fac * 0.5*(MEKE%Kh_diff(i,j)+MEKE%Kh_diff(i,j+1))
        Inv_Kh_max = 2.0*sdt * ((G%dx_Cv(i,J)*G%IdyCv(i,J)) * max(G%IareaT(i,j),G%IareaT(i,j+1)))
        if (Kh_here*Inv_Kh_max > 0.25) Kh_here = 0.25 / Inv_Kh_max
        Kh_v(i,J) = Kh_here

        ! Here the units of MEKE_uflux and MEKE_vflux are [R Z L4 T-3 ~> kg m2 s-3].
        MEKE_vflux(i,J) = ((Kh_here * (G%dx_Cv(i,J)*G%IdyCv(i,J))) * &
            ((2.0*mass(i,j)*mass(i,j+1)) / ((mass(i,j)+mass(i,j+1)) + mass_neglect)) ) * &
            (MEKE%MEKE(i,j) - MEKE%MEKE(i,j+1))
      enddo ; enddo
      if (CS%MEKE_advection_factor>0.) then
        advFac = CS%MEKE_advection_factor / sdt ! [T-1 ~> s-1]
        !$OMP parallel do default(shared)
        do j=js,je ; do I=is-1,ie
          ! Here the units of the quantities added to MEKE_uflux are [R Z L4 T-3 ~> kg m2 s-3].
          if (baroHu(I,j)>0.) then
            MEKE_uflux(I,j) = MEKE_uflux(I,j) + baroHu(I,j)*MEKE%MEKE(i,j)*advFac
          elseif (baroHu(I,j)<0.) then
            MEKE_uflux(I,j) = MEKE_uflux(I,j) + baroHu(I,j)*MEKE%MEKE(i+1,j)*advFac
          endif
        enddo ; enddo
        !$OMP parallel do default(shared)
        do J=js-1,je ; do i=is,ie
          ! Here the units of the quantities added to MEKE_vflux are [R Z L4 T-3 ~> kg m2 s-3].
          if (baroHv(i,J)>0.) then
            MEKE_vflux(i,J) = MEKE_vflux(i,J) + baroHv(i,J)*MEKE%MEKE(i,j)*advFac
          elseif (baroHv(i,J)<0.) then
            MEKE_vflux(i,J) = MEKE_vflux(i,J) + baroHv(i,J)*MEKE%MEKE(i,j+1)*advFac
          endif
        enddo ; enddo
      endif

      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        MEKE%MEKE(i,j) = MEKE%MEKE(i,j) + (sdt*(G%IareaT(i,j)*I_mass(i,j))) * &
            ((MEKE_uflux(I-1,j) - MEKE_uflux(I,j)) + &
             (MEKE_vflux(i,J-1) - MEKE_vflux(i,J)))
      enddo ; enddo
    endif ! MEKE_KH>0

    ! Add on bi-harmonic tendency
    if (CS%MEKE_K4 >= 0.0) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        MEKE%MEKE(i,j) = MEKE%MEKE(i,j) + del4MEKE(i,j)
      enddo ; enddo
    endif

    ! Second stage of Strang splitting
    if (CS%MEKE_KH >= 0.0 .or. CS%MEKE_K4 >= 0.0) then
      if (sdt>sdt_damp) then
        ! Recalculate the drag rate, since MEKE has changed.
        if (use_drag_rate) then
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            drag_rate(i,j) = (US%L_to_Z*Rho0 * I_mass(i,j)) * sqrt( drag_rate_visc(i,j)**2 + &
                   cdrag2 * ( max(0.0, 2.0*bottomFac2(i,j)*MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2 ) )
          enddo ; enddo
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            ldamping = CS%MEKE_damping + drag_rate(i,j) * bottomFac2(i,j)
            if (MEKE%MEKE(i,j) < 0.) ldamping = 0.
            ! notice that the above line ensures a damping only if MEKE is positive,
            ! while leaving MEKE unchanged if it is negative
            MEKE%MEKE(i,j) =  MEKE%MEKE(i,j) / (1.0 + sdt_damp*ldamping)
            MEKE_decay(i,j) = ldamping*G%mask2dT(i,j)
          enddo ; enddo
        endif
      endif
    endif ! MEKE_KH>=0

 !   do j=js,je ; do i=is,ie
 !     MEKE%MEKE(i,j) =  MAX(MEKE%MEKE(i,j),0.0)
 !   enddo ; enddo

    call cpu_clock_begin(CS%id_clock_pass)
    call do_group_pass(CS%pass_MEKE, G%Domain)
    call cpu_clock_end(CS%id_clock_pass)

    ! Calculate diffusivity for main model to use
    if (CS%MEKE_KhCoeff>0.) then
      if (.not.CS%MEKE_GEOMETRIC) then
        if (CS%use_old_lscale) then
          if (CS%Rd_as_max_scale) then
            !$OMP parallel do default(shared)
            do j=js,je ; do i=is,ie
              MEKE%Kh(i,j) = (CS%MEKE_KhCoeff * &
                         sqrt(2.*max(0.,barotrFac2(i,j)*MEKE%MEKE(i,j))*G%areaT(i,j)) ) * &
                         min(MEKE%Rd_dx_h(i,j), 1.0)
            enddo ; enddo
          else
            !$OMP parallel do default(shared)
            do j=js,je ; do i=is,ie
              MEKE%Kh(i,j) = CS%MEKE_KhCoeff * &
                  sqrt(2.*max(0., barotrFac2(i,j)*MEKE%MEKE(i,j))*G%areaT(i,j))
            enddo ; enddo
          endif
        else
          !$OMP parallel do default(shared)
          do j=js,je ; do i=is,ie
            MEKE%Kh(i,j) = CS%MEKE_KhCoeff * &
                sqrt(2.*max(0., barotrFac2(i,j)*MEKE%MEKE(i,j))) * LmixScale(i,j)
          enddo ; enddo
        endif
      endif
    endif

    ! Calculate viscosity for the main model to use
    if (CS%viscosity_coeff_Ku /=0.) then
      do j=js,je ; do i=is,ie
        MEKE%Ku(i,j) = CS%viscosity_coeff_Ku * sqrt(2.*max(0.,MEKE%MEKE(i,j))) * LmixScale(i,j)
      enddo ; enddo
    endif

    if (CS%viscosity_coeff_Au /=0.) then
      do j=js,je ; do i=is,ie
        MEKE%Au(i,j) = CS%viscosity_coeff_Au * sqrt(2.*max(0.,MEKE%MEKE(i,j))) * LmixScale(i,j)**3
      enddo ; enddo
    endif

    if (associated(MEKE%Kh) .or. associated(MEKE%Ku) .or. associated(MEKE%Au)) then
      call cpu_clock_begin(CS%id_clock_pass)
      call do_group_pass(CS%pass_Kh, G%Domain)
      call cpu_clock_end(CS%id_clock_pass)
    endif

    ! Offer fields for averaging.
    if (any([CS%id_Ue, CS%id_Ub, CS%id_Ut] > 0)) &
      tmp(:,:) = 0.
    if (CS%id_MEKE>0) call post_data(CS%id_MEKE, MEKE%MEKE, CS%diag)
    if (CS%id_Ue>0) then
      do j=js,je ; do i=is,ie
        tmp(i,j) = sqrt(max(0., 2. * MEKE%MEKE(i,j)))
      enddo ; enddo
      call post_data(CS%id_Ue, tmp, CS%diag)
    endif
    if (CS%id_Ub>0) then
      do j=js,je ; do i=is,ie
        tmp(i,j) = sqrt(max(0., 2. * MEKE%MEKE(i,j) * bottomFac2(i,j)))
      enddo ; enddo
      call post_data(CS%id_Ub, tmp, CS%diag)
    endif
    if (CS%id_Ut>0) then
      do j=js,je ; do i=is,ie
        tmp(i,j) = sqrt(max(0., 2. * MEKE%MEKE(i,j) * barotrFac2(i,j)))
      enddo ; enddo
      call post_data(CS%id_Ut, tmp, CS%diag)
    endif
    if (CS%id_Kh>0) call post_data(CS%id_Kh, MEKE%Kh, CS%diag)
    if (CS%id_Ku>0) call post_data(CS%id_Ku, MEKE%Ku, CS%diag)
    if (CS%id_Au>0) call post_data(CS%id_Au, MEKE%Au, CS%diag)
    if (CS%id_KhMEKE_u>0) call post_data(CS%id_KhMEKE_u, Kh_u, CS%diag)
    if (CS%id_KhMEKE_v>0) call post_data(CS%id_KhMEKE_v, Kh_v, CS%diag)
    if (CS%id_src>0) call post_data(CS%id_src, src, CS%diag)
    if (CS%id_decay>0) call post_data(CS%id_decay, MEKE_decay, CS%diag)
    if (CS%id_GM_src>0) call post_data(CS%id_GM_src, MEKE%GM_src, CS%diag)
    if (CS%id_mom_src>0) call post_data(CS%id_mom_src, MEKE%mom_src, CS%diag)
    if (CS%id_GME_snk>0) call post_data(CS%id_GME_snk, MEKE%GME_snk, CS%diag)
    if (CS%id_Le>0) call post_data(CS%id_Le, LmixScale, CS%diag)
    if (CS%id_gamma_b>0) then
      do j=js,je ; do i=is,ie
        bottomFac2(i,j) = sqrt(bottomFac2(i,j))
      enddo ; enddo
      call post_data(CS%id_gamma_b, bottomFac2, CS%diag)
    endif
    if (CS%id_gamma_t>0) then
      do j=js,je ; do i=is,ie
        barotrFac2(i,j) = sqrt(barotrFac2(i,j))
      enddo ; enddo
      call post_data(CS%id_gamma_t, barotrFac2, CS%diag)
    endif

end subroutine step_forward_MEKE

!> Calculates the equilibrium solutino where the source depends only on MEKE diffusivity
!! and there is no lateral diffusion of MEKE.
!! Results is in MEKE%MEKE.
subroutine MEKE_equilibrium(CS, MEKE, G, GV, US, SN_u, SN_v, drag_rate_visc, I_mass)
  type(ocean_grid_type),             intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),           intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type
  type(MEKE_CS),                     pointer       :: CS   !< MEKE control structure.
  type(MEKE_type),                   pointer       :: MEKE !< A structure with MEKE data.
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: drag_rate_visc !< Mean flow velocity contribution
                                                           !! to the MEKE drag rate [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: I_mass  !< Inverse of column mass [R-1 Z-1 ~> m2 kg-1].
  ! Local variables
  real :: beta ! Combined topograpic and planetary vorticity gradient [T-1 L-1 ~> s-1 m-1]
  real :: SN   ! The local Eady growth rate [T-1 ~> s-1]
  real :: bottomFac2, barotrFac2    ! Vertical structure factors [nondim]
  real :: LmixScale, LRhines, LEady ! Various mixing length scales [L ~> m]
  real :: I_H, KhCoeff
  real :: Kh    ! A lateral diffusivity [L2 T-1 ~> m2 s-1]
  real :: Ubg2  ! Background (tidal?) velocity squared [L2 T-2 ~> m2 s-2]
  real :: cd2
  real :: drag_rate ! The MEKE spindown timescale due to bottom drag [T-1 ~> s-1].
  real :: src   ! The sum of MEKE sources [L2 T-3 ~> W kg-1]
  real :: ldamping  ! The MEKE damping rate [T-1 ~> s-1].
  real :: EKE, EKEmin, EKEmax, EKEerr ! [L2 T-2 ~> m2 s-2]
  real :: resid, ResMin, ResMax ! Residuals [L2 T-3 ~> W kg-1]
  real :: FatH    ! Coriolis parameter at h points; to compute topographic beta [T-1 ~> s-1]
  real :: beta_topo_x, beta_topo_y    ! Topographic PV gradients in x and y [T-1 L-1 ~> s-1 m-1]
  integer :: i, j, is, ie, js, je, n1, n2
  real :: tolerance ! Width of EKE bracket [L2 T-2 ~> m2 s-2].
  logical :: useSecant, debugIteration

  real :: Lgrid, Ldeform, Lfrict

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  debugIteration = .false.
  KhCoeff = CS%MEKE_KhCoeff
  Ubg2 = CS%MEKE_Uscale**2
  cd2 = CS%cdrag**2
  tolerance = 1.0e-12*US%m_s_to_L_T**2

!$OMP do
  do j=js,je ; do i=is,ie
    ! SN = 0.25*max( (SN_u(I,j) + SN_u(I-1,j)) + (SN_v(i,J) + SN_v(i,J-1)), 0.)
    ! This avoids extremes values in equilibrium solution due to bad values in SN_u, SN_v
    SN = min(SN_u(I,j), SN_u(I-1,j), SN_v(i,J), SN_v(i,J-1))

    if (CS%MEKE_equilibrium_alt) then
      MEKE%MEKE(i,j) = (CS%MEKE_GEOMETRIC_alpha * SN * US%Z_to_m*G%bathyT(i,j))**2 / cd2
    else
      FatH = 0.25*((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1)) + &
                   (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J-1))) ! Coriolis parameter at h points

      ! Since zero-bathymetry cells are masked, this avoids calculations on land
      if (CS%MEKE_topographic_beta == 0. .or. G%bathyT(i,j) == 0.) then
        beta_topo_x = 0. ; beta_topo_y = 0.
      else
        !### Consider different combinations of these estimates of topographic beta, and the use
        !    of the water column thickness instead of the bathymetric depth.
        beta_topo_x = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (G%bathyT(i+1,j)-G%bathyT(i,j)) * G%IdxCu(I,j)  &
                  / max(G%bathyT(i+1,j),G%bathyT(i,j), GV%H_subroundoff) &
              +       (G%bathyT(i,j)-G%bathyT(i-1,j)) * G%IdxCu(I-1,j) &
                  / max(G%bathyT(i,j),G%bathyT(i-1,j), GV%H_subroundoff) )
        beta_topo_y = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (G%bathyT(i,j+1)-G%bathyT(i,j)) * G%IdyCv(i,J)  &
                  / max(G%bathyT(i,j+1),G%bathyT(i,j), GV%H_subroundoff) + &
                      (G%bathyT(i,j)-G%bathyT(i,j-1)) * G%IdyCv(i,J-1) &
                  / max(G%bathyT(i,j),G%bathyT(i,j-1), GV%H_subroundoff) )
      endif
      beta =  sqrt((G%dF_dx(i,j) + beta_topo_x)**2 + &
                   (G%dF_dy(i,j) + beta_topo_y)**2 )

      I_H = US%L_to_Z*GV%Rho0 * I_mass(i,j)

      if (KhCoeff*SN*I_H>0.) then
        ! Solve resid(E) = 0, where resid = Kh(E) * (SN)^2 - damp_rate(E) E
        EKEmin = 0.   ! Use the trivial root as the left bracket
        ResMin = 0.   ! Need to detect direction of left residual
        EKEmax = 0.01*US%m_s_to_L_T**2 ! First guess at right bracket
        useSecant = .false. ! Start using a bisection method

        ! First find right bracket for which resid<0
        resid = 1.0*US%m_to_L**2*US%T_to_s**3 ; n1 = 0
        do while (resid>0.)
          n1 = n1 + 1
          EKE = EKEmax
          call MEKE_lengthScales_0d(CS, US, G%areaT(i,j), beta, G%bathyT(i,j), &
                                    MEKE%Rd_dx_h(i,j), SN, EKE, &
                                    bottomFac2, barotrFac2, LmixScale, LRhines, LEady)
          ! TODO: Should include resolution function in Kh
          Kh = (KhCoeff * sqrt(2.*barotrFac2*EKE) * LmixScale)
          src = Kh * (SN * SN)
          drag_rate = I_H * sqrt(drag_rate_visc(i,j)**2 + cd2 * ( 2.0*bottomFac2*EKE + Ubg2 ) )
          ldamping = CS%MEKE_damping + drag_rate * bottomFac2
          resid = src - ldamping * EKE
          ! if (debugIteration) then
          !   write(0,*) n1, 'EKE=',EKE,'resid=',resid
          !   write(0,*) 'EKEmin=',EKEmin,'ResMin=',ResMin
          !   write(0,*) 'src=',src,'ldamping=',ldamping
          !   write(0,*) 'gamma-b=',bottomFac2,'gamma-t=',barotrFac2
          !   write(0,*) 'drag_visc=',drag_rate_visc(i,j),'Ubg2=',Ubg2
          ! endif
          if (resid>0.) then    ! EKE is to the left of the root
            EKEmin = EKE        ! so we move the left bracket here
            EKEmax = 10. * EKE  ! and guess again for the right bracket
            if (resid<ResMin) useSecant = .true.
            ResMin = resid
            if (EKEmax > 2.e17*US%m_s_to_L_T**2) then
              if (debugIteration) stop 'Something has gone very wrong'
              debugIteration = .true.
              resid = 1. ; n1 = 0
              EKEmin = 0. ; ResMin = 0.
              EKEmax = 0.01*US%m_s_to_L_T**2
              useSecant = .false.
            endif
          endif
        enddo ! while(resid>0.) searching for right bracket
        ResMax = resid

        ! Bisect the bracket
        n2 = 0 ; EKEerr = EKEmax - EKEmin
        do while (EKEerr > tolerance)
          n2 = n2 + 1
          if (useSecant) then
            EKE = EKEmin + (EKEmax - EKEmin) * (ResMin / (ResMin - ResMax))
          else
            EKE = 0.5 * (EKEmin + EKEmax)
          endif
          EKEerr = min( EKE-EKEmin, EKEmax-EKE )
          ! TODO: Should include resolution function in Kh
          Kh = (KhCoeff * sqrt(2.*barotrFac2*EKE) * LmixScale)
          src = Kh * (SN * SN)
          drag_rate = I_H * sqrt( drag_rate_visc(i,j)**2 + cd2 * ( 2.0*bottomFac2*EKE + Ubg2 ) )
          ldamping = CS%MEKE_damping + drag_rate * bottomFac2
          resid = src - ldamping * EKE
          if (useSecant .and. resid>ResMin) useSecant = .false.
          if (resid>0.) then              ! EKE is to the left of the root
            EKEmin = EKE                  ! so we move the left bracket here
            if (resid<ResMin) useSecant = .true.
            ResMin = resid                ! Save this for the secant method
          elseif (resid<0.) then          ! EKE is to the right of the root
            EKEmax = EKE                  ! so we move the right bracket here
            ResMax = resid                ! Save this for the secant method
          else
            exit                          ! resid=0 => EKE is exactly at the root
          endif
          if (n2>200) stop 'Failing to converge?'
        enddo ! while(EKEmax-EKEmin>tolerance)

      else
        EKE = 0.
      endif
      MEKE%MEKE(i,j) = EKE
    endif
  enddo ; enddo

end subroutine MEKE_equilibrium


!< This subroutine calculates a new equilibrium value for MEKE at each time step. This is not copied into
!! MEKE%MEKE; rather, it is used as a restoring term to nudge MEKE%MEKE back to an equilibrium value
subroutine MEKE_equilibrium_restoring(CS, G, US, SN_u, SN_v)
  type(ocean_grid_type),             intent(inout) :: G    !< Ocean grid.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type.
  type(MEKE_CS),                     pointer       :: CS   !< MEKE control structure.
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  ! Local variables
  real :: SN                      ! The local Eady growth rate [T-1 ~> s-1]
  integer :: i, j, is, ie, js, je ! local indices
  real :: cd2                     ! bottom drag

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  cd2 = CS%cdrag**2

  if (.not. associated(CS%equilibrium_value)) allocate(CS%equilibrium_value(SZI_(G),SZJ_(G)))
  CS%equilibrium_value(:,:) = 0.0

!$OMP do
  do j=js,je ; do i=is,ie
    ! SN = 0.25*max( (SN_u(I,j) + SN_u(I-1,j)) + (SN_v(i,J) + SN_v(i,J-1)), 0.)
    ! This avoids extremes values in equilibrium solution due to bad values in SN_u, SN_v
    SN = min(SN_u(I,j), SN_u(I-1,j), SN_v(i,J), SN_v(i,J-1))
    CS%equilibrium_value(i,j) = (CS%MEKE_GEOMETRIC_alpha * SN * US%Z_to_m*G%bathyT(i,j))**2 / cd2
  enddo ; enddo

  if (CS%id_MEKE_equilibrium>0) call post_data(CS%id_MEKE_equilibrium, CS%equilibrium_value, CS%diag)

end subroutine MEKE_equilibrium_restoring

!> Calculates the eddy mixing length scale and \f$\gamma_b\f$ and \f$\gamma_t\f$
!! functions that are ratios of either bottom or barotropic eddy energy to the
!! column eddy energy, respectively.  See \ref section_MEKE_equations.
subroutine MEKE_lengthScales(CS, MEKE, G, GV, US, SN_u, SN_v, &
            EKE, bottomFac2, barotrFac2, LmixScale)
  type(MEKE_CS),                     pointer       :: CS   !< MEKE control structure.
  type(MEKE_type),                   pointer       :: MEKE !< MEKE data.
  type(ocean_grid_type),             intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),           intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: EKE  !< Eddy kinetic energy [L2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: bottomFac2 !< gamma_b^2
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: barotrFac2 !< gamma_t^2
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: LmixScale !< Eddy mixing length [L ~> m].
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: LRhines, LEady  ! Possible mixing length scales [L ~> m]
  real :: beta ! Combined topograpic and planetary vorticity gradient [T-1 L-1 ~> s-1 m-1]
  real :: SN   ! The local Eady growth rate [T-1 ~> s-1]
  real :: FatH ! Coriolis parameter at h points [T-1 ~> s-1]
  real :: beta_topo_x, beta_topo_y  ! Topographic PV gradients in x and y [T-1 L-1 ~> s-1 m-1]
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

!$OMP do
  do j=js,je ; do i=is,ie
    if (.not.CS%use_old_lscale) then
      if (CS%aEady > 0.) then
        SN = 0.25 * ( (SN_u(I,j) + SN_u(I-1,j)) + (SN_v(i,J) + SN_v(i,J-1)) )
      else
        SN = 0.
      endif
      FatH = 0.25* ( ( G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1) ) + &
                     ( G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J-1) ) )  ! Coriolis parameter at h points

      ! If bathyT is zero, then a division by zero FPE will be raised.  In this
      ! case, we apply Adcroft's rule of reciprocals and set the term to zero.
      ! Since zero-bathymetry cells are masked, this should not affect values.
      if (CS%MEKE_topographic_beta == 0. .or. G%bathyT(i,j) == 0.0) then
        beta_topo_x = 0. ; beta_topo_y = 0.
      else
        !### Consider different combinations of these estimates of topographic beta, and the use
        !    of the water column thickness instead of the bathymetric depth.
        beta_topo_x = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (G%bathyT(i+1,j)-G%bathyT(i,j)) * G%IdxCu(I,j)  &
                 / max(G%bathyT(i+1,j),G%bathyT(i,j), GV%H_subroundoff) &
              +       (G%bathyT(i,j)-G%bathyT(i-1,j)) * G%IdxCu(I-1,j) &
                 / max(G%bathyT(i,j),G%bathyT(i-1,j), GV%H_subroundoff) )
        beta_topo_y = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (G%bathyT(i,j+1)-G%bathyT(i,j)) * G%IdyCv(i,J)  &
                 / max(G%bathyT(i,j+1),G%bathyT(i,j), GV%H_subroundoff) + &
                      (G%bathyT(i,j)-G%bathyT(i,j-1)) * G%IdyCv(i,J-1) &
                 / max(G%bathyT(i,j),G%bathyT(i,j-1), GV%H_subroundoff) )
      endif
      beta =  sqrt((G%dF_dx(i,j) + beta_topo_x)**2 + &
                   (G%dF_dy(i,j) + beta_topo_y)**2 )

    else
      beta = 0.
    endif
    ! Returns bottomFac2, barotrFac2 and LmixScale
    call MEKE_lengthScales_0d(CS, US, G%areaT(i,j), beta, G%bathyT(i,j),  &
                              MEKE%Rd_dx_h(i,j), SN, MEKE%MEKE(i,j), &
                              bottomFac2(i,j), barotrFac2(i,j), LmixScale(i,j), &
                              LRhines(i,j), LEady(i,j))
  enddo ; enddo
  if (CS%id_Lrhines>0) call post_data(CS%id_LRhines, LRhines, CS%diag)
  if (CS%id_Leady>0) call post_data(CS%id_LEady, LEady, CS%diag)

end subroutine MEKE_lengthScales

!> Calculates the eddy mixing length scale and \f$\gamma_b\f$ and \f$\gamma_t\f$
!! functions that are ratios of either bottom or barotropic eddy energy to the
!! column eddy energy, respectively.  See \ref section_MEKE_equations.
subroutine MEKE_lengthScales_0d(CS, US, area, beta, depth, Rd_dx, SN, EKE, & ! Z_to_L, &
                                bottomFac2, barotrFac2, LmixScale, Lrhines, Leady)
  type(MEKE_CS), pointer       :: CS         !< MEKE control structure.
  type(unit_scale_type), intent(in) :: US    !< A dimensional unit scaling type
  real,          intent(in)    :: area       !< Grid cell area [L2 ~> m2]
  real,          intent(in)    :: beta       !< Planetary beta = \f$ \nabla f\f$  [T-1 L-1 ~> s-1 m-1]
  real,          intent(in)    :: depth      !< Ocean depth [Z ~> m]
  real,          intent(in)    :: Rd_dx      !< Resolution Ld/dx [nondim].
  real,          intent(in)    :: SN         !< Eady growth rate [T-1 ~> s-1].
  real,          intent(in)    :: EKE        !< Eddy kinetic energy [L2 T-2 ~> m2 s-2].
!  real,          intent(in)    :: Z_to_L     !< A conversion factor from depth units (Z) to
!                                             !! the units for lateral distances (L).
  real,          intent(out)   :: bottomFac2 !< gamma_b^2
  real,          intent(out)   :: barotrFac2 !< gamma_t^2
  real,          intent(out)   :: LmixScale  !< Eddy mixing length [L ~> m].
  real,          intent(out)   :: Lrhines    !< Rhines length scale [L ~> m].
  real,          intent(out)   :: Leady      !< Eady length scale [L ~> m].
  ! Local variables
  real :: Lgrid, Ldeform, Lfrict ! Length scales [L ~> m]
  real :: Ue  ! An eddy velocity [L T-1 ~> m s-1]

  ! Length scale for MEKE derived diffusivity
  Lgrid = sqrt(area)               ! Grid scale
  Ldeform = Lgrid * Rd_dx          ! Deformation scale
  Lfrict = (US%Z_to_L * depth) / CS%cdrag  ! Frictional arrest scale
  ! gamma_b^2 is the ratio of bottom eddy energy to mean column eddy energy
  ! used in calculating bottom drag
  bottomFac2 = CS%MEKE_CD_SCALE**2
  if (Lfrict*CS%MEKE_Cb>0.) bottomFac2 = bottomFac2 + 1./( 1. + CS%MEKE_Cb*(Ldeform/Lfrict) )**0.8
  bottomFac2 = max(bottomFac2, CS%MEKE_min_gamma)
  ! gamma_t^2 is the ratio of barotropic eddy energy to mean column eddy energy
  ! used in the velocity scale for diffusivity
  barotrFac2 = 1.
  if (Lfrict*CS%MEKE_Ct>0.) barotrFac2 = 1. / ( 1. + CS%MEKE_Ct*(Ldeform/Lfrict) )**0.25
  barotrFac2 = max(barotrFac2, CS%MEKE_min_gamma)
  if (CS%use_old_lscale) then
    if (CS%Rd_as_max_scale) then
      LmixScale = min(Ldeform, Lgrid) ! The smaller of Ld or dx
    else
      LmixScale = Lgrid
    endif
  else
    Ue = sqrt( 2.0 * max( 0., barotrFac2*EKE ) ) ! Barotropic eddy flow scale
    Lrhines = sqrt( Ue / max( beta, 1.e-30*US%T_to_s*US%L_to_m ) )       ! Rhines scale
    if (CS%aEady > 0.) then
      Leady = Ue / max( SN, 1.e-15*US%T_to_s ) ! Bound Eady time-scale < 1e15 seconds
    else
      Leady = 0.
    endif
    if (CS%use_min_lscale) then
      LmixScale = 1.e7
      if (CS%aDeform*Ldeform > 0.) LmixScale = min(LmixScale,CS%aDeform*Ldeform)
      if (CS%aFrict *Lfrict  > 0.) LmixScale = min(LmixScale,CS%aFrict *Lfrict)
      if (CS%aRhines*Lrhines > 0.) LmixScale = min(LmixScale,CS%aRhines*Lrhines)
      if (CS%aEady  *Leady   > 0.) LmixScale = min(LmixScale,CS%aEady  *Leady)
      if (CS%aGrid  *Lgrid   > 0.) LmixScale = min(LmixScale,CS%aGrid  *Lgrid)
      if (CS%Lfixed          > 0.) LmixScale = min(LmixScale,CS%Lfixed)
    else
      LmixScale = 0.
      if (CS%aDeform*Ldeform > 0.) LmixScale = LmixScale + 1./(CS%aDeform*Ldeform)
      if (CS%aFrict *Lfrict  > 0.) LmixScale = LmixScale + 1./(CS%aFrict *Lfrict)
      if (CS%aRhines*Lrhines > 0.) LmixScale = LmixScale + 1./(CS%aRhines*Lrhines)
      if (CS%aEady  *Leady   > 0.) LmixScale = LmixScale + 1./(CS%aEady  *Leady)
      if (CS%aGrid  *Lgrid   > 0.) LmixScale = LmixScale + 1./(CS%aGrid  *Lgrid)
      if (CS%Lfixed          > 0.) LmixScale = LmixScale + 1./CS%Lfixed
      if (LmixScale > 0.) LmixScale = 1. / LmixScale
    endif
  endif

end subroutine MEKE_lengthScales_0d

!> Initializes the MOM_MEKE module and reads parameters.
!! Returns True if module is to be used, otherwise returns False.
logical function MEKE_init(Time, G, US, param_file, diag, CS, MEKE, restart_CS)
  type(time_type),         intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),   intent(inout) :: G          !< The ocean's grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file parser structure.
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(MEKE_CS),           pointer       :: CS         !< MEKE control structure.
  type(MEKE_type),         pointer       :: MEKE       !< MEKE-related fields.
  type(MOM_restart_CS),    pointer       :: restart_CS !< Restart control structure for MOM_MEKE.

  ! Local variables
  real    :: I_T_rescale   ! A rescaling factor for time from the internal representation in this
                           ! run to the representation in a restart file.
  real    :: L_rescale     ! A rescaling factor for length from the internal representation in this
                           ! run to the representation in a restart file.
  real    :: MEKE_restoring_timescale ! The timescale used to nudge MEKE toward its equilibrium value.
  real :: cdrag            ! The default bottom drag coefficient [nondim].
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  logical :: laplacian, biharmonic, useVarMix, coldStart
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_MEKE" ! This module's name.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! Determine whether this module will be used
  call get_param(param_file, mdl, "USE_MEKE", MEKE_init, default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, "", all_default=.not.MEKE_init)
  call get_param(param_file, mdl, "USE_MEKE", MEKE_init, &
                 "If true, turns on the MEKE scheme which calculates "// &
                 "a sub-grid mesoscale eddy kinetic energy budget.", &
                 default=.false.)
  if (.not. MEKE_init) return

  if (.not. associated(MEKE)) then
    ! The MEKE structure should have been allocated in MEKE_alloc_register_restart()
    call MOM_error(WARNING, "MEKE_init called with NO associated "// &
                            "MEKE-type structure.")
    return
  endif
  if (associated(CS)) then
    call MOM_error(WARNING, &
      "MEKE_init called with an associated control structure.")
    return
  else ; allocate(CS) ; endif

  call MOM_mesg("MEKE_init: reading parameters ", 5)

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "MEKE_DAMPING", CS%MEKE_damping, &
                 "The local depth-independent MEKE dissipation rate.", &
                 units="s-1", default=0.0, scale=US%T_to_s)
  call get_param(param_file, mdl, "MEKE_CD_SCALE", CS%MEKE_Cd_scale, &
                 "The ratio of the bottom eddy velocity to the column mean "//&
                 "eddy velocity, i.e. sqrt(2*MEKE). This should be less than 1 "//&
                 "to account for the surface intensification of MEKE.", &
                 units="nondim", default=0.)
  call get_param(param_file, mdl, "MEKE_CB", CS%MEKE_Cb, &
                 "A coefficient in the expression for the ratio of bottom projected "//&
                 "eddy energy and mean column energy (see Jansen et al. 2015).",&
                 units="nondim", default=25.)
  call get_param(param_file, mdl, "MEKE_MIN_GAMMA2", CS%MEKE_min_gamma, &
                 "The minimum allowed value of gamma_b^2.",&
                 units="nondim", default=0.0001)
  call get_param(param_file, mdl, "MEKE_CT", CS%MEKE_Ct, &
                 "A coefficient in the expression for the ratio of barotropic "//&
                 "eddy energy and mean column energy (see Jansen et al. 2015).",&
                 units="nondim", default=50.)
  call get_param(param_file, mdl, "MEKE_GMCOEFF", CS%MEKE_GMcoeff, &
                 "The efficiency of the conversion of potential energy "//&
                 "into MEKE by the thickness mixing parameterization. "//&
                 "If MEKE_GMCOEFF is negative, this conversion is not "//&
                 "used or calculated.", units="nondim", default=-1.0)
  call get_param(param_file, mdl, "MEKE_GEOMETRIC", CS%MEKE_GEOMETRIC, &
                 "If MEKE_GEOMETRIC is true, uses the GM coefficient formulation "//&
                 "from the GEOMETRIC framework (Marshall et al., 2012).", default=.false.)
  call get_param(param_file, mdl, "MEKE_GEOMETRIC_ALPHA", CS%MEKE_GEOMETRIC_alpha, &
                 "The nondimensional coefficient governing the efficiency of the GEOMETRIC \n"//&
                 "thickness diffusion.", units="nondim", default=0.05)
  call get_param(param_file, mdl, "MEKE_EQUILIBRIUM_ALT", CS%MEKE_equilibrium_alt, &
                 "If true, use an alternative formula for computing the (equilibrium)"//&
                 "initial value of MEKE.", default=.false.)
  call get_param(param_file, mdl, "MEKE_EQUILIBRIUM_RESTORING", CS%MEKE_equilibrium_restoring, &
                 "If true, restore MEKE back to its equilibrium value, which is calculated at "//&
                 "each time step.", default=.false.)
  if (CS%MEKE_equilibrium_restoring) then
    call get_param(param_file, mdl, "MEKE_RESTORING_TIMESCALE", MEKE_restoring_timescale, &
                   "The timescale used to nudge MEKE toward its equilibrium value.", units="s", &
                   default=1e6, scale=US%T_to_s)
    CS%MEKE_restoring_rate = 1.0 / MEKE_restoring_timescale
  endif

  call get_param(param_file, mdl, "MEKE_FRCOEFF", CS%MEKE_FrCoeff, &
                 "The efficiency of the conversion of mean energy into "//&
                 "MEKE.  If MEKE_FRCOEFF is negative, this conversion "//&
                 "is not used or calculated.", units="nondim", default=-1.0)
  call get_param(param_file, mdl, "MEKE_GMECOEFF", CS%MEKE_GMECoeff, &
                 "The efficiency of the conversion of MEKE into mean energy "//&
                 "by GME.  If MEKE_GMECOEFF is negative, this conversion "//&
                 "is not used or calculated.", units="nondim", default=-1.0)
  call get_param(param_file, mdl, "MEKE_BGSRC", CS%MEKE_BGsrc, &
                 "A background energy source for MEKE.", units="W kg-1", &
                 default=0.0, scale=US%m_to_L**2*US%T_to_s**3)
  call get_param(param_file, mdl, "MEKE_KH", CS%MEKE_Kh, &
                 "A background lateral diffusivity of MEKE. "//&
                 "Use a negative value to not apply lateral diffusion to MEKE.", &
                 units="m2 s-1", default=-1.0, scale=US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "MEKE_K4", CS%MEKE_K4, &
                 "A lateral bi-harmonic diffusivity of MEKE. "//&
                 "Use a negative value to not apply bi-harmonic diffusion to MEKE.", &
                 units="m4 s-1", default=-1.0, scale=US%m_to_L**4*US%T_to_s)
  call get_param(param_file, mdl, "MEKE_DTSCALE", CS%MEKE_dtScale, &
                 "A scaling factor to accelerate the time evolution of MEKE.", &
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "MEKE_KHCOEFF", CS%MEKE_KhCoeff, &
                 "A scaling factor in the expression for eddy diffusivity "//&
                 "which is otherwise proportional to the MEKE velocity- "//&
                 "scale times an eddy mixing-length. This factor "//&
                 "must be >0 for MEKE to contribute to the thickness/ "//&
                 "and tracer diffusivity in the rest of the model.", &
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "MEKE_USCALE", CS%MEKE_Uscale, &
                 "The background velocity that is combined with MEKE to "//&
                 "calculate the bottom drag.", units="m s-1", default=0.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "MEKE_GM_SRC_ALT", CS%GM_src_alt, &
                 "If true, use the GM energy conversion form S^2*N^2*kappa rather "//&
                 "than the streamfunction for the MEKE GM source term.", default=.false.)
  call get_param(param_file, mdl, "MEKE_VISC_DRAG", CS%visc_drag, &
                 "If true, use the vertvisc_type to calculate the bottom "//&
                 "drag acting on MEKE.", default=.true.)
  call get_param(param_file, mdl, "MEKE_KHTH_FAC", MEKE%KhTh_fac, &
                 "A factor that maps MEKE%Kh to KhTh.", units="nondim", &
                 default=0.0)
  call get_param(param_file, mdl, "MEKE_KHTR_FAC", MEKE%KhTr_fac, &
                 "A factor that maps MEKE%Kh to KhTr.", units="nondim", &
                 default=0.0)
  call get_param(param_file, mdl, "MEKE_KHMEKE_FAC", CS%KhMEKE_Fac, &
                 "A factor that maps MEKE%Kh to Kh for MEKE itself.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_OLD_LSCALE", CS%use_old_lscale, &
                 "If true, use the old formula for length scale which is "//&
                 "a function of grid spacing and deformation radius.",  &
                 default=.false.)
  call get_param(param_file, mdl, "MEKE_MIN_LSCALE", CS%use_min_lscale, &
                 "If true, use a strict minimum of provided length scales "//&
                 "rather than harmonic mean.",  &
                 default=.false.)
  call get_param(param_file, mdl, "MEKE_RD_MAX_SCALE", CS%Rd_as_max_scale, &
                 "If true, the length scale used by MEKE is the minimum of "//&
                 "the deformation radius or grid-spacing. Only used if "//&
                 "MEKE_OLD_LSCALE=True", units="nondim", default=.false.)
  call get_param(param_file, mdl, "MEKE_VISCOSITY_COEFF_KU", CS%viscosity_coeff_Ku, &
                 "If non-zero, is the scaling coefficient in the expression for"//&
                 "viscosity used to parameterize harmonic lateral momentum mixing by"//&
                 "unresolved eddies represented by MEKE. Can be negative to"//&
                 "represent backscatter from the unresolved eddies.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_VISCOSITY_COEFF_AU", CS%viscosity_coeff_Au, &
                 "If non-zero, is the scaling coefficient in the expression for"//&
                 "viscosity used to parameterize biharmonic lateral momentum mixing by"//&
                 "unresolved eddies represented by MEKE. Can be negative to"//&
                 "represent backscatter from the unresolved eddies.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_FIXED_MIXING_LENGTH", CS%Lfixed, &
                 "If positive, is a fixed length contribution to the expression "//&
                 "for mixing length used in MEKE-derived diffusivity.", &
                 units="m", default=0.0, scale=US%m_to_L)
  call get_param(param_file, mdl, "MEKE_ALPHA_DEFORM", CS%aDeform, &
                 "If positive, is a coefficient weighting the deformation scale "//&
                 "in the expression for mixing length used in MEKE-derived diffusivity.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ALPHA_RHINES", CS%aRhines, &
                 "If positive, is a coefficient weighting the Rhines scale "//&
                 "in the expression for mixing length used in MEKE-derived diffusivity.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ALPHA_EADY", CS%aEady, &
                 "If positive, is a coefficient weighting the Eady length scale "//&
                 "in the expression for mixing length used in MEKE-derived diffusivity.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ALPHA_FRICT", CS%aFrict, &
                 "If positive, is a coefficient weighting the frictional arrest scale "//&
                 "in the expression for mixing length used in MEKE-derived diffusivity.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ALPHA_GRID", CS%aGrid, &
                 "If positive, is a coefficient weighting the grid-spacing as a scale "//&
                 "in the expression for mixing length used in MEKE-derived diffusivity.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_COLD_START", coldStart, &
                 "If true, initialize EKE to zero. Otherwise a local equilibrium solution "//&
                 "is used as an initial condition for EKE.", default=.false.)
  call get_param(param_file, mdl, "MEKE_BACKSCAT_RO_C", MEKE%backscatter_Ro_c, &
                 "The coefficient in the Rossby number function for scaling the biharmonic "//&
                 "frictional energy source. Setting to non-zero enables the Rossby number function.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_BACKSCAT_RO_POW", MEKE%backscatter_Ro_pow, &
                 "The power in the Rossby number function for scaling the biharmonic "//&
                 "frictional energy source.", units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_ADVECTION_FACTOR", CS%MEKE_advection_factor, &
                 "A scale factor in front of advection of eddy energy. Zero turns advection off. "//&
                 "Using unity would be normal but other values could accommodate a mismatch "//&
                 "between the advecting barotropic flow and the vertical structure of MEKE.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_TOPOGRAPHIC_BETA", CS%MEKE_topographic_beta, &
                 "A scale factor to determine how much topographic beta is weighed in " //&
                 "computing beta in the expression of Rhines scale. Use 1 if full "//&
                 "topographic beta effect is considered; use 0 if it's completely ignored.", &
                 units="nondim", default=0.0)

  ! Nonlocal module parameters
  call get_param(param_file, mdl, "CDRAG", cdrag, &
                 "CDRAG is the drag coefficient relating the magnitude of "//&
                 "the velocity field to the bottom stress.", units="nondim", &
                 default=0.003)
  call get_param(param_file, mdl, "MEKE_CDRAG", CS%cdrag, &
                 "Drag coefficient relating the magnitude of the velocity "//&
                 "field to the bottom stress in MEKE.", units="nondim", &
                 default=cdrag)
  call get_param(param_file, mdl, "LAPLACIAN", laplacian, default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "BIHARMONIC", biharmonic, default=.false., do_not_log=.true.)

  if (CS%viscosity_coeff_Ku/=0. .and. .not. laplacian) call MOM_error(FATAL, &
                 "LAPLACIAN must be true if MEKE_VISCOSITY_COEFF_KU is true.")

  if (CS%viscosity_coeff_Au/=0. .and. .not. biharmonic) call MOM_error(FATAL, &
                 "BIHARMONIC must be true if MEKE_VISCOSITY_COEFF_AU is true.")

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false., do_not_log=.true.)

  ! Identify if any lateral diffusive processes are active
  CS%kh_flux_enabled = .false.
  if ((CS%MEKE_KH >= 0.0)  .or. (CS%KhMEKE_FAC > 0.0) .or. (CS%MEKE_advection_factor > 0.0)) &
    CS%kh_flux_enabled = .true.

! Register fields for output from this module.
  CS%diag => diag
  CS%id_MEKE = register_diag_field('ocean_model', 'MEKE', diag%axesT1, Time, &
     'Mesoscale Eddy Kinetic Energy', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  if (.not. associated(MEKE%MEKE)) CS%id_MEKE = -1
  CS%id_Kh = register_diag_field('ocean_model', 'MEKE_KH', diag%axesT1, Time, &
     'MEKE derived diffusivity', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  if (.not. associated(MEKE%Kh)) CS%id_Kh = -1
  CS%id_Ku = register_diag_field('ocean_model', 'MEKE_KU', diag%axesT1, Time, &
     'MEKE derived lateral viscosity', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  if (.not. associated(MEKE%Ku)) CS%id_Ku = -1
  CS%id_Au = register_diag_field('ocean_model', 'MEKE_AU', diag%axesT1, Time, &
     'MEKE derived lateral biharmonic viscosity', 'm4 s-1', conversion=US%L_to_m**4*US%s_to_T)
  if (.not. associated(MEKE%Au)) CS%id_Au = -1
  CS%id_Ue = register_diag_field('ocean_model', 'MEKE_Ue', diag%axesT1, Time, &
     'MEKE derived eddy-velocity scale', 'm s-1', conversion=US%L_T_to_m_s)
  if (.not. associated(MEKE%MEKE)) CS%id_Ue = -1
  CS%id_Ub = register_diag_field('ocean_model', 'MEKE_Ub', diag%axesT1, Time, &
     'MEKE derived bottom eddy-velocity scale', 'm s-1', conversion=US%L_T_to_m_s)
  if (.not. associated(MEKE%MEKE)) CS%id_Ub = -1
  CS%id_Ut = register_diag_field('ocean_model', 'MEKE_Ut', diag%axesT1, Time, &
     'MEKE derived barotropic eddy-velocity scale', 'm s-1', conversion=US%L_T_to_m_s)
  if (.not. associated(MEKE%MEKE)) CS%id_Ut = -1
  CS%id_src = register_diag_field('ocean_model', 'MEKE_src', diag%axesT1, Time, &
     'MEKE energy source', 'm2 s-3', conversion=(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_decay = register_diag_field('ocean_model', 'MEKE_decay', diag%axesT1, Time, &
     'MEKE decay rate', 's-1', conversion=US%s_to_T)
  CS%id_GM_src = register_diag_field('ocean_model', 'MEKE_GM_src', diag%axesT1, Time, &
     'MEKE energy available from thickness mixing', &
     'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  if (.not. associated(MEKE%GM_src)) CS%id_GM_src = -1
  CS%id_mom_src = register_diag_field('ocean_model', 'MEKE_mom_src',diag%axesT1, Time, &
     'MEKE energy available from momentum', &
     'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  if (.not. associated(MEKE%mom_src)) CS%id_mom_src = -1
  CS%id_GME_snk = register_diag_field('ocean_model', 'MEKE_GME_snk',diag%axesT1, Time, &
     'MEKE energy lost to GME backscatter', &
     'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  if (.not. associated(MEKE%GME_snk)) CS%id_GME_snk = -1
  CS%id_Le = register_diag_field('ocean_model', 'MEKE_Le', diag%axesT1, Time, &
     'Eddy mixing length used in the MEKE derived eddy diffusivity', 'm', conversion=US%L_to_m)
  CS%id_Lrhines = register_diag_field('ocean_model', 'MEKE_Lrhines', diag%axesT1, Time, &
     'Rhines length scale used in the MEKE derived eddy diffusivity', 'm', conversion=US%L_to_m)
  CS%id_Leady = register_diag_field('ocean_model', 'MEKE_Leady', diag%axesT1, Time, &
     'Eady length scale used in the MEKE derived eddy diffusivity', 'm', conversion=US%L_to_m)
  CS%id_gamma_b = register_diag_field('ocean_model', 'MEKE_gamma_b', diag%axesT1, Time, &
     'Ratio of bottom-projected eddy velocity to column-mean eddy velocity', 'nondim')
  CS%id_gamma_t = register_diag_field('ocean_model', 'MEKE_gamma_t', diag%axesT1, Time, &
     'Ratio of barotropic eddy velocity to column-mean eddy velocity', 'nondim')

  if (CS%kh_flux_enabled) then
    CS%id_KhMEKE_u = register_diag_field('ocean_model', 'KHMEKE_u', diag%axesCu1, Time, &
     'Zonal diffusivity of MEKE', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
    CS%id_KhMEKE_v = register_diag_field('ocean_model', 'KHMEKE_v', diag%axesCv1, Time, &
     'Meridional diffusivity of MEKE', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  endif

  if (CS%MEKE_equilibrium_restoring) then
    CS%id_MEKE_equilibrium = register_diag_field('ocean_model', 'MEKE_equilibrium', diag%axesT1, Time, &
     'Equilibrated Mesoscale Eddy Kinetic Energy', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  endif

  CS%id_clock_pass = cpu_clock_id('(Ocean continuity halo updates)', grain=CLOCK_ROUTINE)

  ! Detect whether this instance of MEKE_init() is at the beginning of a run
  ! or after a restart. If at the beginning, we will initialize MEKE to a local
  ! equilibrium.
  CS%initialize = .not.query_initialized(MEKE%MEKE, "MEKE", restart_CS)
  if (coldStart) CS%initialize = .false.
  if (CS%initialize) call MOM_error(WARNING, &
                       "MEKE_init: Initializing MEKE with a local equilibrium balance.")

  ! Account for possible changes in dimensional scaling for variables that have been
  ! read from a restart file.
  I_T_rescale = 1.0
  if ((US%s_to_T_restart /= 0.0) .and. (US%s_to_T_restart /= US%s_to_T)) &
    I_T_rescale = US%s_to_T_restart / US%s_to_T
  L_rescale = 1.0
  if ((US%m_to_L_restart /= 0.0) .and. (US%m_to_L_restart /= US%m_to_L)) &
    L_rescale = US%m_to_L / US%m_to_L_restart

  if (L_rescale*I_T_rescale /= 1.0) then
    if (associated(MEKE%MEKE)) then ; if (query_initialized(MEKE%MEKE, "MEKE_MEKE", restart_CS)) then
      do j=js,je ; do i=is,ie
        MEKE%MEKE(i,j) = L_rescale*I_T_rescale * MEKE%MEKE(i,j)
      enddo ; enddo
    endif ; endif
  endif
  if (L_rescale**2*I_T_rescale /= 1.0) then
    if (associated(MEKE%Kh)) then ; if (query_initialized(MEKE%Kh, "MEKE_Kh", restart_CS)) then
      do j=js,je ; do i=is,ie
        MEKE%Kh(i,j) = L_rescale**2*I_T_rescale * MEKE%Kh(i,j)
      enddo ; enddo
    endif ; endif
    if (associated(MEKE%Ku)) then ; if (query_initialized(MEKE%Ku, "MEKE_Ku", restart_CS)) then
      do j=js,je ; do i=is,ie
        MEKE%Ku(i,j) = L_rescale**2*I_T_rescale * MEKE%Ku(i,j)
      enddo ; enddo
    endif ; endif
    if (associated(MEKE%Kh_diff)) then ; if (query_initialized(MEKE%Kh, "MEKE_Kh_diff", restart_CS)) then
      do j=js,je ; do i=is,ie
        MEKE%Kh_diff(i,j) = L_rescale**2*I_T_rescale * MEKE%Kh_diff(i,j)
      enddo ; enddo
    endif ; endif
  endif
  if (L_rescale**4*I_T_rescale /= 1.0) then
    if (associated(MEKE%Au)) then ; if (query_initialized(MEKE%Au, "MEKE_Au", restart_CS)) then
      do j=js,je ; do i=is,ie
        MEKE%Au(i,j) = L_rescale**4*I_T_rescale * MEKE%Au(i,j)
      enddo ; enddo
    endif ; endif
  endif

  ! Set up group passes.  In the case of a restart, these fields need a halo update now.
  if (associated(MEKE%MEKE)) then
    call create_group_pass(CS%pass_MEKE, MEKE%MEKE, G%Domain)
    if (associated(MEKE%Kh_diff)) call create_group_pass(CS%pass_MEKE, MEKE%Kh_diff, G%Domain)
    if (.not.CS%initialize) call do_group_pass(CS%pass_MEKE, G%Domain)
  endif
  if (associated(MEKE%Kh)) call create_group_pass(CS%pass_Kh, MEKE%Kh, G%Domain)
  if (associated(MEKE%Ku)) call create_group_pass(CS%pass_Kh, MEKE%Ku, G%Domain)
  if (associated(MEKE%Au)) call create_group_pass(CS%pass_Kh, MEKE%Au, G%Domain)

  if (associated(MEKE%Kh) .or. associated(MEKE%Ku) .or. associated(MEKE%Au)) &
    call do_group_pass(CS%pass_Kh, G%Domain)

end function MEKE_init

!> Allocates memory and register restart fields for the MOM_MEKE module.
subroutine MEKE_alloc_register_restart(HI, param_file, MEKE, restart_CS)
! Arguments
  type(hor_index_type),  intent(in)    :: HI         !< Horizontal index structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file parser structure.
  type(MEKE_type),       pointer       :: MEKE       !< A structure with MEKE-related fields.
  type(MOM_restart_CS),  pointer       :: restart_CS !< Restart control structure for MOM_MEKE.
! Local variables
  type(vardesc) :: vd
  real :: MEKE_GMcoeff, MEKE_FrCoeff, MEKE_GMECoeff, MEKE_KHCoeff, MEKE_viscCoeff_Ku, MEKE_viscCoeff_Au
  logical :: Use_KH_in_MEKE
  logical :: useMEKE
  integer :: isd, ied, jsd, jed

! Determine whether this module will be used
  useMEKE = .false.; call read_param(param_file,"USE_MEKE",useMEKE)

! Read these parameters to determine what should be in the restarts
  MEKE_GMcoeff =-1.; call read_param(param_file,"MEKE_GMCOEFF",MEKE_GMcoeff)
  MEKE_FrCoeff =-1.; call read_param(param_file,"MEKE_FRCOEFF",MEKE_FrCoeff)
  MEKE_GMEcoeff =-1.; call read_param(param_file,"MEKE_GMECOEFF",MEKE_GMEcoeff)
  MEKE_KhCoeff =1.; call read_param(param_file,"MEKE_KHCOEFF",MEKE_KhCoeff)
  MEKE_viscCoeff_Ku =0.; call read_param(param_file,"MEKE_VISCOSITY_COEFF_KU",MEKE_viscCoeff_Ku)
  MEKE_viscCoeff_Au =0.; call read_param(param_file,"MEKE_VISCOSITY_COEFF_AU",MEKE_viscCoeff_Au)
  Use_KH_in_MEKE = .false.; call read_param(param_file,"USE_KH_IN_MEKE", Use_KH_in_MEKE)
! Allocate control structure
  if (associated(MEKE)) then
    call MOM_error(WARNING, "MEKE_alloc_register_restart called with an associated "// &
                             "MEKE type.")
    return
  else; allocate(MEKE); endif

  if (.not. useMEKE) return

! Allocate memory
  call MOM_mesg("MEKE_alloc_register_restart: allocating and registering", 5)
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed
  allocate(MEKE%MEKE(isd:ied,jsd:jed)) ; MEKE%MEKE(:,:) = 0.0
  vd = var_desc("MEKE", "m2 s-2", hor_grid='h', z_grid='1', &
           longname="Mesoscale Eddy Kinetic Energy")
  call register_restart_field(MEKE%MEKE, vd, .false., restart_CS)
  if (MEKE_GMcoeff>=0.) then
    allocate(MEKE%GM_src(isd:ied,jsd:jed)) ; MEKE%GM_src(:,:) = 0.0
  endif
  if (MEKE_FrCoeff>=0. .or. MEKE_GMECoeff>=0.)  then
    allocate(MEKE%mom_src(isd:ied,jsd:jed)) ; MEKE%mom_src(:,:) = 0.0
  endif
  if (MEKE_GMECoeff>=0.) then
    allocate(MEKE%GME_snk(isd:ied,jsd:jed)) ; MEKE%GME_snk(:,:) = 0.0
  endif
  if (MEKE_KhCoeff>=0.) then
    allocate(MEKE%Kh(isd:ied,jsd:jed)) ; MEKE%Kh(:,:) = 0.0
    vd = var_desc("MEKE_Kh", "m2 s-1", hor_grid='h', z_grid='1', &
             longname="Lateral diffusivity from Mesoscale Eddy Kinetic Energy")
    call register_restart_field(MEKE%Kh, vd, .false., restart_CS)
  endif
  allocate(MEKE%Rd_dx_h(isd:ied,jsd:jed)) ; MEKE%Rd_dx_h(:,:) = 0.0
  if (MEKE_viscCoeff_Ku/=0.) then
    allocate(MEKE%Ku(isd:ied,jsd:jed)) ; MEKE%Ku(:,:) = 0.0
    vd = var_desc("MEKE_Ku", "m2 s-1", hor_grid='h', z_grid='1', &
             longname="Lateral viscosity from Mesoscale Eddy Kinetic Energy")
    call register_restart_field(MEKE%Ku, vd, .false., restart_CS)
  endif
  if (Use_Kh_in_MEKE) then
    allocate(MEKE%Kh_diff(isd:ied,jsd:jed)) ; MEKE%Kh_diff(:,:) = 0.0
    vd = var_desc("MEKE_Kh_diff", "m2 s-1",hor_grid='h',z_grid='1', &
             longname="Copy of thickness diffusivity for diffusing MEKE")
    call register_restart_field(MEKE%Kh_diff, vd, .false., restart_CS)
  endif

  if (MEKE_viscCoeff_Au/=0.) then
    allocate(MEKE%Au(isd:ied,jsd:jed)) ; MEKE%Au(:,:) = 0.0
    vd = var_desc("MEKE_Au", "m4 s-1", hor_grid='h', z_grid='1', &
             longname="Lateral biharmonic viscosity from Mesoscale Eddy Kinetic Energy")
    call register_restart_field(MEKE%Au, vd, .false., restart_CS)
  endif

end subroutine MEKE_alloc_register_restart

!> Deallocates any variables allocated in MEKE_init or
!! MEKE_alloc_register_restart.
subroutine MEKE_end(MEKE, CS)
  type(MEKE_type), pointer :: MEKE !< A structure with MEKE-related fields.
  type(MEKE_CS),   pointer :: CS   !< The control structure for MOM_MEKE.

  if (associated(CS)) deallocate(CS)

  if (.not.associated(MEKE)) return

  if (associated(MEKE%MEKE)) deallocate(MEKE%MEKE)
  if (associated(MEKE%GM_src)) deallocate(MEKE%GM_src)
  if (associated(MEKE%mom_src)) deallocate(MEKE%mom_src)
  if (associated(MEKE%GME_snk)) deallocate(MEKE%GME_snk)
  if (associated(MEKE%Kh)) deallocate(MEKE%Kh)
  if (associated(MEKE%Kh_diff)) deallocate(MEKE%Kh_diff)
  if (associated(MEKE%Ku)) deallocate(MEKE%Ku)
  if (associated(MEKE%Au)) deallocate(MEKE%Au)
  deallocate(MEKE)

end subroutine MEKE_end

!> \namespace mom_meke
!!
!! \section section_MEKE The Mesoscale Eddy Kinetic Energy (MEKE) framework
!!
!! The MEKE framework accounts for the mean potential energy removed by
!! the first order closures used to parameterize mesoscale eddies.
!! It requires closure at the second order, namely dissipation and transport
!! of eddy energy.
!!
!! Monitoring the sub-grid scale eddy energy budget provides a means to predict
!! a sub-grid eddy-velocity scale which can be used in the lower order closures.
!!
!! \subsection section_MEKE_equations MEKE equations
!!
!! The eddy kinetic energy equation is:
!! \f[ \partial_{\tilde{t}} E =
!!   \overbrace{ \dot{E}_b + \gamma_\eta \dot{E}_\eta + \gamma_v \dot{E}_v
!!             }^\text{sources}
!! - \overbrace{ ( \lambda + C_d | U_d | \gamma_b^2 ) E
!!             }^\text{local dissipation}
!! + \overbrace{ \nabla \cdot ( ( \kappa_E + \gamma_M \kappa_M ) \nabla E
!!                              - \kappa_4 \nabla^3 E )
!!             }^\text{smoothing}
!! \f]
!! where \f$ E \f$ is the eddy kinetic energy (variable <code>MEKE</code>) with units of
!! m<sup>2</sup>s<sup>-2</sup>,
!! and \f$\tilde{t} = a t\f$ is a scaled time. The non-dimensional factor
!! \f$ a\geq 1 \f$ is used to accelerate towards equilibrium.
!!
!! The MEKE equation is two-dimensional and obtained by depth averaging the
!! the three-dimensional eddy energy equation. In the following expressions
!! \f$ \left< \phi \right> = \frac{1}{H} \int^\eta_{-D} \phi \, dz \f$ maps
!! three dimensional terms into the two-dimensional quantities needed.
!!
!! \subsubsection section_MEKE_source_terms MEKE source terms
!!
!! The source term \f$ \dot{E}_b \f$ is a constant background source
!! of energy intended to avoid the limit \f$E\rightarrow 0\f$.
!!
!! The "GM" source term
!! \f[ \dot{E}_\eta = - \left< \overline{w^\prime b^\prime} \right>
!! = \left< \kappa_h N^2S^2 \right>
!! \approx \left< \kappa_h g\prime |\nabla_\sigma \eta|^2 \right>\f]
!! equals the mean potential energy removed by the Gent-McWilliams closure,
!! and is excluded/included in the MEKE budget by the efficiency parameter
!! \f$ \gamma_\eta \in [0,1] \f$.
!!
!! The "frictional" source term
!! \f[ \dot{E}_{v} = \left<  \partial_i u_j \tau_{ij} \right> \f]
!! equals the mean kinetic energy removed by lateral viscous fluxes, and
!! is excluded/included in the MEKE budget by the efficiency parameter
!! \f$ \gamma_v \in [0,1] \f$.
!!
!! \subsubsection section_MEKE_dissipation_terms MEKE dissipation terms
!!
!! The local dissipation of \f$ E \f$ is parameterized through a linear
!! damping, \f$\lambda\f$, and bottom drag, \f$ C_d | U_d | \gamma_b^2 \f$.
!! The \f$ \gamma_b \f$ accounts for the weak projection of the column-mean
!! eddy velocty to the bottom. In other words, the bottom velocity is
!! estimated as \f$ \gamma_b U_e \f$.
!! The bottom drag coefficient, \f$ C_d \f$ is the same as that used in the bottom
!! friction in the mean model equations.
!!
!! The bottom drag velocity scale, \f$ U_d \f$, has contributions from the
!! resolved state and \f$ E \f$:
!! \f[ U_d = \sqrt{ U_b^2 + |u|^2_{z=-D} + |\gamma_b U_e|^2 } .\f]
!! where the eddy velocity scale, \f$ U_e \f$, is given by:
!! \f[ U_e = \sqrt{ 2 E } .\f]
!! \f$ U_b \f$ is a constant background bottom velocity scale and is
!! typically not used (i.e. set to zero).
!!
!! Following Jansen et al., 2015, the projection of eddy energy on to the bottom
!! is given by the ratio of bottom energy to column mean energy:
!! \f[
!! \gamma_b^2  = \frac{E_b}{E} = \gamma_{d0}
!!    + \left( 1 + c_{b} \frac{L_d}{L_f} \right)^{-\frac{4}{5}}
!! ,
!! \f]
!! \f[
!! \gamma_b^2  \leftarrow  \max{\left( \gamma_b^2, \gamma_{min}^2 \right)}
!! .
!! \f]
!!
!! \subsection section_MEKE_smoothing MEKE smoothing terms
!!
!! \f$ E \f$ is laterally diffused by a diffusivity \f$ \kappa_E + \gamma_M
!! \kappa_M \f$ where \f$ \kappa_E \f$ is a constant diffusivity and the term
!! \f$ \gamma_M \kappa_M \f$ is a "self diffusion" using the diffusivity
!! calculated in the section \ref section_MEKE_diffusivity.
!! \f$ \kappa_4 \f$ is a constant bi-harmonic diffusivity.
!!
!! \subsection section_MEKE_diffusivity Diffusivity derived from MEKE
!!
!! The predicted eddy velocity scale, \f$ U_e \f$, can be combined with a
!! mixing length scale to form a diffusivity.
!! The primary use of a MEKE derived diffusivity is for use in thickness
!! diffusion (module mom_thickness_diffuse) and optionally in along
!! isopycnal mixing of tracers (module mom_tracer_hor_diff).
!! The original form used (enabled with MEKE_OLD_LSCALE=True):
!!
!! \f[  \kappa_M = \gamma_\kappa \sqrt{ \gamma_t^2 U_e^2 A_\Delta } \f]
!!
!! where \f$ A_\Delta \f$ is the area of the grid cell.
!! Following Jansen et al., 2015, we now use
!!
!! \f[  \kappa_M = \gamma_\kappa l_M \sqrt{ \gamma_t^2 U_e^2 } \f]
!!
!! where \f$ \gamma_\kappa \in [0,1] \f$ is a non-dimensional factor and,
!! following Jansen et al., 2015, \f$\gamma_t^2\f$ is the ratio of barotropic
!! eddy energy to column mean eddy energy given by
!! \f[
!! \gamma_t^2  = \frac{E_t}{E} = \left( 1 + c_{t} \frac{L_d}{L_f} \right)^{-\frac{1}{4}}
!! ,
!! \f]
!! \f[
!! \gamma_t^2  \leftarrow  \max{\left( \gamma_t^2, \gamma_{min}^2 \right)}
!! .
!! \f]
!!
!! The length-scale is a configurable combination of multiple length scales:
!!
!! \f[
!! l_M = \left(
!!       \frac{\alpha_d}{L_d}
!!     + \frac{\alpha_f}{L_f}
!!     + \frac{\alpha_R}{L_R}
!!     + \frac{\alpha_e}{L_e}
!!     + \frac{\alpha_\Delta}{L_\Delta}
!!     + \frac{\delta[L_c]}{L_c}
!!       \right)^{-1}
!! \f]
!!
!! where
!!
!! \f{eqnarray*}{
!! L_d & = & \sqrt{\frac{c_g^2}{f^2+2\beta c_g}} \sim \frac{ c_g }{f} \\\\
!! L_R & = & \sqrt{\frac{U_e}{\beta^*}} \\\\
!! L_e & = & \frac{U_e}{|S| N} \\\\
!! L_f & = & \frac{H}{c_d} \\\\
!! L_\Delta & = & \sqrt{A_\Delta} .
!! \f}
!!
!! \f$L_c\f$ is a constant and \f$\delta[L_c]\f$ is the impulse function so that the term
!! \f$\frac{\delta[L_c]}{L_c}\f$ evaluates to \f$\frac{1}{L_c}\f$ when \f$L_c\f$ is non-zero
!! but is dropped if \f$L_c=0\f$.
!!
!! \f$\beta^*\f$ is the effective \f$\beta\f$ that combines both the planetary vorticity
!! gradient (i.e. \f$\beta=\nabla f\f$) and the topographic \f$\beta\f$ effect,
!! with the latter weighed by a weighting constant, \f$c_\beta\f$, that varies
!! from 0 to 1, so that \f$c_\beta=0\f$ means the topographic \f$\beta\f$ effect is ignored,
!! while \f$c_\beta=1\f$ means it is fully considered. The new \f$\beta^*\f$ therefore
!! takes the form of
!!
!! \f[
!! \beta^* = \sqrt{( \partial_xf - c_\beta\frac{f}{D}\partial_xD )^2 +
!!           ( \partial_yf - c_\beta\frac{f}{D}\partial_yD )^2}
!! \f]
!! where \f$D\f$ is water column depth at T points.
!!
!! \subsection section_MEKE_viscosity Viscosity derived from MEKE
!!
!! As for \f$ \kappa_M \f$, the predicted eddy velocity scale can be
!! used to form a harmonic eddy viscosity,
!!
!! \f[  \kappa_u = \gamma_u \sqrt{ U_e^2 A_\Delta }  \f]
!!
!! as well as a biharmonic eddy viscosity,
!!
!! \f[  \kappa_4 = \gamma_4 \sqrt{ U_e^2 A_\Delta^3 }  \f]
!!
!! \subsection section_MEKE_limit_case Limit cases for local source-dissipative balance
!!
!! Note that in steady-state (or when \f$ a>>1 \f$) and there is no
!! diffusion of \f$ E \f$ then
!! \f[ \overline{E} \approx \frac{ \dot{E}_b + \gamma_\eta \dot{E}_\eta +
!!               \gamma_v \dot{E}_v }{ \lambda + C_d|U_d|\gamma_b^2 } . \f]
!!
!! In the linear drag limit, where
!! \f$ U_e << \min(U_b, |u|_{z=-D}, C_d^{-1}\lambda) \f$, the equilibrium becomes
!! \f$ \overline{E} \approx \frac{ \dot{E}_b + \gamma_\eta \dot{E}_\eta +
!!               \gamma_v \dot{E}_v }{ \lambda + C_d \sqrt{ U_b^2 + |u|^2_{z=-D} } } \f$.
!!
!! In the nonlinear drag limit, where \f$ U_e >> \max(U_b, |u|_{z=-D}, C_d^{-1}\lambda) \f$,
!! the equilibrium becomes
!! \f$ \overline{E} \approx \left( \frac{ \dot{E}_b + \gamma_\eta \dot{E}_\eta +
!!               \gamma_v \dot{E}_v }{ \sqrt{2} C_d \gamma_b^3 } \right)^\frac{2}{3} \f$.
!!
!! \subsubsection section_MEKE_module_parameters MEKE module parameters
!!
!! | Symbol                | Module parameter |
!! | ------                | --------------- |
!! | -                     | <code>USE_MEKE</code> |
!! | \f$ a \f$             | <code>MEKE_DTSCALE</code> |
!! | \f$ \dot{E}_b \f$     | <code>MEKE_BGSRC</code> |
!! | \f$ \gamma_\eta \f$   | <code>MEKE_GMCOEFF</code> |
!! | \f$ \gamma_v \f$      | <code>MEKE_FrCOEFF</code> |
!! | \f$ \lambda \f$       | <code>MEKE_DAMPING</code> |
!! | \f$ U_b \f$           | <code>MEKE_USCALE</code> |
!! | \f$ \gamma_{d0} \f$   | <code>MEKE_CD_SCALE</code> |
!! | \f$ c_{b} \f$         | <code>MEKE_CB</code> |
!! | \f$ c_{t} \f$         | <code>MEKE_CT</code> |
!! | \f$ \kappa_E \f$      | <code>MEKE_KH</code> |
!! | \f$ \kappa_4 \f$      | <code>MEKE_K4</code> |
!! | \f$ \gamma_\kappa \f$ | <code>MEKE_KHCOEFF</code> |
!! | \f$ \gamma_M \f$      | <code>MEKE_KHMEKE_FAC</code> |
!! | \f$ \gamma_u \f$      | <code>MEKE_VISCOSITY_COEFF_KU</code> |
!! | \f$ \gamma_4 \f$      | <code>MEKE_VISCOSITY_COEFF_AU</code> |
!! | \f$ \gamma_{min}^2 \f$| <code>MEKE_MIN_GAMMA2</code> |
!! | \f$ \alpha_d \f$      | <code>MEKE_ALPHA_DEFORM</code> |
!! | \f$ \alpha_f \f$      | <code>MEKE_ALPHA_FRICT</code> |
!! | \f$ \alpha_R \f$      | <code>MEKE_ALPHA_RHINES</code> |
!! | \f$ \alpha_e \f$      | <code>MEKE_ALPHA_EADY</code> |
!! | \f$ \alpha_\Delta \f$ | <code>MEKE_ALPHA_GRID</code> |
!! | \f$ L_c \f$           | <code>MEKE_FIXED_MIXING_LENGTH</code> |
!! | \f$ c_\beta \f$       | <code>MEKE_TOPOGRAPHIC_BETA</code> |
!! | -                     | <code>MEKE_KHTH_FAC</code> |
!! | -                     | <code>MEKE_KHTR_FAC</code> |
!!
!! | Symbol                | Model parameter |
!! | ------                | --------------- |
!! | \f$ C_d \f$           | <code>CDRAG</code> |
!!
!! \subsection section_MEKE_references References
!!
!! Jansen, M. F., A. J. Adcroft, R. Hallberg, and I. M. Held, 2015: Parameterization of eddy fluxes based on a
!! mesoscale energy budget. Ocean Modelling, 92, 28--41, http://doi.org/10.1016/j.ocemod.2015.05.007 .
!!
!! Marshall, D. P., and A. J. Adcroft, 2010: Parameterization of ocean eddies: Potential vorticity mixing, energetics
!! and Arnold first stability theorem. Ocean Modelling, 32, 188--204, http://doi.org/10.1016/j.ocemod.2010.02.001 .

end module MOM_MEKE

