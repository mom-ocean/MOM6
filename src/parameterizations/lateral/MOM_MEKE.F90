!> Implements the Mesoscale Eddy Kinetic Energy framework
!! with topographic beta effect included in computing beta in Rhines scale

module MOM_MEKE

! This file is part of MOM6. See LICENSE.md for the license.
use iso_fortran_env,       only : real32

use MOM_coms,              only : PE_here
use MOM_database_comms,    only : dbclient_type, dbcomms_CS_type
use MOM_debugging,         only : hchksum, uvchksum
use MOM_cpu_clock,         only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator,     only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,     only : diag_ctrl, time_type
use MOM_domains,           only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,           only : pass_vector, pass_var
use MOM_error_handler,     only : MOM_error, FATAL, WARNING, NOTE, MOM_mesg, is_root_pe
use MOM_file_parser,       only : read_param, get_param, log_version, param_file_type
use MOM_grid,              only : ocean_grid_type
use MOM_hor_index,         only : hor_index_type
use MOM_interface_heights, only : find_eta
use MOM_interpolate,       only : init_external_field, time_interp_external
use MOM_interpolate,       only : time_interp_external_init
use MOM_interpolate,       only : external_field
use MOM_io,                only : vardesc, var_desc, slasher
use MOM_isopycnal_slopes,  only : calc_isoneutral_slopes
use MOM_restart,           only : MOM_restart_CS, register_restart_field, query_initialized
use MOM_string_functions,  only : lowercase
use MOM_time_manager,      only : time_type_to_real
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : vertvisc_type, thermo_var_ptrs
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_MEKE_types,        only : MEKE_type


implicit none ; private

#include <MOM_memory.h>

public step_forward_MEKE, MEKE_init, MEKE_alloc_register_restart, MEKE_end

! Constants for this module
integer, parameter :: NUM_FEATURES = 4 !< How many features used to predict EKE
integer, parameter :: MKE_IDX = 1     !< Index of mean kinetic energy in the feature array
integer, parameter :: SLOPE_Z_IDX = 2 !< Index of vertically averaged isopycnal slope in the feature array
integer, parameter :: RV_IDX = 3      !< Index of surface relative vorticity in the feature array
integer, parameter :: RD_DX_Z_IDX = 4 !< Index of the radius of deformation over the grid size in the feature array

integer, parameter :: EKE_PROG = 1     !< Use prognostic equation to calculate EKE
integer, parameter :: EKE_FILE = 2     !< Read in EKE from a file
integer, parameter :: EKE_DBCLIENT = 3 !< Infer EKE using a neural network

!> Control structure that contains MEKE parameters and diagnostics handles
type, public :: MEKE_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  ! Parameters
  real :: MEKE_FrCoeff  !< Efficiency of conversion of ME into MEKE [nondim]
  real :: MEKE_GMcoeff  !< Efficiency of conversion of PE into MEKE [nondim]
  real :: MEKE_GMECoeff !< Efficiency of conversion of MEKE into ME by GME [nondim]
  real :: MEKE_damping  !< Local depth-independent MEKE dissipation rate [T-1 ~> s-1].
  real :: MEKE_Cd_scale !< The ratio of the bottom eddy velocity to the column mean
                        !! eddy velocity, i.e. sqrt(2*MEKE), [nondim]. This should be less than 1
                        !! to account for the surface intensification of MEKE.
  real :: MEKE_Cb       !< Coefficient in the \f$\gamma_{bot}\f$ expression [nondim]
  real :: MEKE_min_gamma!< Minimum value of gamma_b^2 allowed [nondim]
  real :: MEKE_Ct       !< Coefficient in the \f$\gamma_{bt}\f$ expression [nondim]
  logical :: visc_drag  !< If true use the vertvisc_type to calculate bottom drag.
  logical :: MEKE_GEOMETRIC !< If true, uses the GM coefficient formulation from the GEOMETRIC
                        !! framework (Marshall et al., 2012)
  real    :: MEKE_GEOMETRIC_alpha !< The nondimensional coefficient governing the efficiency of the
                        !! GEOMETRIC thickness diffusion [nondim].
  logical :: MEKE_equilibrium_alt !< If true, use an alternative calculation for the
                        !! equilibrium value of MEKE.
  logical :: MEKE_equilibrium_restoring !< If true, restore MEKE back to its equilibrium value,
                        !!  which is calculated at each time step.
  logical :: GM_src_alt !< If true, use the GM energy conversion form S^2*N^2*kappa rather
                        !! than the streamfunction for the MEKE GM source term.
  real    :: MEKE_min_depth_tot  !< The minimum total thickness over which to distribute MEKE energy
                        !! sources from GM energy conversion [H ~> m or kg m-2].  When the total
                        !! thickness is less than this, the sources are scaled away.
  logical :: Rd_as_max_scale !< If true the length scale can not exceed the
                        !! first baroclinic deformation radius.
  logical :: use_old_lscale !< Use the old formula for mixing length scale.
  logical :: use_min_lscale !< Use simple minimum for mixing length scale.
  real :: lscale_maxval !< The ceiling on the MEKE mixing length scale when use_min_lscale is true [L ~> m].
  real :: cdrag         !< The bottom drag coefficient for MEKE, times rescaling factors [H L-1 ~> nondim or kg m-3]
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
                        !! by unresolved eddies represented by MEKE [nondim].
  real :: viscosity_coeff_Au !< The scaling coefficient in the expression for
                        !! viscosity used to parameterize lateral biharmonic momentum mixing
                        !! by unresolved eddies represented by MEKE [nondim].
  real :: Lfixed        !< Fixed mixing length scale [L ~> m].
  real :: aDeform       !< Weighting towards deformation scale of mixing length [nondim]
  real :: aRhines       !< Weighting towards Rhines scale of mixing length [nondim]
  real :: aFrict        !< Weighting towards frictional arrest scale of mixing length [nondim]
  real :: aEady         !< Weighting towards Eady scale of mixing length [nondim]
  real :: aGrid         !< Weighting towards grid scale of mixing length [nondim]
  real :: MEKE_advection_factor !< A scaling in front of the advection of MEKE [nondim]
  real :: MEKE_topographic_beta !< Weight for how much topographic beta is considered
                                !! when computing beta in Rhines scale [nondim]
  real :: MEKE_restoring_rate !< Inverse of the timescale used to nudge MEKE toward its
                        !! equilibrium value [T-1 ~> s-1].
  logical :: MEKE_advection_bug !< If true, recover a bug in the calculation of the barotropic
                        !! transport for the advection of MEKE, wherein only the transports in the
                        !! deepest layer are used.
  logical :: fixed_total_depth  !< If true, use the nominal bathymetric depth as the estimate of
                        !! the time-varying ocean depth.  Otherwise base the depth on the total
                        !! ocean mass per unit area.
  real :: rho_fixed_total_depth !< A density used to translate the nominal bathymetric depth into an
                        !! estimate of the total ocean mass per unit area when MEKE_FIXED_TOTAL_DEPTH
                        !! is true [R ~> kg m-3]
  logical :: kh_flux_enabled !< If true, lateral diffusive MEKE flux is enabled.
  logical :: initialize !< If True, invokes a steady state solver to calculate MEKE.
  logical :: debug      !< If true, write out checksums of data for debugging
  integer :: eke_src !< Enum specifying whether EKE is stepped forward prognostically (default),
                     !! read in from a file, or inferred via a neural network
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
  type(external_field) :: eke_handle   !< Handle for reading in EKE from a file
  ! Infrastructure
  integer :: id_clock_pass !< Clock for group pass calls
  type(group_pass_type) :: pass_MEKE !< Group halo pass handle for MEKE%MEKE and maybe MEKE%Kh_diff
  type(group_pass_type) :: pass_Kh   !< Group halo pass handle for MEKE%Kh, MEKE%Ku, and/or MEKE%Au

  ! MEKE via Machine Learning
  type(dbclient_type), pointer :: client => NULL() !< Pointer to the database client

  logical :: online_analysis !< If true, post the EKE used in MOM6 at every timestep
  character(len=5) :: model_key  = 'mleke'  !< Key where the ML-model is stored
  character(len=7) :: key_suffix !< Suffix appended to every key sent to Redis
  real :: eke_max !< The maximum value of EKE considered physically reasonable [L2 T-2 ~> m2 s-2]

  ! Clock ids
  integer :: id_client_init   !< Clock id to time initialization of the client
  integer :: id_put_tensor    !< Clock id to time put_tensor routine
  integer :: id_run_model     !< Clock id to time running of the ML model
  integer :: id_unpack_tensor !< Clock id to time retrieval of EKE prediction

  ! Diagnostic ids
  integer :: id_mke     = -1 !< Diagnostic id for surface mean kinetic energy
  integer :: id_slope_z = -1 !< Diagnostic id for vertically averaged horizontal slope magnitude
  integer :: id_slope_x = -1 !< Diagnostic id for isopycnal slope in the x-direction
  integer :: id_slope_y = -1 !< Diagnostic id for isopycnal slope in the y-direction
  integer :: id_rv      = -1 !< Diagnostic id for surface relative vorticity

end type MEKE_CS

contains

!> Integrates forward-in-time the MEKE eddy energy equation.
!! See \ref section_MEKE_equations.
subroutine step_forward_MEKE(MEKE, h, SN_u, SN_v, visc, dt, G, GV, US, CS, hu, hv, u, v, tv, Time)
  type(MEKE_type),                          intent(inout) :: MEKE !< MEKE data.
  type(ocean_grid_type),                    intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),                  intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),                    intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)   :: h    !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)),        intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)),        intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  type(vertvisc_type),                      intent(in)    :: visc !< The vertical viscosity type.
  real,                                     intent(in)    :: dt   !< Model(baroclinic) time-step [T ~> s].
  type(MEKE_CS),                            intent(inout) :: CS   !< MEKE control structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in)  :: hu   !< Accumulated zonal mass flux [H L2 ~> m3 or kg].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in)  :: hv   !< Accumulated meridional mass flux [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v  !< Meridional velocity [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                    intent(in)    :: tv   !< Type containing thermodynamic variables
  type(time_type),                          intent(in)    :: Time !< The time used for interpolating EKE

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    data_eke, &     ! EKE from file [L2 T-2 ~> m2 s-2]
    mass, &         ! The total mass of the water column [R Z ~> kg m-2].
    I_mass, &       ! The inverse of mass [R-1 Z-1 ~> m2 kg-1].
    depth_tot, &    ! The depth of the water column [H ~> m or kg m-2].
    src, &          ! The sum of all MEKE sources [L2 T-3 ~> W kg-1] (= m2 s-3).
    MEKE_decay, &   ! A diagnostic of the MEKE decay timescale [T-1 ~> s-1].
    drag_rate_visc, & ! Near-bottom velocity contribution to bottom drag [H T-1 ~> m s-1 or kg m-2 s-1]
    drag_rate, &    ! The MEKE spindown timescale due to bottom drag [T-1 ~> s-1].
    del2MEKE, &     ! Laplacian of MEKE, used for bi-harmonic diffusion [T-2 ~> s-2].
    del4MEKE, &     ! Time-integrated MEKE tendency arising from the biharmonic of MEKE [L2 T-2 ~> m2 s-2].
    LmixScale, &    ! Eddy mixing length [L ~> m].
    barotrFac2, &   ! Ratio of EKE_barotropic / EKE [nondim]
    bottomFac2, &   ! Ratio of EKE_bottom / EKE [nondim]
    tmp, &          ! Temporary variable for computation of diagnostic velocities [L T-1 ~> m s-1]
    equilibrium_value ! The equilibrium value of MEKE to be calculated at each
                    ! time step [L2 T-2 ~> m2 s-2]

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    MEKE_uflux, &   ! The zonal advective and diffusive flux of MEKE with units of [R Z L4 T-3 ~> kg m2 s-3].
                    ! In one place, MEKE_uflux is used as temporary work space with units of [L2 T-2 ~> m2 s-2].
    Kh_u, &         ! The zonal diffusivity that is actually used [L2 T-1 ~> m2 s-1].
    baroHu, &       ! Depth integrated accumulated zonal mass flux [R Z L2 ~> kg].
    drag_vel_u      ! A piston velocity associated with bottom drag at u-points [H T-1 ~> m s-1 or kg m-2 s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    MEKE_vflux, &   ! The meridional advective and diffusive flux of MEKE with units of [R Z L4 T-3 ~> kg m2 s-3].
                    ! In one place, MEKE_vflux is used as temporary work space with units of [L2 T-2 ~> m2 s-2].
    Kh_v, &         ! The meridional diffusivity that is actually used [L2 T-1 ~> m2 s-1].
    baroHv, &       ! Depth integrated accumulated meridional mass flux [R Z L2 ~> kg].
    drag_vel_v      ! A piston velocity associated with bottom drag at v-points [H T-1 ~> m s-1 or kg m-2 s-1]
  real :: Kh_here   ! The local horizontal viscosity [L2 T-1 ~> m2 s-1]
  real :: Inv_Kh_max ! The inverse of the local horizontal viscosity [T L-2 ~> s m-2]
  real :: K4_here   ! The local horizontal biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: Inv_K4_max ! The inverse of the local horizontal biharmonic viscosity [T L-4 ~> s m-4]
  real :: cdrag2    ! The square of the drag coefficient times unit conversion factors [H2 L-2 ~> nondim or kg2 m-6]
  real :: advFac    ! The product of the advection scaling factor and 1/dt [T-1 ~> s-1]
  real :: mass_neglect ! A negligible mass [R Z ~> kg m-2].
  real :: ldamping  ! The MEKE damping rate [T-1 ~> s-1].
  real :: sdt       ! dt to use locally [T ~> s] (could be scaled to accelerate)
  real :: sdt_damp  ! dt for damping [T ~> s] (sdt could be split).
  logical :: use_drag_rate ! Flag to indicate drag_rate is finite
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  real(kind=real32), dimension(size(MEKE%MEKE),NUM_FEATURES) :: features_array ! The array of features
                                        ! needed for the machine learning inference, with different
                                        ! units for the various subarrays [various]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_MEKE: Module must be initialized before it is used.")

  if ((CS%MEKE_Cd_scale > 0.0) .or. (CS%MEKE_Cb>0.) .or. CS%visc_drag) then
    use_drag_rate = .true.
  else
    use_drag_rate = .false.
  endif

  ! Only integrate the MEKE equations if MEKE is required.
  if (.not. allocated(MEKE%MEKE)) then
!   call MOM_error(FATAL, "MOM_MEKE: MEKE%MEKE is not associated!")
    return
  endif

  select case(CS%eke_src)
  case(EKE_PROG)
    if (CS%debug) then
      if (allocated(MEKE%mom_src)) &
        call hchksum(MEKE%mom_src, 'MEKE mom_src', G%HI, unscale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      if (allocated(MEKE%GME_snk)) &
        call hchksum(MEKE%GME_snk, 'MEKE GME_snk', G%HI, unscale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      if (allocated(MEKE%GM_src)) &
        call hchksum(MEKE%GM_src, 'MEKE GM_src', G%HI, unscale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
      if (allocated(MEKE%MEKE)) &
        call hchksum(MEKE%MEKE, 'MEKE MEKE', G%HI, unscale=US%L_T_to_m_s**2)
      call uvchksum("MEKE SN_[uv]", SN_u, SN_v, G%HI, unscale=US%s_to_T, &
                    scalar_pair=.true.)
      call uvchksum("MEKE h[uv]", hu, hv, G%HI, haloshift=0, symmetric=.true., &
                    unscale=GV%H_to_m*US%L_to_m**2)
    endif

    sdt = dt*CS%MEKE_dtScale ! Scaled dt to use for time-stepping
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
          baroHu(I,j) = baroHu(I,j) + hu(I,j,k) * GV%H_to_RZ
        enddo ; enddo
      enddo
      do J=js-1,je ; do i=is,ie
        baroHv(i,J) = 0.
      enddo ; enddo
      do k=1,nz
        do J=js-1,je ; do i=is,ie
          baroHv(i,J) = baroHv(i,J) + hv(i,J,k) * GV%H_to_RZ
        enddo ; enddo
      enddo
      if (CS%MEKE_advection_bug) then
        ! This obviously incorrect code reproduces a bug in the original implementation of
        ! the MEKE advection.
        do j=js,je ; do I=is-1,ie
          baroHu(I,j) = hu(I,j,nz) * GV%H_to_RZ
        enddo ; enddo
        do J=js-1,je ; do i=is,ie
          baroHv(i,J) = hv(i,J,nz) * GV%H_to_RZ
        enddo ; enddo
      endif
    endif

    ! Calculate drag_rate_visc(i,j) which accounts for the model bottom mean flow
    if (CS%visc_drag .and. allocated(visc%Kv_bbl_u) .and. allocated(visc%Kv_bbl_v)) then
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
        drag_rate_visc(i,j) = (0.25*G%IareaT(i,j) * &
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

    if (CS%fixed_total_depth) then
      if (GV%Boussinesq) then
        !$OMP parallel do default(shared)
        do j=js-1,je+1 ; do i=is-1,ie+1
          depth_tot(i,j) = (G%bathyT(i,j) + G%Z_ref) * GV%Z_to_H
        enddo ; enddo
      else
        !$OMP parallel do default(shared)
        do j=js-1,je+1 ; do i=is-1,ie+1
          depth_tot(i,j) = (G%bathyT(i,j) + G%Z_ref) * CS%rho_fixed_total_depth * GV%RZ_to_H
        enddo ; enddo
      endif
    else
      !$OMP parallel do default(shared)
      do j=js-1,je+1 ; do i=is-1,ie+1
        depth_tot(i,j) = mass(i,j) * GV%RZ_to_H
      enddo ; enddo
    endif

    if (CS%initialize) then
      call MEKE_equilibrium(CS, MEKE, G, GV, US, SN_u, SN_v, drag_rate_visc, I_mass, depth_tot)
      CS%initialize = .false.
    endif

    ! Calculates bottomFac2, barotrFac2 and LmixScale
    call MEKE_lengthScales(CS, MEKE, G, GV, US, SN_u, SN_v, MEKE%MEKE, depth_tot, bottomFac2, barotrFac2, LmixScale)
    if (CS%debug) then
      if (CS%visc_drag) &
        call uvchksum("MEKE drag_vel_[uv]", drag_vel_u, drag_vel_v, G%HI, &
                      unscale=GV%H_to_mks*US%s_to_T, scalar_pair=.true.)
      call hchksum(mass, 'MEKE mass',G%HI,haloshift=1, unscale=US%RZ_to_kg_m2)
      call hchksum(drag_rate_visc, 'MEKE drag_rate_visc', G%HI, unscale=GV%H_to_mks*US%s_to_T)
      call hchksum(bottomFac2, 'MEKE bottomFac2', G%HI)
      call hchksum(barotrFac2, 'MEKE barotrFac2', G%HI)
      call hchksum(LmixScale, 'MEKE LmixScale', G%HI, unscale=US%L_to_m)
    endif

    ! Aggregate sources of MEKE (background, frictional and GM)
    !$OMP parallel do default(shared)
    do j=js,je ; do i=is,ie
      src(i,j) = CS%MEKE_BGsrc
    enddo ; enddo

    if (allocated(MEKE%mom_src)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_FrCoeff*I_mass(i,j)*MEKE%mom_src(i,j)
      enddo ; enddo
    endif

    if (allocated(MEKE%GME_snk)) then
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_GMECoeff*I_mass(i,j)*MEKE%GME_snk(i,j)
      enddo ; enddo
    endif

    if (allocated(MEKE%GM_src)) then
      if (CS%GM_src_alt) then
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src(i,j) = src(i,j) - CS%MEKE_GMcoeff*MEKE%GM_src(i,j) / &
                     (GV%H_to_RZ * MAX(CS%MEKE_min_depth_tot, depth_tot(i,j)))
        enddo ; enddo
      else
        !$OMP parallel do default(shared)
        do j=js,je ; do i=is,ie
          src(i,j) = src(i,j) - CS%MEKE_GMcoeff*I_mass(i,j)*MEKE%GM_src(i,j)
        enddo ; enddo
      endif
    endif

    if (CS%MEKE_equilibrium_restoring) then
      call MEKE_equilibrium_restoring(CS, G, GV, US, SN_u, SN_v, depth_tot, &
                                      equilibrium_value)
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_restoring_rate*(MEKE%MEKE(i,j) - equilibrium_value(i,j))
      enddo ; enddo
    endif

    if (CS%debug) then
      call hchksum(src, "MEKE src", G%HI, haloshift=0, unscale=US%L_to_m**2*US%s_to_T**3)
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
        drag_rate(i,j) = (GV%H_to_RZ * I_mass(i,j)) * sqrt( drag_rate_visc(i,j)**2 + &
                 cdrag2 * ( max(0.0, 2.0*bottomFac2(i,j)*MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2 ) )
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=js,je ; do i=is,ie
        drag_rate(i,j) = 0.
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
        MEKE_uflux(I,j) = ((G%dy_Cu(I,j)*G%IdxCu(I,j)) * G%OBCmaskCu(I,j)) * &
            (MEKE%MEKE(i+1,j) - MEKE%MEKE(i,j))
      ! This would have units of [R Z L2 T-2 ~> kg s-2]
      ! MEKE_uflux(I,j) = ((G%dy_Cu(I,j)*G%IdxCu(I,j)) * &
      !     ((2.0*mass(i,j)*mass(i+1,j)) / ((mass(i,j)+mass(i+1,j)) + mass_neglect)) ) * &
      !     (MEKE%MEKE(i+1,j) - MEKE%MEKE(i,j))
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=js-2,je+1 ; do i=is-1,ie+1
        ! MEKE_vflux is used here as workspace with units of [L2 T-2 ~> m2 s-2].
        MEKE_vflux(i,J) = ((G%dx_Cv(i,J)*G%IdyCv(i,J)) * G%OBCmaskCv(i,J)) * &
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
        if (allocated(MEKE%Kh)) &
          Kh_here = max(0., CS%MEKE_Kh) + &
              CS%KhMEKE_Fac*0.5*(MEKE%Kh(i,j)+MEKE%Kh(i+1,j))
        if (allocated(MEKE%Kh_diff)) &
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
        if (allocated(MEKE%Kh)) &
          Kh_here = max(0.,CS%MEKE_Kh) + CS%KhMEKE_Fac * 0.5*(MEKE%Kh(i,j)+MEKE%Kh(i,j+1))
        if (allocated(MEKE%Kh_diff)) &
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
            drag_rate(i,j) = (GV%H_to_RZ * I_mass(i,j)) * sqrt( drag_rate_visc(i,j)**2 + &
                   cdrag2 * ( max(0.0, 2.0*bottomFac2(i,j)*MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2 ) )
          enddo ; enddo
        endif
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
    endif ! MEKE_KH>=0

    if (CS%debug) then
      call hchksum(MEKE%MEKE, "MEKE post-update MEKE", G%HI, haloshift=0, unscale=US%L_T_to_m_s**2)
    endif

  case(EKE_FILE)
    call time_interp_external(CS%eke_handle, Time, data_eke, scale=US%m_s_to_L_T**2)
    do j=js,je ; do i=is,ie
      MEKE%MEKE(i,j) = data_eke(i,j) * G%mask2dT(i,j)
    enddo; enddo
    call MEKE_lengthScales(CS, MEKE, G, GV, US, SN_u, SN_v, MEKE%MEKE, depth_tot, bottomFac2, barotrFac2, LmixScale)
  case(EKE_DBCLIENT)
    call pass_vector(u, v, G%Domain)
    call MEKE_lengthScales(CS, MEKE, G, GV, US, SN_u, SN_v, MEKE%MEKE, depth_tot, bottomFac2, barotrFac2, LmixScale)
    call ML_MEKE_calculate_features(G, GV, US, CS, MEKE%Rd_dx_h, u, v, tv, h, dt, features_array)
    call predict_MEKE(G, CS, SIZE(h), Time, features_array, MEKE%MEKE)
  case default
    call MOM_error(FATAL,"Invalid method specified for calculating EKE")
  end select

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

  if (allocated(MEKE%Kh) .or. allocated(MEKE%Ku) .or. allocated(MEKE%Au)) then
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

!> Calculates the equilibrium solution where the source depends only on MEKE diffusivity
!! and there is no lateral diffusion of MEKE.
!! Results is in MEKE%MEKE.
subroutine MEKE_equilibrium(CS, MEKE, G, GV, US, SN_u, SN_v, drag_rate_visc, I_mass, depth_tot)
  type(ocean_grid_type),             intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),           intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type
  type(MEKE_CS),                     intent(in)    :: CS   !< MEKE control structure.
  type(MEKE_type),                   intent(inout) :: MEKE !< MEKE fields
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: drag_rate_visc !< Mean flow velocity contribution
                                                           !! to the MEKE drag rate [H T-1 ~> m s-1 or kg m-2 s-1]
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: I_mass  !< Inverse of column mass [R-1 Z-1 ~> m2 kg-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: depth_tot !< The thickness of the water column [H ~> m or kg m-2].

  ! Local variables
  real :: beta ! Combined topographic and planetary vorticity gradient [T-1 L-1 ~> s-1 m-1]
  real :: SN   ! The local Eady growth rate [T-1 ~> s-1]
  real :: bottomFac2, barotrFac2    ! Vertical structure factors [nondim]
  real :: LmixScale, LRhines, LEady ! Various mixing length scales [L ~> m]
  real :: KhCoeff ! A copy of MEKE_KhCoeff from the control structure [nondim]
  real :: Kh    ! A lateral diffusivity [L2 T-1 ~> m2 s-1]
  real :: Ubg2  ! Background (tidal?) velocity squared [L2 T-2 ~> m2 s-2]
  real :: cd2   ! The square of the drag coefficient times unit conversion factors [H2 L-2 ~> nondim or kg2 m-6]
  real :: drag_rate ! The MEKE spindown timescale due to bottom drag [T-1 ~> s-1].
  real :: src   ! The sum of MEKE sources [L2 T-3 ~> W kg-1]
  real :: ldamping  ! The MEKE damping rate [T-1 ~> s-1].
  real :: EKE, EKEmin, EKEmax, EKEerr ! [L2 T-2 ~> m2 s-2]
  real :: resid, ResMin, ResMax ! Residuals [L2 T-3 ~> W kg-1]
  real :: FatH    ! Coriolis parameter at h points; to compute topographic beta [T-1 ~> s-1]
  real :: beta_topo_x, beta_topo_y    ! Topographic PV gradients in x and y [T-1 L-1 ~> s-1 m-1]
  real :: h_neglect ! A negligible thickness [H ~> m or kg m-2]
  integer :: i, j, is, ie, js, je, n1, n2
  real :: tolerance ! Width of EKE bracket [L2 T-2 ~> m2 s-2].
  logical :: useSecant, debugIteration

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  debugIteration = .false.
  KhCoeff = CS%MEKE_KhCoeff
  Ubg2 = CS%MEKE_Uscale**2
  cd2 = CS%cdrag**2
  tolerance = 1.0e-12*US%m_s_to_L_T**2
  h_neglect = GV%H_subroundoff

!$OMP do
  do j=js,je ; do i=is,ie
    ! SN = 0.25*max( (SN_u(I,j) + SN_u(I-1,j)) + (SN_v(i,J) + SN_v(i,J-1)), 0.)
    ! This avoids extremes values in equilibrium solution due to bad values in SN_u, SN_v
    SN = min(SN_u(I,j), SN_u(I-1,j), SN_v(i,J), SN_v(i,J-1))

    if (CS%MEKE_equilibrium_alt) then
      MEKE%MEKE(i,j) = (CS%MEKE_GEOMETRIC_alpha * SN * depth_tot(i,j))**2 / cd2
    else
      FatH = 0.25*((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1)) + &
                   (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J-1))) ! Coriolis parameter at h points

      ! Since zero-bathymetry cells are masked, this avoids calculations on land
      if (CS%MEKE_topographic_beta == 0. .or. (depth_tot(i,j) == 0.0)) then
        beta_topo_x = 0. ; beta_topo_y = 0.
      else
        !### Consider different combinations of these estimates of topographic beta.
        beta_topo_x = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (depth_tot(i+1,j)-depth_tot(i,j)) * G%IdxCu(I,j)  &
                  / max(depth_tot(i+1,j), depth_tot(i,j), h_neglect) &
              +       (depth_tot(i,j)-depth_tot(i-1,j)) * G%IdxCu(I-1,j) &
                  / max(depth_tot(i,j), depth_tot(i-1,j), h_neglect) )
        beta_topo_y = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (depth_tot(i,j+1)-depth_tot(i,j)) * G%IdyCv(i,J)  &
                  / max(depth_tot(i,j+1), depth_tot(i,j), h_neglect) + &
                      (depth_tot(i,j)-depth_tot(i,j-1)) * G%IdyCv(i,J-1) &
                  / max(depth_tot(i,j), depth_tot(i,j-1), h_neglect) )
      endif
      beta =  sqrt((G%dF_dx(i,j) + beta_topo_x)**2 + &
                   (G%dF_dy(i,j) + beta_topo_y)**2 )

      if (KhCoeff*SN*I_mass(i,j)>0.) then
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
          call MEKE_lengthScales_0d(CS, US, G%areaT(i,j), beta, depth_tot(i,j), &
                                    MEKE%Rd_dx_h(i,j), SN, EKE, &
                                    bottomFac2, barotrFac2, LmixScale, LRhines, LEady)
          ! TODO: Should include resolution function in Kh
          Kh = (KhCoeff * sqrt(2.*barotrFac2*EKE) * LmixScale)
          src = Kh * (SN * SN)
          drag_rate = (GV%H_to_RZ * I_mass(i,j)) * sqrt(drag_rate_visc(i,j)**2 + cd2 * ( 2.0*bottomFac2*EKE + Ubg2 ) )
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
          drag_rate = (GV%H_to_RZ * I_mass(i,j)) * sqrt( drag_rate_visc(i,j)**2 + cd2 * ( 2.0*bottomFac2*EKE + Ubg2 ) )
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
subroutine MEKE_equilibrium_restoring(CS, G, GV, US, SN_u, SN_v, depth_tot, &
                                      equilibrium_value)
  type(ocean_grid_type),             intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),           intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type.
  type(MEKE_CS),                     intent(in)    :: CS   !< MEKE control structure.
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: depth_tot !< The thickness of the water column [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: equilibrium_value
      !< Equilbrium value of MEKE to be calculated at each time step [L2 T-2 ~> m2 s-2]

  ! Local variables
  real :: SN                      ! The local Eady growth rate [T-1 ~> s-1]
  integer :: i, j, is, ie, js, je ! local indices
  real :: cd2                     ! The square of the drag coefficient [nondim]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  cd2 = CS%cdrag**2
  equilibrium_value(:,:) = 0.0

!$OMP do
  do j=js,je ; do i=is,ie
    ! SN = 0.25*max( (SN_u(I,j) + SN_u(I-1,j)) + (SN_v(i,J) + SN_v(i,J-1)), 0.)
    ! This avoids extremes values in equilibrium solution due to bad values in SN_u, SN_v
    SN = min(SN_u(I,j), SN_u(I-1,j), SN_v(i,J), SN_v(i,J-1))
    equilibrium_value(i,j) = (CS%MEKE_GEOMETRIC_alpha * SN * depth_tot(i,j))**2 / cd2
  enddo ; enddo

  if (CS%id_MEKE_equilibrium>0) call post_data(CS%id_MEKE_equilibrium, equilibrium_value, CS%diag)
end subroutine MEKE_equilibrium_restoring

!> Calculates the eddy mixing length scale and \f$\gamma_b\f$ and \f$\gamma_t\f$
!! functions that are ratios of either bottom or barotropic eddy energy to the
!! column eddy energy, respectively.  See \ref section_MEKE_equations.
subroutine MEKE_lengthScales(CS, MEKE, G, GV, US, SN_u, SN_v, EKE, depth_tot, &
                             bottomFac2, barotrFac2, LmixScale)
  type(MEKE_CS),                     intent(in)    :: CS   !< MEKE control structure.
  type(MEKE_type),                   intent(in)    :: MEKE !< MEKE field
  type(ocean_grid_type),             intent(inout) :: G    !< Ocean grid.
  type(verticalGrid_type),           intent(in)    :: GV   !< Ocean vertical grid structure.
  type(unit_scale_type),             intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: SN_u !< Eady growth rate at u-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: SN_v !< Eady growth rate at v-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: EKE  !< Eddy kinetic energy [L2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: depth_tot !< The thickness of the water column [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: bottomFac2 !< gamma_b^2 [nondim]
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: barotrFac2 !< gamma_t^2 [nondim]
  real, dimension(SZI_(G),SZJ_(G)),  intent(out)   :: LmixScale !< Eddy mixing length [L ~> m].
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: LRhines, LEady  ! Possible mixing length scales [L ~> m]
  real :: beta ! Combined topographic and planetary vorticity gradient [T-1 L-1 ~> s-1 m-1]
  real :: SN   ! The local Eady growth rate [T-1 ~> s-1]
  real :: FatH ! Coriolis parameter at h points [T-1 ~> s-1]
  real :: beta_topo_x, beta_topo_y  ! Topographic PV gradients in x and y [T-1 L-1 ~> s-1 m-1]
  real :: h_neglect ! A negligible thickness [H ~> m or kg m-2]
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  h_neglect = GV%H_subroundoff

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

      ! If depth_tot is zero, then a division by zero FPE will be raised.  In this
      ! case, we apply Adcroft's rule of reciprocals and set the term to zero.
      ! Since zero-bathymetry cells are masked, this should not affect values.
      if (CS%MEKE_topographic_beta == 0. .or. (depth_tot(i,j) == 0.0)) then
        beta_topo_x = 0. ; beta_topo_y = 0.
      else
        !### Consider different combinations of these estimates of topographic beta.
        beta_topo_x = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (depth_tot(i+1,j)-depth_tot(i,j)) * G%IdxCu(I,j)  &
                 / max(depth_tot(i+1,j), depth_tot(i,j), h_neglect) &
              +       (depth_tot(i,j)-depth_tot(i-1,j)) * G%IdxCu(I-1,j) &
                 / max(depth_tot(i,j), depth_tot(i-1,j), h_neglect) )
        beta_topo_y = -CS%MEKE_topographic_beta * FatH * 0.5 * ( &
                      (depth_tot(i,j+1)-depth_tot(i,j)) * G%IdyCv(i,J)  &
                 / max(depth_tot(i,j+1), depth_tot(i,j), h_neglect) + &
                      (depth_tot(i,j)-depth_tot(i,j-1)) * G%IdyCv(i,J-1) &
                 / max(depth_tot(i,j), depth_tot(i,j-1), h_neglect) )
      endif
      beta =  sqrt((G%dF_dx(i,j) + beta_topo_x)**2 + &
                   (G%dF_dy(i,j) + beta_topo_y)**2 )

    else
      beta = 0.
    endif
    ! Returns bottomFac2, barotrFac2 and LmixScale
    call MEKE_lengthScales_0d(CS, US, G%areaT(i,j), beta, depth_tot(i,j),  &
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
subroutine MEKE_lengthScales_0d(CS, US, area, beta, depth_tot, Rd_dx, SN, EKE, &
                                bottomFac2, barotrFac2, LmixScale, Lrhines, Leady)
  type(MEKE_CS), intent(in)    :: CS         !< MEKE control structure.
  type(unit_scale_type), intent(in) :: US    !< A dimensional unit scaling type
  real,          intent(in)    :: area       !< Grid cell area [L2 ~> m2]
  real,          intent(in)    :: beta       !< Planetary beta = \f$ \nabla f\f$  [T-1 L-1 ~> s-1 m-1]
  real,          intent(in)    :: depth_tot  !< The total thickness of the water column [H ~> m or kg m-2]
  real,          intent(in)    :: Rd_dx      !< Resolution Ld/dx [nondim].
  real,          intent(in)    :: SN         !< Eady growth rate [T-1 ~> s-1].
  real,          intent(in)    :: EKE        !< Eddy kinetic energy [L2 T-2 ~> m2 s-2].
  real,          intent(out)   :: bottomFac2 !< gamma_b^2 [nondim]
  real,          intent(out)   :: barotrFac2 !< gamma_t^2 [nondim]
  real,          intent(out)   :: LmixScale  !< Eddy mixing length [L ~> m].
  real,          intent(out)   :: Lrhines    !< Rhines length scale [L ~> m].
  real,          intent(out)   :: Leady      !< Eady length scale [L ~> m].
  ! Local variables
  real :: Lgrid, Ldeform, Lfrict ! Length scales [L ~> m]
  real :: Ue  ! An eddy velocity [L T-1 ~> m s-1]

  ! Length scale for MEKE derived diffusivity
  Lgrid = sqrt(area)               ! Grid scale
  Ldeform = Lgrid * Rd_dx          ! Deformation scale
  Lfrict = depth_tot / CS%cdrag    ! Frictional arrest scale
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
      LmixScale = CS%lscale_maxval
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
logical function MEKE_init(Time, G, GV, US, param_file, diag, dbcomms_CS, CS, MEKE, restart_CS, meke_in_dynamics)
  type(time_type),         intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),   intent(inout) :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Ocean vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file parser structure.
  type(dbcomms_CS_type),   intent(in)    :: dbcomms_CS !< Database communications control structure
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(MEKE_CS),           intent(inout) :: CS         !< MEKE control structure.
  type(MEKE_type),         intent(inout) :: MEKE       !< MEKE fields
  type(MOM_restart_CS),    intent(in)    :: restart_CS !< MOM restart control structure
  logical,                 intent(  out) :: meke_in_dynamics !< If true, MEKE is stepped forward in dynamics
                                                             !! otherwise in tracer dynamics

  ! Local variables
  real :: MEKE_restoring_timescale ! The timescale used to nudge MEKE toward its equilibrium value [T ~> s]
  real :: cdrag            ! The default bottom drag coefficient [nondim].
  character(len=200) :: eke_filename, eke_varname, inputdir
  character(len=16) :: eke_source_str
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  logical :: laplacian, biharmonic, coldStart
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
  CS%initialized = .true.
  call get_param(param_file, mdl, "MEKE_IN_DYNAMICS", meke_in_dynamics, &
                 "If true, step MEKE forward with the dynamics"// &
                 "otherwise with the tracer timestep.", &
                 default=.true.)

  call get_param(param_file, mdl, "EKE_SOURCE", eke_source_str, &
                 "Determine the where EKE comes from:\n" // &
                 "  'prog': Calculated solving EKE equation\n"// &
                 "  'file': Read in from a file\n"            // &
                 "  'dbclient': Retrieved from ML-database", default='prog')

  call MOM_mesg("MEKE_init: reading parameters ", 5)

  select case (lowercase(eke_source_str))
  case("file")
    CS%eke_src = EKE_FILE
    call time_interp_external_init
    call get_param(param_file, mdl, "EKE_FILE", eke_filename, &
                 "A file in which to find the eddy kineteic energy variable.", &
                 default="eke_file.nc")
    call get_param(param_file, mdl, "EKE_VARIABLE", eke_varname, &
                 "The name of the eddy kinetic energy variable to read from "//&
                 "EKE_FILE to use in MEKE.", &
                 default="eke")
    call get_param(param_file, mdl, "INPUTDIR", inputdir, &
                 "The directory in which all input files are found.", &
                 default=".", do_not_log=.true.)
    inputdir = slasher(inputdir)

    eke_filename = trim(inputdir) // trim(eke_filename)
    CS%eke_handle = init_external_field(eke_filename, eke_varname, domain=G%Domain%mpp_domain)
  case("prog")
    CS%eke_src = EKE_PROG
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
                     "The timescale used to nudge MEKE toward its equilibrium value.", &
                     units="s", default=1e6, scale=US%s_to_T)
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
                   "A background energy source for MEKE.", &
                   units="W kg-1", default=0.0, scale=US%m_to_L**2*US%T_to_s**3)
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
  case("dbclient")
    CS%eke_src = EKE_DBCLIENT
    call ML_MEKE_init(diag, G, US, Time, param_file, dbcomms_CS, CS)
  case default
    call MOM_error(FATAL, "Invalid method selected for calculating EKE")
  end select
  ! GMM, make sure all parameters used to calculated MEKE are within the above if

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
  call get_param(param_file, mdl, "MEKE_MIN_DEPTH_TOT", CS%MEKE_min_depth_tot, &
                 "The minimum total depth over which to distribute MEKE energy sources.  "//&
                 "When the total depth is less than this, the sources are scaled away.", &
                 units="m", default=1.0, scale=GV%m_to_H, do_not_log=.not.CS%GM_src_alt)
  call get_param(param_file, mdl, "MEKE_VISC_DRAG", CS%visc_drag, &
                 "If true, use the vertvisc_type to calculate the bottom "//&
                 "drag acting on MEKE.", default=.true.)
  call get_param(param_file, mdl, "MEKE_KHTH_FAC", MEKE%KhTh_fac, &
                 "A factor that maps MEKE%Kh to KhTh.", units="nondim", default=0.0)
  call get_param(param_file, mdl, "MEKE_KHTR_FAC", MEKE%KhTr_fac, &
                 "A factor that maps MEKE%Kh to KhTr.", units="nondim", default=0.0)
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
  call get_param(param_file, mdl, "MEKE_LSCALE_MAX_VAL", CS%lscale_maxval, &
                 "The ceiling on the value of the MEKE length scale when MEKE_MIN_LSCALE=True.  "//&
                 "The default is the distance from the equator to the pole on Earth, as "//&
                 "estimated by enlightenment era scientists, but should probably scale with RAD_EARTH.", &
                 units="m", default=1.0e7, scale=US%m_to_L, do_not_log=.not.CS%use_min_lscale)
  call get_param(param_file, mdl, "MEKE_RD_MAX_SCALE", CS%Rd_as_max_scale, &
                 "If true, the length scale used by MEKE is the minimum of "//&
                 "the deformation radius or grid-spacing. Only used if "//&
                 "MEKE_OLD_LSCALE=True", default=.false.)
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
  call get_param(param_file, mdl, "MEKE_FIXED_TOTAL_DEPTH", CS%fixed_total_depth, &
                 "If true, use the nominal bathymetric depth as the estimate of the "//&
                 "time-varying ocean depth.  Otherwise base the depth on the total ocean mass"//&
                 "per unit area.", default=.true.)
  call get_param(param_file, mdl, "MEKE_TOTAL_DEPTH_RHO", CS%rho_fixed_total_depth, &
                 "A density used to translate the nominal bathymetric depth into an estimate "//&
                 "of the total ocean mass per unit area when MEKE_FIXED_TOTAL_DEPTH is true.", &
                 units="kg m-3", default=GV%Rho0*US%R_to_kg_m3, scale=US%kg_m3_to_R, &
                 do_not_log=(GV%Boussinesq.or.(.not.CS%fixed_total_depth)))

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
  call get_param(param_file, mdl, "MEKE_ADVECTION_BUG", CS%MEKE_advection_bug, &
                 "If true, recover a bug in the calculation of the barotropic transport for "//&
                 "the advection of MEKE.  With the bug, only the transports in the deepest "//&
                 "layer are used.", default=.false., do_not_log=(CS%MEKE_advection_factor<=0.))
  call get_param(param_file, mdl, "MEKE_TOPOGRAPHIC_BETA", CS%MEKE_topographic_beta, &
                 "A scale factor to determine how much topographic beta is weighed in " //&
                 "computing beta in the expression of Rhines scale. Use 1 if full "//&
                 "topographic beta effect is considered; use 0 if it's completely ignored.", &
                 units="nondim", default=0.0)

  ! Nonlocal module parameters
  call get_param(param_file, mdl, "CDRAG", cdrag, &
                 "CDRAG is the drag coefficient relating the magnitude of the velocity "//&
                 "field to the bottom stress.", units="nondim", default=0.003)
  call get_param(param_file, mdl, "MEKE_CDRAG", CS%cdrag, &
                 "Drag coefficient relating the magnitude of the velocity "//&
                 "field to the bottom stress in MEKE.", units="nondim", default=cdrag, scale=US%L_to_m*GV%m_to_H)
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
  if (.not. allocated(MEKE%MEKE)) CS%id_MEKE = -1
  CS%id_Kh = register_diag_field('ocean_model', 'MEKE_KH', diag%axesT1, Time, &
     'MEKE derived diffusivity', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  if (.not. allocated(MEKE%Kh)) CS%id_Kh = -1
  CS%id_Ku = register_diag_field('ocean_model', 'MEKE_KU', diag%axesT1, Time, &
     'MEKE derived lateral viscosity', 'm2 s-1', conversion=US%L_to_m**2*US%s_to_T)
  if (.not. allocated(MEKE%Ku)) CS%id_Ku = -1
  CS%id_Au = register_diag_field('ocean_model', 'MEKE_AU', diag%axesT1, Time, &
     'MEKE derived lateral biharmonic viscosity', 'm4 s-1', conversion=US%L_to_m**4*US%s_to_T)
  if (.not. allocated(MEKE%Au)) CS%id_Au = -1
  CS%id_Ue = register_diag_field('ocean_model', 'MEKE_Ue', diag%axesT1, Time, &
     'MEKE derived eddy-velocity scale', 'm s-1', conversion=US%L_T_to_m_s)
  if (.not. allocated(MEKE%MEKE)) CS%id_Ue = -1
  CS%id_Ub = register_diag_field('ocean_model', 'MEKE_Ub', diag%axesT1, Time, &
     'MEKE derived bottom eddy-velocity scale', 'm s-1', conversion=US%L_T_to_m_s)
  if (.not. allocated(MEKE%MEKE)) CS%id_Ub = -1
  CS%id_Ut = register_diag_field('ocean_model', 'MEKE_Ut', diag%axesT1, Time, &
     'MEKE derived barotropic eddy-velocity scale', 'm s-1', conversion=US%L_T_to_m_s)
  if (.not. allocated(MEKE%MEKE)) CS%id_Ut = -1
  CS%id_src = register_diag_field('ocean_model', 'MEKE_src', diag%axesT1, Time, &
     'MEKE energy source', 'm2 s-3', conversion=(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_decay = register_diag_field('ocean_model', 'MEKE_decay', diag%axesT1, Time, &
     'MEKE decay rate', 's-1', conversion=US%s_to_T)
  CS%id_GM_src = register_diag_field('ocean_model', 'MEKE_GM_src', diag%axesT1, Time, &
     'MEKE energy available from thickness mixing', &
     'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  if (.not. allocated(MEKE%GM_src)) CS%id_GM_src = -1
  CS%id_mom_src = register_diag_field('ocean_model', 'MEKE_mom_src',diag%axesT1, Time, &
     'MEKE energy available from momentum', &
     'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  if (.not. allocated(MEKE%mom_src)) CS%id_mom_src = -1
  CS%id_GME_snk = register_diag_field('ocean_model', 'MEKE_GME_snk',diag%axesT1, Time, &
     'MEKE energy lost to GME backscatter', &
     'W m-2', conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  if (.not. allocated(MEKE%GME_snk)) CS%id_GME_snk = -1
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

  ! Set up group passes.  In the case of a restart, these fields need a halo update now.
  if (allocated(MEKE%MEKE)) then
    call create_group_pass(CS%pass_MEKE, MEKE%MEKE, G%Domain)
    if (allocated(MEKE%Kh_diff)) call create_group_pass(CS%pass_MEKE, MEKE%Kh_diff, G%Domain)
    if (.not.CS%initialize) call do_group_pass(CS%pass_MEKE, G%Domain)
  endif
  if (allocated(MEKE%Kh)) call create_group_pass(CS%pass_Kh, MEKE%Kh, G%Domain)
  if (allocated(MEKE%Ku)) call create_group_pass(CS%pass_Kh, MEKE%Ku, G%Domain)
  if (allocated(MEKE%Au)) call create_group_pass(CS%pass_Kh, MEKE%Au, G%Domain)

  if (allocated(MEKE%Kh) .or. allocated(MEKE%Ku) .or. allocated(MEKE%Au)) &
    call do_group_pass(CS%pass_Kh, G%Domain)

end function MEKE_init

!> Initializer for the variant of MEKE that uses ML to predict eddy kinetic energy
subroutine ML_MEKE_init(diag, G, US, Time, param_file, dbcomms_CS, CS)
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(ocean_grid_type),         intent(inout) :: G           !< The ocean's grid structure.
  type(unit_scale_type),         intent(in)    :: US          !< A dimensional unit scaling type
  type(time_type),               intent(in)    :: Time        !< The current model time.
  type(param_file_type),         intent(in)    :: param_file  !< Parameter file parser structure.
  type(dbcomms_CS_type),         intent(in)    :: dbcomms_CS  !< Control structure for database communication
  type(MEKE_CS),                 intent(inout) :: CS          !< Control structure for this module

  character(len=200)  :: inputdir, backend, model_filename
  integer :: db_return_code, batch_size
  character(len=40) :: mdl = "MOM_ML_MEKE"

  ! Store pointers in control structure
  write(CS%key_suffix, '(A,I6.6)') '_', PE_here()
  ! Put some basic information into the database
  db_return_code = 0
  db_return_code = CS%client%put_tensor("meta"//CS%key_suffix, &
    REAL([G%isd_global, G%idg_offset, G%jsd_global, G%jdg_offset]),[4]) + db_return_code
  db_return_code = CS%client%put_tensor("geolat"//CS%key_suffix, G%geoLatT, shape(G%geoLatT)) + db_return_code
  db_return_code = CS%client%put_tensor("geolon"//CS%key_suffix, G%geoLonT, shape(G%geoLonT)) + db_return_code
  db_return_code = CS%client%put_tensor("EKE_shape"//CS%key_suffix, shape(G%geolonT), [2]) + db_return_code

  if (CS%client%SR_error_parser(db_return_code)) call MOM_error(FATAL, "Putting metadata into the database failed")

  call read_param(param_file, "INPUTDIR", inputdir)
  inputdir = slasher(inputdir)

  call get_param(param_file, mdl, "BATCH_SIZE", batch_size, "Batch size to use for inference", default=1)
  call get_param(param_file, mdl, "EKE_BACKEND", backend, &
                 "The computational backend to use for EKE inference (CPU or GPU)", default="GPU")
  call get_param(param_file, mdl, "EKE_MODEL", model_filename, &
                 "Filename of the a saved pyTorch model to use", fail_if_missing = .true.)
  call get_param(param_file, mdl, "EKE_MAX", CS%eke_max, &
                 "Maximum value of EKE allowed when inferring EKE", &
                 units="m2 s-2", default=2., scale=US%L_T_to_m_s**2)

  ! Set the machine learning model
  if (dbcomms_CS%colocated) then
    if (modulo(PE_here(),dbcomms_CS%colocated_stride) == 0) then
      db_return_code = CS%client%set_model_from_file(CS%model_key, trim(inputdir)//trim(model_filename), &
                                                  "TORCH", backend, batch_size=batch_size)
    endif
  else
    if (is_root_pe()) then
      db_return_code = CS%client%set_model_from_file(CS%model_key, trim(inputdir)//trim(model_filename), &
                                                  "TORCH", backend, batch_size=batch_size)
    endif
  endif
  if (CS%client%SR_error_parser(db_return_code)) then
    call MOM_error(FATAL, "MEKE: set_model failed")
  endif

  call get_param(param_file, mdl, "ONLINE_ANALYSIS", CS%online_analysis, &
               "If true, post EKE used in MOM6 to the database for analysis", default=.true.)

  ! Set various clock ids
  CS%id_client_init   = cpu_clock_id('(ML_MEKE client init)', grain=CLOCK_ROUTINE)
  CS%id_put_tensor    = cpu_clock_id('(ML_MEKE put tensor)', grain=CLOCK_ROUTINE)
  CS%id_run_model     = cpu_clock_id('(ML_MEKE run model)', grain=CLOCK_ROUTINE)
  CS%id_unpack_tensor = cpu_clock_id('(ML_MEKE unpack tensor )', grain=CLOCK_ROUTINE)

  ! Diagnostics for ML_MEKE
  CS%id_mke = register_diag_field('ocean_model', 'MEKE_MKE', diag%axesT1, Time, &
     'Surface mean (resolved) kinetic energy used in MEKE', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_slope_z= register_diag_field('ocean_model', 'MEKE_slope_z', diag%axesT1, Time, &
     'Vertically averaged isopyncal slope magnitude used in MEKE', 'nondim', conversion=US%Z_to_L)
  CS%id_slope_x= register_diag_field('ocean_model', 'MEKE_slope_x', diag%axesCui, Time, &
     'Isopycnal slope in the x-direction used in MEKE', 'nondim', conversion=US%Z_to_L)
  CS%id_slope_y= register_diag_field('ocean_model', 'MEKE_slope_y', diag%axesCvi, Time, &
     'Isopycnal slope in the y-direction used in MEKE', 'nondim', conversion=US%Z_to_L)
  CS%id_rv = register_diag_field('ocean_model', 'MEKE_RV', diag%axesT1, Time, &
     'Surface relative vorticity used in MEKE', 's-1', conversion=US%s_to_T)

end subroutine ML_MEKE_init

!> Calculate the various features used for the machine learning prediction
subroutine ML_MEKE_calculate_features(G, GV, US, CS, Rd_dx_h, u, v, tv, h, dt, features_array)
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  type(MEKE_CS),                             intent(in)    :: CS !< Control structure for MEKE
  real, dimension(SZI_(G),SZJ_(G)),          intent(in   ) :: Rd_dx_h !< Rossby radius of deformation over
                                                                 !! the grid length scale [nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v  !< Meridional velocity [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Type containing thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt !< Model(baroclinic) time-step [T ~> s].
  real(kind=real32), dimension(SIZE(h),num_features), intent(  out) :: features_array
                                                                 !< The array of features needed for machine
                                                                 !! learning inference, with different units
                                                                 !! for the various subarrays [various]

  real, dimension(SZI_(G),SZJ_(G)) :: mke      ! Surface kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJ_(G)) :: slope_z  ! Vertically averaged isoneutral slopes [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJB_(G)) :: rv_z   ! Surface relative vorticity [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G)) :: rv_z_t ! Surface relative vorticity interpolated to tracer points [T-1 ~> s-1]

  real, dimension(SZIB_(G),SZJ_(G), SZK_(G)) :: h_u ! Thickness at u point [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJB_(G), SZK_(G)) :: h_v ! Thickness at v point [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1) :: slope_x ! Isoneutral slope at U point [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1) :: slope_y ! Isoneutral slope at V point [Z L-1 ~> nondim]
  real, dimension(SZIB_(G),SZJ_(G)) :: slope_x_vert_avg ! Isoneutral slope at U point [Z L-1 ~> nondim]
  real, dimension(SZI_(G),SZJB_(G)) :: slope_y_vert_avg ! Isoneutral slope at V point [Z L-1 ~> nondim]
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) ::  e ! The interface heights relative to mean sea level [Z ~> m].
  real :: slope_t  ! Slope interpolated to thickness points [Z L-1 ~> nondim]
  real :: u_t, v_t ! u and v interpolated to thickness points [L T-1 ~> m s-1]
  real :: dvdx, dudy ! Components of relative vorticity [T-1 ~> s-1]
  real :: a_e, a_w, a_n, a_s ! Fractional areas of neighboring cells for interpolating velocities [nondim]
  real :: Idenom    ! A normalizing factor in calculating weighted averages of areas [L-2 ~> m-2]
  real :: sum_area  ! A sum of adjacent cell areas [L2 ~> m2]

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! Calculate various features for used to infer eddy kinetic energy
  ! Linear interpolation to estimate thickness at a velocity points
  do k=1,nz; do j=js-1,je+1; do i=is-1,ie+1
    h_u(I,j,k) = 0.5*(h(i,j,k)*G%mask2dT(i,j) + h(i+1,j,k)*G%mask2dT(i+1,j)) + GV%Angstrom_H
    h_v(i,J,k) = 0.5*(h(i,j,k)*G%mask2dT(i,j) + h(i,j+1,k)*G%mask2dT(i,j+1)) + GV%Angstrom_H
  enddo; enddo; enddo;
  call find_eta(h, tv, G, GV, US, e, halo_size=2)
  ! Note the hard-coded dimenisional constant in the following line.
  call calc_isoneutral_slopes(G, GV, US, h, e, tv, dt*1.e-7*GV%m2_s_to_HZ_T, .false., slope_x, slope_y)
  call pass_vector(slope_x, slope_y, G%Domain)
  do j=js-1,je+1; do i=is-1,ie+1
    slope_x_vert_avg(I,j) = vertical_average_interface(slope_x(i,j,:), h_u(i,j,:), GV%H_subroundoff)
    slope_y_vert_avg(i,J) = vertical_average_interface(slope_y(i,j,:), h_v(i,j,:), GV%H_subroundoff)
  enddo; enddo
  slope_z(:,:) = 0.

  call pass_vector(slope_x_vert_avg, slope_y_vert_avg, G%Domain)
  do j=js,je; do i=is,ie
    ! Calculate weights for interpolation from velocity points to h points
    sum_area = G%areaCu(I-1,j) + G%areaCu(I,j)
    if (sum_area>0.0) then
      Idenom = sqrt(0.5*G%IareaT(i,j) / sum_area)
      a_w = G%areaCu(I-1,j) * Idenom
      a_e = G%areaCu(I,j) * Idenom
    else
      a_w = 0.0 ; a_e = 0.0
    endif

    sum_area = G%areaCv(i,J-1) + G%areaCv(i,J)
    if (sum_area>0.0) then
      Idenom = sqrt(0.5*G%IareaT(i,j) / sum_area)
      a_s = G%areaCv(i,J-1) * Idenom
      a_n = G%areaCv(i,J) * Idenom
    else
      a_s = 0.0 ; a_n = 0.0
    endif

    ! Calculate mean kinetic energy
    u_t = a_e*u(I,j,1)+a_w*u(I-1,j,1)
    v_t = a_n*v(i,J,1)+a_s*v(i,J-1,1)
    mke(i,j) = 0.5*( u_t*u_t + v_t*v_t )

    ! Calculate the magnitude of the slope
    slope_t = slope_x_vert_avg(I,j)*a_e+slope_x_vert_avg(I-1,j)*a_w
    slope_z(i,j) = sqrt(slope_t*slope_t)
    slope_t = slope_y_vert_avg(i,J)*a_n+slope_y_vert_avg(i,J-1)*a_s
    slope_z(i,j) = 0.5*(slope_z(i,j) + sqrt(slope_t*slope_t))*G%mask2dT(i,j)
  enddo; enddo
  call pass_var(slope_z, G%Domain)

  ! Calculate relative vorticity
  do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
    dvdx = (v(i+1,J,1)*G%dyCv(i+1,J) - v(i,J,1)*G%dyCv(i,J))
    dudy = (u(I,j+1,1)*G%dxCu(I,j+1) - u(I,j,1)*G%dxCu(I,j))
    ! Assumed no slip
    rv_z(I,J) = (2.0-G%mask2dBu(I,J)) * (dvdx - dudy) * G%IareaBu(I,J)
  enddo; enddo
  ! Interpolate RV to t-point, revisit this calculation to include metrics
  do j=js,je; do i=is,ie
    rv_z_t(i,j) = 0.25*(rv_z(i-1,j) + rv_z(i,j) + rv_z(i-1,j-1) + rv_z(i,j-1))
  enddo; enddo


  ! Construct the feature array
  features_array(:,mke_idx) = pack(mke,.true.)
  features_array(:,slope_z_idx) = pack(slope_z,.true.)
  features_array(:,rd_dx_z_idx) = pack(Rd_dx_h,.true.)
  features_array(:,rv_idx) = pack(rv_z_t,.true.)

  if (CS%id_rv>0) call post_data(CS%id_rv, rv_z, CS%diag)
  if (CS%id_mke>0) call post_data(CS%id_mke, mke, CS%diag)
  if (CS%id_slope_z>0) call post_data(CS%id_slope_z, slope_z, CS%diag)
  if (CS%id_slope_x>0) call post_data(CS%id_slope_x, slope_x, CS%diag)
  if (CS%id_slope_y>0) call post_data(CS%id_slope_y, slope_y, CS%diag)
end subroutine ML_MEKE_calculate_features

!> Use the machine learning interface to predict EKE
subroutine predict_MEKE(G, CS, npts, Time, features_array, MEKE)
  type(ocean_grid_type),                                 intent(inout) :: G  !< Ocean grid
  type(MEKE_CS),                                         intent(in   ) :: CS !< Control structure for MEKE
  integer,                                               intent(in   ) :: npts !< Number of T-grid cells on the local
                                                                               !! domain
  type(time_type),                                       intent(in   ) :: Time !< The current model time
  real(kind=real32), dimension(npts,num_features),       intent(in   ) :: features_array
                                                                          !< The array of features needed for machine
                                                                          !! learning inference, with different units
                                                                          !! for the various subarrays [various]
  real, dimension(SZI_(G),SZJ_(G)),                      intent(  out) :: MEKE !< Eddy kinetic energy [L2 T-2 ~> m2 s-2]
  integer :: db_return_code
  character(len=255), dimension(1) :: model_out, model_in
  character(len=255) :: time_suffix
  real(kind=real32), dimension(SIZE(MEKE)) :: MEKE_vec ! A one-dimensional array of eddy kinetic
                                                       ! energy [L2 T-2 ~> m2 s-2]

  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
!> Use the database client to call a machine learning model to predict eddy kinetic energy
  call cpu_clock_begin(CS%id_put_tensor)
  db_return_code = CS%client%put_tensor("features"//CS%key_suffix, features_array, shape(features_array))
  call cpu_clock_end(CS%id_put_tensor)

  ! Run the ML model to predict EKE and return the result
  model_out(1) = "EKE"//CS%key_suffix
  model_in(1) = "features"//CS%key_suffix
  call cpu_clock_begin(CS%id_run_model)
  db_return_code = CS%client%run_model(CS%model_key, model_in, model_out)
  call cpu_clock_end(CS%id_run_model)
  if (CS%client%SR_error_parser(db_return_code)) then
    call MOM_error(FATAL, "MEKE: run_model failed")
  endif
  call cpu_clock_begin(CS%id_unpack_tensor)
  db_return_code = CS%client%unpack_tensor( model_out(1), MEKE_vec, shape(MEKE_vec) )
  call cpu_clock_end(CS%id_unpack_tensor)

  !### Does MEKE_vec need to be rescaled from [m2 s-2] to [L2 T-2 ~> m2 s-2] by
  !    multiplying MEKE_vec by US%m_s_to_L_T**2 here?
  MEKE = reshape(MEKE_vec, shape(MEKE))
  do j=js,je; do i=is,ie
    MEKE(i,j) = MIN(MAX(exp(MEKE(i,j)),0.),CS%eke_max)
  enddo; enddo
  call pass_var(MEKE,G%Domain)

  if (CS%online_analysis) then
    write(time_suffix,"(F16.0)") time_type_to_real(Time)
    db_return_code = CS%client%put_tensor(trim("EKE_")//trim(adjustl(time_suffix))//CS%key_suffix, MEKE, shape(MEKE))
  endif
end subroutine predict_MEKE

!> Compute average of interface quantities weighted by the thickness of the surrounding layers
real function vertical_average_interface(h, w, h_min)

  real, dimension(:), intent(in) :: h  !< Layer Thicknesses [H ~> m or kg m-2]
  real, dimension(:), intent(in) :: w  !< Quantity to average [arbitrary]
  real, intent(in) :: h_min !< The vanishingly small layer thickness [H ~> m or kg m-2]

  real :: htot  ! Twice the sum of the layer thicknesses interpolated to interior interfaces [H ~> m or kg m-2]
  real :: inv_htot ! The inverse of htot  [H-1 ~> m-1 or m2 kg-1]
  integer :: k, nk

  nk = size(h)
  htot = h_min
  do k=2,nk
    htot = htot + (h(k-1)+h(k))
  enddo
  inv_htot = 1./htot

  vertical_average_interface = 0.
  do K=2,nk
    vertical_average_interface = vertical_average_interface + (w(k)*(h(k-1)+h(k)))*inv_htot
  enddo
end function vertical_average_interface

!> Allocates memory and register restart fields for the MOM_MEKE module.
subroutine MEKE_alloc_register_restart(HI, US, param_file, MEKE, restart_CS)
! Arguments
  type(hor_index_type),  intent(in)    :: HI         !< Horizontal index structure
  type(unit_scale_type), intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type), intent(in)    :: param_file !< Parameter file parser structure.
  type(MEKE_type),       intent(inout) :: MEKE       !< MEKE fields
  type(MOM_restart_CS),  intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables
  real :: MEKE_GMcoeff, MEKE_FrCoeff, MEKE_GMECoeff  ! Coefficients for various terms [nondim]
  real :: MEKE_KHCoeff, MEKE_viscCoeff_Ku, MEKE_viscCoeff_Au  ! Coefficients for various terms [nondim]
  logical :: Use_KH_in_MEKE
  logical :: useMEKE
  integer :: isd, ied, jsd, jed

! Determine whether this module will be used
  useMEKE = .false.; call read_param(param_file,"USE_MEKE",useMEKE)

! Read these parameters to determine what should be in the restarts
  MEKE_GMcoeff = -1. ; call read_param(param_file,"MEKE_GMCOEFF",MEKE_GMcoeff)
  MEKE_FrCoeff = -1. ; call read_param(param_file,"MEKE_FRCOEFF",MEKE_FrCoeff)
  MEKE_GMEcoeff = -1. ; call read_param(param_file,"MEKE_GMECOEFF",MEKE_GMEcoeff)
  MEKE_KhCoeff = 1. ; call read_param(param_file,"MEKE_KHCOEFF",MEKE_KhCoeff)
  MEKE_viscCoeff_Ku = 0. ; call read_param(param_file,"MEKE_VISCOSITY_COEFF_KU",MEKE_viscCoeff_Ku)
  MEKE_viscCoeff_Au = 0. ; call read_param(param_file,"MEKE_VISCOSITY_COEFF_AU",MEKE_viscCoeff_Au)
  Use_KH_in_MEKE = .false. ; call read_param(param_file,"USE_KH_IN_MEKE", Use_KH_in_MEKE)

  if (.not. useMEKE) return

! Allocate memory
  call MOM_mesg("MEKE_alloc_register_restart: allocating and registering", 5)
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed
  allocate(MEKE%MEKE(isd:ied,jsd:jed), source=0.0)
  call register_restart_field(MEKE%MEKE, "MEKE", .false., restart_CS, &
           longname="Mesoscale Eddy Kinetic Energy", units="m2 s-2", conversion=US%L_T_to_m_s**2)

  if (MEKE_GMcoeff>=0.) allocate(MEKE%GM_src(isd:ied,jsd:jed), source=0.0)
  if (MEKE_FrCoeff>=0. .or. MEKE_GMECoeff>=0.) &
    allocate(MEKE%mom_src(isd:ied,jsd:jed), source=0.0)
  if (MEKE_GMECoeff>=0.) allocate(MEKE%GME_snk(isd:ied,jsd:jed), source=0.0)
  if (MEKE_KhCoeff>=0.) then
    allocate(MEKE%Kh(isd:ied,jsd:jed), source=0.0)
    call register_restart_field(MEKE%Kh, "MEKE_Kh", .false., restart_CS, &
             longname="Lateral diffusivity from Mesoscale Eddy Kinetic Energy", &
             units="m2 s-1", conversion=US%L_to_m**2*US%s_to_T)
  endif
  allocate(MEKE%Rd_dx_h(isd:ied,jsd:jed), source=0.0)
  if (MEKE_viscCoeff_Ku/=0.) then
    allocate(MEKE%Ku(isd:ied,jsd:jed), source=0.0)
    call register_restart_field(MEKE%Ku, "MEKE_Ku", .false., restart_CS, &
             longname="Lateral viscosity from Mesoscale Eddy Kinetic Energy", &
             units="m2 s-1", conversion=US%L_to_m**2*US%s_to_T)
  endif
  if (Use_Kh_in_MEKE) then
    allocate(MEKE%Kh_diff(isd:ied,jsd:jed), source=0.0)
    call register_restart_field(MEKE%Kh_diff, "MEKE_Kh_diff", .false., restart_CS, &
             longname="Copy of thickness diffusivity for diffusing MEKE", &
             units="m2 s-1", conversion=US%L_to_m**2*US%s_to_T)
  endif

  if (MEKE_viscCoeff_Au/=0.) then
    allocate(MEKE%Au(isd:ied,jsd:jed), source=0.0)
    call register_restart_field(MEKE%Au, "MEKE_Au", .false., restart_CS, &
             longname="Lateral biharmonic viscosity from Mesoscale Eddy Kinetic Energy", &
             units="m4 s-1", conversion=US%L_to_m**4*US%s_to_T)
  endif

end subroutine MEKE_alloc_register_restart

!> Deallocates any variables allocated in MEKE_alloc_register_restart.
subroutine MEKE_end(MEKE)
  type(MEKE_type), intent(inout) :: MEKE !< A structure with MEKE-related fields.

  ! NOTE: MEKE will always be allocated by MEKE_init, even if MEKE is disabled.
  !  So these must all be conditional, even though MEKE%MEKE and MEKE%Rd_dx_h
  !  are always allocated (when MEKE is enabled)

  if (allocated(MEKE%Au)) deallocate(MEKE%Au)
  if (allocated(MEKE%Kh_diff)) deallocate(MEKE%Kh_diff)
  if (allocated(MEKE%Ku)) deallocate(MEKE%Ku)
  if (allocated(MEKE%Rd_dx_h)) deallocate(MEKE%Rd_dx_h)
  if (allocated(MEKE%Kh)) deallocate(MEKE%Kh)
  if (allocated(MEKE%GME_snk)) deallocate(MEKE%GME_snk)
  if (allocated(MEKE%mom_src)) deallocate(MEKE%mom_src)
  if (allocated(MEKE%GM_src)) deallocate(MEKE%GM_src)
  if (allocated(MEKE%MEKE)) deallocate(MEKE%MEKE)
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
!! eddy velocity to the bottom. In other words, the bottom velocity is
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

