!> This routine drives the diabatic/dianeutral physics for MOM
module MOM_diabatic_driver

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_bulk_mixed_layer,    only : bulkmixedlayer, bulkmixedlayer_init, bulkmixedlayer_CS
use MOM_debugging,           only : hchksum
use MOM_checksum_packages,   only : MOM_state_chksum, MOM_state_stats
use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_CVMix_shear,         only : CVMix_shear_is_used
use MOM_CVMix_ddiff,         only : CVMix_ddiff_is_used
use MOM_diabatic_aux,        only : diabatic_aux_init, diabatic_aux_end, diabatic_aux_CS
use MOM_diabatic_aux,        only : make_frazil, adjust_salt, differential_diffuse_T_S, triDiagTS
use MOM_diabatic_aux,        only : triDiagTS_Eulerian, find_uv_at_h, diagnoseMLDbyDensityDifference
use MOM_diabatic_aux,        only : applyBoundaryFluxesInOut, diagnoseMLDbyEnergy, set_pen_shortwave
use MOM_diag_mediator,       only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,       only : diag_ctrl, time_type, diag_update_remap_grids
use MOM_diag_mediator,       only : diag_ctrl, query_averaging_enabled, enable_averages, disable_averaging
use MOM_diag_mediator,       only : diag_grid_storage, diag_grid_storage_init, diag_grid_storage_end
use MOM_diag_mediator,       only : diag_copy_diag_to_storage, diag_copy_storage_to_diag
use MOM_diag_mediator,       only : diag_save_grids, diag_restore_grids
use MOM_diapyc_energy_req,   only : diapyc_energy_req_init, diapyc_energy_req_end
use MOM_diapyc_energy_req,   only : diapyc_energy_req_calc, diapyc_energy_req_test, diapyc_energy_req_CS
use MOM_CVMix_conv,          only : CVMix_conv_init, CVMix_conv_cs
use MOM_CVMix_conv,          only : CVMix_conv_end, calculate_CVMix_conv
use MOM_domains,             only : pass_var, To_West, To_South, To_All, Omit_Corners
use MOM_domains,             only : create_group_pass, do_group_pass, group_pass_type
use MOM_energetic_PBL,       only : energetic_PBL, energetic_PBL_init
use MOM_energetic_PBL,       only : energetic_PBL_end, energetic_PBL_CS
use MOM_energetic_PBL,       only : energetic_PBL_get_MLD
use MOM_entrain_diffusive,   only : entrainment_diffusive, entrain_diffusive_init
use MOM_entrain_diffusive,   only : entrain_diffusive_end, entrain_diffusive_CS
use MOM_EOS,                 only : calculate_density, calculate_TFreeze, EOS_domain
use MOM_error_handler,       only : MOM_error, FATAL, WARNING, callTree_showQuery,MOM_mesg
use MOM_error_handler,       only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,         only : get_param, log_version, param_file_type, read_param
use MOM_forcing_type,        only : forcing, MOM_forcing_chksum
use MOM_forcing_type,        only : calculateBuoyancyFlux2d, forcing_SinglePointPrint
use MOM_geothermal,          only : geothermal_entraining, geothermal_in_place
use MOM_geothermal,          only : geothermal_init, geothermal_end, geothermal_CS
use MOM_grid,                only : ocean_grid_type
use MOM_int_tide_input,      only : set_int_tide_input, int_tide_input_init
use MOM_int_tide_input,      only : int_tide_input_end, int_tide_input_CS, int_tide_input_type
use MOM_interface_heights,   only : find_eta
use MOM_internal_tides,      only : propagate_int_tide
use MOM_internal_tides,      only : internal_tides_init, internal_tides_end, int_tide_CS
use MOM_kappa_shear,         only : kappa_shear_is_used
use MOM_CVMix_KPP,           only : KPP_CS, KPP_init, KPP_compute_BLD, KPP_calculate
use MOM_CVMix_KPP,           only : KPP_end, KPP_get_BLD
use MOM_CVMix_KPP,           only : KPP_NonLocalTransport_temp, KPP_NonLocalTransport_saln
use MOM_opacity,             only : opacity_init, opacity_end, opacity_CS
use MOM_opacity,             only : absorbRemainingSW, optics_type, optics_nbands
use MOM_open_boundary,       only : ocean_OBC_type
use MOM_regularize_layers,   only : regularize_layers, regularize_layers_init, regularize_layers_CS
use MOM_set_diffusivity,     only : set_diffusivity, set_BBL_TKE
use MOM_set_diffusivity,     only : set_diffusivity_init, set_diffusivity_end
use MOM_set_diffusivity,     only : set_diffusivity_CS
use MOM_sponge,              only : apply_sponge, sponge_CS
use MOM_ALE_sponge,          only : apply_ALE_sponge, ALE_sponge_CS
use MOM_time_manager,        only : time_type, real_to_time, operator(-), operator(<=)
use MOM_tracer_flow_control, only : call_tracer_column_fns, tracer_flow_control_CS
use MOM_tracer_diabatic,     only : tracer_vertdiff, tracer_vertdiff_Eulerian
use MOM_unit_scaling,        only : unit_scale_type
use MOM_variables,           only : thermo_var_ptrs, vertvisc_type, accel_diag_ptrs
use MOM_variables,           only : cont_diag_ptrs, MOM_thermovar_chksum, p3d
use MOM_verticalGrid,        only : verticalGrid_type, get_thickness_units
use MOM_wave_speed,          only : wave_speeds
use MOM_wave_interface,      only : wave_parameters_CS

implicit none ; private

#include <MOM_memory.h>

public diabatic
public diabatic_driver_init
public diabatic_driver_end
public extract_diabatic_member
public adiabatic
public adiabatic_driver_init
! public legacy_diabatic

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for this module
type, public:: diabatic_CS; private

  logical :: use_legacy_diabatic     !< If true (default), use the a legacy version of the diabatic
                                     !! algorithm. This is temporary and is needed to avoid change
                                     !! in answers.
  logical :: bulkmixedlayer          !< If true, a refined bulk mixed layer is used with
                                     !! nkml sublayers (and additional buffer layers).
  logical :: use_energetic_PBL       !< If true, use the implicit energetics planetary
                                     !! boundary layer scheme to determine the diffusivity
                                     !! in the surface boundary layer.
  logical :: use_KPP                 !< If true, use CVMix/KPP boundary layer scheme to determine the
                                     !! OBLD and the diffusivities within this layer.
  logical :: use_kappa_shear         !< If true, use the kappa_shear module to find the
                                     !! shear-driven diapycnal diffusivity.
  logical :: use_CVMix_shear         !< If true, use the CVMix module to find the
                                     !! shear-driven diapycnal diffusivity.
  logical :: use_CVMix_ddiff         !< If true, use the CVMix double diffusion module.
  logical :: use_CVMix_conv          !< If true, use the CVMix module to get enhanced
                                     !! mixing due to convection.
  logical :: double_diffuse          !< If true, some form of double-diffusive mixing is used.
  logical :: use_sponge              !< If true, sponges may be applied anywhere in the
                                     !! domain.  The exact location and properties of
                                     !! those sponges are set by calls to
                                     !! initialize_sponge and set_up_sponge_field.
  logical :: use_geothermal          !< If true, apply geothermal heating.
  logical :: use_int_tides           !< If true, use the code that advances a separate set
                                     !! of equations for the internal tide energy density.
  logical :: ePBL_is_additive        !< If true, the diffusivity from ePBL is added to all
                                     !! other diffusivities. Otherwise, the larger of kappa-
                                     !! shear and ePBL diffusivities are used.
  real    :: ePBL_Prandtl            !< The Prandtl number used by ePBL to convert vertical
                                     !! diffusivities into viscosities.
  integer :: nMode = 1               !< Number of baroclinic modes to consider
  real    :: uniform_test_cg         !< Uniform group velocity of internal tide
                                     !! for testing internal tides [L T-1 ~> m s-1]
  logical :: useALEalgorithm         !< If true, use the ALE algorithm rather than layered
                                     !! isopycnal/stacked shallow water mode. This logical
                                     !! passed by argument to diabatic_driver_init.
  logical :: aggregate_FW_forcing    !< Determines whether net incoming/outgoing surface
                                     !! FW fluxes are applied separately or combined before
                                     !! being applied.
  real    :: ML_mix_first            !< The nondimensional fraction of the mixed layer
                                     !! algorithm that is applied before diffusive mixing.
                                     !! The default is 0, while 0.5 gives Strang splitting
                                     !! and 1 is a sensible value too.  Note that if there
                                     !! are convective instabilities in the initial state,
                                     !! the first call may do much more than the second.
  integer :: NKBL                    !< The number of buffer layers (if bulk_mixed_layer)
  logical :: massless_match_targets  !< If true (the default), keep the T & S
                                     !! consistent with the target values.
  logical :: mix_boundary_tracers    !< If true, mix the passive tracers in massless layers at the
                                     !! bottom into the interior as though a diffusivity of
                                     !! Kd_min_tr (see below) were operating.
  logical :: mix_boundary_tracer_ALE !< If true, in ALE mode mix the passive tracers in massless
                                     !! layers at the bottom into the interior as though a
                                     !! diffusivity of Kd_min_tr (see below) were operating.
  real    :: Kd_BBL_tr               !< A bottom boundary layer tracer diffusivity that
                                     !! will allow for explicitly specified bottom fluxes
                                     !! [Z2 T-1 ~> m2 s-1].  The entrainment at the bottom is at
                                     !! least sqrt(Kd_BBL_tr*dt) over the same distance.
  real    :: Kd_min_tr               !< A minimal diffusivity that should always be
                                     !! applied to tracers, especially in massless layers
                                     !! near the bottom [Z2 T-1 ~> m2 s-1].
  real    :: minimum_forcing_depth   !< The smallest depth over which heat and freshwater
                                     !! fluxes are applied [H ~> m or kg m-2].
  real    :: evap_CFL_limit = 0.8    !< The largest fraction of a layer that can be
                                     !! evaporated in one time-step [nondim].
  integer :: halo_TS_diff = 0        !< The temperature, salinity and thickness halo size that
                                     !! must be valid for the diffusivity calculations.
  logical :: useKPP = .false.        !< use CVMix/KPP diffusivities and non-local transport
  logical :: KPPisPassive            !< If true, KPP is in passive mode, not changing answers.
  logical :: debug                   !< If true, write verbose checksums for debugging purposes.
  logical :: debugConservation       !< If true, monitor conservation and extrema.
  logical :: tracer_tridiag          !< If true, use tracer_vertdiff instead of tridiagTS for
                                     !< vertical diffusion of T and S
  logical :: debug_energy_req        !< If true, test the mixing energy requirement code.
  type(diag_ctrl), pointer :: diag   !< structure used to regulate timing of diagnostic output
  real    :: MLDdensityDifference    !< Density difference used to determine MLD_user [R ~> kg m-3]
  real    :: dz_subML_N2             !< The distance over which to calculate a diagnostic of the
                                     !! average stratification at the base of the mixed layer [Z ~> m].
  real    :: MLD_EN_VALS(3)          !< Energy values for energy mixed layer diagnostics

  !>@{ Diagnostic IDs
  integer :: id_cg1      = -1                 ! diag handle for mode-1 speed
  integer, allocatable, dimension(:) :: id_cn ! diag handle for all mode speeds
  integer :: id_ea       = -1, id_eb       = -1 ! used by layer diabatic
  integer :: id_ea_t     = -1, id_eb_t     = -1, id_ea_s   = -1, id_eb_s     = -1
  integer :: id_Kd_heat  = -1, id_Kd_salt  = -1, id_Kd_int = -1, id_Kd_ePBL  = -1
  integer :: id_Tdif     = -1, id_Sdif     = -1, id_Tadv   = -1, id_Sadv     = -1
  ! These are handles to diagnostics related to the mixed layer properties.
  integer :: id_MLD_003 = -1, id_MLD_0125 = -1, id_MLD_user = -1, id_mlotstsq = -1
  integer :: id_MLD_EN1 = -1, id_MLD_EN2  = -1, id_MLD_EN3  = -1, id_subMLN2  = -1

  ! These are handles to diatgnostics that are only available in non-ALE layered mode.
  integer :: id_wd       = -1
  integer :: id_dudt_dia = -1, id_dvdt_dia = -1
  integer :: id_hf_dudt_dia_2d = -1, id_hf_dvdt_dia_2d = -1

  ! diagnostic for fields prior to applying diapycnal physics
  integer :: id_u_predia = -1, id_v_predia = -1, id_h_predia = -1
  integer :: id_T_predia = -1, id_S_predia = -1, id_e_predia = -1

  integer :: id_diabatic_diff_temp_tend     = -1
  integer :: id_diabatic_diff_saln_tend     = -1
  integer :: id_diabatic_diff_heat_tend     = -1
  integer :: id_diabatic_diff_salt_tend     = -1
  integer :: id_diabatic_diff_heat_tend_2d  = -1
  integer :: id_diabatic_diff_salt_tend_2d  = -1
  integer :: id_diabatic_diff_h = -1

  integer :: id_boundary_forcing_h       = -1
  integer :: id_boundary_forcing_h_tendency   = -1
  integer :: id_boundary_forcing_temp_tend    = -1
  integer :: id_boundary_forcing_saln_tend    = -1
  integer :: id_boundary_forcing_heat_tend    = -1
  integer :: id_boundary_forcing_salt_tend    = -1
  integer :: id_boundary_forcing_heat_tend_2d = -1
  integer :: id_boundary_forcing_salt_tend_2d = -1

  integer :: id_frazil_h    = -1
  integer :: id_frazil_temp_tend    = -1
  integer :: id_frazil_heat_tend    = -1
  integer :: id_frazil_heat_tend_2d = -1
  !>@}

  logical :: diabatic_diff_tendency_diag = .false. !< If true calculate diffusive tendency diagnostics
  logical :: boundary_forcing_tendency_diag = .false. !< If true calculate frazil diagnostics
  logical :: frazil_tendency_diag = .false. !< If true calculate frazil tendency diagnostics

  type(diabatic_aux_CS),        pointer :: diabatic_aux_CSp      => NULL() !< Control structure for a child module
  type(entrain_diffusive_CS),   pointer :: entrain_diffusive_CSp => NULL() !< Control structure for a child module
  type(bulkmixedlayer_CS),      pointer :: bulkmixedlayer_CSp    => NULL() !< Control structure for a child module
  type(energetic_PBL_CS),       pointer :: energetic_PBL_CSp     => NULL() !< Control structure for a child module
  type(regularize_layers_CS),   pointer :: regularize_layers_CSp => NULL() !< Control structure for a child module
  type(geothermal_CS),          pointer :: geothermal_CSp        => NULL() !< Control structure for a child module
  type(int_tide_CS),            pointer :: int_tide_CSp          => NULL() !< Control structure for a child module
  type(int_tide_input_CS),      pointer :: int_tide_input_CSp    => NULL() !< Control structure for a child module
  type(int_tide_input_type),    pointer :: int_tide_input        => NULL() !< Control structure for a child module
  type(opacity_CS),             pointer :: opacity_CSp           => NULL() !< Control structure for a child module
  type(set_diffusivity_CS),     pointer :: set_diff_CSp          => NULL() !< Control structure for a child module
  type(sponge_CS),              pointer :: sponge_CSp            => NULL() !< Control structure for a child module
  type(ALE_sponge_CS),          pointer :: ALE_sponge_CSp        => NULL() !< Control structure for a child module
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp       => NULL() !< Control structure for a child module
  type(optics_type),            pointer :: optics                => NULL() !< Control structure for a child module
  type(KPP_CS),                 pointer :: KPP_CSp               => NULL() !< Control structure for a child module
  type(CVMix_conv_cs),          pointer :: CVMix_conv_csp        => NULL() !< Control structure for a child module
  type(diapyc_energy_req_CS),   pointer :: diapyc_en_rec_CSp     => NULL() !< Control structure for a child module

  type(group_pass_type) :: pass_hold_eb_ea !< For group halo pass
  type(group_pass_type) :: pass_Kv         !< For group halo pass
  type(diag_grid_storage) :: diag_grids_prev!< Stores diagnostic grids at some previous point in the algorithm
  ! Data arrays for communicating between components
  real, allocatable, dimension(:,:,:) :: KPP_NLTheat    !< KPP non-local transport for heat [m s-1]
  real, allocatable, dimension(:,:,:) :: KPP_NLTscalar  !< KPP non-local transport for scalars [m s-1]
  real, allocatable, dimension(:,:,:) :: KPP_buoy_flux  !< KPP forcing buoyancy flux [L2 T-3 ~> m2 s-3]
  real, allocatable, dimension(:,:)   :: KPP_temp_flux  !< KPP effective temperature flux [degC m s-1]
  real, allocatable, dimension(:,:)   :: KPP_salt_flux  !< KPP effective salt flux [ppt m s-1]

  type(time_type), pointer :: Time !< Pointer to model time (needed for sponges)
end type diabatic_CS

!>@{ clock ids
integer :: id_clock_entrain, id_clock_mixedlayer, id_clock_set_diffusivity
integer :: id_clock_tracers, id_clock_tridiag, id_clock_pass, id_clock_sponge
integer :: id_clock_geothermal, id_clock_differential_diff, id_clock_remap
integer :: id_clock_kpp
!>@}

contains

!>  This subroutine imposes the diapycnal mass fluxes and the
!!  accompanying diapycnal advection of momentum and tracers.
subroutine diabatic(u, v, h, tv, Hml, fluxes, visc, ADp, CDp, dt, Time_end, &
                    G, GV, US, CS, OBC, Waves)
  type(ocean_grid_type),                      intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV       !< ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h        !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv       !< points to thermodynamic fields
                                                                        !! unused have NULL ptrs
  real, dimension(:,:),                       pointer       :: Hml      !< Active mixed layer depth [Z ~> m]
  type(forcing),                              intent(inout) :: fluxes   !< points to forcing fields
                                                                        !! unused fields have NULL ptrs
  type(vertvisc_type),                        intent(inout) :: visc     !< vertical viscosities, BBL properies, and
  type(accel_diag_ptrs),                      intent(inout) :: ADp      !< related points to accelerations in momentum
                                                                        !! equations, to enable the later derived
                                                                        !! diagnostics, like energy budgets
  type(cont_diag_ptrs),                       intent(inout) :: CDp      !< points to terms in continuity equations
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
  type(unit_scale_type),                      intent(in)    :: US       !< A dimensional unit scaling type
  type(diabatic_CS),                          pointer       :: CS       !< module control structure
  type(ocean_OBC_type),             optional, pointer       :: OBC      !< Open boundaries control structure.
  type(Wave_parameters_CS),         optional, pointer       :: Waves    !< Surface gravity waves

  ! local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    eta      ! Interface heights before diapycnal mixing [m].
  real, dimension(SZI_(G),SZJ_(G),CS%nMode) :: &
    cn_IGW   ! baroclinic internal gravity wave speeds [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: temp_diag  ! Previous temperature for diagnostics [degC]
  integer :: i, j, k, m, is, ie, js, je, nz
  logical :: showCallTree ! If true, show the call tree

  if (GV%ke == 1) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_diabatic_driver: "// &
         "Module must be initialized before it is used.")
  if (dt == 0.0) call MOM_error(FATAL, "MOM_diabatic_driver: "// &
        "diabatic was called with a zero length timestep.")
  if (dt < 0.0) call MOM_error(FATAL, "MOM_diabatic_driver: "// &
        "diabatic was called with a negative timestep.")

  showCallTree = callTree_showQuery()

  ! Offer diagnostics of various state varables at the start of diabatic
  ! these are mostly for debugging purposes.
  if (CS%id_u_predia > 0) call post_data(CS%id_u_predia, u, CS%diag)
  if (CS%id_v_predia > 0) call post_data(CS%id_v_predia, v, CS%diag)
  if (CS%id_h_predia > 0) call post_data(CS%id_h_predia, h, CS%diag)
  if (CS%id_T_predia > 0) call post_data(CS%id_T_predia, tv%T, CS%diag)
  if (CS%id_S_predia > 0) call post_data(CS%id_S_predia, tv%S, CS%diag)
  if (CS%id_e_predia > 0) then
    call find_eta(h, tv, G, GV, US, eta, eta_to_m=1.0)
    call post_data(CS%id_e_predia, eta, CS%diag)
  endif

  if (CS%debug) then
    call MOM_state_chksum("Start of diabatic ", u, v, h, G, GV, US, haloshift=0)
    call MOM_forcing_chksum("Start of diabatic", fluxes, G, US, haloshift=0)
  endif
  if (CS%debugConservation) call MOM_state_stats('Start of diabatic', u, v, h, tv%T, tv%S, G, GV, US)

  if (CS%debug_energy_req) &
    call diapyc_energy_req_test(h, dt, tv, G, GV, US, CS%diapyc_en_rec_CSp)

  call cpu_clock_begin(id_clock_set_diffusivity)
  call set_BBL_TKE(u, v, h, fluxes, visc, G, GV, US, CS%set_diff_CSp, OBC=OBC)
  call cpu_clock_end(id_clock_set_diffusivity)

  ! Frazil formation keeps the temperature above the freezing point.
  ! make_frazil is deliberately called at both the beginning and at
  ! the end of the diabatic processes.
  if (associated(tv%T) .AND. associated(tv%frazil)) then
    ! For frazil diagnostic, the first call covers the first half of the time step
    call enable_averages(0.5*dt, Time_end - real_to_time(0.5*US%T_to_s*dt), CS%diag)
    if (CS%frazil_tendency_diag) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        temp_diag(i,j,k) = tv%T(i,j,k)
      enddo ; enddo ; enddo
    endif

    if (associated(fluxes%p_surf_full)) then
      call make_frazil(h, tv, G, GV, US, CS%diabatic_aux_CSp, fluxes%p_surf_full, halo=CS%halo_TS_diff)
    else
      call make_frazil(h, tv, G, GV, US, CS%diabatic_aux_CSp, halo=CS%halo_TS_diff)
    endif
    if (showCallTree) call callTree_waypoint("done with 1st make_frazil (diabatic)")

    if (CS%frazil_tendency_diag) then
      call diagnose_frazil_tendency(tv, h, temp_diag, 0.5*dt, G, GV, US, CS)
      if (CS%id_frazil_h > 0) call post_data(CS%id_frazil_h, h, CS%diag)
    endif
    call disable_averaging(CS%diag)
  endif ! associated(tv%T) .AND. associated(tv%frazil)
  if (CS%debugConservation) call MOM_state_stats('1st make_frazil', u, v, h, tv%T, tv%S, G, GV, US)

  if (CS%use_int_tides) then
    ! This block provides an interface for the unresolved low-mode internal tide module.
    call set_int_tide_input(u, v, h, tv, fluxes, CS%int_tide_input, dt, G, GV, US, &
                            CS%int_tide_input_CSp)
    cn_IGW(:,:,:) = 0.0
    if (CS%uniform_test_cg > 0.0) then
      do m=1,CS%nMode ; cn_IGW(:,:,m) = CS%uniform_test_cg ; enddo
    else
      call wave_speeds(h, tv, G, GV, US, CS%nMode, cn_IGW, full_halos=.true.)
    endif

    call propagate_int_tide(h, tv, cn_IGW, CS%int_tide_input%TKE_itidal_input, CS%int_tide_input%tideamp, &
                            CS%int_tide_input%Nb, dt, G, GV, US, CS%int_tide_CSp)
    if (showCallTree) call callTree_waypoint("done with propagate_int_tide (diabatic)")
  endif ! end CS%use_int_tides

  if (CS%useALEalgorithm .and. CS%use_legacy_diabatic) then
    call diabatic_ALE_legacy(u, v, h, tv, Hml, fluxes, visc, ADp, CDp, dt, Time_end, &
                      G, GV, US, CS, Waves)
  elseif (CS%useALEalgorithm) then
    call diabatic_ALE(u, v, h, tv, Hml, fluxes, visc, ADp, CDp, dt, Time_end, &
                      G, GV, US, CS, Waves)
  else
    call layered_diabatic(u, v, h, tv, Hml, fluxes, visc, ADp, CDp, dt, Time_end, &
                          G, GV, US, CS, Waves)
  endif


  call cpu_clock_begin(id_clock_pass)
  if (associated(visc%Kv_shear)) &
    call pass_var(visc%Kv_shear, G%Domain, To_All+Omit_Corners, halo=1)
  call cpu_clock_end(id_clock_pass)

  call disable_averaging(CS%diag)
  ! Frazil formation keeps temperature above the freezing point.
  ! make_frazil is deliberately called at both the beginning and at
  ! the end of the diabatic processes.
  if (associated(tv%T) .AND. associated(tv%frazil)) then
    call enable_averages(0.5*dt, Time_end, CS%diag)
    if (CS%frazil_tendency_diag) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        temp_diag(i,j,k) = tv%T(i,j,k)
      enddo ; enddo ; enddo
    endif

    if (associated(fluxes%p_surf_full)) then
      call make_frazil(h, tv, G, GV, US, CS%diabatic_aux_CSp, fluxes%p_surf_full)
    else
      call make_frazil(h, tv, G, GV, US, CS%diabatic_aux_CSp)
    endif

    if (CS%frazil_tendency_diag) then
      call diagnose_frazil_tendency(tv, h, temp_diag, 0.5*dt, G, GV, US, CS)
      if (CS%id_frazil_h > 0 ) call post_data(CS%id_frazil_h, h, CS%diag)
    endif

    if (showCallTree) call callTree_waypoint("done with 2nd make_frazil (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('2nd make_frazil', u, v, h, tv%T, tv%S, G, GV, US)
    call disable_averaging(CS%diag)

  endif  ! endif for frazil


  ! Diagnose mixed layer depths.
  call enable_averages(dt, Time_end, CS%diag)
  if (CS%id_MLD_003 > 0 .or. CS%id_subMLN2 > 0 .or. CS%id_mlotstsq > 0) then
    call diagnoseMLDbyDensityDifference(CS%id_MLD_003, h, tv, 0.03*US%kg_m3_to_R, G, GV, US, CS%diag, &
                                        id_N2subML=CS%id_subMLN2, id_MLDsq=CS%id_mlotstsq, dz_subML=CS%dz_subML_N2)
  endif
  if (CS%id_MLD_0125 > 0) then
    call diagnoseMLDbyDensityDifference(CS%id_MLD_0125, h, tv, 0.125*US%kg_m3_to_R, G, GV, US, CS%diag)
  endif
  if (CS%id_MLD_user > 0) then
    call diagnoseMLDbyDensityDifference(CS%id_MLD_user, h, tv, CS%MLDdensityDifference, G, GV, US, CS%diag)
  endif
  if ((CS%id_MLD_EN1 > 0) .or. (CS%id_MLD_EN2 > 0) .or. (CS%id_MLD_EN3 > 0)) then
    call diagnoseMLDbyEnergy((/CS%id_MLD_EN1, CS%id_MLD_EN2, CS%id_MLD_EN3/),&
                             h, tv, G, GV, US, CS%MLD_EN_VALS, CS%diag)
  endif
  if (CS%use_int_tides) then
    if (CS%id_cg1 > 0) call post_data(CS%id_cg1, cn_IGW(:,:,1),CS%diag)
    do m=1,CS%nMode ; if (CS%id_cn(m) > 0) call post_data(CS%id_cn(m), cn_IGW(:,:,m), CS%diag) ; enddo
  endif
  call disable_averaging(CS%diag)

  if (CS%debugConservation) call MOM_state_stats('leaving diabatic', u, v, h, tv%T, tv%S, G, GV, US)

end subroutine diabatic


!> Applies diabatic forcing and diapycnal mixing of temperature, salinity and other tracers for use
!! with an ALE algorithm.  This version uses an older set of algorithms compared with diabatic_ALE.
subroutine diabatic_ALE_legacy(u, v, h, tv, Hml, fluxes, visc, ADp, CDp, dt, Time_end, &
                           G, GV, US, CS, Waves)
  type(ocean_grid_type),                      intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV       !< ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US       !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h        !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv       !< points to thermodynamic fields
                                                                        !! unused have NULL ptrs
  real, dimension(:,:),                       pointer       :: Hml      !< Active mixed layer depth [Z ~> m]
  type(forcing),                              intent(inout) :: fluxes   !< points to forcing fields
                                                                        !! unused fields have NULL ptrs
  type(vertvisc_type),                        intent(inout) :: visc     !< vertical viscosities, BBL properies, and
  type(accel_diag_ptrs),                      intent(inout) :: ADp      !< related points to accelerations in momentum
                                                                        !! equations, to enable the later derived
                                                                        !! diagnostics, like energy budgets
  type(cont_diag_ptrs),                       intent(inout) :: CDp      !< points to terms in continuity equations
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
  type(diabatic_CS),                          pointer       :: CS       !< module control structure
  type(Wave_parameters_CS),         optional, pointer       :: Waves    !< Surface gravity waves

  ! local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    h_orig, &    ! Initial layer thicknesses [H ~> m or kg m-2]
    dSV_dT, &    ! The partial derivative of specific volume with temperature [R-1 degC-1 ~> m3 kg-1 degC-1]
    dSV_dS, &    ! The partial derivative of specific volume with salinity [R-1 ppt-1 ~> m3 kg-1 ppt-1].
    cTKE,   &    ! convective TKE requirements for each layer [R Z3 T-2 ~> J m-2].
    u_h,    &    ! Zonal velocities interpolated to thickness points [L T-1 ~> m s-1]
    v_h,    &    ! Meridional velocities interpolated to thickness points [L T-1 ~> m s-1]
    temp_diag, & ! Diagnostic array of previous temperatures [degC]
    saln_diag    ! Diagnostic array of previous salinity [ppt]

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    ent_s,    & ! The diffusive coupling across interfaces within one time step for
                ! salinity and passive tracers [H ~> m or kg m-2]
    ent_t,    & ! The diffusive coupling across interfaces within one time step for
                ! temperature [H ~> m or kg m-2]
    Kd_int,   & ! diapycnal diffusivity of interfaces [Z2 T-1 ~> m2 s-1]
    Kd_heat,  & ! diapycnal diffusivity of heat [Z2 T-1 ~> m2 s-1]
    Kd_salt,  & ! diapycnal diffusivity of salt and passive tracers [Z2 T-1 ~> m2 s-1]
    Kd_extra_T , & ! The extra diffusivity of temperature due to double diffusion relative to
                ! Kd_int [Z2 T-1 ~> m2 s-1].
    Kd_extra_S , & !  The extra diffusivity of salinity due to double diffusion relative to
                ! Kd_int [Z2 T-1 ~> m2 s-1].
    Kd_ePBL,  & ! test array of diapycnal diffusivities at interfaces [Z2 T-1 ~> m2 s-1]
    Tdif_flx, & ! diffusive diapycnal heat flux across interfaces [degC H T-1 ~> degC m s-1 or degC kg m-2 s-1]
    Sdif_flx    ! diffusive diapycnal salt flux across interfaces [ppt H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]

  real, dimension(SZI_(G),SZJ_(G)) :: &
    SkinBuoyFlux ! 2d surface buoyancy flux [Z2 T-3 ~> m2 s-3], used by ePBL

  logical, dimension(SZI_(G)) :: &
    in_boundary  ! True if there are no massive layers below, where massive is defined as
                 ! sufficiently thick that the no-flux boundary conditions have not restricted
                 ! the entrainment - usually sqrt(Kd*dt).

  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected [H ~> m or kg m-2]
  real :: h_neglect2   ! h_neglect^2 [H2 ~> m2 or kg2 m-4]
  real :: add_ent      ! Entrainment that needs to be added when mixing tracers [H ~> m or kg m-2]
  real :: I_hval       ! The inverse of the thicknesses averaged to interfaces [H-1 ~> m-1 or m2 kg-1]
  real :: h_tr         ! h_tr is h at tracer points with a tiny thickness
                       ! added to ensure positive definiteness [H ~> m or kg m-2]
  real :: Tr_ea_BBL    ! The diffusive tracer thickness in the BBL that is
                       ! coupled to the bottom within a timestep [H ~> m or kg m-2]
  real :: Kd_add_here    ! An added diffusivity [Z2 T-1 ~> m2 s-1].
  real :: htot(SZIB_(G)) ! The summed thickness from the bottom [H ~> m or kg m-2].

  real :: Ent_int ! The diffusive entrainment rate at an interface [H ~> m or kg m-2]
  real :: Idt     ! The inverse time step [T-1 ~> s-1]

  integer :: dir_flag     ! An integer encoding the directions in which to do halo updates.
  logical :: showCallTree ! If true, show the call tree
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, m

  integer :: ig, jg      ! global indices for testing testing itide point source (BDM)

  is   = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq  = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  h_neglect = GV%H_subroundoff ; h_neglect2 = h_neglect*h_neglect
  Kd_heat(:,:,:) = 0.0 ; Kd_salt(:,:,:) = 0.0

  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("diabatic_ALE_legacy(), MOM_diabatic_driver.F90")

  ! For all other diabatic subroutines, the averaging window should be the entire diabatic timestep
  call enable_averages(dt, Time_end, CS%diag)

  if (CS%use_geothermal) then
    call cpu_clock_begin(id_clock_geothermal)
    call geothermal_in_place(h, tv, dt, G, GV, US, CS%geothermal_CSp, halo=CS%halo_TS_diff)
    call cpu_clock_end(id_clock_geothermal)
    if (showCallTree) call callTree_waypoint("geothermal (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('geothermal', u, v, h, tv%T, tv%S, G, GV, US)
  endif

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  ! Set_pen_shortwave estimates the optical properties of the water column.
  ! It will need to be modified later to include information about the
  ! biological properties and layer thicknesses.
  if (associated(CS%optics)) &
    call set_pen_shortwave(CS%optics, fluxes, G, GV, US, CS%diabatic_aux_CSp, CS%opacity_CSp, CS%tracer_flow_CSp)

  if (CS%debug) call MOM_state_chksum("before find_uv_at_h", u, v, h, G, GV, US, haloshift=0)

  if (CS%use_kappa_shear .or. CS%use_CVMix_shear) then
    if (CS%use_geothermal) then
      call find_uv_at_h(u, v, h, u_h, v_h, G, GV, US, zero_mix=.true.)
    else
      call find_uv_at_h(u, v, h, u_h, v_h, G, GV, US)
    endif
    if (showCallTree) call callTree_waypoint("done with find_uv_at_h (diabatic)")
  endif

  call cpu_clock_begin(id_clock_set_diffusivity)
  ! Sets: Kd_int, Kd_extra_T, Kd_extra_S and visc%TKE_turb
  ! Also changes: visc%Kd_shear, visc%Kv_shear and visc%Kv_slow
  if (CS%debug) &
    call MOM_state_chksum("before set_diffusivity", u, v, h, G, GV, US, haloshift=CS%halo_TS_diff)
  if (CS%double_diffuse) then
    call set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, CS%optics, visc, dt, G, GV, US, CS%set_diff_CSp, &
                         Kd_int=Kd_int, Kd_extra_T=Kd_extra_T, Kd_extra_S=Kd_extra_S)
  else
    call set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, CS%optics, visc, dt, G, GV, US, &
                         CS%set_diff_CSp, Kd_int=Kd_int)
  endif
  call cpu_clock_end(id_clock_set_diffusivity)
  if (showCallTree) call callTree_waypoint("done with set_diffusivity (diabatic)")


  if (CS%debug) then
    call MOM_state_chksum("after set_diffusivity ", u, v, h, G, GV, US, haloshift=0)
    call MOM_forcing_chksum("after set_diffusivity ", fluxes, G, US, haloshift=0)
    call MOM_thermovar_chksum("after set_diffusivity ", tv, G)
    call hchksum(Kd_Int, "after set_diffusivity Kd_Int", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
  endif

  ! Set diffusivities for heat and salt separately
  if (CS%useKPP) then
    ! Add contribution from double diffusion
    if (CS%double_diffuse) then
      !$OMP parallel do default(shared)
      do K=1,nz+1 ; do j=js,je ; do i=is,ie
        Kd_salt(i,j,K) = Kd_int(i,j,K) + Kd_extra_S(i,j,K)
        Kd_heat(i,j,K) = Kd_int(i,j,K) + Kd_extra_T(i,j,K)
      enddo ; enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do K=1,nz+1 ; do j=js,je ; do i=is,ie
        Kd_salt(i,j,K) = Kd_int(i,j,K)
        Kd_heat(i,j,K) = Kd_int(i,j,K)
      enddo ; enddo ; enddo
    endif

    if (CS%debug) then
      call hchksum(Kd_heat, "after set_diffusivity Kd_heat", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
      call hchksum(Kd_salt, "after set_diffusivity Kd_salt", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
    endif

    call cpu_clock_begin(id_clock_kpp)
    ! total vertical viscosity in the interior is represented via visc%Kv_shear

    ! KPP needs the surface buoyancy flux but does not update state variables.
    ! We could make this call higher up to avoid a repeat unpacking of the surface fluxes.
    ! Sets: CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux
    ! NOTE: CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux are returned as rates (i.e. stuff per second)
    ! unlike other instances where the fluxes are integrated in time over a time-step.
    call calculateBuoyancyFlux2d(G, GV, US, fluxes, CS%optics, h, tv%T, tv%S, tv, &
                                 CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux)
    ! The KPP scheme calculates boundary layer diffusivities and non-local transport.

    call KPP_compute_BLD(CS%KPP_CSp, G, GV, US, h, tv%T, tv%S, u, v, tv, &
                         fluxes%ustar, CS%KPP_buoy_flux, Waves=Waves)

    call KPP_calculate(CS%KPP_CSp, G, GV, US, h, fluxes%ustar, CS%KPP_buoy_flux, Kd_heat, &
                       Kd_salt, visc%Kv_shear, CS%KPP_NLTheat, CS%KPP_NLTscalar, Waves=Waves)

    if (associated(Hml)) then
      call KPP_get_BLD(CS%KPP_CSp, Hml(:,:), G, US)
      call pass_var(Hml, G%domain, halo=1)
      ! If visc%MLD exists, copy KPP's BLD into it
      if (associated(visc%MLD)) visc%MLD(:,:) = Hml(:,:)
    endif

    if (.not.CS%KPPisPassive) then
      !$OMP parallel do default(shared)
      do K=1,nz+1 ; do j=js,je ; do i=is,ie
        Kd_int(i,j,K) = min( Kd_salt(i,j,K),  Kd_heat(i,j,K) )
      enddo ; enddo ; enddo
      if (CS%double_diffuse) then
        !$OMP parallel do default(shared)
        do K=1,nz+1 ; do j=js,je ; do i=is,ie
          Kd_extra_S(i,j,K) = (Kd_salt(i,j,K) - Kd_int(i,j,K))
          Kd_extra_T(i,j,K) = (Kd_heat(i,j,K) - Kd_int(i,j,K))
        enddo ; enddo ; enddo
      endif
    endif ! not passive

    if (showCallTree) call callTree_waypoint("done with KPP_calculate (diabatic)")
    if (CS%debug) then
      call MOM_state_chksum("after KPP", u, v, h, G, GV, US, haloshift=0)
      call MOM_forcing_chksum("after KPP", fluxes, G, US, haloshift=0)
      call MOM_thermovar_chksum("after KPP", tv, G)
      call hchksum(Kd_heat, "after KPP Kd_heat", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
      call hchksum(Kd_salt, "after KPP Kd_salt", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
      call hchksum(CS%KPP_temp_flux, "before KPP_applyNLT netHeat", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(CS%KPP_salt_flux, "before KPP_applyNLT netSalt", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(CS%KPP_NLTheat, "before KPP_applyNLT NLTheat", G%HI, haloshift=0)
      call hchksum(CS%KPP_NLTscalar, "before KPP_applyNLT NLTscalar", G%HI, haloshift=0)
    endif
    ! Apply non-local transport of heat and salt
    ! Changes: tv%T, tv%S
    call KPP_NonLocalTransport_temp(CS%KPP_CSp, G, GV, h, CS%KPP_NLTheat,   CS%KPP_temp_flux, &
                                    US%T_to_s*dt, tv%T, US%Q_to_J_kg*tv%C_p)
    call KPP_NonLocalTransport_saln(CS%KPP_CSp, G, GV, h, CS%KPP_NLTscalar, CS%KPP_salt_flux, &
                                    US%T_to_s*dt, tv%S)
    call cpu_clock_end(id_clock_kpp)
    if (showCallTree) call callTree_waypoint("done with KPP_applyNonLocalTransport (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('KPP_applyNonLocalTransport', u, v, h, tv%T, tv%S, G, GV, US)

    if (CS%debug) then
      call MOM_state_chksum("after KPP_applyNLT ", u, v, h, G, GV, US, haloshift=0)
      call MOM_forcing_chksum("after KPP_applyNLT ", fluxes, G, US, haloshift=0)
      call MOM_thermovar_chksum("after KPP_applyNLT ", tv, G)
    endif
  endif ! endif for KPP

  ! This is the "old" method for applying differential diffusion.
  ! Changes: tv%T, tv%S
  if (CS%double_diffuse .and. associated(tv%T)) then

    call cpu_clock_begin(id_clock_differential_diff)
    call differential_diffuse_T_S(h, tv%T, tv%S, Kd_extra_T, Kd_extra_S, dt, G, GV)
    call cpu_clock_end(id_clock_differential_diff)

    if (showCallTree) call callTree_waypoint("done with differential_diffuse_T_S (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('differential_diffuse_T_S', u, v, h, tv%T, tv%S, G, GV, US)

    ! increment heat and salt diffusivity.
    ! CS%useKPP==.true. already has extra_T and extra_S included
    if (.not. CS%useKPP) then
      !$OMP parallel do default(shared)
      do K=2,nz ; do j=js,je ; do i=is,ie
        Kd_heat(i,j,K) = Kd_heat(i,j,K) + Kd_extra_T(i,j,K)
        Kd_salt(i,j,K) = Kd_salt(i,j,K) + Kd_extra_S(i,j,K)
      enddo ; enddo ; enddo
    endif

  endif

  ! Calculate vertical mixing due to convection (computed via CVMix)
  if (CS%use_CVMix_conv) then
    ! Increment vertical diffusion and viscosity due to convection
    call calculate_CVMix_conv(h, tv, G, GV, US, CS%CVMix_conv_csp, Hml, Kd=Kd_int, Kv=visc%Kv_slow)
  endif

  ! This block sets ent_t and ent_s from h and Kd_int.
  do j=js,je ; do i=is,ie
    ent_s(i,j,1) = 0.0 ; ent_s(i,j,nz+1) = 0.0
    ent_t(i,j,1) = 0.0 ; ent_t(i,j,nz+1) = 0.0
  enddo ; enddo
  !$OMP parallel do default(shared)  private(I_hval)
  do K=2,nz ; do j=js,je ; do i=is,ie
    I_hval = 1.0 / (h_neglect + 0.5*(h(i,j,k-1) + h(i,j,k)))
    ent_s(i,j,K) = (GV%Z_to_H**2) * dt * I_hval * Kd_int(i,j,K)
    ent_t(i,j,K) = ent_s(i,j,K)
  enddo ; enddo ; enddo
  if (showCallTree) call callTree_waypoint("done setting ent_s and ent_t from Kd_int (diabatic)")

  if (CS%debug) then
    call MOM_forcing_chksum("after calc_entrain ", fluxes, G, US, haloshift=0)
    call MOM_thermovar_chksum("after calc_entrain ", tv, G)
    call MOM_state_chksum("after calc_entrain ", u, v, h, G, GV, US, haloshift=0)
    call hchksum(ent_s, "after calc_entrain ent_s", G%HI, haloshift=0, scale=GV%H_to_m)
  endif

  ! Save fields before boundary forcing is applied for tendency diagnostics
  do k=1,nz ; do j=js,je ; do i=is,ie
    h_orig(i,j,k) = h(i,j,k)
  enddo ; enddo ; enddo
  if (CS%boundary_forcing_tendency_diag) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      temp_diag(i,j,k) = tv%T(i,j,k)
      saln_diag(i,j,k) = tv%S(i,j,k)
    enddo ; enddo ; enddo
  endif

  ! Apply forcing
  ! Changes made to following fields:  h, tv%T and tv%S.
  call cpu_clock_begin(id_clock_remap)

  if (CS%use_energetic_PBL) then

    skinbuoyflux(:,:) = 0.0
    call applyBoundaryFluxesInOut(CS%diabatic_aux_CSp, G, GV, US, dt, fluxes, CS%optics, &
            optics_nbands(CS%optics), h, tv, CS%aggregate_FW_forcing, CS%evap_CFL_limit,                         &
            CS%minimum_forcing_depth, cTKE, dSV_dT, dSV_dS, SkinBuoyFlux=SkinBuoyFlux)

    if (CS%debug) then
      call hchksum(ent_t, "after applyBoundaryFluxes ent_t", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(ent_s, "after applyBoundaryFluxes ent_s", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(cTKE, "after applyBoundaryFluxes cTKE", G%HI, haloshift=0, &
                   scale=US%RZ3_T3_to_W_m2*US%T_to_s)
      call hchksum(dSV_dT, "after applyBoundaryFluxes dSV_dT", G%HI, haloshift=0, scale=US%kg_m3_to_R)
      call hchksum(dSV_dS, "after applyBoundaryFluxes dSV_dS", G%HI, haloshift=0, scale=US%kg_m3_to_R)
    endif

    call find_uv_at_h(u, v, h, u_h, v_h, G, GV, US)
    call energetic_PBL(h, u_h, v_h, tv, fluxes, dt, Kd_ePBL, G, GV, US, &
                       CS%energetic_PBL_CSp, dSV_dT, dSV_dS, cTKE, SkinBuoyFlux, waves=waves)

    if (associated(Hml)) then
      call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, Hml(:,:), G, US)
      call pass_var(Hml, G%domain, halo=1)
      ! If visc%MLD exists, copy ePBL's MLD into it
      if (associated(visc%MLD)) visc%MLD(:,:) = Hml(:,:)
    elseif (associated(visc%MLD)) then
      call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, visc%MLD, G, US)
      call pass_var(visc%MLD, G%domain, halo=1)
    endif

    ! Augment the diffusivities and viscosity due to those diagnosed in energetic_PBL.
    do K=2,nz ; do j=js,je ; do i=is,ie
      if (CS%ePBL_is_additive) then
        Kd_add_here = Kd_ePBL(i,j,K)
        visc%Kv_shear(i,j,K) = visc%Kv_shear(i,j,K) + CS%ePBL_Prandtl*Kd_ePBL(i,j,K)
      else
        Kd_add_here = max(Kd_ePBL(i,j,K) - visc%Kd_shear(i,j,K), 0.0)
        visc%Kv_shear(i,j,K) = max(visc%Kv_shear(i,j,K), CS%ePBL_Prandtl*Kd_ePBL(i,j,K))
      endif

      Ent_int = Kd_add_here * (GV%Z_to_H**2 * dt) / (0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect)
      ent_s(i,j,K) = ent_s(i,j,K) + Ent_int
      Kd_int(i,j,K) = Kd_int(i,j,K) + Kd_add_here

      ! for diagnostics
      Kd_heat(i,j,K) = Kd_heat(i,j,K) + Kd_int(i,j,K)
      Kd_salt(i,j,K) = Kd_salt(i,j,K) + Kd_int(i,j,K)

    enddo ; enddo ; enddo

    if (CS%debug) then
      call hchksum(ent_t, "after ePBL ent_t", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(ent_s, "after ePBL ent_s", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(Kd_ePBL, "after ePBL Kd_ePBL", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
    endif

  else
    call applyBoundaryFluxesInOut(CS%diabatic_aux_CSp, G, GV, US, dt, fluxes, CS%optics, &
                                  optics_nbands(CS%optics), h, tv, CS%aggregate_FW_forcing, &
                                  CS%evap_CFL_limit, CS%minimum_forcing_depth)

  endif   ! endif for CS%use_energetic_PBL

  ! diagnose the tendencies due to boundary forcing
  ! At this point, the diagnostic grids have not been updated since the call to the boundary layer scheme
  !  so all tendency diagnostics need to be posted on h_orig, and grids rebuilt afterwards
  if (CS%boundary_forcing_tendency_diag) then
    call diagnose_boundary_forcing_tendency(tv, h, temp_diag, saln_diag, h_orig, dt, G, GV, US, CS)
    if (CS%id_boundary_forcing_h > 0) call post_data(CS%id_boundary_forcing_h, h, CS%diag, alt_h=h_orig)
  endif
  ! Boundary fluxes may have changed T, S, and h
  call diag_update_remap_grids(CS%diag)
  call cpu_clock_end(id_clock_remap)
  if (CS%debug) then
    call MOM_forcing_chksum("after applyBoundaryFluxes ", fluxes, G, US, haloshift=0)
    call MOM_thermovar_chksum("after applyBoundaryFluxes ", tv, G)
    call MOM_state_chksum("after applyBoundaryFluxes ", u, v, h, G, GV, US, haloshift=0)
  endif
  if (showCallTree) call callTree_waypoint("done with applyBoundaryFluxes (diabatic)")
  if (CS%debugConservation)  call MOM_state_stats('applyBoundaryFluxes', u, v, h, tv%T, tv%S, G, GV, US)

  if (CS%debug) then
    call MOM_state_chksum("after negative check ", u, v, h, G, GV, US, haloshift=0)
    call MOM_forcing_chksum("after negative check ", fluxes, G, US, haloshift=0)
    call MOM_thermovar_chksum("after negative check ", tv, G)
  endif
  if (showCallTree) call callTree_waypoint("done with h=ea-eb (diabatic)")
  if (CS%debugConservation) call MOM_state_stats('h=ea-eb', u, v, h, tv%T, tv%S, G, GV, US)

  ! calculate change in temperature & salinity due to dia-coordinate surface diffusion
  if (associated(tv%T)) then

    if (CS%debug) then
      call hchksum(ent_t, "before triDiagTS ent_t ", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(ent_s, "before triDiagTS ent_s ", G%HI, haloshift=0, scale=GV%H_to_m)
    endif

    call cpu_clock_begin(id_clock_tridiag)
    !  Keep salinity from falling below a small but positive threshold.
    !  This constraint is needed for SIS1 ice model, which can extract
    !  more salt than is present in the ocean. SIS2 does not suffer
    !  from this limitation, in which case we can let salinity=0 and still
    !  have salt conserved with SIS2 ice. So for SIS2, we can run with
    !  BOUND_SALINITY=False in MOM.F90.
    if (associated(tv%S) .and. associated(tv%salt_deficit)) &
      call adjust_salt(h, tv, G, GV, CS%diabatic_aux_CSp)

    if (CS%diabatic_diff_tendency_diag) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        temp_diag(i,j,k) = tv%T(i,j,k)
        saln_diag(i,j,k) = tv%S(i,j,k)
      enddo ; enddo ; enddo
    endif

    ! Changes T and S via the tridiagonal solver; no change to h
    do k=1,nz ; do j=js,je ; do i=is,ie
      ent_t(i,j,K) = ent_s(i,j,K) ; ent_t(i,j,K+1) = ent_s(i,j,K+1)
    enddo ; enddo ; enddo
    if (CS%tracer_tridiag) then
      call tracer_vertdiff_Eulerian(h, ent_t, dt, tv%T, G, GV)
      call tracer_vertdiff_Eulerian(h, ent_s, dt, tv%S, G, GV)
    else
      call triDiagTS_Eulerian(G, GV, is, ie, js, je, h, ent_s, tv%T, tv%S)
    endif

    ! diagnose temperature, salinity, heat, and salt tendencies
    if (CS%diabatic_diff_tendency_diag) then
      call diagnose_diabatic_diff_tendency(tv, h, temp_diag, saln_diag, dt, G, GV, US, CS)
      if (CS%id_diabatic_diff_h > 0) call post_data(CS%id_diabatic_diff_h, h, CS%diag, alt_h=h)
    endif

    call cpu_clock_end(id_clock_tridiag)

    if (showCallTree) call callTree_waypoint("done with triDiagTS (diabatic)")

  endif  ! endif corresponding to if (associated(tv%T))

  if (CS%debugConservation) call MOM_state_stats('triDiagTS', u, v, h, tv%T, tv%S, G, GV, US)

  if (CS%debug) then
    call MOM_state_chksum("after mixed layer ", u, v, h, G, GV, US, haloshift=0)
    call MOM_thermovar_chksum("after mixed layer ", tv, G)
  endif

  ! Whenever thickness changes let the diag manager know, as the
  ! target grids for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  ! Diagnose the diapycnal diffusivities and other related quantities.
  if (CS%id_Kd_int  > 0) call post_data(CS%id_Kd_int,  Kd_int,  CS%diag)
  if (CS%id_Kd_heat > 0) call post_data(CS%id_Kd_heat, Kd_heat, CS%diag)
  if (CS%id_Kd_salt > 0) call post_data(CS%id_Kd_salt, Kd_salt, CS%diag)
  if (CS%id_Kd_ePBL > 0) call post_data(CS%id_Kd_ePBL, Kd_ePBL, CS%diag)

  if (CS%id_ea   > 0) call post_data(CS%id_ea,   ent_s(:,:,1:nz), CS%diag)
  if (CS%id_eb   > 0) call post_data(CS%id_eb,   ent_s(:,:,2:nz+1), CS%diag)
  if (CS%id_ea_t > 0) call post_data(CS%id_ea_t, ent_t(:,:,1:nz), CS%diag)
  if (CS%id_eb_t > 0) call post_data(CS%id_eb_t, ent_t(:,:,2:nz+1), CS%diag)
  if (CS%id_ea_s > 0) call post_data(CS%id_ea_s, ent_s(:,:,1:nz), CS%diag)
  if (CS%id_eb_s > 0) call post_data(CS%id_eb_s, ent_s(:,:,2:nz+1), CS%diag)
  Idt = 1.0 / dt
  if (CS%id_Tdif > 0) then
    do j=js,je ; do i=is,ie
      Tdif_flx(i,j,1) = 0.0 ; Tdif_flx(i,j,nz+1) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do K=2,nz ; do j=js,je ; do i=is,ie
      Tdif_flx(i,j,K) = (Idt * ent_t(i,j,K)) * (tv%T(i,j,k-1) - tv%T(i,j,k))
    enddo ; enddo ; enddo
    if (CS%id_Tdif > 0) call post_data(CS%id_Tdif, Tdif_flx, CS%diag)
  endif
  if (CS%id_Sdif > 0) then
    do j=js,je ; do i=is,ie
      Sdif_flx(i,j,1) = 0.0 ; Sdif_flx(i,j,nz+1) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do K=2,nz ; do j=js,je ; do i=is,ie
      Sdif_flx(i,j,K) = (Idt * ent_s(i,j,K)) * (tv%S(i,j,k-1) - tv%S(i,j,k))
    enddo ; enddo ; enddo
    if (CS%id_Sdif > 0) call post_data(CS%id_Sdif, Sdif_flx, CS%diag)
  endif

  ! mixing of passive tracers from massless boundary layers to interior
  call cpu_clock_begin(id_clock_tracers)

  if (CS%mix_boundary_tracer_ALE) then
    Tr_ea_BBL = GV%Z_to_H * sqrt(dt*CS%Kd_BBL_tr)

    !$OMP parallel do default(shared) private(htot,in_boundary,add_ent)
    do j=js,je
      do i=is,ie
        htot(i) = 0.0
        in_boundary(i) = (G%mask2dT(i,j) > 0.0)
      enddo
      do k=nz,2,-1 ; do i=is,ie
        if (in_boundary(i)) then
          htot(i) = htot(i) + h(i,j,k)
          !   If diapycnal mixing has been suppressed because this is a massless
          ! layer near the bottom, add some mixing of tracers between these
          ! layers.  This flux is based on the harmonic mean of the two
          ! thicknesses, as this corresponds pretty closely (to within
          ! differences in the density jumps between layers) with what is done
          ! in the calculation of the fluxes in the first place.  Kd_min_tr
          ! should be much less than the values that have been set in Kd_int,
          ! perhaps a molecular diffusivity.
          add_ent = ((dt * CS%Kd_min_tr) * GV%Z_to_H**2) * &
                    ((h(i,j,k-1)+h(i,j,k)+h_neglect) /  (h(i,j,k-1)*h(i,j,k)+h_neglect2)) - &
                    0.5*(ent_s(i,j,K) + ent_s(i,j,K))
          if (htot(i) < Tr_ea_BBL) then
            add_ent = max(0.0, add_ent, (Tr_ea_BBL - htot(i)) - ent_s(i,j,K))
          elseif (add_ent < 0.0) then
            add_ent = 0.0 ; in_boundary(i) = .false.
          endif

          ent_s(i,j,K) = ent_s(i,j,K) + add_ent
        endif

        if (CS%double_diffuse) then ; if (Kd_extra_S(i,j,k) > 0.0) then
          add_ent = ((dt * Kd_extra_S(i,j,k)) * GV%Z_to_H**2) / &
                    (0.5 * (h(i,j,k-1) + h(i,j,k)) +  h_neglect)
          ent_s(i,j,K) = ent_s(i,j,K) + add_ent
        endif ; endif
      enddo ; enddo

    enddo
  elseif (CS%double_diffuse .and. .not.CS%mix_boundary_tracers) then  ! extra diffusivity for passive tracers
    !$OMP parallel do default(shared) private(add_ent)
    do k=nz,2,-1 ; do j=js,je ; do i=is,ie
      if (Kd_extra_S(i,j,k) > 0.0) then
        add_ent = ((dt * Kd_extra_S(i,j,k)) * GV%Z_to_H**2) / &
                  (0.5 * (h(i,j,k-1) + h(i,j,k)) + h_neglect)
      else
        add_ent = 0.0
      endif
      ent_s(i,j,K) = ent_s(i,j,K) + add_ent
    enddo ; enddo ; enddo
  endif  ! (CS%mix_boundary_tracers)

  ! For passive tracers, the changes in thickness due to boundary fluxes has yet to be applied
  call call_tracer_column_fns(h_orig, h, ent_s(:,:,1:nz), ent_s(:,:,2:nz+1), fluxes, Hml, dt, &
                              G, GV, US, tv, CS%optics, CS%tracer_flow_CSp, CS%debug, &
                              evap_CFL_limit = CS%evap_CFL_limit, &
                              minimum_forcing_depth=CS%minimum_forcing_depth)

  call cpu_clock_end(id_clock_tracers)

  ! Apply ALE sponge
  if (CS%use_sponge .and. associated(CS%ALE_sponge_CSp)) then
    call cpu_clock_begin(id_clock_sponge)
    call apply_ALE_sponge(h, dt, G, GV, US, CS%ALE_sponge_CSp, CS%Time)
    call cpu_clock_end(id_clock_sponge)
    if (CS%debug) then
      call MOM_state_chksum("apply_sponge ", u, v, h, G, GV, US, haloshift=0)
      call MOM_thermovar_chksum("apply_sponge ", tv, G)
    endif
  endif ! CS%use_sponge

  call disable_averaging(CS%diag)

  if (showCallTree) call callTree_leave("diabatic_ALE_legacy()")

end subroutine diabatic_ALE_legacy


!>  This subroutine imposes the diapycnal mass fluxes and the
!!  accompanying diapycnal advection of momentum and tracers.
subroutine diabatic_ALE(u, v, h, tv, Hml, fluxes, visc, ADp, CDp, dt, Time_end, &
                        G, GV, US, CS, Waves)
  type(ocean_grid_type),                      intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV       !< ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US       !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h        !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv       !< points to thermodynamic fields
                                                                        !! unused have NULL ptrs
  real, dimension(:,:),                       pointer       :: Hml      !< Active mixed layer depth [Z ~> m]
  type(forcing),                              intent(inout) :: fluxes   !< points to forcing fields
                                                                        !! unused fields have NULL ptrs
  type(vertvisc_type),                        intent(inout) :: visc     !< vertical viscosities, BBL properies, and
  type(accel_diag_ptrs),                      intent(inout) :: ADp      !< related points to accelerations in momentum
                                                                        !! equations, to enable the later derived
                                                                        !! diagnostics, like energy budgets
  type(cont_diag_ptrs),                       intent(inout) :: CDp      !< points to terms in continuity equations
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
  type(diabatic_CS),                          pointer       :: CS       !< module control structure
  type(Wave_parameters_CS),         optional, pointer       :: Waves    !< Surface gravity waves

  ! local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    h_orig, &    ! Initial layer thicknesses [H ~> m or kg m-2]
    dSV_dT, &    ! The partial derivative of specific volume with temperature [R-1 degC-1 ~> m3 kg-1 degC-1]
    dSV_dS, &    ! The partial derivative of specific volume with salinity [R-1 ppt-1 ~> m3 kg-1 ppt-1].
    cTKE,   &    ! convective TKE requirements for each layer [R Z3 T-2 ~> J m-2].
    u_h,    &    ! Zonal velocities interpolated to thickness points [L T-1 ~> m s-1]
    v_h,    &    ! Meridional velocities interpolated to thickness points [L T-1 ~> m s-1]
    temp_diag, & ! Diagnostic array of previous temperatures [degC]
    saln_diag    ! Diagnostic array of previous salinity [ppt]

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    ent_s,    & ! The diffusive coupling across interfaces within one time step for
                ! salinity and passive tracers [H ~> m or kg m-2]
    ent_t,    & ! The diffusive coupling across interfaces within one time step for
                ! temperature [H ~> m or kg m-2]
    Kd_heat,  & ! diapycnal diffusivity of heat or the smaller of the diapycnal diffusivities of
                ! heat and salt [Z2 T-1 ~> m2 s-1]
    Kd_salt,  & ! diapycnal diffusivity of salt and passive tracers [Z2 T-1 ~> m2 s-1]
    Kd_extra_T , & ! The extra diffusivity of temperature due to double diffusion relative to
                ! Kd_int returned from set_diffusivity [Z2 T-1 ~> m2 s-1].
    Kd_extra_S , & !  The extra diffusivity of salinity due to double diffusion relative to
                ! Kd_int returned from set_diffusivity [Z2 T-1 ~> m2 s-1].
    Kd_ePBL,  & ! boundary layer or convective diapycnal diffusivities at interfaces [Z2 T-1 ~> m2 s-1]
    Tdif_flx, & ! diffusive diapycnal heat flux across interfaces [degC H T-1 ~> degC m s-1 or degC kg m-2 s-1]
    Sdif_flx    ! diffusive diapycnal salt flux across interfaces [ppt H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]

  real, dimension(SZI_(G),SZJ_(G)) :: &
    SkinBuoyFlux ! 2d surface buoyancy flux [Z2 T-3 ~> m2 s-3], used by ePBL

  logical, dimension(SZI_(G)) :: &
    in_boundary  ! True if there are no massive layers below, where massive is defined as
                 ! sufficiently thick that the no-flux boundary conditions have not restricted
                 ! the entrainment - usually sqrt(Kd*dt).

  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected [H ~> m or kg m-2]
  real :: h_neglect2   ! h_neglect^2 [H2 ~> m2 or kg2 m-4]
  real :: add_ent      ! Entrainment that needs to be added when mixing tracers [H ~> m or kg m-2]
  real :: I_hval       ! The inverse of the thicknesses averaged to interfaces [H-1 ~> m-1 or m2 kg-1]
  real :: h_tr         ! h_tr is h at tracer points with a tiny thickness
                       ! added to ensure positive definiteness [H ~> m or kg m-2]
  real :: Tr_ea_BBL    ! The diffusive tracer thickness in the BBL that is
                       ! coupled to the bottom within a timestep [H ~> m or kg m-2]
  real :: htot(SZIB_(G)) ! The summed thickness from the bottom [H ~> m or kg m-2].
  real :: Kd_add_here    ! An added diffusivity [Z2 T-1 ~> m2 s-1].
  real :: Idt     ! The inverse time step [T-1 ~> s-1]

  integer :: dir_flag     ! An integer encoding the directions in which to do halo updates.
  logical :: showCallTree ! If true, show the call tree
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, m

  is   = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq  = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  h_neglect = GV%H_subroundoff ; h_neglect2 = h_neglect*h_neglect
  Kd_heat(:,:,:) = 0.0 ; Kd_salt(:,:,:) = 0.0
  ent_s(:,:,:) = 0.0 ; ent_t(:,:,:) = 0.0

  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("diabatic_ALE(), MOM_diabatic_driver.F90")

  if (.not. (CS%useALEalgorithm)) call MOM_error(FATAL, "MOM_diabatic_driver: "// &
         "The ALE algorithm must be enabled when using MOM_diabatic_driver.")

  ! For all other diabatic subroutines, the averaging window should be the entire diabatic timestep
  call enable_averages(dt, Time_end, CS%diag)

  if (CS%use_geothermal) then
    call cpu_clock_begin(id_clock_geothermal)
    call geothermal_in_place(h, tv, dt, G, GV, US, CS%geothermal_CSp, halo=CS%halo_TS_diff)
    call cpu_clock_end(id_clock_geothermal)
    if (showCallTree) call callTree_waypoint("geothermal (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('geothermal', u, v, h, tv%T, tv%S, G, GV, US)
  endif

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  ! Set_pen_shortwave estimates the optical properties of the water column.
  ! It will need to be modified later to include information about the
  ! biological properties and layer thicknesses.
  if (associated(CS%optics)) &
    call set_pen_shortwave(CS%optics, fluxes, G, GV, US, CS%diabatic_aux_CSp, CS%opacity_CSp, CS%tracer_flow_CSp)

  if (CS%debug) call MOM_state_chksum("before find_uv_at_h", u, v, h, G, GV, US, haloshift=0)

  if (CS%use_kappa_shear .or. CS%use_CVMix_shear) then
    if (CS%use_geothermal) then
      call find_uv_at_h(u, v, h, u_h, v_h, G, GV, US, zero_mix=.true.)
    else
      call find_uv_at_h(u, v, h, u_h, v_h, G, GV, US)
    endif
    if (showCallTree) call callTree_waypoint("done with find_uv_at_h (diabatic)")
  endif

  call cpu_clock_begin(id_clock_set_diffusivity)
  ! Sets: Kd_heat, Kd_extra_T, Kd_extra_S and visc%TKE_turb
  ! Also changes: visc%Kd_shear, visc%Kv_shear and visc%Kv_slow
  if (CS%debug) &
    call MOM_state_chksum("before set_diffusivity", u, v, h, G, GV, US, haloshift=CS%halo_TS_diff)
  if (CS%double_diffuse) then
    call set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, CS%optics, visc, dt, G, GV, US, CS%set_diff_CSp, &
                         Kd_int=Kd_heat, Kd_extra_T=Kd_extra_T, Kd_extra_S=Kd_extra_S)
  else
    call set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, CS%optics, visc, dt, G, GV, US, &
                         CS%set_diff_CSp, Kd_int=Kd_heat)
  endif
  call cpu_clock_end(id_clock_set_diffusivity)
  if (showCallTree) call callTree_waypoint("done with set_diffusivity (diabatic)")

  if (CS%debug) then
    call MOM_state_chksum("after set_diffusivity ", u, v, h, G, GV, US, haloshift=0)
    call MOM_forcing_chksum("after set_diffusivity ", fluxes, G, US, haloshift=0)
    call MOM_thermovar_chksum("after set_diffusivity ", tv, G)
    call hchksum(Kd_heat, "after set_diffusivity Kd_heat", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
  endif

  ! Store the diagnosed typical diffusivity at interfaces.
  if (CS%id_Kd_int > 0) call post_data(CS%id_Kd_int, Kd_heat,  CS%diag)

  ! Set diffusivities for heat and salt separately, and possibly change the meaning of Kd_heat.
  if (CS%double_diffuse) then
    ! Add contributions from double diffusion
    !$OMP parallel do default(shared)
    do K=1,nz+1 ; do j=js,je ; do i=is,ie
      Kd_salt(i,j,K) = Kd_heat(i,j,K) + Kd_extra_S(i,j,K)
      Kd_heat(i,j,K) = Kd_heat(i,j,K) + Kd_extra_T(i,j,K)
    enddo ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do K=1,nz+1 ; do j=js,je ; do i=is,ie
      Kd_salt(i,j,K) = Kd_heat(i,j,K)
    enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call hchksum(Kd_heat, "after double diffuse Kd_heat", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
    call hchksum(Kd_salt, "after double diffuse Kd_salt", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
  endif

  if (CS%useKPP) then
    call cpu_clock_begin(id_clock_kpp)
    ! total vertical viscosity in the interior is represented via visc%Kv_shear
    do k=1,nz+1 ; do j=js,je ; do i=is,ie
      visc%Kv_shear(i,j,k) = visc%Kv_shear(i,j,k) + visc%Kv_slow(i,j,k)
    enddo ; enddo ; enddo

    ! KPP needs the surface buoyancy flux but does not update state variables.
    ! We could make this call higher up to avoid a repeat unpacking of the surface fluxes.
    ! Sets: CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux
    ! NOTE: CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux are returned as rates (i.e. stuff per second)
    ! unlike other instances where the fluxes are integrated in time over a time-step.
    call calculateBuoyancyFlux2d(G, GV, US, fluxes, CS%optics, h, tv%T, tv%S, tv, &
                                 CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux)

    ! The KPP scheme calculates boundary layer diffusivities and non-local transport.
    call KPP_compute_BLD(CS%KPP_CSp, G, GV, US, h, tv%T, tv%S, u, v, tv, &
                         fluxes%ustar, CS%KPP_buoy_flux, Waves=Waves)

    call KPP_calculate(CS%KPP_CSp, G, GV, US, h, fluxes%ustar, CS%KPP_buoy_flux, Kd_heat, &
                       Kd_salt, visc%Kv_shear, CS%KPP_NLTheat, CS%KPP_NLTscalar, Waves=Waves)

    if (associated(Hml)) then
      call KPP_get_BLD(CS%KPP_CSp, Hml(:,:), G, US)
      call pass_var(Hml, G%domain, halo=1)
      ! If visc%MLD exists, copy KPP's BLD into it
      if (associated(visc%MLD)) visc%MLD(:,:) = Hml(:,:)
    endif

    if (showCallTree) call callTree_waypoint("done with KPP_calculate (diabatic)")
    if (CS%debug) then
      call MOM_state_chksum("after KPP", u, v, h, G, GV, US, haloshift=0)
      call MOM_forcing_chksum("after KPP", fluxes, G, US, haloshift=0)
      call MOM_thermovar_chksum("after KPP", tv, G)
      call hchksum(Kd_heat, "after KPP Kd_heat", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
      call hchksum(Kd_salt, "after KPP Kd_salt", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
      call hchksum(CS%KPP_temp_flux, "before KPP_applyNLT netHeat", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(CS%KPP_salt_flux, "before KPP_applyNLT netSalt", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(CS%KPP_NLTheat, "before KPP_applyNLT NLTheat", G%HI, haloshift=0)
      call hchksum(CS%KPP_NLTscalar, "before KPP_applyNLT NLTscalar", G%HI, haloshift=0)
    endif
    ! Apply non-local transport of heat and salt
    ! Changes: tv%T, tv%S
    call KPP_NonLocalTransport_temp(CS%KPP_CSp, G, GV, h, CS%KPP_NLTheat,   CS%KPP_temp_flux, &
                                    US%T_to_s*dt, tv%T, US%Q_to_J_kg*tv%C_p)
    call KPP_NonLocalTransport_saln(CS%KPP_CSp, G, GV, h, CS%KPP_NLTscalar, CS%KPP_salt_flux, &
                                    US%T_to_s*dt, tv%S)
    call cpu_clock_end(id_clock_kpp)
    if (showCallTree) call callTree_waypoint("done with KPP_applyNonLocalTransport (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('KPP_applyNonLocalTransport', u, v, h, tv%T, tv%S, G, GV, US)

    if (CS%debug) then
      call MOM_state_chksum("after KPP_applyNLT ", u, v, h, G, GV, US, haloshift=0)
      call MOM_forcing_chksum("after KPP_applyNLT ", fluxes, G, US, haloshift=0)
      call MOM_thermovar_chksum("after KPP_applyNLT ", tv, G)
    endif
  endif ! endif for KPP

  ! Calculate vertical mixing due to convection (computed via CVMix)
  if (CS%use_CVMix_conv) then
    ! Increment vertical diffusion and viscosity due to convection
    if (CS%useKPP) then
      call calculate_CVMix_conv(h, tv, G, GV, US, CS%CVMix_conv_csp, Hml, Kd=Kd_heat, Kv=visc%Kv_shear, Kd_aux=Kd_salt)
    else
      call calculate_CVMix_conv(h, tv, G, GV, US, CS%CVMix_conv_csp, Hml, Kd=Kd_heat, Kv=visc%Kv_slow, Kd_aux=Kd_salt)
    endif
  endif

  ! Save fields before boundary forcing is applied for tendency diagnostics
  do k=1,nz ; do j=js,je ; do i=is,ie
    h_orig(i,j,k) = h(i,j,k)
  enddo ; enddo ; enddo
  if (CS%boundary_forcing_tendency_diag) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      temp_diag(i,j,k) = tv%T(i,j,k)
      saln_diag(i,j,k) = tv%S(i,j,k)
    enddo ; enddo ; enddo
  endif

  ! Apply forcing
  ! Changes made to following fields:  h, tv%T and tv%S.
  call cpu_clock_begin(id_clock_remap)

  if (CS%use_energetic_PBL) then

    skinbuoyflux(:,:) = 0.0
    call applyBoundaryFluxesInOut(CS%diabatic_aux_CSp, G, GV, US, dt, fluxes, CS%optics, &
            optics_nbands(CS%optics), h, tv, CS%aggregate_FW_forcing, CS%evap_CFL_limit, &
            CS%minimum_forcing_depth, cTKE, dSV_dT, dSV_dS, SkinBuoyFlux=SkinBuoyFlux)

    if (CS%debug) then
      call hchksum(ent_t, "after applyBoundaryFluxes ent_t", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(ent_s, "after applyBoundaryFluxes ent_s", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(cTKE, "after applyBoundaryFluxes cTKE", G%HI, haloshift=0, &
                   scale=US%RZ3_T3_to_W_m2*US%T_to_s)
      call hchksum(dSV_dT, "after applyBoundaryFluxes dSV_dT", G%HI, haloshift=0, scale=US%kg_m3_to_R)
      call hchksum(dSV_dS, "after applyBoundaryFluxes dSV_dS", G%HI, haloshift=0, scale=US%kg_m3_to_R)
    endif

    call find_uv_at_h(u, v, h, u_h, v_h, G, GV, US)
    call energetic_PBL(h, u_h, v_h, tv, fluxes, dt, Kd_ePBL, G, GV, US, &
                       CS%energetic_PBL_CSp, dSV_dT, dSV_dS, cTKE, SkinBuoyFlux, waves=waves)

    if (associated(Hml)) then
      call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, Hml(:,:), G, US)
      call pass_var(Hml, G%domain, halo=1)
      ! If visc%MLD exists, copy ePBL's MLD into it
      if (associated(visc%MLD)) visc%MLD(:,:) = Hml(:,:)
    elseif (associated(visc%MLD)) then
      call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, visc%MLD, G, US)
      call pass_var(visc%MLD, G%domain, halo=1)
    endif

    ! Augment the diffusivities and viscosity due to those diagnosed in energetic_PBL.
    do K=2,nz ; do j=js,je ; do i=is,ie
      if (CS%ePBL_is_additive) then
        Kd_add_here = Kd_ePBL(i,j,K)
        visc%Kv_shear(i,j,K) = visc%Kv_shear(i,j,K) + CS%ePBL_Prandtl*Kd_ePBL(i,j,K)
      else
        Kd_add_here = max(Kd_ePBL(i,j,K) - visc%Kd_shear(i,j,K), 0.0)
        visc%Kv_shear(i,j,K) = max(visc%Kv_shear(i,j,K), CS%ePBL_Prandtl*Kd_ePBL(i,j,K))
      endif

      Kd_heat(i,j,K) = Kd_heat(i,j,K) + Kd_add_here
      Kd_salt(i,j,K) = Kd_salt(i,j,K) + Kd_add_here
    enddo ; enddo ; enddo

    if (CS%debug) then
      call hchksum(ent_t, "after ePBL ent_t", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(ent_s, "after ePBL ent_s", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(Kd_ePBL, "after ePBL Kd_ePBL", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
    endif

  else
    call applyBoundaryFluxesInOut(CS%diabatic_aux_CSp, G, GV, US, dt, fluxes, CS%optics, &
                                  optics_nbands(CS%optics), h, tv, CS%aggregate_FW_forcing, &
                                  CS%evap_CFL_limit, CS%minimum_forcing_depth)

  endif   ! endif for CS%use_energetic_PBL

  ! diagnose the tendencies due to boundary forcing
  ! At this point, the diagnostic grids have not been updated since the call to the boundary layer scheme
  !  so all tendency diagnostics need to be posted on h_orig, and grids rebuilt afterwards
  if (CS%boundary_forcing_tendency_diag) then
    call diagnose_boundary_forcing_tendency(tv, h, temp_diag, saln_diag, h_orig, dt, G, GV, US, CS)
    if (CS%id_boundary_forcing_h > 0) call post_data(CS%id_boundary_forcing_h, h, CS%diag, alt_h=h_orig)
  endif
  ! Boundary fluxes may have changed T, S, and h
  call diag_update_remap_grids(CS%diag)
  call cpu_clock_end(id_clock_remap)
  if (CS%debug) then
    call MOM_forcing_chksum("after applyBoundaryFluxes ", fluxes, G, US, haloshift=0)
    call MOM_thermovar_chksum("after applyBoundaryFluxes ", tv, G)
    call MOM_state_chksum("after applyBoundaryFluxes ", u, v, h, G, GV, US, haloshift=0)
  endif
  if (showCallTree) call callTree_waypoint("done with applyBoundaryFluxes (diabatic)")
  if (CS%debugConservation)  call MOM_state_stats('applyBoundaryFluxes', u, v, h, tv%T, tv%S, G, GV, US)

  if (showCallTree) call callTree_waypoint("done with h=ea-eb (diabatic)")
  if (CS%debugConservation) call MOM_state_stats('h=ea-eb', u, v, h, tv%T, tv%S, G, GV, US)

  ! calculate change in temperature & salinity due to dia-coordinate surface diffusion
  if (associated(tv%T)) then

    if (CS%debug) then
      call hchksum(ent_t, "before triDiagTS ent_t ", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(ent_s, "before triDiagTS ent_s ", G%HI, haloshift=0, scale=GV%H_to_m)
    endif

    call cpu_clock_begin(id_clock_tridiag)
    !  Keep salinity from falling below a small but positive threshold.
    !  This constraint is needed for SIS1 ice model, which can extract
    !  more salt than is present in the ocean. SIS2 does not suffer
    !  from this limitation, in which case we can let salinity=0 and still
    !  have salt conserved with SIS2 ice. So for SIS2, we can run with
    !  BOUND_SALINITY=False in MOM.F90.
    if (associated(tv%S) .and. associated(tv%salt_deficit)) &
      call adjust_salt(h, tv, G, GV, CS%diabatic_aux_CSp)

    if (CS%diabatic_diff_tendency_diag) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        temp_diag(i,j,k) = tv%T(i,j,k)
        saln_diag(i,j,k) = tv%S(i,j,k)
      enddo ; enddo ; enddo
    endif

    ! set ent_t=dt*Kd_heat/h_int and est_s=dt*Kd_salt/h_int on interfaces for use in the tridiagonal solver.
    do j=js,je ; do i=is,ie
      ent_t(i,j,1) = 0. ; ent_t(i,j,nz+1) = 0.
      ent_s(i,j,1) = 0. ; ent_s(i,j,nz+1) = 0.
    enddo ; enddo

    !$OMP parallel do default(shared) private(I_hval)
    do K=2,nz ; do j=js,je ; do i=is,ie
      I_hval = 1.0 / (h_neglect + 0.5*(h(i,j,k-1) + h(i,j,k)))
      ent_t(i,j,K) = (GV%Z_to_H**2) * dt * I_hval * Kd_heat(i,j,k)
      ent_s(i,j,K) = (GV%Z_to_H**2) * dt * I_hval * Kd_salt(i,j,k)
    enddo ; enddo ; enddo
    if (showCallTree) call callTree_waypoint("done setting ent_t and ent_t from Kd_heat and " //&
                                             "Kd_salt (diabatic_ALE)")

    ! Changes T and S via the tridiagonal solver; no change to h
    call tracer_vertdiff_Eulerian(h, ent_t, dt, tv%T, G, GV)
    call tracer_vertdiff_Eulerian(h, ent_s, dt, tv%S, G, GV)

    ! In ALE-mode, layer thicknesses do not change. Therefore, we can use h below
    if (CS%diabatic_diff_tendency_diag) then
      call diagnose_diabatic_diff_tendency(tv, h, temp_diag, saln_diag, dt, G, GV, US, CS)
    endif
    call cpu_clock_end(id_clock_tridiag)

    if (showCallTree) call callTree_waypoint("done with triDiagTS (diabatic)")

  endif  ! endif corresponding to if (associated(tv%T))

  if (CS%debugConservation) call MOM_state_stats('triDiagTS', u, v, h, tv%T, tv%S, G, GV, US)

  if (CS%debug) then
    call MOM_state_chksum("after mixed layer ", u, v, h, G, GV, US, haloshift=0)
    call MOM_thermovar_chksum("after mixed layer ", tv, G)
  endif

  ! Whenever thickness changes let the diag manager know, as the
  ! target grids for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  ! Diagnose the diapycnal diffusivities and other related quantities.
  if (CS%id_Kd_heat      > 0) call post_data(CS%id_Kd_heat,      Kd_heat, CS%diag)
  if (CS%id_Kd_salt      > 0) call post_data(CS%id_Kd_salt,      Kd_salt, CS%diag)
  if (CS%id_Kd_ePBL      > 0) call post_data(CS%id_Kd_ePBL,      Kd_ePBL, CS%diag)

  if (CS%id_ea_t       > 0) call post_data(CS%id_ea_t, ent_t(:,:,1:nz), CS%diag)
  if (CS%id_eb_t       > 0) call post_data(CS%id_eb_t, ent_t(:,:,2:nz+1), CS%diag)
  if (CS%id_ea_s       > 0) call post_data(CS%id_ea_s, ent_s(:,:,1:nz), CS%diag)
  if (CS%id_eb_s       > 0) call post_data(CS%id_eb_s, ent_s(:,:,2:nz+1), CS%diag)

  Idt = 1.0 / dt
  if (CS%id_Tdif > 0) then
    do j=js,je ; do i=is,ie
      Tdif_flx(i,j,1) = 0.0 ; Tdif_flx(i,j,nz+1) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do K=2,nz ; do j=js,je ; do i=is,ie
      Tdif_flx(i,j,K) = (Idt * ent_t(i,j,K)) * (tv%T(i,j,k-1) - tv%T(i,j,k))
    enddo ; enddo ; enddo
    if (CS%id_Tdif > 0) call post_data(CS%id_Tdif, Tdif_flx, CS%diag)
  endif
  if (CS%id_Sdif > 0) then
    do j=js,je ; do i=is,ie
      Sdif_flx(i,j,1) = 0.0 ; Sdif_flx(i,j,nz+1) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do K=2,nz ; do j=js,je ; do i=is,ie
      Sdif_flx(i,j,K) = (Idt * ent_s(i,j,K)) * (tv%S(i,j,k-1) - tv%S(i,j,k))
    enddo ; enddo ; enddo
    if (CS%id_Sdif > 0) call post_data(CS%id_Sdif, Sdif_flx, CS%diag)
  endif

  ! mixing of passive tracers from massless boundary layers to interior
  call cpu_clock_begin(id_clock_tracers)

  if (CS%mix_boundary_tracer_ALE) then
    Tr_ea_BBL = GV%Z_to_H * sqrt(dt*CS%Kd_BBL_tr)
    !$OMP parallel do default(shared) private(htot,in_boundary,add_ent)
    do j=js,je
      do i=is,ie
        htot(i) = 0.0
        in_boundary(i) = (G%mask2dT(i,j) > 0.0)
      enddo
      do k=nz,2,-1 ; do i=is,ie
        if (in_boundary(i)) then
          htot(i) = htot(i) + h(i,j,k)
          !   If diapycnal mixing has been suppressed because this is a massless layer near the
          ! bottom, add some mixing of tracers between these layers.  This flux is based on the
          ! harmonic mean of the two thicknesses, following what is done in layered mode. Kd_min_tr
          ! should be much less than the values in Kd_salt, perhaps a molecular diffusivity.
          add_ent = ((dt * CS%Kd_min_tr) * GV%Z_to_H**2) * &
                    ((h(i,j,k-1)+h(i,j,k) + h_neglect) /  (h(i,j,k-1)*h(i,j,k) + h_neglect2)) - &
                    ent_s(i,j,K)
          if (htot(i) < Tr_ea_BBL) then
            add_ent = max(0.0, add_ent, (Tr_ea_BBL - htot(i)) - ent_s(i,j,K))
          elseif (add_ent < 0.0) then
            add_ent = 0.0 ; in_boundary(i) = .false.
          endif

          ent_s(i,j,K) = ent_s(i,j,K) + add_ent
        endif
      enddo ; enddo
    enddo
  endif  ! (CS%mix_boundary_tracer_ALE)

  ! For passive tracers, the changes in thickness due to boundary fluxes has yet to be applied
  call call_tracer_column_fns(h_orig, h, ent_s(:,:,1:nz), ent_s(:,:,2:nz+1), fluxes, Hml, dt, &
                              G, GV, US, tv, CS%optics, CS%tracer_flow_CSp, CS%debug, &
                              evap_CFL_limit=CS%evap_CFL_limit, &
                              minimum_forcing_depth=CS%minimum_forcing_depth)

  call cpu_clock_end(id_clock_tracers)

  ! Apply ALE sponge
  if (CS%use_sponge .and. associated(CS%ALE_sponge_CSp)) then
    call cpu_clock_begin(id_clock_sponge)
    call apply_ALE_sponge(h, dt, G, GV, US, CS%ALE_sponge_CSp, CS%Time)
    call cpu_clock_end(id_clock_sponge)
    if (CS%debug) then
      call MOM_state_chksum("apply_sponge ", u, v, h, G, GV, US, haloshift=0)
      call MOM_thermovar_chksum("apply_sponge ", tv, G)
    endif
  endif ! CS%use_sponge

  call cpu_clock_begin(id_clock_pass)
  ! visc%Kv_slow is not in the group pass because it has larger vertical extent.
  if (associated(visc%Kv_slow)) &
    call pass_var(visc%Kv_slow, G%Domain, To_All+Omit_Corners, halo=1)
  call cpu_clock_end(id_clock_pass)

  call disable_averaging(CS%diag)

  if (showCallTree) call callTree_leave("diabatic_ALE()")

end subroutine diabatic_ALE

!> Imposes the diapycnal mass fluxes and the accompanying diapycnal advection of momentum and tracers
!! using the original MOM6 algorithms.
subroutine layered_diabatic(u, v, h, tv, Hml, fluxes, visc, ADp, CDp, dt, Time_end, &
                            G, GV, US, CS, Waves)
  type(ocean_grid_type),                      intent(inout) :: G        !< ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV       !< ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US       !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: u        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: v        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h        !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(inout) :: tv       !< points to thermodynamic fields
                                                                        !! unused have NULL ptrs
  real, dimension(:,:),                       pointer       :: Hml      !< Active mixed layer depth [Z ~> m]
  type(forcing),                              intent(inout) :: fluxes   !< points to forcing fields
                                                                        !! unused fields have NULL ptrs
  type(vertvisc_type),                        intent(inout) :: visc     !< vertical viscosities, BBL properies, and
  type(accel_diag_ptrs),                      intent(inout) :: ADp      !< related points to accelerations in momentum
                                                                        !! equations, to enable the later derived
                                                                        !! diagnostics, like energy budgets
  type(cont_diag_ptrs),                       intent(inout) :: CDp      !< points to terms in continuity equations
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
  type(diabatic_CS),                          pointer       :: CS       !< module control structure
  type(Wave_parameters_CS),         optional, pointer       :: Waves    !< Surface gravity waves

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    ea,     &    ! amount of fluid entrained from the layer above within
                 ! one time step [H ~> m or kg m-2]
    eb,     &    ! amount of fluid entrained from the layer below within
                 ! one time step [H ~> m or kg m-2]
    Kd_lay, &    ! diapycnal diffusivity of layers [Z2 T-1 ~> m2 s-1]
    h_orig, &    ! initial layer thicknesses [H ~> m or kg m-2]
    hold,   &    ! layer thickness before diapycnal entrainment, and later the initial
                 ! layer thicknesses (if a mixed layer is used) [H ~> m or kg m-2]
    dSV_dT, &    ! The partial derivative of specific volume with temperature [R-1 degC-1 ~> m3 kg-1 degC-1]
    dSV_dS, &    ! The partial derivative of specific volume with salinity [R-1 ppt-1 ~> m3 kg-1 ppt-1].
    u_h,    &    ! Zonal velocities at thickness points after entrainment [L T-1 ~> m s-1]
    v_h,    &    ! Meridional velocities at thickness points after entrainment [L T-1 ~> m s-1]
    temp_diag, & ! Diagnostic array of previous temperatures [degC]
    saln_diag    ! Diagnostic array of previous salinity [ppt]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    Rcv_ml, &    ! Coordinate density of mixed layer [R ~> kg m-3], used for applying sponges
    SkinBuoyFlux ! 2d surface buoyancy flux [Z2 T-3 ~> m2 s-3], used by ePBL

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target :: &
             ! These are targets so that the space can be shared with eaml & ebml.
    eatr, &  ! The equivalent of ea and eb for tracers, which differ from ea and
    ebtr     ! eb in that they tend to homogenize tracers in massless layers
             ! near the boundaries [H ~> m or kg m-2] (for Bous or non-Bouss)

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    Kd_int,   & ! diapycnal diffusivity of interfaces [Z2 T-1 ~> m2 s-1]
    Kd_heat,  & ! diapycnal diffusivity of heat [Z2 T-1 ~> m2 s-1]
    Kd_salt,  & ! diapycnal diffusivity of salt and passive tracers [Z2 T-1 ~> m2 s-1]
    Kd_extra_T , & ! The extra diffusivity of temperature due to double diffusion relative to
                ! Kd_int [Z2 T-1 ~> m2 s-1].
    Kd_extra_S , & !  The extra diffusivity of salinity due to double diffusion relative to
                ! Kd_int [Z2 T-1 ~> m2 s-1].
    Tdif_flx, & ! diffusive diapycnal heat flux across interfaces [degC H T-1 ~> degC m s-1 or degC kg m-2 s-1]
    Tadv_flx, & ! advective diapycnal heat flux across interfaces [degC H T-1 ~> degC m s-1 or degC kg m-2 s-1]
    Sdif_flx, & ! diffusive diapycnal salt flux across interfaces [ppt H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]
    Sadv_flx    ! advective diapycnal salt flux across interfaces [ppt H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]

  real, allocatable, dimension(:,:) :: &
    hf_dudt_dia_2d, hf_dvdt_dia_2d ! Depth sum of diapycnal mixing accelaration * fract. thickness [L T-2 ~> m s-2].

  ! The following 3 variables are only used with a bulk mixed layer.
  real, pointer, dimension(:,:,:) :: &
    eaml, &  ! The equivalent of ea due to mixed layer processes [H ~> m or kg m-2].
    ebml     ! The equivalent of eb due to mixed layer processes [H ~> m or kg m-2].
             ! eaml and ebml are pointers to eatr and ebtr so as to reuse the memory as
             ! the arrays are not needed at the same time.

  integer :: kb(SZI_(G),SZJ_(G)) ! index of the lightest layer denser
                                 ! than the buffer layer [nondim]

  real :: p_ref_cv(SZI_(G))      ! Reference pressure for the potential density that defines the
                                 ! coordinate variable, set to P_Ref [R L2 T-2 ~> Pa].

  logical :: in_boundary(SZI_(G)) ! True if there are no massive layers below,
                                  ! where massive is defined as sufficiently thick that
                                  ! the no-flux boundary conditions have not restricted
                                  ! the entrainment - usually sqrt(Kd*dt).

  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected [H ~> m or kg m-2]
  real :: h_neglect2   ! h_neglect^2 [H2 ~> m2 or kg2 m-4]
  real :: net_ent      ! The net of ea-eb at an interface [H ~> m or kg m-2]
  real :: add_ent      ! Entrainment that needs to be added when mixing tracers [H ~> m or kg m-2]
  real :: eaval        ! eaval is 2*ea at velocity grid points [H ~> m or kg m-2]
  real :: hval         ! hval is 2*h at velocity grid points [H ~> m or kg m-2]
  real :: h_tr         ! h_tr is h at tracer points with a tiny thickness
                       ! added to ensure positive definiteness [H ~> m or kg m-2]
  real :: Tr_ea_BBL    ! The diffusive tracer thickness in the BBL that is
                       ! coupled to the bottom within a timestep [H ~> m or kg m-2]

  real :: htot(SZIB_(G))        ! The summed thickness from the bottom [H ~> m or kg m-2].
  real :: b1(SZIB_(G))          ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1]
  real :: b_denom_1             ! The first term in the denominator of b1 [H ~> m or kg m-2]
  real :: d1(SZIB_(G))          ! A variable used by the tridiagonal solver [nondim]
  real :: c1(SZIB_(G),SZK_(GV)) ! A variable used by the tridiagonal solver [nondim]

  real :: dt_mix  ! The amount of time over which to apply mixing [T ~> s]
  real :: Idt     ! The inverse time step [T-1 ~> s-1]

  integer :: dir_flag     ! An integer encoding the directions in which to do halo updates.
  logical :: showCallTree ! If true, show the call tree
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb, m, halo

  is   = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq  = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nkmb = GV%nk_rho_varies
  h_neglect = GV%H_subroundoff ; h_neglect2 = h_neglect*h_neglect
  Kd_heat(:,:,:) = 0.0 ; Kd_salt(:,:,:) = 0.0


  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("layered_diabatic(), MOM_diabatic_driver.F90")

  ! set equivalence between the same bits of memory for these arrays
  eaml => eatr ; ebml => ebtr

  ! For all other diabatic subroutines, the averaging window should be the entire diabatic timestep
  call enable_averages(dt, Time_end, CS%diag)

  if ((CS%ML_mix_first > 0.0) .or. CS%use_geothermal) then
    halo = CS%halo_TS_diff
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=js-halo,je+halo ; do i=is-halo,ie+halo
      h_orig(i,j,k) = h(i,j,k) ; eaml(i,j,k) = 0.0 ; ebml(i,j,k) = 0.0
    enddo ; enddo ; enddo
  endif

  if (CS%use_geothermal) then
    call cpu_clock_begin(id_clock_geothermal)
    call geothermal_entraining(h, tv, dt, eaml, ebml, G, GV, US, CS%geothermal_CSp, halo=CS%halo_TS_diff)
    call cpu_clock_end(id_clock_geothermal)
    if (showCallTree) call callTree_waypoint("geothermal (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('geothermal', u, v, h, tv%T, tv%S, G, GV, US)
  endif

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  ! Set_pen_shortwave estimates the optical properties of the water column.
  ! It will need to be modified later to include information about the
  ! biological properties and layer thicknesses.
  if (associated(CS%optics)) &
    call set_pen_shortwave(CS%optics, fluxes, G, GV, US, CS%diabatic_aux_CSp, CS%opacity_CSp, CS%tracer_flow_CSp)

  if (CS%bulkmixedlayer) then
    if (CS%debug) call MOM_forcing_chksum("Before mixedlayer", fluxes, G, US, haloshift=0)

    if (CS%ML_mix_first > 0.0) then
!  This subroutine
!    (1) Cools the mixed layer.
!    (2) Performs convective adjustment by mixed layer entrainment.
!    (3) Heats the mixed layer and causes it to detrain to
!        Monin-Obukhov depth or minimum mixed layer depth.
!    (4) Uses any remaining TKE to drive mixed layer entrainment.
!    (5) Possibly splits buffer layer into two isopycnal layers (when using isopycnal coordinate)
      call find_uv_at_h(u, v, h, u_h, v_h, G, GV, US)

      call cpu_clock_begin(id_clock_mixedlayer)
      if (CS%ML_mix_first < 1.0) then
        ! Changes: h, tv%T, tv%S, eaml and ebml  (G is also inout???)
        call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt*CS%ML_mix_first, &
                            eaml,ebml, G, GV, US, CS%bulkmixedlayer_CSp, CS%optics, &
                            Hml, CS%aggregate_FW_forcing, dt, last_call=.false.)
        ! Changes: h, tv%T, tv%S, eaml and ebml  (G is also inout???)
        call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt, eaml, ebml, &
                            G, GV, US, CS%bulkmixedlayer_CSp, CS%optics, &
                            Hml, CS%aggregate_FW_forcing, dt, last_call=.true.)
      endif

      !  Keep salinity from falling below a small but positive threshold.
      !  This constraint is needed for SIS1 ice model, which can extract
      !  more salt than is present in the ocean. SIS2 does not suffer
      !  from this limitation, in which case we can let salinity=0 and still
      !  have salt conserved with SIS2 ice. So for SIS2, we can run with
      !  BOUND_SALINITY=False in MOM.F90.
      if (associated(tv%S) .and. associated(tv%salt_deficit)) &
        call adjust_salt(h, tv, G, GV, CS%diabatic_aux_CSp)
      call cpu_clock_end(id_clock_mixedlayer)
      if (CS%debug) then
        call MOM_state_chksum("After mixedlayer ", u, v, h, G, GV, US, haloshift=0)
        call MOM_forcing_chksum("After mixedlayer", fluxes, G, US, haloshift=0)
      endif
      if (showCallTree) call callTree_waypoint("done with 1st bulkmixedlayer (diabatic)")
      if (CS%debugConservation) call MOM_state_stats('1st bulkmixedlayer', u, v, h, tv%T, tv%S, G, GV, US)
    endif
  endif

  if (CS%debug) &
    call MOM_state_chksum("before find_uv_at_h", u, v, h, G, GV, US, haloshift=0)
  if (CS%use_kappa_shear .or. CS%use_CVMix_shear) then
    if ((CS%ML_mix_first > 0.0) .or. CS%use_geothermal) then
      call find_uv_at_h(u, v, h_orig, u_h, v_h, G, GV, US, eaml, ebml)
      if (CS%debug) then
        call hchksum(eaml, "after find_uv_at_h eaml", G%HI, scale=GV%H_to_m)
        call hchksum(ebml, "after find_uv_at_h ebml", G%HI, scale=GV%H_to_m)
      endif
    else
      call find_uv_at_h(u, v, h, u_h, v_h, G, GV, US)
    endif
    if (showCallTree) call callTree_waypoint("done with find_uv_at_h (diabatic)")
  endif

  call cpu_clock_begin(id_clock_set_diffusivity)
  ! Sets: Kd_lay, Kd_int, Kd_extra_T, Kd_extra_S and visc%TKE_turb
  ! Also changes: visc%Kd_shear and visc%Kv_shear
  if ((CS%halo_TS_diff > 0) .and. (CS%ML_mix_first > 0.0)) then
    if (associated(tv%T)) call pass_var(tv%T, G%Domain, halo=CS%halo_TS_diff, complete=.false.)
    if (associated(tv%T)) call pass_var(tv%S, G%Domain, halo=CS%halo_TS_diff, complete=.false.)
    call pass_var(h, G%domain, halo=CS%halo_TS_diff, complete=.true.)
  endif
  if (CS%debug) &
    call MOM_state_chksum("before set_diffusivity", u, v, h, G, GV, US, haloshift=CS%halo_TS_diff)
  if (CS%double_diffuse) then
    call set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, CS%optics, visc, dt, G, GV, US, CS%set_diff_CSp, &
                         Kd_lay=Kd_lay, Kd_int=Kd_int, Kd_extra_T=Kd_extra_T, Kd_extra_S=Kd_extra_S)
  else
    call set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, CS%optics, visc, dt, G, GV, US, &
                         CS%set_diff_CSp, Kd_lay=Kd_lay, Kd_int=Kd_int)
  endif
  call cpu_clock_end(id_clock_set_diffusivity)
  if (showCallTree) call callTree_waypoint("done with set_diffusivity (diabatic)")

  if (CS%debug) then
    call MOM_state_chksum("after set_diffusivity ", u, v, h, G, GV, US, haloshift=0)
    call MOM_forcing_chksum("after set_diffusivity ", fluxes, G, US, haloshift=0)
    call MOM_thermovar_chksum("after set_diffusivity ", tv, G)
    call hchksum(Kd_lay, "after set_diffusivity Kd_lay", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
    call hchksum(Kd_Int, "after set_diffusivity Kd_Int", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
  endif


  if (CS%useKPP) then
    call cpu_clock_begin(id_clock_kpp)
    ! KPP needs the surface buoyancy flux but does not update state variables.
    ! We could make this call higher up to avoid a repeat unpacking of the surface fluxes.
    ! Sets: CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux
    ! NOTE: CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux are returned as rates (i.e. stuff per second)
    ! unlike other instances where the fluxes are integrated in time over a time-step.
    call calculateBuoyancyFlux2d(G, GV, US, fluxes, CS%optics, h, tv%T, tv%S, tv, &
                                 CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux)
    ! The KPP scheme calculates boundary layer diffusivities and non-local transport.

    ! Set diffusivities for heat and salt separately

    if (CS%double_diffuse) then
      ! Add contribution from double diffusion
      !$OMP parallel do default(shared)
      do K=1,nz+1 ; do j=js,je ; do i=is,ie
        Kd_salt(i,j,K) = Kd_int(i,j,K) + Kd_extra_S(i,j,K)
        Kd_heat(i,j,K) = Kd_int(i,j,K) + Kd_extra_T(i,j,K)
      enddo ; enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do K=1,nz+1 ; do j=js,je ; do i=is,ie
        Kd_salt(i,j,K) = Kd_int(i,j,K)
        Kd_heat(i,j,K) = Kd_int(i,j,K)
      enddo ; enddo ; enddo
    endif

    call KPP_compute_BLD(CS%KPP_CSp, G, GV, US, h, tv%T, tv%S, u, v, tv, &
                         fluxes%ustar, CS%KPP_buoy_flux, Waves=Waves)

    call KPP_calculate(CS%KPP_CSp, G, GV, US, h, fluxes%ustar, CS%KPP_buoy_flux, Kd_heat, &
                       Kd_salt, visc%Kv_shear, CS%KPP_NLTheat, CS%KPP_NLTscalar, Waves=Waves)

    if (associated(Hml)) then
      call KPP_get_BLD(CS%KPP_CSp, Hml(:,:), G, US)
      call pass_var(Hml, G%domain, halo=1)
      ! If visc%MLD exists, copy KPP's BLD into it
      if (associated(visc%MLD)) visc%MLD(:,:) = Hml(:,:)
    endif

    if (.not. CS%KPPisPassive) then
      !$OMP parallel do default(shared)
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        Kd_int(i,j,K) = min( Kd_salt(i,j,k),  Kd_heat(i,j,k) )
      enddo ; enddo ; enddo
      if (CS%double_diffuse) then
        !$OMP parallel do default(shared)
        do k=1,nz+1 ; do j=js,je ; do i=is,ie
          Kd_extra_S(i,j,k) = (Kd_salt(i,j,k) - Kd_int(i,j,K))
          Kd_extra_T(i,j,k) = (Kd_heat(i,j,k) - Kd_int(i,j,K))
        enddo ; enddo ; enddo
      endif
    endif ! not passive

    call cpu_clock_end(id_clock_kpp)
    if (showCallTree) call callTree_waypoint("done with KPP_calculate (diabatic)")
    if (CS%debug) then
      call MOM_state_chksum("after KPP", u, v, h, G, GV, US, haloshift=0)
      call MOM_forcing_chksum("after KPP", fluxes, G, US, haloshift=0)
      call MOM_thermovar_chksum("after KPP", tv, G)
      call hchksum(Kd_lay, "after KPP Kd_lay", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
      call hchksum(Kd_Int, "after KPP Kd_Int", G%HI, haloshift=0, scale=US%Z2_T_to_m2_s)
    endif

  endif  ! endif for KPP

  ! Add vertical diff./visc. due to convection (computed via CVMix)
  if (CS%use_CVMix_conv) then
    call calculate_CVMix_conv(h, tv, G, GV, US, CS%CVMix_conv_csp, Hml, Kd=Kd_int, Kv=visc%Kv_slow)
  endif

  if (CS%useKPP) then
    call cpu_clock_begin(id_clock_kpp)
    if (CS%debug) then
      call hchksum(CS%KPP_temp_flux, "before KPP_applyNLT netHeat", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(CS%KPP_salt_flux, "before KPP_applyNLT netSalt", G%HI, haloshift=0, scale=GV%H_to_m)
      call hchksum(CS%KPP_NLTheat, "before KPP_applyNLT NLTheat", G%HI, haloshift=0)
      call hchksum(CS%KPP_NLTscalar, "before KPP_applyNLT NLTscalar", G%HI, haloshift=0)
    endif
    ! Apply non-local transport of heat and salt
    ! Changes: tv%T, tv%S
    call KPP_NonLocalTransport_temp(CS%KPP_CSp, G, GV, h, CS%KPP_NLTheat,   CS%KPP_temp_flux, &
                                    US%T_to_s*dt, tv%T, US%Q_to_J_kg*tv%C_p)
    call KPP_NonLocalTransport_saln(CS%KPP_CSp, G, GV, h, CS%KPP_NLTscalar, CS%KPP_salt_flux, &
                                    US%T_to_s*dt, tv%S)
    call cpu_clock_end(id_clock_kpp)
    if (showCallTree) call callTree_waypoint("done with KPP_applyNonLocalTransport (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('KPP_applyNonLocalTransport', u, v, h, tv%T, tv%S, G, GV, US)

    if (CS%debug) then
      call MOM_state_chksum("after KPP_applyNLT ", u, v, h, G, GV, US, haloshift=0)
      call MOM_forcing_chksum("after KPP_applyNLT ", fluxes, G, US, haloshift=0)
      call MOM_thermovar_chksum("after KPP_applyNLT ", tv, G)
    endif
  endif ! endif for KPP

  ! Differential diffusion done here.
  ! Changes: tv%T, tv%S
  if (CS%double_diffuse .and. associated(tv%T)) then

    call cpu_clock_begin(id_clock_differential_diff)
    call differential_diffuse_T_S(h, tv%T, tv%S, Kd_extra_T, Kd_extra_S, dt, G, GV)
    call cpu_clock_end(id_clock_differential_diff)
    if (showCallTree) call callTree_waypoint("done with differential_diffuse_T_S (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('differential_diffuse_T_S', u, v, h, tv%T, tv%S, G, GV, US)

    ! increment heat and salt diffusivity.
    ! CS%useKPP==.true. already has extra_T and extra_S included
    if (.not. CS%useKPP) then
      !$OMP parallel do default(shared)
      do K=2,nz ; do j=js,je ; do i=is,ie
        Kd_heat(i,j,K) = Kd_heat(i,j,K) + Kd_extra_T(i,j,K)
        Kd_salt(i,j,K) = Kd_salt(i,j,K) + Kd_extra_S(i,j,K)
      enddo ; enddo ; enddo
    endif

  endif

  ! Calculate layer entrainments and detrainments from diffusivities and differences between
  ! layer and target densities (i.e. do remapping as well as diffusion).
  call cpu_clock_begin(id_clock_entrain)
  ! Calculate appropriately limited diapycnal mass fluxes to account
  ! for diapycnal diffusion and advection.  Sets: ea, eb. Changes: kb
  call Entrainment_diffusive(h, tv, fluxes, dt, G, GV, US, CS%entrain_diffusive_CSp, &
                             ea, eb, kb, Kd_lay=Kd_lay, Kd_int=Kd_int)
  call cpu_clock_end(id_clock_entrain)
  if (showCallTree) call callTree_waypoint("done with Entrainment_diffusive (diabatic)")

  if (CS%debug) then
    call MOM_forcing_chksum("after calc_entrain ", fluxes, G, US, haloshift=0)
    call MOM_thermovar_chksum("after calc_entrain ", tv, G)
    call MOM_state_chksum("after calc_entrain ", u, v, h, G, GV, US, haloshift=0)
    call hchksum(ea, "after calc_entrain ea", G%HI, haloshift=0, scale=GV%H_to_m)
    call hchksum(eb, "after calc_entrain eb", G%HI, haloshift=0, scale=GV%H_to_m)
  endif

  ! Save fields before boundary forcing is applied for tendency diagnostics
  if (CS%boundary_forcing_tendency_diag) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      temp_diag(i,j,k) = tv%T(i,j,k)
      saln_diag(i,j,k) = tv%S(i,j,k)
    enddo ; enddo ; enddo
  endif

  ! Update h according to divergence of the difference between
  ! ea and eb. We keep a record of the original h in hold.
  ! In the following, the checks for negative values are to guard
  ! against instances where entrainment drives a layer to
  ! negative thickness.  This situation will never happen if
  ! enough iterations are permitted in Calculate_Entrainment.
  ! Even if too few iterations are allowed, it is still guarded
  ! against.  In other words the checks are probably unnecessary.
  !$OMP parallel do default(shared)
  do j=js,je
    do i=is,ie
      hold(i,j,1) = h(i,j,1)
      h(i,j,1) = h(i,j,1) + (eb(i,j,1) - ea(i,j,2))
      hold(i,j,nz) = h(i,j,nz)
      h(i,j,nz) = h(i,j,nz) + (ea(i,j,nz) - eb(i,j,nz-1))
      if (h(i,j,1) <= 0.0) then
        h(i,j,1) = GV%Angstrom_H
      endif
      if (h(i,j,nz) <= 0.0) then
        h(i,j,nz) = GV%Angstrom_H
      endif
    enddo
    do k=2,nz-1 ; do i=is,ie
      hold(i,j,k) = h(i,j,k)
      h(i,j,k) = h(i,j,k) + ((ea(i,j,k) - eb(i,j,k-1)) + &
                    (eb(i,j,k) - ea(i,j,k+1)))
      if (h(i,j,k) <= 0.0) then
        h(i,j,k) = GV%Angstrom_H
      endif
    enddo ; enddo
  enddo
  ! Checks for negative thickness may have changed layer thicknesses
  call diag_update_remap_grids(CS%diag)

  if (CS%debug) then
    call MOM_state_chksum("after negative check ", u, v, h, G, GV, US, haloshift=0)
    call MOM_forcing_chksum("after negative check ", fluxes, G, US, haloshift=0)
    call MOM_thermovar_chksum("after negative check ", tv, G)
  endif
  if (showCallTree) call callTree_waypoint("done with h=ea-eb (diabatic)")
  if (CS%debugConservation) call MOM_state_stats('h=ea-eb', u, v, h, tv%T, tv%S, G, GV, US)

  ! Here, T and S are updated according to ea and eb.
  ! If using the bulk mixed layer, T and S are also updated
  ! by surface fluxes (in fluxes%*).
  ! This is a very long block.
  if (CS%bulkmixedlayer) then

    if (associated(tv%T)) then
      call cpu_clock_begin(id_clock_tridiag)
      ! Temperature and salinity (as state variables) are treated
      ! differently from other tracers to insure massless layers that
      ! are lighter than the mixed layer have temperatures and salinities
      ! that correspond to their prescribed densities.
      if (CS%massless_match_targets) then
        !$OMP parallel do default (shared) private(h_tr,b1,d1,c1,b_denom_1)
        do j=js,je
          do i=is,ie
            h_tr = hold(i,j,1) + h_neglect
            b1(i) = 1.0 / (h_tr + eb(i,j,1))
            d1(i) = h_tr * b1(i)
            tv%T(i,j,1) = b1(i) * (h_tr*tv%T(i,j,1))
            tv%S(i,j,1) = b1(i) * (h_tr*tv%S(i,j,1))
          enddo
          do k=2,nkmb ; do i=is,ie
            c1(i,k) = eb(i,j,k-1) * b1(i)
            h_tr = hold(i,j,k) + h_neglect
            b_denom_1 = h_tr + d1(i)*ea(i,j,k)
            b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
            if (k<nkmb) d1(i) = b_denom_1 * b1(i)
            tv%T(i,j,k) = b1(i) * (h_tr*tv%T(i,j,k) + ea(i,j,k)*tv%T(i,j,k-1))
            tv%S(i,j,k) = b1(i) * (h_tr*tv%S(i,j,k) + ea(i,j,k)*tv%S(i,j,k-1))
          enddo ; enddo

          do k=nkmb+1,nz ; do i=is,ie
            if (k == kb(i,j)) then
              c1(i,k) = eb(i,j,k-1) * b1(i)
              d1(i) = (((eb(i,j,nkmb)-eb(i,j,k-1)) + hold(i,j,nkmb) + h_neglect) + &
                       d1(i)*ea(i,j,nkmb)) * b1(i)
              h_tr = hold(i,j,k) + h_neglect
              b_denom_1 = h_tr + d1(i)*ea(i,j,k)
              b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
              d1(i) = b_denom_1 * b1(i)
              tv%T(i,j,k) = b1(i) * (h_tr*tv%T(i,j,k) + ea(i,j,k)*tv%T(i,j,nkmb))
              tv%S(i,j,k) = b1(i) * (h_tr*tv%S(i,j,k) + ea(i,j,k)*tv%S(i,j,nkmb))
            elseif (k > kb(i,j)) then
              c1(i,k) = eb(i,j,k-1) * b1(i)
              h_tr = hold(i,j,k) + h_neglect
              b_denom_1 = h_tr + d1(i)*ea(i,j,k)
              b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
              d1(i) = b_denom_1 * b1(i)
              tv%T(i,j,k) = b1(i) * (h_tr*tv%T(i,j,k) + ea(i,j,k)*tv%T(i,j,k-1))
              tv%S(i,j,k) = b1(i) * (h_tr*tv%S(i,j,k) + ea(i,j,k)*tv%S(i,j,k-1))
            elseif (eb(i,j,k) < eb(i,j,k-1)) then ! (note that k < kb(i,j))
              !   The bottommost buffer layer might entrain all the mass from some
              ! of the interior layers that are thin and lighter in the coordinate
              ! density than that buffer layer.  The T and S of these newly
              ! massless interior layers are unchanged.
              tv%T(i,j,nkmb) = tv%T(i,j,nkmb) + b1(i) * (eb(i,j,k-1) - eb(i,j,k)) * tv%T(i,j,k)
              tv%S(i,j,nkmb) = tv%S(i,j,nkmb) + b1(i) * (eb(i,j,k-1) - eb(i,j,k)) * tv%S(i,j,k)
            endif
          enddo ; enddo

          do k=nz-1,nkmb,-1 ; do i=is,ie
            if (k >= kb(i,j)) then
              tv%T(i,j,k) = tv%T(i,j,k) + c1(i,k+1)*tv%T(i,j,k+1)
              tv%S(i,j,k) = tv%S(i,j,k) + c1(i,k+1)*tv%S(i,j,k+1)
            endif
          enddo ; enddo
          do i=is,ie ; if (kb(i,j) <= nz) then
            tv%T(i,j,nkmb) = tv%T(i,j,nkmb) + c1(i,kb(i,j))*tv%T(i,j,kb(i,j))
            tv%S(i,j,nkmb) = tv%S(i,j,nkmb) + c1(i,kb(i,j))*tv%S(i,j,kb(i,j))
          endif ; enddo
          do k=nkmb-1,1,-1 ; do i=is,ie
            tv%T(i,j,k) = tv%T(i,j,k) + c1(i,k+1)*tv%T(i,j,k+1)
            tv%S(i,j,k) = tv%S(i,j,k) + c1(i,k+1)*tv%S(i,j,k+1)
          enddo ; enddo
        enddo ! end of j loop
      else ! .not. massless_match_targets
        ! This simpler form allows T & S to be too dense for the layers
        ! between the buffer layers and the interior.
        ! Changes: T, S
        if (CS%tracer_tridiag) then
          call tracer_vertdiff(hold, ea, eb, dt, tv%T, G, GV)
          call tracer_vertdiff(hold, ea, eb, dt, tv%S, G, GV)
        else
          call triDiagTS(G, GV, is, ie, js, je, hold, ea, eb, tv%T, tv%S)
        endif
      endif ! massless_match_targets
      call cpu_clock_end(id_clock_tridiag)

    endif ! endif for associated(T)
    if (CS%debugConservation) call MOM_state_stats('BML tridiag', u, v, h, tv%T, tv%S, G, GV, US)

    if ((CS%ML_mix_first > 0.0) .or. CS%use_geothermal) then
      ! The mixed layer code has already been called, but there is some needed
      ! bookkeeping.
      !$OMP parallel do default(shared)
      do k=1,nz ; do j=js,je ; do i=is,ie
        hold(i,j,k) = h_orig(i,j,k)
        ea(i,j,k) = ea(i,j,k) + eaml(i,j,k)
        eb(i,j,k) = eb(i,j,k) + ebml(i,j,k)
      enddo ; enddo ; enddo
      if (CS%debug) then
        call hchksum(ea, "after ea = ea + eaml", G%HI, haloshift=0, scale=GV%H_to_m)
        call hchksum(eb, "after eb = eb + ebml", G%HI, haloshift=0, scale=GV%H_to_m)
      endif
    endif

    if (CS%ML_mix_first < 1.0) then
    !  Call the mixed layer code now, perhaps for a second time.
    !  This subroutine (1)  Cools the mixed layer.
    !    (2) Performs convective adjustment by mixed layer entrainment.
    !    (3) Heats the mixed layer and causes it to detrain to
    !        Monin-Obukhov depth or minimum mixed layer depth.
    !    (4) Uses any remaining TKE to drive mixed layer entrainment.
    !    (5) Possibly splits the buffer layer into two isopycnal layers.

      call find_uv_at_h(u, v, hold, u_h, v_h, G, GV, US, ea, eb)
      if (CS%debug) call MOM_state_chksum("find_uv_at_h1 ", u, v, h, G, GV, US, haloshift=0)

      dt_mix = min(dt, dt*(1.0 - CS%ML_mix_first))
      call cpu_clock_begin(id_clock_mixedlayer)
      ! Changes: h, tv%T, tv%S, ea and eb  (G is also inout???)
      call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt_mix, ea, eb, &
                          G, GV, US, CS%bulkmixedlayer_CSp, CS%optics, &
                          Hml, CS%aggregate_FW_forcing, dt, last_call=.true.)

      !  Keep salinity from falling below a small but positive threshold.
      !  This constraint is needed for SIS1 ice model, which can extract
      !  more salt than is present in the ocean. SIS2 does not suffer
      !  from this limitation, in which case we can let salinity=0 and still
      !  have salt conserved with SIS2 ice. So for SIS2, we can run with
      !  BOUND_SALINITY=False in MOM.F90.
      if (associated(tv%S) .and. associated(tv%salt_deficit)) &
        call adjust_salt(h, tv, G, GV, CS%diabatic_aux_CSp)

      call cpu_clock_end(id_clock_mixedlayer)
      if (showCallTree) call callTree_waypoint("done with 2nd bulkmixedlayer (diabatic)")
      if (CS%debugConservation) call MOM_state_stats('2nd bulkmixedlayer', u, v, h, tv%T, tv%S, G, GV, US)
    endif

  else  ! following block for when NOT using BULKMIXEDLAYER

    ! calculate change in temperature & salinity due to dia-coordinate surface diffusion
    if (associated(tv%T)) then

      if (CS%debug) then
        call hchksum(ea, "before triDiagTS ea ", G%HI, haloshift=0, scale=GV%H_to_m)
        call hchksum(eb, "before triDiagTS eb ", G%HI, haloshift=0, scale=GV%H_to_m)
      endif
      call cpu_clock_begin(id_clock_tridiag)

      !  Keep salinity from falling below a small but positive threshold.
      !  This constraint is needed for SIS1 ice model, which can extract
      !  more salt than is present in the ocean. SIS2 does not suffer
      !  from this limitation, in which case we can let salinity=0 and still
      !  have salt conserved with SIS2 ice. So for SIS2, we can run with
      !  BOUND_SALINITY=False in MOM.F90.
      if (associated(tv%S) .and. associated(tv%salt_deficit)) &
        call adjust_salt(h, tv, G, GV, CS%diabatic_aux_CSp)

      if (CS%diabatic_diff_tendency_diag) then
        do k=1,nz ; do j=js,je ; do i=is,ie
          temp_diag(i,j,k) = tv%T(i,j,k)
          saln_diag(i,j,k) = tv%S(i,j,k)
        enddo ; enddo ; enddo
      endif

      ! Changes T and S via the tridiagonal solver; no change to h
      if (CS%tracer_tridiag) then
        call tracer_vertdiff(hold, ea, eb, dt, tv%T, G, GV)
        call tracer_vertdiff(hold, ea, eb, dt, tv%S, G, GV)
      else
        call triDiagTS(G, GV, is, ie, js, je, hold, ea, eb, tv%T, tv%S)
      endif

      ! diagnose temperature, salinity, heat, and salt tendencies
      ! Note: hold here refers to the thicknesses from before the dual-entraintment when using
      ! the bulk mixed layer scheme, so tendencies should be posted on hold.
      if (CS%diabatic_diff_tendency_diag) then
        call diagnose_diabatic_diff_tendency(tv, hold, temp_diag, saln_diag, dt, G, GV, US, CS)
        if (CS%id_diabatic_diff_h > 0) call post_data(CS%id_diabatic_diff_h, hold, CS%diag, alt_h=hold)
      endif

      call cpu_clock_end(id_clock_tridiag)
      if (showCallTree) call callTree_waypoint("done with triDiagTS (diabatic)")

    endif  ! endif corresponding to if (associated(tv%T))
    if (CS%debugConservation) call MOM_state_stats('triDiagTS', u, v, h, tv%T, tv%S, G, GV, US)

  endif  ! endif for the BULKMIXEDLAYER block

  if (CS%debug) then
    call MOM_state_chksum("after mixed layer ", u, v, h, G, GV, US, haloshift=0)
    call MOM_thermovar_chksum("after mixed layer ", tv, G)
    call hchksum(ea, "after mixed layer ea", G%HI, scale=GV%H_to_m)
    call hchksum(eb, "after mixed layer eb", G%HI, scale=GV%H_to_m)
  endif

  call cpu_clock_begin(id_clock_remap)
  call regularize_layers(h, tv, dt, ea, eb, G, GV, US, CS%regularize_layers_CSp)
  call cpu_clock_end(id_clock_remap)
  if (showCallTree) call callTree_waypoint("done with regularize_layers (diabatic)")
  if (CS%debugConservation) call MOM_state_stats('regularize_layers', u, v, h, tv%T, tv%S, G, GV, US)

  ! Whenever thickness changes let the diag manager know, as the
  ! target grids for vertical remapping may need to be regenerated.
  if (associated(ADp%du_dt_dia) .or. associated(ADp%dv_dt_dia)) &
    ! Remapped d[uv]dt_dia require east/north halo updates of h
    call pass_var(h, G%domain, To_West+To_South+Omit_Corners, halo=1)
  call diag_update_remap_grids(CS%diag)

  ! diagnostics
  Idt = 1.0 / dt
  if ((CS%id_Tdif > 0) .or. (CS%id_Tadv > 0)) then
    do j=js,je ; do i=is,ie
      Tdif_flx(i,j,1) = 0.0 ; Tdif_flx(i,j,nz+1) = 0.0
      Tadv_flx(i,j,1) = 0.0 ; Tadv_flx(i,j,nz+1) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do K=2,nz ; do j=js,je ; do i=is,ie
      Tdif_flx(i,j,K) = (Idt * 0.5*(ea(i,j,k) + eb(i,j,k-1))) * &
                        (tv%T(i,j,k-1) - tv%T(i,j,k))
      Tadv_flx(i,j,K) = (Idt * (ea(i,j,k) - eb(i,j,k-1))) * &
                    0.5*(tv%T(i,j,k-1) + tv%T(i,j,k))
    enddo ; enddo ; enddo
  endif
  if ((CS%id_Sdif > 0) .or. (CS%id_Sadv > 0)) then
    do j=js,je ; do i=is,ie
      Sdif_flx(i,j,1) = 0.0 ; Sdif_flx(i,j,nz+1) = 0.0
      Sadv_flx(i,j,1) = 0.0 ; Sadv_flx(i,j,nz+1) = 0.0
    enddo ; enddo
    !$OMP parallel do default(shared)
    do K=2,nz ; do j=js,je ; do i=is,ie
      Sdif_flx(i,j,K) = (Idt * 0.5*(ea(i,j,k) + eb(i,j,k-1))) * &
                        (tv%S(i,j,k-1) - tv%S(i,j,k))
      Sadv_flx(i,j,K) = (Idt * (ea(i,j,k) - eb(i,j,k-1))) * &
                    0.5*(tv%S(i,j,k-1) + tv%S(i,j,k))
    enddo ; enddo ; enddo
  endif

  ! mixing of passive tracers from massless boundary layers to interior
  call cpu_clock_begin(id_clock_tracers)
  if (CS%mix_boundary_tracers) then
    Tr_ea_BBL = GV%Z_to_H * sqrt(dt*CS%Kd_BBL_tr)
    !$OMP parallel do default(shared) private(htot,in_boundary,add_ent)
    do j=js,je
      do i=is,ie
        ebtr(i,j,nz) = eb(i,j,nz)
        htot(i) = 0.0
        in_boundary(i) = (G%mask2dT(i,j) > 0.0)
      enddo
      do k=nz,2,-1 ; do i=is,ie
        if (in_boundary(i)) then
          htot(i) = htot(i) + h(i,j,k)
          !   If diapycnal mixing has been suppressed because this is a massless
          ! layer near the bottom, add some mixing of tracers between these
          ! layers.  This flux is based on the harmonic mean of the two
          ! thicknesses, as this corresponds pretty closely (to within
          ! differences in the density jumps between layers) with what is done
          ! in the calculation of the fluxes in the first place.  Kd_min_tr
          ! should be much less than the values that have been set in Kd_lay,
          ! perhaps a molecular diffusivity.
          add_ent = ((dt * CS%Kd_min_tr) * GV%Z_to_H**2) * &
                    ((h(i,j,k-1)+h(i,j,k)+h_neglect) / &
                     (h(i,j,k-1)*h(i,j,k)+h_neglect2)) - &
                    0.5*(ea(i,j,k) + eb(i,j,k-1))
          if (htot(i) < Tr_ea_BBL) then
            add_ent = max(0.0, add_ent, &
                          (Tr_ea_BBL - htot(i)) - min(ea(i,j,k), eb(i,j,k-1)))
          elseif (add_ent < 0.0) then
            add_ent = 0.0 ; in_boundary(i) = .false.
          endif

          ebtr(i,j,k-1) = eb(i,j,k-1) + add_ent
          eatr(i,j,k) = ea(i,j,k) + add_ent
        else
          ebtr(i,j,k-1) = eb(i,j,k-1) ; eatr(i,j,k) = ea(i,j,k)
        endif
        if (CS%double_diffuse) then ; if (Kd_extra_S(i,j,K) > 0.0) then
          add_ent = ((dt * Kd_extra_S(i,j,K)) * GV%Z_to_H**2) / &
             (0.25 * ((h(i,j,k-1) + h(i,j,k)) + (hold(i,j,k-1) + hold(i,j,k))) + &
              h_neglect)
          ebtr(i,j,k-1) = ebtr(i,j,k-1) + add_ent
          eatr(i,j,k) = eatr(i,j,k) + add_ent
        endif ; endif
      enddo ; enddo
      do i=is,ie ; eatr(i,j,1) = ea(i,j,1) ; enddo

    enddo

    call call_tracer_column_fns(hold, h, eatr, ebtr, fluxes, Hml, dt, G, GV, US, tv, &
                              CS%optics, CS%tracer_flow_CSp, CS%debug)

  elseif (CS%double_diffuse) then  ! extra diffusivity for passive tracers

    do j=js,je ; do i=is,ie
      ebtr(i,j,nz) = eb(i,j,nz) ; eatr(i,j,1) = ea(i,j,1)
    enddo ; enddo
    !$OMP parallel do default(shared) private(add_ent)
    do k=nz,2,-1 ; do j=js,je ; do i=is,ie
      if (Kd_extra_S(i,j,K) > 0.0) then
        add_ent = ((dt * Kd_extra_S(i,j,K)) * GV%Z_to_H**2) / &
           (0.25 * ((h(i,j,k-1) + h(i,j,k)) + (hold(i,j,k-1) + hold(i,j,k))) + &
            h_neglect)
      else
        add_ent = 0.0
      endif
      ebtr(i,j,k-1) = eb(i,j,k-1) + add_ent
      eatr(i,j,k) = ea(i,j,k) + add_ent
    enddo ; enddo ; enddo

    call call_tracer_column_fns(hold, h, eatr, ebtr, fluxes, Hml, dt, G, GV, US, tv, &
                                CS%optics, CS%tracer_flow_CSp, CS%debug)

  else
    call call_tracer_column_fns(hold, h, ea, eb, fluxes, Hml, dt, G, GV, US, tv, &
                                CS%optics, CS%tracer_flow_CSp, CS%debug)

  endif  ! (CS%mix_boundary_tracers)

  call cpu_clock_end(id_clock_tracers)

  ! sponges
  if (CS%use_sponge) then
    call cpu_clock_begin(id_clock_sponge)
    ! Layer mode sponge
    if (CS%bulkmixedlayer .and. associated(tv%eqn_of_state)) then
      do i=is,ie ; p_ref_cv(i) = tv%P_Ref ; enddo
      EOSdom(:) = EOS_domain(G%HI)
      !$OMP parallel do default(shared)
      do j=js,je
        call calculate_density(tv%T(:,j,1), tv%S(:,j,1), p_ref_cv, Rcv_ml(:,j), &
                               tv%eqn_of_state, EOSdom)
      enddo
      call apply_sponge(h, dt, G, GV, US, ea, eb, CS%sponge_CSp, Rcv_ml)
    else
      call apply_sponge(h, dt, G, GV, US, ea, eb, CS%sponge_CSp)
    endif
    call cpu_clock_end(id_clock_sponge)
    if (CS%debug) then
      call MOM_state_chksum("apply_sponge ", u, v, h, G, GV, US, haloshift=0)
      call MOM_thermovar_chksum("apply_sponge ", tv, G)
    endif
  endif ! CS%use_sponge

!   Save the diapycnal mass fluxes as a diagnostic field.
  if (associated(CDp%diapyc_vel)) then
    !$OMP parallel do default(shared)
    do j=js,je
      do K=2,nz ; do i=is,ie
        CDp%diapyc_vel(i,j,K) = US%s_to_T*Idt * (ea(i,j,k) - eb(i,j,k-1))
      enddo ; enddo
      do i=is,ie
        CDp%diapyc_vel(i,j,1) = 0.0
        CDp%diapyc_vel(i,j,nz+1) = 0.0
      enddo
    enddo
  endif

! For momentum, it is only the net flux that homogenizes within
! the mixed layer.  Vertical viscosity that is proportional to the
! mixed layer turbulence is applied elsewhere.
  if (CS%bulkmixedlayer) then
    if (CS%debug) then
      call hchksum(ea, "before net flux rearrangement ea", G%HI, scale=GV%H_to_m)
      call hchksum(eb, "before net flux rearrangement eb", G%HI, scale=GV%H_to_m)
    endif
    !$OMP parallel do default(shared) private(net_ent)
    do j=js,je
      do K=2,GV%nkml ; do i=is,ie
        net_ent = ea(i,j,k) - eb(i,j,k-1)
        ea(i,j,k) = max(net_ent, 0.0)
        eb(i,j,k-1) = max(-net_ent, 0.0)
      enddo ; enddo
    enddo
    if (CS%debug) then
      call hchksum(ea, "after net flux rearrangement ea", G%HI, scale=GV%H_to_m)
      call hchksum(eb, "after net flux rearrangement eb", G%HI, scale=GV%H_to_m)
    endif
  endif

! Initialize halo regions of ea, eb, and hold to default values.
  !$OMP parallel do default(shared)
  do k=1,nz
    do i=is-1,ie+1
      hold(i,js-1,k) = GV%Angstrom_H ; ea(i,js-1,k) = 0.0 ; eb(i,js-1,k) = 0.0
      hold(i,je+1,k) = GV%Angstrom_H ; ea(i,je+1,k) = 0.0 ; eb(i,je+1,k) = 0.0
    enddo
    do j=js,je
      hold(is-1,j,k) = GV%Angstrom_H ; ea(is-1,j,k) = 0.0 ; eb(is-1,j,k) = 0.0
      hold(ie+1,j,k) = GV%Angstrom_H ; ea(ie+1,j,k) = 0.0 ; eb(ie+1,j,k) = 0.0
    enddo
  enddo

  call cpu_clock_begin(id_clock_pass)
  if (G%symmetric) then ; dir_flag = To_All+Omit_Corners
  else ; dir_flag = To_West+To_South+Omit_Corners ; endif
  call create_group_pass(CS%pass_hold_eb_ea, hold, G%Domain, dir_flag, halo=1)
  call create_group_pass(CS%pass_hold_eb_ea, eb, G%Domain, dir_flag, halo=1)
  call create_group_pass(CS%pass_hold_eb_ea, ea, G%Domain, dir_flag, halo=1)
  call do_group_pass(CS%pass_hold_eb_ea, G%Domain)
  call cpu_clock_end(id_clock_pass)

  !  Use a tridiagonal solver to determine effect of the diapycnal
  !  advection on velocity field. It is assumed that water leaves
  !  or enters the ocean with the surface velocity.
  if (CS%debug) then
    call MOM_state_chksum("before u/v tridiag ", u, v, h, G, GV, US, haloshift=0)
    call hchksum(ea, "before u/v tridiag ea", G%HI, scale=GV%H_to_m)
    call hchksum(eb, "before u/v tridiag eb", G%HI, scale=GV%H_to_m)
    call hchksum(hold, "before u/v tridiag hold", G%HI, scale=GV%H_to_m)
  endif
  call cpu_clock_begin(id_clock_tridiag)

  !$OMP parallel do default(shared) private(hval,b1,d1,c1,eaval)
  do j=js,je
    do I=Isq,Ieq
      if (associated(ADp%du_dt_dia)) ADp%du_dt_dia(I,j,1) = u(I,j,1)
      hval = (hold(i,j,1) + hold(i+1,j,1)) + (ea(i,j,1) + ea(i+1,j,1)) + h_neglect
      b1(I) = 1.0 / (hval + (eb(i,j,1) + eb(i+1,j,1)))
      d1(I) = hval * b1(I)
      u(I,j,1) = b1(I) * (hval * u(I,j,1))
    enddo
    do k=2,nz ; do I=Isq,Ieq
      if (associated(ADp%du_dt_dia)) ADp%du_dt_dia(I,j,k) = u(I,j,k)
      c1(I,k) = (eb(i,j,k-1)+eb(i+1,j,k-1)) * b1(I)
      eaval = ea(i,j,k) + ea(i+1,j,k)
      hval = hold(i,j,k) + hold(i+1,j,k) + h_neglect
      b1(I) = 1.0 / ((eb(i,j,k) + eb(i+1,j,k)) + (hval + d1(I)*eaval))
      d1(I) = (hval + d1(I)*eaval) * b1(I)
      u(I,j,k) = (hval*u(I,j,k) + eaval*u(I,j,k-1))*b1(I)
    enddo ; enddo
    do k=nz-1,1,-1 ; do I=Isq,Ieq
      u(I,j,k) = u(I,j,k) + c1(I,k+1)*u(I,j,k+1)
      if (associated(ADp%du_dt_dia)) &
        ADp%du_dt_dia(I,j,k) = (u(I,j,k) - ADp%du_dt_dia(I,j,k)) * Idt
    enddo ; enddo
    if (associated(ADp%du_dt_dia)) then
      do I=Isq,Ieq
        ADp%du_dt_dia(I,j,nz) = (u(I,j,nz)-ADp%du_dt_dia(I,j,nz)) * Idt
      enddo
    endif
  enddo
  if (CS%debug) then
    call MOM_state_chksum("aft 1st loop tridiag ", u, v, h, G, GV, US, haloshift=0)
  endif
  !$OMP parallel do default(shared) private(hval,b1,d1,c1,eaval)
  do J=Jsq,Jeq
    do i=is,ie
      if (associated(ADp%dv_dt_dia)) ADp%dv_dt_dia(i,J,1) = v(i,J,1)
      hval = (hold(i,j,1) + hold(i,j+1,1)) + (ea(i,j,1) + ea(i,j+1,1)) + h_neglect
      b1(i) = 1.0 / (hval + (eb(i,j,1) + eb(i,j+1,1)))
      d1(I) = hval * b1(I)
      v(i,J,1) = b1(i) * (hval * v(i,J,1))
    enddo
    do k=2,nz ; do i=is,ie
      if (associated(ADp%dv_dt_dia)) ADp%dv_dt_dia(i,J,k) = v(i,J,k)
      c1(i,k) = (eb(i,j,k-1)+eb(i,j+1,k-1)) * b1(i)
      eaval = ea(i,j,k) + ea(i,j+1,k)
      hval = hold(i,j,k) + hold(i,j+1,k) + h_neglect
      b1(i) = 1.0 / ((eb(i,j,k) + eb(i,j+1,k)) + (hval + d1(i)*eaval))
      d1(i) = (hval + d1(i)*eaval) * b1(i)
      v(i,J,k) = (hval*v(i,J,k) + eaval*v(i,J,k-1))*b1(i)
    enddo ; enddo
    do k=nz-1,1,-1 ; do i=is,ie
      v(i,J,k) = v(i,J,k) + c1(i,k+1)*v(i,J,k+1)
      if (associated(ADp%dv_dt_dia)) &
        ADp%dv_dt_dia(i,J,k) = (v(i,J,k) - ADp%dv_dt_dia(i,J,k)) * Idt
    enddo ; enddo
    if (associated(ADp%dv_dt_dia)) then
      do i=is,ie
        ADp%dv_dt_dia(i,J,nz) = (v(i,J,nz)-ADp%dv_dt_dia(i,J,nz)) * Idt
      enddo
    endif
  enddo
  call cpu_clock_end(id_clock_tridiag)
  if (CS%debug) then
    call MOM_state_chksum("after u/v tridiag ", u, v, h, G, GV, US, haloshift=0)
  endif

  ! Diagnose the diapycnal diffusivities and other related quantities.
  if (CS%id_Kd_int  > 0) call post_data(CS%id_Kd_int,  Kd_int,  CS%diag)
  if (CS%id_Kd_heat > 0) call post_data(CS%id_Kd_heat, Kd_heat, CS%diag)
  if (CS%id_Kd_salt > 0) call post_data(CS%id_Kd_salt, Kd_salt, CS%diag)

  if (CS%id_ea > 0) call post_data(CS%id_ea, ea, CS%diag)
  if (CS%id_eb > 0) call post_data(CS%id_eb, eb, CS%diag)

  if (CS%id_dudt_dia > 0) call post_data(CS%id_dudt_dia, ADp%du_dt_dia,  CS%diag)
  if (CS%id_dvdt_dia > 0) call post_data(CS%id_dvdt_dia, ADp%dv_dt_dia,  CS%diag)
  if (CS%id_wd       > 0) call post_data(CS%id_wd,       CDp%diapyc_vel, CS%diag)

  if (CS%id_Tdif > 0) call post_data(CS%id_Tdif, Tdif_flx, CS%diag)
  if (CS%id_Tadv > 0) call post_data(CS%id_Tadv, Tadv_flx, CS%diag)
  if (CS%id_Sdif > 0) call post_data(CS%id_Sdif, Sdif_flx, CS%diag)
  if (CS%id_Sadv > 0) call post_data(CS%id_Sadv, Sadv_flx, CS%diag)

  !! Diagnostics for terms multiplied by fractional thicknesses
  if (CS%id_hf_dudt_dia_2d > 0) then
    allocate(hf_dudt_dia_2d(G%IsdB:G%IedB,G%jsd:G%jed))
    hf_dudt_dia_2d(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      hf_dudt_dia_2d(I,j) = hf_dudt_dia_2d(I,j) + ADp%du_dt_dia(I,j,k) * ADp%diag_hfrac_u(I,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_hf_dudt_dia_2d, hf_dudt_dia_2d, CS%diag)
    deallocate(hf_dudt_dia_2d)
  endif

  if (CS%id_hf_dvdt_dia_2d > 0) then
    allocate(hf_dvdt_dia_2d(G%isd:G%ied,G%JsdB:G%JedB))
    hf_dvdt_dia_2d(:,:) = 0.0
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      hf_dvdt_dia_2d(i,J) = hf_dvdt_dia_2d(i,J) + ADp%dv_dt_dia(i,J,k) * ADp%diag_hfrac_v(i,J,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_hf_dvdt_dia_2d, hf_dvdt_dia_2d, CS%diag)
    deallocate(hf_dvdt_dia_2d)
  endif

  call disable_averaging(CS%diag)

  if (showCallTree) call callTree_leave("layered_diabatic()")

end subroutine layered_diabatic

!> Returns pointers or values of members within the diabatic_CS type. For extensibility,
!! each returned argument is an optional argument
subroutine extract_diabatic_member(CS, opacity_CSp, optics_CSp, evap_CFL_limit, minimum_forcing_depth, &
                                   KPP_CSp, energetic_PBL_CSp, diabatic_aux_CSp, diabatic_halo)
  type(diabatic_CS), intent(in   )           :: CS !< module control structure
  ! All output arguments are optional
  type(opacity_CS),  optional, pointer       :: opacity_CSp !< A pointer to be set to the opacity control structure
  type(optics_type), optional, pointer       :: optics_CSp  !< A pointer to be set to the optics control structure
  type(KPP_CS),      optional, pointer       :: KPP_CSp     !< A pointer to be set to the KPP CS
  type(energetic_PBL_CS), optional, pointer  :: energetic_PBL_CSp !< A pointer to be set to the ePBL CS
  real,              optional, intent(  out) :: evap_CFL_limit !<The largest fraction of a layer that can be
                                                            !! evaporated in one time-step [nondim].
  real,              optional, intent(  out) :: minimum_forcing_depth !< The smallest depth over which heat
                                                            !! and freshwater fluxes are applied [H ~> m or kg m-2].
  type(diabatic_aux_CS), optional, pointer   :: diabatic_aux_CSp !< A pointer to be set to the diabatic_aux
                                                            !! control structure
  integer,           optional, intent(  out) :: diabatic_halo !< The halo size where the diabatic algorithms
                                                            !! assume thermodynamics properties are valid.

  ! Pointers to control structures
  if (present(opacity_CSp))       opacity_CSp => CS%opacity_CSp
  if (present(optics_CSp))        optics_CSp  => CS%optics
  if (present(KPP_CSp))           KPP_CSp     => CS%KPP_CSp
  if (present(energetic_PBL_CSp)) energetic_PBL_CSp => CS%energetic_PBL_CSp

  ! Constants within diabatic_CS
  if (present(evap_CFL_limit))        evap_CFL_limit = CS%evap_CFL_limit
  if (present(minimum_forcing_depth)) minimum_forcing_depth = CS%minimum_forcing_depth
  if (present(diabatic_halo)) diabatic_halo = CS%halo_TS_diff

end subroutine extract_diabatic_member

!> Routine called for adiabatic physics
subroutine adiabatic(h, tv, fluxes, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(inout) :: G      !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h      !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv     !< points to thermodynamic fields
  type(forcing),           intent(inout) :: fluxes !< boundary fluxes
  real,                    intent(in)    :: dt     !< time step [T ~> s]
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(diabatic_CS),       pointer       :: CS     !< module control structure

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: zeros  ! An array of zeros.

  zeros(:,:,:) = 0.0

  call call_tracer_column_fns(h, h, zeros, zeros, fluxes, zeros(:,:,1), dt, G, GV, US, tv, &
                              CS%optics, CS%tracer_flow_CSp, CS%debug)

end subroutine adiabatic


!> This routine diagnoses tendencies from application of diabatic diffusion
!! using ALE algorithm. Note that layer thickness is not altered by
!! diabatic diffusion.
subroutine diagnose_diabatic_diff_tendency(tv, h, temp_old, saln_old, dt, G, GV, US, CS)
  type(ocean_grid_type),                      intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type),                    intent(in) :: GV       !< ocean vertical grid structure
  type(thermo_var_ptrs),                      intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: h        !< thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: temp_old !< temperature prior to diabatic physics
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in) :: saln_old !< salinity prior to diabatic physics [ppt]
  real,                                       intent(in) :: dt       !< time step [T ~> s]
  type(unit_scale_type),                      intent(in) :: US       !< A dimensional unit scaling type
  type(diabatic_CS),                          pointer    :: CS       !< module control structure

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: work_3d
  real, dimension(SZI_(G),SZJ_(G))          :: work_2d
  real :: Idt  ! The inverse of the timestep [T-1 ~> s-1]
  real :: ppt2mks = 0.001  ! Conversion factor from g/kg to kg/kg.
  integer :: i, j, k, is, ie, js, je, nz
  logical :: do_saln_tend   ! Calculate salinity-based tendency diagnosics

  is  = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Idt = 0.0 ; if (dt > 0.0) Idt = 1. / dt
  work_3d(:,:,:) = 0.0
  work_2d(:,:)   = 0.0


  ! temperature tendency
  do k=1,nz ; do j=js,je ; do i=is,ie
    work_3d(i,j,k) = (tv%T(i,j,k)-temp_old(i,j,k))*Idt
  enddo ; enddo ; enddo
  if (CS%id_diabatic_diff_temp_tend > 0) then
    call post_data(CS%id_diabatic_diff_temp_tend, work_3d, CS%diag, alt_h=h)
  endif

  ! heat tendency
  if (CS%id_diabatic_diff_heat_tend > 0 .or. CS%id_diabatic_diff_heat_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = h(i,j,k)*GV%H_to_RZ * tv%C_p * work_3d(i,j,k)
    enddo ; enddo ; enddo
    if (CS%id_diabatic_diff_heat_tend > 0) then
      call post_data(CS%id_diabatic_diff_heat_tend, work_3d, CS%diag, alt_h=h)
    endif
    if (CS%id_diabatic_diff_heat_tend_2d > 0) then
      work_2d(:,:) = 0.0
      do k=1,nz ; do j=js,je ; do i=is,ie
        work_2d(i,j) = work_2d(i,j) + work_3d(i,j,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_diabatic_diff_heat_tend_2d, work_2d, CS%diag)
    endif
  endif

  ! salinity tendency
  do_saln_tend = CS%id_diabatic_diff_saln_tend > 0 &
    .or. CS%id_diabatic_diff_salt_tend > 0 &
    .or. CS%id_diabatic_diff_salt_tend_2d > 0

  if (do_saln_tend) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = (tv%S(i,j,k) - saln_old(i,j,k)) * Idt
    enddo ; enddo ; enddo

    if (CS%id_diabatic_diff_saln_tend > 0) &
      call post_data(CS%id_diabatic_diff_saln_tend, work_3d, CS%diag, alt_h=h)

    ! salt tendency
    if (CS%id_diabatic_diff_salt_tend > 0 .or. CS%id_diabatic_diff_salt_tend_2d > 0) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        work_3d(i,j,k) = h(i,j,k)*GV%H_to_RZ * ppt2mks * work_3d(i,j,k)
      enddo ; enddo ; enddo
      if (CS%id_diabatic_diff_salt_tend > 0) then
        call post_data(CS%id_diabatic_diff_salt_tend, work_3d, CS%diag, alt_h=h)
      endif
      if (CS%id_diabatic_diff_salt_tend_2d > 0) then
        work_2d(:,:) = 0.0
        do k=1,nz ; do j=js,je ; do i=is,ie
          work_2d(i,j) = work_2d(i,j) + work_3d(i,j,k)
        enddo ; enddo ; enddo
        call post_data(CS%id_diabatic_diff_salt_tend_2d, work_2d, CS%diag)
      endif
    endif
  endif

end subroutine diagnose_diabatic_diff_tendency


!> This routine diagnoses tendencies from application of boundary fluxes.
!! These impacts are generally 3d, in particular for penetrative shortwave.
!! Other fluxes contribute 3d in cases when the layers vanish or are very thin,
!! in which case we distribute the flux into k > 1 layers.
subroutine diagnose_boundary_forcing_tendency(tv, h, temp_old, saln_old, h_old, &
                                              dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type), intent(in) :: GV       !< ocean vertical grid structure
  type(thermo_var_ptrs),   intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h        !< thickness after boundary flux application [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: temp_old !< temperature prior to boundary flux application [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: saln_old !< salinity prior to boundary flux application [ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old    !< thickness prior to boundary flux application [H ~> m or kg m-2]
  real,                    intent(in) :: dt       !< time step [T ~> s]
  type(unit_scale_type),   intent(in) :: US       !< A dimensional unit scaling type
  type(diabatic_CS),       pointer    :: CS       !< module control structure

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: work_3d
  real, dimension(SZI_(G),SZJ_(G))          :: work_2d
  real :: Idt  ! The inverse of the timestep [T-1 ~> s-1]
  real :: ppt2mks = 0.001  ! Conversion factor from g/kg to kg/kg.
  integer :: i, j, k, is, ie, js, je, nz

  is  = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Idt = 0.0 ; if (dt > 0.0) Idt = 1. / dt
  work_3d(:,:,:) = 0.0
  work_2d(:,:)   = 0.0

  ! Thickness tendency
  if (CS%id_boundary_forcing_h_tendency > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = (h(i,j,k) - h_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_boundary_forcing_h_tendency, work_3d, CS%diag, alt_h=h_old)
  endif

  ! temperature tendency
  if (CS%id_boundary_forcing_temp_tend > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = (tv%T(i,j,k)-temp_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_boundary_forcing_temp_tend, work_3d, CS%diag, alt_h=h_old)
  endif

  ! heat tendency
  if (CS%id_boundary_forcing_heat_tend > 0 .or. CS%id_boundary_forcing_heat_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = GV%H_to_RZ * tv%C_p * Idt * (h(i,j,k) * tv%T(i,j,k) - h_old(i,j,k) * temp_old(i,j,k))
    enddo ; enddo ; enddo
    if (CS%id_boundary_forcing_heat_tend > 0) then
      call post_data(CS%id_boundary_forcing_heat_tend, work_3d, CS%diag, alt_h=h_old)
    endif
    if (CS%id_boundary_forcing_heat_tend_2d > 0) then
      work_2d(:,:) = 0.0
      do k=1,nz ; do j=js,je ; do i=is,ie
        work_2d(i,j) = work_2d(i,j) + work_3d(i,j,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_boundary_forcing_heat_tend_2d, work_2d, CS%diag)
    endif
  endif

  ! salinity tendency
  if (CS%id_boundary_forcing_saln_tend > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = (tv%S(i,j,k)-saln_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_boundary_forcing_saln_tend, work_3d, CS%diag, alt_h=h_old)
  endif

  ! salt tendency
  if (CS%id_boundary_forcing_salt_tend > 0 .or. CS%id_boundary_forcing_salt_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = GV%H_to_RZ * ppt2mks * Idt * (h(i,j,k) * tv%S(i,j,k) - h_old(i,j,k) * saln_old(i,j,k))
    enddo ; enddo ; enddo
    if (CS%id_boundary_forcing_salt_tend > 0) then
      call post_data(CS%id_boundary_forcing_salt_tend, work_3d, CS%diag, alt_h=h_old)
    endif
    if (CS%id_boundary_forcing_salt_tend_2d > 0) then
      work_2d(:,:) = 0.0
      do k=1,nz ; do j=js,je ; do i=is,ie
        work_2d(i,j) = work_2d(i,j) + work_3d(i,j,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_boundary_forcing_salt_tend_2d, work_2d, CS%diag)
    endif
  endif

end subroutine diagnose_boundary_forcing_tendency


!> This routine diagnoses tendencies for temperature and heat from frazil formation.
!! This routine is called twice from within subroutine diabatic; at start and at
!! end of the diabatic processes. The impacts from frazil are generally a function
!! of depth.  Hence, when checking heat budget, be sure to remove HFSIFRAZIL from HFDS in k=1.
subroutine diagnose_frazil_tendency(tv, h, temp_old, dt, G, GV, US, CS)
  type(ocean_grid_type),                     intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV       !< ocean vertical grid structure
  type(thermo_var_ptrs),                     intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h        !< thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: temp_old !< temperature prior to frazil formation [degC]
  real,                                      intent(in) :: dt       !< time step [T ~> s]
  type(unit_scale_type),                     intent(in) :: US       !< A dimensional unit scaling type
  type(diabatic_CS),                         pointer    :: CS       !< module control structure

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: work_3d
  real, dimension(SZI_(G),SZJ_(G))          :: work_2d
  real    :: Idt ! The inverse of the timestep [T-1 ~> s-1]
  integer :: i, j, k, is, ie, js, je, nz

  is  = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Idt = 0.0 ; if (dt > 0.0) Idt = 1. / dt

  ! temperature tendency
  if (CS%id_frazil_temp_tend > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = Idt * (tv%T(i,j,k)-temp_old(i,j,k))
    enddo ; enddo ; enddo
    call post_data(CS%id_frazil_temp_tend, work_3d, CS%diag)
  endif

  ! heat tendency
  if (CS%id_frazil_heat_tend > 0 .or. CS%id_frazil_heat_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = GV%H_to_RZ * tv%C_p * h(i,j,k) * Idt * (tv%T(i,j,k)-temp_old(i,j,k))
    enddo ; enddo ; enddo
    if (CS%id_frazil_heat_tend > 0) call post_data(CS%id_frazil_heat_tend, work_3d, CS%diag)

    ! As a consistency check, we must have
    ! FRAZIL_HEAT_TENDENCY_2d = HFSIFRAZIL
    if (CS%id_frazil_heat_tend_2d > 0) then
      work_2d(:,:) = 0.0
      do k=1,nz ; do j=js,je ; do i=is,ie
        work_2d(i,j) = work_2d(i,j) + work_3d(i,j,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_frazil_heat_tend_2d, work_2d, CS%diag)
    endif
  endif

end subroutine diagnose_frazil_tendency


!> A simplified version of diabatic_driver_init that will allow
!! tracer column functions to be called without allowing any
!! of the diabatic processes to be used.
subroutine adiabatic_driver_init(Time, G, param_file, diag, CS, &
                                tracer_flow_CSp)
  type(time_type),         intent(in)    :: Time             !< current model time
  type(ocean_grid_type),   intent(in)    :: G                !< model grid structure
  type(param_file_type),   intent(in)    :: param_file       !< the file to parse for parameter values
  type(diag_ctrl), target, intent(inout) :: diag             !< regulates diagnostic output
  type(diabatic_CS),       pointer       :: CS               !< module control structure
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp  !< pointer to control structure of the
                                                             !! tracer flow control module

! This "include" declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_diabatic_driver" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "adiabatic_driver_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag
  if (associated(tracer_flow_CSp)) CS%tracer_flow_CSp => tracer_flow_CSp

! Set default, read and log parameters
  call log_version(param_file, mdl, version, &
                   "The following parameters are used for diabatic processes.")

end subroutine adiabatic_driver_init


!> This routine initializes the diabatic driver module.
subroutine diabatic_driver_init(Time, G, GV, US, param_file, useALEalgorithm, diag, &
                                ADp, CDp, CS, tracer_flow_CSp, sponge_CSp, &
                                ALE_sponge_CSp)
  type(time_type), target                :: Time             !< model time
  type(ocean_grid_type),   intent(inout) :: G                !< model grid structure
  type(verticalGrid_type), intent(in)    :: GV               !< model vertical grid structure
  type(unit_scale_type),   intent(in)    :: US               !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file       !< file to parse for parameter values
  logical,                 intent(in)    :: useALEalgorithm  !< logical for whether to use ALE remapping
  type(diag_ctrl), target, intent(inout) :: diag             !< structure to regulate diagnostic output
  type(accel_diag_ptrs),   intent(inout) :: ADp              !< pointers to accelerations in momentum equations,
                                                             !! to enable diagnostics, like energy budgets
  type(cont_diag_ptrs),    intent(inout) :: CDp              !< pointers to terms in continuity equations
  type(diabatic_CS),       pointer       :: CS               !< module control structure
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp  !< pointer to control structure of the
                                                             !! tracer flow control module
  type(sponge_CS),         pointer       :: sponge_CSp       !< pointer to the sponge module control structure
  type(ALE_sponge_CS),     pointer       :: ALE_sponge_CSp   !< pointer to the ALE sponge module control structure

  real    :: Kd  ! A diffusivity used in the default for other tracer diffusivities, in MKS units [m2 s-1]
  integer :: num_mode
  logical :: use_temperature
  character(len=20) :: EN1, EN2, EN3

! This "include" declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_diabatic_driver" ! This module's name.
  character(len=48)  :: thickness_units
  character(len=40)  :: var_name
  character(len=160) :: var_descript
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz, nbands, m
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "diabatic_driver_init called with an "// &
                            "associated control structure.")
    return
  else
    allocate(CS)
  endif

  CS%diag => diag
  CS%Time => Time

  if (associated(tracer_flow_CSp)) CS%tracer_flow_CSp => tracer_flow_CSp
  if (associated(sponge_CSp))      CS%sponge_CSp      => sponge_CSp
  if (associated(ALE_sponge_CSp))  CS%ALE_sponge_CSp  => ALE_sponge_CSp

  CS%useALEalgorithm = useALEalgorithm
  CS%bulkmixedlayer = (GV%nkml > 0)

  ! Set default, read and log parameters
  call log_version(param_file, mdl, version, &
                   "The following parameters are used for diabatic processes.", &
                   log_to_all=.true., debugging=.true.)
  call get_param(param_file, mdl, "USE_LEGACY_DIABATIC_DRIVER", CS%use_legacy_diabatic, &
                 "If true, use a legacy version of the diabatic subroutine. "//&
                 "This is temporary and is needed to avoid change in answers.", &
                 default=.true.)
  call get_param(param_file, mdl, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified via calls to initialize_sponge and possibly "//&
                 "set_up_sponge_field.", default=.false.)
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", use_temperature, &
                 "If true, temperature and salinity are used as state "//&
                 "variables.", default=.true.)
  call get_param(param_file, mdl, "ENERGETICS_SFC_PBL", CS%use_energetic_PBL, &
                 "If true, use an implied energetics planetary boundary "//&
                 "layer scheme to determine the diffusivity and viscosity "//&
                 "in the surface boundary layer.", default=.false.)
  call get_param(param_file, mdl, "EPBL_IS_ADDITIVE", CS%ePBL_is_additive, &
                 "If true, the diffusivity from ePBL is added to all "//&
                 "other diffusivities. Otherwise, the larger of kappa-shear "//&
                 "and ePBL diffusivities are used.", default=.true.)
  call get_param(param_file, mdl, "PRANDTL_EPBL", CS%ePBL_Prandtl, &
                 "The Prandtl number used by ePBL to convert vertical diffusivities into "//&
                 "viscosities.", default=1.0, units="nondim", do_not_log=.not.CS%use_energetic_PBL)
  call get_param(param_file, mdl, "USE_KPP", CS%use_KPP, &
                 "If true, turns on the [CVMix] KPP scheme of Large et al., 1994, "//&
                 "to calculate diffusivities and non-local transport in the OBL.", &
                 default=.false., do_not_log=.true.)
  CS%use_CVMix_ddiff = CVMix_ddiff_is_used(param_file)

  CS%use_kappa_shear = kappa_shear_is_used(param_file)
  CS%use_CVMix_shear = CVMix_shear_is_used(param_file)

  if (CS%bulkmixedlayer) then
    call get_param(param_file, mdl, "ML_MIX_FIRST", CS%ML_mix_first, &
                 "The fraction of the mixed layer mixing that is applied "//&
                 "before interior diapycnal mixing.  0 by default.", &
                 units="nondim", default=0.0)
    call get_param(param_file, mdl, "NKBL", CS%nkbl, default=2, do_not_log=.true.)
  else
    CS%ML_mix_first = 0.0
  endif
  if (use_temperature) then
    call get_param(param_file, mdl, "DO_GEOTHERMAL", CS%use_geothermal, &
                 "If true, apply geothermal heating.", default=.false.)
  else
    CS%use_geothermal = .false.
  endif
  call get_param(param_file, mdl, "INTERNAL_TIDES", CS%use_int_tides, &
                 "If true, use the code that advances a separate set of "//&
                 "equations for the internal tide energy density.", default=.false.)
  CS%nMode = 1
  if (CS%use_int_tides) then
    call get_param(param_file, mdl, "INTERNAL_TIDE_MODES", CS%nMode, &
                 "The number of distinct internal tide modes "//&
                 "that will be calculated.", default=1, do_not_log=.true.)
    call get_param(param_file, mdl, "UNIFORM_TEST_CG", CS%uniform_test_cg, &
                 "If positive, a uniform group velocity of internal tide for test case", &
                 default=-1., units="m s-1", scale=US%m_s_to_L_T)
  endif

  call get_param(param_file, mdl, "MASSLESS_MATCH_TARGETS", &
                                CS%massless_match_targets, &
                 "If true, the temperature and salinity of massless layers "//&
                 "are kept consistent with their target densities. "//&
                 "Otherwise the properties of massless layers evolve "//&
                 "diffusively to match massive neighboring layers.", &
                 default=.true.)

  call get_param(param_file, mdl, "AGGREGATE_FW_FORCING", CS%aggregate_FW_forcing, &
                 "If true, the net incoming and outgoing fresh water fluxes are combined "//&
                 "and applied as either incoming or outgoing depending on the sign of the net. "//&
                 "If false, the net incoming fresh water flux is added to the model and "//&
                 "thereafter the net outgoing is removed from the topmost non-vanished "//&
                 "layers of the updated state.", default=.true.)

  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_CONSERVATION", CS%debugConservation, &
                 "If true, monitor conservation and extrema.", &
                 default=.false., debuggingParam=.true.)

  call get_param(param_file, mdl, "DEBUG_ENERGY_REQ", CS%debug_energy_req, &
                 "If true, debug the energy requirements.", default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "MIX_BOUNDARY_TRACERS", CS%mix_boundary_tracers, &
                 "If true, mix the passive tracers in massless layers at "//&
                 "the bottom into the interior as though a diffusivity of "//&
                 "KD_MIN_TR were operating.", default=.true.)
  call get_param(param_file, mdl, "MIX_BOUNDARY_TRACER_ALE", CS%mix_boundary_tracer_ALE, &
                 "If true and in ALE mode, mix the passive tracers in massless layers at "//&
                 "the bottom into the interior as though a diffusivity of "//&
                 "KD_MIN_TR were operating.", default=.false., do_not_log=.not.CS%useALEalgorithm)

  if (CS%mix_boundary_tracers .or. CS%mix_boundary_tracer_ALE) then
    call get_param(param_file, mdl, "KD", Kd, default=0.0)
    call get_param(param_file, mdl, "KD_MIN_TR", CS%Kd_min_tr, &
                 "A minimal diffusivity that should always be applied to "//&
                 "tracers, especially in massless layers near the bottom. "//&
                 "The default is 0.1*KD.", units="m2 s-1", default=0.1*Kd, scale=US%m2_s_to_Z2_T)
    call get_param(param_file, mdl, "KD_BBL_TR", CS%Kd_BBL_tr, &
                 "A bottom boundary layer tracer diffusivity that will "//&
                 "allow for explicitly specified bottom fluxes. The "//&
                 "entrainment at the bottom is at least sqrt(Kd_BBL_tr*dt) "//&
                 "over the same distance.", units="m2 s-1", default=0., scale=US%m2_s_to_Z2_T)
  endif

  call get_param(param_file, mdl, "TRACER_TRIDIAG", CS%tracer_tridiag, &
                 "If true, use the passive tracer tridiagonal solver for T and S", &
                 default=.false.)

  call get_param(param_file, mdl, "MINIMUM_FORCING_DEPTH", CS%minimum_forcing_depth, &
                 "The smallest depth over which forcing can be applied. This "//&
                 "only takes effect when near-surface layers become thin "//&
                 "relative to this scale, in which case the forcing tendencies "//&
                 "scaled down by distributing the forcing over this depth scale.", &
                 units="m", default=0.001, scale=GV%m_to_H)
  call get_param(param_file, mdl, "EVAP_CFL_LIMIT", CS%evap_CFL_limit, &
                 "The largest fraction of a layer than can be lost to forcing "//&
                 "(e.g. evaporation, sea-ice formation) in one time-step. The unused "//&
                 "mass loss is passed down through the column.", &
                 units="nondim", default=0.8)

  if (CS%use_energetic_PBL .and. .not.CS%useALEalgorithm) &
    call MOM_error(FATAL, "diabatic_driver_init: "//&
                   "ENERGETICS_SFC_PBL = True is only coded to work when USE_REGRIDDING = True.")

  ! Register all available diagnostics for this module.
  thickness_units = get_thickness_units(GV)

  CS%id_ea_t = register_diag_field('ocean_model', 'ea_t', diag%axesTL, Time, &
      'Layer (heat) entrainment from above per timestep', 'm', conversion=GV%H_to_m)
  CS%id_eb_t = register_diag_field('ocean_model', 'eb_t', diag%axesTL, Time, &
      'Layer (heat) entrainment from below per timestep', 'm', conversion=GV%H_to_m)
  CS%id_ea_s = register_diag_field('ocean_model', 'ea_s', diag%axesTL, Time, &
      'Layer (salt) entrainment from above per timestep', 'm', conversion=GV%H_to_m)
  CS%id_eb_s = register_diag_field('ocean_model', 'eb_s', diag%axesTL, Time, &
      'Layer (salt) entrainment from below per timestep', 'm', conversion=GV%H_to_m)
  ! used by layer diabatic
  CS%id_ea = register_diag_field('ocean_model', 'ea', diag%axesTL, Time, &
      'Layer entrainment from above per timestep', 'm', conversion=GV%H_to_m)
  CS%id_eb = register_diag_field('ocean_model', 'eb', diag%axesTL, Time, &
      'Layer entrainment from below per timestep', 'm', conversion=GV%H_to_m)
  if (.not.CS%useALEalgorithm) then
    CS%id_wd = register_diag_field('ocean_model', 'wd', diag%axesTi, Time, &
      'Diapycnal velocity', 'm s-1', conversion=GV%H_to_m)
    if (CS%id_wd > 0) call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)

    CS%id_dudt_dia = register_diag_field('ocean_model', 'dudt_dia', diag%axesCuL, Time, &
        'Zonal Acceleration from Diapycnal Mixing', 'm s-2', conversion=US%L_T2_to_m_s2)
    CS%id_dvdt_dia = register_diag_field('ocean_model', 'dvdt_dia', diag%axesCvL, Time, &
        'Meridional Acceleration from Diapycnal Mixing', 'm s-2', conversion=US%L_T2_to_m_s2)

    CS%id_hf_dudt_dia_2d = register_diag_field('ocean_model', 'hf_dudt_dia_2d', diag%axesCu1, Time, &
        'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Diapycnal Mixing', &
        'm s-2', conversion=US%L_T2_to_m_s2)
    if (CS%id_hf_dudt_dia_2d > 0) call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)

    CS%id_hf_dvdt_dia_2d = register_diag_field('ocean_model', 'hf_dvdt_dia_2d', diag%axesCv1, Time, &
        'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Diapycnal Mixing', &
        'm s-2', conversion=US%L_T2_to_m_s2)
    if (CS%id_hf_dvdt_dia_2d > 0)  call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,Jsd,JedB,nz)

    if ((CS%id_dudt_dia > 0) .or. (CS%id_hf_dudt_dia_2d > 0)) &
      call safe_alloc_ptr(ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
    if ((CS%id_dvdt_dia > 0) .or. (CS%id_hf_dvdt_dia_2d > 0)) &
      call safe_alloc_ptr(ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
  endif

  if (CS%use_int_tides) then
    CS%id_cg1 = register_diag_field('ocean_model', 'cn1', diag%axesT1, &
                 Time, 'First baroclinic mode (eigen) speed', 'm s-1', conversion=US%L_T_to_m_s)
    allocate(CS%id_cn(CS%nMode)) ; CS%id_cn(:) = -1
    do m=1,CS%nMode
      write(var_name, '("cn_mode",i1)') m
      write(var_descript, '("Baroclinic (eigen) speed of mode ",i1)') m
      CS%id_cn(m) = register_diag_field('ocean_model',var_name, diag%axesT1, &
                   Time, var_descript, 'm s-1', conversion=US%L_T_to_m_s)
      call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    enddo
  endif

  if (use_temperature) then
    CS%id_Tdif = register_diag_field('ocean_model',"Tflx_dia_diff", diag%axesTi, &
        Time, "Diffusive diapycnal temperature flux across interfaces", &
        "degC m s-1", conversion=GV%H_to_m*US%s_to_T)
    if (.not.CS%useALEalgorithm) then
      CS%id_Tadv = register_diag_field('ocean_model',"Tflx_dia_adv", diag%axesTi, &
          Time, "Advective diapycnal temperature flux across interfaces", &
          "degC m s-1", conversion=GV%H_to_m*US%s_to_T)
    endif
    CS%id_Sdif = register_diag_field('ocean_model',"Sflx_dia_diff", diag%axesTi, &
        Time, "Diffusive diapycnal salnity flux across interfaces", &
        "psu m s-1", conversion=GV%H_to_m*US%s_to_T)
    if (.not.CS%useALEalgorithm) then
      CS%id_Sadv = register_diag_field('ocean_model',"Sflx_dia_adv", diag%axesTi, &
          Time, "Advective diapycnal salnity flux across interfaces", &
          "psu m s-1", conversion=GV%H_to_m*US%s_to_T)
    endif
    CS%id_MLD_003 = register_diag_field('ocean_model', 'MLD_003', diag%axesT1, Time, &
        'Mixed layer depth (delta rho = 0.03)', 'm', conversion=US%Z_to_m, &
        cmor_field_name='mlotst', cmor_long_name='Ocean Mixed Layer Thickness Defined by Sigma T', &
        cmor_standard_name='ocean_mixed_layer_thickness_defined_by_sigma_t')
    CS%id_mlotstsq = register_diag_field('ocean_model', 'mlotstsq', diag%axesT1, Time, &
        long_name='Square of Ocean Mixed Layer Thickness Defined by Sigma T', &
        standard_name='square_of_ocean_mixed_layer_thickness_defined_by_sigma_t', &
        units='m2', conversion=US%Z_to_m**2)
    CS%id_MLD_0125 = register_diag_field('ocean_model', 'MLD_0125', diag%axesT1, Time, &
        'Mixed layer depth (delta rho = 0.125)', 'm', conversion=US%Z_to_m)
    call get_param(param_file, mdl, "MLD_EN_VALS", CS%MLD_EN_VALS, &
         "The energy values used to compute MLDs.  If not set (or all set to 0.), the "//&
         "default will overwrite to 25., 2500., 250000.",units='J/m2', default=0., &
         scale=US%kg_m3_to_R*US%m_to_Z**3*US%T_to_s**2)
    if ((CS%MLD_EN_VALS(1)==0.).and.(CS%MLD_EN_VALS(2)==0.).and.(CS%MLD_EN_VALS(3)==0.)) then
      CS%MLD_EN_VALS = (/25.*US%kg_m3_to_R*US%m_to_Z*US%m_to_L**2*US%T_to_s**2,&
           2500.*US%kg_m3_to_R*US%m_to_Z*US%m_to_L**2*US%T_to_s**2,&
           250000.*US%kg_m3_to_R*US%m_to_Z*US%m_to_L**2*US%T_to_s**2/)
    endif
    write(EN1,'(F10.2)') CS%MLD_EN_VALS(1)*US%R_to_kg_m3*US%Z_to_m*US%L_to_m**2*US%s_to_T**2
    write(EN2,'(F10.2)') CS%MLD_EN_VALS(2)*US%R_to_kg_m3*US%Z_to_m*US%L_to_m**2*US%s_to_T**2
    write(EN3,'(F10.2)') CS%MLD_EN_VALS(3)*US%R_to_kg_m3*US%Z_to_m*US%L_to_m**2*US%s_to_T**2
    CS%id_MLD_EN1 = register_diag_field('ocean_model', 'MLD_EN1', diag%axesT1, Time, &
         'Mixed layer depth for energy value set to '//trim(EN1)//' J/m2 (Energy set by 1st MLD_EN_VALS)', &
         'm', conversion=US%Z_to_m)
    CS%id_MLD_EN2 = register_diag_field('ocean_model', 'MLD_EN2', diag%axesT1, Time, &
         'Mixed layer depth for energy value set to '//trim(EN2)//' J/m2 (Energy set by 2nd MLD_EN_VALS)', &
         'm', conversion=US%Z_to_m)
    CS%id_MLD_EN3 = register_diag_field('ocean_model', 'MLD_EN3', diag%axesT1, Time, &
         'Mixed layer depth for energy value set to '//trim(EN3)//' J/m2 (Energy set by 3rd MLD_EN_VALS)', &
         'm', conversion=US%Z_to_m)
    CS%id_subMLN2  = register_diag_field('ocean_model', 'subML_N2', diag%axesT1, Time, &
        'Squared buoyancy frequency below mixed layer', 's-2', conversion=US%s_to_T**2)
    CS%id_MLD_user = register_diag_field('ocean_model', 'MLD_user', diag%axesT1, Time, &
        'Mixed layer depth (used defined)', 'm', conversion=US%Z_to_m)
  endif
  call get_param(param_file, mdl, "DIAG_MLD_DENSITY_DIFF", CS%MLDdensityDifference, &
                 "The density difference used to determine a diagnostic mixed "//&
                 "layer depth, MLD_user, following the definition of Levitus 1982. "//&
                 "The MLD is the depth at which the density is larger than the "//&
                 "surface density by the specified amount.", &
                 units='kg/m3', default=0.1, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "DIAG_DEPTH_SUBML_N2", CS%dz_subML_N2, &
                 "The distance over which to calculate a diagnostic of the "//&
                 "stratification at the base of the mixed layer.", &
                 units='m', default=50.0, scale=US%m_to_Z)

  ! diagnostics for values prior to diabatic and prior to ALE
  CS%id_u_predia = register_diag_field('ocean_model', 'u_predia', diag%axesCuL, Time, &
      'Zonal velocity before diabatic forcing', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_v_predia = register_diag_field('ocean_model', 'v_predia', diag%axesCvL, Time, &
      'Meridional velocity before diabatic forcing', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_h_predia = register_diag_field('ocean_model', 'h_predia', diag%axesTL, Time, &
      'Layer Thickness before diabatic forcing', &
      trim(thickness_units), conversion=GV%H_to_MKS, v_extensive=.true.)
  CS%id_e_predia = register_diag_field('ocean_model', 'e_predia', diag%axesTi, Time, &
      'Interface Heights before diabatic forcing', 'm')
  if (use_temperature) then
    CS%id_T_predia = register_diag_field('ocean_model', 'temp_predia', diag%axesTL, Time, &
        'Potential Temperature', 'degC')
    CS%id_S_predia = register_diag_field('ocean_model', 'salt_predia', diag%axesTL, Time, &
        'Salinity', 'PSU')
  endif


  !call set_diffusivity_init(Time, G, param_file, diag, CS%set_diff_CSp, CS%int_tide_CSp)
  CS%id_Kd_int = register_diag_field('ocean_model', 'Kd_interface', diag%axesTi, Time, &
      'Total diapycnal diffusivity at interfaces', 'm2 s-1', conversion=US%Z2_T_to_m2_s)
  if (CS%use_energetic_PBL) then
      CS%id_Kd_ePBL = register_diag_field('ocean_model', 'Kd_ePBL', diag%axesTi, Time, &
          'ePBL diapycnal diffusivity at interfaces', 'm2 s-1', conversion=US%Z2_T_to_m2_s)
  endif

  CS%id_Kd_heat = register_diag_field('ocean_model', 'Kd_heat', diag%axesTi, Time, &
      'Total diapycnal diffusivity for heat at interfaces', 'm2 s-1', conversion=US%Z2_T_to_m2_s, &
       cmor_field_name='difvho',                                                   &
       cmor_standard_name='ocean_vertical_heat_diffusivity',                       &
       cmor_long_name='Ocean vertical heat diffusivity')
  CS%id_Kd_salt = register_diag_field('ocean_model', 'Kd_salt', diag%axesTi, Time, &
      'Total diapycnal diffusivity for salt at interfaces', 'm2 s-1',  conversion=US%Z2_T_to_m2_s, &
       cmor_field_name='difvso',                                                   &
       cmor_standard_name='ocean_vertical_salt_diffusivity',                       &
       cmor_long_name='Ocean vertical salt diffusivity')

  ! CS%useKPP is set to True if KPP-scheme is to be used, False otherwise.
  ! KPP_init() allocated CS%KPP_Csp and also sets CS%KPPisPassive
  CS%useKPP = KPP_init(param_file, G, GV, US, diag, Time, CS%KPP_CSp, passive=CS%KPPisPassive)
  if (CS%useKPP) then
    allocate( CS%KPP_NLTheat(isd:ied,jsd:jed,nz+1) )   ; CS%KPP_NLTheat(:,:,:)   = 0.
    allocate( CS%KPP_NLTscalar(isd:ied,jsd:jed,nz+1) ) ; CS%KPP_NLTscalar(:,:,:) = 0.
  endif
  if (CS%useKPP) then
    allocate( CS%KPP_buoy_flux(isd:ied,jsd:jed,nz+1) ) ; CS%KPP_buoy_flux(:,:,:) = 0.
    allocate( CS%KPP_temp_flux(isd:ied,jsd:jed) )      ; CS%KPP_temp_flux(:,:)   = 0.
    allocate( CS%KPP_salt_flux(isd:ied,jsd:jed) )      ; CS%KPP_salt_flux(:,:)   = 0.
  endif


  ! diagnostics for tendencies of temp and saln due to diabatic processes
  ! available only for ALE algorithm.
  ! diagnostics for tendencies of temp and heat due to frazil
  CS%id_diabatic_diff_h = register_diag_field('ocean_model', 'diabatic_diff_h', diag%axesTL, Time, &
      long_name='Cell thickness used during diabatic diffusion', &
      units='m', conversion=GV%H_to_m, v_extensive=.true.)
  if (CS%useALEalgorithm) then
    CS%id_diabatic_diff_temp_tend = register_diag_field('ocean_model', &
        'diabatic_diff_temp_tendency', diag%axesTL, Time,              &
        'Diabatic diffusion temperature tendency', 'degC s-1', conversion=US%s_to_T)
    if (CS%id_diabatic_diff_temp_tend > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    CS%id_diabatic_diff_saln_tend = register_diag_field('ocean_model',&
        'diabatic_diff_saln_tendency', diag%axesTL, Time,             &
        'Diabatic diffusion salinity tendency', 'psu s-1', conversion=US%s_to_T)
    if (CS%id_diabatic_diff_saln_tend > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    CS%id_diabatic_diff_heat_tend = register_diag_field('ocean_model',                             &
        'diabatic_heat_tendency', diag%axesTL, Time,                                               &
        'Diabatic diffusion heat tendency',                                                        &
        'W m-2', conversion=US%QRZ_T_to_W_m2, cmor_field_name='opottempdiff',            &
        cmor_standard_name='tendency_of_sea_water_potential_temperature_expressed_as_heat_content_'// &
                           'due_to_parameterized_dianeutral_mixing',                               &
        cmor_long_name='Tendency of sea water potential temperature expressed as heat content '//  &
                       'due to parameterized dianeutral mixing', &
        v_extensive=.true.)
    if (CS%id_diabatic_diff_heat_tend > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    CS%id_diabatic_diff_salt_tend = register_diag_field('ocean_model',                   &
        'diabatic_salt_tendency', diag%axesTL, Time,                                     &
        'Diabatic diffusion of salt tendency',                                           &
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, cmor_field_name='osaltdiff', &
        cmor_standard_name='tendency_of_sea_water_salinity_expressed_as_salt_content_'// &
                           'due_to_parameterized_dianeutral_mixing',                     &
        cmor_long_name='Tendency of sea water salinity expressed as salt content '//     &
                       'due to parameterized dianeutral mixing', &
        v_extensive=.true.)
    if (CS%id_diabatic_diff_salt_tend > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    ! This diagnostic should equal to roundoff if all is working well.
    CS%id_diabatic_diff_heat_tend_2d = register_diag_field('ocean_model',                        &
        'diabatic_heat_tendency_2d', diag%axesT1, Time,                                          &
        'Depth integrated diabatic diffusion heat tendency',                                     &
        'W m-2', conversion=US%QRZ_T_to_W_m2, cmor_field_name='opottempdiff_2d',      &
        cmor_standard_name='tendency_of_sea_water_potential_temperature_expressed_as_heat_content_'//&
                           'due_to_parameterized_dianeutral_mixing_depth_integrated',            &
        cmor_long_name='Tendency of sea water potential temperature expressed as heat content '//&
                       'due to parameterized dianeutral mixing depth integrated')
    if (CS%id_diabatic_diff_heat_tend_2d > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    ! This diagnostic should equal to roundoff if all is working well.
    CS%id_diabatic_diff_salt_tend_2d = register_diag_field('ocean_model',                &
        'diabatic_salt_tendency_2d', diag%axesT1, Time,                                  &
        'Depth integrated diabatic diffusion salt tendency',                             &
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, cmor_field_name='osaltdiff_2d',              &
        cmor_standard_name='tendency_of_sea_water_salinity_expressed_as_salt_content_'// &
                           'due_to_parameterized_dianeutral_mixing_depth_integrated',    &
        cmor_long_name='Tendency of sea water salinity expressed as salt content '//     &
                       'due to parameterized dianeutral mixing depth integrated')
    if (CS%id_diabatic_diff_salt_tend_2d > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    ! diagnostics for tendencies of thickness temp and saln due to boundary forcing
    ! available only for ALE algorithm.
  ! diagnostics for tendencies of temp and heat due to frazil
    CS%id_boundary_forcing_h = register_diag_field('ocean_model', 'boundary_forcing_h', diag%axesTL, Time, &
        long_name='Cell thickness after applying boundary forcing', &
        units='m', conversion=GV%H_to_m, v_extensive=.true.)
    CS%id_boundary_forcing_h_tendency = register_diag_field('ocean_model',   &
        'boundary_forcing_h_tendency', diag%axesTL, Time,                &
        'Cell thickness tendency due to boundary forcing', 'm s-1', &
        conversion=GV%H_to_m*US%s_to_T, v_extensive=.true.)
    if (CS%id_boundary_forcing_h_tendency > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    CS%id_boundary_forcing_temp_tend = register_diag_field('ocean_model',&
        'boundary_forcing_temp_tendency', diag%axesTL, Time,             &
        'Boundary forcing temperature tendency', 'degC s-1', conversion=US%s_to_T)
    if (CS%id_boundary_forcing_temp_tend > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    CS%id_boundary_forcing_saln_tend = register_diag_field('ocean_model',&
        'boundary_forcing_saln_tendency', diag%axesTL, Time,             &
        'Boundary forcing saln tendency', 'psu s-1', conversion=US%s_to_T)
    if (CS%id_boundary_forcing_saln_tend > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    CS%id_boundary_forcing_heat_tend = register_diag_field('ocean_model',&
        'boundary_forcing_heat_tendency', diag%axesTL, Time,             &
        'Boundary forcing heat tendency', &
        'W m-2', conversion=US%QRZ_T_to_W_m2, v_extensive = .true.)
    if (CS%id_boundary_forcing_heat_tend > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    CS%id_boundary_forcing_salt_tend = register_diag_field('ocean_model',&
        'boundary_forcing_salt_tendency', diag%axesTL, Time,             &
        'Boundary forcing salt tendency', 'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, &
        v_extensive = .true.)
    if (CS%id_boundary_forcing_salt_tend > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    ! This diagnostic should equal to surface heat flux if all is working well.
    CS%id_boundary_forcing_heat_tend_2d = register_diag_field('ocean_model',&
        'boundary_forcing_heat_tendency_2d', diag%axesT1, Time,             &
        'Depth integrated boundary forcing of ocean heat', &
        'W m-2', conversion=US%QRZ_T_to_W_m2)
    if (CS%id_boundary_forcing_heat_tend_2d > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    ! This diagnostic should equal to surface salt flux if all is working well.
    CS%id_boundary_forcing_salt_tend_2d = register_diag_field('ocean_model',&
        'boundary_forcing_salt_tendency_2d', diag%axesT1, Time,             &
        'Depth integrated boundary forcing of ocean salt', &
        'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)
    if (CS%id_boundary_forcing_salt_tend_2d > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif
  endif

  ! diagnostics for tendencies of temp and heat due to frazil
  CS%id_frazil_h = register_diag_field('ocean_model', 'frazil_h', diag%axesTL, Time, &
      long_name='Cell Thickness', standard_name='cell_thickness', &
      units='m', conversion=GV%H_to_m, v_extensive=.true.)

  ! diagnostic for tendency of temp due to frazil
  CS%id_frazil_temp_tend = register_diag_field('ocean_model',&
      'frazil_temp_tendency', diag%axesTL, Time,             &
      'Temperature tendency due to frazil formation', 'degC s-1', conversion=US%s_to_T)
  if (CS%id_frazil_temp_tend > 0) then
    CS%frazil_tendency_diag = .true.
  endif

  ! diagnostic for tendency of heat due to frazil
  CS%id_frazil_heat_tend = register_diag_field('ocean_model',&
      'frazil_heat_tendency', diag%axesTL, Time,             &
      'Heat tendency due to frazil formation', &
      'W m-2', conversion=US%QRZ_T_to_W_m2, v_extensive=.true.)
  if (CS%id_frazil_heat_tend > 0) then
    CS%frazil_tendency_diag = .true.
  endif

  ! if all is working propertly, this diagnostic should equal to hfsifrazil
  CS%id_frazil_heat_tend_2d = register_diag_field('ocean_model',&
      'frazil_heat_tendency_2d', diag%axesT1, Time,             &
      'Depth integrated heat tendency due to frazil formation', &
      'W m-2', conversion=US%QRZ_T_to_W_m2)
  if (CS%id_frazil_heat_tend_2d > 0) then
    CS%frazil_tendency_diag = .true.
  endif

  ! CS%use_CVMix_conv is set to True if CVMix convection will be used, otherwise it is False.
  CS%use_CVMix_conv = CVMix_conv_init(Time, G, GV, US, param_file, diag, CS%CVMix_conv_csp)

  call entrain_diffusive_init(Time, G, GV, US, param_file, diag, CS%entrain_diffusive_CSp, &
                              just_read_params=CS%useALEalgorithm)

  ! initialize the geothermal heating module
  if (CS%use_geothermal) &
    call geothermal_init(Time, G, GV, US, param_file, diag, CS%geothermal_CSp, useALEalgorithm)

  ! initialize module for internal tide induced mixing
  if (CS%use_int_tides) then
    call int_tide_input_init(Time, G, GV, US, param_file, diag, CS%int_tide_input_CSp, &
                             CS%int_tide_input)
    call internal_tides_init(Time, G, GV, US, param_file, diag, CS%int_tide_CSp)
  endif

  ! initialize module for setting diffusivities
  call set_diffusivity_init(Time, G, GV, US, param_file, diag, CS%set_diff_CSp, CS%int_tide_CSp, &
                            halo_TS=CS%halo_TS_diff, double_diffuse=CS%double_diffuse)

  if (CS%useKPP .and. (CS%double_diffuse .and. .not.CS%use_CVMix_ddiff)) &
    call MOM_error(FATAL, 'diabatic_driver_init: DOUBLE_DIFFUSION (old method) does not work '//&
                          'with KPP. Please set DOUBLE_DIFFUSION=False and USE_CVMIX_DDIFF=True.')

  ! set up the clocks for this module
  id_clock_entrain = cpu_clock_id('(Ocean diabatic entrain)', grain=CLOCK_MODULE)
  if (CS%bulkmixedlayer) &
    id_clock_mixedlayer = cpu_clock_id('(Ocean mixed layer)', grain=CLOCK_MODULE)
  id_clock_remap = cpu_clock_id('(Ocean vert remap)', grain=CLOCK_MODULE)
  if (CS%use_geothermal) &
    id_clock_geothermal = cpu_clock_id('(Ocean geothermal)', grain=CLOCK_ROUTINE)
  id_clock_set_diffusivity = cpu_clock_id('(Ocean set_diffusivity)', grain=CLOCK_MODULE)
  id_clock_kpp = cpu_clock_id('(Ocean KPP)', grain=CLOCK_MODULE)
  id_clock_tracers = cpu_clock_id('(Ocean tracer_columns)', grain=CLOCK_MODULE_DRIVER+5)
  if (CS%use_sponge) &
    id_clock_sponge = cpu_clock_id('(Ocean sponges)', grain=CLOCK_MODULE)
  id_clock_tridiag = cpu_clock_id('(Ocean diabatic tridiag)', grain=CLOCK_ROUTINE)
  id_clock_pass = cpu_clock_id('(Ocean diabatic message passing)', grain=CLOCK_ROUTINE)
  id_clock_differential_diff = -1 ; if (CS%double_diffuse .and. .not.CS%use_CVMix_ddiff) &
    id_clock_differential_diff = cpu_clock_id('(Ocean differential diffusion)', grain=CLOCK_ROUTINE)

  ! initialize the auxiliary diabatic driver module
  call diabatic_aux_init(Time, G, GV, US, param_file, diag, CS%diabatic_aux_CSp, &
                         CS%useALEalgorithm, CS%use_energetic_PBL)

  ! initialize the boundary layer modules
  if (CS%bulkmixedlayer) &
    call bulkmixedlayer_init(Time, G, GV, US, param_file, diag, CS%bulkmixedlayer_CSp)
  if (CS%use_energetic_PBL) &
    call energetic_PBL_init(Time, G, GV, US, param_file, diag, CS%energetic_PBL_CSp)

  call regularize_layers_init(Time, G, GV, param_file, diag, CS%regularize_layers_CSp)

  if (CS%debug_energy_req) &
    call diapyc_energy_req_init(Time, G, GV, US, param_file, diag, CS%diapyc_en_rec_CSp)

  ! obtain information about the number of bands for penetrative shortwave
  if (use_temperature) then
    call get_param(param_file, mdl, "PEN_SW_NBANDS", nbands, default=1)
    if (nbands > 0) then
      allocate(CS%optics)
      call opacity_init(Time, G, GV, US, param_file, diag, CS%opacity_CSp, CS%optics)
    endif
  endif

  ! Initialize the diagnostic grid storage
  call diag_grid_storage_init(CS%diag_grids_prev, G, diag)

end subroutine diabatic_driver_init


!> Routine to close the diabatic driver module
subroutine diabatic_driver_end(CS)
  type(diabatic_CS), pointer :: CS    !< module control structure

  if (.not.associated(CS)) return

  call diabatic_aux_end(CS%diabatic_aux_CSp)

  call entrain_diffusive_end(CS%entrain_diffusive_CSp)
  call set_diffusivity_end(CS%set_diff_CSp)

  if (CS%useKPP) then
    deallocate( CS%KPP_buoy_flux )
    deallocate( CS%KPP_temp_flux )
    deallocate( CS%KPP_salt_flux )
    deallocate( CS%KPP_NLTheat )
    deallocate( CS%KPP_NLTscalar )
    call KPP_end(CS%KPP_CSp)
  endif

  if (CS%use_CVMix_conv) call CVMix_conv_end(CS%CVMix_conv_csp)

  if (CS%use_energetic_PBL) &
    call energetic_PBL_end(CS%energetic_PBL_CSp)
  if (CS%debug_energy_req) &
    call diapyc_energy_req_end(CS%diapyc_en_rec_CSp)

  if (associated(CS%optics)) then
    call opacity_end(CS%opacity_CSp, CS%optics)
    deallocate(CS%optics)
  endif

  ! GMM, the following is commented out because arrays in
  ! CS%diag_grids_prev are neither pointers or allocatables
  ! and, therefore, cannot be deallocated.

  !call diag_grid_storage_end(CS%diag_grids_prev)

  deallocate(CS)

end subroutine diabatic_driver_end


!> \namespace mom_diabatic_driver
!!
!!  By Robert Hallberg, Alistair Adcroft, and Stephen Griffies
!!
!!    This program contains the subroutine that, along with the
!!  subroutines that it calls, implements diapycnal mass and momentum
!!  fluxes and a bulk mixed layer.  The diapycnal diffusion can be
!!  used without the bulk mixed layer.
!!
!!  \section section_diabatic Outline of MOM diabatic
!!
!!  * diabatic first determines the (diffusive) diapycnal mass fluxes
!!  based on the convergence of the buoyancy fluxes within each layer.
!!
!!  * The dual-stream entrainment scheme of MacDougall and Dewar (JPO,
!!  1997) is used for combined diapycnal advection and diffusion,
!!  calculated implicitly and potentially with the Richardson number
!!  dependent mixing, as described by Hallberg (MWR, 2000).
!!
!!  * Diapycnal advection is the residual of diapycnal diffusion,
!!  so the fully implicit upwind differencing scheme that is used is
!!  entirely appropriate.
!!
!!  * The downward buoyancy flux in each layer is determined from
!!  an implicit calculation based on the previously
!!  calculated flux of the layer above and an estimated flux in the
!!  layer below.  This flux is subject to the following conditions:
!!  (1) the flux in the top and bottom layers are set by the boundary
!!  conditions, and (2) no layer may be driven below an Angstrom thick-
!!  ness.  If there is a bulk mixed layer, the buffer layer is treated
!!  as a fixed density layer with vanishingly small diffusivity.
!!
!!    diabatic takes 5 arguments:  the two velocities (u and v), the
!!  thicknesses (h), a structure containing the forcing fields, and
!!  the length of time over which to act (dt).  The velocities and
!!  thickness are taken as inputs and modified within the subroutine.
!!  There is no limit on the time step.

end module MOM_diabatic_driver
