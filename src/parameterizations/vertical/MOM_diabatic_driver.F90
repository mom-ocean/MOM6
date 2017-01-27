!> This routine drives the diabatic/dianeutral physics for MOM
module MOM_diabatic_driver

! This file is part of MOM6. See LICENSE.md for the license.


use MOM_bulk_mixed_layer,    only : bulkmixedlayer, bulkmixedlayer_init, bulkmixedlayer_CS
use MOM_debugging,           only : hchksum
use MOM_checksum_packages,   only : MOM_state_chksum, MOM_state_stats
use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_CVMix_shear,         only : CVMix_shear_is_used
use MOM_diabatic_aux,        only : diabatic_aux_init, diabatic_aux_end, diabatic_aux_CS
use MOM_diabatic_aux,        only : make_frazil, adjust_salt, insert_brine, differential_diffuse_T_S, triDiagTS
use MOM_diabatic_aux,        only : find_uv_at_h, diagnoseMLDbyDensityDifference, applyBoundaryFluxesInOut
use MOM_diag_mediator,       only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,       only : diag_ctrl, time_type, diag_update_remap_grids
use MOM_diag_mediator,       only : diag_ctrl, query_averaging_enabled
use MOM_diag_to_Z,           only : diag_to_Z_CS, register_Zint_diag, calc_Zint_diags
use MOM_diapyc_energy_req,   only : diapyc_energy_req_init, diapyc_energy_req_end
use MOM_diapyc_energy_req,   only : diapyc_energy_req_calc, diapyc_energy_req_test, diapyc_energy_req_CS
use MOM_diffConvection,      only : diffConvection_CS, diffConvection_init
use MOM_diffConvection,      only : diffConvection_calculate, diffConvection_end
use MOM_domains,             only : pass_var, To_West, To_South
use MOM_domains,             only : create_group_pass, do_group_pass, group_pass_type
use MOM_energetic_PBL,       only : energetic_PBL, energetic_PBL_init
use MOM_energetic_PBL,       only : energetic_PBL_end, energetic_PBL_CS
use MOM_energetic_PBL,       only : energetic_PBL_get_MLD
use MOM_entrain_diffusive,   only : entrainment_diffusive, entrain_diffusive_init
use MOM_entrain_diffusive,   only : entrain_diffusive_end, entrain_diffusive_CS
use MOM_EOS,                 only : calculate_density, calculate_TFreeze
use MOM_EOS,                 only : calculate_specific_vol_derivs
use MOM_error_handler,       only : MOM_error, FATAL, WARNING, callTree_showQuery,MOM_mesg
use MOM_error_handler,       only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,         only : get_param, log_version, param_file_type, read_param
use MOM_forcing_type,        only : forcing, MOM_forcing_chksum
use MOM_forcing_type,        only : calculateBuoyancyFlux2d, forcing_SinglePointPrint
use MOM_geothermal,          only : geothermal, geothermal_init, geothermal_end, geothermal_CS
use MOM_grid,                only : ocean_grid_type
use MOM_io,                  only : vardesc, var_desc
use MOM_int_tide_input,      only : set_int_tide_input, int_tide_input_init
use MOM_int_tide_input,      only : int_tide_input_end, int_tide_input_CS, int_tide_input_type
use MOM_internal_tides,      only : propagate_int_tide
use MOM_internal_tides,      only : internal_tides_init, internal_tides_end, int_tide_CS
use MOM_kappa_shear,         only : kappa_shear_is_used
use MOM_KPP,                 only : KPP_CS, KPP_init, KPP_calculate, KPP_end
use MOM_KPP,                 only : KPP_NonLocalTransport_temp, KPP_NonLocalTransport_saln
use MOM_opacity,             only : opacity_init, set_opacity, opacity_end, opacity_CS
use MOM_regularize_layers,   only : regularize_layers, regularize_layers_init, regularize_layers_CS
use MOM_set_diffusivity,     only : set_diffusivity, set_BBL_TKE
use MOM_set_diffusivity,     only : set_diffusivity_init, set_diffusivity_end
use MOM_set_diffusivity,     only : set_diffusivity_CS
use MOM_shortwave_abs,       only : absorbRemainingSW, optics_type
use MOM_sponge,              only : apply_sponge, sponge_CS
use MOM_ALE_sponge,          only : apply_ALE_sponge, ALE_sponge_CS
use MOM_time_manager,        only : operator(<=), time_type ! for testing itides (BDM)
use MOM_tracer_flow_control, only : call_tracer_column_fns, tracer_flow_control_CS
use MOM_tracer_diabatic,     only : tracer_vertdiff
use MOM_variables,           only : thermo_var_ptrs, vertvisc_type, accel_diag_ptrs
use MOM_variables,           only : cont_diag_ptrs, MOM_thermovar_chksum, p3d
use MOM_verticalGrid,        only : verticalGrid_type
use MOM_wave_speed,          only : wave_speeds
use time_manager_mod,        only : increment_time ! for testing itides (BDM)


implicit none ; private

#include <MOM_memory.h>

public diabatic
public diabatic_driver_init
public diabatic_driver_end
public adiabatic
public adiabatic_driver_init

!> Control structure for this module
type, public :: diabatic_CS ;
  logical :: bulkmixedlayer          !< If true, a refined bulk mixed layer is used with
                                     !! nkml sublayers (and additional buffer layers).
  logical :: use_energetic_PBL       !< If true, use the implicit energetics planetary
                                     !! boundary layer scheme to determine the diffusivity
                                     !! in the surface boundary layer.
  logical :: use_kappa_shear         !< If true, use the kappa_shear module to find the
                                     !! shear-driven diapycnal diffusivity.
  logical :: use_cvmix_shear         !< If true, use the CVMix module to find the
                                     !! shear-driven diapycnal diffusivity.
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
  integer :: nMode = 1               !< Number of baroclinic modes to consider
  logical :: int_tide_source_test    !< If true, apply an arbitrary generation site
                                     !! for internal tide testing (BDM)
  real    :: int_tide_source_x       !< X Location of generation site
                                     !! for internal tide for testing (BDM)
  real    :: int_tide_source_y       !< Y Location of generation site
                                     !! for internal tide for testing (BDM)
  integer :: tlen_days               !< Time interval from start for adding wave source
                                     !! for testing internal tides (BDM)
  logical :: uniform_cg              !< If true, set cg = cg_test everywhere
                                     !! for testing internal tides (BDM)
  real    :: cg_test                 !< Uniform group velocity of internal tide
                                     !! for testing internal tides (BDM)
  type(time_type) :: time_max_source !< For use in testing internal tides (BDM)
  type(time_type) :: time_end        !< For use in testing internal tides (BDM)
  logical :: useALEalgorithm         !< If true, use the ALE algorithm rather than layered
                                     !! isopycnal/stacked shallow water mode. This logical
                                     !! passed by argument to diabatic_driver_init.
  logical :: aggregate_FW_forcing    !< Determines whether net incoming/outgoing surface
                                     !! FW fluxes are applied separately or combined before
                                     !! being applied.
  real :: ML_mix_first               !< The nondimensional fraction of the mixed layer
                                     !! algorithm that is applied before diffusive mixing.
                                     !! The default is 0, while 0.5 gives Strang splitting
                                     !! and 1 is a sensible value too.  Note that if there
                                     !! are convective instabilities in the initial state,
                                     !! the first call may do much more than the second.
  integer :: NKBL                    !< The number of buffer layers (if bulk_mixed_layer)
  logical :: massless_match_targets  !< If true (the default), keep the T & S
                                     !! consistent with the target values.
  logical :: mix_boundary_tracers    !< If true, mix the passive tracers in massless
                                     !! layers at the bottom into the interior as though
                                     !! a diffusivity of Kd_min_tr (see below) were
                                     !! operating.
  real    :: Kd_BBL_tr               !< A bottom boundary layer tracer diffusivity that
                                     !! will allow for explicitly specified bottom fluxes
                                     !! in m2 s-1.  The entrainment at the bottom is at
                                     !! least sqrt(Kd_BBL_tr*dt) over the same distance.
  real    :: Kd_min_tr               !< A minimal diffusivity that should always be
                                     !! applied to tracers, especially in massless layers
                                     !! near the bottom, in m2 s-1.

  logical :: useKPP                  !< use CVmix/KPP diffusivities and non-local transport
  logical :: salt_reject_below_ML    !< If true, add salt below mixed layer (layer mode only)
  logical :: KPPisPassive            !< If true, KPP is in passive mode, not changing answers.
  logical :: useConvection           !< If true, calculate large diffusivities when column
                                     !! is statically unstable.
  logical :: debug                   !< If true, write verbose checksums for debugging purposes.
  logical :: debugConservation       !< If true, monitor conservation and extrema.
  logical :: tracer_tridiag          !< If true, use tracer_vertdiff instead of tridiagTS for
                                     !< vertical diffusion of T and S
  logical :: debug_energy_req        !  If true, test the mixing energy requirement code.
  type(diag_ctrl), pointer :: diag   !< structure used to regulate timing of diagnostic output
  real :: MLDdensityDifference       !< Density difference used to determine MLD_user
  integer :: nsw                     !< SW_NBANDS

  integer :: id_cg1      = -1                 ! diag handle for mode-1 speed (BDM)
  integer, allocatable, dimension(:) :: id_cn ! diag handle for all mode speeds (BDM)
  integer :: id_dudt_dia = -1, id_dvdt_dia = -1, id_wd           = -1
  integer :: id_ea       = -1, id_eb       = -1, id_Kd_z         = -1
  integer :: id_Kd_heat  = -1, id_Kd_salt  = -1, id_Kd_interface = -1, id_Kd_ePBL  = -1
  integer :: id_Tdif_z   = -1, id_Tadv_z   = -1, id_Sdif_z       = -1, id_Sadv_z   = -1
  integer :: id_Tdif     = -1, id_Tadv     = -1, id_Sdif         = -1, id_Sadv     = -1
  integer :: id_MLD_003  = -1, id_MLD_0125  = -1, id_MLD_user     = -1, id_mlotstsq = -1
  integer :: id_subMLN2  = -1, id_brine_lay = -1

  integer :: id_diabatic_diff_temp_tend     = -1
  integer :: id_diabatic_diff_saln_tend     = -1
  integer :: id_diabatic_diff_heat_tend     = -1
  integer :: id_diabatic_diff_salt_tend     = -1
  integer :: id_diabatic_diff_heat_tend_2d  = -1
  integer :: id_diabatic_diff_salt_tend_2d  = -1
  logical :: diabatic_diff_tendency_diag    = .false.

  integer :: id_boundary_forcing_temp_tend    = -1
  integer :: id_boundary_forcing_saln_tend    = -1
  integer :: id_boundary_forcing_heat_tend    = -1
  integer :: id_boundary_forcing_salt_tend    = -1
  integer :: id_boundary_forcing_heat_tend_2d = -1
  integer :: id_boundary_forcing_salt_tend_2d = -1
  logical :: boundary_forcing_tendency_diag   = .false.

  integer :: id_frazil_temp_tend    = -1
  integer :: id_frazil_heat_tend    = -1
  integer :: id_frazil_heat_tend_2d = -1
  logical :: frazil_tendency_diag   = .false.
  real, allocatable, dimension(:,:,:) :: frazil_heat_diag !< diagnose 3d heat tendency from frazil
  real, allocatable, dimension(:,:,:) :: frazil_temp_diag !< diagnose 3d temp tendency from frazil

  real    :: ppt2mks = 0.001

  type(diabatic_aux_CS),        pointer :: diabatic_aux_CSp      => NULL()
  type(entrain_diffusive_CS),   pointer :: entrain_diffusive_CSp => NULL()
  type(bulkmixedlayer_CS),      pointer :: bulkmixedlayer_CSp    => NULL()
  type(energetic_PBL_CS),       pointer :: energetic_PBL_CSp     => NULL()
  type(regularize_layers_CS),   pointer :: regularize_layers_CSp => NULL()
  type(geothermal_CS),          pointer :: geothermal_CSp        => NULL()
  type(int_tide_CS),            pointer :: int_tide_CSp          => NULL()
  type(int_tide_input_CS),      pointer :: int_tide_input_CSp    => NULL()
  type(int_tide_input_type),    pointer :: int_tide_input        => NULL()
  type(opacity_CS),             pointer :: opacity_CSp           => NULL()
  type(set_diffusivity_CS),     pointer :: set_diff_CSp          => NULL()
  type(sponge_CS),              pointer :: sponge_CSp            => NULL()
  type(ALE_sponge_CS),          pointer :: ALE_sponge_CSp        => NULL()
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp       => NULL()
  type(optics_type),            pointer :: optics                => NULL()
  type(diag_to_Z_CS),           pointer :: diag_to_Z_CSp         => NULL()
  type(KPP_CS),                 pointer :: KPP_CSp               => NULL()
  type(diffConvection_CS),      pointer :: Conv_CSp              => NULL()
  type(diapyc_energy_req_CS),   pointer :: diapyc_en_rec_CSp     => NULL()

  type(group_pass_type) :: pass_hold_eb_ea !< For group halo pass

  ! Data arrays for communicating between components
  real, allocatable, dimension(:,:,:) :: KPP_NLTheat    !< KPP non-local transport for heat (m/s)
  real, allocatable, dimension(:,:,:) :: KPP_NLTscalar  !< KPP non-local transport for scalars (m/s)
  real, allocatable, dimension(:,:,:) :: KPP_buoy_flux  !< KPP forcing buoyancy flux (m^2/s^3)
  real, allocatable, dimension(:,:)   :: KPP_temp_flux  !< KPP effective temperature flux (K m/s)
  real, allocatable, dimension(:,:)   :: KPP_salt_flux  !< KPP effective salt flux (ppt m/s)

end type diabatic_CS

! clock ids
integer :: id_clock_entrain, id_clock_mixedlayer, id_clock_set_diffusivity
integer :: id_clock_tracers, id_clock_tridiag, id_clock_pass, id_clock_sponge
integer :: id_clock_geothermal, id_clock_differential_diff, id_clock_remap
integer :: id_clock_kpp

contains

!>  This subroutine imposes the diapycnal mass fluxes and the
!!  accompanying diapycnal advection of momentum and tracers.
subroutine diabatic(u, v, h, tv, fluxes, visc, ADp, CDp, dt, G, GV, CS)
  type(ocean_grid_type),                     intent(inout) :: G      !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u      !< zonal velocity (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v      !< meridional velocity (m/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h      !< thickness (m for Bouss / kg/m2 for non-Bouss)
  type(thermo_var_ptrs),                     intent(inout) :: tv     !< points to thermodynamic fields; unused have NULL ptrs
  type(forcing),                             intent(inout) :: fluxes !< points to forcing fields; unused fields have NULL ptrs
  type(vertvisc_type),                       intent(inout) :: visc   !< vertical viscosities, BBL properies, and related
  type(accel_diag_ptrs),                     intent(inout) :: ADp    !< points to accelerations in momentum equations,
                                                                     !! to enable the later derived diagn, like energy budgets
  type(cont_diag_ptrs),                      intent(inout) :: CDp    !< points to terms in continuity equations
  real,                                      intent(in)    :: dt     !< time increment (seconds)
  type(diabatic_CS),                         pointer       :: CS     !< module control structure

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    ea,     &    ! amount of fluid entrained from the layer above within
                 ! one time step  (m for Bouss, kg/m^2 for non-Bouss)
    eb,     &    ! amount of fluid entrained from the layer below within
                 ! one time step  (m for Bouss, kg/m^2 for non-Bouss)
    Kd,     &    ! diapycnal diffusivity of layers (m^2/sec)
    h_orig, &    ! initial layer thicknesses (m for Bouss, kg/m^2 for non-Bouss)
    h_prebound, &    ! initial layer thicknesses (m for Bouss, kg/m^2 for non-Bouss)
    hold,   &    ! layer thickness before diapycnal entrainment, and later
                 ! the initial layer thicknesses (if a mixed layer is used),
                 ! (m for Bouss, kg/m^2 for non-Bouss)
    dSV_dT, &    ! The partial derivatives of specific volume with temperature
    dSV_dS, &    ! and salinity in m^3/(kg K) and m^3/(kg ppt).
    cTKE,   &    ! convective TKE requirements for each layer in J/m^2.
    u_h,    &    ! zonal and meridional velocities at thickness points after
    v_h          ! entrainment (m/s)

  real, dimension(SZI_(G),SZJ_(G),CS%nMode) :: &
    cn       ! baroclinic gravity wave speeds (formerly cg1 - BDM)

  real, dimension(SZI_(G),SZJ_(G)) :: &
    Rcv_ml   ! coordinate density of mixed layer, used for applying sponges

  real, dimension(SZI_(G),SZJ_(G),G%ke) :: h_diag                ! diagnostic array for thickness
  real, dimension(SZI_(G),SZJ_(G),G%ke) :: temp_diag             ! diagnostic array for temp
  real, dimension(SZI_(G),SZJ_(G),G%ke) :: saln_diag             ! diagnostic array for salinity
  real, dimension(SZI_(G),SZJ_(G))      :: tendency_2d           ! depth integrated content tendency for diagn
  real, dimension(SZI_(G),SZJ_(G))      :: TKE_itidal_input_test ! override of energy input for testing (BDM)

  real :: net_ent  ! The net of ea-eb at an interface.

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target :: &
             ! These are targets so that the space can be shared with eaml & ebml.
    eatr, &  ! The equivalent of ea and eb for tracers, which differ from ea and
    ebtr     ! eb in that they tend to homogenize tracers in massless layers
             ! near the boundaries (m for Bouss and kg/m^2 for non-Bouss)

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), target :: &
    Kd_int,   & ! diapycnal diffusivity of interfaces (m^2/s)
    Kd_heat,  & ! diapycnal diffusivity of heat (m^2/s)
    Kd_salt,  & ! diapycnal diffusivity of salt and passive tracers (m^2/s)
    Kd_ePBL,  & ! test array of diapycnal diffusivities at interfaces (m^2/s)
    Tdif_flx, & ! diffusive diapycnal heat flux across interfaces (K m/s)
    Tadv_flx, & ! advective diapycnal heat flux across interfaces (K m/s)
    Sdif_flx, & ! diffusive diapycnal salt flux across interfaces (ppt m/s)
    Sadv_flx    ! advective diapycnal salt flux across interfaces (ppt m/s)

  ! The following 5 variables are only used with a bulk mixed layer.
  real, pointer, dimension(:,:,:) :: &
    eaml, &  ! The equivalent of ea and eb due to mixed layer processes,
    ebml     ! (m for Bouss and kg/m^2 for non-Bouss).  These will be
             ! pointers to eatr and ebtr so as to reuse the memory as
             ! the arrays are not needed at the same time.

  integer :: kb(SZI_(G),SZJ_(G)) ! index of the lightest layer denser
                                 ! than the buffer laye (nondimensional)

  real :: p_ref_cv(SZI_(G))      ! Reference pressure for the potential
                                 ! density which defines the coordinate
                                 ! variable, set to P_Ref, in Pa.

  logical :: in_boundary(SZI_(G)) ! True if there are no massive layers below,
                                  ! where massive is defined as sufficiently thick that
                                  ! the no-flux boundary conditions have not restricted
                                  ! the entrainment - usually sqrt(Kd*dt).

  real :: b_denom_1    ! The first term in the denominator of b1
                       ! (m for Bouss, kg/m^2 for non-Bouss)
  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected
                       ! (m for Bouss and kg/m^2 for non-Bouss)
  real :: h_neglect2   ! h_neglect^2  (m^2 for Bouss, kg^2/m^4 for non-Bouss)
  real :: add_ent      ! Entrainment that needs to be added when mixing tracers
                       ! (m for Bouss and kg/m^2 for non-Bouss)
  real :: eaval        ! eaval is 2*ea at velocity grid points (m for Bouss, kg/m^2 for non-Bouss)
  real :: hval         ! hval is 2*h at velocity grid points (m for Bouss, kg/m^2 for non-Bouss)
  real :: h_tr         ! h_tr is h at tracer points with a tiny thickness
                       ! added to ensure positive definiteness (m for Bouss, kg/m^2 for non-Bouss)
  real :: Tr_ea_BBL    ! The diffusive tracer thickness in the BBL that is
                       ! coupled to the bottom within a timestep (m)

  real :: htot(SZIB_(G))             ! The summed thickness from the bottom, in m.
  real :: b1(SZIB_(G)), d1(SZIB_(G)) ! b1, c1, and d1 are variables used by the
  real :: c1(SZIB_(G),SZK_(G))       ! tridiagonal solver.

  real :: Ent_int ! The diffusive entrainment rate at an interface
                  ! (H units = m for Bouss, kg/m^2 for non-Bouss).
  real :: dt_mix  ! amount of time over which to apply mixing (seconds)
  real :: Idt     ! inverse time step (1/s)

  type(p3d) :: z_ptrs(7)  ! pointers to diagnostics to be interpolated to depth
  integer :: num_z_diags  ! number of diagnostics to be interpolated to depth
  integer :: z_ids(7)     ! id numbers of diagnostics to be interpolated to depth
  logical :: showCallTree ! If true, show the call tree
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb, m

  integer :: ig, jg      ! global indices for testing testing itide point source (BDM)
  logical :: avg_enabled ! for testing internal tides (BDM)
  real :: Kd_add_here    ! An added diffusivity in m2/s

  is   = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = G%ke
  Isq  = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nkmb = GV%nk_rho_varies
  h_neglect = GV%H_subroundoff ; h_neglect2 = h_neglect*h_neglect
  Kd_heat(:,:,:) = 0.0 ; Kd_salt(:,:,:) = 0.0


  if (nz == 1) return
  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("diabatic(), MOM_diabatic_driver.F90")

  ! set equivalence between the same bits of memory for these arrays
  eaml => eatr ; ebml => ebtr

  ! inverse time step
  Idt = 1.0 / dt

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_diabatic_driver: "// &
         "Module must be initialized before it is used.")

  if (CS%debug) then
    call MOM_state_chksum("Start of diabatic ", u, v, h, G, GV, haloshift=0)
    call MOM_forcing_chksum("Start of diabatic", fluxes, G, haloshift=0)
  endif
  if (CS%debugConservation) call MOM_state_stats('Start of diabatic', u, v, h, tv%T, tv%S, G)

  if (CS%debug_energy_req) &
    call diapyc_energy_req_test(h, dt, tv, G, GV, CS%diapyc_en_rec_CSp)


  call cpu_clock_begin(id_clock_set_diffusivity)
  call set_BBL_TKE(u, v, h, fluxes, visc, G, GV, CS%set_diff_CSp)
  call cpu_clock_end(id_clock_set_diffusivity)

  ! Frazil formation keeps the temperature above the freezing point.
  ! make_frazil is deliberately called at both the beginning and at
  ! the end of the diabatic processes.
  if (ASSOCIATED(tv%T) .AND. ASSOCIATED(tv%frazil)) then

    if(CS%frazil_tendency_diag) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        temp_diag(i,j,k) = tv%T(i,j,k)
      enddo ; enddo ; enddo
    endif

    if (ASSOCIATED(fluxes%p_surf_full)) then
      call make_frazil(h, tv, G, GV, CS%diabatic_aux_CSp, fluxes%p_surf_full)
    else
      call make_frazil(h, tv, G, GV, CS%diabatic_aux_CSp)
    endif
    if (showCallTree) call callTree_waypoint("done with 1st make_frazil (diabatic)")

    if (CS%frazil_tendency_diag) then
      call diagnose_frazil_tendency(tv, h, temp_diag, dt, G, GV, CS, 1)
    endif

  endif

  if (CS%debugConservation) call MOM_state_stats('1st make_frazil', u, v, h, tv%T, tv%S, G)

  if ((CS%ML_mix_first > 0.0) .or. CS%use_geothermal) then
!$OMP parallel do default(none) shared(is,ie,js,je,nz,h_orig,h,eaml,ebml)
    do k=1,nz ; do j=js,je ; do i=is,ie
      h_orig(i,j,k) = h(i,j,k) ; eaml(i,j,k) = 0.0 ; ebml(i,j,k) = 0.0
    enddo ; enddo ; enddo
  endif

  if (CS%use_geothermal) then
    call cpu_clock_begin(id_clock_geothermal)
    call geothermal(h, tv, dt, eaml, ebml, G, GV, CS%geothermal_CSp)
    call cpu_clock_end(id_clock_geothermal)
    if (showCallTree) call callTree_waypoint("geothermal (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('geothermal', u, v, h, tv%T, tv%S, G)
  endif


  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  ! Set_opacity estimates the optical properties of the water column.
  ! It will need to be modified later to include information about the
  ! biological properties and layer thicknesses.
  if (associated(CS%optics)) &
    call set_opacity(CS%optics, fluxes, G, GV, CS%opacity_CSp)

  if (CS%bulkmixedlayer) then
    if (CS%debug) then
      call MOM_forcing_chksum("Before mixedlayer", fluxes, G, haloshift=0)
    endif

    if (CS%ML_mix_first > 0.0) then
!  This subroutine
!    (1) Cools the mixed layer.
!    (2) Performs convective adjustment by mixed layer entrainment.
!    (3) Heats the mixed layer and causes it to detrain to
!        Monin-Obukhov depth or minimum mixed layer depth.
!    (4) Uses any remaining TKE to drive mixed layer entrainment.
!    (5) Possibly splits buffer layer into two isopycnal layers (when using isopycnal coordinate)
      call find_uv_at_h(u, v, h, u_h, v_h, G, GV)

      call cpu_clock_begin(id_clock_mixedlayer)
      if (CS%ML_mix_first < 1.0) then
        ! Changes: h, tv%T, tv%S, eaml and ebml  (G is also inout???)
        call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt*CS%ML_mix_first, &
                            eaml,ebml, G, GV, CS%bulkmixedlayer_CSp, CS%optics, &
                            CS%aggregate_FW_forcing, dt, last_call=.false.)
        if (CS%salt_reject_below_ML) &
          call insert_brine(h, tv, G, GV, fluxes, nkmb, CS%diabatic_aux_CSp, &
                            dt*CS%ML_mix_first, CS%id_brine_lay)
      else
        ! Changes: h, tv%T, tv%S, eaml and ebml  (G is also inout???)
        call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt, eaml, ebml, &
                        G, GV, CS%bulkmixedlayer_CSp, CS%optics, &
                        CS%aggregate_FW_forcing, dt, last_call=.true.)
      endif

      !  Keep salinity from falling below a small but positive threshold.
      !  This constraint is needed for SIS1 ice model, which can extract
      !  more salt than is present in the ocean. SIS2 does not suffer
      !  from this limitation, in which case we can let salinity=0 and still
      !  have salt conserved with SIS2 ice. So for SIS2, we can run with
      !  BOUND_SALINITY=False in MOM.F90.
      if (ASSOCIATED(tv%S) .and. ASSOCIATED(tv%salt_deficit)) &
        call adjust_salt(h, tv, G, GV, CS%diabatic_aux_CSp)
      call cpu_clock_end(id_clock_mixedlayer)
      if (CS%debug) then
        call MOM_state_chksum("After mixedlayer ", u, v, h, G, GV, haloshift=0)
        call MOM_forcing_chksum("After mixedlayer", fluxes, G, haloshift=0)
      endif
      if (showCallTree) call callTree_waypoint("done with 1st bulkmixedlayer (diabatic)")
      if (CS%debugConservation) call MOM_state_stats('1st bulkmixedlayer', u, v, h, tv%T, tv%S, G)
    endif
  endif

  if (CS%debug) then
    call MOM_state_chksum("before find_uv_at_h", u, v, h, G, GV, haloshift=0)
  endif
  if (CS%use_kappa_shear .or. CS%use_CVMix_shear) then
    if ((CS%ML_mix_first > 0.0) .or. CS%use_geothermal) then
      call find_uv_at_h(u, v, h_orig, u_h, v_h, G, GV, eaml, ebml)
      if (CS%debug) then
        call hchksum(eaml, "after find_uv_at_h eaml",G%HI)
        call hchksum(ebml, "after find_uv_at_h ebml",G%HI)
      endif
    else
      call find_uv_at_h(u, v, h, u_h, v_h, G, GV)
    endif
    if (showCallTree) call callTree_waypoint("done with find_uv_at_h (diabatic)")
  endif

  if (CS%use_int_tides) then
    !   This block provides an interface for the unresolved low-mode internal
    ! tide module (BDM).

    ! PROVIDE ENERGY DISTRIBUTION (calculate time-varying energy source)
    call set_int_tide_input(u, v, h, tv, fluxes, CS%int_tide_input, dt, G, GV, &
                            CS%int_tide_input_CSp)
    ! CALCULATE MODAL VELOCITIES
    cn(:,:,:) = 0.0
    if (CS%uniform_cg) then
       ! SET TO CONSTANT VALUE TO TEST PROPAGATE CODE
       do m=1,CS%nMode ; cn(:,:,m) = CS%cg_test ; enddo
    else
       call wave_speeds(h, tv, G, GV, CS%nMode, cn, full_halos=.true.)
       ! uncomment the lines below for a hard-coded cn that changes linearly with latitude
       !do j=G%jsd,G%jed ; do i=G%isd,G%ied
       !  cn(i,j,:) = ((7.-1.)/14000000.)*G%geoLatBu(i,j) + (1.-((7.-1.)/14000000.)*-7000000.)
       !enddo ; enddo
    endif

    if (CS%int_tide_source_test) then
      ! BUILD 2D ARRAY WITH POINT SOURCE FOR TESTING
      !  This block of code should be moved into set_int_tide_input. -RWH
      TKE_itidal_input_test(:,:) = 0.0
      avg_enabled = query_averaging_enabled(CS%diag,time_end=CS%time_end)
      if (CS%time_end <= CS%time_max_source) then
        do j=G%jsc,G%jec ; do i=G%isc,G%iec
          !INPUT ARBITRARY ENERGY POINT SOURCE
          if ((G%idg_offset + i == CS%int_tide_source_x) .and. &
              (G%jdg_offset + j == CS%int_tide_source_y)) then
            TKE_itidal_input_test(i,j) = 1.0
          endif
        enddo ; enddo
      endif
      ! CALL ROUTINE USING PRESCRIBED KE FOR TESTING
      call propagate_int_tide(h, tv, cn, TKE_itidal_input_test, &
                            CS%int_tide_input%tideamp, CS%int_tide_input%Nb, dt, G, GV, CS%int_tide_CSp)
    else
      ! CALL ROUTINE USING CALCULATED KE INPUT
      call propagate_int_tide(h, tv, cn, CS%int_tide_input%TKE_itidal_input, &
                              CS%int_tide_input%tideamp, CS%int_tide_input%Nb, dt, G, GV, CS%int_tide_CSp)
    endif
    if (showCallTree) call callTree_waypoint("done with propagate_int_tide (diabatic)")
  endif

  call cpu_clock_begin(id_clock_set_diffusivity)
  ! Sets: Kd, Kd_int, visc%Kd_extra_T, visc%Kd_extra_S
  ! Also changes: visc%Kd_turb, visc%TKE_turb (not clear that TKE_turb is used as input ????)
  ! And sets visc%Kv_turb
  call set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, CS%optics, visc, dt, G, GV, CS%set_diff_CSp, Kd, Kd_int)
  call cpu_clock_end(id_clock_set_diffusivity)
  if (showCallTree) call callTree_waypoint("done with set_diffusivity (diabatic)")

  if (CS%debug) then
    call MOM_state_chksum("after set_diffusivity ", u, v, h, G, GV, haloshift=0)
    call MOM_forcing_chksum("after set_diffusivity ", fluxes, G, haloshift=0)
    call MOM_thermovar_chksum("after set_diffusivity ", tv, G)
    call hchksum(Kd, "after set_diffusivity Kd",G%HI,haloshift=0)
    call hchksum(Kd_Int, "after set_diffusivity Kd_Int",G%HI,haloshift=0)
  endif


  if (CS%useKPP) then
    call cpu_clock_begin(id_clock_kpp)
    ! KPP needs the surface buoyancy flux but does not update state variables.
    ! We could make this call higher up to avoid a repeat unpacking of the surface fluxes.
    ! Sets: CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux
    ! NOTE: CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux are returned as rates (i.e. stuff per second)
    ! unlike other instances where the fluxes are integrated in time over a time-step.
    call calculateBuoyancyFlux2d(G, GV, fluxes, CS%optics, h, tv%T, tv%S, tv, &
                                 CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux)
    ! The KPP scheme calculates boundary layer diffusivities and non-local transport.
    ! MOM6 implementation of KPP matches the boundary layer to zero interior diffusivity,
    ! since the matching to nonzero interior diffusivity can be problematic.
    ! Changes: Kd_int. Sets: KPP_NLTheat, KPP_NLTscalar

!$OMP parallel default(none) shared(is,ie,js,je,nz,Kd_salt,Kd_int,visc,CS,Kd_heat)
!$OMP do
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        Kd_salt(i,j,k) = Kd_int(i,j,k)
        Kd_heat(i,j,k) = Kd_int(i,j,k)
      enddo ; enddo ; enddo
    if (associated(visc%Kd_extra_S)) then
!$OMP do
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        Kd_salt(i,j,k) = Kd_salt(i,j,k) + visc%Kd_extra_S(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (associated(visc%Kd_extra_T)) then
!$OMP do
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        Kd_heat(i,j,k) = Kd_heat(i,j,k) + visc%Kd_extra_T(i,j,k)
      enddo ; enddo ; enddo
    endif
!$OMP end parallel

    call KPP_calculate(CS%KPP_CSp, G, GV, h, tv%T, tv%S, u, v, tv%eqn_of_state, &
      fluxes%ustar, CS%KPP_buoy_flux, Kd_heat, Kd_salt, visc%Kv_turb, CS%KPP_NLTheat, CS%KPP_NLTscalar)
!$OMP parallel default(none) shared(is,ie,js,je,nz,Kd_salt,Kd_int,visc,CS,Kd_heat)

    if (.not. CS%KPPisPassive) then
!$OMP do
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        Kd_int(i,j,k) = min( Kd_salt(i,j,k),  Kd_heat(i,j,k) )
      enddo ; enddo ; enddo
      if (associated(visc%Kd_extra_S)) then
!$OMP do
        do k=1,nz+1 ; do j=js,je ; do i=is,ie
          visc%Kd_extra_S(i,j,k) = Kd_salt(i,j,k) - Kd_int(i,j,k)
        enddo ; enddo ; enddo
      endif
      if (associated(visc%Kd_extra_T)) then
!$OMP do
        do k=1,nz+1 ; do j=js,je ; do i=is,ie
          visc%Kd_extra_T(i,j,k) = Kd_heat(i,j,k) - Kd_int(i,j,k)
        enddo ; enddo ; enddo
      endif
    endif ! not passive
!$OMP end parallel
    call cpu_clock_end(id_clock_kpp)
    if (showCallTree) call callTree_waypoint("done with KPP_calculate (diabatic)")
    if (CS%debug) then
      call MOM_state_chksum("after KPP", u, v, h, G, GV, haloshift=0)
      call MOM_forcing_chksum("after KPP", fluxes, G, haloshift=0)
      call MOM_thermovar_chksum("after KPP", tv, G)
      call hchksum(Kd, "after KPP Kd",G%HI,haloshift=0)
      call hchksum(Kd_Int, "after KPP Kd_Int",G%HI,haloshift=0)
    endif

  endif  ! endif for KPP

  ! Check for static instabilities and increase Kd_int where unstable
  if (CS%useConvection) call diffConvection_calculate(CS%Conv_CSp, &
         G, GV, h, tv%T, tv%S, tv%eqn_of_state, Kd_int)

  if (CS%useKPP) then

    call cpu_clock_begin(id_clock_kpp)
    if (CS%debug) then
      call hchksum(CS%KPP_temp_flux*GV%H_to_m, "before KPP_applyNLT netHeat",G%HI,haloshift=0)
      call hchksum(CS%KPP_salt_flux*GV%H_to_m, "before KPP_applyNLT netSalt",G%HI,haloshift=0)
      call hchksum(CS%KPP_NLTheat, "before KPP_applyNLT NLTheat",G%HI,haloshift=0)
      call hchksum(CS%KPP_NLTscalar, "before KPP_applyNLT NLTscalar",G%HI,haloshift=0)
    endif
    ! Apply non-local transport of heat and salt
    ! Changes: tv%T, tv%S
    call KPP_NonLocalTransport_temp(CS%KPP_CSp, G, GV, h, CS%KPP_NLTheat,   CS%KPP_temp_flux, dt, tv%T, tv%C_p)
    call KPP_NonLocalTransport_saln(CS%KPP_CSp, G, GV, h, CS%KPP_NLTscalar, CS%KPP_salt_flux, dt, tv%S)
    call cpu_clock_end(id_clock_kpp)
    if (showCallTree) call callTree_waypoint("done with KPP_applyNonLocalTransport (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('KPP_applyNonLocalTransport', u, v, h, tv%T, tv%S, G)

    if (CS%debug) then
      call MOM_state_chksum("after KPP_applyNLT ", u, v, h, G, GV, haloshift=0)
      call MOM_forcing_chksum("after KPP_applyNLT ", fluxes, G, haloshift=0)
      call MOM_thermovar_chksum("after KPP_applyNLT ", tv, G)
    endif

  endif ! endif for KPP

  ! Differential diffusion done here.
  ! Changes: tv%T, tv%S
  ! If using matching within the KPP scheme, then this step needs to provide
  ! a diffusivity and happen before KPP.  But generally in MOM, we do not match
  ! KPP boundary layer to interior, so this diffusivity can be computed when convenient.
  if (associated(visc%Kd_extra_T) .and. associated(visc%Kd_extra_S) .and. associated(tv%T)) then
    call cpu_clock_begin(id_clock_differential_diff)

    call differential_diffuse_T_S(h, tv, visc, dt, G, GV)
    call cpu_clock_end(id_clock_differential_diff)
    if (showCallTree) call callTree_waypoint("done with differential_diffuse_T_S (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('differential_diffuse_T_S', u, v, h, tv%T, tv%S, G)

    ! increment heat and salt diffusivity.
    ! CS%useKPP==.true. already has extra_T and extra_S included
    if(.not. CS%useKPP) then
      do K=2,nz ; do j=js,je ; do i=is,ie
        Kd_heat(i,j,K) = Kd_heat(i,j,K) + visc%Kd_extra_T(i,j,K)
        Kd_salt(i,j,K) = Kd_salt(i,j,K) + visc%Kd_extra_S(i,j,K)
      enddo ; enddo ; enddo
    endif


  endif


  ! This block sets ea, eb from Kd or Kd_int.
  ! If using ALE algorithm, set ea=eb=Kd_int on interfaces for
  ! use in the tri-diagonal solver.
  ! Otherwise, call entrainment_diffusive() which sets ea and eb
  ! based on KD and target densities (ie. does remapping as well).
  if (CS%useALEalgorithm) then

    do j=js,je ; do i=is,ie
      ea(i,j,1) = 0.
    enddo ; enddo
!$OMP parallel do default(none) shared(is,ie,js,je,nz,h_neglect,h,ea,GV,dt,Kd_int,eb) &
!$OMP                          private(hval)
    do k=2,nz ; do j=js,je ; do i=is,ie
      hval=1.0/(h_neglect + 0.5*(h(i,j,k-1) + h(i,j,k)))
      ea(i,j,k) = (GV%m_to_H**2) * dt * hval * Kd_int(i,j,k)
      eb(i,j,k-1) = ea(i,j,k)
    enddo ; enddo ; enddo
    do j=js,je ; do i=is,ie
      eb(i,j,nz) = 0.
    enddo ; enddo
    if (showCallTree) call callTree_waypoint("done setting ea,eb from Kd_int (diabatic)")

  else ! .not. CS%useALEalgorithm
    ! When not using ALE, calculate layer entrainments/detrainments from
    ! diffusivities and differences between layer and target densities
    call cpu_clock_begin(id_clock_entrain)
    ! Calculate appropriately limited diapycnal mass fluxes to account
    ! for diapycnal diffusion and advection.  Sets: ea, eb. Changes: kb
    call Entrainment_diffusive(u, v, h, tv, fluxes, dt, G, GV, CS%entrain_diffusive_CSp, &
                               ea, eb, kb, Kd_Lay=Kd, Kd_int=Kd_int)
    call cpu_clock_end(id_clock_entrain)
    if (showCallTree) call callTree_waypoint("done with Entrainment_diffusive (diabatic)")

  endif ! endif for (CS%useALEalgorithm)

  if (CS%debug) then
    call MOM_forcing_chksum("after calc_entrain ", fluxes, G, haloshift=0)
    call MOM_thermovar_chksum("after calc_entrain ", tv, G)
    call MOM_state_chksum("after calc_entrain ", u, v, h, G, GV, haloshift=0)
    call hchksum(GV%H_to_m*ea, "after calc_entrain ea",G%HI,haloshift=0)
    call hchksum(GV%H_to_m*eb, "after calc_entrain eb",G%HI,haloshift=0)
  endif

  ! Apply forcing when using the ALE algorithm
  if (CS%useALEalgorithm) then
    call cpu_clock_begin(id_clock_remap)

    ! Changes made to following fields:  h, tv%T and tv%S.

    ! save prior values for diagnostics
    if(CS%boundary_forcing_tendency_diag) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        h_diag(i,j,k)    = h(i,j,k)
        temp_diag(i,j,k) = tv%T(i,j,k)
        saln_diag(i,j,k) = tv%S(i,j,k)
      enddo ; enddo ; enddo
    endif

    do k=1,nz ; do j=js,je ; do i=is,ie
        h_prebound(i,j,k) = h(i,j,k)
    enddo ; enddo ; enddo
    if (CS%use_energetic_PBL) then
      call applyBoundaryFluxesInOut(CS%diabatic_aux_CSp, G, GV, dt, fluxes, CS%optics, &
                          h, tv, CS%aggregate_FW_forcing, cTKE, dSV_dT, dSV_dS)
      call calculateBuoyancyFlux2d(G, GV, fluxes, CS%optics, h, tv%T, tv%S, tv, &
           CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux)

      if (CS%debug) then
        call hchksum(ea, "after applyBoundaryFluxes ea",G%HI,haloshift=0)
        call hchksum(eb, "after applyBoundaryFluxes eb",G%HI,haloshift=0)
        call hchksum(cTKE, "after applyBoundaryFluxes cTKE",G%HI,haloshift=0)
        call hchksum(dSV_dT, "after applyBoundaryFluxes dSV_dT",G%HI,haloshift=0)
        call hchksum(dSV_dS, "after applyBoundaryFluxes dSV_dS",G%HI,haloshift=0)
      endif

      call find_uv_at_h(u, v, h, u_h, v_h, G, GV)
      call energetic_PBL(h, u_h, v_h, tv, fluxes, dt, Kd_ePBL, G, GV, &
                         CS%energetic_PBL_CSp, dSV_dT, dSV_dS, cTKE, CS%KPP_Buoy_Flux)

      ! If visc%MLD exists, copy the ePBL's MLD into it
      if (associated(visc%MLD)) then
        call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, visc%MLD, G)
        call pass_var(visc%MLD, G%domain)
      endif

      ! Augment the diffusivities due to those diagnosed in energetic_PBL.
      do K=2,nz ; do j=js,je ; do i=is,ie

        if (CS%ePBL_is_additive) then
          Kd_add_here = Kd_ePBL(i,j,K)
          visc%Kv_turb(i,j,K) = visc%Kv_turb(i,j,K) + Kd_ePBL(i,j,K)
        else
          Kd_add_here = max(Kd_ePBL(i,j,K) - visc%Kd_turb(i,j,K), 0.0)
          visc%Kv_turb(i,j,K) = max(visc%Kv_turb(i,j,K), Kd_ePBL(i,j,K))
        endif
        Ent_int = Kd_add_here * (GV%m_to_H**2 * dt) / &
                    (0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect)
        eb(i,j,k-1) = eb(i,j,k-1) + Ent_int
        ea(i,j,k) = ea(i,j,k) + Ent_int
        Kd_int(i,j,K)  = Kd_int(i,j,K) + Kd_add_here

        ! for diagnostics
        Kd_heat(i,j,K) = Kd_heat(i,j,K) + Kd_int(i,j,K)
        Kd_salt(i,j,K) = Kd_salt(i,j,K) + Kd_int(i,j,K)

      enddo ; enddo ; enddo

      if (CS%debug) then
        call hchksum(ea, "after ePBL ea",G%HI,haloshift=0)
        call hchksum(eb, "after ePBL eb",G%HI,haloshift=0)
        call hchksum(Kd_ePBL, "after ePBL Kd_ePBL",G%HI,haloshift=0)
      endif

    else
      call applyBoundaryFluxesInOut(CS%diabatic_aux_CSp, G, GV, dt, fluxes, CS%optics, &
                                    h, tv, CS%aggregate_FW_forcing)

    endif   ! endif for CS%use_energetic_PBL

    ! diagnose the tendencies due to boundary forcing
    if(CS%boundary_forcing_tendency_diag) then
      call diagnose_boundary_forcing_tendency(tv, h, temp_diag, saln_diag, h_diag, dt, G, GV, CS)
    endif

    call cpu_clock_end(id_clock_remap)
    if (CS%debug) then
      call MOM_forcing_chksum("after applyBoundaryFluxes ", fluxes, G, haloshift=0)
      call MOM_thermovar_chksum("after applyBoundaryFluxes ", tv, G)
      call MOM_state_chksum("after applyBoundaryFluxes ", u, v, h, G, GV, haloshift=0)
    endif
    if (showCallTree) call callTree_waypoint("done with applyBoundaryFluxes (diabatic)")
    if (CS%debugConservation)  call MOM_state_stats('applyBoundaryFluxes', u, v, h, tv%T, tv%S, G)

  endif   ! endif for (CS%useALEalgorithm)

  ! Update h according to divergence of the difference between
  ! ea and eb. We keep a record of the original h in hold.
  ! In the following, the checks for negative values are to guard
  ! against instances where entrainment drives a layer to
  ! negative thickness.  This situation will never happen if
  ! enough iterations are permitted in Calculate_Entrainment.
  ! Even if too few iterations are allowed, it is still guarded
  ! against.  In other words the checks are probably unnecessary.
!$OMP parallel do default(none) shared(is,ie,js,je,nz,hold,h,eb,ea,GV)
  do j=js,je
    do i=is,ie
      hold(i,j,1) = h(i,j,1)
      h(i,j,1) = h(i,j,1) + (eb(i,j,1) - ea(i,j,2))
      hold(i,j,nz) = h(i,j,nz)
      h(i,j,nz) = h(i,j,nz) + (ea(i,j,nz) - eb(i,j,nz-1))
      if (h(i,j,1) <= 0.0) then
        h(i,j,1) = GV%Angstrom
      endif
      if (h(i,j,nz) <= 0.0) then
        h(i,j,nz) = GV%Angstrom
      endif
    enddo
    do k=2,nz-1 ; do i=is,ie
      hold(i,j,k) = h(i,j,k)
      h(i,j,k) = h(i,j,k) + ((ea(i,j,k) - eb(i,j,k-1)) + &
                    (eb(i,j,k) - ea(i,j,k+1)))
      if (h(i,j,k) <= 0.0) then
        h(i,j,k) = GV%Angstrom
      endif
    enddo ; enddo
  enddo
  if (CS%debug) then
    call MOM_state_chksum("after negative check ", u, v, h, G, GV, haloshift=0)
    call MOM_forcing_chksum("after negative check ", fluxes, G, haloshift=0)
    call MOM_thermovar_chksum("after negative check ", tv, G)
  endif
  if (showCallTree) call callTree_waypoint("done with h=ea-eb (diabatic)")
  if (CS%debugConservation) call MOM_state_stats('h=ea-eb', u, v, h, tv%T, tv%S, G)


  ! Here, T and S are updated according to ea and eb.
  ! If using the bulk mixed layer, T and S are also updated
  ! by surface fluxes (in fluxes%*).
  ! This is a very long block.
  if (CS%bulkmixedlayer) then

    if (ASSOCIATED(tv%T)) then
      call cpu_clock_begin(id_clock_tridiag)
      ! Temperature and salinity (as state variables) are treated
      ! differently from other tracers to insure massless layers that
      ! are lighter than the mixed layer have temperatures and salinities
      ! that correspond to their prescribed densities.
      if (CS%massless_match_targets) then
!$OMP parallel do default (none) shared(is,ie,js,je,nkmb,hold,h_neglect,eb,ea,nz,kb,tv) &
!$OMP                           private(h_tr,b1,d1,c1,b_denom_1)
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

    endif ! endif for ASSOCIATED(T)
    if (CS%debugConservation) call MOM_state_stats('BML tridiag', u, v, h, tv%T, tv%S, G)

    if ((CS%ML_mix_first > 0.0) .or. CS%use_geothermal) then
      ! The mixed layer code has already been called, but there is some needed
      ! bookkeeping.
!$OMP parallel do default(none) shared(is,ie,js,je,nz,hold,h_orig,ea,eaml,eb,ebml)
      do k=1,nz ; do j=js,je ; do i=is,ie
        hold(i,j,k) = h_orig(i,j,k)
        ea(i,j,k) = ea(i,j,k) + eaml(i,j,k)
        eb(i,j,k) = eb(i,j,k) + ebml(i,j,k)
      enddo ; enddo ; enddo
      if (CS%debug) then
        call hchksum(GV%H_to_m*ea, "after ea = ea + eaml",G%HI,haloshift=0)
        call hchksum(GV%H_to_m*eb, "after eb = eb + ebml",G%HI,haloshift=0)
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

      call find_uv_at_h(u, v, hold, u_h, v_h, G, GV, ea, eb)
      if (CS%debug) call MOM_state_chksum("find_uv_at_h1 ", u, v, h, G, GV, haloshift=0)

      dt_mix = min(dt,dt*(1.0 - CS%ML_mix_first))
      call cpu_clock_begin(id_clock_mixedlayer)
      ! Changes: h, tv%T, tv%S, ea and eb  (G is also inout???)
      call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt_mix, ea, eb, &
                      G, GV, CS%bulkmixedlayer_CSp, CS%optics, &
                      CS%aggregate_FW_forcing, dt, last_call=.true.)

      if (CS%salt_reject_below_ML) &
        call insert_brine(h, tv, G, GV, fluxes, nkmb, CS%diabatic_aux_CSp, dt_mix, &
                          CS%id_brine_lay)

      !  Keep salinity from falling below a small but positive threshold.
      !  This constraint is needed for SIS1 ice model, which can extract
      !  more salt than is present in the ocean. SIS2 does not suffer
      !  from this limitation, in which case we can let salinity=0 and still
      !  have salt conserved with SIS2 ice. So for SIS2, we can run with
      !  BOUND_SALINITY=False in MOM.F90.
      if (ASSOCIATED(tv%S) .and. ASSOCIATED(tv%salt_deficit)) &
        call adjust_salt(h, tv, G, GV, CS%diabatic_aux_CSp)

      call cpu_clock_end(id_clock_mixedlayer)
      if (showCallTree) call callTree_waypoint("done with 2nd bulkmixedlayer (diabatic)")
      if (CS%debugConservation) call MOM_state_stats('2nd bulkmixedlayer', u, v, h, tv%T, tv%S, G)
    endif

  else  ! following block for when NOT using BULKMIXEDLAYER


    ! calculate change in temperature & salinity due to dia-coordinate surface diffusion
    if (ASSOCIATED(tv%T)) then

      if (CS%debug) then
        call hchksum(GV%H_to_m*ea, "before triDiagTS ea ",G%HI,haloshift=0)
        call hchksum(GV%H_to_m*eb, "before triDiagTS eb ",G%HI,haloshift=0)
      endif
      call cpu_clock_begin(id_clock_tridiag)

      if(CS%diabatic_diff_tendency_diag) then
        do k=1,nz ; do j=js,je ; do i=is,ie
          temp_diag(i,j,k) = tv%T(i,j,k)
          saln_diag(i,j,k) = tv%S(i,j,k)
        enddo ; enddo ; enddo
      endif

      ! Changes T and S via the tridiagonal solver; no change to h
      if(CS%tracer_tridiag) then
          call tracer_vertdiff(hold, ea, eb, dt, tv%T, G, GV)
          call tracer_vertdiff(hold, ea, eb, dt, tv%S, G, GV)
      else
        call triDiagTS(G, GV, is, ie, js, je, hold, ea, eb, tv%T, tv%S)
      endif

      ! diagnose temperature, salinity, heat, and salt tendencies
      if(CS%diabatic_diff_tendency_diag) then
         call diagnose_diabatic_diff_tendency(tv, hold, temp_diag, saln_diag, dt, G, GV, CS)
      endif

      call cpu_clock_end(id_clock_tridiag)
      if (showCallTree) call callTree_waypoint("done with triDiagTS (diabatic)")

    endif  ! endif corresponding to if (ASSOCIATED(tv%T))
    if (CS%debugConservation) call MOM_state_stats('triDiagTS', u, v, h, tv%T, tv%S, G)


  endif  ! endif for the BULKMIXEDLAYER block


  if (CS%debug) then
    call MOM_state_chksum("after mixed layer ", u, v, h, G, GV, haloshift=0)
    call MOM_thermovar_chksum("after mixed layer ", tv, G)
    call hchksum(ea, "after mixed layer ea", G%HI)
    call hchksum(eb, "after mixed layer eb", G%HI)
  endif

  if (.not. CS%useALEalgorithm) then
    call cpu_clock_begin(id_clock_remap)
    call regularize_layers(h, tv, dt, ea, eb, G, GV, CS%regularize_layers_CSp)
    call cpu_clock_end(id_clock_remap)
    if (showCallTree) call callTree_waypoint("done with regularize_layers (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('regularize_layers', u, v, h, tv%T, tv%S, G)
  endif

  ! Whenever thickness changes let the diag manager know, as the
  ! target grids for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  ! diagnostics
  if ((CS%id_Tdif > 0) .or. (CS%id_Tdif_z > 0) .or. &
      (CS%id_Tadv > 0) .or. (CS%id_Tadv_z > 0)) then
    do j=js,je ; do i=is,ie
      Tdif_flx(i,j,1) = 0.0 ; Tdif_flx(i,j,nz+1) = 0.0
      Tadv_flx(i,j,1) = 0.0 ; Tadv_flx(i,j,nz+1) = 0.0
    enddo ; enddo
!$OMP parallel do default(none) shared(is,ie,js,je,nz,Tdif_flx,Idt,ea,eb,Tadv_flx,tv)
    do K=2,nz ; do j=js,je ; do i=is,ie
      Tdif_flx(i,j,K) = (Idt * 0.5*(ea(i,j,k) + eb(i,j,k-1))) * &
                        (tv%T(i,j,k-1) - tv%T(i,j,k))
      Tadv_flx(i,j,K) = (Idt * (ea(i,j,k) - eb(i,j,k-1))) * &
                    0.5*(tv%T(i,j,k-1) + tv%T(i,j,k))
    enddo ; enddo ; enddo
  endif
  if ((CS%id_Sdif > 0) .or. (CS%id_Sdif_z > 0) .or. &
      (CS%id_Sadv > 0) .or. (CS%id_Sadv_z > 0)) then
    do j=js,je ; do i=is,ie
      Sdif_flx(i,j,1) = 0.0 ; Sdif_flx(i,j,nz+1) = 0.0
      Sadv_flx(i,j,1) = 0.0 ; Sadv_flx(i,j,nz+1) = 0.0
    enddo ; enddo
!$OMP parallel do default(none) shared(is,ie,js,je,nz,Sdif_flx,Idt,ea,eb,Sadv_flx,tv)
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
    Tr_ea_BBL = sqrt(dt*CS%Kd_BBL_tr)
!$OMP parallel do default(none) shared(is,ie,js,je,ebtr,nz,G,GV,h,dt,CS,h_neglect,  &
!$OMP                                  ea,eb,Tr_ea_BBL,eatr,visc,hold,h_neglect2 )  &
!$OMP                          private(htot,in_boundary,add_ent)
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
          ! should be much less than the values that have been set in Kd,
          ! perhaps a molecular diffusivity.
          add_ent = ((dt * CS%Kd_min_tr) * GV%m_to_H**2) * &
                    ((h(i,j,k-1)+h(i,j,k)+h_neglect) / &
                     (h(i,j,k-1)*h(i,j,k)+h_neglect2)) - &
                    0.5*(ea(i,j,k) + eb(i,j,k-1))
          if (htot(i) < Tr_ea_BBL) then
            add_ent = max(0.0, add_ent, &
                          (Tr_ea_BBL - htot(i)) - min(ea(i,j,k),eb(i,j,k-1)))
          elseif (add_ent < 0.0) then
            add_ent = 0.0 ; in_boundary(i) = .false.
          endif

          ebtr(i,j,k-1) = eb(i,j,k-1) + add_ent
          eatr(i,j,k) = ea(i,j,k) + add_ent
        else
          ebtr(i,j,k-1) = eb(i,j,k-1) ; eatr(i,j,k) = ea(i,j,k)
        endif
        if (associated(visc%Kd_extra_S)) then ; if (visc%Kd_extra_S(i,j,k) > 0.0) then
          add_ent = ((dt * visc%Kd_extra_S(i,j,k)) * GV%m_to_H**2) / &
             (0.25 * ((h(i,j,k-1) + h(i,j,k)) + (hold(i,j,k-1) + hold(i,j,k))) + &
              h_neglect)
          ebtr(i,j,k-1) = ebtr(i,j,k-1) + add_ent
          eatr(i,j,k) = eatr(i,j,k) + add_ent
        endif ; endif
      enddo ; enddo
      do i=is,ie ; eatr(i,j,1) = ea(i,j,1) ; enddo

    enddo

    if (CS%useALEalgorithm) then
    ! For passive tracers, the changes in thickness due to boundary fluxes has yet to be applied
    ! so hold should be h_orig
      call call_tracer_column_fns(h_prebound, h, ea, eb, fluxes, dt, G, GV, tv, &
                                CS%optics, CS%tracer_flow_CSp, CS%debug, &
                                evap_CFL_limit = CS%diabatic_aux_CSp%evap_CFL_limit, &
                                minimum_forcing_depth = CS%diabatic_aux_CSp%minimum_forcing_depth)
    else
      call call_tracer_column_fns(hold, h, eatr, ebtr, fluxes, dt, G, GV, tv, &
                                CS%optics, CS%tracer_flow_CSp, CS%debug)
    endif

  elseif (associated(visc%Kd_extra_S)) then  ! extra diffusivity for passive tracers

    do j=js,je ; do i=is,ie
      ebtr(i,j,nz) = eb(i,j,nz) ; eatr(i,j,1) = ea(i,j,1)
    enddo ; enddo
!$OMP parallel do default(none) shared(nz,is,ie,js,je,visc,dt,GV,h,hold,h_neglect,&
!$OMP                                  ebtr,eb,eatr,ea )                          &
!$OMP                          private(add_ent)
    do k=nz,2,-1 ; do j=js,je ; do i=is,ie
      if (visc%Kd_extra_S(i,j,k) > 0.0) then
        add_ent = ((dt * visc%Kd_extra_S(i,j,k)) * GV%m_to_H**2) / &
           (0.25 * ((h(i,j,k-1) + h(i,j,k)) + (hold(i,j,k-1) + hold(i,j,k))) + &
            h_neglect)
      else
        add_ent = 0.0
      endif
      ebtr(i,j,k-1) = eb(i,j,k-1) + add_ent
      eatr(i,j,k) = ea(i,j,k) + add_ent
    enddo ; enddo ; enddo

    if (CS%useALEalgorithm) then
    ! For passive tracers, the changes in thickness due to boundary fluxes has yet to be applied
      call call_tracer_column_fns(h_prebound, h, eatr, ebtr, fluxes, dt, G, GV, tv, &
                                  CS%optics, CS%tracer_flow_CSp, CS%debug,&
                                  evap_CFL_limit = CS%diabatic_aux_CSp%evap_CFL_limit, &
                                  minimum_forcing_depth = CS%diabatic_aux_CSp%minimum_forcing_depth)
    else
      call call_tracer_column_fns(hold, h, eatr, ebtr, fluxes, dt, G, GV, tv, &
                                  CS%optics, CS%tracer_flow_CSp, CS%debug)
    endif

  else
    if (CS%useALEalgorithm) then
    ! For passive tracers, the changes in thickness due to boundary fluxes has yet to be applied
      call call_tracer_column_fns(h_prebound, h, eatr, ebtr, fluxes, dt, G, GV, tv, &
                                  CS%optics, CS%tracer_flow_CSp, CS%debug, &
                                  evap_CFL_limit = CS%diabatic_aux_CSp%evap_CFL_limit, &
                                  minimum_forcing_depth = CS%diabatic_aux_CSp%minimum_forcing_depth)
    else
      call call_tracer_column_fns(hold, h, ea, eb, fluxes, dt, G, GV, tv, &
                                  CS%optics, CS%tracer_flow_CSp, CS%debug)
    endif

  endif  ! (CS%mix_boundary_tracers)



  call cpu_clock_end(id_clock_tracers)


  ! sponges
  if (CS%use_sponge) then
    call cpu_clock_begin(id_clock_sponge)
    if (associated(CS%ALE_sponge_CSp)) then
      ! ALE sponge
      call apply_ALE_sponge(h, dt, G, CS%ALE_sponge_CSp)
    else
      ! Layer mode sponge
      if (CS%bulkmixedlayer .and. ASSOCIATED(tv%eqn_of_state)) then
        do i=is,ie ; p_ref_cv(i) = tv%P_Ref ; enddo
!$OMP parallel do default(none) shared(js,je,p_ref_cv,Rcv_ml,is,ie,tv)
        do j=js,je
           call calculate_density(tv%T(:,j,1), tv%S(:,j,1), p_ref_cv, Rcv_ml(:,j), &
                               is, ie-is+1, tv%eqn_of_state)
        enddo
        call apply_sponge(h, dt, G, GV, ea, eb, CS%sponge_CSp, Rcv_ml)
      else
        call apply_sponge(h, dt, G, GV, ea, eb, CS%sponge_CSp)
      endif
    endif
    call cpu_clock_end(id_clock_sponge)
    if (CS%debug) then
      call MOM_state_chksum("apply_sponge ", u, v, h, G, GV, haloshift=0)
      call MOM_thermovar_chksum("apply_sponge ", tv, G)
    endif
  endif ! CS%use_sponge


!$OMP parallel default(none) shared(is,ie,js,je,nz,CDp,Idt,G,GV,ea,eb,CS,hold) private(net_ent)
!   Save the diapycnal mass fluxes as a diagnostic field.
  if (ASSOCIATED(CDp%diapyc_vel)) then
!$OMP do
    do j=js,je
      do K=2,nz ; do i=is,ie
        CDp%diapyc_vel(i,j,K) = Idt * (GV%H_to_m * (ea(i,j,k) - eb(i,j,k-1)))
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
      call hchksum(ea, "before net flux rearrangement ea",G%HI)
      call hchksum(eb, "before net flux rearrangement eb",G%HI)
    endif
!$OMP do
    do j=js,je
      do K=2,GV%nkml ; do i=is,ie
        net_ent = ea(i,j,k) - eb(i,j,k-1)
        ea(i,j,k) = max(net_ent, 0.0)
        eb(i,j,k-1) = max(-net_ent, 0.0)
      enddo ; enddo
    enddo
    if (CS%debug) then
      call hchksum(ea, "after net flux rearrangement ea",G%HI)
      call hchksum(eb, "after net flux rearrangement eb",G%HI)
    endif
  endif


! Initialize halo regions of ea, eb, and hold to default values.
!$OMP do
  do k=1,nz
    do i=is-1,ie+1
      hold(i,js-1,k) = GV%Angstrom ; ea(i,js-1,k) = 0.0 ; eb(i,js-1,k) = 0.0
      hold(i,je+1,k) = GV%Angstrom ; ea(i,je+1,k) = 0.0 ; eb(i,je+1,k) = 0.0
    enddo
    do j=js,je
      hold(is-1,j,k) = GV%Angstrom ; ea(is-1,j,k) = 0.0 ; eb(is-1,j,k) = 0.0
      hold(ie+1,j,k) = GV%Angstrom ; ea(ie+1,j,k) = 0.0 ; eb(ie+1,j,k) = 0.0
    enddo
  enddo
!$OMP end parallel

  call cpu_clock_begin(id_clock_pass)
  if (G%symmetric) then
    call create_group_pass(CS%pass_hold_eb_ea,hold,G%Domain)
    call create_group_pass(CS%pass_hold_eb_ea,eb,G%Domain)
    call create_group_pass(CS%pass_hold_eb_ea,ea,G%Domain)
  else
    call create_group_pass(CS%pass_hold_eb_ea,hold,G%Domain,To_West+To_South)
    call create_group_pass(CS%pass_hold_eb_ea,eb,G%Domain,To_West+To_South)
    call create_group_pass(CS%pass_hold_eb_ea,ea,G%Domain,To_West+To_South)
  endif
  call do_group_pass(CS%pass_hold_eb_ea,G%Domain)
  call cpu_clock_end(id_clock_pass)

  if (.not. CS%useALEalgorithm) then
    !  Use a tridiagonal solver to determine effect of the diapycnal
    !  advection on velocity field. It is assumed that water leaves
    !  or enters the ocean with the surface velocity.
    if (CS%debug) then
      call MOM_state_chksum("before u/v tridiag ", u, v, h, G, GV, haloshift=0)
      call hchksum(ea, "before u/v tridiag ea",G%HI)
      call hchksum(eb, "before u/v tridiag eb",G%HI)
      call hchksum(hold, "before u/v tridiag hold",G%HI)
    endif
    call cpu_clock_begin(id_clock_tridiag)
!$OMP parallel do default(none) shared(js,je,Isq,Ieq,ADp,u,hold,ea,h_neglect,eb,nz,Idt) &
!$OMP                          private(hval,b1,d1,c1,eaval)
    do j=js,je
      do I=Isq,Ieq
        if (ASSOCIATED(ADp%du_dt_dia)) ADp%du_dt_dia(I,j,1) = u(I,j,1)
        hval = (hold(i,j,1) + hold(i+1,j,1)) + (ea(i,j,1) + ea(i+1,j,1)) + h_neglect
        b1(I) = 1.0 / (hval + (eb(i,j,1) + eb(i+1,j,1)))
        d1(I) = hval * b1(I)
        u(I,j,1) = b1(I) * (hval * u(I,j,1))
      enddo
      do k=2,nz ; do I=Isq,Ieq
        if (ASSOCIATED(ADp%du_dt_dia)) ADp%du_dt_dia(I,j,k) = u(I,j,k)
        c1(I,k) = (eb(i,j,k-1)+eb(i+1,j,k-1)) * b1(I)
        eaval = ea(i,j,k) + ea(i+1,j,k)
        hval = hold(i,j,k) + hold(i+1,j,k) + h_neglect
        b1(I) = 1.0 / ((eb(i,j,k) + eb(i+1,j,k)) + (hval + d1(I)*eaval))
        d1(I) = (hval + d1(I)*eaval) * b1(I)
        u(I,j,k) = (hval*u(I,j,k) + eaval*u(I,j,k-1))*b1(I)
      enddo ; enddo
      do k=nz-1,1,-1 ; do I=Isq,Ieq
        u(I,j,k) = u(I,j,k) + c1(I,k+1)*u(I,j,k+1)
        if (ASSOCIATED(ADp%du_dt_dia)) &
          ADp%du_dt_dia(I,j,k) = (u(I,j,k) - ADp%du_dt_dia(I,j,k)) * Idt
      enddo ; enddo
      if (ASSOCIATED(ADp%du_dt_dia)) then
        do I=Isq,Ieq
          ADp%du_dt_dia(I,j,nz) = (u(I,j,nz)-ADp%du_dt_dia(I,j,nz)) * Idt
        enddo
      endif
    enddo
    if (CS%debug) then
      call MOM_state_chksum("aft 1st loop tridiag ", u, v, h, G, GV, haloshift=0)
    endif
!$OMP parallel do default(none) shared(Jsq,Jeq,is,ie,ADp,v,hold,ea,h_neglect,eb,nz,Idt) &
!$OMP                          private(hval,b1,d1,c1,eaval)
    do J=Jsq,Jeq
      do i=is,ie
        if (ASSOCIATED(ADp%dv_dt_dia)) ADp%dv_dt_dia(i,J,1) = v(i,J,1)
        hval = (hold(i,j,1) + hold(i,j+1,1)) + (ea(i,j,1) + ea(i,j+1,1)) + h_neglect
        b1(i) = 1.0 / (hval + (eb(i,j,1) + eb(i,j+1,1)))
        d1(I) = hval * b1(I)
        v(i,J,1) = b1(i) * (hval * v(i,J,1))
      enddo
      do k=2,nz ; do i=is,ie
        if (ASSOCIATED(ADp%dv_dt_dia)) ADp%dv_dt_dia(i,J,k) = v(i,J,k)
        c1(i,k) = (eb(i,j,k-1)+eb(i,j+1,k-1)) * b1(i)
        eaval = ea(i,j,k) + ea(i,j+1,k)
        hval = hold(i,j,k) + hold(i,j+1,k) + h_neglect
        b1(i) = 1.0 / ((eb(i,j,k) + eb(i,j+1,k)) + (hval + d1(i)*eaval))
        d1(i) = (hval + d1(i)*eaval) * b1(i)
        v(i,J,k) = (hval*v(i,J,k) + eaval*v(i,J,k-1))*b1(i)
      enddo ; enddo
      do k=nz-1,1,-1 ; do i=is,ie
        v(i,J,k) = v(i,J,k) + c1(i,k+1)*v(i,J,k+1)
        if (ASSOCIATED(ADp%dv_dt_dia)) &
          ADp%dv_dt_dia(i,J,k) = (v(i,J,k) - ADp%dv_dt_dia(i,J,k)) * Idt
      enddo ; enddo
      if (ASSOCIATED(ADp%dv_dt_dia)) then
        do i=is,ie
          ADp%dv_dt_dia(i,J,nz) = (v(i,J,nz)-ADp%dv_dt_dia(i,J,nz)) * Idt
        enddo
      endif
    enddo
    call cpu_clock_end(id_clock_tridiag)
    if (CS%debug) then
      call MOM_state_chksum("after u/v tridiag ", u, v, h, G, GV, haloshift=0)
    endif
  endif ! useALEalgorithm

  ! Frazil formation keeps temperature above the freezing point.
  ! make_frazil is deliberately called at both the beginning and at
  ! the end of the diabatic processes.
  if (ASSOCIATED(tv%T) .AND. ASSOCIATED(tv%frazil)) then

    if(CS%frazil_tendency_diag) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        temp_diag(i,j,k) = tv%T(i,j,k)
      enddo ; enddo ; enddo
    endif

    if (ASSOCIATED(fluxes%p_surf_full)) then
      call make_frazil(h, tv, G, GV, CS%diabatic_aux_CSp, fluxes%p_surf_full)
    else
      call make_frazil(h, tv, G, GV, CS%diabatic_aux_CSp)
    endif

    if (CS%frazil_tendency_diag) then
      call diagnose_frazil_tendency(tv, h, temp_diag, dt, G, GV, CS, 2)
    endif

    if (showCallTree) call callTree_waypoint("done with 2nd make_frazil (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('2nd make_frazil', u, v, h, tv%T, tv%S, G)

  endif  ! endif for frazil


  ! Diagnose the diapycnal diffusivities and other related quantities.
  if (CS%id_Kd_interface > 0) call post_data(CS%id_Kd_interface, Kd_int,  CS%diag)
  if (CS%id_Kd_heat      > 0) call post_data(CS%id_Kd_heat,      Kd_heat, CS%diag)
  if (CS%id_Kd_salt      > 0) call post_data(CS%id_Kd_salt,      Kd_salt, CS%diag)
  if (CS%id_Kd_ePBL      > 0) call post_data(CS%id_Kd_ePBL,      Kd_ePBL, CS%diag)

  if (CS%id_ea       > 0) call post_data(CS%id_ea,       ea, CS%diag)
  if (CS%id_eb       > 0) call post_data(CS%id_eb,       eb, CS%diag)

  if (CS%id_dudt_dia > 0) call post_data(CS%id_dudt_dia, ADp%du_dt_dia,  CS%diag)
  if (CS%id_dvdt_dia > 0) call post_data(CS%id_dvdt_dia, ADp%dv_dt_dia,  CS%diag)
  if (CS%id_wd       > 0) call post_data(CS%id_wd,       CDp%diapyc_vel, CS%diag)

  if (CS%id_MLD_003 > 0 .or. CS%id_subMLN2 > 0 .or. CS%id_mlotstsq > 0) then
    call diagnoseMLDbyDensityDifference(CS%id_MLD_003, h, tv, 0.03, G, GV, CS%diag, &
                                        id_N2subML=CS%id_subMLN2, id_MLDsq=CS%id_mlotstsq)
  endif
  if (CS%id_MLD_0125 > 0) then
    call diagnoseMLDbyDensityDifference(CS%id_MLD_0125, h, tv, 0.125, G, GV, CS%diag)
  endif
  if (CS%id_MLD_user > 0) then
    call diagnoseMLDbyDensityDifference(CS%id_MLD_user, h, tv, CS%MLDdensityDifference, G, GV, CS%diag)
  endif

  if (CS%id_Tdif > 0) call post_data(CS%id_Tdif, Tdif_flx, CS%diag)
  if (CS%id_Tadv > 0) call post_data(CS%id_Tadv, Tadv_flx, CS%diag)
  if (CS%id_Sdif > 0) call post_data(CS%id_Sdif, Sdif_flx, CS%diag)
  if (CS%id_Sadv > 0) call post_data(CS%id_Sadv, Sadv_flx, CS%diag)
  if (CS%use_int_tides) then
    if (CS%id_cg1 > 0) call post_data(CS%id_cg1, cn(:,:,1),CS%diag)
    do m=1,CS%nMode
      if (CS%id_cn(m) > 0) call post_data(CS%id_cn(m),cn(:,:,m),CS%diag)
    enddo
  endif

  num_z_diags = 0
  if (CS%id_Kd_z > 0) then
    num_z_diags = num_z_diags + 1
    z_ids(num_z_diags) = CS%id_Kd_z ; z_ptrs(num_z_diags)%p => Kd_int
  endif
  if (CS%id_Tdif_z > 0) then
    num_z_diags = num_z_diags + 1
    z_ids(num_z_diags) = CS%id_Tdif_z ; z_ptrs(num_z_diags)%p => Tdif_flx
  endif
  if (CS%id_Tadv_z > 0) then
    num_z_diags = num_z_diags + 1
    z_ids(num_z_diags) = CS%id_Tadv_z ; z_ptrs(num_z_diags)%p => Tadv_flx
  endif
  if (CS%id_Sdif_z > 0) then
    num_z_diags = num_z_diags + 1
    z_ids(num_z_diags) = CS%id_Sdif_z ; z_ptrs(num_z_diags)%p => Sdif_flx
  endif
  if (CS%id_Sadv_z > 0) then
    num_z_diags = num_z_diags + 1
    z_ids(num_z_diags) = CS%id_Sadv_z ; z_ptrs(num_z_diags)%p => Sadv_flx
  endif

  if (num_z_diags > 0) &
    call calc_Zint_diags(h, z_ptrs, z_ids, num_z_diags, G, GV, CS%diag_to_Z_CSp)

  if (CS%debugConservation) call MOM_state_stats('leaving diabatic', u, v, h, tv%T, tv%S, G)
  if (showCallTree) call callTree_leave("diabatic()")

end subroutine diabatic


!> Routine called for adiabatic physics
subroutine adiabatic(h, tv, fluxes, dt, G, GV, CS)
  type(ocean_grid_type),                    intent(inout) :: G      !< ocean grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h      !< thickness (m for Bouss or kg/m2 for non-Bouss)
  type(thermo_var_ptrs),                    intent(inout) :: tv     !< points to thermodynamic fields
  type(forcing),                            intent(inout) :: fluxes !< boundary fluxes
  real,                                     intent(in)    :: dt     !< time step (seconds)
  type(verticalGrid_type),                  intent(in)    :: GV     !< ocean vertical grid structure
  type(diabatic_CS),                        pointer       :: CS     !< module control structure

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: zeros  ! An array of zeros.

  zeros(:,:,:) = 0.0

  call call_tracer_column_fns(h, h, zeros, zeros, fluxes, dt, G, GV, tv, &
                              CS%optics, CS%tracer_flow_CSp, CS%debug)

end subroutine adiabatic


!> This routine diagnoses tendencies from application of diabatic diffusion
!! using ALE algorithm. Note that layer thickness is not altered by
!! diabatic diffusion.
subroutine diagnose_diabatic_diff_tendency(tv, h, temp_old, saln_old, dt, G, GV, CS)
  type(ocean_grid_type),                     intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV       !< ocean vertical grid structure
  type(thermo_var_ptrs),                     intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: h        !< thickness (m or kg/m2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: temp_old !< temperature prior to diabatic physics
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: saln_old !< salinity prior to diabatic physics (PPT)
  real,                                      intent(in) :: dt       !< time step (sec)
  type(diabatic_CS),                         pointer    :: CS       !< module control structure

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: work_3d
  real, dimension(SZI_(G),SZJ_(G))         :: work_2d
  real    :: Idt
  integer :: i, j, k, is, ie, js, je, nz

  is  = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Idt = 1/dt
  work_3d(:,:,:) = 0.0
  work_2d(:,:)   = 0.0


  ! temperature tendency
  do k=1,nz ; do j=js,je ; do i=is,ie
    work_3d(i,j,k) = (tv%T(i,j,k)-temp_old(i,j,k))*Idt
  enddo ; enddo ; enddo
  if(CS%id_diabatic_diff_temp_tend > 0) then
    call post_data(CS%id_diabatic_diff_temp_tend, work_3d, CS%diag)
  endif

  ! heat tendency
  if(CS%id_diabatic_diff_heat_tend > 0 .or. CS%id_diabatic_diff_heat_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = h(i,j,k) * GV%H_to_kg_m2 * tv%C_p * work_3d(i,j,k)
    enddo ; enddo ; enddo
    if(CS%id_diabatic_diff_heat_tend > 0) then
      call post_data(CS%id_diabatic_diff_heat_tend, work_3d, CS%diag)
    endif
    if(CS%id_diabatic_diff_heat_tend_2d > 0) then
      do j=js,je ; do i=is,ie
        work_2d(i,j) = 0.0
        do k=1,nz
          work_2d(i,j) = work_2d(i,j) + work_3d(i,j,k)
        enddo
      enddo ; enddo
      call post_data(CS%id_diabatic_diff_heat_tend_2d, work_2d, CS%diag)
    endif
  endif

  ! salinity tendency
  do k=1,nz ; do j=js,je ; do i=is,ie
    work_3d(i,j,k) = (tv%S(i,j,k)-saln_old(i,j,k))*Idt
  enddo ; enddo ; enddo
  if(CS%id_diabatic_diff_saln_tend > 0) then
    call post_data(CS%id_diabatic_diff_saln_tend, work_3d, CS%diag)
  endif

  ! salt tendency
  if(CS%id_diabatic_diff_salt_tend > 0 .or. CS%id_diabatic_diff_salt_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = h(i,j,k) * GV%H_to_kg_m2 * CS%ppt2mks * work_3d(i,j,k)
    enddo ; enddo ; enddo
    if(CS%id_diabatic_diff_salt_tend > 0) then
      call post_data(CS%id_diabatic_diff_salt_tend, work_3d, CS%diag)
    endif
    if(CS%id_diabatic_diff_salt_tend_2d > 0) then
      do j=js,je ; do i=is,ie
        work_2d(i,j) = 0.0
        do k=1,nz
          work_2d(i,j) = work_2d(i,j) + work_3d(i,j,k)
        enddo
      enddo ; enddo
      call post_data(CS%id_diabatic_diff_salt_tend_2d, work_2d, CS%diag)
    endif
  endif

end subroutine diagnose_diabatic_diff_tendency


!> This routine diagnoses tendencies from application of boundary fluxes.
!! These impacts are generally 3d, in particular for penetrative shortwave.
!! Other fluxes contribute 3d in cases when the layers vanish or are very thin,
!! in which case we distribute the flux into k > 1 layers.
subroutine diagnose_boundary_forcing_tendency(tv, h, temp_old, saln_old, h_old, &
                                              dt, G, GV, CS)
  type(ocean_grid_type),                    intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type),                  intent(in) :: GV       !< ocean vertical grid structure
  type(thermo_var_ptrs),                    intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h        !< thickness after boundary flux application (m or kg/m2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: temp_old !< temperature prior to boundary flux application
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: saln_old !< salinity prior to boundary flux application (PPT)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h_old    !< thickness prior to boundary flux application (m or kg/m2)
  real,                                     intent(in) :: dt       !< time step (sec)
  type(diabatic_CS),                        pointer    :: CS       !< module control structure

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: work_3d
  real, dimension(SZI_(G),SZJ_(G))         :: work_2d
  real    :: Idt
  integer :: i, j, k, is, ie, js, je, nz

  is  = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Idt = 1/dt
  work_3d(:,:,:) = 0.0
  work_2d(:,:)   = 0.0

  ! temperature tendency
  if(CS%id_boundary_forcing_temp_tend > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = (tv%T(i,j,k)-temp_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_boundary_forcing_temp_tend, work_3d, CS%diag)
  endif

  ! heat tendency
  if(CS%id_boundary_forcing_heat_tend > 0 .or. CS%id_boundary_forcing_heat_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = GV%H_to_kg_m2 * tv%C_p * Idt * (h(i,j,k) * tv%T(i,j,k) - h_old(i,j,k) * temp_old(i,j,k))
    enddo ; enddo ; enddo
    if(CS%id_boundary_forcing_heat_tend > 0) then
      call post_data(CS%id_boundary_forcing_heat_tend, work_3d, CS%diag)
    endif
    if(CS%id_boundary_forcing_heat_tend_2d > 0) then
      do j=js,je ; do i=is,ie
        work_2d(i,j) = 0.0
        do k=1,nz
          work_2d(i,j) = work_2d(i,j) + work_3d(i,j,k)
        enddo
      enddo ; enddo
      call post_data(CS%id_boundary_forcing_heat_tend_2d, work_2d, CS%diag)
    endif
  endif


  ! salinity tendency
  if(CS%id_boundary_forcing_saln_tend > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = (tv%S(i,j,k)-saln_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_boundary_forcing_saln_tend, work_3d, CS%diag)
  endif

  ! salt tendency
  if(CS%id_boundary_forcing_salt_tend > 0 .or. CS%id_boundary_forcing_salt_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = GV%H_to_kg_m2 * CS%ppt2mks * Idt * (h(i,j,k) * tv%S(i,j,k) - h_old(i,j,k) * saln_old(i,j,k))
    enddo ; enddo ; enddo
    if(CS%id_boundary_forcing_salt_tend > 0) then
      call post_data(CS%id_boundary_forcing_salt_tend, work_3d, CS%diag)
    endif
    if(CS%id_boundary_forcing_salt_tend_2d > 0) then
      do j=js,je ; do i=is,ie
        work_2d(i,j) = 0.0
        do k=1,nz
          work_2d(i,j) = work_2d(i,j) + work_3d(i,j,k)
        enddo
      enddo ; enddo
      call post_data(CS%id_boundary_forcing_salt_tend_2d, work_2d, CS%diag)
    endif
  endif

end subroutine diagnose_boundary_forcing_tendency


!> This routine diagnoses tendencies for temperature and heat from frazil formation.
!! This routine is called twice from within subroutine diabatic; at start and at
!! end of the diabatic processes. The impacts from frazil are generally a function
!! of depth.  Hence, when checking heat budget, be sure to remove HFSIFRAZIL from HFDS in k=1.
subroutine diagnose_frazil_tendency(tv, h, temp_old, dt, G, GV, CS, ncall)
  type(ocean_grid_type),                    intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type),                  intent(in) :: GV       !< ocean vertical grid structure
  type(thermo_var_ptrs),                    intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h        !< thickness (m or kg/m2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: temp_old !< temperature prior to frazil formation
  real,                                     intent(in) :: dt       !< time step (sec)
  integer,                                  intent(in) :: ncall    !< the first or second call of this routine
  type(diabatic_CS),                        pointer    :: CS       !< module control structure

  real, dimension(SZI_(G),SZJ_(G))         :: work_2d
  real    :: Idt
  integer :: i, j, k, is, ie, js, je, nz

  is  = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Idt = 1/dt
  work_2d(:,:) = 0.0

  ! zero the tendencies at start of first call
  if(ncall == 1) then
    CS%frazil_heat_diag(:,:,:) = 0.0
    CS%frazil_temp_diag(:,:,:) = 0.0
  endif

  ! temperature tendency
  if(CS%id_frazil_temp_tend > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%frazil_temp_diag(i,j,k) = CS%frazil_temp_diag(i,j,k) + Idt * (tv%T(i,j,k)-temp_old(i,j,k))
    enddo ; enddo ; enddo
    if(ncall == 2) then
      call post_data(CS%id_frazil_temp_tend, CS%frazil_temp_diag(:,:,:), CS%diag)
    endif
  endif

  ! heat tendency
  if(CS%id_frazil_heat_tend > 0 .or. CS%id_frazil_heat_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%frazil_heat_diag(i,j,k) = CS%frazil_heat_diag(i,j,k) + &
                                   GV%H_to_kg_m2 * tv%C_p * h(i,j,k) * Idt * (tv%T(i,j,k)-temp_old(i,j,k))
    enddo ; enddo ; enddo
    if(CS%id_frazil_heat_tend  > 0 .and. ncall == 2) then
      call post_data(CS%id_frazil_heat_tend, CS%frazil_heat_diag(:,:,:), CS%diag)
    endif

    ! As a consistency check, we must have
    ! FRAZIL_HEAT_TENDENCY_2d = HFSIFRAZIL
    if(CS%id_frazil_heat_tend_2d > 0 .and. ncall == 2) then
      do j=js,je ; do i=is,ie
        work_2d(i,j) = 0.0
        do k=1,nz
          work_2d(i,j) = work_2d(i,j) + CS%frazil_heat_diag(i,j,k)
        enddo
      enddo ; enddo
      call post_data(CS%id_frazil_heat_tend_2d, work_2d, CS%diag)
    endif
  endif


end subroutine diagnose_frazil_tendency


!> A simplified version of diabatic_driver_init that will allow
!! tracer column functions to be called without allowing any
!! of the diabatic processes to be used.
subroutine adiabatic_driver_init(Time, G, param_file, diag, CS, &
                                tracer_flow_CSp, diag_to_Z_CSp)
  type(time_type),         intent(in)    :: Time              !< current model time
  type(ocean_grid_type),   intent(in)    :: G                 !< model grid structure
  type(param_file_type),   intent(in)    :: param_file        !< the file to parse for parameter values
  type(diag_ctrl), target, intent(inout) :: diag              !< regulates diagnostic output
  type(diabatic_CS),       pointer       :: CS                !< module control structure
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp   !< points to control structure of tracer flow control module
  type(diag_to_Z_CS),      pointer       :: diag_to_Z_CSp     !< pointer to Z-diagnostics control structure

! This "include" declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod  = "MOM_diabatic_driver" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "adiabatic_driver_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag
  if (associated(tracer_flow_CSp)) CS%tracer_flow_CSp => tracer_flow_CSp
  if (associated(diag_to_Z_CSp)) CS%diag_to_Z_CSp => diag_to_Z_CSp

! Set default, read and log parameters
  call log_version(param_file, mod, version, &
                   "The following parameters are used for diabatic processes.")

end subroutine adiabatic_driver_init


!> This routine initializes the diabatic driver module.
subroutine diabatic_driver_init(Time, G, GV, param_file, useALEalgorithm, diag, &
                                ADp, CDp, CS, tracer_flow_CSp, sponge_CSp, &
                                ALE_sponge_CSp, diag_to_Z_CSp)
  type(time_type),         intent(in)    :: Time             !< model time
  type(ocean_grid_type),   intent(inout) :: G                !< model grid structure
  type(verticalGrid_type), intent(in)    :: GV               !< model vertical grid structure
  type(param_file_type),   intent(in)    :: param_file       !< file to parse for parameter values
  logical,                 intent(in)    :: useALEalgorithm  !< logical for whether to use ALE remapping
  type(diag_ctrl), target, intent(inout) :: diag             !< structure to regulate diagnostic output
  type(accel_diag_ptrs),   intent(inout) :: ADp              !< pointers to accelerations in momentum equations,
                                                             !! to enable diagnostics, like energy budgets
  type(cont_diag_ptrs),    intent(inout) :: CDp              !< pointers to terms in continuity equations
  type(diabatic_CS),       pointer       :: CS               !< module control structure
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp  !< pointer to control structure of tracer flow control module
  type(sponge_CS),         pointer       :: sponge_CSp       !< pointer to the sponge module control structure
  type(ALE_sponge_CS),     pointer       :: ALE_sponge_CSp   !< pointer to the ALE sponge module control structure
  type(diag_to_Z_CS),      pointer       :: diag_to_Z_CSp    !< pointer to the Z-diagnostics control structure

  real    :: Kd
  integer :: num_mode
  logical :: use_temperature, differentialDiffusion
  type(vardesc) :: vd

! This "include" declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod  = "MOM_diabatic_driver" ! This module's name.
  character(len=48)  :: thickness_units
  character(len=40)  :: var_name
  character(len=160) :: var_descript
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz, nbands, m
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "diabatic_driver_init called with an "// &
                            "associated control structure.")
    return
  else
    allocate(CS)
  endif

  CS%diag => diag
  if (associated(tracer_flow_CSp)) CS%tracer_flow_CSp => tracer_flow_CSp
  if (associated(sponge_CSp))      CS%sponge_CSp      => sponge_CSp
  if (associated(ALE_sponge_CSp))  CS%ALE_sponge_CSp  => ALE_sponge_CSp
  if (associated(diag_to_Z_CSp))   CS%diag_to_Z_CSp   => diag_to_Z_CSp

  CS%useALEalgorithm = useALEalgorithm
  CS%bulkmixedlayer = (GV%nkml > 0)

  ! Set default, read and log parameters
  call log_version(param_file, mod, version, &
                   "The following parameters are used for diabatic processes.")

  call get_param(param_file, mod, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified via calls to initialize_sponge and possibly \n"//&
                 "set_up_sponge_field.", default=.false.)
  call get_param(param_file, mod, "ENABLE_THERMODYNAMICS", use_temperature, &
                 "If true, temperature and salinity are used as state \n"//&
                 "variables.", default=.true.)
  call get_param(param_file, mod, "ENERGETICS_SFC_PBL", CS%use_energetic_PBL, &
                 "If true, use an implied energetics planetary boundary \n"//&
                 "layer scheme to determine the diffusivity and viscosity \n"//&
                 "in the surface boundary layer.", default=.false.)
  call get_param(param_file, mod, "EPBL_IS_ADDITIVE", CS%ePBL_is_additive, &
                 "If true, the diffusivity from ePBL is added to all\n"//&
                 "other diffusivities. Otherwise, the larger of kappa-\n"//&
                 "shear and ePBL diffusivities are used.", default=.true.)
  call get_param(param_file, mod, "DOUBLE_DIFFUSION", differentialDiffusion, &
                 "If true, apply parameterization of double-diffusion.", &
                 default=.false. )
  CS%use_kappa_shear = kappa_shear_is_used(param_file)
  CS%use_CVMix_shear = cvmix_shear_is_used(param_file)
  if (CS%bulkmixedlayer) then
    call get_param(param_file, mod, "ML_MIX_FIRST", CS%ML_mix_first, &
                 "The fraction of the mixed layer mixing that is applied \n"//&
                 "before interior diapycnal mixing.  0 by default.", &
                 units="nondim", default=0.0)
    call get_param(param_file, mod, "NKBL", CS%nkbl, default=2, do_not_log=.true.)
  else
    CS%ML_mix_first = 0.0
  endif
  if (use_temperature) then
    call get_param(param_file, mod, "DO_GEOTHERMAL", CS%use_geothermal, &
                 "If true, apply geothermal heating.", default=.false.)
  else
    CS%use_geothermal = .false.
  endif
  call get_param(param_file, mod, "INTERNAL_TIDES", CS%use_int_tides, &
                 "If true, use the code that advances a separate set of \n"//&
                 "equations for the internal tide energy density.", default=.false.)
  CS%nMode = 1
  if (CS%use_int_tides) then
    ! SET NUMBER OF MODES TO CONSIDER
    call get_param(param_file, mod, "INTERNAL_TIDE_MODES", CS%nMode, &
                 "The number of distinct internal tide modes \n"//&
                 "that will be calculated.", default=1, do_not_log=.true.)

    ! The following parameters are used in testing the internal tide code.
    ! GET LOCATION AND DURATION OF ENERGY POINT SOURCE FOR TESTING (BDM)
    call get_param(param_file, mod, "INTERNAL_TIDE_SOURCE_TEST", CS%int_tide_source_test, &
                 "If true, apply an arbitrary generation site for internal tide testing", &
                 default=.false.)
    if(CS%int_tide_source_test)then
      call get_param(param_file, mod, "INTERNAL_TIDE_SOURCE_X", CS%int_tide_source_x, &
                 "X Location of generation site for internal tide", default=1.)
      call get_param(param_file, mod, "INTERNAL_TIDE_SOURCE_Y", CS%int_tide_source_y, &
                 "Y Location of generation site for internal tide", default=1.)
      call get_param(param_file, mod, "INTERNAL_TIDE_SOURCE_TLEN_DAYS", CS%tlen_days, &
                 "Time interval from start of experiment for adding wave source", &
                 units="days", default=0)
      CS%time_max_source = increment_time(Time,0,days=CS%tlen_days)
    endif
    ! GET UNIFORM MODE VELOCITY FOR TESTING (BDM)
    call get_param(param_file, mod, "UNIFORM_CG", CS%uniform_cg, &
                 "If true, set cg = cg_test everywhere for test case", default=.false.)
    if(CS%uniform_cg)then
      call get_param(param_file, mod, "CG_TEST", CS%cg_test, &
                 "Uniform group velocity of internal tide for test case", default=1.)
    endif
  endif

  call get_param(param_file, mod, "MASSLESS_MATCH_TARGETS", &
                                CS%massless_match_targets, &
                 "If true, the temperature and salinity of massless layers \n"//&
                 "are kept consistent with their target densities. \n"//&
                 "Otherwise the properties of massless layers evolve \n"//&
                 "diffusively to match massive neighboring layers.", &
                 default=.true.)

  call get_param(param_file, mod, "AGGREGATE_FW_FORCING", CS%aggregate_FW_forcing, &
                 "If true, the net incoming and outgoing fresh water fluxes are combined\n"//&
                 "and applied as either incoming or outgoing depending on the sign of the net.\n"//&
                 "If false, the net incoming fresh water flux is added to the model and\n"//&
                 "thereafter the net outgoing is removed from the updated state."//&
                 "into the first non-vanished layer for which the column remains stable", &
                 default=.true.)

  call get_param(param_file, mod, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)
  call get_param(param_file, mod, "DEBUG_CONSERVATION", CS%debugConservation, &
                 "If true, monitor conservation and extrema.", default=.false.)

  call get_param(param_file, mod, "DEBUG_ENERGY_REQ", CS%debug_energy_req, &
                 "If true, debug the energy requirements.", default=.false., do_not_log=.true.)
  call get_param(param_file, mod, "MIX_BOUNDARY_TRACERS", CS%mix_boundary_tracers, &
                 "If true, mix the passive tracers in massless layers at \n"//&
                 "the bottom into the interior as though a diffusivity of \n"//&
                 "KD_MIN_TR were operating.", default=.true.)

  if (CS%mix_boundary_tracers) then
    call get_param(param_file, mod, "KD", Kd, fail_if_missing=.true.)
    call get_param(param_file, mod, "KD_MIN_TR", CS%Kd_min_tr, &
                 "A minimal diffusivity that should always be applied to \n"//&
                 "tracers, especially in massless layers near the bottom. \n"//&
                 "The default is 0.1*KD.", units="m2 s-1", default=0.1*Kd)
    call get_param(param_file, mod, "KD_BBL_TR", CS%Kd_BBL_tr, &
                 "A bottom boundary layer tracer diffusivity that will \n"//&
                 "allow for explicitly specified bottom fluxes. The \n"//&
                 "entrainment at the bottom is at least sqrt(Kd_BBL_tr*dt) \n"//&
                 "over the same distance.", units="m2 s-1", default=0.)
  endif

  call get_param(param_file, mod, "TRACER_TRIDIAG", CS%tracer_tridiag, &
                 "If true, use the passive tracer tridiagonal solver for T and S\n", &
                 default=.false.)


  ! Register all available diagnostics for this module.
  if (GV%Boussinesq) then ; thickness_units = "meter"
  else ; thickness_units = "kilogram meter-2" ; endif

  CS%id_ea = register_diag_field('ocean_model','ea',diag%axesTL,Time, &
      'Layer entrainment from above per timestep','meter')
  CS%id_eb = register_diag_field('ocean_model','eb',diag%axesTL,Time, &
      'Layer entrainment from below per timestep', 'meter')
  CS%id_dudt_dia = register_diag_field('ocean_model','dudt_dia',diag%axesCuL,Time, &
      'Zonal Acceleration from Diapycnal Mixing', 'meter second-2')
  CS%id_dvdt_dia = register_diag_field('ocean_model','dvdt_dia',diag%axesCvL,Time, &
      'Meridional Acceleration from Diapycnal Mixing', 'meter second-2')
  CS%id_wd = register_diag_field('ocean_model','wd',diag%axesTi,Time, &
      'Diapycnal Velocity', 'meter second-1')
  if (CS%use_int_tides) then
    CS%id_cg1 = register_diag_field('ocean_model','cn1', diag%axesT1, &
                 Time, 'First baroclinic mode (eigen) speed', 'm s-1')
    allocate(CS%id_cn(CS%nMode)) ; CS%id_cn(:) = -1
    do m=1,CS%nMode
      write(var_name, '("cn_mode",i1)') m
      write(var_descript, '("Baroclinic (eigen) speed of mode ",i1)') m
      CS%id_cn(m) = register_diag_field('ocean_model',var_name, diag%axesT1, &
                   Time, var_descript, 'm s-1')
      call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    enddo
  endif

  CS%id_Tdif = register_diag_field('ocean_model',"Tflx_dia_diff",diag%axesTi, &
      Time, "Diffusive diapycnal temperature flux across interfaces", &
      "degC meter second-1")
  CS%id_Tadv = register_diag_field('ocean_model',"Tflx_dia_adv",diag%axesTi, &
      Time, "Advective diapycnal temperature flux across interfaces", &
      "degC meter second-1")
  CS%id_Sdif = register_diag_field('ocean_model',"Sflx_dia_diff",diag%axesTi, &
      Time, "Diffusive diapycnal salnity flux across interfaces", &
      "PSU meter second-1")
  CS%id_Sadv = register_diag_field('ocean_model',"Sflx_dia_adv",diag%axesTi, &
      Time, "Advective diapycnal salnity flux across interfaces", &
      "PSU meter second-1")
  CS%id_MLD_003 = register_diag_field('ocean_model','MLD_003',diag%axesT1,Time,        &
      'Mixed layer depth (delta rho = 0.03)', 'meter', cmor_field_name='mlotst',       &
      cmor_long_name='Ocean Mixed Layer Thickness Defined by Sigma T', cmor_units='m', &
      cmor_standard_name='ocean_mixed_layer_thickness_defined_by_sigma_t')
  CS%id_mlotstsq = register_diag_field('ocean_model','mlotstsq',diag%axesT1,Time,      &
      long_name='Square of Ocean Mixed Layer Thickness Defined by Sigma T',            &
      standard_name='square_of_ocean_mixed_layer_thickness_defined_by_sigma_t',units='m2')
  CS%id_MLD_0125 = register_diag_field('ocean_model','MLD_0125',diag%axesT1,Time, &
      'Mixed layer depth (delta rho = 0.125)', 'meter')
  CS%id_subMLN2  = register_diag_field('ocean_model','subML_N2',diag%axesT1,Time, &
      'Squared buoyancy frequency below mixed layer', 's-2')
  CS%id_MLD_user = register_diag_field('ocean_model','MLD_user',diag%axesT1,Time, &
      'Mixed layer depth (used defined)', 'meter')
  call get_param(param_file, mod, "DIAG_MLD_DENSITY_DIFF", CS%MLDdensityDifference, &
                 "The density difference used to determine a diagnostic mixed\n"//&
                 "layer depth, MLD_user, following the definition of Levitus 1982. \n"//&
                 "The MLD is the depth at which the density is larger than the\n"//&
                 "surface density by the specified amount.", units='kg/m3', default=0.1)

  ! diagnostics making use of the z-gridding code
  if (associated(diag_to_Z_CSp)) then
    vd = var_desc("Kd_interface", "meter2 second-1", &
                  "Diapycnal diffusivity at interfaces, interpolated to z", z_grid='z')
    CS%id_Kd_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = var_desc("Tflx_dia_diff", "degC meter second-1", &
                  "Diffusive diapycnal temperature flux across interfaces, interpolated to z", &
                  z_grid='z')
    CS%id_Tdif_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = var_desc("Tflx_dia_adv", "degC meter second-1", &
                  "Advective diapycnal temperature flux across interfaces, interpolated to z",&
                  z_grid='z')
    CS%id_Tadv_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = var_desc("Sflx_dia_diff", "PSU meter second-1", &
                  "Diffusive diapycnal salinity flux across interfaces, interpolated to z",&
                  z_grid='z')
    CS%id_Sdif_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = var_desc("Sflx_dia_adv", "PSU meter second-1", &
                  "Advective diapycnal salinity flux across interfaces, interpolated to z",&
                  z_grid='z')
    CS%id_Sadv_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
  endif

  if (CS%id_dudt_dia > 0) call safe_alloc_ptr(ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
  if (CS%id_dvdt_dia > 0) call safe_alloc_ptr(ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
  if (CS%id_wd > 0)       call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)

  !call set_diffusivity_init(Time, G, param_file, diag, CS%set_diff_CSp, diag_to_Z_CSp, CS%int_tide_CSp)
  CS%id_Kd_interface = register_diag_field('ocean_model', 'Kd_interface', diag%axesTi, Time, &
      'Total diapycnal diffusivity at interfaces', 'meter2 second-1')
  if (CS%use_energetic_PBL) then
      CS%id_Kd_ePBL = register_diag_field('ocean_model', 'Kd_ePBL', diag%axesTi, Time, &
          'ePBL diapycnal diffusivity at interfaces', 'meter2 second-1')
  endif

  CS%id_Kd_heat = register_diag_field('ocean_model', 'Kd_heat', diag%axesTi, Time, &
      'Total diapycnal diffusivity for heat at interfaces', 'meter2 second-1',     &
       cmor_field_name='difvho', cmor_units='m2 s-1',                              &
       cmor_standard_name='ocean_vertical_heat_diffusivity',                       &
       cmor_long_name='Ocean vertical heat diffusivity')
  CS%id_Kd_salt = register_diag_field('ocean_model', 'Kd_salt', diag%axesTi, Time, &
      'Total diapycnal diffusivity for salt at interfaces', 'meter2 second-1',     &
       cmor_field_name='difvso', cmor_units='m2 s-1',                              &
       cmor_standard_name='ocean_vertical_salt_diffusivity',                       &
       cmor_long_name='Ocean vertical salt diffusivity')

  ! CS%useKPP is set to True if KPP-scheme is to be used, False otherwise.
  ! KPP_init() allocated CS%KPP_Csp and also sets CS%KPPisPassive
  CS%useKPP = KPP_init(param_file, G, diag, Time, CS%KPP_CSp, passive=CS%KPPisPassive)
  if (CS%useKPP .or. CS%use_energetic_PBL) then
    allocate( CS%KPP_NLTheat(isd:ied,jsd:jed,nz+1) )   ; CS%KPP_NLTheat(:,:,:)   = 0.
    allocate( CS%KPP_NLTscalar(isd:ied,jsd:jed,nz+1) ) ; CS%KPP_NLTscalar(:,:,:) = 0.
    allocate( CS%KPP_buoy_flux(isd:ied,jsd:jed,nz+1) ) ; CS%KPP_buoy_flux(:,:,:) = 0.
    allocate( CS%KPP_temp_flux(isd:ied,jsd:jed) )      ; CS%KPP_temp_flux(:,:)   = 0.
    allocate( CS%KPP_salt_flux(isd:ied,jsd:jed) )      ; CS%KPP_salt_flux(:,:)   = 0.
 endif

  call get_param(param_file, mod, "SALT_REJECT_BELOW_ML", CS%salt_reject_below_ML, &
                 "If true, place salt from brine rejection below the mixed layer,\n"// &
                 "into the first non-vanished layer for which the column remains stable", &
                 default=.false.)

  if (CS%salt_reject_below_ML) then
      CS%id_brine_lay = register_diag_field('ocean_model','brine_layer',diag%axesT1,Time, &
      'Brine insertion layer','none')
  endif


  ! diagnostics for tendencies of temp and saln due to diabatic processes;
  ! available only for ALE algorithm.
  if (CS%useALEalgorithm) then
    CS%id_diabatic_diff_temp_tend = register_diag_field('ocean_model', &
        'diabatic_diff_temp_tendency', diag%axesTL, Time,              &
        'Diabatic diffusion temperature tendency', 'Degree C per second')
    if (CS%id_diabatic_diff_temp_tend > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    CS%id_diabatic_diff_saln_tend = register_diag_field('ocean_model',&
        'diabatic_diff_saln_tendency', diag%axesTL, Time,             &
        'Diabatic diffusion salinity tendency', 'PPT per second')
    if (CS%id_diabatic_diff_saln_tend > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    CS%id_diabatic_diff_heat_tend = register_diag_field('ocean_model',                                                 &
        'diabatic_heat_tendency', diag%axesTL, Time,                                                                   &
        'Diabatic diffusion heat tendency',                                                                            &
        'Watts/m2',cmor_field_name='opottempdiff', cmor_units='W m-2',                                                 &
        cmor_standard_name=                                                                                            &
        'tendency_of_sea_water_potential_temperature_expressed_as_heat_content_due_to_parameterized_dianeutral_mixing',&
        cmor_long_name =                                                                                               &
        'Tendency of sea water potential temperature expressed as heat content due to parameterized dianeutral mixing',&
        v_extensive=.true.)
    if (CS%id_diabatic_diff_heat_tend > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    CS%id_diabatic_diff_salt_tend = register_diag_field('ocean_model',                                     &
        'diabatic_salt_tendency', diag%axesTL, Time,                                                       &
        'Diabatic diffusion of salt tendency',                                                             &
        'kg/(m2 * s)',cmor_field_name='osaltdiff', cmor_units='kg m-2 s-1',                                &
        cmor_standard_name=                                                                                &
        'tendency_of_sea_water_salinity_expressed_as_salt_content_due_to_parameterized_dianeutral_mixing', &
        cmor_long_name =                                                                                   &
        'Tendency of sea water salinity expressed as salt content due to parameterized dianeutral mixing', &
        v_extensive=.true.)
    if (CS%id_diabatic_diff_salt_tend > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    ! This diagnostic should equal to roundoff if all is working well.
    CS%id_diabatic_diff_heat_tend_2d = register_diag_field('ocean_model',                                                               &
        'diabatic_heat_tendency_2d', diag%axesT1, Time,                                                                                 &
        'Depth integrated diabatic diffusion heat tendency',                                                                            &
        'Watts/m2',cmor_field_name='opottempdiff_2d', cmor_units='W m-2',                                                               &
        cmor_standard_name=                                                                                                             &
        'tendency_of_sea_water_potential_temperature_expressed_as_heat_content_due_to_parameterized_dianeutral_mixing_depth_integrated',&
        cmor_long_name =                                                                                                                &
        'Tendency of sea water potential temperature expressed as heat content due to parameterized dianeutral mixing depth integrated')
    if (CS%id_diabatic_diff_heat_tend_2d > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    ! This diagnostic should equal to roundoff if all is working well.
    CS%id_diabatic_diff_salt_tend_2d = register_diag_field('ocean_model',                                                  &
        'diabatic_salt_tendency_2d', diag%axesT1, Time,                                                                    &
        'Depth integrated diabatic diffusion salt tendency',                                                               &
        'kg/(m2 * s)',cmor_field_name='osaltdiff_2d', cmor_units='kg m-2 s-1',                                             &
        cmor_standard_name=                                                                                                &
        'tendency_of_sea_water_salinity_expressed_as_salt_content_due_to_parameterized_dianeutral_mixing_depth_integrated',&
        cmor_long_name =                                                                                                   &
        'Tendency of sea water salinity expressed as salt content due to parameterized dianeutral mixing depth integrated')
    if (CS%id_diabatic_diff_salt_tend_2d > 0) then
      CS%diabatic_diff_tendency_diag = .true.
    endif

    ! diagnostics for tendencies of temp and saln due to boundary forcing;
    ! available only for ALE algorithm.
    CS%id_boundary_forcing_temp_tend = register_diag_field('ocean_model',&
        'boundary_forcing_temp_tendency', diag%axesTL, Time,             &
        'Boundary forcing temperature tendency', 'Degree C per second')
    if (CS%id_boundary_forcing_temp_tend > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    CS%id_boundary_forcing_saln_tend = register_diag_field('ocean_model',&
        'boundary_forcing_saln_tendency', diag%axesTL, Time,             &
        'Boundary forcing saln tendency', 'PPT per second')
    if (CS%id_boundary_forcing_saln_tend > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    CS%id_boundary_forcing_heat_tend = register_diag_field('ocean_model',&
        'boundary_forcing_heat_tendency', diag%axesTL, Time,             &
        'Boundary forcing heat tendency','Watts/m2')
    if (CS%id_boundary_forcing_heat_tend > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    CS%id_boundary_forcing_salt_tend = register_diag_field('ocean_model',&
        'boundary_forcing_salt_tendency', diag%axesTL, Time,             &
        'Boundary forcing salt tendency','kg m-2 s-1')
    if (CS%id_boundary_forcing_salt_tend > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    ! This diagnostic should equal to surface heat flux if all is working well.
    CS%id_boundary_forcing_heat_tend_2d = register_diag_field('ocean_model',&
        'boundary_forcing_heat_tendency_2d', diag%axesT1, Time,             &
        'Depth integrated boundary forcing of ocean heat','Watts/m2')
    if (CS%id_boundary_forcing_heat_tend_2d > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif

    ! This diagnostic should equal to surface salt flux if all is working well.
    CS%id_boundary_forcing_salt_tend_2d = register_diag_field('ocean_model',&
        'boundary_forcing_salt_tendency_2d', diag%axesT1, Time,             &
        'Depth integrated boundary forcing of ocean salt','kg m-2 s-1')
    if (CS%id_boundary_forcing_salt_tend_2d > 0) then
      CS%boundary_forcing_tendency_diag = .true.
    endif
  endif

  ! diagnostics for tendencies of temp and heat due to frazil

  ! diagnostic for tendency of temp due to frazil
  CS%id_frazil_temp_tend = register_diag_field('ocean_model',&
      'frazil_temp_tendency', diag%axesTL, Time,             &
      'Temperature tendency due to frazil formation', 'Degree C per second')
  if (CS%id_frazil_temp_tend > 0) then
    CS%frazil_tendency_diag = .true.
  endif

  ! diagnostic for tendency of heat due to frazil
  CS%id_frazil_heat_tend = register_diag_field('ocean_model',&
      'frazil_heat_tendency', diag%axesTL, Time,             &
      'Heat tendency due to frazil formation','Watts/m2')
  if (CS%id_frazil_heat_tend > 0) then
    CS%frazil_tendency_diag = .true.
  endif

  ! if all is working propertly, this diagnostic should equal to hfsifrazil
  CS%id_frazil_heat_tend_2d = register_diag_field('ocean_model',&
      'frazil_heat_tendency_2d', diag%axesT1, Time,             &
      'Depth integrated heat tendency due to frazil formation','Watts/m2')
  if (CS%id_frazil_heat_tend_2d > 0) then
    CS%frazil_tendency_diag = .true.
  endif

  if (CS%frazil_tendency_diag) then
    allocate(CS%frazil_temp_diag(isd:ied,jsd:jed,nz) ) ; CS%frazil_temp_diag(:,:,:) = 0.
    allocate(CS%frazil_heat_diag(isd:ied,jsd:jed,nz) ) ; CS%frazil_heat_diag(:,:,:) = 0.
  endif


  ! CS%useConvection is set to True IF convection will be used, otherwise False.
  ! CS%Conv_CSp is allocated by diffConvection_init()
  CS%useConvection = diffConvection_init(param_file, G, diag, Time, CS%Conv_CSp)

  call entrain_diffusive_init(Time, G, GV, param_file, diag, CS%entrain_diffusive_CSp)

  ! initialize the geothermal heating module
  if (CS%use_geothermal) &
    call geothermal_init(Time, G, param_file, diag, CS%geothermal_CSp)

  ! initialize module for internal tide induced mixing
  if (CS%use_int_tides) then
    call int_tide_input_init(Time, G, GV, param_file, diag, CS%int_tide_input_CSp, &
                             CS%int_tide_input)
    call internal_tides_init(Time, G, GV, param_file, diag, CS%int_tide_CSp)
  endif

  ! initialize module for setting diffusivities
  call set_diffusivity_init(Time, G, GV, param_file, diag, CS%set_diff_CSp, diag_to_Z_CSp, CS%int_tide_CSp)


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
  id_clock_differential_diff = -1 ; if (differentialDiffusion) &
    id_clock_differential_diff = cpu_clock_id('(Ocean differential diffusion)', grain=CLOCK_ROUTINE)

  ! initialize the auxiliary diabatic driver module
  call diabatic_aux_init(Time, G, GV, param_file, diag, CS%diabatic_aux_CSp, &
                         CS%useALEalgorithm, CS%use_energetic_PBL)

  ! initialize the boundary layer modules
  if (CS%bulkmixedlayer) &
    call bulkmixedlayer_init(Time, G, GV, param_file, diag, CS%bulkmixedlayer_CSp)
  if (CS%use_energetic_PBL) &
    call energetic_PBL_init(Time, G, GV, param_file, diag, CS%energetic_PBL_CSp)

  call regularize_layers_init(Time, G, param_file, diag, CS%regularize_layers_CSp)

  if (CS%debug_energy_req) &
    call diapyc_energy_req_init(Time, G, param_file, diag, CS%diapyc_en_rec_CSp)

  ! obtain information about the number of bands for penetrative shortwave
  if (use_temperature) then
    call get_param(param_file, mod, "PEN_SW_NBANDS", nbands, default=1)
    if (nbands > 0) then
      allocate(CS%optics)
      call opacity_init(Time, G, param_file, diag, CS%tracer_flow_CSp, CS%opacity_CSp, CS%optics)
    endif
  endif
  CS%nsw = 0
  if (ASSOCIATED(CS%optics)) CS%nsw = CS%optics%nbands


end subroutine diabatic_driver_init


!> Routine to close the diabatic driver module
subroutine diabatic_driver_end(CS)
  type(diabatic_CS), pointer :: CS    !< module control structure

  if (.not.associated(CS)) return

  call diabatic_aux_end(CS%diabatic_aux_CSp)

  call entrain_diffusive_end(CS%entrain_diffusive_CSp)
  call set_diffusivity_end(CS%set_diff_CSp)
  if (CS%useKPP .or. CS%use_energetic_PBL) then
    deallocate( CS%KPP_NLTheat )
    deallocate( CS%KPP_NLTscalar )
    deallocate( CS%KPP_buoy_flux )
    deallocate( CS%KPP_temp_flux )
    deallocate( CS%KPP_salt_flux )
    call KPP_end(CS%KPP_CSp)
  endif
  if (CS%useConvection) call diffConvection_end(CS%Conv_CSp)
  if (CS%use_energetic_PBL) &
    call energetic_PBL_end(CS%energetic_PBL_CSp)
  if (CS%debug_energy_req) &
    call diapyc_energy_req_end(CS%diapyc_en_rec_CSp)

  if (associated(CS%optics)) then
    call opacity_end(CS%opacity_CSp, CS%optics)
    deallocate(CS%optics)
  endif
  if (associated(CS)) deallocate(CS)

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
!!
!!
!!  \section section_gridlayout MOM grid layout
!!
!!  A small fragment of the grid is shown below:
!!
!! \verbatim
!!    j+1  x ^ x ^ x
!!
!!    j+1  > o > o >
!!
!!    j    x ^ x ^ x
!!
!!    j    > o > o >
!!
!!    j-1  x ^ x ^ x
!!
!!        i-1  i  i+1
!!
!!           i  i+1
!!
!! \endverbatim
!!
!!  Fields at each point
!!  * x =  q, CoriolisBu
!!  * ^ =  v, PFv, CAv, vh, diffv, tauy, vbt, vhtr
!!  * > =  u, PFu, CAu, uh, diffu, taux, ubt, uhtr
!!  * o =  h, bathyT, eta, T, S, tr
!!
!!  The boundaries always run through q grid points (x).
!!

end module MOM_diabatic_driver
