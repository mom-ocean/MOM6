!> This routine drives the diabatic/dianeutral physics for MOM.
!! This is a legacy module that will be deleted in the near future.
module MOM_legacy_diabatic_driver

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
use MOM_diag_mediator,       only : diag_ctrl, query_averaging_enabled, enable_averaging, disable_averaging
use MOM_diag_mediator,       only : diag_grid_storage, diag_grid_storage_init, diag_grid_storage_end
use MOM_diag_mediator,       only : diag_copy_diag_to_storage, diag_copy_storage_to_diag
use MOM_diag_mediator,       only : diag_save_grids, diag_restore_grids
use MOM_diag_to_Z,           only : diag_to_Z_CS, register_Zint_diag, calc_Zint_diags
use MOM_diapyc_energy_req,   only : diapyc_energy_req_init, diapyc_energy_req_end
use MOM_diapyc_energy_req,   only : diapyc_energy_req_calc, diapyc_energy_req_test, diapyc_energy_req_CS
use MOM_CVMix_conv,          only : CVMix_conv_init, CVMix_conv_cs
use MOM_CVMix_conv,          only : CVMix_conv_end, calculate_CVMix_conv
use MOM_domains,             only : pass_var, To_West, To_South, To_All, Omit_Corners
use MOM_domains,             only : create_group_pass, do_group_pass, group_pass_type
use MOM_tidal_mixing,        only : tidal_mixing_init, tidal_mixing_cs
use MOM_tidal_mixing,        only : tidal_mixing_end
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
use MOM_interface_heights,   only : find_eta
use MOM_internal_tides,      only : propagate_int_tide
use MOM_internal_tides,      only : internal_tides_init, internal_tides_end, int_tide_CS
use MOM_kappa_shear,         only : kappa_shear_is_used
use MOM_KPP,                 only : KPP_CS, KPP_init, KPP_compute_BLD, KPP_calculate
use MOM_KPP,                 only : KPP_end, KPP_get_BLD
use MOM_KPP,                 only : KPP_NonLocalTransport_temp, KPP_NonLocalTransport_saln
use MOM_opacity,             only : opacity_init, set_opacity, opacity_end, opacity_CS
use MOM_regularize_layers,   only : regularize_layers, regularize_layers_init, regularize_layers_CS
use MOM_set_diffusivity,     only : set_diffusivity, set_BBL_TKE
use MOM_set_diffusivity,     only : set_diffusivity_init, set_diffusivity_end
use MOM_set_diffusivity,     only : set_diffusivity_CS
use MOM_shortwave_abs,       only : absorbRemainingSW, optics_type
use MOM_sponge,              only : apply_sponge, sponge_CS
use MOM_ALE_sponge,          only : apply_ALE_sponge, ALE_sponge_CS
use MOM_time_manager,        only : operator(-), set_time
use MOM_time_manager,        only : operator(<=), time_type ! for testing itides (BDM)
use MOM_tracer_flow_control, only : call_tracer_column_fns, tracer_flow_control_CS
use MOM_tracer_diabatic,     only : tracer_vertdiff
use MOM_variables,           only : thermo_var_ptrs, vertvisc_type, accel_diag_ptrs
use MOM_variables,           only : cont_diag_ptrs, MOM_thermovar_chksum, p3d
use MOM_verticalGrid,        only : verticalGrid_type
use MOM_wave_speed,          only : wave_speeds
use time_manager_mod,        only : increment_time ! for testing itides (BDM)
use MOM_wave_interface,      only : wave_parameters_CS
use MOM_diabatic_driver,     only : diabatic_CS

implicit none ; private

#include <MOM_memory.h>

public legacy_diabatic

! clock ids
integer :: id_clock_entrain, id_clock_mixedlayer, id_clock_set_diffusivity
integer :: id_clock_tracers, id_clock_tridiag, id_clock_pass, id_clock_sponge
integer :: id_clock_geothermal, id_clock_differential_diff, id_clock_remap
integer :: id_clock_kpp

contains

!>  This subroutine imposes the diapycnal mass fluxes and the
!!  accompanying diapycnal advection of momentum and tracers.
subroutine legacy_diabatic(u, v, h, tv, Hml, fluxes, visc, ADp, CDp, dt, Time_end, &
                    G, GV, CS, WAVES)
  type(ocean_grid_type),                     intent(inout) :: G         !< ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV        !< ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u         !< zonal velocity (m/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v         !< meridional velocity (m/s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h         !< thickness (m for Bouss / kg/m2 for non-Bouss)
  type(thermo_var_ptrs),                     intent(inout) :: tv        !< points to thermodynamic fields
                                                                        !! unused have NULL ptrs
  real, dimension(:,:),                      pointer       :: Hml       !< active mixed layer depth
  type(forcing),                             intent(inout) :: fluxes    !< points to forcing fields
                                                                        !! unused fields have NULL ptrs
  type(vertvisc_type),                       intent(inout) :: visc      !< vertical viscosities, BBL properies, and
  type(accel_diag_ptrs),                     intent(inout) :: ADp       !< related points to accelerations in momentum
                                                                        !! equations, to enable the later derived
                                                                        !! diagnostics, like energy budgets
  type(cont_diag_ptrs),                      intent(inout) :: CDp       !< points to terms in continuity equations
  real,                                      intent(in)    :: dt        !< time increment (seconds)
  type(time_type),                           intent(in)    :: Time_end  !< Time at the end of the interval
  type(diabatic_CS),                         pointer       :: CS        !< module control structure
  type(Wave_parameters_CS),        optional, pointer       :: Waves     !< Surface gravity waves

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
    Rcv_ml, &   ! coordinate density of mixed layer, used for applying sponges
    SkinBuoyFlux! 2d surface buoyancy flux (m2/s3), used by ePBL
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
    eta, &      ! Interface heights before diapycnal mixing, in m.
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
  integer :: dir_flag     ! An integer encoding the directions in which to do halo updates.
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


  ! Offer diagnostics of various state varables at the start of diabatic
  ! these are mostly for debugging purposes.
  if (CS%id_u_predia > 0) call post_data(CS%id_u_predia, u, CS%diag)
  if (CS%id_v_predia > 0) call post_data(CS%id_v_predia, v, CS%diag)
  if (CS%id_h_predia > 0) call post_data(CS%id_h_predia, h, CS%diag)
  if (CS%id_T_predia > 0) call post_data(CS%id_T_predia, tv%T, CS%diag)
  if (CS%id_S_predia > 0) call post_data(CS%id_S_predia, tv%S, CS%diag)
  if (CS%id_e_predia > 0) then
    call find_eta(h, tv, GV%g_Earth, G, GV, eta)
    call post_data(CS%id_e_predia, eta, CS%diag)
  endif


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
  if (associated(tv%T) .AND. associated(tv%frazil)) then
    ! For frazil diagnostic, the first call covers the first half of the time step
    call enable_averaging(0.5*dt, Time_end - set_time(int(floor(0.5*dt+0.5))), CS%diag)
    if (CS%frazil_tendency_diag) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        temp_diag(i,j,k) = tv%T(i,j,k)
      enddo ; enddo ; enddo
    endif

    if (associated(fluxes%p_surf_full)) then
      call make_frazil(h, tv, G, GV, CS%diabatic_aux_CSp, fluxes%p_surf_full)
    else
      call make_frazil(h, tv, G, GV, CS%diabatic_aux_CSp)
    endif
    if (showCallTree) call callTree_waypoint("done with 1st make_frazil (diabatic)")

    if (CS%frazil_tendency_diag) then
      call diagnose_frazil_tendency(tv, h, temp_diag, 0.5*dt, G, GV, CS)
      if (CS%id_frazil_h > 0) call post_data(CS%id_frazil_h, h, CS%diag)
    endif
    call disable_averaging(CS%diag)
  endif
  ! For all other diabatic subroutines, the averaging window should be the entire diabatic timestep
  call enable_averaging(dt, Time_end, CS%diag)
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
                            Hml, CS%aggregate_FW_forcing, dt, last_call=.false.)
        if (CS%salt_reject_below_ML) &
          call insert_brine(h, tv, G, GV, fluxes, nkmb, CS%diabatic_aux_CSp, &
                            dt*CS%ML_mix_first, CS%id_brine_lay)
      else
        ! Changes: h, tv%T, tv%S, eaml and ebml  (G is also inout???)
        call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt, eaml, ebml, &
                        G, GV, CS%bulkmixedlayer_CSp, CS%optics, &
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
        call hchksum(eaml, "after find_uv_at_h eaml",G%HI, scale=GV%H_to_m)
        call hchksum(ebml, "after find_uv_at_h ebml",G%HI, scale=GV%H_to_m)
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
  ! Also changes: visc%Kd_shear, visc%TKE_turb (not clear that TKE_turb is used as input ????
  ! And sets visc%Kv_shear
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

    call KPP_compute_BLD(CS%KPP_CSp, G, GV, h, tv%T, tv%S, u, v, tv%eqn_of_state, &
      fluxes%ustar, CS%KPP_buoy_flux)

    call KPP_calculate(CS%KPP_CSp, G, GV, h, tv%T, tv%S, u, v, tv%eqn_of_state, &
      fluxes%ustar, CS%KPP_buoy_flux, Kd_heat, Kd_salt, visc%Kv_shear, CS%KPP_NLTheat, &
      CS%KPP_NLTscalar, Waves=Waves)
!$OMP parallel default(none) shared(is,ie,js,je,nz,Kd_salt,Kd_int,visc,CS,G,Kd_heat,Hml)

    if (associated(Hml)) then
      call KPP_get_BLD(CS%KPP_CSp, Hml(:,:), G)
      call pass_var(Hml, G%domain, halo=1)
    endif

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

  ! Add vertical diff./visc. due to convection (computed via CVMix)
  if (CS%use_CVMix_conv) then
    call calculate_CVMix_conv(h, tv, G, GV, CS%CVMix_conv_csp, Hml)

      !!!!!!!! GMM, the following needs to be checked !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k=1,nz ; do j=js,je ; do i=is,ie
        Kd_int(i,j,k) = Kd_int(i,j,k) + CS%CVMix_conv_csp%kd_conv(i,j,k)
        visc%Kv_slow(i,j,k) = visc%Kv_slow(i,j,k) + CS%CVMix_conv_csp%kv_conv(i,j,k)
      enddo ; enddo ; enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  endif

  if (CS%useKPP) then

    call cpu_clock_begin(id_clock_kpp)
    if (CS%debug) then
      call hchksum(CS%KPP_temp_flux, "before KPP_applyNLT netHeat",G%HI,haloshift=0, scale=GV%H_to_m)
      call hchksum(CS%KPP_salt_flux, "before KPP_applyNLT netSalt",G%HI,haloshift=0, scale=GV%H_to_m)
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
    if (.not. CS%useKPP) then
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
    call hchksum(ea, "after calc_entrain ea", G%HI, haloshift=0, scale=GV%H_to_m)
    call hchksum(eb, "after calc_entrain eb", G%HI, haloshift=0, scale=GV%H_to_m)
  endif

  ! Save fields before boundary forcing is applied for tendency diagnostics
  if (CS%boundary_forcing_tendency_diag) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h_diag(i,j,k)    = h(i,j,k)
      temp_diag(i,j,k) = tv%T(i,j,k)
      saln_diag(i,j,k) = tv%S(i,j,k)
    enddo ; enddo ; enddo
  endif

  ! Apply forcing when using the ALE algorithm
  if (CS%useALEalgorithm) then
    call cpu_clock_begin(id_clock_remap)

    ! Changes made to following fields:  h, tv%T and tv%S.

    do k=1,nz ; do j=js,je ; do i=is,ie
        h_prebound(i,j,k) = h(i,j,k)
    enddo ; enddo ; enddo
    if (CS%use_energetic_PBL) then

      skinbuoyflux(:,:) = 0.0
      call applyBoundaryFluxesInOut(CS%diabatic_aux_CSp, G, GV, dt, fluxes, CS%optics, &
              h, tv, CS%aggregate_FW_forcing, CS%evap_CFL_limit,                         &
              CS%minimum_forcing_depth, cTKE, dSV_dT, dSV_dS, SkinBuoyFlux=SkinBuoyFlux)

      if (CS%debug) then
        call hchksum(ea, "after applyBoundaryFluxes ea",G%HI,haloshift=0, scale=GV%H_to_m)
        call hchksum(eb, "after applyBoundaryFluxes eb",G%HI,haloshift=0, scale=GV%H_to_m)
        call hchksum(cTKE, "after applyBoundaryFluxes cTKE",G%HI,haloshift=0)
        call hchksum(dSV_dT, "after applyBoundaryFluxes dSV_dT",G%HI,haloshift=0)
        call hchksum(dSV_dS, "after applyBoundaryFluxes dSV_dS",G%HI,haloshift=0)
      endif

      call find_uv_at_h(u, v, h, u_h, v_h, G, GV)
      call energetic_PBL(h, u_h, v_h, tv, fluxes, dt, Kd_ePBL, G, GV, &
           CS%energetic_PBL_CSp, dSV_dT, dSV_dS, cTKE, SkinBuoyFlux, waves=waves)

      ! If visc%MLD exists, copy the ePBL's MLD into it
      if (associated(visc%MLD)) then
        call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, visc%MLD, G)
        call pass_var(visc%MLD, G%domain, halo=1)
        Hml(:,:) = visc%MLD(:,:)
      endif

      ! Augment the diffusivities due to those diagnosed in energetic_PBL.
      do K=2,nz ; do j=js,je ; do i=is,ie

        if (CS%ePBL_is_additive) then
          Kd_add_here = Kd_ePBL(i,j,K)
          visc%Kv_shear(i,j,K) = visc%Kv_shear(i,j,K) + Kd_ePBL(i,j,K)
        else
          Kd_add_here = max(Kd_ePBL(i,j,K) - visc%Kd_shear(i,j,K), 0.0)
          visc%Kv_shear(i,j,K) = max(visc%Kv_shear(i,j,K), Kd_ePBL(i,j,K))
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
        call hchksum(ea, "after ePBL ea",G%HI,haloshift=0, scale=GV%H_to_m)
        call hchksum(eb, "after ePBL eb",G%HI,haloshift=0, scale=GV%H_to_m)
        call hchksum(Kd_ePBL, "after ePBL Kd_ePBL",G%HI,haloshift=0)
      endif

    else
      call applyBoundaryFluxesInOut(CS%diabatic_aux_CSp, G, GV, dt, fluxes, CS%optics, &
                                    h, tv, CS%aggregate_FW_forcing, &
                                    CS%evap_CFL_limit, CS%minimum_forcing_depth)

    endif   ! endif for CS%use_energetic_PBL

    ! diagnose the tendencies due to boundary forcing
    ! At this point, the diagnostic grids have not been updated since the call to the boundary layer scheme
    !  so all tendency diagnostics need to be posted on h_diag, and grids rebuilt afterwards
    if (CS%boundary_forcing_tendency_diag) then
      call diagnose_boundary_forcing_tendency(tv, h, temp_diag, saln_diag, h_diag, dt, G, GV, CS)
      if (CS%id_boundary_forcing_h > 0) call post_data(CS%id_boundary_forcing_h, h, CS%diag, alt_h = h_diag)
    endif
    ! Boundary fluxes may have changed T, S, and h
    call diag_update_remap_grids(CS%diag)

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
  !$OMP parallel do default(shared)
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
  ! Checks for negative thickness may have changed layer thicknesses
  call diag_update_remap_grids(CS%diag)

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
    if (CS%debugConservation) call MOM_state_stats('BML tridiag', u, v, h, tv%T, tv%S, G)

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
        call hchksum(ea, "after ea = ea + eaml",G%HI,haloshift=0, scale=GV%H_to_m)
        call hchksum(eb, "after eb = eb + ebml",G%HI,haloshift=0, scale=GV%H_to_m)
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
                      Hml, CS%aggregate_FW_forcing, dt, last_call=.true.)

      if (CS%salt_reject_below_ML) &
        call insert_brine(h, tv, G, GV, fluxes, nkmb, CS%diabatic_aux_CSp, dt_mix, &
                          CS%id_brine_lay)

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
      if (CS%debugConservation) call MOM_state_stats('2nd bulkmixedlayer', u, v, h, tv%T, tv%S, G)
    endif

  else  ! following block for when NOT using BULKMIXEDLAYER


    ! calculate change in temperature & salinity due to dia-coordinate surface diffusion
    if (associated(tv%T)) then

      if (CS%debug) then
        call hchksum(ea, "before triDiagTS ea ",G%HI,haloshift=0, scale=GV%H_to_m)
        call hchksum(eb, "before triDiagTS eb ",G%HI,haloshift=0, scale=GV%H_to_m)
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
      ! the bulk mixed layer scheme. Otherwise in ALE-mode, layer thicknesses will have changed
      ! In either case, tendencies should be posted on hold
      if (CS%diabatic_diff_tendency_diag) then
        call diagnose_diabatic_diff_tendency(tv, hold, temp_diag, saln_diag, dt, G, GV, CS)
        if (CS%id_diabatic_diff_h > 0) call post_data(CS%id_diabatic_diff_h, hold, CS%diag, alt_h = hold)
      endif

      call cpu_clock_end(id_clock_tridiag)
      if (showCallTree) call callTree_waypoint("done with triDiagTS (diabatic)")

    endif  ! endif corresponding to if (associated(tv%T))
    if (CS%debugConservation) call MOM_state_stats('triDiagTS', u, v, h, tv%T, tv%S, G)


  endif  ! endif for the BULKMIXEDLAYER block


  if (CS%debug) then
    call MOM_state_chksum("after mixed layer ", u, v, h, G, GV, haloshift=0)
    call MOM_thermovar_chksum("after mixed layer ", tv, G)
    call hchksum(ea, "after mixed layer ea", G%HI, scale=GV%H_to_m)
    call hchksum(eb, "after mixed layer eb", G%HI, scale=GV%H_to_m)
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
    !$OMP parallel do default(shared)
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
    Tr_ea_BBL = sqrt(dt*CS%Kd_BBL_tr)
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
      call call_tracer_column_fns(h_prebound, h, ea, eb, fluxes, Hml, dt, G, GV, tv, &
                                CS%optics, CS%tracer_flow_CSp, CS%debug, &
                                evap_CFL_limit = CS%evap_CFL_limit, &
                                minimum_forcing_depth = CS%minimum_forcing_depth)
    else
      call call_tracer_column_fns(hold, h, eatr, ebtr, fluxes, Hml, dt, G, GV, tv, &
                                CS%optics, CS%tracer_flow_CSp, CS%debug)
    endif

  elseif (associated(visc%Kd_extra_S)) then  ! extra diffusivity for passive tracers

    do j=js,je ; do i=is,ie
      ebtr(i,j,nz) = eb(i,j,nz) ; eatr(i,j,1) = ea(i,j,1)
    enddo ; enddo
    !$OMP parallel do default(shared) private(add_ent)
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
      call call_tracer_column_fns(h_prebound, h, eatr, ebtr, fluxes, Hml, dt, G, GV, tv, &
                                  CS%optics, CS%tracer_flow_CSp, CS%debug,&
                                  evap_CFL_limit = CS%evap_CFL_limit, &
                                  minimum_forcing_depth = CS%minimum_forcing_depth)
    else
      call call_tracer_column_fns(hold, h, eatr, ebtr, fluxes, Hml, dt, G, GV, tv, &
                                  CS%optics, CS%tracer_flow_CSp, CS%debug)
    endif

  else
    if (CS%useALEalgorithm) then
    ! For passive tracers, the changes in thickness due to boundary fluxes has yet to be applied
      call call_tracer_column_fns(h_prebound, h, eatr, ebtr, fluxes, Hml, dt, G, GV, tv, &
                                  CS%optics, CS%tracer_flow_CSp, CS%debug, &
                                  evap_CFL_limit = CS%evap_CFL_limit, &
                                  minimum_forcing_depth = CS%minimum_forcing_depth)
    else
      call call_tracer_column_fns(hold, h, ea, eb, fluxes, Hml, dt, G, GV, tv, &
                                  CS%optics, CS%tracer_flow_CSp, CS%debug)
    endif

  endif  ! (CS%mix_boundary_tracers)



  call cpu_clock_end(id_clock_tracers)


  ! sponges
  if (CS%use_sponge) then
    call cpu_clock_begin(id_clock_sponge)
    if (associated(CS%ALE_sponge_CSp)) then
      ! ALE sponge
      call apply_ALE_sponge(h, dt, G, CS%ALE_sponge_CSp, CS%Time)
    else
      ! Layer mode sponge
      if (CS%bulkmixedlayer .and. associated(tv%eqn_of_state)) then
        do i=is,ie ; p_ref_cv(i) = tv%P_Ref ; enddo
        !$OMP parallel do default(shared)
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


!   Save the diapycnal mass fluxes as a diagnostic field.
  if (associated(CDp%diapyc_vel)) then
    !$OMP parallel do default(shared)
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
      call hchksum(ea, "before net flux rearrangement ea",G%HI, scale=GV%H_to_m)
      call hchksum(eb, "before net flux rearrangement eb",G%HI, scale=GV%H_to_m)
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
      call hchksum(ea, "after net flux rearrangement ea",G%HI, scale=GV%H_to_m)
      call hchksum(eb, "after net flux rearrangement eb",G%HI, scale=GV%H_to_m)
    endif
  endif

! Initialize halo regions of ea, eb, and hold to default values.
  !$OMP parallel do default(shared)
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

  call cpu_clock_begin(id_clock_pass)
  if (G%symmetric) then ; dir_flag = To_All+Omit_Corners
  else ; dir_flag = To_West+To_South+Omit_Corners ; endif
  call create_group_pass(CS%pass_hold_eb_ea, hold, G%Domain, dir_flag, halo=1)
  call create_group_pass(CS%pass_hold_eb_ea, eb, G%Domain, dir_flag, halo=1)
  call create_group_pass(CS%pass_hold_eb_ea, ea, G%Domain, dir_flag, halo=1)
  call do_group_pass(CS%pass_hold_eb_ea, G%Domain)
  ! visc%Kv_shear is not in the group pass because it has larger vertical extent.
  if (associated(visc%Kv_shear)) &
    call pass_var(visc%Kv_shear, G%Domain, To_All+Omit_Corners, halo=1)
  call cpu_clock_end(id_clock_pass)

  if (.not. CS%useALEalgorithm) then
    !  Use a tridiagonal solver to determine effect of the diapycnal
    !  advection on velocity field. It is assumed that water leaves
    !  or enters the ocean with the surface velocity.
    if (CS%debug) then
      call MOM_state_chksum("before u/v tridiag ", u, v, h, G, GV, haloshift=0)
      call hchksum(ea, "before u/v tridiag ea",G%HI, scale=GV%H_to_m)
      call hchksum(eb, "before u/v tridiag eb",G%HI, scale=GV%H_to_m)
      call hchksum(hold, "before u/v tridiag hold",G%HI, scale=GV%H_to_m)
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
      call MOM_state_chksum("aft 1st loop tridiag ", u, v, h, G, GV, haloshift=0)
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
      call MOM_state_chksum("after u/v tridiag ", u, v, h, G, GV, haloshift=0)
    endif
  endif ! useALEalgorithm

  call disable_averaging(CS%diag)
  ! Frazil formation keeps temperature above the freezing point.
  ! make_frazil is deliberately called at both the beginning and at
  ! the end of the diabatic processes.
  if (associated(tv%T) .AND. associated(tv%frazil)) then
    call enable_averaging(0.5*dt, Time_end, CS%diag)
    if (CS%frazil_tendency_diag) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        temp_diag(i,j,k) = tv%T(i,j,k)
      enddo ; enddo ; enddo
    endif

    if (associated(fluxes%p_surf_full)) then
      call make_frazil(h, tv, G, GV, CS%diabatic_aux_CSp, fluxes%p_surf_full)
    else
      call make_frazil(h, tv, G, GV, CS%diabatic_aux_CSp)
    endif

    if (CS%frazil_tendency_diag) then
      call diagnose_frazil_tendency(tv, h, temp_diag, 0.5*dt, G, GV, CS)
      if (CS%id_frazil_h > 0 ) call post_data(CS%id_frazil_h, h, CS%diag)
    endif

    if (showCallTree) call callTree_waypoint("done with 2nd make_frazil (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('2nd make_frazil', u, v, h, tv%T, tv%S, G)
    call disable_averaging(CS%diag)

  endif  ! endif for frazil

  ! Diagnose the diapycnal diffusivities and other related quantities.
  call enable_averaging(dt, Time_end, CS%diag)

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

  call disable_averaging(CS%diag)

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

end subroutine legacy_diabatic

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
  if (CS%id_diabatic_diff_temp_tend > 0) then
    call post_data(CS%id_diabatic_diff_temp_tend, work_3d, CS%diag, alt_h = h)
  endif

  ! heat tendency
  if (CS%id_diabatic_diff_heat_tend > 0 .or. CS%id_diabatic_diff_heat_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = h(i,j,k) * GV%H_to_kg_m2 * tv%C_p * work_3d(i,j,k)
    enddo ; enddo ; enddo
    if (CS%id_diabatic_diff_heat_tend > 0) then
      call post_data(CS%id_diabatic_diff_heat_tend, work_3d, CS%diag, alt_h = h)
    endif
    if (CS%id_diabatic_diff_heat_tend_2d > 0) then
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
  if (CS%id_diabatic_diff_saln_tend > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = (tv%S(i,j,k)-saln_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_diabatic_diff_saln_tend, work_3d, CS%diag, alt_h = h)
  endif

  ! salt tendency
  if (CS%id_diabatic_diff_salt_tend > 0 .or. CS%id_diabatic_diff_salt_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = h(i,j,k) * GV%H_to_kg_m2 * CS%ppt2mks * work_3d(i,j,k)
    enddo ; enddo ; enddo
    if (CS%id_diabatic_diff_salt_tend > 0) then
      call post_data(CS%id_diabatic_diff_salt_tend, work_3d, CS%diag, alt_h = h)
    endif
    if (CS%id_diabatic_diff_salt_tend_2d > 0) then
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
  type(ocean_grid_type),   intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type), intent(in) :: GV       !< ocean vertical grid structure
  type(thermo_var_ptrs),   intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h        !< thickness after boundary flux application (m or kg/m2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: temp_old !< temperature prior to boundary flux application
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: saln_old !< salinity prior to boundary flux application (PPT)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_old    !< thickness prior to boundary flux application (m or kg/m2)
  real,                    intent(in) :: dt       !< time step (sec)
  type(diabatic_CS),       pointer    :: CS       !< module control structure

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: work_3d
  real, dimension(SZI_(G),SZJ_(G))         :: work_2d
  real    :: Idt
  integer :: i, j, k, is, ie, js, je, nz

  is  = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Idt = 1/dt
  work_3d(:,:,:) = 0.0
  work_2d(:,:)   = 0.0

  ! Thickness tendency
  if (CS%id_boundary_forcing_h_tendency > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = (h(i,j,k) - h_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_boundary_forcing_h_tendency, work_3d, CS%diag, alt_h = h_old)
  endif

  ! temperature tendency
  if (CS%id_boundary_forcing_temp_tend > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = (tv%T(i,j,k)-temp_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_boundary_forcing_temp_tend, work_3d, CS%diag, alt_h = h_old)
  endif

  ! heat tendency
  if (CS%id_boundary_forcing_heat_tend > 0 .or. CS%id_boundary_forcing_heat_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = GV%H_to_kg_m2 * tv%C_p * Idt * (h(i,j,k) * tv%T(i,j,k) - h_old(i,j,k) * temp_old(i,j,k))
    enddo ; enddo ; enddo
    if (CS%id_boundary_forcing_heat_tend > 0) then
      call post_data(CS%id_boundary_forcing_heat_tend, work_3d, CS%diag, alt_h = h_old)
    endif
    if (CS%id_boundary_forcing_heat_tend_2d > 0) then
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
  if (CS%id_boundary_forcing_saln_tend > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = (tv%S(i,j,k)-saln_old(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(CS%id_boundary_forcing_saln_tend, work_3d, CS%diag, alt_h = h_old)
  endif

  ! salt tendency
  if (CS%id_boundary_forcing_salt_tend > 0 .or. CS%id_boundary_forcing_salt_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = GV%H_to_kg_m2 * CS%ppt2mks * Idt * (h(i,j,k) * tv%S(i,j,k) - h_old(i,j,k) * saln_old(i,j,k))
    enddo ; enddo ; enddo
    if (CS%id_boundary_forcing_salt_tend > 0) then
      call post_data(CS%id_boundary_forcing_salt_tend, work_3d, CS%diag, alt_h = h_old)
    endif
    if (CS%id_boundary_forcing_salt_tend_2d > 0) then
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
subroutine diagnose_frazil_tendency(tv, h, temp_old, dt, G, GV, CS)
  type(ocean_grid_type),                    intent(in) :: G        !< ocean grid structure
  type(verticalGrid_type),                  intent(in) :: GV       !< ocean vertical grid structure
  type(diabatic_CS),                        pointer    :: CS       !< module control structure
  type(thermo_var_ptrs),                    intent(in) :: tv       !< points to updated thermodynamic fields
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h        !< thickness (m or kg/m2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: temp_old !< temperature prior to frazil formation
  real,                                     intent(in) :: dt       !< time step (sec)

  real, dimension(SZI_(G),SZJ_(G))         :: work_2d
  real    :: Idt
  integer :: i, j, k, is, ie, js, je, nz

  is  = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Idt = 1/dt

  ! temperature tendency
  if (CS%id_frazil_temp_tend > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%frazil_temp_diag(i,j,k) = Idt * (tv%T(i,j,k)-temp_old(i,j,k))
    enddo ; enddo ; enddo
    call post_data(CS%id_frazil_temp_tend, CS%frazil_temp_diag(:,:,:), CS%diag)
  endif

  ! heat tendency
  if (CS%id_frazil_heat_tend > 0 .or. CS%id_frazil_heat_tend_2d > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%frazil_heat_diag(i,j,k) = GV%H_to_kg_m2 * tv%C_p * h(i,j,k) * Idt * (tv%T(i,j,k)-temp_old(i,j,k))
    enddo ; enddo ; enddo
    if (CS%id_frazil_heat_tend > 0) call post_data(CS%id_frazil_heat_tend, CS%frazil_heat_diag(:,:,:), CS%diag)

    ! As a consistency check, we must have
    ! FRAZIL_HEAT_TENDENCY_2d = HFSIFRAZIL
    if (CS%id_frazil_heat_tend_2d > 0) then
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

end module MOM_legacy_diabatic_driver
