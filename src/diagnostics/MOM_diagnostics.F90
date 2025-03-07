!> Calculates any requested diagnostic quantities
!! that are not calculated in the various subroutines.
!! Diagnostic quantities are requested by allocating them memory.
module MOM_diagnostics

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,              only : reproducing_sum
use MOM_coupler_types,     only : coupler_type_send_data
use MOM_density_integrals, only : int_density_dz
use MOM_diag_mediator,     only : post_data, get_diag_time_end
use MOM_diag_mediator,     only : post_product_u, post_product_sum_u
use MOM_diag_mediator,     only : post_product_v, post_product_sum_v
use MOM_diag_mediator,     only : register_diag_field, register_scalar_field
use MOM_diag_mediator,     only : register_static_field, diag_register_area_ids
use MOM_diag_mediator,     only : diag_ctrl, time_type, safe_alloc_ptr
use MOM_diag_mediator,     only : diag_get_volume_cell_measure_dm_id
use MOM_diag_mediator,     only : diag_grid_storage
use MOM_diag_mediator,     only : diag_save_grids, diag_restore_grids, diag_copy_storage_to_diag
use MOM_domains,           only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,           only : To_North, To_East
use MOM_EOS,               only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_EOS,               only : cons_temp_to_pot_temp, abs_saln_to_prac_saln
use MOM_error_handler,     only : MOM_error, FATAL, WARNING
use MOM_file_parser,       only : get_param, log_version, param_file_type
use MOM_grid,              only : ocean_grid_type
use MOM_interface_heights, only : find_eta
use MOM_spatial_means,     only : global_area_mean, global_layer_mean
use MOM_spatial_means,     only : global_volume_mean, global_area_integral
use MOM_tracer_registry,   only : tracer_registry_type, post_tracer_transport_diagnostics
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs, ocean_internal_state, p3d
use MOM_variables,         only : accel_diag_ptrs, cont_diag_ptrs, surface
use MOM_verticalGrid,      only : verticalGrid_type, get_thickness_units, get_flux_units
use MOM_wave_speed,        only : wave_speed, wave_speed_CS, wave_speed_init

implicit none ; private

#include <MOM_memory.h>

public calculate_diagnostic_fields, register_time_deriv, write_static_fields
public find_eta
public register_surface_diags, post_surface_dyn_diags, post_surface_thermo_diags
public register_transport_diags, post_transport_diagnostics
public MOM_diagnostics_init, MOM_diagnostics_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure for the MOM_diagnostics module
type, public :: diagnostics_CS ; private
  logical :: initialized = .false.     !< True if this control structure has been initialized.
  real :: mono_N2_column_fraction = 0. !< The lower fraction of water column over which N2 is limited as
                                       !! monotonic for the purposes of calculating the equivalent
                                       !! barotropic wave speed [nondim].
  real :: mono_N2_depth = -1.          !< The depth below which N2 is limited as monotonic for the purposes of
                                       !! calculating the equivalent barotropic wave speed [H ~> m or kg m-2].

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                       !! regulate the timing of diagnostic output.

  ! following arrays store diagnostics calculated here and unavailable outside.

  ! following fields have nz layers.
  real, allocatable :: du_dt(:,:,:) !< net i-acceleration [L T-2 ~> m s-2]
  real, allocatable :: dv_dt(:,:,:) !< net j-acceleration [L T-2 ~> m s-2]
  real, allocatable :: dh_dt(:,:,:) !< thickness rate of change [H T-1 ~> m s-1 or kg m-2 s-1]

  logical :: KE_term_on !< If true, at least one diagnostic term in the KE budget is in use.

  !>@{ Diagnostic IDs
  integer :: id_u   = -1,   id_v   = -1, id_h = -1
  integer :: id_usq = -1,   id_vsq = -1, id_uv = -1
  integer :: id_e              = -1, id_e_D            = -1
  integer :: id_du_dt          = -1, id_dv_dt          = -1
  ! integer :: id_hf_du_dt       = -1, id_hf_dv_dt       = -1
  integer :: id_h_du_dt       = -1, id_h_dv_dt       = -1
  integer :: id_hf_du_dt_2d    = -1, id_hf_dv_dt_2d    = -1
  integer :: id_col_ht         = -1, id_dh_dt          = -1
  integer :: id_KE             = -1, id_dKEdt          = -1
  integer :: id_PE_to_KE       = -1, id_KE_BT          = -1
  integer :: id_KE_SAL         = -1, id_KE_TIDES       = -1
  integer :: id_KE_BT_PF       = -1, id_KE_BT_CF       = -1
  integer :: id_KE_BT_WD       = -1
  integer :: id_PE_to_KE_btbc  = -1, id_KE_Coradv_btbc = -1
  integer :: id_KE_Coradv      = -1, id_KE_adv         = -1
  integer :: id_KE_visc        = -1, id_KE_stress      = -1
  integer :: id_KE_visc_gl90   = -1
  integer :: id_KE_horvisc     = -1, id_KE_dia         = -1
  integer :: id_uh_Rlay        = -1, id_vh_Rlay        = -1
  integer :: id_uhGM_Rlay      = -1, id_vhGM_Rlay      = -1
  integer :: id_h_Rlay         = -1, id_Rd1            = -1
  integer :: id_Rml            = -1, id_Rcv            = -1
  integer :: id_cg1            = -1, id_cfl_cg1        = -1
  integer :: id_cfl_cg1_x      = -1, id_cfl_cg1_y      = -1
  integer :: id_cg_ebt         = -1, id_Rd_ebt         = -1
  integer :: id_p_ebt          = -1
  integer :: id_temp_int       = -1, id_salt_int       = -1
  integer :: id_mass_wt        = -1, id_col_mass       = -1
  integer :: id_masscello      = -1, id_masso          = -1
  integer :: id_volcello       = -1
  integer :: id_Tpot           = -1, id_Sprac          = -1
  integer :: id_tob            = -1, id_sob            = -1
  integer :: id_thetaoga       = -1, id_soga           = -1
  integer :: id_bigthetaoga    = -1, id_abssoga        = -1
  integer :: id_sosga          = -1, id_tosga          = -1
  integer :: id_abssosga       = -1, id_bigtosga       = -1
  integer :: id_temp_layer_ave = -1, id_salt_layer_ave = -1
  integer :: id_bigtemp_layer_ave = -1, id_abssalt_layer_ave = -1
  integer :: id_pbo            = -1
  integer :: id_thkcello       = -1, id_rhoinsitu      = -1
  integer :: id_rhopot0        = -1, id_rhopot2        = -1
  integer :: id_drho_dT        = -1, id_drho_dS        = -1
  integer :: id_h_pre_sync     = -1
  integer :: id_tosq           = -1, id_sosq           = -1

  !>@}
  type(wave_speed_CS) :: wave_speed  !< Wave speed control struct

  type(p3d) :: var_ptr(MAX_FIELDS_)  !< pointers to variables used in the calculation
                                     !! of time derivatives
  type(p3d) :: deriv(MAX_FIELDS_)    !< Time derivatives of various fields
  type(p3d) :: prev_val(MAX_FIELDS_) !< Previous values of variables used in the calculation
                                     !! of time derivatives
  !< previous values of variables used in calculation of time derivatives
  integer   :: nlay(MAX_FIELDS_) !< The number of layers in each diagnostics
  integer   :: num_time_deriv = 0 !< The number of time derivative diagnostics

  type(group_pass_type) :: pass_KE_uv !< A handle used for group halo passes

end type diagnostics_CS


!> A structure with diagnostic IDs of the surface and integrated variables
type, public :: surface_diag_IDs ; private
  !>@{ Diagnostic IDs for 2-d surface and bottom flux and state fields
  !Diagnostic IDs for 2-d surface and bottom fields
  integer :: id_zos  = -1, id_zossq  = -1
  integer :: id_volo = -1, id_speed  = -1
  integer :: id_ssh  = -1, id_ssh_ga = -1
  integer :: id_sst  = -1, id_sst_sq = -1, id_sstcon = -1
  integer :: id_sss  = -1, id_sss_sq = -1, id_sssabs = -1
  integer :: id_ssu  = -1, id_ssv    = -1
  integer :: id_ssu_east = -1, id_ssv_north = -1

  ! Diagnostic IDs for  heat and salt flux fields
  integer :: id_fraz         = -1
  integer :: id_salt_deficit = -1
  integer :: id_Heat_PmE     = -1
  integer :: id_intern_heat  = -1
  !>@}
end type surface_diag_IDs


!> A structure with diagnostic IDs of mass transport related diagnostics
type, public :: transport_diag_IDs ; private
  !>@{  Diagnostics for tracer horizontal transport
  integer :: id_uhtr = -1, id_umo = -1, id_umo_2d = -1
  integer :: id_vhtr = -1, id_vmo = -1, id_vmo_2d = -1
  integer :: id_dynamics_h = -1, id_dynamics_h_tendency = -1
  !>@}
end type transport_diag_IDs


contains
!> Diagnostics not more naturally calculated elsewhere are computed here.
subroutine calculate_diagnostic_fields(u, v, h, uh, vh, tv, ADp, CDp, p_surf, &
                                       dt, diag_pre_sync, G, GV, US, CS)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: uh   !< Transport through zonal faces = u*h*dy,
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: vh   !< Transport through meridional faces = v*h*dx,
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various
                                                 !! thermodynamic variables.
  type(accel_diag_ptrs),   intent(in)    :: ADp  !< structure with pointers to
                                                 !! accelerations in momentum equation.
  type(cont_diag_ptrs),    intent(in)    :: CDp  !< structure with pointers to
                                                 !! terms in continuity equation.
  real, dimension(:,:),    pointer       :: p_surf !< A pointer to the surface pressure [R L2 T-2 ~> Pa].
                                                 !! If p_surf is not associated, it is the same
                                                 !! as setting the surface pressure to 0.
  real,                    intent(in)    :: dt   !< The time difference since the last
                                                 !! call to this subroutine [T ~> s].
  type(diag_grid_storage), intent(in)    :: diag_pre_sync !< Target grids from previous timestep
  type(diagnostics_CS),    intent(inout) :: CS   !< Control structure returned by a
                                                 !! previous call to diagnostics_init.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G))  :: uv  ! u x v at h-points          [L2 T-2 ~> m2 s-2]

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb

  real :: eta(SZI_(G),SZJ_(G),SZK_(GV)+1)  ! Interface heights, either relative to a reference
                                           ! geopotential or the seafloor [Z ~> m].
  real :: Rcv(SZI_(G),SZJ_(G),SZK_(GV)) ! Coordinate variable potential density [R ~> kg m-3].
  real :: work_3d(SZI_(G),SZJ_(G),SZK_(GV)) ! A 3-d temporary work array in various units
                                            ! including [nondim] and [H ~> m or kg m-2].
  real :: uh_tmp(SZIB_(G),SZJ_(G),SZK_(GV)) ! A temporary zonal transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vh_tmp(SZI_(G),SZJB_(G),SZK_(GV)) ! A temporary meridional transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: mass_cell(SZI_(G),SZJ_(G))       ! The vertically integrated mass in a grid cell [R Z L2 ~> kg]
  real :: rho_in_situ(SZI_(G))             ! In situ density [R ~> kg m-3]
  real :: cg1(SZI_(G),SZJ_(G))             ! First baroclinic gravity wave speed [L T-1 ~> m s-1]
  real :: Rd1(SZI_(G),SZJ_(G))             ! First baroclinic deformation radius [L ~> m]
  real :: CFL_cg1(SZI_(G),SZJ_(G))         ! CFL for first baroclinic gravity wave speed, either based on the
                                           ! overall grid spacing or just one direction [nondim]


  ! tmp array for surface properties
  real :: pressure_1d(SZI_(G)) ! Temporary array for pressure when calling EOS [R L2 T-2 ~> Pa]
  real :: wt, wt_p ! The fractional weights of two successive values when interpolating from
                   ! a list [nondim], scaled so that wt + wt_p = 1.
  real :: f2_h     ! Squared Coriolis parameter at to h-points [T-2 ~> s-2]
  real :: mag_beta ! Magnitude of the gradient of f [T-1 L-1 ~> s-1 m-1]
  real :: absurdly_small_freq2 ! Frequency squared used to avoid division by 0 [T-2 ~> s-2]

  integer :: k_list

  real, dimension(SZK_(GV)) :: temp_layer_ave ! The average temperature in a layer [C ~> degC]
  real, dimension(SZK_(GV)) :: salt_layer_ave ! The average salinity in a layer [S ~> ppt]
  real :: thetaoga  ! The volume mean potential temperature [C ~> degC]
  real :: soga      ! The volume mean ocean salinity [S ~> ppt]
  real :: masso     ! The total mass of the ocean [R Z L2 ~> kg]
  real :: tosga     ! The area mean sea surface temperature [C ~> degC]
  real :: sosga     ! The area mean sea surface salinity [S ~> ppt]

  is  = G%isc  ; ie   = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq  = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nz  = GV%ke   ; nkmb = GV%nk_rho_varies

  ! This value is roughly (pi / (the age of the universe) )^2.
  absurdly_small_freq2 = 1e-34*US%T_to_s**2

  if (.not. CS%initialized) call MOM_error(FATAL, &
         "calculate_diagnostic_fields: Module must be initialized before used.")

  call calculate_derivs(dt, G, CS)

  if (dt > 0.0) then
    call diag_save_grids(CS%diag)
    call diag_copy_storage_to_diag(CS%diag, diag_pre_sync)

    if (CS%id_h_pre_sync > 0) &
        call post_data(CS%id_h_pre_sync, diag_pre_sync%h_state, CS%diag, alt_h=diag_pre_sync%h_state)

    if (CS%id_du_dt>0) call post_data(CS%id_du_dt, CS%du_dt, CS%diag, alt_h=diag_pre_sync%h_state)

    if (CS%id_dv_dt>0) call post_data(CS%id_dv_dt, CS%dv_dt, CS%diag, alt_h=diag_pre_sync%h_state)

    if (CS%id_dh_dt>0) call post_data(CS%id_dh_dt, CS%dh_dt, CS%diag, alt_h=diag_pre_sync%h_state)

    !! Diagnostics for terms multiplied by fractional thicknesses

    ! 3D diagnostics hf_du(dv)_dt are commented because there is no clarity on proper remapping grid option.
    ! The code is retained for debugging purposes in the future.
    !if (CS%id_hf_du_dt > 0) then
    !  call post_product_u(CS%id_hf_du_dt, CS%du_dt, ADp%diag_hfrac_u, G, nz, CS%diag, alt_h=diag_pre_sync%h_state)
    !if (CS%id_hf_dv_dt > 0) &
    !  call post_product_v(CS%id_hf_dv_dt, CS%dv_dt, ADp%diag_hfrac_v, G, nz, CS%diag, alt_h=diag_pre_sync%h_state)

    if (CS%id_hf_du_dt_2d > 0) &
      call post_product_sum_u(CS%id_hf_du_dt_2d, CS%du_dt, ADp%diag_hfrac_u, G, nz, CS%diag)
    if (CS%id_hf_dv_dt_2d > 0) &
      call post_product_sum_v(CS%id_hf_dv_dt_2d, CS%dv_dt, ADp%diag_hfrac_v, G, nz, CS%diag)

    if (CS%id_h_du_dt > 0) &
      call post_product_u(CS%id_h_du_dt, CS%du_dt, ADp%diag_hu, G, nz, CS%diag)
    if (CS%id_h_dv_dt > 0) &
      call post_product_v(CS%id_h_dv_dt, CS%dv_dt, ADp%diag_hv, G, nz, CS%diag)

    call diag_restore_grids(CS%diag)

    call calculate_energy_diagnostics(u, v, h, uh, vh, ADp, CDp, G, GV, US, CS)
  endif

  ! smg: is the following robust to ALE? It seems a bit opaque.
  ! If the model is NOT in isopycnal mode then nkmb=0. But we need all the
  ! following diagnostics to treat all layers as variable density, so we set
  ! nkmb = nz, on the expectation that loops nkmb+1,nz will not iterate.
  ! This behavior is ANSI F77 but some compiler options can force at least
  ! one iteration that would break the following one-line workaround!
  if (nkmb==0 .and. nz > 1) nkmb = nz

  if (CS%id_u > 0) call post_data(CS%id_u, u, CS%diag)

  if (CS%id_v > 0) call post_data(CS%id_v, v, CS%diag)

  if (CS%id_h > 0) call post_data(CS%id_h, h, CS%diag)

  if (CS%id_usq > 0) call post_product_u(CS%id_usq, u, u, G, nz, CS%diag)

  if (CS%id_vsq > 0) call post_product_v(CS%id_vsq, v, v, G, nz, CS%diag)

  if (CS%id_uv > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      uv(i,j,k) = (0.5*(u(I-1,j,k) + u(I,j,k))) * &
                  (0.5*(v(i,J-1,k) + v(i,J,k)))
    enddo ; enddo ; enddo
    call post_data(CS%id_uv, uv, CS%diag)
  endif

  ! Find the interface heights, relative either to a reference height or to the bottom [Z ~> m].
  if (CS%id_e > 0) then
    call find_eta(h, tv, G, GV, US, eta, dZref=G%Z_ref)
    if (CS%id_e > 0) call post_data(CS%id_e, eta, CS%diag)
    if (CS%id_e_D > 0) then
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        eta(i,j,k) = eta(i,j,k) + (G%bathyT(i,j) + G%Z_ref)
      enddo ; enddo ; enddo
      call post_data(CS%id_e_D, eta, CS%diag)
    endif
  elseif (CS%id_e_D > 0) then
    call find_eta(h, tv, G, GV, US, eta)
    do k=1,nz+1 ; do j=js,je ; do i=is,ie
      eta(i,j,k) = eta(i,j,k) + G%bathyT(i,j)
    enddo ; enddo ; enddo
    call post_data(CS%id_e_D, eta, CS%diag)
  endif

  ! mass per area of grid cell (for Boussinesq, use Rho0)
  if (CS%id_masscello > 0) then
    call post_data(CS%id_masscello, h, CS%diag)
  endif

  ! mass of liquid ocean (for Bouss, use Rho0). The reproducing sum requires the use of MKS units.
  if (CS%id_masso > 0) then
    mass_cell(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      mass_cell(i,j) = mass_cell(i,j) + (GV%H_to_RZ*h(i,j,k)) * G%areaT(i,j)
    enddo ; enddo ; enddo
    masso = reproducing_sum(mass_cell, unscale=US%RZL2_to_kg)
    call post_data(CS%id_masso, masso, CS%diag)
  endif

  ! diagnose thickness/volumes of grid cells [Z ~> m] and [m3]
  if (CS%id_thkcello>0 .or. CS%id_volcello>0) then
    if (GV%Boussinesq) then ! thkcello = h for Boussinesq
      if (CS%id_thkcello > 0) then ; if (GV%H_to_Z == 1.0) then
        call post_data(CS%id_thkcello, h, CS%diag)
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          work_3d(i,j,k) = GV%H_to_Z*h(i,j,k)
        enddo ; enddo ; enddo
        call post_data(CS%id_thkcello, work_3d, CS%diag)
      endif ; endif
      if (CS%id_volcello > 0) then ! volcello = h*area for Boussinesq
        do k=1,nz ; do j=js,je ; do i=is,ie
          work_3d(i,j,k) = ( GV%H_to_Z*h(i,j,k) ) * US%Z_to_m*US%L_to_m**2*G%areaT(i,j)
        enddo ; enddo ; enddo
        call post_data(CS%id_volcello, work_3d, CS%diag)
      endif
    else ! thkcello = dp/(rho*g) for non-Boussinesq
      EOSdom(:) = EOS_domain(G%HI)
      do j=js,je
        if (associated(p_surf)) then ! Pressure loading at top of surface layer [R L2 T-2 ~> Pa]
          do i=is,ie
            pressure_1d(i) = p_surf(i,j)
          enddo
        else
          do i=is,ie
            pressure_1d(i) = 0.0
          enddo
        endif
        do k=1,nz ! Integrate vertically downward for pressure
          do i=is,ie ! Pressure for EOS at the layer center [R L2 T-2 ~> Pa]
            pressure_1d(i) = pressure_1d(i) + 0.5*(GV%H_to_RZ*GV%g_Earth)*h(i,j,k)
          enddo
          ! Store in-situ density [R ~> kg m-3] in work_3d
          call calculate_density(tv%T(:,j,k), tv%S(:,j,k),  pressure_1d, rho_in_situ, &
                                 tv%eqn_of_state, EOSdom)
          do i=is,ie ! Cell thickness = dz = dp/(g*rho) (meter); store in work_3d
            work_3d(i,j,k) = (GV%H_to_RZ*h(i,j,k)) / rho_in_situ(i)
          enddo
          do i=is,ie ! Pressure for EOS at the bottom interface [R L2 T-2 ~> Pa]
            pressure_1d(i) = pressure_1d(i) + 0.5*(GV%H_to_RZ*GV%g_Earth)*h(i,j,k)
          enddo
        enddo ! k
      enddo ! j
      if (CS%id_thkcello > 0) call post_data(CS%id_thkcello, work_3d, CS%diag)
      if (CS%id_volcello > 0) then
        do k=1,nz ; do j=js,je ; do i=is,ie ! volcello = dp/(rho*g)*area for non-Boussinesq
          work_3d(i,j,k) = US%Z_to_m*US%L_to_m**2*G%areaT(i,j) * work_3d(i,j,k)
        enddo ; enddo ; enddo
        call post_data(CS%id_volcello, work_3d, CS%diag)
      endif
    endif
  endif

  ! Calculate additional, potentially derived temperature diagnostics
  if (tv%T_is_conT) then
    ! Internal T&S variables are conservative temperature & absolute salinity,
    ! so they need to converted to potential temperature and practical salinity
    ! for some diagnostics using TEOS-10 function calls.
    if ((CS%id_Tpot > 0) .or. (CS%id_tob > 0) .or. (CS%id_tosq > 0)) then
      EOSdom(:) = EOS_domain(G%HI)
      do k=1,nz ; do j=js,je
        call cons_temp_to_pot_temp(tv%T(:,j,k), tv%S(:,j,k), work_3d(:,j,k), tv%eqn_of_state, EOSdom)
      enddo ; enddo
      if (CS%id_Tpot > 0) call post_data(CS%id_Tpot, work_3d, CS%diag)
      if (CS%id_tob > 0) call post_data(CS%id_tob, work_3d(:,:,nz), CS%diag, mask=G%mask2dT)
      ! volume mean potential temperature
      if (CS%id_thetaoga>0) then
        thetaoga = global_volume_mean(work_3d, h, G, GV, tmp_scale=US%C_to_degC)
        call post_data(CS%id_thetaoga, thetaoga, CS%diag)
      endif
      ! volume mean conservative temperature
      if (CS%id_bigthetaoga>0) then
        thetaoga = global_volume_mean(tv%T, h, G, GV, tmp_scale=US%C_to_degC)
        call post_data(CS%id_bigthetaoga, thetaoga, CS%diag)
      endif
      ! area mean potential SST
      if (CS%id_tosga > 0) then
        tosga = global_area_mean(work_3d(:,:,1), G, tmp_scale=US%C_to_degC)
        call post_data(CS%id_tosga, tosga, CS%diag)
      endif
      ! area mean conservative SST
      if (CS%id_bigtosga > 0) then
        tosga = global_area_mean(tv%T(:,:,1), G, tmp_scale=US%C_to_degC)
        call post_data(CS%id_bigtosga, tosga, CS%diag)
      endif
      ! layer mean potential temperature
      if (CS%id_temp_layer_ave>0) then
        temp_layer_ave = global_layer_mean(work_3d, h, G, GV, tmp_scale=US%C_to_degC)
        call post_data(CS%id_temp_layer_ave, temp_layer_ave, CS%diag)
      endif
      ! layer mean conservative temperature
      if (CS%id_bigtemp_layer_ave>0) then
        temp_layer_ave = global_layer_mean(tv%T, h, G, GV, tmp_scale=US%C_to_degC)
        call post_data(CS%id_bigtemp_layer_ave, temp_layer_ave, CS%diag)
      endif
      if (CS%id_tosq > 0) then
         do k=1,nz ; do j=js,je ; do i=is,ie
           work_3d(i,j,k) = work_3d(i,j,k)*work_3d(i,j,k)
         enddo ; enddo ; enddo
         call post_data(CS%id_tosq, work_3d, CS%diag)
      endif
    endif
  else
    ! Internal T&S variables are potential temperature & practical salinity
    if (CS%id_tob > 0) call post_data(CS%id_tob, tv%T(:,:,nz), CS%diag, mask=G%mask2dT)
    if (CS%id_tosq > 0) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        work_3d(i,j,k) = tv%T(i,j,k)*tv%T(i,j,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_tosq, work_3d, CS%diag)
    endif
    ! volume mean potential temperature
    if (CS%id_thetaoga>0) then
      thetaoga = global_volume_mean(tv%T, h, G, GV, tmp_scale=US%C_to_degC)
      call post_data(CS%id_thetaoga, thetaoga, CS%diag)
    endif
    ! area mean SST
    if (CS%id_tosga > 0) then
      tosga = global_area_mean(tv%T(:,:,1), G, tmp_scale=US%C_to_degC)
      call post_data(CS%id_tosga, tosga, CS%diag)
    endif
    ! layer mean potential temperature
    if (CS%id_temp_layer_ave>0) then
      temp_layer_ave = global_layer_mean(tv%T, h, G, GV, tmp_scale=US%C_to_degC)
      call post_data(CS%id_temp_layer_ave, temp_layer_ave, CS%diag)
    endif
  endif


  ! Calculate additional, potentially derived salinity diagnostics
  if (tv%S_is_absS) then
    ! Internal T&S variables are conservative temperature & absolute salinity,
    ! so they need to converted to potential temperature and practical salinity
    ! for some diagnostics using TEOS-10 function calls.
    if ((CS%id_Sprac > 0) .or. (CS%id_sob > 0) .or. (CS%id_sosq >0)) then
      EOSdom(:) = EOS_domain(G%HI)
      do k=1,nz ; do j=js,je
        call abs_saln_to_prac_saln(tv%S(:,j,k), work_3d(:,j,k), tv%eqn_of_state, EOSdom)
      enddo ; enddo
      if (CS%id_Sprac > 0) call post_data(CS%id_Sprac, work_3d, CS%diag)
      if (CS%id_sob > 0) call post_data(CS%id_sob, work_3d(:,:,nz), CS%diag, mask=G%mask2dT)
      ! volume mean salinity
      if (CS%id_soga>0) then
        soga = global_volume_mean(work_3d, h, G, GV, tmp_scale=US%S_to_ppt)
        call post_data(CS%id_soga, soga, CS%diag)
      endif
      ! volume mean absolute salinity
      if (CS%id_abssoga>0) then
        soga = global_volume_mean(tv%S, h, G, GV, tmp_scale=US%S_to_ppt)
        call post_data(CS%id_abssoga, soga, CS%diag)
      endif
      ! area mean practical SSS
      if (CS%id_sosga > 0) then
        sosga = global_area_mean(work_3d(:,:,1), G, tmp_scale=US%S_to_ppt)
        call post_data(CS%id_sosga, sosga, CS%diag)
      endif
      ! area mean absolute SSS
      if (CS%id_abssosga > 0) then
        sosga = global_area_mean(tv%S(:,:,1), G, tmp_scale=US%S_to_ppt)
        call post_data(CS%id_abssosga, sosga, CS%diag)
      endif
      ! layer mean practical salinity
      if (CS%id_salt_layer_ave>0) then
        salt_layer_ave = global_layer_mean(work_3d, h, G, GV, tmp_scale=US%S_to_ppt)
        call post_data(CS%id_salt_layer_ave, salt_layer_ave, CS%diag)
      endif
      ! layer mean absolute salinity
      if (CS%id_abssalt_layer_ave>0) then
        salt_layer_ave = global_layer_mean(tv%S, h, G, GV, tmp_scale=US%S_to_ppt)
        call post_data(CS%id_abssalt_layer_ave, salt_layer_ave, CS%diag)
      endif
      if (CS%id_sosq > 0) then
        do k=1,nz ; do j=js,je ; do i=is,ie
           work_3d(i,j,k) = work_3d(i,j,k)*work_3d(i,j,k)
        enddo ; enddo ; enddo
        call post_data(CS%id_sosq, work_3d, CS%diag)
      endif
    endif
  else
    ! Internal T&S variables are potential temperature & practical salinity
    if (CS%id_sob > 0) call post_data(CS%id_sob, tv%S(:,:,nz), CS%diag, mask=G%mask2dT)
    if (CS%id_sosq > 0) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        work_3d(i,j,k) = tv%S(i,j,k)*tv%S(i,j,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_sosq, work_3d, CS%diag)
    endif
    ! volume mean salinity
    if (CS%id_soga>0) then
      soga = global_volume_mean(tv%S, h, G, GV, tmp_scale=US%S_to_ppt)
      call post_data(CS%id_soga, soga, CS%diag)
    endif
    ! area mean SSS
    if (CS%id_sosga > 0) then
      sosga = global_area_mean(tv%S(:,:,1), G, tmp_scale=US%S_to_ppt)
      call post_data(CS%id_sosga, sosga, CS%diag)
    endif
    ! layer mean salinity
    if (CS%id_salt_layer_ave>0) then
      salt_layer_ave = global_layer_mean(tv%S, h, G, GV, tmp_scale=US%S_to_ppt)
      call post_data(CS%id_salt_layer_ave, salt_layer_ave, CS%diag)
    endif
  endif

  call calculate_vertical_integrals(h, tv, p_surf, G, GV, US, CS)

  if ((CS%id_Rml > 0) .or. (CS%id_Rcv > 0) .or. (CS%id_h_Rlay > 0) .or. &
      (CS%id_uh_Rlay > 0) .or. (CS%id_vh_Rlay > 0) .or. &
      (CS%id_uhGM_Rlay > 0) .or. (CS%id_vhGM_Rlay > 0)) then

    if (associated(tv%eqn_of_state)) then
      EOSdom(:) = EOS_domain(G%HI, halo=1)
      pressure_1d(:) = tv%P_Ref
      !$OMP parallel do default(shared)
      do k=1,nz ; do j=js-1,je+1
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k),  pressure_1d, Rcv(:,j,k), tv%eqn_of_state, &
                               EOSdom)
      enddo ; enddo
    else ! Rcv should not be used much in this case, so fill in sensible values.
      do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
        Rcv(i,j,k) = GV%Rlay(k)
      enddo ; enddo ; enddo
    endif
    if (CS%id_Rml > 0) call post_data(CS%id_Rml, Rcv, CS%diag)
    if (CS%id_Rcv > 0) call post_data(CS%id_Rcv, Rcv, CS%diag)

    if (CS%id_h_Rlay > 0) then
      ! Here work_3d is used for the layer thicknesses in potential density coordinates [H ~> m or kg m-2].
      k_list = nz/2
      !$OMP parallel do default(shared) private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do i=is,ie
          work_3d(i,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          work_3d(i,j,k) = h(i,j,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, Rcv(i,j,k), k_list, nz, wt, wt_p)
          work_3d(i,j,k_list)   = work_3d(i,j,k_list)   + h(i,j,k)*wt
          work_3d(i,j,k_list+1) = work_3d(i,j,k_list+1) + h(i,j,k)*wt_p
        enddo ; enddo
      enddo

      call post_data(CS%id_h_Rlay, work_3d, CS%diag)
    endif

    if (CS%id_uh_Rlay > 0) then
      ! Calculate zonal transports in potential density coordinates [H L2 T-1 ~> m3 s-1 or kg s-1].
      k_list = nz/2
      !$OMP parallel do default(shared) private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do I=Isq,Ieq
          uh_tmp(I,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do I=Isq,Ieq
          uh_tmp(I,j,k) = uh(I,j,k)
        enddo ; enddo
        k_list = nz/2
        do k=1,nkmb ; do I=Isq,Ieq
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i+1,j,k)), k_list, nz, wt, wt_p)
          uh_tmp(I,j,k_list)   = uh_tmp(I,j,k_list)   + uh(I,j,k)*wt
          uh_tmp(I,j,k_list+1) = uh_tmp(I,j,k_list+1) + uh(I,j,k)*wt_p
        enddo ; enddo
      enddo

      call post_data(CS%id_uh_Rlay, uh_tmp, CS%diag)
    endif

    if (CS%id_vh_Rlay > 0) then
      ! Calculate meridional transports in potential density coordinates [H L2 T-1 ~> m3 s-1 or kg s-1].
      k_list = nz/2
      !$OMP parallel do default(shared) private(wt,wt_p) firstprivate(k_list)
      do J=Jsq,Jeq
        do k=1,nkmb ; do i=is,ie
          vh_tmp(i,J,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          vh_tmp(i,J,k) = vh(i,J,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i,j+1,k)), k_list, nz, wt, wt_p)
          vh_tmp(i,J,k_list)   = vh_tmp(i,J,k_list)   + vh(i,J,k)*wt
          vh_tmp(i,J,k_list+1) = vh_tmp(i,J,k_list+1) + vh(i,J,k)*wt_p
        enddo ; enddo
      enddo

      call post_data(CS%id_vh_Rlay, vh_tmp, CS%diag)
    endif

    if ((CS%id_uhGM_Rlay > 0) .and. associated(CDp%uhGM)) then
      ! Calculate zonal Gent-McWilliams transports in potential density
      ! coordinates [H L2 T-1 ~> m3 s-1 or kg s-1].
      k_list = nz/2
      !$OMP parallel do default(shared) private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do I=Isq,Ieq
          uh_tmp(I,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do I=Isq,Ieq
          uh_tmp(I,j,k) = CDp%uhGM(I,j,k)
        enddo ; enddo
        do k=1,nkmb ; do I=Isq,Ieq
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i+1,j,k)), k_list, nz, wt, wt_p)
          uh_tmp(I,j,k_list)   = uh_tmp(I,j,k_list)   + CDp%uhGM(I,j,k)*wt
          uh_tmp(I,j,k_list+1) = uh_tmp(I,j,k_list+1) + CDp%uhGM(I,j,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_uhGM_Rlay > 0) call post_data(CS%id_uhGM_Rlay, uh_tmp, CS%diag)
    endif

    if ((CS%id_vhGM_Rlay > 0) .and. associated(CDp%vhGM)) then
      ! Calculate meridional Gent-McWilliams transports in potential density
      ! coordinates [H L2 T-1 ~> m3 s-1 or kg s-1].
      k_list = nz/2
      !$OMP parallel do default(shared) private(wt,wt_p) firstprivate(k_list)
      do J=Jsq,Jeq
        do k=1,nkmb ; do i=is,ie
          vh_tmp(i,J,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          vh_tmp(i,J,k) = CDp%vhGM(i,J,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i,j+1,k)), k_list, nz, wt, wt_p)
          vh_tmp(i,J,k_list)   = vh_tmp(i,J,k_list)   + CDp%vhGM(i,J,k)*wt
          vh_tmp(i,J,k_list+1) = vh_tmp(i,J,k_list+1) + CDp%vhGM(i,J,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_vhGM_Rlay > 0) call post_data(CS%id_vhGM_Rlay, vh_tmp, CS%diag)
    endif
  endif

  if (associated(tv%eqn_of_state)) then
    EOSdom(:) = EOS_domain(G%HI)
    if (CS%id_rhopot0 > 0) then
      pressure_1d(:) = 0.
      !$OMP parallel do default(shared)
      do k=1,nz ; do j=js,je
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k),  pressure_1d, Rcv(:,j,k), &
                                tv%eqn_of_state, EOSdom)
      enddo ; enddo
      if (CS%id_rhopot0 > 0) call post_data(CS%id_rhopot0, Rcv, CS%diag)
    endif
    if (CS%id_rhopot2 > 0) then
      pressure_1d(:) = 2.0e7*US%Pa_to_RL2_T2 ! 2000 dbars
      !$OMP parallel do default(shared)
      do k=1,nz ; do j=js,je
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k),  pressure_1d, Rcv(:,j,k), &
                                tv%eqn_of_state, EOSdom)
      enddo ; enddo
      if (CS%id_rhopot2 > 0) call post_data(CS%id_rhopot2, Rcv, CS%diag)
    endif
    if (CS%id_rhoinsitu > 0) then
      !$OMP parallel do default(shared) private(pressure_1d)
      do j=js,je
        pressure_1d(:) = 0. ! Start at p=0 Pa at surface
        do k=1,nz
          pressure_1d(:) =  pressure_1d(:) + 0.5 * h(:,j,k) * (GV%H_to_RZ*GV%g_Earth) ! Pressure in middle of layer k
          call calculate_density(tv%T(:,j,k), tv%S(:,j,k),  pressure_1d, Rcv(:,j,k), &
                                 tv%eqn_of_state, EOSdom)
          pressure_1d(:) =  pressure_1d(:) + 0.5 * h(:,j,k) * (GV%H_to_RZ*GV%g_Earth) ! Pressure at bottom of layer k
        enddo
      enddo
      if (CS%id_rhoinsitu > 0) call post_data(CS%id_rhoinsitu, Rcv, CS%diag)
    endif

    if (CS%id_drho_dT > 0 .or. CS%id_drho_dS > 0) then
      !$OMP parallel do default(shared) private(pressure_1d)
      do j=js,je
        pressure_1d(:) = 0. ! Start at p=0 Pa at surface
        do k=1,nz
          pressure_1d(:) =  pressure_1d(:) + 0.5 * h(:,j,k) * (GV%H_to_RZ*GV%g_Earth) ! Pressure in middle of layer k
          ! To avoid storing more arrays, put drho_dT into Rcv, and drho_dS into work3d
          call calculate_density_derivs(tv%T(:,j,k), tv%S(:,j,k), pressure_1d, &
                                        Rcv(:,j,k), work_3d(:,j,k), tv%eqn_of_state, EOSdom)
          pressure_1d(:) =  pressure_1d(:) + 0.5 * h(:,j,k) * (GV%H_to_RZ*GV%g_Earth) ! Pressure at bottom of layer k
        enddo
      enddo
      if (CS%id_drho_dT > 0) call post_data(CS%id_drho_dT, Rcv, CS%diag)
      if (CS%id_drho_dS > 0) call post_data(CS%id_drho_dS, work_3d, CS%diag)
    endif
  endif

  if ((CS%id_cg1>0) .or. (CS%id_Rd1>0) .or. (CS%id_cfl_cg1>0) .or. &
      (CS%id_cfl_cg1_x>0) .or. (CS%id_cfl_cg1_y>0)) then
    call wave_speed(h, tv, G, GV, US, cg1, CS%wave_speed)
    if (CS%id_cg1>0) call post_data(CS%id_cg1, cg1, CS%diag)
    if (CS%id_Rd1>0) then
      !$OMP parallel do default(shared) private(f2_h,mag_beta)
      do j=js,je ; do i=is,ie
        ! Blend the equatorial deformation radius with the standard one.
        f2_h = absurdly_small_freq2 + 0.25 * &
            ((G%Coriolis2Bu(I,J) + G%Coriolis2Bu(I-1,J-1)) + &
             (G%Coriolis2Bu(I-1,J) + G%Coriolis2Bu(I,J-1)))
        mag_beta = sqrt(0.5 * ( &
            ((((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2) + &
             (((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2)) + &
            ((((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2) + &
             (((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2)) ))
        Rd1(i,j) = cg1(i,j) / sqrt(f2_h + cg1(i,j) * mag_beta)

      enddo ; enddo
      call post_data(CS%id_Rd1, Rd1, CS%diag)
    endif
    if (CS%id_cfl_cg1>0) then
      do j=js,je ; do i=is,ie
        CFL_cg1(i,j) = (dt*cg1(i,j)) * (G%IdxT(i,j) + G%IdyT(i,j))
      enddo ; enddo
      call post_data(CS%id_cfl_cg1, CFL_cg1, CS%diag)
    endif
    if (CS%id_cfl_cg1_x>0) then
      do j=js,je ; do i=is,ie
        CFL_cg1(i,j) = (dt*cg1(i,j)) * G%IdxT(i,j)
      enddo ; enddo
      call post_data(CS%id_cfl_cg1_x, CFL_cg1, CS%diag)
    endif
    if (CS%id_cfl_cg1_y>0) then
      do j=js,je ; do i=is,ie
        CFL_cg1(i,j) = (dt*cg1(i,j)) * G%IdyT(i,j)
      enddo ; enddo
      call post_data(CS%id_cfl_cg1_y, CFL_cg1, CS%diag)
    endif
  endif
  if ((CS%id_cg_ebt>0) .or. (CS%id_Rd_ebt>0) .or. (CS%id_p_ebt>0)) then
    if (CS%id_p_ebt>0) then
      ! Here work_3d is used for the equivalent barotropic modal structure [nondim].
      work_3d(:,:,:) = 0.0
      call wave_speed(h, tv, G, GV, US, cg1, CS%wave_speed, use_ebt_mode=.true., &
                      mono_N2_column_fraction=CS%mono_N2_column_fraction, &
                      mono_N2_depth=CS%mono_N2_depth, modal_structure=work_3d)
      call post_data(CS%id_p_ebt, work_3d, CS%diag)
    else
      call wave_speed(h, tv, G, GV, US, cg1, CS%wave_speed, use_ebt_mode=.true., &
                      mono_N2_column_fraction=CS%mono_N2_column_fraction, &
                      mono_N2_depth=CS%mono_N2_depth)
    endif
    if (CS%id_cg_ebt>0) call post_data(CS%id_cg_ebt, cg1, CS%diag)
    if (CS%id_Rd_ebt>0) then
      !$OMP parallel do default(shared) private(f2_h,mag_beta)
      do j=js,je ; do i=is,ie
        ! Blend the equatorial deformation radius with the standard one.
        f2_h = absurdly_small_freq2 + 0.25 * &
            ((G%Coriolis2Bu(I,J) + G%Coriolis2Bu(I-1,J-1)) + &
             (G%Coriolis2Bu(I-1,J) + G%Coriolis2Bu(I,J-1)))
        mag_beta = sqrt(0.5 * ( &
            ((((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2) + &
             (((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2)) + &
            ((((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2) + &
             (((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2)) ))
        Rd1(i,j) = cg1(i,j) / sqrt(f2_h + cg1(i,j) * mag_beta)

      enddo ; enddo
      call post_data(CS%id_Rd_ebt, Rd1, CS%diag)
    endif
  endif

end subroutine calculate_diagnostic_fields

!> This subroutine finds the location of R_in in an increasing ordered
!! list, Rlist, returning as k the element such that
!! Rlist(k) <= R_in < Rlist(k+1), and where wt and wt_p are the linear
!! weights that should be assigned to elements k and k+1.
subroutine find_weights(Rlist, R_in, k, nz, wt, wt_p)
  real, dimension(:), &
            intent(in)    :: Rlist !< The list of target densities [R ~> kg m-3]
  real,     intent(in)    :: R_in !< The density being inserted into Rlist [R ~> kg m-3]
  integer,  intent(inout) :: k    !< The value of k such that Rlist(k) <= R_in < Rlist(k+1)
                                  !! The input value is a first guess
  integer,  intent(in)    :: nz   !< The number of layers in Rlist
  real,     intent(out)   :: wt   !< The weight of layer k for interpolation [nondim]
  real,     intent(out)   :: wt_p !< The weight of layer k+1 for interpolation [nondim]

  ! This subroutine finds location of R_in in an increasing ordered
  ! list, Rlist, returning as k the element such that
  !  Rlist(k) <= R_in < Rlist(k+1), and where wt and wt_p are the linear
  ! weights that should be assigned to elements k and k+1.

  integer :: k_upper, k_lower, k_new, inc

  ! First, bracket the desired point.
  if ((k < 1) .or. (k > nz)) k = nz/2

  k_upper = k ; k_lower = k ; inc = 1
  if (R_in < Rlist(k)) then
    do
      k_lower = max(k_lower-inc, 1)
      if ((k_lower == 1) .or. (R_in >= Rlist(k_lower))) exit
      k_upper = k_lower
      inc = inc*2
    enddo
  else
    do
      k_upper = min(k_upper+inc, nz)
      if ((k_upper == nz) .or. (R_in < Rlist(k_upper))) exit
      k_lower = k_upper
      inc = inc*2
    enddo
  endif

  if ((k_lower == 1) .and. (R_in <= Rlist(k_lower))) then
    k = 1 ; wt = 1.0 ; wt_p = 0.0
  elseif ((k_upper == nz) .and. (R_in >= Rlist(k_upper))) then
    k = nz-1 ; wt = 0.0 ; wt_p = 1.0
  else
    do
      if (k_upper <= k_lower+1) exit
      k_new = (k_upper + k_lower) / 2
      if (R_in < Rlist(k_new)) then
        k_upper = k_new
      else
        k_lower = k_new
      endif
    enddo

!   Uncomment this as a code check
!    if ((R_in < Rlist(k_lower)) .or. (R_in >= Rlist(k_upper)) .or. (k_upper-k_lower /= 1)) &
!      write (*,*) "Error: ",R_in," is not between R(",k_lower,") = ", &
!        Rlist(k_lower)," and R(",k_upper,") = ",Rlist(k_upper),"."
    k = k_lower
    wt = (Rlist(k_upper) - R_in) / (Rlist(k_upper) - Rlist(k_lower))
    wt_p = 1.0 - wt

  endif

end subroutine find_weights

!> This subroutine calculates vertical integrals of several tracers, along
!! with the mass-weight of these tracers, the total column mass, and the
!! carefully calculated column height.
subroutine calculate_vertical_integrals(h, tv, p_surf, G, GV, US, CS)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various
                                                 !! thermodynamic variables.
  real, dimension(:,:),    pointer       :: p_surf !< A pointer to the surface pressure [R L2 T-2 ~> Pa].
                                                 !! If p_surf is not associated, it is the same
                                                 !! as setting the surface pressure to 0.
  type(diagnostics_CS),    intent(inout) :: CS   !< Control structure returned by a
                                                 !! previous call to diagnostics_init.

  real, dimension(SZI_(G),SZJ_(G)) :: &
    z_top, &  ! Height of the top of a layer or the ocean [Z ~> m].
    z_bot, &  ! Height of the bottom of a layer (for id_mass) or the
              ! (positive) depth of the ocean (for id_col_ht) [Z ~> m].
    mass, &   ! integrated mass of the water column [R Z ~> kg m-2].  For
              ! non-Boussinesq models this is rho*dz. For Boussinesq
              ! models, this is either the integral of in-situ density
              ! (rho*dz for col_mass) or reference density (Rho_0*dz for mass_wt).
    btm_pres,&! The pressure at the ocean bottom, or CMIP variable 'pbo'.
              ! This is the column mass multiplied by gravity plus the pressure
              ! at the ocean surface [R L2 T-2 ~> Pa].
    dpress, & ! Change in hydrostatic pressure across a layer [R L2 T-2 ~> Pa].
    tr_int    ! vertical integral of a tracer times density,
              ! (Rho_0 in a Boussinesq model) [Conc R Z ~> Conc kg m-2].
  real    :: IG_Earth  ! Inverse of gravitational acceleration [T2 Z L-2 ~> s2 m-1].

  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (CS%id_mass_wt > 0) then
    do j=js,je ; do i=is,ie ; mass(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      mass(i,j) = mass(i,j) + GV%H_to_RZ*h(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_mass_wt, mass, CS%diag)
  endif

  if (CS%id_temp_int > 0) then
    do j=js,je ; do i=is,ie ; tr_int(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      tr_int(i,j) = tr_int(i,j) + (GV%H_to_RZ*h(i,j,k))*tv%T(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_temp_int, tr_int, CS%diag)
  endif

  if (CS%id_salt_int > 0) then
    do j=js,je ; do i=is,ie ; tr_int(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      tr_int(i,j) = tr_int(i,j) + (GV%H_to_RZ*h(i,j,k))*tv%S(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_salt_int, tr_int, CS%diag)
  endif

  if (CS%id_col_ht > 0) then
    call find_eta(h, tv, G, GV, US, z_top)
    do j=js,je ; do i=is,ie
      z_bot(i,j) = z_top(i,j) + G%bathyT(i,j)
    enddo ; enddo
    call post_data(CS%id_col_ht, z_bot, CS%diag)
  endif

  ! NOTE: int_density_z expects z_top and z_btm values from [ij]sq to [ij]eq+1
  if (CS%id_col_mass > 0 .or. CS%id_pbo > 0) then
    do j=js,je ; do i=is,ie ; mass(i,j) = 0.0 ; enddo ; enddo
    if (GV%Boussinesq) then
      if (associated(tv%eqn_of_state)) then
        IG_Earth = 1.0 / GV%g_Earth
        do j=G%jscB,G%jecB+1 ; do i=G%iscB,G%iecB+1
          z_bot(i,j) = 0.0
        enddo ; enddo
        do k=1,nz
          do j=G%jscB,G%jecB+1 ; do i=G%iscB,G%iecB+1
            z_top(i,j) = z_bot(i,j)
            z_bot(i,j) = z_top(i,j) - GV%H_to_Z*h(i,j,k)
          enddo ; enddo
          call int_density_dz(tv%T(:,:,k), tv%S(:,:,k), z_top, z_bot, 0.0, GV%Rho0, GV%g_Earth, &
                              G%HI, tv%eqn_of_state, US, dpress)
          do j=js,je ; do i=is,ie
            mass(i,j) = mass(i,j) + dpress(i,j) * IG_Earth
          enddo ; enddo
        enddo
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          mass(i,j) = mass(i,j) + (GV%H_to_Z*GV%Rlay(k))*h(i,j,k)
        enddo ; enddo ; enddo
      endif
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        mass(i,j) = mass(i,j) + GV%H_to_RZ*h(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (CS%id_col_mass > 0) then
      call post_data(CS%id_col_mass, mass, CS%diag)
    endif
    if (CS%id_pbo > 0) then
      do j=js,je ; do i=is,ie ; btm_pres(i,j) = 0.0 ; enddo ; enddo
      ! 'pbo' is defined as the sea water pressure at the sea floor
      !     pbo = (mass * g) + p_surf
      ! where p_surf is the sea water pressure at sea water surface.
      do j=js,je ; do i=is,ie
        btm_pres(i,j) = GV%g_Earth * mass(i,j)
        if (associated(p_surf)) then
          btm_pres(i,j) = btm_pres(i,j) + p_surf(i,j)
        endif
      enddo ; enddo
      call post_data(CS%id_pbo, btm_pres, CS%diag)
    endif
  endif

end subroutine calculate_vertical_integrals

!> This subroutine calculates terms in the mechanical energy budget.
subroutine calculate_energy_diagnostics(u, v, h, uh, vh, ADp, CDp, G, GV, US, CS)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: uh   !< Transport through zonal faces=u*h*dy,
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: vh   !< Transport through merid faces=v*h*dx,
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(accel_diag_ptrs),   intent(in)    :: ADp  !< Structure pointing to accelerations in momentum equation.
  type(cont_diag_ptrs),    intent(in)    :: CDp  !< Structure pointing to terms in continuity equations.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(diagnostics_CS),    intent(inout) :: CS   !< Control structure returned by a previous call to
                                                 !! diagnostics_init.

  ! Local variables
  real :: KE(SZI_(G),SZJ_(G),SZK_(GV)) ! Kinetic energy per unit mass [L2 T-2 ~> m2 s-2]
  real :: KE_term(SZI_(G),SZJ_(G),SZK_(GV)) ! A term in the kinetic energy budget
                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]
  real :: KE_u(SZIB_(G),SZJ_(G)) ! The area integral of a KE term in a layer at u-points
                                 ! [H L4 T-3 ~> m5 s-3 or W]
  real :: KE_v(SZI_(G),SZJB_(G)) ! The area integral of a KE term in a layer at v-points
                                 ! [H L4 T-3 ~> m5 s-3 or W]
  real :: KE_h(SZI_(G),SZJ_(G))  ! A KE term contribution at tracer points
                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.(CS%KE_term_on .or. (CS%id_KE > 0))) return

  KE_u(:,:) = 0. ; KE_v(:,:) = 0.

  do k=1,nz ; do j=js,je ; do i=is,ie
    KE(i,j,k) = (((u(I,j,k) * u(I,j,k)) + (u(I-1,j,k) * u(I-1,j,k))) &
               + ((v(i,J,k) * v(i,J,k)) + (v(i,J-1,k) * v(i,J-1,k)))) * 0.25
  enddo ; enddo ; enddo
  if (CS%id_KE > 0) call post_data(CS%id_KE, KE, CS%diag)

  if (CS%KE_term_on .and. .not.G%symmetric) then
    call create_group_pass(CS%pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)
  endif

  if (CS%id_dKEdt > 0) then
    ! Calculate the time derivative of the layer KE [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * CS%du_dt(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * CS%dv_dt(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = KE(i,j,k) * CS%dh_dt(i,j,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_dKEdt, KE_term, CS%diag)
  endif

  if (CS%id_PE_to_KE > 0) then
    ! Calculate the potential energy to KE term [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%PFu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%PFv(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_PE_to_KE, KE_term, CS%diag)
  endif

  if (CS%id_KE_SAL > 0) then
    ! Calculate the KE source from self-attraction and loading [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%sal_u(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%sal_v(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_SAL, KE_term, CS%diag)
  endif

  if (CS%id_KE_TIDES > 0) then
    ! Calculate the KE source from astronomical tidal forcing [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%tides_u(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%tides_v(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_TIDES, KE_term, CS%diag)
  endif

  if (CS%id_KE_BT > 0) then
    ! Calculate the barotropic contribution to KE term [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%u_accel_bt(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%v_accel_bt(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_BT, KE_term, CS%diag)
  endif

  if (CS%id_PE_to_KE_btbc > 0) then
    ! Calculate the potential energy to KE term including barotropic solver contribution
    ! [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * (ADp%PFu(I,j,k) + ADp%bt_pgf_u(I,j,k))
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * (ADp%PFv(i,J,k) + ADp%bt_pgf_v(i,J,k))
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_PE_to_KE_btbc, KE_term, CS%diag)
  endif

  if (CS%id_KE_Coradv_btbc > 0) then
    ! Calculate the KE source from Coriolis and advection terms including barotropic solver contribution
    ! [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * (ADp%CAu(I,j,k) + ADp%bt_cor_u(I,j))
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * (ADp%CAv(i,J,k) + ADp%bt_cor_v(i,J))
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = -KE(i,j,k) * G%IareaT(i,j) &
            * ((uh(I,j,k) - uh(I-1,j,k)) + (vh(i,J,k) - vh(i,J-1,k)))
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_Coradv_btbc, KE_term, CS%diag)
  endif

  if (CS%id_KE_BT_PF > 0) then
    ! Calculate the anomalous pressure gradient force contribution to KE term [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%bt_pgf_u(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%bt_pgf_v(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_BT_PF, KE_term, CS%diag)
  endif

  if (CS%id_KE_BT_CF > 0) then
    ! Calculate the anomalous Coriolis force contribution to KE term [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%bt_cor_u(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%bt_cor_v(i,J)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_BT_CF, KE_term, CS%diag)
  endif

  if (CS%id_KE_BT_WD > 0) then
    ! Calculate the barotropic linear wave drag contribution to KE term [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%bt_lwd_u(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%bt_lwd_v(i,J)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_BT_WD, KE_term, CS%diag)
  endif

  if (CS%id_KE_Coradv > 0) then
    ! Calculate the KE source from the combined Coriolis and advection terms [H L2 T-3 ~> m3 s-3 or W m-2].
    ! The Coriolis source should be zero, but is not due to truncation errors.  There should be
    ! near-cancellation of the global integral of this spurious Coriolis source.
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%CAu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%CAv(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = -KE(i,j,k) * G%IareaT(i,j) &
            * ((uh(I,j,k) - uh(I-1,j,k)) + (vh(i,J,k) - vh(i,J-1,k)))
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_Coradv, KE_term, CS%diag)
  endif

  if (CS%id_KE_adv > 0) then
    ! Calculate the KE source from along-layer advection [H L2 T-3 ~> m3 s-3 or W m-2].
    ! NOTE: All terms in KE_adv are multiplied by -1, which can easily produce
    ! negative zeros and may signal a reproducibility issue over land.
    ! We resolve this by re-initializing and only evaluating over water points.
    KE_u(:,:) = 0. ; KE_v(:,:) = 0.
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        if (G%mask2dCu(i,j) /= 0.) &
          KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%gradKEu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        if (G%mask2dCv(i,j) /= 0.) &
          KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%gradKEv(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = -KE(i,j,k) * G%IareaT(i,j) &
            * ((uh(I,j,k) - uh(I-1,j,k)) + (vh(i,J,k) - vh(i,J-1,k)))
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_adv, KE_term, CS%diag)
  endif

  if (CS%id_KE_visc > 0) then
    ! Calculate the KE source from vertical viscosity [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%du_dt_visc(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%dv_dt_visc(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_visc, KE_term, CS%diag)
  endif

  if (CS%id_KE_visc_gl90 > 0) then
    ! Calculate the KE source from GL90 vertical viscosity [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%du_dt_visc_gl90(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%dv_dt_visc_gl90(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_visc_gl90, KE_term, CS%diag)
  endif

  if (CS%id_KE_stress > 0) then
    ! Calculate the KE source from surface stress (included in KE_visc) [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%du_dt_str(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%dv_dt_str(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) * &
            ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_stress, KE_term, CS%diag)
  endif

  if (CS%id_KE_horvisc > 0) then
    ! Calculate the KE source from horizontal viscosity [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%diffu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%diffv(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_horvisc, KE_term, CS%diag)
  endif

  if (CS%id_KE_dia > 0) then
    ! Calculate the KE source from diapycnal diffusion [H L2 T-3 ~> m3 s-3 or W m-2].
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%du_dt_dia(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%dv_dt_dia(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = KE(i,j,k) * (CDp%diapyc_vel(i,j,k) - CDp%diapyc_vel(i,j,k+1))
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) &
            * ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    call post_data(CS%id_KE_dia, KE_term, CS%diag)
  endif

end subroutine calculate_energy_diagnostics

!> This subroutine registers fields to calculate a diagnostic time derivative.
subroutine register_time_deriv(lb, f_ptr, deriv_ptr, CS)
  integer, intent(in), dimension(3) :: lb     !< Lower index bound of f_ptr
  real, dimension(lb(1):,lb(2):,:), target :: f_ptr
                                              !< Time derivative operand, in arbitrary units [A ~> a]
  real, dimension(lb(1):,lb(2):,:), target :: deriv_ptr
                                              !< Time derivative of f_ptr, in units derived from
                                              !! the arbitrary units of f_ptr [A T-1 ~> a s-1]
  type(diagnostics_CS), intent(inout) :: CS   !< Control structure returned by previous call to
                                              !! diagnostics_init.

  ! This subroutine registers fields to calculate a diagnostic time derivative.
  ! NOTE: Lower bound is required for grid indexing in calculate_derivs().
  !       We assume that the vertical axis is 1-indexed.

  integer :: m      !< New index of deriv_ptr in CS%deriv
  integer :: ub(3)  !< Upper index bound of f_ptr, based on shape.

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "register_time_deriv: Module must be initialized before it is used.")

  if (CS%num_time_deriv >= MAX_FIELDS_) then
    call MOM_error(WARNING,"MOM_diagnostics:  Attempted to register more than " // &
                   "MAX_FIELDS_ diagnostic time derivatives via register_time_deriv.")
    return
  endif

  m = CS%num_time_deriv+1 ; CS%num_time_deriv = m

  ub(:) = lb(:) + shape(f_ptr) - 1

  CS%nlay(m) = size(f_ptr, 3)
  CS%deriv(m)%p => deriv_ptr
  allocate(CS%prev_val(m)%p(lb(1):ub(1), lb(2):ub(2), CS%nlay(m)))

  CS%var_ptr(m)%p => f_ptr
  CS%prev_val(m)%p(:,:,:) = f_ptr(:,:,:)

end subroutine register_time_deriv

!> This subroutine calculates all registered time derivatives.
subroutine calculate_derivs(dt, G, CS)
  real,                  intent(in)    :: dt   !< The time interval over which differences occur [T ~> s].
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  type(diagnostics_CS),  intent(inout) :: CS   !< Control structure returned by previous call to
                                               !! diagnostics_init.

! This subroutine calculates all registered time derivatives.
  real :: Idt  ! The inverse timestep [T-1 ~> s-1]
  integer :: i, j, k, m

  if (dt > 0.0) then ; Idt = 1.0/dt
  else ; return ; endif

  ! Because the field is unknown, its grid index bounds are also unknown.
  ! Additionally, two of the fields (dudt, dvdt) require calculation of spatial
  ! derivatives when computing d(KE)/dt.  This raises issues in non-symmetric
  ! mode, where the symmetric boundaries (west, south) may not be updated.

  ! For this reason, we explicitly loop from isc-1:iec and jsc-1:jec, in order
  ! to force boundary value updates, even though it may not be strictly valid
  ! for all fields.  Note this assumes a halo, and that it has been updated.

  do m=1,CS%num_time_deriv
    do k=1,CS%nlay(m) ; do j=G%jsc-1,G%jec ; do i=G%isc-1,G%iec
      CS%deriv(m)%p(i,j,k) = (CS%var_ptr(m)%p(i,j,k) - CS%prev_val(m)%p(i,j,k)) * Idt
      CS%prev_val(m)%p(i,j,k) = CS%var_ptr(m)%p(i,j,k)
    enddo ; enddo ; enddo
  enddo

end subroutine calculate_derivs

!> This routine posts diagnostics of various dynamic ocean surface quantities,
!! including velocities, speed and sea surface height, at the time the ocean
!! state is reported back to the caller
subroutine post_surface_dyn_diags(IDs, G, diag, sfc_state, ssh)
  type(surface_diag_IDs),   intent(in) :: IDs !< A structure with the diagnostic IDs.
  type(ocean_grid_type),    intent(in) :: G   !< ocean grid structure
  type(diag_ctrl),          intent(in) :: diag !< regulates diagnostic output
  type(surface),            intent(in) :: sfc_state !< structure describing the ocean surface state
  real, dimension(SZI_(G),SZJ_(G)), &
                            intent(in) :: ssh !< Time mean surface height without corrections
                                              !! for ice displacement [Z ~> m]

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: speed  ! The surface speed [L T-1 ~> m s-1]
  real :: ssu_east(SZI_(G),SZJ_(G))        ! Surface velocity due east component [L T-1 ~> m s-1]
  real :: ssv_north(SZI_(G),SZJ_(G))       ! Surface velocity due north component [L T-1 ~> m s-1]
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (IDs%id_ssh > 0) &
    call post_data(IDs%id_ssh, ssh, diag, mask=G%mask2dT)

  if (IDs%id_ssu > 0) &
    call post_data(IDs%id_ssu, sfc_state%u, diag, mask=G%mask2dCu)

  if (IDs%id_ssv > 0) &
    call post_data(IDs%id_ssv, sfc_state%v, diag, mask=G%mask2dCv)

  if (IDs%id_speed > 0) then
    do j=js,je ; do i=is,ie
      speed(i,j) = sqrt(0.5*((sfc_state%u(I-1,j)**2) + (sfc_state%u(I,j)**2)) + &
                        0.5*((sfc_state%v(i,J-1)**2) + (sfc_state%v(i,J)**2)))
    enddo ; enddo
    call post_data(IDs%id_speed, speed, diag, mask=G%mask2dT)
  endif

  if (IDs%id_ssu_east > 0 .or. IDs%id_ssv_north > 0) then
    do j=js,je ; do i=is,ie
      ssu_east(i,j) = ((0.5*(sfc_state%u(I-1,j) + sfc_state%u(I,j))) * G%cos_rot(i,j)) + &
                      ((0.5*(sfc_state%v(i,J-1) + sfc_state%v(i,J))) * G%sin_rot(i,j))
      ssv_north(i,j) = ((0.5*(sfc_state%v(i,J-1) + sfc_state%v(i,J))) * G%cos_rot(i,j)) - &
                       ((0.5*(sfc_state%u(I-1,j) + sfc_state%u(I,j))) * G%sin_rot(i,j))
    enddo ; enddo
    if (IDs%id_ssu_east > 0 ) call post_data(IDs%id_ssu_east, ssu_east, diag, mask=G%mask2dT)
    if (IDs%id_ssv_north > 0 ) call post_data(IDs%id_ssv_north, ssv_north, diag, mask=G%mask2dT)
  endif

end subroutine post_surface_dyn_diags


!> This routine posts diagnostics of various ocean surface and integrated
!! quantities at the time the ocean state is reported back to the caller
subroutine post_surface_thermo_diags(IDs, G, GV, US, diag, dt_int, sfc_state, tv, &
                                    ssh, ssh_ibc)
  type(surface_diag_IDs),   intent(in) :: IDs !< A structure with the diagnostic IDs.
  type(ocean_grid_type),    intent(in) :: G   !< ocean grid structure
  type(verticalGrid_type),  intent(in) :: GV  !< ocean vertical grid structure
  type(unit_scale_type),    intent(in) :: US  !< A dimensional unit scaling type
  type(diag_ctrl),          intent(in) :: diag  !< regulates diagnostic output
  real,                     intent(in) :: dt_int !< total time step associated with these diagnostics [T ~> s].
  type(surface),            intent(in) :: sfc_state !< structure describing the ocean surface state
  type(thermo_var_ptrs),    intent(in) :: tv  !< A structure pointing to various thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: ssh !< Time mean surface height without corrections
                                              !! for ice displacement [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: ssh_ibc !< Time mean surface height with corrections
                                              !! for ice displacement and the inverse barometer [Z ~> m]

  real, dimension(SZI_(G),SZJ_(G)) :: work_2d  ! A 2-d work array [various]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    zos  ! dynamic sea lev (zero area mean) from inverse-barometer adjusted ssh [Z ~> m]
  real :: I_time_int    ! The inverse of the time interval [T-1 ~> s-1].
  real :: zos_area_mean ! Global area mean sea surface height [Z ~> m]
  real :: volo          ! Total volume of the ocean [Z L2 ~> m3]
  real :: ssh_ga        ! Global ocean area weighted mean sea seaface height [Z ~> m]
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! area mean SSH
  if (IDs%id_ssh_ga > 0) then
    ssh_ga = global_area_mean(ssh, G, tmp_scale=US%Z_to_m)
    call post_data(IDs%id_ssh_ga, ssh_ga, diag)
  endif

  ! post the dynamic sea level, zos, and zossq.
  ! zos is ave_ssh with sea ice inverse barometer removed, and with zero global area mean.
  if (IDs%id_zos > 0 .or. IDs%id_zossq > 0) then
    zos_area_mean = global_area_mean(ssh_ibc, G, tmp_scale=US%Z_to_m)
    do j=js,je ; do i=is,ie
      zos(i,j) = ssh_ibc(i,j) - G%mask2dT(i,j)*zos_area_mean
    enddo ; enddo
    if (IDs%id_zos > 0) call post_data(IDs%id_zos, zos, diag, mask=G%mask2dT)
    if (IDs%id_zossq > 0) then
      do j=js,je ; do i=is,ie
        work_2d(i,j) = zos(i,j)*zos(i,j)
      enddo ; enddo
      call post_data(IDs%id_zossq, work_2d, diag, mask=G%mask2dT)
    endif
  endif

  ! post total volume of the liquid ocean
  if (IDs%id_volo > 0) then
    do j=js,je ; do i=is,ie
      work_2d(i,j) = G%mask2dT(i,j) * (ssh(i,j) + G%bathyT(i,j))
    enddo ; enddo
    volo = global_area_integral(work_2d, G, tmp_scale=US%Z_to_m)
    call post_data(IDs%id_volo, volo, diag)
  endif

  ! Use Adcroft's rule of reciprocals; it does the right thing here.
  I_time_int = 0.0 ; if (dt_int > 0.0) I_time_int = 1.0 / dt_int

  ! post time-averaged rate of frazil formation
  if (associated(tv%frazil) .and. (IDs%id_fraz > 0)) then
    do j=js,je ; do i=is,ie
      work_2d(i,j) = tv%frazil(i,j) * I_time_int
    enddo ; enddo
    call post_data(IDs%id_fraz, work_2d, diag, mask=G%mask2dT)
  endif

  ! post time-averaged salt deficit
  if (associated(tv%salt_deficit) .and. (IDs%id_salt_deficit > 0)) then
    do j=js,je ; do i=is,ie
      work_2d(i,j) = tv%salt_deficit(i,j) * I_time_int
    enddo ; enddo
    call post_data(IDs%id_salt_deficit, work_2d, diag, mask=G%mask2dT)
  endif

  ! post temperature of P-E+R
  if (associated(tv%TempxPmE) .and. (IDs%id_Heat_PmE > 0)) then
    do j=js,je ; do i=is,ie
      work_2d(i,j) = tv%TempxPmE(i,j) * (tv%C_p * I_time_int)
    enddo ; enddo
    call post_data(IDs%id_Heat_PmE, work_2d, diag, mask=G%mask2dT)
  endif

  ! post geothermal heating or internal heat source/sinks
  if (associated(tv%internal_heat) .and. (IDs%id_intern_heat > 0)) then
    do j=js,je ; do i=is,ie
      work_2d(i,j) = tv%internal_heat(i,j) * (tv%C_p * I_time_int)
    enddo ; enddo
    call post_data(IDs%id_intern_heat, work_2d, diag, mask=G%mask2dT)
  endif

  if (tv%T_is_conT) then
    ! Internal T&S variables are conservative temperature & absolute salinity
    if (IDs%id_sstcon > 0) call post_data(IDs%id_sstcon, sfc_state%SST, diag, mask=G%mask2dT)
    ! Use TEOS-10 function calls convert T&S diagnostics from conservative temp
    ! to potential temperature.
    EOSdom(:) = EOS_domain(G%HI)
    do j=js,je
      call cons_temp_to_pot_temp(sfc_state%SST(:,j), sfc_state%SSS(:,j), work_2d(:,j), tv%eqn_of_state, EOSdom)
    enddo
    if (IDs%id_sst > 0) call post_data(IDs%id_sst, work_2d, diag, mask=G%mask2dT)
  else
    ! Internal T&S variables are potential temperature & practical salinity
    if (IDs%id_sst > 0) call post_data(IDs%id_sst, sfc_state%SST, diag, mask=G%mask2dT)
  endif

  if (tv%S_is_absS) then
    ! Internal T&S variables are conservative temperature & absolute salinity
    if (IDs%id_sssabs > 0) call post_data(IDs%id_sssabs, sfc_state%SSS, diag, mask=G%mask2dT)
    ! Use TEOS-10 function calls convert T&S diagnostics from absolute salinity
    ! to practical salinity.
    EOSdom(:) = EOS_domain(G%HI)
    do j=js,je
      call abs_saln_to_prac_saln(sfc_state%SSS(:,j), work_2d(:,j), tv%eqn_of_state, EOSdom)
    enddo
    if (IDs%id_sss > 0) call post_data(IDs%id_sss, work_2d, diag, mask=G%mask2dT)
  else
    ! Internal T&S variables are potential temperature & practical salinity
    if (IDs%id_sss > 0) call post_data(IDs%id_sss, sfc_state%SSS, diag, mask=G%mask2dT)
  endif

  if (IDs%id_sst_sq > 0) then
    do j=js,je ; do i=is,ie
      work_2d(i,j) = sfc_state%SST(i,j)*sfc_state%SST(i,j)
    enddo ; enddo
    call post_data(IDs%id_sst_sq, work_2d, diag, mask=G%mask2dT)
  endif
  if (IDs%id_sss_sq > 0) then
    do j=js,je ; do i=is,ie
      work_2d(i,j) = sfc_state%SSS(i,j)*sfc_state%SSS(i,j)
    enddo ; enddo
    call post_data(IDs%id_sss_sq, work_2d, diag, mask=G%mask2dT)
  endif

  call coupler_type_send_data(sfc_state%tr_fields, get_diag_time_end(diag))

end subroutine post_surface_thermo_diags


!> This routine posts diagnostics of the transports, including the subgridscale
!! contributions.
subroutine post_transport_diagnostics(G, GV, US, uhtr, vhtr, h, IDs, diag_pre_dyn, &
                                      diag, dt_trans, Reg)
  type(ocean_grid_type),    intent(inout) :: G   !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV  !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uhtr !< Accumulated zonal thickness fluxes
                                                 !! used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: vhtr !< Accumulated meridional thickness fluxes
                                                 !! used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                            intent(in)    :: h   !< The updated layer thicknesses [H ~> m or kg m-2]
  type(transport_diag_IDs), intent(in)    :: IDs !< A structure with the diagnostic IDs.
  type(diag_grid_storage),  intent(inout) :: diag_pre_dyn !< Stored grids from before dynamics
  type(diag_ctrl),          intent(inout) :: diag !< regulates diagnostic output
  real,                     intent(in)    :: dt_trans !< total time step associated with the transports [T ~> s].
  type(tracer_registry_type), pointer     :: Reg !< Pointer to the tracer registry

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: umo2d ! Diagnostics of integrated mass transport [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: vmo2d ! Diagnostics of integrated mass transport [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: umo ! Diagnostics of layer mass transport [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vmo ! Diagnostics of layer mass transport [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))   :: h_tend ! Change in layer thickness due to dynamics
                          ! [H T-1 ~> m s-1 or kg m-2 s-1].
  real :: Idt             ! The inverse of the time interval [T-1 ~> s-1]
  real :: H_to_RZ_dt      ! A conversion factor from accumulated transports to fluxes
                          ! [R Z H-1 T-1 ~> kg m-3 s-1 or s-1].
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  Idt = 1. / dt_trans
  H_to_RZ_dt = GV%H_to_RZ * Idt

  call diag_save_grids(diag)
  call diag_copy_storage_to_diag(diag, diag_pre_dyn)

  if (IDs%id_umo_2d > 0) then
    umo2d(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      umo2d(I,j) = umo2d(I,j) + uhtr(I,j,k) * H_to_RZ_dt
    enddo ; enddo ; enddo
    call post_data(IDs%id_umo_2d, umo2d, diag)
  endif
  if (IDs%id_umo > 0) then
    ! Convert to kg/s.
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      umo(I,j,k) = uhtr(I,j,k) * H_to_RZ_dt
    enddo ; enddo ; enddo
    call post_data(IDs%id_umo, umo, diag, alt_h=diag_pre_dyn%h_state)
  endif
  if (IDs%id_vmo_2d > 0) then
    vmo2d(:,:) = 0.0
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      vmo2d(i,J) = vmo2d(i,J) + vhtr(i,J,k) * H_to_RZ_dt
    enddo ; enddo ; enddo
    call post_data(IDs%id_vmo_2d, vmo2d, diag)
  endif
  if (IDs%id_vmo > 0) then
    ! Convert to kg/s.
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      vmo(i,J,k) = vhtr(i,J,k) * H_to_RZ_dt
    enddo ; enddo ; enddo
    call post_data(IDs%id_vmo, vmo, diag, alt_h=diag_pre_dyn%h_state)
  endif

  if (IDs%id_uhtr > 0) call post_data(IDs%id_uhtr, uhtr, diag, alt_h=diag_pre_dyn%h_state)
  if (IDs%id_vhtr > 0) call post_data(IDs%id_vhtr, vhtr, diag, alt_h=diag_pre_dyn%h_state)
  if (IDs%id_dynamics_h > 0) call post_data(IDs%id_dynamics_h, diag_pre_dyn%h_state, diag, &
                                            alt_h=diag_pre_dyn%h_state)
  ! Post the change in thicknesses
  if (IDs%id_dynamics_h_tendency > 0) then
    h_tend(:,:,:) = 0.
    do k=1,nz ; do j=js,je ; do i=is,ie
      h_tend(i,j,k) = (h(i,j,k) - diag_pre_dyn%h_state(i,j,k))*Idt
    enddo ; enddo ; enddo
    call post_data(IDs%id_dynamics_h_tendency, h_tend, diag, alt_h=diag_pre_dyn%h_state)
  endif

  call post_tracer_transport_diagnostics(G, GV, Reg, diag_pre_dyn%h_state, diag)

  call diag_restore_grids(diag)

end subroutine post_transport_diagnostics

!> This subroutine registers various diagnostics and allocates space for fields
!! that other diagnostics depend upon.
subroutine MOM_diagnostics_init(MIS, ADp, CDp, Time, G, GV, US, param_file, diag, CS, tv)
  type(ocean_internal_state), intent(in)    :: MIS  !< For "MOM Internal State" a set of pointers to
                                                    !! the fields and accelerations that make up the
                                                    !! ocean's internal physical state.
  type(accel_diag_ptrs),      intent(inout) :: ADp  !< Structure with pointers to momentum equation
                                                    !! terms.
  type(cont_diag_ptrs),       intent(inout) :: CDp  !< Structure with pointers to continuity
                                                    !! equation terms.
  type(time_type),            intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),      intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file !< A structure to parse for run-time
                                                    !! parameters.
  type(diag_ctrl), target,    intent(inout) :: diag !< Structure to regulate diagnostic output.
  type(diagnostics_CS),       intent(inout) :: CS   !< Diagnostic control struct
  type(thermo_var_ptrs),      intent(in)    :: tv   !< A structure pointing to various
                                                    !! thermodynamic variables.

  ! Local variables
  real :: wave_speed_min      ! A floor in the first mode speed below which 0 is returned [L T-1 ~> m s-1]
  real :: wave_speed_tol      ! The fractional tolerance for finding the wave speeds [nondim]
  real :: convert_H           ! A conversion factor from internal thickness units to the appropriate
                              ! MKS units (m or kg m-2) for thicknesses depending on whether the
                              ! Boussinesq approximation is being made [m H-1 ~> 1] or [kg m-2 H-1 ~> 1]
  logical :: better_speed_est ! If true, use a more robust estimate of the first
                              ! mode wave speed as the starting point for iterations.
  logical :: split            ! True if using the barotropic-baroclinic split algorithm
  logical :: calc_tides       ! True if using tidal forcing
  logical :: calc_sal         ! True if using self-attraction and loading
  logical :: om4_remap_via_sub_cells ! Use the OM4-era ramap_via_sub_cells for calculating the EBT structure
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_diagnostics" ! This module's name.
  character(len=48) :: thickness_units, flux_units
  logical :: use_temperature, adiabatic
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  integer :: remap_answer_date    ! The vintage of the order of arithmetic and expressions to use
                                  ! for remapping.  Values below 20190101 recover the remapping
                                  ! answers from 2018, while higher values use more robust
                                  ! forms of the same remapping expressions.

  CS%initialized = .true.

  CS%diag => diag
  use_temperature = associated(tv%T)
  call get_param(param_file, mdl, "ADIABATIC", adiabatic, default=.false., do_not_log=.true.)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DIAG_EBT_MONO_N2_COLUMN_FRACTION", CS%mono_N2_column_fraction, &
                 "The lower fraction of water column over which N2 is limited as monotonic "// &
                 "for the purposes of calculating the equivalent barotropic wave speed.", &
                 units='nondim', default=0.)
  call get_param(param_file, mdl, "DIAG_EBT_MONO_N2_DEPTH", CS%mono_N2_depth, &
                 "The depth below which N2 is limited as monotonic for the "// &
                 "purposes of calculating the equivalent barotropic wave speed.", &
                 units='m', scale=GV%m_to_H, default=-1.)
  call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_TOL", wave_speed_tol, &
                 "The fractional tolerance for finding the wave speeds.", &
                 units="nondim", default=0.001)
  !### Set defaults so that wave_speed_min*wave_speed_tol >= 1e-9 m s-1
  call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_MIN", wave_speed_min, &
                 "A floor in the first mode speed below which 0 used instead.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_BETTER_EST", better_speed_est, &
                 "If true, use a more robust estimate of the first mode wave speed as the "//&
                 "starting point for iterations.", default=.true.)
  call get_param(param_file, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 do_not_log=.true., default=.true.)

  call get_param(param_file, mdl, "INTWAVE_REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 "If true, use the OM4 remapping-via-subcells algorithm for calculating EBT structure. "//&
                 "See REMAPPING_USE_OM4_SUBCELLS for details. "//&
                 "We recommend setting this option to false.", default=om4_remap_via_sub_cells)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "REMAPPING_ANSWER_DATE", remap_answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for remapping.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) remap_answer_date = max(remap_answer_date, 20230701)

  call get_param(param_file, mdl, "SPLIT", split, default=.true., do_not_log=.true.)
  call get_param(param_file, mdl, "TIDES", calc_tides, default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "CALCULATE_SAL", calc_sal, default=calc_tides, do_not_log=.true.)

  thickness_units = get_thickness_units(GV)
  flux_units = get_flux_units(GV)
  convert_H = GV%H_to_MKS

  CS%id_masscello = register_diag_field('ocean_model', 'masscello', diag%axesTL, &
      Time, 'Mass per unit area of liquid ocean grid cell', 'kg m-2', conversion=GV%H_to_kg_m2, &
      standard_name='sea_water_mass_per_unit_area', v_extensive=.true.)

  CS%id_masso = register_scalar_field('ocean_model', 'masso', Time, diag, &
      'Mass of liquid ocean', units='kg', conversion=US%RZL2_to_kg, &
      standard_name='sea_water_mass')

  CS%id_thkcello = register_diag_field('ocean_model', 'thkcello', diag%axesTL, Time, &
      long_name='Cell Thickness', standard_name='cell_thickness', &
      units='m', conversion=US%Z_to_m, v_extensive=.true.)
  CS%id_h_pre_sync = register_diag_field('ocean_model', 'h_pre_sync', diag%axesTL, Time, &
      long_name='Cell thickness from the previous timestep', &
      units=thickness_units, conversion=GV%H_to_MKS, v_extensive=.true.)

  ! Note that CS%id_volcello would normally be registered here but because it is a "cell measure" and
  ! must be registered first. We earlier stored the handle of volcello but need it here for posting
  ! by this module.
  CS%id_volcello = diag_get_volume_cell_measure_dm_id(diag)

  if (use_temperature) then
    if (tv%T_is_conT) then
      CS%id_Tpot = register_diag_field('ocean_model', 'temp', diag%axesTL, &
          Time, 'Potential Temperature', 'degC', conversion=US%C_to_degC, cmor_field_name="thetao")
    endif
    if (tv%S_is_absS) then
      CS%id_Sprac = register_diag_field('ocean_model', 'salt', diag%axesTL, &
          Time, 'Salinity', 'psu', conversion=US%S_to_ppt, cmor_field_name='so')
    endif

    CS%id_tob = register_diag_field('ocean_model','tob', diag%axesT1, Time, &
        long_name='Sea Water Potential Temperature at Sea Floor', &
        standard_name='sea_water_potential_temperature_at_sea_floor', &
        units='degC', conversion=US%C_to_degC)
    CS%id_sob = register_diag_field('ocean_model','sob',diag%axesT1, Time, &
        long_name='Sea Water Salinity at Sea Floor', &
        standard_name='sea_water_salinity_at_sea_floor', &
        units='psu', conversion=US%S_to_ppt)

    CS%id_tosq = register_diag_field('ocean_model', 'tosq', diag%axesTL, &
        Time, 'Square of Potential Temperature', 'degC2', conversion=US%C_to_degC**2, &
        standard_name='Potential Temperature Squared')
    CS%id_sosq = register_diag_field('ocean_model', 'sosq', diag%axesTL, &
        Time, 'Square of Salinity', 'psu2', conversion=US%S_to_ppt**2, &
        standard_name='Salinity Squared')

    CS%id_temp_layer_ave = register_diag_field('ocean_model', 'temp_layer_ave', &
        diag%axesZL, Time, 'Layer Average Ocean Temperature', units='degC', conversion=US%C_to_degC)
    CS%id_bigtemp_layer_ave = register_diag_field('ocean_model', 'contemp_layer_ave', &
        diag%axesZL, Time, 'Layer Average Ocean Conservative Temperature', units='Celsius', conversion=US%C_to_degC)
    CS%id_salt_layer_ave = register_diag_field('ocean_model', 'salt_layer_ave', &
        diag%axesZL, Time, 'Layer Average Ocean Salinity', units='psu', conversion=US%S_to_ppt)
    CS%id_abssalt_layer_ave = register_diag_field('ocean_model', 'abssalt_layer_ave', &
        diag%axesZL, Time, 'Layer Average Ocean Absolute Salinity', units='g kg-1', conversion=US%S_to_ppt)

    CS%id_thetaoga = register_scalar_field('ocean_model', 'thetaoga', &
        Time, diag, 'Global Mean Ocean Potential Temperature', units='degC', conversion=US%C_to_degC, &
        standard_name='sea_water_potential_temperature')
    CS%id_bigthetaoga = register_scalar_field('ocean_model', 'bigthetaoga', &
        Time, diag, 'Global Mean Ocean Conservative Temperature', units='Celsius', conversion=US%C_to_degC, &
        standard_name='sea_water_conservative_temperature')
    CS%id_soga = register_scalar_field('ocean_model', 'soga', &
        Time, diag, 'Global Mean Ocean Salinity', units='psu', conversion=US%S_to_ppt, &
        standard_name='sea_water_salinity')
    CS%id_abssoga = register_scalar_field('ocean_model', 'abssoga', &
        Time, diag, 'Global Mean Ocean Absolute Salinity', units='g kg-1', conversion=US%S_to_ppt, &
        standard_name='sea_water_absolute_salinity')

    ! The CMIP convention is potential temperature, but not indicated in the CMIP long name.
    CS%id_tosga = register_scalar_field('ocean_model', 'sst_global', Time, diag, &
        long_name='Global Area Average Sea Surface Temperature', &
        units='degC', conversion=US%C_to_degC, standard_name='sea_surface_temperature', &
        cmor_field_name='tosga', cmor_standard_name='sea_surface_temperature', &
        cmor_long_name='Sea Surface Temperature')
    CS%id_bigtosga = register_scalar_field('ocean_model', 'sscont_global', Time, diag, &
        long_name='Global Area Average Sea Surface Conservative Temperature', &
        units='Celsius', conversion=US%C_to_degC, standard_name='sea_surface_temperature')
    ! The CMIP convention is practical salinity, but not indicated in the CMIP long name.
    CS%id_sosga = register_scalar_field('ocean_model', 'sss_global', Time, diag, &
        long_name='Global Area Average Sea Surface Salinity', &
        units='psu', conversion=US%S_to_ppt, standard_name='sea_surface_salinity', &
        cmor_field_name='sosga', cmor_standard_name='sea_surface_salinity', &
        cmor_long_name='Sea Surface Salinity')
    CS%id_abssosga = register_scalar_field('ocean_model', 'ssabss_global', Time, diag, &
        long_name='Global Area Average Sea Surface Absolute Salinity', &
        units='psu', conversion=US%S_to_ppt, standard_name='sea_surface_absolute_salinity')
  endif

  CS%id_u = register_diag_field('ocean_model', 'u', diag%axesCuL, Time, &
      'Zonal velocity', 'm s-1', conversion=US%L_T_to_m_s, cmor_field_name='uo', &
      cmor_standard_name='sea_water_x_velocity', cmor_long_name='Sea Water X Velocity')
  CS%id_v = register_diag_field('ocean_model', 'v', diag%axesCvL, Time, &
      'Meridional velocity', 'm s-1', conversion=US%L_T_to_m_s, cmor_field_name='vo', &
      cmor_standard_name='sea_water_y_velocity', cmor_long_name='Sea Water Y Velocity')
  CS%id_usq = register_diag_field('ocean_model', 'usq', diag%axesCuL, Time, &
      'Zonal velocity squared', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_vsq = register_diag_field('ocean_model', 'vsq', diag%axesCvL, Time, &
      'Meridional velocity squared', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_uv = register_diag_field('ocean_model', 'uv', diag%axesTL, Time, &
      'Product between zonal and meridional velocities at h-points', &
      'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_h = register_diag_field('ocean_model', 'h', diag%axesTL, Time, &
      'Layer Thickness', thickness_units, v_extensive=.true., conversion=convert_H)

  CS%id_e = register_diag_field('ocean_model', 'e', diag%axesTi, Time, &
      'Interface Height Relative to Mean Sea Level', 'm', conversion=US%Z_to_m)
  CS%id_e_D = register_diag_field('ocean_model', 'e_D', diag%axesTi, Time, &
      'Interface Height above the Seafloor', 'm', conversion=US%Z_to_m)

  CS%id_Rml = register_diag_field('ocean_model', 'Rml', diag%axesTL, Time, &
      'Mixed Layer Coordinate Potential Density', 'kg m-3', conversion=US%R_to_kg_m3)

  CS%id_Rcv = register_diag_field('ocean_model', 'Rho_cv', diag%axesTL, Time, &
      'Coordinate Potential Density', 'kg m-3', conversion=US%R_to_kg_m3)

  CS%id_rhopot0 = register_diag_field('ocean_model', 'rhopot0', diag%axesTL, Time, &
      'Potential density referenced to surface', 'kg m-3', conversion=US%R_to_kg_m3)
  CS%id_rhopot2 = register_diag_field('ocean_model', 'rhopot2', diag%axesTL, Time, &
      'Potential density referenced to 2000 dbar', 'kg m-3', conversion=US%R_to_kg_m3)
  CS%id_rhoinsitu = register_diag_field('ocean_model', 'rhoinsitu', diag%axesTL, Time, &
      'In situ density', 'kg m-3', conversion=US%R_to_kg_m3)
  CS%id_drho_dT = register_diag_field('ocean_model', 'drho_dT', diag%axesTL, Time, &
      'Partial derivative of rhoinsitu with respect to temperature (alpha)', &
      'kg m-3 degC-1', conversion=US%R_to_kg_m3*US%degC_to_C)
  CS%id_drho_dS = register_diag_field('ocean_model', 'drho_dS', diag%axesTL, Time, &
      'Partial derivative of rhoinsitu with respect to salinity (beta)', &
      'kg^2 g-1 m-3', conversion=US%R_to_kg_m3*US%ppt_to_S)

  CS%id_du_dt = register_diag_field('ocean_model', 'dudt', diag%axesCuL, Time, &
      'Zonal Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_dv_dt = register_diag_field('ocean_model', 'dvdt', diag%axesCvL, Time, &
      'Meridional Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_dh_dt = register_diag_field('ocean_model', 'dhdt', diag%axesTL, Time, &
      'Thickness tendency', trim(thickness_units)//" s-1", conversion=convert_H*US%s_to_T, v_extensive=.true.)

  !CS%id_hf_du_dt = register_diag_field('ocean_model', 'hf_dudt', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Zonal Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2, &
  !    v_extensive=.true.)

  !CS%id_hf_dv_dt = register_diag_field('ocean_model', 'hf_dvdt', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2, &
  !    v_extensive=.true.)

  CS%id_hf_du_dt_2d = register_diag_field('ocean_model', 'hf_dudt_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_hf_dv_dt_2d = register_diag_field('ocean_model', 'hf_dvdt_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)

  CS%id_h_du_dt = register_diag_field('ocean_model', 'h_du_dt', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Acceleration', 'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  CS%id_h_dv_dt = register_diag_field('ocean_model', 'h_dv_dt', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Acceleration', 'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)

  ! layer thickness variables
  !if (GV%nk_rho_varies > 0) then
    CS%id_h_Rlay = register_diag_field('ocean_model', 'h_rho', diag%axesTL, Time, &
        'Layer thicknesses in pure potential density coordinates', &
        thickness_units, conversion=convert_H)

    CS%id_uh_Rlay = register_diag_field('ocean_model', 'uh_rho', diag%axesCuL, Time, &
        'Zonal volume transport in pure potential density coordinates', &
        flux_units, conversion=US%L_to_m**2*US%s_to_T*convert_H)

    CS%id_vh_Rlay = register_diag_field('ocean_model', 'vh_rho', diag%axesCvL, Time, &
        'Meridional volume transport in pure potential density coordinates', &
        flux_units, conversion=US%L_to_m**2*US%s_to_T*convert_H)

    CS%id_uhGM_Rlay = register_diag_field('ocean_model', 'uhGM_rho', diag%axesCuL, Time, &
        'Zonal volume transport due to interface height diffusion in pure potential '//&
        'density coordinates', flux_units, conversion=US%L_to_m**2*US%s_to_T*convert_H)

    CS%id_vhGM_Rlay = register_diag_field('ocean_model', 'vhGM_rho', diag%axesCvL, Time, &
        'Meridional volume transport due to interface height diffusion in pure potential '//&
        'density coordinates', flux_units, conversion=US%L_to_m**2*US%s_to_T*convert_H)
  !endif


  ! terms in the kinetic energy budget
  CS%id_KE = register_diag_field('ocean_model', 'KE', diag%axesTL, Time, &
      'Layer kinetic energy per unit mass', &
      'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_dKEdt = register_diag_field('ocean_model', 'dKE_dt', diag%axesTL, Time, &
      'Kinetic Energy Tendency of Layer', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_PE_to_KE = register_diag_field('ocean_model', 'PE_to_KE', diag%axesTL, Time, &
      'Potential to Kinetic Energy Conversion of Layer', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (calc_sal) &
    CS%id_KE_SAL = register_diag_field('ocean_model', 'KE_SAL', diag%axesTL, Time, &
        'Kinetic Energy Source from Self-Attraction and Loading', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (calc_tides) &
    CS%id_KE_TIDES = register_diag_field('ocean_model', 'KE_tides', diag%axesTL, Time, &
        'Kinetic Energy Source from Astronomical Tidal Forcing', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (split) then
    CS%id_KE_BT = register_diag_field('ocean_model', 'KE_BT', diag%axesTL, Time, &
        'Barotropic contribution to Kinetic Energy', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
    CS%id_PE_to_KE_btbc = register_diag_field('ocean_model', 'PE_to_KE_btbc', diag%axesTL, Time, &
        'Potential to Kinetic Energy Conversion of Layer (including barotropic solver contribution)', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
    CS%id_KE_Coradv_btbc = register_diag_field('ocean_model', 'KE_Coradv_btbc', diag%axesTL, Time, &
        'Kinetic Energy Source from Coriolis and Advection (including barotropic solver contribution)', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
    CS%id_KE_BT_PF = register_diag_field('ocean_model', 'KE_BTPF', diag%axesTL, Time, &
        'Kinetic Energy Source from Barotropic Pressure Gradient Force.', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
    CS%id_KE_BT_CF = register_diag_field('ocean_model', 'KE_BTCF', diag%axesTL, Time, &
        'Kinetic Energy Source from Barotropic Coriolis Force.', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
    CS%id_KE_BT_WD = register_diag_field('ocean_model', 'KE_BTWD', diag%axesTL, Time, &
        'Kinetic Energy Source from Barotropic Linear Wave Drag.', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  endif
  CS%id_KE_Coradv = register_diag_field('ocean_model', 'KE_Coradv', diag%axesTL, Time, &
      'Kinetic Energy Source from Coriolis and Advection', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_KE_adv = register_diag_field('ocean_model', 'KE_adv', diag%axesTL, Time, &
      'Kinetic Energy Source from Advection', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_KE_visc = register_diag_field('ocean_model', 'KE_visc', diag%axesTL, Time, &
      'Kinetic Energy Source from Vertical Viscosity and Stresses', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_KE_visc_gl90 = register_diag_field('ocean_model', 'KE_visc_gl90', diag%axesTL, Time, &
      'Kinetic Energy Source from GL90 Vertical Viscosity', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_KE_stress = register_diag_field('ocean_model', 'KE_stress', diag%axesTL, Time, &
      'Kinetic Energy Source from Surface Stresses or Body Wind Stress', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_KE_horvisc = register_diag_field('ocean_model', 'KE_horvisc', diag%axesTL, Time, &
      'Kinetic Energy Source from Horizontal Viscosity', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (.not. adiabatic) then
    CS%id_KE_dia = register_diag_field('ocean_model', 'KE_dia', diag%axesTL, Time, &
        'Kinetic Energy Source from Diapycnal Diffusion', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  endif

  ! gravity wave CFLs
  CS%id_cg1 = register_diag_field('ocean_model', 'cg1', diag%axesT1, Time, &
      'First baroclinic gravity wave speed', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_Rd1 = register_diag_field('ocean_model', 'Rd1', diag%axesT1, Time, &
      'First baroclinic deformation radius', 'm', conversion=US%L_to_m)
  CS%id_cfl_cg1 = register_diag_field('ocean_model', 'CFL_cg1', diag%axesT1, Time, &
      'CFL of first baroclinic gravity wave = dt*cg1*(1/dx+1/dy)', 'nondim')
  CS%id_cfl_cg1_x = register_diag_field('ocean_model', 'CFL_cg1_x', diag%axesT1, Time, &
      'i-component of CFL of first baroclinic gravity wave = dt*cg1*/dx', 'nondim')
  CS%id_cfl_cg1_y = register_diag_field('ocean_model', 'CFL_cg1_y', diag%axesT1, Time, &
      'j-component of CFL of first baroclinic gravity wave = dt*cg1*/dy', 'nondim')
  CS%id_cg_ebt = register_diag_field('ocean_model', 'cg_ebt', diag%axesT1, Time, &
      'Equivalent barotropic gravity wave speed', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_Rd_ebt = register_diag_field('ocean_model', 'Rd_ebt', diag%axesT1, Time, &
      'Equivalent barotropic deformation radius', 'm', conversion=US%L_to_m)
  CS%id_p_ebt = register_diag_field('ocean_model', 'p_ebt', diag%axesTL, Time, &
      'Equivalent barotropic modal strcuture', 'nondim')

  if ((CS%id_cg1>0) .or. (CS%id_Rd1>0) .or. (CS%id_cfl_cg1>0) .or. &
      (CS%id_cfl_cg1_x>0) .or. (CS%id_cfl_cg1_y>0) .or. &
      (CS%id_cg_ebt>0) .or. (CS%id_Rd_ebt>0) .or. (CS%id_p_ebt>0)) then
    call wave_speed_init(CS%wave_speed, GV, remap_answer_date=remap_answer_date, &
                         better_speed_est=better_speed_est, min_speed=wave_speed_min, &
                         wave_speed_tol=wave_speed_tol, om4_remap_via_sub_cells=om4_remap_via_sub_cells)
  endif

  CS%id_mass_wt = register_diag_field('ocean_model', 'mass_wt', diag%axesT1, Time, &
      'The column mass for calculating mass-weighted average properties', 'kg m-2', conversion=US%RZ_to_kg_m2)

  if (use_temperature) then
    CS%id_temp_int = register_diag_field('ocean_model', 'temp_int', diag%axesT1, Time, &
        'Density weighted column integrated potential temperature', &
        'degC kg m-2', conversion=US%C_to_degC*US%RZ_to_kg_m2, &
        cmor_field_name='opottempmint', &
        cmor_long_name='integral_wrt_depth_of_product_of_sea_water_density_and_potential_temperature', &
        cmor_standard_name='Depth integrated density times potential temperature')

    CS%id_salt_int = register_diag_field('ocean_model', 'salt_int', diag%axesT1, Time, &
        'Density weighted column integrated salinity', &
        'psu kg m-2', conversion=US%S_to_ppt*US%RZ_to_kg_m2, &
        cmor_field_name='somint', &
        cmor_long_name='integral_wrt_depth_of_product_of_sea_water_density_and_salinity', &
        cmor_standard_name='Depth integrated density times salinity')
  endif

  CS%id_col_mass = register_diag_field('ocean_model', 'col_mass', diag%axesT1, Time, &
      'The column integrated in situ density', 'kg m-2', conversion=US%RZ_to_kg_m2)

  CS%id_col_ht = register_diag_field('ocean_model', 'col_height', diag%axesT1, Time, &
      'The height of the water column', 'm', conversion=US%Z_to_m)
  CS%id_pbo = register_diag_field('ocean_model', 'pbo', diag%axesT1, Time, &
      long_name='Sea Water Pressure at Sea Floor', standard_name='sea_water_pressure_at_sea_floor', &
      units='Pa', conversion=US%RL2_T2_to_Pa)

  ! Register time derivatives and allocate memory for diagnostics that need
  ! access from across several modules.
  call set_dependent_diagnostics(MIS, ADp, CDp, G, GV, CS)

end subroutine MOM_diagnostics_init


!> Register diagnostics of the surface state and integrated quantities
subroutine register_surface_diags(Time, G, US, IDs, diag, tv)
  type(time_type),         intent(in)    :: Time  !< current model time
  type(ocean_grid_type),   intent(in)    :: G     !< ocean grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_diag_IDs),  intent(inout) :: IDs   !< A structure with the diagnostic IDs.
  type(diag_ctrl),         intent(inout) :: diag  !< regulates diagnostic output
  type(thermo_var_ptrs),   intent(in)    :: tv    !< A structure pointing to various thermodynamic variables

  ! Vertically integrated, budget, and surface state diagnostics
  IDs%id_volo = register_scalar_field('ocean_model', 'volo', Time, diag, &
      long_name='Total volume of liquid ocean', units='m3', conversion=US%Z_to_m*US%L_to_m**2, &
      standard_name='sea_water_volume')
  IDs%id_zos = register_diag_field('ocean_model', 'zos', diag%axesT1, Time, &
      standard_name = 'sea_surface_height_above_geoid', &
      long_name= 'Sea surface height above geoid', units='m', conversion=US%Z_to_m)
  IDs%id_zossq = register_diag_field('ocean_model', 'zossq', diag%axesT1, Time, &
      standard_name='square_of_sea_surface_height_above_geoid', &
      long_name='Square of sea surface height above geoid', units='m2', conversion=US%Z_to_m**2)
  IDs%id_ssh = register_diag_field('ocean_model', 'SSH', diag%axesT1, Time, &
      'Sea Surface Height', 'm', conversion=US%Z_to_m)
  IDs%id_ssh_ga = register_scalar_field('ocean_model', 'ssh_ga', Time, diag, &
      long_name='Area averaged sea surface height', units='m', conversion=US%Z_to_m, &
      standard_name='area_averaged_sea_surface_height')
  IDs%id_ssu = register_diag_field('ocean_model', 'SSU', diag%axesCu1, Time, &
      'Sea Surface Zonal Velocity', 'm s-1', conversion=US%L_T_to_m_s)
  IDs%id_ssv = register_diag_field('ocean_model', 'SSV', diag%axesCv1, Time, &
      'Sea Surface Meridional Velocity', 'm s-1', conversion=US%L_T_to_m_s)
  IDs%id_speed = register_diag_field('ocean_model', 'speed', diag%axesT1, Time, &
      'Sea Surface Speed', 'm s-1', conversion=US%L_T_to_m_s)
  IDs%id_ssu_east = register_diag_field('ocean_model', 'ssu_east', diag%axesT1, Time, &
      'Eastward velocity', 'm s-1', conversion=US%L_T_to_m_s)
  IDs%id_ssv_north = register_diag_field('ocean_model', 'ssv_north', diag%axesT1, Time, &
      'Northward velocity', 'm s-1', conversion=US%L_T_to_m_s)

  if (associated(tv%T)) then
    IDs%id_sst = register_diag_field('ocean_model', 'SST', diag%axesT1, Time, &
        'Sea Surface Temperature', 'degC', conversion=US%C_to_degC, &
        cmor_field_name='tos', cmor_long_name='Sea Surface Temperature', &
        cmor_standard_name='sea_surface_temperature')
    IDs%id_sst_sq = register_diag_field('ocean_model', 'SST_sq', diag%axesT1, Time, &
        'Sea Surface Temperature Squared', 'degC2', conversion=US%C_to_degC**2, &
        cmor_field_name='tossq', cmor_long_name='Square of Sea Surface Temperature ', &
        cmor_standard_name='square_of_sea_surface_temperature')
    IDs%id_sss = register_diag_field('ocean_model', 'SSS', diag%axesT1, Time, &
        'Sea Surface Salinity', 'psu', conversion=US%S_to_ppt, &
        cmor_field_name='sos', cmor_long_name='Sea Surface Salinity', &
        cmor_standard_name='sea_surface_salinity')
    IDs%id_sss_sq = register_diag_field('ocean_model', 'SSS_sq', diag%axesT1, Time, &
        'Sea Surface Salinity Squared', 'psu2', conversion=US%S_to_ppt**2, &
        cmor_field_name='sossq', cmor_long_name='Square of Sea Surface Salinity ', &
        cmor_standard_name='square_of_sea_surface_salinity')
    if (tv%T_is_conT) then
      IDs%id_sstcon = register_diag_field('ocean_model', 'conSST', diag%axesT1, Time, &
          'Sea Surface Conservative Temperature', 'Celsius', conversion=US%C_to_degC)
    endif
    if (tv%S_is_absS) then
      IDs%id_sssabs = register_diag_field('ocean_model', 'absSSS', diag%axesT1, Time, &
          'Sea Surface Absolute Salinity', 'g kg-1', conversion=US%S_to_ppt)
    endif
    if (associated(tv%frazil)) then
      IDs%id_fraz = register_diag_field('ocean_model', 'frazil', diag%axesT1, Time, &
            'Heat from frazil formation', 'W m-2', conversion=US%QRZ_T_to_W_m2, &
            cmor_field_name='hfsifrazil', &
            cmor_standard_name='heat_flux_into_sea_water_due_to_frazil_ice_formation', &
            cmor_long_name='Heat Flux into Sea Water due to Frazil Ice Formation')
    endif
  endif

  IDs%id_salt_deficit = register_diag_field('ocean_model', 'salt_deficit', diag%axesT1, Time, &
         'Salt source in ocean required to supply excessive ice salt fluxes', &
         'ppt kg m-2 s-1', conversion=US%S_to_ppt*US%RZ_T_to_kg_m2s)
  IDs%id_Heat_PmE = register_diag_field('ocean_model', 'Heat_PmE', diag%axesT1, Time, &
         'Heat flux into ocean from mass flux into ocean', &
         'W m-2', conversion=US%QRZ_T_to_W_m2)
  IDs%id_intern_heat = register_diag_field('ocean_model', 'internal_heat', diag%axesT1, Time, &
         'Heat flux into ocean from geothermal or other internal sources', &
         'W m-2', conversion=US%QRZ_T_to_W_m2)

end subroutine register_surface_diags

!> Register certain diagnostics related to transports
subroutine register_transport_diags(Time, G, GV, US, IDs, diag)
  type(time_type),          intent(in)    :: Time  !< current model time
  type(ocean_grid_type),    intent(in)    :: G     !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV    !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US    !< A dimensional unit scaling type
  type(transport_diag_IDs), intent(inout) :: IDs   !< A structure with the diagnostic IDs.
  type(diag_ctrl),          intent(inout) :: diag  !< regulates diagnostic output

  character(len=48) :: thickness_units, accum_flux_units

  thickness_units = get_thickness_units(GV)
  if (GV%Boussinesq) then
    accum_flux_units = "m3"
  else
    accum_flux_units = "kg"
  endif

  ! Diagnostics related to tracer and mass transport
  IDs%id_uhtr = register_diag_field('ocean_model', 'uhtr', diag%axesCuL, Time, &
      'Accumulated zonal thickness fluxes to advect tracers', &
      accum_flux_units, y_cell_method='sum', v_extensive=.true., conversion=GV%H_to_MKS*US%L_to_m**2)
  IDs%id_vhtr = register_diag_field('ocean_model', 'vhtr', diag%axesCvL, Time, &
      'Accumulated meridional thickness fluxes to advect tracers', &
      accum_flux_units, x_cell_method='sum', v_extensive=.true., conversion=GV%H_to_MKS*US%L_to_m**2)
  IDs%id_umo = register_diag_field('ocean_model', 'umo', &
      diag%axesCuL, Time, 'Ocean Mass X Transport', &
      'kg s-1', conversion=US%RZ_T_to_kg_m2s*US%L_to_m**2, &
      standard_name='ocean_mass_x_transport', y_cell_method='sum', v_extensive=.true.)
  IDs%id_vmo = register_diag_field('ocean_model', 'vmo', &
      diag%axesCvL, Time, 'Ocean Mass Y Transport', &
      'kg s-1', conversion=US%RZ_T_to_kg_m2s*US%L_to_m**2, &
      standard_name='ocean_mass_y_transport', x_cell_method='sum', v_extensive=.true.)
  IDs%id_umo_2d = register_diag_field('ocean_model', 'umo_2d', &
      diag%axesCu1, Time, 'Ocean Mass X Transport Vertical Sum', &
      'kg s-1', conversion=US%RZ_T_to_kg_m2s*US%L_to_m**2, &
      standard_name='ocean_mass_x_transport_vertical_sum', y_cell_method='sum')
  IDs%id_vmo_2d = register_diag_field('ocean_model', 'vmo_2d', &
      diag%axesCv1, Time, 'Ocean Mass Y Transport Vertical Sum', &
      'kg s-1', conversion=US%RZ_T_to_kg_m2s*US%L_to_m**2, &
      standard_name='ocean_mass_y_transport_vertical_sum', x_cell_method='sum')
  IDs%id_dynamics_h = register_diag_field('ocean_model','dynamics_h', &
      diag%axesTl, Time, 'Layer thicknesses prior to horizontal dynamics', &
      thickness_units, conversion=GV%H_to_MKS, v_extensive=.true.)
  IDs%id_dynamics_h_tendency = register_diag_field('ocean_model','dynamics_h_tendency', &
      diag%axesTl, Time, 'Change in layer thicknesses due to horizontal dynamics', &
      trim(thickness_units)//" s-1", conversion=GV%H_to_MKS*US%s_to_T, v_extensive=.true.)

end subroutine register_transport_diags

!> Offers the static fields in the ocean grid type for output via the diag_manager.
subroutine write_static_fields(G, GV, US, tv, diag)
  type(ocean_grid_type),   intent(in)    :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various thermodynamic variables
  type(diag_ctrl), target, intent(inout) :: diag !< regulates diagnostic output

  ! Local variables
  real :: work_2d(SZI_(G),SZJ_(G))         ! A 2-d temporary work array [Z ~> m]
  integer :: id, i, j
  logical :: use_temperature

  id = register_static_field('ocean_model', 'geolat', diag%axesT1, &
        'Latitude of tracer (T) points', 'degrees_north')
  if (id > 0) call post_data(id, G%geoLatT, diag, .true.)

  id = register_static_field('ocean_model', 'geolon', diag%axesT1, &
        'Longitude of tracer (T) points', 'degrees_east')
  if (id > 0) call post_data(id, G%geoLonT, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_c', diag%axesB1, &
        'Latitude of corner (Bu) points', 'degrees_north', interp_method='none')
  if (id > 0) call post_data(id, G%geoLatBu, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_c', diag%axesB1, &
        'Longitude of corner (Bu) points', 'degrees_east', interp_method='none')
  if (id > 0) call post_data(id, G%geoLonBu, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_v', diag%axesCv1, &
        'Latitude of meridional velocity (Cv) points', 'degrees_north', interp_method='none')
  if (id > 0) call post_data(id, G%geoLatCv, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_v', diag%axesCv1, &
        'Longitude of meridional velocity (Cv) points', 'degrees_east', interp_method='none')
  if (id > 0) call post_data(id, G%geoLonCv, diag, .true.)

  id = register_static_field('ocean_model', 'geolat_u', diag%axesCu1, &
        'Latitude of zonal velocity (Cu) points', 'degrees_north', interp_method='none')
  if (id > 0) call post_data(id, G%geoLatCu, diag, .true.)

  id = register_static_field('ocean_model', 'geolon_u', diag%axesCu1, &
        'Longitude of zonal velocity (Cu) points', 'degrees_east', interp_method='none')
  if (id > 0) call post_data(id, G%geoLonCu, diag, .true.)

  id = register_static_field('ocean_model', 'area_t', diag%axesT1, &
        'Surface area of tracer (T) cells', 'm2', conversion=US%L_to_m**2, &
        cmor_field_name='areacello', cmor_standard_name='cell_area', &
        cmor_long_name='Ocean Grid-Cell Area', &
        x_cell_method='sum', y_cell_method='sum', area_cell_method='sum')
  if (id > 0) then
    call post_data(id, G%areaT, diag, .true.)
    call diag_register_area_ids(diag, id_area_t=id)
  endif

  id = register_static_field('ocean_model', 'area_u', diag%axesCu1, &
        'Surface area of x-direction flow (U) cells', 'm2', conversion=US%L_to_m**2, &
        cmor_field_name='areacello_cu', cmor_standard_name='cell_area', &
        cmor_long_name='Ocean Grid-Cell Area', &
        x_cell_method='sum', y_cell_method='sum', area_cell_method='sum')
  if (id > 0) call post_data(id, G%areaCu, diag, .true.)

  id = register_static_field('ocean_model', 'area_v', diag%axesCv1, &
        'Surface area of y-direction flow (V) cells', 'm2', conversion=US%L_to_m**2, &
        cmor_field_name='areacello_cv', cmor_standard_name='cell_area', &
        cmor_long_name='Ocean Grid-Cell Area', &
        x_cell_method='sum', y_cell_method='sum', area_cell_method='sum')
  if (id > 0) call post_data(id, G%areaCv, diag, .true.)

  id = register_static_field('ocean_model', 'area_q', diag%axesB1, &
        'Surface area of B-grid flow (Q) cells', 'm2', conversion=US%L_to_m**2, &
        cmor_field_name='areacello_bu', cmor_standard_name='cell_area', &
        cmor_long_name='Ocean Grid-Cell Area', &
        x_cell_method='sum', y_cell_method='sum', area_cell_method='sum')
  if (id > 0) call post_data(id, G%areaBu, diag, .true.)

  id = register_static_field('ocean_model', 'depth_ocean', diag%axesT1, &
        'Depth of the ocean at tracer points', 'm', conversion=US%Z_to_m, &
        standard_name='sea_floor_depth_below_geoid', &
        cmor_field_name='deptho', cmor_long_name='Sea Floor Depth', &
        cmor_standard_name='sea_floor_depth_below_geoid', area=diag%axesT1%id_area, &
        x_cell_method='mean', y_cell_method='mean', area_cell_method='mean')
  if (id > 0) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec ; work_2d(i,j) = G%bathyT(i,j)+G%Z_ref ; enddo ; enddo
    call post_data(id, work_2d, diag, .true., mask=G%mask2dT)
  endif

  id = register_static_field('ocean_model', 'wet', diag%axesT1, &
        '0 if land, 1 if ocean at tracer points', 'none', area=diag%axesT1%id_area)
  if (id > 0) call post_data(id, G%mask2dT, diag, .true.)

  id = register_static_field('ocean_model', 'wet_c', diag%axesB1, &
        '0 if land, 1 if ocean at corner (Bu) points', 'none', interp_method='none')
  if (id > 0) call post_data(id, G%mask2dBu, diag, .true.)

  id = register_static_field('ocean_model', 'wet_u', diag%axesCu1, &
        '0 if land, 1 if ocean at zonal velocity (Cu) points', 'none', interp_method='none')
  if (id > 0) call post_data(id, G%mask2dCu, diag, .true.)

  id = register_static_field('ocean_model', 'wet_v', diag%axesCv1, &
        '0 if land, 1 if ocean at meridional velocity (Cv) points', 'none', interp_method='none')
  if (id > 0) call post_data(id, G%mask2dCv, diag, .true.)

  id = register_static_field('ocean_model', 'Coriolis', diag%axesB1, &
        'Coriolis parameter at corner (Bu) points', 's-1', interp_method='none', conversion=US%s_to_T)
  if (id > 0) call post_data(id, G%CoriolisBu, diag, .true.)

  id = register_static_field('ocean_model', 'dxt', diag%axesT1, &
        'Delta(x) at thickness/tracer points (meter)', 'm', interp_method='none', conversion=US%L_to_m)
  if (id > 0) call post_data(id, G%dxT, diag, .true.)

  id = register_static_field('ocean_model', 'dyt', diag%axesT1, &
        'Delta(y) at thickness/tracer points (meter)', 'm', interp_method='none', conversion=US%L_to_m)
  if (id > 0) call post_data(id, G%dyT, diag, .true.)

  id = register_static_field('ocean_model', 'dxCu', diag%axesCu1, &
        'Delta(x) at u points (meter)', 'm', interp_method='none', conversion=US%L_to_m)
  if (id > 0) call post_data(id, G%dxCu, diag, .true.)

  id = register_static_field('ocean_model', 'dyCu', diag%axesCu1, &
        'Delta(y) at u points (meter)', 'm', interp_method='none', conversion=US%L_to_m)
  if (id > 0) call post_data(id, G%dyCu, diag, .true.)

  id = register_static_field('ocean_model', 'dxCv', diag%axesCv1, &
        'Delta(x) at v points (meter)', 'm', interp_method='none', conversion=US%L_to_m)
  if (id > 0) call post_data(id, G%dxCv, diag, .true.)

  id = register_static_field('ocean_model', 'dyCv', diag%axesCv1, &
        'Delta(y) at v points (meter)', 'm', interp_method='none', conversion=US%L_to_m)
  if (id > 0) call post_data(id, G%dyCv, diag, .true.)

  id = register_static_field('ocean_model', 'dyCuo', diag%axesCu1, &
        'Open meridional grid spacing at u points (meter)', 'm', interp_method='none', conversion=US%L_to_m)
  if (id > 0) call post_data(id, G%dy_Cu, diag, .true.)

  id = register_static_field('ocean_model', 'dxCvo', diag%axesCv1, &
        'Open zonal grid spacing at v points (meter)', 'm', interp_method='none', conversion=US%L_to_m)
  if (id > 0) call post_data(id, G%dx_Cv, diag, .true.)

  id = register_static_field('ocean_model', 'sin_rot', diag%axesT1, &
        'sine of the clockwise angle of the ocean grid north to true north', 'none')
  if (id > 0) call post_data(id, G%sin_rot, diag, .true.)

  id = register_static_field('ocean_model', 'cos_rot', diag%axesT1, &
        'cosine of the clockwise angle of the ocean grid north to true north', 'none')
  if (id > 0) call post_data(id, G%cos_rot, diag, .true.)


  ! This static diagnostic is from CF 1.8, and is the fraction of a cell
  ! covered by ocean, given as a percentage (poorly named).
  id = register_static_field('ocean_model', 'area_t_percent', diag%axesT1, &
        'Percentage of cell area covered by ocean', '%', conversion=100.0, &
        cmor_field_name='sftof', cmor_standard_name='SeaAreaFraction', &
        cmor_long_name='Sea Area Fraction', &
        x_cell_method='mean', y_cell_method='mean', area_cell_method='mean')
  if (id > 0) call post_data(id, G%mask2dT, diag, .true.)


  id = register_static_field('ocean_model','Rho_0', diag%axesNull, &
       'mean ocean density used with the Boussinesq approximation', &
       'kg m-3', conversion=US%R_to_kg_m3, cmor_field_name='rhozero', &
       cmor_standard_name='reference_sea_water_density_for_boussinesq_approximation', &
       cmor_long_name='reference sea water density for boussinesq approximation')
  if (id > 0) call post_data(id, GV%Rho0, diag, .true.)

  use_temperature = associated(tv%T)
  if (use_temperature) then
    id = register_static_field('ocean_model','C_p', diag%axesNull, &
         'heat capacity of sea water', 'J kg-1 K-1', conversion=US%Q_to_J_kg*US%degC_to_C, &
         cmor_field_name='cpocean', &
         cmor_standard_name='specific_heat_capacity_of_sea_water', &
         cmor_long_name='specific_heat_capacity_of_sea_water')
    if (id > 0) call post_data(id, tv%C_p, diag, .true.)
  endif

end subroutine write_static_fields


!> This subroutine sets up diagnostics upon which other diagnostics depend.
subroutine set_dependent_diagnostics(MIS, ADp, CDp, G, GV, CS)
  type(ocean_internal_state), intent(in)    :: MIS !< For "MOM Internal State" a set of pointers to
                                                   !! the fields and accelerations making up ocean
                                                   !! internal physical state.
  type(accel_diag_ptrs),      intent(inout) :: ADp !< Structure pointing to accelerations in
                                                   !! momentum equation.
  type(cont_diag_ptrs),       intent(inout) :: CDp !< Structure pointing to terms in continuity
                                                   !! equation.
  type(ocean_grid_type),      intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV   !< ocean vertical grid structure
  type(diagnostics_CS),       intent(inout) :: CS  !< Pointer to the control structure for this
                                                   !! module.

  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! Allocate and register time derivatives.
  if ( ( (CS%id_du_dt>0) .or. (CS%id_dKEdt > 0) .or. &
       ! (CS%id_hf_du_dt > 0) .or. &
         (CS%id_h_du_dt > 0) .or. (CS%id_hf_du_dt_2d > 0) ) .and. &
       (.not. allocated(CS%du_dt)) ) then
    allocate(CS%du_dt(IsdB:IedB,jsd:jed,nz), source=0.)
    call register_time_deriv(lbound(MIS%u), MIS%u, CS%du_dt, CS)
  endif
  if ( ( (CS%id_dv_dt>0) .or. (CS%id_dKEdt > 0) .or. &
       ! (CS%id_hf_dv_dt > 0) .or. &
         (CS%id_h_dv_dt > 0) .or. (CS%id_hf_dv_dt_2d > 0) ) .and. &
       (.not. allocated(CS%dv_dt)) ) then
    allocate(CS%dv_dt(isd:ied,JsdB:JedB,nz), source=0.)
    call register_time_deriv(lbound(MIS%v), MIS%v, CS%dv_dt, CS)
  endif
  if ( ( (CS%id_dh_dt>0) .or. (CS%id_dKEdt > 0) ) .and. &
       (.not. allocated(CS%dh_dt)) ) then
    allocate(CS%dh_dt(isd:ied,jsd:jed,nz), source=0.)
    call register_time_deriv(lbound(MIS%h), MIS%h, CS%dh_dt, CS)
  endif

  ! Allocate memory for other dependent diagnostics.
  if (CS%id_KE_adv > 0) then
    call safe_alloc_ptr(ADp%gradKEu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%gradKEv,isd,ied,JsdB,JedB,nz)
  endif
  if (CS%id_KE_visc > 0) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  endif
  if (CS%id_KE_visc_gl90 > 0) then
    call safe_alloc_ptr(ADp%du_dt_visc_gl90,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_visc_gl90,isd,ied,JsdB,JedB,nz)
  endif
  if (CS%id_KE_stress > 0) then
    call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)
  endif

  if (CS%id_KE_dia > 0) then
    call safe_alloc_ptr(ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)
  endif

  if ((CS%id_PE_to_KE_btbc > 0) .or. (CS%id_KE_BT_PF > 0)) then
    call safe_alloc_ptr(ADp%bt_pgf_u, IsdB, IedB, jsd, jed, nz)
    call safe_alloc_ptr(ADp%bt_pgf_v, isd, ied, JsdB, JedB, nz)
  endif

  if ((CS%id_KE_Coradv_btbc > 0) .or. (CS%id_KE_BT_CF > 0)) then
    call safe_alloc_ptr(ADp%bt_cor_u, IsdB, IedB, jsd, jed)
    call safe_alloc_ptr(ADp%bt_cor_v, isd, ied, JsdB, JedB)
  endif

  if (CS%id_KE_BT_WD > 0) then
    call safe_alloc_ptr(ADp%bt_lwd_u, IsdB, IedB, jsd, jed)
    call safe_alloc_ptr(ADp%bt_lwd_v, isd, ied, JsdB, JedB)
  endif

  if (CS%id_KE_SAL > 0) then
    call safe_alloc_ptr(ADp%sal_u, IsdB, IedB, jsd, jed, nz)
    call safe_alloc_ptr(ADp%sal_v, isd, ied, JsdB, JedB, nz)
  endif

  if (CS%id_KE_TIDES > 0) then
    call safe_alloc_ptr(ADp%tides_u, IsdB, IedB, jsd, jed, nz)
    call safe_alloc_ptr(ADp%tides_v, isd, ied, JsdB, JedB, nz)
  endif

  CS%KE_term_on = ((CS%id_dKEdt > 0) .or. (CS%id_PE_to_KE > 0) .or. (CS%id_KE_BT > 0) .or. &
                   (CS%id_KE_Coradv > 0) .or. (CS%id_KE_adv > 0) .or. (CS%id_KE_visc > 0) .or. &
                   (CS%id_KE_visc_gl90 > 0) .or. (CS%id_KE_stress > 0) .or. (CS%id_KE_horvisc > 0) .or. &
                   (CS%id_KE_dia > 0) .or. (CS%id_PE_to_KE_btbc > 0) .or. (CS%id_KE_BT_PF > 0) .or. &
                   (CS%id_KE_Coradv_btbc > 0) .or. (CS%id_KE_BT_CF > 0) .or. (CS%id_KE_BT_WD > 0) .or. &
                   (CS%id_KE_SAL > 0) .or. (CS%id_KE_TIDES > 0))

  if (CS%id_h_du_dt > 0) call safe_alloc_ptr(ADp%diag_hu,IsdB,IedB,jsd,jed,nz)
  if (CS%id_h_dv_dt > 0) call safe_alloc_ptr(ADp%diag_hv,isd,ied,JsdB,JedB,nz)

  if (CS%id_hf_du_dt_2d > 0) call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  if (CS%id_hf_dv_dt_2d > 0) call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)
  ! if (CS%id_hf_du_dt > 0) call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  ! if (CS%id_hf_dv_dt > 0) call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)

  if (CS%id_uhGM_Rlay > 0) call safe_alloc_ptr(CDp%uhGM,IsdB,IedB,jsd,jed,nz)
  if (CS%id_vhGM_Rlay > 0) call safe_alloc_ptr(CDp%vhGM,isd,ied,JsdB,JedB,nz)

end subroutine set_dependent_diagnostics

!> Deallocate memory associated with the diagnostics module
subroutine MOM_diagnostics_end(CS, ADp, CDp)
  type(diagnostics_CS),  intent(inout) :: CS  !< Control structure returned by a
                                              !! previous call to diagnostics_init.
  type(accel_diag_ptrs), intent(inout) :: ADp !< structure with pointers to
                                              !! accelerations in momentum equation.
  type(cont_diag_ptrs),  intent(inout) :: CDp !< Structure pointing to terms in continuity
                                              !! equation.
  integer :: m

  if (allocated(CS%dh_dt))      deallocate(CS%dh_dt)
  if (allocated(CS%dv_dt))      deallocate(CS%dv_dt)
  if (allocated(CS%du_dt))      deallocate(CS%du_dt)

  if (associated(ADp%gradKEu))    deallocate(ADp%gradKEu)
  if (associated(ADp%gradKEv))    deallocate(ADp%gradKEv)
  if (associated(ADp%du_dt_visc)) deallocate(ADp%du_dt_visc)
  if (associated(ADp%dv_dt_visc)) deallocate(ADp%dv_dt_visc)
  if (associated(ADp%du_dt_str))  deallocate(ADp%du_dt_str)
  if (associated(ADp%dv_dt_str))  deallocate(ADp%dv_dt_str)
  if (associated(ADp%du_dt_dia))  deallocate(ADp%du_dt_dia)
  if (associated(ADp%dv_dt_dia))  deallocate(ADp%dv_dt_dia)
  if (associated(ADp%du_other))   deallocate(ADp%du_other)
  if (associated(ADp%dv_other))   deallocate(ADp%dv_other)

  if (associated(ADp%bt_pgf_u)) deallocate(ADp%bt_pgf_u)
  if (associated(ADp%bt_pgf_v)) deallocate(ADp%bt_pgf_v)
  if (associated(ADp%bt_cor_u)) deallocate(ADp%bt_cor_u)
  if (associated(ADp%bt_cor_v)) deallocate(ADp%bt_cor_v)
  if (associated(ADp%bt_lwd_u)) deallocate(ADp%bt_lwd_u)
  if (associated(ADp%bt_lwd_v)) deallocate(ADp%bt_lwd_v)

  ! NOTE: sal_[uv] and tide_[uv] may be allocated either here (KE budget diagnostics) or
  ! PressureForce module (momentum acceleration diagnostics)
  if (associated(ADp%sal_u))  deallocate(ADp%sal_u)
  if (associated(ADp%sal_v))  deallocate(ADp%sal_v)
  if (associated(ADp%tides_u)) deallocate(ADp%tides_u)
  if (associated(ADp%tides_v)) deallocate(ADp%tides_v)

  if (associated(ADp%diag_hfrac_u)) deallocate(ADp%diag_hfrac_u)
  if (associated(ADp%diag_hfrac_v)) deallocate(ADp%diag_hfrac_v)

  ! NOTE: [uv]hGM may be allocated either here or the thickness diffuse module
  if (associated(CDp%uhGM)) deallocate(CDp%uhGM)
  if (associated(CDp%vhGM)) deallocate(CDp%vhGM)
  if (associated(CDp%diapyc_vel)) deallocate(CDp%diapyc_vel)

  do m=1,CS%num_time_deriv ; deallocate(CS%prev_val(m)%p) ; enddo
end subroutine MOM_diagnostics_end

end module MOM_diagnostics
