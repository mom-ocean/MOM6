!> Calculates any requested diagnostic quantities
!! that are not calculated in the various subroutines.
!! Diagnostic quantities are requested by allocating them memory.
module MOM_diagnostics

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,              only : reproducing_sum
use MOM_coupler_types,     only : coupler_type_send_data
use MOM_density_integrals, only : int_density_dz
use MOM_diag_mediator,     only : post_data, get_diag_time_end
use MOM_diag_mediator,     only : register_diag_field, register_scalar_field
use MOM_diag_mediator,     only : register_static_field, diag_register_area_ids
use MOM_diag_mediator,     only : diag_ctrl, time_type, safe_alloc_ptr
use MOM_diag_mediator,     only : diag_get_volume_cell_measure_dm_id
use MOM_diag_mediator,     only : diag_grid_storage
use MOM_diag_mediator,     only : diag_save_grids, diag_restore_grids, diag_copy_storage_to_diag
use MOM_domains,           only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,           only : To_North, To_East
use MOM_EOS,               only : calculate_density, calculate_density_derivs, EOS_domain
use MOM_EOS,               only : gsw_sp_from_sr, gsw_pt_from_ct
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
use MOM_verticalGrid,      only : verticalGrid_type, get_thickness_units
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
  real :: mono_N2_column_fraction = 0. !< The lower fraction of water column over which N2 is limited as
                                       !! monotonic for the purposes of calculating the equivalent
                                       !! barotropic wave speed.
  real :: mono_N2_depth = -1.          !< The depth below which N2 is limited as monotonic for the purposes of
                                       !! calculating the equivalent barotropic wave speed [Z ~> m].

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                       !! regulate the timing of diagnostic output.

  ! following arrays store diagnostics calculated here and unavailable outside.

  ! following fields have nz+1 levels.
  real, pointer, dimension(:,:,:) :: &
    e => NULL(), &   !< interface height [Z ~> m]
    e_D => NULL()    !< interface height above bottom [Z ~> m]

  ! following fields have nz layers.
  real, pointer, dimension(:,:,:) :: &
    du_dt => NULL(), & !< net i-acceleration [L T-2 ~> m s-2]
    dv_dt => NULL(), & !< net j-acceleration [L T-2 ~> m s-2]
    dh_dt => NULL(), & !< thickness rate of change [H T-1 ~> m s-1 or kg m-2 s-1]
    p_ebt => NULL()    !< Equivalent barotropic modal structure [nondim]
    ! hf_du_dt => NULL(), hf_dv_dt => NULL() !< du_dt, dv_dt x fract. thickness [L T-2 ~> m s-2].
    ! 3D diagnostics hf_du(dv)_dt are commented because there is no clarity on proper remapping grid option.
    ! The code is retained for degugging purposes in the future.

  real, pointer, dimension(:,:,:) :: h_Rlay => NULL() !< Layer thicknesses in potential density
                                              !! coordinates [H ~> m or kg m-2]
  real, pointer, dimension(:,:,:) :: uh_Rlay => NULL() !< Zonal transports in potential density
                                              !! coordinates [H m2 s-1 ~> m3 s-1 or kg s-1]
  real, pointer, dimension(:,:,:) :: vh_Rlay => NULL() !< Meridional transports in potential density
                                              !! coordinates [H m2 s-1 ~> m3 s-1 or kg s-1]
  real, pointer, dimension(:,:,:) :: uhGM_Rlay => NULL() !< Zonal Gent-McWilliams transports in potential density
                                              !! coordinates [H m2 s-1 ~> m3 s-1 or kg s-1]
  real, pointer, dimension(:,:,:) :: vhGM_Rlay => NULL() !< Meridional Gent-McWilliams transports in potential density
                                              !! coordinates [H m2 s-1 ~> m3 s-1 or kg s-1]

  ! following fields are 2-D.
  real, pointer, dimension(:,:) :: &
    cg1 => NULL(),       & !< First baroclinic gravity wave speed [L T-1 ~> m s-1]
    Rd1 => NULL(),       & !< First baroclinic deformation radius [L ~> m]
    cfl_cg1 => NULL(),   & !< CFL for first baroclinic gravity wave speed [nondim]
    cfl_cg1_x => NULL(), & !< i-component of CFL for first baroclinic gravity wave speed [nondim]
    cfl_cg1_y => NULL()    !< j-component of CFL for first baroclinic gravity wave speed [nondim]

  ! The following arrays hold diagnostics in the layer-integrated energy budget.
  real, pointer, dimension(:,:,:) :: &
    KE        => NULL(), &  !< KE per unit mass [L2 T-2 ~> m2 s-2]
    dKE_dt    => NULL(), &  !< time derivative of the layer KE [H L2 T-3 ~> m3 s-3]
    PE_to_KE  => NULL(), &  !< potential energy to KE term [m3 s-3]
    KE_BT     => NULL(), &  !< barotropic contribution to KE term [m3 s-3]
    KE_CorAdv => NULL(), &  !< KE source from the combined Coriolis and
                            !! advection terms [H L2 T-3 ~> m3 s-3].
                            !! The Coriolis source should be zero, but is not due to truncation
                            !! errors.  There should be near-cancellation of the global integral
                            !! of this spurious Coriolis source.
    KE_adv     => NULL(), & !< KE source from along-layer advection [H L2 T-3 ~> m3 s-3]
    KE_visc    => NULL(), & !< KE source from vertical viscosity [H L2 T-3 ~> m3 s-3]
    KE_stress  => NULL(), & !< KE source from surface stress (included in KE_visc) [H L2 T-3 ~> m3 s-3]
    KE_horvisc => NULL(), & !< KE source from horizontal viscosity [H L2 T-3 ~> m3 s-3]
    KE_dia     => NULL()    !< KE source from diapycnal diffusion [H L2 T-3 ~> m3 s-3]

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
  integer :: id_KE_Coradv      = -1, id_KE_adv         = -1
  integer :: id_KE_visc        = -1, id_KE_stress      = -1
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
  integer :: id_sosga          = -1, id_tosga          = -1
  integer :: id_temp_layer_ave = -1, id_salt_layer_ave = -1
  integer :: id_pbo            = -1
  integer :: id_thkcello       = -1, id_rhoinsitu      = -1
  integer :: id_rhopot0        = -1, id_rhopot2        = -1
  integer :: id_drho_dT        = -1, id_drho_dS        = -1
  integer :: id_h_pre_sync     = -1
  !>@}
  !> The control structure for calculating wave speed.
  type(wave_speed_CS), pointer :: wave_speed_CSp => NULL()

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
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: usq ! squared eastward velocity  [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: vsq ! squared northward velocity [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G))  :: uv  ! u x v at h-points          [L2 T-2 ~> m2 s-2]

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb

  real :: Rcv(SZI_(G),SZJ_(G),SZK_(GV)) ! Coordinate variable potential density [R ~> kg m-3].
  real :: work_3d(SZI_(G),SZJ_(G),SZK_(GV)) ! A 3-d temporary work array.
  real :: work_2d(SZI_(G),SZJ_(G))         ! A 2-d temporary work array.
  real :: rho_in_situ(SZI_(G))             ! In situ density [R ~> kg m-3]

  real, allocatable, dimension(:,:) :: &
    hf_du_dt_2d, hf_dv_dt_2d ! z integeral of hf_du_dt, hf_dv_dt [L T-2 ~> m s-2].

  real, allocatable, dimension(:,:,:) :: h_du_dt ! h x dudt [H L T-2 ~> m2 s-2]
  real, allocatable, dimension(:,:,:) :: h_dv_dt ! h x dvdt [H L T-2 ~> m2 s-2]

  ! tmp array for surface properties
  real :: surface_field(SZI_(G),SZJ_(G))
  real :: pressure_1d(SZI_(G)) ! Temporary array for pressure when calling EOS [R L2 T-2 ~> Pa]
  real :: wt, wt_p

  real :: f2_h     ! Squared Coriolis parameter at to h-points [T-2 ~> s-2]
  real :: mag_beta ! Magnitude of the gradient of f [T-1 L-1 ~> s-1 m-1]
  real :: absurdly_small_freq2 ! Frequency squared used to avoid division by 0 [T-2 ~> s-2]

  integer :: k_list

  real, dimension(SZK_(GV)) :: temp_layer_ave, salt_layer_ave
  real :: thetaoga, soga, masso, tosga, sosga

  is  = G%isc  ; ie   = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq  = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nz  = GV%ke   ; nkmb = GV%nk_rho_varies

  ! This value is roughly (pi / (the age of the universe) )^2.
  absurdly_small_freq2 = 1e-34*US%T_to_s**2

  if (loc(CS)==0) call MOM_error(FATAL, &
         "calculate_diagnostic_fields: Module must be initialized before used.")

  call calculate_derivs(dt, G, CS)

  if (dt > 0.0) then
    call diag_save_grids(CS%diag)
    call diag_copy_storage_to_diag(CS%diag, diag_pre_sync)

    if (CS%id_h_pre_sync > 0) &
        call post_data(CS%id_h_pre_sync, diag_pre_sync%h_state, CS%diag, alt_h = diag_pre_sync%h_state)

    if (CS%id_du_dt>0) call post_data(CS%id_du_dt, CS%du_dt, CS%diag, alt_h = diag_pre_sync%h_state)

    if (CS%id_dv_dt>0) call post_data(CS%id_dv_dt, CS%dv_dt, CS%diag, alt_h = diag_pre_sync%h_state)

    if (CS%id_dh_dt>0) call post_data(CS%id_dh_dt, CS%dh_dt, CS%diag, alt_h = diag_pre_sync%h_state)

    !! Diagnostics for terms multiplied by fractional thicknesses

    ! 3D diagnostics hf_du(dv)_dt are commented because there is no clarity on proper remapping grid option.
    ! The code is retained for degugging purposes in the future.
    !if (CS%id_hf_du_dt > 0) then
    !  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    !    CS%hf_du_dt(I,j,k) = CS%du_dt(I,j,k) * ADp%diag_hfrac_u(I,j,k)
    !  enddo ; enddo ; enddo
    !  call post_data(CS%id_hf_du_dt, CS%hf_du_dt, CS%diag, alt_h = diag_pre_sync%h_state)
    !endif

    !if (CS%id_hf_dv_dt > 0) then
    !  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    !    CS%hf_dv_dt(i,J,k) = CS%dv_dt(i,J,k) * ADp%diag_hfrac_v(i,J,k)
    !  enddo ; enddo ; enddo
    !  call post_data(CS%id_hf_dv_dt, CS%hf_dv_dt, CS%diag, alt_h = diag_pre_sync%h_state)
    !endif

    if (CS%id_hf_du_dt_2d > 0) then
      allocate(hf_du_dt_2d(G%IsdB:G%IedB,G%jsd:G%jed))
      hf_du_dt_2d(:,:) = 0.0
      do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        hf_du_dt_2d(I,j) = hf_du_dt_2d(I,j) + CS%du_dt(I,j,k) * ADp%diag_hfrac_u(I,j,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_hf_du_dt_2d, hf_du_dt_2d, CS%diag)
      deallocate(hf_du_dt_2d)
    endif

    if (CS%id_hf_dv_dt_2d > 0) then
      allocate(hf_dv_dt_2d(G%isd:G%ied,G%JsdB:G%JedB))
      hf_dv_dt_2d(:,:) = 0.0
      do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        hf_dv_dt_2d(i,J) = hf_dv_dt_2d(i,J) + CS%dv_dt(i,J,k) * ADp%diag_hfrac_v(i,J,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_hf_dv_dt_2d, hf_dv_dt_2d, CS%diag)
      deallocate(hf_dv_dt_2d)
    endif

    if (CS%id_h_du_dt > 0) then
      allocate(h_du_dt(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke))
      h_du_dt(:,:,:) = 0.0
      do k=1,nz ; do j=js,je ; do I=Isq,Ieq
        h_du_dt(I,j,k) = CS%du_dt(I,j,k) * ADp%diag_hu(I,j,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_h_du_dt, h_du_dt, CS%diag)
      deallocate(h_du_dt)
    endif
    if (CS%id_h_dv_dt > 0) then
      allocate(h_dv_dt(G%isd:G%ied,G%JsdB:G%JedB,GV%ke))
      h_dv_dt(:,:,:) = 0.0
      do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        h_dv_dt(i,J,k) = CS%dv_dt(i,J,k) * ADp%diag_hv(i,J,k)
      enddo ; enddo ; enddo
      call post_data(CS%id_h_dv_dt, h_dv_dt, CS%diag)
      deallocate(h_dv_dt)
    endif

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

  if (CS%id_usq > 0) then
    do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      usq(I,j,k) = u(I,j,k) * u(I,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_usq, usq, CS%diag)
  endif

  if (CS%id_vsq > 0) then
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      vsq(i,J,k) = v(i,J,k) * v(i,J,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_vsq, vsq, CS%diag)
  endif

  if (CS%id_uv > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      uv(i,j,k) = (0.5*(u(I-1,j,k) + u(I,j,k))) * &
                (0.5*(v(i,J-1,k) + v(i,J,k)))
    enddo ; enddo ; enddo
    call post_data(CS%id_uv, uv, CS%diag)
  endif

  if (associated(CS%e)) then
    call find_eta(h, tv, G, GV, US, CS%e)
    if (CS%id_e > 0) call post_data(CS%id_e, CS%e, CS%diag)
  endif

  if (associated(CS%e_D)) then
    if (associated(CS%e)) then
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        CS%e_D(i,j,k) = CS%e(i,j,k) + G%bathyT(i,j)
      enddo ; enddo ; enddo
    else
      call find_eta(h, tv, G, GV, US, CS%e_D)
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        CS%e_D(i,j,k) = CS%e_D(i,j,k) + G%bathyT(i,j)
      enddo ; enddo ; enddo
    endif

    if (CS%id_e_D > 0) call post_data(CS%id_e_D, CS%e_D, CS%diag)
  endif

  ! mass per area of grid cell (for Bouss, use Rho0)
  if (CS%id_masscello > 0) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_3d(i,j,k) = GV%H_to_kg_m2*h(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_masscello, work_3d, CS%diag)
  endif

  ! mass of liquid ocean (for Bouss, use Rho0). The reproducing sum requires the use of MKS units.
  if (CS%id_masso > 0) then
    work_2d(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      work_2d(i,j) = work_2d(i,j) + (GV%H_to_kg_m2*h(i,j,k)) * US%L_to_m**2*G%areaT(i,j)
    enddo ; enddo ; enddo
    masso = reproducing_sum(work_2d)
    call post_data(CS%id_masso, masso, CS%diag)
  endif

  ! diagnose thickness/volumes of grid cells [m]
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
          call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pressure_1d, rho_in_situ, &
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
    if ((CS%id_Tpot > 0) .or. (CS%id_tob > 0)) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        work_3d(i,j,k) = gsw_pt_from_ct(tv%S(i,j,k),tv%T(i,j,k))
      enddo ; enddo ; enddo
      if (CS%id_Tpot > 0) call post_data(CS%id_Tpot, work_3d, CS%diag)
      if (CS%id_tob > 0) call post_data(CS%id_tob, work_3d(:,:,nz), CS%diag, mask=G%mask2dT)
    endif
  else
    ! Internal T&S variables are potential temperature & practical salinity
    if (CS%id_tob > 0) call post_data(CS%id_tob, tv%T(:,:,nz), CS%diag, mask=G%mask2dT)
  endif

  ! Calculate additional, potentially derived salinity diagnostics
  if (tv%S_is_absS) then
    ! Internal T&S variables are conservative temperature & absolute salinity,
    ! so they need to converted to potential temperature and practical salinity
    ! for some diagnostics using TEOS-10 function calls.
    if ((CS%id_Sprac > 0) .or. (CS%id_sob > 0)) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        work_3d(i,j,k) = gsw_sp_from_sr(tv%S(i,j,k))
      enddo ; enddo ; enddo
      if (CS%id_Sprac > 0) call post_data(CS%id_Sprac, work_3d, CS%diag)
      if (CS%id_sob > 0) call post_data(CS%id_sob, work_3d(:,:,nz), CS%diag, mask=G%mask2dT)
    endif
  else
    ! Internal T&S variables are potential temperature & practical salinity
    if (CS%id_sob > 0) call post_data(CS%id_sob, tv%S(:,:,nz), CS%diag, mask=G%mask2dT)
  endif

  ! volume mean potential temperature
  if (CS%id_thetaoga>0) then
    thetaoga = global_volume_mean(tv%T, h, G, GV)
    call post_data(CS%id_thetaoga, thetaoga, CS%diag)
  endif

  ! area mean SST
  if (CS%id_tosga > 0) then
    do j=js,je ; do i=is,ie
      surface_field(i,j) = tv%T(i,j,1)
    enddo ; enddo
    tosga = global_area_mean(surface_field, G)
    call post_data(CS%id_tosga, tosga, CS%diag)
  endif

  ! volume mean salinity
  if (CS%id_soga>0) then
    soga = global_volume_mean(tv%S, h, G, GV)
    call post_data(CS%id_soga, soga, CS%diag)
  endif

  ! area mean SSS
  if (CS%id_sosga > 0) then
    do j=js,je ; do i=is,ie
       surface_field(i,j) = tv%S(i,j,1)
    enddo ; enddo
    sosga = global_area_mean(surface_field, G)
    call post_data(CS%id_sosga, sosga, CS%diag)
  endif

  ! layer mean potential temperature
  if (CS%id_temp_layer_ave>0) then
    temp_layer_ave = global_layer_mean(tv%T, h, G, GV)
    call post_data(CS%id_temp_layer_ave, temp_layer_ave, CS%diag)
  endif

  ! layer mean salinity
  if (CS%id_salt_layer_ave>0) then
    salt_layer_ave = global_layer_mean(tv%S, h, G, GV)
    call post_data(CS%id_salt_layer_ave, salt_layer_ave, CS%diag)
  endif

  call calculate_vertical_integrals(h, tv, p_surf, G, GV, US, CS)

  if ((CS%id_Rml > 0) .or. (CS%id_Rcv > 0) .or. associated(CS%h_Rlay) .or. &
      associated(CS%uh_Rlay) .or. associated(CS%vh_Rlay) .or. &
      associated(CS%uhGM_Rlay) .or. associated(CS%vhGM_Rlay)) then

    if (associated(tv%eqn_of_state)) then
      EOSdom(:) = EOS_domain(G%HI, halo=1)
      pressure_1d(:) = tv%P_Ref
      !$OMP parallel do default(shared)
      do k=1,nz ; do j=js-1,je+1
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pressure_1d, Rcv(:,j,k), tv%eqn_of_state, &
                               EOSdom)
      enddo ; enddo
    else ! Rcv should not be used much in this case, so fill in sensible values.
      do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
        Rcv(i,j,k) = GV%Rlay(k)
      enddo ; enddo ; enddo
    endif
    if (CS%id_Rml > 0) call post_data(CS%id_Rml, Rcv, CS%diag)
    if (CS%id_Rcv > 0) call post_data(CS%id_Rcv, Rcv, CS%diag)

    if (associated(CS%h_Rlay)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(is,ie,js,je,nz,nkmb,CS,Rcv,h,GV) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do i=is,ie
          CS%h_Rlay(i,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          CS%h_Rlay(i,j,k) = h(i,j,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, Rcv(i,j,k), k_list, nz, wt, wt_p)
          CS%h_Rlay(i,j,k_list)   = CS%h_Rlay(i,j,k_list)   + h(i,j,k)*wt
          CS%h_Rlay(i,j,k_list+1) = CS%h_Rlay(i,j,k_list+1) + h(i,j,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_h_Rlay > 0) call post_data(CS%id_h_Rlay, CS%h_Rlay, CS%diag)
    endif

    if (associated(CS%uh_Rlay)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(Isq,Ieq,js,je,nz,nkmb,Rcv,CS,GV,uh) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do I=Isq,Ieq
          CS%uh_Rlay(I,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do I=Isq,Ieq
          CS%uh_Rlay(I,j,k) = uh(I,j,k)
        enddo ; enddo
        k_list = nz/2
        do k=1,nkmb ; do I=Isq,Ieq
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i+1,j,k)), k_list, nz, wt, wt_p)
          CS%uh_Rlay(I,j,k_list)   = CS%uh_Rlay(I,j,k_list)   + uh(I,j,k)*wt
          CS%uh_Rlay(I,j,k_list+1) = CS%uh_Rlay(I,j,k_list+1) + uh(I,j,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_uh_Rlay > 0) call post_data(CS%id_uh_Rlay, CS%uh_Rlay, CS%diag)
    endif

    if (associated(CS%vh_Rlay)) then
      k_list = nz/2
!$OMP parallel do default(none)  shared(Jsq,Jeq,is,ie,nz,nkmb,Rcv,CS,GV,vh) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do J=Jsq,Jeq
        do k=1,nkmb ; do i=is,ie
          CS%vh_Rlay(i,J,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          CS%vh_Rlay(i,J,k) = vh(i,J,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i,j+1,k)), k_list, nz, wt, wt_p)
          CS%vh_Rlay(i,J,k_list)   = CS%vh_Rlay(i,J,k_list)   + vh(i,J,k)*wt
          CS%vh_Rlay(i,J,k_list+1) = CS%vh_Rlay(i,J,k_list+1) + vh(i,J,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_vh_Rlay > 0) call post_data(CS%id_vh_Rlay, CS%vh_Rlay, CS%diag)
    endif

    if (associated(CS%uhGM_Rlay) .and. associated(CDp%uhGM)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(Isq,Ieq,js,je,nz,nkmb,Rcv,CDP,CS,GV) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do I=Isq,Ieq
          CS%uhGM_Rlay(I,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do I=Isq,Ieq
          CS%uhGM_Rlay(I,j,k) = CDp%uhGM(I,j,k)
        enddo ; enddo
        do k=1,nkmb ; do I=Isq,Ieq
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i+1,j,k)), k_list, nz, wt, wt_p)
          CS%uhGM_Rlay(I,j,k_list)   = CS%uhGM_Rlay(I,j,k_list)   + CDp%uhGM(I,j,k)*wt
          CS%uhGM_Rlay(I,j,k_list+1) = CS%uhGM_Rlay(I,j,k_list+1) + CDp%uhGM(I,j,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_uhGM_Rlay > 0) call post_data(CS%id_uhGM_Rlay, CS%uhGM_Rlay, CS%diag)
    endif

    if (associated(CS%vhGM_Rlay) .and. associated(CDp%vhGM)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(is,ie,Jsq,Jeq,nz,nkmb,CS,CDp,Rcv,GV) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do J=Jsq,Jeq
        do k=1,nkmb ; do i=is,ie
          CS%vhGM_Rlay(i,J,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          CS%vhGM_Rlay(i,J,k) = CDp%vhGM(i,J,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i,j+1,k)), k_list, nz, wt, wt_p)
          CS%vhGM_Rlay(i,J,k_list)   = CS%vhGM_Rlay(i,J,k_list)   + CDp%vhGM(i,J,k)*wt
          CS%vhGM_Rlay(i,J,k_list+1) = CS%vhGM_Rlay(i,J,k_list+1) + CDp%vhGM(i,J,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_vhGM_Rlay > 0) call post_data(CS%id_vhGM_Rlay, CS%vhGM_Rlay, CS%diag)
    endif
  endif

  if (associated(tv%eqn_of_state)) then
    EOSdom(:) = EOS_domain(G%HI)
    if (CS%id_rhopot0 > 0) then
      pressure_1d(:) = 0.
      !$OMP parallel do default(shared)
      do k=1,nz ; do j=js,je
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pressure_1d, Rcv(:,j,k), &
                                tv%eqn_of_state, EOSdom)
      enddo ; enddo
      if (CS%id_rhopot0 > 0) call post_data(CS%id_rhopot0, Rcv, CS%diag)
    endif
    if (CS%id_rhopot2 > 0) then
      pressure_1d(:) = 2.0e7*US%kg_m3_to_R*US%m_s_to_L_T**2 ! 2000 dbars
      !$OMP parallel do default(shared)
      do k=1,nz ; do j=js,je
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pressure_1d, Rcv(:,j,k), &
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
          call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pressure_1d, Rcv(:,j,k), &
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
          pressure_1d(:) =  pressure_1d(:) + 0.5 * h(:,j,k) * GV%H_to_Pa ! Pressure in middle of layer k
          ! To avoid storing more arrays, put drho_dT into Rcv, and drho_dS into work3d
          call calculate_density_derivs(tv%T(:,j,k),tv%S(:,j,k),pressure_1d, &
                                 Rcv(:,j,k),work_3d(:,j,k),is,ie-is+1, tv%eqn_of_state)
          pressure_1d(:) =  pressure_1d(:) + 0.5 * h(:,j,k) * GV%H_to_Pa ! Pressure at bottom of layer k
        enddo
      enddo
      if (CS%id_drho_dT > 0) call post_data(CS%id_drho_dT, Rcv, CS%diag)
      if (CS%id_drho_dS > 0) call post_data(CS%id_drho_dS, work_3d, CS%diag)
    endif
  endif

  if ((CS%id_cg1>0) .or. (CS%id_Rd1>0) .or. (CS%id_cfl_cg1>0) .or. &
      (CS%id_cfl_cg1_x>0) .or. (CS%id_cfl_cg1_y>0)) then
    call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed_CSp)
    if (CS%id_cg1>0) call post_data(CS%id_cg1, CS%cg1, CS%diag)
    if (CS%id_Rd1>0) then
      !$OMP parallel do default(shared) private(f2_h,mag_beta)
      do j=js,je ; do i=is,ie
        ! Blend the equatorial deformation radius with the standard one.
        f2_h = absurdly_small_freq2 + 0.25 * &
            ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
             (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2))
        mag_beta = sqrt(0.5 * ( &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2) ))
        CS%Rd1(i,j) = CS%cg1(i,j) / sqrt(f2_h + CS%cg1(i,j) * mag_beta)

      enddo ; enddo
      call post_data(CS%id_Rd1, CS%Rd1, CS%diag)
    endif
    if (CS%id_cfl_cg1>0) then
      do j=js,je ; do i=is,ie
        CS%cfl_cg1(i,j) = (dt*CS%cg1(i,j)) * (G%IdxT(i,j) + G%IdyT(i,j))
      enddo ; enddo
      call post_data(CS%id_cfl_cg1, CS%cfl_cg1, CS%diag)
    endif
    if (CS%id_cfl_cg1_x>0) then
      do j=js,je ; do i=is,ie
        CS%cfl_cg1_x(i,j) = (dt*CS%cg1(i,j)) * G%IdxT(i,j)
      enddo ; enddo
      call post_data(CS%id_cfl_cg1_x, CS%cfl_cg1_x, CS%diag)
    endif
    if (CS%id_cfl_cg1_y>0) then
      do j=js,je ; do i=is,ie
        CS%cfl_cg1_y(i,j) = (dt*CS%cg1(i,j)) * G%IdyT(i,j)
      enddo ; enddo
      call post_data(CS%id_cfl_cg1_y, CS%cfl_cg1_y, CS%diag)
    endif
  endif
  if ((CS%id_cg_ebt>0) .or. (CS%id_Rd_ebt>0) .or. (CS%id_p_ebt>0)) then
    if (CS%id_p_ebt>0) then
      call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed_CSp, use_ebt_mode=.true., &
                      mono_N2_column_fraction=CS%mono_N2_column_fraction, &
                      mono_N2_depth=CS%mono_N2_depth, modal_structure=CS%p_ebt)
      call post_data(CS%id_p_ebt, CS%p_ebt, CS%diag)
    else
      call wave_speed(h, tv, G, GV, US, CS%cg1, CS%wave_speed_CSp, use_ebt_mode=.true., &
                      mono_N2_column_fraction=CS%mono_N2_column_fraction, &
                      mono_N2_depth=CS%mono_N2_depth)
    endif
    if (CS%id_cg_ebt>0) call post_data(CS%id_cg_ebt, CS%cg1, CS%diag)
    if (CS%id_Rd_ebt>0) then
      !$OMP parallel do default(shared) private(f2_h,mag_beta)
      do j=js,je ; do i=is,ie
        ! Blend the equatorial deformation radius with the standard one.
        f2_h = absurdly_small_freq2 + 0.25 * &
            ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
             (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2))
        mag_beta = sqrt(0.5 * ( &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2) ))
        CS%Rd1(i,j) = CS%cg1(i,j) / sqrt(f2_h + CS%cg1(i,j) * mag_beta)

      enddo ; enddo
      call post_data(CS%id_Rd_ebt, CS%Rd1, CS%diag)
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
  real,     intent(out)   :: wt   !< The weight of layer k for interpolation, nondim
  real,     intent(out)   :: wt_p !< The weight of layer k+1 for interpolation, nondim

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

  real, dimension(SZI_(G), SZJ_(G)) :: &
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
              ! (Rho_0 in a Boussinesq model) [TR kg m-2].
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
!       do j=js,je ; do i=is,ie ; z_bot(i,j) = -P_SURF(i,j)/GV%H_to_Pa ; enddo ; enddo
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

! This subroutine calculates terms in the mechanical energy budget.

  ! Local variables
  real :: KE_u(SZIB_(G),SZJ_(G))
  real :: KE_v(SZI_(G),SZJB_(G))
  real :: KE_h(SZI_(G),SZJ_(G))

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  do j=js-1,je ; do i=is-1,ie
    KE_u(I,j) = 0.0 ; KE_v(i,J) = 0.0
  enddo ; enddo

  if (associated(CS%KE)) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%KE(i,j,k) = ((u(I,j,k) * u(I,j,k) + u(I-1,j,k) * u(I-1,j,k)) &
          + (v(i,J,k) * v(i,J,k) + v(i,J-1,k) * v(i,J-1,k))) * 0.25
      ! DELETE THE FOLLOWING...  Make this 0 to test the momentum balance,
      ! or a huge number to test the continuity balance.
      ! CS%KE(i,j,k) *= 1e20
    enddo ; enddo ; enddo
    if (CS%id_KE > 0) call post_data(CS%id_KE, CS%KE, CS%diag)
  endif

  if (.not.G%symmetric) then
    if (associated(CS%dKE_dt) .OR. associated(CS%PE_to_KE) .OR. associated(CS%KE_BT) .OR. &
        associated(CS%KE_CorAdv) .OR. associated(CS%KE_adv) .OR. associated(CS%KE_visc) .OR. &
        associated(CS%KE_horvisc) .OR. associated(CS%KE_dia) ) then
      call create_group_pass(CS%pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)
    endif
  endif

  if (associated(CS%dKE_dt)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * CS%du_dt(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * CS%dv_dt(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = CS%KE(i,j,k) * CS%dh_dt(i,j,k)
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%dKE_dt(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) &
            * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_dKEdt > 0) call post_data(CS%id_dKEdt, CS%dKE_dt, CS%diag)
  endif

  if (associated(CS%PE_to_KE)) then
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
        CS%PE_to_KE(i,j,k) = 0.5 * G%IareaT(i,j) &
            * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_PE_to_KE > 0) call post_data(CS%id_PE_to_KE, CS%PE_to_KE, CS%diag)
  endif

  if (associated(CS%KE_BT)) then
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
        CS%KE_BT(i,j,k) = 0.5 * G%IareaT(i,j) &
            * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_BT > 0) call post_data(CS%id_KE_BT, CS%KE_BT, CS%diag)
  endif

  if (associated(CS%KE_CorAdv)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%CAu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%CAv(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = -CS%KE(i,j,k) * G%IareaT(i,j) &
            * (uh(I,j,k) - uh(I-1,j,k) + vh(i,J,k) - vh(i,J-1,k))
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_CorAdv(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) &
            * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_Coradv > 0) call post_data(CS%id_KE_Coradv, CS%KE_Coradv, CS%diag)
  endif

  if (associated(CS%KE_adv)) then
    ! NOTE: All terms in KE_adv are multipled by -1, which can easily produce
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
        KE_h(i,j) = -CS%KE(i,j,k) * G%IareaT(i,j) &
            * (uh(I,j,k) - uh(I-1,j,k) + vh(i,J,k) - vh(i,J-1,k))
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_adv(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) &
            * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_adv > 0) call post_data(CS%id_KE_adv, CS%KE_adv, CS%diag)
  endif

  if (associated(CS%KE_visc)) then
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
        CS%KE_visc(i,j,k) = 0.5 * G%IareaT(i,j) &
            * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_visc > 0) call post_data(CS%id_KE_visc, CS%KE_visc, CS%diag)
  endif

  if (associated(CS%KE_stress)) then
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
        CS%KE_stress(i,j,k) = 0.5 * G%IareaT(i,j) * &
            ((KE_u(I,j) + KE_u(I-1,j)) + (KE_v(i,J) + KE_v(i,J-1)))
      enddo ; enddo
    enddo
    if (CS%id_KE_stress > 0) call post_data(CS%id_KE_stress, CS%KE_stress, CS%diag)
  endif

  if (associated(CS%KE_horvisc)) then
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
        CS%KE_horvisc(i,j,k) = 0.5 * G%IareaT(i,j) &
            * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_horvisc > 0) call post_data(CS%id_KE_horvisc, CS%KE_horvisc, CS%diag)
  endif

  if (associated(CS%KE_dia)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k) * G%dxCu(I,j) * ADp%du_dt_dia(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k) * G%dyCv(i,J) * ADp%dv_dt_dia(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = CS%KE(i,j,k) * (CDp%diapyc_vel(i,j,k) - CDp%diapyc_vel(i,j,k+1))
      enddo ; enddo
      if (.not.G%symmetric) &
        call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_dia(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) &
            * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_dia > 0) call post_data(CS%id_KE_dia, CS%KE_dia, CS%diag)
  endif

end subroutine calculate_energy_diagnostics

!> This subroutine registers fields to calculate a diagnostic time derivative.
subroutine register_time_deriv(lb, f_ptr, deriv_ptr, CS)
  integer, intent(in), dimension(3) :: lb     !< Lower index bound of f_ptr
  real, dimension(lb(1):,lb(2):,:), target :: f_ptr
                                              !< Time derivative operand
  real, dimension(lb(1):,lb(2):,:), target :: deriv_ptr
                                              !< Time derivative of f_ptr
  type(diagnostics_CS),  pointer :: CS        !< Control structure returned by previous call to
                                              !! diagnostics_init.

  ! This subroutine registers fields to calculate a diagnostic time derivative.
  ! NOTE: Lower bound is required for grid indexing in calculate_derivs().
  !       We assume that the vertical axis is 1-indexed.

  integer :: m      !< New index of deriv_ptr in CS%deriv
  integer :: ub(3)  !< Upper index bound of f_ptr, based on shape.

  if (.not.associated(CS)) call MOM_error(FATAL, &
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
                            intent(in) :: ssh !< Time mean surface height without corrections for ice displacement [m]

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: speed  ! The surface speed [L T-1 ~> m s-1]
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
      speed(i,j) = sqrt(0.5*(sfc_state%u(I-1,j)**2 + sfc_state%u(I,j)**2) + &
                        0.5*(sfc_state%v(i,J-1)**2 + sfc_state%v(i,J)**2))
    enddo ; enddo
    call post_data(IDs%id_speed, speed, diag, mask=G%mask2dT)
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
  real, dimension(SZI_(G),SZJ_(G)), &
                            intent(in) :: ssh !< Time mean surface height without corrections for ice displacement [m]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: ssh_ibc !< Time mean surface height with corrections
                                                  !! for ice displacement and the inverse barometer [m]

  real, dimension(SZI_(G),SZJ_(G)) :: work_2d  ! A 2-d work array
  real, dimension(SZI_(G),SZJ_(G)) :: &
    zos  ! dynamic sea lev (zero area mean) from inverse-barometer adjusted ssh [m]
  real :: I_time_int    ! The inverse of the time interval [T-1 ~> s-1].
  real :: zos_area_mean, volo, ssh_ga
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! area mean SSH
  if (IDs%id_ssh_ga > 0) then
    ssh_ga = global_area_mean(ssh, G)
    call post_data(IDs%id_ssh_ga, ssh_ga, diag)
  endif

  ! post the dynamic sea level, zos, and zossq.
  ! zos is ave_ssh with sea ice inverse barometer removed,
  ! and with zero global area mean.
  if (IDs%id_zos > 0 .or. IDs%id_zossq > 0) then
    zos(:,:) = 0.0
    do j=js,je ; do i=is,ie
      zos(i,j) = ssh_ibc(i,j)
    enddo ; enddo
    zos_area_mean = global_area_mean(zos, G)
    do j=js,je ; do i=is,ie
      zos(i,j) = zos(i,j) - G%mask2dT(i,j)*zos_area_mean
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
      work_2d(i,j) = G%mask2dT(i,j)*(ssh(i,j) + US%Z_to_m*G%bathyT(i,j))
    enddo ; enddo
    volo = global_area_integral(work_2d, G)
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
    do j=js,je ; do i=is,ie
      work_2d(i,j) = gsw_pt_from_ct(sfc_state%SSS(i,j), sfc_state%SST(i,j))
    enddo ; enddo
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
    do j=js,je ; do i=is,ie
      work_2d(i,j) = gsw_sp_from_sr(sfc_state%SSS(i,j))
    enddo ; enddo
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
  real, dimension(SZIB_(G), SZJ_(G)) :: umo2d ! Diagnostics of integrated mass transport [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZI_(G), SZJB_(G)) :: vmo2d ! Diagnostics of integrated mass transport [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZIB_(G), SZJ_(G),SZK_(GV)) :: umo ! Diagnostics of layer mass transport [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZI_(G), SZJB_(G),SZK_(GV)) :: vmo ! Diagnostics of layer mass transport [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))   :: h_tend ! Change in layer thickness due to dynamics
                          ! [H s-1 ~> m s-1 or kg m-2 s-1].
  real :: Idt             ! The inverse of the time interval [T-1 ~> s-1]
  real :: H_to_RZ_dt   ! A conversion factor from accumulated transports to fluxes
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
!! that other diagnostis depend upon.
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
  type(diagnostics_CS),       pointer       :: CS   !< Pointer set to point to control structure
                                                    !! for this module.
  type(thermo_var_ptrs),      intent(in)    :: tv   !< A structure pointing to various
                                                    !! thermodynamic variables.

  ! Local variables
  real :: omega, f2_min, convert_H
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_diagnostics" ! This module's name.
  character(len=48) :: thickness_units, flux_units
  real :: wave_speed_min      ! A floor in the first mode speed below which 0 is returned [L T-1 ~> m s-1]
  real :: wave_speed_tol      ! The fractional tolerance for finding the wave speeds [nondim]
  logical :: better_speed_est ! If true, use a more robust estimate of the first
                              ! mode wave speed as the starting point for iterations.
  logical :: split            ! True if using the barotropic-baroclinic split algorithm
  logical :: use_temperature, adiabatic
  logical :: default_2018_answers, remap_answers_2018
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz, nkml, nkbl
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j

  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "MOM_diagnostics_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag
  use_temperature = associated(tv%T)
  call get_param(param_file, mdl, "ADIABATIC", adiabatic, default=.false., &
                 do_not_log=.true.)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DIAG_EBT_MONO_N2_COLUMN_FRACTION", CS%mono_N2_column_fraction, &
                 "The lower fraction of water column over which N2 is limited as monotonic "// &
                 "for the purposes of calculating the equivalent barotropic wave speed.", &
                 units='nondim', default=0.)
  call get_param(param_file, mdl, "DIAG_EBT_MONO_N2_DEPTH", CS%mono_N2_depth, &
                 "The depth below which N2 is limited as monotonic for the "// &
                 "purposes of calculating the equivalent barotropic wave speed.", &
                 units='m', scale=US%m_to_Z, default=-1.)
  call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_TOL", wave_speed_tol, &
                 "The fractional tolerance for finding the wave speeds.", &
                 units="nondim", default=0.001)
  !### Set defaults so that wave_speed_min*wave_speed_tol >= 1e-9 m s-1
  call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_MIN", wave_speed_min, &
                 "A floor in the first mode speed below which 0 used instead.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "INTERNAL_WAVE_SPEED_BETTER_EST", better_speed_est, &
                 "If true, use a more robust estimate of the first mode wave speed as the "//&
                 "starting point for iterations.", default=.false.) !### Change the default.
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", remap_answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call get_param(param_file, mdl, "SPLIT", split, default=.true., do_not_log=.true.)

  if (GV%Boussinesq) then
    thickness_units = "m" ; flux_units = "m3 s-1" ; convert_H = GV%H_to_m
  else
    thickness_units = "kg m-2" ; flux_units = "kg s-1" ; convert_H = GV%H_to_kg_m2
  endif

  CS%id_masscello = register_diag_field('ocean_model', 'masscello', diag%axesTL,&
      Time, 'Mass per unit area of liquid ocean grid cell', 'kg m-2',           &
      standard_name='sea_water_mass_per_unit_area', v_extensive=.true.)

  CS%id_masso = register_scalar_field('ocean_model', 'masso', Time,  &
      diag, 'Mass of liquid ocean', 'kg', standard_name='sea_water_mass')

  CS%id_thkcello = register_diag_field('ocean_model', 'thkcello', diag%axesTL, Time, &
      long_name = 'Cell Thickness', standard_name='cell_thickness', &
      units='m', conversion=US%Z_to_m, v_extensive=.true.)
  CS%id_h_pre_sync = register_diag_field('ocean_model', 'h_pre_sync', diag%axesTL, Time, &
      long_name = 'Cell thickness from the previous timestep', &
      units='m', conversion=GV%H_to_m, v_extensive=.true.)

  ! Note that CS%id_volcello would normally be registered here but because it is a "cell measure" and
  ! must be registered first. We earlier stored the handle of volcello but need it here for posting
  ! by this module.
  CS%id_volcello = diag_get_volume_cell_measure_dm_id(diag)

  if (use_temperature) then
    if (tv%T_is_conT) then
      CS%id_Tpot = register_diag_field('ocean_model', 'temp', diag%axesTL, &
          Time, 'Potential Temperature', 'degC')
    endif
    if (tv%S_is_absS) then
      CS%id_Sprac = register_diag_field('ocean_model', 'salt', diag%axesTL, &
          Time, 'Salinity', 'psu')
    endif

    CS%id_tob = register_diag_field('ocean_model','tob', diag%axesT1, Time, &
        long_name='Sea Water Potential Temperature at Sea Floor', &
        standard_name='sea_water_potential_temperature_at_sea_floor', units='degC')
    CS%id_sob = register_diag_field('ocean_model','sob',diag%axesT1, Time, &
        long_name='Sea Water Salinity at Sea Floor', &
        standard_name='sea_water_salinity_at_sea_floor', units='psu')

    CS%id_temp_layer_ave = register_diag_field('ocean_model', 'temp_layer_ave', &
        diag%axesZL, Time, 'Layer Average Ocean Temperature', 'degC')
    CS%id_salt_layer_ave = register_diag_field('ocean_model', 'salt_layer_ave', &
        diag%axesZL, Time, 'Layer Average Ocean Salinity', 'psu')

    CS%id_thetaoga = register_scalar_field('ocean_model', 'thetaoga', &
        Time, diag, 'Global Mean Ocean Potential Temperature', 'degC',&
        standard_name='sea_water_potential_temperature')
    CS%id_soga = register_scalar_field('ocean_model', 'soga', &
        Time, diag, 'Global Mean Ocean Salinity', 'psu', &
        standard_name='sea_water_salinity')

    CS%id_tosga = register_scalar_field('ocean_model', 'sst_global', Time, diag,&
        long_name='Global Area Average Sea Surface Temperature',                &
        units='degC', standard_name='sea_surface_temperature',                  &
        cmor_field_name='tosga', cmor_standard_name='sea_surface_temperature',  &
        cmor_long_name='Sea Surface Temperature')
    CS%id_sosga = register_scalar_field('ocean_model', 'sss_global', Time, diag,&
        long_name='Global Area Average Sea Surface Salinity',                   &
        units='psu', standard_name='sea_surface_salinity',                      &
        cmor_field_name='sosga', cmor_standard_name='sea_surface_salinity',     &
        cmor_long_name='Sea Surface Salinity')
  endif

  CS%id_u = register_diag_field('ocean_model', 'u', diag%axesCuL, Time,              &
      'Zonal velocity', 'm s-1', conversion=US%L_T_to_m_s, cmor_field_name='uo', &
      cmor_standard_name='sea_water_x_velocity', cmor_long_name='Sea Water X Velocity')
  CS%id_v = register_diag_field('ocean_model', 'v', diag%axesCvL, Time,                  &
      'Meridional velocity', 'm s-1', conversion=US%L_T_to_m_s, cmor_field_name='vo', &
      cmor_standard_name='sea_water_y_velocity', cmor_long_name='Sea Water Y Velocity')
  CS%id_usq = register_diag_field('ocean_model', 'usq', diag%axesCuL, Time,              &
      'Zonal velocity squared', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_vsq = register_diag_field('ocean_model', 'vsq', diag%axesCvL, Time,                  &
      'Meridional velocity squared', 'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_uv = register_diag_field('ocean_model', 'uv', diag%axesTL, Time, &
      'Product between zonal and meridional velocities at h-points', &
      'm2 s-2', conversion=US%L_T_to_m_s**2)
  CS%id_h = register_diag_field('ocean_model', 'h', diag%axesTL, Time, &
      'Layer Thickness', thickness_units, v_extensive=.true., conversion=convert_H)

  CS%id_e = register_diag_field('ocean_model', 'e', diag%axesTi, Time, &
      'Interface Height Relative to Mean Sea Level', 'm', conversion=US%Z_to_m)
  if (CS%id_e>0) call safe_alloc_ptr(CS%e,isd,ied,jsd,jed,nz+1)

  CS%id_e_D = register_diag_field('ocean_model', 'e_D', diag%axesTi, Time, &
      'Interface Height above the Seafloor', 'm', conversion=US%Z_to_m)
  if (CS%id_e_D>0) call safe_alloc_ptr(CS%e_D,isd,ied,jsd,jed,nz+1)

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
      'Partial derivative of rhoinsitu with respect to temperature (alpha)', 'kg m-3 degC-1')
  CS%id_drho_dS = register_diag_field('ocean_model', 'drho_dS', diag%axesTL, Time, &
      'Partial derivative of rhoinsitu with respect to salinity (beta)', 'kg^2 g-1 m-3')

  CS%id_du_dt = register_diag_field('ocean_model', 'dudt', diag%axesCuL, Time, &
      'Zonal Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_du_dt>0) .and. .not.associated(CS%du_dt)) then
    call safe_alloc_ptr(CS%du_dt,IsdB,IedB,jsd,jed,nz)
    call register_time_deriv(lbound(MIS%u), MIS%u, CS%du_dt, CS)
  endif

  CS%id_dv_dt = register_diag_field('ocean_model', 'dvdt', diag%axesCvL, Time, &
      'Meridional Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_dv_dt>0) .and. .not.associated(CS%dv_dt)) then
    call safe_alloc_ptr(CS%dv_dt,isd,ied,JsdB,JedB,nz)
    call register_time_deriv(lbound(MIS%v), MIS%v, CS%dv_dt, CS)
  endif

  CS%id_dh_dt = register_diag_field('ocean_model', 'dhdt', diag%axesTL, Time, &
      'Thickness tendency', trim(thickness_units)//" s-1", conversion=convert_H*US%s_to_T, v_extensive=.true.)
  if ((CS%id_dh_dt>0) .and. .not.associated(CS%dh_dt)) then
    call safe_alloc_ptr(CS%dh_dt,isd,ied,jsd,jed,nz)
    call register_time_deriv(lbound(MIS%h), MIS%h, CS%dh_dt, CS)
  endif

  !CS%id_hf_du_dt = register_diag_field('ocean_model', 'hf_dudt', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Zonal Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2, &
  !    v_extensive=.true.)
  !if (CS%id_hf_du_dt > 0) then
  !  call safe_alloc_ptr(CS%hf_du_dt,IsdB,IedB,jsd,jed,nz)
  !  if (.not.associated(CS%du_dt)) then
  !    call safe_alloc_ptr(CS%du_dt,IsdB,IedB,jsd,jed,nz)
  !    call register_time_deriv(lbound(MIS%u), MIS%u, CS%du_dt, CS)
  !  endif
  !  call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  !endif

  !CS%id_hf_dv_dt = register_diag_field('ocean_model', 'hf_dvdt', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2, &
  !    v_extensive=.true.)
  !if (CS%id_hf_dv_dt > 0) then
  !  call safe_alloc_ptr(CS%hf_dv_dt,isd,ied,JsdB,JedB,nz)
  !  if (.not.associated(CS%dv_dt)) then
  !    call safe_alloc_ptr(CS%dv_dt,isd,ied,JsdB,JedB,nz)
  !    call register_time_deriv(lbound(MIS%v), MIS%v, CS%dv_dt, CS)
  !  endif
  !  call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)
  !endif

  CS%id_hf_du_dt_2d = register_diag_field('ocean_model', 'hf_dudt_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_du_dt_2d > 0) then
    if (.not.associated(CS%du_dt)) then
      call safe_alloc_ptr(CS%du_dt,IsdB,IedB,jsd,jed,nz)
      call register_time_deriv(lbound(MIS%u), MIS%u, CS%du_dt, CS)
    endif
    call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_hf_dv_dt_2d = register_diag_field('ocean_model', 'hf_dvdt_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Acceleration', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_dv_dt_2d > 0) then
    if (.not.associated(CS%dv_dt)) then
      call safe_alloc_ptr(CS%dv_dt,isd,ied,JsdB,JedB,nz)
      call register_time_deriv(lbound(MIS%v), MIS%v, CS%dv_dt, CS)
    endif
    call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)
  endif

  CS%id_h_du_dt = register_diag_field('ocean_model', 'h_du_dt', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Acceleration', 'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_du_dt > 0) then
    if (.not.associated(CS%du_dt)) then
      call safe_alloc_ptr(CS%du_dt,IsdB,IedB,jsd,jed,nz)
      call register_time_deriv(lbound(MIS%u), MIS%u, CS%du_dt, CS)
    endif
    call safe_alloc_ptr(ADp%diag_hu,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_h_dv_dt = register_diag_field('ocean_model', 'h_dv_dt', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Acceleration', 'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_dv_dt > 0) then
    if (.not.associated(CS%dv_dt)) then
      call safe_alloc_ptr(CS%dv_dt,isd,ied,JsdB,JedB,nz)
      call register_time_deriv(lbound(MIS%v), MIS%v, CS%dv_dt, CS)
    endif
    call safe_alloc_ptr(ADp%diag_hv,isd,ied,JsdB,JedB,nz)
  endif

  ! layer thickness variables
  !if (GV%nk_rho_varies > 0) then
    CS%id_h_Rlay = register_diag_field('ocean_model', 'h_rho', diag%axesTL, Time, &
        'Layer thicknesses in pure potential density coordinates', &
        thickness_units, conversion=convert_H)
    if (CS%id_h_Rlay>0) call safe_alloc_ptr(CS%h_Rlay,isd,ied,jsd,jed,nz)

    CS%id_uh_Rlay = register_diag_field('ocean_model', 'uh_rho', diag%axesCuL, Time, &
        'Zonal volume transport in pure potential density coordinates', &
        flux_units, conversion=US%L_to_m**2*US%s_to_T*convert_H)
    if (CS%id_uh_Rlay>0) call safe_alloc_ptr(CS%uh_Rlay,IsdB,IedB,jsd,jed,nz)

    CS%id_vh_Rlay = register_diag_field('ocean_model', 'vh_rho', diag%axesCvL, Time, &
        'Meridional volume transport in pure potential density coordinates', &
        flux_units, conversion=US%L_to_m**2*US%s_to_T*convert_H)
    if (CS%id_vh_Rlay>0) call safe_alloc_ptr(CS%vh_Rlay,isd,ied,JsdB,JedB,nz)

    CS%id_uhGM_Rlay = register_diag_field('ocean_model', 'uhGM_rho', diag%axesCuL, Time, &
        'Zonal volume transport due to interface height diffusion in pure potential '//&
        'density coordinates', flux_units, conversion=US%L_to_m**2*US%s_to_T*convert_H)
    if (CS%id_uhGM_Rlay>0) call safe_alloc_ptr(CS%uhGM_Rlay,IsdB,IedB,jsd,jed,nz)

    CS%id_vhGM_Rlay = register_diag_field('ocean_model', 'vhGM_rho', diag%axesCvL, Time, &
        'Meridional volume transport due to interface height diffusion in pure potential '//&
        'density coordinates', flux_units, conversion=US%L_to_m**2*US%s_to_T*convert_H)
    if (CS%id_vhGM_Rlay>0) call safe_alloc_ptr(CS%vhGM_Rlay,isd,ied,JsdB,JedB,nz)
  !endif


  ! terms in the kinetic energy budget
  CS%id_KE = register_diag_field('ocean_model', 'KE', diag%axesTL, Time, &
      'Layer kinetic energy per unit mass', &
      'm2 s-2', conversion=US%L_T_to_m_s**2)
  if (CS%id_KE>0) call safe_alloc_ptr(CS%KE,isd,ied,jsd,jed,nz)

  CS%id_dKEdt = register_diag_field('ocean_model', 'dKE_dt', diag%axesTL, Time, &
      'Kinetic Energy Tendency of Layer', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (CS%id_dKEdt>0) call safe_alloc_ptr(CS%dKE_dt,isd,ied,jsd,jed,nz)

  CS%id_PE_to_KE = register_diag_field('ocean_model', 'PE_to_KE', diag%axesTL, Time, &
      'Potential to Kinetic Energy Conversion of Layer', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (CS%id_PE_to_KE>0) call safe_alloc_ptr(CS%PE_to_KE,isd,ied,jsd,jed,nz)

  if (split) then
    CS%id_KE_BT = register_diag_field('ocean_model', 'KE_BT', diag%axesTL, Time, &
        'Barotropic contribution to Kinetic Energy', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
    if (CS%id_KE_BT>0) call safe_alloc_ptr(CS%KE_BT,isd,ied,jsd,jed,nz)
  endif

  CS%id_KE_Coradv = register_diag_field('ocean_model', 'KE_Coradv', diag%axesTL, Time, &
      'Kinetic Energy Source from Coriolis and Advection', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (CS%id_KE_Coradv>0) call safe_alloc_ptr(CS%KE_Coradv,isd,ied,jsd,jed,nz)

  CS%id_KE_adv = register_diag_field('ocean_model', 'KE_adv', diag%axesTL, Time, &
      'Kinetic Energy Source from Advection', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (CS%id_KE_adv>0) call safe_alloc_ptr(CS%KE_adv,isd,ied,jsd,jed,nz)

  CS%id_KE_visc = register_diag_field('ocean_model', 'KE_visc', diag%axesTL, Time, &
      'Kinetic Energy Source from Vertical Viscosity and Stresses', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (CS%id_KE_visc>0) call safe_alloc_ptr(CS%KE_visc,isd,ied,jsd,jed,nz)

  CS%id_KE_stress = register_diag_field('ocean_model', 'KE_stress', diag%axesTL, Time, &
      'Kinetic Energy Source from Surface Stresses or Body Wind Stress', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (CS%id_KE_stress>0) call safe_alloc_ptr(CS%KE_stress,isd,ied,jsd,jed,nz)

  CS%id_KE_horvisc = register_diag_field('ocean_model', 'KE_horvisc', diag%axesTL, Time, &
      'Kinetic Energy Source from Horizontal Viscosity', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  if (CS%id_KE_horvisc>0) call safe_alloc_ptr(CS%KE_horvisc,isd,ied,jsd,jed,nz)

  if (.not. adiabatic) then
    CS%id_KE_dia = register_diag_field('ocean_model', 'KE_dia', diag%axesTL, Time, &
        'Kinetic Energy Source from Diapycnal Diffusion', &
        'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
    if (CS%id_KE_dia>0) call safe_alloc_ptr(CS%KE_dia,isd,ied,jsd,jed,nz)
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
    call wave_speed_init(CS%wave_speed_CSp, remap_answers_2018=remap_answers_2018, &
                         better_speed_est=better_speed_est, min_speed=wave_speed_min, &
                         wave_speed_tol=wave_speed_tol)
!###    call wave_speed_init(CS%wave_speed_CSp, remap_answers_2018=remap_answers_2018)
    call safe_alloc_ptr(CS%cg1,isd,ied,jsd,jed)
    if (CS%id_Rd1>0)       call safe_alloc_ptr(CS%Rd1,isd,ied,jsd,jed)
    if (CS%id_Rd_ebt>0)    call safe_alloc_ptr(CS%Rd1,isd,ied,jsd,jed)
    if (CS%id_cfl_cg1>0)   call safe_alloc_ptr(CS%cfl_cg1,isd,ied,jsd,jed)
    if (CS%id_cfl_cg1_x>0) call safe_alloc_ptr(CS%cfl_cg1_x,isd,ied,jsd,jed)
    if (CS%id_cfl_cg1_y>0) call safe_alloc_ptr(CS%cfl_cg1_y,isd,ied,jsd,jed)
    if (CS%id_p_ebt>0) call safe_alloc_ptr(CS%p_ebt,isd,ied,jsd,jed,nz)
  endif

  CS%id_mass_wt = register_diag_field('ocean_model', 'mass_wt', diag%axesT1, Time,  &
      'The column mass for calculating mass-weighted average properties', 'kg m-2', conversion=US%RZ_to_kg_m2)

  if (use_temperature) then
    CS%id_temp_int = register_diag_field('ocean_model', 'temp_int', diag%axesT1, Time,                &
        'Density weighted column integrated potential temperature', 'degC kg m-2', conversion=US%RZ_to_kg_m2, &
        cmor_field_name='opottempmint',                                                               &
        cmor_long_name='integral_wrt_depth_of_product_of_sea_water_density_and_potential_temperature',&
        cmor_standard_name='Depth integrated density times potential temperature')

    CS%id_salt_int = register_diag_field('ocean_model', 'salt_int', diag%axesT1, Time,   &
        'Density weighted column integrated salinity', 'psu kg m-2', conversion=US%RZ_to_kg_m2, &
        cmor_field_name='somint',                                                        &
        cmor_long_name='integral_wrt_depth_of_product_of_sea_water_density_and_salinity',&
        cmor_standard_name='Depth integrated density times salinity')
  endif

  CS%id_col_mass = register_diag_field('ocean_model', 'col_mass', diag%axesT1, Time, &
      'The column integrated in situ density', 'kg m-2', conversion=US%RZ_to_kg_m2)

  CS%id_col_ht = register_diag_field('ocean_model', 'col_height', diag%axesT1, Time, &
      'The height of the water column', 'm', conversion=US%Z_to_m)
  CS%id_pbo = register_diag_field('ocean_model', 'pbo', diag%axesT1, Time, &
      long_name='Sea Water Pressure at Sea Floor', standard_name='sea_water_pressure_at_sea_floor', &
      units='Pa', conversion=US%RL2_T2_to_Pa)

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
  IDs%id_volo = register_scalar_field('ocean_model', 'volo', Time, diag,&
      long_name='Total volume of liquid ocean', units='m3',            &
      standard_name='sea_water_volume')
  IDs%id_zos = register_diag_field('ocean_model', 'zos', diag%axesT1, Time,&
      standard_name = 'sea_surface_height_above_geoid',                   &
      long_name= 'Sea surface height above geoid', units='m')
  IDs%id_zossq = register_diag_field('ocean_model', 'zossq', diag%axesT1, Time,&
      standard_name='square_of_sea_surface_height_above_geoid',             &
      long_name='Square of sea surface height above geoid', units='m2')
  IDs%id_ssh = register_diag_field('ocean_model', 'SSH', diag%axesT1, Time, &
      'Sea Surface Height', 'm')
  IDs%id_ssh_ga = register_scalar_field('ocean_model', 'ssh_ga', Time, diag,&
      long_name='Area averaged sea surface height', units='m',            &
      standard_name='area_averaged_sea_surface_height')
  IDs%id_ssu = register_diag_field('ocean_model', 'SSU', diag%axesCu1, Time, &
      'Sea Surface Zonal Velocity', 'm s-1', conversion=US%L_T_to_m_s)
  IDs%id_ssv = register_diag_field('ocean_model', 'SSV', diag%axesCv1, Time, &
      'Sea Surface Meridional Velocity', 'm s-1', conversion=US%L_T_to_m_s)
  IDs%id_speed = register_diag_field('ocean_model', 'speed', diag%axesT1, Time, &
      'Sea Surface Speed', 'm s-1', conversion=US%L_T_to_m_s)

  if (associated(tv%T)) then
    IDs%id_sst = register_diag_field('ocean_model', 'SST', diag%axesT1, Time,     &
        'Sea Surface Temperature', 'degC', cmor_field_name='tos', &
        cmor_long_name='Sea Surface Temperature',                                &
        cmor_standard_name='sea_surface_temperature')
    IDs%id_sst_sq = register_diag_field('ocean_model', 'SST_sq', diag%axesT1, Time, &
        'Sea Surface Temperature Squared', 'degC2', cmor_field_name='tossq', &
        cmor_long_name='Square of Sea Surface Temperature ',                      &
        cmor_standard_name='square_of_sea_surface_temperature')
    IDs%id_sss = register_diag_field('ocean_model', 'SSS', diag%axesT1, Time, &
        'Sea Surface Salinity', 'psu', cmor_field_name='sos', &
        cmor_long_name='Sea Surface Salinity',                            &
        cmor_standard_name='sea_surface_salinity')
    IDs%id_sss_sq = register_diag_field('ocean_model', 'SSS_sq', diag%axesT1, Time, &
        'Sea Surface Salinity Squared', 'psu', cmor_field_name='sossq', &
        cmor_long_name='Square of Sea Surface Salinity ',                     &
        cmor_standard_name='square_of_sea_surface_salinity')
    if (tv%T_is_conT) then
      IDs%id_sstcon = register_diag_field('ocean_model', 'conSST', diag%axesT1, Time,     &
          'Sea Surface Conservative Temperature', 'Celsius')
    endif
    if (tv%S_is_absS) then
      IDs%id_sssabs = register_diag_field('ocean_model', 'absSSS', diag%axesT1, Time,     &
          'Sea Surface Absolute Salinity', 'g kg-1')
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
         'ppt kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s)
  IDs%id_Heat_PmE = register_diag_field('ocean_model', 'Heat_PmE', diag%axesT1, Time, &
         'Heat flux into ocean from mass flux into ocean', &
         'W m-2', conversion=US%QRZ_T_to_W_m2)
  IDs%id_intern_heat = register_diag_field('ocean_model', 'internal_heat', diag%axesT1, Time,&
         'Heat flux into ocean from geothermal or other internal sources', &
         'W m-2', conversion=US%QRZ_T_to_W_m2)

end subroutine register_surface_diags

!> Register certain diagnostics related to transports
subroutine register_transport_diags(Time, G, GV, US, IDs, diag)
  type(time_type),          intent(in)    :: Time  !< current model time
  type(ocean_grid_type),    intent(in)    :: G     !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV    !< ocean vertical grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(transport_diag_IDs), intent(inout) :: IDs   !< A structure with the diagnostic IDs.
  type(diag_ctrl),          intent(inout) :: diag  !< regulates diagnostic output

  real :: H_convert
  character(len=48) :: thickness_units, accum_flux_units

  thickness_units = get_thickness_units(GV)
  if (GV%Boussinesq) then
    H_convert = GV%H_to_m ; accum_flux_units = "m3"
  else
    H_convert = GV%H_to_kg_m2 ; accum_flux_units = "kg"
  endif

  ! Diagnostics related to tracer and mass transport
  IDs%id_uhtr = register_diag_field('ocean_model', 'uhtr', diag%axesCuL, Time, &
      'Accumulated zonal thickness fluxes to advect tracers', &
      accum_flux_units, y_cell_method='sum', v_extensive=.true., conversion=H_convert*US%L_to_m**2)
  IDs%id_vhtr = register_diag_field('ocean_model', 'vhtr', diag%axesCvL, Time, &
      'Accumulated meridional thickness fluxes to advect tracers', &
      accum_flux_units, x_cell_method='sum', v_extensive=.true., conversion=H_convert*US%L_to_m**2)
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
  IDs%id_dynamics_h = register_diag_field('ocean_model','dynamics_h',  &
      diag%axesTl, Time, 'Layer thicknesses prior to horizontal dynamics', &
      'm', v_extensive=.true., conversion=GV%H_to_m)
  IDs%id_dynamics_h_tendency = register_diag_field('ocean_model','dynamics_h_tendency',  &
      diag%axesTl, Time, 'Change in layer thicknesses due to horizontal dynamics', &
      'm s-1', v_extensive=.true., conversion=GV%H_to_m*US%s_to_T)

end subroutine register_transport_diags

!> Offers the static fields in the ocean grid type for output via the diag_manager.
subroutine write_static_fields(G, GV, US, tv, diag)
  type(ocean_grid_type),   intent(in)    :: G    !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in)    :: tv   !< A structure pointing to various thermodynamic variables
  type(diag_ctrl), target, intent(inout) :: diag !< regulates diagnostic output

  ! Local variables
  integer :: id
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

  id = register_static_field('ocean_model', 'depth_ocean', diag%axesT1,  &
        'Depth of the ocean at tracer points', 'm', conversion=US%Z_to_m, &
        standard_name='sea_floor_depth_below_geoid',                     &
        cmor_field_name='deptho', cmor_long_name='Sea Floor Depth',      &
        cmor_standard_name='sea_floor_depth_below_geoid', area=diag%axesT1%id_area, &
        x_cell_method='mean', y_cell_method='mean', area_cell_method='mean')
  if (id > 0) call post_data(id, G%bathyT, diag, .true., mask=G%mask2dT)

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
         'heat capacity of sea water', 'J kg-1 K-1', conversion=US%Q_to_J_kg, &
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
  type(diagnostics_CS),       pointer       :: CS  !< Pointer to the control structure for this
                                                   !! module.

! This subroutine sets up diagnostics upon which other diagnostics depend.
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS%dKE_dt) .or. associated(CS%PE_to_KE) .or. &
      associated(CS%KE_BT) .or. associated(CS%KE_CorAdv) .or. &
      associated(CS%KE_adv) .or. associated(CS%KE_visc) .or. associated(CS%KE_stress) .or. &
      associated(CS%KE_horvisc) .or. associated(CS%KE_dia)) &
    call safe_alloc_ptr(CS%KE,isd,ied,jsd,jed,nz)

  if (associated(CS%dKE_dt)) then
    if (.not.associated(CS%du_dt)) then
      call safe_alloc_ptr(CS%du_dt,IsdB,IedB,jsd,jed,nz)
      call register_time_deriv(lbound(MIS%u), MIS%u, CS%du_dt, CS)
    endif
    if (.not.associated(CS%dv_dt)) then
      call safe_alloc_ptr(CS%dv_dt,isd,ied,JsdB,JedB,nz)
      call register_time_deriv(lbound(MIS%v), MIS%v, CS%dv_dt, CS)
    endif
    if (.not.associated(CS%dh_dt)) then
      call safe_alloc_ptr(CS%dh_dt,isd,ied,jsd,jed,nz)
      call register_time_deriv(lbound(MIS%h), MIS%h, CS%dh_dt, CS)
    endif
  endif

  if (associated(CS%KE_adv)) then
    call safe_alloc_ptr(ADp%gradKEu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%gradKEv,isd,ied,JsdB,JedB,nz)
  endif

  if (associated(CS%KE_visc)) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  endif

  if (associated(CS%KE_stress)) then
    call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)
  endif

  if (associated(CS%KE_dia)) then
    call safe_alloc_ptr(ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)
  endif

  if (associated(CS%uhGM_Rlay)) call safe_alloc_ptr(CDp%uhGM,IsdB,IedB,jsd,jed,nz)
  if (associated(CS%vhGM_Rlay)) call safe_alloc_ptr(CDp%vhGM,isd,ied,JsdB,JedB,nz)

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

  if (associated(CS%e))          deallocate(CS%e)
  if (associated(CS%e_D))        deallocate(CS%e_D)
  if (associated(CS%KE))         deallocate(CS%KE)
  if (associated(CS%dKE_dt))     deallocate(CS%dKE_dt)
  if (associated(CS%PE_to_KE))   deallocate(CS%PE_to_KE)
  if (associated(CS%KE_BT))      deallocate(CS%KE_BT)
  if (associated(CS%KE_Coradv))  deallocate(CS%KE_Coradv)
  if (associated(CS%KE_adv))     deallocate(CS%KE_adv)
  if (associated(CS%KE_visc))    deallocate(CS%KE_visc)
  if (associated(CS%KE_stress))  deallocate(CS%KE_stress)
  if (associated(CS%KE_horvisc)) deallocate(CS%KE_horvisc)
  if (associated(CS%KE_dia))     deallocate(CS%KE_dia)
  if (associated(CS%dv_dt))      deallocate(CS%dv_dt)
  if (associated(CS%dh_dt))      deallocate(CS%dh_dt)
  if (associated(CS%du_dt))      deallocate(CS%du_dt)
  if (associated(CS%h_Rlay))     deallocate(CS%h_Rlay)
  if (associated(CS%uh_Rlay))    deallocate(CS%uh_Rlay)
  if (associated(CS%vh_Rlay))    deallocate(CS%vh_Rlay)
  if (associated(CS%uhGM_Rlay))  deallocate(CS%uhGM_Rlay)
  if (associated(CS%vhGM_Rlay))  deallocate(CS%vhGM_Rlay)

  if (associated(ADp%gradKEu))    deallocate(ADp%gradKEu)
  if (associated(ADp%gradKEu))    deallocate(ADp%gradKEu)
  if (associated(ADp%du_dt_visc)) deallocate(ADp%du_dt_visc)
  if (associated(ADp%dv_dt_visc)) deallocate(ADp%dv_dt_visc)
  if (associated(ADp%du_dt_str))  deallocate(ADp%du_dt_str)
  if (associated(ADp%dv_dt_str))  deallocate(ADp%dv_dt_str)
  if (associated(ADp%du_dt_dia))  deallocate(ADp%du_dt_dia)
  if (associated(ADp%dv_dt_dia))  deallocate(ADp%dv_dt_dia)
  if (associated(ADp%du_other))   deallocate(ADp%du_other)
  if (associated(ADp%dv_other))   deallocate(ADp%dv_other)

  if (associated(ADp%diag_hfrac_u)) deallocate(ADp%diag_hfrac_u)
  if (associated(ADp%diag_hfrac_v)) deallocate(ADp%diag_hfrac_v)

  ! NOTE: [uv]hGM may be allocated either here or the thickness diffuse module
  if (associated(CDp%uhGM)) deallocate(CDp%uhGM)
  if (associated(CDp%vhGM)) deallocate(CDp%vhGM)
  if (associated(CDp%diapyc_vel)) deallocate(CDp%diapyc_vel)

  do m=1,CS%num_time_deriv ; deallocate(CS%prev_val(m)%p) ; enddo
end subroutine MOM_diagnostics_end

end module MOM_diagnostics
