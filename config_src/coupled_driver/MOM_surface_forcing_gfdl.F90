module MOM_surface_forcing_gfdl

! This file is part of MOM6. See LICENSE.md for the license.

!#CTRL# use MOM_controlled_forcing, only : apply_ctrl_forcing, register_ctrl_forcing_restarts
!#CTRL# use MOM_controlled_forcing, only : controlled_forcing_init, controlled_forcing_end
!#CTRL# use MOM_controlled_forcing, only : ctrl_forcing_CS
use MOM_coms,             only : reproducing_sum
use MOM_constants,        only : hlv, hlf
use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_SUBCOMPONENT
use MOM_diag_mediator,    only : diag_ctrl
use MOM_diag_mediator,    only : safe_alloc_ptr, time_type
use MOM_domains,          only : pass_vector, pass_var, fill_symmetric_edges
use MOM_domains,          only : global_field_sum, BITWISE_EXACT_SUM
use MOM_domains,          only : AGRID, BGRID_NE, CGRID_NE, To_All
use MOM_domains,          only : To_North, To_East, Omit_Corners
use MOM_error_handler,    only : MOM_error, WARNING, FATAL, is_root_pe, MOM_mesg
use MOM_file_parser,      only : get_param, log_version, param_file_type
use MOM_forcing_type,     only : forcing, mech_forcing
use MOM_forcing_type,     only : forcing_diags, mech_forcing_diags, register_forcing_type_diags
use MOM_forcing_type,     only : allocate_forcing_type, deallocate_forcing_type
use MOM_forcing_type,     only : allocate_mech_forcing, deallocate_mech_forcing
use MOM_get_input,        only : Get_MOM_Input, directories
use MOM_grid,             only : ocean_grid_type
use MOM_io,               only : slasher, write_version_number, MOM_read_data
use MOM_restart,          only : register_restart_field, restart_init, MOM_restart_CS
use MOM_restart,          only : restart_init_end, save_restart, restore_state
use MOM_string_functions, only : uppercase
use MOM_spatial_means,    only : adjust_area_mean_to_zero
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : surface
use user_revise_forcing,  only : user_alter_forcing, user_revise_forcing_init
use user_revise_forcing,  only : user_revise_forcing_CS

use coupler_types_mod,    only : coupler_2d_bc_type, coupler_type_write_chksums
use coupler_types_mod,    only : coupler_type_initialized, coupler_type_spawn
use coupler_types_mod,    only : coupler_type_copy_data
use data_override_mod,    only : data_override_init, data_override
use fms_mod,              only : stdout
use mpp_mod,              only : mpp_chksum
use time_interp_external_mod, only : init_external_field, time_interp_external
use time_interp_external_mod, only : time_interp_external_init

implicit none ; private

#include <MOM_memory.h>

public convert_IOB_to_fluxes, convert_IOB_to_forces
public surface_forcing_init
public ice_ocn_bnd_type_chksum
public forcing_save_restart

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> surface_forcing_CS is a structure containing pointers to the forcing fields
!! which may be used to drive MOM.  All fluxes are positive downward.
type, public :: surface_forcing_CS ; private
  integer :: wind_stagger       !< AGRID, BGRID_NE, or CGRID_NE (integer values
                                !! from MOM_domains) to indicate the staggering of
                                !! the winds that are being provided in calls to
                                !! update_ocean_model.
  logical :: use_temperature    !< If true, temp and saln used as state variables
  real :: wind_stress_multiplier !< A multiplier applied to incoming wind stress [nondim].

  real :: Rho0                  !< Boussinesq reference density [R ~> kg m-3]
  real :: area_surf = -1.0      !< Total ocean surface area [m2]
  real :: latent_heat_fusion    !< Latent heat of fusion [J kg-1]
  real :: latent_heat_vapor     !< Latent heat of vaporization [J kg-1]

  real :: max_p_surf            !< The maximum surface pressure that can be
                                !! exerted by the atmosphere and floating sea-ice [Pa].
                                !! This is needed because the FMS coupling structure
                                !! does not limit the water that can be frozen out
                                !! of the ocean and the ice-ocean heat fluxes are
                                !! treated explicitly.
  logical :: use_limited_P_SSH  !< If true, return the sea surface height with
                                !! the correction for the atmospheric (and sea-ice)
                                !! pressure limited by max_p_surf instead of the
                                !! full atmospheric pressure.  The default is true.
  logical :: approx_net_mass_src !< If true, use the net mass sources from the ice-ocean boundary
                                !! type without any further adjustments to drive the ocean dynamics.
                                !! The actual net mass source may differ due to corrections.

  real :: gust_const            !< Constant unresolved background gustiness for ustar [R L Z T-1 ~> Pa]
  logical :: read_gust_2d       !< If true, use a 2-dimensional gustiness supplied from an input file.
  real, pointer, dimension(:,:) :: &
    TKE_tidal => NULL()         !< Turbulent kinetic energy introduced to the bottom boundary layer
                                !! by drag on the tidal flows [R Z3 T-3 ~> W m-2].
  real, pointer, dimension(:,:) :: &
    gust => NULL()              !< A spatially varying unresolved background gustiness that
                                !! contributes to ustar [R L Z T-1 ~> Pa].  gust is used when read_gust_2d is true.
  real, pointer, dimension(:,:) :: &
    ustar_tidal => NULL()       !< Tidal contribution to the bottom friction velocity [Z T-1 ~> m s-1]
  real :: cd_tides              !< Drag coefficient that applies to the tides (nondimensional)
  real :: utide                 !< Constant tidal velocity to use if read_tideamp is false [Z T-1 ~> m s-1].
  logical :: read_tideamp       !< If true, spatially varying tidal amplitude read from a file.

  logical :: rigid_sea_ice      !< If true, sea-ice exerts a rigidity that acts to damp surface
                                !! deflections (especially surface gravity waves).  The default is false.
  real    :: G_Earth            !< Gravitational acceleration [m s-2]
  real    :: Kv_sea_ice         !< Viscosity in sea-ice that resists sheared vertical motions [m2 s-1]
  real    :: density_sea_ice    !< Typical density of sea-ice [kg m-3]. The value is only used to convert
                                !! the ice pressure into appropriate units for use with Kv_sea_ice.
  real    :: rigid_sea_ice_mass !< A mass per unit area of sea-ice beyond which sea-ice viscosity
                                !! becomes effective [kg m-2], typically of order 1000 kg m-2.
  logical :: allow_flux_adjustments !< If true, use data_override to obtain flux adjustments

  logical :: restore_salt       !< If true, the coupled MOM driver adds a term to restore surface
                                !! salinity to a specified value.
  logical :: restore_temp       !< If true, the coupled MOM driver adds a term to restore sea
                                !! surface temperature to a specified value.
  real    :: Flux_const                     !< Piston velocity for surface restoring [Z T-1 ~> m s-1]
  logical :: salt_restore_as_sflux          !< If true, SSS restore as salt flux instead of water flux
  logical :: adjust_net_srestore_to_zero    !< Adjust srestore to zero (for both salt_flux or vprec)
  logical :: adjust_net_srestore_by_scaling !< Adjust srestore w/o moving zero contour
  logical :: adjust_net_fresh_water_to_zero !< Adjust net surface fresh-water (with restoring) to zero
  logical :: use_net_FW_adjustment_sign_bug !< Use the wrong sign when adjusting net FW
  logical :: adjust_net_fresh_water_by_scaling !< Adjust net surface fresh-water w/o moving zero contour
  logical :: mask_srestore_under_ice        !< If true, use an ice mask defined by frazil criteria
                                            !! for salinity restoring.
  real    :: ice_salt_concentration         !< Salt concentration for sea ice [kg/kg]
  logical :: mask_srestore_marginal_seas    !< If true, then mask SSS restoring in marginal seas
  real    :: max_delta_srestore             !< Maximum delta salinity used for restoring
  real    :: max_delta_trestore             !< Maximum delta sst used for restoring
  real, pointer, dimension(:,:) :: basin_mask => NULL() !< Mask for surface salinity restoring by basin
  logical :: answers_2018       !< If true, use the order of arithmetic and expressions that recover
                                !! the answers from the end of 2018.  Otherwise, use a simpler
                                !! expression to calculate gustiness.
  logical :: check_no_land_fluxes           !< Return warning if IOB flux over land is non-zero

  type(diag_ctrl), pointer :: diag => NULL()  !< Structure to regulate diagnostic output timing
  character(len=200) :: inputdir              !< Directory where NetCDF input files are
  character(len=200) :: salt_restore_file     !< Filename for salt restoring data
  character(len=30)  :: salt_restore_var_name !< Name of surface salinity in salt_restore_file
  logical            :: mask_srestore         !< If true, apply a 2-dimensional mask to the surface
                                              !! salinity restoring fluxes. The masking file should be
                                              !! in inputdir/salt_restore_mask.nc and the field should
                                              !! be named 'mask'
  real, pointer, dimension(:,:) :: srestore_mask => NULL() !< mask for SSS restoring
  character(len=200) :: temp_restore_file     !< Filename for sst restoring data
  character(len=30)  :: temp_restore_var_name !< Name of surface temperature in temp_restore_file
  logical            :: mask_trestore         !< If true, apply a 2-dimensional mask to the surface
                                              !! temperature restoring fluxes. The masking file should be
                                              !! in inputdir/temp_restore_mask.nc and the field should
                                              !! be named 'mask'
  real, pointer, dimension(:,:) :: trestore_mask => NULL() !< Mask for SST restoring
  integer :: id_srestore = -1  !< An id number for time_interp_external.
  integer :: id_trestore = -1  !< An id number for time_interp_external.

  type(forcing_diags), public :: handles !< Diagnostics handles

!#CTRL#  type(ctrl_forcing_CS), pointer :: ctrl_forcing_CSp => NULL()
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure
  type(user_revise_forcing_CS), pointer :: urf_CS => NULL() !< A control structure for user forcing revisions
end type surface_forcing_CS


!> ice_ocean_boundary_type is a structure corresponding to forcing, but with the elements, units,
!! and conventions that exactly conform to the use for MOM6-based coupled models.
type, public :: ice_ocean_boundary_type
  real, pointer, dimension(:,:) :: u_flux          =>NULL() !< i-direction wind stress [Pa]
  real, pointer, dimension(:,:) :: v_flux          =>NULL() !< j-direction wind stress [Pa]
  real, pointer, dimension(:,:) :: t_flux          =>NULL() !< sensible heat flux [W m-2]
  real, pointer, dimension(:,:) :: q_flux          =>NULL() !< specific humidity flux [kg m-2 s-1]
  real, pointer, dimension(:,:) :: salt_flux       =>NULL() !< salt flux [kg m-2 s-1]
  real, pointer, dimension(:,:) :: lw_flux         =>NULL() !< long wave radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_vis_dir =>NULL() !< direct visible sw radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_vis_dif =>NULL() !< diffuse visible sw radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_nir_dir =>NULL() !< direct Near InfraRed sw radiation [W m-2]
  real, pointer, dimension(:,:) :: sw_flux_nir_dif =>NULL() !< diffuse Near InfraRed sw radiation [W m-2]
  real, pointer, dimension(:,:) :: lprec           =>NULL() !< mass flux of liquid precip [kg m-2 s-1]
  real, pointer, dimension(:,:) :: fprec           =>NULL() !< mass flux of frozen precip [kg m-2 s-1]
  real, pointer, dimension(:,:) :: runoff          =>NULL() !< mass flux of liquid runoff [kg m-2 s-1]
  real, pointer, dimension(:,:) :: calving         =>NULL() !< mass flux of frozen runoff [kg m-2 s-1]
  real, pointer, dimension(:,:) :: stress_mag      =>NULL() !< The time-mean magnitude of the stress on the ocean [Pa]
  real, pointer, dimension(:,:) :: ustar_berg      =>NULL() !< frictional velocity beneath icebergs [m s-1]
  real, pointer, dimension(:,:) :: area_berg       =>NULL() !< fractional area covered by icebergs [m2 m-2]
  real, pointer, dimension(:,:) :: mass_berg       =>NULL() !< mass of icebergs per unit ocean area [kg m-2]
  real, pointer, dimension(:,:) :: runoff_hflx     =>NULL() !< heat content of liquid runoff [W m-2]
  real, pointer, dimension(:,:) :: calving_hflx    =>NULL() !< heat content of frozen runoff [W m-2]
  real, pointer, dimension(:,:) :: p               =>NULL() !< pressure of overlying ice and atmosphere
                                                            !< on ocean surface [Pa]
  real, pointer, dimension(:,:) :: mi              =>NULL() !< mass of ice per unit ocean area [kg m-2]
  real, pointer, dimension(:,:) :: ice_rigidity    =>NULL() !< rigidity of the sea ice, sea-ice and
                                                            !! ice-shelves, expressed as a coefficient
                                                            !! for divergence damping, as determined
                                                            !! outside of the ocean model [m3 s-1]
  integer :: xtype                    !< The type of the exchange - REGRID, REDIST or DIRECT
  type(coupler_2d_bc_type) :: fluxes  !< A structure that may contain an array of named fields
                                      !! used for passive tracer fluxes.
  integer :: wind_stagger = -999      !< A flag indicating the spatial discretization of wind stresses.
                                      !! This flag may be set by the flux-exchange code, based on what
                                      !! the sea-ice model is providing.  Otherwise, the value from
                                      !! the surface_forcing_CS is used.
end type ice_ocean_boundary_type

integer :: id_clock_forcing !< A CPU time clock

contains

!> This subroutine translates the Ice_ocean_boundary_type into a MOM
!! thermodynamic forcing type, including changes of units, sign conventions,
!! and putting the fields into arrays with MOM-standard halos.
subroutine convert_IOB_to_fluxes(IOB, fluxes, index_bounds, Time, valid_time, G, US, CS, sfc_state)
  type(ice_ocean_boundary_type), &
                   target, intent(in)    :: IOB    !< An ice-ocean boundary type with fluxes to drive
                                                   !! the ocean in a coupled model
  type(forcing),           intent(inout) :: fluxes !< A structure containing pointers to all
                                                   !! possible mass, heat or salt flux forcing fields.
                                                   !!  Unused fields have NULL ptrs.
  integer, dimension(4),   intent(in)    :: index_bounds !< The i- and j- size of the arrays in IOB.
  type(time_type),         intent(in)    :: Time   !< The time of the fluxes, used for interpolating the
                                                   !! salinity to the right time, when it is being restored.
  real,                    intent(in)    :: valid_time !< The amount of time over which these fluxes
                                                   !! should be applied [s].
  type(ocean_grid_type),   intent(inout) :: G      !< The ocean's grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS),pointer       :: CS     !< A pointer to the control structure returned by a
                                                   !! previous call to surface_forcing_init.
  type(surface),           intent(in)    :: sfc_state !< A structure containing fields that describe the
                                                   !! surface state of the ocean.

  real, dimension(SZI_(G),SZJ_(G)) :: &
    data_restore,  & ! The surface value toward which to restore [ppt] or [degC]
    SST_anom,      & ! Instantaneous sea surface temperature anomalies from a target value [degC]
    SSS_anom,      & ! Instantaneous sea surface salinity anomalies from a target value [ppt]
    SSS_mean,      & ! A (mean?) salinity about which to normalize local salinity
                     ! anomalies when calculating restorative precipitation anomalies [ppt]
    PmE_adj,       & ! The adjustment to PminusE that will cause the salinity
                     ! to be restored toward its target value [kg m-1 s-1]
    net_FW,        & ! The area integrated net freshwater flux into the ocean [kg s-1]
    net_FW2,       & ! The net freshwater flux into the ocean [kg m-2 s-1]
    work_sum,      & ! A 2-d array that is used as the work space for global sums [m2] or [kg s-1]
    open_ocn_mask    ! a binary field indicating where ice is present based on frazil criteria [nondim]

  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isr, ier, jsr, jer
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd

  real :: delta_sss           ! temporary storage for sss diff from restoring value [ppt]
  real :: delta_sst           ! temporary storage for sst diff from restoring value [degC]

  real :: kg_m2_s_conversion  ! A combination of unit conversion factors for rescaling
                              ! mass fluxes [R Z s m2 kg-1 T-1 ~> 1].
  real :: rhoXcp              ! Reference density times heat capacity times unit scaling
                              ! factors [J T s-1 Z-1 m-2 degC-1 ~> J m-3 degC-1]
  real :: sign_for_net_FW_bug ! Should be +1. but an old bug can be recovered by using -1.

  call cpu_clock_begin(id_clock_forcing)

  isc_bnd = index_bounds(1) ; iec_bnd = index_bounds(2)
  jsc_bnd = index_bounds(3) ; jec_bnd = index_bounds(4)
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  Isq  = G%IscB  ; Ieq  = G%IecB   ; Jsq  = G%JscB  ; Jeq  = G%JecB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB
  isr = is-isd+1 ; ier  = ie-isd+1 ; jsr = js-jsd+1 ; jer = je-jsd+1

  kg_m2_s_conversion = US%kg_m3_to_R*US%m_to_Z*US%T_to_s
  if (CS%restore_temp) rhoXcp = US%R_to_kg_m3*US%Z_to_m*US%s_to_T * CS%Rho0 * fluxes%C_p
  open_ocn_mask(:,:)     = 1.0
  pme_adj(:,:)           = 0.0
  fluxes%vPrecGlobalAdj  = 0.0
  fluxes%vPrecGlobalScl  = 0.0
  fluxes%saltFluxGlobalAdj = 0.0
  fluxes%saltFluxGlobalScl = 0.0
  fluxes%netFWGlobalAdj = 0.0
  fluxes%netFWGlobalScl = 0.0

  ! allocation and initialization if this is the first time that this
  ! flux type has been used.
  if (fluxes%dt_buoy_accum < 0) then
    call allocate_forcing_type(G, fluxes, water=.true., heat=.true., &
                               ustar=.true., press=.true.)

    call safe_alloc_ptr(fluxes%sw_vis_dir,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_vis_dif,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_nir_dir,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_nir_dif,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%p_surf,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%p_surf_full,isd,ied,jsd,jed)
    if (CS%use_limited_P_SSH) then
      fluxes%p_surf_SSH => fluxes%p_surf
    else
      fluxes%p_surf_SSH => fluxes%p_surf_full
    endif

    call safe_alloc_ptr(fluxes%salt_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%salt_flux_in,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%salt_flux_added,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%TKE_tidal,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%ustar_tidal,isd,ied,jsd,jed)

    if (CS%allow_flux_adjustments) then
      call safe_alloc_ptr(fluxes%heat_added,isd,ied,jsd,jed)
      call safe_alloc_ptr(fluxes%salt_flux_added,isd,ied,jsd,jed)
    endif

    do j=js-2,je+2 ; do i=is-2,ie+2
      fluxes%TKE_tidal(i,j)   = CS%TKE_tidal(i,j)
      fluxes%ustar_tidal(i,j) = CS%ustar_tidal(i,j)
    enddo ; enddo

    if (CS%restore_temp) call safe_alloc_ptr(fluxes%heat_added,isd,ied,jsd,jed)

  endif   ! endif for allocation and initialization


  if (((associated(IOB%ustar_berg) .and. (.not.associated(fluxes%ustar_berg))) &
    .or. (associated(IOB%area_berg) .and. (.not.associated(fluxes%area_berg)))) &
    .or. (associated(IOB%mass_berg) .and. (.not.associated(fluxes%mass_berg)))) &
    call allocate_forcing_type(G, fluxes, iceberg=.true.)

  if ((.not.coupler_type_initialized(fluxes%tr_fluxes)) .and. &
      coupler_type_initialized(IOB%fluxes)) &
    call coupler_type_spawn(IOB%fluxes, fluxes%tr_fluxes, &
                            (/is,is,ie,ie/), (/js,js,je,je/))
  !   It might prove valuable to use the same array extents as the rest of the
  ! ocean model, rather than using haloless arrays, in which case the last line
  ! would be: (             (/isd,is,ie,ied/), (/jsd,js,je,jed/))

  ! allocation and initialization on first call to this routine
  if (CS%area_surf < 0.0) then
    do j=js,je ; do i=is,ie
      work_sum(i,j) = US%L_to_m**2*G%areaT(i,j) * G%mask2dT(i,j)
    enddo ; enddo
    CS%area_surf = reproducing_sum(work_sum, isr, ier, jsr, jer)
  endif    ! endif for allocation and initialization


  ! Indicate that there are new unused fluxes.
  fluxes%fluxes_used = .false.
  fluxes%dt_buoy_accum = US%s_to_T*valid_time

  if (CS%allow_flux_adjustments) then
    fluxes%heat_added(:,:)=0.0
    fluxes%salt_flux_added(:,:)=0.0
  endif

  do j=js,je ; do i=is,ie
    fluxes%salt_flux(i,j) = 0.0
    fluxes%vprec(i,j) = 0.0
  enddo ; enddo

  ! Salinity restoring logic
  if (CS%restore_salt) then
    call time_interp_external(CS%id_srestore,Time,data_restore)
    ! open_ocn_mask indicates where to restore salinity (1 means restore, 0 does not)
    open_ocn_mask(:,:) = 1.0
    if (CS%mask_srestore_under_ice) then ! Do not restore under sea-ice
      do j=js,je ; do i=is,ie
        if (sfc_state%SST(i,j) <= -0.0539*sfc_state%SSS(i,j)) open_ocn_mask(i,j)=0.0
      enddo ; enddo
    endif
    if (CS%salt_restore_as_sflux) then
      do j=js,je ; do i=is,ie
        delta_sss = data_restore(i,j)- sfc_state%SSS(i,j)
        delta_sss = sign(1.0,delta_sss)*min(abs(delta_sss),CS%max_delta_srestore)
        fluxes%salt_flux(i,j) = 1.e-3*G%mask2dT(i,j) * (CS%Rho0*CS%Flux_const)* &
             (CS%basin_mask(i,j)*open_ocn_mask(i,j)*CS%srestore_mask(i,j)) *delta_sss  ! R Z T-1 ~> kg Salt m-2 s-1
      enddo ; enddo
      if (CS%adjust_net_srestore_to_zero) then
        if (CS%adjust_net_srestore_by_scaling) then
          call adjust_area_mean_to_zero(fluxes%salt_flux, G, fluxes%saltFluxGlobalScl, &
                          unit_scale=US%R_to_kg_m3*US%Z_to_m*US%s_to_T)
          fluxes%saltFluxGlobalAdj = 0.
        else
          work_sum(is:ie,js:je) = US%L_to_m**2*US%R_to_kg_m3*US%Z_to_m*US%s_to_T * &
                  G%areaT(is:ie,js:je)*fluxes%salt_flux(is:ie,js:je)
          fluxes%saltFluxGlobalAdj = reproducing_sum(work_sum(:,:), isr,ier, jsr,jer)/CS%area_surf
          fluxes%salt_flux(is:ie,js:je) = fluxes%salt_flux(is:ie,js:je) - kg_m2_s_conversion * fluxes%saltFluxGlobalAdj
        endif
      endif
      fluxes%salt_flux_added(is:ie,js:je) = fluxes%salt_flux(is:ie,js:je) ! Diagnostic
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0.5) then
          delta_sss = sfc_state%SSS(i,j) - data_restore(i,j)
          delta_sss = sign(1.0,delta_sss)*min(abs(delta_sss),CS%max_delta_srestore)
          fluxes%vprec(i,j) = (CS%basin_mask(i,j)*open_ocn_mask(i,j)*CS%srestore_mask(i,j))* &
                      (CS%Rho0*CS%Flux_const) * &
                      delta_sss / (0.5*(sfc_state%SSS(i,j) + data_restore(i,j)))
        endif
      enddo ; enddo
      if (CS%adjust_net_srestore_to_zero) then
        if (CS%adjust_net_srestore_by_scaling) then
          call adjust_area_mean_to_zero(fluxes%vprec, G, fluxes%vPrecGlobalScl, &
                       unit_scale=US%R_to_kg_m3*US%Z_to_m*US%s_to_T)
          fluxes%vPrecGlobalAdj = 0.
        else
          work_sum(is:ie,js:je) = US%L_to_m**2*G%areaT(is:ie,js:je) * &
                                  US%R_to_kg_m3*US%Z_to_m*US%s_to_T*fluxes%vprec(is:ie,js:je)
          fluxes%vPrecGlobalAdj = reproducing_sum(work_sum(:,:), isr, ier, jsr, jer) / CS%area_surf
          do j=js,je ; do i=is,ie
            fluxes%vprec(i,j) = ( fluxes%vprec(i,j) - kg_m2_s_conversion*fluxes%vPrecGlobalAdj ) * G%mask2dT(i,j)
          enddo ; enddo
        endif
      endif
    endif
  endif

  ! SST restoring logic
  if (CS%restore_temp) then
    call time_interp_external(CS%id_trestore,Time,data_restore)
    do j=js,je ; do i=is,ie
      delta_sst = data_restore(i,j)- sfc_state%SST(i,j)
      delta_sst = sign(1.0,delta_sst)*min(abs(delta_sst),CS%max_delta_trestore)
      fluxes%heat_added(i,j) = G%mask2dT(i,j) * CS%trestore_mask(i,j) * &
                               rhoXcp * delta_sst * CS%Flux_const  ! W m-2
    enddo ; enddo
  endif


  ! obtain fluxes from IOB; note the staggering of indices
  i0 = is - isc_bnd ; j0 = js - jsc_bnd
  do j=js,je ; do i=is,ie

    if (associated(IOB%lprec)) then
      fluxes%lprec(i,j) = kg_m2_s_conversion * IOB%lprec(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%lprec(i-i0,j-j0), G%mask2dT(i,j), i, j, 'lprec', G)
    endif

    if (associated(IOB%fprec)) then
      fluxes%fprec(i,j) = kg_m2_s_conversion * IOB%fprec(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%fprec(i-i0,j-j0), G%mask2dT(i,j), i, j, 'fprec', G)
    endif

    if (associated(IOB%q_flux)) then
      fluxes%evap(i,j) = - kg_m2_s_conversion * IOB%q_flux(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%q_flux(i-i0,j-j0), G%mask2dT(i,j), i, j, 'q_flux', G)
    endif

    if (associated(IOB%runoff)) then
      fluxes%lrunoff(i,j) = kg_m2_s_conversion * IOB%runoff(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%runoff(i-i0,j-j0), G%mask2dT(i,j), i, j, 'runoff', G)
    endif

    if (associated(IOB%calving)) then
      fluxes%frunoff(i,j) = kg_m2_s_conversion * IOB%calving(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%calving(i-i0,j-j0), G%mask2dT(i,j), i, j, 'calving', G)
    endif

    if (associated(IOB%ustar_berg)) then
      fluxes%ustar_berg(i,j) = US%m_to_Z*US%T_to_s * IOB%ustar_berg(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%ustar_berg(i-i0,j-j0), G%mask2dT(i,j), i, j, 'ustar_berg', G)
    endif

    if (associated(IOB%area_berg)) then
      fluxes%area_berg(i,j) = IOB%area_berg(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%area_berg(i-i0,j-j0), G%mask2dT(i,j), i, j, 'area_berg', G)
    endif

    if (associated(IOB%mass_berg)) then
      fluxes%mass_berg(i,j) = IOB%mass_berg(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%mass_berg(i-i0,j-j0), G%mask2dT(i,j), i, j, 'mass_berg', G)
    endif

    if (associated(IOB%runoff_hflx)) then
      fluxes%heat_content_lrunoff(i,j) = kg_m2_s_conversion * IOB%runoff_hflx(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%runoff_hflx(i-i0,j-j0), G%mask2dT(i,j), i, j, 'runoff_hflx', G)
    endif

    if (associated(IOB%calving_hflx)) then
      fluxes%heat_content_frunoff(i,j) = kg_m2_s_conversion * IOB%calving_hflx(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%calving_hflx(i-i0,j-j0), G%mask2dT(i,j), i, j, 'calving_hflx', G)
    endif

    if (associated(IOB%lw_flux)) then
      fluxes%LW(i,j) = IOB%lw_flux(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%lw_flux(i-i0,j-j0), G%mask2dT(i,j), i, j, 'lw_flux', G)
    endif

    if (associated(IOB%t_flux)) then
      fluxes%sens(i,j) = - IOB%t_flux(i-i0,j-j0) * G%mask2dT(i,j)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%t_flux(i-i0,j-j0), G%mask2dT(i,j), i, j, 't_flux', G)
    endif

    fluxes%latent(i,j) = 0.0
    if (associated(IOB%fprec)) then
      fluxes%latent(i,j)            = fluxes%latent(i,j) - IOB%fprec(i-i0,j-j0)*CS%latent_heat_fusion
      fluxes%latent_fprec_diag(i,j) = -G%mask2dT(i,j) * IOB%fprec(i-i0,j-j0)*CS%latent_heat_fusion
    endif
    if (associated(IOB%calving)) then
      fluxes%latent(i,j)              = fluxes%latent(i,j) - IOB%calving(i-i0,j-j0)*CS%latent_heat_fusion
      fluxes%latent_frunoff_diag(i,j) = -G%mask2dT(i,j) * IOB%calving(i-i0,j-j0)*CS%latent_heat_fusion
    endif
    if (associated(IOB%q_flux)) then
      fluxes%latent(i,j)           = fluxes%latent(i,j) - IOB%q_flux(i-i0,j-j0)*CS%latent_heat_vapor
      fluxes%latent_evap_diag(i,j) = -G%mask2dT(i,j) * IOB%q_flux(i-i0,j-j0)*CS%latent_heat_vapor
    endif

    fluxes%latent(i,j) = G%mask2dT(i,j) * fluxes%latent(i,j)

    if (associated(IOB%sw_flux_vis_dir)) then
      fluxes%sw_vis_dir(i,j) = G%mask2dT(i,j) * IOB%sw_flux_vis_dir(i-i0,j-j0)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%sw_flux_vis_dir(i-i0,j-j0), G%mask2dT(i,j), i, j, 'sw_flux_vis_dir', G)
    endif
    if (associated(IOB%sw_flux_vis_dif)) then
      fluxes%sw_vis_dif(i,j) = G%mask2dT(i,j) * IOB%sw_flux_vis_dif(i-i0,j-j0)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%sw_flux_vis_dif(i-i0,j-j0), G%mask2dT(i,j), i, j, 'sw_flux_vis_dif', G)
    endif
    if (associated(IOB%sw_flux_nir_dir)) then
      fluxes%sw_nir_dir(i,j) = G%mask2dT(i,j) * IOB%sw_flux_nir_dir(i-i0,j-j0)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%sw_flux_nir_dir(i-i0,j-j0), G%mask2dT(i,j), i, j, 'sw_flux_nir_dir', G)
    endif
    if (associated(IOB%sw_flux_nir_dif)) then
      fluxes%sw_nir_dif(i,j) = G%mask2dT(i,j) * IOB%sw_flux_nir_dif(i-i0,j-j0)
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%sw_flux_nir_dif(i-i0,j-j0), G%mask2dT(i,j), i, j, 'sw_flux_nir_dif', G)
    endif
    fluxes%sw(i,j) = fluxes%sw_vis_dir(i,j) + fluxes%sw_vis_dif(i,j) + &
                     fluxes%sw_nir_dir(i,j) + fluxes%sw_nir_dif(i,j)

  enddo ; enddo

  ! applied surface pressure from atmosphere and cryosphere
  if (associated(IOB%p)) then
    if (CS%max_p_surf >= 0.0) then
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * IOB%p(i-i0,j-j0)
        fluxes%p_surf(i,j) = MIN(fluxes%p_surf_full(i,j),CS%max_p_surf)
        if (CS%check_no_land_fluxes) &
          call check_mask_val_consistency(IOB%p(i-i0,j-j0), G%mask2dT(i,j), i, j, 'p', G)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * IOB%p(i-i0,j-j0)
        fluxes%p_surf(i,j) = fluxes%p_surf_full(i,j)
        if (CS%check_no_land_fluxes) &
          call check_mask_val_consistency(IOB%p(i-i0,j-j0), G%mask2dT(i,j), i, j, 'p', G)
      enddo ; enddo
    endif
    fluxes%accumulate_p_surf = .true. ! Multiple components may contribute to surface pressure.
  endif

  ! more salt restoring logic
  if (associated(IOB%salt_flux)) then
    do j=js,je ; do i=is,ie
      fluxes%salt_flux(i,j)    = G%mask2dT(i,j)*(fluxes%salt_flux(i,j) - kg_m2_s_conversion*IOB%salt_flux(i-i0,j-j0))
      fluxes%salt_flux_in(i,j) = G%mask2dT(i,j)*( -kg_m2_s_conversion*IOB%salt_flux(i-i0,j-j0) )
      if (CS%check_no_land_fluxes) &
        call check_mask_val_consistency(IOB%salt_flux(i-i0,j-j0), G%mask2dT(i,j), i, j, 'salt_flux', G)
    enddo ; enddo
  endif

!#CTRL# if (associated(CS%ctrl_forcing_CSp)) then
!#CTRL#   do j=js,je ; do i=is,ie
!#CTRL#     SST_anom(i,j) = sfc_state%SST(i,j) - CS%T_Restore(i,j)
!#CTRL#     SSS_anom(i,j) = sfc_state%SSS(i,j) - CS%S_Restore(i,j)
!#CTRL#     SSS_mean(i,j) = 0.5*(sfc_state%SSS(i,j) + CS%S_Restore(i,j))
!#CTRL#   enddo ; enddo
!#CTRL#   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_restore, &
!#CTRL#                           fluxes%vprec, day, dt, G, CS%ctrl_forcing_CSp)
!#CTRL# endif

  ! adjust the NET fresh-water flux to zero, if flagged
  if (CS%adjust_net_fresh_water_to_zero) then
    sign_for_net_FW_bug = 1.
    if (CS%use_net_FW_adjustment_sign_bug) sign_for_net_FW_bug = -1.
    do j=js,je ; do i=is,ie
      net_FW(i,j) = US%R_to_kg_m3*US%Z_to_m*US%s_to_T* &
                    (((fluxes%lprec(i,j)   + fluxes%fprec(i,j)) + &
                      (fluxes%lrunoff(i,j) + fluxes%frunoff(i,j))) + &
                      (fluxes%evap(i,j)    + fluxes%vprec(i,j)) ) * US%L_to_m**2*G%areaT(i,j)
      !   The following contribution appears to be calculating the volume flux of sea-ice
      ! melt. This calculation is clearly WRONG if either sea-ice has variable
      ! salinity or the sea-ice is completely fresh.
      !   Bob thinks this is trying ensure the net fresh-water of the ocean + sea-ice system
      ! is constant.
      !   To do this correctly we will need a sea-ice melt field added to IOB. -AJA
      if (associated(IOB%salt_flux) .and. (CS%ice_salt_concentration>0.0)) &
        net_FW(i,j) = net_FW(i,j) + sign_for_net_FW_bug * US%L_to_m**2*G%areaT(i,j) * &
                     (IOB%salt_flux(i-i0,j-j0) / CS%ice_salt_concentration)
      net_FW2(i,j) = net_FW(i,j) / (US%L_to_m**2*G%areaT(i,j))
    enddo ; enddo

    if (CS%adjust_net_fresh_water_by_scaling) then
      call adjust_area_mean_to_zero(net_FW2, G, fluxes%netFWGlobalScl)
      do j=js,je ; do i=is,ie
        fluxes%vprec(i,j) = fluxes%vprec(i,j) + US%kg_m3_to_R*US%m_to_Z*US%T_to_s * &
            (net_FW2(i,j) - net_FW(i,j)/(US%L_to_m**2*G%areaT(i,j))) * G%mask2dT(i,j)
      enddo ; enddo
    else
      fluxes%netFWGlobalAdj = reproducing_sum(net_FW(:,:), isr, ier, jsr, jer) / CS%area_surf
      do j=js,je ; do i=is,ie
        fluxes%vprec(i,j) = ( fluxes%vprec(i,j) - kg_m2_s_conversion * fluxes%netFWGlobalAdj ) * G%mask2dT(i,j)
      enddo ; enddo
    endif

  endif

  ! Set the wind stresses and ustar.
  if (associated(fluxes%ustar) .and. associated(fluxes%ustar_gustless)) then
    call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, ustar=fluxes%ustar, &
                              gustless_ustar=fluxes%ustar_gustless)
  elseif (associated(fluxes%ustar)) then
    call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, ustar=fluxes%ustar)
  elseif (associated(fluxes%ustar_gustless)) then
    call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, gustless_ustar=fluxes%ustar_gustless)
  endif

  if (coupler_type_initialized(fluxes%tr_fluxes) .and. &
      coupler_type_initialized(IOB%fluxes)) &
    call coupler_type_copy_data(IOB%fluxes, fluxes%tr_fluxes)

  if (CS%allow_flux_adjustments) then
    ! Apply adjustments to fluxes
    call apply_flux_adjustments(G, US, CS, Time, fluxes)
  endif

  ! Allow for user-written code to alter fluxes after all the above
  call user_alter_forcing(sfc_state, fluxes, Time, G, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)

end subroutine convert_IOB_to_fluxes

!> This subroutine translates the Ice_ocean_boundary_type into a MOM
!! mechanical forcing type, including changes of units, sign conventions,
!! and putting the fields into arrays with MOM-standard halos.
subroutine convert_IOB_to_forces(IOB, forces, index_bounds, Time, G, US, CS, dt_forcing, reset_avg)
  type(ice_ocean_boundary_type), &
                   target, intent(in)    :: IOB    !< An ice-ocean boundary type with fluxes to drive
                                                   !! the ocean in a coupled model
  type(mech_forcing),      intent(inout) :: forces !< A structure with the driving mechanical forces
  integer, dimension(4),   intent(in)    :: index_bounds !< The i- and j- size of the arrays in IOB.
  type(time_type),         intent(in)    :: Time   !< The time of the fluxes, used for interpolating the
                                                   !! salinity to the right time, when it is being restored.
  type(ocean_grid_type),   intent(inout) :: G      !< The ocean's grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(surface_forcing_CS),pointer       :: CS     !< A pointer to the control structure returned by a
                                                   !! previous call to surface_forcing_init.
  real,          optional, intent(in)    :: dt_forcing !< A time interval over which to apply the
                                                   !! current value of ustar as a weighted running
                                                   !! average [s], or if 0 do not average ustar.
                                                   !! Missing is equivalent to 0.
  logical,       optional, intent(in)    :: reset_avg !< If true, reset the time average.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    rigidity_at_h, &  ! Ice rigidity at tracer points [m3 s-1]
    net_mass_src, &   ! A temporary of net mass sources [kg m-2 s-1].
    ustar_tmp         ! A temporary array of ustar values [Z T-1 ~> m s-1].

  real :: I_GEarth      ! 1.0 / G_Earth [s2 m-1]
  real :: Kv_rho_ice    ! (CS%kv_sea_ice / CS%density_sea_ice) [m5 s-1 kg-1]
  real :: mass_ice      ! mass of sea ice at a face [kg m-2]
  real :: mass_eff      ! effective mass of sea ice for rigidity [kg m-2]
  real :: wt1, wt2      ! Relative weights of previous and current values of ustar, ND.

  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isr, ier, jsr, jer
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd

  call cpu_clock_begin(id_clock_forcing)

  isc_bnd = index_bounds(1) ; iec_bnd = index_bounds(2)
  jsc_bnd = index_bounds(3) ; jec_bnd = index_bounds(4)
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  Isq  = G%IscB  ; Ieq  = G%IecB   ; Jsq  = G%JscB  ; Jeq  = G%JecB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB
  isr = is-isd+1 ; ier  = ie-isd+1 ; jsr = js-jsd+1 ; jer = je-jsd+1
  i0 = is - isc_bnd ; j0 = js - jsc_bnd

  ! allocation and initialization if this is the first time that this
  ! mechanical forcing type has been used.
  if (.not.forces%initialized) then
    call allocate_mech_forcing(G, forces, stress=.true., ustar=.true., &
                               press=.true.)

    call safe_alloc_ptr(forces%p_surf,isd,ied,jsd,jed)
    call safe_alloc_ptr(forces%p_surf_full,isd,ied,jsd,jed)
    if (CS%use_limited_P_SSH) then
      forces%p_surf_SSH => forces%p_surf
    else
      forces%p_surf_SSH => forces%p_surf_full
    endif

    if (CS%rigid_sea_ice) then
      call safe_alloc_ptr(forces%rigidity_ice_u,IsdB,IedB,jsd,jed)
      call safe_alloc_ptr(forces%rigidity_ice_v,isd,ied,JsdB,JedB)
    endif

    forces%initialized = .true.
  endif

  if ( (associated(IOB%area_berg) .and. (.not. associated(forces%area_berg))) .or. &
       (associated(IOB%mass_berg) .and. (.not. associated(forces%mass_berg))) ) &
    call allocate_mech_forcing(G, forces, iceberg=.true.)

  if (associated(IOB%ice_rigidity)) then
    rigidity_at_h(:,:) = 0.0
    call safe_alloc_ptr(forces%rigidity_ice_u,IsdB,IedB,jsd,jed)
    call safe_alloc_ptr(forces%rigidity_ice_v,isd,ied,JsdB,JedB)
  endif

  forces%accumulate_rigidity = .true. ! Multiple components may contribute to rigidity.
  if (associated(forces%rigidity_ice_u)) forces%rigidity_ice_u(:,:) = 0.0
  if (associated(forces%rigidity_ice_v)) forces%rigidity_ice_v(:,:) = 0.0

  ! Set the weights for forcing fields that use running time averages.
  if (present(reset_avg)) then ; if (reset_avg) forces%dt_force_accum = 0.0 ; endif
  wt1 = 0.0 ; wt2 = 1.0
  if (present(dt_forcing)) then
    if ((forces%dt_force_accum > 0.0) .and. (dt_forcing > 0.0)) then
      wt1 = forces%dt_force_accum / (forces%dt_force_accum + dt_forcing)
      wt2 = 1.0 - wt1
    endif
    if (dt_forcing > 0.0) then
      forces%dt_force_accum = max(forces%dt_force_accum, 0.0) + dt_forcing
    else
      forces%dt_force_accum = 0.0 ! Reset the averaging time interval.
    endif
  else
    forces%dt_force_accum = 0.0 ! Reset the averaging time interval.
  endif

  ! applied surface pressure from atmosphere and cryosphere
  if (associated(IOB%p)) then
    if (CS%max_p_surf >= 0.0) then
      do j=js,je ; do i=is,ie
        forces%p_surf_full(i,j) = G%mask2dT(i,j) * IOB%p(i-i0,j-j0)
        forces%p_surf(i,j) = MIN(forces%p_surf_full(i,j),CS%max_p_surf)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        forces%p_surf_full(i,j) = G%mask2dT(i,j) * IOB%p(i-i0,j-j0)
        forces%p_surf(i,j) = forces%p_surf_full(i,j)
      enddo ; enddo
    endif
  else
    do j=js,je ; do i=is,ie
      forces%p_surf_full(i,j) = 0.0
      forces%p_surf(i,j) = 0.0
    enddo ; enddo
  endif
  forces%accumulate_p_surf = .true. ! Multiple components may contribute to surface pressure.

  ! Set the wind stresses and ustar.
  if (wt1 <= 0.0) then
    call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, taux=forces%taux, tauy=forces%tauy, &
                              ustar=forces%ustar, tau_halo=1)
  else
    call extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, taux=forces%taux, tauy=forces%tauy, &
                              ustar=ustar_tmp, tau_halo=1)
    do j=js,je ; do i=is,ie
      forces%ustar(i,j) = wt1*forces%ustar(i,j) + wt2*ustar_tmp(i,j)
    enddo ; enddo
  endif

  ! Find the net mass source in the input forcing without other adjustments.
  if (CS%approx_net_mass_src .and. associated(forces%net_mass_src)) then
    net_mass_src(:,:) = 0.0
    i0 = is - isc_bnd ; j0 = js - jsc_bnd
    do j=js,je ; do i=is,ie ; if (G%mask2dT(i,j) > 0.0) then
      if (associated(IOB%lprec)) &
        net_mass_src(i,j) = net_mass_src(i,j) + IOB%lprec(i-i0,j-j0)
      if (associated(IOB%fprec)) &
        net_mass_src(i,j) = net_mass_src(i,j) + IOB%fprec(i-i0,j-j0)
      if (associated(IOB%runoff)) &
        net_mass_src(i,j) = net_mass_src(i,j) + IOB%runoff(i-i0,j-j0)
      if (associated(IOB%calving)) &
        net_mass_src(i,j) = net_mass_src(i,j) + IOB%calving(i-i0,j-j0)
      if (associated(IOB%q_flux)) &
        net_mass_src(i,j) = net_mass_src(i,j) - IOB%q_flux(i-i0,j-j0)
    endif ; enddo ; enddo
    if (wt1 <= 0.0) then
      do j=js,je ; do i=is,ie
        forces%net_mass_src(i,j) = wt2*net_mass_src(i,j)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        forces%net_mass_src(i,j) = wt1*forces%net_mass_src(i,j) + wt2*net_mass_src(i,j)
      enddo ; enddo
    endif
    forces%net_mass_src_set = .true.
  else
    forces%net_mass_src_set = .false.
  endif

  ! Obtain optional ice-berg related fluxes from the IOB type:
  if (associated(IOB%area_berg)) then ; do j=js,je ; do i=is,ie
    forces%area_berg(i,j) = IOB%area_berg(i-i0,j-j0) * G%mask2dT(i,j)
  enddo ; enddo ; endif

  if (associated(IOB%mass_berg)) then ; do j=js,je ; do i=is,ie
    forces%mass_berg(i,j) = IOB%mass_berg(i-i0,j-j0) * G%mask2dT(i,j)
  enddo ; enddo ; endif

  ! Obtain sea ice related dynamic fields
  if (associated(IOB%ice_rigidity)) then
    do j=js,je ; do i=is,ie
      rigidity_at_h(i,j) = IOB%ice_rigidity(i-i0,j-j0) * G%mask2dT(i,j)
    enddo ; enddo
    call pass_var(rigidity_at_h, G%Domain, halo=1)
    do I=is-1,ie ; do j=js,je
      forces%rigidity_ice_u(I,j) = forces%rigidity_ice_u(I,j) + &
              min(rigidity_at_h(i,j), rigidity_at_h(i+1,j))
    enddo ; enddo
    do i=is,ie ; do J=js-1,je
      forces%rigidity_ice_v(i,J) = forces%rigidity_ice_v(i,J) + &
              min(rigidity_at_h(i,j), rigidity_at_h(i,j+1))
    enddo ; enddo
  endif

  if (CS%rigid_sea_ice) then
    call pass_var(forces%p_surf_full, G%Domain, halo=1)
    I_GEarth = 1.0 / CS%G_Earth
    Kv_rho_ice = (CS%kv_sea_ice / CS%density_sea_ice)
    do I=is-1,ie ; do j=js,je
      mass_ice = min(forces%p_surf_full(i,j), forces%p_surf_full(i+1,j)) * I_GEarth
      mass_eff = 0.0
      if (mass_ice > CS%rigid_sea_ice_mass) then
        mass_eff = (mass_ice - CS%rigid_sea_ice_mass) **2 / &
                   (mass_ice + CS%rigid_sea_ice_mass)
      endif
      forces%rigidity_ice_u(I,j) = forces%rigidity_ice_u(I,j) + Kv_rho_ice * mass_eff
    enddo ; enddo
    do i=is,ie ; do J=js-1,je
      mass_ice = min(forces%p_surf_full(i,j), forces%p_surf_full(i,j+1)) * I_GEarth
      mass_eff = 0.0
      if (mass_ice > CS%rigid_sea_ice_mass) then
        mass_eff = (mass_ice - CS%rigid_sea_ice_mass) **2 / &
                   (mass_ice + CS%rigid_sea_ice_mass)
      endif
      forces%rigidity_ice_v(i,J) = forces%rigidity_ice_v(i,J) + Kv_rho_ice * mass_eff
    enddo ; enddo
  endif

  if (CS%allow_flux_adjustments) then
    ! Apply adjustments to forces
    call apply_force_adjustments(G, US, CS, Time, forces)
  endif

!###  ! Allow for user-written code to alter fluxes after all the above
!###  call user_alter_mech_forcing(forces, Time, G, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end subroutine convert_IOB_to_forces


!> This subroutine extracts the wind stresses and related fields like ustar from an
!! Ice_ocean_boundary_type into optional argument arrays, including changes of units, sign
!! conventions, and putting the fields into arrays with MOM-standard sized halos.
subroutine extract_IOB_stresses(IOB, index_bounds, Time, G, US, CS, taux, tauy, ustar, &
                                gustless_ustar, tau_halo)
  type(ice_ocean_boundary_type), &
                   target, intent(in)    :: IOB  !< An ice-ocean boundary type with fluxes to drive
                                                 !! the ocean in a coupled model
  integer, dimension(4),   intent(in)    :: index_bounds !< The i- and j- size of the arrays in IOB.
  type(time_type),         intent(in)    :: Time !< The time of the fluxes, used for interpolating the
                                                 !! salinity to the right time, when it is being restored.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(surface_forcing_CS),pointer       :: CS   !< A pointer to the control structure returned by a
                                                 !! previous call to surface_forcing_init.
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(inout) :: taux !< The zonal wind stresses on a C-grid [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(inout) :: tauy !< The meridional wind stresses on a C-grid [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(inout) :: ustar !< The surface friction velocity [Z T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(out)   :: gustless_ustar !< The surface friction velocity without
                                                 !! any contributions from gustiness [Z T-1 ~> m s-1].
  integer,       optional, intent(in)    :: tau_halo !< The halo size of wind stresses to set, 0 by default.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: taux_in_A   ! Zonal wind stresses [R Z L T-2 ~> Pa] at h points
  real, dimension(SZI_(G),SZJ_(G)) :: tauy_in_A   ! Meridional wind stresses [R Z L T-2 ~> Pa] at h points
  real, dimension(SZIB_(G),SZJ_(G)) :: taux_in_C  ! Zonal wind stresses [Pa] at u points
  real, dimension(SZI_(G),SZJB_(G)) :: tauy_in_C  ! Meridional wind stresses [R Z L T-2 ~> Pa] at v points
  real, dimension(SZIB_(G),SZJB_(G)) :: taux_in_B ! Zonal wind stresses [Pa] at q points
  real, dimension(SZIB_(G),SZJB_(G)) :: tauy_in_B ! Meridional wind stresses [R Z L T-2 ~> Pa] at q points

  real :: gustiness     ! unresolved gustiness that contributes to ustar [R Z L T-2 ~> Pa]
  real :: Irho0         ! Inverse of the mean density rescaled to [Z L-1 R-1 ~> m3 kg-1]
  real :: taux2, tauy2  ! squared wind stresses [R2 Z2 L2 T-4 ~> Pa2]
  real :: tau_mag       ! magnitude of the wind stress [R Z L T-2 ~> Pa]
  real :: Pa_conversion ! A unit conversion factor from Pa to the internal wind stress units [R Z L T-2 Pa-1 ~> 1]
  real :: stress_conversion ! A unit conversion factor from Pa times any stress multiplier [R Z L T-2 Pa-1 ~> 1]

  logical :: do_ustar, do_gustless
  integer :: wind_stagger  ! AGRID, BGRID_NE, or CGRID_NE (integers from MOM_domains)
  integer :: i, j, is, ie, js, je, ish, ieh, jsh, jeh, Isqh, Ieqh, Jsqh, Jeqh, i0, j0, halo

  halo = 0 ; if (present(tau_halo)) halo = tau_halo
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  ish  = G%isc-halo  ; ieh   = G%iec+halo  ; jsh  = G%jsc-halo  ; jeh  = G%jec+halo
  Isqh = G%IscB-halo ; Ieqh  = G%IecB+halo ; Jsqh = G%JscB-halo ; Jeqh = G%JecB+halo
  i0 = is - index_bounds(1) ; j0 = js - index_bounds(3)

  IRho0 = US%L_to_Z / CS%Rho0
  Pa_conversion = US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z
  stress_conversion = Pa_conversion * CS%wind_stress_multiplier

  do_ustar = present(ustar) ; do_gustless = present(gustless_ustar)

  wind_stagger = CS%wind_stagger
  if ((IOB%wind_stagger == AGRID) .or. (IOB%wind_stagger == BGRID_NE) .or. &
      (IOB%wind_stagger == CGRID_NE)) wind_stagger = IOB%wind_stagger

  if (associated(IOB%u_flux).neqv.associated(IOB%v_flux)) call MOM_error(FATAL,"extract_IOB_stresses: "//&
            "associated(IOB%u_flux) /= associated(IOB%v_flux !!!")
  if (present(taux).neqv.present(tauy)) call MOM_error(FATAL,"extract_IOB_stresses: "//&
            "present(taux) /= present(tauy) !!!")

  ! Set surface momentum stress related fields as a function of staggering.
  if (present(taux) .or. present(tauy) .or. &
      ((do_ustar.or.do_gustless) .and. .not.associated(IOB%stress_mag)) ) then

    if (wind_stagger == BGRID_NE) then
      taux_in_B(:,:) = 0.0 ; tauy_in_B(:,:) = 0.0
      if (associated(IOB%u_flux).and.associated(IOB%v_flux)) then
        do J=js,je ; do I=is,ie
          taux_in_B(I,J) = IOB%u_flux(i-i0,j-j0) * stress_conversion
          tauy_in_B(I,J) = IOB%v_flux(i-i0,j-j0) * stress_conversion
        enddo ; enddo
      endif

      if (G%symmetric) call fill_symmetric_edges(taux_in_B, tauy_in_B, G%Domain, stagger=BGRID_NE)
      call pass_vector(taux_in_B, tauy_in_B, G%Domain, stagger=BGRID_NE, halo=max(1,halo))

      if (present(taux).and.present(tauy)) then
        do j=jsh,jeh ; do I=Isqh,Ieqh
          taux(I,j) = 0.0
          if ((G%mask2dBu(I,J) + G%mask2dBu(I,J-1)) > 0) &
            taux(I,j) = (G%mask2dBu(I,J)*taux_in_B(I,J) + G%mask2dBu(I,J-1)*taux_in_B(I,J-1)) / &
                        (G%mask2dBu(I,J) + G%mask2dBu(I,J-1))
        enddo ; enddo
        do J=Jsqh,Jeqh ; do i=ish,ieh
          tauy(i,J) = 0.0
          if ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J)) > 0) &
            tauy(i,J) = (G%mask2dBu(I,J)*tauy_in_B(I,J) + G%mask2dBu(I-1,J)*tauy_in_B(I-1,J)) / &
                        (G%mask2dBu(I,J) + G%mask2dBu(I-1,J))
        enddo ; enddo
      endif
    elseif (wind_stagger == AGRID) then
      taux_in_A(:,:) = 0.0 ; tauy_in_A(:,:) = 0.0
      if (associated(IOB%u_flux).and.associated(IOB%v_flux)) then
        do j=js,je ; do i=is,ie
          taux_in_A(i,j) = IOB%u_flux(i-i0,j-j0) * stress_conversion
          tauy_in_A(i,j) = IOB%v_flux(i-i0,j-j0) * stress_conversion
        enddo ; enddo
      endif

      if (halo == 0) then
        call pass_vector(taux_in_A, tauy_in_A, G%Domain, To_All+Omit_Corners, stagger=AGRID, halo=1)
      else
        call pass_vector(taux_in_A, tauy_in_A, G%Domain, stagger=AGRID, halo=max(1,halo))
      endif

      if (present(taux)) then ; do j=jsh,jeh ; do I=Isqh,Ieqh
        taux(I,j) = 0.0
        if ((G%mask2dT(i,j) + G%mask2dT(i+1,j)) > 0) &
          taux(I,j) = (G%mask2dT(i,j)*taux_in_A(i,j) + G%mask2dT(i+1,j)*taux_in_A(i+1,j)) / &
                      (G%mask2dT(i,j) + G%mask2dT(i+1,j))
      enddo ; enddo ; endif

      if (present(tauy)) then ; do J=Jsqh,Jeqh ; do i=ish,ieh
        tauy(i,J) = 0.0
        if ((G%mask2dT(i,j) + G%mask2dT(i,j+1)) > 0) &
          tauy(i,J) = (G%mask2dT(i,j)*tauy_in_A(i,j) + G%mask2dT(i,J+1)*tauy_in_A(i,j+1)) / &
                      (G%mask2dT(i,j) + G%mask2dT(i,j+1))
      enddo ; enddo ; endif

    else ! C-grid wind stresses.
      taux_in_C(:,:) = 0.0 ; tauy_in_C(:,:) = 0.0
      if (associated(IOB%u_flux).and.associated(IOB%v_flux)) then
        do j=js,je ; do i=is,ie
          taux_in_C(I,j) = IOB%u_flux(i-i0,j-j0) * stress_conversion
          tauy_in_C(i,J) = IOB%v_flux(i-i0,j-j0) * stress_conversion
        enddo ; enddo
      endif

      if (G%symmetric) call fill_symmetric_edges(taux_in_C, tauy_in_C, G%Domain)
      call pass_vector(taux_in_C, tauy_in_C, G%Domain, halo=max(1,halo))

      if (present(taux).and.present(tauy)) then
        do j=jsh,jeh ; do I=Isqh,Ieqh
          taux(I,j) = G%mask2dCu(I,j)*taux_in_C(I,j)
        enddo ; enddo
        do J=Jsqh,Jeqh ; do i=ish,ieh
          tauy(i,J) = G%mask2dCv(i,J)*tauy_in_C(i,J)
        enddo ; enddo
      endif
    endif   ! endif for extracting wind stress fields with various staggerings
  endif

  if (do_ustar .or. do_gustless) then
    ! Set surface friction velocity directly or as a function of staggering.
    ! ustar is required for the bulk mixed layer formulation and other turbulent mixing
    ! parametizations. The background gustiness (for example with a relatively small value
    ! of 0.02 Pa) is intended to give reasonable behavior in regions of very weak winds.
    if (associated(IOB%stress_mag)) then
      if (do_ustar) then ; do j=js,je ; do i=is,ie
        gustiness = CS%gust_const
        if (CS%read_gust_2d) then
          if ((wind_stagger == CGRID_NE) .or. &
              ((wind_stagger == AGRID) .and. (G%mask2dT(i,j) > 0)) .or. &
              ((wind_stagger == BGRID_NE) .and. &
               (((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + &
                (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) > 0)) ) &
            gustiness = CS%gust(i,j)
        endif
        ustar(i,j) = sqrt(gustiness*IRho0 + IRho0*Pa_conversion*IOB%stress_mag(i-i0,j-j0))
      enddo ; enddo ; endif
      if (CS%answers_2018) then
        if (do_gustless) then ; do j=js,je ; do i=is,ie
          gustless_ustar(i,j) = sqrt(Pa_conversion*US%L_to_Z*IOB%stress_mag(i-i0,j-j0) / CS%Rho0)
        enddo ; enddo ; endif
      else
        if (do_gustless) then ; do j=js,je ; do i=is,ie
          gustless_ustar(i,j) = sqrt(IRho0 * Pa_conversion*IOB%stress_mag(i-i0,j-j0))
        enddo ; enddo ; endif
      endif
    elseif (wind_stagger == BGRID_NE) then
      do j=js,je ; do i=is,ie
        tau_mag = 0.0 ; gustiness = CS%gust_const
        if (((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + &
             (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) > 0) then
          tau_mag = sqrt(((G%mask2dBu(I,J)*(taux_in_B(I,J)**2 + tauy_in_B(I,J)**2) + &
              G%mask2dBu(I-1,J-1)*(taux_in_B(I-1,J-1)**2 + tauy_in_B(I-1,J-1)**2)) + &
             (G%mask2dBu(I,J-1)*(taux_in_B(I,J-1)**2 + tauy_in_B(I,J-1)**2) + &
              G%mask2dBu(I-1,J)*(taux_in_B(I-1,J)**2 + tauy_in_B(I-1,J)**2)) ) / &
            ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) )
          if (CS%read_gust_2d) gustiness = CS%gust(i,j)
        endif
        if (do_ustar) ustar(i,j) = sqrt(gustiness*IRho0 + IRho0 * tau_mag)
        if (CS%answers_2018) then
          if (do_gustless) gustless_ustar(i,j) = sqrt(US%L_to_Z*tau_mag / CS%Rho0)
        else
          if (do_gustless) gustless_ustar(i,j) = sqrt(IRho0 * tau_mag)
        endif
      enddo ; enddo
    elseif (wind_stagger == AGRID) then
      do j=js,je ; do i=is,ie
        tau_mag = G%mask2dT(i,j) * sqrt(taux_in_A(i,j)**2 + tauy_in_A(i,j)**2)
        gustiness = CS%gust_const
        if (CS%read_gust_2d .and. (G%mask2dT(i,j) > 0)) gustiness = CS%gust(i,j)
        if (do_ustar) ustar(i,j) = sqrt(gustiness*IRho0 + IRho0 * tau_mag)
        if (CS%answers_2018) then
          if (do_gustless) gustless_ustar(i,j) = sqrt(US%L_to_Z*tau_mag / CS%Rho0)
        else
          if (do_gustless) gustless_ustar(i,j) = sqrt(IRho0 * tau_mag)
        endif
      enddo ; enddo
    else  ! C-grid wind stresses.
      do j=js,je ; do i=is,ie
        taux2 = 0.0 ; tauy2 = 0.0
        if ((G%mask2dCu(I-1,j) + G%mask2dCu(I,j)) > 0) &
          taux2 = (G%mask2dCu(I-1,j)*taux_in_C(I-1,j)**2 + G%mask2dCu(I,j)*taux_in_C(I,j)**2) / &
                  (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))
        if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0) &
          tauy2 = (G%mask2dCv(i,J-1)*tauy_in_C(i,J-1)**2 + G%mask2dCv(i,J)*tauy_in_C(i,J)**2) / &
                  (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))
        tau_mag = sqrt(taux2 + tauy2)

        gustiness = CS%gust_const
        if (CS%read_gust_2d) gustiness = CS%gust(i,j)

        if (do_ustar) ustar(i,j) = sqrt(gustiness*IRho0 + IRho0 * tau_mag)
        if (CS%answers_2018) then
          if (do_gustless) gustless_ustar(i,j) = sqrt(US%L_to_Z*tau_mag / CS%Rho0)
        else
          if (do_gustless) gustless_ustar(i,j) = sqrt(IRho0 * tau_mag)
        endif
      enddo ; enddo
    endif ! endif for wind friction velocity fields
  endif

end subroutine extract_IOB_stresses


!> Adds thermodynamic flux adjustments obtained via data_override
!! Component name is 'OCN'
!! Available adjustments are:
!! - hflx_adj (Heat flux into the ocean [W m-2])
!! - sflx_adj (Salt flux into the ocean [kg salt m-2 s-1])
!! - prcme_adj (Fresh water flux into the ocean [kg m-2 s-1])
subroutine apply_flux_adjustments(G, US, CS, Time, fluxes)
  type(ocean_grid_type),    intent(inout) :: G  !< Ocean grid structure
  type(unit_scale_type),    intent(in)    :: US !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS !< Surface forcing control structure
  type(time_type),          intent(in)    :: Time !< Model time structure
  type(forcing),            intent(inout) :: fluxes !< Surface fluxes structure

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: temp_at_h ! Various fluxes at h points [W m-2] or [kg m-2 s-1]

  integer :: isc, iec, jsc, jec, i, j
  logical :: overrode_h

  isc = G%isc; iec = G%iec ; jsc = G%jsc; jec = G%jec

  overrode_h = .false.
  call data_override('OCN', 'hflx_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then ; do j=jsc,jec ; do i=isc,iec
    fluxes%heat_added(i,j) = fluxes%heat_added(i,j) + temp_at_h(i,j)* G%mask2dT(i,j)
  enddo ; enddo ; endif
  ! Not needed? ! if (overrode_h) call pass_var(fluxes%heat_added, G%Domain)

  overrode_h = .false.
  call data_override('OCN', 'sflx_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then ; do j=jsc,jec ; do i=isc,iec
    fluxes%salt_flux_added(i,j) = fluxes%salt_flux_added(i,j) + &
        US%kg_m3_to_R*US%m_to_Z*US%T_to_s * temp_at_h(i,j)* G%mask2dT(i,j)
  enddo ; enddo ; endif
  ! Not needed? ! if (overrode_h) call pass_var(fluxes%salt_flux_added, G%Domain)

  overrode_h = .false.
  call data_override('OCN', 'prcme_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then ; do j=jsc,jec ; do i=isc,iec
    fluxes%vprec(i,j) = fluxes%vprec(i,j) + US%kg_m3_to_R*US%m_to_Z*US%T_to_s * temp_at_h(i,j)* G%mask2dT(i,j)
  enddo ; enddo ; endif
  ! Not needed? ! if (overrode_h) call pass_var(fluxes%vprec, G%Domain)
end subroutine apply_flux_adjustments

!> Adds mechanical forcing adjustments obtained via data_override
!! Component name is 'OCN'
!! Available adjustments are:
!! - taux_adj (Zonal wind stress delta, positive to the east [Pa])
!! - tauy_adj (Meridional wind stress delta, positive to the north [Pa])
subroutine apply_force_adjustments(G, US, CS, Time, forces)
  type(ocean_grid_type),    intent(inout) :: G  !< Ocean grid structure
  type(unit_scale_type),    intent(in)    :: US !< A dimensional unit scaling type
  type(surface_forcing_CS), pointer       :: CS !< Surface forcing control structure
  type(time_type),          intent(in)    :: Time !< Model time structure
  type(mech_forcing),       intent(inout) :: forces !< A structure with the driving mechanical forces

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: tempx_at_h ! Delta to zonal wind stress at h points [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)) :: tempy_at_h ! Delta to meridional wind stress at h points [R Z L T-2 ~> Pa]

  integer :: isc, iec, jsc, jec, i, j
  real :: dLonDx, dLonDy, rDlon, cosA, sinA, zonal_tau, merid_tau
  real :: Pa_conversion ! A unit conversion factor from Pa to the internal units [R Z L T-2 Pa-1 ~> 1]
  logical :: overrode_x, overrode_y

  isc = G%isc; iec = G%iec ; jsc = G%jsc; jec = G%jec
  Pa_conversion = US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z

  tempx_at_h(:,:) = 0.0 ; tempy_at_h(:,:) = 0.0
  ! Either reads data or leaves contents unchanged
  overrode_x = .false. ; overrode_y = .false.
  call data_override('OCN', 'taux_adj', tempx_at_h(isc:iec,jsc:jec), Time, override=overrode_x)
  call data_override('OCN', 'tauy_adj', tempy_at_h(isc:iec,jsc:jec), Time, override=overrode_y)

  if (overrode_x .or. overrode_y) then
    if (.not. (overrode_x .and. overrode_y)) call MOM_error(FATAL,"apply_flux_adjustments: "//&
            "Both taux_adj and tauy_adj must be specified, or neither, in data_table")

    ! Rotate winds
    call pass_vector(tempx_at_h, tempy_at_h, G%Domain, To_All, AGRID, halo=1)
    do j=jsc-1,jec+1 ; do i=isc-1,iec+1
      dLonDx = G%geoLonCu(I,j) - G%geoLonCu(I-1,j)
      dLonDy = G%geoLonCv(i,J) - G%geoLonCv(i,J-1)
      rDlon = sqrt( dLonDx * dLonDx + dLonDy * dLonDy )
      if (rDlon > 0.) rDlon = 1. / rDlon
      cosA = dLonDx * rDlon
      sinA = dLonDy * rDlon
      zonal_tau = Pa_conversion * tempx_at_h(i,j)
      merid_tau = Pa_conversion * tempy_at_h(i,j)
      tempx_at_h(i,j) = cosA * zonal_tau - sinA * merid_tau
      tempy_at_h(i,j) = sinA * zonal_tau + cosA * merid_tau
    enddo ; enddo

    ! Average to C-grid locations
    do j=jsc,jec ; do I=isc-1,iec
      forces%taux(I,j) = forces%taux(I,j) + 0.5 * ( tempx_at_h(i,j) + tempx_at_h(i+1,j) )
    enddo ; enddo

    do J=jsc-1,jec ; do i=isc,iec
      forces%tauy(i,J) = forces%tauy(i,J) + 0.5 * ( tempy_at_h(i,j) + tempy_at_h(i,j+1) )
    enddo ; enddo
  endif ! overrode_x .or. overrode_y

end subroutine apply_force_adjustments

!> Save any restart files associated with the surface forcing.
subroutine forcing_save_restart(CS, G, Time, directory, time_stamped, &
                                filename_suffix)
  type(surface_forcing_CS),   pointer       :: CS   !< A pointer to the control structure returned
                                                    !! by a previous call to surface_forcing_init
  type(ocean_grid_type),      intent(inout) :: G    !< The ocean's grid structure
  type(time_type),            intent(in)    :: Time !< The current model time
  character(len=*),           intent(in)    :: directory !< The directory into which to write the
                                                    !! restart files
  logical,          optional, intent(in)    :: time_stamped !< If true, the restart file names include
                                                    !! a unique time stamp.  The default is false.
  character(len=*), optional, intent(in)    :: filename_suffix !< An optional suffix (e.g., a time-
                                                    !! stamp) to append to the restart file names.

  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return
  call save_restart(directory, Time, G, CS%restart_CSp, time_stamped)

end subroutine forcing_save_restart

!> Initialize the surface forcing, including setting parameters and allocating permanent memory.
subroutine surface_forcing_init(Time, G, US, param_file, diag, CS)
  type(time_type),          intent(in)    :: Time !< The current model time
  type(ocean_grid_type),    intent(in)    :: G    !< The ocean's grid structure
  type(unit_scale_type),    intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),    intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target,  intent(inout) :: diag !< A structure that is used to regulate
                                                  !! diagnostic output
  type(surface_forcing_CS), pointer       :: CS   !< A pointer that is set to point to the control
                                                  !! structure for this module

  ! Local variables
  real :: utide  ! The RMS tidal velocity [Z T-1 ~> m s-1].
  type(directories)  :: dirs
  logical            :: new_sim, iceberg_flux_diags
  logical            :: default_2018_answers
  type(time_type)    :: Time_frc
  character(len=200) :: TideAmp_file, gust_file, salt_file, temp_file ! Input file names.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_surface_forcing"  ! This module's name.
  character(len=48)  :: stagger
  character(len=48)  :: flnam
  character(len=240) :: basin_file
  integer :: i, j, isd, ied, jsd, jed

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
    call MOM_error(WARNING, "surface_forcing_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  id_clock_forcing=cpu_clock_id('Ocean surface forcing', grain=CLOCK_SUBCOMPONENT)
  call cpu_clock_begin(id_clock_forcing)

  CS%diag => diag

  call write_version_number(version)
  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, &
                 "The directory in which all input files are found.", &
                 default=".")
  CS%inputdir = slasher(CS%inputdir)
  call get_param(param_file, mdl, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state "//&
                 "variables.", default=.true.)
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", CS%latent_heat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf)
  call get_param(param_file, mdl, "LATENT_HEAT_VAPORIZATION", CS%latent_heat_vapor, &
                 "The latent heat of fusion.", units="J/kg", default=hlv)
  call get_param(param_file, mdl, "MAX_P_SURF", CS%max_p_surf, &
                 "The maximum surface pressure that can be exerted by the "//&
                 "atmosphere and floating sea-ice or ice shelves. This is "//&
                 "needed because the FMS coupling structure does not "//&
                 "limit the water that can be frozen out of the ocean and "//&
                 "the ice-ocean heat fluxes are treated explicitly.  No "//&
                 "limit is applied if a negative value is used.", units="Pa", &
                 default=-1.0)
  call get_param(param_file, mdl, "RESTORE_SALINITY", CS%restore_salt, &
                 "If true, the coupled driver will add a globally-balanced "//&
                 "fresh-water flux that drives sea-surface salinity "//&
                 "toward specified values.", default=.false.)
  call get_param(param_file, mdl, "RESTORE_TEMPERATURE", CS%restore_temp, &
                 "If true, the coupled driver will add a  "//&
                 "heat flux that drives sea-surface temperature "//&
                 "toward specified values.", default=.false.)
  call get_param(param_file, mdl, "ADJUST_NET_SRESTORE_TO_ZERO", &
                 CS%adjust_net_srestore_to_zero, &
                 "If true, adjusts the salinity restoring seen to zero "//&
                 "whether restoring is via a salt flux or virtual precip.",&
                 default=CS%restore_salt)
  call get_param(param_file, mdl, "ADJUST_NET_SRESTORE_BY_SCALING", &
                 CS%adjust_net_srestore_by_scaling, &
                 "If true, adjustments to salt restoring to achieve zero net are "//&
                 "made by scaling values without moving the zero contour.",&
                 default=.false.)
  call get_param(param_file, mdl, "ADJUST_NET_FRESH_WATER_TO_ZERO", &
                 CS%adjust_net_fresh_water_to_zero, &
                 "If true, adjusts the net fresh-water forcing seen "//&
                 "by the ocean (including restoring) to zero.", default=.false.)
  if (CS%adjust_net_fresh_water_to_zero) &
    call get_param(param_file, mdl, "USE_NET_FW_ADJUSTMENT_SIGN_BUG", &
                 CS%use_net_FW_adjustment_sign_bug, &
                   "If true, use the wrong sign for the adjustment to "//&
                   "the net fresh-water.", default=.true.)
  call get_param(param_file, mdl, "ADJUST_NET_FRESH_WATER_BY_SCALING", &
                 CS%adjust_net_fresh_water_by_scaling, &
                 "If true, adjustments to net fresh water to achieve zero net are "//&
                 "made by scaling values without moving the zero contour.",&
                 default=.false.)
  call get_param(param_file, mdl, "ICE_SALT_CONCENTRATION", &
                 CS%ice_salt_concentration, &
                 "The assumed sea-ice salinity needed to reverse engineer the "//&
                 "melt flux (or ice-ocean fresh-water flux).", &
                 units="kg/kg", default=0.005)
  call get_param(param_file, mdl, "USE_LIMITED_PATM_SSH", CS%use_limited_P_SSH, &
                 "If true, return the sea surface height with the "//&
                 "correction for the atmospheric (and sea-ice) pressure "//&
                 "limited by max_p_surf instead of the full atmospheric "//&
                 "pressure.", default=.true.)
  call get_param(param_file, mdl, "APPROX_NET_MASS_SRC", CS%approx_net_mass_src, &
                 "If true, use the net mass sources from the ice-ocean "//&
                 "boundary type without any further adjustments to drive "//&
                 "the ocean dynamics.  The actual net mass source may differ "//&
                 "due to internal corrections.", default=.false.)

  call get_param(param_file, mdl, "WIND_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the "//&
                 "staggering of the input wind stress field.  Valid "//&
                 "values are 'A', 'B', or 'C'.", default="C")
  if (uppercase(stagger(1:1)) == 'A') then ; CS%wind_stagger = AGRID
  elseif (uppercase(stagger(1:1)) == 'B') then ; CS%wind_stagger = BGRID_NE
  elseif (uppercase(stagger(1:1)) == 'C') then ; CS%wind_stagger = CGRID_NE
  else ; call MOM_error(FATAL,"surface_forcing_init: WIND_STAGGER = "// &
                        trim(stagger)//" is invalid.") ; endif
  call get_param(param_file, mdl, "WIND_STRESS_MULTIPLIER", CS%wind_stress_multiplier, &
                 "A factor multiplying the wind-stress given to the ocean by the "//&
                 "coupler. This is used for testing and should be =1.0 for any "//&
                 "production runs.", default=1.0)

  if (CS%restore_salt) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes "//&
                 "to the relative surface anomalies (akin to a piston "//&
                 "velocity).  Note the non-MKS units.", &
                 units="m day-1", scale=US%m_to_Z*US%T_to_s, fail_if_missing=.true.)
    call get_param(param_file, mdl, "SALT_RESTORE_FILE", CS%salt_restore_file, &
                 "A file in which to find the surface salinity to use for restoring.", &
                 default="salt_restore.nc")
    call get_param(param_file, mdl, "SALT_RESTORE_VARIABLE", CS%salt_restore_var_name, &
                 "The name of the surface salinity variable to read from "//&
                 "SALT_RESTORE_FILE for restoring salinity.", &
                 default="salt")
! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0

    call get_param(param_file, mdl, "SRESTORE_AS_SFLUX", CS%salt_restore_as_sflux, &
                 "If true, the restoring of salinity is applied as a salt "//&
                 "flux instead of as a freshwater flux.", default=.false.)
    call get_param(param_file, mdl, "MAX_DELTA_SRESTORE", CS%max_delta_srestore, &
                 "The maximum salinity difference used in restoring terms.", &
                 units="PSU or g kg-1", default=999.0)
    call get_param(param_file, mdl, "MASK_SRESTORE_UNDER_ICE", &
                 CS%mask_srestore_under_ice, &
                 "If true, disables SSS restoring under sea-ice based on a frazil "//&
                 "criteria (SST<=Tf). Only used when RESTORE_SALINITY is True.",      &
                 default=.false.)
    call get_param(param_file, mdl, "MASK_SRESTORE_MARGINAL_SEAS", &
                 CS%mask_srestore_marginal_seas, &
                 "If true, disable SSS restoring in marginal seas. Only used when "//&
                 "RESTORE_SALINITY is True.", default=.false.)
    call get_param(param_file, mdl, "BASIN_FILE", basin_file, &
                 "A file in which to find the basin masks, in variable 'basin'.", &
                 default="basin.nc")
    basin_file = trim(CS%inputdir) // trim(basin_file)
    call safe_alloc_ptr(CS%basin_mask,isd,ied,jsd,jed) ; CS%basin_mask(:,:) = 1.0
    if (CS%mask_srestore_marginal_seas) then
      call MOM_read_data(basin_file,'basin',CS%basin_mask,G%domain, timelevel=1)
      do j=jsd,jed ; do i=isd,ied
        if (CS%basin_mask(i,j) >= 6.0) then ; CS%basin_mask(i,j) = 0.0
        else ; CS%basin_mask(i,j) = 1.0 ; endif
      enddo ; enddo
    endif
    call get_param(param_file, mdl, "MASK_SRESTORE", CS%mask_srestore, &
                 "If true, read a file (salt_restore_mask) containing "//&
                 "a mask for SSS restoring.", default=.false.)
  endif

  if (CS%restore_temp) then
    call get_param(param_file, mdl, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes "//&
                 "to the relative surface anomalies (akin to a piston "//&
                 "velocity).  Note the non-MKS units.", &
                 units="m day-1", scale=US%m_to_Z*US%T_to_s, fail_if_missing=.true.)
    call get_param(param_file, mdl, "SST_RESTORE_FILE", CS%temp_restore_file, &
                 "A file in which to find the surface temperature to use for restoring.", &
                 default="temp_restore.nc")
    call get_param(param_file, mdl, "SST_RESTORE_VARIABLE", CS%temp_restore_var_name, &
                 "The name of the surface temperature variable to read from "//&
                 "SST_RESTORE_FILE for restoring sst.", &
                 default="temp")
  ! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0

    call get_param(param_file, mdl, "MAX_DELTA_TRESTORE", CS%max_delta_trestore, &
                 "The maximum sst difference used in restoring terms.", &
                 units="degC ", default=999.0)
    call get_param(param_file, mdl, "MASK_TRESTORE", CS%mask_trestore, &
                 "If true, read a file (temp_restore_mask) containing "//&
                 "a mask for SST restoring.", default=.false.)

  endif

! Optionally read tidal amplitude from input file [m s-1] on model grid.
! Otherwise use default tidal amplitude for bottom frictionally-generated
! dissipation. Default cd_tides is chosen to yield approx 1 TWatt of
! work done against tides globally using OSU tidal amplitude.
  call get_param(param_file, mdl, "CD_TIDES", CS%cd_tides, &
                 "The drag coefficient that applies to the tides.", &
                 units="nondim", default=1.0e-4)
  call get_param(param_file, mdl, "READ_TIDEAMP", CS%read_TIDEAMP, &
                 "If true, read a file (given by TIDEAMP_FILE) containing "//&
                 "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
  if (CS%read_TIDEAMP) then
    call get_param(param_file, mdl, "TIDEAMP_FILE", TideAmp_file, &
                 "The path to the file containing the spatially varying "//&
                 "tidal amplitudes with INT_TIDE_DISSIPATION.", &
                 default="tideamp.nc")
    CS%utide=0.0
  else
    call get_param(param_file, mdl, "UTIDE", CS%utide, &
                 "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
                 units="m s-1", default=0.0, scale=US%m_to_Z*US%T_to_s)
  endif

  call safe_alloc_ptr(CS%TKE_tidal,isd,ied,jsd,jed)
  call safe_alloc_ptr(CS%ustar_tidal,isd,ied,jsd,jed)

  if (CS%read_TIDEAMP) then
    TideAmp_file = trim(CS%inputdir) // trim(TideAmp_file)
    call MOM_read_data(TideAmp_file,'tideamp',CS%TKE_tidal,G%domain,timelevel=1, scale=US%m_to_Z*US%T_to_s)
    do j=jsd, jed; do i=isd, ied
      utide = CS%TKE_tidal(i,j)
      CS%TKE_tidal(i,j) = G%mask2dT(i,j)*CS%Rho0*CS%cd_tides*(utide*utide*utide)
      CS%ustar_tidal(i,j) = sqrt(CS%cd_tides)*utide
    enddo ; enddo
  else
    do j=jsd,jed; do i=isd,ied
      utide = CS%utide
      CS%TKE_tidal(i,j) = CS%Rho0*CS%cd_tides*(utide*utide*utide)
      CS%ustar_tidal(i,j) = sqrt(CS%cd_tides)*utide
    enddo ; enddo
  endif

  call time_interp_external_init

  ! Optionally read a x-y gustiness field in place of a global constant.
  call get_param(param_file, mdl, "READ_GUST_2D", CS%read_gust_2d, &
                 "If true, use a 2-dimensional gustiness supplied from "//&
                 "an input file", default=.false.)
  call get_param(param_file, mdl, "GUST_CONST", CS%gust_const, &
               "The background gustiness in the winds.", &
               units="Pa", default=0.02, scale=US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z)
  if (CS%read_gust_2d) then
    call get_param(param_file, mdl, "GUST_2D_FILE", gust_file, &
                 "The file in which the wind gustiness is found in "//&
                 "variable gustiness.")

    call safe_alloc_ptr(CS%gust,isd,ied,jsd,jed)
    gust_file = trim(CS%inputdir) // trim(gust_file)
    call MOM_read_data(gust_file, 'gustiness', CS%gust, G%domain, timelevel=1, &
               scale=US%kg_m3_to_R*US%m_s_to_L_T**2*US%L_to_Z) ! units in file should be Pa
  endif
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.true.)
  call get_param(param_file, mdl, "SURFACE_FORCING_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the answers "//&
                 "from the end of 2018.  Otherwise, use a simpler expression to calculate gustiness.", &
                 default=default_2018_answers)

! See whether sufficiently thick sea ice should be treated as rigid.
  call get_param(param_file, mdl, "USE_RIGID_SEA_ICE", CS%rigid_sea_ice, &
                 "If true, sea-ice is rigid enough to exert a "//&
                 "nonhydrostatic pressure that resist vertical motion.", &
                 default=.false.)
  if (CS%rigid_sea_ice) then
    call get_param(param_file, mdl, "G_EARTH", CS%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)
    call get_param(param_file, mdl, "SEA_ICE_MEAN_DENSITY", CS%density_sea_ice, &
                 "A typical density of sea ice, used with the kinematic "//&
                 "viscosity, when USE_RIGID_SEA_ICE is true.", units="kg m-3", &
                 default=900.0)
    call get_param(param_file, mdl, "SEA_ICE_VISCOSITY", CS%Kv_sea_ice, &
                 "The kinematic viscosity of sufficiently thick sea ice "//&
                 "for use in calculating the rigidity of sea ice.", &
                 units="m2 s-1", default=1.0e9)
    call get_param(param_file, mdl, "SEA_ICE_RIGID_MASS", CS%rigid_sea_ice_mass, &
                 "The mass of sea-ice per unit area at which the sea-ice "//&
                 "starts to exhibit rigidity", units="kg m-2", default=1000.0)
  endif

  call get_param(param_file, mdl, "ALLOW_ICEBERG_FLUX_DIAGNOSTICS", iceberg_flux_diags, &
                 "If true, makes available diagnostics of fluxes from icebergs "//&
                 "as seen by MOM6.", default=.false.)
  call register_forcing_type_diags(Time, diag, US, CS%use_temperature, CS%handles, &
                                   use_berg_fluxes=iceberg_flux_diags)

  call get_param(param_file, mdl, "ALLOW_FLUX_ADJUSTMENTS", CS%allow_flux_adjustments, &
                 "If true, allows flux adjustments to specified via the "//&
                 "data_table using the component name 'OCN'.", default=.false.)

  call get_param(param_file, mdl, "CHECK_NO_LAND_FLUXES", CS%check_no_land_fluxes, &
                 "If true, checks that values from IOB fluxes are zero "//&
                 "above land points (i.e. G%mask2dT = 0).", default=.false., &
                 debuggingParam=.true.)

  call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)

  if (CS%restore_salt) then
    salt_file = trim(CS%inputdir) // trim(CS%salt_restore_file)
    CS%id_srestore = init_external_field(salt_file, CS%salt_restore_var_name, domain=G%Domain%mpp_domain)
    call safe_alloc_ptr(CS%srestore_mask,isd,ied,jsd,jed); CS%srestore_mask(:,:) = 1.0
    if (CS%mask_srestore) then ! read a 2-d file containing a mask for restoring fluxes
      flnam = trim(CS%inputdir) // 'salt_restore_mask.nc'
      call MOM_read_data(flnam,'mask', CS%srestore_mask, G%domain, timelevel=1)
    endif
  endif

  if (CS%restore_temp) then
    temp_file = trim(CS%inputdir) // trim(CS%temp_restore_file)
    CS%id_trestore = init_external_field(temp_file, CS%temp_restore_var_name, domain=G%Domain%mpp_domain)
    call safe_alloc_ptr(CS%trestore_mask,isd,ied,jsd,jed); CS%trestore_mask(:,:) = 1.0
    if (CS%mask_trestore) then  ! read a 2-d file containing a mask for restoring fluxes
      flnam = trim(CS%inputdir) // 'temp_restore_mask.nc'
      call MOM_read_data(flnam, 'mask', CS%trestore_mask, G%domain, timelevel=1)
    endif
  endif

  ! Set up any restart fields associated with the forcing.
  call restart_init(param_file, CS%restart_CSp, "MOM_forcing.res")
!#CTRL#  call register_ctrl_forcing_restarts(G, param_file, CS%ctrl_forcing_CSp, &
!#CTRL#                                      CS%restart_CSp)
  call restart_init_end(CS%restart_CSp)

  if (associated(CS%restart_CSp)) then
    call Get_MOM_Input(dirs=dirs)

    new_sim = .false.
    if ((dirs%input_filename(1:1) == 'n') .and. &
        (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.
    if (.not.new_sim) then
      call restore_state(dirs%input_filename, dirs%restart_input_dir, Time_frc, &
                         G, CS%restart_CSp)
    endif
  endif

!#CTRL#  call controlled_forcing_init(Time, G, param_file, diag, CS%ctrl_forcing_CSp)

  call user_revise_forcing_init(param_file, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end subroutine surface_forcing_init

!> Clean up and deallocate any memory associated with this module and its children.
subroutine surface_forcing_end(CS, fluxes)
  type(surface_forcing_CS), pointer       :: CS !< A pointer to the control structure returned by
                                                !! a previous call to surface_forcing_init, it will
                                                !! be deallocated here.
  type(forcing), optional,  intent(inout) :: fluxes !< A structure containing pointers to all
                                                !! possible mass, heat or salt flux forcing fields.
                                                !! If present, it will be deallocated here.

  if (present(fluxes)) call deallocate_forcing_type(fluxes)

!#CTRL#  call controlled_forcing_end(CS%ctrl_forcing_CSp)

  if (associated(CS)) deallocate(CS)
  CS => NULL()

end subroutine surface_forcing_end

!> Write out a set of messages with checksums of the fields in an ice_ocen_boundary type
subroutine ice_ocn_bnd_type_chksum(id, timestep, iobt)

  character(len=*), intent(in) :: id     !< An identifying string for this call
  integer,          intent(in) :: timestep !< The number of elapsed timesteps
  type(ice_ocean_boundary_type), &
                    intent(in) :: iobt   !< An ice-ocean boundary type with fluxes to drive the
                                         !! ocean in a coupled model whose checksums are reported
  integer ::   n,m, outunit

  outunit = stdout()

  write(outunit,*) "BEGIN CHECKSUM(ice_ocean_boundary_type):: ", id, timestep
  write(outunit,100) 'iobt%u_flux         ', mpp_chksum( iobt%u_flux         )
  write(outunit,100) 'iobt%v_flux         ', mpp_chksum( iobt%v_flux         )
  write(outunit,100) 'iobt%t_flux         ', mpp_chksum( iobt%t_flux         )
  write(outunit,100) 'iobt%q_flux         ', mpp_chksum( iobt%q_flux         )
  write(outunit,100) 'iobt%salt_flux      ', mpp_chksum( iobt%salt_flux      )
  write(outunit,100) 'iobt%lw_flux        ', mpp_chksum( iobt%lw_flux        )
  write(outunit,100) 'iobt%sw_flux_vis_dir', mpp_chksum( iobt%sw_flux_vis_dir)
  write(outunit,100) 'iobt%sw_flux_vis_dif', mpp_chksum( iobt%sw_flux_vis_dif)
  write(outunit,100) 'iobt%sw_flux_nir_dir', mpp_chksum( iobt%sw_flux_nir_dir)
  write(outunit,100) 'iobt%sw_flux_nir_dif', mpp_chksum( iobt%sw_flux_nir_dif)
  write(outunit,100) 'iobt%lprec          ', mpp_chksum( iobt%lprec          )
  write(outunit,100) 'iobt%fprec          ', mpp_chksum( iobt%fprec          )
  write(outunit,100) 'iobt%runoff         ', mpp_chksum( iobt%runoff         )
  write(outunit,100) 'iobt%calving        ', mpp_chksum( iobt%calving        )
  write(outunit,100) 'iobt%p              ', mpp_chksum( iobt%p              )
  if (associated(iobt%ustar_berg)) &
    write(outunit,100) 'iobt%ustar_berg     ', mpp_chksum( iobt%ustar_berg )
  if (associated(iobt%area_berg)) &
    write(outunit,100) 'iobt%area_berg      ', mpp_chksum( iobt%area_berg  )
  if (associated(iobt%mass_berg)) &
    write(outunit,100) 'iobt%mass_berg      ', mpp_chksum( iobt%mass_berg  )
100 FORMAT("   CHECKSUM::",A20," = ",Z20)

  call coupler_type_write_chksums(iobt%fluxes, outunit, 'iobt%')

end subroutine ice_ocn_bnd_type_chksum

!> Check the values passed by IOB over land are zero
subroutine check_mask_val_consistency(val, mask, i, j, varname, G)

  real, intent(in) :: val  !< value of flux/variable passed by IOB
  real, intent(in) :: mask !< value of ocean mask
  integer, intent(in) :: i, j !< model grid cell indices
  character(len=*), intent(in) :: varname !< variable name
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure
  ! Local variables
  character(len=48) :: ci, cj !< model local grid cell indices as strings
  character(len=48) :: ciglo, cjglo !< model global grid cell indices as strings
  character(len=48) :: cval !< value to be displayed
  character(len=256) :: error_message !< error message to be displayed

  if ((mask == 0.) .and. (val /= 0.)) then
    write(ci, '(I8)') i
    write(cj, '(I8)') j
    write(ciglo, '(I8)') i + G%HI%idg_offset
    write(cjglo, '(I8)') j + G%HI%jdg_offset
    write(cval, '(E22.16)') val
    error_message = "MOM_surface_forcing: found non-zero value (="//trim(cval)//") over land "//&
                    "for variable "//trim(varname)//" at local point (i, j) = ("//trim(ci)//", "//trim(cj)//&
                    ", global point (iglo, jglo) = ("//trim(ciglo)//", "//trim(cjglo)//")"
    call MOM_error(WARNING, error_message)
  endif

end subroutine

end module MOM_surface_forcing_gfdl
