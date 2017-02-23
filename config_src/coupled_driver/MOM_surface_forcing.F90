module MOM_surface_forcing
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, May 2004                                       *
!*  Edited by Stephen Griffies June 2014                               *
!*                                                                     *
!*    This program contains the subroutines that transform the surface *
!*  forcing fields for the coupled model and manage the output of      *
!*  these fields.                                                      *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, tauy                                  *
!*    j    x ^ x ^ x   At >:  u, taux                                  *
!*    j    > o > o >   At o:  h, fluxes.                               *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
!### use MOM_controlled_forcing, only : apply_ctrl_forcing, register_ctrl_forcing_restarts
!### use MOM_controlled_forcing, only : controlled_forcing_init, controlled_forcing_end
!### use MOM_controlled_forcing, only : ctrl_forcing_CS
use MOM_coms,             only : reproducing_sum
use MOM_constants,        only : hlv, hlf
use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_SUBCOMPONENT
use MOM_diag_mediator,    only : diag_ctrl
use MOM_diag_mediator,    only : safe_alloc_ptr, time_type
use MOM_domains,          only : pass_vector, pass_var, fill_symmetric_edges
use MOM_domains,          only : global_field_sum, BITWISE_EXACT_SUM
use MOM_domains,          only : AGRID, BGRID_NE, CGRID_NE, To_All
use MOM_error_handler,    only : MOM_error, WARNING, FATAL, is_root_pe, MOM_mesg
use MOM_file_parser,      only : get_param, log_version, param_file_type
use MOM_forcing_type,     only : forcing, forcing_diags, register_forcing_type_diags
use MOM_forcing_type,     only : allocate_forcing_type, deallocate_forcing_type
use MOM_get_input,        only : Get_MOM_Input, directories
use MOM_grid,             only : ocean_grid_type
use MOM_io,               only : slasher, write_version_number
use MOM_restart,          only : register_restart_field, restart_init, MOM_restart_CS
use MOM_restart,          only : restart_init_end, save_restart, restore_state
use MOM_string_functions, only : uppercase
use MOM_spatial_means,    only : adjust_area_mean_to_zero
use MOM_variables,        only : surface
use user_revise_forcing,  only : user_alter_forcing, user_revise_forcing_init, user_revise_forcing_CS

use coupler_types_mod,        only : coupler_2d_bc_type
use data_override_mod,        only : data_override_init, data_override
use fms_mod,                  only : read_data, stdout
use mpp_mod,                  only : mpp_chksum
use time_interp_external_mod, only : init_external_field, time_interp_external
use time_interp_external_mod, only : time_interp_external_init

implicit none ; private

#include <MOM_memory.h>

public convert_IOB_to_fluxes
public surface_forcing_init
public ice_ocn_bnd_type_chksum
public forcing_save_restart


! surface_forcing_CS is a structure containing pointers to the forcing fields
! which may be used to drive MOM.  All fluxes are positive downward.
type, public :: surface_forcing_CS ; private
  integer :: wind_stagger       ! AGRID, BGRID_NE, or CGRID_NE (integer values
                                ! from MOM_domains) to indicate the staggering of
                                ! the winds that are being provided in calls to
                                ! update_ocean_model.
  logical :: use_temperature    ! If true, temp and saln used as state variables
  real :: wind_stress_multiplier !< A multiplier applied to incoming wind stress (nondim).

  ! smg: remove when have A=B code reconciled
  logical :: bulkmixedlayer     ! If true, model based on bulk mixed layer code

  real :: Rho0                  ! Boussinesq reference density (kg/m^3)
  real :: area_surf = -1.0      ! total ocean surface area (m^2)
  real :: latent_heat_fusion    ! latent heat of fusion (J/kg)
  real :: latent_heat_vapor     ! latent heat of vaporization (J/kg)

  real :: max_p_surf            ! maximum surface pressure that can be
                                ! exerted by the atmosphere and floating sea-ice,
                                ! in Pa.  This is needed because the FMS coupling
                                ! structure does not limit the water that can be
                                ! frozen out of the ocean and the ice-ocean heat
                                ! fluxes are treated explicitly.
  logical :: use_limited_P_SSH  ! If true, return the sea surface height with
                                ! the correction for the atmospheric (and sea-ice)
                                ! pressure limited by max_p_surf instead of the
                                ! full atmospheric pressure.  The default is true.

  real :: gust_const            ! constant unresolved background gustiness for ustar (Pa)
  logical :: read_gust_2d       ! If true, use a 2-dimensional gustiness supplied
                                ! from an input file.
  real, pointer, dimension(:,:) :: &
    TKE_tidal => NULL(), &      ! turbulent kinetic energy introduced to the
                                ! bottom boundary layer by drag on the tidal flows,
                                ! in W m-2.
    gust => NULL(), &           ! spatially varying unresolved background
                                ! gustiness that contributes to ustar (Pa).
                                ! gust is used when read_gust_2d is true.
    ustar_tidal => NULL()       ! tidal contribution to the bottom friction velocity (m/s)
  real :: cd_tides              ! drag coefficient that applies to the tides (nondimensional)
  real :: utide                 ! constant tidal velocity to use if read_tideamp
                                ! is false, in m s-1.
  logical :: read_tideamp       ! If true, spatially varying tidal amplitude read from a file.

  logical :: rigid_sea_ice      ! If true, sea-ice exerts a rigidity that acts
                                ! to damp surface deflections (especially surface
                                ! gravity waves).  The default is false.
  real    :: Kv_sea_ice         ! viscosity in sea-ice that resists sheared vertical motions (m^2/s)
  real    :: density_sea_ice    ! typical density of sea-ice (kg/m^3). The value is
                                ! only used to convert the ice pressure into
                                ! appropriate units for use with Kv_sea_ice.
  real    :: rigid_sea_ice_mass ! A mass per unit area of sea-ice beyond which
                                ! sea-ice viscosity becomes effective, in kg m-2,
                                ! typically of order 1000 kg m-2.
  logical :: allow_flux_adjustments ! If true, use data_override to obtain flux adjustments

  real    :: Flux_const                     ! piston velocity for surface restoring (m/s)
  logical :: salt_restore_as_sflux          ! If true, SSS restore as salt flux instead of water flux
  logical :: adjust_net_srestore_to_zero    ! adjust srestore to zero (for both salt_flux or vprec)
  logical :: adjust_net_srestore_by_scaling ! adjust srestore w/o moving zero contour
  logical :: adjust_net_fresh_water_to_zero ! adjust net surface fresh-water (w/ restoring) to zero
  logical :: adjust_net_fresh_water_by_scaling ! adjust net surface fresh-water  w/o moving zero contour
  logical :: mask_srestore_under_ice        ! If true, use an ice mask defined by frazil
                                            ! criteria for salinity restoring.
  real    :: ice_salt_concentration         ! salt concentration for sea ice (kg/kg)
  logical :: mask_srestore_marginal_seas    ! if true, then mask SSS restoring in marginal seas
  real    :: max_delta_srestore             ! maximum delta salinity used for restoring
  real    :: max_delta_trestore             ! maximum delta sst used for restoring
  real, pointer, dimension(:,:) :: basin_mask => NULL() ! mask for SSS restoring

  type(diag_ctrl), pointer :: diag                  ! structure to regulate diagnostic output timing
  character(len=200)       :: inputdir              ! directory where NetCDF input files are
  character(len=200)       :: salt_restore_file     ! filename for salt restoring data
  character(len=30)        :: salt_restore_var_name ! name of surface salinity in salt_restore_file
  character(len=200)       :: temp_restore_file     ! filename for sst restoring data
  character(len=30)        :: temp_restore_var_name ! name of surface temperature in temp_restore_file

  integer :: id_srestore = -1     ! id number for time_interp_external.
  integer :: id_trestore = -1     ! id number for time_interp_external.

  ! Diagnostics handles
  type(forcing_diags), public :: handles

!###  type(ctrl_forcing_CS), pointer :: ctrl_forcing_CSp => NULL()
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()
  type(user_revise_forcing_CS), pointer :: urf_CS => NULL()
end type surface_forcing_CS


! ice_ocean_boundary_type is a structure corresponding to forcing, but with
! the elements, units, and conventions that exactly conform to the use for
! MOM-based coupled models.
type, public :: ice_ocean_boundary_type
  real, pointer, dimension(:,:) :: u_flux               =>NULL() ! i-direction wind stress (Pa)
  real, pointer, dimension(:,:) :: v_flux               =>NULL() ! j-direction wind stress (Pa)
  real, pointer, dimension(:,:) :: t_flux               =>NULL() ! sensible heat flux (W/m2)
  real, pointer, dimension(:,:) :: q_flux               =>NULL() ! specific humidity flux (kg/m2/s)
  real, pointer, dimension(:,:) :: salt_flux            =>NULL() ! salt flux (kg/m2/s)
  real, pointer, dimension(:,:) :: lw_flux              =>NULL() ! long wave radiation (W/m2)
  real, pointer, dimension(:,:) :: sw_flux_vis_dir      =>NULL() ! direct visible sw radiation (W/m2)
  real, pointer, dimension(:,:) :: sw_flux_vis_dif      =>NULL() ! diffuse visible sw radiation (W/m2)
  real, pointer, dimension(:,:) :: sw_flux_nir_dir      =>NULL() ! direct Near InfraRed sw radiation (W/m2)
  real, pointer, dimension(:,:) :: sw_flux_nir_dif      =>NULL() ! diffuse Near InfraRed sw radiation (W/m2)
  real, pointer, dimension(:,:) :: lprec                =>NULL() ! mass flux of liquid precip (kg/m2/s)
  real, pointer, dimension(:,:) :: fprec                =>NULL() ! mass flux of frozen precip (kg/m2/s)
  real, pointer, dimension(:,:) :: runoff               =>NULL() ! mass flux of liquid runoff (kg/m2/s)
  real, pointer, dimension(:,:) :: calving              =>NULL() ! mass flux of frozen runoff (kg/m2/s)
  real, pointer, dimension(:,:) :: ustar_berg           =>NULL() ! frictional velocity beneath icebergs (m/s)
  real, pointer, dimension(:,:) :: area_berg            =>NULL() ! area covered by icebergs(m2/m2)
  real, pointer, dimension(:,:) :: mass_berg            =>NULL() ! mass of icebergs(kg/m2)
  real, pointer, dimension(:,:) :: runoff_hflx          =>NULL() ! heat content of liquid runoff (W/m2)
  real, pointer, dimension(:,:) :: calving_hflx         =>NULL() ! heat content of frozen runoff (W/m2)
  real, pointer, dimension(:,:) :: p                    =>NULL() ! pressure of overlying ice and atmosphere
                                                                 ! on ocean surface (Pa)
  real, pointer, dimension(:,:) :: mi                   =>NULL() ! mass of ice (kg/m2)
  integer :: xtype                                   ! REGRID, REDIST or DIRECT
  type(coupler_2d_bc_type)      :: fluxes            ! A structure that may contain an
                                                     ! array of named fields used for
                                                     ! passive tracer fluxes.
  integer :: wind_stagger = -999                     ! A flag indicating the spatial discretization of
                                                     ! wind stresses.  This flag may be set by the
                                                     ! flux-exchange code, based on what the sea-ice
                                                     ! model is providing.  Otherwise, the value from
                                                     ! the surface_forcing_CS is used.
end type ice_ocean_boundary_type

integer :: id_clock_forcing

contains

subroutine convert_IOB_to_fluxes(IOB, fluxes, index_bounds, Time, G, CS, state, restore_salt, restore_temp)
  type(ice_ocean_boundary_type), intent(in), target :: IOB
  type(forcing),              intent(inout)         :: fluxes
  integer, dimension(4),      intent(in)            :: index_bounds
  type(time_type),            intent(in)            :: Time
  type(ocean_grid_type),      intent(inout)         :: G
  type(surface_forcing_CS),   pointer               :: CS
  type(surface),              intent(in)            :: state
  logical, optional,          intent(in)            :: restore_salt, restore_temp

! This subroutine translates the Ice_ocean_boundary_type into a
! MOM forcing type, including changes of units, sign conventions,
! and puting the fields into arrays with MOM-standard halos.

! Arguments:
!   IOB  ice-ocean boundary type w/ fluxes to drive ocean in a coupled model
!  (out) fluxes - A structure containing pointers to any possible
!                 forcing fields.  Unused fields have NULL ptrs.
!  (in)  index_bounds - the i- and j- size of the arrays in IOB.
!  (in)  Time - The time of the fluxes, used for interpolating the salinity
!               to the right time, when it is being restored.
!  (in)  G - The ocean's grid structure.
!  (in)  CS - A pointer to the control structure returned by a previous
!             call to surface_forcing_init.
!  (in)  state - A structure containing fields that describe the
!                surface state of the ocean.
!  (in)  restore_salt - if true, salinity is restored to a target value.
!  (in)  restore_temp - if true, temperature is restored to a target value.

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    taux_at_q, &     ! Zonal wind stresses at q points (Pa)
    tauy_at_q        ! Meridional wind stresses at q points (Pa)

  real, dimension(SZI_(G),SZJ_(G)) :: &
    taux_at_h,     & ! Zonal wind stresses at h points (Pa)
    tauy_at_h,     & ! Meridional wind stresses at h points (Pa)
    data_restore,  & ! The surface value toward which to restore (g/kg or degC)
    SST_anom,      & ! Instantaneous sea surface temperature anomalies from a target value (deg C)
    SSS_anom,      & ! Instantaneous sea surface salinity anomalies from a target value (g/kg)
    SSS_mean,      & ! A (mean?) salinity about which to normalize local salinity
                     ! anomalies when calculating restorative precipitation anomalies (g/kg)
    PmE_adj,       & ! The adjustment to PminusE that will cause the salinity
                     ! to be restored toward its target value (kg/(m^2 * s))
    net_FW,        & ! The area integrated net freshwater flux into the ocean (kg/s)
    net_FW2,       & ! The area integrated net freshwater flux into the ocean (kg/s)
    work_sum,      & ! A 2-d array that is used as the work space for a global
                     ! sum, used with units of m2 or (kg/s)
    open_ocn_mask    ! a binary field indicating where ice is present based on frazil criteria

  real :: gustiness     ! unresolved gustiness that contributes to ustar (Pa)
  real :: Irho0         ! inverse of the mean density in (m^3/kg)
  real :: taux2, tauy2  ! squared wind stresses (Pa^2)
  real :: tau_mag       ! magnitude of the wind stress (Pa)
  real :: I_GEarth      ! 1.0 / G%G_Earth  (s^2/m)
  real :: Kv_rho_ice    ! (CS%kv_sea_ice / CS%density_sea_ice) ( m^5/(s*kg) )
  real :: mass_ice      ! mass of sea ice at a face (kg/m^2)
  real :: mass_eff      ! effective mass of sea ice for rigidity (kg/m^2)

  integer :: wind_stagger  ! AGRID, BGRID_NE, or CGRID_NE (integers from MOM_domains)
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isr, ier, jsr, jer
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd

  logical :: restore_salinity ! local copy of the argument restore_salt, if it
                              ! is present, or false (no restoring) otherwise.
  logical :: restore_sst      ! local copy of the argument restore_temp, if it
                              ! is present, or false (no restoring) otherwise.
  real :: delta_sss           ! temporary storage for sss diff from restoring value
  real :: delta_sst           ! temporary storage for sst diff from restoring value

  real :: C_p                 ! heat capacity of seawater ( J/(K kg) )

  call cpu_clock_begin(id_clock_forcing)

  isc_bnd = index_bounds(1) ; iec_bnd = index_bounds(2)
  jsc_bnd = index_bounds(3) ; jec_bnd = index_bounds(4)
  is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
  Isq  = G%IscB  ; Ieq  = G%IecB   ; Jsq  = G%JscB  ; Jeq  = G%JecB
  isd  = G%isd   ; ied  = G%ied    ; jsd  = G%jsd   ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB   ; JsdB = G%JsdB  ; JedB = G%JedB
  isr = is-isd+1 ; ier  = ie-isd+1 ; jsr = js-jsd+1 ; jer = je-jsd+1

  C_p                    = fluxes%C_p
  Irho0                  = 1.0/CS%Rho0
  open_ocn_mask(:,:)     = 1.0
  pme_adj(:,:)           = 0.0
  fluxes%vPrecGlobalAdj  = 0.0
  fluxes%vPrecGlobalScl  = 0.0
  fluxes%saltFluxGlobalAdj = 0.0
  fluxes%saltFluxGlobalScl = 0.0
  fluxes%netFWGlobalAdj = 0.0
  fluxes%netFWGlobalScl = 0.0

  restore_salinity = .false.
  if (present(restore_salt)) restore_salinity = restore_salt
  restore_sst = .false.
  if (present(restore_temp)) restore_sst = restore_temp

  ! allocation and initialization if this is the first time that this
  ! flux type has been used.
  if (fluxes%dt_buoy_accum < 0) then
    call allocate_forcing_type(G, fluxes, stress=.true., ustar=.true., &
                               water=.true., heat=.true.)

    call safe_alloc_ptr(fluxes%sw_vis_dir,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_vis_dif,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_nir_dir,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%sw_nir_dif,isd,ied,jsd,jed)

    call safe_alloc_ptr(fluxes%p_surf,isd,ied,jsd,jed)
    call safe_alloc_ptr(fluxes%p_surf_full,isd,ied,jsd,jed)

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
    enddo; enddo

    if (CS%rigid_sea_ice) then
      call safe_alloc_ptr(fluxes%rigidity_ice_u,IsdB,IedB,jsd,jed)
      call safe_alloc_ptr(fluxes%rigidity_ice_v,isd,ied,JsdB,JedB)
    endif

    if (restore_temp) call safe_alloc_ptr(fluxes%heat_added,isd,ied,jsd,jed)

    fluxes%dt_buoy_accum = 0.0
  endif   ! endif for allocation and initialization

  if (CS%allow_flux_adjustments) then
   fluxes%heat_added(:,:)=0.0
   fluxes%salt_flux_added(:,:)=0.0
  endif

  ! allocation and initialization on first call to this routine
  if (CS%area_surf < 0.0) then
    do j=js,je ; do i=is,ie
      work_sum(i,j) = G%areaT(i,j) * G%mask2dT(i,j)
    enddo ; enddo
    CS%area_surf = reproducing_sum(work_sum, isr, ier, jsr, jer)
  endif    ! endif for allocation and initialization

  do j=js,je ; do i=is,ie
    fluxes%salt_flux(i,j) = 0.0
    fluxes%vprec(i,j) = 0.0
  enddo ; enddo

  ! Salinity restoring logic
  if (restore_salinity) then
    call time_interp_external(CS%id_srestore,Time,data_restore)
    ! open_ocn_mask indicates where to restore salinity (1 means restore, 0 does not)
    open_ocn_mask(:,:) = 1.0
    if (CS%mask_srestore_under_ice) then ! Do not restore under sea-ice
      do j=js,je ; do i=is,ie
        if (state%SST(i,j) .le. -0.0539*state%SSS(i,j)) open_ocn_mask(i,j)=0.0
      enddo; enddo
    endif
    if (CS%salt_restore_as_sflux) then
      do j=js,je ; do i=is,ie
        delta_sss = data_restore(i,j)- state%SSS(i,j)
        delta_sss = sign(1.0,delta_sss)*min(abs(delta_sss),CS%max_delta_srestore)
        fluxes%salt_flux(i,j) = 1.e-3*G%mask2dT(i,j) * (CS%Rho0*CS%Flux_const)* &
                  (CS%basin_mask(i,j)*open_ocn_mask(i,j)) *delta_sss  ! kg Salt m-2 s-1
      enddo; enddo
      if (CS%adjust_net_srestore_to_zero) then
        if (CS%adjust_net_srestore_by_scaling) then
          call adjust_area_mean_to_zero(fluxes%salt_flux, G, fluxes%saltFluxGlobalScl)
          fluxes%saltFluxGlobalAdj = 0.
        else
          work_sum(is:ie,js:je) = G%areaT(is:ie,js:je)*fluxes%salt_flux(is:ie,js:je)
          fluxes%saltFluxGlobalAdj = reproducing_sum(work_sum(:,:), isr,ier, jsr,jer)/CS%area_surf
          fluxes%salt_flux(is:ie,js:je) = fluxes%salt_flux(is:ie,js:je) - fluxes%saltFluxGlobalAdj
        endif
      endif
      fluxes%salt_flux_added(is:ie,js:je) = fluxes%salt_flux(is:ie,js:je) ! Diagnostic
    else
      do j=js,je ; do i=is,ie
        if (G%mask2dT(i,j) > 0.5) then
          delta_sss = state%SSS(i,j) - data_restore(i,j)
          delta_sss = sign(1.0,delta_sss)*min(abs(delta_sss),CS%max_delta_srestore)
          fluxes%vprec(i,j) = (CS%basin_mask(i,j)*open_ocn_mask(i,j))* &
                      (CS%Rho0*CS%Flux_const) * &
                      delta_sss / (0.5*(state%SSS(i,j) + data_restore(i,j)))
        endif
      enddo; enddo
      if (CS%adjust_net_srestore_to_zero) then
        if (CS%adjust_net_srestore_by_scaling) then
          call adjust_area_mean_to_zero(fluxes%vprec, G, fluxes%vPrecGlobalScl)
          fluxes%vPrecGlobalAdj = 0.
        else
          work_sum(is:ie,js:je) = G%areaT(is:ie,js:je)*fluxes%vprec(is:ie,js:je)
          fluxes%vPrecGlobalAdj = reproducing_sum(work_sum(:,:), isr, ier, jsr, jer) / CS%area_surf
          do j=js,je ; do i=is,ie
            fluxes%vprec(i,j) = ( fluxes%vprec(i,j) - fluxes%vPrecGlobalAdj ) * G%mask2dT(i,j)
          enddo; enddo
        endif
      endif
    endif
  endif



  ! SST restoring logic
  if (restore_sst) then
    call time_interp_external(CS%id_trestore,Time,data_restore)
    do j=js,je ; do i=is,ie
       delta_sst = data_restore(i,j)- state%SST(i,j)
       delta_sst = sign(1.0,delta_sst)*min(abs(delta_sst),CS%max_delta_trestore)
       fluxes%heat_added(i,j) = G%mask2dT(i,j) * (CS%Rho0*fluxes%C_p) * delta_sst * CS%Flux_const   ! W m-2
    enddo; enddo
  endif

  wind_stagger = CS%wind_stagger
  if ((IOB%wind_stagger == AGRID) .or. (IOB%wind_stagger == BGRID_NE) .or. &
      (IOB%wind_stagger == CGRID_NE)) wind_stagger = IOB%wind_stagger
  if (wind_stagger == BGRID_NE) then
    ! This is necessary to fill in the halo points.
    taux_at_q(:,:) = 0.0 ; tauy_at_q(:,:) = 0.0
  endif
  if (wind_stagger == AGRID) then
    ! This is necessary to fill in the halo points.
    taux_at_h(:,:) = 0.0 ; tauy_at_h(:,:) = 0.0
  endif


  ! obtain fluxes from IOB; note the staggering of indices
  i0 = is - isc_bnd ; j0 = js - jsc_bnd
  do j=js,je ; do i=is,ie

    if (wind_stagger == BGRID_NE) then
      if (ASSOCIATED(IOB%u_flux)) taux_at_q(I,J) = IOB%u_flux(i-i0,j-j0) * CS%wind_stress_multiplier
      if (ASSOCIATED(IOB%v_flux)) tauy_at_q(I,J) = IOB%v_flux(i-i0,j-j0) * CS%wind_stress_multiplier
    elseif (wind_stagger == AGRID) then
      if (ASSOCIATED(IOB%u_flux)) taux_at_h(i,j) = IOB%u_flux(i-i0,j-j0) * CS%wind_stress_multiplier
      if (ASSOCIATED(IOB%v_flux)) tauy_at_h(i,j) = IOB%v_flux(i-i0,j-j0) * CS%wind_stress_multiplier
    else ! C-grid wind stresses.
      if (ASSOCIATED(IOB%u_flux)) fluxes%taux(I,j) = IOB%u_flux(i-i0,j-j0) * CS%wind_stress_multiplier
      if (ASSOCIATED(IOB%v_flux)) fluxes%tauy(i,J) = IOB%v_flux(i-i0,j-j0) * CS%wind_stress_multiplier
    endif

    if (ASSOCIATED(IOB%lprec)) &
      fluxes%lprec(i,j) =  IOB%lprec(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%fprec)) &
      fluxes%fprec(i,j) = IOB%fprec(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%q_flux)) &
      fluxes%evap(i,j) = - IOB%q_flux(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%runoff)) &
      fluxes%lrunoff(i,j) = IOB%runoff(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%calving)) &
      fluxes%frunoff(i,j) = IOB%calving(i-i0,j-j0) * G%mask2dT(i,j)

    if (((ASSOCIATED(IOB%ustar_berg) .and. (.not. ASSOCIATED(fluxes%ustar_berg)))   &
      .or. (ASSOCIATED(IOB%area_berg) .and. (.not. ASSOCIATED(fluxes%area_berg)))) &
      .or. (ASSOCIATED(IOB%mass_berg) .and. (.not. ASSOCIATED(fluxes%mass_berg)))) &
      call allocate_forcing_type(G, fluxes, iceberg=.true.)

    if (ASSOCIATED(IOB%ustar_berg)) &
      fluxes%ustar_berg(i,j) = IOB%ustar_berg(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%area_berg)) &
      fluxes%area_berg(i,j) = IOB%area_berg(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%mass_berg)) &
      fluxes%mass_berg(i,j) = IOB%mass_berg(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%runoff_hflx)) &
      fluxes%heat_content_lrunoff(i,j) = IOB%runoff_hflx(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%calving_hflx)) &
      fluxes%heat_content_frunoff(i,j) = IOB%calving_hflx(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%lw_flux)) &
      fluxes%LW(i,j) = IOB%lw_flux(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%t_flux)) &
      fluxes%sens(i,j) = - IOB%t_flux(i-i0,j-j0) * G%mask2dT(i,j)

    fluxes%latent(i,j) = 0.0
    if (ASSOCIATED(IOB%fprec)) then
      fluxes%latent(i,j)            = fluxes%latent(i,j) - IOB%fprec(i-i0,j-j0)*CS%latent_heat_fusion
      fluxes%latent_fprec_diag(i,j) = -G%mask2dT(i,j) * IOB%fprec(i-i0,j-j0)*CS%latent_heat_fusion
    endif
    if (ASSOCIATED(IOB%calving)) then
      fluxes%latent(i,j)              = fluxes%latent(i,j) - IOB%calving(i-i0,j-j0)*CS%latent_heat_fusion
      fluxes%latent_frunoff_diag(i,j) = -G%mask2dT(i,j) * IOB%calving(i-i0,j-j0)*CS%latent_heat_fusion
    endif
    if (ASSOCIATED(IOB%q_flux)) then
      fluxes%latent(i,j)           = fluxes%latent(i,j) - IOB%q_flux(i-i0,j-j0)*CS%latent_heat_vapor
      fluxes%latent_evap_diag(i,j) = -G%mask2dT(i,j) * IOB%q_flux(i-i0,j-j0)*CS%latent_heat_vapor
    endif

    fluxes%latent(i,j) = G%mask2dT(i,j) * fluxes%latent(i,j)

    if (ASSOCIATED(IOB%sw_flux_vis_dir)) &
      fluxes%sw_vis_dir(i,j) = G%mask2dT(i,j) * IOB%sw_flux_vis_dir(i-i0,j-j0)
    if (ASSOCIATED(IOB%sw_flux_vis_dif)) &
      fluxes%sw_vis_dif(i,j) = G%mask2dT(i,j) * IOB%sw_flux_vis_dif(i-i0,j-j0)
    if (ASSOCIATED(IOB%sw_flux_nir_dir)) &
      fluxes%sw_nir_dir(i,j) = G%mask2dT(i,j) * IOB%sw_flux_nir_dir(i-i0,j-j0)
    if (ASSOCIATED(IOB%sw_flux_nir_dif)) &
      fluxes%sw_nir_dif(i,j) = G%mask2dT(i,j) * IOB%sw_flux_nir_dif(i-i0,j-j0)
    fluxes%sw(i,j) = fluxes%sw_vis_dir(i,j) + fluxes%sw_vis_dif(i,j) + &
                     fluxes%sw_nir_dir(i,j) + fluxes%sw_nir_dif(i,j)

  enddo ; enddo

  ! more salt restoring logic
  if (ASSOCIATED(IOB%salt_flux)) then
    do j=js,je ; do i=is,ie
      fluxes%salt_flux(i,j)    = G%mask2dT(i,j)*(fluxes%salt_flux(i,j) - IOB%salt_flux(i-i0,j-j0))
      fluxes%salt_flux_in(i,j) = G%mask2dT(i,j)*( -IOB%salt_flux(i-i0,j-j0) )
    enddo ; enddo
  endif

!### if (associated(CS%ctrl_forcing_CSp)) then
!###   do j=js,je ; do i=is,ie
!###     SST_anom(i,j) = state%SST(i,j) - CS%T_Restore(i,j)
!###     SSS_anom(i,j) = state%SSS(i,j) - CS%S_Restore(i,j)
!###     SSS_mean(i,j) = 0.5*(state%SSS(i,j) + CS%S_Restore(i,j))
!###   enddo ; enddo
!###   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_restore, &
!###                           fluxes%vprec, day, dt, G, CS%ctrl_forcing_CSp)
!### endif

  ! adjust the NET fresh-water flux to zero, if flagged
  if (CS%adjust_net_fresh_water_to_zero) then
    do j=js,je ; do i=is,ie
      net_FW(i,j) = (((fluxes%lprec(i,j)   + fluxes%fprec(i,j)) + &
                      (fluxes%lrunoff(i,j) + fluxes%frunoff(i,j))) + &
                      (fluxes%evap(i,j)    + fluxes%vprec(i,j)) ) * G%areaT(i,j)
      !   The following contribution appears to be calculating the volume flux of sea-ice
      ! melt. This calculation is clearly WRONG if either sea-ice has variable
      ! salinity or the sea-ice is completely fresh.
      !   Bob thinks this is trying ensure the net fresh-water of the ocean + sea-ice system
      ! is constant.
      !   To do this correctly we will need a sea-ice melt field added to IOB. -AJA
      if (ASSOCIATED(IOB%salt_flux) .and. (CS%ice_salt_concentration>0.0)) &
        net_FW(i,j) = net_FW(i,j) - G%areaT(i,j) * &
                     (IOB%salt_flux(i-i0,j-j0) / CS%ice_salt_concentration)
      net_FW2(i,j) = net_FW(i,j)
    enddo ; enddo

    if (CS%adjust_net_fresh_water_by_scaling) then
      call adjust_area_mean_to_zero(net_FW2, G, fluxes%netFWGlobalScl)
      do j=js,je ; do i=is,ie
        fluxes%vprec(i,j) = fluxes%vprec(i,j) + (net_FW2(i,j) - net_FW(i,j)) * G%mask2dT(i,j)
      enddo; enddo
    else
      fluxes%netFWGlobalAdj = reproducing_sum(net_FW(:,:), isr, ier, jsr, jer) / CS%area_surf
      do j=js,je ; do i=is,ie
        fluxes%vprec(i,j) = ( fluxes%vprec(i,j) - fluxes%netFWGlobalAdj ) * G%mask2dT(i,j)
      enddo; enddo
    endif

  endif

  ! applied surface pressure from atmosphere and cryosphere
  if (ASSOCIATED(IOB%p)) then
    if (CS%max_p_surf >= 0.0) then
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * IOB%p(i-i0,j-j0)
        fluxes%p_surf(i,j) = MIN(fluxes%p_surf_full(i,j),CS%max_p_surf)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * IOB%p(i-i0,j-j0)
        fluxes%p_surf(i,j) = fluxes%p_surf_full(i,j)
      enddo ; enddo
    endif
    if (CS%use_limited_P_SSH) then
      fluxes%p_surf_SSH => fluxes%p_surf
    else
      fluxes%p_surf_SSH => fluxes%p_surf_full
    endif
  endif

  ! surface momentum stress related fields as function of staggering
  if (wind_stagger == BGRID_NE) then
    if (G%symmetric) &
      call fill_symmetric_edges(taux_at_q, tauy_at_q, G%Domain, stagger=BGRID_NE)
    call pass_vector(taux_at_q, tauy_at_q, G%Domain, stagger=BGRID_NE)

    do j=js,je ; do I=Isq,Ieq
      fluxes%taux(I,j) = 0.0
      if ((G%mask2dBu(I,J) + G%mask2dBu(I,J-1)) > 0) &
        fluxes%taux(I,j) = (G%mask2dBu(I,J)*taux_at_q(I,J) + &
                            G%mask2dBu(I,J-1)*taux_at_q(I,J-1)) / &
                           (G%mask2dBu(I,J) + G%mask2dBu(I,J-1))
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie
      fluxes%tauy(i,J) = 0.0
      if ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J)) > 0) &
        fluxes%tauy(i,J) = (G%mask2dBu(I,J)*tauy_at_q(I,J) + &
                            G%mask2dBu(I-1,J)*tauy_at_q(I-1,J)) / &
                           (G%mask2dBu(I,J) + G%mask2dBu(I-1,J))
    enddo ; enddo

    ! ustar is required for the bulk mixed layer formulation. The background value
    ! of 0.02 Pa is a relatively small value intended to give reasonable behavior
    ! in regions of very weak winds.

    do j=js,je ; do i=is,ie
      tau_mag = 0.0 ; gustiness = CS%gust_const
      if (((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + &
           (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) > 0) then
        tau_mag = sqrt(((G%mask2dBu(I,J)*(taux_at_q(I,J)**2 + tauy_at_q(I,J)**2) + &
            G%mask2dBu(I-1,J-1)*(taux_at_q(I-1,J-1)**2 + tauy_at_q(I-1,J-1)**2)) + &
           (G%mask2dBu(I,J-1)*(taux_at_q(I,J-1)**2 + tauy_at_q(I,J-1)**2) + &
            G%mask2dBu(I-1,J)*(taux_at_q(I-1,J)**2 + tauy_at_q(I-1,J)**2)) ) / &
          ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) )
        if (CS%read_gust_2d) gustiness = CS%gust(i,j)
      endif
      fluxes%ustar(i,j) = sqrt(gustiness*Irho0 + Irho0*tau_mag)
    enddo ; enddo

  elseif (wind_stagger == AGRID) then
    call pass_vector(taux_at_h, tauy_at_h, G%Domain,stagger=AGRID)

    do j=js,je ; do I=Isq,Ieq
      fluxes%taux(I,j) = 0.0
      if ((G%mask2dT(i,j) + G%mask2dT(i+1,j)) > 0) &
        fluxes%taux(I,j) = (G%mask2dT(i,j)*taux_at_h(i,j) + &
                            G%mask2dT(i+1,j)*taux_at_h(i+1,j)) / &
                           (G%mask2dT(i,j) + G%mask2dT(i+1,j))
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie
      fluxes%tauy(i,J) = 0.0
      if ((G%mask2dT(i,j) + G%mask2dT(i,j+1)) > 0) &
        fluxes%tauy(i,J) = (G%mask2dT(i,j)*tauy_at_h(i,j) + &
                            G%mask2dT(i,J+1)*tauy_at_h(i,j+1)) / &
                           (G%mask2dT(i,j) + G%mask2dT(i,j+1))
    enddo ; enddo

    do j=js,je ; do i=is,ie
      gustiness = CS%gust_const
      if (CS%read_gust_2d .and. (G%mask2dT(i,j) > 0)) gustiness = CS%gust(i,j)
      fluxes%ustar(i,j) = sqrt(gustiness*Irho0 + Irho0 * G%mask2dT(i,j) * &
                               sqrt(taux_at_h(i,j)**2 + tauy_at_h(i,j)**2))
    enddo ; enddo

  else ! C-grid wind stresses.
    if (G%symmetric) &
      call fill_symmetric_edges(fluxes%taux, fluxes%tauy, G%Domain)
    call pass_vector(fluxes%taux, fluxes%tauy, G%Domain)

    do j=js,je ; do i=is,ie
      taux2 = 0.0
      if ((G%mask2dCu(I-1,j) + G%mask2dCu(I,j)) > 0) &
        taux2 = (G%mask2dCu(I-1,j)*fluxes%taux(I-1,j)**2 + &
                 G%mask2dCu(I,j)*fluxes%taux(I,j)**2) / (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))

      tauy2 = 0.0
      if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0) &
        tauy2 = (G%mask2dCv(i,J-1)*fluxes%tauy(i,J-1)**2 + &
                 G%mask2dCv(i,J)*fluxes%tauy(i,J)**2) / (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))

      if (CS%read_gust_2d) then
        fluxes%ustar(i,j) = sqrt(CS%gust(i,j)*Irho0 + Irho0*sqrt(taux2 + tauy2))
      else
        fluxes%ustar(i,j) = sqrt(CS%gust_const*Irho0 + Irho0*sqrt(taux2 + tauy2))
      endif
    enddo ; enddo

  endif   ! endif for wind related fields


  ! sea ice related fields
  if (CS%rigid_sea_ice) then
    ! The commented out code here and in the following lines is the correct
    ! version, but the incorrect version is being retained temporarily to avoid
    ! changing answers.
    call pass_var(fluxes%p_surf_full, G%Domain)
    I_GEarth = 1.0 / G%G_Earth
    Kv_rho_ice = (CS%kv_sea_ice / CS%density_sea_ice)
    do I=isd,ied-1 ; do j=jsd,jed
      mass_ice = min(fluxes%p_surf_full(i,j), fluxes%p_surf_full(i+1,j)) * I_GEarth
      mass_eff = 0.0
      if (mass_ice > CS%rigid_sea_ice_mass) then
        mass_eff = (mass_ice - CS%rigid_sea_ice_mass) **2 / &
                   (mass_ice + CS%rigid_sea_ice_mass)
      endif
      ! CAUTION: with both rigid_sea_ice and ice shelves, we will need to make this
      ! a maximum for the second call.
      fluxes%rigidity_ice_u(I,j) = Kv_rho_ice * mass_eff
    enddo ; enddo
    do i=isd,ied ; do J=jsd,jed-1
      mass_ice = min(fluxes%p_surf_full(i,j), fluxes%p_surf_full(i,j+1)) * I_GEarth
      mass_eff = 0.0
      if (mass_ice > CS%rigid_sea_ice_mass) then
        mass_eff = (mass_ice - CS%rigid_sea_ice_mass) **2 / &
                   (mass_ice + CS%rigid_sea_ice_mass)
      endif
      fluxes%rigidity_ice_v(i,J) = Kv_rho_ice * mass_eff
    enddo ; enddo
  endif

  !   At a later time, it might prove valuable to translate this array into the
  ! index space of the ocean model, rather than leaving it in the (haloless)
  ! arrays that come in from the surface forcing.
  fluxes%tr_fluxes => IOB%fluxes

  if (CS%allow_flux_adjustments) then
    ! Apply adjustments to fluxes
    call apply_flux_adjustments(G, CS, Time, fluxes)
  endif

  ! Allow for user-written code to alter fluxes after all the above
  call user_alter_forcing(state, fluxes, Time, G, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end subroutine convert_IOB_to_fluxes

!> Adds flux adjustments obtained via data_override
!! Component name is 'OCN'
!! Available adjustments are:
!! - taux_adj (Zonal wind stress delta, positive to the east, in Pa)
!! - tauy_adj (Meridional wind stress delta, positive to the north, in Pa)
subroutine apply_flux_adjustments(G, CS, Time, fluxes)
  type(ocean_grid_type),    intent(inout) :: G !< Ocean grid structure
  type(surface_forcing_CS), pointer       :: CS !< Surface forcing control structure
  type(time_type),          intent(in)    :: Time !< Model time structure
  type(forcing), optional,  intent(inout) :: fluxes !< Surface fluxes structure
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: tempx_at_h ! Delta to zonal wind stress at h points (Pa)
  real, dimension(SZI_(G),SZJ_(G)) :: tempy_at_h ! Delta to meridional wind stress at h points (Pa)
  real, dimension(SZI_(G),SZJ_(G)) :: temp_at_h ! Fluxes at h points (W m-2 or kg m-2 s-1)

  integer :: isc, iec, jsc, jec, i, j
  real :: dLonDx, dLonDy, rDlon, cosA, sinA, zonal_tau, merid_tau
  logical :: overrode_x, overrode_y, overrode_h

  isc = G%isc; iec = G%iec
  jsc = G%jsc; jec = G%jec

  overrode_h = .false.
  call data_override('OCN', 'hflx_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      fluxes%heat_added(i,j) = fluxes%heat_added(i,j) + temp_at_h(i,j)* G%mask2dT(i,j)
    enddo; enddo
  endif

  call pass_var(fluxes%heat_added, G%Domain)

  overrode_h = .false.
  call data_override('OCN', 'sflx_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      fluxes%salt_flux_added(i,j) = fluxes%salt_flux_added(i,j) + temp_at_h(i,j)* G%mask2dT(i,j)
    enddo; enddo
  endif

  call pass_var(fluxes%salt_flux_added, G%Domain)
  overrode_h = .false.

  call data_override('OCN', 'prcme_adj', temp_at_h(isc:iec,jsc:jec), Time, override=overrode_h)

  if (overrode_h) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      fluxes%vprec(i,j) = fluxes%vprec(i,j) + temp_at_h(i,j)* G%mask2dT(i,j)
    enddo; enddo
  endif

  call pass_var(fluxes%vprec, G%Domain)


  tempx_at_h(:,:) = 0.0 ; tempy_at_h(:,:) = 0.0
  ! Either reads data or leaves contents unchanged
  overrode_x = .false. ; overrode_y = .false.
  call data_override('OCN', 'taux_adj', tempx_at_h(isc:iec,jsc:jec), Time, override=overrode_x)
  call data_override('OCN', 'tauy_adj', tempy_at_h(isc:iec,jsc:jec), Time, override=overrode_y)

  if (overrode_x .or. overrode_y) then
    if (.not. (overrode_x .and. overrode_y)) call MOM_error(FATAL,"apply_flux_adjustments: "//&
            "Both taux_adj and tauy_adj must be specified, or neither, in data_table")

    ! Rotate winds
    call pass_vector(tempx_at_h, tempy_at_h, G%Domain, To_All, AGRID)
    do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      dLonDx = G%geoLonCu(I,j) - G%geoLonCu(I-1,j)
      dLonDy = G%geoLonCv(i,J) - G%geoLonCv(i,J-1)
      rDlon = sqrt( dLonDx * dLonDx + dLonDy * dLonDy )
      if (rDlon > 0.) rDlon = 1. / rDlon
      cosA = dLonDx * rDlon
      sinA = dLonDy * rDlon
      zonal_tau = tempx_at_h(i,j)
      merid_tau = tempy_at_h(i,j)
      tempx_at_h(i,j) = cosA * zonal_tau - sinA * merid_tau
      tempy_at_h(i,j) = sinA * zonal_tau + cosA * merid_tau
    enddo ; enddo

    ! Average to C-grid locations
    do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
      fluxes%taux(I,j) = fluxes%taux(I,j) + 0.5 * ( tempx_at_h(i,j) + tempx_at_h(i+1,j) )
    enddo ; enddo

    do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
      fluxes%tauy(i,J) = fluxes%tauy(i,J) + 0.5 * ( tempy_at_h(i,j) + tempy_at_h(i,j+1) )
    enddo ; enddo
  endif ! overrode_x .or. overrode_y

end subroutine apply_flux_adjustments

subroutine forcing_save_restart(CS, G, Time, directory, time_stamped, &
                                filename_suffix)
  type(surface_forcing_CS),   pointer       :: CS
  type(ocean_grid_type),      intent(inout) :: G
  type(time_type),            intent(in)    :: Time
  character(len=*),           intent(in)    :: directory
  logical,          optional, intent(in)    :: time_stamped
  character(len=*), optional, intent(in)    :: filename_suffix
! Arguments: CS - A pointer to the control structure returned by a previous
!                 call to surface_forcing_init.
!  (in)      G - The ocean's grid structure.
!  (in)      Time - The model time at this call.  This is needed for mpp_write calls.
!  (in, opt) directory - An optional directory into which to write these restart files.
!  (in, opt) time_stamped - If true, the restart file names include
!                           a unique time stamp.  The default is false.
!  (in, opt) filename_suffix - An optional suffix (e.g., a time-stamp) to append
!                              to the restart file names.

  if (.not.associated(CS)) return
  if (.not.associated(CS%restart_CSp)) return
  call save_restart(directory, Time, G, CS%restart_CSp, time_stamped)

end subroutine forcing_save_restart

subroutine surface_forcing_init(Time, G, param_file, diag, CS, restore_salt, restore_temp)
  type(time_type),          intent(in)    :: Time
  type(ocean_grid_type),    intent(in)    :: G
  type(param_file_type),    intent(in)    :: param_file
  type(diag_ctrl), target,  intent(inout) :: diag
  type(surface_forcing_CS), pointer       :: CS
  logical, optional,        intent(in)    :: restore_salt, restore_temp
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      restore_salt - If present and true, salinity restoring will be
!                           applied in this model.
  real :: utide  ! The RMS tidal velocity, in m s-1.
  type(directories)  :: dirs
  logical            :: new_sim, iceberg_flux_diags
  type(time_type)    :: Time_frc
  character(len=200) :: TideAmp_file, gust_file, salt_file, temp_file ! Input file names.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_surface_forcing"  ! This module's name.
  character(len=48)  :: stagger
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

  call write_version_number (version)
  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")

  call get_param(param_file, mod, "INPUTDIR", CS%inputdir, &
                 "The directory in which all input files are found.", &
                 default=".")
  CS%inputdir = slasher(CS%inputdir)
  call get_param(param_file, mod, "ENABLE_THERMODYNAMICS", CS%use_temperature, &
                 "If true, Temperature and salinity are used as state \n"//&
                 "variables.", default=.true.)
  call get_param(param_file, mod, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mod, "LATENT_HEAT_FUSION", CS%latent_heat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf)
  call get_param(param_file, mod, "LATENT_HEAT_VAPORIZATION", CS%latent_heat_vapor, &
                 "The latent heat of fusion.", units="J/kg", default=hlv)
  call get_param(param_file, mod, "MAX_P_SURF", CS%max_p_surf, &
                 "The maximum surface pressure that can be exerted by the \n"//&
                 "atmosphere and floating sea-ice or ice shelves. This is \n"//&
                 "needed because the FMS coupling structure does not \n"//&
                 "limit the water that can be frozen out of the ocean and \n"//&
                 "the ice-ocean heat fluxes are treated explicitly.  No \n"//&
                 "limit is applied if a negative value is used.", units="Pa", &
                 default=-1.0)
  call get_param(param_file, mod, "ADJUST_NET_SRESTORE_TO_ZERO", &
                 CS%adjust_net_srestore_to_zero, &
                 "If true, adjusts the salinity restoring seen to zero\n"//&
                 "whether restoring is via a salt flux or virtual precip.",&
                 default=restore_salt)
  call get_param(param_file, mod, "ADJUST_NET_SRESTORE_BY_SCALING", &
                 CS%adjust_net_srestore_by_scaling, &
                 "If true, adjustments to salt restoring to achieve zero net are\n"//&
                 "made by scaling values without moving the zero contour.",&
                 default=.false.)
  call get_param(param_file, mod, "ADJUST_NET_FRESH_WATER_TO_ZERO", &
                 CS%adjust_net_fresh_water_to_zero, &
                 "If true, adjusts the net fresh-water forcing seen \n"//&
                 "by the ocean (including restoring) to zero.", default=.false.)
  call get_param(param_file, mod, "ADJUST_NET_FRESH_WATER_BY_SCALING", &
                 CS%adjust_net_fresh_water_by_scaling, &
                 "If true, adjustments to net fresh water to achieve zero net are\n"//&
                 "made by scaling values without moving the zero contour.",&
                 default=.false.)
  call get_param(param_file, mod, "ICE_SALT_CONCENTRATION", &
                 CS%ice_salt_concentration, &
                 "The assumed sea-ice salinity needed to reverse engineer the \n"//&
                 "melt flux (or ice-ocean fresh-water flux).", &
                 units="kg/kg", default=0.005)
  call get_param(param_file, mod, "USE_LIMITED_PATM_SSH", CS%use_limited_P_SSH, &
                 "If true, return the sea surface height with the \n"//&
                 "correction for the atmospheric (and sea-ice) pressure \n"//&
                 "limited by max_p_surf instead of the full atmospheric \n"//&
                 "pressure.", default=.true.)

! smg: should get_param call should be removed when have A=B code reconciled.
! this param is used to distinguish how to diagnose surface heat content from water.
  call get_param(param_file, mod, "BULKMIXEDLAYER", CS%bulkmixedlayer, &
                 default=CS%use_temperature,do_not_log=.true.)

  call get_param(param_file, mod, "WIND_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the \n"//&
                 "staggering of the input wind stress field.  Valid \n"//&
                 "values are 'A', 'B', or 'C'.", default="C")
  if (uppercase(stagger(1:1)) == 'A') then ; CS%wind_stagger = AGRID
  elseif (uppercase(stagger(1:1)) == 'B') then ; CS%wind_stagger = BGRID_NE
  elseif (uppercase(stagger(1:1)) == 'C') then ; CS%wind_stagger = CGRID_NE
  else ; call MOM_error(FATAL,"surface_forcing_init: WIND_STAGGER = "// &
                        trim(stagger)//" is invalid.") ; endif
  call get_param(param_file, mod, "WIND_STRESS_MULTIPLIER", CS%wind_stress_multiplier, &
                 "A factor multiplying the wind-stress given to the ocean by the\n"//&
                 "coupler. This is used for testing and should be =1.0 for any\n"//&
                 "production runs.", default=1.0)

  if (restore_salt) then
    call get_param(param_file, mod, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes \n"//&
                 "to the relative surface anomalies (akin to a piston \n"//&
                 "velocity).  Note the non-MKS units.", units="m day-1", &
                 fail_if_missing=.true.)
    call get_param(param_file, mod, "SALT_RESTORE_FILE", CS%salt_restore_file, &
                 "A file in which to find the surface salinity to use for restoring.", &
                 default="salt_restore.nc")
    call get_param(param_file, mod, "SALT_RESTORE_VARIABLE", CS%salt_restore_var_name, &
                 "The name of the surface salinity variable to read from "//&
                 "SALT_RESTORE_FILE for restoring salinity.", &
                 default="salt")
! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0

    call get_param(param_file, mod, "SRESTORE_AS_SFLUX", CS%salt_restore_as_sflux, &
                 "If true, the restoring of salinity is applied as a salt \n"//&
                 "flux instead of as a freshwater flux.", default=.false.)
    call get_param(param_file, mod, "MAX_DELTA_SRESTORE", CS%max_delta_srestore, &
                 "The maximum salinity difference used in restoring terms.", &
                 units="PSU or g kg-1", default=999.0)
    call get_param(param_file, mod, "MASK_SRESTORE_UNDER_ICE", &
                 CS%mask_srestore_under_ice, &
                 "If true, disables SSS restoring under sea-ice based on a frazil\n"//&
                 "criteria (SST<=Tf). Only used when RESTORE_SALINITY is True.",      &
                 default=.false.)
    call get_param(param_file, mod, "MASK_SRESTORE_MARGINAL_SEAS", &
                 CS%mask_srestore_marginal_seas, &
                 "If true, disable SSS restoring in marginal seas. Only used when\n"//&
                 "RESTORE_SALINITY is True.", default=.false.)
    call get_param(param_file, mod, "BASIN_FILE", basin_file, &
                 "A file in which to find the basin masks, in variable 'basin'.", &
                 default="basin.nc")
    basin_file = trim(CS%inputdir) // trim(basin_file)
    call safe_alloc_ptr(CS%basin_mask,isd,ied,jsd,jed) ; CS%basin_mask(:,:) = 1.0
    if (CS%mask_srestore_marginal_seas) then
      call read_data(basin_file,'basin',CS%basin_mask,domain=G%domain%mpp_domain,timelevel=1)
      do j=jsd,jed ; do i=isd,ied
        if (CS%basin_mask(i,j) >= 6.0) then ; CS%basin_mask(i,j) = 0.0
        else ; CS%basin_mask(i,j) = 1.0 ; endif
      enddo ; enddo
    endif
  endif

  if (restore_temp) then
    call get_param(param_file, mod, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes \n"//&
                 "to the relative surface anomalies (akin to a piston \n"//&
                 "velocity).  Note the non-MKS units.", units="m day-1", &
                 fail_if_missing=.true.)
    call get_param(param_file, mod, "SST_RESTORE_FILE", CS%temp_restore_file, &
                 "A file in which to find the surface temperature to use for restoring.", &
                 default="temp_restore.nc")
    call get_param(param_file, mod, "SST_RESTORE_VARIABLE", CS%temp_restore_var_name, &
                 "The name of the surface temperature variable to read from "//&
                 "SST_RESTORE_FILE for restoring sst.", &
                 default="temp")
! Convert CS%Flux_const from m day-1 to m s-1.
    CS%Flux_const = CS%Flux_const / 86400.0

    call get_param(param_file, mod, "MAX_DELTA_TRESTORE", CS%max_delta_trestore, &
                 "The maximum sst difference used in restoring terms.", &
                 units="degC ", default=999.0)

  endif

! Optionally read tidal amplitude from input file (m s-1) on model grid.
! Otherwise use default tidal amplitude for bottom frictionally-generated
! dissipation. Default cd_tides is chosen to yield approx 1 TWatt of
! work done against tides globally using OSU tidal amplitude.
  call get_param(param_file, mod, "CD_TIDES", CS%cd_tides, &
                 "The drag coefficient that applies to the tides.", &
                 units="nondim", default=1.0e-4)
  call get_param(param_file, mod, "READ_TIDEAMP", CS%read_TIDEAMP, &
                 "If true, read a file (given by TIDEAMP_FILE) containing \n"//&
                 "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)
  if (CS%read_TIDEAMP) then
    call get_param(param_file, mod, "TIDEAMP_FILE", TideAmp_file, &
                 "The path to the file containing the spatially varying \n"//&
                 "tidal amplitudes with INT_TIDE_DISSIPATION.", &
                 default="tideamp.nc")
    CS%utide=0.0
  else
    call get_param(param_file, mod, "UTIDE", CS%utide, &
                 "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
                 units="m s-1", default=0.0)
  endif

  call safe_alloc_ptr(CS%TKE_tidal,isd,ied,jsd,jed)
  call safe_alloc_ptr(CS%ustar_tidal,isd,ied,jsd,jed)

  if (CS%read_TIDEAMP) then
    TideAmp_file = trim(CS%inputdir) // trim(TideAmp_file)
    call read_data(TideAmp_file,'tideamp',CS%TKE_tidal,domain=G%domain%mpp_domain,timelevel=1)
    do j=jsd, jed; do i=isd, ied
      utide = CS%TKE_tidal(i,j)
      CS%TKE_tidal(i,j) = G%mask2dT(i,j)*CS%Rho0*CS%cd_tides*(utide*utide*utide)
      CS%ustar_tidal(i,j)=sqrt(CS%cd_tides)*utide
    enddo ; enddo
  else
    do j=jsd,jed; do i=isd,ied
      utide=CS%utide
      CS%TKE_tidal(i,j) = CS%Rho0*CS%cd_tides*(utide*utide*utide)
      CS%ustar_tidal(i,j)=sqrt(CS%cd_tides)*utide
    enddo ; enddo
  endif

  call time_interp_external_init

! Optionally read a x-y gustiness field in place of a global
! constant.

  call get_param(param_file, mod, "READ_GUST_2D", CS%read_gust_2d, &
                 "If true, use a 2-dimensional gustiness supplied from \n"//&
                 "an input file", default=.false.)
  call get_param(param_file, mod, "GUST_CONST", CS%gust_const, &
               "The background gustiness in the winds.", units="Pa", &
               default=0.02)
  if (CS%read_gust_2d) then
    call get_param(param_file, mod, "GUST_2D_FILE", gust_file, &
                 "The file in which the wind gustiness is found in \n"//&
                 "variable gustiness.")

    call safe_alloc_ptr(CS%gust,isd,ied,jsd,jed)
    gust_file = trim(CS%inputdir) // trim(gust_file)
    call read_data(gust_file,'gustiness',CS%gust,domain=G%domain%mpp_domain, &
                   timelevel=1) ! units should be Pa
  endif

! See whether sufficiently thick sea ice should be treated as rigid.
  call get_param(param_file, mod, "USE_RIGID_SEA_ICE", CS%rigid_sea_ice, &
                 "If true, sea-ice is rigid enough to exert a \n"//&
                 "nonhydrostatic pressure that resist vertical motion.", &
                 default=.false.)
  if (CS%rigid_sea_ice) then
    call get_param(param_file, mod, "SEA_ICE_MEAN_DENSITY", CS%density_sea_ice, &
                 "A typical density of sea ice, used with the kinematic \n"//&
                 "viscosity, when USE_RIGID_SEA_ICE is true.", units="kg m-3", &
                 default=900.0)
    call get_param(param_file, mod, "SEA_ICE_VISCOSITY", CS%Kv_sea_ice, &
                 "The kinematic viscosity of sufficiently thick sea ice \n"//&
                 "for use in calculating the rigidity of sea ice.", &
                 units="m2 s-1", default=1.0e9)
    call get_param(param_file, mod, "SEA_ICE_RIGID_MASS", CS%rigid_sea_ice_mass, &
                 "The mass of sea-ice per unit area at which the sea-ice \n"//&
                 "starts to exhibit rigidity", units="kg m-2", default=1000.0)
  endif

  call get_param(param_file, mod, "ALLOW_ICEBERG_FLUX_DIAGNOSTICS", iceberg_flux_diags, &
                 "If true, makes available diagnostics of fluxes from icebergs\n"//&
                 "as seen by MOM6.", default=.false.)
  call register_forcing_type_diags(Time, diag, CS%use_temperature, CS%handles, &
                                   use_berg_fluxes=iceberg_flux_diags)

  call get_param(param_file, mod, "ALLOW_FLUX_ADJUSTMENTS", CS%allow_flux_adjustments, &
                 "If true, allows flux adjustments to specified via the \n"//&
                 "data_table using the component name 'OCN'.", default=.false.)
  if (CS%allow_flux_adjustments) then
    call data_override_init(Ocean_domain_in=G%Domain%mpp_domain)
  endif

  if (present(restore_salt)) then ; if (restore_salt) then
    salt_file = trim(CS%inputdir) // trim(CS%salt_restore_file)
    CS%id_srestore = init_external_field(salt_file, CS%salt_restore_var_name, domain=G%Domain%mpp_domain)
  endif ; endif

  if (present(restore_temp)) then ; if (restore_temp) then
    temp_file = trim(CS%inputdir) // trim(CS%temp_restore_file)
    CS%id_trestore = init_external_field(temp_file, CS%temp_restore_var_name, domain=G%Domain%mpp_domain)
  endif ; endif

  ! Set up any restart fields associated with the forcing.
  call restart_init(param_file, CS%restart_CSp, "MOM_forcing.res")
!###  call register_ctrl_forcing_restarts(G, param_file, CS%ctrl_forcing_CSp, &
!###                                      CS%restart_CSp)
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

!###  call controlled_forcing_init(Time, G, param_file, diag, CS%ctrl_forcing_CSp)

  call user_revise_forcing_init(param_file, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end subroutine surface_forcing_init

subroutine surface_forcing_end(CS, fluxes)
  type(surface_forcing_CS), pointer       :: CS
  type(forcing), optional,  intent(inout) :: fluxes
! Arguments:  CS - A pointer to the control structure returned by a previous
!                  call to surface_forcing_init, it will be deallocated here.
!  (inout)    fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.

  if (present(fluxes)) call deallocate_forcing_type(fluxes)

!###  call controlled_forcing_end(CS%ctrl_forcing_CSp)

  if (associated(CS)) deallocate(CS)
  CS => NULL()

end subroutine surface_forcing_end

subroutine ice_ocn_bnd_type_chksum(id, timestep, iobt)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_ocean_boundary_type), intent(in) :: iobt
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
    if (ASSOCIATED(iobt%ustar_berg)) &
      write(outunit,100) 'iobt%ustar_berg     ', mpp_chksum( iobt%ustar_berg     )
    if (ASSOCIATED(iobt%area_berg)) &
      write(outunit,100) 'iobt%area_berg      ', mpp_chksum( iobt%area_berg      )
    if (ASSOCIATED(iobt%mass_berg)) &
      write(outunit,100) 'iobt%mass_berg      ', mpp_chksum( iobt%mass_berg      )

100 FORMAT("   CHECKSUM::",A20," = ",Z20)
    do n = 1, iobt%fluxes%num_bcs  !{
       do m = 1, iobt%fluxes%bc(n)%num_fields  !{
          write(outunit,101) 'iobt%',trim(iobt%fluxes%bc(n)%name), &
               trim(iobt%fluxes%bc(n)%field(m)%name), &
               mpp_chksum(iobt%fluxes%bc(n)%field(m)%values)
       enddo  !} m
    enddo  !} n
101 FORMAT("   CHECKSUM::",A6,a,'%',a," = ",Z20)

end subroutine ice_ocn_bnd_type_chksum

end module MOM_surface_forcing
