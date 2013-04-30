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
use MOM_coms, only : reproducing_sum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_SUBCOMPONENT
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ptrs
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_domains, only : pass_vector, pass_var, global_field_sum, BITWISE_EXACT_SUM
use MOM_error_handler, only : MOM_error, WARNING, FATAL, is_root_pe, MOM_mesg
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_get_input, only : Get_MOM_Input, directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : slasher, write_version_number
use MOM_restart, only : register_restart_field, restart_init, MOM_restart_CS
use MOM_restart, only : restart_init_end, save_restart, restore_state
use MOM_variables, only : surface
use user_revise_forcing, only : user_alter_forcing, user_revise_forcing_init, &
                                user_revise_forcing_CS
!   Forcing is a structure containing pointers to the forcing fields
! which may be used to drive MOM.  All fluxes are positive downward.
!   Surface is a structure containing pointers to various fields that
! may be used describe the surface state of MOM.
!   ice_ocean_boundary_type is a structure corresponding to forcing, but with
! the elements, units, and conventions that exactly conform to the use for
! MOM-based coupled models.

use coupler_types_mod, only : coupler_2d_bc_type
use time_interp_external_mod, only : init_external_field, time_interp_external, &
                                     time_interp_external_init
use fms_mod, only : read_data, stdout
use mpp_mod, only : mpp_chksum

implicit none ; private

#include <MOM_memory.h>

public convert_IOB_to_fluxes, surface_forcing_init, average_forcing, ice_ocn_bnd_type_chksum
public forcing_save_restart

type, public :: surface_forcing_CS ; private
  integer :: wind_stagger    !   A_GRID, B_GRID, or C_GRID (integer module
                             ! parameters defined below) to indicate the
                             ! staggering of the winds that are being
                             ! provided in calls to update_ocean_model.
  logical :: use_temperature !   If true, temperature and salinity are used as
                             ! state variables.
  real :: Rho0               !   The density used in the Boussinesq
                             ! approximation, in kg m-3.
  real :: max_p_surf         !   The maximum surface pressure that can be
                             ! exerted by the atmosphere and floating sea-ice,
                             ! in Pa.  This is needed because the FMS coupling
                             ! structure does not limit the water that can be
                             ! frozen out of the ocean and the ice-ocean heat
                             ! fluxes are treated explicitly.
  logical :: use_limited_P_SSH ! If true, return the the sea surface height with
                             ! the correction for the atmospheric (and sea-ice)
                             ! pressure limited by max_p_surf instead of the
                             ! full atmospheric pressure.  The default is true.
  real :: Flux_const         !   The restoring rate at the surface, in m s-1.
  real :: gust_const         !   A constant unresolved background gustiness
                             ! that contributes to ustar, in Pa. 
  real :: area_surf          ! The total ocean surface area, in m2.
  logical :: read_gust_2d    !   If true, use a 2-dimensional gustiness supplied
                             ! from an input file.
  real, pointer, dimension(:,:) :: &
    TKE_tidal => NULL(), &   !   The turbulent kinetic energy introduced to the
                             ! bottom boundary layer by drag on the tidal flows,
                             ! in W m-2.
    gust => NULL(), &        !   A spatially varying unresolved background
                             ! gustiness that contributes to ustar, in Pa.
                             ! gust is used when read_gust_2d is true.
    ustar_tidal => NULL()    !   A tidal contribution to the bottom friction
                             ! velocity, in m s-1.
  real :: cd_tides           !   The drag coefficient that applies to the tides,
                             ! nondimensional.
  real :: utide              !   A constant tidal velocity to use if read_tideamp
                             ! is false, in m s-1.
  logical :: rigid_sea_ice   !   If true, sea-ice exerts a rigidity that acts
                             ! to damp surface deflections (especially surface
                             ! gravity waves).  The default is false.
  real    :: Kv_sea_ice      !   A viscosity in sea-ice that resists sheared
                             ! vertical motions, in m2 s-1.                           
  real    :: density_sea_ice !   A typical density of sea-ice, in kg m-3.  THis
                             ! is only used to convert the ice pressure into
                             ! appropriate units for use with Kv_sea_ice.
  real    :: rigid_sea_ice_mass ! A mass per unit area of sea-ice beyond which
                             ! sea-ice viscosity becomes effective, in kg m-2,
                             ! typically of order 1000 kg m-2.       
  logical :: read_tideamp    !   If true, the spatially varying tidal amplitude
                             ! is read from a file.
  logical :: salt_restore_as_sflux ! If true, the restoring of salinity is applied
                             ! as a salt flux in place of a freshwater flux.
  logical :: adjust_net_fresh_water_to_zero ! Adjust the NET fresh-water forcing 
                             ! to zero, including restoring
  logical :: mask_srestore_under_ice  ! If true, use an ice mask defined by
                             ! frazil criteria for salinity restoring.
  real    :: ice_salt_concentration  ! salt concentration for sea ice in kg kg-1
  logical :: mask_srestore_marginal_seas   ! if true, then mask sss restoring in marginal seas
  real :: max_delta_srestore ! maximum delta salinity used for restoring (duplicates mom4 option)
  real, pointer, dimension(:,:) :: basin_mask => NULL() ! mask for sss restoring
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
                             ! ocean diagnostic fields and control variables.
  character(len=200) :: inputdir ! The directory where NetCDF input files are.

  logical :: first_call = .true. ! True if convert_IOB_to_fluxes has not been
                                 ! called yet.

  integer :: id_taux = -1, id_tauy = -1, id_ustar = -1
  integer :: id_PminusE = -1, id_evap = -1, id_precip = -1
  integer :: id_liq_precip = -1, id_froz_precip = -1, id_virt_precip = -1
  integer :: id_liq_runoff = -1, id_froz_runoff = -1
  integer :: id_runoff_hflx = -1, id_calving_hflx = -1
  integer :: id_Net_Heating = -1, id_sw = -1, id_LwLatSens = -1, id_buoy = -1
  integer :: id_LW = -1, id_lat = -1, id_sens = -1
  integer :: id_psurf = -1, id_saltflux = -1, id_TKE_tidal = -1
  integer :: id_srestore = -1  ! An id number for time_interp_external.
!###  type(ctrl_forcing_CS), pointer :: ctrl_forcing_CSp => NULL()
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()
  type(user_revise_forcing_CS), pointer :: urf_CS => NULL() 
end type surface_forcing_CS

type, public :: ice_ocean_boundary_type
  real, pointer, dimension(:,:) :: u_flux =>NULL()   ! i-direction wind stress (Pa)
  real, pointer, dimension(:,:) :: v_flux =>NULL()   ! j-direction wind stress (Pa)
  real, pointer, dimension(:,:) :: t_flux =>NULL()   ! sensible heat flux (W/m2)
  real, pointer, dimension(:,:) :: q_flux =>NULL()   ! specific humidity flux (kg/m2/s)
  real, pointer, dimension(:,:) :: salt_flux =>NULL()! salt flux (kg/m2/s)
  real, pointer, dimension(:,:) :: lw_flux =>NULL()  ! long wave radiation (w/m2)
  real, pointer, dimension(:,:) :: sw_flux_vis_dir => NULL() ! direct visible sw radiation (w/m2)
  real, pointer, dimension(:,:) :: sw_flux_vis_dif => NULL() ! diffuse visible sw radiation (w/m2)
  real, pointer, dimension(:,:) :: sw_flux_nir_dir => NULL() ! direct Near InfraRed sw radiation (w/m2)
  real, pointer, dimension(:,:) :: sw_flux_nir_dif => NULL() ! diffuse Near InfraRed sw radiation (w/m2)
  real, pointer, dimension(:,:) :: lprec =>NULL()    ! mass flux of liquid precip (kg/m2/s)
  real, pointer, dimension(:,:) :: fprec =>NULL()    ! mass flux of frozen precip (kg/m2/s)
  real, pointer, dimension(:,:) :: runoff =>NULL()   ! mass flux of liquid runoff (kg/m2/s)
  real, pointer, dimension(:,:) :: calving =>NULL()  ! mass flux of frozen runoff (kg/m2/s)
  real, pointer, dimension(:,:) :: runoff_hflx =>NULL() ! heat flux associated with liquid runoff (w/m2)
  real, pointer, dimension(:,:) :: calving_hflx =>NULL()! heat flux associated with frozen runoff (w/m2)
  real, pointer, dimension(:,:) :: p =>NULL()        ! pressure of overlying ice and atmosphere
                                                     ! on ocean surface (Pa)
  real, pointer, dimension(:,:) :: mi =>NULL()       ! mass of ice (kg/m2)
  integer :: xtype                                   ! REGRID, REDIST or DIRECT
  type(coupler_2d_bc_type)      :: fluxes            ! A structure that may contain an
                                                     ! array of named fields used for
                                                     ! passive tracer fluxes.
end type ice_ocean_boundary_type

integer, parameter :: A_GRID = 1, B_GRID = 2, C_GRID = 3
integer :: id_clock_forcing

contains

subroutine convert_IOB_to_fluxes(IOB, fluxes, index_bounds, Time, G, CS, state, restore_salt)
  use MOM_constants, only : hlv, hlf
  use MOM_domains, only : BGRID_NE
  type(ice_ocean_boundary_type), intent(in), target :: IOB
  type(forcing),              intent(inout) :: fluxes
  integer, dimension(4),      intent(in)    :: index_bounds
  type(time_type),            intent(in)    :: Time
  type(ocean_grid_type),      intent(inout) :: G
  type(surface_forcing_CS),   pointer       :: CS
  type(surface),              intent(in)    :: state
  logical, optional,          intent(in)    :: restore_salt
! This subroutine translates the Ice_ocean_boundary_type into a
! MOM forcing type, including changes of units, sign conventions,
! and puting the fields into arrays with MOM-standard halos.
! Arguments: IOB - the ice-ocean boundary type, containing all the fluxes used
!                  to drive the ocean in a coupled model.
!  (out)     fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      index_bounds - the i- and j- size of the arrays in IOB.
!  (in)      Time - The time of the fluxes, used for interpolating the salinity
!                   to the right time, when it is being restored.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to surface_forcing_init.
!  (in)      state - A structure containing fields that describe the
!                    surface state of the ocean.
!  (in)      restore_salt - if true, salinity is restored to a target value.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    taux_at_q, &     ! Zonal wind stresses at q points in Pa.
    tauy_at_q, &     ! Meridional wind stresses at q points in Pa.
    data_srestore, & ! The surface salinity toward which to restore, in PSU.
    SST_anom, &      ! Instantaneous sea surface temperature anomalies from a
                     ! target (observed) value, in deg C.
    SSS_anom, &      ! Instantaneous sea surface salinity anomalies from a target
                     ! (observed) value, in g kg-1.
    SSS_mean, &      ! A (mean?) salinity about which to normalize local salinity
                     ! anomalies when calculating restorative precipitation
                     ! anomalies, in g kg-1.
    PmE_adj, &       ! The adjustment to PminusE that will cause the salinity
                     ! to be restored toward its target value, in kg m-2 s-1.
    net_FW, &        ! The area integrated net freshwater flux into the ocean, in kg s-1.
    work_sum, &      ! A 2-d array that is used as the work space for a global
                     ! sum, used with units of m2 or kg s-1.
    open_ocn_mask    ! a binary field  indicating where ice is present based on
                     ! frazil criteria (per MOM4)
  real :: gustiness     ! The unresolved gustiness that contributes to ustar in Pa.
  real :: PmE_adj_total ! The globally area integrated PmE_adj, in kg s-1.
  real :: net_FW_avg    ! The globally averaged net fresh water input to the 
                        ! ocean/sea-ice system, in kg m-2 s-1.
  real :: Sflux_adj_total  ! The globally area integrated salt flux adjustment, in kg s-1.  
  real :: Irho0         ! The inverse of the mean density in m3 kg-1.
  real :: taux2, tauy2  ! The squared wind stresses in Pa2.
  real :: tau_mag       ! The magnitude of the wind stress, in Pa.
  real :: I_GEarth      ! 1.0 / G%G_Earth, in s2 m-1.
  real :: Kv_rho_ice    ! (CS%kv_sea_ice / CS%density_sea_ice), in m5 s-1 kg-1.
  real :: mass_ice      ! The mass of sea ice at a face, in kg m-2.
  real :: mass_eff      ! The effective mass of sea ice for rigidity, in kg m-2.
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq, i0, j0
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isr, ier, jsr, jer
  integer :: isc_bnd, iec_bnd, jsc_bnd, jec_bnd
  logical :: restore_salinity ! A local copy of the argument restore_salt, if it
                              ! is present, or false (no restoring) otherwise.
  real :: delta_sss   ! temporary storage for sss difference from restoring value
  
  call cpu_clock_begin(id_clock_forcing)

  isc_bnd = index_bounds(1) ; iec_bnd = index_bounds(2)
  jsc_bnd = index_bounds(3) ; jec_bnd = index_bounds(4)
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  isr = is-isd+1 ; ier = ie-isd+1 ; jsr = js-jsd+1 ; jer = je-jsd+1

  Irho0 = 1.0/CS%Rho0

  open_ocn_mask(:,:) = 1.0
  pme_adj(:,:) = 0.0
  PmE_adj_total = 0.0
  Sflux_adj_total = 0.0
  
  restore_salinity = .false.
  if (present(restore_salt)) restore_salinity = restore_salt

  if (CS%first_call) then
    call safe_alloc_ptr(fluxes%taux,IsdB,IedB,jsd,jed)      ; fluxes%taux(:,:) = 0.0
    call safe_alloc_ptr(fluxes%tauy,isd,ied,JsdB,JedB)      ; fluxes%tauy(:,:) = 0.0
    call safe_alloc_ptr(fluxes%ustar,isd,ied,jsd,jed)       ; fluxes%ustar(:,:) = 0.0
    call safe_alloc_ptr(fluxes%evap,isd,ied,jsd,jed)        ; fluxes%evap(:,:) = 0.0
    call safe_alloc_ptr(fluxes%liq_precip,isd,ied,jsd,jed)  ; fluxes%liq_precip(:,:) = 0.0
    call safe_alloc_ptr(fluxes%froz_precip,isd,ied,jsd,jed) ; fluxes%froz_precip(:,:) = 0.0
    call safe_alloc_ptr(fluxes%virt_precip,isd,ied,jsd,jed) ; fluxes%virt_precip(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw,isd,ied,jsd,jed)          ; fluxes%sw(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw_vis_dir,isd,ied,jsd,jed)  ; fluxes%sw_vis_dir(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw_vis_dif,isd,ied,jsd,jed)  ; fluxes%sw_vis_dif(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw_nir_dir,isd,ied,jsd,jed)  ; fluxes%sw_nir_dir(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sw_nir_dif,isd,ied,jsd,jed)  ; fluxes%sw_nir_dif(:,:) = 0.0
    call safe_alloc_ptr(fluxes%lw,isd,ied,jsd,jed)          ; fluxes%lw(:,:) = 0.0
    call safe_alloc_ptr(fluxes%latent,isd,ied,jsd,jed)      ; fluxes%latent(:,:) = 0.0
    call safe_alloc_ptr(fluxes%sens,isd,ied,jsd,jed)        ; fluxes%sens(:,:) = 0.0
    call safe_alloc_ptr(fluxes%p_surf,isd,ied,jsd,jed)      ; fluxes%p_surf(:,:) = 0.0
    call safe_alloc_ptr(fluxes%p_surf_full,isd,ied,jsd,jed) ; fluxes%p_surf_full(:,:) = 0.0
    call safe_alloc_ptr(fluxes%salt_flux,isd,ied,jsd,jed)   ; fluxes%salt_flux(:,:) = 0.0
    call safe_alloc_ptr(fluxes%TKE_tidal,isd,ied,jsd,jed)   ; fluxes%TKE_tidal(:,:) = 0.0
    call safe_alloc_ptr(fluxes%ustar_tidal,isd,ied,jsd,jed) ; fluxes%ustar_tidal(:,:) = 0.0
    call safe_alloc_ptr(fluxes%liq_runoff,isd,ied,jsd,jed)  ; fluxes%liq_runoff(:,:) = 0.0        
    call safe_alloc_ptr(fluxes%froz_runoff,isd,ied,jsd,jed) ; fluxes%froz_runoff(:,:) = 0.0        
    if (ASSOCIATED(IOB%calving_hflx)) then
      call safe_alloc_ptr(fluxes%calving_hflx,isd,ied,jsd,jed) ; fluxes%calving_hflx(:,:) = 0.0        
    endif
    if (ASSOCIATED(IOB%runoff_hflx)) then
      call safe_alloc_ptr(fluxes%runoff_hflx,isd,ied,jsd,jed) ; fluxes%runoff_hflx(:,:) = 0.0        
    endif
    if (CS%rigid_sea_ice) then
      call safe_alloc_ptr(fluxes%rigidity_ice_u,IsdB,IedB,jsd,jed)
      call safe_alloc_ptr(fluxes%rigidity_ice_v,isd,ied,JsdB,JedB)
      fluxes%rigidity_ice_u(:,:) = 0.0 ; fluxes%rigidity_ice_v(:,:) = 0.0  
    endif

    CS%first_call = .false.

    do j=js,je ; do i=is,ie
      work_sum(i,j) = G%areaT(i,j) * G%mask2dT(i,j)
    enddo ; enddo
    CS%area_surf = reproducing_sum(work_sum, isr, ier, jsr, jer)

    do j=js-2,je+2 ; do i=is-2,ie+2
      fluxes%TKE_tidal(i,j) = CS%TKE_tidal(i,j)
      fluxes%ustar_tidal(i,j) = CS%ustar_tidal(i,j)
    enddo; enddo

  endif

  if (restore_salinity) then
      call time_interp_external(CS%id_srestore,Time,data_srestore)
      open_ocn_mask = 1.0
      if (CS%mask_srestore_under_ice) then
        do j=js,je ; do i=is,ie
          if (state%SST(i,j) .le. -0.0539*state%SSS(i,j)) open_ocn_mask(i,j)=0.0
        enddo; enddo
      endif
      if (CS%salt_restore_as_sflux) then
        do j=js,je ; do i=is,ie
          delta_sss = data_srestore(i,j)- state%SSS(i,j)
          delta_sss = sign(1.0,delta_sss)*min(abs(delta_sss),CS%max_delta_srestore)
          fluxes%salt_flux(i,j) = 1.e-3*G%mask2dT(i,j) * (CS%Rho0*CS%Flux_const)* &
                    (CS%basin_mask(i,j)*open_ocn_mask(i,j)) *delta_sss  ! kg Salt m-2 s-1
          work_sum(i,j) = G%areaT(i,j)*fluxes%salt_flux(i,j)
        enddo; enddo
        Sflux_adj_total = reproducing_sum(work_sum(:,:), isr,ier, jsr,jer) / &
                          CS%area_surf
      else
        do j=js,je ; do i=is,ie
          if (G%mask2dT(i,j) > 0.5) then
            delta_sss = state%SSS(i,j) - data_srestore(i,j)
            delta_sss = sign(1.0,delta_sss)*min(abs(delta_sss),CS%max_delta_srestore)
            pme_adj(i,j) = (CS%basin_mask(i,j)*open_ocn_mask(i,j))* &
                      (CS%Rho0*CS%Flux_const) * &
                      delta_sss / (0.5*(state%SSS(i,j) + data_srestore(i,j)))
          else
            pme_adj(i,j) = 0.0
          endif
          work_sum(i,j) = G%areaT(i,j) * pme_adj(i,j)
        enddo; enddo
        PmE_adj_total = reproducing_sum(work_sum(:,:), isr, ier, jsr, jer) / &
                           CS%area_surf
        ! Note that when CS%adjust_net_fresh_water_to_zero is true, this adjustment
        ! of the net salt-restoring to zero is redundant but has been left here
        ! for backward compatibility. See section below where
        ! CS%adjust_net_fresh_water_to_zero is tested to be true.
      endif
  endif

  if (CS%wind_stagger == B_GRID) then
    ! This is necessary to fill in the halo points.
    taux_at_q(:,:) = 0.0 ; tauy_at_q(:,:) = 0.0
  endif

  i0 = is - isc_bnd ; j0 = js - jsc_bnd
  do j=js,je ; do i=is,ie
    if (CS%wind_stagger == B_GRID) then
      if (ASSOCIATED(IOB%u_flux)) taux_at_q(i,j) = IOB%u_flux(i-i0,j-j0)
      if (ASSOCIATED(IOB%v_flux)) tauy_at_q(i,j) = IOB%v_flux(i-i0,j-j0)
    else ! C-grid wind stresses.
      if (ASSOCIATED(IOB%u_flux)) fluxes%taux(i,j) = IOB%u_flux(i-i0,j-j0)
      if (ASSOCIATED(IOB%v_flux)) fluxes%tauy(i,j) = IOB%v_flux(i-i0,j-j0)
    endif

    if (ASSOCIATED(IOB%lprec)) &
      fluxes%liq_precip(i,j) =  IOB%lprec(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%fprec)) &
      fluxes%froz_precip(i,j) = IOB%fprec(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%q_flux)) &
      fluxes%evap(i,j) = - IOB%q_flux(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%runoff)) &
      fluxes%liq_runoff(i,j) = IOB%runoff(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%calving)) &
      fluxes%froz_runoff(i,j) = IOB%calving(i-i0,j-j0) * G%mask2dT(i,j)

    if (restore_salinity) &
      fluxes%virt_precip(i,j) = (pme_adj(i,j) - PmE_adj_total) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%calving_hflx)) &
      fluxes%calving_hflx(i,j) = IOB%calving_hflx(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%runoff_hflx)) &
      fluxes%runoff_hflx(i,j) = IOB%runoff_hflx(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%lw_flux)) &
      fluxes%LW(i,j) = IOB%lw_flux(i-i0,j-j0) * G%mask2dT(i,j)

    if (ASSOCIATED(IOB%t_flux)) &
      fluxes%sens(i,j) = - IOB%t_flux(i-i0,j-j0) * G%mask2dT(i,j)

    fluxes%latent(i,j) = 0.0
    if (ASSOCIATED(IOB%fprec)) &
      fluxes%latent(i,j) = fluxes%latent(i,j) - IOB%fprec(i-i0,j-j0)*hlf
    if (ASSOCIATED(IOB%calving)) &
      fluxes%latent(i,j) = fluxes%latent(i,j) - IOB%calving(i-i0,j-j0)*hlf
    if (ASSOCIATED(IOB%q_flux)) &
      fluxes%latent(i,j) = fluxes%latent(i,j) - IOB%q_flux(i-i0,j-j0)*hlv
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


    if (restore_salinity .and. CS%salt_restore_as_sflux) then
      fluxes%salt_flux(i,j) = G%mask2dT(i,j)*(fluxes%salt_flux(i,j)-Sflux_adj_total)
    else
      fluxes%salt_flux(i,j) = 0.0
    endif

    if (ASSOCIATED(IOB%salt_flux)) then
      fluxes%salt_flux(i,j) = G%mask2dT(i,j)*(fluxes%salt_flux(i,j) - IOB%salt_flux(i-i0,j-j0))
    endif
  enddo ; enddo

!### if (associated(CS%ctrl_forcing_CSp)) then
!###   do j=js,je ; do i=is,ie
!###     SST_anom(i,j) = state%SST(i,j) - CS%T_Restore(i,j)
!###     SSS_anom(i,j) = state%SSS(i,j) - CS%S_Restore(i,j)
!###     SSS_mean(i,j) = 0.5*(state%SSS(i,j) + CS%S_Restore(i,j))
!###   enddo ; enddo
!###   call apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, fluxes%heat_restore, &
!###                           fluxes%virt_precip, day, dt, G, CS%ctrl_forcing_CSp)
!### endif

  ! Adjust the NET fresh-water flux to zero, if flagged
  if (CS%adjust_net_fresh_water_to_zero) then
    do j=js,je ; do i=is,ie
      net_FW(i,j) = (((fluxes%liq_precip(i,j) + fluxes%froz_precip(i,j)) + &
                      (fluxes%liq_runoff(i,j) + fluxes%froz_runoff(i,j))) + &
                      (fluxes%evap(i,j) + fluxes%virt_precip(i,j)) ) * G%areaT(i,j)
      if (ASSOCIATED(IOB%salt_flux) .and. (CS%ice_salt_concentration>0.0)) &
        net_FW(i,j) = net_FW(i,j) - G%areaT(i,j) * &
                     (IOB%salt_flux(i-i0,j-j0) / CS%ice_salt_concentration)
     
    enddo ; enddo
  
    net_FW_avg = reproducing_sum(net_FW(:,:), isr, ier, jsr, jer) / &
                 CS%area_surf
    do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j) > 0.5) fluxes%virt_precip(i,j) = &
             fluxes%virt_precip(i,j) - net_FW_avg
    enddo; enddo
  endif

  if (ASSOCIATED(IOB%p) .and. (CS%max_p_surf >= 0.0)) then
    if (CS%use_limited_P_SSH) then
      do j=js,je ; do i=is,ie
        fluxes%p_surf(i,j) = G%mask2dT(i,j) * MIN(IOB%p(i-i0,j-j0),CS%max_p_surf)
        fluxes%p_surf_full(i,j) = fluxes%p_surf(i,j)
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * IOB%p(i-i0,j-j0)
        fluxes%p_surf(i,j) = MIN(fluxes%p_surf_full(i,j),CS%max_p_surf)
      enddo ; enddo
    endif
  elseif (ASSOCIATED(IOB%p)) then
    do j=js,je ; do i=is,ie
      fluxes%p_surf_full(i,j) = G%mask2dT(i,j) * IOB%p(i-i0,j-j0)
      fluxes%p_surf(i,j) = fluxes%p_surf_full(i,j)
    enddo ; enddo
  endif

  if (CS%wind_stagger == B_GRID) then
    call pass_vector(taux_at_q,tauy_at_q,G%Domain,stagger=BGRID_NE)

    do j=js,je ; do I=Isq,Ieq
      fluxes%taux(I,j) = 0.0
      If ((G%mask2dBu(I,J) + G%mask2dBu(I,J-1)) > 0) &
        fluxes%taux(I,j) = (G%mask2dBu(I,J)*taux_at_q(I,J) + G%mask2dBu(I,J-1)*taux_at_q(I,J-1)) / &
            (G%mask2dBu(I,J) + G%mask2dBu(I,J-1))
    enddo ; enddo

    do J=Jsq,Jeq ; do i=is,ie
      fluxes%tauy(i,J) = 0.0
      if ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J)) > 0) &
        fluxes%tauy(i,J) = (G%mask2dBu(I,J)*tauy_at_q(I,J) + G%mask2dBu(I-1,J)*tauy_at_q(I-1,J)) / &
            (G%mask2dBu(I,J) + G%mask2dBu(I-1,J))
    enddo ; enddo

    ! ustar is required for MOM's mixed layer formulation.  The background value
    ! of 0.02 Pa is a relatively small value intended to give reasonable behavior
    ! in regions of very weak winds.

    do j=js,je ; do i=is,ie
      tau_mag = 0.0 ; gustiness = CS%gust_const
      if (((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) > 0) then
        tau_mag = sqrt(((G%mask2dBu(I,J)*(taux_at_q(I,J)**2 + tauy_at_q(I,J)**2) + &
            G%mask2dBu(I-1,J-1)*(taux_at_q(I-1,J-1)**2 + tauy_at_q(I-1,J-1)**2)) + &
           (G%mask2dBu(I,J-1)*(taux_at_q(I,J-1)**2 + tauy_at_q(I,J-1)**2) + &
            G%mask2dBu(I-1,J)*(taux_at_q(I-1,J)**2 + tauy_at_q(I-1,J)**2)) ) / &
          ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) )
        if (CS%read_gust_2d) gustiness = CS%gust(i,j)
      endif
      fluxes%ustar(i,j) = sqrt(gustiness*Irho0 + Irho0*tau_mag)
    enddo ; enddo
  else ! C-grid wind stresses.
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
  endif

  if (CS%rigid_sea_ice) then
    call pass_var(fluxes%p_surf_full, G%Domain)
    I_GEarth = 1.0 / G%G_Earth
    Kv_rho_ice = (CS%kv_sea_ice / CS%density_sea_ice)
    do I=isd,ied-1 ; do j=isd,jed
      mass_ice = min(fluxes%p_surf_full(i,j), fluxes%p_surf_full(i+1,j)) * I_GEarth
      mass_eff = 0.0
      if (mass_ice > CS%rigid_sea_ice_mass) then
        mass_eff = (mass_ice - CS%rigid_sea_ice_mass) **2 / &
                   (mass_ice + CS%rigid_sea_ice_mass)
        ! Alistar thinks that for simplicity this should be
        !   mass_eff = (mass_ice - CS%rigid_sea_ice_mass)
        ! but Bob thinks it should vary smoothly, like (m-m1)^2/(m+m1) does.
      endif
      ! CAUTION: with both rigid_sea_ice and ice shelves, we will need to make this
      ! a maximum for the second call.
      fluxes%rigidity_ice_u(I,j) = Kv_rho_ice * mass_eff
    enddo ; enddo
    do i=isd,ied ; do J=isd,jed-1
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

  ! Allow for user-written code to alter fluxes after all the above
  call user_alter_forcing(state, fluxes, Time, G, CS%urf_CS)

  call cpu_clock_end(id_clock_forcing)
end subroutine convert_IOB_to_fluxes

subroutine average_forcing(fluxes, dt, G, CS)
  type(forcing),         intent(in) :: fluxes
  real,                  intent(in) :: dt
  type(ocean_grid_type), intent(in) :: G
  type(surface_forcing_CS), pointer :: CS
!   This subroutine offers forcing fields for time averaging.  These
! fields must first be registered in surface_forcing_init (below).
! This subroutine will typically not be modified, except when new
! forcing fields are added.
!
! Arguments: fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields are unallocated.
!  (in)      dt - The amount of time over which to average.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to surface_forcing_init.

  real, dimension(SZI_(G),SZJ_(G)) :: sum

  call cpu_clock_begin(id_clock_forcing)

  if (query_averaging_enabled(CS%diag)) then
    if ((CS%id_taux > 0) .and. ASSOCIATED(fluxes%taux)) &
      call post_data(CS%id_taux, fluxes%taux, CS%diag)
    if ((CS%id_tauy > 0) .and. ASSOCIATED(fluxes%tauy)) &
      call post_data(CS%id_tauy, fluxes%tauy, CS%diag)
    if ((CS%id_ustar > 0) .and. ASSOCIATED(fluxes%ustar)) &
      call post_data(CS%id_ustar, fluxes%ustar, CS%diag)

    if (CS%id_PminusE > 0) then
      sum(:,:) = 0.0
      if (ASSOCIATED(fluxes%liq_precip)) sum(:,:) = sum(:,:)+fluxes%liq_precip(:,:)
      if (ASSOCIATED(fluxes%froz_precip)) sum(:,:) = sum(:,:)+fluxes%froz_precip(:,:)
      if (ASSOCIATED(fluxes%evap)) sum(:,:) = sum(:,:)+fluxes%evap(:,:)
      if (ASSOCIATED(fluxes%liq_runoff)) sum(:,:) = sum(:,:)+fluxes%liq_runoff(:,:)
      if (ASSOCIATED(fluxes%froz_runoff)) sum(:,:) = sum(:,:)+fluxes%froz_runoff(:,:)
      if (ASSOCIATED(fluxes%virt_precip)) sum(:,:) = sum(:,:)+fluxes%virt_precip(:,:)
      call post_data(CS%id_PminusE, sum, CS%diag)
    endif

    if ((CS%id_evap > 0) .and. ASSOCIATED(fluxes%evap)) &
      call post_data(CS%id_evap, fluxes%evap, CS%diag)
    if ((CS%id_precip > 0) .and. ASSOCIATED(fluxes%liq_precip) &
         .and. ASSOCIATED(fluxes%froz_precip)) then
      sum(:,:) = fluxes%liq_precip(:,:) + fluxes%froz_precip(:,:)
      call post_data(CS%id_precip, sum, CS%diag)
    endif

    if ((CS%id_liq_precip > 0) .and. ASSOCIATED(fluxes%liq_precip)) &
      call post_data(CS%id_liq_precip, fluxes%liq_precip, CS%diag)
    if ((CS%id_froz_precip > 0) .and. ASSOCIATED(fluxes%froz_precip)) &
      call post_data(CS%id_froz_precip, fluxes%froz_precip, CS%diag)
    if ((CS%id_virt_precip > 0) .and. ASSOCIATED(fluxes%virt_precip)) &
      call post_data(CS%id_virt_precip, fluxes%virt_precip, CS%diag)
    if ((CS%id_liq_runoff > 0) .and. ASSOCIATED(fluxes%liq_runoff)) &
      call post_data(CS%id_liq_runoff, fluxes%liq_runoff, CS%diag)
    if ((CS%id_froz_runoff > 0) .and. ASSOCIATED(fluxes%froz_runoff)) &
      call post_data(CS%id_froz_runoff, fluxes%froz_runoff, CS%diag)

    if ((CS%id_runoff_hflx > 0) .and. ASSOCIATED(fluxes%runoff_hflx)) &
      call post_data(CS%id_runoff_hflx, fluxes%runoff_hflx, CS%diag)
    if ((CS%id_calving_hflx > 0) .and. ASSOCIATED(fluxes%calving_hflx)) &
      call post_data(CS%id_calving_hflx, fluxes%calving_hflx, CS%diag)

    if (CS%id_Net_Heating > 0) then
      sum(:,:) = 0.0
      if (ASSOCIATED(fluxes%LW)) sum(:,:) = sum(:,:) + fluxes%LW(:,:)
      if (ASSOCIATED(fluxes%latent)) sum(:,:) = sum(:,:) + fluxes%latent(:,:)
      if (ASSOCIATED(fluxes%sens)) sum(:,:) = sum(:,:) + fluxes%sens(:,:)
      if (ASSOCIATED(fluxes%SW)) sum(:,:) = sum(:,:) + fluxes%SW(:,:)
      call post_data(CS%id_Net_Heating, sum, CS%diag)
    endif
    if ((CS%id_LwLatSens > 0) .and. ASSOCIATED(fluxes%lw) .and. &
         ASSOCIATED(fluxes%latent) .and. ASSOCIATED(fluxes%sens)) then
      sum(:,:) = (fluxes%lw(:,:) + fluxes%latent(:,:)) + fluxes%sens(:,:)
      call post_data(CS%id_LwLatSens, sum, CS%diag)
    endif

    if ((CS%id_sw > 0) .and. ASSOCIATED(fluxes%sw)) &
      call post_data(CS%id_sw, fluxes%sw, CS%diag)
    if ((CS%id_LW > 0) .and. ASSOCIATED(fluxes%LW)) &
      call post_data(CS%id_LW, fluxes%LW, CS%diag)
    if ((CS%id_lat > 0) .and. ASSOCIATED(fluxes%latent)) &
      call post_data(CS%id_lat, fluxes%latent, CS%diag)
    if ((CS%id_sens > 0) .and. ASSOCIATED(fluxes%sens)) &
      call post_data(CS%id_sens, fluxes%sens, CS%diag)

    if ((CS%id_psurf > 0) .and. ASSOCIATED(fluxes%p_surf)) &
      call post_data(CS%id_psurf, fluxes%p_surf, CS%diag)
    if ((CS%id_saltflux > 0) .and. ASSOCIATED(fluxes%salt_flux)) &
      call post_data(CS%id_saltflux, fluxes%salt_flux, CS%diag)
    if ((CS%id_TKE_tidal > 0) .and. ASSOCIATED(fluxes%TKE_tidal)) &
      call post_data(CS%id_TKE_tidal, fluxes%TKE_tidal, CS%diag)

    if ((CS%id_buoy > 0) .and. ASSOCIATED(fluxes%buoy)) &
      call post_data(CS%id_buoy, fluxes%buoy, CS%diag)
  endif

  call cpu_clock_end(id_clock_forcing)
end subroutine average_forcing

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

subroutine surface_forcing_init(Time, G, param_file, diag, CS, restore_salt)
  type(time_type),          intent(in) :: Time
  type(ocean_grid_type),    intent(in) :: G
  type(param_file_type),    intent(in) :: param_file
  type(diag_ptrs), target,  intent(in) :: diag
  type(surface_forcing_CS), pointer    :: CS
  logical, optional,       intent(in) :: restore_salt
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      restore_salt - If present and true, salinity restoring will be
!                           applied in this model.
  real :: utide  ! The RMS tidal velocity, in m s-1.
  type(directories)  :: dirs
  logical            :: new_sim
  type(time_type)    :: Time_frc
  character(len=200) :: TideAmp_file, gust_file, salt_file ! Input file names.
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'
  character(len=40)  :: mod = "MOM_surface_forcing"  ! This module's name.
  character(len=48)  :: stagger
  character(len=128) :: basin_file
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

  call write_version_number (version, tagname)
  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, tagname, "")

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
  call get_param(param_file, mod, "MAX_P_SURF", CS%max_p_surf, &
                 "The maximum surface pressure that can be exerted by the \n"//&
                 "atmosphere and floating sea-ice or ice shelves. This is \n"//&
                 "needed because the FMS coupling structure does not \n"//&
                 "limit the water that can be frozen out of the ocean and \n"//&
                 "the ice-ocean heat fluxes are treated explicitly.  No \n"//&
                 "limit is applied if a negative value is used.", units="Pa", &
                 default=-1.0)
  call get_param(param_file, mod, "ADJUST_NET_FRESH_WATER_TO_ZERO", &
                 CS%adjust_net_fresh_water_to_zero, &
                 "If true, adjusts the net fresh-water forcing seen \n"//&
                 "by the ocean (including restoring) to zero.", default=.false.)
  call get_param(param_file, mod, "ICE_SALT_CONCENTRATION", &
                 CS%ice_salt_concentration, &
                 "The assumed sea-ice salinity needed to reverse engineer the \n"//&
                 "melt flux (or ice-ocean fresh-water flux).", &
                 units="kg/kg", default=0.005)
  call get_param(param_file, mod, "USE_LIMITED_PATM_SSH", CS%use_limited_P_SSH, &
                 "If true, return the the sea surface height with the \n"//&
                 "correction for the atmospheric (and sea-ice) pressure \n"//&
                 "limited by max_p_surf instead of the full atmospheric \n"//&
                 "pressure.", default=.true.)

  call get_param(param_file, mod, "WIND_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the \n"//&
                 "staggering of the input wind stress field.  Valid \n"//&
                 "values are 'A', 'B', or 'C'.", default="C")
  if ((stagger(1:1) == 'a') .or. (stagger(1:1) == 'A')) then
    CS%wind_stagger = A_GRID
    call MOM_error(FATAL,"surface_forcing_init: A-grid input wind stagger "// &
                          "is not supported yet.")
  elseif ((stagger(1:1) == 'b') .or. (stagger(1:1) == 'B')) then
    CS%wind_stagger = B_GRID
  elseif ((stagger(1:1) == 'c') .or. (stagger(1:1) == 'C')) then
    CS%wind_stagger = C_GRID
  else
    call MOM_error(FATAL,"surface_forcing_init: #define WIND_STAGGER "// &
                    trim(stagger)//" is invalid.")
  endif

  if (restore_salt) then
    call get_param(param_file, mod, "FLUXCONST", CS%Flux_const, &
                 "The constant that relates the restoring surface fluxes \n"//&
                 "to the relative surface anomalies (akin to a piston \n"//&
                 "velocity).  Note the non-MKS units.", units="m day-1", &
                 fail_if_missing=.true.)
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
                 " If true, use an ice mask defined by frazil criteria to \n"//&
                 "determine where to apply salinity restoring.", default=.false.)
    call get_param(param_file, mod, "MASK_SRESTORE_MARGINAL_SEAS", &
                 CS%mask_srestore_marginal_seas, &
                 "If true, mask sss restoring in marginal seas.", default=.false.)
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

  call safe_alloc_ptr(CS%TKE_tidal,isd,ied,jsd,jed)   ; CS%TKE_tidal(:,:) = 0.0
  call safe_alloc_ptr(CS%ustar_tidal,isd,ied,jsd,jed) ; CS%ustar_tidal(:,:) = 0.0
  
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

    call safe_alloc_ptr(CS%gust,isd,ied,jsd,jed) ; CS%gust(:,:) = 0.0
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

  CS%id_taux = register_diag_field('ocean_model', 'taux', G%axesCu1, Time, &
        'Zonal Wind Stress', 'Pascal', standard_name='surface_downward_x_stress')
  CS%id_tauy = register_diag_field('ocean_model', 'tauy', G%axesCv1, Time, &
        'Meridional Wind Stress', 'Pascal', standard_name='surface_downward_y_stress')
  CS%id_ustar = register_diag_field('ocean_model', 'ustar', G%axesT1, Time, &
      'Surface friction velocity', 'meter second-1')

  if (CS%use_temperature) then
    CS%id_PminusE = register_diag_field('ocean_model', 'PmE', G%axesT1, Time, &
          'Net fresh water flux (P-E+C+R)', 'kilogram meter-2 second-1')
    CS%id_evap = register_diag_field('ocean_model', 'evap', G%axesT1, Time, &
          'Evaporation at ocean surface (usually negative)', 'kilogram meter-2 second-1')
    CS%id_precip = register_diag_field('ocean_model', 'precip', G%axesT1, Time, &
          'Precipitation into ocean', 'kilogram meter-2 second-1')
    CS%id_froz_precip = register_diag_field('ocean_model', 'froz_precip', G%axesT1, Time, &
          'Frozen Precipitation into ocean', 'kilogram meter-2 second-1', &
          standard_name='snowfall_flux')
    CS%id_liq_precip = register_diag_field('ocean_model', 'liq_precip', G%axesT1, Time, &
          'Liquid Precipitation into ocean', 'kilogram meter-2 second-1', &
          standard_name='rainfall_flux')
    CS%id_virt_precip = register_diag_field('ocean_model', 'virt_precip', G%axesT1, Time, &
          'Virtual Precipitation into ocean (due to salinity restoring)', 'kilogram meter-2 second-1')
    CS%id_froz_runoff = register_diag_field('ocean_model', 'froz_runoff', G%axesT1, Time, &
          'Frozen runoff (calving) into ocean', 'kilogram meter-2 second-1', &
          standard_name='water_flux_into_sea_water_from_icebergs')
    CS%id_liq_runoff = register_diag_field('ocean_model', 'liq_runoff', G%axesT1, Time, &
          'Liquid runoff (rivers) into ocean', 'kilogram meter-2 second-1', &
          standard_name='water_flux_into_sea_water_from_rivers')
    CS%id_calving_hflx = register_diag_field('ocean_model', 'calving_hflx', G%axesT1, Time, &
          'Heat content of frozen runoff (calving) into ocean', 'Watt meter-2')
    CS%id_runoff_hflx = register_diag_field('ocean_model', 'runoff_hflx', G%axesT1, Time, &
          'Heat content of liquid runoff (rivers) into ocean', 'Watt meter-2')

    CS%id_Net_Heating = register_diag_field('ocean_model', 'Net_Heat', G%axesT1, Time, &
          'Net Surface Heating of Ocean', 'Watt meter-2')
    CS%id_sw = register_diag_field('ocean_model', 'SW', G%axesT1, Time, &
        'Shortwave radiation flux into ocean', 'Watt meter-2', &
        standard_name='surface_net_downward_shortwave_flux')
    CS%id_LwLatSens = register_diag_field('ocean_model', 'LwLatSens', G%axesT1, Time, &
          'Combined longwave, latent, and sensible heating', 'Watt meter-2')
    CS%id_lw = register_diag_field('ocean_model', 'LW', G%axesT1, Time, &
        'Longwave radiation flux into ocean', 'Watt meter-2', &
        standard_name='surface_net_downward_longwave_flux')
    CS%id_lat = register_diag_field('ocean_model', 'latent', G%axesT1, Time, &
        'Latent heat flux into ocean due to fusion and evaporation', 'Watt meter-2')
    CS%id_sens = register_diag_field('ocean_model', 'sensible', G%axesT1, Time, &
        'Sensible heat flux into ocean', 'Watt meter-2', &
        standard_name='surface_downward_sensible_heat_flux')

    CS%id_psurf = register_diag_field('ocean_model', 'p_surf', G%axesT1, Time, &
          'Pressure at ice-ocean or atmosphere-ocean interface', 'Pascal')
    CS%id_saltflux = register_diag_field('ocean_model', 'salt_flux', G%axesT1, Time, &
          'Salt flux into ocean at surface', 'kilogram meter-2 second-1')
    CS%id_TKE_tidal = register_diag_field('ocean_model', 'TKE_tidal', G%axesT1, Time, &
          'Tidal source of BBL mixing', 'Watt meter-2')
  else
    CS%id_buoy = register_diag_field('ocean_model', 'buoy', G%axesT1, Time, &
          'Buoyancy forcing', 'meter2 second-3')
  endif

  if (present(restore_salt)) then ; if (restore_salt) then
    salt_file = trim(CS%inputdir) // "salt_restore.nc"
    CS%id_srestore = init_external_field(salt_file,'salt',domain=G%Domain%mpp_domain)
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

  if (present(fluxes)) then
    if (associated(fluxes%taux))        deallocate(fluxes%taux)
    if (associated(fluxes%tauy))        deallocate(fluxes%tauy)
    if (associated(fluxes%ustar))       deallocate(fluxes%ustar)
    if (associated(fluxes%buoy))        deallocate(fluxes%buoy)
    if (associated(fluxes%sw))          deallocate(fluxes%sw)
    if (associated(fluxes%sw_vis_dir))  deallocate(fluxes%sw_vis_dir)
    if (associated(fluxes%sw_vis_dif))  deallocate(fluxes%sw_vis_dif)
    if (associated(fluxes%sw_nir_dir))  deallocate(fluxes%sw_nir_dir)
    if (associated(fluxes%sw_nir_dif))  deallocate(fluxes%sw_nir_dif)
    if (associated(fluxes%lw))          deallocate(fluxes%lw)
    if (associated(fluxes%latent))      deallocate(fluxes%latent)
    if (associated(fluxes%sens))        deallocate(fluxes%sens)
    if (associated(fluxes%evap))        deallocate(fluxes%evap)
    if (associated(fluxes%liq_precip))  deallocate(fluxes%liq_precip)
    if (associated(fluxes%froz_precip)) deallocate(fluxes%froz_precip)
    if (associated(fluxes%virt_precip)) deallocate(fluxes%virt_precip)
    if (associated(fluxes%runoff_hflx)) deallocate(fluxes%runoff_hflx)
    if (associated(fluxes%calving_hflx)) deallocate(fluxes%calving_hflx)
    if (associated(fluxes%p_surf))      deallocate(fluxes%p_surf)
    if (associated(fluxes%salt_flux))   deallocate(fluxes%salt_flux)
    if (associated(fluxes%TKE_tidal))   deallocate(fluxes%TKE_tidal)
    if (associated(fluxes%ustar_tidal)) deallocate(fluxes%ustar_tidal)
    ! Deallocate any elements of fluxes%tr_fluxes.
    if (associated(fluxes%tr_fluxes))   deallocate(fluxes%tr_fluxes)
  endif

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
