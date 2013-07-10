module MOM_diabatic_driver
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
!*  By Robert Hallberg, April 1994 - July 2000                         *
!*                                                                     *
!*    This program contains the subroutine that, along with the        *
!*  subroutines that it calls, implements diapycnal mass and momentum  *
!*  fluxes and a bulk mixed layer.  The diapycnal diffusion can be     *
!*  used without the bulk mixed layer.                                 *
!*                                                                     *
!*    diabatic first determines the (diffusive) diapycnal mass fluxes  *
!*  based on the convergence of the buoyancy fluxes within each layer. *
!*  The dual-stream entrainment scheme of MacDougall and Dewar (JPO,   *
!*  1997) is used for combined diapycnal advection and diffusion,      *
!*  calculated implicitly and potentially with the Richardson number   *
!*  dependent mixing, as described by Hallberg (MWR, 2000). Diapycnal  *
!*  advection is fundamentally the residual of diapycnal diffusion,    *
!*  so the fully implicit upwind differencing scheme that is used is   *
!*  entirely appropriate.  The downward buoyancy flux in each layer    *
!*  is determined from an implicit calculation based on the previously *
!*  calculated flux of the layer above and an estimated flux in the    *
!*  layer below.  This flux is subject to the following conditions:    *
!*  (1) the flux in the top and bottom layers are set by the boundary  *
!*  conditions, and (2) no layer may be driven below an Angstrom thick-*
!*  ness.  If there is a bulk mixed layer, the buffer layer is treat-  *
!*  ed as a fixed density layer with vanishingly small diffusivity.    *
!*                                                                     *
!*    diabatic takes 5 arguments:  the two velocities (u and v), the   *
!*  thicknesses (h), a structure containing the forcing fields, and    *
!*  the length of time over which to act (dt).  The velocities and     *
!*  thickness are taken as inputs and modified within the subroutine.  *
!*  There is no limit on the time step.                                *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, T, S, buoy, ustar, ea, eb, etc.       *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_entrain_diffusive, only : entrainment_diffusive, entrain_diffusive_init
use MOM_entrain_diffusive, only : entrain_diffusive_end, entrain_diffusive_CS
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_diag_to_Z, only : diag_to_Z_CS, register_Zint_diag, calc_Zint_diags
use MOM_domains, only : pass_var, To_West, To_South
use MOM_checksums, only : hchksum, uchksum, vchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, optics_type, MOM_forcing_chksum
use MOM_forcing_type, only : extractFluxes1d, absorbRemainingSW
use MOM_geothermal, only : geothermal, geothermal_init, geothermal_end, geothermal_CS
use MOM_grid, only : ocean_grid_type
use MOM_io, only : vardesc
use MOM_int_tide_input, only : set_int_tide_input, int_tide_input_init
use MOM_int_tide_input, only : int_tide_input_end, int_tide_input_CS, int_tide_input_type
use MOM_internal_tides, only : propagate_int_tide, register_int_tide_restarts
use MOM_internal_tides, only : internal_tides_init, internal_tides_end, int_tide_CS
use MOM_kappa_shear, only : Calculate_kappa_shear, kappa_shear_init, Kappa_shear_CS
use MOM_bulk_mixed_layer, only : bulkmixedlayer, bulkmixedlayer_init, bulkmixedlayer_CS
use MOM_opacity, only : opacity_init, set_opacity, opacity_end, opacity_CS
use MOM_set_diffusivity, only : set_diffusivity, set_BBL_diffusivity
use MOM_set_diffusivity, only : set_diffusivity_init, set_diffusivity_end
use MOM_set_diffusivity, only : set_diffusivity_CS
use MOM_sponge, only : apply_sponge, sponge_CS
use MOM_tracer_flow_control, only : call_tracer_column_fns, tracer_flow_control_CS
use MOM_variables, only : thermo_var_ptrs, vertvisc_type, accel_diag_ptrs
use MOM_variables, only : cont_diag_ptrs, MOM_thermovar_chksum, p3d
use MOM_regularize_layers, only : regularize_layers, regularize_layers_init, regularize_layers_CS
use MOM_wave_speed, only : wave_speed
use MOM_EOS, only : calculate_density, calculate_2_densities, calculate_TFreeze

implicit none ; private

#include <MOM_memory.h>

public diabatic, diabatic_driver_init, diabatic_driver_end
public adiabatic, adiabatic_driver_init

type, public :: diabatic_CS ; private
  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used with
                             ! nkml sublayers (and additional buffer layers).
  logical :: use_kappa_shear ! If true, use the kappa_shear module to find the
                             ! shear-driven diapycnal diffusivity.
  logical :: use_sponge      ! If true, sponges may be applied anywhere in the
                             ! domain.  The exact location and properties of
                             ! those sponges are set by calls to
                             ! initialize_sponge and set_up_sponge_field.
  logical :: use_geothermal  ! If true, apply geothermal heating.
  logical :: use_int_tides
  logical :: useALEalgorithm ! If true, use the ALE algorithm rather than layered
                             ! isopycnal/stacked shallow water mode. This logical
                             ! passed by argument to diabatic_driver_init.
  real :: ML_mix_first       !   The nondimensional fraction of the mixed layer
                             ! algorithm that is applied before diffusive mixing.
                             ! The default is 0, while 0.5 gives Strang splitting
                             ! and 1 is a sensible value too.  Note that if there
                             ! are convective instabilities in the initial state,
                             ! the first call may do much more than the second.
  logical :: massless_match_targets  ! If true (the default), keep the T & S
                             ! consistent with the target values.
  logical :: mix_boundary_tracers ! If true, mix the passive tracers in massless
                             ! layers at the bottom into the interior as though
                             ! a diffusivity of Kd_min_tr (see below) were
                             ! operating.
  real    :: Kd_BBL_tr       !   A bottom boundary layer tracer diffusivity that
                             ! will allow for explicitly specified bottom fluxes
                             ! in m2 s-1.  The entrainment at the bottom is at
                             ! least sqrt(Kd_BBL_tr*dt) over the same distance.
  real    :: Kd_min_tr       !   A minimal diffusivity that should always be
                             ! applied to tracers, especially in massless layers
                             ! near the bottom, in m2 s-1.
  logical :: reclaim_frazil  !   If true, try to use any frazil heat deficit to
                             ! to cool the topmost layer down to the freezing
                             ! point.  The default is false.                           
  logical :: pressure_dependent_frazil  ! If true, use a pressure dependent
                             ! freezing temperature when making frazil.  The
                             ! default is false, which will be faster but is
                             ! inappropriate with ice-shelf cavities.
  logical :: debug           ! If true, write verbose checksums for debugging purposes.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer :: id_dudt_dia = -1, id_dvdt_dia = -1, id_wd = -1
  integer :: id_ea = -1 , id_eb = -1, id_Kd_z = -1
  integer :: id_Tdif_z = -1, id_Tadv_z = -1, id_Sdif_z = -1, id_Sadv_z = -1
  integer :: id_Tdif = -1, id_Tadv = -1, id_Sdif = -1, id_Sadv = -1

  type(entrain_diffusive_CS), pointer :: entrain_diffusive_CSp => NULL()
  type(bulkmixedlayer_CS),    pointer :: bulkmixedlayer_CSp => NULL()
  type(regularize_layers_CS), pointer :: regularize_layers_CSp => NULL()
  type(geothermal_CS),        pointer :: geothermal_CSp => NULL()
  type(Kappa_shear_CS),       pointer :: kappa_shear_CSp => NULL()
  type(int_tide_CS),          pointer :: int_tide_CSp => NULL()
  type(int_tide_input_CS),    pointer :: int_tide_input_CSp => NULL()
  type(int_tide_input_type),  pointer :: int_tide_input => NULL()
  type(opacity_CS),           pointer :: opacity_CSp => NULL()
  type(set_diffusivity_CS),   pointer :: set_diff_CSp => NULL()
  type(sponge_CS),            pointer :: sponge_CSp => NULL()
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp => NULL()
  type(optics_type),          pointer :: optics => NULL()
  type(diag_to_Z_CS),         pointer :: diag_to_Z_CSp => NULL()
end type diabatic_CS

integer :: id_clock_entrain, id_clock_mixedlayer, id_clock_set_diffusivity
integer :: id_clock_uv_at_h, id_clock_frazil, id_clock_kappa_shear
integer :: id_clock_tracers, id_clock_tridiag, id_clock_pass, id_clock_sponge
integer :: id_clock_geothermal, id_clock_double_diff, id_clock_remap

contains

subroutine diabatic(u, v, h, tv, fluxes, visc, ADp, CDp, dt, G, CS)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(inout) :: h
  type(thermo_var_ptrs),                  intent(inout) :: tv
  type(forcing),                          intent(inout) :: fluxes
  type(vertvisc_type),                    intent(inout) :: visc
  type(accel_diag_ptrs),                  intent(inout) :: ADp
  type(cont_diag_ptrs),                   intent(inout) :: CDp
  real,                                   intent(in)    :: dt
  type(ocean_grid_type),                  intent(inout) :: G
  type(diabatic_CS),                      pointer       :: CS
!    This subroutine imposes the diapycnal mass fluxes and the
!  accompanying diapycnal advection of momentum and tracers.

! Arguments: u - Zonal velocity, in m s-1. (Intent in/out)
!  (in/out)  v - Meridional velocity, in m s-1.
!  (in/out)  h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in/out)  visc - A structure containing vertical viscosities, bottom boundary
!                   layer properies, and related fields.
!  (inout)   ADp - A structure with pointers to the various accelerations in
!                  the momentum equations, to enable the later calculation
!                  of derived diagnostics, like energy budgets.
!  (inout)   CDp - A structure with pointers to various terms in the continuity
!                  equations.
!  (in)      dt - Time increment, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 diabatic_driver_init.

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    ea, &    ! ea is the amount of fluid entrained from the layer above within
             ! one time step, in units of m or kg m-2.
    eb, &    ! eb is the amount of fluid entrained from the layer below within
             ! one time step, in units of m or kg m-2.
    Kd, &    ! The diapycnal diffusivity of layers, in m2 s-1.
    h_orig, &  ! The initial layer thicknesses, in m or kg m-2.
    hold, &  ! The layer thickness before diapycnal entrainment, and later
             ! the initial layer thicknesses (if a mixed layer is used),
             ! in units of m or kg m-2.
    u_h, &   ! Zonal and meridional velocities at thickness points after
    v_h      ! entrainment, in m s-1.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    cg1, &   ! The first baroclinic gravity wave speed.
    Rcv_ml   ! The coordinate density of the mixed layer, used only for applying
             ! sponges.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target :: &
    ! These are targets so that the space can be shared with eaml & ebml.
    eatr, &  ! The equivalent of ea and eb for tracers, which differ from ea and
    ebtr     ! eb in that they tend to homogenize tracers in massless layers
             ! near the boundaries, in units of m or kg m-2.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), target :: &
    Kd_int, &  ! The diapycnal diffusivity of interfaces, in m2 s-1.
    Tdif_flx, &! The diffusive diapycnal heat flux across interfaces, in K m s-1.
    Tadv_flx, &! The advective diapycnal heat flux across interfaces, in K m s-1.
    Sdif_flx, &! The diffusive diapycnal salt flux across interfaces, in PSU m s-1.
    Sadv_flx   ! The advective diapycnal salt flux across interfaces, in PSU m s-1.

  ! The following 5 variables are only used with a bulk mixed layer.
  real, pointer, dimension(:,:,:) :: &
    eaml, &  ! The equivalent of ea and eb due to mixed layer processes,
    ebml     ! in m or kg m-2.  These will be pointers to eatr and ebtr so
             ! as to reuse the memory as the arrays are not needed at the
             ! same time.

  integer :: kb(SZI_(G),SZJ_(G)) ! The index of the lightest layer denser
                              ! than the buffer layer.  Nondimensional.
  real :: p_ref_cv(SZI_(G))   ! Reference pressure for the potential
                              ! density which defines the coordinate
                              ! variable, set to P_Ref, in Pa.

  logical :: in_boundary(SZI_(G)) ! True if there are no massive layers below,
                       ! where massive is defined as sufficiently thick that
                       ! the no-flux boundary conditions have not restricted
                       ! the entrainment - usually sqrt(Kd*dt).
  real :: b_denom_1    ! The first term in the denominator of b1 in m or kg m-2.
  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected, in m or kg m-2.
  real :: h_neglect2   ! h_neglect^2, in m2 or kg2 m-4.
  real :: add_ent      ! Entrainment that needs to be added when mixing tracers,
                       ! in units of m or kg m-2.
  real :: eaval        ! eaval is 2*ea at velocity grid points, in m or kg m-2.
  real :: hval         ! hval is 2*h at velocity grid points, in m or kg m-2.
  real :: h_tr         ! h_tr is h at tracer points with a tiny thickness
                       ! added to ensure positive definiteness, in m or kg m-2.
  real :: Tr_ea_BBL    ! The diffusive tracer thickness in the BBL that is
                       ! coupled to the bottom within a timestep, in m.
  real :: htot(SZIB_(G))  ! The summed thickness from the bottom, in m.
  real :: b1(SZIB_(G)), d1(SZIB_(G)) ! b1, c1, and d1 are variables used by the
  real :: c1(SZIB_(G),SZK_(G))       ! tridiagonal solver.
  real :: dt_mix       ! The amount of time over which to apply mixing, in s.
  real :: Idt          ! The inverse of the time step, in s-1.
  type(p3d) :: z_ptrs(7) ! Pointers to the diagnostics that are to be
                       ! interpolated into depth space.
  integer :: num_z_diags  ! The number of diagnostics that are to be
                       ! interpolated into depth space.
  integer :: z_ids(7)  ! The id numbers of the diagnostics that are to be
                       ! interpolated into depth space.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  real, pointer :: T(:,:,:), S(:,:,:)
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nkmb = G%nk_rho_varies
  h_neglect = G%H_subroundoff ; h_neglect2 = h_neglect*h_neglect

  if (nz == 1) return

  T => tv%T
  S => tv%S

  ! This sets the equivalence between the same bits of memory for these arrays.
  eaml => eatr ; ebml => ebtr

  Idt = 1.0 / dt

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_diabatic_driver: "// &
         "Module must be initialized before it is used.")

  if (CS%debug) then
    call MOM_state_chksum("Start of diabatic ", u(:,:,:), v(:,:,:), h(:,:,:), G)
  endif

  call cpu_clock_begin(id_clock_set_diffusivity)
  call set_BBL_diffusivity(u, v, h, fluxes, visc, G, CS%set_diff_CSp)
  call cpu_clock_end(id_clock_set_diffusivity)

  !   Frazil formation keeps the temperature above the freezing point.
  ! make_frazil is deliberately called at both the beginning and at
  ! the end of the diabatic processes.
  if (ASSOCIATED(T) .AND. ASSOCIATED(tv%frazil)) call make_frazil(h,tv,G,CS)

  if ((CS%ML_mix_first > 0.0) .or. CS%use_geothermal) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h_orig(i,j,k) = h(i,j,k) ; eaml(i,j,k) = 0.0 ; ebml(i,j,k) = 0.0
    enddo ; enddo ; enddo
  endif
  if (CS%use_geothermal) then
    call cpu_clock_begin(id_clock_geothermal)
    call geothermal(h, tv, dt, eaml, ebml, G, CS%geothermal_CSp)
    call cpu_clock_end(id_clock_geothermal)
  endif

  ! Set_opacity estimates the optical properties of the water column.
  !   It will need to be modified later to include information about the
  ! biological properties and layer thicknesses.
  if (associated(CS%optics)) &
    call set_opacity(CS%optics, fluxes, G, CS%opacity_CSp)
  
  if (CS%bulkmixedlayer) then
    if (CS%ML_mix_first > 0.0) then
!  This subroutine (1)  Cools the mixed layer.
!    (2) Performs convective adjustment by mixed layer entrainment.
!    (3) Heats the mixed layer and causes it to detrain to
!        Monin-Obukhov depth or minimum mixed layer depth.
!    (4) Uses any remaining TKE to drive mixed layer entrainment.
!    (5) Possibly splits the buffer layer into two isopycnal layers.
      call find_uv_at_h(u, v, h, u_h, v_h, G)

      call cpu_clock_begin(id_clock_mixedlayer)
      if (CS%ML_mix_first < 1.0) then
        call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt*CS%ML_mix_first, &
                            eaml,ebml, G, CS%bulkmixedlayer_CSp, CS%optics, dt, &
                            last_call=.false.)
      else
        call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt, eaml, ebml, &
                        G, CS%bulkmixedlayer_CSp, CS%optics, dt, last_call=.true.)
      endif
!  Keep salinity from falling below a small but positive threshold.
!  This occurs when the ice model attempts to extract more salt than
!  is actually present in the ocean.  
      if (ASSOCIATED(S) .and. ASSOCIATED(tv%salt_deficit)) &
        call adjust_salt(h, tv, G, CS)
      call cpu_clock_end(id_clock_mixedlayer)
      if (CS%debug) call MOM_state_chksum("After mixedlayer ", u, v, h, G)
    endif

  endif

!    This following subroutine applies diapycnal diffusion and the
!  surface buoyancy forcing.  No layer is less than an Angstrom thick at
!  the end, and total density (heat or salt) is conserved unless all
!  of the mass is in the lightest or densest layer.
!  When the bulk mixed layer is used, the flux from the buffer layer
!  is applied to the next denser layer, while all lighter layers are
!  effectively removed from the calculation.  Also, buoyancy forcing
!  is applied to the bulk mixed layer elsewhere.

!   Calculate appropriately limited diapycnal mass fluxes to account
! for diapycnal diffusion and advection.
  if (CS%debug) then
    call MOM_state_chksum("before use_kappa_shear ", u(:,:,:), v(:,:,:), h(:,:,:), G)
  endif
  if (CS%use_kappa_shear) then
    if ((CS%ML_mix_first > 0.0) .or. CS%use_geothermal) then
      call find_uv_at_h(u, v, h_orig, u_h, v_h, G, eaml, ebml)
      if (CS%debug) then
        call hchksum(eaml, "aft find_uv_at_h eaml",G)
        call hchksum(ebml, "aft find_uv_at_h ebml",G)
      endif
    else
      call find_uv_at_h(u, v, h, u_h, v_h, G)
    endif
    if (CS%debug) then
      call MOM_state_chksum("aft find_uv_at_h KS ", u(:,:,:), v(:,:,:), h(:,:,:), G)
      call MOM_forcing_chksum("before calc_KS ", fluxes, G, haloshift=0)
      call MOM_thermovar_chksum("before calc_KS ", tv, G)
      call uchksum(u_h, "before calc_KS u_h",G)
      call vchksum(v_h, "before calc_KS v_h",G)
    endif

    call cpu_clock_begin(id_clock_kappa_shear)
    call Calculate_kappa_shear(u_h, v_h, h, tv, fluxes%p_surf, visc%Kd_turb, visc%TKE_turb, &
                               dt, G, CS%kappa_shear_CSp)
    call cpu_clock_end(id_clock_kappa_shear)
    if (CS%debug) then
      call MOM_state_chksum("after Calc_KS ", u(:,:,:), v(:,:,:), h(:,:,:), G)
      call MOM_forcing_chksum("after calc_KS ", fluxes, G, haloshift=0)
      call MOM_thermovar_chksum("after calc_KS ", tv, G)
      call hchksum(visc%Kd_turb, "after calc_KS visc%Kd_turb",G)
      call hchksum(visc%TKE_turb, "after calc_KS visc%TKE_turb",G)
    endif
  endif

  if (CS%use_int_tides) then
    call set_int_tide_input(u, v, h, tv, fluxes, CS%int_tide_input, dt, G, &
                            CS%int_tide_input_CSp)
    cg1(:,:) = 0.0
    call wave_speed(h, tv, G, cg1, full_halos=.true.)
    call propagate_int_tide(cg1, CS%int_tide_input%TKE_itidal_input, &
                            CS%int_tide_input%tideamp, dt, G, CS%int_tide_CSp)
  endif

  call cpu_clock_begin(id_clock_set_diffusivity)
  call set_diffusivity(u, v, h, tv, fluxes, visc, dt, G, CS%set_diff_CSp, Kd, Kd_int)
  call cpu_clock_end(id_clock_set_diffusivity)
  if (CS%debug) then
    call MOM_state_chksum("after set_diffusivity ", u(:,:,:), v(:,:,:), h(:,:,:), G)
    call MOM_forcing_chksum("after set_diffusivity ", fluxes, G, haloshift=0)
    call MOM_thermovar_chksum("after set_diffusivity ", tv, G)
    call hchksum(Kd, "after set_diffusivity Kd",G,haloshift=0)
    call hchksum(Kd_Int, "after set_diffusivity Kd_Int",G,haloshift=0)
  endif

  if (associated(visc%Kd_extra_T) .and. associated(visc%Kd_extra_S) .and. &
      associated(T)) then
    call cpu_clock_begin(id_clock_double_diff)
    call double_diffuse_T_S(h, tv, visc, dt, G)
    call cpu_clock_end(id_clock_double_diff)
  endif

  ! If using the ALE algorithm, set ea=eb=Kd on interfaces for
  ! use in the tri-diagonal solver.
  ! Otherwise, call entrainment_diffusive() which sets ea and eb
  ! based on KD and target densities (ie. does remapping as well).
  if (CS%useALEalgorithm) then
    do j=js,je ; do i=is,ie
      ea(i,j,1) = 0.
    enddo ; enddo
    do k=2,nz ; do j=js,je ; do i=is,ie
      hval=1.0/(h_neglect + 0.5*(h(i,j,k-1) + h(i,j,k)))
      ea(i,j,k) = (G%m_to_H**2) * dt * hval * Kd_int(i,j,k) 
    ! Alternative to use Kd rather than Kd_int:
    ! ea(i,j,k) = (G%m_to_H**2) * dt * hval * &
    !      0.5*(h(i,j,k)*Kd(i,j,k-1)+h(i,j,k-1)*Kd(i,j,k))*hval
      eb(i,j,k-1) = ea(i,j,k)
    enddo ; enddo ; enddo
    do j=js,je ; do i=is,ie
      eb(i,j,nz) = 0.
    enddo ; enddo
!   This block could replace the more general block below ... - AJA ?
!   do k=1,nz ; do j=js,je ; do i=is,ie
!     hold(i,j,k) = h(i,j,k)
!     if (h(i,j,k) <= 0.0) h(i,j,k) = G%Angstrom
!   enddo ; enddo ; enddo
  else ! .not. CS%useALEalgorithm
    ! If not useing ALE, then calculate layer entrainments/detrainments from
    ! diffusivities and differences between layer and target densities
    call cpu_clock_begin(id_clock_entrain)
    call Entrainment_diffusive(u, v, h, tv, fluxes, dt, G, CS%entrain_diffusive_CSp, &
                               ea, eb, kb, Kd_Lay=Kd, Kd_int=Kd_int)
    call cpu_clock_end(id_clock_entrain)
  endif ! (CS%useALEalgorithm)

  if (CS%debug) then
    call MOM_forcing_chksum("after calc_entrain ", fluxes, G, haloshift=0)
    call MOM_thermovar_chksum("after calc_entrain ", tv, G)
    call MOM_state_chksum("after calc_entrain ", u(:,:,:), v(:,:,:), h(:,:,:), G)
    call hchksum(G%H_to_m*ea, "after calc_entrain ea",G,haloshift=0)
    call hchksum(G%H_to_m*eb, "after calc_entrain eb",G,haloshift=0)
  endif

  ! Update h according to divergence of the difference between
  ! ea and eb. We keep a record of the original h in hold.
  ! In the following, the checks for negative values are to guard
  ! against against unforeseen instances where the entrainment
  ! drives a layer to negative thickness.  This will never happen if
  ! enough iterations are permitted in Calculate_Entrainment, and
  ! even if too few iterations are allowed, it is still guarded
  ! against.  In other words the checks are probably unnecessary.
  do j=js,je
    do i=is,ie
      hold(i,j,1) = h(i,j,1)
      h(i,j,1) = h(i,j,1) + (eb(i,j,1) - ea(i,j,2))
      hold(i,j,nz) = h(i,j,nz)
      h(i,j,nz) = h(i,j,nz) + (ea(i,j,nz) - eb(i,j,nz-1))
      if (h(i,j,1) <= 0.0) h(i,j,1) = G%Angstrom
      if (h(i,j,nz) <= 0.0) h(i,j,nz) = G%Angstrom
    enddo
    do k=2,nz-1 ; do i=is,ie
      hold(i,j,k) = h(i,j,k)
      h(i,j,k) = h(i,j,k) + ((ea(i,j,k) - eb(i,j,k-1)) + &
                    (eb(i,j,k) - ea(i,j,k+1)))
      if (h(i,j,k) <= 0.0) h(i,j,k) = G%Angstrom
    enddo ; enddo
  enddo
  if (CS%debug) then
    call MOM_state_chksum("after negative check ", u(:,:,:), v(:,:,:), h(:,:,:), G)
  endif

  ! Apply forcing when using the ALE algorithm
  if (CS%useALEalgorithm) then
    call cpu_clock_begin(id_clock_remap)
    call applyBoundaryFluxes(G, dt, fluxes, CS%optics, ea, eb, h, tv)
    call cpu_clock_end(id_clock_remap)
  endif
  
  ! Here, T and S are updated according to ea and eb.
  ! If using the bulk mixed layer, T and S are also updated
  ! by surface fluxes (in fluxes%*).
  if (CS%bulkmixedlayer) then
    if (ASSOCIATED(T)) then
      call cpu_clock_begin(id_clock_tridiag)
      ! Temperature and salinity (as state variables) are treated slightly
      ! differently from other tracers to insure that massless layers that
      ! are lighter than the mixed layer have temperatures and salinities
      ! that correspond to their prescribed densities.
      do j=js,je

        if (CS%massless_match_targets) then
          do i=is,ie
            h_tr = hold(i,j,1) + h_neglect
            b1(i) = 1.0 / (h_tr + eb(i,j,1))
            d1(i) = h_tr * b1(i)
            T(i,j,1) = b1(i) * (h_tr*T(i,j,1))
            S(i,j,1) = b1(i) * (h_tr*S(i,j,1))
          enddo
          do k=2,nkmb ; do i=is,ie
            c1(i,k) = eb(i,j,k-1) * b1(i)
            h_tr = hold(i,j,k) + h_neglect
            b_denom_1 = h_tr + d1(i)*ea(i,j,k)
            b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
            if (k<nkmb) d1(i) = b_denom_1 * b1(i)
            T(i,j,k) = b1(i) * (h_tr*T(i,j,k) + ea(i,j,k)*T(i,j,k-1))
            S(i,j,k) = b1(i) * (h_tr*S(i,j,k) + ea(i,j,k)*S(i,j,k-1))
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
              T(i,j,k) = b1(i) * (h_tr*T(i,j,k) + ea(i,j,k)*T(i,j,nkmb))
              S(i,j,k) = b1(i) * (h_tr*S(i,j,k) + ea(i,j,k)*S(i,j,nkmb))
            elseif (k > kb(i,j)) then
              c1(i,k) = eb(i,j,k-1) * b1(i)
              h_tr = hold(i,j,k) + h_neglect
              b_denom_1 = h_tr + d1(i)*ea(i,j,k)
              b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
              d1(i) = b_denom_1 * b1(i)
              T(i,j,k) = b1(i) * (h_tr*T(i,j,k) + ea(i,j,k)*T(i,j,k-1))
              S(i,j,k) = b1(i) * (h_tr*S(i,j,k) + ea(i,j,k)*S(i,j,k-1))
            elseif (eb(i,j,k) < eb(i,j,k-1)) then ! (note that k < kb(i,j))
              !   The bottommost buffer layer might entrain all the mass from some
              ! of the interior layers that are thin and lighter in the coordinate
              ! density than that buffer layer.  The T and S of these newly
              ! massless interior layers are unchanged.
              T(i,j,nkmb) = T(i,j,nkmb) + b1(i) * (eb(i,j,k-1) - eb(i,j,k)) * T(i,j,k)
              S(i,j,nkmb) = S(i,j,nkmb) + b1(i) * (eb(i,j,k-1) - eb(i,j,k)) * S(i,j,k)
            endif
          enddo ; enddo

          do k=nz-1,nkmb,-1 ; do i=is,ie
            if (k >= kb(i,j)) then
              T(i,j,k) = T(i,j,k) + c1(i,k+1)*T(i,j,k+1)
              S(i,j,k) = S(i,j,k) + c1(i,k+1)*S(i,j,k+1)
            endif
          enddo ; enddo
          do i=is,ie ; if (kb(i,j) <= nz) then
            T(i,j,nkmb) = T(i,j,nkmb) + c1(i,kb(i,j))*T(i,j,kb(i,j))
            S(i,j,nkmb) = S(i,j,nkmb) + c1(i,kb(i,j))*S(i,j,kb(i,j))
          endif ; enddo
          do k=nkmb-1,1,-1 ; do i=is,ie
            T(i,j,k) = T(i,j,k) + c1(i,k+1)*T(i,j,k+1)
            S(i,j,k) = S(i,j,k) + c1(i,k+1)*S(i,j,k+1)
          enddo ; enddo
        else
          ! This simpler form allows T & S to be too dense for the layers
          ! between the buffer layers and the interior.
          do i=is,ie
            h_tr = hold(i,j,1) + h_neglect
            b1(i) = 1.0 / (h_tr + eb(i,j,1))
            d1(i) = h_tr * b1(i)
            T(i,j,1) = b1(i) * (h_tr*T(i,j,1))
            S(i,j,1) = b1(i) * (h_tr*S(i,j,1))
          enddo
          do k=2,nz ; do i=is,ie
            c1(i,k) = eb(i,j,k-1) * b1(i)
            h_tr = hold(i,j,k) + h_neglect
            b_denom_1 = h_tr + d1(i)*ea(i,j,k)
            b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
            d1(i) = b_denom_1 * b1(i)
            T(i,j,k) = b1(i) * (h_tr*T(i,j,k) + ea(i,j,k)*T(i,j,k-1))
            S(i,j,k) = b1(i) * (h_tr*S(i,j,k) + ea(i,j,k)*S(i,j,k-1))
          enddo ; enddo
          do k=nz-1,1,-1 ; do i=is,ie
            T(i,j,k) = T(i,j,k) + c1(i,k+1)*T(i,j,k+1)
            S(i,j,k) = S(i,j,k) + c1(i,k+1)*S(i,j,k+1)
          enddo ; enddo
        endif
      enddo ! end of j loop
      call cpu_clock_end(id_clock_tridiag)
    endif ! end of ASSOCIATED(T)

    if ((CS%ML_mix_first > 0.0) .or. CS%use_geothermal) then
      ! The mixed layer code has already been called, but there is some needed
      ! bookkeeping.
      do k=1,nz ; do j=js,je ; do i=is,ie
        hold(i,j,k) = h_orig(i,j,k)
        ea(i,j,k) = ea(i,j,k) + eaml(i,j,k)
        eb(i,j,k) = eb(i,j,k) + ebml(i,j,k)
      enddo ; enddo ; enddo
    endif

    if (CS%ML_mix_first < 1.0) then
  !  Call the mixed layer code now, perhaps for a second time.
  !  This subroutine (1)  Cools the mixed layer.
  !    (2) Performs convective adjustment by mixed layer entrainment.
  !    (3) Heats the mixed layer and causes it to detrain to
  !        Monin-Obukhov depth or minimum mixed layer depth.
  !    (4) Uses any remaining TKE to drive mixed layer entrainment.
  !    (5) Possibly splits the buffer layer into two isopycnal layers.

      call find_uv_at_h(u, v, hold, u_h, v_h, G, ea, eb)
      if (CS%debug) call MOM_state_chksum("find_uv_at_h1 ", u, v, h, G)

      dt_mix = min(dt,dt*(1.0 - CS%ML_mix_first))
      call cpu_clock_begin(id_clock_mixedlayer)
      call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt_mix, ea, eb, &
                      G, CS%bulkmixedlayer_CSp, CS%optics, dt, last_call=.true.)

!  Keep salinity from falling below a small but positive threshold.
!  This occurs when the ice model attempts to extract more salt than
!  is actually present in the ocean.  
      if (ASSOCIATED(S) .and. ASSOCIATED(tv%salt_deficit)) &
        call adjust_salt(h, tv, G, CS)

      call cpu_clock_end(id_clock_mixedlayer)
      if (CS%debug) call MOM_state_chksum("after Mixedlayer ", u, v, h, G)
    endif

  else                                             ! Not BULKMIXEDLAYER.

! Calculate the change in temperature & salinity due to entrainment.
    if (ASSOCIATED(T)) then
      call cpu_clock_begin(id_clock_tridiag)
      do j=js,je
        do i=is,ie
          h_tr = hold(i,j,1) + h_neglect
          b1(i) = 1.0 / (h_tr + eb(i,j,1))
          d1(i) = h_tr * b1(i)
          T(i,j,1) = (b1(i)*h_tr)*T(i,j,1)
          S(i,j,1) = (b1(i)*h_tr)*S(i,j,1)
        enddo
        do k=2,nz ; do i=is,ie
          c1(i,k) = eb(i,j,k-1) * b1(i)
          h_tr = hold(i,j,k) + h_neglect
          b_denom_1 = h_tr + d1(i)*ea(i,j,k)
          b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
          d1(i) = b_denom_1 * b1(i)
          T(i,j,k) = b1(i) * (h_tr*T(i,j,k) + ea(i,j,k)*T(i,j,k-1))
          S(i,j,k) = b1(i) * (h_tr*S(i,j,k) + ea(i,j,k)*S(i,j,k-1))
        enddo ; enddo
        do k=nz-1,1,-1 ; do i=is,ie
          T(i,j,k) = T(i,j,k) + c1(i,k+1)*T(i,j,k+1)
          S(i,j,k) = S(i,j,k) + c1(i,k+1)*S(i,j,k+1)
        enddo ; enddo
      enddo
      call cpu_clock_end(id_clock_tridiag)
    endif

  endif                                          ! end BULKMIXEDLAYER

  if (.not. CS%useALEalgorithm) then
    call cpu_clock_begin(id_clock_remap)
    call regularize_layers(h, tv, dt, ea, eb, G, CS%regularize_layers_CSp)
    call cpu_clock_end(id_clock_remap)
  endif

  if ((CS%id_Tdif > 0) .or. (CS%id_Tdif_z > 0) .or. &
      (CS%id_Tadv > 0) .or. (CS%id_Tadv_z > 0)) then
    do j=js,je ; do i=is,ie
      Tdif_flx(i,j,1) = 0.0 ; Tdif_flx(i,j,nz+1) = 0.0
      Tadv_flx(i,j,1) = 0.0 ; Tadv_flx(i,j,nz+1) = 0.0
    enddo ; enddo
    do K=2,nz ; do j=js,je ; do i=is,ie
      Tdif_flx(i,j,K) = (Idt * 0.5*(ea(i,j,k) + eb(i,j,k-1))) * &
                        (T(i,j,k-1) - T(i,j,k))
      Tadv_flx(i,j,K) = (Idt * (ea(i,j,k) - eb(i,j,k-1))) * &
                    0.5*(T(i,j,k-1) + T(i,j,k))
    enddo ; enddo ; enddo
  endif
  if ((CS%id_Sdif > 0) .or. (CS%id_Sdif_z > 0) .or. &
      (CS%id_Sadv > 0) .or. (CS%id_Sadv_z > 0)) then
    do j=js,je ; do i=is,ie
      Sdif_flx(i,j,1) = 0.0 ; Sdif_flx(i,j,nz+1) = 0.0
      Sadv_flx(i,j,1) = 0.0 ; Sadv_flx(i,j,nz+1) = 0.0
    enddo ; enddo
    do K=2,nz ; do j=js,je ; do i=is,ie
      Sdif_flx(i,j,K) = (Idt * 0.5*(ea(i,j,k) + eb(i,j,k-1))) * &
                        (S(i,j,k-1) - S(i,j,k))
      Sadv_flx(i,j,K) = (Idt * (ea(i,j,k) - eb(i,j,k-1))) * &
                    0.5*(S(i,j,k-1) + S(i,j,k))
    enddo ; enddo ; enddo
  endif

  call cpu_clock_begin(id_clock_tracers)
  if (CS%mix_boundary_tracers) then
    Tr_ea_BBL = sqrt(dt*CS%Kd_BBL_tr)

    do j=js,je
      do i=is,ie
        ebtr(i,j,nz) = eb(i,j,nz)
        htot(i) = 0.0
        in_boundary(i) = (G%mask2dT(i,j) > 0.0)
      enddo
      do k=nz,2,-1 ; do i=is,ie
        if (in_boundary(i)) then
          htot(i) = htot(i) + h(i,j,k)
          !   If diapyncal mixing has been suppressed because this is a massless
          ! layer near the bottom, add some mixing of tracers between these
          ! layers.  This flux is based on the harmonic mean of the two
          ! thicknesses, as this corresponds pretty closely (to within
          ! differences in the density jumps between layers) with what is done
          ! in the calculation of the fluxes in the first place.  Kd_min_tr
          ! should be much less than the values that have been set in Kd,
          ! perhaps a molecular diffusivity.
          add_ent = ((dt * CS%Kd_min_tr) * G%m_to_H**2) * &
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
          add_ent = ((dt * visc%Kd_extra_S(i,j,k)) * G%m_to_H**2) / &
             (0.25 * ((h(i,j,k-1) + h(i,j,k)) + (hold(i,j,k-1) + hold(i,j,k))) + &
              h_neglect)
          ebtr(i,j,k-1) = ebtr(i,j,k-1) + add_ent
          eatr(i,j,k) = eatr(i,j,k) + add_ent
        endif ; endif
      enddo ; enddo
      do i=is,ie ; eatr(i,j,1) = ea(i,j,1) ; enddo
    enddo
    call call_tracer_column_fns(hold, h, eatr, ebtr, fluxes, dt, G, tv, &
                                CS%optics, CS%tracer_flow_CSp)
  elseif (associated(visc%Kd_extra_S)) then
    do j=js,je ; do i=is,ie
      ebtr(i,j,nz) = eb(i,j,nz) ; eatr(i,j,1) = ea(i,j,1)
    enddo ; enddo
    do k=nz,2,-1 ; do j=js,je ; do i=is,ie
      if (visc%Kd_extra_S(i,j,k) > 0.0) then
        add_ent = ((dt * visc%Kd_extra_S(i,j,k)) * G%m_to_H**2) / &
           (0.25 * ((h(i,j,k-1) + h(i,j,k)) + (hold(i,j,k-1) + hold(i,j,k))) + &
            h_neglect)
      else
        add_ent = 0.0
      endif
      ebtr(i,j,k-1) = eb(i,j,k-1) + add_ent
      eatr(i,j,k) = ea(i,j,k) + add_ent
    enddo ; enddo ; enddo
    call call_tracer_column_fns(hold, h, eatr, ebtr, fluxes, dt, G, tv, &
                                CS%optics, CS%tracer_flow_CSp)
  else
    call call_tracer_column_fns(hold, h, ea, eb, fluxes, dt, G, tv, &
                                CS%optics, CS%tracer_flow_CSp)
  endif
  call cpu_clock_end(id_clock_tracers)

  if (CS%use_sponge) then
    call cpu_clock_begin(id_clock_sponge)
    if (CS%bulkmixedlayer .and. ASSOCIATED(tv%eqn_of_state)) then
      do i=is,ie ; p_ref_cv(i) = tv%P_Ref ; enddo
      do j=js,je
        call calculate_density(T(:,j,1), S(:,j,1), p_ref_cv, Rcv_ml(:,j), &
                               is, ie-is+1, tv%eqn_of_state)
      enddo
      call apply_sponge(h, dt, G, ea, eb, CS%sponge_CSp, Rcv_ml)
    else
      call apply_sponge(h, dt, G, ea, eb, CS%sponge_CSp)
    endif
    call cpu_clock_end(id_clock_sponge)
    if (CS%debug) then
      call MOM_state_chksum("apply_sponge ", u, v, h, G)
    endif
  endif

!   Save the diapycnal mass fluxes as a diagnostic field.
  if (ASSOCIATED(CDp%diapyc_vel)) then
    do K=2,nz ; do j=js,je ; do i=is,ie
      CDp%diapyc_vel(i,j,K) = Idt * (G%H_to_m * (ea(i,j,k) - eb(i,j,k-1)))
    enddo ; enddo ; enddo
    do j=js,je ; do i=is,ie
      CDp%diapyc_vel(i,j,1) = 0.0
      CDp%diapyc_vel(i,j,nz+1) = 0.0
    enddo ; enddo
  endif

!   For momentum, it is only the net flux that homogenizes within
!  the mixed layer.  Vertical viscosity that is proportional to the
!  mixed layer turbulence is applied elsewhere.
  if (CS%bulkmixedlayer) then
    do k=2,G%nkml ; do j=js,je ; do i=is,ie
      if (ea(i,j,k) >= eb(i,j,k-1)) then
        ea(i,j,k) = ea(i,j,k) - eb(i,j,k-1)
        eb(i,j,k-1) = 0.0
      else
        eb(i,j,k-1) = eb(i,j,k-1) - ea(i,j,k)
        ea(i,j,k) = 0.0
      endif
    enddo ; enddo ; enddo
  endif


!   This initializes the halo regions of ea, eb, and hold to default values.
  do k=1,nz
    do i=is-1,ie+1
      hold(i,js-1,k) = G%Angstrom ; ea(i,js-1,k) = 0.0 ; eb(i,js-1,k) = 0.0
      hold(i,je+1,k) = G%Angstrom ; ea(i,je+1,k) = 0.0 ; eb(i,je+1,k) = 0.0
    enddo
    do j=js,je
      hold(is-1,j,k) = G%Angstrom ; ea(is-1,j,k) = 0.0 ; eb(is-1,j,k) = 0.0
      hold(ie+1,j,k) = G%Angstrom ; ea(ie+1,j,k) = 0.0 ; eb(ie+1,j,k) = 0.0
    enddo
  enddo
  
  call cpu_clock_begin(id_clock_pass)
  if (G%symmetric) then
    call pass_var(hold,G%Domain,complete=.false.)
    call pass_var(eb,G%Domain,complete=.false.)
    call pass_var(ea,G%Domain)
  else
    call pass_var(hold,G%Domain,To_West+To_South,complete=.false.)
    call pass_var(eb,G%Domain,To_West+To_South,complete=.false.)
    call pass_var(ea,G%Domain,To_West+To_South)
  endif  
  call cpu_clock_end(id_clock_pass)

!  Use a tridiagonal solver to determine the effect of the diapycnal
!  advection on the velocity field.   It is assumed that water leaves
!  or enters the ocean with the surface velocity.

  if (CS%debug) then
    call MOM_state_chksum("before tridiag ", u, v, h, G)
  endif
  call cpu_clock_begin(id_clock_tridiag)
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
    call MOM_state_chksum("aft 1st loop tridiag ", u, v, h, G)
  endif
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
    call MOM_state_chksum("after tridiag ", u, v, h, G)
  endif

!   Frazil formation keeps the temperature above the freezing point.
! make_frazil is deliberately called at both the beginning and at
! the end of the diabatic processes.
  if (ASSOCIATED(T) .AND. ASSOCIATED(tv%frazil)) call make_frazil(h,tv,G,CS)

  if (CS%id_ea > 0) call post_data(CS%id_ea, ea, CS%diag)
  if (CS%id_eb > 0) call post_data(CS%id_eb, eb, CS%diag)
  if (CS%id_dudt_dia > 0) call post_data(CS%id_dudt_dia, ADp%du_dt_dia, CS%diag)
  if (CS%id_dvdt_dia > 0) call post_data(CS%id_dvdt_dia, ADp%dv_dt_dia, CS%diag)
  if (CS%id_wd > 0) call post_data(CS%id_wd, CDp%diapyc_vel, CS%diag)

  if (CS%id_Tdif > 0) call post_data(CS%id_Tdif, Tdif_flx, CS%diag)
  if (CS%id_Tadv > 0) call post_data(CS%id_Tadv, Tadv_flx, CS%diag)
  if (CS%id_Sdif > 0) call post_data(CS%id_Sdif, Sdif_flx, CS%diag)
  if (CS%id_Sadv > 0) call post_data(CS%id_Sadv, Sadv_flx, CS%diag)

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
    call calc_Zint_diags(h, z_ptrs, z_ids, num_z_diags, G, CS%diag_to_Z_CSp)

end subroutine diabatic

subroutine adiabatic(h, tv, fluxes, dt, G, CS)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(inout) :: h
  type(thermo_var_ptrs),                  intent(inout) :: tv
  type(forcing),                          intent(inout) :: fluxes
  real,                                   intent(in)    :: dt
  type(ocean_grid_type),                  intent(inout) :: G
  type(diabatic_CS),                      pointer       :: CS

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    zeros    ! An array of zeros.

  zeros(:,:,:) = 0.0

  call call_tracer_column_fns(h, h, zeros, zeros, fluxes, dt, G, tv, &
                              CS%optics, CS%tracer_flow_CSp)

end subroutine adiabatic

subroutine make_frazil(h, tv, G, CS)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h
  type(thermo_var_ptrs),                 intent(inout) :: tv
  type(ocean_grid_type),                 intent(in)    :: G
  type(diabatic_CS),                     intent(in)    :: CS

!   Frazil formation keeps the temperature above the freezing point.
! This subroutine warms any water that is colder than the (currently
! surface) freezing point up to the freezing point and accumulates
! the required heat (in J m-2) in tv%frazil.
!   The expression, below, for the freezing point of sea water comes
! from Millero (1978) via Appendix A of Gill, 1982.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 diabatic_driver_init.
  real, dimension(SZI_(G)) :: &
    fraz_col, & ! The accumulated heat requirement due to frazil, in J.
    T_freeze    ! The freezing potential temperature at the current salinity, C.
  real, dimension(SZI_(G),SZK_(G)) :: &
    pressure    ! The pressure at the middle of each layer in Pa.
  real :: hc    ! A layer's heat capacity in J m-2 K-1.
  logical :: T_fr_set  ! True if the freezing point has been calculated for a
                       ! row of points.
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call cpu_clock_begin(id_clock_frazil)

  if (.not.CS%pressure_dependent_frazil) then
    do k=1,nz ; do i=is,ie ; pressure(i,k) = 0.0 ; enddo ; enddo
  endif

  do j=js,je
    do i=is,ie ; fraz_col(:) = 0.0 ; enddo

    if (CS%pressure_dependent_frazil) then
      do i=is,ie
        pressure(i,1) = (0.5*G%H_to_Pa)*h(i,j,1)
      enddo
      do k=2,nz ; do i=is,ie
        pressure(i,k) = pressure(i,k-1) + &
          (0.5*G%H_to_Pa) * (h(i,j,k) + h(i,j,k-1))
      enddo ; enddo
    endif

    if (CS%reclaim_frazil) then
      T_fr_set = .false.
      do i=is,ie ; if (tv%frazil(i,j) > 0.0) then
        if (.not.T_fr_set) then
          call calculate_TFreeze(tv%S(i:,j,1), pressure(i:,1), T_freeze(i:), &
                                 1, ie-i+1, tv%eqn_of_state)
          T_fr_set = .true.
        endif

        if (tv%T(i,j,1) > T_freeze(i)) then
    ! If frazil had previously been formed, but the surface temperature is now
    ! above freezing, cool the surface layer with the frazil heat deficit.
          hc = (tv%C_p*G%H_to_kg_m2) * h(i,j,1)
          if (tv%frazil(i,j) - hc * (tv%T(i,j,1) - T_freeze(i)) <= 0.0) then
            tv%T(i,j,1) = tv%T(i,j,1) - tv%frazil(i,j)/hc
            tv%frazil(i,j) = 0.0
          else
            tv%frazil(i,j) = tv%frazil(i,j) - hc * (tv%T(i,j,1) - T_freeze(i))
            tv%T(i,j,1) = T_freeze(i)
          endif
        endif
      endif ; enddo
    endif

    do k=nz,1,-1
      T_fr_set = .false.
      do i=is,ie
        if ((G%mask2dT(i,j) > 0.0) .and. &
            ((tv%T(i,j,k) < 0.0) .or. (fraz_col(i) > 0.0))) then
          if (.not.T_fr_set) then
            call calculate_TFreeze(tv%S(i:,j,k), pressure(i:,k), T_freeze(i:), &
                                   1, ie-i+1, tv%eqn_of_state)
            T_fr_set = .true.
          endif

          hc = (tv%C_p*G%H_to_kg_m2) * h(i,j,k)
          if (h(i,j,k) <= 10.0*G%Angstrom) then
            ! Very thin layers should not be cooled by the frazil flux.
            if (tv%T(i,j,k) < T_freeze(i)) then
              fraz_col(i) = fraz_col(i) + hc * (T_freeze(i) - tv%T(i,j,k))
              tv%T(i,j,k) = T_freeze(i)
            endif
          else
            if (fraz_col(i) + hc * (T_freeze(i) - tv%T(i,j,k)) <= 0.0) then
              tv%T(i,j,k) = tv%T(i,j,k) - fraz_col(i)/hc
              fraz_col(i) = 0.0
            else
              fraz_col(i) = fraz_col(i) + hc * (T_freeze(i) - tv%T(i,j,k))
              tv%T(i,j,k) = T_freeze(i)
            endif
          endif
        endif
      enddo
    enddo
    do i=is,ie
      tv%frazil(i,j) = tv%frazil(i,j) + fraz_col(i)
    enddo
  enddo
  call cpu_clock_end(id_clock_frazil)

end subroutine make_frazil

subroutine double_diffuse_T_S(h, tv, visc, dt, G)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h
  type(thermo_var_ptrs),                 intent(inout) :: tv
  type(vertvisc_type),                   intent(in)    :: visc
  real,                                  intent(in)    :: dt
  type(ocean_grid_type),                 intent(in)    :: G

! This subroutine applies double diffusion to T & S, assuming no diapycal mass
! fluxes, using a simple triadiagonal solver.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      visc - A structure containing vertical viscosities, bottom boundary
!                   layer properies, and related fields.
!  (in)      dt - Time increment, in s.
!  (in)      G - The ocean's grid structure.

  real, dimension(SZI_(G)) :: &
    b1_T, b1_S, &  !  Variables used by the tridiagonal solvers of T & S, in H.
    d1_T, d1_S     !  Variables used by the tridiagonal solvers, nondim.
  real, dimension(SZI_(G),SZK_(G)) :: &
    c1_T, c1_S     !  Variables used by the tridiagonal solvers, in m or kg m-2.
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    mix_T, mix_S   !  Mixing distances in both directions across each
                   !  interface, in m or kg m-2.
  real :: h_tr         ! h_tr is h at tracer points with a tiny thickness
                       ! added to ensure positive definiteness, in m or kg m-2.
  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected, in m or kg m-2.
  real :: I_h_int      ! The inverse of the thickness associated with an
                       ! interface, in m-1 or m2 kg-1.
  real :: b_denom_T    ! The first term in the denominators for the expressions
  real :: b_denom_S    ! for b1_T and b1_S, both in m or kg m-2.

  integer :: i, j, k, is, ie, js, je, nz
  real, pointer :: T(:,:,:), S(:,:,:), Kd_T(:,:,:), Kd_S(:,:,:)
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  h_neglect = G%H_subroundoff

  if (.not.associated(tv%T)) call MOM_error(FATAL, &
      "double_diffuse_T_S: Called with an unassociated tv%T")
  if (.not.associated(tv%S)) call MOM_error(FATAL, &
      "double_diffuse_T_S: Called with an unassociated tv%S")
  if (.not.associated(visc%Kd_extra_T)) call MOM_error(FATAL, &
      "double_diffuse_T_S: Called with an unassociated visc%Kd_extra_T")
  if (.not.associated(visc%Kd_extra_S)) call MOM_error(FATAL, &
      "double_diffuse_T_S: Called with an unassociated visc%Kd_extra_S")

  T => tv%T ; S => tv%S
  Kd_T => visc%Kd_extra_T ; Kd_S => visc%Kd_extra_S
  do j=js,je
    do i=is,ie
      I_h_int = 1.0 / (0.5 * (h(i,j,1) + h(i,j,2)) + h_neglect)
      mix_T(i,2) = ((dt * Kd_T(i,j,2)) * G%m_to_H**2) * I_h_int
      mix_S(i,2) = ((dt * Kd_S(i,j,2)) * G%m_to_H**2) * I_h_int
  
      h_tr = h(i,j,1) + h_neglect
      b1_T(i) = 1.0 / (h_tr + mix_T(i,2))
      b1_S(i) = 1.0 / (h_tr + mix_S(i,2))
      d1_T(i) = h_tr * b1_T(i)
      d1_S(i) = h_tr * b1_S(i)
      T(i,j,1) = (b1_T(i)*h_tr)*T(i,j,1)
      S(i,j,1) = (b1_S(i)*h_tr)*S(i,j,1)
    enddo
    do k=2,nz-1 ; do i=is,ie
      ! Calculate the mixing across the interface below this layer.
      I_h_int = 1.0 / (0.5 * (h(i,j,k) + h(i,j,k+1)) + h_neglect)
      mix_T(i,K+1) = ((dt * Kd_T(i,j,K+1)) * G%m_to_H**2) * I_h_int
      mix_S(i,K+1) = ((dt * Kd_S(i,j,K+1)) * G%m_to_H**2) * I_h_int

      c1_T(i,k) = mix_T(i,K) * b1_T(i)
      c1_S(i,k) = mix_S(i,K) * b1_S(i)

      h_tr = h(i,j,k) + h_neglect
      b_denom_T = h_tr + d1_T(i)*mix_T(i,K)
      b_denom_S = h_tr + d1_S(i)*mix_S(i,K)
      b1_T(i) = 1.0 / (b_denom_T + mix_T(i,K+1))
      b1_S(i) = 1.0 / (b_denom_S + mix_S(i,K+1))
      d1_T(i) = b_denom_T * b1_T(i)
      d1_S(i) = b_denom_S * b1_S(i)

      T(i,j,k) = b1_T(i) * (h_tr*T(i,j,k) + mix_T(i,K)*T(i,j,k-1))
      S(i,j,k) = b1_S(i) * (h_tr*S(i,j,k) + mix_S(i,K)*S(i,j,k-1))
    enddo ; enddo
    do i=is,ie
      c1_T(i,nz) = mix_T(i,nz) * b1_T(i)
      c1_S(i,nz) = mix_S(i,nz) * b1_S(i)

      h_tr = h(i,j,nz) + h_neglect
      b1_T(i) = 1.0 / (h_tr + d1_T(i)*mix_T(i,nz))
      b1_S(i) = 1.0 / (h_tr + d1_S(i)*mix_S(i,nz))

      T(i,j,nz) = b1_T(i) * (h_tr*T(i,j,nz) + mix_T(i,nz)*T(i,j,nz-1))
      S(i,j,nz) = b1_S(i) * (h_tr*S(i,j,nz) + mix_S(i,nz)*S(i,j,nz-1))
    enddo
    do k=nz-1,1,-1 ; do i=is,ie
      T(i,j,k) = T(i,j,k) + c1_T(i,k+1)*T(i,j,k+1)
      S(i,j,k) = S(i,j,k) + c1_S(i,k+1)*S(i,j,k+1)
    enddo ; enddo
  enddo

end subroutine double_diffuse_T_S

subroutine adjust_salt(h, tv, G, CS)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h
  type(thermo_var_ptrs),                 intent(inout) :: tv
  type(ocean_grid_type),                 intent(in)    :: G
  type(diabatic_CS),                     intent(in)    :: CS
  
!  Keep salinity from falling below a small but positive threshold
!  This occurs when the ice model attempts to extract more salt then
!  is actually available to it from the ocean.  

! Arguments: h - Layer thickness, in m.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 diabatic_driver_init.
  real :: salt_add_col(SZI_(G),SZJ_(G)) ! The accumulated salt requirement
  real :: S_min      ! The minimum salinity
  real :: mc         ! A layer's mass kg  m-2 .
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

!  call cpu_clock_begin(id_clock_adjust_salt)

  S_min = 0.01
  
  salt_add_col(:,:) = 0.0

  do k=nz,1,-1 ; do j=js,je ; do i=is,ie
    if ((G%mask2dT(i,j) > 0.0) .and. &
         ((tv%S(i,j,k) < S_min) .or. (salt_add_col(i,j) > 0.0))) then
      mc = G%H_to_kg_m2 * h(i,j,k)
      if (h(i,j,k) <= 10.0*G%Angstrom) then
        ! Very thin layers should not be adjusted by the salt flux
        if (tv%S(i,j,k) < S_min) then
          salt_add_col(i,j) = salt_add_col(i,j) +  mc * (S_min - tv%S(i,j,k))
          tv%S(i,j,k) = S_min
        endif
      else
        if (salt_add_col(i,j) + mc * (S_min - tv%S(i,j,k)) <= 0.0) then
          tv%S(i,j,k) = tv%S(i,j,k) - salt_add_col(i,j)/mc
          salt_add_col(i,j) = 0.0
        else
          salt_add_col(i,j) = salt_add_col(i,j) + mc * (S_min - tv%S(i,j,k))
          tv%S(i,j,k) = S_min
        endif
      endif
    endif
  enddo ; enddo ; enddo
  do j=js,je ; do i=is,ie
    tv%salt_deficit(i,j) = tv%salt_deficit(i,j) + salt_add_col(i,j)
  enddo ; enddo
!  call cpu_clock_end(id_clock_adjust_salt)

end subroutine adjust_salt

subroutine adiabatic_driver_init(Time, G, param_file, diag, CS, &
                                tracer_flow_CSp, diag_to_Z_CSp)
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(diabatic_CS),       pointer       :: CS
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp
  type(diag_to_Z_CS),      pointer       :: diag_to_Z_CSp
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      tracer_flow_CSp - A pointer to the control structure of the tracer
!                              flow control module.
!  (in)      diag_to_Z_CSp - A pointer to the Z-diagnostics control structure.

!   This is a highly simplified version of diabatic_driver_init that will allow
! the tracer column functions to be called without allowing any of the diabatic
! processes to be used.
! This include declares and sets the variable "version".
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

subroutine diabatic_driver_init(Time, G, param_file, useALEalgorithm, diag, &
                     ADp, CDp, CS, tracer_flow_CSp, sponge_CSp, diag_to_Z_CSp)
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  logical,                 intent(in)    :: useALEalgorithm
  type(diag_ctrl), target, intent(inout) :: diag
  type(accel_diag_ptrs),   intent(inout) :: ADp
  type(cont_diag_ptrs),    intent(inout) :: CDp
  type(diabatic_CS),       pointer       :: CS
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp
  type(sponge_CS),         pointer       :: sponge_CSp
  type(diag_to_Z_CS),      pointer       :: diag_to_Z_CSp
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (inout)   ADp - A structure with pointers to the various accelerations in
!                  the momentum equations, to enable the later calculation
!                  of derived diagnostics, like energy budgets.
!  (inout)   CDp - A structure with pointers to various terms in the continuity
!                  equations.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      tracer_flow_CSp - A pointer to the control structure of the tracer
!                              flow control module.
!  (in)      sponge_CSp - A pointer to the sponge module control structure.
!  (in)      diag_to_Z_CSp - A pointer to the Z-diagnostics control structure.
  real :: Kd
  logical :: use_temperature, double_diffusion
  type(vardesc) :: vd
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod  = "MOM_diabatic_driver" ! This module's name.
  character(len=48)  :: thickness_units
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz, nbands
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "diabatic_driver_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag
  if (associated(tracer_flow_CSp)) CS%tracer_flow_CSp => tracer_flow_CSp
  if (associated(sponge_CSp)) CS%sponge_CSp => sponge_CSp
  if (associated(diag_to_Z_CSp)) CS%diag_to_Z_CSp => diag_to_Z_CSp

  CS%useALEalgorithm = useALEalgorithm
  CS%bulkmixedlayer = (G%nkml > 0)

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
  call get_param(param_file, mod, "DOUBLE_DIFFUSION", double_diffusion, &
                 "If true, apply parameterization of double-diffusion.", &
                 default=.false. )
  call get_param(param_file, mod, "USE_JACKSON_PARAM", CS%use_kappa_shear, &
                 "If true, use the Jackson-Hallberg-Legg (JPO 2008) \n"//& 
                 "shear mixing parameterization.", default=.false.)
  if (CS%bulkmixedlayer) then
    call get_param(param_file, mod, "ML_MIX_FIRST", CS%ML_mix_first, &
                 "The fraction of the mixed layer mixing that is applied \n"//&
                 "before interior diapycnal mixing.  0 by default.", &
                 units="nondim", default=0.0)
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
                 "If true, use the code that advances as separate set of \n"//&
                 "equations for the internal tide energy density.", default=.false.)
  call get_param(param_file, mod, "MASSLESS_MATCH_TARGETS", &
                                CS%massless_match_targets, &
                 "If true, the temperature and salinity of massless layers \n"//&
                 "are kept consistent with their target densities. \n"//&
                 "Otherwise the properties of massless layers evolve \n"//&
                 "diffusively to match massive neighboring layers.", &
                 default=.true.)
  call get_param(param_file, mod, "RECLAIM_FRAZIL", CS%reclaim_frazil, &
                 "If true, try to use any frazil heat deficit to cool any\n"//&
                 "overlying layers down to the freezing point, thereby \n"//&
                 "avoiding the creation of thin ice when the SST is above \n"//&
                 "the freezing point.", default=.true.)
  call get_param(param_file, mod, "PRESSURE_DEPENDENT_FRAZIL", &
                                CS%pressure_dependent_frazil, &
                 "If true, use a pressure dependent freezing temperature \n"//&
                 "when making frazil. The default is false, which will be \n"//&
                 "faster but is inappropriate with ice-shelf cavities.", &
                 default=.false.)
  call get_param(param_file, mod, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)
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

  if (G%Boussinesq) then ; thickness_units = "meter"
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

  if (associated(diag_to_Z_CSp)) then
    vd = vardesc("Kd_z","Diapycnal diffusivity at interfaces, interpolated to z",&
                 'h','z','s',"meter2 second-1")
    CS%id_Kd_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = vardesc("Tflx_dia_dif_z","Diffusive diapycnal temperature flux across interfaces, interpolated to z",&
                 'h','z','s',"degC meter second-1")
    CS%id_Tdif_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = vardesc("Tflx_dia_adv_z","Advective diapycnal temperature flux across interfaces, interpolated to z",&
                 'h','z','s',"degC meter second-1")
    CS%id_Tadv_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = vardesc("Sflx_dia_dif_z","Diffusive diapycnal salinity flux across interfaces, interpolated to z",&
                 'h','z','s',"PSU meter second-1")
    CS%id_Sdif_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = vardesc("Sflx_dia_adv_z","Advective diapycnal salinity flux across interfaces, interpolated to z",&
                 'h','z','s',"PSU meter second-1")
    CS%id_Sadv_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
  endif

  if (CS%id_dudt_dia > 0) call safe_alloc_ptr(ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
  if (CS%id_dvdt_dia > 0) call safe_alloc_ptr(ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
  if (CS%id_wd > 0) call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)

  call set_diffusivity_init(Time, G, param_file, diag, CS%set_diff_CSp, diag_to_Z_CSp)
  call entrain_diffusive_init(Time, G, param_file, diag, CS%entrain_diffusive_CSp)
  if (CS%use_geothermal) &
    call geothermal_init(Time, G, param_file, diag, CS%geothermal_CSp)
  if (CS%use_kappa_shear) &
    call kappa_shear_init(Time, G, param_file, diag, CS%kappa_shear_CSp)

  if (CS%use_int_tides) then
    call int_tide_input_init(Time, G, param_file, diag, CS%int_tide_input_CSp, &
                             CS%int_tide_input)
    call internal_tides_init(Time, G, param_file, diag, CS%int_tide_CSp)
  endif

  id_clock_entrain = cpu_clock_id('(Ocean diabatic entrain)', grain=CLOCK_MODULE)
  if (CS%bulkmixedlayer) &
    id_clock_mixedlayer = cpu_clock_id('(Ocean mixed layer)', grain=CLOCK_MODULE)
  id_clock_uv_at_h = cpu_clock_id('(Ocean find_uv_at_h)', grain=CLOCK_ROUTINE)
  id_clock_remap = cpu_clock_id('(Ocean vert remap)', grain=CLOCK_MODULE)
  if (use_temperature) &
    id_clock_frazil = cpu_clock_id('(Ocean frazil)', grain=CLOCK_ROUTINE)
  if (CS%use_geothermal) &
    id_clock_geothermal = cpu_clock_id('(Ocean geothermal)', grain=CLOCK_ROUTINE)
  if (CS%use_kappa_shear) &
    id_clock_kappa_shear = cpu_clock_id('(Ocean kappa_shear)', grain=CLOCK_MODULE)
  id_clock_set_diffusivity = cpu_clock_id('(Ocean set_diffusivity)', grain=CLOCK_MODULE)
  id_clock_tracers = cpu_clock_id('(Ocean tracer_columns)', grain=CLOCK_MODULE_DRIVER+5)
  if (CS%use_sponge) &
    id_clock_sponge = cpu_clock_id('(Ocean sponges)', grain=CLOCK_MODULE)
  id_clock_tridiag = cpu_clock_id('(Ocean diabatic tridiag)', grain=CLOCK_ROUTINE)
  id_clock_pass = cpu_clock_id('(Ocean diabatic message passing)', grain=CLOCK_ROUTINE)
  id_clock_double_diff = -1 ; if (double_diffusion) &
    id_clock_double_diff = cpu_clock_id('(Ocean double diffusion)', grain=CLOCK_ROUTINE)

  if (CS%bulkmixedlayer) &
    call bulkmixedlayer_init(Time, G, param_file, diag, CS%bulkmixedlayer_CSp)

  call regularize_layers_init(Time, G, param_file, diag, CS%regularize_layers_CSp)

  if (use_temperature) then
    call get_param(param_file, mod, "PEN_SW_NBANDS", nbands, default=1)
    if (nbands > 0) then
      allocate(CS%optics)
      call opacity_init(Time, G, param_file, diag, CS%tracer_flow_CSp, CS%opacity_CSp, CS%optics)
    endif
  endif

end subroutine diabatic_driver_init

subroutine diabatic_driver_end(CS)
  type(diabatic_CS), pointer :: CS

  call entrain_diffusive_end(CS%entrain_diffusive_CSp)
  call set_diffusivity_end(CS%set_diff_CSp)

  if (associated(CS%optics)) then
    call opacity_end(CS%opacity_CSp, CS%optics)
    deallocate(CS%optics)
  endif
  if (associated(CS)) deallocate(CS)

end subroutine diabatic_driver_end

subroutine find_uv_at_h(u, v, h, u_h, v_h, G, ea, eb)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)  :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)  :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)  :: h
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(out) :: u_h, v_h
  type(ocean_grid_type),                  intent(in)  :: G
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in), optional  :: ea, eb
!   This subroutine calculates u_h and v_h (velocities at thickness
! points), optionally using the entrainments (in m) passed in as arguments.

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (out)     u_h - The zonal velocity at thickness points after
!                  entrainment, in m s-1.
!  (out)     v_h - The meridional velocity at thickness points after
!                  entrainment, in m s-1.
!  (in)      G - The ocean's grid structure.
!  (in, opt) ea - The amount of fluid entrained from the layer above within
!                 this time step, in units of m or kg m-2.  Omitting ea is the
!                 same as setting it to 0.
!  (in, opt) eb - The amount of fluid entrained from the layer below within
!                 this time step, in units of m or kg m-2.  Omitting eb is the
!                 same as setting it to 0.  ea and eb must either be both
!                 present or both absent.

  real :: b_denom_1    ! The first term in the denominator of b1 in m or kg m-2.
  real :: h_neglect    ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected, in m or kg m-2.
  real :: b1(SZI_(G)), d1(SZI_(G)), c1(SZI_(G),SZK_(G))
  real :: a_n(SZI_(G)), a_s(SZI_(G))  ! Fractional weights of the neighboring
  real :: a_e(SZI_(G)), a_w(SZI_(G))  ! velocity points, ~1/2 in the open
                                      ! ocean, nondimensional.
  real :: s, Idenom
  logical :: mix_vertically
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  call cpu_clock_begin(id_clock_uv_at_h)
  h_neglect = G%H_subroundoff

  mix_vertically = present(ea)
  if (present(ea) .neqv. present(eb)) call MOM_error(FATAL, &
      "find_uv_at_h: Either both ea and eb or neither one must be present "// &
      "in call to find_uv_at_h.")

  do j=js,je
    do i=is,ie
      s = G%areaCu(i-1,j)+G%areaCu(i,j)
      if (s>0.0) then
        Idenom = sqrt(0.5*G%IareaT(i,j)/s)
        a_w(i) = G%areaCu(i-1,j)*Idenom
        a_e(i) = G%areaCu(i,j)*Idenom
      else
        a_w(i) = 0.0 ; a_e(i) = 0.0
      endif

      s = G%areaCv(i,j-1)+G%areaCv(i,j)
      if (s>0.0) then
        Idenom = sqrt(0.5*G%IareaT(i,j)/s)
        a_s(i) = G%areaCv(i,j-1)*Idenom
        a_n(i) = G%areaCv(i,j)*Idenom
      else
        a_s(i) = 0.0 ; a_n(i) = 0.0
      endif
    enddo

    if (mix_vertically) then
      do i=is,ie
        b_denom_1 = h(i,j,1) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,1))
        d1(i) = b_denom_1 * b1(i)
        u_h(i,j,1) = (h(i,j,1)*b1(i)) * (a_e(i)*u(i,j,1) + a_w(i)*u(i-1,j,1))
        v_h(i,j,1) = (h(i,j,1)*b1(i)) * (a_n(i)*v(i,j,1) + a_s(i)*v(i,j-1,1))
      enddo
      do k=2,nz ; do i=is,ie
        c1(i,k) = eb(i,j,k-1) * b1(i)
        b_denom_1 = h(i,j,k) + d1(i)*ea(i,j,k) + h_neglect
        b1(i) = 1.0 / (b_denom_1 + eb(i,j,k))
        d1(i) = b_denom_1 * b1(i)
        u_h(i,j,k) = (h(i,j,k) * (a_e(i)*u(i,j,k) + a_w(i)*u(i-1,j,k)) + &
                      ea(i,j,k)*u_h(i,j,k-1))*b1(i)
        v_h(i,j,k) = (h(i,j,k) * (a_n(i)*v(i,j,k) + a_s(i)*v(i,j-1,k)) + &
                      ea(i,j,k)*v_h(i,j,k-1))*b1(i)
      enddo ; enddo
      do k=nz-1,1,-1 ; do i=is,ie
        u_h(i,j,k) = u_h(i,j,k) + c1(i,k+1)*u_h(i,j,k+1)
        v_h(i,j,k) = v_h(i,j,k) + c1(i,k+1)*v_h(i,j,k+1)
      enddo ; enddo
    else
      do k=1,nz ; do i=is,ie
        u_h(i,j,k) = a_e(i)*u(i,j,k) + a_w(i)*u(i-1,j,k)
        v_h(i,j,k) = a_n(i)*v(i,j,k) + a_s(i)*v(i,j-1,k)
      enddo ; enddo
    endif
  enddo

  call cpu_clock_end(id_clock_uv_at_h)
end subroutine find_uv_at_h

subroutine MOM_state_chksum(mesg, u, v, h, G)
  character(len=*),                       intent(in) :: mesg
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: v
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in) :: h
  type(ocean_grid_type),                  intent(in) :: G
!   This subroutine writes out chksums for the model's basic state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m.
!  (in)      uh - Volume flux through zonal faces = u*h*dy, m3 s-1.
!  (in)      vh - Volume flux through meridional faces = v*h*dx, in m3 s-1.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  call uchksum(u,mesg//" u",G,haloshift=0)
  call vchksum(v,mesg//" v",G,haloshift=0)
  call hchksum(G%H_to_m*h, mesg//" h",G,haloshift=0)
end subroutine MOM_state_chksum

subroutine applyBoundaryFluxes(G, dt, fluxes, optics, ea, eb, h, tv)
  type(ocean_grid_type),                 intent(in)    :: G
  real,                                  intent(in)    :: dt
  type(forcing),                         intent(in)    :: fluxes
  type(optics_type),                     pointer       :: optics
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: ea, eb, h
  type(thermo_var_ptrs),                 intent(inout) :: tv
  !   This subroutine applies thermodynamic forcing (contained in fluxes type)
  ! to h, tv%T and tv%S. If the event the P-E+R is negatively in excess of the
  ! available layer thickness, ea and eb might be adjusted to exchange (entrain)
  ! volume between layers to maintain non-negative values.
  real :: Irho0, I_Cp, Irho_cp, H_limit_fluxes, IforcingDepthScale
  real :: dThickness, dTemp, dSalt, fractionOfForcing, hOld, Ithickness
  real, dimension(SZI_(G)) :: netThickness, netHeat, netSalt, htot, Ttot
  real, dimension(SZI_(G), SZK_(G)) :: h2d, T2d, eps
  integer, dimension(SZI_(G), SZK_(G)) :: ksort
  real, allocatable, dimension(:,:) :: Pen_SW_bnd
  real, allocatable, dimension(:,:,:) :: opacityBand
  logical :: use_riverHeatContent, useCalvingHeatContent
  integer :: i, j, is, ie, js, je, k, nz, n, nsw
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  ! Skip applying forcing if fluxes%sw is not associated.
  if (.not.ASSOCIATED(fluxes%sw)) return

  nsw = 0
  if (ASSOCIATED(optics)) nsw = optics%nbands
  allocate(Pen_SW_bnd(max(nsw,1),SZI_(G)))
  allocate(opacityBand(max(nsw,1),SZI_(G),SZK_(G)))

  Irho0 = 1.0 / G%Rho0
  I_Cp = 1.0 / fluxes%C_p
  Irho_cp = 1.0 / (G%H_to_kg_m2 * fluxes%C_p)
  ! H_limit_fluxes is used by extractFluxes1d to scale away fluxes if the total
  ! depth of the ocean is vanishing. It does not (yet) handle a value of zero!
  H_limit_fluxes = max(G%Angstrom, 1.E-30) ! This is a hack to avoid division by zero
  ! To accomodate vanishing upper layers, we need to allow for an instantaneous
  ! distribution of forcing over some finite vertical extent. The bulk mixed layer
  ! code handled this properly. The inverse scale, IforcingDepthScale, is a hack which 
  ! should not be tickled in Eulerian mode. It stops all the forcing going into a
  ! vanish(ed/ing) layer.
  IforcingDepthScale = 1000. ! Use 1 mm to distribute the surface fluxes uniformly

  use_riverHeatContent = .false.
  useCalvingHeatContent = .false.

  ! s/r aborbRemaining uses an indirect indexing in the vertical, a hold over for use
  ! with the bulk mixed layer
  do k=1,nz ; do i=is,ie
      ksort(i,k) = k
  enddo ; enddo

  do j=js,je ! Work in vertical slices (this is a hold over from the routines called with a j argument)
    ! Copy state into 2D-slice arrays
    do k=1,nz
      do i=is,ie
        h2d(i,k) = h(i,j,k)
        T2d(i,k) = tv%T(i,j,k)
        do n=1,nsw
          opacityBand(n,i,k) = G%H_to_m*optics%opacity_band(n,i,j,k)
        enddo
        eps(i,k) = 0.
      enddo
    enddo
    do i=is,ie ; htot(i) = 0. ; enddo

    ! The surface forcing is contained in the fluxes type.
    ! Here, we "unpack" it and aggregate all the thermodynamic forcing into
    ! netThickness, netHeat, netSalt and the SW penetrative componennts, Pen_SW_bnd.
    call extractFluxes1d(G, fluxes, optics, nsw, j, dt, &
                  H_limit_fluxes, use_riverHeatContent, useCalvingHeatContent, &
                  h2d, T2d, netThickness, netHeat, netSalt, Pen_SW_bnd, tv)
 
    do i=is,ie
      do k=1,nz
        ! Fraction of forcing that we put into this layer is normally 100%, unless the
        ! layer is thin relative to 1/IforcingDepthScale
        fractionOfForcing = min(1.0, h2d(i,k)*IforcingDepthScale)

        ! Change in state due to forcing
        ! (Limit mass loss to the available mass in layer)
        dThickness = max( fractionOfForcing*netThickness(i), -h2d(i,k) )
        dTemp = fractionOfForcing*netHeat(i)
        dSalt = fractionOfForcing*netSalt(i)

        ! Update the forcing by the part to be consumed
        netThickness(i) = netThickness(i) - dThickness
        netHeat(i) = netHeat(i) - dTemp
        netSalt(i) = netSalt(i) - dSalt

        ! Adjust heating by the temperature of rain/water vapor
        dTemp = dTemp + dThickness*tv%T(i,j,k)

        ! Update state by the appropriate delta (change in state calculated above)
        hOld = h2d(i,k) ! Need to keep original thickness in hand
        h2d(i,k) = h2d(i,k) + dThickness
        Ithickness = 1./h2d(i,k)
        tv%T(i,j,k) = (hOld*tv%T(i,j,k) + dTemp)*Ithickness
        tv%S(i,j,k) = (hOld*tv%S(i,j,k) + dSalt)*Ithickness

      enddo ! k
    enddo ! i

    ! Heat by the divergence of penetrating SW (this uses the updated thicknesses)
    call absorbRemainingSW(G, h2d, eps, htot, opacityBand, nsw, j, dt, &
                           H_limit_fluxes, .true., .true., &
                           ksort, T2d, Ttot, Pen_SW_bnd)

    ! Copy slice back into model state
    do k=1,nz
      do i=is,ie
        h(i,j,k) = h2d(i,k)
        tv%T(i,j,k) = T2d(i,k)
      enddo
    enddo
  enddo ! j

  deallocate(Pen_SW_bnd)
  deallocate(opacityBand)
  
end subroutine applyBoundaryFluxes

end module MOM_diabatic_driver
