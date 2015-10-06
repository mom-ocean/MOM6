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
!*  By Robert Hallberg, Alistair Adcroft, and Stephen Griffies         *
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

use MOM_bulk_mixed_layer,    only : bulkmixedlayer, bulkmixedlayer_init, bulkmixedlayer_CS
use MOM_checksums,           only : hchksum, uchksum, vchksum
use MOM_checksum_packages,   only : MOM_state_chksum, MOM_state_stats
use MOM_cpu_clock,           only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,           only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diabatic_aux,        only : diabatic_aux_init, diabatic_aux_end, diabatic_aux_CS
use MOM_diabatic_aux,        only : make_frazil, adjust_salt, insert_brine, differential_diffuse_T_S, triDiagTS
use MOM_diabatic_aux,        only : find_uv_at_h, diagnoseMLDbyDensityDifference, applyBoundaryFluxesInOut
use MOM_diag_mediator,       only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,       only : diag_ctrl, time_type, diag_update_target_grids
use MOM_diag_mediator,       only : diag_ctrl, query_averaging_enabled !BDM: see src/framework/MOM_diag_mediator.F90
use MOM_diag_to_Z,           only : diag_to_Z_CS, register_Zint_diag, calc_Zint_diags
use MOM_diffConvection,      only : diffConvection_CS, diffConvection_init
use MOM_diffConvection,      only : diffConvection_calculate, diffConvection_end
use MOM_domains,             only : pass_var, To_West, To_South
use MOM_domains,             only : create_group_pass, do_group_pass, group_pass_type
use MOM_energetic_PBL,       only : energetic_PBL, energetic_PBL_init
use MOM_energetic_PBL,       only : energetic_PBL_end, energetic_PBL_CS
use MOM_entrain_diffusive,   only : entrainment_diffusive, entrain_diffusive_init
use MOM_entrain_diffusive,   only : entrain_diffusive_end, entrain_diffusive_CS
use MOM_EOS,                 only : calculate_density, calculate_2_densities, calculate_TFreeze
use MOM_EOS,                 only : calculate_specific_vol_derivs
use MOM_error_handler,       only : MOM_error, FATAL, WARNING, callTree_showQuery, is_root_pe
use MOM_error_handler,       only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,         only : get_param, log_version, param_file_type
use MOM_forcing_type,        only : forcing, MOM_forcing_chksum
use MOM_forcing_type,        only : extractFluxes1d, calculateBuoyancyFlux2d
use MOM_forcing_type,        only : forcing_SinglePointPrint
use MOM_geothermal,          only : geothermal, geothermal_init, geothermal_end, geothermal_CS
use MOM_grid,                only : ocean_grid_type
use MOM_io,                  only : vardesc
use MOM_int_tide_input,      only : set_int_tide_input, int_tide_input_init
use MOM_int_tide_input,      only : int_tide_input_end, int_tide_input_CS, int_tide_input_type
use MOM_internal_tides,      only : propagate_int_tide, register_int_tide_restarts
use MOM_internal_tides,      only : internal_tides_init, internal_tides_end, int_tide_CS
use MOM_kappa_shear,         only : kappa_shear_is_used
use MOM_KPP,                 only : KPP_CS, KPP_init, KPP_calculate, KPP_end
use MOM_KPP,                 only : KPP_NonLocalTransport_temp, KPP_NonLocalTransport_saln
use MOM_opacity,             only : opacity_init, set_opacity, opacity_end, opacity_CS
use MOM_set_diffusivity,     only : set_diffusivity, set_BBL_TKE
use MOM_set_diffusivity,     only : set_diffusivity_init, set_diffusivity_end
use MOM_set_diffusivity,     only : set_diffusivity_CS
use MOM_shortwave_abs,       only : absorbRemainingSW, optics_type
use MOM_sponge,              only : apply_sponge, sponge_CS
use MOM_thickness_diffuse,   only : vert_fill_TS ! BDM
use MOM_time_manager,        only : operator(<=), time_type !BDM: see src/framework/MOM_time_manager.F90
use MOM_tracer_flow_control, only : call_tracer_column_fns, tracer_flow_control_CS
use MOM_variables,           only : thermo_var_ptrs, vertvisc_type, accel_diag_ptrs
use MOM_variables,           only : cont_diag_ptrs, MOM_thermovar_chksum, p3d
use MOM_regularize_layers,   only : regularize_layers, regularize_layers_init, regularize_layers_CS
use time_manager_mod,        only : increment_time !BDM: see shared/time_manager/time_manager.F90
use MOM_wave_speed,          only : wave_speed

implicit none ; private

#include <MOM_memory.h>

public diabatic, diabatic_driver_init, diabatic_driver_end
public adiabatic, adiabatic_driver_init

type, public :: diabatic_CS ; private
  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used with
                             ! nkml sublayers (and additional buffer layers).
  logical :: use_energetic_PBL ! If true, use the implicit energetics planetary
                             ! boundary layer scheme to determine the diffusivity
                             ! in the surface boundary layer.
  logical :: use_kappa_shear ! If true, use the kappa_shear module to find the
                             ! shear-driven diapycnal diffusivity.
  logical :: use_sponge      ! If true, sponges may be applied anywhere in the
                             ! domain.  The exact location and properties of
                             ! those sponges are set by calls to
                             ! initialize_sponge and set_up_sponge_field.
  logical :: use_geothermal  ! If true, apply geothermal heating.
  logical :: use_int_tides   ! If true, use the code that advances a separate set
                             ! of equations for the internal tide energy density.
  logical :: int_tide_source_test
                             ! If true, apply an arbitrary generation site 
                             ! for internal tide testing (BDM) 
  real    :: int_tide_source_x 
                             ! X Location of generation site 
                             ! for internal tide for testing (BDM)
  real    :: int_tide_source_y 
                             ! Y Location of generation site 
                             ! for internal tide for testing (BDM)
  integer :: tlen_days         
                             ! Time interval from start for adding wave source for testing (BDM)
  logical :: uniform_cg	     ! If true, set cg = cg_test everywhere for test case (BDM)
  real    :: cg_test         ! Uniform group velocity of internal tide for test case (BDM)
  type(time_type) :: time_end, time_max_source 
                             ! For use in testing internal tides (BDM)
  logical :: useALEalgorithm ! If true, use the ALE algorithm rather than layered
                             ! isopycnal/stacked shallow water mode. This logical
                             ! passed by argument to diabatic_driver_init.
  logical :: aggregate_FW_forcing ! Determines whether net incoming/outgoing surface
                             ! FW fluxes are applied separately or combined before
                             ! being applied.
  real :: ML_mix_first       ! The nondimensional fraction of the mixed layer
                             ! algorithm that is applied before diffusive mixing.
                             ! The default is 0, while 0.5 gives Strang splitting
                             ! and 1 is a sensible value too.  Note that if there
                             ! are convective instabilities in the initial state,
                             ! the first call may do much more than the second.
  integer :: NKBL            ! The number of buffer layers (if bulk_mixed_layer)
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

  logical :: useKPP          ! use CVmix/KPP diffusivities and non-local transport
  logical :: salt_reject_below_ML ! If true, add salt below mixed layer (layer mode only)
  logical :: KPPisPassive    ! If true, KPP is in passive mode, not changing answers.
  logical :: useConvection   ! If true, calculate large diffusivities when column
                             ! is statically unstable.
  logical :: debug                 ! If true, write verbose checksums for debugging purposes.
  logical :: debugConservation     ! If true, monitor conservation and extrema.
  type(diag_ctrl), pointer :: diag ! structure used to regulate timing of diagnostic output
  real :: MLDdensityDifference ! Density difference used to determine MLD_user
  integer :: nsw               ! SW_NBANDS

  integer :: id_dudt_dia = -1, id_dvdt_dia = -1, id_wd           = -1
  integer :: id_cg1      = -1  !(BDM)
  integer :: id_ea       = -1, id_eb       = -1, id_Kd_z         = -1
  integer :: id_Kd_heat  = -1, id_Kd_salt  = -1, id_Kd_interface = -1, id_Kd_ePBL  = -1
  integer :: id_Tdif_z   = -1, id_Tadv_z   = -1, id_Sdif_z       = -1, id_Sadv_z   = -1
  integer :: id_Tdif     = -1, id_Tadv     = -1, id_Sdif         = -1, id_Sadv     = -1
  integer :: id_subMLN2  = -1, id_brine_lay = -1
  integer :: id_MLD_003  = -1, id_MLD_0125 = -1, id_MLD_user     = -1, id_mlotstsq = -1

  type(diabatic_aux_CS),        pointer :: diabatic_aux_CSp => NULL()
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
  type(tracer_flow_control_CS), pointer :: tracer_flow_CSp       => NULL()
  type(optics_type),            pointer :: optics                => NULL()
  type(diag_to_Z_CS),           pointer :: diag_to_Z_CSp         => NULL()
  type(KPP_CS),                 pointer :: KPP_CSp               => NULL()
  type(diffConvection_CS),      pointer :: Conv_CSp              => NULL()

  type(group_pass_type) :: pass_hold_eb_ea ! For group halo pass

  ! Data arrays for communicating between components
  real, allocatable, dimension(:,:,:) :: KPP_NLTheat    ! KPP non-local transport for heat (m/s)
  real, allocatable, dimension(:,:,:) :: KPP_NLTscalar  ! KPP non-local transport for scalars (m/s)
  real, allocatable, dimension(:,:,:) :: KPP_buoy_flux  ! KPP forcing buoyancy flux (m^2/s^3)
  real, allocatable, dimension(:,:)   :: KPP_temp_flux  ! KPP effective temperature flux (K m/s)
  real, allocatable, dimension(:,:)   :: KPP_salt_flux  ! KPP effective salt flux (ppt m/s)

end type diabatic_CS

integer :: id_clock_entrain, id_clock_mixedlayer, id_clock_set_diffusivity
integer :: id_clock_tracers, id_clock_tridiag, id_clock_pass, id_clock_sponge
integer :: id_clock_geothermal, id_clock_differential_diff, id_clock_remap
integer :: id_clock_kpp

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

!  This subroutine imposes the diapycnal mass fluxes and the
!  accompanying diapycnal advection of momentum and tracers.

! Arguments:
!  (in/out)  u      = Zonal velocity (m/s)
!  (in/out)  v      = Meridional velocity (m/s)
!  (in/out)  h      = Layer thickness (m for Bouss and kg/m^2 for non-Bouss)
!  (in)      tv     = A structure containing pointers to any available
!                     thermodynamic fields; unused fields have NULL ptrs.
!  (in)      fluxes = A structure containing pointers to any possible
!                     forcing fields; unused fields have NULL ptrs.
!  (in/out)  visc   = A structure containing vertical viscosities, bottom boundary
!                     layer properies, and related fields
!  (inout)   ADp    = A structure with pointers to the various accelerations in
!                     the momentum equations, to enable the later calculation
!                     of derived diagnostics, like energy budgets
!  (inout)   CDp    = structure with pointers to terms in continuity equations
!  (in)      dt     = Time increment (seconds)
!  (in)      G      = ocean grid structure
!  (in)      CS     = control structure returned by a previous diabatic_driver_init call

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    ea,     &    ! amount of fluid entrained from the layer above within
                 ! one time step  (m for Bouss, kg/m^2 for non-Bouss)
    eb,     &    ! amount of fluid entrained from the layer below within
                 ! one time step  (m for Bouss, kg/m^2 for non-Bouss)
    Kd,     &    ! diapycnal diffusivity of layers (m^2/sec)
    h_orig, &    ! initial layer thicknesses (m for Bouss, kg/m^2 for non-Bouss)
    hold,   &    ! layer thickness before diapycnal entrainment, and later
                 ! the initial layer thicknesses (if a mixed layer is used),
                 ! (m for Bouss, kg/m^2 for non-Bouss)
    dSV_dT, &    ! The partial derivatives of specific volume with temperature
    dSV_dS, &    ! and salinity in m^3/(kg K) and m^3/(kg ppt).
    cTKE, &      ! convective TKE requirements for each layer in J/m^2.
    u_h,   &     ! zonal and meridional velocities at thickness points after
    v_h          ! entrainment (m/s)

  real, dimension(SZI_(G),SZJ_(G)) :: &
    cg1, &   ! first baroclinic gravity wave speed.
    TKE_itidal_input_test, & 
             ! OVERRIDE OF ENERGY INPUT TO INT TIDE (BDM)
    Rcv_ml   ! coordinate density of mixed layer, used for applying sponges.

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

  real :: htot(SZIB_(G))  ! The summed thickness from the bottom, in m.
  real :: b1(SZIB_(G)), d1(SZIB_(G)) ! b1, c1, and d1 are variables used by the
  real :: c1(SZIB_(G),SZK_(G))       ! tridiagonal solver.

  real :: Ent_int      ! The diffusive entrainment rate at an interface, in H.
  real :: dt_mix       ! amount of time over which to apply mixing (seconds)
  real :: Idt          ! inverse time step (1/s)

  type(p3d) :: z_ptrs(7)  ! pointers to diagnostics to be interpolated to depth
  integer :: num_z_diags  ! number of diagnostics to be interpolated to depth
  integer :: z_ids(7)     ! id numbers of diagnostics to be interpolated to depth
  logical :: showCallTree ! If true, show the call tree
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  
  integer :: ig, jg      ! (BDM)
  logical :: avg_enabled ! (BDM)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: T_f, S_f      
                         ! The temperature and salinity in C and PSU with the values in
                         ! the massless layers filled vertically by diffusion. (BDM)
  real, dimension(SZI_(G),SZJ_(G)) :: N2_bot        
                         ! The bottom squared buoyancy frequency, in s-2. (BDM)

  is   = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = G%ke
  Isq  = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nkmb = G%nk_rho_varies
  h_neglect = G%H_subroundoff ; h_neglect2 = h_neglect*h_neglect
  Kd_heat(:,:,:) = 0.0 ; Kd_salt(:,:,:) = 0.0


  if (nz == 1) return
  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("diabatic(), MOM_diabatic_driver.F90")

  ! set equivalence between the same bits of memory for these arrays
  eaml => eatr ; ebml => ebtr

  Idt = 1.0 / dt

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_diabatic_driver: "// &
         "Module must be initialized before it is used.")

  if (CS%debug) then
    call MOM_state_chksum("Start of diabatic ", u(:,:,:), v(:,:,:), h(:,:,:), G)
    call MOM_forcing_chksum("Start of diabatic", fluxes, G, haloshift=0)
  endif
  if (CS%debugConservation) call MOM_state_stats('Start of diabatic', u, v, h, tv%T, tv%S, G)

  call cpu_clock_begin(id_clock_set_diffusivity)
  call set_BBL_TKE(u, v, h, fluxes, visc, G, CS%set_diff_CSp)
  call cpu_clock_end(id_clock_set_diffusivity)

  ! Frazil formation keeps the temperature above the freezing point.
  ! make_frazil is deliberately called at both the beginning and at
  ! the end of the diabatic processes.
  if (ASSOCIATED(tv%T) .AND. ASSOCIATED(tv%frazil)) then
    if (ASSOCIATED(fluxes%p_surf_full)) then
        call make_frazil(h, tv, G, CS%diabatic_aux_CSp, fluxes%p_surf_full)
    else
        call make_frazil(h, tv, G, CS%diabatic_aux_CSp)
    endif
    if (showCallTree) call callTree_waypoint("done with 1st make_frazil (diabatic)")
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
    call geothermal(h, tv, dt, eaml, ebml, G, CS%geothermal_CSp)
    call cpu_clock_end(id_clock_geothermal)
    if (showCallTree) call callTree_waypoint("geothermal (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('geothermal', u, v, h, tv%T, tv%S, G)
  endif

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_target_grids(CS%diag)

  ! Set_opacity estimates the optical properties of the water column.
  ! It will need to be modified later to include information about the
  ! biological properties and layer thicknesses.
  if (associated(CS%optics)) &
    call set_opacity(CS%optics, fluxes, G, CS%opacity_CSp)

  if (CS%bulkmixedlayer) then
    if (CS%debug) then
      call MOM_forcing_chksum("Before mixedlayer", fluxes, G, haloshift=0)
    endif

    if (CS%ML_mix_first > 0.0) then
!  This subroutine (1)  Cools the mixed layer.
!    (2) Performs convective adjustment by mixed layer entrainment.
!    (3) Heats the mixed layer and causes it to detrain to
!        Monin-Obukhov depth or minimum mixed layer depth.
!    (4) Uses any remaining TKE to drive mixed layer entrainment.
!    (5) Possibly splits buffer layer into two isopycnal layers (when using isopycnal coordinate)
      call find_uv_at_h(u, v, h, u_h, v_h, G)

      call cpu_clock_begin(id_clock_mixedlayer)
      if (CS%ML_mix_first < 1.0) then
        ! Changes: h, tv%T, tv%S, eaml and ebml  (G is also inout???)
        call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt*CS%ML_mix_first, &
                            eaml,ebml, G, CS%bulkmixedlayer_CSp, CS%optics, &
                            CS%aggregate_FW_forcing, dt, last_call=.false.)
        if (CS%salt_reject_below_ML) &
          call insert_brine(h, tv, G, fluxes, nkmb, CS%diabatic_aux_CSp, &
                            dt*CS%ML_mix_first, CS%id_brine_lay)
      else
        ! Changes: h, tv%T, tv%S, eaml and ebml  (G is also inout???)
        call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt, eaml, ebml, &
                        G, CS%bulkmixedlayer_CSp, CS%optics, &
                        CS%aggregate_FW_forcing, dt, last_call=.true.)
      endif

!  Keep salinity from falling below a small but positive threshold.
!  This occurs when the ice model attempts to extract more salt than
!  is actually present in the ocean.
      if (ASSOCIATED(tv%S) .and. ASSOCIATED(tv%salt_deficit)) &
        call adjust_salt(h, tv, G, CS%diabatic_aux_CSp)
      call cpu_clock_end(id_clock_mixedlayer)
      if (CS%debug) then
        call MOM_state_chksum("After mixedlayer ", u, v, h, G)
        call MOM_forcing_chksum("After mixedlayer", fluxes, G, haloshift=0)
      endif
      if (showCallTree) call callTree_waypoint("done with 1st bulkmixedlayer (diabatic)")
      if (CS%debugConservation) call MOM_state_stats('1st bulkmixedlayer', u, v, h, tv%T, tv%S, G)
    endif
  endif

!  This subroutine applies diapycnal diffusion and the
!  surface buoyancy forcing.  No layer is less than an Angstrom thick at
!  the end, and total density (heat or salt) is conserved unless all
!  of the mass is in the lightest or densest layer.
!  When the bulk mixed layer is used, the flux from the buffer layer
!  is applied to the next denser layer, while all lighter layers are
!  effectively removed from the calculation.  Also, buoyancy forcing
!  is applied to the bulk mixed layer elsewhere.

! Calculate appropriately limited diapycnal mass fluxes to account
! for diapycnal diffusion and advection.
  if (CS%debug) then
    call MOM_state_chksum("before find_uv_at_h", u(:,:,:), v(:,:,:), h(:,:,:), G)
  endif
  if (CS%use_kappa_shear) then
    if ((CS%ML_mix_first > 0.0) .or. CS%use_geothermal) then
      call find_uv_at_h(u, v, h_orig, u_h, v_h, G, eaml, ebml)
      if (CS%debug) then
        call hchksum(eaml, "after find_uv_at_h eaml",G)
        call hchksum(ebml, "after find_uv_at_h ebml",G)
      endif
    else
      call find_uv_at_h(u, v, h, u_h, v_h, G)
    endif
    if (showCallTree) call callTree_waypoint("done with find_uv_at_h (diabatic)")
  endif

  if (CS%use_int_tides) then
    !   This block provides an interface for the unresolved low-mode internal
    ! tide module. It will eventually be used to provide an energy input to
    ! set_diffusivity.
    
    ! PROVIDE ENERGY DISTRIBUTION (calculate time-varying energy source)
    call set_int_tide_input(u, v, h, tv, fluxes, CS%int_tide_input, dt, G, &
                            CS%int_tide_input_CSp)
    ! CALCULATE MODAL VELOCITY
    cg1(:,:) = 0.0
    if (CS%uniform_cg) then
       ! SET TO CONSTANT VALUE TO TEST PROPAGATE CODE (BDM)
       cg1(:,:) = CS%cg_test
    else
       ! CALCULATE cg1 OR OVERRIDE WITH HARD-CODED DISTRIBUTION (BDM)
       call wave_speed(h, tv, G, cg1, full_halos=.true.)
       !do j=G%jsd,G%jed        
       !  cg1(:,j) = ((7.-1.)/14000000.)*G%geoLatBu(:,j) + (1.-((7.-1.)/14000000.)*-7000000.)
       !enddo
    endif

    if (CS%int_tide_source_test) then
      ! BUILD 2D ARRAY WITH POINT SOURCE FOR TESTING (BDM)
      TKE_itidal_input_test(:,:) = 0.0
      avg_enabled = query_averaging_enabled(CS%diag,time_end=CS%time_end)
      if (CS%time_end <= CS%time_max_source) then
        do j=G%jsc,G%jec; do i=G%isc,G%iec
         ig = G%isd_global + i - 1.0
         jg = G%jsd_global + j - 1.0
         !INPUT ARBITRARY ENERGY POINT SOURCE (BDM)
         if (ig .eq. CS%int_tide_source_x .and. jg .eq. CS%int_tide_source_y) then
             TKE_itidal_input_test(i,j) = 1.0
        endif
        enddo; enddo
      endif
      ! CALL ROUTINE USING PRESCRIBED KE FOR TESTING (BDM)
      call propagate_int_tide(h, tv, cg1, TKE_itidal_input_test, &
                            CS%int_tide_input%tideamp, CS%int_tide_input%Nb, dt, G, CS%int_tide_CSp)
    else    
      ! CALL ROUTINE USING CALCULATED KE INPUT
      call propagate_int_tide(h, tv, cg1, CS%int_tide_input%TKE_itidal_input, &
                              CS%int_tide_input%tideamp, CS%int_tide_input%Nb, dt, G, CS%int_tide_CSp)    
    endif
    if (showCallTree) call callTree_waypoint("done with propagate_int_tide (diabatic)")
  endif

  call cpu_clock_begin(id_clock_set_diffusivity)
  ! Sets: Kd, Kd_int, visc%Kd_extra_T, visc%Kd_extra_S
  ! Also changes: visc%Kd_turb, visc%TKE_turb (not clear that TKE_turb is used as input ????)
  ! And sets visc%Kv_turb
  call set_diffusivity(u, v, h, u_h, v_h, tv, fluxes, CS%optics, visc, dt, G, CS%set_diff_CSp, Kd, Kd_int)
  call cpu_clock_end(id_clock_set_diffusivity)
  if (showCallTree) call callTree_waypoint("done with set_diffusivity (diabatic)")

  if (CS%debug) then
    call MOM_state_chksum("after set_diffusivity ", u(:,:,:), v(:,:,:), h(:,:,:), G)
    call MOM_forcing_chksum("after set_diffusivity ", fluxes, G, haloshift=0)
    call MOM_thermovar_chksum("after set_diffusivity ", tv, G)
    call hchksum(Kd, "after set_diffusivity Kd",G,haloshift=0)
    call hchksum(Kd_Int, "after set_diffusivity Kd_Int",G,haloshift=0)
  endif


  if (CS%useKPP) then
    call cpu_clock_begin(id_clock_kpp)
    ! KPP needs the surface buoyancy flux but does not update state variables.
    ! We could make this call higher up to avoid a repeat unpacking of the surface fluxes. 
    ! Sets: CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux
    ! NOTE: CS%KPP_buoy_flux, CS%KPP_temp_flux, CS%KPP_salt_flux are returned as rates (i.e. stuff per second)
    ! unlike other instances where the fluxes are integrated in time over a time-step.
    call calculateBuoyancyFlux2d(G, fluxes, CS%optics, h, tv%T, tv%S, tv, &
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

    call KPP_calculate(CS%KPP_CSp, G, h, tv%T, tv%S, u, v, tv%eqn_of_state, &
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
      call MOM_state_chksum("after KPP", u(:,:,:), v(:,:,:), h(:,:,:), G)
      call MOM_forcing_chksum("after KPP", fluxes, G, haloshift=0)
      call MOM_thermovar_chksum("after KPP", tv, G)
      call hchksum(Kd, "after KPP Kd",G,haloshift=0)
      call hchksum(Kd_Int, "after KPP Kd_Int",G,haloshift=0)
    endif

  endif  ! endif for KPP

  ! Check for static instabilities and increase Kd_int where unstable
  if (CS%useConvection) call diffConvection_calculate(CS%Conv_CSp, &
         G, h, tv%T, tv%S, tv%eqn_of_state, Kd_int)

  if (CS%useKPP) then
    call cpu_clock_begin(id_clock_kpp)
    if (CS%debug) then
      call hchksum(CS%KPP_temp_flux*G%H_to_m, "before KPP_applyNLT netHeat",G,haloshift=0)
      call hchksum(CS%KPP_salt_flux*G%H_to_m, "before KPP_applyNLT netSalt",G,haloshift=0)
      call hchksum(CS%KPP_NLTheat, "before KPP_applyNLT NLTheat",G,haloshift=0)
      call hchksum(CS%KPP_NLTscalar, "before KPP_applyNLT NLTscalar",G,haloshift=0)
    endif
    ! Apply non-local transport of heat and salt
    ! Changes: tv%T, tv%S
    call KPP_NonLocalTransport_temp(CS%KPP_CSp, G, h, CS%KPP_NLTheat,   CS%KPP_temp_flux, dt, tv%T, tv%C_p)
    call KPP_NonLocalTransport_saln(CS%KPP_CSp, G, h, CS%KPP_NLTscalar, CS%KPP_salt_flux, dt, tv%S)
    call cpu_clock_end(id_clock_kpp)
    if (showCallTree) call callTree_waypoint("done with KPP_applyNonLocalTransport (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('KPP_applyNonLocalTransport', u, v, h, tv%T, tv%S, G)

    if (CS%debug) then
      call MOM_state_chksum("after KPP_applyNLT ", u(:,:,:), v(:,:,:), h(:,:,:), G)
      call MOM_forcing_chksum("after KPP_applyNLT ", fluxes, G, haloshift=0)
      call MOM_thermovar_chksum("after KPP_applyNLT ", tv, G)
    endif
  endif

  ! if using matching within the KPP scheme, then this step needs to provide
  ! a diffusivity and happen before KPP.  But generally in MOM, we do not match
  ! KPP boundary layer to interior, so this diffusivity can be computed when convenient.
  if (associated(visc%Kd_extra_T) .and. associated(visc%Kd_extra_S) .and. associated(tv%T)) then
    call cpu_clock_begin(id_clock_differential_diff)
    ! Changes: tv%T, tv%S
    call differential_diffuse_T_S(h, tv, visc, dt, G)
    call cpu_clock_end(id_clock_differential_diff)
    if (showCallTree) call callTree_waypoint("done with differential_diffuse_T_S (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('differential_diffuse_T_S', u, v, h, tv%T, tv%S, G)
  endif

  ! This block sets ea, eb from Kd or Kd_int.
  !   If using the ALE algorithm, set ea=eb=Kd_int on interfaces for
  ! use in the tri-diagonal solver.
  !   Otherwise, call entrainment_diffusive() which sets ea and eb
  ! based on KD and target densities (ie. does remapping as well).
  if (CS%useALEalgorithm) then
    do j=js,je ; do i=is,ie
      ea(i,j,1) = 0.
    enddo ; enddo
!$OMP parallel do default(none) shared(is,ie,js,je,nz,h_neglect,h,ea,G,dt,Kd_int,eb) &
!$OMP                          private(hval)
    do k=2,nz ; do j=js,je ; do i=is,ie
      hval=1.0/(h_neglect + 0.5*(h(i,j,k-1) + h(i,j,k)))
      ea(i,j,k) = (G%m_to_H**2) * dt * hval * Kd_int(i,j,k)
      eb(i,j,k-1) = ea(i,j,k)
    enddo ; enddo ; enddo
    do j=js,je ; do i=is,ie
      eb(i,j,nz) = 0.
    enddo ; enddo
    if (showCallTree) call callTree_waypoint("done setting ea,eb from Kd_int (diabatic)")
  else ! .not. CS%useALEalgorithm
    ! If not useing ALE, then calculate layer entrainments/detrainments from
    ! diffusivities and differences between layer and target densities
    call cpu_clock_begin(id_clock_entrain)
    ! Sets: ea, eb. Changes: kb
    call Entrainment_diffusive(u, v, h, tv, fluxes, dt, G, CS%entrain_diffusive_CSp, &
                               ea, eb, kb, Kd_Lay=Kd, Kd_int=Kd_int)
    call cpu_clock_end(id_clock_entrain)
    if (showCallTree) call callTree_waypoint("done with Entrainment_diffusive (diabatic)")
  endif ! (CS%useALEalgorithm)

  if (CS%debug) then
    call MOM_forcing_chksum("after calc_entrain ", fluxes, G, haloshift=0)
    call MOM_thermovar_chksum("after calc_entrain ", tv, G)
    call MOM_state_chksum("after calc_entrain ", u(:,:,:), v(:,:,:), h(:,:,:), G)
    call hchksum(G%H_to_m*ea, "after calc_entrain ea",G,haloshift=0)
    call hchksum(G%H_to_m*eb, "after calc_entrain eb",G,haloshift=0)
  endif

  if (CS%id_Kd_ePBL > -1) Kd_ePBL(:,:,:) = 0.0

  ! Apply forcing when using the ALE algorithm
  if (CS%useALEalgorithm) then
    call cpu_clock_begin(id_clock_remap)
    ! Changes: ea(:,:,1), h, tv%T and tv%S.

    if (CS%use_energetic_PBL) then
      call applyBoundaryFluxesInOut(CS%diabatic_aux_CSp, G, dt, fluxes, CS%optics, &
                          ea, h, tv, CS%aggregate_FW_forcing, cTKE, dSV_dT, dSV_dS)
      if (CS%debug) then
        call hchksum(ea, "after applyBoundaryFluxes ea",G,haloshift=0)
        call hchksum(eb, "after applyBoundaryFluxes eb",G,haloshift=0)
        call hchksum(cTKE, "after applyBoundaryFluxes cTKE",G,haloshift=0)
        call hchksum(dSV_dT, "after applyBoundaryFluxes dSV_dT",G,haloshift=0)
        call hchksum(dSV_dS, "after applyBoundaryFluxes dSV_dS",G,haloshift=0)
      endif
      call find_uv_at_h(u, v, h, u_h, v_h, G)
      call energetic_PBL(h, u_h, v_h, tv, fluxes, dt, Kd_ePBL, G, &
                         CS%energetic_PBL_CSp, dSV_dT, dSV_dS, cTKE)
      ! Augment the diffusivities due to those diagnosed in energetic_PBL.
      do K=2,nz ; do j=js,je ; do i=is,ie
        Ent_int = Kd_ePBL(i,j,K) * (G%m_to_H**2 * dt) / &
                  (0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect)
        eb(i,j,k-1) = eb(i,j,k-1) + Ent_int
        ea(i,j,k) = ea(i,j,k) + Ent_int
        visc%Kv_turb(i,j,K) = visc%Kv_turb(i,j,K) + Kd_ePBL(i,j,K)
        Kd_int(i,j,K)  = Kd_int(i,j,K)  + Kd_ePBL(i,j,K)
        Kd_heat(i,j,K) = Kd_heat(i,j,K) + Kd_ePBL(i,j,K)
        Kd_salt(i,j,K) = Kd_salt(i,j,K) + Kd_ePBL(i,j,K)
      enddo ; enddo ; enddo
    if (CS%debug) then
      call hchksum(ea, "after ePBL ea",G,haloshift=0)
      call hchksum(eb, "after ePBL eb",G,haloshift=0)
      call hchksum(Kd_ePBL, "after ePBL Kd_ePBL",G,haloshift=0)
    endif
    else
      call applyBoundaryFluxesInOut(CS%diabatic_aux_CSp, G, dt, fluxes, CS%optics, &
                                    ea, h, tv, CS%aggregate_FW_forcing)
    endif

    call cpu_clock_end(id_clock_remap)
    if (CS%debug) then
      call MOM_forcing_chksum("after applyBoundaryFluxes ", fluxes, G, haloshift=0)
      call MOM_thermovar_chksum("after applyBoundaryFluxes ", tv, G)
      call MOM_state_chksum("after applyBoundaryFluxes ", u(:,:,:), v(:,:,:), h(:,:,:), G)
    endif
    if (showCallTree) call callTree_waypoint("done with applyBoundaryFluxes (diabatic)")
    if (CS%debugConservation)  call MOM_state_stats('applyBoundaryFluxes', u, v, h, tv%T, tv%S, G)
  endif

  ! Update h according to divergence of the difference between
  ! ea and eb. We keep a record of the original h in hold.
  ! In the following, the checks for negative values are to guard
  ! against against unforeseen instances where the entrainment
  ! drives a layer to negative thickness.  This will never happen if
  ! enough iterations are permitted in Calculate_Entrainment, and
  ! even if too few iterations are allowed, it is still guarded
  ! against.  In other words the checks are probably unnecessary.

!$OMP parallel do default(none) shared(is,ie,js,je,nz,hold,h,eb,ea,G)
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
    call MOM_forcing_chksum("after negative check ", fluxes, G, haloshift=0)
    call MOM_thermovar_chksum("after negative check ", tv, G)
  endif
  if (showCallTree) call callTree_waypoint("done with h=ea-eb (diabatic)")
  if (CS%debugConservation) call MOM_state_stats('h=ea-eb', u, v, hold, tv%T, tv%S, G)

  ! Here, T and S are updated according to ea and eb.
  ! If using the bulk mixed layer, T and S are also updated
  ! by surface fluxes (in fluxes%*).
  if (CS%bulkmixedlayer) then
    if (ASSOCIATED(tv%T)) then
      call cpu_clock_begin(id_clock_tridiag)
      ! Temperature and salinity (as state variables) are treated slightly
      ! differently from other tracers to insure that massless layers that
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
        call triDiagTS(G, is, ie, js, je, hold, ea, eb, tv%T, tv%S)
      endif ! massless_match_targets
      call cpu_clock_end(id_clock_tridiag)
      if (CS%debugConservation) call MOM_state_stats('BML tridiag', u, v, h, tv%T, tv%S, G)
    endif ! end of ASSOCIATED(T)

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
        call hchksum(G%H_to_m*ea, "after ea = ea + eaml",G,haloshift=0)
        call hchksum(G%H_to_m*eb, "after eb = eb + ebml",G,haloshift=0)
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

      call find_uv_at_h(u, v, hold, u_h, v_h, G, ea, eb)
      if (CS%debug) call MOM_state_chksum("find_uv_at_h1 ", u, v, h, G)

      dt_mix = min(dt,dt*(1.0 - CS%ML_mix_first))
      call cpu_clock_begin(id_clock_mixedlayer)
      ! Changes: h, tv%T, tv%S, ea and eb  (G is also inout???)
      call bulkmixedlayer(h, u_h, v_h, tv, fluxes, dt_mix, ea, eb, &
                      G, CS%bulkmixedlayer_CSp, CS%optics, &
                      CS%aggregate_FW_forcing, dt, last_call=.true.)

      if (CS%salt_reject_below_ML) &
        call insert_brine(h, tv, G, fluxes, nkmb, CS%diabatic_aux_CSp, dt_mix, &
                          CS%id_brine_lay)
!  Keep salinity from falling below a small but positive threshold.
!  This occurs when the ice model attempts to extract more salt than
!  is actually present in the ocean.
      if (ASSOCIATED(tv%S) .and. ASSOCIATED(tv%salt_deficit)) &
        call adjust_salt(h, tv, G, CS%diabatic_aux_CSp)

      call cpu_clock_end(id_clock_mixedlayer)
      if (showCallTree) call callTree_waypoint("done with 2nd bulkmixedlayer (diabatic)")
      if (CS%debugConservation) call MOM_state_stats('2nd bulkmixedlayer', u, v, h, tv%T, tv%S, G)
    endif

  else                                             ! Not BULKMIXEDLAYER.

! Calculate the change in temperature & salinity due to entrainment.
    if (ASSOCIATED(tv%T)) then
      if (CS%debug) then
        call hchksum(G%H_to_m*ea, "before triDiagTS ea ",G,haloshift=0)
        call hchksum(G%H_to_m*eb, "before triDiagTS eb ",G,haloshift=0)
      endif
      call cpu_clock_begin(id_clock_tridiag)
      ! Changes: T, S
      call triDiagTS(G, is, ie, js, je, hold, ea, eb, tv%T, tv%S)
      call cpu_clock_end(id_clock_tridiag)
      if (showCallTree) call callTree_waypoint("done with triDiagTS (diabatic)")
      if (CS%debugConservation) call MOM_state_stats('triDiagTS', u, v, h, tv%T, tv%S, G)
    endif

  endif                                          ! end BULKMIXEDLAYER
  if (CS%debug) then
    call MOM_state_chksum("after mixed layer ", u, v, h, G)
    call MOM_thermovar_chksum("after mixed layer ", tv, G)
    call hchksum(ea, "after mixed layer ea",G)
    call hchksum(eb, "after mixed layer eb",G)
  endif

  if (.not. CS%useALEalgorithm) then
    call cpu_clock_begin(id_clock_remap)
    call regularize_layers(h, tv, dt, ea, eb, G, CS%regularize_layers_CSp)
    call cpu_clock_end(id_clock_remap)
    if (showCallTree) call callTree_waypoint("done with regularize_layers (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('regularize_layers', u, v, h, tv%T, tv%S, G)
  endif

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_target_grids(CS%diag)

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

  call cpu_clock_begin(id_clock_tracers)
  if (CS%mix_boundary_tracers) then
    Tr_ea_BBL = sqrt(dt*CS%Kd_BBL_tr)
!$OMP parallel do default(none) shared(is,ie,js,je,ebtr,nz,G,h,dt,CS,h_neglect,     &
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
!$OMP parallel do default(none) shared(nz,is,ie,js,je,visc,dt,G,h,hold,h_neglect, &
!$OMP                                  ebtr,eb,eatr,ea )                          &
!$OMP                          private(add_ent)
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
!$OMP parallel do default(none) shared(js,je,p_ref_cv,Rcv_ml,is,ie,tv)
      do j=js,je
        call calculate_density(tv%T(:,j,1), tv%S(:,j,1), p_ref_cv, Rcv_ml(:,j), &
                               is, ie-is+1, tv%eqn_of_state)
      enddo
      call apply_sponge(h, dt, G, ea, eb, CS%sponge_CSp, Rcv_ml)
    else
      call apply_sponge(h, dt, G, ea, eb, CS%sponge_CSp)
    endif
    call cpu_clock_end(id_clock_sponge)
    if (CS%debug) then
      call MOM_state_chksum("apply_sponge ", u, v, h, G)
      call MOM_thermovar_chksum("apply_sponge ", tv, G)
    endif
  endif

!$OMP parallel default(none) shared(is,ie,js,je,nz,CDp,Idt,G,ea,eb,CS,hold)
!   Save the diapycnal mass fluxes as a diagnostic field.
  if (ASSOCIATED(CDp%diapyc_vel)) then
!$OMP do
    do j=js,je
      do K=2,nz ;  do i=is,ie
        CDp%diapyc_vel(i,j,K) = Idt * (G%H_to_m * (ea(i,j,k) - eb(i,j,k-1)))
      enddo ; enddo
      do i=is,ie
        CDp%diapyc_vel(i,j,1) = 0.0
        CDp%diapyc_vel(i,j,nz+1) = 0.0
      enddo
    enddo
  endif

!   For momentum, it is only the net flux that homogenizes within
!  the mixed layer.  Vertical viscosity that is proportional to the
!  mixed layer turbulence is applied elsewhere.
  if (CS%bulkmixedlayer) then
    if (CS%debug) then
      call hchksum(ea, "before net flux rearrangement ea",G)
      call hchksum(eb, "before net flux rearrangement eb",G)
    endif
!$OMP do
    do k=2,G%nkml ; do j=js,je ; do i=is,ie
      if (ea(i,j,k) >= eb(i,j,k-1)) then
        ea(i,j,k) = ea(i,j,k) - eb(i,j,k-1)
        eb(i,j,k-1) = 0.0
      else
        eb(i,j,k-1) = eb(i,j,k-1) - ea(i,j,k)
        ea(i,j,k) = 0.0
      endif
    enddo ; enddo ; enddo
    if (CS%debug) then
      call hchksum(ea, "after net flux rearrangement ea",G)
      call hchksum(eb, "after net flux rearrangement eb",G)
    endif
  endif


!   This initializes the halo regions of ea, eb, and hold to default values.
!$OMP do
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

!  Use a tridiagonal solver to determine the effect of the diapycnal
!  advection on the velocity field.   It is assumed that water leaves
!  or enters the ocean with the surface velocity.

  if (CS%debug) then
    call MOM_state_chksum("before u/v tridiag ", u, v, h, G)
    call hchksum(ea, "before u/v tridiag ea",G)
    call hchksum(eb, "before u/v tridiag eb",G)
    call hchksum(hold, "before u/v tridiag hold",G)
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
    call MOM_state_chksum("aft 1st loop tridiag ", u, v, h, G)
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
    call MOM_state_chksum("after u/v tridiag ", u, v, h, G)
  endif

!   Frazil formation keeps the temperature above the freezing point.
! make_frazil is deliberately called at both the beginning and at
! the end of the diabatic processes.
  if (ASSOCIATED(tv%T) .AND. ASSOCIATED(tv%frazil)) then
    if (ASSOCIATED(fluxes%p_surf_full)) then
       call make_frazil(h, tv, G, CS%diabatic_aux_CSp, fluxes%p_surf_full)
    else
       call make_frazil(h, tv, G, CS%diabatic_aux_CSp)
    endif

    if (showCallTree) call callTree_waypoint("done with 2nd make_frazil (diabatic)")
    if (CS%debugConservation) call MOM_state_stats('2nd make_frazil', u, v, h, tv%T, tv%S, G)
  endif

  ! Diagnose the diapycnal diffusivities and other quantities related to diapycnal mixing.
  if (CS%id_Kd_interface > 0) call post_data(CS%id_Kd_interface, Kd_int, CS%diag)
  if (CS%id_Kd_heat > 0)      call post_data(CS%id_Kd_heat,     Kd_heat, CS%diag)
  if (CS%id_Kd_salt > 0)      call post_data(CS%id_Kd_salt,     Kd_salt, CS%diag)
  if (CS%id_Kd_ePBL > 0)      call post_data(CS%id_Kd_ePBL,     Kd_ePBL, CS%diag)

  if (CS%id_ea > 0) call post_data(CS%id_ea, ea, CS%diag)
  if (CS%id_eb > 0) call post_data(CS%id_eb, eb, CS%diag)
  if (CS%id_dudt_dia > 0) call post_data(CS%id_dudt_dia, ADp%du_dt_dia, CS%diag)
  if (CS%id_dvdt_dia > 0) call post_data(CS%id_dvdt_dia, ADp%dv_dt_dia, CS%diag)
  if (CS%id_wd > 0) call post_data(CS%id_wd, CDp%diapyc_vel, CS%diag)
  if (CS%id_MLD_003 > 0 .or. CS%id_subMLN2 > 0 .or. CS%id_mlotstsq > 0) then
    call diagnoseMLDbyDensityDifference(CS%id_MLD_003, h, tv, 0.03, G, CS%diag, &
                                        id_N2subML=CS%id_subMLN2, id_MLDsq=CS%id_mlotstsq)
  endif
  if (CS%id_MLD_0125 > 0) then
    call diagnoseMLDbyDensityDifference(CS%id_MLD_0125, h, tv, 0.125, G, CS%diag)
  endif
  if (CS%id_MLD_user > 0) then
    call diagnoseMLDbyDensityDifference(CS%id_MLD_user, h, tv, CS%MLDdensityDifference, G, CS%diag)
  endif

  if (CS%id_Tdif > 0) call post_data(CS%id_Tdif, Tdif_flx, CS%diag)
  if (CS%id_Tadv > 0) call post_data(CS%id_Tadv, Tadv_flx, CS%diag)
  if (CS%id_Sdif > 0) call post_data(CS%id_Sdif, Sdif_flx, CS%diag)
  if (CS%id_Sadv > 0) call post_data(CS%id_Sadv, Sadv_flx, CS%diag)
  if (CS%id_cg1  > 0) call post_data(CS%id_cg1,  cg1,      CS%diag) !(BDM)
  
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

  if (CS%debugConservation) call MOM_state_stats('leaving diabatic', u, v, h, tv%T, tv%S, G)
  if (showCallTree) call callTree_leave("diabatic()")
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


subroutine adiabatic_driver_init(Time, G, param_file, diag, CS, &
                                tracer_flow_CSp, diag_to_Z_CSp)
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(diabatic_CS),       pointer       :: CS
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp
  type(diag_to_Z_CS),      pointer       :: diag_to_Z_CSp

! Arguments:
!  (in)      Time       = current model time
!  (in)      G          = ocean grid structure
!  (in)      param_file = structure indicating the open file to parse for parameter values
!  (in)      diag       = structure used to regulate diagnostic output
!  (in/out)  CS         = pointer set to point to the control structure for this module
!  (in)      tracer_flow_CSp = A pointer to the control structure of the tracer
!                              flow control module.
!  (in)      diag_to_Z_CSp   = A pointer to the Z-diagnostics control structure.

! This is a simplified version of diabatic_driver_init that will allow
! tracer column functions to be called without allowing any of the diabatic
! processes to be used.

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


subroutine diabatic_driver_init(Time, G, param_file, useALEalgorithm, diag, &
                                ADp, CDp, CS, tracer_flow_CSp, sponge_CSp, diag_to_Z_CSp)
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(inout) :: G
  type(param_file_type),   intent(in)    :: param_file
  logical,                 intent(in)    :: useALEalgorithm
  type(diag_ctrl), target, intent(inout) :: diag
  type(accel_diag_ptrs),   intent(inout) :: ADp
  type(cont_diag_ptrs),    intent(inout) :: CDp
  type(diabatic_CS),       pointer       :: CS
  type(tracer_flow_control_CS), pointer  :: tracer_flow_CSp
  type(sponge_CS),         pointer       :: sponge_CSp
  type(diag_to_Z_CS),      pointer       :: diag_to_Z_CSp

! Arguments:
!  (in)      Time            = current model time
!  (in)      G               = ocean grid structure
!  (in)      param_file      = structure indicating the open file to parse for parameter values
!  (in)      useALEalgorithm = logical determining whether to use ALE vertical remapping
!  (in)      diag            = structure used to regulate diagnostic output
!  (inout)   ADp             = structure with pointers to accelerations in momentum equations,
!                              to enable later calculation of diagnostics, like energy budgets
!  (inout)   CDp             = structure with pointers to terms in continuity equations
!  (in/out)  CS              = pointer set to point to the control structure for this module
!  (in)      tracer_flow_CSp = pointer to control structure of tracer flow control module
!  (in)      sponge_CSp      = pointer to the sponge module control structure
!  (in)      diag_to_Z_CSp   = A pointer to the Z-diagnostics control structure.

  real :: Kd
  logical :: use_temperature, differentialDiffusion
  type(vardesc) :: vd

! This "include" declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod  = "MOM_diabatic_driver" ! This module's name.
  character(len=48)  :: thickness_units
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz, nbands
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
  if (associated(diag_to_Z_CSp))   CS%diag_to_Z_CSp   => diag_to_Z_CSp

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
  call get_param(param_file, mod, "ENERGETICS_SFC_PBL", CS%use_energetic_PBL, &
                 "If true, use an implied energetics planetary boundary \n"//&
                 "layer scheme to determine the diffusivity and viscosity \n"//&
                 "in the surface boundary layer.", default=.false.)
  call get_param(param_file, mod, "DOUBLE_DIFFUSION", differentialDiffusion, &
                 "If true, apply parameterization of double-diffusion.", &
                 default=.false. )
  CS%use_kappa_shear = kappa_shear_is_used(param_file)
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
  if (CS%use_int_tides) then
    ! GET LOCATION AND DURATION OF ENERGY POINT SOURCE FOR TESTING (BDM)
    call get_param(param_file, mod, "INTERNAL_TIDE_SOURCE_TEST", CS%int_tide_source_test, &
                 "If true, apply an arbitrary generation site for internal tide testing", &
                 default=.false.)
    call get_param(param_file, mod, "INTERNAL_TIDE_SOURCE_X", CS%int_tide_source_x, &
                 "X Location of generation site for internal tide", default=1.)
    call get_param(param_file, mod, "INTERNAL_TIDE_SOURCE_Y", CS%int_tide_source_y, &
                 "Y Location of generation site for internal tide", default=1.)
    call get_param(param_file, mod, "INTERNAL_TIDE_SOURCE_TLEN_DAYS", CS%tlen_days, &
                 "Time interval from start of experiment for adding wave source", &
                 units="days", default=0)
    CS%time_max_source = increment_time(Time,0,days=CS%tlen_days)
    ! GET UNIFORM MODE VELOCITY FOR TESTING (BDM)  
    call get_param(param_file, mod, "UNIFORM_CG", CS%uniform_cg, &
                 "If true, set cg = cg_test everywhere for test case", default=.false.)
    call get_param(param_file, mod, "CG_TEST", CS%cg_test, &
                 "Uniform group velocity of internal tide for test case", default=1.)    
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
  CS%id_cg1 = register_diag_field('ocean_model','cn1', diag%axesT1, &
                 Time, 'First baroclinic mode (eigen) speed', 'm s-1') ! (BDM)

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

  if (associated(diag_to_Z_CSp)) then
    vd = vardesc("Kd","Diapycnal diffusivity at interfaces, interpolated to z",&
                 'h','z','s',"meter2 second-1")
    CS%id_Kd_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = vardesc("Tflx_dia_dif","Diffusive diapycnal temperature flux across interfaces, interpolated to z",&
                 'h','z','s',"degC meter second-1")
    CS%id_Tdif_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = vardesc("Tflx_dia_adv","Advective diapycnal temperature flux across interfaces, interpolated to z",&
                 'h','z','s',"degC meter second-1")
    CS%id_Tadv_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = vardesc("Sflx_dia_dif","Diffusive diapycnal salinity flux across interfaces, interpolated to z",&
                 'h','z','s',"PSU meter second-1")
    CS%id_Sdif_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
    vd = vardesc("Sflx_dia_adv","Advective diapycnal salinity flux across interfaces, interpolated to z",&
                 'h','z','s',"PSU meter second-1")
    CS%id_Sadv_z = register_Zint_diag(vd, CS%diag_to_Z_CSp, Time)
  endif

  if (CS%id_dudt_dia > 0) call safe_alloc_ptr(ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
  if (CS%id_dvdt_dia > 0) call safe_alloc_ptr(ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
  if (CS%id_wd > 0)       call safe_alloc_ptr(CDp%diapyc_vel,isd,ied,jsd,jed,nz+1)
  
  if (CS%use_int_tides) then !BDM
    call int_tide_input_init(Time, G, param_file, diag, CS%int_tide_input_CSp, &
                             CS%int_tide_input)
    call internal_tides_init(Time, G, param_file, diag, CS%int_tide_CSp)
  endif

  call set_diffusivity_init(Time, G, param_file, diag, CS%set_diff_CSp, diag_to_Z_CSp, CS%int_tide_CSp)
  CS%id_Kd_interface = register_diag_field('ocean_model', 'Kd_interface', diag%axesTi, Time, &
      'Total diapycnal diffusivity at interfaces', 'meter2 second-1')
  CS%id_Kd_ePBL = register_diag_field('ocean_model', 'Kd_ePBL', diag%axesTi, Time, &
      'ePBL diapycnal diffusivity at interfaces', 'meter2 second-1')

  CS%id_Kd_heat = register_diag_field('ocean_model', 'Kd_heat', diag%axesTi, Time, &
      'Total diapycnal diffusivity for heat at interfaces', 'meter2 second-1',     &
       cmor_field_name='difvho', cmor_units='m2 s-1',                              &
       cmor_standard_name='ocean_vertical_heat_diffusivity',                       &
       cmor_long_name='Net diapycnal diffusivity for ocean heat')

  CS%id_Kd_salt = register_diag_field('ocean_model', 'Kd_salt', diag%axesTi, Time, &
      'Total diapycnal diffusivity for salt at interfaces', 'meter2 second-1',     &
       cmor_field_name='difvso', cmor_units='m2 s-1',                              &
       cmor_standard_name='ocean_vertical_salt_diffusivity',                       &
       cmor_long_name='Net diapycnal diffusivity for ocean salt')

  ! CS%useKPP is set to True if KPP-scheme is to be used, False otherwise.
  ! KPP_init() allocated CS%KPP_Csp and also sets CS%KPPisPassive
  CS%useKPP = KPP_init(param_file, G, diag, Time, CS%KPP_CSp, passive=CS%KPPisPassive)
  if (CS%useKPP) then
    allocate( CS%KPP_NLTheat(isd:ied,jsd:jed,nz+1) ) ; CS%KPP_NLTheat(:,:,:) = 0.
    allocate( CS%KPP_NLTscalar(isd:ied,jsd:jed,nz+1) ) ; CS%KPP_NLTscalar(:,:,:) = 0.
    allocate( CS%KPP_buoy_flux(isd:ied,jsd:jed,nz+1) ) ; CS%KPP_buoy_flux(:,:,:) = 0.
    allocate( CS%KPP_temp_flux(isd:ied,jsd:jed) ) ; CS%KPP_temp_flux(:,:) = 0.
    allocate( CS%KPP_salt_flux(isd:ied,jsd:jed) ) ; CS%KPP_salt_flux(:,:) = 0.
  endif

  call get_param(param_file, mod, "SALT_REJECT_BELOW_ML", CS%salt_reject_below_ML, &
                 "If true, place salt from brine rejection below the mixed layer,\n"// &
                 "into the first non-vanished layer for which the column remains stable", &
                 default=.false.)

  if (CS%salt_reject_below_ML) then
     CS%id_brine_lay = register_diag_field('ocean_model','brine_layer',diag%axesT1,Time, &
      'Brine insertion layer','none')
  endif

  ! CS%useConvection is set to True IF convection will be used, otherwise False.
  ! CS%Conv_CSp is allocated by diffConvection_init()
  CS%useConvection = diffConvection_init(param_file, G, diag, Time, CS%Conv_CSp)

  call entrain_diffusive_init(Time, G, param_file, diag, CS%entrain_diffusive_CSp)
  if (CS%use_geothermal) &
    call geothermal_init(Time, G, param_file, diag, CS%geothermal_CSp)

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

  call diabatic_aux_init(Time, G, param_file, diag, CS%diabatic_aux_CSp, &
                         CS%use_energetic_PBL)

  if (CS%bulkmixedlayer) &
    call bulkmixedlayer_init(Time, G, param_file, diag, CS%bulkmixedlayer_CSp)
  if (CS%use_energetic_PBL) &
    call energetic_PBL_init(Time, G, param_file, diag, CS%energetic_PBL_CSp)

  call regularize_layers_init(Time, G, param_file, diag, CS%regularize_layers_CSp)

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


subroutine diabatic_driver_end(CS)
  type(diabatic_CS), pointer :: CS

  if (.not.associated(CS)) return

  call diabatic_aux_end(CS%diabatic_aux_CSp)

  call entrain_diffusive_end(CS%entrain_diffusive_CSp)
  call set_diffusivity_end(CS%set_diff_CSp)
  if (CS%useKPP) then
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
    
  !if (CS%use_int_tides) &
  !  call internal_tides_end(CS%int_tide_CSp) !BDM  
  !  
  !!--------------------check----------------------------------------------
  !if (is_root_pe()) then
  !  print *,'diabatic_driver_end: internal_tides_end called!' !BDM
  !endif
  !!-----------------------------------------------------------------------

  if (associated(CS%optics)) then
    call opacity_end(CS%opacity_CSp, CS%optics)
    deallocate(CS%optics)
  endif
  if (associated(CS)) deallocate(CS)

end subroutine diabatic_driver_end


end module MOM_diabatic_driver
