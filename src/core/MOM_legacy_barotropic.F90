module MOM_legacy_barotropic
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
!*  By Robert Hallberg, April 1994 - January 2007                      *
!*                                                                     *
!*    This program contains the subroutines that time steps the        *
!*  linearized barotropic equations.  btstep is used to actually       *
!*  time step the barotropic equations, and contains most of the       *
!*  substance of this module.                                          *
!*                                                                     *
!*    btstep uses a forwards-backwards based scheme to time step       *
!*  the barotropic equations, returning the layers' accelerations due  *
!*  to the barotropic changes in the ocean state, the final free       *
!*  surface height (or column mass), and the volume (or mass) fluxes   *
!*  summed through the layers and averaged over the baroclinic time    *
!*  step.  As input, btstep takes the initial 3-D velocities, the      *
!*  inital free surface height, the 3-D accelerations of the layers,   *
!*  and the external forcing.  Everything in btstep is cast in terms   *
!*  of anomalies, so if everything is in balance, there is explicitly  *
!*  no acceleration due to btstep.                                     *
!*                                                                     *
!*    The spatial discretization of the continuity equation is second  *
!*  order accurate.  A flux conservative form is used to guarantee     *
!*  global conservation of volume.  The spatial discretization of the  *
!*  momentum equation is second order accurate.  The Coriolis force    *
!*  is written in a form which does not contribute to the energy       *
!*  tendency and which conserves linearized potential vorticity, f/D.  *
!*  These terms are exactly removed from the baroclinic momentum       *
!*  equations, so the linearization of vorticity advection will not    *
!*  degrade the overall solution.                                      *
!*                                                                     *
!*    btcalc calculates the fractional thickness of each layer at the  *
!*  velocity points, for later use in calculating the barotropic       *
!*  velocities and the averaged accelerations.  Harmonic mean          *
!*  thicknesses (i.e. 2*h_L*h_R/(h_L + h_R)) are used to avoid overly  *
!*  strong weighting of overly thin layers.  This may later be relaxed *
!*  to use thicknesses determined from the continuity equations.       *
!*                                                                     *
!*    bt_mass_source determines the real mass sources for the          *
!*  barotropic solver, along with the corrective pseudo-fluxes that    *
!*  keep the barotropic and baroclinic estimates of the free surface   *
!*  height close to each other.  Given the layer thicknesses and the   *
!*  free surface height that correspond to each other, it calculates   *
!*  a corrective mass source that is added to the barotropic continuity*
!*  equation, and optionally adjusts a slowly varying correction rate. *
!*  Newer algorithmic changes have deemphasized the need for this, but *
!*  it is still here to add net water sources to the barotropic solver.*
!*                                                                     *
!*    barotropic_init allocates and initializes any barotropic arrays  *
!*  that have not been read from a restart file, reads parameters from *
!*  the inputfile, and sets up diagnostic fields.                      *
!*                                                                     *
!*    barotropic_end deallocates anything allocated in barotropic_init *
!*  or register_barotropic_restarts.                                   *
!*                                                                     *
!*    register_barotropic_restarts is used to indicate any fields that *
!*  are private to the barotropic solver that need to be included in   *
!*  the restart files, and to ensure that they are read.               *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v_in, vbt, accel_layer_v, vbtav          *
!*    j    x ^ x ^ x   At >:  u_in, ubt, accel_layer_u, ubtav, amer    *
!*    j    > o > o >   At o:  eta, h, bathyT, pbce                     *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1                                                  *
!*           i  i+1                                                    *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_debugging, only : hchksum, uchksum, vchksum
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator, only : safe_alloc_ptr, diag_ctrl, enable_averaging
use MOM_domains, only : pass_var, pass_vector, min_across_PEs, clone_MOM_domain
use MOM_domains, only : pass_var_start, pass_var_complete
use MOM_domains, only : pass_vector_start, pass_vector_complete
use MOM_domains, only : To_All, Scalar_Pair, AGRID, CORNER, MOM_domain_type
use MOM_domains, only : pass_var_start, pass_var_complete
use MOM_domains, only : pass_vector_start, pass_vector_complete
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : vardesc, var_desc
use MOM_open_boundary, only : ocean_OBC_type, OBC_SIMPLE, OBC_NONE, OBC_FLATHER
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W
use MOM_open_boundary, only : OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_restart, only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_tidal_forcing, only : tidal_forcing_sensitivity, tidal_forcing_CS
use MOM_time_manager, only : time_type, set_time, operator(+), operator(-)
use MOM_variables, only : BT_cont_type, alloc_bt_cont_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>
#ifdef STATIC_MEMORY_
#  ifndef BTHALO_
#    define BTHALO_ 0
#  endif
#  define WHALOI_ MAX(BTHALO_-NIHALO_,0)
#  define WHALOJ_ MAX(BTHALO_-NJHALO_,0)
#  define NIMEMW_   1-WHALOI_:NIMEM_+WHALOI_
#  define NJMEMW_   1-WHALOJ_:NJMEM_+WHALOJ_
#  define NIMEMBW_  -WHALOI_:NIMEM_+WHALOI_
#  define NJMEMBW_  -WHALOJ_:NJMEM_+WHALOJ_
#  define SZIW_(G)  NIMEMW_
#  define SZJW_(G)  NJMEMW_
#  define SZIBW_(G) NIMEMBW_
#  define SZJBW_(G) NJMEMBW_
#else
#  define NIMEMW_   :
#  define NJMEMW_   :
#  define NIMEMBW_  :
#  define NJMEMBW_  :
#  define SZIW_(G)  G%isdw:G%iedw
#  define SZJW_(G)  G%jsdw:G%jedw
#  define SZIBW_(G) G%isdw-1:G%iedw
#  define SZJBW_(G) G%jsdw-1:G%jedw
#endif

public legacy_btcalc, legacy_bt_mass_source, legacy_btstep, legacy_barotropic_init
public legacy_barotropic_end, register_legacy_barotropic_restarts, legacy_set_dtbt

type, public :: legacy_barotropic_CS ; private
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: frhatu
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: frhatv
      ! frhatu and frhatv are the fraction of the total column thickness
      ! interpolated to u or v grid points in each layer, nondimensional.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    IDatu, &        ! Inverse of the basin depth at u grid points, in m-1.
    uhbt_IC, &      ! The barotropic solver's estimate of the zonal
                    ! transport as the initla condition for the next call
                    ! to btstep, in H m2 s-1.
    ubt_IC, &       ! The barotropic solver's estimate of the zonal velocity
                    ! that will be the initial condition for the next call
                    ! to btstep, in m s-1.
    ubtav           ! The barotropic zonal velocity averaged over the
                    ! baroclinic time step, m s-1.
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    IDatv, &        ! Inverse of the basin depth at v grid points, in m-1.
    vhbt_IC, &      ! The barotropic solver's estimate of the zonal
                    ! transport as the initla condition for the next call
                    ! to btstep, in H m2 s-1.
    vbt_IC, &       ! The barotropic solver's estimate of the zonal velocity
                    ! that will be the initial condition for the next call
                    ! to btstep, in m s-1.
    vbtav           ! The barotropic meridional velocity averaged over the
                    ! baroclinic time step, m s-1.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    eta_source, &   ! The net mass source to be applied within the
                    ! barotropic solver, in H s-1.
    eta_cor, &      ! The difference between the free surface height from
                    ! the barotropic calculation and the sum of the layer
                    ! thicknesses. This difference is imposed as a forcing
                    ! term in the barotropic calculation over a baroclinic
                    ! timestep, in H (m or kg m-2).
    eta_cor_bound   ! A limit on the rate at which eta_cor can be applied
                    ! while avoiding instability, in units of H s-1. This
                    ! is only used if CS%bound_BT_corr is true.
  real ALLOCABLE_, dimension(NIMEMW_,NJMEMW_) :: &
    ua_polarity, &  ! Test vector components for checking grid polarity.
    va_polarity, &  ! Test vector components for checking grid polarity.
    bathyT          !   A copy of bathyT (ocean bottom depth) with wide halos.
  real ALLOCABLE_, dimension(NIMEMW_,NJMEMW_) :: &
    IareaT          !   This is a copy of G%IareaT with wide halos, but will
                    ! still utilize the macro IareaT when referenced, m-2.
  real ALLOCABLE_, dimension(NIMEMBW_,NJMEMW_) :: &
    Datu_res, &     ! A nondimensional factor by which the zonal face areas
                    ! are to be rescaled to account for the effective face
                    ! areas of the layers, if rescale_D_bt is true.
                    ! Datu_res is set in btcalc.
    D_u_Cor, &      !   A simply averaged depth at u points, in m.
    dy_Cu, &         !   A copy of G%dy_Cu with wide halos, in m.
    IdxCu            !   A copy of G%IdxCu with wide halos, in m-1.
  real ALLOCABLE_, dimension(NIMEMW_,NJMEMBW_) :: &
    Datv_res, &     ! A nondimensional factor by which the meridional face
                    ! areas are to be rescaled to account for the effective
                    ! face areas of the layers, if rescale_D_bt is true.
                    ! Datv_res is set in btcalc.
    D_v_Cor, &      !   A simply averaged depth at v points, in m.
    dx_Cv, &         !   A copy of G%dx_Cv with wide halos, in m.
    IdyCv            !   A copy of G%IdyCv with wide halos, in m-1.
  real ALLOCABLE_, dimension(NIMEMBW_,NJMEMBW_) :: &
    q_D             ! f / D at PV points, in m-1 s-1.

  real, pointer, dimension(:,:,:) :: frhatu1, frhatv1 ! Predictor values.

  real    :: Rho0            !   The density used in the Boussinesq
                             ! approximation, in kg m-3.
  real    :: dtbt            ! The barotropic time step, in s.
  real    :: dtbt_fraction   !   The fraction of the maximum time-step that
                             ! should used.  The default is 0.98.
  real    :: dtbt_max        !   The maximum stable barotropic time step, in s.
  real    :: dt_bt_filter    !   The time-scale over which the barotropic mode
                             ! solutions are filtered, in s.  This can never
                             ! be taken to be longer than 2*dt.  The default, 0,
                             ! applies no filtering.
  integer :: nstep_last = 0  ! The number of barotropic timesteps per baroclinic
                             ! time step the last time btstep was called.
  real    :: bebt            ! A nondimensional number, from 0 to 1, that
                             ! determines the gravity wave time stepping scheme.
                             ! 0.0 gives a forward-backward scheme, while 1.0
                             ! give backward Euler. In practice, bebt should be
                             ! of order 0.2 or greater.
  logical :: split           ! If true, use the split time stepping scheme.
  real    :: eta_source_limit  ! The fraction of the initial depth of the ocean
                             ! that can be added or removed to the bartropic
                             ! solution within a thermodynamic time step.  By
                             ! default this is 0 (i.e., no correction). Nondim.
  logical :: bound_BT_corr   ! If true, the magnitude of the fake mass source
                             ! in the barotropic equation that drives the two
                             ! estimates of the free surface height toward each
                             ! other is bounded to avoid driving corrective
                             ! velocities that exceed 0.1*MAXVEL.
  logical :: gradual_BT_ICs  ! If true, adjust the initial conditions for the
                             ! barotropic solver to the values from the layered
                             ! solution over a whole timestep instead of
                             ! instantly.  This is a decent approximation to the
                             ! inclusion of sum(u dh_dt) while also correcting
                             ! for truncation errors.
  logical :: Sadourny        ! If true, the Coriolis terms are discretized
                             ! with Sadourny's energy conserving scheme,
                             ! otherwise the Arakawa & Hsu scheme is used.  If
                             ! the deformation radius is not resolved Sadourny's
                             ! scheme should probably be used.
  logical :: Nonlinear_continuity ! If true, the barotropic continuity equation
                             ! uses the full ocean thickness for transport.
  integer :: Nonlin_cont_update_period ! The number of barotropic time steps
                             ! between updates to the face area, or 0 only to
                             ! update at the start of a call to btstep.  The
                             ! default is 1.
  logical :: BT_project_velocity ! If true, step the barotropic velocity first
                             ! and project out the velocity tendancy by 1+BEBT
                             ! when calculating the transport.  The default
                             ! (false) is to use a predictor continuity step to
                             ! find the pressure field, and then do a corrector
                             ! continuity step using a weighted average of the
                             ! old and new velocities, with weights of (1-BEBT)
                             ! and BEBT.
  logical :: dynamic_psurf   ! If true, add a dynamic pressure due to a viscous
                             ! ice shelf, for instance.
  real    :: Dmin_dyn_psurf  ! The minimum depth to use in limiting the size
                             ! of the dynamic surface pressure for stability,
                             ! in m.
  real    :: ice_strength_length  ! The length scale at which the damping rate
                             ! due to the ice strength should be the same as if
                             ! a Laplacian were applied, in m.
  real    :: const_dyn_psurf ! The constant that scales the dynamic surface
                             ! pressure, nondim.  Stable values are < ~1.0.
                             ! The default is 0.9.
  logical :: tides           ! If true, apply tidal momentum forcing.
  real    :: G_extra         ! A nondimensional factor by which gtot is enhanced.
  real    :: drag_amp        ! A nondimensional value (presumably between 0 and
                             ! 1) scaling the magnitude of the bottom drag
                             ! applied in the barotropic solver, 1 by default.
  integer :: hvel_scheme     ! An integer indicating how the thicknesses at
                             ! velocity points are calculated. Valid values are
                             ! given by the parameters defined below:
                             !   HARMONIC, ARITHMETIC, HYBRID, and FROM_BT_CONT
  logical :: strong_drag     ! If true, use a stronger estimate of the retarding
                             ! effects of strong bottom drag.
  logical :: rescale_D_bt    ! If true, the open face areas are rescaled by a
                             ! function of the ratio of the summed harmonic
                             ! mean thicknesses to the harmonic mean of the
                             ! summed thicknesses.
  logical :: linearized_BT_PV  ! If true, the PV and interface thicknesses used
                             ! in the barotropic Coriolis calculation is time
                             ! invariant and linearized.
  logical :: use_wide_halos  ! If true, use wide halos and march in during the
                             ! barotropic time stepping for efficiency.
  logical :: clip_velocity   ! If true, limit any velocity components that
                             ! exceed maxvel.  This should only be used as a
                             ! desperate debugging measure.
  logical :: debug           ! If true, write verbose checksums for debugging purposes.
  logical :: debug_bt        ! If true, write verbose checksums for debugging purposes.
  real    :: maxvel          ! Velocity components greater than maxvel are
                             ! truncated to maxvel, in m s-1.
  real    :: maxCFL_BT_cont  ! The maximum permitted CFL number associated with the
                             ! barotropic accelerations from the summed velocities
                             ! times the time-derivatives of thicknesses.  The
                             ! default is 0.1, and there will probably be real
                             ! problems if this were set close to 1.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  type(MOM_domain_type), pointer :: BT_Domain => NULL()
  type(hor_index_type), pointer :: debug_BT_HI ! debugging copy of horizontal index_type
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL()
  logical :: module_is_initialized = .false.

  integer :: isdw, iedw, jsdw, jedw ! The memory limits of the wide halo arrays.

  integer :: id_PFu_bt = -1, id_PFv_bt = -1, id_Coru_bt = -1, id_Corv_bt = -1
  integer :: id_ubtforce = -1, id_vbtforce = -1, id_uaccel = -1, id_vaccel = -1
  integer :: id_visc_rem_u = -1, id_visc_rem_v = -1, id_eta_cor = -1
  integer :: id_ubt = -1, id_vbt = -1, id_eta_bt = -1, id_ubtav = -1, id_vbtav = -1
  integer :: id_ubt_st = -1, id_vbt_st = -1, id_eta_st = -1
  integer :: id_ubt_hifreq = -1, id_vbt_hifreq = -1, id_eta_hifreq = -1
  integer :: id_uhbt_hifreq = -1, id_vhbt_hifreq = -1, id_eta_pred_hifreq = -1
  integer :: id_gtotn = -1, id_gtots = -1, id_gtote = -1, id_gtotw = -1
  integer :: id_Datu_res = -1, id_Datv_res = -1
  integer :: id_uhbt = -1, id_frhatu = -1, id_vhbt = -1, id_frhatv = -1
  integer :: id_frhatu1 = -1, id_frhatv1 = -1
end type legacy_barotropic_CS

type, private :: local_BT_cont_u_type
  real :: FA_u_EE, FA_u_E0, FA_u_W0, FA_u_WW
  real :: ubt_EE, ubt_WW
  real :: uh_crvE, uh_crvW
  real :: uh_EE, uh_WW
end type local_BT_cont_u_type
type, private :: local_BT_cont_v_type
  real :: FA_v_NN, FA_v_N0, FA_v_S0, FA_v_SS
  real :: vbt_NN, vbt_SS
  real :: vh_crvN, vh_crvS
  real :: vh_NN, vh_SS
end type local_BT_cont_v_type

type, private :: memory_size_type
  integer :: isdw, iedw, jsdw, jedw ! The memory limits of the wide halo arrays.
end type memory_size_type


type, private :: BT_OBC_type
  real, dimension(:,:), pointer :: &
    Cg_u => NULL(), &     ! The external wave speed at u-points, in m s-1.
    Cg_v => NULL(), &     ! The external wave speed at u-points, in m s-1.
    H_u => NULL(), &      ! The total thickness at the u-points, in m or kg m-2.
    H_v => NULL(), &      ! The total thickness at the v-points, in m or kg m-2.
    uhbt => NULL(), &     ! The zonal and meridional barotropic thickness fluxes
    vhbt => NULL(), &     ! specified for open boundary conditions (if any),
                          ! in units of m3 s-1.
    ubt_outer => NULL(), & ! The zonal and meridional velocities just outside
    vbt_outer => NULL(), & ! the domain, as set by the open boundary conditions,
                           ! in units of m s-1.
    eta_outer_u => NULL(), & ! The surface height outside of the domain at a
    eta_outer_v => NULL()    ! u- or v- point with an open boundary condition,
                             ! in units of m or kg m-2.
  integer :: is_u_obc, ie_u_obc, js_u_obc, je_u_obc
  integer :: is_v_obc, ie_v_obc, js_v_obc, je_v_obc
end type BT_OBC_type

integer :: id_clock_sync=-1, id_clock_calc=-1
integer :: id_clock_calc_pre=-1, id_clock_calc_post=-1
integer :: id_clock_pass_step=-1, id_clock_pass_pre=-1, id_clock_pass_post=-1
logical :: apply_u_OBCs, apply_v_OBCs

! Enumeration values for various schemes
integer, parameter :: HARMONIC        = 1
integer, parameter :: ARITHMETIC      = 2
integer, parameter :: HYBRID          = 3
integer, parameter :: FROM_BT_CONT    = 4
integer, parameter :: HYBRID_BT_CONT  = 5
character*(20), parameter :: HYBRID_STRING = "HYBRID"
character*(20), parameter :: HARMONIC_STRING = "HARMONIC"
character*(20), parameter :: ARITHMETIC_STRING = "ARITHMETIC"
character*(20), parameter :: BT_CONT_STRING = "FROM_BT_CONT"

contains

subroutine legacy_btstep(use_fluxes, U_in, V_in, eta_in, dt, bc_accel_u, bc_accel_v, &
                  fluxes, pbce, eta_PF_in, U_Cor, V_Cor, &
                  accel_layer_u, accel_layer_v, eta_out, uhbtav, vhbtav, G, GV, CS, &
                  visc_rem_u, visc_rem_v, etaav, uhbt_out, vhbt_out, OBC, &
                  BT_cont, eta_PF_start, &
                  taux_bot, tauy_bot, uh0, vh0, u_uh0, v_vh0)
  type(ocean_grid_type),                   intent(inout) :: G
  type(verticalGrid_type),                 intent(in)    :: GV
  logical,                                 intent(in)    :: use_fluxes
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)  :: U_in
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)  :: V_in
  real, dimension(SZI_(G),SZJ_(G)),        intent(in)    :: eta_in
  real,                                    intent(in)    :: dt
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)  :: bc_accel_u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)  :: bc_accel_v
  type(forcing),                           intent(in)    :: fluxes
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)   :: pbce
  real, dimension(SZI_(G),SZJ_(G)),        intent(in)    :: eta_PF_in
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)  :: U_Cor
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)  :: V_Cor
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out) :: accel_layer_u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out) :: accel_layer_v
  real, dimension(SZI_(G),SZJ_(G)),        intent(inout) :: eta_out
  real, dimension(SZIB_(G),SZJ_(G)),       intent(out)   :: uhbtav
  real, dimension(SZI_(G),SZJB_(G)),       intent(out)   :: vhbtav
  type(legacy_barotropic_CS),              pointer       :: CS
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in), optional :: visc_rem_u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in), optional :: visc_rem_v
  real, dimension(SZI_(G),SZJ_(G)),    intent(out), optional :: etaav
  real, dimension(SZIB_(G),SZJ_(G)),   intent(out), optional :: uhbt_out
  real, dimension(SZI_(G),SZJB_(G)),   intent(out), optional :: vhbt_out
  type(ocean_OBC_type),                pointer,     optional :: OBC
  type(BT_cont_type),                  pointer,     optional :: BT_cont
  real, dimension(:,:),                pointer,     optional :: eta_PF_start
  real, dimension(:,:),                pointer,     optional :: taux_bot
  real, dimension(:,:),                pointer,     optional :: tauy_bot
  real, dimension(:,:,:),              pointer,     optional :: uh0, u_uh0
  real, dimension(:,:,:),              pointer,     optional :: vh0, v_vh0

! Arguments: use_fluxes - A logical indicating whether velocities (false) or
!                         fluxes (true) are used to initialize the barotropic
!                         velocities.
!  (in)      U_in - The initial (3-D) zonal velocity or volume or mass fluxes,
!                   depending on flux_form, in m s-1 or m3 s-1 or kg s-1.
!  (in)      V_in - The initial (3-D) meridional velocity or volume/mass fluxes,
!                   depending on flux_form, in m s-1 or m3 s-1 or kg s-1.
!  (in)      eta_in - The initial barotropic free surface height anomaly or
!                     column mass anomaly, in m or kg m-2.
!  (in)      dt - The time increment to integrate over.
!  (in)      bc_accel_u - The zonal baroclinic accelerations, in m s-2.
!  (in)      bc_accel_v - The meridional baroclinic accelerations, in m s-2.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      pbce - The baroclinic pressure anomaly in each layer
!                   due to free surface height anomalies, in m2 H-1 s-2.
!  (in)      eta_PF_in - The 2-D eta field (either SSH anomaly or column mass
!                        anomaly) that was used to calculate the input pressure
!                        gradient accelerations (or its final value if
!                        eta_PF_start is provided, in m or kg m-2.
!      Note: eta_in, pbce, and eta_PF_in must have up-to-date halos.
!  (in)      U_Cor - The (3-D) zonal- and meridional- velocities or volume or
!  (in)      V_Cor   mass fluxes used to calculate the Coriolis terms in
!                    bc_accel_u and bc_accel_v, in m s-1 or m3 s-1 or kg s-1.
!  (out)     accel_layer_u - The accelerations of each layer due to the
!  (out)     accel_layer_v - barotropic calculation, in m s-2.
!  (out)     eta_out - The final barotropic free surface height anomaly or
!                      column mass anomaly, in m or kg m-2.
!  (out)     uhbtav - the barotropic zonal volume or mass fluxes averaged
!                     through the barotropic steps, in m3 s-1 or kg s-1.
!  (out)     vhbtav - the barotropic meridional volume or mass fluxes averaged
!                     through the barotropic steps, in m3 s-1 or kg s-1.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 barotropic_init.
!  (in,opt)  visc_rem_u - Both the fraction of the momentum originally in a
!  (in,opt)  visc_rem_v - layer that remains after a time-step of viscosity,
!                         and the fraction of a time-step's worth of a
!                         barotropic acceleration that a layer experiences
!                         after viscosity is applied, in the zonal (_u) and
!                         meridional (_v) directions.  Nondimensional between
!                         0 (at the bottom) and 1 (far above the bottom).
!  (out,opt) etaav - The free surface height or column mass averaged over the
!                    barotropic integration, in m or kg m-2.
!  (out,opt) uhbt_out - The barotropic zonal volume or mass fluxes at the end
!                     of the barotropic steps, in m3 s-1 or kg s-1.
!  (out,opt) vhbt_out - The barotropic meridional volume or mass fluxes at the
!                     end of the barotropic steps, in m3 s-1 or kg s-1.
!  (in,opt)  OBC - An open boundary condition type, which contains the values
!                  associated with open boundary conditions.
!  (in,opt)  BT_cont - A structure with elements that describe the effective
!                      open face areas as a function of barotropic flow.
!  (in,opt)  eta_PF_start - The eta field consistent with the pressure gradient
!                      at the start of the barotropic stepping, in m or kg m-2.
!  (in,opt)  taux_bot - The zonal bottom frictional stress from ocean to the
!                       seafloor, in Pa.
!  (in,opt)  tauy_bot - The meridional bottom frictional stress from ocean to
!                       the seafloor, in Pa.

!    This subroutine time steps the barotropic equations explicitly.
!  For gravity waves, anything between a forwards-backwards scheme
!  and a simulated backwards Euler scheme is used, with bebt between
!  0.0 and 1.0 determining the scheme.  In practice, bebt must be of
!  order 0.2 or greater.  A forwards-backwards treatment of the
!  Coriolis terms is always used.
!    Depending on the value of use_fluxes, the initial velocities are determined
!  from input velocites or volume (Boussinesq) or mass (non-Boussinesq) fluxes.

  real :: ubt_Cor(SZIB_(G),SZJ_(G)) ! The barotropic velocities that had been
  real :: vbt_Cor(SZI_(G),SZJB_(G)) ! use to calculate the input Coriolis
                                    ! terms, in m s-1.
  real :: wt_u(SZIB_(G),SZJ_(G),SZK_(G)) ! wt_u and wt_v are the
  real :: wt_v(SZI_(G),SZJB_(G),SZK_(G)) ! normalized weights to
                ! be used in calculating barotropic velocities, possibly with
                ! sums less than one due to viscous losses.  Nondimensional.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    av_rem_u, &   ! The weighted average of visc_rem_u, if it is present. ND.
    tmp_u         ! A temporary array at u points.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    av_rem_v, &   ! The weighted average of visc_rem_u, if it is present. ND.
    tmp_v         ! A temporary array at v points.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    e_anom        ! The anomaly in the sea surface height or column mass
                  ! averaged between the beginning and end of the time step,
                  ! relative to eta_PF, with SAL effects included, in units
                  ! of H (m or kg m-2, the same as eta and h).

  ! These are always allocated with symmetric memory and wide halos.
  real :: q(SZIBW_(CS),SZJBW_(CS))  ! A pseudo potential vorticity in s-1 m-1.
  real, dimension(SZIBW_(CS),SZJW_(CS)) :: &
    ubt, &        ! The zonal barotropic velocity in m s-1.
    bt_rem_u, &   ! The fraction of the barotropic zonal velocity that remains
                  ! after a time step, the remainder being lost to bottom drag.
                  ! bt_rem_u is a nondimensional number between 0 and 1.
    BT_force_u, & ! The vertical average of all of the u-accelerations that are
                  ! not explicitly included in the barotropic equation, m s-2.
    u_accel_bt, & ! The difference between the zonal acceleration from the
                  ! barotropic calculation and BT_force_u, in m s-2.
    uhbt, &       ! The zonal barotropic thickness fluxes, in H m2 s-1.
    uhbt0, &      ! The difference between the sum of the layer zonal thickness
                  ! fluxes and the barotropic thickness flux using the same
                  ! velocity, in H m2 s-1.
    ubt_old, &    ! The starting value of ubt in a barotropic step, in m s-1.
    ubt_sum, &    ! The sum of ubt over the time steps, in m s-1.
    uhbt_sum, &   ! The sum of uhbt over the time steps, in H m2 s-1.
    ubt_wtd, &    ! A weighted sum used to find the filtered final ubt, in m s-1.
    ubt_trans, &  ! The latest value of ubt used for a transport, in m s-1.
    azon, bzon, & ! _zon & _mer are the values of the Coriolis force which
    czon, dzon, & ! are applied to the neighboring values of vbtav & ubtav,
    amer, bmer, & ! respectively to get the barotropic inertial rotation,
    cmer, dmer, & ! in units of s-1.
    Cor_ref_u, &  ! The zonal barotropic Coriolis acceleration due
                  ! to the reference velocities, in m s-2.
    PFu_bt_sum, & ! The summed zonal barotropic pressure gradient force, in m s-2.
    Coru_bt_sum, & ! The summed zonal barotropic Coriolis acceleration, in m s-2.
    DCor_u, &     ! A simply averaged depth at u points, in m.
    Datu          ! Basin depth at u-velocity grid points times the y-grid
                  ! spacing, in H m.
  real, dimension(SZIW_(CS),SZJBW_(CS)) :: &
    vbt, &        ! The meridional barotropic velocity in m s-1.
    bt_rem_v, &   ! The fraction of the barotropic meridional velocity that
                  ! remains after a time step, the rest being lost to bottom
                  ! drag.  bt_rem_v is a nondimensional number between 0 and 1.
    BT_force_v, & ! The vertical average of all of the v-accelerations that are
                  ! not explicitly included in the barotropic equation, m s-2.
    v_accel_bt, & ! The difference between the meridional acceleration from the
                  ! barotropic calculation and BT_force_v, in m s-2.
    vhbt, &       ! The meridional barotropic thickness fluxes, in H m2 s-1.
    vhbt0, &      ! The difference between the sum of the layer meridional
                  ! thickness fluxes and the barotropic thickness flux using
                  ! the same velocities, in H m2 s-1.
    vbt_old, &    ! The starting value of vbt in a barotropic step, in m s-1.
    vbt_sum, &    ! The sum of vbt over the time steps, in m s-1.
    vhbt_sum, &   ! The sum of vhbt over the time steps, in H m2 s-1.
    vbt_wtd, &    ! A weighted sum used to find the filtered final vbt, in m s-1.
    vbt_trans, &  ! The latest value of vbt used for a transport, in m s-1.
    Cor_ref_v, &  ! The meridional barotropic Coriolis acceleration due
                  ! to the reference velocities, in m s-2.
    PFv_bt_sum, & ! The summed meridional barotropic pressure gradient force,
                  ! in m s-2.
    Corv_bt_sum, & ! The summed meridional barotropic Coriolis acceleration,
                  ! in m s-2.
    DCor_v, &     ! A simply averaged depth at v points, in m.
    Datv          ! Basin depth at v-velocity grid points times the x-grid
                  ! spacing, in H m.
  real, target, dimension(SZIW_(CS),SZJW_(CS)) :: &
    eta, &        ! The barotropic free surface height anomaly or column mass
                  ! anomaly, in m or kg m-2.
    eta_pred      ! A predictor value of eta, in m or kg m-2 like eta.
  real, pointer, dimension(:,:) :: &
    eta_PF_BT     ! A pointer to the eta array (either eta or eta_pred) that
                  ! determines the barotropic pressure force, in m or kg m-2.
  real, dimension(SZIW_(CS),SZJW_(CS)) :: &
    eta_sum, &    ! eta summed across the timesteps, in m or kg m-2.
    eta_wtd, &    ! A weighted estimate used to calculate eta_out, in m or kg m-2.
    eta_PF, &     ! A local copy of the 2-D eta field (either SSH anomaly or
                  ! column mass anomaly) that was used to calculate the input
                  ! pressure gradient accelerations, in m or kg m-2.
    eta_PF_1, &   ! The initial value of eta_PF, when interp_eta_PF is
                  ! true, in m or kg m-2.
    d_eta_PF, &   ! The change in eta_PF over the barotropic time stepping when
                  ! interp_eta_PF is true, in m or kg m-2.
    gtot_E, &     ! gtot_X is the effective total reduced gravity used to relate
    gtot_W, &     ! free surface height deviations to pressure forces (including
    gtot_N, &     ! GFS and baroclinic  contributions) in the barotropic momentum
    gtot_S, &     ! equations half a grid-point in the X-direction (X is N, S,
                  ! E, or W) from the thickness point. gtot_X has units of m2 H-1 s-2.
                  ! (See Hallberg, J Comp Phys 1997 for a discussion.)
    eta_src, &    ! The source of eta per barotropic timestep, in m or kg m-2.
    dyn_coef_eta, & ! The coefficient relating the changes in eta to the
                  ! dynamic surface pressure under rigid ice, in m2 s-2 H-1.
    p_surf_dyn    ! A dynamic surface pressure under rigid ice, in m2 s-2.
  type(local_BT_cont_u_type), dimension(SZIBW_(CS),SZJW_(CS)) :: &
    BTCL_u        ! A repackaged version of the u-point information in BT_cont.
  type(local_BT_cont_v_type), dimension(SZIW_(CS),SZJBW_(CS)) :: &
    BTCL_v        ! A repackaged version of the v-point information in BT_cont.
  ! End of wide-sized variables.

  real :: I_Rho0      ! The inverse of the mean density (Rho0), in m3 kg-1.
  real :: visc_rem    ! A work variable that may equal visc_rem_[uv].  Nondim.
  real :: vel_prev    ! The previous velocity in m s-1.
  real :: vel_trans   ! The combination of the previous and current velocity
                      ! that does the mass transport, in m s-1.
  real :: dtbt        ! The barotropic time step in s.
  real :: bebt        ! A copy of CS%bebt.
  real :: be_proj     ! The fractional amount by which velocities are projected
                      ! when project_velocity is true. For now be_proj is set
                      ! to equal bebt, as they have similar roles and meanings.
  real :: Idt         ! The inverse of dt, in s-1.
  real :: det_de      ! The partial derivative due to self-attraction and loading
                      ! of the reference geopotential with the sea surface height.
                      ! This is typically ~0.09 or less.
  real :: dgeo_de     ! The constant of proportionality between geopotential and
                      ! sea surface height.  It is a nondimensional number of
                      ! order 1.  For stability, this may be made larger
                      ! than physical problem would suggest.
  real :: Instep      ! The inverse of the number of barotropic time steps
                      ! to take.
  real :: Cor, gradP  ! The Coriolis and pressure gradient accelerations, m s-1.
  real :: wt_end      ! The weighting of the final value of eta_PF, ND.
  integer :: nstep    ! The number of barotropic time steps to take.
  type(time_type) :: &
    time_bt_start, &  ! The starting time of the barotropic steps.
    time_step_end, &  ! The end time of a barotropic step.
    time_end_in       ! The end time for diagnostics when this routine started.
  real :: time_int_in ! The diagnostics' time interval when this routine started.
  logical :: do_hifreq_output  ! If true, output occurs every barotropic step.
  logical :: use_visc_rem, use_BT_cont
  logical :: do_ave, find_etaav, find_PF, find_Cor
  logical :: ice_is_rigid, nonblock_setup, interp_eta_PF
  logical :: project_velocity, add_uh0

  real :: dyn_coef_max ! The maximum stable value of dyn_coef_eta, in m2 s-2 H-1.
  real :: ice_strength = 0.0  ! The effective strength of the ice in m s-2.
  real :: Idt_max2    ! The squared inverse of the local maximum stable
                      ! barotropic time step, in s-2.
  real :: H_min_dyn   ! The minimum depth to use in limiting the size of the
                      ! dynamic surface pressure for stability, in H.
  real :: H_eff_dx2   ! The effective total thickness divided by the grid spacing
                      ! squared, in H m-2.
  real :: vel_tmp     ! A temporary velocity, in m s-1.

  real, allocatable, dimension(:) :: wt_vel, wt_eta, wt_accel, wt_trans, wt_accel2
  real :: sum_wt_vel, sum_wt_eta, sum_wt_accel, sum_wt_trans
  real :: I_sum_wt_vel, I_sum_wt_eta, I_sum_wt_accel, I_sum_wt_trans
  real :: dt_filt     ! The half-width of the barotropic filter, in s.
  integer :: nfilter

  logical :: apply_OBCs, apply_OBC_flather
  type(BT_OBC_type) :: BT_OBC  ! A structure with all of this module's fields
                               ! for applying open boundary conditions.
  type(memory_size_type) :: MS
  character(len=200) :: mesg
  integer :: pid_ubt, pid_eta, pid_e_anom, pid_etaav, pid_uhbtav, pid_ubtav
  integer :: pid_q, pid_eta_PF, pid_dyn_coef_eta, pid_eta_src
  integer :: pid_DCor_u, pid_Datu_res, pid_tmp_u, pid_gtot_E, pid_gtot_W
  integer :: pid_bt_rem_u, pid_Datu, pid_BT_force_u, pid_Cor_ref
  integer :: pid_eta_PF_1, pid_d_eta_PF, pid_uhbt0
  integer :: isv, iev, jsv, jev ! The valid array size at the end of a step.
  integer :: stencil  ! The stencil size of the algorithm, often 1 or 2.
  integer :: isvf, ievf, jsvf, jevf, num_cycles
  integer :: i, j, k, n
  integer :: is, ie, js, je, nz, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  if (.not.associated(CS)) call MOM_error(FATAL, &
      "legacy_btstep: Module MOM_legacy_barotropic must be initialized before it is used.")
  if (.not.CS%split) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  MS%isdw = CS%isdw ; MS%iedw = CS%iedw ; MS%jsdw = CS%jsdw ; MS%jedw = CS%jedw
  Idt = 1.0 / dt

  use_BT_cont = .false.
  if (present(BT_cont)) use_BT_cont = (associated(BT_cont))

  interp_eta_PF = .false.
  if (present(eta_PF_start)) interp_eta_PF = (associated(eta_PF_start))

  project_velocity = CS%BT_project_velocity

  ! Figure out the fullest arrays that could be updated.
  stencil = 1
  if ((.not.use_BT_cont) .and. CS%Nonlinear_continuity .and. &
      (CS%Nonlin_cont_update_period > 0)) stencil = 2

  num_cycles = 1
  if (CS%use_wide_halos) &
    num_cycles = min((is-CS%isdw) / stencil, (js-CS%jsdw) / stencil)
  isvf = is - (num_cycles-1)*stencil ; ievf = ie + (num_cycles-1)*stencil
  jsvf = js - (num_cycles-1)*stencil ; jevf = je + (num_cycles-1)*stencil

  do_ave = query_averaging_enabled(CS%diag)
  find_etaav = present(etaav)
  use_visc_rem = present(visc_rem_u)
  if ((use_visc_rem) .neqv. present(visc_rem_v)) call MOM_error(FATAL, &
      "btstep: Either both visc_rem_u and visc_rem_v or neither"// &
       " one must be present in call to btstep.")
  find_PF = (do_ave .and. ((CS%id_PFu_bt > 0) .or. (CS%id_PFv_bt > 0)))
  find_Cor = (do_ave .and. ((CS%id_Coru_bt > 0) .or. (CS%id_Corv_bt > 0)))

  add_uh0 = .false.
  if (present(uh0)) add_uh0 = associated(uh0)
  if (add_uh0 .and. .not.(present(vh0) .and. present(u_uh0) .and. &
                          present(v_vh0))) call MOM_error(FATAL, &
      "legacy_btstep: vh0, u_uh0, and v_vh0 must be present if uh0 is used.")
  if (add_uh0 .and. .not.(associated(vh0) .and. associated(u_uh0) .and. &
                          associated(v_vh0))) call MOM_error(FATAL, &
      "legacy_btstep: vh0, u_uh0, and v_vh0 must be associated if uh0 is used.")
  if (add_uh0 .and. use_fluxes) call MOM_error(WARNING, &
      "legacy_btstep: with use_fluxes, add_uh0 does nothing!")
  if (use_fluxes) add_uh0 = .false.

  ! This can be changed to try to optimize the performance.
  nonblock_setup = G%nonblocking_updates

  if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)

  apply_OBCs = .false. ; apply_u_OBCs = .false. ; apply_v_OBCs = .false.
  apply_OBC_flather = .false.
  if (present(OBC)) then ; if (associated(OBC)) then
    apply_u_OBCs = OBC%Flather_u_BCs_exist_globally .or. OBC%specified_u_BCs_exist_globally
    apply_v_OBCs = OBC%Flather_v_BCs_exist_globally .or. OBC%specified_v_BCs_exist_globally
    apply_OBC_flather = OBC%Flather_u_BCs_exist_globally .or. OBC%Flather_v_BCs_exist_globally
    apply_OBCs = OBC%specified_u_BCs_exist_globally .or. OBC%specified_v_BCs_exist_globally .or. apply_OBC_flather

    if (apply_OBC_flather .and. .not.GV%Boussinesq) call MOM_error(FATAL, &
      "legacy_btstep: Flather open boundary conditions have not yet been "// &
      "implemented for a non-Boussinesq model.")
  endif ; endif

  nstep = CEILING(dt/CS%dtbt - 0.0001)
  if (is_root_PE() .and. (nstep /= CS%nstep_last)) then
    write(mesg,'("legacy_btstep is using a dynamic barotropic timestep of ", ES12.6, &
               & " seconds, max ", ES12.6, ".")') (dt/nstep), CS%dtbt_max
    call MOM_mesg(mesg, 3)
  endif
  CS%nstep_last = nstep

  ! Set the actual barotropic time step.
  Instep = 1.0 / real(nstep)
  dtbt = dt * Instep
  bebt = CS%bebt
  be_proj = CS%bebt
  I_Rho0 = 1.0/GV%Rho0
  do_ave = query_averaging_enabled(CS%diag)

  do_hifreq_output = .false.
  if ((CS%id_ubt_hifreq > 0) .or. (CS%id_vbt_hifreq > 0) .or. &
      (CS%id_eta_hifreq > 0) .or. (CS%id_eta_pred_hifreq > 0) .or. &
      (CS%id_uhbt_hifreq > 0) .or. (CS%id_vhbt_hifreq > 0)) then
    do_hifreq_output = query_averaging_enabled(CS%diag, time_int_in, time_end_in)
    if (do_hifreq_output) &
      time_bt_start = time_end_in - set_time(int(floor(dt+0.5)))
  endif

!   Calculate the constant coefficients for the Coriolis force terms in the
! barotropic momentum equations.  This has to be done quite early to start
! the halo update that needs to be completed before the next calculations.
  if (CS%linearized_BT_PV) then
!$OMP parallel do default(none) shared(jsvf,jevf,isvf,ievf,q,CS)
    do J=jsvf-2,jevf+1 ; do I=isvf-2,ievf+1
      q(I,J) = CS%q_D(I,j)
    enddo ; enddo
!$OMP parallel do default(none) shared(jsvf,jevf,isvf,ievf,DCor_u,CS)
    do j=jsvf-1,jevf+1 ; do I=isvf-2,ievf+1
      DCor_u(I,j) = CS%D_u_Cor(I,j)
    enddo ; enddo
!$OMP parallel do default(none) shared(jsvf,jevf,isvf,ievf,DCor_v,CS)
    do J=jsvf-2,jevf+1 ; do i=isvf-1,ievf+1
      DCor_v(i,J) = CS%D_v_Cor(i,J)
    enddo ; enddo
  else
    q(:,:) = 0.0 ; DCor_u(:,:) = 0.0 ; DCor_v(:,:) = 0.0
    !  This option has not yet been written properly.
    !  D here should be replaced with D+eta(Bous) or eta(non-Bous).
!$OMP parallel do default(none) shared(js,je,is,ie,DCor_u,G)
    do j=js,je ; do I=is-1,ie
      DCor_u(I,j) = 0.5 * (G%bathyT(i+1,j) + G%bathyT(i,j))
    enddo ; enddo
!$OMP parallel do default(none) shared(js,je,is,ie,DCor_v,G)
    do J=js-1,je ; do i=is,ie
      DCor_v(i,J) = 0.5 * (G%bathyT(i,j+1) + G%bathyT(i,j))
    enddo ; enddo
!$OMP parallel do default(none) shared(js,je,is,ie,q,G)
    do J=js-1,je ; do I=is-1,ie
      q(I,J) = 0.25 * G%CoriolisBu(I,J) * &
           ((G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i+1,j) + G%areaT(i,j+1))) / &
           ((G%areaT(i,j) * G%bathyT(i,j) + G%areaT(i+1,j+1) * G%bathyT(i+1,j+1)) + &
            (G%areaT(i+1,j) * G%bathyT(i+1,j) + G%areaT(i,j+1) * G%bathyT(i,j+1)))
    enddo ; enddo
    ! With very wide halos, q and D need to be calculated on the available data
    ! domain and then updated onto the full computational domain.
    ! These calculations can be done almost immediately, but the halo updates
    ! must be done before the [abcd]mer and [abcd]zon are calculated.
    if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
    if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
    if (nonblock_setup) then
      pid_q = pass_var_start(q, CS%BT_Domain, To_All, position=CORNER)
      pid_DCor_u = pass_vector_start(DCor_u, DCor_v, CS%BT_Domain, To_All+Scalar_Pair)
    else
      call pass_var(q, CS%BT_Domain, To_All, position=CORNER)
      call pass_vector(DCor_u, DCor_v, CS%BT_Domain, To_All+Scalar_Pair)
    endif
    if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
    if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  endif

  if (nonblock_setup) then
    ! Start all halo updates that are ready to go.
  !###   if (use_BT_cont) call start_set_local_BT_cont_types( ... )
    if ((.not.use_BT_cont) .and. CS%rescale_D_bt .and. (ievf>ie)) then
      if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
      if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
      pid_Datu_res = pass_vector_start(CS%Datu_res, CS%Datv_res, CS%BT_Domain, &
                                       To_All+Scalar_Pair)
      if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
      if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
    endif
  endif

  ! Zero out various wide-halo arrays.
  do j=CS%jsdw,CS%jedw ; do i=CS%isdw,CS%iedw
    gtot_E(i,j) = 0.0 ; gtot_W(i,j) = 0.0
    gtot_N(i,j) = 0.0 ; gtot_S(i,j) = 0.0
    eta(i,j) = 0.0
    eta_PF(i,j) = 0.0
    if (interp_eta_PF) then
      eta_PF_1(i,j) = 0.0 ; d_eta_PF(i,j) = 0.0
    endif
    p_surf_dyn(i,j) = 0.0
    if (CS%dynamic_psurf) dyn_coef_eta(i,j) = 0.0
  enddo ; enddo
  !   The halo regions of various arrays need to be initialized to
  ! non-NaNs in case the neighboring domains are not part of the ocean.
  ! Otherwise a halo update later on fills in the correct values.
  do j=CS%jsdw,CS%jedw ; do I=CS%isdw-1,CS%iedw
    Cor_ref_u(I,j) = 0.0 ; BT_force_u(I,j) = 0.0 ; ubt(I,j) = 0.0
    Datu(I,j) = 0.0 ; bt_rem_u(I,j) = 0.0 ; uhbt0(I,j) = 0.0
  enddo ; enddo
  do J=CS%jsdw-1,CS%jedw ; do i=CS%isdw,CS%iedw
    Cor_ref_v(i,J) = 0.0 ; BT_force_v(i,J) = 0.0 ; vbt(i,J) = 0.0
    Datv(i,J) = 0.0 ; bt_rem_v(i,J) = 0.0 ; vhbt0(I,j) = 0.0
  enddo ; enddo

  ! Copy input arrays into their wide-halo counterparts.  eta_in and eta_PF_in
  ! have the correct values in their halo regions.
  if (interp_eta_PF) then
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      eta(i,j) = eta_in(i,j)
      eta_PF_1(i,j) = eta_PF_start(i,j)
      d_eta_PF(i,j) = eta_PF_in(i,j) - eta_PF_start(i,j)
    enddo ; enddo
  else
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      eta(i,j) = eta_in(i,j)
      eta_PF(i,j) = eta_PF_in(i,j)
    enddo ; enddo
  endif

  if (use_visc_rem) then
!$OMP parallel do default(none) shared(Isq,Ieq,js,je,nz,visc_rem_u,Instep,wt_u,CS) private(visc_rem)
    do k=1,nz ; do j=js-1,je+1 ; do I=Isq-1,Ieq+1
      ! rem needs greater than visc_rem_u and 1-Instep/visc_rem_u.
      ! The 0.5 below is just for safety.
      if (visc_rem_u(I,j,k) <= 0.0) then ; visc_rem = 0.0
      elseif (visc_rem_u(I,j,k) >= 1.0) then ; visc_rem = 1.0
      elseif (visc_rem_u(I,j,k)**2 > visc_rem_u(I,j,k) - 0.5*Instep) then
        visc_rem = visc_rem_u(I,j,k)
      else ; visc_rem = 1.0 - 0.5*Instep/visc_rem_u(I,j,k) ; endif
      wt_u(I,j,k) = CS%frhatu(I,j,k) * visc_rem
    enddo ; enddo ; enddo
!$OMP parallel do default(none) shared(is,ie,Jsq,Jeq,nz,visc_rem_v,Instep,wt_v,CS) private(visc_rem)
    do k=1,nz ; do J=Jsq-1,Jeq+1 ; do i=is-1,ie+1
      ! rem needs greater than visc_rem_v and 1-Instep/visc_rem_v.
      if (visc_rem_v(i,J,k) <= 0.0) then ; visc_rem = 0.0
      elseif (visc_rem_v(i,J,k) >= 1.0) then ; visc_rem = 1.0
      elseif (visc_rem_v(i,J,k)**2 > visc_rem_v(i,J,k) - 0.5*Instep) then
        visc_rem = visc_rem_v(i,J,k)
      else ; visc_rem = 1.0 - 0.5*Instep/visc_rem_v(i,J,k) ; endif
      wt_v(i,J,k) = CS%frhatv(i,J,k) * visc_rem
    enddo ; enddo ; enddo
  else
    do k=1,nz ; do j=js-1,je+1 ; do I=Isq-1,Ieq+1
      wt_u(I,j,k) = CS%frhatu(I,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do J=Jsq-1,Jeq+1 ; do i=is-1,ie+1
      wt_v(i,J,k) = CS%frhatv(i,J,k)
    enddo ; enddo ; enddo
  endif

  ! The gtot arrays are the effective layer-weighted reduced gravities for
  ! accelerations across the various faces, with names for the relative
  ! locations of the faces to the pressure point.  They will have their halos
  ! updated later on.
  do k=1,nz
    do j=js,je ; do I=is-1,ie
      gtot_E(i,j)   = gtot_E(i,j)   + pbce(i,j,k)   * wt_u(I,j,k)
      gtot_W(i+1,j) = gtot_W(i+1,j) + pbce(i+1,j,k) * wt_u(I,j,k)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      gtot_N(i,j)   = gtot_N(i,j)   + pbce(i,j,k)   * wt_v(i,J,k)
      gtot_S(i,j+1) = gtot_S(i,j+1) + pbce(i,j+1,k) * wt_v(i,J,k)
    enddo ; enddo
  enddo

  if (CS%tides) then
    call tidal_forcing_sensitivity(G, CS%tides_CSp, det_de)
    dgeo_de = 1.0 + det_de + CS%G_extra
  else
    dgeo_de = 1.0 + CS%G_extra
  endif

  if (nonblock_setup .and. .not.CS%linearized_BT_PV) then
    if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
    if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
      call pass_var_complete(pid_q, q, CS%BT_Domain, To_All, position=CORNER)
      call pass_vector_complete(pid_DCor_u, DCor_u, DCor_v, CS%BT_Domain, To_All+Scalar_Pair)
    if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
    if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  endif

  ! Calculate the open areas at the velocity points.
  ! The halo updates are needed before Datu is first used, either in set_up_BT_OBC or ubt_Cor.
  if (use_BT_cont) then
    call set_local_BT_cont_types(BT_cont, BTCL_u, BTCL_v, G, MS, CS%BT_Domain, 1+ievf-ie)
  else
    if (CS%rescale_D_bt .and. (ievf>ie)) then
      ! Datu_res was previously calculated in btcalc, and will be used in find_face_areas.
      ! This halo update is needed for wider halos than 1.  The complete goes here.
      if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
      if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
      if (nonblock_setup) then
        call pass_vector_complete(pid_Datu_res, CS%Datu_res, CS%Datv_res, &
                                  CS%BT_Domain, To_All+Scalar_Pair)
      else
        call pass_vector(CS%Datu_res, CS%Datv_res, CS%BT_Domain, To_All+Scalar_Pair)
      endif
      if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
      if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
    endif
    if (CS%Nonlinear_continuity) then
      call find_face_areas(Datu, Datv, G, GV, CS, MS, CS%rescale_D_bt, eta, 1)
    else
      call find_face_areas(Datu, Datv, G, GV, CS, MS, CS%rescale_D_bt, halo=1)
    endif
  endif

  ! Set up fields related to the open boundary conditions.
  if (apply_OBCs) then
    call set_up_BT_OBC(OBC, eta, BT_OBC, G, GV, MS, ievf-ie, use_BT_cont, &
                                     Datu, Datv, BTCL_u, BTCL_v)
  endif

!   Here the vertical average accelerations due to the Coriolis, advective,
! pressure gradient and horizontal viscous terms in the layer momentum
! equations are calculated.  These will be used to determine the difference
! between the accelerations due to the average of the layer equations and the
! barotropic calculation.
  !  ### Should IDatu here be replaced with 1/D+eta(Bous) or 1/eta(non-Bous)?
  if (use_visc_rem) then
    do j=js,je ; do I=is-1,ie
      BT_force_u(I,j) = fluxes%taux(I,j) * I_rho0*CS%IDatu(I,j)*visc_rem_u(I,j,1)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      BT_force_v(i,J) = fluxes%tauy(i,J) * I_rho0*CS%IDatv(i,J)*visc_rem_v(i,J,1)
    enddo ; enddo
  else
    do j=js,je ; do I=is-1,ie
      BT_force_u(I,j) = fluxes%taux(I,j) * I_rho0 * CS%IDatu(I,j)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      BT_force_v(i,J) = fluxes%tauy(i,J) * I_rho0 * CS%IDatv(i,J)
    enddo ; enddo
  endif
  if (present(taux_bot) .and. present(tauy_bot)) then
    if (associated(taux_bot) .and. associated(tauy_bot)) then
      do j=js,je ; do I=is-1,ie
        BT_force_u(I,j) = BT_force_u(I,j) - taux_bot(I,j) * I_rho0 * CS%IDatu(I,j)
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        BT_force_v(i,J) = BT_force_v(i,J) - tauy_bot(i,J) * I_rho0 * CS%IDatv(i,J)
      enddo ; enddo
    endif
  endif

  ! bc_accel_u & bc_accel_v are only available on the potentially
  ! non-symmetric computational domain.
  do k=1,nz ; do j=js,je ; do I=Isq,Ieq
    BT_force_u(I,j) = BT_force_u(I,j) + wt_u(I,j,k) * bc_accel_u(I,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
    BT_force_v(i,J) = BT_force_v(i,J) + wt_v(i,J,k) * bc_accel_v(i,J,k)
  enddo ; enddo ; enddo

  ! Determine the difference between the sum of the layer fluxes and the
  ! barotropic fluxes found from the same input velocities.
  if (add_uh0) then
    do j=js,je ; do I=is-1,ie ; uhbt(I,j) = 0.0 ; ubt(I,j) = 0.0 ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; vhbt(i,J) = 0.0 ; vbt(i,J) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      uhbt(I,j) = uhbt(I,j) + uh0(I,j,k)
      ubt(I,j) = ubt(I,j) + wt_u(I,j,k) * u_uh0(I,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      vhbt(i,J) = vhbt(i,J) + vh0(i,J,k)
      vbt(i,J) = vbt(i,J) + wt_v(i,J,k) * v_vh0(i,J,k)
    enddo ; enddo ; enddo
    if (use_BT_cont) then
      do j=js,je ; do I=is-1,ie
        uhbt0(I,j) = uhbt(I,j) - find_uhbt(ubt(I,j),BTCL_u(I,j))
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        vhbt0(i,J) = vhbt(i,J) - find_vhbt(vbt(i,J),BTCL_v(i,J))
      enddo ; enddo
    else
      do j=js,je ; do I=is-1,ie
        uhbt0(I,j) = uhbt(I,j) - Datu(I,j)*ubt(I,j)
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        vhbt0(i,J) = vhbt(i,J) - Datv(i,J)*vbt(i,J)
      enddo ; enddo
    endif
  endif

! Calculate the initial barotropic velocities from the layer's velocities.
  do j=jsvf-1,jevf+1 ; do I=isvf-2,ievf+1
    ubt(I,j) = 0.0 ; uhbt(I,j) = 0.0 ; u_accel_bt(I,j) = 0.0
  enddo ; enddo
  do J=jsvf-2,jevf+1 ; do i=isvf-1,ievf+1
    vbt(i,J) = 0.0 ; vhbt(i,J) = 0.0 ; v_accel_bt(i,J) = 0.0
  enddo ; enddo
  if (use_fluxes) then
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      uhbt(I,j) = uhbt(I,j) + U_in(I,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      vhbt(i,J) = vhbt(i,J) + V_in(i,J,k)
    enddo ; enddo ; enddo
    if (use_BT_cont) then
      do j=js,je ; do I=is-1,ie
        ubt(I,j) = uhbt_to_ubt(uhbt(I,j),BTCL_u(I,j), guess=CS%ubt_IC(I,j))
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        vbt(i,J) = vhbt_to_vbt(vhbt(i,J),BTCL_v(i,J), guess=CS%vbt_IC(i,J))
      enddo ; enddo
    else
      do j=js,je ; do I=is-1,ie
        ubt(I,j) = 0.0 ; if (Datu(I,j) > 0.0) ubt(I,j) = uhbt(I,j) / Datu(I,j)
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        vbt(i,J) = 0.0 ; if (Datv(i,J) > 0.0) vbt(i,J) = vhbt(i,J) / Datv(i,J)
      enddo ; enddo
    endif
  else
    do k=1,nz ; do j=js,je ; do I=is-1,ie
      ubt(I,j) = ubt(I,j) + wt_u(I,j,k) * U_in(I,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-1,je ; do i=is,ie
      vbt(i,J) = vbt(i,J) + wt_v(i,J,k) * V_in(i,J,k)
    enddo ; enddo ; enddo
  endif

  if (CS%gradual_BT_ICs) then
    if (use_fluxes) then
      if (use_BT_cont) then
        do j=js,je ; do I=is-1,ie
          vel_tmp = uhbt_to_ubt(CS%uhbt_IC(I,j),BTCL_u(I,j), guess=CS%ubt_IC(I,j))
          BT_force_u(I,j) = BT_force_u(I,j) + (ubt(I,j) - vel_tmp) * Idt
          ubt(I,j) = vel_tmp
        enddo ; enddo
        do J=js-1,je ; do i=is,ie
          vel_tmp = vhbt_to_vbt(CS%vhbt_IC(i,J),BTCL_v(i,J), guess=CS%vbt_IC(i,J))
          BT_force_v(i,J) = BT_force_v(i,J) + (vbt(i,J) - vel_tmp) * Idt
          vbt(i,J) = vel_tmp
        enddo ; enddo
      else
        do j=js,je ; do I=is-1,ie
          vel_tmp = 0.0 ; if (Datu(I,j) > 0.0) vel_tmp = CS%uhbt_IC(I,j) / Datu(I,j)
          BT_force_u(I,j) = BT_force_u(I,j) + (ubt(I,j) - vel_tmp) * Idt
          ubt(I,j) = vel_tmp
        enddo ; enddo
        do J=js-1,je ; do i=is,ie
          vel_tmp = 0.0 ; if (Datv(i,J) > 0.0) vel_tmp = CS%vhbt_IC(i,J) / Datv(i,J)
          BT_force_v(i,J) = BT_force_v(i,J) + (vbt(i,J) - vel_tmp) * Idt
          vbt(i,J) = vel_tmp
        enddo ; enddo
      endif
    else
      do j=js,je ; do I=is-1,ie
        BT_force_u(I,j) = BT_force_u(I,j) + (ubt(I,j) - CS%ubt_IC(I,j)) * Idt
        ubt(I,j) = CS%ubt_IC(I,j)
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        BT_force_v(i,J) = BT_force_v(i,J) + (vbt(i,J) - CS%vbt_IC(i,J)) * Idt
        vbt(i,J) = CS%vbt_IC(i,J)
      enddo ; enddo
    endif
  endif

  if ((Isq > is-1) .or. (Jsq > js-1)) then
    ! Non-symmetric memory is being used, so the edge values need to be
    ! filled in with a halo update of a non-symmetric array.
    if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
    if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
    tmp_u(:,:) = 0.0 ; tmp_v(:,:) = 0.0
    do j=js,je ; do I=Isq,Ieq ; tmp_u(I,j) = BT_force_u(I,j) ; enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie ; tmp_v(i,J) = BT_force_v(i,J) ; enddo ; enddo
    if (nonblock_setup) then
      pid_tmp_u = pass_vector_start(tmp_u, tmp_v, G%Domain)
    else
      call pass_vector(tmp_u, tmp_v, G%Domain, complete=.true.)
      do j=jsd,jed ; do I=IsdB,IedB ; BT_force_u(I,j) = tmp_u(I,j) ; enddo ; enddo
      do J=JsdB,JedB ; do i=isd,ied ; BT_force_v(i,J) = tmp_v(i,J) ; enddo ; enddo
    endif
    if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
    if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  endif

  if (nonblock_setup) then
    if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
    if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
    pid_gtot_E = pass_vector_start(gtot_E, gtot_N, CS%BT_Domain, &
                     To_All+Scalar_Pair, AGRID, complete=.false.)
    pid_gtot_W = pass_vector_start(gtot_W, gtot_S, CS%BT_Domain, &
                     To_All+Scalar_Pair, AGRID)
    if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
    if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  endif

  ! Determine the weighted Coriolis parameters for the neighboring velocities.
!$OMP parallel do default(none) shared(isvf,ievf,jsvf,jevf,amer,bmer,cmer,dmer,DCor_u,q,CS)
  do J=jsvf-1,jevf ; do i=isvf-1,ievf+1
    if (CS%Sadourny) then
      amer(I-1,j) = DCor_u(I-1,j) * q(I-1,J)
      bmer(I,j) = DCor_u(I,j) * q(I,J)
      cmer(I,j+1) = DCor_u(I,j+1) * q(I,J)
      dmer(I-1,j+1) = DCor_u(I-1,j+1) * q(I-1,J)
    else
      amer(I-1,j) = DCor_u(I-1,j) * &
                    ((q(I,J) + q(I-1,J-1)) + q(I-1,J)) / 3.0
      bmer(I,j) = DCor_u(I,j) * &
                  (q(I,J) + (q(I-1,J) + q(I,J-1))) / 3.0
      cmer(I,j+1) = DCor_u(I,j+1) * &
                    (q(I,J) + (q(I-1,J) + q(I,J+1))) / 3.0
      dmer(I-1,j+1) = DCor_u(I-1,j+1) * &
                      ((q(I,J) + q(I-1,J+1)) + q(I-1,J)) / 3.0
    endif
  enddo ; enddo

!$OMP parallel do default(none) shared(isvf,ievf,jsvf,jevf,azon,bzon,czon,dzon,DCor_v,q,CS)
  do j=jsvf-1,jevf+1 ; do I=isvf-1,ievf
    if (CS%Sadourny) then
      azon(I,j) = DCor_v(i+1,J) * q(I,J)
      bzon(I,j) = DCor_v(i,J) * q(I,J)
      czon(I,j) = DCor_v(i,J-1) * q(I,J-1)
      dzon(I,j) = DCor_v(i+1,J-1) * q(I,J-1)
    else
      azon(I,j) = DCor_v(i+1,J) * &
                  (q(I,J) + (q(I+1,J) + q(I,J-1))) / 3.0
      bzon(I,j) = DCor_v(i,J) * &
                  (q(I,J) + (q(I-1,J) + q(I,J-1))) / 3.0
      czon(I,j) = DCor_v(i,J-1) * &
                  ((q(I,J) + q(I-1,J-1)) + q(I,J-1)) / 3.0
      dzon(I,j) = DCor_v(i+1,J-1) * &
                  ((q(I,J) + q(I+1,J-1)) + q(I,J-1)) / 3.0
    endif
  enddo ; enddo

  !   If they are present, use u_Cor and v_Cor as the reference values for the
  ! Coriolis terms, including the viscous remnant if it is present.
  if (use_fluxes) then
    do j=js-1,je+1 ; do I=is-1,ie ; uhbt(I,j) = 0.0 ; enddo ; enddo
    do J=js-1,je ; do i=is-1,ie+1 ; vhbt(i,J) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js-1,je+1 ; do I=is-1,ie
      uhbt(I,j) = uhbt(I,j) + U_Cor(I,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-1,je ; do i=is-1,ie+1
      vhbt(i,J) = vhbt(i,J) + V_Cor(i,J,k)
    enddo ; enddo ; enddo
    if (use_BT_cont) then
      do j=js-1,je+1 ; do I=is-1,ie
        ubt_Cor(I,j) = uhbt_to_ubt(uhbt(I,j), BTCL_u(I,j), guess=CS%ubtav(I,j))
      enddo ; enddo
      do J=js-1,je ; do i=is-1,ie+1
        vbt_Cor(i,J) = vhbt_to_vbt(vhbt(i,J), BTCL_v(i,J), guess=CS%vbtav(i,J))
      enddo ; enddo
    else
      do j=js-1,je+1 ; do I=is-1,ie
        ubt_Cor(I,j) = 0.0 ; if (Datu(I,j)>0.0) ubt_Cor(I,j) = uhbt(I,j)/Datu(I,j)
      enddo ; enddo
      do J=js-1,je ; do i=is-1,ie+1
        vbt_Cor(i,J) = 0.0 ; if (Datv(i,J)>0.0) vbt_Cor(i,J) = vhbt(i,J)/Datv(i,J)
      enddo ; enddo
    endif
  else
    do j=js-1,je+1 ; do I=is-1,ie ; ubt_Cor(I,j) = 0.0 ; enddo ; enddo
    do J=js-1,je ; do i=is-1,ie+1 ; vbt_Cor(i,J) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js-1,je+1 ; do I=is-1,ie
      ubt_Cor(I,j) = ubt_Cor(I,j) + wt_u(I,j,k) * U_Cor(I,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-1,je ; do i=is-1,ie+1
      vbt_Cor(i,J) = vbt_Cor(i,J) + wt_v(i,J,k) * V_Cor(i,J,k)
    enddo ; enddo ; enddo
  endif

!$OMP parallel do default(none) shared(is,ie,js,je,Cor_ref_u,azon,bzon,czon,dzon,vbt_Cor)
  do j=js,je ; do I=is-1,ie
    Cor_ref_u(I,j) =  &
        ((azon(I,j) * vbt_Cor(i+1,j) + czon(I,j) * vbt_Cor(i  ,j-1)) + &
         (bzon(I,j) * vbt_Cor(i  ,j) + dzon(I,j) * vbt_Cor(i+1,j-1)))
  enddo ; enddo
!$OMP parallel do default(none) shared(is,ie,js,je,Cor_ref_v,amer,bmer,cmer,dmer,ubt_Cor)
  do J=js-1,je ; do i=is,ie
    Cor_ref_v(i,J) = -1.0 * &
        ((amer(I-1,j) * ubt_Cor(I-1,j) + cmer(I  ,j+1) * ubt_Cor(I  ,j+1)) + &
         (bmer(I  ,j) * ubt_Cor(I  ,j) + dmer(I-1,j+1) * ubt_Cor(I-1,j+1)))
  enddo ; enddo

! Complete the previously initiated message passing.
  if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
  if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
  if (nonblock_setup) then
    if ((Isq > is-1) .or. (Jsq > js-1)) then
      call pass_vector_complete(pid_tmp_u, tmp_u, tmp_v, G%Domain)
      do j=jsd,jed ; do I=IsdB,IedB ; BT_force_u(I,j) = tmp_u(I,j) ; enddo ; enddo
      do J=JsdB,JedB ; do i=isd,ied ; BT_force_v(i,J) = tmp_v(i,J) ; enddo ; enddo
    endif
    call pass_vector_complete(pid_gtot_E, gtot_E, gtot_N, CS%BT_Domain, To_All+Scalar_Pair, AGRID)
    call pass_vector_complete(pid_gtot_W, gtot_W, gtot_S, CS%BT_Domain, To_All+Scalar_Pair, AGRID)
  else
    call pass_vector(gtot_E, gtot_N, CS%BT_Domain, To_All+Scalar_Pair, AGRID, complete=.false.)
    call pass_vector(gtot_W, gtot_S, CS%BT_Domain, To_All+Scalar_Pair, AGRID, complete=.true.)
  endif
  ! The various elements of gtot are positive definite but directional, so use
  ! the polarity arrays to sort out when the directions have shifted.
  do j=jsvf-1,jevf+1 ; do i=isvf-1,ievf+1
    if (CS%ua_polarity(i,j) < 0.0) call swap(gtot_E(i,j), gtot_W(i,j))
    if (CS%va_polarity(i,j) < 0.0) call swap(gtot_N(i,j), gtot_S(i,j))
  enddo ; enddo

  ! Now start new halo updates.
  if (nonblock_setup) then
    if (.not.use_BT_cont) &
      pid_Datu = pass_vector_start(Datu, Datv, CS%BT_Domain, To_All+Scalar_Pair)
    ! The following halo update is not needed without wide halos.  RWH
    if (((G%isd > CS%isdw) .or. (G%jsd > CS%jsdw)) .or. (Isq <= is-1) .or. (Jsq <= js-1)) &
      pid_BT_force_u = pass_vector_start(BT_force_u, BT_force_v, CS%BT_Domain, &
                                         complete=.false.)
    if (add_uh0) pid_uhbt0 = pass_vector_start(uhbt0, vhbt0, CS%BT_Domain)
    pid_Cor_ref = pass_vector_start(Cor_ref_u, Cor_ref_v, CS%BT_Domain)
  endif
  if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
  if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)

  do j=js-1,je+1 ; do I=is-1,ie ; bt_rem_u(I,j) = G%mask2dCu(I,j) ; enddo ; enddo
  do J=js-1,je ; do i=is-1,ie+1 ; bt_rem_v(i,J) = G%mask2dCv(i,J) ; enddo ; enddo
  if ((use_visc_rem) .and. (CS%drag_amp > 0.0)) then
    do j=js-1,je+1 ; do I=is-1,ie ; av_rem_u(I,j) = 0.0 ; enddo ; enddo
    do J=js-1,je ; do i=is-1,ie+1 ; av_rem_v(i,J) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js-1,je+1 ; do I=is-1,ie
      av_rem_u(I,j) = av_rem_u(I,j) + CS%frhatu(I,j,k) * visc_rem_u(I,j,k)
    enddo ; enddo ; enddo
    do k=1,nz ; do J=js-1,je ; do i=is-1,ie+1
      av_rem_v(i,J) = av_rem_v(i,J) + CS%frhatv(i,J,k) * visc_rem_v(i,J,k)
    enddo ; enddo ; enddo
    if (CS%strong_drag) then
      do j=js-1,je+1 ; do I=is-1,ie
        bt_rem_u(I,j) = G%mask2dCu(I,j) * CS%drag_amp * &
           ((nstep * av_rem_u(I,j)) / (1.0 + (nstep-1)*av_rem_u(I,j)))
      enddo ; enddo
      do J=js-1,je ; do i=is-1,ie+1
        bt_rem_v(i,J) = G%mask2dCv(i,J) * CS%drag_amp * &
           ((nstep * av_rem_v(i,J)) / (1.0 + (nstep-1)*av_rem_v(i,J)))
      enddo ; enddo
    else
      do j=js-1,je+1 ; do I=is-1,ie
        bt_rem_u(I,j) = 0.0
        if (G%mask2dCu(I,j) * av_rem_u(I,j) > 0.0) &
          bt_rem_u(I,j) = G%mask2dCu(I,j) * CS%drag_amp * (av_rem_u(I,j)**Instep)
      enddo ; enddo
      do J=js-1,je ; do i=is-1,ie+1
        bt_rem_v(i,J) = 0.0
        if (G%mask2dCv(i,J) * av_rem_v(i,J) > 0.0) &
          bt_rem_v(i,J) = G%mask2dCv(i,J) * CS%drag_amp * (av_rem_v(i,J)**Instep)
      enddo ; enddo
    endif
  endif

  ! Zero out the arrays for various time-averaged quantities.
  if (find_etaav) then ; do j=jsvf-1,jevf+1 ; do i=isvf-1,ievf+1
    eta_sum(i,j) = 0.0 ; eta_wtd(i,j) = 0.0
  enddo ; enddo ; else ; do j=jsvf-1,jevf+1 ; do i=isvf-1,ievf+1
    eta_wtd(i,j) = 0.0
  enddo ; enddo ; endif
  do j=jsvf-1,jevf+1 ; do I=isvf-1,ievf
    ubt_sum(I,j) = 0.0 ; uhbt_sum(I,j) = 0.0
    PFu_bt_sum(I,j) = 0.0 ; Coru_bt_sum(I,j) = 0.0
    ubt_wtd(I,j) = 0.0 ; ubt_trans(I,j) = 0.0
  enddo ; enddo
  do J=jsvf-1,jevf ; do i=isvf-1,ievf+1
    vbt_sum(i,J) = 0.0 ; vhbt_sum(i,J) = 0.0
    PFv_bt_sum(i,J) = 0.0 ; Corv_bt_sum(i,J) = 0.0
    vbt_wtd(i,J) = 0.0 ; vbt_trans(I,j) = 0.0
  enddo ; enddo

  ! Set the mass source, after first initializing the halos to 0.
  do j=jsvf-1,jevf+1; do i=isvf-1,ievf+1 ; eta_src(i,j) = 0.0 ; enddo ; enddo
  if (CS%bound_BT_corr) then ; do j=js,je ; do i=is,ie
    if (abs(CS%eta_cor(i,j)) > dt*CS%eta_cor_bound(i,j)) &
      CS%eta_cor(i,j) = sign(dt*CS%eta_cor_bound(i,j),CS%eta_cor(i,j))
  enddo ; enddo ; endif
  do j=js,je ; do i=is,ie
    eta_src(i,j) = G%mask2dT(i,j) * (Instep * CS%eta_cor(i,j) + dtbt * CS%eta_source(i,j))
  enddo ; enddo

  if (CS%dynamic_psurf) then
    ice_is_rigid = (associated(fluxes%rigidity_ice_u) .and. &
                    associated(fluxes%rigidity_ice_v))
    H_min_dyn = GV%m_to_H * CS%Dmin_dyn_psurf
    if (ice_is_rigid .and. use_BT_cont) &
      call BT_cont_to_face_areas(BT_cont, Datu, Datv, G, MS, 0, .true.)
    if (ice_is_rigid) then ; do j=js,je ; do i=is,ie
      ! First determine the maximum stable value for dyn_coef_eta.

      !   This estimate of the maximum stable time step is pretty accurate for
      ! gravity waves, but it is a conservative estimate since it ignores the
      ! stabilizing effect of the bottom drag.
      Idt_max2 = 0.5 * (dgeo_de * (1.0 + 2.0*bebt)) * (G%IareaT(i,j) * &
            ((gtot_E(i,j) * (Datu(I,j)*G%IdxCu(I,j)) + &
              gtot_W(i,j) * (Datu(I-1,j)*G%IdxCu(I-1,j))) + &
             (gtot_N(i,j) * (Datv(i,J)*G%IdyCv(i,J)) + &
              gtot_S(i,j) * (Datv(i,J-1)*G%IdyCv(i,J-1)))) + &
            ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
             (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2)))
      H_eff_dx2 = max(H_min_dyn * (G%IdxT(i,j)**2 + G%IdyT(i,j)**2), &
                      G%IareaT(i,j) * &
                        ((Datu(I,j)*G%IdxCu(I,j) + Datu(I-1,j)*G%IdxCu(I-1,j)) + &
                         (Datv(i,J)*G%IdyCv(i,J) + Datv(i,J-1)*G%IdyCv(i,J-1)) ) )
      dyn_coef_max = CS%const_dyn_psurf * max(0.0, 1.0 - dtbt**2 * Idt_max2) / &
                     (dtbt**2 * H_eff_dx2)

      ! ice_strength has units of m s-2. rigidity_ice_[uv] has units of m3 s-1.
      ice_strength = ((fluxes%rigidity_ice_u(I,j) + fluxes%rigidity_ice_u(I-1,j)) + &
                      (fluxes%rigidity_ice_v(i,J) + fluxes%rigidity_ice_v(i,J-1))) / &
                      (CS%ice_strength_length**2 * dtbt)

      ! Units of dyn_coef: m2 s-2 H-1
      dyn_coef_eta(I,j) = min(dyn_coef_max, ice_strength * GV%H_to_m)
    enddo ; enddo ; endif
  endif

  if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
  if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
  if (nonblock_setup) then
    if (CS%dynamic_psurf) pid_dyn_coef_eta = &
      pass_var_start(dyn_coef_eta, CS%BT_Domain, complete=.false.)
    if (interp_eta_PF) then
      pid_eta_PF_1 = pass_var_start(eta_PF_1, CS%BT_Domain, complete=.false.)
      pid_d_eta_PF = pass_var_start(d_eta_PF, CS%BT_Domain, complete=.false.)
    else
      if ((G%isd > CS%isdw) .or. (G%jsd > CS%jsdw)) &
        pid_eta_PF = pass_var_start(eta_PF, CS%BT_Domain, complete=.false.)
    endif
    pid_eta_src = pass_var_start(eta_src, CS%BT_Domain)

    ! The following halo update is not needed without wide halos.  RWH
    if (ievf > ie) &
      pid_bt_rem_u = pass_vector_start(bt_rem_u, bt_rem_v, CS%BT_Domain, &
                                       To_All+Scalar_Pair)
  else
    if (interp_eta_PF) then
      call pass_var(eta_PF_1, CS%BT_Domain, complete=.false.)
      call pass_var(d_eta_PF, CS%BT_Domain, complete=.false.)
    else
      !   eta_PF_in had correct values in its halos, so only update eta_PF with
      ! extra-wide halo arrays.  This could have started almost immediately.
      if ((G%isd > CS%isdw) .or. (G%jsd > CS%jsdw)) &
        call pass_var(eta_PF, CS%BT_Domain, complete=.false.)
    endif
    if (CS%dynamic_psurf) call pass_var(dyn_coef_eta, CS%BT_Domain, complete=.false.)
    ! The following halo update is not needed without wide halos.  RWH
    if (ievf > ie) &
      call pass_vector(bt_rem_u, bt_rem_v, CS%BT_Domain, To_All+Scalar_Pair)
    if (.not.use_BT_cont) &
      call pass_vector(Datu, Datv, CS%BT_Domain, To_All+Scalar_Pair)
    if (((G%isd > CS%isdw) .or. (G%jsd > CS%jsdw)) .or. (Isq <= is-1) .or. (Jsq <= js-1)) &
      call pass_vector(BT_force_u, BT_force_v, CS%BT_Domain, complete=.false.)
    if (G%nonblocking_updates) then ! Passing needs to be completed now.
      call pass_var(eta_src, CS%BT_Domain, complete=.true.)
      if (add_uh0) call pass_vector(uhbt0, vhbt0, CS%BT_Domain, complete=.false.)
      call pass_vector(Cor_ref_u, Cor_ref_v, CS%BT_Domain, complete=.true.)
    else
      call pass_var(eta_src, CS%BT_Domain, complete=.false.)
      if (add_uh0) call pass_vector(uhbt0, vhbt0, CS%BT_Domain, complete=.false.)
      call pass_vector(Cor_ref_u, Cor_ref_v, CS%BT_Domain, complete=.false.)
    endif
  endif
  if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
  if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)

  ! Complete all of the outstanding halo updates.
  if (nonblock_setup) then
    if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
    if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)

    if (.not.use_BT_cont) &  !### IS THIS OK HERE?
      call pass_vector_complete(pid_Datu, Datu, Datv, CS%BT_Domain, To_All+Scalar_Pair)
    ! The following halo update is not needed without wide halos.  RWH
    if (((G%isd > CS%isdw) .or. (G%jsd > CS%jsdw)) .or. (Isq <= is-1) .or. (Jsq <= js-1)) &
      call pass_vector_complete(pid_BT_force_u, BT_force_u, BT_force_v, CS%BT_Domain)
    if (add_uh0) call pass_vector_complete(pid_uhbt0, uhbt0, vhbt0, CS%BT_Domain)
    call pass_vector_complete(pid_Cor_ref, Cor_ref_u, Cor_ref_v, CS%BT_Domain)

    if (CS%dynamic_psurf) &
      call pass_var_complete(pid_dyn_coef_eta, dyn_coef_eta, CS%BT_Domain)
    if (interp_eta_PF) then
      call pass_var_complete(pid_eta_PF_1, eta_PF_1, CS%BT_Domain)
      call pass_var_complete(pid_d_eta_PF, d_eta_PF, CS%BT_Domain)
    else
      if ((G%isd > CS%isdw) .or. (G%jsd > CS%jsdw)) &
        call pass_var_complete(pid_eta_PF, eta_PF, CS%BT_Domain)
    endif
    call pass_var_complete(pid_eta_src, eta_src, CS%BT_Domain)

    if (ievf > ie) &
      call pass_vector_complete(pid_bt_rem_u, bt_rem_u, bt_rem_v, CS%BT_Domain,&
      To_All+Scalar_Pair)

    if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
    if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)
  endif

  if (CS%debug) then
    call uchksum(uhbt, "BT uhbt",CS%debug_BT_HI,haloshift=0)
    call vchksum(vhbt, "BT vhbt",CS%debug_BT_HI,haloshift=0)
    call uchksum(ubt, "BT Initial ubt",CS%debug_BT_HI,haloshift=0)
    call vchksum(vbt, "BT Initial vbt",CS%debug_BT_HI,haloshift=0)
    call hchksum(GV%H_to_kg_m2*eta, "BT Initial eta",CS%debug_BT_HI,haloshift=0)
    call uchksum(BT_force_u, "BT BT_force_u",CS%debug_BT_HI,haloshift=0)
    call vchksum(BT_force_v, "BT BT_force_v",CS%debug_BT_HI,haloshift=0)
    if (interp_eta_PF) then
      call hchksum(GV%H_to_kg_m2*eta_PF_1, "BT eta_PF_1",CS%debug_BT_HI,haloshift=0)
      call hchksum(GV%H_to_kg_m2*d_eta_PF, "BT d_eta_PF",CS%debug_BT_HI,haloshift=0)
    else
      call hchksum(GV%H_to_kg_m2*eta_PF, "BT eta_PF",CS%debug_BT_HI,haloshift=0)
      call hchksum(GV%H_to_kg_m2*eta_PF_in, "BT eta_PF_in",G%HI,haloshift=0)
    endif
    call uchksum(Cor_ref_u, "BT Cor_ref_u",CS%debug_BT_HI,haloshift=0)
    call vchksum(Cor_ref_v, "BT Cor_ref_v",CS%debug_BT_HI,haloshift=0)
    call uchksum(uhbt0, "BT uhbt0",CS%debug_BT_HI,haloshift=0)
    call vchksum(vhbt0, "BT vhbt0",CS%debug_BT_HI,haloshift=0)
    if (.not. use_BT_cont) then
      call uchksum(GV%H_to_m*Datu, "BT Datu",CS%debug_BT_HI,haloshift=1)
      call vchksum(GV%H_to_m*Datv, "BT Datv",CS%debug_BT_HI,haloshift=1)
    endif
    call uchksum(wt_u, "BT wt_u",G%HI,haloshift=1)
    call vchksum(wt_v, "BT wt_v",G%HI,haloshift=1)
    call uchksum(CS%frhatu, "BT frhatu",G%HI,haloshift=1)
    call vchksum(CS%frhatv, "BT frhatv",G%HI,haloshift=1)
    call uchksum(bc_accel_u, "BT bc_accel_u",G%HI,haloshift=0)
    call vchksum(bc_accel_v, "BT bc_accel_v",G%HI,haloshift=0)
    call uchksum(CS%IDatu, "BT IDatu",G%HI,haloshift=0)
    call vchksum(CS%IDatv, "BT IDatv",G%HI,haloshift=0)
    if (use_visc_rem) then
      call uchksum(visc_rem_u, "BT visc_rem_u",G%HI,haloshift=1)
      call vchksum(visc_rem_v, "BT visc_rem_v",G%HI,haloshift=1)
    endif
  endif

  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_eta_st > 0) call post_data(CS%id_eta_st, eta(isd:ied,jsd:jed), CS%diag)
    if (CS%id_ubt_st > 0) call post_data(CS%id_ubt_st, ubt(IsdB:IedB,jsd:jed), CS%diag)
    if (CS%id_vbt_st > 0) call post_data(CS%id_vbt_st, vbt(isd:ied,JsdB:JedB), CS%diag)
  endif

  if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
  if (id_clock_calc > 0) call cpu_clock_begin(id_clock_calc)

  if (project_velocity) then ; eta_PF_BT => eta ; else ; eta_PF_BT => eta_pred ; endif

  if (CS%dt_bt_filter >= 0.0) then
    dt_filt = 0.5 * max(0.0, min(CS%dt_bt_filter, 2.0*dt))
  else
    dt_filt = 0.5 * max(0.0, dt * min(-CS%dt_bt_filter, 2.0))
  endif
  nfilter = ceiling(dt_filt / dtbt)

  ! Set up the normalized weights for the filtered velocity.
  sum_wt_vel = 0.0 ; sum_wt_eta = 0.0 ; sum_wt_accel = 0.0 ; sum_wt_trans = 0.0
  allocate(wt_vel(nstep+nfilter)) ; allocate(wt_eta(nstep+nfilter))
  allocate(wt_trans(nstep+nfilter+1)) ; allocate(wt_accel(nstep+nfilter+1))
  allocate(wt_accel2(nstep+nfilter+1))
  do n=1,nstep+nfilter
    ! Modify this to use a different filter...
    if ( (n==nstep) .or. (dt_filt - abs(n-nstep)*dtbt >= 0.0)) then
      wt_vel(n) = 1.0  ; wt_eta(n) = 1.0
    elseif (dtbt + dt_filt - abs(n-nstep)*dtbt > 0.0) then
      wt_vel(n) = 1.0 + (dt_filt / dtbt) - abs(n-nstep) ; wt_eta(n) = wt_vel(n)
    else
      wt_vel(n) = 0.0  ; wt_eta(n) = 0.0
    endif
!###    if (n < nstep-nfilter) then ; wt_vel(n) = 0.0 ; else ; wt_vel(n) = 1.0 ; endif
!###    if (n < nstep-nfilter) then ; wt_eta(n) = 0.0 ; else ; wt_eta(n) = 1.0 ; endif

    ! The rest should not be changed.
    sum_wt_vel = sum_wt_vel + wt_vel(n) ; sum_wt_eta = sum_wt_eta + wt_eta(n)
  enddo
  wt_trans(nstep+nfilter+1) = 0.0 ; wt_accel(nstep+nfilter+1) = 0.0
  do n=nstep+nfilter,1,-1
    wt_trans(n) = wt_trans(n+1) + wt_eta(n)
    wt_accel(n) = wt_accel(n+1) + wt_vel(n)
    sum_wt_accel = sum_wt_accel + wt_accel(n) ; sum_wt_trans = sum_wt_trans + wt_trans(n)
  enddo
  ! Normalize the weights.
  I_sum_wt_vel = 1.0 / sum_wt_vel ; I_sum_wt_accel = 1.0 / sum_wt_accel
  I_sum_wt_eta = 1.0 / sum_wt_eta ; I_sum_wt_trans = 1.0 / sum_wt_trans
  do n=1,nstep+nfilter
    wt_vel(n) = wt_vel(n) * I_sum_wt_vel
    wt_accel2(n) = wt_accel(n)
    wt_accel(n) = wt_accel(n) * I_sum_wt_accel
    wt_eta(n) = wt_eta(n) * I_sum_wt_eta
!    wt_trans(n) = wt_trans(n) * I_sum_wt_trans
  enddo


  sum_wt_vel = 0.0 ; sum_wt_eta = 0.0 ; sum_wt_accel = 0.0 ; sum_wt_trans = 0.0

  ! The following loop contains all of the time steps.
  isv=is ; iev=ie ; jsv=js ; jev=je
  do n=1,nstep+nfilter

    sum_wt_vel = sum_wt_vel + wt_vel(n)
    sum_wt_eta = sum_wt_eta + wt_eta(n)
    sum_wt_accel = sum_wt_accel + wt_accel2(n)
    sum_wt_trans = sum_wt_trans + wt_trans(n)

    if (CS%clip_velocity) then
      do j=jsv,jev ; do I=isv-1,iev
        if (abs(ubt(I,j)) > CS%maxvel) then
          ! Add some error handling later.
          ubt(I,j) = SIGN(CS%maxvel,ubt(I,j))
        endif
      enddo ; enddo
      do J=jsv-1,jev ; do i=isv,iev
        if (abs(vbt(i,J)) > CS%maxvel) then
          ! Add some error handling later.
          vbt(i,J) = SIGN(CS%maxvel,vbt(i,J))
        endif
      enddo ; enddo
    endif

    if ((iev - stencil < ie) .or. (jev - stencil < je)) then
      if (id_clock_calc > 0) call cpu_clock_end(id_clock_calc)
      if (id_clock_pass_step > 0) call cpu_clock_begin(id_clock_pass_step)
      if (G%nonblocking_updates) then
        pid_ubt = pass_vector_start(ubt, vbt, CS%BT_Domain)
        pid_eta = pass_var_start(eta, CS%BT_Domain)
        call pass_vector_complete(pid_ubt, ubt, vbt, CS%BT_Domain)
        call pass_var_complete(pid_eta, eta, CS%BT_Domain)
      else
        call pass_var(eta, CS%BT_Domain)
        call pass_vector(ubt, vbt, CS%BT_Domain)
      endif
      isv = isvf ; iev = ievf ; jsv = jsvf ; jev = jevf
      if (id_clock_pass_step > 0) call cpu_clock_end(id_clock_pass_step)
      if (id_clock_calc > 0) call cpu_clock_begin(id_clock_calc)
    else
      isv = isv+stencil ; iev = iev-stencil
      jsv = jsv+stencil ; jev = jev-stencil
    endif

    if ((.not.use_BT_cont) .and. CS%Nonlinear_continuity .and. &
        (CS%Nonlin_cont_update_period > 0)) then
      if ((n>1) .and. (mod(n-1,CS%Nonlin_cont_update_period) == 0)) &
        call find_face_areas(Datu, Datv, G, GV, CS, MS, CS%rescale_D_bt, eta, 1+iev-ie)
    endif

    if (CS%dynamic_psurf .or. .not.project_velocity) then
      if (use_BT_cont) then
!$OMP parallel do default(none) shared(isv,iev,jsv,jev,uhbt,ubt,BTCL_u,uhbt0)
        do j=jsv-1,jev+1 ; do I=isv-2,iev+1
          uhbt(I,j) = find_uhbt(ubt(I,j),BTCL_u(I,j)) + uhbt0(I,j)
        enddo ; enddo
        do J=jsv-2,jev+1 ; do i=isv-1,iev+1
          vhbt(i,J) = find_vhbt(vbt(i,J),BTCL_v(i,J)) + vhbt0(i,J)
        enddo ; enddo
        do j=jsv-1,jev+1 ; do i=isv-1,iev+1
          eta_pred(i,j) = (eta(i,j) + eta_src(i,j)) + (dtbt * CS%IareaT(i,j)) * &
                     ((uhbt(I-1,j) - uhbt(I,j)) + (vhbt(i,J-1) - vhbt(i,J)))
        enddo ; enddo
      else
!$OMP parallel do default(none) shared(isv,iev,jsv,jev,eta_pred,eta,eta_src,dtbt,&
!$OMP                                  CS,Datu,ubt,uhbt0,Datv,vhbt0,vbt)
        do j=jsv-1,jev+1 ; do i=isv-1,iev+1
          eta_pred(i,j) = (eta(i,j) + eta_src(i,j)) + (dtbt * CS%IareaT(i,j)) * &
              (((Datu(I-1,j)*ubt(I-1,j) + uhbt0(I-1,j)) - &
                (Datu(I,j)*ubt(I,j) + uhbt0(I,j))) + &
               ((Datv(i,J-1)*vbt(i,J-1) + vhbt0(i,J-1)) - &
                (Datv(i,J)*vbt(i,J) + vhbt0(i,J))))
        enddo ; enddo
      endif

      if (CS%dynamic_psurf) then ; do j=jsv-1,jev+1 ; do i=isv-1,iev+1
        p_surf_dyn(i,j) = dyn_coef_eta(I,j) * (eta_pred(i,j) - eta(i,j))
      enddo ; enddo ; endif
    endif

    ! Recall that just outside the do n loop, there is code like...
    !  eta_PF_BT => eta_pred ; if (project_velocity) eta_PF_BT => eta

    if (find_etaav) then ; do j=js,je ; do i=is,ie
      eta_sum(i,j) = eta_sum(i,j) + wt_accel2(n) * eta_PF_BT(i,j)
    enddo ; enddo ; endif

    if (interp_eta_PF) then
      wt_end = n*Instep  ! This could be (n-0.5)*Instep.
      do j=jsv-1,jev+1 ; do i=isv-1,iev+1
        eta_PF(i,j) = eta_PF_1(i,j) + wt_end*d_eta_PF(i,j)
      enddo ; enddo
    endif

    if (apply_OBC_flather) then
!$OMP parallel do default(none) shared(isv,iev,jsv,jev,ubt_old,ubt)
      do j=jsv,jev ; do I=isv-2,iev+1
        ubt_old(I,j) = ubt(I,j)
      enddo; enddo
!$OMP parallel do default(none) shared(isv,iev,jsv,jev,vbt_old,vbt)
      do J=jsv-2,jev+1 ; do i=isv,iev
        vbt_old(i,J) = vbt(i,J)
      enddo ; enddo
    endif

!$OMP parallel default(none) shared(isv,iev,jsv,jev,G,amer,ubt,cmer,bmer,dmer,eta_PF_BT, &
!$OMP                               eta_PF,gtot_N,gtot_S,dgeo_de,CS,p_surf_dyn,n,        &
!$OMP                               v_accel_bt,wt_accel,PFv_bt_sum,wt_accel2,            &
!$OMP                               Corv_bt_sum,BT_OBC,vbt,bt_rem_v,BT_force_v,vhbt,     &
!$OMP                               Cor_ref_v,find_PF,find_Cor,apply_v_OBCs,dtbt,        &
!$OMP                               project_velocity,be_proj,bebt,use_BT_cont,BTCL_v,    &
!$OMP                               vhbt0,Datv,vbt_sum,wt_trans,vhbt_sum,vbt_wtd,wt_vel, &
!$OMP                               azon,bzon,czon,dzon,Cor_ref_u,gtot_E,gtot_W,OBC,     &
!$OMP                               u_accel_bt,PFu_bt_sum,Coru_bt_sum,apply_u_OBCs,      &
!$OMP                               bt_rem_u,BT_force_u,uhbt,BTCL_u,uhbt0,Datu,ubt_sum,  &
!$OMP                               uhbt_sum,ubt_wtd)                                    &
!$OMP                       private(Cor, gradP, vel_prev, vel_trans )
    if (MOD(n+G%first_direction,2)==1) then
      ! On odd-steps, update v first.
!$OMP do
      do J=jsv-1,jev ; do i=isv-1,iev+1
        Cor = -1.0*((amer(I-1,j) * ubt(I-1,j) + cmer(I,j+1) * ubt(I,j+1)) + &
               (bmer(I,j) * ubt(I,j) + dmer(I-1,j+1) * ubt(I-1,j+1))) - Cor_ref_v(i,J)
        gradP = ((eta_PF_BT(i,j)-eta_PF(i,j))*gtot_N(i,j) - &
                 (eta_PF_BT(i,j+1)-eta_PF(i,j+1))*gtot_S(i,j+1)) * &
                dgeo_de * CS%IdyCv(i,J)
        if (CS%dynamic_psurf) &
          gradP = gradP + (p_surf_dyn(i,j) - p_surf_dyn(i,j+1)) * CS%IdyCv(i,J)
        v_accel_bt(i,J) = v_accel_bt(i,J) + wt_accel(n) * (Cor + gradP)
        if (find_PF)  PFv_bt_sum(i,J)  = PFv_bt_sum(i,J) + wt_accel2(n) * gradP
        if (find_Cor) Corv_bt_sum(i,J) = Corv_bt_sum(i,J) + wt_accel2(n) * Cor

        if (apply_v_OBCs) then ; if (OBC%OBC_segment_v(i,J) /= OBC_NONE) cycle ; endif

        vel_prev = vbt(i,J)
        vbt(i,J) = bt_rem_v(i,J) * (vbt(i,J) + &
                        dtbt * ((BT_force_v(i,J) + Cor) + gradP))
        if (project_velocity) then
          vel_trans = (1.0 + be_proj)*vbt(i,J) - be_proj*vel_prev
        else
          vel_trans = (1.0-bebt)*vel_prev + bebt*vbt(i,J)
        endif
        if (use_BT_cont) then
          vhbt(i,J) = find_vhbt(vel_trans,BTCL_v(i,J)) + vhbt0(i,J)
        else
          vhbt(i,J) = Datv(i,J)*vel_trans + vhbt0(i,J)
        endif

        vbt_sum(i,J) = vbt_sum(i,J) + wt_trans(n) * vel_trans
        vhbt_sum(i,J) = vhbt_sum(i,J) + wt_trans(n) * vhbt(i,J)
        vbt_wtd(i,J) = vbt_wtd(i,J) + wt_vel(n) * vbt(i,J)
      enddo ; enddo

!$OMP do
      do j=jsv,jev ; do I=isv-1,iev
        Cor = ((azon(I,j) * vbt(i+1,J) + czon(I,j) * vbt(i,J-1)) + &
               (bzon(I,j) * vbt(i,J) + dzon(I,j) * vbt(i+1,J-1))) - Cor_ref_u(I,j)
        gradP = ((eta_PF_BT(i,j)-eta_PF(i,j))*gtot_E(i,j) - &
                 (eta_PF_BT(i+1,j)-eta_PF(i+1,j))*gtot_W(i+1,j)) * &
                dgeo_de * CS%IdxCu(I,j)
        if (CS%dynamic_psurf) &
          gradP = gradP + (p_surf_dyn(i,j) - p_surf_dyn(i+1,j)) * CS%IdxCu(I,j)
        u_accel_bt(I,j) = u_accel_bt(I,j) + wt_accel(n) * (Cor + gradP)
        if (find_PF)  PFu_bt_sum(I,j)  = PFu_bt_sum(I,j) + wt_accel2(n) * gradP
        if (find_Cor) Coru_bt_sum(I,j) = Coru_bt_sum(I,j) + wt_accel2(n) * Cor

        if (apply_u_OBCs) then ; if (OBC%OBC_segment_u(I,j) /= OBC_NONE) cycle ; endif

        vel_prev = ubt(I,j)
        ubt(I,j) = bt_rem_u(I,j) * (ubt(I,j) + &
                        dtbt * ((BT_force_u(I,j) + Cor) + gradP))
        if (project_velocity) then
          vel_trans = (1.0 + be_proj)*ubt(I,j) - be_proj*vel_prev
        else
          vel_trans = (1.0-bebt)*vel_prev + bebt*ubt(I,j)
        endif
        if (use_BT_cont) then
          uhbt(I,j) = find_uhbt(vel_trans, BTCL_u(I,j)) + uhbt0(I,j)
        else
          uhbt(I,j) = Datu(I,j)*vel_trans + uhbt0(I,j)
        endif

        ubt_sum(I,j) = ubt_sum(I,j) + wt_trans(n) * vel_trans
        uhbt_sum(I,j) = uhbt_sum(I,j) + wt_trans(n) * uhbt(I,j)
        ubt_wtd(I,j) = ubt_wtd(I,j) + wt_vel(n) * ubt(I,j)
      enddo ; enddo

    else
      ! On even steps, update u first.
!$OMP do
      do j=jsv-1,jev+1 ; do I=isv-1,iev
        Cor = ((azon(I,j) * vbt(i+1,J) + czon(I,j) * vbt(i,J-1)) + &
               (bzon(I,j) * vbt(i,J) +  dzon(I,j) * vbt(i+1,J-1))) - Cor_ref_u(I,j)
        gradP = ((eta_PF_BT(i,j)-eta_PF(i,j))*gtot_E(i,j) - &
                 (eta_PF_BT(i+1,j)-eta_PF(i+1,j))*gtot_W(i+1,j)) * &
                dgeo_de * CS%IdxCu(I,j)
        if (CS%dynamic_psurf) &
          gradP = gradP + (p_surf_dyn(i,j) - p_surf_dyn(i+1,j)) * CS%IdxCu(I,j)
        u_accel_bt(I,j) = u_accel_bt(I,j) + wt_accel(n) * (Cor + gradP)
        if (find_PF)  PFu_bt_sum(I,j)  = PFu_bt_sum(I,j) + wt_accel2(n) * gradP
        if (find_Cor) Coru_bt_sum(I,j) = Coru_bt_sum(I,j) + wt_accel2(n) * Cor

        if (apply_u_OBCs) then ; if (OBC%OBC_segment_u(I,j) /= OBC_NONE) cycle ; endif

        vel_prev = ubt(I,j)
        ubt(I,j) = bt_rem_u(I,j) * (ubt(I,j) + &
                        dtbt * ((BT_force_u(I,j) + Cor) + gradP))
        if (project_velocity) then
          vel_trans = (1.0 + be_proj)*ubt(I,j) - be_proj*vel_prev
        else
          vel_trans = (1.0-bebt)*vel_prev + bebt*ubt(I,j)
        endif
        if (use_BT_cont) then
          uhbt(I,j) = find_uhbt(vel_trans,BTCL_u(I,j)) + uhbt0(I,j)
        else
          uhbt(I,j) = Datu(I,j)*vel_trans + uhbt0(I,j)
        endif

        ubt_sum(I,j) = ubt_sum(I,j) + wt_trans(n) * vel_trans
        uhbt_sum(I,j) = uhbt_sum(I,j) + wt_trans(n) * uhbt(I,j)
        ubt_wtd(I,j) = ubt_wtd(I,j) + wt_vel(n) * ubt(I,j)
      enddo ; enddo

!$OMP do
      do J=jsv-1,jev ; do i=isv,iev
        Cor = -1.0*((amer(I-1,j) * ubt(I-1,j) + bmer(I,j) * ubt(I,j)) + &
                (cmer(I,j+1) * ubt(I,j+1) + dmer(I-1,j+1) * ubt(I-1,j+1))) - Cor_ref_v(i,J)
        gradP = ((eta_PF_BT(i,j)-eta_PF(i,j))*gtot_N(i,j) - &
                 (eta_PF_BT(i,j+1)-eta_PF(i,j+1))*gtot_S(i,j+1)) * &
                dgeo_de * CS%IdyCv(i,J)
        if (CS%dynamic_psurf) &
          gradP = gradP + (p_surf_dyn(i,j) - p_surf_dyn(i,j+1)) * CS%IdyCv(i,J)
        v_accel_bt(I,j) = v_accel_bt(I,j) + wt_accel(n) * (Cor + gradP)
        if (find_PF)  PFv_bt_sum(i,J)  = PFv_bt_sum(i,J) + wt_accel2(n) * gradP
        if (find_Cor) Corv_bt_sum(i,J) = Corv_bt_sum(i,J) + wt_accel2(n) * Cor

        if (apply_v_OBCs) then ; if (OBC%OBC_segment_v(i,J) /= OBC_NONE) cycle ; endif

        vel_prev = vbt(i,J)
        vbt(i,J) = bt_rem_v(i,J) * (vbt(i,J) + &
                        dtbt * ((BT_force_v(i,J) + Cor) + gradP))
        if (project_velocity) then
          vel_trans = (1.0 + be_proj)*vbt(i,J) - be_proj*vel_prev
        else
          vel_trans = (1.0-bebt)*vel_prev + bebt*vbt(i,J)
        endif
         if (use_BT_cont) then
          vhbt(i,J) = find_vhbt(vel_trans,BTCL_v(i,J)) + vhbt0(i,J)
        else
          vhbt(i,J) = Datv(i,J)*vel_trans + vhbt0(i,J)
        endif

        vbt_sum(i,J) = vbt_sum(i,J) + wt_trans(n) * vel_trans
        vhbt_sum(i,J) = vhbt_sum(i,J) + wt_trans(n) * vhbt(i,J)
        vbt_wtd(i,J) = vbt_wtd(i,J) + wt_vel(n) * vbt(i,J)
      enddo ; enddo
    endif
!$OMP end parallel

    if (apply_OBCs) then
      call apply_velocity_OBCs(OBC, ubt, vbt, uhbt, vhbt, &
           ubt_trans, vbt_trans, eta, ubt_old, vbt_old, BT_OBC, &
           G, MS, iev-ie, dtbt, bebt, use_BT_cont, Datu, Datv, BTCL_u, BTCL_v, &
           uhbt0, vhbt0)
      if (apply_u_OBCs) then ; do j=js,je ; do I=is-1,ie
        if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
          ubt_sum(I,j) = ubt_sum(I,j) + wt_trans(n) * ubt_trans(I,j)
          uhbt_sum(I,j) = uhbt_sum(I,j) + wt_trans(n) * uhbt(I,j)
          ubt_wtd(I,j) = ubt_wtd(I,j) + wt_vel(n) * ubt(I,j)
        endif
      enddo ; enddo ; endif
      if (apply_v_OBCs) then ; do J=js-1,je ; do i=is,ie
        if (OBC%OBC_segment_v(i,J) /= OBC_NONE) then
          vbt_sum(i,J) = vbt_sum(i,J) + wt_trans(n) * vbt_trans(i,J)
          vhbt_sum(i,J) = vhbt_sum(i,J) + wt_trans(n) * vhbt(i,J)
          vbt_wtd(i,J) = vbt_wtd(i,J) + wt_vel(n) * vbt(i,J)
        endif
      enddo ; enddo ; endif
    endif

    if (CS%debug_bt) then
      call uchksum(uhbt, "BT uhbt just after OBC",CS%debug_BT_HI,haloshift=iev-ie)
      call vchksum(vhbt, "BT vhbt just after OBC",CS%debug_BT_HI,haloshift=iev-ie)
    endif

!$OMP parallel do default(none) shared(isv,iev,jsv,jev,n,eta,eta_src,dtbt,CS,uhbt,vhbt,eta_wtd,wt_eta)
    do j=jsv,jev ; do i=isv,iev
      eta(i,j) = (eta(i,j) + eta_src(i,j)) + (dtbt * CS%IareaT(i,j)) * &
                 ((uhbt(I-1,j) - uhbt(I,j)) + (vhbt(i,J-1) - vhbt(i,J)))
      eta_wtd(i,j) = eta_wtd(i,j) + eta(i,j) * wt_eta(n)
    enddo ; enddo
    if (apply_OBCs) call apply_eta_OBCs(OBC, eta, ubt_old, vbt_old, BT_OBC, &
                                        G, MS, iev-ie, dtbt)

    if (do_hifreq_output) then
      time_step_end = time_bt_start + set_time(int(floor(n*dtbt+0.5)))
      call enable_averaging(dtbt, time_step_end, CS%diag)
      if (CS%id_ubt_hifreq > 0) call post_data(CS%id_ubt_hifreq, ubt(IsdB:IedB,jsd:jed), CS%diag)
      if (CS%id_vbt_hifreq > 0) call post_data(CS%id_vbt_hifreq, vbt(isd:ied,JsdB:JedB), CS%diag)
      if (CS%id_eta_hifreq > 0) call post_data(CS%id_eta_hifreq, eta(isd:ied,jsd:jed), CS%diag)
      if (CS%id_uhbt_hifreq > 0) call post_data(CS%id_uhbt_hifreq, uhbt(IsdB:IedB,jsd:jed), CS%diag)
      if (CS%id_vhbt_hifreq > 0) call post_data(CS%id_vhbt_hifreq, vhbt(isd:ied,JsdB:JedB), CS%diag)
      if (CS%id_eta_pred_hifreq > 0) call post_data(CS%id_eta_pred_hifreq, eta_PF_BT(isd:ied,jsd:jed), CS%diag)
    endif

    if (CS%debug_bt) then
      write(mesg,'("BT step ",I4)') n
      call uchksum(ubt, trim(mesg)//" ubt",CS%debug_BT_HI,haloshift=iev-ie)
      call vchksum(vbt, trim(mesg)//" vbt",CS%debug_BT_HI,haloshift=iev-ie)
      call hchksum(GV%H_to_kg_m2*eta, trim(mesg)//" eta",CS%debug_BT_HI,haloshift=iev-ie)
    endif

  enddo ! end of do n=1,ntimestep
  if (id_clock_calc > 0) call cpu_clock_end(id_clock_calc)
  if (id_clock_calc_post > 0) call cpu_clock_begin(id_clock_calc_post)

  ! Reset the time information in the diag type.
  if (do_hifreq_output) call enable_averaging(time_int_in, time_end_in, CS%diag)

  I_sum_wt_vel = 1.0 / sum_wt_vel ; I_sum_wt_eta = 1.0 / sum_wt_eta
  I_sum_wt_accel = 1.0 / sum_wt_accel ; I_sum_wt_trans = 1.0 / sum_wt_trans

  if (find_etaav) then ; do j=js,je ; do i=is,ie
    etaav(i,j) = eta_sum(i,j) * I_sum_wt_accel
  enddo ; enddo ; endif
  do j=js-1,je+1 ; do i=is-1,ie+1 ; e_anom(i,j) = 0.0 ; enddo ; enddo
  if (interp_eta_PF) then
    do j=js,je ; do i=is,ie
      e_anom(i,j) = dgeo_de * (0.5 * (eta(i,j) + eta_in(i,j)) - &
                               (eta_PF_1(i,j) + 0.5*d_eta_PF(i,j)))
    enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      e_anom(i,j) = dgeo_de * (0.5 * (eta(i,j) + eta_in(i,j)) - eta_PF(i,j))
    enddo ; enddo
  endif

  ! It is possible that eta_out and eta_in are the same.
  do j=js,je ; do i=is,ie
    eta_out(i,j) = eta_wtd(i,j) * I_sum_wt_eta
  enddo ; enddo

  if (id_clock_calc_post > 0) call cpu_clock_end(id_clock_calc_post)
  if (id_clock_pass_post > 0) call cpu_clock_begin(id_clock_pass_post)
  if (G%nonblocking_updates) then
    pid_e_anom = pass_var_start(e_anom, G%Domain)
  else
    if (find_etaav) call pass_var(etaav, G%Domain, complete=.false.)
    call pass_var(e_anom, G%Domain)
  endif
  if (id_clock_pass_post > 0) call cpu_clock_end(id_clock_pass_post)
  if (id_clock_calc_post > 0) call cpu_clock_begin(id_clock_calc_post)

  do j=js,je ; do I=is-1,ie
    CS%ubtav(I,j) = ubt_sum(I,j) * I_sum_wt_trans
    uhbtav(I,j) = uhbt_sum(I,j) * I_sum_wt_trans
 !###   u_accel_bt(I,j) = u_accel_bt(I,j) * I_sum_wt_accel
    ubt_wtd(I,j) = ubt_wtd(I,j) * I_sum_wt_vel
  enddo ; enddo

  do J=js-1,je ; do i=is,ie
    CS%vbtav(i,J) = vbt_sum(i,J) * I_sum_wt_trans
    vhbtav(i,J) = vhbt_sum(i,J) * I_sum_wt_trans
 !###   v_accel_bt(i,J) = v_accel_bt(i,J)  * I_sum_wt_accel
    vbt_wtd(i,J) = vbt_wtd(i,J) * I_sum_wt_vel
  enddo ; enddo

  if (present(uhbt_out)) then
    uhbt_out(:,:) = 0.0
    if (use_BT_cont) then ; do j=js,je ; do I=is-1,ie
      uhbt_out(I,j) = find_uhbt(ubt_wtd(I,j), BTCL_u(I,j)) + uhbt0(I,j)
    enddo ; enddo ; else ; do j=js,je ; do I=is-1,ie
      uhbt_out(I,j) = ubt_wtd(I,j) * Datu(I,j) + uhbt0(I,j)
    enddo ; enddo ; endif

    vhbt_out(:,:) = 0.0
    if (use_BT_cont) then ; do J=js-1,je ; do i=is,ie
      vhbt_out(i,J) = find_vhbt(vbt_wtd(i,J), BTCL_v(i,J)) + vhbt0(i,J)
    enddo ; enddo ; else ; do J=js-1,je ; do i=is,ie
      vhbt_out(i,J) = vbt_wtd(i,J) * Datv(i,J) + vhbt0(i,J)
    enddo ; enddo ; endif
  endif

  if (id_clock_calc_post > 0) call cpu_clock_end(id_clock_calc_post)
  if (id_clock_pass_post > 0) call cpu_clock_begin(id_clock_pass_post)
  if (G%nonblocking_updates) then
    call pass_var_complete(pid_e_anom, e_anom, G%Domain)

    if (find_etaav) pid_etaav = pass_var_start(etaav, G%Domain)
    pid_uhbtav = pass_vector_start(uhbtav, vhbtav, G%Domain, complete=.false.)
    pid_ubtav = pass_vector_start(CS%ubtav, CS%vbtav, G%Domain)
  else
    call pass_vector(CS%ubtav, CS%vbtav, G%Domain, complete=.false.)
    call pass_vector(uhbtav, vhbtav, G%Domain)
  endif
  if (id_clock_pass_post > 0) call cpu_clock_end(id_clock_pass_post)
  if (id_clock_calc_post > 0) call cpu_clock_begin(id_clock_calc_post)

! Now calculate each layer's accelerations.
!$OMP parallel do default(none) shared(is,ie,js,je,nz,accel_layer_u,u_accel_bt,pbce,gtot_W, &
!$OMP                                  e_anom,gtot_E,CS,accel_layer_v,v_accel_bt,      &
!$OMP                                  gtot_S,gtot_N)
  do k=1,nz
    do j=js,je ; do I=is-1,ie
      accel_layer_u(I,j,k) = u_accel_bt(I,j) - &
           ((pbce(i+1,j,k) - gtot_W(i+1,j)) * e_anom(i+1,j) - &
            (pbce(i,j,k) - gtot_E(i,j)) * e_anom(i,j)) * CS%IdxCu(I,j)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      accel_layer_v(i,J,k) = v_accel_bt(i,J) - &
           ((pbce(i,j+1,k) - gtot_S(i,j+1))*e_anom(i,j+1) - &
            (pbce(i,j,k) - gtot_N(i,j))*e_anom(i,j)) * CS%IdyCv(i,J)
    enddo ; enddo
  enddo

  if (id_clock_calc_post > 0) call cpu_clock_end(id_clock_calc_post)

  ! Calculate diagnostic quantities.
  if (query_averaging_enabled(CS%diag)) then

    do j=js,je ; do I=is-1,ie ; CS%ubt_IC(I,j) = ubt_wtd(I,j) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; CS%vbt_IC(i,J) = vbt_wtd(i,J) ; enddo ; enddo
    if (present(uhbt_out)) then
      do j=js,je ; do I=is-1,ie ; CS%uhbt_IC(I,j) = uhbt_out(I,j) ; enddo ; enddo
      do J=js-1,je ; do i=is,ie ; CS%vhbt_IC(i,J) = vhbt_out(i,J) ; enddo ; enddo
    elseif (use_BT_cont) then
      do j=js,je ; do I=is-1,ie
        CS%uhbt_IC(I,j) = find_uhbt(ubt_wtd(I,j), BTCL_u(I,j)) + uhbt0(I,j)
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        CS%vhbt_IC(i,J) = find_vhbt(vbt_wtd(i,J), BTCL_v(i,J)) + vhbt0(i,J)
      enddo ; enddo
    else
      do j=js,je ; do I=is-1,ie
        CS%uhbt_IC(I,j) = ubt_wtd(I,j) * Datu(I,j) + uhbt0(I,j)
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        CS%vhbt_IC(i,J) = vbt_wtd(i,J) * Datv(i,J) + vhbt0(i,J)
      enddo ; enddo
    endif

!  Offer various barotropic terms for averaging.
    if (CS%id_PFu_bt > 0) then
      do j=js,je ; do I=is-1,ie
        PFu_bt_sum(I,j) = PFu_bt_sum(I,j) * I_sum_wt_accel
      enddo ; enddo
      call post_data(CS%id_PFu_bt, PFu_bt_sum(IsdB:IedB,jsd:jed), CS%diag)
    endif
    if (CS%id_PFv_bt > 0) then
      do J=js-1,je ; do i=is,ie
        PFv_bt_sum(i,J) = PFv_bt_sum(i,J) * I_sum_wt_accel
      enddo ; enddo
      call post_data(CS%id_PFv_bt, PFv_bt_sum(isd:ied,JsdB:JedB), CS%diag)
    endif
    if (CS%id_Coru_bt > 0) then
      do j=js,je ; do I=is-1,ie
        Coru_bt_sum(I,j) = Coru_bt_sum(I,j) * I_sum_wt_accel
      enddo ; enddo
      call post_data(CS%id_Coru_bt, Coru_bt_sum(IsdB:IedB,jsd:jed), CS%diag)
    endif
    if (CS%id_Corv_bt > 0) then
      do J=js-1,je ; do i=is,ie
        Corv_bt_sum(i,J) = Corv_bt_sum(i,J) * I_sum_wt_accel
      enddo ; enddo
      call post_data(CS%id_Corv_bt, Corv_bt_sum(isd:ied,JsdB:JedB), CS%diag)
    endif
    if (CS%id_ubtforce > 0) call post_data(CS%id_ubtforce, BT_force_u(IsdB:IedB,jsd:jed), CS%diag)
    if (CS%id_vbtforce > 0) call post_data(CS%id_vbtforce, BT_force_v(isd:ied,JsdB:JedB), CS%diag)
    if (CS%id_uaccel > 0) call post_data(CS%id_uaccel, u_accel_bt(IsdB:IedB,jsd:jed), CS%diag)
    if (CS%id_vaccel > 0) call post_data(CS%id_vaccel, v_accel_bt(isd:ied,JsdB:JedB), CS%diag)
    if (CS%id_Datu_res > 0) call post_data(CS%id_Datu_res, CS%Datu_res(IsdB:IedB,jsd:jed), CS%diag)
    if (CS%id_Datv_res > 0) call post_data(CS%id_Datv_res, CS%Datv_res(isd:ied,JsdB:JedB), CS%diag)

    if (CS%id_eta_cor > 0) call post_data(CS%id_eta_cor, CS%eta_cor, CS%diag)
    if (CS%id_eta_bt > 0) call post_data(CS%id_eta_bt, eta_out, CS%diag)
    if (CS%id_gtotn > 0) call post_data(CS%id_gtotn, gtot_N(isd:ied,jsd:jed), CS%diag)
    if (CS%id_gtots > 0) call post_data(CS%id_gtots, gtot_S(isd:ied,jsd:jed), CS%diag)
    if (CS%id_gtote > 0) call post_data(CS%id_gtote, gtot_E(isd:ied,jsd:jed), CS%diag)
    if (CS%id_gtotw > 0) call post_data(CS%id_gtotw, gtot_W(isd:ied,jsd:jed), CS%diag)
    if (CS%id_ubt > 0) call post_data(CS%id_ubt, ubt_wtd(IsdB:IedB,jsd:jed), CS%diag)
    if (CS%id_vbt > 0) call post_data(CS%id_vbt, vbt_wtd(isd:ied,JsdB:JedB), CS%diag)
    if (CS%id_ubtav > 0) call post_data(CS%id_ubtav, CS%ubtav, CS%diag)
    if (CS%id_vbtav > 0) call post_data(CS%id_vbtav, CS%vbtav, CS%diag)
    if (use_visc_rem) then
      if (CS%id_visc_rem_u > 0) call post_data(CS%id_visc_rem_u, visc_rem_u, CS%diag)
      if (CS%id_visc_rem_v > 0) call post_data(CS%id_visc_rem_v, visc_rem_v, CS%diag)
    endif

    if (CS%id_frhatu > 0) call post_data(CS%id_frhatu, CS%frhatu, CS%diag)
    if (CS%id_uhbt > 0) call post_data(CS%id_uhbt, uhbtav, CS%diag)
    if (CS%id_frhatv > 0) call post_data(CS%id_frhatv, CS%frhatv, CS%diag)
    if (CS%id_vhbt > 0) call post_data(CS%id_vhbt, vhbtav, CS%diag)

    if (CS%id_frhatu1 > 0) call post_data(CS%id_frhatu1, CS%frhatu1, CS%diag)
    if (CS%id_frhatv1 > 0) call post_data(CS%id_frhatv1, CS%frhatv1, CS%diag)
  else
    if (CS%id_frhatu1 > 0) CS%frhatu1(:,:,:) = CS%frhatu(:,:,:)
    if (CS%id_frhatv1 > 0) CS%frhatv1(:,:,:) = CS%frhatv(:,:,:)
  endif

  if (G%nonblocking_updates) then
    if (find_etaav) call pass_var_complete(pid_etaav, etaav, G%Domain)
    call pass_vector_complete(pid_uhbtav, uhbtav, vhbtav, G%Domain)
    call pass_vector_complete(pid_ubtav, CS%ubtav, CS%vbtav, G%Domain)
  endif

  if (apply_OBCs) call destroy_BT_OBC(BT_OBC)

end subroutine legacy_btstep

subroutine legacy_set_dtbt(G, GV, CS, eta, pbce, BT_cont, gtot_est, SSH_add)
  type(ocean_grid_type),                    intent(inout) :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  type(legacy_barotropic_CS),               pointer       :: CS
  real, dimension(SZI_(G),SZJ_(G)),         intent(in), optional :: eta
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in), optional :: pbce
  type(BT_cont_type),                       pointer,    optional :: BT_cont
  real,                                     intent(in), optional :: gtot_est
  real,                                     intent(in), optional :: SSH_add
! Arguments: G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 barotropic_init.
!  (in,opt)  eta - The barotropic free surface height anomaly or
!                  column mass anomaly, in m or kg m-2.
!  (in,opt)  pbce - The baroclinic pressure anomaly in each layer
!                   due to free surface height anomalies, in m2 H-1 s-2.
!  (in,opt)  BT_cont - A structure with elements that describe the effective
!                      open face areas as a function of barotropic flow.
!  (in,opt)  gtot_est - An estimate of the total gravitational acceleration,
!                       in m s-2.
!  (in,opt)  SSH_add - An additional contribution to SSH to provide a margin of
!                      error when calculating the external wave speed, in m.
  ! This subroutine automatically determines an optimal value for dtbt based
  ! on some state of the ocean.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    gtot_E, &     ! gtot_X is the effective total reduced gravity used to relate
    gtot_W, &     ! free surface height deviations to pressure forces (including
    gtot_N, &     ! GFS and baroclinic  contributions) in the barotropic momentum
    gtot_S        ! equations half a grid-point in the X-direction (X is N, S,
                  ! E, or W) from the thickness point. gtot_X has units of m2 H-1 s-2.
                  ! (See Hallberg, J Comp Phys 1997 for a discussion.)
  real, dimension(SZIBS_(G),SZJ_(G)) :: &
    Datu          ! Basin depth at u-velocity grid points times the y-grid
                  ! spacing, in m2.
  real, dimension(SZI_(G),SZJBS_(G)) :: &
    Datv          ! Basin depth at v-velocity grid points times the x-grid
                  ! spacing, in m2.
  real :: det_de  ! The partial derivative due to self-attraction and loading
                  ! of the reference geopotential with the sea surface height.
                  ! This is typically ~0.09 or less.
  real :: dgeo_de ! The constant of proportionality between geopotential and
                  ! sea surface height.  It is a nondimensional number of
                  ! order 1.  For stability, this may be made larger
                  ! than physical problem would suggest.
  real :: add_SSH ! An additional contribution to SSH to provide a margin of error
                  ! when calculating the external wave speed, in m.
  real :: min_max_dt2, Idt_max2, dtbt_max
  logical :: use_BT_cont
  type(memory_size_type) :: MS

  character(len=200) :: mesg
  integer :: i, j, k, is, ie, js, je, nz

  if (.not.associated(CS)) call MOM_error(FATAL, &
      "legacy_set_dtbt: Module MOM_barotropic must be initialized before it is used.")
  if (.not.CS%split) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  MS%isdw = G%isd ; MS%iedw = G%ied ; MS%jsdw = G%jsd ; MS%jedw = G%jed

  if (.not.(present(pbce) .or. present(gtot_est))) call MOM_error(FATAL, &
      "legacy_set_dtbt: Either pbce or gtot_est must be present.")

  add_SSH = 0.0 ; if (present(SSH_add)) add_SSH = SSH_add

  use_BT_cont = .false.
  if (present(BT_cont)) use_BT_cont = (associated(BT_cont))

  if (use_BT_cont) then
    call BT_cont_to_face_areas(BT_cont, Datu, Datv, G, MS, 0, .true.)
  elseif (CS%Nonlinear_continuity .and. present(eta)) then
    call find_face_areas(Datu, Datv, G, GV, CS, MS, eta=eta, halo=0)
  else
    call find_face_areas(Datu, Datv, G, GV, CS, MS, halo=0, add_max=add_SSH)
  endif

  det_de = 0.0
  if (CS%tides) call tidal_forcing_sensitivity(G, CS%tides_CSp, det_de)
  dgeo_de = 1.0 + max(0.0, det_de + CS%G_extra)
  if (present(pbce)) then
    do j=js,je ; do i=is,ie
      gtot_E(i,j) = 0.0 ; gtot_W(i,j) = 0.0
      gtot_N(i,j) = 0.0 ; gtot_S(i,j) = 0.0
    enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      gtot_E(i,j) = gtot_E(i,j) + pbce(i,j,k) * CS%frhatu(I,j,k)
      gtot_W(i,j) = gtot_W(i,j) + pbce(i,j,k) * CS%frhatu(I-1,j,k)
      gtot_N(i,j) = gtot_N(i,j) + pbce(i,j,k) * CS%frhatv(i,J,k)
      gtot_S(i,j) = gtot_S(i,j) + pbce(i,j,k) * CS%frhatv(i,J-1,k)
    enddo ; enddo ; enddo
  else
    do j=js,je ; do i=is,ie
      gtot_E(i,j) = gtot_est * GV%H_to_m ; gtot_W(i,j) = gtot_est * GV%H_to_m
      gtot_N(i,j) = gtot_est * GV%H_to_m ; gtot_S(i,j) = gtot_est * GV%H_to_m
    enddo ; enddo
  endif

  min_max_dt2 = 1.0e38  ! A huge number.
  do j=js,je ; do i=is,ie
    !   This is pretty accurate for gravity waves, but it is a conservative
    ! estimate since it ignores the stabilizing effect of the bottom drag.
    Idt_max2 = 0.5 * (1.0 + 2.0*CS%bebt) * (G%IareaT(i,j) * &
      ((gtot_E(i,j)*Datu(I,j)*G%IdxCu(I,j) + gtot_W(i,j)*Datu(I-1,j)*G%IdxCu(I-1,j)) + &
       (gtot_N(i,j)*Datv(i,J)*G%IdyCv(i,J) + gtot_S(i,j)*Datv(i,J-1)*G%IdyCv(i,J-1))) + &
      ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
       (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2)))
    if (Idt_max2 * min_max_dt2 > 1.0) min_max_dt2 = 1.0 / Idt_max2
  enddo ; enddo
  dtbt_max = sqrt(min_max_dt2 / dgeo_de)
  if (id_clock_sync > 0) call cpu_clock_begin(id_clock_sync)
  call min_across_PEs(dtbt_max)
  if (id_clock_sync > 0) call cpu_clock_end(id_clock_sync)

  CS%dtbt = CS%dtbt_fraction * dtbt_max
  CS%dtbt_max = dtbt_max
end subroutine legacy_set_dtbt

! The following 4 subroutines apply the open boundary conditions.
subroutine apply_velocity_OBCs(OBC, ubt, vbt, uhbt, vhbt, ubt_trans, vbt_trans, &
                               eta, ubt_old, vbt_old, BT_OBC, &
                               G, MS, halo, dtbt, bebt, use_BT_cont, Datu, Datv, &
                               BTCL_u, BTCL_v, uhbt0, vhbt0)
  type(ocean_OBC_type),              pointer       :: OBC
  type(ocean_grid_type),             intent(inout) :: G
  type(memory_size_type),            intent(in)    :: MS
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(inout) :: ubt, uhbt, ubt_trans
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(inout) :: vbt, vhbt, vbt_trans
  real, dimension(SZIW_(MS),SZJW_(MS)),   intent(in)    :: eta
  real, dimension(SZIBW_(MS),SZJW_(MS)),  intent(in)    :: ubt_old
  real, dimension(SZIW_(MS),SZJBW_(MS)),  intent(in)    :: vbt_old
  type(BT_OBC_type),                 intent(in)    :: BT_OBC
  integer,                           intent(in)    :: halo
  real,                              intent(in)    :: dtbt, bebt
  logical,                           intent(in)    :: use_BT_cont
  real, dimension(SZIBW_(MS),SZJW_(MS)),  intent(in)    :: Datu
  real, dimension(SZIW_(MS),SZJBW_(MS)),  intent(in)    :: Datv
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), intent(in) :: BTCL_u
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), intent(in) :: BTCL_v
  real, dimension(SZIBW_(MS),SZJW_(MS)),  intent(in)    :: uhbt0
  real, dimension(SZIW_(MS),SZJBW_(MS)),  intent(in)    :: vhbt0
!   This subroutine applies the open boundary conditions on barotropic
! velocities and mass transports, as developed by Mehmet Ilicak.

! Arguments: OBC - an associated pointer to an OBC type.
!  (inout)   ubt - the zonal barotropic velocity, in m s-1.
!  (inout)   vbt - the meridional barotropic velocity, in m s-1.
!  (inout)   uhbt - the zonal barotropic transport, in H m2 s-1.
!  (inout)   vhbt - the meridional barotropic transport, in H m2 s-1.
!  (inout)   ubt_trans - the zonal barotropic velocity used in transport, m s-1.
!  (inout)   vbt_trans - the meridional BT velocity used in transports, m s-1.
!  (in)      eta - The barotropic free surface height anomaly or
!                  column mass anomaly, in m or kg m-2.
!  (in)      ubt_old - The starting value of ubt in a barotropic step, m s-1.
!  (in)      vbt_old - The starting value of ubt in a barotropic step, m s-1.
!  (in)      BT_OBC - A structure with the private barotropic arrays related
!                     to the open boundary conditions, set by set_up_BT_OBC.
!  (in)      G - The ocean's grid structure.
!  (in)      MS - A type that describes the memory sizes of the argument arrays.
!  (in)      halo - The extra halo size to use here.
!  (in)      dtbt - The time step, in s.
!  (in)      bebt - The fractional weighting of the future velocity in
!                   determining the transport.
!  (in)      use_BT_cont - If true, use the BT_cont_types to calculate transports.
!  (in)      Datu - A fixed estimate of the face areas at u points.
!  (in)      Datv - A fixed estimate of the face areas at u points.
!  (in)      BCTL_u - Structures of information used for a dynamic estimate
!  (in)      BCTL_v - of the face areas at u- and v- points.
  real :: vel_prev    ! The previous velocity in m s-1.
  real :: vel_trans   ! The combination of the previous and current velocity
                      ! that does the mass transport, in m s-1.
  real :: H_u         ! The total thickness at the u-point, in m or kg m-2.
  real :: H_v         ! The total thickness at the v-point, in m or kg m-2.
  real :: cfl         ! The CFL number at the point in question, ND.
  real :: u_inlet
  real :: v_inlet
  real :: h_in
  integer :: i, j, is, ie, js, je
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo

  if (apply_u_OBCs) then
    do j=js,je ; do I=is-1,ie ; if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
      if (OBC%segment(OBC%OBC_segment_u(I,j))%specified) then
        uhbt(I,j) = BT_OBC%uhbt(I,j)
        ubt(I,j) = BT_OBC%ubt_outer(I,j)
        vel_trans = ubt(I,j)
      elseif (OBC%segment(OBC%OBC_segment_u(I,j))%Flather) then
        if (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_E) then
          cfl = dtbt * BT_OBC%Cg_u(I,j) * G%IdxCu(I,j)            ! CFL
          u_inlet = cfl*ubt_old(I-1,j) + (1.0-cfl)*ubt_old(I,j)  ! Valid for cfl<1
        !  h_in = 2.0*cfl*eta(i,j) + (1.0-2.0*cfl)*eta(i+1,j)    ! external
          h_in = eta(i,j) + (0.5-cfl)*(eta(i,j)-eta(i-1,j))      ! internal

          H_u = BT_OBC%H_u(I,j)
          vel_prev = ubt(I,j)
          ubt(I,j) = 0.5*((u_inlet + BT_OBC%ubt_outer(I,j)) + &
              (BT_OBC%Cg_u(I,j)/H_u) * (h_in-BT_OBC%eta_outer_u(I,j)))

          vel_trans = (1.0-bebt)*vel_prev + bebt*ubt(I,j)
        elseif (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_W) then
          cfl = dtbt * BT_OBC%Cg_u(I,j) * G%IdxCu(I,j)            ! CFL
          u_inlet = cfl*ubt_old(I+1,j) + (1.0-cfl)*ubt_old(I,j)  ! Valid for cfl<1
        !  h_in = 2.0*cfl*eta(i+1,j) + (1.0-2.0*cfl)*eta(i,j)    ! external
          h_in = eta(i+1,j) + (0.5-cfl)*(eta(i+1,j)-eta(i+2,j))  ! internal

          H_u = BT_OBC%H_u(I,j)
          vel_prev = ubt(I,j)
          ubt(I,j) = 0.5*((u_inlet+BT_OBC%ubt_outer(I,j)) + &
              (BT_OBC%Cg_u(I,j)/H_u) * (BT_OBC%eta_outer_u(I,j)-h_in))

          vel_trans = (1.0-bebt)*vel_prev + bebt*ubt(I,j)
        elseif (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_N) then
          if ((vbt(i,J-1)+vbt(i+1,J-1)) > 0.0) then
            ubt(I,j) = 2.0*ubt(I,j-1)-ubt(I,j-2)
          else
            ubt(I,j) = BT_OBC%ubt_outer(I,j)
          endif
          vel_trans = ubt(I,j)
        elseif (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_S) then
          if ((vbt(i,J)+vbt(i+1,J)) < 0.0) then
            ubt(I,j) = 2.0*ubt(I,j+1)-ubt(I,j+2)
          else
            ubt(I,j) = BT_OBC%ubt_outer(I,j)
          endif
          vel_trans = ubt(I,j)
        endif
      endif

      if (.not. OBC%segment(OBC%OBC_segment_u(I,j))%specified) then
        if (use_BT_cont) then
          uhbt(I,j) = find_uhbt(vel_trans,BTCL_u(I,j)) + uhbt0(I,j)
        else
          uhbt(I,j) = Datu(I,j)*vel_trans + uhbt0(I,j)
        endif
      endif

      ubt_trans(I,j) = vel_trans
    endif ; enddo ; enddo
  endif

  if (apply_v_OBCs) then
    do J=js-1,je ; do i=is,ie ; if (OBC%OBC_segment_v(i,J) /= OBC_NONE) then
      if (OBC%segment(OBC%OBC_segment_v(i,J))%specified) then
        vhbt(i,J) = BT_OBC%vhbt(i,J)
        vbt(i,J) = BT_OBC%vbt_outer(i,J)
        vel_trans = vbt(i,J)
      elseif (OBC%segment(OBC%OBC_segment_v(i,J))%Flather) then
        if (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_N) then
          cfl = dtbt * BT_OBC%Cg_v(i,J) * G%IdyCv(I,j)            ! CFL
          v_inlet = cfl*vbt_old(i,J-1) + (1.0-cfl)*vbt_old(i,J)  ! Valid for cfl<1
        !  h_in = 2.0*cfl*eta(i,j) + (1.0-2.0*cfl)*eta(i,j+1)    ! external
          h_in = eta(i,j) + (0.5-cfl)*(eta(i,j)-eta(i,j-1))      ! internal

          H_v = BT_OBC%H_v(i,J)
          vel_prev = vbt(i,J)
          vbt(i,J) = 0.5*((v_inlet+BT_OBC%vbt_outer(i,J)) + &
              (BT_OBC%Cg_v(i,J)/H_v) * (h_in-BT_OBC%eta_outer_v(i,J)))

          vel_trans = (1.0-bebt)*vel_prev + bebt*vbt(i,J)
        elseif (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_S) then
          cfl = dtbt * BT_OBC%Cg_v(i,J) * G%IdyCv(I,j)            ! CFL
          v_inlet = cfl*vbt_old(i,J+1) + (1.0-cfl)*vbt_old(i,J)  ! Valid for cfl <1
        !  h_in = 2.0*cfl*eta(i,j+1) + (1.0-2.0*cfl)*eta(i,j)    ! external
          h_in = eta(i,j+1) + (0.5-cfl)*(eta(i,j+1)-eta(i,j+2))  ! internal

          H_v = BT_OBC%H_v(i,J)
          vel_prev = vbt(i,J)
          vbt(i,J) = 0.5*((v_inlet+BT_OBC%vbt_outer(i,J)) + &
              (BT_OBC%Cg_v(i,J)/H_v) * (BT_OBC%eta_outer_v(i,J)-h_in))

          vel_trans = (1.0-bebt)*vel_prev + bebt*vbt(i,J)
        elseif (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_E) then
          if ((ubt(I-1,j)+ubt(I-1,j+1)) > 0.0) then
            vbt(i,J) = 2.0*vbt(i-1,J)-vbt(i-2,J)
          else
            vbt(i,J) = BT_OBC%vbt_outer(i,J)
          endif
          vel_trans = vbt(i,J)
!!!!!!!!!!!!!!!!!!! CLAMPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       cfl = dtbt * BT_OBC%Cg_v(i,J) * G%IdyCv(i,J)           !
!       vbt(i,J) = (vbt(i-1,J) + CFL*vbt(i,J)) / (1.0 + CFL)  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_W) then
          if ((ubt(I,j)+ubt(I,j+1)) < 0.0) then
            vbt(i,J) = 2.0*vbt(i+1,J)-vbt(i+2,J)
          else
            vbt(i,J) = BT_OBC%vbt_outer(i,J)
          endif
          vel_trans = vbt(i,J)
!!!!!!!!!!!!!!!!!! CLAMPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       cfl = dtbt * BT_OBC%Cg_v(i,J) * G%IdyCv(i,J)           !
!       vbt(i,J) = (vbt(i-1,J) + CFL*vbt(i,J)) / (1.0 + CFL)  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif
      endif

      if (OBC%OBC_segment_v(i,J) /= OBC_SIMPLE) then
        if (use_BT_cont) then
           vhbt(i,J) = find_vhbt(vel_trans,BTCL_v(i,J)) + vhbt0(i,J)
        else
           vhbt(i,J) = vel_trans*Datv(i,J) + vhbt0(i,J)
        endif
      endif

      vbt_trans(i,J) = vel_trans
    endif ; enddo ; enddo
  endif

end subroutine apply_velocity_OBCs

subroutine apply_eta_OBCs(OBC, eta, ubt, vbt, BT_OBC, G, MS, halo, dtbt)
  type(ocean_OBC_type),                  pointer       :: OBC
  type(memory_size_type),                intent(in)    :: MS
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(inout) :: eta
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: ubt
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: vbt
  type(BT_OBC_type),                     intent(in)    :: BT_OBC
  type(ocean_grid_type),                 intent(inout) :: G
  integer,                               intent(in)    :: halo
  real,                                  intent(in)    :: dtbt
!   This subroutine applies the open boundary conditions on the free surface
! height, as coded by Mehmet Ilicak.

! Arguments: OBC - an associated pointer to an OBC type.
!  (inout)   eta - The barotropic free surface height anomaly or
!                  column mass anomaly, in m or kg m-2.
!  (inout)   ubt - the zonal barotropic velocity, in m s-1.
!  (inout)   vbt - the meridional barotropic velocity, in m s-1.
!  (in)      G - The ocean's grid structure.
!  (in)      MS - A type that describes the memory sizes of the argument arrays.
!  (in)      halo - The extra halo size to use here.
!  (in)      dtbt - The time step, in s.

  real :: H_u         ! The total thickness at the u-point, in m or kg m-2.
  real :: H_v         ! The total thickness at the v-point, in m or kg m-2.
  real :: cfl         ! The CFL number at the point in question, ND.
  real :: u_inlet
  real :: v_inlet
  real :: h_in
  integer :: i, j, is, ie, js, je
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo

  if ((OBC%Flather_u_BCs_exist_globally) .and. apply_u_OBCs) then
    do j=js,je ; do I=is-1,ie ; if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
      if (OBC%segment(OBC%OBC_segment_u(I,j))%Flather) then
        if (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_E) then
          cfl = dtbt * BT_OBC%Cg_u(I,j) * G%IdxCu(I,j)            ! CFL
          u_inlet = cfl*ubt(I-1,j) + (1.0-cfl)*ubt(I,j)          ! Valid for cfl <1
!          h_in = 2.0*cfl*eta(i,j) + (1.0-2.0*cfl)*eta(i+1,j)    ! external
          h_in = eta(i,j) + (0.5-cfl)*(eta(i,j)-eta(i-1,j))      ! internal

          H_u = BT_OBC%H_u(I,j)
          eta(i+1,j) = 2.0 * 0.5*((BT_OBC%eta_outer_u(I,j)+h_in) + &
              (H_u/BT_OBC%Cg_u(I,j))*(u_inlet-BT_OBC%ubt_outer(I,j))) - eta(i,j)
        elseif (OBC%segment(OBC%OBC_segment_u(I,j))%direction == OBC_DIRECTION_W) then
          cfl = dtbt*BT_OBC%Cg_u(I,j)*G%IdxCu(I,j)                ! CFL
          u_inlet = cfl*ubt(I+1,j) + (1.0-cfl)*ubt(I,j)          ! Valid for cfl <1
!          h_in = 2.0*cfl*eta(i+1,j) + (1.0-2.0*cfl)*eta(i,j)    ! external
          h_in = eta(i+1,j) + (0.5-cfl)*(eta(i+1,j)-eta(i+2,j))  ! internal

          H_u = BT_OBC%H_u(I,j)
          eta(i,j) = 2.0 * 0.5*((BT_OBC%eta_outer_u(I,j)+h_in) + &
              (H_u/BT_OBC%Cg_u(I,j))*(BT_OBC%ubt_outer(I,j)-u_inlet)) - eta(i+1,j)
        endif
      endif
    endif ; enddo ; enddo
  endif

  if ((OBC%Flather_v_BCs_exist_globally) .and. apply_v_OBCs) then
    do J=js-1,je ; do i=is,ie ; if (OBC%OBC_segment_v(i,J) /= OBC_NONE) then
      if (OBC%segment(OBC%OBC_segment_v(i,J))%Flather) then
        if (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_N) then
          cfl = dtbt*BT_OBC%Cg_v(i,J)*G%IdyCv(i,J)                ! CFL
          v_inlet = cfl*vbt(i,J-1) + (1.0-cfl)*vbt(i,J)          ! Valid for cfl <1
!          h_in = 2.0*cfl*eta(i,j) + (1.0-2.0*cfl)*eta(i,j+1)    ! external
          h_in = eta(i,j) + (0.5-cfl)*(eta(i,j)-eta(i,j-1))      ! internal

          H_v = BT_OBC%H_v(i,J)
          eta(i,j+1) = 2.0 * 0.5*((BT_OBC%eta_outer_v(i,J)+h_in) + &
              (H_v/BT_OBC%Cg_v(i,J))*(v_inlet-BT_OBC%vbt_outer(i,J))) - eta(i,j)
        elseif (OBC%segment(OBC%OBC_segment_v(i,J))%direction == OBC_DIRECTION_S) then
          cfl = dtbt*BT_OBC%Cg_v(i,J)*G%IdyCv(i,J)                ! CFL
          v_inlet = cfl*vbt(i,J+1) + (1.0-cfl)*vbt(i,J)          ! Valid for cfl <1
!          h_in = 2.0*cfl*eta(i,j+1) + (1.0-2.0*cfl)*eta(i,j)    ! external
          h_in = eta(i,j+1) + (0.5-cfl)*(eta(i,j+1)-eta(i,j+2))  ! internal

          H_v = BT_OBC%H_v(i,J)
          eta(i,j) = 2.0 * 0.5*((BT_OBC%eta_outer_v(i,J)+h_in) + &
              (H_v/BT_OBC%Cg_v(i,J))*(BT_OBC%vbt_outer(i,J)-v_inlet)) - eta(i,j+1)
        endif
      endif
    endif ; enddo ; enddo
  endif

end subroutine apply_eta_OBCs

subroutine set_up_BT_OBC(OBC, eta, BT_OBC, G, GV, MS, halo, use_BT_cont, Datu, Datv, BTCL_u, BTCL_v)
  type(ocean_OBC_type),                  pointer       :: OBC
  type(memory_size_type),                intent(in)    :: MS
  real, dimension(SZIW_(MS),SZJW_(MS)),  intent(in)    :: eta
  type(BT_OBC_type),                     intent(inout) :: BT_OBC
  type(ocean_grid_type),                 intent(inout) :: G
  type(verticalGrid_type),               intent(in)    :: GV
  integer,                               intent(in)    :: halo
  logical,                               intent(in)    :: use_BT_cont
  real, dimension(SZIBW_(MS),SZJW_(MS)), intent(in)    :: Datu
  real, dimension(SZIW_(MS),SZJBW_(MS)), intent(in)    :: Datv
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), intent(in) :: BTCL_u
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), intent(in) :: BTCL_v
!   This subroutine sets up the private structure used to apply the open
! boundary conditions, as developed by Mehmet Ilicak.

! Arguments: OBC - an associated pointer to an OBC type.
!  (in)      eta - The barotropic free surface height anomaly or
!                  column mass anomaly, in m or kg m-2.
!  (in)      BT_OBC - A structure with the private barotropic arrays related
!                     to the open boundary conditions, set by set_up_BT_OBC.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      MS - A type that describes the memory sizes of the argument arrays.
!  (in)      halo - The extra halo size to use here.
!  (in)      dtbt - The time step, in s.
!  (in)      use_BT_cont - If true, use the BT_cont_types to calculate transports.
!  (in)      Datu - A fixed estimate of the face areas at u points.
!  (in)      Datv - A fixed estimate of the face areas at u points.
!  (in)      BCTL_u - Structures of information used for a dynamic estimate
!  (in)      BCTL_v - of the face areas at u- and v- points.

  integer :: i, j, k, is, ie, js, je, nz, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: isdw, iedw, jsdw, jedw
  logical :: OBC_used
  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  isdw = MS%isdw ; iedw = MS%iedw ; jsdw = MS%jsdw ; jedw = MS%jedw

  if ((isdw < isd) .or. (jsdw < jsd)) then
    call MOM_error(FATAL, "set_up_BT_OBC: Open boundary conditions are not "//&
                           "yet fully implemented with wide barotropic halos.")
  endif

  allocate(BT_OBC%Cg_u(isdw-1:iedw,jsdw:jedw))        ; BT_OBC%Cg_u(:,:) = 0.0
  allocate(BT_OBC%H_u(isdw-1:iedw,jsdw:jedw))         ; BT_OBC%H_u(:,:) = 0.0
  allocate(BT_OBC%uhbt(isdw-1:iedw,jsdw:jedw))        ; BT_OBC%uhbt(:,:) = 0.0
  allocate(BT_OBC%ubt_outer(isdw-1:iedw,jsdw:jedw))   ; BT_OBC%ubt_outer(:,:) = 0.0
  allocate(BT_OBC%eta_outer_u(isdw-1:iedw,jsdw:jedw)) ; BT_OBC%eta_outer_u(:,:) = 0.0

  allocate(BT_OBC%Cg_v(isdw:iedw,jsdw-1:jedw))        ; BT_OBC%Cg_v(:,:) = 0.0
  allocate(BT_OBC%H_v(isdw:iedw,jsdw-1:jedw))         ; BT_OBC%H_v(:,:) = 0.0
  allocate(BT_OBC%vhbt(isdw:iedw,jsdw-1:jedw))        ; BT_OBC%vhbt(:,:) = 0.0
  allocate(BT_OBC%vbt_outer(isdw:iedw,jsdw-1:jedw))   ; BT_OBC%vbt_outer(:,:) = 0.0
  allocate(BT_OBC%eta_outer_v(isdw:iedw,jsdw-1:jedw)) ; BT_OBC%eta_outer_v(:,:)=0.0

  if (apply_u_OBCs) then
    if (OBC%specified_u_BCs_exist_globally) then
      do k=1,nz ; do j=js,je ; do I=is-1,ie
        BT_OBC%uhbt(I,j) = BT_OBC%uhbt(I,j) + OBC%uh(I,j,k)
      enddo ; enddo ; enddo
    endif
    do j=js,je ; do I=is-1,ie ; if (OBC%OBC_segment_u(I,j) /= OBC_NONE) then
      if (OBC%segment(OBC%OBC_segment_u(I,j))%specified) then
        if (use_BT_cont) then
          BT_OBC%ubt_outer(I,j) = uhbt_to_ubt(BT_OBC%uhbt(I,j),BTCL_u(I,j))
        else
          if (Datu(I,j) > 0.0) BT_OBC%ubt_outer(I,j) = BT_OBC%uhbt(I,j) / Datu(I,j)
        endif
      else
        BT_OBC%Cg_u(I,j) = SQRT(GV%g_prime(1)*(0.5* &
                                (G%bathyT(i,j) + G%bathyT(i+1,j))))
        if (GV%Boussinesq) then
          BT_OBC%H_u(I,j) = 0.5*((G%bathyT(i,j) + eta(i,j)) + &
                                 (G%bathyT(i+1,j) + eta(i+1,j)))
        else
          BT_OBC%H_u(I,j) = 0.5*(eta(i,j) + eta(i+1,j))
        endif
      endif
    endif ; enddo ; enddo
    if (associated(OBC%ubt_outer)) then ; do j=js,je ; do I=is-1,ie
      BT_OBC%ubt_outer(I,j) = OBC%ubt_outer(I,j)
    enddo ; enddo ; endif
    if (associated(OBC%eta_outer_u)) then ; do j=js,je ; do I=is-1,ie
      BT_OBC%eta_outer_u(I,j) = OBC%eta_outer_u(I,j)
    enddo ; enddo ; endif
  endif
  if (apply_v_OBCs) then
    if (OBC%specified_v_BCs_exist_globally) then
      do k=1,nz ; do J=js-1,je ; do i=is,ie
        BT_OBC%vhbt(i,J) = BT_OBC%vhbt(i,J) + OBC%vh(i,J,k)
      enddo ; enddo ; enddo
    endif

    do J=js-1,je ; do i=is,ie ; if (OBC%OBC_segment_v(i,J) /= OBC_NONE) then
      if (OBC%OBC_segment_v(i,J) == OBC_SIMPLE) then
        if (use_BT_cont) then
          BT_OBC%vbt_outer(i,J) = vhbt_to_vbt(BT_OBC%vhbt(i,J),BTCL_v(i,J))
        else
          if (Datv(i,J) > 0.0) BT_OBC%vbt_outer(i,J) = BT_OBC%vhbt(i,J) / Datv(i,J)
        endif
      else
        BT_OBC%Cg_v(i,J) = SQRT(GV%g_prime(1)*(0.5* &
                                (G%bathyT(i,j) + G%bathyT(i,j+1))))
        if (GV%Boussinesq) then
          BT_OBC%H_v(i,J) = 0.5*((G%bathyT(i,j) + eta(i,j)) + &
                                 (G%bathyT(i,j+1) + eta(i,j+1)))
        else
          BT_OBC%H_v(i,J) = 0.5*(eta(i,j) + eta(i,j+1))
        endif
      endif
    endif ; enddo ; enddo
    if (associated(OBC%vbt_outer)) then ; do J=js-1,je ; do i=is,ie
      BT_OBC%vbt_outer(i,J) = OBC%vbt_outer(i,J)
    enddo ; enddo ; endif
    if (associated(OBC%eta_outer_v)) then ; do J=js-1,je ; do i=is,ie
      BT_OBC%eta_outer_v(i,J) = OBC%eta_outer_v(i,J)
    enddo ; enddo ; endif
  endif
end subroutine set_up_BT_OBC

subroutine destroy_BT_OBC(BT_OBC)
  type(BT_OBC_type), intent(inout) :: BT_OBC

  deallocate(BT_OBC%Cg_u)
  deallocate(BT_OBC%H_u)
  deallocate(BT_OBC%uhbt)
  deallocate(BT_OBC%ubt_outer)
  deallocate(BT_OBC%eta_outer_u)

  deallocate(BT_OBC%Cg_v)
  deallocate(BT_OBC%H_v)
  deallocate(BT_OBC%vhbt)
  deallocate(BT_OBC%vbt_outer)
  deallocate(BT_OBC%eta_outer_v)
end subroutine destroy_BT_OBC


subroutine legacy_btcalc(h, G, GV, CS, h_u, h_v, may_use_default)
  type(ocean_grid_type),                  intent(inout) :: G
  type(verticalGrid_type),                intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: h
  type(legacy_barotropic_CS),             pointer       :: CS
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in), optional :: h_u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in), optional :: h_v
  logical,                                   intent(in), optional :: may_use_default
!   btcalc calculates the barotropic velocities from the full velocity and
! thickness fields, determines the fraction of the total water column in each
! layer at velocity points, and determines a corrective fictitious mass source
! that will drive the barotropic estimate of the free surface height toward the
! baroclinic estimate.

! Arguments: h - Layer thickness, in m or kg m-2 (H in later comments).
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 barotropic_init.
!  (in,opt)  h_u, h_v - The specified thicknesses at u- and v- points, in m or kg m-2.

! All of these variables are in the same units as h - usually m or kg m-2.
  real :: hatutot(SZIB_(G))    ! The sum of the layer thicknesses
  real :: hatvtot(SZI_(G))     ! interpolated to the u & v grid points.
  real :: Ihatutot(SZIB_(G))   ! Ihatutot and Ihatvtot are the inverses
  real :: Ihatvtot(SZI_(G))    ! of hatutot and hatvtot, both in H-1.
  real :: h_arith              ! The arithmetic mean thickness, in H.
  real :: h_harm               ! The harmonic mean thicknesses, in H.
  real :: h_neglect            ! A thickness that is so small it is usually lost
                               ! in roundoff and can be neglected, in H.
  real :: wt_arith             ! The nondimensional weight for the arithmetic
                               ! mean thickness.  The harmonic mean uses
                               ! a weight of (1 - wt_arith).
  real :: Rh                   ! A ratio of summed thicknesses, nondim.
  real :: htot(SZI_(G),SZJ_(G)) !   The sum of the layer thicknesses, in H.
  real :: e_u(SZIB_(G),SZK_(G)+1) !   The interface heights at u-velocity and
  real :: e_v(SZI_(G),SZK_(G)+1)  ! v-velocity points in H.
  real :: D_shallow_u(SZI_(G)) ! The shallower of the adjacent depths in H.
  real :: D_shallow_v(SZIB_(G))! The shallower of the adjacent depths in H.

  logical :: use_default, test_dflt
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, i, j, k

!    This section interpolates thicknesses onto u & v grid points with the
! second order accurate estimate h = 2*(h+ * h-)/(h+ + h-).
  if (.not.associated(CS)) call MOM_error(FATAL, &
      "legacy_btcalc: Module MOM_legacy_barotropic must be initialized before it is used.")
  if (.not.CS%split) return

  use_default = .false.
  test_dflt = .false. ; if (present(may_use_default)) test_dflt = may_use_default

  if (test_dflt) then
    if (.not.((present(h_u) .and. present(h_v)) .or. &
              (CS%hvel_scheme == HARMONIC) .or. (CS%hvel_scheme == HYBRID) .or.&
              (CS%hvel_scheme == ARITHMETIC))) use_default = .true.
  else
    if (.not.((present(h_u) .and. present(h_v)) .or. &
              (CS%hvel_scheme == HARMONIC) .or. (CS%hvel_scheme == HYBRID) .or.&
              (CS%hvel_scheme == ARITHMETIC))) call MOM_error(FATAL, &
        "legacy_btcalc: Inconsistent settings of optional arguments and hvel_scheme.")
  endif

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  h_neglect = GV%H_subroundoff

  !   This estimates the fractional thickness of each layer at the velocity
  ! points, using a harmonic mean estimate.
!$OMP parallel do default(none) shared(is,ie,js,je,nz,h_u,CS,h_neglect,h,use_default,G,GV) &
!$OMP                          private(hatutot,Ihatutot,e_u,D_shallow_u,h_arith,h_harm,wt_arith)
  do j=js-1,je+1
    if (present(h_u)) then
      do I=is-2,ie+1 ; hatutot(I) = h_u(I,j,1) ; enddo
      do k=2,nz ; do I=is-2,ie+1
        hatutot(I) = hatutot(I) + h_u(I,j,k)
      enddo ; enddo
      do I=is-2,ie+1 ; Ihatutot(I) = 1.0 / (hatutot(I) + h_neglect) ; enddo
      do k=1,nz ; do I=is-2,ie+1
        CS%frhatu(I,j,k) = h_u(I,j,k) * Ihatutot(I)
      enddo ; enddo
    else
      if (CS%hvel_scheme == ARITHMETIC) then
        do I=is-2,ie+1
          CS%frhatu(I,j,1) = 0.5 * (h(i+1,j,1) + h(i,j,1))
          hatutot(I) = CS%frhatu(I,j,1)
        enddo
        do k=2,nz ; do I=is-2,ie+1
          CS%frhatu(I,j,k) = 0.5 * (h(i+1,j,k) + h(i,j,k))
          hatutot(I) = hatutot(I) + CS%frhatu(I,j,k)
        enddo ; enddo
      elseif (CS%hvel_scheme == HYBRID .or. use_default) then
        do I=is-2,ie+1
          e_u(I,nz+1) = -0.5 * GV%m_to_H * (G%bathyT(i+1,j) + G%bathyT(i,j))
          D_shallow_u(I) = -GV%m_to_H * min(G%bathyT(i+1,j), G%bathyT(i,j))
          hatutot(I) = 0.0
        enddo
        do k=nz,1,-1 ; do I=is-2,ie+1
          e_u(I,K) = e_u(I,K+1) + 0.5 * (h(i+1,j,k) + h(i,j,k))
          h_arith = 0.5 * (h(i+1,j,k) + h(i,j,k))
          if (e_u(I,K+1) >= D_shallow_u(I)) then
            CS%frhatu(I,j,k) = h_arith
          else
            h_harm = (h(i+1,j,k) * h(i,j,k)) / (h_arith + h_neglect)
            if (e_u(I,K) <= D_shallow_u(I)) then
              CS%frhatu(I,j,k) = h_harm
            else
              wt_arith = (e_u(I,K) - D_shallow_u(I)) / (h_arith + h_neglect)
              CS%frhatu(I,j,k) = wt_arith*h_arith + (1.0-wt_arith)*h_harm
            endif
          endif
          hatutot(I) = hatutot(I) + CS%frhatu(I,j,k)
        enddo ; enddo
      elseif (CS%hvel_scheme == HARMONIC) then
        do I=is-2,ie+1
          CS%frhatu(I,j,1) = 2.0*(h(i+1,j,1) * h(i,j,1)) / &
                             ((h(i+1,j,1) + h(i,j,1)) + h_neglect)
          hatutot(I) = CS%frhatu(I,j,1)
        enddo
        do k=2,nz ; do I=is-2,ie+1
          CS%frhatu(I,j,k) = 2.0*(h(i+1,j,k) * h(i,j,k)) / &
                             ((h(i+1,j,k) + h(i,j,k)) + h_neglect)
          hatutot(I) = hatutot(I) + CS%frhatu(I,j,k)
        enddo ; enddo
      endif
      do I=is-2,ie+1 ; Ihatutot(I) = 1.0 / (hatutot(I) + h_neglect) ; enddo
      do k=1,nz ; do I=is-2,ie+1
        CS%frhatu(I,j,k) = CS%frhatu(I,j,k) * Ihatutot(I)
      enddo ; enddo
    endif
  enddo

!$OMP parallel do default(none) shared(is,ie,js,je,nz,CS,G,GV,h_v,h_neglect,h,use_default) &
!$OMP                          private(hatvtot,Ihatvtot,e_v,D_shallow_v,h_arith,h_harm,wt_arith)
  do J=js-2,je+1
    if (present(h_v)) then
      do i=is-1,ie+1 ; hatvtot(i) = h_v(i,J,1) ; enddo
      do k=2,nz ; do i=is-1,ie+1
        hatvtot(i) = hatvtot(i) + h_v(i,J,k)
      enddo ; enddo
      do i=is-1,ie+1 ; Ihatvtot(i) = 1.0 / (hatvtot(i) + h_neglect) ; enddo
      do k=1,nz ; do i=is-1,ie+1
        CS%frhatv(i,J,k) = h_v(i,J,k) * Ihatvtot(i)
      enddo ; enddo
    else
      if (CS%hvel_scheme == ARITHMETIC) then
        do i=is-1,ie+1
          CS%frhatv(i,J,1) = 0.5 * (h(i,j+1,1) + h(i,j,1))
          hatvtot(i) = CS%frhatv(i,J,1)
        enddo
        do k=2,nz ; do i=is-1,ie+1
          CS%frhatv(i,J,k) = 0.5 * (h(i,j+1,k) + h(i,j,k))
          hatvtot(i) = hatvtot(i) + CS%frhatv(i,J,k)
        enddo ; enddo
      elseif (CS%hvel_scheme == HYBRID .or. use_default) then
        do i=is-1,ie+1
          e_v(i,nz+1) = -0.5 * GV%m_to_H * (G%bathyT(i,j+1) + G%bathyT(i,j))
          D_shallow_v(I) = -GV%m_to_H * min(G%bathyT(i,j+1), G%bathyT(i,j))
          hatvtot(I) = 0.0
        enddo
        do k=nz,1,-1 ; do i=is-1,ie+1
          e_v(i,K) = e_v(i,K+1) + 0.5 * (h(i,j+1,k) + h(i,j,k))
          h_arith = 0.5 * (h(i,j+1,k) + h(i,j,k))
          if (e_v(i,K+1) >= D_shallow_v(i)) then
            CS%frhatv(i,J,k) = h_arith
          else
            h_harm = (h(i,j+1,k) * h(i,j,k)) / (h_arith + h_neglect)
            if (e_v(i,K) <= D_shallow_v(i)) then
              CS%frhatv(i,J,k) = h_harm
            else
              wt_arith = (e_v(i,K) - D_shallow_v(i)) / (h_arith + h_neglect)
              CS%frhatv(i,J,k) = wt_arith*h_arith + (1.0-wt_arith)*h_harm
            endif
          endif
          hatvtot(i) = hatvtot(i) + CS%frhatv(i,J,k)
        enddo ; enddo
      elseif (CS%hvel_scheme == HARMONIC) then
        do i=is-1,ie+1
          CS%frhatv(i,J,1) = 2.0*(h(i,j+1,1) * h(i,j,1)) / &
                             ((h(i,j+1,1) + h(i,j,1)) + h_neglect)
          hatvtot(i) = CS%frhatv(i,J,1)
        enddo
        do k=2,nz ; do i=is-1,ie+1
          CS%frhatv(i,J,k) = 2.0*(h(i,j+1,k) * h(i,j,k)) / &
                             ((h(i,j+1,k) + h(i,j,k)) + h_neglect)
          hatvtot(i) = hatvtot(i) + CS%frhatv(i,J,k)
        enddo ; enddo
      endif
      do i=is-1,ie+1 ; Ihatvtot(i) = 1.0 / (hatvtot(i) + h_neglect) ; enddo
      do k=1,nz ; do i=is-1,ie+1
        CS%frhatv(i,J,k) = CS%frhatv(i,J,k) * Ihatvtot(i)
      enddo ; enddo
    endif
  enddo

  if (CS%rescale_D_bt) then
    do j=js-2,je+2 ; do i=is-2,ie+2 ; htot(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js-2,je+2 ; do i=is-2,ie+2
      htot(i,j) = htot(i,j) + h(i,j,k)
    enddo ; enddo ; enddo

!$OMP parallel do default(none) shared(is,ie,js,je,nz,h,htot,h_neglect,CS) &
!$OMP                          private(hatutot, h_harm, Rh)
    do j=js-1,je+1
      do I=is-2,ie+1 ; hatutot(I) = 0.0 ; enddo
      do k=1,nz ; do I=is-2,ie+1
        h_harm = 2.0*(h(i+1,j,k) * h(i,j,k)) / &
                 ((h(i+1,j,k) + h(i,j,k)) + h_neglect)
        hatutot(I) = hatutot(I) + h_harm
      enddo ; enddo
      do I=is-2,ie+1
        Rh = hatutot(I) * (htot(i+1,j) + htot(i,j)) / &
             (2.0*(htot(i+1,j) * htot(i,j) + h_neglect**2))
        CS%Datu_res(I,j) = 1.0
        ! This curve satisfies F(1/2) = 1; F'(1/2) = 0; F(x) ~ x for x<<1.
        if (Rh < 0.5) CS%Datu_res(I,j) = Rh * ((4.0-7.0*Rh) / (2.0-3.0*Rh)**2)
      enddo
    enddo

!$OMP parallel do default(none) shared(is,ie,js,je,nz,h,htot,h_neglect,CS) &
!$OMP                          private(hatvtot, h_harm, Rh)
    do J=js-2,je+1
      do i=is-1,ie+1 ; hatvtot(i) = 0.0 ; enddo
      do k=1,nz ; do i=is-1,ie+1
        h_harm = 2.0*(h(i,j+1,k) * h(i,j,k)) / &
                 ((h(i,j+1,k) + h(i,j,k)) + h_neglect)
        hatvtot(i) = hatvtot(i) + h_harm
      enddo ; enddo
      do i=is-1,ie+1
        Rh = hatvtot(i) * (htot(i,j+1) + htot(i,j)) / &
             (2.0*(htot(i,j+1) * htot(i,j) + h_neglect**2))
        CS%Datv_res(i,J) = 1.0
        ! This curve satisfies F(1/2) = 1; F'(1/2) = 0; F(x) ~ x for x<<1.
        if (Rh < 0.5) CS%Datv_res(i,J) = Rh * ((4.0-7.0*Rh) / (2.0-3.0*Rh)**2)
      enddo
    enddo
  endif

  if (CS%debug) then
    call uchksum(CS%frhatu, "btcalc frhatu",G%HI,haloshift=1)
    call vchksum(CS%frhatv, "btcalc frhatv",G%HI,haloshift=1)
    call hchksum(GV%H_to_m*h, "btcalc h",G%HI,haloshift=1)
  endif

end subroutine legacy_btcalc

function find_uhbt(u, BTC) result(uhbt)
  real, intent(in) :: u
  type(local_BT_cont_u_type), intent(in) :: BTC
  real :: uhbt ! The result
  ! This function evaluates the zonal transport function.

  if (u == 0.0) then
    uhbt = 0.0
  elseif (u < BTC%uBT_EE) then
    uhbt = (u - BTC%uBT_EE) * BTC%FA_u_EE + BTC%uh_EE
  elseif (u < 0.0) then
    uhbt = u * (BTC%FA_u_E0 + BTC%uh_crvE * u**2)
  elseif (u <= BTC%uBT_WW) then
    uhbt = u * (BTC%FA_u_W0 + BTC%uh_crvW * u**2)
  else ! (u > BTC%uBT_WW)
    uhbt = (u - BTC%uBT_WW) * BTC%FA_u_WW + BTC%uh_WW
  endif
end function find_uhbt

function uhbt_to_ubt(uhbt, BTC, guess) result(ubt)
  real, intent(in) :: uhbt
  type(local_BT_cont_u_type), intent(in) :: BTC
  real, optional, intent(in) :: guess
  real :: ubt ! The result
  !   This function inverts the transport function to determine the barotopic
  ! velocity that is consistent with a given transport.
  ! Arguments: uhbt - The barotropic zonal transport that should be inverted
  !                   for, in units of H m2 s-1.
  !  (in)      BTC - A structure containing various fields that allow the
  !                  barotropic transports to be calculated consistently with
  !                  the layers' continuity equations.
  !  (in,opt)  FA_rat_EE - The current value of the far-east face area divided
  !                        by its value when uhbt was originally calculated, ND.
  !  (in,opt)  FA_rat_WW - The current value of the far-west face area divided
  !                        by its value when uhbt was originally calculated, ND.
  !  (in, opt) guess - A guess at what ubt will be.  The result is not allowed
  !                    to be dramatically larger than guess.
  ! result: ubt - The velocity that gives uhbt transport, in m s-1.
  real :: ubt_min, ubt_max, uhbt_err, derr_du
  real :: uherr_min, uherr_max
  real, parameter :: tol = 1.0e-10
  real :: dvel, vsr  ! Temporary variables used in the limiting the velocity.
  real, parameter :: vs1 = 1.25  ! Nondimensional parameters used in limiting
  real, parameter :: vs2 = 2.0   ! the velocity, starting at vs1, with the
                                 ! maximum increase of vs2, both nondim.
  integer :: itt, max_itt = 20

  ! Find the value of ubt that gives uhbt.
  if (uhbt == 0.0) then
    ubt = 0.0
  elseif (uhbt < BTC%uh_EE) then
    ubt = BTC%uBT_EE + (uhbt - BTC%uh_EE) / BTC%FA_u_EE
  elseif (uhbt < 0.0) then
    ! Iterate to convergence with Newton's method (when bounded) and the
    ! false position method otherwise.  ubt will be negative.
    ubt_min = BTC%uBT_EE ; uherr_min = BTC%uh_EE - uhbt
    ubt_max = 0.0 ; uherr_max = -uhbt
    ! Use a false-position method first guess.
    ubt = BTC%uBT_EE * (uhbt / BTC%uh_EE)
    do itt = 1, max_itt
      uhbt_err = ubt * (BTC%FA_u_E0 + BTC%uh_crvE * ubt**2) - uhbt

      if (abs(uhbt_err) < tol*abs(uhbt)) exit
      if (uhbt_err > 0.0) then ; ubt_max = ubt ; uherr_max = uhbt_err ; endif
      if (uhbt_err < 0.0) then ; ubt_min = ubt ; uherr_min = uhbt_err ; endif

      derr_du = BTC%FA_u_E0 + 3.0 * BTC%uh_crvE * ubt**2
      if ((uhbt_err >= derr_du*(ubt - ubt_min)) .or. &
          (-uhbt_err >= derr_du*(ubt_max - ubt)) .or. (derr_du <= 0.0)) then
        ! Use a false-position method guess.
        ubt = ubt_max + (ubt_min-ubt_max) * (uherr_max / (uherr_max-uherr_min))
      else ! Use Newton's method.
        ubt = ubt - uhbt_err / derr_du
        if (abs(uhbt_err) < (0.01*tol)*abs(ubt_min*derr_du)) exit
      endif
    enddo
  elseif (uhbt <= BTC%uh_WW) then
    ! Iterate to convergence with Newton's method.  ubt will be positive.
    ubt_min = 0.0 ; uherr_min = -uhbt
    ubt_max = BTC%uBT_WW ; uherr_max = BTC%uh_WW - uhbt
    ! Use a false-position method first guess.
    ubt = BTC%uBT_WW * (uhbt / BTC%uh_WW)
    do itt = 1, max_itt
      uhbt_err = ubt * (BTC%FA_u_W0 + BTC%uh_crvW * ubt**2) - uhbt

      if (abs(uhbt_err) < tol*abs(uhbt)) exit
      if (uhbt_err > 0.0) then ; ubt_max = ubt ; uherr_max = uhbt_err ; endif
      if (uhbt_err < 0.0) then ; ubt_min = ubt ; uherr_min = uhbt_err ; endif

      derr_du = BTC%FA_u_W0 + 3.0 * BTC%uh_crvW * ubt**2
      if ((uhbt_err >= derr_du*(ubt - ubt_min)) .or. &
          (-uhbt_err >= derr_du*(ubt_max - ubt)) .or. (derr_du <= 0.0)) then
        ! Use a false-position method guess.
        ubt = ubt_min + (ubt_max-ubt_min) * (-uherr_min / (uherr_max-uherr_min))
      else ! Use Newton's method.
        ubt = ubt - uhbt_err / derr_du
        if (abs(uhbt_err) < (0.01*tol)*(ubt_max*derr_du)) exit
      endif
    enddo
  else ! (uhbt > BTC%uh_WW)
    ubt = BTC%uBT_WW + (uhbt - BTC%uh_WW) / BTC%FA_u_WW
  endif

  if (present(guess)) then
    dvel = abs(ubt) - vs1*abs(guess)
    if (dvel > 0.0) then ! Limit the velocity
      if (dvel < 40.0 * (abs(guess)*(vs2-vs1)) ) then
        vsr = vs2 - (vs2-vs1) * exp(-dvel / (abs(guess)*(vs2-vs1)))
      else  ! The exp be less than 4e-18 anyway in this case, so neglect it.
        vsr = vs2
      endif
      ubt = SIGN(vsr * guess, ubt)
    endif
  endif

end function uhbt_to_ubt

function find_vhbt(v, BTC) result(vhbt)
  real, intent(in) :: v
  type(local_BT_cont_v_type), intent(in) :: BTC
  real :: vhbt ! The result
  ! This function evaluates the meridional transport function.

  if (v == 0.0) then
    vhbt = 0.0
  elseif (v < BTC%vBT_NN) then
    vhbt = (v - BTC%vBT_NN) * BTC%FA_v_NN + BTC%vh_NN
  elseif (v < 0.0) then
    vhbt = v * (BTC%FA_v_N0 + BTC%vh_crvN * v**2)
  elseif (v <= BTC%vBT_SS) then
    vhbt = v * (BTC%FA_v_S0 + BTC%vh_crvS * v**2)
  else ! (v > BTC%vBT_SS)
    vhbt = (v - BTC%vBT_SS) * BTC%FA_v_SS + BTC%vh_SS
  endif
end function find_vhbt

function vhbt_to_vbt(vhbt, BTC, guess) result(vbt)
  real, intent(in) :: vhbt
  type(local_BT_cont_v_type), intent(in) :: BTC
  real, optional, intent(in) :: guess
  real :: vbt ! The result
  !   This function inverts the transport function to determine the barotopic
  ! velocity that is consistent with a given transport.
  ! Arguments: vhbt_in - The barotropic meridional transport that should be
  !                      inverted for, in units of H m2 s-1.
  !  (in)      BTC - A structure containing various fields that allow the
  !                  barotropic transports to be calculated consistently with
  !                  the layers' continuity equations.
  !  (in,opt)  FA_rat_NN - The current value of the far-north face area divided
  !                        by its value when vhbt was originally calculated, ND.
  !  (in,opt)  FA_rat_SS - The current value of the far-south face area divided
  !                        by its value when vhbt was originally calculated, ND.
  !  (in, opt) guess - A guess at what vbt will be.  The result is not allowed
  !                    to be dramatically larger than guess.
  ! result: vbt - The velocity that gives vhbt transport, in m s-1.
  real :: vbt_min, vbt_max, vhbt_err, derr_dv
  real :: vherr_min, vherr_max
  real, parameter :: tol = 1.0e-10
  real :: dvel, vsr  ! Temporary variables used in the limiting the velocity.
  real, parameter :: vs1 = 1.25  ! Nondimensional parameters used in limiting
  real, parameter :: vs2 = 2.0   ! the velocity, starting at vs1, with the
                                 ! maximum increase of vs2, both nondim.
  integer :: itt, max_itt = 20

  ! Find the value of vbt that gives vhbt.
  if (vhbt == 0.0) then
    vbt = 0.0
  elseif (vhbt < BTC%vh_NN) then
    vbt = BTC%vBT_NN + (vhbt - BTC%vh_NN) / BTC%FA_v_NN
  elseif (vhbt < 0.0) then
    ! Iterate to convergence with Newton's method (when bounded) and the
    ! false position method otherwise.  vbt will be negative.
    vbt_min = BTC%vBT_NN ; vherr_min = BTC%vh_NN - vhbt
    vbt_max = 0.0 ; vherr_max = -vhbt
    ! Use a false-position method first guess.
    vbt = BTC%vBT_NN * (vhbt / BTC%vh_NN)
    do itt = 1, max_itt
      vhbt_err = vbt * (BTC%FA_v_N0 + BTC%vh_crvN * vbt**2) - vhbt

      if (abs(vhbt_err) < tol*abs(vhbt)) exit
      if (vhbt_err > 0.0) then ; vbt_max = vbt ; vherr_max = vhbt_err ; endif
      if (vhbt_err < 0.0) then ; vbt_min = vbt ; vherr_min = vhbt_err ; endif

      derr_dv = BTC%FA_v_N0 + 3.0 * BTC%vh_crvN * vbt**2
      if ((vhbt_err >= derr_dv*(vbt - vbt_min)) .or. &
          (-vhbt_err >= derr_dv*(vbt_max - vbt)) .or. (derr_dv <= 0.0)) then
        ! Use a false-position method guess.
        vbt = vbt_max + (vbt_min-vbt_max) * (vherr_max / (vherr_max-vherr_min))
      else ! Use Newton's method.
        vbt = vbt - vhbt_err / derr_dv
        if (abs(vhbt_err) < (0.01*tol)*abs(derr_dv*vbt_min)) exit
      endif
    enddo
  elseif (vhbt <= BTC%vh_SS) then
    ! Iterate to convergence with Newton's method.  vbt will be positive.
    vbt_min = 0.0 ; vherr_min = -vhbt
    vbt_max = BTC%vBT_SS ; vherr_max = BTC%vh_SS - vhbt
    ! Use a false-position method first guess.
    vbt = BTC%vBT_SS * (vhbt / BTC%vh_SS)
    do itt = 1, max_itt
      vhbt_err = vbt * (BTC%FA_v_S0 + BTC%vh_crvS * vbt**2) - vhbt

      if (abs(vhbt_err) < tol*abs(vhbt)) exit
      if (vhbt_err > 0.0) then ; vbt_max = vbt ; vherr_max = vhbt_err ; endif
      if (vhbt_err < 0.0) then ; vbt_min = vbt ; vherr_min = vhbt_err ; endif

      derr_dv = BTC%FA_v_S0 + 3.0 * BTC%vh_crvS * vbt**2
      if ((vhbt_err >= derr_dv*(vbt - vbt_min)) .or. &
          (-vhbt_err >= derr_dv*(vbt_max - vbt)) .or. (derr_dv <= 0.0)) then
        ! Use a false-position method guess.
        vbt = vbt_min + (vbt_max-vbt_min) * (-vherr_min / (vherr_max-vherr_min))
      else ! Use Newton's method.
        vbt = vbt - vhbt_err / derr_dv
        if (abs(vhbt_err) < (0.01*tol)*(vbt_max*derr_dv)) exit
      endif
    enddo
  else ! (vhbt > BTC%vh_SS)
    vbt = BTC%vBT_SS + (vhbt - BTC%vh_SS) / BTC%FA_v_SS
  endif

  if (present(guess)) then
    dvel = abs(vbt) - vs1*abs(guess)
    if (dvel > 0.0) then ! Limit the velocity
      if (dvel < 40.0 * (abs(guess)*(vs2-vs1)) ) then
        vsr = vs2 - (vs2-vs1) * exp(-dvel / (abs(guess)*(vs2-vs1)))
      else  ! The exp be less than 4e-18 anyway in this case, so neglect it.
        vsr = vs2
      endif
      vbt = SIGN(guess * vsr, vbt)
    endif
  endif

end function vhbt_to_vbt


subroutine set_local_BT_cont_types(BT_cont, BTCL_u, BTCL_v, G, MS, BT_Domain, halo)
  type(BT_cont_type),                                    intent(inout) :: BT_cont
  type(memory_size_type),                                intent(in)    :: MS
  type(local_BT_cont_u_type), dimension(SZIBW_(MS),SZJW_(MS)), intent(out) :: BTCL_u
  type(local_BT_cont_v_type), dimension(SZIW_(MS),SZJBW_(MS)), intent(out) :: BTCL_v
  type(ocean_grid_type),                                 intent(inout) :: G
  type(MOM_domain_type),                                intent(inout) :: BT_Domain
  integer,                                     optional, intent(in)    :: halo
!   This subroutine sets up reordered versions of the BT_cont type in the
! local_BT_cont types, which have wide halos properly filled in.
! Arguments: BT_cont - The BT_cont_type input to the barotropic solver.
!  (out)     BTCL_u - A structure with the u information from BT_cont.
!  (out)     BTCL_v - A structure with the v information from BT_cont.
!  (in)      G - The ocean's grid structure.
!  (in)      MS - A type that describes the memory sizes of the argument arrays.
!  (in)      BT_Domain - The domain to use for updating the halos of wide arrays.
!  (in)      halo - The extra halo size to use here.

  real, dimension(SZIBW_(MS),SZJW_(MS)) :: &
    u_polarity, uBT_EE, uBT_WW, FA_u_EE, FA_u_E0, FA_u_W0, FA_u_WW
  real, dimension(SZIW_(MS),SZJBW_(MS)) :: &
    v_polarity, vBT_NN, vBT_SS, FA_v_NN, FA_v_N0, FA_v_S0, FA_v_SS
  real, parameter :: C1_3 = 1.0/3.0
  integer :: i, j, is, ie, js, je, hs
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  hs = 1 ; if (present(halo)) hs = max(halo,0)

  ! Copy the BT_cont arrays into symmetric, potentially wide haloed arrays.
  u_polarity(:,:) = 1.0
  uBT_EE(:,:) = 0.0 ; uBT_WW(:,:) = 0.0
  FA_u_EE(:,:) = 0.0 ; FA_u_E0(:,:) = 0.0 ; FA_u_W0(:,:) = 0.0 ; FA_u_WW(:,:) = 0.0
  do I=is-1,ie ; do j=js,je
    uBT_EE(I,j) = BT_cont%uBT_EE(I,j) ; uBT_WW(I,j) = BT_cont%uBT_WW(I,j)
    FA_u_EE(I,j) = BT_cont%FA_u_EE(I,j) ; FA_u_E0(I,j) = BT_cont%FA_u_E0(I,j)
    FA_u_W0(I,j) = BT_cont%FA_u_W0(I,j) ; FA_u_WW(I,j) = BT_cont%FA_u_WW(I,j)
  enddo ; enddo

  v_polarity(:,:) = 1.0
  vBT_NN(:,:) = 0.0 ; vBT_SS(:,:) = 0.0
  FA_v_NN(:,:) = 0.0 ; FA_v_N0(:,:) = 0.0 ; FA_v_S0(:,:) = 0.0 ; FA_v_SS(:,:) = 0.0
  do i=is,ie ; do J=js-1,je
    vBT_NN(i,J) = BT_cont%vBT_NN(i,J) ; vBT_SS(i,J) = BT_cont%vBT_SS(i,J)
    FA_v_NN(i,J) = BT_cont%FA_v_NN(i,J) ; FA_v_N0(i,J) = BT_cont%FA_v_N0(i,J)
    FA_v_S0(i,J) = BT_cont%FA_v_S0(i,J) ; FA_v_SS(i,J) = BT_cont%FA_v_SS(i,J)
  enddo ; enddo

  if (id_clock_calc_pre > 0) call cpu_clock_end(id_clock_calc_pre)
  if (id_clock_pass_pre > 0) call cpu_clock_begin(id_clock_pass_pre)
  ! Do halo updates on BT_cont.
  call pass_vector(u_polarity, v_polarity, BT_Domain, complete=.false.)
  call pass_vector(uBT_EE, vBT_NN, BT_Domain, complete=.false.)
  call pass_vector(uBT_WW, vBT_SS, BT_Domain, complete=.true.)

  call pass_vector(FA_u_EE, FA_v_NN, BT_Domain, To_All+Scalar_Pair, complete=.false.)
  call pass_vector(FA_u_E0, FA_v_N0, BT_Domain, To_All+Scalar_Pair, complete=.false.)
  call pass_vector(FA_u_W0, FA_v_S0, BT_Domain, To_All+Scalar_Pair, complete=.false.)
  call pass_vector(FA_u_WW, FA_v_SS, BT_Domain, To_All+Scalar_Pair, complete=.true.)
  if (id_clock_pass_pre > 0) call cpu_clock_end(id_clock_pass_pre)
  if (id_clock_calc_pre > 0) call cpu_clock_begin(id_clock_calc_pre)

  do j=js-hs,je+hs ; do I=is-hs-1,ie+hs
    BTCL_u(I,j)%FA_u_EE = FA_u_EE(I,j) ; BTCL_u(I,j)%FA_u_E0 = FA_u_E0(I,j)
    BTCL_u(I,j)%FA_u_W0 = FA_u_W0(I,j) ; BTCL_u(I,j)%FA_u_WW = FA_u_WW(I,j)
    BTCL_u(I,j)%uBT_EE = uBT_EE(I,j)   ; BTCL_u(I,j)%uBT_WW = uBT_WW(I,j)
    ! Check for reversed polarity in the tripolar halo regions.
    if (u_polarity(I,j) < 0.0) then
      call swap(BTCL_u(I,j)%FA_u_EE, BTCL_u(I,j)%FA_u_WW)
      call swap(BTCL_u(I,j)%FA_u_E0, BTCL_u(I,j)%FA_u_W0)
      call swap(BTCL_u(I,j)%uBT_EE,  BTCL_u(I,j)%uBT_WW)
    endif

    BTCL_u(I,j)%uh_EE = BTCL_u(I,j)%uBT_EE * &
        (C1_3 * (2.0*BTCL_u(I,j)%FA_u_E0 + BTCL_u(I,j)%FA_u_EE))
    BTCL_u(I,j)%uh_WW = BTCL_u(I,j)%uBT_WW * &
        (C1_3 * (2.0*BTCL_u(I,j)%FA_u_W0 + BTCL_u(I,j)%FA_u_WW))

    BTCL_u(I,j)%uh_crvE = 0.0 ; BTCL_u(I,j)%uh_crvW = 0.0
    if (abs(BTCL_u(I,j)%uBT_WW) > 0.0) BTCL_u(I,j)%uh_crvW = &
      (C1_3 * (BTCL_u(I,j)%FA_u_WW - BTCL_u(I,j)%FA_u_W0)) / BTCL_u(I,j)%uBT_WW**2
    if (abs(BTCL_u(I,j)%uBT_EE) > 0.0) BTCL_u(I,j)%uh_crvE = &
      (C1_3 * (BTCL_u(I,j)%FA_u_EE - BTCL_u(I,j)%FA_u_E0)) / BTCL_u(I,j)%uBT_EE**2
  enddo ; enddo

  do J=js-hs-1,je+hs ; do i=is-hs,ie+hs
    BTCL_v(i,J)%FA_v_NN = FA_v_NN(i,J) ; BTCL_v(i,J)%FA_v_N0 = FA_v_N0(i,J)
    BTCL_v(i,J)%FA_v_S0 = FA_v_S0(i,J) ; BTCL_v(i,J)%FA_v_SS = FA_v_SS(i,J)
    BTCL_v(i,J)%vBT_NN = vBT_NN(i,J)   ; BTCL_v(i,J)%vBT_SS = vBT_SS(i,J)
    ! Check for reversed polarity in the tripolar halo regions.
    if (v_polarity(i,J) < 0.0) then
      call swap(BTCL_v(i,J)%FA_v_NN, BTCL_v(i,J)%FA_v_SS)
      call swap(BTCL_v(i,J)%FA_v_N0, BTCL_v(i,J)%FA_v_S0)
      call swap(BTCL_v(i,J)%vBT_NN,  BTCL_v(i,J)%vBT_SS)
    endif

    BTCL_v(i,J)%vh_NN = BTCL_v(i,J)%vBT_NN * &
        (C1_3 * (2.0*BTCL_v(i,J)%FA_v_N0 + BTCL_v(i,J)%FA_v_NN))
    BTCL_v(i,J)%vh_SS = BTCL_v(i,J)%vBT_SS * &
        (C1_3 * (2.0*BTCL_v(i,J)%FA_v_S0 + BTCL_v(i,J)%FA_v_SS))

    BTCL_v(i,J)%vh_crvN = 0.0 ; BTCL_v(i,J)%vh_crvS = 0.0
    if (abs(BTCL_v(i,J)%vBT_SS) > 0.0) BTCL_v(i,J)%vh_crvS = &
      (C1_3 * (BTCL_v(i,J)%FA_v_SS - BTCL_v(i,J)%FA_v_S0)) / BTCL_v(i,J)%vBT_SS**2
    if (abs(BTCL_v(i,J)%vBT_NN) > 0.0) BTCL_v(i,J)%vh_crvN = &
      (C1_3 * (BTCL_v(i,J)%FA_v_NN - BTCL_v(i,J)%FA_v_N0)) / BTCL_v(i,J)%vBT_NN**2
  enddo ; enddo

end subroutine set_local_BT_cont_types

subroutine BT_cont_to_face_areas(BT_cont, Datu, Datv, G, MS, halo, maximize)
  type(BT_cont_type),                         intent(inout) :: BT_cont
  type(memory_size_type),                     intent(in)    :: MS
  real, dimension(MS%isdw-1:MS%iedw,MS%jsdw:MS%jedw), intent(out)   :: Datu
  real, dimension(MS%isdw:MS%iedw,MS%jsdw-1:MS%jedw), intent(out)   :: Datv
  type(ocean_grid_type),                      intent(in)  :: G
  integer,                          optional, intent(in)  :: halo
  logical,                          optional, intent(in)  :: maximize
  !   This subroutine uses the BTCL types to find typical or maximum face
  ! areas, which can then be used for finding wave speeds, etc.
  logical :: find_max
  integer :: i, j, is, ie, js, je, hs
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  hs = 1 ; if (present(halo)) hs = max(halo,0)
  find_max = .false. ; if (present(maximize)) find_max = maximize

  if (find_max) then
    do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
      Datu(I,j) = max(BT_cont%FA_u_EE(I,j), BT_cont%FA_u_E0(I,j), &
                      BT_cont%FA_u_W0(I,j), BT_cont%FA_u_WW(I,j))
    enddo ; enddo
    do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
      Datv(i,J) = max(BT_cont%FA_v_NN(i,J), BT_cont%FA_v_N0(i,J), &
                      BT_cont%FA_v_S0(i,J), BT_cont%FA_v_SS(i,J))
    enddo ; enddo
  else
    do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
      Datu(I,j) = 0.5 * (BT_cont%FA_u_E0(I,j) + BT_cont%FA_u_W0(I,j))
    enddo ; enddo
    do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
      Datv(i,J) = 0.5 * (BT_cont%FA_v_N0(i,J) + BT_cont%FA_v_S0(i,J))
    enddo ; enddo
  endif
end subroutine BT_cont_to_face_areas

subroutine swap(a,b)
  real, intent(inout) :: a, b
  real :: tmp
  tmp = a ; a = b ; b = tmp
end subroutine swap

subroutine find_face_areas(Datu, Datv, G, GV, CS, MS, rescale_faces, eta, halo, add_max)
  type(memory_size_type),                   intent(in) :: MS
  real, dimension(MS%isdw-1:MS%iedw,MS%jsdw:MS%jedw), intent(out)   :: Datu
  real, dimension(MS%isdw:MS%iedw,MS%jsdw-1:MS%jedw), intent(out)   :: Datv
  type(ocean_grid_type),                    intent(in) :: G
  type(verticalGrid_type),                  intent(in) :: GV
  type(legacy_barotropic_CS),               pointer    :: CS
  logical,                        optional, intent(in) :: rescale_faces
  real, dimension(MS%isdw:MS%iedw,MS%jsdw:MS%jedw), optional, intent(in) :: eta
  integer,                        optional, intent(in) :: halo
  real,                           optional, intent(in) :: add_max
! Arguments: Datu - The open zonal face area, in H m (m2 or kg m-1).
!  (out)     Datv - The open meridional face area, in H m (m2 or kg m-1).
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 barotropic_init.
!  (in)      MS - A type that describes the memory sizes of the argument arrays.
!  (in, opt) rescale_faces - If true, rescale the face areas by Datu_res, etc.
!  (in, opt) eta - The barotropic free surface height anomaly or
!                  column mass anomaly, in m or kg m-2.
!  (in, opt) halo - The halo size to use, default = 1.
!  (in, opt) add_max - A value to add to the maximum depth (used to overestimate
!                      the external wave speed) in m.


!   This subroutine determines the open face areas of cells for calculating
! the barotropic transport.
  real :: H1, H2      ! Temporary total thicknesses, in m or kg m-2.
  logical :: rescale
  integer :: i, j, is, ie, js, je, hs
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  hs = 1 ; if (present(halo)) hs = max(halo,0)
  rescale = .false. ; if (present(rescale_faces)) rescale = rescale_faces

!$OMP parallel default(none) shared(is,ie,js,je,hs,eta,GV,CS,Datu,Datv,add_max,rescale) &
!$OMP                       private(H1,H2)
  if (present(eta)) then
    ! The use of harmonic mean thicknesses ensure positive definiteness.
    if (GV%Boussinesq) then
!$OMP do
      do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
        H1 = CS%bathyT(i,j) + eta(i,j) ; H2 = CS%bathyT(i+1,j) + eta(i+1,j)
        Datu(I,j) = 0.0 ; if ((H1 > 0.0) .and. (H2 > 0.0)) &
        Datu(I,j) = CS%dy_Cu(I,j) * (2.0 * H1 * H2) / (H1 + H2)
!       Datu(I,j) = CS%dy_Cu(I,j) * 0.5 * (H1 + H2)
      enddo; enddo
!$OMP do
      do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
        H1 = CS%bathyT(i,j) + eta(i,j) ; H2 = CS%bathyT(i,j+1) + eta(i,j+1)
        Datv(i,J) = 0.0 ; if ((H1 > 0.0) .and. (H2 > 0.0)) &
        Datv(i,J) = CS%dx_Cv(i,J) * (2.0 * H1 * H2) / (H1 + H2)
!       Datv(i,J) = CS%dy_v(i,J) * 0.5 * (H1 + H2)
      enddo; enddo
    else
!$OMP do
      do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
        Datu(I,j) = 0.0 ; if ((eta(i,j) > 0.0) .and. (eta(i+1,j) > 0.0)) &
        Datu(I,j) = CS%dy_Cu(I,j) * (2.0 * eta(i,j) * eta(i+1,j)) / &
                                  (eta(i,j) + eta(i+1,j))
        ! Datu(I,j) = CS%dy_Cu(I,j) * 0.5 * (eta(i,j) + eta(i+1,j))
      enddo; enddo
!$OMP do
      do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
        Datv(i,J) = 0.0 ; if ((eta(i,j) > 0.0) .and. (eta(i,j+1) > 0.0)) &
        Datv(i,J) = CS%dx_Cv(i,J) * (2.0 * eta(i,j) * eta(i,j+1)) / &
                                  (eta(i,j) + eta(i,j+1))
        ! Datv(i,J) = CS%dy_v(i,J) * 0.5 * (eta(i,j) + eta(i,j+1))
      enddo; enddo
    endif
  elseif (present(add_max)) then
!$OMP do
    do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
      Datu(I,j) = CS%dy_Cu(I,j) * GV%m_to_H * &
                  (max(CS%bathyT(i+1,j), CS%bathyT(i,j)) + add_max)
    enddo ; enddo
!$OMP do
    do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
      Datv(i,J) = CS%dx_Cv(i,J) * GV%m_to_H * &
                  (max(CS%bathyT(i,j+1), CS%bathyT(i,j)) + add_max)
    enddo ; enddo
  else
!$OMP do
    do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
      Datu(I,j) = 2.0*CS%dy_Cu(I,j) * GV%m_to_H * &
                  (CS%bathyT(i+1,j) * CS%bathyT(i,j)) / &
                  (CS%bathyT(i+1,j) + CS%bathyT(i,j))
    enddo ; enddo
!$OMP do
    do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
      Datv(i,J) = 2.0*CS%dx_Cv(i,J) * GV%m_to_H * &
                  (CS%bathyT(i,j+1) * CS%bathyT(i,j)) / &
                  (CS%bathyT(i,j+1) + CS%bathyT(i,j))
    enddo ; enddo
  endif

  if (rescale) then
!$OMP do
    do j=js-hs,je+hs ; do I=is-1-hs,ie+hs
      Datu(I,j) = Datu(I,j) * CS%Datu_res(I,j)
    enddo ; enddo
!$OMP do
    do J=js-1-hs,je+hs ; do i=is-hs,ie+hs
      Datv(i,J) = Datv(i,J) * CS%Datv_res(i,J)
    enddo ; enddo
  endif
!$OMP end parallel

end subroutine find_face_areas

subroutine legacy_bt_mass_source(h, eta, fluxes, set_cor, dt_therm, &
                                 dt_since_therm, G, GV, CS)
  type(ocean_grid_type),                intent(in) :: G
  type(verticalGrid_type),              intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h
  real, dimension(SZI_(G),SZJ_(G)),     intent(in) :: eta
  type(forcing),                        intent(in) :: fluxes
  logical,                              intent(in) :: set_cor
  real,                                 intent(in) :: dt_therm, dt_since_therm
  type(legacy_barotropic_CS),           pointer    :: CS
!   bt_mass_source determines the appropriately limited mass source for
! the barotropic solver, along with a corrective fictitious mass source that
! will drive the barotropic estimate of the free surface height toward the
! baroclinic estimate.

! Arguments: h - Layer thickness, in m or kg m-2 (H).
!  (in)      eta - The free surface height that is to be corrected, in m.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      set_cor - A flag to indicate whether to set the corrective fluxes
!                      (and update the slowly varying part of eta_cor) (.true.)
!                      or whether to incrementally update the corrective fluxes.
!  (in)      dt_therm - The thermodynamic time step, in s.
!  (in)      dt_since_therm - The elapsed time since mass forcing was applied, s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 barotropic_init.
  real :: h_tot(SZI_(G))      ! The sum of the layer thicknesses, in H.
  real :: eta_h(SZI_(G))      ! The free surface height determined from
                              ! the sum of the layer thicknesses, in H.
  real :: d_eta               ! The difference between estimates of the total
                              ! thicknesses, in H.
  real :: limit_dt            ! The fractional mass-source limit divided by the
                              ! thermodynamic time step, in s-1.
  integer :: is, ie, js, je, nz, i, j, k
  real, parameter :: frac_cor = 0.25
  real, parameter :: slow_rate = 0.125

  if (.not.associated(CS)) call MOM_error(FATAL, "bt_mass_source: "// &
        "Module MOM_barotropic must be initialized before it is used.")
  if (.not.CS%split) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

!$OMP parallel do default(none) shared(is,ie,js,je,nz,G,GV,h,set_cor,CS,dt_therm, &
!$OMP                                  fluxes,eta,dt_since_therm)              &
!$OMP                          private(eta_h,h_tot,limit_dt,d_eta)
  do j=js,je
    do i=is,ie ; h_tot(i) = h(i,j,1) ; enddo
    if (GV%Boussinesq) then
      do i=is,ie ; eta_h(i) = h(i,j,1) - G%bathyT(i,j) ; enddo
    else
      do i=is,ie ; eta_h(i) = h(i,j,1) ; enddo
    endif
    do k=2,nz ; do i=is,ie
      eta_h(i) = eta_h(i) + h(i,j,k)
      h_tot(i) = h_tot(i) + h(i,j,k)
    enddo ; enddo

    if (set_cor) then
      do i=is,ie ; CS%eta_source(i,j) = 0.0 ; enddo
      if (CS%eta_source_limit > 0.0) then
        limit_dt = CS%eta_source_limit/dt_therm
        if (associated(fluxes%lprec)) then ; do i=is,ie
          CS%eta_source(i,j) = CS%eta_source(i,j) + fluxes%lprec(i,j)
        enddo ; endif
        if (associated(fluxes%fprec)) then ; do i=is,ie
          CS%eta_source(i,j) = CS%eta_source(i,j) + fluxes%fprec(i,j)
        enddo ; endif
        if (associated(fluxes%vprec)) then ; do i=is,ie
          CS%eta_source(i,j) = CS%eta_source(i,j) + fluxes%vprec(i,j)
        enddo ; endif
        if (associated(fluxes%lrunoff)) then ; do i=is,ie
          CS%eta_source(i,j) = CS%eta_source(i,j) + fluxes%lrunoff(i,j)
        enddo ; endif
        if (associated(fluxes%frunoff)) then ; do i=is,ie
          CS%eta_source(i,j) = CS%eta_source(i,j) + fluxes%frunoff(i,j)
        enddo ; endif
        if (associated(fluxes%evap)) then ; do i=is,ie
          CS%eta_source(i,j) = CS%eta_source(i,j) + fluxes%evap(i,j)
        enddo ; endif
        do i=is,ie
          CS%eta_source(i,j) = CS%eta_source(i,j)*GV%kg_m2_to_H
          if (abs(CS%eta_source(i,j)) > limit_dt * h_tot(i)) then
            CS%eta_source(i,j) = SIGN(limit_dt * h_tot(i), CS%eta_source(i,j))
          endif
        enddo
      endif
    endif

    if (set_cor) then
      do i=is,ie
        d_eta = eta_h(i) - (eta(i,j) - dt_since_therm*CS%eta_source(i,j))
        CS%eta_cor(i,j) = d_eta
      enddo
    else
      do i=is,ie
        d_eta = eta_h(i) - (eta(i,j) - dt_since_therm*CS%eta_source(i,j))
        CS%eta_cor(i,j) = CS%eta_cor(i,j) + d_eta
      enddo
    endif
  enddo

end subroutine legacy_bt_mass_source

subroutine legacy_barotropic_init(u, v, h, eta, Time, G, GV, param_file, diag, CS, &
                           restart_CS, BT_cont, tides_CSp)
  type(ocean_grid_type),              intent(inout) :: G
  type(verticalGrid_type),            intent(in)    :: GV
  real, intent(in), dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: u
  real, intent(in), dimension(SZI_(G),SZJB_(G),SZK_(G)) :: v
  real, intent(in), dimension(SZI_(G),SZJ_(G),SZK_(G))  :: h
  real, intent(in), dimension(SZI_(G),SZJ_(G))      :: eta
  type(time_type), target,            intent(in)    :: Time
  type(param_file_type),              intent(in)    :: param_file
  type(diag_ctrl), target,            intent(inout) :: diag
  type(legacy_barotropic_CS),         pointer       :: CS
  type(MOM_restart_CS),               pointer       :: restart_CS
  type(BT_cont_type),       optional, pointer       :: BT_cont
  type(tidal_forcing_CS),   optional, pointer       :: tides_CSp
!   barotropic_init initializes a number of time-invariant fields used in the
! barotropic calculation and initializes any barotropic fields that have not
! already been initialized.

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      eta - Free surface height or column mass anomaly, in m or kg m-2.
!  (in)      Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer to the control structure for this module that is
!                 set in register_barotropic_restarts.
!  (in)      restart_CS - A pointer to the restart control structure.
!  (in,opt)  BT_cont - A structure with elements that describe the effective
!                      open face areas as a function of barotropic flow.
!  (in)      tides_CSp - a pointer to the control structure of the tide module.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_barotropic"  ! This module's name.
  real :: Datu(SZIBS_(G),SZJ_(G)), Datv(SZI_(G),SZJBS_(G))
  real :: gtot_estimate ! Summing GV%g_prime gives an upper-bound estimate for pbce.
  real :: SSH_extra     ! An estimate of how much higher SSH might get, for use
                        ! in calculating the safe external wave speed.
  real :: dtbt_input
  type(memory_size_type) :: MS
  logical :: apply_bt_drag, use_BT_cont_type
  character(len=48) :: thickness_units, flux_units
  character*(40) :: hvel_str
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: isdw, iedw, jsdw, jedw
  integer :: i, j, k
  integer :: wd_halos(2), bt_halo_sz
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  MS%isdw = G%isd ; MS%iedw = G%ied ; MS%jsdw = G%jsd ; MS%jedw = G%jed

  if (CS%module_is_initialized) then
    call MOM_error(WARNING, "barotropic_init called with a control structure "// &
                            "that has already been initialized.")
    return
  endif
  CS%module_is_initialized = .true.

  CS%diag => diag ; CS%Time => Time
  if (present(tides_CSp)) then
    if (associated(tides_CSp)) CS%tides_CSp => tides_CSp
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "SPLIT", CS%split, &
                 "Use the split time stepping if true.", default=.true.)
  if (.not.CS%split) return

  ! ### USE SOMETHING OTHER THAN MAXVEL FOR THIS...
  call get_param(param_file, mod, "BOUND_BT_CORRECTION", CS%bound_BT_corr, &
                 "If true, the corrective pseudo mass-fluxes into the \n"//&
                 "barotropic solver are limited to values that require \n"//&
                 "less than 0.1*MAXVEL to be accommodated.",default=.false.)
  call get_param(param_file, mod, "GRADUAL_BT_ICS", CS%gradual_BT_ICs, &
                 "If true, adjust the initial conditions for the \n"//&
                 "barotropic solver to the values from the layered \n"//&
                 "solution over a whole timestep instead of instantly. \n"//&
                 "This is a decent approximation to the inclusion of \n"//&
                 "sum(u dh_dt) while also correcting for truncation errors.", &
                 default=.false.)
  call get_param(param_file, mod, "BT_USE_WIDE_HALOS", CS%use_wide_halos, &
                 "If true, use wide halos and march in during the \n"//&
                 "barotropic time stepping for efficiency.", default=.true., &
                 layoutParam=.true.)
  call get_param(param_file, mod, "BTHALO", bt_halo_sz, &
                 "The minimum halo size for the barotropic solver.", default=0, &
                 layoutParam=.true.)
#ifdef STATIC_MEMORY_
  if ((bt_halo_sz > 0) .and. (bt_halo_sz /= BTHALO_)) call MOM_error(FATAL, &
      "barotropic_init: Run-time values of BTHALO must agree with the \n"//&
      "macro BTHALO_ with STATIC_MEMORY_.")
  wd_halos(1) = WHALOI_+NIHALO_ ; wd_halos(2) = WHALOJ_+NJHALO_
#else
  wd_halos(1) = bt_halo_sz; wd_halos(2) =  bt_halo_sz
#endif
  call log_param(param_file, mod, "!BT x-halo", wd_halos(1), &
                 "The barotropic x-halo size that is actually used.", &
                 layoutParam=.true.)
  call log_param(param_file, mod, "!BT y-halo", wd_halos(2), &
                 "The barotropic y-halo size that is actually used.", &
                 layoutParam=.true.)

  call get_param(param_file, mod, "USE_BT_CONT_TYPE", use_BT_cont_type, &
               "If true, use a structure with elements that describe \n"//&
               "effective face areas from the summed continuity solver \n"//&
               "as a function the barotropic flow in coupling between \n"//&
               "the barotropic and baroclinic flow.  This is only used \n"//&
               "if SPLIT is true. \n", default=.true.)
  call get_param(param_file, mod, "NONLINEAR_BT_CONTINUITY", &
                                CS%Nonlinear_continuity, &
                 "If true, use nonlinear transports in the barotropic \n"//&
                 "continuity equation.  This does not apply if \n"//&
                 "USE_BT_CONT_TYPE is true.", default=.false.)
  CS%Nonlin_cont_update_period = 1
  if (CS%Nonlinear_continuity) &
    call get_param(param_file, mod, "NONLIN_BT_CONT_UPDATE_PERIOD", &
                                  CS%Nonlin_cont_update_period, &
                 "If NONLINEAR_BT_CONTINUITY is true, this is the number \n"//&
                 "of barotropic time steps between updates to the face \n"//&
                 "areas, or 0 to update only before the barotropic stepping.",&
                 units="nondim", default=1)
  call get_param(param_file, mod, "RESCALE_BT_FACE_AREAS", CS%rescale_D_bt, &
                 "If true, the face areas used by the barotropic solver \n"//&
                 "are rescaled to approximately reflect the open face \n"//&
                 "areas of the interior layers.  This probably requires \n"//&
                 "FLUX_BT_COUPLING to work, and should not be used with \n"//&
                 "USE_BT_CONT_TYPE.", default=.false.)
  call get_param(param_file, mod, "BT_MASS_SOURCE_LIMIT", CS%eta_source_limit, &
                 "The fraction of the initial depth of the ocean that can \n"//&
                 "be added to or removed from the bartropic solution \n"//&
                 "within a thermodynamic time step.  By default this is 0 \n"//&
                 "for no correction.", units="nondim", default=0.0)
  call get_param(param_file, mod, "BT_PROJECT_VELOCITY", CS%BT_project_velocity,&
                 "If true, step the barotropic velocity first and project \n"//&
                 "out the velocity tendancy by 1+BEBT when calculating the \n"//&
                 "transport.  The default (false) is to use a predictor \n"//&
                 "continuity step to find the pressure field, and then \n"//&
                 "to do a corrector continuity step using a weighted \n"//&
                 "average of the old and new velocities, with weights \n"//&
                 "of (1-BEBT) and BEBT.", default=.false.)

  call get_param(param_file, mod, "DYNAMIC_SURFACE_PRESSURE", CS%dynamic_psurf, &
                 "If true, add a dynamic pressure due to a viscous ice \n"//&
                 "shelf, for instance.", default=.false.)
  if (CS%dynamic_psurf) then
    call get_param(param_file, mod, "ICE_LENGTH_DYN_PSURF", CS%ice_strength_length, &
                 "The length scale at which the Rayleigh damping rate due \n"//&
                 "to the ice strength should be the same as if a Laplacian \n"//&
                 "were applied, if DYNAMIC_SURFACE_PRESSURE is true.", &
                 units="m", default=1.0e4)
    call get_param(param_file, mod, "DEPTH_MIN_DYN_PSURF", CS%Dmin_dyn_psurf, &
                  "The minimum depth to use in limiting the size of the \n"//&
                  "dynamic surface pressure for stability, if \n"//&
                  "DYNAMIC_SURFACE_PRESSURE is true..", units="m", &
                  default=1.0e-6)
    call get_param(param_file, mod, "CONST_DYN_PSURF", CS%const_dyn_psurf, &
                 "The constant that scales the dynamic surface pressure, \n"//&
                 "if DYNAMIC_SURFACE_PRESSURE is true.  Stable values \n"//&
                 "are < ~1.0.", units="nondim", default=0.9)
  endif

  call get_param(param_file, mod, "TIDES", CS%tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  call get_param(param_file, mod, "SADOURNY", CS%Sadourny, &
                 "If true, the Coriolis terms are discretized with the \n"//&
                 "Sadourny (1975) energy conserving scheme, otherwise \n"//&
                 "the Arakawa & Hsu scheme is used.  If the internal \n"//&
                 "deformation radius is not resolved, the Sadourny scheme \n"//&
                 "should probably be used.", default=.true.)

  call get_param(param_file, mod, "BT_THICK_SCHEME", hvel_str, &
                 "A string describing the scheme that is used to set the \n"//&
                 "open face areas used for barotropic transport and the \n"//&
                 "relative weights of the accelerations. Valid values are:\n"//&
                 "\t ARITHMETIC - arithmetic mean layer thicknesses \n"//&
                 "\t HARMONIC - harmonic mean layer thicknesses \n"//&
                 "\t HYBRID (the default) - use arithmetic means for \n"//&
                 "\t    layers above the shallowest bottom, the harmonic \n"//&
                 "\t    mean for layers below, and a weighted average for \n"//&
                 "\t    layers that straddle that depth \n"//&
                 "\t FROM_BT_CONT - use the average thicknesses kept \n"//&
                 "\t    in the h_u and h_v fields of the BT_cont_type", &
                 default=BT_CONT_STRING)
  select case (hvel_str)
    case (HYBRID_STRING) ; CS%hvel_scheme = HYBRID
    case (HARMONIC_STRING) ; CS%hvel_scheme = HARMONIC
    case (ARITHMETIC_STRING) ; CS%hvel_scheme = ARITHMETIC
    case (BT_CONT_STRING) ; CS%hvel_scheme = FROM_BT_CONT
    case default
      call MOM_mesg('barotropic_init: BT_THICK_SCHEME ="'//trim(hvel_str)//'"', 0)
      call MOM_error(FATAL, "barotropic_init: Unrecognized setting "// &
            "#define BT_THICK_SCHEME "//trim(hvel_str)//" found in input file.")
  end select
  if ((CS%hvel_scheme == FROM_BT_CONT) .and. .not.use_BT_cont_type) &
    call MOM_error(FATAL, "barotropic_init: BT_THICK_SCHEME FROM_BT_CONT "//&
                           "can only be used if USE_BT_CONT_TYPE is defined.")

  call get_param(param_file, mod, "APPLY_BT_DRAG", apply_bt_drag, &
                 "If defined, bottom drag is applied within the \n"//&
                 "barotropic solver.", default=.true.)
  call get_param(param_file, mod, "BT_STRONG_DRAG", CS%strong_drag, &
                 "If true, use a stronger estimate of the retarding \n"//&
                 "effects of strong bottom drag, by making it implicit \n"//&
                 "with the barotropic time-step instead of implicit with \n"//&
                 "the baroclinic time-step and dividing by the number of \n"//&
                 "barotropic steps.", default=.true.)

  call get_param(param_file, mod, "CLIP_BT_VELOCITY", CS%clip_velocity, &
                 "If true, limit any velocity components that exceed \n"//&
                 "MAXVEL.  This should only be used as a desperate \n"//&
                 "debugging measure.", default=.false.)
  call get_param(param_file, mod, "MAXVEL", CS%maxvel, &
                 "The maximum velocity allowed before the velocity \n"//&
                 "components are truncated.", units="m s-1", default=3.0e8, &
                 do_not_log=.not.CS%clip_velocity)
  call get_param(param_file, mod, "MAXCFL_BT_CONT", CS%maxCFL_BT_cont, &
                 "The maximum permitted CFL number associated with the \n"//&
                 "barotropic accelerations from the summed velocities \n"//&
                 "times the time-derivatives of thicknesses.", units="nondim", &
                 default=0.1)

  call get_param(param_file, mod, "DT_BT_FILTER", CS%dt_bt_filter, &
                 "A time-scale over which the barotropic mode solutions \n"//&
                 "are filtered, in seconds if positive, or as a fraction \n"//&
                 "of DT if negative. When used this can never be taken to \n"//&
                 "be longer than 2*dt.  Set this to 0 to apply no filtering.", &
                 units="sec or nondim", default=-0.25)
  call get_param(param_file, mod, "G_BT_EXTRA", CS%G_extra, &
                 "A nondimensional factor by which gtot is enhanced.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mod, "SSH_EXTRA", SSH_extra, &
                 "An estimate of how much higher SSH might get, for use \n"//&
                 "in calculating the safe external wave speed. The \n"//&
                 "default is the minimum of 10 m or 5% of MAXIMUM_DEPTH.", &
                 units="m", default=min(10.0,0.05*G%max_depth))

  call get_param(param_file, mod, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)
  call get_param(param_file, mod, "DEBUG_BT", CS%debug_bt, &
                 "If true, write out verbose debugging data within the \n"//&
                 "barotropic time-stepping loop. The data volume can be \n"//&
                 "quite large if this is true.", default=CS%debug)

  CS%linearized_BT_PV = .true.
  call get_param(param_file, mod, "BEBT", CS%bebt, &
                 "BEBT determines whether the barotropic time stepping \n"//&
                 "uses the forward-backward time-stepping scheme or a \n"//&
                 "backward Euler scheme. BEBT is valid in the range from \n"//&
                 "0 (for a forward-backward treatment of nonrotating \n"//&
                 "gravity waves) to 1 (for a backward Euler treatment). \n"//&
                 "In practice, BEBT must be greater than about 0.05.", &
                 units="nondim", default=0.1)
  call get_param(param_file, mod, "DTBT", CS%dtbt, &
                 "The barotropic time step, in s. DTBT is only used with \n"//&
                 "the split explicit time stepping. To set the time step \n"//&
                 "automatically based the maximum stable value use 0, or \n"//&
                 "a negative value gives the fraction of the stable value. \n"//&
                 "Setting DTBT to 0 is the same as setting it to -0.98. \n"//&
                 "The value of DTBT that will actually be used is an \n"//&
                 "integer fraction of DT, rounding down.", units="s or nondim",&
                 default = -0.98)


  if (apply_bt_drag) then ; CS%drag_amp = 1.0 ; else ; CS%drag_amp = 0.0 ; endif

  ! Initialize a version of the MOM domain that is specific to the barotropic solver.
  call clone_MOM_domain(G%Domain, CS%BT_Domain, min_halo=wd_halos, symmetric=.true.)
#ifdef STATIC_MEMORY_
  if (wd_halos(1) /= WHALOI_+NIHALO_) call MOM_error(FATAL, "barotropic_init: "//&
          "Barotropic x-halo sizes are incorrectly resized with STATIC_MEMORY_.")
  if (wd_halos(2) /= WHALOJ_+NJHALO_) call MOM_error(FATAL, "barotropic_init: "//&
          "Barotropic y-halo sizes are incorrectly resized with STATIC_MEMORY_.")
#else
  if (bt_halo_sz > 0) then
    if (wd_halos(1) > bt_halo_sz) &
      call MOM_mesg("barotropic_init: barotropic x-halo size increased.", 3)
    if (wd_halos(2) > bt_halo_sz) &
      call MOM_mesg("barotropic_init: barotropic y-halo size increased.", 3)
  endif
#endif

  CS%isdw = G%isc-wd_halos(1) ; CS%iedw = G%iec+wd_halos(1)
  CS%jsdw = G%jsc-wd_halos(2) ; CS%jedw = G%jec+wd_halos(2)
  isdw = CS%isdw ; iedw = CS%iedw ; jsdw = CS%jsdw ; jedw = CS%jedw

  ALLOC_(CS%frhatu(IsdB:IedB,jsd:jed,nz)) ; ALLOC_(CS%frhatv(isd:ied,JsdB:JedB,nz))
  ALLOC_(CS%eta_source(isd:ied,jsd:jed)) ; ALLOC_(CS%eta_cor(isd:ied,jsd:jed))
  if (CS%bound_BT_corr) then
    ALLOC_(CS%eta_cor_bound(isd:ied,jsd:jed)) ; CS%eta_cor_bound(:,:) = 0.0
  endif
  ALLOC_(CS%IDatu(IsdB:IedB,jsd:jed)) ; ALLOC_(CS%IDatv(isd:ied,JsdB:JedB))

  ALLOC_(CS%Datu_res(isdw-1:iedw,jsdw:jedw))
  ALLOC_(CS%Datv_res(isdw:iedw,jsdw-1:jedw))
  ALLOC_(CS%ua_polarity(isdw:iedw,jsdw:jedw))
  ALLOC_(CS%va_polarity(isdw:iedw,jsdw:jedw))

  CS%frhatu(:,:,:) = 0.0 ; CS%frhatv(:,:,:) = 0.0
  CS%eta_source(:,:) = 0.0 ; CS%eta_cor(:,:) = 0.0
  CS%IDatu(:,:) = 0.0 ; CS%IDatv(:,:) = 0.0
  CS%Datu_res(:,:) = 1.0 ; CS%Datv_res(:,:) = 1.0

  CS%ua_polarity(:,:) = 1.0 ; CS%va_polarity(:,:) = 1.0
  call pass_vector(CS%ua_polarity, CS%va_polarity, CS%BT_domain, To_All, AGRID)

  if (use_BT_cont_type) &
    call alloc_BT_cont_type(BT_cont, G, (CS%hvel_scheme == FROM_BT_CONT))

  if (CS%debug) then ! Make a local copy of loop ranges for chksum calls
    allocate(CS%debug_BT_HI)
    CS%debug_BT_HI%isc=G%isc
    CS%debug_BT_HI%iec=G%iec
    CS%debug_BT_HI%jsc=G%jsc
    CS%debug_BT_HI%jec=G%jec
    CS%debug_BT_HI%IscB=G%isc-1
    CS%debug_BT_HI%IecB=G%iec
    CS%debug_BT_HI%JscB=G%jsc-1
    CS%debug_BT_HI%JecB=G%jec
    CS%debug_BT_HI%isd=CS%isdw
    CS%debug_BT_HI%ied=CS%iedw
    CS%debug_BT_HI%jsd=CS%jsdw
    CS%debug_BT_HI%jed=CS%jedw
    CS%debug_BT_HI%IsdB=CS%isdw-1
    CS%debug_BT_HI%IedB=CS%iedw
    CS%debug_BT_HI%JsdB=CS%jsdw-1
    CS%debug_BT_HI%JedB=CS%jedw

  endif

  ! IareaT, IdxCu, and IdyCv need to be allocated with wide halos.
  ALLOC_(CS%IareaT(CS%isdw:CS%iedw,CS%jsdw:CS%jedw)) ; CS%IareaT(:,:) = 0.0
  ALLOC_(CS%bathyT(CS%isdw:CS%iedw,CS%jsdw:CS%jedw)) ; CS%bathyT(:,:) = GV%Angstrom_z !### Should this be 0 instead?
  ALLOC_(CS%IdxCu(CS%isdw-1:CS%iedw,CS%jsdw:CS%jedw)) ; CS%IdxCu(:,:) = 0.0
  ALLOC_(CS%IdyCv(CS%isdw:CS%iedw,CS%jsdw-1:CS%jedw)) ; CS%IdyCv(:,:) = 0.0
  ALLOC_(CS%dy_Cu(CS%isdw-1:CS%iedw,CS%jsdw:CS%jedw)) ; CS%dy_Cu(:,:) = 0.0
  ALLOC_(CS%dx_Cv(CS%isdw:CS%iedw,CS%jsdw-1:CS%jedw)) ; CS%dx_Cv(:,:) = 0.0
  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    CS%IareaT(i,j) = G%IareaT(i,j)
    CS%bathyT(i,j) = G%bathyT(i,j)
  enddo ; enddo
  ! Note: G%IdxCu & G%IdyCv may be smaller than CS%IdxCu & CS%IdyCv, even without
  !   wide halos.
  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
    CS%IdxCu(I,j) = G%IdxCu(I,j) ; CS%dy_Cu(I,j) = G%dy_Cu(I,j)
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
    CS%IdyCv(I,j) = G%IdyCv(I,j) ; CS%dx_Cv(i,J) = G%dx_Cv(i,J)
  enddo ; enddo
  call pass_var(CS%IareaT, CS%BT_domain, To_All)
  call pass_var(CS%bathyT, CS%BT_domain, To_All)
  call pass_vector(CS%IdxCu, CS%IdyCv, CS%BT_domain, To_All+Scalar_Pair)
  call pass_vector(CS%dy_Cu, CS%dx_Cv, CS%BT_domain, To_All+Scalar_Pair)

  if (CS%linearized_BT_PV) then
    ALLOC_(CS%q_D(CS%isdw-1:CS%iedw,CS%jsdw-1:CS%jedw))
    ALLOC_(CS%D_u_Cor(CS%isdw-1:CS%iedw,CS%jsdw:CS%jedw))
    ALLOC_(CS%D_v_Cor(CS%isdw:CS%iedw,CS%jsdw-1:CS%jedw))
    CS%q_D(:,:) = 0.0 ; CS%D_u_Cor(:,:) = 0.0 ; CS%D_v_Cor(:,:) = 0.0
    do j=js,je ; do I=is-1,ie
      CS%D_u_Cor(I,j) = 0.5 * (G%bathyT(i+1,j) + G%bathyT(i,j))
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      CS%D_v_Cor(i,J) = 0.5 * (G%bathyT(i,j+1) + G%bathyT(i,j))
    enddo ; enddo
    do J=js-1,je ; do I=is-1,ie
      CS%q_D(I,J) = 0.25 * G%CoriolisBu(I,J) * &
           ((G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i+1,j) + G%areaT(i,j+1))) / &
           ((G%areaT(i,j) * G%bathyT(i,j) + G%areaT(i+1,j+1) * G%bathyT(i+1,j+1)) + &
            (G%areaT(i+1,j) * G%bathyT(i+1,j) + G%areaT(i,j+1) * G%bathyT(i,j+1)))
    enddo ; enddo
    ! With very wide halos, q and D need to be calculated on the available data
    ! domain and then updated onto the full computational domain.
    call pass_var(CS%q_D, CS%BT_Domain, To_All, position=CORNER)
    call pass_vector(CS%D_u_Cor, CS%D_v_Cor, CS%BT_Domain, To_All+Scalar_Pair)
  endif

  ! Estimate the maximum stable barotropic time step.
  dtbt_input = CS%dtbt
  CS%dtbt_fraction = 0.98 ; if (CS%dtbt < 0.0) CS%dtbt_fraction = -CS%dtbt
  gtot_estimate = 0.0
  do k=1,G%ke ; gtot_estimate = gtot_estimate + GV%g_prime(K) ; enddo
  call legacy_set_dtbt(G, GV, CS, gtot_est = gtot_estimate, SSH_add = SSH_extra)
  if (dtbt_input > 0.0) CS%dtbt = dtbt_input

  call log_param(param_file, mod, "!DTBT as used", CS%dtbt)
  call log_param(param_file, mod, "!estimated maximum DTBT", CS%dtbt_max)

  ! ubtav, vbtav, ubt_IC, vbt_IC, uhbt_IC, and vhbt_IC are allocated and
  ! initialized in register_barotropic_restarts.

  if (GV%Boussinesq) then
    thickness_units = "meter" ; flux_units = "meter3 second-1"
  else
    thickness_units = "kilogram meter-2" ; flux_units = "kilogram second-1"
  endif

  CS%id_PFu_bt = register_diag_field('ocean_model', 'PFuBT', diag%axesCu1, Time, &
      'Zonal Anomalous Barotropic Pressure Force Force Acceleration', 'meter second-2')
  CS%id_PFv_bt = register_diag_field('ocean_model', 'PFvBT', diag%axesCv1, Time, &
      'Meridional Anomalous Barotropic Pressure Force Acceleration', 'meter second-2')
  CS%id_Coru_bt = register_diag_field('ocean_model', 'CoruBT', diag%axesCu1, Time, &
      'Zonal Barotropic Coriolis Acceleration', 'meter second-2')
  CS%id_Corv_bt = register_diag_field('ocean_model', 'CorvBT', diag%axesCv1, Time, &
      'Meridional Barotropic Coriolis Acceleration', 'meter second-2')
  CS%id_uaccel = register_diag_field('ocean_model', 'u_accel_bt', diag%axesCu1, Time, &
      'Barotropic zonal acceleration', 'meter second-2')
  CS%id_vaccel = register_diag_field('ocean_model', 'v_accel_bt', diag%axesCv1, Time, &
      'Barotropic meridional acceleration', 'meter second-2')
  CS%id_ubtforce = register_diag_field('ocean_model', 'ubtforce', diag%axesCu1, Time, &
      'Barotropic zonal acceleration from baroclinic terms', 'meter second-2')
  CS%id_vbtforce = register_diag_field('ocean_model', 'vbtforce', diag%axesCv1, Time, &
      'Barotropic meridional acceleration from baroclinic terms', 'meter second-2')

  CS%id_eta_bt = register_diag_field('ocean_model', 'eta_bt', diag%axesT1, Time, &
      'Barotropic end SSH', thickness_units)
  CS%id_ubt = register_diag_field('ocean_model', 'ubt', diag%axesCu1, Time, &
      'Barotropic end zonal velocity', 'meter second-1')
  CS%id_vbt = register_diag_field('ocean_model', 'vbt', diag%axesCv1, Time, &
      'Barotropic end meridional velocity', 'meter second-1')
  CS%id_eta_st = register_diag_field('ocean_model', 'eta_st', diag%axesT1, Time, &
      'Barotropic start SSH', thickness_units)
  CS%id_ubt_st = register_diag_field('ocean_model', 'ubt_st', diag%axesCu1, Time, &
      'Barotropic start zonal velocity', 'meter second-1')
  CS%id_vbt_st = register_diag_field('ocean_model', 'vbt_st', diag%axesCv1, Time, &
      'Barotropic start meridional velocity', 'meter second-1')
  CS%id_ubtav = register_diag_field('ocean_model', 'ubtav', diag%axesCu1, Time, &
      'Barotropic time-average zonal velocity', 'meter second-1')
  CS%id_vbtav = register_diag_field('ocean_model', 'vbtav', diag%axesCv1, Time, &
      'Barotropic time-average meridional velocity', 'meter second-1')
  CS%id_eta_cor = register_diag_field('ocean_model', 'eta_cor', diag%axesT1, Time, &
      'Corrective mass flux', 'meter second-1')
  CS%id_visc_rem_u = register_diag_field('ocean_model', 'visc_rem_u', diag%axesCuL, Time, &
      'Viscous remnant at u', 'Nondim')
  CS%id_visc_rem_v = register_diag_field('ocean_model', 'visc_rem_v', diag%axesCvL, Time, &
      'Viscous remnant at v', 'Nondim')
  CS%id_gtotn = register_diag_field('ocean_model', 'gtot_n', diag%axesT1, Time, &
      'gtot to North', 'm s-2')
  CS%id_gtots = register_diag_field('ocean_model', 'gtot_s', diag%axesT1, Time, &
      'gtot to South', 'm s-2')
  CS%id_gtote = register_diag_field('ocean_model', 'gtot_e', diag%axesT1, Time, &
      'gtot to East', 'm s-2')
  CS%id_gtotw = register_diag_field('ocean_model', 'gtot_w', diag%axesT1, Time, &
      'gtot to West', 'm s-2')
  CS%id_eta_hifreq = register_diag_field('ocean_model', 'eta_hifreq', diag%axesT1, Time, &
      'High Frequency Barotropic SSH', thickness_units)
  CS%id_ubt_hifreq = register_diag_field('ocean_model', 'ubt_hifreq', diag%axesCu1, Time, &
      'High Frequency Barotropic zonal velocity', 'meter second-1')
  CS%id_vbt_hifreq = register_diag_field('ocean_model', 'vbt_hifreq', diag%axesCv1, Time, &
      'High Frequency Barotropic meridional velocity', 'meter second-1')
  CS%id_eta_pred_hifreq = register_diag_field('ocean_model', 'eta_pred_hifreq', diag%axesT1, Time, &
      'High Frequency Predictor Barotropic SSH', thickness_units)
  CS%id_uhbt_hifreq = register_diag_field('ocean_model', 'uhbt_hifreq', diag%axesCu1, Time, &
      'High Frequency Barotropic zonal transport', 'meter3 second-1')
  CS%id_vhbt_hifreq = register_diag_field('ocean_model', 'vhbt_hifreq', diag%axesCv1, Time, &
      'High Frequency Barotropic meridional transport', 'meter3 second-1')
  if (CS%rescale_D_bt) then
    CS%id_Datu_res = register_diag_field('ocean_model', 'Datu_res', diag%axesCu1, Time, &
      'Rescaling for zonal face area in barotropic continuity', 'Nondimensional')
    CS%id_Datv_res = register_diag_field('ocean_model', 'Datv_res', diag%axesCv1, Time, &
      'Rescaling for meridional face area in barotropic continuity', 'Nondimensional')
  endif
  CS%id_frhatu = register_diag_field('ocean_model', 'frhatu', diag%axesCuL, Time, &
      'Fractional thickness of layers in u-columns', 'Nondim')
  CS%id_frhatv = register_diag_field('ocean_model', 'frhatv', diag%axesCvL, Time, &
      'Fractional thickness of layers in v-columns', 'Nondim')
  CS%id_frhatu1 = register_diag_field('ocean_model', 'frhatu1', diag%axesCuL, Time, &
      'Predictor Fractional thickness of layers in u-columns', 'Nondim')
  CS%id_frhatv1 = register_diag_field('ocean_model', 'frhatv1', diag%axesCvL, Time, &
      'Predictor Fractional thickness of layers in v-columns', 'Nondim')
  CS%id_uhbt = register_diag_field('ocean_model', 'uhbt', diag%axesCu1, Time, &
      'Barotropic zonal transport averaged over a baroclinic step', 'meter3 second-1')
  CS%id_vhbt = register_diag_field('ocean_model', 'vhbt', diag%axesCv1, Time, &
      'Barotropic meridional transport averaged over a baroclinic step', 'meter3 second-1')

  if (CS%id_frhatu1 > 0) call safe_alloc_ptr(CS%frhatu1, IsdB,IedB,jsd,jed,nz)
  if (CS%id_frhatv1 > 0) call safe_alloc_ptr(CS%frhatv1, isd,ied,JsdB,JedB,nz)

  if (.NOT.query_initialized(CS%ubtav,"ubtav",restart_CS) .or. &
      .NOT.query_initialized(CS%vbtav,"vbtav",restart_CS)) then
    call legacy_btcalc(h, G, GV, CS, may_use_default=.true.)
    do j=js-1,je+1 ; do i=is-1,ie+1
      CS%ubtav(I,j) = 0.0 ; CS%vbtav(i,J) = 0.0
    enddo ; enddo
    do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      CS%ubtav(I,j) = CS%ubtav(I,j) + CS%frhatu(I,j,k) * u(I,j,k)
      CS%vbtav(I,j) = CS%vbtav(I,j) + CS%frhatv(i,J,k) * v(i,J,k)
    enddo ; enddo ; enddo
  endif

  if (.NOT.query_initialized(CS%ubt_IC,"ubt_IC",restart_CS) .or. &
      .NOT.query_initialized(CS%vbt_IC,"vbt_IC",restart_CS)) then
    do j=js,je ; do I=is-1,ie ; CS%ubt_IC(I,j) = CS%ubtav(I,j) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; CS%vbt_IC(i,J) = CS%vbtav(i,J) ; enddo ; enddo
  endif

!   Calculate other constants which are used for btstep.

  ! The following is only valid with the Boussinesq approximation.
! if (GV%Boussinesq) then
    do j=js,je ; do I=is-1,ie
      CS%IDatu(I,j) = G%mask2dCu(I,j) * 2.0 / (G%bathyT(i+1,j) + G%bathyT(i,j))
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      CS%IDatv(i,J) = G%mask2dCv(i,J) * 2.0 / (G%bathyT(i,j+1) + G%bathyT(i,j))
    enddo ; enddo
! else
!   do j=js,je ; do I=is-1,ie
!     CS%IDatu(I,j) = G%mask2dCu(I,j) * 2.0 / (GV%Rho0*(G%bathyT(i+1,j) + G%bathyT(i,j)))
!   enddo ; enddo
!   do J=js-1,je ; do i=is,ie
!     CS%IDatv(i,J) = G%mask2dCv(i,J) * 2.0 / (GV%Rho0*(G%bathyT(i,j+1) + G%bathyT(i,j)))
!   enddo ; enddo
! endif

  call find_face_areas(Datu, Datv, G, GV, CS, MS, halo=1)
  if (CS%bound_BT_corr) then
    do j=js,je ; do i=is,ie
      CS%eta_cor_bound(i,j) = GV%m_to_H * G%IareaT(i,j) * 0.1 * CS%maxvel * &
         ((Datu(I-1,j) + Datu(I,j)) + (Datv(i,J) + Datv(i,J-1)))
    enddo ; enddo
  endif

  if (.NOT.query_initialized(CS%uhbt_IC,"uhbt_IC",restart_CS) .or. &
      .NOT.query_initialized(CS%vhbt_IC,"vhbt_IC",restart_CS)) then
    do j=js,je ; do I=is-1,ie ; CS%uhbt_IC(I,j) = CS%ubtav(I,j) * Datu(I,j) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; CS%vhbt_IC(i,J) = CS%vbtav(i,J) * Datv(i,J) ; enddo ; enddo
  endif

  call pass_vector(CS%ubt_IC, CS%vbt_IC, G%Domain, complete=.false.)
  call pass_vector(CS%uhbt_IC, CS%vhbt_IC, G%Domain, complete=.false.)
  call pass_vector(CS%ubtav, CS%vbtav, G%Domain)


!  id_clock_pass = cpu_clock_id('(Ocean BT halo updates)', grain=CLOCK_ROUTINE)
  id_clock_calc_pre  = cpu_clock_id('(Ocean BT pre-calcs only)', grain=CLOCK_ROUTINE)
  id_clock_pass_pre = cpu_clock_id('(Ocean BT pre-step halo updates)', grain=CLOCK_ROUTINE)
  id_clock_calc = cpu_clock_id('(Ocean BT stepping calcs only)', grain=CLOCK_ROUTINE)
  id_clock_pass_step = cpu_clock_id('(Ocean BT stepping halo updates)', grain=CLOCK_ROUTINE)
  id_clock_calc_post = cpu_clock_id('(Ocean BT post-calcs only)', grain=CLOCK_ROUTINE)
  id_clock_pass_post = cpu_clock_id('(Ocean BT post-step halo updates)', grain=CLOCK_ROUTINE)
  if (dtbt_input <= 0.0) &
    id_clock_sync = cpu_clock_id('(Ocean BT global synch)', grain=CLOCK_ROUTINE)

end subroutine legacy_barotropic_init

subroutine legacy_barotropic_end(CS)
  type(legacy_barotropic_CS), pointer :: CS
  DEALLOC_(CS%frhatu)   ; DEALLOC_(CS%frhatv)
  DEALLOC_(CS%IDatu)    ; DEALLOC_(CS%IDatv)
  DEALLOC_(CS%Datu_res) ; DEALLOC_(CS%Datv_res)
  DEALLOC_(CS%ubtav)    ; DEALLOC_(CS%vbtav)
  DEALLOC_(CS%eta_cor)  ; DEALLOC_(CS%eta_source)
  DEALLOC_(CS%ua_polarity) ; DEALLOC_(CS%va_polarity)
  if (CS%bound_BT_corr) then
    DEALLOC_(CS%eta_cor_bound)
  endif

  deallocate(CS)
end subroutine legacy_barotropic_end

subroutine register_legacy_barotropic_restarts(HI, GV, param_file, CS, restart_CS)
  type(hor_index_type),    intent(in) :: HI
  type(verticalGrid_type), intent(in) :: GV
  type(param_file_type),   intent(in) :: param_file
  type(legacy_barotropic_CS), pointer :: CS
  type(MOM_restart_CS),    pointer    :: restart_CS
! This subroutine is used to register any fields from MOM_barotropic.F90
! that should be written to or read from the restart file.
! Arguments: HI - A horizontal index type structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in)      restart_CS - A pointer to the restart control structure.
  type(vardesc) :: vd(3)
  real :: slow_rate
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed
  IsdB = HI%IsdB ; IedB = HI%IedB ; JsdB = HI%JsdB ; JedB = HI%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "register_barotropic_restarts called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  ALLOC_(CS%ubtav(IsdB:IedB,jsd:jed))      ; CS%ubtav(:,:) = 0.0
  ALLOC_(CS%vbtav(isd:ied,JsdB:JedB))      ; CS%vbtav(:,:) = 0.0
  ALLOC_(CS%ubt_IC(IsdB:IedB,jsd:jed))     ; CS%ubt_IC(:,:) = 0.0
  ALLOC_(CS%vbt_IC(isd:ied,JsdB:JedB))     ; CS%vbt_IC(:,:) = 0.0
  ALLOC_(CS%uhbt_IC(IsdB:IedB,jsd:jed))    ; CS%uhbt_IC(:,:) = 0.0
  ALLOC_(CS%vhbt_IC(isd:ied,JsdB:JedB))    ; CS%vhbt_IC(:,:) = 0.0

  vd(2) = var_desc("ubtav","meter second-1","Time mean barotropic zonal velocity", &
                hor_grid='u', z_grid='1')
  vd(3) = var_desc("vbtav","meter second-1","Time mean barotropic meridional velocity",&
                hor_grid='v', z_grid='1')
  call register_restart_field(CS%ubtav, vd(2), .false., restart_CS)
  call register_restart_field(CS%vbtav, vd(3), .false., restart_CS)

  vd(2) = var_desc("ubt_IC", "meter second-1", &
              longname="Next initial condition for the barotropic zonal velocity", &
              hor_grid='u', z_grid='1')
  vd(3) = var_desc("vbt_IC", "meter second-1", &
              longname="Next initial condition for the barotropic meridional velocity",&
              hor_grid='v', z_grid='1')
  call register_restart_field(CS%ubt_IC, vd(2), .false., restart_CS)
  call register_restart_field(CS%vbt_IC, vd(3), .false., restart_CS)

  if (GV%Boussinesq) then
    vd(2) = var_desc("uhbt_IC", "meter3 second-1", &
                longname="Next initial condition for the barotropic zonal transport", &
                hor_grid='u', z_grid='1')
    vd(3) = var_desc("vhbt_IC", "meter3 second-1", &
                longname="Next initial condition for the barotropic meridional transport",&
                hor_grid='v', z_grid='1')
  else
    vd(2) = var_desc("uhbt_IC", "kg second-1", &
                longname="Next initial condition for the barotropic zonal transport", &
                hor_grid='u', z_grid='1')
    vd(3) = var_desc("vhbt_IC", "kg second-1", &
                longname="Next initial condition for the barotropic meridional transport",&
                hor_grid='v', z_grid='1')
  endif
  call register_restart_field(CS%uhbt_IC, vd(2), .false., restart_CS)
  call register_restart_field(CS%vhbt_IC, vd(3), .false., restart_CS)

end subroutine register_legacy_barotropic_restarts

end module MOM_legacy_barotropic
