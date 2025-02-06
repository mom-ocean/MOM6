!> Subroutines that use the ray-tracing equations to propagate the internal tide energy density.
!!
!! \author Benjamin Mater & Robert Hallberg, 2015
module MOM_internal_tides

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums,     only : hchksum
use MOM_debugging,     only : is_NaN
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_axis_init
use MOM_diag_mediator, only : disable_averaging, enable_averages
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_diag_mediator, only : axes_grp, define_axes_group
use MOM_domains, only       : AGRID, To_South, To_West, To_All, CGRID_NE
use MOM_domains, only       : create_group_pass, pass_var, pass_vector
use MOM_domains, only       : group_pass_type, start_group_pass, complete_group_pass
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only   : read_param, get_param, log_param, log_version, param_file_type
use MOM_forcing_type,only   : forcing
use MOM_grid, only          : ocean_grid_type
use MOM_int_tide_input, only: int_tide_input_CS, get_input_TKE, get_barotropic_tidal_vel
use MOM_io, only            : slasher, MOM_read_data, file_exists, axis_info
use MOM_io, only            : set_axis_info, get_axis_info, stdout
use MOM_restart, only       : register_restart_field, MOM_restart_CS, restart_init, save_restart
use MOM_restart, only       : lock_check, restart_registry_lock
use MOM_spatial_means, only : global_area_integral
use MOM_string_functions, only: extract_real
use MOM_time_manager, only  : time_type, time_type_to_real, operator(+), operator(/), operator(-)
use MOM_unit_scaling, only  : unit_scale_type
use MOM_variables, only     : surface, thermo_var_ptrs, vertvisc_type
use MOM_verticalGrid, only  : verticalGrid_type
use MOM_wave_speed, only    : wave_speeds, wave_speed_CS, wave_speed_init

implicit none ; private

#include <MOM_memory.h>

public propagate_int_tide, register_int_tide_restarts
public internal_tides_init, internal_tides_end
public get_lowmode_loss, get_lowmode_diffusivity

!> This control structure has parameters for the MOM_internal_tides module
type, public :: int_tide_CS ; private
  logical :: initialized = .false.   !< True if this control structure has been initialized.
  logical :: do_int_tides    !< If true, use the internal tide code.
  integer :: nFreq = 0       !< The number of internal tide frequency bands
  integer :: nMode = 1       !< The number of internal tide vertical modes
  integer :: nAngle = 24     !< The number of internal tide angular orientations
  integer :: energized_angle = -1 !< If positive, only this angular band is energized for debugging purposes
  real    :: uniform_test_cg !< Uniform group velocity of internal tide
                             !! for testing internal tides [L T-1 ~> m s-1]
  logical :: corner_adv      !< If true, use a corner advection rather than PPM.
  logical :: upwind_1st      !< If true, use a first-order upwind scheme.
  logical :: simple_2nd      !< If true, use a simple second order (arithmetic mean) interpolation
                             !! of the edge values instead of the higher order interpolation.
  logical :: vol_CFL         !< If true, use the ratio of the open face lengths to the tracer cell
                             !! areas when estimating CFL numbers.  Without aggress_adjust,
                             !! the default is false; it is always true with aggress_adjust.
  logical :: use_PPMang      !< If true, use PPM for advection of energy in angular space.
  logical :: update_Kd       !< If true, the scheme will modify the diffusivities seen by the dynamics
  logical :: apply_refraction  !< If false, skip refraction (for debugging)
  logical :: apply_propagation !< If False, do not propagate energy (for debugging)
  logical :: debug             !< If true, use debugging prints
  logical :: init_forcing_only !< if True, add TKE forcing only at first step (for debugging)
  logical :: force_posit_En    !< if True, remove subroundoff negative values (needs enhancement)
  logical :: add_tke_forcing = .true. !< Whether to add forcing, used by init_forcing_only

  real, allocatable, dimension(:,:) :: fraction_tidal_input
                        !< how the energy from one tidal component is distributed
                        !! over the various vertical modes, 2d in frequency and mode [nondim]
  real, allocatable, dimension(:,:) :: refl_angle
                        !< local coastline/ridge/shelf angles read from file [rad]
                        ! (could be in G control structure)
  real :: nullangle = -999.9 !< placeholder value in cells with no reflection [rad]
  real, allocatable, dimension(:,:) :: refl_pref
                        !< partial reflection coeff for each "coast cell" [nondim]
                        ! (could be in G control structure)
  logical, allocatable, dimension(:,:) :: refl_pref_logical
                        !< true if reflecting cell with partial reflection
                        ! (could be in G control structure)
  logical, allocatable, dimension(:,:) :: refl_dbl
                        !< identifies reflection cells where double reflection
                        !! is possible (i.e. ridge cells)
                        ! (could be in G control structure)
  real, allocatable, dimension(:,:) :: trans
                        !< partial transmission coeff for each "coast cell" [nondim]
  real, allocatable, dimension(:,:) :: residual
                        !< residual of reflection and transmission coeff for each "coast cell" [nondim]
  real, allocatable, dimension(:,:,:,:) :: cp
                        !< horizontal phase speed [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_leak_loss
                        !< energy lost due to misc background processes [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_quad_loss
                        !< energy lost due to quadratic bottom drag [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_Froude_loss
                        !< energy lost due to wave breaking [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: TKE_itidal_loss_fixed
                        !< Fixed part of the energy lost due to small-scale drag [H Z2 L-2 ~> kg m-2] here;
                        !! This will be multiplied by N and the squared near-bottom velocity (and by
                        !! the near-bottom density in non-Boussinesq mode) to get the energy losses
                        !! in [R Z4 H-1 L-2 ~> kg m-2 or m]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_itidal_loss
                        !< energy lost due to small-scale wave drag [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_residual_loss
                        !< internal tide energy loss due to the residual at slopes [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_slope_loss
                        !< internal tide energy loss due to the residual at slopes [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: TKE_input_glo_dt
                        !< The integrated energy input to the internal waves [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: TKE_leak_loss_glo_dt
                        !< Integrated energy lost due to misc background processes [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: TKE_quad_loss_glo_dt
                        !< Integrated energy lost due to quadratic bottom drag [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: TKE_Froude_loss_glo_dt
                        !< Integrated energy lost due to wave breaking [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: TKE_itidal_loss_glo_dt
                        !< energy lost due to small-scale wave drag [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable, dimension(:,:) :: TKE_residual_loss_glo_dt
                        !< internal tide energy loss due to the residual at slopes [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: error_mode
                        !< internal tide energy budget error for each mode [H Z2 L2 T-2 ~> m5 s-2 or J]
  real, allocatable, dimension(:,:) :: tot_leak_loss !< Energy loss rates due to misc background processes,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: tot_quad_loss !< Energy loss rates due to quadratic bottom drag,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: tot_itidal_loss !< Energy loss rates due to small-scale drag,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~>  m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: tot_Froude_loss !< Energy loss rates due to wave breaking,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: tot_residual_loss !< Energy loss rates due to residual on slopes,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:) :: tot_allprocesses_loss !< Energy loss rates due to all processes,
                        !! summed over angle, frequency and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:,:) :: w_struct !< Vertical structure of vertical velocity (normalized)
                        !! for each frequency and each mode [nondim]
  real, allocatable, dimension(:,:,:,:) :: u_struct !< Vertical structure of horizontal velocity (normalized and
                        !! divided by layer thicknesses) for each frequency and each mode [Z-1 ~> m-1]
  real, allocatable, dimension(:,:,:) :: u_struct_max !< Maximum of u_struct,
                        !! for each mode [Z-1 ~> m-1]
  real, allocatable, dimension(:,:,:) :: u_struct_bot !< Bottom value of u_struct,
                        !! for each mode [Z-1 ~> m-1]
  real, allocatable, dimension(:,:,:) :: int_w2 !< Vertical integral of w_struct squared,
                        !! for each mode [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:) :: int_U2 !< Vertical integral of u_struct squared,
                        !! for each mode [H Z-2 ~> m-1 or kg m-4]
  real, allocatable, dimension(:,:,:) :: int_N2w2 !< Depth-integrated Brunt Vaissalla freqency times
                        !! vertical profile squared, for each mode [H T-2 ~> m s-2 or kg m-2 s-2]
  real :: q_itides      !< fraction of local dissipation [nondim]
  real :: mixing_effic  !< mixing efficiency [nondim]
  real :: En_sum        !< global sum of energy for use in debugging, in MKS units [m5 s-2 or J]
  real :: En_underflow  !< A minuscule amount of energy [H Z2 T-2 ~> m3 s-2 or J m-2]
  integer :: En_restart_power !< A power factor of 2 by which to multiply the energy in restart [nondim]
  type(time_type), pointer :: Time => NULL() !< A pointer to the model's clock.
  character(len=200) :: inputdir !< directory to look for coastline angle file
  real, allocatable, dimension(:,:,:,:) :: decay_rate_2d !< rate at which internal tide energy is
                                                         !! lost to the interior ocean internal wave field
                                                         !! as a function of longitude, latitude, frequency
                                                         !! and vertical mode [T-1 ~> s-1].
  real :: cdrag         !< The bottom drag coefficient [nondim].
  real :: drag_min_depth !< The minimum total ocean thickness that will be used in the denominator
                        !! of the quadratic drag terms for internal tides when
                        !! INTERNAL_TIDE_QUAD_DRAG is true [H ~> m or kg m-2]
  real :: gamma_osborn  !< Mixing efficiency from Osborn 1980 [nondim]
  real :: Kd_min        !< The minimum diapycnal diffusivity. [L2 T-1 ~> m2 s-1]
  real :: max_TKE_to_Kd !< Maximum allowed value for TKE_to_kd [H Z2 T-3 ~> m3 s-3 or W m-2]
  real :: min_thick_layer_Kd !< minimum layer thickness allowed to use with TKE_to_kd [H ~> m or kg m-2]
  logical :: apply_background_drag
                        !< If true, apply a drag due to background processes as a sink.
  logical :: apply_bottom_drag
                        !< If true, apply a quadratic bottom drag as a sink.
  logical :: apply_wave_drag
                        !< If true, apply scattering due to small-scale roughness as a sink.
  logical :: apply_Froude_drag
                        !< If true, apply wave breaking as a sink.
  real :: En_check_tol  !< An energy density tolerance for flagging points with small negative
                        !! internal tide energy [H Z2 T-2 ~> m3 s-2 or J m-2]
  logical :: apply_residual_drag
                        !< If true, apply sink from residual term of reflection/transmission.
  logical :: use_2d_decay_rate
                        !< If true, use a spatially varying decay rate for each harmonic.
  real, allocatable :: En(:,:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,frequency,mode)
                        !! integrated within an angular and frequency band [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable :: En_ini_glo(:,:)
                        !< The internal wave energy density as a function of (frequency,mode) spatially
                        !! integrated within an angular and frequency band [H Z2 L2 T-2 ~> m5 s-2 or J]
                        !! only at the start of the routine (for diags)
  real, allocatable :: En_end_glo(:,:)
                        !< The internal wave energy density as a function of (frequency,mode) spatially
                        !! integrated within an angular and frequency band [H Z2 L2 T-2 ~> m5 s-2 or J]
                        !! only at the end of the routine (for diags)
  real, allocatable :: En_restart_mode1(:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,freq)
                        !! for mode 1 [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable :: En_restart_mode2(:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,freq)
                        !! for mode 2 [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable :: En_restart_mode3(:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,freq)
                        !! for mode 3 [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable :: En_restart_mode4(:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,freq)
                        !! for mode 4 [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, allocatable :: En_restart_mode5(:,:,:,:)
                        !< The internal wave energy density as a function of (i,j,angle,freq)
                        !! for mode 5 [H Z2 T-2 ~> m3 s-2 or J m-2]

  real, allocatable, dimension(:) :: frequency  !< The frequency of each band [T-1 ~> s-1].
  real :: Int_tide_decay_scale  !< vertical decay scale for St Laurent profile [Z ~> m]
  real :: Int_tide_decay_scale_slope  !< vertical decay scale for St Laurent profile on slopes [Z ~> m]

  type(wave_speed_CS) :: wave_speed  !< Wave speed control structure
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                        !! timing of diagnostic output.

  !>@{ Diag handles
  ! Diag handles relevant to all modes, frequencies, and angles
  integer :: id_cg1      = -1                 ! diagnostic handle for mode-1 speed
  integer, allocatable, dimension(:) :: id_cn ! diagnostic handle for all mode speeds
  integer :: id_tot_En = -1
  integer :: id_refl_pref = -1, id_refl_ang = -1, id_land_mask = -1
  integer :: id_trans = -1, id_residual = -1
  integer :: id_dx_Cv = -1, id_dy_Cu = -1
  ! Diag handles considering: sums over all modes, frequencies, and angles
  integer :: id_tot_leak_loss = -1, id_tot_quad_loss = -1, id_tot_itidal_loss = -1
  integer :: id_tot_Froude_loss = -1, id_tot_residual_loss = -1, id_tot_allprocesses_loss = -1
  ! Diag handles considering: all modes & frequencies; summed over angles
  integer, allocatable, dimension(:,:) :: &
             id_En_mode, &
             id_itidal_loss_mode, &
             id_leak_loss_mode, &
             id_quad_loss_mode, &
             id_Froude_loss_mode, &
             id_residual_loss_mode, &
             id_allprocesses_loss_mode, &
             id_itide_drag, &
             id_Ub_mode, &
             id_cp_mode
  ! Diag handles considering: all modes, frequencies, and angles
  integer, allocatable, dimension(:,:) :: &
             id_En_ang_mode, &
             id_itidal_loss_ang_mode
  integer, allocatable, dimension(:) :: &
             id_TKE_itidal_input, &
             id_Ustruct_mode, &
             id_Wstruct_mode, &
             id_int_w2_mode, &
             id_int_U2_mode, &
             id_int_N2w2_mode
  !>@}

end type int_tide_CS

!> A structure with the active energy loop bounds.
type :: loop_bounds_type ; private
  !>@{ The active loop bounds
  integer :: ish, ieh, jsh, jeh
  !>@}
end type loop_bounds_type

contains

!> Calls subroutines in this file that are needed to refract, propagate,
!! and dissipate energy density of the internal tide.
subroutine propagate_int_tide(h, tv, Nb, Rho_bot, dt, G, GV, US, inttide_input_CSp, CS)
  type(ocean_grid_type),            intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in)    :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),            intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)    :: tv !< Pointer to thermodynamic variables
                                                        !! (needed for wave structure).
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: Nb !< Near-bottom buoyancy frequency [T-1 ~> s-1].
                                                        !! In some cases the input values are used, but in
                                                        !! others this is set along with the wave speeds.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: Rho_bot !< Near-bottom density or the Boussinesq
                                                        !! reference density [R ~> kg m-3].
  real,                             intent(in)    :: dt !< Length of time over which to advance
                                                        !! the internal tides [T ~> s].
  type(int_tide_input_CS),          intent(in)    :: inttide_input_CSp !< Internal tide input control structure
  type(int_tide_CS),                intent(inout) :: CS !< Internal tide control structure

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),CS%nFreq) :: &
    TKE_itidal_input, & !< The energy input to the internal waves [H Z2 T-3 ~> m3 s-3 or W m-2].
    vel_btTide !< Barotropic velocity read from file [L T-1 ~> m s-1].

  real, dimension(SZI_(G),SZJ_(G),2) :: &
    test           ! A test unit vector used to determine grid rotation in halos [nondim]
  real, dimension(SZI_(G),SZJ_(G),CS%nMode) :: &
    cn             ! baroclinic internal gravity wave speeds for each mode [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),CS%nFreq,CS%nMode) :: &
    tot_En_mode, & ! energy summed over angles only [H Z2 T-2 ~> m3 s-2 or J m-2]
    Ub, &          ! near-bottom horizontal velocity of wave (modal) [L T-1 ~> m s-1]
    Umax           ! Maximum horizontal velocity of wave (modal) [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),CS%nFreq,CS%nMode) :: &
    drag_scale     ! bottom drag scale [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    tot_vel_btTide2, & ! [L2 T-2 ~> m2 s-2]
    tot_En, &      ! energy summed over angles, modes, frequencies [H Z2 T-2 ~> m3 s-2 or J m-2]
    tot_leak_loss, tot_quad_loss, tot_itidal_loss, tot_Froude_loss, tot_residual_loss, tot_allprocesses_loss, &
                   ! energy loss rates summed over angle, freq, and mode [H Z2 T-3 ~> m3 s-3 or W m-2]
    htot, &        ! The vertical sum of the layer thicknesses [H ~> m or kg m-2]
    itidal_loss_mode, & ! Energy lost due to small-scale wave drag, summed over angles [H Z2 T-3 ~> m3 s-3 or W m-2]
    leak_loss_mode, &
    quad_loss_mode, &
    Froude_loss_mode, &
    residual_loss_mode, &
    allprocesses_loss_mode  ! Total energy loss rates for a given mode and frequency (summed over
                            ! all angles) [H Z2 T-3 ~> m3 s-3 or W m-2]
  real :: frac_per_sector ! The inverse of the number of angular, modal and frequency bins [nondim]
  real :: f2       ! The squared Coriolis parameter interpolated to a tracer point [T-2 ~> s-2]
  real :: Kmag2    ! A squared horizontal wavenumber [L-2 ~> m-2]
  real :: I_D_here ! The inverse of the local water column thickness [H-1 ~> m-1 or m2 kg-1]
  real :: I_mass   ! The inverse of the local water mass [R-1 Z-1 ~> m2 kg-1]
  real :: I_dt     ! The inverse of the timestep [T-1 ~> s-1]
  real :: En_restart_factor ! A multiplicative factor of the form 2**En_restart_power [nondim]
  real :: I_En_restart_factor ! The inverse of the restart mult factor [nondim]
  real :: freq2    ! The frequency squared [T-2 ~> s-2]
  real :: PE_term  ! total potential energy of profile [R Z ~> kg m-2]
  real :: KE_term  ! total kinetic energy of profile [R Z ~> kg m-2]
  real :: U_mag    ! rescaled magnitude of horizontal profile [L Z T-1 ~> m2 s-1]
  real :: W0       ! rescaled magnitude of vertical profile [Z T-1 ~> m s-1]
  real :: c_phase  ! The phase speed [L T-1 ~> m s-1]
  real :: loss_rate  ! An energy loss rate [T-1 ~> s-1]
  real :: Fr2_max    ! The column maximum internal wave Froude number squared [nondim]
  real :: cn_subRO        ! A tiny wave speed to prevent division by zero [L T-1 ~> m s-1]
  real :: en_subRO        ! A tiny energy to prevent division by zero [H Z2 T-2 ~> m3 s-2 or J m-2]
  real :: En_a, En_b                                 ! Energies for time stepping [H Z2 T-2 ~> m3 s-2 or J m-2]
  real :: En_new, En_check                           ! Energies for debugging [H Z2 T-2 ~> m3 s-2 or J m-2]
  real :: En_sumtmp                                  ! Energies for debugging [H Z2 L2 T-2 ~> m5 s-2 or J]
  real :: En_initial, Delta_E_check                  ! Energies for debugging [H Z2 T-2 ~> m3 s-2 or J m-2]
  real :: TKE_Froude_loss_check, TKE_Froude_loss_tot ! Energy losses for debugging [H Z2 T-3 ~> m3 s-3 or W m-2]
  real :: HZ2_T2_to_J_m2                             ! unit conversion factor for Energy from internal units
                                                     ! to mks [T2 kg H-1 Z-2 s-2 ~> kg m-3 or 1]
  real :: J_m2_to_HZ2_T2                             ! unit conversion factor for Energy from mks to internal
                                                     ! units [H Z2 s2 T-2 kg-1 ~> m3 kg-1 or 1]
  character(len=160) :: mesg  ! The text of an error message
  integer :: En_halo_ij_stencil ! The halo size needed for energy advection
  integer :: a, m, fr, i, j, k, is, ie, js, je, isd, ied, jsd, jed, nAngle
  integer :: id_g, jd_g         ! global (decomp-invar) indices (for debugging)
  type(group_pass_type), save :: pass_test, pass_En
  type(time_type) :: time_end
  logical:: avg_enabled

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nAngle = CS%NAngle

  HZ2_T2_to_J_m2 = GV%H_to_kg_m2*(US%Z_to_m**2)*(US%s_to_T**2)
  J_m2_to_HZ2_T2 = GV%kg_m2_to_H*(US%m_to_Z**2)*(US%T_to_s**2)

  cn_subRO = 1e-30*US%m_s_to_L_T
  en_subRO = 1e-30*J_m2_to_HZ2_T2

  I_dt = 1.0 / dt
  En_restart_factor = 2**CS%En_restart_power
  I_En_restart_factor = 1.0 / En_restart_factor

  ! initialize local arrays
  TKE_itidal_input(:,:,:) = 0.
  vel_btTide(:,:,:) = 0.
  tot_vel_btTide2(:,:) = 0.
  drag_scale(:,:,:,:) = 0.
  Ub(:,:,:,:) = 0.
  Umax(:,:,:,:) = 0.

  cn(:,:,:) = 0.

  ! Rebuild energy density array from multiple restarts
  do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
    CS%En(i,j,a,fr,1) = CS%En_restart_mode1(i,j,a,fr) * I_En_restart_factor
  enddo ; enddo ; enddo ; enddo

  if (CS%nMode >= 2) then
    do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
      CS%En(i,j,a,fr,2) = CS%En_restart_mode2(i,j,a,fr) * I_En_restart_factor
    enddo ; enddo ; enddo ; enddo
  endif

  if (CS%nMode >= 3) then
    do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
      CS%En(i,j,a,fr,3) = CS%En_restart_mode3(i,j,a,fr) * I_En_restart_factor
    enddo ; enddo ; enddo ; enddo
  endif

  if (CS%nMode >= 4) then
    do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
      CS%En(i,j,a,fr,4) = CS%En_restart_mode4(i,j,a,fr) * I_En_restart_factor
    enddo ; enddo ; enddo ; enddo
  endif

  if (CS%nMode >= 5) then
    do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
      CS%En(i,j,a,fr,5) = CS%En_restart_mode5(i,j,a,fr) * I_En_restart_factor
    enddo ; enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    ! save initial energy for online budget
    do m=1,CS%nMode ; do fr=1,CS%nFreq
      En_sumtmp = 0.
      do a=1,CS%nAngle
        En_sumtmp = En_sumtmp + global_area_integral(CS%En(:,:,a,fr,m), G, tmp_scale=HZ2_T2_to_J_m2)
      enddo
      CS%En_ini_glo(fr,m) = En_sumtmp
    enddo ; enddo
  endif

  ! Set properties related to the internal tides, such as the wave speeds, storing some
  ! of them in the control structure for this module.
  if (CS%uniform_test_cg > 0.0) then
    do m=1,CS%nMode ; cn(:,:,m) = CS%uniform_test_cg ; enddo
  else
    call wave_speeds(h, tv, G, GV, US, CS%nMode, cn, CS%wave_speed, &
                     CS%w_struct, CS%u_struct, CS%u_struct_max, CS%u_struct_bot, &
                     Nb, CS%int_w2, CS%int_U2, CS%int_N2w2, halo_size=2)
                   ! The value of halo_size above would have to be larger if there were
                   ! not a halo update between the calls to propagate_x and propagate_y.
                   ! It can be 1 point smaller if teleport is not used.
  endif

  call pass_var(cn,G%Domain)

  if (CS%debug) then
    call hchksum(cn(:,:,1), "CN mode 1", G%HI, haloshift=0, unscale=US%L_to_m*US%s_to_T)
    call hchksum(CS%w_struct(:,:,:,1), "Wstruct mode 1", G%HI, haloshift=0)
    call hchksum(CS%u_struct(:,:,:,1), "Ustruct mode 1", G%HI, haloshift=0, unscale=US%m_to_Z)
    call hchksum(CS%u_struct_bot(:,:,1), "Ustruct_bot mode 1", G%HI, haloshift=0, unscale=US%m_to_Z)
    call hchksum(CS%u_struct_max(:,:,1), "Ustruct_max mode 1", G%HI, haloshift=0, unscale=US%m_to_Z)
    call hchksum(CS%int_w2(:,:,1),   "int_w2", G%HI, haloshift=0, unscale=GV%H_to_MKS)
    call hchksum(CS%int_U2(:,:,1),   "int_U2", G%HI, haloshift=0, unscale=GV%H_to_mks*US%m_to_Z**2)
    call hchksum(CS%int_N2w2(:,:,1), "int_N2w2", G%HI, haloshift=0, unscale=GV%H_to_mks*US%s_to_T**2)
  endif

  ! Set the wave speeds for the modes, using cg(n) ~ cg(1)/n.**********************
  ! This is wrong, of course, but it works reasonably in some cases.
  ! Uncomment if wave_speed is not used to calculate the true values (BDM).
  !do m=1,CS%nMode ; do j=js-2,je+2 ; do i=is-2,ie+2
  !  cn(i,j,m) = cn(i,j,1) / real(m)
  !enddo ; enddo ; enddo

  ! Add the forcing.***************************************************************

  if (CS%add_tke_forcing) then

    call get_input_TKE(G, TKE_itidal_input, CS%nFreq, inttide_input_CSp)

    if (CS%debug) then
      call hchksum(TKE_itidal_input(:,:,1), "TKE_itidal_input", G%HI, haloshift=0, &
                   unscale=GV%H_to_mks*(US%Z_to_m**2)*(US%s_to_T)**3)
      call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides bf input", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    endif

    if (CS%energized_angle <= 0) then
      frac_per_sector = 1.0 / real(CS%nAngle)
      do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
        f2 = 0.25*((G%Coriolis2Bu(I,J) + G%Coriolis2Bu(I-1,J-1)) + &
                   (G%Coriolis2Bu(I-1,J) + G%Coriolis2Bu(I,J-1)))
        if (CS%frequency(fr)**2 > f2) then
          CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) + (dt*frac_per_sector*(1.0-CS%q_itides) * &
                              CS%fraction_tidal_input(fr,m) * TKE_itidal_input(i,j,fr))
        else
          ! zero out input TKE value to get correct diagnostics
          TKE_itidal_input(i,j,fr) = 0.
        endif
      enddo ; enddo ; enddo ; enddo ; enddo
    elseif (CS%energized_angle <= CS%nAngle) then
      frac_per_sector = 1.0
      a = CS%energized_angle
      do m=1,CS%nMode ; do fr=1,CS%nFreq ; do j=js,je ; do i=is,ie
        f2 = 0.25*((G%Coriolis2Bu(I,J) + G%Coriolis2Bu(I-1,J-1)) + &
                   (G%Coriolis2Bu(I-1,J) + G%Coriolis2Bu(I,J-1)))
        if (CS%frequency(fr)**2 > f2) then
          CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) + (dt*frac_per_sector*(1.0-CS%q_itides) * &
                              CS%fraction_tidal_input(fr,m) * TKE_itidal_input(i,j,fr))
        else
          ! zero out input TKE value to get correct diagnostics
          TKE_itidal_input(i,j,fr) = 0.
        endif
      enddo ; enddo ; enddo ; enddo
    else
      call MOM_error(WARNING, "Internal tide energy is being put into a angular "//&
                              "band that does not exist.")
    endif
  endif ! add tke forcing

  if (CS%init_forcing_only) CS%add_tke_forcing=.false.

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides af input", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    ! save forcing for online budget
    do m=1,CS%nMode ; do fr=1,CS%nFreq
      En_sumtmp = 0.
      do a=1,CS%nAngle
        En_sumtmp = En_sumtmp + global_area_integral(dt*frac_per_sector*(1.0-CS%q_itides)* &
                                                     CS%fraction_tidal_input(fr,m)*TKE_itidal_input(:,:,fr), &
                                                     G, tmp_scale=HZ2_T2_to_J_m2)
      enddo
      CS%TKE_input_glo_dt(fr,m) = En_sumtmp
    enddo ; enddo
  endif

  ! Pass a test vector to check for grid rotation in the halo updates.
  do j=jsd,jed ; do i=isd,ied ; test(i,j,1) = 1.0 ; test(i,j,2) = 0.0 ; enddo ; enddo
  call create_group_pass(pass_test, test(:,:,1), test(:,:,2), G%domain, stagger=AGRID)
  call start_group_pass(pass_test, G%domain)

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides af halo", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'prop_int_tide: after forcing')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'prop_int_tide: after forcing', CS%En_sum
    enddo ; enddo
  endif

  ! Apply half the refraction.
  if (CS%apply_refraction) then
    do m=1,CS%nMode ; do fr=1,CS%nFreq
      call refract(CS%En(:,:,:,fr,m), cn(:,:,m), CS%frequency(fr), 0.5*dt, &
                   G, US, CS%nAngle, CS%use_PPMang)
    enddo ; enddo
  endif
  ! A this point, CS%En is only valid on the computational domain.

  if (CS%force_posit_En) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=jsd,jed ; do i=isd,ied
        if (CS%En(i,j,a,fr,m)<0.0) then
          CS%En(i,j,a,fr,m) = 0.0
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides af refr", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'prop_int_tide: after 1/2 refraction')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'prop_int_tide: after 1/2 refraction', CS%En_sum
    enddo ; enddo
    ! Check for En<0 - for debugging, delete later
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=js,je ; do i=is,ie
        if (CS%En(i,j,a,fr,m)<0.0) then
          id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
          write(mesg,*) 'After first refraction: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                        'En=', HZ2_T2_to_J_m2*CS%En(i,j,a,fr,m)
          call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg))
          !call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  call complete_group_pass(pass_test, G%domain)

  ! Set the halo size to work on, using similar logic to that used in propagate.  This may need
  ! to be adjusted depending on the advection scheme and whether teleport is used.
  if (CS%upwind_1st) then ; En_halo_ij_stencil = 2
  else ; En_halo_ij_stencil = 3 ; endif

  ! Rotate points in the halos as necessary.
  call correct_halo_rotation(CS%En, test, G, CS%nAngle, halo=En_halo_ij_stencil)

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides af halo R", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'prop_int_tide: after correct halo rotation')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'prop_int_tide: after correct halo rotation', CS%En_sum
    enddo ; enddo
  endif

  ! Propagate the waves.
  do m=1,CS%nMode ; do fr=1,CS%Nfreq

    ! initialize residual loss, will be computed in propagate
    CS%TKE_residual_loss(:,:,:,fr,m) = 0.
    CS%TKE_slope_loss(:,:,:,fr,m) = 0.

    if (CS%apply_propagation) then
      call propagate(CS%En(:,:,:,fr,m), cn(:,:,m), CS%frequency(fr), dt, &
                     G, GV, US, CS, CS%NAngle, CS%TKE_slope_loss(:,:,:,fr,m))
      endif
  enddo ; enddo

  if (CS%force_posit_En) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=jsd,jed ; do i=isd,ied
        if (CS%En(i,j,a,fr,m)<0.0) then
          CS%En(i,j,a,fr,m) = 0.0
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides af prop", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'prop_int_tide: after propagate')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'prop_int_tide: after propagate', CS%En_sum
    enddo ; enddo
    ! Check for En<0 - for debugging, delete later
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=js,je ; do i=is,ie
        if (CS%En(i,j,a,fr,m)<0.0) then
          id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset
          if (abs(CS%En(i,j,a,fr,m))>CS%En_check_tol) then ! only print if large
            write(mesg,*)  'After propagation: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                           'En=', HZ2_T2_to_J_m2*CS%En(i,j,a,fr,m)
            call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg))
            ! RD propagate produces very little negative energy (diff 2 large numbers), needs fix
            !call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
          endif
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  if (CS%apply_refraction) then
    ! Apply the other half of the refraction.
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call refract(CS%En(:,:,:,fr,m), cn(:,:,m), CS%frequency(fr), 0.5*dt, &
                   G, US, CS%NAngle, CS%use_PPMang)
    enddo ; enddo
    ! A this point, CS%En is only valid on the computational domain.
  endif

  if (CS%force_posit_En) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=jsd,jed ; do i=isd,ied
        if (CS%En(i,j,a,fr,m)<0.0) then
          CS%En(i,j,a,fr,m) = 0.0
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides af refr2", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'prop_int_tide: after 2/2 refraction')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'prop_int_tide: after 2/2 refraction', CS%En_sum
    enddo ; enddo
    ! Check for En<0 - for debugging, delete later
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=js,je ; do i=is,ie
        if (CS%En(i,j,a,fr,m)<0.0) then
          id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
          write(mesg,*) 'After second refraction: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                        'En=', HZ2_T2_to_J_m2*CS%En(i,j,a,fr,m)
          call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg))
          !call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  ! Apply various dissipation mechanisms.
  if (CS%apply_background_drag .or. CS%apply_bottom_drag &
      .or. CS%apply_wave_drag .or. CS%apply_Froude_drag &
      .or. (CS%id_tot_En > 0)) then
    tot_En(:,:) = 0.0
    tot_En_mode(:,:,:,:) = 0.0
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      do j=js,je ; do i=is,ie ; do a=1,CS%nAngle
        tot_En(i,j) = tot_En(i,j) + CS%En(i,j,a,fr,m)
        tot_En_mode(i,j,fr,m) = tot_En_mode(i,j,fr,m) + CS%En(i,j,a,fr,m)
      enddo ; enddo ; enddo
    enddo ; enddo
  endif

  ! Extract the energy for mixing due to misc. processes (background leakage)------
  if (CS%apply_background_drag) then
    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
      ! Calculate loss rate and apply loss over the time step ; apply the same drag timescale
      ! to each En component (technically not correct; fix later)
      En_b = CS%En(i,j,a,fr,m) ! save previous value
      En_a = CS%En(i,j,a,fr,m) / (1.0 + (dt * CS%decay_rate_2d(i,j,fr,m))) ! implicit update
      CS%TKE_leak_loss(i,j,a,fr,m) = (En_b - En_a) * I_dt ! compute exact loss rate [H Z2 T-3 ~> m3 s-3 or W m-2]
      CS%En(i,j,a,fr,m) = En_a ! update value
    enddo ; enddo ; enddo ; enddo ; enddo
  endif

  if (CS%force_posit_En) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=jsd,jed ; do i=isd,ied
        if (CS%En(i,j,a,fr,m)<0.0) then
          CS%En(i,j,a,fr,m) = 0.0
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides after leak", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'prop_int_tide: after background drag')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'prop_int_tide: after background drag', CS%En_sum
      call sum_En(G, GV, US, CS, CS%TKE_leak_loss(:,:,:,fr,m) * dt, 'prop_int_tide: loss after background drag')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'prop_int_tide: loss after background drag', CS%En_sum
    enddo ; enddo
    ! Check for En<0 - for debugging, delete later
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=js,je ; do i=is,ie
        if (CS%En(i,j,a,fr,m)<0.0) then
          id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
          write(mesg,*) 'After leak loss: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                        'En=', HZ2_T2_to_J_m2*CS%En(i,j,a,fr,m)
          call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg), all_print=.true.)
          !call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
    ! save loss term for online budget
    do m=1,CS%nMode ; do fr=1,CS%nFreq
      En_sumtmp = 0.
      do a=1,CS%nAngle
        En_sumtmp = En_sumtmp + global_area_integral(CS%TKE_leak_loss(:,:,a,fr,m)*dt, G, &
                                                     tmp_scale=HZ2_T2_to_J_m2)
      enddo
      CS%TKE_leak_loss_glo_dt(fr,m) = En_sumtmp
    enddo ; enddo
  endif

  ! Extract the energy for mixing due to bottom drag-------------------------------
  if (CS%apply_bottom_drag) then
    do j=jsd,jed ; do i=isd,ied ; htot(i,j) = 0.0 ; enddo ; enddo

    call get_barotropic_tidal_vel(G, vel_btTide, CS%nFreq, inttide_input_CSp)

    do fr=1,CS%Nfreq ; do j=jsd,jed ; do i=isd,ied
      tot_vel_btTide2(i,j) =  tot_vel_btTide2(i,j) + (vel_btTide(i,j,fr)**2)
    enddo ; enddo ; enddo

    do k=1,GV%ke ; do j=jsd,jed ; do i=isd,ied
      htot(i,j) = htot(i,j) + h(i,j,k)
    enddo ; enddo ; enddo
    if (GV%Boussinesq) then
      ! This is mathematically equivalent to the form in the option below, but they differ at roundoff.
      do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do j=jsd,jed ; do i=isd,ied
        I_D_here = 1.0 / (max(htot(i,j), CS%drag_min_depth))
        drag_scale(i,j,fr,m) = CS%cdrag * sqrt(max(0.0, US%L_to_Z**2*tot_vel_btTide2(i,j) + &
                             (tot_En_mode(i,j,fr,m) * I_D_here))) * GV%Z_to_H*I_D_here
      enddo ; enddo ; enddo ; enddo
    else
      do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do j=jsd,jed ; do i=isd,ied
        I_D_here = 1.0 / (max(htot(i,j), CS%drag_min_depth))
        I_mass = GV%RZ_to_H * I_D_here
        drag_scale(i,j,fr,m) = (CS%cdrag * (Rho_bot(i,j)*I_mass)) * &
                              sqrt(max(0.0, US%L_to_Z**2*tot_vel_btTide2(i,j) + &
                                            (tot_En_mode(i,j,fr,m) * I_D_here)))
      enddo ; enddo ; enddo ; enddo
    endif

    if (CS%debug) call hchksum(drag_scale(:,:,1,1), "dragscale", G%HI, haloshift=0, unscale=US%s_to_T)
    if (CS%debug) call hchksum(tot_vel_btTide2(:,:), "tot_vel_btTide2", G%HI, haloshift=0, unscale=US%L_T_to_m_s**2)

    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
      ! Calculate loss rate and apply loss over the time step ; apply the same drag timescale
      ! to each En component (technically not correct; fix later)
      En_b = CS%En(i,j,a,fr,m)
      En_a = CS%En(i,j,a,fr,m) / (1.0 + (dt * drag_scale(i,j,fr,m))) ! implicit update
      CS%TKE_quad_loss(i,j,a,fr,m)  = (En_b - En_a) * I_dt
      CS%En(i,j,a,fr,m) = En_a
    enddo ; enddo ; enddo ; enddo ; enddo
  endif

  if (CS%force_posit_En) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=jsd,jed ; do i=isd,ied
        if (CS%En(i,j,a,fr,m)<0.0) then
          CS%En(i,j,a,fr,m) = 0.0
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides after quad", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    ! save loss term for online budget
    do m=1,CS%nMode ; do fr=1,CS%nFreq
      En_sumtmp = 0.
      do a=1,CS%nAngle
        En_sumtmp = En_sumtmp + global_area_integral(CS%TKE_quad_loss(:,:,a,fr,m)*dt, G, &
                                                     tmp_scale=HZ2_T2_to_J_m2)
      enddo
      CS%TKE_quad_loss_glo_dt(fr,m) = En_sumtmp
    enddo ; enddo
    ! Check for En<0 - for debugging, delete later
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=js,je ; do i=is,ie
        if (CS%En(i,j,a,fr,m)<0.0) then
          id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
          write(mesg,*) 'After bottom loss: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                        'En=', HZ2_T2_to_J_m2*CS%En(i,j,a,fr,m)
          call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg), all_print=.true.)
          !call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  ! Extract the energy for mixing due to scattering (wave-drag)--------------------
  ! still need to allow a portion of the extracted energy to go to higher modes.
  ! First, find velocity profiles
  if (CS%apply_wave_drag .or. CS%apply_Froude_drag) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq

      ! compute near-bottom and max horizontal baroclinic velocity values at each point
      do j=js,je ; do i=is,ie
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging

        ! Calculate wavenumber magnitude
        freq2 = CS%frequency(fr)**2

        f2 = (0.25*(G%CoriolisBu(I,J) + G%CoriolisBu(max(I-1,1),max(J-1,1)) + &
                    G%CoriolisBu(I,max(J-1,1)) + G%CoriolisBu(max(I-1,1),J)))**2
        Kmag2 = (freq2 - f2) / ((cn(i,j,m)**2) + (cn_subRO**2))


        ! Back-calculate amplitude from energy equation
        if ( (G%mask2dT(i,j) > 0.5) .and. (freq2*Kmag2 > 0.0)) then
          ! Units here are [R Z ~> kg m-2]
          KE_term = 0.25*GV%H_to_RZ*( (((freq2 + f2) / (freq2*Kmag2))*US%L_to_Z**2*CS%int_U2(i,j,m)) + &
                                   CS%int_w2(i,j,m) )
          PE_term = 0.25*GV%H_to_RZ*( CS%int_N2w2(i,j,m) / freq2 )

          if (KE_term + PE_term > 0.0) then
            W0 = sqrt( GV%H_to_RZ * tot_En_mode(i,j,fr,m) / (KE_term + PE_term) )
          else
            !call MOM_error(WARNING, "MOM internal tides: KE + PE <= 0.0; setting to W0 to 0.0")
            W0 = 0.0
          endif

          U_mag = W0 * sqrt((freq2 + f2) / (2.0*freq2*Kmag2))
          ! scaled maximum tidal velocity
          Umax(i,j,fr,m) = abs(U_mag * CS%u_struct_max(i,j,m))
          ! scaled bottom tidal velocity
          Ub(i,j,fr,m) = abs(U_mag * CS%u_struct_bot(i,j,m))
        else
          Umax(i,j,fr,m) = 0.
          Ub(i,j,fr,m) = 0.
        endif

      enddo ; enddo ! i-loop, j-loop
    enddo ; enddo ! fr-loop, m-loop
  endif ! apply_wave or _Froude_drag (Ub or Umax needed)
  ! Finally, apply loss
  if (CS%apply_wave_drag) then
    ! Calculate loss rate and apply loss over the time step
    call itidal_lowmode_loss(G, GV, US, CS, Nb, Rho_bot, Ub, CS%En, CS%TKE_itidal_loss_fixed, &
                             CS%TKE_itidal_loss, dt, halo_size=0)
  endif

  if (CS%force_posit_En) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=jsd,jed ; do i=isd,ied
        if (CS%En(i,j,a,fr,m)<0.0) then
          CS%En(i,j,a,fr,m) = 0.0
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides after wave", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'prop_int_tide: before Froude drag')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'prop_int_tide: before Froude drag', CS%En_sum
    enddo ; enddo
    ! save loss term for online budget, may want to add a debug flag later
    do m=1,CS%nMode ; do fr=1,CS%nFreq
      En_sumtmp = 0.
      do a=1,CS%nAngle
        En_sumtmp = En_sumtmp + global_area_integral(CS%TKE_itidal_loss(:,:,a,fr,m)*dt, G, &
                                                     tmp_scale=HZ2_T2_to_J_m2)
      enddo
      CS%TKE_itidal_loss_glo_dt(fr,m) = En_sumtmp
    enddo ; enddo
    ! Check for En<0 - for debugging, delete later
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=js,je ; do i=is,ie
        if (CS%En(i,j,a,fr,m)<0.0) then
          id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
          write(mesg,*) 'After wave drag loss: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                        'En=', HZ2_T2_to_J_m2*CS%En(i,j,a,fr,m)
          call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg), all_print=.true.)
          !call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  ! Extract the energy for mixing due to wave breaking-----------------------------
  if (CS%apply_Froude_drag) then
    ! Pick out maximum baroclinic velocity values; calculate Fr=max(u)/cg
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      freq2 = CS%frequency(fr)**2
      do j=js,je ; do i=is,ie
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
        ! Calculate horizontal phase velocity magnitudes
        f2 = 0.25*((G%Coriolis2Bu(I,J) + G%Coriolis2Bu(I-1,J-1)) + &
                   (G%Coriolis2Bu(I-1,J) + G%Coriolis2Bu(I,J-1)))
        Kmag2 = (freq2 - f2) / ((cn(i,j,m)**2) + (cn_subRO**2))
        c_phase = 0.0
        CS%TKE_Froude_loss(i,j,:,fr,m) = 0. ! init for all angles
        if (Kmag2 > 0.0) then
          c_phase = sqrt(freq2/Kmag2)
          Fr2_max = (Umax(i,j,fr,m) / c_phase)**2
          ! Dissipate energy if Fr>1; done here with an arbitrary time scale
          if (Fr2_max > 1.0) then
            ! Calculate effective decay rate [T-1 ~> s-1] if breaking occurs over a time step
            !loss_rate = (1.0 - Fr2_max) / (Fr2_max * dt)
            do a=1,CS%nAngle
              ! Determine effective dissipation rate (Wm-2)
              !CS%TKE_Froude_loss(i,j,a,fr,m) = CS%En(i,j,a,fr,m) * abs(loss_rate)
              En_b = CS%En(i,j,a,fr,m)
              En_a = CS%En(i,j,a,fr,m)/Fr2_max
              CS%TKE_Froude_loss(i,j,a,fr,m) = (En_b - En_a) * I_dt
              CS%En(i,j,a,fr,m) = En_a
            enddo
          endif ! Fr2>1
        endif ! Kmag2>0
        CS%cp(i,j,fr,m) = c_phase
      enddo ; enddo
    enddo ; enddo
  endif

  if (CS%force_posit_En) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=jsd,jed ; do i=isd,ied
        if (CS%En(i,j,a,fr,m)<0.0) then
          CS%En(i,j,a,fr,m) = 0.0
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides after froude", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'prop_int_tide: after Froude drag')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'prop_int_tide: after Froude drag', CS%En_sum
      call sum_En(G, GV, US, CS, CS%TKE_Froude_loss(:,:,:,fr,m) * dt, 'prop_int_tide: loss after Froude drag')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'prop_int_tide: loss after Froude drag', CS%En_sum
    enddo ; enddo
    ! save loss term for online budget, may want to add a debug flag later
    do m=1,CS%nMode ; do fr=1,CS%nFreq
      En_sumtmp = 0.
      do a=1,CS%nAngle
        En_sumtmp = En_sumtmp + global_area_integral(CS%TKE_Froude_loss(:,:,a,fr,m)*dt, G, &
                                                     tmp_scale=HZ2_T2_to_J_m2)
      enddo
      CS%TKE_Froude_loss_glo_dt(fr,m) = En_sumtmp
    enddo ; enddo
    ! Check for En<0 - for debugging, delete later
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=js,je ; do i=is,ie
        if (CS%En(i,j,a,fr,m)<0.0) then
          id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset
          write(mesg,*) 'After Froude loss: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                        'En=', HZ2_T2_to_J_m2*CS%En(i,j,a,fr,m)
          call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg), all_print=.true.)
          !call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
        endif
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  ! loss from residual of reflection/transmission coefficients
  if (CS%apply_residual_drag) then

    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
      ! implicit form is rewritten to minimize number of divisions
      !CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) / (1.0 + dt * CS%TKE_residual_loss(i,j,a,fr,m) / &
      !                                        (CS%En(i,j,a,fr,m) + en_subRO))
      ! only compute when partial reflection is present not to create noise elsewhere
      if (CS%refl_pref_logical(i,j)) then
        En_b = CS%En(i,j,a,fr,m)
        En_a = (CS%En(i,j,a,fr,m) * (CS%En(i,j,a,fr,m) + en_subRO)) / &
               ((CS%En(i,j,a,fr,m) + en_subRO) + (dt * CS%TKE_slope_loss(i,j,a,fr,m)))
        CS%TKE_residual_loss(i,j,a,fr,m) = (En_b - En_a) * I_dt
        CS%En(i,j,a,fr,m) = En_a
      endif
    enddo ; enddo ; enddo ; enddo ; enddo

  else
    ! zero out the residual loss term so it does not count towards diagnostics
    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
       CS%TKE_residual_loss(i,j,a,fr,m) = 0.
    enddo ; enddo ; enddo ; enddo ; enddo
  endif

  if (CS%debug) then
    call hchksum(CS%En(:,:,:,1,1), "EnergyIntTides after slope", G%HI, haloshift=0, unscale=HZ2_T2_to_J_m2)
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'prop_int_tide')
    enddo ; enddo
    ! save loss term for online budget
    do m=1,CS%nMode ; do fr=1,CS%nFreq
      En_sumtmp = 0.
      do a=1,CS%nAngle
        En_sumtmp = En_sumtmp + global_area_integral(CS%TKE_residual_loss(:,:,a,fr,m)*dt, G, &
                                                     tmp_scale=HZ2_T2_to_J_m2)
      enddo
      CS%TKE_residual_loss_glo_dt(fr,m) = En_sumtmp
    enddo ; enddo
  endif

  !---- energy budget ----
  if (CS%debug) then
    ! save final energy for online budget
    do m=1,CS%nMode ; do fr=1,CS%nFreq
      En_sumtmp = 0.
      do a=1,CS%nAngle
        En_sumtmp = En_sumtmp + global_area_integral(CS%En(:,:,a,fr,m), G, tmp_scale=HZ2_T2_to_J_m2)
      enddo
      CS%En_end_glo(fr,m) = En_sumtmp
    enddo ; enddo

    do m=1,CS%nMode ; do fr=1,CS%nFreq
      CS%error_mode(fr,m) = CS%En_ini_glo(fr,m) + CS%TKE_input_glo_dt(fr,m) - CS%TKE_leak_loss_glo_dt(fr,m) - &
                            CS%TKE_quad_loss_glo_dt(fr,m) - CS%TKE_itidal_loss_glo_dt(fr,m) - &
                            CS%TKE_Froude_loss_glo_dt(fr,m) - CS%TKE_residual_loss_glo_dt(fr,m) - &
                            CS%En_end_glo(fr,m)
      if (is_root_pe()) write(stdout,'(A,F18.10)') &
          "error in Energy budget", US%L_to_m**2*HZ2_T2_to_J_m2*CS%error_mode(fr,m)
    enddo ; enddo
  endif

  ! Output diagnostics.************************************************************
  avg_enabled = query_averaging_enabled(CS%diag, time_end=time_end)
  call enable_averages(dt, time_end, CS%diag)

  if (query_averaging_enabled(CS%diag)) then
    ! Output internal wave modal wave speeds
    if (CS%id_cg1 > 0) call post_data(CS%id_cg1, cn(:,:,1),CS%diag)
    do m=1,CS%nMode ; if (CS%id_cn(m) > 0) call post_data(CS%id_cn(m), cn(:,:,m), CS%diag) ; enddo

    ! Output two-dimensional diagnostics
    if (CS%id_tot_En > 0)     call post_data(CS%id_tot_En, tot_En, CS%diag)
    do fr=1,CS%nFreq
      if (CS%id_TKE_itidal_input(fr) > 0) call post_data(CS%id_TKE_itidal_input(fr), &
                                                         TKE_itidal_input(:,:,fr), CS%diag)
    enddo

    do m=1,CS%nMode ; do fr=1,CS%nFreq
      if (CS%id_itide_drag(fr,m) > 0) call post_data(CS%id_itide_drag(fr,m), drag_scale(:,:,fr,m), CS%diag)
    enddo ; enddo

    ! Output 2-D energy density (summed over angles) for each frequency and mode
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; if (CS%id_En_mode(fr,m) > 0) then
      tot_En(:,:) = 0.0
      do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
        tot_En(i,j) = tot_En(i,j) + CS%En(i,j,a,fr,m)
      enddo ; enddo ; enddo
      call post_data(CS%id_En_mode(fr,m), tot_En, CS%diag)
    endif ; enddo ; enddo

    ! split energy array into multiple restarts
    do fr=1,CS%Nfreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
      CS%En_restart_mode1(i,j,a,fr) = CS%En(i,j,a,fr,1) * En_restart_factor
    enddo ; enddo ; enddo ; enddo

    if (CS%nMode >= 2) then
      do fr=1,CS%Nfreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
        CS%En_restart_mode2(i,j,a,fr) = CS%En(i,j,a,fr,2) * En_restart_factor
      enddo ; enddo ; enddo ; enddo
    endif

    if (CS%nMode >= 3) then
      do fr=1,CS%Nfreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
        CS%En_restart_mode3(i,j,a,fr) = CS%En(i,j,a,fr,3) * En_restart_factor
      enddo ; enddo ; enddo ; enddo
    endif

    if (CS%nMode >= 4) then
      do fr=1,CS%Nfreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
        CS%En_restart_mode4(i,j,a,fr) = CS%En(i,j,a,fr,4) * En_restart_factor
      enddo ; enddo ; enddo ; enddo
    endif

    if (CS%nMode >= 5) then
      do fr=1,CS%Nfreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
        CS%En_restart_mode5(i,j,a,fr) = CS%En(i,j,a,fr,5) * En_restart_factor
      enddo ; enddo ; enddo ; enddo
    endif

    ! Output 3-D (i,j,a) energy density for each frequency and mode
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; if (CS%id_En_ang_mode(fr,m) > 0) then
      call post_data(CS%id_En_ang_mode(fr,m), CS%En(:,:,:,fr,m) , CS%diag)
    endif ; enddo ; enddo

    ! Output 2-D energy loss (summed over angles, freq, modes)
    tot_leak_loss(:,:)   = 0.0
    tot_quad_loss(:,:)   = 0.0
    tot_itidal_loss(:,:) = 0.0
    tot_Froude_loss(:,:) = 0.0
    tot_residual_loss(:,:) = 0.0
    tot_allprocesses_loss(:,:) = 0.0
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
      tot_leak_loss(i,j)   = tot_leak_loss(i,j)   + CS%TKE_leak_loss(i,j,a,fr,m)
      tot_quad_loss(i,j)   = tot_quad_loss(i,j)   + CS%TKE_quad_loss(i,j,a,fr,m)
      tot_itidal_loss(i,j) = tot_itidal_loss(i,j) + CS%TKE_itidal_loss(i,j,a,fr,m)
      tot_Froude_loss(i,j) = tot_Froude_loss(i,j) + CS%TKE_Froude_loss(i,j,a,fr,m)
      tot_residual_loss(i,j) = tot_residual_loss(i,j) + CS%TKE_residual_loss(i,j,a,fr,m)
    enddo ; enddo ; enddo ; enddo ; enddo
    do j=js,je ; do i=is,ie
      tot_allprocesses_loss(i,j) = ((((tot_leak_loss(i,j) + tot_quad_loss(i,j)) + &
                          tot_itidal_loss(i,j)) + tot_Froude_loss(i,j)) + &
                          tot_residual_loss(i,j))
    enddo ; enddo
    CS%tot_leak_loss         = tot_leak_loss
    CS%tot_quad_loss         = tot_quad_loss
    CS%tot_itidal_loss       = tot_itidal_loss
    CS%tot_Froude_loss       = tot_Froude_loss
    CS%tot_residual_loss     = tot_residual_loss
    CS%tot_allprocesses_loss = tot_allprocesses_loss
    if (CS%id_tot_leak_loss > 0) then
      call post_data(CS%id_tot_leak_loss, tot_leak_loss, CS%diag)
    endif
    if (CS%id_tot_quad_loss > 0) then
      call post_data(CS%id_tot_quad_loss, tot_quad_loss, CS%diag)
    endif
    if (CS%id_tot_itidal_loss > 0) then
      call post_data(CS%id_tot_itidal_loss, tot_itidal_loss, CS%diag)
    endif
    if (CS%id_tot_Froude_loss > 0) then
      call post_data(CS%id_tot_Froude_loss, tot_Froude_loss, CS%diag)
    endif
    if (CS%id_tot_residual_loss > 0) then
      call post_data(CS%id_tot_residual_loss, tot_residual_loss, CS%diag)
    endif
    if (CS%id_tot_allprocesses_loss > 0) then
      call post_data(CS%id_tot_allprocesses_loss, tot_allprocesses_loss, CS%diag)
    endif

    ! Output 2-D energy loss (summed over angles) for each frequency and mode
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
    if (CS%id_itidal_loss_mode(fr,m) > 0 .or. CS%id_allprocesses_loss_mode(fr,m) > 0) then
      itidal_loss_mode(:,:) = 0.0 ! wave-drag processes (could do others as well)
      leak_loss_mode(:,:) = 0.0
      quad_loss_mode(:,:) = 0.0
      Froude_loss_mode(:,:) = 0.0
      residual_loss_mode(:,:) = 0.0
      allprocesses_loss_mode(:,:) = 0.0 ! all processes summed together
      do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
        itidal_loss_mode(i,j)       = itidal_loss_mode(i,j) + CS%TKE_itidal_loss(i,j,a,fr,m)
        leak_loss_mode(i,j) =  leak_loss_mode(i,j) + CS%TKE_leak_loss(i,j,a,fr,m)
        quad_loss_mode(i,j) = quad_loss_mode(i,j) + CS%TKE_quad_loss(i,j,a,fr,m)
        Froude_loss_mode(i,j) = Froude_loss_mode(i,j) + CS%TKE_Froude_loss(i,j,a,fr,m)
        residual_loss_mode(i,j) = residual_loss_mode(i,j) +  CS%TKE_residual_loss(i,j,a,fr,m)
        allprocesses_loss_mode(i,j) = allprocesses_loss_mode(i,j) + &
                                 ((((CS%TKE_leak_loss(i,j,a,fr,m) + CS%TKE_quad_loss(i,j,a,fr,m)) + &
                                 CS%TKE_itidal_loss(i,j,a,fr,m)) + CS%TKE_Froude_loss(i,j,a,fr,m)) + &
                                 CS%TKE_residual_loss(i,j,a,fr,m))
      enddo ; enddo ; enddo
      call post_data(CS%id_itidal_loss_mode(fr,m), itidal_loss_mode, CS%diag)
      call post_data(CS%id_leak_loss_mode(fr,m), leak_loss_mode, CS%diag)
      call post_data(CS%id_quad_loss_mode(fr,m), quad_loss_mode, CS%diag)
      call post_data(CS%id_Froude_loss_mode(fr,m), Froude_loss_mode, CS%diag)
      call post_data(CS%id_residual_loss_mode(fr,m), residual_loss_mode, CS%diag)
      call post_data(CS%id_allprocesses_loss_mode(fr,m), allprocesses_loss_mode, CS%diag)
    endif ; enddo ; enddo

    ! Output 3-D (i,j,a) energy loss for each frequency and mode
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; if (CS%id_itidal_loss_ang_mode(fr,m) > 0) then
      call post_data(CS%id_itidal_loss_ang_mode(fr,m), CS%TKE_itidal_loss(:,:,:,fr,m) , CS%diag)
    endif ; enddo ; enddo

    ! Output 2-D period-averaged horizontal near-bottom mode velocity for each frequency and mode
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; if (CS%id_Ub_mode(fr,m) > 0) then
      call post_data(CS%id_Ub_mode(fr,m), Ub(:,:,fr,m), CS%diag)
    endif ; enddo ; enddo

    do m=1,CS%nMode ; if (CS%id_Ustruct_mode(m) > 0) then
      call post_data(CS%id_Ustruct_mode(m), CS%u_struct(:,:,:,m), CS%diag)
    endif ; enddo

    do m=1,CS%nMode ; if (CS%id_Wstruct_mode(m) > 0) then
      call post_data(CS%id_Wstruct_mode(m), CS%w_struct(:,:,:,m), CS%diag)
    endif ; enddo

    do m=1,CS%nMode ; if (CS%id_int_w2_mode(m) > 0) then
      call post_data(CS%id_int_w2_mode(m), CS%int_w2(:,:,m), CS%diag)
    endif ; enddo

    do m=1,CS%nMode ; if (CS%id_int_U2_mode(m) > 0) then
      call post_data(CS%id_int_U2_mode(m), CS%int_U2(:,:,m), CS%diag)
    endif ; enddo

    do m=1,CS%nMode ; if (CS%id_int_N2w2_mode(m) > 0) then
      call post_data(CS%id_int_N2w2_mode(m), CS%int_N2w2(:,:,m), CS%diag)
    endif ; enddo

    ! Output 2-D horizontal phase velocity for each frequency and mode
    do m=1,CS%nMode ; do fr=1,CS%Nfreq ; if (CS%id_cp_mode(fr,m) > 0) then
      call post_data(CS%id_cp_mode(fr,m), CS%cp(:,:,fr,m), CS%diag)
    endif ; enddo ; enddo

  endif

  call disable_averaging(CS%diag)

end subroutine propagate_int_tide

!> Checks for energy conservation on computational domain
subroutine sum_En(G, GV, US, CS, En, label)
  type(ocean_grid_type),  intent(in) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),intent(in) :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),  intent(in) :: US !< A dimensional unit scaling type
  type(int_tide_CS),      intent(inout) :: CS !< Internal tide control structure
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle), &
                          intent(in) :: En !< The energy density of the internal tides [H Z2 T-2 ~> m3 s-2 or J m-2].
  character(len=*),       intent(in) :: label !< A label to use in error messages
  ! Local variables
  real :: En_sum   ! The total energy in MKS units for potential output [m5 s-2 or J]
  integer :: a
  ! real :: En_sum_diff  ! Change in energy from the expected value [m5 s-2 or J]
  ! real :: En_sum_pdiff ! Percentage change in energy from the expected value [nondim]
  ! character(len=160) :: mesg  ! The text of an error message
  ! real :: days          ! The time in days for use in output messages [days]

  En_sum = 0.0
  do a=1,CS%nAngle
    En_sum = En_sum + global_area_integral(En(:,:,a), G, unscale=GV%H_to_mks*(US%Z_to_m**2)*(US%s_to_T)**2)
  enddo
  CS%En_sum = En_sum
  !En_sum_diff = En_sum - CS%En_sum
  !if (CS%En_sum /= 0.0) then
  !  En_sum_pdiff= (En_sum_diff/CS%En_sum)*100.0
  !else
  !  En_sum_pdiff= 0.0
  !endif
  !! Print to screen
  !if (is_root_pe()) then
  !  days = time_type_to_real(CS%Time) / 86400.0
  !  write(mesg,*) trim(label)//': days =', days, ', En_sum=', En_sum, &
  !                ', En_sum_diff=', En_sum_diff, ', Percent change=', En_sum_pdiff, '%'
  !  call MOM_mesg(mesg)
  !if (is_root_pe() .and. (abs(En_sum_pdiff) > 1.0)) &
  !  call MOM_error(FATAL, "Run stopped due to excessive internal tide energy change.")
  !endif

end subroutine sum_En

!> Calculates the energy lost from the propagating internal tide due to
!! scattering over small-scale roughness along the lines of Jayne & St. Laurent (2001).
subroutine itidal_lowmode_loss(G, GV, US, CS, Nb, Rho_bot, Ub, En, TKE_loss_fixed, TKE_loss, dt, halo_size)
  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  type(int_tide_CS),         intent(in)    :: CS !< Internal tide control structure
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                             intent(in)    :: Nb !< Near-bottom stratification [T-1 ~> s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                             intent(in)    :: Rho_bot !< Near-bottom density [R ~> kg m-3].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%nFreq,CS%nMode), &
                             intent(inout) :: Ub !< RMS (over one period) near-bottom horizontal
                                                 !! mode velocity [L T-1 ~> m s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                             intent(in) :: TKE_loss_fixed !< Fixed part of energy loss [R Z4 H-1 L-2 ~> kg m-2 or m]
                                                 !! (rho*kappa*h^2) or (kappa*h^2).
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle,CS%nFreq,CS%nMode), &
                             intent(inout) :: En !< Energy density of the internal waves [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle,CS%nFreq,CS%nMode), &
                             intent(out)   :: TKE_loss    !< Energy loss rate [H Z2 T-3 ~> m3 s-3 or W m-2]
                                                 !! (q*rho*kappa*h^2*N*U^2).
  real,                      intent(in)    :: dt !< Time increment [T ~> s].
  integer, optional,         intent(in)    :: halo_size !< The halo size over which to do the calculations
  ! Local variables
  integer :: j, i, m, fr, a, is, ie, js, je, halo
  real    :: En_tot          ! energy for a given mode, frequency
                             ! and point summed over angles [H Z2 T-2 ~> m3 s-2 or J m-2]
  real    :: TKE_loss_tot    ! dissipation for a given mode, frequency
                             ! and point summed over angles [H Z2 T-3 ~> m3 s-3 or W m-2]
  real    :: frac_per_sector ! fraction of energy in each wedge [nondim]
  real    :: q_itides        ! fraction of energy actually lost to mixing (remainder, 1-q, is
                             ! assumed to stay in propagating mode for now - BDM) [nondim]
  real    :: loss_rate       ! approximate loss rate for implicit calc [T-1 ~> s-1]
  real    :: En_negl         ! negligibly small number to prevent division by zero [H Z2 T-2 ~> m3 s-2 or J m-2]
  real    :: En_a, En_b      ! energy before and after timestep [H Z2 T-2 ~> m3 s-2 or J m-2]
  real    :: I_dt            ! The inverse of the timestep [T-1 ~> s-1]
  real    :: J_m2_to_HZ2_T2  ! unit conversion factor for Energy from mks to internal
                             ! units [H Z2 s2 T-2 kg-1 ~> m3 kg-1 or 1]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  J_m2_to_HZ2_T2 = GV%m_to_H*(US%m_to_Z**2)*(US%T_to_s**2)

  I_dt = 1.0 / dt
  q_itides = CS%q_itides
  En_negl = 1e-30*J_m2_to_HZ2_T2

  if (present(halo_size)) then
    halo = halo_size
    is = G%isc - halo ; ie = G%iec + halo ; js = G%jsc - halo ; je = G%jec + halo
  endif

  if (CS%debug) then
    call hchksum(TKE_loss_fixed, "TKE loss fixed", G%HI, haloshift=0, &
                 unscale=US%RZ_to_kg_m2*(US%Z_to_m**3)*GV%m_to_H*(US%m_to_L**2))
    call hchksum(Nb(:,:), "Nbottom", G%HI, haloshift=0, unscale=US%s_to_T)
    call hchksum(Ub(:,:,1,1), "Ubottom", G%HI, haloshift=0, unscale=US%L_to_m*US%s_to_T)
  endif

  do j=js,je ; do i=is,ie ; do m=1,CS%nMode ; do fr=1,CS%nFreq

    ! Sum energy across angles
    En_tot = 0.0
    do a=1,CS%nAngle
      En_tot = En_tot + En(i,j,a,fr,m)
    enddo

    ! Calculate TKE loss rate; units of [H Z2 T-3 ~> m3 s-3 or W m-2] here.
    if (GV%Boussinesq .or. GV%semi_Boussinesq) then
      TKE_loss_tot = q_itides * GV%RZ_to_H*GV%Z_to_H*TKE_loss_fixed(i,j)*Nb(i,j)*Ub(i,j,fr,m)**2
    else
      TKE_loss_tot = q_itides * (GV%RZ_to_H*GV%RZ_to_H*Rho_bot(i,j))*TKE_loss_fixed(i,j)*Nb(i,j)*Ub(i,j,fr,m)**2
    endif

    ! Update energy remaining (this is a pseudo implicit calc)
    ! (E(t+1)-E(t))/dt = -TKE_loss(E(t+1)/E(t)), which goes to zero as E(t+1) goes to zero
    if (En_tot > 0.0) then
      do a=1,CS%nAngle
        frac_per_sector = En(i,j,a,fr,m)/En_tot
        TKE_loss(i,j,a,fr,m) = frac_per_sector*TKE_loss_tot           ! [H Z2 T-3  ~> m3 s-3 or W m-2]
        loss_rate = TKE_loss(i,j,a,fr,m) / (En(i,j,a,fr,m) + En_negl) ! [T-1 ~> s-1]
        En_b = En(i,j,a,fr,m)
        En_a = En(i,j,a,fr,m) / (1.0 + (dt*loss_rate))
        TKE_loss(i,j,a,fr,m) = (En_b - En_a) * I_dt ! overwrite with exact value
        En(i,j,a,fr,m) = En_a
      enddo
    else
      ! no loss if no energy
      do a=1,CS%nAngle
        TKE_loss(i,j,a,fr,m) = 0.0
      enddo
    endif

  enddo ; enddo ; enddo ; enddo

end subroutine itidal_lowmode_loss

!> This subroutine extracts the energy lost from the propagating internal which has
!> been summed across all angles, frequencies, and modes for a given mechanism and location.
!!
!> It can be called from another module to get values from this module's (private) CS.
subroutine get_lowmode_loss(i,j,G,CS,mechanism,TKE_loss_sum)
  integer,               intent(in)  :: i   !< The i-index of the value to be reported.
  integer,               intent(in)  :: j   !< The j-index of the value to be reported.
  type(ocean_grid_type), intent(in)  :: G   !< The ocean's grid structure
  type(int_tide_CS),     intent(in)  :: CS  !< Internal tide control structure
  character(len=*),      intent(in)  :: mechanism    !< The named mechanism of loss to return
  real,                  intent(out) :: TKE_loss_sum !< Total energy loss rate due to specified
                                                     !! mechanism [H Z2 T-3 ~> m3 s-3 or W m-2].

  if (mechanism == 'LeakDrag')  TKE_loss_sum = CS%tot_leak_loss(i,j)
  if (mechanism == 'QuadDrag')  TKE_loss_sum = CS%tot_quad_loss(i,j)
  if (mechanism == 'WaveDrag')  TKE_loss_sum = CS%tot_itidal_loss(i,j)
  if (mechanism == 'Froude')    TKE_loss_sum = CS%tot_Froude_loss(i,j)
  if (mechanism == 'SlopeDrag') TKE_loss_sum = CS%tot_residual_loss(i,j)

end subroutine get_lowmode_loss


!> Returns the values of diffusivity corresponding to various mechanisms
subroutine get_lowmode_diffusivity(G, GV, h, tv, US, h_bot, k_bot, j, N2_lay, N2_int, TKE_to_Kd, Kd_max, CS, &
                                   Kd_leak, Kd_quad, Kd_itidal, Kd_Froude, Kd_slope, &
                                   Kd_lay, Kd_int, profile_leak, profile_quad, profile_itidal, &
                                   profile_Froude, profile_slope)

  type(ocean_grid_type),               intent(in)  :: G       !< The ocean's grid structure
  type(verticalGrid_type),             intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                       intent(in)  :: h       !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),               intent(in)  :: tv      !< Structure containing pointers to any available
  type(unit_scale_type),               intent(in)  :: US      !< A dimensional unit scaling type
  real, dimension(SZI_(G)),            intent(in)  :: h_bot   !< Bottom boundary layer thickness [H ~> m or kg m-2]
  integer, dimension(SZI_(G)),         intent(in)  :: k_bot   !< Bottom boundary layer top layer index
  integer,                             intent(in)  :: j       !< The j-index to work on
  real, dimension(SZI_(G),SZK_(GV)),   intent(in)  :: N2_lay  !< The squared buoyancy frequency of the
                                                              !! layers [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)+1), intent(in)  :: N2_int  !< The squared buoyancy frequency of the
                                                              !! interfaces [T-2 ~> s-2].
  real, dimension(SZI_(G),SZK_(GV)),   intent(in)  :: TKE_to_Kd !< The conversion rate between the TKE
                                                              !! dissipated within a layer and the
                                                              !! diapycnal diffusivity within that layer,
                                                              !! usually (~Rho_0 / (G_Earth * dRho_lay))
                                                              !! [H Z T-1 / H Z2 T-3 = T2 Z-1 ~> s2 m-1]
  real,                                 intent(in) :: Kd_max  !< The maximum increment for diapycnal
                                                              !! diffusivity due to TKE-based processes
                                                              !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
                                                              !! Set this to a negative value to have no limit.
                                                              !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  type(int_tide_cs),                    intent(in)    :: CS   !< The control structure for this module

  real, dimension(SZI_(G),SZK_(GV)+1),  intent(out) :: Kd_leak        !< Diffusivity due to background drag
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),  intent(out) :: Kd_quad        !< Diffusivity due to bottom drag
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),  intent(out) :: Kd_itidal      !< Diffusivity due to wave drag
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),  intent(out) :: Kd_Froude      !< Diffusivity due to high Froude breaking
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),  intent(out) :: Kd_slope       !< Diffusivity due to critical slopes
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)),    intent(inout) :: Kd_lay       !< The diapycnal diffusivity in layers
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G),SZK_(GV)+1),  intent(inout) :: Kd_int       !< The diapycnal diffusivity at interfaces
                                                                      !! [H Z T-1 ~> m2 s-1 or kg m-1 s-1].
  real, dimension(SZI_(G), SZK_(GV)),   intent(out) :: profile_leak   !< Normalized profile for background drag
                                                                      !! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZI_(G), SZK_(GV)),   intent(out) :: profile_quad   !< Normalized profile for  bottom drag
                                                                      !! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZI_(G), SZK_(GV)),   intent(out) :: profile_itidal !< Normalized profile for wave drag
                                                                      !! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZI_(G), SZK_(GV)),   intent(out) :: profile_Froude !< Normalized profile for Froude drag
                                                                      !! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZI_(G), SZK_(GV)),   intent(out) :: profile_slope  !< Normalized profile for critical slopes
                                                                      !! [H-1 ~> m-1 or m2 kg-1]

  ! local variables
  real :: TKE_loss          ! temp variable to pass value of internal tides TKE loss [R Z-3 T-3 ~> W m-2]
  real :: renorm_N          ! renormalization for N profile [H T-1 ~> m s-1 or kg m-2 s-1]
  real :: renorm_N2         ! renormalization for N2 profile [H T-2 ~> m s-2 or kg m-2 s-2]
  real :: tmp_StLau         ! tmp var for renormalization for StLaurent profile [nondim]
  real :: tmp_StLau_slope   ! tmp var for renormalization for StLaurent profile [nondim]
  real :: renorm_StLau      ! renormalization for StLaurent profile [nondim]
  real :: renorm_StLau_slope! renormalization for StLaurent profile [nondim]
  real :: htot              ! total depth of water column [H ~> m or kg m-2]
  real :: htmp              ! local value of thickness in layers [H ~> m or kg m-2]
  real :: h_d               ! expomential decay length scale [H ~> m or kg m-2]
  real :: h_s               ! expomential decay length scale on the slope [H ~> m or kg m-2]
  real :: I_h_d             ! inverse of expomential decay length scale [H-1 ~> m-1 or m2 kg-1]
  real :: I_h_s             ! inverse of expomential decay length scale on the slope [H-1 ~> m-1 or m2 kg-1]
  real :: TKE_to_Kd_lim     ! limited version of TKE_to_Kd [T2 Z-1 ~> s2 m-1]

  ! vertical profiles have units Z-1 for conversion to Kd to be dim correct (see eq 2 of St Laurent GRL 2002)
  real, dimension(SZK_(GV)) :: profile_N               ! vertical profile varying with N [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZK_(GV)) :: profile_N2              ! vertical profile varying with N2 [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZK_(GV)) :: profile_StLaurent       ! vertical profile according to St Laurent 2002
                                                       ! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZK_(GV)) :: profile_StLaurent_slope ! vertical profile according to St Laurent 2002
                                                       ! [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZK_(GV)) :: profile_BBL             ! vertical profile Heavyside BBL  [H-1 ~> m-1 or m2 kg-1]
  real, dimension(SZK_(GV)) :: Kd_leak_lay   ! Diffusivity due to background drag [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZK_(GV)) :: Kd_quad_lay   ! Diffusivity due to bottom drag [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZK_(GV)) :: Kd_itidal_lay ! Diffusivity due to wave drag [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZK_(GV)) :: Kd_Froude_lay ! Diffusivity due to high Froude breaking [H Z T-1 ~> m2 s-1 or kg m-1 s-1]
  real, dimension(SZK_(GV)) :: Kd_slope_lay  ! Diffusivity due to critical slopes [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  real :: hmin                  ! A minimum allowable thickness [H ~> m or kg m-2]
  real :: h_rmn                 ! Remaining thickness in k-loop [H ~> m or kg m-2]
  real :: frac                  ! A fraction of thicknesses [nondim]
  real :: verif_N,   &          ! profile verification [nondim]
          verif_N2,  &          ! profile verification [nondim]
          verif_bbl, &          ! profile verification [nondim]
          verif_stl1,&          ! profile verification [nondim]
          verif_stl2,&          ! profile verification [nondim]
          threshold_renorm_N2,& ! Maximum allowable error on N2 profile [H T-2 ~> m s-2 or kg m-2 s-2]
          threshold_renorm_N, & ! Maximum allowable error on N profile [H T-1 ~> m s-1 or kg m-2 s-1]
          threshold_verif       ! Maximum allowable error on verification [nondim]

  logical :: non_Bous ! fully Non-Boussinesq
  integer :: i, k, is, ie, nz

  is=G%isc ; ie=G%iec ; nz=GV%ke

  non_Bous = .not.(GV%Boussinesq .or. GV%semi_Boussinesq)

  h_d = CS%Int_tide_decay_scale
  h_s = CS%Int_tide_decay_scale_slope
  I_h_d = 1 / h_d
  I_h_s = 1 / h_s

  hmin = 1.0e-6*GV%m_to_H
  threshold_renorm_N2 = 1.0e-13 * GV%m_to_H * US%T_to_s**2
  threshold_renorm_N  = 1.0e-13 * GV%m_to_H * US%T_to_s
  threshold_verif = 1.0e-13

  ! init output arrays
  profile_leak(:,:) = 0.0
  profile_quad(:,:) = 0.0
  profile_slope(:,:) = 0.0
  profile_itidal(:,:) = 0.0
  profile_Froude(:,:) = 0.0

  Kd_leak_lay(:) = 0.0
  Kd_quad_lay(:) = 0.0
  Kd_itidal_lay(:) = 0.0
  Kd_Froude_lay(:) = 0.0
  Kd_slope_lay(:) = 0.0

  Kd_leak(:,:) = 0.0
  Kd_quad(:,:) = 0.0
  Kd_itidal(:,:) = 0.0
  Kd_Froude(:,:) = 0.0
  Kd_slope(:,:) = 0.0

  do i=is,ie

    ! create vertical profiles for diffusivites in layers
    renorm_N = 0.0
    renorm_N2 = 0.0
    renorm_StLau = 0.0
    renorm_StLau_slope = 0.0
    tmp_StLau = 0.0
    tmp_StLau_slope = 0.0
    htot = 0.0
    htmp = 0.0

    do k=1,nz
      ! N-profile
      if (N2_lay(i,k) < 0.) call MOM_error(WARNING, "negative buoyancy freq")
      renorm_N = renorm_N + (sqrt(max(N2_lay(i,k), 0.)) * h(i,j,k))
      ! N2-profile
      renorm_N2 = renorm_N2 + (max(N2_lay(i,k), 0.) * h(i,j,k))
      ! total depth
      htot = htot + h(i,j,k)
    enddo

    profile_N2(:) = 0.0
    profile_N(:) = 0.0
    profile_BBL(:) = 0.0
    profile_StLaurent(:) = 0.0
    profile_StLaurent_slope(:) = 0.0

    ! BBL-profile
    h_rmn = h_bot(i)
    do k=nz,1,-1
      if (G%mask2dT(i,j) > 0.0) then
        profile_BBL(k) = 0.0
        if (h(i,j,k) <= h_rmn) then
          profile_BBL(k) = 1.0 / h_bot(i)
          h_rmn = h_rmn - h(i,j,k)
        else
          if (h_rmn > 0.0) then
            frac = h_rmn / h(i,j,k)
            profile_BBL(k) = frac / h_bot(i)
            h_rmn = h_rmn - frac*h(i,j,k)
          endif
        endif
      endif
    enddo

    do k=1,nz
      if (G%mask2dT(i,j) > 0.0) then
        ! N - profile
        if (renorm_N > threshold_renorm_N) then
           profile_N(k) = sqrt(max(N2_lay(i,k), 0.)) / renorm_N
        else
           profile_N(k) = 1 / htot
        endif

        ! N2 - profile
        if (renorm_N2 > threshold_renorm_N2) then
           profile_N2(k) = max(N2_lay(i,k), 0.) / renorm_N2
        else
           profile_N2(k) = 1 / htot
        endif

        ! slope intensified (St Laurent GRL 2002) - profile
        ! in paper, z is defined positive upwards, range 0 to -H
        ! here depth positive downwards
        ! profiles are almost normalized but differ from a few percent
        ! so we add a second renormalization factor

        ! add first half of layer: get to the layer center
        htmp = htmp + 0.5*h(i,j,k)

        profile_StLaurent(k) = exp(-I_h_d*(htot-htmp)) / &
                              (h_d*(1 - exp(-I_h_d*htot)))

        profile_StLaurent_slope(k) = exp(-I_h_s*(htot-htmp)) / &
                                    (h_s*(1 - exp(-I_h_s*htot)))

        tmp_StLau = tmp_StLau + (profile_StLaurent(k) * h(i,j,k))
        tmp_StLau_slope = tmp_StLau_slope + (profile_StLaurent_slope(k) * h(i,j,k))

        ! add second half of layer: get to the next interface
        htmp = htmp + 0.5*h(i,j,k)
      endif
    enddo

    if (G%mask2dT(i,j) > 0.0) then
      ! allow for difference less than verification threshold
      renorm_StLau = 1.0
      renorm_StLau_slope = 1.0
      if (abs(tmp_StLau -1.0) > threshold_verif) renorm_StLau = 1.0 / tmp_StLau
      if (abs(tmp_StLau_slope -1.0) > threshold_verif) renorm_StLau_slope = 1.0 / tmp_StLau_slope

      do k=1,nz
        profile_StLaurent(k) = profile_StLaurent(k) * renorm_StLau
        profile_StLaurent_slope(k) = profile_StLaurent_slope(k) * renorm_StLau_slope
      enddo
    endif

    ! verif integrals
    if (CS%debug) then
      if (G%mask2dT(i,j) > 0.0) then
         verif_N = 0.0
         verif_N2 = 0.0
         verif_bbl = 0.0
         verif_stl1 = 0.0
         verif_stl2 = 0.0
         do k=1,nz
           verif_N = verif_N + (profile_N(k) * h(i,j,k))
           verif_N2 = verif_N2 + (profile_N2(k) * h(i,j,k))
           verif_bbl = verif_bbl + (profile_BBL(k) * h(i,j,k))
           verif_stl1 = verif_stl1 + (profile_StLaurent(k) * h(i,j,k))
           verif_stl2 = verif_stl2 + (profile_StLaurent_slope(k) * h(i,j,k))
         enddo

         if (abs(verif_N -1.0) > threshold_verif) then
           write(stdout,'(I5,I5,F18.10)') i, j, verif_N
           call MOM_error(FATAL, "mismatch integral for N profile")
         endif
         if (abs(verif_N2 -1.0) > threshold_verif) then
           write(stdout,'(I5,I5,F18.10)') i, j, verif_N2
           call MOM_error(FATAL, "mismatch integral for N2 profile")
         endif
         if (abs(verif_bbl -1.0) > threshold_verif) then
           write(stdout,'(I5,I5,F18.10)') i, j, verif_bbl
           call MOM_error(FATAL, "mismatch integral for bbl profile")
         endif
         if (abs(verif_stl1 -1.0) > threshold_verif) then
           write(stdout,'(I5,I5,F18.10)') i, j, verif_stl1
           call MOM_error(FATAL, "mismatch integral for stl1 profile")
         endif
         if (abs(verif_stl2 -1.0) > threshold_verif) then
           write(stdout,'(I5,I5,F18.10)') i, j, verif_stl2
           call MOM_error(FATAL, "mismatch integral for stl2 profile")
         endif

      endif
    endif

    ! note on units: TKE_to_Kd = 1 / ((g/rho0) * drho) Z-1 T2
    ! mult by dz gives -1/N2 in T2

    ! get TKE loss value and compute diffusivites in layers
    if (CS%apply_background_drag) then
      call get_lowmode_loss(i, j, G, CS, "LeakDrag", TKE_loss)
      ! insert logic to switch between profiles here
      ! if trim(CS%leak_profile) == "N2" then
      profile_leak(i,:) = profile_N2(:)
      ! elseif trim(CS%leak_profile) == "N" then
      ! profile_leak(:) = profile_N(:)
      ! something else
      ! endif
      Kd_leak_lay(:) = 0.
      do k=1,nz
        ! layer diffusivity for processus
        if (h(i,j,k) >= CS%min_thick_layer_Kd) then
          TKE_to_Kd_lim = min(TKE_to_Kd(i,k), CS%max_TKE_to_Kd)
          Kd_leak_lay(k) = CS%mixing_effic * TKE_loss * TKE_to_Kd_lim * profile_leak(i,k) * h(i,j,k)
        else
          Kd_leak_lay(k) = 0.
        endif
        ! add to total Kd in layer
        if (CS%update_Kd) Kd_lay(i,k) = Kd_lay(i,k) + min(Kd_leak_lay(k), Kd_max)
      enddo
    endif

    if (CS%apply_Froude_drag) then
      call get_lowmode_loss(i, j, G, CS, "Froude", TKE_loss)
      ! insert logic to switch between profiles here
      ! if trim(CS%Froude_profile) == "N" then
      profile_Froude(i,:) = profile_N(:)
      ! elseif trim(CS%Froude_profile) == "N2" then
      ! profile_Froude(:) = profile_N2(:)
      ! something else
      ! endif
      do k=1,nz
        ! layer diffusivity for processus
        if (h(i,j,k) >= CS%min_thick_layer_Kd) then
          TKE_to_Kd_lim = min(TKE_to_Kd(i,k), CS%max_TKE_to_Kd)
          Kd_Froude_lay(k) = CS%mixing_effic * TKE_loss * TKE_to_Kd_lim * profile_Froude(i,k) * h(i,j,k)
        else
          Kd_Froude_lay(k) = 0.
        endif
        ! add to total Kd in layer
        if (CS%update_Kd) Kd_lay(i,k) = Kd_lay(i,k) + min(Kd_Froude_lay(k), Kd_max)
      enddo
    endif

    if (CS%apply_wave_drag) then
      call get_lowmode_loss(i, j, G, CS, "WaveDrag", TKE_loss)
      ! insert logic to switch between profiles here
      ! if trim(CS%wave_profile) == "StLaurent" then
      profile_itidal(i,:) = profile_StLaurent(:)
      ! elseif trim(CS%Froude_profile) == "N2" then
      ! profile_itidal(:) = profile_N2(:)
      ! something else
      ! endif
      do k=1,nz
        ! layer diffusivity for processus
        if (h(i,j,k) >= CS%min_thick_layer_Kd) then
          TKE_to_Kd_lim = min(TKE_to_Kd(i,k), CS%max_TKE_to_Kd)
          Kd_itidal_lay(k) = CS%mixing_effic * TKE_loss * TKE_to_Kd_lim * profile_itidal(i,k) * h(i,j,k)
        else
          Kd_itidal_lay(k) = 0.
        endif
        ! add to total Kd in layer
        if (CS%update_Kd) Kd_lay(i,k) = Kd_lay(i,k) + min(Kd_itidal_lay(k), Kd_max)
      enddo
    endif

    if (CS%apply_residual_drag) then
      call get_lowmode_loss(i, j, G, CS, "SlopeDrag", TKE_loss)
      ! insert logic to switch between profiles here
      ! if trim(CS%wave_profile) == "StLaurent" then
      profile_slope(i,:) = profile_StLaurent_slope(:)
      ! elseif trim(CS%Froude_profile) == "N2" then
      ! profile_itidal(:) = profile_N2(:)
      ! something else
      ! endif
      do k=1,nz
        ! layer diffusivity for processus
        if (h(i,j,k) >= CS%min_thick_layer_Kd) then
          TKE_to_Kd_lim = min(TKE_to_Kd(i,k), CS%max_TKE_to_Kd)
          Kd_slope_lay(k) = CS%mixing_effic * TKE_loss * TKE_to_Kd_lim * profile_slope(i,k) * h(i,j,k)
        else
          Kd_slope_lay(k) = 0.
        endif
        ! add to total Kd in layer
        if (CS%update_Kd) Kd_lay(i,k) = Kd_lay(i,k) + min(Kd_slope_lay(k), Kd_max)
      enddo
    endif

    if (CS%apply_bottom_drag) then
      call get_lowmode_loss(i, j, G, CS, "QuadDrag", TKE_loss)
      ! insert logic to switch between profiles here
      ! if trim(CS%bottom_profile) == "BBL" then
      profile_quad(i,:) = profile_BBL(:)
      ! elseif trim(CS%bottom_profile) == "N2" then
      ! profile_quad(:) = profile_N2(:)
      ! something else
      ! endif
      do k=1,nz
        ! layer diffusivity for processus
        if (h(i,j,k) >= CS%min_thick_layer_Kd) then
          TKE_to_Kd_lim = min(TKE_to_Kd(i,k), CS%max_TKE_to_Kd)
          Kd_quad_lay(k) = CS%mixing_effic * TKE_loss * TKE_to_Kd_lim * profile_quad(i,k) * h(i,j,k)
        else
          Kd_quad_lay(k) = 0.
        endif
        ! add to total Kd in layer
        if (CS%update_Kd) Kd_lay(i,k) = Kd_lay(i,k) + min(Kd_quad_lay(k), Kd_max)
      enddo
    endif

    ! interpolate Kd_[] to interfaces and add to Kd_int
    if (CS%apply_background_drag) then
      do k=1,nz+1
        if (k>1)    Kd_leak(i,K) = 0.5*Kd_leak_lay(k-1)
        if (k<nz+1) Kd_leak(i,K) = Kd_leak(i,K) + 0.5*Kd_leak_lay(k)
        ! add to Kd_int
        if (CS%update_Kd) Kd_int(i,K) = min(Kd_int(i,K) + Kd_leak(i,K), Kd_max)
      enddo
    endif

    if (CS%apply_wave_drag) then
      do k=1,nz+1
        if (k>1)    Kd_itidal(i,K) = 0.5*Kd_itidal_lay(k-1)
        if (k<nz+1) Kd_itidal(i,K) = Kd_itidal(i,K) + 0.5*Kd_itidal_lay(k)
        ! add to Kd_int
        if (CS%update_Kd) Kd_int(i,K) = min(Kd_int(i,K) + Kd_itidal(i,K), Kd_max)
      enddo
    endif

    if (CS%apply_Froude_drag) then
      do k=1,nz+1
        if (k>1)    Kd_Froude(i,K) = 0.5*Kd_Froude_lay(k-1)
        if (k<nz+1) Kd_Froude(i,K) = Kd_Froude(i,K) + 0.5*Kd_Froude_lay(k)
        ! add to Kd_int
        if (CS%update_Kd) Kd_int(i,K) = min(Kd_int(i,K) + Kd_Froude(i,K), Kd_max)
      enddo
    endif

    if (CS%apply_residual_drag) then
      do k=1,nz+1
        if (k>1)    Kd_slope(i,K) = 0.5*Kd_slope_lay(k-1)
        if (k<nz+1) Kd_slope(i,K) = Kd_slope(i,K) + 0.5*Kd_slope_lay(k)
        ! add to Kd_int
        if (CS%update_Kd) Kd_int(i,K) = min(Kd_int(i,K) + Kd_slope(i,K), Kd_max)
      enddo
    endif

    if (CS%apply_bottom_drag) then
      do k=1,nz+1
        if (k>1)    Kd_quad(i,K) = 0.5*Kd_quad_lay(k-1)
        if (k<nz+1) Kd_quad(i,K) = Kd_quad(i,K) + 0.5*Kd_quad_lay(k)
        ! add to Kd_int
        if (CS%update_Kd) Kd_int(i,K) = min(Kd_int(i,K) + Kd_quad(i,K), Kd_max)
      enddo
    endif
  enddo ! i-loop

end subroutine get_lowmode_diffusivity

!> Implements refraction on the internal waves at a single frequency.
subroutine refract(En, cn, freq, dt, G, US, NAngle, use_PPMang)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure.
  integer,               intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                         intent(inout) :: En   !< The internal gravity wave energy density as a
                                               !! function of space and angular resolution,
                                               !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed),        &
                         intent(in)    :: cn   !< Baroclinic mode speed [L T-1 ~> m s-1].
  real,                  intent(in)    :: freq !< Wave frequency [T-1 ~> s-1].
  real,                  intent(in)    :: dt   !< Time step [T ~> s].
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  logical,               intent(in)    :: use_PPMang !< If true, use PPM for advection rather
                                               !! than upwind.
  ! Local variables
  integer, parameter :: stencil = 2
  real, dimension(SZI_(G),1-stencil:NAngle+stencil) :: &
    En2d                  ! The internal gravity wave energy density in zonal slices [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(1-stencil:NAngle+stencil) :: &
    cos_angle, sin_angle  ! The cosine and sine of each angle [nondim]
  real, dimension(SZI_(G)) :: &
    Dk_Dt_Kmag, Dl_Dt_Kmag ! Rates of angular refraction [T-1 ~> s-1]
  real, dimension(SZI_(G),0:nAngle) :: &
    Flux_E                ! The flux of energy between successive angular wedges
                          ! within a timestep [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G),1-stencil:NAngle+stencil) :: &
    CFL_ang               ! The CFL number of angular refraction [nondim]
  real, dimension(G%IsdB:G%IedB,G%jsd:G%jed) :: cn_u !< Internal wave group velocity at U-point [L T-1 ~> m s-1]
  real, dimension(G%isd:G%ied,G%JsdB:G%JedB) :: cn_v !< Internal wave group velocity at V-point [L T-1 ~> m s-1]
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: cnmask !< Local mask for group velocity [nondim]
  real :: f2              ! The squared Coriolis parameter [T-2 ~> s-2].
  real :: favg            ! The average Coriolis parameter at a point [T-1 ~> s-1].
  real :: df_dy, df_dx    ! The x- and y- gradients of the Coriolis parameter [T-1 L-1 ~> s-1 m-1].
  real :: dlnCn_dx        ! The x-gradient of the wave speed divided by itself [L-1 ~> m-1].
  real :: dlnCn_dy        ! The y-gradient of the wave speed divided by itself [L-1 ~> m-1].
  real :: Angle_size      ! The size of each wedge of angles [rad]
  real :: dt_Angle_size   ! The time step divided by the angle size [T rad-1 ~> s rad-1]
  real :: angle           ! The central angle of each wedge [rad]
  real :: Ifreq           ! The inverse of the wave frequency [T ~> s]
  real :: Kmag2           ! A squared horizontal wavenumber [L-2 ~> m-2]
  real :: I_Kmag          ! The inverse of the magnitude of the horizontal wavenumber [L ~> m]
  real :: cn_subRO        ! A tiny wave speed to prevent division by zero [L T-1 ~> m s-1]
  integer :: is, ie, js, je, asd, aed, na
  integer :: i, j, a
  real :: wgt1, wgt2      ! Weights in an average, both of which may be 0 [nondim]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; na = size(En,3)
  asd = 1-stencil ; aed = NAngle+stencil

  cnmask(:,:) = merge(0., 1., cn(:,:) == 0.)

  do j=js,je ; do I=is-1,ie
    ! wgt = 0 if local cn == 0, wgt = 0.5 if both contiguous values != 0
    ! and wgt = 1 if neighbour cn == 0
    wgt1 = cnmask(i,j) - (0.5 * cnmask(i,j) * cnmask(i+1,j))
    wgt2 = cnmask(i+1,j) - (0.5 * cnmask(i,j) * cnmask(i+1,j))
    cn_u(I,j) = (wgt1*cn(i,j)) + (wgt2*cn(i+1,j))
  enddo ; enddo

  do J=js-1,je ; do i=is,ie
    wgt1 = cnmask(i,j) - (0.5 * cnmask(i,j) * cnmask(i,j+1))
    wgt2 = cnmask(i,j+1) - (0.5 * cnmask(i,j) * cnmask(i,j+1))
    cn_v(i,J) = (wgt1*cn(i,j)) + (wgt2*cn(i,j+1))
  enddo ; enddo

  Ifreq = 1.0 / freq
  cn_subRO = 1e-30*US%m_s_to_L_T
  Angle_size = (8.0*atan(1.0)) / (real(NAngle))
  dt_Angle_size = dt / Angle_size

  do A=asd,aed
    angle = (real(A) - 0.5) * Angle_size
    cos_angle(A) = cos(angle) ; sin_angle(A) = sin(angle)
  enddo

  !### There should also be refraction due to cn.grad(grid_orientation).
  CFL_ang(:,:,:) = 0.0
  do j=js,je
  ! Copy En into angle space with halos.
    do a=1,na ; do i=is,ie
      En2d(i,a) = En(i,j,a)
    enddo ; enddo
    do a=asd,0 ; do i=is,ie
      En2d(i,a) = En2d(i,a+NAngle)
      En2d(i,NAngle+stencil+a) = En2d(i,stencil+a)
    enddo ; enddo

  ! Do the refraction.
    do i=is,ie
      f2 = 0.25* ((G%Coriolis2Bu(I,J) + G%Coriolis2Bu(I-1,J-1)) + &
                 (G%Coriolis2Bu(I,J-1) + G%Coriolis2Bu(I-1,J)))
      favg = 0.25*((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1)) + &
                   (G%CoriolisBu(I,J-1) + G%CoriolisBu(I-1,J)))
      df_dx = 0.5*G%IdxT(i,j)*((G%CoriolisBu(I,J) - G%CoriolisBu(I-1,J-1)) + &
                               (G%CoriolisBu(I,J-1) - G%CoriolisBu(I-1,J)))
      df_dy = 0.5*G%IdyT(i,j)*((G%CoriolisBu(I,J) - G%CoriolisBu(I-1,J-1)) + &
                               (G%CoriolisBu(I-1,J) - G%CoriolisBu(I,J-1)))

      dlnCn_dx = G%IdxT(i,j) * (cn_u(I,j) - cn_u(I-1,j)) / (0.5 * (cn_u(I,j) + cn_u(I-1,j)) + cn_subRO)
      dlnCn_dy = G%IdyT(i,j) * (cn_v(i,J) - cn_v(i,J-1)) / (0.5 * (cn_v(i,J) + cn_v(i,J-1)) + cn_subRO)

      Kmag2 = (freq**2 - f2) / (cn(i,j)**2 + cn_subRO**2)
      if (Kmag2 > 0.0) then
        I_Kmag = 1.0 / sqrt(Kmag2)
        Dk_Dt_Kmag(i) = -Ifreq * (favg*df_dx + (freq**2 - f2) * dlnCn_dx) * I_Kmag
        Dl_Dt_Kmag(i) = -Ifreq * (favg*df_dy + (freq**2 - f2) * dlnCn_dy) * I_Kmag
      else
        Dk_Dt_Kmag(i) = 0.0
        Dl_Dt_Kmag(i) = 0.0
      endif
    enddo

    ! Determine the energy fluxes in angular orientation space.
    do A=asd,aed ; do i=is,ie
      CFL_ang(i,j,A) = ((cos_angle(A) * Dl_Dt_Kmag(i)) - (sin_angle(A) * Dk_Dt_Kmag(i))) * dt_Angle_size
      if (abs(CFL_ang(i,j,A)) > 1.0) then
        call MOM_error(WARNING, "refract: CFL exceeds 1.", .true.)
        if (CFL_ang(i,j,A) > 1.0) then ; CFL_ang(i,j,A) = 1.0 ; else ; CFL_ang(i,j,A) = -1.0 ; endif
      endif
    enddo ; enddo

    ! Advect in angular space
    if (.not.use_PPMang) then
      ! Use simple upwind
      do  A=0,na ; do i=is,ie
        if (CFL_ang(i,j,A) > 0.0) then
          Flux_E(i,A) = CFL_ang(i,j,A) * En2d(i,A)
        else
          Flux_E(i,A) = CFL_ang(i,j,A) * En2d(i,A+1)
        endif
      enddo ; enddo
    else
      ! Use PPM
      do i=is,ie
        call PPM_angular_advect(En2d(i,:),CFL_ang(i,j,:),Flux_E(i,:),NAngle,dt,stencil)
      enddo
    endif

  ! Update and copy back to En.
    do a=1,na ; do i=is,ie
      !if (En2d(i,a)+(Flux_E(i,A-1)-Flux_E(i,A)) < 0.0) then ! for debugging
      !  call MOM_error(FATAL, "refract: OutFlux>Available")
      !endif
      En(i,j,a) = En2d(i,a) + (Flux_E(i,A-1) - Flux_E(i,A))
    enddo ; enddo
  enddo ! j-loop
end subroutine refract

!> This subroutine calculates the 1-d flux for advection in angular space using a monotonic
!! piecewise parabolic scheme. This needs to be called from within i and j spatial loops.
subroutine PPM_angular_advect(En2d, CFL_ang, Flux_En, NAngle, dt, halo_ang)
  integer,                   intent(in)    :: NAngle  !< The number of wave orientations in the
                                                      !! discretized wave energy spectrum [nondim]
  real,                      intent(in)    :: dt      !< Time increment [T ~> s].
  integer,                   intent(in)    :: halo_ang !< The halo size in angular space
  real, dimension(1-halo_ang:NAngle+halo_ang),   &
                             intent(in)    :: En2d    !< The internal gravity wave energy density as a
                                                      !! function of angular resolution [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(1-halo_ang:NAngle+halo_ang),   &
                             intent(in)    :: CFL_ang !< The CFL number of the energy advection across angles [nondim]
  real, dimension(0:NAngle), intent(out)   :: Flux_En !< The time integrated internal wave energy flux
                                                      !! across angles  [H Z2 T-2 ~> m3 s-2 or J m-2].
  ! Local variables
  real :: flux         ! The internal wave energy flux across angles  [H Z2 T-3 ~> m3 s-3 or W m-2].
  real :: u_ang        ! Angular propagation speed [Rad T-1 ~> Rad s-1]
  real :: Angle_size   ! The size of each orientation wedge in radians [Rad]
  real :: I_Angle_size ! The inverse of the orientation wedges [Rad-1]
  real :: I_dt         ! The inverse of the timestep [T-1 ~> s-1]
  real :: aR, aL       ! Left and right edge estimates of energy density [H Z2 T-2 rad-1 ~> m3 s-2 rad-1 or J m-2 rad-1]
  real :: Ep, Ec, Em   ! Mean angular energy density for three successive wedges in angular
                       ! orientation [H Z2 T-2 rad-1 ~> m3 s-2 rad-1 or J m-2 rad-1]
  real :: dA, curv_3   ! Difference and curvature of energy density [H Z2 T-2 rad-1 ~> m3 s-2 rad-1 or J m-2 rad-1]
  real, parameter :: oneSixth = 1.0/6.0  ! One sixth [nondim]
  integer :: a

  I_dt = 1.0 / dt
  Angle_size = (8.0*atan(1.0)) / (real(NAngle))
  I_Angle_size = 1 / Angle_size
  Flux_En(:) = 0

  do A=0,NAngle
    u_ang = CFL_ang(A)*Angle_size*I_dt
    if (u_ang >= 0.0) then
      ! Implementation of PPM-H3
      ! Convert wedge-integrated energy density into angular energy densities for three successive
      ! wedges around the source wedge for this flux [H Z2 T-2 rad-1 ~> m3 s-2 rad-1 or J m-2 rad-1].
      Ep = En2d(a+1)*I_Angle_size
      Ec = En2d(a)  *I_Angle_size
      Em = En2d(a-1)*I_Angle_size
      ! Calculate and bound edge values of energy density.
      aL = ( 5.*Ec + ( 2.*Em - Ep ) ) * oneSixth ! H3 estimate
      aL = max( min(Ec,Em), aL) ; aL = min( max(Ec,Em), aL) ! Bound
      aR = ( 5.*Ec + ( 2.*Ep - Em ) ) * oneSixth ! H3 estimate
      aR = max( min(Ec,Ep), aR) ; aR = min( max(Ec,Ep), aR) ! Bound
      dA = aR - aL
      if ((Ep-Ec)*(Ec-Em) <= 0.) then
        aL = Ec ; aR = Ec    ! use PCM for local extremum
      elseif ( 3.0*dA*(2.*Ec - (aR + aL)) > (dA*dA) ) then
        aL = 3.*Ec - 2.*aR   ! Flatten the profile to move the extremum to the left edge
      elseif ( 3.0*dA*(2.*Ec - (aR + aL)) < - (dA*dA) ) then
        aR = 3.*Ec - 2.*aL   ! Flatten the profile to move the extremum to the right edge
      endif
      curv_3 = (aR + aL) - 2.0*Ec ! Curvature
      ! Calculate angular flux rate [H Z2 T-3 ~> m3 s-3 or W m-2]
      flux = u_ang*( aR + CFL_ang(A) * ( 0.5*(aL - aR) + curv_3 * (CFL_ang(A) - 1.5) ) )
      ! Calculate amount of energy fluxed between wedges [H Z2 T-2 ~> m3 s-2 or J m-2]
      Flux_En(A) = dt * flux
      !Flux_En(A) = (dt * I_Angle_size) * flux
    else
      ! Implementation of PPM-H3
      ! Convert wedge-integrated energy density into angular energy densities for three successive
      ! wedges around the source wedge for this flux [H Z2 T-2 rad-1 ~> m3 s-2 rad-1 or J m-2 rad-1].
      Ep = En2d(a+2)*I_Angle_size
      Ec = En2d(a+1)*I_Angle_size
      Em = En2d(a)  *I_Angle_size
      ! Calculate and bound edge values of energy density.
      aL = ( 5.*Ec + ( 2.*Em - Ep ) ) * oneSixth ! H3 estimate
      aL = max( min(Ec,Em), aL) ; aL = min( max(Ec,Em), aL) ! Bound
      aR = ( 5.*Ec + ( 2.*Ep - Em ) ) * oneSixth ! H3 estimate
      aR = max( min(Ec,Ep), aR) ; aR = min( max(Ec,Ep), aR) ! Bound
      dA = aR - aL
      if ((Ep-Ec)*(Ec-Em) <= 0.) then
        aL = Ec ; aR = Ec    ! use PCM for local extremum
      elseif ( 3.0*dA*(2.*Ec - (aR + aL)) > (dA*dA) ) then
        aL = 3.*Ec - 2.*aR   ! Flatten the profile to move the extremum to the left edge
      elseif ( 3.0*dA*(2.*Ec - (aR + aL)) < - (dA*dA) ) then
        aR = 3.*Ec - 2.*aL   ! Flatten the profile to move the extremum to the right edge
      endif
      curv_3 = (aR + aL) - 2.0*Ec ! Curvature
      ! Calculate angular flux rate [H Z2 T-3 ~> m3 s-3 or W m-2]
      ! Note that CFL_ang is negative here, so it looks odd compared with equivalent expressions.
      flux = u_ang*( aL - CFL_ang(A) * ( 0.5*(aR - aL) + curv_3 * (-CFL_ang(A) - 1.5) ) )
      ! Calculate amount of energy fluxed between wedges [H Z2 T-2 ~> m3 s-2 or J m-2]
      Flux_En(A) = dt * flux
      !Flux_En(A) = (dt * I_Angle_size) * flux
    endif
  enddo
end subroutine PPM_angular_advect

!> Propagates internal waves at a single frequency.
subroutine propagate(En, cn, freq, dt, G, GV, US, CS, NAngle, residual_loss)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure.
  integer,               intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                         intent(inout) :: En   !< The internal gravity wave energy density as a
                                               !! function of space and angular resolution,
                                               !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed),        &
                         intent(in)    :: cn   !< Baroclinic mode speed [L T-1 ~> m s-1].
  real,                  intent(in)    :: freq !< Wave frequency [T-1 ~> s-1].
  real,                  intent(in)    :: dt   !< Time step [T ~> s].
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(int_tide_CS),     intent(inout)    :: CS   !< Internal tide control structure
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                         intent(inout) :: residual_loss !< internal tide energy loss due
                                                        !! to the residual at slopes [H Z2 T-3 ~> m3 s-3 or W m-2].
  ! Local variables
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB) :: &
    speed  ! The magnitude of the group velocity at the q points for corner adv [L T-1 ~> m s-1].
  integer, parameter :: stencil = 2
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    speed_x  ! The magnitude of the group velocity at the Cu points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    speed_y  ! The magnitude of the group velocity at the Cv points [L T-1 ~> m s-1].
  real, dimension(0:NAngle) :: &
    cos_angle, sin_angle  ! The cosine and sine of each angle [nondim]
  real, dimension(NAngle) :: &
    Cgx_av, &  ! The average projection of the wedge into the x-direction [nondim]
    Cgy_av, &  ! The average projection of the wedge into the y-direction [nondim]
    dCgx, &    ! The difference in x-projections between the edges of each angular band [nondim].
    dCgy       ! The difference in y-projections between the edges of each angular band [nondim].
  real :: f2   ! The squared Coriolis parameter [T-2 ~> s-2].
  real :: Angle_size      ! The size of each wedge of angles [rad]
  real :: I_Angle_size    ! The inverse of the size of each wedge of angles [rad-1]
  real :: angle           ! The central angle of each wedge [rad]
  real :: Ifreq ! The inverse of the frequency [T ~> s]
  real :: freq2 ! The frequency squared [T-2 ~> s-2]
  type(loop_bounds_type) :: LB
  logical :: x_first
  integer :: is, ie, js, je, asd, aed, na
  integer :: ish, ieh, jsh, jeh
  integer :: i, j, a, fr, m

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; na = size(En,3)
  asd = 1-stencil ; aed = NAngle+stencil

  Ifreq = 1.0 / freq
  freq2 = freq**2

  ! Define loop bounds: Need extensions on j-loop so propagate_y
  ! (done after propagate_x) will have updated values in the halo
  ! for correct PPM reconstruction. Use if no teleporting and
  ! no pass_var between propagate_x and propagate_y.
  !jsh = js-3 ; jeh = je+3 ; ish = is ; ieh = ie

  ! Define loop bounds: Need 1-pt extensions on loops because
  ! teleporting eats up a halo point. Use if teleporting.
  ! Also requires pass_var before propagate_y.
  jsh = js-1 ; jeh = je+1 ; ish = is-1 ; ieh = ie+1

  Angle_size = (8.0*atan(1.0)) / real(NAngle)
  I_Angle_size = 1.0 / Angle_size

  x_first = .true. ! x_first = (MOD(G%first_direction,2) == 0)

  if (CS%debug) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'propagate: top of routine')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'propagate: top of routine', CS%En_sum
    enddo ; enddo
  endif

  ! IMPLEMENT PPM ADVECTION IN HORIZONTAL-----------------------
  ! These could be in the control structure, as they do not vary.
  do A=0,na
    ! These are the angles at the cell edges...
    angle = (real(A) - 0.5) * Angle_size
    cos_angle(A) = cos(angle) ; sin_angle(A) = sin(angle)
  enddo

  do a=1,na
    Cgx_av(a) = (sin_angle(A) - sin_angle(A-1)) * I_Angle_size
    Cgy_av(a) = -(cos_angle(A) - cos_angle(A-1)) * I_Angle_size
    dCgx(a) = sqrt(0.5 + 0.5*(sin_angle(A)*cos_angle(A) - &
                              sin_angle(A-1)*cos_angle(A-1)) * I_Angle_size - &
                   Cgx_av(a)**2)
    dCgy(a) = sqrt(0.5 - 0.5*(sin_angle(A)*cos_angle(A) - &
                              sin_angle(A-1)*cos_angle(A-1)) * I_Angle_size - &
                   Cgy_av(a)**2)
  enddo

  speed_x(:,:) = 0.
  do j=jsh,jeh ; do I=ish-1,ieh
    f2 = 0.5 * (G%Coriolis2Bu(I,J) + G%Coriolis2Bu(I,J-1))
    speed_x(I,j) = 0.5*(cn(i,j) + cn(i+1,j)) * G%mask2dCu(I,j) * &
                   sqrt(max(freq2 - f2, 0.0)) * Ifreq
  enddo ; enddo

  speed_y(:,:) = 0.
  do J=jsh-1,jeh ; do i=ish,ieh
    f2 = 0.5 * (G%Coriolis2Bu(I,J) + G%Coriolis2Bu(I-1,J))
    speed_y(i,J) = 0.5*(cn(i,j) + cn(i,j+1)) * G%mask2dCv(i,J) * &
                   sqrt(max(freq2 - f2, 0.0)) * Ifreq
  enddo ; enddo

  call pass_vector(speed_x, speed_y, G%Domain, stagger=CGRID_NE)
  call pass_var(En, G%domain)

  ! Apply propagation in the first direction (reflection included)
  LB%jsh = jsh ; LB%jeh = jeh ; LB%ish = ish ; LB%ieh = ieh
  if (x_first) then
    call propagate_x(En, speed_x, Cgx_av, dCgx, dt, G, US, CS%nAngle, CS, LB, residual_loss)
  else
    call propagate_y(En, speed_y, Cgy_av, dCgy, dt, G, US, CS%nAngle, CS, LB, residual_loss)
  endif

  ! fix underflows
  do a=1,na ; do j=jsh,jeh ; do i=ish,ieh
    if (abs(En(i,j,a)) < CS%En_underflow) En(i,j,a) = 0.0
  enddo ; enddo ; enddo

  if (CS%debug) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'propagate: after propagate_x')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'propagate: after propagate_x', CS%En_sum
    enddo ; enddo
  endif

  ! Update halos
  call pass_var(En, G%domain)
  call pass_var(residual_loss, G%domain)

  if (CS%debug) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'propagate: after halo update')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'propagate: after halo update', CS%En_sum
    enddo ; enddo
  endif
  ! Apply propagation in the second direction (reflection included)
  ! LB%jsh = js ; LB%jeh = je ; LB%ish = is ; LB%ieh = ie ! Use if no teleport
  LB%jsh = jsh ; LB%jeh = jeh ; LB%ish = ish ; LB%ieh = ieh
  if (x_first) then
    call propagate_y(En, speed_y, Cgy_av, dCgy, dt, G, US, CS%nAngle, CS, LB, residual_loss)
  else
    call propagate_x(En, speed_x, Cgx_av, dCgx, dt, G, US, CS%nAngle, CS, LB, residual_loss)
  endif

  ! fix underflows
  do a=1,na ; do j=jsh,jeh ; do i=ish,ieh
    if (abs(En(i,j,a)) < CS%En_underflow) En(i,j,a) = 0.0
  enddo ; enddo ; enddo

  call pass_var(En, G%domain)
  call pass_var(residual_loss, G%domain)

  if (CS%debug) then
    do m=1,CS%nMode ; do fr=1,CS%Nfreq
      call sum_En(G, GV, US, CS, CS%En(:,:,:,fr,m), 'propagate: bottom of routine')
      if (is_root_pe()) write(stdout,'(A,E18.10)') 'propagate: bottom of routine', CS%En_sum
    enddo ; enddo
  endif

end subroutine propagate


!> Propagates the internal wave energy in the logical x-direction.
subroutine propagate_x(En, speed_x, Cgx_av, dCgx, dt, G, US, Nangle, CS, LB, residual_loss)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  integer,                 intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle),   &
                           intent(inout) :: En !< The energy density integrated over an angular
                                               !! band [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(G%IsdB:G%IedB,G%jsd:G%jed),        &
                           intent(in)    :: speed_x !< The magnitude of the group velocity at the
                                               !! Cu points [L T-1 ~> m s-1].
  real, dimension(Nangle), intent(in)    :: Cgx_av !< The average x-projection in each angular band [nondim]
  real, dimension(Nangle), intent(in)    :: dCgx !< The difference in x-projections between the
                                               !! edges of each angular band [nondim].
  real,                    intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(int_tide_CS),       intent(in)    :: CS !< Internal tide control structure
  type(loop_bounds_type),  intent(in)    :: LB !< A structure with the active energy loop bounds.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle),   &
                           intent(inout) :: residual_loss !< internal tide energy loss due
                                                          !! to the residual at slopes [H Z2 T-3 ~> m3 s-3 or W m-2].
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    EnL, EnR    ! Left and right face energy densities [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    flux_x      ! The internal wave energy flux [H Z2 L2 T-3 ~> m5 s-3 or J s-1].
  real, dimension(SZIB_(G)) :: &
    cg_p, &     ! The x-direction group velocity [L T-1 ~> m s-1]
    flux1       ! A 1-d copy of the x-direction internal wave energy flux [H Z2 L2 T-3 ~> m5 s-3 or J s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle) :: &
    Fdt_m, Fdt_p! Left and right energy fluxes [H Z2 L2 T-2 ~> m5 s-2 or J]
  integer :: i, j, ish, ieh, jsh, jeh, a

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh
  do a=1,Nangle
  ! This sets EnL and EnR.
    if (CS%upwind_1st) then
      do j=jsh,jeh ; do i=ish-1,ieh+1
        EnL(i,j) = En(i,j,a) ; EnR(i,j) = En(i,j,a)
      enddo ; enddo
    else
      call PPM_reconstruction_x(En(:,:,a), EnL, EnR, G, LB, simple_2nd=CS%simple_2nd)
    endif

    do j=jsh,jeh
      ! This is done once with single speed (GARDEN SPRINKLER EFFECT POSSIBLE)
      do I=ish-1,ieh
        cg_p(I) = speed_x(I,j) * (Cgx_av(a))
      enddo
      call zonal_flux_En(cg_p, En(:,j,a), EnL(:,j), EnR(:,j), flux1, &
                         dt, G, US, j, ish, ieh, CS%vol_CFL)
      do I=ish-1,ieh ; flux_x(I,j) = flux1(I); enddo
    enddo

    do j=jsh,jeh ; do i=ish,ieh
      Fdt_m(i,j,a) = dt*flux_x(I-1,j) ! left face influx  [H Z2 L2 T-2 ~> m5 s-2 or J]
      Fdt_p(i,j,a) = -dt*flux_x(I,j)  ! right face influx [H Z2 L2 T-2 ~> m5 s-2 or J]

      ! only compute residual loss on partial reflection cells, remove numerical noise elsewhere
      if (CS%refl_pref_logical(i,j)) then
        residual_loss(i,j,a) = residual_loss(i,j,a) + &
                              ((abs(flux_x(I-1,j)) * CS%residual(i,j) * G%IareaT(i,j)) + &
                               (abs(flux_x(I,j)) * CS%residual(i,j) * G%IareaT(i,j)))
      endif
    enddo ; enddo

  enddo ! a-loop

  ! Only reflect newly arrived energy; existing energy in incident wedge is not reflected
  ! and will eventually propagate out of cell. (This code only reflects if En > 0.)
  call reflect(Fdt_m, Nangle, CS, G, LB)
  !call teleport(Fdt_m, Nangle, CS, G, LB)
  call reflect(Fdt_p, Nangle, CS, G, LB)
  !call teleport(Fdt_p, Nangle, CS, G, LB)

  ! Update reflected energy [H Z2 T-2 ~> m3 s-2 or J m-2]
  do a=1,Nangle ; do j=jsh,jeh ; do i=ish,ieh
    En(i,j,a) = En(i,j,a) + (G%IareaT(i,j)*(Fdt_m(i,j,a) + Fdt_p(i,j,a)))
  enddo ; enddo ; enddo

end subroutine propagate_x

!> Propagates the internal wave energy in the logical y-direction.
subroutine propagate_y(En, speed_y, Cgy_av, dCgy, dt, G, US, Nangle, CS, LB, residual_loss)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  integer,                 intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle), &
                           intent(inout) :: En !< The energy density integrated over an angular
                                               !! band [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(G%isd:G%ied,G%JsdB:G%JedB),      &
                           intent(in)    :: speed_y !< The magnitude of the group velocity at the
                                               !! Cv points [L T-1 ~> m s-1].
  real, dimension(Nangle), intent(in)    :: Cgy_av !< The average y-projection in each angular band [nondim]
  real, dimension(Nangle), intent(in)    :: dCgy !< The difference in y-projections between the
                                               !! edges of each angular band [nondim]
  real,                    intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(int_tide_CS),       intent(in)    :: CS !< Internal tide control structure
  type(loop_bounds_type),  intent(in)    :: LB !< A structure with the active energy loop bounds.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle),   &
                           intent(inout) :: residual_loss !< internal tide energy loss due
                                                          !! to the residual at slopes [H Z2 T-3 ~> m3 s-3 or W m-2].
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    EnL, EnR    ! South and north face energy densities [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    flux_y      ! The internal wave energy flux [H Z2 L2 T-3 ~> m5 s-3 or J s-1].
  real, dimension(SZI_(G)) :: &
    cg_p, &     ! The y-direction group velocity [L T-1 ~> m s-1]
    flux1       ! A 1-d copy of the y-direction internal wave energy flux [H Z2 L2 T-3 ~> m5 s-3 or J s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle) :: &
    Fdt_m, Fdt_p! South and north energy fluxes [H Z2 L2 T-2 ~> m5 s-2 or J]
  integer :: i, j, ish, ieh, jsh, jeh, a

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh
  do a=1,Nangle
    ! This sets EnL and EnR.
    if (CS%upwind_1st) then
      do j=jsh-1,jeh+1 ; do i=ish,ieh
        EnL(i,j) = En(i,j,a) ; EnR(i,j) = En(i,j,a)
      enddo ; enddo
    else
      call PPM_reconstruction_y(En(:,:,a), EnL, EnR, G, LB, simple_2nd=CS%simple_2nd)
    endif

    do J=jsh-1,jeh
      !   This is done once with single speed (GARDEN SPRINKLER EFFECT POSSIBLE)
      do i=ish,ieh
        cg_p(i) = speed_y(i,J) * (Cgy_av(a))
      enddo
      call merid_flux_En(cg_p, En(:,:,a), EnL(:,:), EnR(:,:), flux1, &
                         dt, G, US, J, ish, ieh, CS%vol_CFL)
      do i=ish,ieh ; flux_y(i,J) = flux1(i); enddo
    enddo

    do j=jsh,jeh ; do i=ish,ieh
      Fdt_m(i,j,a) = dt*flux_y(i,J-1) ! south face influx [H Z2 L2 T-2 ~> m5 s-2 or J]
      Fdt_p(i,j,a) = -dt*flux_y(i,J)  ! north face influx [H Z2 L2 T-2 ~> m5 s-2 or J]

      ! only compute residual loss on partial reflection cells, remove numerical noise elsewhere
      if (CS%refl_pref_logical(i,j)) then
        residual_loss(i,j,a) = residual_loss(i,j,a) + &
                              ((abs(flux_y(i,J-1)) * CS%residual(i,j) * G%IareaT(i,j)) + &
                               (abs(flux_y(i,J)) * CS%residual(i,j) * G%IareaT(i,j)))
      endif
    enddo ; enddo

  enddo ! a-loop

  ! Only reflect newly arrived energy; existing energy in incident wedge is not reflected
  ! and will eventually propagate out of cell. (This code only reflects if En > 0.)
  call reflect(Fdt_m, Nangle, CS, G, LB)
  !call teleport(Fdt_m, Nangle, CS, G, LB)
  call reflect(Fdt_p, Nangle, CS, G, LB)
  !call teleport(Fdt_p, Nangle, CS, G, LB)

  ! Update reflected energy [H Z2 T-2 ~> m3 s-2 or J m-2]
  do a=1,Nangle ; do j=jsh,jeh ; do i=ish,ieh
    En(i,j,a) = En(i,j,a) + G%IareaT(i,j)*(Fdt_m(i,j,a) + Fdt_p(i,j,a))
  enddo ; enddo ; enddo

end subroutine propagate_y

!> Evaluates the zonal mass or volume fluxes in a layer.
subroutine zonal_flux_En(u, h, hL, hR, uh, dt, G, US, j, ish, ieh, vol_CFL)
  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure.
  real, dimension(SZIB_(G)), intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G)),  intent(in)    :: h  !< Energy density used to calculate the fluxes
                                                 !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G)),  intent(in)    :: hL !< Left- Energy densities in the reconstruction
                                                 !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G)),  intent(in)    :: hR !< Right- Energy densities in the reconstruction
                                                 !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZIB_(G)), intent(out) :: uh !< The zonal energy transport [H Z2 L2 T-3 ~> m5 s-3 or J s-1].
  real,                      intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  integer,                   intent(in)    :: j  !< The j-index to work on.
  integer,                   intent(in)    :: ish !< The start i-index range to work on.
  integer,                   intent(in)    :: ieh !< The end i-index range to work on.
  logical,                   intent(in)    :: vol_CFL !< If true, rescale the ratio of face areas to
                                                 !! the cell areas when estimating the CFL number.
  ! Local variables
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim].
  real :: curv_3 ! A measure of the energy density curvature over a grid length [H Z2 T-2 ~> m3 s-2 or J m-2]
  integer :: i

  do I=ish-1,ieh
    ! Set new values of uh and duhdu.
    if (u(I) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      else ; CFL = u(I) * dt * G%IdxT(i,j) ; endif
      curv_3 = (hL(i) + hR(i)) - 2.0*h(i)
      uh(I) = G%dy_Cu(I,j) * u(I) * &
          (hR(i) + CFL * (0.5*(hL(i) - hR(i)) + curv_3*(CFL - 1.5)))
    elseif (u(I) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I) * dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else ; CFL = -u(I) * dt * G%IdxT(i+1,j) ; endif
      curv_3 = (hL(i+1) + hR(i+1)) - 2.0*h(i+1)
      uh(I) = G%dy_Cu(I,j) * u(I) * &
          (hL(i+1) + CFL * (0.5*(hR(i+1)-hL(i+1)) + curv_3*(CFL - 1.5)))
    else
      uh(I) = 0.0
    endif
  enddo
end subroutine zonal_flux_En

!> Evaluates the meridional mass or volume fluxes in a layer.
subroutine merid_flux_En(v, h, hL, hR, vh, dt, G, US, J, ish, ieh, vol_CFL)
  type(ocean_grid_type),            intent(in)    :: G  !< The ocean's grid structure.
  real, dimension(SZI_(G)),         intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: h  !< Energy density used to calculate the
                                                        !! fluxes [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: hL !< Left- Energy densities in the
                                                        !! reconstruction [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: hR !< Right- Energy densities in the
                                                        !! reconstruction [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G)),         intent(out) :: vh !< The meridional energy transport
                                                      !! [H Z2 L2 T-3 ~> m5 s-3 or J s-1].
  real,                             intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),            intent(in)    :: US !< A dimensional unit scaling type
  integer,                          intent(in)    :: J  !< The j-index to work on.
  integer,                          intent(in)    :: ish !< The start i-index range to work on.
  integer,                          intent(in)    :: ieh !< The end i-index range to work on.
  logical,                          intent(in)    :: vol_CFL !< If true, rescale the ratio of face
                                                        !! areas to the cell areas when estimating
                                                        !! the CFL number.
  ! Local variables
  real :: CFL ! The CFL number based on the local velocity and grid spacing [nondim].
  real :: curv_3 ! A measure of the energy density curvature over a grid length [H Z2 T-2 ~> m3 s-2 or J m-2]
  integer :: i

  do i=ish,ieh
    if (v(i) > 0.0) then
      if (vol_CFL) then ; CFL = (v(i) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
      else ; CFL = v(i) * dt * G%IdyT(i,j) ; endif
      curv_3 = (hL(i,j) + hR(i,j)) - 2.0*h(i,j)
      vh(i) = G%dx_Cv(i,J) * v(i) * ( hR(i,j) + CFL * &
          (0.5*(hL(i,j) - hR(i,j)) + curv_3*(CFL - 1.5)) )
    elseif (v(i) < 0.0) then
      if (vol_CFL) then ; CFL = (-v(i) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
      else ; CFL = -v(i) * dt * G%IdyT(i,j+1) ; endif
      curv_3 = (hL(i,j+1) + hR(i,j+1)) - 2.0*h(i,j+1)
      vh(i) = G%dx_Cv(i,J) * v(i) * ( hL(i,j+1) + CFL * &
          (0.5*(hR(i,j+1)-hL(i,j+1)) + curv_3*(CFL - 1.5)) )
    else
      vh(i) = 0.0
    endif
  enddo
end subroutine merid_flux_En

!> Reflection of the internal waves at a single frequency.
subroutine reflect(En, NAngle, CS, G, LB)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure
  integer,                intent(in)    :: NAngle !< The number of wave orientations in the
                                              !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                          intent(inout) :: En !< The internal gravity wave energy density as a
                                              !! function of space and angular resolution
                                              !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  type(int_tide_CS),      intent(in)    :: CS !< Internal tide control structure
  type(loop_bounds_type), intent(in)    :: LB !< A structure with the active energy loop bounds.

  ! Local variables
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: angle_c
                                           ! angle of boundary wrt equator [rad]
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: part_refl
                                           ! fraction of wave energy reflected
                                           ! values should collocate with angle_c [nondim]
  logical, dimension(G%isd:G%ied,G%jsd:G%jed) :: ridge
                                           ! tags of cells with double reflection
  real, dimension(1:Nangle) :: En_reflected ! Energy reflected [H Z2 T-2 ~> m3 s-2 or J m-2].

  real    :: TwoPi                         ! 2*pi = 6.2831853... [nondim]
  real    :: Angle_size                    ! size of beam wedge [rad]
  integer :: angle_wall                    ! angle-bin of coast/ridge/shelf wrt equator
  integer :: angle_wall0                   ! angle-bin of coast/ridge/shelf wrt equator
  integer :: angle_r                       ! angle-bin of reflected ray wrt equator
  integer :: angle_r0                      ! angle-bin of reflected ray wrt equator
  integer :: angle_to_wall                 ! angle-bin relative to wall
  integer :: a, a0                         ! loop index for angles
  integer :: i, j
  integer :: Nangle_d2            ! Nangle / 2
  integer :: isc, iec, jsc, jec   ! start and end local indices on PE
                                  ! (values exclude halos)
  integer :: ish, ieh, jsh, jeh   ! start and end local indices on data domain
                                  ! leaving out outdated halo points (march in)

  isc = G%isc  ; iec = G%iec  ; jsc = G%jsc  ; jec = G%jec
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh

  TwoPi = 8.0*atan(1.0)
  Angle_size = TwoPi / (real(NAngle))
  Nangle_d2 = (Nangle / 2)

  ! init local arrays
  angle_c(:,:) = CS%nullangle
  part_refl(:,:) = 0.
  ridge(:,:) = .false.

  do j=jsh,jeh ; do i=ish,ieh
    if (CS%refl_angle(i,j) /= CS%nullangle) then
      angle_c(i,j) = mod(CS%refl_angle(i,j) + TwoPi, TwoPi)
    endif
    part_refl(i,j) = CS%refl_pref(i,j)
    ridge(i,j)     = CS%refl_dbl(i,j)
  enddo ; enddo
  En_reflected(:) = 0.0

  do j=jsh,jeh ; do i=ish,ieh
    ! redistribute energy in angular space if ray will hit boundary
    ! i.e., if energy is in a reflecting cell
    if (angle_c(i,j) /= CS%nullangle) then
      ! refection angle is given in rad, convert to the discrete angle
      angle_wall = nint(angle_c(i,j)/Angle_size) + 1
      do a=1,NAngle ; if (En(i,j,a) > 0.0) then
        ! reindex to 0 -> Nangle-1 for trig
        a0 = a - 1
        angle_wall0 = angle_wall - 1
        ! compute relative angle from wall and use cyclic properties
        ! to ensure it is bounded by 0 -> Nangle-1
        angle_to_wall = mod((a0 - angle_wall0) + Nangle, Nangle)

        if (ridge(i,j)) then
          ! if ray is not incident but in ridge cell, use complementary angle
          if ((Nangle_d2 < angle_to_wall) .and. (angle_to_wall < Nangle)) then
            angle_wall0 = mod(angle_wall0 + (Nangle_d2 + Nangle), Nangle)
          endif
        endif

        ! do reflection
        if ((0 < angle_to_wall) .and. (angle_to_wall < Nangle_d2)) then
          angle_r0 = mod(2*angle_wall0 - a0 + Nangle, Nangle)
          angle_r = angle_r0 + 1 !re-index to 1 -> Nangle
          if (a /= angle_r) then
            En_reflected(angle_r) = part_refl(i,j)*En(i,j,a)
            En(i,j,a) = (1.0-part_refl(i,j))*En(i,j,a)
          endif
        endif
      endif ; enddo ! a-loop
      do a=1,NAngle
        En(i,j,a) = En(i,j,a) + En_reflected(a)
        En_reflected(a) = 0.0  ! reset values
      enddo ! a-loop
    endif
  enddo ; enddo ! i- and j-loops

  ! Check to make sure no energy gets onto land (only run for debugging)
  ! do a=1,NAngle ; do j=jsc,jec ; do i=isc,iec
  !   if (En(i,j,a) > 0.001 .and. G%mask2dT(i,j) == 0) then
  !     write (mesg,*) 'En=', HZ2_T2_to_J_m2*En(i,j,a), 'a=', a, 'ig_g=',i+G%idg_offset, 'jg_g=',j+G%jdg_offset
  !     call MOM_error(FATAL, "reflect: Energy detected out of bounds: "//trim(mesg), .true.)
  !   endif
  ! enddo ; enddo ; enddo

end subroutine reflect

!> Moves energy across lines of partial reflection to prevent
!! reflection of energy that is supposed to get across.
subroutine teleport(En, NAngle, CS, G, LB)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure.
  integer,                intent(in)    :: NAngle !< The number of wave orientations in the
                                              !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                          intent(inout) :: En !< The internal gravity wave energy density as a
                                              !! function of space and angular resolution
                                              !! [H Z2 T-2 ~> m3 s-2 or J m-2].
  type(int_tide_CS),      intent(in)    :: CS !< Internal tide control structure
  type(loop_bounds_type), intent(in)    :: LB !< A structure with the active energy loop bounds.
  ! Local variables
  real, dimension(G%isd:G%ied,G%jsd:G%jed)    :: angle_c
                                              ! angle of boundary wrt equator [rad]
  real, dimension(G%isd:G%ied,G%jsd:G%jed)    :: part_refl
                                              ! fraction of wave energy reflected
                                              ! values should collocate with angle_c [nondim]
  logical, dimension(G%isd:G%ied,G%jsd:G%jed) :: pref_cell
                                              ! flag for partial reflection
  logical, dimension(G%isd:G%ied,G%jsd:G%jed) :: ridge
                                              ! tags of cells with double reflection
  real                        :: TwoPi      ! 2*pi = 6.2831853... [nondim]
  real                        :: Angle_size ! size of beam wedge [rad]
  real, dimension(1:NAngle)   :: angle_i    ! angle of incident ray wrt equator [rad]
  real, dimension(1:NAngle)   :: cos_angle  ! Cosine of the beam angle relative to eastward [nondim]
  real, dimension(1:NAngle)   :: sin_angle  ! Sine of the beam angle relative to eastward [nondim]
  real                        :: En_tele    ! energy to be "teleported" [H Z2 T-2 ~> m3 s-2 or J m-2]
  character(len=160) :: mesg  ! The text of an error message
  integer :: i, j, a
  integer :: ish, ieh, jsh, jeh     ! start and end local indices on data domain
                                    ! leaving out outdated halo points (march in)
  integer :: id_g, jd_g             ! global (decomposition-invariant) indices
  integer :: jos, ios               ! offsets
  real    :: cos_normal, sin_normal ! cos/sin of cross-ridge normal direction [nondim]
  real    :: angle_wall             ! The coastline angle or the complementary angle [radians]

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh

  TwoPi = 8.0*atan(1.0)
  Angle_size = TwoPi / (real(NAngle))

  do a=1,Nangle
    ! These are the angles at the cell centers
    ! (should do this elsewhere since doesn't change with time)
    angle_i(a) = Angle_size * real(a - 1) ! for a=1 aligned with x-axis
    cos_angle(a) = cos(angle_i(a)) ; sin_angle(a) = sin(angle_i(a))
  enddo

  angle_c   = CS%refl_angle
  part_refl = CS%refl_pref
  pref_cell = CS%refl_pref_logical
  ridge     = CS%refl_dbl

  do j=jsh,jeh
    do i=ish,ieh
      id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset
      if (pref_cell(i,j)) then
        do a=1,Nangle
          if (En(i,j,a) > 0) then
            ! if ray is incident, keep specified boundary angle
            if (sin(angle_i(a) - angle_c(i,j)) >= 0.0) then
              angle_wall = angle_c(i,j)
            ! if ray is not incident but in ridge cell, use complementary angle
            elseif (ridge(i,j)) then
              angle_wall = angle_c(i,j) + 0.5*TwoPi
            ! if ray is not incident and not in a ridge cell, keep specified angle
            else
              angle_wall = angle_c(i,j)
            endif
            ! teleport if incident
            if (sin(angle_i(a) - angle_wall) >= 0.0) then
              En_tele = En(i,j,a)
              cos_normal = cos(angle_wall + 0.25*TwoPi)
              sin_normal = sin(angle_wall + 0.25*TwoPi)
              ! find preferred zonal offset based on shelf/ridge angle
              ios = int(sign(1.,cos_normal))
              ! find preferred meridional offset based on shelf/ridge angle
              jos = int(sign(1.,sin_normal))
              ! find receptive ocean cell in direction of offset
              if (.not. pref_cell(i+ios,j+jos)) then
                En(i,j,a) = En(i,j,a) - En_tele
                En(i+ios,j+jos,a) = En(i+ios,j+jos,a) + En_tele
              else
                write(mesg,*) 'idg=',id_g,'jd_g=',jd_g,'a=',a
                call MOM_error(FATAL, "teleport: no receptive ocean cell at "//trim(mesg), .true.)
              endif
            endif ! incidence check
          endif ! energy check
        enddo ! a-loop
      endif ! pref check
    enddo ! i-loop
  enddo ! j-loop

end subroutine teleport

!> Rotates points in the halos where required to accommodate
!! changes in grid orientation, such as at the tripolar fold.
subroutine correct_halo_rotation(En, test, G, NAngle, halo)
  type(ocean_grid_type),      intent(in)    :: G    !< The ocean's grid structure
  real, dimension(:,:,:,:,:), intent(inout) :: En   !< The internal gravity wave energy density as a
                                       !! function of space, angular orientation, frequency,
                                       !! and vertical mode [H Z2 T-2 ~> m3 s-2 or J m-2].
  real, dimension(SZI_(G),SZJ_(G),2), &
                              intent(in)    :: test !< An x-unit vector that has been passed through
                                       !! the halo updates, to enable the rotation of the
                                       !! wave energies in the halo region to be corrected [nondim].
  integer,                    intent(in)    :: NAngle !< The number of wave orientations in the
                                                      !! discretized wave energy spectrum.
  integer,                    intent(in)    :: halo   !< The halo size over which to do the calculations
  ! Local variables
  real, dimension(G%isd:G%ied,NAngle) :: En2d ! A zonal row of the internal gravity wave energy density
                                              ! in a frequency band and mode [H Z2 T-2 ~> m3 s-2 or J m-2].
  integer, dimension(G%isd:G%ied) :: a_shift
  integer :: i_first, i_last, a_new
  integer :: a, i, j, ish, ieh, jsh, jeh, m, fr
  character(len=160) :: mesg  ! The text of an error message
  ish = G%isc-halo ; ieh = G%iec+halo ; jsh = G%jsc-halo ; jeh = G%jec+halo

  do j=jsh,jeh
    i_first = ieh+1 ; i_last = ish-1
    do i=ish,ieh
      a_shift(i) = 0
      if (test(i,j,1) /= 1.0) then
        if (i<i_first) i_first = i
        if (i>i_last) i_last = i

        if (test(i,j,1) == -1.0) then ; a_shift(i) = nAngle/2
        elseif (test(i,j,2) == 1.0) then ; a_shift(i) = -nAngle/4
        elseif (test(i,j,2) == -1.0) then ; a_shift(i) = nAngle/4
        else
          write(mesg,'("Unrecognized rotation test vector ",2ES9.2," at ",F7.2," E, ",&
                       &F7.2," N; i,j=",2i4)') &
                test(i,j,1), test(i,j,2), G%GeoLonT(i,j), G%GeoLatT(i,j), i, j
          call MOM_error(FATAL, mesg)
        endif
      endif
    enddo

    if (i_first <= i_last) then
      ! At least one point in this row needs to be rotated.
      do m=1,size(En,5) ; do fr=1,size(En,4)
        do a=1,nAngle ; do i=i_first,i_last ; if (a_shift(i) /= 0) then
          a_new = a + a_shift(i)
          if (a_new < 1) a_new = a_new + nAngle
          if (a_new > nAngle) a_new = a_new - nAngle
          En2d(i,a_new) = En(i,j,a,fr,m)
        endif ; enddo ; enddo
        do a=1,nAngle ; do i=i_first,i_last ; if (a_shift(i) /= 0) then
          En(i,j,a,fr,m) = En2d(i,a)
        endif ; enddo ; enddo
      enddo ; enddo
    endif
  enddo
end subroutine correct_halo_rotation

!> Calculates left/right edge values for PPM reconstruction in x-direction.
subroutine PPM_reconstruction_x(h_in, h_l, h_r, G, LB, simple_2nd)
  type(ocean_grid_type),            intent(in)  :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Energy density in a sector (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l  !< Left edge value of reconstruction (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r  !< Right edge value of reconstruction (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  type(loop_bounds_type),           intent(in)  :: LB   !< A structure with the active loop bounds.
  logical,                          intent(in)  :: simple_2nd !< If true, use the arithmetic mean
                                                        !! energy densities as default edge values
                                                        !! for a simple 2nd order scheme.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slope in energy density times the cell width
                                           ! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, parameter :: oneSixth = 1./6. ! One sixth [nondim]
  real :: h_ip1, h_im1 ! The energy densities at adjacent points [H Z2 T-2 ~> m3 s-2 or J m-2]
  real :: dMx, dMn ! The maximum and minimum of values of energy density at adjacent points
                   ! relative to the center point [H Z2 T-2 ~> m3 s-2 or J m-2]
  character(len=256) :: mesg  ! The text of an error message
  integer :: i, j, isl, iel, jsl, jel, stencil

  isl = LB%ish-1 ; iel = LB%ieh+1 ; jsl = LB%jsh ; jel = LB%jeh

  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (simple_2nd) stencil = 1

  if ((isl-stencil < G%isd) .or. (iel+stencil > G%ied)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_x called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               stencil + max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl < G%jsd) .or. (jel > G%jed)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_x called with a ", &
               & "y-halo that needs to be increased by ",i2,".")') &
               max(G%jsd-jsl,jel-G%jed)
    call MOM_error(FATAL,mesg)
  endif

  if (simple_2nd) then
    do j=jsl,jel ; do i=isl,iel
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_im1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) )
    enddo ; enddo
  else
    do j=jsl,jel ; do i=isl-1,iel+1
      if ((G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j)) == 0.0) then
        slp(i,j) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j) = 0.5 * (h_in(i+1,j) - h_in(i-1,j))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i+1,j), h_in(i-1,j), h_in(i,j)) - h_in(i,j)
        dMn = h_in(i,j) - min(h_in(i+1,j), h_in(i-1,j), h_in(i,j))
        slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j))
      endif
    enddo ; enddo

    do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
    ! h_im1 = h_in(i-1,j,k) ; if (G%mask2dT(i-1,j) < 0.5) h_im1 = h_in(i,j)
    ! h_ip1 = h_in(i+1,j,k) ; if (G%mask2dT(i+1,j) < 0.5) h_ip1 = h_in(i,j)
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_l(i,j) = 0.5*( h_im1 + h_in(i,j) ) + oneSixth*( slp(i-1,j) - slp(i,j) )
      h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i+1,j) )
    enddo ; enddo
  endif

  call PPM_limit_pos(h_in, h_l, h_r, 0.0, G, isl, iel, jsl, jel)
end subroutine PPM_reconstruction_x

!> Calculates left/right edge valus for PPM reconstruction in y-direction.
subroutine PPM_reconstruction_y(h_in, h_l, h_r, G, LB, simple_2nd)
  type(ocean_grid_type),            intent(in)  :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Energy density in a sector (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l  !< Left edge value of reconstruction (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r  !< Right edge value of reconstruction (2D)
                                                        !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  type(loop_bounds_type),           intent(in)  :: LB   !< A structure with the active loop bounds.
  logical,                          intent(in)  :: simple_2nd !< If true, use the arithmetic mean
                                                        !! energy densities as default edge values
                                                        !! for a simple 2nd order scheme.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slope in energy density times the cell width
                                           ! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, parameter :: oneSixth = 1./6. ! One sixth [nondim]
  real :: h_jp1, h_jm1 ! The energy densities at adjacent points [H Z2 T-2 ~> m3 s-2 or J m-2]
  real :: dMx, dMn ! The maximum and minimum of values of energy density at adjacent points
                   ! relative to the center point [H Z2 T-2 ~> m3 s-2 or J m-2]
  character(len=256) :: mesg  ! The text of an error message
  integer :: i, j, isl, iel, jsl, jel, stencil

  isl = LB%ish ; iel = LB%ieh ; jsl = LB%jsh-1 ; jel = LB%jeh+1

  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (simple_2nd) stencil = 1

  if ((isl < G%isd) .or. (iel > G%ied)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_y called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl-stencil < G%jsd) .or. (jel+stencil > G%jed)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_y called with a ", &
                 & "y-halo that needs to be increased by ",i2,".")') &
                 stencil + max(G%jsd-jsl,jel-G%jed)
    call MOM_error(FATAL,mesg)
  endif

  if (simple_2nd) then
    do j=jsl,jel ; do i=isl,iel
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) )
    enddo ; enddo
  else
    do j=jsl-1,jel+1 ; do i=isl,iel
      if ((G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1)) == 0.0) then
        slp(i,j) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j) = 0.5 * (h_in(i,j+1) - h_in(i,j-1))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i,j+1), h_in(i,j-1), h_in(i,j)) - h_in(i,j)
        dMn = h_in(i,j) - min(h_in(i,j+1), h_in(i,j-1), h_in(i,j))
        slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1))
      endif
    enddo ; enddo

    do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) ) + oneSixth*( slp(i,j-1) - slp(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i,j+1) )
    enddo ; enddo
  endif

  call PPM_limit_pos(h_in, h_l, h_r, 0.0, G, isl, iel, jsl, jel)
end subroutine PPM_reconstruction_y

!> Limits the left/right edge values of the PPM reconstruction
!! to give a reconstruction that is positive-definite.  Here this is
!! reinterpreted as giving a constant value if the mean value is less
!! than h_min, with a minimum of h_min otherwise.
subroutine PPM_limit_pos(h_in, h_L, h_R, h_min, G, iis, iie, jis, jie)
  type(ocean_grid_type),            intent(in)     :: G     !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)     :: h_in  !< Energy density in each sector (2D)
                                                            !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout)  :: h_L   !< Left edge value of reconstruction
                                                            !!  [H Z2 T-2 ~> m3 s-2 or J m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout)  :: h_R   !< Right edge value of reconstruction
                                                            !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  real,                             intent(in)     :: h_min !< The minimum value that can be
                                                            !! obtained by a concave parabolic fit
                                                            !! [H Z2 T-2 ~> m3 s-2 or J m-2]
  integer,                          intent(in)     :: iis   !< Start i-index for computations
  integer,                          intent(in)     :: iie   !< End i-index for computations
  integer,                          intent(in)     :: jis   !< Start j-index for computations
  integer,                          intent(in)     :: jie   !< End j-index for computations
  ! Local variables
  real    :: curv    ! The cell-area normalized curvature [H Z2 T-2 ~> m3 s-2 or J m-2]
  real    :: dh      ! The difference between the edge values [H Z2 T-2 ~> m3 s-2 or J m-2]
  real    :: scale   ! A rescaling factor used to give a minimum cell value of at least h_min [nondim]
  integer :: i, j

  do j=jis,jie ; do i=iis,iie
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    curv = 3.0*((h_L(i,j) + h_R(i,j)) - 2.0*h_in(i,j))
    if (curv > 0.0) then ! Only minima are limited.
      dh = h_R(i,j) - h_L(i,j)
      if (abs(dh) < curv) then ! The parabola's minimum is within the cell.
        if (h_in(i,j) <= h_min) then
          h_L(i,j) = h_in(i,j) ; h_R(i,j) = h_in(i,j)
        elseif (12.0*curv*(h_in(i,j) - h_min) < (curv**2 + 3.0*dh**2)) then
          ! The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
          ! be limited in this case.  0 < scale < 1.
          scale = 12.0*curv*(h_in(i,j) - h_min) / (curv**2 + 3.0*dh**2)
          h_L(i,j) = h_in(i,j) + scale*(h_L(i,j) - h_in(i,j))
          h_R(i,j) = h_in(i,j) + scale*(h_R(i,j) - h_in(i,j))
        endif
      endif
    endif
  enddo ; enddo
end subroutine PPM_limit_pos

subroutine register_int_tide_restarts(G, GV, US, param_file, CS, restart_CS)
  type(ocean_grid_type), intent(in) :: G          !< The ocean's grid structure
  type(verticalGrid_type),intent(in):: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type), intent(in) :: US         !< A dimensional unit scaling type
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(int_tide_CS),     pointer    :: CS         !< Internal tide control structure
  type(MOM_restart_CS),  pointer    :: restart_CS !< MOM restart control structure

  ! This subroutine is used to allocate and register any fields in this module
  ! that should be written to or read from the restart file.
  logical :: non_Bous          ! If true, this run is fully non-Boussinesq
  logical :: Boussinesq        ! If true, this run is fully Boussinesq
  logical :: semi_Boussinesq   ! If true, this run is partially non-Boussinesq
  logical :: use_int_tides
  integer :: num_freq, num_angle , num_mode, period_1
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, i, j, a, fr
  character(64) :: var_name, cfr, units

  type(axis_info) :: axes_inttides(2)
  real, dimension(:), allocatable :: angles, freqs ! Lables for angles and frequencies [nondim]
  real :: HZ2_T2_to_J_m2  ! unit conversion factor for Energy from internal to mks [H Z2 T-2 ~> m3 s-2 or J m-2]

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  HZ2_T2_to_J_m2 = GV%H_to_MKS*(US%Z_to_m**2)*(US%s_to_T**2)

  if (associated(CS)) then
    call MOM_error(WARNING, "register_int_tide_restarts called "//&
                             "with an associated control structure.")
    return
  endif

  allocate(CS)

  ! write extra axes
  call get_param(param_file, "MOM", "INTERNAL_TIDE_ANGLES", num_angle, default=24)
  call get_param(param_file, "MOM", "INTERNAL_TIDE_FREQS", num_freq, default=1)
  call get_param(param_file, "MOM", "INTERNAL_TIDE_MODES", num_mode, default=1)

  ! define restart units depemding on Boussinesq
  call get_param(param_file, "MOM", "BOUSSINESQ", Boussinesq, &
                 "If true, make the Boussinesq approximation.", default=.true., do_not_log=.true.)
  call get_param(param_file, "MOM", "SEMI_BOUSSINESQ", semi_Boussinesq, &
                 "If true, do non-Boussinesq pressure force calculations and use mass-based "//&
                 "thicknesses, but use RHO_0 to convert layer thicknesses into certain "//&
                 "height changes.  This only applies if BOUSSINESQ is false.", &
                 default=.true., do_not_log=.true.)
  non_Bous = .not.(Boussinesq .or. semi_Boussinesq)

  units = "J m-2"
  if (Boussinesq) units = "m3 s-2"

  allocate (angles(num_angle))
  allocate (freqs(num_freq))

  do a=1,num_angle ; angles(a) = a ; enddo
  do fr=1,num_freq ; freqs(fr) = fr ; enddo

  call set_axis_info(axes_inttides(1), "angle", "", "angle direction", num_angle, angles, "N", 1)
  call set_axis_info(axes_inttides(2), "freq", "", "wave frequency", num_freq, freqs, "N", 1)

  ! full energy array
  allocate(CS%En(isd:ied, jsd:jed, num_angle, num_freq, num_mode), source=0.0)

  ! restart strategy: support for 5d restart is not yet available so we split into
  ! 4d restarts. Vertical modes >= 6 are dissipated locally and do not propagate
  ! so we only allow for 5 vertical modes and each has its own variable

  ! allocate restart arrays
  allocate(CS%En_restart_mode1(isd:ied, jsd:jed, num_angle, num_freq), source=0.0)
  if (num_mode >= 2) allocate(CS%En_restart_mode2(isd:ied, jsd:jed, num_angle, num_freq), source=0.0)
  if (num_mode >= 3) allocate(CS%En_restart_mode3(isd:ied, jsd:jed, num_angle, num_freq), source=0.0)
  if (num_mode >= 4) allocate(CS%En_restart_mode4(isd:ied, jsd:jed, num_angle, num_freq), source=0.0)
  if (num_mode >= 5) allocate(CS%En_restart_mode5(isd:ied, jsd:jed, num_angle, num_freq), source=0.0)

  ! register all 4d restarts and copy into full Energy array when restarting from previous state
  call register_restart_field(CS%En_restart_mode1(:,:,:,:), "IW_energy_mode1", .false., restart_CS, &
                              longname="The internal wave energy density f(i,j,angle,freq) for mode 1", &
                              units=units, conversion=HZ2_T2_to_J_m2, z_grid='1', t_grid="s", &
                              extra_axes=axes_inttides)

  do fr=1,num_freq ; do a=1,num_angle ; do j=jsd,jed ; do i=isd,ied
    CS%En(i,j,a,fr,1) = CS%En_restart_mode1(i,j,a,fr)
  enddo ; enddo ; enddo ; enddo

  if (num_mode >= 2) then
    call register_restart_field(CS%En_restart_mode2(:,:,:,:), "IW_energy_mode2", .false., restart_CS, &
                                longname="The internal wave energy density f(i,j,angle,freq) for mode 2", &
                                units=units, conversion=HZ2_T2_to_J_m2, z_grid='1', t_grid="s", &
                                extra_axes=axes_inttides)

    do fr=1,num_freq ; do a=1,num_angle ; do j=jsd,jed ; do i=isd,ied
      CS%En(i,j,a,fr,2) = CS%En_restart_mode2(i,j,a,fr)
    enddo ; enddo ; enddo ; enddo

  endif

  if (num_mode >= 3) then
    call register_restart_field(CS%En_restart_mode3(:,:,:,:), "IW_energy_mode3", .false., restart_CS, &
                                longname="The internal wave energy density f(i,j,angle,freq) for mode 3", &
                                units=units, conversion=HZ2_T2_to_J_m2, z_grid='1', t_grid="s", &
                                extra_axes=axes_inttides)

    do fr=1,num_freq ; do a=1,num_angle ; do j=jsd,jed ; do i=isd,ied
      CS%En(i,j,a,fr,3) = CS%En_restart_mode3(i,j,a,fr)
    enddo ; enddo ; enddo ; enddo

  endif

  if (num_mode >= 4) then
    call register_restart_field(CS%En_restart_mode4(:,:,:,:), "IW_energy_mode4", .false., restart_CS, &
                                longname="The internal wave energy density f(i,j,angle,freq) for mode 4", &
                                units=units, conversion=HZ2_T2_to_J_m2, z_grid='1', t_grid="s", &
                                extra_axes=axes_inttides)

    do fr=1,num_freq ; do a=1,num_angle ; do j=jsd,jed ; do i=isd,ied
      CS%En(i,j,a,fr,4) = CS%En_restart_mode4(i,j,a,fr)
    enddo ; enddo ; enddo ; enddo

  endif

  if (num_mode >= 5) then
    call register_restart_field(CS%En_restart_mode5(:,:,:,:), "IW_energy_mode5", .false., restart_CS, &
                                longname="The internal wave energy density f(i,j,angle,freq) for mode 5", &
                                units=units, conversion=HZ2_T2_to_J_m2, z_grid='1', t_grid="s", &
                                extra_axes=axes_inttides)

    do fr=1,num_freq ; do a=1,num_angle ; do j=jsd,jed ; do i=isd,ied
      CS%En(i,j,a,fr,5) = CS%En_restart_mode5(i,j,a,fr)
    enddo ; enddo ; enddo ; enddo

  endif

end subroutine register_int_tide_restarts

!> This subroutine initializes the internal tides module.
subroutine internal_tides_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target,   intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),     intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),     intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),     intent(in)    :: param_file !< A structure to parse for run-time
                                                   !! parameters.
  type(diag_ctrl), target,   intent(in)    :: diag !< A structure that is used to regulate
                                                   !! diagnostic output.
  type(int_tide_CS),         pointer :: CS         !< Internal tide control structure

  ! Local variables
  real                              :: Angle_size ! size of wedges [rad]
  real, allocatable                 :: angles(:)  ! orientations of wedge centers [rad]
  real, dimension(:,:), allocatable :: h2         ! topographic roughness scale squared [Z2 ~> m2]
  real                              :: kappa_itides ! characteristic topographic wave number [L-1 ~> m-1]
  real, dimension(:,:), allocatable :: ridge_temp ! array for temporary storage of flags
                                                  ! of cells with double-reflecting ridges [nondim]
  real, dimension(:,:), allocatable :: tmp_decay ! a temp array to store decay rates [T-1 ~> s-1]
  real :: decay_rate                             ! A constant rate at which internal tide energy is
                                                 ! lost to the interior ocean internal wave field [T-1 ~> s-1].
  logical :: use_int_tides, use_temperature
  logical :: om4_remap_via_sub_cells ! Use the OM4-era ramap_via_sub_cells for calculating the EBT structure
  real    :: IGW_c1_thresh ! A threshold first mode internal wave speed below which all higher
                 ! mode speeds are not calculated but simply assigned a speed of 0 [L T-1 ~> m s-1].
  real    :: kappa_h2_factor    ! A roughness scaling factor [nondim]
  real    :: RMS_roughness_frac ! The maximum RMS topographic roughness as a fraction of the
                                ! nominal ocean depth, or a negative value for no limit [nondim]
  real    :: period_1           ! The period of the gravest modeled mode [T ~> s]
  real    :: period             ! A tidal period read from namelist [T ~> s]
  real    :: HZ2_T2_to_J_m2     ! unit conversion factor for Energy from internal units
                                ! to mks [T2 kg H-1 Z-2 s-2 ~> kg m-3 or 1]
  real    :: HZ2_T3_to_W_m2     ! unit conversion factor for TKE from internal units
                                ! to mks [T3 kg H-1 Z-2 s-3 ~> kg m-3 or 1]
  real    :: J_m2_to_HZ2_T2     ! unit conversion factor for Energy from mks to internal
                                ! units [H Z2 s2 T-2 kg-1 ~> m3 kg-1 or 1]
  integer :: num_angle, num_freq, num_mode, m, fr
  integer :: isd, ied, jsd, jed, a, id_ang, i, j, nz
  type(axes_grp) :: axes_ang
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_internal_tides" ! This module's name.
  character(len=16), dimension(8) :: freq_name
  character(len=40)  :: var_name
  character(len=160) :: var_descript
  character(len=200) :: filename
  character(len=200) :: refl_angle_file
  character(len=200) :: refl_pref_file, refl_dbl_file, trans_file
  character(len=200) :: h2_file, decay_file
  character(len=80)  :: rough_var ! Input file variable names

  character(len=240), dimension(:), allocatable :: energy_fractions
  character(len=240) :: periods

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  nz = GV%ke

  HZ2_T2_to_J_m2 = GV%H_to_kg_m2*(US%Z_to_m**2)*(US%s_to_T**2)
  HZ2_T3_to_W_m2 = GV%H_to_kg_m2*(US%Z_to_m**2)*(US%s_to_T**3)
  J_m2_to_HZ2_T2 = GV%kg_m2_to_H*(US%m_to_Z**2)*(US%T_to_s**2)

  CS%initialized = .true.

  use_int_tides = .false.
  call read_param(param_file, "INTERNAL_TIDES", use_int_tides)
  CS%do_int_tides = use_int_tides
  if (.not.use_int_tides) return

  use_temperature = .true.
  call read_param(param_file, "ENABLE_THERMODYNAMICS", use_temperature)
  if (.not.use_temperature) call MOM_error(FATAL, &
      "internal_tides_init: internal_tides only works with ENABLE_THERMODYNAMICS defined.")

  ! Set number of frequencies, angles, and modes to consider
  num_freq = 1 ; num_angle = 24 ; num_mode = 1
  call read_param(param_file, "INTERNAL_TIDE_FREQS", num_freq)
  call read_param(param_file, "INTERNAL_TIDE_ANGLES", num_angle)
  call read_param(param_file, "INTERNAL_TIDE_MODES", num_mode)
  if (.not.((num_freq > 0) .and. (num_angle > 0) .and. (num_mode > 0))) return
  CS%nFreq = num_freq ; CS%nAngle = num_angle ; CS%nMode = num_mode

  allocate(energy_fractions(num_freq))
  allocate(CS%fraction_tidal_input(num_freq,num_mode))

  call read_param(param_file, "ENERGY_FRACTION_PER_MODE", energy_fractions)

  do fr=1,num_freq ; do m=1,num_mode
    CS%fraction_tidal_input(fr,m) = extract_real(energy_fractions(fr), " ,", m, 0.)
  enddo ; enddo

  ! Allocate phase speed array
  allocate(CS%cp(isd:ied, jsd:jed, num_freq, num_mode), source=0.0)

  ! Allocate and populate frequency array (each a multiple of first for now)
  allocate(CS%frequency(num_freq))


  call get_param(param_file, mdl, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)

  ! The periods of the tidal constituents for internal tides raytracing
  call read_param(param_file, "TIDAL_PERIODS", periods)

  do fr=1,num_freq
    period = US%s_to_T*extract_real(periods, " ,", fr, 0.)
    if (period == 0.) call MOM_error(FATAL, "MOM_internal_tides: invalid tidal period")
    CS%frequency(fr) = 8.0*atan(1.0)/period
  enddo

  ! Read all relevant parameters and write them to the model log.

  CS%Time => Time ! direct a pointer to the current model time target

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".")
                 CS%inputdir = slasher(CS%inputdir)

  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "INTERNAL_TIDE_FREQS", num_freq, &
                 "The number of distinct internal tide frequency bands "//&
                 "that will be calculated.", default=1)
  call get_param(param_file, mdl, "INTERNAL_TIDE_MODES", num_mode, &
                 "The number of distinct internal tide modes "//&
                 "that will be calculated.", default=1)
  call get_param(param_file, mdl, "INTERNAL_TIDE_ANGLES", num_angle, &
                 "The number of angular resolution bands for the internal "//&
                 "tide calculations.", default=24)

  if (use_int_tides) then
    if ((num_freq <= 0) .and. (num_mode <= 0) .and. (num_angle <= 0)) then
      call MOM_error(WARNING, "Internal tides were enabled, but the number "//&
             "of requested frequencies, modes and angles were not all positive.")
      return
    endif
  else
    if ((num_freq > 0) .and. (num_mode > 0) .and. (num_angle > 0)) then
      call MOM_error(WARNING, "Internal tides were not enabled, even though "//&
             "the number of requested frequencies, modes and angles were all "//&
             "positive.")
      return
    endif
  endif

  if (CS%NFreq /= num_freq) call MOM_error(FATAL, "Internal_tides_init: "//&
      "Inconsistent number of frequencies.")
  if (CS%NAngle /= num_angle) call MOM_error(FATAL, "Internal_tides_init: "//&
      "Inconsistent number of angles.")
  if (CS%nMode /= num_mode) call MOM_error(FATAL, "Internal_tides_init: "//&
      "Inconsistent number of modes.")
  if (4*(num_angle/4) /= num_angle) call MOM_error(FATAL, &
    "Internal_tides_init: INTERNAL_TIDE_ANGLES must be a multiple of 4.")

  CS%diag => diag

  call get_param(param_file, mdl, "INTERNAL_TIDES_UPDATE_KD", CS%update_Kd, &
                 "If true, internal tides ray tracing changes Kd for dynamics.", &
                 default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDES_REFRACTION", CS%apply_refraction, &
                 "If true, internal tides ray tracing does refraction.", &
                 default=.true.)
  call get_param(param_file, mdl, "INTERNAL_TIDES_PROPAGATION", CS%apply_propagation, &
                 "If true, internal tides ray tracing does propagate.", &
                 default=.true.)
  call get_param(param_file, mdl, "INTERNAL_TIDES_ONLY_INIT_FORCING", CS%init_forcing_only, &
                 "If true, internal tides ray tracing only applies forcing at first step (debugging).", &
                 default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDES_FORCE_POS_EN", CS%force_posit_En, &
                 "If true, force energy to be positive by removing subroundoff negative values.", &
                 default=.true.)
  call get_param(param_file, mdl, "KD_MIN", CS%Kd_min, &
                 "The minimum diapycnal diffusivity.", &
                 units="m2 s-1", default=2e-6, scale=GV%m2_s_to_HZ_T)
  call get_param(param_file, mdl, "MINTHICK_TKE_TO_KD", CS%min_thick_layer_Kd, &
                 "The minimum thickness allowed with TKE_to_Kd.", &
                 units="m", default=1e-6, scale=GV%m_to_H)
  call get_param(param_file, mdl, "ITIDES_MIXING_EFFIC", CS%mixing_effic, &
                 "Mixing efficiency for internal tides raytracing", &
                 units="nondim", default=0.2)
  call get_param(param_file, mdl, "MAX_TKE_TO_KD", CS%max_TKE_to_Kd, &
                 "Limiter for TKE_to_Kd.", &
                 units="s2 m-1", default=1e9, scale=US%Z_to_m*US%s_to_T**2)
  call get_param(param_file, mdl, "INTERNAL_TIDE_DECAY_RATE", decay_rate, &
                 "The rate at which internal tide energy is lost to the "//&
                 "interior ocean internal wave field.", &
                 units="s-1", default=0.0, scale=US%T_to_s)
  call get_param(param_file, mdl, "USE_2D_INTERNAL_TIDE_DECAY_RATE", CS%use_2d_decay_rate, &
                 "If true, use a spatially varying decay rate for leakage loss in the "// &
                 "internal tide code.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_VOLUME_BASED_CFL", CS%vol_CFL, &
                 "If true, use the ratio of the open face lengths to the "//&
                 "tracer cell areas when estimating CFL numbers in the "//&
                 "internal tide code.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_SIMPLE_2ND_PPM", CS%simple_2nd, &
                 "If true, CONTINUITY_PPM uses a simple 2nd order "//&
                 "(arithmetic mean) interpolation of the edge values. "//&
                 "This may give better PV conservation properties. While "//&
                 "it formally reduces the accuracy of the continuity "//&
                 "solver itself in the strongly advective limit, it does "//&
                 "not reduce the overall order of accuracy of the dynamic "//&
                 "core.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_UPWIND_1ST", CS%upwind_1st, &
                 "If true, the internal tide ray-tracing advection uses "//&
                 "1st-order upwind advection.  This scheme is highly "//&
                 "continuity solver.  This scheme is highly "//&
                 "diffusive but may be useful for debugging.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_BACKGROUND_DRAG", CS%apply_background_drag, &
                 "If true, the internal tide ray-tracing advection uses a background drag "//&
                 "term as a sink.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_QUAD_DRAG", CS%apply_bottom_drag, &
                 "If true, the internal tide ray-tracing advection uses "//&
                 "a quadratic bottom drag term as a sink.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_WAVE_DRAG", CS%apply_wave_drag, &
                 "If true, apply scattering due to small-scale roughness as a sink.", &
                 default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_RESIDUAL_DRAG", CS%apply_residual_drag, &
                 "If true, apply drag due to critical slopes", &
                 default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_DRAG_MIN_DEPTH", CS%drag_min_depth, &
                 "The minimum total ocean thickness that will be used in the denominator "//&
                 "of the quadratic drag terms for internal tides.", &
                 units="m", default=1.0, scale=GV%m_to_H, do_not_log=.not.CS%apply_bottom_drag)
  CS%drag_min_depth = MAX(CS%drag_min_depth, GV%H_subroundoff)
  call get_param(param_file, mdl, "INTERNAL_TIDE_FROUDE_DRAG", CS%apply_Froude_drag, &
                 "If true, apply wave breaking as a sink.", &
                 default=.false.)
  call get_param(param_file, mdl, "EN_CHECK_TOLERANCE", CS%En_check_tol, &
                 "An energy density tolerance for flagging points with small negative "//&
                 "internal tide energy.", &
                 units="J m-2", default=1.0, scale=J_m2_to_HZ2_T2, &
                 do_not_log=.not.CS%apply_Froude_drag)
  call get_param(param_file, mdl, "EN_UNDERFLOW", CS%En_underflow, &
                 "A small energy density below which Energy is set to zero.", &
                 units="J m-2", default=1.0e-100, scale=J_m2_to_HZ2_T2)
  call get_param(param_file, mdl, "EN_RESTART_POWER", CS%En_restart_power, &
                 "A power factor to save larger values x 2**(power) in restart files.", &
                 units="nondim", default=0)
  call get_param(param_file, mdl, "CDRAG", CS%cdrag, &
                 "CDRAG is the drag coefficient relating the magnitude of "//&
                 "the velocity field to the bottom stress.", &
                 units="nondim", default=0.003)
  call get_param(param_file, mdl, "INTERNAL_WAVE_CG1_THRESH", IGW_c1_thresh, &
                 "A minimal value of the first mode internal wave speed below which all higher "//&
                 "mode speeds are not calculated but are simply reported as 0.  This must be "//&
                 "non-negative for the wave_speeds routine to be used.", &
                 units="m s-1", default=0.01, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 do_not_log=.true., default=.true.)
  call get_param(param_file, mdl, "INTWAVE_REMAPPING_USE_OM4_SUBCELLS", om4_remap_via_sub_cells, &
                 "If true, use the OM4 remapping-via-subcells algorithm for calculating EBT structure. "//&
                 "See REMAPPING_USE_OM4_SUBCELLS for details. "//&
                 "We recommend setting this option to false.", default=om4_remap_via_sub_cells)
  call get_param(param_file, mdl, "UNIFORM_TEST_CG", CS%uniform_test_cg, &
                 "If positive, a uniform group velocity of internal tide for test case", &
                 default=-1., units="m s-1", scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "INTERNAL_TIDE_ENERGIZED_ANGLE", CS%energized_angle, &
                 "If positive, only one angular band of the internal tides "//&
                 "gets all of the energy.  (This is for debugging.)", default=-1)
  call get_param(param_file, mdl, "USE_PPM_ANGULAR", CS%use_PPMang, &
                 "If true, use PPM for advection of energy in angular space.", &
                 default=.false.)
  call get_param(param_file, mdl, "GAMMA_ITIDES", CS%q_itides, &
                 "The fraction of the internal tidal energy that is "//&
                 "dissipated locally with INT_TIDE_DISSIPATION. "//&
                 "THIS NAME COULD BE BETTER.", &
                 units="nondim", default=0.3333)
  call get_param(param_file, mdl, "KAPPA_ITIDES", kappa_itides, &
               "A topographic wavenumber used with INT_TIDE_DISSIPATION. "//&
               "The default is 2pi/10 km, as in St.Laurent et al. 2002.", &
               units="m-1", default=8.e-4*atan(1.0), scale=US%L_to_m)
  call get_param(param_file, mdl, "KAPPA_H2_FACTOR", kappa_h2_factor, &
               "A scaling factor for the roughness amplitude with "//&
               "INT_TIDE_DISSIPATION.",  units="nondim", default=1.0)
  call get_param(param_file, mdl, "GAMMA_OSBORN", CS%gamma_osborn, &
               "The mixing efficiency for internan tides from Osborn 1980 ", &
               units="nondim", default=0.2)
  call get_param(param_file, mdl, "INT_TIDE_DECAY_SCALE", CS%Int_tide_decay_scale, &
                 "The decay scale away from the bottom for tidal TKE with "//&
                 "the new coding when INT_TIDE_DISSIPATION is used.", &
                 units="m", default=500.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "INT_TIDE_DECAY_SCALE_SLOPES", CS%Int_tide_decay_scale_slope, &
                 "The slope decay scale away from the bottom for tidal TKE with "//&
                 "the new coding when INT_TIDE_DISSIPATION is used.", &
                 units="m", default=100.0, scale=GV%m_to_H)

  ! Allocate various arrays needed for loss rates
  allocate(h2(isd:ied,jsd:jed), source=0.0)
  allocate(CS%TKE_itidal_loss_fixed(isd:ied,jsd:jed), source=0.0)
  allocate(CS%TKE_leak_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode), source=0.0)
  allocate(CS%TKE_quad_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode), source=0.0)
  allocate(CS%TKE_itidal_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode), source=0.0)
  allocate(CS%TKE_Froude_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode), source=0.0)
  allocate(CS%TKE_residual_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode), source=0.0)
  allocate(CS%TKE_slope_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode), source=0.0)
  allocate(CS%tot_leak_loss(isd:ied,jsd:jed), source=0.0)
  allocate(CS%tot_quad_loss(isd:ied,jsd:jed), source=0.0)
  allocate(CS%tot_itidal_loss(isd:ied,jsd:jed), source=0.0)
  allocate(CS%tot_Froude_loss(isd:ied,jsd:jed), source=0.0)
  allocate(CS%tot_residual_loss(isd:ied,jsd:jed), source=0.0)
  allocate(CS%u_struct_bot(isd:ied,jsd:jed,num_mode), source=0.0)
  allocate(CS%u_struct_max(isd:ied,jsd:jed,num_mode), source=0.0)
  allocate(CS%int_w2(isd:ied,jsd:jed,num_mode), source=0.0)
  allocate(CS%int_U2(isd:ied,jsd:jed,num_mode), source=0.0)
  allocate(CS%int_N2w2(isd:ied,jsd:jed,num_mode), source=0.0)
  allocate(CS%w_struct(isd:ied,jsd:jed,1:nz+1,num_mode), source=0.0)
  allocate(CS%u_struct(isd:ied,jsd:jed,1:nz,num_mode), source=0.0)
  allocate(CS%error_mode(num_freq,num_mode), source=0.0)
  allocate(CS%En_ini_glo(num_freq,num_mode), source=0.0)
  allocate(CS%En_end_glo(num_freq,num_mode), source=0.0)
  allocate(CS%TKE_leak_loss_glo_dt(num_freq,num_mode), source=0.0)
  allocate(CS%TKE_quad_loss_glo_dt(num_freq,num_mode), source=0.0)
  allocate(CS%TKE_Froude_loss_glo_dt(num_freq,num_mode), source=0.0)
  allocate(CS%TKE_itidal_loss_glo_dt(num_freq,num_mode), source=0.0)
  allocate(CS%TKE_residual_loss_glo_dt(num_freq,num_mode), source=0.0)
  allocate(CS%TKE_input_glo_dt(num_freq,num_mode), source=0.0)
  allocate(CS%decay_rate_2d(isd:ied,jsd:jed,num_freq,num_mode), source=0.0)
  allocate(tmp_decay(isd:ied,jsd:jed), source=0.0)

  if (CS%use_2d_decay_rate) then
    call get_param(param_file, mdl, "ITIDES_DECAY_FILE", decay_file, &
            "The path to the file containing the decay rates "//&
            "for internal tides with USE_2D_INTERNAL_TIDE_DECAY_RATE.", &
            fail_if_missing=.true.)
    do m=1,num_mode ; do fr=1,num_freq
      ! read 2d field for each harmonic
      filename = trim(CS%inputdir) // trim(decay_file)
      write(var_name, '("decay_rate_freq",i1,"_mode",i1)') fr, m
      call MOM_read_data(filename, var_name, tmp_decay, G%domain, scale=US%T_to_s)
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        CS%decay_rate_2d(i,j,fr,m) = tmp_decay(i,j)
      enddo ; enddo
    enddo ; enddo
  else
    do m=1,num_mode ; do fr=1,num_freq ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      CS%decay_rate_2d(i,j,fr,m) = decay_rate
    enddo ; enddo ; enddo ; enddo
  endif

  do m=1,num_mode
    call pass_var(CS%decay_rate_2d(:,:,:,m), G%domain)
  enddo

  ! Compute the fixed part of the bottom drag loss from baroclinic modes
  call get_param(param_file, mdl, "H2_FILE", h2_file, &
          "The path to the file containing the sub-grid-scale "//&
          "topographic roughness amplitude with INT_TIDE_DISSIPATION.", &
          fail_if_missing=.true.)
  filename = trim(CS%inputdir) // trim(h2_file)
  call log_param(param_file, mdl, "INPUTDIR/H2_FILE", filename)
  call get_param(param_file, mdl, "ROUGHNESS_VARNAME", rough_var, &
                 "The name in the input file of the squared sub-grid-scale "//&
                 "topographic roughness amplitude variable.", default="h2")
  call get_param(param_file, mdl, "INTERNAL_TIDE_ROUGHNESS_FRAC", RMS_roughness_frac, &
                 "The maximum RMS topographic roughness as a fraction of the nominal ocean depth, "//&
                 "or a negative value for no limit.",  units="nondim", default=0.1)

  call MOM_read_data(filename, rough_var, h2, G%domain, scale=US%m_to_Z**2)
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! Restrict RMS topographic roughness to a fraction (10 percent by default) of the column depth.
    if (RMS_roughness_frac >= 0.0) then
      h2(i,j) = max(min((RMS_roughness_frac*(G%bathyT(i,j)+G%Z_ref))**2, h2(i,j)), 0.0)
    else
      h2(i,j) = max(h2(i,j), 0.0)
    endif
    ! Compute the fixed part; units are [R Z4 H-1 L-2 ~> kg m-2 or m] here
    ! will be multiplied by N and the squared near-bottom velocity (and by the
    ! near-bottom density in non-Boussinesq mode) to get into [H Z2 T-3 ~> m3 s-3 or W m-2]
    CS%TKE_itidal_loss_fixed(i,j) = 0.5*kappa_h2_factor* GV%H_to_RZ * US%L_to_Z*kappa_itides * h2(i,j)
  enddo ; enddo

  deallocate(h2)

  ! Read in prescribed coast/ridge/shelf angles from file
  call get_param(param_file, mdl, "REFL_ANGLE_FILE", refl_angle_file, &
               "The path to the file containing the local angle of "//&
               "the coastline/ridge/shelf with respect to the equator.", &
               fail_if_missing=.false., default='')
  filename = trim(CS%inputdir) // trim(refl_angle_file)
  allocate(CS%refl_angle(isd:ied,jsd:jed), source=CS%nullangle)
  if (file_exists(filename, G%domain)) then
    call log_param(param_file, mdl, "INPUTDIR/REFL_ANGLE_FILE", filename)
    call MOM_read_data(filename, 'refl_angle', CS%refl_angle, G%domain)
  else
    if (trim(refl_angle_file) /= '' ) call MOM_error(FATAL, &
                                                     "REFL_ANGLE_FILE: "//trim(filename)//" not found")
  endif
  ! replace NaNs with null value
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if (is_NaN(CS%refl_angle(i,j))) CS%refl_angle(i,j) = CS%nullangle
  enddo ; enddo
  call pass_var(CS%refl_angle, G%domain)

  ! Read in prescribed partial reflection coefficients from file
  call get_param(param_file, mdl, "REFL_PREF_FILE", refl_pref_file, &
               "The path to the file containing the reflection coefficients.", &
               fail_if_missing=.false., default='')
  filename = trim(CS%inputdir) // trim(refl_pref_file)
  allocate(CS%refl_pref(isd:ied,jsd:jed), source=1.0)
  if (file_exists(filename, G%domain)) then
    call log_param(param_file, mdl, "INPUTDIR/REFL_PREF_FILE", filename)
    call MOM_read_data(filename, 'refl_pref', CS%refl_pref, G%domain)
  else
    if (trim(refl_pref_file) /= '' ) call MOM_error(FATAL, &
                                                    "REFL_PREF_FILE: "//trim(filename)//" not found")
  endif
  !CS%refl_pref = CS%refl_pref*1 ! adjust partial reflection if desired
  call pass_var(CS%refl_pref, G%domain)

  ! Tag reflection cells with partial reflection (done here for speed)
  allocate(CS%refl_pref_logical(isd:ied,jsd:jed), source=.false.)
  do j=jsd,jed ; do i=isd,ied
    ! flag cells with partial reflection
    if ((CS%refl_angle(i,j) /= CS%nullangle) .and. &
        (CS%refl_pref(i,j) < 1.0) .and. (CS%refl_pref(i,j) > 0.0)) then
      CS%refl_pref_logical(i,j) = .true.
    endif
  enddo ; enddo

  ! Read in double-reflective (ridge) tags from file
  call get_param(param_file, mdl, "REFL_DBL_FILE", refl_dbl_file, &
               "The path to the file containing the double-reflective ridge tags.", &
               fail_if_missing=.false., default='')
  filename = trim(CS%inputdir) // trim(refl_dbl_file)
  allocate(ridge_temp(isd:ied,jsd:jed), source=0.0)
  if (file_exists(filename, G%domain)) then
    call log_param(param_file, mdl, "INPUTDIR/REFL_DBL_FILE", filename)
    call MOM_read_data(filename, 'refl_dbl', ridge_temp, G%domain)
  else
    if (trim(refl_dbl_file) /= '' ) call MOM_error(FATAL, &
                                                   "REFL_DBL_FILE: "//trim(filename)//" not found")
  endif
  call pass_var(ridge_temp, G%domain)
  allocate(CS%refl_dbl(isd:ied,jsd:jed), source=.false.)
  do j=jsd,jed ; do i=isd,ied
    CS%refl_dbl(i,j) = (ridge_temp(i,j) == 1)
  enddo ; enddo

  ! Read in the transmission coefficient and infer the residual
  call get_param(param_file, mdl, "TRANS_FILE", trans_file, &
               "The path to the file containing the transmission coefficent for internal tides.", &
               fail_if_missing=.false., default='')
  filename = trim(CS%inputdir) // trim(trans_file)
  allocate(CS%trans(isd:ied,jsd:jed), source=0.0)
  if (file_exists(filename, G%domain)) then
    call log_param(param_file, mdl, "INPUTDIR/TRANS_FILE", filename)
    call MOM_read_data(filename, 'trans', CS%trans, G%domain)
  else
    if (trim(trans_file) /= '' ) call MOM_error(FATAL, &
                                                "TRANS_FILE: "//trim(filename)//" not found")
  endif

  call pass_var(CS%trans, G%domain)

  ! residual
  allocate(CS%residual(isd:ied,jsd:jed), source=0.0)
  if (CS%apply_residual_drag) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (CS%refl_pref_logical(i,j)) then
        CS%residual(i,j) = 1. - (CS%refl_pref(i,j) - CS%trans(i,j))
      endif
    enddo ; enddo
    call pass_var(CS%residual, G%domain)
  else
    ! report residual of transmission/reflection onto reflection
    ! this ensure energy budget is conserved
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (CS%refl_pref_logical(i,j)) then
        CS%refl_pref(i,j) = 1. - CS%trans(i,j)
      endif
    enddo ; enddo
    call pass_var(CS%refl_pref, G%domain)
  endif

  CS%id_cg1 = register_diag_field('ocean_model', 'cn1', diag%axesT1, &
               Time, 'First baroclinic mode (eigen) speed', 'm s-1', conversion=US%L_T_to_m_s)
  allocate(CS%id_cn(CS%nMode), source=-1)
  do m=1,CS%nMode
    write(var_name, '("cn_mode",i1)') m
    write(var_descript, '("Baroclinic (eigen) speed of mode ",i1)') m
    CS%id_cn(m) = register_diag_field('ocean_model',var_name, diag%axesT1, &
                 Time, var_descript, 'm s-1', conversion=US%L_T_to_m_s)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
  enddo

  ! Register maps of reflection parameters
  CS%id_refl_ang = register_diag_field('ocean_model', 'refl_angle', diag%axesT1, &
                 Time, 'Local angle of coastline/ridge/shelf with respect to equator', 'rad')
  CS%id_refl_pref = register_diag_field('ocean_model', 'refl_pref', diag%axesT1, &
                 Time, 'Partial reflection coefficients', '')
  CS%id_trans = register_diag_field('ocean_model', 'trans', diag%axesT1, &
                 Time, 'Partial transmission coefficients', '')
  CS%id_residual = register_diag_field('ocean_model', 'residual', diag%axesT1, &
                 Time, 'Residual of reflection and transmission coefficients', '')
  CS%id_dx_Cv = register_diag_field('ocean_model', 'dx_Cv', diag%axesT1, &
                 Time, 'North face unblocked width', 'm', conversion=US%L_to_m)
  CS%id_dy_Cu = register_diag_field('ocean_model', 'dy_Cu', diag%axesT1, &
                 Time, 'East face unblocked width', 'm', conversion=US%L_to_m)
  CS%id_land_mask = register_diag_field('ocean_model', 'land_mask', diag%axesT1, &
                 Time, 'Land mask', 'nondim')
  ! Output reflection parameters as diagnostics here (not needed every timestep)
  if (CS%id_refl_ang > 0)   call post_data(CS%id_refl_ang, CS%refl_angle, CS%diag)
  if (CS%id_refl_pref > 0)  call post_data(CS%id_refl_pref, CS%refl_pref, CS%diag)
  if (CS%id_trans > 0)      call post_data(CS%id_trans, CS%trans, CS%diag)
  if (CS%id_residual > 0)   call post_data(CS%id_residual, CS%residual, CS%diag)
  if (CS%id_dx_Cv > 0)      call post_data(CS%id_dx_Cv, G%dx_Cv, CS%diag)
  if (CS%id_dy_Cu > 0)      call post_data(CS%id_dy_Cu, G%dy_Cu, CS%diag)
  if (CS%id_land_mask > 0)  call post_data(CS%id_land_mask, G%mask2dT, CS%diag)

  ! Register 2-D energy density (summed over angles, freq, modes)
  CS%id_tot_En = register_diag_field('ocean_model', 'ITide_tot_En', diag%axesT1, &
                 Time, 'Internal tide total energy density', &
                 'J m-2', conversion=HZ2_T2_to_J_m2)

  allocate(CS%id_itide_drag(CS%nFreq, CS%nMode), source=-1)
  allocate(CS%id_TKE_itidal_input(CS%nFreq), source=-1)
  do fr=1,CS%nFreq
    ! Register 2-D energy input into internal tides for each frequency
    write(var_name, '("TKE_itidal_input_freq",i1)') fr
    write(var_descript, '("a fraction of which goes into rays in frequency ",i1)') fr

    CS%id_TKE_itidal_input(fr) = register_diag_field('ocean_model', var_name, diag%axesT1, &
                                                     Time, 'Conversion from barotropic to baroclinic tide, '//&
                                                     var_descript, 'W m-2', conversion=HZ2_T3_to_W_m2)
  enddo
  ! Register 2-D energy losses (summed over angles, freq, modes)
  CS%id_tot_leak_loss = register_diag_field('ocean_model', 'ITide_tot_leak_loss', diag%axesT1, &
                Time, 'Internal tide energy loss to background drag', &
                'W m-2', conversion=HZ2_T3_to_W_m2)
  CS%id_tot_quad_loss = register_diag_field('ocean_model', 'ITide_tot_quad_loss', diag%axesT1, &
                Time, 'Internal tide energy loss to bottom drag', &
                'W m-2', conversion=HZ2_T3_to_W_m2)
  CS%id_tot_itidal_loss = register_diag_field('ocean_model', 'ITide_tot_itidal_loss', diag%axesT1, &
                Time, 'Internal tide energy loss to wave drag', &
                'W m-2', conversion=HZ2_T3_to_W_m2)
  CS%id_tot_Froude_loss = register_diag_field('ocean_model', 'ITide_tot_Froude_loss', diag%axesT1, &
                Time, 'Internal tide energy loss to wave breaking', &
                'W m-2', conversion=HZ2_T3_to_W_m2)
  CS%id_tot_residual_loss = register_diag_field('ocean_model', 'ITide_tot_residual_loss', diag%axesT1, &
                Time, 'Internal tide energy loss to residual on slopes', &
                'W m-2', conversion=HZ2_T3_to_W_m2)
  CS%id_tot_allprocesses_loss = register_diag_field('ocean_model', 'ITide_tot_allprocesses_loss', diag%axesT1, &
                Time, 'Internal tide energy loss summed over all processes', &
                'W m-2', conversion=HZ2_T3_to_W_m2)

  allocate(CS%id_En_mode(CS%nFreq,CS%nMode), source=-1)
  allocate(CS%id_En_ang_mode(CS%nFreq,CS%nMode), source=-1)
  allocate(CS%id_itidal_loss_mode(CS%nFreq,CS%nMode), source=-1)
  allocate(CS%id_leak_loss_mode(CS%nFreq,CS%nMode), source=-1)
  allocate(CS%id_quad_loss_mode(CS%nFreq,CS%nMode), source=-1)
  allocate(CS%id_Froude_loss_mode(CS%nFreq,CS%nMode), source=-1)
  allocate(CS%id_residual_loss_mode(CS%nFreq,CS%nMode), source=-1)
  allocate(CS%id_allprocesses_loss_mode(CS%nFreq,CS%nMode), source=-1)
  allocate(CS%id_itidal_loss_ang_mode(CS%nFreq,CS%nMode), source=-1)
  allocate(CS%id_Ub_mode(CS%nFreq,CS%nMode), source=-1)
  allocate(CS%id_Ustruct_mode(CS%nMode), source=-1)
  allocate(CS%id_Wstruct_mode(CS%nMode), source=-1)
  allocate(CS%id_int_w2_mode(CS%nMode), source=-1)
  allocate(CS%id_int_U2_mode(CS%nMode), source=-1)
  allocate(CS%id_int_N2w2_mode(CS%nMode), source=-1)
  allocate(CS%id_cp_mode(CS%nFreq,CS%nMode), source=-1)

  allocate(angles(CS%NAngle), source=0.0)
  Angle_size = (8.0*atan(1.0)) / (real(num_angle))
  do a=1,num_angle ; angles(a) = (real(a) - 1) * Angle_size ; enddo

  id_ang = diag_axis_init("angle", angles, "Radians", "N", "Angular Orientation of Fluxes")
  call define_axes_group(diag, (/ diag%axesT1%handles(1), diag%axesT1%handles(2), id_ang /), &
                         axes_ang, is_h_point=.true.)
  do fr=1,CS%nFreq ; write(freq_name(fr), '("freq",i1)') fr ; enddo
  do m=1,CS%nMode ; do fr=1,CS%nFreq
    ! Register 2-D energy density (summed over angles) for each frequency and mode
    write(var_name, '("Itide_En_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy density in frequency ",i1," mode ",i1)') fr, m
    CS%id_En_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'J m-2', conversion=HZ2_T2_to_J_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 3-D (i,j,a) energy density for each frequency and mode
    write(var_name, '("Itide_En_ang_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide angular energy density in frequency ",i1," mode ",i1)') fr, m
    CS%id_En_ang_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 axes_ang, Time, var_descript, 'J m-2 band-1', conversion=HZ2_T2_to_J_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 2-D energy loss (summed over angles) for each frequency and mode
    ! wave-drag only
    write(var_name, '("Itide_wavedrag_loss_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy loss due to wave-drag from frequency ",i1," mode ",i1)') fr, m
    CS%id_itidal_loss_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'W m-2', conversion=HZ2_T3_to_W_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    ! Leakage loss
    write(var_name, '("Itide_leak_loss_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy loss due to leakage from frequency ",i1," mode ",i1)') fr, m
    CS%id_leak_loss_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'W m-2', conversion=HZ2_T3_to_W_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    ! Quad loss
    write(var_name, '("Itide_quad_loss_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy quad loss from frequency ",i1," mode ",i1)') fr, m
    CS%id_quad_loss_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'W m-2', conversion=HZ2_T3_to_W_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    ! Froude loss
    write(var_name, '("Itide_froude_loss_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy Froude loss from frequency ",i1," mode ",i1)') fr, m
    CS%id_froude_loss_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'W m-2', conversion=HZ2_T3_to_W_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    ! residual losses
    write(var_name, '("Itide_residual_loss_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy residual loss from frequency ",i1," mode ",i1)') fr, m
    CS%id_residual_loss_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'W m-2', conversion=HZ2_T3_to_W_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    ! all loss processes
    write(var_name, '("Itide_allprocesses_loss_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy loss due to all processes from frequency ",i1," mode ",i1)') fr, m
    CS%id_allprocesses_loss_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'W m-2', conversion=HZ2_T3_to_W_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 3-D (i,j,a) energy loss for each frequency and mode
    ! wave-drag only
    write(var_name, '("Itide_wavedrag_loss_ang_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy loss due to wave-drag from frequency ",i1," mode ",i1)') fr, m
    CS%id_itidal_loss_ang_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 axes_ang, Time, var_descript, 'W m-2 band-1', conversion=HZ2_T3_to_W_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 2-D period-averaged near-bottom horizontal velocity for each frequency and mode
    write(var_name, '("Itide_Ub_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Near-bottom horizontal velocity for frequency ",i1," mode ",i1)') fr, m
    CS%id_Ub_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'm s-1', conversion=US%L_T_to_m_s)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 2-D horizontal phase velocity for each frequency and mode
    write(var_name, '("Itide_cp_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Horizontal phase velocity for frequency ",i1," mode ",i1)') fr, m
    CS%id_cp_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'm s-1', conversion=US%L_T_to_m_s)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 2-D drag scale used for quadratic bottom drag for each frequency and mode
    write(var_name, '("ITide_drag_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Interior and bottom drag int tide decay timescale in frequency ",i1, " mode ",i1)') fr, m

    CS%id_itide_drag(fr,m) = register_diag_field('ocean_model', var_name, diag%axesT1, Time, &
                                                 's-1', conversion=US%s_to_T)
  enddo ; enddo


  do m=1,CS%nMode

    ! Register 3-D internal tide horizonal velocity profile for each mode
    write(var_name, '("Itide_Ustruct","_mode",i1)') m
    write(var_descript, '("horizonal velocity profile for mode ",i1)') m
    CS%id_Ustruct_mode(m) = register_diag_field('ocean_model', var_name, &
                 diag%axesTl, Time, var_descript, 'm-1', conversion=US%m_to_L)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 3-D internal tide vertical velocity profile for each mode
    write(var_name, '("Itide_Wstruct","_mode",i1)') m
    write(var_descript, '("vertical velocity profile for mode ",i1)') m
    CS%id_Wstruct_mode(m) = register_diag_field('ocean_model', var_name, &
                 diag%axesTi, Time, var_descript, 'nondim')
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    write(var_name, '("Itide_int_w2","_mode",i1)') m
    write(var_descript, '("integral of w2 for mode ",i1)') m
    CS%id_int_w2_mode(m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'm', conversion=GV%H_to_m)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    write(var_name, '("Itide_int_U2","_mode",i1)') m
    write(var_descript, '("integral of U2 for mode ",i1)') m
    CS%id_int_U2_mode(m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'm-1', conversion=US%m_to_Z*GV%H_to_Z)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    write(var_name, '("Itide_int_N2w2","_mode",i1)') m
    write(var_descript, '("integral of N2w2 for mode ",i1)') m
    CS%id_int_N2w2_mode(m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'm s-2', conversion=GV%H_to_m*US%s_to_T**2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

  enddo

  ! Initialize the module that calculates the wave speeds.
  call wave_speed_init(CS%wave_speed, GV, c1_thresh=IGW_c1_thresh, &
                       om4_remap_via_sub_cells=om4_remap_via_sub_cells)

end subroutine internal_tides_init

!> This subroutine deallocates the memory associated with the internal tides control structure
subroutine internal_tides_end(CS)
  type(int_tide_CS), intent(inout) :: CS  !<  Internal tide control structure

  if (allocated(CS%En)) deallocate(CS%En)
  if (allocated(CS%frequency)) deallocate(CS%frequency)
  if (allocated(CS%id_En_mode)) deallocate(CS%id_En_mode)
  if (allocated(CS%id_Ub_mode)) deallocate(CS%id_Ub_mode)
  if (allocated(CS%id_cp_mode)) deallocate(CS%id_cp_mode)
  if (allocated(CS%id_Ustruct_mode)) deallocate(CS%id_Ustruct_mode)
  if (allocated(CS%id_Wstruct_mode)) deallocate(CS%id_Wstruct_mode)
  if (allocated(CS%id_int_w2_mode)) deallocate(CS%id_int_w2_mode)
  if (allocated(CS%id_int_U2_mode)) deallocate(CS%id_int_U2_mode)
  if (allocated(CS%id_int_N2w2_mode)) deallocate(CS%id_int_N2w2_mode)

end subroutine internal_tides_end

end module MOM_internal_tides
