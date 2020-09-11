!> Implements the thermodynamic aspects of ocean / ice-shelf interactions,
!!  along with a crude placeholder for a later implementation of full
!!  ice shelf dynamics, all using the MOM framework and coding style.
module MOM_ice_shelf

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_constants, only : hlf
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_mediator_init, set_diag_mediator_grid, diag_ctrl, time_type
use MOM_diag_mediator, only : enable_averages, enable_averaging, disable_averaging
use MOM_domains, only : MOM_domains_init, clone_MOM_domain
use MOM_domains, only : pass_var, pass_vector, TO_ALL, CGRID_NE, BGRID_NE, CORNER
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_dyn_horgrid, only : rescale_dyn_horgrid_bathymetry
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : MOM_grid_init, ocean_grid_type
use MOM_grid_initialize, only : set_grid_metrics
use MOM_fixed_initialization, only : MOM_initialize_topography
use MOM_fixed_initialization, only : MOM_initialize_rotation
use user_initialization, only : user_initialize_topography
use MOM_io, only : field_exists, file_exists, MOM_read_data, write_version_number
use MOM_io, only : slasher, fieldtype
use MOM_io, only : write_field, close_file, SINGLE_FILE, MULTIPLE
use MOM_restart, only : register_restart_field, query_initialized, save_restart
use MOM_restart, only : restart_init, restore_state, MOM_restart_CS
use MOM_time_manager, only : time_type, time_type_to_real, real_to_time, operator(>), operator(-)
use MOM_transcribe_grid, only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
use MOM_unit_scaling, only : unit_scale_type, unit_scaling_init, fix_restart_unit_scaling
use MOM_variables, only : surface
use MOM_forcing_type, only : forcing, allocate_forcing_type, MOM_forcing_chksum
use MOM_forcing_type, only : mech_forcing, allocate_mech_forcing, MOM_mech_forcing_chksum
use MOM_forcing_type, only : copy_common_forcing_fields
use MOM_get_input, only : directories, Get_MOM_input
use MOM_EOS, only : calculate_density, calculate_density_derivs, calculate_TFreeze, EOS_domain
use MOM_EOS, only : EOS_type, EOS_init
use MOM_ice_shelf_dynamics, only : ice_shelf_dyn_CS, update_ice_shelf
use MOM_ice_shelf_dynamics, only : register_ice_shelf_dyn_restarts, initialize_ice_shelf_dyn
use MOM_ice_shelf_dynamics, only : ice_shelf_min_thickness_calve
use MOM_ice_shelf_dynamics, only : ice_time_step_CFL, ice_shelf_dyn_end
use MOM_ice_shelf_initialize, only : initialize_ice_thickness
!MJH use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary
use MOM_ice_shelf_state, only : ice_shelf_state, ice_shelf_state_end, ice_shelf_state_init
use user_shelf_init, only : USER_initialize_shelf_mass, USER_update_shelf_mass
use user_shelf_init, only : user_ice_shelf_CS
use MOM_coms, only : reproducing_sum
use MOM_spatial_means, only : global_area_integral
use MOM_checksums, only : hchksum, qchksum, chksum, uchksum, vchksum, uvchksum
use time_interp_external_mod, only : init_external_field, time_interp_external
use time_interp_external_mod, only : time_interp_external_init
implicit none ; private

#include <MOM_memory.h>
#ifdef SYMMETRIC_LAND_ICE
#  define GRID_SYM_ .true.
#else
#  define GRID_SYM_ .false.
#endif

public shelf_calc_flux, add_shelf_flux, initialize_ice_shelf, ice_shelf_end
public ice_shelf_save_restart, solo_step_ice_shelf, add_shelf_forces

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure that contains ice shelf parameters and diagnostics handles
type, public :: ice_shelf_CS ; private
  ! Parameters
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control
                                          !! structure for the ice shelves
  type(ocean_grid_type) :: grid           !< Grid for the ice-shelf model
  type(unit_scale_type), pointer :: &
    US => NULL()       !< A structure containing various unit conversion factors
  !type(dyn_horgrid_type), pointer :: dG  !< Dynamic grid for the ice-shelf model
  type(ocean_grid_type), pointer :: ocn_grid => NULL() !< A pointer to the ocean model grid
                                          !! The rest is private
  real ::   flux_factor = 1.0             !< A factor that can be used to turn off ice shelf
                                          !! melting (flux_factor = 0) [nondim].
  character(len=128) :: restart_output_dir = ' ' !< The directory in which to write restart files
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
                                          !! the ice-shelf state
  type(ice_shelf_dyn_CS), pointer :: dCS => NULL() !< The control structure for the ice-shelf dynamics.

  real, pointer, dimension(:,:) :: &
    utide   => NULL()  !< An unresolved tidal velocity [L T-1 ~> m s-1]

  real :: ustar_bg     !< A minimum value for ustar under ice shelves [Z T-1 ~> m s-1].
  real :: cdrag        !< drag coefficient under ice shelves [nondim].
  real :: g_Earth      !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: Cp           !< The heat capacity of sea water [Q degC-1 ~> J kg-1 degC-1].
  real :: Rho_ocn      !< A reference ocean density [R ~> kg m-3].
  real :: Cp_ice       !< The heat capacity of fresh ice [Q degC-1 ~> J kg-1 degC-1].
  real :: gamma_t      !< The (fixed) turbulent exchange velocity in the
                       !< 2-equation formulation [Z T-1 ~> m s-1].
  real :: Salin_ice    !< The salinity of shelf ice [ppt].
  real :: Temp_ice     !< The core temperature of shelf ice [degC].
  real :: kv_ice       !< The viscosity of ice [L4 Z-2 T-1 ~> m2 s-1].
  real :: density_ice  !< A typical density of ice [R ~> kg m-3].
  real :: kv_molec     !< The molecular kinematic viscosity of sea water [Z2 T-1 ~> m2 s-1].
  real :: kd_molec_salt!< The molecular diffusivity of salt [Z2 T-1 ~> m2 s-1].
  real :: kd_molec_temp!< The molecular diffusivity of heat [Z2 T-1 ~> m2 s-1].
  real :: Lat_fusion   !< The latent heat of fusion [Q ~> J kg-1].
  real :: Gamma_T_3EQ  !<  Nondimensional heat-transfer coefficient, used in the 3Eq. formulation
  real :: Gamma_S_3EQ  !<  Nondimensional salt-transfer coefficient, used in the 3Eq. formulation
                       !<  This number should be specified by the user.
  real :: col_mass_melt_threshold !< An ocean column mass below the iceshelf below which melting
                       !! does not occur [R Z ~> kg m-2]
  logical :: mass_from_file !< Read the ice shelf mass from a file every dt

  !!!! PHYSICAL AND NUMERICAL PARAMETERS FOR ICE DYNAMICS !!!!!!

  real :: time_step    !< this is the shortest timestep that the ice shelf sees, and
                       !! is equal to the forcing timestep (it is passed in when the shelf
                       !! is initialized - so need to reorganize MOM driver.
                       !! it will be the prognistic timestep ... maybe.

  logical :: solo_ice_sheet !< whether the ice model is running without being
                            !! coupled to the ocean
  logical :: GL_regularize  !< whether to regularize the floatation condition
                            !! at the grounding line a la Goldberg Holland Schoof 2009
  logical :: GL_couple      !< whether to let the floatation condition be
                            !!determined by ocean column thickness means update_OD_ffrac
                            !! will be called (note: GL_regularize and GL_couple
                            !! should be exclusive)
  logical :: calve_to_mask  !< If true, calve any ice that passes outside of a masked area
  real :: min_thickness_simple_calve !< min. ice shelf thickness criteria for calving [Z ~> m].
  real :: T0                !< temperature at ocean surface in the restoring region [degC]
  real :: S0                !< Salinity at ocean surface in the restoring region [ppt].
  real :: input_flux        !< Ice volume flux at an upstream open boundary [m3 s-1].
  real :: input_thickness   !< Ice thickness at an upstream open boundary [m].

  type(time_type) :: Time                !< The component's time.
  type(EOS_type), pointer :: eqn_of_state => NULL() !< Type that indicates the
                                         !! equation of state to use.
  logical :: active_shelf_dynamics       !< True if the ice shelf mass changes as a result
                                         !! the dynamic ice-shelf model.
  logical :: override_shelf_movement     !< If true, user code specifies the shelf movement
                                         !! instead of using the dynamic ice-shelf mode.
  logical :: isthermo                    !< True if the ice shelf can exchange heat and
                                         !! mass with the underlying ocean.
  logical :: threeeq                     !< If true, the 3 equation consistency equations are
                                         !! used to calculate the flux at the ocean-ice
                                         !! interface.
  logical :: insulator                   !< If true, ice shelf is a perfect insulator
  logical :: const_gamma                 !< If true, gamma_T is specified by the user.
  logical :: constant_sea_level          !< if true, apply an evaporative, heat and salt
                                         !! fluxes. It will avoid large increase in sea level.
  real    :: min_ocean_mass_float        !< The minimum ocean mass per unit area before the ice
                                         !! shelf is considered to float when constant_sea_level
                                         !! is used [R Z ~> kg m-2]
  real    :: cutoff_depth                !< Depth above which melt is set to zero (>= 0) [Z ~> m].
  logical :: find_salt_root              !< If true, if true find Sbdry using a quadratic eq.
  real    :: TFr_0_0                     !< The freezing point at 0 pressure and 0 salinity [degC]
  real    :: dTFr_dS                     !< Partial derivative of freezing temperature with salinity [degC ppt-1]
  real    :: dTFr_dp                     !< Partial derivative of freezing temperature with
                                         !! pressure [degC T2 R-1 L-2 ~> degC Pa-1]
  !>@{ Diagnostic handles
  integer :: id_melt = -1, id_exch_vel_s = -1, id_exch_vel_t = -1, &
             id_tfreeze = -1, id_tfl_shelf = -1, &
             id_thermal_driving = -1, id_haline_driving = -1, &
             id_u_ml = -1, id_v_ml = -1, id_sbdry = -1, &
             id_h_shelf = -1, id_h_mask = -1, &
             id_surf_elev = -1, id_bathym = -1, &
             id_area_shelf_h = -1, &
             id_ustar_shelf = -1, id_shelf_mass = -1, id_mass_flux = -1
  !>@}

  integer :: id_read_mass !< An integer handle used in time interpolation of
                          !! the ice shelf mass read from a file
  integer :: id_read_area !< An integer handle used in time interpolation of
                          !! the ice shelf mass read from a file

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to control diagnostic output.
  type(user_ice_shelf_CS), pointer :: user_CS => NULL() !< A pointer to the control structure for
                                  !! user-supplied modifications to the ice shelf code.

  logical :: debug                !< If true, write verbose checksums for debugging purposes
                                  !! and use reproducible sums
end type ice_shelf_CS

integer :: id_clock_shelf !< CPU Clock for the ice shelf code
integer :: id_clock_pass !< CPU Clock for group pass calls

contains

!> Calculates fluxes between the ocean and ice-shelf using the three-equations
!! formulation (optional to use just two equations).
!! See \ref section_ICE_SHELF_equations
subroutine shelf_calc_flux(sfc_state, fluxes, Time, time_step, CS, forces)
  type(surface),         intent(inout) :: sfc_state !< A structure containing fields that
                                                !! describe the surface state of the ocean.  The
                                                !! intent is only inout to allow for halo updates.
  type(forcing),         intent(inout) :: fluxes !< structure containing pointers to any possible
                                                !! thermodynamic or mass-flux forcing fields.
  type(time_type),       intent(in)    :: Time  !< Start time of the fluxes.
  real,                  intent(in)    :: time_step !< Length of time over which these fluxes
                                                !! will be applied [s].
  type(ice_shelf_CS),    pointer       :: CS    !< A pointer to the control structure returned
                                                !! by a previous call to initialize_ice_shelf.
  type(mech_forcing), optional, intent(inout) :: forces !< A structure with the driving mechanical forces

  ! Local variables
  type(ocean_grid_type), pointer :: G => NULL()  !< The grid structure used by the ice shelf.
  type(unit_scale_type), pointer :: US => NULL() !< Pointer to a structure containing
                                                 !! various unit conversion factors
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
                                                 !! the ice-shelf state

  real, dimension(SZI_(CS%grid)) :: &
    Rhoml, &   !< Ocean mixed layer density [R ~> kg m-3].
    dR0_dT, &  !< Partial derivative of the mixed layer density
               !< with temperature [R degC-1 ~> kg m-3 degC-1].
    dR0_dS, &  !< Partial derivative of the mixed layer density
               !< with salinity [R ppt-1 ~> kg m-3 ppt-1].
    p_int      !< The pressure at the ice-ocean interface [R L2 T-2 ~> Pa].

  real, dimension(SZI_(CS%grid),SZJ_(CS%grid)) :: &
    exch_vel_t, &  !< Sub-shelf thermal exchange velocity [Z T-1 ~> m s-1]
    exch_vel_s     !< Sub-shelf salt exchange velocity [Z T-1 ~> m s-1]

  real, dimension(SZDI_(CS%grid),SZDJ_(CS%grid)) :: &
    mass_flux  !< Total mass flux of freshwater across the ice-ocean interface. [R Z L2 T-1 ~> kg/s]
  real, dimension(SZDI_(CS%grid),SZDJ_(CS%grid)) :: &
    haline_driving !< (SSS - S_boundary) ice-ocean
               !! interface, positive for melting and negative for freezing.
               !! This is computed as part of the ISOMIP diagnostics.
  real, parameter :: VK    = 0.40 !< Von Karman's constant - dimensionless
  real :: ZETA_N = 0.052 !> The fraction of the boundary layer over which the
               !! viscosity is linearly increasing [nondim]. (Was 1/8. Why?)
  real, parameter :: RC    = 0.20     ! critical flux Richardson number.
  real :: I_ZETA_N !< The inverse of ZETA_N [nondim].
  real :: I_LF     !< The inverse of the latent heat of fusion [Q-1 ~> kg J-1].
  real :: I_VK     !< The inverse of the Von Karman constant [nondim].
  real :: PR, SC   !< The Prandtl number and Schmidt number [nondim].

  ! 3 equations formulation variables
  real, dimension(SZDI_(CS%grid),SZDJ_(CS%grid)) :: &
    Sbdry     !< Salinities in the ocean at the interface with the ice shelf [ppt].
  real :: Sbdry_it
  real :: Sbdry1, Sbdry2
  real :: S_a, S_b, S_c  ! Variables used to find salt roots
  real :: dS_it    !< The interface salinity change during an iteration [ppt].
  real :: hBL_neut !< The neutral boundary layer thickness [Z ~> m].
  real :: hBL_neut_h_molec !< The ratio of the neutral boundary layer thickness
                   !! to the molecular boundary layer thickness [nondim].
  real :: wT_flux !< The downward vertical flux of heat just inside the ocean [degC Z T-1 ~> degC m s-1].
  real :: wB_flux !< The downward vertical flux of buoyancy just inside the ocean [Z2 T-3 ~> m2 s-3].
  real :: dB_dS   !< The derivative of buoyancy with salinity [Z T-2 ppt-1 ~> m s-2 ppt-1].
  real :: dB_dT   !< The derivative of buoyancy with temperature [Z T-2 degC-1 ~> m s-2 degC-1].
  real :: I_n_star ! [nondim]
  real :: n_star_term ! A term in the expression for nstar [T3 Z-2 ~> s3 m-2]
  real :: absf     ! The absolute value of the Coriolis parameter [T-1 ~> s-1]
  real :: dIns_dwB !< The partial derivative of I_n_star with wB_flux, in [T3 Z-2 ~> s3 m-2]
  real :: dT_ustar ! The difference between the the freezing point and the ocean boundary layer
                   ! temperature times the friction velocity [degC Z T-1 ~> degC m s-1]
  real :: dS_ustar ! The difference between the salinity at the ice-ocean interface and the ocean
                   ! boundary layer salinity times the friction velocity [ppt Z T-1 ~> ppt m s-1]
  real :: ustar_h  ! The friction velocity in the water below the ice shelf [Z T-1 ~> m s-1]
  real :: Gam_turb ! [nondim]
  real :: Gam_mol_t, Gam_mol_s ! Relative coefficients of molecular diffusivites [nondim]
  real :: RhoCp     ! A typical ocean density times the heat capacity of water [Q R ~> J m-3]
  real :: ln_neut
  real :: mass_exch ! A mass exchange rate [R Z T-1 ~> kg m-2 s-1]
  real :: Sb_min, Sb_max
  real :: dS_min, dS_max
  ! Variables used in iterating for wB_flux.
  real :: wB_flux_new, dDwB_dwB_in
  real :: I_Gam_T, I_Gam_S
  real :: dG_dwB   ! The derivative of Gam_turb with wB [T3 Z-2 ~> s3 m-2]
  real :: taux2, tauy2 ! The squared surface stresses [R2 L2 Z2 T-4 ~> Pa2].
  real :: u2_av, v2_av ! The ice-area weighted average squared ocean velocities [L2 T-2 ~> m2 s-2]
  real :: asu1, asu2   ! Ocean areas covered by ice shelves at neighboring u-
  real :: asv1, asv2   ! and v-points [L2 ~> m2].
  real :: I_au, I_av   ! The Adcroft reciprocals of the ice shelf areas at adjacent points [L-2 ~> m-2]
  real :: Irho0        ! The inverse of the mean density times a unit conversion factor [R-1 L Z-1 ~> m3 kg-1]
  logical :: Sb_min_set, Sb_max_set
  logical :: update_ice_vel ! If true, it is time to update the ice shelf velocities.
  logical :: coupled_GL     ! If true, the grouding line position is determined based on
                            ! coupled ice-ocean dynamics.

  real, parameter :: c2_3 = 2.0/3.0
  character(len=160) :: mesg  ! The text of an error message
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, is, ie, js, je, ied, jed, it1, it3

  if (.not. associated(CS)) call MOM_error(FATAL, "shelf_calc_flux: "// &
       "initialize_ice_shelf must be called before shelf_calc_flux.")
  call cpu_clock_begin(id_clock_shelf)

  G => CS%grid ; US => CS%US
  ISS => CS%ISS

  ! useful parameters
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; ied = G%ied ; jed = G%jed
  I_ZETA_N = 1.0 / ZETA_N
  I_LF = 1.0 / CS%Lat_fusion
  SC = CS%kv_molec/CS%kd_molec_salt
  PR = CS%kv_molec/CS%kd_molec_temp
  I_VK = 1.0/VK
  RhoCp = CS%Rho_ocn * CS%Cp

  !first calculate molecular component
  Gam_mol_t = 12.5 * (PR**c2_3) - 6.0
  Gam_mol_s = 12.5 * (SC**c2_3) - 6.0

  ! GMM, zero some fields of the ice shelf structure (ice_shelf_CS)
  ! these fields are already set to zero during initialization
  ! However, they seem to be changed somewhere and, for diagnostic
  ! reasons, it is better to set them to zero again.
  exch_vel_t(:,:) = 0.0 ; exch_vel_s(:,:) = 0.0
  ISS%tflux_shelf(:,:) = 0.0 ; ISS%water_flux(:,:) = 0.0
  ISS%salt_flux(:,:) = 0.0 ; ISS%tflux_ocn(:,:) = 0.0 ; ISS%tfreeze(:,:) = 0.0
  ! define Sbdry to avoid Run-Time Check Failure, when melt is not computed.
  haline_driving(:,:) = 0.0
  Sbdry(:,:) = sfc_state%sss(:,:)

  !update time
  CS%Time = Time

  if (CS%override_shelf_movement) then
    CS%time_step = time_step
    ! update shelf mass
    if (CS%mass_from_file) call update_shelf_mass(G, US, CS, ISS, Time)
  endif

  if (CS%debug) then
    call hchksum(fluxes%frac_shelf_h, "frac_shelf_h before apply melting", G%HI, haloshift=0)
    call hchksum(sfc_state%sst, "sst before apply melting", G%HI, haloshift=0)
    call hchksum(sfc_state%sss, "sss before apply melting", G%HI, haloshift=0)
    call hchksum(sfc_state%u, "u_ml before apply melting", G%HI, haloshift=0, scale=US%L_T_to_m_s)
    call hchksum(sfc_state%v, "v_ml before apply melting", G%HI, haloshift=0, scale=US%L_T_to_m_s)
    call hchksum(sfc_state%ocean_mass, "ocean_mass before apply melting", G%HI, haloshift=0, &
                 scale=US%RZ_to_kg_m2)
  endif

  ! Calculate the friction velocity under ice shelves, using taux_shelf and tauy_shelf if possible.
  if (allocated(sfc_state%taux_shelf) .and. allocated(sfc_state%tauy_shelf)) then
    call pass_vector(sfc_state%taux_shelf, sfc_state%tauy_shelf, G%domain, TO_ALL, CGRID_NE)
  endif
  Irho0 = US%Z_to_L / CS%Rho_ocn
  do j=js,je ; do i=is,ie ; if (fluxes%frac_shelf_h(i,j) > 0.0) then
    taux2 = 0.0 ; tauy2 = 0.0 ; u2_av = 0.0 ; v2_av = 0.0
    asu1 = (ISS%area_shelf_h(i-1,j) + ISS%area_shelf_h(i,j))
    asu2 = (ISS%area_shelf_h(i,j) + ISS%area_shelf_h(i+1,j))
    asv1 = (ISS%area_shelf_h(i,j-1) + ISS%area_shelf_h(i,j))
    asv2 = (ISS%area_shelf_h(i,j) + ISS%area_shelf_h(i,j+1))
    I_au = 0.0 ; if (asu1 + asu2 > 0.0) I_au = 1.0 / (asu1 + asu2)
    I_av = 0.0 ; if (asv1 + asv2 > 0.0) I_av = 1.0 / (asv1 + asv2)
    if (allocated(sfc_state%taux_shelf) .and. allocated(sfc_state%tauy_shelf)) then
      taux2 = (asu1 * sfc_state%taux_shelf(I-1,j)**2 + asu2 * sfc_state%taux_shelf(I,j)**2  ) * I_au
      tauy2 = (asv1 * sfc_state%tauy_shelf(i,J-1)**2 + asv2 * sfc_state%tauy_shelf(i,J)**2  ) * I_av
    endif
    u2_av = (asu1 * sfc_state%u(I-1,j)**2 + asu2 * sfc_state%u(I,j)**2) * I_au
    v2_av = (asv1 * sfc_state%v(i,J-1)**2 + asu2 * sfc_state%v(i,J)**2) * I_av

    if (taux2 + tauy2 > 0.0) then
      fluxes%ustar_shelf(i,j) = MAX(CS%ustar_bg, US%L_to_Z * &
          sqrt(Irho0 * sqrt(taux2 + tauy2) + CS%cdrag*CS%utide(i,j)**2))
    else   ! Take care of the cases when taux_shelf is not set or not allocated.
      fluxes%ustar_shelf(i,j) = MAX(CS%ustar_bg, US%L_TO_Z * &
          sqrt(CS%cdrag*((u2_av + v2_av) + CS%utide(i,j)**2)))
    endif
  else ! There is no shelf here.
    fluxes%ustar_shelf(i,j) = 0.0
  endif ; enddo ; enddo

  EOSdom(:) = EOS_domain(G%HI)
  do j=js,je
    ! Find the pressure at the ice-ocean interface, averaged only over the
    ! part of the cell covered by ice shelf.
    do i=is,ie ; p_int(i) = CS%g_Earth * ISS%mass_shelf(i,j) ; enddo

    ! Calculate insitu densities and expansion coefficients
    call calculate_density(sfc_state%sst(:,j), sfc_state%sss(:,j), p_int, Rhoml(:), &
                                 CS%eqn_of_state, EOSdom)
    call calculate_density_derivs(sfc_state%sst(:,j), sfc_state%sss(:,j), p_int, dR0_dT, dR0_dS, &
                                 CS%eqn_of_state, EOSdom)

    do i=is,ie
      if ((sfc_state%ocean_mass(i,j) > CS%col_mass_melt_threshold) .and. &
          (ISS%area_shelf_h(i,j) > 0.0) .and. CS%isthermo) then

        if (CS%threeeq) then
          !   Iteratively determine a self-consistent set of fluxes, with the ocean
          ! salinity just below the ice-shelf as the variable that is being
          ! iterated for.

          ustar_h = fluxes%ustar_shelf(i,j)

          ! Estimate the neutral ocean boundary layer thickness as the minimum of the
          ! reported ocean mixed layer thickness and the neutral Ekman depth.
          absf = 0.25*((abs(G%CoriolisBu(I,J)) + abs(G%CoriolisBu(I-1,J-1))) + &
                                 (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I-1,J))))
          if (absf*sfc_state%Hml(i,j) <= VK*ustar_h) then ; hBL_neut = sfc_state%Hml(i,j)
          else ; hBL_neut = (VK*ustar_h) / absf ; endif
          hBL_neut_h_molec = ZETA_N * ((hBL_neut * ustar_h) / (5.0 * CS%kv_molec))

          ! Determine the mixed layer buoyancy flux, wB_flux.
          dB_dS = (US%L_to_Z**2*CS%g_Earth / Rhoml(i)) * dR0_dS(i)
          dB_dT = (US%L_to_Z**2*CS%g_Earth / Rhoml(i)) * dR0_dT(i)
          ln_neut = 0.0 ; if (hBL_neut_h_molec > 1.0) ln_neut = log(hBL_neut_h_molec)

          if (CS%find_salt_root) then
            ! Solve for the skin salinity using the linearized liquidus parameters and
            ! balancing the turbulent fresh water flux in the near-boundary layer with
            ! the net fresh water or salt added by melting:
            ! (Cp/Lat_fusion)*Gamma_T_3Eq*(TFr_skin-T_ocn) = Gamma_S_3Eq*(S_skin-S_ocn)/S_skin

            ! S_a is always < 0.0 with a realistic expression for the freezing point.
            S_a = CS%dTFr_dS * CS%Gamma_T_3EQ * CS%Cp
            S_b = CS%Gamma_T_3EQ*CS%Cp*(CS%TFr_0_0 + CS%dTFr_dp*p_int(i) - sfc_state%sst(i,j)) - &
                  CS%Lat_fusion * CS%Gamma_S_3EQ    ! S_b Can take either sign, but is usually negative.
            S_c = CS%Lat_fusion * CS%Gamma_S_3EQ * sfc_state%sss(i,j) ! Always >= 0

            if (S_c == 0.0) then  ! The solution for fresh water.
              Sbdry(i,j) = 0.0
            elseif (S_a < 0.0) then ! This is the usual ocean case
              if (S_b < 0.0) then ! This is almost always the case
                Sbdry(i,j) = 2.0*S_c / (-S_b + SQRT(S_b*S_b - 4.*S_a*S_c))
              else
                Sbdry(i,j) = (S_b + SQRT(S_b*S_b - 4.*S_a*S_c)) / (-2.*S_a)
              endif
            elseif ((S_a == 0.0) .and. (S_b < 0.0)) then ! It should be the case that S_b < 0.
              Sbdry(i,j) = -S_c / S_b
            else
              call MOM_error(FATAL, "Impossible conditions found in 3-equation skin salinity calculation.")
            endif

            ! Safety check
            if (Sbdry(i,j) < 0.) then
              write(mesg,*) 'sfc_state%sss(i,j) = ',sfc_state%sss(i,j), 'S_a, S_b, S_c', S_a, S_b, S_c
              call MOM_error(WARNING, mesg, .true.)
              write(mesg,*) 'I,J,Sbdry1,Sbdry2',i,j,Sbdry1,Sbdry2
              call MOM_error(WARNING, mesg, .true.)
              call MOM_error(FATAL, "shelf_calc_flux: Negative salinity (Sbdry).")
            endif
          else
            ! Guess sss as the iteration starting point for the boundary salinity.
            Sbdry(i,j) = sfc_state%sss(i,j) ; Sb_max_set = .false.
            Sb_min_set = .false.
          endif !find_salt_root

          do it1 = 1,20
            ! Determine the potential temperature at the ice-ocean interface.
            call calculate_TFreeze(Sbdry(i,j), p_int(i), ISS%tfreeze(i,j), CS%eqn_of_state, &
                                   pres_scale=US%RL2_T2_to_Pa)

            dT_ustar = (ISS%tfreeze(i,j) - sfc_state%sst(i,j)) * ustar_h
            dS_ustar = (Sbdry(i,j) - sfc_state%sss(i,j)) * ustar_h

            ! First, determine the buoyancy flux assuming no effects of stability
            ! on the turbulence.  Following H & J '99, this limit also applies
            ! when the buoyancy flux is destabilizing.

            if (CS%const_gamma) then ! if using a constant gamma_T
              ! note the different form, here I_Gam_T is NOT 1/Gam_T!
              I_Gam_T = CS%Gamma_T_3EQ
              I_Gam_S = CS%Gamma_S_3EQ
            else
              Gam_turb = I_VK * (ln_neut + (0.5 * I_ZETA_N - 1.0))
              I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
              I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
            endif

            wT_flux = dT_ustar * I_Gam_T
            wB_flux = dB_dS * (dS_ustar * I_Gam_S) + dB_dT * wT_flux

            if (wB_flux < 0.0) then
              ! The buoyancy flux is stabilizing and will reduce the tubulent
              ! fluxes, and iteration is required.
              n_star_term = (ZETA_N/RC) * (hBL_neut * VK) / (ustar_h)**3
              do it3 = 1,30
               ! n_star <= 1.0 is the ratio of working boundary layer thickness
               ! to the neutral thickness.
               ! hBL = n_star*hBL_neut ; hSub = 1/8*n_star*hBL

                I_n_star = sqrt(1.0 - n_star_term * wB_flux)
                dIns_dwB = 0.5 * n_star_term / I_n_star
                if (hBL_neut_h_molec > I_n_star**2) then
                  Gam_turb = I_VK * ((ln_neut - 2.0*log(I_n_star)) + &
                                    (0.5*I_ZETA_N*I_n_star - 1.0))
                  dG_dwB =  I_VK * ( -2.0 / I_n_star + (0.5 * I_ZETA_N)) * dIns_dwB
                else
                  !   The layer dominated by molecular viscosity is smaller than
                  ! the assumed boundary layer.  This should be rare!
                  Gam_turb = I_VK * (0.5 * I_ZETA_N*I_n_star - 1.0)
                  dG_dwB = I_VK * (0.5 * I_ZETA_N) * dIns_dwB
                endif

                if (CS%const_gamma) then ! if using a constant gamma_T
                  ! note the different form, here I_Gam_T is NOT 1/Gam_T!
                  I_Gam_T = CS%Gamma_T_3EQ
                  I_Gam_S = CS%Gamma_S_3EQ
                else
                  I_Gam_T = 1.0 / (Gam_mol_t + Gam_turb)
                  I_Gam_S = 1.0 / (Gam_mol_s + Gam_turb)
                endif

                wT_flux = dT_ustar * I_Gam_T
                wB_flux_new = dB_dS * (dS_ustar * I_Gam_S) + dB_dT * wT_flux

                ! Find the root where wB_flux_new = wB_flux.  Make the 1.0e-4 below into a parameter?
                if (abs(wB_flux_new - wB_flux) < 1.0e-4*(abs(wB_flux_new) + abs(wB_flux))) exit

                dDwB_dwB_in = dG_dwB * (dB_dS * (dS_ustar * I_Gam_S**2) + &
                                        dB_dT * (dT_ustar * I_Gam_T**2)) - 1.0
                ! This is Newton's method without any bounds.  Should bounds be needed?
                wB_flux_new = wB_flux - (wB_flux_new - wB_flux) / dDwB_dwB_in
              enddo !it3
            endif

            ISS%tflux_ocn(i,j)  = RhoCp * wT_flux
            exch_vel_t(i,j) = ustar_h * I_Gam_T
            exch_vel_s(i,j) = ustar_h * I_Gam_S

            ! Calculate the heat flux inside the ice shelf.
            ! Vertical adv/diff as in H+J 1999, eqns (26) & approx from (31).
            !   Q_ice = density_ice * CS%Cp_ice * K_ice * dT/dz (at interface)
            ! vertical adv/diff as in H+J 1999, eqs (31) & (26)...
            !   dT/dz ~= min( (lprec/(density_ice*K_ice))*(CS%Temp_Ice-T_freeze) , 0.0 )
            ! If this approximation is not made, iterations are required... See H+J Fig 3.

            if (ISS%tflux_ocn(i,j) >= 0.0) then
              ! Freezing occurs due to downward ocean heat flux, so zero iout ce heat flux.
              ISS%water_flux(i,j) = -I_LF * ISS%tflux_ocn(i,j)
              ISS%tflux_shelf(i,j) = 0.0
            else
              if (CS%insulator) then
                !no conduction/perfect insulator
                ISS%tflux_shelf(i,j) = 0.0
                ISS%water_flux(i,j) = I_LF * (ISS%tflux_shelf(i,j) - ISS%tflux_ocn(i,j))

              else
                ! With melting, from H&J 1999, eqs (31) & (26)...
                !   Q_ice ~= Cp_ice * (CS%Temp_Ice-T_freeze) * lprec
                !   RhoLF*lprec = Q_ice - ISS%tflux_ocn(i,j)
                !   lprec = -(ISS%tflux_ocn(i,j)) / (CS%Lat_fusion + Cp_ice * (T_freeze-CS%Temp_Ice))
                ISS%water_flux(i,j) = -ISS%tflux_ocn(i,j) / &
                     (CS%Lat_fusion + CS%Cp_ice * (ISS%tfreeze(i,j) - CS%Temp_Ice))

                ISS%tflux_shelf(i,j) = ISS%tflux_ocn(i,j) + CS%Lat_fusion*ISS%water_flux(i,j)
              endif

            endif
            !other options: dTi/dz linear through shelf, with draft in [Z ~> m], KTI in [Z2 T-1 ~> m2 s-1]
            !    dTi_dz = (CS%Temp_Ice - ISS%tfreeze(i,j)) / draft(i,j)
            !    ISS%tflux_shelf(i,j) = Rho_Ice * CS%Cp_ice * KTI * dTi_dz


            if (CS%find_salt_root) then
              exit ! no need to do interaction, so exit loop
            else

              mass_exch = exch_vel_s(i,j) * CS%Rho_ocn
              Sbdry_it = (sfc_state%sss(i,j) * mass_exch + CS%Salin_ice * ISS%water_flux(i,j)) / &
                         (mass_exch + ISS%water_flux(i,j))
              dS_it = Sbdry_it - Sbdry(i,j)
              if (abs(dS_it) < 1.0e-4*(0.5*(sfc_state%sss(i,j) + Sbdry(i,j) + 1.0e-10))) exit


              if (dS_it < 0.0) then ! Sbdry is now the upper bound.
                if (Sb_max_set .and. (Sbdry(i,j) > Sb_max)) &
                  call MOM_error(FATAL,"shelf_calc_flux: Irregular iteration for Sbdry (max).")
                Sb_max = Sbdry(i,j) ; dS_max = dS_it ; Sb_max_set = .true.
              else ! Sbdry is now the lower bound.
                if (Sb_min_set .and. (Sbdry(i,j) < Sb_min)) &
                   call MOM_error(FATAL, "shelf_calc_flux: Irregular iteration for Sbdry (min).")
                 Sb_min = Sbdry(i,j) ; dS_min = dS_it ; Sb_min_set = .true.
              endif ! dS_it < 0.0

              if (Sb_min_set .and. Sb_max_set) then
                ! Use the false position method for the next iteration.
                Sbdry(i,j) = Sb_min + (Sb_max-Sb_min) * (dS_min / (dS_min - dS_max))
              else
                Sbdry(i,j) = Sbdry_it
              endif ! Sb_min_set

              Sbdry(i,j) = Sbdry_it
            endif ! CS%find_salt_root

          enddo !it1
          ! Check for non-convergence and/or non-boundedness?

        else
          !   In the 2-equation form, the mixed layer turbulent exchange velocity
          ! is specified and large enough that the ocean salinity at the interface
          ! is about the same as the boundary layer salinity.

          call calculate_TFreeze(sfc_state%sss(i,j), p_int(i), ISS%tfreeze(i,j), CS%eqn_of_state, &
                                 pres_scale=US%RL2_T2_to_Pa)

          exch_vel_t(i,j) = CS%gamma_t
          ISS%tflux_ocn(i,j) = RhoCp * exch_vel_t(i,j) * (ISS%tfreeze(i,j) - sfc_state%sst(i,j))
          ISS%tflux_shelf(i,j) = 0.0
          ISS%water_flux(i,j) = -I_LF * ISS%tflux_ocn(i,j)
          Sbdry(i,j) = 0.0
        endif
      elseif (ISS%area_shelf_h(i,j) > 0.0) then ! This is an ice-sheet, not a floating shelf.
        ISS%tflux_ocn(i,j) = 0.0
      else ! There is no ice shelf or sheet here.
        ISS%tflux_ocn(i,j) = 0.0
      endif

!      haline_driving(i,j) = sfc_state%sss(i,j) - Sbdry(i,j)

    enddo ! i-loop
  enddo ! j-loop


  do j=js,je ; do i=is,ie
    ! ISS%water_flux = net liquid water into the ocean [R Z T-1 ~> kg m-2 s-1]
    fluxes%iceshelf_melt(i,j) = ISS%water_flux(i,j) * CS%flux_factor

    if ((sfc_state%ocean_mass(i,j) > CS%col_mass_melt_threshold) .and. &
        (ISS%area_shelf_h(i,j) > 0.0) .and.  (CS%isthermo)) then

      ! Set melt to zero above a cutoff pressure (CS%Rho_ocn*CS%cutoff_depth*CS%g_Earth).
      ! This is needed for the ISOMIP test case.
      if (ISS%mass_shelf(i,j) < CS%Rho_ocn*CS%cutoff_depth) then
        ISS%water_flux(i,j) = 0.0
        fluxes%iceshelf_melt(i,j) = 0.0
      endif
      ! Compute haline driving, which is one of the diags. used in ISOMIP
      haline_driving(i,j) = (ISS%water_flux(i,j) * Sbdry(i,j)) / (CS%Rho_ocn * exch_vel_s(i,j))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!Safety checks !!!!!!!!!!!!!!!!!!!!!!!!!
      !1)Check if haline_driving computed above is consistent with
      ! haline_driving = sfc_state%sss - Sbdry
      !if (fluxes%iceshelf_melt(i,j) /= 0.0) then
      !   if (haline_driving(i,j) /= (sfc_state%sss(i,j) - Sbdry(i,j))) then
      !     write(mesg,*) 'at i,j=',i,j,' haline_driving, sss-Sbdry',haline_driving(i,j), &
      !                   (sfc_state%sss(i,j) - Sbdry(i,j))
      !     call MOM_error(FATAL, &
      !            "shelf_calc_flux: Inconsistency in melt and haline_driving"//trim(mesg))
      !   endif
      !endif

      ! 2) check if |melt| > 0 when ustar_shelf = 0.
      ! this should never happen
      if ((abs(fluxes%iceshelf_melt(i,j))>0.0) .and. (fluxes%ustar_shelf(i,j) == 0.0)) then
        write(mesg,*) "|melt| = ",fluxes%iceshelf_melt(i,j)," > 0 and ustar_shelf = 0. at i,j", i, j
        call MOM_error(FATAL, "shelf_calc_flux: "//trim(mesg))
      endif
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!End of safety checks !!!!!!!!!!!!!!!!!!!
    elseif (ISS%area_shelf_h(i,j) > 0.0) then
      ! This is grounded ice, that could be modified to melt if a geothermal heat flux were used.
      haline_driving(i,j) = 0.0
      ISS%water_flux(i,j) = 0.0
      fluxes%iceshelf_melt(i,j) = 0.0
    endif ! area_shelf_h

    ! mass flux [R Z L2 T-1 ~> kg s-1], part of ISOMIP diags.
    mass_flux(i,j) = ISS%water_flux(i,j) * ISS%area_shelf_h(i,j)
  enddo ; enddo ! i- and j-loops

  if (CS%active_shelf_dynamics .or. CS%override_shelf_movement) then
    call cpu_clock_begin(id_clock_pass)
    call pass_var(ISS%area_shelf_h, G%domain, complete=.false.)
    call pass_var(ISS%mass_shelf, G%domain)
    call cpu_clock_end(id_clock_pass)
  endif

  ! Melting has been computed, now is time to update thickness and mass
  if ( CS%override_shelf_movement .and. (.not.CS%mass_from_file)) then
    call change_thickness_using_melt(ISS, G, US, US%s_to_T*time_step, fluxes, CS%density_ice, CS%debug)

    if (CS%debug) then
      call hchksum(ISS%h_shelf, "h_shelf after change thickness using melt", G%HI, haloshift=0, scale=US%Z_to_m)
      call hchksum(ISS%mass_shelf, "mass_shelf after change thickness using melt", G%HI, haloshift=0, &
                   scale=US%RZ_to_kg_m2)
    endif
  endif

  if (CS%debug) call MOM_forcing_chksum("Before add shelf flux", fluxes, G, CS%US, haloshift=0)

  call add_shelf_flux(G, US, CS, sfc_state, fluxes)

  ! now the thermodynamic data is passed on... time to update the ice dynamic quantities

  if (CS%active_shelf_dynamics) then
    update_ice_vel = .false.
    coupled_GL = (CS%GL_couple .and. .not.CS%solo_ice_sheet)

    ! advect the ice shelf, and advance the front. Calving will be in here somewhere as well..
    ! when we decide on how to do it
    call update_ice_shelf(CS%dCS, ISS, G, US, US%s_to_T*time_step, Time, &
                          sfc_state%ocean_mass, coupled_GL)

  endif

  call enable_averaging(time_step,Time,CS%diag)
  if (CS%id_shelf_mass > 0) call post_data(CS%id_shelf_mass, ISS%mass_shelf, CS%diag)
  if (CS%id_area_shelf_h > 0) call post_data(CS%id_area_shelf_h, ISS%area_shelf_h, CS%diag)
  if (CS%id_ustar_shelf > 0) call post_data(CS%id_ustar_shelf, fluxes%ustar_shelf, CS%diag)
  if (CS%id_melt > 0) call post_data(CS%id_melt, fluxes%iceshelf_melt, CS%diag)
  if (CS%id_thermal_driving > 0) call post_data(CS%id_thermal_driving, (sfc_state%sst-ISS%tfreeze), CS%diag)
  if (CS%id_Sbdry > 0) call post_data(CS%id_Sbdry, Sbdry, CS%diag)
  if (CS%id_haline_driving > 0) call post_data(CS%id_haline_driving, haline_driving, CS%diag)
  if (CS%id_mass_flux > 0) call post_data(CS%id_mass_flux, mass_flux, CS%diag)
  if (CS%id_u_ml > 0) call post_data(CS%id_u_ml, sfc_state%u, CS%diag)
  if (CS%id_v_ml > 0) call post_data(CS%id_v_ml, sfc_state%v, CS%diag)
  if (CS%id_tfreeze > 0) call post_data(CS%id_tfreeze, ISS%tfreeze, CS%diag)
  if (CS%id_tfl_shelf > 0) call post_data(CS%id_tfl_shelf, ISS%tflux_shelf, CS%diag)
  if (CS%id_exch_vel_t > 0) call post_data(CS%id_exch_vel_t, exch_vel_t, CS%diag)
  if (CS%id_exch_vel_s > 0) call post_data(CS%id_exch_vel_s, exch_vel_s, CS%diag)
  if (CS%id_h_shelf > 0) call post_data(CS%id_h_shelf, ISS%h_shelf, CS%diag)
  if (CS%id_h_mask > 0) call post_data(CS%id_h_mask,ISS%hmask,CS%diag)
  call disable_averaging(CS%diag)

  if (present(forces)) then
    call add_shelf_forces(G, US, CS, forces, do_shelf_area=(CS%active_shelf_dynamics .or. &
                                                        CS%override_shelf_movement))
  endif

  call cpu_clock_end(id_clock_shelf)

  if (CS%debug) call MOM_forcing_chksum("End of shelf calc flux", fluxes, G, CS%US, haloshift=0)

end subroutine shelf_calc_flux

!> Changes the thickness (mass) of the ice shelf based on sub-ice-shelf melting
subroutine change_thickness_using_melt(ISS, G, US, time_step, fluxes, density_ice, debug)
  type(ocean_grid_type), intent(inout) :: G  !< The ocean's grid structure.
  type(ice_shelf_state), intent(inout) :: ISS !< A structure with elements that describe
                                              !! the ice-shelf state
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  real,                  intent(in)    :: time_step !< The time step for this update [T ~> s].
  type(forcing),         intent(inout) :: fluxes !< structure containing pointers to any possible
                                                 !! thermodynamic or mass-flux forcing fields.
  real,                  intent(in)    :: density_ice !< The density of ice-shelf ice [R ~> kg m-3].
  logical,     optional, intent(in)    :: debug !< If present and true, write chksums

  ! locals
  real :: I_rho_ice ! Ice specific volume [R-1 ~> m3 kg-1]
  integer :: i, j

  I_rho_ice = 1.0 / density_ice


  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then
      ! first, zero out fluxes applied during previous time step
      if (associated(fluxes%lprec)) fluxes%lprec(i,j) = 0.0
      if (associated(fluxes%sens)) fluxes%sens(i,j) = 0.0
      if (associated(fluxes%frac_shelf_h)) fluxes%frac_shelf_h(i,j) = 0.0
      if (associated(fluxes%salt_flux)) fluxes%salt_flux(i,j) = 0.0

      if (ISS%water_flux(i,j) * time_step / density_ice < ISS%h_shelf(i,j)) then
        ISS%h_shelf(i,j) = ISS%h_shelf(i,j) - ISS%water_flux(i,j) * time_step / density_ice
      else
        ! the ice is about to melt away, so set thickness, area, and mask to zero
        ! NOTE: this is not mass conservative should maybe scale salt & heat flux for this cell
        ISS%h_shelf(i,j) = 0.0
        ISS%hmask(i,j) = 0.0
        ISS%area_shelf_h(i,j) = 0.0
      endif
      ISS%mass_shelf(i,j) = ISS%h_shelf(i,j) * density_ice
    endif
  enddo ; enddo

  call pass_var(ISS%area_shelf_h, G%domain)
  call pass_var(ISS%h_shelf, G%domain)
  call pass_var(ISS%hmask, G%domain)
  call pass_var(ISS%mass_shelf, G%domain)

end subroutine change_thickness_using_melt

!> This subroutine adds the mechanical forcing fields and perhaps shelf areas, based on
!! the ice state in ice_shelf_CS.
subroutine add_shelf_forces(G, US, CS, forces, do_shelf_area)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(ice_shelf_CS),    pointer       :: CS   !< This module's control structure.
  type(mech_forcing),    intent(inout) :: forces !< A structure with the driving mechanical forces
  logical, optional,     intent(in)    :: do_shelf_area !< If true find the shelf-covered areas.

  real :: kv_rho_ice ! The viscosity of ice divided by its density [L4 T-1 R-1 Z-2 ~> m5 kg-1 s-1].
  real :: press_ice  ! The pressure of the ice shelf per unit area of ocean (not ice) [R L2 T-2 ~> Pa].
  logical :: find_area ! If true find the shelf areas at u & v points.
  type(ice_shelf_state), pointer :: ISS => NULL() ! A structure with elements that describe
                                          ! the ice-shelf state

  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed

  if ((CS%grid%isc /= G%isc) .or. (CS%grid%iec /= G%iec) .or. &
      (CS%grid%jsc /= G%jsc) .or. (CS%grid%jec /= G%jec)) &
    call MOM_error(FATAL,"add_shelf_forces: Incompatible ocean and ice shelf grids.")

  ISS => CS%ISS

  find_area = .true. ; if (present(do_shelf_area)) find_area = do_shelf_area

  if (find_area) then
    ! The frac_shelf is set over the widest possible area. Could it be smaller?
    do j=jsd,jed ; do I=isd,ied-1
      forces%frac_shelf_u(I,j) = 0.0
      if ((G%areaT(i,j) + G%areaT(i+1,j) > 0.0)) & ! .and. (G%areaCu(I,j) > 0.0)) &
        forces%frac_shelf_u(I,j) = (ISS%area_shelf_h(i,j) + ISS%area_shelf_h(i+1,j)) / &
                                   (G%areaT(i,j) + G%areaT(i+1,j))
    enddo ; enddo
    do J=jsd,jed-1 ; do i=isd,ied
      forces%frac_shelf_v(i,J) = 0.0
      if ((G%areaT(i,j) + G%areaT(i,j+1) > 0.0)) & ! .and. (G%areaCv(i,J) > 0.0)) &
        forces%frac_shelf_v(i,J) = (ISS%area_shelf_h(i,j) + ISS%area_shelf_h(i,j+1)) / &
                                   (G%areaT(i,j) + G%areaT(i,j+1))
    enddo ; enddo
    call pass_vector(forces%frac_shelf_u, forces%frac_shelf_v, G%domain, TO_ALL, CGRID_NE)
  endif

  do j=js,je ; do i=is,ie
    press_ice = (ISS%area_shelf_h(i,j) * G%IareaT(i,j)) * (CS%g_Earth * ISS%mass_shelf(i,j))
    if (associated(forces%p_surf)) then
      if (.not.forces%accumulate_p_surf) forces%p_surf(i,j) = 0.0
      forces%p_surf(i,j) = forces%p_surf(i,j) + press_ice
    endif
    if (associated(forces%p_surf_full)) then
      if (.not.forces%accumulate_p_surf) forces%p_surf_full(i,j) = 0.0
      forces%p_surf_full(i,j) = forces%p_surf_full(i,j) + press_ice
    endif
  enddo ; enddo

  ! For various reasons, forces%rigidity_ice_[uv] is always updated here. Note
  ! that it may have been zeroed out where IOB is translated to forces and
  ! contributions from icebergs and the sea-ice pack added subsequently.
  !### THE RIGIDITY SHOULD ALSO INCORPORATE AREAL-COVERAGE INFORMATION.
  kv_rho_ice = CS%kv_ice / CS%density_ice
  do j=js,je ; do I=is-1,ie
    if (.not.forces%accumulate_rigidity) forces%rigidity_ice_u(I,j) = 0.0
    forces%rigidity_ice_u(I,j) = forces%rigidity_ice_u(I,j) + &
            kv_rho_ice * min(ISS%mass_shelf(i,j), ISS%mass_shelf(i+1,j))
  enddo ; enddo
  do J=js-1,je ; do i=is,ie
    if (.not.forces%accumulate_rigidity) forces%rigidity_ice_v(i,J) = 0.0
    forces%rigidity_ice_v(i,J) = forces%rigidity_ice_v(i,J) + &
            kv_rho_ice * min(ISS%mass_shelf(i,j), ISS%mass_shelf(i,j+1))
  enddo ; enddo

  if (CS%debug) then
    call uvchksum("rigidity_ice_[uv]", forces%rigidity_ice_u, forces%rigidity_ice_v, &
                  G%HI, symmetric=.true., scale=US%L_to_m**3*US%L_to_Z*US%s_to_T)
    call uvchksum("frac_shelf_[uv]", forces%frac_shelf_u, forces%frac_shelf_v, &
                  G%HI, symmetric=.true.)
  endif

end subroutine add_shelf_forces

!> This subroutine adds the ice shelf pressure to the fluxes type.
subroutine add_shelf_pressure(G, US, CS, fluxes)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(ice_shelf_CS),    intent(in)    :: CS   !< This module's control structure.
  type(forcing),         intent(inout) :: fluxes  !< A structure of surface fluxes that may be updated.

  real :: press_ice       !< The pressure of the ice shelf per unit area of ocean (not ice) [R L2 T-2 ~> Pa].
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if ((CS%grid%isc /= G%isc) .or. (CS%grid%iec /= G%iec) .or. &
      (CS%grid%jsc /= G%jsc) .or. (CS%grid%jec /= G%jec)) &
    call MOM_error(FATAL,"add_shelf_pressure: Incompatible ocean and ice shelf grids.")

  do j=js,je ; do i=is,ie
    press_ice = (CS%ISS%area_shelf_h(i,j) * G%IareaT(i,j)) * (CS%g_Earth * CS%ISS%mass_shelf(i,j))
    if (associated(fluxes%p_surf)) then
      if (.not.fluxes%accumulate_p_surf) fluxes%p_surf(i,j) = 0.0
      fluxes%p_surf(i,j) = fluxes%p_surf(i,j) + press_ice
    endif
    if (associated(fluxes%p_surf_full)) then
      if (.not.fluxes%accumulate_p_surf) fluxes%p_surf_full(i,j) = 0.0
      fluxes%p_surf_full(i,j) = fluxes%p_surf_full(i,j) + press_ice
    endif
  enddo ; enddo

end subroutine add_shelf_pressure

!> Updates surface fluxes that are influenced by sub-ice-shelf melting
subroutine add_shelf_flux(G, US, CS, sfc_state, fluxes)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(ice_shelf_CS),    pointer       :: CS   !< This module's control structure.
  type(surface),         intent(inout) :: sfc_state !< Surface ocean state
  type(forcing),         intent(inout) :: fluxes  !< A structure of surface fluxes that may be used/updated.

  ! local variables
  real :: frac_shelf       !< The fractional area covered by the ice shelf [nondim].
  real :: frac_open        !< The fractional area of the ocean that is not covered by the ice shelf [nondim].
  real :: delta_mass_shelf !< Change in ice shelf mass over one time step [R Z m2 T-1 ~> kg s-1]
  real :: balancing_flux   !< The fresh water flux that balances the integrated melt flux [R Z T-1 ~> kg m-2 s-1]
  real :: balancing_area   !< total area where the balancing flux is applied [m2]
  type(time_type) :: dTime !< The time step as a time_type
  type(time_type) :: Time0 !< The previous time (Time-dt)
  real, dimension(SZDI_(G),SZDJ_(G)) :: bal_frac  !< Fraction of the cel1 where the mass flux
                          !! balancing the net melt flux occurs, 0 to 1 [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)) :: last_mass_shelf !< Ice shelf mass
                          !! at at previous time (Time-dt) [R Z ~> kg m-2]
  real, dimension(SZDI_(G),SZDJ_(G)) :: delta_float_mass   !< The change in the floating mass between
                          !! the two timesteps at (Time) and (Time-dt) [R Z ~> kg m-2].
  real, dimension(SZDI_(G),SZDJ_(G))  :: last_h_shelf !< Ice shelf thickness [Z ~> m]
                          !! at at previous time (Time-dt)
  real, dimension(SZDI_(G),SZDJ_(G))  :: last_hmask !< Ice shelf mask [nondim]
                          !! at at previous time (Time-dt)
  real, dimension(SZDI_(G),SZDJ_(G))  :: last_area_shelf_h !< Ice shelf area [L2 ~> m2]
                          !! at at previous time (Time-dt)
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
                                          !! the ice-shelf state

  character(len=160) :: mesg  ! The text of an error message
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed

  if ((CS%grid%isc /= G%isc) .or. (CS%grid%iec /= G%iec) .or. &
      (CS%grid%jsc /= G%jsc) .or. (CS%grid%jec /= G%jec)) &
    call MOM_error(FATAL,"add_shelf_flux: Incompatible ocean and ice shelf grids.")

  ISS => CS%ISS


  call add_shelf_pressure(G, US, CS, fluxes)

  ! Determine ustar and the square magnitude of the velocity in the
  ! bottom boundary layer. Together these give the TKE source and
  ! vertical decay scale.

  if (CS%debug) then
    if (allocated(sfc_state%taux_shelf) .and. allocated(sfc_state%tauy_shelf)) then
      call uvchksum("tau[xy]_shelf", sfc_state%taux_shelf, sfc_state%tauy_shelf, &
                    G%HI, haloshift=0, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    endif
  endif

  if (CS%active_shelf_dynamics .or. CS%override_shelf_movement) then
    do j=jsd,jed ; do i=isd,ied
      if (G%areaT(i,j) > 0.0) &
        fluxes%frac_shelf_h(i,j) = min(1.0, ISS%area_shelf_h(i,j) * G%IareaT(i,j))
    enddo ; enddo
  endif

  if (CS%debug) then
    call MOM_forcing_chksum("Before adding shelf fluxes", fluxes, G, CS%US, haloshift=0)
  endif

  do j=js,je ; do i=is,ie ; if (ISS%area_shelf_h(i,j) > 0.0) then
    ! Replace fluxes intercepted by the ice shelf with fluxes from the ice shelf
    frac_shelf = min(1.0, ISS%area_shelf_h(i,j) * G%IareaT(i,j))
    frac_open = max(0.0, 1.0 - frac_shelf)

    if (associated(fluxes%sw)) fluxes%sw(i,j) = frac_open * fluxes%sw(i,j)
    if (associated(fluxes%sw_vis_dir)) fluxes%sw_vis_dir(i,j) = frac_open * fluxes%sw_vis_dir(i,j)
    if (associated(fluxes%sw_vis_dif)) fluxes%sw_vis_dif(i,j) = frac_open * fluxes%sw_vis_dif(i,j)
    if (associated(fluxes%sw_nir_dir)) fluxes%sw_nir_dir(i,j) = frac_open * fluxes%sw_nir_dir(i,j)
    if (associated(fluxes%sw_nir_dif)) fluxes%sw_nir_dif(i,j) = frac_open * fluxes%sw_nir_dif(i,j)
    if (associated(fluxes%lw)) fluxes%lw(i,j) = frac_open * fluxes%lw(i,j)
    if (associated(fluxes%latent)) fluxes%latent(i,j) = frac_open * fluxes%latent(i,j)
    if (associated(fluxes%evap)) fluxes%evap(i,j) = frac_open * fluxes%evap(i,j)
    if (associated(fluxes%lprec)) then
      if (ISS%water_flux(i,j) > 0.0) then
        fluxes%lprec(i,j) =  frac_shelf*ISS%water_flux(i,j)*CS%flux_factor + frac_open * fluxes%lprec(i,j)
      else
        fluxes%lprec(i,j) = frac_open * fluxes%lprec(i,j)
        fluxes%evap(i,j) = fluxes%evap(i,j) + frac_shelf*ISS%water_flux(i,j)*CS%flux_factor
      endif
    endif

    if (associated(fluxes%sens)) &
      fluxes%sens(i,j) = frac_shelf*ISS%tflux_ocn(i,j)*CS%flux_factor + frac_open * fluxes%sens(i,j)
    ! The salt flux should be mostly from sea ice, so perhaps none should be intercepted and this should be changed.
    if (associated(fluxes%salt_flux)) &
      fluxes%salt_flux(i,j) = frac_shelf * ISS%salt_flux(i,j)*CS%flux_factor + frac_open * fluxes%salt_flux(i,j)
  endif ; enddo ; enddo

  if (CS%debug) then
    call hchksum(ISS%water_flux, "water_flux add shelf fluxes", G%HI, haloshift=0, scale=US%RZ_T_to_kg_m2s)
    call hchksum(ISS%tflux_ocn, "tflux_ocn add shelf fluxes", G%HI, haloshift=0, scale=US%QRZ_T_to_W_m2)
    call MOM_forcing_chksum("After adding shelf fluxes", fluxes, G, CS%US, haloshift=0)
  endif

  ! Keep sea level constant by removing mass via a balancing flux that might be applied
  ! in the open ocean or the sponge region (via virtual precip, vprec). Apply additional
  ! salt/heat fluxes so that the resultant surface buoyancy forcing is ~ 0.
  ! This is needed for some of the ISOMIP+ experiments.

  if (CS%constant_sea_level) then
    if (.not. associated(fluxes%salt_flux)) allocate(fluxes%salt_flux(ie,je))
    if (.not. associated(fluxes%vprec)) allocate(fluxes%vprec(ie,je))
    fluxes%salt_flux(:,:) = 0.0 ; fluxes%vprec(:,:) = 0.0

    ! take into account changes in mass (or thickness) when imposing ice shelf mass
    if (CS%override_shelf_movement .and. CS%mass_from_file) then
      dTime = real_to_time(CS%time_step)

      ! Compute changes in mass after at least one full time step
      if (CS%Time > dTime) then
        Time0 = CS%Time - dTime
        do j=js,je ; do i=is,ie
          last_hmask(i,j) = ISS%hmask(i,j) ; last_area_shelf_h(i,j) = ISS%area_shelf_h(i,j)
        enddo ; enddo
        call time_interp_external(CS%id_read_mass, Time0, last_mass_shelf)
        do j=js,je ; do i=is,ie
        ! This should only be done if time_interp_external did an update.
          last_mass_shelf(i,j) = US%kg_m3_to_R*US%m_to_Z * last_mass_shelf(i,j) ! Rescale after time_interp
          last_h_shelf(i,j) = last_mass_shelf(i,j) / CS%density_ice
        enddo ; enddo

        ! apply calving
        if (CS%min_thickness_simple_calve > 0.0) then
          call ice_shelf_min_thickness_calve(G, last_h_shelf, last_area_shelf_h, last_hmask, &
                                       CS%min_thickness_simple_calve, halo=0)
          ! convert to mass again
          do j=js,je ; do i=is,ie
            last_mass_shelf(i,j) = last_h_shelf(i,j) * CS%density_ice
          enddo ; enddo
        endif

        ! get total ice shelf mass at (Time-dt) and (Time), in kg
        do j=js,je ; do i=is,ie
          ! Just consider the change in the mass of the floating shelf.
          if ((sfc_state%ocean_mass(i,j) > CS%min_ocean_mass_float) .and. &
              (ISS%area_shelf_h(i,j) > 0.0)) then
            delta_float_mass(i,j) = ISS%mass_shelf(i,j) - last_mass_shelf(i,j)
          else
            delta_float_mass(i,j) = 0.0
          endif
        enddo ; enddo
        delta_mass_shelf = US%kg_m2s_to_RZ_T*(global_area_integral(delta_float_mass, G, scale=US%RZ_to_kg_m2, &
                                                                   area=ISS%area_shelf_h) / CS%time_step)
      else! first time step
        delta_mass_shelf = 0.0
      endif
    else ! ice shelf mass does not change
      delta_mass_shelf = 0.0
    endif

    ! average total melt flux over sponge area
    do j=js,je ; do i=is,ie
      if ((G%mask2dT(i,j) > 0.0) .AND. (ISS%area_shelf_h(i,j) * G%IareaT(i,j) < 1.0)) then
         ! Uncomment this for some ISOMIP cases:
         !  .AND. (G%geoLonT(i,j) >= 790.0) .AND. (G%geoLonT(i,j) <= 800.0)) then
        bal_frac(i,j) = max(1.0 - ISS%area_shelf_h(i,j) * G%IareaT(i,j), 0.0)
      else
        bal_frac(i,j) = 0.0
      endif
    enddo ; enddo

    balancing_area = global_area_integral(bal_frac, G)
    if (balancing_area > 0.0) then
      balancing_flux = ( US%kg_m2s_to_RZ_T*global_area_integral(ISS%water_flux, G, scale=US%RZ_T_to_kg_m2s, &
                                                                area=ISS%area_shelf_h) + &
                         delta_mass_shelf ) / balancing_area
    else
      balancing_flux = 0.0
    endif

    ! apply fluxes
    do j=js,je ; do i=is,ie
      if (bal_frac(i,j) > 0.0) then
        ! evap is negative, and vprec has units of [R Z T-1 ~> kg m-2 s-1]
        fluxes%vprec(i,j) = -balancing_flux
        fluxes%sens(i,j) = fluxes%vprec(i,j) * CS%Cp * CS%T0 ! [ Q R Z T-1 ~> W /m^2 ]
        fluxes%salt_flux(i,j) = fluxes%vprec(i,j) * CS%S0*1.0e-3 ! [kgSalt/kg R Z T-1 ~> kgSalt m-2 s-1]
      endif
    enddo ; enddo

    if (CS%debug) then
      write(mesg,*) 'Balancing flux (kg/(m^2 s)), dt = ', balancing_flux*US%RZ_T_to_kg_m2s, CS%time_step
      call MOM_mesg(mesg)
      call MOM_forcing_chksum("After constant sea level", fluxes, G, CS%US, haloshift=0)
    endif

  endif ! constant_sea_level

end subroutine add_shelf_flux


!> Initializes shelf model data, parameters and diagnostics
subroutine initialize_ice_shelf(param_file, ocn_grid, Time, CS, diag, forces, fluxes, Time_in, solo_ice_sheet_in)
  type(param_file_type),        intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(ocean_grid_type),        pointer       :: ocn_grid   !< The calling ocean model's horizontal grid structure
  type(time_type),              intent(inout) :: Time !< The clock that that will indicate the model time
  type(ice_shelf_CS),           pointer       :: CS   !< A pointer to the ice shelf control structure
  type(diag_ctrl),    target,   intent(in)    :: diag !< A structure that is used to regulate the diagnostic output.
  type(forcing),      optional, intent(inout) :: fluxes !< A structure containing pointers to any possible
                                                   !! thermodynamic or mass-flux forcing fields.
  type(mech_forcing), optional, intent(inout) :: forces !< A structure with the driving mechanical forces
  type(time_type),    optional, intent(in)    :: Time_in !< The time at initialization.
  logical,            optional, intent(in)    :: solo_ice_sheet_in !< If present, this indicates whether
                                                   !! a solo ice-sheet driver.

  type(ocean_grid_type), pointer :: G  => NULL(), OG  => NULL() ! Pointers to grids for convenience.
  type(unit_scale_type), pointer :: US => NULL() ! Pointer to a structure containing
                                                 ! various unit conversion factors
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
                                          !! the ice-shelf state
  type(directories)  :: dirs
  type(dyn_horgrid_type), pointer :: dG => NULL()
  real    :: Z_rescale  ! A rescaling factor for heights from the representation in
                        ! a restart file to the internal representation in this run.
  real    :: RZ_rescale ! A rescaling factor for mass loads from the representation in
                        ! a restart file to the internal representation in this run.
  real    :: L_rescale  ! A rescaling factor for horizontal lengths from the representation in
                        ! a restart file to the internal representation in this run.
  real :: meltrate_conversion ! The conversion factor to use for in the melt rate diagnostic.
  real :: dz_ocean_min_float ! The minimum ocean thickness above which the ice shelf is considered
                        ! to be floating when CONST_SEA_LEVEL = True [Z ~> m].
  real :: cdrag, drag_bg_vel
  logical :: new_sim, save_IC, var_force
  !This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=200) :: config
  character(len=200) :: IC_file,filename,inputdir
  character(len=40)  :: mdl = "MOM_ice_shelf"  ! This module's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq
  integer :: wd_halos(2)
  logical :: read_TideAmp, shelf_mass_is_dynamic, debug
  character(len=240) :: Tideamp_file
  real    :: utide  ! A tidal velocity [L T-1 ~> m s-1]
  real    :: col_thick_melt_thresh ! An ocean column thickness below which iceshelf melting
                       ! does not occur [Z ~> m]
  if (associated(CS)) then
    call MOM_error(FATAL, "MOM_ice_shelf.F90, initialize_ice_shelf: "// &
                          "called with an associated control structure.")
    return
  endif
  allocate(CS)

  !   Go through all of the infrastructure initialization calls, since this is
  ! being treated as an independent component that just happens to use the
  ! MOM's grid and infrastructure.
  call Get_MOM_Input(dirs=dirs)

  ! Determining the internal unit scaling factors for this run.
  call unit_scaling_init(param_file, CS%US)

  ! Set up the ice-shelf domain and grid
  wd_halos(:)=0
  call MOM_domains_init(CS%grid%domain, param_file, min_halo=wd_halos, symmetric=GRID_SYM_)
  ! call diag_mediator_init(CS%grid,param_file,CS%diag)
  ! this needs to be fixed - will probably break when not using coupled driver 0
  call MOM_grid_init(CS%grid, param_file, CS%US)

  call create_dyn_horgrid(dG, CS%grid%HI)
  call clone_MOM_domain(CS%grid%Domain, dG%Domain)

  call set_grid_metrics(dG, param_file, CS%US)
  ! call set_diag_mediator_grid(CS%grid, CS%diag)

  ! The ocean grid possibly uses different symmetry.
  if (associated(ocn_grid)) then ; CS%ocn_grid => ocn_grid
  else ; CS%ocn_grid => CS%grid ; endif

  ! Convenience pointers
  G => CS%grid
  OG => CS%ocn_grid
  US => CS%US

  if (is_root_pe()) then
    write(0,*) 'OG: ', OG%isd, OG%isc, OG%iec, OG%ied, OG%jsd, OG%jsc, OG%jsd, OG%jed
    write(0,*) 'IG: ', G%isd, G%isc, G%iec, G%ied, G%jsd, G%jsc, G%jsd, G%jed
  endif

  CS%diag => diag

  ! Are we being called from the solo ice-sheet driver? When called by the ocean
  ! model solo_ice_sheet_in is not preset.
  CS%solo_ice_sheet = .false.
  if (present(solo_ice_sheet_in)) CS%solo_ice_sheet = solo_ice_sheet_in

  if (present(Time_in)) Time = Time_in

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB

  CS%override_shelf_movement = .false. ; CS%active_shelf_dynamics = .false.

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DEBUG", debug, default=.false.)
  call get_param(param_file, mdl, "DEBUG_IS", CS%debug, &
                 "If true, write verbose debugging messages for the ice shelf.", &
                 default=debug)
  call get_param(param_file, mdl, "DYNAMIC_SHELF_MASS", shelf_mass_is_dynamic, &
                 "If true, the ice sheet mass can evolve with time.", &
                 default=.false.)
  if (shelf_mass_is_dynamic) then
    call get_param(param_file, mdl, "OVERRIDE_SHELF_MOVEMENT", CS%override_shelf_movement, &
                 "If true, user provided code specifies the ice-shelf "//&
                 "movement instead of the dynamic ice model.", default=.false.)
    CS%active_shelf_dynamics = .not.CS%override_shelf_movement
    call get_param(param_file, mdl, "GROUNDING_LINE_INTERPOLATE", CS%GL_regularize, &
                 "If true, regularize the floatation condition at the "//&
                 "grounding line as in Goldberg Holland Schoof 2009.", default=.false.)
    call get_param(param_file, mdl, "GROUNDING_LINE_COUPLE", CS%GL_couple, &
                 "If true, let the floatation condition be determined by "//&
                 "ocean column thickness. This means that update_OD_ffrac "//&
                 "will be called.  GL_REGULARIZE and GL_COUPLE are exclusive.", &
                 default=.false., do_not_log=CS%GL_regularize)
    if (CS%GL_regularize) CS%GL_couple = .false.
  endif

  call get_param(param_file, mdl, "SHELF_THERMO", CS%isthermo, &
                 "If true, use a thermodynamically interactive ice shelf.", &
                 default=.false.)
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", CS%Lat_fusion, &
                 "The latent heat of fusion.", units="J/kg", default=hlf, scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "SHELF_THREE_EQN", CS%threeeq, &
                 "If true, use the three equation expression of "//&
                 "consistency to calculate the fluxes at the ice-ocean "//&
                 "interface.", default=.true.)
  call get_param(param_file, mdl, "SHELF_INSULATOR", CS%insulator, &
                 "If true, the ice shelf is a perfect insulatior "//&
                 "(no conduction).", default=.false.)
  call get_param(param_file, mdl, "MELTING_CUTOFF_DEPTH", CS%cutoff_depth, &
                 "Depth above which the melt is set to zero (it must be >= 0) "//&
                 "Default value won't affect the solution.", units="m", default=0.0, scale=US%m_to_Z)
  if (CS%cutoff_depth < 0.) &
    call MOM_error(WARNING,"Initialize_ice_shelf: MELTING_CUTOFF_DEPTH must be >= 0.")

  call get_param(param_file, mdl, "CONST_SEA_LEVEL", CS%constant_sea_level, &
                 "If true, apply evaporative, heat and salt fluxes in "//&
                 "the sponge region. This will avoid a large increase "//&
                 "in sea level. This option is needed for some of the "//&
                 "ISOMIP+ experiments (Ocean3 and Ocean4). "//&
                 "IMPORTANT: it is not currently possible to do "//&
                 "prefect restarts using this flag.", default=.false.)
  call get_param(param_file, mdl, "MIN_OCEAN_FLOAT_THICK", dz_ocean_min_float, &
                 "The minimum ocean thickness above which the ice shelf is considered to be "//&
                 "floating when CONST_SEA_LEVEL = True.", &
                 default=0.1, units="m", scale=US%m_to_Z, do_not_log=.not.CS%constant_sea_level)

  call get_param(param_file, mdl, "ISOMIP_S_SUR_SPONGE", CS%S0, &
                 "Surface salinity in the restoring region.", &
                default=33.8, units='ppt', do_not_log=.true.)

  call get_param(param_file, mdl, "ISOMIP_T_SUR_SPONGE", CS%T0, &
                "Surface temperature in the restoring region.", &
                default=-1.9, units='degC', do_not_log=.true.)

  call get_param(param_file, mdl, "SHELF_3EQ_GAMMA", CS%const_gamma, &
                 "If true, user specifies a constant nondimensional heat-transfer coefficient "//&
                 "(GAMMA_T_3EQ), from which the default salt-transfer coefficient is set "//&
                 "as GAMMA_T_3EQ/35. This is used with SHELF_THREE_EQN.", default=.false.)
  if (CS%threeeq) then
    call get_param(param_file, mdl, "SHELF_S_ROOT", CS%find_salt_root, &
                 "If SHELF_S_ROOT = True, salinity at the ice/ocean interface (Sbdry) "//&
                 "is computed from a quadratic equation. Otherwise, the previous "//&
                 "interactive method to estimate Sbdry is used.", default=.false.)
  else
    call get_param(param_file, mdl, "SHELF_2EQ_GAMMA_T", CS%gamma_t, &
                 "If SHELF_THREE_EQN is false, this the fixed turbulent "//&
                 "exchange velocity at the ice-ocean interface.", &
                 units="m s-1", scale=US%m_to_Z*US%T_to_s, fail_if_missing=.true.)
  endif
  if (CS%const_gamma .or. CS%find_salt_root) then
    call get_param(param_file, mdl, "SHELF_3EQ_GAMMA_T", CS%Gamma_T_3EQ, &
                 "Nondimensional heat-transfer coefficient.", &
                  units="nondim", default=2.2e-2)
    call get_param(param_file, mdl, "SHELF_3EQ_GAMMA_S", CS%Gamma_S_3EQ, &
                 "Nondimensional salt-transfer coefficient.", &
                 default=CS%Gamma_T_3EQ/35.0, units="nondim")
  endif

  call get_param(param_file, mdl, "ICE_SHELF_MASS_FROM_FILE", &
                 CS%mass_from_file, "Read the mass of the "//&
                 "ice shelf (every time step) from a file.", default=.false.)

  if (CS%find_salt_root) then ! read liquidus coeffs.
     call get_param(param_file, mdl, "TFREEZE_S0_P0", CS%TFr_0_0, &
                 "this is the freezing potential temperature at "//&
                 "S=0, P=0.", units="degC", default=0.0, do_not_log=.true.)
    call get_param(param_file, mdl, "DTFREEZE_DS", CS%dTFr_dS, &
                 "this is the derivative of the freezing potential temperature with salinity.", &
                 units="degC psu-1", default=-0.054, do_not_log=.true.)
    call get_param(param_file, mdl, "DTFREEZE_DP", CS%dTFr_dp, &
                 "this is the derivative of the freezing potential temperature with pressure.", &
                 units="degC Pa-1", default=0.0, scale=US%RL2_T2_to_Pa, do_not_log=.true.)
  endif

  call get_param(param_file, mdl, "G_EARTH", CS%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "C_P", CS%Cp, &
                 "The heat capacity of sea water, approximated as a constant. "//&
                 "The default value is from the TEOS-10 definition of conservative temperature.", &
                 units="J kg-1 K-1", default=3991.86795711963, scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "RHO_0", CS%Rho_ocn, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "C_P_ICE", CS%Cp_ice, &
                 "The heat capacity of ice.", units="J kg-1 K-1", scale=US%J_kg_to_Q, &
                 default=2.10e3)
  if (CS%constant_sea_level) CS%min_ocean_mass_float = dz_ocean_min_float*CS%Rho_ocn

  call get_param(param_file, mdl, "ICE_SHELF_FLUX_FACTOR", CS%flux_factor, &
                 "Non-dimensional factor applied to shelf thermodynamic "//&
                 "fluxes.", units="none", default=1.0)

  call get_param(param_file, mdl, "KV_ICE", CS%kv_ice, &
                 "The viscosity of the ice.", &
                 units="m2 s-1", default=1.0e10, scale=US%Z_to_L**2*US%m_to_L**2*US%T_to_s)
  call get_param(param_file, mdl, "KV_MOLECULAR", CS%kv_molec, &
                 "The molecular kinimatic viscosity of sea water at the "//&
                 "freezing temperature.", units="m2 s-1", default=1.95e-6, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "ICE_SHELF_SALINITY", CS%Salin_ice, &
                 "The salinity of the ice inside the ice shelf.", units="psu", &
                 default=0.0)
  call get_param(param_file, mdl, "ICE_SHELF_TEMPERATURE", CS%Temp_ice, &
                 "The temperature at the center of the ice shelf.", &
                 units = "degC", default=-15.0)
  call get_param(param_file, mdl, "KD_SALT_MOLECULAR", CS%kd_molec_salt, &
                 "The molecular diffusivity of salt in sea water at the "//&
                 "freezing point.", units="m2 s-1", default=8.02e-10, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "KD_TEMP_MOLECULAR", CS%kd_molec_temp, &
                 "The molecular diffusivity of heat in sea water at the "//&
                 "freezing point.", units="m2 s-1", default=1.41e-7, scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "DT_FORCING", CS%time_step, &
                 "The time step for changing forcing, coupling with other "//&
                 "components, or potentially writing certain diagnostics. "//&
                 "The default value is given by DT.", units="s", default=0.0)

  call get_param(param_file, mdl, "COL_THICK_MELT_THRESHOLD", col_thick_melt_thresh, &
                 "The minimum ocean column thickness where melting is allowed.", &
                 units="m", scale=US%m_to_Z, default=0.0)
  CS%col_mass_melt_threshold =  CS%Rho_ocn * col_thick_melt_thresh

  call get_param(param_file, mdl, "READ_TIDEAMP", read_TIDEAMP, &
                 "If true, read a file (given by TIDEAMP_FILE) containing "//&
                 "the tidal amplitude with INT_TIDE_DISSIPATION.", default=.false.)

  call safe_alloc_ptr(CS%utide,isd,ied,jsd,jed)   ; CS%utide(:,:) = 0.0

  if (read_TIDEAMP) then
    call get_param(param_file, mdl, "TIDEAMP_FILE", TideAmp_file, &
                 "The path to the file containing the spatially varying "//&
                 "tidal amplitudes.", &
                 default="tideamp.nc")
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    TideAmp_file = trim(inputdir) // trim(TideAmp_file)
    call MOM_read_data(TideAmp_file, 'tideamp', CS%utide, G%domain, timelevel=1, scale=US%m_s_to_L_T)
  else
    call get_param(param_file, mdl, "UTIDE", utide, &
                 "The constant tidal amplitude used with INT_TIDE_DISSIPATION.", &
                 units="m s-1", default=0.0 , scale=US%m_s_to_L_T)
    CS%utide(:,:) = utide
  endif

  call EOS_init(param_file, CS%eqn_of_state)

  !! new parameters that need to be in MOM_input

  if (CS%active_shelf_dynamics) then

    call get_param(param_file, mdl, "DENSITY_ICE", CS%density_ice, &
                 "A typical density of ice.", units="kg m-3", default=917.0, scale=US%kg_m3_to_R)

    call get_param(param_file, mdl, "INPUT_FLUX_ICE_SHELF", CS%input_flux, &
                 "volume flux at upstream boundary", units="m2 s-1", default=0.)
    call get_param(param_file, mdl, "INPUT_THICK_ICE_SHELF", CS%input_thickness, &
                 "flux thickness at upstream boundary", units="m", default=1000.)
  else
    ! This is here because of inconsistent defaults.  I don't know why.  RWH
    call get_param(param_file, mdl, "DENSITY_ICE", CS%density_ice, &
                 "A typical density of ice.", units="kg m-3", default=900.0, scale=US%kg_m3_to_R)
  endif
  call get_param(param_file, mdl, "MIN_THICKNESS_SIMPLE_CALVE", &
                CS%min_thickness_simple_calve, &
                 "Min thickness rule for the very simple calving law",&
                 units="m", default=0.0, scale=US%m_to_Z)

  call get_param(param_file, mdl, "USTAR_SHELF_BG", CS%ustar_bg, &
                 "The minimum value of ustar under ice shelves.", &
                 units="m s-1", default=0.0, scale=US%m_to_Z*US%T_to_s)
  call get_param(param_file, mdl, "CDRAG_SHELF", cdrag, &
       "CDRAG is the drag coefficient relating the magnitude of "//&
       "the velocity field to the surface stress.", units="nondim", &
       default=0.003)
  CS%cdrag = cdrag
  if (CS%ustar_bg <= 0.0) then
    call get_param(param_file, mdl, "DRAG_BG_VEL_SHELF", drag_bg_vel, &
                 "DRAG_BG_VEL is either the assumed bottom velocity (with "//&
                 "LINEAR_DRAG) or an unresolved  velocity that is "//&
                 "combined with the resolved velocity to estimate the "//&
                 "velocity magnitude.", units="m s-1", default=0.0, scale=US%m_to_Z*US%T_to_s)
    if (CS%cdrag*drag_bg_vel > 0.0) CS%ustar_bg = sqrt(CS%cdrag)*drag_bg_vel
  endif

  ! Allocate and initialize state variables to default values
  call ice_shelf_state_init(CS%ISS, CS%grid)
  ISS => CS%ISS

  ! Allocate the arrays for passing ice-shelf data through the forcing type.
  if (.not. CS%solo_ice_sheet) then
    call MOM_mesg("MOM_ice_shelf.F90, initialize_ice_shelf: allocating fluxes.")
     ! GMM: the following assures that water/heat fluxes are just allocated
     ! when SHELF_THERMO = True. These fluxes are necessary if one wants to
     ! use either ENERGETICS_SFC_PBL (ALE mode) or BULKMIXEDLAYER (layer mode).
    if (present(fluxes)) &
      call allocate_forcing_type(CS%ocn_grid, fluxes, ustar=.true., shelf=.true., &
                                 press=.true., water=CS%isthermo, heat=CS%isthermo)
    if (present(forces)) &
      call allocate_mech_forcing(CS%ocn_grid, forces, ustar=.true., shelf=.true., press=.true.)
  else
    call MOM_mesg("MOM_ice_shelf.F90, initialize_ice_shelf: allocating fluxes in solo mode.")
    if (present(fluxes)) &
      call allocate_forcing_type(G, fluxes, ustar=.true., shelf=.true., press=.true.)
    if (present(forces)) &
      call allocate_mech_forcing(G, forces, ustar=.true., shelf=.true., press=.true.)
  endif

  ! Set up the bottom depth, G%D either analytically or from file
  call MOM_initialize_topography(dG%bathyT, G%max_depth, dG, param_file)
  call rescale_dyn_horgrid_bathymetry(dG, US%Z_to_m)

  ! Set up the Coriolis parameter, G%f, usually analytically.
  call MOM_initialize_rotation(dG%CoriolisBu, dG, param_file, US)
  ! This copies grid elements, including bathyT and CoriolisBu from dG to CS%grid.
  call copy_dyngrid_to_MOM_grid(dG, CS%grid, US)

  call destroy_dyn_horgrid(dG)

  ! Set up the restarts.
  call restart_init(param_file, CS%restart_CSp, "Shelf.res")
  call register_restart_field(ISS%mass_shelf, "shelf_mass", .true., CS%restart_CSp, &
                              "Ice shelf mass", "kg m-2")
  call register_restart_field(ISS%area_shelf_h, "shelf_area", .true., CS%restart_CSp, &
                              "Ice shelf area in cell", "m2")
  call register_restart_field(ISS%h_shelf, "h_shelf", .true., CS%restart_CSp, &
                              "ice sheet/shelf thickness", "m")
  call register_restart_field(US%m_to_Z_restart, "m_to_Z", .false., CS%restart_CSp, &
                              "Height unit conversion factor", "Z meter-1")
  call register_restart_field(US%m_to_L_restart, "m_to_L", .false., CS%restart_CSp, &
                              "Length unit conversion factor", "L meter-1")
  call register_restart_field(US%kg_m3_to_R_restart, "kg_m3_to_R", .false., CS%restart_CSp, &
                              "Density unit conversion factor", "R m3 kg-1")
  if (CS%active_shelf_dynamics) then
    call register_restart_field(ISS%hmask, "h_mask", .true., CS%restart_CSp, &
                                "ice sheet/shelf thickness mask" ,"none")
  endif

  if (CS%active_shelf_dynamics) then
    ! Allocate CS%dCS and specify additional restarts for ice shelf dynamics
    call register_ice_shelf_dyn_restarts(G, param_file, CS%dCS, CS%restart_CSp)
  endif

  !GMM - I think we do not need to save ustar_shelf and iceshelf_melt in the restart file
  !if (.not. CS%solo_ice_sheet) then
  !  call register_restart_field(fluxes%ustar_shelf, "ustar_shelf", .false., CS%restart_CSp, &
  !                              "Friction velocity under ice shelves", "m s-1")
  !endif

  CS%restart_output_dir = dirs%restart_output_dir

  new_sim = .false.
  if ((dirs%input_filename(1:1) == 'n') .and. &
      (LEN_TRIM(dirs%input_filename) == 1)) new_sim = .true.

  if (CS%override_shelf_movement .and. CS%mass_from_file) then

    ! initialize the ids for reading shelf mass from a netCDF
    call initialize_shelf_mass(G, param_file, CS, ISS)

    if (new_sim) then
      ! new simulation, initialize ice thickness as in the static case
      call initialize_ice_thickness(ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, G, US, param_file)

    ! next make sure mass is consistent with thickness
      do j=G%jsd,G%jed ; do i=G%isd,G%ied
        if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then
          ISS%mass_shelf(i,j) = ISS%h_shelf(i,j)*CS%density_ice
        endif
      enddo ; enddo

      if (CS%min_thickness_simple_calve > 0.0) &
        call ice_shelf_min_thickness_calve(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, &
                                           CS%min_thickness_simple_calve)
    endif
  endif

  if (CS%active_shelf_dynamics) then
    ! the only reason to initialize boundary conds is if the shelf is dynamic - MJH

    ! call initialize_ice_shelf_boundary ( CS%u_face_mask_bdry, CS%v_face_mask_bdry, &
    !                                      CS%u_flux_bdry_val, CS%v_flux_bdry_val, &
    !                                      CS%u_bdry_val, CS%v_bdry_val, CS%h_bdry_val, &
    !                                      ISS%hmask, G, param_file)

  endif

  if (new_sim .and. (.not. (CS%override_shelf_movement .and. CS%mass_from_file))) then

    ! This model is initialized internally or from a file.
    call initialize_ice_thickness(ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, G, US, param_file)

    ! next make sure mass is consistent with thickness
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then
        ISS%mass_shelf(i,j) = ISS%h_shelf(i,j)*CS%density_ice
      endif
    enddo ; enddo

  ! else ! Previous block for new_sim=.T., this block restores the state.
  elseif (.not.new_sim) then
    ! This line calls a subroutine that reads the initial conditions from a restart file.
    call MOM_mesg("MOM_ice_shelf.F90, initialize_ice_shelf: Restoring ice shelf from file.")
    call restore_state(dirs%input_filename, dirs%restart_input_dir, Time, &
                       G, CS%restart_CSp)

    if ((US%m_to_Z_restart /= 0.0) .and. (US%m_to_Z_restart /= US%m_to_Z)) then
      Z_rescale = US%m_to_Z / US%m_to_Z_restart
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        ISS%h_shelf(i,j) = Z_rescale * ISS%h_shelf(i,j)
      enddo ; enddo
    endif

    if ((US%m_to_Z_restart*US%kg_m3_to_R_restart /= 0.0) .and. &
        (US%m_to_Z*US%kg_m3_to_R /= US%m_to_Z_restart * US%kg_m3_to_R_restart)) then
      RZ_rescale = US%m_to_Z*US%kg_m3_to_R / (US%m_to_Z_restart * US%kg_m3_to_R_restart)
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        ISS%mass_shelf(i,j) = RZ_rescale * ISS%mass_shelf(i,j)
      enddo ; enddo
    endif

    if ((US%m_to_L_restart /= 0.0) .and. (US%m_to_L_restart /= US%m_to_L)) then
      L_rescale = US%m_to_L / US%m_to_L_restart
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        ISS%area_shelf_h(i,j) = L_rescale**2 * ISS%area_shelf_h(i,j)
      enddo ; enddo
    endif

  endif ! .not. new_sim

  CS%Time = Time

  call cpu_clock_begin(id_clock_pass)
  call pass_var(ISS%area_shelf_h, G%domain)
  call pass_var(ISS%h_shelf, G%domain)
  call pass_var(ISS%mass_shelf, G%domain)
  call pass_var(ISS%hmask, G%domain)
  call pass_var(G%bathyT, G%domain)
  call cpu_clock_end(id_clock_pass)

  do j=jsd,jed ; do i=isd,ied
    if (ISS%area_shelf_h(i,j) > G%areaT(i,j)) then
      call MOM_error(WARNING,"Initialize_ice_shelf: area_shelf_h exceeds G%areaT.")
      ISS%area_shelf_h(i,j) = G%areaT(i,j)
    endif
  enddo ; enddo
  if (present(fluxes)) then ; do j=jsd,jed ; do i=isd,ied
    if (G%areaT(i,j) > 0.0) fluxes%frac_shelf_h(i,j) = ISS%area_shelf_h(i,j) / G%areaT(i,j)
  enddo ; enddo ; endif

  if (CS%debug) then
    call hchksum(fluxes%frac_shelf_h, "IS init: frac_shelf_h", G%HI, haloshift=0)
  endif

  if (present(forces)) &
    call add_shelf_forces(G, US, CS, forces, do_shelf_area=.not.CS%solo_ice_sheet)

  if (present(fluxes)) call add_shelf_pressure(G, US, CS, fluxes)

  if (CS%active_shelf_dynamics .and. .not.CS%isthermo) then
    ISS%water_flux(:,:) = 0.0
  endif

  if (shelf_mass_is_dynamic) &
    call initialize_ice_shelf_dyn(param_file, Time, ISS, CS%dCS, G, US, diag, new_sim, solo_ice_sheet_in)

  call fix_restart_unit_scaling(US)

  call get_param(param_file, mdl, "SAVE_INITIAL_CONDS", save_IC, &
                 "If true, save the ice shelf initial conditions.", &
                 default=.false.)
  if (save_IC) call get_param(param_file, mdl, "SHELF_IC_OUTPUT_FILE", IC_file,&
                 "The name-root of the output file for the ice shelf "//&
                 "initial conditions.", default="MOM_Shelf_IC")

  if (save_IC .and. .not.((dirs%input_filename(1:1) == 'r') .and. &
                          (LEN_TRIM(dirs%input_filename) == 1))) then
    call save_restart(dirs%output_directory, CS%Time, G, &
                      CS%restart_CSp, filename=IC_file)
  endif


  CS%id_area_shelf_h = register_diag_field('ocean_model', 'area_shelf_h', CS%diag%axesT1, CS%Time, &
     'Ice Shelf Area in cell', 'meter-2', conversion=US%L_to_m**2)
  CS%id_shelf_mass = register_diag_field('ocean_model', 'shelf_mass', CS%diag%axesT1, CS%Time, &
     'mass of shelf', 'kg/m^2', conversion=US%RZ_to_kg_m2)
  CS%id_h_shelf = register_diag_field('ocean_model', 'h_shelf', CS%diag%axesT1, CS%Time, &
       'ice shelf thickness', 'm', conversion=US%Z_to_m)
  CS%id_mass_flux = register_diag_field('ocean_model', 'mass_flux', CS%diag%axesT1,&
     CS%Time, 'Total mass flux of freshwater across the ice-ocean interface.', &
     'kg/s', conversion=US%RZ_T_to_kg_m2s*US%L_to_m**2)

  if (CS%const_gamma) then ! use ISOMIP+ eq. with rho_fw = 1000. kg m-3
    meltrate_conversion = 86400.0*365.0*US%Z_to_m*US%s_to_T / (1000.0*US%kg_m3_to_R)
  else ! use original eq.
    meltrate_conversion = 86400.0*365.0*US%Z_to_m*US%s_to_T / CS%density_ice
  endif
  CS%id_melt = register_diag_field('ocean_model', 'melt', CS%diag%axesT1, CS%Time, &
     'Ice Shelf Melt Rate', 'm yr-1', conversion= meltrate_conversion)
  CS%id_thermal_driving = register_diag_field('ocean_model', 'thermal_driving', CS%diag%axesT1, CS%Time, &
     'pot. temp. in the boundary layer minus freezing pot. temp. at the ice-ocean interface.', 'Celsius')
  CS%id_haline_driving = register_diag_field('ocean_model', 'haline_driving', CS%diag%axesT1, CS%Time, &
     'salinity in the boundary layer minus salinity at the ice-ocean interface.', 'psu')
  CS%id_Sbdry = register_diag_field('ocean_model', 'sbdry', CS%diag%axesT1, CS%Time, &
     'salinity at the ice-ocean interface.', 'psu')
  CS%id_u_ml = register_diag_field('ocean_model', 'u_ml', CS%diag%axesCu1, CS%Time, &
     'Eastward vel. in the boundary layer (used to compute ustar)', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_v_ml = register_diag_field('ocean_model', 'v_ml', CS%diag%axesCv1, CS%Time, &
     'Northward vel. in the boundary layer (used to compute ustar)', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_exch_vel_s = register_diag_field('ocean_model', 'exch_vel_s', CS%diag%axesT1, CS%Time, &
     'Sub-shelf salinity exchange velocity', 'm s-1', conversion=US%Z_to_m*US%s_to_T)
  CS%id_exch_vel_t = register_diag_field('ocean_model', 'exch_vel_t', CS%diag%axesT1, CS%Time, &
     'Sub-shelf thermal exchange velocity', 'm s-1' , conversion=US%Z_to_m*US%s_to_T)
  CS%id_tfreeze = register_diag_field('ocean_model', 'tfreeze', CS%diag%axesT1, CS%Time, &
     'In Situ Freezing point at ice shelf interface', 'degC')
  CS%id_tfl_shelf = register_diag_field('ocean_model', 'tflux_shelf', CS%diag%axesT1, CS%Time, &
     'Heat conduction into ice shelf', 'W m-2', conversion=-US%QRZ_T_to_W_m2)
  CS%id_ustar_shelf = register_diag_field('ocean_model', 'ustar_shelf', CS%diag%axesT1, CS%Time, &
     'Fric vel under shelf', 'm/s', conversion=US%Z_to_m*US%s_to_T)
  if (CS%active_shelf_dynamics) then
    CS%id_h_mask = register_diag_field('ocean_model', 'h_mask', CS%diag%axesT1, CS%Time, &
       'ice shelf thickness mask', 'none')
  endif

  id_clock_shelf = cpu_clock_id('Ice shelf', grain=CLOCK_COMPONENT)
  id_clock_pass = cpu_clock_id(' Ice shelf halo updates', grain=CLOCK_ROUTINE)

end subroutine initialize_ice_shelf

!> Initializes shelf mass based on three options (file, zero and user)
subroutine initialize_shelf_mass(G, param_file, CS, ISS, new_sim)

  type(ocean_grid_type), intent(in) :: G   !< The ocean's grid structure.
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(ice_shelf_CS),    pointer    :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state), intent(inout) :: ISS !< The ice shelf state type that is being updated
  logical,     optional, intent(in) :: new_sim !< If present and false, this run is being restarted

  integer :: i, j, is, ie, js, je
  logical :: read_shelf_area, new_sim_2
  character(len=240) :: config, inputdir, shelf_file, filename
  character(len=120) :: shelf_mass_var  ! The name of shelf mass in the file.
  character(len=120) :: shelf_area_var ! The name of shelf area in the file.
  character(len=40)  :: mdl = "MOM_ice_shelf"
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  new_sim_2 = .true. ; if (present(new_sim)) new_sim_2 = new_sim

  call get_param(param_file, mdl, "ICE_SHELF_CONFIG", config, &
                 "A string that specifies how the ice shelf is "//&
                 "initialized. Valid options include:\n"//&
                 " \tfile\t Read from a file.\n"//&
                 " \tzero\t Set shelf mass to 0 everywhere.\n"//&
                 " \tUSER\t Call USER_initialize_shelf_mass.\n", &
                 fail_if_missing=.true.)

  select case ( trim(config) )
    case ("file")

      call time_interp_external_init()

      call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
      inputdir = slasher(inputdir)

      call get_param(param_file, mdl, "SHELF_FILE", shelf_file, &
              "If DYNAMIC_SHELF_MASS = True, OVERRIDE_SHELF_MOVEMENT = True "//&
              "and ICE_SHELF_MASS_FROM_FILE = True, this is the file from "//&
              "which to read the shelf mass and area.", &
               default="shelf_mass.nc")
      call get_param(param_file, mdl, "SHELF_MASS_VAR", shelf_mass_var, &
                 "The variable in SHELF_FILE with the shelf mass.", &
                 default="shelf_mass")
      call get_param(param_file, mdl, "READ_SHELF_AREA", read_shelf_area, &
                 "If true, also read the area covered by ice-shelf from SHELF_FILE.", &
                 default=.false.)

      filename = trim(slasher(inputdir))//trim(shelf_file)
      call log_param(param_file, mdl, "INPUTDIR/SHELF_FILE", filename)

      CS%id_read_mass = init_external_field(filename, shelf_mass_var, &
                          domain=G%Domain%mpp_domain, verbose=CS%debug)

      if (read_shelf_area) then
         call get_param(param_file, mdl, "SHELF_AREA_VAR", shelf_area_var, &
                  "The variable in SHELF_FILE with the shelf area.", &
                  default="shelf_area")

         CS%id_read_area = init_external_field(filename,shelf_area_var, &
                             domain=G%Domain%mpp_domain)
      endif

      if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
           " initialize_shelf_mass: Unable to open "//trim(filename))

    case ("zero")
      do j=js,je ; do i=is,ie
        ISS%mass_shelf(i,j) = 0.0
        ISS%area_shelf_h(i,j) = 0.0
      enddo ; enddo

    case ("USER")
      call USER_initialize_shelf_mass(ISS%mass_shelf, ISS%area_shelf_h, &
                   ISS%h_shelf, ISS%hmask, G, CS%US, CS%user_CS, param_file, new_sim_2)

    case default ;  call MOM_error(FATAL,"initialize_ice_shelf: "// &
      "Unrecognized ice shelf setup "//trim(config))
  end select

end subroutine initialize_shelf_mass

!> Updates the ice shelf mass using data from a file.
subroutine update_shelf_mass(G, US, CS, ISS, Time)
  type(ocean_grid_type), intent(inout) :: G   !< The ocean's grid structure.
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(ice_shelf_CS),    intent(in)    :: CS  !< A pointer to the ice shelf control structure
  type(ice_shelf_state), intent(inout) :: ISS !< The ice shelf state type that is being updated
  type(time_type),       intent(in)    :: Time !< The current model time

  ! local variables
  integer :: i, j, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  call time_interp_external(CS%id_read_mass, Time, ISS%mass_shelf)
  ! This should only be done if time_interp_external did an update.
  do j=js,je ; do i=is,ie
    ISS%mass_shelf(i,j) = US%kg_m3_to_R*US%m_to_Z * ISS%mass_shelf(i,j) ! Rescale after time_interp
  enddo ; enddo

  do j=js,je ; do i=is,ie
    ISS%area_shelf_h(i,j) = 0.0
    ISS%hmask(i,j) = 0.
    if (ISS%mass_shelf(i,j) > 0.0) then
      ISS%area_shelf_h(i,j) = G%areaT(i,j)
      ISS%h_shelf(i,j) = ISS%mass_shelf(i,j) / CS%density_ice
      ISS%hmask(i,j) = 1.
    endif
  enddo ; enddo

  !call USER_update_shelf_mass(ISS%mass_shelf, ISS%area_shelf_h, ISS%h_shelf, &
  !                            ISS%hmask, CS%grid, CS%user_CS, Time, .true.)

  if (CS%min_thickness_simple_calve > 0.0) then
    call ice_shelf_min_thickness_calve(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, &
                                       CS%min_thickness_simple_calve, halo=0)
  endif

  call pass_var(ISS%area_shelf_h, G%domain)
  call pass_var(ISS%h_shelf, G%domain)
  call pass_var(ISS%hmask, G%domain)
  call pass_var(ISS%mass_shelf, G%domain)

end subroutine update_shelf_mass

!> Save the ice shelf restart file
subroutine ice_shelf_save_restart(CS, Time, directory, time_stamped, filename_suffix)
  type(ice_shelf_CS),         pointer    :: CS !< ice shelf control structure
  type(time_type),            intent(in) :: Time !< model time at this call
  character(len=*), optional, intent(in) :: directory !< An optional directory into which to write
                                               !! these restart files.
  logical,          optional, intent(in) :: time_stamped !< f true, the restart file names include
                                               !! a unique time stamp.  The default is false.
  character(len=*), optional, intent(in) :: filename_suffix !< An optional suffix (e.g., a
                                               !! time-stamp) to append to the restart file names.
  ! local variables
  type(ocean_grid_type), pointer :: G => NULL()
  character(len=200) :: restart_dir

  G => CS%grid

  if (present(directory)) then ; restart_dir = directory
  else ; restart_dir = CS%restart_output_dir ; endif

  call save_restart(restart_dir, Time, CS%grid, CS%restart_CSp, time_stamped)

end subroutine ice_shelf_save_restart

!> Deallocates all memory associated with this module
subroutine ice_shelf_end(CS)
  type(ice_shelf_CS), pointer   :: CS !< A pointer to the ice shelf control structure

  if (.not.associated(CS)) return

  call ice_shelf_state_end(CS%ISS)

  if (CS%active_shelf_dynamics) call ice_shelf_dyn_end(CS%dCS)

  deallocate(CS)

end subroutine ice_shelf_end

!> This routine is for stepping a stand-alone ice shelf model without an ocean.
subroutine solo_step_ice_shelf(CS, time_interval, nsteps, Time, min_time_step_in)
  type(ice_shelf_CS), pointer    :: CS      !< A pointer to the ice shelf control structure
  type(time_type), intent(in)    :: time_interval !< The time interval for this update [s].
  integer,         intent(inout) :: nsteps  !< The running number of ice shelf steps.
  type(time_type), intent(inout) :: Time    !< The current model time
  real,  optional, intent(in)    :: min_time_step_in !< The minimum permitted time step [T ~> s].

  type(ocean_grid_type), pointer :: G => NULL()  ! A pointer to the ocean's grid structure
  type(unit_scale_type), pointer :: US => NULL() ! Pointer to a structure containing
                                                 ! various unit conversion factors
  type(ice_shelf_state), pointer :: ISS => NULL() !< A structure with elements that describe
                                          !! the ice-shelf state
  real :: remaining_time    ! The remaining time in this call [T ~> s]
  real :: time_step         ! The internal time step during this call [T ~> s]
  real :: min_time_step     ! The minimal required timestep that would indicate a fatal problem [T ~> s]
  character(len=240) :: mesg
  logical :: update_ice_vel ! If true, it is time to update the ice shelf velocities.
  logical :: coupled_GL     ! If true the grouding line position is determined based on
                            ! coupled ice-ocean dynamics.
  integer :: is, iec, js, jec, i, j

  G => CS%grid
  US => CS%US
  ISS => CS%ISS
  is = G%isc ; iec = G%iec ; js = G%jsc ; jec = G%jec

  remaining_time = US%s_to_T*time_type_to_real(time_interval)

  if (present (min_time_step_in)) then
    min_time_step = min_time_step_in
  else
    min_time_step = 1000.0*US%s_to_T ! At 1 km resolution this would imply ice is moving at ~1 meter per second
  endif

  write (mesg,*) "TIME in ice shelf call, yrs: ", time_type_to_real(Time)/(365. * 86400.)
  call MOM_mesg("solo_step_ice_shelf: "//mesg, 5)

  do while (remaining_time > 0.0)
    nsteps = nsteps+1

    ! If time_interval is not too long, this is unnecessary.
    time_step = min(ice_time_step_CFL(CS%dCS, ISS, G), remaining_time)

    write (mesg,*) "Ice model timestep = ", US%T_to_s*time_step, " seconds"
    if ((time_step < min_time_step) .and. (time_step < remaining_time))  then
      call MOM_error(FATAL, "MOM_ice_shelf:solo_step_ice_shelf: abnormally small timestep "//mesg)
    else
      call MOM_mesg("solo_step_ice_shelf: "//mesg, 5)
    endif

    remaining_time = remaining_time - time_step

    ! If the last mini-timestep is a day or less, we cannot expect velocities to change by much.
    ! Do not update the velocities if the last step is very short.
    update_ice_vel = ((time_step > min_time_step) .or. (remaining_time > 0.0))
    coupled_GL = .false.

    call update_ice_shelf(CS%dCS, ISS, G, US, time_step, Time, must_update_vel=update_ice_vel)

    call enable_averages(time_step, Time, CS%diag)
    if (CS%id_area_shelf_h > 0) call post_data(CS%id_area_shelf_h, ISS%area_shelf_h, CS%diag)
    if (CS%id_h_shelf > 0) call post_data(CS%id_h_shelf, ISS%h_shelf, CS%diag)
    if (CS%id_h_mask > 0) call post_data(CS%id_h_mask, ISS%hmask, CS%diag)
    call disable_averaging(CS%diag)

  enddo

end subroutine solo_step_ice_shelf

!> \namespace mom_ice_shelf
!!
!! \section section_ICE_SHELF
!!
!! This module implements the thermodynamic aspects of ocean/ice-shelf
!! inter-actions using the MOM framework and coding style.
!!
!! Derived from code by Chris Little, early 2010.
!!
!!   The ice-sheet dynamics subroutines do the following:
!!  initialize_shelf_mass - Initializes the ice shelf mass distribution.
!!      - Initializes h_shelf, h_mask, area_shelf_h
!!      - CURRENTLY: initializes mass_shelf as well, but this is unnecessary, as mass_shelf is initialized based on
!!             h_shelf and density_ice immediately afterwards. Possibly subroutine should be renamed
!!  update_shelf_mass - updates ice shelf mass via netCDF file
!!                      USER_update_shelf_mass (TODO).
!!    solo_step_ice_shelf - called only in ice-only mode.
!!    shelf_calc_flux - after melt rate & fluxes are calculated, ice dynamics are done. currently mass_shelf is
!! updated immediately after ice_shelf_advect in fully dynamic mode.
!!
!!   NOTES: be aware that hmask(:,:) has a number of functions; it is used for front advancement,
!! for subroutines in the velocity solve, and for thickness boundary conditions (this last one may be removed).
!! in other words, interfering with its updates will have implications you might not expect.
!!
!!  Overall issues: Many variables need better documentation and units and the
!!                  subgrid on which they are discretized.
!!
!! \subsection section_ICE_SHELF_equations ICE_SHELF equations
!!
!! The three fundamental equations are:
!! Heat flux
!! \f[ \qquad \rho_w  C_{pw} \gamma_T (T_w - T_b) = \rho_i  \dot{m}  L_f \f]
!! Salt flux
!! \f[  \qquad \rho_w \gamma_s (S_w - S_b) =  \rho_i \dot{m} S_b \f]
!! Freezing temperature
!! \f[  \qquad T_b = a S_b + b + c P \f]
!!
!! where ....
!!
!! \subsection section_ICE_SHELF_references References
!!
!! Asay-Davis, Xylar S., Stephen L. Cornford, Benjamin K. Galton-Fenzi, Rupert M. Gladstone, G. Hilmar Gudmundsson,
!! David M. Holland, Paul R. Holland, and Daniel F. Martin. Experimental design for three interrelated marine ice sheet
!! and ocean model intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1).
!! Geoscientific Model Development 9, no. 7 (2016): 2471.
!!
!! Goldberg, D. N., et al. Investigation of land ice-ocean interaction with a fully coupled ice-ocean model: 1.
!!  Model description and behavior. Journal of Geophysical Research: Earth Surface 117.F2 (2012).
!!
!! Goldberg, D. N., et al. Investigation of land ice-ocean interaction with a fully coupled ice-ocean model: 2.
!! Sensitivity to external forcings. Journal of Geophysical Research: Earth Surface 117.F2 (2012).
!!
!! Holland, David M., and Adrian Jenkins. Modeling thermodynamic ice-ocean interactions at the base of an ice shelf.
!! Journal of Physical Oceanography 29.8 (1999): 1787-1800.

end module MOM_ice_shelf
