!> \brief Parameterization of mixed layer restratification by unresolved mixed-layer eddies.
module MOM_mixed_layer_restrat

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_diag_mediator, only : diag_update_remap_grids
use MOM_domains,       only : pass_var, To_West, To_South, Omit_Corners
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_file_parser,   only : openParameterBlock, closeParameterBlock
use MOM_forcing_type,  only : mech_forcing, find_ustar
use MOM_grid,          only : ocean_grid_type
use MOM_hor_index,     only : hor_index_type
use MOM_interpolate,   only : init_external_field, time_interp_external, time_interp_external_init
use MOM_interpolate,   only : external_field
use MOM_intrinsic_functions, only : cuberoot
use MOM_io,            only : slasher, MOM_read_data
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_restart,       only : register_restart_field, query_initialized, MOM_restart_CS
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type, get_thickness_units
use MOM_EOS,           only : calculate_density, calculate_spec_vol, EOS_domain

implicit none ; private

#include <MOM_memory.h>

public mixedlayer_restrat
public mixedlayer_restrat_init
public mixedlayer_restrat_register_restarts
public mixedlayer_restrat_unit_tests

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for mom_mixed_layer_restrat
type, public :: mixedlayer_restrat_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: ml_restrat_coef       !< A non-dimensional factor by which the instability is enhanced
                                   !! over what would be predicted based on the resolved gradients
                                   !! [nondim].  This increases with grid spacing^2, up to something
                                   !! of order 500.
  real    :: ml_restrat_coef2      !< As for ml_restrat_coef but using the slow filtered MLD [nondim].
  real    :: front_length          !< If non-zero, is the frontal-length scale [L ~> m] used to calculate the
                                   !! upscaling of buoyancy gradients that is otherwise represented
                                   !! by the parameter FOX_KEMPER_ML_RESTRAT_COEF. If MLE_FRONT_LENGTH is
                                   !! non-zero, it is recommended to set FOX_KEMPER_ML_RESTRAT_COEF=1.0.
  logical :: fl_from_file          !< If true, read the MLE front-length scale from a netCDF file.
  logical :: MLE_use_PBL_MLD       !< If true, use the MLD provided by the PBL parameterization.
                                   !! if false, MLE will calculate a MLD based on a density difference
                                   !! based on the parameter MLE_DENSITY_DIFF.
  logical :: Bodner_detect_MLD        !< If true, detect the MLD based on given density difference criterion
                                   !! (MLE_DENSITY_DIFF) in the Bodner et al. parameterization.
  real    :: vonKar                !< The von Karman constant as used for mixed layer viscosity [nondim]
  real    :: MLE_MLD_decay_time    !< Time-scale to use in a running-mean when MLD is retreating [T ~> s].
  real    :: MLE_MLD_decay_time2   !< Time-scale to use in a running-mean when filtered MLD is retreating [T ~> s].
  real    :: MLE_density_diff      !< Density difference used in detecting mixed-layer depth [R ~> kg m-3].
  real    :: MLE_tail_dh           !< Fraction by which to extend the mixed-layer restratification
                                   !! depth used for a smoother stream function at the base of
                                   !! the mixed-layer [nondim].
  real    :: MLE_MLD_stretch       !< A scaling coefficient for stretching/shrinking the MLD used in
                                   !! the MLE scheme [nondim]. This simply multiplies MLD wherever used.

  ! The following parameters are used in the Bodner et al., 2023, parameterization
  logical :: use_Bodner = .false.  !< If true, use the Bodner et al., 2023, parameterization.
  real    :: Cr                    !< Efficiency coefficient from Bodner et al., 2023 [nondim]
  real    :: mstar                 !< The m* value used to estimate the turbulent vertical momentum flux [nondim]
  real    :: nstar                 !< The n* value used to estimate the turbulent vertical momentum flux [nondim]
  real    :: min_wstar2            !< The minimum lower bound to apply to the vertical momentum flux,
                                   !! w'u', in the Bodner et al., restratification parameterization
                                   !! [Z2 T-2 ~> m2 s-2].  This avoids a division-by-zero in the limit when u*
                                   !! and the buoyancy flux are zero.
  real    :: BLD_growing_Tfilt     !< The time-scale for a running-mean filter applied to the boundary layer
                                   !! depth (BLD) when the BLD is deeper than the running mean [T ~> s].
                                   !! A value of 0 instantaneously sets the running mean to the current value of BLD.
  real    :: BLD_decaying_Tfilt    !< The time-scale for a running-mean filter applied to the boundary layer
                                   !! depth (BLD) when the BLD is shallower than the running mean [T ~> s].
                                   !! A value of 0 instantaneously sets the running mean to the current value of BLD.
  real    :: MLD_decaying_Tfilt    !< The time-scale for a running-mean filter applied to the time-filtered
                                   !! MLD, when the latter is shallower than the running mean [T ~> s].
                                   !! A value of 0 instantaneously sets the running mean to the current value of MLD.
  real    :: MLD_growing_Tfilt     !< The time-scale for a running-mean filter applied to the time-filtered
                                   !! MLD, when the latter is deeper than the running mean [T ~> s].
                                   !! A value of 0 instantaneously sets the running mean to the current value of MLD.
  integer :: answer_date           !< The vintage of the order of arithmetic and expressions in the
                                   !! mixed layer restrat calculations.  Values below 20240201 recover
                                   !! the answers from the end of 2023, while higher values use the new
                                   !! cuberoot function in the Bodner code to avoid needing to undo
                                   !! dimensional rescaling.

  logical :: debug = .false.       !< If true, calculate checksums of fields for debugging.

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  type(external_field) :: sbc_fl   !< A handle used in time interpolation of
                                   !! front-length scales read from a file.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  logical :: use_Stanley_ML        !< If true, use the Stanley parameterization of SGS T variance
  real    :: ustar_min             !< A minimum value of ustar in thickness units to avoid numerical
                                   !! problems [H T-1 ~> m s-1 or kg m-2 s-1]
  real    :: Kv_restrat            !< A viscosity that sets a floor on the momentum mixing rate
                                   !! during restratification, rescaled into thickness-based
                                   !! units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1]
  logical :: MLD_grid              !< If true, read a spacially varying field for MLD_decaying_Tfilt
  logical :: Cr_grid               !< If true, read a spacially varying field for Cr

  real, dimension(:,:), allocatable :: &
         MLD_filtered, &           !< Time-filtered MLD [H ~> m or kg m-2]
         MLD_filtered_slow, &      !< Slower time-filtered MLD [H ~> m or kg m-2]
         wpup_filtered, &          !< Time-filtered vertical momentum flux [H L T-2 ~> m2 s-2 or kg m-1 s-2]
         MLD_Tfilt_space, &        !< Spatially varying time scale for MLD filter [T ~> s]
         Cr_space                  !< Spatially varying Cr coefficient [nondim]

  !>@{
  !! Diagnostic identifier
  integer :: id_urestrat_time = -1
  integer :: id_vrestrat_time = -1
  integer :: id_uhml = -1
  integer :: id_vhml = -1
  integer :: id_MLD  = -1
  integer :: id_BLD  = -1
  integer :: id_Rml  = -1
  integer :: id_uDml = -1
  integer :: id_vDml = -1
  integer :: id_uml  = -1
  integer :: id_vml  = -1
  integer :: id_wpup = -1
  integer :: id_ustar = -1
  integer :: id_bflux = -1
  integer :: id_lfbod = -1
  integer :: id_mle_fl = -1
  !>@}

end type mixedlayer_restrat_CS

character(len=40)  :: mdl = "MOM_mixed_layer_restrat" !< This module's name.

contains

!> Driver for the mixed-layer restratification parameterization.
!! The code branches between two different implementations depending
!! on whether the bulk-mixed layer or a general coordinate are in use.
subroutine mixedlayer_restrat(h, uhtr, vhtr, tv, forces, dt, MLD, h_MLD, bflux, VarMix, G, GV, US, CS)
  type(ocean_grid_type),                      intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !! [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !! [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                         intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  real, dimension(:,:),                       pointer       :: MLD    !< Mixed layer depth provided by the
                                                                      !! planetary boundary layer scheme [Z ~> m]
  real, dimension(:,:),                       pointer       :: h_MLD  !< Mixed layer thickness provided
                                                                      !! by the planetary boundary layer
                                                                      !! scheme [H ~> m or kg m-2]
  real, dimension(:,:),                       pointer       :: bflux  !< Surface buoyancy flux provided by the
                                                                      !! PBL scheme [Z2 T-3 ~> m2 s-3]
  type(VarMix_CS),                            intent(in)    :: VarMix !< Variable mixing control structure
  type(mixedlayer_restrat_CS),                intent(inout) :: CS     !< Module control structure


  if (.not. CS%initialized) call MOM_error(FATAL, "mixedlayer_restrat: "// &
         "Module must be initialized before it is used.")

  if (GV%nkml>0) then
    ! Original form, written for the isopycnal model with a bulk mixed layer
    call mixedlayer_restrat_BML(h, uhtr, vhtr, tv, forces, dt, G, GV, US, CS)
  elseif (CS%use_Bodner) then
    ! Implementation of Bodner et al., 2023
    call mixedlayer_restrat_Bodner(CS, G, GV, US, h, uhtr, vhtr, tv, forces, dt, MLD, h_MLD, bflux)
  else
    ! Implementation of Fox-Kemper et al., 2008, to work in general coordinates
    call mixedlayer_restrat_OM4(h, uhtr, vhtr, tv, forces, dt, h_MLD, VarMix, G, GV, US, CS)
  endif

end subroutine mixedlayer_restrat

!> Calculates a restratifying flow in the mixed layer, following the formulation used in OM4
subroutine mixedlayer_restrat_OM4(h, uhtr, vhtr, tv, forces, dt, h_MLD, VarMix, G, GV, US, CS)
  ! Arguments
  type(ocean_grid_type),                      intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                         intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  real, dimension(:,:),                       pointer       :: h_MLD  !< Thickness of water within the
                                                                      !! mixed layer depth provided by
                                                                      !!  the PBL scheme [H ~> m or kg m-2]
  type(VarMix_CS),                            intent(in)    :: VarMix !< Variable mixing control structure
  type(mixedlayer_restrat_CS),                intent(inout) :: CS     !< Module control structure

  ! Local variables
  real :: uhml(SZIB_(G),SZJ_(G),SZK_(GV)) ! Restratifying zonal thickness transports [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vhml(SZI_(G),SZJB_(G),SZK_(GV)) ! Restratifying meridional thickness transports [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    h_avail               ! The volume available for diffusion out of each face of each
                          ! sublayer of the mixed layer, divided by dt [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    U_star_2d, &          ! The wind friction velocity in thickness-based units, calculated using
                          ! the Boussinesq reference density or the time-evolving surface density
                          ! in non-Boussinesq mode [H T-1 ~> m s-1 or kg m-2 s-1]
    MLD_fast, &           ! Mixed layer depth actually used in MLE restratification parameterization [H ~> m or kg m-2]
    htot_fast, &          ! The sum of the thicknesses of layers in the mixed layer [H ~> m or kg m-2]
    Rml_av_fast, &        ! Negative g_Rho0 times the average mixed layer density or G_Earth
                          ! times the average specific volume [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2]
    MLD_slow,  &          ! Mixed layer depth actually used in MLE restratification parameterization [H ~> m or kg m-2]
    mle_fl_2d, &          ! MLE frontal length-scale                                [L ~> m]
    htot_slow, &          ! The sum of the thicknesses of layers in the mixed layer [H ~> m or kg m-2]
    Rml_av_slow           ! Negative g_Rho0 times the average mixed layer density or G_Earth
                          ! times the average specific volume [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2]
  real :: g_Rho0          ! G_Earth/Rho0 times a thickness conversion factor
                          ! [L2 H-1 T-2 R-1 ~> m4 s-2 kg-1 or m7 s-2 kg-2]
  real :: rho_ml(SZI_(G)) ! Potential density relative to the surface [R ~> kg m-3]
  real :: rml_int_fast(SZI_(G)) ! The integral of density over the mixed layer depth [R H ~> kg m-2 or kg2 m-3]
  real :: rml_int_slow(SZI_(G)) ! The integral of density over the mixed layer depth [R H ~> kg m-2 or kg2 m-3]
  real :: SpV_ml(SZI_(G)) ! Specific volume evaluated at the surface pressure [R-1 ~> m3 kg-1]
  real :: SpV_int_fast(SZI_(G)) ! Specific volume integrated through the mixed layer [H R-1 ~> m4 kg-1 or m]
  real :: SpV_int_slow(SZI_(G)) ! Specific volume integrated through the mixed layer [H R-1 ~> m4 kg-1 or m]
  real :: p0(SZI_(G))     ! A pressure of 0 [R L2 T-2 ~> Pa]

  real :: h_vel           ! htot interpolated onto velocity points [H ~> m or kg m-2]
  real :: absf            ! absolute value of f, interpolated to velocity points [T-1 ~> s-1]
  real :: u_star          ! surface friction velocity, interpolated to velocity points and recast into
                          ! thickness-based units [H T-1 ~> m s-1 or kg m-2 s-1].
  real :: mom_mixrate     ! rate at which momentum is homogenized within mixed layer [T-1 ~> s-1]
  real :: timescale       ! mixing growth timescale [T ~> s]
  real :: h_min           ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  real :: h_neglect       ! tiny thickness usually lost in roundoff so can be neglected [H ~> m or kg m-2]
  real :: I4dt            ! 1/(4 dt) [T-1 ~> s-1]
  real :: Ihtot,Ihtot_slow! Inverses of the total mixed layer thickness [H-1 ~> m-1 or m2 kg-1]
  real :: a(SZK_(GV))     ! A non-dimensional value relating the overall flux
                          ! magnitudes (uDml & vDml) to the realized flux in a
                          ! layer [nondim].  The vertical sum of a() through the pieces of
                          ! the mixed layer must be 0.
  real :: b(SZK_(GV))     ! As for a(k) but for the slow-filtered MLD [nondim]
  real :: uDml(SZIB_(G))  ! Zonal volume fluxes in the upper half of the mixed layer [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vDml(SZI_(G))   ! Meridional volume fluxes in the upper half of the mixed layer [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: uDml_slow(SZIB_(G)) ! Zonal volume fluxes in the upper half of the boundary layer to
                          ! restratify the time-filtered boundary layer depth [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vDml_slow(SZI_(G))  ! Meridional volume fluxes in the upper half of the boundary layer to
                          ! restratify the time-filtered boundary layer depth [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: utimescale_diag(SZIB_(G),SZJ_(G)) ! Zonal restratification timescale [T ~> s], stored for diagnostics.
  real :: vtimescale_diag(SZI_(G),SZJB_(G)) ! Meridional restratification timescale [T ~> s], stored for diagnostics.
  real :: uDml_diag(SZIB_(G),SZJ_(G))  ! A 2D copy of uDml for diagnostics [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vDml_diag(SZI_(G),SZJB_(G))  ! A 2D copy of vDml for diagnostics [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G)) :: covTS, & ! SGS TS covariance in Stanley param; currently 0 [C S ~> degC ppt]
                              varS     ! SGS S variance in Stanley param; currently 0    [S2 ~> ppt2]
  real :: aFac, bFac ! Nondimensional ratios [nondim]
  real :: hAtVel    ! Thickness at the velocity points [H ~> m or kg m-2]
  real :: zpa       ! Fractional position within the mixed layer of the interface above a layer [nondim]
  real :: zpb       ! Fractional position within the mixed layer of the interface below a layer [nondim]
  real :: dh        ! Portion of the layer thickness that is in the mixed layer [H ~> m or kg m-2]
  real :: res_scaling_fac ! The resolution-dependent scaling factor [nondim]
  real :: lfront    ! Frontal length scale at velocity points [L ~> m]
  real :: I_LFront  ! The inverse of the frontal length scale [L-1 ~> m-1]
  real :: vonKar_x_pi2    ! A scaling constant that is approximately the von Karman constant times
                          ! pi squared [nondim]
  character(len=128) :: mesg
  logical :: line_is_empty, keep_going, res_upscale
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  h_min = 0.5*GV%Angstrom_H ! This should be GV%Angstrom_H, but that value would change answers.
  covTS(:) = 0.0 !!Functionality not implemented yet; in future, should be passed in tv
  varS(:) = 0.0
  mle_fl_2d(:,:) = 0.0
  vonKar_x_pi2 = CS%vonKar * 9.8696

  if (.not.associated(tv%eqn_of_state)) call MOM_error(FATAL, "mixedlayer_restrat_OM4: "// &
         "An equation of state must be used with this module.")
  if (.not. allocated(VarMix%Rd_dx_h) .and. CS%front_length > 0.) &
    call MOM_error(FATAL, "mixedlayer_restrat_OM4: "// &
         "The resolution argument, Rd/dx, was not associated.")
  if (CS%use_Stanley_ML .and. .not.GV%Boussinesq) call MOM_error(FATAL, &
       "MOM_mixedlayer_restrat: The Stanley parameterization is not"//&
       "available without the Boussinesq approximation.")

  ! Extract the friction velocity from the forcing type.
  call find_ustar(forces, tv, U_star_2d, G, GV, US, halo=1, H_T_units=.true.)

  if (CS%MLE_density_diff > 0.) then ! We need to calculate a mixed layer depth, MLD.
    call detect_mld(h, tv, MLD_fast, G, GV, CS)
  elseif (CS%MLE_use_PBL_MLD) then
    do j=js-1,je+1 ; do i=is-1,ie+1
      MLD_fast(i,j) = CS%MLE_MLD_stretch * h_MLD(i,j)
    enddo ; enddo
  else
    call MOM_error(FATAL, "mixedlayer_restrat_OM4: "// &
         "No MLD to use for MLE parameterization.")
  endif

  ! Apply time filter (to remove diurnal cycle)
  if (CS%MLE_MLD_decay_time>0.) then
    if (CS%debug) then
      call hchksum(CS%MLD_filtered, 'mixed_layer_restrat: MLD_filtered', G%HI, haloshift=1, unscale=GV%H_to_mks)
      call hchksum(h_MLD, 'mixed_layer_restrat: MLD in', G%HI, haloshift=1, unscale=GV%H_to_mks)
    endif
    aFac = CS%MLE_MLD_decay_time / ( dt + CS%MLE_MLD_decay_time )
    bFac = dt / ( dt + CS%MLE_MLD_decay_time )
    do j=js-1,je+1 ; do i=is-1,ie+1
      ! Expression bFac*MLD_fast(i,j) + aFac*CS%MLD_filtered(i,j) is the time-filtered
      ! (running mean) of MLD. The max() allows the "running mean" to be reset
      ! instantly to a deeper MLD.
      CS%MLD_filtered(i,j) = max( MLD_fast(i,j), bFac*MLD_fast(i,j) + aFac*CS%MLD_filtered(i,j) )
      MLD_fast(i,j) = CS%MLD_filtered(i,j)
    enddo ; enddo
  endif

  ! Apply slower time filter (to remove seasonal cycle) on already filtered MLD_fast
  if (CS%MLE_MLD_decay_time2>0.) then
    if (CS%debug) then
      call hchksum(CS%MLD_filtered_slow, 'mixed_layer_restrat: MLD_filtered_slow', G%HI, &
                   haloshift=1, unscale=GV%H_to_mks)
      call hchksum(MLD_fast, 'mixed_layer_restrat: MLD fast', G%HI, haloshift=1, unscale=GV%H_to_mks)
    endif
    aFac = CS%MLE_MLD_decay_time2 / ( dt + CS%MLE_MLD_decay_time2 )
    bFac = dt / ( dt + CS%MLE_MLD_decay_time2 )
    do j=js-1,je+1 ; do i=is-1,ie+1
      ! Expression bFac*MLD_fast(i,j) + aFac*CS%MLD_filtered(i,j) is the time-filtered
      ! (running mean) of MLD. The max() allows the "running mean" to be reset
      ! instantly to a deeper MLD.
      CS%MLD_filtered_slow(i,j) = max( MLD_fast(i,j), bFac*MLD_fast(i,j) + aFac*CS%MLD_filtered_slow(i,j) )
      MLD_slow(i,j) = CS%MLD_filtered_slow(i,j)
    enddo ; enddo
  else
    do j=js-1,je+1 ; do i=is-1,ie+1
      MLD_slow(i,j) = MLD_fast(i,j)
    enddo ; enddo
  endif

  uDml(:) = 0.0 ; vDml(:) = 0.0
  uDml_slow(:) = 0.0 ; vDml_slow(:) = 0.0
  I4dt = 0.25 / dt
  g_Rho0 = GV%H_to_Z * GV%g_Earth / GV%Rho0
  h_neglect = GV%H_subroundoff
  if (CS%front_length>0.) then
    res_upscale = .true.
    do j=js-1,je+1 ; do i=is-1,ie+1
      mle_fl_2d(i,j) = CS%front_length
    enddo ; enddo
  elseif (CS%front_length == 0. .and. CS%fl_from_file) then
    res_upscale = .true.
    call time_interp_external(CS%sbc_fl, CS%Time, mle_fl_2d, turns=G%HI%turns, scale=US%m_to_L)
    call pass_var(mle_fl_2d, G%domain, halo=1)
    do j=js,je ; do i=is,ie
      if ((G%mask2dT(i,j) > 0.0) .and. (mle_fl_2d(i,j) < 0.0)) then
        write(mesg,'(" Time_interp negative MLE frontal-length scale of ",(1pe12.4)," at i,j = ",&
                  & 2(i3), "lon/lat = ",(1pe12.4)," E ", (1pe12.4), " N.")') &
                  mle_fl_2d(i,j), i, j, G%geoLonT(i,j), G%geoLatT(i,j)
        call MOM_error(FATAL, "MOM_mixed_layer_restrat mixedlayer_restrat_OM4: "//trim(mesg))
      endif
    enddo ; enddo
  else
    res_upscale = .false.
  endif

  p0(:) = 0.0
  EOSdom(:) = EOS_domain(G%HI, halo=1)
  !$OMP parallel default(shared) private(rho_ml,h_vel,u_star,absf,mom_mixrate,timescale, &
  !$OMP                                SpV_ml,SpV_int_fast,SpV_int_slow,Rml_int_fast,Rml_int_slow, &
  !$OMP                                line_is_empty,keep_going,res_scaling_fac, &
  !$OMP                                a,IhTot,b,Ihtot_slow,zpb,hAtVel,zpa,dh,lfront,I_LFront)         &
  !$OMP                        firstprivate(uDml,vDml,uDml_slow,vDml_slow)

  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    !$OMP do
    do j=js-1,je+1
      do i=is-1,ie+1
        htot_fast(i,j) = 0.0 ; Rml_int_fast(i) = 0.0
        htot_slow(i,j) = 0.0 ; Rml_int_slow(i) = 0.0
      enddo
      keep_going = .true.
      do k=1,nz
        do i=is-1,ie+1
          h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom_H),0.0)
        enddo
        if (keep_going) then
          if (CS%use_Stanley_ML) then
            call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p0, tv%varT(:,j,k), covTS, varS, &
                                   rho_ml(:), tv%eqn_of_state, EOSdom)
          else
            call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p0, rho_ml(:), tv%eqn_of_state, EOSdom)
          endif
          line_is_empty = .true.
          do i=is-1,ie+1
            if (htot_fast(i,j) < MLD_fast(i,j)) then
              dh = min( h(i,j,k), MLD_fast(i,j)-htot_fast(i,j) )
              Rml_int_fast(i) = Rml_int_fast(i) + dh*rho_ml(i)
              htot_fast(i,j) = htot_fast(i,j) + dh
              line_is_empty = .false.
            endif
            if (htot_slow(i,j) < MLD_slow(i,j)) then
              dh = min( h(i,j,k), MLD_slow(i,j)-htot_slow(i,j) )
              Rml_int_slow(i) = Rml_int_slow(i) + dh*rho_ml(i)
              htot_slow(i,j) = htot_slow(i,j) + dh
              line_is_empty = .false.
            endif
          enddo
          if (line_is_empty) keep_going=.false.
        endif
      enddo

      do i=is-1,ie+1
        Rml_av_fast(i,j) = -(g_Rho0*Rml_int_fast(i)) / (htot_fast(i,j) + h_neglect)
        Rml_av_slow(i,j) = -(g_Rho0*Rml_int_slow(i)) / (htot_slow(i,j) + h_neglect)
      enddo
    enddo
  else  ! This is only used in non-Boussinesq mode.
    !$OMP do
    do j=js-1,je+1
      do i=is-1,ie+1
        htot_fast(i,j) = 0.0 ; SpV_int_fast(i) = 0.0
        htot_slow(i,j) = 0.0 ; SpV_int_slow(i) = 0.0
      enddo
      keep_going = .true.
      do k=1,nz
        do i=is-1,ie+1
          h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom_H),0.0)
        enddo
        if (keep_going) then
          ! if (CS%use_Stanley_ML) then  ! This is not implemented yet in the EoS code.
          !   call calculate_spec_vol(tv%T(:,j,k), tv%S(:,j,k), p0, tv%varT(:,j,k), covTS, varS, &
          !                          rho_ml(:), tv%eqn_of_state, EOSdom)
          ! else
            call calculate_spec_vol(tv%T(:,j,k), tv%S(:,j,k), p0, SpV_ml, tv%eqn_of_state, EOSdom)
          ! endif
          line_is_empty = .true.
          do i=is-1,ie+1
            if (htot_fast(i,j) < MLD_fast(i,j)) then
              dh = min( h(i,j,k), MLD_fast(i,j)-htot_fast(i,j) )
              SpV_int_fast(i) = SpV_int_fast(i) + dh*SpV_ml(i)
              htot_fast(i,j) = htot_fast(i,j) + dh
              line_is_empty = .false.
            endif
            if (htot_slow(i,j) < MLD_slow(i,j)) then
              dh = min( h(i,j,k), MLD_slow(i,j)-htot_slow(i,j) )
              SpV_int_slow(i) = SpV_int_slow(i) + dh*SpV_ml(i)
              htot_slow(i,j) = htot_slow(i,j) + dh
              line_is_empty = .false.
            endif
          enddo
          if (line_is_empty) keep_going=.false.
        endif
      enddo

      ! Convert the vertically integrated specific volume into a positive variable with units of density.
      do i=is-1,ie+1
        Rml_av_fast(i,j) = (GV%H_to_RZ*GV%g_Earth * SpV_int_fast(i)) / (htot_fast(i,j) + h_neglect)
        Rml_av_slow(i,j) = (GV%H_to_RZ*GV%g_Earth * SpV_int_slow(i)) / (htot_slow(i,j) + h_neglect)
      enddo
    enddo
  endif

  if (CS%debug) then
    call hchksum(h, 'mixed_layer_restrat: h', G%HI, haloshift=1, unscale=GV%H_to_mks)
    call hchksum(U_star_2d, 'mixed_layer_restrat: u*', G%HI, haloshift=1, unscale=GV%H_to_m*US%s_to_T)
    call hchksum(MLD_fast, 'mixed_layer_restrat: MLD', G%HI, haloshift=1, unscale=GV%H_to_mks)
    call hchksum(Rml_av_fast, 'mixed_layer_restrat: rml', G%HI, haloshift=1, &
                 unscale=GV%m_to_H*US%L_T_to_m_s**2)
  endif

! TO DO:
!   1. Mixing extends below the mixing layer to the mixed layer.  Find it!
!   2. Add exponential tail to stream-function?

!   U - Component
  !$OMP do
  do j=js,je ; do I=is-1,ie
    u_star = max(CS%ustar_min, 0.5*(U_star_2d(i,j) + U_star_2d(i+1,j)))

    absf = 0.5*(abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))
    ! Compute I_LFront = 1 / (frontal length scale) [m-1]
    lfront = 0.5 * (mle_fl_2d(i,j) + mle_fl_2d(i+1,j))
    ! Adcroft reciprocal
    I_LFront = 0.0 ; if (lfront /= 0.0) I_LFront = 1.0/lfront
    ! If needed, res_scaling_fac = min( ds, L_d ) / l_f
    if (res_upscale) res_scaling_fac = &
          ( sqrt( 0.5 * ( (G%dxCu(I,j)**2) + (G%dyCu(I,j)**2) ) ) * I_LFront ) &
          * min( 1., 0.5*( VarMix%Rd_dx_h(i,j) + VarMix%Rd_dx_h(i+1,j) ) )

    ! peak ML visc: u_star * von_Karman * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2
    h_vel = 0.5*((htot_fast(i,j) + htot_fast(i+1,j)) + h_neglect)

    ! NOTE: growth_time changes answers on some systems, see below.
    ! timescale = growth_time(u_star, h_vel, absf, h_neglect, CS%vonKar, CS%Kv_restrat, CS%ml_restrat_coef)

    mom_mixrate = vonKar_x_pi2*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+h_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)
    timescale = timescale * CS%ml_restrat_coef

    if (res_upscale) timescale = timescale * res_scaling_fac
    uDml(I) = timescale * G%OBCmaskCu(I,j)*G%dyCu(I,j)*G%IdxCu(I,j) * &
        (Rml_av_fast(i+1,j)-Rml_av_fast(i,j)) * (h_vel**2)

    ! As above but using the slow filtered MLD
    h_vel = 0.5*((htot_slow(i,j) + htot_slow(i+1,j)) + h_neglect)

    ! NOTE: growth_time changes answers on some systems, see below.
    ! timescale = growth_time(u_star, h_vel, absf, h_neglect, CS%vonKar, CS%Kv_restrat, CS%ml_restrat_coef2)

    mom_mixrate = vonKar_x_pi2*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+h_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)
    timescale = timescale * CS%ml_restrat_coef2

    if (res_upscale) timescale = timescale * res_scaling_fac
    uDml_slow(I) = timescale * G%OBCmaskCu(I,j)*G%dyCu(I,j)*G%IdxCu(I,j) * &
        (Rml_av_slow(i+1,j)-Rml_av_slow(i,j)) * (h_vel**2)

    if (uDml(I) + uDml_slow(I) == 0.) then
      do k=1,nz ; uhml(I,j,k) = 0.0 ; enddo
    else
      IhTot = 2.0 / ((htot_fast(i,j) + htot_fast(i+1,j)) + h_neglect)
      IhTot_slow = 2.0 / ((htot_slow(i,j) + htot_slow(i+1,j)) + h_neglect)
      zpa = 0.0 ; zpb = 0.0
      ! a(k) relates the sublayer transport to uDml with a linear profile.
      ! The sum of a(k) through the mixed layers must be 0.
      do k=1,nz
        hAtVel = 0.5*(h(i,j,k) + h(i+1,j,k))
        a(k) = mu(zpa, CS%MLE_tail_dh)        ! mu(z/MLD) for upper interface
        zpa = zpa - (hAtVel * IhTot)          ! z/H for lower interface
        a(k) = a(k) - mu(zpa, CS%MLE_tail_dh) ! Transport profile
        ! Limit magnitude (uDml) if it would violate CFL
        if (a(k)*uDml(I) > 0.0) then
          if (a(k)*uDml(I) > h_avail(i,j,k)) uDml(I) = h_avail(i,j,k) / a(k)
        elseif (a(k)*uDml(I) < 0.0) then
          if (-a(k)*uDml(I) > h_avail(i+1,j,k)) uDml(I) = -h_avail(i+1,j,k) / a(k)
        endif
      enddo
      do k=1,nz
        ! Transport for slow-filtered MLD
        hAtVel = 0.5*(h(i,j,k) + h(i+1,j,k))
        b(k) = mu(zpb, CS%MLE_tail_dh)        ! mu(z/MLD) for upper interface
        zpb = zpb - (hAtVel * IhTot_slow)     ! z/H for lower interface
        b(k) = b(k) - mu(zpb, CS%MLE_tail_dh) ! Transport profile
        ! Limit magnitude (uDml_slow) if it would violate CFL when added to uDml
        if (b(k)*uDml_slow(I) > 0.0) then
          if (b(k)*uDml_slow(I) > h_avail(i,j,k) - a(k)*uDml(I)) &
             uDml_slow(I) = max( 0., h_avail(i,j,k) - a(k)*uDml(I) ) / b(k)
        elseif (b(k)*uDml_slow(I) < 0.0) then
          if (-b(k)*uDml_slow(I) > h_avail(i+1,j,k) + a(k)*uDml(I)) &
             uDml_slow(I) = -max( 0., h_avail(i+1,j,k) + a(k)*uDml(I) ) / b(k)
        endif
      enddo
      do k=1,nz
        uhml(I,j,k) = a(k)*uDml(I) + b(k)*uDml_slow(I)
        uhtr(I,j,k) = uhtr(I,j,k) + uhml(I,j,k)*dt
      enddo
    endif

    utimescale_diag(I,j) = timescale
    uDml_diag(I,j) = uDml(I)
  enddo ; enddo

!  V- component
  !$OMP do
  do J=js-1,je ; do i=is,ie
    u_star = max(CS%ustar_min, 0.5*(U_star_2d(i,j) + U_star_2d(i,j+1)))
    ! Compute I_LFront = 1 / (frontal length scale) [m-1]
    lfront = 0.5 * (mle_fl_2d(i,j) + mle_fl_2d(i,j+1))
    ! Adcroft reciprocal
    I_LFront = 0.0 ; if (lfront /= 0.0) I_LFront = 1.0/lfront
    absf = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
    ! If needed, res_scaling_fac = min( ds, L_d ) / l_f
    if (res_upscale) res_scaling_fac = &
          ( sqrt( 0.5 * ( (G%dxCv(i,J)**2) + (G%dyCv(i,J)**2) ) ) * I_LFront ) &
          * min( 1., 0.5*( VarMix%Rd_dx_h(i,j) + VarMix%Rd_dx_h(i,j+1) ) )

    ! peak ML visc: u_star * von_Karman * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2
    h_vel = 0.5*((htot_fast(i,j) + htot_fast(i,j+1)) + h_neglect)

    ! NOTE: growth_time changes answers on some systems, see below.
    ! timescale = growth_time(u_star, h_vel, absf, h_neglect, CS%vonKar, CS%Kv_restrat, CS%ml_restrat_coef)

    mom_mixrate = vonKar_x_pi2*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+h_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)
    timescale = timescale * CS%ml_restrat_coef

    if (res_upscale) timescale = timescale * res_scaling_fac
    vDml(i) = timescale * G%OBCmaskCv(i,J)*G%dxCv(i,J)*G%IdyCv(i,J) * &
        (Rml_av_fast(i,j+1)-Rml_av_fast(i,j)) * (h_vel**2)

    ! As above but using the slow filtered MLD
    h_vel = 0.5*((htot_slow(i,j) + htot_slow(i,j+1)) + h_neglect)

    ! NOTE: growth_time changes answers on some systems, see below.
    ! timescale = growth_time(u_star, h_vel, absf, h_neglect, CS%vonKar, CS%Kv_restrat, CS%ml_restrat_coef2)

    mom_mixrate = vonKar_x_pi2*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+h_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)
    timescale = timescale * CS%ml_restrat_coef2

    if (res_upscale) timescale = timescale * res_scaling_fac
    vDml_slow(i) = timescale * G%OBCmaskCv(i,J)*G%dxCv(i,J)*G%IdyCv(i,J) * &
        (Rml_av_slow(i,j+1)-Rml_av_slow(i,j)) * (h_vel**2)

    if (vDml(i) + vDml_slow(i) == 0.) then
      do k=1,nz ; vhml(i,J,k) = 0.0 ; enddo
    else
      IhTot = 2.0 / ((htot_fast(i,j) + htot_fast(i,j+1)) + h_neglect)
      IhTot_slow = 2.0 / ((htot_slow(i,j) + htot_slow(i,j+1)) + h_neglect)
      zpa = 0.0 ; zpb = 0.0
      ! a(k) relates the sublayer transport to vDml with a linear profile.
      ! The sum of a(k) through the mixed layers must be 0.
      do k=1,nz
        hAtVel = 0.5*(h(i,j,k) + h(i,j+1,k))
        a(k) = mu(zpa, CS%MLE_tail_dh)        ! mu(z/MLD) for upper interface
        zpa = zpa - (hAtVel * IhTot)          ! z/H for lower interface
        a(k) = a(k) - mu(zpa, CS%MLE_tail_dh) ! Transport profile
        ! Limit magnitude (vDml) if it would violate CFL
        if (a(k)*vDml(i) > 0.0) then
          if (a(k)*vDml(i) > h_avail(i,j,k)) vDml(i) = h_avail(i,j,k) / a(k)
        elseif (a(k)*vDml(i) < 0.0) then
          if (-a(k)*vDml(i) > h_avail(i,j+1,k)) vDml(i) = -h_avail(i,j+1,k) / a(k)
        endif
      enddo
      do k=1,nz
        ! Transport for slow-filtered MLD
        hAtVel = 0.5*(h(i,j,k) + h(i,j+1,k))
        b(k) = mu(zpb, CS%MLE_tail_dh)        ! mu(z/MLD) for upper interface
        zpb = zpb - (hAtVel * IhTot_slow)     ! z/H for lower interface
        b(k) = b(k) - mu(zpb, CS%MLE_tail_dh) ! Transport profile
        ! Limit magnitude (vDml_slow) if it would violate CFL when added to vDml
        if (b(k)*vDml_slow(i) > 0.0) then
          if (b(k)*vDml_slow(i) > h_avail(i,j,k) - a(k)*vDml(i)) &
             vDml_slow(i) = max( 0., h_avail(i,j,k) - a(k)*vDml(i) ) / b(k)
        elseif (b(k)*vDml_slow(i) < 0.0) then
          if (-b(k)*vDml_slow(i) > h_avail(i,j+1,k) + a(k)*vDml(i)) &
             vDml_slow(i) = -max( 0., h_avail(i,j+1,k) + a(k)*vDml(i) ) / b(k)
        endif
      enddo
      do k=1,nz
        vhml(i,J,k) = a(k)*vDml(i) + b(k)*vDml_slow(i)
        vhtr(i,J,k) = vhtr(i,J,k) + vhml(i,J,k)*dt
      enddo
    endif

    vtimescale_diag(i,J) = timescale
    vDml_diag(i,J) = vDml(i)
  enddo ; enddo

  !$OMP do
  do j=js,je ; do k=1,nz ; do i=is,ie
    h(i,j,k) = h(i,j,k) - dt*G%IareaT(i,j) * &
        ((uhml(I,j,k) - uhml(I-1,j,k)) + (vhml(i,J,k) - vhml(i,J-1,k)))
    if (h(i,j,k) < h_min) h(i,j,k) = h_min
  enddo ; enddo ; enddo
  !$OMP end parallel

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  if (CS%id_uhml > 0 .or. CS%id_vhml > 0) &
    ! Remapped uhml and vhml require east/north halo updates of h
    call pass_var(h, G%domain, To_West+To_South+Omit_Corners, halo=1)
  call diag_update_remap_grids(CS%diag)

  ! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_urestrat_time > 0) call post_data(CS%id_urestrat_time, utimescale_diag, CS%diag)
    if (CS%id_vrestrat_time > 0) call post_data(CS%id_vrestrat_time, vtimescale_diag, CS%diag)
    if (CS%id_uhml          > 0) call post_data(CS%id_uhml, uhml, CS%diag)
    if (CS%id_vhml          > 0) call post_data(CS%id_vhml, vhml, CS%diag)
    if (CS%id_BLD           > 0) call post_data(CS%id_BLD, MLD_fast, CS%diag)
    if (CS%id_MLD           > 0) call post_data(CS%id_MLD, MLD_slow, CS%diag)
    if (CS%id_Rml           > 0) call post_data(CS%id_Rml, Rml_av_fast, CS%diag)
    if (CS%id_uDml          > 0) call post_data(CS%id_uDml, uDml_diag, CS%diag)
    if (CS%id_vDml          > 0) call post_data(CS%id_vDml, vDml_diag, CS%diag)
    if (CS%id_mle_fl        > 0) call post_data(CS%id_mle_fl, mle_fl_2d, CS%diag)

    if (CS%id_uml > 0) then
      do j=js,je ; do I=is-1,ie
        h_vel = 0.5*((htot_fast(i,j) + htot_fast(i+1,j)) + h_neglect)
        uDml_diag(I,j) = uDml_diag(I,j) / (0.01*h_vel) * G%IdyCu(I,j) * (mu(0.,0.)-mu(-.01,0.))
      enddo ; enddo
      call post_data(CS%id_uml, uDml_diag, CS%diag)
    endif
    if (CS%id_vml > 0) then
      do J=js-1,je ; do i=is,ie
        h_vel = 0.5*((htot_fast(i,j) + htot_fast(i,j+1)) + h_neglect)
        vDml_diag(i,J) = vDml_diag(i,J) / (0.01*h_vel) * G%IdxCv(i,J) * (mu(0.,0.)-mu(-.01,0.))
      enddo ; enddo
      call post_data(CS%id_vml, vDml_diag, CS%diag)
    endif
  endif
  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  ! This needs to happen after the H update and before the next post_data.
  call diag_update_remap_grids(CS%diag)

end subroutine mixedlayer_restrat_OM4

!> Stream function shape as a function of non-dimensional position within mixed-layer [nondim]
real function mu(sigma, dh)
  real, intent(in) :: sigma !< Fractional position within mixed layer [nondim]
                            !! z=0 is surface, z=-1 is the bottom of the mixed layer
  real, intent(in) :: dh    !< Non-dimensional distance over which to extend stream
                            !! function to smooth transport at base [nondim]
  ! Local variables
  real :: xp            !< A linear function from mid-point of the mixed-layer
                        !! to the extended mixed-layer bottom [nondim]
  real :: bottop        !< A mask, 0 in upper half of mixed layer, 1 otherwise [nondim]
  real :: dd            !< A cubic(-ish) profile in lower half of extended mixed
                        !! layer to smooth out the parameterized transport [nondim]

  ! Lower order shape (not used), see eq 10 from FK08b.
  ! Apparently used in CM2G, see eq 14 of FK11.
  !mu = max(0., (1. - (2.*sigma + 1.)**2))

  ! Second order, in Rossby number, shape. See eq 21 from FK08a, eq 9 from FK08b, eq 5 FK11
  mu = max(0., (1. - (2.*sigma + 1.)**2) * (1. + (5./21.)*(2.*sigma + 1.)**2))

  ! -0.5    < sigma           : xp(sigma)=0      (upper half of mixed layer)
  ! -1.0+dh < sigma < -0.5    : xp(sigma)=linear (lower half +dh of mixed layer)
  !           sigma < -1.0+dh : xp(sigma)=1      (below mixed layer + dh)
  xp = max(0., min(1., (-sigma - 0.5)*2. / (1. + 2.*dh)))

  ! -0.5    < sigma           : dd(sigma)=1      (upper half of mixed layer)
  ! -1.0+dh < sigma < -0.5    : dd(sigma)=cubic  (lower half +dh of mixed layer)
  !           sigma < -1.0+dh : dd(sigma)=0      (below mixed layer + dh)
  dd = (max(1. - xp**2 * (3. - 2.*xp), 0.))**(1. + 2.*dh)

  ! -0.5    < sigma           : bottop(sigma)=0  (upper half of mixed layer)
  !           sigma < -0.5    : bottop(sigma)=1  (below upper half)
  bottop = 0.5*(1. - sign(1., sigma + 0.5))  ! =0 for sigma>-0.5, =1 for sigma<-0.5

  mu = max(mu, dd*bottop)  ! Combines original psi1 with tail
end function mu

!> Calculates a restratifying flow in the mixed layer, following the formulation
!! used in Bodner et al., 2023 (B22)
subroutine mixedlayer_restrat_Bodner(CS, G, GV, US, h, uhtr, vhtr, tv, forces, dt, BLD, h_MLD, bflux)
  ! Arguments
  type(mixedlayer_restrat_CS),                intent(inout) :: CS     !< Module control structure
  type(ocean_grid_type),                      intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                         intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  real, dimension(:,:),                       pointer       :: BLD    !< Active boundary layer depth provided by the
                                                                      !! PBL scheme [Z ~> m] (not H)
  real, dimension(:,:),                       pointer       :: h_MLD  !< Thickness of water within the
                                                                      !! active boundary layer depth provided by
                                                                      !! the PBL scheme [H ~> m or kg m-2]
  real, dimension(:,:),                       pointer       :: bflux  !< Surface buoyancy flux provided by the
                                                                      !! PBL scheme [Z2 T-3 ~> m2 s-3]
  ! Local variables
  real :: uhml(SZIB_(G),SZJ_(G),SZK_(GV)) ! zonal mixed layer transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vhml(SZI_(G),SZJB_(G),SZK_(GV)) ! merid mixed layer transport [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vol_dt_avail(SZI_(G),SZJ_(G),SZK_(GV)) ! The volume available for exchange out of each face of
                          ! each layer, divided by dt [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    little_h, &           ! "Little h" representing active mixing layer depth [H ~> m or kg m-2]
    big_H, &              ! "Big H" representing the mixed layer depth [H ~> m or kg m-2]
    mld, &                ! The mixed layer depth returned by detect_mld [H ~> m or kg m-2]
    htot, &               ! The sum of the thicknesses of layers in the mixed layer [H ~> m or kg m-2]
    buoy_av, &            ! g_Rho0 times the average mixed layer density or G_Earth
                          ! times the average specific volume [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2]
    wpup                  ! Turbulent vertical momentum [L H T-2 ~> m2 s-2 or kg m-1 s-2]
  real :: uDml_diag(SZIB_(G),SZJ_(G))  ! A 2D copy of uDml for diagnostics [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vDml_diag(SZI_(G),SZJB_(G))  ! A 2D copy of vDml for diagnostics [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: lf_bodner_diag(SZI_(G),SZJ_(G)) ! Front width as in Bodner et al., 2023 (B22), eq 24 [L ~> m]
  real :: U_star_2d(SZI_(G),SZJ_(G))   ! The wind friction velocity, calculated using the Boussinesq
                          ! reference density or the time-evolving surface density in non-Boussinesq
                          ! mode [Z T-1 ~> m s-1]
  real :: covTS(SZI_(G))  ! SGS TS covariance in Stanley param; currently 0 [C S ~> degC ppt]
  real :: varS(SZI_(G))   ! SGS S variance in Stanley param; currently 0 [S2 ~> ppt2]
  real :: dmu(SZK_(GV))   ! Change in mu(z) across layer k [nondim]
  real :: Rml_int(SZI_(G)) ! Potential density integrated through the mixed layer [R H ~> kg m-2 or kg2 m-5]
  real :: SpV_ml(SZI_(G)) ! Specific volume evaluated at the surface pressure [R-1 ~> m3 kg-1]
  real :: SpV_int(SZI_(G)) ! Specific volume integrated through the mixed layer [H R-1 ~> m4 kg-1 or m]
  real :: rho_ml(SZI_(G)) ! Potential density relative to the surface [R ~> kg m-3]
  real :: p0(SZI_(G))     ! A pressure of 0 [R L2 T-2 ~> Pa]
  real :: g_Rho0          ! G_Earth/Rho0 times a thickness conversion factor
                          ! [L2 H-1 T-2 R-1 ~> m4 s-2 kg-1 or m7 s-2 kg-2]
  real :: h_vel           ! htot interpolated onto velocity points [H ~> m or kg m-2]
  real :: w_star3         ! Cube of turbulent convective velocity [Z3 T-3 ~> m3 s-3]
  real :: u_star3         ! Cube of surface friction velocity [Z3 T-3 ~> m3 s-3]
  real :: r_wpup          ! reciprocal of vertical momentum flux [T2 L-1 H-1 ~> s2 m-2 or m s2 kg-1]
  real :: absf            ! absolute value of f, interpolated to velocity points [T-1 ~> s-1]
  real :: f_h             ! Coriolis parameter at h-points [T-1 ~> s-1]
  real :: f2_h            ! Coriolis parameter at h-points squared [T-2 ~> s-2]
  real :: absurdly_small_freq2 ! Frequency squared used to avoid division by 0 [T-2 ~> s-2]
  real :: grid_dsd        ! combination of grid scales [L2 ~> m2]
  real :: h_sml           ! "Little h", the active mixing depth with diurnal cycle removed [H ~> m or kg m-2]
  real :: h_big           ! "Big H", the mixed layer depth based on a time filtered "little h" [H ~> m or kg m-2]
  real :: grd_b           ! The vertically average gradient of buoyancy [L H-1 T-2 ~> s-2 or m-3 kg-1 s-2]
  real :: psi_mag         ! Magnitude of stream function [L2 H T-1 ~> m3 s-1 or kg s-1]
  real :: h_neglect       ! tiny thickness usually lost in roundoff so can be neglected [H ~> m or kg m-2]
  real :: I4dt            ! 1/(4 dt) [T-1 ~> s-1]
  real :: Ihtot,Ihtot_slow! Inverses of the total mixed layer thickness [H-1 ~> m-1 or m2 kg-1]
  real :: hAtVel          ! Thickness at the velocity points [H ~> m or kg m-2]
  real :: sigint          ! Fractional position within the mixed layer of the interface above a layer [nondim]
  real :: muzb            ! mu(z) at bottom of the layer [nondim]
  real :: muza            ! mu(z) at top of the layer [nondim]
  real :: dh              ! Portion of the layer thickness that is in the mixed layer [H ~> m or kg m-2]
  real :: res_scaling_fac ! The resolution-dependent scaling factor [nondim]
  real :: Z3_T3_to_m3_s3  ! Conversion factors to undo scaling and permit terms to be raised to a
                          ! fractional power [T3 m3 Z-3 s-3 ~> 1]
  real :: m2_s2_to_Z2_T2  ! Conversion factors to restore scaling after a term is raised to a
                          ! fractional power [Z2 s2 T-2 m-2 ~> 1]
  real, parameter :: two_thirds = 2./3.  ! [nondim]
  logical :: line_is_empty, keep_going
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  I4dt = 0.25 / dt
  g_Rho0 = GV%H_to_Z * GV%g_Earth / GV%Rho0
  h_neglect = GV%H_subroundoff

  covTS(:) = 0.0 ! Might be in tv% in the future. Not implemented for the time being.
  varS(:) = 0.0  ! Ditto.

 ! This value is roughly (pi / (the age of the universe) )^2.
  absurdly_small_freq2 = 1e-34*US%T_to_s**2

  if (.not.associated(tv%eqn_of_state)) call MOM_error(FATAL, "mixedlayer_restrat_Bodner: "// &
         "An equation of state must be used with this module.")
  if (CS%MLE_use_PBL_MLD) then
    if (CS%MLE_density_diff > 0.) call MOM_error(FATAL, "mixedlayer_restrat_Bodner: "// &
           "MLE_density_diff is +ve and should not be in mixedlayer_restrat_Bodner.")
    if (.not.associated(bflux)) call MOM_error(FATAL, "mixedlayer_restrat_Bodner: "// &
           "Surface buoyancy flux was not associated.")
  else
    if (.not.CS%Bodner_detect_MLD) call MOM_error(FATAL, "mixedlayer_restrat_Bodner: "// &
           "To use the Bodner et al., 2023, MLE parameterization, either MLE_USE_PBL_MLD or "// &
           "Bodner_detect_MLD must be True.")
  endif

  if (associated(bflux)) &
    call pass_var(bflux, G%domain, halo=1)

  ! Extract the friction velocity from the forcing type.
  call find_ustar(forces, tv, U_star_2d, G, GV, US, halo=1)

  if (CS%debug) then
    call hchksum(h,'mixed_Bodner: h', G%HI, haloshift=1, unscale=GV%H_to_mks)
    call hchksum(BLD, 'mle_Bodner: BLD', G%HI, haloshift=1, unscale=US%Z_to_m)
    call hchksum(h_MLD, 'mle_Bodner: h_MLD', G%HI, haloshift=1, unscale=GV%H_to_mks)
    if (associated(bflux)) &
      call hchksum(bflux, 'mle_Bodner: bflux', G%HI, haloshift=1, unscale=US%Z_to_m**2*US%s_to_T**3)
    call hchksum(U_star_2d, 'mle_Bodner: u*', G%HI, haloshift=1, unscale=US%Z_to_m*US%s_to_T)
    call hchksum(CS%MLD_filtered, 'mle_Bodner: MLD_filtered 1', &
                 G%HI, haloshift=1, unscale=GV%H_to_mks)
    call hchksum(CS%MLD_filtered_slow,'mle_Bodner: MLD_filtered_slow 1', &
                 G%HI, haloshift=1, unscale=GV%H_to_mks)
  endif

  ! Apply time filter to h_MLD (to remove diurnal cycle) to obtain "little h".
  ! "little h" is representative of the active mixing layer depth, used in B22 formula (eq 27).
  do j=js-1,je+1 ; do i=is-1,ie+1
    little_h(i,j) = rmean2ts(h_MLD(i,j), CS%MLD_filtered(i,j), &
                             CS%BLD_growing_Tfilt, CS%BLD_decaying_Tfilt, dt)
    CS%MLD_filtered(i,j) = little_h(i,j)
  enddo ; enddo

  ! Calculate "big H", representative of the mixed layer depth, used in B22 formula (eq 27).
  if (CS%MLD_grid) then
    do j=js-1,je+1 ; do i=is-1,ie+1
      big_H(i,j) = rmean2ts(little_h(i,j), CS%MLD_filtered_slow(i,j), &
                            CS%MLD_growing_Tfilt, CS%MLD_Tfilt_space(i,j), dt)
    enddo ; enddo
  elseif (CS%Bodner_detect_MLD) then
    call detect_mld(h, tv, MLD, G, GV, CS)
    do j=js-1,je+1 ; do i=is-1,ie+1
      big_H(i,j) = rmean2ts(MLD(i,j), CS%MLD_filtered_slow(i,j), &
                            CS%MLD_growing_Tfilt, CS%MLD_decaying_Tfilt, dt)
    enddo ; enddo
  else
    do j=js-1,je+1 ; do i=is-1,ie+1
      big_H(i,j) = rmean2ts(little_h(i,j), CS%MLD_filtered_slow(i,j), &
                            CS%MLD_growing_Tfilt, CS%MLD_decaying_Tfilt, dt)
    enddo ; enddo
  endif
  do j=js-1,je+1 ; do i=is-1,ie+1
    CS%MLD_filtered_slow(i,j) = big_H(i,j)
  enddo ; enddo

  ! Estimate w'u' at h-points, with a floor to avoid division by zero later.
  if (allocated(tv%SpV_avg) .and. .not.(GV%Boussinesq .or. GV%semi_Boussinesq)) then
    do j=js-1,je+1 ; do i=is-1,ie+1
      ! This expression differs by a factor of 1. / (Rho_0 * SpV_avg) compared with the other
      ! expressions below, and it is invariant to the value of Rho_0 in non-Boussinesq mode.
      wpup(i,j) = max((cuberoot( CS%mstar * U_star_2d(i,j)**3 + &
                                 CS%nstar * max(0., -bflux(i,j)) * BLD(i,j) ))**2, CS%min_wstar2) &
                  * (US%Z_to_L * GV%RZ_to_H / tv%SpV_avg(i,j,1))
                ! The final line above converts from [Z2 T-2 ~> m2 s-2] to [L H T-2 ~> m2 s-2 or Pa].
                ! Some rescaling factors and the division by specific volume compensating for other
                ! factors that are in find_ustar_mech, and others effectively converting the wind
                ! stresses from [R L Z T-2 ~> Pa] to [L H T-2 ~> m2 s-2 or Pa].  The rescaling factors
                ! and density being applied to the buoyancy flux are not so neatly explained because
                ! fractional powers cancel out or combine with terms in the definitions of BLD and
                ! bflux (such as SpV_avg**-2/3 combining with other terms in bflux to give the thermal
                ! expansion coefficient) and because the specific volume does vary within the mixed layer.
    enddo ; enddo
  elseif (CS%answer_date < 20240201) then
    Z3_T3_to_m3_s3 = (US%Z_to_m * US%s_to_T)**3
    m2_s2_to_Z2_T2 = (US%m_to_Z * US%T_to_s)**2
    do j=js-1,je+1 ; do i=is-1,ie+1
      w_star3 = max(0., -bflux(i,j)) * BLD(i,j)    ! In [Z3 T-3 ~> m3 s-3]
      u_star3 = U_star_2d(i,j)**3                  ! In [Z3 T-3 ~> m3 s-3]
      wpup(i,j) = max(m2_s2_to_Z2_T2 * (Z3_T3_to_m3_s3 * ( CS%mstar * u_star3 + CS%nstar * w_star3 ) )**two_thirds, &
          CS%min_wstar2) * US%Z_to_L * GV%Z_to_H ! In [L H T-2 ~> m2 s-2 or kg m-1 s-2]
    enddo ; enddo
  else
    do j=js-1,je+1 ; do i=is-1,ie+1
      w_star3 = max(0., -bflux(i,j)) * BLD(i,j)    ! In [Z3 T-3 ~> m3 s-3]
      wpup(i,j) = max( (cuberoot(CS%mstar * U_star_2d(i,j)**3 + CS%nstar * w_star3))**2, CS%min_wstar2 ) &
          * US%Z_to_L * GV%Z_to_H ! In [L H T-2 ~> m2 s-2 or kg m-1 s-2]
    enddo ; enddo
  endif

  ! We filter w'u' with the same time scales used for "little h"
  do j=js-1,je+1 ; do i=is-1,ie+1
    wpup(i,j) = rmean2ts(wpup(i,j), CS%wpup_filtered(i,j), &
                         CS%BLD_growing_Tfilt, CS%BLD_decaying_Tfilt, dt)
    CS%wpup_filtered(i,j) = wpup(i,j)
  enddo ; enddo

  if (CS%id_lfbod > 0) then
    do j=js-1,je+1 ; do i=is-1,ie+1
      ! Calculate front length used in B22 formula (eq 24).
      w_star3 = max(0., -bflux(i,j)) * BLD(i,j)
      u_star3 = U_star_2d(i,j)**3

      ! Include an absurdly_small_freq2 to prevent division by zero.
      f_h = 0.25 * ((G%CoriolisBu(I,J)  + G%CoriolisBu(I-1,J-1)) &
          + (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J-1)))
      f2_h = max(f_h**2, absurdly_small_freq2)

      lf_bodner_diag(i,j) = &
          0.25 * cuberoot(CS%mstar * u_star3 + CS%nstar * w_star3)**2 &
            / (f2_h * max(little_h(i,j), GV%Angstrom_H))
    enddo ; enddo

    ! Rescale from [Z2 H-1 to L]
    if (allocated(tv%SpV_avg) .and. .not.(GV%Boussinesq .or. GV%semi_Boussinesq)) then
      do j=js-1,je+1 ; do i=is-1,ie+1
        lf_bodner_diag(i,j) = lf_bodner_diag(i,j) &
            * (US%Z_to_L * GV%RZ_to_H / tv%SpV_avg(i,j,1))
      enddo ; enddo
    else
      do j=js-1,je+1 ; do i=is-1,ie+1
        lf_bodner_diag(i,j) = lf_bodner_diag(i,j) * US%Z_to_L * GV%Z_to_H
      enddo ; enddo
    endif
  endif

  if (CS%debug) then
    call hchksum(little_h,'mle_Bodner: little_h', G%HI, haloshift=1, unscale=GV%H_to_mks)
    call hchksum(big_H,'mle_Bodner: big_H', G%HI, haloshift=1, unscale=GV%H_to_mks)
    call hchksum(CS%MLD_filtered,'mle_Bodner: MLD_filtered 2', &
                 G%HI, haloshift=1, unscale=GV%H_to_mks)
    call hchksum(CS%MLD_filtered_slow,'mle_Bodner: MLD_filtered_slow 2', &
                 G%HI, haloshift=1, unscale=GV%H_to_mks)
    call hchksum(wpup,'mle_Bodner: wpup', G%HI, haloshift=1, unscale=US%L_to_m*GV%H_to_mks*US%s_to_T**2)
  endif

  ! Calculate the average density in the "mixed layer".
  ! Notice we use p=0 (sigma_0) since horizontal differences of vertical averages of
  ! in-situ density would contain the MLD gradient (through the pressure dependence).
  p0(:) = 0.0
  EOSdom(:) = EOS_domain(G%HI, halo=1)
  !$OMP parallel &
  !$OMP default(shared) &
  !$OMP private(i, j, k, keep_going, line_is_empty, dh, &
  !$OMP   grid_dsd, absf, h_sml, h_big, grd_b, r_wpup, psi_mag, IhTot, &
  !$OMP   sigint, muzb, muza, hAtVel, Rml_int, SpV_int)

  !$OMP do
  do j=js-1,je+1
    rho_ml(:) = 0.0 ; SpV_ml(:) = 0.0
    do i=is-1,ie+1
      htot(i,j) = 0.0 ; Rml_int(i) = 0.0 ; SpV_int(i) = 0.0
    enddo
    keep_going = .true.
    do k=1,nz
      do i=is-1,ie+1
        vol_dt_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom_H),0.0)
      enddo
      if (keep_going) then
        if (GV%Boussinesq .or. GV%semi_Boussinesq) then
          if (CS%use_Stanley_ML) then
            call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p0, tv%varT(:,j,k), covTS, varS, &
                                   rho_ml, tv%eqn_of_state, EOSdom)
          else
            call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p0, rho_ml, tv%eqn_of_state, EOSdom)
          endif
        else
          call calculate_spec_vol(tv%T(:,j,k), tv%S(:,j,k), p0, SpV_ml, tv%eqn_of_state, EOSdom)
        endif
        line_is_empty = .true.
        do i=is-1,ie+1
          if (htot(i,j) < big_H(i,j)) then
            dh = min( h(i,j,k), big_H(i,j) - htot(i,j) )
            Rml_int(i) = Rml_int(i) + dh*rho_ml(i) ! Rml_int has units of [R H ~> kg m-2]
            SpV_int(i) = SpV_int(i) + dh*SpV_ml(i) ! SpV_int has units of [H R-1 ~> m4 kg-1 or m]
            htot(i,j) = htot(i,j) + dh
            line_is_empty = .false.
          endif
        enddo
        if (line_is_empty) keep_going=.false.
      endif
    enddo

    if (GV%Boussinesq .or. GV%semi_Boussinesq) then
      do i=is-1,ie+1
        ! Buoy_av has units (L2 H-1 T-2 R-1) * (R H) * H-1 = L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2
        buoy_av(i,j) = -( g_Rho0 * Rml_int(i) ) / (htot(i,j) + h_neglect)
      enddo
    else
      do i=is-1,ie+1
        ! Buoy_av has units (R L2 H-1 T-2) * (R-1 H) * H-1 = L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2
        buoy_av(i,j) = (GV%H_to_RZ*GV%g_Earth * SpV_int(i)) / (htot(i,j) + h_neglect)
      enddo
    endif
  enddo

  if (CS%debug) then
    call hchksum(htot,'mle_Bodner: htot', G%HI, haloshift=1, unscale=GV%H_to_mks)
    call hchksum(vol_dt_avail,'mle_Bodner: vol_dt_avail', G%HI, haloshift=1, &
                 unscale=US%L_to_m**2*GV%H_to_mks*US%s_to_T)
    call hchksum(buoy_av,'mle_Bodner: buoy_av', G%HI, haloshift=1, unscale=GV%m_to_H*US%L_T_to_m_s**2)
  endif

  ! U - Component
  !$OMP do
  do j=js,je ; do I=is-1,ie
    if (G%OBCmaskCu(I,j) > 0.) then
      grid_dsd = sqrt(0.5*( G%dxCu(I,j)**2 + G%dyCu(I,j)**2 )) * G%dyCu(I,j) ! L2 ~> m2
      absf = 0.5*(abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))  ! T-1 ~> s-1
      h_sml = 0.5*( little_h(i,j) + little_h(i+1,j) )                 ! H ~> m or kg m-3
      h_big = 0.5*( big_H(i,j) + big_H(i+1,j) )                       ! H ~> m or kg m-3
      grd_b = ( buoy_av(i+1,j) - buoy_av(i,j) ) * G%IdxCu(I,j)        ! L H-1 T-2 ~> s-2 or m3 kg-1 s-2
      r_wpup = 2. / ( wpup(i,j) + wpup(i+1,j) )                       ! T2 L-1 H-1 ~> s2 m-2 or m s2 kg-1
      psi_mag = ( ( ( (0.5*(CS%Cr_space(i,j) + CS%Cr_space(i+1,j))) * grid_dsd ) & ! L2 H T-1 ~> m3 s-1 or kg s-1
                  * ( absf * h_sml ) ) * ( ( h_big**2 ) * grd_b ) ) * r_wpup
    else  ! There is no flux on land and no gradient at open boundary points.
      psi_mag = 0.0
    endif

    IhTot = 2.0 / ((htot(i,j) + htot(i+1,j)) + h_neglect) ! [H-1]
    sigint = 0.0
    muzb = 0.0 ! This will be the first value of muza = mu(z=0)
    do k=1,nz
      muza = muzb                           ! mu(z/MLD) for upper interface [nondim]
      hAtVel = 0.5*(h(i,j,k) + h(i+1,j,k))  ! Thickness at velocity point [H]
      sigint = sigint - (hAtVel * IhTot)    ! z/H for lower interface [nondim]
      muzb = mu(sigint, CS%MLE_tail_dh)     ! mu(z/MLD) for lower interface [nondim]
      dmu(k) = muza - muzb                  ! Change in mu(z) across layer [nondim]
      ! dmu(k)*psi_mag is the transport in this layer [L2 H T-1 ~> m3 s-1]
      ! Limit magnitude (psi_mag) if it would violate CFL
      if (dmu(k)*psi_mag > 0.0) then
        if (dmu(k)*psi_mag > vol_dt_avail(i,j,k)) psi_mag = vol_dt_avail(i,j,k) / dmu(k)
      elseif (dmu(k)*psi_mag < 0.0) then
        if (-dmu(k)*psi_mag > vol_dt_avail(i+1,j,k)) psi_mag = -vol_dt_avail(i+1,j,k) / dmu(k)
      endif
    enddo ! These loops cannot be fused because psi_mag applies to the whole column
    do k=1,nz
      uhml(I,j,k) = dmu(k) * psi_mag ! [ L2 H T-1 ]
      uhtr(I,j,k) = uhtr(I,j,k) + uhml(I,j,k) * dt ! [ L2 H ]
    enddo

    uDml_diag(I,j) = psi_mag
  enddo ; enddo

  ! V- component
  !$OMP do
  do J=js-1,je ; do i=is,ie
    if (G%OBCmaskCv(i,J) > 0.) then
      grid_dsd = sqrt(0.5*( G%dxCv(i,J)**2 + G%dyCv(i,J)**2 )) * G%dxCv(i,J) ! L2 ~> m2
      absf = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))  ! T-1 ~> s-1
      h_sml = 0.5*( little_h(i,j) + little_h(i,j+1) )                 ! H ~> m or kg m-3
      h_big = 0.5*( big_H(i,j) + big_H(i,j+1) )                       ! H ~> m or kg m-3
      grd_b = ( buoy_av(i,j+1) - buoy_av(i,j) ) * G%IdyCv(I,j)        ! L H-1 T-2 ~> s-2 or m3 kg-1 s-2
      r_wpup = 2. / ( wpup(i,j) + wpup(i,j+1) )                       ! T2 L-1 H-1 ~> s2 m-2 or m s2 kg-1
      psi_mag = ( ( ( (0.5*(CS%Cr_space(i,j) + CS%Cr_space(i,j+1))) * grid_dsd ) & ! L2 H T-1 ~> m3 s-1 or kg s-1
                  * ( absf * h_sml ) ) * ( ( h_big**2 ) * grd_b ) ) * r_wpup
    else  ! There is no flux on land and no gradient at open boundary points.
      psi_mag = 0.0
    endif

    IhTot = 2.0 / ((htot(i,j) + htot(i,j+1)) + h_neglect) ! [H-1]
    sigint = 0.0
    muzb = 0.0 ! This will be the first value of muza = mu(z=0)
    do k=1,nz
      muza = muzb                           ! mu(z/MLD) for upper interface [nondim]
      hAtVel = 0.5*(h(i,j,k) + h(i,j+1,k))  ! Thickness at velocity point [H]
      sigint = sigint - (hAtVel * IhTot)    ! z/H for lower interface [nondim]
      muzb = mu(sigint, CS%MLE_tail_dh)     ! mu(z/MLD) for lower interface [nondim]
      dmu(k) = muza - muzb                  ! Change in mu(z) across layer [nondim]
      ! dmu(k)*psi_mag is the transport in this layer [L2 H T-1 ~> m3 s-1]
      ! Limit magnitude (psi_mag) if it would violate CFL
      if (dmu(k)*psi_mag > 0.0) then
        if (dmu(k)*psi_mag > vol_dt_avail(i,j,k)) psi_mag = vol_dt_avail(i,j,k) / dmu(k)
      elseif (dmu(k)*psi_mag < 0.0) then
        if (-dmu(k)*psi_mag > vol_dt_avail(i,j+1,k)) psi_mag = -vol_dt_avail(i,j+1,k) / dmu(k)
      endif
    enddo ! These loops cannot be fused because psi_mag applies to the whole column
    do k=1,nz
      vhml(i,J,k) = dmu(k) * psi_mag ! [ L2 H T-1 ]
      vhtr(i,J,k) = vhtr(i,J,k) + vhml(i,J,k) * dt ! [ L2 H ]
    enddo

    vDml_diag(i,J) = psi_mag
  enddo ; enddo

  !$OMP do
  do j=js,je ; do k=1,nz ; do i=is,ie
    h(i,j,k) = h(i,j,k) - dt*G%IareaT(i,j) * &
        ((uhml(I,j,k) - uhml(I-1,j,k)) + (vhml(i,J,k) - vhml(i,J-1,k)))
  enddo ; enddo ; enddo
  !$OMP end parallel

  if (CS%id_uhml > 0 .or. CS%id_vhml > 0) &
    ! Remapped uhml and vhml require east/north halo updates of h
    call pass_var(h, G%domain, To_West+To_South+Omit_Corners, halo=1)
  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  call diag_update_remap_grids(CS%diag)

  ! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_ustar > 0) call post_data(CS%id_ustar, U_star_2d, CS%diag)
    if (CS%id_bflux > 0) call post_data(CS%id_bflux, bflux, CS%diag)
    if (CS%id_wpup  > 0) call post_data(CS%id_wpup, wpup, CS%diag)
    if (CS%id_Rml   > 0) call post_data(CS%id_Rml, buoy_av, CS%diag)
    if (CS%id_BLD   > 0) call post_data(CS%id_BLD, little_h, CS%diag)
    if (CS%id_MLD   > 0) call post_data(CS%id_MLD, big_H, CS%diag)
    if (CS%id_uhml  > 0) call post_data(CS%id_uhml, uhml, CS%diag)
    if (CS%id_vhml  > 0) call post_data(CS%id_vhml, vhml, CS%diag)
    if (CS%id_uDml  > 0) call post_data(CS%id_uDml, uDml_diag, CS%diag)
    if (CS%id_vDml  > 0) call post_data(CS%id_vDml, vDml_diag, CS%diag)
    if (CS%id_lfbod  > 0) call post_data(CS%id_lfbod, lf_bodner_diag, CS%diag)

    if (CS%id_uml > 0) then
      do j=js,je ; do I=is-1,ie
        h_vel = 0.5*((htot(i,j) + htot(i+1,j)) + h_neglect)
        uDml_diag(I,j) = uDml_diag(I,j) / (0.01*h_vel) * G%IdyCu(I,j) * (mu(0.,0.)-mu(-.01,0.))
      enddo ; enddo
      call post_data(CS%id_uml, uDml_diag, CS%diag)
    endif
    if (CS%id_vml > 0) then
      do J=js-1,je ; do i=is,ie
        h_vel = 0.5*((htot(i,j) + htot(i,j+1)) + h_neglect)
        vDml_diag(i,J) = vDml_diag(i,J) / (0.01*h_vel) * G%IdxCv(i,J) * (mu(0.,0.)-mu(-.01,0.))
      enddo ; enddo
      call post_data(CS%id_vml, vDml_diag, CS%diag)
    endif
  endif

end subroutine mixedlayer_restrat_Bodner

!> Two time-scale running mean [units of "signal" and "filtered"]
!!
!! If signal > filtered, returns running-mean with time scale "tau_growing".
!! If signal <= filtered, returns running-mean with time scale "tau_decaying".
!!
!! The running mean of \f$ s \f$ with time scale "of \f$ \tau \f$ is:
!! \f[
!!   \bar{s} <- ( \Delta t * s + \tau * \bar{s} ) / ( \Delta t + \tau )
!! \f]
!!
!! Note that if \f$ tau=0 \f$, then the running mean equals the signal. Thus,
!! rmean2ts with tau_growing=0 recovers the "resetting running mean" used in OM4.
real elemental function rmean2ts(signal, filtered, tau_growing, tau_decaying, dt)
  ! Arguments
  real, intent(in) :: signal       ! Unfiltered signal [arbitrary units]
  real, intent(in) :: filtered     ! Current value of running mean [arbitrary units]
  real, intent(in) :: tau_growing  ! Time scale for growing signal [T ~> s]
  real, intent(in) :: tau_decaying ! Time scale for decaying signal [T ~> s]
  real, intent(in) :: dt           ! Time step [T ~> s]
  ! Local variables
  real :: afac, bfac ! Non-dimensional fractional weights [nondim]
  real :: rt ! Reciprocal time scale [T-1 ~> s-1]

  if (signal>=filtered) then
    rt = 1.0 / ( dt + tau_growing )
    aFac = tau_growing * rt
    bFac = 1. - aFac
  else
    rt = 1.0 / ( dt + tau_decaying )
    aFac = tau_decaying * rt
    bFac = 1. - aFac
  endif

  rmean2ts = aFac * filtered + bFac * signal

end function rmean2ts

!> Calculates a restratifying flow assuming a 2-layer bulk mixed layer.
subroutine mixedlayer_restrat_BML(h, uhtr, vhtr, tv, forces, dt, G, GV, US, CS)
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr   !< Accumulated zonal mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhtr   !< Accumulated meridional mass flux
                                                                      !!   [H L2 ~> m3 or kg]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamic variables structure
  type(mech_forcing),                         intent(in)    :: forces !< A structure with the driving mechanical forces
  real,                                       intent(in)    :: dt     !< Time increment [T ~> s]
  type(mixedlayer_restrat_CS),                intent(inout) :: CS     !< Module control structure

  ! Local variables
  real :: uhml(SZIB_(G),SZJ_(G),SZK_(GV)) ! Restratifying zonal thickness transports [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vhml(SZI_(G),SZJB_(G),SZK_(GV)) ! Restratifying meridional thickness transports [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    h_avail               ! The volume available for diffusion out of each face of each
                          ! sublayer of the mixed layer, divided by dt [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    U_star_2d, &          ! The wind friction velocity in thickness-based units, calculated using
                          ! the Boussinesq reference density or the time-evolving surface density
                          ! in non-Boussinesq mode [H T-1 ~> m s-1 or kg m-2 s-1]
    htot, &               ! The sum of the thicknesses of layers in the mixed layer [H ~> m or kg m-2]
    Rml_av                ! g_Rho0 times the average mixed layer density or negative G_Earth
                          ! times the average specific volume [L2 H-1 T-2 ~> m s-2 or m4 kg-1 s-2]
  real :: g_Rho0          ! G_Earth/Rho0 times a thickness conversion factor
                          ! [L2 H-1 T-2 R-1 ~> m4 s-2 kg-1 or m7 s-2 kg-2]
  real :: Rho_ml(SZI_(G)) ! Potential density relative to the surface [R ~> kg m-3]
  real :: rho_int(SZI_(G)) ! The integral of density over the mixed layer depth [R H ~> kg m-2 or kg2 m-3]
  real :: SpV_ml(SZI_(G)) ! Specific volume evaluated at the surface pressure [R-1 ~> m3 kg-1]
  real :: SpV_int(SZI_(G)) ! Specific volume integrated through the surface layer [H R-1 ~> m4 kg-1 or m]
  real :: p0(SZI_(G))     ! A pressure of 0 [R L2 T-2 ~> Pa]

  real :: h_vel           ! htot interpolated onto velocity points [H ~> m or kg m-2]
  real :: absf            ! absolute value of f, interpolated to velocity points [T-1 ~> s-1]
  real :: u_star          ! surface friction velocity, interpolated to velocity points and recast into
                          ! thickness-based units [H T-1 ~> m s-1 or kg m-2 s-1].
  real :: vonKar_x_pi2    ! A scaling constant that is approximately the von Karman constant times
                          ! pi squared [nondim]
  real :: mom_mixrate     ! rate at which momentum is homogenized within mixed layer [T-1 ~> s-1]
  real :: timescale       ! mixing growth timescale [T ~> s]
  real :: h_min           ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  real :: h_neglect       ! tiny thickness usually lost in roundoff and can be neglected [H ~> m or kg m-2]
  real :: I4dt            ! 1/(4 dt) [T-1 ~> s-1]
  real :: I2htot          ! Twice the total mixed layer thickness at velocity points [H ~> m or kg m-2]
  real :: z_topx2         ! depth of the top of a layer at velocity points [H ~> m or kg m-2]
  real :: hx2             ! layer thickness at velocity points [H ~> m or kg m-2]
  real :: a(SZK_(GV))     ! A non-dimensional value relating the overall flux magnitudes (uDml & vDml)
                          ! to the realized flux in a layer [nondim].  The vertical sum of a()
                          ! through the pieces of the mixed layer must be 0.
  real :: uDml(SZIB_(G))  ! Zonal volume fluxes in the upper half of the mixed layer [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vDml(SZI_(G))   ! Meridional volume fluxes in the upper half of the mixed layer [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: utimescale_diag(SZIB_(G),SZJ_(G)) ! Zonal restratification timescale [T ~> s], stored for diagnostics.
  real :: vtimescale_diag(SZI_(G),SZJB_(G)) ! Meridional restratification timescale [T ~> s], stored for diagnostics.
  real :: uDml_diag(SZIB_(G),SZJ_(G))  ! A 2D copy of uDml for diagnostics [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: vDml_diag(SZI_(G),SZJB_(G))  ! A 2D copy of vDml for diagnostics [H L2 T-1 ~> m3 s-1 or kg s-1]
  logical :: use_EOS    ! If true, density is calculated from T & S using an equation of state.

  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkml
  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nkml = GV%nkml

  if (.not. CS%initialized) call MOM_error(FATAL, "mixedlayer_restrat_BML: "// &
         "Module must be initialized before it is used.")

  if ((nkml<2) .or. (CS%ml_restrat_coef<=0.0)) return


  h_min = 0.5*GV%Angstrom_H ! This should be GV%Angstrom_H, but that value would change answers.
  uDml(:)    = 0.0 ; vDml(:) = 0.0
  I4dt       = 0.25 / dt
  g_Rho0     = GV%H_to_Z * GV%g_Earth / GV%Rho0
  vonKar_x_pi2 = CS%vonKar * 9.8696
  use_EOS    = associated(tv%eqn_of_state)
  h_neglect  = GV%H_subroundoff

  if (.not.use_EOS) call MOM_error(FATAL, "mixedlayer_restrat_BML: "// &
         "An equation of state must be used with this module.")

  if (CS%use_Stanley_ML) call MOM_error(FATAL, "mixedlayer_restrat_BML: "// &
         "The Stanley parameterization is not available with the BML.")

  ! Extract the friction velocity from the forcing type.
  call find_ustar(forces, tv, U_star_2d, G, GV, US, halo=1, H_T_units=.true.)

  ! Fix this later for nkml >= 3.

  p0(:) = 0.0
  EOSdom(:) = EOS_domain(G%HI, halo=1)
  !$OMP parallel default(shared) private(Rho_ml,rho_int,h_vel,u_star,absf,mom_mixrate,timescale, &
  !$OMP                                  SpV_ml,SpV_int,I2htot,z_topx2,hx2,a) &
  !$OMP                       firstprivate(uDml,vDml)

  if (GV%Boussinesq .or. GV%semi_Boussinesq) then
    !$OMP do
    do j=js-1,je+1
      do i=is-1,ie+1
        htot(i,j) = 0.0 ; rho_int(i) = 0.0
      enddo
      do k=1,nkml
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p0, Rho_ml(:), tv%eqn_of_state, EOSdom)
        do i=is-1,ie+1
          rho_int(i) = rho_int(i) + h(i,j,k)*Rho_ml(i)
          htot(i,j) = htot(i,j) + h(i,j,k)
          h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom_H),0.0)
        enddo
      enddo

      do i=is-1,ie+1
        Rml_av(i,j) = (g_Rho0*rho_int(i)) / (htot(i,j) + h_neglect)
      enddo
    enddo
  else  ! This is only used in non-Boussinesq mode.
    !$OMP do
    do j=js-1,je+1
      do i=is-1,ie+1
        htot(i,j) = 0.0 ; SpV_int(i) = 0.0
      enddo
      do k=1,nkml
        call calculate_spec_vol(tv%T(:,j,k), tv%S(:,j,k), p0, SpV_ml, tv%eqn_of_state, EOSdom)
        do i=is-1,ie+1
          SpV_int(i) = SpV_int(i) + h(i,j,k)*SpV_ml(i)  ! [H R-1 ~> m4 kg-1 or m]
          htot(i,j) = htot(i,j) + h(i,j,k)
          h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom_H),0.0)
        enddo
      enddo

      ! Convert the vertically integrated specific volume into a negative variable with units of density.
      do i=is-1,ie+1
        Rml_av(i,j) = (-GV%H_to_RZ*GV%g_Earth * SpV_int(i)) / (htot(i,j) + h_neglect)
      enddo
    enddo
  endif

! TO DO:
!   1. Mixing extends below the mixing layer to the mixed layer.  Find it!
!   2. Add exponential tail to stream-function?

!   U - Component
  !$OMP do
  do j=js,je ; do I=is-1,ie
    h_vel = 0.5*(htot(i,j) + htot(i+1,j))

    u_star = max(CS%ustar_min, 0.5*(U_star_2d(i,j) + U_star_2d(i+1,j)))

    absf = 0.5*(abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))

    ! NOTE: growth_time changes answers on some systems, see below.
    ! timescale = growth_time(u_star, h_vel, absf, h_neglect, CS%vonKar, CS%Kv_restrat, CS%ml_restrat_coef)

    ! peak ML visc: u_star * von_Karman * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2
    mom_mixrate = vonKar_x_pi2*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+h_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)

    timescale = timescale * CS%ml_restrat_coef
!      timescale = timescale*(2?)*(L_def/L_MLI) * min(EKE/MKE,1.0 + (G%dyCv(i,j)/L_def)**2)

    uDml(I) = timescale * G%OBCmaskCu(I,j)*G%dyCu(I,j)*G%IdxCu(I,j) * &
        (Rml_av(i+1,j)-Rml_av(i,j)) * (h_vel**2)

    if (uDml(I) == 0) then
      do k=1,nkml ; uhml(I,j,k) = 0.0 ; enddo
    else
      I2htot = 1.0 / (htot(i,j) + htot(i+1,j) + h_neglect)
      z_topx2 = 0.0
      ! a(k) relates the sublayer transport to uDml with a linear profile.
      ! The sum of a(k) through the mixed layers must be 0.
      do k=1,nkml
        hx2 = (h(i,j,k) + h(i+1,j,k) + h_neglect)
        a(k) = (hx2 * I2htot) * (2.0 - 4.0*(z_topx2+0.5*hx2)*I2htot)
        z_topx2 = z_topx2 + hx2
        if (a(k)*uDml(I) > 0.0) then
          if (a(k)*uDml(I) > h_avail(i,j,k)) uDml(I) = h_avail(i,j,k) / a(k)
        else
          if (-a(k)*uDml(I) > h_avail(i+1,j,k)) uDml(I) = -h_avail(i+1,j,k)/a(k)
        endif
      enddo
      do k=1,nkml
        uhml(I,j,k) = a(k)*uDml(I)
        uhtr(I,j,k) = uhtr(I,j,k) + uhml(I,j,k)*dt
      enddo
    endif

    uDml_diag(I,j) = uDml(I)
    utimescale_diag(I,j) = timescale
  enddo ; enddo

!  V- component
  !$OMP do
  do J=js-1,je ; do i=is,ie
    h_vel = 0.5*(htot(i,j) + htot(i,j+1))

    u_star = max(CS%ustar_min, 0.5*(U_star_2d(i,j) + U_star_2d(i,j+1)))

    absf = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))

    ! NOTE: growth_time changes answers on some systems, see below.
    ! timescale = growth_time(u_star, h_vel, absf, h_neglect, CS%vonKar, CS%Kv_restrat, CS%ml_restrat_coef)

    ! peak ML visc: u_star * von_Karman * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2
    mom_mixrate = vonKar_x_pi2*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+h_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)

    timescale = timescale * CS%ml_restrat_coef
!     timescale = timescale*(2?)*(L_def/L_MLI) * min(EKE/MKE,1.0 + (G%dyCv(i,j)/L_def)**2)

    vDml(i) = timescale * G%OBCmaskCv(i,J)*G%dxCv(i,J)*G%IdyCv(i,J) * &
        (Rml_av(i,j+1)-Rml_av(i,j)) * (h_vel**2)
    if (vDml(i) == 0) then
      do k=1,nkml ; vhml(i,J,k) = 0.0 ; enddo
    else
      I2htot = 1.0 / (htot(i,j) + htot(i,j+1) + h_neglect)
      z_topx2 = 0.0
      ! a(k) relates the sublayer transport to vDml with a linear profile.
      ! The sum of a(k) through the mixed layers must be 0.
      do k=1,nkml
        hx2 = (h(i,j,k) + h(i,j+1,k) + h_neglect)
        a(k) = (hx2 * I2htot) * (2.0 - 4.0*(z_topx2+0.5*hx2)*I2htot)
        z_topx2 = z_topx2 + hx2
        if (a(k)*vDml(i) > 0.0) then
          if (a(k)*vDml(i) > h_avail(i,j,k)) vDml(i) = h_avail(i,j,k) / a(k)
        else
          if (-a(k)*vDml(i) > h_avail(i,j+1,k)) vDml(i) = -h_avail(i,j+1,k)/a(k)
        endif
      enddo
      do k=1,nkml
        vhml(i,J,k) = a(k)*vDml(i)
        vhtr(i,J,k) = vhtr(i,J,k) + vhml(i,J,k)*dt
      enddo
    endif

    vtimescale_diag(i,J) = timescale
    vDml_diag(i,J) = vDml(i)
  enddo ; enddo

  !$OMP do
  do j=js,je ; do k=1,nkml ; do i=is,ie
    h(i,j,k) = h(i,j,k) - dt*G%IareaT(i,j) * &
        ((uhml(I,j,k) - uhml(I-1,j,k)) + (vhml(i,J,k) - vhml(i,J-1,k)))
    if (h(i,j,k) < h_min) h(i,j,k) = h_min
  enddo ; enddo ; enddo
  !$OMP end parallel

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  if (CS%id_uhml > 0 .or. CS%id_vhml > 0) &
    ! Remapped uhml and vhml require east/north halo updates of h
    call pass_var(h, G%domain, To_West+To_South+Omit_Corners, halo=1)
  call diag_update_remap_grids(CS%diag)

  ! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag) .and. &
    ((CS%id_urestrat_time > 0)  .or. (CS%id_vrestrat_time > 0))) then
    call post_data(CS%id_urestrat_time, utimescale_diag, CS%diag)
    call post_data(CS%id_vrestrat_time, vtimescale_diag, CS%diag)
  endif
  if (query_averaging_enabled(CS%diag) .and. &
      ((CS%id_uhml>0) .or. (CS%id_vhml>0))) then
    do k=nkml+1,nz
      do j=js,je ; do I=Isq,Ieq ; uhml(I,j,k) = 0.0 ; enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie ; vhml(i,J,k) = 0.0 ; enddo ; enddo
    enddo
    if (CS%id_uhml > 0) call post_data(CS%id_uhml, uhml,      CS%diag)
    if (CS%id_vhml > 0) call post_data(CS%id_vhml, vhml,      CS%diag)
    if (CS%id_MLD  > 0) call post_data(CS%id_MLD,  htot,      CS%diag)
    if (CS%id_Rml  > 0) call post_data(CS%id_Rml,  Rml_av,    CS%diag)
    if (CS%id_uDml > 0) call post_data(CS%id_uDml, uDml_diag, CS%diag)
    if (CS%id_vDml > 0) call post_data(CS%id_vDml, vDml_diag, CS%diag)
  endif

end subroutine mixedlayer_restrat_BML

!> Detects the mixed layer depth using a density difference criterion (MLE_DENSITY_DIFF)
subroutine detect_mld(h, tv, MLD_fast, G, GV, CS)
  type(mixedlayer_restrat_CS),                intent(inout) :: CS     !< Module control structure
  type(ocean_grid_type),                      intent(inout) :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),           intent(out)   :: MLD_fast !< detected mixed layer depth [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamic variables structure

  ! Local variables
  real, dimension(SZI_(G)) :: pRef_MLD ! A reference pressure for calculating the mixed layer
                                       ! densities [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G)) :: rhoSurf, deltaRhoAtKm1, deltaRhoAtK ! Densities and density differences [R ~> kg m-3]
  real, dimension(SZI_(G)) :: dK, dKm1 ! Depths of layer centers [H ~> m or kg m-2].
  real :: ddRho     ! A density difference [R ~> kg m-3]
  real :: aFac  ! A nondimensional ratio [nondim]
  real :: covTS(SZI_(G))  ! SGS TS covariance in Stanley param; currently 0 [C S ~> degC ppt]
  real :: varS(SZI_(G))   ! SGS S variance in Stanley param; currently 0 [S2 ~> ppt2]
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, nz

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke

  covTS(:) = 0.0 ! Might be in tv% in the future. Not implemented for the time being.
  varS(:) = 0.0  ! Ditto.

  !! TODO: use derivatives and mid-MLD pressure. Currently this is sigma-0. -AJA
  pRef_MLD(:) = 0.
  EOSdom(:) = EOS_domain(G%HI, halo=1)
  do j=js-1,je+1
    dK(:) = 0.5 * h(:,j,1) ! Depth of center of surface layer
    if (CS%use_Stanley_ML) then
      call calculate_density(tv%T(:,j,1), tv%S(:,j,1), pRef_MLD, tv%varT(:,j,1), covTS, varS, &
        rhoSurf, tv%eqn_of_state, EOSdom)
    else
      call calculate_density(tv%T(:,j,1), tv%S(:,j,1), pRef_MLD, rhoSurf, tv%eqn_of_state, EOSdom)
    endif
    deltaRhoAtK(:) = 0.
    MLD_fast(:,j) = 0.
    do k=2,nz
      dKm1(:) = dK(:) ! Depth of center of layer K-1
      dK(:) = dK(:) + 0.5 * ( h(:,j,k) + h(:,j,k-1) ) ! Depth of center of layer K
      ! Mixed-layer depth, using sigma-0 (surface reference pressure)
      deltaRhoAtKm1(:) = deltaRhoAtK(:) ! Store value from previous iteration of K
      if (CS%use_Stanley_ML) then
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pRef_MLD, tv%varT(:,j,k), covTS, varS, &
          deltaRhoAtK, tv%eqn_of_state, EOSdom)
      else
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pRef_MLD, deltaRhoAtK, tv%eqn_of_state, EOSdom)
      endif
      do i=is-1,ie+1
        deltaRhoAtK(i) = deltaRhoAtK(i) - rhoSurf(i) ! Density difference between layer K and surface
      enddo
      do i=is-1,ie+1
        ddRho = deltaRhoAtK(i) - deltaRhoAtKm1(i)
        if ((MLD_fast(i,j)==0.) .and. (ddRho>0.) .and. &
            (deltaRhoAtKm1(i)<CS%MLE_density_diff) .and. (deltaRhoAtK(i)>=CS%MLE_density_diff)) then
          aFac = ( CS%MLE_density_diff - deltaRhoAtKm1(i) ) / ddRho
          MLD_fast(i,j) = dK(i) * aFac + dKm1(i) * (1. - aFac)
        endif
      enddo ! i-loop
    enddo ! k-loop
    do i=is-1,ie+1
      MLD_fast(i,j) = CS%MLE_MLD_stretch * MLD_fast(i,j)
      if ((MLD_fast(i,j)==0.) .and. (deltaRhoAtK(i)<CS%MLE_density_diff)) &
        MLD_fast(i,j) = dK(i) ! Assume mixing to the bottom
    enddo
  enddo ! j-loop
end subroutine detect_mld

! NOTE: This function appears to change answers on some platforms, so it is
! currently unused in the model, but we intend to introduce it in the future.

!> Return the growth timescale for the submesoscale mixed layer eddies in [T ~> s]
real function growth_time(u_star, hBL, absf, h_neg, vonKar, Kv_rest, restrat_coef)
  real, intent(in) :: u_star   !< Surface friction velocity in thickness-based units [H T-1 ~> m s-1 or kg m-2 s-1]
  real, intent(in) :: hBL      !< Boundary layer thickness including at least a negligible
                               !! value to keep it positive definite [H ~> m or kg m-2]
  real, intent(in) :: absf     !< Absolute value of the Coriolis parameter [T-1 ~> s-1]
  real, intent(in) :: h_neg    !< A tiny thickness that is usually lost in roundoff so can be
                               !! neglected [H ~> m or kg m-2]
  real, intent(in) :: Kv_rest  !< The background laminar vertical viscosity used for restratification,
                               !! rescaled into thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1]
  real, intent(in) :: vonKar   !< The von Karman constant, used to scale the turbulent limits
                               !! on the restratification timescales [nondim]
  real, intent(in) :: restrat_coef !< An overall scaling factor for the restratification timescale [nondim]

  ! Local variables
  real :: mom_mixrate  ! rate at which momentum is homogenized within mixed layer [T-1 ~> s-1]
  real :: Kv_eff       ! An effective overall viscosity in thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1]
  real :: pi2          ! A scaling constant that is approximately pi^2 [nondim]

  ! peak ML visc: u_star * von_Karman * (h_ml*u_star)/(absf*h_ml + 4.0*u_star) + Kv_water
  ! momentum mixing rate: pi^2*visc/h_ml^2
  pi2 = 9.8696  ! Approximately pi^2.  This is more accurate than the overall uncertainty of the
                ! scheme, with a value that is chosen to reproduce previous answers.
  if (Kv_rest <= 0.0) then
    ! This case reproduces the previous answers, but the extra h_neg is otherwise unnecessary.
    mom_mixrate = (pi2*vonKar)*u_star**2 / (absf*hBL**2 + 4.0*(hBL + h_neg)*u_star)
    growth_time = restrat_coef * (0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2))
  else
    ! Set the mixing rate to the sum of a turbulent mixing rate and a laminar viscous rate.
    ! mom_mixrate = pi2*vonKar*u_star**2 / (absf*hBL**2 + 4.0*hBL*u_star) + pi2*Kv_rest / hBL**2
    if (absf*hBL <= 4.0e-16*u_star) then
      Kv_eff = pi2 * (Kv_rest + 0.25*vonKar*hBL*u_star)
    else
      Kv_eff = pi2 * (Kv_rest + vonKar*u_star**2*hBL / (absf*hBL + 4.0*u_star))
    endif
    growth_time = (restrat_coef*0.0625) * ((hBL**2*(hBL**2*absf + 2.0*Kv_eff)) / ((hBL**2*absf)**2 + Kv_eff**2))
  endif

end function growth_time

!> Initialize the mixed layer restratification module
logical function mixedlayer_restrat_init(Time, G, GV, US, param_file, diag, CS, restart_CS)
  type(time_type), target,     intent(in)    :: Time       !< Current model time
  type(ocean_grid_type),       intent(inout) :: G          !< Ocean grid structure
  type(verticalGrid_type),     intent(in)    :: GV         !< Ocean vertical grid structure
  type(unit_scale_type),       intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),       intent(in)    :: param_file !< Parameter file to parse
  type(diag_ctrl), target,     intent(inout) :: diag       !< Regulate diagnostics
  type(mixedlayer_restrat_CS), intent(inout) :: CS         !< Module control structure
  type(MOM_restart_CS),        intent(in)    :: restart_CS !< MOM restart control structure

  ! Local variables
  real :: flux_to_kg_per_s ! A unit conversion factor for fluxes. [kg T s-1 H-1 L-2 ~> kg m-3 or 1]
  real :: omega            ! The Earth's rotation rate [T-1 ~> s-1].
  real :: ustar_min_dflt   ! The default value for RESTRAT_USTAR_MIN [Z T-1 ~> m s-1]
  real :: Stanley_coeff    ! Coefficient relating the temperature gradient and sub-gridscale
                           ! temperature variance [nondim]
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags
  ! This include declares and sets the variable "version".
  character(len=200) :: inputdir   ! The directory where NetCDF input files
  character(len=240) :: mle_fl_filename ! A file from which chl_a concentrations are to be read.
  character(len=128) :: mle_fl_file ! Data containing MLE front-length scale. Used
                                    ! when reading from file.
  character(len=32)  :: fl_varname ! Name of front-length scale variable in mle_fl_file.

# include "version_variable.h"
  integer :: i, j
  character(len=200) :: filename, varname

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "MIXEDLAYER_RESTRAT", mixedlayer_restrat_init, &
             default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, "", all_default=.not.mixedlayer_restrat_init)
  call get_param(param_file, mdl, "MIXEDLAYER_RESTRAT", mixedlayer_restrat_init, &
             "If true, a density-gradient dependent re-stratifying "//&
             "flow is imposed in the mixed layer. Can be used in ALE mode "//&
             "without restriction but in layer mode can only be used if "//&
             "BULKMIXEDLAYER is true.", default=.false.)
  if (.not. mixedlayer_restrat_init) return

  CS%initialized = .true.
  CS%Time => Time

  ! Nonsense values to cause problems when these parameters are not used
  CS%MLE_MLD_decay_time = -9.e9*US%s_to_T
  CS%MLE_density_diff = -9.e9*US%kg_m3_to_R
  CS%MLE_tail_dh = -9.e9
  CS%MLE_use_PBL_MLD = .false.
  CS%MLE_MLD_stretch = -9.e9
  CS%use_Stanley_ML = .false.
  CS%use_Bodner = .false.
  CS%fl_from_file = .false.
  CS%MLD_grid = .false.
  CS%Cr_grid = .false.

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
      "This sets the default value for the various _ANSWER_DATE parameters.", &
      default=99991231, do_not_log=.true.)
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  call openParameterBlock(param_file,'MLE') ! Prepend MLE% to all parameters
  if (GV%nkml==0) then
    call get_param(param_file, mdl, "USE_BODNER23", CS%use_Bodner, &
             "If true, use the Bodner et al., 2023, formulation of the re-stratifying "//&
             "mixed-layer restratification parameterization. This only works in ALE mode.", &
             default=.false.)
  endif
  if (CS%use_Bodner) then
    call get_param(param_file, mdl, "CR", CS%Cr, &
             "The efficiency coefficient in eq 27 of Bodner et al., 2023.", &
             units="nondim", default=0.0)
    call get_param(param_file, mdl, "BODNER_NSTAR", CS%Nstar, &
             "The n* value used to estimate the turbulent vertical momentum flux "//&
             "in Bodner et al., 2023, eq. 18. This is independent of the value used in "//&
             "the PBL scheme but should be set to be the same for consistency.", &
             units="nondim", default=0.066)
    call get_param(param_file, mdl, "BODNER_MSTAR", CS%Mstar, &
             "The m* value used to estimate the turbulent vertical momentum flux "//&
             "in Bodner et al., 2023, eq. 18. This is independent of the value used in "//&
             "the PBL scheme but should be set to be the same for consistency.", &
             units="nondim", default=0.5)
    call get_param(param_file, mdl, "BLD_GROWING_TFILTER", CS%BLD_growing_Tfilt, &
             "The time-scale for a running-mean filter applied to the boundary layer "//&
             "depth (BLD) when the BLD is deeper than the running mean. A value of 0 "//&
             "instantaneously sets the running mean to the current value of BLD.", &
             units="s", default=0., scale=US%s_to_T)
    call get_param(param_file, mdl, "BLD_DECAYING_TFILTER", CS%BLD_decaying_Tfilt, &
             "The time-scale for a running-mean filter applied to the boundary layer "//&
             "depth (BLD) when the BLD is shallower than the running mean. A value of 0 "//&
             "instantaneously sets the running mean to the current value of BLD.", &
             units="s", default=0., scale=US%s_to_T)
    call get_param(param_file, mdl, "MLD_GROWING_TFILTER", CS%MLD_growing_Tfilt, &
             "The time-scale for a running-mean filter applied to the time-filtered "//&
             "BLD, when the latter is deeper than the running mean. A value of 0 "//&
             "instantaneously sets the running mean to the current value filtered BLD.", &
             units="s", default=0., scale=US%s_to_T)
    call get_param(param_file, mdl, "MLD_DECAYING_TFILTER", CS%MLD_decaying_Tfilt, &
             "The time-scale for a running-mean filter applied to the time-filtered "//&
             "BLD, when the latter is shallower than the running mean. A value of 0 "//&
             "instantaneously sets the running mean to the current value filtered BLD.", &
             units="s", default=0., scale=US%s_to_T)
    call get_param(param_file, mdl, "ML_RESTRAT_ANSWER_DATE", CS%answer_date, &
             "The vintage of the order of arithmetic and expressions in the mixed layer "//&
             "restrat calculations.  Values below 20240201 recover the answers from the end "//&
             "of 2023, while higher values use the new cuberoot function in the Bodner code "//&
             "to avoid needing to undo dimensional rescaling.", &
             default=default_answer_date, &
             do_not_log=.not.(CS%use_Bodner.and.(GV%Boussinesq.or.GV%semi_Boussinesq)))
    call get_param(param_file, mdl, "MIN_WSTAR2", CS%min_wstar2, &
             "The minimum lower bound to apply to the vertical momentum flux, w'u', "//&
             "in the Bodner et al., restratification parameterization. This avoids "//&
             "a division-by-zero in the limit when u* and the buoyancy flux are zero.  "//&
             "The default is less than the molecular viscosity of water times the Coriolis "//&
             "parameter a micron away from the equator.", &
             units="m2 s-2", default=1.0e-24, scale=US%m_to_Z**2*US%T_to_s**2)
    call get_param(param_file, mdl, "TAIL_DH", CS%MLE_tail_dh, &
             "Fraction by which to extend the mixed-layer restratification "//&
             "depth used for a smoother stream function at the base of "//&
             "the mixed-layer.", units="nondim", default=0.0)
    call get_param(param_file, mdl, "USE_STANLEY_TVAR", CS%use_Stanley_ML, &
             "If true, turn on Stanley SGS T variance parameterization "// &
             "in ML restrat code.", default=.false.)
    call get_param(param_file, mdl, "USE_CR_GRID", CS%Cr_grid, &
             "If true, read in a spatially varying Cr field.", default=.false.)
    call get_param(param_file, mdl, "USE_MLD_GRID", CS%MLD_grid, &
             "If true, read in a spatially varying MLD_decaying_Tfilt field.", default=.false.)
    if (CS%MLD_grid) then
      call get_param(param_file, mdl, "MLD_TFILT_FILE", filename, &
             "The path to the file containing the MLD_decaying_Tfilt fields.", &
             default="")
      call get_param(param_file, mdl, "MLD_TFILT_VAR", varname, &
              "The variable name for MLD_decaying_Tfilt field.", &
              default="MLD_tfilt")
      filename = trim(inputdir) // "/" // trim(filename)
      allocate(CS%MLD_Tfilt_space(G%isd:G%ied,G%jsd:G%jed), source=0.0)
      call MOM_read_data(filename, varname, CS%MLD_Tfilt_space, G%domain, scale=US%s_to_T)
      call pass_var(CS%MLD_Tfilt_space, G%domain)
    endif
    allocate(CS%Cr_space(G%isd:G%ied,G%jsd:G%jed), source=CS%Cr)
    if (CS%Cr_grid) then
      call get_param(param_file, mdl, "CR_FILE", filename, &
             "The path to the file containing the Cr fields.", &
             default="")
      call get_param(param_file, mdl, "CR_VAR", varname, &
              "The variable name for Cr field.", &
              default="Cr")
      filename = trim(inputdir) // "/" // trim(filename)
      call MOM_read_data(filename, varname, CS%Cr_space, G%domain)
      call pass_var(CS%Cr_space, G%domain)
    endif
    call closeParameterBlock(param_file) ! The remaining parameters do not have MLE% prepended
    call get_param(param_file, mdl, "MLE_USE_PBL_MLD", CS%MLE_use_PBL_MLD, &
             "If true, the MLE parameterization will use the mixed-layer "//&
             "depth provided by the active PBL parameterization. If false, "//&
             "MLE will estimate a MLD based on a density difference with the "//&
             "surface using the parameter MLE_DENSITY_DIFF, unless "//&
              "BODNER_DETECT_MLD is true.", default=.false.)
    call get_param(param_file, mdl, "BODNER_DETECT_MLD", CS%Bodner_detect_MLD, &
             "If true, the Bodner parameterization will use the mixed-layer depth "//&
             "detected via the density difference criterion MLE_DENSITY_DIFF.", default=.false.)
    if (.not.(CS%MLE_use_PBL_MLD.or.CS%Bodner_detect_MLD)) call MOM_error(FATAL, "mixedlayer_restrat_init: "// &
             "To use MLE%USE_BODNER23=True then MLE_USE_PBL_MLD or BODNER_DETECT_MLD must be true.")
    if (CS%MLE_use_PBL_MLD.and.CS%Bodner_detect_MLD) call MOM_error(FATAL, "mixedlayer_restrat_init: "// &
             "MLE_USE_PBL_MLD and BODNER_DETECT_MLD cannot both be true.")
  else
    call closeParameterBlock(param_file) ! The remaining parameters do not have MLE% prepended
  endif

  if (.not.CS%use_Bodner) then
    ! This coefficient is used in both layered and ALE versions of Fox-Kemper but not Bodner
    call get_param(param_file, mdl, "FOX_KEMPER_ML_RESTRAT_COEF", CS%ml_restrat_coef, &
             "A nondimensional coefficient that is proportional to "//&
             "the ratio of the deformation radius to the dominant "//&
             "lengthscale of the submesoscale mixed layer "//&
             "instabilities, times the minimum of the ratio of the "//&
             "mesoscale eddy kinetic energy to the large-scale "//&
             "geostrophic kinetic energy or 1 plus the square of the "//&
             "grid spacing over the deformation radius, as detailed "//&
             "by Fox-Kemper et al. (2011)", units="nondim", default=0.0)
    ! These parameters are only used in the OM4-era version of Fox-Kemper
    call get_param(param_file, mdl, "USE_STANLEY_ML", CS%use_Stanley_ML, &
                   "If true, turn on Stanley SGS T variance parameterization "// &
                   "in ML restrat code.", default=.false.)
    if (CS%use_Stanley_ML) then
      call get_param(param_file, mdl, "STANLEY_COEFF", Stanley_coeff, &
                   "Coefficient correlating the temperature gradient and SGS T variance.", &
                   units="nondim", default=-1.0, do_not_log=.true.)
      if (Stanley_coeff < 0.0) call MOM_error(FATAL, &
               "STANLEY_COEFF must be set >= 0 if USE_STANLEY_ML is true.")
    endif
    call get_param(param_file, mdl, 'VON_KARMAN_CONST', CS%vonKar, &
                   'The value the von Karman constant as used for mixed layer viscosity.', &
                   units='nondim', default=0.41)
    ! We use GV%nkml to distinguish between the old and new implementation of MLE.
    ! The old implementation only works for the layer model with nkml>0.
    if (GV%nkml==0) then
      call get_param(param_file, mdl, "FOX_KEMPER_ML_RESTRAT_COEF2", CS%ml_restrat_coef2, &
             "As for FOX_KEMPER_ML_RESTRAT_COEF but used in a second application "//&
             "of the MLE restratification parameterization.", units="nondim", default=0.0)
      call get_param(param_file, mdl, "MLE_FRONT_LENGTH", CS%front_length, &
             "If non-zero, is the frontal-length scale used to calculate the "//&
             "upscaling of buoyancy gradients that is otherwise represented "//&
             "by the parameter FOX_KEMPER_ML_RESTRAT_COEF. If MLE_FRONT_LENGTH is "//&
             "non-zero, it is recommended to set FOX_KEMPER_ML_RESTRAT_COEF=1.0.",&
             units="m", default=0.0, scale=US%m_to_L)
      call get_param(param_file, mdl, "MLE_FRONT_LENGTH_FROM_FILE", CS%fl_from_file, &
                   "If true, the MLE front-length scale is read from a file.", default=.false.)
      if (CS%fl_from_file) then
        call time_interp_external_init()

        call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
        call get_param(param_file, mdl, "MLE_FL_FILE", mle_fl_file, &
                   "MLE_FL_FILE is the file containing MLE front-length scale. "//&
                   "It is used when MLE_FRONT_LENGTH_FROM_FILE is true.", fail_if_missing=.true.)
        mle_fl_filename = trim(slasher(inputdir))//trim(mle_fl_file)
        call log_param(param_file, mdl, "INPUTDIR/MLE_FL_FILE", mle_fl_filename)
        call get_param(param_file, mdl, "FL_VARNAME", fl_varname, &
                   "Name of MLE front-length scale variable in MLE_FL_FILE.", default='mle_fl')
        if (modulo(G%Domain%turns, 4) /= 0) then
          CS%sbc_fl = init_external_field(mle_fl_filename, trim(fl_varname), MOM_domain=G%Domain%domain_in)
        else
          CS%sbc_fl = init_external_field(mle_fl_filename, trim(fl_varname), MOM_domain=G%Domain)
        endif
      endif
      if (CS%fl_from_file .and. CS%front_length>0.0) call MOM_error(FATAL, "mixedlayer_restrat_init: "// &
             "MLE_FRONT_LENGTH_FROM_FILE cannot be true when MLE_FRONT_LENGTH > 0.0. "// &
             "If you want to use MLE_FRONT_LENGTH, set MLE_FRONT_LENGTH_FROM_FILE to false." // &
             "If you want to use MLE_FRONT_LENGTH_FROM_FILE, set MLE_FRONT_LENGTH to 0.0.")
      call get_param(param_file, mdl, "MLE_USE_PBL_MLD", CS%MLE_use_PBL_MLD, &
             "If true, the MLE parameterization will use the mixed-layer "//&
             "depth provided by the active PBL parameterization. If false, "//&
             "MLE will estimate a MLD based on a density difference with the "//&
             "surface using the parameter MLE_DENSITY_DIFF.", default=.false.)
      call get_param(param_file, mdl, "MLE_MLD_DECAY_TIME", CS%MLE_MLD_decay_time, &
             "The time-scale for a running-mean filter applied to the mixed-layer "//&
             "depth used in the MLE restratification parameterization. When "//&
             "the MLD deepens below the current running-mean the running-mean "//&
             "is instantaneously set to the current MLD.", units="s", default=0., scale=US%s_to_T)
      call get_param(param_file, mdl, "MLE_MLD_DECAY_TIME2", CS%MLE_MLD_decay_time2, &
             "The time-scale for a running-mean filter applied to the filtered "//&
             "mixed-layer depth used in a second MLE restratification parameterization. "//&
             "When the MLD deepens below the current running-mean the running-mean "//&
             "is instantaneously set to the current MLD.", units="s", default=0., scale=US%s_to_T)
      if (.not. CS%MLE_use_PBL_MLD) then
        call get_param(param_file, mdl, "MLE_DENSITY_DIFF", CS%MLE_density_diff, &
             "Density difference used to detect the mixed-layer "//&
             "depth used for the mixed-layer eddy parameterization "//&
             "by Fox-Kemper et al. (2011)", units="kg/m3", default=0.03, scale=US%kg_m3_to_R)
      endif
      call get_param(param_file, mdl, "MLE_TAIL_DH", CS%MLE_tail_dh, &
             "Fraction by which to extend the mixed-layer restratification "//&
             "depth used for a smoother stream function at the base of "//&
             "the mixed-layer.", units="nondim", default=0.0)
      call get_param(param_file, mdl, "MLE_MLD_STRETCH", CS%MLE_MLD_stretch, &
             "A scaling coefficient for stretching/shrinking the MLD "//&
             "used in the MLE scheme. This simply multiplies MLD wherever used.",&
             units="nondim", default=1.0)
    endif
    call get_param(param_file, mdl, "KV_RESTRAT", CS%Kv_restrat, &
                 "A small viscosity that sets a floor on the momentum mixing rate during "//&
                 "restratification.  If this is positive, it will prevent some possible "//&
                 "divisions by zero even if ustar, RESTRAT_USTAR_MIN, and f are all 0.", &
                 units="m2 s-1", default=0.0, scale=GV%m2_s_to_HZ_T*(US%Z_to_m*GV%m_to_H))
    call get_param(param_file, mdl, "OMEGA", omega, &
                 "The rotation rate of the earth.", &
                 units="s-1", default=7.2921e-5, scale=US%T_to_s)
    ustar_min_dflt = 2.0e-4 * omega * (GV%Angstrom_Z + GV%dZ_subroundoff)
    call get_param(param_file, mdl, "RESTRAT_USTAR_MIN", CS%ustar_min, &
                 "The minimum value of ustar that will be used by the mixed layer "//&
                 "restratification module.  This can be tiny, but if this is greater than 0, "//&
                 "it will prevent divisions by zero when f and KV_RESTRAT are zero.", &
                 units="m s-1", default=US%Z_to_m*US%s_to_T*ustar_min_dflt, scale=GV%m_to_H*US%T_to_s)
  elseif (CS%Bodner_detect_MLD) then
    call get_param(param_file, mdl, "MLE_DENSITY_DIFF", CS%MLE_density_diff, &
           "Density difference used to detect the mixed-layer "//&
           "depth used for the mixed-layer eddy parameterization "//&
           "by Fox-Kemper et al. (2010)", units="kg/m3", default=0.03, scale=US%kg_m3_to_R)
    call get_param(param_file, mdl, "MLE_MLD_STRETCH", CS%MLE_MLD_stretch, &
           "A scaling coefficient for stretching/shrinking the MLD "//&
           "used in the MLE scheme. This simply multiplies MLD wherever used.",&
           units="nondim", default=1.0)
  endif

  CS%diag => diag

  flux_to_kg_per_s = GV%H_to_kg_m2 * US%L_to_m**2 * US%s_to_T

  CS%id_uhml = register_diag_field('ocean_model', 'uhml', diag%axesCuL, Time, &
      'Zonal Thickness Flux to Restratify Mixed Layer', &
      'kg s-1', conversion=flux_to_kg_per_s, y_cell_method='sum', v_extensive=.true.)
  CS%id_vhml = register_diag_field('ocean_model', 'vhml', diag%axesCvL, Time, &
      'Meridional Thickness Flux to Restratify Mixed Layer', &
      'kg s-1', conversion=flux_to_kg_per_s, x_cell_method='sum', v_extensive=.true.)
  CS%id_urestrat_time = register_diag_field('ocean_model', 'MLu_restrat_time', diag%axesCu1, Time, &
      'Mixed Layer Zonal Restratification Timescale', 's', conversion=US%T_to_s)
  CS%id_vrestrat_time = register_diag_field('ocean_model', 'MLv_restrat_time', diag%axesCv1, Time, &
      'Mixed Layer Meridional Restratification Timescale', 's', conversion=US%T_to_s)
  CS%id_MLD = register_diag_field('ocean_model', 'MLD_restrat', diag%axesT1, Time, &
      'Mixed Layer Depth as used in the mixed-layer restratification parameterization', &
      'm', conversion=GV%H_to_m)
  CS%id_BLD = register_diag_field('ocean_model', 'BLD_restrat', diag%axesT1, Time, &
      'Boundary Layer Depth as used in the mixed-layer restratification parameterization', &
      'm', conversion=GV%H_to_m)
  CS%id_Rml = register_diag_field('ocean_model', 'ML_buoy_restrat', diag%axesT1, Time, &
      'Mixed Layer Buoyancy as used in the mixed-layer restratification parameterization', &
      'm s-2', conversion=GV%m_to_H*(US%L_T_to_m_s**2))
  CS%id_uDml = register_diag_field('ocean_model', 'udml_restrat', diag%axesCu1, Time, &
      'Transport stream function amplitude for zonal restratification of mixed layer', &
      'm3 s-1', conversion=GV%H_to_m*(US%L_to_m**2)*US%s_to_T)
  CS%id_vDml = register_diag_field('ocean_model', 'vdml_restrat', diag%axesCv1, Time, &
      'Transport stream function amplitude for meridional restratification of mixed layer', &
      'm3 s-1', conversion=GV%H_to_m*(US%L_to_m**2)*US%s_to_T)
  CS%id_uml = register_diag_field('ocean_model', 'uml_restrat', diag%axesCu1, Time, &
      'Surface zonal velocity component of mixed layer restratification', &
      'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vml = register_diag_field('ocean_model', 'vml_restrat', diag%axesCv1, Time, &
      'Surface meridional velocity component of mixed layer restratification', &
      'm s-1', conversion=US%L_T_to_m_s)
  if (CS%use_Bodner) then
    CS%id_wpup = register_diag_field('ocean_model', 'MLE_wpup', diag%axesT1, Time, &
        'Vertical turbulent momentum flux in Bodner mixed layer restratification parameterization', &
        'm2 s-2', conversion=US%L_to_m*GV%H_to_m*US%s_to_T**2)
    CS%id_ustar = register_diag_field('ocean_model', 'MLE_ustar', diag%axesT1, Time, &
        'Surface turbulent friction velocity, u*, in Bodner mixed layer restratification parameterization', &
        'm s-1', conversion=(US%Z_to_m*US%s_to_T))
    CS%id_bflux = register_diag_field('ocean_model', 'MLE_bflux', diag%axesT1, Time, &
        'Surface buoyancy flux, B0, in Bodner mixed layer restratification parameterization', &
        'm2 s-3', conversion=(US%Z_to_m**2*US%s_to_T**3))
    CS%id_lfbod = register_diag_field('ocean_model', 'lf_bodner', diag%axesT1, Time, &
        'Front length in Bodner mixed layer restratificiation parameterization', &
        'm', conversion=US%L_to_m)
  else
    CS%id_mle_fl = register_diag_field('ocean_model', 'mle_fl', diag%axesT1, Time, &
        'Frontal length scale used in the mixed layer restratificiation parameterization', &
        'm', conversion=US%L_to_m)
  endif

  ! If MLD_filtered is being used, we need to update halo regions after a restart
  if (allocated(CS%MLD_filtered)) call pass_var(CS%MLD_filtered, G%domain)
  if (allocated(CS%MLD_filtered_slow)) call pass_var(CS%MLD_filtered_slow, G%domain)
  if (allocated(CS%wpup_filtered)) call pass_var(CS%wpup_filtered, G%domain)

end function mixedlayer_restrat_init

!> Allocate and register fields in the mixed layer restratification structure for restarts
subroutine mixedlayer_restrat_register_restarts(HI, GV, US, param_file, CS, restart_CS)
  ! Arguments
  type(hor_index_type),        intent(in)    :: HI         !< Horizontal index structure
  type(verticalGrid_type),     intent(in)    :: GV         !< Ocean vertical grid structure
  type(unit_scale_type),       intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),       intent(in)    :: param_file !< Parameter file to parse
  type(mixedlayer_restrat_CS), intent(inout) :: CS         !< Module control structure
  type(MOM_restart_CS),        intent(inout) :: restart_CS !< MOM restart control structure

  ! Local variables
  character(len=64) :: mom_flux_units
  logical :: mixedlayer_restrat_init, use_Bodner

  ! Check to see if this module will be used
  call get_param(param_file, mdl, "MIXEDLAYER_RESTRAT", mixedlayer_restrat_init, &
             default=.false., do_not_log=.true.)
  if (.not. mixedlayer_restrat_init) return

  call get_param(param_file, mdl, "MLE_MLD_DECAY_TIME", CS%MLE_MLD_decay_time, &
                 units="s", default=0., scale=US%s_to_T, do_not_log=.true.)
  call get_param(param_file, mdl, "MLE_MLD_DECAY_TIME2", CS%MLE_MLD_decay_time2, &
                 units="s", default=0., scale=US%s_to_T, do_not_log=.true.)
  call openParameterBlock(param_file, 'MLE', do_not_log=.true.)
  call get_param(param_file, mdl, "USE_BODNER23", use_Bodner, &
                 default=.false., do_not_log=.true.)
  call closeParameterBlock(param_file)
  if (CS%MLE_MLD_decay_time>0. .or. CS%MLE_MLD_decay_time2>0. .or. use_Bodner) then
    ! CS%MLD_filtered is used to keep a running mean of the PBL's actively mixed MLD.
    allocate(CS%MLD_filtered(HI%isd:HI%ied,HI%jsd:HI%jed), source=0.)
    call register_restart_field(CS%MLD_filtered, "MLD_MLE_filtered", .false., restart_CS, &
                                longname="Time-filtered MLD for use in MLE", &
                                units=get_thickness_units(GV), conversion=GV%H_to_MKS)
  endif
  if (CS%MLE_MLD_decay_time2>0. .or. use_Bodner) then
    ! CS%MLD_filtered_slow is used to keep a running mean of the PBL's seasonal or winter MLD.
    allocate(CS%MLD_filtered_slow(HI%isd:HI%ied,HI%jsd:HI%jed), source=0.)
    call register_restart_field(CS%MLD_filtered_slow, "MLD_MLE_filtered_slow", .false., restart_CS, &
                                longname="Slower time-filtered MLD for use in MLE", &
                                units=get_thickness_units(GV), conversion=GV%H_to_MKS)
  endif
  if (use_Bodner) then
    ! CS%MLD_filtered_slow is used to keep a running mean of the PBL's seasonal or winter MLD.
    mom_flux_units = "m2 s-2" ; if (.not.GV%Boussinesq) mom_flux_units = "kg m-1 s-2"
    allocate(CS%wpup_filtered(HI%isd:HI%ied,HI%jsd:HI%jed), source=0.)
    call register_restart_field(CS%wpup_filtered, "MLE_Bflux", .false., restart_CS, &
                                longname="Time-filtered vertical turbulent momentum flux for use in MLE", &
                                units=mom_flux_units, conversion=US%L_to_m*GV%H_to_mks*US%s_to_T**2 )
  endif

end subroutine mixedlayer_restrat_register_restarts

!> Returns true if a unit test of functions in MOM_mixedlayer_restrat fail.
!! Returns false otherwise.
logical function mixedlayer_restrat_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables
  type(mixedlayer_restrat_CS) :: CS ! Control structure
  logical :: this_test

  print *,'===== mixedlayer_restrat: mixedlayer_restrat_unit_tests =================='

  ! Tests of the shape function mu(z)
  this_test = &
    test_answer(verbose, mu(3.,0.), 0., 'mu(3)=0')
  this_test = this_test .or. &
    test_answer(verbose, mu(0.,0.), 0., 'mu(0)=0')
  this_test = this_test .or. &
    test_answer(verbose, mu(-0.25,0.), 0.7946428571428572, 'mu(-0.25)=0.7946...', tol=epsilon(1.))
  this_test = this_test .or. &
    test_answer(verbose, mu(-0.5,0.), 1., 'mu(-0.5)=1')
  this_test = this_test .or. &
    test_answer(verbose, mu(-0.75,0.), 0.7946428571428572, 'mu(-0.75)=0.7946...', tol=epsilon(1.))
  this_test = this_test .or. &
    test_answer(verbose, mu(-1.,0.), 0., 'mu(-1)=0')
  this_test = this_test .or. &
    test_answer(verbose, mu(-3.,0.), 0., 'mu(-3)=0')
  this_test = this_test .or. &
    test_answer(verbose, mu(-0.5,0.5), 1., 'mu(-0.5,0.5)=1')
  this_test = this_test .or. &
    test_answer(verbose, mu(-1.,0.5), 0.25, 'mu(-1,0.5)=0.25')
  this_test = this_test .or. &
    test_answer(verbose, mu(-1.5,0.5), 0., 'mu(-1.5,0.5)=0')
  if (.not. this_test) print '(a)','  Passed tests of mu(z)'
  mixedlayer_restrat_unit_tests = this_test

  ! Tests of the two time-scale running mean function
  this_test = &
    test_answer(verbose, rmean2ts(3.,2.,0.,0.,3.), 3., 'rmean2ts(3,2,0,0,3)=3')
  this_test = this_test .or. &
    test_answer(verbose, rmean2ts(1.,2.,0.,0.,3.), 1., 'rmean2ts(1,2,0,0,3)=1')
  this_test = this_test .or. &
    test_answer(verbose, rmean2ts(4.,0.,3.,0.,1.), 1., 'rmean2ts(4,0,3,0,1)=1')
  this_test = this_test .or. &
    test_answer(verbose, rmean2ts(0.,4.,0.,3.,1.), 3., 'rmean2ts(0,4,0,3,1)=3')
  if (.not. this_test) print '(a)','  Passed tests of rmean2ts(s,f,g,d,dt)'
  mixedlayer_restrat_unit_tests = mixedlayer_restrat_unit_tests .or. this_test

end function mixedlayer_restrat_unit_tests

!> Returns true if any cell of u and u_true are not identical. Returns false otherwise.
logical function test_answer(verbose, u, u_true, label, tol)
  logical,            intent(in) :: verbose !< If true, write results to stdout
  real,               intent(in) :: u      !< Values to test in arbitrary units [A]
  real,               intent(in) :: u_true !< Values to test against (correct answer) [A]
  character(len=*),   intent(in) :: label  !< Message
  real, optional,     intent(in) :: tol    !< The tolerance for differences between u and u_true [A]
  ! Local variables
  real :: tolerance ! The tolerance for differences between u and u_true [A]
  integer :: k

  tolerance = 0.0 ; if (present(tol)) tolerance = tol
  test_answer = .false.

  if (abs(u - u_true) > tolerance) test_answer = .true.
  if (test_answer .or. verbose) then
    if (test_answer) then
      print '(3(a,1pe24.16),1x,a,1x,a)','computed =',u,' correct =',u_true, &
            ' err=',u-u_true,' < wrong',label
    else
      print '(2(a,1pe24.16),1x,a)','computed =',u,' correct =',u_true,label
    endif
  endif

end function test_answer

!> \namespace mom_mixed_layer_restrat
!!
!! \section section_mle Mixed-layer eddy parameterization module
!!
!! The subroutines in this module implement a parameterization of unresolved viscous
!! mixed layer restratification of the mixed layer as described in \cite fox-kemper2008,
!! and whose impacts are described in \cite fox-kemper2011.
!! This is derived in part from the older parameterization that is described in
!! \cite Hallberg2003, which this new parameterization surpasses, which
!! in turn is based on the sub-inertial mixed layer theory of \cite Young1994.
!! There is no net horizontal volume transport due to this parameterization, and
!! no direct effect below the mixed layer. A revised version of the parameterization by
!! \cite Bodner2023 is also available as an option.
!!
!! This parameterization sets the restratification timescale to agree with
!! high-resolution studies of mixed layer restratification.
!!
!! The run-time parameter <code>FOX_KEMPER_ML_RESTRAT_COEF</code> is a non-dimensional number of
!! order a few tens, proportional to the ratio of the deformation radius or the
!! grid scale (whichever is smaller) to the dominant horizontal length-scale of the
!! sub-meso-scale mixed layer instabilities.
!!
!! \subsection section_mle_nutshell "Sub-meso" in a nutshell
!!
!! The parameterization is colloquially referred to as "sub-meso".
!!
!! The original \cite fox-kemper2008-2 paper proposed a quasi-Stokes
!! advection described by the stream function (eq. 5 of \cite fox-kemper2011):
!! \f[
!!    {\bf \Psi}_o = C_e \frac{ H^2 \nabla \bar{b} \times \hat{\bf z} }{ |f| } \mu(z)
!! \f]
!!
!! where the vertical profile function is
!! \f[
!!    \mu(z) = \max \left\{ 0, \left[ 1 - \left(\frac{2z}{H}+1\right)^2 \right]
!!                            \left[ 1 + \frac{5}{21} \left(\frac{2z}{H}+1\right)^2 \right] \right\}
!! \f]
!! and \f$ H \f$ is the mixed-layer depth, \f$ f \f$ is the local Coriolis parameter, \f$ C_e \sim 0.06-0.08 \f$ and
!! \f$ \nabla \bar{b} \f$ is a depth mean buoyancy gradient averaged over the mixed layer.
!!
!! For use in coarse-resolution models, an upscaling of the buoyancy gradients and adaption for the equator
!! leads to the following parameterization (eq. 6 of \cite fox-kemper2011):
!! \f[
!!    {\bf \Psi} = C_e \Gamma_\Delta \frac{\Delta s}{l_f} \frac{ H^2 \nabla \bar{b} \times \hat{\bf z} }
!!                 { \sqrt{ f^2 + \tau^{-2}} } \mu(z)
!! \f]
!! where \f$ \Delta s \f$ is the minimum of grid-scale and deformation radius,
!! \f$ l_f \f$ is the width of the mixed-layer fronts, and \f$ \Gamma_\Delta=1 \f$.
!! \f$ \tau \f$ is a time-scale for mixing momentum across the mixed layer.
!! \f$ l_f \f$ is thought to be of order hundreds of meters.
!!
!! The upscaling factor \f$ \frac{\Delta s}{l_f} \f$ can be a global constant, model parameter
!! <code>FOX_KEMPER_ML_RESTRAT</code>, so that in practice the parameterization is:
!! \f[
!!    {\bf \Psi} = C_e \Gamma_\Delta \frac{ H^2 \nabla \bar{b} \times \hat{\bf z} }{ \sqrt{ f^2 + \tau^{-2}} } \mu(z)
!! \f]
!! with non-unity \f$ \Gamma_\Delta \f$.
!!
!! \f$ C_e \f$ is hard-coded as 0.0625. \f$ \tau \f$ is calculated from the surface friction velocity \f$ u^* \f$.
!!
!! \todo Explain expression for momentum mixing time-scale.
!!
!! | Symbol                       | Module parameter      |
!! | ---------------------------- | --------------------- |
!! | \f$ \Gamma_\Delta \f$        | FOX_KEMPER_ML_RESTRAT |
!! | \f$ l_f \f$                  | MLE_FRONT_LENGTH      |
!! | \f$ \Delta \rho \f$          | MLE_DENSITY_DIFF      |
!!
!! \subsection section_mle_filtering Time-filtering of mixed-layer depth
!!
!! Using the instantaneous mixed-layer depth is inconsistent with the finite life-time of
!! mixed-layer instabilities. We provide a one-sided running-mean filter of mixed-layer depth, \f$ H \f$, of the form:
!! \f[
!!    \bar{H} \leftarrow \max \left( H, \frac{ \Delta t H + \tau_h \bar{H} }{ \Delta t + \tau_h } \right)
!! \f]
!! which allows the effective mixed-layer depth seen by the parameterization, \f$\bar{H}\f$, to instantaneously deepen
!! but to decay with time-scale \f$ \tau_h \f$.
!! \f$ \bar{H} \f$ is substituted for \f$ H \f$ in the above equations.
!!
!! | Symbol                       | Module parameter      |
!! | ---------------------------- | --------------------- |
!! | \f$ \tau_h \f$               | MLE_MLD_DECAY_TIME    |
!!
!! \subsection section_mle_mld Defining the mixed-layer-depth
!!
!! If the parameter MLE_USE_PBL_MLD=True then the mixed-layer depth is defined/diagnosed by the
!! boundary-layer parameterization (e.g. ePBL, KPP, etc.).
!!
!! If the parameter MLE_USE_PBL_MLD=False then the mixed-layer depth is diagnosed in this module
!! as the depth of a given density difference, \f$ \Delta \rho \f$, with the surface where the
!! density difference is the parameter MLE_DENSITY_DIFF.
!!
!! \subsection The Bodner (2023) modification
!!
!! To use this variant of the parameterization, set MLE\%USE_BODNER23=True which then changes the
!! available parameters.
!! MLE_USE_PBL_MLD must be True to use the B23 modification.
!!
!! \cite Bodner2023, (B23) use an expression for the frontal width which changes the scaling from \f$ H^2 \f$
!! to \f$ h H^2 \f$:
!! \f[
!!    {\bf \Psi} = C_r \frac{\Delta s |f| \bar{h} \bar{H}^2 \nabla \bar{b} \times \hat{\bf z} }
!!                 { \left( m_*u_*^3 + n_* w_*^3 \right)^{2/3} } \mu(z)
!! \f]
!! (see eq. 27 of B23).
!! Here, the \f$h\f$ is the activate boundary layer depth, and \f$H\f$ is the mixed layer depth.
!! The denominator is an approximation of the vertical turbulent momentum flux \f$\overline{w'u'}\f$ (see
!! eq. 18 of B23) calculated from the surface friction velocity \f$u_*\f$, and from the surface buoyancy flux,
!! \f$B\f$, using the relation \f$ w_*^3 \sim -B h \f$.
!! An advantage of this form of "sub-meso" is the denominator is well behaved at the equator but we apply a
!! lower bound of \f$w_{min}^2\f$ to avoid division by zero under zero forcing.
!! As for the original Fox-Kemper parameterization, \f$\nabla \bar{b}\f$ is the buoyancy gradient averaged
!! over the mixed-layer.
!!
!! The instantaneous boundary layer depth, \f$h\f$, is time filtered primarily to remove the diurnal cycle:
!! \f[
!!    \bar{h} \leftarrow \max \left(
!!        \min \left( h, \frac{ \Delta t h + \tau_{h+} \bar{h} }{ \Delta t + \tau_{h+} } \right),
!!                   \frac{ \Delta t h + \tau_{h-} \bar{h} }{ \Delta t + \tau_{h-} } \right)
!! \f]
!! Setting \f$ \tau_{h+}=0 \f$ means that when \f$ h>\bar{h} \f$ then \f$\bar{h}\leftarrow h\f$, i.e. the
!! effective (filtered) depth, \f$\bar{h}\f$, is instantly deepened. When \f$h<\bar{h}\f$ then the effective
!! depth shoals with time-scale \f$\tau_{h-}\f$.
!!
!! A second filter is applied to \f$\bar{h}\f$ to yield and effective "mixed layer depth", \f$\bar{H}\f$,
!! defined as the deepest the boundary layer over some time-scale \f$\tau_{H-}\f$:
!! \f[
!!    \bar{H} \leftarrow \max \left(
!!        \min \left( \bar{h}, \frac{ \Delta t \bar{h} + \tau_{H+} \bar{H} }{ \Delta t + \tau_{H+} } \right),
!!                   \frac{ \Delta t \bar{h} + \tau_{h-} \bar{H} }{ \Delta t + \tau_{H-} } \right)
!! \f]
!! Again, setting \f$ \tau_{H+}=0 \f$ allows the effective mixed layer to instantly deepend to \f$ \bar{h} \f$.
!!
!! | Symbol                       | Module parameter          |
!! | ---------------------------- | ------------------------- |
!! | \f$ C_r \f$                  | MLE\%CR                   |
!! | \f$ n_* \f$                  | MLE\%BODNER_NSTAR         |
!! | \f$ m_* \f$                  | MLE\%BODNER_MSTAR         |
!! | \f$ w_* \f$                  | MLE\%BODNER_MSTAR         |
!! | \f$ w_{min}^2 \f$            | MLE\%MIN_WSTAR2           |
!! | \f$ \tau_{h+} \f$            | MLE\%BLD_GROWING_TFILTER  |
!! | \f$ \tau_{h-} \f$            | MLE\%BLD_DECAYING_TFILTER |
!! | \f$ \tau_{H+} \f$            | MLE\%MLD_GROWING_TFILTER  |
!! | \f$ \tau_{H-} \f$            | MLE\%MLD_DECAYING_TFILTER |
!!
!! \subsection section_mle_ref References
!!
!! Fox-Kemper, B., Ferrari, R. and Hallberg, R., 2008:
!! Parameterization of Mixed Layer Eddies. Part I: Theory and Diagnosis
!! J. Phys. Oceangraphy, 38 (6), p1145-1165.
!! https://doi.org/10.1175/2007JPO3792.1
!!
!! Fox-Kemper, B. and Ferrari, R. 2008:
!! Parameterization of Mixed Layer Eddies. Part II: Prognosis and Impact
!! J. Phys. Oceangraphy, 38 (6), p1166-1179.
!! https://doi.org/10.1175/2007JPO3788.1
!!
!! B. Fox-Kemper, G. Danabasoglu, R. Ferrari, S.M. Griffies, R.W. Hallberg, M.M. Holland, M.E. Maltrud,
!! S. Peacock, and B.L. Samuels, 2011: Parameterization of mixed layer eddies. III: Implementation and impact
!! in global ocean climate simulations. Ocean Modell., 39(1), p61-78.
!! https://doi.org/10.1016/j.ocemod.2010.09.002
!!
!! A.S. Bodner, B. Fox-Kemper, L. Johnson, L. P. Van Roekel, J. C. McWilliams, P. P. Sullivan, P. S. Hall,
!! and J. Dong, 2023: Modifying the Mixed Layer Eddy Parameterization to Include Frontogenesis Arrest by
!! Boundary Layer Turbulence. J. Phys. Oceanogr., 53(1), p323-339.
!! https://doi.org/10.1175/JPO-D-21-0297.1

end module MOM_mixed_layer_restrat
