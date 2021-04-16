!> A column-wise toolbox for implementing neutral diffusion
module MOM_neutral_diffusion

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,             only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_domains,               only : pass_var
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_diag_mediator,         only : post_data, register_diag_field
use MOM_EOS,                   only : EOS_type, EOS_manual_init, EOS_domain
use MOM_EOS,                   only : calculate_density, calculate_density_derivs
use MOM_EOS,                   only : extract_member_EOS, EOS_LINEAR, EOS_TEOS10, EOS_WRIGHT
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_file_parser,           only : openParameterBlock, closeParameterBlock
use MOM_grid,                  only : ocean_grid_type
use MOM_remapping,             only : remapping_CS, initialize_remapping
use MOM_remapping,             only : extract_member_remapping_CS, build_reconstructions_1d
use MOM_remapping,             only : average_value_ppoly, remappingSchemesDoc, remappingDefaultScheme
use MOM_tracer_registry,       only : tracer_registry_type, tracer_type
use MOM_unit_scaling,          only : unit_scale_type
use MOM_verticalGrid,          only : verticalGrid_type
use polynomial_functions,      only : evaluation_polynomial, first_derivative_polynomial
use PPM_functions,             only : PPM_reconstruction, PPM_boundary_extrapolation
use regrid_edge_values,        only : edge_values_implicit_h4
use MOM_CVMix_KPP,             only : KPP_get_BLD, KPP_CS
use MOM_energetic_PBL,         only : energetic_PBL_get_MLD, energetic_PBL_CS
use MOM_diabatic_driver,       only : diabatic_CS, extract_diabatic_member
use MOM_lateral_boundary_diffusion, only : boundary_k_range, SURFACE, BOTTOM
use MOM_io,                    only : stdout, stderr

implicit none ; private

#include <MOM_memory.h>

public neutral_diffusion, neutral_diffusion_init, neutral_diffusion_end
public neutral_diffusion_calc_coeffs
public neutral_diffusion_unit_tests

!> The control structure for the MOM_neutral_diffusion module
type, public :: neutral_diffusion_CS ; private
  integer :: nkp1     !< Number of interfaces for a column = nk + 1
  integer :: nsurf    !< Number of neutral surfaces
  integer :: deg = 2  !< Degree of polynomial used for reconstructions
  logical :: continuous_reconstruction = .true. !< True if using continuous PPM reconstruction at interfaces
  logical :: debug = .false. !< If true, write verbose debugging messages
  logical :: hard_fail_heff !< Bring down the model if a problem with heff is detected
  integer :: max_iter !< Maximum number of iterations if refine_position is defined
  real :: drho_tol    !< Convergence criterion representing density difference from true neutrality [R ~> kg m-3]
  real :: x_tol       !< Convergence criterion for how small an update of the position can be
  real :: ref_pres    !< Reference pressure, negative if using locally referenced neutral
                      !! density [R L2 T-2 ~> Pa]
  logical :: interior_only !< If true, only applies neutral diffusion in the ocean interior.
                      !! That is, the algorithm will exclude the surface and bottom boundary layers.
  logical :: use_unmasked_transport_bug !< If true, use an older form for the accumulation of
                      !! neutral-diffusion transports that were unmasked, as used prior to Jan 2018.
  ! Positions of neutral surfaces in both the u, v directions
  real,    allocatable, dimension(:,:,:) :: uPoL  !< Non-dimensional position with left layer uKoL-1, u-point
  real,    allocatable, dimension(:,:,:) :: uPoR  !< Non-dimensional position with right layer uKoR-1, u-point
  integer, allocatable, dimension(:,:,:) :: uKoL  !< Index of left interface corresponding to neutral surface,
                                                  !! at a u-point
  integer, allocatable, dimension(:,:,:) :: uKoR  !< Index of right interface corresponding to neutral surface,
                                                  !! at a u-point
  real,    allocatable, dimension(:,:,:) :: uHeff !< Effective thickness at u-point [H ~> m or kg m-2]
  real,    allocatable, dimension(:,:,:) :: vPoL  !< Non-dimensional position with left layer uKoL-1, v-point
  real,    allocatable, dimension(:,:,:) :: vPoR  !< Non-dimensional position with right layer uKoR-1, v-point
  integer, allocatable, dimension(:,:,:) :: vKoL  !< Index of left interface corresponding to neutral surface,
                                                  !! at a v-point
  integer, allocatable, dimension(:,:,:) :: vKoR  !< Index of right interface corresponding to neutral surface,
                                                  !! at a v-point
  real,    allocatable, dimension(:,:,:) :: vHeff !< Effective thickness at v-point [H ~> m or kg m-2]
  ! Coefficients of polynomial reconstructions for temperature and salinity
  real,    allocatable, dimension(:,:,:,:) :: ppoly_coeffs_T !< Polynomial coefficients for temperature
  real,    allocatable, dimension(:,:,:,:) :: ppoly_coeffs_S !< Polynomial coefficients for salinity
  ! Variables needed for continuous reconstructions
  real,    allocatable, dimension(:,:,:) :: dRdT !< dRho/dT [R degC-1 ~> kg m-3 degC-1] at interfaces
  real,    allocatable, dimension(:,:,:) :: dRdS !< dRho/dS [R ppt-1 ~> kg m-3 ppt-1] at interfaces
  real,    allocatable, dimension(:,:,:) :: Tint !< Interface T [degC]
  real,    allocatable, dimension(:,:,:) :: Sint !< Interface S [ppt]
  real,    allocatable, dimension(:,:,:) :: Pint !< Interface pressure [R L2 T-2 ~> Pa]
  ! Variables needed for discontinuous reconstructions
  real,    allocatable, dimension(:,:,:,:) :: T_i    !< Top edge reconstruction of temperature [degC]
  real,    allocatable, dimension(:,:,:,:) :: S_i    !< Top edge reconstruction of salinity [ppt]
  real,    allocatable, dimension(:,:,:,:) :: P_i    !< Interface pressures [R L2 T-2 ~> Pa]
  real,    allocatable, dimension(:,:,:,:) :: dRdT_i !< dRho/dT [R degC-1 ~> kg m-3 degC-1] at top edge
  real,    allocatable, dimension(:,:,:,:) :: dRdS_i !< dRho/dS [R ppt-1 ~> kg m-3 ppt-1] at top edge
  integer, allocatable, dimension(:,:)     :: ns     !< Number of interfacs in a column
  logical, allocatable, dimension(:,:,:) :: stable_cell !< True if the cell is stably stratified wrt to the next cell
  real :: R_to_kg_m3 = 1.0                   !< A rescaling factor translating density to kg m-3 for
                                             !! use in diagnostic messages [kg m-3 R-1 ~> 1].
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.
  integer :: neutral_pos_method              !< Method to find the position of a neutral surface within the layer
  character(len=40)  :: delta_rho_form       !< Determine which (if any) approximation is made to the
                                             !! equation describing the difference in density

  integer :: id_uhEff_2d = -1 !< Diagnostic IDs
  integer :: id_vhEff_2d = -1 !< Diagnostic IDs

  type(EOS_type), pointer :: EOS => NULL()  !< Equation of state parameters
  type(remapping_CS) :: remap_CS   !< Remapping control structure used to create sublayers
  logical :: remap_answers_2018    !< If true, use the order of arithmetic and expressions that
                                   !! recover the answers for remapping from the end of 2018.
                                   !! Otherwise, use more robust forms of the same expressions.
  type(KPP_CS),           pointer :: KPP_CSp => NULL()          !< KPP control structure needed to get BLD
  type(energetic_PBL_CS), pointer :: energetic_PBL_CSp => NULL()!< ePBL control structure needed to get MLD
end type neutral_diffusion_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40)  :: mdl = "MOM_neutral_diffusion" !< module name

contains

!> Read parameters and allocate control structure for neutral_diffusion module.
logical function neutral_diffusion_init(Time, G, GV, US, param_file, diag, EOS, diabatic_CSp, CS)
  type(time_type), target,    intent(in)    :: Time       !< Time structure
  type(ocean_grid_type),      intent(in)    :: G          !< Grid structure
  type(verticalGrid_type),    intent(in)    :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in)    :: US         !< A dimensional unit scaling type
  type(diag_ctrl), target,    intent(inout) :: diag       !< Diagnostics control structure
  type(param_file_type),      intent(in)    :: param_file !< Parameter file structure
  type(EOS_type),  target,    intent(in)    :: EOS        !< Equation of state
  type(diabatic_CS),          pointer       :: diabatic_CSp!< KPP control structure needed to get BLD
  type(neutral_diffusion_CS), pointer       :: CS         !< Neutral diffusion control structure

  ! Local variables
  character(len=256) :: mesg    ! Message for error messages.
  character(len=80)  :: string  ! Temporary strings
  logical :: default_2018_answers
  logical :: boundary_extrap

  if (associated(CS)) then
    call MOM_error(FATAL, "neutral_diffusion_init called with associated control structure.")
    return
  endif

  ! Log this module and master switch for turning it on/off
  call get_param(param_file, mdl, "USE_NEUTRAL_DIFFUSION", neutral_diffusion_init, &
                 default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, &
           "This module implements neutral diffusion of tracers", &
           all_default=.not.neutral_diffusion_init)
  call get_param(param_file, mdl, "USE_NEUTRAL_DIFFUSION", neutral_diffusion_init, &
                 "If true, enables the neutral diffusion module.", &
                 default=.false.)

  if (.not.neutral_diffusion_init) return

  allocate(CS)
  CS%diag => diag
  CS%EOS => EOS
 ! call openParameterBlock(param_file,'NEUTRAL_DIFF')

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "NDIFF_CONTINUOUS", CS%continuous_reconstruction, &
                 "If true, uses a continuous reconstruction of T and S when "//&
                 "finding neutral surfaces along which diffusion will happen. "//&
                 "If false, a PPM discontinuous reconstruction of T and S "//&
                 "is done which results in a higher order routine but exacts "//&
                 "a higher computational cost.", default=.true.)
  call get_param(param_file, mdl, "NDIFF_REF_PRES", CS%ref_pres,                    &
                 "The reference pressure (Pa) used for the derivatives of "//&
                 "the equation of state. If negative (default), local pressure is used.", &
                 units="Pa", default = -1., scale=US%kg_m3_to_R*US%m_s_to_L_T**2)
  call get_param(param_file, mdl, "NDIFF_INTERIOR_ONLY", CS%interior_only, &
                 "If true, only applies neutral diffusion in the ocean interior."//&
                 "That is, the algorithm will exclude the surface and bottom"//&
                 "boundary layers.", default = .false.)
  call get_param(param_file, mdl, "NDIFF_USE_UNMASKED_TRANSPORT_BUG", CS%use_unmasked_transport_bug, &
                 "If true, use an older form for the accumulation of neutral-diffusion "//&
                 "transports that were unmasked, as used prior to Jan 2018. This is not "//&
                 "recommended.", default = .false.)

  ! Initialize and configure remapping
  if ( .not.CS%continuous_reconstruction ) then
    call get_param(param_file, mdl, "NDIFF_BOUNDARY_EXTRAP", boundary_extrap, &
                   "Extrapolate at the top and bottommost cells, otherwise   \n"//  &
                   "assume boundaries are piecewise constant",                      &
                   default=.false.)
    call get_param(param_file, mdl, "NDIFF_REMAPPING_SCHEME", string, &
                   "This sets the reconstruction scheme used "//&
                   "for vertical remapping for all variables. "//&
                   "It can be one of the following schemes: "//&
                   trim(remappingSchemesDoc), default=remappingDefaultScheme)
    call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
    call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", CS%remap_answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
    call initialize_remapping( CS%remap_CS, string, boundary_extrapolation=boundary_extrap, &
                               answers_2018=CS%remap_answers_2018 )
    call extract_member_remapping_CS(CS%remap_CS, degree=CS%deg)
    call get_param(param_file, mdl, "NEUTRAL_POS_METHOD", CS%neutral_pos_method,   &
                   "Method used to find the neutral position                 \n"// &
                   "1. Delta_rho varies linearly, find 0 crossing            \n"// &
                   "2. Alpha and beta vary linearly from top to bottom,      \n"// &
                   "   Newton's method for neutral position                  \n"// &
                   "3. Full nonlinear equation of state, use regula falsi    \n"// &
                   "   for neutral position", default=3)
    if (CS%neutral_pos_method > 4 .or. CS%neutral_pos_method < 0) then
      call MOM_error(FATAL,"Invalid option for NEUTRAL_POS_METHOD")
    endif

    call get_param(param_file, mdl, "DELTA_RHO_FORM", CS%delta_rho_form,           &
                   "Determine how the difference in density is calculated    \n"// &
                   "  full       : Difference of in-situ densities           \n"// &
                   "  no_pressure: Calculated from dRdT, dRdS, but no        \n"// &
                   "               pressure dependence",                           &
                   default="mid_pressure")
    if (CS%neutral_pos_method > 1) then
      call get_param(param_file, mdl, "NDIFF_DRHO_TOL", CS%drho_tol,            &
                     "Sets the convergence criterion for finding the neutral\n"// &
                     "position within a layer in kg m-3.",                        &
                     default=1.e-10, scale=US%kg_m3_to_R)
      call get_param(param_file, mdl, "NDIFF_X_TOL", CS%x_tol,            &
                     "Sets the convergence criterion for a change in nondim\n"// &
                     "position within a layer.",                        &
                     default=0.)
      call get_param(param_file, mdl, "NDIFF_MAX_ITER", CS%max_iter,              &
                    "The maximum number of iterations to be done before \n"//     &
                     "exiting the iterative loop to find the neutral surface",    &
                     default=10)
    endif
    call get_param(param_file, mdl, "NDIFF_DEBUG", CS%debug,             &
                   "Turns on verbose output for discontinuous neutral "//&
                   "diffusion routines.", &
                   default = .false.)
    call get_param(param_file, mdl, "HARD_FAIL_HEFF", CS%hard_fail_heff, &
                  "Bring down the model if a problem with heff is detected",&
                   default = .true.)
  endif

  if (CS%interior_only) then
    call extract_diabatic_member(diabatic_CSp, KPP_CSp=CS%KPP_CSp)
    call extract_diabatic_member(diabatic_CSp, energetic_PBL_CSp=CS%energetic_PBL_CSp)
    if ( .not. ASSOCIATED(CS%energetic_PBL_CSp) .and. .not. ASSOCIATED(CS%KPP_CSp) ) then
      call MOM_error(FATAL,"NDIFF_INTERIOR_ONLY is true, but no valid boundary layer scheme was found")
    endif
  endif
  ! Store a rescaling factor for use in diagnostic messages.
  CS%R_to_kg_m3 = US%R_to_kg_m3

! call get_param(param_file, mdl, "KHTR", CS%KhTr, &
!                "The background along-isopycnal tracer diffusivity.", &
!                units="m2 s-1", default=0.0)
!  call closeParameterBlock(param_file)
  if (CS%continuous_reconstruction) then
    CS%nsurf = 2*GV%ke+2 ! Continuous reconstruction means that every interface has two connections
    allocate(CS%dRdT(SZI_(G),SZJ_(G),SZK_(GV)+1)) ; CS%dRdT(:,:,:) = 0.
    allocate(CS%dRdS(SZI_(G),SZJ_(G),SZK_(GV)+1)) ; CS%dRdS(:,:,:) = 0.
  else
    CS%nsurf = 4*GV%ke   ! Discontinuous means that every interface has four connections
    allocate(CS%T_i(SZI_(G),SZJ_(G),SZK_(GV),2))    ; CS%T_i(:,:,:,:) = 0.
    allocate(CS%S_i(SZI_(G),SZJ_(G),SZK_(GV),2))    ; CS%S_i(:,:,:,:) = 0.
    allocate(CS%P_i(SZI_(G),SZJ_(G),SZK_(GV),2))    ; CS%P_i(:,:,:,:) = 0.
    allocate(CS%dRdT_i(SZI_(G),SZJ_(G),SZK_(GV),2)) ; CS%dRdT_i(:,:,:,:) = 0.
    allocate(CS%dRdS_i(SZI_(G),SZJ_(G),SZK_(GV),2)) ; CS%dRdS_i(:,:,:,:) = 0.
    allocate(CS%ppoly_coeffs_T(SZI_(G),SZJ_(G),SZK_(GV),CS%deg+1)) ; CS%ppoly_coeffs_T(:,:,:,:) = 0.
    allocate(CS%ppoly_coeffs_S(SZI_(G),SZJ_(G),SZK_(GV),CS%deg+1)) ; CS%ppoly_coeffs_S(:,:,:,:) = 0.
    allocate(CS%ns(SZI_(G),SZJ_(G)))    ; CS%ns(:,:) = 0.
  endif
  ! T-points
  allocate(CS%Tint(SZI_(G),SZJ_(G),SZK_(GV)+1)) ; CS%Tint(:,:,:) = 0.
  allocate(CS%Sint(SZI_(G),SZJ_(G),SZK_(GV)+1)) ; CS%Sint(:,:,:) = 0.
  allocate(CS%Pint(SZI_(G),SZJ_(G),SZK_(GV)+1)) ; CS%Pint(:,:,:) = 0.
  allocate(CS%stable_cell(SZI_(G),SZJ_(G),SZK_(GV))) ; CS%stable_cell(:,:,:) = .true.
  ! U-points
  allocate(CS%uPoL(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%uPoL(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0.
  allocate(CS%uPoR(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%uPoR(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0.
  allocate(CS%uKoL(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%uKoL(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0
  allocate(CS%uKoR(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%uKoR(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0
  allocate(CS%uHeff(G%isd:G%ied,G%jsd:G%jed,CS%nsurf-1)); CS%uHeff(G%isc-1:G%iec,G%jsc:G%jec,:) = 0
  ! V-points
  allocate(CS%vPoL(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%vPoL(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0.
  allocate(CS%vPoR(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%vPoR(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0.
  allocate(CS%vKoL(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%vKoL(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0
  allocate(CS%vKoR(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%vKoR(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0
  allocate(CS%vHeff(G%isd:G%ied,G%jsd:G%jed,CS%nsurf-1)); CS%vHeff(G%isc:G%iec,G%jsc-1:G%jec,:) = 0

end function neutral_diffusion_init

!> Calculate remapping factors for u/v columns used to map adjoining columns to
!! a shared coordinate space.
subroutine neutral_diffusion_calc_coeffs(G, GV, US, h, T, S, CS, p_surf)
  type(ocean_grid_type),                     intent(in) :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in) :: GV  !< ocean vertical grid structure
  type(unit_scale_type),                     intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h   !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: T   !< Potential temperature [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: S   !< Salinity [ppt]
  type(neutral_diffusion_CS),                pointer    :: CS  !< Neutral diffusion control structure
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: p_surf !< Surface pressure to include in pressures used
                                                              !! for equation of state calculations [R L2 T-2 ~> Pa]

  ! Local variables
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k
  ! Variables used for reconstructions
  real, dimension(SZK_(GV),2) :: ppoly_r_S      ! Reconstruction slopes
  real, dimension(SZI_(G), SZJ_(G)) :: hEff_sum ! Summed effective face thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G))  :: hbl      ! Boundary layer depth [H ~> m or kg m-2]
  integer :: iMethod
  real, dimension(SZI_(G)) :: ref_pres ! Reference pressure used to calculate alpha/beta [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G)) :: rho_tmp  ! Routine to calculate drho_dp, returns density which is not used
  real :: h_neglect, h_neglect_edge    ! Negligible thicknesses [H ~> m or kg m-2]
  integer, dimension(SZI_(G), SZJ_(G)) :: k_top  ! Index of the first layer within the boundary
  real,    dimension(SZI_(G), SZJ_(G)) :: zeta_top ! Distance from the top of a layer to the intersection of the
                                                   ! top extent of the boundary layer (0 at top, 1 at bottom) [nondim]
  integer, dimension(SZI_(G), SZJ_(G)) :: k_bot    ! Index of the last layer within the boundary
  real,    dimension(SZI_(G), SZJ_(G)) :: zeta_bot ! Distance of the lower layer to the boundary layer depth
  real :: pa_to_H                      ! A conversion factor from pressure to H units [H T2 R-1 Z-2 ~> m Pa-1 or s2 m-2]

  pa_to_H = 1. / (GV%H_to_RZ * GV%g_Earth)

  k_top(:,:) = 1     ; k_bot(:,:) = 1
  zeta_top(:,:) = 0. ; zeta_bot(:,:) = 0.

  ! Check if hbl needs to be extracted
  if (CS%interior_only) then
    if (ASSOCIATED(CS%KPP_CSp)) call KPP_get_BLD(CS%KPP_CSp, hbl, G, US, m_to_BLD_units=GV%m_to_H)
    if (ASSOCIATED(CS%energetic_PBL_CSp)) call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, hbl, G, US, &
                                                                     m_to_MLD_units=GV%m_to_H)
    call pass_var(hbl,G%Domain)
    ! get k-indices and zeta
    do j=G%jsc-1, G%jec+1 ; do i=G%isc-1,G%iec+1
      call boundary_k_range(SURFACE, GV%ke, h(i,j,:), hbl(i,j), k_top(i,j), zeta_top(i,j), k_bot(i,j), zeta_bot(i,j))
    enddo ; enddo
    ! TODO: add similar code for BOTTOM boundary layer
  endif

  if (.not.CS%remap_answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  ! If doing along isopycnal diffusion (as opposed to neutral diffusion, set the reference pressure)
  if (CS%ref_pres>=0.) then
    ref_pres(:) = CS%ref_pres
  endif

  if (CS%continuous_reconstruction) then
    CS%dRdT(:,:,:) = 0.
    CS%dRdS(:,:,:) = 0.
  else
    CS%T_i(:,:,:,:) = 0.
    CS%S_i(:,:,:,:) = 0.
    CS%dRdT_i(:,:,:,:) = 0.
    CS%dRdS_i(:,:,:,:) = 0.
    CS%ns(:,:) = 0.
    CS%stable_cell(:,:,:) = .true.
  endif

  ! Calculate pressure at interfaces and layer averaged alpha/beta
  if (present(p_surf)) then
    do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      CS%Pint(i,j,1) = p_surf(i,j)
    enddo ; enddo
  else
    CS%Pint(:,:,1) = 0.
  endif
  do k=1,GV%ke ; do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
    CS%Pint(i,j,k+1) = CS%Pint(i,j,k) + h(i,j,k)*(GV%g_Earth*GV%H_to_RZ)
  enddo ; enddo ; enddo

  ! Pressures at the interfaces, this is redundant as P_i(k,1) = P_i(k-1,2) however retain this
  ! for now to ensure consitency of indexing for diiscontinuous reconstructions
  if (.not. CS%continuous_reconstruction) then
    if (present(p_surf)) then
      do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
        CS%P_i(i,j,1,1) = p_surf(i,j)
        CS%P_i(i,j,1,2) = p_surf(i,j) + h(i,j,1)*(GV%H_to_RZ*GV%g_Earth)
      enddo ; enddo
    else
      do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
        CS%P_i(i,j,1,1) = 0.
        CS%P_i(i,j,1,2) = h(i,j,1)*(GV%H_to_RZ*GV%g_Earth)
      enddo ; enddo
    endif
    do k=2,GV%ke ; do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      CS%P_i(i,j,k,1) = CS%P_i(i,j,k-1,2)
      CS%P_i(i,j,k,2) = CS%P_i(i,j,k-1,2) + h(i,j,k)*(GV%H_to_RZ*GV%g_Earth)
    enddo ; enddo ; enddo
  endif

  EOSdom(:) = EOS_domain(G%HI, halo=1)
  do j = G%jsc-1, G%jec+1
    ! Interpolate state to interface
    do i = G%isc-1, G%iec+1
      if (CS%continuous_reconstruction) then
        call interface_scalar(GV%ke, h(i,j,:), T(i,j,:), CS%Tint(i,j,:), 2, h_neglect)
        call interface_scalar(GV%ke, h(i,j,:), S(i,j,:), CS%Sint(i,j,:), 2, h_neglect)
      else
        call build_reconstructions_1d( CS%remap_CS, GV%ke, h(i,j,:), T(i,j,:), CS%ppoly_coeffs_T(i,j,:,:), &
                                       CS%T_i(i,j,:,:), ppoly_r_S, iMethod, h_neglect, h_neglect_edge )
        call build_reconstructions_1d( CS%remap_CS, GV%ke, h(i,j,:), S(i,j,:), CS%ppoly_coeffs_S(i,j,:,:), &
                                       CS%S_i(i,j,:,:), ppoly_r_S, iMethod, h_neglect, h_neglect_edge )
        ! In the current ALE formulation, interface values are not exactly at the 0. or 1. of the
        ! polynomial reconstructions
        do k=1,GV%ke
           CS%T_i(i,j,k,1) = evaluation_polynomial( CS%ppoly_coeffs_T(i,j,k,:), CS%deg+1, 0. )
           CS%T_i(i,j,k,2) = evaluation_polynomial( CS%ppoly_coeffs_T(i,j,k,:), CS%deg+1, 1. )
           CS%S_i(i,j,k,1) = evaluation_polynomial( CS%ppoly_coeffs_S(i,j,k,:), CS%deg+1, 0. )
           CS%S_i(i,j,k,2) = evaluation_polynomial( CS%ppoly_coeffs_S(i,j,k,:), CS%deg+1, 1. )
        enddo
      endif
    enddo

    ! Continuous reconstruction
    if (CS%continuous_reconstruction) then
      do k = 1, GV%ke+1
        if (CS%ref_pres<0) ref_pres(:) = CS%Pint(:,j,k)
        call calculate_density_derivs(CS%Tint(:,j,k), CS%Sint(:,j,k), ref_pres, CS%dRdT(:,j,k), &
                                      CS%dRdS(:,j,k), CS%EOS, EOSdom)
      enddo
    else ! Discontinuous reconstruction
      do k = 1, GV%ke
        if (CS%ref_pres<0) ref_pres(:) = CS%Pint(:,j,k)
        ! Calculate derivatives for the top interface
        call calculate_density_derivs(CS%T_i(:,j,k,1), CS%S_i(:,j,k,1), ref_pres, CS%dRdT_i(:,j,k,1), &
                                      CS%dRdS_i(:,j,k,1), CS%EOS, EOSdom)
        if (CS%ref_pres<0) ref_pres(:) = CS%Pint(:,j,k+1)
        ! Calculate derivatives at the bottom interface
        call calculate_density_derivs(CS%T_i(:,j,k,2), CS%S_i(:,j,k,2), ref_pres, CS%dRdT_i(:,j,k,2), &
                                      CS%dRdS_i(:,j,k,2), CS%EOS, EOSdom)
      enddo
    endif
  enddo

  if (.not. CS%continuous_reconstruction) then
    do j = G%jsc-1, G%jec+1 ; do i = G%isc-1, G%iec+1
      call mark_unstable_cells( CS, GV%ke, CS%T_i(i,j,:,:), CS%S_i(i,j,:,:), CS%P_i(i,j,:,:), CS%stable_cell(i,j,:) )
      if (CS%interior_only) then
        if (.not. CS%stable_cell(i,j,k_bot(i,j))) zeta_bot(i,j) = -1.
        ! set values in the surface and bottom boundary layer to false.
        do k = 1, k_bot(i,j)
          CS%stable_cell(i,j,k) = .false.
        enddo
      endif
    enddo ; enddo
  endif

  CS%uhEff(:,:,:) = 0.
  CS%vhEff(:,:,:) = 0.
  CS%uPoL(:,:,:) = 0.
  CS%vPoL(:,:,:) = 0.
  CS%uPoR(:,:,:) = 0.
  CS%vPoR(:,:,:) = 0.
  CS%uKoL(:,:,:) = 1
  CS%vKoL(:,:,:) = 1
  CS%uKoR(:,:,:) = 1
  CS%vKoR(:,:,:) = 1

  ! Neutral surface factors at U points
  do j = G%jsc, G%jec ; do I = G%isc-1, G%iec
    if (G%mask2dCu(I,j) > 0.) then
      if (CS%continuous_reconstruction) then
        call find_neutral_surface_positions_continuous(GV%ke,                                    &
                CS%Pint(i,j,:), CS%Tint(i,j,:), CS%Sint(i,j,:), CS%dRdT(i,j,:), CS%dRdS(i,j,:),            &
                CS%Pint(i+1,j,:), CS%Tint(i+1,j,:), CS%Sint(i+1,j,:), CS%dRdT(i+1,j,:), CS%dRdS(i+1,j,:),  &
                CS%uPoL(I,j,:), CS%uPoR(I,j,:), CS%uKoL(I,j,:), CS%uKoR(I,j,:), CS%uhEff(I,j,:),           &
                k_bot(I,j), k_bot(I+1,j), zeta_bot(I,j), zeta_bot(I+1,j))
      else
        call find_neutral_surface_positions_discontinuous(CS, GV%ke, &
            CS%P_i(i,j,:,:), h(i,j,:), CS%T_i(i,j,:,:), CS%S_i(i,j,:,:), CS%ppoly_coeffs_T(i,j,:,:),           &
            CS%ppoly_coeffs_S(i,j,:,:),CS%stable_cell(i,j,:),                                                  &
            CS%P_i(i+1,j,:,:), h(i+1,j,:), CS%T_i(i+1,j,:,:), CS%S_i(i+1,j,:,:), CS%ppoly_coeffs_T(i+1,j,:,:), &
            CS%ppoly_coeffs_S(i+1,j,:,:), CS%stable_cell(i+1,j,:),                                             &
            CS%uPoL(I,j,:), CS%uPoR(I,j,:), CS%uKoL(I,j,:), CS%uKoR(I,j,:), CS%uhEff(I,j,:),                   &
            hard_fail_heff = CS%hard_fail_heff)
      endif
    endif
  enddo ; enddo

  ! Neutral surface factors at V points
  do J = G%jsc-1, G%jec ; do i = G%isc, G%iec
    if (G%mask2dCv(i,J) > 0.) then
      if (CS%continuous_reconstruction) then
        call find_neutral_surface_positions_continuous(GV%ke,                                              &
                CS%Pint(i,j,:), CS%Tint(i,j,:), CS%Sint(i,j,:), CS%dRdT(i,j,:), CS%dRdS(i,j,:),           &
                CS%Pint(i,j+1,:), CS%Tint(i,j+1,:), CS%Sint(i,j+1,:), CS%dRdT(i,j+1,:), CS%dRdS(i,j+1,:), &
                CS%vPoL(i,J,:), CS%vPoR(i,J,:), CS%vKoL(i,J,:), CS%vKoR(i,J,:), CS%vhEff(i,J,:),          &
                k_bot(i,J), k_bot(i,J+1), zeta_bot(i,J), zeta_bot(i,J+1))
      else
        call find_neutral_surface_positions_discontinuous(CS, GV%ke, &
            CS%P_i(i,j,:,:), h(i,j,:), CS%T_i(i,j,:,:), CS%S_i(i,j,:,:), CS%ppoly_coeffs_T(i,j,:,:),           &
            CS%ppoly_coeffs_S(i,j,:,:),CS%stable_cell(i,j,:),                                                  &
            CS%P_i(i,j+1,:,:), h(i,j+1,:), CS%T_i(i,j+1,:,:), CS%S_i(i,j+1,:,:), CS%ppoly_coeffs_T(i,j+1,:,:), &
            CS%ppoly_coeffs_S(i,j+1,:,:), CS%stable_cell(i,j+1,:),                                             &
            CS%vPoL(I,j,:), CS%vPoR(I,j,:), CS%vKoL(I,j,:), CS%vKoR(I,j,:), CS%vhEff(I,j,:),                   &
            hard_fail_heff = CS%hard_fail_heff)
      endif
    endif
  enddo ; enddo

  ! Continuous reconstructions calculate hEff as the difference between the pressures of the
  ! neutral surfaces which need to be reconverted to thickness units. The discontinuous version
  ! calculates hEff from the nondimensional fraction of the layer spanned by adjacent neutral
  ! surfaces, so hEff is already in thickness units.
  if (CS%continuous_reconstruction) then
    if (CS%use_unmasked_transport_bug) then
      ! This option is not recommended but needed to recover answers prior to Jan 2018.
      ! It is independent of the other 2018 answers flags.
      do k = 1, CS%nsurf-1 ; do j = G%jsc, G%jec ; do I = G%isc-1, G%iec
        CS%uhEff(I,j,k) = CS%uhEff(I,j,k) / GV%H_to_pa
      enddo ; enddo ; enddo
      do k = 1, CS%nsurf-1 ; do J = G%jsc-1, G%jec ; do i = G%isc, G%iec
        CS%vhEff(I,j,k) = CS%vhEff(I,j,k) / GV%H_to_pa
      enddo ; enddo ; enddo
    else
      do k = 1, CS%nsurf-1 ; do j = G%jsc, G%jec ; do I = G%isc-1, G%iec
        if (G%mask2dCu(I,j) > 0.) CS%uhEff(I,j,k) = CS%uhEff(I,j,k) * pa_to_H
      enddo ; enddo ; enddo
      do k = 1, CS%nsurf-1 ; do J = G%jsc-1, G%jec ; do i = G%isc, G%iec
        if (G%mask2dCv(i,J) > 0.) CS%vhEff(i,J,k) = CS%vhEff(i,J,k) * pa_to_H
      enddo ; enddo ; enddo
    endif
  endif

  if (CS%id_uhEff_2d>0) then
    hEff_sum(:,:) = 0.
    do k = 1,CS%nsurf-1 ; do j=G%jsc,G%jec ; do i=G%isc-1,G%iec
      hEff_sum(i,j) = hEff_sum(i,j) + CS%uhEff(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_uhEff_2d, hEff_sum, CS%diag)
  endif
  if (CS%id_vhEff_2d>0) then
    hEff_sum(:,:) = 0.
    do k = 1,CS%nsurf-1 ; do j=G%jsc-1,G%jec ; do i=G%isc,G%iec
      hEff_sum(i,j) = hEff_sum(i,j) + CS%vhEff(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_vhEff_2d, hEff_sum, CS%diag)
  endif

end subroutine neutral_diffusion_calc_coeffs

!> Update tracer concentration due to neutral diffusion; layer thickness unchanged by this update.
subroutine neutral_diffusion(G, GV, h, Coef_x, Coef_y, dt, Reg, US, CS)
  type(ocean_grid_type),                     intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: Coef_x !< dt * Kh * dy / dx at u-points [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: Coef_y !< dt * Kh * dx / dy at v-points [L2 ~> m2]
  real,                                      intent(in)    :: dt     !< Tracer time step * I_numitts [T ~> s]
                                                                     !! (I_numitts in tracer_hordiff)
  type(tracer_registry_type),                pointer       :: Reg    !< Tracer registry
  type(unit_scale_type),                     intent(in)    :: US     !< A dimensional unit scaling type
  type(neutral_diffusion_CS),                pointer       :: CS     !< Neutral diffusion control structure

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),CS%nsurf-1) :: uFlx        ! Zonal flux of tracer [H conc ~> m conc or conc kg m-2]
  real, dimension(SZI_(G),SZJB_(G),CS%nsurf-1) :: vFlx        ! Meridional flux of tracer
                                                              ! [H conc ~> m conc or conc kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))    :: tendency    ! tendency array for diagnostics
                                                              ! [H conc T-1 ~> m conc s-1 or kg m-2 conc s-1]
  real, dimension(SZI_(G),SZJ_(G))             :: tendency_2d ! depth integrated content tendency for diagn
                                                              ! [H conc T-1 ~> m conc s-1 or kg m-2 conc s-1]
  real, dimension(SZIB_(G),SZJ_(G))            :: trans_x_2d  ! depth integrated diffusive tracer x-transport diagn
  real, dimension(SZI_(G),SZJB_(G))            :: trans_y_2d  ! depth integrated diffusive tracer y-transport diagn
  real, dimension(SZK_(GV))                    :: dTracer     ! change in tracer concentration due to ndiffusion
                                                              ! [H L2 conc ~> m3 conc or kg conc]

  type(tracer_type), pointer                   :: Tracer => NULL() ! Pointer to the current tracer

  integer :: i, j, k, m, ks, nk
  real :: Idt  ! The inverse of the time step [T-1 ~> s-1]
  real :: h_neglect, h_neglect_edge

  if (.not.CS%remap_answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  else
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  endif

  nk = GV%ke

  do m = 1,Reg%ntr ! Loop over tracer registry

    tracer => Reg%Tr(m)

    ! for diagnostics
    if (tracer%id_dfxy_conc > 0 .or. tracer%id_dfxy_cont > 0 .or. tracer%id_dfxy_cont_2d > 0 .or. &
        tracer%id_dfx_2d > 0 .or. tracer%id_dfy_2d > 0) then
      Idt = 1.0 / dt
      tendency(:,:,:)  = 0.0
    endif

    uFlx(:,:,:) = 0.
    vFlx(:,:,:) = 0.

    ! x-flux
    do j = G%jsc,G%jec ; do I = G%isc-1,G%iec
      if (G%mask2dCu(I,j)>0.) then
        call neutral_surface_flux(nk, CS%nsurf, CS%deg, h(i,j,:), h(i+1,j,:),       &
                                  tracer%t(i,j,:), tracer%t(i+1,j,:), &
                                  CS%uPoL(I,j,:), CS%uPoR(I,j,:), &
                                  CS%uKoL(I,j,:), CS%uKoR(I,j,:), &
                                  CS%uhEff(I,j,:), uFlx(I,j,:), &
                                  CS%continuous_reconstruction, h_neglect, CS%remap_CS, h_neglect_edge)
      endif
    enddo ; enddo

    ! y-flux
    do J = G%jsc-1,G%jec ; do i = G%isc,G%iec
      if (G%mask2dCv(i,J)>0.) then
        call neutral_surface_flux(nk, CS%nsurf, CS%deg, h(i,j,:), h(i,j+1,:),       &
                                  tracer%t(i,j,:), tracer%t(i,j+1,:), &
                                  CS%vPoL(i,J,:), CS%vPoR(i,J,:), &
                                  CS%vKoL(i,J,:), CS%vKoR(i,J,:), &
                                  CS%vhEff(i,J,:), vFlx(i,J,:),   &
                                  CS%continuous_reconstruction, h_neglect, CS%remap_CS, h_neglect_edge)
      endif
    enddo ; enddo

    ! Update the tracer concentration from divergence of neutral diffusive flux components
    do j = G%jsc,G%jec ; do i = G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then

        dTracer(:) = 0.
        do ks = 1,CS%nsurf-1
          k = CS%uKoL(I,j,ks)
          dTracer(k) = dTracer(k) + Coef_x(I,j)   * uFlx(I,j,ks)
          k = CS%uKoR(I-1,j,ks)
          dTracer(k) = dTracer(k) - Coef_x(I-1,j) * uFlx(I-1,j,ks)
          k = CS%vKoL(i,J,ks)
          dTracer(k) = dTracer(k) + Coef_y(i,J)   * vFlx(i,J,ks)
          k = CS%vKoR(i,J-1,ks)
          dTracer(k) = dTracer(k) - Coef_y(i,J-1) * vFlx(i,J-1,ks)
        enddo
        do k = 1, GV%ke
          tracer%t(i,j,k) = tracer%t(i,j,k) + dTracer(k) * &
                          ( G%IareaT(i,j) / ( h(i,j,k) + GV%H_subroundoff ) )
        enddo

        if (tracer%id_dfxy_conc > 0  .or. tracer%id_dfxy_cont > 0 .or. tracer%id_dfxy_cont_2d > 0 ) then
          do k = 1, GV%ke
            tendency(i,j,k) = dTracer(k) * G%IareaT(i,j) * Idt
          enddo
        endif

      endif
    enddo ; enddo

    ! Diagnose vertically summed zonal flux, giving zonal tracer transport from ndiff.
    ! Note sign corresponds to downgradient flux convention.
    if (tracer%id_dfx_2d > 0) then
      do j = G%jsc,G%jec ; do I = G%isc-1,G%iec
        trans_x_2d(I,j) = 0.
        if (G%mask2dCu(I,j)>0.) then
          do ks = 1,CS%nsurf-1
            trans_x_2d(I,j) = trans_x_2d(I,j) - Coef_x(I,j) * uFlx(I,j,ks)
          enddo
          trans_x_2d(I,j) = trans_x_2d(I,j) * Idt
        endif
      enddo ; enddo
      call post_data(tracer%id_dfx_2d, trans_x_2d(:,:), CS%diag)
    endif

    ! Diagnose vertically summed merid flux, giving meridional tracer transport from ndiff.
    ! Note sign corresponds to downgradient flux convention.
    if (tracer%id_dfy_2d > 0) then
      do J = G%jsc-1,G%jec ; do i = G%isc,G%iec
        trans_y_2d(i,J) = 0.
        if (G%mask2dCv(i,J)>0.) then
          do ks = 1,CS%nsurf-1
            trans_y_2d(i,J) = trans_y_2d(i,J) - Coef_y(i,J) * vFlx(i,J,ks)
          enddo
          trans_y_2d(i,J) = trans_y_2d(i,J) * Idt
        endif
      enddo ; enddo
      call post_data(tracer%id_dfy_2d, trans_y_2d(:,:), CS%diag)
    endif

    ! post tendency of layer-integrated tracer content
    if (tracer%id_dfxy_cont > 0) then
      call post_data(tracer%id_dfxy_cont, tendency(:,:,:), CS%diag)
    endif

    ! post depth summed tendency for tracer content
    if (tracer%id_dfxy_cont_2d > 0) then
      tendency_2d(:,:) = 0.
      do j = G%jsc,G%jec ; do i = G%isc,G%iec
        do k = 1, GV%ke
          tendency_2d(i,j) = tendency_2d(i,j) + tendency(i,j,k)
        enddo
      enddo ; enddo
      call post_data(tracer%id_dfxy_cont_2d, tendency_2d(:,:), CS%diag)
    endif

    ! post tendency of tracer concentration; this step must be
    ! done after posting tracer content tendency, since we alter
    ! the tendency array.
    if (tracer%id_dfxy_conc > 0) then
      do k = 1, GV%ke ; do j = G%jsc,G%jec ; do i = G%isc,G%iec
        tendency(i,j,k) =  tendency(i,j,k) / ( h(i,j,k) + GV%H_subroundoff )
      enddo ; enddo ; enddo
      call post_data(tracer%id_dfxy_conc, tendency, CS%diag)
    endif
  enddo ! Loop over tracer registry

end subroutine neutral_diffusion

!> Returns interface scalar, Si, for a column of layer values, S.
subroutine interface_scalar(nk, h, S, Si, i_method, h_neglect)
  integer,               intent(in)    :: nk       !< Number of levels
  real, dimension(nk),   intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2]
  real, dimension(nk),   intent(in)    :: S        !< Layer scalar (conc, e.g. ppt)
  real, dimension(nk+1), intent(inout) :: Si       !< Interface scalar (conc, e.g. ppt)
  integer,               intent(in)    :: i_method !< =1 use average of PLM edges
                                                   !! =2 use continuous PPM edge interpolation
  real,                  intent(in)    :: h_neglect !< A negligibly small thickness [H ~> m or kg m-2]
  ! Local variables
  integer :: k, km2, kp1
  real, dimension(nk) :: diff
  real :: Sb, Sa

  call PLM_diff(nk, h, S, 2, 1, diff)
  Si(1) = S(1) - 0.5 * diff(1)
  if (i_method==1) then
    do k = 2, nk
      ! Average of the two edge values (will be bounded and,
      ! when slopes are unlimited, notionally second-order accurate)
      Sa = S(k-1) + 0.5 * diff(k-1) ! Lower edge value of a PLM reconstruction for layer above
      Sb = S(k) - 0.5 * diff(k) ! Upper edge value of a PLM reconstruction for layer below
      Si(k) = 0.5 * ( Sa + Sb )
    enddo
  elseif (i_method==2) then
    do k = 2, nk
      ! PPM quasi-fourth order interpolation for edge values following
      ! equation 1.6 in Colella & Woodward, 1984: JCP 54, 174-201.
      km2 = max(1, k-2)
      kp1 = min(nk, k+1)
      Si(k) = ppm_edge(h(km2), h(k-1), h(k), h(kp1),  S(k-1), S(k), diff(k-1), diff(k), h_neglect)
    enddo
  endif
  Si(nk+1) = S(nk) + 0.5 * diff(nk)

end subroutine interface_scalar

!> Returns the PPM quasi-fourth order edge value at k+1/2 following
!! equation 1.6 in Colella & Woodward, 1984: JCP 54, 174-201.
real function ppm_edge(hkm1, hk, hkp1, hkp2,  Ak, Akp1, Pk, Pkp1, h_neglect)
  real, intent(in) :: hkm1 !< Width of cell k-1
  real, intent(in) :: hk   !< Width of cell k
  real, intent(in) :: hkp1 !< Width of cell k+1
  real, intent(in) :: hkp2 !< Width of cell k+2
  real, intent(in) :: Ak   !< Average scalar value of cell k
  real, intent(in) :: Akp1 !< Average scalar value of cell k+1
  real, intent(in) :: Pk   !< PLM slope for cell k
  real, intent(in) :: Pkp1 !< PLM slope for cell k+1
  real, intent(in) :: h_neglect !< A negligibly small thickness [H ~> m or kg m-2]

  ! Local variables
  real :: R_hk_hkp1, R_2hk_hkp1, R_hk_2hkp1, f1, f2, f3, f4

  R_hk_hkp1 = hk + hkp1
  if (R_hk_hkp1 <= 0.) then
    ppm_edge = 0.5 * ( Ak + Akp1 )
    return
  endif
  R_hk_hkp1 = 1. / R_hk_hkp1
  if (hk<hkp1) then
    ppm_edge = Ak + ( hk * R_hk_hkp1 ) * ( Akp1 - Ak )
  else
    ppm_edge = Akp1 + ( hkp1 * R_hk_hkp1 ) * ( Ak - Akp1 )
  endif

  R_2hk_hkp1 = 1. / ( ( 2. * hk + hkp1 ) + h_neglect )
  R_hk_2hkp1 = 1. / ( ( hk + 2. * hkp1 ) + h_neglect )
  f1 = 1./ ( ( hk + hkp1) + ( hkm1 + hkp2 ) )
  f2 = 2. * ( hkp1 * hk ) * R_hk_hkp1 * &
            ( ( hkm1 + hk ) * R_2hk_hkp1  - ( hkp2 + hkp1 ) * R_hk_2hkp1 )
  f3 = hk * ( hkm1 + hk ) * R_2hk_hkp1
  f4 = hkp1 * ( hkp1 + hkp2 ) * R_hk_2hkp1

  ppm_edge = ppm_edge + f1 * ( f2 * ( Akp1 - Ak ) - ( f3 * Pkp1 - f4 * Pk ) )

end function ppm_edge

!> Returns the average of a PPM reconstruction between two
!! fractional positions.
real function ppm_ave(xL, xR, aL, aR, aMean)
  real, intent(in) :: xL    !< Fraction position of left bound (0,1)
  real, intent(in) :: xR    !< Fraction position of right bound (0,1)
  real, intent(in) :: aL    !< Left edge scalar value, at x=0
  real, intent(in) :: aR    !< Right edge scalar value, at x=1
  real, intent(in) :: aMean !< Average scalar value of cell

  ! Local variables
  real :: dx, xave, a6, a6o3

  dx = xR - xL
  xave = 0.5 * ( xR + xL )
  a6o3 = 2. * aMean - ( aL + aR ) ! a6 / 3.
  a6 = 3. * a6o3

  if (dx<0.) then
    stop 'ppm_ave: dx<0 should not happend!'
  elseif (dx>1.) then
    stop 'ppm_ave: dx>1 should not happend!'
  elseif (dx==0.) then
    ppm_ave = aL + ( aR - aL ) * xR + a6 * xR * ( 1. - xR )
  else
    ppm_ave = ( aL + xave * ( ( aR - aL ) + a6 ) )  - a6o3 * ( xR**2 + xR * xL + xL**2 )
  endif
end function ppm_ave

!> A true signum function that returns either -abs(a), when x<0; or abs(a) when x>0; or 0 when x=0.
real function signum(a,x)
  real, intent(in) :: a !< The magnitude argument
  real, intent(in) :: x !< The sign (or zero) argument

  signum = sign(a,x)
  if (x==0.) signum = 0.

end function signum

!> Returns PLM slopes for a column where the slopes are the difference in value across each cell.
!! The limiting follows equation 1.8 in Colella & Woodward, 1984: JCP 54, 174-201.
subroutine PLM_diff(nk, h, S, c_method, b_method, diff)
  integer,             intent(in)    :: nk       !< Number of levels
  real, dimension(nk), intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2]
  real, dimension(nk), intent(in)    :: S        !< Layer salinity (conc, e.g. ppt)
  integer,             intent(in)    :: c_method !< Method to use for the centered difference
  integer,             intent(in)    :: b_method !< =1, use PCM in first/last cell, =2 uses linear extrapolation
  real, dimension(nk), intent(inout) :: diff     !< Scalar difference across layer (conc, e.g. ppt)
                                                 !! determined by the following values for c_method:
                                                 !!   1. Second order finite difference (not recommended)
                                                 !!   2. Second order finite volume (used in original PPM)
                                                 !!   3. Finite-volume weighted least squares linear fit
                                                 !! \todo  The use of c_method to choose a scheme is inefficient
                                                 !! and should eventually be moved up the call tree.

  ! Local variables
  integer :: k
  real :: hkm1, hk, hkp1, Skm1, Sk, Skp1, diff_l, diff_r, diff_c

  do k = 2, nk-1
    hkm1 = h(k-1)
    hk = h(k)
    hkp1 = h(k+1)

    if ( ( hkp1 + hk ) * ( hkm1 + hk ) > 0.) then
      Skm1 = S(k-1)
      Sk = S(k)
      Skp1 = S(k+1)
      if (c_method==1) then
        ! Simple centered diff (from White)
        if ( hk + 0.5 * (hkm1 + hkp1) /= 0. ) then
          diff_c = ( Skp1 - Skm1 ) * ( hk / ( hk + 0.5 * (hkm1 + hkp1) ) )
        else
          diff_c = 0.
        endif
      elseif (c_method==2) then
        ! Second order accurate centered FV slope (from Colella and Woodward, JCP 1984)
        diff_c = fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
      elseif (c_method==3) then
        ! Second order accurate finite-volume least squares slope
        diff_c = hk * fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
      endif
      ! Limit centered slope by twice the side differenced slopes
      diff_l = 2. * ( Sk - Skm1 )
      diff_r = 2. * ( Skp1 - Sk )
      if ( signum(1., diff_l) * signum(1., diff_r) <= 0. ) then
        diff(k) = 0. ! PCM for local extrema
      else
        diff(k) = sign( min( abs(diff_l), abs(diff_c), abs(diff_r) ), diff_c )
      endif
    else
      diff(k) = 0. ! PCM next to vanished layers
    endif
  enddo
  if (b_method==1) then ! PCM for top and bottom layer
    diff(1) = 0.
    diff(nk) = 0.
  elseif (b_method==2) then ! Linear extrapolation for top and bottom interfaces
    diff(1) = ( S(2) - S(1) ) * 2. * ( h(1) / ( h(1) + h(2) ) )
    diff(nk) = S(nk) - S(nk-1) * 2. * ( h(nk) / ( h(nk-1) + h(nk) ) )
  endif

end subroutine PLM_diff

!> Returns the cell-centered second-order finite volume (unlimited PLM) slope
!! using three consecutive cell widths and average values. Slope is returned
!! as a difference across the central cell (i.e. units of scalar S).
!! Discretization follows equation 1.7 in Colella & Woodward, 1984: JCP 54, 174-201.
real function fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  real, intent(in) :: hkm1 !< Left cell width
  real, intent(in) :: hk   !< Center cell width
  real, intent(in) :: hkp1 !< Right cell width
  real, intent(in) :: Skm1 !< Left cell average value
  real, intent(in) :: Sk   !< Center cell average value
  real, intent(in) :: Skp1 !< Right cell average value

  ! Local variables
  real :: h_sum, hp, hm

  h_sum = ( hkm1 + hkp1 ) + hk
  if (h_sum /= 0.) h_sum = 1./ h_sum
  hm =  hkm1 + hk
  if (hm /= 0.) hm = 1./ hm
  hp =  hkp1 + hk
  if (hp /= 0.) hp = 1./ hp
  fv_diff = ( hk * h_sum ) * &
            (   ( 2. * hkm1 + hk ) * hp * ( Skp1 - Sk ) &
              + ( 2. * hkp1 + hk ) * hm * ( Sk - Skm1 ) )
end function fv_diff


!> Returns the cell-centered second-order weighted least squares slope
!! using three consecutive cell widths and average values. Slope is returned
!! as a gradient (i.e. units of scalar S over width units).
real function fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  real, intent(in) :: hkm1 !< Left cell width
  real, intent(in) :: hk   !< Center cell width
  real, intent(in) :: hkp1 !< Right cell width
  real, intent(in) :: Skm1 !< Left cell average value
  real, intent(in) :: Sk   !< Center cell average value
  real, intent(in) :: Skp1 !< Right cell average value

  ! Local variables
  real :: xkm1, xkp1
  real :: h_sum, hx_sum, hxsq_sum, hxy_sum, hy_sum, det

  xkm1 = -0.5 * ( hk + hkm1 )
  xkp1 = 0.5 * ( hk + hkp1 )
  h_sum = ( hkm1 + hkp1 ) + hk
  hx_sum = hkm1*xkm1 + hkp1*xkp1
  hxsq_sum = hkm1*(xkm1**2) + hkp1*(xkp1**2)
  hxy_sum = hkm1*xkm1*Skm1 + hkp1*xkp1*Skp1
  hy_sum = ( hkm1*Skm1 + hkp1*Skp1 ) + hk*Sk
  det = h_sum * hxsq_sum - hx_sum**2
  if (det /= 0.) then
    !a = ( hxsq_sum * hy_sum - hx_sum*hxy_sum ) / det ! a would be mean of straight line fit
    fvlsq_slope = ( h_sum * hxy_sum - hx_sum*hy_sum ) / det ! Gradient of straight line fit
  else
    fvlsq_slope = 0. ! Adcroft's reciprocal rule
  endif
end function fvlsq_slope


!> Returns positions within left/right columns of combined interfaces using continuous reconstructions of T/S
subroutine find_neutral_surface_positions_continuous(nk, Pl, Tl, Sl, dRdTl, dRdSl, Pr, Tr, Sr, &
                                                     dRdTr, dRdSr, PoL, PoR, KoL, KoR, hEff, bl_kl, bl_kr, bl_zl, bl_zr)
  integer,                    intent(in)    :: nk    !< Number of levels
  real, dimension(nk+1),      intent(in)    :: Pl    !< Left-column interface pressure [R L2 T-2 ~> Pa] or other units
  real, dimension(nk+1),      intent(in)    :: Tl    !< Left-column interface potential temperature [degC]
  real, dimension(nk+1),      intent(in)    :: Sl    !< Left-column interface salinity [ppt]
  real, dimension(nk+1),      intent(in)    :: dRdTl !< Left-column dRho/dT [R degC-1 ~> kg m-3 degC-1]
  real, dimension(nk+1),      intent(in)    :: dRdSl !< Left-column dRho/dS [R ppt-1 ~> kg m-3 ppt-1]
  real, dimension(nk+1),      intent(in)    :: Pr    !< Right-column interface pressure [R L2 T-2 ~> Pa] or other units
  real, dimension(nk+1),      intent(in)    :: Tr    !< Right-column interface potential temperature [degC]
  real, dimension(nk+1),      intent(in)    :: Sr    !< Right-column interface salinity [ppt]
  real, dimension(nk+1),      intent(in)    :: dRdTr !< Left-column dRho/dT [R degC-1 ~> kg m-3 degC-1]
  real, dimension(nk+1),      intent(in)    :: dRdSr !< Left-column dRho/dS [R ppt-1 ~> kg m-3 ppt-1]
  real, dimension(2*nk+2),    intent(inout) :: PoL   !< Fractional position of neutral surface within
                                                     !! layer KoL of left column
  real, dimension(2*nk+2),    intent(inout) :: PoR   !< Fractional position of neutral surface within
                                                     !! layer KoR of right column
  integer, dimension(2*nk+2), intent(inout) :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(2*nk+2), intent(inout) :: KoR   !< Index of first right interface above neutral surface
  real, dimension(2*nk+1),    intent(inout) :: hEff  !< Effective thickness between two neutral surfaces
                                                     !! [R L2 T-2 ~> Pa] or other units following Pl and Pr.
  integer, optional,          intent(in)    :: bl_kl !< Layer index of the boundary layer (left)
  integer, optional,          intent(in)    :: bl_kr !< Layer index of the boundary layer (right)
  real, optional,             intent(in)    :: bl_zl !< Nondimensional position of the boundary layer (left)
  real, optional,             intent(in)    :: bl_zr !< Nondimensional position of the boundary layer (right)

  ! Local variables
  integer :: ns                     ! Number of neutral surfaces
  integer :: k_surface              ! Index of neutral surface
  integer :: kl                     ! Index of left interface
  integer :: kr                     ! Index of right interface
  real    :: dRdT, dRdS             ! dRho/dT [kg m-3 degC-1] and dRho/dS [kg m-3 ppt-1] for the neutral surface
  logical :: searching_left_column  ! True if searching for the position of a right interface in the left column
  logical :: searching_right_column ! True if searching for the position of a left interface in the right column
  logical :: reached_bottom         ! True if one of the bottom-most interfaces has been used as the target
  integer :: krm1, klm1
  real    :: dRho, dRhoTop, dRhoBot ! Potential density differences at various points [R ~> kg m-3]
  real    :: hL, hR                 ! Pressure thicknesses [R L2 T-2 ~> Pa]
  integer :: lastK_left, lastK_right ! Layers used during the last iteration
  real    :: lastP_left, lastP_right ! Fractional positions during the last iteration [nondim]
  logical :: interior_limit

  ns = 2*nk+2

  ! Initialize variables for the search
  kr = 1 ;
  kl = 1 ;
  lastP_right = 0.
  lastP_left = 0.
  lastK_right = 1
  lastK_left  = 1
  reached_bottom = .false.

  ! Check to see if we should limit the diffusion to the interior
  interior_limit = PRESENT(bl_kl) .and. PRESENT(bl_kr) .and. PRESENT(bl_zr) .and. PRESENT(bl_zl)

  ! Loop over each neutral surface, working from top to bottom
  neutral_surfaces: do k_surface = 1, ns
    klm1 = max(kl-1, 1)
    if (klm1>nk) stop 'find_neutral_surface_positions(): klm1 went out of bounds!'
    krm1 = max(kr-1, 1)
    if (krm1>nk) stop 'find_neutral_surface_positions(): krm1 went out of bounds!'

    ! Potential density difference, rho(kr) - rho(kl)
    dRho = 0.5 * ( ( dRdTr(kr) + dRdTl(kl) ) * ( Tr(kr) - Tl(kl) ) &
                 + ( dRdSr(kr) + dRdSl(kl) ) * ( Sr(kr) - Sl(kl) ) )
    ! Which column has the lighter surface for the current indexes, kr and kl
    if (.not. reached_bottom) then
      if (dRho < 0.) then
        searching_left_column = .true.
        searching_right_column = .false.
      elseif (dRho > 0.) then
        searching_right_column = .true.
        searching_left_column = .false.
      else ! dRho == 0.
        if (kl + kr == 2) then ! Still at surface
          searching_left_column = .true.
          searching_right_column = .false.
        else ! Not the surface so we simply change direction
          searching_left_column = .not.  searching_left_column
          searching_right_column = .not.  searching_right_column
        endif
      endif
    endif

    if (searching_left_column) then
      ! Interpolate for the neutral surface position within the left column, layer klm1
      ! Potential density difference, rho(kl-1) - rho(kr) (should be negative)
      dRhoTop = 0.5 * ( ( dRdTl(klm1) + dRdTr(kr) ) * ( Tl(klm1) - Tr(kr) ) &
                     + ( dRdSl(klm1) + dRdSr(kr) ) * ( Sl(klm1) - Sr(kr) ) )
      ! Potential density difference, rho(kl) - rho(kr) (will be positive)
      dRhoBot = 0.5 * ( ( dRdTl(klm1+1) + dRdTr(kr) ) * ( Tl(klm1+1) - Tr(kr) ) &
                      + ( dRdSl(klm1+1) + dRdSr(kr) ) * ( Sl(klm1+1) - Sr(kr) ) )

      ! Because we are looking left, the right surface, kr, is lighter than klm1+1 and should be denser than klm1
      ! unless we are still at the top of the left column (kl=1)
      if (dRhoTop > 0. .or. kr+kl==2) then
        PoL(k_surface) = 0. ! The right surface is lighter than anything in layer klm1
      elseif (dRhoTop >= dRhoBot) then ! Left layer is unstratified
        PoL(k_surface) = 1.
      else
        ! Linearly interpolate for the position between Pl(kl-1) and Pl(kl) where the density difference
        ! between right and left is zero. The Pl here are only used to handle massless layers.
        PoL(k_surface) = interpolate_for_nondim_position( dRhoTop, Pl(klm1), dRhoBot, Pl(klm1+1) )
      endif
      if (PoL(k_surface)>=1. .and. klm1<nk) then ! >= is really ==, when PoL==1 we point to the bottom of the cell
        klm1 = klm1 + 1
        PoL(k_surface) = PoL(k_surface) - 1.
      endif
      if (real(klm1-lastK_left)+(PoL(k_surface)-lastP_left)<0.) then
        PoL(k_surface) = lastP_left
        klm1 = lastK_left
      endif
      KoL(k_surface) = klm1
      if (kr <= nk) then
        PoR(k_surface) = 0.
        KoR(k_surface) = kr
      else
        PoR(k_surface) = 1.
        KoR(k_surface) = nk
      endif
      if (kr <= nk) then
        kr = kr + 1
      else
        reached_bottom = .true.
        searching_right_column = .true.
        searching_left_column = .false.
      endif
    elseif (searching_right_column) then
      ! Interpolate for the neutral surface position within the right column, layer krm1
      ! Potential density difference, rho(kr-1) - rho(kl) (should be negative)
      dRhoTop = 0.5 * ( ( dRdTr(krm1) + dRdTl(kl) ) * ( Tr(krm1) - Tl(kl) ) + &
                        ( dRdSr(krm1) + dRdSl(kl) ) * ( Sr(krm1) - Sl(kl) ) )
      ! Potential density difference, rho(kr) - rho(kl) (will be positive)
      dRhoBot = 0.5 * ( ( dRdTr(krm1+1) + dRdTl(kl) ) * ( Tr(krm1+1) - Tl(kl) ) + &
                        ( dRdSr(krm1+1) + dRdSl(kl) ) * ( Sr(krm1+1) - Sl(kl) ) )

      ! Because we are looking right, the left surface, kl, is lighter than krm1+1 and should be denser than krm1
      ! unless we are still at the top of the right column (kr=1)
      if (dRhoTop >= 0. .or. kr+kl==2) then
        PoR(k_surface) = 0. ! The left surface is lighter than anything in layer krm1
      elseif (dRhoTop >= dRhoBot) then ! Right layer is unstratified
        PoR(k_surface) = 1.
      else
        ! Linearly interpolate for the position between Pr(kr-1) and Pr(kr) where the density difference
        ! between right and left is zero. The Pr here are only used to handle massless layers.
        PoR(k_surface) = interpolate_for_nondim_position( dRhoTop, Pr(krm1), dRhoBot, Pr(krm1+1) )
      endif
      if (PoR(k_surface)>=1. .and. krm1<nk) then ! >= is really ==, when PoR==1 we point to the bottom of the cell
        krm1 = krm1 + 1
        PoR(k_surface) = PoR(k_surface) - 1.
      endif
      if (real(krm1-lastK_right)+(PoR(k_surface)-lastP_right)<0.) then
        PoR(k_surface) = lastP_right
        krm1 = lastK_right
      endif
      KoR(k_surface) = krm1
      if (kl <= nk) then
        PoL(k_surface) = 0.
        KoL(k_surface) = kl
      else
        PoL(k_surface) = 1.
        KoL(k_surface) = nk
      endif
      if (kl <= nk) then
        kl = kl + 1
      else
        reached_bottom = .true.
        searching_right_column = .false.
        searching_left_column = .true.
      endif
    else
      stop 'Else what?'
    endif
    if (interior_limit) then
      if (KoL(k_surface)<=bl_kl) then
        KoL(k_surface) = bl_kl
        if (PoL(k_surface)<bl_zl) then
          PoL(k_surface) = bl_zl
        endif
      endif
      if (KoR(k_surface)<=bl_kr) then
        KoR(k_surface) = bl_kr
        if (PoR(k_surface)<bl_zr) then
          PoR(k_surface) = bl_zr
        endif
      endif
    endif

    lastK_left = KoL(k_surface) ; lastP_left = PoL(k_surface)
    lastK_right = KoR(k_surface) ; lastP_right = PoR(k_surface)
    ! Effective thickness
    ! NOTE: This would be better expressed in terms of the layers thicknesses rather
    ! than as differences of position - AJA
    if (k_surface>1) then
      hL = absolute_position(nk,ns,Pl,KoL,PoL,k_surface) - absolute_position(nk,ns,Pl,KoL,PoL,k_surface-1)
      hR = absolute_position(nk,ns,Pr,KoR,PoR,k_surface) - absolute_position(nk,ns,Pr,KoR,PoR,k_surface-1)
      if ( hL + hR > 0.) then
        hEff(k_surface-1) = 2. * hL * hR / ( hL + hR ) ! Harmonic mean of layer thicknesses
      else
        hEff(k_surface-1) = 0.
      endif
    endif

  enddo neutral_surfaces

end subroutine find_neutral_surface_positions_continuous

!> Returns the non-dimensional position between Pneg and Ppos where the
!! interpolated density difference equals zero.
!! The result is always bounded to be between 0 and 1.
real function interpolate_for_nondim_position(dRhoNeg, Pneg, dRhoPos, Ppos)
  real, intent(in) :: dRhoNeg !< Negative density difference [R ~> kg m-3]
  real, intent(in) :: Pneg    !< Position of negative density difference [R L2 T-2 ~> Pa] or [nondim]
  real, intent(in) :: dRhoPos !< Positive density difference [R ~> kg m-3]
  real, intent(in) :: Ppos    !< Position of positive density difference [R L2 T-2 ~> Pa] or [nondim]

  character(len=120) :: mesg

  if (Ppos < Pneg) then
    call MOM_error(FATAL, 'interpolate_for_nondim_position: Houston, we have a problem! Ppos<Pneg')
  elseif (dRhoNeg>dRhoPos) then
    write(stderr,*) 'dRhoNeg, Pneg, dRhoPos, Ppos=',dRhoNeg, Pneg, dRhoPos, Ppos
    write(mesg,*) 'dRhoNeg, Pneg, dRhoPos, Ppos=', dRhoNeg, Pneg, dRhoPos, Ppos
    call MOM_error(WARNING, 'interpolate_for_nondim_position: '//trim(mesg))
  elseif (dRhoNeg>dRhoPos) then !### Does this duplicated test belong here?
    call MOM_error(FATAL, 'interpolate_for_nondim_position: Houston, we have a problem! dRhoNeg>dRhoPos')
  endif
  if (Ppos<=Pneg) then ! Handle vanished or inverted layers
    interpolate_for_nondim_position = 0.5
  elseif ( dRhoPos - dRhoNeg > 0. ) then
    interpolate_for_nondim_position = min( 1., max( 0., -dRhoNeg / ( dRhoPos - dRhoNeg ) ) )
  elseif ( dRhoPos - dRhoNeg == 0) then
    if (dRhoNeg>0.) then
      interpolate_for_nondim_position = 0.
    elseif (dRhoNeg<0.) then
      interpolate_for_nondim_position = 1.
    else ! dRhoPos = dRhoNeg = 0
      interpolate_for_nondim_position = 0.5
    endif
  else ! dRhoPos - dRhoNeg < 0
    interpolate_for_nondim_position = 0.5
  endif
  if ( interpolate_for_nondim_position < 0. ) &
    call MOM_error(FATAL, 'interpolate_for_nondim_position: Houston, we have a problem! Pint < Pneg')
  if ( interpolate_for_nondim_position > 1. ) &
    call MOM_error(FATAL, 'interpolate_for_nondim_position: Houston, we have a problem! Pint > Ppos')
end function interpolate_for_nondim_position

!> Higher order version of find_neutral_surface_positions. Returns positions within left/right columns
!! of combined interfaces using intracell reconstructions of T/S. Note that the polynomial reconstrcutions
!! of T and S are optional to aid with unit testing, but will always be passed otherwise
subroutine find_neutral_surface_positions_discontinuous(CS, nk, &
                   Pres_l, hcol_l, Tl, Sl, ppoly_T_l, ppoly_S_l, stable_l, &
                   Pres_r, hcol_r, Tr, Sr, ppoly_T_r, ppoly_S_r, stable_r, &
                   PoL, PoR, KoL, KoR, hEff, zeta_bot_L, zeta_bot_R, k_bot_L, k_bot_R, hard_fail_heff)

  type(neutral_diffusion_CS),     intent(inout) :: CS        !< Neutral diffusion control structure
  integer,                        intent(in)    :: nk        !< Number of levels
  real, dimension(nk,2),          intent(in)    :: Pres_l    !< Left-column interface pressure [R L2 T-2 ~> Pa]
  real, dimension(nk),            intent(in)    :: hcol_l    !< Left-column layer thicknesses [H ~> m or kg m-2]
                                                             !! or other units
  real, dimension(nk,2),          intent(in)    :: Tl        !< Left-column top interface potential temperature [degC]
  real, dimension(nk,2),          intent(in)    :: Sl        !< Left-column top interface salinity [ppt]
  real, dimension(:,:),           intent(in)    :: ppoly_T_l !< Left-column coefficients of T reconstruction [degC]
  real, dimension(:,:),           intent(in)    :: ppoly_S_l !< Left-column coefficients of S reconstruction [ppt]
  logical, dimension(nk),         intent(in)    :: stable_l  !< True where the left-column is stable
  real, dimension(nk,2),          intent(in)    :: Pres_r    !< Right-column interface pressure [R L2 T-2 ~> Pa]
  real, dimension(nk),            intent(in)    :: hcol_r    !< Left-column layer thicknesses [H ~> m or kg m-2]
                                                             !! or other units
  real, dimension(nk,2),          intent(in)    :: Tr        !< Right-column top interface potential temperature [degC]
  real, dimension(nk,2),          intent(in)    :: Sr        !< Right-column top interface salinity [ppt]
  real, dimension(:,:),           intent(in)    :: ppoly_T_r !< Right-column coefficients of T reconstruction [degC]
  real, dimension(:,:),           intent(in)    :: ppoly_S_r !< Right-column coefficients of S reconstruction [ppt]
  logical, dimension(nk),         intent(in)    :: stable_r  !< True where the right-column is stable
  real, dimension(4*nk),          intent(inout) :: PoL       !< Fractional position of neutral surface within
                                                             !! layer KoL of left column [nondim]
  real, dimension(4*nk),          intent(inout) :: PoR       !< Fractional position of neutral surface within
                                                             !! layer KoR of right column [nondim]
  integer, dimension(4*nk),       intent(inout) :: KoL       !< Index of first left interface above neutral surface
  integer, dimension(4*nk),       intent(inout) :: KoR       !< Index of first right interface above neutral surface
  real, dimension(4*nk-1),        intent(inout) :: hEff      !< Effective thickness between two neutral surfaces
                                                             !! [H ~> m or kg m-2] or other units taken from hcol_l
  real, optional,                 intent(in)    :: zeta_bot_L!< Non-dimensional distance to where the boundary layer
                                                             !! intersetcs the cell (left) [nondim]
  real, optional,                 intent(in)    :: zeta_bot_R!< Non-dimensional distance to where the boundary layer
                                                             !! intersetcs the cell (right) [nondim]

  integer, optional,              intent(in)    :: k_bot_L   !< k-index for the boundary layer (left) [nondim]
  integer, optional,              intent(in)    :: k_bot_R   !< k-index for the boundary layer (right) [nondim]
  logical, optional,              intent(in)    :: hard_fail_heff !< If true (default) bring down the model if the
                                                             !! neutral surfaces ever cross [logical]
  ! Local variables
  integer :: ns                     ! Number of neutral surfaces
  integer :: k_surface              ! Index of neutral surface
  integer :: kl_left, kl_right      ! Index of layers on the left/right
  integer :: ki_left, ki_right      ! Index of interfaces on the left/right
  logical :: searching_left_column  ! True if searching for the position of a right interface in the left column
  logical :: searching_right_column ! True if searching for the position of a left interface in the right column
  logical :: reached_bottom         ! True if one of the bottom-most interfaces has been used as the target
  logical :: search_layer
  logical :: fail_heff              ! Fail if negative thickness are encountered.  By default this
                                    ! is true, but it can take its value from hard_fail_heff.
  real    :: dRho                   ! A density difference between columns [R ~> kg m-3]
  real    :: hL, hR                 ! Left and right layer thicknesses [H ~> m or kg m-2] or units from hcol_l
  real    :: lastP_left, lastP_right ! Previous positions for left and right [nondim]
  integer :: k_init_L, k_init_R      ! Starting indices layers for left and right
  real    :: p_init_L, p_init_R      ! Starting positions for left and right [nondim]
  ! Initialize variables for the search
  ns = 4*nk
  ki_right = 1
  ki_left = 1
  kl_left = 1
  kl_right = 1
  lastP_left = 0.
  lastP_right = 0.
  reached_bottom = .false.
  searching_left_column = .false.
  searching_right_column = .false.

  fail_heff = .true.
  if (PRESENT(hard_fail_heff)) fail_heff = hard_fail_heff

  if (PRESENT(k_bot_L) .and. PRESENT(k_bot_R) .and. PRESENT(zeta_bot_L) .and. PRESENT(zeta_bot_R)) then
    k_init_L = k_bot_L; k_init_R = k_bot_R
    p_init_L = zeta_bot_L; p_init_R = zeta_bot_R
    lastP_left = zeta_bot_L; lastP_right = zeta_bot_R
    kl_left = k_bot_L; kl_right = k_bot_R
  else
    k_init_L = 1  ; k_init_R = 1
    p_init_L = 0. ; p_init_R = 0.
  endif
  ! Loop over each neutral surface, working from top to bottom
  neutral_surfaces: do k_surface = 1, ns

    if (k_surface == ns) then
        PoL(k_surface) = 1.
        PoR(k_surface) = 1.
        KoL(k_surface) = nk
        KoR(k_surface) = nk
    ! If the layers are unstable, then simply point the surface to the previous location
    elseif (.not. stable_l(kl_left)) then
      if (k_surface > 1) then
        PoL(k_surface) = ki_left - 1 ! Top interface is at position = 0., Bottom is at position = 1
        KoL(k_surface) = kl_left
        PoR(k_surface) = PoR(k_surface-1)
        KoR(k_surface) = KoR(k_surface-1)
      else
        PoR(k_surface) = p_init_R
        KoR(k_surface) = k_init_R
        PoL(k_surface) = p_init_L
        KoL(k_Surface) = k_init_L
      endif
      call increment_interface(nk, kl_left, ki_left, reached_bottom, searching_left_column, searching_right_column)
      searching_left_column = .true.
      searching_right_column = .false.
    elseif (.not. stable_r(kl_right)) then ! Check the right layer for stability
      if (k_surface > 1) then
        PoR(k_surface) = ki_right - 1 ! Top interface is at position = 0., Bottom is at position = 1
        KoR(k_surface) = kl_right
        PoL(k_surface) = PoL(k_surface-1)
        KoL(k_surface) = KoL(k_surface-1)
      else
        PoR(k_surface) = 0.
        KoR(k_surface) = 1
        PoL(k_surface) = 0.
        KoL(k_surface) = 1
      endif
      call increment_interface(nk, kl_right, ki_right, reached_bottom, searching_right_column, searching_left_column)
      searching_left_column = .false.
      searching_right_column = .true.
    else ! Layers are stable so need to figure out whether we need to search right or left
      ! For convenience, the left column uses the searched "from" interface variables, and the right column
      ! uses the searched 'to'. These will get reset in subsequent calc_delta_rho calls

      call calc_delta_rho_and_derivs(CS,                                                                        &
                                     Tr(kl_right, ki_right), Sr(kl_right, ki_right), Pres_r(kl_right,ki_right), &
                                     Tl(kl_left, ki_left),   Sl(kl_left, ki_left)  , Pres_l(kl_left,ki_left),   &
                                     dRho)
      if (CS%debug) write(stdout,'(A,I2,A,E12.4,A,I2,A,I2,A,I2,A,I2)') &
          "k_surface=",k_surface, "  dRho=",CS%R_to_kg_m3*dRho, &
          "kl_left=",kl_left, "  ki_left=",ki_left, "  kl_right=",kl_right, "  ki_right=",ki_right
      ! Which column has the lighter surface for the current indexes, kr and kl
      if (.not. reached_bottom) then
        if (dRho < 0.) then
          searching_left_column  = .true.
          searching_right_column = .false.
        elseif (dRho > 0.) then
          searching_left_column  = .false.
          searching_right_column = .true.
        else ! dRho == 0.
          if (  ( kl_left + kl_right == 2 ) .and. (ki_left + ki_right == 2) ) then ! Still at surface
            searching_left_column  = .true.
            searching_right_column = .false.
          else ! Not the surface so we simply change direction
            searching_left_column = .not. searching_left_column
            searching_right_column = .not. searching_right_column
          endif
        endif
      endif
      if (searching_left_column) then
        ! Position of the right interface is known and all quantities are fixed
        PoR(k_surface) = ki_right - 1.
        KoR(k_surface) = kl_right
        PoL(k_surface) = search_other_column(CS, k_surface, lastP_left,                                &
                           Tr(kl_right, ki_right), Sr(kl_right, ki_right), Pres_r(kl_right, ki_right), &
                           Tl(kl_left,1),          Sl(kl_left,1),          Pres_l(kl_left,1),          &
                           Tl(kl_left,2),          Sl(kl_left,2),          Pres_l(kl_left,2),          &
                           ppoly_T_l(kl_left,:), ppoly_S_l(kl_left,:))
        KoL(k_surface) = kl_left

        if (CS%debug) then
          write(stdout,'(A,I2)') "Searching left layer ", kl_left
          write(stdout,'(A,I2,X,I2)') "Searching from right: ", kl_right, ki_right
          write(stdout,*) "Temp/Salt Reference: ", Tr(kl_right,ki_right), Sr(kl_right,ki_right)
          write(stdout,*) "Temp/Salt Top L: ", Tl(kl_left,1), Sl(kl_left,1)
          write(stdout,*) "Temp/Salt Bot L: ", Tl(kl_left,2), Sl(kl_left,2)
        endif
        call increment_interface(nk, kl_right, ki_right, reached_bottom, searching_right_column, searching_left_column)
        lastP_left = PoL(k_surface)
        ! If the right layer increments, then we need to reset the last position on the right
        if ( kl_right == (KoR(k_surface) + 1) ) lastP_right = 0.

      elseif (searching_right_column) then
        ! Position of the right interface is known and all quantities are fixed
        PoL(k_surface) = ki_left - 1.
        KoL(k_surface) = kl_left
        PoR(k_surface) = search_other_column(CS, k_surface, lastP_right,                         &
                           Tl(kl_left, ki_left), Sl(kl_left, ki_left), Pres_l(kl_left, ki_left), &
                           Tr(kl_right,1),       Sr(kl_right,1),       Pres_r(kl_right,1),       &
                           Tr(kl_right,2),       Sr(kl_right,2),       Pres_r(kl_right,2),       &
                           ppoly_T_r(kl_right,:), ppoly_S_r(kl_right,:))
        KoR(k_surface) = kl_right

        if (CS%debug) then
          write(stdout,'(A,I2)') "Searching right layer ", kl_right
          write(stdout,'(A,I2,X,I2)') "Searching from left: ", kl_left, ki_left
          write(stdout,*) "Temp/Salt Reference: ", Tl(kl_left,ki_left), Sl(kl_left,ki_left)
          write(stdout,*) "Temp/Salt Top L: ", Tr(kl_right,1), Sr(kl_right,1)
          write(stdout,*) "Temp/Salt Bot L: ", Tr(kl_right,2), Sr(kl_right,2)
        endif
        call increment_interface(nk, kl_left, ki_left, reached_bottom, searching_left_column, searching_right_column)
        lastP_right = PoR(k_surface)
        ! If the right layer increments, then we need to reset the last position on the right
        if ( kl_left == (KoL(k_surface) + 1) ) lastP_left = 0.
      else
        stop 'Else what?'
      endif
      if (CS%debug)  write(stdout,'(A,I3,A,ES16.6,A,I2,A,ES16.6)') "KoL:", KoL(k_surface), " PoL:", PoL(k_surface), &
                     "     KoR:", KoR(k_surface), " PoR:", PoR(k_surface)
    endif
    ! Effective thickness
    if (k_surface>1) then
      if ( KoL(k_surface) == KoL(k_surface-1) .and. KoR(k_surface) == KoR(k_surface-1) ) then
        hL = (PoL(k_surface) - PoL(k_surface-1))*hcol_l(KoL(k_surface))
        hR = (PoR(k_surface) - PoR(k_surface-1))*hcol_r(KoR(k_surface))
        if (hL < 0. .or. hR < 0.) then
          if (fail_heff) then
            call MOM_error(FATAL,"Negative thicknesses in neutral diffusion")
          else
            if (searching_left_column) then
              PoL(k_surface) = PoL(k_surface-1)
              KoL(k_surface) = KoL(k_surface-1)
            elseif (searching_right_column) then
              PoR(k_surface) = PoR(k_surface-1)
              KoR(k_surface) = KoR(k_surface-1)
            endif
          endif
        elseif ( hL + hR == 0. ) then
          hEff(k_surface-1) = 0.
        else
          hEff(k_surface-1) = 2. * ( (hL * hR) / ( hL + hR ) )! Harmonic mean
          if ( KoL(k_surface) /= KoL(k_surface-1) ) then
            call MOM_error(FATAL,"Neutral sublayer spans multiple layers")
          endif
          if ( KoR(k_surface) /= KoR(k_surface-1) ) then
            call MOM_error(FATAL,"Neutral sublayer spans multiple layers")
          endif
        endif
      else
        hEff(k_surface-1) = 0.
      endif
    endif
  enddo neutral_surfaces
end subroutine find_neutral_surface_positions_discontinuous

!> Sweep down through the column and mark as stable if the bottom interface of a cell is denser than the top
subroutine mark_unstable_cells(CS, nk, T, S, P, stable_cell)
  type(neutral_diffusion_CS), intent(inout) :: CS      !< Neutral diffusion control structure
  integer,                intent(in)    :: nk          !< Number of levels in a column
  real, dimension(nk,2),  intent(in)    :: T           !< Temperature at interfaces [degC]
  real, dimension(nk,2),  intent(in)    :: S           !< Salinity at interfaces [ppt]
  real, dimension(nk,2),  intent(in)    :: P           !< Pressure at interfaces [R L2 T-2 ~> Pa]
  logical, dimension(nk), intent(  out) :: stable_cell !< True if this cell is unstably stratified

  integer :: k, first_stable, prev_stable
  real :: delta_rho ! A density difference [R ~> kg m-3]

  do k = 1,nk
    call calc_delta_rho_and_derivs( CS, T(k,2), S(k,2), max(P(k,2), CS%ref_pres), &
                                        T(k,1), S(k,1), max(P(k,1), CS%ref_pres), delta_rho )
    stable_cell(k) = (delta_rho > 0.)
  enddo
end subroutine mark_unstable_cells

!> Searches the "other" (searched) column for the position of the neutral surface
real function search_other_column(CS, ksurf, pos_last, T_from, S_from, P_from, T_top, S_top, P_top, &
                                  T_bot, S_bot, P_bot, T_poly, S_poly ) result(pos)
  type(neutral_diffusion_CS), intent(in   ) :: CS       !< Neutral diffusion control structure
  integer,                    intent(in   ) :: ksurf    !< Current index of neutral surface
  real,                       intent(in   ) :: pos_last !< Last position within the current layer, used as the lower
                                                        !! bound in the root finding algorithm [nondim]
  real,                       intent(in   ) :: T_from   !< Temperature at the searched from interface [degC]
  real,                       intent(in   ) :: S_from   !< Salinity    at the searched from interface [ppt]
  real,                       intent(in   ) :: P_from   !< Pressure at the searched from interface [R L2 T-2 ~> Pa]
  real,                       intent(in   ) :: T_top    !< Temperature at the searched to top interface [degC]
  real,                       intent(in   ) :: S_top    !< Salinity    at the searched to top interface [ppt]
  real,                       intent(in   ) :: P_top    !< Pressure at the searched to top interface [R L2 T-2 ~> Pa]
                                                        !! interface [R L2 T-2 ~> Pa]
  real,                       intent(in   ) :: T_bot    !< Temperature at the searched to bottom interface [degC]
  real,                       intent(in   ) :: S_bot    !< Salinity    at the searched to bottom interface [ppt]
  real,                       intent(in   ) :: P_bot    !< Pressure at the searched to bottom
                                                        !! interface [R L2 T-2 ~> Pa]
  real, dimension(:),         intent(in   ) :: T_poly   !< Temperature polynomial reconstruction coefficients [degC]
  real, dimension(:),         intent(in   ) :: S_poly   !< Salinity    polynomial reconstruction coefficients [ppt]
  ! Local variables
  real :: dRhotop, dRhobot ! Density differences [R ~> kg m-3]
  real :: dRdT_top, dRdT_bot, dRdT_from ! Partial derivatives of density with temperature [R degC-1 ~> kg m-3 degC-1]
  real :: dRdS_top, dRdS_bot, dRdS_from ! Partial derivatives of density with salinity [R ppt-1 ~> kg m-3 ppt-1]

  ! Calculate the differencei in density at the tops or the bottom
  if (CS%neutral_pos_method == 1 .or. CS%neutral_pos_method == 3) then
    call calc_delta_rho_and_derivs(CS, T_top, S_top, P_top, T_from, S_from, P_from, dRhoTop)
    call calc_delta_rho_and_derivs(CS, T_bot, S_bot, P_bot, T_from, S_from, P_from, dRhoBot)
  elseif (CS%neutral_pos_method == 2) then
    call calc_delta_rho_and_derivs(CS, T_top, S_top, P_top, T_from, S_from, P_from, dRhoTop, &
                                   dRdT_top, dRdS_top, dRdT_from, dRdS_from)
    call calc_delta_rho_and_derivs(CS, T_bot, S_bot, P_bot, T_from, S_from, P_from, dRhoBot, &
                                   dRdT_bot, dRdS_bot, dRdT_from, dRdS_from)
  endif

  ! Handle all the special cases EXCEPT if it connects within the layer
  if ( (dRhoTop > 0.) .or. (ksurf == 1) ) then      ! First interface or lighter than anything in layer
    pos = pos_last
  elseif ( dRhoTop > dRhoBot ) then                 ! Unstably stratified
    pos = 1.
  elseif ( dRhoTop < 0. .and. dRhoBot < 0.) then    ! Denser than anything in layer
    pos = 1.
  elseif ( dRhoTop == 0. .and. dRhoBot == 0. ) then ! Perfectly unstratified
    pos = 1.
  elseif ( dRhoBot == 0. ) then                     ! Matches perfectly at the Top
    pos = 1.
  elseif ( dRhoTop == 0. ) then                     ! Matches perfectly at the Bottom
    pos = pos_last
  else                                              ! Neutral surface within layer
    pos = -1
  endif

  ! Can safely return if position is >= 0 otherwise will need to find the position within the layer
  if (pos>=0) return

  if (CS%neutral_pos_method==1) then
    pos = interpolate_for_nondim_position( dRhoTop, P_top, dRhoBot, P_bot )
  ! For the 'Linear' case of finding the neutral position, the fromerence pressure to use is the average
  ! of the midpoint of the layer being searched and the interface being searched from
  elseif (CS%neutral_pos_method == 2) then
    pos = find_neutral_pos_linear( CS, pos_last, T_from, S_from, dRdT_from, dRdS_from, &
                                   dRdT_top, dRdS_top, dRdT_bot, dRdS_bot, T_poly, S_poly )
  elseif (CS%neutral_pos_method == 3) then
    pos = find_neutral_pos_full( CS, pos_last, T_from, S_from, P_from, P_top, P_bot, T_poly, S_poly)
  endif

end function search_other_column

!> Increments the interface which was just connected and also set flags if the bottom is reached
subroutine increment_interface(nk, kl, ki, reached_bottom, searching_this_column, searching_other_column)
  integer, intent(in   )                :: nk                     !< Number of vertical levels
  integer, intent(inout)                :: kl                     !< Current layer (potentially updated)
  integer, intent(inout)                :: ki                     !< Current interface
  logical, intent(inout)                :: reached_bottom         !< Updated when kl == nk and ki == 2
  logical, intent(inout)                :: searching_this_column  !< Updated when kl == nk and ki == 2
  logical, intent(inout)                :: searching_other_column !< Updated when kl == nk and ki == 2
  integer :: k

  reached_bottom = .false.
  if (ki == 2) then ! At the bottom interface
    if ((ki == 2) .and. (kl < nk) ) then ! Not at the bottom so just go to the next layer
      kl = kl+1
      ki = 1
    elseif ((kl == nk) .and. (ki==2)) then
      reached_bottom = .true.
      searching_this_column = .false.
      searching_other_column = .true.
    endif
  elseif (ki==1) then ! At the top interface
    ki = 2 ! Next interface is same layer, but bottom interface
  else
    call MOM_error(FATAL,"Unanticipated eventuality in increment_interface")
  endif
end subroutine increment_interface

!> Search a layer to find where delta_rho = 0 based on a linear interpolation of alpha and beta of the top and bottom
!! being searched and polynomial reconstructions of T and S. Compressibility is not needed because either, we are
!! assuming incompressibility in the equation of state for this module or alpha and beta are calculated having been
!! displaced to the average pressures of the two pressures We need Newton's method because the T and S reconstructions
!! make delta_rho a polynomial function of z if using PPM or higher. If Newton's method would search fall out of the
!! interval [0,1], a bisection step would be taken instead. Also this linearization of alpha, beta means that second
!! derivatives of the EOS are not needed. Note that delta in variable names below refers to horizontal differences and
!! 'd' refers to vertical differences
function find_neutral_pos_linear( CS, z0, T_ref, S_ref, dRdT_ref, dRdS_ref, &
                                  dRdT_top, dRdS_top, dRdT_bot, dRdS_bot, ppoly_T, ppoly_S ) result( z )
  type(neutral_diffusion_CS),intent(in) :: CS        !< Control structure with parameters for this module
  real,                      intent(in) :: z0        !< Lower bound of position, also serves as the
                                                     !! initial guess [nondim]
  real,                      intent(in) :: T_ref     !< Temperature at the searched from interface [degC]
  real,                      intent(in) :: S_ref     !< Salinity at the searched from interface [ppt]
  real,                      intent(in) :: dRdT_ref  !< dRho/dT at the searched from interface
                                                     !! [R degC-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_ref  !< dRho/dS at the searched from interface
                                                     !! [R ppt-1 ~> kg m-3 ppt-1]
  real,                      intent(in) :: dRdT_top  !< dRho/dT at top of layer being searched
                                                     !! [R degC-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_top  !< dRho/dS at top of layer being searched
                                                     !! [R ppt-1 ~> kg m-3 ppt-1]
  real,                      intent(in) :: dRdT_bot  !< dRho/dT at bottom of layer being searched
                                                     !! [R degC-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_bot  !< dRho/dS at bottom of layer being searched
                                                     !! [R ppt-1 ~> kg m-3 ppt-1]
  real, dimension(:),        intent(in) :: ppoly_T   !< Coefficients of the polynomial reconstruction of T within
                                                     !! the layer to be searched [degC].
  real, dimension(:),        intent(in) :: ppoly_S   !< Coefficients of the polynomial reconstruction of S within
                                                     !! the layer to be searched [ppt].
  real                                  :: z         !< Position where drho = 0 [nondim]
  ! Local variables
  real :: dRdT_diff  ! Difference in the partial derivative of density with temperature across the
                     ! layer [R degC-1 ~> kg m-3 degC-1]
  real :: dRdS_diff  ! Difference in the partial derivative of density with salinity across the
                     ! layer [R ppt-1 ~> kg m-3 ppt-1]
  real :: drho, drho_dz ! Density anomaly and its derivative with fracitonal position [R ~> kg m-3]
  real :: dRdT_z     ! Partial derivative of density with temperature at a point [R degC-1 ~> kg m-3 degC-1]
  real :: dRdS_z     ! Partial derivative of density with salinity at a point [R ppt-1 ~> kg m-3 ppt-1]
  real :: T_z, dT_dz ! Temperature at a point and its derivative with fractional position [degC]
  real :: S_z, dS_dz ! Salinity at a point and its derivative with fractional position [ppt]
  real :: drho_min, drho_max ! Bounds on density differences [R ~> kg m-3]
  real :: ztest, zmin, zmax ! Fractional positions in the cell [nondim]
  real :: dz         ! Change in position in the cell [nondim]
  real :: a1, a2     ! Fractional weights of the top and bottom values [nondim]
  integer :: iter
  integer :: nterm

  nterm = SIZE(ppoly_T)

  ! Position independent quantities
  dRdT_diff = dRdT_bot - dRdT_top
  dRdS_diff = dRdS_bot - dRdS_top
  ! Initial starting drho (used for bisection)
  zmin = z0        ! Lower bounding interval
  zmax = 1.        ! Maximum bounding interval (bottom of layer)
  a1 = 1. - zmin
  a2 = zmin
  T_z = evaluation_polynomial( ppoly_T, nterm, zmin )
  S_z = evaluation_polynomial( ppoly_S, nterm, zmin )
  dRdT_z = a1*dRdT_top + a2*dRdT_bot
  dRdS_z = a1*dRdS_top + a2*dRdS_bot
  drho_min = 0.5*((dRdT_z+dRdT_ref)*(T_z-T_ref) + (dRdS_z+dRdS_ref)*(S_z-S_ref))

  T_z = evaluation_polynomial( ppoly_T, nterm, 1. )
  S_z = evaluation_polynomial( ppoly_S, nterm, 1. )
  drho_max = 0.5*((dRdT_bot+dRdT_ref)*(T_z-T_ref) + (dRdS_bot+dRdS_ref)*(S_z-S_ref))

  if (drho_min >= 0.) then
    z = z0
    return
  elseif (drho_max == 0.) then
    z = 1.
    return
  endif
  if ( SIGN(1.,drho_min) == SIGN(1.,drho_max) ) then
    call MOM_error(FATAL, "drho_min is the same sign as dhro_max")
  endif

  z = z0
  ztest = z0
  do iter = 1, CS%max_iter
    ! Calculate quantities at the current nondimensional position
    a1 = 1.-z
    a2 = z
    dRdT_z    = a1*dRdT_top + a2*dRdT_bot
    dRdS_z    = a1*dRdS_top + a2*dRdS_bot
    T_z       = evaluation_polynomial( ppoly_T, nterm, z )
    S_z       = evaluation_polynomial( ppoly_S, nterm, z )
    drho = 0.5*((dRdT_z+dRdT_ref)*(T_z-T_ref) + (dRdS_z+dRdS_ref)*(S_z-S_ref))

    ! Check for convergence
    if (ABS(drho) <= CS%drho_tol) exit
    ! Update bisection bracketing intervals
    if (drho < 0. .and. drho > drho_min) then
      drho_min = drho
      zmin = z
    elseif (drho > 0. .and. drho < drho_max) then
      drho_max = drho
      zmax = z
    endif

    ! Calculate a Newton step
    dT_dz = first_derivative_polynomial( ppoly_T, nterm, z )
    dS_dz = first_derivative_polynomial( ppoly_S, nterm, z )
    drho_dz = 0.5*( (dRdT_diff*(T_z - T_ref) + (dRdT_ref+dRdT_z)*dT_dz) + &
                    (dRdS_diff*(S_z - S_ref) + (dRdS_ref+dRdS_z)*dS_dz) )

    ztest = z - drho/drho_dz
    ! Take a bisection if z falls out of [zmin,zmax]
    if (ztest < zmin .or. ztest > zmax) then
      if ( drho < 0. ) then
        ztest = 0.5*(z + zmax)
      else
        ztest = 0.5*(zmin + z)
      endif
    endif

    ! Test to ensure we haven't stalled out
    if ( abs(z-ztest) <= CS%x_tol ) exit
    ! Reset for next iteration
    z = ztest
  enddo

end function find_neutral_pos_linear

!> Use the full equation of state to calculate the difference in locally referenced potential density. The derivatives
!! in this case are not trivial to calculate, so instead we use a regula falsi method
function find_neutral_pos_full( CS, z0, T_ref, S_ref, P_ref, P_top, P_bot, ppoly_T, ppoly_S ) result( z )
  type(neutral_diffusion_CS),intent(in) :: CS        !< Control structure with parameters for this module
  real,                      intent(in) :: z0        !< Lower bound of position, also serves as the
                                                     !! initial guess [nondim]
  real,                      intent(in) :: T_ref     !< Temperature at the searched from interface [degC]
  real,                      intent(in) :: S_ref     !< Salinity at the searched from interface [ppt]
  real,                      intent(in) :: P_ref     !< Pressure at the searched from interface [R L2 T-2 ~> Pa]
  real,                      intent(in) :: P_top     !< Pressure at top of layer being searched [R L2 T-2 ~> Pa]
  real,                      intent(in) :: P_bot     !< Pressure at bottom of layer being searched [R L2 T-2 ~> Pa]
  real, dimension(:),        intent(in) :: ppoly_T   !< Coefficients of the polynomial reconstruction of T within
                                                     !! the layer to be searched [degC]
  real, dimension(:),        intent(in) :: ppoly_S   !< Coefficients of the polynomial reconstruction of T within
                                                     !! the layer to be searched [ppt]
  real                                  :: z         !< Position where drho = 0 [nondim]
  ! Local variables
  integer :: iter
  integer :: nterm

  real :: drho_a, drho_b, drho_c ! Density differences [R ~> kg m-3]
  real :: a, b, c     ! Fractional positions [nondim]
  real :: Ta, Tb, Tc  ! Temperatures [degC]
  real :: Sa, Sb, Sc  ! Salinities [ppt]
  real :: Pa, Pb, Pc  ! Pressures [R L2 T-2 ~> Pa]
  integer :: side

  side = 0
  ! Set the first two evaluation to the endpoints of the interval
  b = z0 ; c = 1
  nterm = SIZE(ppoly_T)

  ! Calculate drho at the minimum bound
  Tb = evaluation_polynomial( ppoly_T, nterm, b )
  Sb = evaluation_polynomial( ppoly_S, nterm, b )
  Pb = P_top*(1.-b) + P_bot*b
  call calc_delta_rho_and_derivs(CS, Tb, Sb, Pb, T_ref, S_ref, P_ref, drho_b)

  ! Calculate drho at the maximum bound
  Tc = evaluation_polynomial( ppoly_T, nterm, 1. )
  Sc = evaluation_polynomial( ppoly_S, nterm, 1. )
  Pc = P_Bot
  call calc_delta_rho_and_derivs(CS, Tc, Sc, Pc, T_ref, S_ref, P_ref, drho_c)

  if (drho_b >= 0.) then
    z = z0
    return
  elseif (drho_c == 0.) then
    z = 1.
    return
  endif
  if ( SIGN(1.,drho_b) == SIGN(1.,drho_c) ) then
    z = z0
    return
  endif

  do iter = 1, CS%max_iter
    ! Calculate new position and evaluate if we have converged
    a = (drho_b*c - drho_c*b)/(drho_b-drho_c)
    Ta = evaluation_polynomial( ppoly_T, nterm, a )
    Sa = evaluation_polynomial( ppoly_S, nterm, a )
    Pa = P_top*(1.-a) + P_bot*a
    call calc_delta_rho_and_derivs(CS, Ta, Sa, Pa, T_ref, S_ref, P_ref, drho_a)
    if (ABS(drho_a) < CS%drho_tol) then
      z = a
      return
    endif

    if (drho_a*drho_c > 0.) then
      if ( ABS(a-c)<CS%x_tol) then
        z = a
        return
      endif
      c = a ; drho_c = drho_a;
      if (side == -1) drho_b = 0.5*drho_b
      side = -1
    elseif ( drho_b*drho_a > 0 ) then
      if ( ABS(a-b)<CS%x_tol) then
        z = a
        return
      endif
      b = a ; drho_b = drho_a
      if (side == 1) drho_c = 0.5*drho_c
      side = 1
    else
      z = a
      return
    endif
  enddo

  z = a

end function find_neutral_pos_full

!> Calculate the difference in density between two points in a variety of ways
subroutine calc_delta_rho_and_derivs(CS, T1, S1, p1_in, T2, S2, p2_in, drho, &
                                     drdt1_out, drds1_out, drdt2_out, drds2_out )
  type(neutral_diffusion_CS)    :: CS        !< Neutral diffusion control structure
  real,           intent(in   ) :: T1        !< Temperature at point 1 [degC]
  real,           intent(in   ) :: S1        !< Salinity at point 1 [ppt]
  real,           intent(in   ) :: p1_in     !< Pressure at point 1 [R L2 T-2 ~> Pa]
  real,           intent(in   ) :: T2        !< Temperature at point 2 [degC]
  real,           intent(in   ) :: S2        !< Salinity at point 2 [ppt]
  real,           intent(in   ) :: p2_in     !< Pressure at point 2 [R L2 T-2 ~> Pa]
  real,           intent(  out) :: drho      !< Difference in density between the two points [R ~> kg m-3]
  real, optional, intent(  out) :: dRdT1_out !< drho_dt at point 1 [R degC-1 ~> kg m-3 degC-1]
  real, optional, intent(  out) :: dRdS1_out !< drho_ds at point 1 [R ppt-1 ~> kg m-3 ppt-1]
  real, optional, intent(  out) :: dRdT2_out !< drho_dt at point 2 [R degC-1 ~> kg m-3 degC-1]
  real, optional, intent(  out) :: dRdS2_out !< drho_ds at point 2 [R ppt-1 ~> kg m-3 ppt-1]
  ! Local variables
  real :: rho1, rho2   ! Densities [R ~> kg m-3]
  real :: p1, p2, pmid ! Pressures [R L2 T-2 ~> Pa]
  real :: drdt1, drdt2 ! Partial derivatives of density with temperature [R degC-1 ~> kg m-3 degC-1]
  real :: drds1, drds2 ! Partial derivatives of density with salinity [R ppt-1 ~> kg m-3 ppt-1]
  real :: drdp1, drdp2 ! Partial derivatives of density with pressure [T2 L-2 ~> s2 m-2]

  ! Use the same reference pressure or the in-situ pressure
  if (CS%ref_pres > 0.) then
    p1 = CS%ref_pres
    p2 = CS%ref_pres
  else
    p1 = p1_in
    p2 = p2_in
  endif

  ! Use the full linear equation of state to calculate the difference in density (expensive!)
  if     (TRIM(CS%delta_rho_form) == 'full') then
    pmid = 0.5 * (p1 + p2)
    call calculate_density( T1, S1, pmid, rho1, CS%EOS)
    call calculate_density( T2, S2, pmid, rho2, CS%EOS)
    drho = rho1 - rho2
  ! Use the density derivatives at the average of pressures and the differentces int temperature
  elseif (TRIM(CS%delta_rho_form) == 'mid_pressure') then
    pmid = 0.5 * (p1 + p2)
    if (CS%ref_pres>=0) pmid = CS%ref_pres
    call calculate_density_derivs(T1, S1, pmid, drdt1, drds1, CS%EOS)
    call calculate_density_derivs(T2, S2, pmid, drdt2, drds2, CS%EOS)
    drho = delta_rho_from_derivs( T1, S1, p1, drdt1, drds1, T2, S2, p2, drdt2, drds2)
  elseif (TRIM(CS%delta_rho_form) == 'local_pressure') then
    call calculate_density_derivs(T1, S1, p1, drdt1, drds1, CS%EOS)
    call calculate_density_derivs(T2, S2, p2, drdt2, drds2, CS%EOS)
    drho = delta_rho_from_derivs( T1, S1, p1, drdt1, drds1, T2, S2, p2, drdt2, drds2)
  else
    call MOM_error(FATAL, "delta_rho_form is not recognized")
  endif

  if (PRESENT(drdt1_out)) drdt1_out = drdt1
  if (PRESENT(drds1_out)) drds1_out = drds1
  if (PRESENT(drdt2_out)) drdt2_out = drdt2
  if (PRESENT(drds2_out)) drds2_out = drds2

end subroutine calc_delta_rho_and_derivs

!> Calculate delta rho from derivatives and gradients of properties
!! \f$ \Delta \rho = \frac{1}{2}\left[ (\alpha_1 + \alpha_2)*(T_1-T_2) +
!!                                   (\beta_1 + \beta_2)*(S_1-S_2) +
!!                                   (\gamma^{-1}_1 + \gamma^{-1}_2)*(P_1-P_2) \right] \f$
function delta_rho_from_derivs( T1, S1, P1, dRdT1, dRdS1, &
                                T2, S2, P2, dRdT2, dRdS2  ) result (drho)
  real :: T1    !< Temperature at point 1 [degC]
  real :: S1    !< Salinity at point 1 [ppt]
  real :: P1    !< Pressure at point 1 [R L2 T-2 ~> Pa]
  real :: dRdT1 !< The partial derivative of density with temperature at point 1 [R degC-1 ~> kg m-3 degC-1]
  real :: dRdS1 !< The partial derivative of density with salinity at point 1 [R ppt-1 ~> kg m-3 ppt-1]
  real :: T2    !< Temperature at point 2 [degC]
  real :: S2    !< Salinity at point 2 [ppt]
  real :: P2    !< Pressure at point 2 [R L2 T-2 ~> Pa]
  real :: dRdT2 !< The partial derivative of density with temperature at point 2 [R degC-1 ~> kg m-3 degC-1]
  real :: dRdS2 !< The partial derivative of density with salinity at point 2 [R ppt-1 ~> kg m-3 ppt-1]
  ! Local variables
  real :: drho  ! The density difference [R ~> kg m-3]

  drho = 0.5 * ( (dRdT1+dRdT2)*(T1-T2) + (dRdS1+dRdS2)*(S1-S2))

end function delta_rho_from_derivs

!> Converts non-dimensional position within a layer to absolute position (for debugging)
function absolute_position(n,ns,Pint,Karr,NParr,k_surface)
  integer, intent(in) :: n            !< Number of levels
  integer, intent(in) :: ns           !< Number of neutral surfaces
  real,    intent(in) :: Pint(n+1)    !< Position of interfaces [R L2 T-2 ~> Pa] or other units
  integer, intent(in) :: Karr(ns)     !< Index of interface above position
  real,    intent(in) :: NParr(ns)    !< Non-dimensional position within layer Karr(:) [nondim]
  integer, intent(in) :: k_surface    !< k-interface to query
  real                :: absolute_position !< The absolute position of a location [R L2 T-2 ~> Pa]
                                      !! or other units following Pint
  ! Local variables
  integer :: k

  k = Karr(k_surface)
  if (k>n) stop 'absolute_position: k>nk is out of bounds!'
  absolute_position = Pint(k) + NParr(k_surface) * ( Pint(k+1) - Pint(k) )

end function absolute_position

!> Converts non-dimensional positions within layers to absolute positions (for debugging)
function absolute_positions(n,ns,Pint,Karr,NParr)
  integer, intent(in) :: n         !< Number of levels
  integer, intent(in) :: ns        !< Number of neutral surfaces
  real,    intent(in) :: Pint(n+1) !< Position of interface [R L2 T-2 ~> Pa] or other units
  integer, intent(in) :: Karr(ns)  !< Indexes of interfaces about positions
  real,    intent(in) :: NParr(ns) !< Non-dimensional positions within layers Karr(:)

  real,  dimension(ns) :: absolute_positions !< Absolute positions [R L2 T-2 ~> Pa]
                                   !! or other units following Pint

  ! Local variables
  integer :: k_surface, k

  do k_surface = 1, ns
    absolute_positions(k_surface) = absolute_position(n,ns,Pint,Karr,NParr,k_surface)
  enddo

end function absolute_positions

!> Returns a single column of neutral diffusion fluxes of a tracer.
subroutine neutral_surface_flux(nk, nsurf, deg, hl, hr, Tl, Tr, PiL, PiR, KoL, KoR, &
                                hEff, Flx, continuous, h_neglect, remap_CS, h_neglect_edge)
  integer,                      intent(in)    :: nk    !< Number of levels
  integer,                      intent(in)    :: nsurf !< Number of neutral surfaces
  integer,                      intent(in)    :: deg   !< Degree of polynomial reconstructions
  real, dimension(nk),          intent(in)    :: hl    !< Left-column layer thickness [H ~> m or kg m-2]
  real, dimension(nk),          intent(in)    :: hr    !< Right-column layer thickness [H ~> m or kg m-2]
  real, dimension(nk),          intent(in)    :: Tl    !< Left-column layer tracer (conc, e.g. degC)
  real, dimension(nk),          intent(in)    :: Tr    !< Right-column layer tracer (conc, e.g. degC)
  real, dimension(nsurf),       intent(in)    :: PiL   !< Fractional position of neutral surface
                                                       !! within layer KoL of left column
  real, dimension(nsurf),       intent(in)    :: PiR   !< Fractional position of neutral surface
                                                       !! within layer KoR of right column
  integer, dimension(nsurf),    intent(in)    :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(nsurf),    intent(in)    :: KoR   !< Index of first right interface above neutral surface
  real, dimension(nsurf-1),     intent(in)    :: hEff  !< Effective thickness between two neutral
                                                       !! surfaces [H ~> m or kg m-2]
  real, dimension(nsurf-1),     intent(inout) :: Flx   !< Flux of tracer between pairs of neutral layers (conc H)
  logical,                      intent(in)    :: continuous !< True if using continuous reconstruction
  real,                         intent(in)    :: h_neglect !< A negligibly small width for the
                                             !! purpose of cell reconstructions [H ~> m or kg m-2]
  type(remapping_CS), optional, intent(in)    :: remap_CS !< Remapping control structure used
                                             !! to create sublayers
  real,               optional, intent(in)    :: h_neglect_edge !< A negligibly small width used for
                                             !! edge value calculations if continuous is false [H ~> m or kg m-2]
  ! Local variables
  integer :: k_sublayer, klb, klt, krb, krt, k
  real :: T_right_top, T_right_bottom, T_right_layer, T_right_sub, T_right_top_int, T_right_bot_int
  real :: T_left_top, T_left_bottom, T_left_layer, T_left_sub, T_left_top_int, T_left_bot_int
  real :: dT_top, dT_bottom, dT_layer, dT_ave, dT_sublayer, dT_top_int, dT_bot_int
  real, dimension(nk+1) :: Til !< Left-column interface tracer (conc, e.g. degC)
  real, dimension(nk+1) :: Tir !< Right-column interface tracer (conc, e.g. degC)
  real, dimension(nk) :: aL_l !< Left-column left edge value of tracer (conc, e.g. degC)
  real, dimension(nk) :: aR_l !< Left-column right edge value of tracer (conc, e.g. degC)
  real, dimension(nk) :: aL_r !< Right-column left edge value of tracer (conc, e.g. degC)
  real, dimension(nk) :: aR_r !< Right-column right edge value of tracer (conc, e.g. degC)
  ! Discontinuous reconstruction
  integer               :: iMethod
  real, dimension(nk,2) :: Tid_l !< Left-column interface tracer (conc, e.g. degC)
  real, dimension(nk,2) :: Tid_r !< Right-column interface tracer (conc, e.g. degC)
  real, dimension(nk,deg+1) :: ppoly_r_coeffs_l
  real, dimension(nk,deg+1) :: ppoly_r_coeffs_r
  real, dimension(nk,deg+1) :: ppoly_r_S_l
  real, dimension(nk,deg+1) :: ppoly_r_S_r
  logical :: down_flux
  ! Setup reconstruction edge values
  if (continuous) then
    call interface_scalar(nk, hl, Tl, Til, 2, h_neglect)
    call interface_scalar(nk, hr, Tr, Tir, 2, h_neglect)
    call ppm_left_right_edge_values(nk, Tl, Til, aL_l, aR_l)
    call ppm_left_right_edge_values(nk, Tr, Tir, aL_r, aR_r)
  else
    ppoly_r_coeffs_l(:,:) = 0.
    ppoly_r_coeffs_r(:,:) = 0.
    Tid_l(:,:) = 0.
    Tid_r(:,:) = 0.

    call build_reconstructions_1d( remap_CS, nk, hl, Tl, ppoly_r_coeffs_l, Tid_l, &
                                   ppoly_r_S_l, iMethod, h_neglect, h_neglect_edge )
    call build_reconstructions_1d( remap_CS, nk, hr, Tr, ppoly_r_coeffs_r, Tid_r, &
                                   ppoly_r_S_r, iMethod, h_neglect, h_neglect_edge )
  endif

  do k_sublayer = 1, nsurf-1
    if (hEff(k_sublayer) == 0.) then
      Flx(k_sublayer) = 0.
    else
      if (continuous) then
        klb = KoL(k_sublayer+1)
        T_left_bottom = ( 1. - PiL(k_sublayer+1) ) * Til(klb) + PiL(k_sublayer+1) * Til(klb+1)
        klt = KoL(k_sublayer)
        T_left_top = ( 1. - PiL(k_sublayer) ) * Til(klt) + PiL(k_sublayer) * Til(klt+1)
        T_left_layer = ppm_ave(PiL(k_sublayer), PiL(k_sublayer+1) + real(klb-klt), &
                               aL_l(klt), aR_l(klt), Tl(klt))

        krb = KoR(k_sublayer+1)
        T_right_bottom = ( 1. - PiR(k_sublayer+1) ) * Tir(krb) + PiR(k_sublayer+1) * Tir(krb+1)
        krt = KoR(k_sublayer)
        T_right_top = ( 1. - PiR(k_sublayer) ) * Tir(krt) + PiR(k_sublayer) * Tir(krt+1)
        T_right_layer = ppm_ave(PiR(k_sublayer), PiR(k_sublayer+1) + real(krb-krt), &
                                aL_r(krt), aR_r(krt), Tr(krt))
        dT_top = T_right_top - T_left_top
        dT_bottom = T_right_bottom - T_left_bottom
        dT_ave = 0.5 * ( dT_top + dT_bottom )
        dT_layer = T_right_layer - T_left_layer
        if (signum(1.,dT_top) * signum(1.,dT_bottom) <= 0. .or. signum(1.,dT_ave) * signum(1.,dT_layer) <= 0.) then
          dT_ave = 0.
        else
          dT_ave = dT_layer
        endif
        Flx(k_sublayer) = dT_ave * hEff(k_sublayer)
      else ! Discontinuous reconstruction
        ! Calculate tracer values on left and right side of the neutral surface
        call neutral_surface_T_eval(nk, nsurf, k_sublayer, KoL, PiL, Tl, Tid_l, deg, iMethod, &
                                    ppoly_r_coeffs_l, T_left_top, T_left_bottom, T_left_sub, &
                                    T_left_top_int, T_left_bot_int, T_left_layer)
        call neutral_surface_T_eval(nk, nsurf, k_sublayer, KoR, PiR, Tr, Tid_r, deg, iMethod, &
                                    ppoly_r_coeffs_r, T_right_top, T_right_bottom, T_right_sub, &
                                    T_right_top_int, T_right_bot_int, T_right_layer)

        dT_top      = T_right_top     - T_left_top
        dT_bottom   = T_right_bottom  - T_left_bottom
        dT_sublayer = T_right_sub     - T_left_sub
        dT_top_int  = T_right_top_int - T_left_top_int
        dT_bot_int  = T_right_bot_int - T_left_bot_int
        ! Enforcing the below criterion incorrectly zero out fluxes
        !dT_layer = T_right_layer - T_left_layer

        down_flux = dT_top <= 0. .and. dT_bottom <= 0. .and.       &
                    dT_sublayer <= 0. .and. dT_top_int <= 0. .and. &
                    dT_bot_int <= 0.
        down_flux = down_flux .or.                                 &
                    (dT_top >= 0. .and. dT_bottom >= 0. .and.      &
                    dT_sublayer >= 0. .and. dT_top_int >= 0. .and. &
                    dT_bot_int >= 0.)
        if (down_flux) then
          Flx(k_sublayer) = dT_sublayer * hEff(k_sublayer)
        else
          Flx(k_sublayer) = 0.
        endif
      endif
    endif
  enddo

end subroutine neutral_surface_flux

!> Evaluate various parts of the reconstructions to calculate gradient-based flux limter
subroutine neutral_surface_T_eval(nk, ns, k_sub, Ks, Ps, T_mean, T_int, deg, iMethod, T_poly, &
                                  T_top, T_bot, T_sub, T_top_int, T_bot_int, T_layer)
  integer,                   intent(in   ) :: nk        !< Number of cell everages
  integer,                   intent(in   ) :: ns        !< Number of neutral surfaces
  integer,                   intent(in   ) :: k_sub     !< Index of current neutral layer
  integer, dimension(ns),    intent(in   ) :: Ks        !< List of the layers associated with each neutral surface
  real, dimension(ns),       intent(in   ) :: Ps        !< List of the positions within a layer of each surface
  real, dimension(nk),       intent(in   ) :: T_mean    !< Cell average of tracer
  real, dimension(nk,2),     intent(in   ) :: T_int     !< Cell interface values of tracer from reconstruction
  integer,                   intent(in   ) :: deg       !< Degree of reconstruction polynomial (e.g. 1 is linear)
  integer,                   intent(in   ) :: iMethod   !< Method of integration to use
  real, dimension(nk,deg+1), intent(in   ) :: T_poly    !< Coefficients of polynomial reconstructions
  real,                      intent(  out) :: T_top     !< Tracer value at top (across discontinuity if necessary)
  real,                      intent(  out) :: T_bot     !< Tracer value at bottom (across discontinuity if necessary)
  real,                      intent(  out) :: T_sub     !< Average of the tracer value over the sublayer
  real,                      intent(  out) :: T_top_int !< Tracer value at top interface of neutral layer
  real,                      intent(  out) :: T_bot_int !< Tracer value at bottom interface of neutral layer
  real,                      intent(  out) :: T_layer   !< Cell-average that the the reconstruction belongs to

  integer :: kl, ks_top, ks_bot

  ks_top = k_sub
  ks_bot = k_sub + 1
  if ( Ks(ks_top) /= Ks(ks_bot) ) then
    call MOM_error(FATAL, "Neutral surfaces span more than one layer")
  endif
  kl = Ks(k_sub)
  ! First if the neutral surfaces spans the entirety of a cell, then do not search across the discontinuity
  if ( (Ps(ks_top) == 0.) .and. (Ps(ks_bot) == 1.)) then
    T_top = T_int(kl,1)
    T_bot = T_int(kl,2)
  else
    ! Search across potential discontinuity at top
    if ( (kl > 1) .and. (Ps(ks_top) == 0.)  ) then
      T_top = T_int(kl-1,2)
    else
      T_top = evaluation_polynomial( T_poly(kl,:), deg+1, Ps(ks_top) )
    endif
    ! Search across potential discontinuity at bottom
    if ( (kl < nk) .and. (Ps(ks_bot) == 1.) ) then
      T_bot = T_int(kl+1,1)
    else
      T_bot = evaluation_polynomial( T_poly(kl,:), deg+1, Ps(ks_bot) )
    endif
  endif
  T_sub = average_value_ppoly(nk, T_mean, T_int, T_poly, iMethod, kl, Ps(ks_top), Ps(ks_bot))
  T_top_int = evaluation_polynomial( T_poly(kl,:), deg+1, Ps(ks_top))
  T_bot_int = evaluation_polynomial( T_poly(kl,:), deg+1, Ps(ks_bot))
  T_layer = T_mean(kl)

end subroutine neutral_surface_T_eval

!> Discontinuous PPM reconstructions of the left/right edge values within a cell
subroutine ppm_left_right_edge_values(nk, Tl, Ti, aL, aR)
  integer,                    intent(in)    :: nk !< Number of levels
  real, dimension(nk),        intent(in)    :: Tl !< Layer tracer (conc, e.g. degC)
  real, dimension(nk+1),      intent(in)    :: Ti !< Interface tracer (conc, e.g. degC)
  real, dimension(nk),        intent(inout) :: aL !< Left edge value of tracer (conc, e.g. degC)
  real, dimension(nk),        intent(inout) :: aR !< Right edge value of tracer (conc, e.g. degC)

  integer :: k
  ! Setup reconstruction edge values
  do k = 1, nk
    aL(k) = Ti(k)
    aR(k) = Ti(k+1)
    if ( signum(1., aR(k) - Tl(k))*signum(1., Tl(k) - aL(k)) <= 0.0 ) then
      aL(k) = Tl(k)
      aR(k) = Tl(k)
    elseif ( sign(3., aR(k) - aL(k)) * ( (Tl(k) - aL(k)) + (Tl(k) - aR(k))) > abs(aR(k) - aL(k)) ) then
      aL(k) = Tl(k) + 2.0 * ( Tl(k) - aR(k) )
    elseif ( sign(3., aR(k) - aL(k)) * ( (Tl(k) - aL(k)) + (Tl(k) - aR(k))) < -abs(aR(k) - aL(k)) ) then
      aR(k) = Tl(k) + 2.0 * ( Tl(k) - aL(k) )
    endif
  enddo
end subroutine ppm_left_right_edge_values

!> Returns true if unit tests of neutral_diffusion functions fail. Otherwise returns false.
logical function neutral_diffusion_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout

  neutral_diffusion_unit_tests = .false. .or. &
    ndiff_unit_tests_continuous(verbose) .or. ndiff_unit_tests_discontinuous(verbose)

end function neutral_diffusion_unit_tests

!> Returns true if unit tests of neutral_diffusion functions fail. Otherwise returns false.
logical function ndiff_unit_tests_continuous(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables
  integer, parameter         :: nk = 4
  real, dimension(nk+1)      :: TiL, TiR1, TiR2, TiR4, Tio ! Test interface temperatures
  real, dimension(nk)        :: TL                         ! Test layer temperatures
  real, dimension(nk+1)      :: SiL                        ! Test interface salinities
  real, dimension(nk+1)      :: PiL, PiR4                  ! Test interface positions
  real, dimension(2*nk+2)    :: PiLRo, PiRLo               ! Test positions
  integer, dimension(2*nk+2) :: KoL, KoR                   ! Test indexes
  real, dimension(2*nk+1)    :: hEff                       ! Test positions
  real, dimension(2*nk+1)    :: Flx                        ! Test flux
  integer :: k
  logical :: v
  real :: h_neglect

  h_neglect = 1.0e-30

  v = verbose

  ndiff_unit_tests_continuous = .false. ! Normally return false
  write(stdout,*) '==== MOM_neutral_diffusion: ndiff_unit_tests_continuous ='

  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,1.,1.,1., 0.,1.,2., 1., 'FV: Straight line on uniform grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,1.,1.,0., 0.,4.,8., 7., 'FV: Vanished right cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,0.,1.,1., 0.,4.,8., 7., 'FV: Vanished left cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,1.,2.,4., 0.,3.,9., 4., 'FV: Stretched grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,2.,0.,2., 0.,1.,2., 0., 'FV: Vanished middle cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,0.,1.,0., 0.,1.,2., 2., 'FV: Vanished on both sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,1.,0.,0., 0.,1.,2., 0., 'FV: Two vanished cell sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,0.,0.,0., 0.,1.,2., 0., 'FV: All vanished cells')

  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,1.,1.,1., 0.,1.,2., 1., 'LSQ: Straight line on uniform grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,1.,1.,0., 0.,1.,2., 1., 'LSQ: Vanished right cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,0.,1.,1., 0.,1.,2., 1., 'LSQ: Vanished left cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,1.,2.,4., 0.,3.,9., 2., 'LSQ: Stretched grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,1.,0.,1., 0.,1.,2., 2., 'LSQ: Vanished middle cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,0.,1.,0., 0.,1.,2., 0., 'LSQ: Vanished on both sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,1.,0.,0., 0.,1.,2., 0., 'LSQ: Two vanished cell sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,0.,0.,0., 0.,1.,2., 0., 'LSQ: All vanished cells')

  call interface_scalar(4, (/10.,10.,10.,10./), (/24.,18.,12.,6./), Tio, 1, h_neglect)
  !ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
  !  test_data1d(5, Tio, (/27.,21.,15.,9.,3./), 'Linear profile, interface temperatures')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_data1d(v,5, Tio, (/24.,22.5,15.,7.5,6./), 'Linear profile, linear interface temperatures')
  call interface_scalar(4, (/10.,10.,10.,10./), (/24.,18.,12.,6./), Tio, 2, h_neglect)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_data1d(v,5, Tio, (/24.,22.,15.,8.,6./), 'Linear profile, PPM interface temperatures')

  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-1.0, 0.,  1.0, 1.0, 0.5, 'Check mid-point')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v, 0.0, 0.,  1.0, 1.0, 0.0, 'Check bottom')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v, 0.1, 0.,  1.1, 1.0, 0.0, 'Check below')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-1.0, 0.,  0.0, 1.0, 1.0, 'Check top')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-1.0, 0., -0.1, 1.0, 1.0, 'Check above')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-1.0, 0.,  3.0, 1.0, 0.25, 'Check 1/4')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-3.0, 0.,  1.0, 1.0, 0.75, 'Check 3/4')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v, 1.0, 0.,  1.0, 1.0, 0.0, 'Check dRho=0 below')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-1.0, 0., -1.0, 1.0, 1.0, 'Check dRho=0 above')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v, 0.0, 0.,  0.0, 1.0, 0.5, 'Check dRho=0 mid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-2.0, .5,  5.0, 0.5, 0.5, 'Check dP=0')

  ! Identical columns
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,2,2,3,3,3,3/), & ! KoL
                                   (/1,1,2,2,3,3,3,3/), & ! KoR
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pR
                                   (/0.,10.,0.,10.,0.,10.,0./), & ! hEff
                                   'Identical columns')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 8, &
                                   absolute_positions(3, 8, (/0.,10.,20.,30./), KoL, PiLRo), &
                                   (/0.,0.,10.,10.,20.,20.,30.,30./), '... left positions')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 8, &
                                   absolute_positions(3, 8, (/0.,10.,20.,30./), KoR, PiRLo), &
                                   (/0.,0.,10.,10.,20.,20.,30.,30./), '... right positions')
  call neutral_surface_flux(3, 2*3+2, 2, (/10.,10.,10./), (/10.,10.,10./), & ! nk, hL, hR
                               (/20.,16.,12./), (/20.,16.,12./), & ! Tl, Tr
                               PiLRo, PiRLo, KoL, KoR, hEff, Flx, .true., h_neglect)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 7, Flx, &
              (/0.,0.,0.,0.,0.,0.,0./), 'Identical columns, rho flux (=0)')
  call neutral_surface_flux(3, 2*3+2, 2, (/10.,10.,10./), (/10.,10.,10./), & ! nk, hL, hR
                               (/-1.,-1.,-1./), (/1.,1.,1./), & ! Sl, Sr
                               PiLRo, PiRLo, KoL, KoR, hEff, Flx, .true., h_neglect)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 7, Flx, &
              (/0.,20.,0.,20.,0.,20.,0./), 'Identical columns, S flux')

  ! Right column slightly cooler than left
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/20.,16.,12.,8./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,2,2,3,3,3,3/), & ! kL
                                   (/1,1,1,2,2,3,3,3/), & ! kR
                                   (/0.,0.5,0.,0.5,0.,0.5,1.,1./), & ! pL
                                   (/0.,0.,0.5,0.,0.5,0.,0.5,1./), & ! pR
                                   (/0.,5.,5.,5.,5.,5.,0./), & ! hEff
                                   'Right column slightly cooler')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 8, &
                                   absolute_positions(3, 8, (/0.,10.,20.,30./), KoL, PiLRo), &
                                   (/0.,5.,10.,15.,20.,25.,30.,30./), '... left positions')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 8, &
                                   absolute_positions(3, 8, (/0.,10.,20.,30./), KoR, PiRLo), &
                                   (/0.,0.,5.,10.,15.,20.,25.,30./), '... right positions')

  ! Right column slightly warmer than left
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/24.,20.,16.,12./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,1,2,2,3,3,3/), & ! kL
                                   (/1,1,2,2,3,3,3,3/), & ! kR
                                   (/0.,0.,0.5,0.,0.5,0.,0.5,1./), & ! pL
                                   (/0.,0.5,0.,0.5,0.,0.5,1.,1./), & ! pR
                                   (/0.,5.,5.,5.,5.,5.,0./), & ! hEff
                                   'Right column slightly warmer')

  ! Right column somewhat cooler than left
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/16.,12.,8.,4./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,2,3,3,3,3,3/), & ! kL
                                   (/1,1,1,1,2,2,3,3/), & ! kR
                                   (/0.,0.,0.5,0.,0.5,1.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.5,0.,0.5,0.,1./), & ! pR
                                   (/0.,0.,5.,5.,5.,0.,0./), & ! hEff
                                   'Right column somewhat cooler')

  ! Right column much colder than left with no overlap
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/9.,7.,5.,3./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,3,3,3,3,3,3/), & ! kL
                                   (/1,1,1,1,1,2,3,3/), & ! kR
                                   (/0.,0.,0.,1.,1.,1.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,0.,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,0.,0./), & ! hEff
                                   'Right column much cooler')

  ! Right column with mixed layer
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,3,3,3,3,3,3/), & ! kL
                                   (/1,1,1,1,2,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,1.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,0.,1./), & ! pR
                                   (/0.,0.,0.,0.,10.,0.,0./), & ! hEff
                                   'Right column with mixed layer')

  ! Identical columns with mixed layer
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,2,2,3,3,3,3/), & ! kL
                                   (/1,1,2,2,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pR
                                   (/0.,10.,0.,10.,0.,10.,0./), & ! hEff
                                   'Identical columns with mixed layer')

  ! Right column with unstable mixed layer
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/10.,14.,12.,4./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,3,3,3,3,3,3/), & ! kL
                                   (/1,1,1,2,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.,.75,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.25,1.,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,7.5,0./), & ! hEff
                                   'Right column with unstable mixed layer')

  ! Left column with unstable mixed layer
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/10.,14.,12.,4./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,1,2,3,3,3,3/), & ! kL
                                   (/1,2,3,3,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.25,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,.75,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,7.5,0./), & ! hEff
                                   'Left column with unstable mixed layer')

  ! Two unstable mixed layers
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/8.,12.,10.,2./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/10.,14.,12.,4./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,1,1,2,3,3,3/), & ! kL
                                   (/1,2,3,3,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.,0.75,1./), & ! pL
                                   (/0.,0.,0.,0.5,0.5,0.5,1.,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,6.,0./), & ! hEff
                                   'Two unstable mixed layers')

  if (.not. ndiff_unit_tests_continuous) write(stdout,*) 'Pass'

end function ndiff_unit_tests_continuous

logical function ndiff_unit_tests_discontinuous(verbose)
  logical, intent(in) :: verbose !< It true, write results to stdout
  ! Local variables
  integer, parameter          :: nk = 3
  integer, parameter          :: ns = nk*4
  real, dimension(nk)         :: Sl, Sr, Tl, Tr ! Salinities [ppt] and temperatures [degC]
  real, dimension(nk)         :: hl, hr    ! Thicknesses in pressure units [R L2 T-2 ~> Pa]
  real, dimension(nk,2)       :: TiL, SiL, TiR, SiR ! Cell edge salinities [ppt] and temperatures [degC]
  real, dimension(nk,2)       :: Pres_l, Pres_r ! Interface pressures [R L2 T-2 ~> Pa]
  integer, dimension(ns)      :: KoL, KoR
  real, dimension(ns)         :: PoL, PoR
  real, dimension(ns-1)       :: hEff, Flx
  type(neutral_diffusion_CS)  :: CS        !< Neutral diffusion control structure
  type(EOS_type),     pointer :: EOS       !< Structure for linear equation of state
  type(remapping_CS), pointer :: remap_CS  !< Remapping control structure (PLM)
  real, dimension(nk,2)       :: ppoly_T_l, ppoly_T_r ! Linear reconstruction for T
  real, dimension(nk,2)       :: ppoly_S_l, ppoly_S_r ! Linear reconstruction for S
  real, dimension(nk,2)       :: dRdT      !< Partial derivative of density with temperature at
                                           !! cell edges [R degC-1 ~> kg m-3 degC-1]
  real, dimension(nk,2)       :: dRdS      !< Partial derivative of density with salinity at
                                           !! cell edges [R ppt-1 ~> kg m-3 ppt-1]
  logical, dimension(nk)      :: stable_l, stable_r
  integer                     :: iMethod
  integer                     :: ns_l, ns_r
  integer :: k
  logical :: v

  v = verbose
  ndiff_unit_tests_discontinuous = .false. ! Normally return false
  write(stdout,*) '==== MOM_neutral_diffusion: ndiff_unit_tests_discontinuous ='

  ! Unit tests for find_neutral_surface_positions_discontinuous
  ! Salinity is 0 for all these tests
  allocate(CS%EOS)
  call EOS_manual_init(CS%EOS, form_of_EOS=EOS_LINEAR, dRho_dT=-1., dRho_dS=0.)
  Sl(:) = 0. ; Sr(:) = 0. ; ; SiL(:,:) = 0. ; SiR(:,:) = 0.
  ppoly_T_l(:,:) = 0.; ppoly_T_r(:,:) = 0.
  ppoly_S_l(:,:) = 0.; ppoly_S_r(:,:) = 0.
  ! Intialize any control structures needed for unit tests
  CS%ref_pres = -1.

  hL = (/10.,10.,10./) ; hR = (/10.,10.,10./)
  Pres_l(1,1) = 0. ; Pres_l(1,2) = hL(1) ; Pres_r(1,1) = 0. ; Pres_r(1,2) = hR(1)
  do k = 2,nk
    Pres_l(k,1) = Pres_l(k-1,2)
    Pres_l(k,2) = Pres_l(k,1) + hL(k)
    Pres_r(k,1) = Pres_r(k-1,2)
    Pres_r(k,2) = Pres_r(k,1) + hR(k)
  enddo
  CS%delta_rho_form = 'mid_pressure'
  CS%neutral_pos_method = 1

  TiL(1,:) = (/ 22.00, 18.00 /); TiL(2,:) = (/ 18.00, 14.00 /); TiL(3,:) = (/ 14.00, 10.00 /);
  TiR(1,:) = (/ 22.00, 18.00 /); TiR(2,:) = (/ 18.00, 14.00 /); TiR(3,:) = (/ 14.00, 10.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 1.00, 1.00 /),  & ! PoL
    (/ 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00, 1.00 /),  & ! PoR
    (/ 0.00, 10.00, 0.00, 0.00, 0.00, 10.00, 0.00, 0.00, 0.00, 10.00, 0.00 /),  & ! hEff
    'Identical Columns')

  TiL(1,:) = (/ 22.00, 18.00 /); TiL(2,:) = (/ 18.00, 14.00 /); TiL(3,:) = (/ 14.00, 10.00 /);
  TiR(1,:) = (/ 20.00, 16.00 /); TiR(2,:) = (/ 16.00, 12.00 /); TiR(3,:) = (/ 12.00, 8.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 1.00 /),  & ! PoL
    (/ 0.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 1.00 /),  & ! PoR
    (/ 0.00, 5.00, 0.00, 5.00, 0.00, 5.00, 0.00, 5.00, 0.00, 5.00, 0.00 /),  & ! hEff
    'Right slightly cooler')

  TiL(1,:) = (/ 20.00, 16.00 /); TiL(2,:) = (/ 16.00, 12.00 /); TiL(3,:) = (/ 12.00, 8.00 /);
  TiR(1,:) = (/ 22.00, 18.00 /); TiR(2,:) = (/ 18.00, 14.00 /); TiR(3,:) = (/ 14.00, 10.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 1.00 /),  & ! PoL
    (/ 0.00, 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 1.00 /),  & ! PoR
    (/ 0.00, 5.00, 0.00, 5.00, 0.00, 5.00, 0.00, 5.00, 0.00, 5.00, 0.00 /),  & ! hEff
    'Left slightly cooler')

  TiL(1,:) = (/ 22.00, 20.00 /); TiL(2,:) = (/ 18.00, 16.00 /); TiL(3,:) = (/ 14.00, 12.00 /);
  TiR(1,:) = (/ 32.00, 24.00 /); TiR(2,:) = (/ 22.00, 14.00 /); TiR(3,:) = (/ 12.00, 4.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 1.00, 1.00, 1.00 /),  & ! PoL
    (/ 0.00, 1.00, 0.00, 0.00, 0.25, 0.50, 0.75, 1.00, 0.00, 0.00, 0.00, 1.00 /),  & ! PoR
    (/ 0.00, 0.00, 0.00, 4.00, 0.00, 4.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Right more strongly stratified')

  TiL(1,:) = (/ 22.00, 18.00 /); TiL(2,:) = (/ 18.00, 14.00 /); TiL(3,:) = (/ 14.00, 10.00 /);
  TiR(1,:) = (/ 14.00, 14.00 /); TiR(2,:) = (/ 14.00, 14.00 /); TiR(3,:) = (/ 12.00, 8.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.50, 1.00, 1.00 /),  & ! PoL
    (/ 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 1.00 /),  & ! PoR
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 5.00, 0.00 /),  & ! hEff
    'Deep Mixed layer on the right')

  TiL(1,:) = (/ 14.00, 14.00 /); TiL(2,:) = (/ 14.00, 12.00 /); TiL(3,:) = (/ 10.00, 8.00 /);
  TiR(1,:) = (/ 14.00, 14.00 /); TiR(2,:) = (/ 14.00, 14.00 /); TiR(3,:) = (/ 14.00, 14.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3 /),  & ! KoL
    (/ 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /),  & ! PoL
    (/ 0.00, 0.00, 0.00, 1.00, 0.00, 1.00, 0.00, 1.00, 1.00, 1.00, 1.00, 1.00 /),  & ! PoR
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Right unstratified column')

  TiL(1,:) = (/ 14.00, 14.00 /); TiL(2,:) = (/ 14.00, 12.00 /); TiL(3,:) = (/ 10.00, 8.00 /);
  TiR(1,:) = (/ 14.00, 14.00 /); TiR(2,:) = (/ 14.00, 14.00 /); TiR(3,:) = (/ 12.00, 4.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.00, 1.00, 1.00, 0.00, 1.00, 1.00 /),  & ! PoL
    (/ 0.00, 0.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.25, 0.50, 1.00 /),  & ! PoR
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 4.00, 0.00 /),  & ! hEff
    'Right unstratified column')

  TiL(1,:) = (/ 14.00, 14.00 /); TiL(2,:) = (/ 14.00, 10.00 /); TiL(3,:) = (/ 10.00, 2.00 /);
  TiR(1,:) = (/ 14.00, 14.00 /); TiR(2,:) = (/ 14.00, 10.00 /); TiR(3,:) = (/ 10.00, 2.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 1.00, 1.00, 1.00, 0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 1.00, 1.00 /),  & ! PoL
    (/ 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 1.00, 1.00 /),  & ! PoR
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 10.00, 0.00, 0.00, 0.00, 10.00, 0.00 /),  & ! hEff
    'Identical columns with mixed layer')

  TiL(1,:) = (/ 14.00, 12.00 /); TiL(2,:) = (/ 10.00, 10.00 /); TiL(3,:) = (/ 8.00, 2.00 /);
  TiR(1,:) = (/ 14.00, 12.00 /); TiR(2,:) = (/ 12.00, 8.00 /); TiR(3,:) = (/ 8.00, 2.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 0.00, 1.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 1.00, 1.00 /),  & ! PoL
    (/ 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.00, 1.00, 1.00 /),  & ! PoR
    (/ 0.00, 10.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 10.00, 0.00 /),  & ! hEff
    'Left interior unstratified')

  TiL(1,:) = (/ 12.00, 12.00 /); TiL(2,:) = (/ 12.00, 10.00 /); TiL(3,:) = (/ 10.00, 6.00 /);
  TiR(1,:) = (/ 12.00, 10.00 /); TiR(2,:) = (/ 10.00, 12.00 /); TiR(3,:) = (/ 8.00, 4.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 1.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.50, 1.00, 1.00 /),  & ! PoL
    (/ 0.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.50, 1.00 /),  & ! PoR
    (/ 0.00, 0.00, 0.00, 10.00, 0.00, 0.00, 0.00, 0.00, 0.00, 5.00, 0.00 /),  & ! hEff
    'Left mixed layer, Right unstable interior')

  TiL(1,:) = (/ 14.00, 14.00 /); TiL(2,:) = (/ 10.00, 10.00 /); TiL(3,:) = (/ 8.00, 6.00 /);
  TiR(1,:) = (/ 10.00, 14.00 /); TiR(2,:) = (/ 16.00, 16.00 /); TiR(3,:) = (/ 12.00, 4.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 1.00, 0.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.00, 0.00, 1.00, 1.00 /),  & ! PoL
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.50, 0.75, 1.00 /),  & ! PoR
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 4.00, 0.00 /),  & ! hEff
    'Left thick mixed layer, Right unstable mixed')

  TiL(1,:) = (/ 8.00, 12.00 /); TiL(2,:) = (/ 12.00, 10.00 /); TiL(3,:) = (/ 8.00, 4.00 /);
  TiR(1,:) = (/ 10.00, 14.00 /); TiR(2,:) = (/ 14.00, 12.00 /); TiR(3,:) = (/ 10.00, 6.00 /);
  call mark_unstable_cells( CS, nk, Til, Sil, Pres_l, stable_l )
  call mark_unstable_cells( CS, nk, Tir, Sir, Pres_r, stable_r )
  call find_neutral_surface_positions_discontinuous(CS, nk, Pres_l, hL, TiL, SiL, ppoly_T_l, ppoly_S_l, stable_l, &
           Pres_r, hR, TiR, SiR, ppoly_T_r, ppoly_S_r, stable_r, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 /),  & ! KoL
    (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3 /),  & ! KoR
    (/ 0.00, 1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.50, 1.00 /),  & ! PoL
    (/ 0.00, 0.00, 0.00, 1.00, 0.00, 1.00, 1.00, 0.00, 0.00, 0.50, 1.00, 1.00 /),  & ! PoR
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 5.00, 0.00 /),  & ! hEff
    'Unstable mixed layers, left cooler')

  call EOS_manual_init(CS%EOS, form_of_EOS = EOS_LINEAR, dRho_dT = -1., dRho_dS = 2.)
  ! Tests for linearized version of searching the layer for neutral surface position
  ! EOS linear in T, uniform alpha
  CS%max_iter = 10
  ! Unit tests require explicit initialization of tolerance
  CS%Drho_tol = 0.
  CS%x_tol = 0.
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5, &
             find_neutral_pos_linear(CS, 0., 10., 35., -0.2, 0., &
                                     -0.2, 0., -0.2, 0.,                     &
                                     (/12.,-4./), (/34.,0./)), "Temp Uniform Linearized Alpha/Beta"))
  ! EOS linear in S, uniform beta
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5, &
             find_neutral_pos_linear(CS, 0., 10., 35., 0., 0.8, &
                                     0., 0.8, 0., 0.8,                &
                                    (/12.,0./), (/34.,2./)), "Salt Uniform Linearized Alpha/Beta"))
  ! EOS linear in T/S, uniform alpha/beta
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5,   &
             find_neutral_pos_linear(CS, 0., 10., 35., -0.5, 0.5,                &
                                     -0.5, 0.5, -0.5, 0.5,  &
                                     (/12.,-4./), (/34.,2./)), "Temp/salt Uniform Linearized Alpha/Beta"))
  ! EOS linear in T, insensitive to So
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5, &
             find_neutral_pos_linear(CS, 0., 10., 35., -0.2, 0., &
                                     -0.4, 0., -0.6, 0.,  &
                                     (/12.,-4./), (/34.,0./)), "Temp stratified Linearized Alpha/Beta"))
  ! EOS linear in S, insensitive to T
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5, &
             find_neutral_pos_linear(CS, 0., 10., 35., 0., 0.8,  &
                                      0., 1.0,  0., 0.5,  &
                                     (/12.,0./), (/34.,2./)), "Salt stratified Linearized Alpha/Beta"))
  if (.not. ndiff_unit_tests_discontinuous) write(stdout,*) 'Pass'

end function ndiff_unit_tests_discontinuous

!> Returns true if a test of fv_diff() fails, and conditionally writes results to stream
logical function test_fv_diff(verbose, hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
  logical,          intent(in) :: verbose !< If true, write results to stdout
  real,             intent(in) :: hkm1  !< Left cell width [nondim]
  real,             intent(in) :: hk    !< Center cell width [nondim]
  real,             intent(in) :: hkp1  !< Right cell width [nondim]
  real,             intent(in) :: Skm1  !< Left cell average value
  real,             intent(in) :: Sk    !< Center cell average value
  real,             intent(in) :: Skp1  !< Right cell average value
  real,             intent(in) :: Ptrue !< True answer [nondim]
  character(len=*), intent(in) :: title !< Title for messages

  ! Local variables
  integer :: stdunit
  real :: Pret

  Pret = fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  test_fv_diff = (Pret /= Ptrue)

  if (test_fv_diff .or. verbose) then
    stdunit = stdout
    if (test_fv_diff) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    if (test_fv_diff) then
      write(stdunit,'(2(x,a,f20.16),x,a)') 'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
    else
      write(stdunit,'(2(x,a,f20.16))') 'pRet=',Pret,'pTrue=',Ptrue
    endif
  endif

end function test_fv_diff

!> Returns true if a test of fvlsq_slope() fails, and conditionally writes results to stream
logical function test_fvlsq_slope(verbose, hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
  logical,          intent(in) :: verbose !< If true, write results to stdout
  real,             intent(in) :: hkm1  !< Left cell width
  real,             intent(in) :: hk    !< Center cell width
  real,             intent(in) :: hkp1  !< Right cell width
  real,             intent(in) :: Skm1  !< Left cell average value
  real,             intent(in) :: Sk    !< Center cell average value
  real,             intent(in) :: Skp1  !< Right cell average value
  real,             intent(in) :: Ptrue !< True answer
  character(len=*), intent(in) :: title !< Title for messages

  ! Local variables
  integer :: stdunit
  real :: Pret

  Pret = fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  test_fvlsq_slope = (Pret /= Ptrue)

  if (test_fvlsq_slope .or. verbose) then
    stdunit = stdout
    if (test_fvlsq_slope) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    if (test_fvlsq_slope) then
      write(stdunit,'(2(x,a,f20.16),x,a)') 'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
    else
      write(stdunit,'(2(x,a,f20.16))') 'pRet=',Pret,'pTrue=',Ptrue
    endif
  endif

end function test_fvlsq_slope

!> Returns true if a test of interpolate_for_nondim_position() fails, and conditionally writes results to stream
logical function test_ifndp(verbose, rhoNeg, Pneg, rhoPos, Ppos, Ptrue, title)
  logical,          intent(in) :: verbose !< If true, write results to stdout
  real,             intent(in) :: rhoNeg !< Lighter density [R ~> kg m-3]
  real,             intent(in) :: Pneg   !< Interface position of lighter density [nondim]
  real,             intent(in) :: rhoPos !< Heavier density [R ~> kg m-3]
  real,             intent(in) :: Ppos   !< Interface position of heavier density [nondim]
  real,             intent(in) :: Ptrue  !< True answer [nondim]
  character(len=*), intent(in) :: title  !< Title for messages

  ! Local variables
  integer :: stdunit
  real :: Pret

  Pret = interpolate_for_nondim_position(rhoNeg, Pneg, rhoPos, Ppos)
  test_ifndp = (Pret /= Ptrue)

  if (test_ifndp .or. verbose) then
    stdunit = stdout
    if (test_ifndp) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    if (test_ifndp) then
      write(stdunit,'(4(x,a,f20.16),2(x,a,1pe22.15),x,a)') &
            'r1=',rhoNeg,'p1=',Pneg,'r2=',rhoPos,'p2=',Ppos,'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
    else
      write(stdunit,'(4(x,a,f20.16),2(x,a,1pe22.15))') &
            'r1=',rhoNeg,'p1=',Pneg,'r2=',rhoPos,'p2=',Ppos,'pRet=',Pret,'pTrue=',Ptrue
    endif
  endif

end function test_ifndp

!> Returns true if comparison of Po and Ptrue fails, and conditionally writes results to stream
logical function test_data1d(verbose, nk, Po, Ptrue, title)
  logical,             intent(in) :: verbose !< If true, write results to stdout
  integer,             intent(in) :: nk    !< Number of layers
  real, dimension(nk), intent(in) :: Po    !< Calculated answer
  real, dimension(nk), intent(in) :: Ptrue !< True answer
  character(len=*),    intent(in) :: title !< Title for messages

  ! Local variables
  integer :: k, stdunit

  test_data1d = .false.
  do k = 1,nk
    if (Po(k) /= Ptrue(k)) test_data1d = .true.
  enddo

  if (test_data1d .or. verbose) then
    stdunit = stdout
    if (test_data1d) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    do k = 1,nk
      if (Po(k) /= Ptrue(k)) then
        test_data1d = .true.
        write(stdunit,'(a,i2,2(x,a,f20.16),x,a,1pe22.15,x,a)') &
              'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k),'WRONG!'
      else
        if (verbose) &
          write(stdunit,'(a,i2,2(x,a,f20.16),x,a,1pe22.15)') &
                'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k)
      endif
    enddo
  endif

end function test_data1d

!> Returns true if comparison of Po and Ptrue fails, and conditionally writes results to stream
logical function test_data1di(verbose, nk, Po, Ptrue, title)
  logical,                intent(in) :: verbose !< If true, write results to stdout
  integer,                intent(in) :: nk    !< Number of layers
  integer, dimension(nk), intent(in) :: Po    !< Calculated answer
  integer, dimension(nk), intent(in) :: Ptrue !< True answer
  character(len=*),       intent(in) :: title !< Title for messages

  ! Local variables
  integer :: k, stdunit

  test_data1di = .false.
  do k = 1,nk
    if (Po(k) /= Ptrue(k)) test_data1di = .true.
  enddo

  if (test_data1di .or. verbose) then
    stdunit = stdout
    if (test_data1di) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    do k = 1,nk
      if (Po(k) /= Ptrue(k)) then
        test_data1di = .true.
        write(stdunit,'(a,i2,2(x,a,i5),x,a)') 'k=',k,'Io=',Po(k),'Itrue=',Ptrue(k),'WRONG!'
      else
        if (verbose) &
          write(stdunit,'(a,i2,2(x,a,i5))') 'k=',k,'Io=',Po(k),'Itrue=',Ptrue(k)
      endif
    enddo
  endif

end function test_data1di

!> Returns true if output of find_neutral_surface_positions() does not match correct values,
!! and conditionally writes results to stream
logical function test_nsp(verbose, ns, KoL, KoR, pL, pR, hEff, KoL0, KoR0, pL0, pR0, hEff0, title)
  logical,                intent(in) :: verbose !< If true, write results to stdout
  integer,                intent(in) :: ns    !< Number of surfaces
  integer, dimension(ns), intent(in) :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(ns), intent(in) :: KoR   !< Index of first right interface above neutral surface
  real, dimension(ns),    intent(in) :: pL    !< Fractional position of neutral surface within layer KoL of left column
  real, dimension(ns),    intent(in) :: pR    !< Fractional position of neutral surface within layer KoR of right column
  real, dimension(ns-1),  intent(in) :: hEff  !< Effective thickness between two neutral surfaces [R L2 T-2 ~> Pa]
  integer, dimension(ns), intent(in) :: KoL0  !< Correct value for KoL
  integer, dimension(ns), intent(in) :: KoR0  !< Correct value for KoR
  real, dimension(ns),    intent(in) :: pL0   !< Correct value for pL
  real, dimension(ns),    intent(in) :: pR0   !< Correct value for pR
  real, dimension(ns-1),  intent(in) :: hEff0 !< Correct value for hEff
  character(len=*),       intent(in) :: title !< Title for messages

  ! Local variables
  integer :: k, stdunit
  logical :: this_row_failed

  test_nsp = .false.
  do k = 1,ns
    test_nsp = test_nsp .or. compare_nsp_row(KoL(k), KoR(k), pL(k), pR(k), KoL0(k), KoR0(k), pL0(k), pR0(k))
    if (k < ns) then
      if (hEff(k) /= hEff0(k)) test_nsp = .true.
    endif
  enddo

  if (test_nsp .or. verbose) then
    stdunit = stdout
    if (test_nsp) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    do k = 1,ns
      this_row_failed = compare_nsp_row(KoL(k), KoR(k), pL(k), pR(k), KoL0(k), KoR0(k), pL0(k), pR0(k))
      if (this_row_failed) then
        write(stdunit,10) k,KoL(k),pL(k),KoR(k),pR(k),' <-- WRONG!'
        write(stdunit,10) k,KoL0(k),pL0(k),KoR0(k),pR0(k),' <-- should be this'
      else
        write(stdunit,10) k,KoL(k),pL(k),KoR(k),pR(k)
      endif
      if (k < ns) then
        if (hEff(k) /= hEff0(k)) then
          write(stdunit,'(i3,8x,"layer hEff =",2(f20.16,a))') k,hEff(k)," .neq. ",hEff0(k),' <-- WRONG!'
        else
          write(stdunit,'(i3,8x,"layer hEff =",f20.16)') k,hEff(k)
        endif
      endif
    enddo
  endif
  if (test_nsp) call MOM_error(FATAL,"test_nsp failed")

10 format("ks=",i3," kL=",i3," pL=",f20.16," kR=",i3," pR=",f20.16,a)
end function test_nsp

!> Compares a single row, k, of output from find_neutral_surface_positions()
logical function compare_nsp_row(KoL, KoR, pL, pR, KoL0, KoR0, pL0, pR0)
  integer,  intent(in) :: KoL   !< Index of first left interface above neutral surface
  integer,  intent(in) :: KoR   !< Index of first right interface above neutral surface
  real,     intent(in) :: pL    !< Fractional position of neutral surface within layer KoL of left column
  real,     intent(in) :: pR    !< Fractional position of neutral surface within layer KoR of right column
  integer,  intent(in) :: KoL0  !< Correct value for KoL
  integer,  intent(in) :: KoR0  !< Correct value for KoR
  real,     intent(in) :: pL0   !< Correct value for pL
  real,     intent(in) :: pR0   !< Correct value for pR

  compare_nsp_row = .false.
  if (KoL /= KoL0) compare_nsp_row = .true.
  if (KoR /= KoR0) compare_nsp_row = .true.
  if (pL /= pL0) compare_nsp_row = .true.
  if (pR /= pR0) compare_nsp_row = .true.
end function compare_nsp_row

!> Compares output position from refine_nondim_position with an expected value
logical function test_rnp(expected_pos, test_pos, title)
  real,             intent(in) :: expected_pos !< The expected position
  real,             intent(in) :: test_pos !< The position returned by the code
  character(len=*), intent(in) :: title    !< A label for this test
  ! Local variables
  integer :: stdunit

  stdunit = stdout ! Output to standard error
  test_rnp = ABS(expected_pos - test_pos) > 2*EPSILON(test_pos)
  if (test_rnp) then
    write(stdunit,'(A, f20.16, " .neq. ", f20.16, " <-- WRONG")') title, expected_pos, test_pos
  else
    write(stdunit,'(A, f20.16, " ==  ", f20.16)') title, expected_pos, test_pos
  endif
end function test_rnp
!> Deallocates neutral_diffusion control structure
subroutine neutral_diffusion_end(CS)
  type(neutral_diffusion_CS), pointer :: CS  !< Neutral diffusion control structure

  if (associated(CS)) deallocate(CS)

end subroutine neutral_diffusion_end

end module MOM_neutral_diffusion
