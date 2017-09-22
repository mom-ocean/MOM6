!> A column-wise toolbox for implementing neutral diffusion
module MOM_neutral_diffusion

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,        only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,    only : diag_ctrl, time_type
use MOM_diag_mediator,    only : post_data, register_diag_field
use MOM_EOS,              only : EOS_type, EOS_manual_init, calculate_compress, calculate_density_derivs
use MOM_EOS,              only : calculate_density_second_derivs
use MOM_EOS,              only : extract_member_EOS, EOS_LINEAR, EOS_TEOS10, EOS_WRIGHT
use MOM_error_handler,    only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,      only : get_param, log_version, param_file_type
use MOM_file_parser,      only : openParameterBlock, closeParameterBlock
use MOM_grid,             only : ocean_grid_type
use MOM_remapping,        only : remapping_CS, initialize_remapping
use MOM_remapping,        only : extract_member_remapping_CS, build_reconstructions_1d
use MOM_remapping,        only : average_value_ppoly, remappingSchemesDoc, remappingDefaultScheme
use MOM_tracer_registry,  only : tracer_registry_type
use MOM_verticalGrid,     only : verticalGrid_type
use polynomial_functions, only : evaluation_polynomial, first_derivative_polynomial
use PPM_functions,        only : PPM_reconstruction, PPM_boundary_extrapolation
use regrid_edge_values,   only : edge_values_implicit_h4

implicit none ; private

#include <MOM_memory.h>

public neutral_diffusion
public neutral_diffusion_init
public neutral_diffusion_diag_init
public neutral_diffusion_end
public neutral_diffusion_calc_coeffs
public neutral_diffusion_unit_tests

type, public :: neutral_diffusion_CS ; private
  integer :: nkp1   ! Number of interfaces for a column = nk + 1
  integer :: nsurf  ! Number of neutral surfaces
  integer :: ppoly_deg = 2 ! Degree of polynomial used for reconstructions
  logical :: continuous_reconstruction = .true.   ! True if using continuous PPM reconstruction at interfaces
  logical :: boundary_extrap = .true.
  logical :: refine_position = .false.
  integer :: max_iter ! Maximum number of iterations if refine_position is defined
  real :: tolerance   ! Convergence criterion representing difference from true neutrality
  real :: ref_pres    ! Reference pressure, negative if using locally referenced neutral density

  ! Positions of neutral surfaces in both the u, v directions
  real,    allocatable, dimension(:,:,:) :: uPoL  ! Non-dimensional position with left layer uKoL-1, u-point
  real,    allocatable, dimension(:,:,:) :: uPoR  ! Non-dimensional position with right layer uKoR-1, u-point
  integer, allocatable, dimension(:,:,:) :: uKoL  ! Index of left interface corresponding to neutral surface, u-point
  integer, allocatable, dimension(:,:,:) :: uKoR  ! Index of right interface corresponding to neutral surface, u-point
  real,    allocatable, dimension(:,:,:) :: uHeff ! Effective thickness at u-point (H units)
  real,    allocatable, dimension(:,:,:) :: vPoL  ! Non-dimensional position with left layer uKoL-1, v-point
  real,    allocatable, dimension(:,:,:) :: vPoR  ! Non-dimensional position with right layer uKoR-1, v-point
  integer, allocatable, dimension(:,:,:) :: vKoL  ! Index of left interface corresponding to neutral surface, v-point
  integer, allocatable, dimension(:,:,:) :: vKoR  ! Index of right interface corresponding to neutral surface, v-point
  real,    allocatable, dimension(:,:,:) :: vHeff ! Effective thickness at v-point (H units)
  ! Coefficients of polynomial reconstructions for temperature and salinity
  real,    allocatable, dimension(:,:,:,:) :: ppoly_coeffs_T !< Polynomial coefficients for temperature
  real,    allocatable, dimension(:,:,:,:) :: ppoly_coeffs_S !< Polynomial coefficients for temperature
  ! Variables needed for continuous reconstructions
  real,    allocatable, dimension(:,:,:) :: dRdT ! dRho/dT (kg/m3/degC) at interfaces
  real,    allocatable, dimension(:,:,:) :: dRdS ! dRho/dS (kg/m3/ppt) at interfaces
  real,    allocatable, dimension(:,:,:) :: Tint ! Interface T (degC)
  real,    allocatable, dimension(:,:,:) :: Sint ! Interface S (ppt)
  real,    allocatable, dimension(:,:,:) :: Pint ! Interface pressure (Pa)
  ! Variables needed for discontinuous reconstructions
  real,    allocatable, dimension(:,:,:,:) :: T_i     ! Top edge reconstruction of temperature (degC)
  real,    allocatable, dimension(:,:,:,:) :: S_i     ! Top edge reconstruction of salinity (ppt)
  real,    allocatable, dimension(:,:,:,:) :: dRdT_i     ! dRho/dT (kg/m3/degC) at top edge
  real,    allocatable, dimension(:,:,:,:) :: dRdS_i     ! dRho/dS (kg/m3/ppt) at top edge

  type(diag_ctrl), pointer :: diag ! structure to regulate output
  integer, allocatable, dimension(:) :: id_neutral_diff_tracer_conc_tend    ! tracer concentration tendency
  integer, allocatable, dimension(:) :: id_neutral_diff_tracer_cont_tend    ! tracer content tendency
  integer, allocatable, dimension(:) :: id_neutral_diff_tracer_cont_tend_2d ! k-summed tracer content tendency
  integer, allocatable, dimension(:) :: id_neutral_diff_tracer_trans_x_2d   ! k-summed ndiff zonal tracer transport
  integer, allocatable, dimension(:) :: id_neutral_diff_tracer_trans_y_2d   ! k-summed ndiff merid tracer transport

  real    :: C_p ! heat capacity of seawater (J kg-1 K-1)

  type(remapping_CS) :: remap_CS
end type neutral_diffusion_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40)  :: mdl = "MOM_neutral_diffusion" ! module name

logical :: debug_this_module = .false. ! If true, verbose output of find neutral position

contains

!> Read parameters and allocate control structure for neutral_diffusion module.
logical function neutral_diffusion_init(Time, G, param_file, diag, CS)
  type(time_type), target,    intent(in)    :: Time       !< Time structure
  type(ocean_grid_type),      intent(in)    :: G          !< Grid structure
  type(diag_ctrl), target,    intent(inout) :: diag       !< Diagnostics control structure
  type(param_file_type),      intent(in)    :: param_file !< Parameter file structure
  type(neutral_diffusion_CS), pointer       :: CS         !< Neutral diffusion control structure

  ! Local variables
  character(len=256) :: mesg    ! Message for error messages.
  character(len=80)  :: string  ! Temporary strings

  if (associated(CS)) then
    call MOM_error(FATAL, "neutral_diffusion_init called with associated control structure.")
    return
  endif

  ! Log this module and master switch for turning it on/off
  call log_version(param_file, mdl, version, &
       "This module implements neutral diffusion of tracers")
  call get_param(param_file, mdl, "USE_NEUTRAL_DIFFUSION", neutral_diffusion_init, &
                 "If true, enables the neutral diffusion module.", &
                 default=.false.)

  if (.not.neutral_diffusion_init) then
    return
  endif

  allocate(CS)
  CS%diag => diag
 ! call openParameterBlock(param_file,'NEUTRAL_DIFF')

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "NDIFF_CONTINUOUS", CS%continuous_reconstruction, &
                 "If true, uses a continuous reconstruction of T and S when  \n"//  &
                 "finding neutral surfaces along which diffusion will happen.\n"//  &
                 "If false, a PPM discontinuous reconstruction of T and S    \n"//  &
                 "is done which results in a higher order routine but exacts \n"//  &
                 "a higher computational cost.", default=.true.)
  call get_param(param_file, mdl, "NDIFF_REF_PRES", CS%ref_pres,                    &
                 "The reference pressure (Pa) used for the derivatives of    \n"//  &
                 "the equation of state. If negative (default), local        \n"//  &
                 "pressure is used.", &
                 default = -1.)
  ! Initialize and configure remapping
  if (CS%continuous_reconstruction .eqv. .false.) then
    call get_param(param_file, mdl, "NDIFF_BOUNDARY_EXTRAP", CS%boundary_extrap, &
                   "Uses a rootfinding approach to find the position of a\n"//   &
                   "neutral surface within a layer taking into account the\n"//  &
                   "nonlinearity of the equation of state and the\n"//           &
                   "polynomial reconstructions of T/S.",                         &
                   default=.false.)
    call get_param(param_file, mdl, "NDIFF_REMAPPING_SCHEME", string, &
                   "This sets the reconstruction scheme used\n"//&
                   "for vertical remapping for all variables.\n"//&
                   "It can be one of the following schemes:\n"//&
                   trim(remappingSchemesDoc), default=remappingDefaultScheme)
    call initialize_remapping( CS%remap_CS, string, boundary_extrapolation = CS%boundary_extrap )
    call extract_member_remapping_CS(CS%remap_CS, degree=CS%ppoly_deg)
    call get_param(param_file, mdl, "NDIFF_REFINE_POSITION", CS%refine_position, &
                   "Uses a rootfinding approach to find the position of a\n"//   &
                   "neutral surface within a layer taking into account the\n"//  &
                   "nonlinearity of the equation of state and the\n"//           &
                   "polynomial reconstructions of T/S.",                         &
                   default=.false.)
    if (CS%refine_position) then
      call get_param(param_file, mdl, "NDIFF_TOLERANCE", CS%tolerance,            &
                     "Sets the convergence criterion for finding the neutral\n"// &
                     "position within a layer in kg m-3.",                        &
                     default=1.e-10)
      call get_param(param_file, mdl, "NDIFF_MAX_ITER", CS%max_iter,              &
                    "The maximum number of iterations to be done before \n"//     &
                     "exiting the iterative loop to find the neutral surface",    &
                     default=10)
    endif
    call get_param(param_file, mdl, "NDIFF_DEBUG", debug_this_module,             &
                   "Turns on verbose output for discontinuous neutral \n"//      &
                   "diffusion routines.", &
                   default = .false.)
  endif

! call get_param(param_file, mdl, "KHTR", CS%KhTr, &
!                "The background along-isopycnal tracer diffusivity.", &
!                units="m2 s-1", default=0.0)
!  call closeParameterBlock(param_file)
  if (CS%continuous_reconstruction) then
    CS%nsurf = 2*G%ke+2 ! Continuous reconstruction means that every interface has two connections
    allocate(CS%dRdT(SZI_(G),SZJ_(G),SZK_(G)+1)) ; CS%dRdT(:,:,:) = 0.
    allocate(CS%dRdS(SZI_(G),SZJ_(G),SZK_(G)+1)) ; CS%dRdS(:,:,:) = 0.
  else
    CS%nsurf = 4*G%ke   ! Discontinuous means that every interface has four connections
    allocate(CS%T_i(SZI_(G),SZJ_(G),SZK_(G),2))    ; CS%T_i(:,:,:,:) = 0.
    allocate(CS%S_i(SZI_(G),SZJ_(G),SZK_(G),2))    ; CS%S_i(:,:,:,:) = 0.
    allocate(CS%dRdT_i(SZI_(G),SZJ_(G),SZK_(G),2)) ; CS%dRdT_i(:,:,:,:) = 0.
    allocate(CS%dRdS_i(SZI_(G),SZJ_(G),SZK_(G),2)) ; CS%dRdS_i(:,:,:,:) = 0.
    allocate(CS%ppoly_coeffs_T(SZI_(G),SZJ_(G),SZK_(G),CS%ppoly_deg+1)) ; CS%ppoly_coeffs_T(:,:,:,:) = 0.
    allocate(CS%ppoly_coeffs_S(SZI_(G),SZJ_(G),SZK_(G),CS%ppoly_deg+1)) ; CS%ppoly_coeffs_S(:,:,:,:) = 0.
  endif
  ! T-points
  allocate(CS%Tint(SZI_(G),SZJ_(G),SZK_(G)+1)) ; CS%Tint(:,:,:) = 0.
  allocate(CS%Sint(SZI_(G),SZJ_(G),SZK_(G)+1)) ; CS%Sint(:,:,:) = 0.
  allocate(CS%Pint(SZI_(G),SZJ_(G),SZK_(G)+1)) ; CS%Pint(:,:,:) = 0.
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

!> Diagnostic handles for neutral diffusion tendencies.
subroutine neutral_diffusion_diag_init(Time, G, diag, C_p, Reg, CS)
  type(time_type),target,     intent(in)  :: Time   !< Time structure
  type(ocean_grid_type),      intent(in)  :: G      !< Grid structure
  type(diag_ctrl),            intent(in)  :: diag   !< Diagnostics control structure
  type(tracer_registry_type), intent(in)  :: Reg    !< Tracer structure
  real,                       intent(in)  :: C_p    !< Seawater heat capacity
  type(neutral_diffusion_CS), pointer     :: CS     !< Neutral diffusion control structure

  ! local
  integer :: n,ntr

  if(.not. associated(CS)) return

  ntr    = Reg%ntr
  CS%C_p = C_p

  allocate(CS%id_neutral_diff_tracer_conc_tend(ntr))
  allocate(CS%id_neutral_diff_tracer_cont_tend(ntr))
  allocate(CS%id_neutral_diff_tracer_cont_tend_2d(ntr))
  allocate(CS%id_neutral_diff_tracer_trans_x_2d(ntr))
  allocate(CS%id_neutral_diff_tracer_trans_y_2d(ntr))
  CS%id_neutral_diff_tracer_conc_tend(:)    = -1
  CS%id_neutral_diff_tracer_cont_tend(:)    = -1
  CS%id_neutral_diff_tracer_cont_tend_2d(:) = -1
  CS%id_neutral_diff_tracer_trans_x_2d(:)   = -1
  CS%id_neutral_diff_tracer_trans_y_2d(:)   = -1

  do n=1,ntr

    if(trim(Reg%Tr(n)%name) == 'T') then

      CS%id_neutral_diff_tracer_conc_tend(n) = register_diag_field('ocean_model',  &
      'ndiff_tracer_conc_tendency_'//trim(Reg%Tr(n)%name), diag%axesTL, Time,      &
      'Neutral diffusion tracer concentration tendency for '//trim(Reg%Tr(n)%name),&
      'degC s-1')

      CS%id_neutral_diff_tracer_cont_tend(n) = register_diag_field('ocean_model',                                      &
      'ndiff_tracer_cont_tendency_'//trim(Reg%Tr(n)%name), diag%axesTL, Time,                                          &
      'Neutral diffusion tracer content tendency for '//trim(Reg%Tr(n)%name),                                          &
      'W m-2',cmor_field_name='opottemppmdiff',                                                                     &
      cmor_standard_name=                                                                                              &
      'tendency_of_sea_water_potential_temperature_expressed_as_heat_content_due_to_parameterized_mesocale_diffusion', &
      cmor_long_name =                                                                                                 &
      'Tendency of sea water potential temperature expressed as heat content due to parameterized mesocale diffusion', &
      v_extensive=.true.)

      CS%id_neutral_diff_tracer_cont_tend_2d(n) = register_diag_field('ocean_model',                                                   &
      'ndiff_tracer_cont_tendency_2d_'//trim(Reg%Tr(n)%name), diag%axesT1, Time,                                                       &
      'Depth integrated neutral diffusion tracer content tendency for '//trim(Reg%Tr(n)%name),                                         &
      'W m-2',cmor_field_name='opottemppmdiff_2d',                                                                                  &
      cmor_standard_name=                                                                                                              &
      'tendency_of_sea_water_potential_temperature_expressed_as_heat_content_due_to_parameterized_mesocale_diffusion_depth_integrated',&
      cmor_long_name =                                                                                                                 &
      'Tendency of sea water potential temperature expressed as heat content due to parameterized mesocale diffusion depth integrated')

      CS%id_neutral_diff_tracer_trans_x_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_x_2d_'//trim(Reg%Tr(n)%name), diag%axesCu1, Time,                  &
      'Depth integrated neutral diffusion zonal tracer transport for '//trim(Reg%Tr(n)%name),&
      'W')

      CS%id_neutral_diff_tracer_trans_y_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_y_2d_'//trim(Reg%Tr(n)%name), diag%axesCv1, Time,                  &
      'Depth integrated neutral diffusion merid tracer transport for '//trim(Reg%Tr(n)%name),&
      'W')

    elseif(trim(Reg%Tr(n)%name) == 'S') then

      CS%id_neutral_diff_tracer_conc_tend(n) = register_diag_field('ocean_model',  &
      'ndiff_tracer_conc_tendency_'//trim(Reg%Tr(n)%name), diag%axesTL, Time,      &
      'Neutral diffusion tracer concentration tendency for '//trim(Reg%Tr(n)%name),&
      'tracer concentration * s-1')

      CS%id_neutral_diff_tracer_cont_tend(n) = register_diag_field('ocean_model',                         &
      'ndiff_tracer_cont_tendency_'//trim(Reg%Tr(n)%name), diag%axesTL, Time,                             &
      'Neutral diffusion tracer content tendency for '//trim(Reg%Tr(n)%name),                             &
      'kg m-2 s-1',cmor_field_name='osaltpmdiff',                                                         &
      cmor_standard_name=                                                                                 &
      'tendency_of_sea_water_salinity_expressed_as_salt_content_due_to_parameterized_mesocale_diffusion', &
      cmor_long_name =                                                                                    &
      'Tendency of sea water salinity expressed as salt content due to parameterized mesocale diffusion', &
      v_extensive=.true.)

      CS%id_neutral_diff_tracer_cont_tend_2d(n) = register_diag_field('ocean_model',                                      &
      'ndiff_tracer_cont_tendency_2d_'//trim(Reg%Tr(n)%name), diag%axesT1, Time,                                          &
      'Depth integrated neutral diffusion tracer content tendency for '//trim(Reg%Tr(n)%name),                            &
      'kg m-2 s-1',cmor_field_name='osaltpmdiff_2d',                                                                      &
      cmor_standard_name=                                                                                                 &
      'tendency_of_sea_water_salinity_expressed_as_salt_content_due_to_parameterized_mesocale_diffusion_depth_integrated',&
      cmor_long_name =                                                                                                    &
      'Tendency of sea water salinity expressed as salt content due to parameterized mesocale diffusion depth integrated')

      CS%id_neutral_diff_tracer_trans_x_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_x_2d_'//trim(Reg%Tr(n)%name), diag%axesCu1, Time,                  &
      'Depth integrated neutral diffusion zonal tracer transport for '//trim(Reg%Tr(n)%name),&
      'kg s-1')

      CS%id_neutral_diff_tracer_trans_y_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_y_2d_'//trim(Reg%Tr(n)%name), diag%axesCv1, Time,                  &
      'Depth integrated neutral diffusion merid tracer transport for '//trim(Reg%Tr(n)%name),&
      'kg s-1')

    else

      CS%id_neutral_diff_tracer_conc_tend(n) = register_diag_field('ocean_model',  &
      'ndiff_tracer_conc_tendency_'//trim(Reg%Tr(n)%name), diag%axesTL, Time,      &
      'Neutral diffusion tracer concentration tendency for '//trim(Reg%Tr(n)%name),&
       'tracer concentration * m-2 s-1')

      CS%id_neutral_diff_tracer_cont_tend(n) = register_diag_field('ocean_model',&
      'ndiff_tracer_cont_tendency_'//trim(Reg%Tr(n)%name), diag%axesTL, Time,    &
      'Neutral diffusion tracer content tendency for '//trim(Reg%Tr(n)%name),    &
      'tracer content * m-2 s-1', v_extensive=.true.)

      CS%id_neutral_diff_tracer_cont_tend_2d(n) = register_diag_field('ocean_model',          &
      'ndiff_tracer_cont_tendency_2d_'//trim(Reg%Tr(n)%name), diag%axesTL, Time,              &
      'Depth integrated neutral diffusion tracer content tendency for '//trim(Reg%Tr(n)%name),&
      'tracer content * m-2 s-1')

      CS%id_neutral_diff_tracer_trans_x_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_x_2d_'//trim(Reg%Tr(n)%name), diag%axesCu1, Time,                  &
      'Depth integrated neutral diffusion zonal tracer transport for '//trim(Reg%Tr(n)%name),&
      'kg s-1')

      CS%id_neutral_diff_tracer_trans_y_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_y_2d_'//trim(Reg%Tr(n)%name), diag%axesCv1, Time,                  &
      'Depth integrated neutral diffusion merid tracer transport for '//trim(Reg%Tr(n)%name),&
      'kg s-1')

    endif

  enddo

end subroutine neutral_diffusion_diag_init


!> Calculate remapping factors for u/v columns used to map adjoining columns to
!! a shared coordinate space.
subroutine neutral_diffusion_calc_coeffs(G, GV, h, T, S, EOS, CS)
  type(ocean_grid_type),                    intent(in) :: G   !< Ocean grid structure
  type(verticalGrid_type),                  intent(in) :: GV  !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h   !< Layer thickness (H units)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: T   !< Potential temperature (degC)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: S   !< Salinity (ppt)
  type(EOS_type),                           pointer    :: EOS !< Equation of state structure
  type(neutral_diffusion_CS),               pointer    :: CS  !< Neutral diffusion control structure

  ! Local variables
  integer :: i, j, k
  ! Variables used for reconstructions
  real, dimension(SZK_(G),2) :: ppoly_r_S            ! Reconstruction slopes
  integer :: iMethod
  real, dimension(SZI_(G)) :: ref_pres ! Reference pressure used to calculate alpha/beta

  ! If doing along isopycnal diffusion (as opposed to neutral diffusion, set the reference pressure)
  if (CS%ref_pres>=0.) ref_pres(:) = CS%ref_pres

  if (CS%continuous_reconstruction) then
    CS%dRdT(:,:,:) = 0.
    CS%dRdS(:,:,:) = 0.
  else
    CS%T_i(:,:,:,:) = 0.
    CS%S_i(:,:,:,:) = 0.
    CS%dRdT_i(:,:,:,:) = 0.
    CS%dRdS_i(:,:,:,:) = 0.
  endif

  ! Calculate pressure at interfaces
  CS%Pint(:,:,1) = 0.
  do k=1,G%ke ; do j=G%jsc-1, G%jec+1 ; do i=G%isc-1,G%iec+1
    CS%Pint(i,j,k+1) = CS%Pint(i,j,k) + h(i,j,k)*GV%H_to_Pa
  enddo ; enddo ; enddo

  do j = G%jsc-1, G%jec+1
    ! Interpolate state to interface
    do i = G%isc-1, G%iec+1
      if (CS%continuous_reconstruction) then
        call interface_scalar(G%ke, h(i,j,:), T(i,j,:), CS%Tint(i,j,:), 2)
        call interface_scalar(G%ke, h(i,j,:), S(i,j,:), CS%Sint(i,j,:), 2)
      else
        call build_reconstructions_1d( CS%remap_CS, G%ke, h(i,j,:), T(i,j,:), CS%ppoly_coeffs_T(i,j,:,:), &
                                       CS%T_i(i,j,:,:), ppoly_r_S, iMethod )
        call build_reconstructions_1d( CS%remap_CS, G%ke, h(i,j,:), S(i,j,:), CS%ppoly_coeffs_S(i,j,:,:), &
                                       CS%S_i(i,j,:,:), ppoly_r_S, iMethod )
      endif
    enddo

    ! Continuous reconstruction
    if (CS%continuous_reconstruction) then
      do k = 1, G%ke+1
        if (CS%ref_pres<0) ref_pres(:) = CS%Pint(:,j,k)
        call calculate_density_derivs(CS%Tint(:,j,k), CS%Sint(:,j,k), ref_pres, &
                                      CS%dRdT(:,j,k), CS%dRdS(:,j,k), G%isc-1, G%iec-G%isc+3, EOS)
      enddo
    else ! Discontinuous reconstruction
      do k = 1, G%ke
        if (CS%ref_pres<0) ref_pres(:) = CS%Pint(:,j,k)
        call calculate_density_derivs(CS%T_i(:,j,k,1), CS%S_i(:,j,k,1), ref_pres, &
                                      CS%dRdT_i(:,j,k,1), CS%dRdS_i(:,j,k,1), G%isc-1, G%iec-G%isc+3, EOS)
        if (CS%ref_pres<0) ref_pres(:) = CS%Pint(:,j,k+1)
        call calculate_density_derivs(CS%T_i(:,j,k,2), CS%S_i(:,j,k,2), ref_pres, &
                                         CS%dRdT_i(:,j,k,2), CS%dRdS_i(:,j,k,2), G%isc-1, G%iec-G%isc+3, EOS)
      enddo
    endif
  enddo

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
        call find_neutral_surface_positions_continuous(G%ke,                                    &
                CS%Pint(i,j,:), CS%Tint(i,j,:), CS%Sint(i,j,:), CS%dRdT(i,j,:), CS%dRdS(i,j,:),            &
                CS%Pint(i+1,j,:), CS%Tint(i+1,j,:), CS%Sint(i+1,j,:), CS%dRdT(i+1,j,:), CS%dRdS(i+1,j,:),  &
                CS%uPoL(I,j,:), CS%uPoR(I,j,:), CS%uKoL(I,j,:), CS%uKoR(I,j,:), CS%uhEff(I,j,:) )
      else
        call find_neutral_surface_positions_discontinuous(G%ke, CS%ppoly_deg,                                   &
            CS%Pint(i,j,:), CS%T_i(i,j,:,:), CS%S_i(i,j,:,:), CS%dRdT_i(i,j,:,:), CS%dRdS_i(i,j,:,:),           &
            CS%Pint(i+1,j,:), CS%T_i(i+1,j,:,:), CS%S_i(i+1,j,:,:), CS%dRdT_i(i+1,j,:,:), CS%dRdS_i(i+1,j,:,:), &
            CS%uPoL(I,j,:), CS%uPoR(I,j,:), CS%uKoL(I,j,:), CS%uKoR(I,j,:), CS%uhEff(I,j,:),                    &
            CS%refine_position, CS%ppoly_coeffs_T(i,j,:,:), CS%ppoly_coeffs_S(i,j,:,:),                         &
            CS%ppoly_coeffs_T(i+1,j,:,:), CS%ppoly_coeffs_S(i+1,j,:,:), EOS, CS%max_iter, CS%tolerance, CS%ref_pres)
      endif
    endif
  enddo ; enddo

  ! Neutral surface factors at V points
  do J = G%jsc-1, G%jec ; do i = G%isc, G%iec
    if (G%mask2dCv(i,J) > 0.) then
      if (CS%continuous_reconstruction) then
        call find_neutral_surface_positions_continuous(G%ke,                                  &
                CS%Pint(i,j,:), CS%Tint(i,j,:), CS%Sint(i,j,:), CS%dRdT(i,j,:), CS%dRdS(i,j,:),           &
                CS%Pint(i,j+1,:), CS%Tint(i,j+1,:), CS%Sint(i,j+1,:), CS%dRdT(i,j+1,:), CS%dRdS(i,j+1,:), &
                CS%vPoL(i,J,:), CS%vPoR(i,J,:), CS%vKoL(i,J,:), CS%vKoR(i,J,:), CS%vhEff(i,J,:) )
      else
        call find_neutral_surface_positions_discontinuous(G%ke, CS%ppoly_deg,                                   &
            CS%Pint(i,j,:), CS%T_i(i,j,:,:), CS%S_i(i,j,:,:), CS%dRdT_i(i,j,:,:), CS%dRdS_i(i,j,:,:),           &
            CS%Pint(i,j+1,:), CS%T_i(i,j+1,:,:), CS%S_i(i,j+1,:,:), CS%dRdT_i(i,j+1,:,:), CS%dRdS_i(i,j+1,:,:), &
            CS%vPoL(I,j,:), CS%vPoR(I,j,:), CS%vKoL(I,j,:), CS%vKoR(I,j,:), CS%vhEff(I,j,:),                    &
            CS%refine_position, CS%ppoly_coeffs_T(i,j,:,:), CS%ppoly_coeffs_S(i,j,:,:),                         &
            CS%ppoly_coeffs_T(i,j+1,:,:), CS%ppoly_coeffs_S(i,j+1,:,:), EOS, CS%max_iter, CS%tolerance, CS%ref_pres)
      endif
    endif
  enddo ; enddo

  CS%uhEff(:,:,:) = CS%uhEff(:,:,:) / GV%H_to_pa
  CS%vhEff(:,:,:) = CS%vhEff(:,:,:) / GV%H_to_pa

end subroutine neutral_diffusion_calc_coeffs

!> Update tracer concentration due to neutral diffusion; layer thickness unchanged by this update.
subroutine neutral_diffusion(G, GV, h, Coef_x, Coef_y, Tracer, m, dt, name, CS)
  type(ocean_grid_type),                     intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h      !< Layer thickness (H units)
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: Coef_x !< dt * Kh * dy / dx at u-points (m^2)
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: Coef_y !< dt * Kh * dx / dy at u-points (m^2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: Tracer !< Tracer concentration
  integer,                                   intent(in)    :: m      !< Tracer number
  real,                                      intent(in)    :: dt     !< Tracer time step * I_numitts (I_numitts in tracer_hordiff)
  character(len=32),                         intent(in)    :: name   !< Tracer name
  type(neutral_diffusion_CS),                pointer       :: CS     !< Neutral diffusion control structure

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),CS%nsurf-1) :: uFlx        ! Zonal flux of tracer      (concentration * H)
  real, dimension(SZI_(G),SZJB_(G),CS%nsurf-1) :: vFlx        ! Meridional flux of tracer (concentration * H)
  real, dimension(SZI_(G),SZJ_(G),G%ke)        :: tendency    ! tendency array for diagn
  real, dimension(SZI_(G),SZJ_(G))             :: tendency_2d ! depth integrated content tendency for diagn
  real, dimension(SZIB_(G),SZJ_(G))            :: trans_x_2d  ! depth integrated diffusive tracer x-transport diagn
  real, dimension(SZI_(G),SZJB_(G))            :: trans_y_2d  ! depth integrated diffusive tracer y-transport diagn
  real, dimension(G%ke)                        :: dTracer     ! change in tracer concentration due to ndiffusion
  integer :: i, j, k, ks, nk
  real :: ppt2mks, Idt, convert

  nk = GV%ke

  ! for diagnostics
  if(CS%id_neutral_diff_tracer_conc_tend(m)    > 0  .or.  &
     CS%id_neutral_diff_tracer_cont_tend(m)    > 0  .or.  &
     CS%id_neutral_diff_tracer_cont_tend_2d(m) > 0  .or.  &
     CS%id_neutral_diff_tracer_trans_x_2d(m)   > 0  .or.  &
     CS%id_neutral_diff_tracer_trans_y_2d(m)   > 0) then
     ppt2mks          = 0.001
     Idt              = 1.0/dt
     tendency(:,:,:)  = 0.0
     tendency_2d(:,:) = 0.0
     trans_x_2d(:,:)  = 0.0
     trans_y_2d(:,:)  = 0.0
     convert          = 1.0
     if(trim(name) == 'T') convert = CS%C_p  * GV%H_to_kg_m2
     if(trim(name) == 'S') convert = ppt2mks * GV%H_to_kg_m2
  endif

  uFlx(:,:,:) = 0.
  vFlx(:,:,:) = 0.

  ! x-flux
  do j = G%jsc,G%jec ; do I = G%isc-1,G%iec
    if (G%mask2dCu(I,j)>0.) then
      call neutral_surface_flux(nk, CS%nsurf, CS%ppoly_deg, h(i,j,:), h(i+1,j,:),       &
                                Tracer(i,j,:), Tracer(i+1,j,:), &
                                CS%uPoL(I,j,:), CS%uPoR(I,j,:), &
                                CS%uKoL(I,j,:), CS%uKoR(I,j,:), &
                                CS%uhEff(I,j,:), uFlx(I,j,:), &
                                CS%continuous_reconstruction, CS%remap_CS)
    endif
  enddo ; enddo

  ! y-flux
  do J = G%jsc-1,G%jec ; do i = G%isc,G%iec
    if (G%mask2dCv(i,J)>0.) then
      call neutral_surface_flux(nk, CS%nsurf, CS%ppoly_deg, h(i,j,:), h(i,j+1,:),       &
                                Tracer(i,j,:), Tracer(i,j+1,:), &
                                CS%vPoL(i,J,:), CS%vPoR(i,J,:), &
                                CS%vKoL(i,J,:), CS%vKoR(i,J,:), &
                                CS%vhEff(i,J,:), vFlx(i,J,:),   &
                                CS%continuous_reconstruction, CS%remap_CS)
    endif
  enddo ; enddo

  ! Update the tracer concentration from divergence of neutral diffusive flux components
  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    if (G%mask2dT(i,j)>0.) then

      dTracer(:) = 0.
      do ks = 1,CS%nsurf-1 ;
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
        Tracer(i,j,k) = Tracer(i,j,k) + dTracer(k) * &
                        ( G%IareaT(i,j) / ( h(i,j,k) + GV%H_subroundoff ) )
      enddo

      if(CS%id_neutral_diff_tracer_conc_tend(m)    > 0  .or.  &
         CS%id_neutral_diff_tracer_cont_tend(m)    > 0  .or.  &
         CS%id_neutral_diff_tracer_cont_tend_2d(m) > 0 ) then
        do k = 1, GV%ke
          tendency(i,j,k) = dTracer(k) * G%IareaT(i,j) * Idt
        enddo
      endif

    endif
  enddo ; enddo


  ! Diagnose vertically summed zonal flux, giving zonal tracer transport from ndiff.
  ! Note sign corresponds to downgradient flux convention.
  if(CS%id_neutral_diff_tracer_trans_x_2d(m) > 0) then
    do j = G%jsc,G%jec ; do I = G%isc-1,G%iec
      trans_x_2d(I,j) = 0.
      if (G%mask2dCu(I,j)>0.) then
        do ks = 1,CS%nsurf-1 ;
          trans_x_2d(I,j) = trans_x_2d(I,j) - Coef_x(I,j) * uFlx(I,j,ks)
        enddo
        trans_x_2d(I,j) = trans_x_2d(I,j) * Idt * convert
      endif
    enddo ; enddo
    call post_data(CS%id_neutral_diff_tracer_trans_x_2d(m), trans_x_2d(:,:), CS%diag)
  endif

  ! Diagnose vertically summed merid flux, giving meridional tracer transport from ndiff.
  ! Note sign corresponds to downgradient flux convention.
  if(CS%id_neutral_diff_tracer_trans_y_2d(m) > 0) then
    do J = G%jsc-1,G%jec ; do i = G%isc,G%iec
      trans_y_2d(i,J) = 0.
      if (G%mask2dCv(i,J)>0.) then
        do ks = 1,CS%nsurf-1 ;
          trans_y_2d(i,J) = trans_y_2d(i,J) - Coef_y(i,J) * vFlx(i,J,ks)
        enddo
        trans_y_2d(i,J) = trans_y_2d(i,J) * Idt * convert
      endif
    enddo ; enddo
    call post_data(CS%id_neutral_diff_tracer_trans_y_2d(m), trans_y_2d(:,:), CS%diag)
  endif

  ! post tendency of tracer content
  if(CS%id_neutral_diff_tracer_cont_tend(m) > 0) then
    call post_data(CS%id_neutral_diff_tracer_cont_tend(m), tendency(:,:,:)*convert, CS%diag)
  endif

  ! post depth summed tendency for tracer content
  if(CS%id_neutral_diff_tracer_cont_tend_2d(m) > 0) then
    do j = G%jsc,G%jec ; do i = G%isc,G%iec
      do k = 1, GV%ke
        tendency_2d(i,j) = tendency_2d(i,j) + tendency(i,j,k)
      enddo
    enddo ; enddo
    call post_data(CS%id_neutral_diff_tracer_cont_tend_2d(m), tendency_2d(:,:)*convert, CS%diag)
  endif

  ! post tendency of tracer concentration; this step must be
  ! done after posting tracer content tendency, since we alter
  ! the tendency array.
  if(CS%id_neutral_diff_tracer_conc_tend(m) > 0) then
    do k = 1, GV%ke ; do j = G%jsc,G%jec ; do i = G%isc,G%iec
      tendency(i,j,k) =  tendency(i,j,k) / ( h(i,j,k) + GV%H_subroundoff )
    enddo ; enddo ; enddo
    call post_data(CS%id_neutral_diff_tracer_conc_tend(m), tendency, CS%diag)
  endif


end subroutine neutral_diffusion

!> Returns interface scalar, Si, for a column of layer values, S.
subroutine interface_scalar(nk, h, S, Si, i_method)
  integer,               intent(in)    :: nk       !< Number of levels
  real, dimension(nk),   intent(in)    :: h        !< Layer thickness (H units)
  real, dimension(nk),   intent(in)    :: S        !< Layer scalar (conc, e.g. ppt)
  real, dimension(nk+1), intent(inout) :: Si       !< Interface scalar (conc, e.g. ppt)
  integer,               intent(in)    :: i_method !< =1 use average of PLM edges
                                                   !! =2 use continuous PPM edge interpolation
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
      Si(k) = ppm_edge(h(km2), h(k-1), h(k), h(kp1),  S(k-1), S(k), diff(k-1), diff(k))
    enddo
  endif
  Si(nk+1) = S(nk) + 0.5 * diff(nk)

end subroutine interface_scalar

!> Returns the PPM quasi-fourth order edge value at k+1/2 following
!! equation 1.6 in Colella & Woodward, 1984: JCP 54, 174-201.
real function ppm_edge(hkm1, hk, hkp1, hkp2,  Ak, Akp1, Pk, Pkp1)
  real, intent(in) :: hkm1 !< Width of cell k-1
  real, intent(in) :: hk   !< Width of cell k
  real, intent(in) :: hkp1 !< Width of cell k+1
  real, intent(in) :: hkp2 !< Width of cell k+2
  real, intent(in) :: Ak   !< Average scalar value of cell k
  real, intent(in) :: Akp1 !< Average scalar value of cell k+1
  real, intent(in) :: Pk   !< PLM slope for cell k
  real, intent(in) :: Pkp1 !< PLM slope for cell k+1

  ! Local variables
  real :: R_hk_hkp1, R_2hk_hkp1, R_hk_2hkp1, f1, f2, f3, f4
  real, parameter :: h_neglect = 1.e-30

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
  real, dimension(nk), intent(in)    :: h        !< Layer thickness (H units)
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
subroutine find_neutral_surface_positions_continuous(nk, Pl, Tl, Sl, dRdTl, dRdSl, Pr, Tr, Sr, dRdTr, dRdSr, PoL, &
                                                     PoR, KoL, KoR, hEff)
  integer,                    intent(in)    :: nk    !< Number of levels
  real, dimension(nk+1),      intent(in)    :: Pl    !< Left-column interface pressure (Pa)
  real, dimension(nk+1),      intent(in)    :: Tl    !< Left-column interface potential temperature (degC)
  real, dimension(nk+1),      intent(in)    :: Sl    !< Left-column interface salinity (ppt)
  real, dimension(nk+1),      intent(in)    :: dRdTl !< Left-column dRho/dT (kg/m3/degC)
  real, dimension(nk+1),      intent(in)    :: dRdSl !< Left-column dRho/dS (kg/m3/ppt)
  real, dimension(nk+1),      intent(in)    :: Pr    !< Right-column interface pressure (Pa)
  real, dimension(nk+1),      intent(in)    :: Tr    !< Right-column interface potential temperature (degC)
  real, dimension(nk+1),      intent(in)    :: Sr    !< Right-column interface salinity (ppt)
  real, dimension(nk+1),      intent(in)    :: dRdTr !< Left-column dRho/dT (kg/m3/degC)
  real, dimension(nk+1),      intent(in)    :: dRdSr !< Left-column dRho/dS (kg/m3/ppt)
  real, dimension(2*nk+2),    intent(inout) :: PoL   !< Fractional position of neutral surface within layer KoL of left column
  real, dimension(2*nk+2),    intent(inout) :: PoR   !< Fractional position of neutral surface within layer KoR of right column
  integer, dimension(2*nk+2), intent(inout) :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(2*nk+2), intent(inout) :: KoR   !< Index of first right interface above neutral surface
  real, dimension(2*nk+1),    intent(inout) :: hEff  !< Effective thickness between two neutral surfaces (Pa)

  ! Local variables
  integer :: ns                     ! Number of neutral surfaces
  integer :: k_surface              ! Index of neutral surface
  integer :: kl                     ! Index of left interface
  integer :: kr                     ! Index of right interface
  real    :: dRdT, dRdS             ! dRho/dT and dRho/dS for the neutral surface
  logical :: searching_left_column  ! True if searching for the position of a right interface in the left column
  logical :: searching_right_column ! True if searching for the position of a left interface in the right column
  logical :: reached_bottom         ! True if one of the bottom-most interfaces has been used as the target
  integer :: krm1, klm1
  real    :: dRho, dRhoTop, dRhoBot, hL, hR
  integer :: lastK_left, lastK_right
  real    :: lastP_left, lastP_right

  ns = 2*nk+2
  ! Initialize variables for the search
  kr = 1 ; lastK_right = 1 ; lastP_right = 0.
  kl = 1 ; lastK_left = 1 ; lastP_left = 0.
  reached_bottom = .false.

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
        ! between right and left is zero.
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
      dRhoTop = 0.5 * ( ( dRdTr(krm1) + dRdTl(kl) ) * ( Tr(krm1) - Tl(kl) ) &
                     + ( dRdSr(krm1) + dRdSl(kl) ) * ( Sr(krm1) - Sl(kl) ) )
      ! Potential density difference, rho(kr) - rho(kl) (will be positive)
      dRhoBot = 0.5 * ( ( dRdTr(krm1+1) + dRdTl(kl) ) * ( Tr(krm1+1) - Tl(kl) ) &
                   + ( dRdSr(krm1+1) + dRdSl(kl) ) * ( Sr(krm1+1) - Sl(kl) ) )

      ! Because we are looking right, the left surface, kl, is lighter than krm1+1 and should be denser than krm1
      ! unless we are still at the top of the right column (kr=1)
      if (dRhoTop >= 0. .or. kr+kl==2) then
        PoR(k_surface) = 0. ! The left surface is lighter than anything in layer krm1
      elseif (dRhoTop >= dRhoBot) then ! Right layer is unstratified
        PoR(k_surface) = 1.
      else
        ! Linearly interpolate for the position between Pr(kr-1) and Pr(kr) where the density difference
        ! between right and left is zero.
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

!> Higher order version of find_neutral_surface_positions. Returns positions within left/right columns
!! of combined interfaces using intracell reconstructions of T/S
subroutine find_neutral_surface_positions_discontinuous(nk, deg,                                                   &
                Pres_l, Tl, Sl, dRdT_l, dRdS_l, Pres_r, Tr, Sr, dRdT_r, dRdS_r, PoL, PoR, KoL, KoR, hEff,          &
                refine_pos_in, ppoly_T_l, ppoly_S_l, ppoly_T_r, ppoly_S_r, EOS, max_iter, tolerance, ref_pres)
  integer,                    intent(in)    :: nk        !< Number of levels
  integer,                    intent(in)    :: deg       !< Degree of polynomial used for reconstructions
  real, dimension(nk+1),      intent(in)    :: Pres_l    !< Left-column interface pressure (Pa)
  real, dimension(nk,2),      intent(in)    :: Tl        !< Left-column top interface potential temperature (degC)
  real, dimension(nk,2),      intent(in)    :: Sl        !< Left-column top interface salinity (ppt)
  real, dimension(nk,2),      intent(in)    :: dRdT_l    !< Left-column, top interface dRho/dT (kg/m3/degC)
  real, dimension(nk,2),      intent(in)    :: dRdS_l    !< Left-column, top interface dRho/dS (kg/m3/ppt)
  real, dimension(nk+1),      intent(in)    :: Pres_r    !< Right-column interface pressure (Pa)
  real, dimension(nk,2),      intent(in)    :: Tr        !< Right-column top interface potential temperature (degC)
  real, dimension(nk,2),      intent(in)    :: Sr        !< Right-column top interface salinity (ppt)
  real, dimension(nk,2),      intent(in)    :: dRdT_r    !< Right-column, top interface dRho/dT (kg/m3/degC)
  real, dimension(nk,2),      intent(in)    :: dRdS_r    !< Right-column, top interface dRho/dS (kg/m3/ppt)
  real, dimension(4*nk),      intent(inout) :: PoL       !< Fractional position of neutral surface within
                                                         !! layer KoL of left column
  real, dimension(4*nk),      intent(inout) :: PoR       !< Fractional position of neutral surface within
                                                         !! layer KoR of right column
  integer, dimension(4*nk),   intent(inout) :: KoL       !< Index of first left interface above neutral surface
  integer, dimension(4*nk),   intent(inout) :: KoR       !< Index of first right interface above neutral surface
  real, dimension(4*nk-1),    intent(inout) :: hEff      !< Effective thickness between two neutral surfaces (Pa)
  logical,                    optional, intent(in) :: refine_pos_in !< True if rootfinding is used for position
  real, dimension(nk,deg+1),  optional, intent(in) :: ppoly_T_l !< Left-column coefficients of T reconstruction
  real, dimension(nk,deg+1),  optional, intent(in) :: ppoly_S_l !< Left-column coefficients of S reconstruction
  real, dimension(nk,deg+1),  optional, intent(in) :: ppoly_T_r !< Right-column coefficients of T reconstruction
  real, dimension(nk,deg+1),  optional, intent(in) :: ppoly_S_r !< Right-column coefficients of S reconstruction
  type(EOS_type),             optional, pointer    :: EOS       !< Equation of state structure
  integer,                    optional, intent(in) :: max_iter  !< Maximum number of iterations in refine_position
  real,                       optional, intent(in) :: tolerance !< Convergence criterion for refine_position
  real,                       optional, intent(in) :: ref_pres  !< Reference pressure to use for deriviative calculation

  ! Local variables
  integer :: ns                     ! Number of neutral surfaces
  integer :: k_surface              ! Index of neutral surface
  integer :: kl_left, kl_right      ! Index of layers on the left/right
  integer :: ki_left, ki_right      ! Index of interfaces on the left/right
  logical :: searching_left_column  ! True if searching for the position of a right interface in the left column
  logical :: searching_right_column ! True if searching for the position of a left interface in the right column
  logical :: reached_bottom         ! True if one of the bottom-most interfaces has been used as the target
  logical :: refine_pos             ! Use rootfinding to find the true neutral surface position
  real    :: dRho, dRhoTop, dRhoBot, dRhoTopm1, hL, hR
  integer :: lastK_left, lastK_right, maxK_left, maxK_right
  real    :: lastP_left, lastP_right
  real    :: min_bound
  logical, dimension(nk) :: top_connected_l, top_connected_r
  logical, dimension(nk) :: bot_connected_l, bot_connected_r

  ns = 4*nk
  top_connected_l(:) = .false. ; top_connected_r(:) = .false.
  bot_connected_l(:) = .false. ; bot_connected_r(:) = .false.
  maxK_left = -1 ; maxK_right = -1
  ! Vectors with all the values of the discontinuous reconstruction.
  ! Dimensions are [number of layers x number of interfaces]. Second dimension = 1 for top interface, = 2 for bottom
!  real, dimension(nk,2) :: Sl, Sr, Tl, Tr, dRdT_l, dRdS_l, dRdT_r, dRdS_r

! Check to make sure that polynomial reconstructions were passed if refine_pos defined)
  refine_pos = .false.
  if (present(refine_pos_in)) then
    refine_pos = refine_pos_in
    if (refine_pos .and. (.not. ( present(ppoly_T_l) .and. present(ppoly_S_l) .and.  &
                                  present(ppoly_T_r) .and. present(ppoly_S_r) .and.  &
                                  present(tolerance) .and. present(max_iter) .and. present(ref_pres) ) )) &
        call MOM_error(FATAL, "fine_neutral_surface_positions_discontinuous: refine_pos is requested, but polynomial"// &
                              "coefficients not available for T and S")
  endif

  ! Initialize variables for the search
  kl_right = 1 ; ki_right = 1 ; lastK_right = 1 ; lastP_right = -1.
  kl_left = 1  ; ki_left = 1  ; lastK_left = 1  ; lastP_left = -1.

  reached_bottom = .false.
  searching_left_column = .false.
  searching_right_column = .false.

  ! Loop over each neutral surface, working from top to bottom
  neutral_surfaces: do k_surface = 1, 4*nk
    ! Potential density difference, rho(kr) - rho(kl)
    dRho = 0.5 * &
      ( ( dRdT_r(kl_right,ki_right) + dRdT_l(kl_left,ki_left) ) * ( Tr(kl_right,ki_right) - Tl(kl_left,ki_left) ) &
      + ( dRdS_r(kl_right,ki_right) + dRdS_l(kl_left,ki_left) ) * ( Sr(kl_right,ki_right) - Sl(kl_left,ki_left) ) )
    if (debug_this_module)  write(*,'(A,I2,A,E12.4,A,I2,A,I2,A,I2,A,I2)') "k_surface=",k_surface,"  dRho=",dRho,"  kl_left=",kl_left,   &
      "  ki_left=",ki_left,"  kl_right=",kl_right, "  ki_right=",ki_right
    ! Which column has the lighter surface for the current indexes, kr and kl
    if (.not. reached_bottom) then
      if (dRho < 0.) then
        searching_left_column = .true.
        searching_right_column = .false.
      elseif (dRho > 0.) then
        searching_right_column = .true.
        searching_left_column = .false.
      else ! dRho == 0.
        if ((kl_left + kl_left == 2) .and. (ki_left + ki_right == 2)) then ! Still at surface
          searching_left_column = .true.
          searching_right_column = .false.
        else ! Not the surface so we simply change direction
          searching_left_column = .not. searching_left_column
          searching_right_column = .not. searching_right_column
        endif
      endif
    endif

    if (searching_left_column) then
      ! Determine differences between right column interface and potentially three different parts of the left
      ! Potential density difference, rho(kl-1) - rho(kr) (should be negative)
      dRhoTop = 0.5 * &
        ( ( dRdT_l(kl_left,1) + dRdT_r(kl_right,ki_right) ) * ( Tl(kl_left,1) - Tr(kl_right,ki_right) ) &
        + ( dRdS_l(kl_left,1) + dRdS_r(kl_right,ki_right) ) * ( Sl(kl_left,1) - Sr(kl_right,ki_right) ) )
      ! Potential density difference, rho(kl) - rho(kl_right,ki_right) (will be positive)
      dRhoBot = 0.5 * &
        ( ( dRdT_l(kl_left,2) + dRdT_r(kl_right,ki_right) ) * ( Tl(kl_left,2) - Tr(kl_right,ki_right) ) &
        + ( dRdS_l(kl_left,2) + dRdS_r(kl_right,ki_right) ) * ( Sl(kl_left,2) - Sr(kl_right,ki_right) ) )
      if (kl_left>1) then ! Calculate the density difference at top of discontinuity
        dRhoTopm1 = 0.5 * &
          ( ( dRdT_l(kl_left-1,2) + dRdT_r(kl_right,ki_right) ) * ( Tl(kl_left-1,2) - Tr(kl_right,ki_right) ) &
          + ( dRdS_l(kl_left-1,2) + dRdS_r(kl_right,ki_right) ) * ( Sl(kl_left-1,2) - Sr(kl_right,ki_right) ) )
      else
        dRhoTopm1 = dRhoTop
      endif
      if (debug_this_module) then
        write(*,'(A,I2,A,E12.4,A,E12.4,A,E12.4)') "Searching left layer ", kl_left, ":  dRhoTopm1=", dRhoTopm1, &
                                                  "  dRhoTop=", dRhoTop, "  dRhoBot=", dRhoBot
        write(*,'(A,I2,X,I2)') "Searching from right: ", kl_right, ki_right
        write(*,*) "Temp/Salt Reference: ", Tr(kl_right,ki_right), Sr(kl_right,ki_right)
        write(*,*) "Temp/Salt Top L: ", Tl(kl_left,1), Sl(kl_left,1)
        write(*,*) "Temp/Salt Bot L: ", Tl(kl_left,2), Sl(kl_left,2)
      endif

      ! Set the position within the starting column
      PoR(k_surface) = REAL(ki_right-1)
      KoR(k_surface) = REAL(kl_right)

      ! Set position within the searched column
      call search_other_column_discontinuous(dRhoTopm1, dRhoTop, dRhoBot, Pres_l(kl_left), Pres_l(kl_left+1),          &
             lastP_left, lastK_left, kl_left, ki_left, top_connected_l, bot_connected_l, PoL(k_surface), KoL(k_surface))
      if ( refine_pos .and. (PoL(k_surface) > 0.) .and. (PoL(k_surface) < 1.) ) then
        min_bound = 0.
        if ( (k_surface > 1) .and. ( KoL(k_surface) == KoL(k_surface-1) ) ) min_bound = PoL(k_surface-1)
        PoL(k_surface) = refine_nondim_position(max_iter, tolerance, Tr(kl_right,ki_right), Sr(kl_right,ki_right),     &
                    dRdT_r(kl_right,ki_right), dRdS_r(kl_right,ki_right), Pres_l(kl_left), Pres_l(kl_left+1),          &
                    deg, ppoly_T_l(kl_left,:), ppoly_S_l(kl_left,:), EOS, PoL(k_surface), dRhoTop, dRhoBot, min_bound, &
                    ref_pres)
      endif
      if (PoL(k_surface) == 0.) top_connected_l(KoL(k_surface)) = .true.
      if (PoL(k_surface) == 1.) bot_connected_l(KoL(k_surface)) = .true.
      call increment_interface(nk, kl_right, ki_right, reached_bottom, searching_right_column, searching_left_column)
      lastK_left = KoL(k_surface)  ; lastP_left = PoL(k_surface)

    elseif (searching_right_column) then
      ! Interpolate for the neutral surface position within the right column, layer krm1
      ! Potential density difference, rho(kr-1) - rho(kl) (should be negative)
      dRhoTop = 0.5 * &
        ( ( dRdT_r(kl_right,1) + dRdT_l(kl_left,ki_left) ) * ( Tr(kl_right,1) - Tl(kl_left,ki_left) ) &
        + ( dRdS_r(kl_right,1) + dRdS_l(kl_left,ki_left) ) * ( Sr(kl_right,1) - Sl(kl_left,ki_left) ) )
      dRhoBot = 0.5 * &
        ( ( dRdT_r(kl_right,2) + dRdT_l(kl_left,ki_left) ) * ( Tr(kl_right,2) - Tl(kl_left,ki_left) ) &
        + ( dRdS_r(kl_right,2) + dRdS_l(kl_left,ki_left) ) * ( Sr(kl_right,2) - Sl(kl_left,ki_left) ) )
      if (kl_right>1) then
        dRhoTopm1 = 0.5 * &
          ( ( dRdT_r(kl_right-1,2) + dRdT_l(kl_left,ki_left) ) * ( Tr(kl_right-1,2) - Tl(kl_left,ki_left) ) &
          + ( dRdS_r(kl_right-1,2) + dRdS_l(kl_left,ki_left) ) * ( Sr(kl_right-1,2) - Sl(kl_left,ki_left) ) )
      else
        dRhoTopm1 = dRhoTop
      endif
      if (debug_this_module) then
        write(*,'(A,I2,A,E12.4,A,E12.4,A,E12.4)') "Searching right layer ", kl_right, ":  dRhoTopm1=", dRhoTopm1, &
                                                  "  dRhoTop=", dRhoTop, "  dRhoBot=", dRhoBot
        write(*,'(A,I2,X,I2)') "Searching from left: ", kl_left, ki_left
        write(*,*) "Temp/Salt Reference: ", Tl(kl_left,ki_left), Sl(kl_left,ki_left)
        write(*,*) "Temp/Salt Top R: ", Tr(kl_right,1), Sr(kl_right,1)
        write(*,*) "Temp/Salt Bot R: ", Tr(kl_right,2), Sr(kl_right,2)
      endif
      ! Set the position within the starting column
      PoL(k_surface) = REAL(ki_left-1)
      KoL(k_surface) = REAL(kl_left)

      ! Set position within the searched column
      call search_other_column_discontinuous(dRhoTopm1, dRhoTop, dRhoBot, Pres_r(kl_right), Pres_r(kl_right+1),        &
         lastP_right, lastK_right, kl_right, ki_right, top_connected_r, bot_connected_r, PoR(k_surface), KoR(k_surface))
      if ( refine_pos .and. (PoR(k_surface) > 0. .and. PoR(k_surface) < 1.) ) then
        min_bound = 0.
        if ( (k_surface > 1) .and. ( KoR(k_surface) == KoR(k_surface-1) ) ) min_bound = PoR(k_surface-1)
        PoR(k_surface) = refine_nondim_position(max_iter, tolerance, Tl(kl_left,ki_left), Sl(kl_left,ki_left),         &
                  dRdT_l(kl_left,ki_left), dRdS_l(kl_left,ki_left), Pres_r(kl_right), Pres_r(kl_right+1),              &
                  deg, ppoly_T_r(kl_right,:), ppoly_S_r(kl_right,:), EOS, PoR(k_surface), dRhoTop, dRhoBot, min_bound, &
                  ref_pres)
      endif
      if (PoR(k_surface) == 0.) top_connected_r(KoR(k_surface)) = .true.
      if (PoR(k_surface) == 1.) bot_connected_r(KoR(k_surface)) = .true.
      call increment_interface(nk, kl_left, ki_left, reached_bottom, searching_left_column, searching_right_column)
      lastK_right = KoR(k_surface) ; lastP_right = PoR(k_surface)

    else
      stop 'Else what?'
    endif
    if (debug_this_module)  write(*,'(A,I2,A,F6.2,A,I2,A,F6.2)') "KoL:", KoL(k_surface), " PoL:", PoL(k_surface), "     KoR:", &
      KoR(k_surface), " PoR:", PoR(k_surface)
    maxK_left= MAX(KoL(k_surface), maxK_left)
    maxK_right= MAX(KoR(k_surface), maxK_right)
    ! Effective thickness
    ! NOTE: This would be better expressed in terms of the layers thicknesses rather
    ! than as differences of position - AJA
    if (k_surface>1) then
      hL = absolute_position(nk,ns,Pres_l,KoL,PoL,k_surface) - absolute_position(nk,ns,Pres_l,KoL,PoL,k_surface-1)
      hR = absolute_position(nk,ns,Pres_r,KoR,PoR,k_surface) - absolute_position(nk,ns,Pres_r,KoR,PoR,k_surface-1)
      ! In the case of a layer being unstably stratified, may get a negative thickness. Set the previous position
      ! to the current location
      if (hL < 0.) then
        if ( (KoL(k_surface)<maxK_left) .and. PoL(k_surface)==1.) then
          KoL(k_surface) = maxK_left ; PoL(k_surface) = 0.
        endif
        PoL(k_surface-1) = PoL(k_surface) ; KoL(k_surface-1) = KoL(k_surface)
        ! Make sure that the right point across a continuity is chose
        hL = 0.
      endif
      if (hR < 0.) then
        if ( (KoR(k_surface)<maxK_right) .and. PoR(k_surface)==1.) then
          KoR(k_surface) = maxK_right ; PoR(k_surface) = 0.
        endif
        PoR(k_surface-1) = PoR(k_surface) ; KoR(k_surface-1) = KoR(k_surface)
        hR = 0.
      endif
      if ( hL + hR > 0.) then
        hEff(k_surface-1) = 2. * hL * hR / ( hL + hR ) ! Harmonic mean
      else
        hEff(k_surface-1) = 0.
      endif
    endif

  enddo neutral_surfaces

end subroutine find_neutral_surface_positions_discontinuous

!> Increments the interface which was just connected and also set flags if the bottom is reached
subroutine increment_interface(nk, kl, ki, reached_bottom, searching_this_column, searching_other_column)
  integer, intent(in   ) :: nk                     !< Number of vertical levels
  integer, intent(inout) :: kl                     !< Current layer (potentially updated)
  integer, intent(inout) :: ki                     !< Current interface
  logical, intent(inout) :: reached_bottom         !< Updated when kl == nk and ki == 2
  logical, intent(inout) :: searching_this_column  !< Updated when kl == nk and ki == 2
  logical, intent(inout) :: searching_other_column !< Updated when kl == nk and ki == 2

  if (ki == 1) then
    ki = 2
  elseif ((ki == 2) .and. (kl < nk)) then
    ki = 1
    kl = kl + 1
  elseif ((kl == nk) .and. (ki==2)) then
    reached_bottom = .true.
    searching_this_column = .true.
    searching_other_column = .false.
  else
    call MOM_error(FATAL,"Unanticipated eventuality in increment_interface")
  endif


end subroutine increment_interface

!> Searches the "other" (searched) column for the position of the neutral surface
subroutine search_other_column_discontinuous(dRhoTopm1, dRhoTop, dRhoBot, Ptop, Pbot, lastP, lastK, kl, ki, &
                                             top_connected, bot_connected, out_P, out_K)
  real,                  intent(in   ) :: dRhoTopm1      !< Density difference across previous interface
  real,                  intent(in   ) :: dRhoTop        !< Density difference across top interface
  real,                  intent(in   ) :: dRhoBot        !< Density difference across top interface
  real,                  intent(in   ) :: Ptop           !< Pressure at top interface
  real,                  intent(in   ) :: Pbot           !< Pressure at bottom interface
  real,                  intent(in   ) :: lastP          !< Last position connected in the searched column
  integer,               intent(in   ) :: lastK          !< Last layer connected in the searched column
  integer,               intent(in   ) :: kl             !< Layer in the searched column
  integer,               intent(in   ) :: ki             !< Interface of the searched column
  logical, dimension(:), intent(inout) :: top_connected  !< True if the top interface was pointed to
  logical, dimension(:), intent(inout) :: bot_connected  !< True if the top interface was pointed to
  real,                  intent(  out) :: out_P          !< Position within searched column
  integer,               intent(  out) :: out_K          !< Layer within searched column
  ! Local variables
  logical :: search_layer

  search_layer = .true.
  ! Bad values to make sure that the particular setup has been processed
  out_P = -1. ; out_K = -1
  ! Check if everything in this layer is denser than neutral surface or if at the top of the water column
  if ((kl==1 .and. ki==1)) then
    if (debug_this_module)  write(*,*) "At surface"
    out_P = 0. ; out_K = kl
    search_layer = .false.
  ! Deal with the case where reconstruction is continuous
  elseif ( kl>1 ) then
    if (dRhoTopm1==dRhoTop .and. dRhoTopm1 == 0. .and. (.not.bot_connected(kl-1)) ) then
      out_P = 1 ; out_K = kl-1;
      search_layer = .false.
    elseif ( (dRhoTopm1<dRhoTop .and. dRhoTop > 0.) ) then
      out_P = 1. ; out_K = kl-1
      search_layer = .false.
    endif
  endif

  if (search_layer) then
    if (dRhoTop > 0.) then
      out_P = 0. ; out_K = kl
    elseif ( dRhoTop == 0. .and. (.not. top_connected(kl)) ) then
      out_P = 0. ; out_K = kl
    elseif (dRhoTop >= dRhoBot) then
      out_P = 1. ; out_K = kl
    else
      if (debug_this_module)  write(*,*) "Zero crossing point within layer"
      out_P = interpolate_for_nondim_position(dRhoTop, Ptop, dRhoBot, Pbot)
      out_K = kl
    endif
  endif

  if ( (out_P < 0.) .and. (out_K < 0) ) then
    call MOM_error(WARNING, "Unanticipated case in search_other_column_discontinuous")
  endif
  ! Check to make sure that the layer index is always increasing
  if ( (out_K < lastK) .and. lastP==0. .and. out_P == 1. ) then
    out_K = lastK ; out_P = 0.
  endif

end subroutine search_other_column_discontinuous
!> Converts non-dimensional position within a layer to absolute position (for debugging)
real function absolute_position(n,ns,Pint,Karr,NParr,k_surface)
  integer, intent(in) :: n            !< Number of levels
  integer, intent(in) :: ns           !< Number of neutral surfaces
  real,    intent(in) :: Pint(n+1)    !< Position of interfaces (Pa)
  integer, intent(in) :: Karr(ns)     !< Index of interface above position
  real,    intent(in) :: NParr(ns)    !< Non-dimensional position within layer Karr(:)

  ! Local variables
  integer :: k_surface, k

  k = Karr(k_surface)
  if (k>n) stop 'absolute_position: k>nk is out of bounds!'
  absolute_position = Pint(k) + NParr(k_surface) * ( Pint(k+1) - Pint(k) )

end function absolute_position

!> Converts non-dimensional positions within layers to absolute positions (for debugging)
function absolute_positions(n,ns,Pint,Karr,NParr)
  integer, intent(in) :: n            !< Number of levels
  integer, intent(in) :: ns           !< Number of neutral surfaces
  real,    intent(in) :: Pint(n+1)    !< Position of interface (Pa)
  integer, intent(in) :: Karr(ns)  !< Indexes of interfaces about positions
  real,    intent(in) :: NParr(ns) !< Non-dimensional positions within layers Karr(:)

  real,  dimension(ns) :: absolute_positions ! Absolute positions (Pa)

  ! Local variables
  integer :: k_surface, k

  do k_surface = 1, ns
    absolute_positions(k_surface) = absolute_position(n,ns,Pint,Karr,NParr,k_surface)
  enddo

end function absolute_positions

!> Returns the non-dimensional position between Pneg and Ppos where the
!! interpolated density difference equals zero.
!! The result is always bounded to be between 0 and 1.
real function interpolate_for_nondim_position(dRhoNeg, Pneg, dRhoPos, Ppos)
  real, intent(in) :: dRhoNeg !< Negative density difference
  real, intent(in) :: Pneg    !< Position of negative density difference
  real, intent(in) :: dRhoPos !< Positive density difference
  real, intent(in) :: Ppos    !< Position of positive density difference

  if (Ppos<Pneg) stop 'interpolate_for_nondim_position: Houston, we have a problem! Ppos<Pneg'
  if (dRhoNeg>dRhoPos) write(0,*) 'dRhoNeg, Pneg, dRhoPos, Ppos=',dRhoNeg, Pneg, dRhoPos, Ppos
  if (dRhoNeg>dRhoPos) stop 'interpolate_for_nondim_position: Houston, we have a problem! dRhoNeg>dRhoPos'
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
  if ( interpolate_for_nondim_position < 0. ) stop 'interpolate_for_nondim_position: Houston, we have a problem! Pint < Pneg'
  if ( interpolate_for_nondim_position > 1. ) stop 'interpolate_for_nondim_position: Houston, we have a problem! Pint > Ppos'
end function interpolate_for_nondim_position

!> Use root-finding methods to find where dRho = 0, based on the equation of state and the polynomial
!! reconstructions of temperature, salinity. Initial guess is based on the zero crossing of based on linear
!! profiles of dRho, T, and S, between the top and bottom interface. If second derivatives of the EOS are available,
!! it starts with a Newton's method. However, Newton's method is not guaranteed to be bracketed, a check is performed
!! to see if it it diverges outside the interval. In that case (or in the case that second derivatives are not
!! available), Brent's method is used following the implementation found at
!! https://people.sc.fsu.edu/~jburkardt/f_src/brent/brent.f90
real function refine_nondim_position(max_iter, tolerance, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, deg, &
      ppoly_T, ppoly_S, EOS, x0, drho_top, drho_bot, min_bound, ref_pres, force_brent)
  integer,            intent(in) :: max_iter    !< Number of maximum iterations to use
  real,               intent(in) :: tolerance   !< Convergence criterion for delta_rho
  real,               intent(in) :: T_ref       !< Temperature of the neutral surface at the searched from interface
  real,               intent(in) :: S_ref       !< Salinity of the neutral surface at the searched from interface
  real,               intent(in) :: alpha_ref   !< dRho/dT of the neutral surface at the searched from interface
  real,               intent(in) :: beta_ref    !< dRho/dS of the neutral surface at the searched from interface
  real,               intent(in) :: P_top       !< Pressure at the top interface in the layer to be searched
  real,               intent(in) :: P_bot       !< Pressure at the bottom interface in the layer to be searched
  integer,            intent(in) :: deg         !< Order of the polynomimal used for reconstructions
  real, dimension(:), intent(in) :: ppoly_T     !< Coefficients of the order N polynomial reconstruction of T within
                                                !! the layer to be searched.
  real, dimension(:), intent(in) :: ppoly_S     !< Coefficients of the order N polynomial reconstruction of T within
                                                !! the layer to be searched.
  real,               intent(in) :: x0          !< Nondimensional position within the layer where the neutral
                                                !! surface connects. If interpolate_for_nondim_position was
                                                !! previously called, this would be based on linear profile of dRho
  real,               intent(in) :: drho_top, drho_bot, min_bound
  real,               intent(in) :: ref_pres    !< Optionally use a different reference pressure other than local
  type(EOS_type),     pointer    :: EOS         !< Equation of state structure
  logical, optional,  intent(in) :: force_brent !< Forces the use of Brent's method instead of Newton's method to find
                                                !! position of neutral surface

  ! Local variables
  integer :: form_of_EOS
  integer :: iter
  logical :: do_newton, do_brent

  real :: delta_rho, d_delta_rho_dP ! Terms for the Newton iteration
  real :: P_int, P_min, P_ref ! Interpolated pressure
  real :: delta_rho_init, delta_rho_final, x_init
  real :: T, S, alpha, beta, alpha_avg, beta_avg
  ! Newton's Method variables
  real :: dT_dP, dS_dP, delta_T, delta_S, delta_P
  real :: dbeta_dS, dbeta_dT, dalpha_dT, dalpha_dS, dbeta_dP, dalpha_dP
  ! Brent's Method variables
  real :: a, b, c, d, e, f, fa, fb, fc, m, p, q, r, s0, sa, sb, tol, machep

  real :: P_last
  logical :: debug = .false.
  if (ref_pres>=0.) P_ref = ref_pres
  delta_P = P_bot-P_top
  refine_nondim_position = min_bound
  x_init = refine_nondim_position

  call extract_member_EOS(EOS, form_of_EOS = form_of_EOS)
  do_newton = (form_of_EOS == EOS_LINEAR) .or. (form_of_EOS == EOS_TEOS10) .or. (form_of_EOS == EOS_WRIGHT)
  do_brent = .not. do_newton
  if (present(force_brent)) then
    do_newton = .not. force_brent
    do_brent = force_brent
  endif

  ! Check to make sure that a root exists between the minimum bound and the bottom of the layer
  call calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, refine_nondim_position, &
                      ref_pres, EOS, delta_rho)
  delta_rho_init = delta_rho
  if ( SIGN(1.,delta_rho) == SIGN(1.,drho_bot) ) then
    ! Return the position of min_bound if closer to 0 than drho_bot
    if (ABS(delta_rho) < ABS(drho_bot)) then
      refine_nondim_position = min_bound
    else
      refine_nondim_position = 1.
    endif
    do_newton = .false. ; do_brent = .false.
  endif

  if (debug) then
    write (*,*) "------"
    write (*,*) "Starting delta_rho: ", delta_rho
  endif

  ! For now only linear, Wright, and TEOS-10 equations of state have functions providing second derivatives and
  ! thus can use Newton's method for the equation of state
  if (do_newton) then
    ! Set lower bound of pressure
    P_min = P_top*(1.-min_bound) + P_bot*(min_bound)
    ! Iterate over Newton's method for the function: x0 = x0 - delta_rho/d_delta_rho_dP
    do iter = 1, max_iter
      ! Evaluate delta_rho(x0)
      call calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S,   &
                          refine_nondim_position, ref_pres, EOS, delta_rho, P_int, T, S, alpha_avg, &
                          beta_avg, delta_T, delta_S)
      ! Check for convergence
      if (ABS(delta_rho) <= tolerance) then
        do_brent = .false.
        exit
      endif
      ! Evaluate total derivative of delta_rho
      if (ref_pres<0.) P_ref = P_int
      call calculate_density_second_derivs( T, S, P_ref, dbeta_dS, dbeta_dT, dalpha_dT, dbeta_dP, dalpha_dP, EOS )
      ! In the case of a constant reference pressure, no dependence on neutral direction with pressure
      if (ref_pres>=0.) then
        dalpha_dP = 0. ; dbeta_dP = 0.
      endif
      dalpha_dS = dbeta_dT ! Cross derivatives are identicial
      ! By chain rule dT_dP= (dT_dz)*(dz/dP) = dT_dz / (Pbot-Ptop)
      dT_dP = first_derivative_polynomial( ppoly_T, deg+1, refine_nondim_position ) / delta_P
      dS_dP = first_derivative_polynomial( ppoly_S, deg+1, refine_nondim_position ) / delta_P
      ! Total derivative of d_delta_rho wrt P
      d_delta_rho_dP = 0.5*( delta_S*(dS_dP*dbeta_dS + dT_dP*dbeta_dT + dbeta_dP) +     &
                             ( delta_T*(dS_dP*dalpha_dS + dT_dP*dalpha_dT + dalpha_dP))) + &
                             dS_dP*beta_avg + dT_dP*alpha_avg
      if (d_delta_rho_dP == 0.) then
        do_brent = .true.
        exit
      endif
      ! Newton step update
      P_last = P_int
      P_int = P_int - (delta_rho / d_delta_rho_dP)
      if (P_int < P_min .or. P_int > P_bot) then
        if (debug) then
          write (*,*) "Iteration: ", iter
          write (*,*) "delta_rho, d_delta_rho_dP: ", delta_rho, d_delta_rho_dP
          write (*,*) "T, T Poly Coeffs: ", T, ppoly_T
          write (*,*) "S, S Poly Coeffs: ", S, ppoly_S
          write (*,*) "T_ref, alpha_ref: ", T_ref, alpha_ref
          write (*,*) "S_ref, beta_ref : ", S_ref, beta_ref
          write (*,*) "P, dT_dP, dS_dP:", P_int, dT_dP, dS_dP
          write (*,*) "dRhoTop, dRhoBot:", drho_top, drho_bot
          write (*,*) "x0: ", x0
          write (*,*) "refine_nondim_position: ", refine_nondim_position
          write (*,*)
        endif
!        call MOM_error(WARNING, "Step went out of bounds")
        ! Switch to Brent's method by setting the converged flag to false
        do_brent = .true.
        ! Reset to first guess if already diverged
        if (ABS(delta_rho_init)<ABS(delta_rho)) then
          refine_nondim_position = x_init
        endif
        exit
      endif
      refine_nondim_position = (P_top-P_int)/(P_top-P_bot)
      ! Check to see if the updated position is too small
!      if ( ABS(P_last-P_int) < 2*EPSILON(P_int)*MAX(P_last,P_int) ) exit
    enddo
  endif
  ! Do Brent if Newton's method kicked out or if second derivatives don't exist
  if (do_brent) then
    machep = EPSILON(sa)
    ! Find the bracketing interval based on where the sign changes
    call calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, &
                        refine_nondim_position, ref_pres, EOS, delta_rho)
    ! Bracketed between min_bound and 1.
    if ( SIGN(1.,delta_rho)*SIGN(1.,drho_bot)==-1. ) then
      sa = refine_nondim_position ; fa = delta_rho
      sb = 1. ; fb = drho_bot
    elseif ( SIGN(1.,delta_rho)*SIGN(1.,drho_top)==-1. ) then
      sa = 0. ; fa = drho_top
      sb = refine_nondim_position ; fb = delta_rho
    else
      call MOM_error(FATAL, "refine_nondim_position: No root exists in Brent's method interval")
    endif
    c = sa ; fc = fa ; e = sb - sa; d = e

    ! This is from https://people.sc.fsu.edu/~jburkardt/f_src/brent/brent.f90
    do iter = 1,max_iter
      if ( abs ( fc ) < abs ( fb ) ) then
        sa = sb
        sb = c
        c = sa
        fa = fb
        fb = fc
        fc = fa
      end if
      tol = 2. * machep * abs ( sb ) + tolerance
      m = 0.5 * ( c - sb )
      if ( abs ( m ) <= tol .or. fb == 0. ) then
        exit
      end if
      if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then
        e = m
        d = e
      else
        s0 = fb / fa
        if ( sa == c ) then
          p = 2. * m * s0
          q = 1. - s0
        else
          q = fa / fc
          r = fb / fc
          p = s0 * ( 2. * m * q * ( q - r ) - ( sb - sa ) * ( r - 1. ) )
          q = ( q - 1. ) * ( r - 1. ) * ( s0 - 1. )
        end if
        if ( 0. < p ) then
          q = - q
        else
          p = - p
        end if
        s0 = e
        e = d
        if ( 2. * p < 3. * m * q - abs ( tol * q ) .and. &
          p < abs ( 0.5 * s0 * q ) ) then
          d = p / q
        else
          e = m
          d = e
        end if
      end if
      sa = sb
      fa = fb
      if ( tol < abs ( d ) ) then
        sb = sb + d
      else if ( 0. < m ) then
        sb = sb + tol
      else
        sb = sb - tol
      end if
      call calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, &
                        sb, ref_pres, EOS, fb)
      if ( ( 0. < fb .and. 0. < fc ) .or. &
           ( fb <= 0. .and. fc <= 0. ) ) then
        c = sa
        fc = fa
        e = sb - sa
        d = e
      end if
    enddo
    ! Modified from original to ensure that the minimum position is found
    fa = ABS(fa) ; fb = ABS(fb) ; fc = ABS(fc)
    delta_rho = MIN(fa, fb, fc)
    if (fb==delta_rho) then
      refine_nondim_position = sb
    elseif (fa==delta_rho) then
      refine_nondim_position = sa
    elseif (fc==delta_rho) then
      refine_nondim_position = c
    endif
  endif

  ! Make sure that the result is bounded between 0 and 1
  if (refine_nondim_position>1.) then
    if (debug) then
      write (*,*) "T, T Poly Coeffs: ", T, ppoly_T
      write (*,*) "S, S Poly Coeffs: ", S, ppoly_S
      write (*,*) "T_ref, alpha_ref: ", T_ref, alpha_ref
      write (*,*) "S_ref, beta_ref : ", S_ref, beta_ref
      write (*,*) "P, dT_dP, dS_dP:", P_int, dT_dP, dS_dP
      write (*,*) "x0: ", x0
      write (*,*) "refine_nondim_position: ", refine_nondim_position
    endif
    call MOM_error(WARNING, "refine_nondim_position>1.")
    refine_nondim_position = MAX(x0,min_bound)
  endif

  if (refine_nondim_position<0.) then
    if (debug) then
      write (*,*) "T, T Poly Coeffs: ", T, ppoly_T
      write (*,*) "S, S Poly Coeffs: ", S, ppoly_S
      write (*,*) "T_ref, alpha_ref: ", T_ref, alpha_ref
      write (*,*) "S_ref, beta_ref : ", S_ref, beta_ref
      write (*,*) "dT_dP, dS_dP:", dT_dP, dS_dP
      write (*,*) "x0: ", x0
      write (*,*) "refine_nondim_position: ", refine_nondim_position
    endif
    call MOM_error(WARNING, "refine_nondim_position<0.")
    refine_nondim_position = MAX(x0,min_bound)
  endif

  if (debug) then
    call calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, &
                        refine_nondim_position, ref_pres, EOS, delta_rho)
    write (*,*) "End delta_rho: ", delta_rho
    write (*,*) "x0, delta_x: ", x0, refine_nondim_position-x0
    write (*,*) "Iterations: ", iter
    write (*,*) "******"
  endif

end function refine_nondim_position

!> Calculate the difference in neutral density between a reference T, S, alpha, and beta
!! and a point on the polynomial reconstructions of T, S
subroutine calc_delta_rho(deg, T_ref, S_ref, alpha_ref, beta_ref, P_top, P_bot, ppoly_T, ppoly_S, x0, ref_pres, EOS, &
                          delta_rho, P_out, T_out, S_out, alpha_avg_out, beta_avg_out, delta_T_out, delta_S_out)
  integer,                intent(in)  :: deg       !< Degree of polynomial reconstruction
  real,                   intent(in)  :: T_ref     !< Temperature at reference surface
  real,                   intent(in)  :: S_ref     !< Salinity at reference surface
  real,                   intent(in)  :: alpha_ref !< dRho/dT at reference surface
  real,                   intent(in)  :: beta_ref  !< dRho/dS at reference surface
  real,                   intent(in)  :: P_top     !< Pressure (Pa) at top interface of layer to be searched
  real,                   intent(in)  :: P_bot     !< Pressure (Pa) at bottom interface
  real, dimension(deg+1), intent(in)  :: ppoly_T   !< Coefficients of T reconstruction
  real, dimension(deg+1), intent(in)  :: ppoly_S   !< Coefficients of S reconstruciton
  real,                   intent(in)  :: x0        !< Nondimensional position to evaluate
  real,                   intent(in)  :: ref_pres  !< Reference pressure
  type(EOS_type),         pointer     :: EOS       !< Equation of state structure
  real,                   intent(out) :: delta_rho
  real,         optional, intent(out) :: P_out         !< Pressure at point x0
  real,         optional, intent(out) :: T_out         !< Temperature at point x0
  real,         optional, intent(out) :: S_out         !< Salinity at point x0
  real,         optional, intent(out) :: alpha_avg_out !< Average of alpha between reference and x0
  real,         optional, intent(out) :: beta_avg_out  !< Average of beta between reference and x0
  real,         optional, intent(out) :: delta_T_out   !< Difference in temperature between reference and x0
  real,         optional, intent(out) :: delta_S_out   !< Difference in salinity between reference and x0

  real :: alpha, beta, alpha_avg, beta_avg, P_int, T, S, delta_T, delta_S

  P_int = (1. - x0)*P_top + x0*P_bot
  T = evaluation_polynomial( ppoly_T, deg+1, x0 )
  S = evaluation_polynomial( ppoly_S, deg+1, x0 )
  ! Interpolated pressure if using locally referenced neutral density
  if (ref_pres<0.) then
    call calculate_density_derivs( T, S, P_int, alpha, beta, EOS )
  else
  ! Constant reference pressure (isopycnal)
    call calculate_density_derivs( T, S, ref_pres, alpha, beta, EOS )
  endif

  ! Calculate the f(P) term for Newton's method
  alpha_avg = 0.5*( alpha + alpha_ref )
  beta_avg = 0.5*( beta + beta_ref )
  delta_T = T - T_ref
  delta_S = S - S_ref
  delta_rho = alpha_avg*delta_T + beta_avg*delta_S

  ! If doing a Newton step, these quantities are needed, otherwise they can just be optional
  if (present(P_out)) P_out = P_int
  if (present(T_out)) T_out = T
  if (present(S_out)) S_out = S
  if (present(alpha_avg_out)) alpha_avg_out = alpha_avg
  if (present(beta_avg_out))  beta_avg_out = beta_avg
  if (present(delta_T_out)) delta_T_out = delta_T
  if (present(delta_S_out)) delta_S_out = delta_S

end subroutine calc_delta_rho

!> Returns a single column of neutral diffusion fluxes of a tracer.
subroutine neutral_surface_flux(nk, nsurf, deg, hl, hr, Tl, Tr, PiL, PiR, KoL, KoR, hEff, Flx, continuous, remap_CS)
  integer,                      intent(in)    :: nk    !< Number of levels
  integer,                      intent(in)    :: nsurf !< Number of neutral surfaces
  integer,                      intent(in)    :: deg   !< Degree of polynomial reconstructions
  real, dimension(nk),          intent(in)    :: hl    !< Left-column layer thickness (Pa)
  real, dimension(nk),          intent(in)    :: hr    !< Right-column layer thickness (Pa)
  real, dimension(nk),          intent(in)    :: Tl    !< Left-column layer tracer (conc, e.g. degC)
  real, dimension(nk),          intent(in)    :: Tr    !< Right-column layer tracer (conc, e.g. degC)
  real, dimension(nsurf),       intent(in)    :: PiL   !< Fractional position of neutral surface
                                                       !! within layer KoL of left column
  real, dimension(nsurf),       intent(in)    :: PiR   !< Fractional position of neutral surface
                                                       !! within layer KoR of right column
  integer, dimension(nsurf),    intent(in)    :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(nsurf),    intent(in)    :: KoR   !< Index of first right interface above neutral surface
  real, dimension(nsurf-1),     intent(in)    :: hEff  !< Effective thickness between two neutral surfaces (Pa)
  real, dimension(nsurf-1),     intent(inout) :: Flx   !< Flux of tracer between pairs of neutral layers (conc H)
  logical,                      intent(in)    :: continuous !< True if using continuous reconstruction
  type(remapping_CS), optional, intent(in)    :: remap_CS
  ! Local variables
  integer :: k_sublayer, klb, klt, krb, krt, k
  real :: T_right_top, T_right_bottom, T_right_layer
  real :: T_left_top, T_left_bottom, T_left_layer
  real :: dT_top, dT_bottom, dT_layer, dT_ave
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

  ! Setup reconstruction edge values
  if (continuous) then
    call interface_scalar(nk, hl, Tl, Til, 2)
    call interface_scalar(nk, hr, Tr, Tir, 2)
    call ppm_left_right_edge_values(nk, Tl, Til, aL_l, aR_l)
    call ppm_left_right_edge_values(nk, Tr, Tir, aL_r, aR_r)
  else
    ppoly_r_coeffs_l(:,:) = 0.
    ppoly_r_coeffs_r(:,:) = 0.
    Tid_l(:,:) = 0.
    Tid_r(:,:) = 0.

    call build_reconstructions_1d( remap_CS, nk, hl, Tl, ppoly_r_coeffs_l, Tid_l, ppoly_r_S_l, iMethod )
    call build_reconstructions_1d( remap_CS, nk, hr, Tr, ppoly_r_coeffs_r, Tid_r, ppoly_r_S_r, iMethod )
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
      else ! Discontinuous reconstruction
        klb = KoL(k_sublayer+1)
        klt = KoL(k_sublayer)
        if (klt .ne. klb) then
          call MOM_error(WARNING, "Neutral surfaces span more than one layer")
          Flx(k_sublayer) = 0.
          cycle
        endif
        T_left_bottom = evaluation_polynomial( ppoly_r_coeffs_l(klb,:), deg+1, PiL(k_sublayer+1))
        T_left_top = evaluation_polynomial( ppoly_r_coeffs_l(klt,:), deg+1, PiL(k_sublayer))
        T_left_layer = average_value_ppoly(nk, Tl, Tid_l, ppoly_r_coeffs_l, iMethod, klb, &
                                           PiL(k_sublayer), PiL(k_sublayer+1))
        krb = KoR(k_sublayer+1)
        krt = KoR(k_sublayer)
        if (krt .ne. krb) then
          call MOM_error(WARNING, "Neutral surfaces span more than one layer")
          Flx(k_sublayer) = 0.
          cycle
        endif
        T_right_bottom = evaluation_polynomial( ppoly_r_coeffs_r(krb,:), deg+1, PiR(k_sublayer+1))
        T_right_top = evaluation_polynomial( ppoly_r_coeffs_r(krt,:), deg+1, PiR(k_sublayer))
        T_right_layer = average_value_ppoly(nk, Tr, Tid_r, ppoly_r_coeffs_r, iMethod, krb, &
                                            PiR(k_sublayer), PiR(k_sublayer+1))
      endif
      dT_top = T_right_top - T_left_top
      dT_bottom = T_right_bottom - T_left_bottom
      dT_ave = 0.5 * ( dT_top + dT_bottom )
      dT_layer = T_right_layer - T_left_layer
      if (signum(1.,dT_top) * signum(1.,dT_bottom) <= 0. .or. signum(1.,dT_ave) * signum(1.,dT_layer) <= 0. ) then
        dT_ave = 0.
      else
        dT_ave = dT_layer
      endif
      Flx(k_sublayer) = dT_ave * hEff(k_sublayer)
    endif
  enddo

end subroutine neutral_surface_flux

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
  logical, intent(in) :: verbose

  neutral_diffusion_unit_tests = .false. .or. &
    ndiff_unit_tests_continuous(verbose) .or. ndiff_unit_tests_discontinuous(verbose)


end function neutral_diffusion_unit_tests

!> Returns true if unit tests of neutral_diffusion functions fail. Otherwise returns false.
logical function ndiff_unit_tests_continuous(verbose)
  logical, intent(in) :: verbose !< It true, write results to stdout
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

  v = verbose

  ndiff_unit_tests_continuous = .false. ! Normally return false
  write(*,*) '==== MOM_neutral_diffusion: ndiff_unit_tests_continuous ='

  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fv_diff(v,1.,1.,1., 0.,1.,2., 1., 'FV: Straight line on uniform grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fv_diff(v,1.,1.,0., 0.,4.,8., 7., 'FV: Vanished right cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fv_diff(v,0.,1.,1., 0.,4.,8., 7., 'FV: Vanished left cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fv_diff(v,1.,2.,4., 0.,3.,9., 4., 'FV: Stretched grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fv_diff(v,2.,0.,2., 0.,1.,2., 0., 'FV: Vanished middle cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fv_diff(v,0.,1.,0., 0.,1.,2., 2., 'FV: Vanished on both sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fv_diff(v,1.,0.,0., 0.,1.,2., 0., 'FV: Two vanished cell sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fv_diff(v,0.,0.,0., 0.,1.,2., 0., 'FV: All vanished cells')

  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fvlsq_slope(v,1.,1.,1., 0.,1.,2., 1., 'LSQ: Straight line on uniform grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fvlsq_slope(v,1.,1.,0., 0.,1.,2., 1., 'LSQ: Vanished right cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fvlsq_slope(v,0.,1.,1., 0.,1.,2., 1., 'LSQ: Vanished left cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fvlsq_slope(v,1.,2.,4., 0.,3.,9., 2., 'LSQ: Stretched grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fvlsq_slope(v,1.,0.,1., 0.,1.,2., 2., 'LSQ: Vanished middle cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fvlsq_slope(v,0.,1.,0., 0.,1.,2., 0., 'LSQ: Vanished on both sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fvlsq_slope(v,1.,0.,0., 0.,1.,2., 0., 'LSQ: Two vanished cell sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_fvlsq_slope(v,0.,0.,0., 0.,1.,2., 0., 'LSQ: All vanished cells')

  call interface_scalar(4, (/10.,10.,10.,10./), (/24.,18.,12.,6./), Tio, 1)
  !ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(5, Tio, (/27.,21.,15.,9.,3./), 'Linear profile, interface temperatures')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v,5, Tio, (/24.,22.5,15.,7.5,6./), 'Linear profile, linear interface temperatures')
  call interface_scalar(4, (/10.,10.,10.,10./), (/24.,18.,12.,6./), Tio, 2)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v,5, Tio, (/24.,22.,15.,8.,6./), 'Linear profile, PPM interface temperatures')

  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v,-1.0, 0.,  1.0, 1.0, 0.5, 'Check mid-point')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v, 0.0, 0.,  1.0, 1.0, 0.0, 'Check bottom')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v, 0.1, 0.,  1.1, 1.0, 0.0, 'Check below')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v,-1.0, 0.,  0.0, 1.0, 1.0, 'Check top')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v,-1.0, 0., -0.1, 1.0, 1.0, 'Check above')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v,-1.0, 0.,  3.0, 1.0, 0.25, 'Check 1/4')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v,-3.0, 0.,  1.0, 1.0, 0.75, 'Check 3/4')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v, 1.0, 0.,  1.0, 1.0, 0.0, 'Check dRho=0 below')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v,-1.0, 0., -1.0, 1.0, 1.0, 'Check dRho=0 above')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v, 0.0, 0.,  0.0, 1.0, 0.5, 'Check dRho=0 mid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_ifndp(v,-2.0, .5,  5.0, 0.5, 0.5, 'Check dP=0')

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
                               PiLRo, PiRLo, KoL, KoR, hEff, Flx, .true.)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 7, Flx, &
              (/0.,0.,0.,0.,0.,0.,0./), 'Identical columns, rho flux (=0)')
  call neutral_surface_flux(3, 2*3+2, 2, (/10.,10.,10./), (/10.,10.,10./), & ! nk, hL, hR
                               (/-1.,-1.,-1./), (/1.,1.,1./), & ! Sl, Sr
                               PiLRo, PiRLo, KoL, KoR, hEff, Flx, .true.)
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
                                   'Indentical columns with mixed layer')

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

  if (.not. ndiff_unit_tests_continuous) write(*,*) 'Pass'

end function ndiff_unit_tests_continuous

logical function ndiff_unit_tests_discontinuous(verbose)
  logical, intent(in) :: verbose !< It true, write results to stdout
  ! Local variables
  integer, parameter          :: nk = 3
  integer, parameter          :: ns = nk*4
  real, dimension(nk)         :: Sl, Sr, Tl, Tr, hl, hr
  real, dimension(nk,2)       :: TiL, SiL, TiR, SiR
  real, dimension(nk+1)       :: Pres_l, Pres_R
  integer, dimension(ns)      :: KoL, KoR
  real, dimension(ns)         :: PoL, PoR
  real, dimension(ns-1)       :: hEff, Flx
  type(EOS_type),     pointer :: EOS       ! Structure for linear equation of state
  type(remapping_CS), pointer :: remap_CS  ! Remapping control structure (PLM)
  real, dimension(nk,2)       :: poly_T_l, poly_T_r, poly_S,  poly_slope    ! Linear reconstruction for T
  real, dimension(nk,2)       :: dRdT, dRdS
  integer                     :: iMethod
  integer :: k
  logical :: v

  v = verbose

  ndiff_unit_tests_discontinuous = .false. ! Normally return false
  write(*,*) '==== MOM_neutral_diffusion: ndiff_unit_tests_discontinuous ='

  ! Unit tests for find_neutral_surface_positions_discontinuous
  ! Salinity is 0 for all these tests
  Sl(:) = 0. ; Sr(:) = 0. ; poly_S(:,:) = 0. ; SiL(:,:) = 0. ; SiR(:,:) = 0.
  dRdT(:,:) = -1. ; dRdS(:,:) = 0.
  allocate(remap_CS)
  call initialize_remapping( remap_CS, "PLM", boundary_extrapolation = .true. )

  hL = (/10.,10.,10./) ; hR = (/10.,10.,10./) ; Pres_l(1) = 0. ; Pres_r(1) = 0.
  do k = 1,nk ; Pres_l(k+1) = Pres_l(k) + hL(k) ; Pres_r(k+1) = Pres_r(k) + hR(k) ; enddo
  ! Identical columns
  Tl = (/20.,16.,12./) ; Tr = (/20.,16.,12./)
  call build_reconstructions_1d( remap_CS, nk, hL, Tl, poly_T_l, TiL, poly_slope, iMethod )
  call build_reconstructions_1d( remap_CS, nk, hR, Tr, poly_T_r, TiR, poly_slope, iMethod )
  call find_neutral_surface_positions_discontinuous(nk, 1, Pres_l, TiL, SiL, dRdT, dRdS, &
            Pres_r, TiR, SiR, dRdT, dRdS, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
                                   (/1,1,1,1,2,2,2,2,3,3,3,3/), & ! KoL
                                   (/1,1,1,1,2,2,2,2,3,3,3,3/), & ! KoR
                                   (/0.,0.,1.,1.,0.,0.,1.,1.,0.,0.,1.,1./), & ! pL
                                   (/0.,0.,1.,1.,0.,0.,1.,1.,0.,0.,1.,1./), & ! pR
                                   (/0.,10.,0.,0.,0.,10.,0.,0.,0.,10.,0.,0./), & ! hEff
                                   'Identical columns')
  Tl = (/20.,16.,12./) ; Tr = (/18.,14.,10./)
  call build_reconstructions_1d( remap_CS, nk, hL, Tl, poly_T_l, TiL, poly_slope, iMethod )
  call build_reconstructions_1d( remap_CS, nk, hR, Tr, poly_T_r, TiR, poly_slope, iMethod )
  call find_neutral_surface_positions_discontinuous(nk, 1, Pres_l, TiL, SiL, dRdT, dRdS, &
            Pres_r, TiR, SiR, dRdT, dRdS, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
                                   (/1,1,1,2,2,2,2,3,3,3,3,3/), & ! KoL
                                   (/1,1,1,1,1,2,2,2,2,3,3,3/),  & ! KoR
                                   (/0.0, 0.5, 1.0, 0.0, 0.5, 0.5, 1.0, 0.0, 0.5, 0.5, 1.0, 1.0/), & ! pL
                                   (/0.0, 0.0, 0.5, 0.5, 1.0, 0.0, 0.5, 0.5, 1.0, 0.0, 0.5, 1.0/), & ! pR
                                   (/0.0, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0/), & ! hEff
                                   'Right column slightly cooler')
  Tl = (/18.,14.,10./) ; Tr = (/20.,16.,12./) ;
  call build_reconstructions_1d( remap_CS, nk, hL, Tl, poly_T_l, TiL, poly_slope, iMethod )
  call build_reconstructions_1d( remap_CS, nk, hR, Tr, poly_T_r, TiR, poly_slope, iMethod )
  call find_neutral_surface_positions_discontinuous(nk, 1, Pres_l, TiL, SiL, dRdT, dRdS, &
            Pres_r, TiR, SiR, dRdT, dRdS, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
                                   (/1,1,1,1,1,2,2,2,2,3,3,3/),  & ! KoL
                                   (/1,1,1,2,2,2,2,3,3,3,3,3/), & ! KoR
                                   (/0.0, 0.0, 0.5, 0.5, 1.0, 0.0, 0.5, 0.5, 1.0, 0.0, 0.5, 1.0/), & ! pL
                                   (/0.0, 0.5, 1.0, 0.0, 0.5, 0.5, 1.0, 0.0, 0.5, 0.5, 1.0, 1.0/), & ! pR
                                   (/0.0, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0/), & ! hEff
                                   'Left column slightly cooler')
  Tl = (/20.,16.,12./) ; Tr = (/14.,10.,6./)
  call build_reconstructions_1d( remap_CS, nk, hL, Tl, poly_T_l, TiL, poly_slope, iMethod )
  call build_reconstructions_1d( remap_CS, nk, hR, Tr, poly_T_r, TiR, poly_slope, iMethod )
  call find_neutral_surface_positions_discontinuous(nk, 1, Pres_l, TiL, SiL, dRdT, dRdS, &
            Pres_r, TiR, SiR, dRdT, dRdS, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
                                   (/1,1,2,2,2,3,3,3,3,3,3,3/),  & ! KoL
                                   (/1,1,1,1,1,1,1,2,2,2,3,3/), & ! KoR
                                   (/0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0/), & ! pL
                                   (/0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 1.0/), & ! pR
                                   (/0.0, 0.0, 0.0, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0, 0.0, 0.0/), & ! hEff
                                   'Right column somewhat cooler')
  Tl = (/20.,16.,12./) ; Tr = (/8.,6.,4./)
  call build_reconstructions_1d( remap_CS, nk, hL, Tl, poly_T_l, TiL, poly_slope, iMethod )
  call build_reconstructions_1d( remap_CS, nk, hR, Tr, poly_T_r, TiR, poly_slope, iMethod )
  call find_neutral_surface_positions_discontinuous(nk, 1, Pres_l, TiL, SiL, dRdT, dRdS, &
            Pres_r, TiR, SiR, dRdT, dRdS, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
                                   (/1,1,2,2,3,3,3,3,3,3,3,3/),  & ! KoL
                                   (/1,1,1,1,1,1,1,1,2,2,3,3/), & ! KoR
                                   (/0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/), & ! pL
                                   (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0/), & ! pR
                                   (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/), & ! hEff
                                   'Right column much cooler')
  Tl = (/14.,14.,10./) ; Tr = (/14.,14.,10./)
  call build_reconstructions_1d( remap_CS, nk, hL, Tl, poly_T_l, TiL, poly_slope, iMethod )
  call build_reconstructions_1d( remap_CS, nk, hR, Tr, poly_T_r, TiR, poly_slope, iMethod )
  call find_neutral_surface_positions_discontinuous(nk, 1, Pres_l, TiL, SiL, dRdT, dRdS, &
            Pres_r, TiR, SiR, dRdT, dRdS, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
                                   (/1,1,1,1,2,2,2,2,3,3,3,3/), & ! KoL
                                   (/1,1,1,1,2,2,2,2,3,3,3,3/), & ! KoR
                                   (/0.,0.,1.,1.,0.,0.,1.,1.,0.,0.,1.,1./), & ! pL
                                   (/0.,0.,1.,1.,0.,0.,1.,1.,0.,0.,1.,1./), & ! pR
                                   (/0.,10.,0.,0.,0.,10.,0.,0.,0.,10.,0.,0./), & ! hEff
                                   'Identical columns with mixed layer')
  Tl = (/20.,16.,12./) ; Tr = (/14.,14.,10./)
  call build_reconstructions_1d( remap_CS, nk, hL, Tl, poly_T_l, TiL, poly_slope, iMethod )
  call build_reconstructions_1d( remap_CS, nk, hR, Tr, poly_T_r, TiR, poly_slope, iMethod )
  call find_neutral_surface_positions_discontinuous(nk, 1, Pres_l, TiL, SiL, dRdT, dRdS, &
            Pres_r, TiR, SiR, dRdT, dRdS, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
                                   (/1,1,2,2,2,3,3,3,3,3,3,3/),  & ! KoL
                                   (/1,1,1,1,1,1,2,2,2,3,3,3/), & ! KoR
                                   (/0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0/), & ! pL
                                   (/0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 1.0/), & ! pR
                                   (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0/), & ! hEff
                                   'Right column with mixed layer')
  Tl = (/14.,14.,6./) ; Tr = (/12.,16.,8./)
  call build_reconstructions_1d( remap_CS, nk, hL, Tl, poly_T_l, TiL, poly_slope, iMethod )
  call build_reconstructions_1d( remap_CS, nk, hR, Tr, poly_T_r, TiR, poly_slope, iMethod )
  call find_neutral_surface_positions_discontinuous(nk, 1, Pres_l, TiL, SiL, dRdT, dRdS, &
            Pres_r, TiR, SiR, dRdT, dRdS, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
                                   (/1,1,2,2,3,3,3,3,3,3,3,3/), & ! KoL
                                   (/1,1,1,1,1,1,2,2,3,3,3,3/), & ! KoR
                                   (/0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .75, 1.0/), & ! pL
                                   (/0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, .25, 1.0, 1.0/), & ! pR
                                   (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.5, 0.0/), & ! hEff
                                   'Left mixed layer, right unstable mixed layer')
  Til(:,1) = (/8.,12.,10./) ; Til(:,2) = (/12.,10.,2./)
  Tir(:,1) = (/10.,14.,12./) ; TiR(:,2) = (/14.,12.,4./)
  call find_neutral_surface_positions_discontinuous(nk, 1, Pres_l, TiL, SiL, dRdT, dRdS, &
            Pres_r, TiR, SiR, dRdT, dRdS, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
                                   (/1,1,1,1,1,1,1,2,2,3,3,3/), & ! KoL
                                   (/1,1,2,2,3,3,3,3,3,3,3,3/), & ! KoR
                                   (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, .75, 1.0/), & ! pL
                                   (/0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, .25, .25, 1.0, 1.0/), & ! pR
                                   (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 7.5, 0.0/), & ! hEff
                                   'Two unstable mixed layers')
  deallocate(remap_CS)

  allocate(EOS)
  call EOS_manual_init(EOS, form_of_EOS = EOS_LINEAR, dRho_dT = -1., dRho_dS = 2.)
  ! Unit tests for refine_nondim_position
  ! Tests using Newton's method
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5,refine_nondim_position( &
            100, 0., 20., 35., -1., 2., 0., 1., 1, (/21., -2./), (/35., 0./), EOS, 0., -1., 1., 0., -1.), &
            "Temperature stratified (Newton) "))
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5,refine_nondim_position( &
            100, 0., 20., 35., -1., 2., 0., 1., 1, (/20., 0./), (/34., 2./), EOS, 0., -2., 2., 0., -1.), &
            "Salinity stratified    (Newton) "))
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5,refine_nondim_position( &
            100, 0., 20., 35., -1., 2., 0., 1., 1, (/21., -2./), (/34., 2./), EOS, 0., -1., 1., 0., -1.), &
            "Temp/Salt stratified   (Newton) "))
  ! Tests using Brent's method
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5,refine_nondim_position( &
            100, 0., 20., 35., -1., 2., 0., 1., 1, (/21., -2./), (/35., 0./), EOS, 0., -1., 1., 0., -1., force_brent = .true.), &
            "Temperature stratified (Brent)  "))
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5,refine_nondim_position( &
            100, 0., 20., 35., -1., 2., 0., 1., 1, (/20., 0./), (/34., 2./), EOS, 0., -2., 2., 0., -1., force_brent = .true.), &
            "Salinity stratified    (Brent)  "))
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5,refine_nondim_position( &
            100, 0., 20., 35., -1., 2., 0., 1., 1, (/21., -2./), (/34., 2./), EOS, 0., -1., 1., 0., -1., force_brent = .true.), &
            "Temp/Salt stratified   (Brent)  "))
  deallocate(EOS)

  if (.not. ndiff_unit_tests_discontinuous) write(*,*) 'Pass'

end function ndiff_unit_tests_discontinuous

!> Returns true if a test of fv_diff() fails, and conditionally writes results to stream
logical function test_fv_diff(verbose, hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
  logical,          intent(in) :: verbose !< If true, write results to stdout
  real,             intent(in) :: hkm1  !< Left cell width
  real,             intent(in) :: hk    !< Center cell width
  real,             intent(in) :: hkp1  !< Right cell width
  real,             intent(in) :: Skm1  !< Left cell average value
  real,             intent(in) :: Sk    !< Center cell average value
  real,             intent(in) :: Skp1  !< Right cell average value
  real,             intent(in) :: Ptrue !< True answer (Pa)
  character(len=*), intent(in) :: title !< Title for messages

  ! Local variables
  integer :: stdunit
  real :: Pret

  Pret = fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  test_fv_diff = (Pret /= Ptrue)

  if (test_fv_diff .or. verbose) then
    stdunit = 6
    if (test_fv_diff) stdunit = 0 ! In case of wrong results, write to error stream
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
  real,             intent(in) :: Ptrue !< True answer (Pa)
  character(len=*), intent(in) :: title !< Title for messages

  ! Local variables
  integer :: stdunit
  real :: Pret

  Pret = fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  test_fvlsq_slope = (Pret /= Ptrue)

  if (test_fvlsq_slope .or. verbose) then
    stdunit = 6
    if (test_fvlsq_slope) stdunit = 0 ! In case of wrong results, write to error stream
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
  real,             intent(in) :: rhoNeg !< Lighter density (kg/m3)
  real,             intent(in) :: Pneg   !< Interface position of lighter density (Pa)
  real,             intent(in) :: rhoPos !< Heavier density (kg/m3)
  real,             intent(in) :: Ppos   !< Interface position of heavier density (Pa)
  real,             intent(in) :: Ptrue  !< True answer (Pa)
  character(len=*), intent(in) :: title  !< Title for messages

  ! Local variables
  integer :: stdunit
  real :: Pret

  Pret = interpolate_for_nondim_position(rhoNeg, Pneg, rhoPos, Ppos)
  test_ifndp = (Pret /= Ptrue)

  if (test_ifndp .or. verbose) then
    stdunit = 6
    if (test_ifndp) stdunit = 0 ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    if (test_ifndp) then
      write(stdunit,'(4(x,a,f20.16),2(x,a,1pe22.15),x,a)') 'r1=',rhoNeg,'p1=',Pneg,'r2=',rhoPos,'p2=',Ppos,'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
    else
      write(stdunit,'(4(x,a,f20.16),2(x,a,1pe22.15))') 'r1=',rhoNeg,'p1=',Pneg,'r2=',rhoPos,'p2=',Ppos,'pRet=',Pret,'pTrue=',Ptrue
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
    stdunit = 6
    if (test_data1d) stdunit = 0 ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    do k = 1,nk
      if (Po(k) /= Ptrue(k)) then
        test_data1d = .true.
        write(stdunit,'(a,i2,2(x,a,f20.16),x,a,1pe22.15,x,a)') 'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k),'WRONG!'
      else
        if (verbose) &
          write(stdunit,'(a,i2,2(x,a,f20.16),x,a,1pe22.15)') 'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k)
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
    stdunit = 6
    if (test_data1di) stdunit = 0 ! In case of wrong results, write to error stream
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

!> Returns true if output of find_neutral_surface_positions() does not match correct values, and conditionally writes results to stream
logical function test_nsp(verbose, ns, KoL, KoR, pL, pR, hEff, KoL0, KoR0, pL0, pR0, hEff0, title)
  logical,                    intent(in) :: verbose !< If true, write results to stdout
  integer,                    intent(in) :: ns    !< Number of surfaces
  integer, dimension(ns), intent(in) :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(ns), intent(in) :: KoR   !< Index of first right interface above neutral surface
  real, dimension(ns),    intent(in) :: pL    !< Fractional position of neutral surface within layer KoL of left column
  real, dimension(ns),    intent(in) :: pR    !< Fractional position of neutral surface within layer KoR of right column
  real, dimension(ns-1),    intent(in) :: hEff  !< Effective thickness between two neutral surfaces (Pa)
  integer, dimension(ns), intent(in) :: KoL0  !< Correct value for KoL
  integer, dimension(ns), intent(in) :: KoR0  !< Correct value for KoR
  real, dimension(ns),    intent(in) :: pL0   !< Correct value for pL
  real, dimension(ns),    intent(in) :: pR0   !< Correct value for pR
  real, dimension(ns-1),    intent(in) :: hEff0 !< Correct value for hEff
  character(len=*),           intent(in) :: title !< Title for messages

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
    stdunit = 6
    if (test_nsp) stdunit = 0 ! In case of wrong results, write to error stream
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
  real,             intent(in) :: expected_pos
  real,             intent(in) :: test_pos
  character(len=*), intent(in) :: title
  ! Local variables
  integer :: stdunit = 6 ! Output to standard error
  test_rnp = expected_pos /= test_pos
  if (test_rnp) then
    write(stdunit,'(A, f20.16, " .neq. ", f20.16, " <-- WRONG")') title, expected_pos, test_pos
  else
    write(stdunit,'(A, f20.16, " .eq.  ", f20.16)') title, expected_pos, test_pos
  endif
end function test_rnp
!> Deallocates neutral_diffusion control structure
subroutine neutral_diffusion_end(CS)
  type(neutral_diffusion_CS), pointer :: CS

  if (associated(CS)) deallocate(CS)

end subroutine neutral_diffusion_end

end module MOM_neutral_diffusion
