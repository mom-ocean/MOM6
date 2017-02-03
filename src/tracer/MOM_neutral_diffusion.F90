!> A column-wise toolbox for implementing neutral diffusion
module MOM_neutral_diffusion

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,      only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,      only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,  only : diag_ctrl, time_type
use MOM_diag_mediator,  only : post_data, register_diag_field
use MOM_EOS,            only : EOS_type, calculate_compress, calculate_density_derivs
use MOM_error_handler,  only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_error_handler,  only : MOM_get_verbosity
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_file_parser,    only : openParameterBlock, closeParameterBlock
use MOM_grid,           only : ocean_grid_type
use MOM_tracer_registry,only : tracer_registry_type
use MOM_verticalGrid,   only : verticalGrid_type

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
  integer :: nkp1X2 ! Number of intersecting interfaces between columns = 2 * nkp1

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

  type(diag_ctrl), pointer :: diag ! structure to regulate output
  integer, allocatable, dimension(:) :: id_neutral_diff_tracer_conc_tend    ! tracer concentration tendency
  integer, allocatable, dimension(:) :: id_neutral_diff_tracer_cont_tend    ! tracer content tendency
  integer, allocatable, dimension(:) :: id_neutral_diff_tracer_cont_tend_2d ! k-summed tracer content tendency
  integer, allocatable, dimension(:) :: id_neutral_diff_tracer_trans_x_2d   ! k-summed ndiff zonal tracer transport
  integer, allocatable, dimension(:) :: id_neutral_diff_tracer_trans_y_2d   ! k-summed ndiff merid tracer transport

  real    :: C_p ! heat capacity of seawater (J kg-1 K-1)

end type neutral_diffusion_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40)  :: mod = "MOM_neutral_diffusion" ! module name

logical, parameter :: debug_this_module = .false.

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

  if (associated(CS)) then
    call MOM_error(FATAL, "neutral_diffusion_init called with associated control structure.")
    return
  endif

  ! Log this module and master switch for turning it on/off
  call log_version(param_file, mod, version, &
       "This module implements neutral diffusion of tracers")
  call get_param(param_file, mod, "USE_NEUTRAL_DIFFUSION", neutral_diffusion_init, &
                 "If true, enables the neutral diffusion module.", &
                 default=.false.)

  if (.not.neutral_diffusion_init) then
    return
  endif

  allocate(CS)
  CS%diag => diag


  ! Read all relevant parameters and write them to the model log.
! call openParameterBlock(param_file,'NEUTRAL_DIFF')
! call get_param(param_file, mod, "KHTR", CS%KhTr, &
!                "The background along-isopycnal tracer diffusivity.", &
!                units="m2 s-1", default=0.0)
! call closeParameterBlock(param_file)

  ! U-points
  allocate(CS%uPoL(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%uPoL(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0.
  allocate(CS%uPoR(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%uPoR(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0.
  allocate(CS%uKoL(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%uKoL(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0
  allocate(CS%uKoR(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%uKoR(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0
  allocate(CS%uHeff(G%isd:G%ied,G%jsd:G%jed,2*G%ke+1)); CS%uHeff(G%isc-1:G%iec,G%jsc:G%jec,:) = 0
  ! V-points
  allocate(CS%vPoL(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%vPoL(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0.
  allocate(CS%vPoR(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%vPoR(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0.
  allocate(CS%vKoL(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%vKoL(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0
  allocate(CS%vKoR(G%isd:G%ied,G%jsd:G%jed,2*G%ke+2)); CS%vKoR(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0
  allocate(CS%vHeff(G%isd:G%ied,G%jsd:G%jed,2*G%ke+1)); CS%vHeff(G%isc:G%iec,G%jsc-1:G%jec,:) = 0

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
       'Degree C per second')

      CS%id_neutral_diff_tracer_cont_tend(n) = register_diag_field('ocean_model',                                      &
      'ndiff_tracer_cont_tendency_'//trim(Reg%Tr(n)%name), diag%axesTL, Time,                                          &
      'Neutral diffusion tracer content tendency for '//trim(Reg%Tr(n)%name),                                          &
      'Watts/m2',cmor_field_name='opottemppmdiff', cmor_units='W m-2',                                                 &
      cmor_standard_name=                                                                                              &
      'tendency_of_sea_water_potential_temperature_expressed_as_heat_content_due_to_parameterized_mesocale_diffusion', &
      cmor_long_name =                                                                                                 &
      'Tendency of sea water potential temperature expressed as heat content due to parameterized mesocale diffusion', &
      v_extensive=.true.)

      CS%id_neutral_diff_tracer_cont_tend_2d(n) = register_diag_field('ocean_model',                                                   &
      'ndiff_tracer_cont_tendency_2d_'//trim(Reg%Tr(n)%name), diag%axesT1, Time,                                                       &
      'Depth integrated neutral diffusion tracer content tendency for '//trim(Reg%Tr(n)%name),                                         &
      'Watts/m2',cmor_field_name='opottemppmdiff_2d', cmor_units='W m-2',                                                              &
      cmor_standard_name=                                                                                                              &
      'tendency_of_sea_water_potential_temperature_expressed_as_heat_content_due_to_parameterized_mesocale_diffusion_depth_integrated',&
      cmor_long_name =                                                                                                                 &
      'Tendency of sea water potential temperature expressed as heat content due to parameterized mesocale diffusion depth integrated')

      CS%id_neutral_diff_tracer_trans_x_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_x_2d_'//trim(Reg%Tr(n)%name), diag%axesCu1, Time,                  &
      'Depth integrated neutral diffusion zonal tracer transport for '//trim(Reg%Tr(n)%name),&
      'Watts')

      CS%id_neutral_diff_tracer_trans_y_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_y_2d_'//trim(Reg%Tr(n)%name), diag%axesCv1, Time,                  &
      'Depth integrated neutral diffusion merid tracer transport for '//trim(Reg%Tr(n)%name),&
      'Watts')

    elseif(trim(Reg%Tr(n)%name) == 'S') then

      CS%id_neutral_diff_tracer_conc_tend(n) = register_diag_field('ocean_model',  &
      'ndiff_tracer_conc_tendency_'//trim(Reg%Tr(n)%name), diag%axesTL, Time,      &
      'Neutral diffusion tracer concentration tendency for '//trim(Reg%Tr(n)%name),&
       'tracer concentration units per second')

      CS%id_neutral_diff_tracer_cont_tend(n) = register_diag_field('ocean_model',                         &
      'ndiff_tracer_cont_tendency_'//trim(Reg%Tr(n)%name), diag%axesTL, Time,                             &
      'Neutral diffusion tracer content tendency for '//trim(Reg%Tr(n)%name),                             &
      'kg m-2 s-1',cmor_field_name='osaltpmdiff', cmor_units='kg m-2 s-1',                                &
      cmor_standard_name=                                                                                 &
      'tendency_of_sea_water_salinity_expressed_as_salt_content_due_to_parameterized_mesocale_diffusion', &
      cmor_long_name =                                                                                    &
      'Tendency of sea water salinity expressed as salt content due to parameterized mesocale diffusion', &
      v_extensive=.true.)

      CS%id_neutral_diff_tracer_cont_tend_2d(n) = register_diag_field('ocean_model',                                      &
      'ndiff_tracer_cont_tendency_2d_'//trim(Reg%Tr(n)%name), diag%axesT1, Time,                                          &
      'Depth integrated neutral diffusion tracer content tendency for '//trim(Reg%Tr(n)%name),                            &
      'kg m-2 s-1',cmor_field_name='osaltpmdiff_2d', cmor_units='kg m-2 s-1',                                             &
      cmor_standard_name=                                                                                                 &
      'tendency_of_sea_water_salinity_expressed_as_salt_content_due_to_parameterized_mesocale_diffusion_depth_integrated',&
      cmor_long_name =                                                                                                    &
      'Tendency of sea water salinity expressed as salt content due to parameterized mesocale diffusion depth integrated')

      CS%id_neutral_diff_tracer_trans_x_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_x_2d_'//trim(Reg%Tr(n)%name), diag%axesCu1, Time,                  &
      'Depth integrated neutral diffusion zonal tracer transport for '//trim(Reg%Tr(n)%name),&
      'kg/s')

      CS%id_neutral_diff_tracer_trans_y_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_y_2d_'//trim(Reg%Tr(n)%name), diag%axesCv1, Time,                  &
      'Depth integrated neutral diffusion merid tracer transport for '//trim(Reg%Tr(n)%name),&
      'kg/s')

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
      'kg/s')

      CS%id_neutral_diff_tracer_trans_y_2d(n) = register_diag_field('ocean_model',           &
      'ndiff_tracer_trans_y_2d_'//trim(Reg%Tr(n)%name), diag%axesCv1, Time,                  &
      'Depth integrated neutral diffusion merid tracer transport for '//trim(Reg%Tr(n)%name),&
      'kg/s')

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
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: Tint ! Interface T (degC)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: Sint ! Interface S (ppt)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: Pint ! Interface pressure (Pa)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: dRdT ! Interface thermal expansion coefficient (kg/m3/degC)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: dRdS ! Interface haline expansion coefficient (kg/m3/ppt)

  do j = G%jsc-1, G%jec+1
    ! Interpolate state to interface
    do i = G%isc-1, G%iec+1
      call interface_scalar(G%ke, h(i,j,:), T(i,j,:), Tint(i,j,:), 2)
      call interface_scalar(G%ke, h(i,j,:), S(i,j,:), Sint(i,j,:), 2)
    enddo

    ! Calculate interface properties
    Pint(:,j,1) = 0. ! Assume P=0 (Pa) at surface - needs correcting for atmospheric and ice loading - AJA
    do k = 1, G%ke+1
      call calculate_density_derivs(Tint(:,j,k), Sint(:,j,k), Pint(:,j,k), &
                                    dRdT(:,j,k), dRdS(:,j,k), G%isc-1, G%iec-G%isc+3, EOS)
      if (k<=G%ke) Pint(:,j,k+1) = Pint(:,j,k) + h(:,j,k) * GV%H_to_Pa ! Pressure at next interface, k+1 (Pa)
    enddo
  enddo

  ! Neutral surface factors at U points
  do j = G%jsc, G%jec
    do I = G%isc-1, G%iec
      call find_neutral_surface_positions(G%ke,                                          &
               Pint(i,j,:), Tint(i,j,:), Sint(i,j,:), dRdT(i,j,:), dRdS(i,j,:),          &
               Pint(i+1,j,:), Tint(i+1,j,:), Sint(i+1,j,:), dRdT(i+1,j,:), dRdS(i+1,j,:),&
               CS%uPoL(I,j,:), CS%uPoR(I,j,:), CS%uKoL(I,j,:), CS%uKoR(I,j,:), CS%uhEff(I,j,:) )
    enddo
  enddo

  ! Neutral surface factors at V points
  do J = G%jsc-1, G%jec
    do i = G%isc, G%iec
      call find_neutral_surface_positions(G%ke,                                          &
               Pint(i,j,:), Tint(i,j,:), Sint(i,j,:), dRdT(i,j,:), dRdS(i,j,:),          &
               Pint(i,j+1,:), Tint(i,j+1,:), Sint(i,j+1,:), dRdT(i,j+1,:), dRdS(i,j+1,:),&
               CS%vPoL(i,J,:), CS%vPoR(i,J,:), CS%vKoL(i,J,:), CS%vKoR(i,J,:), CS%vhEff(i,J,:) )
    enddo
  enddo

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
  real, dimension(SZIB_(G),SZJ_(G),2*G%ke+1) :: uFlx        ! Zonal flux of tracer      (concentration * H)
  real, dimension(SZI_(G),SZJB_(G),2*G%ke+1) :: vFlx        ! Meridional flux of tracer (concentration * H)
  real, dimension(SZI_(G),SZJ_(G),G%ke)      :: tendency    ! tendency array for diagn
  real, dimension(SZI_(G),SZJ_(G))           :: tendency_2d ! depth integrated content tendency for diagn
  real, dimension(SZIB_(G),SZJ_(G))          :: trans_x_2d  ! depth integrated diffusive tracer x-transport diagn
  real, dimension(SZI_(G),SZJB_(G))          :: trans_y_2d  ! depth integrated diffusive tracer y-transport diagn
  real, dimension(G%ke)                      :: dTracer     ! change in tracer concentration due to ndiffusion
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


  ! x-flux
  do j = G%jsc,G%jec ; do I = G%isc-1,G%iec
    if (G%mask2dCu(I,j)>0.) then
      call neutral_surface_flux(nk, h(i,j,:), h(i+1,j,:),       &
                                Tracer(i,j,:), Tracer(i+1,j,:), &
                                CS%uPoL(I,j,:), CS%uPoR(I,j,:), &
                                CS%uKoL(I,j,:), CS%uKoR(I,j,:), &
                                CS%uhEff(I,j,:), uFlx(I,j,:))
    else
      uFlx(I,j,:) = 0.
    endif
  enddo ; enddo

  ! y-flux
  do J = G%jsc-1,G%jec ; do i = G%isc,G%iec
    if (G%mask2dCv(i,J)>0.) then
      call neutral_surface_flux(nk, h(i,j,:), h(i,j+1,:),       &
                                Tracer(i,j,:), Tracer(i,j+1,:), &
                                CS%vPoL(i,J,:), CS%vPoR(i,J,:), &
                                CS%vKoL(i,J,:), CS%vKoR(i,J,:), &
                                CS%vhEff(i,J,:), vFlx(i,J,:))
    else
      vFlx(I,j,:) = 0.
    endif
  enddo ; enddo

  ! Update the tracer concentration from divergence of neutral diffusive flux components
  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    if (G%mask2dT(i,j)>0.) then

      dTracer(:) = 0.
      do ks = 1,2*nk+1 ;
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
        do ks = 1,2*nk+1 ;
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
        do ks = 1,2*nk+1 ;
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


!> Returns positions within left/right columns of combined interfaces
subroutine find_neutral_surface_positions(nk, Pl, Tl, Sl, dRdTl, dRdSl, Pr, Tr, Sr, dRdTr, dRdSr, PoL, PoR, KoL, KoR, hEff)
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

  ! Initialize variables for the search
  kr = 1 ; lastK_right = 1 ; lastP_right = 0.
  kl = 1 ; lastK_left = 1 ; lastP_left = 0.
  reached_bottom = .false.

  ! Loop over each neutral surface, working from top to bottom
  neutral_surfaces: do k_surface = 1, 2*nk+2
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
      hL = absolute_position(nk,Pl,KoL,PoL,k_surface) - absolute_position(nk,Pl,KoL,PoL,k_surface-1)
      hR = absolute_position(nk,Pr,KoR,PoR,k_surface) - absolute_position(nk,Pr,KoR,PoR,k_surface-1)
      if ( hL + hR > 0.) then
        hEff(k_surface-1) = 2. * hL * hR / ( hL + hR )
      else
        hEff(k_surface-1) = 0.
      endif
    endif

  enddo neutral_surfaces

end subroutine find_neutral_surface_positions

!> Converts non-dimensional position within a layer to absolute position (for debugging)
real function absolute_position(n,Pint,Karr,NParr,k_surface)
  integer, intent(in) :: n            !< Number of levels
  real,    intent(in) :: Pint(n+1)    !< Position of interfaces (Pa)
  integer, intent(in) :: Karr(2*n+2)  !< Index of interface above position
  real,    intent(in) :: NParr(2*n+2) !< Non-dimensional position within layer Karr(:)

  ! Local variables
  integer :: k_surface, k

  k = Karr(k_surface)
  if (k>n) stop 'absolute_position: k>nk is out of bounds!'
  absolute_position = Pint(k) + NParr(k_surface) * ( Pint(k+1) - Pint(k) )

end function absolute_position

!> Converts non-dimensional positions within layers to absolute positions (for debugging)
function absolute_positions(n,Pint,Karr,NParr)
  integer, intent(in) :: n            !< Number of levels
  real,    intent(in) :: Pint(n+1)    !< Position of interface (Pa)
  integer, intent(in) :: Karr(2*n+2)  !< Indexes of interfaces about positions
  real,    intent(in) :: NParr(2*n+2) !< Non-dimensional positions within layers Karr(:)

  real,  dimension(2*n+2) :: absolute_positions ! Absolute positions (Pa)

  ! Local variables
  integer :: k_surface, k

  do k_surface = 1, 2*n+2
    absolute_positions(k_surface) = absolute_position(n,Pint,Karr,NParr,k_surface)
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


!> Returns a single column of neutral diffusion fluxes of a tracer.
subroutine neutral_surface_flux(nk, hl, hr, Tl, Tr, PiL, PiR, KoL, KoR, hEff, Flx)
  integer,                    intent(in)    :: nk    !< Number of levels
  real, dimension(nk),        intent(in)    :: hl    !< Left-column layer thickness (Pa)
  real, dimension(nk),        intent(in)    :: hr    !< Right-column layer thickness (Pa)
  real, dimension(nk),        intent(in)    :: Tl    !< Left-column layer tracer (conc, e.g. degC)
  real, dimension(nk),        intent(in)    :: Tr    !< Right-column layer tracer (conc, e.g. degC)
  real, dimension(2*nk+2),    intent(in)    :: PiL   !< Fractional position of neutral surface within layer KoL of left column
  real, dimension(2*nk+2),    intent(in)    :: PiR   !< Fractional position of neutral surface within layer KoR of right column
  integer, dimension(2*nk+2), intent(in)    :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(2*nk+2), intent(in)    :: KoR   !< Index of first right interface above neutral surface
  real, dimension(2*nk+1),    intent(in)    :: hEff  !< Effective thickness between two neutral surfaces (Pa)
  real, dimension(2*nk+1),    intent(inout) :: Flx   !< Flux of tracer between pairs of neutral layers (conc H)

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

  call interface_scalar(nk, hl, Tl, Til, 2)
  call interface_scalar(nk, hr, Tr, Tir, 2)

  ! Setup reconstruction edge values
  do k = 1, nk
    aL_l(k) = Til(k)
    aR_l(k) = Til(k+1)
    if ( signum(1., aR_l(k) - Tl(k))*signum(1., Tl(k) - aL_l(k)) <= 0.0 ) then
      aL_l(k) = Tl(k)
      aR_l(k) = Tl(k)
    elseif ( sign(3., aR_l(k) - aL_l(k)) * ( (Tl(k) - aL_l(k)) + (Tl(k) - aR_l(k))) > abs(aR_l(k) - aL_l(k)) ) then
      aL_l(k) = Tl(k) + 2.0 * ( Tl(k) - aR_l(k) )
    elseif ( sign(3., aR_l(k) - aL_l(k)) * ( (Tl(k) - aL_l(k)) + (Tl(k) - aR_l(k))) < -abs(aR_l(k) - aL_l(k)) ) then
      aR_l(k) = Tl(k) + 2.0 * ( Tl(k) - aL_l(k) )
    endif
    aL_r(k) = Tir(k)
    aR_r(k) = Tir(k+1)
    if ( signum(1., aR_r(k) - Tr(k))*signum(1., Tr(k) - aL_r(k)) <= 0.0 ) then
      aL_r(k) = Tr(k)
      aR_r(k) = Tr(k)
    elseif ( sign(3., aR_r(k) - aL_r(k)) * ( (Tr(k) - aL_r(k)) + (Tr(k) - aR_r(k))) > abs(aR_r(k) - aL_r(k)) ) then
      aL_r(k) = Tr(k) + 2.0 * ( Tr(k) - aR_r(k) )
    elseif ( sign(3., aR_r(k) - aL_r(k)) * ( (Tr(k) - aL_r(k)) + (Tr(k) - aR_r(k))) < -abs(aR_r(k) - aL_r(k)) ) then
      aR_r(k) = Tr(k) + 2.0 * ( Tr(k) - aL_r(k) )
    endif
  enddo

  do k_sublayer = 1, 2*nk+1
    if (hEff(k_sublayer) == 0.) then
      Flx(k_sublayer) = 0.
    else

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
      if (signum(1.,dT_top) * signum(1.,dT_bottom) <= 0. .or. signum(1.,dT_ave) * signum(1.,dT_layer) <= 0. ) then
        dT_ave = 0.
      else
        dT_ave = dT_layer
      endif
      Flx(k_sublayer) = dT_ave * hEff(k_sublayer)
    endif
  enddo

end subroutine neutral_surface_flux

!> Returns true if unit tests of neutral_diffusion functions fail. Otherwise returns false.
logical function neutral_diffusion_unit_tests()
  integer, parameter         :: nk = 4
  real, dimension(nk+1)      :: TiL, TiR1, TiR2, TiR4, Tio ! Test interface temperatures
  real, dimension(nk)        :: TL                         ! Test layer temperatures
  real, dimension(nk+1)      :: SiL                        ! Test interface salinities
  real, dimension(nk+1)      :: PiL, PiR4                  ! Test interface positions
  real, dimension(2*nk+2)    :: PiLRo, PiRLo               ! Test positions
  integer, dimension(2*nk+2) :: KoL, KoR                   ! Test indexes
  real, dimension(2*nk+1)    :: hEff                       ! Test positions
  real, dimension(2*nk+1)    :: Flx                        ! Test flux

  integer :: k, verbosity

  verbosity = MOM_get_verbosity()

  neutral_diffusion_unit_tests = .false. ! Normally return false
  write(*,'(a)') '===== MOM_neutral_diffusion: neutral_diffusion_unit_tests ==============='

  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fv_diff(1.,1.,1., 0.,1.,2., 1., 'FV: Straight line on uniform grid')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fv_diff(1.,1.,0., 0.,4.,8., 7., 'FV: Vanished right cell')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fv_diff(0.,1.,1., 0.,4.,8., 7., 'FV: Vanished left cell')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fv_diff(1.,2.,4., 0.,3.,9., 4., 'FV: Stretched grid')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fv_diff(2.,0.,2., 0.,1.,2., 0., 'FV: Vanished middle cell')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fv_diff(0.,1.,0., 0.,1.,2., 2., 'FV: Vanished on both sides')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fv_diff(1.,0.,0., 0.,1.,2., 0., 'FV: Two vanished cell sides')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fv_diff(0.,0.,0., 0.,1.,2., 0., 'FV: All vanished cells')

  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fvlsq_slope(1.,1.,1., 0.,1.,2., 1., 'LSQ: Straight line on uniform grid')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fvlsq_slope(1.,1.,0., 0.,1.,2., 1., 'LSQ: Vanished right cell')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fvlsq_slope(0.,1.,1., 0.,1.,2., 1., 'LSQ: Vanished left cell')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fvlsq_slope(1.,2.,4., 0.,3.,9., 2., 'LSQ: Stretched grid')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fvlsq_slope(1.,0.,1., 0.,1.,2., 2., 'LSQ: Vanished middle cell')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fvlsq_slope(0.,1.,0., 0.,1.,2., 0., 'LSQ: Vanished on both sides')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fvlsq_slope(1.,0.,0., 0.,1.,2., 0., 'LSQ: Two vanished cell sides')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_fvlsq_slope(0.,0.,0., 0.,1.,2., 0., 'LSQ: All vanished cells')

  call interface_scalar(4, (/10.,10.,10.,10./), (/24.,18.,12.,6./), Tio, 1)
  !neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_data1d(5, Tio, (/27.,21.,15.,9.,3./), 'Linear profile, interface temperatures')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_data1d(5, Tio, (/24.,22.5,15.,7.5,6./), 'Linear profile, linear interface temperatures')
  call interface_scalar(4, (/10.,10.,10.,10./), (/24.,18.,12.,6./), Tio, 2)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_data1d(5, Tio, (/24.,22.,15.,8.,6./), 'Linear profile, PPM interface temperatures')

  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp(-1.0, 0.,  1.0, 1.0, 0.5, 'Check mid-point')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp( 0.0, 0.,  1.0, 1.0, 0.0, 'Check bottom')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp( 0.1, 0.,  1.1, 1.0, 0.0, 'Check below')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp(-1.0, 0.,  0.0, 1.0, 1.0, 'Check top')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp(-1.0, 0., -0.1, 1.0, 1.0, 'Check above')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp(-1.0, 0.,  3.0, 1.0, 0.25, 'Check 1/4')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp(-3.0, 0.,  1.0, 1.0, 0.75, 'Check 3/4')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp( 1.0, 0.,  1.0, 1.0, 0.0, 'Check dRho=0 below')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp(-1.0, 0., -1.0, 1.0, 1.0, 'Check dRho=0 above')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp( 0.0, 0.,  0.0, 1.0, 0.5, 'Check dRho=0 mid')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_ifndp(-2.0, .5,  5.0, 0.5, 0.5, 'Check dP=0')

  ! Identical columns
  call find_neutral_surface_positions(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or.  test_nsp(3, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,2,2,3,3,3,3/), & ! KoL
                                   (/1,1,2,2,3,3,3,3/), & ! KoR
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pR
                                   (/0.,10.,0.,10.,0.,10.,0./), & ! hEff
                                   'Identical columns')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_data1d(8, &
                                   absolute_positions(3, (/0.,10.,20.,30./), KoL, PiLRo), &
                                   (/0.,0.,10.,10.,20.,20.,30.,30./), '... left positions')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_data1d(8, &
                                   absolute_positions(3, (/0.,10.,20.,30./), KoR, PiRLo), &
                                   (/0.,0.,10.,10.,20.,20.,30.,30./), '... right positions')
  call neutral_surface_flux(3, (/10.,10.,10./), (/10.,10.,10./), & ! nk, hL, hR
                               (/20.,16.,12./), (/20.,16.,12./), & ! Tl, Tr
                               PiLRo, PiRLo, KoL, KoR, hEff, Flx)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_data1d(7, Flx, &
              (/0.,0.,0.,0.,0.,0.,0./), 'Identical columns, rho flux (=0)')
  call neutral_surface_flux(3, (/10.,10.,10./), (/10.,10.,10./), & ! nk, hL, hR
                               (/-1.,-1.,-1./), (/1.,1.,1./), & ! Sl, Sr
                               PiLRo, PiRLo, KoL, KoR, hEff, Flx)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_data1d(7, Flx, &
              (/0.,20.,0.,20.,0.,20.,0./), 'Identical columns, S flux')

  ! Right column slightly cooler than left
  call find_neutral_surface_positions(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/20.,16.,12.,8./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or.  test_nsp(3, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,2,2,3,3,3,3/), & ! kL
                                   (/1,1,1,2,2,3,3,3/), & ! kR
                                   (/0.,0.5,0.,0.5,0.,0.5,1.,1./), & ! pL
                                   (/0.,0.,0.5,0.,0.5,0.,0.5,1./), & ! pR
                                   (/0.,5.,5.,5.,5.,5.,0./), & ! hEff
                                   'Right column slightly cooler')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_data1d(8, &
                                   absolute_positions(3, (/0.,10.,20.,30./), KoL, PiLRo), &
                                   (/0.,5.,10.,15.,20.,25.,30.,30./), '... left positions')
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or. test_data1d(8, &
                                   absolute_positions(3, (/0.,10.,20.,30./), KoR, PiRLo), &
                                   (/0.,0.,5.,10.,15.,20.,25.,30./), '... right positions')

  ! Right column slightly warmer than left
  call find_neutral_surface_positions(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/24.,20.,16.,12./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or.  test_nsp(3, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,1,2,2,3,3,3/), & ! kL
                                   (/1,1,2,2,3,3,3,3/), & ! kR
                                   (/0.,0.,0.5,0.,0.5,0.,0.5,1./), & ! pL
                                   (/0.,0.5,0.,0.5,0.,0.5,1.,1./), & ! pR
                                   (/0.,5.,5.,5.,5.,5.,0./), & ! hEff
                                   'Right column slightly warmer')

  ! Right column somewhat cooler than left
  call find_neutral_surface_positions(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/16.,12.,8.,4./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or.  test_nsp(3, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,2,3,3,3,3,3/), & ! kL
                                   (/1,1,1,1,2,2,3,3/), & ! kR
                                   (/0.,0.,0.5,0.,0.5,1.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.5,0.,0.5,0.,1./), & ! pR
                                   (/0.,0.,5.,5.,5.,0.,0./), & ! hEff
                                   'Right column somewhat cooler')

  ! Right column much colder than left with no overlap
  call find_neutral_surface_positions(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/9.,7.,5.,3./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or.  test_nsp(3, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,3,3,3,3,3,3/), & ! kL
                                   (/1,1,1,1,1,2,3,3/), & ! kR
                                   (/0.,0.,0.,1.,1.,1.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,0.,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,0.,0./), & ! hEff
                                   'Right column much cooler')

  ! Right column with mixed layer
  call find_neutral_surface_positions(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or.  test_nsp(3, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,3,3,3,3,3,3/), & ! kL
                                   (/1,1,1,1,2,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,1.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,0.,1./), & ! pR
                                   (/0.,0.,0.,0.,10.,0.,0./), & ! hEff
                                   'Right column with mixed layer')

  ! Identical columns with mixed layer
  call find_neutral_surface_positions(3, &
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or.  test_nsp(3, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,2,2,3,3,3,3/), & ! kL
                                   (/1,1,2,2,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pR
                                   (/0.,10.,0.,10.,0.,10.,0./), & ! hEff
                                   'Indentical columns with mixed layer')

  ! Right column with unstable mixed layer
  call find_neutral_surface_positions(3, &
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/10.,14.,12.,4./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or.  test_nsp(3, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,3,3,3,3,3,3/), & ! kL
                                   (/1,1,1,2,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.,.75,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.25,1.,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,7.5,0./), & ! hEff
                                   'Right column with unstable mixed layer')

  ! Left column with unstable mixed layer
  call find_neutral_surface_positions(3, &
             (/0.,10.,20.,30./), (/10.,14.,12.,4./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or.  test_nsp(3, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,1,2,3,3,3,3/), & ! kL
                                   (/1,2,3,3,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.25,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,.75,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,7.5,0./), & ! hEff
                                   'Left column with unstable mixed layer')

  ! Two unstable mixed layers
  call find_neutral_surface_positions(3, &
             (/0.,10.,20.,30./), (/8.,12.,10.,2./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/10.,14.,12.,4./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  neutral_diffusion_unit_tests = neutral_diffusion_unit_tests .or.  test_nsp(3, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,1,1,2,3,3,3/), & ! kL
                                   (/1,2,3,3,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.,0.75,1./), & ! pL
                                   (/0.,0.,0.,0.5,0.5,0.5,1.,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,6.,0./), & ! hEff
                                   'Two unstable mixed layers')

  write(*,'(a)') '=========================================================='

  contains

  !> Returns true if a test of fv_diff() fails, and conditionally writes results to stream
  logical function test_fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
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

    if (test_fv_diff .or. verbosity>5) then
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
  logical function test_fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
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

    if (test_fvlsq_slope .or. verbosity>5) then
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
  logical function test_ifndp(rhoNeg, Pneg, rhoPos, Ppos, Ptrue, title)
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

    if (test_ifndp .or. verbosity>5) then
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
  logical function test_data1d(nk, Po, Ptrue, title)
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

    if (test_data1d .or. verbosity>5) then
      stdunit = 6
      if (test_data1d) stdunit = 0 ! In case of wrong results, write to error stream
      write(stdunit,'(a)') title
      do k = 1,nk
        if (Po(k) /= Ptrue(k)) then
          test_data1d = .true.
          write(stdunit,'(a,i2,2(x,a,f20.16),x,a,1pe22.15,x,a)') 'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k),'WRONG!'
        else
          if (verbosity>5) &
            write(stdunit,'(a,i2,2(x,a,f20.16),x,a,1pe22.15)') 'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k)
        endif
      enddo
    endif

  end function test_data1d

  !> Returns true if comparison of Po and Ptrue fails, and conditionally writes results to stream
  logical function test_data1di(nk, Po, Ptrue, title)
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

    if (test_data1di .or. verbosity>5) then
      stdunit = 6
      if (test_data1di) stdunit = 0 ! In case of wrong results, write to error stream
      write(stdunit,'(a)') title
      do k = 1,nk
        if (Po(k) /= Ptrue(k)) then
          test_data1di = .true.
          write(stdunit,'(a,i2,2(x,a,i5),x,a)') 'k=',k,'Io=',Po(k),'Itrue=',Ptrue(k),'WRONG!'
        else
          if (verbosity>5) &
            write(stdunit,'(a,i2,2(x,a,i5))') 'k=',k,'Io=',Po(k),'Itrue=',Ptrue(k)
        endif
      enddo
    endif

  end function test_data1di

  !> Returns true if output of find_neutral_surface_positions() does not match correct values, and conditionally writes results to stream
  logical function test_nsp(nk, KoL, KoR, pL, pR, hEff, KoL0, KoR0, pL0, pR0, hEff0, title)
    integer,                    intent(in) :: nk    !< Number of layers
    integer, dimension(2*nk+2), intent(in) :: KoL   !< Index of first left interface above neutral surface
    integer, dimension(2*nk+2), intent(in) :: KoR   !< Index of first right interface above neutral surface
    real, dimension(2*nk+2),    intent(in) :: pL    !< Fractional position of neutral surface within layer KoL of left column
    real, dimension(2*nk+2),    intent(in) :: pR    !< Fractional position of neutral surface within layer KoR of right column
    real, dimension(2*nk+1),    intent(in) :: hEff  !< Effective thickness between two neutral surfaces (Pa)
    integer, dimension(2*nk+2), intent(in) :: KoL0  !< Correct value for KoL
    integer, dimension(2*nk+2), intent(in) :: KoR0  !< Correct value for KoR
    real, dimension(2*nk+2),    intent(in) :: pL0   !< Correct value for pL
    real, dimension(2*nk+2),    intent(in) :: pR0   !< Correct value for pR
    real, dimension(2*nk+1),    intent(in) :: hEff0 !< Correct value for hEff
    character(len=*),           intent(in) :: title !< Title for messages

    ! Local variables
    integer :: k, stdunit
    logical :: this_row_failed

    test_nsp = .false.
    do k = 1,2*nk+2
      test_nsp = test_nsp .or. compare_nsp_row(KoL(k), KoR(k), pL(k), pR(k), KoL0(k), KoR0(k), pL0(k), pR0(k))
      if (k < 2*nk+2) then
        if (hEff(k) /= hEff0(k)) test_nsp = .true.
      endif
    enddo

    if (test_nsp .or. verbosity>5) then
      stdunit = 6
      if (test_nsp) stdunit = 0 ! In case of wrong results, write to error stream
      write(stdunit,'(a)') title
      do k = 1,2*nk+2
        this_row_failed = compare_nsp_row(KoL(k), KoR(k), pL(k), pR(k), KoL0(k), KoR0(k), pL0(k), pR0(k))
        if (this_row_failed) then
          write(stdunit,10) k,KoL(k),pL(k),KoR(k),pR(k),' <-- WRONG!'
          write(stdunit,10) k,KoL0(k),pL0(k),KoR0(k),pR0(k),' <-- should be this'
        else
          write(stdunit,10) k,KoL(k),pL(k),KoR(k),pR(k)
        endif
        if (k < 2*nk+2) then
          if (hEff(k) /= hEff0(k)) then
            write(stdunit,'(i3,8x,"layer hEff =",2(f20.16,a))') k,hEff(k)," .neq. ",hEff0(k),' <-- WRONG!'
          else
            write(stdunit,'(i3,8x,"layer hEff =",f20.16)') k,hEff(k)
          endif
        endif
      enddo
    endif

10  format("ks=",i3," kL=",i3," pL=",f20.16," kR=",i3," pR=",f20.16,a)
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

end function neutral_diffusion_unit_tests

!> Deallocates neutral_diffusion control structure
subroutine neutral_diffusion_end(CS)
  type(neutral_diffusion_CS), pointer :: CS

  if (associated(CS)) deallocate(CS)

end subroutine neutral_diffusion_end

end module MOM_neutral_diffusion
