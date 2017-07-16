module MOM_diffConvection

use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : query_averaging_enabled, register_diag_field
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_PE
use MOM_EOS, only : EOS_type, calculate_density
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_file_parser, only : openParameterBlock, closeParameterBlock
use MOM_grid, only : ocean_grid_type, isPointInCell
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include "MOM_memory.h"

public :: diffConvection_init, diffConvection_calculate, diffConvection_end

! Control structure for containing KPP parameters/data
type, public :: diffConvection_CS ; private

  ! Parameters
  real    :: Kd_convection ! The value of diffusivity to add at statically unstable interfaces (m2/s)
  logical :: debug         ! If true, turn on debugging
  logical :: passiveMode   ! If true, make the motions but go nowhere

  ! Daignostic handles and pointers
  type(diag_ctrl), pointer :: diag => NULL()
  integer :: id_N2 = -1, id_Kd_conv = -1

  ! Diagnostics arrays
  real, allocatable, dimension(:,:,:) :: N2      ! Brunt-Vaisala frequency (1/s2)
  real, allocatable, dimension(:,:,:) :: Kd_conv ! Diffusivity added by convection (m2/s)

end type diffConvection_CS

! Module data used for debugging only
logical, parameter :: verbose = .False.

contains

logical function diffConvection_init(paramFile, G, diag, Time, CS)
! Initialize the CVmix KPP module and set up diagnostics
! Returns True if module is to be used, otherwise returns False.

! Arguments
  type(param_file_type),   intent(in)    :: paramFile ! File parser
  type(ocean_grid_type),   intent(in)    :: G         ! Ocean grid
  type(diag_ctrl), target, intent(in)    :: diag      ! Diagnostics
  type(time_type),         intent(in)    :: Time      ! Time
  type(diffConvection_CS), pointer       :: CS        ! Control structure
! Local variables
#include "version_variable.h"
  character(len=40) :: mod = 'MOM_diffConvection' ! This module's name.

  if (associated(CS)) call MOM_error(FATAL, 'MOM_diffConvection, diffConvection_init: '// &
           'Control structure has already been initialized')
  allocate(CS)

! Read parameters
  call log_version(paramFile, mod, version, &
            'This module implements enhanced diffusivity as a\n' // &
            'function of static stability, N^2.')
  call get_param(paramFile, mod, "USE_CONVECTION", diffConvection_init, &
                 "If true, turns on the diffusive convection scheme that\n"// &
                 "increases diapycnal diffusivities at statically unstable\n"// &
                 "interfaces. Relevant parameters are contained in the\n"// &
                 "CONVECTION% parameter block.", &
                 default=.false.)

  call openParameterBlock(paramFile,'CONVECTION')
  call get_param(paramFile, mod, 'PASSIVE', CS%passiveMode,           &
                 'If True, puts KPP into a passive-diagnostic mode.', &
                 default=.False.)
  call get_param(paramFile, mod, 'KD_CONV', CS%Kd_convection,                  &
                 'DIffusivity used in statically unstable regions of column.', &
                 units='m2/s', default=1.00)
  call closeParameterBlock(paramFile)
  call get_param(paramFile, mod, 'DEBUG', CS%debug, default=.False., do_not_log=.True.)

! Forego remainder of initialization if not using this scheme
  if (.not. diffConvection_init) return

! Register diagnostics
  CS%diag => diag
  CS%id_N2 = register_diag_field('ocean_model', 'Conv_N2', diag%axesTi, Time, &
      'Square of Brunt-Vaisala frequency used by diffConvection module', '1/s2')
  if (CS%id_N2 > 0) allocate( CS%N2( SZI_(G), SZJ_(G), SZK_(G)+1 ) )
  CS%id_Kd_conv = register_diag_field('ocean_model', 'Conv_Kd', diag%axesTi, Time, &
      'Additional diffusivity added by diffConvection module', 'm2/s')
  if (CS%id_Kd_conv > 0) allocate( CS%Kd_conv( SZI_(G), SZJ_(G), SZK_(G)+1 ) )

  if (CS%id_N2 > 0) CS%N2(:,:,:) = 0.
  if (CS%id_Kd_conv > 0) CS%Kd_conv(:,:,:) = 0.

end function diffConvection_init


subroutine diffConvection_calculate(CS, G, GV, h, Temp, Salt, EOS, Kd_int)
! Calculates diffusivity and non-local transport for KPP parameterization

! Arguments
  type(diffConvection_CS),                   pointer       :: CS    ! Control structure
  type(ocean_grid_type),                     intent(in)    :: G     ! Ocean grid
  type(verticalGrid_type),                   intent(in)    :: GV    ! Ocean vertical grid
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h     ! Layer/level thicknesses (units of H)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: Temp  ! Pot. temperature (degrees C)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: Salt  ! Salinity (ppt)
  type(EOS_type),                            pointer       :: EOS   ! Equation of state
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kd_int ! (in) Vertical diffusivity on interfaces (m2/s)
                                                                 ! (out) Modified vertical diffusivity (m2/s)
! Local variables
  integer :: i, j, k
  real, dimension( G%ke+1 ) :: N2_1d ! Brunt-Vaisala frequency squared, at interfaces (1/s2)
  real, dimension( G%ke+1 ) :: Kd_1d ! Vertical diffusivity at interfaces (m2/s)
  real :: GoRho, pRef, rhoK, rhoKm1

  GoRho = GV%g_Earth / GV%Rho0

  N2_1d( 1 ) = 0.
  N2_1d( G%ke+1 ) = 0.
  Kd_1d( 1 ) = 0.
  Kd_1d( G%ke+1 ) = 0.
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! This k-loop calculates external quantities independent of any iterations
    ! Start at bottom of top level
    pRef = 0. ! Ignore atmospheric pressure
    do K = 2, G%ke
      ! Pressure at interface K is incremented by mass of level above
      pRef = pRef + GV%g_Earth * GV%Rho0 * h(i,j,k-1) * GV%H_to_m ! Boussinesq approximation!!!! ?????
      ! Compute Brunt-Vaisala frequency (static stability) on interfaces
      call calculate_density(Temp(i,j,k),   Salt(i,j,k),   pRef, rhoK,   EOS)
      call calculate_density(Temp(i,j,k-1), Salt(i,j,k-1), pRef, rhoKm1, EOS)
      N2_1d(K) = GoRho * (rhoK - rhoKm1) / &
              (0.5*(h(i,j,k-1) + h(i,j,k)) + GV%H_subroundoff) ! Can be negative
      Kd_1d(K) = 0.
      if (N2_1d(K) < 0.) Kd_1d(K) = CS%Kd_convection
    enddo ! k

    if (.not. CS%passiveMode) Kd_int(i,j,:) = Kd_int(i,j,:) + Kd_1d(:)

    if (CS%id_N2 > 0) CS%N2(i,j,:) = N2_1d(:)
    if (CS%id_Kd_conv > 0) CS%Kd_conv(i,j,:) = Kd_1d(:)

  enddo ; enddo ! j

  if (CS%id_N2 > 0) call post_data(CS%id_N2, CS%N2, CS%diag)
  if (CS%id_Kd_conv > 0) call post_data(CS%id_Kd_conv, CS%Kd_conv, CS%diag)

end subroutine diffConvection_calculate


subroutine diffConvection_end(CS)
! Clear pointers, dealocate memory
  type(diffConvection_CS), pointer :: CS ! Control structure

  if (CS%id_N2 > 0) deallocate(CS%N2, CS%diag)
  if (CS%id_Kd_conv > 0) deallocate(CS%Kd_conv, CS%diag)
  deallocate(CS)
end subroutine diffConvection_end

end module MOM_diffConvection
