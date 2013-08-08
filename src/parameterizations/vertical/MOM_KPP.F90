module MOM_KPP

use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : query_averaging_enabled, register_diag_field
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_EOS, only : EOS_type, calculate_density
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_file_parser, only : openParameterBlock, closeParameterBlock
use MOM_grid, only : ocean_grid_type

use cvmix_kpp, only : cvmix_init_kpp, cvmix_put_kpp, cvmix_get_kpp_real
use cvmix_kpp, only : cvmix_coeffs_kpp, cvmix_kpp_compute_OBL_depth
use cvmix_kpp, only : cvmix_kpp_params_type

implicit none ; private

#include "MOM_memory.h"

public :: KPP_init, KPP_calculate, KPP_end

! Control structure for containing KPP parameters/data
type, public :: KPP_CS ; private

  ! Parameters
  real    :: Ri_crit    ! Critical Richardson number (defines OBL depth)
  real    :: vonKarman  ! von Karman constant
! real    :: zeta_m     ! parameter for computing vel scale func
! real    :: zeta_s     ! parameter for computing vel scale func
! real    :: a_m        ! parameter for computing vel scale func
! real    :: a_s        ! parameter for computing vel scale func
! real    :: c_m        ! parameter for computing vel scale func
! real    :: c_s        ! parameter for computing vel scale func
! real    :: eps        ! small non-negative val (rec 1e-10)
  character(len=10) :: interpType ! Type of iterpolation to use in determining OBL
  logical :: computeEkman ! If True, compute Ekman depth limit
  logical :: computeMoninObukhov ! If True, compute Monin-Obukhov limit

  ! CVmix parameters
  type(cvmix_kpp_params_type), pointer :: KPP_params => NULL()

  ! Daignostic handles and pointers
  type(diag_ctrl), pointer :: diag => NULL()
  integer :: id_OBL = -1, id_BulkRi = -1

end type KPP_CS

contains

subroutine KPP_init(paramFile, G, diag, Time, CS)
! Initialize the CVmix KPP module and set up diagnostics

! Arguments
  type(param_file_type),   intent(in)    :: paramFile ! File parser
  type(ocean_grid_type),   intent(in)    :: G         ! Ocean grid
  type(diag_ctrl), target, intent(in)    :: diag      ! Diagnostics
  type(time_type),         intent(in)    :: Time      ! Time
  type(KPP_CS),            pointer       :: CS        ! Control structure
! Local variables
#include "version_variable.h"
  character(len=40) :: mod = 'MOM_KPP' ! This module's name.

  if (associated(CS)) call MOM_error(FATAL, 'MOM_KPP, KPP_init: '// &
           'Control structure has already been initialized')
  allocate(CS)

! Read parameters
  call log_version(paramFile, mod, version, '')
  call openParameterBlock(paramFile,'KPP')
  call get_param(paramFile, mod, 'RI_CRIT', CS%Ri_crit,                       &
                 'Critical Richardson number used to define depth of the\n'// &
                 'Oceab Boundary Layer (OBL).',                               &
                 units='nondim', default=0.3)
  call get_param(paramFile, mod, 'VON_KARMAN', CS%vonKarman, &
                 'von Karman constant.',                     &
                 units='nondim', default=0.41)
  call get_param(paramFile, mod, 'INTERP_TYPE', CS%interpType,                  &
                 'Type of interpolation to use to determine the OBL depth.\n'// &
                 'Allowed types are: linear, quadratic, cubic.',                &
                 default='quadratic')
  call get_param(paramFile, mod, 'COMPUTE_EKMAN', CS%computeEkman,                     &
                 'If True, limit the OBL depth to be shallower than the Ekman depth.', &
                 default=.False.)
  call get_param(paramFile, mod, 'COMPUTE_MONIN_OBUKHOV', CS%computeMoninObukhov,  &
                 'If True, limit the OBL depth to be shallower than the\n'//       &
                 'Monin-Obukhov depth.',                                           &
                 default=.False.)
  call closeParameterBlock(paramFile)

  call cvmix_init_kpp( Ri_crit=CS%Ri_crit,                 &
                       vonKarman=CS%vonKarman,             &
                       interp_type=CS%interpType,          &
                       lEkman=CS%computeEkman,             &
                       lMonOb=CS%computeMoninObukhov,      &
                       CVmix_kpp_params_user=CS%KPP_params )

! Register diagnostics
  CS%diag => diag
  CS%id_OBL = register_diag_field('ocean_model', 'KPP_OBLdepth', diag%axesT1, Time, &
      'Thickness of the surface Ocean Boundary Layer calculated by [CVmix] KPP', 'meter')
  CS%id_BulkRi = register_diag_field('ocean_model', 'KPP_BulkRi', diag%axesTL, Time, &
      'Bulk Richardson number used to find the OBL depth used by [CVmix] KPP', 'nondim')

end subroutine KPP_init

subroutine KPP_calculate(CS, G, h, Temp, Salt, u, v, EOS, uStar, bFlux, Kv)
! Calculates diffusivity and non-local transport for KPP parameterization 

! Arguments
  type(KPP_CS),                           intent(in)    :: CS    ! Control structure
  type(ocean_grid_type),                  intent(in)    :: G     ! Ocean grid
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h     ! layer/level thicknesses (units of H)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: Temp  ! Pot. temperature (degrees C)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: Salt  ! Salinity (ppt)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)    :: u     ! Velocity components (m/s)
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)    :: v     ! Velocity components (m/s)
  type(EOS_type),                         pointer       :: EOS   ! Equation of state
  real, dimension(NIMEM_,NJMEM_),         intent(in)    :: uStar ! Piston velocity (m/s)
  real, dimension(NIMEM_,NJMEM_),         intent(in)    :: bFlux ! Buoyancy flux (m2/s3)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(inout) :: Kv ! Vertical diffusivity due to KPP

! Local variables
  real :: BulkRi( SZI_(G), SZJ_(G), SZK_(G) ) ! Bulk Richardson number for each layer
  real :: cellHeight( G%ke ) ! Cell center heights referenced to surface (m)
  real :: OBLdepth( SZI_(G), SZJ_(G) ) ! Depth (positive) of OBL (m)
  real :: OBLheight, surfFricVel, surfBuoyFlux, iFaceHeight, Coriolis
  integer :: i, j, k, kOBL

  call calculateBulkRichardson(CS, G, h, Temp, Salt, u, v, EOS, BulkRi )

  do j = G%jsc, G%jec
    do i = G%isc, G%iec

      ! Compute heights, referenced to the surface (z=0)
      iFaceHeight = 0.
      do k = 1, G%ke
        cellHeight(k) = iFaceHeight - 0.5 * h(i,j,k) * G%H_to_m ! in metres 
        iFaceHeight = iFaceHeight - h(i,j,k) * G%H_to_m ! in metres
      enddo ! k

      Coriolis = 0.25*( (G%CoriolisBu(i,j) + G%CoriolisBu(i-1,j-1)) &
                       +(G%CoriolisBu(i-1,j) + G%CoriolisBu(i,j-1)) )
      surfFricVel = uStar(i,j)
      surfBuoyFlux = bFlux(i,j)

      ! Compute the OBL thickness
      call cvmix_kpp_compute_OBL_depth( &
           BulkRi(i,j,:),  & ! (in) Bulk Richardson number
           cellHeight,     & ! (in) Height of cells (m)         ???? or interfaces????
           OBLheight,      & ! (out) OBL height (m)
           kOBL,           & ! (out) level of OBL extent
           surfFricVel,    & ! (in) Turbulent friction velocity at surface (m/s)
           surfBuoyFlux,   & ! (in) Buoyancy flux at surface (m2/s3)                  SIGNS???
           Coriolis,       & ! (in) Coriolis parameter (1/s)
           CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters
      OBLdepth(i,j) = -OBLheight ! Change sign for depth/thickness

    enddo ! i
  enddo ! j

  if (CS%id_OBL > 0) call post_data(CS%id_OBL, OBLdepth, CS%diag)
  if (CS%id_BulkRi > 0) call post_data(CS%id_BulkRi, BulkRi, CS%diag)

end subroutine KPP_calculate


subroutine calculateBulkRichardson(CS, G, h, Temp, Salt, u, v, EOS, BulkRi )
! Calculates Bulk richardson number

! Arguments
  type(KPP_CS),                           intent(in)    :: CS     ! Control structure
  type(ocean_grid_type),                  intent(in)    :: G      ! Ocean grid
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h     ! layer/level thicknesses (units of H)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: Temp  ! Pot. temperature (degrees C)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: Salt  ! Salinity (ppt)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)    :: u     ! Velocity components (m/s)
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)    :: v     ! Velocity components (m/s)
  type(EOS_type),                         pointer       :: EOS    ! Equation of state
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(inout) :: BulkRi ! Bulk Richardson number (nondim)

! Local variables
  integer :: i, j, k
  real :: rho1( SZI_(G) ) ! Density of surface properties at depth z
  real :: rho2( SZI_(G) ) ! In-situ density at depth z
  real :: pRef( SZI_(G) ) ! Pressure at top of rho2 cell
  real :: dRef( SZI_(G) ) ! Depth of top of rho2 cell (positve in m)
  real :: dLev( SZI_(G) ) ! Depth of center of rho2 cell (positve in m)
  real :: GoRho, Ut2, Uk, Vk

  GoRho = G%g_Earth / G%Rho0
  Ut2 = 1.e-15 ! Ut2 = Ut**2, Ut is an unresolved vertical shear.  Add to arguments/parameters ???????

  do j = G%jsc, G%jec
    pRef(:) = 0. ! Ignore atmospheric loading in this calculation ????
    dRef(:) = 0.
    do k = 1, G%ke

      ! rho1 is meant to be the average over some [Monin-Obukhov] scale at the surface
      ! In z-mode, this will typically just be the top level. but a proper integral
      ! will be needed for fine vertical resolution or arbitray coordinates.   ???????
      call calculate_density(Temp(:,j,1), Salt(:,j,1), pRef, &
                      rho1, G%isc, G%iec, EOS)
      call calculate_density(Temp(:,j,k), Salt(:,j,k), pRef, &
                      rho2, G%isc, G%iec, EOS)

      dLev(:) = dRef(:) + 0.5 * h(:,j,k) * G%H_to_m ! Depth of center of level k

      do i = G%isc, G%iec
        Uk = 0.5 * ( abs( u(i,j,k) - u(i,j,1) ) + abs( u(i-1,j,k) - u(i-1,j,1) ) ) ! delta_k U
        Vk = 0.5 * ( abs( v(i,j,k) - v(i,j,1) ) + abs( v(i,j-1,k) - v(i,j-1,1) ) ) ! delta_k V
        BulkRi(i,j,k) = ( GoRho * dLev(i) ) * ( rho2(i) - rho1(i) ) / ( ( Uk**2 + Vk**2 ) + Ut2 )
        ! Notes:
        ! o Using dLev includes an extra half layer thickness from surface for all levels
        ! o BulRi(k=1)=0 because rho1=rho2
      enddo ! i

      ! Pressure at bottom of level k will become pressure at top of level on next iteration
      pRef(:) = pRef(:) + G%g_Earth * G%Rho0 * h(:,j,k) ! Boussinesq approximation!!!! ?????
      dRef(:) = dRef(:) + h(:,j,k) * G%H_to_m ! Depth of bottom of level k

    enddo ! k
  enddo ! j

end subroutine calculateBulkRichardson


subroutine KPP_end(CS)
! Clear pointers, dealocate memory
  type(KPP_CS), pointer :: CS ! Control structure

  deallocate(CS)
end subroutine KPP_end

end module MOM_KPP
