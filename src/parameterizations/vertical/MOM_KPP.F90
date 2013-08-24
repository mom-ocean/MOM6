module MOM_KPP

use MOM_coms, only : max_across_PEs
use MOM_checksums, only : hchksum, is_NaN
use MOM_diag_mediator, only : time_type, diag_ctrl, safe_alloc_ptr, post_data
use MOM_diag_mediator, only : query_averaging_enabled, register_diag_field
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_PE
use MOM_EOS, only : EOS_type, calculate_density
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_file_parser, only : openParameterBlock, closeParameterBlock
use MOM_grid, only : ocean_grid_type, isPointInCell

use CVmix_kpp, only : CVmix_init_kpp, CVmix_put_kpp, CVmix_get_kpp_real
use CVmix_kpp, only : CVmix_coeffs_kpp
use CVmix_kpp, only : CVmix_kpp_compute_OBL_depth
use CVmix_kpp, only : CVmix_kpp_compute_turbulent_scales
use CVmix_kpp, only : CVmix_kpp_params_type

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
  real    :: cs         ! Parameter for computing velocity scale function
! real    :: eps        ! small non-negative val (rec 1e-10)
  character(len=10) :: interpType ! Type of iterpolation to use in determining OBL
  logical :: computeEkman ! If True, compute Ekman depth limit
  logical :: computeMoninObukhov ! If True, compute Monin-Obukhov limit
  logical :: passiveMode ! If True, makes KPP passive meaning it does NOT alter the diffusivity
  logical :: debug      ! If True, calculate checksums and write debugging information

  ! CVmix parameters
  type(CVmix_kpp_params_type), pointer :: KPP_params => NULL()

  ! Daignostic handles and pointers
  type(diag_ctrl), pointer :: diag => NULL()
  integer :: id_OBL = -1, id_BulkRi = -1, id_Ws = -1, id_N = -1
  integer :: id_Ut2 = -1, id_BulkUz2 = -1, id_BulkDrho = -1
  integer :: id_uStar = -1, id_buoyFlux = -1
  integer :: id_Kt_KPP = -1, id_Ks_KPP = -1
  integer :: id_NLt_KPP = -1, id_NLs_KPP = -1

end type KPP_CS

! Module data used for debugging only
logical, parameter :: verbose = .True.
#define __DO_SAFETY_CHECKS__

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
  call log_version(paramFile, mod, version, 'This is the MOM wrapper to CVmix:KPP\n' // &
            'See http://code.google.com/p/cvmix/')
  call get_param(paramFile, mod, 'DEBUG', CS%debug, default=.False., do_not_log=.True.)
  call openParameterBlock(paramFile,'KPP')
  call get_param(paramFile, mod, 'PASSIVE', CS%passiveMode,           &
                 'If True, puts KPP into a passive-diagnostic mode.', &
                 default=.False.)
  call get_param(paramFile, mod, 'RI_CRIT', CS%Ri_crit,                       &
                 'Critical Richardson number used to define depth of the\n'// &
                 'Oceab Boundary Layer (OBL).',                               &
                 units='nondim', default=0.3)
  call get_param(paramFile, mod, 'VON_KARMAN', CS%vonKarman, &
                 'von Karman constant.',                     &
                 units='nondim', default=0.40)
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
  call get_param(paramFile, mod, 'CS', CS%cs, &
                 'Parameter for computing velocity scale function.', &
                 units='nondim', default=98.96)
  call closeParameterBlock(paramFile)

  call CVmix_init_kpp( Ri_crit=CS%Ri_crit,                 &
                       vonKarman=CS%vonKarman,             &
                       interp_type=CS%interpType,          &
                       lEkman=CS%computeEkman,             &
                       lMonOb=CS%computeMoninObukhov,      &
                       c_s=CS%cs,                          &
                       CVmix_kpp_params_user=CS%KPP_params )

! Register diagnostics
  CS%diag => diag
  CS%id_OBL = register_diag_field('ocean_model', 'KPP_OBLdepth', diag%axesT1, Time, &
      'Thickness of the surface Ocean Boundary Layer calculated by [CVmix] KPP', 'meter')
  CS%id_BulkRi = register_diag_field('ocean_model', 'KPP_BulkRi', diag%axesTL, Time, &
      'Bulk Richardson number used to find the OBL depth used by [CVmix] KPP', 'nondim')
  CS%id_Ws = register_diag_field('ocean_model', 'KPP_Ws', diag%axesTi, Time, &
      'Turbulent vertical velocity scale for scalars used by [CVmix] KPP', 'm/s')
  CS%id_N = register_diag_field('ocean_model', 'KPP_N', diag%axesTi, Time, &
      'Brunt-Vaisala frequency used by [CVmix] KPP', '1/s')
  CS%id_Ut2 = register_diag_field('ocean_model', 'KPP_Ut2', diag%axesTi, Time, &
      'Unresolved shear turbulence used by [CVmix] KPP', '1/s2')
  CS%id_BulkUz2 = register_diag_field('ocean_model', 'KPP_BulkUz2', diag%axesTL, Time, &
      'Square of bulk difference in resolved velocity used in Bulk Richardson number, as used by [CVmix] KPP', 'm2/s2')
  CS%id_BulkDrho = register_diag_field('ocean_model', 'KPP_BulkDrho', diag%axesTL, Time, &
      'Bulk difference in density used in Bulk Richardson number, as used by [CVmix] KPP', 'kg/m3')
  CS%id_uStar = register_diag_field('ocean_model', 'KPP_uStar', diag%axesT1, Time, &
      'Frictional velocity, u*, as used by [CVmix] KPP', 'm/s')
  CS%id_buoyFlux = register_diag_field('ocean_model', 'KPP_buoyFlux', diag%axesT1, Time, &
      'Buoyancy flux, as used by [CVmix] KPP', 'm2/s3')
  CS%id_Kt_KPP = register_diag_field('ocean_model', 'KPP_Kheat', diag%axesTi, Time, &
      'Heat diffusivity due to KPP, as calculated by [CVmix] KPP', 'm2/s')
  CS%id_Ks_KPP = register_diag_field('ocean_model', 'KPP_Ksalt', diag%axesTi, Time, &
      'Salt diffusivity due to KPP, as calculated by [CVmix] KPP', 'm2/s')
  CS%id_NLt_KPP = register_diag_field('ocean_model', 'KPP_NLtransport_heat', diag%axesTi, Time, &
      'Non-local transport for heat, as calculated by [CVmix] KPP', 'm/s')
  CS%id_NLs_KPP = register_diag_field('ocean_model', 'KPP_NLtransport_salt', diag%axesTi, Time, &
      'Non-local tranpsort for salt, as calculated by [CVmix] KPP', 'm/s')

end subroutine KPP_init


subroutine KPP_calculate(CS, G, h, Temp, Salt, u, v, EOS, uStar, bFlux, Kv)
! Calculates diffusivity and non-local transport for KPP parameterization 

! Arguments
  type(KPP_CS),                           intent(in)    :: CS    ! Control structure
  type(ocean_grid_type),                  intent(in)    :: G     ! Ocean grid
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h     ! Layer/level thicknesses (units of H)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: Temp  ! Pot. temperature (degrees C)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: Salt  ! Salinity (ppt)
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in)    :: u     ! Velocity components (m/s)
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in)    :: v     ! Velocity components (m/s)
  type(EOS_type),                         pointer       :: EOS   ! Equation of state
  real, dimension(NIMEM_,NJMEM_),         intent(in)    :: uStar ! Piston velocity (m/s)
  real, dimension(NIMEM_,NJMEM_),         intent(in)    :: bFlux ! Buoyancy flux (m2/s3)
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(inout) :: Kv ! (in) Vertical diffusivity in interior (m2/s)
                                                                 ! (out) Vertical diffusivity including KPP (m2/s)

! Diagnostics arrays                   should these become allocatables ??????????
  real, dimension( SZI_(G), SZJ_(G), SZK_(G) ) :: BulkRi ! Bulk Richardson number for each layer
  real, dimension( SZI_(G), SZJ_(G) ) :: OBLdepth ! Depth (positive) of OBL (m)
  real, dimension( SZI_(G), SZJ_(G), SZK_(G)+1 ) :: Ws ! Turbulent velocity scale for scalars (m/s)
  real, dimension( SZI_(G), SZJ_(G), SZK_(G)+1 ) :: N ! Brunt-Vaisala frequency (1/s)
  real, dimension( SZI_(G), SZJ_(G), SZK_(G) ) :: dRho ! Bulk difference in density (kg/m3)
  real, dimension( SZI_(G), SZJ_(G), SZK_(G) ) :: Uz2 ! Square of bulk difference in resolved velocity (m2/s2)
  real, dimension( SZI_(G), SZJ_(G), SZK_(G)+1 ) :: Ut2 ! Unresolved shear turbulence (1/s2)
  real, dimension( SZI_(G), SZJ_(G), SZK_(G)+1 ) :: Kt_KPP, Ks_KPP ! Temp/alt diffusivity due to KPP (m2/s)
  real, dimension( SZI_(G), SZJ_(G), SZK_(G)+1 ) :: NLt_KPP, NLs_KPP ! Temp/alt non-local transport (m/s)

! Local variables
  integer :: i, j, k, km1, iteration, largestIterationCount
  real, dimension( G%ke ) :: cellHeight ! Cell center heights referenced to surface (m)
  real, dimension( G%ke+1 ) :: iFaceHeight ! Interface heights referenced to surface (m)
  real, dimension( G%ke+1 ) :: sigmaCoord ! Normalized coordiante, =0 at surface, =1 at z=-OBLd
  real, dimension( G%ke+1 ) :: Ws_1d, Wm_1d ! Profiles of vertical velocity scale for scalars/momentum (m/s)
  real, dimension( G%ke+1 ) :: N_1d ! Brunt-Vaisala frequency, at interfaces (1/s)
  real, dimension( G%ke+1 ) :: Ut2_1d ! Unresolved shear turbulence, at interfaces (1/s2)
  real, dimension( G%ke ) :: BulkRi_1d ! Bulk Richardson number for each layer
  real, dimension( G%ke ) :: deltaRho ! delta Rho as appears in numerator of Bulk Richardson number
  real, dimension( G%ke ) :: deltaU2 ! square of delta U (shear) as appears in denominator of Bulk Richardson number (m2/s2)
  real, dimension( G%ke+1, 2) :: Kdiffusivity ! Vertical diffusivity at interfaces (m2/s)
  real, dimension( G%ke+1 ) :: Kviscosity ! Vertical viscosity at interfaces (m2/s)
  real, dimension( G%ke+1, 2) :: nonLocalTrans ! Non-local transport for heat/salt at interfaces (m/s)
  real :: kOBL, OBLdepth_0d, surfFricVel, surfBuoyFlux, Coriolis, lastOBLdepth, penulOBLdepth
  real :: correction, largestCorrection
  real :: GoRho, pRef, rho1, rhoK, rhoKm1, Uk, Vk, const1, Cv
  real, parameter :: negligibleShear = 1.e-15 ! A small number added to (un)resolved shears to avoid divide by zero
  integer, parameter :: maxIterations = 30 ! Number of iteration on OBL depth to make
  real, parameter :: tolerance = 1.e-4 ! (m) What change in OBL depth is acceptably accurate to stop iterating
  real, parameter :: eps = 0.1 ! Nondimensional extent of surface layer. Used for const1 below.
  real, parameter :: BetaT = -0.2 ! Ratio of entrainment flux to surface buoyancy flux. Used for const1 below.

  BulkRi(:,:,:) = 0.
  OBLdepth(:,:) = 0.
  Ws(:,:,:) = 0.
  N(:,:,:) = 0.
  Ut2(:,:,:) = 0.
  Uz2(:,:,:) = 0.
  Kt_KPP(:,:,:) = 0.
  Ks_KPP(:,:,:) = 0.
  NLt_KPP(:,:,:) = 0.
  NLs_KPP(:,:,:) = 0.

#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
    call hchksum(h, "KPP in: h",G,haloshift=0)
    call hchksum(Temp, "KPP in: T",G,haloshift=0)
    call hchksum(Salt, "KPP in: S",G,haloshift=0)
    call hchksum(u, "KPP in: u",G,haloshift=0)
    call hchksum(v, "KPP in: v",G,haloshift=0)
    call hchksum(uStar, "KPP in: uStar",G,haloshift=0)
    call hchksum(bFlux, "KPP in: bFlux",G,haloshift=0)
    call hchksum(Kv, "KPP in: Kv",G,haloshift=0)
  endif
#endif

  GoRho = G%g_Earth / G%Rho0
  ! const1 is a constant factor in the equation for unresolved shear, Ut (eq. 23 in LMD94)
  const1 = sqrt( abs(BetaT) / (CS%cs * eps) )/( CS%Ri_crit * (CS%vonKarman**2) )

  largestIterationCount = 0
  largestCorrection = 0.
  do j = G%jsc, G%jec
    do i = G%isc, G%iec

      iFaceHeight(1) = 0.
      pRef = 0.
      do k = 1, G%ke
        ! Compute heights, referenced to the surface (z=0)
        cellHeight(k) = iFaceHeight(k) - 0.5 * h(i,j,k) * G%H_to_m ! cell center in metres 
        iFaceHeight(k+1) = iFaceHeight(k) - h(i,j,k) * G%H_to_m ! cell bottom in metres

        ! Compute Bulk Richardson number
        ! rho1 is meant to be the average over some [Monin-Obukhov] scale at the surface
        ! In z-mode, this will typically just be the top level. but a proper integral
        ! will be needed for fine vertical resolution or arbitray coordinates.   ???????
        km1 = max(1, k-1)
        call calculate_density(Temp(i,j,1), Salt(i,j,1), pRef, rho1, EOS)
        call calculate_density(Temp(i,j,k), Salt(i,j,k), pRef, rhoK, EOS)
        call calculate_density(Temp(i,j,km1), Salt(i,j,km1), pRef, rhoKm1, EOS)
        Uk = 0.5 * ( abs( u(i,j,k) - u(i,j,1) ) + abs( u(i-1,j,k) - u(i-1,j,1) ) ) ! delta_k U
        Vk = 0.5 * ( abs( v(i,j,k) - v(i,j,1) ) + abs( v(i,j-1,k) - v(i,j-1,1) ) ) ! delta_k V
        deltaRho(k) = rhoK - rho1
        deltaU2(k) = ( Uk**2 + Vk**2 ) + negligibleShear
        N_1d(k) = sqrt( GoRho * max(rhoK - rhoKm1, 0.) / (0.5*(h(i,j,km1) + h(i,j,k))+G%H_subroundoff) )
        BulkRi_1d(k) = ( GoRho * ( -cellHeight(k) ) ) * deltaRho(k) / deltaU2(k)
        ! Notes:
        ! o Using cellHeight includes an extra half layer thickness from surface for all levels   ????
        ! o BulRi(k=1)=0 because rho1=rhoK

        ! Pressure at bottom of level k will become pressure at top of level on next iteration
        pRef = pRef + G%g_Earth * G%Rho0 * h(i,j,k) ! Boussinesq approximation!!!! ?????
      enddo ! k
      N_1d( G%ke+1 ) = 0.
      Ut2_1d( G%ke+1 ) = 0.

      Coriolis = 0.25*( (G%CoriolisBu(i,j) + G%CoriolisBu(i-1,j-1)) &
                       +(G%CoriolisBu(i-1,j) + G%CoriolisBu(i,j-1)) )
      surfFricVel = uStar(i,j)
      surfBuoyFlux = bFlux(i,j)
! if (isPointInCell(G,i,j,-235.5,13.3)) then
!   print *,'u*,bFlux',surfFricVel,surfBuoyFlux
!   do k=1,G%ke
!     print *,'k,h,z,S,T=',k,h(i,j,k),cellHeight(k),Temp(i,j,k),Salt(i,j,k)
!   enddo
! endif

      OBLdepth_0d = 1.e10 ! Silly initial value
      OBLiterater: do iteration = 0, maxIterations ! Iterate of the estimates of Bulk Ri, Ws and OBL depth

        ! Compute the OBL thickness
        penulOBLdepth = lastOBLdepth ! Store penultimate estimate to catch oscillations in iterator
        lastOBLdepth = OBLdepth_0d ! Record last estimate to measure convergence
        call CVmix_kpp_compute_OBL_depth( &
           BulkRi_1d,              & ! (in) Bulk Richardson number
           iFaceHeight,            & ! (in) Height of interfaces (m)
           OBLdepth_0d,            & ! (out) OBL depth (m)
           kOBL,                   & ! (out) level (+fraction) of OBL extent
           zt_cntr=cellHeight,     & ! (in) Height of cell centers (m)
           surf_fric=surfFricVel,  & ! (in) Turbulent friction velocity at surface (m/s)
           surf_buoy=surfBuoyFlux, & ! (in) Buoyancy flux at surface (m2/s3)
           Coriolis=Coriolis,      & ! (in) Coriolis parameter (1/s)
           CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters
        !if (isPointInCell(G,i,j,-30.5,60.)) print *,'iter,OBL',iteration, OBLdepth_0d
!       OBLdepth_0d = max( OBLdepth_0d, h(i,j,1) ) ! Limit OBL to thicker than top layer ?????

        ! Exit loop if converged, one iteration is guaranteed because initial value
        ! of OBLdepth_0d is silly
        if (abs(OBLdepth_0d - penulOBLdepth) < tolerance .and. & ! Detect oscillatory state
            lastOBLdepth < OBLdepth_0d) then ! Select deeper solution to smooth over problem areas????
          lastOBLdepth = penulOBLdepth ! Force exit through simple criteria
        endif
        correction = abs(OBLdepth_0d - lastOBLdepth)
        if (correction < tolerance) exit OBLiterater ! Simple exit criteria

!       if (iteration > 7 .and. mod(iteration,4) == 0) then ! Slow convergence
!         if ((OBLdepth_0d-lastOBLdepth)*(lastOBLdepth-penulOBLdepth) < 0.) & ! Oscillating
!           OBLdepth_0d = 0.25*( OBLdepth_0d + 2.*lastOBLdepth + penulOBLdepth ) ! Filter guess
!       endif

        ! Now calculate the unresolved turbulence velocity scales
        sigmaCoord(:) = -iFaceHeight/OBLdepth_0d ! =0 at surface, =1 at z=-OBLd
        call CVmix_kpp_compute_turbulent_scales( &
           sigmaCoord,     & ! (in) Normalized boundary layer coordinate (at interfaces)
           OBLdepth_0d,    & ! (in) OBL depth (m)
           surfBuoyFlux,   & ! (in) Buoyancy flux at surface (m2/s3)
           surfFricVel,    & ! (in) Turbulent friction velocity at surface (m/s)
           w_s=Ws_1d,      & ! (out) Turbulent velocity scale profile, at interfaces (m/s)
           CVmix_kpp_params_user=CS%KPP_params ) ! KPP parameters

        ! Re-calculate the Bulk Richardson number adding the turbulent velocity scale
        if ( iteration < maxIterations ) then
          do k = 1, G%ke
            ! Unresolved turbulence shear
            Cv = max( 1.7, 2.1 - 200. * N_1d(k) )
            Ut2_1d(k) = const1 * Cv * (-cellHeight(k)) * N_1d(k) * Ws_1d(k)
            ! Note upward-biased used of Ws since Ws is at interfaces
            BulkRi_1d(k) = ( GoRho * ( -cellHeight(k) ) ) * deltaRho(k) / ( deltaU2(k) + Ut2_1d(k) )
          enddo ! k
        endif

      enddo OBLiterater ! iteration

      if (maxIterations > 0 .and. correction >= tolerance) then
        write(0,*) 'i,j,x,y',i,j,G%geoLonT(i,j),G%geoLatT(i,j)
        write(0,*) 'iters',iteration
        write(0,*) 'penul,last, current OBL',penulOBLdepth,lastOBLdepth,OBLdepth_0d
        write(0,*) 'u*, bFlux',surfFricVel,surfBuoyFlux
        do k = 1, G%ke
          write(0,*) 'k,zw,h,T,S',k,iFaceHeight(k),h(i,j,k),temp(i,j,k),salt(i,j,k)
        enddo
        do k = 1, G%ke
          write(0,*) 'k,h,dRho,dU,Ri_b',k,h(i,j,k),deltaRho(k),sqrt(deltaU2(k)), &
              ( GoRho * ( -cellHeight(k) ) ) * deltaRho(k) / ( deltaU2(k) + Ut2_1d(k) )
        enddo
        call MOM_error(FATAL, 'MOM_KPP, KPP_calculate: '// &
              'The OBL depth iteration failed to converge!!!')
      endif

      if (verbose .and. maxIterations > 0) then
        largestCorrection = max( largestCorrection, correction )
        largestIterationCount = max( largestIterationCount, iteration )
      endif

      ! Now call KPP proper to obtain BL diffusivities, viscosities and non-local transports
      Kdiffusivity(:,1) = Kv(i,j,:) ! Diffusivty for heat ????
      Kdiffusivity(:,2) = Kv(i,j,:) ! Diffusivity for salt ????
      Kviscosity(:) = Kv(i,j,:) ! Viscosity    ???????
      call cvmix_coeffs_kpp(Kdiffusivity, Kviscosity, iFaceHeight, cellHeight, OBLdepth_0d, &
                            kOBL, nonLocalTrans, surfFricVel, surfBuoyFlux,                 &
                            CVmix_kpp_params_user=CS%KPP_params )
#ifdef __DO_SAFETY_CHECKS__
!     if (is_NaN(Kdiffusivity(:,2),skip_mpp=.True.)) then
!       write(0,*) 'i,j=',i,j
!       write(0,*) 'u*,bFlux',surfFricVel, surfBuoyFlux
!       write(0,*) 'OBLd,kOBL',OBLdepth_0d, kOBL
!       do k = 1, G%ke+1
!         write(0,*) 'k,zw,Kin,Kout',k,iFaceHeight(k),Kv(i,j,k),Kdiffusivity(k,2)
!       enddo
!       call MOM_error(FATAL, 'MOM_KPP, KPP_calculate: '// &
!             'NaN detected on return from KPP!!!')
!     endif
#endif

      if (CS%id_OBL > 0) OBLdepth(i,j) = OBLdepth_0d ! Change sign for depth/thickness
      if (CS%id_BulkRi > 0) BulkRi(i,j,:) = BulkRi_1d(:)
      if (CS%id_Ws > 0) Ws(i,j,:) = Ws_1d(:)
      if (CS%id_N > 0) N(i,j,:) = N_1d(:)
      if (CS%id_Ut2 > 0) Ut2(i,j,:) = Ut2_1d(:)
      if (CS%id_BulkUz2 > 0) Uz2(i,j,:) = deltaU2(:)
      if (CS%id_BulkDrho > 0) dRho(i,j,:) = deltaRho(:)
      if (CS%id_Kt_KPP > 0) Kt_KPP(i,j,:) = Kdiffusivity(:,1) - Kv(i,j,:) ! Heat diffusivity due to KPP  (correct index ???)
      if (CS%id_Ks_KPP > 0) Ks_KPP(i,j,:) = Kdiffusivity(:,2) - Kv(i,j,:) ! Salt diffusivity due to KPP  (correct index ???)
      if (CS%id_NLt_KPP > 0) NLt_KPP(i,j,:) = nonLocalTrans(:,1) ! correct index ???
      if (CS%id_NLs_KPP > 0) NLs_KPP(i,j,:) = nonLocalTrans(:,2) ! correct index ???
      if (.not. CS%passiveMode) Kv(i,j,:) = Kdiffusivity(:,2)
    enddo ! i
  enddo ! j

  if (verbose) then
    call max_across_PEs( largestIterationCount )
    call max_across_PEs( largestCorrection )
    if (is_root_PE()) &
      write(*,'("MOM_KPP: max(iter, correction)=",i3,es10.2," m")') largestIterationCount, largestCorrection
  endif

#ifdef __DO_SAFETY_CHECKS__
  if (CS%debug) then
    call hchksum(Kv, "KPP out: Kv",G,haloshift=0)
  endif
#endif

  if (CS%id_OBL > 0) call post_data(CS%id_OBL, OBLdepth, CS%diag)
  if (CS%id_BulkRi > 0) call post_data(CS%id_BulkRi, BulkRi, CS%diag)
  if (CS%id_Ws > 0) call post_data(CS%id_Ws, Ws, CS%diag)
  if (CS%id_N > 0) call post_data(CS%id_N, N, CS%diag)
  if (CS%id_Ut2 > 0) call post_data(CS%id_Ut2, Ut2, CS%diag)
  if (CS%id_BulkUz2 > 0) call post_data(CS%id_BulkUz2, Uz2, CS%diag)
  if (CS%id_BulkDrho > 0) call post_data(CS%id_BulkDrho, dRho, CS%diag)
  if (CS%id_uStar > 0) call post_data(CS%id_uStar, uStar, CS%diag)
  if (CS%id_buoyFlux > 0) call post_data(CS%id_buoyFlux, bFlux, CS%diag)
  if (CS%id_Kt_KPP > 0) call post_data(CS%id_Kt_KPP, Kt_KPP, CS%diag)
  if (CS%id_Ks_KPP > 0) call post_data(CS%id_Ks_KPP, Ks_KPP, CS%diag)
  if (CS%id_NLt_KPP > 0) call post_data(CS%id_NLt_KPP, NLt_KPP, CS%diag)
  if (CS%id_NLs_KPP > 0) call post_data(CS%id_NLs_KPP, NLs_KPP, CS%diag)

end subroutine KPP_calculate


subroutine KPP_end(CS)
! Clear pointers, dealocate memory
  type(KPP_CS), pointer :: CS ! Control structure

  deallocate(CS)
end subroutine KPP_end

end module MOM_KPP
