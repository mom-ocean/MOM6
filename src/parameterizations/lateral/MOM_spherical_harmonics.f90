!> Laplace's spherical harmonics transforms (SHT)
module MOM_spherical_harmonics
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
  CLOCK_MODULE
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_coms_infra, only : sum_across_PEs

implicit none ; private

public spherical_harmonics_init, spherical_harmonics_end, order2index, calc_lmax
public spherical_harmonics_forward, spherical_harmonics_inverse

#include <MOM_memory.h>

type, public :: sht_CS ; private
  logical :: initialized = .False. !< True if this control structure has been initialized.
  integer :: nOrder !< Maximum degree of the spherical harmonics [nodim].
  integer :: lmax !< Number of associated Legendre polynomials of nonnegative m [=(nOrder+1)*(nOrder+2)/2] [nodim].
  real, allocatable :: cosCoLatT(:,:) !< Precomputed cosine of colatitude at the t-cells [nondim].
  real, allocatable :: Pmm(:,:,:) !< Precomputed associated Legendre polynomials of which m=n at the t-cells [nondim].
  real, allocatable :: complexFactorRe(:,:,:), complexFactorIm(:,:,:), & !< Precomputed exponential factors
                       complexExpRe(:,:,:), complexExpIm(:,:,:)          !! at the t-cells [nondim].
  real, allocatable :: aRecurrenceCoeff(:,:), bRecurrenceCoeff(:,:) !< Precomputed recurrennce coefficients [nondim].
  logical :: bfb !< True if use reproducable global sums
end type sht_CS

contains
subroutine spherical_harmonics_forward(G, CS, var, SnmRe, SnmIm, Nd)
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure.
  type(sht_CS), intent(in) :: CS !< Control structure for spherical harmonics trasnforms
  real, intent(in)  :: var(:,:) !< Input 2-D variable
  real, intent(out) :: SnmRe(:), SnmIm(:) !< Real and imaginary SHT coefficients
  integer, intent(in), optional :: Nd !< Maximum degree of the spherical harmonics, overriding nOrder
                                   !! in the control structure.
  ! local variables
  integer :: Nmax ! Local copy of the maximum degree of the spherical harmonics
  integer :: i, j, k
  integer :: is, ie, js, je
  integer :: m, n, l
  real, allocatable :: Snm_local(:), SnmRe_local(:), SnmIm_local(:)
  real :: pmn,   & ! Current associated Legendre polynomials of degree n and order m
          pmnm1, & ! Associated Legendre polynomials of degree n-1 and order m
          pmnm2    ! Associated Legendre polynomials of degree n-2 and order m

  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_spherical_harmonics " // &
    "spherical_harmonics_forward: Module must be initialized before it is used.")

  Nmax = CS%nOrder; if (present(Nd)) Nmax = Nd

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (CS%bfb) then
    ! allocate(Snm_local_reproSum(2*CS%lmax)); Snm_local_reproSum = 0.0
    ! allocate(SnmRe_local_reproSum(CS%lmax)); SnmRe_local_reproSum = 0.0
    ! allocate(SnmIm_local_reproSum(CS%lmax)); SnmIm_local_reproSum = 0.0
  else
    allocate(Snm_local(2*CS%lmax)); Snm_local = 0.0
    allocate(SnmRe_local(CS%lmax)); SnmRe_local = 0.0
    allocate(SnmIm_local(CS%lmax)); SnmIm_local = 0.0
  endif

  do m=0,Nmax
    l = order2index(m, Nmax)

    do j=js,je ; do i=is,ie
      pmn = CS%Pmm(i,j,m+1)
      SnmRe_local(l) = SnmRe_local(l) + var(i,j) * pmn * CS%complexFactorRe(i,j,m+1)
      SnmIm_local(l) = SnmIm_local(l) + var(i,j) * pmn * CS%complexFactorIm(i,j,m+1)

      pmnm2 = 0.0; pmnm1 = pmn
      do n = m+1, Nmax
        pmn =  CS%aRecurrenceCoeff(n+1,m+1) * CS%cosCoLatT(i,j) * pmnm1 - CS%bRecurrenceCoeff(n+1,m+1) * pmnm2
        SnmRe_local(l+n-m) = SnmRe_local(l+n-m) + var(i,j) * pmn * CS%complexFactorRe(i,j,m+1)
        SnmIm_local(l+n-m) = SnmIm_local(l+n-m) + var(i,j) * pmn * CS%complexFactorIm(i,j,m+1)
        pmnm2 = pmnm1; pmnm1 = pmn
      enddo
    enddo ; enddo
  enddo

  ! call mpas_timer_stop('Parallel SAL: Forward Transform')

  ! call mpas_timer_start('Parallel SAL: Communication')
  if (CS%bfb) then
    ! do m = 1,lmax
    !     do iCell = 1,nCellsOwned 
    !       Snm_local_reproSum(iCell,m) = SnmRe_local_reproSum(iCell,m)
    !       Snm_local_reproSum(iCell,lmax+m) = SnmIm_local_reproSum(iCell,m)
    !     enddo
    ! enddo
  else
    do m = 1,CS%lmax
      Snm_local(m) = SnmRe_local(m)
      Snm_local(CS%lmax+m) = SnmIm_local(m)
    enddo
  endif

  ! Compute global integral by summing local contributions
  if (CS%bfb) then
     !threadNum = mpas_threading_get_thread_num()
     !if ( threadNum == 0 ) then
        !  Snm = mpas_global_sum_nfld(Snm_local_reproSum,dminfo%comm)
     !endif
  else
    call sum_across_PEs(Snm_local, 2*CS%lmax)
  endif

  do m=1,CS%lmax
    SnmRe(m) = Snm_local(m)
    SnmIm(m) = Snm_local(CS%lmax+m)
  enddo
end subroutine spherical_harmonics_forward

subroutine spherical_harmonics_inverse(G, CS, SnmRe, SnmIm, var, Nd)
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure.
  type(sht_CS), intent(in) :: CS !< Control structure for spherical harmonics trasnforms
  real, intent(out) :: var(:,:) !< Output 2-D variable
  real, intent(in)  :: SnmRe(:), SnmIm(:) !< Real and imaginary SHT coefficients including
                                          !! any additional scaling factors such as Love numbers
  integer, intent(in), optional :: Nd !< Maximum degree of the spherical harmonics, overriding nOrder
                                   !! in the control structure.
  ! local variables
  integer :: i, j, k
  integer :: is, ie, js, je
  integer :: m, n, l
  integer :: Nmax ! Local copy of the maximum degree of the spherical harmonics
  real :: mFac ! A constant multiplier. mFac = 1 (if m==0) or 2 (if m>0)
  real :: pmn,   & ! Current associated Legendre polynomials of degree n and order m
          pmnm1, & ! Associated Legendre polynomials of degree n-1 and order m
          pmnm2    ! Associated Legendre polynomials of degree n-2 and order m

  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_spherical_harmonics " // &
    "spherical_harmonics_inverse: Module must be initialized before it is used.")
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  Nmax = CS%nOrder; if (present(Nd)) Nmax = Nd

  var = 0.0
  do m=0,Nmax
    mFac = sign(1.0, m-0.5)*0.5 + 1.5
    l = order2index(m, Nmax)

    do j=js,je ; do i=is,ie
      pmn = CS%Pmm(i,j,m+1)
      var(i,j) = var(i,j) &
        + mFac * pmn * (SnmRe(l) * CS%complexExpRe(i,j,m+1) + SnmIm(l) * CS%complexExpIm(i,j,m+1))

      pmnm2 = 0.0; pmnm1 = pmn
      do n=m+1,Nmax
        pmn =  CS%aRecurrenceCoeff(n+1,m+1) * CS%cosCoLatT(i,j) * pmnm1 - CS%bRecurrenceCoeff(n+1,m+1) * pmnm2
        var(i,j) = var(i,j) &
          + mFac * pmn * (SnmRe(l+n-m) * CS%complexExpRe(i,j,m+1) + SnmIm(l+n-m) * CS%complexExpIm(i,j,m+1))
        pmnm2 = pmnm1; pmnm1 = pmn
      enddo
    enddo ; enddo
  enddo
end subroutine spherical_harmonics_inverse

subroutine spherical_harmonics_init(G, param_file, CS)
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure.
  type(param_file_type), intent(in) :: param_file !< A structure indicating
  type(sht_CS), intent(inout) :: CS !< Control structure for spherical harmonics trasnforms

  ! local variables
  real, parameter :: PI = 4.0*atan(1.0) ! 3.1415926... calculated as 4*atan(1) [nodim]
  real, parameter :: RADIAN = PI / 180.0 ! Degree to Radian constant [rad/degree]
  real, dimension(SZI_(G),SZJ_(G)) :: sinCoLatT ! sine of colatitude at the t-cells [nondim].
  real :: Pmm_coef ! = sqrt{ 1.0/(4.0*PI) * prod[(2k+1)/2k)] } [nondim].
  integer :: is, ie, js, je
  integer :: i, j, k
  integer :: m, n
  integer :: Nd_tidal_SAL
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40) :: mdl = "MOM_spherical_harmonics" ! This module's name.

  CS%initialized = .True.

  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "TIDAL_SAL_SHT_DEGREE", Nd_tidal_SAL, &
                 "The maximum degree of the spherical harmonics transformation used for "// &
                 "calculating the self-attraction and loading term for tides.", &
                 default=0, do_not_log=.true.)
  CS%nOrder = Nd_tidal_SAL
  CS%lmax = calc_lmax(CS%nOrder)
  call get_param(param_file, mdl, "SHT_BFB", CS%bfb, &
                 "If true, use bfb sum. Default is False.", default=.False.)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  ! Recurrence relationship coefficients
  !  a has an additional row of zero to accommodate the nonexistent element P(m+2,m+1) when m=n=CS%nOrder
  allocate(CS%aRecurrenceCoeff(CS%nOrder+2,CS%nOrder+1)); CS%aRecurrenceCoeff(:,:) = 0.0
  allocate(CS%bRecurrenceCoeff(CS%nOrder+1,CS%nOrder+1)); CS%bRecurrenceCoeff(:,:) = 0.0
  do m = 0, CS%nOrder
    do n = m, CS%nOrder
      if (m /= n) then
        CS%aRecurrenceCoeff(n+1,m+1) = sqrt(real((2*n-1)*(2*n+1)) / real((n-m)*(n+m)))
        CS%bRecurrenceCoeff(n+1,m+1) = sqrt(real((2*n+1)*(n+m-1)*(n-m-1)) / real((n-m)*(n+m)*(2*n-3)))
      endif
    enddo
  enddo

  ! Complex exponential factors
  allocate(CS%complexFactorRe(is:ie, js:je, CS%nOrder+1)); CS%complexFactorRe(:,:,:) = 0.0
  allocate(CS%complexFactorIm(is:ie, js:je, CS%nOrder+1)); CS%complexFactorIm(:,:,:) = 0.0
  allocate(CS%complexExpRe(is:ie, js:je, CS%nOrder+1)); CS%complexExpRe(:,:,:) = 0.0
  allocate(CS%complexExpIm(is:ie, js:je, CS%nOrder+1)); CS%complexExpIm(:,:,:) = 0.0
  do m=0,CS%nOrder
    do j=js,je ; do i=is,ie
      CS%complexExpRe(i, j, m+1)    = cos(real(m) * (G%geolonT(i,j)*RADIAN))
      CS%complexExpIm(i, j, m+1)    = sin(real(m) * (G%geolonT(i,j)*RADIAN))
      CS%complexFactorRe(i, j, m+1) = CS%complexExpRe(i, j, m+1) * G%areaT(i,j) / G%Rad_Earth**2
      CS%complexFactorIm(i, j, m+1) = CS%complexExpIm(i, j, m+1) * G%areaT(i,j) / G%Rad_Earth**2
    enddo ; enddo
  enddo

  ! sine and cosine of colatitude
  allocate(CS%cosCoLatT(is:ie,js:je)); CS%cosCoLatT(:,:) = 0.0
  do j=js,je ; do i=is,ie
    CS%cosCoLatT(i,j) = cos(0.5*PI - G%geolatT(i,j)*RADIAN)
    sinCoLatT(i,j)    = sin(0.5*PI - G%geolatT(i,j)*RADIAN)
  enddo ; enddo

  ! The diagonal elements of the associated Legendre polynomials (n=m)
  allocate(CS%Pmm(is:ie,js:je,m+1)); CS%Pmm(:,:,:) = 0.0
  do m=0,CS%nOrder
    ! Pmm_coef = 1.0/(4.0*PI)
    ! do k=1,m ; Pmm_coef = Pmm_coef * real(2*k+1)/real(2*k); enddo
    ! Pmm_coef = sqrt(Pmm_coef)
    ! do j=js,je ; do i=is,ie
    !   CS%Pmm(i,j,m+1) = Pmm_coef * sinCoLatT(i,j)**m
    ! enddo ; enddo
    do j=js,je ; do i=is,ie
      CS%Pmm(i,j,m+1) = sqrt(1.0/(4.0*PI)) * sinCoLatT(i,j)**m
      do k = 1, m
        CS%Pmm(i,j,m+1) = CS%Pmm(i,j,m+1) * sqrt(real(2*k+1)/real(2*k))
      enddo
    enddo ; enddo
  enddo
end subroutine spherical_harmonics_init

subroutine spherical_harmonics_end(CS)
  type(sht_CS), intent(inout) :: CS

  deallocate(CS%cosCoLatT)
  deallocate(CS%Pmm)
  deallocate(CS%complexFactorRe, CS%complexFactorIm, CS%complexExpRe, CS%complexExpIm)
  deallocate(CS%aRecurrenceCoeff, CS%bRecurrenceCoeff)
end subroutine spherical_harmonics_end

!> The function calc_lmax returns the number of real elements (cosine) of the spherical harmonics,
!! given the maximum degree,
function calc_lmax(Nd) result(lmax)
  integer :: lmax
  integer, intent(in) :: Nd

  lmax = (Nd+2) * (Nd+1) / 2
end function calc_lmax

!> The function returns the one-dimension index number at (n=0, m=m), given order (m) and maximum degree (Nd)
!! The one-dimensional array is organized following degree being the faster moving dimension.
function order2index(m, Nd) result(l)
  integer :: l
  integer, intent(in) :: m
  integer, intent(in) :: Nd

  l = ((Nd+1) + (Nd+1-(m-1)))*m/2 + 1
end function order2index

end module MOM_spherical_harmonics