!> Laplace's spherical harmonic transforms (SHT)
module MOM_spherical_harmonics
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
                              CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_coms_infra,    only : sum_across_PEs
use MOM_coms,          only : reproducing_sum

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
  logical :: reprod_sum !< True if use reproducable global sums
end type sht_CS

integer :: id_clock_sht=-1 !< CPU clock for SHT [MODULE]
integer :: id_clock_sht_forward=-1 !< CPU clock for forward transforms [ROUTINE]
integer :: id_clock_sht_inverse=-1  !< CPU clock for inverse transforms [ROUTINE]
integer :: id_clock_sht_global_sum=-1  !< CPU clock for global summation in forward transforms [LOOP]

contains
subroutine spherical_harmonics_forward(G, CS, var, SnmRe, SnmIm, Nd)
  type(ocean_grid_type), intent(in)  :: G           !< The ocean's grid structure.
  type(sht_CS),          intent(in)  :: CS          !< Control structure for SHT
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)  :: var(:,:)    !< Input 2-D variable
  real,                  intent(out) :: SnmRe(:), & !< Output real and imaginary SHT coefficients
                                        SnmIm(:)    !! [nondim]
  integer,     optional, intent(in)  :: Nd          !< Maximum degree of the spherical harmonics
                                                    !! overriding nOrder in the CS [nondim]
  ! local variables
  integer :: Nmax ! Local copy of the maximum degree of the spherical harmonics
  integer :: Ltot ! Local copy of the number of spherical harmonics
  real, dimension(SZI_(G),SZJ_(G)) :: &
    pmn,   & ! Current associated Legendre polynomials of degree n and order m
    pmnm1, & ! Associated Legendre polynomials of degree n-1 and order m
    pmnm2    ! Associated Legendre polynomials of degree n-2 and order m
  real, allocatable :: SnmRe_reproSum(:,:,:), SnmIm_reproSum(:,:,:)
  integer :: i, j, k
  integer :: is, ie, js, je
  integer :: m, n, l

  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_spherical_harmonics " // &
    "spherical_harmonics_forward: Module must be initialized before it is used.")

  if (id_clock_sht>0) call cpu_clock_begin(id_clock_sht)
  if (id_clock_sht_forward>0) call cpu_clock_begin(id_clock_sht_forward)

  Nmax = CS%nOrder; if (present(Nd)) Nmax = Nd
  Ltot = calc_lmax(Nmax)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do l=1,Ltot ; SnmRe(l) = 0.0; SnmIm(l) = 0.0 ; enddo

  if (CS%reprod_sum) then
    allocate(SnmRe_reproSum(is:ie, js:je, CS%lmax)); SnmRe_reproSum = 0.0
    allocate(SnmIm_reproSum(is:ie, js:je, CS%lmax)); SnmIm_reproSum = 0.0

    do m=0,Nmax
      l = order2index(m, Nmax)

      do j=js,je ; do i=is,ie
        SnmRe_reproSum(i,j,l) = var(i,j) * CS%Pmm(i,j,m+1) * CS%complexFactorRe(i,j,m+1)
        SnmIm_reproSum(i,j,l) = var(i,j) * CS%Pmm(i,j,m+1) * CS%complexFactorIm(i,j,m+1)
        pmnm2(i,j) = 0.0
        pmnm1(i,j) = CS%Pmm(i,j,m+1)
      enddo ; enddo

      do n = m+1, Nmax ; do j=js,je ; do i=is,ie
        pmn(i,j) = CS%aRecurrenceCoeff(n+1,m+1) * CS%cosCoLatT(i,j) * pmnm1(i,j) - CS%bRecurrenceCoeff(n+1,m+1) * pmnm2(i,j)
        SnmRe_reproSum(i,j,l+n-m) = var(i,j) * pmn(i,j) * CS%complexFactorRe(i,j,m+1)
        SnmIm_reproSum(i,j,l+n-m) = var(i,j) * pmn(i,j) * CS%complexFactorIm(i,j,m+1)
        pmnm2(i,j) = pmnm1(i,j)
        pmnm1(i,j) = pmn(i,j)
      enddo ; enddo ; enddo
    enddo
  else
    do m=0,Nmax
      l = order2index(m, Nmax)

      do j=js,je ; do i=is,ie
        SnmRe(l) = SnmRe(l) + var(i,j) * CS%Pmm(i,j,m+1) * CS%complexFactorRe(i,j,m+1)
        SnmIm(l) = SnmIm(l) + var(i,j) * CS%Pmm(i,j,m+1) * CS%complexFactorIm(i,j,m+1)
        pmnm2(i,j) = 0.0
        pmnm1(i,j) = CS%Pmm(i,j,m+1)
      enddo ; enddo

      do n=m+1, Nmax ; do j=js,je ; do i=is,ie
        pmn(i,j) = CS%aRecurrenceCoeff(n+1,m+1) * CS%cosCoLatT(i,j) * pmnm1(i,j) - CS%bRecurrenceCoeff(n+1,m+1) * pmnm2(i,j)
        SnmRe(l+n-m) = SnmRe(l+n-m) + var(i,j) * pmn(i,j) * CS%complexFactorRe(i,j,m+1)
        SnmIm(l+n-m) = SnmIm(l+n-m) + var(i,j) * pmn(i,j) * CS%complexFactorIm(i,j,m+1)
        pmnm2(i,j) = pmnm1(i,j)
        pmnm1(i,j) = pmn(i,j)
      enddo ; enddo ; enddo
    enddo
  endif

  if (id_clock_sht_global_sum>0) call cpu_clock_begin(id_clock_sht_global_sum)

  if (CS%reprod_sum) then
    do l=1,Ltot
      SnmRe(l) = reproducing_sum(SnmRe_reproSum(:,:,l))
      SnmIm(l) = reproducing_sum(SnmIm_reproSum(:,:,l))
    enddo
  else
    call sum_across_PEs(SnmRe, Ltot)
    call sum_across_PEs(SnmIm, Ltot)
  endif

  if (id_clock_sht_global_sum>0) call cpu_clock_end(id_clock_sht_global_sum)
  if (id_clock_sht_forward>0) call cpu_clock_end(id_clock_sht_forward)
  if (id_clock_sht>0) call cpu_clock_end(id_clock_sht)
end subroutine spherical_harmonics_forward

subroutine spherical_harmonics_inverse(G, CS, SnmRe, SnmIm, var, Nd)
  type(ocean_grid_type), intent(in)  :: G           !< The ocean's grid structure.
  type(sht_CS),          intent(in)  :: CS          !< Control structure for SHT
  real,                  intent(in)  :: SnmRe(:), & !< Real and imaginary SHT coefficients with
                                        SnmIm(:)    !! any scaling factors such as Love numbers [nondim]
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(out) :: var(:,:)    !< Output 2-D variable
  integer,     optional, intent(in)  :: Nd          !< Maximum degree of the spherical harmonics
                                                    !! overriding nOrder in the CS [nondim]
  ! local variables
  integer :: Nmax ! Local copy of the maximum degree of the spherical harmonics
  real    :: mFac ! A constant multiplier. mFac = 1 (if m==0) or 2 (if m>0)
  real, dimension(SZI_(G),SZJ_(G)) :: &
    pmn,   & ! Current associated Legendre polynomials of degree n and order m
    pmnm1, & ! Associated Legendre polynomials of degree n-1 and order m
    pmnm2    ! Associated Legendre polynomials of degree n-2 and order m
  integer :: i, j, k
  integer :: is, ie, js, je
  integer :: m, n, l

  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_spherical_harmonics " // &
    "spherical_harmonics_inverse: Module must be initialized before it is used.")

  if (id_clock_sht>0) call cpu_clock_begin(id_clock_sht)
  if (id_clock_sht_inverse>0) call cpu_clock_begin(id_clock_sht_inverse)

  Nmax = CS%nOrder; if (present(Nd)) Nmax = Nd

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  var = 0.0
  do m=0,Nmax
    mFac = sign(1.0, m-0.5)*0.5 + 1.5
    l = order2index(m, Nmax)

    do j=js,je ; do i=is,ie
      var(i,j) = var(i,j) &
        + mFac * CS%Pmm(i,j,m+1) * (SnmRe(l) * CS%complexExpRe(i,j,m+1) + SnmIm(l) * CS%complexExpIm(i,j,m+1))
      pmnm2(i,j) = 0.0
      pmnm1(i,j) = CS%Pmm(i,j,m+1)
    enddo ; enddo

    do n=m+1,Nmax ; do j=js,je ; do i=is,ie
      pmn(i,j) = CS%aRecurrenceCoeff(n+1,m+1) * CS%cosCoLatT(i,j) * pmnm1(i,j) - CS%bRecurrenceCoeff(n+1,m+1) * pmnm2(i,j)
      var(i,j) = var(i,j) &
        + mFac * pmn(i,j) * (SnmRe(l+n-m) * CS%complexExpRe(i,j,m+1) + SnmIm(l+n-m) * CS%complexExpIm(i,j,m+1))
      pmnm2(i,j) = pmnm1(i,j)
      pmnm1(i,j) = pmn(i,j)
    enddo ; enddo ; enddo
  enddo

  if (id_clock_sht_inverse>0) call cpu_clock_end(id_clock_sht_inverse)
  if (id_clock_sht>0) call cpu_clock_end(id_clock_sht)
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

  if (CS%initialized) return
  CS%initialized = .True.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "TIDAL_SAL_SHT_DEGREE", Nd_tidal_SAL, &
                 "The maximum degree of the spherical harmonics transformation used for "// &
                 "calculating the self-attraction and loading term for tides.", &
                 default=0, do_not_log=.true.)
  CS%nOrder = Nd_tidal_SAL
  CS%lmax = calc_lmax(CS%nOrder)
  call get_param(param_file, mdl, "SHT_REPRODUCING_SUM", CS%reprod_sum, &
                 "If true, use reproducing sums (invariant to PE layout) in inverse transform "// &
                 "of spherical harmonics. Otherwise use a simple sum of floationg point numbers. ", default=.False.)

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

  id_clock_sht = cpu_clock_id('(Ocean spherical harmonics)', grain=CLOCK_MODULE)
  id_clock_sht_forward = cpu_clock_id('(Ocean SHT forward)', grain=CLOCK_ROUTINE)
  id_clock_sht_inverse = cpu_clock_id('(Ocean SHT inverse)', grain=CLOCK_ROUTINE)
  id_clock_sht_global_sum = cpu_clock_id('(Ocean SHT global sum)', grain=CLOCK_LOOP)

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