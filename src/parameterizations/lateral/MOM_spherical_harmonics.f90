module MOM_spherical_harmonics
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
  CLOCK_MODULE
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type

implicit none ; private

public spherical_harmonics_init, spherical_harmonics_end, associatedLegendrePolynomials, SHOrderDegreeToIndex

#include <MOM_memory.h>

type, public :: sht_CS
  integer :: nOrder, lmax
  real, allocatable :: sinLatCell(:,:), cosLatCell(:,:)
  real, allocatable :: complexFactorRe(:,:,:), complexFactorIm(:,:,:), complexExpRe(:,:,:), complexExpIm(:,:,:)
  real, allocatable :: pmnm2(:,:), pmnm1(:,:), pmn(:,:)
  real, allocatable :: aRecurrenceCoeff(:,:), bRecurrenceCoeff(:,:)
  logical :: bfb
end type sht_CS

contains

subroutine spherical_harmonics_init(G, param_file, CS)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(sht_CS), intent(inout) :: CS
  type(param_file_type),   intent(in)    :: param_file !< A structure indicating
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: RADIAN
  integer :: is, ie, js, je
  integer :: i, j
  integer :: m, n
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  PI = 4.0*atan(1.0)
  RADIAN = PI / 180.0

  call get_param(param_file, '', "SHT_ORDER", CS%nOrder, &
                 "Order of spherical harmonics transformation. ", default=1)
  CS%lmax = (CS%nOrder + 1) * (CS%nOrder + 2) / 2
  call get_param(param_file, '', "SHT_BFB", CS%bfb, &
                 "If true, use bfb sum. Default is False.", default=.False.)

  allocate(CS%pmn(is:ie,js:je))  ; CS%pmn(:,:)   = 0.0
  allocate(CS%pmnm1(is:ie,js:je)); CS%pmnm1(:,:) = 0.0
  allocate(CS%pmnm2(is:ie,js:je)); CS%pmnm2(:,:) = 0.0

  ! Compute recurrence relationship coefficients
  allocate(CS%aRecurrenceCoeff(CS%nOrder+1,CS%nOrder+1)); CS%aRecurrenceCoeff(:,:) = 0.0
  allocate(CS%bRecurrenceCoeff(CS%nOrder+1,CS%nOrder+1)); CS%bRecurrenceCoeff(:,:) = 0.0

  do m = 0, CS%nOrder
    do n = m, CS%nOrder
      if (m /= n) then
        CS%aRecurrenceCoeff(n+1,m+1) = sqrt(real((2*n-1)*(2*n+1)) / real((n-m)*(n+m)))
        CS%bRecurrenceCoeff(n+1,m+1) = sqrt(real((2*n+1)*(n+m-1)*(n-m-1)) / real((n-m)*(n+m)*(2*n-3)))
      endif
    enddo
  enddo

  ! Precompute complex exponential factors
  allocate(CS%complexFactorRe(is:ie, js:je, CS%nOrder+1)); CS%complexFactorRe(:,:,:) = 0.0
  allocate(CS%complexFactorIm(is:ie, js:je, CS%nOrder+1)); CS%complexFactorIm(:,:,:) = 0.0
  allocate(CS%complexExpRe(is:ie, js:je, CS%nOrder+1)); CS%complexExpRe(:,:,:) = 0.0
  allocate(CS%complexExpIm(is:ie, js:je, CS%nOrder+1)); CS%complexExpIm(:,:,:) = 0.0

  do m = 0, CS%nOrder
    do j = js,je ; do i = is,ie
      CS%complexExpRe(i, j, m+1)    = cos(real(m) * (G%geolonT(i,j)*RADIAN))
      CS%complexExpIm(i, j, m+1)    = sin(real(m) * (G%geolonT(i,j)*RADIAN))
      CS%complexFactorRe(i, j, m+1) = CS%complexExpRe(i, j, m+1) * G%areaT(i,j) / G%Rad_Earth**2
      CS%complexFactorIm(i, j, m+1) = CS%complexExpIm(i, j, m+1) * G%areaT(i,j) / G%Rad_Earth**2
    enddo ; enddo
  enddo

  ! allocate(Snm_local(2*lmax),Snm(2*lmax))
  ! allocate(SnmRe_local(lmax),SnmRe(lmax))
  ! allocate(SnmIm_local(lmax),SnmIm(lmax))
  ! allocate(SnmIm_local_reproSum(nCellsOwned,lmax))
  ! allocate(SnmRe_local_reproSum(nCellsOwned,lmax))
  ! allocate(Snm_local_reproSum(nCellsOwned,2*lmax))


  ! Pre-compute sin and cos of latCell (co-latitude) values
  allocate(CS%sinLatCell(is:ie,js:je)); CS%sinLatCell(:,:) = 0.0
  allocate(CS%cosLatCell(is:ie,js:je)); CS%cosLatCell(:,:) = 0.0
  do j = js,je ; do i = is,ie
    CS%sinLatCell(i,j) = sin(0.5*PI - G%geolatT(i,j)*RADIAN)
    CS%cosLatCell(i,j) = cos(0.5*PI - G%geolatT(i,j)*RADIAN)
  enddo ; enddo
endsubroutine spherical_harmonics_init

subroutine spherical_harmonics_end(CS)
  type(sht_CS), intent(inout) :: CS

  deallocate(CS%sinLatCell, CS%cosLatCell)
  deallocate(CS%complexFactorRe, CS%complexFactorIm, CS%complexExpRe, CS%complexExpIm)
  deallocate(CS%pmnm2, CS%pmnm1, CS%pmn)
  deallocate(CS%aRecurrenceCoeff, CS%bRecurrenceCoeff)
end subroutine spherical_harmonics_end

subroutine associatedLegendrePolynomials(n, m, l, CS, G)
  type(ocean_grid_type),   intent(in) :: G
  integer, intent(in) :: n
  integer, intent(in) :: m
  integer, intent(out) :: l
  type(sht_CS), intent(inout) :: CS
  ! real, dimension(:,:), intent(inout) :: pmnm2
  ! real, dimension(:,:), intent(inout) :: pmnm1
  ! real, dimension(:,:), intent(inout) :: pmn

  integer :: i, j, k
  real, parameter :: PI = 4.0*atan(1.0)                 ! 3.1415926... calculated as 4*atan(1)
  integer :: is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  l = SHOrderDegreeToIndex(n,m, CS%nOrder)

  if (n == m) then
    do j = js, je ; do i = is, ie
      CS%pmnm2(i,j) = sqrt(1.0/(4.0*PI)) * CS%sinLatCell(i,j)**m
      do k = 1, m
        CS%pmnm2(i,j) = CS%pmnm2(i,j) * sqrt(real(2*k+1)/real(2*k))
      enddo
    enddo ; enddo
  else if (n == m+1) then
    do j = js, je ; do i = is, ie
      CS%pmnm1(i,j) = CS%aRecurrenceCoeff(n+1,m+1) * CS%cosLatCell(i,j) * CS%pmnm2(i,j)
    enddo ; enddo
  else
    do j = js, je ; do i = is, ie
      CS%pmn(i,j) =  CS%aRecurrenceCoeff(n+1,m+1) * CS%cosLatCell(i,j) * CS%pmnm1(i,j) &
                   - CS%bRecurrenceCoeff(n+1,m+1) * CS%pmnm2(i,j)
    enddo ; enddo
  endif
end subroutine associatedLegendrePolynomials

function SHOrderDegreeToIndex(n,m, nOrder) result(l)!{{{

  integer :: l
  integer :: n
  integer :: m
  integer :: nOrder

  l = (nOrder+1)*m - m*(m+1)/2 + n+1

end function SHOrderDegreeToIndex

end module MOM_spherical_harmonics