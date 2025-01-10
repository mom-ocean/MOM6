!> Laplace's spherical harmonic transforms (SHT)
module MOM_spherical_harmonics
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
                              CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_coms_infra,    only : sum_across_PEs
use MOM_coms,          only : reproducing_sum

implicit none ; private

public spherical_harmonics_init, spherical_harmonics_end, order2index, calc_lmax
public spherical_harmonics_forward, spherical_harmonics_inverse

#include <MOM_memory.h>

!> Control structure for spherical harmonic transforms
type, public :: sht_CS ; private
  logical :: initialized = .False. !< True if this control structure has been initialized.
  integer :: ndegree !< Maximum degree of the spherical harmonics [nondim].
  integer :: lmax !< Number of associated Legendre polynomials of nonnegative m
                  !! [lmax=(ndegree+1)*(ndegree+2)/2] [nondim].
  real, allocatable :: cos_clatT(:,:) !< Precomputed cosine of colatitude at the t-cells [nondim].
  real, allocatable :: Pmm(:,:,:) !< Precomputed associated Legendre polynomials (m=n) at the t-cells [nondim].
  real, allocatable :: cos_lonT(:,:,:), & !< Precomputed cosine factors at the t-cells [nondim].
                       sin_lonT(:,:,:)    !< Precomputed sine factors at the t-cells [nondim].
  real, allocatable :: cos_lonT_wtd(:,:,:), & !< Precomputed area-weighted cosine factors at the t-cells [nondim]
                       sin_lonT_wtd(:,:,:)    !< Precomputed area-weighted sine factors at the t-cells [nondim]
  real, allocatable :: a_recur(:,:), & !< Precomputed recurrence coefficients a [nondim].
                       b_recur(:,:)    !< Precomputed recurrence coefficients b [nondim].
  real, allocatable :: Snm_Re_raw(:,:,:), & !< Array to store un-summed SHT coefficients
                       Snm_Im_raw(:,:,:)    !< at the t-cells for reproducing sums [same as input variable]
  logical :: reprod_sum !< True if use reproducible global sums
end type sht_CS

integer :: id_clock_sht=-1 !< CPU clock for SHT [MODULE]
integer :: id_clock_sht_forward=-1 !< CPU clock for forward transforms [ROUTINE]
integer :: id_clock_sht_inverse=-1  !< CPU clock for inverse transforms [ROUTINE]
integer :: id_clock_sht_global_sum=-1  !< CPU clock for global summation in forward transforms [LOOP]

contains

!> Calculates forward spherical harmonics transforms
subroutine spherical_harmonics_forward(G, CS, var, Snm_Re, Snm_Im, Nd, tmp_scale)
  type(ocean_grid_type), intent(in)    :: G            !< The ocean's grid structure.
  type(sht_CS),          intent(inout) :: CS           !< Control structure for SHT
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(in)    :: var          !< Input 2-D variable [A]
  real,                  intent(out)   :: Snm_Re(:)    !< SHT coefficients for the real modes (cosine) [A]
  real,                  intent(out)   :: Snm_Im(:)    !< SHT coefficients for the imaginary modes (sine) [A]
  integer,     optional, intent(in)    :: Nd           !< Maximum degree of the spherical harmonics
                                                       !! overriding ndegree in the CS [nondim]
  real,        optional, intent(in)    :: tmp_scale    !< A temporary rescaling factor to convert
                                                       !! var to MKS units during the reproducing
                                                       !! sums [a A-1 ~> 1]
  ! local variables
  integer :: Nmax ! Local copy of the maximum degree of the spherical harmonics
  integer :: Ltot ! Local copy of the number of spherical harmonics
  real, dimension(SZI_(G),SZJ_(G)) :: &
    pmn,   & ! Current associated Legendre polynomials of degree n and order m [nondim]
    pmnm1, & ! Associated Legendre polynomials of degree n-1 and order m [nondim]
    pmnm2    ! Associated Legendre polynomials of degree n-2 and order m [nondim]
  real :: scale ! A rescaling factor to temporarily convert var to MKS units during the
                ! reproducing sums [a A-1 ~> 1]
  real :: I_scale ! The inverse of scale [A a-1 ~> 1]
  real :: sum_tot ! The total of all components output by the reproducing sum in arbitrary units [a]
  integer :: i, j, k
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: m, n, l

  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_spherical_harmonics " // &
    "spherical_harmonics_forward: Module must be initialized before it is used.")

  if (id_clock_sht>0) call cpu_clock_begin(id_clock_sht)
  if (id_clock_sht_forward>0) call cpu_clock_begin(id_clock_sht_forward)

  Nmax = CS%ndegree; if (present(Nd)) Nmax = Nd
  Ltot = calc_lmax(Nmax)

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed ; do i=isd,ied
    pmn(i,j) = 0.0; pmnm1(i,j) = 0.0; pmnm2(i,j) = 0.0
  enddo ; enddo

  do l=1,Ltot ; Snm_Re(l) = 0.0; Snm_Im(l) = 0.0 ; enddo

  if (CS%reprod_sum) then
    scale = 1.0 ; if (present(tmp_scale)) scale = tmp_scale
    do m=0,Nmax
      l = order2index(m, Nmax)

      do j=js,je ; do i=is,ie
        CS%Snm_Re_raw(i,j,l) = (scale*var(i,j)) * CS%Pmm(i,j,m+1) * CS%cos_lonT_wtd(i,j,m+1)
        CS%Snm_Im_raw(i,j,l) = (scale*var(i,j)) * CS%Pmm(i,j,m+1) * CS%sin_lonT_wtd(i,j,m+1)
        pmnm2(i,j) = 0.0
        pmnm1(i,j) = CS%Pmm(i,j,m+1)
      enddo ; enddo

      do n = m+1, Nmax ; do j=js,je ; do i=is,ie
        pmn(i,j) = &
          CS%a_recur(n+1,m+1) * CS%cos_clatT(i,j) * pmnm1(i,j) - CS%b_recur(n+1,m+1) * pmnm2(i,j)
        CS%Snm_Re_raw(i,j,l+n-m) = (scale*var(i,j)) * pmn(i,j) * CS%cos_lonT_wtd(i,j,m+1)
        CS%Snm_Im_raw(i,j,l+n-m) = (scale*var(i,j)) * pmn(i,j) * CS%sin_lonT_wtd(i,j,m+1)
        pmnm2(i,j) = pmnm1(i,j)
        pmnm1(i,j) = pmn(i,j)
      enddo ; enddo ; enddo
    enddo
  else
    do m=0,Nmax
      l = order2index(m, Nmax)

      do j=js,je ; do i=is,ie
        Snm_Re(l) = Snm_Re(l) + var(i,j) * CS%Pmm(i,j,m+1) * CS%cos_lonT_wtd(i,j,m+1)
        Snm_Im(l) = Snm_Im(l) + var(i,j) * CS%Pmm(i,j,m+1) * CS%sin_lonT_wtd(i,j,m+1)
        pmnm2(i,j) = 0.0
        pmnm1(i,j) = CS%Pmm(i,j,m+1)
      enddo ; enddo

      do n=m+1, Nmax ; do j=js,je ; do i=is,ie
        pmn(i,j) = &
          CS%a_recur(n+1,m+1) * CS%cos_clatT(i,j) * pmnm1(i,j) - CS%b_recur(n+1,m+1) * pmnm2(i,j)
        Snm_Re(l+n-m) = Snm_Re(l+n-m) + var(i,j) * pmn(i,j) * CS%cos_lonT_wtd(i,j,m+1)
        Snm_Im(l+n-m) = Snm_Im(l+n-m) + var(i,j) * pmn(i,j) * CS%sin_lonT_wtd(i,j,m+1)
        pmnm2(i,j) = pmnm1(i,j)
        pmnm1(i,j) = pmn(i,j)
      enddo ; enddo ; enddo
    enddo
  endif

  if (id_clock_sht_global_sum>0) call cpu_clock_begin(id_clock_sht_global_sum)

  if (CS%reprod_sum) then
    sum_tot = reproducing_sum(CS%Snm_Re_raw(:,:,1:Ltot), sums=Snm_Re(1:Ltot))
    sum_tot = reproducing_sum(CS%Snm_Im_raw(:,:,1:Ltot), sums=Snm_Im(1:Ltot))
    if (scale /= 1.0) then
      I_scale = 1.0 / scale
      do l=1,Ltot
        Snm_Re(l) = I_scale * Snm_Re(l)
        Snm_Im(l) = I_scale * Snm_Im(l)
      enddo
    endif
  else
    call sum_across_PEs(Snm_Re, Ltot)
    call sum_across_PEs(Snm_Im, Ltot)
  endif

  if (id_clock_sht_global_sum>0) call cpu_clock_end(id_clock_sht_global_sum)
  if (id_clock_sht_forward>0) call cpu_clock_end(id_clock_sht_forward)
  if (id_clock_sht>0) call cpu_clock_end(id_clock_sht)
end subroutine spherical_harmonics_forward

!> Calculates inverse spherical harmonics transforms
subroutine spherical_harmonics_inverse(G, CS, Snm_Re, Snm_Im, var, Nd)
  type(ocean_grid_type), intent(in)  :: G            !< The ocean's grid structure.
  type(sht_CS),          intent(in)  :: CS           !< Control structure for SHT
  real,                  intent(in)  :: Snm_Re(:)    !< SHT coefficients for the real modes (cosine) [A]
  real,                  intent(in)  :: Snm_Im(:)    !< SHT coefficients for the imaginary modes (sine) [A]
  real, dimension(SZI_(G),SZJ_(G)), &
                         intent(out) :: var          !< Output 2-D variable [A]
  integer,     optional, intent(in)  :: Nd           !< Maximum degree of the spherical harmonics
                                                     !! overriding ndegree in the CS [nondim]
  ! local variables
  integer :: Nmax ! Local copy of the maximum degree of the spherical harmonics [nondim]
  real    :: mFac ! A constant multiplier. mFac = 1 (if m==0) or 2 (if m>0) [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    pmn,   & ! Current associated Legendre polynomials of degree n and order m [nondim]
    pmnm1, & ! Associated Legendre polynomials of degree n-1 and order m [nondim]
    pmnm2    ! Associated Legendre polynomials of degree n-2 and order m [nondim]
  integer :: i, j, k
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: m, n, l

  if (.not.CS%initialized) call MOM_error(FATAL, "MOM_spherical_harmonics " // &
    "spherical_harmonics_inverse: Module must be initialized before it is used.")

  if (id_clock_sht>0) call cpu_clock_begin(id_clock_sht)
  if (id_clock_sht_inverse>0) call cpu_clock_begin(id_clock_sht_inverse)

  Nmax = CS%ndegree; if (present(Nd)) Nmax = Nd

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed ; do i=isd,ied
    pmn(i,j) = 0.0; pmnm1(i,j) = 0.0; pmnm2(i,j) = 0.0
    var(i,j) = 0.0
  enddo ; enddo

  do m=0,Nmax
    mFac = sign(1.0, m-0.5)*0.5 + 1.5
    l = order2index(m, Nmax)

    do j=js,je ; do i=is,ie
      var(i,j) = var(i,j) &
        + mFac * CS%Pmm(i,j,m+1) * (  Snm_Re(l) * CS%cos_lonT(i,j,m+1) &
                                    + Snm_Im(l) * CS%sin_lonT(i,j,m+1))
      pmnm2(i,j) = 0.0
      pmnm1(i,j) = CS%Pmm(i,j,m+1)
    enddo ; enddo

    do n=m+1,Nmax ; do j=js,je ; do i=is,ie
      pmn(i,j) = &
        CS%a_recur(n+1,m+1) * CS%cos_clatT(i,j) * pmnm1(i,j) - CS%b_recur(n+1,m+1) * pmnm2(i,j)
      var(i,j) = var(i,j) &
        + mFac * pmn(i,j) * (  Snm_Re(l+n-m) * CS%cos_lonT(i,j,m+1) &
                             + Snm_Im(l+n-m) * CS%sin_lonT(i,j,m+1))
      pmnm2(i,j) = pmnm1(i,j)
      pmnm1(i,j) = pmn(i,j)
    enddo ; enddo ; enddo
  enddo

  if (id_clock_sht_inverse>0) call cpu_clock_end(id_clock_sht_inverse)
  if (id_clock_sht>0) call cpu_clock_end(id_clock_sht)
end subroutine spherical_harmonics_inverse

!> Calculate precomputed coefficients
subroutine spherical_harmonics_init(G, param_file, CS)
  type(ocean_grid_type), intent(in) :: G !< The ocean's grid structure.
  type(param_file_type), intent(in) :: param_file !< A structure indicating
  type(sht_CS), intent(inout)       :: CS !< Control structure for spherical harmonic transforms

  ! local variables
  real, parameter :: PI = 4.0*atan(1.0) ! 3.1415926... calculated as 4*atan(1) [nondim]
  real, parameter :: RADIAN = PI / 180.0 ! Degree to Radian constant [rad/degree]
  real, dimension(SZI_(G),SZJ_(G)) :: sin_clatT ! sine of colatitude at the t-cells [nondim].
  real :: Pmm_coef ! = sqrt{ 1.0/(4.0*PI) * prod[(2k+1)/2k)] } [nondim].
  integer :: is, ie, js, je
  integer :: i, j, k
  integer :: m, n
  integer :: Nd_SAL ! Maximum degree for SAL
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40) :: mdl = "MOM_spherical_harmonics" ! This module's name.

  if (CS%initialized) return
  CS%initialized = .True.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SAL_HARMONICS_DEGREE", Nd_SAL, "", default=0, do_not_log=.true.)
  CS%ndegree = Nd_SAL
  CS%lmax = calc_lmax(CS%ndegree)
  call get_param(param_file, mdl, "SHT_REPRODUCING_SUM", CS%reprod_sum, &
                 "If true, use reproducing sums (invariant to PE layout) in inverse transform "// &
                 "of spherical harmonics. Otherwise use a simple sum of floating point numbers. ", &
                 default=.False.)

  ! Calculate recurrence relationship coefficients
  allocate(CS%a_recur(CS%ndegree+1, CS%ndegree+1)); CS%a_recur(:,:) = 0.0
  allocate(CS%b_recur(CS%ndegree+1, CS%ndegree+1)); CS%b_recur(:,:) = 0.0
  do m=0,CS%ndegree ; do n=m+1,CS%ndegree
    ! These expressione will give NaNs with 32-bit integers for n > 23170, but this is trapped elsewhere.
    CS%a_recur(n+1,m+1) = sqrt(real((2*n-1) * (2*n+1)) / real((n-m) * (n+m)))
    CS%b_recur(n+1,m+1) = sqrt((real(2*n+1) * real((n+m-1) * (n-m-1))) / (real((n-m) * (n+m)) * real(2*n-3)))
  enddo ; enddo

  ! Calculate complex exponential factors
  allocate(CS%cos_lonT_wtd(is:ie, js:je, CS%ndegree+1)); CS%cos_lonT_wtd(:,:,:) = 0.0
  allocate(CS%sin_lonT_wtd(is:ie, js:je, CS%ndegree+1)); CS%sin_lonT_wtd(:,:,:) = 0.0
  allocate(CS%cos_lonT(is:ie, js:je, CS%ndegree+1)); CS%cos_lonT(:,:,:) = 0.0
  allocate(CS%sin_lonT(is:ie, js:je, CS%ndegree+1)); CS%sin_lonT(:,:,:) = 0.0
  do m=0,CS%ndegree
    do j=js,je ; do i=is,ie
      CS%cos_lonT(i,j,m+1)     = cos(real(m) * (G%geolonT(i,j)*RADIAN))
      CS%sin_lonT(i,j,m+1)     = sin(real(m) * (G%geolonT(i,j)*RADIAN))
      CS%cos_lonT_wtd(i,j,m+1) = CS%cos_lonT(i,j,m+1) * G%areaT(i,j) / G%Rad_Earth_L**2
      CS%sin_lonT_wtd(i,j,m+1) = CS%sin_lonT(i,j,m+1) * G%areaT(i,j) / G%Rad_Earth_L**2
    enddo ; enddo
  enddo

  ! Calculate sine and cosine of colatitude
  allocate(CS%cos_clatT(is:ie, js:je)); CS%cos_clatT(:,:) = 0.0
  do j=js,je ; do i=is,ie
    CS%cos_clatT(i,j) = cos(0.5*PI - G%geolatT(i,j)*RADIAN)
    sin_clatT(i,j)    = sin(0.5*PI - G%geolatT(i,j)*RADIAN)
  enddo ; enddo

  ! Calculate the diagonal elements of the associated Legendre polynomials (n=m)
  allocate(CS%Pmm(is:ie,js:je,m+1)); CS%Pmm(:,:,:) = 0.0
  do m=0,CS%ndegree
    Pmm_coef = 1.0/(4.0*PI)
    do k=1,m ; Pmm_coef = Pmm_coef * (real(2*k+1) / real(2*k)); enddo
    Pmm_coef = sqrt(Pmm_coef)
    do j=js,je ; do i=is,ie
      CS%Pmm(i,j,m+1) = Pmm_coef * (sin_clatT(i,j)**m)
    enddo ; enddo
  enddo

  if (CS%reprod_sum) then
    allocate(CS%Snm_Re_raw(is:ie, js:je, CS%lmax)); CS%Snm_Re_raw = 0.0
    allocate(CS%Snm_Im_raw(is:ie, js:je, CS%lmax)); CS%Snm_Im_raw = 0.0
  endif

  id_clock_sht = cpu_clock_id('(Ocean spherical harmonics)', grain=CLOCK_MODULE)
  id_clock_sht_forward = cpu_clock_id('(Ocean SHT forward)', grain=CLOCK_ROUTINE)
  id_clock_sht_inverse = cpu_clock_id('(Ocean SHT inverse)', grain=CLOCK_ROUTINE)
  id_clock_sht_global_sum = cpu_clock_id('(Ocean SHT global sum)', grain=CLOCK_LOOP)

end subroutine spherical_harmonics_init

!> Deallocate any variables allocated in spherical_harmonics_init
subroutine spherical_harmonics_end(CS)
  type(sht_CS), intent(inout) :: CS !< Control structure for spherical harmonic transforms

  deallocate(CS%cos_clatT)
  deallocate(CS%Pmm)
  deallocate(CS%cos_lonT_wtd, CS%sin_lonT_wtd, CS%cos_lonT, CS%sin_lonT)
  deallocate(CS%a_recur, CS%b_recur)
  if (CS%reprod_sum) &
    deallocate(CS%Snm_Re_raw, CS%Snm_Im_raw)
end subroutine spherical_harmonics_end

!> Calculates the number of real elements (cosine) of spherical harmonics given maximum degree Nd.
function calc_lmax(Nd) result(lmax)
  integer :: lmax           !< Number of real spherical harmonic modes [nondim]
  integer, intent(in) :: Nd !< Maximum degree [nondim]

  lmax = (Nd+2) * (Nd+1) / 2
end function calc_lmax

!> Calculates the one-dimensional index number at (n=0, m=m), given order m and maximum degree Nd.
!! It is sequenced with degree (n) changing first and order (m) changing second.
function order2index(m, Nd) result(l)
  integer :: l              !< One-dimensional index number [nondim]
  integer, intent(in) :: m  !< Current order number [nondim]
  integer, intent(in) :: Nd !< Maximum degree [nondim]

  l = ((Nd+1) + (Nd+1-(m-1)))*m/2 + 1
end function order2index

!> \namespace mom_spherical_harmonics
!!
!! \section section_spherical_harmonics Spherical harmonics
!!
!! This module contains the subroutines to calculate spherical harmonic transforms (SHT), namely, forward transform
!! of a two-dimensional field into a given number of spherical harmonic modes and its inverse transform.  This module
!! is primarily used to but not limited to calculate self-attraction and loading (SAL) term, which is mostly relevant to
!! high frequency motions such as tides.  Should other needs arise in the future, this API can be easily modified.
!! Currently, the transforms are for t-cell fields only.
!!
!! This module is stemmed from SAL calculation in Model for Prediction Across Scales (MPAS)-Ocean developed by Los
!! Alamos National Laboratory and University of Michigan [\cite Barton2022 and \cite Brus2023]. The algorithm
!! for forward and inverse transforms loosely follows \cite Schaeffer2013.
!!
!! In forward transform, a two-dimensional physical field can be projected into a series of spherical harmonics. The
!! spherical harmonic coefficient of degree n and order m for a field \f$f(\theta, \phi)\f$ is calculated as follows:
!! \f[
!!   f^m_n = \int^{2\pi}_{0}\int^{\pi}_{0}f(\theta,\phi)Y^m_n(\theta,\phi)\sin\theta d\theta d\phi
!! \f]
!! and
!! \f[
!!  Y^m_n(\theta,\phi) = P^m_n(\cos\theta)\exp(im\phi)
!! \f]
!! where \f$P^m_n(\cos \theta)\f$ is the normalized associated Legendre polynomial of degree n and order m. \f$\phi\f$
!! is the longitude and \f$\theta\f$ is the colatitude.
!! Or, written in the discretized form:
!! \f[
!!  f^m_n = \sum^{Nj}_{0}\sum^{Ni}_{0}f(i,j)Y^m_n(i,j)A(i,j)/r_e^2
!! \f]
!! where \f$A\f$ is the area of the cell and \f$r_e\f$ is the radius of the Earth.
!!
!! In inverse transform, the first N degree spherical harmonic coefficients are used to reconstruct a two-dimensional
!! physical field:
!! \f[
!!   f(\theta,\phi) = \sum^N_{n=0}\sum^{n}_{m=-n}f^m_nY^m_n(\theta,\phi)
!! \f]
!!
!! The exponential coefficients are pre-computed and stored in the memory. The associated Legendre polynomials are
!! computed "on-the-fly", using the recurrence relationships to avoid large memory usage and take the advantage of
!! array vectorization.
!!
!! The maximum degree of the spherical harmonics is a runtime parameter and the maximum used by all SHT applications.
!! At the moment, it is only decided by <code>SAL_HARMONICS_DEGREE</code>.
!!
!! The forward transforms involve a global summation. Runtime flag <code>SHT_REPRODUCING_SUM</code> controls
!! whether this is done in a bit-wise reproducing way or not.
!!
!! References:
!!
!! Barton, K.N., Pal, N., Brus, S.R., Petersen, M.R., Arbic, B.K., Engwirda, D., Roberts, A.F., Westerink, J.J.,
!! Wirasaet, D. and Schindelegger, M., 2022. Global Barotropic Tide Modeling Using Inline Self‐Attraction and Loading in
!! MPAS‐Ocean. Journal of Advances in Modeling Earth Systems, 14(11), p.e2022MS003207.
!! https://doi.org/10.1029/2022MS003207
!!
!! Brus, S.R., Barton, K.N., Pal, N., Roberts, A.F., Engwirda, D., Petersen, M.R., Arbic, B.K., Wirasaet, D.,
!! Westerink, J.J. and Schindelegger, M., 2023. Scalable self attraction and loading calculations for unstructured ocean
!! tide models. Ocean Modelling, p.102160.
!! https://doi.org/10.1016/j.ocemod.2023.102160
!!
!! Schaeffer, N., 2013. Efficient spherical harmonic transforms aimed at pseudospectral numerical simulations.
!! Geochemistry, Geophysics, Geosystems, 14(3), pp.751-758.
!! https://doi.org/10.1002/ggge.20071
end module MOM_spherical_harmonics
