module MOM_self_attr_load

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_MODULE
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_spherical_harmonics, only : spherical_harmonics_init, spherical_harmonics_end, order2index, calc_lmax
use MOM_spherical_harmonics, only : spherical_harmonics_forward, spherical_harmonics_inverse
use MOM_spherical_harmonics, only : sht_CS
use MOM_load_love_numbers, only : Love_Data

implicit none ; private

public calc_SAL, scalar_SAL_sensitivity, SAL_init, SAL_end

#include <MOM_memory.h>

!> The control structure for the MOM_self_attr_load module
type, public :: SAL_CS ; private
  logical :: use_sal_scalar !< If true, use the scalar approximation when
                      !! calculating self-attraction and loading.
  real    :: sal_scalar !< The constant of proportionality between sea surface
                      !! height (really it should be bottom pressure) anomalies
                      !! and bottom geopotential anomalies [nondim].
  logical :: use_prev_tides !< If true, use the SAL from the previous iteration of the tides
                            !! to facilitate convergence.
  logical :: use_sal_sht !< If true, use online spherical harmonics to calculate SAL
  type(sht_CS) :: sht !< Spherical harmonic transforms (SHT) for SAL
  integer :: sal_sht_Nd !< Maximum degree for SHT [nodim]
  real, allocatable :: Love_Scaling(:) !< Love number for each SHT mode [nodim]
  real, allocatable :: Snm_Re(:), & !< Real and imaginary SHT coefficient for SHT SAL
                       Snm_Im(:)    !< [Z ~> m]
end type SAL_CS

integer :: id_clock_SAL   !< CPU clock for self-attraction and loading

contains

!> This subroutine calculates seawater self-attraction and loading based on sea surface height. This should
!! be changed into bottom pressure anomaly in the future. Note that the SAL calculation applies to all motions
!! across the spectrum. Tidal-specific methods that assume periodicity, i.e. iterative and read-in SAL, are
!! stored in MOM_tidal_forcing module.
subroutine calc_SAL(eta, eta_sal, G, CS)
  type(ocean_grid_type), intent(in)  :: G  !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: eta     !< The sea surface height anomaly from
                                                           !! a time-mean geoid [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: eta_sal !< The sea surface height anomaly from
                                                           !! self-attraction and loading [Z ~> m].
  type(SAL_CS), intent(inout) :: CS !< The control structure returned by a previous call to SAL_init.

  ! Local variables
  integer :: n, m, l
  integer :: Isq, Ieq, Jsq, Jeq
  integer :: i, j
  real :: eta_prop

  call cpu_clock_begin(id_clock_SAL)

  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  ! use the scalar approximation, iterative tidal SAL or no SAL
  call scalar_SAL_sensitivity(CS, eta_prop)
  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    eta_sal(i,j) = eta_prop*eta(i,j)
  enddo ; enddo

  if (CS%use_sal_sht) then ! use the spherical harmonics method
    call spherical_harmonics_forward(G, CS%sht, eta, CS%Snm_Re, CS%Snm_Im, CS%sal_sht_Nd)

    ! Multiply scaling factors to each mode
    do m = 0,CS%sal_sht_Nd
      l = order2index(m, CS%sal_sht_Nd)
      do n = m,CS%sal_sht_Nd
        CS%Snm_Re(l+n-m) = CS%Snm_Re(l+n-m) * CS%Love_Scaling(l+n-m)
        CS%Snm_Im(l+n-m) = CS%Snm_Im(l+n-m) * CS%Love_Scaling(l+n-m)
      enddo
    enddo

    call spherical_harmonics_inverse(G, CS%sht, CS%Snm_Re, CS%Snm_Im, eta_sal, CS%sal_sht_Nd)

    call pass_var(eta_sal, G%domain)
  endif

  call cpu_clock_end(id_clock_SAL)
end subroutine calc_SAL

!>   This subroutine calculates the partial derivative of the local geopotential height with the input
!! sea surface height due to the scalar approximation of self-attraction and loading.
subroutine scalar_SAL_sensitivity(CS, deta_sal_deta)
  type(SAL_CS), intent(in)  :: CS !< The control structure returned by a previous call to SAL_init.
  real,         intent(out) :: deta_sal_deta !< The partial derivative of eta_sal with
                                             !! the local value of eta [nondim].

  if (CS%USE_SAL_SCALAR .and. CS%USE_PREV_TIDES) then
    deta_sal_deta = 2.0*CS%SAL_SCALAR
  elseif (CS%USE_SAL_SCALAR .or. CS%USE_PREV_TIDES) then
    deta_sal_deta = CS%SAL_SCALAR
  else
    deta_sal_deta = 0.0
  endif
end subroutine scalar_SAL_sensitivity

!> This subroutine calculates coefficients of the spherical harmonic modes for self-attraction and loading.
!! The algorithm is based on the SAL implementation in MPAS-ocean, which was modified by Kristin Barton from
!! routine written by K. Quinn (March 2010) and modified by M. Schindelegger (May 2017).
subroutine calc_love_scaling(nlm, rhoW, rhoE, Love_Scaling)
  integer, intent(in) :: nlm  !< Maximum spherical harmonics degree [nondim]
  real,    intent(in) :: rhoW !< The average density of sea water [R ~> kg m-3]
  real,    intent(in) :: rhoE !< The average density of Earth [R ~> kg m-3]
  real, dimension(:), intent(out) :: Love_Scaling !< Scaling factors for inverse SHT [nondim]

  ! Local variables
  real, dimension(:), allocatable :: HDat, LDat, KDat ! Love numbers converted in CF reference frames
  real :: H1, L1, K1 ! Temporary variables to store degree 1 Love numbers
  integer :: n_tot ! Size of the stored Love numbers
  integer :: n, m, l

  n_tot = size(Love_Data, dim=2)

  if (nlm+1 > n_tot) call MOM_error(FATAL, "MOM_tidal_forcing " // &
    "calc_love_scaling: maximum spherical harmonics degree is larger than " // &
    "the size of the stored Love numbers in MOM_load_love_number.")

  allocate(HDat(nlm+1), LDat(nlm+1), KDat(nlm+1))
  HDat(:) = Love_Data(2,1:nlm+1) ; LDat(:) = Love_Data(3,1:nlm+1) ; KDat(:) = Love_Data(4,1:nlm+1)

  ! Convert reference frames from CM to CF
  if (nlm > 0) then
    H1 = HDat(2) ; L1 = LDat(2) ;  K1 = KDat(2)
    HDat(2) = ( 2.0 / 3.0) * (H1 - L1)
    LDat(2) = (-1.0 / 3.0) * (H1 - L1)
    KDat(2) = (-1.0 / 3.0) * H1 - (2.0 / 3.0) * L1 - 1.0
  endif

  do m=0,nlm ; do n=m,nlm
    l = order2index(m,nlm)
    Love_Scaling(l+n-m) = (3.0 / real(2*n+1)) * (rhoW / rhoE) * (1.0 + KDat(n+1) - HDat(n+1))
  enddo ; enddo
end subroutine calc_love_scaling

!> This subroutine initializeds the self-attraction and loading control structure.
subroutine SAL_init(G, US, param_file, CS)
  type(ocean_grid_type),  intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type),  intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(SAL_CS), intent(inout) :: CS   !< Self-attraction and loading control structure

# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_self_attr_load" ! This module's name.
  integer :: lmax ! Total modes of the real spherical harmonics [nondim]
  real :: rhoW    ! The average density of sea water [R ~> kg m-3].
  real :: rhoE    ! The average density of Earth [R ~> kg m-3].

  logical :: calculate_sal
  logical :: tides, tidal_sal_from_file

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, '', "TIDES", tides, default=.false., do_not_log=.True.)

  CS%use_prev_tides = .false.
  tidal_sal_from_file = .false.
  if (tides) then
    call get_param(param_file, '', "USE_PREVIOUS_TIDES", CS%use_prev_tides,&
                   default=.false., do_not_log=.True.)
    call get_param(param_file, '', "TIDAL_SAL_FROM_FILE", tidal_sal_from_file,&
                   default=.false., do_not_log=.True.)
  endif

  call get_param(param_file, mdl, "TIDE_USE_SAL_SCALAR", CS%use_sal_scalar, &
                "If true and TIDES is true, use the scalar approximation "//&
                "when calculating self-attraction and loading.", &
                default=.not.tidal_sal_from_file)
  if (CS%use_sal_scalar .or. CS%use_prev_tides) &
    call get_param(param_file, mdl, "TIDE_SAL_SCALAR_VALUE", CS%sal_scalar, &
                  "The constant of proportionality between sea surface "//&
                  "height (really it should be bottom pressure) anomalies "//&
                  "and bottom geopotential anomalies. This is only used if "//&
                  "TIDES and TIDE_USE_SAL_SCALAR are true.", units="m m-1", &
                  fail_if_missing=.true.)

  call get_param(param_file, mdl, "TIDAL_SAL_SHT", CS%use_sal_sht, &
                 "If true, use the online spherical harmonics method to calculate "//&
                 "self-attraction and loading term in tides.", default=.false.)

  call get_param(param_file, mdl, "CALCULATE_SAL", calculate_sal, &
                 "If true, calculate self-attraction and loading.", default=tides)

  ! ! Default USE_SAL is TRUE for now to keep backward compatibility with old MOM_INPUT files. It should be changed to
  ! ! FALSE in the future (mostly to avoid the SSH calculations in MOM_PressureForce). In that case, the following check
  ! ! informs prior tidal experiments that use scalar or iterative SAL to include USE_SAL flag, as the USE_SAL flag
  ! ! overrules the option flags.
  ! if ((.not. calculate_sal) .and. (CS%use_prev_tides .or. CS%use_sal_scalar .or. CS%use_sal_sht)) &
  !   call MOM_error(FATAL, trim(mdl)//": USE_SAL is False but one of the options is True. Nothing will happen.")

  if (CS%use_sal_sht) then
    call get_param(param_file, mdl, "TIDAL_SAL_SHT_DEGREE", CS%sal_sht_Nd, &
                   "The maximum degree of the spherical harmonics transformation used for "// &
                   "calculating the self-attraction and loading term.", &
                   default=0, do_not_log=.not. CS%use_sal_sht)
    call get_param(param_file, mdl, "RHO_0", rhoW, default=1035.0, scale=US%kg_m3_to_R, do_not_log=.True.)
    call get_param(param_file, mdl, "RHO_E", rhoE, &
                   "The mean solid earth density.  This is used for calculating the "// &
                   "self-attraction and loading term.", units="kg m-3", &
                   default=5517.0, scale=US%kg_m3_to_R, do_not_log=.not. CS%use_sal_sht)
    lmax = calc_lmax(CS%sal_sht_Nd)
    allocate(CS%Snm_Re(lmax)); CS%Snm_Re(:) = 0.0
    allocate(CS%Snm_Im(lmax)); CS%Snm_Im(:) = 0.0

    allocate(CS%Love_Scaling(lmax)); CS%Love_Scaling(:) = 0.0
    call calc_love_scaling(CS%sal_sht_Nd, rhoW, rhoE, CS%Love_Scaling)
    call spherical_harmonics_init(G, param_file, CS%sht)
  endif

  id_clock_SAL = cpu_clock_id('(Ocean SAL)', grain=CLOCK_MODULE)

end subroutine SAL_init

!> This subroutine deallocates memory associated with the SAL module.
subroutine SAL_end(CS)
  type(SAL_CS), intent(inout) :: CS !< The control structure returned by a previous call
                                    !! to SAL_init; it is deallocated here.
  if (CS%use_sal_sht) then
    if (allocated(CS%Love_Scaling)) deallocate(CS%Love_Scaling)
    if (allocated(CS%Snm_Re)) deallocate(CS%Snm_Re)
    if (allocated(CS%Snm_Im)) deallocate(CS%Snm_Im)
    call spherical_harmonics_end(CS%sht)
  endif
end subroutine SAL_end

!> \namespace self_attr_load
!!
!! This module contains methods to calculate self-attraction and loading (SAL) as a function of sea surface height (SSH)
!! (rather, it should be bottom pressure anomaly). SAL is primarily used for fast evolving processes like tides or
!! storm surges, but the effect applys to all motions.
!!
!!     If TIDE_USE_SAL_SCALAR is true, a scalar approximiation is applied (Accad and Pekeris 1978) and the SAL is simply
!! a fraction (set by TIDE_SAL_SCALAR_VALUE, usualy around 10% for global tides) of local SSH . For the tides, the
!! scalar approximation can also be used to iterate the SAL to convergence [see USE_PREVIOUS_TIDES in MOM_tidal_forcing,
!! Arbic et al. (2004)].
!!
!!    If TIDAL_SAL_SHT is true, a more accurate online spherical harmonic transforms are used to calculate SAL.
!! Subroutines in module MOM_spherical_harmonics are called and the degree of spherical harmonic transforms is set by
!! TIDAL_SAL_SHT_DEGREE. The algorithm is based on SAL calculation in Model for Prediction Across Scales (MPAS)-Ocean
!! developed by Los Alamos National Laboratory and University of Michigan [Barton et al. (2022) and Brus et al. (2023)].
!!
!! References:
!!
!! Accad, Y. and Pekeris, C.L., 1978. Solution of the tidal equations for the M2 and S2 tides in the world oceans from a
!! knowledge of the tidal potential alone. Philosophical Transactions of the Royal Society of London. Series A,
!! Mathematical and Physical Sciences, 290(1368), pp.235-266.
!! https://doi.org/10.1098/rsta.1978.0083
!!
!! Arbic, B.K., Garner, S.T., Hallberg, R.W. and Simmons, H.L., 2004. The accuracy of surface elevations in forward
!! global barotropic and baroclinic tide models. Deep Sea Research Part II: Topical Studies in Oceanography, 51(25-26),
!! pp.3069-3101.
!! https://doi.org/10.1016/j.dsr2.2004.09.014
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
end module MOM_self_attr_load
