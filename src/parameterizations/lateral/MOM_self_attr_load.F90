module MOM_self_attr_load

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_MODULE
use MOM_domains, only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : slasher, MOM_read_data
use MOM_load_love_numbers, only : Love_Data
use MOM_obsolete_params, only : obsolete_logical, obsolete_int
use MOM_spherical_harmonics, only : spherical_harmonics_init, spherical_harmonics_end
use MOM_spherical_harmonics, only : spherical_harmonics_forward, spherical_harmonics_inverse
use MOM_spherical_harmonics, only : sht_CS, order2index, calc_lmax
use MOM_unit_scaling, only : unit_scale_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

public calc_SAL, scalar_SAL_sensitivity, SAL_init, SAL_end

#include <MOM_memory.h>

!> The control structure for the MOM_self_attr_load module
type, public :: SAL_CS ; private
  logical :: use_sal_scalar = .false.
    !< If true, use the scalar approximation to calculate SAL.
  logical :: use_sal_sht = .false.
    !< If true, use online spherical harmonics to calculate SAL
  logical :: use_tidal_sal_prev = .false.
    !< If true, read the tidal SAL from the previous iteration of the tides to
    !! facilitate convergence.
  logical :: use_bpa = .false.
    !< If true, use bottom pressure anomaly instead of SSH to calculate SAL.
  real :: eta_prop
    !< The partial derivative of eta_sal with the local value of eta [nondim].
  real :: linear_scaling
    !< Dimensional coefficients for scalar SAL [nondim or Z T2 L-2 R-1 ~> m Pa-1]
  type(sht_CS), allocatable :: sht
    !< Spherical harmonic transforms (SHT) control structure
  integer :: sal_sht_Nd
    !< Maximum degree for spherical harmonic transforms [nondim]
  real, allocatable :: ebot_ref(:,:)
    !< Reference bottom pressure scaled by Rho_0 and G_Earth[Z ~> m]
  real, allocatable :: Love_scaling(:)
    !< Dimensional coefficients for harmonic SAL, which are functions of Love numbers
    !! [nondim] or [Z T2 L-2 R-1 ~> m Pa-1], depending on the value of use_ppa.
  real, allocatable :: Snm_Re(:), &    !< Real SHT coefficient for SHT SAL [Z ~> m]
                       Snm_Im(:)       !< Imaginary SHT coefficient for SHT SAL [Z ~> m]
end type SAL_CS

integer :: id_clock_SAL   !< CPU clock for self-attraction and loading

contains

!> This subroutine calculates seawater self-attraction and loading based on either sea surface height (SSH) or bottom
!! pressure anomaly.  Note that the SAL calculation applies to all motions across the spectrum. Tidal-specific methods
!! that assume periodicity, i.e. iterative and read-in SAL, are stored in MOM_tidal_forcing module.
!! The input field can be either SSH [Z ~> m] or total bottom pressure [R L2 T-2 ~> Pa].  If total bottom pressure is
!! used, bottom pressure anomaly is first calculated by subtracting a reference bottom pressure from an input file.
!! The output field is expressed as geopotential height anomaly, and therefore has the unit of [Z ~> m].
subroutine calc_SAL(eta, eta_sal, G, CS, tmp_scale)
  type(ocean_grid_type), intent(in) :: G  !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: eta     !< The sea surface height anomaly from
              !! a time-mean geoid or total bottom pressure [Z ~> m] or [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: eta_sal !< The geopotential height anomaly from
              !! self-attraction and loading [Z ~> m].
  type(SAL_CS), intent(inout) :: CS !< The control structure returned by a previous call to SAL_init.
  real, optional, intent(in)  :: tmp_scale !< A rescaling factor to temporarily convert eta
              !! to MKS units in reproducing sumes [m Z-1 ~> 1]

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: bpa ! SSH or bottom pressure anomaly [Z ~> m] or [R L2 T-2 ~> Pa]
  integer :: n, m, l
  integer :: Isq, Ieq, Jsq, Jeq
  integer :: i, j

  call cpu_clock_begin(id_clock_SAL)

  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (CS%use_bpa) then ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    bpa(i,j) = eta(i,j) - CS%ebot_ref(i,j)
  enddo ; enddo ; else ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    bpa(i,j) = eta(i,j)
  enddo ; enddo ; endif

  ! use the scalar approximation and/or iterative tidal SAL
  if (CS%use_sal_scalar .or. CS%use_tidal_sal_prev) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_sal(i,j) = CS%linear_scaling * bpa(i,j)
    enddo ; enddo

  ! use the spherical harmonics method
  elseif (CS%use_sal_sht) then
    call spherical_harmonics_forward(G, CS%sht, bpa, CS%Snm_Re, CS%Snm_Im, CS%sal_sht_Nd, tmp_scale=tmp_scale)

    ! Multiply scaling factors to each mode
    do m = 0,CS%sal_sht_Nd
      l = order2index(m, CS%sal_sht_Nd)
      do n = m,CS%sal_sht_Nd
        CS%Snm_Re(l+n-m) = CS%Snm_Re(l+n-m) * CS%Love_scaling(l+n-m)
        CS%Snm_Im(l+n-m) = CS%Snm_Im(l+n-m) * CS%Love_scaling(l+n-m)
      enddo
    enddo

    call spherical_harmonics_inverse(G, CS%sht, CS%Snm_Re, CS%Snm_Im, eta_sal, CS%sal_sht_Nd)
    ! Halo was not calculated in spherical harmonic transforms.
    call pass_var(eta_sal, G%domain)

  else
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_sal(i,j) = 0.0
    enddo ; enddo
  endif

  call cpu_clock_end(id_clock_SAL)
end subroutine calc_SAL

!>   This subroutine returns eta_prop member of SAL_CS type, which is the non-dimensional partial
!! derivative of the local geopotential height with the input sea surface height due to the scalar
!! approximation of self-attraction and loading.
subroutine scalar_SAL_sensitivity(CS, deta_sal_deta)
  type(SAL_CS), intent(in)  :: CS !< The control structure returned by a previous call to SAL_init.
  real,         intent(out) :: deta_sal_deta !< The partial derivative of eta_sal with
                                             !! the local value of eta [nondim].
  deta_sal_deta = CS%eta_prop
end subroutine scalar_SAL_sensitivity

!> This subroutine calculates coefficients of the spherical harmonic modes for self-attraction and loading.
!! The algorithm is based on the SAL implementation in MPAS-ocean, which was modified by Kristin Barton from
!! routine written by K. Quinn (March 2010) and modified by M. Schindelegger (May 2017).
subroutine calc_love_scaling(rhoW, rhoE, grav, CS)
  real, intent(in) :: rhoW !< The average density of sea water [R ~> kg m-3]
  real, intent(in) :: rhoE !< The average density of Earth [R ~> kg m-3]
  real, intent(in) :: grav !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  type(SAL_CS), intent(inout) :: CS !< The control structure returned by a previous call to SAL_init.

  ! Local variables
  real :: coef_rhoE ! A scaling coefficient of solid Earth density. coef_rhoE = rhoW / rhoE with USE_BPA=False
      ! and coef_rhoE = 1.0 / (rhoE * grav) with USE_BPA=True. [nondim] or [Z T2 L-2 R-1 ~> m Pa-1]
  real, dimension(:), allocatable :: HDat, LDat, KDat ! Love numbers converted in CF reference frames [nondim]
  real :: H1, L1, K1 ! Temporary variables to store degree 1 Love numbers [nondim]
  integer :: n_tot ! Size of the stored Love numbers [nondim]
  integer :: nlm  ! Maximum spherical harmonics degree [nondim]
  integer :: n, m, l

  n_tot = size(Love_Data, dim=2)
  nlm = CS%sal_sht_Nd

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

  if (CS%use_bpa) then
    coef_rhoE = 1.0 / (rhoE * grav) ! [Z T2 L-2 R-1 ~> m Pa-1]
  else
    coef_rhoE = rhoW / rhoE ! [nondim]
  endif

  do m=0,nlm ; do n=m,nlm
    l = order2index(m, nlm)
    ! Love_scaling has the same as coef_rhoE.
    CS%Love_scaling(l+n-m) = (3.0 / real(2*n+1)) * coef_rhoE * (1.0 + KDat(n+1) - HDat(n+1))
  enddo ; enddo
end subroutine calc_love_scaling

!> This subroutine initializes the self-attraction and loading control structure.
subroutine SAL_init(G, GV, US, param_file, CS)
  type(ocean_grid_type),   intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(SAL_CS), intent(inout) :: CS   !< Self-attraction and loading control structure

  ! Local variables
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_self_attr_load" ! This module's name.
  integer :: lmax ! Total modes of the real spherical harmonics [nondim]
  real :: rhoE    ! The average density of Earth [R ~> kg m-3].
  character(len=200) :: filename, ebot_ref_file, inputdir ! Strings for file/path
  character(len=200) :: ebot_ref_varname                  ! Variable name in file
  logical :: calculate_sal, tides, use_tidal_sal_file
  integer :: tides_answer_date ! Recover old answers with tides
  real :: sal_scalar_value ! Scaling SAL factors [nondim]
  integer :: isd, ied, jsd, jed

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, '', "TIDES", tides, default=.false., do_not_log=.True.)
  call get_param(param_file, '', "CALCULATE_SAL", calculate_sal, default=tides, do_not_log=.True.)
  if (.not. calculate_sal) return

  if (tides) then
    call get_param(param_file, '', "USE_PREVIOUS_TIDES", CS%use_tidal_sal_prev, &
                   default=.false., do_not_log=.True.)
    call get_param(param_file, '', "TIDAL_SAL_FROM_FILE", use_tidal_sal_file, &
                   default=.false., do_not_log=.True.)
    call get_param(param_file, '', "TIDES_ANSWER_DATE", tides_answer_date, &
                   default=20230630, do_not_log=.True.)
  endif

  call get_param(param_file, mdl, "SAL_USE_BPA", CS%use_bpa, &
                 "If true, use bottom pressure anomaly to calculate self-attraction and "// &
                 "loading (SAL). Otherwise sea surface height anomaly is used, which is "// &
                 "only correct for homogenous flow.", default=.False.)
  if (CS%use_bpa) then
    call get_param(param_file, '', "INPUTDIR", inputdir, default=".", do_not_log=.True.)
    inputdir = slasher(inputdir)
    call get_param(param_file, mdl, "REF_BOT_PRES_FILE", ebot_ref_file, &
                   "Reference bottom pressure file used by self-attraction and loading (SAL).", &
                   default="pbot.nc")
    call get_param(param_file, mdl, "REF_BOT_PRES_VARNAME", ebot_ref_varname, &
                   "The name of the variable in REF_BOT_PRES_FILE with reference bottom "//&
                   "pressure.  The variable should have the unit of Pa.", &
                   default="pbot")
    filename = trim(inputdir)//trim(ebot_ref_file)
    call log_param(param_file, mdl, "INPUTDIR/REF_BOT_PRES_FILE", filename)

    allocate(CS%ebot_ref(isd:ied, jsd:jed), source=0.0)
    call MOM_read_data(filename, trim(ebot_ref_varname), CS%ebot_ref, G%Domain,&
                       scale=US%Pa_to_RL2_T2)
    call pass_var(CS%ebot_ref, G%Domain)
  endif
  if (tides_answer_date<=20250131 .and. CS%use_bpa) &
    call MOM_error(FATAL, trim(mdl) // ", SAL_init: SAL_USE_BPA needs to be false to recover "//&
                   "tide answers before 20250131.")
  call get_param(param_file, mdl, "SAL_SCALAR_APPROX", CS%use_sal_scalar, &
                 "If true, use the scalar approximation to calculate self-attraction and "//&
                 "loading.", default=tides .and. (.not. use_tidal_sal_file))
  if (CS%use_sal_scalar .and. CS%use_bpa) &
    call MOM_error(WARNING, trim(mdl) // ", SAL_init: Using bottom pressure anomaly for scalar "//&
                   "approximation SAL is unsubstantiated.")
  call get_param(param_file, mdl, "SAL_SCALAR_VALUE", sal_scalar_value, "The constant of "//&
                 "proportionality between self-attraction and loading (SAL) geopotential "//&
                 "anomaly and barotropic geopotential anomaly. This is only used if "//&
                 "SAL_SCALAR_APPROX is true or USE_PREVIOUS_TIDES is true.", default=0.0, &
                 units="m m-1", do_not_log=.not.(CS%use_sal_scalar .or. CS%use_tidal_sal_prev), &
                 old_name='TIDE_SAL_SCALAR_VALUE')
  call get_param(param_file, mdl, "SAL_HARMONICS", CS%use_sal_sht, &
                 "If true, use the online spherical harmonics method to calculate "//&
                 "self-attraction and loading.", default=.false.)
  call get_param(param_file, mdl, "SAL_HARMONICS_DEGREE", CS%sal_sht_Nd, &
                 "The maximum degree of the spherical harmonics transformation used for "// &
                 "calculating the self-attraction and loading term.", &
                 default=0, do_not_log=(.not. CS%use_sal_sht))
  call get_param(param_file, mdl, "RHO_SOLID_EARTH", rhoE, &
                 "The mean solid earth density.  This is used for calculating the "// &
                 "self-attraction and loading term.", units="kg m-3", &
                 default=5517.0, scale=US%kg_m3_to_R, do_not_log=(.not. CS%use_sal_sht))

  ! Set scaling coefficients for scalar approximation
  if (CS%use_sal_scalar .or. CS%use_tidal_sal_prev) then
    if (CS%use_sal_scalar .and. CS%use_tidal_sal_prev) then
      CS%eta_prop = 2.0 * sal_scalar_value
    else
      CS%eta_prop = sal_scalar_value
    endif
    if (CS%use_bpa) then
      CS%linear_scaling = CS%eta_prop / (GV%Rho0 * GV%g_Earth)
    else
      CS%linear_scaling = CS%eta_prop
    endif
  else
    CS%eta_prop = 0.0 ; CS%linear_scaling = 0.0
  endif

  ! Set scaling coefficients for spherical harmonics
  if (CS%use_sal_sht) then
    lmax = calc_lmax(CS%sal_sht_Nd)
    allocate(CS%Snm_Re(lmax), source=0.0)
    allocate(CS%Snm_Im(lmax), source=0.0)

    allocate(CS%Love_scaling(lmax), source=0.0)
    call calc_love_scaling(GV%Rho0, rhoE, GV%g_Earth, CS)

    allocate(CS%sht)
    call spherical_harmonics_init(G, param_file, CS%sht)
  endif

  id_clock_SAL = cpu_clock_id('(Ocean SAL)', grain=CLOCK_MODULE)

end subroutine SAL_init

!> This subroutine deallocates memory associated with the SAL module.
subroutine SAL_end(CS)
  type(SAL_CS), intent(inout) :: CS !< The control structure returned by a previous call
                                    !! to SAL_init; it is deallocated here.

  if (allocated(CS%ebot_ref)) deallocate(CS%ebot_ref)

  if (CS%use_sal_sht) then
    if (allocated(CS%Love_scaling)) deallocate(CS%Love_scaling)
    if (allocated(CS%Snm_Re)) deallocate(CS%Snm_Re)
    if (allocated(CS%Snm_Im)) deallocate(CS%Snm_Im)
    call spherical_harmonics_end(CS%sht)
    deallocate(CS%sht)
  endif
end subroutine SAL_end

!> \namespace self_attr_load
!!
!! \section section_SAL Self attraction and loading
!!
!! This module contains methods to calculate self-attraction and loading (SAL) as a function of sea surface height or
!! bottom pressure anomaly. SAL is primarily used for fast evolving processes like tides or storm surges, but the
!! effect applies to all motions.
!!
!! If <code>SAL_SCALAR_APPROX</code> is true, a scalar approximation is applied (\cite Accad1978) and the SAL is simply
!! a fraction (set by <code>SAL_SCALAR_VALUE</code>, usually around 10% for global tides) of local SSH. For tides, the
!! scalar approximation can also be used to iterate the SAL to convergence [see <code>USE_PREVIOUS_TIDES</code> in
!! MOM_tidal_forcing, \cite Arbic2004].
!!
!! If <code>SAL_HARMONICS</code> is true, a more accurate online spherical harmonic transforms are used to calculate
!! SAL. Subroutines in module MOM_spherical_harmonics are called and the degree of spherical harmonic transforms is set
!! by <code>SAL_HARMONICS_DEGREE</code>. The algorithm is based on SAL calculation in Model for Prediction Across
!! Scales (MPAS)-Ocean developed by Los Alamos National Laboratory and University of Michigan
!! [\cite Barton2022 and \cite Brus2023].
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
