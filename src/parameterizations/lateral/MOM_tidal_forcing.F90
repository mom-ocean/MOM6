!> Tidal contributions to geopotential
module MOM_tidal_forcing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
                              CLOCK_MODULE
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : field_exists, file_exists, MOM_read_data
use MOM_time_manager,  only : set_date, time_type, time_type_to_real, operator(-)
use MOM_unit_scaling,  only : unit_scale_type
use MOM_spherical_harmonics, only : spherical_harmonics_init, spherical_harmonics_end
use MOM_spherical_harmonics, only : associatedLegendrePolynomials, SHOrderDegreeToIndex
use MOM_spherical_harmonics, only : sht_CS
use MOM_coms_infra, only : sum_across_PEs

implicit none ; private

public calc_tidal_forcing, tidal_forcing_init, tidal_forcing_end
public tidal_forcing_sensitivity
! MOM_open_boundary uses the following to set tides on the boundary.
public astro_longitudes_init, eq_phase, nodal_fu, tidal_frequency

#include <MOM_memory.h>

integer, parameter :: MAX_CONSTITUENTS = 10 !< The maximum number of tidal
                                            !! constituents that could be used.
!> Simple type to store astronomical longitudes used to calculate tidal phases.
type, public :: astro_longitudes
  real :: &
    s, &  !< Mean longitude of moon [rad]
    h, &  !< Mean longitude of sun [rad]
    p, &  !< Mean longitude of lunar perigee [rad]
    N     !< Longitude of ascending node [rad]
end type astro_longitudes

!> The control structure for the MOM_tidal_forcing module
type, public :: tidal_forcing_CS ; private
  logical :: use_sal_scalar !< If true, use the scalar approximation when
                      !! calculating self-attraction and loading.
  logical :: tidal_sal_from_file !< If true, Read the tidal self-attraction
                      !! and loading from input files, specified
                      !! by TIDAL_INPUT_FILE.
  logical :: use_prev_tides !< If true, use the SAL from the previous
                      !! iteration of the tides to facilitate convergence.
  logical :: use_eq_phase !< If true, tidal forcing is phase-shifted to match
                      !! equilibrium tide. Set to false if providing tidal phases
                      !! that have already been shifted by the
                      !! astronomical/equilibrium argument.
  logical :: tidal_sal_sht
  real    :: sal_scalar !< The constant of proportionality between sea surface
                      !! height (really it should be bottom pressure) anomalies
                      !! and bottom geopotential anomalies [nondim].
  integer :: nc       !< The number of tidal constituents in use.
  real, dimension(MAX_CONSTITUENTS) :: &
    freq, &           !< The frequency of a tidal constituent [T-1 ~> s-1].
    phase0, &         !< The phase of a tidal constituent at time 0 [rad].
    amp, &            !< The amplitude of a tidal constituent at time 0 [Z ~> m].
    love_no           !< The Love number of a tidal constituent at time 0 [nondim].
  integer :: struct(MAX_CONSTITUENTS) !< An encoded spatial structure for each constituent
  character (len=16) :: const_name(MAX_CONSTITUENTS) !< The name of each constituent

  type(time_type) :: time_ref !< Reference time (t = 0) used to calculate tidal forcing.
  type(astro_longitudes) :: tidal_longitudes !< Astronomical longitudes used to calculate
                                   !! tidal phases at t = 0.
  real, allocatable :: &
    sin_struct(:,:,:), &    !< The sine and cosine based structures that can
    cos_struct(:,:,:), &    !< be associated with the astronomical forcing [nondim].
    cosphasesal(:,:,:), &   !< The cosine and sine of the phase of the
    sinphasesal(:,:,:), &   !< self-attraction and loading amphidromes.
    ampsal(:,:,:), &        !< The amplitude of the SAL [Z ~> m].
    cosphase_prev(:,:,:), & !< The cosine and sine of the phase of the
    sinphase_prev(:,:,:), & !< amphidromes in the previous tidal solutions.
    amp_prev(:,:,:)         !< The amplitude of the previous tidal solution [Z ~> m].
  type(sht_CS) :: sht
end type tidal_forcing_CS

integer :: id_clock_tides !< CPU clock for tides
integer :: id_clock_SAL   !< CPU clock for inline self-attration and loading

contains

!> Finds astronomical longitudes s, h, p, and N,
!! the mean longitude of the moon, sun, lunar perigee, and ascending node, respectively,
!! at the specified reference time time_ref.
!! These formulas were obtained from
!! Kowalik and Luick, "Modern Theory and Practice of Tide Analysis and Tidal Power", 2019
!! (their Equation I.71), which are based on Schureman, 1958.
!! For simplicity, the time associated with time_ref should
!! be at midnight. These formulas also only make sense if
!! the calendar is gregorian.
subroutine astro_longitudes_init(time_ref, longitudes)
  type(time_type), intent(in) :: time_ref            !> Time to calculate longitudes for.
  type(astro_longitudes), intent(out) :: longitudes  !> Lunar and solar longitudes at time_ref.
  real :: D                                          !> Time since the reference date [days]
  real :: T                                          !> Time in Julian centuries [centuries]
  real, parameter :: PI = 4.0 * atan(1.0)            !> 3.14159... [nondim]
  ! Find date at time_ref in days since 1900-01-01
  D = time_type_to_real(time_ref - set_date(1900, 1, 1)) / (24.0 * 3600.0)
  ! Time since 1900-01-01 in Julian centuries
  ! Kowalik and Luick use 36526, but Schureman uses 36525 which I think is correct.
  T = D / 36525.0
  ! Calculate longitudes, including converting to radians on [0, 2pi)
  ! s: Mean longitude of moon
  longitudes%s = mod((277.0248 + 481267.8906 * T) + 0.0011 * (T**2), 360.0) * PI / 180.0
  ! h: Mean longitude of sun
  longitudes%h = mod((280.1895 + 36000.7689 * T) + 3.0310e-4 * (T**2), 360.0) * PI / 180.0
  ! p: Mean longitude of lunar perigee
  longitudes%p = mod((334.3853 + 4069.0340 * T) - 0.0103 * (T**2), 360.0) * PI / 180.0
  ! n: Longitude of ascending node
  longitudes%N = mod((259.1568 - 1934.142 * T) + 0.0021 * (T**2), 360.0) * PI / 180.0
end subroutine astro_longitudes_init

!> Calculates the equilibrium phase argument for the given tidal
!! constituent constit and the astronomical longitudes and the reference time.
!! These formulas follow Table I.4 of Kowalik and Luick,
!! "Modern Theory and Practice of Tide Analysis and Tidal Power", 2019.
function eq_phase(constit, longitudes)
  character (len=2), intent(in) :: constit !> Name of constituent (e.g., M2).
  type(astro_longitudes), intent(in) :: longitudes   !> Mean longitudes calculated using astro_longitudes_init
  real, parameter :: PI = 4.0 * atan(1.0)  !> 3.14159...
  real :: eq_phase                         !> The equilibrium phase argument for the constituent [rad].

  select case (constit)
    case ("M2")
      eq_phase = 2 * (longitudes%h - longitudes%s)
    case ("S2")
      eq_phase = 0.0
    case ("N2")
      eq_phase = (- 3 * longitudes%s + 2 * longitudes%h) + longitudes%p
    case ("K2")
      eq_phase = 2 * longitudes%h
    case ("K1")
      eq_phase = longitudes%h + PI / 2.0
    case ("O1")
      eq_phase = (- 2 * longitudes%s + longitudes%h) - PI / 2.0
    case ("P1")
      eq_phase = - longitudes%h - PI / 2.0
    case ("Q1")
      eq_phase = ((- 3 * longitudes%s + longitudes%h) + longitudes%p) - PI / 2.0
    case ("MF")
      eq_phase = 2 * longitudes%s
    case ("MM")
      eq_phase = longitudes%s - longitudes%p
    case default
      call MOM_error(FATAL, "eq_phase: unrecognized constituent")
  end select
end function eq_phase

!> Looks up angular frequencies for the main tidal constituents.
!! Values used here are from previous versions of MOM.
function tidal_frequency(constit)
  character (len=2), intent(in) :: constit !> Constituent to look up
  real :: tidal_frequency                  !> Angular frequency [s-1]

  select case (constit)
    case ("M2")
      tidal_frequency = 1.4051890e-4
    case ("S2")
      tidal_frequency = 1.4544410e-4
    case ("N2")
      tidal_frequency = 1.3787970e-4
    case ("K2")
      tidal_frequency = 1.4584234e-4
    case ("K1")
      tidal_frequency = 0.7292117e-4
    case ("O1")
      tidal_frequency = 0.6759774e-4
    case ("P1")
      tidal_frequency = 0.7252295e-4
    case ("Q1")
      tidal_frequency = 0.6495854e-4
    case ("MF")
      tidal_frequency = 0.053234e-4
    case ("MM")
      tidal_frequency = 0.026392e-4
    case default
      call MOM_error(FATAL, "tidal_frequency: unrecognized constituent")
  end select
end function tidal_frequency

!> Find amplitude (f) and phase (u) modulation of tidal constituents by the 18.6
!! year nodal cycle. Values here follow Table I.6 in Kowalik and Luick,
!! "Modern Theory and Practice of Tide Analysis and Tidal Power", 2019.
subroutine nodal_fu(constit, nodelon, fn, un)
  character (len=2), intent(in)  :: constit !> Tidal constituent to find modulation for.
  real,              intent(in)  :: nodelon !> Longitude of ascending node [rad], which
                                            !! can be calculated using astro_longitudes_init.
  real,              intent(out) :: fn      !> Amplitude modulation [nondim]
  real,              intent(out) :: un      !> Phase modulation [rad]

  real, parameter :: RADIANS = 4.0 * atan(1.0) / 180.0  !> Converts degrees to radians [nondim]

  select case (constit)
    case ("M2")
      fn = 1.0 - 0.037 * cos(nodelon)
      un = -2.1 * RADIANS * sin(nodelon)
    case ("S2")
      fn = 1.0  ! Solar S2 has no amplitude modulation.
      un = 0.0  ! S2 has no phase modulation.
    case ("N2")
      fn = 1.0 - 0.037 * cos(nodelon)
      un = -2.1 * RADIANS * sin(nodelon)
    case ("K2")
      fn = 1.024 + 0.286 * cos(nodelon)
      un = -17.7 * RADIANS * sin(nodelon)
    case ("K1")
      fn = 1.006 + 0.115 * cos(nodelon)
      un = -8.9 * RADIANS * sin(nodelon)
    case ("O1")
      fn = 1.009 + 0.187 * cos(nodelon)
      un = 10.8 * RADIANS * sin(nodelon)
    case ("P1")
      fn = 1.0  ! P1 has no amplitude modulation.
      un = 0.0  ! P1 has no phase modulation.
    case ("Q1")
      fn = 1.009 + 0.187 * cos(nodelon)
      un = 10.8 * RADIANS * sin(nodelon)
    case ("MF")
      fn = 1.043 + 0.414 * cos(nodelon)
      un = -23.7 * RADIANS * sin(nodelon)
    case ("MM")
      fn = 1.0 - 0.130 * cos(nodelon)
      un = 0.0  ! MM has no phase modulation.
    case default
      call MOM_error(FATAL, "nodal_fu: unrecognized constituent")
  end select

end subroutine nodal_fu

!> This subroutine allocates space for the static variables used
!! by this module.  The metrics may be effectively 0, 1, or 2-D arrays,
!! while fields like the background viscosities are 2-D arrays.
!! ALLOC is a macro defined in MOM_memory.h for allocate or nothing with
!! static memory.
subroutine tidal_forcing_init(Time, G, US, param_file, CS)
  type(time_type),        intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),  intent(inout) :: G    !< The ocean's grid structure.
  type(unit_scale_type),  intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(tidal_forcing_CS), intent(inout) :: CS   !< Tidal forcing control struct

  ! Local variables
  real, dimension(SZI_(G), SZJ_(G)) :: &
    phase, &          ! The phase of some tidal constituent.
    lat_rad, lon_rad  ! Latitudes and longitudes of h-points in radians.
  real :: deg_to_rad
  real, dimension(MAX_CONSTITUENTS) :: freq_def ! Default frequency for each tidal constituent [s-1]
  real, dimension(MAX_CONSTITUENTS) :: phase0_def ! Default reference phase for each tidal constituent [rad]
  real, dimension(MAX_CONSTITUENTS) :: amp_def  ! Default amplitude for each tidal constituent [m]
  real, dimension(MAX_CONSTITUENTS) :: love_def ! Default love number for each constituent [nondim]
  integer, dimension(3) :: tide_ref_date !< Reference date (t = 0) for tidal forcing.
  logical :: use_M2, use_S2, use_N2, use_K2, use_K1, use_O1, use_P1, use_Q1
  logical :: use_MF, use_MM
  logical :: tides      ! True if a tidal forcing is to be used.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_tidal_forcing" ! This module's name.
  character(len=128) :: mesg
  character(len=200) :: tidal_input_files(4*MAX_CONSTITUENTS)
  integer :: i, j, c, is, ie, js, je, isd, ied, jsd, jed, nc
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd; jed = G%jed

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "TIDES", tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)

  if (.not.tides) return

  ! Set up the spatial structure functions for the diurnal, semidiurnal, and
  ! low-frequency tidal components.
  allocate(CS%sin_struct(isd:ied,jsd:jed,3), source=0.0)
  allocate(CS%cos_struct(isd:ied,jsd:jed,3), source=0.0)
  deg_to_rad = 4.0*ATAN(1.0)/180.0
  do j=js-1,je+1 ; do i=is-1,ie+1
    lat_rad(i,j) = G%geoLatT(i,j)*deg_to_rad
    lon_rad(i,j) = G%geoLonT(i,j)*deg_to_rad
  enddo ; enddo
  do j=js-1,je+1 ; do i=is-1,ie+1
    CS%sin_struct(i,j,1) = -sin(2.0*lat_rad(i,j)) * sin(lon_rad(i,j))
    CS%cos_struct(i,j,1) =  sin(2.0*lat_rad(i,j)) * cos(lon_rad(i,j))
    CS%sin_struct(i,j,2) = -cos(lat_rad(i,j))**2 * sin(2.0*lon_rad(i,j))
    CS%cos_struct(i,j,2) =  cos(lat_rad(i,j))**2 * cos(2.0*lon_rad(i,j))
    CS%sin_struct(i,j,3) =  0.0
    CS%cos_struct(i,j,3) = (0.5-1.5*sin(lat_rad(i,j))**2)
  enddo ; enddo

  call get_param(param_file, mdl, "TIDE_M2", use_M2, &
                 "If true, apply tidal momentum forcing at the M2 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_S2", use_S2, &
                 "If true, apply tidal momentum forcing at the S2 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_N2", use_N2, &
                 "If true, apply tidal momentum forcing at the N2 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_K2", use_K2, &
                 "If true, apply tidal momentum forcing at the K2 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_K1", use_K1, &
                 "If true, apply tidal momentum forcing at the K1 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_O1", use_O1, &
                 "If true, apply tidal momentum forcing at the O1 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_P1", use_P1, &
                 "If true, apply tidal momentum forcing at the P1 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_Q1", use_Q1, &
                 "If true, apply tidal momentum forcing at the Q1 "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_MF", use_MF, &
                 "If true, apply tidal momentum forcing at the MF "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)
  call get_param(param_file, mdl, "TIDE_MM", use_MM, &
                 "If true, apply tidal momentum forcing at the MM "//&
                 "frequency. This is only used if TIDES is true.", &
                 default=.false.)

  ! Determine how many tidal components are to be used.
  nc = 0
  if (use_M2) nc=nc+1 ; if (use_S2) nc=nc+1
  if (use_N2) nc=nc+1 ; if (use_K2) nc=nc+1
  if (use_K1) nc=nc+1 ; if (use_O1) nc=nc+1
  if (use_P1) nc=nc+1 ; if (use_Q1) nc=nc+1
  if (use_MF) nc=nc+1 ; if (use_MM) nc=nc+1
  CS%nc = nc

  if (nc == 0) then
    call MOM_error(FATAL, "tidal_forcing_init: "// &
        "TIDES are defined, but no tidal constituents are used.")
    return
  endif

  call get_param(param_file, mdl, "TIDAL_SAL_FROM_FILE", CS%tidal_sal_from_file, &
                 "If true, read the tidal self-attraction and loading "//&
                 "from input files, specified by TIDAL_INPUT_FILE. "//&
                 "This is only used if TIDES is true.", default=.false.)
  call get_param(param_file, mdl, "USE_PREVIOUS_TIDES", CS%use_prev_tides, &
                 "If true, use the SAL from the previous iteration of the "//&
                 "tides to facilitate convergent iteration. "//&
                 "This is only used if TIDES is true.", default=.false.)
  call get_param(param_file, mdl, "TIDE_USE_SAL_SCALAR", CS%use_sal_scalar, &
                 "If true and TIDES is true, use the scalar approximation "//&
                 "when calculating self-attraction and loading.", &
                 default=.not.CS%tidal_sal_from_file)
  ! If it is being used, sal_scalar MUST be specified in param_file.
  if (CS%use_sal_scalar .or. CS%use_prev_tides) &
    call get_param(param_file, mdl, "TIDE_SAL_SCALAR_VALUE", CS%sal_scalar, &
                 "The constant of proportionality between sea surface "//&
                 "height (really it should be bottom pressure) anomalies "//&
                 "and bottom geopotential anomalies. This is only used if "//&
                 "TIDES and TIDE_USE_SAL_SCALAR are true.", units="m m-1", &
                 fail_if_missing=.true.)

  call get_param(param_file, mdl, "TIDAL_SAL_SHT", CS%tidal_sal_sht, &
                 "If true, use inline SAL.", default=.false.)

  if (nc > MAX_CONSTITUENTS) then
    write(mesg,'("Increase MAX_CONSTITUENTS in MOM_tidal_forcing.F90 to at least",I3, &
                &"to accommodate all the registered tidal constituents.")') nc
    call MOM_error(FATAL, "MOM_tidal_forcing"//mesg)
  endif

  do c=1,4*MAX_CONSTITUENTS ; tidal_input_files(c) = "" ; enddo

  if (CS%tidal_sal_from_file .or. CS%use_prev_tides) then
    call get_param(param_file, mdl, "TIDAL_INPUT_FILE", tidal_input_files, &
                   "A list of input files for tidal information.",         &
                   default = "", fail_if_missing=.true.)
  endif

  call get_param(param_file, mdl, "TIDE_REF_DATE", tide_ref_date, &
                 "Year,month,day to use as reference date for tidal forcing. "//&
                 "If not specified, defaults to 0.", &
                 default=0)

  call get_param(param_file, mdl, "TIDE_USE_EQ_PHASE", CS%use_eq_phase, &
                 "Correct phases by calculating equilibrium phase arguments for TIDE_REF_DATE. ", &
                 default=.false., fail_if_missing=.false.)

  if (sum(tide_ref_date) == 0) then  ! tide_ref_date defaults to 0.
    CS%time_ref = set_date(1, 1, 1)
  else
    if (.not. CS%use_eq_phase) then
      ! Using a reference date but not using phase relative to equilibrium.
      ! This makes sense as long as either phases are overridden, or
      ! correctly simulating tidal phases is not desired.
      call MOM_mesg('Tidal phases will *not* be corrected with equilibrium arguments.')
    endif
    CS%time_ref = set_date(tide_ref_date(1), tide_ref_date(2), tide_ref_date(3))
  endif

  ! Initialize reference time for tides and find relevant lunar and solar
  ! longitudes at the reference time.
  if (CS%use_eq_phase) call astro_longitudes_init(CS%time_ref, CS%tidal_longitudes)

  ! Set the parameters for all components that are in use.
  c=0
  if (use_M2) then
    c=c+1 ; CS%const_name(c) = "M2" ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.242334 ! Default amplitude in m.
  endif

  if (use_S2) then
    c=c+1 ; CS%const_name(c) = "S2" ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.112743 ! Default amplitude in m.
  endif

  if (use_N2) then
    c=c+1 ; CS%const_name(c) = "N2" ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.046397 ! Default amplitude in m.
  endif

  if (use_K2) then
    c=c+1 ; CS%const_name(c) = "K2" ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.030684 ! Default amplitude in m.
  endif

  if (use_K1) then
    c=c+1 ; CS%const_name(c) = "K1" ; CS%struct(c) = 1
    CS%love_no(c) = 0.736 ; amp_def(c) = 0.141565 ! Default amplitude in m.
  endif

  if (use_O1) then
    c=c+1 ; CS%const_name(c) = "O1" ; CS%struct(c) = 1
    CS%love_no(c) = 0.695 ; amp_def(c) = 0.100661 ! Default amplitude in m.
  endif

  if (use_P1) then
    c=c+1 ; CS%const_name(c) = "P1" ; CS%struct(c) = 1
    CS%love_no(c) = 0.706 ; amp_def(c) = 0.046848 ! Default amplitude in m.
  endif

  if (use_Q1) then
    c=c+1 ; CS%const_name(c) = "Q1" ; CS%struct(c) = 1
    CS%love_no(c) = 0.695 ; amp_def(c) = 0.019273 ! Default amplitude in m.
  endif

  if (use_MF) then
    c=c+1 ; CS%const_name(c) = "MF" ; CS%struct(c) = 3
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.042041 ! Default amplitude in m.
  endif

  if (use_MM) then
    c=c+1 ; CS%const_name(c) = "MM" ; CS%struct(c) = 3
    CS%love_no(c) = 0.693 ; amp_def(c) = 0.022191 ! Default amplitude in m.
  endif

  ! Set defaults for all included constituents
  ! and things that can be set by functions
  do c=1,nc
    freq_def(c) = tidal_frequency(CS%const_name(c))
    love_def(c) = CS%love_no(c)
    CS%phase0(c) = 0.0
    if (CS%use_eq_phase) then
      phase0_def(c) = eq_phase(CS%const_name(c), CS%tidal_longitudes)
    else
      phase0_def(c) = 0.0
    endif
  enddo

  !  Parse the input file to potentially override the default values for the
  ! frequency, amplitude and initial phase of each constituent, and log the
  ! values that are actually used.
  do c=1,nc
    call get_param(param_file, mdl, "TIDE_"//trim(CS%const_name(c))//"_FREQ", CS%freq(c), &
                   "Frequency of the "//trim(CS%const_name(c))//" tidal constituent. "//&
                   "This is only used if TIDES and TIDE_"//trim(CS%const_name(c))// &
                   " are true, or if OBC_TIDE_N_CONSTITUENTS > 0 and "//trim(CS%const_name(c))// &
                   " is in OBC_TIDE_CONSTITUENTS.", units="s-1", default=freq_def(c), scale=US%T_to_s)
    call get_param(param_file, mdl, "TIDE_"//trim(CS%const_name(c))//"_AMP", CS%amp(c), &
                   "Amplitude of the "//trim(CS%const_name(c))//" tidal constituent. "//&
                   "This is only used if TIDES and TIDE_"//trim(CS%const_name(c))// &
                   " are true.", units="m", default=amp_def(c), scale=US%m_to_Z)
    call get_param(param_file, mdl, "TIDE_"//trim(CS%const_name(c))//"_PHASE_T0", CS%phase0(c), &
                   "Phase of the "//trim(CS%const_name(c))//" tidal constituent at time 0. "//&
                   "This is only used if TIDES and TIDE_"//trim(CS%const_name(c))// &
                   " are true.", units="radians", default=phase0_def(c))
  enddo

  if (CS%tidal_sal_from_file) then
    allocate(CS%cosphasesal(isd:ied,jsd:jed,nc))
    allocate(CS%sinphasesal(isd:ied,jsd:jed,nc))
    allocate(CS%ampsal(isd:ied,jsd:jed,nc))
    do c=1,nc
      ! Read variables with names like PHASE_SAL_M2 and AMP_SAL_M2.
      call find_in_files(tidal_input_files, "PHASE_SAL_"//trim(CS%const_name(c)), phase, G)
      call find_in_files(tidal_input_files, "AMP_SAL_"//trim(CS%const_name(c)), CS%ampsal(:,:,c), &
                         G, scale=US%m_to_Z)
      call pass_var(phase,           G%domain,complete=.false.)
      call pass_var(CS%ampsal(:,:,c),G%domain,complete=.true.)
      do j=js-1,je+1 ; do i=is-1,ie+1
        CS%cosphasesal(i,j,c) = cos(phase(i,j)*deg_to_rad)
        CS%sinphasesal(i,j,c) = sin(phase(i,j)*deg_to_rad)
      enddo ; enddo
    enddo
  endif

  if (CS%USE_PREV_TIDES) then
    allocate(CS%cosphase_prev(isd:ied,jsd:jed,nc))
    allocate(CS%sinphase_prev(isd:ied,jsd:jed,nc))
    allocate(CS%amp_prev(isd:ied,jsd:jed,nc))
    do c=1,nc
      ! Read variables with names like PHASE_PREV_M2 and AMP_PREV_M2.
      call find_in_files(tidal_input_files, "PHASE_PREV_"//trim(CS%const_name(c)), phase, G)
      call find_in_files(tidal_input_files, "AMP_PREV_"//trim(CS%const_name(c)), CS%amp_prev(:,:,c), &
                         G, scale=US%m_to_Z)
      call pass_var(phase,             G%domain,complete=.false.)
      call pass_var(CS%amp_prev(:,:,c),G%domain,complete=.true.)
      do j=js-1,je+1 ; do i=is-1,ie+1
        CS%cosphase_prev(i,j,c) = cos(phase(i,j)*deg_to_rad)
        CS%sinphase_prev(i,j,c) = sin(phase(i,j)*deg_to_rad)
      enddo ; enddo
    enddo
  endif

  if (CS%tidal_sal_sht) then
    call spherical_harmonics_init(G, param_file, CS%sht)
    id_clock_SAL = cpu_clock_id('(Ocean SAL)', grain=CLOCK_MODULE)
  endif

  id_clock_tides = cpu_clock_id('(Ocean tides)', grain=CLOCK_MODULE)

end subroutine tidal_forcing_init

subroutine getloadLoveNums(nlm, LoveScaling) !{{{

  integer, intent(in) :: nlm
  real, dimension(:), intent(out) :: LoveScaling

  real, dimension(:), allocatable :: H, L, K
  real, dimension(:,:), allocatable :: LoveDat
  real :: H1, L1, K1
  integer :: i, j, n, m
  real, parameter :: rhoE=5517.0 ! Average density of Earth (kg/m^3)
  real, parameter :: rhoW=1035.0 ! Density of water (kg/m^3)
  integer, parameter :: lmax=1440

  allocate(LoveDat(4,lmax+1))

  LoveDat(1:4,1) = (/ 0.0, 0.0000000000, 0.0000000000, -1.0000000000 /) !{{{
  LoveDat(1:4,2) = (/ 1.0, -1.2858777580, -8.9608179370e-1, -1.0000000000 /)
  LoveDat(1:4,3) = (/ 2.0, -9.9079949000e-1, 2.3286695000e-2, -3.0516104000e-1 /)
  LoveDat(1:4,4) = (/ 3.0, -1.0499631000, 6.9892136000e-2, -1.9585733000e-1 /)
  LoveDat(1:4,5) = (/ 4.0, -1.0526477000, 5.8670467000e-2, -1.3352284000e-1 /)
  LoveDat(1:4,6) = (/ 5.0, -1.0855918000, 4.6165153000e-2, -1.0456531000e-1 /)
  LoveDat(1:4,7) = (/ 6.0, -1.1431163000, 3.8586926000e-2, -9.0184841000e-2 /)
  LoveDat(1:4,8) = (/ 7.0, -1.2116273000, 3.4198827000e-2, -8.1906787000e-2 /)
  LoveDat(1:4,9) = (/ 8.0, -1.2831157000, 3.1474998000e-2, -7.6379141000e-2 /)
  LoveDat(1:4,10) = (/ 9.0, -1.3538554000, 2.9624407000e-2, -7.2250183000e-2 /)
  LoveDat(1:4,11) = (/ 10.0, -1.4223516000, 2.8273961000e-2, -6.8934145000e-2 /)
  LoveDat(1:4,12) = (/ 11.0, -1.4881117000, 2.7242278000e-2, -6.6147992000e-2 /)
  LoveDat(1:4,13) = (/ 12.0, -1.5510428000, 2.6431124000e-2, -6.3736253000e-2 /)
  LoveDat(1:4,14) = (/ 13.0, -1.6111895000, 2.5779507000e-2, -6.1602870000e-2 /)
  LoveDat(1:4,15) = (/ 14.0, -1.6686329000, 2.5245139000e-2, -5.9683159000e-2 /)
  LoveDat(1:4,16) = (/ 15.0, -1.7234569000, 2.4796803000e-2, -5.7931180000e-2 /)
  LoveDat(1:4,17) = (/ 16.0, -1.7757418000, 2.4410861000e-2, -5.6313294000e-2 /)
  LoveDat(1:4,18) = (/ 17.0, -1.8255646000, 2.4069336000e-2, -5.4804452000e-2 /)
  LoveDat(1:4,19) = (/ 18.0, -1.8730019000, 2.3758645000e-2, -5.3385807000e-2 /)
  LoveDat(1:4,20) = (/ 19.0, -1.9181321000, 2.3468646000e-2, -5.2043088000e-2 /)
  LoveDat(1:4,21) = (/ 20.0, -1.9610366000, 2.3191893000e-2, -5.0765423000e-2 /)
  LoveDat(1:4,22) = (/ 21.0, -2.0018000000, 2.2923032000e-2, -4.9544487000e-2 /)
  LoveDat(1:4,23) = (/ 22.0, -2.0405101000, 2.2658321000e-2, -4.8373866000e-2 /)
  LoveDat(1:4,24) = (/ 23.0, -2.0772571000, 2.2395242000e-2, -4.7248575000e-2 /)
  LoveDat(1:4,25) = (/ 24.0, -2.1121328000, 2.2132200000e-2, -4.6164708000e-2 /)
  LoveDat(1:4,26) = (/ 25.0, -2.1452296000, 2.1868280000e-2, -4.5119160000e-2 /)
  LoveDat(1:4,27) = (/ 26.0, -2.1766398000, 2.1603063000e-2, -4.4109431000e-2 /)
  LoveDat(1:4,28) = (/ 27.0, -2.2064546000, 2.1336479000e-2, -4.3133464000e-2 /)
  LoveDat(1:4,29) = (/ 28.0, -2.2347634000, 2.1068700000e-2, -4.2189540000e-2 /)
  LoveDat(1:4,30) = (/ 29.0, -2.2616531000, 2.0800053000e-2, -4.1276184000e-2 /)
  LoveDat(1:4,31) = (/ 30.0, -2.2872080000, 2.0530962000e-2, -4.0392105000e-2 /)
  LoveDat(1:4,32) = (/ 31.0, -2.3115088000, 2.0261897000e-2, -3.9536148000e-2 /)
  LoveDat(1:4,33) = (/ 32.0, -2.3346328000, 1.9993346000e-2, -3.8707260000e-2 /)
  LoveDat(1:4,34) = (/ 33.0, -2.3566536000, 1.9725790000e-2, -3.7904463000e-2 /)
  LoveDat(1:4,35) = (/ 34.0, -2.3776409000, 1.9459686000e-2, -3.7126837000e-2 /)
  LoveDat(1:4,36) = (/ 35.0, -2.3976605000, 1.9195459000e-2, -3.6373510000e-2 /)
  LoveDat(1:4,37) = (/ 36.0, -2.4167746000, 1.8933494000e-2, -3.5643644000e-2 /)
  LoveDat(1:4,38) = (/ 37.0, -2.4350414000, 1.8674136000e-2, -3.4936432000e-2 /)
  LoveDat(1:4,39) = (/ 38.0, -2.4525156000, 1.8417687000e-2, -3.4251094000e-2 /)
  LoveDat(1:4,40) = (/ 39.0, -2.4692484000, 1.8164407000e-2, -3.3586873000e-2 /)
  LoveDat(1:4,41) = (/ 40.0, -2.4852876000, 1.7914518000e-2, -3.2943035000e-2 /)
  LoveDat(1:4,42) = (/ 41.0, -2.5006779000, 1.7668203000e-2, -3.2318866000e-2 /)
  LoveDat(1:4,43) = (/ 42.0, -2.5154609000, 1.7425613000e-2, -3.1713675000e-2 /)
  LoveDat(1:4,44) = (/ 43.0, -2.5296755000, 1.7186866000e-2, -3.1126789000e-2 /)
  LoveDat(1:4,45) = (/ 44.0, -2.5433577000, 1.6952053000e-2, -3.0557557000e-2 /)
  LoveDat(1:4,46) = (/ 45.0, -2.5565412000, 1.6721240000e-2, -3.0005347000e-2 /)
  LoveDat(1:4,47) = (/ 46.0, -2.5692574000, 1.6494470000e-2, -2.9469547000e-2 /)
  LoveDat(1:4,48) = (/ 47.0, -2.5815353000, 1.6271769000e-2, -2.8949568000e-2 /)
  LoveDat(1:4,49) = (/ 48.0, -2.5934022000, 1.6053144000e-2, -2.8444838000e-2 /)
  LoveDat(1:4,50) = (/ 49.0, -2.6048833000, 1.5838586000e-2, -2.7954806000e-2 /)
  LoveDat(1:4,51) = (/ 50.0, -2.6160021000, 1.5628077000e-2, -2.7478940000e-2 /)
  LoveDat(1:4,52) = (/ 51.0, -2.6267805000, 1.5421585000e-2, -2.7016729000e-2 /)
  LoveDat(1:4,53) = (/ 52.0, -2.6372389000, 1.5219071000e-2, -2.6567679000e-2 /)
  LoveDat(1:4,54) = (/ 53.0, -2.6473964000, 1.5020486000e-2, -2.6131317000e-2 /)
  LoveDat(1:4,55) = (/ 54.0, -2.6572706000, 1.4825779000e-2, -2.5707185000e-2 /)
  LoveDat(1:4,56) = (/ 55.0, -2.6668781000, 1.4634888000e-2, -2.5294846000e-2 /)
  LoveDat(1:4,57) = (/ 56.0, -2.6762345000, 1.4447752000e-2, -2.4893877000e-2 /)
  LoveDat(1:4,58) = (/ 57.0, -2.6853540000, 1.4264303000e-2, -2.4503874000e-2 /)
  LoveDat(1:4,59) = (/ 58.0, -2.6942503000, 1.4084474000e-2, -2.4124449000e-2 /)
  LoveDat(1:4,60) = (/ 59.0, -2.7029358000, 1.3908192000e-2, -2.3755228000e-2 /)
  LoveDat(1:4,61) = (/ 60.0, -2.7114225000, 1.3735386000e-2, -2.3395852000e-2 /)
  LoveDat(1:4,62) = (/ 61.0, -2.7197214000, 1.3565983000e-2, -2.3045980000e-2 /)
  LoveDat(1:4,63) = (/ 62.0, -2.7278428000, 1.3399909000e-2, -2.2705280000e-2 /)
  LoveDat(1:4,64) = (/ 63.0, -2.7357965000, 1.3237092000e-2, -2.2373437000e-2 /)
  LoveDat(1:4,65) = (/ 64.0, -2.7435916000, 1.3077458000e-2, -2.2050147000e-2 /)
  LoveDat(1:4,66) = (/ 65.0, -2.7512366000, 1.2920935000e-2, -2.1735119000e-2 /)
  LoveDat(1:4,67) = (/ 66.0, -2.7587397000, 1.2767451000e-2, -2.1428073000e-2 /)
  LoveDat(1:4,68) = (/ 67.0, -2.7661083000, 1.2616936000e-2, -2.1128742000e-2 /)
  LoveDat(1:4,69) = (/ 68.0, -2.7733496000, 1.2469319000e-2, -2.0836869000e-2 /)
  LoveDat(1:4,70) = (/ 69.0, -2.7804703000, 1.2324532000e-2, -2.0552206000e-2 /)
  LoveDat(1:4,71) = (/ 70.0, -2.7874767000, 1.2182508000e-2, -2.0274516000e-2 /)
  LoveDat(1:4,72) = (/ 71.0, -2.7943748000, 1.2043181000e-2, -2.0003572000e-2 /)
  LoveDat(1:4,73) = (/ 72.0, -2.8011702000, 1.1906487000e-2, -1.9739156000e-2 /)
  LoveDat(1:4,74) = (/ 73.0, -2.8078682000, 1.1772362000e-2, -1.9481058000e-2 /)
  LoveDat(1:4,75) = (/ 74.0, -2.8144738000, 1.1640746000e-2, -1.9229076000e-2 /)
  LoveDat(1:4,76) = (/ 75.0, -2.8209918000, 1.1511578000e-2, -1.8983017000e-2 /)
  LoveDat(1:4,77) = (/ 76.0, -2.8274266000, 1.1384799000e-2, -1.8742695000e-2 /)
  LoveDat(1:4,78) = (/ 77.0, -2.8337824000, 1.1260352000e-2, -1.8507931000e-2 /)
  LoveDat(1:4,79) = (/ 78.0, -2.8400633000, 1.1138183000e-2, -1.8278553000e-2 /)
  LoveDat(1:4,80) = (/ 79.0, -2.8462730000, 1.1018236000e-2, -1.8054395000e-2 /)
  LoveDat(1:4,81) = (/ 80.0, -2.8524152000, 1.0900460000e-2, -1.7835300000e-2 /)
  LoveDat(1:4,82) = (/ 81.0, -2.8584932000, 1.0784802000e-2, -1.7621113000e-2 /)
  LoveDat(1:4,83) = (/ 82.0, -2.8645103000, 1.0671213000e-2, -1.7411688000e-2 /)
  LoveDat(1:4,84) = (/ 83.0, -2.8704696000, 1.0559645000e-2, -1.7206882000e-2 /)
  LoveDat(1:4,85) = (/ 84.0, -2.8763739000, 1.0450051000e-2, -1.7006560000e-2 /)
  LoveDat(1:4,86) = (/ 85.0, -2.8822260000, 1.0342384000e-2, -1.6810590000e-2 /)
  LoveDat(1:4,87) = (/ 86.0, -2.8880285000, 1.0236599000e-2, -1.6618845000e-2 /)
  LoveDat(1:4,88) = (/ 87.0, -2.8937839000, 1.0132655000e-2, -1.6431203000e-2 /)
  LoveDat(1:4,89) = (/ 88.0, -2.8994945000, 1.0030508000e-2, -1.6247547000e-2 /)
  LoveDat(1:4,90) = (/ 89.0, -2.9051627000, 9.9301169000e-3, -1.6067762000e-2 /)
  LoveDat(1:4,91) = (/ 90.0, -2.9107905000, 9.8314429000e-3, -1.5891741000e-2 /)
  LoveDat(1:4,92) = (/ 91.0, -2.9163799000, 9.7344467000e-3, -1.5719376000e-2 /)
  LoveDat(1:4,93) = (/ 92.0, -2.9219330000, 9.6390907000e-3, -1.5550567000e-2 /)
  LoveDat(1:4,94) = (/ 93.0, -2.9274514000, 9.5453383000e-3, -1.5385215000e-2 /)
  LoveDat(1:4,95) = (/ 94.0, -2.9329370000, 9.4531538000e-3, -1.5223225000e-2 /)
  LoveDat(1:4,96) = (/ 95.0, -2.9383913000, 9.3625026000e-3, -1.5064506000e-2 /)
  LoveDat(1:4,97) = (/ 96.0, -2.9438161000, 9.2733509000e-3, -1.4908968000e-2 /)
  LoveDat(1:4,98) = (/ 97.0, -2.9492127000, 9.1856660000e-3, -1.4756526000e-2 /)
  LoveDat(1:4,99) = (/ 98.0, -2.9545826000, 9.0994159000e-3, -1.4607099000e-2 /)
  LoveDat(1:4,100) = (/ 99.0, -2.9599272000, 9.0145695000e-3, -1.4460604000e-2 /)
  LoveDat(1:4,101) = (/ 100.0, -2.9652476000, 8.9310967000e-3, -1.4316967000e-2 /)
  LoveDat(1:4,102) = (/ 101.0, -2.9705453000, 8.8489681000e-3, -1.4176111000e-2 /)
  LoveDat(1:4,103) = (/ 102.0, -2.9758213000, 8.7681548000e-3, -1.4037965000e-2 /)
  LoveDat(1:4,104) = (/ 103.0, -2.9810767000, 8.6886292000e-3, -1.3902458000e-2 /)
  LoveDat(1:4,105) = (/ 104.0, -2.9863125000, 8.6103640000e-3, -1.3769523000e-2 /)
  LoveDat(1:4,106) = (/ 105.0, -2.9915299000, 8.5333328000e-3, -1.3639094000e-2 /)
  LoveDat(1:4,107) = (/ 106.0, -2.9967298000, 8.4575097000e-3, -1.3511108000e-2 /)
  LoveDat(1:4,108) = (/ 107.0, -3.0019129000, 8.3828699000e-3, -1.3385503000e-2 /)
  LoveDat(1:4,109) = (/ 108.0, -3.0070803000, 8.3093886000e-3, -1.3262220000e-2 /)
  LoveDat(1:4,110) = (/ 109.0, -3.0122328000, 8.2370423000e-3, -1.3141201000e-2 /)
  LoveDat(1:4,111) = (/ 110.0, -3.0173710000, 8.1658076000e-3, -1.3022390000e-2 /)
  LoveDat(1:4,112) = (/ 111.0, -3.0224958000, 8.0956619000e-3, -1.2905734000e-2 /)
  LoveDat(1:4,113) = (/ 112.0, -3.0276079000, 8.0265832000e-3, -1.2791179000e-2 /)
  LoveDat(1:4,114) = (/ 113.0, -3.0327080000, 7.9585500000e-3, -1.2678675000e-2 /)
  LoveDat(1:4,115) = (/ 114.0, -3.0377966000, 7.8915413000e-3, -1.2568172000e-2 /)
  LoveDat(1:4,116) = (/ 115.0, -3.0428744000, 7.8255367000e-3, -1.2459622000e-2 /)
  LoveDat(1:4,117) = (/ 116.0, -3.0479420000, 7.7605163000e-3, -1.2352979000e-2 /)
  LoveDat(1:4,118) = (/ 117.0, -3.0529999000, 7.6964606000e-3, -1.2248198000e-2 /)
  LoveDat(1:4,119) = (/ 118.0, -3.0580486000, 7.6333507000e-3, -1.2145235000e-2 /)
  LoveDat(1:4,120) = (/ 119.0, -3.0630887000, 7.5711680000e-3, -1.2044048000e-2 /)
  LoveDat(1:4,121) = (/ 120.0, -3.0681205000, 7.5098946000e-3, -1.1944594000e-2 /)
  LoveDat(1:4,122) = (/ 121.0, -3.0731446000, 7.4495128000e-3, -1.1846835000e-2 /)
  LoveDat(1:4,123) = (/ 122.0, -3.0781614000, 7.3900054000e-3, -1.1750732000e-2 /)
  LoveDat(1:4,124) = (/ 123.0, -3.0831713000, 7.3313557000e-3, -1.1656245000e-2 /)
  LoveDat(1:4,125) = (/ 124.0, -3.0881747000, 7.2735474000e-3, -1.1563340000e-2 /)
  LoveDat(1:4,126) = (/ 125.0, -3.0931718000, 7.2165644000e-3, -1.1471980000e-2 /)
  LoveDat(1:4,127) = (/ 126.0, -3.0981632000, 7.1603911000e-3, -1.1382130000e-2 /)
  LoveDat(1:4,128) = (/ 127.0, -3.1031490000, 7.1050124000e-3, -1.1293757000e-2 /)
  LoveDat(1:4,129) = (/ 128.0, -3.1081296000, 7.0504134000e-3, -1.1206828000e-2 /)
  LoveDat(1:4,130) = (/ 129.0, -3.1131054000, 6.9965795000e-3, -1.1121311000e-2 /)
  LoveDat(1:4,131) = (/ 130.0, -3.1180765000, 6.9434967000e-3, -1.1037175000e-2 /)
  LoveDat(1:4,132) = (/ 131.0, -3.1230433000, 6.8911509000e-3, -1.0954391000e-2 /)
  LoveDat(1:4,133) = (/ 132.0, -3.1280059000, 6.8395288000e-3, -1.0872928000e-2 /)
  LoveDat(1:4,134) = (/ 133.0, -3.1329647000, 6.7886171000e-3, -1.0792758000e-2 /)
  LoveDat(1:4,135) = (/ 134.0, -3.1379199000, 6.7384029000e-3, -1.0713853000e-2 /)
  LoveDat(1:4,136) = (/ 135.0, -3.1428716000, 6.6888735000e-3, -1.0636187000e-2 /)
  LoveDat(1:4,137) = (/ 136.0, -3.1478201000, 6.6400168000e-3, -1.0559733000e-2 /)
  LoveDat(1:4,138) = (/ 137.0, -3.1527656000, 6.5918206000e-3, -1.0484466000e-2 /)
  LoveDat(1:4,139) = (/ 138.0, -3.1577082000, 6.5442732000e-3, -1.0410360000e-2 /)
  LoveDat(1:4,140) = (/ 139.0, -3.1626481000, 6.4973631000e-3, -1.0337392000e-2 /)
  LoveDat(1:4,141) = (/ 140.0, -3.1675855000, 6.4510790000e-3, -1.0265537000e-2 /)
  LoveDat(1:4,142) = (/ 141.0, -3.1725205000, 6.4054099000e-3, -1.0194773000e-2 /)
  LoveDat(1:4,143) = (/ 142.0, -3.1774533000, 6.3603452000e-3, -1.0125078000e-2 /)
  LoveDat(1:4,144) = (/ 143.0, -3.1823840000, 6.3158742000e-3, -1.0056429000e-2 /)
  LoveDat(1:4,145) = (/ 144.0, -3.1873127000, 6.2719868000e-3, -9.9888045000e-3 /)
  LoveDat(1:4,146) = (/ 145.0, -3.1922396000, 6.2286729000e-3, -9.9221850000e-3 /)
  LoveDat(1:4,147) = (/ 146.0, -3.1971648000, 6.1859227000e-3, -9.8565496000e-3 /)
  LoveDat(1:4,148) = (/ 147.0, -3.2020883000, 6.1437265000e-3, -9.7918788000e-3 /)
  LoveDat(1:4,149) = (/ 148.0, -3.2070102000, 6.1020749000e-3, -9.7281532000e-3 /)
  LoveDat(1:4,150) = (/ 149.0, -3.2119308000, 6.0609589000e-3, -9.6653542000e-3 /)
  LoveDat(1:4,151) = (/ 150.0, -3.2168500000, 6.0203693000e-3, -9.6034635000e-3 /)
  LoveDat(1:4,152) = (/ 151.0, -3.2217679000, 5.9802974000e-3, -9.5424633000e-3 /)
  LoveDat(1:4,153) = (/ 152.0, -3.2266847000, 5.9407346000e-3, -9.4823362000e-3 /)
  LoveDat(1:4,154) = (/ 153.0, -3.2316003000, 5.9016724000e-3, -9.4230652000e-3 /)
  LoveDat(1:4,155) = (/ 154.0, -3.2365149000, 5.8631026000e-3, -9.3646338000e-3 /)
  LoveDat(1:4,156) = (/ 155.0, -3.2414284000, 5.8250172000e-3, -9.3070259000e-3 /)
  LoveDat(1:4,157) = (/ 156.0, -3.2463411000, 5.7874081000e-3, -9.2502257000e-3 /)
  LoveDat(1:4,158) = (/ 157.0, -3.2512529000, 5.7502678000e-3, -9.1942178000e-3 /)
  LoveDat(1:4,159) = (/ 158.0, -3.2561639000, 5.7135886000e-3, -9.1389873000e-3 /)
  LoveDat(1:4,160) = (/ 159.0, -3.2610741000, 5.6773630000e-3, -9.0845194000e-3 /)
  LoveDat(1:4,161) = (/ 160.0, -3.2659835000, 5.6415839000e-3, -9.0308000000e-3 /)
  LoveDat(1:4,162) = (/ 161.0, -3.2708923000, 5.6062442000e-3, -8.9778149000e-3 /)
  LoveDat(1:4,163) = (/ 162.0, -3.2758004000, 5.5713368000e-3, -8.9255506000e-3 /)
  LoveDat(1:4,164) = (/ 163.0, -3.2807079000, 5.5368550000e-3, -8.8739938000e-3 /)
  LoveDat(1:4,165) = (/ 164.0, -3.2856148000, 5.5027920000e-3, -8.8231314000e-3 /)
  LoveDat(1:4,166) = (/ 165.0, -3.2905211000, 5.4691413000e-3, -8.7729507000e-3 /)
  LoveDat(1:4,167) = (/ 166.0, -3.2954269000, 5.4358966000e-3, -8.7234394000e-3 /)
  LoveDat(1:4,168) = (/ 167.0, -3.3003322000, 5.4030515000e-3, -8.6745852000e-3 /)
  LoveDat(1:4,169) = (/ 168.0, -3.3052370000, 5.3705998000e-3, -8.6263763000e-3 /)
  LoveDat(1:4,170) = (/ 169.0, -3.3101414000, 5.3385356000e-3, -8.5788012000e-3 /)
  LoveDat(1:4,171) = (/ 170.0, -3.3150452000, 5.3068529000e-3, -8.5318484000e-3 /)
  LoveDat(1:4,172) = (/ 171.0, -3.3199486000, 5.2755459000e-3, -8.4855070000e-3 /)
  LoveDat(1:4,173) = (/ 172.0, -3.3248516000, 5.2446089000e-3, -8.4397661000e-3 /)
  LoveDat(1:4,174) = (/ 173.0, -3.3297541000, 5.2140364000e-3, -8.3946150000e-3 /)
  LoveDat(1:4,175) = (/ 174.0, -3.3346563000, 5.1838229000e-3, -8.3500435000e-3 /)
  LoveDat(1:4,176) = (/ 175.0, -3.3395580000, 5.1539630000e-3, -8.3060415000e-3 /)
  LoveDat(1:4,177) = (/ 176.0, -3.3444593000, 5.1244515000e-3, -8.2625990000e-3 /)
  LoveDat(1:4,178) = (/ 177.0, -3.3493602000, 5.0952833000e-3, -8.2197063000e-3 /)
  LoveDat(1:4,179) = (/ 178.0, -3.3542607000, 5.0664532000e-3, -8.1773539000e-3 /)
  LoveDat(1:4,180) = (/ 179.0, -3.3591609000, 5.0379563000e-3, -8.1355327000e-3 /)
  LoveDat(1:4,181) = (/ 180.0, -3.3640606000, 5.0097879000e-3, -8.0942335000e-3 /)
  LoveDat(1:4,182) = (/ 181.0, -3.3689599000, 4.9819430000e-3, -8.0534474000e-3 /)
  LoveDat(1:4,183) = (/ 182.0, -3.3738588000, 4.9544170000e-3, -8.0131658000e-3 /)
  LoveDat(1:4,184) = (/ 183.0, -3.3787572000, 4.9272053000e-3, -7.9733801000e-3 /)
  LoveDat(1:4,185) = (/ 184.0, -3.3836553000, 4.9003034000e-3, -7.9340821000e-3 /)
  LoveDat(1:4,186) = (/ 185.0, -3.3885529000, 4.8737069000e-3, -7.8952635000e-3 /)
  LoveDat(1:4,187) = (/ 186.0, -3.3934501000, 4.8474114000e-3, -7.8569164000e-3 /)
  LoveDat(1:4,188) = (/ 187.0, -3.3983469000, 4.8214127000e-3, -7.8190330000e-3 /)
  LoveDat(1:4,189) = (/ 188.0, -3.4032432000, 4.7957066000e-3, -7.7816057000e-3 /)
  LoveDat(1:4,190) = (/ 189.0, -3.4081390000, 4.7702889000e-3, -7.7446269000e-3 /)
  LoveDat(1:4,191) = (/ 190.0, -3.4130344000, 4.7451557000e-3, -7.7080893000e-3 /)
  LoveDat(1:4,192) = (/ 191.0, -3.4179292000, 4.7203030000e-3, -7.6719857000e-3 /)
  LoveDat(1:4,193) = (/ 192.0, -3.4228236000, 4.6957268000e-3, -7.6363091000e-3 /)
  LoveDat(1:4,194) = (/ 193.0, -3.4277174000, 4.6714235000e-3, -7.6010526000e-3 /)
  LoveDat(1:4,195) = (/ 194.0, -3.4326107000, 4.6473891000e-3, -7.5662095000e-3 /)
  LoveDat(1:4,196) = (/ 195.0, -3.4375035000, 4.6236200000e-3, -7.5317730000e-3 /)
  LoveDat(1:4,197) = (/ 196.0, -3.4423957000, 4.6001126000e-3, -7.4977367000e-3 /)
  LoveDat(1:4,198) = (/ 197.0, -3.4472873000, 4.5768634000e-3, -7.4640943000e-3 /)
  LoveDat(1:4,199) = (/ 198.0, -3.4521783000, 4.5538688000e-3, -7.4308395000e-3 /)
  LoveDat(1:4,200) = (/ 199.0, -3.4570687000, 4.5311254000e-3, -7.3979662000e-3 /)
  LoveDat(1:4,201) = (/ 200.0, -3.4619585000, 4.5086298000e-3, -7.3654685000e-3 /)
  LoveDat(1:4,202) = (/ 201.0, -3.4668476000, 4.4863788000e-3, -7.3333403000e-3 /)
  LoveDat(1:4,203) = (/ 202.0, -3.4717360000, 4.4643689000e-3, -7.3015761000e-3 /)
  LoveDat(1:4,204) = (/ 203.0, -3.4766237000, 4.4425971000e-3, -7.2701701000e-3 /)
  LoveDat(1:4,205) = (/ 204.0, -3.4815107000, 4.4210601000e-3, -7.2391168000e-3 /)
  LoveDat(1:4,206) = (/ 205.0, -3.4863970000, 4.3997550000e-3, -7.2084108000e-3 /)
  LoveDat(1:4,207) = (/ 206.0, -3.4912825000, 4.3786785000e-3, -7.1780467000e-3 /)
  LoveDat(1:4,208) = (/ 207.0, -3.4961672000, 4.3578278000e-3, -7.1480193000e-3 /)
  LoveDat(1:4,209) = (/ 208.0, -3.5010512000, 4.3371999000e-3, -7.1183236000e-3 /)
  LoveDat(1:4,210) = (/ 209.0, -3.5059343000, 4.3167918000e-3, -7.0889544000e-3 /)
  LoveDat(1:4,211) = (/ 210.0, -3.5108165000, 4.2966008000e-3, -7.0599068000e-3 /)
  LoveDat(1:4,212) = (/ 211.0, -3.5156979000, 4.2766239000e-3, -7.0311760000e-3 /)
  LoveDat(1:4,213) = (/ 212.0, -3.5205784000, 4.2568586000e-3, -7.0027573000e-3 /)
  LoveDat(1:4,214) = (/ 213.0, -3.5254580000, 4.2373019000e-3, -6.9746460000e-3 /)
  LoveDat(1:4,215) = (/ 214.0, -3.5303366000, 4.2179514000e-3, -6.9468375000e-3 /)
  LoveDat(1:4,216) = (/ 215.0, -3.5352143000, 4.1988043000e-3, -6.9193272000e-3 /)
  LoveDat(1:4,217) = (/ 216.0, -3.5400909000, 4.1798580000e-3, -6.8921109000e-3 /)
  LoveDat(1:4,218) = (/ 217.0, -3.5449666000, 4.1611101000e-3, -6.8651842000e-3 /)
  LoveDat(1:4,219) = (/ 218.0, -3.5498412000, 4.1425580000e-3, -6.8385428000e-3 /)
  LoveDat(1:4,220) = (/ 219.0, -3.5547147000, 4.1241992000e-3, -6.8121826000e-3 /)
  LoveDat(1:4,221) = (/ 220.0, -3.5595871000, 4.1060313000e-3, -6.7860995000e-3 /)
  LoveDat(1:4,222) = (/ 221.0, -3.5644584000, 4.0880520000e-3, -6.7602894000e-3 /)
  LoveDat(1:4,223) = (/ 222.0, -3.5693286000, 4.0702588000e-3, -6.7347484000e-3 /)
  LoveDat(1:4,224) = (/ 223.0, -3.5741976000, 4.0526495000e-3, -6.7094726000e-3 /)
  LoveDat(1:4,225) = (/ 224.0, -3.5790654000, 4.0352217000e-3, -6.6844583000e-3 /)
  LoveDat(1:4,226) = (/ 225.0, -3.5839320000, 4.0179733000e-3, -6.6597016000e-3 /)
  LoveDat(1:4,227) = (/ 226.0, -3.5887973000, 4.0009020000e-3, -6.6351989000e-3 /)
  LoveDat(1:4,228) = (/ 227.0, -3.5936613000, 3.9840057000e-3, -6.6109466000e-3 /)
  LoveDat(1:4,229) = (/ 228.0, -3.5985240000, 3.9672821000e-3, -6.5869411000e-3 /)
  LoveDat(1:4,230) = (/ 229.0, -3.6033854000, 3.9507293000e-3, -6.5631791000e-3 /)
  LoveDat(1:4,231) = (/ 230.0, -3.6082455000, 3.9343450000e-3, -6.5396569000e-3 /)
  LoveDat(1:4,232) = (/ 231.0, -3.6131041000, 3.9181273000e-3, -6.5163713000e-3 /)
  LoveDat(1:4,233) = (/ 232.0, -3.6179613000, 3.9020742000e-3, -6.4933190000e-3 /)
  LoveDat(1:4,234) = (/ 233.0, -3.6228171000, 3.8861836000e-3, -6.4704966000e-3 /)
  LoveDat(1:4,235) = (/ 234.0, -3.6276714000, 3.8704536000e-3, -6.4479012000e-3 /)
  LoveDat(1:4,236) = (/ 235.0, -3.6325242000, 3.8548822000e-3, -6.4255293000e-3 /)
  LoveDat(1:4,237) = (/ 236.0, -3.6373754000, 3.8394677000e-3, -6.4033781000e-3 /)
  LoveDat(1:4,238) = (/ 237.0, -3.6422252000, 3.8242080000e-3, -6.3814445000e-3 /)
  LoveDat(1:4,239) = (/ 238.0, -3.6470733000, 3.8091013000e-3, -6.3597254000e-3 /)
  LoveDat(1:4,240) = (/ 239.0, -3.6519198000, 3.7941458000e-3, -6.3382179000e-3 /)
  LoveDat(1:4,241) = (/ 240.0, -3.6567647000, 3.7793398000e-3, -6.3169193000e-3 /)
  LoveDat(1:4,242) = (/ 241.0, -3.6616079000, 3.7646814000e-3, -6.2958265000e-3 /)
  LoveDat(1:4,243) = (/ 242.0, -3.6664494000, 3.7501690000e-3, -6.2749370000e-3 /)
  LoveDat(1:4,244) = (/ 243.0, -3.6712891000, 3.7358007000e-3, -6.2542478000e-3 /)
  LoveDat(1:4,245) = (/ 244.0, -3.6761271000, 3.7215749000e-3, -6.2337563000e-3 /)
  LoveDat(1:4,246) = (/ 245.0, -3.6809634000, 3.7074899000e-3, -6.2134599000e-3 /)
  LoveDat(1:4,247) = (/ 246.0, -3.6857978000, 3.6935441000e-3, -6.1933559000e-3 /)
  LoveDat(1:4,248) = (/ 247.0, -3.6906303000, 3.6797359000e-3, -6.1734419000e-3 /)
  LoveDat(1:4,249) = (/ 248.0, -3.6954610000, 3.6660636000e-3, -6.1537152000e-3 /)
  LoveDat(1:4,250) = (/ 249.0, -3.7002898000, 3.6525257000e-3, -6.1341734000e-3 /)
  LoveDat(1:4,251) = (/ 250.0, -3.7051167000, 3.6391206000e-3, -6.1148140000e-3 /)
  LoveDat(1:4,252) = (/ 251.0, -3.7099416000, 3.6258468000e-3, -6.0956346000e-3 /)
  LoveDat(1:4,253) = (/ 252.0, -3.7147645000, 3.6127027000e-3, -6.0766330000e-3 /)
  LoveDat(1:4,254) = (/ 253.0, -3.7195854000, 3.5996869000e-3, -6.0578067000e-3 /)
  LoveDat(1:4,255) = (/ 254.0, -3.7244043000, 3.5867979000e-3, -6.0391534000e-3 /)
  LoveDat(1:4,256) = (/ 255.0, -3.7292211000, 3.5740342000e-3, -6.0206710000e-3 /)
  LoveDat(1:4,257) = (/ 256.0, -3.7340357000, 3.5613944000e-3, -6.0023572000e-3 /)
  LoveDat(1:4,258) = (/ 257.0, -3.7388483000, 3.5488772000e-3, -5.9842098000e-3 /)
  LoveDat(1:4,259) = (/ 258.0, -3.7436587000, 3.5364810000e-3, -5.9662266000e-3 /)
  LoveDat(1:4,260) = (/ 259.0, -3.7484669000, 3.5242045000e-3, -5.9484056000e-3 /)
  LoveDat(1:4,261) = (/ 260.0, -3.7532729000, 3.5120464000e-3, -5.9307447000e-3 /)
  LoveDat(1:4,262) = (/ 261.0, -3.7580766000, 3.5000053000e-3, -5.9132419000e-3 /)
  LoveDat(1:4,263) = (/ 262.0, -3.7628780000, 3.4880799000e-3, -5.8958950000e-3 /)
  LoveDat(1:4,264) = (/ 263.0, -3.7676772000, 3.4762689000e-3, -5.8787022000e-3 /)
  LoveDat(1:4,265) = (/ 264.0, -3.7724740000, 3.4645710000e-3, -5.8616614000e-3 /)
  LoveDat(1:4,266) = (/ 265.0, -3.7772685000, 3.4529849000e-3, -5.8447709000e-3 /)
  LoveDat(1:4,267) = (/ 266.0, -3.7820605000, 3.4415093000e-3, -5.8280285000e-3 /)
  LoveDat(1:4,268) = (/ 267.0, -3.7868501000, 3.4301431000e-3, -5.8114326000e-3 /)
  LoveDat(1:4,269) = (/ 268.0, -3.7916373000, 3.4188851000e-3, -5.7949812000e-3 /)
  LoveDat(1:4,270) = (/ 269.0, -3.7964220000, 3.4077339000e-3, -5.7786726000e-3 /)
  LoveDat(1:4,271) = (/ 270.0, -3.8012042000, 3.3966884000e-3, -5.7625050000e-3 /)
  LoveDat(1:4,272) = (/ 271.0, -3.8059839000, 3.3857475000e-3, -5.7464766000e-3 /)
  LoveDat(1:4,273) = (/ 272.0, -3.8107610000, 3.3749099000e-3, -5.7305857000e-3 /)
  LoveDat(1:4,274) = (/ 273.0, -3.8155355000, 3.3641746000e-3, -5.7148305000e-3 /)
  LoveDat(1:4,275) = (/ 274.0, -3.8203074000, 3.3535404000e-3, -5.6992095000e-3 /)
  LoveDat(1:4,276) = (/ 275.0, -3.8250766000, 3.3430061000e-3, -5.6837210000e-3 /)
  LoveDat(1:4,277) = (/ 276.0, -3.8298432000, 3.3325707000e-3, -5.6683633000e-3 /)
  LoveDat(1:4,278) = (/ 277.0, -3.8346070000, 3.3222331000e-3, -5.6531348000e-3 /)
  LoveDat(1:4,279) = (/ 278.0, -3.8393682000, 3.3119922000e-3, -5.6380340000e-3 /)
  LoveDat(1:4,280) = (/ 279.0, -3.8441265000, 3.3018470000e-3, -5.6230593000e-3 /)
  LoveDat(1:4,281) = (/ 280.0, -3.8488821000, 3.2917964000e-3, -5.6082092000e-3 /)
  LoveDat(1:4,282) = (/ 281.0, -3.8536348000, 3.2818393000e-3, -5.5934822000e-3 /)
  LoveDat(1:4,283) = (/ 282.0, -3.8583847000, 3.2719748000e-3, -5.5788767000e-3 /)
  LoveDat(1:4,284) = (/ 283.0, -3.8631317000, 3.2622018000e-3, -5.5643913000e-3 /)
  LoveDat(1:4,285) = (/ 284.0, -3.8678759000, 3.2525193000e-3, -5.5500246000e-3 /)
  LoveDat(1:4,286) = (/ 285.0, -3.8726170000, 3.2429264000e-3, -5.5357752000e-3 /)
  LoveDat(1:4,287) = (/ 286.0, -3.8773553000, 3.2334221000e-3, -5.5216416000e-3 /)
  LoveDat(1:4,288) = (/ 287.0, -3.8820905000, 3.2240054000e-3, -5.5076224000e-3 /)
  LoveDat(1:4,289) = (/ 288.0, -3.8868227000, 3.2146753000e-3, -5.4937164000e-3 /)
  LoveDat(1:4,290) = (/ 289.0, -3.8915519000, 3.2054310000e-3, -5.4799221000e-3 /)
  LoveDat(1:4,291) = (/ 290.0, -3.8962780000, 3.1962715000e-3, -5.4662383000e-3 /)
  LoveDat(1:4,292) = (/ 291.0, -3.9010010000, 3.1871958000e-3, -5.4526635000e-3 /)
  LoveDat(1:4,293) = (/ 292.0, -3.9057209000, 3.1782032000e-3, -5.4391967000e-3 /)
  LoveDat(1:4,294) = (/ 293.0, -3.9104377000, 3.1692926000e-3, -5.4258363000e-3 /)
  LoveDat(1:4,295) = (/ 294.0, -3.9151512000, 3.1604632000e-3, -5.4125813000e-3 /)
  LoveDat(1:4,296) = (/ 295.0, -3.9198616000, 3.1517142000e-3, -5.3994305000e-3 /)
  LoveDat(1:4,297) = (/ 296.0, -3.9245687000, 3.1430446000e-3, -5.3863824000e-3 /)
  LoveDat(1:4,298) = (/ 297.0, -3.9292725000, 3.1344537000e-3, -5.3734361000e-3 /)
  LoveDat(1:4,299) = (/ 298.0, -3.9339731000, 3.1259405000e-3, -5.3605902000e-3 /)
  LoveDat(1:4,300) = (/ 299.0, -3.9386704000, 3.1175043000e-3, -5.3478437000e-3 /)
  LoveDat(1:4,301) = (/ 300.0, -3.9433643000, 3.1091442000e-3, -5.3351954000e-3 /)
  LoveDat(1:4,302) = (/ 301.0, -3.9480548000, 3.1008594000e-3, -5.3226441000e-3 /)
  LoveDat(1:4,303) = (/ 302.0, -3.9527420000, 3.0926491000e-3, -5.3101888000e-3 /)
  LoveDat(1:4,304) = (/ 303.0, -3.9574257000, 3.0845126000e-3, -5.2978283000e-3 /)
  LoveDat(1:4,305) = (/ 304.0, -3.9621060000, 3.0764490000e-3, -5.2855615000e-3 /)
  LoveDat(1:4,306) = (/ 305.0, -3.9667828000, 3.0684575000e-3, -5.2733874000e-3 /)
  LoveDat(1:4,307) = (/ 306.0, -3.9714561000, 3.0605375000e-3, -5.2613050000e-3 /)
  LoveDat(1:4,308) = (/ 307.0, -3.9761259000, 3.0526881000e-3, -5.2493131000e-3 /)
  LoveDat(1:4,309) = (/ 308.0, -3.9807921000, 3.0449085000e-3, -5.2374107000e-3 /)
  LoveDat(1:4,310) = (/ 309.0, -3.9854548000, 3.0371982000e-3, -5.2255969000e-3 /)
  LoveDat(1:4,311) = (/ 310.0, -3.9901138000, 3.0295562000e-3, -5.2138707000e-3 /)
  LoveDat(1:4,312) = (/ 311.0, -3.9947693000, 3.0219820000e-3, -5.2022310000e-3 /)
  LoveDat(1:4,313) = (/ 312.0, -3.9994210000, 3.0144747000e-3, -5.1906768000e-3 /)
  LoveDat(1:4,314) = (/ 313.0, -4.0040691000, 3.0070337000e-3, -5.1792073000e-3 /)
  LoveDat(1:4,315) = (/ 314.0, -4.0087135000, 2.9996584000e-3, -5.1678215000e-3 /)
  LoveDat(1:4,316) = (/ 315.0, -4.0133542000, 2.9923479000e-3, -5.1565183000e-3 /)
  LoveDat(1:4,317) = (/ 316.0, -4.0179911000, 2.9851016000e-3, -5.1452970000e-3 /)
  LoveDat(1:4,318) = (/ 317.0, -4.0226242000, 2.9779189000e-3, -5.1341566000e-3 /)
  LoveDat(1:4,319) = (/ 318.0, -4.0272535000, 2.9707990000e-3, -5.1230962000e-3 /)
  LoveDat(1:4,320) = (/ 319.0, -4.0318790000, 2.9637414000e-3, -5.1121150000e-3 /)
  LoveDat(1:4,321) = (/ 320.0, -4.0365006000, 2.9567453000e-3, -5.1012119000e-3 /)
  LoveDat(1:4,322) = (/ 321.0, -4.0411184000, 2.9498101000e-3, -5.0903863000e-3 /)
  LoveDat(1:4,323) = (/ 322.0, -4.0457322000, 2.9429353000e-3, -5.0796372000e-3 /)
  LoveDat(1:4,324) = (/ 323.0, -4.0503421000, 2.9361201000e-3, -5.0689638000e-3 /)
  LoveDat(1:4,325) = (/ 324.0, -4.0549481000, 2.9293639000e-3, -5.0583652000e-3 /)
  LoveDat(1:4,326) = (/ 325.0, -4.0595501000, 2.9226662000e-3, -5.0478407000e-3 /)
  LoveDat(1:4,327) = (/ 326.0, -4.0641480000, 2.9160263000e-3, -5.0373894000e-3 /)
  LoveDat(1:4,328) = (/ 327.0, -4.0687420000, 2.9094435000e-3, -5.0270106000e-3 /)
  LoveDat(1:4,329) = (/ 328.0, -4.0733319000, 2.9029174000e-3, -5.0167034000e-3 /)
  LoveDat(1:4,330) = (/ 329.0, -4.0779177000, 2.8964474000e-3, -5.0064671000e-3 /)
  LoveDat(1:4,331) = (/ 330.0, -4.0824995000, 2.8900327000e-3, -4.9963009000e-3 /)
  LoveDat(1:4,332) = (/ 331.0, -4.0870771000, 2.8836730000e-3, -4.9862041000e-3 /)
  LoveDat(1:4,333) = (/ 332.0, -4.0916505000, 2.8773676000e-3, -4.9761758000e-3 /)
  LoveDat(1:4,334) = (/ 333.0, -4.0962198000, 2.8711159000e-3, -4.9662155000e-3 /)
  LoveDat(1:4,335) = (/ 334.0, -4.1007850000, 2.8649173000e-3, -4.9563223000e-3 /)
  LoveDat(1:4,336) = (/ 335.0, -4.1053459000, 2.8587715000e-3, -4.9464955000e-3 /)
  LoveDat(1:4,337) = (/ 336.0, -4.1099025000, 2.8526777000e-3, -4.9367344000e-3 /)
  LoveDat(1:4,338) = (/ 337.0, -4.1144549000, 2.8466354000e-3, -4.9270384000e-3 /)
  LoveDat(1:4,339) = (/ 338.0, -4.1190030000, 2.8406442000e-3, -4.9174066000e-3 /)
  LoveDat(1:4,340) = (/ 339.0, -4.1235469000, 2.8347035000e-3, -4.9078386000e-3 /)
  LoveDat(1:4,341) = (/ 340.0, -4.1280863000, 2.8288128000e-3, -4.8983335000e-3 /)
  LoveDat(1:4,342) = (/ 341.0, -4.1326215000, 2.8229715000e-3, -4.8888907000e-3 /)
  LoveDat(1:4,343) = (/ 342.0, -4.1371523000, 2.8171792000e-3, -4.8795095000e-3 /)
  LoveDat(1:4,344) = (/ 343.0, -4.1416786000, 2.8114353000e-3, -4.8701893000e-3 /)
  LoveDat(1:4,345) = (/ 344.0, -4.1462006000, 2.8057394000e-3, -4.8609295000e-3 /)
  LoveDat(1:4,346) = (/ 345.0, -4.1507181000, 2.8000909000e-3, -4.8517295000e-3 /)
  LoveDat(1:4,347) = (/ 346.0, -4.1552312000, 2.7944894000e-3, -4.8425885000e-3 /)
  LoveDat(1:4,348) = (/ 347.0, -4.1597397000, 2.7889344000e-3, -4.8335060000e-3 /)
  LoveDat(1:4,349) = (/ 348.0, -4.1642438000, 2.7834254000e-3, -4.8244814000e-3 /)
  LoveDat(1:4,350) = (/ 349.0, -4.1687434000, 2.7779620000e-3, -4.8155141000e-3 /)
  LoveDat(1:4,351) = (/ 350.0, -4.1732384000, 2.7725436000e-3, -4.8066034000e-3 /)
  LoveDat(1:4,352) = (/ 351.0, -4.1777288000, 2.7671698000e-3, -4.7977488000e-3 /)
  LoveDat(1:4,353) = (/ 352.0, -4.1822147000, 2.7618402000e-3, -4.7889498000e-3 /)
  LoveDat(1:4,354) = (/ 353.0, -4.1866959000, 2.7565543000e-3, -4.7802057000e-3 /)
  LoveDat(1:4,355) = (/ 354.0, -4.1911725000, 2.7513117000e-3, -4.7715160000e-3 /)
  LoveDat(1:4,356) = (/ 355.0, -4.1956445000, 2.7461118000e-3, -4.7628800000e-3 /)
  LoveDat(1:4,357) = (/ 356.0, -4.2001118000, 2.7409544000e-3, -4.7542974000e-3 /)
  LoveDat(1:4,358) = (/ 357.0, -4.2045744000, 2.7358388000e-3, -4.7457675000e-3 /)
  LoveDat(1:4,359) = (/ 358.0, -4.2090323000, 2.7307648000e-3, -4.7372897000e-3 /)
  LoveDat(1:4,360) = (/ 359.0, -4.2134854000, 2.7257319000e-3, -4.7288636000e-3 /)
  LoveDat(1:4,361) = (/ 360.0, -4.2179338000, 2.7207397000e-3, -4.7204886000e-3 /)
  LoveDat(1:4,362) = (/ 361.0, -4.2223775000, 2.7157877000e-3, -4.7121643000e-3 /)
  LoveDat(1:4,363) = (/ 362.0, -4.2268163000, 2.7108756000e-3, -4.7038900000e-3 /)
  LoveDat(1:4,364) = (/ 363.0, -4.2312503000, 2.7060029000e-3, -4.6956653000e-3 /)
  LoveDat(1:4,365) = (/ 364.0, -4.2356795000, 2.7011692000e-3, -4.6874897000e-3 /)
  LoveDat(1:4,366) = (/ 365.0, -4.2401039000, 2.6963742000e-3, -4.6793627000e-3 /)
  LoveDat(1:4,367) = (/ 366.0, -4.2445234000, 2.6916175000e-3, -4.6712838000e-3 /)
  LoveDat(1:4,368) = (/ 367.0, -4.2489380000, 2.6868986000e-3, -4.6632526000e-3 /)
  LoveDat(1:4,369) = (/ 368.0, -4.2533476000, 2.6822172000e-3, -4.6552684000e-3 /)
  LoveDat(1:4,370) = (/ 369.0, -4.2577524000, 2.6775728000e-3, -4.6473310000e-3 /)
  LoveDat(1:4,371) = (/ 370.0, -4.2621522000, 2.6729652000e-3, -4.6394397000e-3 /)
  LoveDat(1:4,372) = (/ 371.0, -4.2665470000, 2.6683940000e-3, -4.6315942000e-3 /)
  LoveDat(1:4,373) = (/ 372.0, -4.2709369000, 2.6638587000e-3, -4.6237940000e-3 /)
  LoveDat(1:4,374) = (/ 373.0, -4.2753218000, 2.6593590000e-3, -4.6160387000e-3 /)
  LoveDat(1:4,375) = (/ 374.0, -4.2797016000, 2.6548946000e-3, -4.6083277000e-3 /)
  LoveDat(1:4,376) = (/ 375.0, -4.2840764000, 2.6504651000e-3, -4.6006607000e-3 /)
  LoveDat(1:4,377) = (/ 376.0, -4.2884462000, 2.6460701000e-3, -4.5930373000e-3 /)
  LoveDat(1:4,378) = (/ 377.0, -4.2928108000, 2.6417093000e-3, -4.5854569000e-3 /)
  LoveDat(1:4,379) = (/ 378.0, -4.2971704000, 2.6373823000e-3, -4.5779192000e-3 /)
  LoveDat(1:4,380) = (/ 379.0, -4.3015249000, 2.6330888000e-3, -4.5704238000e-3 /)
  LoveDat(1:4,381) = (/ 380.0, -4.3058742000, 2.6288285000e-3, -4.5629702000e-3 /)
  LoveDat(1:4,382) = (/ 381.0, -4.3102184000, 2.6246011000e-3, -4.5555581000e-3 /)
  LoveDat(1:4,383) = (/ 382.0, -4.3145575000, 2.6204061000e-3, -4.5481870000e-3 /)
  LoveDat(1:4,384) = (/ 383.0, -4.3188914000, 2.6162432000e-3, -4.5408565000e-3 /)
  LoveDat(1:4,385) = (/ 384.0, -4.3232200000, 2.6121122000e-3, -4.5335663000e-3 /)
  LoveDat(1:4,386) = (/ 385.0, -4.3275435000, 2.6080128000e-3, -4.5263159000e-3 /)
  LoveDat(1:4,387) = (/ 386.0, -4.3318617000, 2.6039445000e-3, -4.5191050000e-3 /)
  LoveDat(1:4,388) = (/ 387.0, -4.3361747000, 2.5999071000e-3, -4.5119331000e-3 /)
  LoveDat(1:4,389) = (/ 388.0, -4.3404824000, 2.5959002000e-3, -4.5048000000e-3 /)
  LoveDat(1:4,390) = (/ 389.0, -4.3447848000, 2.5919236000e-3, -4.4977052000e-3 /)
  LoveDat(1:4,391) = (/ 390.0, -4.3490820000, 2.5879770000e-3, -4.4906484000e-3 /)
  LoveDat(1:4,392) = (/ 391.0, -4.3533738000, 2.5840600000e-3, -4.4836292000e-3 /)
  LoveDat(1:4,393) = (/ 392.0, -4.3576603000, 2.5801724000e-3, -4.4766472000e-3 /)
  LoveDat(1:4,394) = (/ 393.0, -4.3619414000, 2.5763138000e-3, -4.4697021000e-3 /)
  LoveDat(1:4,395) = (/ 394.0, -4.3662172000, 2.5724840000e-3, -4.4627935000e-3 /)
  LoveDat(1:4,396) = (/ 395.0, -4.3704876000, 2.5686827000e-3, -4.4559212000e-3 /)
  LoveDat(1:4,397) = (/ 396.0, -4.3747527000, 2.5649095000e-3, -4.4490846000e-3 /)
  LoveDat(1:4,398) = (/ 397.0, -4.3790123000, 2.5611642000e-3, -4.4422836000e-3 /)
  LoveDat(1:4,399) = (/ 398.0, -4.3832665000, 2.5574466000e-3, -4.4355178000e-3 /)
  LoveDat(1:4,400) = (/ 399.0, -4.3875152000, 2.5537563000e-3, -4.4287868000e-3 /)
  LoveDat(1:4,401) = (/ 400.0, -4.3917586000, 2.5500930000e-3, -4.4220903000e-3 /)
  LoveDat(1:4,402) = (/ 401.0, -4.3959964000, 2.5464565000e-3, -4.4154280000e-3 /)
  LoveDat(1:4,403) = (/ 402.0, -4.4002288000, 2.5428466000e-3, -4.4087995000e-3 /)
  LoveDat(1:4,404) = (/ 403.0, -4.4044556000, 2.5392629000e-3, -4.4022046000e-3 /)
  LoveDat(1:4,405) = (/ 404.0, -4.4086770000, 2.5357051000e-3, -4.3956430000e-3 /)
  LoveDat(1:4,406) = (/ 405.0, -4.4128928000, 2.5321731000e-3, -4.3891142000e-3 /)
  LoveDat(1:4,407) = (/ 406.0, -4.4171031000, 2.5286666000e-3, -4.3826181000e-3 /)
  LoveDat(1:4,408) = (/ 407.0, -4.4213078000, 2.5251852000e-3, -4.3761543000e-3 /)
  LoveDat(1:4,409) = (/ 408.0, -4.4255070000, 2.5217289000e-3, -4.3697225000e-3 /)
  LoveDat(1:4,410) = (/ 409.0, -4.4297006000, 2.5182972000e-3, -4.3633224000e-3 /)
  LoveDat(1:4,411) = (/ 410.0, -4.4338886000, 2.5148899000e-3, -4.3569537000e-3 /)
  LoveDat(1:4,412) = (/ 411.0, -4.4380709000, 2.5115069000e-3, -4.3506162000e-3 /)
  LoveDat(1:4,413) = (/ 412.0, -4.4422477000, 2.5081478000e-3, -4.3443095000e-3 /)
  LoveDat(1:4,414) = (/ 413.0, -4.4464188000, 2.5048125000e-3, -4.3380334000e-3 /)
  LoveDat(1:4,415) = (/ 414.0, -4.4505843000, 2.5015006000e-3, -4.3317876000e-3 /)
  LoveDat(1:4,416) = (/ 415.0, -4.4547441000, 2.4982119000e-3, -4.3255718000e-3 /)
  LoveDat(1:4,417) = (/ 416.0, -4.4588982000, 2.4949463000e-3, -4.3193857000e-3 /)
  LoveDat(1:4,418) = (/ 417.0, -4.4630466000, 2.4917034000e-3, -4.3132290000e-3 /)
  LoveDat(1:4,419) = (/ 418.0, -4.4671894000, 2.4884831000e-3, -4.3071016000e-3 /)
  LoveDat(1:4,420) = (/ 419.0, -4.4713264000, 2.4852851000e-3, -4.3010031000e-3 /)
  LoveDat(1:4,421) = (/ 420.0, -4.4754577000, 2.4821092000e-3, -4.2949332000e-3 /)
  LoveDat(1:4,422) = (/ 421.0, -4.4795832000, 2.4789551000e-3, -4.2888918000e-3 /)
  LoveDat(1:4,423) = (/ 422.0, -4.4837030000, 2.4758227000e-3, -4.2828785000e-3 /)
  LoveDat(1:4,424) = (/ 423.0, -4.4878171000, 2.4727118000e-3, -4.2768931000e-3 /)
  LoveDat(1:4,425) = (/ 424.0, -4.4919253000, 2.4696220000e-3, -4.2709353000e-3 /)
  LoveDat(1:4,426) = (/ 425.0, -4.4960278000, 2.4665532000e-3, -4.2650050000e-3 /)
  LoveDat(1:4,427) = (/ 426.0, -4.5001245000, 2.4635053000e-3, -4.2591017000e-3 /)
  LoveDat(1:4,428) = (/ 427.0, -4.5042153000, 2.4604778000e-3, -4.2532254000e-3 /)
  LoveDat(1:4,429) = (/ 428.0, -4.5083003000, 2.4574708000e-3, -4.2473758000e-3 /)
  LoveDat(1:4,430) = (/ 429.0, -4.5123795000, 2.4544839000e-3, -4.2415526000e-3 /)
  LoveDat(1:4,431) = (/ 430.0, -4.5164529000, 2.4515170000e-3, -4.2357555000e-3 /)
  LoveDat(1:4,432) = (/ 431.0, -4.5205204000, 2.4485699000e-3, -4.2299844000e-3 /)
  LoveDat(1:4,433) = (/ 432.0, -4.5245820000, 2.4456423000e-3, -4.2242391000e-3 /)
  LoveDat(1:4,434) = (/ 433.0, -4.5286377000, 2.4427340000e-3, -4.2185193000e-3 /)
  LoveDat(1:4,435) = (/ 434.0, -4.5326876000, 2.4398450000e-3, -4.2128247000e-3 /)
  LoveDat(1:4,436) = (/ 435.0, -4.5367315000, 2.4369749000e-3, -4.2071552000e-3 /)
  LoveDat(1:4,437) = (/ 436.0, -4.5407695000, 2.4341235000e-3, -4.2015105000e-3 /)
  LoveDat(1:4,438) = (/ 437.0, -4.5448016000, 2.4312908000e-3, -4.1958904000e-3 /)
  LoveDat(1:4,439) = (/ 438.0, -4.5488278000, 2.4284765000e-3, -4.1902947000e-3 /)
  LoveDat(1:4,440) = (/ 439.0, -4.5528480000, 2.4256804000e-3, -4.1847233000e-3 /)
  LoveDat(1:4,441) = (/ 440.0, -4.5568623000, 2.4229023000e-3, -4.1791757000e-3 /)
  LoveDat(1:4,442) = (/ 441.0, -4.5608706000, 2.4201420000e-3, -4.1736520000e-3 /)
  LoveDat(1:4,443) = (/ 442.0, -4.5648729000, 2.4173995000e-3, -4.1681518000e-3 /)
  LoveDat(1:4,444) = (/ 443.0, -4.5688693000, 2.4146744000e-3, -4.1626750000e-3 /)
  LoveDat(1:4,445) = (/ 444.0, -4.5728596000, 2.4119666000e-3, -4.1572213000e-3 /)
  LoveDat(1:4,446) = (/ 445.0, -4.5768440000, 2.4092760000e-3, -4.1517905000e-3 /)
  LoveDat(1:4,447) = (/ 446.0, -4.5808223000, 2.4066023000e-3, -4.1463825000e-3 /)
  LoveDat(1:4,448) = (/ 447.0, -4.5847946000, 2.4039454000e-3, -4.1409971000e-3 /)
  LoveDat(1:4,449) = (/ 448.0, -4.5887608000, 2.4013051000e-3, -4.1356340000e-3 /)
  LoveDat(1:4,450) = (/ 449.0, -4.5927211000, 2.3986813000e-3, -4.1302931000e-3 /)
  LoveDat(1:4,451) = (/ 450.0, -4.5966752000, 2.3960738000e-3, -4.1249742000e-3 /)
  LoveDat(1:4,452) = (/ 451.0, -4.6006234000, 2.3934824000e-3, -4.1196771000e-3 /)
  LoveDat(1:4,453) = (/ 452.0, -4.6045654000, 2.3909070000e-3, -4.1144015000e-3 /)
  LoveDat(1:4,454) = (/ 453.0, -4.6085014000, 2.3883473000e-3, -4.1091474000e-3 /)
  LoveDat(1:4,455) = (/ 454.0, -4.6124313000, 2.3858033000e-3, -4.1039146000e-3 /)
  LoveDat(1:4,456) = (/ 455.0, -4.6163550000, 2.3832748000e-3, -4.0987028000e-3 /)
  LoveDat(1:4,457) = (/ 456.0, -4.6202727000, 2.3807615000e-3, -4.0935118000e-3 /)
  LoveDat(1:4,458) = (/ 457.0, -4.6241843000, 2.3782635000e-3, -4.0883416000e-3 /)
  LoveDat(1:4,459) = (/ 458.0, -4.6280897000, 2.3757804000e-3, -4.0831919000e-3 /)
  LoveDat(1:4,460) = (/ 459.0, -4.6319890000, 2.3733122000e-3, -4.0780626000e-3 /)
  LoveDat(1:4,461) = (/ 460.0, -4.6358822000, 2.3708588000e-3, -4.0729534000e-3 /)
  LoveDat(1:4,462) = (/ 461.0, -4.6397692000, 2.3684198000e-3, -4.0678643000e-3 /)
  LoveDat(1:4,463) = (/ 462.0, -4.6436501000, 2.3659953000e-3, -4.0627950000e-3 /)
  LoveDat(1:4,464) = (/ 463.0, -4.6475249000, 2.3635851000e-3, -4.0577454000e-3 /)
  LoveDat(1:4,465) = (/ 464.0, -4.6513934000, 2.3611889000e-3, -4.0527153000e-3 /)
  LoveDat(1:4,466) = (/ 465.0, -4.6552558000, 2.3588068000e-3, -4.0477046000e-3 /)
  LoveDat(1:4,467) = (/ 466.0, -4.6591120000, 2.3564384000e-3, -4.0427131000e-3 /)
  LoveDat(1:4,468) = (/ 467.0, -4.6629620000, 2.3540838000e-3, -4.0377406000e-3 /)
  LoveDat(1:4,469) = (/ 468.0, -4.6668058000, 2.3517427000e-3, -4.0327870000e-3 /)
  LoveDat(1:4,470) = (/ 469.0, -4.6706434000, 2.3494150000e-3, -4.0278521000e-3 /)
  LoveDat(1:4,471) = (/ 470.0, -4.6744748000, 2.3471006000e-3, -4.0229358000e-3 /)
  LoveDat(1:4,472) = (/ 471.0, -4.6783000000, 2.3447994000e-3, -4.0180379000e-3 /)
  LoveDat(1:4,473) = (/ 472.0, -4.6821189000, 2.3425111000e-3, -4.0131582000e-3 /)
  LoveDat(1:4,474) = (/ 473.0, -4.6859316000, 2.3402357000e-3, -4.0082967000e-3 /)
  LoveDat(1:4,475) = (/ 474.0, -4.6897381000, 2.3379731000e-3, -4.0034532000e-3 /)
  LoveDat(1:4,476) = (/ 475.0, -4.6935383000, 2.3357231000e-3, -3.9986274000e-3 /)
  LoveDat(1:4,477) = (/ 476.0, -4.6973323000, 2.3334855000e-3, -3.9938194000e-3 /)
  LoveDat(1:4,478) = (/ 477.0, -4.7011201000, 2.3312604000e-3, -3.9890289000e-3 /)
  LoveDat(1:4,479) = (/ 478.0, -4.7049015000, 2.3290474000e-3, -3.9842557000e-3 /)
  LoveDat(1:4,480) = (/ 479.0, -4.7086767000, 2.3268466000e-3, -3.9794999000e-3 /)
  LoveDat(1:4,481) = (/ 480.0, -4.7124456000, 2.3246577000e-3, -3.9747611000e-3 /)
  LoveDat(1:4,482) = (/ 481.0, -4.7162083000, 2.3224807000e-3, -3.9700393000e-3 /)
  LoveDat(1:4,483) = (/ 482.0, -4.7199646000, 2.3203154000e-3, -3.9653344000e-3 /)
  LoveDat(1:4,484) = (/ 483.0, -4.7237147000, 2.3181618000e-3, -3.9606461000e-3 /)
  LoveDat(1:4,485) = (/ 484.0, -4.7274585000, 2.3160196000e-3, -3.9559744000e-3 /)
  LoveDat(1:4,486) = (/ 485.0, -4.7311959000, 2.3138889000e-3, -3.9513192000e-3 /)
  LoveDat(1:4,487) = (/ 486.0, -4.7349271000, 2.3117694000e-3, -3.9466802000e-3 /)
  LoveDat(1:4,488) = (/ 487.0, -4.7386519000, 2.3096610000e-3, -3.9420575000e-3 /)
  LoveDat(1:4,489) = (/ 488.0, -4.7423704000, 2.3075637000e-3, -3.9374508000e-3 /)
  LoveDat(1:4,490) = (/ 489.0, -4.7460826000, 2.3054773000e-3, -3.9328600000e-3 /)
  LoveDat(1:4,491) = (/ 490.0, -4.7497885000, 2.3034017000e-3, -3.9282850000e-3 /)
  LoveDat(1:4,492) = (/ 491.0, -4.7534880000, 2.3013368000e-3, -3.9237256000e-3 /)
  LoveDat(1:4,493) = (/ 492.0, -4.7571812000, 2.2992825000e-3, -3.9191818000e-3 /)
  LoveDat(1:4,494) = (/ 493.0, -4.7608681000, 2.2972386000e-3, -3.9146535000e-3 /)
  LoveDat(1:4,495) = (/ 494.0, -4.7645486000, 2.2952052000e-3, -3.9101404000e-3 /)
  LoveDat(1:4,496) = (/ 495.0, -4.7682227000, 2.2931820000e-3, -3.9056425000e-3 /)
  LoveDat(1:4,497) = (/ 496.0, -4.7718905000, 2.2911690000e-3, -3.9011597000e-3 /)
  LoveDat(1:4,498) = (/ 497.0, -4.7755520000, 2.2891660000e-3, -3.8966919000e-3 /)
  LoveDat(1:4,499) = (/ 498.0, -4.7792071000, 2.2871729000e-3, -3.8922389000e-3 /)
  LoveDat(1:4,500) = (/ 499.0, -4.7828558000, 2.2851898000e-3, -3.8878005000e-3 /)
  LoveDat(1:4,501) = (/ 500.0, -4.7864981000, 2.2832163000e-3, -3.8833768000e-3 /)
  LoveDat(1:4,502) = (/ 501.0, -4.7901341000, 2.2812525000e-3, -3.8789676000e-3 /)
  LoveDat(1:4,503) = (/ 502.0, -4.7937636000, 2.2792983000e-3, -3.8745728000e-3 /)
  LoveDat(1:4,504) = (/ 503.0, -4.7973868000, 2.2773535000e-3, -3.8701922000e-3 /)
  LoveDat(1:4,505) = (/ 504.0, -4.8010036000, 2.2754180000e-3, -3.8658258000e-3 /)
  LoveDat(1:4,506) = (/ 505.0, -4.8046141000, 2.2734918000e-3, -3.8614735000e-3 /)
  LoveDat(1:4,507) = (/ 506.0, -4.8082181000, 2.2715748000e-3, -3.8571351000e-3 /)
  LoveDat(1:4,508) = (/ 507.0, -4.8118157000, 2.2696668000e-3, -3.8528105000e-3 /)
  LoveDat(1:4,509) = (/ 508.0, -4.8154069000, 2.2677678000e-3, -3.8484997000e-3 /)
  LoveDat(1:4,510) = (/ 509.0, -4.8189918000, 2.2658777000e-3, -3.8442025000e-3 /)
  LoveDat(1:4,511) = (/ 510.0, -4.8225702000, 2.2639964000e-3, -3.8399188000e-3 /)
  LoveDat(1:4,512) = (/ 511.0, -4.8261422000, 2.2621237000e-3, -3.8356485000e-3 /)
  LoveDat(1:4,513) = (/ 512.0, -4.8297078000, 2.2602597000e-3, -3.8313916000e-3 /)
  LoveDat(1:4,514) = (/ 513.0, -4.8332670000, 2.2584041000e-3, -3.8271479000e-3 /)
  LoveDat(1:4,515) = (/ 514.0, -4.8368197000, 2.2565570000e-3, -3.8229173000e-3 /)
  LoveDat(1:4,516) = (/ 515.0, -4.8403661000, 2.2547183000e-3, -3.8186997000e-3 /)
  LoveDat(1:4,517) = (/ 516.0, -4.8439060000, 2.2528877000e-3, -3.8144951000e-3 /)
  LoveDat(1:4,518) = (/ 517.0, -4.8474395000, 2.2510654000e-3, -3.8103033000e-3 /)
  LoveDat(1:4,519) = (/ 518.0, -4.8509666000, 2.2492511000e-3, -3.8061243000e-3 /)
  LoveDat(1:4,520) = (/ 519.0, -4.8544872000, 2.2474448000e-3, -3.8019578000e-3 /)
  LoveDat(1:4,521) = (/ 520.0, -4.8580014000, 2.2456465000e-3, -3.7978040000e-3 /)
  LoveDat(1:4,522) = (/ 521.0, -4.8615092000, 2.2438560000e-3, -3.7936626000e-3 /)
  LoveDat(1:4,523) = (/ 522.0, -4.8650105000, 2.2420732000e-3, -3.7895335000e-3 /)
  LoveDat(1:4,524) = (/ 523.0, -4.8685054000, 2.2402981000e-3, -3.7854168000e-3 /)
  LoveDat(1:4,525) = (/ 524.0, -4.8719939000, 2.2385305000e-3, -3.7813122000e-3 /)
  LoveDat(1:4,526) = (/ 525.0, -4.8754759000, 2.2367705000e-3, -3.7772197000e-3 /)
  LoveDat(1:4,527) = (/ 526.0, -4.8789515000, 2.2350179000e-3, -3.7731392000e-3 /)
  LoveDat(1:4,528) = (/ 527.0, -4.8824206000, 2.2332727000e-3, -3.7690706000e-3 /)
  LoveDat(1:4,529) = (/ 528.0, -4.8858833000, 2.2315347000e-3, -3.7650139000e-3 /)
  LoveDat(1:4,530) = (/ 529.0, -4.8893395000, 2.2298040000e-3, -3.7609689000e-3 /)
  LoveDat(1:4,531) = (/ 530.0, -4.8927893000, 2.2280804000e-3, -3.7569356000e-3 /)
  LoveDat(1:4,532) = (/ 531.0, -4.8962327000, 2.2263638000e-3, -3.7529139000e-3 /)
  LoveDat(1:4,533) = (/ 532.0, -4.8996696000, 2.2246542000e-3, -3.7489037000e-3 /)
  LoveDat(1:4,534) = (/ 533.0, -4.9031000000, 2.2229515000e-3, -3.7449049000e-3 /)
  LoveDat(1:4,535) = (/ 534.0, -4.9065240000, 2.2212556000e-3, -3.7409174000e-3 /)
  LoveDat(1:4,536) = (/ 535.0, -4.9099415000, 2.2195665000e-3, -3.7369411000e-3 /)
  LoveDat(1:4,537) = (/ 536.0, -4.9133526000, 2.2178841000e-3, -3.7329761000e-3 /)
  LoveDat(1:4,538) = (/ 537.0, -4.9167573000, 2.2162082000e-3, -3.7290221000e-3 /)
  LoveDat(1:4,539) = (/ 538.0, -4.9201554000, 2.2145390000e-3, -3.7250792000e-3 /)
  LoveDat(1:4,540) = (/ 539.0, -4.9235472000, 2.2128762000e-3, -3.7211471000e-3 /)
  LoveDat(1:4,541) = (/ 540.0, -4.9269324000, 2.2112198000e-3, -3.7172260000e-3 /)
  LoveDat(1:4,542) = (/ 541.0, -4.9303112000, 2.2095698000e-3, -3.7133156000e-3 /)
  LoveDat(1:4,543) = (/ 542.0, -4.9336836000, 2.2079261000e-3, -3.7094160000e-3 /)
  LoveDat(1:4,544) = (/ 543.0, -4.9370495000, 2.2062885000e-3, -3.7055269000e-3 /)
  LoveDat(1:4,545) = (/ 544.0, -4.9404089000, 2.2046571000e-3, -3.7016485000e-3 /)
  LoveDat(1:4,546) = (/ 545.0, -4.9437619000, 2.2030318000e-3, -3.6977805000e-3 /)
  LoveDat(1:4,547) = (/ 546.0, -4.9471084000, 2.2014125000e-3, -3.6939229000e-3 /)
  LoveDat(1:4,548) = (/ 547.0, -4.9504485000, 2.1997991000e-3, -3.6900757000e-3 /)
  LoveDat(1:4,549) = (/ 548.0, -4.9537821000, 2.1981917000e-3, -3.6862387000e-3 /)
  LoveDat(1:4,550) = (/ 549.0, -4.9571092000, 2.1965901000e-3, -3.6824120000e-3 /)
  LoveDat(1:4,551) = (/ 550.0, -4.9604299000, 2.1949942000e-3, -3.6785954000e-3 /)
  LoveDat(1:4,552) = (/ 551.0, -4.9637442000, 2.1934040000e-3, -3.6747888000e-3 /)
  LoveDat(1:4,553) = (/ 552.0, -4.9670519000, 2.1918195000e-3, -3.6709922000e-3 /)
  LoveDat(1:4,554) = (/ 553.0, -4.9703533000, 2.1902406000e-3, -3.6672056000e-3 /)
  LoveDat(1:4,555) = (/ 554.0, -4.9736481000, 2.1886671000e-3, -3.6634288000e-3 /)
  LoveDat(1:4,556) = (/ 555.0, -4.9769366000, 2.1870992000e-3, -3.6596618000e-3 /)
  LoveDat(1:4,557) = (/ 556.0, -4.9802185000, 2.1855366000e-3, -3.6559045000e-3 /)
  LoveDat(1:4,558) = (/ 557.0, -4.9834940000, 2.1839795000e-3, -3.6521569000e-3 /)
  LoveDat(1:4,559) = (/ 558.0, -4.9867631000, 2.1824276000e-3, -3.6484189000e-3 /)
  LoveDat(1:4,560) = (/ 559.0, -4.9900257000, 2.1808809000e-3, -3.6446904000e-3 /)
  LoveDat(1:4,561) = (/ 560.0, -4.9932819000, 2.1793394000e-3, -3.6409714000e-3 /)
  LoveDat(1:4,562) = (/ 561.0, -4.9965316000, 2.1778031000e-3, -3.6372617000e-3 /)
  LoveDat(1:4,563) = (/ 562.0, -4.9997749000, 2.1762718000e-3, -3.6335615000e-3 /)
  LoveDat(1:4,564) = (/ 563.0, -5.0030117000, 2.1747455000e-3, -3.6298704000e-3 /)
  LoveDat(1:4,565) = (/ 564.0, -5.0062421000, 2.1732242000e-3, -3.6261887000e-3 /)
  LoveDat(1:4,566) = (/ 565.0, -5.0094660000, 2.1717078000e-3, -3.6225160000e-3 /)
  LoveDat(1:4,567) = (/ 566.0, -5.0126835000, 2.1701963000e-3, -3.6188525000e-3 /)
  LoveDat(1:4,568) = (/ 567.0, -5.0158946000, 2.1686895000e-3, -3.6151980000e-3 /)
  LoveDat(1:4,569) = (/ 568.0, -5.0190992000, 2.1671876000e-3, -3.6115525000e-3 /)
  LoveDat(1:4,570) = (/ 569.0, -5.0222974000, 2.1656903000e-3, -3.6079159000e-3 /)
  LoveDat(1:4,571) = (/ 570.0, -5.0254891000, 2.1641977000e-3, -3.6042882000e-3 /)
  LoveDat(1:4,572) = (/ 571.0, -5.0286744000, 2.1627096000e-3, -3.6006692000e-3 /)
  LoveDat(1:4,573) = (/ 572.0, -5.0318533000, 2.1612262000e-3, -3.5970590000e-3 /)
  LoveDat(1:4,574) = (/ 573.0, -5.0350258000, 2.1597472000e-3, -3.5934575000e-3 /)
  LoveDat(1:4,575) = (/ 574.0, -5.0381918000, 2.1582727000e-3, -3.5898647000e-3 /)
  LoveDat(1:4,576) = (/ 575.0, -5.0413514000, 2.1568026000e-3, -3.5862804000e-3 /)
  LoveDat(1:4,577) = (/ 576.0, -5.0445046000, 2.1553369000e-3, -3.5827047000e-3 /)
  LoveDat(1:4,578) = (/ 577.0, -5.0476514000, 2.1538755000e-3, -3.5791374000e-3 /)
  LoveDat(1:4,579) = (/ 578.0, -5.0507917000, 2.1524183000e-3, -3.5755785000e-3 /)
  LoveDat(1:4,580) = (/ 579.0, -5.0539256000, 2.1509654000e-3, -3.5720280000e-3 /)
  LoveDat(1:4,581) = (/ 580.0, -5.0570532000, 2.1495166000e-3, -3.5684858000e-3 /)
  LoveDat(1:4,582) = (/ 581.0, -5.0601743000, 2.1480720000e-3, -3.5649519000e-3 /)
  LoveDat(1:4,583) = (/ 582.0, -5.0632890000, 2.1466315000e-3, -3.5614262000e-3 /)
  LoveDat(1:4,584) = (/ 583.0, -5.0663973000, 2.1451950000e-3, -3.5579086000e-3 /)
  LoveDat(1:4,585) = (/ 584.0, -5.0694991000, 2.1437625000e-3, -3.5543992000e-3 /)
  LoveDat(1:4,586) = (/ 585.0, -5.0725946000, 2.1423339000e-3, -3.5508978000e-3 /)
  LoveDat(1:4,587) = (/ 586.0, -5.0756837000, 2.1409093000e-3, -3.5474044000e-3 /)
  LoveDat(1:4,588) = (/ 587.0, -5.0787664000, 2.1394885000e-3, -3.5439189000e-3 /)
  LoveDat(1:4,589) = (/ 588.0, -5.0818427000, 2.1380716000e-3, -3.5404414000e-3 /)
  LoveDat(1:4,590) = (/ 589.0, -5.0849126000, 2.1366585000e-3, -3.5369717000e-3 /)
  LoveDat(1:4,591) = (/ 590.0, -5.0879762000, 2.1352491000e-3, -3.5335099000e-3 /)
  LoveDat(1:4,592) = (/ 591.0, -5.0910333000, 2.1338434000e-3, -3.5300557000e-3 /)
  LoveDat(1:4,593) = (/ 592.0, -5.0940841000, 2.1324413000e-3, -3.5266094000e-3 /)
  LoveDat(1:4,594) = (/ 593.0, -5.0971285000, 2.1310429000e-3, -3.5231706000e-3 /)
  LoveDat(1:4,595) = (/ 594.0, -5.1001665000, 2.1296481000e-3, -3.5197395000e-3 /)
  LoveDat(1:4,596) = (/ 595.0, -5.1031982000, 2.1282569000e-3, -3.5163160000e-3 /)
  LoveDat(1:4,597) = (/ 596.0, -5.1062234000, 2.1268691000e-3, -3.5129000000e-3 /)
  LoveDat(1:4,598) = (/ 597.0, -5.1092424000, 2.1254848000e-3, -3.5094915000e-3 /)
  LoveDat(1:4,599) = (/ 598.0, -5.1122549000, 2.1241039000e-3, -3.5060904000e-3 /)
  LoveDat(1:4,600) = (/ 599.0, -5.1152611000, 2.1227265000e-3, -3.5026968000e-3 /)
  LoveDat(1:4,601) = (/ 600.0, -5.1182610000, 2.1213524000e-3, -3.4993105000e-3 /)
  LoveDat(1:4,602) = (/ 601.0, -5.1212545000, 2.1199816000e-3, -3.4959315000e-3 /)
  LoveDat(1:4,603) = (/ 602.0, -5.1242417000, 2.1186141000e-3, -3.4925597000e-3 /)
  LoveDat(1:4,604) = (/ 603.0, -5.1272225000, 2.1172498000e-3, -3.4891952000e-3 /)
  LoveDat(1:4,605) = (/ 604.0, -5.1301970000, 2.1158888000e-3, -3.4858379000e-3 /)
  LoveDat(1:4,606) = (/ 605.0, -5.1331651000, 2.1145309000e-3, -3.4824877000e-3 /)
  LoveDat(1:4,607) = (/ 606.0, -5.1361270000, 2.1131762000e-3, -3.4791445000e-3 /)
  LoveDat(1:4,608) = (/ 607.0, -5.1390825000, 2.1118246000e-3, -3.4758085000e-3 /)
  LoveDat(1:4,609) = (/ 608.0, -5.1420316000, 2.1104761000e-3, -3.4724795000e-3 /)
  LoveDat(1:4,610) = (/ 609.0, -5.1449745000, 2.1091306000e-3, -3.4691574000e-3 /)
  LoveDat(1:4,611) = (/ 610.0, -5.1479111000, 2.1077881000e-3, -3.4658423000e-3 /)
  LoveDat(1:4,612) = (/ 611.0, -5.1508413000, 2.1064486000e-3, -3.4625341000e-3 /)
  LoveDat(1:4,613) = (/ 612.0, -5.1537652000, 2.1051120000e-3, -3.4592327000e-3 /)
  LoveDat(1:4,614) = (/ 613.0, -5.1566829000, 2.1037784000e-3, -3.4559381000e-3 /)
  LoveDat(1:4,615) = (/ 614.0, -5.1595942000, 2.1024476000e-3, -3.4526504000e-3 /)
  LoveDat(1:4,616) = (/ 615.0, -5.1624993000, 2.1011196000e-3, -3.4493693000e-3 /)
  LoveDat(1:4,617) = (/ 616.0, -5.1653981000, 2.0997945000e-3, -3.4460950000e-3 /)
  LoveDat(1:4,618) = (/ 617.0, -5.1682905000, 2.0984722000e-3, -3.4428273000e-3 /)
  LoveDat(1:4,619) = (/ 618.0, -5.1711768000, 2.0971526000e-3, -3.4395663000e-3 /)
  LoveDat(1:4,620) = (/ 619.0, -5.1740567000, 2.0958358000e-3, -3.4363118000e-3 /)
  LoveDat(1:4,621) = (/ 620.0, -5.1769304000, 2.0945216000e-3, -3.4330639000e-3 /)
  LoveDat(1:4,622) = (/ 621.0, -5.1797978000, 2.0932101000e-3, -3.4298226000e-3 /)
  LoveDat(1:4,623) = (/ 622.0, -5.1826589000, 2.0919012000e-3, -3.4265877000e-3 /)
  LoveDat(1:4,624) = (/ 623.0, -5.1855138000, 2.0905950000e-3, -3.4233592000e-3 /)
  LoveDat(1:4,625) = (/ 624.0, -5.1883625000, 2.0892913000e-3, -3.4201372000e-3 /)
  LoveDat(1:4,626) = (/ 625.0, -5.1912049000, 2.0879902000e-3, -3.4169215000e-3 /)
  LoveDat(1:4,627) = (/ 626.0, -5.1940410000, 2.0866915000e-3, -3.4137122000e-3 /)
  LoveDat(1:4,628) = (/ 627.0, -5.1968710000, 2.0853954000e-3, -3.4105092000e-3 /)
  LoveDat(1:4,629) = (/ 628.0, -5.1996947000, 2.0841018000e-3, -3.4073124000e-3 /)
  LoveDat(1:4,630) = (/ 629.0, -5.2025121000, 2.0828105000e-3, -3.4041219000e-3 /)
  LoveDat(1:4,631) = (/ 630.0, -5.2053234000, 2.0815217000e-3, -3.4009376000e-3 /)
  LoveDat(1:4,632) = (/ 631.0, -5.2081285000, 2.0802353000e-3, -3.3977595000e-3 /)
  LoveDat(1:4,633) = (/ 632.0, -5.2109273000, 2.0789512000e-3, -3.3945875000e-3 /)
  LoveDat(1:4,634) = (/ 633.0, -5.2137199000, 2.0776695000e-3, -3.3914216000e-3 /)
  LoveDat(1:4,635) = (/ 634.0, -5.2165064000, 2.0763900000e-3, -3.3882618000e-3 /)
  LoveDat(1:4,636) = (/ 635.0, -5.2192866000, 2.0751129000e-3, -3.3851080000e-3 /)
  LoveDat(1:4,637) = (/ 636.0, -5.2220607000, 2.0738380000e-3, -3.3819602000e-3 /)
  LoveDat(1:4,638) = (/ 637.0, -5.2248286000, 2.0725653000e-3, -3.3788184000e-3 /)
  LoveDat(1:4,639) = (/ 638.0, -5.2275903000, 2.0712949000e-3, -3.3756826000e-3 /)
  LoveDat(1:4,640) = (/ 639.0, -5.2303458000, 2.0700266000e-3, -3.3725527000e-3 /)
  LoveDat(1:4,641) = (/ 640.0, -5.2330952000, 2.0687604000e-3, -3.3694286000e-3 /)
  LoveDat(1:4,642) = (/ 641.0, -5.2358384000, 2.0674964000e-3, -3.3663104000e-3 /)
  LoveDat(1:4,643) = (/ 642.0, -5.2385755000, 2.0662346000e-3, -3.3631981000e-3 /)
  LoveDat(1:4,644) = (/ 643.0, -5.2413064000, 2.0649747000e-3, -3.3600915000e-3 /)
  LoveDat(1:4,645) = (/ 644.0, -5.2440312000, 2.0637170000e-3, -3.3569907000e-3 /)
  LoveDat(1:4,646) = (/ 645.0, -5.2467498000, 2.0624613000e-3, -3.3538957000e-3 /)
  LoveDat(1:4,647) = (/ 646.0, -5.2494624000, 2.0612076000e-3, -3.3508063000e-3 /)
  LoveDat(1:4,648) = (/ 647.0, -5.2521688000, 2.0599559000e-3, -3.3477227000e-3 /)
  LoveDat(1:4,649) = (/ 648.0, -5.2548690000, 2.0587062000e-3, -3.3446446000e-3 /)
  LoveDat(1:4,650) = (/ 649.0, -5.2575632000, 2.0574585000e-3, -3.3415722000e-3 /)
  LoveDat(1:4,651) = (/ 650.0, -5.2602513000, 2.0562126000e-3, -3.3385054000e-3 /)
  LoveDat(1:4,652) = (/ 651.0, -5.2629332000, 2.0549687000e-3, -3.3354442000e-3 /)
  LoveDat(1:4,653) = (/ 652.0, -5.2656091000, 2.0537266000e-3, -3.3323885000e-3 /)
  LoveDat(1:4,654) = (/ 653.0, -5.2682789000, 2.0524865000e-3, -3.3293383000e-3 /)
  LoveDat(1:4,655) = (/ 654.0, -5.2709426000, 2.0512481000e-3, -3.3262936000e-3 /)
  LoveDat(1:4,656) = (/ 655.0, -5.2736002000, 2.0500116000e-3, -3.3232543000e-3 /)
  LoveDat(1:4,657) = (/ 656.0, -5.2762518000, 2.0487769000e-3, -3.3202205000e-3 /)
  LoveDat(1:4,658) = (/ 657.0, -5.2788973000, 2.0475440000e-3, -3.3171921000e-3 /)
  LoveDat(1:4,659) = (/ 658.0, -5.2815367000, 2.0463128000e-3, -3.3141691000e-3 /)
  LoveDat(1:4,660) = (/ 659.0, -5.2841701000, 2.0450834000e-3, -3.3111514000e-3 /)
  LoveDat(1:4,661) = (/ 660.0, -5.2867975000, 2.0438557000e-3, -3.3081390000e-3 /)
  LoveDat(1:4,662) = (/ 661.0, -5.2894188000, 2.0426297000e-3, -3.3051319000e-3 /)
  LoveDat(1:4,663) = (/ 662.0, -5.2920341000, 2.0414054000e-3, -3.3021302000e-3 /)
  LoveDat(1:4,664) = (/ 663.0, -5.2946433000, 2.0401828000e-3, -3.2991336000e-3 /)
  LoveDat(1:4,665) = (/ 664.0, -5.2972466000, 2.0389618000e-3, -3.2961423000e-3 /)
  LoveDat(1:4,666) = (/ 665.0, -5.2998438000, 2.0377425000e-3, -3.2931562000e-3 /)
  LoveDat(1:4,667) = (/ 666.0, -5.3024350000, 2.0365247000e-3, -3.2901753000e-3 /)
  LoveDat(1:4,668) = (/ 667.0, -5.3050203000, 2.0353086000e-3, -3.2871995000e-3 /)
  LoveDat(1:4,669) = (/ 668.0, -5.3075995000, 2.0340940000e-3, -3.2842288000e-3 /)
  LoveDat(1:4,670) = (/ 669.0, -5.3101728000, 2.0328810000e-3, -3.2812633000e-3 /)
  LoveDat(1:4,671) = (/ 670.0, -5.3127401000, 2.0316695000e-3, -3.2783028000e-3 /)
  LoveDat(1:4,672) = (/ 671.0, -5.3153014000, 2.0304596000e-3, -3.2753474000e-3 /)
  LoveDat(1:4,673) = (/ 672.0, -5.3178568000, 2.0292512000e-3, -3.2723970000e-3 /)
  LoveDat(1:4,674) = (/ 673.0, -5.3204062000, 2.0280443000e-3, -3.2694517000e-3 /)
  LoveDat(1:4,675) = (/ 674.0, -5.3229496000, 2.0268388000e-3, -3.2665113000e-3 /)
  LoveDat(1:4,676) = (/ 675.0, -5.3254871000, 2.0256348000e-3, -3.2635759000e-3 /)
  LoveDat(1:4,677) = (/ 676.0, -5.3280187000, 2.0244322000e-3, -3.2606454000e-3 /)
  LoveDat(1:4,678) = (/ 677.0, -5.3305444000, 2.0232311000e-3, -3.2577199000e-3 /)
  LoveDat(1:4,679) = (/ 678.0, -5.3330641000, 2.0220314000e-3, -3.2547992000e-3 /)
  LoveDat(1:4,680) = (/ 679.0, -5.3355779000, 2.0208331000e-3, -3.2518834000e-3 /)
  LoveDat(1:4,681) = (/ 680.0, -5.3380858000, 2.0196361000e-3, -3.2489725000e-3 /)
  LoveDat(1:4,682) = (/ 681.0, -5.3405878000, 2.0184406000e-3, -3.2460664000e-3 /)
  LoveDat(1:4,683) = (/ 682.0, -5.3430840000, 2.0172463000e-3, -3.2431652000e-3 /)
  LoveDat(1:4,684) = (/ 683.0, -5.3455742000, 2.0160534000e-3, -3.2402687000e-3 /)
  LoveDat(1:4,685) = (/ 684.0, -5.3480585000, 2.0148619000e-3, -3.2373770000e-3 /)
  LoveDat(1:4,686) = (/ 685.0, -5.3505370000, 2.0136716000e-3, -3.2344900000e-3 /)
  LoveDat(1:4,687) = (/ 686.0, -5.3530097000, 2.0124826000e-3, -3.2316078000e-3 /)
  LoveDat(1:4,688) = (/ 687.0, -5.3554764000, 2.0112949000e-3, -3.2287303000e-3 /)
  LoveDat(1:4,689) = (/ 688.0, -5.3579373000, 2.0101085000e-3, -3.2258574000e-3 /)
  LoveDat(1:4,690) = (/ 689.0, -5.3603924000, 2.0089233000e-3, -3.2229893000e-3 /)
  LoveDat(1:4,691) = (/ 690.0, -5.3628417000, 2.0077394000e-3, -3.2201258000e-3 /)
  LoveDat(1:4,692) = (/ 691.0, -5.3652851000, 2.0065567000e-3, -3.2172669000e-3 /)
  LoveDat(1:4,693) = (/ 692.0, -5.3677227000, 2.0053752000e-3, -3.2144126000e-3 /)
  LoveDat(1:4,694) = (/ 693.0, -5.3701545000, 2.0041948000e-3, -3.2115629000e-3 /)
  LoveDat(1:4,695) = (/ 694.0, -5.3725805000, 2.0030157000e-3, -3.2087178000e-3 /)
  LoveDat(1:4,696) = (/ 695.0, -5.3750006000, 2.0018377000e-3, -3.2058772000e-3 /)
  LoveDat(1:4,697) = (/ 696.0, -5.3774150000, 2.0006609000e-3, -3.2030412000e-3 /)
  LoveDat(1:4,698) = (/ 697.0, -5.3798237000, 1.9994853000e-3, -3.2002097000e-3 /)
  LoveDat(1:4,699) = (/ 698.0, -5.3822265000, 1.9983108000e-3, -3.1973826000e-3 /)
  LoveDat(1:4,700) = (/ 699.0, -5.3846236000, 1.9971373000e-3, -3.1945601000e-3 /)
  LoveDat(1:4,701) = (/ 700.0, -5.3870149000, 1.9959650000e-3, -3.1917420000e-3 /)
  LoveDat(1:4,702) = (/ 701.0, -5.3894005000, 1.9947938000e-3, -3.1889283000e-3 /)
  LoveDat(1:4,703) = (/ 702.0, -5.3917803000, 1.9936237000e-3, -3.1861191000e-3 /)
  LoveDat(1:4,704) = (/ 703.0, -5.3941544000, 1.9924547000e-3, -3.1833143000e-3 /)
  LoveDat(1:4,705) = (/ 704.0, -5.3965228000, 1.9912867000e-3, -3.1805139000e-3 /)
  LoveDat(1:4,706) = (/ 705.0, -5.3988854000, 1.9901198000e-3, -3.1777178000e-3 /)
  LoveDat(1:4,707) = (/ 706.0, -5.4012423000, 1.9889539000e-3, -3.1749261000e-3 /)
  LoveDat(1:4,708) = (/ 707.0, -5.4035936000, 1.9877890000e-3, -3.1721387000e-3 /)
  LoveDat(1:4,709) = (/ 708.0, -5.4059391000, 1.9866252000e-3, -3.1693556000e-3 /)
  LoveDat(1:4,710) = (/ 709.0, -5.4082790000, 1.9854623000e-3, -3.1665769000e-3 /)
  LoveDat(1:4,711) = (/ 710.0, -5.4106131000, 1.9843005000e-3, -3.1638024000e-3 /)
  LoveDat(1:4,712) = (/ 711.0, -5.4129416000, 1.9831396000e-3, -3.1610322000e-3 /)
  LoveDat(1:4,713) = (/ 712.0, -5.4152645000, 1.9819797000e-3, -3.1582662000e-3 /)
  LoveDat(1:4,714) = (/ 713.0, -5.4175816000, 1.9808208000e-3, -3.1555045000e-3 /)
  LoveDat(1:4,715) = (/ 714.0, -5.4198932000, 1.9796628000e-3, -3.1527469000e-3 /)
  LoveDat(1:4,716) = (/ 715.0, -5.4221991000, 1.9785058000e-3, -3.1499936000e-3 /)
  LoveDat(1:4,717) = (/ 716.0, -5.4244993000, 1.9773497000e-3, -3.1472445000e-3 /)
  LoveDat(1:4,718) = (/ 717.0, -5.4267939000, 1.9761945000e-3, -3.1444995000e-3 /)
  LoveDat(1:4,719) = (/ 718.0, -5.4290830000, 1.9750402000e-3, -3.1417587000e-3 /)
  LoveDat(1:4,720) = (/ 719.0, -5.4313664000, 1.9738869000e-3, -3.1390221000e-3 /)
  LoveDat(1:4,721) = (/ 720.0, -5.4336442000, 1.9727344000e-3, -3.1362895000e-3 /)
  LoveDat(1:4,722) = (/ 721.0, -5.4359164000, 1.9715828000e-3, -3.1335611000e-3 /)
  LoveDat(1:4,723) = (/ 722.0, -5.4381830000, 1.9704321000e-3, -3.1308367000e-3 /)
  LoveDat(1:4,724) = (/ 723.0, -5.4404441000, 1.9692823000e-3, -3.1281164000e-3 /)
  LoveDat(1:4,725) = (/ 724.0, -5.4426996000, 1.9681333000e-3, -3.1254002000e-3 /)
  LoveDat(1:4,726) = (/ 725.0, -5.4449495000, 1.9669852000e-3, -3.1226881000e-3 /)
  LoveDat(1:4,727) = (/ 726.0, -5.4471939000, 1.9658379000e-3, -3.1199799000e-3 /)
  LoveDat(1:4,728) = (/ 727.0, -5.4494328000, 1.9646915000e-3, -3.1172758000e-3 /)
  LoveDat(1:4,729) = (/ 728.0, -5.4516661000, 1.9635458000e-3, -3.1145757000e-3 /)
  LoveDat(1:4,730) = (/ 729.0, -5.4538938000, 1.9624010000e-3, -3.1118796000e-3 /)
  LoveDat(1:4,731) = (/ 730.0, -5.4561161000, 1.9612570000e-3, -3.1091874000e-3 /)
  LoveDat(1:4,732) = (/ 731.0, -5.4583329000, 1.9601138000e-3, -3.1064992000e-3 /)
  LoveDat(1:4,733) = (/ 732.0, -5.4605441000, 1.9589714000e-3, -3.1038149000e-3 /)
  LoveDat(1:4,734) = (/ 733.0, -5.4627499000, 1.9578298000e-3, -3.1011346000e-3 /)
  LoveDat(1:4,735) = (/ 734.0, -5.4649502000, 1.9566889000e-3, -3.0984582000e-3 /)
  LoveDat(1:4,736) = (/ 735.0, -5.4671450000, 1.9555488000e-3, -3.0957857000e-3 /)
  LoveDat(1:4,737) = (/ 736.0, -5.4693343000, 1.9544095000e-3, -3.0931171000e-3 /)
  LoveDat(1:4,738) = (/ 737.0, -5.4715182000, 1.9532709000e-3, -3.0904524000e-3 /)
  LoveDat(1:4,739) = (/ 738.0, -5.4736966000, 1.9521331000e-3, -3.0877915000e-3 /)
  LoveDat(1:4,740) = (/ 739.0, -5.4758696000, 1.9509960000e-3, -3.0851345000e-3 /)
  LoveDat(1:4,741) = (/ 740.0, -5.4780372000, 1.9498596000e-3, -3.0824813000e-3 /)
  LoveDat(1:4,742) = (/ 741.0, -5.4801993000, 1.9487240000e-3, -3.0798319000e-3 /)
  LoveDat(1:4,743) = (/ 742.0, -5.4823560000, 1.9475891000e-3, -3.0771864000e-3 /)
  LoveDat(1:4,744) = (/ 743.0, -5.4845073000, 1.9464549000e-3, -3.0745446000e-3 /)
  LoveDat(1:4,745) = (/ 744.0, -5.4866533000, 1.9453214000e-3, -3.0719066000e-3 /)
  LoveDat(1:4,746) = (/ 745.0, -5.4887938000, 1.9441885000e-3, -3.0692724000e-3 /)
  LoveDat(1:4,747) = (/ 746.0, -5.4909289000, 1.9430564000e-3, -3.0666420000e-3 /)
  LoveDat(1:4,748) = (/ 747.0, -5.4930587000, 1.9419250000e-3, -3.0640153000e-3 /)
  LoveDat(1:4,749) = (/ 748.0, -5.4951831000, 1.9407942000e-3, -3.0613923000e-3 /)
  LoveDat(1:4,750) = (/ 749.0, -5.4973021000, 1.9396641000e-3, -3.0587731000e-3 /)
  LoveDat(1:4,751) = (/ 750.0, -5.4994158000, 1.9385347000e-3, -3.0561575000e-3 /)
  LoveDat(1:4,752) = (/ 751.0, -5.5015242000, 1.9374059000e-3, -3.0535457000e-3 /)
  LoveDat(1:4,753) = (/ 752.0, -5.5036272000, 1.9362778000e-3, -3.0509375000e-3 /)
  LoveDat(1:4,754) = (/ 753.0, -5.5057250000, 1.9351503000e-3, -3.0483331000e-3 /)
  LoveDat(1:4,755) = (/ 754.0, -5.5078174000, 1.9340235000e-3, -3.0457322000e-3 /)
  LoveDat(1:4,756) = (/ 755.0, -5.5099044000, 1.9328973000e-3, -3.0431351000e-3 /)
  LoveDat(1:4,757) = (/ 756.0, -5.5119863000, 1.9317717000e-3, -3.0405415000e-3 /)
  LoveDat(1:4,758) = (/ 757.0, -5.5140628000, 1.9306468000e-3, -3.0379516000e-3 /)
  LoveDat(1:4,759) = (/ 758.0, -5.5161340000, 1.9295224000e-3, -3.0353653000e-3 /)
  LoveDat(1:4,760) = (/ 759.0, -5.5182000000, 1.9283987000e-3, -3.0327826000e-3 /)
  LoveDat(1:4,761) = (/ 760.0, -5.5202607000, 1.9272756000e-3, -3.0302035000e-3 /)
  LoveDat(1:4,762) = (/ 761.0, -5.5223161000, 1.9261531000e-3, -3.0276280000e-3 /)
  LoveDat(1:4,763) = (/ 762.0, -5.5243664000, 1.9250311000e-3, -3.0250561000e-3 /)
  LoveDat(1:4,764) = (/ 763.0, -5.5264113000, 1.9239098000e-3, -3.0224877000e-3 /)
  LoveDat(1:4,765) = (/ 764.0, -5.5284511000, 1.9227890000e-3, -3.0199228000e-3 /)
  LoveDat(1:4,766) = (/ 765.0, -5.5304856000, 1.9216689000e-3, -3.0173615000e-3 /)
  LoveDat(1:4,767) = (/ 766.0, -5.5325150000, 1.9205493000e-3, -3.0148038000e-3 /)
  LoveDat(1:4,768) = (/ 767.0, -5.5345391000, 1.9194303000e-3, -3.0122495000e-3 /)
  LoveDat(1:4,769) = (/ 768.0, -5.5365581000, 1.9183118000e-3, -3.0096987000e-3 /)
  LoveDat(1:4,770) = (/ 769.0, -5.5385718000, 1.9171939000e-3, -3.0071515000e-3 /)
  LoveDat(1:4,771) = (/ 770.0, -5.5405804000, 1.9160766000e-3, -3.0046077000e-3 /)
  LoveDat(1:4,772) = (/ 771.0, -5.5425839000, 1.9149598000e-3, -3.0020674000e-3 /)
  LoveDat(1:4,773) = (/ 772.0, -5.5445822000, 1.9138435000e-3, -2.9995305000e-3 /)
  LoveDat(1:4,774) = (/ 773.0, -5.5465753000, 1.9127278000e-3, -2.9969971000e-3 /)
  LoveDat(1:4,775) = (/ 774.0, -5.5485633000, 1.9116127000e-3, -2.9944671000e-3 /)
  LoveDat(1:4,776) = (/ 775.0, -5.5505462000, 1.9104981000e-3, -2.9919406000e-3 /)
  LoveDat(1:4,777) = (/ 776.0, -5.5525239000, 1.9093840000e-3, -2.9894175000e-3 /)
  LoveDat(1:4,778) = (/ 777.0, -5.5544966000, 1.9082704000e-3, -2.9868978000e-3 /)
  LoveDat(1:4,779) = (/ 778.0, -5.5564641000, 1.9071573000e-3, -2.9843815000e-3 /)
  LoveDat(1:4,780) = (/ 779.0, -5.5584266000, 1.9060448000e-3, -2.9818686000e-3 /)
  LoveDat(1:4,781) = (/ 780.0, -5.5603840000, 1.9049328000e-3, -2.9793591000e-3 /)
  LoveDat(1:4,782) = (/ 781.0, -5.5623363000, 1.9038213000e-3, -2.9768530000e-3 /)
  LoveDat(1:4,783) = (/ 782.0, -5.5642835000, 1.9027103000e-3, -2.9743502000e-3 /)
  LoveDat(1:4,784) = (/ 783.0, -5.5662257000, 1.9015998000e-3, -2.9718507000e-3 /)
  LoveDat(1:4,785) = (/ 784.0, -5.5681628000, 1.9004898000e-3, -2.9693547000e-3 /)
  LoveDat(1:4,786) = (/ 785.0, -5.5700949000, 1.8993803000e-3, -2.9668619000e-3 /)
  LoveDat(1:4,787) = (/ 786.0, -5.5720220000, 1.8982713000e-3, -2.9643725000e-3 /)
  LoveDat(1:4,788) = (/ 787.0, -5.5739440000, 1.8971627000e-3, -2.9618864000e-3 /)
  LoveDat(1:4,789) = (/ 788.0, -5.5758611000, 1.8960547000e-3, -2.9594036000e-3 /)
  LoveDat(1:4,790) = (/ 789.0, -5.5777731000, 1.8949471000e-3, -2.9569240000e-3 /)
  LoveDat(1:4,791) = (/ 790.0, -5.5796802000, 1.8938401000e-3, -2.9544478000e-3 /)
  LoveDat(1:4,792) = (/ 791.0, -5.5815822000, 1.8927334000e-3, -2.9519749000e-3 /)
  LoveDat(1:4,793) = (/ 792.0, -5.5834793000, 1.8916273000e-3, -2.9495052000e-3 /)
  LoveDat(1:4,794) = (/ 793.0, -5.5853715000, 1.8905216000e-3, -2.9470388000e-3 /)
  LoveDat(1:4,795) = (/ 794.0, -5.5872586000, 1.8894164000e-3, -2.9445756000e-3 /)
  LoveDat(1:4,796) = (/ 795.0, -5.5891409000, 1.8883117000e-3, -2.9421157000e-3 /)
  LoveDat(1:4,797) = (/ 796.0, -5.5910182000, 1.8872074000e-3, -2.9396591000e-3 /)
  LoveDat(1:4,798) = (/ 797.0, -5.5928905000, 1.8861036000e-3, -2.9372056000e-3 /)
  LoveDat(1:4,799) = (/ 798.0, -5.5947580000, 1.8850002000e-3, -2.9347554000e-3 /)
  LoveDat(1:4,800) = (/ 799.0, -5.5966205000, 1.8838973000e-3, -2.9323083000e-3 /)
  LoveDat(1:4,801) = (/ 800.0, -5.5984781000, 1.8827948000e-3, -2.9298645000e-3 /)
  LoveDat(1:4,802) = (/ 801.0, -5.6003309000, 1.8816928000e-3, -2.9274239000e-3 /)
  LoveDat(1:4,803) = (/ 802.0, -5.6021787000, 1.8805912000e-3, -2.9249864000e-3 /)
  LoveDat(1:4,804) = (/ 803.0, -5.6040217000, 1.8794901000e-3, -2.9225522000e-3 /)
  LoveDat(1:4,805) = (/ 804.0, -5.6058598000, 1.8783894000e-3, -2.9201211000e-3 /)
  LoveDat(1:4,806) = (/ 805.0, -5.6076931000, 1.8772891000e-3, -2.9176931000e-3 /)
  LoveDat(1:4,807) = (/ 806.0, -5.6095215000, 1.8761892000e-3, -2.9152683000e-3 /)
  LoveDat(1:4,808) = (/ 807.0, -5.6113451000, 1.8750898000e-3, -2.9128467000e-3 /)
  LoveDat(1:4,809) = (/ 808.0, -5.6131638000, 1.8739909000e-3, -2.9104282000e-3 /)
  LoveDat(1:4,810) = (/ 809.0, -5.6149777000, 1.8728923000e-3, -2.9080128000e-3 /)
  LoveDat(1:4,811) = (/ 810.0, -5.6167869000, 1.8717942000e-3, -2.9056005000e-3 /)
  LoveDat(1:4,812) = (/ 811.0, -5.6185912000, 1.8706965000e-3, -2.9031914000e-3 /)
  LoveDat(1:4,813) = (/ 812.0, -5.6203907000, 1.8695992000e-3, -2.9007853000e-3 /)
  LoveDat(1:4,814) = (/ 813.0, -5.6221855000, 1.8685023000e-3, -2.8983824000e-3 /)
  LoveDat(1:4,815) = (/ 814.0, -5.6239754000, 1.8674058000e-3, -2.8959825000e-3 /)
  LoveDat(1:4,816) = (/ 815.0, -5.6257606000, 1.8663098000e-3, -2.8935857000e-3 /)
  LoveDat(1:4,817) = (/ 816.0, -5.6275411000, 1.8652141000e-3, -2.8911920000e-3 /)
  LoveDat(1:4,818) = (/ 817.0, -5.6293168000, 1.8641189000e-3, -2.8888014000e-3 /)
  LoveDat(1:4,819) = (/ 818.0, -5.6310878000, 1.8630241000e-3, -2.8864138000e-3 /)
  LoveDat(1:4,820) = (/ 819.0, -5.6328540000, 1.8619297000e-3, -2.8840292000e-3 /)
  LoveDat(1:4,821) = (/ 820.0, -5.6346155000, 1.8608357000e-3, -2.8816477000e-3 /)
  LoveDat(1:4,822) = (/ 821.0, -5.6363723000, 1.8597420000e-3, -2.8792693000e-3 /)
  LoveDat(1:4,823) = (/ 822.0, -5.6381245000, 1.8586488000e-3, -2.8768939000e-3 /)
  LoveDat(1:4,824) = (/ 823.0, -5.6398719000, 1.8575560000e-3, -2.8745215000e-3 /)
  LoveDat(1:4,825) = (/ 824.0, -5.6416146000, 1.8564636000e-3, -2.8721521000e-3 /)
  LoveDat(1:4,826) = (/ 825.0, -5.6433527000, 1.8553716000e-3, -2.8697857000e-3 /)
  LoveDat(1:4,827) = (/ 826.0, -5.6450861000, 1.8542799000e-3, -2.8674223000e-3 /)
  LoveDat(1:4,828) = (/ 827.0, -5.6468149000, 1.8531887000e-3, -2.8650619000e-3 /)
  LoveDat(1:4,829) = (/ 828.0, -5.6485390000, 1.8520979000e-3, -2.8627045000e-3 /)
  LoveDat(1:4,830) = (/ 829.0, -5.6502584000, 1.8510074000e-3, -2.8603501000e-3 /)
  LoveDat(1:4,831) = (/ 830.0, -5.6519733000, 1.8499173000e-3, -2.8579986000e-3 /)
  LoveDat(1:4,832) = (/ 831.0, -5.6536835000, 1.8488276000e-3, -2.8556502000e-3 /)
  LoveDat(1:4,833) = (/ 832.0, -5.6553891000, 1.8477384000e-3, -2.8533046000e-3 /)
  LoveDat(1:4,834) = (/ 833.0, -5.6570901000, 1.8466494000e-3, -2.8509621000e-3 /)
  LoveDat(1:4,835) = (/ 834.0, -5.6587866000, 1.8455609000e-3, -2.8486224000e-3 /)
  LoveDat(1:4,836) = (/ 835.0, -5.6604784000, 1.8444728000e-3, -2.8462858000e-3 /)
  LoveDat(1:4,837) = (/ 836.0, -5.6621657000, 1.8433850000e-3, -2.8439520000e-3 /)
  LoveDat(1:4,838) = (/ 837.0, -5.6638484000, 1.8422976000e-3, -2.8416212000e-3 /)
  LoveDat(1:4,839) = (/ 838.0, -5.6655266000, 1.8412106000e-3, -2.8392933000e-3 /)
  LoveDat(1:4,840) = (/ 839.0, -5.6672002000, 1.8401239000e-3, -2.8369683000e-3 /)
  LoveDat(1:4,841) = (/ 840.0, -5.6688693000, 1.8390377000e-3, -2.8346462000e-3 /)
  LoveDat(1:4,842) = (/ 841.0, -5.6705338000, 1.8379518000e-3, -2.8323270000e-3 /)
  LoveDat(1:4,843) = (/ 842.0, -5.6721939000, 1.8368663000e-3, -2.8300107000e-3 /)
  LoveDat(1:4,844) = (/ 843.0, -5.6738494000, 1.8357811000e-3, -2.8276973000e-3 /)
  LoveDat(1:4,845) = (/ 844.0, -5.6755004000, 1.8346964000e-3, -2.8253867000e-3 /)
  LoveDat(1:4,846) = (/ 845.0, -5.6771470000, 1.8336120000e-3, -2.8230791000e-3 /)
  LoveDat(1:4,847) = (/ 846.0, -5.6787890000, 1.8325279000e-3, -2.8207743000e-3 /)
  LoveDat(1:4,848) = (/ 847.0, -5.6804266000, 1.8314443000e-3, -2.8184724000e-3 /)
  LoveDat(1:4,849) = (/ 848.0, -5.6820597000, 1.8303610000e-3, -2.8161733000e-3 /)
  LoveDat(1:4,850) = (/ 849.0, -5.6836884000, 1.8292781000e-3, -2.8138771000e-3 /)
  LoveDat(1:4,851) = (/ 850.0, -5.6853127000, 1.8281955000e-3, -2.8115837000e-3 /)
  LoveDat(1:4,852) = (/ 851.0, -5.6869325000, 1.8271133000e-3, -2.8092932000e-3 /)
  LoveDat(1:4,853) = (/ 852.0, -5.6885478000, 1.8260315000e-3, -2.8070055000e-3 /)
  LoveDat(1:4,854) = (/ 853.0, -5.6901588000, 1.8249501000e-3, -2.8047206000e-3 /)
  LoveDat(1:4,855) = (/ 854.0, -5.6917653000, 1.8238690000e-3, -2.8024385000e-3 /)
  LoveDat(1:4,856) = (/ 855.0, -5.6933675000, 1.8227882000e-3, -2.8001593000e-3 /)
  LoveDat(1:4,857) = (/ 856.0, -5.6949653000, 1.8217079000e-3, -2.7978829000e-3 /)
  LoveDat(1:4,858) = (/ 857.0, -5.6965586000, 1.8206279000e-3, -2.7956092000e-3 /)
  LoveDat(1:4,859) = (/ 858.0, -5.6981477000, 1.8195482000e-3, -2.7933384000e-3 /)
  LoveDat(1:4,860) = (/ 859.0, -5.6997323000, 1.8184690000e-3, -2.7910703000e-3 /)
  LoveDat(1:4,861) = (/ 860.0, -5.7013126000, 1.8173900000e-3, -2.7888051000e-3 /)
  LoveDat(1:4,862) = (/ 861.0, -5.7028886000, 1.8163115000e-3, -2.7865426000e-3 /)
  LoveDat(1:4,863) = (/ 862.0, -5.7044602000, 1.8152333000e-3, -2.7842829000e-3 /)
  LoveDat(1:4,864) = (/ 863.0, -5.7060275000, 1.8141555000e-3, -2.7820260000e-3 /)
  LoveDat(1:4,865) = (/ 864.0, -5.7075905000, 1.8130780000e-3, -2.7797718000e-3 /)
  LoveDat(1:4,866) = (/ 865.0, -5.7091492000, 1.8120009000e-3, -2.7775204000e-3 /)
  LoveDat(1:4,867) = (/ 866.0, -5.7107035000, 1.8109241000e-3, -2.7752717000e-3 /)
  LoveDat(1:4,868) = (/ 867.0, -5.7122536000, 1.8098477000e-3, -2.7730258000e-3 /)
  LoveDat(1:4,869) = (/ 868.0, -5.7137995000, 1.8087717000e-3, -2.7707826000e-3 /)
  LoveDat(1:4,870) = (/ 869.0, -5.7153410000, 1.8076960000e-3, -2.7685421000e-3 /)
  LoveDat(1:4,871) = (/ 870.0, -5.7168783000, 1.8066207000e-3, -2.7663044000e-3 /)
  LoveDat(1:4,872) = (/ 871.0, -5.7184113000, 1.8055458000e-3, -2.7640694000e-3 /)
  LoveDat(1:4,873) = (/ 872.0, -5.7199401000, 1.8044712000e-3, -2.7618372000e-3 /)
  LoveDat(1:4,874) = (/ 873.0, -5.7214646000, 1.8033969000e-3, -2.7596076000e-3 /)
  LoveDat(1:4,875) = (/ 874.0, -5.7229850000, 1.8023230000e-3, -2.7573808000e-3 /)
  LoveDat(1:4,876) = (/ 875.0, -5.7245011000, 1.8012495000e-3, -2.7551566000e-3 /)
  LoveDat(1:4,877) = (/ 876.0, -5.7260130000, 1.8001763000e-3, -2.7529352000e-3 /)
  LoveDat(1:4,878) = (/ 877.0, -5.7275207000, 1.7991035000e-3, -2.7507164000e-3 /)
  LoveDat(1:4,879) = (/ 878.0, -5.7290242000, 1.7980311000e-3, -2.7485003000e-3 /)
  LoveDat(1:4,880) = (/ 879.0, -5.7305236000, 1.7969590000e-3, -2.7462870000e-3 /)
  LoveDat(1:4,881) = (/ 880.0, -5.7320187000, 1.7958873000e-3, -2.7440763000e-3 /)
  LoveDat(1:4,882) = (/ 881.0, -5.7335097000, 1.7948159000e-3, -2.7418682000e-3 /)
  LoveDat(1:4,883) = (/ 882.0, -5.7349966000, 1.7937449000e-3, -2.7396629000e-3 /)
  LoveDat(1:4,884) = (/ 883.0, -5.7364793000, 1.7926742000e-3, -2.7374601000e-3 /)
  LoveDat(1:4,885) = (/ 884.0, -5.7379579000, 1.7916039000e-3, -2.7352601000e-3 /)
  LoveDat(1:4,886) = (/ 885.0, -5.7394323000, 1.7905340000e-3, -2.7330627000e-3 /)
  LoveDat(1:4,887) = (/ 886.0, -5.7409027000, 1.7894644000e-3, -2.7308680000e-3 /)
  LoveDat(1:4,888) = (/ 887.0, -5.7423689000, 1.7883951000e-3, -2.7286759000e-3 /)
  LoveDat(1:4,889) = (/ 888.0, -5.7438310000, 1.7873263000e-3, -2.7264864000e-3 /)
  LoveDat(1:4,890) = (/ 889.0, -5.7452891000, 1.7862578000e-3, -2.7242996000e-3 /)
  LoveDat(1:4,891) = (/ 890.0, -5.7467430000, 1.7851896000e-3, -2.7221154000e-3 /)
  LoveDat(1:4,892) = (/ 891.0, -5.7481929000, 1.7841218000e-3, -2.7199338000e-3 /)
  LoveDat(1:4,893) = (/ 892.0, -5.7496387000, 1.7830544000e-3, -2.7177548000e-3 /)
  LoveDat(1:4,894) = (/ 893.0, -5.7510805000, 1.7819873000e-3, -2.7155785000e-3 /)
  LoveDat(1:4,895) = (/ 894.0, -5.7525182000, 1.7809206000e-3, -2.7134047000e-3 /)
  LoveDat(1:4,896) = (/ 895.0, -5.7539519000, 1.7798543000e-3, -2.7112336000e-3 /)
  LoveDat(1:4,897) = (/ 896.0, -5.7553816000, 1.7787883000e-3, -2.7090651000e-3 /)
  LoveDat(1:4,898) = (/ 897.0, -5.7568072000, 1.7777227000e-3, -2.7068991000e-3 /)
  LoveDat(1:4,899) = (/ 898.0, -5.7582289000, 1.7766574000e-3, -2.7047358000e-3 /)
  LoveDat(1:4,900) = (/ 899.0, -5.7596465000, 1.7755925000e-3, -2.7025750000e-3 /)
  LoveDat(1:4,901) = (/ 900.0, -5.7610602000, 1.7745280000e-3, -2.7004169000e-3 /)
  LoveDat(1:4,902) = (/ 901.0, -5.7624698000, 1.7734638000e-3, -2.6982613000e-3 /)
  LoveDat(1:4,903) = (/ 902.0, -5.7638755000, 1.7724000000e-3, -2.6961082000e-3 /)
  LoveDat(1:4,904) = (/ 903.0, -5.7652772000, 1.7713365000e-3, -2.6939578000e-3 /)
  LoveDat(1:4,905) = (/ 904.0, -5.7666750000, 1.7702734000e-3, -2.6918099000e-3 /)
  LoveDat(1:4,906) = (/ 905.0, -5.7680688000, 1.7692107000e-3, -2.6896645000e-3 /)
  LoveDat(1:4,907) = (/ 906.0, -5.7694587000, 1.7681484000e-3, -2.6875218000e-3 /)
  LoveDat(1:4,908) = (/ 907.0, -5.7708447000, 1.7670864000e-3, -2.6853815000e-3 /)
  LoveDat(1:4,909) = (/ 908.0, -5.7722267000, 1.7660247000e-3, -2.6832438000e-3 /)
  LoveDat(1:4,910) = (/ 909.0, -5.7736048000, 1.7649635000e-3, -2.6811087000e-3 /)
  LoveDat(1:4,911) = (/ 910.0, -5.7749791000, 1.7639026000e-3, -2.6789761000e-3 /)
  LoveDat(1:4,912) = (/ 911.0, -5.7763494000, 1.7628421000e-3, -2.6768460000e-3 /)
  LoveDat(1:4,913) = (/ 912.0, -5.7777158000, 1.7617819000e-3, -2.6747185000e-3 /)
  LoveDat(1:4,914) = (/ 913.0, -5.7790784000, 1.7607221000e-3, -2.6725934000e-3 /)
  LoveDat(1:4,915) = (/ 914.0, -5.7804371000, 1.7596627000e-3, -2.6704709000e-3 /)
  LoveDat(1:4,916) = (/ 915.0, -5.7817919000, 1.7586037000e-3, -2.6683510000e-3 /)
  LoveDat(1:4,917) = (/ 916.0, -5.7831429000, 1.7575450000e-3, -2.6662335000e-3 /)
  LoveDat(1:4,918) = (/ 917.0, -5.7844901000, 1.7564867000e-3, -2.6641185000e-3 /)
  LoveDat(1:4,919) = (/ 918.0, -5.7858334000, 1.7554288000e-3, -2.6620060000e-3 /)
  LoveDat(1:4,920) = (/ 919.0, -5.7871729000, 1.7543712000e-3, -2.6598961000e-3 /)
  LoveDat(1:4,921) = (/ 920.0, -5.7885086000, 1.7533140000e-3, -2.6577886000e-3 /)
  LoveDat(1:4,922) = (/ 921.0, -5.7898405000, 1.7522572000e-3, -2.6556836000e-3 /)
  LoveDat(1:4,923) = (/ 922.0, -5.7911686000, 1.7512008000e-3, -2.6535811000e-3 /)
  LoveDat(1:4,924) = (/ 923.0, -5.7924928000, 1.7501447000e-3, -2.6514811000e-3 /)
  LoveDat(1:4,925) = (/ 924.0, -5.7938134000, 1.7490890000e-3, -2.6493836000e-3 /)
  LoveDat(1:4,926) = (/ 925.0, -5.7951301000, 1.7480337000e-3, -2.6472885000e-3 /)
  LoveDat(1:4,927) = (/ 926.0, -5.7964431000, 1.7469788000e-3, -2.6451960000e-3 /)
  LoveDat(1:4,928) = (/ 927.0, -5.7977523000, 1.7459242000e-3, -2.6431058000e-3 /)
  LoveDat(1:4,929) = (/ 928.0, -5.7990578000, 1.7448700000e-3, -2.6410182000e-3 /)
  LoveDat(1:4,930) = (/ 929.0, -5.8003595000, 1.7438162000e-3, -2.6389330000e-3 /)
  LoveDat(1:4,931) = (/ 930.0, -5.8016575000, 1.7427628000e-3, -2.6368502000e-3 /)
  LoveDat(1:4,932) = (/ 931.0, -5.8029518000, 1.7417098000e-3, -2.6347699000e-3 /)
  LoveDat(1:4,933) = (/ 932.0, -5.8042424000, 1.7406571000e-3, -2.6326921000e-3 /)
  LoveDat(1:4,934) = (/ 933.0, -5.8055293000, 1.7396049000e-3, -2.6306167000e-3 /)
  LoveDat(1:4,935) = (/ 934.0, -5.8068125000, 1.7385530000e-3, -2.6285437000e-3 /)
  LoveDat(1:4,936) = (/ 935.0, -5.8080920000, 1.7375015000e-3, -2.6264732000e-3 /)
  LoveDat(1:4,937) = (/ 936.0, -5.8093679000, 1.7364504000e-3, -2.6244051000e-3 /)
  LoveDat(1:4,938) = (/ 937.0, -5.8106400000, 1.7353996000e-3, -2.6223394000e-3 /)
  LoveDat(1:4,939) = (/ 938.0, -5.8119085000, 1.7343493000e-3, -2.6202762000e-3 /)
  LoveDat(1:4,940) = (/ 939.0, -5.8131734000, 1.7332993000e-3, -2.6182153000e-3 /)
  LoveDat(1:4,941) = (/ 940.0, -5.8144346000, 1.7322497000e-3, -2.6161569000e-3 /)
  LoveDat(1:4,942) = (/ 941.0, -5.8156922000, 1.7312006000e-3, -2.6141009000e-3 /)
  LoveDat(1:4,943) = (/ 942.0, -5.8169461000, 1.7301518000e-3, -2.6120473000e-3 /)
  LoveDat(1:4,944) = (/ 943.0, -5.8181965000, 1.7291034000e-3, -2.6099961000e-3 /)
  LoveDat(1:4,945) = (/ 944.0, -5.8194432000, 1.7280553000e-3, -2.6079473000e-3 /)
  LoveDat(1:4,946) = (/ 945.0, -5.8206864000, 1.7270077000e-3, -2.6059010000e-3 /)
  LoveDat(1:4,947) = (/ 946.0, -5.8219259000, 1.7259605000e-3, -2.6038570000e-3 /)
  LoveDat(1:4,948) = (/ 947.0, -5.8231619000, 1.7249137000e-3, -2.6018154000e-3 /)
  LoveDat(1:4,949) = (/ 948.0, -5.8243943000, 1.7238672000e-3, -2.5997761000e-3 /)
  LoveDat(1:4,950) = (/ 949.0, -5.8256231000, 1.7228212000e-3, -2.5977393000e-3 /)
  LoveDat(1:4,951) = (/ 950.0, -5.8268484000, 1.7217755000e-3, -2.5957048000e-3 /)
  LoveDat(1:4,952) = (/ 951.0, -5.8280701000, 1.7207303000e-3, -2.5936728000e-3 /)
  LoveDat(1:4,953) = (/ 952.0, -5.8292883000, 1.7196854000e-3, -2.5916430000e-3 /)
  LoveDat(1:4,954) = (/ 953.0, -5.8305029000, 1.7186410000e-3, -2.5896157000e-3 /)
  LoveDat(1:4,955) = (/ 954.0, -5.8317141000, 1.7175970000e-3, -2.5875907000e-3 /)
  LoveDat(1:4,956) = (/ 955.0, -5.8329217000, 1.7165533000e-3, -2.5855681000e-3 /)
  LoveDat(1:4,957) = (/ 956.0, -5.8341258000, 1.7155101000e-3, -2.5835478000e-3 /)
  LoveDat(1:4,958) = (/ 957.0, -5.8353264000, 1.7144672000e-3, -2.5815299000e-3 /)
  LoveDat(1:4,959) = (/ 958.0, -5.8365235000, 1.7134248000e-3, -2.5795144000e-3 /)
  LoveDat(1:4,960) = (/ 959.0, -5.8377172000, 1.7123828000e-3, -2.5775012000e-3 /)
  LoveDat(1:4,961) = (/ 960.0, -5.8389074000, 1.7113411000e-3, -2.5754903000e-3 /)
  LoveDat(1:4,962) = (/ 961.0, -5.8400941000, 1.7102999000e-3, -2.5734818000e-3 /)
  LoveDat(1:4,963) = (/ 962.0, -5.8412773000, 1.7092591000e-3, -2.5714756000e-3 /)
  LoveDat(1:4,964) = (/ 963.0, -5.8424571000, 1.7082187000e-3, -2.5694717000e-3 /)
  LoveDat(1:4,965) = (/ 964.0, -5.8436335000, 1.7071787000e-3, -2.5674702000e-3 /)
  LoveDat(1:4,966) = (/ 965.0, -5.8448065000, 1.7061391000e-3, -2.5654710000e-3 /)
  LoveDat(1:4,967) = (/ 966.0, -5.8459760000, 1.7051000000e-3, -2.5634741000e-3 /)
  LoveDat(1:4,968) = (/ 967.0, -5.8471421000, 1.7040612000e-3, -2.5614796000e-3 /)
  LoveDat(1:4,969) = (/ 968.0, -5.8483048000, 1.7030229000e-3, -2.5594873000e-3 /)
  LoveDat(1:4,970) = (/ 969.0, -5.8494641000, 1.7019850000e-3, -2.5574974000e-3 /)
  LoveDat(1:4,971) = (/ 970.0, -5.8506200000, 1.7009475000e-3, -2.5555098000e-3 /)
  LoveDat(1:4,972) = (/ 971.0, -5.8517725000, 1.6999104000e-3, -2.5535245000e-3 /)
  LoveDat(1:4,973) = (/ 972.0, -5.8529217000, 1.6988737000e-3, -2.5515415000e-3 /)
  LoveDat(1:4,974) = (/ 973.0, -5.8540675000, 1.6978374000e-3, -2.5495608000e-3 /)
  LoveDat(1:4,975) = (/ 974.0, -5.8552099000, 1.6968016000e-3, -2.5475824000e-3 /)
  LoveDat(1:4,976) = (/ 975.0, -5.8563490000, 1.6957662000e-3, -2.5456062000e-3 /)
  LoveDat(1:4,977) = (/ 976.0, -5.8574847000, 1.6947312000e-3, -2.5436324000e-3 /)
  LoveDat(1:4,978) = (/ 977.0, -5.8586172000, 1.6936966000e-3, -2.5416609000e-3 /)
  LoveDat(1:4,979) = (/ 978.0, -5.8597463000, 1.6926625000e-3, -2.5396916000e-3 /)
  LoveDat(1:4,980) = (/ 979.0, -5.8608720000, 1.6916288000e-3, -2.5377246000e-3 /)
  LoveDat(1:4,981) = (/ 980.0, -5.8619945000, 1.6905955000e-3, -2.5357599000e-3 /)
  LoveDat(1:4,982) = (/ 981.0, -5.8631137000, 1.6895626000e-3, -2.5337975000e-3 /)
  LoveDat(1:4,983) = (/ 982.0, -5.8642296000, 1.6885302000e-3, -2.5318373000e-3 /)
  LoveDat(1:4,984) = (/ 983.0, -5.8653422000, 1.6874982000e-3, -2.5298794000e-3 /)
  LoveDat(1:4,985) = (/ 984.0, -5.8664515000, 1.6864666000e-3, -2.5279238000e-3 /)
  LoveDat(1:4,986) = (/ 985.0, -5.8675576000, 1.6854355000e-3, -2.5259704000e-3 /)
  LoveDat(1:4,987) = (/ 986.0, -5.8686604000, 1.6844048000e-3, -2.5240193000e-3 /)
  LoveDat(1:4,988) = (/ 987.0, -5.8697599000, 1.6833745000e-3, -2.5220704000e-3 /)
  LoveDat(1:4,989) = (/ 988.0, -5.8708562000, 1.6823447000e-3, -2.5201238000e-3 /)
  LoveDat(1:4,990) = (/ 989.0, -5.8719493000, 1.6813153000e-3, -2.5181795000e-3 /)
  LoveDat(1:4,991) = (/ 990.0, -5.8730391000, 1.6802863000e-3, -2.5162374000e-3 /)
  LoveDat(1:4,992) = (/ 991.0, -5.8741258000, 1.6792578000e-3, -2.5142975000e-3 /)
  LoveDat(1:4,993) = (/ 992.0, -5.8752092000, 1.6782297000e-3, -2.5123598000e-3 /)
  LoveDat(1:4,994) = (/ 993.0, -5.8762894000, 1.6772020000e-3, -2.5104244000e-3 /)
  LoveDat(1:4,995) = (/ 994.0, -5.8773664000, 1.6761748000e-3, -2.5084913000e-3 /)
  LoveDat(1:4,996) = (/ 995.0, -5.8784403000, 1.6751480000e-3, -2.5065603000e-3 /)
  LoveDat(1:4,997) = (/ 996.0, -5.8795109000, 1.6741217000e-3, -2.5046316000e-3 /)
  LoveDat(1:4,998) = (/ 997.0, -5.8805784000, 1.6730958000e-3, -2.5027051000e-3 /)
  LoveDat(1:4,999) = (/ 998.0, -5.8816427000, 1.6720704000e-3, -2.5007809000e-3 /)
  LoveDat(1:4,1000) = (/ 999.0, -5.8827039000, 1.6710454000e-3, -2.4988588000e-3 /)
  LoveDat(1:4,1001) = (/ 1000.0, -5.8837619000, 1.6700209000e-3, -2.4969390000e-3 /)
  LoveDat(1:4,1002) = (/ 1001.0, -5.8848168000, 1.6689968000e-3, -2.4950213000e-3 /)
  LoveDat(1:4,1003) = (/ 1002.0, -5.8858686000, 1.6679732000e-3, -2.4931059000e-3 /)
  LoveDat(1:4,1004) = (/ 1003.0, -5.8869172000, 1.6669500000e-3, -2.4911927000e-3 /)
  LoveDat(1:4,1005) = (/ 1004.0, -5.8879627000, 1.6659272000e-3, -2.4892817000e-3 /)
  LoveDat(1:4,1006) = (/ 1005.0, -5.8890051000, 1.6649049000e-3, -2.4873729000e-3 /)
  LoveDat(1:4,1007) = (/ 1006.0, -5.8900444000, 1.6638831000e-3, -2.4854663000e-3 /)
  LoveDat(1:4,1008) = (/ 1007.0, -5.8910807000, 1.6628617000e-3, -2.4835619000e-3 /)
  LoveDat(1:4,1009) = (/ 1008.0, -5.8921138000, 1.6618408000e-3, -2.4816596000e-3 /)
  LoveDat(1:4,1010) = (/ 1009.0, -5.8931439000, 1.6608204000e-3, -2.4797596000e-3 /)
  LoveDat(1:4,1011) = (/ 1010.0, -5.8941708000, 1.6598004000e-3, -2.4778617000e-3 /)
  LoveDat(1:4,1012) = (/ 1011.0, -5.8951948000, 1.6587808000e-3, -2.4759660000e-3 /)
  LoveDat(1:4,1013) = (/ 1012.0, -5.8962156000, 1.6577617000e-3, -2.4740725000e-3 /)
  LoveDat(1:4,1014) = (/ 1013.0, -5.8972335000, 1.6567431000e-3, -2.4721812000e-3 /)
  LoveDat(1:4,1015) = (/ 1014.0, -5.8982483000, 1.6557250000e-3, -2.4702920000e-3 /)
  LoveDat(1:4,1016) = (/ 1015.0, -5.8992600000, 1.6547073000e-3, -2.4684051000e-3 /)
  LoveDat(1:4,1017) = (/ 1016.0, -5.9002688000, 1.6536901000e-3, -2.4665202000e-3 /)
  LoveDat(1:4,1018) = (/ 1017.0, -5.9012745000, 1.6526733000e-3, -2.4646376000e-3 /)
  LoveDat(1:4,1019) = (/ 1018.0, -5.9022772000, 1.6516570000e-3, -2.4627571000e-3 /)
  LoveDat(1:4,1020) = (/ 1019.0, -5.9032769000, 1.6506412000e-3, -2.4608788000e-3 /)
  LoveDat(1:4,1021) = (/ 1020.0, -5.9042737000, 1.6496258000e-3, -2.4590026000e-3 /)
  LoveDat(1:4,1022) = (/ 1021.0, -5.9052674000, 1.6486109000e-3, -2.4571286000e-3 /)
  LoveDat(1:4,1023) = (/ 1022.0, -5.9062582000, 1.6475965000e-3, -2.4552567000e-3 /)
  LoveDat(1:4,1024) = (/ 1023.0, -5.9072460000, 1.6465826000e-3, -2.4533870000e-3 /)
  LoveDat(1:4,1025) = (/ 1024.0, -5.9082308000, 1.6455691000e-3, -2.4515194000e-3 /)
  LoveDat(1:4,1026) = (/ 1025.0, -5.9092127000, 1.6445561000e-3, -2.4496539000e-3 /)
  LoveDat(1:4,1027) = (/ 1026.0, -5.9101917000, 1.6435436000e-3, -2.4477906000e-3 /)
  LoveDat(1:4,1028) = (/ 1027.0, -5.9111677000, 1.6425316000e-3, -2.4459295000e-3 /)
  LoveDat(1:4,1029) = (/ 1028.0, -5.9121408000, 1.6415200000e-3, -2.4440704000e-3 /)
  LoveDat(1:4,1030) = (/ 1029.0, -5.9131109000, 1.6405089000e-3, -2.4422135000e-3 /)
  LoveDat(1:4,1031) = (/ 1030.0, -5.9140782000, 1.6394983000e-3, -2.4403587000e-3 /)
  LoveDat(1:4,1032) = (/ 1031.0, -5.9150425000, 1.6384882000e-3, -2.4385061000e-3 /)
  LoveDat(1:4,1033) = (/ 1032.0, -5.9160039000, 1.6374786000e-3, -2.4366555000e-3 /)
  LoveDat(1:4,1034) = (/ 1033.0, -5.9169625000, 1.6364694000e-3, -2.4348071000e-3 /)
  LoveDat(1:4,1035) = (/ 1034.0, -5.9179181000, 1.6354607000e-3, -2.4329608000e-3 /)
  LoveDat(1:4,1036) = (/ 1035.0, -5.9188709000, 1.6344526000e-3, -2.4311166000e-3 /)
  LoveDat(1:4,1037) = (/ 1036.0, -5.9198208000, 1.6334449000e-3, -2.4292746000e-3 /)
  LoveDat(1:4,1038) = (/ 1037.0, -5.9207678000, 1.6324377000e-3, -2.4274346000e-3 /)
  LoveDat(1:4,1039) = (/ 1038.0, -5.9217120000, 1.6314309000e-3, -2.4255967000e-3 /)
  LoveDat(1:4,1040) = (/ 1039.0, -5.9226533000, 1.6304247000e-3, -2.4237610000e-3 /)
  LoveDat(1:4,1041) = (/ 1040.0, -5.9235918000, 1.6294190000e-3, -2.4219273000e-3 /)
  LoveDat(1:4,1042) = (/ 1041.0, -5.9245275000, 1.6284137000e-3, -2.4200958000e-3 /)
  LoveDat(1:4,1043) = (/ 1042.0, -5.9254603000, 1.6274090000e-3, -2.4182663000e-3 /)
  LoveDat(1:4,1044) = (/ 1043.0, -5.9263904000, 1.6264047000e-3, -2.4164389000e-3 /)
  LoveDat(1:4,1045) = (/ 1044.0, -5.9273176000, 1.6254009000e-3, -2.4146136000e-3 /)
  LoveDat(1:4,1046) = (/ 1045.0, -5.9282420000, 1.6243977000e-3, -2.4127904000e-3 /)
  LoveDat(1:4,1047) = (/ 1046.0, -5.9291636000, 1.6233949000e-3, -2.4109693000e-3 /)
  LoveDat(1:4,1048) = (/ 1047.0, -5.9300824000, 1.6223926000e-3, -2.4091503000e-3 /)
  LoveDat(1:4,1049) = (/ 1048.0, -5.9309984000, 1.6213909000e-3, -2.4073333000e-3 /)
  LoveDat(1:4,1050) = (/ 1049.0, -5.9319117000, 1.6203896000e-3, -2.4055185000e-3 /)
  LoveDat(1:4,1051) = (/ 1050.0, -5.9328222000, 1.6193888000e-3, -2.4037057000e-3 /)
  LoveDat(1:4,1052) = (/ 1051.0, -5.9337299000, 1.6183886000e-3, -2.4018949000e-3 /)
  LoveDat(1:4,1053) = (/ 1052.0, -5.9346349000, 1.6173888000e-3, -2.4000862000e-3 /)
  LoveDat(1:4,1054) = (/ 1053.0, -5.9355371000, 1.6163895000e-3, -2.3982796000e-3 /)
  LoveDat(1:4,1055) = (/ 1054.0, -5.9364366000, 1.6153908000e-3, -2.3964751000e-3 /)
  LoveDat(1:4,1056) = (/ 1055.0, -5.9373334000, 1.6143925000e-3, -2.3946726000e-3 /)
  LoveDat(1:4,1057) = (/ 1056.0, -5.9382274000, 1.6133948000e-3, -2.3928722000e-3 /)
  LoveDat(1:4,1058) = (/ 1057.0, -5.9391187000, 1.6123976000e-3, -2.3910738000e-3 /)
  LoveDat(1:4,1059) = (/ 1058.0, -5.9400074000, 1.6114009000e-3, -2.3892775000e-3 /)
  LoveDat(1:4,1060) = (/ 1059.0, -5.9408933000, 1.6104046000e-3, -2.3874833000e-3 /)
  LoveDat(1:4,1061) = (/ 1060.0, -5.9417765000, 1.6094089000e-3, -2.3856910000e-3 /)
  LoveDat(1:4,1062) = (/ 1061.0, -5.9426570000, 1.6084138000e-3, -2.3839009000e-3 /)
  LoveDat(1:4,1063) = (/ 1062.0, -5.9435349000, 1.6074191000e-3, -2.3821127000e-3 /)
  LoveDat(1:4,1064) = (/ 1063.0, -5.9444100000, 1.6064249000e-3, -2.3803266000e-3 /)
  LoveDat(1:4,1065) = (/ 1064.0, -5.9452825000, 1.6054313000e-3, -2.3785426000e-3 /)
  LoveDat(1:4,1066) = (/ 1065.0, -5.9461524000, 1.6044382000e-3, -2.3767606000e-3 /)
  LoveDat(1:4,1067) = (/ 1066.0, -5.9470196000, 1.6034456000e-3, -2.3749806000e-3 /)
  LoveDat(1:4,1068) = (/ 1067.0, -5.9478841000, 1.6024535000e-3, -2.3732026000e-3 /)
  LoveDat(1:4,1069) = (/ 1068.0, -5.9487460000, 1.6014619000e-3, -2.3714267000e-3 /)
  LoveDat(1:4,1070) = (/ 1069.0, -5.9496053000, 1.6004709000e-3, -2.3696528000e-3 /)
  LoveDat(1:4,1071) = (/ 1070.0, -5.9504619000, 1.5994804000e-3, -2.3678809000e-3 /)
  LoveDat(1:4,1072) = (/ 1071.0, -5.9513160000, 1.5984904000e-3, -2.3661110000e-3 /)
  LoveDat(1:4,1073) = (/ 1072.0, -5.9521674000, 1.5975009000e-3, -2.3643432000e-3 /)
  LoveDat(1:4,1074) = (/ 1073.0, -5.9530162000, 1.5965120000e-3, -2.3625773000e-3 /)
  LoveDat(1:4,1075) = (/ 1074.0, -5.9538624000, 1.5955236000e-3, -2.3608135000e-3 /)
  LoveDat(1:4,1076) = (/ 1075.0, -5.9547061000, 1.5945357000e-3, -2.3590517000e-3 /)
  LoveDat(1:4,1077) = (/ 1076.0, -5.9555471000, 1.5935483000e-3, -2.3572919000e-3 /)
  LoveDat(1:4,1078) = (/ 1077.0, -5.9563856000, 1.5925615000e-3, -2.3555341000e-3 /)
  LoveDat(1:4,1079) = (/ 1078.0, -5.9572215000, 1.5915752000e-3, -2.3537782000e-3 /)
  LoveDat(1:4,1080) = (/ 1079.0, -5.9580548000, 1.5905894000e-3, -2.3520244000e-3 /)
  LoveDat(1:4,1081) = (/ 1080.0, -5.9588856000, 1.5896041000e-3, -2.3502726000e-3 /)
  LoveDat(1:4,1082) = (/ 1081.0, -5.9597138000, 1.5886194000e-3, -2.3485228000e-3 /)
  LoveDat(1:4,1083) = (/ 1082.0, -5.9605395000, 1.5876353000e-3, -2.3467750000e-3 /)
  LoveDat(1:4,1084) = (/ 1083.0, -5.9613627000, 1.5866516000e-3, -2.3450292000e-3 /)
  LoveDat(1:4,1085) = (/ 1084.0, -5.9621833000, 1.5856685000e-3, -2.3432853000e-3 /)
  LoveDat(1:4,1086) = (/ 1085.0, -5.9630014000, 1.5846860000e-3, -2.3415434000e-3 /)
  LoveDat(1:4,1087) = (/ 1086.0, -5.9638170000, 1.5837039000e-3, -2.3398036000e-3 /)
  LoveDat(1:4,1088) = (/ 1087.0, -5.9646301000, 1.5827224000e-3, -2.3380657000e-3 /)
  LoveDat(1:4,1089) = (/ 1088.0, -5.9654407000, 1.5817415000e-3, -2.3363297000e-3 /)
  LoveDat(1:4,1090) = (/ 1089.0, -5.9662488000, 1.5807611000e-3, -2.3345958000e-3 /)
  LoveDat(1:4,1091) = (/ 1090.0, -5.9670544000, 1.5797812000e-3, -2.3328638000e-3 /)
  LoveDat(1:4,1092) = (/ 1091.0, -5.9678575000, 1.5788019000e-3, -2.3311338000e-3 /)
  LoveDat(1:4,1093) = (/ 1092.0, -5.9686582000, 1.5778231000e-3, -2.3294057000e-3 /)
  LoveDat(1:4,1094) = (/ 1093.0, -5.9694563000, 1.5768449000e-3, -2.3276797000e-3 /)
  LoveDat(1:4,1095) = (/ 1094.0, -5.9702521000, 1.5758672000e-3, -2.3259555000e-3 /)
  LoveDat(1:4,1096) = (/ 1095.0, -5.9710453000, 1.5748901000e-3, -2.3242334000e-3 /)
  LoveDat(1:4,1097) = (/ 1096.0, -5.9718361000, 1.5739135000e-3, -2.3225132000e-3 /)
  LoveDat(1:4,1098) = (/ 1097.0, -5.9726245000, 1.5729374000e-3, -2.3207950000e-3 /)
  LoveDat(1:4,1099) = (/ 1098.0, -5.9734105000, 1.5719619000e-3, -2.3190787000e-3 /)
  LoveDat(1:4,1100) = (/ 1099.0, -5.9741940000, 1.5709870000e-3, -2.3173643000e-3 /)
  LoveDat(1:4,1101) = (/ 1100.0, -5.9749751000, 1.5700126000e-3, -2.3156519000e-3 /)
  LoveDat(1:4,1102) = (/ 1101.0, -5.9757538000, 1.5690387000e-3, -2.3139415000e-3 /)
  LoveDat(1:4,1103) = (/ 1102.0, -5.9765300000, 1.5680654000e-3, -2.3122330000e-3 /)
  LoveDat(1:4,1104) = (/ 1103.0, -5.9773039000, 1.5670927000e-3, -2.3105264000e-3 /)
  LoveDat(1:4,1105) = (/ 1104.0, -5.9780754000, 1.5661205000e-3, -2.3088218000e-3 /)
  LoveDat(1:4,1106) = (/ 1105.0, -5.9788445000, 1.5651489000e-3, -2.3071191000e-3 /)
  LoveDat(1:4,1107) = (/ 1106.0, -5.9796112000, 1.5641778000e-3, -2.3054184000e-3 /)
  LoveDat(1:4,1108) = (/ 1107.0, -5.9803755000, 1.5632073000e-3, -2.3037195000e-3 /)
  LoveDat(1:4,1109) = (/ 1108.0, -5.9811375000, 1.5622374000e-3, -2.3020226000e-3 /)
  LoveDat(1:4,1110) = (/ 1109.0, -5.9818971000, 1.5612680000e-3, -2.3003277000e-3 /)
  LoveDat(1:4,1111) = (/ 1110.0, -5.9826543000, 1.5602991000e-3, -2.2986346000e-3 /)
  LoveDat(1:4,1112) = (/ 1111.0, -5.9834092000, 1.5593309000e-3, -2.2969435000e-3 /)
  LoveDat(1:4,1113) = (/ 1112.0, -5.9841618000, 1.5583632000e-3, -2.2952543000e-3 /)
  LoveDat(1:4,1114) = (/ 1113.0, -5.9849120000, 1.5573960000e-3, -2.2935670000e-3 /)
  LoveDat(1:4,1115) = (/ 1114.0, -5.9856599000, 1.5564294000e-3, -2.2918817000e-3 /)
  LoveDat(1:4,1116) = (/ 1115.0, -5.9864054000, 1.5554634000e-3, -2.2901982000e-3 /)
  LoveDat(1:4,1117) = (/ 1116.0, -5.9871487000, 1.5544979000e-3, -2.2885167000e-3 /)
  LoveDat(1:4,1118) = (/ 1117.0, -5.9878896000, 1.5535331000e-3, -2.2868370000e-3 /)
  LoveDat(1:4,1119) = (/ 1118.0, -5.9886282000, 1.5525687000e-3, -2.2851593000e-3 /)
  LoveDat(1:4,1120) = (/ 1119.0, -5.9893646000, 1.5516050000e-3, -2.2834835000e-3 /)
  LoveDat(1:4,1121) = (/ 1120.0, -5.9900986000, 1.5506418000e-3, -2.2818095000e-3 /)
  LoveDat(1:4,1122) = (/ 1121.0, -5.9908304000, 1.5496792000e-3, -2.2801375000e-3 /)
  LoveDat(1:4,1123) = (/ 1122.0, -5.9915599000, 1.5487171000e-3, -2.2784674000e-3 /)
  LoveDat(1:4,1124) = (/ 1123.0, -5.9922871000, 1.5477557000e-3, -2.2767991000e-3 /)
  LoveDat(1:4,1125) = (/ 1124.0, -5.9930120000, 1.5467948000e-3, -2.2751328000e-3 /)
  LoveDat(1:4,1126) = (/ 1125.0, -5.9937347000, 1.5458344000e-3, -2.2734683000e-3 /)
  LoveDat(1:4,1127) = (/ 1126.0, -5.9944551000, 1.5448747000e-3, -2.2718058000e-3 /)
  LoveDat(1:4,1128) = (/ 1127.0, -5.9951733000, 1.5439155000e-3, -2.2701451000e-3 /)
  LoveDat(1:4,1129) = (/ 1128.0, -5.9958892000, 1.5429569000e-3, -2.2684863000e-3 /)
  LoveDat(1:4,1130) = (/ 1129.0, -5.9966029000, 1.5419989000e-3, -2.2668294000e-3 /)
  LoveDat(1:4,1131) = (/ 1130.0, -5.9973144000, 1.5410414000e-3, -2.2651743000e-3 /)
  LoveDat(1:4,1132) = (/ 1131.0, -5.9980236000, 1.5400845000e-3, -2.2635212000e-3 /)
  LoveDat(1:4,1133) = (/ 1132.0, -5.9987307000, 1.5391282000e-3, -2.2618699000e-3 /)
  LoveDat(1:4,1134) = (/ 1133.0, -5.9994355000, 1.5381725000e-3, -2.2602205000e-3 /)
  LoveDat(1:4,1135) = (/ 1134.0, -6.0001381000, 1.5372174000e-3, -2.2585729000e-3 /)
  LoveDat(1:4,1136) = (/ 1135.0, -6.0008385000, 1.5362628000e-3, -2.2569272000e-3 /)
  LoveDat(1:4,1137) = (/ 1136.0, -6.0015367000, 1.5353089000e-3, -2.2552834000e-3 /)
  LoveDat(1:4,1138) = (/ 1137.0, -6.0022328000, 1.5343555000e-3, -2.2536414000e-3 /)
  LoveDat(1:4,1139) = (/ 1138.0, -6.0029266000, 1.5334027000e-3, -2.2520013000e-3 /)
  LoveDat(1:4,1140) = (/ 1139.0, -6.0036183000, 1.5324504000e-3, -2.2503631000e-3 /)
  LoveDat(1:4,1141) = (/ 1140.0, -6.0043078000, 1.5314988000e-3, -2.2487267000e-3 /)
  LoveDat(1:4,1142) = (/ 1141.0, -6.0049952000, 1.5305477000e-3, -2.2470922000e-3 /)
  LoveDat(1:4,1143) = (/ 1142.0, -6.0056804000, 1.5295973000e-3, -2.2454595000e-3 /)
  LoveDat(1:4,1144) = (/ 1143.0, -6.0063635000, 1.5286474000e-3, -2.2438287000e-3 /)
  LoveDat(1:4,1145) = (/ 1144.0, -6.0070444000, 1.5276981000e-3, -2.2421997000e-3 /)
  LoveDat(1:4,1146) = (/ 1145.0, -6.0077231000, 1.5267494000e-3, -2.2405726000e-3 /)
  LoveDat(1:4,1147) = (/ 1146.0, -6.0083998000, 1.5258013000e-3, -2.2389473000e-3 /)
  LoveDat(1:4,1148) = (/ 1147.0, -6.0090743000, 1.5248538000e-3, -2.2373238000e-3 /)
  LoveDat(1:4,1149) = (/ 1148.0, -6.0097467000, 1.5239068000e-3, -2.2357022000e-3 /)
  LoveDat(1:4,1150) = (/ 1149.0, -6.0104170000, 1.5229605000e-3, -2.2340824000e-3 /)
  LoveDat(1:4,1151) = (/ 1150.0, -6.0110851000, 1.5220147000e-3, -2.2324645000e-3 /)
  LoveDat(1:4,1152) = (/ 1151.0, -6.0117512000, 1.5210696000e-3, -2.2308484000e-3 /)
  LoveDat(1:4,1153) = (/ 1152.0, -6.0124152000, 1.5201250000e-3, -2.2292341000e-3 /)
  LoveDat(1:4,1154) = (/ 1153.0, -6.0130771000, 1.5191811000e-3, -2.2276216000e-3 /)
  LoveDat(1:4,1155) = (/ 1154.0, -6.0137369000, 1.5182377000e-3, -2.2260110000e-3 /)
  LoveDat(1:4,1156) = (/ 1155.0, -6.0143946000, 1.5172949000e-3, -2.2244022000e-3 /)
  LoveDat(1:4,1157) = (/ 1156.0, -6.0150502000, 1.5163527000e-3, -2.2227952000e-3 /)
  LoveDat(1:4,1158) = (/ 1157.0, -6.0157038000, 1.5154112000e-3, -2.2211900000e-3 /)
  LoveDat(1:4,1159) = (/ 1158.0, -6.0163553000, 1.5144702000e-3, -2.2195866000e-3 /)
  LoveDat(1:4,1160) = (/ 1159.0, -6.0170048000, 1.5135298000e-3, -2.2179851000e-3 /)
  LoveDat(1:4,1161) = (/ 1160.0, -6.0176522000, 1.5125900000e-3, -2.2163853000e-3 /)
  LoveDat(1:4,1162) = (/ 1161.0, -6.0182976000, 1.5116508000e-3, -2.2147874000e-3 /)
  LoveDat(1:4,1163) = (/ 1162.0, -6.0189409000, 1.5107122000e-3, -2.2131913000e-3 /)
  LoveDat(1:4,1164) = (/ 1163.0, -6.0195822000, 1.5097742000e-3, -2.2115970000e-3 /)
  LoveDat(1:4,1165) = (/ 1164.0, -6.0202215000, 1.5088369000e-3, -2.2100045000e-3 /)
  LoveDat(1:4,1166) = (/ 1165.0, -6.0208588000, 1.5079001000e-3, -2.2084138000e-3 /)
  LoveDat(1:4,1167) = (/ 1166.0, -6.0214940000, 1.5069639000e-3, -2.2068249000e-3 /)
  LoveDat(1:4,1168) = (/ 1167.0, -6.0221273000, 1.5060283000e-3, -2.2052377000e-3 /)
  LoveDat(1:4,1169) = (/ 1168.0, -6.0227585000, 1.5050934000e-3, -2.2036524000e-3 /)
  LoveDat(1:4,1170) = (/ 1169.0, -6.0233877000, 1.5041590000e-3, -2.2020689000e-3 /)
  LoveDat(1:4,1171) = (/ 1170.0, -6.0240150000, 1.5032253000e-3, -2.2004872000e-3 /)
  LoveDat(1:4,1172) = (/ 1171.0, -6.0246402000, 1.5022921000e-3, -2.1989072000e-3 /)
  LoveDat(1:4,1173) = (/ 1172.0, -6.0252635000, 1.5013596000e-3, -2.1973290000e-3 /)
  LoveDat(1:4,1174) = (/ 1173.0, -6.0258848000, 1.5004276000e-3, -2.1957527000e-3 /)
  LoveDat(1:4,1175) = (/ 1174.0, -6.0265042000, 1.4994963000e-3, -2.1941781000e-3 /)
  LoveDat(1:4,1176) = (/ 1175.0, -6.0271215000, 1.4985656000e-3, -2.1926052000e-3 /)
  LoveDat(1:4,1177) = (/ 1176.0, -6.0277369000, 1.4976355000e-3, -2.1910342000e-3 /)
  LoveDat(1:4,1178) = (/ 1177.0, -6.0283504000, 1.4967060000e-3, -2.1894649000e-3 /)
  LoveDat(1:4,1179) = (/ 1178.0, -6.0289619000, 1.4957771000e-3, -2.1878974000e-3 /)
  LoveDat(1:4,1180) = (/ 1179.0, -6.0295715000, 1.4948489000e-3, -2.1863317000e-3 /)
  LoveDat(1:4,1181) = (/ 1180.0, -6.0301791000, 1.4939212000e-3, -2.1847678000e-3 /)
  LoveDat(1:4,1182) = (/ 1181.0, -6.0307848000, 1.4929942000e-3, -2.1832056000e-3 /)
  LoveDat(1:4,1183) = (/ 1182.0, -6.0313886000, 1.4920677000e-3, -2.1816451000e-3 /)
  LoveDat(1:4,1184) = (/ 1183.0, -6.0319905000, 1.4911419000e-3, -2.1800865000e-3 /)
  LoveDat(1:4,1185) = (/ 1184.0, -6.0325904000, 1.4902167000e-3, -2.1785296000e-3 /)
  LoveDat(1:4,1186) = (/ 1185.0, -6.0331885000, 1.4892921000e-3, -2.1769744000e-3 /)
  LoveDat(1:4,1187) = (/ 1186.0, -6.0337846000, 1.4883682000e-3, -2.1754210000e-3 /)
  LoveDat(1:4,1188) = (/ 1187.0, -6.0343789000, 1.4874448000e-3, -2.1738694000e-3 /)
  LoveDat(1:4,1189) = (/ 1188.0, -6.0349712000, 1.4865221000e-3, -2.1723195000e-3 /)
  LoveDat(1:4,1190) = (/ 1189.0, -6.0355617000, 1.4856000000e-3, -2.1707714000e-3 /)
  LoveDat(1:4,1191) = (/ 1190.0, -6.0361503000, 1.4846785000e-3, -2.1692250000e-3 /)
  LoveDat(1:4,1192) = (/ 1191.0, -6.0367370000, 1.4837576000e-3, -2.1676804000e-3 /)
  LoveDat(1:4,1193) = (/ 1192.0, -6.0373218000, 1.4828373000e-3, -2.1661375000e-3 /)
  LoveDat(1:4,1194) = (/ 1193.0, -6.0379048000, 1.4819177000e-3, -2.1645963000e-3 /)
  LoveDat(1:4,1195) = (/ 1194.0, -6.0384859000, 1.4809986000e-3, -2.1630569000e-3 /)
  LoveDat(1:4,1196) = (/ 1195.0, -6.0390651000, 1.4800802000e-3, -2.1615192000e-3 /)
  LoveDat(1:4,1197) = (/ 1196.0, -6.0396426000, 1.4791625000e-3, -2.1599833000e-3 /)
  LoveDat(1:4,1198) = (/ 1197.0, -6.0402181000, 1.4782453000e-3, -2.1584490000e-3 /)
  LoveDat(1:4,1199) = (/ 1198.0, -6.0407919000, 1.4773288000e-3, -2.1569166000e-3 /)
  LoveDat(1:4,1200) = (/ 1199.0, -6.0413638000, 1.4764129000e-3, -2.1553858000e-3 /)
  LoveDat(1:4,1201) = (/ 1200.0, -6.0419338000, 1.4754976000e-3, -2.1538568000e-3 /)
  LoveDat(1:4,1202) = (/ 1201.0, -6.0425021000, 1.4745829000e-3, -2.1523295000e-3 /)
  LoveDat(1:4,1203) = (/ 1202.0, -6.0430685000, 1.4736689000e-3, -2.1508039000e-3 /)
  LoveDat(1:4,1204) = (/ 1203.0, -6.0436331000, 1.4727555000e-3, -2.1492800000e-3 /)
  LoveDat(1:4,1205) = (/ 1204.0, -6.0441959000, 1.4718427000e-3, -2.1477579000e-3 /)
  LoveDat(1:4,1206) = (/ 1205.0, -6.0447570000, 1.4709305000e-3, -2.1462375000e-3 /)
  LoveDat(1:4,1207) = (/ 1206.0, -6.0453162000, 1.4700190000e-3, -2.1447188000e-3 /)
  LoveDat(1:4,1208) = (/ 1207.0, -6.0458736000, 1.4691081000e-3, -2.1432018000e-3 /)
  LoveDat(1:4,1209) = (/ 1208.0, -6.0464292000, 1.4681978000e-3, -2.1416865000e-3 /)
  LoveDat(1:4,1210) = (/ 1209.0, -6.0469831000, 1.4672882000e-3, -2.1401729000e-3 /)
  LoveDat(1:4,1211) = (/ 1210.0, -6.0475352000, 1.4663791000e-3, -2.1386610000e-3 /)
  LoveDat(1:4,1212) = (/ 1211.0, -6.0480855000, 1.4654707000e-3, -2.1371509000e-3 /)
  LoveDat(1:4,1213) = (/ 1212.0, -6.0486341000, 1.4645630000e-3, -2.1356424000e-3 /)
  LoveDat(1:4,1214) = (/ 1213.0, -6.0491809000, 1.4636558000e-3, -2.1341356000e-3 /)
  LoveDat(1:4,1215) = (/ 1214.0, -6.0497259000, 1.4627493000e-3, -2.1326306000e-3 /)
  LoveDat(1:4,1216) = (/ 1215.0, -6.0502692000, 1.4618435000e-3, -2.1311272000e-3 /)
  LoveDat(1:4,1217) = (/ 1216.0, -6.0508107000, 1.4609382000e-3, -2.1296255000e-3 /)
  LoveDat(1:4,1218) = (/ 1217.0, -6.0513505000, 1.4600336000e-3, -2.1281255000e-3 /)
  LoveDat(1:4,1219) = (/ 1218.0, -6.0518886000, 1.4591297000e-3, -2.1266272000e-3 /)
  LoveDat(1:4,1220) = (/ 1219.0, -6.0524250000, 1.4582263000e-3, -2.1251306000e-3 /)
  LoveDat(1:4,1221) = (/ 1220.0, -6.0529596000, 1.4573236000e-3, -2.1236357000e-3 /)
  LoveDat(1:4,1222) = (/ 1221.0, -6.0534925000, 1.4564215000e-3, -2.1221424000e-3 /)
  LoveDat(1:4,1223) = (/ 1222.0, -6.0540237000, 1.4555201000e-3, -2.1206509000e-3 /)
  LoveDat(1:4,1224) = (/ 1223.0, -6.0545531000, 1.4546193000e-3, -2.1191610000e-3 /)
  LoveDat(1:4,1225) = (/ 1224.0, -6.0550809000, 1.4537191000e-3, -2.1176728000e-3 /)
  LoveDat(1:4,1226) = (/ 1225.0, -6.0556070000, 1.4528196000e-3, -2.1161863000e-3 /)
  LoveDat(1:4,1227) = (/ 1226.0, -6.0561314000, 1.4519207000e-3, -2.1147014000e-3 /)
  LoveDat(1:4,1228) = (/ 1227.0, -6.0566541000, 1.4510224000e-3, -2.1132182000e-3 /)
  LoveDat(1:4,1229) = (/ 1228.0, -6.0571751000, 1.4501248000e-3, -2.1117367000e-3 /)
  LoveDat(1:4,1230) = (/ 1229.0, -6.0576944000, 1.4492278000e-3, -2.1102569000e-3 /)
  LoveDat(1:4,1231) = (/ 1230.0, -6.0582120000, 1.4483315000e-3, -2.1087787000e-3 /)
  LoveDat(1:4,1232) = (/ 1231.0, -6.0587280000, 1.4474358000e-3, -2.1073022000e-3 /)
  LoveDat(1:4,1233) = (/ 1232.0, -6.0592424000, 1.4465407000e-3, -2.1058273000e-3 /)
  LoveDat(1:4,1234) = (/ 1233.0, -6.0597550000, 1.4456463000e-3, -2.1043541000e-3 /)
  LoveDat(1:4,1235) = (/ 1234.0, -6.0602660000, 1.4447525000e-3, -2.1028826000e-3 /)
  LoveDat(1:4,1236) = (/ 1235.0, -6.0607754000, 1.4438593000e-3, -2.1014127000e-3 /)
  LoveDat(1:4,1237) = (/ 1236.0, -6.0612831000, 1.4429668000e-3, -2.0999445000e-3 /)
  LoveDat(1:4,1238) = (/ 1237.0, -6.0617892000, 1.4420750000e-3, -2.0984779000e-3 /)
  LoveDat(1:4,1239) = (/ 1238.0, -6.0622936000, 1.4411837000e-3, -2.0970130000e-3 /)
  LoveDat(1:4,1240) = (/ 1239.0, -6.0627964000, 1.4402931000e-3, -2.0955497000e-3 /)
  LoveDat(1:4,1241) = (/ 1240.0, -6.0632976000, 1.4394032000e-3, -2.0940880000e-3 /)
  LoveDat(1:4,1242) = (/ 1241.0, -6.0637972000, 1.4385139000e-3, -2.0926281000e-3 /)
  LoveDat(1:4,1243) = (/ 1242.0, -6.0642951000, 1.4376252000e-3, -2.0911697000e-3 /)
  LoveDat(1:4,1244) = (/ 1243.0, -6.0647914000, 1.4367372000e-3, -2.0897130000e-3 /)
  LoveDat(1:4,1245) = (/ 1244.0, -6.0652862000, 1.4358498000e-3, -2.0882579000e-3 /)
  LoveDat(1:4,1246) = (/ 1245.0, -6.0657793000, 1.4349631000e-3, -2.0868045000e-3 /)
  LoveDat(1:4,1247) = (/ 1246.0, -6.0662708000, 1.4340770000e-3, -2.0853527000e-3 /)
  LoveDat(1:4,1248) = (/ 1247.0, -6.0667608000, 1.4331916000e-3, -2.0839025000e-3 /)
  LoveDat(1:4,1249) = (/ 1248.0, -6.0672491000, 1.4323068000e-3, -2.0824540000e-3 /)
  LoveDat(1:4,1250) = (/ 1249.0, -6.0677359000, 1.4314226000e-3, -2.0810070000e-3 /)
  LoveDat(1:4,1251) = (/ 1250.0, -6.0682211000, 1.4305391000e-3, -2.0795618000e-3 /)
  LoveDat(1:4,1252) = (/ 1251.0, -6.0687047000, 1.4296562000e-3, -2.0781181000e-3 /)
  LoveDat(1:4,1253) = (/ 1252.0, -6.0691867000, 1.4287740000e-3, -2.0766760000e-3 /)
  LoveDat(1:4,1254) = (/ 1253.0, -6.0696672000, 1.4278925000e-3, -2.0752356000e-3 /)
  LoveDat(1:4,1255) = (/ 1254.0, -6.0701462000, 1.4270115000e-3, -2.0737968000e-3 /)
  LoveDat(1:4,1256) = (/ 1255.0, -6.0706235000, 1.4261312000e-3, -2.0723596000e-3 /)
  LoveDat(1:4,1257) = (/ 1256.0, -6.0710993000, 1.4252516000e-3, -2.0709241000e-3 /)
  LoveDat(1:4,1258) = (/ 1257.0, -6.0715736000, 1.4243726000e-3, -2.0694901000e-3 /)
  LoveDat(1:4,1259) = (/ 1258.0, -6.0720464000, 1.4234943000e-3, -2.0680577000e-3 /)
  LoveDat(1:4,1260) = (/ 1259.0, -6.0725175000, 1.4226166000e-3, -2.0666270000e-3 /)
  LoveDat(1:4,1261) = (/ 1260.0, -6.0729872000, 1.4217396000e-3, -2.0651979000e-3 /)
  LoveDat(1:4,1262) = (/ 1261.0, -6.0734554000, 1.4208632000e-3, -2.0637703000e-3 /)
  LoveDat(1:4,1263) = (/ 1262.0, -6.0739220000, 1.4199874000e-3, -2.0623444000e-3 /)
  LoveDat(1:4,1264) = (/ 1263.0, -6.0743871000, 1.4191123000e-3, -2.0609201000e-3 /)
  LoveDat(1:4,1265) = (/ 1264.0, -6.0748507000, 1.4182379000e-3, -2.0594973000e-3 /)
  LoveDat(1:4,1266) = (/ 1265.0, -6.0753127000, 1.4173641000e-3, -2.0580762000e-3 /)
  LoveDat(1:4,1267) = (/ 1266.0, -6.0757733000, 1.4164910000e-3, -2.0566566000e-3 /)
  LoveDat(1:4,1268) = (/ 1267.0, -6.0762324000, 1.4156185000e-3, -2.0552387000e-3 /)
  LoveDat(1:4,1269) = (/ 1268.0, -6.0766899000, 1.4147466000e-3, -2.0538223000e-3 /)
  LoveDat(1:4,1270) = (/ 1269.0, -6.0771460000, 1.4138754000e-3, -2.0524076000e-3 /)
  LoveDat(1:4,1271) = (/ 1270.0, -6.0776006000, 1.4130049000e-3, -2.0509944000e-3 /)
  LoveDat(1:4,1272) = (/ 1271.0, -6.0780537000, 1.4121350000e-3, -2.0495828000e-3 /)
  LoveDat(1:4,1273) = (/ 1272.0, -6.0785054000, 1.4112658000e-3, -2.0481728000e-3 /)
  LoveDat(1:4,1274) = (/ 1273.0, -6.0789555000, 1.4103972000e-3, -2.0467643000e-3 /)
  LoveDat(1:4,1275) = (/ 1274.0, -6.0794042000, 1.4095293000e-3, -2.0453575000e-3 /)
  LoveDat(1:4,1276) = (/ 1275.0, -6.0798515000, 1.4086620000e-3, -2.0439522000e-3 /)
  LoveDat(1:4,1277) = (/ 1276.0, -6.0802972000, 1.4077954000e-3, -2.0425485000e-3 /)
  LoveDat(1:4,1278) = (/ 1277.0, -6.0807415000, 1.4069294000e-3, -2.0411464000e-3 /)
  LoveDat(1:4,1279) = (/ 1278.0, -6.0811844000, 1.4060641000e-3, -2.0397458000e-3 /)
  LoveDat(1:4,1280) = (/ 1279.0, -6.0816258000, 1.4051994000e-3, -2.0383468000e-3 /)
  LoveDat(1:4,1281) = (/ 1280.0, -6.0820658000, 1.4043354000e-3, -2.0369494000e-3 /)
  LoveDat(1:4,1282) = (/ 1281.0, -6.0825043000, 1.4034720000e-3, -2.0355536000e-3 /)
  LoveDat(1:4,1283) = (/ 1282.0, -6.0829414000, 1.4026093000e-3, -2.0341593000e-3 /)
  LoveDat(1:4,1284) = (/ 1283.0, -6.0833771000, 1.4017473000e-3, -2.0327665000e-3 /)
  LoveDat(1:4,1285) = (/ 1284.0, -6.0838114000, 1.4008859000e-3, -2.0313754000e-3 /)
  LoveDat(1:4,1286) = (/ 1285.0, -6.0842442000, 1.4000251000e-3, -2.0299858000e-3 /)
  LoveDat(1:4,1287) = (/ 1286.0, -6.0846756000, 1.3991650000e-3, -2.0285977000e-3 /)
  LoveDat(1:4,1288) = (/ 1287.0, -6.0851056000, 1.3983056000e-3, -2.0272112000e-3 /)
  LoveDat(1:4,1289) = (/ 1288.0, -6.0855342000, 1.3974468000e-3, -2.0258263000e-3 /)
  LoveDat(1:4,1290) = (/ 1289.0, -6.0859614000, 1.3965887000e-3, -2.0244429000e-3 /)
  LoveDat(1:4,1291) = (/ 1290.0, -6.0863871000, 1.3957312000e-3, -2.0230610000e-3 /)
  LoveDat(1:4,1292) = (/ 1291.0, -6.0868115000, 1.3948744000e-3, -2.0216807000e-3 /)
  LoveDat(1:4,1293) = (/ 1292.0, -6.0872345000, 1.3940183000e-3, -2.0203020000e-3 /)
  LoveDat(1:4,1294) = (/ 1293.0, -6.0876561000, 1.3931628000e-3, -2.0189247000e-3 /)
  LoveDat(1:4,1295) = (/ 1294.0, -6.0880764000, 1.3923079000e-3, -2.0175491000e-3 /)
  LoveDat(1:4,1296) = (/ 1295.0, -6.0884952000, 1.3914537000e-3, -2.0161749000e-3 /)
  LoveDat(1:4,1297) = (/ 1296.0, -6.0889127000, 1.3906002000e-3, -2.0148024000e-3 /)
  LoveDat(1:4,1298) = (/ 1297.0, -6.0893288000, 1.3897473000e-3, -2.0134313000e-3 /)
  LoveDat(1:4,1299) = (/ 1298.0, -6.0897435000, 1.3888951000e-3, -2.0120618000e-3 /)
  LoveDat(1:4,1300) = (/ 1299.0, -6.0901569000, 1.3880436000e-3, -2.0106938000e-3 /)
  LoveDat(1:4,1301) = (/ 1300.0, -6.0905689000, 1.3871927000e-3, -2.0093273000e-3 /)
  LoveDat(1:4,1302) = (/ 1301.0, -6.0909796000, 1.3863424000e-3, -2.0079624000e-3 /)
  LoveDat(1:4,1303) = (/ 1302.0, -6.0913889000, 1.3854928000e-3, -2.0065990000e-3 /)
  LoveDat(1:4,1304) = (/ 1303.0, -6.0917969000, 1.3846439000e-3, -2.0052371000e-3 /)
  LoveDat(1:4,1305) = (/ 1304.0, -6.0922035000, 1.3837956000e-3, -2.0038767000e-3 /)
  LoveDat(1:4,1306) = (/ 1305.0, -6.0926088000, 1.3829480000e-3, -2.0025179000e-3 /)
  LoveDat(1:4,1307) = (/ 1306.0, -6.0930127000, 1.3821011000e-3, -2.0011606000e-3 /)
  LoveDat(1:4,1308) = (/ 1307.0, -6.0934154000, 1.3812548000e-3, -1.9998048000e-3 /)
  LoveDat(1:4,1309) = (/ 1308.0, -6.0938167000, 1.3804091000e-3, -1.9984505000e-3 /)
  LoveDat(1:4,1310) = (/ 1309.0, -6.0942166000, 1.3795642000e-3, -1.9970977000e-3 /)
  LoveDat(1:4,1311) = (/ 1310.0, -6.0946153000, 1.3787199000e-3, -1.9957464000e-3 /)
  LoveDat(1:4,1312) = (/ 1311.0, -6.0950127000, 1.3778762000e-3, -1.9943967000e-3 /)
  LoveDat(1:4,1313) = (/ 1312.0, -6.0954087000, 1.3770332000e-3, -1.9930484000e-3 /)
  LoveDat(1:4,1314) = (/ 1313.0, -6.0958034000, 1.3761909000e-3, -1.9917017000e-3 /)
  LoveDat(1:4,1315) = (/ 1314.0, -6.0961969000, 1.3753492000e-3, -1.9903564000e-3 /)
  LoveDat(1:4,1316) = (/ 1315.0, -6.0965890000, 1.3745082000e-3, -1.9890127000e-3 /)
  LoveDat(1:4,1317) = (/ 1316.0, -6.0969798000, 1.3736678000e-3, -1.9876705000e-3 /)
  LoveDat(1:4,1318) = (/ 1317.0, -6.0973694000, 1.3728281000e-3, -1.9863297000e-3 /)
  LoveDat(1:4,1319) = (/ 1318.0, -6.0977577000, 1.3719890000e-3, -1.9849905000e-3 /)
  LoveDat(1:4,1320) = (/ 1319.0, -6.0981446000, 1.3711507000e-3, -1.9836527000e-3 /)
  LoveDat(1:4,1321) = (/ 1320.0, -6.0985303000, 1.3703129000e-3, -1.9823165000e-3 /)
  LoveDat(1:4,1322) = (/ 1321.0, -6.0989148000, 1.3694759000e-3, -1.9809817000e-3 /)
  LoveDat(1:4,1323) = (/ 1322.0, -6.0992979000, 1.3686395000e-3, -1.9796485000e-3 /)
  LoveDat(1:4,1324) = (/ 1323.0, -6.0996798000, 1.3678037000e-3, -1.9783167000e-3 /)
  LoveDat(1:4,1325) = (/ 1324.0, -6.1000605000, 1.3669686000e-3, -1.9769864000e-3 /)
  LoveDat(1:4,1326) = (/ 1325.0, -6.1004399000, 1.3661342000e-3, -1.9756576000e-3 /)
  LoveDat(1:4,1327) = (/ 1326.0, -6.1008180000, 1.3653005000e-3, -1.9743302000e-3 /)
  LoveDat(1:4,1328) = (/ 1327.0, -6.1011948000, 1.3644674000e-3, -1.9730044000e-3 /)
  LoveDat(1:4,1329) = (/ 1328.0, -6.1015705000, 1.3636349000e-3, -1.9716800000e-3 /)
  LoveDat(1:4,1330) = (/ 1329.0, -6.1019449000, 1.3628031000e-3, -1.9703571000e-3 /)
  LoveDat(1:4,1331) = (/ 1330.0, -6.1023180000, 1.3619720000e-3, -1.9690357000e-3 /)
  LoveDat(1:4,1332) = (/ 1331.0, -6.1026899000, 1.3611416000e-3, -1.9677157000e-3 /)
  LoveDat(1:4,1333) = (/ 1332.0, -6.1030606000, 1.3603118000e-3, -1.9663972000e-3 /)
  LoveDat(1:4,1334) = (/ 1333.0, -6.1034300000, 1.3594826000e-3, -1.9650802000e-3 /)
  LoveDat(1:4,1335) = (/ 1334.0, -6.1037983000, 1.3586541000e-3, -1.9637647000e-3 /)
  LoveDat(1:4,1336) = (/ 1335.0, -6.1041653000, 1.3578263000e-3, -1.9624506000e-3 /)
  LoveDat(1:4,1337) = (/ 1336.0, -6.1045311000, 1.3569992000e-3, -1.9611380000e-3 /)
  LoveDat(1:4,1338) = (/ 1337.0, -6.1048956000, 1.3561727000e-3, -1.9598268000e-3 /)
  LoveDat(1:4,1339) = (/ 1338.0, -6.1052590000, 1.3553468000e-3, -1.9585171000e-3 /)
  LoveDat(1:4,1340) = (/ 1339.0, -6.1056212000, 1.3545217000e-3, -1.9572089000e-3 /)
  LoveDat(1:4,1341) = (/ 1340.0, -6.1059821000, 1.3536972000e-3, -1.9559021000e-3 /)
  LoveDat(1:4,1342) = (/ 1341.0, -6.1063419000, 1.3528733000e-3, -1.9545968000e-3 /)
  LoveDat(1:4,1343) = (/ 1342.0, -6.1067005000, 1.3520501000e-3, -1.9532929000e-3 /)
  LoveDat(1:4,1344) = (/ 1343.0, -6.1070578000, 1.3512276000e-3, -1.9519905000e-3 /)
  LoveDat(1:4,1345) = (/ 1344.0, -6.1074140000, 1.3504057000e-3, -1.9506896000e-3 /)
  LoveDat(1:4,1346) = (/ 1345.0, -6.1077690000, 1.3495845000e-3, -1.9493900000e-3 /)
  LoveDat(1:4,1347) = (/ 1346.0, -6.1081229000, 1.3487640000e-3, -1.9480920000e-3 /)
  LoveDat(1:4,1348) = (/ 1347.0, -6.1084755000, 1.3479441000e-3, -1.9467953000e-3 /)
  LoveDat(1:4,1349) = (/ 1348.0, -6.1088270000, 1.3471249000e-3, -1.9455001000e-3 /)
  LoveDat(1:4,1350) = (/ 1349.0, -6.1091773000, 1.3463063000e-3, -1.9442064000e-3 /)
  LoveDat(1:4,1351) = (/ 1350.0, -6.1095265000, 1.3454884000e-3, -1.9429141000e-3 /)
  LoveDat(1:4,1352) = (/ 1351.0, -6.1098744000, 1.3446712000e-3, -1.9416232000e-3 /)
  LoveDat(1:4,1353) = (/ 1352.0, -6.1102213000, 1.3438546000e-3, -1.9403338000e-3 /)
  LoveDat(1:4,1354) = (/ 1353.0, -6.1105669000, 1.3430387000e-3, -1.9390458000e-3 /)
  LoveDat(1:4,1355) = (/ 1354.0, -6.1109115000, 1.3422234000e-3, -1.9377592000e-3 /)
  LoveDat(1:4,1356) = (/ 1355.0, -6.1112548000, 1.3414088000e-3, -1.9364741000e-3 /)
  LoveDat(1:4,1357) = (/ 1356.0, -6.1115971000, 1.3405949000e-3, -1.9351903000e-3 /)
  LoveDat(1:4,1358) = (/ 1357.0, -6.1119382000, 1.3397816000e-3, -1.9339081000e-3 /)
  LoveDat(1:4,1359) = (/ 1358.0, -6.1122781000, 1.3389690000e-3, -1.9326272000e-3 /)
  LoveDat(1:4,1360) = (/ 1359.0, -6.1126170000, 1.3381571000e-3, -1.9313478000e-3 /)
  LoveDat(1:4,1361) = (/ 1360.0, -6.1129546000, 1.3373458000e-3, -1.9300697000e-3 /)
  LoveDat(1:4,1362) = (/ 1361.0, -6.1132912000, 1.3365352000e-3, -1.9287931000e-3 /)
  LoveDat(1:4,1363) = (/ 1362.0, -6.1136267000, 1.3357252000e-3, -1.9275180000e-3 /)
  LoveDat(1:4,1364) = (/ 1363.0, -6.1139610000, 1.3349159000e-3, -1.9262442000e-3 /)
  LoveDat(1:4,1365) = (/ 1364.0, -6.1142942000, 1.3341073000e-3, -1.9249718000e-3 /)
  LoveDat(1:4,1366) = (/ 1365.0, -6.1146263000, 1.3332993000e-3, -1.9237009000e-3 /)
  LoveDat(1:4,1367) = (/ 1366.0, -6.1149573000, 1.3324920000e-3, -1.9224314000e-3 /)
  LoveDat(1:4,1368) = (/ 1367.0, -6.1152872000, 1.3316853000e-3, -1.9211632000e-3 /)
  LoveDat(1:4,1369) = (/ 1368.0, -6.1156160000, 1.3308793000e-3, -1.9198965000e-3 /)
  LoveDat(1:4,1370) = (/ 1369.0, -6.1159437000, 1.3300740000e-3, -1.9186312000e-3 /)
  LoveDat(1:4,1371) = (/ 1370.0, -6.1162702000, 1.3292693000e-3, -1.9173673000e-3 /)
  LoveDat(1:4,1372) = (/ 1371.0, -6.1165957000, 1.3284653000e-3, -1.9161048000e-3 /)
  LoveDat(1:4,1373) = (/ 1372.0, -6.1169202000, 1.3276619000e-3, -1.9148437000e-3 /)
  LoveDat(1:4,1374) = (/ 1373.0, -6.1172435000, 1.3268592000e-3, -1.9135840000e-3 /)
  LoveDat(1:4,1375) = (/ 1374.0, -6.1175657000, 1.3260572000e-3, -1.9123257000e-3 /)
  LoveDat(1:4,1376) = (/ 1375.0, -6.1178869000, 1.3252558000e-3, -1.9110688000e-3 /)
  LoveDat(1:4,1377) = (/ 1376.0, -6.1182070000, 1.3244551000e-3, -1.9098133000e-3 /)
  LoveDat(1:4,1378) = (/ 1377.0, -6.1185260000, 1.3236551000e-3, -1.9085591000e-3 /)
  LoveDat(1:4,1379) = (/ 1378.0, -6.1188440000, 1.3228557000e-3, -1.9073064000e-3 /)
  LoveDat(1:4,1380) = (/ 1379.0, -6.1191609000, 1.3220569000e-3, -1.9060550000e-3 /)
  LoveDat(1:4,1381) = (/ 1380.0, -6.1194767000, 1.3212589000e-3, -1.9048051000e-3 /)
  LoveDat(1:4,1382) = (/ 1381.0, -6.1197915000, 1.3204614000e-3, -1.9035565000e-3 /)
  LoveDat(1:4,1383) = (/ 1382.0, -6.1201052000, 1.3196647000e-3, -1.9023093000e-3 /)
  LoveDat(1:4,1384) = (/ 1383.0, -6.1204179000, 1.3188686000e-3, -1.9010635000e-3 /)
  LoveDat(1:4,1385) = (/ 1384.0, -6.1207295000, 1.3180732000e-3, -1.8998190000e-3 /)
  LoveDat(1:4,1386) = (/ 1385.0, -6.1210401000, 1.3172784000e-3, -1.8985760000e-3 /)
  LoveDat(1:4,1387) = (/ 1386.0, -6.1213496000, 1.3164843000e-3, -1.8973343000e-3 /)
  LoveDat(1:4,1388) = (/ 1387.0, -6.1216581000, 1.3156908000e-3, -1.8960940000e-3 /)
  LoveDat(1:4,1389) = (/ 1388.0, -6.1219656000, 1.3148980000e-3, -1.8948550000e-3 /)
  LoveDat(1:4,1390) = (/ 1389.0, -6.1222720000, 1.3141059000e-3, -1.8936175000e-3 /)
  LoveDat(1:4,1391) = (/ 1390.0, -6.1225774000, 1.3133144000e-3, -1.8923813000e-3 /)
  LoveDat(1:4,1392) = (/ 1391.0, -6.1228818000, 1.3125236000e-3, -1.8911464000e-3 /)
  LoveDat(1:4,1393) = (/ 1392.0, -6.1231851000, 1.3117334000e-3, -1.8899130000e-3 /)
  LoveDat(1:4,1394) = (/ 1393.0, -6.1234875000, 1.3109439000e-3, -1.8886809000e-3 /)
  LoveDat(1:4,1395) = (/ 1394.0, -6.1237888000, 1.3101551000e-3, -1.8874501000e-3 /)
  LoveDat(1:4,1396) = (/ 1395.0, -6.1240891000, 1.3093669000e-3, -1.8862207000e-3 /)
  LoveDat(1:4,1397) = (/ 1396.0, -6.1243884000, 1.3085794000e-3, -1.8849927000e-3 /)
  LoveDat(1:4,1398) = (/ 1397.0, -6.1246867000, 1.3077925000e-3, -1.8837660000e-3 /)
  LoveDat(1:4,1399) = (/ 1398.0, -6.1249840000, 1.3070063000e-3, -1.8825407000e-3 /)
  LoveDat(1:4,1400) = (/ 1399.0, -6.1252803000, 1.3062208000e-3, -1.8813168000e-3 /)
  LoveDat(1:4,1401) = (/ 1400.0, -6.1255756000, 1.3054359000e-3, -1.8800942000e-3 /)
  LoveDat(1:4,1402) = (/ 1401.0, -6.1258699000, 1.3046517000e-3, -1.8788729000e-3 /)
  LoveDat(1:4,1403) = (/ 1402.0, -6.1261632000, 1.3038681000e-3, -1.8776530000e-3 /)
  LoveDat(1:4,1404) = (/ 1403.0, -6.1264555000, 1.3030852000e-3, -1.8764345000e-3 /)
  LoveDat(1:4,1405) = (/ 1404.0, -6.1267469000, 1.3023029000e-3, -1.8752172000e-3 /)
  LoveDat(1:4,1406) = (/ 1405.0, -6.1270372000, 1.3015213000e-3, -1.8740014000e-3 /)
  LoveDat(1:4,1407) = (/ 1406.0, -6.1273266000, 1.3007404000e-3, -1.8727868000e-3 /)
  LoveDat(1:4,1408) = (/ 1407.0, -6.1276150000, 1.2999601000e-3, -1.8715736000e-3 /)
  LoveDat(1:4,1409) = (/ 1408.0, -6.1279024000, 1.2991804000e-3, -1.8703618000e-3 /)
  LoveDat(1:4,1410) = (/ 1409.0, -6.1281889000, 1.2984015000e-3, -1.8691513000e-3 /)
  LoveDat(1:4,1411) = (/ 1410.0, -6.1284744000, 1.2976231000e-3, -1.8679421000e-3 /)
  LoveDat(1:4,1412) = (/ 1411.0, -6.1287590000, 1.2968455000e-3, -1.8667343000e-3 /)
  LoveDat(1:4,1413) = (/ 1412.0, -6.1290425000, 1.2960685000e-3, -1.8655277000e-3 /)
  LoveDat(1:4,1414) = (/ 1413.0, -6.1293252000, 1.2952921000e-3, -1.8643226000e-3 /)
  LoveDat(1:4,1415) = (/ 1414.0, -6.1296068000, 1.2945164000e-3, -1.8631187000e-3 /)
  LoveDat(1:4,1416) = (/ 1415.0, -6.1298876000, 1.2937414000e-3, -1.8619162000e-3 /)
  LoveDat(1:4,1417) = (/ 1416.0, -6.1301673000, 1.2929670000e-3, -1.8607150000e-3 /)
  LoveDat(1:4,1418) = (/ 1417.0, -6.1304462000, 1.2921933000e-3, -1.8595151000e-3 /)
  LoveDat(1:4,1419) = (/ 1418.0, -6.1307241000, 1.2914202000e-3, -1.8583165000e-3 /)
  LoveDat(1:4,1420) = (/ 1419.0, -6.1310010000, 1.2906478000e-3, -1.8571193000e-3 /)
  LoveDat(1:4,1421) = (/ 1420.0, -6.1312770000, 1.2898761000e-3, -1.8559234000e-3 /)
  LoveDat(1:4,1422) = (/ 1421.0, -6.1315521000, 1.2891050000e-3, -1.8547288000e-3 /)
  LoveDat(1:4,1423) = (/ 1422.0, -6.1318263000, 1.2883345000e-3, -1.8535355000e-3 /)
  LoveDat(1:4,1424) = (/ 1423.0, -6.1320995000, 1.2875647000e-3, -1.8523435000e-3 /)
  LoveDat(1:4,1425) = (/ 1424.0, -6.1323718000, 1.2867956000e-3, -1.8511528000e-3 /)
  LoveDat(1:4,1426) = (/ 1425.0, -6.1326432000, 1.2860271000e-3, -1.8499635000e-3 /)
  LoveDat(1:4,1427) = (/ 1426.0, -6.1329136000, 1.2852592000e-3, -1.8487754000e-3 /)
  LoveDat(1:4,1428) = (/ 1427.0, -6.1331832000, 1.2844921000e-3, -1.8475887000e-3 /)
  LoveDat(1:4,1429) = (/ 1428.0, -6.1334518000, 1.2837255000e-3, -1.8464033000e-3 /)
  LoveDat(1:4,1430) = (/ 1429.0, -6.1337196000, 1.2829597000e-3, -1.8452191000e-3 /)
  LoveDat(1:4,1431) = (/ 1430.0, -6.1339864000, 1.2821945000e-3, -1.8440363000e-3 /)
  LoveDat(1:4,1432) = (/ 1431.0, -6.1342523000, 1.2814299000e-3, -1.8428548000e-3 /)
  LoveDat(1:4,1433) = (/ 1432.0, -6.1345173000, 1.2806660000e-3, -1.8416745000e-3 /)
  LoveDat(1:4,1434) = (/ 1433.0, -6.1347815000, 1.2799027000e-3, -1.8404956000e-3 /)
  LoveDat(1:4,1435) = (/ 1434.0, -6.1350447000, 1.2791401000e-3, -1.8393180000e-3 /)
  LoveDat(1:4,1436) = (/ 1435.0, -6.1353070000, 1.2783782000e-3, -1.8381416000e-3 /)
  LoveDat(1:4,1437) = (/ 1436.0, -6.1355685000, 1.2776168000e-3, -1.8369666000e-3 /)
  LoveDat(1:4,1438) = (/ 1437.0, -6.1358290000, 1.2768562000e-3, -1.8357928000e-3 /)
  LoveDat(1:4,1439) = (/ 1438.0, -6.1360887000, 1.2760962000e-3, -1.8346203000e-3 /)
  LoveDat(1:4,1440) = (/ 1439.0, -6.1363475000, 1.2753368000e-3, -1.8334492000e-3 /) !}}}
  LoveDat(1:4,1441) = (/ 1440.0, -6.1366054000, 1.2745781000e-3, -1.8322792000e-3 /)


  ! if (trim(config_ocean_run_mode) .eq. 'init') then
  if (.False.) then

      LoveScaling(:) = 1.0

  else

      allocate(H(lmax+1),L(lmax+1),K(lmax+1))

      H(1:lmax+1) = LoveDat(2,1:lmax+1)
      L(1:lmax+1) = LoveDat(3,1:lmax+1)
      K(1:lmax+1) = LoveDat(4,1:lmax+1)

      ! Convert from CM to CF
      H1 = H(2)
      L1 = L(2)
      K1 = K(2)
      H(2) = 2.0/3.0*(H1 - L1)
      L(2) = -1.0/3.0*(H1 - L1)
      K(2) = -1.0/3.0*H1 - 2.0/3.0*L1 - 1.0

      do m = 0,nlm
          do n = m,nlm
              i = SHOrderDegreeToIndex(n,m,nlm)
              LoveScaling(i) = (1.0 + K(n+1) - H(n+1))/real(2*n+1)
          enddo
      enddo
      LoveScaling = LoveScaling * 3.0 * rhoW / rhoE

  endif

end subroutine!}}}

!> This subroutine finds a named variable in a list of files and reads its
!! values into a domain-decomposed 2-d array
subroutine find_in_files(filenames, varname, array, G, scale)
  character(len=*), dimension(:),   intent(in)  :: filenames !< The names of the files to search for the named variable
  character(len=*),                 intent(in)  :: varname   !< The name of the variable to read
  type(ocean_grid_type),            intent(in)  :: G         !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: array     !< The array to fill with the data
  real,                   optional, intent(in)  :: scale     !< A factor by which to rescale the array.
  ! Local variables
  integer :: nf

  do nf=1,size(filenames)
    if (LEN_TRIM(filenames(nf)) == 0) cycle
    if (field_exists(filenames(nf), varname, MOM_domain=G%Domain)) then
      call MOM_read_data(filenames(nf), varname, array, G%Domain, scale=scale)
      return
    endif
  enddo

  do nf=size(filenames),1,-1
    if (file_exists(filenames(nf), G%Domain)) then
      call MOM_error(FATAL, "MOM_tidal_forcing.F90: Unable to find "// &
         trim(varname)//" in any of the tidal input files, last tried "// &
         trim(filenames(nf)))
    endif
  enddo

  call MOM_error(FATAL, "MOM_tidal_forcing.F90: Unable to find any of the "// &
                  "tidal input files, including "//trim(filenames(1)))

end subroutine find_in_files

!>   This subroutine calculates returns the partial derivative of the local
!! geopotential height with the input sea surface height due to self-attraction
!! and loading.
subroutine tidal_forcing_sensitivity(G, CS, deta_tidal_deta)
  type(ocean_grid_type),  intent(in)  :: G  !< The ocean's grid structure.
  type(tidal_forcing_CS), intent(in)  :: CS !< The control structure returned by a previous call to tidal_forcing_init.
  real,                   intent(out) :: deta_tidal_deta !< The partial derivative of eta_tidal with
                                            !! the local value of eta [nondim].

  if (CS%USE_SAL_SCALAR .and. CS%USE_PREV_TIDES) then
    deta_tidal_deta = 2.0*CS%SAL_SCALAR
  elseif (CS%USE_SAL_SCALAR .or. CS%USE_PREV_TIDES) then
    deta_tidal_deta = CS%SAL_SCALAR
  else
    deta_tidal_deta = 0.0
  endif
end subroutine tidal_forcing_sensitivity

!>   This subroutine calculates the geopotential anomalies that drive the tides,
!! including self-attraction and loading.  Optionally, it also returns the
!! partial derivative of the local geopotential height with the input sea surface
!! height.  For now, eta and eta_tidal are both geopotential heights in depth
!! units, but probably the input for eta should really be replaced with the
!! column mass anomalies.
subroutine calc_tidal_forcing(Time, eta, eta_tidal, G, US, CS)
  type(ocean_grid_type),            intent(in)  :: G         !< The ocean's grid structure.
  type(time_type),                  intent(in)  :: Time      !< The time for the caluculation.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: eta       !< The sea surface height anomaly from
                                                             !! a time-mean geoid [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: eta_tidal !< The tidal forcing geopotential height
                                                             !! anomalies [Z ~> m].
  type(unit_scale_type),            intent(in)  :: US        !< A dimensional unit scaling type
  type(tidal_forcing_CS),           intent(inout)  :: CS        !< The control structure returned by a
                                                             !! previous call to tidal_forcing_init.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))  :: eta_sal      !< SAL
  real :: now       ! The relative time compared with the tidal reference [T ~> s]
  real :: amp_cosomegat, amp_sinomegat ! The tidal amplitudes times the components of phase [Z ~> m]
  real :: cosomegat, sinomegat ! The components of the phase [nondim]
  real :: eta_prop  ! The nondimenional constant of proportionality beteen eta and eta_tidal [nondim]
  integer :: i, j, c, m, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  call cpu_clock_begin(id_clock_tides)

  if (CS%nc == 0) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1 ; eta_tidal(i,j) = 0.0 ; enddo ; enddo
    return
  endif

  now = US%s_to_T * time_type_to_real(Time - cs%time_ref)

  if (CS%USE_SAL_SCALAR .and. CS%USE_PREV_TIDES) then
    eta_prop = 2.0*CS%SAL_SCALAR
  elseif (CS%USE_SAL_SCALAR .or. CS%USE_PREV_TIDES) then
    eta_prop = CS%SAL_SCALAR
  else
    eta_prop = 0.0
  endif

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    eta_tidal(i,j) = eta_prop*eta(i,j)
  enddo ; enddo

  do c=1,CS%nc
    m = CS%struct(c)
    amp_cosomegat = CS%amp(c)*CS%love_no(c) * cos(CS%freq(c)*now + CS%phase0(c))
    amp_sinomegat = CS%amp(c)*CS%love_no(c) * sin(CS%freq(c)*now + CS%phase0(c))
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_tidal(i,j) = eta_tidal(i,j) + (amp_cosomegat*CS%cos_struct(i,j,m) + &
                                         amp_sinomegat*CS%sin_struct(i,j,m))
    enddo ; enddo
  enddo

  if (CS%tidal_sal_from_file) then ; do c=1,CS%nc
    cosomegat = cos(CS%freq(c)*now)
    sinomegat = sin(CS%freq(c)*now)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_tidal(i,j) = eta_tidal(i,j) + CS%ampsal(i,j,c) * &
           (cosomegat*CS%cosphasesal(i,j,c) + sinomegat*CS%sinphasesal(i,j,c))
    enddo ; enddo
  enddo ; endif

  if (CS%USE_PREV_TIDES) then ; do c=1,CS%nc
    cosomegat = cos(CS%freq(c)*now)
    sinomegat = sin(CS%freq(c)*now)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_tidal(i,j) = eta_tidal(i,j) - CS%SAL_SCALAR*CS%amp_prev(i,j,c) * &
          (cosomegat*CS%cosphase_prev(i,j,c) + sinomegat*CS%sinphase_prev(i,j,c))
    enddo ; enddo
  enddo ; endif

  if (CS%TIDAL_SAL_SHT) then
    eta_sal = 0.0
    call calc_tidal_SAL(eta, eta_sal, G, CS%sht)
    call pass_var(eta_sal, G%domain)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_tidal(i,j) = eta_tidal(i,j) + eta_sal(i,j)
    enddo ; enddo
  endif
  call cpu_clock_end(id_clock_tides)

end subroutine calc_tidal_forcing

subroutine calc_tidal_SAL(eta, eta_sal, G, sht)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: eta       !< The sea surface height anomaly from
                                                             !! a time-mean geoid [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: eta_sal   !< The sea surface height anomaly from
                                                             !! a time-mean geoid [Z ~> m].
  type(sht_CS), intent(inout) :: sht
  real, allocatable :: SnmRe_local(:), SnmIm_local(:)
  real, allocatable :: SnmRe_local_reproSum(:), SnmIm_local_reproSum(:)
  real, allocatable :: SnmRe(:), SnmIm(:)
  real, allocatable :: Snm_local(:), Snm_local_reproSum(:), Snm(:)
  real, allocatable :: LoveScaling(:)

  integer :: i, j
  integer :: is, ie, js, je
  integer :: n, m, l
  real :: mFac

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  call cpu_clock_begin(id_clock_SAL)

  allocate(Snm(2*sht%lmax)); Snm = 0.0
  allocate(SnmRe(sht%lmax)); SnmRe = 0.0
  allocate(SnmIm(sht%lmax)); SnmIm = 0.0

  if (sht%bfb) then
    allocate(Snm_local_reproSum(2*sht%lmax)); Snm_local_reproSum = 0.0
    allocate(SnmRe_local_reproSum(sht%lmax)); SnmRe_local_reproSum = 0.0
    allocate(SnmIm_local_reproSum(sht%lmax)); SnmIm_local_reproSum = 0.0
  else
    allocate(Snm_local(2*sht%lmax)); Snm_local = 0.0
    allocate(SnmRe_local(sht%lmax)); SnmRe_local = 0.0
    allocate(SnmIm_local(sht%lmax)); SnmIm_local = 0.0
  endif

  ! Get SAL scaling factors
  ALLOC_(LoveScaling(sht%lmax))
  call getloadLoveNums(sht%nOrder, LoveScaling)

  !!!!!!!!!!!!!!!!!!!!!
  ! Forward Transform
  !!!!!!!!!!!!!!!!!!!!!

  do m = 0, sht%nOrder
    !------------
    ! n = m
    !------------
    n = m
    ! Calculate associated Legendre polynomial for n=m (output pmnm2)
    ! call associatedLegendrePolynomials(n, m, l, sht%pmnm2, sht%pmnm1, sht%pmn)
    call associatedLegendrePolynomials(n, m, l, sht, G)

    ! Compute local integral contribution
    if (sht%bfb) then
      ! do iCell = endIdx, startIdx, -1
      !     SnmRe_local_reproSum(iCell,l) = SnmRe_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmnm2(iCell)*complexFactorRe(iCell,m+1)
      !     SnmIm_local_reproSum(iCell,l) = SnmIm_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmnm2(iCell)*complexFactorIm(iCell,m+1)
      ! enddo
    else
      do j = js,je ; do i = is,ie
        SnmRe_local(l) = SnmRe_local(l) + eta(i,j) * sht%pmnm2(i,j) * sht%complexFactorRe(i,j,m+1)
        SnmIm_local(l) = SnmIm_local(l) + eta(i,j) * sht%pmnm2(i,j) * sht%complexFactorIm(i,j,m+1)
      enddo ; enddo
    endif

    !------------
    ! n = m+1
    !------------
    n = m+1
    if (n <= sht%nOrder) then
      ! Calculate associated Legendre polynomial for n = m+1 using recurrence relationship
      ! call associatedLegendrePolynomials(n, m, l, sht%pmnm2, sht%pmnm1, sht%pmn)
      call associatedLegendrePolynomials(n, m, l, sht, G)

      ! Compute local integral contribution
      if (sht%bfb) then
        ! do iCell = endIdx, startIdx, -1
        !     SnmRe_local_reproSum(iCell,l) = SnmRe_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmnm1(iCell)*complexFactorRe(iCell,m+1)
        !     SnmIm_local_reproSum(iCell,l) = SnmIm_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmnm1(iCell)*complexFactorIm(iCell,m+1)
        ! enddo
      else
        do j = js,je ; do i = is,ie
          SnmRe_local(l) = SnmRe_local(l) + eta(i,j) * sht%pmnm1(i,j) * sht%complexFactorRe(i,j,m+1)
          SnmIm_local(l) = SnmIm_local(l) + eta(i,j) * sht%pmnm1(i,j) * sht%complexFactorIm(i,j,m+1)
        enddo ; enddo
      endif
    endif

    !------------
    ! n > m+1
    !------------
    do n = m+2,sht%nOrder
      ! Calculate associated Legendre polynomial using recurrence relationship
      ! call associatedLegendrePolynomials(n, m, l, sht%pmnm2, sht%pmnm1, sht%pmn)
      call associatedLegendrePolynomials(n, m, l, sht, G)

      ! Update associated Ledgendre polynomial values for next recurrence
      do j = js,je ; do i = is,ie
        sht%pmnm2(i,j) = sht%pmnm1(i,j)
        sht%pmnm1(i,j) = sht%pmn(i,j)
      enddo ; enddo

      ! Compute local integral contribution
      if (sht%bfb) then
        ! do iCell = endIdx, startIdx, -1
        !     SnmRe_local_reproSum(iCell,l) = SnmRe_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmn(iCell)*complexFactorRe(iCell,m+1)
        !     SnmIm_local_reproSum(iCell,l) = SnmIm_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmn(iCell)*complexFactorIm(iCell,m+1)
        ! enddo
      else
        do j = js,je ; do i = is,ie
          SnmRe_local(l) = SnmRe_local(l) + eta(i,j) * sht%pmn(i,j) * sht%complexFactorRe(i,j,m+1)
          SnmIm_local(l) = SnmIm_local(l) + eta(i,j) * sht%pmn(i,j) * sht%complexFactorIm(i,j,m+1)
        enddo ; enddo
      endif
    enddo ! n loop
  enddo ! m loop

  ! call mpas_timer_stop('Parallel SAL: Forward Transform')

  ! call mpas_timer_start('Parallel SAL: Communication')
  if (sht%bfb) then
    ! do m = 1,lmax
    !     do iCell = 1,nCellsOwned 
    !       Snm_local_reproSum(iCell,m) = SnmRe_local_reproSum(iCell,m)
    !       Snm_local_reproSum(iCell,lmax+m) = SnmIm_local_reproSum(iCell,m)
    !     enddo
    ! enddo
  else
    do m = 1,sht%lmax
      Snm_local(m) = SnmRe_local(m)
      Snm_local(sht%lmax+m) = SnmIm_local(m)
    enddo
  endif

  ! Compute global integral by summing local contributions
  if (sht%bfb) then
     !threadNum = mpas_threading_get_thread_num()
     !if ( threadNum == 0 ) then
        !  Snm = mpas_global_sum_nfld(Snm_local_reproSum,dminfo%comm)
     !endif
  else
    call sum_across_PEs(Snm_local, 2*sht%lmax)
  endif

  do m = 1,sht%lmax
    SnmRe(m) = Snm_local(m)
    SnmIm(m) = Snm_local(sht%lmax+m)
  enddo

  ! call mpas_timer_stop('Parallel SAL: Communication')

  ! call mpas_timer_start('Parallel SAL: Inverse Transform')

  !!!!!!!!!!!!!!!!!!!!!
  ! Apply SAL scaling
  !!!!!!!!!!!!!!!!!!!!!

  do m = 0,sht%nOrder
    do n = m,sht%nOrder
      l = SHOrderDegreeToIndex(n,m,sht%norder)
      SnmRe(l) = SnmRe(l)*LoveScaling(l)
      SnmIm(l) = SnmIm(l)*LoveScaling(l)
    enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!
  ! Inverse transform
  !!!!!!!!!!!!!!!!!!!!

  do m = 0,sht%nOrder
    if (m>0) then
      mFac = 2.0
    else
      mFac = 1.0
    endif

    !------------
    ! n = m
    !------------
    n = m
    ! Calculate associated Legendre polynomial using recurrence relationship
    ! call associatedLegendrePolynomials(n, m, l, sht%pmnm2, sht%pmnm1, sht%pmn)
    call associatedLegendrePolynomials(n, m, l, sht, G)

    ! Sum together product of spherical harmonic functions and coefficients
    do j = js,je ; do i = is,ie
      eta_sal(i,j) = eta_sal(i,j) &
        + mFac * sht%pmnm2(i,j) * (SnmRe(l) * sht%complexExpRe(i,j,m+1) + SnmIm(l) * sht%complexExpIm(i,j,m+1))
    enddo ; enddo

    !------------
    ! n = m+1
    !------------
    n = m+1
    if (n <= sht%nOrder) then
      ! Calculate associated Legendre polynomial using recurrence relationship
      ! call associatedLegendrePolynomials(n, m, l, sht%pmnm2, sht%pmnm1, sht%pmn)
      call associatedLegendrePolynomials(n, m, l, sht, G)

      ! Sum together product of spherical harmonic functions and coefficients
      do j = js,je ; do i = is,ie
        eta_sal(i,j) = eta_sal(i,j) &
          + mFac * sht%pmnm1(i,j) * (SnmRe(l) * sht%complexExpRe(i,j,m+1) + SnmIm(l) * sht%complexExpIm(i,j,m+1))
      enddo ; enddo
    endif

    !------------
    ! n > m+1
    !------------
    do n = m+2,sht%nOrder
      ! Calculate associated Legendre polynomial using recurrence relationship
      ! call associatedLegendrePolynomials(n, m, l, sht%pmnm2, sht%pmnm1, sht%pmn)
      call associatedLegendrePolynomials(n, m, l, sht, G)

      ! Update associated Ledgendre polynomial values for next recurrence
      do j = js,je ; do i = is,ie
        sht%pmnm2(i,j) = sht%pmnm1(i,j)
        sht%pmnm1(i,j) = sht%pmn(i,j)
      enddo ; enddo

      ! Sum together product of spherical harmonic functions and coefficients
      do j = js,je ; do i = is,ie
        eta_sal(i,j) = eta_sal(i,j) &
          + mFac * sht%pmn(i,j) * (SnmRe(l) * sht%complexExpRe(i,j,m+1) + SnmIm(l) * sht%complexExpIm(i,j,m+1))
      enddo ; enddo
    enddo ! n loop
  enddo ! m loop

  call cpu_clock_end(id_clock_SAL)
end subroutine calc_tidal_SAL

!> This subroutine deallocates memory associated with the tidal forcing module.
subroutine tidal_forcing_end(CS)
  type(tidal_forcing_CS), intent(inout) :: CS !< The control structure returned by a previous call
                                              !! to tidal_forcing_init; it is deallocated here.

  if (allocated(CS%sin_struct)) deallocate(CS%sin_struct)
  if (allocated(CS%cos_struct)) deallocate(CS%cos_struct)

  if (allocated(CS%cosphasesal)) deallocate(CS%cosphasesal)
  if (allocated(CS%sinphasesal)) deallocate(CS%sinphasesal)
  if (allocated(CS%ampsal))      deallocate(CS%ampsal)

  if (allocated(CS%cosphase_prev)) deallocate(CS%cosphase_prev)
  if (allocated(CS%sinphase_prev)) deallocate(CS%sinphase_prev)
  if (allocated(CS%amp_prev))      deallocate(CS%amp_prev)

  if (CS%tidal_sal_sht) &
    call spherical_harmonics_end(CS%sht)
end subroutine tidal_forcing_end

!> \namespace tidal_forcing
!!
!! Code by Robert Hallberg, August 2005, based on C-code by Harper
!! Simmons, February, 2003, in turn based on code by Brian Arbic.
!!
!!   The main subroutine in this file calculates the total tidal
!! contribution to the geopotential, including self-attraction and
!! loading terms and the astronomical contributions.  All options
!! are selected with entries in a file that is parsed at run-time.
!! Overall tides are enabled with the run-time parameter 'TIDES=True'.
!! Tidal constituents must be individually enabled with lines like
!! 'TIDE_M2=True'.  This file has default values of amplitude,
!! frequency, Love number, and phase at time 0 for the Earth's M2,
!! S2, N2, K2, K1, O1, P1, Q1,  MF, and MM tidal constituents, but
!! the frequency, amplitude and phase ant time 0 for each constituent
!! can be changed at run time by setting variables like TIDE_M2_FREQ,
!! TIDE_M2_AMP and TIDE_M2_PHASE_T0 (for M2).
!!
!!   In addition, the approach to calculating self-attraction and
!! loading is set at run time.  The default is to use the scalar
!! approximation, with a coefficient TIDE_SAL_SCALAR_VALUE that must
!! be set in the run-time file (for global runs, 0.094 is typical).
!! Alternately, TIDAL_SAL_FROM_FILE can be set to read the SAL from
!! a file containing the results of a previous simulation. To iterate
!! the SAL to convergence, USE_PREVIOUS_TIDES may be useful (for
!! details, see Arbic et al., 2004, DSR II). With TIDAL_SAL_FROM_FILE
!! or USE_PREVIOUS_TIDES,a list of input files must be provided to
!! describe each constituent's properties from a previous solution.

end module MOM_tidal_forcing
