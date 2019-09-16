!> Tidal contributions to geopotential
module MOM_tidal_forcing

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
                              CLOCK_MODULE
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : field_exists, file_exists, MOM_read_data
use MOM_time_manager,  only : time_type, time_type_to_real

implicit none ; private

public calc_tidal_forcing, tidal_forcing_init, tidal_forcing_end
public tidal_forcing_sensitivity

#include <MOM_memory.h>

integer, parameter :: MAX_CONSTITUENTS = 10 !< The maximum number of tidal
                                            !! constituents that could be used.

!> The control structure for the MOM_tidal_forcing module
type, public :: tidal_forcing_CS ; private
  logical :: use_sal_scalar !< If true, use the scalar approximation when
                      !! calculating self-attraction and loading.
  logical :: tidal_sal_from_file !< If true, Read the tidal self-attraction
                      !! and loading from input files, specified
                      !! by TIDAL_INPUT_FILE.
  logical :: use_prev_tides !< If true, use the SAL from the previous
                      !! iteration of the tides to facilitate convergence.
  real    :: sal_scalar !< The constant of proportionality between sea surface
                      !! height (really it should be bottom pressure) anomalies
                      !! and bottom geopotential anomalies.
  integer :: nc       !< The number of tidal constituents in use.
  real, dimension(MAX_CONSTITUENTS) :: &
    freq, &           !< The frequency of a tidal constituent [s-1].
    phase0, &         !< The phase of a tidal constituent at time 0, in radians.
    amp, &            !< The amplitude of a tidal constituent at time 0 [m].
    love_no           !< The Love number of a tidal constituent at time 0 [nondim].
  integer :: struct(MAX_CONSTITUENTS) !< An encoded spatial structure for each constituent
  character (len=16) :: const_name(MAX_CONSTITUENTS) !< The name of each constituent

  real, pointer, dimension(:,:,:) :: &
    sin_struct => NULL(), &    !< The sine and cosine based structures that can
    cos_struct => NULL(), &    !< be associated with the astronomical forcing.
    cosphasesal => NULL(), &   !< The cosine and sine of the phase of the
    sinphasesal => NULL(), &   !< self-attraction and loading amphidromes.
    ampsal => NULL(), &        !< The amplitude of the SAL [m].
    cosphase_prev => NULL(), & !< The cosine and sine of the phase of the
    sinphase_prev => NULL(), & !< amphidromes in the previous tidal solutions.
    amp_prev => NULL()         !< The amplitude of the previous tidal solution [m].
end type tidal_forcing_CS

integer :: id_clock_tides !< CPU clock for tides

contains

!> This subroutine allocates space for the static variables used
!! by this module.  The metrics may be effectively 0, 1, or 2-D arrays,
!! while fields like the background viscosities are 2-D arrays.
!! ALLOC is a macro defined in MOM_memory.h for allocate or nothing with
!! static memory.
subroutine tidal_forcing_init(Time, G, param_file, CS)
  type(time_type),       intent(in)    :: Time !< The current model time.
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  type(param_file_type), intent(in)    :: param_file !< A structure to parse for run-time parameters.
  type(tidal_forcing_CS), pointer      :: CS   !< A pointer that is set to point to the control
                                               !! structure for this module.
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G)) :: &
    phase, &          ! The phase of some tidal constituent.
    lat_rad, lon_rad  ! Latitudes and longitudes of h-points in radians.
  real :: deg_to_rad
  real, dimension(MAX_CONSTITUENTS) :: freq_def, phase0_def, amp_def, love_def
  logical :: use_const  ! True if a constituent is being used.
  logical :: use_M2, use_S2, use_N2, use_K2, use_K1, use_O1, use_P1, use_Q1
  logical :: use_MF, use_MM
  logical :: tides      ! True if a tidal forcing is to be used.
  logical :: FAIL_IF_MISSING = .true.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_tidal_forcing" ! This module's name.
  character(len=128) :: mesg
  character(len=200) :: tidal_input_files(4*MAX_CONSTITUENTS)
  integer :: i, j, c, is, ie, js, je, isd, ied, jsd, jed, nc
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd; jed = G%jed

  if (associated(CS)) then
    call MOM_error(WARNING, "tidal_forcing_init called with an associated "// &
                            "control structure.")
    return
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "TIDES", tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)

  if (.not.tides) return

  allocate(CS)

  ! Set up the spatial structure functions for the diurnal, semidiurnal, and
  ! low-frequency tidal components.
  allocate(CS%sin_struct(isd:ied,jsd:jed,3)) ; CS%sin_struct(:,:,:) = 0.0
  allocate(CS%cos_struct(isd:ied,jsd:jed,3)) ; CS%cos_struct(:,:,:) = 0.0
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

  ! Set the parameters for all components that are in use.
  c=0
  if (use_M2) then
    c=c+1 ; CS%const_name(c) = "M2" ; CS%freq(c) = 1.4051890e-4 ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; CS%amp(c) = 0.242334 ; CS%phase0(c) = 0.0
    freq_def(c) = CS%freq(c) ; love_def(c) = CS%love_no(c)
    amp_def(c) = CS%amp(c) ; phase0_def(c) = CS%phase0(c)
  endif

  if (use_S2) then
    c=c+1 ; CS%const_name(c) = "S2" ; CS%freq(c) = 1.4544410e-4 ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; CS%amp(c) = 0.112743 ; CS%phase0(c) = 0.0
    freq_def(c) = CS%freq(c) ; love_def(c) = CS%love_no(c)
    amp_def(c) = CS%amp(c) ; phase0_def(c) = CS%phase0(c)
  endif

  if (use_N2) then
    c=c+1 ; CS%const_name(c) = "N2" ; CS%freq(c) = 1.3787970e-4 ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; CS%amp(c) = 0.046397 ; CS%phase0(c) = 0.0
    freq_def(c) = CS%freq(c) ; love_def(c) = CS%love_no(c)
    amp_def(c) = CS%amp(c) ; phase0_def(c) = CS%phase0(c)
  endif

  if (use_K2) then
    c=c+1 ; CS%const_name(c) = "K2" ; CS%freq(c) = 1.4584234e-4 ; CS%struct(c) = 2
    CS%love_no(c) = 0.693 ; CS%amp(c) = 0.030684 ; CS%phase0(c) = 0.0
    freq_def(c) = CS%freq(c) ; love_def(c) = CS%love_no(c)
    amp_def(c) = CS%amp(c) ; phase0_def(c) = CS%phase0(c)
  endif

  if (use_K1) then
    c=c+1 ; CS%const_name(c) = "K1" ; CS%freq(c) = 0.7292117e-4 ; CS%struct(c) = 1
    CS%love_no(c) = 0.736 ; CS%amp(c) = 0.141565 ; CS%phase0(c) = 0.0
    freq_def(c) = CS%freq(c) ; love_def(c) = CS%love_no(c)
    amp_def(c) = CS%amp(c) ; phase0_def(c) = CS%phase0(c)
  endif

  if (use_O1) then
    c=c+1 ; CS%const_name(c) = "O1" ; CS%freq(c) = 0.6759774e-4 ; CS%struct(c) = 1
    CS%love_no(c) = 0.695 ; CS%amp(c) = 0.100661 ; CS%phase0(c) = 0.0
    freq_def(c) = CS%freq(c) ; love_def(c) = CS%love_no(c)
    amp_def(c) = CS%amp(c) ; phase0_def(c) = CS%phase0(c)
  endif

  if (use_P1) then
    c=c+1 ; CS%const_name(c) = "P1" ; CS%freq(c) = 0.7252295e-4 ; CS%struct(c) = 1
    CS%love_no(c) = 0.706 ; CS%amp(c) = 0.046848 ; CS%phase0(c) = 0.0
    freq_def(c) = CS%freq(c) ; love_def(c) = CS%love_no(c)
    amp_def(c) = CS%amp(c) ; phase0_def(c) = CS%phase0(c)
  endif

  if (use_Q1) then
    c=c+1 ; CS%const_name(c) = "Q1" ; CS%freq(c) = 0.6495854e-4 ; CS%struct(c) = 1
    CS%love_no(c) = 0.695 ; CS%amp(c) = 0.019273 ; CS%phase0(c) = 0.0
    freq_def(c) = CS%freq(c) ; love_def(c) = CS%love_no(c)
    amp_def(c) = CS%amp(c) ; phase0_def(c) = CS%phase0(c)
  endif

  if (use_MF) then
    c=c+1 ; CS%const_name(c) = "MF" ; CS%freq(c) = 0.053234e-4 ; CS%struct(c) = 3
    CS%love_no(c) = 0.693 ; CS%amp(c) = 0.042041 ; CS%phase0(c) = 0.0
    freq_def(c) = CS%freq(c) ; love_def(c) = CS%love_no(c)
    amp_def(c) = CS%amp(c) ; phase0_def(c) = CS%phase0(c)
  endif

  if (use_MM) then
    c=c+1 ; CS%const_name(c) = "MM" ; CS%freq(c) = 0.026392e-4 ; CS%struct(c) = 3
    CS%love_no(c) = 0.693 ; CS%amp(c) = 0.022191 ; CS%phase0(c) = 0.0
    freq_def(c) = CS%freq(c) ; love_def(c) = CS%love_no(c)
    amp_def(c) = CS%amp(c) ; phase0_def(c) = CS%phase0(c)
  endif

  !   Parse the input file to potentially override the default values for the
  ! frequency, amplitude and initial phase of each constituent, and log the
  ! values that are actually used.
  do c=1,nc
    call get_param(param_file, mdl, "TIDE_"//trim(CS%const_name(c))//"_FREQ", CS%freq(c), &
                   "Frequency of the "//trim(CS%const_name(c))//" tidal constituent. "//&
                   "This is only used if TIDES and TIDE_"//trim(CS%const_name(c))// &
                   " are true.", units="s-1", default=freq_def(c))
    call get_param(param_file, mdl, "TIDE_"//trim(CS%const_name(c))//"_AMP", CS%amp(c), &
                   "Amplitude of the "//trim(CS%const_name(c))//" tidal constituent. "//&
                   "This is only used if TIDES and TIDE_"//trim(CS%const_name(c))// &
                   " are true.", units="m", default=amp_def(c))
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
      call find_in_files(tidal_input_files,"PHASE_SAL_"//trim(CS%const_name(c)),phase,G)
      call find_in_files(tidal_input_files,"AMP_SAL_"//trim(CS%const_name(c)),CS%ampsal(:,:,c),G)
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
      call find_in_files(tidal_input_files,"PHASE_PREV_"//trim(CS%const_name(c)),phase,G)
      call find_in_files(tidal_input_files,"AMP_PREV_"//trim(CS%const_name(c)),CS%amp_prev(:,:,c),G)
      call pass_var(phase,             G%domain,complete=.false.)
      call pass_var(CS%amp_prev(:,:,c),G%domain,complete=.true.)
      do j=js-1,je+1 ; do i=is-1,ie+1
        CS%cosphase_prev(i,j,c) = cos(phase(i,j)*deg_to_rad)
        CS%sinphase_prev(i,j,c) = sin(phase(i,j)*deg_to_rad)
      enddo ; enddo
    enddo
  endif

  id_clock_tides = cpu_clock_id('(Ocean tides)', grain=CLOCK_MODULE)

end subroutine tidal_forcing_init

!> This subroutine finds a named variable in a list of files and reads its
!! values into a domain-decomposed 2-d array
subroutine find_in_files(filenames, varname, array, G)
  character(len=*), dimension(:),   intent(in)  :: filenames !< The names of the files to search for the named variable
  character(len=*),                 intent(in)  :: varname   !< The name of the variable to read
  type(ocean_grid_type),            intent(in)  :: G         !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: array     !< The array to fill with the data
  ! Local variables
  integer :: nf

  do nf=1,size(filenames)
    if (LEN_TRIM(filenames(nf)) == 0) cycle
    if (field_exists(filenames(nf), varname, G%Domain%mpp_domain)) then
      call MOM_read_data(filenames(nf), varname, array, G%Domain)
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
  type(tidal_forcing_CS), pointer     :: CS !< The control structure returned by a previous call to tidal_forcing_init.
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
subroutine calc_tidal_forcing(Time, eta, eta_tidal, G, CS, deta_tidal_deta, m_to_Z)
  type(ocean_grid_type),            intent(in)  :: G         !< The ocean's grid structure.
  type(time_type),                  intent(in)  :: Time      !< The time for the caluculation.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: eta       !< The sea surface height anomaly from
                                                             !! a time-mean geoid [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: eta_tidal !< The tidal forcing geopotential height
                                                             !! anomalies [Z ~> m].
  type(tidal_forcing_CS),           pointer     :: CS        !< The control structure returned by a
                                                             !! previous call to tidal_forcing_init.
  real, optional,                   intent(out) :: deta_tidal_deta !< The partial derivative of
                                                             !! eta_tidal with the local value of
                                                             !! eta [nondim].
  real, optional,                   intent(in)  :: m_to_Z    !< A scaling factor from m to the units of eta.

  ! Local variables
  real :: eta_astro(SZI_(G),SZJ_(G))
  real :: eta_SAL(SZI_(G),SZJ_(G))
  real :: now       ! The relative time in seconds.
  real :: amp_cosomegat, amp_sinomegat
  real :: cosomegat, sinomegat
  real :: m_Z       ! A scaling factor from m to depth units.
  real :: eta_prop  ! The nondimenional constant of proportionality beteen eta and eta_tidal.
  integer :: i, j, c, m, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.associated(CS)) return

  call cpu_clock_begin(id_clock_tides)

  if (CS%nc == 0) then
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1 ; eta_tidal(i,j) = 0.0 ; enddo ; enddo
    return
  endif

  now = time_type_to_real(Time)

  if (CS%USE_SAL_SCALAR .and. CS%USE_PREV_TIDES) then
    eta_prop = 2.0*CS%SAL_SCALAR
  elseif (CS%USE_SAL_SCALAR .or. CS%USE_PREV_TIDES) then
    eta_prop = CS%SAL_SCALAR
  else
    eta_prop = 0.0
  endif

  if (present(deta_tidal_deta)) then
    deta_tidal_deta = eta_prop
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1 ; eta_tidal(i,j) = 0.0 ; enddo ; enddo
  else
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_tidal(i,j) = eta_prop*eta(i,j)
    enddo ; enddo
  endif

  m_Z = 1.0 ; if (present(m_to_Z)) m_Z = m_to_Z

  do c=1,CS%nc
    m = CS%struct(c)
    amp_cosomegat = m_Z*CS%amp(c)*CS%love_no(c) * cos(CS%freq(c)*now + CS%phase0(c))
    amp_sinomegat = m_Z*CS%amp(c)*CS%love_no(c) * sin(CS%freq(c)*now + CS%phase0(c))
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_tidal(i,j) = eta_tidal(i,j) + (amp_cosomegat*CS%cos_struct(i,j,m) + &
                                         amp_sinomegat*CS%sin_struct(i,j,m))
    enddo ; enddo
  enddo

  if (CS%tidal_sal_from_file) then ; do c=1,CS%nc
    cosomegat = cos(CS%freq(c)*now)
    sinomegat = sin(CS%freq(c)*now)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_tidal(i,j) = eta_tidal(i,j) + m_Z*CS%ampsal(i,j,c) * &
           (cosomegat*CS%cosphasesal(i,j,c) + sinomegat*CS%sinphasesal(i,j,c))
    enddo ; enddo
  enddo ; endif

  if (CS%USE_PREV_TIDES) then ; do c=1,CS%nc
    cosomegat = cos(CS%freq(c)*now)
    sinomegat = sin(CS%freq(c)*now)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      eta_tidal(i,j) = eta_tidal(i,j) - m_Z*CS%SAL_SCALAR*CS%amp_prev(i,j,c) * &
          (cosomegat*CS%cosphase_prev(i,j,c) + sinomegat*CS%sinphase_prev(i,j,c))
    enddo ; enddo
  enddo ; endif

  call cpu_clock_end(id_clock_tides)

end subroutine calc_tidal_forcing

!> This subroutine deallocates memory associated with the tidal forcing module.
subroutine tidal_forcing_end(CS)
  type(tidal_forcing_CS), pointer :: CS !< The control structure returned by a previous call
                                        !! to tidal_forcing_init; it is deallocated here.

  if (associated(CS%sin_struct)) deallocate(CS%sin_struct)
  if (associated(CS%cos_struct)) deallocate(CS%cos_struct)

  if (associated(CS%cosphasesal)) deallocate(CS%cosphasesal)
  if (associated(CS%sinphasesal)) deallocate(CS%sinphasesal)
  if (associated(CS%ampsal))      deallocate(CS%ampsal)

  if (associated(CS%cosphase_prev)) deallocate(CS%cosphase_prev)
  if (associated(CS%sinphase_prev)) deallocate(CS%sinphase_prev)
  if (associated(CS%amp_prev))      deallocate(CS%amp_prev)

  if (associated(CS)) deallocate(CS)

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
