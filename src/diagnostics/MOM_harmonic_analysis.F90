! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Inline harmonic analysis (conventional)
module MOM_harmonic_analysis

use MOM_time_manager,  only : time_type, real_to_time, time_to_real, time_minus_signed
use MOM_time_manager,  only : set_date, get_date, increment_date
use MOM_time_manager,  only : operator(+), operator(-), operator(<), operator(>), operator(>=)
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_file_parser,   only : param_file_type, get_param
use MOM_io,            only : file_exists, open_ASCII_file, READONLY_FILE, close_file
use MOM_io,            only : MOM_infra_file, vardesc, MOM_field
use MOM_io,            only : var_desc, create_MOM_file, SINGLE_FILE, MOM_write_field
use MOM_error_handler, only : MOM_mesg, MOM_error, NOTE
use MOM_tidal_forcing, only : astro_longitudes, astro_longitudes_init, eq_phase, nodal_fu, tidal_frequency

implicit none ; private

public HA_init, HA_accum

#include <MOM_memory.h>

!> The private control structure for storing the HA info of a particular field
type, private :: HA_type
  character(len=16) :: key = "none"          !< Name of the field of which harmonic analysis is to be performed
  character(len=1)  :: grid                  !< The grid on which the field is defined ('h', 'q', 'u', or 'v')
  real :: old_time = -1.0                    !< The time of the previous accumulating step [T ~> s]
  real, allocatable :: ref(:,:)              !< The initial field in arbitrary units [A]
  real, allocatable :: FtF(:,:)              !< Accumulator of (F' * F) [nondim]
  real, allocatable :: FtSSH(:,:,:)          !< Accumulator of (F' * SSH_in) in arbitrary units [A]
  !>@{ Lower and upper bounds of input data
  integer :: is, ie, js, je
  !>@}
end type HA_type

!> A linked list of control structures that store the HA info of different fields
type, private :: HA_node
  type(HA_type)          :: this             !< Control structure of the current field in the list
  type(HA_node), pointer :: next             !< The list of other fields
end type HA_node

!> The public control structure of the MOM_harmonic_analysis module
type, public :: harmonic_analysis_CS ; private
  logical :: HAready = .false.               !< If true, perform harmonic analysis
  type(time_type) :: &
    time_start, &                            !< Start time of harmonic analysis
    time_end, &                              !< End time of harmonic analysis
    time_ref                                 !< Reference time (t = 0) used to calculate tidal forcing
  real, allocatable, dimension(:) :: &
    freq, &                                  !< The frequency of a tidal constituent [T-1 ~> s-1]
    phase0, &                                !< The phase of a tidal constituent at time 0 [rad]
    tide_fn, &                               !< Amplitude modulation of tides by nodal cycle [nondim].
    tide_un                                  !< Phase modulation of tides by nodal cycle [rad].
  integer :: nc                              !< The number of tidal constituents in use
  integer :: length                          !< Number of fields of which harmonic analysis is to be performed
  character(len=4), allocatable, dimension(:) :: const_name !< The name of each constituent
  character(len=255) :: path                 !< Path to directory where output will be written
  type(unit_scale_type)  :: US               !< A dimensional unit scaling type
  type(HA_node), pointer :: list => NULL()   !< A linked list for storing the HA info of different fields
end type harmonic_analysis_CS

contains

!> This subroutine sets static variables used by this module and initializes CS%list.
!! THIS MUST BE CALLED AT THE END OF tidal_forcing_init.
subroutine HA_init(Time, US, param_file, nc, CS)
  type(time_type),       intent(in)  :: Time        !< The current model time
  type(unit_scale_type), intent(in)  :: US          !< A dimensional unit scaling type
  type(param_file_type), intent(in)  :: param_file  !< A structure to parse for run-time parameters
  integer,               intent(in)  :: nc          !< The number of tidal constituents in use
  type(harmonic_analysis_CS), intent(out) :: CS     !< Control structure of the MOM_harmonic_analysis module

  ! Local variables
  logical :: tides                                  !< True if tidal forcing module is enabled
  logical :: use_eq_phase                           !< If true, tidal forcing is phase-shifted to match
                                                    !! equilibrium tide. Set to false if providing tidal phases
                                                    !! that have already been shifted by the
                                                    !! astronomical/equilibrium argument
  logical :: add_nodal_terms                        !< If true, insert terms for the 18.6 year modulation when
                                                    !! calculating tidal forcing.
  integer, dimension(3)  :: tide_ref_date           !< Reference date (t = 0) for tidal forcing (year, month, day)
  integer, dimension(3)  :: nodal_ref_date          !< Date to calculate nodal modulation for (year, month, day)
  type(time_type)        :: nodal_time              !< Model time to calculate nodal modulation for.
  type(astro_longitudes) :: tidal_longitudes        !< Astronomical longitudes used to calculate
                                                    !! tidal phases at t = 0.
  type(astro_longitudes) :: nodal_longitudes        !< Solar and lunar longitudes for tidal forcing
  character(len=50)      :: const_name              !< Names of all tidal constituents to be harmonically analyzed
  integer :: c

  type(HA_type) :: ha1                              !< A temporary, null field used for initializing CS%list
  real :: HA_start_time                             !< Start time of harmonic analysis [T ~> s]
  real :: HA_end_time                               !< End time of harmonic analysis [T ~> s]
  logical :: HA_ssh, HA_ubt, HA_vbt
  character(len=40)  :: mdl="MOM_harmonic_analysis" !< This module's name
  character(len=255) :: mesg
  integer :: year, month, day, hour, minute, second

  call get_param(param_file, mdl, "TIDES", tides, &
      "If true, apply tidal momentum forcing.", default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "TIDE_USE_EQ_PHASE", use_eq_phase, &
      "If true, add the equilibrium phase argument to the specified tidal phases.", &
      old_name="OBC_TIDE_ADD_EQ_PHASE", default=.false., do_not_log=tides)
  call get_param(param_file, mdl, "TIDE_ADD_NODAL", add_nodal_terms, &
      "If true, include 18.6 year nodal modulation in the boundary tidal forcing.", &
      old_name="OBC_TIDE_ADD_NODAL", default=.false., do_not_log=tides)
  call get_param(param_file, mdl, "TIDE_REF_DATE", tide_ref_date, &
      "Reference date to use for tidal calculations and equilibrium phase.", &
      old_name="OBC_TIDE_REF_DATE", defaults=(/0, 0, 0/), do_not_log=tides)
  call get_param(param_file, mdl, "TIDE_NODAL_REF_DATE", nodal_ref_date, &
      "Fixed reference date to use for nodal modulation.", &
      old_name="OBC_TIDE_NODAL_REF_DATE", defaults=(/0, 0, 0/), do_not_log=tides)
  call get_param(param_file, mdl, "HA_CONSTITUENTS", const_name, &
      "Names of tidal constituents to be harmonically analyzed. "//&
      "They don't have to be the same as those used in MOM_tidal_forcing.", &
      fail_if_missing=.true.)

  if (sum(tide_ref_date) == 0) then  ! tide_ref_date defaults to 0.
    CS%time_ref = set_date(1, 1, 1, 0, 0, 0)
  else
    if (.not. use_eq_phase) then
      ! Using a reference date but not using phase relative to equilibrium.
      ! This makes sense as long as either phases are overridden, or
      ! correctly simulating tidal phases is not desired.
      call MOM_mesg('Tidal phases will *not* be corrected with equilibrium arguments.')
    endif
    CS%time_ref = set_date(tide_ref_date(1), tide_ref_date(2), tide_ref_date(3), 0, 0, 0)
  endif

  ! Initialize reference time for tides and find relevant lunar and solar
  ! longitudes at the reference time.
  if (use_eq_phase) call astro_longitudes_init(CS%time_ref, tidal_longitudes)

  ! If the nodal correction is based on a different time, initialize that.
  ! Otherwise, it can use N from the time reference.
  if (add_nodal_terms) then
    if (sum(nodal_ref_date) /= 0) then
      ! A reference date was provided for the nodal correction
      nodal_time = set_date(nodal_ref_date(1), nodal_ref_date(2), nodal_ref_date(3))
      call astro_longitudes_init(nodal_time, nodal_longitudes)
    elseif (use_eq_phase) then
      ! Astronomical longitudes were already calculated for use in equilibrium phases,
      ! so use nodal longitude from that.
      nodal_longitudes = tidal_longitudes
    else
      ! Tidal reference time is a required parameter, so calculate the longitudes from that.
      call astro_longitudes_init(CS%time_ref, nodal_longitudes)
    endif
  endif

  allocate(CS%const_name(nc))
  allocate(CS%freq(nc))
  allocate(CS%phase0(nc))
  allocate(CS%tide_fn(nc))
  allocate(CS%tide_un(nc))

  ! Tidal constituents for harmonic analysis can be different from those defined in MOM_tidal_forcing
  read(const_name, *) CS%const_name

  ! For major tidal constituents, tidal parameters defined in MOM_tidal_forcing will be used.
  ! For those not available in MOM_tidal_forcing, parameters needs to be defined in MOM_input.
  do c=1,nc
    call get_param(param_file, mdl, "HA_"//trim(CS%const_name(c))//"_FREQ", &
                   CS%freq(c), "Frequency of the "//trim(CS%const_name(c))//&
                   " constituent. This is used if USE_HA is true and "//trim(CS%const_name(c))//&
                   " is in HA_CONSTITUENTS.", units="rad s-1", scale=US%T_to_s, default=0.0)
    if (CS%freq(c)<=0.0) then
      select case (trim(CS%const_name(c)))
        case ('M4')
          CS%freq(c) = tidal_frequency('M2') * 2
        case ('M6')
          CS%freq(c) = tidal_frequency('M2') * 3
        case ('M8')
          CS%freq(c) = tidal_frequency('M2') * 4
        case ('S4')
          CS%freq(c) = tidal_frequency('S2') * 2
        case ('S6')
          CS%freq(c) = tidal_frequency('S2') * 3
        case ('MK3')
          CS%freq(c) = tidal_frequency('M2') + tidal_frequency('K1')
        case ('MS4')
          CS%freq(c) = tidal_frequency('M2') + tidal_frequency('S2')
        case ('MN4')
          CS%freq(c) = tidal_frequency('M2') + tidal_frequency('N2')
        case default
          CS%freq(c) = tidal_frequency(trim(CS%const_name(c)))
      end select
    endif

    call get_param(param_file, mdl, "HA_"//trim(CS%const_name(c))//"_PHASE_T0", CS%phase0(c), &
                   "Phase of the "//trim(CS%const_name(c))//" tidal constituent at time 0. "//&
                   "This is only used if USE_HA is true and "//trim(CS%const_name(c))// &
                   " is in HA_CONSTITUENTS.", units="radians", default=0.0)
    if (use_eq_phase) CS%phase0(c) = eq_phase(trim(CS%const_name(c)), tidal_longitudes)

    ! Nodal modulation should be turned off for tidal constituents not available in MOM_tidal_forcing
    if (add_nodal_terms) then
      call nodal_fu(trim(trim(CS%const_name(c))), nodal_longitudes%N, CS%tide_fn(c), CS%tide_un(c))
    else
      CS%tide_fn(c) = 1.0
      CS%tide_un(c) = 0.0
    endif
  enddo

  ! Determine CS%time_start and CS%time_end
  call get_param(param_file, mdl, "HA_START_TIME", HA_start_time, &
                 "Start time of harmonic analysis, in units of days after "//&
                 "the start of the current run segment. Must be smaller than "//&
                 "HA_END_TIME, otherwise harmonic analysis will not be performed. "//&
                 "If negative, |HA_START_TIME| determines the length of harmonic analysis, "//&
                 "and harmonic analysis will start |HA_START_TIME| days before HA_END_TIME, "//&
                 "or at the beginning of the run segment, whichever occurs later.", &
                 units="days", default=0.0, scale=86400.0*US%s_to_T)
  call get_param(param_file, mdl, "HA_END_TIME", HA_end_time, &
                 "End time of harmonic analysis, in units of days after "//&
                 "the start of the current run segment. Must be positive "//&
                 "and smaller than the length of the currnet run segment, "//&
                 "otherwise harmonic analysis will not be performed.", &
                 units="days", default=0.0, scale=86400.0*US%s_to_T)

  if (HA_end_time <= 0.0) then
    call MOM_mesg('MOM_harmonic_analysis: HA_END_TIME is zero or negative. '//&
                  'Harmonic analysis will not be performed.')
    CS%HAready = .false. ; return
  endif

  if (HA_end_time <= HA_start_time) then
    call MOM_mesg('MOM_harmonic_analysis: HA_END_TIME is smaller than or equal to HA_START_TIME. '//&
                  'Harmonic analysis will not be performed.')
    CS%HAready = .false. ; return
  endif

  CS%HAready = .true.

  if (HA_start_time < 0.0) then
    HA_start_time = HA_end_time + HA_start_time
    if (HA_start_time <= 0.0) HA_start_time = 0.0
  endif

  CS%time_start = Time + real_to_time(HA_start_time, unscale=US%T_to_s)
  CS%time_end = Time + real_to_time(HA_end_time, unscale=US%T_to_s)

  call get_date(Time, year, month, day, hour, minute, second)
  write(mesg,*) "MOM_harmonic_analysis: run segment starts on ", year, month, day, hour, minute, second
  call MOM_error(NOTE, trim(mesg))
  call get_date(CS%time_start, year, month, day, hour, minute, second)
  write(mesg,*) "MOM_harmonic_analysis: harmonic analysis starts on ", year, month, day, hour, minute, second
  call MOM_error(NOTE, trim(mesg))
  call get_date(CS%time_end, year, month, day, hour, minute, second)
  write(mesg,*) "MOM_harmonic_analysis: harmonic analysis ends on ", year, month, day, hour, minute, second
  call MOM_error(NOTE, trim(mesg))

  ! Set path to directory where output will be written
  call get_param(param_file, mdl, "HA_PATH", CS%path, &
                 "Path to output files for runtime harmonic analysis.", default="./")

  ! Populate some parameters of the control structure
  CS%nc         =  nc
  CS%length     =  0
  CS%US         =  US

  ! Initialize CS%list
  allocate(CS%list)
  CS%list%this  =  ha1
  nullify(CS%list%next)

  ! Register variables/fields to be analyzed
  call get_param(param_file, mdl, "HA_SSH", HA_ssh, &
                 "If true, perform harmonic analysis of sea serface height.", default=.false.)
  if (HA_ssh) call HA_register('ssh', 'h', CS)
  call get_param(param_file, mdl, "HA_UBT", HA_ubt, &
                 "If true, perform harmonic analysis of zonal barotropic velocity.", default=.false.)
  if (HA_ubt) call HA_register('ubt', 'u', CS)
  call get_param(param_file, mdl, "HA_VBT", HA_vbt, &
                 "If true, perform harmonic analysis of meridional barotropic velocity.", default=.false.)
  if (HA_vbt) call HA_register('vbt', 'v', CS)

end subroutine HA_init

!> This subroutine registers each of the fields on which HA is to be performed.
subroutine HA_register(key, grid, CS)
  character(len=*),           intent(in)    :: key     !< Name of the current field
  character(len=1),           intent(in)    :: grid    !< The grid on which the key field is defined
  type(harmonic_analysis_CS), intent(inout) :: CS      !< Control structure of the MOM_harmonic_analysis module

  ! Local variables
  type(HA_type)          :: ha1                        !< Control structure for the current field
  type(HA_node), pointer :: tmp                        !< A temporary list to hold the current field

  if (.not. CS%HAready) return

  allocate(tmp)
  ha1%key   =  trim(key)
  ha1%grid  =  trim(grid)
  tmp%this  =  ha1
  tmp%next  => CS%list
  CS%list   => tmp
  CS%length =  CS%length + 1

end subroutine HA_register

!> This subroutine accumulates the temporal basis functions in FtF and FtSSH and then calls HA_write to compute
!! harmonic constants and write results. The tidal constituents are those used in MOM_tidal_forcing, plus the
!! mean (of zero frequency). For FtF, only the main diagonal and entries below it are calculated, which are needed
!! for Cholesky decomposition.
subroutine HA_accum(key, data, Time, G, CS)
  character(len=*),           intent(in) :: key  !< Name of the current field
  real, dimension(:,:),       intent(in) :: data !< Input data of which harmonic analysis is to be performed [A]
  type(time_type),            intent(in) :: Time !< The current model time
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(harmonic_analysis_CS), intent(inout) :: CS   !< Control structure of the MOM_harmonic_analysis module

  ! Local variables
  type(HA_type), pointer :: ha1
  type(HA_node), pointer :: tmp
  real :: now                                    !< The relative time compared with the tidal reference [T ~> s]
  real :: dt                                     !< The current time step size of the accumulator [T ~> s]
  real :: cosomegat, sinomegat, ccosomegat, ssinomegat !< The components of the phase [nondim]
  integer :: nc, i, j, k, c, cc, icos, isin, iccos, issin, is, ie, js, je
  character(len=128) :: mesg

  ! Exit the accumulator in the following cases
  if (.not. CS%HAready) return
  if (CS%length == 0) return
  if (Time < CS%time_start) return
  if (Time > CS%time_end) return

  ! Loop through the full list to find the current field
  tmp => CS%list
  do k=1,CS%length
    ha1 => tmp%this
    if (trim(key) == trim(ha1%key)) exit
    tmp => tmp%next
    if (k == CS%length) return              !< Do not perform harmonic analysis of a field that is not registered
  enddo

  nc  = CS%nc
  now = time_minus_signed(Time, CS%time_ref, scale=CS%US%s_to_T)

  !!! Additional processing at the initial accumulating step !!!
  if (ha1%old_time < 0.0) then
    ha1%old_time = now

    write(mesg,*) "MOM_harmonic_analysis: initializing accumulator, key = ", trim(ha1%key)
    call MOM_error(NOTE, trim(mesg))

    ! Get the lower and upper bounds of input data
    ha1%is = LBOUND(data,1) ; is = ha1%is
    ha1%ie = UBOUND(data,1) ; ie = ha1%ie
    ha1%js = LBOUND(data,2) ; js = ha1%js
    ha1%je = UBOUND(data,2) ; je = ha1%je

    allocate(ha1%ref(is:ie,js:je), source=0.0)
    allocate(ha1%FtF(2*nc+1,2*nc+1), source=0.0)
    allocate(ha1%FtSSH(is:ie,js:je,2*nc+1), source=0.0)
    ha1%ref(:,:) = data(:,:)
  endif

  dt = now - ha1%old_time
  ha1%old_time = now                        !< Keep track of time so we know when Time approaches CS%time_end

  is = ha1%is ; ie = ha1%ie ; js = ha1%js ; je = ha1%je

  !!! Accumulator of FtF !!!
  !< First entry, corresponding to the zero frequency constituent (mean)
  ha1%FtF(1,1) = ha1%FtF(1,1) + 1.0

  do c=1,nc
    icos = 2*c
    isin = 2*c+1
    cosomegat = CS%tide_fn(c) * cos(CS%freq(c) * now + (CS%phase0(c) + CS%tide_un(c)))
    sinomegat = CS%tide_fn(c) * sin(CS%freq(c) * now + (CS%phase0(c) + CS%tide_un(c)))

    ! First column, corresponding to the zero frequency constituent (mean)
    ha1%FtF(icos,1) = ha1%FtF(icos,1) + cosomegat
    ha1%FtF(isin,1) = ha1%FtF(isin,1) + sinomegat

    do cc=1,c
      iccos = 2*cc
      issin = 2*cc+1
      ccosomegat = CS%tide_fn(cc) * cos(CS%freq(cc) * now + (CS%phase0(cc) + CS%tide_un(cc)))
      ssinomegat = CS%tide_fn(cc) * sin(CS%freq(cc) * now + (CS%phase0(cc) + CS%tide_un(cc)))

      ! Interior of the matrix, corresponding to the products of cosine and sine terms
      ha1%FtF(icos,iccos) = ha1%FtF(icos,iccos) + cosomegat * ccosomegat
      ha1%FtF(icos,issin) = ha1%FtF(icos,issin) + cosomegat * ssinomegat
      ha1%FtF(isin,iccos) = ha1%FtF(isin,iccos) + sinomegat * ccosomegat
      ha1%FtF(isin,issin) = ha1%FtF(isin,issin) + sinomegat * ssinomegat
    enddo ! cc=1,c
  enddo ! c=1,nc

  !!! Accumulator of FtSSH !!!
  !< First entry, corresponding to the zero frequency constituent (mean)
  do j=js,je ; do i=is,ie
    ha1%FtSSH(i,j,1) = ha1%FtSSH(i,j,1) + (data(i,j) - ha1%ref(i,j))
  enddo ; enddo

  !< The remaining entries
  do c=1,nc
    icos = 2*c
    isin = 2*c+1
    cosomegat = CS%tide_fn(c) * cos(CS%freq(c) * now + (CS%phase0(c) + CS%tide_un(c)))
    sinomegat = CS%tide_fn(c) * sin(CS%freq(c) * now + (CS%phase0(c) + CS%tide_un(c)))
    do j=js,je ; do i=is,ie
      ha1%FtSSH(i,j,icos) = ha1%FtSSH(i,j,icos) + (data(i,j) - ha1%ref(i,j)) * cosomegat
      ha1%FtSSH(i,j,isin) = ha1%FtSSH(i,j,isin) + (data(i,j) - ha1%ref(i,j)) * sinomegat
    enddo ; enddo
  enddo ! c=1,nc

  !!! Compute harmonic constants and write output as Time approaches CS%time_end !!!
  ! This guarantees that HA_write will be called before Time becomes larger than CS%time_end.
  ! Result of subtracting time types is always >= 0, which is acceptable here.
  if (time_to_real(CS%time_end - Time, scale=CS%US%s_to_T) <= dt) then
    call HA_write(ha1, Time, G, CS)

    write(mesg,*) "MOM_harmonic_analysis: harmonic analysis done, key = ", trim(ha1%key)
    call MOM_error(NOTE, trim(mesg))

    ! De-register the current field and deallocate memory
    ha1%key = 'none'
    deallocate(ha1%ref)
    deallocate(ha1%FtSSH)
  endif

end subroutine HA_accum

!> This subroutine computes the harmonic constants and write output for the current field
subroutine HA_write(ha1, Time, G, CS)
  type(HA_type), pointer,     intent(in) :: ha1    !< Control structure for the current field
  type(time_type),            intent(in) :: Time   !< The current model time
  type(ocean_grid_type),      intent(in) :: G      !< The ocean's grid structure
  type(harmonic_analysis_CS), intent(in) :: CS     !< Control structure of the MOM_harmonic_analysis module

  ! Local variables
  real, dimension(:,:,:), allocatable :: FtSSHw    !< An array containing the harmonic constants [A]
  integer :: year, month, day, hour, minute, second
  integer :: nc, i, j, k, is, ie, js, je

  character(len=255)           :: filename         !< Output file name
  type(MOM_infra_file)         :: cdf              !< The file handle for output harmonic constants
  type(vardesc),   allocatable :: cdf_vars(:)      !< Output variable names
  type(MOM_field), allocatable :: cdf_fields(:)    !< Field type variables for the output fields

  nc = CS%nc ; is = ha1%is ; ie = ha1%ie ; js = ha1%js ; je = ha1%je

  allocate(FtSSHw(is:ie,js:je,2*nc+1), source=0.0)

  ! Compute the harmonic coefficients
  call HA_solver(ha1, nc, ha1%FtF, FtSSHw)

  ! Output file name
  call get_date(Time, year, month, day, hour, minute, second)
  write(filename, '(a,"HA_",a,i0.4,i0.2,i0.2,".nc")') &
      trim(CS%path), trim(ha1%key), year, month, day

  allocate(cdf_vars(2*nc+1))
  allocate(cdf_fields(2*nc+1))

  ! Variable names
  cdf_vars(1) = var_desc("z0", "m" ,"mean value", ha1%grid, '1', '1')
  do k=1,nc
    cdf_vars(2*k  ) = var_desc(trim(CS%const_name(k))//"cos", "m", "cosine coefficient", ha1%grid, '1', '1')
    cdf_vars(2*k+1) = var_desc(trim(CS%const_name(k))//"sin", "m", "sine coefficient", ha1%grid, '1', '1')
  enddo

  ! Create output file
  call create_MOM_file(cdf, trim(filename), cdf_vars, &
                       2*nc+1, cdf_fields, SINGLE_FILE, 86400.0, G=G)

  ! Add the initial field back to the mean state
  do j=js,je ; do i=is,ie
    FtSSHw(i,j,1) = FtSSHw(i,j,1) + ha1%ref(i,j)
  enddo ; enddo

  ! Write data
  call MOM_write_field(cdf, cdf_fields(1), G%domain, FtSSHw(:,:,1), 0.0)
  do k=1,nc
    call MOM_write_field(cdf, cdf_fields(2*k  ), G%domain, FtSSHw(:,:,2*k  ), 0.0)
    call MOM_write_field(cdf, cdf_fields(2*k+1), G%domain, FtSSHw(:,:,2*k+1), 0.0)
  enddo

  call cdf%flush()
  deallocate(cdf_vars)
  deallocate(cdf_fields)
  deallocate(FtSSHw)

end subroutine HA_write

!> This subroutine computes the harmonic constants (stored in x) using the dot products of the temporal
!! basis functions accumulated in FtF, and the dot products of the SSH (or other fields) with the temporal basis
!! functions accumulated in FtSSH. The system is solved by Cholesky decomposition,
!!
!!     FtF * x = FtSSH,    =>    L * (L' * x) = FtSSH,    =>    L * y = FtSSH,
!!
!! where L is the lower triangular matrix, y = L' * x, and x is the solution vector.
!!
subroutine HA_solver(ha1, nc, FtF, x)
  type(HA_type), pointer,              intent(in)  :: ha1    !< Control structure for the current field
  integer,                             intent(in)  :: nc     !< Number of harmonic constituents
  real, dimension(:,:),                intent(in)  :: FtF    !< Accumulator of (F' * F) for all fields [nondim]
  real, dimension(ha1%is:ha1%ie,ha1%js:ha1%je,2*nc+1), &
                                       intent(out) :: x      !< Solution vector of harmonic constants [A]

  ! Local variables
  real :: tmp0                                !< Temporary variable for Cholesky decomposition [nondim]
  real, dimension(2*nc+1,2*nc+1)      :: L    !< Lower triangular matrix of Cholesky decomposition [nondim]
  real, dimension(2*nc+1)             :: tmp1 !< Inverse of the diagonal entries of L [nondim]
  real, dimension(ha1%is:ha1%ie,ha1%js:ha1%je)        :: tmp2 !< 2D temporary array involving FtSSH [A]
  real, dimension(ha1%is:ha1%ie,ha1%js:ha1%je,2*nc+1) :: y    !< 3D temporary array, i.e., L' * x [A]
  integer :: k, m, n

  ! Cholesky decomposition
  do m=1,2*nc+1

    ! First, calculate the diagonal entries
    tmp0 = 0.0
    do k=1,m-1                             ! This loop operates along the m-th row
      tmp0 = tmp0 + L(m,k) * L(m,k)
    enddo
    L(m,m) = sqrt(FtF(m,m) - tmp0)         ! This is the m-th diagonal entry

    ! Now calculate the off-diagonal entries
    tmp1(m) = 1 / L(m,m)
    do k=m+1,2*nc+1                        ! This loop operates along the column below the m-th diagonal entry
      tmp0 = 0.0
      do n=1,m-1
        tmp0 = tmp0 + L(k,n) * L(m,n)
      enddo
      L(k,m) = (FtF(k,m) - tmp0) * tmp1(m) ! This is the k-th off-diagonal entry below the m-th diagonal entry
    enddo
  enddo

  ! Solve for y from L * y = FtSSH
  do k=1,2*nc+1
    tmp2(:,:) = 0.0
    do m=1,k-1
      tmp2(:,:) = tmp2(:,:) + L(k,m) * y(:,:,m)
    enddo
    y(:,:,k) = (ha1%FtSSH(:,:,k) - tmp2(:,:)) * tmp1(k)
  enddo

  ! Solve for x from L' * x = y
  do k=2*nc+1,1,-1
    tmp2(:,:) = 0.0
    do m=k+1,2*nc+1
      tmp2(:,:) = tmp2(:,:) + L(m,k) * x(:,:,m)
    enddo
    x(:,:,k) = (y(:,:,k) - tmp2(:,:)) * tmp1(k)
  enddo

end subroutine HA_solver

!> \namespace harmonic_analysis
!!
!! Major revision (August, 2025)
!!
!! This module is now independent of MOM_tidal_forcing, providing more flexibility for performing harmonic analyses
!! on tidal constituents not available in MOM_tidal_forcing (e.g., MK3, M4), with the following conditions:
!!   1) For tidal constituents not available in MOM_tidal_forcing, the frequencies and equilibrium phases (if used)
!!      must be specified manually in MOM_input.
!!   2) If any tidal constituents not available in MOM_tidal_forcing are used, the nodal modulation cannot be added.
!!      Or, if nodal modulation is added, then harmonic analysis can only be performed on major tidal constituents.
!!
!! Original version (April, 2024)
!!
!! This module computes the harmonic constants which can be used to reconstruct the tidal elevation (or other
!! fields) through SSH = F * x, where F is an nt-by-2*nc matrix (nt is the number of time steps and nc is the
!! number of tidal constituents) containing the cosine/sine functions for each frequency evaluated at each time
!! step, and x is a 2*nc-by-1 vector containing the constant coefficients of the sine/cosine for each constituent
!! (i.e., the harmonic constants). At each grid point, the harmonic constants are computed using least squares,
!!
!!     (F' * F) * x = F' * SSH_in,    =>    FtF * x = FtSSH,
!!
!! where the prime denotes matrix transpose, and SSH_in is the sea surface height (or other fields) determined by
!! the model. The dot products (F' * F) and (F' * SSH_in) are computed by accumulating the sums as the model is
!! running and stored in the arrays FtF and FtSSH, respectively. The FtF matrix is inverted as needed before
!! computing and writing out the harmonic constants.
!!
!! Ed Zaron and William Xu (chengzhu.xu@oregonstate.edu)

end module MOM_harmonic_analysis

