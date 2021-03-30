!> Describes the decomposed MOM domain and has routines for communications across PEs
module MOM_domains

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms_infra,       only : MOM_infra_init, MOM_infra_end
use MOM_coms_infra,       only : PE_here, root_PE, num_PEs, broadcast
use MOM_coms_infra,       only : sum_across_PEs, min_across_PEs, max_across_PEs
use MOM_domain_infra,     only : MOM_domain_type, domain2D, domain1D, group_pass_type
use MOM_domain_infra,     only : create_MOM_domain, clone_MOM_domain, deallocate_MOM_domain
use MOM_domain_infra,     only : get_domain_extent, get_domain_components
use MOM_domain_infra,     only : compute_block_extent, get_global_shape
use MOM_domain_infra,     only : pass_var, pass_vector, fill_symmetric_edges, global_field_sum
use MOM_domain_infra,     only : pass_var_start, pass_var_complete
use MOM_domain_infra,     only : pass_vector_start, pass_vector_complete
use MOM_domain_infra,     only : create_group_pass, do_group_pass
use MOM_domain_infra,     only : start_group_pass, complete_group_pass
use MOM_domain_infra,     only : rescale_comp_data, global_field, redistribute_array, broadcast_domain
use MOM_domain_infra,     only : MOM_thread_affinity_set, set_MOM_thread_affinity
use MOM_domain_infra,     only : AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR, BITWISE_EXACT_SUM
use MOM_domain_infra,     only : CORNER, CENTER, NORTH_FACE, EAST_FACE
use MOM_domain_infra,     only : To_East, To_West, To_North, To_South, To_All, Omit_Corners
use MOM_error_handler,    only : MOM_error, MOM_mesg, NOTE, WARNING, FATAL
use MOM_file_parser,      only : get_param, log_param, log_version, param_file_type
use MOM_io_infra,         only : file_exists
use MOM_string_functions, only : slasher

implicit none ; private

public :: MOM_infra_init, MOM_infra_end
!  Domain types and creation and destruction routines
public :: MOM_domain_type, domain2D, domain1D
public :: MOM_domains_init, create_MOM_domain, clone_MOM_domain, deallocate_MOM_domain
public :: MOM_thread_affinity_set, set_MOM_thread_affinity
!  Domain query routines
public :: get_domain_extent, get_domain_components, compute_block_extent, get_global_shape
public :: PE_here, root_PE, num_PEs
!  Single call communication routines
public :: pass_var, pass_vector, fill_symmetric_edges, broadcast
!  Non-blocking communication routines
public :: pass_var_start, pass_var_complete, pass_vector_start, pass_vector_complete
!  Multi-variable group communication routines and type
public :: create_group_pass, do_group_pass, group_pass_type, start_group_pass, complete_group_pass
!  Global reduction routines
public :: sum_across_PEs, min_across_PEs, max_across_PEs
public :: global_field, redistribute_array, broadcast_domain
!  Simple index-convention-invariant array manipulation routine
public :: rescale_comp_data
!> These encoding constants are used to indicate the staggering of scalars and vectors
public :: AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR
!> These encoding constants are used to indicate the discretization position of a variable
public :: CORNER, CENTER, NORTH_FACE, EAST_FACE
!> These encoding constants indicate communication patterns.  In practice they can be added.
public :: To_East, To_West, To_North, To_South, To_All, Omit_Corners
! These are no longer used by MOM6 because the reproducing sum works so well, but they are
! still referenced by some of the non-GFDL couplers.
public :: global_field_sum, BITWISE_EXACT_SUM

contains

!> MOM_domains_init initializes a MOM_domain_type variable, based on the information
!! read in from a param_file_type, and optionally returns data describing various
!! properties of the domain type.
subroutine MOM_domains_init(MOM_dom, param_file, symmetric, static_memory, &
                            NIHALO, NJHALO, NIGLOBAL, NJGLOBAL, NIPROC, NJPROC, &
                            min_halo, domain_name, include_name, param_suffix)
  type(MOM_domain_type),           pointer       :: MOM_dom      !< A pointer to the MOM_domain_type
                                                                 !! being defined here.
  type(param_file_type),           intent(in)    :: param_file   !< A structure to parse for
                                                                 !! run-time parameters
  logical, optional,               intent(in)    :: symmetric    !< If present, this specifies
                                                  !! whether this domain is symmetric, regardless of
                                                  !! whether the macro SYMMETRIC_MEMORY_ is defined.
  logical, optional,               intent(in)    :: static_memory !< If present and true, this
                                                  !! domain type is set up for static memory and
                                                  !! error checking of various input values is
                                                  !! performed against those in the input file.
  integer, optional,               intent(in)    :: NIHALO       !< Default halo sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NJHALO       !< Default halo sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NIGLOBAL     !< Total domain sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NJGLOBAL     !< Total domain sizes, required
                                                                 !! with static memory.
  integer, optional,               intent(in)    :: NIPROC       !< Processor counts, required with
                                                                 !! static memory.
  integer, optional,               intent(in)    :: NJPROC       !< Processor counts, required with
                                                                 !! static memory.
  integer, dimension(2), optional, intent(inout) :: min_halo     !< If present, this sets the
                                            !! minimum halo size for this domain in the i- and j-
                                            !! directions, and returns the actual halo size used.
  character(len=*),      optional, intent(in)    :: domain_name  !< A name for this domain, "MOM"
                                                                 !! if missing.
  character(len=*),      optional, intent(in)    :: include_name !< A name for model's include file,
                                                                 !! "MOM_memory.h" if missing.
  character(len=*),      optional, intent(in)    :: param_suffix !< A suffix to apply to
                                                                 !! layout-specific parameters.

  ! Local variables
  integer, dimension(2) :: layout    ! The number of logical processors in the i- and j- directions
  integer, dimension(2) :: io_layout ! The layout of logical processors for input and output
  !$ integer :: ocean_nthreads       ! Number of openMP threads
  !$ logical :: ocean_omp_hyper_thread ! If true use openMP hyper-threads
  integer, dimension(2) :: n_global  ! The number of i- and j- points in the global computational domain
  integer, dimension(2) :: n_halo    ! The number of i- and j- points in the halos
  integer :: nihalo_dflt, njhalo_dflt ! The default halo sizes
  integer :: PEs_used       ! The number of processors used
  logical, dimension(2) :: reentrant ! True if the x- and y- directions are periodic.
  logical :: tripolar_N     ! A flag indicating whether there is northern tripolar connectivity
  logical :: is_static      ! If true, static memory is being used for this domain.
  logical :: is_symmetric   ! True if the domain being set up will use symmetric memory.
  logical :: nonblocking    ! If true, nonblocking halo updates will be used.
  logical :: thin_halos     ! If true, If true, optional arguments may be used to specify the
                            ! width of the halos that are updated with each call.
  logical            :: mask_table_exists ! True if there is a mask table file
  character(len=128) :: inputdir   ! The directory in which to find the diag table
  character(len=200) :: mask_table ! The file name and later the full path to the diag table
  character(len=64)  :: inc_nm     ! The name of the memory include file
  character(len=200) :: mesg       ! A string to use for error messages

  integer :: nip_parsed, njp_parsed
  character(len=8) :: char_xsiz, char_ysiz, char_niglobal, char_njglobal
  character(len=40) :: nihalo_nm, njhalo_nm, layout_nm, io_layout_nm, masktable_nm
  character(len=40) :: niproc_nm, njproc_nm
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl ! This module's name.

  PEs_used = num_PEs()

  mdl = "MOM_domains"

  is_symmetric = .true. ; if (present(symmetric)) is_symmetric = symmetric
  if (present(min_halo)) mdl = trim(mdl)//" min_halo"

  inc_nm = "MOM_memory.h" ; if (present(include_name)) inc_nm = trim(include_name)

  nihalo_nm = "NIHALO" ; njhalo_nm = "NJHALO"
  layout_nm = "LAYOUT" ; io_layout_nm = "IO_LAYOUT" ; masktable_nm = "MASKTABLE"
  niproc_nm = "NIPROC" ; njproc_nm = "NJPROC"
  if (present(param_suffix)) then ; if (len(trim(adjustl(param_suffix))) > 0) then
    nihalo_nm = "NIHALO"//(trim(adjustl(param_suffix)))
    njhalo_nm = "NJHALO"//(trim(adjustl(param_suffix)))
    layout_nm = "LAYOUT"//(trim(adjustl(param_suffix)))
    io_layout_nm = "IO_LAYOUT"//(trim(adjustl(param_suffix)))
    masktable_nm = "MASKTABLE"//(trim(adjustl(param_suffix)))
    niproc_nm = "NIPROC"//(trim(adjustl(param_suffix)))
    njproc_nm = "NJPROC"//(trim(adjustl(param_suffix)))
  endif ; endif

  is_static = .false. ; if (present(static_memory)) is_static = static_memory
  if (is_static) then
    if (.not.present(NIHALO)) call MOM_error(FATAL, "NIHALO must be "// &
      "present in the call to MOM_domains_init with static memory.")
    if (.not.present(NJHALO)) call MOM_error(FATAL, "NJHALO must be "// &
      "present in the call to MOM_domains_init with static memory.")
    if (.not.present(NIGLOBAL)) call MOM_error(FATAL, "NIGLOBAL must be "// &
      "present in the call to MOM_domains_init with static memory.")
    if (.not.present(NJGLOBAL)) call MOM_error(FATAL, "NJGLOBAL must be "// &
      "present in the call to MOM_domains_init with static memory.")
    if (.not.present(NIPROC)) call MOM_error(FATAL, "NIPROC must be "// &
      "present in the call to MOM_domains_init with static memory.")
    if (.not.present(NJPROC)) call MOM_error(FATAL, "NJPROC must be "// &
      "present in the call to MOM_domains_init with static memory.")
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "", log_to_all=.true., layout=.true.)
  call get_param(param_file, mdl, "REENTRANT_X", reentrant(1), &
                 "If true, the domain is zonally reentrant.", default=.true.)
  call get_param(param_file, mdl, "REENTRANT_Y", reentrant(2), &
                 "If true, the domain is meridionally reentrant.", &
                 default=.false.)
  call get_param(param_file, mdl, "TRIPOLAR_N", tripolar_N, &
                 "Use tripolar connectivity at the northern edge of the "//&
                 "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", &
                 default=.false.)

# ifndef NOT_SET_AFFINITY
  !$ if (.not.MOM_thread_affinity_set()) then
  !$   call get_param(param_file, mdl, "OCEAN_OMP_THREADS", ocean_nthreads, &
  !$              "The number of OpenMP threads that MOM6 will use.", &
  !$              default = 1, layoutParam=.true.)
  !$   call get_param(param_file, mdl, "OCEAN_OMP_HYPER_THREAD", ocean_omp_hyper_thread, &
  !$              "If True, use hyper-threading.", default = .false., layoutParam=.true.)
  !$   call set_MOM_thread_affinity(ocean_nthreads, ocean_omp_hyper_thread)
  !$ endif
# endif

  call log_param(param_file, mdl, "!SYMMETRIC_MEMORY_", is_symmetric, &
                 "If defined, the velocity point data domain includes every face of the "//&
                 "thickness points. In other words, some arrays are larger than others, "//&
                 "depending on where they are on the staggered grid.  Also, the starting "//&
                 "index of the velocity-point arrays is usually 0, not 1. "//&
                 "This can only be set at compile time.",&
                 layoutParam=.true.)
  call get_param(param_file, mdl, "NONBLOCKING_UPDATES", nonblocking, &
                 "If true, non-blocking halo updates may be used.", &
                 default=.false., layoutParam=.true.)
  call get_param(param_file, mdl, "THIN_HALO_UPDATES", thin_halos, &
                 "If true, optional arguments may be used to specify the width of the "//&
                 "halos that are updated with each call.", &
                 default=.true., layoutParam=.true.)

  nihalo_dflt = 4 ; njhalo_dflt = 4
  if (present(NIHALO)) nihalo_dflt = NIHALO
  if (present(NJHALO)) njhalo_dflt = NJHALO

  call log_param(param_file, mdl, "!STATIC_MEMORY_", is_static, &
                 "If STATIC_MEMORY_ is defined, the principle variables will have sizes that "//&
                 "are statically determined at compile time.  Otherwise the sizes are not "//&
                 "determined until run time. The STATIC option is substantially faster, but "//&
                 "does not allow the PE count to be changed at run time.  This can only be "//&
                 "set at compile time.", layoutParam=.true.)

  if (is_static) then
    call get_param(param_file, mdl, "NIGLOBAL", n_global(1), &
                 "The total number of thickness grid points in the x-direction in the physical "//&
                 "domain. With STATIC_MEMORY_ this is set in "//trim(inc_nm)//" at compile time.", &
                 static_value=NIGLOBAL)
    call get_param(param_file, mdl, "NJGLOBAL", n_global(2), &
                 "The total number of thickness grid points in the y-direction in the physical "//&
                 "domain. With STATIC_MEMORY_ this is set in "//trim(inc_nm)//" at compile time.", &
                 static_value=NJGLOBAL)
    if (n_global(1) /= NIGLOBAL) call MOM_error(FATAL,"MOM_domains_init: " // &
          "static mismatch for NIGLOBAL_ domain size. Header file does not match input namelist")
    if (n_global(2) /= NJGLOBAL) call MOM_error(FATAL,"MOM_domains_init: " // &
          "static mismatch for NJGLOBAL_ domain size. Header file does not match input namelist")

    ! Check the requirement of equal sized compute domains when STATIC_MEMORY_ is used.
    if ((MOD(NIGLOBAL, NIPROC) /= 0) .OR. (MOD(NJGLOBAL, NJPROC) /= 0)) then
      write( char_xsiz, '(i4)' ) NIPROC
      write( char_ysiz, '(i4)' ) NJPROC
      write( char_niglobal, '(i4)' ) NIGLOBAL
      write( char_njglobal, '(i4)' ) NJGLOBAL
      call MOM_error(WARNING, 'MOM_domains: Processor decomposition (NIPROC_,NJPROC_) = ('//&
              trim(char_xsiz)//','//trim(char_ysiz)//') does not evenly divide size '//&
              'set by preprocessor macro ('//trim(char_niglobal)//','//trim(char_njglobal)//').')
      call MOM_error(FATAL,'MOM_domains:  #undef STATIC_MEMORY_ in '//trim(inc_nm)//' to use '//&
              'dynamic allocation, or change processor decomposition to evenly divide the domain.')
    endif
  else
    call get_param(param_file, mdl, "NIGLOBAL", n_global(1), &
                 "The total number of thickness grid points in the x-direction in the physical "//&
                 "domain. With STATIC_MEMORY_ this is set in "//trim(inc_nm)//" at compile time.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "NJGLOBAL", n_global(2), &
                 "The total number of thickness grid points in the y-direction in the physical "//&
                 "domain. With STATIC_MEMORY_ this is set in "//trim(inc_nm)//" at compile time.", &
                 fail_if_missing=.true.)
  endif

  call get_param(param_file, mdl, trim(nihalo_nm), n_halo(1), &
                 "The number of halo points on each side in the x-direction.  How this is set "//&
                 "varies with the calling component and static or dynamic memory configuration.", &
                 default=nihalo_dflt, static_value=nihalo_dflt)
  call get_param(param_file, mdl, trim(njhalo_nm), n_halo(2), &
                 "The number of halo points on each side in the y-direction.  How this is set "//&
                 "varies with the calling component and static or dynamic memory configuration.", &
                 default=njhalo_dflt, static_value=njhalo_dflt)
  if (present(min_halo)) then
    n_halo(1) = max(n_halo(1), min_halo(1))
    min_halo(1) = n_halo(1)
    n_halo(2) = max(n_halo(2), min_halo(2))
    min_halo(2) = n_halo(2)
    ! These are generally used only with static memory, so they are considered layout params.
    call log_param(param_file, mdl, "!NIHALO min_halo", n_halo(1), layoutParam=.true.)
    call log_param(param_file, mdl, "!NJHALO min_halo", n_halo(2), layoutParam=.true.)
  endif
  if (is_static .and. .not.present(min_halo)) then
    if (n_halo(1) /= NIHALO) call MOM_error(FATAL,"MOM_domains_init: " // &
           "static mismatch for "//trim(nihalo_nm)//" domain size")
    if (n_halo(2) /= NJHALO) call MOM_error(FATAL,"MOM_domains_init: " // &
           "static mismatch for "//trim(njhalo_nm)//" domain size")
  endif

  call get_param(param_file, mdl, "INPUTDIR", inputdir, do_not_log=.true., default=".")
  inputdir = slasher(inputdir)

  call get_param(param_file, mdl, trim(masktable_nm), mask_table, &
                 "A text file to specify n_mask, layout and mask_list. This feature masks out "//&
                 "processors that contain only land points. The first line of mask_table is the "//&
                 "number of regions to be masked out. The second line is the layout of the "//&
                 "model and must be consistent with the actual model layout. The following "//&
                 "(n_mask) lines give the logical positions of the processors that are masked "//&
                 "out. The mask_table can be created by tools like check_mask. The following "//&
                 "example of mask_table masks out 2 processors, (1,2) and (3,6), out of the 24 "//&
                 "in a 4x6 layout: \n 2\n 4,6\n 1,2\n 3,6\n", default="MOM_mask_table", &
                 layoutParam=.true.)
  mask_table = trim(inputdir)//trim(mask_table)
  mask_table_exists = file_exists(mask_table)

  if (is_static) then
    layout(1) = NIPROC ; layout(2) = NJPROC
  else
    call get_param(param_file, mdl, trim(layout_nm), layout, &
                 "The processor layout to be used, or 0, 0 to automatically set the layout "//&
                 "based on the number of processors.", default=0, do_not_log=.true.)
    call get_param(param_file, mdl, trim(niproc_nm), nip_parsed, &
                 "The number of processors in the x-direction.", default=-1, do_not_log=.true.)
    call get_param(param_file, mdl, trim(njproc_nm), njp_parsed, &
                 "The number of processors in the y-direction.", default=-1, do_not_log=.true.)
    if (nip_parsed > -1) then
      if ((layout(1) > 0) .and. (layout(1) /= nip_parsed)) &
        call MOM_error(FATAL, trim(layout_nm)//" and "//trim(niproc_nm)//" set inconsistently. "//&
                              "Only LAYOUT should be used.")
      layout(1) = nip_parsed
      call MOM_mesg(trim(niproc_nm)//" used to set "//trim(layout_nm)//" in dynamic mode.  "//&
                    "Shift to using "//trim(layout_nm)//" instead.")
    endif
    if (njp_parsed > -1) then
      if ((layout(2) > 0) .and. (layout(2) /= njp_parsed)) &
        call MOM_error(FATAL, trim(layout_nm)//" and "//trim(njproc_nm)//" set inconsistently. "//&
                              "Only "//trim(layout_nm)//" should be used.")
      layout(2) = njp_parsed
      call MOM_mesg(trim(njproc_nm)//" used to set "//trim(layout_nm)//" in dynamic mode.  "//&
                    "Shift to using "//trim(layout_nm)//" instead.")
    endif

    if ( (layout(1) == 0) .and. (layout(2) == 0) ) &
      call MOM_define_layout(n_global, PEs_used, layout)
    if ( (layout(1) /= 0) .and. (layout(2) == 0) ) layout(2) = PEs_used / layout(1)
    if ( (layout(1) == 0) .and. (layout(2) /= 0) ) layout(1) = PEs_used / layout(2)

    if (layout(1)*layout(2) /= PEs_used .and. (.not. mask_table_exists) ) then
      write(mesg,'("MOM_domains_init: The product of the two components of layout, ", &
            &      2i4,", is not the number of PEs used, ",i5,".")') &
            layout(1), layout(2), PEs_used
      call MOM_error(FATAL, mesg)
    endif
  endif
  call log_param(param_file, mdl, trim(niproc_nm), layout(1), &
                 "The number of processors in the x-direction. With STATIC_MEMORY_ this "//&
                 "is set in "//trim(inc_nm)//" at compile time.", layoutParam=.true.)
  call log_param(param_file, mdl, trim(njproc_nm), layout(2), &
                 "The number of processors in the y-direction. With STATIC_MEMORY_ this "//&
                 "is set in "//trim(inc_nm)//" at compile time.", layoutParam=.true.)
  call log_param(param_file, mdl, trim(layout_nm), layout, &
                 "The processor layout that was actually used.", layoutParam=.true.)

  ! Idiot check that fewer PEs than columns have been requested
  if (layout(1)*layout(2) > n_global(1)*n_global(2))  then
    write(mesg,'(a,2(i5,x,a))') 'You requested to use',layout(1)*layout(2), &
      'PEs but there are only', n_global(1)*n_global(2), 'columns in the model'
    call MOM_error(FATAL, mesg)
  endif

  if (mask_table_exists) &
    call MOM_error(NOTE, 'MOM_domains_init: reading maskmap information from '//trim(mask_table))

  ! Set up the I/O layout, it will be checked later that it uses an even multiple of the number of
  ! PEs in each direction.
  io_layout(:) = (/ 1, 1 /)
  call get_param(param_file, mdl, trim(io_layout_nm), io_layout, &
                 "The processor layout to be used, or 0,0 to automatically set the io_layout "//&
                 "to be the same as the layout.", default=1, layoutParam=.true.)

  call create_MOM_domain(MOM_dom, n_global, n_halo, reentrant, tripolar_N, layout, &
                         io_layout=io_layout, domain_name=domain_name, mask_table=mask_table, &
                         symmetric=symmetric, thin_halos=thin_halos, nonblocking=nonblocking)

end subroutine MOM_domains_init

!> Given a global array size and a number of (logical) processors, provide a layout of the
!! processors in the two directions where the total number of processors is the product of
!! the two layouts and number of points in the partitioned arrays are as close as possible
!! to an aspect ratio of 1.
subroutine MOM_define_layout(n_global, ndivs, layout)
  integer, dimension(2), intent(in)  :: n_global !< The total number of gridpoints in 2 directions
  integer,               intent(in)  :: ndivs    !< The total number of (logical) PEs
  integer, dimension(2), intent(out) :: layout   !< The generated layout of PEs

  ! Local variables
  integer :: isz, jsz, idiv, jdiv

  ! At present, this algorithm is a copy of mpp_define_layout, but it could perhaps be improved?

  isz = n_global(1) ; jsz = n_global(2)
  ! First try to divide ndivs to match the domain aspect ratio.  If this is not an even
  ! divisor of ndivs, reduce idiv until a factor is found.
  idiv = max(nint( sqrt(float(ndivs*isz)/jsz) ), 1)
  do while( mod(ndivs,idiv) /= 0 )
    idiv = idiv - 1
  enddo ! This will terminate at idiv=1 if not before
  jdiv = ndivs / idiv

  layout = (/ idiv, jdiv /)
end subroutine MOM_define_layout

end module MOM_domains
