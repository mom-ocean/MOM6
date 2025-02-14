!> Describes the decomposed MOM domain and has routines for communications across PEs
module MOM_domains

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms_infra,       only : MOM_infra_init, MOM_infra_end
use MOM_coms_infra,       only : PE_here, root_PE, num_PEs, broadcast
use MOM_coms_infra,       only : sum_across_PEs, min_across_PEs, max_across_PEs
use MOM_domain_infra,     only : MOM_domain_type, domain2D, domain1D, group_pass_type
use MOM_domain_infra,     only : create_MOM_domain, clone_MOM_domain, deallocate_MOM_domain
use MOM_domain_infra,     only : get_domain_extent, get_domain_components, same_domain
use MOM_domain_infra,     only : compute_block_extent, get_global_shape
use MOM_domain_infra,     only : pass_var, pass_vector, fill_symmetric_edges
use MOM_domain_infra,     only : pass_var_start, pass_var_complete
use MOM_domain_infra,     only : pass_vector_start, pass_vector_complete
use MOM_domain_infra,     only : create_group_pass, do_group_pass
use MOM_domain_infra,     only : start_group_pass, complete_group_pass
use MOM_domain_infra,     only : rescale_comp_data, global_field, redistribute_array, broadcast_domain
use MOM_domain_infra,     only : MOM_thread_affinity_set, set_MOM_thread_affinity
use MOM_domain_infra,     only : AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR
use MOM_domain_infra,     only : CORNER, CENTER, NORTH_FACE, EAST_FACE
use MOM_domain_infra,     only : To_East, To_West, To_North, To_South, To_All, Omit_Corners
use MOM_domain_infra,     only : compute_extent
use MOM_error_handler,    only : MOM_error, MOM_mesg, NOTE, WARNING, FATAL, is_root_pe
use MOM_file_parser,      only : get_param, log_param, log_version, param_file_type
use MOM_io_infra,         only : file_exists, read_field, open_ASCII_file, close_file, WRITEONLY_FILE
use MOM_string_functions, only : slasher
use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_unit_scaling,     only : unit_scale_type

implicit none ; private

public :: MOM_infra_init, MOM_infra_end
!  Domain types and creation and destruction routines
public :: MOM_domain_type, domain2D, domain1D
public :: MOM_domains_init, create_MOM_domain, clone_MOM_domain, deallocate_MOM_domain
public :: MOM_thread_affinity_set, set_MOM_thread_affinity
!  Domain query routines
public :: get_domain_extent, get_domain_components, get_global_shape, same_domain
public :: PE_here, root_PE, num_PEs
!  Blocks are not actively used in MOM6, so this routine could be deprecated.
public :: compute_block_extent
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

contains

!> MOM_domains_init initializes a MOM_domain_type variable, based on the information
!! read in from a param_file_type, and optionally returns data describing various
!! properties of the domain type.
subroutine MOM_domains_init(MOM_dom, param_file, symmetric, static_memory, &
                            NIHALO, NJHALO, NIGLOBAL, NJGLOBAL, NIPROC, NJPROC, &
                            min_halo, domain_name, include_name, param_suffix, US, MOM_dom_unmasked)
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
  type(unit_scale_type), optional, pointer       :: US           !< A dimensional unit scaling type
  type(MOM_domain_type), optional, pointer       :: MOM_dom_unmasked !< Unmasked MOM domain instance.
                                                                 !! Set to null if masking is not enabled.

  ! Local variables
  integer, dimension(2) :: layout    ! The number of logical processors in the i- and j- directions
  integer, dimension(2) :: auto_layout ! The layout determined by the auto masking routine
  integer, dimension(2) :: layout_unmasked ! A temporary layout for unmasked domain
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
  logical            :: auto_mask_table ! Runtime flag that turns on automatic mask table generator
  integer            :: auto_io_layout_fac ! Used to compute IO layout when auto_mask_table is True.
  logical            :: mask_table_exists ! True if there is a mask table file
  logical :: is_MOM_domain  ! True if this domain is being set for MOM, and not another component like SIS2.
  character(len=128) :: inputdir   ! The directory in which to find the diag table
  character(len=200) :: mask_table ! The file name and later the full path to the diag table
  character(len=64)  :: inc_nm     ! The name of the memory include file
  character(len=200) :: mesg       ! A string to use for error messages

  integer :: nip_parsed, njp_parsed
  character(len=8) :: char_xsiz, char_ysiz, char_niglobal, char_njglobal
  character(len=40) :: nihalo_nm, njhalo_nm, layout_nm, io_layout_nm, masktable_nm
  character(len=40) :: niproc_nm, njproc_nm
  character(len=200) :: topo_config
  integer :: id_clock_auto_mask
  character(len=:), allocatable :: masktable_desc
  character(len=:), allocatable :: auto_mask_table_fname ! Auto-generated mask table file name
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
  !$              default=1, layoutParam=.true.)
  !$   call get_param(param_file, mdl, "OCEAN_OMP_HYPER_THREAD", ocean_omp_hyper_thread, &
  !$              "If True, use hyper-threading.", default=.false., layoutParam=.true.)
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
                 default=NIGLOBAL)
    call get_param(param_file, mdl, "NJGLOBAL", n_global(2), &
                 "The total number of thickness grid points in the y-direction in the physical "//&
                 "domain. With STATIC_MEMORY_ this is set in "//trim(inc_nm)//" at compile time.", &
                 default=NJGLOBAL)
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
                 default=nihalo_dflt)
  call get_param(param_file, mdl, trim(njhalo_nm), n_halo(2), &
                 "The number of halo points on each side in the y-direction.  How this is set "//&
                 "varies with the calling component and static or dynamic memory configuration.", &
                 default=njhalo_dflt)
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

  is_MOM_domain = .true.
  if (present(domain_name)) then
    is_MOM_domain = (index(domain_name, "MOM") > 1)
  endif

  if (is_MOM_domain) then
    call get_param(param_file, mdl, "TOPO_CONFIG", topo_config, do_not_log=.true., fail_if_missing=.true.)
  else ! SIS2 has a default value for TOPO_CONFIG.
    call get_param(param_file, mdl, "TOPO_CONFIG", topo_config, default="file", do_not_log=.true.)
  endif

  auto_mask_table = .false.
  if (.not. present(param_suffix) .and. .not. is_static .and. trim(topo_config) == 'file') then
     call get_param(param_file, mdl, 'AUTO_MASKTABLE', auto_mask_table, &
                 "Turn on automatic mask table generation to eliminate land blocks.", &
                 default=.false., layoutParam=.true.)
  endif

  masktable_desc = "A text file to specify n_mask, layout and mask_list. This feature masks out "//&
      "processors that contain only land points. The first line of mask_table is the "//&
      "number of regions to be masked out. The second line is the layout of the "//&
      "model and must be consistent with the actual model layout. The following "//&
      "(n_mask) lines give the logical positions of the processors that are masked "//&
      "out. The mask_table can be created by tools like check_mask. The following "//&
      "example of mask_table masks out 2 processors, (1,2) and (3,6), out of the 24 "//&
      "in a 4x6 layout: \n 2\n 4,6\n 1,2\n 3,6\n"

  if (auto_mask_table) then
    id_clock_auto_mask = cpu_clock_id('(Ocean gen_auto_mask_table)', grain=CLOCK_ROUTINE)
    auto_mask_table_fname = "MOM_auto_mask_table"

    ! Auto-generate a mask file and determine the layout
    call cpu_clock_begin(id_clock_auto_mask)
    if (is_root_PE()) then
      call gen_auto_mask_table(n_global, reentrant, tripolar_N, PEs_used, param_file, inputdir, &
                               auto_mask_table_fname, auto_layout, US)
    endif
    call broadcast(auto_layout, length=2)
    call cpu_clock_end(id_clock_auto_mask)

    mask_table = auto_mask_table_fname
    call log_param(param_file, mdl, trim(masktable_nm), mask_table, masktable_desc, &
                   default="MOM_mask_table", layoutParam=.true.)
  else
    call get_param(param_file, mdl, trim(masktable_nm), mask_table, masktable_desc, &
                   default="MOM_mask_table", layoutParam=.true.)
  endif

  ! First, check the run directory for the mask_table input file.
  mask_table_exists = file_exists(trim(mask_table))
  ! If not found, check the input directory
  if (.not. mask_table_exists) then
    mask_table = trim(inputdir)//trim(mask_table)
    mask_table_exists = file_exists(mask_table)
  endif

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

    if (auto_mask_table) then
        if (layout(1) /= 0 .and. layout(1) /= auto_layout(1)) then
          call MOM_error(FATAL, "Cannot set LAYOUT or NIPROC when AUTO_MASKTABLE is enabled.")
        endif
        if (layout(2) /= 0 .and. layout(2) /= auto_layout(2)) then
          call MOM_error(FATAL, "Cannot set LAYOUT or NJPROC when AUTO_MASKTABLE is enabled.")
        endif
        layout(:) = auto_layout(:)
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
    write(mesg,'(a,2(i5,1x,a))') 'You requested to use', layout(1)*layout(2), &
      'PEs but there are only', n_global(1)*n_global(2), 'columns in the model'
    call MOM_error(FATAL, mesg)
  endif

  if (mask_table_exists) &
    call MOM_error(NOTE, 'MOM_domains_init: reading maskmap information from '//trim(mask_table))

  ! Set up the I/O layout, it will be checked later that it uses an even multiple of the number of
  ! PEs in each direction.
  io_layout(:) = (/ 1, 1 /)

  ! Compute a valid IO layout if auto_mask_table is on. Otherwise, read in IO_LAYOUT parameter,
  if (auto_mask_table) then
    call get_param(param_file, mdl, "AUTO_IO_LAYOUT_FAC", auto_io_layout_fac, &
            "When AUTO_MASKTABLE is enabled, io layout is calculated by performing integer "//&
            "division of the runtime-determined domain layout with this factor. If the factor "//&
            "is set to 0 (default), the io layout is set to 1,1.", &
            default=0, layoutParam=.true.)
    if (auto_io_layout_fac>0) then
      io_layout(1) = max(layout(1)/auto_io_layout_fac, 1)
      io_layout(2) = max(layout(2)/auto_io_layout_fac, 1)
    elseif (auto_io_layout_fac<0) then
      call MOM_error(FATAL, 'AUTO_IO_LAYOUT_FAC must be a nonnegative integer.')
    endif
    call log_param(param_file, mdl, trim(io_layout_nm), io_layout, &
                   "The processor layout to be used, or 0,0 to automatically set the io_layout "//&
                   "to be the same as the layout.", layoutParam=.true.)
  else
    call get_param(param_file, mdl, trim(io_layout_nm), io_layout, &
                   "The processor layout to be used, or 0,0 to automatically set the io_layout "//&
                   "to be the same as the layout.", default=1, layoutParam=.true.)
  endif

  ! Create an unmasked domain if requested. This is used for writing out unmasked ocean geometry.
  if (present(MOM_dom_unmasked) .and. mask_table_exists) then
    call MOM_define_layout(n_global, PEs_used, layout_unmasked)
    call create_MOM_domain(MOM_dom_unmasked, n_global, n_halo, reentrant, tripolar_N, layout_unmasked, &
                           domain_name=domain_name, symmetric=symmetric, thin_halos=thin_halos, &
                           nonblocking=nonblocking)
  endif

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

!> Given a desired number of active npes, generate a layout and mask_table
subroutine gen_auto_mask_table(n_global, reentrant, tripolar_N, npes, param_file, inputdir, filename, layout, US)
  integer, dimension(2), intent(in)         :: n_global   !< The total number of gridpoints in 2 directions
  logical, dimension(2), intent(in)         :: reentrant  !< True if the x- and y- directions are periodic.
  logical,               intent(in)         :: tripolar_N !< A flag indicating whether there is n. tripolar connectivity
  integer,               intent(in)         :: npes       !< The desired number of active PEs.
  type(param_file_type), intent(in)         :: param_file !< A structure to parse for run-time parameters
  character(len=128),    intent(in)         :: inputdir   !< INPUTDIR parameter
  character(len=:), allocatable, intent(in) :: filename   !< Mask table file path (to be auto-generated.)
  integer, dimension(2), intent(out)        :: layout     !< The generated layout of PEs (incl. masked blocks)
  type(unit_scale_type), optional, pointer  :: US         !< A dimensional unit scaling type

  ! Local variables
  real, dimension(n_global(1), n_global(2)) :: D        ! Bathymetric depth (to be read in from TOPO_FILE) [Z ~> m]
  integer, dimension(:,:), allocatable :: mask          ! Cell masks (based on D and MINIMUM_DEPTH)
  character(len=200) :: topo_filepath, topo_file        ! Strings for file/path
  character(len=200) :: topo_varname                    ! Variable name in file
  character(len=200) :: topo_config
  character(len=40)  :: mdl = "gen_auto_mask_table"      ! This subroutine's name.
  integer :: i, j, p
  real :: Dmask          ! The depth for masking in the same units as D             [Z ~> m]
  real :: min_depth      ! The minimum ocean depth in the same units as D           [Z ~> m]
  real :: mask_depth     ! The depth shallower than which to mask a point as land.  [Z ~> m]
  real :: glob_ocn_frac  ! ratio of ocean points to total number of points          [nondim]
  real :: r_p            ! aspect ratio for division count p.                       [nondim]
  real :: m_to_Z         ! A conversion factor from m to height units           [Z m-1 ~> 1]
  integer :: nx, ny      ! global domain sizes
  integer, parameter :: ibuf=2, jbuf=2
  real, parameter :: r_extreme = 4.0 ! aspect ratio limit (>1) for a layout to be considered [nondim]
  integer :: num_masked_blocks
  integer, allocatable :: mask_table(:,:)

  m_to_Z = 1.0 ; if (present(US)) m_to_Z = US%m_to_Z

  ! Read in params necessary for auto-masking
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 units="m", default=0.0, scale=m_to_Z, do_not_log=.true.)
  call get_param(param_file, mdl, "MASKING_DEPTH", mask_depth, &
                 units="m", default=-9999.0, scale=m_to_Z, do_not_log=.true.)
  call get_param(param_file, mdl, "TOPO_CONFIG", topo_config, default="file", do_not_log=.true.)
  call get_param(param_file, mdl, "TOPO_FILE", topo_file, do_not_log=.true., default="topog.nc")
  call get_param(param_file, mdl, "TOPO_VARNAME", topo_varname, do_not_log=.true., default="depth")
  topo_filepath = trim(inputdir)//trim(topo_file)

  ! Sanity checks
  if (.not. is_root_pe()) then
    call MOM_error(FATAL, 'gen_auto_mask_table should only be called by the root PE.')
  endif
  if (trim(topo_config) /= "file") then
    call MOM_error(FATAL, 'Auto mask table only works with TOPO_CONFIG="file"')
  endif
  if (.not.file_exists(topo_filepath)) then
    call MOM_error(FATAL, " gen_auto_mask_table: Unable to open "//trim(topo_filepath))
  endif

  nx = n_global(1)
  ny = n_global(2)

  ! Read in bathymetric depth.
  D(:,:) = -9.0e30 * m_to_Z ! Initializing to a very large negative depth (tall mountains) everywhere.
  call read_field(topo_filepath, trim(topo_varname), D, start=(/1, 1/), nread=n_global, no_domain=.true., &
                  scale=m_to_Z)

  allocate(mask(nx+2*ibuf, ny+2*jbuf), source=0)

  ! Determine cell masks
  Dmask = mask_depth
  if (mask_depth == -9999.0*m_to_Z) Dmask = min_depth
  do i=1,nx ; do j=1,ny
    if (D(i,j) <= Dmask) then
      mask(i+ibuf,j+jbuf) = 0
    else
      mask(i+ibuf,j+jbuf) = 1
    endif
  enddo ; enddo

  ! fill in buffer cells

  if (reentrant(1)) then ! REENTRANT_X
    mask(1:ibuf, :) = mask(nx+1:nx+ibuf, :)
    mask(ibuf+nx+1:nx+2*ibuf, :) = mask(ibuf+1:2*ibuf, :)
  endif

  if (reentrant(2)) then ! REENTRANT_Y
    mask(:, 1:jbuf) = mask(:, ny+1:ny+jbuf)
    mask(:, jbuf+ny+1:ny+2*jbuf) = mask(:, jbuf+1:2*jbuf)
  endif

  if (tripolar_N) then ! TRIPOLAR_N
    do i=1,nx+2*ibuf
      do j=1,jbuf
        mask(i, jbuf+ny+j) = mask(nx+2*ibuf+1-i, jbuf+ny+1-j)
      enddo
    enddo
  endif

  ! Tripolar Stitch Fix: In cases where masking is asymmetrical across the tripolar stitch, there's a possibility
  ! that certain unmasked blocks won't be able to obtain grid metrics from the halo points. This occurs when the
  ! neighboring block on the opposite side of the tripolar stitch is masked. As a consequence, certain metrics like
  ! dxT and dyT may be calculated through extrapolation (refer to extrapolate_metric), potentially leading to the
  ! generation of non-positive values. This can result in divide-by-zero errors elsewhere, e.g., in MOM_hor_visc.F90.
  ! Currently, the safest and most general solution is to prohibit masking along the tripolar stitch:
  if (tripolar_N) then
    mask(:, jbuf+ny) = 1
  endif

  glob_ocn_frac = real(sum(mask(1+ibuf:nx+ibuf, 1+jbuf:ny+jbuf))) / (nx * ny)

  ! Iteratively check for all possible division counts starting from the upper bound of npes/glob_ocn_frac,
  ! which is over-optimistic for realistic domains, but may be satisfied with idealized domains.
  do p = ceiling(npes/glob_ocn_frac), npes, -1

    ! compute the layout for the current division count, p
    call MOM_define_layout(n_global, p, layout)

    ! don't bother checking this p if the aspect ratio is extreme
    r_p = (real(nx)/layout(1)) / (real(ny)/layout(2))
    if ( r_p * r_extreme < 1 .or. r_extreme < r_p ) cycle

    ! Get the number of masked_blocks for this particular division count
    call determine_land_blocks(mask, nx, ny, layout(1), layout(2), ibuf, jbuf, num_masked_blocks)

    ! If we can eliminate enough blocks to reach the target npes, adopt
    ! this p (and the associated layout) and terminate the iteration.
    if (p-num_masked_blocks <= npes) then
      call MOM_error(NOTE, "Found the optimum layout for auto-masking. Terminating iteration...")
      exit
    endif
  enddo

  if (num_masked_blocks == 0) then
    call MOM_error(FATAL, "Couldn't auto-eliminate any land blocks. Try to increase the number "//&
        "of MOM6 PEs or set AUTO_MASKTABLE to False.")
  endif

  ! Call determine_land_blocks once again, this time to retrieve and write out the mask_table.
  allocate(mask_table(num_masked_blocks,2))
  call determine_land_blocks(mask, nx, ny, layout(1), layout(2), ibuf, jbuf, num_masked_blocks, mask_table)
  call write_auto_mask_file(mask_table, layout, npes, filename)
  deallocate(mask_table)
  deallocate(mask)

end subroutine gen_auto_mask_table

!> Given a number of domain divisions, compute the max number of land blocks that can be eliminated,
!! and return the resulting mask table if requested.
subroutine determine_land_blocks(mask, nx, ny, idiv, jdiv, ibuf, jbuf, num_masked_blocks, mask_table)
  integer, dimension(:,:), intent(in)   :: mask     !< cell masks based on depth and MINIMUM_DEPTH
  integer, intent(in)                   :: nx       !< Total number of gridpoints in x-dir (global)
  integer, intent(in)                   :: ny       !< Total number of gridpoints in y-dir (global)
  integer, intent(in)                   :: idiv     !< number of divisions along x-dir
  integer, intent(in)                   :: jdiv     !< number of divisions along y-dir
  integer, intent(in)                   :: ibuf     !< number of buffer cells in x-dir.
                                                    !! (not necessarily the same as NIHALO)
  integer, intent(in)                   :: jbuf     !< number of buffer cells in y-dir.
                                                    !! (not necessarily the same as NJHALO)
  integer, intent(out)                  :: num_masked_blocks !< the final number of masked blocks
  integer, intent(out), optional        :: mask_table(:,:) !< the resulting array of mask_table
  ! integer
  integer, dimension(idiv) :: ibegin   !< The starting index of each division along x axis
  integer, dimension(idiv) :: iend     !< The ending index of each division along x axis
  integer, dimension(jdiv) :: jbegin   !< The starting index of each division along y axis
  integer, dimension(jdiv) :: jend     !< The ending index of each division along y axis
  integer :: i, j, ib, ie, jb,je

  call compute_extent(1, nx, idiv, ibegin, iend)
  call compute_extent(1, ny, jdiv, jbegin, jend)

  num_masked_blocks = 0

  do i=1,idiv
    ib = ibegin(i)
    ie = iend(i) + 2 * ibuf
    do j=1,jdiv
      jb = jbegin(j)
      je = jend(j) + 2 * jbuf

      if (any(mask(ib:ie,jb:je)==1)) cycle

      num_masked_blocks = num_masked_blocks + 1

      if (present(mask_table)) then
        if ( num_masked_blocks > size(mask_table, dim=1)) then
          call MOM_error(FATAL, "The mask_table argument passed to determine_land_blocks() has insufficient size.")
        endif

        mask_table(num_masked_blocks,1) = i
        mask_table(num_masked_blocks,2) = j
      endif
    enddo
  enddo

end subroutine determine_land_blocks

!> Write out the auto-generated mask information to a file in the run directory.
subroutine write_auto_mask_file(mask_table, layout, npes, filename)
  integer, intent(in) :: mask_table(:,:)      !> mask table array to be written out.
  integer, dimension(2), intent(in) :: layout !> PE layout
  integer, intent(in) :: npes                 !> Number of divisions (incl. eliminated ones)
  character(len=:), allocatable, intent(in) :: filename !> file name for the mask_table to be written
  ! local
  integer :: file_ascii= -1  !< The unit number of the auto-generated mask_file file.
  integer :: true_num_masked_blocks
  integer :: p

  ! Eliminate only enough blocks to ensure that the number of active blocks precisely matches the target npes.
  true_num_masked_blocks = layout(1) * layout(2) - npes

  call open_ASCII_file(file_ascii, trim(filename), action=WRITEONLY_FILE)
  write(file_ascii, '(I0)') true_num_masked_blocks
  write(file_ascii, '(I0,",",I0)') layout(1), layout(2)
  do p = 1, true_num_masked_blocks
    write(file_ascii, '(I0,",",I0)') mask_table(p,1), mask_table(p,2)
  enddo
  call close_file(file_ascii)
end subroutine write_auto_mask_file

end module MOM_domains
