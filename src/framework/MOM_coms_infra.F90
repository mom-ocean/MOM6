!> Thin interfaces to non-domain-oriented mpp communication subroutines
module MOM_coms_infra

! This file is part of MOM6. See LICENSE.md for the license.

use fms_mod, only : fms_end, MOM_infra_init => fms_init
use memutils_mod, only : print_memuse_stats
use mpp_mod, only : PE_here => mpp_pe, root_PE => mpp_root_pe, num_PEs => mpp_npes
use mpp_mod, only : set_rootPE => mpp_set_root_pe
use mpp_mod, only : Set_PElist => mpp_set_current_pelist, Get_PElist => mpp_get_current_pelist
use mpp_mod, only : mpp_broadcast, mpp_sync, mpp_sync_self, field_chksum => mpp_chksum
use mpp_mod, only : sum_across_PEs => mpp_sum, max_across_PEs => mpp_max, min_across_PEs => mpp_min

implicit none ; private

public :: PE_here, root_PE, num_PEs, MOM_infra_init, MOM_infra_end, Set_PElist, Get_PElist
public :: set_rootPE, broadcast, sum_across_PEs, min_across_PEs, max_across_PEs, field_chksum

! This module provides interfaces to the non-domain-oriented communication subroutines.

!> Communicate an array, string or scalar from one PE to others
interface broadcast
  module procedure broadcast_char, broadcast_int0D, broadcast_int1D
  module procedure broadcast_real0D, broadcast_real1D, broadcast_real2D
end interface broadcast

contains

!> Communicate a 1-D array of character strings from one PE to others
subroutine broadcast_char(dat, length, from_PE, PElist, blocking)
  character(len=*),  intent(inout) :: dat(:)    !< The data to communicate and destination
  integer,           intent(in)    :: length    !< The length of each string
  integer, optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer, optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                !! active PE set as previously set via Set_PElist.
  logical, optional, intent(in)    :: blocking  !< If true, barriers are added around the call

  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers

  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, length, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end subroutine broadcast_char

!> Communicate an integer from one PE to others
subroutine broadcast_int0D(dat, from_PE, PElist, blocking)
  integer,               intent(inout) :: dat       !< The data to communicate and destination
  integer,     optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,     optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                    !! active PE set as previously set via Set_PElist.
  logical,     optional, intent(in)    :: blocking  !< If true, barriers are added around the call

  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers

  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end subroutine broadcast_int0D

!> Communicate a 1-D array of integers from one PE to others
subroutine broadcast_int1D(dat, length, from_PE, PElist, blocking)
  integer, dimension(:), intent(inout) :: dat       !< The data to communicate and destination
  integer,               intent(in)    :: length    !< The number of data elements
  integer,     optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,     optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                    !! active PE set as previously set via Set_PElist.
  logical,     optional, intent(in)    :: blocking  !< If true, barriers are added around the call

  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers

  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, length, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end subroutine broadcast_int1D

!> Communicate a real number from one PE to others
subroutine broadcast_real0D(dat, from_PE, PElist, blocking)
  real,                 intent(inout) :: dat       !< The data to communicate and destination
  integer,    optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,    optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                   !! active PE set as previously set via Set_PElist.
  logical,    optional, intent(in)    :: blocking  !< If true, barriers are added around the call

  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers

  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end subroutine broadcast_real0D

!> Communicate a 1-D array of reals from one PE to others
subroutine broadcast_real1D(dat, length, from_PE, PElist, blocking)
  real, dimension(:),   intent(inout) :: dat       !< The data to communicate and destination
  integer,              intent(in)    :: length    !< The number of data elements
  integer,    optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,    optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                   !! active PE set as previously set via Set_PElist.
  logical,    optional, intent(in)    :: blocking  !< If true, barriers are added around the call

  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers

  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, length, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end subroutine broadcast_real1D

!> Communicate a 2-D array of reals from one PE to others
subroutine broadcast_real2D(dat, length, from_PE, PElist, blocking)
  real, dimension(:,:), intent(inout) :: dat       !< The data to communicate and destination
  integer,              intent(in)    :: length    !< The total number of data elements
  integer,    optional, intent(in)    :: from_PE   !< The source PE, by default the root PE
  integer,    optional, intent(in)    :: PElist(:) !< The list of participating PEs, by default the
                                                   !! active PE set as previously set via Set_PElist.
  logical,    optional, intent(in)    :: blocking  !< If true, barriers are added around the call

  integer :: src_PE   ! The processor that is sending the data
  logical :: do_block ! If true add synchronizing barriers

  do_block = .false. ; if (present(blocking)) do_block = blocking
  if (present(from_PE)) then ; src_PE = from_PE ; else ; src_PE = root_PE() ; endif

  if (do_block) call mpp_sync(PElist)
  call mpp_broadcast(dat, length, src_PE, PElist)
  if (do_block) call mpp_sync_self(PElist)

end subroutine broadcast_real2D


!> This subroutine carries out all of the calls required to close out the infrastructure cleanly.
!! This should only be called in ocean-only runs, as the coupler takes care of this in coupled runs.
subroutine MOM_infra_end
  call print_memuse_stats( 'Memory HiWaterMark', always=.TRUE. )
  call fms_end()
end subroutine MOM_infra_end

end module MOM_coms_infra
