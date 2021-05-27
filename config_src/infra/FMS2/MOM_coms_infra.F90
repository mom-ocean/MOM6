!> Thin interfaces to non-domain-oriented mpp communication subroutines
module MOM_coms_infra

! This file is part of MOM6. See LICENSE.md for the license.

use iso_fortran_env, only : int32, int64

use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_npes, mpp_set_root_pe
use mpp_mod, only : mpp_set_current_pelist, mpp_get_current_pelist
use mpp_mod, only : mpp_broadcast, mpp_sync, mpp_sync_self, mpp_chksum
use mpp_mod, only : mpp_sum, mpp_max, mpp_min
use memutils_mod, only : print_memuse_stats
use fms_mod, only : fms_end, fms_init

implicit none ; private

public :: PE_here, root_PE, num_PEs, set_rootPE, Set_PElist, Get_PElist
public :: broadcast, sum_across_PEs, min_across_PEs, max_across_PEs
public :: field_chksum, MOM_infra_init, MOM_infra_end

! This module provides interfaces to the non-domain-oriented communication
! subroutines.

!> Communicate an array, string or scalar from one PE to others
interface broadcast
  module procedure broadcast_char, broadcast_int32_0D, broadcast_int64_0D, broadcast_int1D
  module procedure broadcast_real0D, broadcast_real1D, broadcast_real2D
end interface broadcast

!> Compute a checksum for a field distributed over a PE list.  If no PE list is
!! provided, then the current active PE list is used.
interface field_chksum
  module procedure field_chksum_real_0d
  module procedure field_chksum_real_1d
  module procedure field_chksum_real_2d
  module procedure field_chksum_real_3d
  module procedure field_chksum_real_4d
end interface field_chksum

!> Find the sum of field across PEs, and update PEs with the sums.
interface sum_across_PEs
  module procedure sum_across_PEs_int4_0d
  module procedure sum_across_PEs_int4_1d
  module procedure sum_across_PEs_int8_0d
  module procedure sum_across_PEs_int8_1d
  module procedure sum_across_PEs_int8_2d
  module procedure sum_across_PEs_real_0d
  module procedure sum_across_PEs_real_1d
  module procedure sum_across_PEs_real_2d
end interface sum_across_PEs

!> Find the maximum value of field across PEs, and update PEs with the values.
interface max_across_PEs
  module procedure max_across_PEs_int_0d
  module procedure max_across_PEs_real_0d
  module procedure max_across_PEs_real_1d
end interface max_across_PEs

!> Find the minimum value of field across PEs, and update PEs with the values.
interface min_across_PEs
  module procedure min_across_PEs_int_0d
  module procedure min_across_PEs_real_0d
  module procedure min_across_PEs_real_1d
end interface min_across_PEs

contains

!> Return the ID of the PE for the current process.
function PE_here() result(pe)
  integer :: pe   !< PE ID of the current process
  pe = mpp_pe()
end function PE_here

!> Return the ID of the root PE for the PE list of the current procss.
function root_PE() result(pe)
  integer :: pe   !< root PE ID
  pe = mpp_root_pe()
end function root_PE

!> Return the number of PEs for the current PE list.
function num_PEs() result(npes)
  integer :: npes   !< Number of PEs
  npes = mpp_npes()
end function num_PEs

!> Designate a PE as the root PE
subroutine set_rootPE(pe)
  integer, intent(in) :: pe   !< ID of the PE to be assigned as root
  call mpp_set_root_pe(pe)
end subroutine

!> Set the current PE list.  If no list is provided, then the current PE list
!! is set to the list of all available PEs on the communicator.  Setting the
!! list will trigger a rank synchronization unless the `no_sync` flag is set.
subroutine Set_PEList(pelist, no_sync)
  integer, optional, intent(in) :: pelist(:)  !< List of PEs to set for communication
  logical, optional, intent(in) :: no_sync    !< Do not sync after list update.
  call mpp_set_current_pelist(pelist, no_sync)
end subroutine Set_PEList

!> Retrieve the current PE list and any metadata if requested.
subroutine Get_PEList(pelist, name, commID)
  integer,                    intent(out) :: pelist(:) !< List of PE IDs of the current PE list
  character(len=*), optional, intent(out) :: name   !< Name of PE list
  integer,          optional, intent(out) :: commID !< Communicator ID of PE list

  call mpp_get_current_pelist(pelist, name, commiD)
end subroutine Get_PEList

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
subroutine broadcast_int64_0D(dat, from_PE, PElist, blocking)
  integer(kind=int64),   intent(inout) :: dat       !< The data to communicate and destination
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

end subroutine broadcast_int64_0D


!> Communicate an integer from one PE to others
subroutine broadcast_int32_0D(dat, from_PE, PElist, blocking)
  integer(kind=int32),   intent(inout) :: dat       !< The data to communicate and destination
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

end subroutine broadcast_int32_0D

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

! field_chksum wrappers

!> Compute a checksum for a field distributed over a PE list.  If no PE list is
!! provided, then the current active PE list is used.
function field_chksum_real_0d(field, pelist, mask_val) result(chksum)
  real,              intent(in) :: field      !< Input scalar
  integer, optional, intent(in) :: pelist(:)  !< PE list of ranks to checksum
  real,    optional, intent(in) :: mask_val   !< FMS mask value
  integer(kind=int64) :: chksum               !< checksum of array

  chksum = mpp_chksum(field, pelist, mask_val)
end function field_chksum_real_0d

!> Compute a checksum for a field distributed over a PE list.  If no PE list is
!! provided, then the current active PE list is used.
function field_chksum_real_1d(field, pelist, mask_val) result(chksum)
  real, dimension(:), intent(in) :: field     !< Input array
  integer,  optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,     optional, intent(in) :: mask_val  !< FMS mask value
  integer(kind=int64) :: chksum               !< checksum of array

  chksum = mpp_chksum(field, pelist, mask_val)
end function field_chksum_real_1d

!> Compute a checksum for a field distributed over a PE list.  If no PE list is
!! provided, then the current active PE list is used.
function field_chksum_real_2d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:), intent(in) :: field     !< Unrotated input field
  integer,    optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,       optional, intent(in) :: mask_val  !< FMS mask value
  integer(kind=int64) :: chksum                 !< checksum of array

  chksum = mpp_chksum(field, pelist, mask_val)
end function field_chksum_real_2d

!> Compute a checksum for a field distributed over a PE list.  If no PE list is
!! provided, then the current active PE list is used.
function field_chksum_real_3d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:,:), intent(in) :: field     !< Unrotated input field
  integer,      optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,         optional, intent(in) :: mask_val  !< FMS mask value
  integer(kind=int64) :: chksum               !< checksum of array

  chksum = mpp_chksum(field, pelist, mask_val)
end function field_chksum_real_3d

!> Compute a checksum for a field distributed over a PE list.  If no PE list is
!! provided, then the current active PE list is used.
function field_chksum_real_4d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:,:,:), intent(in) :: field     !< Unrotated input field
  integer,        optional, intent(in) :: pelist(:) !< PE list of ranks to checksum
  real,           optional, intent(in) :: mask_val  !< FMS mask value
  integer(kind=int64) :: chksum               !< checksum of array

  chksum = mpp_chksum(field, pelist, mask_val)
end function field_chksum_real_4d

! sum_across_PEs wrappers

!> Find the sum of field across PEs, and return this sum in field.
subroutine sum_across_PEs_int4_0d(field, pelist)
  integer(kind=int32), intent(inout) :: field     !< Value on this PE, and the sum across PEs upon return
  integer,   optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_sum(field, pelist)
end subroutine sum_across_PEs_int4_0d

!> Find the sum of the values in corresponding positions of field across PEs, and return these sums in field.
subroutine sum_across_PEs_int4_1d(field, length, pelist)
  integer(kind=int32), dimension(:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,                           intent(in)    :: length    !< Number of elements in field to add
  integer,                 optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_sum(field, length, pelist)
end subroutine sum_across_PEs_int4_1d

!> Find the sum of field across PEs, and return this sum in field.
subroutine sum_across_PEs_int8_0d(field, pelist)
  integer(kind=int64), intent(inout) :: field     !< Value on this PE, and the sum across PEs upon return
  integer,   optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_sum(field, pelist)
end subroutine sum_across_PEs_int8_0d

!> Find the sum of the values in corresponding positions of field across PEs, and return these sums in field.
subroutine sum_across_PEs_int8_1d(field, length, pelist)
  integer(kind=int64), dimension(:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,                           intent(in)    :: length    !< Number of elements in field to add
  integer,                 optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_sum(field, length, pelist)
end subroutine sum_across_PEs_int8_1d

!> Find the sum of the values in corresponding positions of field across PEs, and return these sums in field.
subroutine sum_across_PEs_int8_2d(field, length, pelist)
  integer(kind=int64), &
           dimension(:,:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,                 intent(in)    :: length    !< The total number of positions to sum, usually
                                                      !! the product of the array sizes.
  integer,       optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_sum(field, length, pelist)
end subroutine sum_across_PEs_int8_2d

!> Find the sum of field across PEs, and return this sum in field.
subroutine sum_across_PEs_real_0d(field, pelist)
  real,              intent(inout) :: field     !< Value on this PE, and the sum across PEs upon return
  integer, optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_sum(field, pelist)
end subroutine sum_across_PEs_real_0d

!> Find the sum of the values in corresponding positions of field across PEs, and return these sums in field.
subroutine sum_across_PEs_real_1d(field, length, pelist)
  real, dimension(:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,            intent(in)    :: length    !< Number of elements in field to add
  integer,  optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_sum(field, length, pelist)
end subroutine sum_across_PEs_real_1d

!> Find the sum of the values in corresponding positions of field across PEs, and return these sums in field.
subroutine sum_across_PEs_real_2d(field, length, pelist)
  real, dimension(:,:), intent(inout) :: field     !< The values to add, the sums upon return
  integer,              intent(in)    :: length    !< The total number of positions to sum, usually
                                                   !! the product of the array sizes.
  integer,    optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_sum(field, length, pelist)
end subroutine sum_across_PEs_real_2d

! max_across_PEs wrappers

!> Find the maximum value of field across PEs, and store this maximum in field.
subroutine max_across_PEs_int_0d(field, pelist)
  integer,           intent(inout) :: field     !< The values to compare, the maximum upon return
  integer, optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_max(field, pelist)
end subroutine max_across_PEs_int_0d

!> Find the maximum value of field across PEs, and store this maximum in field.
subroutine max_across_PEs_real_0d(field, pelist)
  real,              intent(inout) :: field     !< The values to compare, the maximum upon return
  integer, optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_max(field, pelist)
end subroutine max_across_PEs_real_0d

!> Find the maximum values in each position of field across PEs, and store these minima in field.
subroutine max_across_PEs_real_1d(field, length, pelist)
  real, dimension(:), intent(inout) :: field     !< The list of values being compared, with the
                                                 !! maxima in each position upon return
  integer,            intent(in)    :: length    !< Number of elements in field to compare
  integer,  optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_max(field, length, pelist)
end subroutine max_across_PEs_real_1d

! min_across_PEs wrappers

!> Find the minimum value of field across PEs, and store this minimum in field.
subroutine min_across_PEs_int_0d(field, pelist)
  integer,           intent(inout) :: field     !< The values to compare, the minimum upon return
  integer, optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_min(field, pelist)
end subroutine min_across_PEs_int_0d

!> Find the minimum value of field across PEs, and store this minimum in field.
subroutine min_across_PEs_real_0d(field, pelist)
  real,              intent(inout) :: field     !< The values to compare, the minimum upon return
  integer, optional, intent(in)    :: pelist(:) !< List of PEs to work with
  call mpp_min(field, pelist)
end subroutine min_across_PEs_real_0d

!> Find the minimum values in each position of field across PEs, and store these minima in field.
subroutine min_across_PEs_real_1d(field, length, pelist)
  real, dimension(:), intent(inout) :: field     !< The list of values being compared, with the
                                                 !! minima in each position upon return
  integer,            intent(in)    :: length    !< Number of elements in field to compare
  integer,  optional, intent(in)    :: pelist(:) !< List of PEs to work with

  call mpp_min(field, length, pelist)
end subroutine min_across_PEs_real_1d

!> Initialize the model framework, including PE communication over a designated communicator.
!! If no communicator ID is provided, the framework's default communicator is used.
subroutine MOM_infra_init(localcomm)
  integer, optional, intent(in) :: localcomm  !< Communicator ID to initialize
  call fms_init(localcomm)
end subroutine

!> This subroutine carries out all of the calls required to close out the infrastructure cleanly.
!! This should only be called in ocean-only runs, as the coupler takes care of this in coupled runs.
subroutine MOM_infra_end
  call print_memuse_stats( 'Memory HiWaterMark', always=.TRUE. )
  call fms_end()
end subroutine MOM_infra_end

end module MOM_coms_infra
