!> Describes the decomposed MOM domain and has routines for communications across PEs
module MOM_domain_infra

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms_infra,  only : PE_here, root_PE, num_PEs
use MOM_cpu_clock_infra, only : cpu_clock_begin, cpu_clock_end
use MOM_error_infra, only : MOM_error=>MOM_err, NOTE, WARNING, FATAL

use mpp_domains_mod, only : domain2D, domain1D
use mpp_domains_mod, only : mpp_define_io_domain, mpp_define_domains, mpp_deallocate_domain
use mpp_domains_mod, only : mpp_get_domain_components, mpp_get_domain_extents
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
use mpp_domains_mod, only : mpp_get_boundary, mpp_update_domains
use mpp_domains_mod, only : mpp_start_update_domains, mpp_complete_update_domains
use mpp_domains_mod, only : mpp_create_group_update, mpp_do_group_update
use mpp_domains_mod, only : mpp_reset_group_update_field, mpp_group_update_initialized
use mpp_domains_mod, only : mpp_start_group_update, mpp_complete_group_update
use mpp_domains_mod, only : mpp_compute_block_extent
use mpp_domains_mod, only : mpp_broadcast_domain, mpp_redistribute, mpp_global_field
use mpp_domains_mod, only : AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR, BITWISE_EXACT_SUM
use mpp_domains_mod, only : CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
use mpp_domains_mod, only : To_East => WUPDATE, To_West => EUPDATE, Omit_Corners => EDGEUPDATE
use mpp_domains_mod, only : To_North => SUPDATE, To_South => NUPDATE
use mpp_domains_mod, only : CENTER, CORNER, NORTH_FACE => NORTH, EAST_FACE => EAST
use fms_io_mod,      only : file_exist, parse_mask_table
use fms_affinity_mod, only : fms_affinity_init, fms_affinity_set, fms_affinity_get

! This subroutine is not in MOM6/src but may be required by legacy drivers
use mpp_domains_mod, only : global_field_sum => mpp_global_sum

! The `group_pass_type` fields are never accessed, so we keep it as an FMS type
use mpp_domains_mod, only : group_pass_type => mpp_group_update_type

implicit none ; private

! These types are inherited from mpp, but are treated as opaque here.
public :: domain2D, domain1D, group_pass_type
! These interfaces are actually implemented or have explicit interfaces in this file.
public :: create_MOM_domain, clone_MOM_domain, get_domain_components, get_domain_extent
public :: deallocate_MOM_domain, get_global_shape, compute_block_extent
public :: pass_var, pass_vector, fill_symmetric_edges, rescale_comp_data
public :: pass_var_start, pass_var_complete, pass_vector_start, pass_vector_complete
public :: create_group_pass, do_group_pass, start_group_pass, complete_group_pass
public :: redistribute_array, broadcast_domain, global_field
public :: get_simple_array_i_ind, get_simple_array_j_ind
public :: MOM_thread_affinity_set, set_MOM_thread_affinity
! These are encoding constant parmeters.
public :: To_East, To_West, To_North, To_South, To_All, Omit_Corners
public :: AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR
public :: CORNER, CENTER, NORTH_FACE, EAST_FACE
! These are no longer used by MOM6 because the reproducing sum works so well, but they are
! still referenced by some of the non-GFDL couplers.
public :: global_field_sum, BITWISE_EXACT_SUM

!> Do a halo update on an array
interface pass_var
  module procedure pass_var_3d, pass_var_2d
end interface pass_var

!> Do a halo update on a pair of arrays representing the two components of a vector
interface pass_vector
  module procedure pass_vector_3d, pass_vector_2d
end interface pass_vector

!> Initiate a non-blocking halo update on an array
interface pass_var_start
  module procedure pass_var_start_3d, pass_var_start_2d
end interface pass_var_start

!> Complete a non-blocking halo update on an array
interface pass_var_complete
  module procedure pass_var_complete_3d, pass_var_complete_2d
end interface pass_var_complete

!> Initiate a halo update on a pair of arrays representing the two components of a vector
interface pass_vector_start
  module procedure pass_vector_start_3d, pass_vector_start_2d
end interface pass_vector_start

!> Complete a halo update on a pair of arrays representing the two components of a vector
interface pass_vector_complete
  module procedure pass_vector_complete_3d, pass_vector_complete_2d
end interface pass_vector_complete

!> Set up a group of halo updates
interface create_group_pass
  module procedure create_var_group_pass_2d
  module procedure create_var_group_pass_3d
  module procedure create_vector_group_pass_2d
  module procedure create_vector_group_pass_3d
end interface create_group_pass

!> Do a set of halo updates that fill in the values at the duplicated edges
!! of a staggered symmetric memory domain
interface fill_symmetric_edges
  module procedure fill_vector_symmetric_edges_2d !, fill_vector_symmetric_edges_3d
!   module procedure fill_scalar_symmetric_edges_2d, fill_scalar_symmetric_edges_3d
end interface fill_symmetric_edges

!> Rescale the values of an array in its computational domain by a constant factor
interface rescale_comp_data
  module procedure rescale_comp_data_4d, rescale_comp_data_3d, rescale_comp_data_2d
end interface rescale_comp_data

!> Pass an array from one MOM domain to another
interface redistribute_array
  module procedure redistribute_array_3d, redistribute_array_2d
end interface redistribute_array

!> Copy one MOM_domain_type into another
interface clone_MOM_domain
  module procedure clone_MD_to_MD, clone_MD_to_d2D
end interface clone_MOM_domain

!> Extract the 1-d domain components from a MOM_domain or domain2d
interface get_domain_components
  module procedure get_domain_components_MD, get_domain_components_d2D
end interface get_domain_components

!> Returns the index ranges that have been stored in a MOM_domain_type
interface get_domain_extent
  module procedure get_domain_extent_MD, get_domain_extent_d2D
end interface get_domain_extent


!> The MOM_domain_type contains information about the domain decomposition.
type, public :: MOM_domain_type
  character(len=64) :: name     !< The name of this domain
  type(domain2D), pointer :: mpp_domain => NULL() !< The FMS domain with halos
                                !! on this processor, centered at h points.
  type(domain2D), pointer :: mpp_domain_d2 => NULL() !< A coarse FMS domain with halos
                                !! on this processor, centered at h points.
  integer :: niglobal           !< The total horizontal i-domain size.
  integer :: njglobal           !< The total horizontal j-domain size.
  integer :: nihalo             !< The i-halo size in memory.
  integer :: njhalo             !< The j-halo size in memory.
  logical :: symmetric          !< True if symmetric memory is used with this domain.
  logical :: nonblocking_updates  !< If true, non-blocking halo updates are
                                !! allowed.  The default is .false. (for now).
  logical :: thin_halo_updates  !< If true, optional arguments may be used to
                                !! specify the width of the halos that are
                                !! updated with each call.
  integer :: layout(2)          !< This domain's processor layout.  This is
                                !! saved to enable the construction of related
                                !! new domains with different resolutions or
                                !! other properties.
  integer :: io_layout(2)       !< The IO-layout used with this domain.
  integer :: X_FLAGS            !< Flag that specifies the properties of the
                                !! domain in the i-direction in a define_domain call.
  integer :: Y_FLAGS            !< Flag that specifies the properties of the
                                !! domain in the j-direction in a define_domain call.
  logical, pointer :: maskmap(:,:) => NULL() !< A pointer to an array indicating
                                !! which logical processors are actually used for
                                !! the ocean code. The other logical processors
                                !! would be contain only land points and are not
                                !! assigned to actual processors. This need not be
                                !! assigned if all logical processors are used.
end type MOM_domain_type

integer, parameter :: To_All = To_East + To_West + To_North + To_South !< A flag for passing in all directions

contains

!> pass_var_3d does a halo update for a three-dimensional array.
subroutine pass_var_3d(array, MOM_dom, sideflag, complete, position, halo, &
                       clock)
  real, dimension(:,:,:), intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! sothe halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                                    !! halo updates should be completed before
                                                    !! progress resumes. Omitting complete is the
                                                    !! same as setting complete to .true.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.

  integer :: dirflag
  logical :: block_til_complete

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif
  block_til_complete = .true.
  if (present(complete)) block_til_complete = complete

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                        complete=block_til_complete, position=position, &
                        whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                          complete=block_til_complete, position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine pass_var_3d

!> pass_var_2d does a halo update for a two-dimensional array.
subroutine pass_var_2d(array, MOM_dom, sideflag, complete, position, halo, inner_halo, clock)
  real, dimension(:,:),  intent(inout) :: array    !< The array which is having its halos points
                                                   !! exchanged.
  type(MOM_domain_type), intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                   !! needed to determine where data should be sent.
  integer,     optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  logical,     optional, intent(in)    :: complete !< An optional argument indicating whether the
                                                   !! halo updates should be completed before
                                                   !! progress resumes.  Omitting complete is the
                                                   !! same as setting complete to .true.
  integer,     optional, intent(in)    :: position !< An optional argument indicating the position.
                                                   !! This is CENTER by default and is often CORNER,
                                                   !! but could also be EAST_FACE or NORTH_FACE.
  integer,     optional, intent(in)    :: halo     !< The size of the halo to update - the full halo
                                                   !! by default.
  integer,     optional, intent(in)    :: inner_halo !< The size of an inner halo to avoid updating,
                                                   !! or 0 to avoid updating symmetric memory
                                                   !! computational domain points.  Setting this >=0
                                                   !! also enforces that complete=.true.
  integer,     optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                   !! started then stopped to time this routine.

  ! Local variables
  real, allocatable, dimension(:,:) :: tmp
  integer :: pos, i_halo, j_halo
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, IscB, IecB, JscB, JecB
  integer :: inner, i, j, isfw, iefw, isfe, iefe, jsfs, jefs, jsfn, jefn
  integer :: dirflag
  logical :: block_til_complete

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif
  block_til_complete = .true. ; if (present(complete)) block_til_complete = complete
  pos = CENTER ; if (present(position)) pos = position

  if (present(inner_halo)) then ; if (inner_halo >= 0) then
    ! Store the original values.
    allocate(tmp(size(array,1), size(array,2)))
    tmp(:,:) = array(:,:)
    block_til_complete = .true.
  endif ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                        complete=block_til_complete, position=position, &
                        whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                        complete=block_til_complete, position=position)
  endif

  if (present(inner_halo)) then ; if (inner_halo >= 0) then
    call mpp_get_compute_domain(MOM_dom%mpp_domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(MOM_dom%mpp_domain, isd, ied, jsd, jed)
    ! Convert to local indices for arrays starting at 1.
    isc = isc - (isd-1) ; iec = iec - (isd-1) ; ied = ied - (isd-1) ; isd = 1
    jsc = jsc - (jsd-1) ; jec = jec - (jsd-1) ; jed = jed - (jsd-1) ; jsd = 1
    i_halo = min(inner_halo, isc-1) ; j_halo = min(inner_halo, jsc-1)

    ! Figure out the array index extents of the eastern, western, northern and southern regions to copy.
    if (pos == CENTER) then
      if (size(array,1) == ied) then
        isfw = isc - i_halo ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      else ; call MOM_error(FATAL, "pass_var_2d: wrong i-size for CENTER array.") ; endif
      if (size(array,2) == jed) then
        isfw = isc - i_halo ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      else ; call MOM_error(FATAL, "pass_var_2d: wrong j-size for CENTER array.") ; endif
    elseif (pos == CORNER) then
      if (size(array,1) == ied) then
        isfw = max(isc - (i_halo+1), 1) ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      elseif (size(array,1) == ied+1) then
        isfw = isc - i_halo ; iefw = isc+1 ; isfe = iec+1 ; iefe = min(iec + 1 + i_halo, ied+1)
      else ; call MOM_error(FATAL, "pass_var_2d: wrong i-size for CORNER array.") ; endif
      if (size(array,2) == jed) then
        jsfs = max(jsc - (j_halo+1), 1) ; jefs = jsc ; jsfn = jec ; jefn = jec + j_halo
      elseif (size(array,2) == jed+1) then
        jsfs = jsc - j_halo ; jefs = jsc+1 ; jsfn = jec+1 ; jefn = min(jec + 1 + j_halo, jed+1)
      else ; call MOM_error(FATAL, "pass_var_2d: wrong j-size for CORNER array.") ; endif
    elseif (pos == NORTH_FACE) then
      if (size(array,1) == ied) then
        isfw = isc - i_halo ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      else ; call MOM_error(FATAL, "pass_var_2d: wrong i-size for NORTH_FACE array.") ; endif
      if (size(array,2) == jed) then
        jsfs = max(jsc - (j_halo+1), 1) ; jefs = jsc ; jsfn = jec ; jefn = jec + j_halo
      elseif (size(array,2) == jed+1) then
        jsfs = jsc - j_halo ; jefs = jsc+1 ; jsfn = jec+1 ; jefn = min(jec + 1 + j_halo, jed+1)
      else ; call MOM_error(FATAL, "pass_var_2d: wrong j-size for NORTH_FACE array.") ; endif
    elseif (pos == EAST_FACE) then
      if (size(array,1) == ied) then
        isfw = max(isc - (i_halo+1), 1) ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      elseif (size(array,1) == ied+1) then
        isfw = isc - i_halo ; iefw = isc+1 ; isfe = iec+1 ; iefe = min(iec + 1 + i_halo, ied+1)
      else ; call MOM_error(FATAL, "pass_var_2d: wrong i-size for EAST_FACE array.") ; endif
      if (size(array,2) == jed) then
        isfw = isc - i_halo ; iefw = isc ; isfe = iec ; iefe = iec + i_halo
      else ; call MOM_error(FATAL, "pass_var_2d: wrong j-size for EAST_FACE array.") ; endif
    else
      call MOM_error(FATAL, "pass_var_2d: Unrecognized position")
    endif

    ! Copy back the stored inner halo points
    do j=jsfs,jefn ; do i=isfw,iefw ; array(i,j) = tmp(i,j) ; enddo ; enddo
    do j=jsfs,jefn ; do i=isfe,iefe ; array(i,j) = tmp(i,j) ; enddo ; enddo
    do j=jsfs,jefs ; do i=isfw,iefe ; array(i,j) = tmp(i,j) ; enddo ; enddo
    do j=jsfn,jefn ; do i=isfw,iefe ; array(i,j) = tmp(i,j) ; enddo ; enddo

    deallocate(tmp)
  endif ; endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine pass_var_2d

!> pass_var_start_2d starts a halo update for a two-dimensional array.
function pass_var_start_2d(array, MOM_dom, sideflag, position, complete, halo, &
                           clock)
  real, dimension(:,:),   intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                                    !! halo updates should be completed before
                                                    !! progress resumes.  Omitting complete is the
                                                    !! same as setting complete to .true.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  integer                               :: pass_var_start_2d  !<The integer index for this update.

  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    pass_var_start_2d = mpp_start_update_domains(array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position, &
                            whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    pass_var_start_2d = mpp_start_update_domains(array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end function pass_var_start_2d

!> pass_var_start_3d starts a halo update for a three-dimensional array.
function pass_var_start_3d(array, MOM_dom, sideflag, position, complete, halo, &
                           clock)
  real, dimension(:,:,:), intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                                    !! halo updates should be completed before
                                                    !! progress resumes.  Omitting complete is the
                                                    !! same as setting complete to .true.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  integer                               :: pass_var_start_3d  !< The integer index for this update.

  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    pass_var_start_3d = mpp_start_update_domains(array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position, &
                            whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    pass_var_start_3d = mpp_start_update_domains(array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end function pass_var_start_3d

!> pass_var_complete_2d completes a halo update for a two-dimensional array.
subroutine pass_var_complete_2d(id_update, array, MOM_dom, sideflag, position, halo, &
                                clock)
  integer,                intent(in)    :: id_update !< The integer id of this update which has
                                                    !! been returned from a previous call to
                                                    !! pass_var_start.
  real, dimension(:,:),   intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.

  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_complete_update_domains(id_update, array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position, &
                            whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_complete_update_domains(id_update, array, MOM_dom%mpp_domain, &
                                     flags=dirflag, position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine pass_var_complete_2d

!> pass_var_complete_3d completes a halo update for a three-dimensional array.
subroutine pass_var_complete_3d(id_update, array, MOM_dom, sideflag, position, halo, &
                                clock)
  integer,                intent(in)    :: id_update !< The integer id of this update which has
                                                    !! been returned from a previous call to
                                                    !! pass_var_start.
  real, dimension(:,:,:), intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.

  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_complete_update_domains(id_update, array, MOM_dom%mpp_domain, &
                            flags=dirflag, position=position, &
                            whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_complete_update_domains(id_update, array, MOM_dom%mpp_domain, &
                                     flags=dirflag, position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine pass_var_complete_3d

!> pass_vector_2d does a halo update for a pair of two-dimensional arrays
!! representing the components of a two-dimensional horizontal vector.
subroutine pass_vector_2d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete, halo, &
                          clock)
  real, dimension(:,:),  intent(inout) :: u_cmpt    !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:),  intent(inout) :: v_cmpt    !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type), intent(inout) :: MOM_dom   !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,     optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,     optional, intent(in)    :: stagger   !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  logical,     optional, intent(in)    :: complete  !< An optional argument indicating whether the
                                     !! halo updates should be completed before progress resumes.
                                     !! Omitting complete is the same as setting complete to .true.
  integer,     optional, intent(in)    :: halo      !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,     optional, intent(in)    :: clock     !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.

  ! Local variables
  integer :: stagger_local
  integer :: dirflag
  logical :: block_til_complete

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif
  block_til_complete = .true.
  if (present(complete)) block_til_complete = complete

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_update_domains(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                   gridtype=stagger_local, complete = block_til_complete, &
                   whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_update_domains(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                   gridtype=stagger_local, complete = block_til_complete)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine pass_vector_2d

!> fill_vector_symmetric_edges_2d does an usual set of halo updates that only
!! fill in the values at the edge of a pair of symmetric memory two-dimensional
!! arrays representing the components of a two-dimensional horizontal vector.
!! If symmetric memory is not being used, this subroutine does nothing except to
!! possibly turn optional cpu clocks on or off.
subroutine fill_vector_symmetric_edges_2d(u_cmpt, v_cmpt, MOM_dom, stagger, scalar, &
                                          clock)
  real, dimension(:,:),  intent(inout) :: u_cmpt  !< The nominal zonal (u) component of the vector
                                                  !! pair which is having its halos points
                                                  !! exchanged.
  real, dimension(:,:),  intent(inout) :: v_cmpt  !< The nominal meridional (v) component of the
                                                  !! vector pair which is having its halos points
                                                  !! exchanged.
  type(MOM_domain_type), intent(inout) :: MOM_dom !< The MOM_domain_type containing the mpp_domain
                                                  !! needed to determine where data should be
                                                  !! sent.
  integer,     optional, intent(in)    :: stagger !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  logical,     optional, intent(in)    :: scalar  !< An optional argument indicating whether.
  integer,     optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                   !! started then stopped to time this routine.

  ! Local variables
  integer :: stagger_local
  integer :: dirflag
  integer :: i, j, isc, iec, jsc, jec, isd, ied, jsd, jed, IscB, IecB, JscB, JecB
  real, allocatable, dimension(:) :: sbuff_x, sbuff_y, wbuff_x, wbuff_y
  logical :: block_til_complete

  if (.not. MOM_dom%symmetric) then
      return
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  if (.not.(stagger_local == CGRID_NE .or. stagger_local == BGRID_NE)) return

  call mpp_get_compute_domain(MOM_dom%mpp_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(MOM_dom%mpp_domain, isd, ied, jsd, jed)

  ! Adjust isc, etc., to account for the fact that the input arrays indices all
  ! start at 1 (and are effectively on a SW grid!).
  isc = isc - (isd-1) ; iec = iec - (isd-1)
  jsc = jsc - (jsd-1) ; jec = jec - (jsd-1)
  IscB = isc ; IecB = iec+1 ; JscB = jsc ; JecB = jec+1

  dirflag = To_All ! 60
  if (present(scalar)) then ; if (scalar) dirflag = To_All+SCALAR_PAIR ; endif

  if (stagger_local == CGRID_NE) then
    allocate(wbuff_x(jsc:jec)) ; allocate(sbuff_y(isc:iec))
    wbuff_x(:) = 0.0 ; sbuff_y(:) = 0.0
    call mpp_get_boundary(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                          wbufferx=wbuff_x, sbuffery=sbuff_y, &
                          gridtype=CGRID_NE)
    do i=isc,iec
      v_cmpt(i,JscB) = sbuff_y(i)
    enddo
    do j=jsc,jec
      u_cmpt(IscB,j) = wbuff_x(j)
    enddo
    deallocate(wbuff_x) ; deallocate(sbuff_y)
  elseif  (stagger_local == BGRID_NE) then
    allocate(wbuff_x(JscB:JecB)) ; allocate(sbuff_x(IscB:IecB))
    allocate(wbuff_y(JscB:JecB)) ; allocate(sbuff_y(IscB:IecB))
    wbuff_x(:) = 0.0 ; wbuff_y(:) = 0.0 ; sbuff_x(:) = 0.0 ; sbuff_y(:) = 0.0
    call mpp_get_boundary(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                          wbufferx=wbuff_x, sbufferx=sbuff_x, &
                          wbuffery=wbuff_y, sbuffery=sbuff_y, &
                          gridtype=BGRID_NE)
    do I=IscB,IecB
      u_cmpt(I,JscB) = sbuff_x(I) ; v_cmpt(I,JscB) = sbuff_y(I)
    enddo
    do J=JscB,JecB
      u_cmpt(IscB,J) = wbuff_x(J) ; v_cmpt(IscB,J) = wbuff_y(J)
    enddo
    deallocate(wbuff_x) ; deallocate(sbuff_x)
    deallocate(wbuff_y) ; deallocate(sbuff_y)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine fill_vector_symmetric_edges_2d

!> pass_vector_3d does a halo update for a pair of three-dimensional arrays
!! representing the components of a three-dimensional horizontal vector.
subroutine pass_vector_3d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete, halo, &
                          clock)
  real, dimension(:,:,:), intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:,:), intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                     !! halo updates should be completed before progress resumes.
                                     !! Omitting complete is the same as setting complete to .true.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.

  ! Local variables
  integer :: stagger_local
  integer :: dirflag
  logical :: block_til_complete

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif
  block_til_complete = .true.
  if (present(complete)) block_til_complete = complete

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_update_domains(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                   gridtype=stagger_local, complete = block_til_complete, &
                   whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_update_domains(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                   gridtype=stagger_local, complete = block_til_complete)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine pass_vector_3d

!> pass_vector_start_2d starts a halo update for a pair of two-dimensional arrays
!! representing the components of a two-dimensional horizontal vector.
function pass_vector_start_2d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete, halo, &
                              clock)
  real, dimension(:,:),   intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:),   intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                     !! halo updates should be completed before progress resumes.
                                     !! Omitting complete is the same as setting complete to .true.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  integer                               :: pass_vector_start_2d !< The integer index for this
                                                                !! update.

  ! Local variables
  integer :: stagger_local
  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    pass_vector_start_2d = mpp_start_update_domains(u_cmpt, v_cmpt, &
        MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local, &
        whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    pass_vector_start_2d = mpp_start_update_domains(u_cmpt, v_cmpt, &
        MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end function pass_vector_start_2d

!> pass_vector_start_3d starts a halo update for a pair of three-dimensional arrays
!! representing the components of a three-dimensional horizontal vector.
function pass_vector_start_3d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete, halo, &
                              clock)
  real, dimension(:,:,:), intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:,:), intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  logical,      optional, intent(in)    :: complete !< An optional argument indicating whether the
                                     !! halo updates should be completed before progress resumes.
                                     !! Omitting complete is the same as setting complete to .true.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  integer                               :: pass_vector_start_3d !< The integer index for this
                                                                !! update.
  ! Local variables
  integer :: stagger_local
  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    pass_vector_start_3d = mpp_start_update_domains(u_cmpt, v_cmpt, &
        MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local, &
        whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    pass_vector_start_3d = mpp_start_update_domains(u_cmpt, v_cmpt, &
        MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end function pass_vector_start_3d

!> pass_vector_complete_2d completes a halo update for a pair of two-dimensional arrays
!! representing the components of a two-dimensional horizontal vector.
subroutine pass_vector_complete_2d(id_update, u_cmpt, v_cmpt, MOM_dom, direction, stagger, halo, &
                                   clock)
  integer,                intent(in)    :: id_update !< The integer id of this update which has been
                                                    !! returned from a previous call to
                                                    !! pass_var_start.
  real, dimension(:,:),   intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:),   intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  ! Local variables
  integer :: stagger_local
  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_complete_update_domains(id_update, u_cmpt, v_cmpt, &
             MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local, &
             whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_complete_update_domains(id_update, u_cmpt, v_cmpt, &
             MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine pass_vector_complete_2d

!> pass_vector_complete_3d completes a halo update for a pair of three-dimensional
!! arrays representing the components of a three-dimensional horizontal vector.
subroutine pass_vector_complete_3d(id_update, u_cmpt, v_cmpt, MOM_dom, direction, stagger, halo, &
                                   clock)
  integer,                intent(in)    :: id_update !< The integer id of this update which has been
                                                    !! returned from a previous call to
                                                    !! pass_var_start.
  real, dimension(:,:,:), intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:,:), intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  ! Local variables
  integer :: stagger_local
  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_complete_update_domains(id_update, u_cmpt, v_cmpt, &
             MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local, &
                   whalo=halo, ehalo=halo, shalo=halo, nhalo=halo)
  else
    call mpp_complete_update_domains(id_update, u_cmpt, v_cmpt, &
             MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine pass_vector_complete_3d

!> create_var_group_pass_2d sets up a group of two-dimensional array halo updates.
subroutine create_var_group_pass_2d(group, array, MOM_dom, sideflag, position, &
                                    halo, clock)
  type(group_pass_type),  intent(inout) :: group    !< The data type that store information for
                                                    !! group update. This data will be used in
                                                    !! do_group_pass.
  real, dimension(:,:),   intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  ! Local variables
  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  elseif (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_create_group_update(group, array, MOM_dom%mpp_domain, flags=dirflag, &
                                 position=position, whalo=halo, ehalo=halo, &
                                 shalo=halo, nhalo=halo)
  else
    call mpp_create_group_update(group, array, MOM_dom%mpp_domain, flags=dirflag, &
                                 position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine create_var_group_pass_2d

!> create_var_group_pass_3d sets up a group of three-dimensional array halo updates.
subroutine create_var_group_pass_3d(group, array, MOM_dom, sideflag, position, halo, &
                                    clock)
  type(group_pass_type),  intent(inout) :: group    !< The data type that store information for
                                                    !! group update. This data will be used in
                                                    !! do_group_pass.
  real, dimension(:,:,:), intent(inout) :: array    !< The array which is having its halos points
                                                    !! exchanged.
  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: sideflag !< An optional integer indicating which
      !! directions the data should be sent. It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH.  For example, TO_EAST sends the data to the processor to the east,
      !! so the halos on the western side are filled.  TO_ALL is the default if sideflag is omitted.
  integer,      optional, intent(in)    :: position !< An optional argument indicating the position.
                                                    !! This is CENTER by default and is often CORNER,
                                                    !! but could also be EAST_FACE or NORTH_FACE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  ! Local variables
  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  dirflag = To_All ! 60
  if (present(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  elseif (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_create_group_update(group, array, MOM_dom%mpp_domain, flags=dirflag, &
                                 position=position, whalo=halo, ehalo=halo, &
                                 shalo=halo, nhalo=halo)
  else
    call mpp_create_group_update(group, array, MOM_dom%mpp_domain, flags=dirflag, &
                                 position=position)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine create_var_group_pass_3d

!> create_vector_group_pass_2d sets up a group of two-dimensional vector halo updates.
subroutine create_vector_group_pass_2d(group, u_cmpt, v_cmpt, MOM_dom, direction, stagger, halo, &
                                       clock)
  type(group_pass_type),  intent(inout) :: group    !< The data type that store information for
                                                    !! group update. This data will be used in
                                                    !! do_group_pass.
  real, dimension(:,:),   intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:),   intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.

  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  ! Local variables
  integer :: stagger_local
  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,u_cmpt, v_cmpt)
  elseif (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_create_group_update(group, u_cmpt, v_cmpt, MOM_dom%mpp_domain, &
            flags=dirflag, gridtype=stagger_local, whalo=halo, ehalo=halo, &
            shalo=halo, nhalo=halo)
  else
    call mpp_create_group_update(group, u_cmpt, v_cmpt, MOM_dom%mpp_domain, &
            flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine create_vector_group_pass_2d

!> create_vector_group_pass_3d sets up a group of three-dimensional vector halo updates.
subroutine create_vector_group_pass_3d(group, u_cmpt, v_cmpt, MOM_dom, direction, stagger, halo, &
                                       clock)
  type(group_pass_type),  intent(inout) :: group    !< The data type that store information for
                                                    !! group update. This data will be used in
                                                    !! do_group_pass.
  real, dimension(:,:,:), intent(inout) :: u_cmpt   !< The nominal zonal (u) component of the vector
                                                    !! pair which is having its halos points
                                                    !! exchanged.
  real, dimension(:,:,:), intent(inout) :: v_cmpt   !< The nominal meridional (v) component of the
                                                    !! vector pair which is having its halos points
                                                    !! exchanged.

  type(MOM_domain_type),  intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,      optional, intent(in)    :: direction !< An optional integer indicating which
      !! directions the data should be sent.  It is TO_ALL or the sum of any of TO_EAST, TO_WEST,
      !! TO_NORTH, and TO_SOUTH, possibly plus SCALAR_PAIR if these are paired non-directional
      !! scalars discretized at the typical vector component locations.  For example, TO_EAST sends
      !! the data to the processor to the east, so the halos on the western side are filled. TO_ALL
      !! is the default if omitted.
  integer,      optional, intent(in)    :: stagger  !< An optional flag, which may be one of A_GRID,
                     !! BGRID_NE, or CGRID_NE, indicating where the two components of the vector are
                     !! discretized. Omitting stagger is the same as setting it to CGRID_NE.
  integer,      optional, intent(in)    :: halo     !< The size of the halo to update - the full
                                                    !! halo by default.
  integer,      optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.

  ! Local variables
  integer :: stagger_local
  integer :: dirflag

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (present(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,u_cmpt, v_cmpt)
  elseif (present(halo) .and. MOM_dom%thin_halo_updates) then
    call mpp_create_group_update(group, u_cmpt, v_cmpt, MOM_dom%mpp_domain, &
            flags=dirflag, gridtype=stagger_local, whalo=halo, ehalo=halo, &
            shalo=halo, nhalo=halo)
  else
    call mpp_create_group_update(group, u_cmpt, v_cmpt, MOM_dom%mpp_domain, &
            flags=dirflag, gridtype=stagger_local)
  endif

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine create_vector_group_pass_3d

!> do_group_pass carries out a group halo update.
subroutine do_group_pass(group, MOM_dom, clock)
  type(group_pass_type), intent(inout) :: group     !< The data type that store information for
                                                    !! group update. This data will be used in
                                                    !! do_group_pass.
  type(MOM_domain_type), intent(inout) :: MOM_dom   !< The MOM_domain_type containing the mpp_domain
                                                    !! needed to determine where data should be
                                                    !! sent.
  integer,     optional, intent(in)    :: clock     !< The handle for a cpu time clock that should be
                                                    !! started then stopped to time this routine.
  real :: d_type

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  call mpp_do_group_update(group, MOM_dom%mpp_domain, d_type)

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine do_group_pass

!> start_group_pass starts out a group halo update.
subroutine start_group_pass(group, MOM_dom, clock)
  type(group_pass_type), intent(inout) :: group    !< The data type that store information for
                                                   !! group update. This data will be used in
                                                   !! do_group_pass.
  type(MOM_domain_type), intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                   !! needed to determine where data should be
                                                   !! sent.
  integer,     optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                   !! started then stopped to time this routine.

  real                                 :: d_type

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  call mpp_start_group_update(group, MOM_dom%mpp_domain, d_type)

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine start_group_pass

!> complete_group_pass completes a group halo update.
subroutine complete_group_pass(group, MOM_dom, clock)
  type(group_pass_type), intent(inout) :: group    !< The data type that store information for
                                                   !! group update. This data will be used in
                                                   !! do_group_pass.
  type(MOM_domain_type), intent(inout) :: MOM_dom  !< The MOM_domain_type containing the mpp_domain
                                                   !! needed to determine where data should be
                                                   !! sent.
  integer,     optional, intent(in)    :: clock    !< The handle for a cpu time clock that should be
                                                   !! started then stopped to time this routine.
  real                                 :: d_type

  if (present(clock)) then ; if (clock>0) call cpu_clock_begin(clock) ; endif

  call mpp_complete_group_update(group, MOM_dom%mpp_domain, d_type)

  if (present(clock)) then ; if (clock>0) call cpu_clock_end(clock) ; endif

end subroutine complete_group_pass


!> Pass a 2-D array from one MOM domain to another
subroutine redistribute_array_2d(Domain1, array1, Domain2, array2, complete)
  type(domain2d), &
           intent(in)  :: Domain1 !< The MOM domain from which to extract information.
  real, dimension(:,:), intent(in) :: array1 !< The array from which to extract information.
  type(domain2d), &
           intent(in)  :: Domain2 !< The MOM domain receiving information.
  real, dimension(:,:), intent(out) :: array2 !< The array receiving information.
  logical, optional, intent(in) :: complete  !< If true, finish communication before proceeding.

  ! Local variables
  logical :: do_complete

  do_complete=.true.;if (PRESENT(complete)) do_complete = complete

  call mpp_redistribute(Domain1, array1, Domain2, array2, do_complete)

end subroutine redistribute_array_2d

!> Pass a 3-D array from one MOM domain to another
subroutine redistribute_array_3d(Domain1, array1, Domain2, array2, complete)
  type(domain2d), &
           intent(in)  :: Domain1 !< The MOM domain from which to extract information.
  real, dimension(:,:,:), intent(in) :: array1 !< The array from which to extract information.
  type(domain2d), &
           intent(in)  :: Domain2 !< The MOM domain receiving information.
  real, dimension(:,:,:), intent(out) :: array2 !< The array receiving information.
  logical, optional, intent(in) :: complete  !< If true, finish communication before proceeding.

  ! Local variables
  logical :: do_complete

  do_complete=.true.;if (PRESENT(complete)) do_complete = complete

  call mpp_redistribute(Domain1, array1, Domain2, array2, do_complete)

end subroutine redistribute_array_3d


!> Rescale the values of a 4-D array in its computational domain by a constant factor
subroutine rescale_comp_data_4d(domain, array, scale)
  type(MOM_domain_type),    intent(in)    :: domain !< MOM domain from which to extract information
  real, dimension(:,:,:,:), intent(inout) :: array  !< The array which is having the data in its
                                                    !! computational domain rescaled
  real,                     intent(in)    :: scale  !< A scaling factor by which to multiply the
                                                    !! values in the computational domain of array
  integer :: is, ie, js, je

  if (scale == 1.0) return

  call get_simple_array_i_ind(domain, size(array,1), is, ie)
  call get_simple_array_j_ind(domain, size(array,2), js, je)
  array(is:ie,js:je,:,:) = scale*array(is:ie,js:je,:,:)

end subroutine rescale_comp_data_4d

!> Rescale the values of a 3-D array in its computational domain by a constant factor
subroutine rescale_comp_data_3d(domain, array, scale)
  type(MOM_domain_type),  intent(in)    :: domain !< MOM domain from which to extract information
  real, dimension(:,:,:), intent(inout) :: array  !< The array which is having the data in its
                                                  !! computational domain rescaled
  real,                   intent(in)    :: scale  !< A scaling factor by which to multiply the
                                                  !! values in the computational domain of array
  integer :: is, ie, js, je

  if (scale == 1.0) return

  call get_simple_array_i_ind(domain, size(array,1), is, ie)
  call get_simple_array_j_ind(domain, size(array,2), js, je)
  array(is:ie,js:je,:) = scale*array(is:ie,js:je,:)

end subroutine rescale_comp_data_3d

!> Rescale the values of a 2-D array in its computational domain by a constant factor
subroutine rescale_comp_data_2d(domain, array, scale)
  type(MOM_domain_type), intent(in)    :: domain !< MOM domain from which to extract information
  real, dimension(:,:),  intent(inout) :: array  !< The array which is having the data in its
                                                 !! computational domain rescaled
  real,                  intent(in)    :: scale  !< A scaling factor by which to multiply the
                                                 !! values in the computational domain of array
  integer :: is, ie, js, je

  if (scale == 1.0) return

  call get_simple_array_i_ind(domain, size(array,1), is, ie)
  call get_simple_array_j_ind(domain, size(array,2), js, je)
  array(is:ie,js:je) = scale*array(is:ie,js:je)

end subroutine rescale_comp_data_2d

!> create_MOM_domain creates and initializes a MOM_domain_type variables, based on the information
!! provided in arguments.
subroutine create_MOM_domain(MOM_dom, n_global, n_halo, reentrant, tripolar_N, layout, io_layout, &
                             domain_name, mask_table, symmetric, thin_halos, nonblocking)
  type(MOM_domain_type),      pointer    :: MOM_dom   !< A pointer to the MOM_domain_type being defined here.
  integer, dimension(2),      intent(in) :: n_global  !< The number of points on the global grid in
                                                      !! the i- and j-directions
  integer, dimension(2),      intent(in) :: n_halo    !< The number of halo points on each processor
  logical, dimension(2),      intent(in) :: reentrant !< If true the grid is periodic in the i- and j- directions
  logical,                    intent(in) :: tripolar_N !< If true the grid uses northern tripolar connectivity
  integer, dimension(2),      intent(in) :: layout    !< The layout of logical PEs in the i- and j-directions.
  integer, dimension(2), optional, intent(in) :: io_layout !< The layout for parallel input and output.
  character(len=*), optional, intent(in) :: domain_name !< A name for this domain, "MOM" if missing.
  character(len=*), optional, intent(in) :: mask_table !< The full relative or absolute path to the mask table.
  logical,          optional, intent(in) :: symmetric !< If present, this specifies whether this domain
                                                      !! uses symmetric memory, or true if missing.
  logical,          optional, intent(in) :: thin_halos !< If present, this specifies whether to permit the use of
                                                      !! thin halo updates, or true if missing.
  logical,          optional, intent(in) :: nonblocking !< If present, this specifies whether to permit the use of
                                                      !! nonblocking halo updates, or false if missing.

  ! local variables
  integer, dimension(4) :: global_indices ! The lower and upper global i- and j-index bounds
  integer :: X_FLAGS  ! A combination of integers encoding the x-direction grid connectivity.
  integer :: Y_FLAGS  ! A combination of integers encoding the y-direction grid connectivity.
  integer :: xhalo_d2, yhalo_d2
  character(len=200) :: mesg    ! A string for use in error messages
  logical :: mask_table_exists  ! Mask_table is present and the file it points to exists

  if (.not.associated(MOM_dom)) then
    allocate(MOM_dom)
    allocate(MOM_dom%mpp_domain)
    allocate(MOM_dom%mpp_domain_d2)
  endif

  MOM_dom%name = "MOM" ; if (present(domain_name)) MOM_dom%name = trim(domain_name)

  X_FLAGS = 0 ; Y_FLAGS = 0
  if (reentrant(1)) X_FLAGS = CYCLIC_GLOBAL_DOMAIN
  if (reentrant(2)) Y_FLAGS = CYCLIC_GLOBAL_DOMAIN
  if (tripolar_N) then
    Y_FLAGS = FOLD_NORTH_EDGE
    if (reentrant(2)) call MOM_error(FATAL,"MOM_domains: "// &
      "TRIPOLAR_N and REENTRANT_Y may not be used together.")
  endif

  MOM_dom%nonblocking_updates = nonblocking
  MOM_dom%thin_halo_updates = thin_halos
  MOM_dom%symmetric = .true. ; if (present(symmetric)) MOM_dom%symmetric = symmetric
  MOM_dom%niglobal = n_global(1) ; MOM_dom%njglobal = n_global(2)
  MOM_dom%nihalo = n_halo(1) ; MOM_dom%njhalo = n_halo(2)

  ! Save the extra data for creating other domains of different resolution that overlay this domain.
  MOM_dom%X_FLAGS = X_FLAGS
  MOM_dom%Y_FLAGS = Y_FLAGS
  MOM_dom%layout(:) = layout(:)

  ! Set up the io_layout, with error handling.
  MOM_dom%io_layout(:) = (/ 1, 1 /)
  if (present(io_layout)) then
    if (io_layout(1) == 0) then
      MOM_dom%io_layout(1) = layout(1)
    elseif (io_layout(1) > 1) then
      MOM_dom%io_layout(1) = io_layout(1)
      if (modulo(layout(1), io_layout(1)) /= 0) then
        write(mesg,'("MOM_domains_init: The i-direction I/O-layout, IO_LAYOUT(1)=",i4, &
              &", does not evenly divide the i-direction layout, NIPROC=,",i4,".")') io_layout(1), layout(1)
        call MOM_error(FATAL, mesg)
      endif
    endif

    if (io_layout(2) == 0) then
      MOM_dom%io_layout(2) = layout(2)
    elseif (io_layout(2) > 1) then
      MOM_dom%io_layout(2) = io_layout(2)
      if (modulo(layout(2), io_layout(2)) /= 0) then
        write(mesg,'("MOM_domains_init: The j-direction I/O-layout, IO_LAYOUT(2)=",i4, &
              &", does not evenly divide the j-direction layout, NJPROC=,",i4,".")') io_layout(2), layout(2)
        call MOM_error(FATAL, mesg)
      endif
    endif
  endif

  if (present(mask_table)) then
    mask_table_exists = file_exist(mask_table)
    if (mask_table_exists) then
      allocate(MOM_dom%maskmap(layout(1), layout(2)))
      call parse_mask_table(mask_table, MOM_dom%maskmap, MOM_dom%name)
    endif
  else
    mask_table_exists = .false.
  endif

  call clone_MD_to_d2D(MOM_dom, MOM_dom%mpp_domain)

  !For downsampled domain, recommend a halo of 1 (or 0?) since we're not doing wide-stencil computations.
  !But that does not work because the downsampled field would not have the correct size to pass the checks, e.g., we get
  !error: downsample_diag_indices_get: peculiar size 28 in i-direction\ndoes not match one of 24 25 26 27
  ! call clone_MD_to_d2D(MOM_dom, MOM_dom%mpp_domain_d2, halo_size=(MOM_dom%nihalo/2), coarsen=2)
  call clone_MD_to_d2D(MOM_dom, MOM_dom%mpp_domain_d2, coarsen=2)

end subroutine create_MOM_domain

!> dealloc_MOM_domain deallocates memory associated with a pointer to a MOM_domain_type
!! and potentially all of its contents
subroutine deallocate_MOM_domain(MOM_domain, cursory)
  type(MOM_domain_type), pointer :: MOM_domain !< A pointer to the MOM_domain_type being deallocated
  logical,  optional, intent(in) :: cursory    !< If true do not deallocate fields associated
                                               !! with the underlying infrastructure
  logical :: invasive  ! If true, deallocate fields associated with the underlying infrastructure

  invasive = .true. ; if (present(cursory)) invasive = .not.cursory

  if (associated(MOM_domain)) then
    if (associated(MOM_domain%mpp_domain)) then
      if (invasive) call mpp_deallocate_domain(MOM_domain%mpp_domain)
      deallocate(MOM_domain%mpp_domain)
    endif
    if (associated(MOM_domain%mpp_domain_d2)) then
      if (invasive) call mpp_deallocate_domain(MOM_domain%mpp_domain_d2)
      deallocate(MOM_domain%mpp_domain_d2)
    endif
    if (associated(MOM_domain%maskmap)) deallocate(MOM_domain%maskmap)
    deallocate(MOM_domain)
  endif

end subroutine deallocate_MOM_domain

!> MOM_thread_affinity_set returns true if the number of openMP threads have been set to a value greater than 1.
function MOM_thread_affinity_set()
  ! Local variables
  !$ integer :: ocean_nthreads       ! Number of openMP threads
  !$ integer :: omp_get_num_threads  ! An openMP function that returns the number of threads
  logical :: MOM_thread_affinity_set

  MOM_thread_affinity_set = .false.
  !$ call fms_affinity_init()
  !$OMP PARALLEL
  !$OMP   MASTER
  !$        ocean_nthreads = omp_get_num_threads()
  !$OMP   END MASTER
  !$OMP END PARALLEL
  !$ MOM_thread_affinity_set = (ocean_nthreads > 1 )
end function MOM_thread_affinity_set

!> set_MOM_thread_affinity sets the number of openMP threads to use with the ocean.
subroutine set_MOM_thread_affinity(ocean_nthreads, ocean_hyper_thread)
  integer, intent(in) :: ocean_nthreads     !< Number of openMP threads to use for the ocean model
  logical, intent(in) :: ocean_hyper_thread !< If true, use hyper threading

  ! Local variables
  !$ integer :: omp_get_thread_num, omp_get_num_threads !< These are the results of openMP functions

  !$ call fms_affinity_init()  ! fms_affinity_init can be safely called more than once.
  !$ call fms_affinity_set('OCEAN', ocean_hyper_thread, ocean_nthreads)
  !$ call omp_set_num_threads(ocean_nthreads)
  !$OMP PARALLEL
  !$ write(6,*) "MOM_domains_mod OMPthreading ", fms_affinity_get(), omp_get_thread_num(), omp_get_num_threads()
  !$ flush(6)
  !$OMP END PARALLEL
end subroutine set_MOM_thread_affinity

!> This subroutine retrieves the 1-d domains that make up the 2d-domain in a MOM_domain
subroutine get_domain_components_MD(MOM_dom, x_domain, y_domain)
  type(MOM_domain_type),    intent(in)    :: MOM_dom  !< The MOM_domain whose contents are being extracted
  type(domain1D), optional, intent(inout) :: x_domain !< The 1-d logical x-domain
  type(domain1D), optional, intent(inout) :: y_domain !< The 1-d logical y-domain

  call mpp_get_domain_components(MOM_dom%mpp_domain, x_domain, y_domain)
end subroutine get_domain_components_MD

!> This subroutine retrieves the 1-d domains that make up a 2d-domain
subroutine get_domain_components_d2D(domain, x_domain, y_domain)
  type(domain2D),           intent(in)    :: domain  !< The 2D domain whose contents are being extracted
  type(domain1D), optional, intent(inout) :: x_domain !< The 1-d logical x-domain
  type(domain1D), optional, intent(inout) :: y_domain !< The 1-d logical y-domain

  call mpp_get_domain_components(domain, x_domain, y_domain)
end subroutine get_domain_components_d2D

!> clone_MD_to_MD copies one MOM_domain_type into another, while allowing
!! some properties of the new type to differ from the original one.
subroutine clone_MD_to_MD(MD_in, MOM_dom, min_halo, halo_size, symmetric, domain_name, &
                          turns, refine, extra_halo)
  type(MOM_domain_type), intent(in)    :: MD_in  !< An existing MOM_domain
  type(MOM_domain_type), pointer       :: MOM_dom !< A pointer to a MOM_domain that will be
                                  !! allocated if it is unassociated, and will have data
                                  !! copied from MD_in
  integer, dimension(2), &
               optional, intent(inout) :: min_halo !< If present, this sets the
                                  !! minimum halo size for this domain in the i- and j-
                                  !! directions, and returns the actual halo size used.
  integer,     optional, intent(in)    :: halo_size !< If present, this sets the halo
                                  !! size for the domain in the i- and j-directions.
                                  !! min_halo and halo_size can not both be present.
  logical,     optional, intent(in)    :: symmetric !< If present, this specifies
                                  !! whether the new domain is symmetric, regardless of
                                  !! whether the macro SYMMETRIC_MEMORY_ is defined.
  character(len=*), &
               optional, intent(in)    :: domain_name !< A name for the new domain, copied
                                  !! from MD_in if missing.
  integer, optional, intent(in) :: turns   !< Number of quarter turns
  integer, optional, intent(in) :: refine  !< A factor by which to enhance the grid resolution.
  integer, optional, intent(in) :: extra_halo !< An extra number of points in the halos
                                  !! compared with MD_in

  integer :: global_indices(4)
  logical :: mask_table_exists
  integer, dimension(:), allocatable :: exni ! The extents of the grid for each i-row of the layout.
                                             ! The sum of exni must equal MOM_dom%niglobal.
  integer, dimension(:), allocatable :: exnj ! The extents of the grid for each j-row of the layout.
                                             ! The sum of exni must equal MOM_dom%niglobal.
  integer :: qturns ! The number of quarter turns, restricted to the range of 0 to 3.
  integer :: i, j, nl1, nl2

  qturns = 0
  if (present(turns)) qturns = modulo(turns, 4)

  if (.not.associated(MOM_dom)) then
    allocate(MOM_dom)
    allocate(MOM_dom%mpp_domain)
    allocate(MOM_dom%mpp_domain_d2)
  endif

! Save the extra data for creating other domains of different resolution that overlay this domain
  MOM_dom%symmetric = MD_in%symmetric
  MOM_dom%nonblocking_updates = MD_in%nonblocking_updates
  MOM_dom%thin_halo_updates = MD_in%thin_halo_updates

  if (modulo(qturns, 2) /= 0) then
    MOM_dom%niglobal = MD_in%njglobal ; MOM_dom%njglobal = MD_in%niglobal
    MOM_dom%nihalo = MD_in%njhalo ; MOM_dom%njhalo = MD_in%nihalo
    call get_layout_extents(MD_in, exnj, exni)

    MOM_dom%X_FLAGS = MD_in%Y_FLAGS ; MOM_dom%Y_FLAGS = MD_in%X_FLAGS
    MOM_dom%layout(:) = MD_in%layout(2:1:-1)
    MOM_dom%io_layout(:) = MD_in%io_layout(2:1:-1)
  else
    MOM_dom%niglobal = MD_in%niglobal ; MOM_dom%njglobal = MD_in%njglobal
    MOM_dom%nihalo = MD_in%nihalo ; MOM_dom%njhalo = MD_in%njhalo
    call get_layout_extents(MD_in, exni, exnj)

    MOM_dom%X_FLAGS = MD_in%X_FLAGS ; MOM_dom%Y_FLAGS = MD_in%Y_FLAGS
    MOM_dom%layout(:) = MD_in%layout(:)
    MOM_dom%io_layout(:) = MD_in%io_layout(:)
  endif

  ! Ensure that the points per processor are the same on the source and densitation grids.
  select case (qturns)
    case (1) ; call invert(exni)
    case (2) ; call invert(exni) ; call invert(exnj)
    case (3) ; call invert(exnj)
  end select

  if (associated(MD_in%maskmap)) then
    mask_table_exists = .true.
    allocate(MOM_dom%maskmap(MOM_dom%layout(1), MOM_dom%layout(2)))

    nl1 = MOM_dom%layout(1) ; nl2 = MOM_dom%layout(2)
    select case (qturns)
      case (0)
        do j=1,nl2 ; do i=1,nl1
          MOM_dom%maskmap(i,j) = MD_in%maskmap(i, j)
        enddo ; enddo
      case (1)
        do j=1,nl2 ; do i=1,nl1
          MOM_dom%maskmap(i,j) = MD_in%maskmap(j, nl1+1-i)
        enddo ; enddo
      case (2)
        do j=1,nl2 ; do i=1,nl1
          MOM_dom%maskmap(i,j) = MD_in%maskmap(nl1+1-i, nl2+1-j)
        enddo ; enddo
      case (3)
        do j=1,nl2 ; do i=1,nl1
          MOM_dom%maskmap(i,j) = MD_in%maskmap(nl2+1-j, i)
        enddo ; enddo
    end select
  else
    mask_table_exists = .false.
  endif

  ! Optionally enhance the grid resolution.
  if (present(refine)) then ; if (refine > 1) then
    MOM_dom%niglobal = refine*MOM_dom%niglobal ; MOM_dom%njglobal = refine*MOM_dom%njglobal
    MOM_dom%nihalo = refine*MOM_dom%nihalo ; MOM_dom%njhalo = refine*MOM_dom%njhalo
    do i=1,MOM_dom%layout(1) ; exni(i) = refine*exni(i) ; enddo
    do j=1,MOM_dom%layout(2) ; exnj(j) = refine*exnj(j) ; enddo
  endif ; endif

  ! Optionally enhance the grid resolution.
  if (present(extra_halo)) then ; if (extra_halo > 0) then
    MOM_dom%nihalo = MOM_dom%nihalo + extra_halo ; MOM_dom%njhalo = MOM_dom%njhalo + extra_halo
  endif ; endif

  if (present(halo_size) .and. present(min_halo)) call MOM_error(FATAL, &
      "clone_MOM_domain can not have both halo_size and min_halo present.")

  if (present(min_halo)) then
    MOM_dom%nihalo = max(MOM_dom%nihalo, min_halo(1))
    min_halo(1) = MOM_dom%nihalo
    MOM_dom%njhalo = max(MOM_dom%njhalo, min_halo(2))
    min_halo(2) = MOM_dom%njhalo
  endif

  if (present(halo_size)) then
    MOM_dom%nihalo = halo_size ; MOM_dom%njhalo = halo_size
  endif

  if (present(symmetric)) then ; MOM_dom%symmetric = symmetric ; endif

  if (present(domain_name)) then
    MOM_dom%name = trim(domain_name)
  else
    MOM_dom%name = MD_in%name
  endif

  call clone_MD_to_d2D(MOM_dom, MOM_dom%mpp_domain, xextent=exni, yextent=exnj)

  call clone_MD_to_d2D(MOM_dom, MOM_dom%mpp_domain_d2, domain_name=MOM_dom%name, coarsen=2)

end subroutine clone_MD_to_MD


!> clone_MD_to_d2D uses information from a MOM_domain_type to create a new
!! domain2d type, while allowing some properties of the new type to differ from
!! the original one.
subroutine clone_MD_to_d2D(MD_in, mpp_domain, min_halo, halo_size, symmetric, &
                           domain_name, turns, xextent, yextent, coarsen)
  type(MOM_domain_type), intent(in)    :: MD_in !< An existing MOM_domain to be cloned
  type(domain2d),        intent(inout) :: mpp_domain !< The new mpp_domain to be set up
  integer, dimension(2), &
               optional, intent(inout) :: min_halo !< If present, this sets the
                                  !! minimum halo size for this domain in the i- and j-
                                  !! directions, and returns the actual halo size used.
  integer,     optional, intent(in)    :: halo_size !< If present, this sets the halo
                                  !! size for the domain in the i- and j-directions.
                                  !! min_halo and halo_size can not both be present.
  logical,     optional, intent(in)    :: symmetric !< If present, this specifies
                                  !! whether the new domain is symmetric, regardless of
                                  !! whether the macro SYMMETRIC_MEMORY_ is defined or
                                  !! whether MD_in is symmetric.
  character(len=*), &
               optional, intent(in)    :: domain_name !< A name for the new domain, "MOM"
                                  !! if missing.
  integer, optional, intent(in) :: turns  !< Number of quarter turns - not implemented here.
  integer, optional, intent(in) :: coarsen !< A factor by which to coarsen this grid.
                                  !! The default of 1 is for no coarsening.
  integer, dimension(:), optional, intent(in) :: xextent !< The number of grid points in the
                                  !! tracer computational domain for division of the x-layout.
  integer, dimension(:), optional, intent(in) :: yextent !< The number of grid points in the
                                  !! tracer computational domain for division of the y-layout.

  integer :: global_indices(4)
  integer :: nihalo, njhalo
  logical :: symmetric_dom, do_coarsen
  character(len=64) :: dom_name

  if (present(turns)) &
    call MOM_error(FATAL, "Rotation not supported for MOM_domain to domain2d")

  if (present(halo_size) .and. present(min_halo)) call MOM_error(FATAL, &
      "clone_MOM_domain can not have both halo_size and min_halo present.")

  do_coarsen = .false. ; if (present(coarsen)) then ; do_coarsen = (coarsen > 1) ; endif

  nihalo = MD_in%nihalo ; njhalo = MD_in%njhalo
  if (do_coarsen) then
    nihalo = int(MD_in%nihalo / coarsen) ; njhalo = int(MD_in%njhalo / coarsen)
  endif

  if (present(min_halo)) then
    nihalo = max(nihalo, min_halo(1))
    njhalo = max(njhalo, min_halo(2))
    min_halo(1) = nihalo ; min_halo(2) = njhalo
  endif
  if (present(halo_size)) then
    nihalo = halo_size ; njhalo = halo_size
  endif

  symmetric_dom = MD_in%symmetric
  if (present(symmetric)) then ; symmetric_dom = symmetric ; endif

  dom_name = MD_in%name
  if (do_coarsen) dom_name = trim(MD_in%name)//"c"
  if (present(domain_name)) dom_name = trim(domain_name)

  global_indices(1:4) = (/ 1, MD_in%niglobal, 1, MD_in%njglobal /)
  if (do_coarsen) then
    global_indices(1:4) = (/ 1, (MD_in%niglobal/coarsen), 1, (MD_in%njglobal/coarsen) /)
  endif

  if (associated(MD_in%maskmap)) then
    call mpp_define_domains( global_indices, MD_in%layout, mpp_domain, &
                xflags=MD_in%X_FLAGS, yflags=MD_in%Y_FLAGS, xhalo=nihalo, yhalo=njhalo, &
                xextent=xextent, yextent=yextent, symmetry=symmetric_dom, name=dom_name, &
                maskmap=MD_in%maskmap )
  else
    call mpp_define_domains( global_indices, MD_in%layout, mpp_domain, &
                xflags=MD_in%X_FLAGS, yflags=MD_in%Y_FLAGS, xhalo=nihalo, yhalo=njhalo, &
                symmetry=symmetric_dom, xextent=xextent, yextent=yextent, name=dom_name)
  endif

  if ((MD_in%io_layout(1) + MD_in%io_layout(2) > 0) .and. &
      (MD_in%layout(1)*MD_in%layout(2) > 1)) then
    call mpp_define_io_domain(mpp_domain, MD_in%io_layout)
  else
    call mpp_define_io_domain(mpp_domain, (/ 1, 1 /) )
  endif

end subroutine clone_MD_to_d2D

!> Returns the index ranges that have been stored in a MOM_domain_type
subroutine get_domain_extent_MD(Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, &
                                isg, ieg, jsg, jeg, idg_offset, jdg_offset, &
                                symmetric, local_indexing, index_offset, coarsen)
  type(MOM_domain_type), &
                     intent(in)  :: Domain !< The MOM domain from which to extract information
  integer,           intent(out) :: isc    !< The start i-index of the computational domain
  integer,           intent(out) :: iec    !< The end i-index of the computational domain
  integer,           intent(out) :: jsc    !< The start j-index of the computational domain
  integer,           intent(out) :: jec    !< The end j-index of the computational domain
  integer,           intent(out) :: isd    !< The start i-index of the data domain
  integer,           intent(out) :: ied    !< The end i-index of the data domain
  integer,           intent(out) :: jsd    !< The start j-index of the data domain
  integer,           intent(out) :: jed    !< The end j-index of the data domain
  integer, optional, intent(out) :: isg    !< The start i-index of the global domain
  integer, optional, intent(out) :: ieg    !< The end i-index of the global domain
  integer, optional, intent(out) :: jsg    !< The start j-index of the global domain
  integer, optional, intent(out) :: jeg    !< The end j-index of the global domain
  integer, optional, intent(out) :: idg_offset !< The offset between the corresponding global and
                                           !! data i-index spaces.
  integer, optional, intent(out) :: jdg_offset !< The offset between the corresponding global and
                                           !! data j-index spaces.
  logical, optional, intent(out) :: symmetric  !< True if symmetric memory is used.
  logical, optional, intent(in)  :: local_indexing !< If true, local tracer array indices start at 1,
                                           !! as in most MOM6 code.  The default is true.
  integer, optional, intent(in)  :: index_offset   !< A fixed additional offset to all indices. This
                                           !! can be useful for some types of debugging with
                                           !! dynamic memory allocation.  The default is 0.
  integer, optional, intent(in)  :: coarsen !< A factor by which the grid is coarsened.
                                           !!  The default is 1, for no coarsening.

  ! Local variables
  integer :: isg_, ieg_, jsg_, jeg_
  integer :: ind_off, idg_off, jdg_off, coarsen_lev
  logical :: local

  local = .true. ; if (present(local_indexing)) local = local_indexing
  ind_off = 0 ; if (present(index_offset)) ind_off = index_offset

  coarsen_lev = 1 ; if (present(coarsen)) coarsen_lev = coarsen

  if (coarsen_lev == 1) then
    call mpp_get_compute_domain(Domain%mpp_domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain%mpp_domain, isd, ied, jsd, jed)
    call mpp_get_global_domain(Domain%mpp_domain, isg_, ieg_, jsg_, jeg_)
  elseif (coarsen_lev == 2) then
    if (.not.associated(Domain%mpp_domain_d2)) call MOM_error(FATAL, &
            "get_domain_extent called with coarsen=2, but Domain%mpp_domain_d2 is not associated.")
    call mpp_get_compute_domain(Domain%mpp_domain_d2, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain%mpp_domain_d2, isd, ied, jsd, jed)
    call mpp_get_global_domain(Domain%mpp_domain_d2, isg_, ieg_, jsg_, jeg_)
  else
    call MOM_error(FATAL, "get_domain_extent called with an unsupported level of coarsening.")
  endif

  if (local) then
    ! This code institutes the MOM convention that local array indices start at 1.
    idg_off = isd - 1 ; jdg_off = jsd - 1
    isc = isc - isd + 1 ; iec = iec - isd + 1 ; jsc = jsc - jsd + 1 ; jec = jec - jsd + 1
    ied = ied - isd + 1 ; jed = jed - jsd + 1
    isd = 1 ; jsd = 1
  else
    idg_off = 0 ; jdg_off = 0
  endif
  if (ind_off /= 0) then
    idg_off = idg_off + ind_off ; jdg_off = jdg_off + ind_off
    isc = isc + ind_off ; iec = iec + ind_off
    jsc = jsc + ind_off ; jec = jec + ind_off
    isd = isd + ind_off ; ied = ied + ind_off
    jsd = jsd + ind_off ; jed = jed + ind_off
  endif
  if (present(isg)) isg = isg_
  if (present(ieg)) ieg = ieg_
  if (present(jsg)) jsg = jsg_
  if (present(jeg)) jeg = jeg_
  if (present(idg_offset)) idg_offset = idg_off
  if (present(jdg_offset)) jdg_offset = jdg_off
  if (present(symmetric)) symmetric = Domain%symmetric

end subroutine get_domain_extent_MD

!> Returns the index ranges that have been stored in a domain2D type
subroutine get_domain_extent_d2D(Domain, isc, iec, jsc, jec, isd, ied, jsd, jed)
  type(domain2d),    intent(in)  :: Domain !< The MOM domain from which to extract information
  integer,           intent(out) :: isc    !< The start i-index of the computational domain
  integer,           intent(out) :: iec    !< The end i-index of the computational domain
  integer,           intent(out) :: jsc    !< The start j-index of the computational domain
  integer,           intent(out) :: jec    !< The end j-index of the computational domain
  integer, optional, intent(out) :: isd    !< The start i-index of the data domain
  integer, optional, intent(out) :: ied    !< The end i-index of the data domain
  integer, optional, intent(out) :: jsd    !< The start j-index of the data domain
  integer, optional, intent(out) :: jed    !< The end j-index of the data domain

  ! Local variables
  integer :: isd_, ied_, jsd_, jed_, jsg_, jeg_, isg_, ieg_

  call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(Domain, isd_, ied_, jsd_, jed_)

  if (present(isd)) isd = isd_
  if (present(ied)) ied = ied_
  if (present(jsd)) jsd = jsd_
  if (present(jed)) jed = jed_

end subroutine get_domain_extent_d2D

!> Return the (potentially symmetric) computational domain i-bounds for an array
!! passed without index specifications (i.e. indices start at 1) based on an array size.
subroutine get_simple_array_i_ind(domain, size, is, ie, symmetric)
  type(MOM_domain_type), intent(in)  :: domain !< MOM domain from which to extract information
  integer,               intent(in)  :: size   !< The i-array size
  integer,               intent(out) :: is     !< The computational domain starting i-index.
  integer,               intent(out) :: ie     !< The computational domain ending i-index.
  logical,     optional, intent(in)  :: symmetric !< If present, indicates whether symmetric sizes
                                               !! can be considered.
  ! Local variables
  logical :: sym
  character(len=120) :: mesg, mesg2
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed

  call mpp_get_compute_domain(Domain%mpp_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(Domain%mpp_domain, isd, ied, jsd, jed)

  isc = isc-isd+1 ; iec = iec-isd+1 ; ied = ied-isd+1 ; isd = 1
  sym = Domain%symmetric ; if (present(symmetric)) sym = symmetric

  if (size == ied) then ; is = isc ; ie = iec
  elseif (size == 1+iec-isc) then ; is = 1 ; ie = size
  elseif (sym .and. (size == 1+ied)) then ; is = isc ; ie = iec+1
  elseif (sym .and. (size == 2+iec-isc)) then ; is = 1 ; ie = size+1
  else
    write(mesg,'("Unrecognized size ", i6, "in call to get_simple_array_i_ind.  \")') size
    if (sym) then
      write(mesg2,'("Valid sizes are : ", 2i7)') ied, 1+iec-isc
    else
      write(mesg2,'("Valid sizes are : ", 4i7)') ied, 1+iec-isc, 1+ied, 2+iec-isc
    endif
    call MOM_error(FATAL, trim(mesg)//trim(mesg2))
  endif

end subroutine get_simple_array_i_ind


!> Return the (potentially symmetric) computational domain j-bounds for an array
!! passed without index specifications (i.e. indices start at 1) based on an array size.
subroutine get_simple_array_j_ind(domain, size, js, je, symmetric)
  type(MOM_domain_type), intent(in)  :: domain !< MOM domain from which to extract information
  integer,               intent(in)  :: size   !< The j-array size
  integer,               intent(out) :: js     !< The computational domain starting j-index.
  integer,               intent(out) :: je     !< The computational domain ending j-index.
  logical,     optional, intent(in)  :: symmetric !< If present, indicates whether symmetric sizes
                                               !! can be considered.
  ! Local variables
  logical :: sym
  character(len=120) :: mesg, mesg2
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed

  call mpp_get_compute_domain(Domain%mpp_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(Domain%mpp_domain, isd, ied, jsd, jed)

  jsc = jsc-jsd+1 ; jec = jec-jsd+1 ; jed = jed-jsd+1 ; jsd = 1
  sym = Domain%symmetric ; if (present(symmetric)) sym = symmetric

  if (size == jed) then ; js = jsc ; je = jec
  elseif (size == 1+jec-jsc) then ; js = 1 ; je = size
  elseif (sym .and. (size == 1+jed)) then ; js = jsc ; je = jec+1
  elseif (sym .and. (size == 2+jec-jsc)) then ; js = 1 ; je = size+1
  else
    write(mesg,'("Unrecognized size ", i6, "in call to get_simple_array_j_ind.  \")') size
    if (sym) then
      write(mesg2,'("Valid sizes are : ", 2i7)') jed, 1+jec-jsc
    else
      write(mesg2,'("Valid sizes are : ", 4i7)') jed, 1+jec-jsc, 1+jed, 2+jec-jsc
    endif
    call MOM_error(FATAL, trim(mesg)//trim(mesg2))
  endif

end subroutine get_simple_array_j_ind

!> Invert the contents of a 1-d array
subroutine invert(array)
 integer, dimension(:), intent(inout) :: array !< The 1-d array to invert
 integer :: i, ni, swap
 ni = size(array)
 do i=1,ni
   swap = array(i)
   array(i) = array(ni+1-i)
   array(ni+1-i) = swap
 enddo
end subroutine invert

!> Returns the global shape of h-point arrays
subroutine get_global_shape(domain, niglobal, njglobal)
  type(MOM_domain_type), intent(in)  :: domain   !< MOM domain from which to extract information
  integer,               intent(out) :: niglobal !< i-index global size of h-point arrays
  integer,               intent(out) :: njglobal !< j-index global size of h-point arrays

  niglobal = domain%niglobal
  njglobal = domain%njglobal
end subroutine get_global_shape

!> Get the array ranges in one dimension for the divisions of a global index space
subroutine compute_block_extent(isg, ieg, ndivs, ibegin, iend)
  integer,               intent(in)  :: isg    !< The starting index of the global index space
  integer,               intent(in)  :: ieg    !< The ending index of the global index space
  integer,               intent(in)  :: ndivs  !< The number of divisions
  integer, dimension(:), intent(out) :: ibegin !< The starting index of each division
  integer, dimension(:), intent(out) :: iend   !< The ending index of each division

  call mpp_compute_block_extent(isg, ieg, ndivs, ibegin, iend)
end subroutine compute_block_extent

!> Broadcast a 2-d domain from the root PE to the other PEs
subroutine broadcast_domain(domain)
  type(domain2d),  intent(inout) :: domain !< The domain2d type that will be shared across PEs.

  call mpp_broadcast_domain(domain)
end subroutine broadcast_domain

!> Broadcast an entire 2-d array from the root processor to all others.
subroutine global_field(domain, local, global)
  type(domain2d),       intent(inout) :: domain !< The domain2d type that describes the decomposition
  real, dimension(:,:), intent(in)    :: local  !< The portion of the array on the local PE
  real, dimension(:,:), intent(out)   :: global !< The whole global array

  call mpp_global_field(domain, local, global)
end subroutine global_field

!> Returns arrays of the i- and j- sizes of the h-point computational domains for each
!! element of the grid layout.  Any input values in the extent arrays are discarded, so
!! they are effectively intent out despite their declared intent of inout.
subroutine get_layout_extents(Domain, extent_i, extent_j)
  type(MOM_domain_type), intent(in)  :: domain !< MOM domain from which to extract information
  integer, dimension(:), allocatable, intent(inout) :: extent_i  !< The number of points in the
                                               !! i-direction in each i-row of the layout
  integer, dimension(:), allocatable, intent(inout) :: extent_j  !< The number of points in the
                                               !! j-direction in each j-row of the layout

  if (allocated(extent_i)) deallocate(extent_i)
  if (allocated(extent_j)) deallocate(extent_j)
  allocate(extent_i(domain%layout(1))) ; extent_i(:) = 0
  allocate(extent_j(domain%layout(2))) ; extent_j(:) = 0
  call mpp_get_domain_extents(domain%mpp_domain, extent_i, extent_j)
end subroutine get_layout_extents

end module MOM_domain_infra
