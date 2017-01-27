module MOM_domains

!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

use MOM_coms, only : PE_here, root_PE, num_PEs, MOM_infra_init, MOM_infra_end
use MOM_coms, only : broadcast, sum_across_PEs, min_across_PEs, max_across_PEs
use MOM_error_handler, only : MOM_error, MOM_mesg, NOTE, WARNING, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version
use MOM_file_parser, only : param_file_type
use MOM_string_functions, only : slasher

use mpp_domains_mod, only : mpp_define_layout, mpp_get_boundary
use mpp_domains_mod, only : MOM_define_io_domain => mpp_define_io_domain
use mpp_domains_mod, only : MOM_define_domain => mpp_define_domains
use mpp_domains_mod, only : domain2D, domain1D, mpp_get_data_domain
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_global_domain
use mpp_domains_mod, only : global_field_sum => mpp_global_sum
use mpp_domains_mod, only : mpp_update_domains, CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
use mpp_domains_mod, only : mpp_start_update_domains, mpp_complete_update_domains
use mpp_domains_mod, only : mpp_create_group_update, mpp_do_group_update
use mpp_domains_mod, only : group_pass_type => mpp_group_update_type
use mpp_domains_mod, only : mpp_reset_group_update_field
use mpp_domains_mod, only : mpp_group_update_initialized
use mpp_domains_mod, only : mpp_start_group_update, mpp_complete_group_update
use mpp_domains_mod, only : compute_block_extent => mpp_compute_block_extent
use mpp_parameter_mod, only : AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR, BITWISE_EXACT_SUM, CORNER
use mpp_parameter_mod, only : To_East => WUPDATE, To_West => EUPDATE
use mpp_parameter_mod, only : To_North => SUPDATE, To_South => NUPDATE
use fms_io_mod,        only : file_exist, parse_mask_table

implicit none ; private

! #include <MOM_memory.h>

public :: MOM_domains_init, MOM_infra_init, MOM_infra_end, get_domain_extent
public :: MOM_define_domain, MOM_define_io_domain, clone_MOM_domain
public :: pass_var, pass_vector, broadcast, PE_here, root_PE, num_PEs
public :: pass_var_start, pass_var_complete, fill_symmetric_edges
public :: pass_vector_start, pass_vector_complete
public :: global_field_sum, sum_across_PEs, min_across_PEs, max_across_PEs
public :: AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR, BITWISE_EXACT_SUM, CORNER
public :: To_East, To_West, To_North, To_South, To_All
public :: create_group_pass, do_group_pass, group_pass_type
public :: start_group_pass, complete_group_pass
public :: compute_block_extent

interface pass_var
  module procedure pass_var_3d, pass_var_2d
end interface pass_var

interface pass_vector
  module procedure pass_vector_3d, pass_vector_2d
end interface pass_vector

interface pass_var_start
  module procedure pass_var_start_3d, pass_var_start_2d
end interface pass_var_start

interface pass_var_complete
  module procedure pass_var_complete_3d, pass_var_complete_2d
end interface pass_var_complete

interface pass_vector_start
  module procedure pass_vector_start_3d, pass_vector_start_2d
end interface pass_vector_start

interface pass_vector_complete
  module procedure pass_vector_complete_3d, pass_vector_complete_2d
end interface pass_vector_complete

interface create_group_pass
  module procedure create_var_group_pass_2d
  module procedure create_var_group_pass_3d
  module procedure create_vector_group_pass_2d
  module procedure create_vector_group_pass_3d
end interface create_group_pass

interface fill_symmetric_edges
  module procedure fill_vector_symmetric_edges_2d !, fill_vector_symmetric_edges_3d
!   module procedure fill_scalar_symmetric_edges_2d, fill_scalar_symmetric_edges_3d
end interface fill_symmetric_edges

interface clone_MOM_domain
  module procedure clone_MD_to_MD, clone_MD_to_d2D
end interface clone_MOM_domain

type, public :: MOM_domain_type
  type(domain2D), pointer :: mpp_domain => NULL() ! The domain with halos on
                                        ! this processor, centered at h points.
  integer :: niglobal, njglobal         ! The total horizontal domain sizes.
  integer :: nihalo, njhalo             ! The X- and Y- halo sizes in memory.
  logical :: symmetric                  ! True if symmetric memory is used with
                                        ! this domain.
  logical :: nonblocking_updates        ! If true, non-blocking halo updates are
                                        ! allowed.  The default is .false. (for now).
  integer :: layout(2), io_layout(2)    ! Saved data for sake of constructing
  integer :: X_FLAGS, Y_FLAGS           ! new domains of different resolution.
  logical :: use_io_layout              ! True if an I/O layout is available.
  logical, pointer :: maskmap(:,:) => NULL() ! A pointer to an array indicating
                                ! which logical processors are actually used for
                                ! the ocean code. The other logical processors
                                ! would be all land points and are not assigned
                                ! to actual processors. This need not be
                                ! assigned if all logical processors are used.
end type MOM_domain_type

integer, parameter :: To_All = To_East + To_West + To_North + To_South

contains

subroutine pass_var_3d(array, MOM_dom, sideflag, complete, position)
  real, dimension(:,:,:), intent(inout) :: array
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: sideflag
  logical,      optional, intent(in)    :: complete
  integer,      optional, intent(in)    :: position
! Arguments: array - The array which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      sideflag - An optional integer indicating which directions the
!                       data should be sent.  It is TO_ALL or the sum of any of
!                       TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH.  For example,
!                       TO_EAST sends the data to the processor to the east, so
!                       the halos on the western side are filled.  TO_ALL is
!                       the default if sideflag is omitted.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be completed before progress resumes.  Omitting
!                       complete is the same as setting complete to .true.
!  (in)      position - An optional argument indicating the position.  This is
!                       usally CORNER, but is CENTER by default.
  integer :: dirflag
  logical :: block_til_complete

  dirflag = To_All ! 60
  if (PRESENT(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif
  block_til_complete = .true.
  if (present(complete)) block_til_complete = complete

  call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                          complete=block_til_complete, position=position)

end subroutine pass_var_3d


subroutine pass_var_2d(array, MOM_dom, sideflag, complete, position)
  real, dimension(:,:),  intent(inout) :: array
  type(MOM_domain_type), intent(inout) :: MOM_dom
  integer,     optional, intent(in)    :: sideflag
  logical,     optional, intent(in)    :: complete
  integer,     optional, intent(in)    :: position
! Arguments: array - The array which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      sideflag - An optional integer indicating which directions the
!                       data should be sent.  It is TO_ALL or the sum of any of
!                       TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH.  For example,
!                       TO_EAST sends the data to the processor to the east, so
!                       the halos on the western side are filled.  TO_ALL is
!                       the default if sideflag is omitted.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be completed before progress resumes.  Omitting
!                       complete is the same as setting complete to .true.
!  (in)      position - An optional argument indicating the position.  This is
!                       usally CORNER, but is CENTER by default.

  integer :: dirflag
  logical :: block_til_complete

  dirflag = To_All ! 60
  if (PRESENT(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif
  block_til_complete = .true.
  if (present(complete)) block_til_complete = complete

  call mpp_update_domains(array, MOM_dom%mpp_domain, flags=dirflag, &
                        complete=block_til_complete, position=position)

end subroutine pass_var_2d

function pass_var_start_2d(array, MOM_dom, sideflag, position, complete)
  real, dimension(:,:),   intent(inout) :: array
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: sideflag
  integer,      optional, intent(in)    :: position
  logical,      optional, intent(in)    :: complete
  integer :: pass_var_start_2d
! Arguments: array - The array which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      sideflag - An optional integer indicating which directions the
!                       data should be sent.  It is TO_ALL or the sum of any of
!                       TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH.  For example,
!                       TO_EAST sends the data to the processor to the east, so
!                       the halos on the western side are filled.  TO_ALL is
!                       the default if sideflag is omitted.
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second
!                       pass_..._start call.  Omitting complete is the same as
!                       setting complete to .true.
!  (return value) - The integer index for this update.
  integer :: dirflag

  dirflag = To_All ! 60
  if (PRESENT(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  pass_var_start_2d = mpp_start_update_domains(array, MOM_dom%mpp_domain, &
                          flags=dirflag, position=position)
end function pass_var_start_2d

function pass_var_start_3d(array, MOM_dom, sideflag, position, complete)
  real, dimension(:,:,:), intent(inout) :: array
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: sideflag
  integer,      optional, intent(in)    :: position
  logical,      optional, intent(in)    :: complete
  integer                               :: pass_var_start_3d
! Arguments: array - The array which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      sideflag - An optional integer indicating which directions the
!                       data should be sent.  It is TO_ALL or the sum of any of
!                       TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH.  For example,
!                       TO_EAST sends the data to the processor to the east, so
!                       the halos on the western side are filled.  TO_ALL is
!                       the default if sideflag is omitted.
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second
!                       pass_..._start call.  Omitting complete is the same as
!                       setting complete to .true.
!  (return value) - The integer index for this update.
  integer :: dirflag

  dirflag = To_All ! 60
  if (PRESENT(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  pass_var_start_3d = mpp_start_update_domains(array, MOM_dom%mpp_domain, &
                          flags=dirflag, position=position)
end function pass_var_start_3d

subroutine pass_var_complete_2d(id_update, array, MOM_dom, sideflag, position)
  integer,                intent(in)    :: id_update
  real, dimension(:,:),   intent(inout) :: array
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: sideflag
  integer,      optional, intent(in)    :: position
! Arguments: id_update - The integer id of this update which has been returned
!                        from a previous call to pass_var_start.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      sideflag - An optional integer indicating which directions the
!                       data should be sent.  It is TO_ALL or the sum of any of
!                       TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH.  For example,
!                       TO_EAST sends the data to the processor to the east, so
!                       the halos on the western side are filled.  TO_ALL is
!                       the default if sideflag is omitted.
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
  integer :: dirflag

  dirflag = To_All ! 60
  if (PRESENT(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  call mpp_complete_update_domains(id_update, array, MOM_dom%mpp_domain, &
                                   flags=dirflag, position=position)
end subroutine pass_var_complete_2d

subroutine pass_var_complete_3d(id_update, array, MOM_dom, sideflag, position)
  integer,                intent(in)    :: id_update
  real, dimension(:,:,:), intent(inout) :: array
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: sideflag
  integer,      optional, intent(in)    :: position
! Arguments: id_update - The integer id of this update which has been returned
!                        from a previous call to pass_var_start.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      sideflag - An optional integer indicating which directions the
!                       data should be sent.  It is TO_ALL or the sum of any of
!                       TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH.  For example,
!                       TO_EAST sends the data to the processor to the east, so
!                       the halos on the western side are filled.  TO_ALL is
!                       the default if sideflag is omitted.
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
  integer :: dirflag

  dirflag = To_All ! 60
  if (PRESENT(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  call mpp_complete_update_domains(id_update, array, MOM_dom%mpp_domain, &
                                   flags=dirflag, position=position)
end subroutine pass_var_complete_3d


subroutine pass_vector_2d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete)
  real, dimension(:,:),  intent(inout) :: u_cmpt, v_cmpt
  type(MOM_domain_type), intent(inout) :: MOM_dom
  integer,     optional, intent(in)    :: direction
  integer,     optional, intent(in)    :: stagger
  logical,     optional, intent(in)    :: complete
! Arguments: u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      direction - An optional integer indicating which directions the
!                        data should be sent.  It is TO_ALL or the sum of any of
!                        TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH, possibly
!                        plus SCALAR_PAIR if these are paired non-directional
!                        scalars discretized at the typical vector component
!                        locations.  For example, TO_EAST sends the data to the
!                        processor to the east, so the halos on the western
!                        side are filled.  TO_ALL is the default if omitted.
!  (in)      stagger - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      or CGRID_NE, indicating where the two components of the
!                      vector are discretized.  Omitting stagger is the same as
!                      setting it to CGRID_NE.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be completed before progress resumes.  Omitting
!                       complete is the same as setting complete to .true.

  integer :: stagger_local
  integer :: dirflag
  logical :: block_til_complete

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (PRESENT(direction)) then ; if (direction > 0) dirflag = direction ; endif
  block_til_complete = .true.
  if (present(complete)) block_til_complete = complete

  call mpp_update_domains(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                          gridtype=stagger_local, complete = block_til_complete)
end subroutine pass_vector_2d

subroutine fill_vector_symmetric_edges_2d(u_cmpt, v_cmpt, MOM_dom, stagger, scalar)
  real, dimension(:,:),  intent(inout) :: u_cmpt, v_cmpt
  type(MOM_domain_type), intent(inout) :: MOM_dom
  integer,     optional, intent(in)    :: stagger
  logical,     optional, intent(in)    :: scalar
! Arguments: u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      stagger - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      or CGRID_NE, indicating where the two components of the
!                      vector are discretized.  Omitting stagger is the same as
!                      setting it to CGRID_NE.
!  (in)      scalar -  An optional argument indicating whether

  integer :: stagger_local
  integer :: dirflag
  integer :: i, j, isc, iec, jsc, jec, isd, ied, jsd, jed, IscB, IecB, JscB, JecB
  real, allocatable, dimension(:) :: sbuff_x, sbuff_y, wbuff_x, wbuff_y
  logical :: block_til_complete

  if (.not. MOM_dom%symmetric) return

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
  if (PRESENT(scalar)) then ; if (scalar) dirflag = To_All+SCALAR_PAIR ; endif

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

end subroutine fill_vector_symmetric_edges_2d

subroutine pass_vector_3d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete)
  real, dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: direction
  integer,      optional, intent(in)    :: stagger
  logical,      optional, intent(in)    :: complete
! Arguments: u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      direction - An optional integer indicating which directions the
!                        data should be sent.  It is TO_ALL or the sum of any of
!                        TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH, possibly
!                        plus SCALAR_PAIR if these are paired non-directional
!                        scalars discretized at the typical vector component
!                        locations.  For example, TO_EAST sends the data to the
!                        processor to the east, so the halos on the western
!                        side are filled.  TO_ALL is the default if omitted.
!  (in)      stagger - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      or CGRID_NE, indicating where the two components of the
!                      vector are discretized.  Omitting stagger is the same as
!                      setting it to CGRID_NE.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be completed before progress resumes.  Omitting
!                       complete is the same as setting complete to .true.

  integer :: stagger_local
  integer :: dirflag
  logical :: block_til_complete

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (PRESENT(direction)) then ; if (direction > 0) dirflag = direction ; endif
  block_til_complete = .true.
  if (present(complete)) block_til_complete = complete

  call mpp_update_domains(u_cmpt, v_cmpt, MOM_dom%mpp_domain, flags=dirflag, &
                          gridtype=stagger_local, complete = block_til_complete)
end subroutine pass_vector_3d

function pass_vector_start_2d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete)
  real, dimension(:,:),   intent(inout) :: u_cmpt, v_cmpt
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: direction
  integer,      optional, intent(in)    :: stagger
  logical,      optional, intent(in)    :: complete
  integer                               :: pass_vector_start_2d
! Arguments: u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      direction - An optional integer indicating which directions the
!                        data should be sent.  It is TO_ALL or the sum of any of
!                        TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH, possibly
!                        plus SCALAR_PAIR if these are paired non-directional
!                        scalars discretized at the typical vector component
!                        locations.  For example, TO_EAST sends the data to the
!                        processor to the east, so the halos on the western
!                        side are filled.  TO_ALL is the default if omitted.
!  (in)      stagger - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      or CGRID_NE, indicating where the two components of the
!                      vector are discretized.  Omitting stagger is the same as
!                      setting it to CGRID_NE.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second
!                       pass_..._start call.  Omitting complete is the same as
!                       setting complete to .true.
!  (return value) - The integer index for this update.
  integer :: stagger_local
  integer :: dirflag

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (PRESENT(direction)) then ; if (direction > 0) dirflag = direction ; endif

  pass_vector_start_2d = mpp_start_update_domains(u_cmpt, v_cmpt, &
      MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)

end function pass_vector_start_2d

function pass_vector_start_3d(u_cmpt, v_cmpt, MOM_dom, direction, stagger, complete)
  real, dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: direction
  integer,      optional, intent(in)    :: stagger
  logical,      optional, intent(in)    :: complete
  integer                               :: pass_vector_start_3d
! Arguments: u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      direction - An optional integer indicating which directions the
!                        data should be sent.  It is TO_ALL or the sum of any of
!                        TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH, possibly
!                        plus SCALAR_PAIR if these are paired non-directional
!                        scalars discretized at the typical vector component
!                        locations.  For example, TO_EAST sends the data to the
!                        processor to the east, so the halos on the western
!                        side are filled.  TO_ALL is the default if omitted.
!  (in)      stagger - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      or CGRID_NE, indicating where the two components of the
!                      vector are discretized.  Omitting stagger is the same as
!                      setting it to CGRID_NE.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second
!                       pass_..._start call.  Omitting complete is the same as
!                       setting complete to .true.
!  (return value) - The integer index for this update.
  integer :: stagger_local
  integer :: dirflag

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (PRESENT(direction)) then ; if (direction > 0) dirflag = direction ; endif

  pass_vector_start_3d = mpp_start_update_domains(u_cmpt, v_cmpt, &
      MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)

end function pass_vector_start_3d

subroutine pass_vector_complete_2d(id_update, u_cmpt, v_cmpt, MOM_dom, direction, stagger)
  integer,                intent(in)    :: id_update
  real, dimension(:,:),   intent(inout) :: u_cmpt, v_cmpt
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: direction
  integer,      optional, intent(in)    :: stagger
! Arguments: id_update - The integer id of this update which has been returned
!                        from a previous call to pass_var_start.
!  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      direction - An optional integer indicating which directions the
!                        data should be sent.  It is TO_ALL or the sum of any of
!                        TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH, possibly
!                        plus SCALAR_PAIR if these are paired non-directional
!                        scalars discretized at the typical vector component
!                        locations.  For example, TO_EAST sends the data to the
!                        processor to the east, so the halos on the western
!                        side are filled.  TO_ALL is the default if omitted.
!  (in)      stagger - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      or CGRID_NE, indicating where the two components of the
!                      vector are discretized.  Omitting stagger is the same as
!                      setting it to CGRID_NE.
!  (return value) - The integer index for this update.
  integer :: stagger_local
  integer :: dirflag

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (PRESENT(direction)) then ; if (direction > 0) dirflag = direction ; endif

  call mpp_complete_update_domains(id_update, u_cmpt, v_cmpt, &
           MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)

end subroutine pass_vector_complete_2d

subroutine pass_vector_complete_3d(id_update, u_cmpt, v_cmpt, MOM_dom, direction, stagger)
  integer,                intent(in)    :: id_update
  real, dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: direction
  integer,      optional, intent(in)    :: stagger
! Arguments: id_update - The integer id of this update which has been returned
!                        from a previous call to pass_var_start.
!  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      direction - An optional integer indicating which directions the
!                        data should be sent.  It is TO_ALL or the sum of any of
!                        TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH, possibly
!                        plus SCALAR_PAIR if these are paired non-directional
!                        scalars discretized at the typical vector component
!                        locations.  For example, TO_EAST sends the data to the
!                        processor to the east, so the halos on the western
!                        side are filled.  TO_ALL is the default if omitted.
!  (in)      stagger - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      or CGRID_NE, indicating where the two components of the
!                      vector are discretized.  Omitting stagger is the same as
!                      setting it to CGRID_NE.
!  (return value) - The integer index for this update.
  integer :: stagger_local
  integer :: dirflag

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (PRESENT(direction)) then ; if (direction > 0) dirflag = direction ; endif

  call mpp_complete_update_domains(id_update, u_cmpt, v_cmpt, &
           MOM_dom%mpp_domain, flags=dirflag, gridtype=stagger_local)

end subroutine pass_vector_complete_3d

subroutine create_var_group_pass_2d(group, array, MOM_dom, sideflag, position)
  type(group_pass_type),  intent(inout) :: group
  real, dimension(:,:),   intent(inout) :: array
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: sideflag
  integer,      optional, intent(in)    :: position
! Arguments:
!  (inout)   group - The data type that store information for group update.
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      sideflag - An optional integer indicating which directions the
!                       data should be sent.  It is TO_ALL or the sum of any of
!                       TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH.  For example,
!                       TO_EAST sends the data to the processor to the east, so
!                       the halos on the western side are filled.  TO_ALL is
!                       the default if sideflag is omitted.
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
  integer :: dirflag

  dirflag = To_All ! 60
  if (PRESENT(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  else
    call mpp_create_group_update(group, array, MOM_dom%mpp_domain, flags=dirflag, &
                                 position=position)
  endif

end subroutine create_var_group_pass_2d

subroutine create_var_group_pass_3d(group, array, MOM_dom, sideflag, position)
  type(group_pass_type),  intent(inout) :: group
  real, dimension(:,:,:), intent(inout) :: array
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: sideflag
  integer,      optional, intent(in)    :: position
! Arguments:
!  (inout)   group - The data type that store information for group update.
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      sideflag - An optional integer indicating which directions the
!                       data should be sent.  It is TO_ALL or the sum of any of
!                       TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH.  For example,
!                       TO_EAST sends the data to the processor to the east, so
!                       the halos on the western side are filled.  TO_ALL is
!                       the default if sideflag is omitted.
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
  integer :: dirflag

  dirflag = To_All ! 60
  if (PRESENT(sideflag)) then ; if (sideflag > 0) dirflag = sideflag ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  else
    call mpp_create_group_update(group, array, MOM_dom%mpp_domain, flags=dirflag, &
                                 position=position)
  endif

end subroutine create_var_group_pass_3d


subroutine create_vector_group_pass_2d(group, u_cmpt, v_cmpt, MOM_dom, direction, stagger)
  type(group_pass_type),  intent(inout) :: group
  real, dimension(:,:),   intent(inout) :: u_cmpt, v_cmpt
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: direction
  integer,      optional, intent(in)    :: stagger
! Arguments:
!  (inout)   group - The data type that store information for group update.
!                    This data will be used in do_group_pass.
!  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      direction - An optional integer indicating which directions the
!                        data should be sent.  It is TO_ALL or the sum of any of
!                        TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH, possibly
!                        plus SCALAR_PAIR if these are paired non-directional
!                        scalars discretized at the typical vector component
!                        locations.  For example, TO_EAST sends the data to the
!                        processor to the east, so the halos on the western
!                        side are filled.  TO_ALL is the default if omitted.
!  (in)      stagger - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      or CGRID_NE, indicating where the two components of the
!                      vector are discretized.  Omitting stagger is the same as
!                      setting it to CGRID_NE.
  integer :: stagger_local
  integer :: dirflag

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (PRESENT(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,u_cmpt, v_cmpt)
  else
    call mpp_create_group_update(group, u_cmpt, v_cmpt, MOM_dom%mpp_domain, &
            flags=dirflag, gridtype=stagger_local)
  endif

end subroutine create_vector_group_pass_2d

subroutine create_vector_group_pass_3d(group, u_cmpt, v_cmpt, MOM_dom, direction, stagger)
  type(group_pass_type),  intent(inout) :: group
  real, dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt
  type(MOM_domain_type),  intent(inout) :: MOM_dom
  integer,      optional, intent(in)    :: direction
  integer,      optional, intent(in)    :: stagger
! Arguments:
!  (inout)   group - The data type that store information for group update.
!                    This data will be used in do_group_pass.
!  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.
!  (in)      direction - An optional integer indicating which directions the
!                        data should be sent.  It is TO_ALL or the sum of any of
!                        TO_EAST, TO_WEST, TO_NORTH, and TO_SOUTH, possibly
!                        plus SCALAR_PAIR if these are paired non-directional
!                        scalars discretized at the typical vector component
!                        locations.  For example, TO_EAST sends the data to the
!                        processor to the east, so the halos on the western
!                        side are filled.  TO_ALL is the default if omitted.
!  (in)      stagger - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      or CGRID_NE, indicating where the two components of the
!                      vector are discretized.  Omitting stagger is the same as
!                      setting it to CGRID_NE.

  integer :: stagger_local
  integer :: dirflag

  stagger_local = CGRID_NE ! Default value for type of grid
  if (present(stagger)) stagger_local = stagger

  dirflag = To_All ! 60
  if (PRESENT(direction)) then ; if (direction > 0) dirflag = direction ; endif

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,u_cmpt, v_cmpt)
  else
    call mpp_create_group_update(group, u_cmpt, v_cmpt, MOM_dom%mpp_domain, &
            flags=dirflag, gridtype=stagger_local)
  endif

end subroutine create_vector_group_pass_3d

subroutine do_group_pass(group, MOM_dom)
  type(group_pass_type), intent(inout) :: group
  type(MOM_domain_type), intent(inout) :: MOM_dom
  real                                 :: d_type

! Arguments:
!  (inout)   group - The data type that store information for group update.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.

  call mpp_do_group_update(group, MOM_dom%mpp_domain, d_type)

end subroutine do_group_pass

subroutine start_group_pass(group, MOM_dom)
  type(group_pass_type), intent(inout) :: group
  type(MOM_domain_type), intent(inout) :: MOM_dom
  real                                 :: d_type

! Arguments:
!  (inout)   group - The data type that store information for group update.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.

  call mpp_start_group_update(group, MOM_dom%mpp_domain, d_type)

end subroutine start_group_pass

subroutine complete_group_pass(group, MOM_dom)
  type(group_pass_type), intent(inout) :: group
  type(MOM_domain_type), intent(inout) :: MOM_dom
  real                                 :: d_type

! Arguments:
!  (inout)   group - The data type that store information for group update.
!  (in)      MOM_dom - The MOM_domain_type containing the mpp_domain needed to
!                      determine where data should be sent.

  call mpp_complete_group_update(group, MOM_dom%mpp_domain, d_type)

end subroutine complete_group_pass


subroutine MOM_domains_init(MOM_dom, param_file, symmetric, static_memory, &
                            NIHALO, NJHALO, NIGLOBAL, NJGLOBAL, NIPROC, NJPROC, &
                            min_halo, domain_name, include_name, param_suffix)
  type(MOM_domain_type),           pointer       :: MOM_dom
  type(param_file_type),           intent(in)    :: param_file
  logical, optional,               intent(in)    :: symmetric
  logical, optional,               intent(in)    :: static_memory
  integer, optional,               intent(in)    :: NIHALO, NJHALO
  integer, optional,               intent(in)    :: NIGLOBAL, NJGLOBAL
  integer, optional,               intent(in)    :: NIPROC, NJPROC
  integer, dimension(2), optional, intent(inout) :: min_halo
  character(len=*),      optional, intent(in)    :: domain_name
  character(len=*),      optional, intent(in)    :: include_name
  character(len=*),      optional, intent(in)    :: param_suffix


! Arguments: MOM_dom - A pointer to the MOM_domain_type being defined here.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in,opt)  symmetric - If present, this specifies whether this domain
!                        is symmetric, regardless of whether the macro
!                        SYMMETRIC_MEMORY_ is defined.
!  (in,opt)  static_memory - If present and true, this domain type is set up for
!                            static memory and error checking of various input
!                            values is performed against those in the input file.
!  (in,opt)  NIHALO, NJHALO - Default halo sizes, required with static memory.
!  (in,opt)  NIGLOBAL, NJGLOBAL - Total domain sizes, required with static memory.
!  (in,opt)  NIPROC, NJPROC - Processor counts, required with static memory.
!  (in,opt)  min_halo - If present, this sets the minimum halo size for this
!                       domain in the x- and y- directions, and returns the
!                       actual halo size used.
!  (in,opt)  domain_name - A name for this domain, "MOM" if missing.
!  (in,opt)  include_name - A name for model's include file, "MOM_memory.h" if missing.
!  (in,opt)  param_suffix - A suffix to apply to layout-specific parameters.

  integer, dimension(2) :: layout = (/ 1, 1 /)
  integer, dimension(2) :: io_layout = (/ 0, 0 /)
  integer, dimension(4) :: global_indices
!$ integer :: ocean_nthreads       ! Number of Openmp threads
!$ integer :: get_cpu_affinity, omp_get_thread_num, omp_get_num_threads
!$ integer :: omp_cores_per_node, adder, base_cpu
!$ logical :: ocean_omp_hyper_thread
  integer :: nihalo_dflt, njhalo_dflt
  integer :: pe, proc_used
  integer :: X_FLAGS, Y_FLAGS
  logical :: reentrant_x, reentrant_y, tripolar_N, is_static
  logical            :: mask_table_exists
  character(len=128) :: mask_table, inputdir
  character(len=64)  :: dom_name, inc_nm
  character(len=200) :: mesg

  integer :: xsiz, ysiz, nip_parsed, njp_parsed
  integer :: isc,iec,jsc,jec ! The bounding indices of the computational domain.
  character(len=8) :: char_xsiz, char_ysiz, char_niglobal, char_njglobal
  character(len=40) :: nihalo_nm, njhalo_nm, layout_nm, io_layout_nm, masktable_nm
  character(len=40) :: niproc_nm, njproc_nm

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl ! This module's name.

  if (.not.associated(MOM_dom)) then
    allocate(MOM_dom)
    allocate(MOM_dom%mpp_domain)
  endif

  pe = PE_here()
  proc_used = num_PEs()

  mdl = "MOM_domains"

  MOM_dom%symmetric = .true.
  if (present(symmetric)) then ; MOM_dom%symmetric = symmetric ; endif
  if (present(min_halo)) mdl = trim(mdl)//" min_halo"

  dom_name = "MOM" ; inc_nm = "MOM_memory.h"
  if (present(domain_name)) dom_name = trim(domain_name)
  if (present(include_name)) inc_nm = trim(include_name)

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
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "REENTRANT_X", reentrant_x, &
                 "If true, the domain is zonally reentrant.", default=.true.)
  call get_param(param_file, mdl, "REENTRANT_Y", reentrant_y, &
                 "If true, the domain is meridionally reentrant.", &
                 default=.false.)
  call get_param(param_file, mdl, "TRIPOLAR_N", tripolar_N, &
                 "Use tripolar connectivity at the northern edge of the \n"//&
                 "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", &
                 default=.false.)

#ifndef NOT_SET_AFFINITY
!$ call get_param(param_file, mdl, "OCEAN_OMP_THREADS", ocean_nthreads, &
!$            "The number of OpenMP threads that MOM6 will use.", &
!$            default = 1, layoutParam=.true.)
!$ call get_param(param_file, mdl, "OCEAN_OMP_HYPER_THREAD", ocean_omp_hyper_thread, &
!$            "If True, use hyper-threading.", default = .false., layoutParam=.true.)
!$ if (ocean_omp_hyper_thread) then
!$   call get_param(param_file, mdl, "OMP_CORES_PER_NODE", omp_cores_per_node, &
!$            "Number of cores per node needed for hyper-threading.", &
!$            fail_if_missing=.true., layoutParam=.true.)
!$ endif
!$ call omp_set_num_threads(ocean_nthreads)
!$OMP PARALLEL private(adder)
!$ base_cpu = get_cpu_affinity()
!$ if (ocean_omp_hyper_thread) then
!$   if (mod(omp_get_thread_num(),2) == 0) then
!$     adder = omp_get_thread_num()/2
!$   else
!$     adder = omp_cores_per_node + omp_get_thread_num()/2
!$   endif
!$ else
!$   adder = omp_get_thread_num()
!$ endif
!$ call set_cpu_affinity(base_cpu + adder)
!!$     write(6,*) " ocean  ", omp_get_num_threads(), get_cpu_affinity(), adder, omp_get_thread_num()
!!$     call flush(6)
!$OMP END PARALLEL
#endif

  call log_param(param_file, mdl, "!SYMMETRIC_MEMORY_", MOM_dom%symmetric, &
                 "If defined, the velocity point data domain includes \n"//&
                 "every face of the thickness points. In other words, \n"//&
                 "some arrays are larger than others, depending on where \n"//&
                 "they are on the staggered grid.  Also, the starting \n"//&
                 "index of the velocity-point arrays is usually 0, not 1. \n"//&
                 "This can only be set at compile time.",&
                 layoutParam=.true.)
  call get_param(param_file, mdl, "NONBLOCKING_UPDATES", MOM_dom%nonblocking_updates, &
                 "If true, non-blocking halo updates may be used.", &
                 default=.false., layoutParam=.true.)

  nihalo_dflt = 4 ; njhalo_dflt = 4
  if (present(NIHALO)) nihalo_dflt = NIHALO
  if (present(NJHALO)) njhalo_dflt = NJHALO

  call log_param(param_file, mdl, "!STATIC_MEMORY_", is_static, &
                 "If STATIC_MEMORY_ is defined, the principle variables \n"//&
                 "will have sizes that are statically determined at \n"//&
                 "compile time.  Otherwise the sizes are not determined \n"//&
                 "until run time. The STATIC option is substantially \n"//&
                 "faster, but does not allow the PE count to be changed \n"//&
                 "at run time.  This can only be set at compile time.",&
                 layoutParam=.true.)

  call get_param(param_file, mdl, trim(nihalo_nm), MOM_dom%nihalo, &
                 "The number of halo points on each side in the \n"//&
                 "x-direction.  With STATIC_MEMORY_ this is set as NIHALO_ \n"//&
                 "in "//trim(inc_nm)//" at compile time; without STATIC_MEMORY_ \n"//&
                 "the default is NIHALO_ in "//trim(inc_nm)//" (if defined) or 2.", &
                 default=4, static_value=nihalo_dflt, layoutParam=.true.)
  call get_param(param_file, mdl, trim(njhalo_nm), MOM_dom%njhalo, &
                 "The number of halo points on each side in the \n"//&
                 "y-direction.  With STATIC_MEMORY_ this is set as NJHALO_ \n"//&
                 "in "//trim(inc_nm)//" at compile time; without STATIC_MEMORY_ \n"//&
                 "the default is NJHALO_ in "//trim(inc_nm)//" (if defined) or 2.", &
                 default=4, static_value=njhalo_dflt, layoutParam=.true.)
  if (present(min_halo)) then
    MOM_dom%nihalo = max(MOM_dom%nihalo, min_halo(1))
    min_halo(1) = MOM_dom%nihalo
    MOM_dom%njhalo = max(MOM_dom%njhalo, min_halo(2))
    min_halo(2) = MOM_dom%njhalo
    call log_param(param_file, mdl, "!NIHALO min_halo", MOM_dom%nihalo, layoutParam=.true.)
    call log_param(param_file, mdl, "!NJHALO min_halo", MOM_dom%nihalo, layoutParam=.true.)
  endif
  if (is_static) then
    call get_param(param_file, mdl, "NIGLOBAL", MOM_dom%niglobal, &
                 "The total number of thickness grid points in the \n"//&
                 "x-direction in the physical domain. With STATIC_MEMORY_ \n"//&
                 "this is set in "//trim(inc_nm)//" at compile time.", &
                 static_value=NIGLOBAL)
    call get_param(param_file, mdl, "NJGLOBAL", MOM_dom%njglobal, &
                 "The total number of thickness grid points in the \n"//&
                 "y-direction in the physical domain. With STATIC_MEMORY_ \n"//&
                 "this is set in "//trim(inc_nm)//" at compile time.", &
                 static_value=NJGLOBAL)
    if (MOM_dom%niglobal /= NIGLOBAL) call MOM_error(FATAL,"MOM_domains_init: " // &
     "static mismatch for NIGLOBAL_ domain size. Header file does not match input namelist")
    if (MOM_dom%njglobal /= NJGLOBAL) call MOM_error(FATAL,"MOM_domains_init: " // &
     "static mismatch for NJGLOBAL_ domain size. Header file does not match input namelist")

    if (.not.present(min_halo)) then
      if (MOM_dom%nihalo /= NIHALO) call MOM_error(FATAL,"MOM_domains_init: " // &
             "static mismatch for "//trim(nihalo_nm)//" domain size")
      if (MOM_dom%njhalo /= NJHALO) call MOM_error(FATAL,"MOM_domains_init: " // &
             "static mismatch for "//trim(njhalo_nm)//" domain size")
    endif
  else
    call get_param(param_file, mdl, "NIGLOBAL", MOM_dom%niglobal, &
                 "The total number of thickness grid points in the \n"//&
                 "x-direction in the physical domain. With STATIC_MEMORY_ \n"//&
                 "this is set in "//trim(inc_nm)//" at compile time.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "NJGLOBAL", MOM_dom%njglobal, &
                 "The total number of thickness grid points in the \n"//&
                 "y-direction in the physical domain. With STATIC_MEMORY_ \n"//&
                 "this is set in "//trim(inc_nm)//" at compile time.", &
                 fail_if_missing=.true.)
  endif

  global_indices(1) = 1 ; global_indices(2) = MOM_dom%niglobal
  global_indices(3) = 1 ; global_indices(4) = MOM_dom%njglobal

  call get_param(param_file, mdl, "INPUTDIR", inputdir, do_not_log=.true., default=".")
  inputdir = slasher(inputdir)

  call get_param(param_file, mdl, trim(masktable_nm), mask_table, &
                 "A text file to specify n_mask, layout and mask_list. \n"//&
                 "This feature masks out processors that contain only land points. \n"//&
                 "The first line of mask_table is the number of regions to be masked out.\n"//&
                 "The second line is the layout of the model and must be \n"//&
                 "consistent with the actual model layout.\n"//&
                 "The following (n_mask) lines give the logical positions \n"//&
                 "of the processors that are masked out. The mask_table \n"//&
                 "can be created by tools like check_mask. The \n"//&
                 "following example of mask_table masks out 2 processors, \n"//&
                 "(1,2) and (3,6), out of the 24 in a 4x6 layout: \n"//&
                 " 2\n 4,6\n 1,2\n 3,6\n", default="MOM_mask_table", &
                 layoutParam=.true.)
  mask_table = trim(inputdir)//trim(mask_table)
  mask_table_exists = file_exist(mask_table)

  if (is_static) then
    layout(1) = NIPROC ; layout(2) = NJPROC
  else
    call get_param(param_file, mdl, trim(layout_nm), layout, &
                 "The processor layout to be used, or 0, 0 to automatically \n"//&
                 "set the layout based on the number of processors.", default=0, &
                 do_not_log=.true.)
    call get_param(param_file, mdl, trim(niproc_nm), nip_parsed, &
                 "The number of processors in the x-direction.", default=-1, &
                 do_not_log=.true.)
    call get_param(param_file, mdl, trim(njproc_nm), njp_parsed, &
                 "The number of processors in the y-direction.", default=-1, &
                 do_not_log=.true.)
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

    if ( layout(1)==0 .and. layout(2)==0 ) &
      call mpp_define_layout(global_indices, proc_used, layout)
    if ( layout(1)/=0 .and. layout(2)==0 ) layout(2) = proc_used/layout(1)
    if ( layout(1)==0 .and. layout(2)/=0 ) layout(1) = proc_used/layout(2)

    if (layout(1)*layout(2) /= proc_used .and. (.not. mask_table_exists) ) then
      write(mesg,'("MOM_domains_init: The product of the two components of layout, ", &
            &      2i4,", is not the number of PEs used, ",i5,".")') &
            layout(1),layout(2),proc_used
      call MOM_error(FATAL, mesg)
    endif
  endif
  call log_param(param_file, mdl, trim(niproc_nm), layout(1), &
                 "The number of processors in the x-direction. With \n"//&
                 "STATIC_MEMORY_ this is set in "//trim(inc_nm)//" at compile time.",&
                 layoutParam=.true.)
  call log_param(param_file, mdl, trim(njproc_nm), layout(2), &
                 "The number of processors in the x-direction. With \n"//&
                 "STATIC_MEMORY_ this is set in "//trim(inc_nm)//" at compile time.",&
                 layoutParam=.true.)
  call log_param(param_file, mdl, trim(layout_nm), layout, &
                 "The processor layout that was acutally used.",&
                 layoutParam=.true.)

  ! Idiot check that fewer PEs than columns have been requested
  if (layout(1)*layout(2)>MOM_dom%niglobal*MOM_dom%njglobal)  then
    write(mesg,'(a,2(i5,x,a))') 'You requested to use',layout(1)*layout(2), &
      'PEs but there are only',MOM_dom%niglobal*MOM_dom%njglobal,'columns in the model'
    call MOM_error(FATAL, mesg)
  endif

  if (mask_table_exists) then
    call MOM_error(NOTE, 'MOM_domains_init: reading maskmap information from '//&
                         trim(mask_table))
    allocate(MOM_dom%maskmap(layout(1), layout(2)))
    call parse_mask_table(mask_table, MOM_dom%maskmap, dom_name)
  endif

  !   Set up the I/O lay-out, and check that it uses an even multiple of the
  ! number of PEs in each direction.
  io_layout(:) = (/ 1, 1 /)
  call get_param(param_file, mdl, trim(io_layout_nm), io_layout, &
                 "The processor layout to be used, or 0,0 to automatically \n"//&
                 "set the io_layout to be the same as the layout.", default=1, &
                 layoutParam=.true.)

  if (io_layout(1) < 0) then
    write(mesg,'("MOM_domains_init: IO_LAYOUT(1) = ",i4,".  Negative values "//&
         &"are not allowed in ")') io_layout(1)
    call MOM_error(FATAL, mesg//trim(IO_layout_nm))
  elseif (io_layout(1) > 0) then ; if (modulo(layout(1), io_layout(1)) /= 0) then
    write(mesg,'("MOM_domains_init: The x-direction I/O-layout, IO_LAYOUT(1)=",i4, &
         &", does not evenly divide the x-direction layout, NIPROC=,",i4,".")') &
          io_layout(1),layout(1)
    call MOM_error(FATAL, mesg)
  endif ; endif

  if (io_layout(2) < 0) then
    write(mesg,'("MOM_domains_init: IO_LAYOUT(2) = ",i4,".  Negative values "//&
         &"are not allowed in ")') io_layout(2)
    call MOM_error(FATAL, mesg//trim(IO_layout_nm))
  elseif (io_layout(2) /= 0) then ; if (modulo(layout(2), io_layout(2)) /= 0) then
    write(mesg,'("MOM_domains_init: The y-direction I/O-layout, IO_LAYOUT(2)=",i4, &
         &", does not evenly divide the y-direction layout, NJPROC=,",i4,".")') &
          io_layout(2),layout(2)
    call MOM_error(FATAL, mesg)
  endif ; endif

  if (io_layout(2) == 0) io_layout(2) = layout(2)
  if (io_layout(1) == 0) io_layout(1) = layout(1)

  X_FLAGS = 0 ; Y_FLAGS = 0
  if (reentrant_x) X_FLAGS = CYCLIC_GLOBAL_DOMAIN
  if (reentrant_y) Y_FLAGS = CYCLIC_GLOBAL_DOMAIN
  if (tripolar_N) then
    Y_FLAGS = FOLD_NORTH_EDGE
    if (reentrant_y) call MOM_error(FATAL,"MOM_domains: "// &
      "TRIPOLAR_N and REENTRANT_Y may not be defined together.")
  endif

  global_indices(1) = 1 ; global_indices(2) = MOM_dom%niglobal
  global_indices(3) = 1 ; global_indices(4) = MOM_dom%njglobal

  if (mask_table_exists) then
    call MOM_define_domain( global_indices, layout, MOM_dom%mpp_domain, &
                xflags=X_FLAGS, yflags=Y_FLAGS, &
                xhalo=MOM_dom%nihalo, yhalo=MOM_dom%njhalo, &
                symmetry = MOM_dom%symmetric, name=dom_name, &
                maskmap=MOM_dom%maskmap )
  else
    call MOM_define_domain( global_indices, layout, MOM_dom%mpp_domain, &
                xflags=X_FLAGS, yflags=Y_FLAGS, &
                xhalo=MOM_dom%nihalo, yhalo=MOM_dom%njhalo, &
                symmetry = MOM_dom%symmetric, name=dom_name)
  endif

  if ((io_layout(1) > 0) .and. (io_layout(2) > 0) .and. &
      (layout(1)*layout(2) > 1)) then
    call MOM_define_io_domain(MOM_dom%mpp_domain, io_layout)
  endif

! Save the extra data for creating other domains of different resolution that overlay this domain
  MOM_dom%X_FLAGS = X_FLAGS
  MOM_dom%Y_FLAGS = Y_FLAGS
  MOM_dom%layout = layout
  MOM_dom%io_layout = io_layout
  MOM_dom%use_io_layout = (io_layout(1) + io_layout(2) > 0)

  if (is_static) then
  !   A requirement of equal sized compute domains is necessary when STATIC_MEMORY_
  ! is used.
    call mpp_get_compute_domain(MOM_dom%mpp_domain,isc,iec,jsc,jec)
    xsiz = iec - isc + 1
    ysiz = jec - jsc + 1
    if (xsiz*NIPROC /= MOM_dom%niglobal .OR. ysiz*NJPROC /= MOM_dom%njglobal) then
       write( char_xsiz,'(i4)' ) NIPROC
       write( char_ysiz,'(i4)' ) NJPROC
       write( char_niglobal,'(i4)' ) MOM_dom%niglobal
       write( char_njglobal,'(i4)' ) MOM_dom%njglobal
       call MOM_error(WARNING,'MOM_domains: Processor decomposition (NIPROC_,NJPROC_) = (' &
           //trim(char_xsiz)//','//trim(char_ysiz)// &
           ') does not evenly divide size set by preprocessor macro ('&
           //trim(char_niglobal)//','//trim(char_njglobal)// '). ')
       call MOM_error(FATAL,'MOM_domains:  #undef STATIC_MEMORY_ in "//trim(inc_nm)//" to use &
           &dynamic allocation, or change processor decomposition to evenly divide the domain.')
    endif
  endif

end subroutine MOM_domains_init

subroutine clone_MD_to_MD(MD_in, MOM_dom, min_halo, halo_size, symmetric, &
                          domain_name)
  type(MOM_domain_type),           intent(in)    :: MD_in
  type(MOM_domain_type),           pointer       :: MOM_dom
  integer, dimension(2), optional, intent(inout) :: min_halo
  integer,               optional, intent(in)    :: halo_size
  logical,               optional, intent(in)    :: symmetric
  character(len=*),      optional, intent(in)    :: domain_name

  integer :: global_indices(4)
  logical :: mask_table_exists
  character(len=64) :: dom_name

  if (.not.associated(MOM_dom)) then
    allocate(MOM_dom)
    allocate(MOM_dom%mpp_domain)
  endif

! Save the extra data for creating other domains of different resolution that overlay this domain
  MOM_dom%niglobal = MD_in%niglobal ; MOM_dom%njglobal = MD_in%njglobal
  MOM_dom%nihalo = MD_in%nihalo ; MOM_dom%njhalo = MD_in%njhalo

  MOM_dom%symmetric = MD_in%symmetric
  MOM_dom%nonblocking_updates = MD_in%nonblocking_updates

  MOM_dom%X_FLAGS = MD_in%X_FLAGS ; MOM_dom%Y_FLAGS = MD_in%Y_FLAGS
  MOM_dom%layout(:) = MD_in%layout(:) ; MOM_dom%io_layout(:) = MD_in%io_layout(:)
  MOM_dom%use_io_layout = (MOM_dom%io_layout(1) + MOM_dom%io_layout(2) > 0)

  if (associated(MD_in%maskmap)) then
    mask_table_exists = .true.
    allocate(MOM_dom%maskmap(MOM_dom%layout(1), MOM_dom%layout(2)))
    MOM_dom%maskmap(:,:) = MD_in%maskmap(:,:)
  else
    mask_table_exists = .false.
  endif

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

  dom_name = "MOM"
  if (present(domain_name)) dom_name = trim(domain_name)

  global_indices(1) = 1 ; global_indices(2) = MOM_dom%niglobal
  global_indices(3) = 1 ; global_indices(4) = MOM_dom%njglobal
  if (mask_table_exists) then
    call MOM_define_domain( global_indices, MOM_dom%layout, MOM_dom%mpp_domain, &
                xflags=MOM_dom%X_FLAGS, yflags=MOM_dom%Y_FLAGS, &
                xhalo=MOM_dom%nihalo, yhalo=MOM_dom%njhalo, &
                symmetry = MOM_dom%symmetric, name=dom_name, &
                maskmap=MOM_dom%maskmap )
  else
    call MOM_define_domain( global_indices, MOM_dom%layout, MOM_dom%mpp_domain, &
                xflags=MOM_dom%X_FLAGS, yflags=MOM_dom%Y_FLAGS, &
                xhalo=MOM_dom%nihalo, yhalo=MOM_dom%njhalo, &
                symmetry = MOM_dom%symmetric, name=dom_name)
  endif

  if ((MOM_dom%io_layout(1) + MOM_dom%io_layout(2) > 0) .and. &
      (MOM_dom%layout(1)*MOM_dom%layout(2) > 1)) then
    call MOM_define_io_domain(MOM_dom%mpp_domain, MOM_dom%io_layout)
  endif

end subroutine clone_MD_to_MD

subroutine clone_MD_to_d2D(MD_in, mpp_domain, min_halo, halo_size, symmetric, &
                           domain_name)
  type(MOM_domain_type),           intent(in)    :: MD_in
  type(domain2d),                  intent(inout) :: mpp_domain
  integer, dimension(2), optional, intent(inout) :: min_halo
  integer,               optional, intent(in)    :: halo_size
  logical,               optional, intent(in)    :: symmetric
  character(len=*),      optional, intent(in)    :: domain_name

  integer :: global_indices(4), layout(2), io_layout(2)
  integer :: X_FLAGS, Y_FLAGS, niglobal, njglobal, nihalo, njhalo
  logical :: symmetric_dom
  character(len=64) :: dom_name

! Save the extra data for creating other domains of different resolution that overlay this domain
  niglobal = MD_in%niglobal ; njglobal = MD_in%njglobal
  nihalo = MD_in%nihalo ; njhalo = MD_in%njhalo

  symmetric_dom = MD_in%symmetric

  X_FLAGS = MD_in%X_FLAGS ; Y_FLAGS = MD_in%Y_FLAGS
  layout(:) = MD_in%layout(:) ; io_layout(:) = MD_in%io_layout(:)

  if (present(halo_size) .and. present(min_halo)) call MOM_error(FATAL, &
      "clone_MOM_domain can not have both halo_size and min_halo present.")

  if (present(min_halo)) then
    nihalo = max(nihalo, min_halo(1))
    njhalo = max(njhalo, min_halo(2))
    min_halo(1) = nihalo ; min_halo(2) = njhalo
  endif

  if (present(halo_size)) then
    nihalo = halo_size ; njhalo = halo_size
  endif

  if (present(symmetric)) then ; symmetric_dom = symmetric ; endif

  dom_name = "MOM"
  if (present(domain_name)) dom_name = trim(domain_name)

  global_indices(1) = 1 ; global_indices(2) = niglobal
  global_indices(3) = 1 ; global_indices(4) = njglobal
  if (associated(MD_in%maskmap)) then
    call MOM_define_domain( global_indices, layout, mpp_domain, &
                xflags=X_FLAGS, yflags=Y_FLAGS, &
                xhalo=nihalo, yhalo=njhalo, &
                symmetry = symmetric, name=dom_name, &
                maskmap=MD_in%maskmap )
  else
    call MOM_define_domain( global_indices, layout, mpp_domain, &
                xflags=X_FLAGS, yflags=Y_FLAGS, &
                xhalo=nihalo, yhalo=njhalo, &
                symmetry = symmetric, name=dom_name)
  endif

  if ((io_layout(1) + io_layout(2) > 0) .and. &
      (layout(1)*layout(2) > 1)) then
    call MOM_define_io_domain(mpp_domain, io_layout)
  endif

end subroutine clone_MD_to_d2D

subroutine get_domain_extent(Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, &
                             isg, ieg, jsg, jeg, idg_offset, jdg_offset, &
                             symmetric, local_indexing, index_offset)
  type(MOM_domain_type), intent(in) :: Domain
  integer, intent(out) :: isc, iec, jsc, jec
  integer, intent(out) :: isd, ied, jsd, jed
  integer, intent(out) :: isg, ieg, jsg, jeg
  integer, intent(out) :: idg_offset, jdg_offset
  logical, intent(out) :: symmetric
  logical, optional, intent(in) :: local_indexing
  integer, optional, intent(in) :: index_offset
! Arguments: Domain - The MOM_domain_type from which the indices are extracted.
!  (out)     isc, iec, jsc, jec - the start & end indices of the
!                                 computational domain.
!  (out)     isd, ied, jsd, jed - the start & end indices of the data domain.
!  (out)     isg, ieg, jsg, jeg - the start & end indices of the global domain.
!  (out)     idg_offset, jdg_offset - the offset between the corresponding
!                                     global and data index spaces.
!  (out)     symmetric - true if symmetric memory is used.
!  (in,opt)  local_indexing - if true, local tracer array indices start at 1, as
!                             in most MOM6 or GOLD code.
!  (in,opt)  index_offset - A fixed additional offset to all indices.  This can
!                           be useful for some types of debugging with dynamic
!                           memory allocation.

  integer :: ind_off
  logical :: local
  local = .true. ; if (present(local_indexing)) local = local_indexing
  ind_off = 0 ; if (present(index_offset)) ind_off = index_offset

  call mpp_get_compute_domain(Domain%mpp_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(Domain%mpp_domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(Domain%mpp_domain, isg, ieg, jsg, jeg)

  ! This code institutes the MOM convention that local array indices start at 1.
  if (local) then
    idg_offset = isd-1 ; jdg_offset = jsd-1
    isc = isc-isd+1 ; iec = iec-isd+1 ; jsc = jsc-jsd+1 ; jec = jec-jsd+1
    ied = ied-isd+1 ; jed = jed-jsd+1
    isd = 1 ; jsd = 1
  else
    idg_offset = 0 ; jdg_offset = 0
  endif
  if (ind_off /= 0) then
    idg_offset = idg_offset + ind_off ; jdg_offset = jdg_offset + ind_off
    isc = isc + ind_off ; iec = iec + ind_off
    jsc = jsc + ind_off ; jec = jec + ind_off
    isd = isd + ind_off ; ied = ied + ind_off
    jsd = jsd + ind_off ; jed = jed + ind_off
  endif
  symmetric = Domain%symmetric

end subroutine get_domain_extent

end module MOM_domains
