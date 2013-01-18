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

use MOM_error_handler, only : MOM_error, NOTE, WARNING, FATAL, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version
use MOM_file_parser, only : param_file_type

use MOM_coms, only : PE_here, root_PE, num_PEs, MOM_infra_init, MOM_infra_end
use MOM_coms, only : broadcast, sum_across_PEs, min_across_PEs, max_across_PEs
use mpp_domains_mod, only : mpp_define_layout
use mpp_domains_mod, only : MOM_define_io_domain => mpp_define_io_domain
use mpp_domains_mod, only : MOM_define_domain => mpp_define_domains
use mpp_domains_mod, only : domain2D, domain1D, mpp_get_data_domain
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_global_domain
use mpp_domains_mod, only : global_field_sum => mpp_global_sum
use mpp_domains_mod, only : mpp_update_domains, CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
use mpp_domains_mod, only : mpp_start_update_domains, mpp_complete_update_domains
use mpp_parameter_mod, only : AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR, BITWISE_EXACT_SUM, CORNER
use mpp_parameter_mod, only : To_East => WUPDATE, To_West => EUPDATE
use mpp_parameter_mod, only : To_North => SUPDATE, To_South => NUPDATE

implicit none ; private

#include <MOM_memory.h>

public :: MOM_domains_init, MOM_infra_init, MOM_infra_end, get_domain_extent
public :: MOM_define_domain, MOM_define_io_domain
public :: pass_var, pass_vector, broadcast, PE_here, root_PE, num_PEs
public :: pass_var_start, pass_var_complete
public :: pass_vector_start, pass_vector_complete
public :: global_field_sum, sum_across_PEs, min_across_PEs, max_across_PEs
public :: AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR, BITWISE_EXACT_SUM, CORNER
public :: To_East, To_West, To_North, To_South, To_All

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
  integer :: X_FLAGS,Y_FLAGS            ! new domains of different resolution.
  logical :: use_io_layout              ! True if an I/O layout is available.

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
  integer :: pass_var_start_2d
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
  integer                               :: pass_vector_start_2d
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
  integer                               :: pass_vector_start_2d
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

subroutine MOM_domains_init(MOM_dom, param_file, min_halo, symmetric)
  type(MOM_domain_type),           pointer       :: MOM_dom
  type(param_file_type),           intent(in)    :: param_file
  integer, dimension(2), optional, intent(inout) :: min_halo
  logical, optional,               intent(in)    :: symmetric
! Arguments: MOM_dom - A pointer to the MOM_domain_type being defined here.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in,opt)  min_halo - If present, this sets the minimum halo size for this
!                       domain in the x- and y- directions, and returns the
!                       actual halo size used.
!  (in,opt)  symmetric - If present, this specified whether this domain
!                        is symmetric, regardless of whether the macro
!                        SYMMETRIC_MEMORY_ is defined.

  integer, dimension(2) :: layout = (/ 1, 1 /)
  integer, dimension(2) :: io_layout = (/ 0, 0 /)
  integer, dimension(4) :: global_indices
  integer :: nihalo, njhalo, nihalo_dflt, njhalo_dflt
  integer :: pe, proc_used
  integer :: isc,iec,jsc,jec ! The bounding indices of the computational domain.
  integer :: X_FLAGS, Y_FLAGS
  integer :: i, xsiz, ysiz
  logical :: reentrant_x, reentrant_y, tripolar_N, is_static
  character(len=200) :: mesg
  character(len=8) :: char_xsiz, char_ysiz, char_niglobal, char_njglobal
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'
  character(len=40)  :: mod ! This module's name.

  if (.not.associated(MOM_dom)) then
    allocate(MOM_dom)
    allocate(MOM_dom%mpp_domain)
  endif

  pe = PE_here()
  proc_used = num_PEs()

  if (present(symmetric)) then
    MOM_dom%symmetric = symmetric
    if (symmetric) then ; mod = "MOM_domains symmetric"
    else ; mod = "MOM_domains non-symmetric" ; endif
  else
    mod = "MOM_domains"
#ifdef SYMMETRIC_MEMORY_
    MOM_dom%symmetric = .true.
#else
    MOM_dom%symmetric = .false.
#endif
  endif
  if (present(min_halo)) mod = trim(mod)//" min_halo"
  
  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, tagname)
  call get_param(param_file, mod, "REENTRANT_X", reentrant_x, &
                 "If true, the domain is zonally reentrant.", default=.true.)
  call get_param(param_file, mod, "REENTRANT_Y", reentrant_y, &
                 "If true, the domain is meridionally reentrant.", &
                 default=.false.)
  call get_param(param_file, mod, "TRIPOLAR_N", tripolar_N, &
                 "Use tripolar connectivity at the northern edge of the \n"//&
                 "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", &
                 default=.false.)
  
  call log_param(param_file, mod, "!SYMMETRIC_MEMORY_", MOM_dom%symmetric, &
                 "If defined, the velocity point data domain includes \n"//&
                 "every face of the thickness points. In other words, \n"//&
                 "some arrays are larger than others, depending on where \n"//&
                 "they are on the staggered grid.  Also, the starting \n"//&
                 "index of the velocity-point arrays is usually 0, not 1. \n"//&
                 "This can only be set at compile time.")
  call get_param(param_file, mod, "NONBLOCKING_UPDATES", MOM_dom%nonblocking_updates, &
                 "If true, non-blocking halo updates may be used.", &
                 default=.false.)

  is_static = .false.
  nihalo_dflt = 2 ; njhalo_dflt = 2
#ifdef STATIC_MEMORY_
  is_static = .true.
  nihalo_dflt = NIHALO_ ; njhalo_dflt = NJHALO_
#else
# ifdef NIHALO_
  nihalo_dflt = NIHALO_
# endif
# ifdef NJHALO_
  njhalo_dflt = NJHALO_
# endif
#endif
  call log_param(param_file, mod, "!STATIC_MEMORY_", is_static, &
                 "If STATIC_MEMORY_ is defined, the principle variables \n"//&
                 "will have sizes that are statically determined at \n"//&
                 "compile time.  Otherwise the sizes are not determined \n"//&
                 "until run time. The STATIC option is substantially \n"//&
                 "faster, but does not allow the PE count to be changed \n"//&
                 "at run time.  This can only be set at compile time.")

  call get_param(param_file, mod, "NIHALO", MOM_dom%nihalo, &
                 "The number of halo points on each side in the \n"//&
                 "x-direction.  With STATIC_MEMORY_ this is set as NIHALO_ \n"//&
                 "in MOM_memory.h at compile time; without STATIC_MEMORY_ \n"//&
                 "the default is NIHALO_ in MOM_memory.h (if defined) or 2.", &
                 default=nihalo_dflt)
  call get_param(param_file, mod, "NJHALO", MOM_dom%njhalo, &
                 "The number of halo points on each side in the \n"//&
                 "y-direction.  With STATIC_MEMORY_ this is set as NJHALO_ \n"//&
                 "in MOM_memory.h at compile time; without STATIC_MEMORY_ \n"//&
                 "the default is NJHALO_ in MOM_memory.h (if defined) or 2.", &
                 default=njhalo_dflt)
  if (present(min_halo)) then
    MOM_dom%nihalo = max(MOM_dom%nihalo, min_halo(1))
    min_halo(1) = MOM_dom%nihalo
    MOM_dom%njhalo = max(MOM_dom%njhalo, min_halo(2))
    min_halo(2) = MOM_dom%njhalo
    call log_param(param_file, mod, "!NIHALO min_halo", MOM_dom%nihalo)
    call log_param(param_file, mod, "!NJHALO min_halo", MOM_dom%nihalo)
  endif
#ifdef STATIC_MEMORY_
  call get_param(param_file, mod, "NIGLOBAL", MOM_dom%niglobal, &
                 "The total number of thickness grid points in the \n"//&
                 "x-direction in the physical domain. With STATIC_MEMORY_ \n"//&
                 "this is set in MOM_memory.h at compile time.", default=NIGLOBAL_)
  call get_param(param_file, mod, "NJGLOBAL", MOM_dom%njglobal, &
                 "The total number of thickness grid points in the \n"//&
                 "x-direction in the physical domain. With STATIC_MEMORY_ \n"//&
                 "this is set in MOM_memory.h at compile time.", default=NJGLOBAL_)
  if (MOM_dom%niglobal /= NIGLOBAL_) call MOM_error(FATAL,"MOM_domains_init: " // &
   "static mismatch for NIGLOBAL_ domain size. Header file does not match input namelist")
  if (MOM_dom%njglobal /= NJGLOBAL_) call MOM_error(FATAL,"MOM_domains_init: " // &
   "static mismatch for NJGLOBAL_ domain size. Header file does not match input namelist")

  if (.not.present(min_halo)) then
    if (MOM_dom%nihalo /= NIHALO_) call MOM_error(FATAL,"MOM_domains_init: " // &
           "static mismatch for NIHALO domain size")
    if (MOM_dom%njhalo /= NJHALO_) call MOM_error(FATAL,"MOM_domains_init: " // &
           "static mismatch for NJHALO domain size")
  endif
#else
  call get_param(param_file, mod, "NIGLOBAL", MOM_dom%niglobal, &
                 "The total number of thickness grid points in the \n"//&
                 "x-direction in the physical domain. With STATIC_MEMORY_ \n"//&
                 "this is set in MOM_memory.h at compile time.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "NJGLOBAL", MOM_dom%njglobal, &
                 "The total number of thickness grid points in the \n"//&
                 "x-direction in the physical domain. With STATIC_MEMORY_ \n"//&
                 "this is set in MOM_memory.h at compile time.", &
                 fail_if_missing=.true.)
#endif
  nihalo = MOM_dom%nihalo
  njhalo = MOM_dom%njhalo

  global_indices(1) = 1
  global_indices(2) = MOM_dom%niglobal
  global_indices(3) = 1
  global_indices(4) = MOM_dom%njglobal

#ifdef STATIC_MEMORY_
  layout(1) = NIPROC_ ; layout(2) = NJPROC_
#else
  call mpp_define_layout(global_indices, proc_used, layout)
  call read_param(param_file,"NIPROC",layout(1))
  call read_param(param_file,"NJPROC",layout(2))
  if (layout(1)*layout(2) /= proc_used) then
    write(mesg,'("MOM_domains_init: The product of the two components of layout, ", &
          &      2i4,", is not the number of PEs used, ",i5,".")') &
          layout(1),layout(2),proc_used
    call MOM_error(FATAL, mesg)
  endif
#endif
  call log_param(param_file, mod, "!NIPROC", layout(1), &
                 "The number of processors in the x-direction. With \n"//&
                 "STATIC_MEMORY_ this is set in MOM_memory.h at compile time.")
  call log_param(param_file, mod, "!NJPROC", layout(2), &
                 "The number of processors in the x-direction. With \n"//&
                 "STATIC_MEMORY_ this is set in MOM_memory.h at compile time.")
!  write(*,*) 'layout is now ',layout, global_indices

  !   Set up the I/O lay-out, and check that it uses an even multiple of the
  ! number of PEs in each direction.
  io_layout(:) = (/ 0, 0 /)
  call get_param(param_file, mod, "NIPROC_IO", io_layout(1), &
                 "The number of processors used for I/O in the \n"//&
                 "x-direction, or 0 to equal NIPROC.", default=0)
  call get_param(param_file, mod, "NJPROC_IO", io_layout(2), &
                 "The number of processors used for I/O in the \n"//&
                 "y-direction, or 0 to equal NJPROC.", default=0)
  if (io_layout(1) < 0) then
    write(mesg,'("MOM_domains_init: NIPROC_IO = ",i4,".  Negative values of "//&
         &" of NIPROC_IO are not allowed.")') io_layout(1)
    call MOM_error(FATAL, mesg)
  elseif (io_layout(1) > 0) then ; if (modulo(layout(1), io_layout(1)) /= 0) then
    write(mesg,'("MOM_domains_init: The x-direction I/O-layout, NIPROC_IO=",i4, &
         &", does not evenly divide the x-direction layout, NIPROC=,",i4,".")') &
          io_layout(1),layout(1)
    call MOM_error(FATAL, mesg)
  endif ; endif
  
  if (io_layout(2) < 0) then
    write(mesg,'("MOM_domains_init: NJPROC_IO = ",i4,".  Negative values of "//&
         &" of NJPROC_IO are not allowed.")') io_layout(2)
    call MOM_error(FATAL, mesg)
  elseif (io_layout(2) /= 0) then ; if (modulo(layout(2), io_layout(2)) /= 0) then
    write(mesg,'("MOM_domains_init: The y-direction I/O-layout, NJPROC_IO=",i4, &
         &", does not evenly divide the y-direction layout, NJPROC=,",i4,".")') &
          io_layout(2),layout(2)
    call MOM_error(FATAL, mesg)
  endif ; endif

  X_FLAGS = 0 ; Y_FLAGS = 0
  if (reentrant_x) X_FLAGS = CYCLIC_GLOBAL_DOMAIN
  if (reentrant_y) Y_FLAGS = CYCLIC_GLOBAL_DOMAIN
  if (tripolar_N) then
    Y_FLAGS = FOLD_NORTH_EDGE
    if (reentrant_y) call MOM_error(FATAL,"MOM_domains: "// &
      "TRIPOLAR_N and REENTRANT_Y may not be defined together.")
  endif
  
  call MOM_define_domain((/1+nihalo,MOM_dom%niglobal+nihalo,1+njhalo, &
                 MOM_dom%njglobal+njhalo/), layout, MOM_dom%mpp_domain, &
                 xflags=X_FLAGS, yflags=Y_FLAGS, xhalo=nihalo, yhalo=njhalo, &
                 symmetry = MOM_dom%symmetric, name="MOM")

  if ((io_layout(1) + io_layout(2) > 0)) then
    call MOM_define_io_domain(MOM_dom%mpp_domain, io_layout)
  endif

! Save the extra data for creating other domains of different resolution that overlay this domain
  MOM_dom%X_FLAGS = X_FLAGS
  MOM_dom%Y_FLAGS = Y_FLAGS
  MOM_dom%layout = layout
  MOM_dom%io_layout = io_layout
  MOM_dom%use_io_layout = (io_layout(1) + io_layout(2) > 0)

#ifdef STATIC_MEMORY_
!   A requirement of equal sized compute domains is necessary when STATIC_MEMORY_
! is used.
  call mpp_get_compute_domain(MOM_dom%mpp_domain,isc,iec,jsc,jec)
  xsiz = iec - isc + 1
  ysiz = jec - jsc + 1
  if (xsiz*NIPROC_ /= MOM_dom%niglobal .OR. ysiz*NJPROC_ /= MOM_dom%njglobal) then
     write( char_xsiz,'(i4)' ) NIPROC_
     write( char_ysiz,'(i4)' ) NJPROC_
     write( char_niglobal,'(i4)' ) MOM_dom%niglobal
     write( char_njglobal,'(i4)' ) MOM_dom%njglobal
     call MOM_error(WARNING,'MOM_domains: Processor decomposition (NIPROC_,NJPROC_) = (' &
         //trim(char_xsiz)//','//trim(char_ysiz)// &
         ') does not evenly divide size set by preprocessor macro ('&
         //trim(char_niglobal)//','//trim(char_njglobal)// '). ')
     call MOM_error(FATAL,'MOM_domains:  #undef STATIC_MEMORY_ in MOM_memory.h to use &
         &dynamic allocation, or change processor decomposition to evenly divide the domain.')
  endif
#endif

end subroutine MOM_domains_init

subroutine get_domain_extent(Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, &
                             isg, ieg, jsg, jeg, idg_offset, jdg_offset, symmetric)
  type(MOM_domain_type), intent(in) :: Domain
  integer, intent(out) :: isc, iec, jsc, jec
  integer, intent(out) :: isd, ied, jsd, jed
  integer, intent(out) :: isg, ieg, jsg, jeg
  integer, intent(out) :: idg_offset, jdg_offset
  logical, intent(out) :: symmetric
! Arguments: Domain - The MOM_domain_type from which the indices are extracted.
!  (out)     isc, iec, jsc, jec - the start & end indices of the
!                                 computational domain.
!  (out)     isd, ied, jsd, jed - the start & end indices of the data domain.
!  (out)     isg, ieg, jsg, jeg - the start & end indices of the global domain.
!  (out)     idg_offset, jdg_offset - the offset between the corresponding
!                                     global and data index spaces.
!  (out)     symmetric - true if symmetric memory is used.

  call mpp_get_compute_domain(Domain%mpp_domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(Domain%mpp_domain, isd, ied, jsd, jed)
  call mpp_get_global_domain(Domain%mpp_domain, isg, ieg, jsg, jeg)

  ! This code institutes the MOM convention that local array indices start at 1.
  idg_offset = isd-1 ; jdg_offset = jsd-1
  isc = isc-isd+1 ; iec = iec-isd+1 ; jsc = jsc-jsd+1 ; jec = jec-jsd+1
  ied = ied-isd+1 ; jed = jed-jsd+1
  isd = 1 ; jsd = 1
  symmetric = Domain%symmetric

end subroutine get_domain_extent

end module MOM_domains
