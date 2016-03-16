module MOM_ALE_sponge
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

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, March 1999-June 2000                           *
!*  Modified by Gustavo Marques, March 2016
!*                                                                     *
!*    This program contains the subroutines that implement sponge      *
!*  regions, in which the stratification and water mass properties     *
!*  are damped toward some profiles.  There are three externally       *
!*  callable subroutines in this file.                                 *
!*                                                                     *
!*    initialize_sponge determines the mapping from the model          *
!*  variables into the arrays of damped columns.  This remapping is    *
!*  done for efficiency and to conserve memory.  Only columns which    *
!*  have positive inverse damping times and which are deeper than a    *
!*  supplied depth are placed in sponges.  The inverse damping         *
!*  time is also stored in this subroutine, and memory is allocated    *
!*  for all of the reference profiles which will subsequently be       *
!*  provided through calls to set_up_sponge_field.  The first two      *
!*  arguments are a two-dimensional array containing the damping       *
!*  rates, and the interface heights to damp towards.                  *
!*                                                                     *
!*    set_up_sponge_field is called to provide a reference profile     *
!*  and the location of the field that will be damped back toward      *
!*  that reference profile.  A third argument, the number of layers    *
!*  in the field is also provided, but this should always be nz.       *
!*                                                                     *
!*    Apply_sponge damps all of the fields that have been registered   *
!*  with set_up_sponge_field toward their reference profiles.  The     *
!*  four arguments are the thickness to be damped, the amount of time  *
!*  over which the damping occurs, and arrays to which the movement    *
!*  of fluid into a layer from above and below will be added. The      *
!*  effect on momentum of the sponge may be accounted for later using  *
!*  the movement of water recorded in these later arrays.              *
!*                                                                     *
!*    All of the variables operated upon in this file are defined at   *
!*  the thickness points.                                              *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, T, S, Iresttime, ea, eb               *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_coms, only : sum_across_PEs
use MOM_diag_mediator, only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, NOTE, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_spatial_means, only : global_i_mean
use MOM_time_manager, only : time_type
use MOM_remapping, only : remapping_cs, remapping_core_h, initialize_remapping

! GM, Do we need remapping_CS? I don't hink so...

! March 11 2016, NOTES:
! Need to update MOM6, the code in the wiki page is newer
! Need to add comments to the functions below
! Look at MOM wiki --> doxygen https://github.com/NOAA-GFDL/MOM6/wiki/Doxygen

! Planned extension:  Support for time varying sponge targets.

implicit none ; private

#include <MOM_memory.h>

public set_up_ALE_sponge_field
public initialize_ALE_sponge, apply_ALE_sponge, ALE_sponge_end, init_ALE_sponge_diags

type :: p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d
type :: p2d
  real, dimension(:,:), pointer :: p => NULL()
end type p2d

type, public :: sponge_CS ; private
  integer :: nz              ! The total number of layers.
  integer :: nz_data              ! The total number of arbritary layers.
  integer :: isc, iec, jsc, jec  ! The index ranges of the computational domain.
  integer :: isd, ied, jsd, jed  ! The index ranges of the data domain.
  integer :: num_col         ! The number of sponge points within the
                             ! computational domain.
  integer :: fldno = 0       ! The number of fields which have already been
                             ! registered by calls to set_up_sponge_field
  integer, pointer :: col_i(:) => NULL()  ! Arrays containing the i- and j- indicies
  integer, pointer :: col_j(:) => NULL()  ! of each of the columns being damped.
  real, pointer :: Iresttime_col(:) => NULL()  ! The inverse restoring time of
                             ! each column.
  type(p3d) :: var(MAX_FIELDS_)  ! Pointers to the fields that are being damped.
  type(p2d) :: Ref_val(MAX_FIELDS_)  ! The values to which the fields are damped.
  real, pointer :: Ref_h(:,:) => NULL() ! Grid on which reference data is provided

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  type(remapping_cs) :: remap_cs ! remapping control structure

  integer :: id_w_sponge = -1

end type sponge_CS

contains

subroutine initialize_ALE_sponge(Iresttime, data_h, nz_data , G, param_file, CS)

  integer,  intent(in)  :: nz_data
  real, dimension(NIMEM_,NJMEM_),       intent(in) :: Iresttime
  type(ocean_grid_type),                intent(in) :: G
  real, dimension(SZI_(G),SZJ_(G),nz_data), intent(in) :: data_h
  type(param_file_type),                intent(in) :: param_file
  type(sponge_CS),                      pointer    :: CS

!   This subroutine determines the number of points which are within
! sponges in this computational domain.  Only points that have
! positive values of Iresttime and which mask2dT indicates are ocean
! points are included in the sponges.  It also stores the target interface
! heights.

! Arguments: Iresttime - The inverse of the restoring time, in s-1.
!  (in)      int_height - The interface heights to damp back toward, in m.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_sponge"  ! This module's name.
  logical :: use_sponge
  integer :: i, j, k, col, total_sponge_cols

  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_sponge called with an associated "// &
                            "control structure.")
    return
  endif

! Set default, read and log parameters
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "SPONGE", use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified from MOM_initialization.F90.", default=.false.)

  if (.not.use_sponge) return
  
  allocate(CS)

  CS%nz = G%ke
  CS%isc = G%isc ; CS%iec = G%iec ; CS%jsc = G%jsc ; CS%jec = G%jec
  CS%isd = G%isd ; CS%ied = G%ied ; CS%jsd = G%jsd ; CS%jed = G%jed

  CS%num_col = 0 ; CS%fldno = 0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((Iresttime(i,j)>0.0) .and. (G%mask2dT(i,j)>0)) &
      CS%num_col = CS%num_col + 1
  enddo ; enddo

  if (CS%num_col > 0) then

    allocate(CS%Iresttime_col(CS%num_col)) ; CS%Iresttime_col = 0.0
    allocate(CS%col_i(CS%num_col))         ; CS%col_i = 0
    allocate(CS%col_j(CS%num_col))         ; CS%col_j = 0

    col = 1
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if ((Iresttime(i,j)>0.0) .and. (G%mask2dT(i,j)>0)) then
        CS%col_i(col) = i ; CS%col_j(col) = j
        CS%Iresttime_col(col) = Iresttime(i,j)
        col = col +1
      endif
    enddo ; enddo

    CS%nz_data = nz_data
    allocate(CS%Ref_h(CS%nz_data,CS%num_col))
    do col=1,CS%num_col ; do K=1,CS%nz_data
      CS%Ref_h(K,col) = data_h(CS%col_i(col),CS%col_j(col),K)
    enddo ; enddo

  endif

  total_sponge_cols = CS%num_col
  call sum_across_PEs(total_sponge_cols)

!  call initialize_remapping(CS%remap_cs, 'PPM_H4', boundary_extrapolation=.false.)
   call initialize_remapping(nz_data, 'PPM_H4', CS%remap_cs)

  call log_param(param_file, mod, "!Total sponge columns", total_sponge_cols, &
                 "The total number of columns where sponges are applied.")

end subroutine initialize_ALE_sponge

subroutine init_ALE_sponge_diags(Time, G, diag, CS)
  type(time_type),       target, intent(in)    :: Time
  type(ocean_grid_type),         intent(in)    :: G
  type(diag_ctrl),       target, intent(inout) :: diag
  type(sponge_CS),               pointer       :: CS

!   This subroutine sets up diagnostics for the sponges.  It is separate
! from initialize_sponge because it requires fields that are not readily
! availble where initialize_sponge is called.

! Arguments:  Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer to the control structure for this module that is
!                 set by a previous call to initialize_sponge.

  if (.not.associated(CS)) return

  CS%diag => diag
  CS%id_w_sponge = register_diag_field('ocean_model', 'w_sponge', diag%axesTi, &
      Time, 'The diapycnal motion due to the sponges', 'meter second-1')

end subroutine init_ALE_sponge_diags

subroutine set_up_ALE_sponge_field(sp_val, G, f_ptr, CS)
  type(ocean_grid_type),                intent(in) :: G
  type(sponge_CS),                      pointer    :: CS
  real, dimension(SZI_(G),SZJ_(G),CS%nz_data), intent(in) :: sp_val
  real, dimension(NIMEM_,NJMEM_,NKMEM_), target, intent(in) :: f_ptr

!   This subroutine stores the reference profile for the variable
! whose address is given by f_ptr. nlay is the number of layers in
! this variable.

! Arguments: sp_val - The reference profiles of the quantity being
!                     registered.
!  (in)      f_ptr - a pointer to the field which will be damped.
!  (in)      nlay - the number of layers in this quantity.
!  (in/out)  CS - A pointer to the control structure for this module that is
!                 set by a previous call to initialize_sponge.

  integer :: j, k, col
  character(len=256) :: mesg ! String for error messages

  if (.not.associated(CS)) return

  CS%fldno = CS%fldno + 1

  if (CS%fldno > MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ to at least ",I3," in MOM_memory.h or decrease &
           &the number of fields to be damped in the call to &
           &initialize_sponge." )') CS%fldno
    call MOM_error(FATAL,"set_up_sponge_field: "//mesg)
  endif

  allocate(CS%Ref_val(CS%fldno)%p(CS%nz_data,CS%num_col))
  CS%Ref_val(CS%fldno)%p(:,:) = 0.0
  do col=1,CS%num_col
    do k=1,CS%nz_data
      CS%Ref_val(CS%fldno)%p(k,col) = sp_val(CS%col_i(col),CS%col_j(col),k)
    enddo
  enddo

  CS%var(CS%fldno)%p => f_ptr

end subroutine set_up_ALE_sponge_field

subroutine apply_ALE_sponge(h, dt, G, CS)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: h
  real,                                  intent(in)    :: dt
  type(ocean_grid_type),                 intent(inout) :: G
  type(sponge_CS),                       pointer       :: CS

! This subroutine applies damping to the layers thicknesses, mixed
! layer buoyancy, and a variety of tracers for every column where
! there is damping.

! Arguments: h -  Layer thickness, in m.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!                 from the layer below during this call will be
!                 added, in m.
!  (in)      CS - A pointer to the control structure for this module that is
!                 set by a previous call to initialize_sponge.

  real :: damp     ! The timestep times the local damping  coefficient.  ND.
  real :: I1pdamp  ! I1pdamp is 1/(1 + damp).  Nondimensional.
  real :: damp_1pdamp ! damp_1pdamp is damp/(1 + damp).  Nondimensional.
  real :: Idt      ! 1.0/dt, in s-1.
  real :: tmp_val1(cs%nz)      ! data values remapped to model grid
  real :: tmp_val2(cs%nz_data)      ! data values remapped to model grid

  integer :: c, m, nkmb, i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) return

  do c=1,CS%num_col
! c is an index for the next 3 lines but a multiplier for the rest of the loop
! Therefore we use c as per C code and increment the index where necessary.
    i = CS%col_i(c) ; j = CS%col_j(c)
    damp = dt*CS%Iresttime_col(c)

        do m=1,CS%fldno
            tmp_val2(1:CS%nz_data) = CS%Ref_val(m)%p(1:CS%nz_data,c)
            call remapping_core_h(CS%nz_data, CS%Ref_h(:,c), tmp_val2, &
                  CS%nz, h(i,j,:), tmp_val1, CS%remap_cs)

            CS%var(m)%p(i,j,:) = I1pdamp * &
                                 (CS%var(m)%p(i,j,:) + tmp_val1 * damp)

!            CS%var(m)%p(i,j,k) = I1pdamp * &
!                                 (CS%var(m)%p(i,j,k) + tmp_val1 * damp)
        enddo

  enddo ! end of c loop

end subroutine apply_ALE_sponge

subroutine ALE_sponge_end(CS)
  type(sponge_CS),              pointer       :: CS
!  (in)      CS - A pointer to the control structure for this module that is
!                 set by a previous call to initialize_sponge.
  integer :: m

  if (.not.associated(CS)) return

  if (associated(CS%col_i)) deallocate(CS%col_i)
  if (associated(CS%col_j)) deallocate(CS%col_j)

  if (associated(CS%Iresttime_col)) deallocate(CS%Iresttime_col)

  do m=1,CS%fldno
    if (associated(CS%Ref_val(CS%fldno)%p)) deallocate(CS%Ref_val(CS%fldno)%p)
  enddo

  deallocate(CS)

end subroutine ALE_sponge_end
end module MOM_ALE_sponge
