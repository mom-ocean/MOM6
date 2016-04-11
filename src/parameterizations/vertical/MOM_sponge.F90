module MOM_sponge
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
use MOM_verticalGrid, only : verticalGrid_type

! Planned extension:  Support for time varying sponge targets.

implicit none ; private

#include <MOM_memory.h>

public set_up_sponge_field, set_up_sponge_ML_density
public initialize_sponge, apply_sponge, sponge_end, init_sponge_diags

type :: p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d
type :: p2d
  real, dimension(:,:), pointer :: p => NULL()
end type p2d

type, public :: sponge_CS ; private
  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used with
                             ! nkml sublayers and nkbl buffer layer.
  integer :: nz              ! The total number of layers.
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
  real, pointer :: Rcv_ml_ref(:) => NULL() ! The value toward which the mixed layer
                             ! coordinate-density is being damped, in kg m-3.
  real, pointer :: Ref_eta(:,:) => NULL() ! The value toward which the interface
                             ! heights are being damped.
  type(p3d) :: var(MAX_FIELDS_)  ! Pointers to the fields that are being damped.
  type(p2d) :: Ref_val(MAX_FIELDS_)  ! The values to which the fields are damped.

  logical :: do_i_mean_sponge     ! If true, apply sponges to the i-mean fields.
  real, pointer :: Iresttime_im(:) => NULL()  ! The inverse restoring time of
                             ! each row for i-mean sponges.
  real, pointer :: Rcv_ml_ref_im(:) => NULL() ! The value toward which the i-mean
                             ! mixed layer coordinate-density is being damped,
                             ! in kg m-3.
  real, pointer :: Ref_eta_im(:,:) => NULL() ! The value toward which the i-mean
                             ! interface heights are being damped.
  type(p2d) :: Ref_val_im(MAX_FIELDS_)  ! The values toward which the i-means of
                             ! fields are damped.

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer :: id_w_sponge = -1

end type sponge_CS

contains

subroutine initialize_sponge(Iresttime, int_height, G, param_file, CS, &
                             Iresttime_i_mean, int_height_i_mean)
  type(ocean_grid_type),                  intent(in) :: G
  real, dimension(SZI_(G),SZJ_(G)),       intent(in) :: Iresttime
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(in) :: int_height
  type(param_file_type),                  intent(in) :: param_file
  type(sponge_CS),                        pointer    :: CS
  real, dimension(SZJ_(G)),     optional, intent(in) :: Iresttime_i_mean
  real, dimension(SZJ_(G),SZK_(G)+1), optional, intent(in) :: int_height_i_mean

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

  if (present(Iresttime_i_mean) .neqv. present(int_height_i_mean)) &
    call MOM_error(FATAL, "initialize_sponge:  The optional arguments \n"//&
           "Iresttime_i_mean and int_height_i_mean must both be present \n"//&
           "if either one is.")

  CS%do_i_mean_sponge = present(Iresttime_i_mean)

  CS%nz = G%ke
  CS%isc = G%isc ; CS%iec = G%iec ; CS%jsc = G%jsc ; CS%jec = G%jec
  CS%isd = G%isd ; CS%ied = G%ied ; CS%jsd = G%jsd ; CS%jed = G%jed
  ! CS%bulkmixedlayer may be set later via a call to set_up_sponge_ML_density.
  CS%bulkmixedlayer = .false.

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

    allocate(CS%Ref_eta(CS%nz+1,CS%num_col))
    do col=1,CS%num_col ; do K=1,CS%nz+1
      CS%Ref_eta(K,col) = int_height(CS%col_i(col),CS%col_j(col),K)
    enddo ; enddo

  endif

  if (CS%do_i_mean_sponge) then
    allocate(CS%Iresttime_im(G%jsd:G%jed)) ; CS%Iresttime_im(:) = 0.0
    allocate(CS%Ref_eta_im(G%jsd:G%jed,G%ke+1)) ; CS%Ref_eta_im(:,:) = 0.0

    do j=G%jsc,G%jec
      CS%Iresttime_im(j) = Iresttime_i_mean(j)
    enddo
    do K=1,CS%nz+1 ; do j=G%jsc,G%jec
      CS%Ref_eta_im(j,K) = int_height_i_mean(j,K)
    enddo ; enddo
  endif

  total_sponge_cols = CS%num_col
  call sum_across_PEs(total_sponge_cols)

  call log_param(param_file, mod, "!Total sponge columns", total_sponge_cols, &
                 "The total number of columns where sponges are applied.")

end subroutine initialize_sponge

subroutine init_sponge_diags(Time, G, diag, CS)
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

end subroutine init_sponge_diags

subroutine set_up_sponge_field(sp_val, f_ptr, G, nlay, CS, sp_val_i_mean)
  type(ocean_grid_type),                            intent(in) :: G
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),         intent(in) :: sp_val
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target, intent(in) :: f_ptr
  integer,                                          intent(in) :: nlay
  type(sponge_CS),                                  pointer    :: CS
  real, dimension(SZJ_(G),SZK_(G)),       optional, intent(in) :: sp_val_i_mean
!   This subroutine stores the reference profile for the variable
! whose address is given by f_ptr. nlay is the number of layers in
! this variable.

! Arguments: sp_val - The reference profiles of the quantity being
!                     registered.
!  (in)      f_ptr - a pointer to the field which will be damped.
!  (in)      nlay - the number of layers in this quantity.
!  (in/out)  CS - A pointer to the control structure for this module that is
!                 set by a previous call to initialize_sponge.
!  (in,opt)  sp_val_i_mean - The i-mean reference value for this field with
!                 i-mean sponges.

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

  allocate(CS%Ref_val(CS%fldno)%p(CS%nz,CS%num_col))
  CS%Ref_val(CS%fldno)%p(:,:) = 0.0
  do col=1,CS%num_col
    do k=1,nlay
      CS%Ref_val(CS%fldno)%p(k,col) = sp_val(CS%col_i(col),CS%col_j(col),k)
    enddo
    do k=nlay+1,CS%nz
      CS%Ref_val(CS%fldno)%p(k,col) = 0.0
    enddo
  enddo

  CS%var(CS%fldno)%p => f_ptr

  if (nlay/=CS%nz) then
    write(mesg,'("Danger: Sponge reference fields require nz (",I3,") layers.&
        & A field with ",I3," layers was passed to set_up_sponge_field.")') &
          CS%nz, nlay
    if (is_root_pe()) call MOM_error(WARNING, "set_up_sponge_field: "//mesg)
  endif

  if (CS%do_i_mean_sponge) then
    if (.not.present(sp_val_i_mean)) call MOM_error(FATAL, &
      "set_up_sponge_field: sp_val_i_mean must be present with i-mean sponges.")

    allocate(CS%Ref_val_im(CS%fldno)%p(CS%jsd:CS%jed,CS%nz))
    CS%Ref_val(CS%fldno)%p(:,:) = 0.0
    do k=1,CS%nz ; do j=CS%jsc,CS%jec
      CS%Ref_val_im(CS%fldno)%p(j,k) = sp_val_i_mean(j,k)
    enddo ; enddo
  endif

end subroutine set_up_sponge_field


subroutine set_up_sponge_ML_density(sp_val, G, CS, sp_val_i_mean)
  type(ocean_grid_type),              intent(in) :: G
  real, dimension(SZI_(G),SZJ_(G)),   intent(in) :: sp_val
  type(sponge_CS),                    pointer    :: CS
  real, dimension(SZJ_(G)), optional, intent(in) :: sp_val_i_mean
!   This subroutine stores the reference value for mixed layer density.  It is
! handled differently from other values because it is only used in determining
! which layers can be inflated.

! Arguments: sp_val - The reference values of the mixed layer density.
!  (in/out)  CS - A pointer to the control structure for this module that is
!                 set by a previous call to initialize_sponge.

  integer :: j, col
  character(len=256) :: mesg ! String for error messages

  if (.not.associated(CS)) return

  if (associated(CS%Rcv_ml_ref)) then
    call MOM_error(FATAL, "set_up_sponge_ML_density appears to have been "//&
                           "called twice.")
  endif

  CS%bulkmixedlayer = .true.
  allocate(CS%Rcv_ml_ref(CS%num_col)) ; CS%Rcv_ml_ref(:) = 0.0
  do col=1,CS%num_col
    CS%Rcv_ml_ref(col) = sp_val(CS%col_i(col),CS%col_j(col))
  enddo

  if (CS%do_i_mean_sponge) then
    if (.not.present(sp_val_i_mean)) call MOM_error(FATAL, &
      "set_up_sponge_field: sp_val_i_mean must be present with i-mean sponges.")

    allocate(CS%Rcv_ml_ref_im(CS%jsd:CS%jed)) ; CS%Rcv_ml_ref_im(:) = 0.0
    do j=CS%jsc,CS%jec
      CS%Rcv_ml_ref_im(j) = sp_val_i_mean(j)
    enddo
  endif

end subroutine set_up_sponge_ML_density

subroutine apply_sponge(h, dt, G, GV, ea, eb, CS, Rcv_ml)
  type(ocean_grid_type),                    intent(inout) :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h
  real,                                     intent(in)    :: dt
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out)   :: ea
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out)   :: eb
  type(sponge_CS),                          pointer       :: CS
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(inout) :: Rcv_ml

! This subroutine applies damping to the layers thicknesses, mixed
! layer buoyancy, and a variety of tracers for every column where
! there is damping.

! Arguments: h -  Layer thickness, in m.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (out)     ea - an array to which the amount of fluid entrained
!                 from the layer above during this call will be
!                 added, in m.
!  (out)     eb - an array to which the amount of fluid entrained
!                 from the layer below during this call will be
!                 added, in m.
!  (in)      CS - A pointer to the control structure for this module that is
!                 set by a previous call to initialize_sponge.
!  (inout,opt)  Rcv_ml - The coordinate density of the mixed layer.

  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) :: &
    w_int, &       ! Water moved upward across an interface within a timestep,
                   ! in m.
    e_D            ! Interface heights that are dilated to have a value of 0
                   ! at the surface, in m.
  real, dimension(SZI_(G), SZJ_(G)) :: &
    eta_anom, &    ! Anomalies in the interface height, relative to the i-mean
                   ! target value.
    fld_anom       ! Anomalies in a tracer concentration, relative to the
                   ! i-mean target value.
  real, dimension(SZJ_(G), SZK_(G)+1) :: &
    eta_mean_anom  ! The i-mean interface height anomalies.
  real, allocatable, dimension(:,:,:) :: &
    fld_mean_anom  ! THe i-mean tracer concentration anomalies.
  real, dimension(SZI_(G), SZK_(G)+1) :: &
    h_above, &     ! The total thickness above an interface, in m.
    h_below        ! The total thickness below an interface, in m.
  real, dimension(SZI_(G)) :: &
    dilate         ! A nondimensional factor by which to dilate layers to
                   ! give 0 at the surface.

  real :: e(SZK_(G)+1)  ! The interface heights, in m or kg m-2, usually negative.
  real :: e0       ! The height of the free surface in m.
  real :: e_str    ! A nondimensional amount by which the reference
                   ! profile must be stretched for the free surfaces
                   ! heights in the two profiles to agree.
  real :: w        ! The thickness of water moving upward through an
                   ! interface within 1 timestep, in m.
  real :: wm       ! wm is w if w is negative and 0 otherwise, in m.
  real :: wb       ! w at the interface below a layer, in m.
  real :: wpb      ! wpb is wb if wb is positive and 0 otherwise, m.
  real :: ea_k, eb_k
  real :: damp     ! The timestep times the local damping  coefficient.  ND.
  real :: I1pdamp  ! I1pdamp is 1/(1 + damp).  Nondimensional.
  real :: damp_1pdamp ! damp_1pdamp is damp/(1 + damp).  Nondimensional.
  real :: Idt      ! 1.0/dt, in s-1.
  integer :: c, m, nkmb, i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) return
  if (CS%bulkmixedlayer) nkmb = GV%nk_rho_varies
  if (CS%bulkmixedlayer .and. (.not.present(Rcv_ml))) &
    call MOM_error(FATAL, "Rml must be provided to apply_sponge when using "//&
                           "a bulk mixed layer.")

  if ((CS%id_w_sponge > 0) .or. CS%do_i_mean_sponge) then
    do k=1,nz+1 ; do j=js,je ; do i=is,ie
      w_int(i,j,K) = 0.0
    enddo ; enddo ; enddo
  endif

  if (CS%do_i_mean_sponge) then
    ! Apply forcing to restore the zonal-mean properties to prescribed values.

    if (CS%bulkmixedlayer) call MOM_error(FATAL, "apply_sponge is not yet set up to "//&
                  "work properly with i-mean sponges and a bulk mixed layer.")

    do j=js,je ; do i=is,ie ; e_D(i,j,nz+1) = -G%bathyT(i,j) ; enddo ; enddo
    do k=nz,1,-1 ; do j=js,je ; do i=is,ie
      e_D(i,j,K) = e_D(i,j,K+1) + h(i,j,k)
    enddo ; enddo ; enddo
    do j=js,je
      do i=is,ie
        dilate(i) = G%bathyT(i,j) / (e_D(i,j,1) + G%bathyT(i,j))
      enddo
      do k=1,nz+1 ; do i=is,ie
        e_D(i,j,K) = dilate(i) * (e_D(i,j,K) + G%bathyT(i,j)) - G%bathyT(i,j)
      enddo ; enddo
    enddo

    do k=2,nz
      do j=js,je ; do i=is,ie
        eta_anom(i,j) = e_D(i,j,k) - CS%Ref_eta_im(j,k)
        if (CS%Ref_eta_im(j,K) < -G%bathyT(i,j)) eta_anom(i,j) = 0.0
      enddo ; enddo
      call global_i_mean(eta_anom(:,:), eta_mean_anom(:,K), G)
    enddo

    if (CS%fldno > 0) allocate(fld_mean_anom(G%isd:G%ied,nz,CS%fldno))
    do m=1,CS%fldno
      do j=js,je ; do i=is,ie
        fld_anom(i,j) = CS%var(m)%p(i,j,k) - CS%Ref_val_im(m)%p(j,k)
      enddo ; enddo
      call global_i_mean(fld_anom(:,:), fld_mean_anom(:,k,m), G, h(:,:,k))
    enddo

    do j=js,je ; if (CS%Iresttime_im(j) > 0.0) then
      damp = dt*CS%Iresttime_im(j) ; damp_1pdamp = damp / (1.0 + damp)

      do i=is,ie
        h_above(i,1) = 0.0 ; h_below(i,nz+1) = 0.0
      enddo
      do K=nz,1,-1 ; do i=is,ie
        h_below(i,K) = h_below(i,K+1) + max(h(i,j,k)-GV%Angstrom, 0.0)
      enddo ; enddo
      do K=2,nz+1 ; do i=is,ie
        h_above(i,K) = h_above(i,K-1) + max(h(i,j,k-1)-GV%Angstrom, 0.0)
      enddo ; enddo
      do K=2,nz
        ! w is positive for an upward (lightward) flux of mass, resulting
        ! in the downward movement of an interface.
        w = damp_1pdamp * eta_mean_anom(j,K)
        do i=is,ie
          if (w > 0.0) then
            w_int(i,j,K) = min(w, h_below(i,K))
            eb(i,j,k-1) = eb(i,j,k-1) + w_int(i,j,K)
          else
            w_int(i,j,K) = max(w, -h_above(i,K))
            ea(i,j,k) = ea(i,j,k) - w_int(i,j,K)
          endif
        enddo
      enddo
      do k=1,nz ; do i=is,ie
        ea_k = max(0.0, -w_int(i,j,K))
        eb_k = max(0.0, w_int(i,j,K+1))
        do m=1,CS%fldno
          CS%var(m)%p(i,j,k) = (h(i,j,k)*CS%var(m)%p(i,j,k) + &
              CS%Ref_val_im(m)%p(j,k) * (ea_k + eb_k)) / &
                     (h(i,j,k) + (ea_k + eb_k)) - &
              damp_1pdamp * fld_mean_anom(j,k,m)
        enddo

        h(i,j,k) = max(h(i,j,k) + (w_int(i,j,K+1) - w_int(i,j,K)), &
                       min(h(i,j,k), GV%Angstrom))
      enddo ; enddo
    endif ; enddo

    if (CS%fldno > 0) deallocate(fld_mean_anom)

  endif

  do c=1,CS%num_col
! c is an index for the next 3 lines but a multiplier for the rest of the loop
! Therefore we use c as per C code and increment the index where necessary.
    i = CS%col_i(c) ; j = CS%col_j(c)
    damp = dt*CS%Iresttime_col(c)

    e(1) = 0.0 ; e0 = 0.0
    do K=1,nz
      e(K+1) = e(K) - h(i,j,k)
    enddo
    e_str = e(nz+1) / CS%Ref_eta(nz+1,c)

    if ( CS%bulkmixedlayer ) then
      I1pdamp = 1.0 / (1.0 + damp)
      if (associated(CS%Rcv_ml_ref)) &
        Rcv_ml(i,j) = I1pdamp * (Rcv_ml(i,j) + CS%Rcv_ml_ref(c)*damp)
      do k=1,nkmb
        do m=1,CS%fldno
          CS%var(m)%p(i,j,k) = I1pdamp * &
              (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(k,c)*damp)
        enddo
      enddo

      wpb = 0.0; wb = 0.0
      do k=nz,nkmb+1,-1
        if (GV%Rlay(k) > Rcv_ml(i,j)) then
          w = MIN((((e(K)-e0) - e_str*CS%Ref_eta(K,c)) * damp), &
                    ((wb + h(i,j,k)) - GV%Angstrom))
          wm = 0.5*(w-ABS(w))
          do m=1,CS%fldno
            CS%var(m)%p(i,j,k) = (h(i,j,k)*CS%var(m)%p(i,j,k) + &
                     CS%Ref_val(m)%p(k,c)*(damp*h(i,j,k) + (wpb - wm))) / &
                     (h(i,j,k)*(1.0 + damp) + (wpb - wm))
          enddo
        else
          do m=1,CS%fldno
            CS%var(m)%p(i,j,k) = I1pdamp * &
              (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(k,c)*damp)
          enddo
          w = wb + (h(i,j,k) - GV%Angstrom)
          wm = 0.5*(w-ABS(w))
        endif
        eb(i,j,k) = eb(i,j,k) + wpb
        ea(i,j,k) = ea(i,j,k) - wm
        h(i,j,k)  = h(i,j,k)  + (wb - w)
        wb = w
        wpb = w - wm
      enddo

      if (wb < 0) then
        do k=nkmb,1,-1
          w = MIN((wb + (h(i,j,k) - GV%Angstrom)),0.0)
          h(i,j,k)  = h(i,j,k)  + (wb - w)
          ea(i,j,k) = ea(i,j,k) - w
          wb = w
        enddo
      else
        w = wb
        do k=GV%nkml,nkmb
          eb(i,j,k) = eb(i,j,k) + w
        enddo

        k = GV%nkml
        h(i,j,k) = h(i,j,k) + w
        do m=1,CS%fldno
          CS%var(m)%p(i,j,k) = (CS%var(m)%p(i,j,k)*h(i,j,k) + &
                                CS%Ref_val(m)%p(k,c)*w) / (h(i,j,k) + w)
        enddo
      endif

      do k=1,nkmb
        do m=1,CS%fldno
          CS%var(m)%p(i,j,k) = I1pdamp * &
              (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(GV%nkml,c)*damp)
        enddo
      enddo

    else                                          ! not BULKMIXEDLAYER

      wpb = 0.0
      wb = 0.0
      do k=nz,1,-1
        w = MIN((((e(K)-e0) - e_str*CS%Ref_eta(K,c)) * damp), &
                  ((wb + h(i,j,k)) - GV%Angstrom))
        wm = 0.5*(w - ABS(w))
        do m=1,CS%fldno
          CS%var(m)%p(i,j,k) = (h(i,j,k)*CS%var(m)%p(i,j,k) + &
              CS%Ref_val(m)%p(k,c) * (damp*h(i,j,k) + (wpb - wm))) / &
                     (h(i,j,k)*(1.0 + damp) + (wpb - wm))
        enddo
        eb(i,j,k) = eb(i,j,k) + wpb
        ea(i,j,k) = ea(i,j,k) - wm
        h(i,j,k)  = h(i,j,k)  + (wb - w)
        wb = w
        wpb = w - wm
      enddo

    endif                                         ! end BULKMIXEDLAYER
  enddo ! end of c loop

  if (associated(CS%diag)) then ; if (query_averaging_enabled(CS%diag)) then
    Idt = 1.0 / dt
    if (CS%id_w_sponge > 0) call post_data(CS%id_w_sponge, Idt*w_int, CS%diag)
  endif ; endif

end subroutine apply_sponge

subroutine sponge_end(CS)
  type(sponge_CS),              pointer       :: CS
!  (in)      CS - A pointer to the control structure for this module that is
!                 set by a previous call to initialize_sponge.
  integer :: m

  if (.not.associated(CS)) return

  if (associated(CS%col_i)) deallocate(CS%col_i)
  if (associated(CS%col_j)) deallocate(CS%col_j)

  if (associated(CS%Iresttime_col)) deallocate(CS%Iresttime_col)
  if (associated(CS%Rcv_ml_ref)) deallocate(CS%Rcv_ml_ref)
  if (associated(CS%Ref_eta)) deallocate(CS%Ref_eta)

  if (associated(CS%Iresttime_im)) deallocate(CS%Iresttime_im)
  if (associated(CS%Rcv_ml_ref_im)) deallocate(CS%Rcv_ml_ref_im)
  if (associated(CS%Ref_eta_im)) deallocate(CS%Ref_eta_im)

  do m=1,CS%fldno
    if (associated(CS%Ref_val(CS%fldno)%p)) deallocate(CS%Ref_val(CS%fldno)%p)
    if (associated(CS%Ref_val_im(CS%fldno)%p)) &
      deallocate(CS%Ref_val_im(CS%fldno)%p)
  enddo

  deallocate(CS)

end subroutine sponge_end

end module MOM_sponge
