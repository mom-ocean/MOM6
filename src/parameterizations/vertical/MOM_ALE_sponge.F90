!> This module contains the routines used to apply sponge layers when using
!! the ALE mode.
!! Applying sponges requires the following:
!! (1) initialize_ALE_sponge
!! (2) set_up_ALE_sponge_field
!! (3) apply_ALE_sponge
!! (4) init_ALE_sponge_diags (not being used for now)
!! (5) ALE_sponge_end (not being used for now)
module MOM_ALE_sponge

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms, only : sum_across_PEs
use MOM_diag_mediator, only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, NOTE, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_spatial_means, only : global_i_mean
use MOM_time_manager, only : time_type
use MOM_remapping, only : remapping_cs, remapping_core_h, initialize_remapping

! GMM - Planned extension:  Support for time varying sponge targets.

implicit none ; private

#include <MOM_memory.h>

!< Publicly available functions
public set_up_ALE_sponge_field
public initialize_ALE_sponge, apply_ALE_sponge, ALE_sponge_end, init_ALE_sponge_diags

type :: p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d
type :: p2d
  real, dimension(:,:), pointer :: p => NULL()
end type p2d

!> SPONGE control structure
type, public :: ALE_sponge_CS ; private
  integer :: nz              !< The total number of layers.
  integer :: nz_data              !< The total number of arbritary layers.
  integer :: isc, iec, jsc, jec  !< The index ranges of the computational domain.
  integer :: isd, ied, jsd, jed  !< The index ranges of the data domain.
  integer :: num_col         !< The number of sponge points within the
                             ! computational domain.
  integer :: fldno = 0       !< The number of fields which have already been
                             ! registered by calls to set_up_sponge_field
  integer, pointer :: col_i(:) => NULL()  !< Arrays containing the i- and j- indicies
  integer, pointer :: col_j(:) => NULL()  ! of each of the columns being damped.
  real, pointer :: Iresttime_col(:) => NULL()  !< The inverse restoring time of
                                                ! each column.
  type(p3d) :: var(MAX_FIELDS_)      !< Pointers to the fields that are being damped.
  type(p2d) :: Ref_val(MAX_FIELDS_)  !< The values to which the fields are damped.
  real, pointer :: Ref_h(:,:) => NULL() !< Grid on which reference data is provided

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   ! timing of diagnostic output.

  type(remapping_cs) :: remap_cs   !< Remapping parameters and work arrays

end type ALE_sponge_CS

contains

!> This subroutine determines the number of points which are within
! sponges in this computational domain.  Only points that have
! positive values of Iresttime and which mask2dT indicates are ocean
! points are included in the sponges.  It also stores the target interface
! heights.

subroutine initialize_ALE_sponge(Iresttime, data_h, nz_data, G, param_file, CS)

  integer,                              intent(in) :: nz_data !< The total number of arbritary layers (in).
  type(ocean_grid_type),                intent(in) :: G !< The ocean's grid structure (in).
  real, dimension(SZI_(G),SZJ_(G)),     intent(in) :: Iresttime !< The inverse of the restoring time, in s-1 (in).
  real, dimension(SZI_(G),SZJ_(G),nz_data), intent(in) :: data_h !< The thickness to be used in the sponge. It has arbritary layers (in).
  type(param_file_type),                intent(in) :: param_file !< A structure indicating the open file to parse for model parameter values (in). GM: I don't know what this does!
  type(ALE_sponge_CS),                  pointer    :: CS !< A pointer that is set to point to the control structure for this module (in/out).


! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_sponge"  ! This module's name.
  logical :: use_sponge
  logical :: bndExtrapolation = .true. ! If true, extrapolate boundaries
  integer :: i, j, k, col, total_sponge_cols
  character(len=10)  :: remapScheme
  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_sponge called with an associated "// &
                            "control structure.")
    return
  endif

! Set default, read and log parameters
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "SPONGE", use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified from MOM_initialization.F90.", default=.false.)

  if (.not.use_sponge) return

  allocate(CS)

  call get_param(param_file, mod, "REMAPPING_SCHEME", remapScheme, &
                 "This sets the reconstruction scheme used \n"//&
                 " for vertical remapping for all variables.", &
                 default="PLM", do_not_log=.true.)

  call get_param(param_file, mod, "BOUNDARY_EXTRAPOLATION", bndExtrapolation, &
                 "When defined, a proper high-order reconstruction \n"//&
                 "scheme is used within boundary cells rather \n"// &
                 "than PCM. E.g., if PPM is used for remapping, a \n" //&
                 "PPM reconstruction will also be used within boundary cells.", &
                 default=.false., do_not_log=.true.)

  CS%nz = G%ke
  CS%isc = G%isc ; CS%iec = G%iec ; CS%jsc = G%jsc ; CS%jec = G%jec
  CS%isd = G%isd ; CS%ied = G%ied ; CS%jsd = G%jsd ; CS%jed = G%jed

  ! number of columns to be restored
  CS%num_col = 0 ; CS%fldno = 0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((Iresttime(i,j)>0.0) .and. (G%mask2dT(i,j)>0)) &
      CS%num_col = CS%num_col + 1
  enddo ; enddo


  if (CS%num_col > 0) then

    allocate(CS%Iresttime_col(CS%num_col)) ; CS%Iresttime_col = 0.0
    allocate(CS%col_i(CS%num_col))         ; CS%col_i = 0
    allocate(CS%col_j(CS%num_col))         ; CS%col_j = 0

    ! pass indices, restoring time to the CS structure
    col = 1
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if ((Iresttime(i,j)>0.0) .and. (G%mask2dT(i,j)>0)) then
        CS%col_i(col) = i ; CS%col_j(col) = j
        CS%Iresttime_col(col) = Iresttime(i,j)
        col = col +1
      endif
    enddo ; enddo

    ! same for total number of arbritary layers and correspondent data
    CS%nz_data = nz_data
    allocate(CS%Ref_h(CS%nz_data,CS%num_col))
    do col=1,CS%num_col ; do K=1,CS%nz_data
      CS%Ref_h(K,col) = data_h(CS%col_i(col),CS%col_j(col),K)
    enddo ; enddo

  endif

  total_sponge_cols = CS%num_col
  call sum_across_PEs(total_sponge_cols)

! Call the constructor for remapping control structure
  call initialize_remapping(CS%remap_cs, remapScheme, boundary_extrapolation=bndExtrapolation)

  call log_param(param_file, mod, "!Total sponge columns", total_sponge_cols, &
                 "The total number of columns where sponges are applied.")

end subroutine initialize_ALE_sponge

!> Initialize diagnostics for the ALE_sponge module.
! GMM: this routine is not being used for now.
subroutine init_ALE_sponge_diags(Time, G, diag, CS)
  type(time_type),       target, intent(in)    :: Time
  type(ocean_grid_type),         intent(in)    :: G
  type(diag_ctrl),       target, intent(inout) :: diag
  type(ALE_sponge_CS),           pointer       :: CS

  if (.not.associated(CS)) return

  CS%diag => diag

end subroutine init_ALE_sponge_diags

!> This subroutine stores the reference profile for the variable
! whose address is given by f_ptr.
subroutine set_up_ALE_sponge_field(sp_val, G, f_ptr, CS)
  type(ocean_grid_type),                   intent(in) :: G !< Grid structure (in).
  type(ALE_sponge_CS),                     pointer  :: CS !< Sponge structure (in/out).
  real, dimension(SZI_(G),SZJ_(G),CS%nz_data), intent(in) :: sp_val !< Field to be used in the sponge, it has arbritary number of layers (in).
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target, intent(in) :: f_ptr !< Pointer to the field to be damped (in).

  integer :: j, k, col
  character(len=256) :: mesg ! String for error messages

  if (.not.associated(CS)) return

  CS%fldno = CS%fldno + 1

  if (CS%fldno > MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ to at least ",I3," in MOM_memory.h or decrease &
           &the number of fields to be damped in the call to &
           &initialize_sponge." )') CS%fldno
    call MOM_error(FATAL,"set_up_ALE_sponge_field: "//mesg)
  endif

  ! stores the reference profile
  allocate(CS%Ref_val(CS%fldno)%p(CS%nz_data,CS%num_col))
  CS%Ref_val(CS%fldno)%p(:,:) = 0.0
  do col=1,CS%num_col
    do k=1,CS%nz_data
      CS%Ref_val(CS%fldno)%p(k,col) = sp_val(CS%col_i(col),CS%col_j(col),k)
    enddo
  enddo

  CS%var(CS%fldno)%p => f_ptr

end subroutine set_up_ALE_sponge_field

!> This subroutine applies damping to the layers thicknesses, temp, salt and a variety of tracers for every column where there is damping.
subroutine apply_ALE_sponge(h, dt, G, CS)
  type(ocean_grid_type),                    intent(inout) :: G !< The ocean's grid structure (in).
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h !< Layer thickness, in m (in)
  real,                                     intent(in)    :: dt !< The amount of time covered by this call, in s (in).
  type(ALE_sponge_CS),                      pointer       :: CS !<A pointer to the control structure for this module that is set by a previous call to initialize_sponge (in).

  real :: damp                                               !< The timestep times the local damping  coefficient.  ND.
  real :: I1pdamp                                            !< I1pdamp is 1/(1 + damp).  Nondimensional.
  real :: Idt                                                !< 1.0/dt, in s-1.
  real :: tmp_val1(cs%nz)                                    !< data values remapped to model grid
  real :: tmp_val2(cs%nz_data)                               ! data values remapped to model grid

  integer :: c, m, nkmb, i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not.associated(CS)) return

  do c=1,CS%num_col
! c is an index for the next 3 lines but a multiplier for the rest of the loop
! Therefore we use c as per C code and increment the index where necessary.
    i = CS%col_i(c) ; j = CS%col_j(c)
    damp = dt*CS%Iresttime_col(c)
    I1pdamp = 1.0 / (1.0 + damp)
        do m=1,CS%fldno
            tmp_val2(1:CS%nz_data) = CS%Ref_val(m)%p(1:CS%nz_data,c)
            call remapping_core_h(CS%nz_data, CS%Ref_h(:,c), tmp_val2, &
                  CS%nz, h(i,j,:), tmp_val1, CS%remap_cs)

!Backward Euler method
            CS%var(m)%p(i,j,:) = I1pdamp * &
                                 (CS%var(m)%p(i,j,:) + tmp_val1 * damp)

!            CS%var(m)%p(i,j,k) = I1pdamp * &
!                                 (CS%var(m)%p(i,j,k) + tmp_val1 * damp)
        enddo

  enddo ! end of c loop
  ! for debugging
  !c=CS%num_col
  !do m=1,CS%fldno
  !   write(*,*)'APPLY SPONGE,m,CS%Ref_h(:,c),h(i,j,:),tmp_val2,tmp_val1',m,CS%Ref_h(:,c),h(i,j,:),tmp_val2,tmp_val1
  !enddo
end subroutine apply_ALE_sponge

!> GMM: I could not find where sponge_end is being called, but I am keeping
!  ALE_sponge_end here so we can add that if needed.
subroutine ALE_sponge_end(CS)
  type(ALE_sponge_CS),              pointer       :: CS
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
