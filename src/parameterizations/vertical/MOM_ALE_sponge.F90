!> This module contains the routines used to apply sponge layers when using
!! the ALE mode.
!!
!! Applying sponges requires the following:
!! 1. initialize_ALE_sponge
!! 2. set_up_ALE_sponge_field (tracers) and set_up_ALE_sponge_vel_field (vel)
!! 3. apply_ALE_sponge
!! 4. init_ALE_sponge_diags (not being used for now)
!! 5. ALE_sponge_end (not being used for now)

module MOM_ALE_sponge


! This file is part of MOM6. See LICENSE.md for the license.
use MOM_array_transform, only: rotate_array
use MOM_coms,          only : sum_across_PEs
use MOM_diag_mediator, only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator, only : diag_ctrl
use MOM_domains, only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, NOTE, WARNING, is_root_pe
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_horizontal_regridding, only : horiz_interp_and_extrap_tracer
use MOM_interpolate,   only : init_external_field, get_external_field_info, time_interp_external_init
use MOM_remapping,     only : remapping_cs, remapping_core_h, initialize_remapping
use MOM_spatial_means, only : global_i_mean
use MOM_time_manager,  only : time_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_verticalGrid,  only : verticalGrid_type

use mpp_io_mod, only : mpp_get_axis_length
use mpp_io_mod, only : axistype

implicit none ; private

#include <MOM_memory.h>

!> Store the reference profile at h points for a variable
interface set_up_ALE_sponge_field
  module procedure set_up_ALE_sponge_field_fixed
  module procedure set_up_ALE_sponge_field_varying
end interface

!> This subroutine stores the reference profile at u and v points for a vector
interface set_up_ALE_sponge_vel_field
  module procedure set_up_ALE_sponge_vel_field_fixed
  module procedure set_up_ALE_sponge_vel_field_varying
end interface

!> Ddetermine the number of points which are within sponges in this computational domain.
!!
!! Only points that have positive values of Iresttime and which mask2dT indicates are ocean
!! points are included in the sponges.  It also stores the target interface heights.
interface initialize_ALE_sponge
  module procedure initialize_ALE_sponge_fixed
  module procedure initialize_ALE_sponge_varying
end interface

!  Publicly available functions
public set_up_ALE_sponge_field, set_up_ALE_sponge_vel_field
public get_ALE_sponge_thicknesses, get_ALE_sponge_nz_data
public initialize_ALE_sponge, apply_ALE_sponge, ALE_sponge_end, init_ALE_sponge_diags
public rotate_ALE_sponge, update_ALE_sponge_field

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A structure for creating arrays of pointers to 3D arrays with extra gridding information
type :: p3d
  integer :: id !< id for FMS external time interpolator
  integer :: nz_data !< The number of vertical levels in the input field.
  integer :: num_tlevs !< The number of time records contained in the file
  real, dimension(:,:,:), pointer :: mask_in => NULL() !< pointer to the data mask.
  real, dimension(:,:,:), pointer :: p => NULL() !< pointer to the data.
  real, dimension(:,:,:), pointer :: h => NULL() !< pointer to the data grid.
end type p3d

!> A structure for creating arrays of pointers to 2D arrays with extra gridding information
type :: p2d
  integer :: id !< id for FMS external time interpolator
  integer :: nz_data !< The number of vertical levels in the input field
  integer :: num_tlevs !< The number of time records contained in the file
  real, dimension(:,:), pointer :: mask_in => NULL()!< pointer to the data mask.
  real, dimension(:,:), pointer :: p => NULL() !< pointer the data.
  real, dimension(:,:), pointer :: h => NULL() !< pointer the data grid.
end type p2d

!> ALE sponge control structure
type, public :: ALE_sponge_CS ; private
  integer :: nz        !< The total number of layers.
  integer :: nz_data   !< The total number of arbritary layers (used by older code).
  integer :: isc       !< The starting i-index of the computational domain at h.
  integer :: iec       !< The ending i-index of the computational domain at h.
  integer :: jsc       !< The starting j-index of the computational domain at h.
  integer :: jec       !< The ending j-index of the computational domain at h.
  integer :: IscB      !< The starting I-index of the computational domain at u/v.
  integer :: IecB      !< The ending I-index of the computational domain at u/v.
  integer :: JscB      !< The starting J-index of the computational domain at u/v.
  integer :: JecB      !< The ending J-index of the computational domain at h.
  integer :: isd       !< The starting i-index of the data domain at h.
  integer :: ied       !< The ending i-index of the data domain at h.
  integer :: jsd       !< The starting j-index of the data domain at h.
  integer :: jed       !< The ending j-index of the data domain at h.
  integer :: num_col   !< The number of sponge tracer points within the computational domain.
  integer :: num_col_u !< The number of sponge u-points within the computational domain.
  integer :: num_col_v !< The number of sponge v-points within the computational domain.
  integer :: fldno = 0 !< The number of fields which have already been
                       !! registered by calls to set_up_sponge_field
  logical :: sponge_uv !< Control whether u and v are included in sponge
  integer, pointer :: col_i(:) => NULL()   !< Array of the i-indicies of each tracer columns being damped.
  integer, pointer :: col_j(:) => NULL()   !< Array of the j-indicies of each tracer columns being damped.
  integer, pointer :: col_i_u(:) => NULL() !< Array of the i-indicies of each u-columns being damped.
  integer, pointer :: col_j_u(:) => NULL() !< Array of the j-indicies of each u-columns being damped.
  integer, pointer :: col_i_v(:) => NULL() !< Array of the i-indicies of each v-columns being damped.
  integer, pointer :: col_j_v(:) => NULL() !< Array of the j-indicies of each v-columns being damped.

  real, pointer :: Iresttime_col(:)   => NULL() !< The inverse restoring time of each tracer column [T-1 ~> s-1].
  real, pointer :: Iresttime_col_u(:) => NULL() !< The inverse restoring time of each u-column [T-1 ~> s-1].
  real, pointer :: Iresttime_col_v(:) => NULL() !< The inverse restoring time of each v-column [T-1 ~> s-1].

  type(p3d) :: var(MAX_FIELDS_)      !< Pointers to the fields that are being damped.
  type(p2d) :: Ref_val(MAX_FIELDS_) !< The values to which the fields are damped.
  type(p2d) :: Ref_val_u  !< The values to which the u-velocities are damped.
  type(p2d) :: Ref_val_v  !< The values to which the v-velocities are damped.
  type(p3d) :: var_u  !< Pointer to the u velocities. that are being damped.
  type(p3d) :: var_v  !< Pointer to the v velocities. that are being damped.
  type(p2d) :: Ref_h  !< Grid on which reference data is provided (older code).
  type(p2d) :: Ref_hu !< u-point grid on which reference data is provided (older code).
  type(p2d) :: Ref_hv !< v-point grid on which reference data is provided (older code).

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  type(remapping_cs) :: remap_cs   !< Remapping parameters and work arrays
  logical :: remap_answers_2018    !< If true, use the order of arithmetic and expressions that
                                   !! recover the answers for remapping from the end of 2018.
                                   !! Otherwise, use more robust forms of the same expressions.
  logical :: hor_regrid_answers_2018 !< If true, use the order of arithmetic for horizonal regridding
                                   !! that recovers the answers from the end of 2018.  Otherwise, use
                                   !! rotationally symmetric forms of the same expressions.

  logical :: time_varying_sponges  !< True if using newer sponge code
  logical :: spongeDataOngrid !< True if the sponge data are on the model horizontal grid

  logical :: reentrant_x !< grid is reentrant in the x direction
  logical :: tripolar_N !< grid is folded at its north edge

  !>@{ Diagnostic IDs
  integer, dimension(2) :: id_sp_tendency      !< Diagnostic ids for temperature and salinity
                                               !! tendency due to sponges
  integer :: id_sp_u_tendency                  !< Diagnostic id for zonal momentum tendency due to
                                               !! Rayleigh damping
  integer :: id_sp_v_tendency                  !< Diagnostic id for meridional momentum tendency due to
                                               !! Rayleigh damping
end type ALE_sponge_CS

contains

!> This subroutine determines the number of points which are within sponges in this computational
!! domain.  Only points that have positive values of Iresttime and which mask2dT indicates are ocean
!! points are included in the sponges.  It also stores the target interface heights. This
subroutine initialize_ALE_sponge_fixed(Iresttime, G, GV, param_file, CS, data_h, nz_data, &
                                        Iresttime_u_in, Iresttime_v_in)

  type(ocean_grid_type),            intent(in) :: G !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  integer,                          intent(in) :: nz_data !< The total number of sponge input layers.
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: Iresttime !< The inverse of the restoring time [T-1 ~> s-1].
  type(param_file_type),            intent(in) :: param_file !< A structure indicating the open file
                                                             !! to parse for model parameter values.
  type(ALE_sponge_CS),              pointer    :: CS !< A pointer that is set to point to the control
                                                     !! structure for this module (in/out).
  real, dimension(SZI_(G),SZJ_(G),nz_data), intent(in) :: data_h !< The thicknesses of the sponge
                                                     !! input layers [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), intent(in), optional :: Iresttime_u_in  !< The inverse of the restoring
                                                                             !! time at U-points [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in), optional :: Iresttime_v_in  !< The inverse of the restoring
                                                                             ! time at v-points [T-1 ~> s-1].


! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_sponge"  ! This module's name.
  logical :: use_sponge
  real, allocatable, dimension(:,:) :: Iresttime_u !< inverse of the restoring time at u points [T-1 ~> s-1]
  real, allocatable, dimension(:,:) :: Iresttime_v !< inverse of the restoring time at v points [T-1 ~> s-1]
  logical :: bndExtrapolation = .true. ! If true, extrapolate boundaries
  logical :: default_2018_answers
  integer :: i, j, k, col, total_sponge_cols, total_sponge_cols_u, total_sponge_cols_v
  character(len=10)  :: remapScheme
  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_ALE_sponge_fixed called with an associated "// &
                            "control structure.")
    return
  endif

! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SPONGE", use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  if (.not.use_sponge) return

  allocate(CS)

  call get_param(param_file, mdl, "SPONGE_UV", CS%sponge_uv, &
                 "Apply sponges in u and v, in addition to tracers.", &
                 default=.false.)

  call get_param(param_file, mdl, "REMAPPING_SCHEME", remapScheme, &
                 "This sets the reconstruction scheme used "//&
                 " for vertical remapping for all variables.", &
                 default="PLM", do_not_log=.true.)

  call get_param(param_file, mdl, "BOUNDARY_EXTRAPOLATION", bndExtrapolation, &
                 "When defined, a proper high-order reconstruction "//&
                 "scheme is used within boundary cells rather "//&
                 "than PCM. E.g., if PPM is used for remapping, a "//&
                 "PPM reconstruction will also be used within boundary cells.", &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", CS%remap_answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call get_param(param_file, mdl, "HOR_REGRID_2018_ANSWERS", CS%hor_regrid_answers_2018, &
                 "If true, use the order of arithmetic for horizonal regridding that recovers "//&
                 "the answers from the end of 2018.  Otherwise, use rotationally symmetric "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call get_param(param_file, mdl, "REENTRANT_X", CS%reentrant_x, &
                 "If true, the domain is zonally reentrant.", default=.true.)
  call get_param(param_file, mdl, "TRIPOLAR_N", CS%tripolar_N, &
                 "Use tripolar connectivity at the northern edge of the "//&
                 "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", default=.false.)

  CS%time_varying_sponges = .false.
  CS%nz = GV%ke
  CS%isc = G%isc ; CS%iec = G%iec ; CS%jsc = G%jsc ; CS%jec = G%jec
  CS%isd = G%isd ; CS%ied = G%ied ; CS%jsd = G%jsd ; CS%jed = G%jed
  CS%iscB = G%iscB ; CS%iecB = G%iecB; CS%jscB = G%jscB ; CS%jecB = G%jecB

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
    allocate(CS%Ref_h%p(CS%nz_data,CS%num_col))
    do col=1,CS%num_col ; do K=1,CS%nz_data
      CS%Ref_h%p(K,col) = data_h(CS%col_i(col),CS%col_j(col),K)
    enddo ; enddo
  endif

  total_sponge_cols = CS%num_col
  call sum_across_PEs(total_sponge_cols)

! Call the constructor for remapping control structure
  call initialize_remapping(CS%remap_cs, remapScheme, boundary_extrapolation=bndExtrapolation, &
                            answers_2018=CS%remap_answers_2018)

  call log_param(param_file, mdl, "!Total sponge columns at h points", total_sponge_cols, &
                 "The total number of columns where sponges are applied at h points.", like_default=.true.)

  if (CS%sponge_uv) then
    allocate(Iresttime_u(G%isdB:G%iedB,G%jsd:G%jed)) ; Iresttime_u(:,:) = 0.0
    allocate(Iresttime_v(G%isd:G%ied,G%jsdB:G%jedB)) ; Iresttime_v(:,:) = 0.0

    ! u points
    CS%num_col_u = 0 ;
    if (present(Iresttime_u_in)) then
       Iresttime_u(:,:) = Iresttime_u_in(:,:)
    else
      do j=CS%jsc,CS%jec ; do I=CS%iscB,CS%iecB
        Iresttime_u(I,j) = 0.5 * (Iresttime(i,j) + Iresttime(i+1,j))
      enddo ; enddo
    endif
    do j=CS%jsc,CS%jec ; do I=CS%iscB,CS%iecB
       if ((Iresttime_u(I,j)>0.0) .and. (G%mask2dCu(I,j)>0)) &
          CS%num_col_u = CS%num_col_u + 1
    enddo ; enddo

    if (CS%num_col_u > 0) then

      allocate(CS%Iresttime_col_u(CS%num_col_u)) ; CS%Iresttime_col_u(:) = 0.0
      allocate(CS%col_i_u(CS%num_col_u))         ; CS%col_i_u(:) = 0
      allocate(CS%col_j_u(CS%num_col_u))         ; CS%col_j_u(:) = 0

      ! Store the column indices and restoring rates in the CS structure
      col = 1
      do j=CS%jsc,CS%jec ; do I=CS%iscB,CS%iecB
        if ((Iresttime_u(I,j)>0.0) .and. (G%mask2dCu(I,j)>0)) then
          CS%col_i_u(col) = I ; CS%col_j_u(col) = j
          CS%Iresttime_col_u(col) = Iresttime_u(I,j)
          col = col + 1
        endif
      enddo ; enddo

      ! same for total number of arbritary layers and correspondent data
      allocate(CS%Ref_hu%p(CS%nz_data,CS%num_col_u))
      do col=1,CS%num_col_u
        I = CS%col_i_u(col) ; j = CS%col_j_u(col)
        do k=1,CS%nz_data
          CS%Ref_hu%p(k,col) = 0.5 * (data_h(i,j,k) + data_h(i+1,j,k))
        enddo
      enddo
    endif
    total_sponge_cols_u = CS%num_col_u
    call sum_across_PEs(total_sponge_cols_u)
    call log_param(param_file, mdl, "!Total sponge columns at u points", total_sponge_cols_u, &
                "The total number of columns where sponges are applied at u points.", like_default=.true.)

    ! v points
    CS%num_col_v = 0 ;
    if (present(Iresttime_v_in)) then
      Iresttime_v(:,:) = Iresttime_v_in(:,:)
    else
      do J=CS%jscB,CS%jecB; do i=CS%isc,CS%iec
        Iresttime_v(i,J) = 0.5 * (Iresttime(i,j) + Iresttime(i,j+1))
      enddo ; enddo
    endif
    do J=CS%jscB,CS%jecB; do i=CS%isc,CS%iec
      if ((Iresttime_v(i,J)>0.0) .and. (G%mask2dCv(i,J)>0)) &
        CS%num_col_v = CS%num_col_v + 1
    enddo ; enddo

    if (CS%num_col_v > 0) then

      allocate(CS%Iresttime_col_v(CS%num_col_v)) ; CS%Iresttime_col_v = 0.0
      allocate(CS%col_i_v(CS%num_col_v))         ; CS%col_i_v = 0
      allocate(CS%col_j_v(CS%num_col_v))         ; CS%col_j_v = 0

      ! pass indices, restoring time to the CS structure
      col = 1
      do J=CS%jscB,CS%jecB ; do i=CS%isc,CS%iec
        if ((Iresttime_v(i,J)>0.0) .and. (G%mask2dCv(i,J)>0)) then
          CS%col_i_v(col) = i ; CS%col_j_v(col) = j
          CS%Iresttime_col_v(col) = Iresttime_v(i,j)
          col = col + 1
        endif
      enddo ; enddo

      ! same for total number of arbritary layers and correspondent data
      allocate(CS%Ref_hv%p(CS%nz_data,CS%num_col_v))
      do col=1,CS%num_col_v
        i = CS%col_i_v(col) ; J = CS%col_j_v(col)
        do k=1,CS%nz_data
          CS%Ref_hv%p(k,col) = 0.5 * (data_h(i,j,k) + data_h(i,j+1,k))
        enddo
      enddo
    endif
    total_sponge_cols_v = CS%num_col_v
    call sum_across_PEs(total_sponge_cols_v)
    call log_param(param_file, mdl, "!Total sponge columns at v points", total_sponge_cols_v, &
                 "The total number of columns where sponges are applied at v points.", like_default=.true.)
  endif

end subroutine initialize_ALE_sponge_fixed

!> Return the number of layers in the data with a fixed ALE sponge, or 0 if there are
!! no sponge columns on this PE.
function get_ALE_sponge_nz_data(CS)
  type(ALE_sponge_CS),   pointer       :: CS !< A pointer that is set to point to the control
                                             !! structure for the ALE_sponge module.
  integer :: get_ALE_sponge_nz_data  !< The number of layers in the fixed sponge data.

  if (associated(CS)) then
    get_ALE_sponge_nz_data = CS%nz_data
  else
    get_ALE_sponge_nz_data = 0
  endif
end function get_ALE_sponge_nz_data

!> Return the thicknesses used for the data with a fixed ALE sponge
subroutine get_ALE_sponge_thicknesses(G, data_h, sponge_mask, CS)
  type(ocean_grid_type), intent(in)    :: G !< The ocean's grid structure (in).
  real, allocatable, dimension(:,:,:), &
                         intent(inout) :: data_h !< The thicknesses of the sponge input layers [H ~> m or kg m-2].
  logical, dimension(SZI_(G),SZJ_(G)), &
                         intent(out)   :: sponge_mask !< A logical mask that is true where
                                                 !! sponges are being applied.
  type(ALE_sponge_CS),   pointer       :: CS !< A pointer that is set to point to the control
                                             !! structure for the ALE_sponge module.
  integer :: c, i, j, k

  if (allocated(data_h)) call MOM_error(FATAL, &
    "get_ALE_sponge_thicknesses called with an allocated data_h.")

  if (.not.associated(CS)) then
    ! There are no sponge points on this PE.
    allocate(data_h(G%isd:G%ied,G%jsd:G%jed,1)) ; data_h(:,:,:) = -1.0
    sponge_mask(:,:) = .false.
    return
  endif

  allocate(data_h(G%isd:G%ied,G%jsd:G%jed,CS%nz_data)) ; data_h(:,:,:) = -1.0
  sponge_mask(:,:) = .false.

  do c=1,CS%num_col
    i = CS%col_i(c) ; j = CS%col_j(c)
    sponge_mask(i,j) = .true.
    do k=1,CS%nz_data
      data_h(i,j,k) = CS%Ref_h%p(k,c)
    enddo
  enddo

end subroutine get_ALE_sponge_thicknesses

!> This subroutine determines the number of points which are to be restoref in the computational
!! domain.  Only points that have positive values of Iresttime and which mask2dT indicates are ocean
!! points are included in the sponges.
subroutine initialize_ALE_sponge_varying(Iresttime, G, GV, param_file, CS, Iresttime_u_in, Iresttime_v_in)

  type(ocean_grid_type),            intent(in) :: G !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: Iresttime !< The inverse of the restoring time [T-1 ~> s-1].
  type(param_file_type),            intent(in) :: param_file !< A structure indicating the open file to parse
                                                             !! for model parameter values.
  type(ALE_sponge_CS),              pointer    :: CS !< A pointer that is set to point to the control
                                                     !! structure for this module (in/out).
  real, dimension(SZIB_(G),SZJ_(G)), intent(in), optional :: Iresttime_u_in !< The inverse of the restoring time
                                                                            !! for u [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in), optional :: Iresttime_v_in !< The inverse of the restoring time
                                                                            !! for v [T-1 ~> s-1].

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_sponge"  ! This module's name.
  logical :: use_sponge
  real, allocatable, dimension(:,:) :: Iresttime_u !< inverse of the restoring time at u points [T-1 ~> s-1]
  real, allocatable, dimension(:,:) :: Iresttime_v !< inverse of the restoring time at v points [T-1 ~> s-1]
  logical :: bndExtrapolation = .true. ! If true, extrapolate boundaries
  logical :: default_2018_answers
  logical :: spongeDataOngrid = .false.
  integer :: i, j, k, col, total_sponge_cols, total_sponge_cols_u, total_sponge_cols_v
  character(len=10)  :: remapScheme

  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_ALE_sponge_varying called with an associated "// &
                            "control structure.")
    return
  endif
! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SPONGE", use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)
  if (.not.use_sponge) return
  allocate(CS)
  call get_param(param_file, mdl, "SPONGE_UV", CS%sponge_uv, &
                 "Apply sponges in u and v, in addition to tracers.", &
                 default=.false.)
  call get_param(param_file, mdl, "REMAPPING_SCHEME", remapScheme, &
                 "This sets the reconstruction scheme used "//&
                 " for vertical remapping for all variables.", &
                 default="PLM", do_not_log=.true.)
  call get_param(param_file, mdl, "BOUNDARY_EXTRAPOLATION", bndExtrapolation, &
                 "When defined, a proper high-order reconstruction "//&
                 "scheme is used within boundary cells rather "//&
                 "than PCM. E.g., if PPM is used for remapping, a "//&
                 "PPM reconstruction will also be used within boundary cells.", &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", CS%remap_answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call get_param(param_file, mdl, "HOR_REGRID_2018_ANSWERS", CS%hor_regrid_answers_2018, &
                 "If true, use the order of arithmetic for horizonal regridding that recovers "//&
                 "the answers from the end of 2018 and retain a bug in the 3-dimensional mask "//&
                 "returned in certain cases.  Otherwise, use rotationally symmetric "//&
                 "forms of the same expressions and initialize the mask properly.", &
                 default=default_2018_answers)
  call get_param(param_file, mdl, "SPONGE_DATA_ONGRID", CS%spongeDataOngrid, &
                 "When defined, the incoming sponge data are "//&
                 "assumed to be on the model grid " , &
                 default=.false.)
  call get_param(param_file, mdl, "REENTRANT_X", CS%reentrant_x, &
                 "If true, the domain is zonally reentrant.", default=.true.)
  call get_param(param_file, mdl, "TRIPOLAR_N", CS%tripolar_N, &
                 "Use tripolar connectivity at the northern edge of the "//&
                 "domain.  With TRIPOLAR_N, NIGLOBAL must be even.", default=.false.)

  CS%time_varying_sponges = .true.
  CS%nz = GV%ke
  CS%isc = G%isc ; CS%iec = G%iec ; CS%jsc = G%jsc ; CS%jec = G%jec
  CS%isd = G%isd ; CS%ied = G%ied ; CS%jsd = G%jsd ; CS%jed = G%jed
  CS%iscB = G%iscB ; CS%iecB = G%iecB; CS%jscB = G%jscB ; CS%jecB = G%jecB

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
        col = col + 1
      endif
    enddo ; enddo
  endif
  total_sponge_cols = CS%num_col
  call sum_across_PEs(total_sponge_cols)

! Call the constructor for remapping control structure
  call initialize_remapping(CS%remap_cs, remapScheme, boundary_extrapolation=bndExtrapolation, &
                            answers_2018=CS%remap_answers_2018)
  call log_param(param_file, mdl, "!Total sponge columns at h points", total_sponge_cols, &
                 "The total number of columns where sponges are applied at h points.", like_default=.true.)
  if (CS%sponge_uv) then
    allocate(Iresttime_u(G%isdB:G%iedB,G%jsd:G%jed)) ; Iresttime_u(:,:) = 0.0
    allocate(Iresttime_v(G%isd:G%ied,G%jsdB:G%jedB)) ; Iresttime_v(:,:) = 0.0

    call pass_var(Iresttime,G%Domain)
    ! u points
    if (present(Iresttime_u_in)) then
      Iresttime_u(:,:) = Iresttime_u_in(:,:)
    else
      do j=CS%jsc,CS%jec ; do I=CS%iscB,CS%iecB
        Iresttime_u(I,j) = 0.5 * (Iresttime(i,j) + Iresttime(i+1,j))
      enddo ; enddo
    endif
    CS%num_col_u = 0 ;
    do j=CS%jsc,CS%jec; do I=CS%iscB,CS%iecB
      if ((Iresttime_u(I,j)>0.0) .and. (G%mask2dCu(I,j)>0)) &
        CS%num_col_u = CS%num_col_u + 1
    enddo ; enddo
    if (CS%num_col_u > 0) then
      allocate(CS%Iresttime_col_u(CS%num_col_u)) ; CS%Iresttime_col_u = 0.0
      allocate(CS%col_i_u(CS%num_col_u))         ; CS%col_i_u = 0
      allocate(CS%col_j_u(CS%num_col_u))         ; CS%col_j_u = 0
      ! pass indices, restoring time to the CS structure
      col = 1
      do j=CS%jsc,CS%jec ; do I=CS%iscB,CS%iecB
        if ((Iresttime_u(I,j)>0.0) .and. (G%mask2dCu(I,j)>0)) then
          CS%col_i_u(col) = i ; CS%col_j_u(col) = j
          CS%Iresttime_col_u(col) = Iresttime_u(i,j)
          col = col + 1
        endif
      enddo ; enddo
      ! same for total number of arbritary layers and correspondent data
    endif
    total_sponge_cols_u = CS%num_col_u
    call sum_across_PEs(total_sponge_cols_u)
    call log_param(param_file, mdl, "!Total sponge columns at u points", total_sponge_cols_u, &
                "The total number of columns where sponges are applied at u points.", like_default=.true.)
    ! v points
    if (present(Iresttime_v_in)) then
      Iresttime_v(:,:) = Iresttime_v_in(:,:)
    else
      do J=CS%jscB,CS%jecB; do i=CS%isc,CS%iec
        Iresttime_v(i,J) = 0.5 * (Iresttime(i,j) + Iresttime(i,j+1))
      enddo ; enddo
    endif
    CS%num_col_v = 0 ;
    do J=CS%jscB,CS%jecB; do i=CS%isc,CS%iec
      if ((Iresttime_v(i,J)>0.0) .and. (G%mask2dCv(i,J)>0)) &
        CS%num_col_v = CS%num_col_v + 1
    enddo ; enddo
    if (CS%num_col_v > 0) then
      allocate(CS%Iresttime_col_v(CS%num_col_v)) ; CS%Iresttime_col_v = 0.0
      allocate(CS%col_i_v(CS%num_col_v))         ; CS%col_i_v = 0
      allocate(CS%col_j_v(CS%num_col_v))         ; CS%col_j_v = 0
      ! pass indices, restoring time to the CS structure
      col = 1
      do J=CS%jscB,CS%jecB ; do i=CS%isc,CS%iec
        if ((Iresttime_v(i,J)>0.0) .and. (G%mask2dCv(i,J)>0)) then
          CS%col_i_v(col) = i ; CS%col_j_v(col) = j
          CS%Iresttime_col_v(col) = Iresttime_v(i,j)
          col = col + 1
        endif
      enddo ; enddo
    endif
    total_sponge_cols_v = CS%num_col_v
    call sum_across_PEs(total_sponge_cols_v)
    call log_param(param_file, mdl, "!Total sponge columns at v points", total_sponge_cols_v, &
                "The total number of columns where sponges are applied at v points.", like_default=.true.)
  endif

end subroutine initialize_ALE_sponge_varying

!> Initialize diagnostics for the ALE_sponge module.
! GMM: this routine is not being used for now.
subroutine init_ALE_sponge_diags(Time, G, diag, CS, US)
  type(time_type), target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(ALE_sponge_CS),     pointer       :: CS   !< ALE sponge control structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type

  if (.not.associated(CS)) return

  CS%diag => diag

  CS%id_sp_tendency(1) = -1
  CS%id_sp_tendency(1) = register_diag_field('ocean_model', 'sp_tendency_temp', diag%axesTL, Time, &
       'Time tendency due to temperature restoring', 'degC s-1',conversion=US%s_to_T)
  CS%id_sp_tendency(2) = -1
  CS%id_sp_tendency(2) = register_diag_field('ocean_model', 'sp_tendency_salt', diag%axesTL, Time, &
       'Time tendency due to salinity restoring', 'g kg-1 s-1',conversion=US%s_to_T)
  CS%id_sp_u_tendency = -1
  CS%id_sp_u_tendency = register_diag_field('ocean_model', 'sp_tendency_u', diag%axesCuL, Time, &
       'Zonal acceleration due to sponges', 'm s-2',conversion=US%L_T2_to_m_s2)
  CS%id_sp_v_tendency = -1
  CS%id_sp_v_tendency = register_diag_field('ocean_model', 'sp_tendency_v', diag%axesCvL, Time, &
       'Meridional acceleration due to sponges', 'm s-2',conversion=US%L_T2_to_m_s2)

end subroutine init_ALE_sponge_diags

!> This subroutine stores the reference profile at h points for the variable
!! whose address is given by f_ptr.
subroutine set_up_ALE_sponge_field_fixed(sp_val, G, GV, f_ptr, CS)
  type(ocean_grid_type),   intent(in) :: G  !< Grid structure
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  type(ALE_sponge_CS),     pointer    :: CS !< ALE sponge control structure (in/out).
  real, dimension(SZI_(G),SZJ_(G),CS%nz_data), &
                           intent(in) :: sp_val !< Field to be used in the sponge, it can have an
                                            !! arbitrary number of layers.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                   target, intent(in) :: f_ptr !< Pointer to the field to be damped

  integer :: j, k, col
  character(len=256) :: mesg ! String for error messages

  if (.not.associated(CS)) return

  CS%fldno = CS%fldno + 1
  if (CS%fldno > MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ to at least ",I3," in MOM_memory.h or decrease &
           &the number of fields to be damped in the call to &
           &initialize_ALE_sponge." )') CS%fldno
    call MOM_error(FATAL,"set_up_ALE_sponge_field: "//mesg)
  endif

  ! stores the reference profile
  CS%Ref_val(CS%fldno)%nz_data = CS%nz_data
  allocate(CS%Ref_val(CS%fldno)%p(CS%nz_data,CS%num_col))
  CS%Ref_val(CS%fldno)%p(:,:) = 0.0
  do col=1,CS%num_col
    do k=1,CS%nz_data
      CS%Ref_val(CS%fldno)%p(k,col) = sp_val(CS%col_i(col),CS%col_j(col),k)
    enddo
  enddo

  CS%var(CS%fldno)%p => f_ptr

end subroutine set_up_ALE_sponge_field_fixed

!> This subroutine stores the reference profile at h points for the variable
!! whose address is given by filename and fieldname.
subroutine set_up_ALE_sponge_field_varying(filename, fieldname, Time, G, GV, US, f_ptr, CS)
  character(len=*),        intent(in) :: filename !< The name of the file with the
                                                  !! time varying field data
  character(len=*),        intent(in) :: fieldname !< The name of the field in the file
                                                  !! with the time varying field data
  type(time_type),         intent(in) :: Time  !< The current model time
  type(ocean_grid_type),   intent(in) :: G     !< Grid structure (in).
  type(verticalGrid_type), intent(in) :: GV    !< ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                   target, intent(in) :: f_ptr !< Pointer to the field to be damped (in).
  type(ALE_sponge_CS),     pointer    :: CS    !< Sponge control structure (in/out).

  ! Local variables
  real, allocatable, dimension(:,:,:) :: sp_val !< Field to be used in the sponge
  real, allocatable, dimension(:,:,:) :: mask_z !< Field mask for the sponge data
  real, allocatable, dimension(:), target :: z_in, z_edges_in ! Heights [Z ~> m].
  real :: missing_value
  integer :: j, k, col
  integer :: isd,ied,jsd,jed
  integer :: nPoints
  integer, dimension(4) :: fld_sz
  integer :: nz_data !< the number of vertical levels in this input field
  character(len=256) :: mesg ! String for error messages
  ! Local variables for ALE remapping
  real, dimension(:), allocatable :: tmpT1d
  real :: zTopOfCell, zBottomOfCell ! Heights [Z ~> m].
  type(remapping_CS) :: remapCS ! Remapping parameters and work arrays

  if (.not.associated(CS)) return
  ! initialize time interpolator module
  call time_interp_external_init()
  isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
  CS%fldno = CS%fldno + 1
  if (CS%fldno > MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ to at least ",I3," in MOM_memory.h or decrease &
           &the number of fields to be damped in the call to &
           &initialize_ALE_sponge." )') CS%fldno
    call MOM_error(FATAL,"set_up_ALE_sponge_field: "//mesg)
  endif
  ! get a unique time interp id for this field. If sponge data is ongrid, then setup
  ! to only read on the computational domain
  if (CS%spongeDataOngrid) then
    CS%Ref_val(CS%fldno)%id = init_external_field(filename, fieldname, MOM_domain=G%Domain)
  else
    CS%Ref_val(CS%fldno)%id = init_external_field(filename, fieldname)
  endif
  fld_sz(1:4)=-1
  call get_external_field_info(CS%Ref_val(CS%fldno)%id, size=fld_sz)
  nz_data = fld_sz(3)
  CS%Ref_val(CS%fldno)%nz_data = nz_data !< individual sponge fields may reside on a different vertical grid
  CS%Ref_val(CS%fldno)%num_tlevs = fld_sz(4)
  ! initializes the target profile array for this field
  ! for all columns which will be masked
  allocate(CS%Ref_val(CS%fldno)%p(nz_data,CS%num_col))
  CS%Ref_val(CS%fldno)%p(:,:) = 0.0
  allocate( CS%Ref_val(CS%fldno)%h(nz_data,CS%num_col) )
  CS%Ref_val(CS%fldno)%h(:,:) = 0.0
  CS%var(CS%fldno)%p => f_ptr

end subroutine set_up_ALE_sponge_field_varying

!> This subroutine stores the reference profile at u and v points for the variable
!! whose address is given by u_ptr and v_ptr.
subroutine set_up_ALE_sponge_vel_field_fixed(u_val, v_val, G, GV, u_ptr, v_ptr, CS)
  type(ocean_grid_type),   intent(in) :: G  !< Grid structure (in).
  type(verticalGrid_type), intent(in) :: GV    !< ocean vertical grid structure
  type(ALE_sponge_CS),     pointer    :: CS !< Sponge structure (in/out).
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                          intent(in) :: u_val !< u field to be used in the sponge, it has arbritary number of layers but
                                              !! not to exceed the total number of model layers
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                          intent(in) :: v_val !< v field to be used in the sponge, it has arbritary number of layers but
                                              !! not to exceed the number of model layers
  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u_ptr !< u pointer to the field to be damped
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v_ptr !< v pointer to the field to be damped

  integer :: j, k, col, fld_sz(4)
  character(len=256) :: mesg ! String for error messages

  if (.not.associated(CS)) return

  ! stores the reference profile
  allocate(CS%Ref_val_u%p(CS%nz_data,CS%num_col_u))
  CS%Ref_val_u%p(:,:) = 0.0
  do col=1,CS%num_col_u
    do k=1,CS%nz_data
      CS%Ref_val_u%p(k,col) = u_val(CS%col_i_u(col),CS%col_j_u(col),k)
    enddo
  enddo
  CS%var_u%p => u_ptr
  allocate(CS%Ref_val_v%p(CS%nz_data,CS%num_col_v))
  CS%Ref_val_v%p(:,:) = 0.0
  do col=1,CS%num_col_v
    do k=1,CS%nz_data
      CS%Ref_val_v%p(k,col) = v_val(CS%col_i_v(col),CS%col_j_v(col),k)
    enddo
  enddo
  CS%var_v%p => v_ptr

end subroutine set_up_ALE_sponge_vel_field_fixed

!> This subroutine stores the reference profile at u and v points for the variable
!! whose address is given by u_ptr and v_ptr.
subroutine set_up_ALE_sponge_vel_field_varying(filename_u, fieldname_u, filename_v, fieldname_v, &
                                               Time, G, GV, US, CS, u_ptr, v_ptr)
  character(len=*), intent(in)    :: filename_u  !< File name for u field
  character(len=*), intent(in)    :: fieldname_u !< Name of u variable in file
  character(len=*), intent(in)    :: filename_v  !< File name for v field
  character(len=*), intent(in)    :: fieldname_v !< Name of v variable in file
  type(time_type),  intent(in)    :: Time        !< Model time
  type(ocean_grid_type), intent(in) :: G         !< Ocean grid (in)
  type(verticalGrid_type), intent(in) :: GV    !< ocean vertical grid structure
  type(unit_scale_type), intent(in)    :: US     !< A dimensional unit scaling type
  type(ALE_sponge_CS), pointer    :: CS          !< Sponge structure (in/out).
  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u_ptr !< u pointer to the field to be damped (in).
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v_ptr !< v pointer to the field to be damped (in).
  ! Local variables
  real, allocatable, dimension(:,:,:) :: u_val !< U field to be used in the sponge.
  real, allocatable, dimension(:,:,:) :: v_val !< V field to be used in the sponge.

  real, allocatable, dimension(:), target :: z_in, z_edges_in
  real :: missing_value
  logical :: override

  integer :: j, k, col
  integer :: isd, ied, jsd, jed
  integer :: isdB, iedB, jsdB, jedB
  integer, dimension(4) :: fld_sz
  character(len=256) :: mesg ! String for error messages
  type(axistype), dimension(4) :: axes_data
  integer :: tmp
  integer :: axis_sizes(4)
  if (.not.associated(CS)) return

  override =.true.

  isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
  isdB = G%isdB; iedB = G%iedB; jsdB = G%jsdB; jedB = G%jedB
  ! get a unique id for this field which will allow us to return an array
  ! containing time-interpolated values from an external file corresponding
  ! to the current model date.
  CS%Ref_val_u%id = init_external_field(filename_u, fieldname_u, domain=G%Domain%mpp_domain)
  fld_sz(1:4)=-1
  call get_external_field_info(CS%Ref_val_u%id, size=fld_sz)
  CS%Ref_val_u%nz_data = fld_sz(3)
  CS%Ref_val_u%num_tlevs = fld_sz(4)

  CS%Ref_val_v%id = init_external_field(filename_v, fieldname_v, domain=G%Domain%mpp_domain)
  fld_sz(1:4)=-1
  call get_external_field_info(CS%Ref_val_v%id, size=fld_sz)
  CS%Ref_val_v%nz_data = fld_sz(3)
  CS%Ref_val_v%num_tlevs = fld_sz(4)

  ! stores the reference profile
  allocate(CS%Ref_val_u%p(fld_sz(3),CS%num_col_u))
  CS%Ref_val_u%p(:,:) = 0.0
  allocate(CS%Ref_val_u%h(fld_sz(3),CS%num_col_u) )
  CS%Ref_val_u%h(:,:) = 0.0
  CS%var_u%p => u_ptr
  allocate(CS%Ref_val_v%p(fld_sz(3),CS%num_col_v))
  CS%Ref_val_v%p(:,:) = 0.0
  allocate(CS%Ref_val_v%h(fld_sz(3),CS%num_col_v) )
  CS%Ref_val_v%h(:,:) = 0.0
  CS%var_v%p => v_ptr

end subroutine set_up_ALE_sponge_vel_field_varying

!> This subroutine applies damping to the layers thicknesses, temp, salt and a variety of tracers
!! for every column where there is damping.
subroutine apply_ALE_sponge(h, dt, G, GV, US, CS, Time)
  type(ocean_grid_type),     intent(inout) :: G  !< The ocean's grid structure (in).
  type(verticalGrid_type),   intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2] (in)
  real,                      intent(in)    :: dt !< The amount of time covered by this call [T ~> s].
  type(ALE_sponge_CS),       pointer       :: CS !< A pointer to the control structure for this module
                                                 !! that is set by a previous call to initialize_ALE_sponge (in).
  type(time_type), optional, intent(in)    :: Time !< The current model date

  real :: damp                                  ! The timestep times the local damping coefficient [nondim].
  real :: I1pdamp                               ! I1pdamp is 1/(1 + damp). [nondim].
  real :: m_to_Z                                ! A unit conversion factor from m to Z.
  real, allocatable, dimension(:) :: tmp_val2   ! data values on the original grid
  real, dimension(SZK_(GV)) :: tmp_val1         ! data values remapped to model grid
  real, dimension(SZK_(GV)) :: h_col            ! A column of thicknesses at h, u or v points [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:) :: sp_val ! A temporary array for fields
  real, allocatable, dimension(:,:,:) :: sp_val_u ! A temporary array for fields
  real, allocatable, dimension(:,:,:) :: sp_val_v ! A temporary array for fields
  real, allocatable, dimension(:,:,:) :: mask_z ! A temporary array for field mask at h pts
  real, allocatable, dimension(:,:,:) :: tmp    !< A temporary array for thermodynamic sponge tendency diagnostics,
  real, allocatable, dimension(:,:,:) :: tmp_u  !< A temporary array for u sponge acceleration diagnostics
  real, allocatable, dimension(:,:,:) :: tmp_v  !< A temporary array for v sponge acceleration diagnostics
  real, dimension(:), allocatable :: hsrc       ! Source thicknesses [Z ~> m].
  ! Local variables for ALE remapping
  real, dimension(:), allocatable :: tmpT1d
  integer :: c, m, nkmb, i, j, k, is, ie, js, je, nz, nz_data
  integer :: col, total_sponge_cols
  real, allocatable, dimension(:), target :: z_in, z_edges_in
  real :: missing_value, Idt
  real :: h_neglect, h_neglect_edge
  real :: zTopOfCell, zBottomOfCell ! Heights [Z ~> m].
  integer :: nPoints

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  if (.not.associated(CS)) return

  Idt = 1.0/dt

  if (.not.CS%remap_answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  if (CS%time_varying_sponges) then
    if (.not. present(Time)) &
      call MOM_error(FATAL,"apply_ALE_sponge: No time information provided")
    do m=1,CS%fldno
      nz_data = CS%Ref_val(m)%nz_data
      allocate(sp_val(G%isd:G%ied,G%jsd:G%jed,1:nz_data)); sp_val(:,:,:) = 0.0
      allocate(mask_z(G%isd:G%ied,G%jsd:G%jed,1:nz_data)); mask_z(:,:,:) = 0.0
      call horiz_interp_and_extrap_tracer(CS%Ref_val(m)%id, Time, 1.0, G, sp_val, mask_z, z_in, &
                     z_edges_in,  missing_value, CS%reentrant_x, CS%tripolar_N, .false., &
                     spongeOnGrid=CS%SpongeDataOngrid, m_to_Z=US%m_to_Z, &
                     answers_2018=CS%hor_regrid_answers_2018)
      allocate( hsrc(nz_data) )
      allocate( tmpT1d(nz_data) )
      do c=1,CS%num_col
        i = CS%col_i(c) ; j = CS%col_j(c)
        CS%Ref_val(m)%p(1:nz_data,c) = sp_val(i,j,1:nz_data)
        ! Build the source grid
        zTopOfCell = 0. ; zBottomOfCell = 0. ; nPoints = 0; hsrc(:) = 0.0; tmpT1d(:) = -99.9
        do k=1,nz_data
          if (mask_z(CS%col_i(c),CS%col_j(c),k) == 1.0) then
            zBottomOfCell = -min( z_edges_in(k+1), G%bathyT(CS%col_i(c),CS%col_j(c)) )
            tmpT1d(k) = sp_val(CS%col_i(c),CS%col_j(c),k)
          elseif (k>1) then
            zBottomOfCell = -G%bathyT(CS%col_i(c),CS%col_j(c))
            tmpT1d(k) = tmpT1d(k-1)
          else ! This next block should only ever be reached over land
            tmpT1d(k) = -99.9
          endif
          hsrc(k) = zTopOfCell - zBottomOfCell
          if (hsrc(k)>0.) nPoints = nPoints + 1
          zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
        enddo
        ! In case data is deeper than model
        hsrc(nz_data) = hsrc(nz_data) + ( zTopOfCell + G%bathyT(CS%col_i(c),CS%col_j(c)) )
        CS%Ref_val(m)%h(1:nz_data,c) = GV%Z_to_H*hsrc(1:nz_data)
        CS%Ref_val(m)%p(1:nz_data,c) = tmpT1d(1:nz_data)
        do k=2,nz_data
          if (CS%Ref_val(m)%h(k,c) <= 0.001*GV%m_to_H) &
            ! some confusion here about why the masks are not correct returning from horiz_interp
            ! reverting to using a minimum thickness criteria
            CS%Ref_val(m)%p(k,c) = CS%Ref_val(m)%p(k-1,c)
        enddo
      enddo
      deallocate(sp_val, mask_z, hsrc, tmpT1d)
    enddo
  endif

  tmp_val1(:)=0.0;h_col(:)=0.0
  do m=1,CS%fldno
    nz_data = CS%Ref_val(m)%nz_data
    allocate(tmp_val2(CS%Ref_val(m)%nz_data))
    if (CS%id_sp_tendency(m) > 0) then
      allocate(tmp(G%isd:G%ied,G%jsd:G%jed,nz));tmp(:,:,:) = 0.0
    endif
    do c=1,CS%num_col
      ! c is an index for the next 3 lines but a multiplier for the rest of the loop
      ! Therefore we use c as per C code and increment the index where necessary.
      i = CS%col_i(c) ; j = CS%col_j(c)
      damp = dt * CS%Iresttime_col(c)
      I1pdamp = 1.0 / (1.0 + damp)
      tmp_val2(1:nz_data) = CS%Ref_val(m)%p(1:nz_data,c)
      do k=1,nz
        h_col(k)=h(i,j,k)
      enddo
      if (CS%time_varying_sponges) then

        call remapping_core_h(CS%remap_cs, nz_data, CS%Ref_val(m)%h(1:nz_data,c), tmp_val2, &
             CS%nz, h_col, tmp_val1, h_neglect, h_neglect_edge)
      else
        call remapping_core_h(CS%remap_cs, nz_data, CS%Ref_h%p(1:nz_data,c), tmp_val2, &
             CS%nz, h_col, tmp_val1, h_neglect, h_neglect_edge)
      endif
      !Backward Euler method
      if (CS%id_sp_tendency(m) > 0) tmp(i,j,1:nz) = CS%var(m)%p(i,j,1:nz)
      CS%var(m)%p(i,j,1:nz) = I1pdamp * (CS%var(m)%p(i,j,1:nz) + tmp_val1(1:nz) * damp)
      if (CS%id_sp_tendency(m) > 0) &
           tmp(i,j,1:CS%nz) = Idt*(CS%var(m)%p(i,j,1:nz) - tmp(i,j,1:nz))
    enddo

    if (CS%id_sp_tendency(m) > 0)  then
      call post_data(CS%id_sp_tendency(m), tmp, CS%diag)
      deallocate(tmp)
    endif
    deallocate(tmp_val2)
  enddo

  if (CS%sponge_uv) then

    if (CS%time_varying_sponges) then
      if (.not. present(Time)) &
         call MOM_error(FATAL,"apply_ALE_sponge: No time information provided")

      nz_data = CS%Ref_val_u%nz_data
      allocate(sp_val(G%isd:G%ied,G%jsd:G%jed,1:nz_data))
      allocate(sp_val_u(G%isdB:G%iedB,G%jsd:G%jed,1:nz_data))
      allocate(mask_z(G%isdB:G%iedB,G%jsd:G%jed,1:nz_data))
      sp_val(:,:,:) = 0.0
      sp_val_u(:,:,:) = 0.0
      mask_z(:,:,:) = 0.0
      ! Interpolate from the external horizontal grid and in time
      call horiz_interp_and_extrap_tracer(CS%Ref_val_u%id, Time, 1.0, G, sp_val, mask_z, z_in, &
                                          z_edges_in, missing_value, CS%reentrant_x, CS%tripolar_N, .false., &
                                          spongeOnGrid=CS%SpongeDataOngrid, m_to_Z=US%m_to_Z,&
                                          answers_2018=CS%hor_regrid_answers_2018)

      call pass_var(sp_val,G%Domain)
      do j=CS%jsc,CS%jec; do I=CS%iscB,CS%iecB
       sp_val_u(I,j,1:nz_data) = 0.5*(sp_val(i,j,1:nz_data)+sp_val(i+1,j,1:nz_data))
      enddo ; enddo

      allocate( hsrc(nz_data) )
      allocate( tmpT1d(nz_data) )
      do c=1,CS%num_col_u
        ! c is an index for the next 3 lines but a multiplier for the rest of the loop
        ! Therefore we use c as per C code and increment the index where necessary.
        i = CS%col_i_u(c) ; j = CS%col_j_u(c)
        CS%Ref_val_u%p(1:nz_data,c) = sp_val_u(i,j,1:nz_data)
        ! Build the source grid
        zTopOfCell = 0. ; zBottomOfCell = 0. ; nPoints = 0; hsrc(:) = 0.0; tmpT1d(:) = -99.9
        do k=1,nz_data
          if (mask_z(i,j,k) == 1.0) then
            zBottomOfCell = -min( z_edges_in(k+1), G%bathyT(i,j) )
            tmpT1d(k) = sp_val_u(i,j,k)
          elseif (k>1) then
            zBottomOfCell = -G%bathyT(i,j)
            tmpT1d(k) = tmpT1d(k-1)
          else ! This next block should only ever be reached over land
            tmpT1d(k) = -99.9
          endif
          hsrc(k) = zTopOfCell - zBottomOfCell
          if (hsrc(k)>0.) nPoints = nPoints + 1
            zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
        enddo
        ! In case data is deeper than model
        hsrc(nz_data) = hsrc(nz_data) + ( zTopOfCell + G%bathyT(i,j) )
        CS%Ref_val_u%h(1:nz_data,c) = GV%Z_to_H*hsrc(1:nz_data)
      enddo
      deallocate(sp_val, sp_val_u, mask_z, hsrc, tmpT1d)
      nz_data = CS%Ref_val_v%nz_data
      allocate(sp_val( G%isd:G%ied,G%jsd:G%jed,1:nz_data))
      allocate(sp_val_v(G%isd:G%ied,G%jsdB:G%jedB,1:nz_data))
      allocate(mask_z(G%isd:G%ied,G%jsdB:G%jedB,1:nz_data))
      sp_val(:,:,:) = 0.0
      sp_val_v(:,:,:) = 0.0
      mask_z(:,:,:) = 0.0
      ! Interpolate from the external horizontal grid and in time
      call horiz_interp_and_extrap_tracer(CS%Ref_val_v%id, Time, 1.0, G, sp_val, mask_z, z_in, &
                                          z_edges_in, missing_value, CS%reentrant_x, CS%tripolar_N, .false., &
                                          spongeOnGrid=CS%SpongeDataOngrid, m_to_Z=US%m_to_Z,&
                                          answers_2018=CS%hor_regrid_answers_2018)
      call pass_var(sp_val,G%Domain)
      do J=CS%jscB,CS%jecB; do i=CS%isc,CS%iec
        sp_val_v(i,J,1:nz_data) = 0.5*(sp_val(i,j,1:nz_data)+sp_val(i,j+1,1:nz_data))
      enddo ; enddo
      !call pass_var(mask_z,G%Domain)
      allocate( hsrc(nz_data) )
      allocate( tmpT1d(nz_data) )
      do c=1,CS%num_col_v
        ! c is an index for the next 3 lines but a multiplier for the rest of the loop
        ! Therefore we use c as per C code and increment the index where necessary.
        i = CS%col_i_v(c) ; j = CS%col_j_v(c)
        CS%Ref_val_v%p(1:nz_data,c) = sp_val_v(i,j,1:nz_data)
        ! Build the source grid
        zTopOfCell = 0. ; zBottomOfCell = 0. ; nPoints = 0; hsrc(:) = 0.0; tmpT1d(:) = -99.9
        do k=1,nz_data
          if (mask_z(i,j,k) == 1.0) then
            zBottomOfCell = -min( z_edges_in(k+1), G%bathyT(i,j) )
            tmpT1d(k) = sp_val_v(i,j,k)
          elseif (k>1) then
            zBottomOfCell = -G%bathyT(i,j)
            tmpT1d(k) = tmpT1d(k-1)
          else ! This next block should only ever be reached over land
            tmpT1d(k) = -99.9
          endif
          hsrc(k) = zTopOfCell - zBottomOfCell
          if (hsrc(k)>0.) nPoints = nPoints + 1
            zTopOfCell = zBottomOfCell ! Bottom becomes top for next value of k
        enddo
        ! In case data is deeper than model
        hsrc(nz_data) = hsrc(nz_data) + ( zTopOfCell + G%bathyT(i,j) )
        CS%Ref_val_v%h(1:nz_data,c) = GV%Z_to_H*hsrc(1:nz_data)
      enddo
      deallocate(sp_val, sp_val_v, mask_z, hsrc, tmpT1d)
    endif

    nz_data = CS%Ref_val_u%nz_data
    allocate(tmp_val2(nz_data))
    if (CS%id_sp_u_tendency > 0) then
      allocate(tmp_u(G%isdB:G%iedB,G%jsd:G%jed,nz));tmp_u(:,:,:)=0.0
    endif
    ! u points
    do c=1,CS%num_col_u
      I = CS%col_i_u(c) ; j = CS%col_j_u(c)
      damp = dt * CS%Iresttime_col_u(c)
      I1pdamp = 1.0 / (1.0 + damp)
      tmp_val2(1:nz_data) = CS%Ref_val_u%p(1:nz_data,c)
      do k=1,nz
        h_col(k) =  0.5 * (h(i,j,k) + h(i+1,j,k))
      enddo
      if (CS%time_varying_sponges) then
        call remapping_core_h(CS%remap_cs, nz_data, CS%Ref_val_u%h(1:nz_data,c), tmp_val2, &
                 CS%nz, h_col, tmp_val1, h_neglect, h_neglect_edge)
      else
        call remapping_core_h(CS%remap_cs, nz_data, CS%Ref_hu%p(1:nz_data,c), tmp_val2, &
                 CS%nz, h_col, tmp_val1, h_neglect, h_neglect_edge)
      endif
      if (CS%id_sp_u_tendency > 0) tmp_u(i,j,1:nz) = CS%var_u%p(i,j,1:nz)
      !Backward Euler method
      CS%var_u%p(i,j,1:nz) = I1pdamp * (CS%var_u%p(i,j,1:nz) + tmp_val1 * damp)
      if (CS%id_sp_u_tendency > 0) tmp_u(i,j,1:nz) = Idt*(CS%var_u%p(i,j,1:nz) - tmp_u(i,j,1:nz))
    enddo
    deallocate(tmp_val2)
    if (CS%id_sp_u_tendency > 0) then
      call post_data(CS%id_sp_u_tendency, tmp_u, CS%diag)
      deallocate(tmp_u)
    endif
    ! v points
    if (CS%id_sp_v_tendency > 0) then
      allocate(tmp_v(G%isd:G%ied,G%jsdB:G%jedB,nz));tmp_v(:,:,:)=0.0
    endif
    nz_data = CS%Ref_val_v%nz_data
    allocate(tmp_val2(nz_data))

    do c=1,CS%num_col_v
      i = CS%col_i_v(c) ; j = CS%col_j_v(c)
      damp = dt * CS%Iresttime_col_v(c)
      I1pdamp = 1.0 / (1.0 + damp)
      if (CS%time_varying_sponges) nz_data = CS%Ref_val_v%nz_data
      tmp_val2(1:nz_data) = CS%Ref_val_v%p(1:nz_data,c)
      do k=1,nz
        h_col(k) =  0.5 * (h(i,j,k) + h(i,j+1,k))
      enddo
      if (CS%time_varying_sponges) then
        call remapping_core_h(CS%remap_cs, nz_data, CS%Ref_val_v%h(1:nz_data,c), tmp_val2, &
                 CS%nz, h_col, tmp_val1, h_neglect, h_neglect_edge)
      else
        call remapping_core_h(CS%remap_cs, nz_data, CS%Ref_hv%p(1:nz_data,c), tmp_val2, &
                 CS%nz, h_col, tmp_val1, h_neglect, h_neglect_edge)
      endif
      if (CS%id_sp_v_tendency > 0) tmp_v(i,j,1:nz) = CS%var_v%p(i,j,1:nz)
      !Backward Euler method
      CS%var_v%p(i,j,1:nz) = I1pdamp * (CS%var_v%p(i,j,1:nz) + tmp_val1 * damp)
      if (CS%id_sp_v_tendency > 0) tmp_v(i,j,1:nz) = Idt*(CS%var_v%p(i,j,1:nz) - tmp_v(i,j,1:nz))
    enddo
    if (CS%id_sp_v_tendency > 0) then
      call post_data(CS%id_sp_v_tendency, tmp_v, CS%diag)
      deallocate(tmp_v)
    endif
    deallocate(tmp_val2)
  endif




end subroutine apply_ALE_sponge

!> Rotate the ALE sponge fields from the input to the model index map.
subroutine rotate_ALE_sponge(sponge_in, G_in, sponge, G, GV, turns, param_file)
  type(ALE_sponge_CS),     intent(in) :: sponge_in !< The control structure for this module with the
                                                   !! original grid rotation
  type(ocean_grid_type),   intent(in) :: G_in      !< The ocean's grid structure with the original rotation.
  type(ALE_sponge_CS),     pointer    :: sponge    !< A pointer to the control that will be set up with
                                                   !! the new grid rotation
  type(ocean_grid_type),   intent(in) :: G         !< The ocean's grid structure with the new rotation.
  type(verticalGrid_type), intent(in) :: GV        !< The ocean's vertical grid structure
  integer,                 intent(in) :: turns     !< The number of 90-degree turns between grids
  type(param_file_type),   intent(in) :: param_file !< A structure indicating the open file
                                                   !! to parse for model parameter values.

  ! First part: Index construction
  !   1. Reconstruct Iresttime(:,:) from sponge_in
  !   2. rotate Iresttime(:,:)
  !   3. Call initialize_ALE_sponge using new grid and rotated Iresttime(:,:)
  ! All the index adjustment should follow from the Iresttime rotation

  real, dimension(:,:), allocatable :: Iresttime_in, Iresttime
  real, dimension(:,:,:), allocatable :: data_h_in, data_h
  real, dimension(:,:,:), allocatable :: sp_val_in, sp_val
  real, dimension(:,:,:), pointer :: sp_ptr => NULL()
  integer :: c, c_i, c_j
  integer :: k, nz_data
  integer :: n
  logical :: fixed_sponge

  fixed_sponge = .not. sponge_in%time_varying_sponges
  ! NOTE: nz_data is only conditionally set when fixed_sponge is true.

  allocate(Iresttime_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed))
  allocate(Iresttime(G%isd:G%ied, G%jsd:G%jed))
  Iresttime_in(:,:) = 0.0

  if (fixed_sponge) then
    nz_data = sponge_in%nz_data
    allocate(data_h_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed, nz_data))
    allocate(data_h(G%isd:G%ied, G%jsd:G%jed, nz_data))
    data_h_in(:,:,:) = 0.
  endif

  ! Re-populate the 2D Iresttime and data_h arrays on the original grid
  do c=1,sponge_in%num_col
    c_i = sponge_in%col_i(c)
    c_j = sponge_in%col_j(c)
    Iresttime_in(c_i, c_j) = sponge_in%Iresttime_col(c)
    if (fixed_sponge) then
      do k = 1, nz_data
        data_h(c_i, c_j, k) = sponge_in%Ref_h%p(k,c)
      enddo
    endif
  enddo

  call rotate_array(Iresttime_in, turns, Iresttime)
  if (fixed_sponge) then
    call rotate_array(data_h_in, turns, data_h)
    call initialize_ALE_sponge_fixed(Iresttime, G, GV, param_file, sponge, &
                                     data_h, nz_data)
  else
    call initialize_ALE_sponge_varying(Iresttime, G, GV, param_file, sponge)
  endif

  deallocate(Iresttime_in)
  deallocate(Iresttime)
  if (fixed_sponge) then
    deallocate(data_h_in)
    deallocate(data_h)
  endif

  ! Second part: Provide rotated fields for which relaxation is applied

  sponge%fldno = sponge_in%fldno

  if (fixed_sponge) then
    allocate(sp_val_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed, nz_data))
    allocate(sp_val(G%isd:G%ied, G%jsd:G%jed, nz_data))
  endif

  do n=1,sponge_in%fldno
    ! Assume that tracers are pointers and are remapped in other functions(?)
    sp_ptr => sponge_in%var(n)%p
    if (fixed_sponge) then
      sp_val_in(:,:,:) = 0.0
      do c = 1, sponge_in%num_col
        c_i = sponge_in%col_i(c)
        c_j = sponge_in%col_j(c)
        do k = 1, nz_data
          sp_val_in(c_i, c_j, k) = sponge_in%Ref_val(n)%p(k,c)
        enddo
      enddo

      call rotate_array(sp_val_in, turns, sp_val)

      ! NOTE: This points sp_val with the unrotated field.  See note below.
      call set_up_ALE_sponge_field(sp_val, G, GV, sp_ptr, sponge)

      deallocate(sp_val_in)
    else
      ! We don't want to repeat FMS init in set_up_ALE_sponge_field_varying()
      ! (time_interp_external_init, init_external_field, etc), so we manually
      ! do a portion of this function below.
      sponge%Ref_val(n)%id = sponge_in%Ref_val(n)%id
      sponge%Ref_val(n)%num_tlevs = sponge_in%Ref_val(n)%num_tlevs

      nz_data = sponge_in%Ref_val(n)%nz_data
      sponge%Ref_val(n)%nz_data = nz_data

      allocate(sponge%Ref_val(n)%p(nz_data, sponge_in%num_col))
      allocate(sponge%Ref_val(n)%h(nz_data, sponge_in%num_col))
      sponge%Ref_val(n)%p(:,:) = 0.0
      sponge%Ref_val(n)%h(:,:) = 0.0

      ! TODO: There is currently no way to associate a generic field pointer to
      !   its rotated equivalent without introducing a new data structure which
      !   explicitly tracks the pairing.
      !
      !   As a temporary fix, we store the pointer to the unrotated field in
      !   the rotated sponge, and use this reference to replace the pointer
      !   to the rotated field update_ALE_sponge field.
      !
      !   This makes a lot of unverifiable assumptions, and should not be
      !   considered the final solution.
      sponge%var(n)%p => sp_ptr
    endif
  enddo

  ! TODO: var_u and var_v sponge dampling is not yet supported.
  if (associated(sponge_in%var_u%p) .or. associated(sponge_in%var_v%p)) &
    call MOM_error(FATAL, "Rotation of ALE sponge velocities is not yet " &
      // "implemented.")

  ! Transfer any existing diag_CS reference pointer
  sponge%diag => sponge_in%diag

  ! NOTE: initialize_ALE_sponge_* resolves remap_cs
end subroutine rotate_ALE_sponge


!> Scan the ALE sponge variables and replace a prescribed pointer to a new value.
! TODO: This function solely exists to replace field pointers in the sponge
!   after rotation.  This function is part of a temporary solution until
!   something more robust is developed.
subroutine update_ALE_sponge_field(sponge, p_old, G, GV, p_new)
  type(ALE_sponge_CS),     pointer    :: sponge !< A pointer to the control structure for this module
                                               !! that is set by a previous call to initialize_ALE_sponge.
  real, dimension(:,:,:), &
                   target, intent(in) :: p_old !< The previous array of target values
  type(ocean_grid_type),   intent(in) :: G     !< The updated ocean grid structure
  type(verticalGrid_type), intent(in) :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                   target, intent(in) :: p_new !< The new array of target values

  integer :: n

  do n=1,sponge%fldno
    if (associated(sponge%var(n)%p, p_old)) sponge%var(n)%p => p_new
  enddo

end subroutine update_ALE_sponge_field


! GMM: I could not find where sponge_end is being called, but I am keeping
!  ALE_sponge_end here so we can add that if needed.
!> This subroutine deallocates any memory associated with the ALE_sponge module.
subroutine ALE_sponge_end(CS)
  type(ALE_sponge_CS), pointer :: CS !< A pointer to the control structure that is
                                     !! set by a previous call to initialize_ALE_sponge.

  integer :: m

  if (.not.associated(CS)) return

  if (associated(CS%col_i)) deallocate(CS%col_i)
  if (associated(CS%col_i_u)) deallocate(CS%col_i_u)
  if (associated(CS%col_i_v)) deallocate(CS%col_i_v)
  if (associated(CS%col_j)) deallocate(CS%col_j)
  if (associated(CS%col_j_u)) deallocate(CS%col_j_u)
  if (associated(CS%col_j_v)) deallocate(CS%col_j_v)

  if (associated(CS%Iresttime_col)) deallocate(CS%Iresttime_col)
  if (associated(CS%Iresttime_col_u)) deallocate(CS%Iresttime_col_u)
  if (associated(CS%Iresttime_col_v)) deallocate(CS%Iresttime_col_v)

  do m=1,CS%fldno
    if (associated(CS%Ref_val(m)%p)) deallocate(CS%Ref_val(m)%p)
  enddo

  deallocate(CS)

end subroutine ALE_sponge_end

end module MOM_ALE_sponge
