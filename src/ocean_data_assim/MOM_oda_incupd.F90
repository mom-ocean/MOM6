!> This module contains the routines used to apply incremental updates
!! from data assimilation.
!!
!! Applying incremental updates requires the following:
!! 1. initialize_oda_incupd_fixed and initialize_oda_incupd
!! 2. set_up_oda_incupd_field (tracers) and set_up_oda_incupd_vel_field (vel)
!! 3. calc_oda_increments (if using full fields input)
!! 4. apply_oda_incupd
!! 5. output_oda_incupd_inc (output increment if using full fields input)
!! 6. init_oda_incupd_diags (to output increments in diagnostics)
!! 7. oda_incupd_end (not being used for now)

module MOM_oda_incupd


! This file is part of MOM6. See LICENSE.md for the license.
use MOM_array_transform, only : rotate_array
use MOM_coms,            only : sum_across_PEs
use MOM_diag_mediator,   only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator,   only : diag_ctrl
use MOM_domains,         only : pass_var,pass_vector
use MOM_error_handler,   only : MOM_error, FATAL, NOTE, WARNING, is_root_pe
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_get_input,       only : directories, Get_MOM_input
use MOM_grid,            only : ocean_grid_type
use MOM_io,              only : vardesc, var_desc
use MOM_remapping,       only : remapping_cs, remapping_core_h, initialize_remapping
use MOM_restart,         only : register_restart_field, register_restart_pair, MOM_restart_CS
use MOM_restart,         only : restart_init, save_restart, query_initialized
use MOM_spatial_means,   only : global_i_mean
use MOM_time_manager,    only : time_type
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_verticalGrid,    only : get_thickness_units

use mpp_io_mod, only : mpp_get_axis_length
use mpp_io_mod, only : axistype

implicit none ; private

#include <MOM_memory.h>


!  Publicly available functions
public set_up_oda_incupd_field, set_up_oda_incupd_vel_field
public initialize_oda_incupd_fixed, initialize_oda_incupd, apply_oda_incupd, oda_incupd_end
public init_oda_incupd_diags,calc_oda_increments,output_oda_incupd_inc

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A structure for creating arrays of pointers to 3D arrays with extra gridding information
type :: p3d
  integer :: id !< id for FMS external time interpolator
  integer :: nz_data !< The number of vertical levels in the input field.
  real, dimension(:,:,:), pointer :: mask_in => NULL() !< pointer to the data mask.
  real, dimension(:,:,:), pointer :: p => NULL() !< pointer to the data.
  real, dimension(:,:,:), pointer :: h => NULL() !< pointer to the data grid.
end type p3d

!> oda incupd control structure
type, public :: oda_incupd_CS ; private
  integer :: nz        !< The total number of layers.
  integer :: nz_data   !< The total number of arbritary layers (used by older code).
  integer :: fldno = 0 !< The number of fields which have already been
                       !! registered by calls to set_up_oda_incupd_field

  type(p3d) :: Inc(MAX_FIELDS_)      !< The increments to be applied to the field
  type(p3d) :: Inc_u  !< The increments to be applied to the u-velocities
  type(p3d) :: Inc_v  !< The increments to be applied to the v-velocities
  type(p3d) :: Ref_h  !< Vertical grid on which the increments are provided


  integer :: nstep_incupd          !< number of time step for full update
  real    :: ncount = 0.0            !< increment time step counter
  type(remapping_cs) :: remap_cs   !< Remapping parameters and work arrays
  logical :: incupdDataOngrid  !< True if the incupd data are on the model horizontal grid
  logical :: uv_inc      !< use u and v increments

  ! for diagnostics
  type(diag_ctrl), pointer           :: diag !<structure to regulate output
  ! diagnostic for inc. fields
  integer :: id_u_oda_inc = -1 !< diagnostic id for zonal velocity inc.
  integer :: id_v_oda_inc = -1 !< diagnostic id for meridional velocity inc.
  integer :: id_h_oda_inc = -1 !< diagnostic id for layer thicknesses inc.
  integer :: id_T_oda_inc = -1 !< diagnostic id for temperature inc.
  integer :: id_S_oda_inc = -1 !< diagnostic id for salinity inc.

end type oda_incupd_CS

contains

!> This subroutine defined the control structure of module and register
!the time counter to full update in restart
subroutine initialize_oda_incupd_fixed( G, GV, US, CS, restart_CS)

  type(ocean_grid_type),            intent(in) :: G          !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in) :: GV         !< ocean vertical grid structure
  type(unit_scale_type),            intent(in) :: US         !< A dimensional unit scaling type
  type(oda_incupd_CS),              pointer    :: CS         !< A pointer that is set to point to the control
                                                             !! structure for this module (in/out).
  type(MOM_restart_CS),             pointer    :: restart_CS !< A pointer to the restart control structure.


! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=256) :: mesg
  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_oda_incupd_fixed called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  ! initialize time counter
  CS%ncount = 0.0
  ! register ncount in restart
  call register_restart_field(CS%ncount, "oda_incupd_ncount", .false., restart_CS,&
                              "Number of inc. update already done", "N/A")


end subroutine initialize_oda_incupd_fixed


!> This subroutine defined the number of time step for full update, stores the layer pressure
!! increments and initialize remap structure.
subroutine initialize_oda_incupd( G, GV, US, param_file, CS, data_h,nz_data, restart_CS)

  type(ocean_grid_type),            intent(in) :: G          !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in) :: GV         !< ocean vertical grid structure
  type(unit_scale_type),            intent(in) :: US         !< A dimensional unit scaling type
  integer,                          intent(in) :: nz_data    !< The total number of incr. input layers.
  type(param_file_type),            intent(in) :: param_file !< A structure indicating the open file
                                                             !! to parse for model parameter values.
  type(oda_incupd_CS),              pointer    :: CS         !< A pointer that is set to point to the control
                                                             !! structure for this module (in/out).
  real, dimension(SZI_(G),SZJ_(G),nz_data), intent(in) :: data_h !< The ODA h
                                                                 !! [H ~> m or kg m-2].
  type(MOM_restart_CS),       pointer       :: restart_CS !< A pointer to the restart control
                                                    !! structure.


! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_oda"  ! This module's name.
  logical :: use_oda_incupd
  logical :: bndExtrapolation = .true.   ! If true, extrapolate boundaries
  logical :: reset_ncount
  integer :: i, j, k
  real    :: nhours_incupd, dt, dt_therm
  type(vardesc) :: vd
  character(len=256) :: mesg
  character(len=10)  :: remapScheme
  if (.not.associated(CS)) then
    call MOM_error(WARNING, "initialize_oda_incupd called without an associated "// &
                            "control structure.")
    return
  endif

! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ODA_INCUPD", use_oda_incupd, &
                 "If true, oda incremental updates will be applied "//&
                 "everywhere in the domain.", default=.false.)

  if (.not.use_oda_incupd) return

  call get_param(param_file, mdl, "ODA_INCUPD_NHOURS", nhours_incupd, &
                 "Number of hours for full update (0=direct insertion).", &
                 default=3.0,units="h", scale=US%s_to_T)
  call get_param(param_file, mdl, "ODA_INCUPD_RESET_NCOUNT", reset_ncount, &
                 "If True, reinitialize number of updates already done, ncount.", &
                 default=.true.)
  call get_param(param_file, mdl, "DT", dt, &
                 "The (baroclinic) dynamics time step.  The time-step that "//&
                 "is actually used will be an integer fraction of the "//&
                 "forcing time-step (DT_FORCING in ocean-only mode or the "//&
                 "coupling timestep in coupled mode.)", units="s", scale=US%s_to_T, &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "DT_THERM", dt_therm, &
                 "The thermodynamic and tracer advection time step. "//&
                 "Ideally DT_THERM should be an integer multiple of DT "//&
                 "and less than the forcing or coupling time-step, unless "//&
                 "THERMO_SPANS_COUPLING is true, in which case DT_THERM "//&
                 "can be an integer multiple of the coupling timestep.  By "//&
                 "default DT_THERM is set to DT.", &
                 units="s", scale=US%s_to_T, default=US%T_to_s*dt)
  call get_param(param_file, mdl, "ODA_INCUPD_UV", CS%uv_inc, &
                 "use U,V increments.", &
                 default=.true.)
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
  call get_param(param_file, mdl, "ODA_INCUPD_DATA_ONGRID", CS%incupdDataOngrid, &
                 "When defined, the incoming oda_incupd data are "//&
                 "assumed to be on the model horizontal grid " , &
                 default=.true.)

  CS%nz = GV%ke

  ! increments on horizontal grid
  if (.not. CS%incupdDataOngrid) call MOM_error(FATAL,'initialize_oda_incupd: '// &
           'The oda_incupd code only applies ODA increments on the same horizontal grid. ')

  ! get number of timestep for full update
  if (nhours_incupd == 0) then
     CS%nstep_incupd = 1 !! direct insertion
  else
     CS%nstep_incupd = floor( nhours_incupd * 3600. / dt_therm + 0.001 ) - 1
  endif
  write(mesg,'(i12)') CS%nstep_incupd
  if (is_root_pe()) &
     call MOM_error(NOTE,"initialize_oda_incupd: Number of Timestep of inc. update:"//&
                       trim(mesg))

  ! number of inc. update already done, CS%ncount, either from restart or set to 0.0
  if (query_initialized(CS%ncount, "oda_incupd_ncount", restart_CS) .and. &
      .not.reset_ncount) then
    CS%ncount = CS%ncount
  else
    CS%ncount = 0.0
  endif
  write(mesg,'(f4.1)') CS%ncount
  if (is_root_pe()) &
     call MOM_error(NOTE,"initialize_oda_incupd: Inc. update already done:"//&
                       trim(mesg))

  ! get the vertical grid (h_obs) of the increments
  CS%nz_data = nz_data
  allocate(CS%Ref_h%p(G%isd:G%ied,G%jsd:G%jed,CS%nz_data))
  CS%Ref_h%p(:,:,:) = 0.0 ;
  do j=G%jsc,G%jec; do i=G%isc,G%iec ; do k=1,CS%nz_data
      CS%Ref_h%p(i,j,k) = data_h(i,j,k)
  enddo;  enddo ; enddo

  ! Call the constructor for remapping control structure
  call initialize_remapping(CS%remap_cs, remapScheme, boundary_extrapolation=bndExtrapolation, &
                            answers_2018=.false.)


end subroutine initialize_oda_incupd


!> This subroutine stores the increments at h points for the variable
!! whose address is given by f_ptr.
subroutine set_up_oda_incupd_field(sp_val, G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G      !< Grid structure
  type(verticalGrid_type), intent(in) :: GV     !< ocean vertical grid structure
  type(oda_incupd_CS),     pointer    :: CS     !< oda_incupd control structure (in/out).
  real, dimension(SZI_(G),SZJ_(G),CS%nz_data), &
                           intent(in) :: sp_val !< increment field, it can have an
                                                !! arbitrary number of layers.

  integer :: i, j, k
  character(len=256) :: mesg ! String for error messages

  if (.not.associated(CS)) return

  CS%fldno = CS%fldno + 1
  if (CS%fldno > MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ to at least ",I3," in MOM_memory.h or decrease &
           &the number of fields increments in the call to &
           &initialize_oda_incupd." )') CS%fldno
    call MOM_error(FATAL,"set_up_oda_incupd_field: "//mesg)
  endif

  ! store the increment/full field tracer profiles
  CS%Inc(CS%fldno)%nz_data = CS%nz_data
  allocate(CS%Inc(CS%fldno)%p(G%isd:G%ied,G%jsd:G%jed,CS%nz_data))
  CS%Inc(CS%fldno)%p(:,:,:) = 0.0
  do k=1,CS%nz_data ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
     CS%Inc(CS%fldno)%p(i,j,k) = sp_val(i,j,k)
  enddo ; enddo ; enddo

end subroutine set_up_oda_incupd_field


!> This subroutine stores the increments at u and v points for the variable
!! whose address is given by u_ptr and v_ptr.
subroutine set_up_oda_incupd_vel_field(u_val, v_val, G, GV, CS)
  type(ocean_grid_type),   intent(in) :: G  !< Grid structure (in).
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  type(oda_incupd_CS),     pointer    :: CS !< oda incupd structure (in/out).

  real, dimension(SZIB_(G),SZJ_(G),CS%nz_data), &
                          intent(in) :: u_val !< u increment, it has arbritary number of layers but
                                              !! not to exceed the total number of model layers [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),CS%nz_data), &
                          intent(in) :: v_val !< v increment, it has arbritary number of layers but
                                              !! not to exceed the number of model layers [L T-1 ~> m s-1]
  integer :: i, j, k

  if (.not.associated(CS)) return


  ! store the increment/full field u profile
  allocate(CS%Inc_u%p(G%isdB:G%iedB,G%jsd:G%jed,CS%nz_data))
  CS%Inc_u%p(:,:,:) = 0.0
  do j=G%jsc,G%jec ; do i=G%iscB,G%iecB
    do k=1,CS%nz_data
      CS%Inc_u%p(i,j,k) = u_val(i,j,k)
    enddo
  enddo ; enddo

  ! store the increment/full field v profile
  allocate(CS%Inc_v%p(G%isd:G%ied,G%jsdB:G%jedB,CS%nz_data))
  CS%Inc_v%p(:,:,:) = 0.0
  do j=G%jscB,G%jecB ; do i=G%isc,G%iec
    do k=1,CS%nz_data
      CS%Inc_v%p(i,j,k) = v_val(i,j,k)
    enddo
  enddo ; enddo

end subroutine set_up_oda_incupd_vel_field

! calculation of the increments if using full fields (ODA_INCUPD_INC=.false.) at initialization
subroutine calc_oda_increments(h, tv, u, v, G, GV, US, CS)

  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure (in).
  type(verticalGrid_type),   intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2] (in)
  type(thermo_var_ptrs),     intent(in)    :: tv !< A structure pointing to various thermodynamic variables

  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)      :: u    !< The zonal velocity that is being
                                                   !! initialized [L T-1 ~> m s-1]
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)      :: v    !< The meridional velocity that is being
                                                   !! initialized [L T-1 ~> m s-1]
  type(oda_incupd_CS),       pointer       :: CS !< A pointer to the control structure for this module
                                                 !! that is set by a previous call to initialize_oda_incupd (in).


  real, dimension(SZK_(GV)) :: tmp_val1          ! data values on the model grid
  real, allocatable, dimension(:) :: tmp_val2    ! data values remapped to increment grid
  real, allocatable, dimension(:,:,:) :: h_obs   !< h of increments
  real, allocatable, dimension(:) :: tmp_h       ! temporary array for corrected h_obs
  real, allocatable, dimension(:) :: hu_obs,hv_obs  ! A column of thicknesses at h points [H ~> m or kg m-2]
  real, dimension(SZK_(GV)) :: hu, hv            ! A column of thicknesses at h, u or v points [H ~> m or kg m-2]


  integer ::  i, j, k, is, ie, js, je, nz, nz_data
  integer :: isB, ieB, jsB, jeB
  real :: h_neglect, h_neglect_edge
  real :: sum_h1, sum_h2 !vertical sums of h's
  character(len=256) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isB = G%iscB ; ieB = G%iecB ; jsB = G%jscB ; jeB = G%jecB
  if (.not.associated(CS)) return


  ! increments calculated on if CS%ncount = 0.0
  if (CS%ncount /= 0.0) call MOM_error(FATAL,'calc_oda_increments: '// &
           'CS%ncount should be 0.0 to get accurate increments.')


  if (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  endif

  ! get h_obs
  nz_data = CS%Inc(1)%nz_data
  allocate(h_obs(G%isd:G%ied,G%jsd:G%jed,nz_data)) ; h_obs(:,:,:) = 0.0
  do k=1,nz_data  ; do j=js,je ; do i=is,ie
    h_obs(i,j,k) = CS%Ref_h%p(i,j,k)
  enddo ; enddo ; enddo
  call pass_var(h_obs,G%Domain)


  ! allocate 1-d arrays
  allocate(tmp_h(nz_data)); tmp_h(:) = 0.0
  allocate(tmp_val2(nz_data)) ; tmp_val2(:) = 0.0
  allocate(hu_obs(nz_data)) ; hu_obs(:) = 0.0
  allocate(hv_obs(nz_data)) ; hv_obs(:) = 0.0

  ! remap t,s (on h_init) to h_obs to get increment
  tmp_val1(:) = 0.0
  do j=js,je ; do i=is,ie
     if (G%mask2dT(i,j) == 1) then
        ! account for the different SSH
        sum_h1 = 0.0
        sum_h2 = 0.0
        do k=1,nz
          sum_h1 = sum_h1+h(i,j,k)
        enddo

        do k=1,nz_data
          sum_h2 = sum_h2+h_obs(i,j,k)
        enddo
        do k=1,nz_data
          tmp_h(k)=(sum_h1/sum_h2)*h_obs(i,j,k)
        enddo
        ! get temperature
        do k=1,nz
           tmp_val1(k) = tv%T(i,j,k)
        enddo
        ! remap tracer on h_obs
        call remapping_core_h(CS%remap_cs, nz, h(i,j,1:nz), tmp_val1, &
                              nz_data, tmp_h(1:nz_data), tmp_val2, &
                              h_neglect, h_neglect_edge)
        ! get increment from full field on h_obs
        do k=1,nz_data
           CS%Inc(1)%p(i,j,k) = CS%Inc(1)%p(i,j,k) - tmp_val2(k)
        enddo

        ! get salinity
        do k=1,nz
           tmp_val1(k) = tv%S(i,j,k)
        enddo
        ! remap tracer on h_obs
        call remapping_core_h(CS%remap_cs, nz, h(i,j,1:nz), tmp_val1, &
                              nz_data, tmp_h(1:nz_data), tmp_val2, &
                              h_neglect, h_neglect_edge)
        ! get increment from full field on h_obs
        do k=1,nz_data
           CS%Inc(2)%p(i,j,k) = CS%Inc(2)%p(i,j,k) - tmp_val2(k)
        enddo
     endif
  enddo; enddo

  ! remap u to h_obs to get increment
  if (CS%uv_inc) then
     call pass_var(h, G%Domain)

     hu(:) = 0.0
     do j=js,je ; do i=isB,ieB
        if (G%mask2dCu(i,j) == 1) then
           ! get u-velocity
           do k=1,nz
              tmp_val1(k) = u(i,j,k)
              ! get the h and h_obs at u points
              hu(k) = 0.5*( h(i,j,k)+ h(i+1,j,k))
           enddo
           do k=1,nz_data
              hu_obs(k) = 0.5*(h_obs(i,j,k)+h_obs(i+1,j,k))
           enddo
           ! account for the different SSH
           sum_h1 = 0.0
           do k=1,nz
             sum_h1 = sum_h1+hu(k)
           enddo
           sum_h2 = 0.0
           do k=1,nz_data
             sum_h2 = sum_h2+hu_obs(k)
           enddo
           do k=1,nz_data
             hu_obs(k)=(sum_h1/sum_h2)*hu_obs(k)
           enddo
           ! remap model u on hu_obs
           call remapping_core_h(CS%remap_cs, nz, hu(1:nz), tmp_val1, &
                                 nz_data, hu_obs(1:nz_data), tmp_val2, &
                                 h_neglect, h_neglect_edge)
           ! get increment from full field on h_obs
           do k=1,nz_data
              CS%Inc_u%p(i,j,k) = CS%Inc_u%p(i,j,k) - tmp_val2(k)
           enddo
        endif
     enddo; enddo

     ! remap v to h_obs to get increment
     hv(:) = 0.0;
     do j=jsB,jeB ; do i=is,ie
        if (G%mask2dCv(i,j) == 1) then
           ! get v-velocity
           do k=1,nz
              tmp_val1(k) = v(i,j,k)
              ! get the h and h_obs at v points
              hv(k) = 0.5*(h(i,j,k)+h(i,j+1,k))
           enddo
           do k=1,nz_data
              hv_obs(k) = 0.5*(h_obs(i,j,k)+h_obs(i,j+1,k))
           enddo
           ! account for the different SSH
           sum_h1 = 0.0
           do k=1,nz
             sum_h1 = sum_h1+hv(k)
           enddo
           sum_h2 = 0.0
           do k=1,nz_data
             sum_h2 = sum_h2+hv_obs(k)
           enddo
           do k=1,nz_data
             hv_obs(k)=(sum_h1/sum_h2)*hv_obs(k)
           enddo
           ! remap model v on hv_obs
           call remapping_core_h(CS%remap_cs, nz, hv(1:nz), tmp_val1, &
                                 nz_data, hv_obs(1:nz_data), tmp_val2, &
                                 h_neglect, h_neglect_edge)
           ! get increment from full field on h_obs
           do k=1,nz_data
              CS%Inc_v%p(i,j,k) = CS%Inc_v%p(i,j,k) - tmp_val2(k)
           enddo
        endif
     enddo; enddo
  endif ! uv_inc

  call pass_var(CS%Inc(1)%p, G%Domain)
  call pass_var(CS%Inc(2)%p, G%Domain)
  call pass_vector(CS%Inc_u%p,CS%Inc_v%p,G%Domain)

  ! deallocate  arrays
  deallocate(tmp_h,tmp_val2,hu_obs,hv_obs)
  deallocate(h_obs)

end subroutine calc_oda_increments

!> This subroutine applies oda increments to layers thicknesses, temp, salt, U
!!  and V everywhere .
subroutine apply_oda_incupd(h, tv, u, v, dt, G, GV, US, CS)
  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure (in).
  type(verticalGrid_type),   intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(inout)    :: h  !< Layer thickness [H ~> m or kg m-2] (in)
  type(thermo_var_ptrs),     intent(inout) :: tv !< A structure pointing to various thermodynamic variables

  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout)   :: u    !< The zonal velocity that is being
                                                   !! initialized [L T-1 ~> m s-1]
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout)   :: v    !< The meridional velocity that is being
                                                   !! initialized [L T-1 ~> m s-1]

  real,                      intent(in)    :: dt !< The amount of time covered by this call [T ~> s].
  type(oda_incupd_CS),       pointer       :: CS !< A pointer to the control structure for this module
                                                 !! that is set by a previous call to initialize_oda_incupd (in).

  real :: m_to_Z                                ! A unit conversion factor from m to Z.
  real, allocatable, dimension(:) :: tmp_val2   ! data values on the increment grid
  real, dimension(SZK_(GV)) :: tmp_val1         ! data values remapped to model grid
  real, dimension(SZK_(GV)) :: hu, hv           ! A column of thicknesses at h, u or v points [H ~> m or kg m-2]

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: tmp_t  !< A temporary array for t inc.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: tmp_s  !< A temporary array for s inc.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: tmp_u  !< A temporary array for u inc.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: tmp_v  !< A temporary array for v inc.

  real, allocatable, dimension(:,:,:) :: h_obs     !< h of increments
  real, allocatable, dimension(:) :: tmp_h         !< temporary array for corrected h_obs
  real, allocatable, dimension(:) :: hu_obs,hv_obs  ! A column of thicknesses at h points [H ~> m or kg m-2]

  integer ::  i, j, k, is, ie, js, je, nz, nz_data
  integer :: isB, ieB, jsB, jeB
!  integer :: ncount      ! time step counter
  real :: inc_wt           ! weight of the update for this time-step
  real :: h_neglect, h_neglect_edge
  real :: sum_h1, sum_h2 !vertical sums of h's
  character(len=256) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isB = G%iscB ; ieB = G%iecB ; jsB = G%jscB ; jeB = G%jecB
  if (.not.associated(CS)) return

  ! no assimilation after CS%step_incupd
  if (CS%ncount >= CS%nstep_incupd) then
    if (is_root_pe()) call MOM_error(NOTE,"ended updating fields with increments. ")
    return
  endif !ncount>CS%nstep_incupd

  ! update counter
  CS%ncount = CS%ncount+1.0
  inc_wt = 1.0/CS%nstep_incupd

  ! print out increments
  write(mesg,'(f10.0)') CS%ncount
  if (is_root_pe()) call MOM_error(NOTE,"updating fields with increments ncount:"//trim(mesg))
  write(mesg,'(f10.8)') inc_wt
  if (is_root_pe()) call MOM_error(NOTE,"updating fields with weight inc_wt:"//trim(mesg))

  if (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  endif

  ! get h_obs
  nz_data = CS%Inc(1)%nz_data
  allocate(h_obs(G%isd:G%ied,G%jsd:G%jed,nz_data)) ; h_obs(:,:,:) = 0.0
  do k=1,nz_data  ; do j=js,je ; do i=is,ie
    h_obs(i,j,k) = CS%Ref_h%p(i,j,k)
  enddo ; enddo ; enddo
  call pass_var(h_obs,G%Domain)

  ! allocate 1-d array
  allocate(tmp_h(nz_data)); tmp_h(:) = 0.0
  allocate(tmp_val2(nz_data))
  allocate(hu_obs(nz_data)) ; hu_obs(:) = 0.0
  allocate(hv_obs(nz_data)) ; hv_obs(:) = 0.0

  ! add increments to tracers
  tmp_val1(:) = 0.0
  tmp_t(:,:,:) = 0.0 ; tmp_s(:,:,:) = 0.0 ! diagnostics
  do j=js,je ; do i=is,ie
    ! account for the different SSH
    sum_h1 = 0.0
    do k=1,nz
      sum_h1 = sum_h1+h(i,j,k)
    enddo
    sum_h2 = 0.0
    do k=1,nz_data
      sum_h2 = sum_h2+h_obs(i,j,k)
    enddo
    do k=1,nz_data
      tmp_h(k) = ( sum_h1 / sum_h2 ) * h_obs(i,j,k)
    enddo
    if (G%mask2dT(i,j) == 1) then
    ! get temperature increment
    do k=1,nz_data
       tmp_val2(k) = CS%Inc(1)%p(i,j,k)
    enddo
    ! remap increment profile on model h
    call remapping_core_h(CS%remap_cs, nz_data, tmp_h(1:nz_data), tmp_val2, &
                          nz, h(i,j,1:nz),tmp_val1, h_neglect, h_neglect_edge)
    do k=1,nz
    ! add increment to tracer on model h
       tv%T(i,j,k) = tv%T(i,j,k) + inc_wt * tmp_val1(k)
       tmp_t(i,j,k) = tmp_val1(k) ! store T increment for diagnostics
    enddo

    ! get salinity increment
    do k=1,nz_data
       tmp_val2(k) = CS%Inc(2)%p(i,j,k)
    enddo
    ! remap increment profile on model h
    call remapping_core_h(CS%remap_cs, nz_data, tmp_h(1:nz_data),tmp_val2,&
                          nz, h(i,j,1:nz),tmp_val1, h_neglect, h_neglect_edge)
    ! add increment to tracer on model h
    do k=1,nz
       tv%S(i,j,k) = tv%S(i,j,k) + inc_wt * tmp_val1(k)
       tmp_s(i,j,k) = tmp_val1(k) ! store S increment for diagnostics
    ! bound salinity values ! check if it is correct to do that or if it hides
    ! other problems ...
      tv%S(i,j,k) = max(0.0 , tv%S(i,j,k))
    enddo
    endif
  enddo; enddo


  ! add u and v increments
  if (CS%uv_inc) then

     call pass_var(h,G%Domain) ! to ensure reproducibility

     ! add increments to u
     hu(:) = 0.0
     tmp_u(:,:,:) = 0.0 ! diagnostics
     do j=js,je ; do i=isB,ieB
        if (G%mask2dCu(i,j) == 1) then
           do k=1,nz_data
              ! get u increment
              tmp_val2(k) = CS%Inc_u%p(i,j,k)
              ! get the h and h_obs at u points
              hu_obs(k) = 0.5 * ( h_obs(i,j,k) + h_obs(i+1,j,k) )
           enddo
           do k=1,nz
              hu(k) = 0.5 * ( h(i,j,k) + h(i+1,j,k) )
           enddo
           ! account for different SSH
           sum_h1 = 0.0
           do k=1,nz
             sum_h1 = sum_h1 + hu(k)
           enddo
           sum_h2 = 0.0
           do k=1,nz_data
             sum_h2 = sum_h2 + hu_obs(k)
           enddo
           do k=1,nz_data
             hu_obs(k)=( sum_h1 / sum_h2 ) * hu_obs(k)
           enddo
           ! remap increment profile on hu
           call remapping_core_h(CS%remap_cs, nz_data, hu_obs(1:nz_data), tmp_val2, &
                                 nz, hu(1:nz), tmp_val1, h_neglect, h_neglect_edge)
           ! add increment to u-velocity on hu
           do k=1,nz
              u(i,j,k) = u(i,j,k) + inc_wt * tmp_val1(k)
              ! store increment for diagnostics
              tmp_u(i,j,k) = tmp_val1(k)
           enddo
        endif
     enddo; enddo

     ! add increments to v
     hv(:) = 0.0
     tmp_v(:,:,:) = 0.0 ! diagnostics
     do j=jsB,jeB ; do i=is,ie
        if (G%mask2dCv(i,j) == 1) then
           ! get v increment
           do k=1,nz_data
              tmp_val2(k) = CS%Inc_v%p(i,j,k)
              ! get the h and h_obs at v points
              hv_obs(k) = 0.5 * ( h_obs(i,j,k) + h_obs(i,j+1,k) )
           enddo
           do k=1,nz
              hv(k) = 0.5 * (h(i,j,k) + h(i,j+1,k) )
           enddo
           ! account for different SSH
           sum_h1 = 0.0
           do k=1,nz
             sum_h1 = sum_h1 + hv(k)
           enddo
           sum_h2 = 0.0
           do k=1,nz_data
             sum_h2 = sum_h2 + hv_obs(k)
           enddo
           do k=1,nz_data
             hv_obs(k)=( sum_h1 / sum_h2 ) * hv_obs(k)
           enddo
           ! remap increment profile on hv
           call remapping_core_h(CS%remap_cs, nz_data, hv_obs(1:nz_data), tmp_val2, &
                                 nz, hv(1:nz), tmp_val1, h_neglect, h_neglect_edge)
           ! add increment to v-velocity on hv
           do k=1,nz
              v(i,j,k) = v(i,j,k) + inc_wt * tmp_val1(k)
              ! store increment for diagnostics
              tmp_v(i,j,k) = tmp_val1(k)
           enddo
        endif
     enddo; enddo

  endif ! uv_inc

  call pass_var(tv%T, G%Domain)
  call pass_var(tv%S, G%Domain)
  call pass_vector(u,v,G%Domain)

  ! Diagnostics of increments, mostly used for debugging.
  if (CS%uv_inc) then
    if (CS%id_u_oda_inc > 0) call post_data(CS%id_u_oda_inc, tmp_u, CS%diag)
    if (CS%id_v_oda_inc > 0) call post_data(CS%id_v_oda_inc, tmp_v, CS%diag)
  endif
  if (CS%id_h_oda_inc > 0) call post_data(CS%id_h_oda_inc, h    , CS%diag)
  if (CS%id_T_oda_inc > 0) call post_data(CS%id_T_oda_inc, tmp_t, CS%diag)
  if (CS%id_S_oda_inc > 0) call post_data(CS%id_S_oda_inc, tmp_s, CS%diag)

  ! deallocate arrays
  deallocate(tmp_h,tmp_val2,hu_obs,hv_obs)
  deallocate(h_obs)

  end subroutine apply_oda_incupd

!> Output increment if using full fields for the oda_incupd module.
subroutine output_oda_incupd_inc(Time, G, GV, param_file, CS, US)
  type(time_type), target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV !< ocean vertical grid structure
  type(param_file_type),   intent(in)    :: param_file !< A structure indicating the open file
                                                             !! to parse for
                                                             !model parameter
                                                             !values.
  type(oda_incupd_CS),     pointer       :: CS   !< ODA incupd control structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling

  type(MOM_restart_CS),    pointer :: restart_CSp_tmp => NULL()

  type(directories)  :: dirs
  type(vardesc) :: u_desc, v_desc

  character(len=40)  :: mdl = "MOM_oda"  ! This module's name.
  character(len=200)  :: inc_file  ! name of the increment file

  if (.not.associated(CS)) return
  ! get the output_directory
  call Get_MOM_Input(dirs=dirs)
  if (is_root_pe()) call MOM_error(NOTE,"output increments in output_directory")

  ! get a restart structure
  call restart_init(param_file, restart_CSp_tmp)

  ! register the variables to write
  call register_restart_field(CS%Inc(1)%p, "T_inc", .true., restart_CSp_tmp, &
                              "Pot. T. increment", "degC")
  call register_restart_field(CS%Inc(2)%p, "S_inc", .true., restart_CSp_tmp, &
                              "Salinity increment", "psu")
  call register_restart_field(CS%Ref_h%p, "h_obs", .true., restart_CSp_tmp, &
                              "Observational h", "m")
  if (CS%uv_inc) then
    u_desc = var_desc("u_inc", "m s-1", "U-vel increment", hor_grid='Cu')
    v_desc = var_desc("v_inc", "m s-1", "V-vel increment", hor_grid='Cv')
    call register_restart_pair(CS%Inc_u%p, CS%Inc_v%p, u_desc, v_desc, &
              .false., restart_CSp_tmp)
  endif

  ! get the name of the output file
  call get_param(param_file, mdl, "ODA_INCUPD_OUTPUT_FILE", inc_file,&
                 "The name-root of the output file for the increment if using full fields.", &
                 default="MOM.inc")

  ! write the increments file
  call save_restart(dirs%output_directory, Time, G, restart_CSp_tmp, &
                      filename=inc_file, GV=GV) !, write_ic=.true.)

end subroutine output_oda_incupd_inc



!> Initialize diagnostics for the oda_incupd module.
subroutine init_oda_incupd_diags(Time, G, GV, diag, CS, US)
  type(time_type), target, intent(in)    :: Time !< The current model time
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),   intent(in)  :: GV !< ocean vertical grid structure
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(oda_incupd_CS),     pointer       :: CS   !< ALE sponge control structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling

  if (.not.associated(CS)) return

  CS%diag => diag
  ! These diagnostics of the state variables increments,useful for debugging the
  ! ODA code.
  CS%id_u_oda_inc = register_diag_field('ocean_model', 'u_oda_inc', diag%axesCuL, Time, &
      'Zonal velocity ODA inc.', 'm s-1')
  CS%id_v_oda_inc = register_diag_field('ocean_model', 'v_oda_inc', diag%axesCvL, Time, &
      'Meridional velocity ODA inc.', 'm s-1')
  CS%id_h_oda_inc = register_diag_field('ocean_model', 'h_oda_inc', diag%axesTL, Time, &
      'Layer Thickness ODA inc.', get_thickness_units(GV))
  CS%id_T_oda_inc = register_diag_field('ocean_model', 'T_oda_inc', diag%axesTL, Time, &
      'Temperature ODA inc.', 'degC')
  CS%id_S_oda_inc = register_diag_field('ocean_model', 'S_oda_inc', diag%axesTL, Time, &
      'Salinity ODA inc.', 'PSU')

end subroutine init_oda_incupd_diags

!> This subroutine deallocates any memory associated with the oda_incupd module.
subroutine oda_incupd_end(CS)
  type(oda_incupd_CS), pointer :: CS !< A pointer to the control structure that is
                                     !! set by a previous call to initialize_oda_incupd.

  integer :: m

  if (.not.associated(CS)) return

  do m=1,CS%fldno
    if (associated(CS%Inc(m)%p)) deallocate(CS%Inc(m)%p)
  enddo

  deallocate(CS)

end subroutine oda_incupd_end

end module MOM_oda_incupd
