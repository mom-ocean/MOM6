!> This module contains the routines used to apply incremental updates
!! from data assimilation.
!!
!! Applying incremental updates requires the following:
!! 1. initialize_oda_incupd
!! 2. set_up_oda_incupd_field (tracers) and set_up_oda_incupd_vel_field (vel)
!! 3. apply_oda_incupd
!! 5. oda_incupd_end (not being used for now)

module MOM_oda_incupd


! This file is part of MOM6. See LICENSE.md for the license.
use MOM_array_transform, only: rotate_array
use MOM_coms,          only : sum_across_PEs
use MOM_diag_mediator, only : post_data, query_averaging_enabled, register_diag_field
use MOM_diag_mediator, only : diag_ctrl
use MOM_domains, only : pass_var,pass_vector
use MOM_error_handler, only : MOM_error, FATAL, NOTE, WARNING, is_root_pe
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_remapping,     only : remapping_cs, remapping_core_h, initialize_remapping
use MOM_spatial_means, only : global_i_mean
use MOM_time_manager,  only : time_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_verticalGrid,  only : get_thickness_units

use mpp_io_mod, only : mpp_get_axis_length
use mpp_io_mod, only : axistype

implicit none ; private

#include <MOM_memory.h>


!  Publicly available functions
public set_up_oda_incupd_field, set_up_oda_incupd_vel_field
public initialize_oda_incupd, apply_oda_incupd, oda_incupd_end
public init_oda_incupd_diags,get_oda_increments
!public rotate_oda_incupd, update_oda_incupd_field

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
  integer :: IsdB      !< The starting I-index of the data domain at u/v.
  integer :: IedB      !< The ending I-index of the data domain at u/v.
  integer :: JsdB      !< The starting J-index of the data domain at u/v.
  integer :: JedB      !< The ending J-index of the data domain at h.
  integer :: fldno = 0 !< The number of fields which have already been
                       !! registered by calls to set_up_oda_incupd_field

  type(p3d) :: Var(MAX_FIELDS_)      !< Pointers to the fields that are being updated.
  type(p3d) :: Inc(MAX_FIELDS_)      !< The increments to be applied to the field
  type(p3d) :: Var_u  !< Pointer to the u-velocities. that are being updated.
  type(p3d) :: Inc_u  !< The increments to be applied to the u-velocities
  type(p3d) :: Var_v  !< Pointer to the v-velocities. that are being updated.
  type(p3d) :: Inc_v  !< The increments to be applied to the v-velocities
  type(p3d) :: Ref_h  !< Vertical grid on which the increments are provided 


  integer :: nstep_incupd          !< number of time step for full update
  integer :: ncount                !< increment time step counter
  type(remapping_cs) :: remap_cs   !< Remapping parameters and work arrays
  logical :: remap_answers_2018    !< If true, use the order of arithmetic and expressions that
                                   !! recover the answers for remapping from the end of 2018.
                                   !! Otherwise, use more robust forms of the same expressions.
  logical :: hor_regrid_answers_2018 !< If true, use the order of arithmetic for horizonal regridding
                                   !! that recovers the answers from the end of 2018.  Otherwise, use
                                   !! rotationally symmetric forms of the same expressions.

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

!> This subroutine defined the number of time step for full update, stores the layer pressure 
!! increments and initialize remap structure.
subroutine initialize_oda_incupd( G, GV, US, param_file, CS, data_h,nz_data)

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

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_oda"  ! This module's name.
  logical :: use_oda_incupd
  logical :: bndExtrapolation = .true.   ! If true, extrapolate boundaries
  logical :: default_2018_answers
  integer :: i, j, k 
  real    :: nhours_incupd, dt, dt_therm
  character(len=256) :: mesg 
  character(len=10)  :: remapScheme
  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_oda_incupd_fixed called with an associated "// &
                            "control structure.")
    return
  endif

! Set default, read and log parameters
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "ODA_INCUPD", use_oda_incupd, &
                 "If true, oda incremental updates will be applied "//&
                 "everywhere in the domain.", default=.false.)

  if (.not.use_oda_incupd) return

  allocate(CS)

  call get_param(param_file, mdl, "ODA_INCUPD_NHOURS", nhours_incupd, &
                 "Number of hours for full update (0=direct insertion).", &
                 default=3.0)
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
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
  call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", CS%remap_answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
  call get_param(param_file, mdl, "ODA_INCUPD_DATA_ONGRID", CS%incupdDataOngrid, &
                 "When defined, the incoming oda_incupd data are "//&
                 "assumed to be on the model horizontal grid " , &
                 default=.true.)

  CS%nz = GV%ke
  CS%isc = G%isc ; CS%iec = G%iec ; CS%jsc = G%jsc ; CS%jec = G%jec
  CS%isd = G%isd ; CS%ied = G%ied ; CS%jsd = G%jsd ; CS%jed = G%jed
  CS%isdB = G%isdB ; CS%iedB = G%iedB; CS%jsdB = G%jsdB ; CS%jedB = G%jedB

  ! increments on horizontal grid
  if (.not. CS%incupdDataOngrid) call MOM_error(FATAL,'initialize_oda_incupd: '// &
           'The oda_incupd code only applies ODA increments on the same horizontal grid. ')

  ! get number of timestep for full update
  if (nhours_incupd == 0) then 
     CS%nstep_incupd = 1 !! direct insertion
  else   
     CS%nstep_incupd = floor(nhours_incupd*3600./dt_therm + 0.001 ) - 1
  endif
  write(mesg,'(i12)') CS%nstep_incupd
  if (is_root_pe()) &
     call MOM_error(NOTE,"initialize_oda_incupd: Number of Timestep of inc. update:"//&
                       trim(mesg))

  ! initialize time counter
  CS%ncount=0

  ! get the vertical grid (h_obs) of the increments
  CS%nz_data = nz_data
  allocate(CS%Ref_h%p(G%isd:G%ied,G%jsd:G%jed,CS%nz_data))
  CS%Ref_h%p(:,:,:) = 0.0 ;
  do j=G%jsc,G%jec; do i=G%isc,G%iec ; do k=1,CS%nz_data
      CS%Ref_h%p(i,j,k) = data_h(i,j,k) 
  enddo;  enddo ; enddo

  ! Call the constructor for remapping control structure
  call initialize_remapping(CS%remap_cs, remapScheme, boundary_extrapolation=bndExtrapolation, &
                            answers_2018=CS%remap_answers_2018)


end subroutine initialize_oda_incupd


!> This subroutine stores theincrements at h points for the variable
!! whose address is given by f_ptr.
subroutine set_up_oda_incupd_field(sp_val, G, GV, f_ptr, CS)
  type(ocean_grid_type),   intent(in) :: G      !< Grid structure
  type(verticalGrid_type), intent(in) :: GV     !< ocean vertical grid structure
  type(oda_incupd_CS),     pointer    :: CS     !< oda_incupd control structure (in/out).
  real, dimension(SZI_(G),SZJ_(G),CS%nz_data), &
                           intent(in) :: sp_val !< increment field, it can have an
                                                !! arbitrary number of layers.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                   target, intent(in) :: f_ptr  !< Pointer to the field to be updated

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
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    do k=1,CS%nz_data
      CS%Inc(CS%fldno)%p(i,j,k) = sp_val(i,j,k)
    enddo
  enddo ; enddo
  CS%Var(CS%fldno)%p => f_ptr  ! pointer to model tracers

end subroutine set_up_oda_incupd_field


!> This subroutine stores the increments at u and v points for the variable
!! whose address is given by u_ptr and v_ptr.
subroutine set_up_oda_incupd_vel_field(u_val, v_val, G, GV, u_ptr, v_ptr, CS)
  type(ocean_grid_type),   intent(in) :: G  !< Grid structure (in).
  type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
  type(oda_incupd_CS),     pointer    :: CS !< oda incupd structure (in/out).

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                          intent(in) :: u_val !< u increment, it has arbritary number of layers but
                                              !! not to exceed the total number of model layers
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                          intent(in) :: v_val !< v increment, it has arbritary number of layers but
                                              !! not to exceed the number of model layers
  real, target, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: u_ptr !< u ptr to the field to be assimilated 
  real, target, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(in) :: v_ptr !< v ptr to the field to be assimilated

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
  CS%Var_u%p => u_ptr   ! pointer to model u-velocity

  ! store the increment/full field v profile
  allocate(CS%Inc_v%p(G%isd:G%ied,G%jsdB:G%jedB,CS%nz_data))
  CS%Inc_v%p(:,:,:) = 0.0
  do j=G%jscB,G%jecB ; do i=G%isc,G%iec
    do k=1,CS%nz_data
      CS%Inc_v%p(i,j,k) = v_val(i,j,k)
    enddo
  enddo ; enddo
  CS%Var_v%p => v_ptr  ! pointer to model v-velocity

end subroutine set_up_oda_incupd_vel_field

subroutine get_oda_increments(h, G, GV, US, CS)

  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure (in).
  type(verticalGrid_type),   intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2] (in)
  type(oda_incupd_CS),       pointer       :: CS !< A pointer to the control structure for this module
                                                 !! that is set by a previous call to initialize_oda_incupd (in).


  real, dimension(SZK_(GV)) :: tmp_val1          ! data values on the model grid
  real, allocatable, dimension(:) :: tmp_val2    ! data values remapped to increment grid
  real, allocatable, dimension(:,:,:) :: h_obs   !< h of increments
  real, allocatable, dimension(:) :: tmp_h       ! temporary array for corrected h_obs
  real, allocatable, dimension(:) :: hu_obs,hv_obs  ! A column of thicknesses at h points [H ~> m or kg m-2]
  real, dimension(SZK_(GV)) :: hu, hv            ! A column of thicknesses at h, u or v points [H ~> m or kg m-2]


  integer ::  m, i, j, k, is, ie, js, je, nz, nz_data
  integer :: isB, ieB, jsB, jeB
  real :: h_neglect, h_neglect_edge
  character(len=256) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isB = G%iscB ; ieB = G%iecB ; jsB = G%jscB ; jeB = G%jecB
  if (.not.associated(CS)) return

  if (.not.CS%remap_answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  ! get h_obs
  nz_data = CS%Inc(1)%nz_data
  allocate(h_obs(G%isd:G%ied,G%jsd:G%jed,nz_data)) ; h_obs(:,:,:)=0.0
  do k=1,nz_data  ; do j=js,je ; do i=is,ie
    h_obs(i,j,k) = CS%Ref_h%p(i,j,k)
  enddo ; enddo ; enddo
  call pass_var(h_obs,G%Domain)

  ! allocate 1-d arrays
  allocate(tmp_h(nz_data)); tmp_h(:)=0.0
  allocate(tmp_val2(nz_data)) ; tmp_val2(:) = 0.0
  allocate(hu_obs(nz_data)) ; hu_obs(:)=0.0
  allocate(hv_obs(nz_data)) ; hv_obs(:)=0.0

  ! remap t,s (on h_init) to h_obs to get increment
  tmp_val1(:)=0.0
  do j=js,je ; do i=is,ie
     ! account for the different SSH
     tmp_h(1:nz_data)=(sum(h(i,j,1:nz))/sum(h_obs(i,j,1:nz_data)))*h_obs(i,j,1:nz_data)
     do m=1,CS%fldno
       ! get tracer
       tmp_val1(1:nz) = CS%Var(m)%p(i,j,1:nz)
       ! remap tracer on h_obs
       call remapping_core_h(CS%remap_cs, nz     ,     h(i,j,1:nz     ), tmp_val1, &
                                          nz_data, tmp_h(    1:nz_data), tmp_val2, &
                             h_neglect, h_neglect_edge)
       ! get increment from full field on h_obs
       CS%Inc(m)%p(i,j,1:nz_data) = CS%Inc(m)%p(i,j,1:nz_data) - tmp_val2(1:nz_data)
     enddo
  enddo; enddo

  ! remap u to h_obs to get increment
  if (CS%uv_inc) then 
    hu(:)=0.0
    do j=js,je ; do i=isB,ieB
      ! get u-velocity
      tmp_val1(1:nz) = CS%Var_u%p(i  ,j,1:nz)
      ! get the h and h_obs at u points
      hu(1:nz)          = 0.5*(    h(i,j,1:nz     )+    h(i+1,j,1:nz     ))
      hu_obs(1:nz_data) = 0.5*(h_obs(i,j,1:nz_data)+h_obs(i+1,j,1:nz_data))
      ! account for the different SSH 
      hu_obs(1:nz_data) = (sum(hu(1:nz))/sum(hu_obs(1:nz_data))) * hu_obs(1:nz_data)
      ! remap model u on hu_obs
      call remapping_core_h(CS%remap_cs, nz     ,     hu(1:nz     ), tmp_val1, &
                                         nz_data, hu_obs(1:nz_data), tmp_val2, &
                            h_neglect, h_neglect_edge)
      ! get increment from full field on h_obs
      CS%Inc_u%p(i,j,1:nz_data) = CS%Inc_u%p(i,j,1:nz_data) - tmp_val2(1:nz_data)
    enddo; enddo
  
    ! remap v to h_obs to get increment
    hv(:) = 0.0;
    do j=jsB,jeB ; do i=is,ie
      ! get v-velocity
      tmp_val1(1:nz) = CS%Var_v%p(i  ,j,1:nz)
      ! get the h and h_obs at v points      
      hv(1:nz)          = 0.5*(    h(i,j,1:nz     )+    h(i,j+1,1:nz     ))
      hv_obs(1:nz_data) = 0.5*(h_obs(i,j,1:nz_data)+h_obs(i,j+1,1:nz_data))
      ! account for the different SSH
      hv_obs(1:nz_data) = (sum(hv(1:nz))/sum(hv_obs(1:nz_data))) *hv_obs(1:nz_data)
      ! remap model v on hv_obs 
      call remapping_core_h(CS%remap_cs, nz     ,     hv(1:nz     ), tmp_val1, &
                                         nz_data, hv_obs(1:nz_data), tmp_val2, &
                            h_neglect, h_neglect_edge)
      ! get increment from full field on h_obs
      CS%Inc_v%p(i,j,1:nz_data) = CS%Inc_v%p(i,j,1:nz_data) - tmp_val2(1:nz_data)
    enddo; enddo
  endif ! uv_inc
 
  ! deallocate  arrays
  deallocate(tmp_h,tmp_val2,hu_obs,hv_obs)
  deallocate(h_obs)

end subroutine 

!> This subroutine applies oda increments to layers thicknesses, temp, salt, U
!!  and V everywhere .
subroutine apply_oda_incupd(h, dt, G, GV, US, CS)
  type(ocean_grid_type),     intent(in) :: G  !< The ocean's grid structure (in).
  type(verticalGrid_type),   intent(in)    :: GV !< ocean vertical grid structure
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                             intent(in) :: h  !< Layer thickness [H ~> m or kg m-2] (in)
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

  integer ::  m, i, j, k, is, ie, js, je, nz, nz_data
  integer :: isB, ieB, jsB, jeB
  integer :: ncount ! time step counter
  real :: q
  real :: h_neglect, h_neglect_edge
  character(len=256) :: mesg  

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isB = G%iscB ; ieB = G%iecB ; jsB = G%jscB ; jeB = G%jecB
  if (.not.associated(CS)) return

  ! update counter
  CS%ncount = CS%ncount+1
  q = 1.0/CS%nstep_incupd

  ! no assimilation after CS%step_incupd
  if (CS%ncount > CS%nstep_incupd) then
    if (is_root_pe()) call MOM_error(NOTE,"ended updating fields with increments. ")
    return
  endif !ncount>CS%nstep_incupd


  ! print out increments
  write(mesg,'(i8)') CS%ncount
  if (is_root_pe()) call MOM_error(NOTE,"updating fields with increments ncount:"//trim(mesg))
  write(mesg,'(f10.8)') q
  if (is_root_pe()) call MOM_error(NOTE,"updating fields with weight q:"//trim(mesg))

  if (.not.CS%remap_answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  ! get h_obs
  nz_data = CS%Inc(1)%nz_data
  allocate(h_obs(G%isd:G%ied,G%jsd:G%jed,nz_data)) ; h_obs(:,:,:)=0.0
  do k=1,nz_data  ; do j=js,je ; do i=is,ie
    h_obs(i,j,k) = CS%Ref_h%p(i,j,k)
  enddo ; enddo ; enddo
  call pass_var(h_obs,G%Domain)

  ! allocate 1-d array
  allocate(tmp_h(nz_data)); tmp_h(:)=0.0
  allocate(tmp_val2(nz_data))
  allocate(hu_obs(nz_data)) ; hu_obs(:)=0.0
  allocate(hv_obs(nz_data)) ; hv_obs(:)=0.0

  ! add increments to tracers
  tmp_val1(:)=0.0
  tmp_t(:,:,:) = 0.0 ; tmp_s(:,:,:) = 0.0 ! diagnostics
  do j=js,je ; do i=is,ie
    ! account for the different SSH
    tmp_h(1:nz_data)=(sum(h(i,j,1:nz))/sum(h_obs(i,j,1:nz_data)))*h_obs(i,j,1:nz_data)
    do m=1,CS%fldno
      ! get increment
      tmp_val2(1:nz_data) = CS%Inc(m)%p(i,j,1:nz_data)
      ! remap increment profile on model h
      call remapping_core_h(CS%remap_cs, nz_data, tmp_h(    1:nz_data),tmp_val2, &
                                         nz     ,     h(i,j,1:nz     ),tmp_val1, &
                            h_neglect, h_neglect_edge)
      ! add increment to tracer on model h
      CS%Var(m)%p(i,j,1:nz) = CS%Var(m)%p(i,j,1:nz) + q*tmp_val1(1:nz)
      ! store increment for diagnostics
      if (m == 1) then
        tmp_t(i,j,1:nz) = tmp_val1(1:nz)
      else
        tmp_s(i,j,1:nz) = tmp_val1(1:nz)
      endif
    enddo
    ! bound salinity values ! check if it is correct to do that or if it hides
    ! other problems ...
    do k=1,nz
      CS%Var(2)%p(i,j,k) = max(0.0 , CS%Var(2)%p(i,j,k))
    enddo
  enddo; enddo


  ! add u and v increments
  if (CS%uv_inc) then

    ! add increments to u
    hu(:) = 0.0
    tmp_u(:,:,:) = 0.0 ! diagnostics
    do j=js,je ; do i=isB,ieB
      ! get increment
      tmp_val2(1:nz_data) = CS%Inc_u%p(i,j,1:nz_data)
      ! get the h and h_obs at u points
      hu_obs(1:nz_data) = 0.5*(h_obs(i,j,1:nz_data)+h_obs(i+1,j,1:nz_data))
      hu(1:nz)          = 0.5*(    h(i,j,1:nz     )+    h(i+1,j,1:nz     ))
      ! account for different SSH
      hu_obs(1:nz_data) = (sum(hu(1:nz))/sum(hu_obs(1:nz_data))) * hu_obs(1:nz_data)
      ! remap increment profile on hu
      call remapping_core_h(CS%remap_cs, nz_data, hu_obs(1:nz_data), tmp_val2, &
                                         nz     ,     hu(1:nz     ), tmp_val1, &
                            h_neglect, h_neglect_edge)
      ! add increment to u-velocity on hu
      CS%Var_u%p(i,j,1:nz) = CS%Var_u%p(i,j,1:nz) + q*tmp_val1(1:nz)
      ! store increment for diagnostics
      tmp_u(i,j,1:nz) = tmp_val1(1:nz)
    enddo; enddo

    ! add increments to v
    hv(:) = 0.0
    tmp_v(:,:,:) = 0.0 ! diagnostics
    do j=jsB,jeB ; do i=is,ie
      ! get increment
      tmp_val2(1:nz_data) = CS%Inc_v%p(i,j,1:nz_data)
      ! get the h and h_obs at v points
      hv_obs(1:nz_data) = 0.5 * (h_obs(i,j,1:nz_data)+h_obs(i,j+1,1:nz_data))
      hv(1:nz)          = 0.5 * (    h(i,j,1:nz     )+    h(i,j+1,1:nz     ))
      ! account for different SSH
      hv_obs(1:nz_data) = (sum(hv(1:nz))/sum(hv_obs(1:nz_data))) * hv_obs(1:nz_data)
      ! remap increment profile on hv
      call remapping_core_h(CS%remap_cs, nz_data, hv_obs(1:nz_data), tmp_val2, &
                                         nz     ,     hv(1:nz     ), tmp_val1, &
                            h_neglect, h_neglect_edge)
      ! add increment to v-velocity on hv
      CS%Var_v%p(i,j,1:nz) = CS%Var_v%p(i,j,1:nz) + q*tmp_val1(1:nz)
      ! store increment for diagnostics
      tmp_v(i,j,1:nz) = tmp_val1(1:nz)
    enddo; enddo

  endif ! uv_inc 

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

!> Rotate the incr. fields from the input to the model index map.
!> Not implemented yet
!subroutine rotate_oda_incupd(oda_incupd_in, G_in, oda_incupd, G, GV, turns, param_file)
!  type(oda_incupd_CS),     intent(in) :: oda_incupd_in !< The control structure for this module with the
!                                                   !! original grid rotation
!  type(ocean_grid_type),   intent(in) :: G_in      !< The ocean's grid structure with the original rotation.
!  type(oda_incupd_CS),     pointer    :: oda_incup !< A pointer to the control that will be set up with
!                                                   !! the new grid rotation
!  type(ocean_grid_type),   intent(in) :: G         !< The ocean's grid structure with the new rotation.
!  type(verticalGrid_type), intent(in) :: GV        !< The ocean's vertical grid structure
!  integer,                 intent(in) :: turns     !< The number of 90-degree turns between grids
!  type(param_file_type),   intent(in) :: param_file !< A structure indicating the open file
!                                                   !! to parse for model parameter values.
!
!  ! First part: Index construction
!  !   1. Reconstruct Iresttime(:,:) from sponge_in
!  !   2. rotate Iresttime(:,:)
!  !   3. Call initialize_oda_incupd using new grid and rotated Iresttime(:,:)
!  ! All the index adjustment should follow from the Iresttime rotation
!
!  real, dimension(:,:), allocatable :: Iresttime_in, Iresttime
!  real, dimension(:,:,:), allocatable :: data_h_in, data_h
!  real, dimension(:,:,:), allocatable :: sp_val_in, sp_val
!  real, dimension(:,:,:), pointer :: sp_ptr => NULL()
!  integer :: c, c_i, c_j
!  integer :: k, nz_data
!  integer :: n
!  logical :: fixed_sponge
!
!  fixed_sponge = .not. sponge_in%time_varying_sponges
!  ! NOTE: nz_data is only conditionally set when fixed_sponge is true.
!
!  allocate(Iresttime_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed))
!  allocate(Iresttime(G%isd:G%ied, G%jsd:G%jed))
!  Iresttime_in(:,:) = 0.0
!
!  if (fixed_sponge) then
!    nz_data = sponge_in%nz_data
!    allocate(data_h_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed, nz_data))
!    allocate(data_h(G%isd:G%ied, G%jsd:G%jed, nz_data))
!    data_h_in(:,:,:) = 0.
!  endif
!
!  ! Re-populate the 2D Iresttime and data_h arrays on the original grid
!  do c=1,sponge_in%num_col
!    c_i = sponge_in%col_i(c)
!    c_j = sponge_in%col_j(c)
!    Iresttime_in(c_i, c_j) = sponge_in%Iresttime_col(c)
!    if (fixed_sponge) then
!      do k = 1, nz_data
!        data_h(c_i, c_j, k) = sponge_in%Ref_h%p(k,c)
!      enddo
!    endif
!  enddo
!
!  call rotate_array(Iresttime_in, turns, Iresttime)
!  if (fixed_sponge) then
!    call rotate_array(data_h_in, turns, data_h)
!    call initialize_oda_incupd_fixed(Iresttime, G, GV, param_file, sponge, &
!                                     data_h, nz_data)
!  else
!    call initialize_oda_incupd_varying(Iresttime, G, GV, param_file, sponge)
!  endif
!
!  deallocate(Iresttime_in)
!  deallocate(Iresttime)
!  if (fixed_sponge) then
!    deallocate(data_h_in)
!    deallocate(data_h)
!  endif
!
!  ! Second part: Provide rotated fields for which relaxation is applied
!
!  sponge%fldno = sponge_in%fldno
!
!  if (fixed_sponge) then
!    allocate(sp_val_in(G_in%isd:G_in%ied, G_in%jsd:G_in%jed, nz_data))
!    allocate(sp_val(G%isd:G%ied, G%jsd:G%jed, nz_data))
!  endif
!
!  do n=1,sponge_in%fldno
!    ! Assume that tracers are pointers and are remapped in other functions(?)
!    sp_ptr => sponge_in%var(n)%p
!    if (fixed_sponge) then
!      sp_val_in(:,:,:) = 0.0
!      do c = 1, sponge_in%num_col
!        c_i = sponge_in%col_i(c)
!        c_j = sponge_in%col_j(c)
!        do k = 1, nz_data
!          sp_val_in(c_i, c_j, k) = sponge_in%Ref_val(n)%p(k,c)
!        enddo
!      enddo
!
!      call rotate_array(sp_val_in, turns, sp_val)
!
!      ! NOTE: This points sp_val with the unrotated field.  See note below.
!      call set_up_oda_incupd_field(sp_val, G, GV, sp_ptr, sponge)
!
!      deallocate(sp_val_in)
!    else
!      ! We don't want to repeat FMS init in set_up_oda_incupd_field_varying()
!      ! (time_interp_external_init, init_external_field, etc), so we manually
!      ! do a portion of this function below.
!      sponge%Ref_val(n)%id = sponge_in%Ref_val(n)%id
!      sponge%Ref_val(n)%num_tlevs = sponge_in%Ref_val(n)%num_tlevs
!
!      nz_data = sponge_in%Ref_val(n)%nz_data
!      sponge%Ref_val(n)%nz_data = nz_data
!
!      allocate(sponge%Ref_val(n)%p(nz_data, sponge_in%num_col))
!      allocate(sponge%Ref_val(n)%h(nz_data, sponge_in%num_col))
!      sponge%Ref_val(n)%p(:,:) = 0.0
!      sponge%Ref_val(n)%h(:,:) = 0.0
!
!      ! TODO: There is currently no way to associate a generic field pointer to
!      !   its rotated equivalent without introducing a new data structure which
!      !   explicitly tracks the pairing.
!      !
!      !   As a temporary fix, we store the pointer to the unrotated field in
!      !   the rotated sponge, and use this reference to replace the pointer
!      !   to the rotated field update_oda_incupd field.
!      !
!      !   This makes a lot of unverifiable assumptions, and should not be
!      !   considered the final solution.
!      sponge%var(n)%p => sp_ptr
!    endif
!  enddo
!
!  ! TODO: var_u and var_v sponge dampling is not yet supported.
!  if (associated(sponge_in%var_u%p) .or. associated(sponge_in%var_v%p)) &
!    call MOM_error(FATAL, "Rotation of ALE sponge velocities is not yet " &
!      // "implemented.")
!
!  ! Transfer any existing diag_CS reference pointer
!  sponge%diag => sponge_in%diag
!
!  ! NOTE: initialize_oda_incupd_* resolves remap_cs
!end subroutine rotate_oda_incupd


!> Scan the oda_incupd variables and replace a prescribed pointer to a new value.
!> Not implemented yet
! TODO: This function solely exists to replace field pointers in the sponge
!   after rotation.  This function is part of a temporary solution until
!   something more robust is developed.
!subroutine update_oda_incupd_field(sponge, p_old, G, GV, p_new)
!  type(oda_incupd_CS),     pointer    :: sponge !< A pointer to the control structure for this module
!                                               !! that is set by a previous call to initialize_oda_incupd.
!  real, dimension(:,:,:), &
!                   target, intent(in) :: p_old !< The previous array of target values
!  type(ocean_grid_type),   intent(in) :: G     !< The updated ocean grid structure
!  type(verticalGrid_type), intent(in) :: GV    !< ocean vertical grid structure
!  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
!                   target, intent(in) :: p_new !< The new array of target values
!
!  integer :: n
!
!  do n=1,sponge%fldno
!    if (associated(sponge%var(n)%p, p_old)) sponge%var(n)%p => p_new
!  enddo
!
!end subroutine update_oda_incupd_field


! GMM: I could not find where sponge_end is being called, but I am keeping
!  oda_incupd_end here so we can add that if needed.
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
