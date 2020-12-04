!> This module contains the routines used to set up a
!! dynamically passive tracer.
!! Set up and use passive tracers requires the following:
!! (1) register_RGC_tracer
!! (2) apply diffusion, physics/chemistry and advect the tracer

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Elizabeth Yankovsky, June 2019                                  *
!*********+*********+*********+*********+*********+*********+***********

module RGC_tracer

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_hor_index, only : hor_index_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_restart, only :  MOM_restart_CS
use MOM_ALE_sponge, only : set_up_ALE_sponge_field, ALE_sponge_CS, get_ALE_sponge_nz_data
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_open_boundary, only : ocean_OBC_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!< Publicly available functions
public register_RGC_tracer, initialize_RGC_tracer
public RGC_tracer_column_physics, RGC_tracer_end

integer, parameter :: NTR = 1 !< The number of tracers in this module.

!> tracer control structure
type, public :: RGC_tracer_CS ; private
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  character(len = 200) :: tracer_IC_file !< The full path to the IC file, or " " to initialize internally.
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry.
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this package.
  real, pointer :: tr_aux(:,:,:,:) => NULL() !< The masked tracer concentration.
  real :: land_val(NTR) = -1.0 !< The value of tr used where land is masked out.
  real :: lenlat           !< the latitudinal or y-direction length of the domain.
  real :: lenlon           !< the longitudinal or x-direction length of the domain.
  real :: CSL              !< The length of the continental shelf (x dir, km)
  real :: lensponge        !< the length of the sponge layer.
  logical :: mask_tracers  !< If true, tracers are masked out in massless layers.
  logical :: use_sponge    !< If true, sponges may be applied somewhere in the domain.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the timing of diagnostic output.
  type(vardesc) :: tr_desc(NTR) !< Descriptions and metadata for the tracers.
end type RGC_tracer_CS

contains


!> This subroutine is used to register tracer fields
function register_RGC_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !<A structure indicating the open file to parse
                                                 !! for model parameter values.
  type(RGC_tracer_CS),        pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module (in/out).
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer to the tracer registry.
  type(MOM_restart_CS),       pointer    :: restart_CS !< A pointer to the restart control structure.

  character(len=80)  :: name, longname
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "RGC_tracer" ! This module's name.
  character(len=200) :: inputdir
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_RGC_tracer
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "RGC_register_tracer called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "RGC_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial \n"//&
                 "conditions for the RGC tracers, or blank to initialize \n"//&
                 "them internally.", default=" ")
  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    CS%tracer_IC_file = trim(inputdir)//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/RGC_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mdl, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. \n"//&
                 "The exact location and properties of those sponges are \n"//&
                 "specified from MOM_initialization.F90.", default=.false.)

  call get_param(param_file, mdl, "LENLAT", CS%lenlat, &
                 "The latitudinal or y-direction length of the domain", &
                 fail_if_missing=.true., do_not_log=.true.)

  call get_param(param_file, mdl, "LENLON", CS%lenlon, &
                 "The longitudinal or x-direction length of the domain", &
                 fail_if_missing=.true., do_not_log=.true.)

  call get_param(param_file, mdl, "CONT_SHELF_LENGTH", CS%CSL, &
                 "The length of the continental shelf (x dir, km).", &
                 default=15.0)

  call get_param(param_file, mdl, "LENSPONGE", CS%lensponge, &
                 "The length of the sponge layer (km).", &
                 default=10.0)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR)) ; CS%tr(:,:,:,:) = 0.0
  if (CS%mask_tracers) then
    allocate(CS%tr_aux(isd:ied,jsd:jed,nz,NTR)) ; CS%tr_aux(:,:,:,:) = 0.0
  endif

  do m=1,NTR
    if (m < 10) then ; write(name,'("tr_RGC",I1.1)') m
    else ; write(name,'("tr_RGC",I2.2)') m ; endif
    write(longname,'("Concentration of RGC Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)

    ! This is needed to force the compiler not to do a copy in the registration calls.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for horizontal advection & diffusion.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                         name=name, longname=longname, units="kg kg-1", &
                         registry_diags=.true., flux_units="kg/s", &
                         restart_CS=restart_CS)
  enddo

  CS%tr_Reg => tr_Reg
  register_RGC_tracer = .true.
end function register_RGC_tracer

!> Initializes the NTR tracer fields in tr(:,:,:,:)
!! and it sets up the tracer output.
subroutine initialize_RGC_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                    layer_CSp, sponge_CSp)

  type(ocean_grid_type),   intent(in) :: G !< Grid structure.
  type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure.
  logical,                 intent(in) :: restart !< .true. if the fields have already
                                             !! been read from a restart file.
  type(time_type), target, intent(in) :: day !< Time of the start of the run.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h !< Layer thickness, in m or kg m-2.
  type(diag_ctrl), target, intent(in) :: diag !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),    pointer    :: OBC !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used. This is not being used for now.
  type(RGC_tracer_CS),     pointer    :: CS  !< The control structure returned by a previous
                                             !!   call to RGC_register_tracer.
  type(sponge_CS),         pointer    :: layer_CSp  !< A pointer to the control structure
  type(ALE_sponge_CS),     pointer    :: sponge_CSp !< A pointer to the control structure for the
                                             !! sponges, if they are in use.  Otherwise this may be unassociated.

  real, allocatable :: temp(:,:,:)
  real, pointer, dimension(:,:,:) :: &
    OBC_tr1_u => NULL(), & ! These arrays should be allocated and set to
    OBC_tr1_v => NULL()    ! specify the values of tracer 1 that should come
                           ! in through u- and v- points through the open
                           ! boundary conditions, in the same units as tr.
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for tracer fluxes, usually
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  real :: h_neglect         ! A thickness that is so small it is usually lost
                            ! in roundoff and can be neglected [H ~> m or kg-2].
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB
  integer :: nzdata

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  h_neglect = GV%H_subroundoff

  CS%Time => day
  CS%diag => diag

  if (.not.restart) then
    if (len_trim(CS%tracer_IC_file) >= 1) then
      !  Read the tracer concentrations from a netcdf file.
      if (.not.file_exists(CS%tracer_IC_file, G%Domain)) &
        call MOM_error(FATAL, "RGC_initialize_tracer: Unable to open "// &
                        CS%tracer_IC_file)
      do m=1,NTR
        call query_vardesc(CS%tr_desc(m), name, caller="initialize_RGC_tracer")
        call MOM_read_data(CS%tracer_IC_file, trim(name), CS%tr(:,:,:,m), G%Domain)
      enddo
    else
      do m=1,NTR
        do k=1,nz ; do j=js,je ; do i=is,ie
          CS%tr(i,j,k,m) = 0.0
        enddo ; enddo ; enddo
      enddo
      m=1
      do j=js,je ; do i=is,ie
         !set tracer to 1.0 in the surface of the continental shelf
         if (G%geoLonT(i,j) <= (CS%CSL)) then
            CS%tr(i,j,1,m) = 1.0 !first layer
         endif
      enddo ; enddo

    endif
  endif ! restart

  if ( CS%use_sponge ) then
!  If sponges are used, this damps values to zero in the offshore boundary.
!  For any tracers that are not damped in the sponge, the call
! to set_up_sponge_field can simply be omitted.
    if (associated(sponge_CSp)) then !ALE mode
      nzdata = get_ALE_sponge_nz_data(sponge_CSp)
      if (nzdata>0) then
        allocate(temp(G%isd:G%ied,G%jsd:G%jed,nzdata))
        do k=1,nzdata ; do j=js,je ; do i=is,ie
          if (G%geoLonT(i,j) >= (CS%lenlon - CS%lensponge) .AND. G%geoLonT(i,j) <= CS%lenlon) then
            temp(i,j,k) = 0.0
          endif
        enddo ; enddo ; enddo
        do m=1,1
        ! This is needed to force the compiler not to do a copy in the sponge calls.
          tr_ptr => CS%tr(:,:,:,m)
          call set_up_ALE_sponge_field(temp, G, GV, tr_ptr, sponge_CSp)
        enddo
        deallocate(temp)
      endif

    elseif (associated(layer_CSp)) then !layer mode
      if (nz>0) then
        allocate(temp(G%isd:G%ied,G%jsd:G%jed,nz))
        do k=1,nz ; do j=js,je ; do i=is,ie
          if (G%geoLonT(i,j) >= (CS%lenlon - CS%lensponge) .AND. G%geoLonT(i,j) <= CS%lenlon) then
            temp(i,j,k) = 0.0
          endif
        enddo ; enddo ; enddo
        do m=1,1
          tr_ptr => CS%tr(:,:,:,m)
          call set_up_sponge_field(temp, tr_ptr, G, GV, nz, layer_CSp)
        enddo
        deallocate(temp)
      endif
    else
      call MOM_error(FATAL, "RGC_initialize_tracer: "// &
        "The pointer to sponge_CSp must be associated if SPONGE is defined.")
    endif !selecting mode/calling error if no pointer
  endif !using sponge

end subroutine initialize_RGC_tracer

!> This subroutine applies diapycnal diffusion and any other column
!! tracer physics or chemistry to the tracers from this file.
!! This is a simple example of a set of advected passive tracers.
subroutine RGC_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, US, CS, &
                              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),                 intent(in) :: G !< The ocean's grid structure.
  type(verticalGrid_type),               intent(in) :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to any possible
                                              !! forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s].
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(RGC_tracer_CS),     pointer    :: CS   !< The control structure returned by a previous call.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can be
                                               !! fluxed out of the top layer in a timestep [nondim].
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which fluxes
                                               !! can be applied [H ~> m or kg m-2].

! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  real :: in_flux(SZI_(G),SZJ_(G),2)  ! total amount of tracer to be injected

  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

  in_flux(:,:,:) = 0.0
  m=1
  do j=js,je ; do i=is,ie
    ! set tracer to 1.0 in the surface of the continental shelf
    if (G%geoLonT(i,j) <= (CS%CSL)) then
      CS%tr(i,j,1,m) = 1.0 !first layer
    endif
  enddo ; enddo

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,NTR
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo;
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m) , dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth, in_flux(:,:,m))

      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,NTR
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

end subroutine RGC_tracer_column_physics

subroutine RGC_tracer_end(CS)
  type(RGC_tracer_CS), pointer :: CS !< The control structure returned by a previous call to RGC_register_tracer.
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine RGC_tracer_end

end module RGC_tracer
