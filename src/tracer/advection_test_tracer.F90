!> This tracer package is used to test advection schemes
module advection_test_tracer

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : file_exists, read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : query_initialized, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type

use coupler_types_mod, only : coupler_type_set_data, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_advection_test_tracer, initialize_advection_test_tracer
public advection_test_tracer_surface_state, advection_test_tracer_end
public advection_test_tracer_column_physics, advection_test_stock

integer, parameter :: NTR = 11  !< The number of tracers in this module.

!> The control structure for the advect_test_tracer module
type, public :: advection_test_tracer_CS ; private
  integer :: ntr = NTR                 !< Number of tracers in this module
  logical :: coupled_tracers = .false. !< These tracers are not offered to the coupler.
  character(len=200) :: tracer_IC_file !< The full path to the IC file, or " " to initialize internally.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM tracer registry
  real, pointer :: tr(:,:,:,:) => NULL() !< The array of tracers used in this subroutine, in g m-3?
  real :: land_val(NTR) = -1.0 !< The value of tr used where land is masked out.
  logical :: use_sponge    !< If true, sponges may be applied somewhere in the domain.
  logical :: tracers_may_reinit !< If true, the tracers may be set up via the initialization code if
                           !! they are not found in the restart files.  Otherwise it is a fatal error
                           !! if the tracers are not found in the restart files of a restarted run.
  real :: x_origin !< Parameters describing the test functions
  real :: x_width  !< Parameters describing the test functions
  real :: y_origin !< Parameters describing the test functions
  real :: y_width  !< Parameters describing the test functions

  integer, dimension(NTR) :: ind_tr !< Indices returned by aof_set_coupler_flux if it is used and
                   !! the surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure.

  type(vardesc) :: tr_desc(NTR) !< Descriptions and metadata for the tracers
end type advection_test_tracer_CS

contains

!> Register tracer fields and subroutines to be used with MOM.
function register_advection_test_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),        intent(in) :: HI   !< A horizontal index type structure
  type(verticalGrid_type),     intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),       intent(in) :: param_file !< A structure to parse for run-time parameters
  type(advection_test_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                                !! call to register_advection_test_tracer.
  type(tracer_registry_type),  pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                !! structure for the tracer advection and
                                                !! diffusion module
  type(MOM_restart_CS),        pointer    :: restart_CS !< A pointer to the restart control structure

  ! Local variables
  character(len=80)  :: name, longname
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "advection_test_tracer" ! This module's name.
  character(len=200) :: inputdir
  character(len=48)  :: flux_units ! The units for tracer fluxes, usually
                            ! kg(tracer) kg(water)-1 m3 s-1 or kg(tracer) s-1.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_advection_test_tracer
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_advection_test_tracer called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "ADVECTION_TEST_X_ORIGIN", CS%x_origin, &
        "The x-coorindate of the center of the test-functions.", default=0.)
  call get_param(param_file, mdl, "ADVECTION_TEST_Y_ORIGIN", CS%y_origin, &
        "The y-coorindate of the center of the test-functions.", default=0.)
  call get_param(param_file, mdl, "ADVECTION_TEST_X_WIDTH", CS%x_width, &
        "The x-width of the test-functions.", default=0.)
  call get_param(param_file, mdl, "ADVECTION_TEST_Y_WIDTH", CS%y_width, &
        "The y-width of the test-functions.", default=0.)
  call get_param(param_file, mdl, "ADVECTION_TEST_TRACER_IC_FILE", CS%tracer_IC_file, &
                 "The name of a file from which to read the initial "//&
                 "conditions for the tracers, or blank to initialize "//&
                 "them internally.", default=" ")

  if (len_trim(CS%tracer_IC_file) >= 1) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%tracer_IC_file = trim(slasher(inputdir))//trim(CS%tracer_IC_file)
    call log_param(param_file, mdl, "INPUTDIR/ADVECTION_TEST_TRACER_IC_FILE", &
                   CS%tracer_IC_file)
  endif
  call get_param(param_file, mdl, "SPONGE", CS%use_sponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if the tracers are not found in the "//&
                 "restart files of a restarted run.", default=.false.)


  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR)) ; CS%tr(:,:,:,:) = 0.0

  do m=1,NTR
    if (m < 10) then ; write(name,'("tr",I1.1)') m
    else ; write(name,'("tr",I2.2)') m ; endif
    write(longname,'("Concentration of Tracer ",I2.2)') m
    CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)
    if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
    else ; flux_units = "kg s-1" ; endif


    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                         name=name, longname=longname, units="kg kg-1", &
                         registry_diags=.true., flux_units=flux_units, &
                         restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_advection_test_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_advection_test_tracer = .true.
end function register_advection_test_tracer

!>   Initializes the NTR tracer fields in tr(:,:,:,:) and it sets up the tracer output.
subroutine initialize_advection_test_tracer(restart, day, G, GV, h,diag, OBC, CS, &
                                            sponge_CSp)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(time_type),            target, intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(advection_test_tracer_CS),     pointer    :: CS !< The control structure returned by a previous
                                                       !! call to register_advection_test_tracer.
  type(sponge_CS),                    pointer    :: sponge_CSp !< Pointer to the control structure for the sponges.

  ! Local variables
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
                            ! in roundoff and can be neglected [H ~> m or kg m-2].
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB
  real :: tmpx, tmpy, locx, locy

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  h_neglect = GV%H_subroundoff

  CS%diag => diag
  CS%ntr = NTR
  do m=1,NTR
    call query_vardesc(CS%tr_desc(m), name=name, &
                       caller="initialize_advection_test_tracer")
    if ((.not.restart) .or. (CS%tracers_may_reinit .and. .not. &
        query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then
      do k=1,nz ; do j=js,je ; do i=is,ie
        CS%tr(i,j,k,m) = 0.0
      enddo ; enddo ; enddo
      k=1 ! Square wave
      do j=js,je ; do i=is,ie
        if (abs(G%geoLonT(i,j)-CS%x_origin)<0.5*CS%x_width .and. &
            abs(G%geoLatT(i,j)-CS%y_origin)<0.5*CS%y_width) CS%tr(i,j,k,m) = 1.0
      enddo ; enddo
      k=2 ! Triangle wave
      do j=js,je ; do i=is,ie
        locx=abs(G%geoLonT(i,j)-CS%x_origin)/CS%x_width
        locy=abs(G%geoLatT(i,j)-CS%y_origin)/CS%y_width
        CS%tr(i,j,k,m) = max(0.0, 1.0-locx)*max(0.0, 1.0-locy)
      enddo ; enddo
      k=3 ! Cosine bell
      do j=js,je ; do i=is,ie
        locx=min(1.0, abs(G%geoLonT(i,j)-CS%x_origin)/CS%x_width)*(acos(0.0)*2.)
        locy=min(1.0, abs(G%geoLatT(i,j)-CS%y_origin)/CS%y_width)*(acos(0.0)*2.)
        CS%tr(i,j,k,m) = (1.0+cos(locx))*(1.0+cos(locy))*0.25
      enddo ; enddo
      k=4 ! Cylinder
      do j=js,je ; do i=is,ie
        locx=abs(G%geoLonT(i,j)-CS%x_origin)/CS%x_width
        locy=abs(G%geoLatT(i,j)-CS%y_origin)/CS%y_width
        if (locx**2+locy**2<=1.0) CS%tr(i,j,k,m) = 1.0
      enddo ; enddo
      k=5 ! Cut cylinder
      do j=js,je ; do i=is,ie
        locx=(G%geoLonT(i,j)-CS%x_origin)/CS%x_width
        locy=(G%geoLatT(i,j)-CS%y_origin)/CS%y_width
        if (locx**2+locy**2<=1.0) CS%tr(i,j,k,m) = 1.0
        if (locx>0.0.and.abs(locy)<0.2) CS%tr(i,j,k,m) = 0.0
      enddo ; enddo
    endif ! restart
  enddo


end subroutine initialize_advection_test_tracer


!>   Applies diapycnal diffusion and any other column tracer physics or chemistry to the tracers
!! from this package.  This is a simple example of a set of advected passive tracers.
subroutine advection_test_tracer_column_physics(h_old, h_new,  ea,  eb, fluxes, dt, G, GV, US, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(advection_test_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                              !! call to register_advection_test_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  real :: b1(SZI_(G))          ! b1 and c1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))  ! tridiagonal solver.
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%ntr
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,NTR
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

end subroutine advection_test_tracer_column_physics

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
!! This particular tracer package does not report anything back to the coupler.
subroutine advection_test_tracer_surface_state(state, h, G, CS)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure.
  type(surface),          intent(inout) :: state !< A structure containing fields that
                                              !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                          intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(advection_test_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                              !! call to register_advection_test_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,CS%ntr
      !   This call loads the surface values into the appropriate array in the
      ! coupler-type structure.
      call coupler_type_set_data(CS%tr(:,:,1,m), CS%ind_tr(m), ind_csurf, &
                   state%tr_fields, idim=(/isd, is, ie, ied/), &
                   jdim=(/jsd, js, je, jed/) )
    enddo
  endif

end subroutine advection_test_tracer_surface_state

!> Calculate the mass-weighted integral of all tracer stocks, returning the number of stocks it has calculated.
!!  If the stock_index is present, only the stock corresponding to that coded index is returned.
function advection_test_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G      !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(:),                 intent(out)   :: stocks !< the mass-weighted integrated amount of each
                                                              !! tracer, in kg times concentration units [kg conc].
  type(verticalGrid_type),            intent(in)    :: GV     !< The ocean's vertical grid structure
  type(advection_test_tracer_CS),     pointer       :: CS     !< The control structure returned by a previous
                                                              !! call to register_advection_test_tracer.
  character(len=*), dimension(:),     intent(out)   :: names  !< the names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units  !< the units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< the coded index of a specific stock being sought.
  integer                                           :: advection_test_stock !< the number of stocks calculated here.

  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  advection_test_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="advection_test_stock")
    stocks(m) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      stocks(m) = stocks(m) + CS%tr(i,j,k,m) * &
                             (G%mask2dT(i,j) * G%US%L_to_m**2*G%areaT(i,j) * h(i,j,k))
    enddo ; enddo ; enddo
    stocks(m) = GV%H_to_kg_m2 * stocks(m)
  enddo
  advection_test_stock = CS%ntr

end function advection_test_stock

!> Deallocate memory associated with this module
subroutine advection_test_tracer_end(CS)
  type(advection_test_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                              !! call to register_advection_test_tracer.
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine advection_test_tracer_end

end module advection_test_tracer
