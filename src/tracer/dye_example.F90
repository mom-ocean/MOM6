!> A tracer package for using dyes to diagnose regional flows.
module regional_dyes

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coupler_types,      only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator,      only : diag_ctrl
use MOM_error_handler,      only : MOM_error, FATAL, WARNING
use MOM_file_parser,        only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,       only : forcing
use MOM_grid,               only : ocean_grid_type
use MOM_hor_index,          only : hor_index_type
use MOM_io,                 only : vardesc, var_desc, query_vardesc
use MOM_open_boundary,      only : ocean_OBC_type
use MOM_restart,            only : query_initialized, MOM_restart_CS
use MOM_sponge,             only : set_up_sponge_field, sponge_CS
use MOM_time_manager,       only : time_type
use MOM_tracer_registry,    only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic,    only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init,      only : tracer_Z_init
use MOM_unit_scaling,       only : unit_scale_type
use MOM_variables,          only : surface
use MOM_verticalGrid,       only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_dye_tracer, initialize_dye_tracer
public dye_tracer_column_physics, dye_tracer_surface_state
public dye_stock, regional_dyes_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure for the regional dyes tracer package
type, public :: dye_tracer_CS ; private
  integer :: ntr    !< The number of tracers that are actually used.
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  real, allocatable, dimension(:) :: dye_source_minlon !< Minimum longitude of region dye will be injected.
  real, allocatable, dimension(:) :: dye_source_maxlon !< Maximum longitude of region dye will be injected.
  real, allocatable, dimension(:) :: dye_source_minlat !< Minimum latitude of region dye will be injected.
  real, allocatable, dimension(:) :: dye_source_maxlat !< Maximum latitude of region dye will be injected.
  real, allocatable, dimension(:) :: dye_source_mindepth !< Minimum depth of region dye will be injected [Z ~> m].
  real, allocatable, dimension(:) :: dye_source_maxdepth !< Maximum depth of region dye will be injected [Z ~> m].
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL() !< The array of tracers used in this subroutine, in g m-3?

  integer, allocatable, dimension(:) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux if it is used and the
                                               !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc), allocatable :: tr_desc(:) !< Descriptions and metadata for the tracers
  logical :: tracers_may_reinit = .false. !< If true the tracers may be initialized if not found in a restart file
end type dye_tracer_CS

contains

!> This subroutine is used to register tracer fields and subroutines
!! to be used with MOM.
function register_dye_tracer(HI, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(dye_tracer_CS),        pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                 !! structure for the tracer advection and diffusion module.
  type(MOM_restart_CS),       pointer    :: restart_CS !< A pointer to the restart control structure.

! Local variables
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "regional_dyes" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  character(len=48)  :: desc_name ! The variable's descriptor.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_dye_tracer
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_dye_tracer called with an "// &
                             "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NUM_DYE_TRACERS", CS%ntr, &
                 "The number of dye tracers in this run. Each tracer "//&
                 "should have a separate region.", default=0)
  allocate(CS%dye_source_minlon(CS%ntr), &
           CS%dye_source_maxlon(CS%ntr), &
           CS%dye_source_minlat(CS%ntr), &
           CS%dye_source_maxlat(CS%ntr), &
           CS%dye_source_mindepth(CS%ntr), &
           CS%dye_source_maxdepth(CS%ntr))
  allocate(CS%ind_tr(CS%ntr))
  allocate(CS%tr_desc(CS%ntr))

  CS%dye_source_minlon(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MINLON", CS%dye_source_minlon, &
                 "This is the starting longitude at which we start injecting dyes.", &
                 fail_if_missing=.true.)
  if (minval(CS%dye_source_minlon(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MINLON ")

  CS%dye_source_maxlon(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MAXLON", CS%dye_source_maxlon, &
                 "This is the ending longitude at which we finish injecting dyes.", &
                 fail_if_missing=.true.)
  if (minval(CS%dye_source_maxlon(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MAXLON ")

  CS%dye_source_minlat(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MINLAT", CS%dye_source_minlat, &
                 "This is the starting latitude at which we start injecting dyes.", &
                 fail_if_missing=.true.)
  if (minval(CS%dye_source_minlat(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MINLAT ")

  CS%dye_source_maxlat(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MAXLAT", CS%dye_source_maxlat, &
                 "This is the ending latitude at which we finish injecting dyes.", &
                 fail_if_missing=.true.)
  if (minval(CS%dye_source_maxlat(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MAXLAT ")

  CS%dye_source_mindepth(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MINDEPTH", CS%dye_source_mindepth, &
                 "This is the minimum depth at which we inject dyes.", &
                 units="m", scale=US%m_to_Z, fail_if_missing=.true.)
  if (minval(CS%dye_source_mindepth(:)) < -1.e29*US%m_to_Z) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MINDEPTH")

  CS%dye_source_maxdepth(:) = -1.e30
  call get_param(param_file, mdl, "DYE_SOURCE_MAXDEPTH", CS%dye_source_maxdepth, &
                 "This is the maximum depth at which we inject dyes.", &
                 units="m", scale=US%m_to_Z, fail_if_missing=.true.)
  if (minval(CS%dye_source_maxdepth(:)) < -1.e29*US%m_to_Z) &
    call MOM_error(FATAL, "register_dye_tracer: Not enough values provided for DYE_SOURCE_MAXDEPTH ")

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr)) ; CS%tr(:,:,:,:) = 0.0

  do m = 1, CS%ntr
    write(var_name(:),'(A,I3.3)') "dye",m
    write(desc_name(:),'(A,I3.3)') "Dye Tracer ",m
    CS%tr_desc(m) = var_desc(trim(var_name), "conc", trim(desc_name), caller=mdl)

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, &
                       caller="register_dye_tracer")
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                         tr_desc=CS%tr_desc(m), registry_diags=.true., &
                         restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_dye_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_dye_tracer = .true.
end function register_dye_tracer

!> This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
!! and it sets up the tracer output.
subroutine initialize_dye_tracer(restart, day, G, GV, h, diag, OBC, CS, sponge_CSp)
  logical,                            intent(in) :: restart !< .true. if the fields have already been
                                                            !! read from a restart file.
  type(time_type), target,            intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl), target,            intent(in) :: diag !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(dye_tracer_CS),                pointer    :: CS   !< The control structure returned by a previous
                                                         !! call to register_dye_tracer.
  type(sponge_CS),                    pointer    :: sponge_CSp    !< A pointer to the control structure
                                                                  !! for the sponges, if they are in use.

! Local variables
  character(len=24) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                ! years m3 s-1 or years kg s-1.
  logical :: OK
  integer :: i, j, k, m
  real    :: z_bot, z_center

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  CS%diag => diag

  ! Establish location of source
  do m= 1, CS%ntr
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      ! A dye is set dependent on the center of the cell being inside the rectangular box.
      if (CS%dye_source_minlon(m)<G%geoLonT(i,j) .and. &
          CS%dye_source_maxlon(m)>=G%geoLonT(i,j) .and. &
          CS%dye_source_minlat(m)<G%geoLatT(i,j) .and. &
          CS%dye_source_maxlat(m)>=G%geoLatT(i,j) .and. &
          G%mask2dT(i,j) > 0.0 ) then
        z_bot = -G%bathyT(i,j)
        do k = GV%ke, 1, -1
          z_center = z_bot + 0.5*h(i,j,k)*GV%H_to_Z
          if ( z_center > -CS%dye_source_maxdepth(m) .and. &
               z_center < -CS%dye_source_mindepth(m) ) then
            CS%tr(i,j,k,m) = 1.0
          endif
          z_bot = z_bot + h(i,j,k)*GV%H_to_Z
        enddo
      endif
    enddo ; enddo
  enddo

end subroutine initialize_dye_tracer

!> This subroutine applies diapycnal diffusion and any other column
!! tracer physics or chemistry to the tracers from this file.
!! This is a simple example of a set of advected passive tracers.
!! The arguments to this subroutine are redundant in that
!!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)
subroutine dye_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
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
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(dye_tracer_CS),     pointer    :: CS   !< The control structure returned by a previous
                                              !! call to register_dye_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified
  real :: sfc_val  ! The surface value for the tracers.
  real :: Isecs_per_year  ! The number of seconds in a year.
  real :: year            ! The time in years.
  integer :: secs, days   ! Integer components of the time type.
  integer :: i, j, k, is, ie, js, je, nz, m
  real    :: z_bot, z_center

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

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
    do m=1,CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

  do m=1,CS%ntr
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      ! A dye is set dependent on the center of the cell being inside the rectangular box.
      if (CS%dye_source_minlon(m)<G%geoLonT(i,j) .and. &
          CS%dye_source_maxlon(m)>=G%geoLonT(i,j) .and. &
          CS%dye_source_minlat(m)<G%geoLatT(i,j) .and. &
          CS%dye_source_maxlat(m)>=G%geoLatT(i,j) .and. &
          G%mask2dT(i,j) > 0.0 ) then
        z_bot = -G%bathyT(i,j)
        do k=nz,1,-1
          z_center = z_bot + 0.5*h_new(i,j,k)*GV%H_to_Z
          if ( z_center > -CS%dye_source_maxdepth(m) .and. &
               z_center < -CS%dye_source_mindepth(m) ) then
            CS%tr(i,j,k,m) = 1.0
          endif
          z_bot = z_bot + h_new(i,j,k)*GV%H_to_Z
        enddo
      endif
    enddo ; enddo
  enddo

end subroutine dye_tracer_column_physics

!> This function calculates the mass-weighted integral of all tracer stocks,
!! returning the number of stocks it has calculated.  If the stock_index
!! is present, only the stock corresponding to that coded index is returned.
function dye_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(:),                 intent(out)   :: stocks !< the mass-weighted integrated amount of
                                                            !! each tracer, in kg times concentration units [kg conc].
  type(dye_tracer_CS),                pointer       :: CS   !< The control structure returned by a
                                                            !! previous call to register_dye_tracer.
  character(len=*), dimension(:),     intent(out)   :: names !< the names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units !< the units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< the coded index of a specific stock
                                                                   !! being sought.
  integer                                           :: dye_stock   !< Return value: the number of stocks
                                                                   !! calculated here.

! Local variables
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  dye_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="dye_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      stocks(m) = stocks(m) + CS%tr(i,j,k,m) * &
                             (G%mask2dT(i,j) * G%US%L_to_m**2*G%areaT(i,j) * h(i,j,k))
    enddo ; enddo ; enddo
    stocks(m) = GV%H_to_kg_m2 * stocks(m)
  enddo
  dye_stock = CS%ntr

end function dye_stock

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
!! This particular tracer package does not report anything back to the coupler.
subroutine dye_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(dye_tracer_CS),     pointer       :: CS !< The control structure returned by a previous
                                               !! call to register_dye_tracer.

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
      call set_coupler_type_data(CS%tr(:,:,1,m), CS%ind_tr(m), sfc_state%tr_fields, &
                   idim=(/isd, is, ie, ied/), jdim=(/jsd, js, je, jed/) )
    enddo
  endif

end subroutine dye_tracer_surface_state

!> Clean up any allocated memory after the run.
subroutine regional_dyes_end(CS)
  type(dye_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                     !! call to register_dye_tracer.
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine regional_dyes_end

!> \namespace regional_dyes
!!
!!    This file contains an example of the code that is needed to set
!!  up and use a set (in this case two) of dynamically passive tracers
!!  for diagnostic purposes.  The tracers here are dye tracers which
!!  are set to 1 within the geographical region specified. The depth
!!  which a tracer is set is determined by calculating the depth from
!!  the seafloor upwards through the column.

end module regional_dyes
