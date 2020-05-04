!> A tracer package of ideal age tracers
module ideal_age_example

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : query_initialized, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type, time_type_to_real
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init, only : tracer_Z_init
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_verticalGrid, only : verticalGrid_type

use coupler_types_mod, only : coupler_type_set_data, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_ideal_age_tracer, initialize_ideal_age_tracer
public ideal_age_tracer_column_physics, ideal_age_tracer_surface_state
public ideal_age_stock, ideal_age_example_end

integer, parameter :: NTR_MAX = 3 !< the maximum number of tracers in this module.

!> The control structure for the ideal_age_tracer package
type, public :: ideal_age_tracer_CS ; private
  integer :: ntr    !< The number of tracers that are actually used.
  logical :: coupled_tracers = .false. !< These tracers are not offered to the coupler.
  integer :: nkml   !< The number of layers in the mixed layer.  The ideal
                    !1 age tracers are reset in the top nkml layers.
  character(len=200) :: IC_file !< The file in which the age-tracer initial values
                    !! can be found, or an empty string for internal initialization.
  logical :: Z_IC_file !< If true, the IC_file is in Z-space.  The default is false.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this package, in g m-3?
  real, dimension(NTR_MAX) :: IC_val = 0.0    !< The (uniform) initial condition value.
  real, dimension(NTR_MAX) :: young_val = 0.0 !< The value assigned to tr at the surface.
  real, dimension(NTR_MAX) :: land_val = -1.0 !< The value of tr used where land is masked out.
  real, dimension(NTR_MAX) :: sfc_growth_rate !< The exponential growth rate for the surface value [year-1].
  real, dimension(NTR_MAX) :: tracer_start_year !< The year in which tracers start aging, or at which the
                                              !! surface value equals young_val, in years.
  logical :: tracers_may_reinit  !< If true, these tracers be set up via the initialization code if
                                 !! they are not found in the restart files.
  logical :: tracer_ages(NTR_MAX) !< Indicates whether each tracer ages.

  integer, dimension(NTR_MAX) :: ind_tr !< Indices returned by aof_set_coupler_flux if it is used and the
                                        !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart controls structure

  type(vardesc) :: tr_desc(NTR_MAX) !< Descriptions and metadata for the tracers
end type ideal_age_tracer_CS

contains

!> Register the ideal age tracer fields to be used with MOM.
function register_ideal_age_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(ideal_age_tracer_CS),  pointer    :: CS !< The control structure returned by a previous
                                               !! call to register_ideal_age_tracer.
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                  !! structure for the tracer advection and
                                                  !! diffusion module
  type(MOM_restart_CS),       pointer    :: restart_CS !< A pointer to the restart control structure

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "ideal_age_example" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_ideal_age_tracer
  logical :: do_ideal_age, do_vintage, do_ideal_age_dated
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_ideal_age_tracer called with an "// &
                             "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DO_IDEAL_AGE", do_ideal_age, &
                 "If true, use an ideal age tracer that is set to 0 age "//&
                 "in the mixed layer and ages at unit rate in the interior.", &
                 default=.true.)
  call get_param(param_file, mdl, "DO_IDEAL_VINTAGE", do_vintage, &
                 "If true, use an ideal vintage tracer that is set to an "//&
                 "exponentially increasing value in the mixed layer and "//&
                 "is conserved thereafter.", default=.false.)
  call get_param(param_file, mdl, "DO_IDEAL_AGE_DATED", do_ideal_age_dated, &
                 "If true, use an ideal age tracer that is everywhere 0 "//&
                 "before IDEAL_AGE_DATED_START_YEAR, but the behaves like "//&
                 "the standard ideal age tracer - i.e. is set to 0 age in "//&
                 "the mixed layer and ages at unit rate in the interior.", &
                 default=.false.)


  call get_param(param_file, mdl, "AGE_IC_FILE", CS%IC_file, &
                 "The file in which the age-tracer initial values can be "//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mdl, "INPUTDIR/AGE_IC_FILE", CS%IC_file)
  endif
  call get_param(param_file, mdl, "AGE_IC_FILE_IS_Z", CS%Z_IC_file, &
                 "If true, AGE_IC_FILE is in depth space, not layer space", &
                 default=.false.)
  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code "//&
                 "if they are not found in the restart files.  Otherwise "//&
                 "it is a fatal error if the tracers are not found in the "//&
                 "restart files of a restarted run.", default=.false.)

  CS%ntr = 0
  if (do_ideal_age) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("age", "yr", "Ideal Age Tracer", cmor_field_name="agessc", caller=mdl)
    CS%tracer_ages(m) = .true. ; CS%sfc_growth_rate(m) = 0.0
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = 0.0
  endif

  if (do_vintage) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("vintage", "yr", "Exponential Vintage Tracer", &
                            caller=mdl)
    CS%tracer_ages(m) = .false. ; CS%sfc_growth_rate(m) = 1.0/30.0
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 1e-20 ; CS%tracer_start_year(m) = 0.0
    call get_param(param_file, mdl, "IDEAL_VINTAGE_START_YEAR", CS%tracer_start_year(m), &
                 "The date at which the ideal vintage tracer starts.", &
                 units="years", default=0.0)
  endif

  if (do_ideal_age_dated) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("age_dated","yr","Ideal Age Tracer with a Start Date",&
                            caller=mdl)
    CS%tracer_ages(m) = .true. ; CS%sfc_growth_rate(m) = 0.0
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = 0.0
    call get_param(param_file, mdl, "IDEAL_AGE_DATED_START_YEAR", CS%tracer_start_year(m), &
                 "The date at which the dated ideal age tracer starts.", &
                 units="years", default=0.0)
  endif

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr)) ; CS%tr(:,:,:,:) = 0.0

  do m=1,CS%ntr
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, &
                       caller="register_ideal_age_tracer")
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, tr_desc=CS%tr_desc(m), &
                         registry_diags=.true., restart_CS=restart_CS, &
                         mandatory=.not.CS%tracers_may_reinit, &
                         flux_scale=GV%H_to_m)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_ideal_age_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_ideal_age_tracer = .true.
end function register_ideal_age_tracer

!> Sets the ideal age traces to their initial values and sets up the tracer output
subroutine initialize_ideal_age_tracer(restart, day, G, GV, US, h, diag, OBC, CS, &
                                       sponge_CSp)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(time_type),            target, intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(ideal_age_tracer_CS),          pointer    :: CS !< The control structure returned by a previous
                                                       !! call to register_ideal_age_tracer.
  type(sponge_CS),                    pointer    :: sponge_CSp !< Pointer to the control structure for the sponges.

!   This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.

  ! Local variables
  character(len=24) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                ! years m3 s-1 or years kg s-1.
  character(len=72) :: cmorname ! The CMOR name of that variable.
  logical :: OK
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%Time => day
  CS%diag => diag
  CS%nkml = max(GV%nkml,1)

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=name, &
                       caller="initialize_ideal_age_tracer")
    if ((.not.restart) .or. (CS%tracers_may_reinit .and. .not. &
        query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then

      if (len_trim(CS%IC_file) > 0) then
  !  Read the tracer concentrations from a netcdf file.
        if (.not.file_exists(CS%IC_file, G%Domain)) &
          call MOM_error(FATAL, "initialize_ideal_age_tracer: "// &
                                 "Unable to open "//CS%IC_file)

        if (CS%Z_IC_file) then
          OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, name,&
                             G, US, -1e34, 0.0) ! CS%land_val(m))
          if (.not.OK) then
            OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, &
                     trim(name), G, US, -1e34, 0.0) ! CS%land_val(m))
            if (.not.OK) call MOM_error(FATAL,"initialize_ideal_age_tracer: "//&
                    "Unable to read "//trim(name)//" from "//&
                    trim(CS%IC_file)//".")
          endif
        else
          call MOM_read_data(CS%IC_file, trim(name), CS%tr(:,:,:,m), G%Domain)
        endif
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          if (G%mask2dT(i,j) < 0.5) then
            CS%tr(i,j,k,m) = CS%land_val(m)
          else
            CS%tr(i,j,k,m) = CS%IC_val(m)
          endif
        enddo ; enddo ; enddo
      endif

    endif ! restart
  enddo ! Tracer loop

  if (associated(OBC)) then
  ! Steal from updated DOME in the fullness of time.
  endif

end subroutine initialize_ideal_age_tracer

!> Applies diapycnal diffusion, aging and regeneration at the surface to the ideal age tracers
subroutine ideal_age_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, &
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
  type(ideal_age_tracer_CS), pointer  :: CS   !< The control structure returned by a previous
                                              !! call to register_ideal_age_tracer.
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
  real :: sfc_val  ! The surface value for the tracers.
  real :: Isecs_per_year  ! The inverse of the amount of time in a year [T-1 ~> s-1]
  real :: year            ! The time in years.
  integer :: i, j, k, is, ie, js, je, nz, m
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

  Isecs_per_year = 1.0 / (365.0*86400.0*US%s_to_T)
  !   Set the surface value of tracer 1 to increase exponentially
  ! with a 30 year time scale.
  year = US%s_to_T*time_type_to_real(CS%Time) * Isecs_per_year

  do m=1,CS%ntr
    if (CS%sfc_growth_rate(m) == 0.0) then
      sfc_val = CS%young_val(m)
    else
      sfc_val = CS%young_val(m) * &
          exp((year-CS%tracer_start_year(m)) * CS%sfc_growth_rate(m))
    endif
    do k=1,CS%nkml ; do j=js,je ; do i=is,ie
      if (G%mask2dT(i,j) > 0.5) then
        CS%tr(i,j,k,m) = sfc_val
      else
        CS%tr(i,j,k,m) = CS%land_val(m)
      endif
    enddo ; enddo ; enddo
  enddo
  do m=1,CS%ntr ; if (CS%tracer_ages(m) .and. &
                      (year>=CS%tracer_start_year(m))) then
!$OMP parallel do default(none) shared(is,ie,js,je,CS,nz,G,dt,Isecs_per_year,m)
    do k=CS%nkml+1,nz ; do j=js,je ; do i=is,ie
      CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + G%mask2dT(i,j)*dt*Isecs_per_year
    enddo ; enddo ; enddo
  endif ; enddo

end subroutine ideal_age_tracer_column_physics

!> Calculates the mass-weighted integral of all tracer stocks, returning the number of stocks it
!! has calculated.  If stock_index is present, only the stock corresponding to that coded index is found.
function ideal_age_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                      intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(:),                 intent(out)   :: stocks !< the mass-weighted integrated amount of each
                                                            !! tracer, in kg times concentration units [kg conc].
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  type(ideal_age_tracer_CS),          pointer       :: CS   !< The control structure returned by a previous
                                                            !! call to register_ideal_age_tracer.
  character(len=*), dimension(:),     intent(out)   :: names  !< the names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units  !< the units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< the coded index of a specific stock
                                                                   !! being sought.
  integer                                           :: ideal_age_stock !< The number of stocks calculated here.
! This function calculates the mass-weighted integral of all tracer stocks,
! returning the number of stocks it has calculated.  If the stock_index
! is present, only the stock corresponding to that coded index is returned.

  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  ideal_age_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="ideal_age_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      stocks(m) = stocks(m) + CS%tr(i,j,k,m) * &
                             (G%mask2dT(i,j) * G%US%L_to_m**2*G%areaT(i,j) * h(i,j,k))
    enddo ; enddo ; enddo
    stocks(m) = GV%H_to_kg_m2 * stocks(m)
  enddo
  ideal_age_stock = CS%ntr

end function ideal_age_stock

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
!! This particular tracer package does not report anything back to the coupler.
subroutine ideal_age_tracer_surface_state(sfc_state, h, G, CS)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure.
  type(surface),          intent(inout) :: sfc_state !< A structure containing fields that
                                              !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                          intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(ideal_age_tracer_CS), pointer    :: CS !< The control structure returned by a previous
                                              !! call to register_ideal_age_tracer.

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
                   sfc_state%tr_fields, idim=(/isd, is, ie, ied/), &
                   jdim=(/jsd, js, je, jed/) )
    enddo
  endif

end subroutine ideal_age_tracer_surface_state

!> Deallocate any memory associated with this tracer package
subroutine ideal_age_example_end(CS)
  type(ideal_age_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                           !! call to register_ideal_age_tracer.

  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine ideal_age_example_end

!> \namespace ideal_age_example
!!
!!  Originally by Robert Hallberg, 2002
!!
!!    This file contains an example of the code that is needed to set
!!  up and use a set (in this case two) of dynamically passive tracers
!!  for diagnostic purposes.  The tracers here are an ideal age tracer
!!  that ages at a rate of 1/year once it is isolated from the surface,
!!  and a vintage tracer, whose surface concentration grows exponen-
!!  with time with a 30-year timescale (similar to CFCs).

end module ideal_age_example
