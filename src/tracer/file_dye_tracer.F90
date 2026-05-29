!> A tracer package for injecting dyes from file
module file_dye_tracer

use MOM_coms,            only : EFP_type
use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : MOM_error, FATAL
use MOM_file_parser,     only : get_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_hor_index,       only : hor_index_type
use MOM_io,              only : MOM_read_data, file_exists, query_vardesc
use MOM_io,              only : slasher, vardesc, var_desc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : MOM_restart_CS
use MOM_spatial_means,   only : global_mass_int_EFP
use MOM_sponge,          only : sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_diabatic, only : applyTracerBoundaryFluxesInOut, tracer_vertdiff
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_file_dye_tracer, initialize_file_dye_tracer
public file_dye_tracer_column_physics, file_dye_tracer_surface_state
public file_dye_stock, file_dyes_end

!> The control structure for the file_dye_tracer module
type, public :: file_dye_tracer_CS ; private
  !> The number of tracers that are actually used
  integer :: ntr

  !> Path to the tracer IC file
  character(len=200) :: tracer_file

  !> Input directory
  character(len=200) :: inputdir

  !> Descriptions and metadata for tracers
  type(vardesc), allocatable :: tr_desc(:)

  !> Array of tracers
  real, pointer :: tr(:,:,:,:) => NULL()

  !> Array of tracer source fields
  real, pointer :: tr_source(:,:,:,:) => NULL()

  !> If true, tracers may be initialized if not in a restart file
  logical :: tracers_may_reinit = .false.

  !> Tracer registry (pointer)
  type(tracer_registry_type), pointer :: tr_Reg => NULL()

  !> Restart control structure (pointer)
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()

  !> Diag control (pointer)
  type(diag_ctrl), pointer :: diag => NULL()

end type file_dye_tracer_CS

contains

!> This subroutine is used to register tracer fields and subroutines to be used with MOM.
function register_file_dye_tracer(HI, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type), intent(in) :: HI          !< The domain horizontal indexing
  type(verticalGrid_type), intent(in) :: GV       !< The ocean's vertical grid structure
  type(unit_scale_type), intent(in) :: US         !< A dimensional unit scaling type
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(file_dye_tracer_CS), pointer :: CS         !< A pointer to this module's control structure
  type(tracer_registry_type), pointer :: tr_Reg   !< A pointer to the tracer registry
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

  logical :: register_file_dye_tracer
  character(len=40) :: mdl = "file_dyes"
  character(len=48) :: var_name, desc_name
#include "version_variable.h"
  real, pointer :: tr_ptr(:,:,:) => NULL()
  integer :: m
  integer :: isd, ied, jsd, jed, nz
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  register_file_dye_tracer = .false.

  if (associated(CS)) then
    call MOM_error(FATAL, "register_file_dye_tracer called with an "// &
                   "associated control structure.")
  end if
  allocate(CS)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "FILE_DYE_TRACERS_FILE", CS%tracer_file, &
                 "The name of a file from which to read the initial conditions for the tracers", &
                 default="")
  ! expect a file, maybe make this fatal
  if (len_trim(CS%tracer_file) == 0) return

  call get_param(param_file, mdl, "NUM_FILE_DYE_TRACERS", CS%ntr, &
                 "The number of file-based dye tracers in this run.", &
                 default=0)

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".")
  CS%inputdir = slasher(CS%inputdir)

  allocate(CS%tr_desc(CS%ntr))

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr), source=0.0)
  allocate(CS%tr_source(isd:ied,jsd:jed,nz,CS%ntr), source=0.0)

  do m = 1, CS%ntr
    write(var_name(:), '(A,I3.3)') "dye", m
    write(desc_name(:), '(A,I3.3)') "Dye Tracer ", m
    CS%tr_desc(m) = var_desc(trim(var_name), "conc", trim(desc_name), caller=mdl)
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, &
                       caller="register_file_dye_tracer")
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                         tr_desc=CS%tr_desc(m), registry_diags=.true., &
                         restart_CS=restart_CS, mandatory=.not. CS%tracers_may_reinit)
  end do

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_file_dye_tracer = .true.
end function register_file_dye_tracer

!> This subroutine initializes the NTR tracer fields in tr(:,:,:,:)
!! and it sets up the tracer output.
subroutine initialize_file_dye_tracer(restart, day, G, GV, h, diag, OBC, CS, sponge_CSp)
  logical, intent(in) :: restart              !< fields already read from a restart file
  type(time_type), target, intent(in) :: day  !< Time of start of run
  type(ocean_grid_type), intent(in) :: G      !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, intent(in), &
    dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl), target, intent(in) :: diag !< Diagnostic control structure
  type(ocean_OBC_type), pointer :: OBC        !< Open boundary control structure
  type(file_dye_tracer_CS), pointer :: CS     !< This module's control structure
  type(sponge_CS), pointer :: sponge_CSp      !< Sponge control structure

  integer :: m
  character(len=32) :: name

  if (.not. associated(CS)) return
  if (CS%ntr < 1) return

  CS%diag => diag

  if (.not. file_exists(trim(CS%inputdir)//trim(CS%tracer_file), G%Domain)) then
    call MOM_error(FATAL, "initialize_file_dye_tracer: Unable to open " // &
                   CS%tracer_file)
  end if

  do m = 1, CS%ntr
    call query_vardesc(CS%tr_desc(m), name, caller="initialize_file_dye_tracer")
    call MOM_read_data(trim(CS%inputdir)//trim(CS%tracer_file), trim(name), CS%tr_source(:,:,:,m), G%Domain)
  end do

end subroutine initialize_file_dye_tracer

!> This subroutine applies diapycnal diffusion and any other column
!! tracer physics or chemistry to the tracers from this file.
!! This is a simple example of a set of advected passive tracers.
!! The arguments to this subroutine are redundant in that
!!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)
subroutine file_dye_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, &
                                          evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type), intent(in) :: G       !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV    !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
    intent(in) :: h_old                        !< Layer thickness before entrainment [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
    intent(in) :: h_new                        !< Layer thickness after entrainment [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
    intent(in) :: ea                           !< An array to which the amount of fluid entrained
                                               !! from the layer above during this call will be
                                               !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
    intent(in) :: eb                           !< An array to which the amount of fluid entrained
                                               !! from the layer below during this call will be
                                               !! added [H ~> m or kg m-2].
  type(forcing), intent(in) :: fluxes          !< A structure with pointers to forcing fields
  real, intent(in) :: dt                       !< Amount of time covered by this call [T ~> s]
  type(unit_scale_type), intent(in) :: US      !< A dimensional unit scaling type
  type(file_dye_tracer_CS), pointer :: CS      !< This module's control structure
  real, optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                               !! be fluxed out of the top layer in a timestep [nondim]
  real, optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                                      !! fluxes can be applied [H ~> m or kg m-2]

  integer :: m
  integer :: is, ie, js, je, nz
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) return
  if (CS%ntr < 1) return

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m = 1, CS%ntr
      h_work(is:ie,js:je,1:nz) = h_old(is:ie,js:je,1:nz)

      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    end do
  else
    do m = 1, CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    end do
  end if

  ! add from source every timestep
  do m = 1, CS%ntr
    CS%tr(is:ie,js:je,1:nz,m) = CS%tr(is:ie,js:je,1:nz,m) + CS%tr_source(is:ie,js:je,1:nz,m)
  end do
end subroutine file_dye_tracer_column_physics

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
subroutine file_dye_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type), intent(in) :: G     !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV  !< The ocean's vertical grid structure
  type(surface), intent(inout) :: sfc_state  !< A structure containing surface state fields
  real, intent(in), &
    dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h !< Layer thickness [H ~> m or kg m-2]
  type(file_dye_tracer_CS), pointer :: CS    !< This module's control structure

  if (.not. associated(CS)) return

  ! we would potentially couple the surface to the coupler here
end subroutine file_dye_tracer_surface_state

!> This function calculates the mass-weighted integral of all tracer stocks,
!! returning the number of stocks it has calculated.  If the stock_index
!! is present, only the stock corresponding to that coded index is returned.
function file_dye_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type), intent(in) :: G       !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV    !< The ocean's vertical grid structure
  real, intent(in), &
    dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h   !< Layer thickness [H ~> m or kg m-2]
  type(EFP_type), intent(out) :: stocks(:)     !< Mass-weighted integrated amount of each
                                               !! tracer, in kg times concentration units [kg conc].
  type(file_dye_tracer_CS), pointer :: CS      !< This module's control structure
  character(len=*), intent(out) :: names(:)    !< The names of the stocks calculated
  character(len=*), intent(out) :: units(:)    !< The units of the stocks calculated
  integer, optional, intent(in) :: stock_index !< The number of stocks calculated here

  integer :: file_dye_stock

  integer :: m

  file_dye_stock = 0
  if (.not. associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
      ! check whether this stock is available...
      return
    end if ; end if

  do m = 1, CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="file_dye_stock")
    units(m) = trim(units(m)) // " kg"
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%tr(:,:,:,m), on_PE_only=.true.)
  end do

  file_dye_stock = CS%ntr
end function file_dye_stock

!> Deallocate memory associated with this module
subroutine file_dyes_end(CS)
  type(file_dye_tracer_CS), pointer :: CS !< This module's control structure

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    if (associated(CS%tr_source)) deallocate(CS%tr_source)
    deallocate(CS)
  end if
end subroutine file_dyes_end

end module file_dye_tracer
