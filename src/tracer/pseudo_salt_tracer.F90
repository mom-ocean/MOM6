!> A tracer package that mimics salinity
module pseudo_salt_tracer

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : query_initialized, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init, only : tracer_Z_init
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface, thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_pseudo_salt_tracer, initialize_pseudo_salt_tracer
public pseudo_salt_tracer_column_physics, pseudo_salt_tracer_surface_state
public pseudo_salt_stock, pseudo_salt_tracer_end

!> The control structure for the pseudo-salt tracer
type, public :: pseudo_salt_tracer_CS ; private
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM tracer registry
  real, pointer :: ps(:,:,:) => NULL()   !< The array of pseudo-salt tracer used in this
                                         !! subroutine [ppt}
  real, pointer :: diff(:,:,:) => NULL() !< The difference between the pseudo-salt
                                         !! tracer and the real salt [ppt].
  logical :: pseudo_salt_may_reinit = .true. !< Hard coding since this should not matter

  integer :: id_psd = -1   !< A diagnostic ID

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc) :: tr_desc !< A description and metadata for the pseudo-salt tracer
end type pseudo_salt_tracer_CS

contains

!> Register the pseudo-salt tracer with MOM6
function register_pseudo_salt_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(pseudo_salt_tracer_CS),  pointer  :: CS !< The control structure returned by a previous
                                               !! call to register_pseudo_salt_tracer.
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                  !! structure for the tracer advection and
                                                  !! diffusion module
  type(MOM_restart_CS),       pointer    :: restart_CS !< A pointer to the restart control structure
! This subroutine is used to register tracer fields and subroutines
! to be used with MOM.

  ! Local variables
  character(len=40)  :: mdl = "pseudo_salt_tracer" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  character(len=3)   :: name_tag ! String for creating identifying pseudo_salt
! This include declares and sets the variable "version".
#include "version_variable.h"
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_pseudo_salt_tracer
  integer :: isd, ied, jsd, jed, nz, i, j
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_pseudo_salt_tracer called with an "// &
                             "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  allocate(CS%ps(isd:ied,jsd:jed,nz)) ; CS%ps(:,:,:) = 0.0
  allocate(CS%diff(isd:ied,jsd:jed,nz)) ; CS%diff(:,:,:) = 0.0

  CS%tr_desc = var_desc(trim("pseudo_salt"), "psu", &
                     "Pseudo salt passive tracer", caller=mdl)

  tr_ptr => CS%ps(:,:,:)
  call query_vardesc(CS%tr_desc, name=var_name, caller="register_pseudo_salt_tracer")
  ! Register the tracer for horizontal advection, diffusion, and restarts.
  call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, name="pseudo_salt", &
                       longname="Pseudo salt passive tracer", units="psu", &
                       registry_diags=.true., restart_CS=restart_CS, &
                       mandatory=.not.CS%pseudo_salt_may_reinit)

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_pseudo_salt_tracer = .true.

end function register_pseudo_salt_tracer

!> Initialize the pseudo-salt tracer
subroutine initialize_pseudo_salt_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                  sponge_CSp, tv)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(time_type),            target, intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(pseudo_salt_tracer_CS),        pointer    :: CS !< The control structure returned by a previous
                                                       !! call to register_pseudo_salt_tracer.
  type(sponge_CS),                    pointer    :: sponge_CSp !< Pointer to the control structure for the sponges.
  type(thermo_var_ptrs),              intent(in) :: tv   !< A structure pointing to various thermodynamic variables
!   This subroutine initializes the tracer fields in CS%ps(:,:,:).

  ! Local variables
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                ! years m3 s-1 or years kg s-1.
  logical :: OK
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB

  if (.not.associated(CS)) return
  if (.not.associated(CS%diff)) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%Time => day
  CS%diag => diag
  name = "pseudo_salt"

  call query_vardesc(CS%tr_desc, name=name, caller="initialize_pseudo_salt_tracer")
  if ((.not.restart) .or. (.not.query_initialized(CS%ps, name, CS%restart_CSp))) then
    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      CS%ps(i,j,k) = tv%S(i,j,k)
    enddo ; enddo ; enddo
  endif

  if (associated(OBC)) then
  ! Steal from updated DOME in the fullness of time.
  endif

  CS%id_psd = register_diag_field("ocean_model", "pseudo_salt_diff", CS%diag%axesTL, &
        day, "Difference between pseudo salt passive tracer and salt tracer", "psu")

end subroutine initialize_pseudo_salt_tracer

!> Apply sources, sinks and diapycnal diffusion to the tracers in this package.
subroutine pseudo_salt_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, tv, debug, &
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
  type(pseudo_salt_tracer_CS), pointer :: CS  !< The control structure returned by a previous
                                              !! call to register_pseudo_salt_tracer.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  logical,                 intent(in) :: debug !< If true calculate checksums
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables
  real :: year, h_total, scale, htot, Ih_limit
  integer :: secs, days
  integer :: i, j, k, is, ie, js, je, nz, k_max
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (.not.associated(CS%diff)) return

  if (debug) then
    call hchksum(tv%S,"salt pre pseudo-salt vertdiff", G%HI)
    call hchksum(CS%ps,"pseudo_salt pre pseudo-salt vertdiff", G%HI)
  endif

  ! This uses applyTracerBoundaryFluxesInOut, usually in ALE mode
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo
    call applyTracerBoundaryFluxesInOut(G, GV, CS%ps, dt, fluxes, h_work, &
                                        evap_CFL_limit, minimum_forcing_depth, out_flux_optional=fluxes%netSalt)
    call tracer_vertdiff(h_work, ea, eb, dt, CS%ps, G, GV)
  else
    call tracer_vertdiff(h_old, ea, eb, dt, CS%ps, G, GV)
  endif

  do k=1,nz ; do j=js,je ; do i=is,ie
    CS%diff(i,j,k) = CS%ps(i,j,k)-tv%S(i,j,k)
  enddo ; enddo ; enddo

  if (debug) then
    call hchksum(tv%S,"salt post pseudo-salt vertdiff", G%HI)
    call hchksum(CS%ps,"pseudo_salt post pseudo-salt vertdiff", G%HI)
  endif

  if (CS%id_psd>0) call post_data(CS%id_psd, CS%diff, CS%diag)

end subroutine pseudo_salt_tracer_column_physics


!> Calculates the mass-weighted integral of all tracer stocks, returning the number of stocks it has
!! calculated.  If the stock_index is present, only the stock corresponding to that coded index is returned.
function pseudo_salt_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(:),                 intent(out)   :: stocks !< the mass-weighted integrated amount of each
                                                              !! tracer, in kg times concentration units [kg conc].
  type(pseudo_salt_tracer_CS),        pointer       :: CS !< The control structure returned by a previous
                                                          !! call to register_pseudo_salt_tracer.
  character(len=*), dimension(:),     intent(out)   :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units  !< The units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< The coded index of a specific stock
                                                              !! being sought.
  integer                                           :: pseudo_salt_stock !< Return value: the number of
                                                              !! stocks calculated here.

  ! Local variables
  real :: stock_scale ! The dimensional scaling factor to convert stocks to kg [kg H-1 L-2 ~> kg m-3 or nondim]
  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  pseudo_salt_stock = 0
  if (.not.associated(CS)) return
  if (.not.associated(CS%diff)) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  stock_scale = G%US%L_to_m**2 * GV%H_to_kg_m2
  call query_vardesc(CS%tr_desc, name=names(1), units=units(1), caller="pseudo_salt_stock")
  units(1) = trim(units(1))//" kg"
  stocks(1) = 0.0
  do k=1,nz ; do j=js,je ; do i=is,ie
    stocks(1) = stocks(1) + CS%diff(i,j,k) * (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
  enddo ; enddo ; enddo
  stocks(1) = stock_scale * stocks(1)

  pseudo_salt_stock = 1

end function pseudo_salt_stock

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
!! This particular tracer package does not report anything back to the coupler.
subroutine pseudo_salt_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(pseudo_salt_tracer_CS),  pointer  :: CS !< The control structure returned by a previous
                                               !! call to register_pseudo_salt_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  ! By design, this tracer package does not return any surface states.

end subroutine pseudo_salt_tracer_surface_state

!> Deallocate memory associated with this tracer package
subroutine pseudo_salt_tracer_end(CS)
  type(pseudo_salt_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                              !! call to register_pseudo_salt_tracer.
  integer :: m

  if (associated(CS)) then
    if (associated(CS%ps)) deallocate(CS%ps)
    if (associated(CS%diff)) deallocate(CS%diff)
    deallocate(CS)
  endif
end subroutine pseudo_salt_tracer_end

!> \namespace pseudo_salt_tracer
!!
!!  By Andrew Shao, 2016
!!
!!    This file contains the routines necessary to model a passive
!!  tracer that uses the same boundary fluxes as salinity. At the
!!  beginning of the run, salt is set to the same as tv%S. Any
!!  deviations between this salt-like tracer and tv%S signifies a
!!  difference between how active and passive tracers are treated.

end module pseudo_salt_tracer
