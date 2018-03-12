module pseudo_salt_tracer

! This file is part of MOM6. See LICENSE.md for the license.

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Andrew Shao, 2016                                               *
!*                                                                     *
!*    This file contains the routines necessary to model a passive     *
!*  tracer that uses the same boundary fluxes as salinity. At the      *
!*  beginning of the run, salt is set to the same as tv%S. Any         *
!*  deviations between this salt-like tracer and tv%S signifies a      *
!*  difference between how active and passive tracers are treated.     *
!*    A single subroutine is called from within each file to register  *
!*  each of the tracers for reinitialization and advection and to      *
!*  register the subroutine that initializes the tracers and set up    *
!*  their output and the subroutine that does any tracer physics or    *
!*  chemistry along with diapycnal mixing (included here because some  *
!*  tracers may float or swim vertically or dye diapycnal processes).  *
!*                                                                     *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, tr                                    *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_debugging,     only : hchksum
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl
use MOM_diag_to_Z, only : diag_to_Z_CS
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : file_exists, read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : query_initialized, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type, get_time
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init, only : tracer_Z_init
use MOM_variables, only : surface
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

use coupler_types_mod, only : coupler_type_set_data, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_pseudo_salt_tracer, initialize_pseudo_salt_tracer
public pseudo_salt_tracer_column_physics, pseudo_salt_tracer_surface_state
public pseudo_salt_stock, pseudo_salt_tracer_end

type, public :: pseudo_salt_tracer_CS ; private
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL()
  real, pointer :: ps(:,:,:) => NULL()   ! The array of pseudo-salt tracer used in this
                                         ! subroutine, in psu
  real, pointer :: diff(:,:,:) => NULL() ! The difference between the pseudo-salt
                                         ! tracer and the real salt, in psu.
  logical :: pseudo_salt_may_reinit = .true. ! Hard coding since this should not matter

  integer :: id_psd = -1

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()

  type(vardesc) :: tr_desc
end type pseudo_salt_tracer_CS

contains

function register_pseudo_salt_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(pseudo_salt_tracer_CS),  pointer    :: CS
  type(tracer_registry_type), pointer    :: tr_Reg
  type(MOM_restart_CS),       pointer    :: restart_CS
! This subroutine is used to register tracer fields and subroutines
! to be used with MOM.
! Arguments: HI - A horizontal index type structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in/out)  tr_Reg - A pointer that is set to point to the control structure
!                  for the tracer advection and diffusion module.
!  (in)      restart_CS - A pointer to the restart control structure.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "pseudo_salt_tracer" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  character(len=3)   :: name_tag ! String for creating identifying pseudo_salt
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

subroutine initialize_pseudo_salt_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                  sponge_CSp, diag_to_Z_CSp, tv)
  logical,                            intent(in) :: restart
  type(time_type), target,            intent(in) :: day
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h    !< Layer thicknesses, in H (usually m or kg m-2)
  type(diag_ctrl), target,            intent(in) :: diag
  type(ocean_OBC_type),               pointer    :: OBC
  type(pseudo_salt_tracer_CS),          pointer    :: CS
  type(sponge_CS),                    pointer    :: sponge_CSp
  type(diag_to_Z_CS),                 pointer    :: diag_to_Z_CSp
  type(thermo_var_ptrs),              intent(in) :: tv   !< A structure pointing to various thermodynamic variables
!   This subroutine initializes the tracer fields in CS%ps(:,:,:).

! Arguments: restart - .true. if the fields have already been read from
!                     a restart file.
!  (in)      day - Time of the start of the run.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in)      OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in/out)  CS - The control structure returned by a previous call to
!                 register_pseudo_salt_tracer.
!  (in/out)  sponge_CSp - A pointer to the control structure for the sponges, if
!                         they are in use.  Otherwise this may be unassociated.
!  (in/out)  diag_to_Z_Csp - A pointer to the control structure for diagnostics
!                            in depth space.
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

subroutine pseudo_salt_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, CS, tv, debug, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h_old, h_new, ea, eb
  type(forcing),                      intent(in) :: fluxes
  real,                               intent(in) :: dt   !< The amount of time covered by this call, in s
  type(pseudo_salt_tracer_CS),                pointer    :: CS
  type(thermo_var_ptrs),              intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  logical,                            intent(in) :: debug
  real,                             optional,intent(in)  :: evap_CFL_limit
  real,                             optional,intent(in)  :: minimum_forcing_depth

!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

! Arguments: h_old -  Layer thickness before entrainment, in m or kg m-2.
!  (in)      h_new -  Layer thickness after entrainment, in m or kg m-2.
!  (in)      ea - an array to which the amount of fluid entrained
!                 from the layer above during this call will be
!                 added, in m or kg m-2.
!  (in)      eb - an array to which the amount of fluid entrained
!                 from the layer below during this call will be
!                 added, in m or kg m-2.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_pseudo_salt_tracer.
!  (in)      tv - Thermodynamic structure with T and S
!  (in)      evap_CFL_limit - Limits how much water can be fluxed out of the top layer
!                             Stored previously in diabatic CS.
!  (in)      minimum_forcing_depth - The smallest depth over which fluxes can be applied
!                             Stored previously in diabatic CS.
!  (in)      debug - Calculates checksums
!
! The arguments to this subroutine are redundant in that
!     h_new[k] = h_old[k] + ea[k] - eb[k-1] + eb[k] - ea[k+1]

  real :: Isecs_per_year = 1.0 / (365.0*86400.0)
  real :: year, h_total, scale, htot, Ih_limit
  integer :: secs, days
  integer :: i, j, k, is, ie, js, je, nz, k_max
  real, allocatable :: local_tr(:,:,:)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  real, dimension(:,:), pointer :: net_salt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  net_salt=>fluxes%netSalt

  if (.not.associated(CS)) return
  if (.not.associated(CS%diff)) return

  if (debug) then
    call hchksum(tv%S,"salt pre pseudo-salt vertdiff", G%HI)
    call hchksum(CS%ps,"pseudo_salt pre pseudo-salt vertdiff", G%HI)
  endif

  ! This uses applyTracerBoundaryFluxesInOut, usually in ALE mode
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do k=1,nz ;do j=js,je ; do i=is,ie
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo;
    call applyTracerBoundaryFluxesInOut(G, GV, CS%ps, dt, fluxes, h_work, &
      evap_CFL_limit, minimum_forcing_depth, out_flux_optional=net_salt)
    call tracer_vertdiff(h_work, ea, eb, dt, CS%ps, G, GV)
  else
    call tracer_vertdiff(h_old, ea, eb, dt, CS%ps, G, GV)
  endif

  do k=1,nz ; do j=js,je ; do i=is,ie
    CS%diff(i,j,k) = CS%ps(i,j,k)-tv%S(i,j,k)
  enddo ; enddo ; enddo

  if(debug) then
    call hchksum(tv%S,"salt post pseudo-salt vertdiff", G%HI)
    call hchksum(CS%ps,"pseudo_salt post pseudo-salt vertdiff", G%HI)
  endif

  if (CS%id_psd>0) call post_data(CS%id_psd, CS%diff, CS%diag)

end subroutine pseudo_salt_tracer_column_physics

function pseudo_salt_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h    !< Layer thicknesses, in H (usually m or kg m-2)
  real, dimension(:),                 intent(out)   :: stocks
  type(pseudo_salt_tracer_CS),                pointer       :: CS
  character(len=*), dimension(:),     intent(out)   :: names
  character(len=*), dimension(:),     intent(out)   :: units
  integer, optional,                  intent(in)    :: stock_index
  integer                                           :: pseudo_salt_stock
! This function calculates the mass-weighted integral of all tracer stocks,
! returning the number of stocks it has calculated.  If the stock_index
! is present, only the stock corresponding to that coded index is returned.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (out)     stocks - the mass-weighted integrated amount of each tracer,
!                     in kg times concentration units.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_pseudo_salt_tracer.
!  (out)     names - the names of the stocks calculated.
!  (out)     units - the units of the stocks calculated.
!  (in,opt)  stock_index - the coded index of a specific stock being sought.
! Return value: the number of stocks calculated here.

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

  call query_vardesc(CS%tr_desc, name=names(1), units=units(1), caller="pseudo_salt_stock")
  units(1) = trim(units(1))//" kg"
  stocks(1) = 0.0
  do k=1,nz ; do j=js,je ; do i=is,ie
    stocks(1) = stocks(1) + CS%diff(i,j,k) * &
                         (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
  enddo ; enddo ; enddo
  stocks(1) = GV%H_to_kg_m2 * stocks(1)

  pseudo_salt_stock = 1

end function pseudo_salt_stock

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
!! This particular tracer package does not report anything back to the coupler.
subroutine pseudo_salt_tracer_surface_state(state, h, G, CS)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure.
  type(surface),          intent(inout) :: state !< A structure containing fields that
                                              !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                          intent(in)    :: h  !< Layer thickness, in m or kg m-2.
  type(pseudo_salt_tracer_CS), pointer  :: CS !< The control structure returned by a previous
                                              !! call to register_pseudo_salt_tracer.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  ! By design, this tracer package does not return any surface states.

end subroutine pseudo_salt_tracer_surface_state

subroutine pseudo_salt_tracer_end(CS)
  type(pseudo_salt_tracer_CS), pointer :: CS
  integer :: m

  if (associated(CS)) then
    if (associated(CS%ps)) deallocate(CS%ps)
    if (associated(CS%diff)) deallocate(CS%diff)
    deallocate(CS)
  endif
end subroutine pseudo_salt_tracer_end

end module pseudo_salt_tracer
