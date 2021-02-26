!> Dyed open boundary conditions
module dyed_obcs_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_dyn_horgrid,     only : dyn_horgrid_type
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,     only : get_param, log_version, param_file_type
use MOM_get_input,       only : directories
use MOM_grid,            only : ocean_grid_type
use MOM_open_boundary,   only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_open_boundary,   only : OBC_segment_type, register_segment_tracer
use MOM_tracer_registry, only : tracer_registry_type, tracer_name_lookup
use MOM_tracer_registry, only : tracer_type
use MOM_variables,       only : thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public dyed_obcs_set_OBC_data

integer :: ntr = 0 !< Number of dye tracers
                   !! \todo This is a module variable. Move this variable into the control structure.

contains

!> This subroutine sets the dye properties at open boundary conditions.
subroutine dyed_obcs_set_OBC_data(OBC, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                                                !! to parse for model parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.
  ! Local variables
  character(len=40)  :: mdl = "dyed_obcs_set_OBC_data" ! This subroutine's name.
  character(len=80)  :: name, longname
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, m, n, nz
  integer :: IsdB, IedB, JsdB, JedB
  real :: dye
  type(OBC_segment_type), pointer :: segment => NULL()
  type(tracer_type), pointer      :: tr_ptr => NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return

  call get_param(param_file, mdl, "NUM_DYE_TRACERS", ntr, &
                 "The number of dye tracers in this run. Each tracer "//&
                 "should have a separate boundary segment.", default=0,   &
                 do_not_log=.true.)

  if (OBC%number_of_segments < ntr) then
    call MOM_error(WARNING, "Error in dyed_obc segment setup")
    return   !!! Need a better error message here
  endif

! ! Set the inflow values of the dyes, one per segment.
! ! We know the order: north, south, east, west
  do m=1,ntr
    write(name,'("dye_",I2.2)') m
    write(longname,'("Concentration of dyed_obc Tracer ",I2.2, " on segment ",I2.2)') m, m
    call tracer_name_lookup(tr_Reg, tr_ptr, name)

    do n=1,OBC%number_of_segments
      if (n == m) then
        dye = 1.0
      else
        dye = 0.0
      endif
      call register_segment_tracer(tr_ptr, param_file, GV, &
                                   OBC%segment(n), OBC_scalar=dye)
    enddo
  enddo

end subroutine dyed_obcs_set_OBC_data

!> \namespace dyed_obcs_initialization
!!
!! Setting dyes, one for painting the inflow on each side.
end module dyed_obcs_initialization
