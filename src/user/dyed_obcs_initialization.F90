module dyed_obcs_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_dyn_horgrid,     only : dyn_horgrid_type
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser,     only : get_param, log_version, param_file_type
use MOM_get_input,       only : directories
use MOM_grid,            only : ocean_grid_type
use MOM_io,              only : vardesc, var_desc
use MOM_open_boundary,   only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_open_boundary,   only : OBC_segment_type, register_segment_tracer
use MOM_tracer_registry, only : tracer_registry_type, add_tracer_OBC_values
use MOM_variables,       only : thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public dyed_obcs_set_OBC_data

integer, parameter :: NTR = 4

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
! Don't even need these any more for this problem...
! real, pointer, dimension(:,:,:) :: &
!   OBC_dye_1_v => NULL(), &    ! Specify the dye concentrations at the boundaries,
!   OBC_dye_2_v => NULL(), &    ! at both u and v points.
!   OBC_dye_3_u => NULL(), &
!   OBC_dye_4_u => NULL()

  character(len=40)  :: mdl = "dyed_obcs_set_OBC_data" ! This subroutine's name.
  character(len=80)  :: name, longname
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, m, nz
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment
  type(vardesc) :: tr_desc(NTR)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return

  if (OBC%number_of_segments .ne. 4) then
    print *, 'Error in dyed_obcs segment setup'
    return   !!! Need a better error message here
  endif

! allocate(OBC_dye_1_v(isd:ied,JsdB:JedB,nz)) ; OBC_dye_1_v(:,:,:) = 0.0
! allocate(OBC_dye_2_v(isd:ied,JsdB:JedB,nz)) ; OBC_dye_2_v(:,:,:) = 0.0
! allocate(OBC_dye_3_u(IsdB:IedB,jsd:jed,nz)) ; OBC_dye_3_u(:,:,:) = 0.0
! allocate(OBC_dye_4_u(IsdB:IedB,jsd:jed,nz)) ; OBC_dye_4_u(:,:,:) = 0.0

! ! Set the inflow values of the dyes, one per segment.
! ! We know the order: north, south, east, west
  do m=1,NTR
    write(name,'("dye_",I1.1)') m
    write(longname,'("Concentration of dyed_obc Tracer ",I1.1, " on segment ",I1.1)') m, m
    tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)

    call register_segment_tracer(tr_desc(m), param_file, OBC%segment(m)%HI, GV, &
                                 OBC%segment(m)%Reg, m, OBC_scalar=1.0)
  enddo
! call add_tracer_OBC_values("dye_1", tr_Reg, OBC_in_v=OBC_dye_1_v)
! call add_tracer_OBC_values("dye_2", tr_Reg, OBC_in_v=OBC_dye_2_v)
! call add_tracer_OBC_values("dye_3", tr_Reg, OBC_in_u=OBC_dye_3_u)
! call add_tracer_OBC_values("dye_4", tr_Reg, OBC_in_u=OBC_dye_4_u)

end subroutine dyed_obcs_set_OBC_data

!> \namespace dyed_obcs_initialization
!! Setting dyes, one for painting the inflow on each side.
end module dyed_obcs_initialization
