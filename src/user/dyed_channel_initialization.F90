module dyed_channel_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_dyn_horgrid,     only : dyn_horgrid_type
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
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

public dyed_channel_set_OBC_data

integer :: ntr = 0

contains

!> This subroutine sets the dye and flow properties at open boundary conditions.
subroutine dyed_channel_set_OBC_data(OBC, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                                                !! to parse for model parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.

! Local variables
  character(len=40)  :: mdl = "dyed_channel_set_OBC_data" ! This subroutine's name.
  character(len=80)  :: name, longname
  real :: zonal_flow
  integer :: i, j, k, l, itt, isd, ied, jsd, jed, m, n, nz
  integer :: IsdB, IedB, JsdB, JedB
  real :: dye
  type(OBC_segment_type), pointer :: segment
  type(vardesc), allocatable, dimension(:) :: tr_desc

  nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) call MOM_error(FATAL, 'dyed_channel_initialization.F90: '// &
        'dyed_channel_set_OBC_data() was called but OBC type was not initialized!')

  call get_param(param_file, mdl, "SUPERCRITICAL_ZONAL_FLOW", zonal_flow, &
                 "Constant zonal flow imposed at upstream open boundary.", &
                 units="m/s", default=8.57)

  do l=1, OBC%number_of_segments
    segment => OBC%segment(l)
    if (.not. segment%on_pe) cycle
    if (segment%gradient) cycle
    if (segment%oblique .and. .not. segment%nudged .and. .not. segment%Flather) cycle

    if (segment%is_E_or_W) then
      jsd = segment%HI%jsd ; jed = segment%HI%jed
      IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
      do k=1,G%ke
        do j=jsd,jed ; do I=IsdB,IedB
          if (segment%specified .or. segment%nudged) then
            segment%normal_vel(I,j,k) = zonal_flow
          endif
          if (segment%specified) then
            segment%normal_trans(I,j,k) = zonal_flow * G%dyCu(I,j)
          endif
        enddo ; enddo
      enddo
      do j=jsd,jed ; do I=IsdB,IedB
        segment%normal_vel_bt(I,j) = zonal_flow
      enddo ; enddo
    else
      isd = segment%HI%isd ; ied = segment%HI%ied
      JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
      do J=JsdB,JedB ; do i=isd,ied
        segment%normal_vel_bt(i,J) = 0.0
      enddo ; enddo
    endif
  enddo

  call get_param(param_file, mdl, "NUM_DYE_TRACERS", ntr, &
                 "The number of dye tracers in this run. Each tracer \n"//&
                 "should have a separate boundary segment.", default=0,   &
                 do_not_log=.true.)

  if (OBC%number_of_segments .lt. ntr) then
    call MOM_error(WARNING, "Error in dyed_obc segment setup")
    return   !!! Need a better error message here
  endif
  allocate(tr_desc(ntr))

! ! Set the inflow values of the dyes, one per segment.
! ! We know the order: north, south, east, west
  do m=1,ntr
    write(name,'("dye_",I1.1)') m
    write(longname,'("Concentration of dyed_obc Tracer ",I1.1, " on segment ",I1.1)') m, m
    tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)

    do n=1,OBC%number_of_segments
      if (n == m) then
        dye = 1.0
      else
        dye = 0.0
      endif
      call register_segment_tracer(tr_desc(m), param_file, GV, &
                                   OBC%segment(n), OBC_scalar=dye)
    enddo
  enddo
  deallocate(tr_desc)

end subroutine dyed_channel_set_OBC_data

!> \namespace dyed_channel_initialization
!! Setting dyes, one for painting the inflow on each side.
end module dyed_channel_initialization
