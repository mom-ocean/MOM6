! This file is part of MOM6. See LICENSE.md for the license.
!> Controls where open boundary conditions are applied
module MOM_boundary_update

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator,         only : time_type
use MOM_domains,               only : pass_var, pass_vector
use MOM_domains,               only : To_All, SCALAR_PAIR, CGRID_NE
use MOM_error_handler,         only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser,           only : get_param, log_version, param_file_type, log_param
use MOM_grid,                  only : ocean_grid_type
use MOM_dyn_horgrid,           only : dyn_horgrid_type
use MOM_open_boundary,         only : ocean_obc_type, update_OBC_segment_data
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_tracer_registry,       only : add_tracer_OBC_values, tracer_registry_type
use MOM_variables,             only : thermo_var_ptrs
use tidal_bay_initialization,  only : tidal_bay_set_OBC_data

implicit none ; private

#include <MOM_memory.h>

public update_OBC_data

integer :: id_clock_pass

character(len=40)  :: mod = "MOM_boundary_update" ! This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"

contains

!> Calls appropriate routine to update the open boundary conditions.
subroutine update_OBC_data(OBC, G, GV, tv, h, Time)
  type(ocean_grid_type),          intent(in) :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !<  Ocean vertical grid structure
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h !< layer thickness
  type(ocean_OBC_type),           pointer    :: OBC !< Open boundary structure
  type(time_type),                intent(in) :: Time !< Model time
  ! Local variables
  logical :: read_OBC_eta = .false.
  logical :: read_OBC_uv = .false.
  logical :: read_OBC_TS = .false.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: isd_off, jsd_off
  integer :: IsdB, IedB, JsdB, JedB
  character(len=40)  :: mod = "update_OBC_data" ! This subroutine's name.
  character(len=200) :: filename, OBC_file, inputdir ! Strings for file/path

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (OBC%OBC_user_config == "tidal_bay") then
    call tidal_bay_set_OBC_data(OBC, G, h, Time)
  elseif (OBC%needs_IO_for_data) then
    call update_OBC_segment_data(G, GV, OBC, tv, h, Time)
  endif

end subroutine update_OBC_data

!> \namespace mom_boundary_update
!! This module updates the open boundary arrays when time-varying.
!! It caused a circular dependency with the tidal_bay setup when
!! MOM_open_boundary.
!!
!! A small fragment of the grid is shown below:
!!
!!    j+1  x ^ x ^ x   At x:  q, CoriolisBu
!!    j+1  > o > o >   At ^:  v, tauy
!!    j    x ^ x ^ x   At >:  u, taux
!!    j    > o > o >   At o:  h, bathyT, buoy, tr, T, S, Rml, ustar
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!!
!! The boundaries always run through q grid points (x).

end module MOM_boundary_update
