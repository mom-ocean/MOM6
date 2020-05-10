! This file is part of MOM6. See LICENSE.md for the license.
!> Controls where open boundary conditions are applied
module MOM_boundary_update

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator,         only : time_type
use MOM_error_handler,         only : MOM_mesg, MOM_error, FATAL, WARNING
use MOM_file_parser,           only : get_param, log_version, param_file_type, log_param
use MOM_grid,                  only : ocean_grid_type
use MOM_dyn_horgrid,           only : dyn_horgrid_type
use MOM_open_boundary,         only : ocean_obc_type, update_OBC_segment_data
use MOM_open_boundary,         only : OBC_registry_type, file_OBC_CS
use MOM_open_boundary,         only : register_file_OBC, file_OBC_end
use MOM_unit_scaling,          only : unit_scale_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_tracer_registry,       only : tracer_registry_type
use MOM_variables,             only : thermo_var_ptrs
use tidal_bay_initialization,  only : tidal_bay_set_OBC_data, register_tidal_bay_OBC
use tidal_bay_initialization,  only : tidal_bay_OBC_end, tidal_bay_OBC_CS
use Kelvin_initialization,     only : Kelvin_set_OBC_data, register_Kelvin_OBC
use Kelvin_initialization,     only : Kelvin_OBC_end, Kelvin_OBC_CS
use shelfwave_initialization,  only : shelfwave_set_OBC_data, register_shelfwave_OBC
use shelfwave_initialization,  only : shelfwave_OBC_end, shelfwave_OBC_CS
use dyed_channel_initialization, only : dyed_channel_update_flow, register_dyed_channel_OBC
use dyed_channel_initialization, only : dyed_channel_OBC_end, dyed_channel_OBC_CS

implicit none ; private

#include <MOM_memory.h>

public call_OBC_register, OBC_register_end
public update_OBC_data

!> The control structure for the MOM_boundary_update module
type, public :: update_OBC_CS ; private
  logical :: use_files = .false.        !< If true, use external files for the open boundary.
  logical :: use_Kelvin = .false.       !< If true, use the Kelvin wave open boundary.
  logical :: use_tidal_bay = .false.    !< If true, use the tidal_bay open boundary.
  logical :: use_shelfwave = .false.    !< If true, use the shelfwave open boundary.
  logical :: use_dyed_channel = .false. !< If true, use the dyed channel open boundary.
  !>@{ Pointers to the control structures for named OBC specifications
  type(file_OBC_CS), pointer :: file_OBC_CSp => NULL()
  type(Kelvin_OBC_CS), pointer :: Kelvin_OBC_CSp => NULL()
  type(tidal_bay_OBC_CS), pointer :: tidal_bay_OBC_CSp => NULL()
  type(shelfwave_OBC_CS), pointer :: shelfwave_OBC_CSp => NULL()
  type(dyed_channel_OBC_CS), pointer :: dyed_channel_OBC_CSp => NULL()
  !>@}
end type update_OBC_CS

integer :: id_clock_pass !< A CPU time clock ID

! character(len=40)  :: mdl = "MOM_boundary_update" ! This module's name.

contains

!> The following subroutines and associated definitions provide the
!! machinery to register and call the subroutines that initialize
!! open boundary conditions.
subroutine call_OBC_register(param_file, CS, OBC)
  type(param_file_type),     intent(in) :: param_file !< Parameter file to parse
  type(update_OBC_CS),       pointer    :: CS         !< Control structure for OBCs
  type(ocean_OBC_type),      pointer    :: OBC        !< Open boundary structure

  ! Local variables
  character(len=40)  :: mdl = "MOM_boundary_update" ! This module's name.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  if (associated(CS)) then
    call MOM_error(WARNING, "call_OBC_register called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "USE_FILE_OBC", CS%use_files, &
                 "If true, use external files for the open boundary.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_TIDAL_BAY_OBC", CS%use_tidal_bay, &
                 "If true, use the tidal_bay open boundary.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_KELVIN_WAVE_OBC", CS%use_Kelvin, &
                 "If true, use the Kelvin wave open boundary.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_SHELFWAVE_OBC", CS%use_shelfwave, &
                 "If true, use the shelfwave open boundary.", &
                 default=.false.)
  call get_param(param_file, mdl, "USE_DYED_CHANNEL_OBC", CS%use_dyed_channel, &
                 "If true, use the dyed channel open boundary.", &
                 default=.false.)

  if (CS%use_files) CS%use_files = &
    register_file_OBC(param_file, CS%file_OBC_CSp, &
               OBC%OBC_Reg)
  if (CS%use_tidal_bay) CS%use_tidal_bay = &
    register_tidal_bay_OBC(param_file, CS%tidal_bay_OBC_CSp, &
               OBC%OBC_Reg)
  if (CS%use_Kelvin) CS%use_Kelvin = &
    register_Kelvin_OBC(param_file, CS%Kelvin_OBC_CSp, &
               OBC%OBC_Reg)
  if (CS%use_shelfwave) CS%use_shelfwave = &
    register_shelfwave_OBC(param_file, CS%shelfwave_OBC_CSp, &
               OBC%OBC_Reg)
  if (CS%use_dyed_channel) CS%use_dyed_channel = &
    register_dyed_channel_OBC(param_file, CS%dyed_channel_OBC_CSp, &
               OBC%OBC_Reg)

end subroutine call_OBC_register

!> Calls appropriate routine to update the open boundary conditions.
subroutine update_OBC_data(OBC, G, GV, US, tv, h, CS, Time)
  type(ocean_grid_type),                    intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                  intent(in)    :: GV   !< Ocean vertical grid structure
  type(unit_scale_type),                    intent(in)    :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),                    intent(in)    :: tv   !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h    !< layer thicknesses [H ~> m or kg m-2]
  type(ocean_OBC_type),                     pointer       :: OBC  !< Open boundary structure
  type(update_OBC_CS),                      pointer       :: CS   !< Control structure for OBCs
  type(time_type),                          intent(in)    :: Time !< Model time

  ! Local variables
  logical :: read_OBC_eta = .false.
  logical :: read_OBC_uv = .false.
  logical :: read_OBC_TS = .false.
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: isd_off, jsd_off
  integer :: IsdB, IedB, JsdB, JedB
  character(len=40)  :: mdl = "update_OBC_data" ! This subroutine's name.
  character(len=200) :: filename, OBC_file, inputdir ! Strings for file/path

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

! Something here... with CS%file_OBC_CSp?
! if (CS%use_files) &
!     call update_OBC_segment_data(G, GV, OBC, tv, h, Time)
  if (CS%use_tidal_bay) &
      call tidal_bay_set_OBC_data(OBC, CS%tidal_bay_OBC_CSp, G, h, Time)
  if (CS%use_Kelvin)  &
      call Kelvin_set_OBC_data(OBC, CS%Kelvin_OBC_CSp, G, GV, US, h, Time)
  if (CS%use_shelfwave) &
      call shelfwave_set_OBC_data(OBC, CS%shelfwave_OBC_CSp, G, h, Time)
  if (CS%use_dyed_channel) &
      call dyed_channel_update_flow(OBC, CS%dyed_channel_OBC_CSp, G, Time)
  if (OBC%needs_IO_for_data)  &
      call update_OBC_segment_data(G, GV, US, OBC, tv, h, Time)

end subroutine update_OBC_data

!> Clean up the OBC registry.
subroutine OBC_register_end(CS)
  type(update_OBC_CS),       pointer    :: CS !< Control structure for OBCs

  if (CS%use_files) call file_OBC_end(CS%file_OBC_CSp)
  if (CS%use_tidal_bay) call tidal_bay_OBC_end(CS%tidal_bay_OBC_CSp)
  if (CS%use_Kelvin) call Kelvin_OBC_end(CS%Kelvin_OBC_CSp)

  if (associated(CS)) deallocate(CS)
end subroutine OBC_register_end

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
