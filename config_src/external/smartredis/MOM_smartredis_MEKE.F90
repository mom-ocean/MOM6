!> Contains routines that contain dummy routines for the smart
module MOM_smartredis_meke

use MOM_diag_mediator,     only : diag_ctrl, time_type
use MOM_error_handler,     only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_grid,              only : ocean_grid_type
use MOM_file_parser,       only : param_file_type
use MOM_smartredis,        only : smartredis_CS_type
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs
use MOM_verticalGrid,      only : verticalGrid_type

implicit none; private

#include <MOM_memory.h>

public smartredis_meke_init, infer_meke

type, public :: smartredis_meke_CS_type; private

end type smartredis_meke_CS_type

contains

!> Initializer for the SmartRedis MEKE module that uses ML to predict eddy kinetic energy
subroutine smartredis_meke_init(diag, G, US, Time, param_file, smartredis_CS, CS)
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(ocean_grid_type),         intent(inout) :: G          !< The ocean's grid structure.
  type(unit_scale_type),         intent(in)    :: US         !< A dimensional unit scaling type
  type(time_type),               intent(in)    :: Time       !< The current model time.
  type(param_file_type),         intent(in)    :: param_file !< Parameter file parser structure.
  type(smartredis_CS_type), target,     intent(in)    :: smartredis_CS !< SmartRedis client
  type(smartredis_meke_CS_type), intent(inout) :: CS         !< Control structure for this module

  call MOM_error(FATAL,"smartredis_meke_init was compiled using the dummy module. Recompile"//&
                       "with source code from https://github.com/CrayLabs/MOM6-smartredis")
end subroutine smartredis_meke_init

!> Use the SmartRedis client to call a machine learning to predict eddy kinetic energy
subroutine infer_meke(G, GV, US, CS, Time, MEKE, Rd_dx_h, u, v, tv, h, dt)
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US         !< A dimensional unit scaling type
  type(time_type),                           intent(in)    :: Time       !< The current model time.
  type(smartredis_meke_CS_type),             intent(in)    :: CS !< Control structure for inferring MEKE
                                                                 !! using SmartRedis
  real, dimension(SZI_(G),SZJ_(G)), intent(  out) :: MEKE !< Vertically averaged eddy kinetic energy [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(in   ) :: Rd_dx_h !< Rossby radius of deformation over
                                                             !! the grid length scale [nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in) :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in) :: v  !< Meridional velocity [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Type containing thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt !< Model(baroclinic) time-step [T ~> s].

  call MOM_error(FATAL,"infer_meke was compiled using the dummy module. Recompile"//&
                       "with source code from https://github.com/CrayLabs/MOM6-smartredis")

end subroutine infer_meke

end module MOM_smartredis_meke
