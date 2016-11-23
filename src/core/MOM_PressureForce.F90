!> A thin wrapper for Boussinesq/non-Boussinesq forms of the pressure force calculation.
module MOM_PressureForce

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_PressureForce_AFV, only : PressureForce_AFV_Bouss, PressureForce_AFV_nonBouss
use MOM_PressureForce_AFV, only : PressureForce_AFV_init, PressureForce_AFV_CS
use MOM_PressureForce_Mont, only : PressureForce_Mont_Bouss, PressureForce_Mont_nonBouss
use MOM_PressureForce_Mont, only : PressureForce_Mont_init, PressureForce_Mont_CS
use MOM_tidal_forcing, only : calc_tidal_forcing, tidal_forcing_CS
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_ALE, only: ALE_CS
implicit none ; private

#include <MOM_memory.h>

public PressureForce, PressureForce_init, PressureForce_end

! Pressure force control structure
type, public :: PressureForce_CS ; private
  logical :: Analytic_FV_PGF !< If true, use the analytic finite volume form
                             !! (Adcroft et al., Ocean Mod. 2008) of the PGF.
  !> Control structure for the analytically integrated finite volume pressure force
  type(PressureForce_AFV_CS), pointer :: PressureForce_AFV_CSp => NULL()
  !> Control structure for the Montgomery potential form of pressure force
  type(PressureForce_Mont_CS), pointer :: PressureForce_Mont_CSp => NULL()
end type PressureForce_CS

contains

!> A thin layer between the model and the Boussinesq and non-Boussinesq pressure force routines.
subroutine PressureForce(h, tv, PFu, PFv, G, GV, CS, ALE_CSp, p_atm, pbce, eta)
  type(ocean_grid_type),                     intent(in)  :: G
  type(verticalGrid_type),                   intent(in)  :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)  :: h
  type(thermo_var_ptrs),                     intent(in)  :: tv
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out) :: PFu
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out) :: PFv
  type(PressureForce_CS),                    pointer     :: CS
  type(ALE_CS),                              pointer     :: ALE_CSp
  real, dimension(:,:),                     optional, pointer     :: p_atm
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), optional, intent(out) :: pbce
  real, dimension(SZI_(G),SZJ_(G)),         optional, intent(out) :: eta


  if (CS%Analytic_FV_PGF) then
    if (GV%Boussinesq) then
      call PressureForce_AFV_Bouss(h, tv, PFu, PFv, G, GV, CS%PressureForce_AFV_CSp, &
                                   ALE_CSp, p_atm, pbce, eta)
    else
      call PressureForce_AFV_nonBouss(h, tv, PFu, PFv, G, GV, CS%PressureForce_AFV_CSp, &
                                      p_atm, pbce, eta)
    endif
  else
    if (GV%Boussinesq) then
      call PressureForce_Mont_Bouss(h, tv, PFu, PFv, G, GV, CS%PressureForce_Mont_CSp, &
                                    p_atm, pbce, eta)
    else
      call PressureForce_Mont_nonBouss(h, tv, PFu, PFv, G, GV, CS%PressureForce_Mont_CSp, &
                                       p_atm, pbce, eta)
    endif
  endif

end subroutine Pressureforce

!> Initialize the pressure force control structure
subroutine PressureForce_init(Time, G, GV, param_file, diag, CS, tides_CSp)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< Vertical grid structure
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_CS),  pointer       :: CS   !< Pressure force control structure
  type(tidal_forcing_CS), optional, pointer :: tides_CSp !< Tide control structure
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_PressureForce" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "PressureForce_init called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "ANALYTIC_FV_PGF", CS%Analytic_FV_PGF, &
                 "If true the pressure gradient forces are calculated \n"//&
                 "with a finite volume form that analytically integrates \n"//&
                 "the equations of state in pressure to avoid any \n"//&
                 "possibility of numerical thermobaric instability, as \n"//&
                 "described in Adcroft et al., O. Mod. (2008).", default=.true.)

  if (CS%Analytic_FV_PGF) then
    call PressureForce_AFV_init(Time, G, GV, param_file, diag, &
             CS%PressureForce_AFV_CSp, tides_CSp)
  else
    call PressureForce_Mont_init(Time, G, GV, param_file, diag, &
             CS%PressureForce_Mont_CSp, tides_CSp)
  endif

end subroutine PressureForce_init

!> Deallocate the pressure force control structure
subroutine PressureForce_end(CS)
  type(PressureForce_CS), pointer :: CS !< Pressure force control structure
  if (associated(CS)) deallocate(CS)
end subroutine PressureForce_end

!> \namespace mom_pressureforce
!!
!! This thin module provides a branch to two forms of the horizontal accelerations
!! due to pressure gradients. The two options currently available are a
!! Montgomery potential form (used in traditional isopycnal layer models), and the
!! analytic finite volume form.

end module MOM_PressureForce
