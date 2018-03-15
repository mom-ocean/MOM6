!> Interface to CVMix tidal mixing scheme.
module MOM_cvmix_tidal

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator,  only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,  only : post_data
use MOM_EOS,            only : calculate_density
use MOM_variables,      only : thermo_var_ptrs
use MOM_error_handler,  only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,    only : openParameterBlock, closeParameterBlock
use MOM_debugging,      only : hchksum
use MOM_grid,           only : ocean_grid_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_file_parser,    only : get_param, log_version, param_file_type
use cvmix_kpp,          only : CVmix_kpp_compute_kOBL_depth

implicit none ; private

#include <MOM_memory.h>

logical :: debug = is_root_pe .and. .true.

public cvmix_tidal_init
public calculate_cvmix_tidal
public cvmix_tidal_end

!> Control structure including parameters for CVMix tidal mixing.
type, public :: cvmix_tidal_cs

  ! Parameters
  real    :: kd_conv       !< diffusivity constant used in convective regime (m2/s)

end type cvmix_tidal_cs

character(len=40)  :: mdl = "MOM_cvmix_tidal"     !< This module's name.

contains

!> Initialize the cvmix tidal mixing routine.
logical function cvmix_tidal_init(Time, G, GV, param_file, diag, CS)

  type(time_type),          intent(in)    :: Time       !< The current time.
  type(ocean_grid_type),    intent(in)    :: G          !< Grid structure.
  type(verticalGrid_type),  intent(in)    :: GV         !< Vertical grid structure.
  type(param_file_type),    intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target,  intent(inout) :: diag       !< Diagnostics control structure.
  type(cvmix_tidal_cs),     pointer       :: CS         !< This module's control structure.

  ! Local variables

! This include declares and sets the variable "version".
#include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(WARNING, "cvmix_tidal_init called when control structure "// &
                            "is already associated.")
    return
  endif
  allocate(CS)

  ! Read parameters
  call log_version(param_file, mdl, version, &
    "Parameterization of enhanced mixing due to convection via CVMix")
  call get_param(param_file, mdl, "USE_CVMIX_TIDAL", cvmix_tidal_init, &
                 "If true, turns on tidal mixing scheme via CVMix\n", &
                 default=.false.)

  if (.not. cvmix_conv_init) return

  call closeParameterBlock(param_file)

end function cvmix_tidal_init


!> ....
subroutine calculate_cvmix_tidal()
  continue
end subroutine calculate_cvmix_tidal


!> Clear pointers and deallocate memory
subroutine cvmix_tidal_end(CS)
  type(cvmix_tidal_cs), pointer :: CS ! This module's control structure

  !TODO deallocate all the dynamically allocated members here ...
  deallocate(CS)
end subroutine cvmix_tidal_end


end module MOM_cvmix_tidal
