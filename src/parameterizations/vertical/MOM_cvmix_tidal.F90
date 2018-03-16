!> Interface to CVMix tidal mixing scheme.
module MOM_cvmix_tidal

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator,  only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,  only : post_data
use MOM_EOS,            only : calculate_density
use MOM_variables,      only : thermo_var_ptrs
use MOM_error_handler,  only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_debugging,      only : hchksum
use MOM_grid,           only : ocean_grid_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_file_parser,    only : openParameterBlock, closeParameterBlock
use MOM_file_parser,    only : get_param, log_version, param_file_type
use cvmix_tidal,        only : cvmix_init_tidal

implicit none ; private

#include <MOM_memory.h>

public cvmix_tidal_init
public calculate_cvmix_tidal
public cvmix_tidal_end

!> Control structure including parameters for CVMix tidal mixing.
type, public :: cvmix_tidal_cs
  logical :: debug = .true.

  ! Parameters
  real :: local_mixing_frac !< fraction of wave energy dissipated locally.
  real :: mixing_efficiency !< The efficiency that mechanical energy dissipation translates into mixing
                            !! that can be parameterized by a diffusivity acting on vertical stratification.
  real :: vert_decay_scale  !< zeta in the Simmons paper (to compute the vertical deposition function). [m]
  real :: tidal_max_coef    !< maximum allowable tidel diffusivity. [m^2/s]

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

  CS%debug = CS%debug.and.is_root_pe()

  ! Read parameters
  call log_version(param_file, mdl, version, &
    "Parameterization of tidal mixing via CVMix")
  call get_param(param_file, mdl, "USE_CVMIX_TIDAL", cvmix_tidal_init, &
                 "If true, turns on tidal mixing scheme via CVMix", &
                 default=.false.)
  call openParameterBlock(param_file,'CVMIX_TIDAL')
  call get_param(param_file, mdl, "LOCAL_MIXING_FRAC", CS%local_mixing_frac, &
                 "Fraction of wave energy dissipated locally.", &
                 units="nondim", default=0.33)
  call get_param(param_file, mdl, "MIXING_EFFICIENCY", CS%mixing_efficiency, &
                 "Gamma in Simmons, 2004", &
                 units="nondim", default=0.20)
  !TODO: make sure GAMMA_ITIDES (same as LOCAL_MIXING_FRAC
  call get_param(param_file, mdl, "VERTICAL_DECAY_SCALE", CS%vert_decay_scale, &
                 "zeta in Simmons, 2004. Used to compute the vertical deposition function", &
                 units="m", default=500.0)
  !TODO: make sure int_tide_decay scale (same as VERTICAL_DECAY_SCALE is removed from code).
  call get_param(param_file, mdl, "TIDAL_MAX_COEF", CS%tidal_max_coef, &
                 "largest acceptable value for tidal diffusivity", &
                 units="m^2/s", default=100e-4) ! the default is 50e-4 in CVMIX, 100e-4 in POP.
  call closeParameterBlock(param_file)

  if (.not. cvmix_tidal_init) return

  if (CS%debug) print *, __FILE__, __LINE__, cvmix_tidal_init

  ! Set up CVMix
  call cvmix_init_tidal(mix_scheme            = 'Simmons',            &
                        efficiency            = cs%mixing_efficiency, &
                        vertical_decay_scale  = cs%vert_decay_scale,  &
                        max_coefficient       = cs%tidal_max_coef,    &
                        local_mixing_frac     = cs%local_mixing_frac, &
                        depth_cutoff          = 0.0)
                        



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
