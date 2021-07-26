!> Top-level module for the MOM6 ocean model in coupled mode.
module MOM_stochastics

! This file is part of MOM6. See LICENSE.md for the license.

! This is the top level module for the MOM6 ocean model.  It contains routines
! for initialization, termination and update of ocean model state.  This
! particular version wraps all of the calls for MOM6 in the calls that had
! been used for MOM4.
!
! This code is a stop-gap wrapper of the MOM6 code to enable it to be called
! in the same way as MOM4.

use MOM_diag_mediator,       only : register_diag_field, diag_ctrl, time_type
use MOM_grid,                only : ocean_grid_type
use MOM_verticalGrid,        only : verticalGrid_type
use MOM_error_handler,       only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler,       only : callTree_enter, callTree_leave
use MOM_file_parser,         only : get_param, log_version, close_param_file, param_file_type
use mpp_domains_mod,         only : domain2d, mpp_get_layout, mpp_get_global_domain
use mpp_domains_mod,         only : mpp_define_domains, mpp_get_compute_domain, mpp_get_data_domain
use MOM_domains,             only : root_PE,num_PEs
use MOM_coms,                only : Get_PElist

#include <MOM_memory.h>

implicit none ; private

public stochastics_init, update_stochastics

!> This control structure holds parameters for the MOM_stochastics module
type, public:: stochastic_CS
  logical :: do_sppt         !< If true, stochastically perturb the diabatic
  logical :: pert_epbl       !< If true, then randomly perturb the KE dissipation and genration terms
  !>@{ Diagnostic IDs
  integer :: id_sppt_wts  = -1
  integer :: id_epbl1_wts=-1,id_epbl2_wts=-1
  !>@}
  ! stochastic patterns
  real, allocatable :: sppt_wts(:,:)  !< Random pattern for ocean SPPT
                                     !! tendencies with a number between 0 and 2
  real, allocatable :: epbl1_wts(:,:) !< Random pattern for K.E. generation
  real, allocatable :: epbl2_wts(:,:) !< Random pattern for K.E. dissipation
  type(diag_ctrl), pointer :: diag   !< structure used to regulate timing of diagnostic output
  type(time_type), pointer :: Time !< Pointer to model time (needed for sponges)
end type stochastic_CS

contains

subroutine stochastics_init(dt, grid, GV, CS, param_file, diag, Time)
  real, intent(in)                     :: dt       !< time step [T ~> s]
  type(ocean_grid_type),   intent(in)  :: grid     ! horizontal grid information
  type(verticalGrid_type), intent(in)  :: GV       ! vertical grid structure
  type(stochastic_CS), pointer,     intent(inout):: CS
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag             !< structure to regulate diagnostic output
  type(time_type), target                :: Time             !< model time
  return
end subroutine stochastics_init

subroutine update_stochastics(CS)
  type(stochastic_CS),      intent(inout) :: CS        !< diabatic control structure
  return
end subroutine update_stochastics

end module MOM_stochastics

