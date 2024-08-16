!> Top-level module for the MOM6 ocean model in coupled mode.
module MOM_stochastics

! This file is part of MOM6. See LICENSE.md for the license.

! This is the top level module for the MOM6 ocean model.  It contains routines
! for initialization, update, and writing restart of stochastic physics. This
! particular version wraps all of the calls for MOM6 in the calls that had
! been used for MOM4.
!
use MOM_diag_mediator,       only : register_diag_field, diag_ctrl, time_type
use MOM_grid,                only : ocean_grid_type
use MOM_verticalGrid,        only : verticalGrid_type
use MOM_error_handler,       only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler,       only : callTree_enter, callTree_leave
use MOM_file_parser,         only : get_param, log_version, close_param_file, param_file_type
use mpp_domains_mod,         only : domain2d, mpp_get_layout, mpp_get_global_domain
use mpp_domains_mod,         only : mpp_define_domains, mpp_get_compute_domain, mpp_get_data_domain
use MOM_domains,             only : root_PE, num_PEs
use MOM_coms,                only : Get_PElist
use stochastic_physics,      only : init_stochastic_physics_ocn, run_stochastic_physics_ocn

#include <MOM_memory.h>

implicit none ; private

public stochastics_init, update_stochastics

!> This control structure holds parameters for the MOM_stochastics module
type, public:: stochastic_CS
  logical :: do_sppt         !< If true, stochastically perturb the diabatic
  logical :: pert_epbl       !< If true, then randomly perturb the KE dissipation and genration terms
  integer :: id_sppt_wts  = -1 !< Diagnostic id for SPPT
  integer :: id_epbl1_wts = -1 !< Diagnostic id for epbl generation perturbation
  integer :: id_epbl2_wts = -1 !< Diagnostic id for epbl dissipation perturbation
  ! stochastic patterns
  real, allocatable :: sppt_wts(:,:)  !< Random pattern for ocean SPPT
                                      !! tendencies with a number between 0 and 2 [nondim]
  real, allocatable :: epbl1_wts(:,:) !< Random pattern for K.E. generation [nondim]
  real, allocatable :: epbl2_wts(:,:) !< Random pattern for K.E. dissipation [nondim]
  type(diag_ctrl), pointer :: diag   !< structure used to regulate timing of diagnostic output
  type(time_type), pointer :: Time !< Pointer to model time (needed for sponges)
end type stochastic_CS

contains

!!   This subroutine initializes the stochastics physics control structure.
subroutine stochastics_init(dt, grid, GV, CS, param_file, diag, Time)
  real, intent(in)                       :: dt      !< time step [T ~> s]
  type(ocean_grid_type),   intent(in)    :: grid    !< horizontal grid information
  type(verticalGrid_type), intent(in)    :: GV      !< vertical grid structure
  type(stochastic_CS), pointer, intent(inout) :: CS !< stochastic control structure
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag    !< structure to regulate diagnostic output
  type(time_type), target                :: Time    !< model time

  ! Local variables
  integer, allocatable :: pelist(:) ! list of pes for this instance of the ocean
  integer :: mom_comm          ! list of pes for this instance of the ocean
  integer :: num_procs         ! number of processors to pass to stochastic physics
  integer :: iret              ! return code from stochastic physics
  integer :: pe_zero           !  root pe
  integer :: nx                ! number of x-points including halo
  integer :: ny                ! number of x-points including halo

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "ocean_stochastics_init"  ! This module's name.

  call callTree_enter("ocean_model_stochastic_init(), MOM_stochastics.F90")
  if (associated(CS)) then
    call MOM_error(WARNING, "MOM_stochastics_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag
  CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  ! get number of processors and PE list for stochastic physics initialization
  call get_param(param_file, mdl, "DO_SPPT", CS%do_sppt, &
                 "If true, then stochastically perturb the thermodynamic "//&
                 "tendemcies of T,S, amd h.  Amplitude and correlations are "//&
                 "controlled by the nam_stoch namelist in the UFS model only.", &
                 default=.false.)
  call get_param(param_file, mdl, "PERT_EPBL", CS%pert_epbl, &
                 "If true, then stochastically perturb the kinetic energy "//&
                 "production and dissipation terms.  Amplitude and correlations are "//&
                 "controlled by the nam_stoch namelist in the UFS model only.", &
                 default=.false.)
  if (CS%do_sppt .OR. CS%pert_epbl) then
    num_procs = num_PEs()
    allocate(pelist(num_procs))
    call Get_PElist(pelist,commID = mom_comm)
    pe_zero = root_PE()
    nx = grid%ied - grid%isd + 1
    ny = grid%jed - grid%jsd + 1
    call init_stochastic_physics_ocn(dt,grid%geoLonT,grid%geoLatT,nx,ny,GV%ke, &
                                     CS%pert_epbl,CS%do_sppt,pe_zero,mom_comm,iret)
    if (iret/=0)  then
      call MOM_error(FATAL, "call to init_stochastic_physics_ocn failed")
    endif

    if (CS%do_sppt) allocate(CS%sppt_wts(grid%isd:grid%ied,grid%jsd:grid%jed), source=0.0)
    if (CS%pert_epbl) then
      allocate(CS%epbl1_wts(grid%isd:grid%ied,grid%jsd:grid%jed), source=0.0)
      allocate(CS%epbl2_wts(grid%isd:grid%ied,grid%jsd:grid%jed), source=0.0)
    endif
  endif
  if (CS%do_sppt) then
    CS%id_sppt_wts = register_diag_field('ocean_model', 'sppt_pattern', CS%diag%axesT1, Time, &
         'random pattern for sppt', 'None')
  endif
  if (CS%pert_epbl) then
    CS%id_epbl1_wts = register_diag_field('ocean_model', 'epbl1_wts', CS%diag%axesT1, Time, &
        'random pattern for KE generation', 'None')
    CS%id_epbl2_wts = register_diag_field('ocean_model', 'epbl2_wts', CS%diag%axesT1, Time, &
        'random pattern for KE dissipation', 'None')
  endif

  if (CS%do_sppt .OR. CS%pert_epbl) &
    call MOM_mesg('            === COMPLETED MOM STOCHASTIC INITIALIZATION =====')

  call callTree_leave("ocean_model_init(")

end subroutine stochastics_init

!> update_ocean_model uses the forcing in Ice_ocean_boundary to advance the
!! ocean model's state from the input value of Ocean_state (which must be for
!! time time_start_update) for a time interval of Ocean_coupling_time_step,
!! returning the publicly visible ocean surface properties in Ocean_sfc and
!! storing the new ocean properties in Ocean_state.
subroutine update_stochastics(CS)
  type(stochastic_CS),      intent(inout) :: CS        !< diabatic control structure
  call callTree_enter("update_stochastics(), MOM_stochastics.F90")

! update stochastic physics patterns before running next time-step
  call run_stochastic_physics_ocn(CS%sppt_wts,CS%epbl1_wts,CS%epbl2_wts)

  return
end subroutine update_stochastics

end module MOM_stochastics

