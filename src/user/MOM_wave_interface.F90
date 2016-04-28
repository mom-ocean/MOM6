!> Initial conditions and forcing for the single column model (SCM) CVmix
!! test set.
module MOM_wave_interface

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing, allocate_forcing_type
use MOM_grid, only : ocean_grid_type
use MOM_verticalgrid, only: verticalGrid_type
use MOM_safe_alloc, only : safe_alloc_ptr
use MOM_time_manager, only : time_type, operator(+), operator(/), get_time,&
                             time_type_to_real
use MOM_variables, only : thermo_var_ptrs, surface
implicit none ; private

#include <MOM_memory.h>

public MOM_wave_interface_init


!> Container for wave related parameters
type, public:: wave_parameters_CS ;
private
  logical :: UseWaves  !< True to Compute Wave parameters
  integer :: WaveMode  !< Options for various wave methods
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
       Us_x ! Stokes drift (zonal) 
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
       Us_y ! Stokes drift (meridional) 
end type

! This include declares and sets the variable "version".
#include "version_variable.h"

character(len=40)  :: mod = "MOM_wave_interface" ! This module's name.

contains

!> Initializes parameters related to MOM_wave_interface
subroutine MOM_wave_interface_init(G,GV,param_file, CS)
  type(ocean_grid_type),                  intent(in)  :: G !< Grid structure
  type(verticalGrid_type),                intent(in)  :: GV!< Vertical grid structure
  type(param_file_type),                  intent(in)  :: param_file !< Input parameter structure
  type(wave_parameters_CS),              pointer     :: CS
  ! Local variables
  integer :: is, ie, js, je, isd, ied, jsd, jed, isdB, iedB, jsdB, jedB, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
     call MOM_error(WARNING, "wave_interface_init called with an associated"//&
                             "control structure.")
     return
  endif
  
  allocate(CS)

  call log_version(param_file, mod, version)
  call get_param(param_file,mod,"Use_Waves",CS%UseWaves, &
                 'Main switch to use wave input', units='',default=.false.)
  if (CS%UseWaves) then 
     !allocate and initialize Stokes drift
     ALLOC_ (CS%Us_x(isdB:IedB,jsd:jed,nz)) ; CS%Us_x(:,:,:) = 0.0
     ALLOC_ (CS%Us_y(isd:Ied,jsdB:jedB,nz)) ; CS%Us_y(:,:,:) = 0.0 
  endif

  print*,'Brandon you are here'
  stop

end subroutine MOM_wave_interface_init

end module MOM_wave_interface
