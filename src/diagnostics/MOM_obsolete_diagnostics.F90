!> Provides a mechanism for recording diagnostic variables that are no longer
!! valid, along with their replacement name if appropriate.
module MOM_obsolete_diagnostics

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : diag_ctrl, found_in_diagtable
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,   only : param_file_type, log_version, get_param

implicit none ; private

#include <MOM_memory.h>

public register_obsolete_diagnostics

contains

!> Scan through the diag_table searching for obsolete parameters and issue informational
!! messages and optionallly a FATAL error.
subroutine register_obsolete_diagnostics(param_file, diag)
  type(param_file_type),       intent(in)    :: param_file !< The parameter file handle.
  type(diag_ctrl),             intent(in)    :: diag       !< A structure used to control diagnostics.
! This include declares and sets the variable "version".
#include "version_variable.h"
  ! Local variables
  character(len=40)  :: mdl = "MOM_obsolete_diagnostics"  !< This module's name.
  logical :: foundEntry, causeFatal
  integer :: errType

  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "OBSOLETE_DIAGNOSTIC_IS_FATAL", causeFatal,              &
                 "If an obsolete diagnostic variable appears in the diag_table, "//        &
                 "cause a FATAL error rather than issue a WARNING.", default=.true.)

  foundEntry = .false.
  ! Each obsolete entry, with replacement name is available.
  if (diag_found(diag, 'Net_Heat', 'net_heat_surface or net_heat_coupler')) foundEntry = .true.
  if (diag_found(diag, 'PmE', 'PRCmE'))                                     foundEntry = .true.
  if (diag_found(diag, 'froz_precip', 'fprec'))                             foundEntry = .true.
  if (diag_found(diag, 'liq_precip', 'lprec'))                              foundEntry = .true.
  if (diag_found(diag, 'virt_precip', 'vprec'))                             foundEntry = .true.
  if (diag_found(diag, 'froz_runoff', 'frunoff'))                           foundEntry = .true.
  if (diag_found(diag, 'liq_runoff', 'lrunoff'))                            foundEntry = .true.
  if (diag_found(diag, 'calving_heat_content', 'heat_content_frunoff'))     foundEntry = .true.
  if (diag_found(diag, 'precip_heat_content', 'heat_content_lprec'))        foundEntry = .true.
  if (diag_found(diag, 'evap_heat_content', 'heat_content_massout'))        foundEntry = .true.
  if (diag_found(diag, 'runoff_heat_content', 'heat_content_lrunoff'))      foundEntry = .true.
  if (diag_found(diag, 'latent_fprec'))                                     foundEntry = .true.
  if (diag_found(diag, 'latent_calve'))                                     foundEntry = .true.
  if (diag_found(diag, 'heat_rest', 'heat_restore'))                        foundEntry = .true.
  if (diag_found(diag, 'KPP_dTdt', 'KPP_NLT_dTdt'))                         foundEntry = .true.
  if (diag_found(diag, 'KPP_dSdt', 'KPP_NLT_dSdt'))                         foundEntry = .true.

  if (causeFatal) then; errType = FATAL
  else ; errType = WARNING ; endif
  if (foundEntry .and. is_root_pe()) &
    call MOM_error(errType, 'MOM_obsolete_diagnostics: Obsolete diagnostics found in diag_table.')

end subroutine register_obsolete_diagnostics

!> Determines whether an obsolete parameter appears in the diag_table.
logical function diag_found(diag, varName, newVarName)
  type(diag_ctrl),            intent(in) :: diag       !< A structure used to control diagnostics.
  character(len=*),           intent(in) :: varName    !< The obsolete diagnostic name
  character(len=*), optional, intent(in) :: newVarName !< The valid name of this diagnostic
  ! Local
  integer :: handle ! Integer handle returned from diag_manager

  diag_found = found_in_diagtable(diag, varName)

  if (diag_found .and. is_root_pe()) then
    if (present(newVarName)) then
      call MOM_error(WARNING, 'MOM_obsolete_params: '//'diag_table entry "'// &
          trim(varName)//'" found. Use ''"'//trim(newVarName)//'" instead.' )
    else
      call MOM_error(WARNING, 'MOM_obsolete_params: '//'diag_table entry "'// &
          trim(varName)//'" is obsolete.' )
    endif
  endif

end function diag_found

end module MOM_obsolete_diagnostics
