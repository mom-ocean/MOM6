!> Provides a mechanism for recording diagnostic variables that are no longer
!! valid, along with their replacement name if appropriate.
module MOM_obsolete_diagnostics

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,   only : param_file_type, log_version, get_param
use MOM_diag_mediator, only : diag_ctrl
use diag_manager_mod, only  : register_static_field_fms=>register_static_field

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
  character(len=40)  :: mod = "MOM_obsolete_diagnostics"  !< This module's name.
  logical :: foundEntry, causeFatal
  integer :: errType

  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "OBSOLETE_DIAGNOSTIC_IS_FATAL", causeFatal,              &
                 "If an obsolete diagnostic variable appears in the diag_table\n"//        &
                 "then cause a FATAL error rather than issue a WARNING.", default=.true.)

  foundEntry = .false.
  ! Each obsolete entry, with replacement name is available.
  if (found_in_diagtable(diag, 'Net_Heat', 'net_heat_surface or net_heat_coupler')) foundEntry = .true.
  if (found_in_diagtable(diag, 'PmE', 'PRCmE'))                                     foundEntry = .true.
  if (found_in_diagtable(diag, 'froz_precip', 'fprec'))                             foundEntry = .true.
  if (found_in_diagtable(diag, 'liq_precip', 'lprec'))                              foundEntry = .true.
  if (found_in_diagtable(diag, 'virt_precip', 'vprec'))                             foundEntry = .true.
  if (found_in_diagtable(diag, 'froz_runoff', 'frunoff'))                           foundEntry = .true.
  if (found_in_diagtable(diag, 'liq_runoff', 'lrunoff'))                            foundEntry = .true.
  if (found_in_diagtable(diag, 'calving_heat_content', 'heat_content_frunoff'))     foundEntry = .true.
  if (found_in_diagtable(diag, 'precip_heat_content', 'heat_content_lprec'))        foundEntry = .true.
  if (found_in_diagtable(diag, 'evap_heat_content', 'heat_content_massout'))        foundEntry = .true.
  if (found_in_diagtable(diag, 'runoff_heat_content', 'heat_content_lrunoff'))      foundEntry = .true.
  if (found_in_diagtable(diag, 'latent_fprec'))                                     foundEntry = .true.
  if (found_in_diagtable(diag, 'latent_calve'))                                     foundEntry = .true.
  if (found_in_diagtable(diag, 'heat_rest', 'heat_restore'))                        foundEntry = .true.
  if (found_in_diagtable(diag, 'KPP_dTdt', 'KPP_NLT_dTdt'))                         foundEntry = .true.
  if (found_in_diagtable(diag, 'KPP_dSdt', 'KPP_NLT_dSdt'))                         foundEntry = .true.

  if (causeFatal) then; errType = FATAL
  else ; errType = WARNING ; endif
  if (foundEntry .and. is_root_pe()) &
    call MOM_error(errType, 'MOM_obsolete_diagnostics: '//&
                            'Obsolete diagnostics found in diag_table')

end subroutine register_obsolete_diagnostics

!> Fakes a register of a diagnostic to find out if an obsolete
!! parameter appears in the diag_table.
logical function found_in_diagtable(diag, varName, newVarName)
  type(diag_ctrl),            intent(in) :: diag
  character(len=*),            intent(in) :: varName
  character(len=*), optional, intent(in) :: newVarName
  ! Local
  integer :: handle ! Integer handle returned from diag_manager

  ! We use register_static_field_fms() instead of register_static_field() so
  ! that the diagnostic does not appear in the available diagnostics list.
  handle = register_static_field_fms('ocean_model', varName, &
            diag%axesT1%handles, 'Obsolete parameter', 'N/A')

  found_in_diagtable = (handle>0)

  if (handle>0 .and. is_root_pe()) then
    if (present(newVarName)) then
      call MOM_error(WARNING, 'MOM_obsolete_params: '//                        &
          'diag_table entry "'//trim(varName)//'" found. Use '// &
          '"'//trim(newVarName)//'" instead.' )
    else
      call MOM_error(WARNING, 'MOM_obsolete_params: '//                        &
          'diag_table entry "'//trim(varName)//'" is obsolete.' )
    endif
  endif

end function found_in_diagtable

end module MOM_obsolete_diagnostics
