!> A simple (very thin) wrapper for the FMS diag_manager routines, with some name changes
module MOM_diag_manager

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_manager_infra,    only : diag_axis_init, get_diag_axis_name, EAST, NORTH
use MOM_diag_manager_infra,    only : null_axis_id
use MOM_diag_manager_infra, only : diag_manager_init, diag_manager_end
use MOM_diag_manager_infra, only : send_data, diag_field_add_attribute, DIAG_FIELD_NOT_FOUND
use MOM_diag_manager_infra, only : register_diag_field_fms
use MOM_diag_manager_infra, only : register_static_field_fms
use MOM_diag_manager_infra, only : get_diag_field_id_fms
use MOM_domain_infra, only : MOM_domain_type, domain2d
use MOM_error_infra,  only : MOM_error=>MOM_err, FATAL
use MOM_time_manager, only : time_type

implicit none ; private

public :: diag_manager_init, diag_manager_end
public :: diag_axis_init, get_diag_axis_name, EAST, NORTH
public :: send_data, diag_field_add_attribute, DIAG_FIELD_NOT_FOUND
public :: register_diag_field_fms, register_static_field_fms,  get_diag_field_id_fms

!> \namespace mom_diag_manager
!!
!! This module simply wraps register_diag_field() from FMS's diag_manager_mod.
!! We used to be able to import register_diag_field and rename it to register_diag_field_fms
!! with a simple "use, only : register_diag_field_fms => register_diag_field" but PGI 16.5
!! has a bug that refuses to compile this - earlier versions did work.

end module MOM_diag_manager
