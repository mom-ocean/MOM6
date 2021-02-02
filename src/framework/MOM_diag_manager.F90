!> A module for MOM6 diagnostics which interfaces with shared infrastructure
!! (currently FMS) imported from MOM_diag_manager_infra. APIs published
!! in this module are available to the diag_mediators (MOM,MOM_Ice_Shelf,SIS2).
!! send_data_fms_wrapper and register_diag_field_fms_wrapper are the only
!! infrastructure specific modules being published for legacy reasons and to
!! avoid potential namespace conflicts with some compilers.
module MOM_diag_manager

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_manager_infra, only : MOM_diag_axis_init, get_MOM_diag_axis_name
use MOM_diag_manager_infra, only : EAST, NORTH
use MOM_diag_manager_infra, only : null_axis_id
use MOM_diag_manager_infra, only : MOM_diag_manager_init, MOM_diag_manager_end
use MOM_diag_manager_infra, only : MOM_diag_field_add_attribute
use MOM_diag_manager_infra, only : DIAG_FIELD_NOT_FOUND
use MOM_diag_manager_infra, only : send_data_fms_wrapper
use MOM_diag_manager_infra, only : register_diag_field_fms_wrapper
use MOM_diag_manager_infra, only : register_static_field_fms_wrapper
use MOM_diag_manager_infra, only : get_MOM_diag_field_id
use MOM_domain_infra, only : MOM_domain_type, domain2d
use MOM_error_infra,  only : MOM_error=>MOM_err, FATAL
use MOM_time_manager, only : time_type

implicit none ; private

public :: MOM_diag_manager_init, MOM_diag_manager_end
public :: MOM_diag_axis_init, get_MOM_diag_axis_name, EAST, NORTH
public :: send_data_fms_wrapper, MOM_diag_field_add_attribute, DIAG_FIELD_NOT_FOUND
public :: register_diag_field_fms_wrapper, register_static_field_fms_wrapper
public :: get_MOM_diag_field_id

!> \namespace mom_diag_manager
!!
!! This module simply wraps register_diag_field() from FMS's diag_manager_mod.
!! We used to be able to import register_diag_field and rename it to register_diag_field_fms
!! with a simple "use, only : register_diag_field_fms => register_diag_field" but PGI 16.5
!! has a bug that refuses to compile this - earlier versions did work.

end module MOM_diag_manager
