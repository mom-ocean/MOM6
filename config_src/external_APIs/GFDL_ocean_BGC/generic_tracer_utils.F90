module g_tracer_utils
#include <fms_platform.h>

  use coupler_types_mod, only: coupler_2d_bc_type
  use time_manager_mod, only : time_type

  use field_manager_mod, only: fm_string_len

  use MOM_diag_mediator, only : register_diag_field_MOM=>register_diag_field
  use MOM_diag_mediator, only : post_data_MOM=>post_data, post_data_1d_k
  use MOM_diag_mediator, only : g_diag_ctrl=>diag_ctrl


implicit none ; private

  type g_tracer_type
!   ! Tracer concentration field in space (and time)
!   ! MOM keeps the prognostic tracer fields at 3 time levels, hence 4D.
    real, pointer, dimension(:,:,:,:) :: field  => NULL()
!   ! Tracer concentration in river runoff
    real, _ALLOCATABLE, dimension(:,:)    :: trunoff _NULL
    logical :: requires_restart = .true.
    character(len=fm_string_len) :: src_file, src_var_name, src_var_unit, src_var_gridspec
    integer :: src_var_record
    logical :: requires_src_info = .false.
    real    :: src_var_unit_conversion = 1.0 !This factor depends on the tracer. Ask  Jasmin
    real    :: src_var_valid_min = 0.0
  end type g_tracer_type


  type g_diag_type
    integer :: dummy
  end type g_diag_type

  type g_tracer_common
    type(g_diag_ctrl) :: diag_CS
    !Domain extents
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk
  end type g_tracer_common

  type(g_tracer_common), target, save :: g_tracer_com

  public :: g_tracer_type
  public :: g_tracer_find
  public :: g_tracer_add
  public :: g_tracer_init
  public :: g_tracer_flux_init
  public :: g_tracer_column_int
  public :: g_tracer_flux_at_depth
  public :: g_tracer_add_param
  public :: g_tracer_set_values
  public :: g_tracer_get_values
  public :: g_tracer_get_pointer
  public :: g_tracer_get_common
  public :: g_tracer_set_common
  public :: g_tracer_set_csdiag
  public :: g_tracer_set_files
  public :: g_tracer_coupler_set
  public :: g_tracer_coupler_get
  public :: g_tracer_send_diag
  public :: g_tracer_diag
  public :: g_tracer_get_name
  public :: g_tracer_get_alias
  public :: g_tracer_get_next
  public :: g_tracer_register_diag
  public :: g_tracer_is_prog
  public :: g_tracer_vertdiff_G
  public :: g_tracer_vertdiff_M
  public :: g_tracer_start_param_list
  public :: g_tracer_end_param_list
  public :: g_diag_type
  public :: g_diag_field_add
  public :: g_tracer_set_pointer
  public :: g_tracer_print_info
  public :: g_tracer_coupler_accumulate
  public :: g_tracer_get_src_info
  public :: g_register_diag_field
  public :: g_send_data

  interface g_tracer_add_param
     module procedure g_tracer_add_param_real
     module procedure g_tracer_add_param_logical
     module procedure g_tracer_add_param_integer
     module procedure g_tracer_add_param_string
  end interface

  interface g_tracer_set_pointer
    module procedure g_tracer_set_pointer_3d
    module procedure g_tracer_set_pointer_4d
  end interface g_tracer_set_pointer

  interface g_send_data
     module procedure g_send_data_0d
     module procedure g_send_data_1d
     module procedure g_send_data_2d
     module procedure g_send_data_3d
  end interface

  interface g_tracer_set_values
     module procedure g_tracer_set_real
     module procedure g_tracer_set_2D
     module procedure g_tracer_set_3D
     module procedure g_tracer_set_4D
  end interface

  interface g_tracer_get_values
     module procedure g_tracer_get_4D_val
     module procedure g_tracer_get_3D_val
     module procedure g_tracer_get_2D_val
     module procedure g_tracer_get_real
     module procedure g_tracer_get_string
  end interface

  interface g_tracer_get_pointer
     module procedure g_tracer_get_4D
     module procedure g_tracer_get_3D
     module procedure g_tracer_get_2D
  end interface

contains

  subroutine g_tracer_start_param_list(package_name)
    character(len=fm_string_len), intent(in) :: package_name
    character(len=fm_string_len), parameter  :: sub_name = 'g_tracer_start_param_list'
    character(len=fm_string_len) :: list_path
    integer                      :: list_index
  end subroutine g_tracer_start_param_list

  subroutine g_tracer_end_param_list(package_name)
    character(len=fm_string_len) :: package_name
  end subroutine g_tracer_end_param_list

  subroutine g_tracer_add_param_real(name, var,  value)
    character(len=*), intent(in)  :: name
    real,             intent(in)  :: value
    real,             intent(out) :: var
  end subroutine g_tracer_add_param_real

  subroutine g_tracer_add_param_logical(name, var,  value)
    character(len=*), intent(in)  :: name
    logical,          intent(in)  :: value
    logical,          intent(out) :: var
  end subroutine g_tracer_add_param_logical

  subroutine g_tracer_add_param_integer(name, var,  value)
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: value
    integer,          intent(out) :: var
  end subroutine g_tracer_add_param_integer

  subroutine g_tracer_add_param_string(name, var,  value)
    character(len=*), intent(in)  :: name
    character(len=*), intent(in)  :: value
    character(len=*), intent(out) :: var
  end subroutine g_tracer_add_param_string

  subroutine g_tracer_add(node_ptr, package, name, longname, units,  prog, const_init_value,init_value,&
       flux_gas, flux_gas_name, flux_runoff, flux_wetdep, flux_drydep, flux_gas_molwt, flux_gas_param, &
       flux_param, flux_bottom, btm_reservoir, move_vertical, diff_vertical, sink_rate, flux_gas_restart_file, &
       flux_gas_type, requires_src_info, standard_name,diag_name, diag_field_units,diag_field_scaling_factor, &
       implementation)
    type(g_tracer_type), pointer :: node_ptr
    character(len=*),   intent(in) :: package,name,longname,units
    logical,            intent(in) :: prog
    real,               intent(in), optional :: const_init_value
    real,               intent(in), optional :: init_value
    real,               intent(in), optional :: sink_rate
    logical,            intent(in), optional :: flux_gas
    logical,            intent(in), optional :: flux_runoff
    logical,            intent(in), optional :: flux_wetdep
    logical,            intent(in), optional :: flux_drydep
    logical,            intent(in), optional :: flux_bottom
    logical,            intent(in), optional :: btm_reservoir
    logical,            intent(in), optional :: move_vertical
    logical,            intent(in), optional :: diff_vertical
    real,               intent(in), optional :: flux_gas_molwt
    real, dimension(:), intent(in), optional :: flux_gas_param
    real, dimension(:), intent(in), optional :: flux_param
    character(len=*),   intent(in), optional :: flux_gas_name
    character(len=*),   intent(in), optional :: implementation
    character(len=*),   intent(in), optional :: flux_gas_type
    character(len=*),   intent(in), optional :: flux_gas_restart_file
    logical,            intent(in), optional :: requires_src_info
    character(len=*),   intent(in), optional :: standard_name
    character(len=*),   intent(in), optional :: diag_name
    character(len=*),   intent(in), optional :: diag_field_units
    real,               intent(in), optional :: diag_field_scaling_factor
  end subroutine g_tracer_add

  function remap_bounds(ilb, jlb, klb, array) result(ptr)
  real, dimension(:,:,:),          pointer              :: ptr
  integer,                                 intent(in)   :: ilb
  integer,                                 intent(in)   :: jlb
  integer,                                 intent(in)   :: klb
  real, dimension(ilb:,jlb:,klb:), target, intent(in)   :: array
  end function remap_bounds

  subroutine g_tracer_init(g_tracer)
    type(g_tracer_type), pointer :: g_tracer
    integer :: isc,iec,jsc,jec,isd,ied,jsd,jed, nk,ntau,axes(3)
  end subroutine g_tracer_init

  subroutine g_tracer_flux_init(g_tracer)
    type(g_tracer_type), pointer :: g_tracer
  end subroutine g_tracer_flux_init

  subroutine g_tracer_register_diag(g_tracer)
    type(g_tracer_type), pointer :: g_tracer
  end subroutine g_tracer_register_diag

  subroutine g_tracer_coupler_set(g_tracer_list,IOB_struc,value)
    type(g_tracer_type), pointer :: g_tracer_list,g_tracer
    type(coupler_2d_bc_type), intent(inout) :: IOB_struc
    real, optional :: value
  end subroutine g_tracer_coupler_set

  subroutine g_tracer_coupler_get(g_tracer_list,IOB_struc, weight, model_time)
    type(g_tracer_type),          pointer :: g_tracer_list, g_tracer
    type(coupler_2d_bc_type),    intent(in) :: IOB_struc
    type(time_type),    optional,intent(in) :: model_time
    real,               optional,intent(in) :: weight
  end subroutine g_tracer_coupler_get

  subroutine g_tracer_coupler_accumulate(g_tracer_list,IOB_struc, weight, model_time)
    type(g_tracer_type),          pointer    :: g_tracer_list, g_tracer
    type(coupler_2d_bc_type),    intent(in)  :: IOB_struc
    real,                        intent(in)  :: weight
    type(time_type), optional,   intent(in)  :: model_time
  end subroutine g_tracer_coupler_accumulate

  subroutine g_tracer_set_csdiag(diag_CS)
    type(g_diag_ctrl),  target,intent(in) :: diag_CS
    g_tracer_com%diag_CS = diag_CS
  end subroutine g_tracer_set_csdiag

  subroutine g_tracer_set_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,grid_kmt,init_time)
    integer,                     intent(in) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes(3)
    real, dimension(isd:,jsd:,:),intent(in) :: grid_tmask
    integer,dimension(isd:,jsd:),intent(in) :: grid_kmt
    type(time_type),             intent(in) :: init_time
  end subroutine g_tracer_set_common

  subroutine g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
       axes,grid_tmask,grid_mask_coast,grid_kmt,init_time,diag_CS)
    integer,               intent(out) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau
    integer,optional,      intent(out) :: axes(3)
    type(time_type), optional,      intent(out) :: init_time
    real, optional, dimension(:,:,:),pointer    :: grid_tmask
    integer, optional, dimension(:,:),  pointer :: grid_mask_coast
    integer, optional, dimension(:,:),  pointer :: grid_kmt
    type(g_diag_ctrl), optional,        pointer :: diag_CS
  end subroutine g_tracer_get_common

  subroutine g_tracer_get_diagCS(diag_CS)
    type(g_diag_ctrl),        pointer :: diag_CS
  end subroutine g_tracer_get_diagCS

  subroutine g_tracer_set_files(ice_restart_file,ocean_restart_file)
    character(len=*),   intent(in) :: ice_restart_file
    character(len=*),   intent(in) :: ocean_restart_file
  end subroutine g_tracer_set_files

  subroutine g_tracer_get_4D(g_tracer_list,name,member,array_ptr)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    real, dimension(:,:,:,:), pointer    :: array_ptr
  end subroutine g_tracer_get_4D

  subroutine g_tracer_get_3D(g_tracer_list,name,member,array_ptr)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    real, dimension(:,:,:), pointer    :: array_ptr
  end subroutine g_tracer_get_3D

  subroutine g_tracer_get_2D(g_tracer_list,name,member,array_ptr)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    real, dimension(:,:), pointer    :: array_ptr
  end subroutine g_tracer_get_2D

  subroutine g_tracer_get_4D_val(g_tracer_list,name,member,array,isd,jsd)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    integer,                  intent(in) :: isd,jsd
    real, dimension(isd:,jsd:,:,:), intent(out):: array
  end subroutine g_tracer_get_4D_val

  subroutine g_tracer_get_3D_val(g_tracer_list,name,member,array,isd,jsd,ntau,positive)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    integer,                  intent(in) :: isd,jsd
    integer, optional,        intent(in) :: ntau
    logical, optional,        intent(in) :: positive
    real, dimension(isd:,jsd:,:), intent(out):: array
    integer :: tau
    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_3D_val'
  end subroutine g_tracer_get_3D_val

  subroutine g_tracer_get_2D_val(g_tracer_list,name,member,array,isd,jsd)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    integer,                  intent(in) :: isd,jsd
    real, dimension(isd:,jsd:), intent(out):: array
  end subroutine g_tracer_get_2D_val

  subroutine g_tracer_get_real(g_tracer_list,name,member,value)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    real,                     intent(out):: value
  end subroutine g_tracer_get_real

  subroutine g_tracer_get_string(g_tracer_list,name,member,string)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    character(len=fm_string_len), intent(out) :: string
    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_string'
  end subroutine g_tracer_get_string

  subroutine g_tracer_set_2D(g_tracer_list,name,member,array,isd,jsd,weight)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),      pointer    :: g_tracer_list, g_tracer
    integer,                   intent(in) :: isd,jsd
    real, dimension(isd:,jsd:),intent(in) :: array
    real, optional            ,intent(in) :: weight
  end subroutine g_tracer_set_2D

  subroutine g_tracer_set_3D(g_tracer_list,name,member,array,isd,jsd,ntau)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    integer,                  intent(in) :: isd,jsd
    integer, optional,        intent(in) :: ntau
    real, dimension(isd:,jsd:,:), intent(in)       :: array
  end subroutine g_tracer_set_3D

  subroutine g_tracer_set_4D(g_tracer_list,name,member,array,isd,jsd)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    integer,                  intent(in) :: isd,jsd
    real, dimension(isd:,jsd:,:,:), intent(in)       :: array
  end subroutine g_tracer_set_4D

  subroutine g_tracer_set_real(g_tracer_list,name,member,value)
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: member
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    real,                     intent(in) :: value
  end subroutine g_tracer_set_real

  subroutine g_tracer_set_pointer_4D(g_tracer_list,name,member,array,ilb,jlb)
    character(len=*),               intent(in) :: name
    character(len=*),               intent(in) :: member
    type(g_tracer_type),            pointer    :: g_tracer_list, g_tracer
    integer,                        intent(in) :: ilb,jlb
    real, dimension(ilb:,jlb:,:,:), target, intent(in) :: array
  end subroutine g_tracer_set_pointer_4D

  subroutine g_tracer_set_pointer_3D(g_tracer_list,name,member,array,ilb,jlb)
    character(len=*),               intent(in) :: name
    character(len=*),               intent(in) :: member
    type(g_tracer_type),            pointer    :: g_tracer_list, g_tracer
    integer,                        intent(in) :: ilb,jlb
    real, dimension(ilb:,jlb:,:), target, intent(in) :: array
  end subroutine g_tracer_set_pointer_3D

  subroutine g_tracer_find(g_tracer,name)
    character(len=*),         intent(in) :: name
    type(g_tracer_type),    pointer    :: g_tracer
  end subroutine g_tracer_find

  subroutine g_tracer_column_int(depth, ilb, jlb, var, dzt, rho_dzt, rd, k_level, integral, caller)
    real,                         intent(in)            :: depth
    integer,                      intent(in)            :: ilb
    integer,                      intent(in)            :: jlb
    real, dimension(ilb:,jlb:,:), intent(in)            :: var
    real, dimension(ilb:,jlb:,:), intent(in)            :: dzt
    real, dimension(ilb:,jlb:,:), intent(in)            :: rho_dzt
    real, dimension(ilb:,jlb:,:), intent(inout)         :: rd
    integer,                      intent(inout)         :: k_level
    real, dimension(ilb:,jlb:),   intent(out)           :: integral
    character(len=*),             intent(in), optional  :: caller
  end subroutine g_tracer_column_int

  subroutine g_tracer_flux_at_depth(depth, ilb, jlb, var, dzt, k_level, frac, initialized, flux, caller)
    real,                            intent(in)                 :: depth
    integer,                         intent(in)                 :: ilb
    integer,                         intent(in)                 :: jlb
    real,    dimension(ilb:,jlb:,:), intent(in)                 :: var
    real,    dimension(ilb:,jlb:,:), intent(in)                 :: dzt
    integer, dimension(ilb:,jlb:),   intent(inout)              :: k_level
    real,    dimension(ilb:,jlb:),   intent(inout)              :: frac
    logical,                         intent(inout)              :: initialized
    real,    dimension(ilb:,jlb:),   intent(out)                :: flux
    character(len=*),                intent(in),    optional    :: caller
  end subroutine g_tracer_flux_at_depth

  subroutine g_tracer_send_diag(g_tracer_list,model_time,tau)
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    type(time_type),        intent(in) :: model_time
    integer,                intent(in) :: tau
  end subroutine g_tracer_send_diag

  subroutine g_tracer_diag(g_tracer_list, ilb, jlb, rho_dzt_tau, rho_dzt_taup1, model_time, tau, taup1, dtts)
    type(g_tracer_type),    pointer    :: g_tracer_list
    integer,                  intent(in) :: ilb
    integer,                  intent(in) :: jlb
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt_tau
    real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt_taup1
    type(time_type),        intent(in) :: model_time
    integer,                  intent(in) :: tau
    integer,                  intent(in) :: taup1
    real,                     intent(in) :: dtts
  end subroutine g_tracer_diag

  subroutine g_tracer_traverse(g_tracer_list)
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
  end subroutine g_tracer_traverse

  subroutine g_tracer_get_name(g_tracer,string)
    type(g_tracer_type),    pointer    :: g_tracer
    character(len=*),        intent(out) :: string
  end subroutine g_tracer_get_name

  subroutine g_tracer_get_alias(g_tracer,string)
    type(g_tracer_type),    pointer    :: g_tracer
    character(len=*),        intent(out) :: string
  end subroutine g_tracer_get_alias

  function g_tracer_is_prog(g_tracer)
    logical :: g_tracer_is_prog
    type(g_tracer_type),    pointer    :: g_tracer
  end function g_tracer_is_prog

  subroutine g_tracer_get_next(g_tracer,g_tracer_next)
    type(g_tracer_type),    pointer    :: g_tracer,g_tracer_next
  end subroutine g_tracer_get_next

  subroutine g_tracer_vertdiff_G(g_tracer, h_old, ea, eb, dt, kg_m2_to_H, m_to_H, tau, mom)
    type(g_tracer_type),    pointer  :: g_tracer
    real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,:), intent(in) :: h_old, ea, eb
    real,                   intent(in) :: dt, kg_m2_to_H, m_to_H
    integer,                intent(in) :: tau
    logical,                                                intent(in), optional :: mom
  end subroutine g_tracer_vertdiff_G

  subroutine g_tracer_vertdiff_M(g_tracer,dh, dhw, diff_cbt, dt, rho0,tau)
    type(g_tracer_type),    pointer  :: g_tracer
    real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,:), intent(in) :: dh, diff_cbt
    real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,0:), intent(in) :: dhw
    real,                   intent(in) :: dt,rho0
    integer,                intent(in) :: tau
  end subroutine g_tracer_vertdiff_M

  subroutine g_diag_field_add(node_ptr, diag_id, package_name, name, axes, init_time, longname, units, &
                            missing_value, Z_diag, field_ptr, Zname, Zlongname, Zunits)
    type(g_diag_type), pointer :: node_ptr
    integer, intent(inout) :: diag_id
    character(len=*), intent(in) :: package_name, name
    integer, intent(in) :: axes(:)
    type(time_type), intent(in) :: init_time
    character(len=*), intent(in) :: longname, units
    real, optional, intent(in) :: missing_value
    integer, optional, intent(in) :: Z_diag
    character(len=*), optional, intent(in) :: Zname, Zlongname, Zunits
    real, optional, pointer :: field_ptr(:,:,:)
  end subroutine g_diag_field_add

  subroutine g_tracer_print_info(g_tracer_list)
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
    integer               :: num_prog,num_diag
  end subroutine g_tracer_print_info

  subroutine g_tracer_get_src_info(g_tracer_list,name,src_file, src_var_name, src_var_unit, src_var_gridspec,&
                                   src_var_record, src_var_valid_min, src_var_valid_max)
    type(g_tracer_type),      pointer    :: g_tracer_list,g_tracer
    character(len=*),         intent(in) :: name
    character(len=*),         intent(out):: src_file, src_var_name, src_var_unit, src_var_gridspec
    integer,                  intent(out):: src_var_record
    real,                     intent(out):: src_var_valid_min, src_var_valid_max
  end subroutine g_tracer_get_src_info

  function g_register_diag_field(module_name, field_name, axes, init_time,         &
       long_name, units, missing_value, range, mask_variant, standard_name,      &
       verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
       cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
       x_cell_method, y_cell_method, v_cell_method, diag_CS)
    integer :: g_register_diag_field !< An integer handle for a diagnostic array.
    character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model" or "ice_shelf_model"
    character(len=*), intent(in) :: field_name !< Name of the diagnostic field
    type(time_type),intent(in)  :: init_time !< Time at which a field is first available?
    type(g_diag_ctrl),optional, pointer :: diag_CS
    integer,          optional, intent(in) :: axes(:)
    character(len=*), optional, intent(in) :: long_name !< Long name of a field.
    character(len=*), optional, intent(in) :: units !< Units of a field.
    character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
    real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
    real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
    logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with post_data calls (not used in MOM?)
    logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
    logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
    character(len=*), optional, intent(out):: err_msg !< String into which an error message might be placed (not used in MOM?)
    character(len=*), optional, intent(in) :: interp_method !< no clue (not used in MOM?)
    integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
    character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
    character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
    character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
    character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
    character(len=*), optional, intent(in) :: cell_methods !< String to append as cell_methods attribute. Use '' to have no attribute.
    !! If present, this overrides the default constructed from the default for
    !! each individual axis direction.
    character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction. Use '' have no method.
    character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction. Use '' have no method.
    character(len=*), optional, intent(in) :: v_cell_method !< Specifies the cell method for the vertical direction. Use '' have no method.
  end function g_register_diag_field

  logical function g_send_data_0d(diag_field_id, field, time, err_msg, diag_CS)
    integer, intent(in) :: diag_field_id
    real, intent(in) :: field
    type(time_type), intent(in), optional :: time
    character(len=*), intent(out), optional :: err_msg
    type(g_diag_ctrl),optional, pointer :: diag_CS
  end function g_send_data_0d

  logical function g_send_data_1d(diag_field_id, field, time, is_in, mask, rmask, ie_in, weight, err_msg, diag_CS)
    integer, intent(in) :: diag_field_id
    real, dimension(:), intent(in) :: field
    real, intent(in), optional :: weight
    real, intent(in), dimension(:), optional :: rmask
    type (time_type), intent(in), optional :: time
    integer, intent(in), optional :: is_in, ie_in
    logical, intent(in), dimension(:), optional :: mask
    character(len=*), intent(out), optional :: err_msg
    type(g_diag_ctrl),optional, pointer :: diag_CS
  end function g_send_data_1d

  logical function g_send_data_2d(diag_field_id, field, time, is_in, js_in, &
       & mask, rmask, ie_in, je_in, weight, err_msg, diag_CS)
    integer, intent(in) :: diag_field_id
    real, intent(in), dimension(:,:) :: field
    real, intent(in), optional :: weight
    type (time_type), intent(in), optional :: time
    integer, intent(in), optional :: is_in, js_in, ie_in, je_in
    logical, intent(in), dimension(:,:), optional :: mask
    real, intent(in), dimension(:,:),optional :: rmask
    character(len=*), intent(out), optional :: err_msg
    type(g_diag_ctrl),optional, pointer :: diag_CS
  end function g_send_data_2d

  logical function g_send_data_3d(diag_field_id, field, time, is_in, js_in, ks_in, &
             & mask, rmask, ie_in, je_in, ke_in, weight, err_msg, diag_CS)
    integer, intent(in) :: diag_field_id
    real, dimension(:,:,:), intent(in) :: field
    real, intent(in), optional :: weight
    type (time_type), intent(in), optional :: time
    integer, intent(in), optional :: is_in, js_in, ks_in,ie_in,je_in, ke_in
    logical, dimension(:,:,:), intent(in), optional :: mask
    real, dimension(:,:,:), intent(in), optional :: rmask
    character(len=*), intent(out), optional :: err_msg
    type(g_diag_ctrl),optional, pointer :: diag_CS
  end function g_send_data_3d

end module g_tracer_utils
