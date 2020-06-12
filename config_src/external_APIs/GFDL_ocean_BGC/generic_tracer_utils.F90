!> g_tracer_utils module consists of core utility subroutines to be used by
!! all generic tracer modules.  These include the lowest level functions
!! for adding, allocating memory, and record keeping of individual generic
!! tracers irrespective of their physical/chemical nature.
module g_tracer_utils
#include <fms_platform.h>

  use coupler_types_mod, only: coupler_2d_bc_type
  use time_manager_mod, only : time_type
  use field_manager_mod, only: fm_string_len
  use MOM_diag_mediator, only : g_diag_ctrl=>diag_ctrl

implicit none ; private

  !> Each generic tracer node is an instant of a FORTRAN type with the following member variables.
  !! These member fields are supposed to uniquely define an individual tracer.
  !! One such type shall be instantiated for EACH individual tracer.
  type g_tracer_type
    !> Tracer concentration field in space (and time)
    !! MOM keeps the prognostic tracer fields at 3 time levels, hence 4D.
    real, pointer, dimension(:,:,:,:) :: field  => NULL()
    !> Tracer concentration in river runoff
    real, _ALLOCATABLE, dimension(:,:) :: trunoff _NULL
    logical :: requires_restart = .true. !< Unknown
    !> Tracer source: filename, type, var name, units, record, gridfile
    character(len=fm_string_len) :: src_file, src_var_name, src_var_unit, src_var_gridspec
    integer :: src_var_record !< Unknown
    logical :: requires_src_info = .false. !< Unknown
    real    :: src_var_unit_conversion = 1.0 !< This factor depends on the tracer. Ask Jasmin
    real    :: src_var_valid_min = 0.0 !< Unknown
  end type g_tracer_type

  !> Unknown
  type g_diag_type
    integer :: dummy !< A dummy member, not part of the API
  end type g_diag_type

  !> The following type fields are common to ALL generic tracers and hence has to be instantiated only once
  type g_tracer_common
!   type(g_diag_ctrl) :: diag_CS !< Unknown
    !> Domain extents
    integer :: isd,jsd
  end type g_tracer_common

  !> Unknown dangerous module data!
  type(g_tracer_common), target, save :: g_tracer_com

  public :: g_tracer_type
! public :: g_tracer_find
! public :: g_tracer_add
! public :: g_tracer_init
  public :: g_tracer_flux_init
! public :: g_tracer_column_int
! public :: g_tracer_flux_at_depth
! public :: g_tracer_add_param
  public :: g_tracer_set_values
  public :: g_tracer_get_values
  public :: g_tracer_get_pointer
  public :: g_tracer_get_common
  public :: g_tracer_set_common
  public :: g_tracer_set_csdiag
! public :: g_tracer_set_files
! public :: g_tracer_coupler_set
! public :: g_tracer_coupler_get
  public :: g_tracer_send_diag
! public :: g_tracer_diag
  public :: g_tracer_get_name
  public :: g_tracer_get_alias
  public :: g_tracer_get_next
! public :: g_tracer_register_diag
  public :: g_tracer_is_prog
! public :: g_tracer_vertdiff_G
! public :: g_tracer_vertdiff_M
! public :: g_tracer_start_param_list
! public :: g_tracer_end_param_list
  public :: g_diag_type
! public :: g_diag_field_add
! public :: g_tracer_set_pointer
! public :: g_tracer_print_info
! public :: g_tracer_coupler_accumulate
! public :: g_tracer_get_src_info
! public :: g_register_diag_field
! public :: g_send_data

! !> Add a new parameter for the generic tracer package
! !!
! !!  This subroutine is used to add a new parameter by the calling tracer package.
! !!  It provides a mechanism for parameter overwrite through the field_table.
! !! For each tracer package there is a field called namelists and there
! !! the parameters can be modified from their value set by this method.
! !! E.g., we may have the following in the field_table
! !!
! !! "namelists","ocean_mod","generic_topaz"
! !! init = t
! !!  /
! !!
! !! This will overwrite the parameter topaz%init to be .true. at the run time
! !! even though generic_topaz package had in the code
! !! `call g_tracer_add_param('init', topaz%init, .false. )`
! !!
! !! For the parameters overwrite mechanism to work all calls
! !! for adding new parameters (refer to description for subroutine g_tracer_add_param)
! !! should happen between a `call g_tracer_start_param_list(package_name)`
! !! and a `call g_tracer_end_param_list(package_name)`
! interface g_tracer_add_param
!   module procedure g_tracer_add_param_real
!   module procedure g_tracer_add_param_logical
!   module procedure g_tracer_add_param_integer
!   module procedure g_tracer_add_param_string
! end interface

! !> Unknown
! interface g_tracer_set_pointer
!   module procedure g_tracer_set_pointer_3d
!   module procedure g_tracer_set_pointer_4d
! end interface g_tracer_set_pointer

! !> Unknown
! interface g_send_data
!   module procedure g_send_data_0d
!   module procedure g_send_data_1d
!   module procedure g_send_data_2d
!   module procedure g_send_data_3d
! end interface

  !> Set the values of various (array) members of the tracer node g_tracer_type
  !!
  !! This function is overloaded to set the values of the following member variables
  interface g_tracer_set_values
    module procedure g_tracer_set_real
    module procedure g_tracer_set_2D
    module procedure g_tracer_set_3D
    module procedure g_tracer_set_4D
  end interface

  !> Reverse of interface g_tracer_set_values for getting the tracer member arrays  in the argument value
  !!
  !! This means "get the values of array %field_name for tracer tracer_name and put them in argument array_out"
  interface g_tracer_get_values
    module procedure g_tracer_get_4D_val
    module procedure g_tracer_get_3D_val
    module procedure g_tracer_get_2D_val
    module procedure g_tracer_get_real
    module procedure g_tracer_get_string
  end interface

  !> Return the pointer to the requested field of a particular tracer
  !!
  !! This means "get the pointer of array %field_name for tracer tracer_name in argument array_ptr"
  interface g_tracer_get_pointer
    module procedure g_tracer_get_4D
    module procedure g_tracer_get_3D
    module procedure g_tracer_get_2D
  end interface

contains

! !> Mark the start of adding new parameters for a package
! !! For the parameters override mechanism to work all calls
! !! for adding new parameters (refer to description for subroutine g_tracer_add_param)
! !! should happen between a `call g_tracer_start_param_list(package_name)`
! !! and a `call g_tracer_end_param_list(package_name)`
! subroutine g_tracer_start_param_list(package_name)
!   !> Name of the generic tracer package that is adding the parameters (e.g., "generic_cfc")
!   character(len=fm_string_len), intent(in) :: package_name
! end subroutine g_tracer_start_param_list

! !> Mark the start of adding new parameters for a package
! subroutine g_tracer_end_param_list(package_name)
!   !> Name of the generic tracer package that is adding the parameters (e.g., "generic_cfc")
!   character(len=fm_string_len) :: package_name
! end subroutine g_tracer_end_param_list

! !> Overload interface g_tracer_add_param for real parameter
! subroutine g_tracer_add_param_real(name, var,  value)
!   character(len=*), intent(in)  :: name !< Unknown
!   real,             intent(in)  :: value !< Unknown
!   real,             intent(out) :: var !< Unknown
! end subroutine g_tracer_add_param_real

! !> Overload interface g_tracer_add_param for logical parameter
! subroutine g_tracer_add_param_logical(name, var,  value)
!   character(len=*), intent(in)  :: name !< Unknown
!   logical,          intent(in)  :: value !< Unknown
!   logical,          intent(out) :: var !< Unknown
! end subroutine g_tracer_add_param_logical

! !> Overload interface g_tracer_add_param for integer parameter
! subroutine g_tracer_add_param_integer(name, var,  value)
!   character(len=*), intent(in)  :: name !< Unknown
!   integer,          intent(in)  :: value !< Unknown
!   integer,          intent(out) :: var !< Unknown
! end subroutine g_tracer_add_param_integer

! !> Overload interface g_tracer_add_param for string parameter
! subroutine g_tracer_add_param_string(name, var,  value)
!   character(len=*), intent(in)  :: name !< Unknown
!   character(len=*), intent(in)  :: value !< Unknown
!   character(len=*), intent(out) :: var !< Unknown
! end subroutine g_tracer_add_param_string

! !> Add a new tracer (node) at the top of the list of generic tracers
! !! This subroutine call adds an individual new tracer to the growing list of generic tracers.
! !! It then allocates all the necessary arrays for using this tracer in the Ocean model that requested it.
! !! The information passed into this subroutine should be enough to fully describe the individual tracer
! subroutine g_tracer_add(node_ptr, package, name, longname, units,  prog, const_init_value,init_value,&
!      flux_gas, flux_gas_name, flux_runoff, flux_wetdep, flux_drydep, flux_gas_molwt, flux_gas_param, &
!      flux_param, flux_bottom, btm_reservoir, move_vertical, diff_vertical, sink_rate, flux_gas_restart_file, &
!      flux_gas_type, requires_src_info, standard_name,diag_name, diag_field_units,diag_field_scaling_factor, &
!      implementation)
!   !> Pointer to the head node of the tracer list. This is also going to be the pointer to the node being added after the call
!   type(g_tracer_type), pointer :: node_ptr
!   character(len=*),   intent(in) :: package !< Name of tracer package adding this node
!   character(len=*),   intent(in) :: name !< Name of this tracer
!   character(len=*),   intent(in) :: longname !< Descriptive name of this tracer
!   character(len=*),   intent(in) :: units !< Concentration units (units of array %field)
!   logical,            intent(in) :: prog !< .true. for prognastic , .false. for diagnostic tracer
!   real,               intent(in), optional :: const_init_value !< Initial value of concenteration if constant
!   real,               intent(in), optional :: init_value !< Unknown
!   real,               intent(in), optional :: sink_rate !< Sinking rate if non-zero
!   logical,            intent(in), optional :: flux_gas !< .true. if there is gas flux exchange with atmos
!   logical,            intent(in), optional :: flux_runoff !< .true. if there is runoff flux
!   logical,            intent(in), optional :: flux_wetdep !< .true. if there is wetdep flux
!   logical,            intent(in), optional :: flux_drydep !< .true. if there is drydep flux
!   logical,            intent(in), optional :: flux_bottom !< .true. if there is bottom flux
!   logical,            intent(in), optional :: btm_reservoir !< .true. if there is bottom reservoir
!   logical,            intent(in), optional :: move_vertical !< .true. if there is active vertical movement
!   logical,            intent(in), optional :: diff_vertical !< Unknown
!   real,               intent(in), optional :: flux_gas_molwt !< Unknown
!   real, dimension(:), intent(in), optional :: flux_gas_param !< Array of parameters for gas flux (refer to documentation for subroutine aof_set_coupler_flux() )
!   real, dimension(:), intent(in), optional :: flux_param !< Array of parameters for non-gas flux (refer to documentation for subroutine aof_set_coupler_flux() )
!   character(len=*),   intent(in), optional :: flux_gas_name !< Name of the atmospheric tracer to exchange flux with (if flux_gas=.true.)
!   character(len=*),   intent(in), optional :: implementation !< Unknown
!   character(len=*),   intent(in), optional :: flux_gas_type !< Unknown
!   character(len=*),   intent(in), optional :: flux_gas_restart_file !< Unknown
!   logical,            intent(in), optional :: requires_src_info !< Unknown
!   character(len=*),   intent(in), optional :: standard_name !< Unknown
!   character(len=*),   intent(in), optional :: diag_name !< Unknown
!   character(len=*),   intent(in), optional :: diag_field_units !< Unknown
!   real,               intent(in), optional :: diag_field_scaling_factor !< Unknown
! end subroutine g_tracer_add

! !> Unknown
! subroutine g_tracer_init(g_tracer)
!   type(g_tracer_type), pointer :: g_tracer !< Pointer to this tracer node
! end subroutine g_tracer_init

  !> Unknown
  subroutine g_tracer_flux_init(g_tracer)
    type(g_tracer_type), pointer :: g_tracer !< Pointer to this tracer node
  end subroutine g_tracer_flux_init

! !> Diag-register all the internal fields that were _ALLOCATED for a tracer
! !!
! !! Use diag_manager register_diag_field for each of the field arrays that were _ALLOCATED for a tracer node.
! !! These include %field,  %tendency, %stf, %stf_gas, %deltap, %kw, %btf, %trunoff, %alpha, %csurf, %sc_no, %btm_reservoir.
! subroutine g_tracer_register_diag(g_tracer)
!   type(g_tracer_type), pointer :: g_tracer !< Pointer to this tracer node
! end subroutine g_tracer_register_diag

! !> Set coupler values only for tracers that have _ALLOCATED %alpha, %csurf and %sc_no
! !!
! !! Use coupler_util subroutine set_coupler_values() to set the coupler values
! !! for fluxes to be exchanged with Ice for the requested fluxes.
! !! NOTE:
! !! This is a collective subroutine and will transverse the list of generic tracers and 
! !! set the coupler values for each tracer node accordingly.
! subroutine g_tracer_coupler_set(g_tracer_list,IOB_struc,value)
!   type(g_tracer_type), pointer :: g_tracer_list !< Pointer to the head of the generic tracer list
!   type(g_tracer_type), pointer :: g_tracer !< Pointer to this tracer node
!   type(coupler_2d_bc_type), intent(inout) :: IOB_struc !< The coupler flux IOB structure
!   real, optional :: value !< Set the coupler values to a constant (particularly 0) is desired
! end subroutine g_tracer_coupler_set

! !> Get coupler values only for tracers that have _ALLOCATED arrays for the fluxes
! !!
! !! Use coupler_util subroutine extract_coupler_values() to get the coupler values
! !! for fluxes to be exchanged with Ice for the requested fluxes only.
! !! NOTE:
! !! This is a collective subroutine and will transverse the list of generic tracers and 
! !! get the coupler values for each tracer node accordingly
! subroutine g_tracer_coupler_get(g_tracer_list,IOB_struc, weight, model_time)
!   type(g_tracer_type),         pointer    :: g_tracer_list !< Pointer to the head of the generic tracer list
!   type(g_tracer_type),         pointer    :: g_tracer !< Pointer to this tracer node
!   type(coupler_2d_bc_type),    intent(in) :: IOB_struc !< The coupler flux IOB structure
!   real,               optional,intent(in) :: weight !< Unknown
!   type(time_type),    optional,intent(in) :: model_time !< Time
! end subroutine g_tracer_coupler_get

! !> Unknown
! subroutine g_tracer_coupler_accumulate(g_tracer_list,IOB_struc, weight, model_time)
!   type(g_tracer_type),         pointer    :: g_tracer_list !< Pointer to the head of the generic tracer list
!   type(g_tracer_type),         pointer    :: g_tracer !< Pointer to this tracer node
!   type(coupler_2d_bc_type),    intent(in) :: IOB_struc !< The coupler flux IOB structure
!   real,               optional,intent(in) :: weight !< Unknown
!   type(time_type),    optional,intent(in) :: model_time !< Time
! end subroutine g_tracer_coupler_accumulate

  !> Unknown
  subroutine g_tracer_set_csdiag(diag_CS)
    type(g_diag_ctrl),  target,intent(in) :: diag_CS !< Unknown
  end subroutine g_tracer_set_csdiag

  subroutine g_tracer_set_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes,grid_tmask,grid_kmt,init_time)
    integer,                     intent(in) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,axes(3) !< Unknown
    real, dimension(isd:,jsd:,:),intent(in) :: grid_tmask !< Unknown
    integer,dimension(isd:,jsd:),intent(in) :: grid_kmt !< Unknown
    type(time_type),             intent(in) :: init_time !< Unknown
  end subroutine g_tracer_set_common

  subroutine g_tracer_get_common(isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau,&
       axes,grid_tmask,grid_mask_coast,grid_kmt,init_time,diag_CS)
    integer,               intent(out) :: isc,iec,jsc,jec,isd,ied,jsd,jed,nk,ntau !< Unknown
    integer,optional,      intent(out) :: axes(3) !< Unknown
    type(time_type), optional,      intent(out) :: init_time !< Unknown
    real, optional, dimension(:,:,:),pointer    :: grid_tmask !< Unknown
    integer, optional, dimension(:,:),  pointer :: grid_mask_coast !< Unknown
    integer, optional, dimension(:,:),  pointer :: grid_kmt !< Unknown
    type(g_diag_ctrl), optional,        pointer :: diag_CS !< Unknown
  end subroutine g_tracer_get_common

! subroutine g_tracer_get_diagCS(diag_CS)
!   type(g_diag_ctrl),        pointer :: diag_CS
! end subroutine g_tracer_get_diagCS

! subroutine g_tracer_set_files(ice_restart_file,ocean_restart_file)
!   character(len=*),   intent(in) :: ice_restart_file
!   character(len=*),   intent(in) :: ocean_restart_file
! end subroutine g_tracer_set_files

  !> Unknown
  subroutine g_tracer_get_4D(g_tracer_list,name,member,array_ptr)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    real, dimension(:,:,:,:), pointer    :: array_ptr
  end subroutine g_tracer_get_4D

  !> Unknown
  subroutine g_tracer_get_3D(g_tracer_list,name,member,array_ptr)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    real, dimension(:,:,:), pointer    :: array_ptr !< Unknown
  end subroutine g_tracer_get_3D

  !> Unknown
  subroutine g_tracer_get_2D(g_tracer_list,name,member,array_ptr)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    real, dimension(:,:), pointer    :: array_ptr !< Unknown
  end subroutine g_tracer_get_2D

  !> Unknown
  subroutine g_tracer_get_4D_val(g_tracer_list,name,member,array,isd,jsd)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    integer,                  intent(in) :: isd,jsd !< Unknown
    real, dimension(isd:,jsd:,:,:), intent(out):: array !< Unknown
  end subroutine g_tracer_get_4D_val

  !> Unknown
  subroutine g_tracer_get_3D_val(g_tracer_list,name,member,array,isd,jsd,ntau,positive)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    integer,                  intent(in) :: isd,jsd !< Unknown
    integer, optional,        intent(in) :: ntau !< Unknown
    logical, optional,        intent(in) :: positive !< Unknown
    real, dimension(isd:,jsd:,:), intent(out):: array !< Unknown
    integer :: tau
    character(len=fm_string_len), parameter :: sub_name = 'g_tracer_get_3D_val'
  end subroutine g_tracer_get_3D_val

  !> Unknown
  subroutine g_tracer_get_2D_val(g_tracer_list,name,member,array,isd,jsd)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    integer,                  intent(in) :: isd,jsd !< Unknown
    real, dimension(isd:,jsd:), intent(out):: array !< Unknown
  end subroutine g_tracer_get_2D_val

  !> Unknown
  subroutine g_tracer_get_real(g_tracer_list,name,member,value)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    real,                     intent(out):: value
  end subroutine g_tracer_get_real

  !> Unknown
  subroutine g_tracer_get_string(g_tracer_list,name,member,string)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    character(len=fm_string_len), intent(out) :: string !< Unknown
  end subroutine g_tracer_get_string

  !> Unknown
  subroutine g_tracer_set_2D(g_tracer_list,name,member,array,isd,jsd,weight)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),      pointer    :: g_tracer_list, g_tracer !< Unknown
    integer,                   intent(in) :: isd,jsd !< Unknown
    real, dimension(isd:,jsd:),intent(in) :: array !< Unknown
    real, optional            ,intent(in) :: weight !< Unknown
  end subroutine g_tracer_set_2D

  !> Unknown
  subroutine g_tracer_set_3D(g_tracer_list,name,member,array,isd,jsd,ntau)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    integer,                  intent(in) :: isd,jsd !< Unknown
    integer, optional,        intent(in) :: ntau !< Unknown
    real, dimension(isd:,jsd:,:), intent(in)       :: array !< Unknown
  end subroutine g_tracer_set_3D

  !> Unknown
  subroutine g_tracer_set_4D(g_tracer_list,name,member,array,isd,jsd)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    integer,                  intent(in) :: isd,jsd !< Unknown
    real, dimension(isd:,jsd:,:,:), intent(in)       :: array !< Unknown
  end subroutine g_tracer_set_4D

  !> Unknown
  subroutine g_tracer_set_real(g_tracer_list,name,member,value)
    character(len=*),         intent(in) :: name !< Unknown
    character(len=*),         intent(in) :: member !< Unknown
    type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer !< Unknown
    real,                     intent(in) :: value !< Unknown
  end subroutine g_tracer_set_real

! subroutine g_tracer_set_pointer_4D(g_tracer_list,name,member,array,ilb,jlb)
!   character(len=*),               intent(in) :: name
!   character(len=*),               intent(in) :: member
!   type(g_tracer_type),            pointer    :: g_tracer_list, g_tracer
!   integer,                        intent(in) :: ilb,jlb
!   real, dimension(ilb:,jlb:,:,:), target, intent(in) :: array
! end subroutine g_tracer_set_pointer_4D

! subroutine g_tracer_set_pointer_3D(g_tracer_list,name,member,array,ilb,jlb)
!   character(len=*),               intent(in) :: name
!   character(len=*),               intent(in) :: member
!   type(g_tracer_type),            pointer    :: g_tracer_list, g_tracer
!   integer,                        intent(in) :: ilb,jlb
!   real, dimension(ilb:,jlb:,:), target, intent(in) :: array
! end subroutine g_tracer_set_pointer_3D

! subroutine g_tracer_find(g_tracer,name)
!   character(len=*),         intent(in) :: name
!   type(g_tracer_type),    pointer    :: g_tracer
! end subroutine g_tracer_find

! subroutine g_tracer_column_int(depth, ilb, jlb, var, dzt, rho_dzt, rd, k_level, integral, caller)
!   real,                         intent(in)            :: depth
!   integer,                      intent(in)            :: ilb
!   integer,                      intent(in)            :: jlb
!   real, dimension(ilb:,jlb:,:), intent(in)            :: var
!   real, dimension(ilb:,jlb:,:), intent(in)            :: dzt
!   real, dimension(ilb:,jlb:,:), intent(in)            :: rho_dzt
!   real, dimension(ilb:,jlb:,:), intent(inout)         :: rd
!   integer,                      intent(inout)         :: k_level
!   real, dimension(ilb:,jlb:),   intent(out)           :: integral
!   character(len=*),             intent(in), optional  :: caller
! end subroutine g_tracer_column_int

! subroutine g_tracer_flux_at_depth(depth, ilb, jlb, var, dzt, k_level, frac, initialized, flux, caller)
!   real,                            intent(in)                 :: depth
!   integer,                         intent(in)                 :: ilb
!   integer,                         intent(in)                 :: jlb
!   real,    dimension(ilb:,jlb:,:), intent(in)                 :: var
!   real,    dimension(ilb:,jlb:,:), intent(in)                 :: dzt
!   integer, dimension(ilb:,jlb:),   intent(inout)              :: k_level
!   real,    dimension(ilb:,jlb:),   intent(inout)              :: frac
!   logical,                         intent(inout)              :: initialized
!   real,    dimension(ilb:,jlb:),   intent(out)                :: flux
!   character(len=*),                intent(in),    optional    :: caller
! end subroutine g_tracer_flux_at_depth

  subroutine g_tracer_send_diag(g_tracer_list,model_time,tau)
    type(g_tracer_type), pointer    :: g_tracer_list !< pointer to the head of the generic tracer list
    type(g_tracer_type), pointer    :: g_tracer !< Pointer to tracer node
    type(time_type),     intent(in) :: model_time !< Time
    integer,             intent(in) :: tau !< The time step for the %field 4D field to be reported
  end subroutine g_tracer_send_diag

! subroutine g_tracer_diag(g_tracer_list, ilb, jlb, rho_dzt_tau, rho_dzt_taup1, model_time, tau, taup1, dtts)
!   type(g_tracer_type),    pointer    :: g_tracer_list
!   integer,                  intent(in) :: ilb
!   integer,                  intent(in) :: jlb
!   real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt_tau
!   real, dimension(ilb:,jlb:,:),   intent(in) :: rho_dzt_taup1
!   type(time_type),        intent(in) :: model_time
!   integer,                  intent(in) :: tau
!   integer,                  intent(in) :: taup1
!   real,                     intent(in) :: dtts
! end subroutine g_tracer_diag

! subroutine g_tracer_traverse(g_tracer_list)
!   type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
! end subroutine g_tracer_traverse

  !> Unknown
  subroutine g_tracer_get_name(g_tracer,string)
    type(g_tracer_type),    pointer    :: g_tracer !< Unknown
    character(len=*),        intent(out) :: string !< Unknown
  end subroutine g_tracer_get_name

  !> Unknown
  subroutine g_tracer_get_alias(g_tracer,string)
    type(g_tracer_type), pointer  :: g_tracer !< Unknown
    character(len=*), intent(out) :: string !< Unknown
  end subroutine g_tracer_get_alias

  !> Is the tracer prognostic?
  function g_tracer_is_prog(g_tracer)
    logical :: g_tracer_is_prog
    type(g_tracer_type), pointer :: g_tracer !< Pointer to tracer node
  end function g_tracer_is_prog

  !> get the next tracer in the list
  subroutine g_tracer_get_next(g_tracer,g_tracer_next)
    type(g_tracer_type), pointer :: g_tracer !< Pointer to tracer node
    type(g_tracer_type), pointer :: g_tracer_next !< Pointer to the next tracer node in the list
  end subroutine g_tracer_get_next

  !>Vertical Diffusion of a tracer node
  !!
  !! This subroutine solves a tridiagonal equation to find and set values of vertically diffused field for a tracer node.
  !! This is ported from GOLD (vertdiff) and simplified
  !! Since the surface flux from the atmosphere (%stf) has the units of mol/m^2/sec the resulting tracer concentration
  !! has units of mol/Kg
  subroutine g_tracer_vertdiff_G(g_tracer, h_old, ea, eb, dt, kg_m2_to_H, m_to_H, tau, mom)
    type(g_tracer_type),    pointer  :: g_tracer
    !> Layer thickness before entrainment, in m or kg m-2.
    real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,:), intent(in) :: h_old
    !> The amount of fluid entrained from the layer above, in H.
    real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,:), intent(in) :: ea
    !> The amount of fluid entrained from the layer below, in H.
    real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,:), intent(in) :: eb
    real,     intent(in) :: dt !< The amount of time covered by this call, in s.
    real,     intent(in) :: kg_m2_to_H !< A conversion factor that translates kg m-2 into
                                       !! the units of h_old (H)
    real,     intent(in) :: m_to_H !< A conversion factor that translates m into the units
                                   !! of h_old (H).
    integer,  intent(in) :: tau !< Unknown
    logical,  intent(in), optional :: mom
  end subroutine g_tracer_vertdiff_G

! subroutine g_tracer_vertdiff_M(g_tracer,dh, dhw, diff_cbt, dt, rho0,tau)
!   type(g_tracer_type),    pointer  :: g_tracer
!   real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,:), intent(in) :: dh, diff_cbt
!   real, dimension(g_tracer_com%isd:,g_tracer_com%jsd:,0:), intent(in) :: dhw
!   real,                   intent(in) :: dt,rho0
!   integer,                intent(in) :: tau
! end subroutine g_tracer_vertdiff_M

! subroutine g_diag_field_add(node_ptr, diag_id, package_name, name, axes, init_time, longname, units, &
!                           missing_value, Z_diag, field_ptr, Zname, Zlongname, Zunits)
!   type(g_diag_type), pointer :: node_ptr
!   integer, intent(inout) :: diag_id
!   character(len=*), intent(in) :: package_name, name
!   integer, intent(in) :: axes(:)
!   type(time_type), intent(in) :: init_time
!   character(len=*), intent(in) :: longname, units
!   real, optional, intent(in) :: missing_value
!   integer, optional, intent(in) :: Z_diag
!   character(len=*), optional, intent(in) :: Zname, Zlongname, Zunits
!   real, optional, pointer :: field_ptr(:,:,:)
! end subroutine g_diag_field_add

! subroutine g_tracer_print_info(g_tracer_list)
!   type(g_tracer_type),    pointer    :: g_tracer_list, g_tracer
!   integer               :: num_prog,num_diag
! end subroutine g_tracer_print_info

! subroutine g_tracer_get_src_info(g_tracer_list,name,src_file, src_var_name, src_var_unit, src_var_gridspec,&
!                                  src_var_record, src_var_valid_min, src_var_valid_max)
!   type(g_tracer_type),      pointer    :: g_tracer_list,g_tracer
!   character(len=*),         intent(in) :: name
!   character(len=*),         intent(out):: src_file, src_var_name, src_var_unit, src_var_gridspec
!   integer,                  intent(out):: src_var_record
!   real,                     intent(out):: src_var_valid_min, src_var_valid_max
! end subroutine g_tracer_get_src_info

! function g_register_diag_field(module_name, field_name, axes, init_time,         &
!      long_name, units, missing_value, range, mask_variant, standard_name,      &
!      verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
!      cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
!      x_cell_method, y_cell_method, v_cell_method, diag_CS)
!   integer :: g_register_diag_field !< An integer handle for a diagnostic array.
!   character(len=*), intent(in) :: module_name !< Name of this module, usually "ocean_model" or "ice_shelf_model"
!   character(len=*), intent(in) :: field_name !< Name of the diagnostic field
!   type(time_type),intent(in)  :: init_time !< Time at which a field is first available?
!   type(g_diag_ctrl),optional, pointer :: diag_CS
!   integer,          optional, intent(in) :: axes(:)
!   character(len=*), optional, intent(in) :: long_name !< Long name of a field.
!   character(len=*), optional, intent(in) :: units !< Units of a field.
!   character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
!   real,             optional, intent(in) :: missing_value !< A value that indicates missing values.
!   real,             optional, intent(in) :: range(2) !< Valid range of a variable (not used in MOM?)
!   logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with post_data calls (not used in MOM?)
!   logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in MOM?)
!   logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in MOM?)
!   character(len=*), optional, intent(out):: err_msg !< String into which an error message might be placed (not used in MOM?)
!   character(len=*), optional, intent(in) :: interp_method !< no clue (not used in MOM?)
!   integer,          optional, intent(in) :: tile_count !< no clue (not used in MOM?)
!   character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
!   character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
!   character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
!   character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
!   character(len=*), optional, intent(in) :: cell_methods !< String to append as cell_methods attribute. Use '' to have no attribute.
!   !! If present, this overrides the default constructed from the default for
!   !! each individual axis direction.
!   character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction. Use '' have no method.
!   character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction. Use '' have no method.
!   character(len=*), optional, intent(in) :: v_cell_method !< Specifies the cell method for the vertical direction. Use '' have no method.
! end function g_register_diag_field

! logical function g_send_data_0d(diag_field_id, field, time, err_msg, diag_CS)
!   integer, intent(in) :: diag_field_id
!   real, intent(in) :: field
!   type(time_type), intent(in), optional :: time
!   character(len=*), intent(out), optional :: err_msg
!   type(g_diag_ctrl),optional, pointer :: diag_CS
! end function g_send_data_0d

! logical function g_send_data_1d(diag_field_id, field, time, is_in, mask, rmask, ie_in, weight, err_msg, diag_CS)
!   integer, intent(in) :: diag_field_id
!   real, dimension(:), intent(in) :: field
!   real, intent(in), optional :: weight
!   real, intent(in), dimension(:), optional :: rmask
!   type (time_type), intent(in), optional :: time
!   integer, intent(in), optional :: is_in, ie_in
!   logical, intent(in), dimension(:), optional :: mask
!   character(len=*), intent(out), optional :: err_msg
!   type(g_diag_ctrl),optional, pointer :: diag_CS
! end function g_send_data_1d

! logical function g_send_data_2d(diag_field_id, field, time, is_in, js_in, &
!      & mask, rmask, ie_in, je_in, weight, err_msg, diag_CS)
!   integer, intent(in) :: diag_field_id
!   real, intent(in), dimension(:,:) :: field
!   real, intent(in), optional :: weight
!   type (time_type), intent(in), optional :: time
!   integer, intent(in), optional :: is_in, js_in, ie_in, je_in
!   logical, intent(in), dimension(:,:), optional :: mask
!   real, intent(in), dimension(:,:),optional :: rmask
!   character(len=*), intent(out), optional :: err_msg
!   type(g_diag_ctrl),optional, pointer :: diag_CS
! end function g_send_data_2d

! logical function g_send_data_3d(diag_field_id, field, time, is_in, js_in, ks_in, &
!            & mask, rmask, ie_in, je_in, ke_in, weight, err_msg, diag_CS)
!   integer, intent(in) :: diag_field_id
!   real, dimension(:,:,:), intent(in) :: field
!   real, intent(in), optional :: weight
!   type (time_type), intent(in), optional :: time
!   integer, intent(in), optional :: is_in, js_in, ks_in,ie_in,je_in, ke_in
!   logical, dimension(:,:,:), intent(in), optional :: mask
!   real, dimension(:,:,:), intent(in), optional :: rmask
!   character(len=*), intent(out), optional :: err_msg
!   type(g_diag_ctrl),optional, pointer :: diag_CS
! end function g_send_data_3d

end module g_tracer_utils
