!> g_tracer_utils module consists of core utility subroutines to be used by
!! all generic tracer modules.  These include the lowest level functions
!! for adding, allocating memory, and record keeping of individual generic
!! tracers irrespective of their physical/chemical nature.
module g_tracer_utils

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
    real, allocatable, dimension(:,:) :: trunoff
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
  public :: g_tracer_flux_init
  public :: g_tracer_set_values
  public :: g_tracer_get_values
  public :: g_tracer_get_pointer
  public :: g_tracer_get_common
  public :: g_tracer_set_common
  public :: g_tracer_set_csdiag
  public :: g_tracer_send_diag
  public :: g_tracer_get_name
  public :: g_tracer_get_alias
  public :: g_tracer_get_next
  public :: g_tracer_is_prog
  public :: g_diag_type

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

  !> Unknown
  subroutine g_tracer_flux_init(g_tracer)
    type(g_tracer_type), pointer :: g_tracer !< Pointer to this tracer node
  end subroutine g_tracer_flux_init

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

  subroutine g_tracer_send_diag(g_tracer_list,model_time,tau)
    type(g_tracer_type), pointer    :: g_tracer_list !< pointer to the head of the generic tracer list
    type(g_tracer_type), pointer    :: g_tracer !< Pointer to tracer node
    type(time_type),     intent(in) :: model_time !< Time
    integer,             intent(in) :: tau !< The time step for the %field 4D field to be reported
  end subroutine g_tracer_send_diag

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
  !! This subroutine solves a tridiagonal equation to find and set values of vertically diffused field
  !!  for a tracer node.This is ported from GOLD (vertdiff) and simplified
  !! Since the surface flux from the atmosphere (%stf) has the units of mol/m^2/sec the resulting
  !!  tracer concentration has units of mol/Kg
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

end module g_tracer_utils
