!> The MOM6 facility for reading and writing restart files, and querying what has been read.
module MOM_restart

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains, only : pe_here, num_PEs, MOM_domain_type
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_string_functions, only : lowercase
use MOM_grid, only : ocean_grid_type
use MOM_io, only : create_file, fieldtype, file_exists, open_file, close_file
use MOM_io, only : write_field, MOM_read_data, read_data, get_filename_appendix
use MOM_io, only : get_file_info, get_file_atts, get_file_fields, get_file_times
use MOM_io, only : vardesc, var_desc, query_vardesc, modify_vardesc
use MOM_io, only : MULTIPLE, NETCDF_FILE, READONLY_FILE, SINGLE_FILE
use MOM_io, only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_time_manager, only : time_type, time_type_to_real, real_to_time
use MOM_time_manager, only : days_in_month, get_date, set_date
use MOM_verticalGrid, only : verticalGrid_type
use mpp_mod,         only:  mpp_chksum,mpp_pe
use mpp_domains_mod, only : domain1d, domain2d
use mpp_io_mod,      only: mpp_attribute_exist, mpp_get_atts
use mpp_io_mod,      only: axistype
use fms2_io_mod,     only: fms2_register_restart_field => register_restart_field, &
                           fms2_register_axis => register_axis, &
                           fms2_register_field => register_field, &
                           fms2_write_data => write_data, &
                           fms2_read_data => read_data, &
                           fms2_open_file => open_file, &
                           fms2_close_file => close_file, &
                           fms2_get_variable_attribute => get_variable_attribute, &
                           fms2_attribute_exists => variable_att_exists, &
                           fms2_get_variable_names => get_variable_names, &
                           fms2_register_variable_attribute => register_variable_attribute, &
                           fms2_get_dimension_size => get_dimension_size, &
                           fms2_get_num_variables => get_num_variables, &
                           fms2_variable_exists => variable_exists, &
                           fms2_dimension_exists => dimension_exists, &
                           fms2_file_exists => file_exists, &
                           FmsNetcdfDomainFile_t, unlimited
#include <fms_platform.h>
!!                           
implicit none ; private

public restart_init, restart_end, restore_state, register_restart_field
public save_restart, query_initialized, restart_init_end, vardesc
public restart_files_exist, determine_is_new_run, is_new_run

!> A type for making arrays of pointers to 4-d arrays
type p4d
  real, dimension(:,:,:,:), pointer :: p => NULL() !< A pointer to a 4d array
end type p4d

!> A type for making arrays of pointers to 3-d arrays
type p3d
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3d array
end type p3d

!> A type for making arrays of pointers to 2-d arrays
type p2d
  real, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2d array
end type p2d

!> A type for making arrays of pointers to 1-d arrays
type p1d
  real, dimension(:), pointer :: p => NULL() !< A pointer to a 1d array
end type p1d

!> A type for making arrays of pointers to scalars
type p0d
  real, pointer :: p => NULL() !< A pointer to a scalar
end type p0d

!> A structure with information about a single restart field
type field_restart
  type(vardesc) :: vars         !< Description of a field that is to be read from or written
                                !! to the restart file.
  logical :: mand_var           !< If .true. the run will abort if this field is not successfully
                                !! read from the restart file.
  logical :: initialized        !< .true. if this field has been read from the restart file.
  character(len=32) :: var_name !< A name by which a variable may be queried.
end type field_restart

!> A restart registry and the control structure for restarts
type, public :: MOM_restart_CS ; private
  logical :: restart    !< restart is set to .true. if the run has been started from a full restart
                        !! file.  Otherwise some fields must be initialized approximately.
  integer :: novars = 0 !< The number of restart fields that have been registered.
  logical :: parallel_restartfiles  !< If true, each PE writes its own restart file,
                                    !! otherwise they are combined internally.
  logical :: large_file_support     !< If true, NetCDF 3.6 or later is being used
                                    !! and large-file-support is enabled.
  logical :: new_run                !< If true, the input filenames and restart file existence will
                                    !! result in a new run that is not initialized from restart files.
  logical :: new_run_set = .false.  !< If true, new_run has been determined for this restart_CS.
  logical :: checksum_required      !< If true, require the restart checksums to match and error out otherwise.
                                    !! Users may want to avoid this comparison if for example the restarts are
                                    !! made from a run with a different mask_table than the current run,
                                    !! in which case the checksums will not match and cause crash.
  character(len=240) :: restartfile !< The name or name root for MOM restart files.

  !> An array of descriptions of the registered fields
  type(field_restart), pointer :: restart_field(:) => NULL()

  !>@{ Pointers to the fields that have been registered for restarts
  type(p0d), pointer :: var_ptr0d(:) => NULL()
  type(p1d), pointer :: var_ptr1d(:) => NULL()
  type(p2d), pointer :: var_ptr2d(:) => NULL()
  type(p3d), pointer :: var_ptr3d(:) => NULL()
  type(p4d), pointer :: var_ptr4d(:) => NULL()
  !!@}
  integer :: max_fields !< The maximum number of restart fields
  type(FmsNetcdfDomainFile_t) :: fileObj
  
end type MOM_restart_CS

!> Register fields for restarts
interface register_restart_field
  module procedure register_restart_field_ptr4d, register_restart_field_4d
  module procedure register_restart_field_ptr3d, register_restart_field_3d
  module procedure register_restart_field_ptr2d, register_restart_field_2d
  module procedure register_restart_field_ptr1d, register_restart_field_1d
  module procedure register_restart_field_ptr0d, register_restart_field_0d
end interface

!> Indicate whether a field has been read from a restart file
interface query_initialized
  module procedure query_initialized_name
  module procedure query_initialized_0d, query_initialized_0d_name
  module procedure query_initialized_1d, query_initialized_1d_name
  module procedure query_initialized_2d, query_initialized_2d_name
  module procedure query_initialized_3d, query_initialized_3d_name
  module procedure query_initialized_4d, query_initialized_4d_name
end interface

contains

!> Register a 3-d field for restarts, providing the metadata in a structure
subroutine register_restart_field_ptr3d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  CS%novars = CS%novars+1! remove .res. from the file name since fms read automatically appends it to
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr3d")

  CS%var_ptr3d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr3d

!> Register a 4-d field for restarts, providing the metadata in a structure
subroutine register_restart_field_ptr4d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:,:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr4d")

  CS%var_ptr4d(CS%novars)%p => f_ptr
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr4d

!> Register a 2-d field for restarts, providing the metadata in a structure
subroutine register_restart_field_ptr2d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr2d")

  CS%var_ptr2d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr2d

!> Register a 1-d field for restarts, providing the metadata in a structure
subroutine register_restart_field_ptr1d(f_ptr, var_desc, mandatory, CS)
  real, dimension(:), target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr1d")

  CS%var_ptr1d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr1d

!> Register a 0-d field for restarts, providing the metadata in a structure
subroutine register_restart_field_ptr0d(f_ptr, var_desc, mandatory, CS)
  real,               target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  type(vardesc),              intent(in) :: var_desc  !< A structure with metadata about this variable
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "register_restart_field: Module must be initialized before it is used.")

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                     ! once the total number of fields is known.

  CS%restart_field(CS%novars)%vars = var_desc
  CS%restart_field(CS%novars)%mand_var = mandatory
  CS%restart_field(CS%novars)%initialized = .false.
  call query_vardesc(CS%restart_field(CS%novars)%vars, &
                     name=CS%restart_field(CS%novars)%var_name, &
                     caller="register_restart_field_ptr0d")

  CS%var_ptr0d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()

end subroutine register_restart_field_ptr0d

! The following provide alternate interfaces to register restarts.

!> Register a 4-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_4d(f_ptr, name, mandatory, CS, G, GV, longname, units, &
                                     hor_grid, z_grid, t_grid)
  real, dimension(:,:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                      !! required if the new file uses any vertical grid axes.
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  ! local
  type(vardesc) :: vd
  type(MOM_restart_CS) :: fileObj
  logical :: file_open_success = .false.
  logical :: file_exists = .false.
  character(len=200) :: file_Name_1, file_Name_2, file_name_join
  integer :: str_split_index = 1, str_end_index = 1
  character(len=50) :: dimNames(4), dimNameStr
  character(len=8) :: ncAction 
  character(len=200) :: mesg
  logical :: use_lath = .false., use_lonh = .false., &
             use_latq = .false., use_lonq = .false., &
             use_layer = .false., use_int = .false., &
             use_time = .false., use_periodic = .false.
  integer :: horgrid_position = 1
  integer :: num_axes
  integer        :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  real, pointer, dimension(:) :: &
    gridLatT => NULL(), & ! The latitude or longitude of T or B points for
    gridLatB => NULL(), & ! the purpose of labeling the output axes.
    gridLonT => NULL(), gridLonB => NULL()
  integer :: axis_length, name_length
          
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: "//&
      "register_restart_field_4d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  WRITE(mpp_pe()+2000,*) "register_restart_field_2d: registering restart variable ", trim(name)
  call flush(mpp_pe()+2000)

  gridLatT => G%gridLatT
  gridLatB => G%gridLatB
  gridLonT => G%gridLonT
  gridLonB => G%gridLonB
  isg = G%isg
  ieg = G%ieg 
  jsg = G%jsg
  jeg = G%jeg
  IsgB = G%IsgB
  IegB = G%IegB
  JsgB = G%JsgB
  JegB = G%JegB

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  ! remove '.res.' from the file name if present since fms read automatically appends it to the file name
  file_Name_1 = ''
  file_Name_2 = ''
  file_name_join = ''
  name_length = 0
  file_Name_1(1:len_trim(CS%restartfile)) = trim(CS%restartfile)
  str_split_index = INDEX(file_Name_1,'.res')
  str_end_index = INDEX(file_Name_1,'.nc') 
  if (str_split_index > 1) then 
     if (str_end_index > 1 ) then
        file_name_join = trim(file_Name_1(1:str_split_index-1)//&
                         file_Name_1(str_split_index+4:str_end_index-1))//'.nc'
     else
        file_name_join = trim(file_Name_1(1:str_split_index-1)//& 
                         file_Name_1(str_split_index+4:len_trim(file_Name_1)))//'.nc'   
     endif          
  else 
      file_name_join = trim(file_Name_1)
  endif

  name_length = len_trim(file_name_join)
  file_Name_2(1:name_length) = trim(file_name_join)
  ! check whether restart file exists
  ! note: file_exists function does not automatically insert .res., so use full file name
  ncAction = ' '
  file_exists = fms2_file_exists(file_Name_1)
  if (file_exists) then
     ncAction = "append"
  else
     ncAction = "write"      
  endif

  ! open the restart file for domain-decomposed write
  file_open_success=fms2_open_file(CS%fileObj, file_Name_2, ncAction, G%Domain%mpp_domain, is_restart = .true.)
  if (.not. file_open_success) then 
     write(mesg,'( "ERROR, unable to open restart file ",A) ') trim(file_Name_2)
     call MOM_error(FATAL,"MOM_restart:register_restart_field_4d: "//mesg)
  endif
                                          
  ! check if variable axes are registered in file, and register them if they are not
  ! 4d variables are lon x lat x vertical level x time
  axis_length = 1
  num_axes = 0
  dimNameStr = ''
 
  ! longitude
  call get_horizontal_grid_coordinates(vd%hor_grid,use_lath,use_lonh,use_latq,use_lonq,horgrid_position)
  if (use_lonh) then
     num_axes = num_axes+1
     axis_length = size(gridLonT(isg:ieg))
     call check_for_restart_axis(CS%fileObj,'lonh', axis_length)
     dimNameStr = 'lonh'
  elseif (use_lonq) then
     num_axes = num_axes+1
     axis_length = size(gridLonB(IsgB:IegB))
     call check_for_restart_axis(CS%fileObj,'lonq', axis_length)
     dimNameStr ='lonq'
  endif
  
  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  
  WRITE(mpp_pe()+2000,*) "register_restart_field_4d: dimNameStr, dimNames 1 ", trim(dimNameStr), dimNames(1)
  call flush(mpp_pe()+2000)
  ! latitude
  
  if (use_lath) then
     num_axes = num_axes+1
     axis_length = size(gridLatT(jsg:jeg))
     call check_for_restart_axis(CS%fileObj,'lath', axis_length)
     dimNameStr = 'lath'
  elseif (use_latq) then
     num_axes = num_axes+1
     axis_length = size(gridLatB(JsgB:JegB))
     call check_for_restart_axis(CS%fileObj,'latq', axis_length)
     dimNameStr = 'latq'
  endif

  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  
  WRITE(mpp_pe()+2000,*) "register_restart_field_4d: dimNameStr, dimNames 2 ", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)

  ! Vertical (z) grid 
  call get_vertical_grid_coordinates(vd%z_grid,use_layer,use_int)
  if ((use_layer .or. use_int) .and. .not.(present(GV))) then
      call MOM_error(FATAL, "MOM_restart:: register_restart_field_4d: "//&
                     " Vertical (z) grid is defined for the variable "// trim(name) //& 
                     " but GV is missing from arguments.")
  endif
     
  if (use_layer) then
     num_axes = num_axes+1
     axis_length = size(GV%sLayer(1:GV%ke))
     call check_for_restart_axis(CS%fileObj,'Layer', axis_length)
     dimNameStr = 'Layer'
  elseif (use_int) then
     num_axes = num_axes+1
     axis_length = size(GV%sInterface(1:GV%ke+1))
     call check_for_restart_axis(CS%fileObj,'Interface', axis_length)
     dimNameStr = 'Interface'
  endif

  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  
  WRITE(mpp_pe()+2000,*) "register_restart_field_4d: dimNameStr, dimNames 3 ", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)


  ! time (t) grid  
  call get_time_coordinates(vd%t_grid,use_time,use_periodic)
     
  if (use_time) then
     num_axes = num_axes+1
     call check_for_restart_axis(CS%fileObj,'Time',unlimited)
     dimNameStr = 'Time'
  elseif (use_periodic) then
     num_axes = num_axes+1
     call check_for_restart_axis(CS%fileObj,'Period',unlimited)
     dimNameStr =  'Period'
  endif

  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  
  WRITE(mpp_pe()+2000,*) "register_restart_field_4d: dimNameStr, dimNames 4 ", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)
 
  ! register the restart field
  call register_restart_field_ptr4d(f_ptr, vd, mandatory, CS)
  ! Need to get the dimension names to register the domain-decomposed variables
  ! The NUMBER of dimensions is defined in MOM_io::create_file
  ! 1. axis metadata are written to file, returns a an axis structure with the metadaa
  ! 2. axes for a given variable are determined
  ! 3. axes metadata are associated with variable by referencing the axis structure
    
  call fms2_register_restart_field(CS%fileObj, name, f_ptr, dimensions=dimNames(1:num_axes), domain_position=horgrid_position)
  ! register variable attributes
  if (present(units)) call fms2_register_variable_attribute(CS%fileObj,name,'units',vd%units)
  if (present(longname)) call fms2_register_variable_attribute(CS%fileObj,name,'long_name',vd%longname)

!  call fms2_register_variable_attribute(CS%fileObj, name, t_grid, vd%t_grid)
!  call fms2_register_variable_attribute(CS%fileObj, name, t_grid, vd%hor_grid)
  call fms2_close_file(CS%fileObj)

end subroutine register_restart_field_4d

!> Register a 3-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_3d(f_ptr, name, mandatory, CS, G, GV, longname, units, &
                                     hor_grid, z_grid, t_grid)
  real, dimension(:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                      !! required if the new file uses any vertical grid axes.
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  ! local
  type(vardesc) :: vd
  type(MOM_restart_CS) :: fileObj
  logical :: file_open_success = .false.
  logical :: file_exists = .false.
  character(len=200) :: file_Name_1, file_Name_2, file_name_join
  integer :: str_split_index = 1, str_end_index = 1
  character(len=50) :: dimNames(4), dimNameStr
  character(len=8) :: ncAction 
  character(len=200) :: mesg
  logical :: use_lath = .false., use_lonh = .false., &
             use_latq = .false., use_lonq = .false., &
             use_layer = .false., use_int = .false., &
             use_time = .false., use_periodic = .false.
  integer :: horgrid_position = 1
  integer :: num_axes
  integer        :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  real, pointer, dimension(:) :: &
    gridLatT => NULL(), & ! The latitude or longitude of T or B points for
    gridLatB => NULL(), & ! the purpose of labeling the output axes.
    gridLonT => NULL(), gridLonB => NULL()
  integer :: axis_length, name_length
          
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: "//&
      "register_restart_field_3d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  WRITE(mpp_pe()+2000,*) "register_restart_field_3d: registering restart variable ", trim(name)
  call flush(mpp_pe()+2000)

  gridLatT => G%gridLatT
  gridLatB => G%gridLatB
  gridLonT => G%gridLonT
  gridLonB => G%gridLonB
  isg = G%isg
  ieg = G%ieg 
  jsg = G%jsg
  jeg = G%jeg
  IsgB = G%IsgB
  IegB = G%IegB
  JsgB = G%JsgB
  JegB = G%JegB

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=z_grid, t_grid=t_grid)

  ! remove '.res.' from the file name if present since fms read automatically appends it to the file name
  file_Name_1 = ''
  file_Name_2 = ''
  file_name_join = ''
  name_length = 0
  file_Name_1(1:len_trim(CS%restartfile)) = trim(CS%restartfile)
  str_split_index = INDEX(file_Name_1,'.res')
  str_end_index = INDEX(file_Name_1,'.nc') 
  if (str_split_index > 1) then 
     if (str_end_index > 1 ) then
        file_name_join = trim(file_Name_1(1:str_split_index-1)//&
                         file_Name_1(str_split_index+4:str_end_index-1))//'.nc'
     else
        file_name_join = trim(file_Name_1(1:str_split_index-1)//& 
                         file_Name_1(str_split_index+4:len_trim(file_Name_1)))//'.nc'   
     endif          
  else 
      file_name_join = trim(file_Name_1)
  endif

  name_length = len_trim(file_name_join)
  file_Name_2(1:name_length) = trim(file_name_join)
  ! Debug problems with string parsing
  !WRITE(mpp_pe()+2000,*) "register_restart_field_3d"
  !call flush(mpp_pe()+2000)
  !WRITE(mpp_pe()+2000,*) "string split index is: ",str_split_index
  !call flush(mpp_pe()+2000)
  !WRITE(mpp_pe()+2000,*) "file_Name_1: ",file_Name_1
  !call flush(mpp_pe()+2000)
  !WRITE(mpp_pe()+2000,*) "file_name_join: ",file_name_join
  !call flush(mpp_pe()+2000)
  !WRITE(mpp_pe()+2000,*) "file_Name_2: ",file_Name_2
  !call flush(mpp_pe()+2000)
  ! check whether restart file exists
  ! note: file_exists function does not automatically insert .res., so use full file name
  ncAction = ' '
  file_exists = fms2_file_exists(file_Name_1)
  if (file_exists) then
     ncAction = "append"
  else
     ncAction = "write"      
  endif

  ! open the restart file for domain-decomposed write
  file_open_success=fms2_open_file(CS%fileObj, file_Name_2, ncAction, G%Domain%mpp_domain, is_restart = .true.)
  if (.not. file_open_success) then 
     write(mesg,'( "ERROR, unable to open restart file ",A) ') trim(file_Name_2)
     call MOM_error(FATAL,"MOM_restart:register_restart_field_3d: "//mesg)
  endif
 
  ! check if variable axes are registered in file, and register them if they are not
  axis_length = 0
  num_axes = 0
  dimNameStr = ''                                      
  ! longitude
  call get_horizontal_grid_coordinates(vd%hor_grid,use_lath,use_lonh,use_latq,use_lonq,horgrid_position)

  if (use_lonh) then
     num_axes = num_axes+1
     axis_length = size(gridLonT(isg:ieg))
     call check_for_restart_axis(CS%fileObj,'lonh', axis_length)
     dimNameStr = 'lonh'
  elseif (use_lonq) then
     num_axes = num_axes+1
     axis_length = size(gridLonB(IsgB:IegB))
     call check_for_restart_axis(CS%fileObj,'lonq', axis_length)
     dimNameStr ='lonq'
  endif
  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  
  WRITE(mpp_pe()+2000,*) "register_restart_field_3d: dimNameStr, dimNames 1 ", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)
  ! latitude
  if (use_lath) then
     num_axes = num_axes+1
     axis_length = size(gridLatT(jsg:jeg))
     call check_for_restart_axis(CS%fileObj,'lath', axis_length)
     dimNameStr = 'lath'
  elseif (use_latq) then
     num_axes = num_axes+1
     axis_length = size(gridLatB(JsgB:JegB))
     call check_for_restart_axis(CS%fileObj,'latq', axis_length)
     dimNameStr = 'latq'
  endif
 
  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  ! Check that dim Name 2 is defined
  WRITE(mpp_pe()+2000,*) "register_restart_field_3d: dimNameStr, dimNames 2 ", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)

  ! Vertical (z) grid 
  call get_vertical_grid_coordinates(vd%z_grid,use_layer,use_int)
  if ((use_layer .or. use_int) .and. .not.(present(GV))) then
      call MOM_error(FATAL, "MOM_restart:: register_restart_field_3d: "//&
                     " Vertical (z) grid is defined for the variable "// trim(name) //& 
                     " but GV is missing from arguments.")
  endif
        
  if (use_layer) then
     num_axes = num_axes+1
     axis_length = size(GV%sLayer(1:GV%ke))
     call check_for_restart_axis(CS%fileObj,'Layer', axis_length)
     dimNameStr = 'Layer'
  elseif (use_int) then
     num_axes = num_axes+1
     axis_length = size(GV%sInterface(1:GV%ke+1))
     call check_for_restart_axis(CS%fileObj,'Interface', axis_length)
     dimNameStr = 'Interface'
  endif
 
  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)

  ! time (t) grid  
  call get_time_coordinates(vd%t_grid, use_time, use_periodic)
     
  if (use_time) then
     num_axes = num_axes+1
     call check_for_restart_axis(CS%fileObj,'Time',unlimited)
     dimNameStr = 'Time'
  elseif (use_periodic) then
     num_axes = num_axes+1
     call check_for_restart_axis(CS%fileObj,'Period',unlimited)
     dimNameStr = 'Period'
  endif

  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  ! Check that dim Name 3 is defined
  WRITE(mpp_pe()+2000,*) "register_restart_field_3d: dimNameStr, dimNames 3 ", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)

  ! register the restart field
  call register_restart_field_ptr3d(f_ptr, vd, mandatory, CS) 
  call fms2_register_restart_field(CS%fileObj, name, f_ptr, dimensions=dimNames(1:num_axes), domain_position=horgrid_position)
  ! register variable attributes
  if (present(units)) call fms2_register_variable_attribute(CS%fileObj,name,'units',vd%units)
  if (present(longname)) call fms2_register_variable_attribute(CS%fileObj,name,'long_name',vd%longname)

  call fms2_close_file(CS%fileObj)

end subroutine register_restart_field_3d

!> Register a 2-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_2d(f_ptr, name, mandatory, CS, G, GV, longname, units, &
                                     hor_grid, z_grid, t_grid)
  real, dimension(:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure which is
                                                      !! required if the new file uses any vertical grid axes.
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, 'h' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, '1' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  ! local
  type(vardesc) :: vd
  character(len=8) :: Zgrid
  type(MOM_restart_CS) :: fileObj
  logical :: file_open_success = .false.
  logical :: file_exists = .false.
  character(len=200) :: file_Name_1, file_Name_2, file_name_join
  integer :: str_split_index = 1, str_end_index = 1
  character(len=50) :: dimNames(4), dimNameStr
  character(len=8) ::  ncAction
  character(len=200) :: mesg
  logical :: use_lath = .false., use_lonh = .false., &
             use_latq = .false., use_lonq = .false., &
             use_layer = .false., use_int = .false., &
             use_time = .false., use_periodic = .false.
  integer :: horgrid_position = 1
  integer :: num_axes=0
  integer        :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  real, pointer, dimension(:) :: &
    gridLatT => NULL(), & ! The latitude or longitude of T or B points for
    gridLatB => NULL(), & ! the purpose of labeling the output axes.
    gridLonT => NULL(), gridLonB => NULL()
  integer :: axis_length, name_length
          
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: "//&
      "register_restart_field_2d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  WRITE(mpp_pe()+2000,*) "register_restart_field_2d: registering restart variable ", trim(Name)
  call flush(mpp_pe()+2000)

  gridLatT => G%gridLatT
  gridLatB => G%gridLatB
  gridLonT => G%gridLonT
  gridLonB => G%gridLonB
  isg = G%isg
  ieg = G%ieg 
  jsg = G%jsg
  jeg = G%jeg
  IsgB = G%IsgB
  IegB = G%IegB
  JsgB = G%JsgB
  JegB = G%JegB
  
  if (present(z_grid)) then 
     Zgrid = z_grid
  else
     Zgrid = '1' ; 
  endif

  vd = var_desc(name, units=units, longname=longname, hor_grid=hor_grid, &
                z_grid=Zgrid, t_grid=t_grid)

  ! remove '.res.' from the file name if present since fms read automatically appends it to the file name
  file_Name_1 = ''
  file_Name_2 = ''
  file_name_join = ''
  name_length = 0
  file_Name_1(1:len_trim(CS%restartfile)) = trim(CS%restartfile)
  str_split_index = INDEX(file_Name_1,'.res')
  str_end_index = INDEX(file_Name_1,'.nc') 
  if (str_split_index > 1) then 
     if (str_end_index > 1 ) then
        file_name_join = trim(file_Name_1(1:str_split_index-1)//&
                         file_Name_1(str_split_index+4:str_end_index-1))//'.nc'
     else
        file_name_join = trim(file_Name_1(1:str_split_index-1)//& 
                         file_Name_1(str_split_index+4:len_trim(file_Name_1)))//'.nc'   
     endif          
  else 
      file_name_join = trim(file_Name_1)
  endif

  name_length = len_trim(file_name_join)
  file_Name_2(1:name_length) = trim(file_name_join)
  ! check whether restart file exists
  ! note: file_exists function does not automatically insert .res., so use full file name
  ncAction = ' '
  file_exists = fms2_file_exists(file_Name_1)
  if (file_exists) then
     ncAction = "append"
  else
     ncAction = "write"      
  endif

  ! open the restart file for domain-decomposed write
  file_open_success=fms2_open_file(CS%fileObj, trim(file_Name_2), ncAction, G%Domain%mpp_domain, is_restart = .true.)
  if (.not. file_open_success) then 
     write(mesg,'( "ERROR, unable to open restart file ",A) ') trim(file_Name_2)
     call MOM_error(FATAL,"MOM_restart:register_restart_field_2d: "//mesg)
  endif
                                          
  ! check if global axis variables are registered in file, and register them if they are not
  axis_length = 1
  num_axes = 0
  dimNameStr = ''
  
  ! longitude
  call get_horizontal_grid_coordinates(vd%hor_grid,use_lath,use_lonh,use_latq,use_lonq,horgrid_position)
  if (use_lonh) then
     num_axes = num_axes+1
     axis_length = size(gridLonT(isg:ieg))
     call check_for_restart_axis(CS%fileObj,'lonh', axis_length)
     dimNameStr = 'lonh'
  elseif (use_lonq) then
     num_axes = num_axes+1
     axis_length = size(gridLonB(IsgB:IegB))
     call check_for_restart_axis(CS%fileObj,'lonq', axis_length)
     dimNameStr ='lonq'
  endif

  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  
  WRITE(mpp_pe()+2000,*) "register_restart_field_2d: dimNameStr, dimNames 1 ", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)
  ! latitude
  if (use_lath) then
     num_axes = num_axes+1
     axis_length = size(gridLatT(jsg:jeg))
     call check_for_restart_axis(CS%fileObj,'lath', axis_length)
     dimNameStr = 'lath'
  elseif (use_latq) then
     num_axes = num_axes+1
     axis_length = size(gridLatB(JsgB:JegB))
     call check_for_restart_axis(CS%fileObj,'latq', axis_length)
     dimNameStr = 'latq'
  endif
  dimNames(num_axes) = ''        
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  
  WRITE(mpp_pe()+2000,*) "register_restart_field_2d: dimNameStr, dimNames 2 ", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)

  ! Vertical (z) grid 
  call get_vertical_grid_coordinates(vd%z_grid,use_layer,use_int)

  !if (use_layer .or. use_int .and. .not.(present(GV))) then
  !   call MOM_error(FATAL, "MOM_restart:: register_restart_field_2d: "//&
  !                   " Vertical (z) grid is defined, but GV is missing from arguments.")
  !endif

  if (use_layer) then
     num_axes = num_axes+1
     axis_length = size(GV%sLayer(1:GV%ke))
     call check_for_restart_axis(CS%fileObj,'Layer', axis_length)
     dimNameStr = 'Layer'
  elseif (use_int) then
     num_axes = num_axes+1
     axis_length = size(GV%sInterface(1:GV%ke+1))
     call check_for_restart_axis(CS%fileObj,'Interface', axis_length)
     dimNameStr = 'Interface'
  endif

  ! time (t) grid  
  call get_time_coordinates(vd%t_grid,use_time,use_periodic)
     
  if (use_time) then
     num_axes = num_axes+1
     call check_for_restart_axis(CS%fileObj,'Time',unlimited)
     dimNameStr = 'Time'
  elseif (use_periodic) then
     num_axes = num_axes+1
     call check_for_restart_axis(CS%fileObj,'Period',unlimited)
     dimNameStr = 'Period'
  endif

  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  ! Check that dim Name 3 is defined
  WRITE(mpp_pe()+2000,*) "register_restart_field_2d: dimNameStr, dimNames 3 ", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)
 
  ! register the restart field
  call register_restart_field_ptr2d(f_ptr, vd, mandatory, CS)
  call fms2_register_restart_field(CS%fileObj, name, f_ptr, dimensions=dimNames(1:num_axes), domain_position=horgrid_position)
  ! register variable attributes
  if (present(units)) call fms2_register_variable_attribute(CS%fileObj,name,'units',vd%units)
  if (present(longname)) call fms2_register_variable_attribute(CS%fileObj,name,'long_name',vd%longname)

  call fms2_close_file(CS%fileObj)

end subroutine register_restart_field_2d

!> Register a 1-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_1d(f_ptr, name, mandatory, CS, G, GV, longname, units, &
                                     hor_grid, z_grid, t_grid)
  real, dimension(:), target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure which is
                                                      !! required if the new file uses any vertical grid axes.
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: hor_grid  !< variable horizonal staggering, '1' if absent
  character(len=*), optional, intent(in) :: z_grid    !< variable vertical staggering, 'L' if absent
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  ! local 
  type(vardesc) :: vd
  character(len=8) :: hgrid
  type(MOM_restart_CS) :: fileObj
  logical :: file_open_success = .false.
  logical :: file_exists = .false.
  character(len=50) :: dimNames(4), dimNameStr
  character(len=200) :: file_Name_1, file_Name_2, file_name_join
  character(len=8) :: ncAction
  integer :: str_split_index = 1, str_end_index = 1
  character(len=200) :: mesg
  logical :: use_layer = .false., use_int = .false., &
             use_time = .false., use_periodic = .false.
  integer :: horgrid_position = 1
  integer :: num_axes
  integer :: axis_length, name_length

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_3d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  WRITE(mpp_pe()+2000,*) "register_restart_field_1d: registering restart variable ", trim(Name)
  call flush(mpp_pe()+2000)

  if (present(hor_grid)) then 
     hgrid = hor_grid
  else 
     hgrid = '1'
  endif 

  vd = var_desc(name, units=units, longname=longname, hor_grid=hgrid, &
                z_grid=z_grid, t_grid=t_grid)

  ! remove '.res.' from the file name if present since fms read automatically appends it to the file name
  file_Name_1 = ''
  file_Name_2 = ''
  file_name_join = ''
  name_length = 0
  file_Name_1(1:len_trim(CS%restartfile)) = trim(CS%restartfile)
  str_split_index = INDEX(file_Name_1,'.res')
  str_end_index = INDEX(file_Name_1,'.nc') 
  if (str_split_index > 1) then 
     if (str_end_index > 1 ) then
        file_name_join = trim(file_Name_1(1:str_split_index-1)//&
                         file_Name_1(str_split_index+4:str_end_index-1))//'.nc'
     else
        file_name_join = trim(file_Name_1(1:str_split_index-1)//& 
                         file_Name_1(str_split_index+4:len_trim(file_Name_1)))//'.nc'   
     endif          
  else 
      file_name_join = trim(file_Name_1)
  endif

  name_length = len_trim(file_name_join)
  file_Name_2(1:name_length) = trim(file_name_join)
  ! check whether restart file exists
  ! note: file_exists function does not automatically insert .res., so use full file name
  ncAction = ' '
  file_exists = fms2_file_exists(file_Name_1)
  if (file_exists) then
     ncAction = "append"
  else
     ncAction = "write"      
  endif

  ! open the restart file for domain-decomposed write
  file_open_success=fms2_open_file(CS%fileObj, file_Name_2, ncAction, G%Domain%mpp_domain, is_restart = .true.)
  if (.not. file_open_success) then 
     write(mesg,'( "ERROR, unable to open restart file ",A) ') trim(file_Name_2)
     call MOM_error(FATAL,"MOM_restart:register_restart_field_4d: "//mesg)
  endif
                                          
  ! Vertical (z) grid
  axis_length = 1
  num_axes = 0
  dimNameStr = '' 
  call get_vertical_grid_coordinates(vd%z_grid,use_layer,use_int)

  !if ((use_layer .or. use_int) .and. .not.(present(GV))) then
  !    call MOM_error(FATAL, "MOM_restart:: register_restart_field_1d: "//&
  !                   " Vertical (z) grid is defined, but GV is missing from arguments.")
  !endif
     
  if (use_layer) then
     num_axes = num_axes+1
     axis_length = size(GV%sLayer(1:GV%ke))
     call check_for_restart_axis(CS%fileObj,'Layer', axis_length)
     dimNameStr = 'Layer'
  elseif (use_int) then
     num_axes = num_axes+1
     axis_length = size(GV%sInterface(1:GV%ke+1))
     call check_for_restart_axis(CS%fileObj,'Interface', axis_length)
     dimNameStr = 'Interface'
  endif

  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)

  ! time (t) grid  
  call get_time_coordinates(vd%t_grid,use_time,use_periodic)
     
  if (use_time) then
     num_axes = num_axes+1
     call check_for_restart_axis(CS%fileObj,'Time',unlimited)
     dimNameStr = 'Time'
  elseif (use_periodic) then
     num_axes = num_axes+1
     call check_for_restart_axis(CS%fileObj,'Period',unlimited)
     dimNameStr = 'Period'
  endif

  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)
  ! Check that dim Name is defined
  WRITE(mpp_pe()+2000,*) "register_restart_field_1d: dimNameStr, dimNames", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)

  ! register the restart field
  call register_restart_field_ptr1d(f_ptr, vd, mandatory, CS)
  call fms2_register_restart_field(CS%fileObj, name, f_ptr, dimensions=dimNames(1:num_axes), domain_position=horgrid_position)
  ! register variable attributes
  if (present(units)) call fms2_register_variable_attribute(CS%fileObj,name,'units',vd%units)
  if (present(longname)) call fms2_register_variable_attribute(CS%fileObj,name,'long_name',vd%longname)

  call fms2_close_file(CS%fileObj)

end subroutine register_restart_field_1d

!> Register a 0-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_0d(f_ptr, name, mandatory, CS, G, GV, longname, units, &
                                     t_grid)
  real,               target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*),           intent(in) :: name      !< variable name to be used in the restart file
  logical,                    intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.
  type(MOM_restart_CS),       pointer    :: CS        !< A pointer to a MOM_restart_CS object (intent in/out)
  type(ocean_grid_type),      intent(in) :: G         !< The ocean's grid structure
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure which is
                                                      !! required if the new file uses any vertical grid axes.
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: t_grid    !< time description: s, p, or 1, 's' if absent
  ! local
  logical :: file_open_success = .false. ! returned by call to fms2_open_file 
  logical :: file_exists = .false.
  type(MOM_restart_CS) :: fileObj ! fms2 data structure
  type(vardesc) :: vd
  character(len=50) :: dimNames(1), dimNameStr
  character(len=8) :: ncAction
  character(len=200) :: file_Name_1, file_Name_2, file_name_join
  integer :: str_split_index = 1, str_end_index = 1
  character(len=200) :: mesg
  logical :: use_time = .false., use_periodic = .false.
  integer :: horgrid_position = 1
  integer :: num_axes
  integer :: axis_length, name_length
          
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart: " // &
      "register_restart_field_0d: Module must be initialized before "//&
      "it is used to register "//trim(name))

  WRITE(mpp_pe()+2000,*) "register_restart_field_0d: registering restart variable ", trim(Name)
  call flush(mpp_pe()+2000)
  
  vd = var_desc(name, units=units, longname=longname, hor_grid='1', &
                z_grid='1', t_grid=t_grid)

  ! remove '.res.' from the file name if present since fms read automatically appends it to the file name
  file_Name_1 = ''
  file_Name_2 = ''
  file_name_join = ''
  name_length = 0
  file_Name_1(1:len_trim(CS%restartfile)) = trim(CS%restartfile)
  str_split_index = INDEX(file_Name_1,'.res')
  str_end_index = INDEX(file_Name_1,'.nc') 
  if (str_split_index > 1) then 
     if (str_end_index > 1 ) then
        file_name_join = trim(file_Name_1(1:str_split_index-1)//&
                         file_Name_1(str_split_index+4:str_end_index-1))//'.nc'
     else
        file_name_join = trim(file_Name_1(1:str_split_index-1)//& 
                         file_Name_1(str_split_index+4:len_trim(file_Name_1)))//'.nc'   
     endif          
  else 
      file_name_join = trim(file_Name_1)
  endif

  name_length = len_trim(file_name_join)
  file_Name_2(1:name_length) = trim(file_name_join)
  ! check whether restart file exists
  ! note: file_exists function does not automatically insert .res., so use full file name
  ncAction = ' '
  file_exists = fms2_file_exists(file_Name_1)
  if (file_exists) then
     ncAction = "append"
  else
     ncAction = "write"      
  endif

  ! open the restart file for domain-decomposed write
  file_open_success=fms2_open_file(CS%fileObj, file_Name_2, ncAction, G%Domain%mpp_domain, is_restart = .true.)
  if (.not. file_open_success) then 
     write(mesg,'( "ERROR, unable to open restart file ",A) ') trim(file_Name_2)
     call MOM_error(FATAL,"MOM_restart:register_restart_field_0d: "//mesg)
  endif
                                          
  ! check if global axis variables are registered in file, and register them if they are not
  ! Vertical (z) grid
  axis_length = 1
  num_axes = 0
  dimNameStr = '' 
 
  ! time (t) grid  
  call get_time_coordinates(vd%t_grid, use_time, use_periodic)
     
  if (use_time) then
     num_axes = num_axes+1
     call check_for_restart_axis(CS%fileObj,'Time', axis_length)
     dimNameStr = 'Time'
  elseif (use_periodic) then
     num_axes = num_axes+1
     call check_for_restart_axis(CS%fileObj,'Period', axis_length)
     dimNameStr = 'Period'
  endif

  dimNames(num_axes) = ''
  name_length = len_trim(dimNameStr)
  dimNames(num_axes)(1:name_length) = trim(dimNameStr)

  ! Check that dim Name is defined
  WRITE(mpp_pe()+2000,*) "register_restart_field_0d: dimNameStr, dimNames", trim(dimNameStr), dimNames(num_axes)
  call flush(mpp_pe()+2000)
  ! register the restart field
  call register_restart_field_ptr0d(f_ptr, vd, mandatory, CS)
  call fms2_register_restart_field(CS%fileObj, name, f_ptr, dimensions=dimNames, domain_position=horgrid_position)
  ! register variable attributes
  if (present(units)) call fms2_register_variable_attribute(CS%fileObj,name,'units',vd%units)
  if (present(longname)) call fms2_register_variable_attribute(CS%fileObj,name,'long_name',vd%longname)

  call fms2_close_file(CS%fileObj)

end subroutine register_restart_field_0d

!> query_initialized_name determines whether a named field has been successfully
!! read from a restart file yet.
function query_initialized_name(name, CS) result(query_initialized)
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine returns .true. if the field referred to by name has
! initialized from a restart file, and .false. otherwise.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (trim(name) == CS%restart_field(m)%var_name) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if ((n==CS%novars+1) .and. (is_root_pe())) &
    call MOM_error(NOTE,"MOM_restart: Unknown restart variable "//name// &
                        " queried for initialization.")

  if ((is_root_pe()) .and. query_initialized) &
    call MOM_error(NOTE,"MOM_restart: "//name// &
                         " initialization confirmed by name.")

end function query_initialized_name

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_0d(f_ptr, CS) result(query_initialized)
  real,         target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr0d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_0d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_1d(f_ptr, CS) result(query_initialized)
  real, dimension(:), target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS),       pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr1d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_1d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_2d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:), &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr2d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_2d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_3d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:,:), &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr3d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_3d

!> Indicate whether the field pointed to by f_ptr has been initialized from a restart file.
function query_initialized_4d(f_ptr, CS) result(query_initialized)
  real, dimension(:,:,:,:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr has
! been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr4d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.

end function query_initialized_4d

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_0d_name(f_ptr, name, CS) result(query_initialized)
  real,         target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr0d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_0d_name

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_1d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr1d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_1d_name

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_2d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.

  integer :: m,n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr2d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE,"MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "probably because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_2d_name

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_3d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:,:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.

  integer :: m, n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr3d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE, "MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "possibly because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_3d_name

!> Indicate whether the field pointed to by f_ptr or with the specified variable
!! name has been initialized from a restart file.
function query_initialized_4d_name(f_ptr, name, CS) result(query_initialized)
  real, dimension(:,:,:,:),  &
                target, intent(in) :: f_ptr !< A pointer to the field that is being queried
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  type(MOM_restart_CS), pointer    :: CS !< A pointer to a MOM_restart_CS object (intent in)
  logical :: query_initialized
!   This subroutine tests whether the field pointed to by f_ptr or with the
! specified variable name has been initialized from a restart file.

  integer :: m, n
  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (associated(CS%var_ptr4d(m)%p,f_ptr)) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
! Assume that you are going to initialize it now, so set flag to initialized if
! queried again.
  if (n<=CS%novars) CS%restart_field(n)%initialized = .true.
  if (n==CS%novars+1) then
    if (is_root_pe()) &
      call MOM_error(NOTE, "MOM_restart: Unable to find "//name//" queried by pointer, "//&
        "possibly because of the suspect comparison of pointers by ASSOCIATED.")
    query_initialized = query_initialized_name(name, CS)
  endif

end function query_initialized_4d_name

!> save_restart saves all registered variables to restart files.
subroutine save_restart(directory, time, G, CS, time_stamped, filename, GV)
  character(len=*),        intent(in)    :: directory !< The directory where the restart files
                                                  !! are to be written
  type(time_type),         intent(in)    :: time  !< The current model time
  type(ocean_grid_type),   intent(inout) :: G     !< The ocean's grid structure
  type(MOM_restart_CS),    pointer       :: CS    !< The control structure returned by a previous
                                                  !! call to restart_init.
  logical,          optional, intent(in) :: time_stamped !< If present and true, add time-stamp
                                                  !! to the restart file names.
  character(len=*), optional, intent(in) :: filename !< A filename that overrides the name in CS%restartfile.
  type(verticalGrid_type), optional, intent(in) :: GV   !< The ocean's vertical grid structure

  ! Local variables
  type(vardesc) :: vars(CS%max_fields)  ! Descriptions of the fields that
                                        ! are to be read from the restart file.
  type(fieldtype) :: fields(CS%max_fields) !
  character(len=512) :: restartpath     ! The restart file path (dir/file).
  character(len=256) :: restartname     ! The restart file name (no dir).
  character(len=8)   :: suffix          ! A suffix (like _2) that is appended
                                        ! to the name of files after the first.
  integer(kind=8) :: var_sz, size_in_file ! The size in bytes of each variable
                                        ! and the variables already in a file.
  integer(kind=8) :: max_file_size = 2147483647_8 ! The maximum size in bytes
                                        ! for any one file.  With NetCDF3,
                                        ! this should be 2 Gb or less.
  integer :: start_var, next_var        ! The starting variables of the
                                        ! current and next files.
  integer :: unit                       ! The mpp unit of the open file.
  integer :: m, nz, num_files, var_periods
  integer :: seconds, days, year, month, hour, minute
  character(len=8) :: hor_grid, z_grid, t_grid ! Variable grid info.
  character(len=8) :: t_grid_read
  character(len=64) :: var_name         ! A variable's name.
  real :: restart_time
  character(len=32) :: filename_appendix = '' !fms appendix to filename for ensemble runs
  integer :: length
  integer(kind=8) :: check_val(CS%max_fields,1)
  integer :: isL, ieL, jsL, jeL, pos
  integer :: str_split_index, ext_index
  type(MOM_restart_CS) :: fileObj
  logical :: file_open_success
  logical :: use_latq, use_lonq, use_lath,use_lonh, &
             use_layer, use_int, use_time, use_period

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "save_restart: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  ! With parallel read & write, it is possible to disable the following...

  ! The maximum file size is 4294967292, according to the NetCDF documentation.
  if (CS%large_file_support) max_file_size = 4294967292_8

  num_files = 0
  next_var = 0
  nz = 1 ; if (present(GV)) nz = GV%ke

  restart_time = time_type_to_real(time) / 86400.0

  use_lath = .false.
  use_lonh = .false.
  use_latq = .false.
  use_lonq = .false.
  use_layer = .false.
  use_int = .false.
  use_time = .false.
  use_period = .false.
  pos = 1
  restartname = " "
  restartpath = " "
  str_split_index = 1
  ext_index = 1

  if (present(filename)) then 
     ! remove .res. from name because fms2 adds it
     str_split_index = INDEX(filename,'.res.')
     ext_index = INDEX(filename,'.nc')
     if ((str_split_index > 1) .and. (ext_index >1)) then 
        restartname = trim(filename(1:str_split_index-1))&
        //trim(filename(ext_index:ext_index+2))
     elseif ((str_split_index > 1) .and. (ext_index <= 1)) then
        restartname = trim(filename(1:str_split_index-1))
     else
        restartname = trim(filename)
     endif    
  else
     str_split_index = INDEX(filename,'.res.')
     ext_index = INDEX(filename,'.nc')
     if ((str_split_index > 1) .and. (ext_index >1)) then 
        restartname=trim(CS%restartfile(1:str_split_index-1))&
        //trim(filename(ext_index:ext_index+2))
     elseif ((str_split_index > 1) .and. (ext_index <= 1)) then
        restartname=trim(CS%restartfile(1:str_split_index-1))
     else
        restartname=trim(CS%restartfile)
     endif
  endif

  if (PRESENT(time_stamped)) then ; if (time_stamped) then
    call get_date(time,year,month,days,hour,minute,seconds)
    ! Compute the year-day, because I don't like months. - RWH
    do m=1,month-1
      days = days + days_in_month(set_date(year,m,2,0,0,0))
    enddo
    seconds = seconds + 60*minute + 3600*hour
    if (year <= 9999) then
      write(restartname,'("_Y",I4.4,"_D",I3.3,"_S",I5.5)') year, days, seconds
    elseif (year <= 99999) then
      write(restartname,'("_Y",I5.5,"_D",I3.3,"_S",I5.5)') year, days, seconds
    else
      write(restartname,'("_Y",I10.10,"_D",I3.3,"_S",I5.5)') year, days, seconds
    endif
    restartname = trim(CS%restartfile)//trim(restartname)
  endif ; endif
 
  next_var = 1
  do while (next_var <= CS%novars )
     start_var = next_var
     size_in_file = 8*(2*G%Domain%niglobal+2*G%Domain%njglobal+2*nz+1000)
     
     do m=start_var,CS%novars
        call query_vardesc(CS%restart_field(m)%vars, hor_grid=hor_grid, &
                           z_grid=z_grid, t_grid=t_grid, caller="save_restart")
        if (hor_grid == '1') then
           var_sz = 8
        else
           var_sz = 8*(G%Domain%niglobal+1)*(G%Domain%njglobal+1)
        endif
        select case (z_grid)
           case ('L') ; var_sz = var_sz * nz
           case ('i') ; var_sz = var_sz * (nz+1)
        end select
        t_grid = adjustl(t_grid)
        if (t_grid(1:1) == 'p') then
           if (len_trim(t_grid(2:8)) > 0) then
              var_periods = -1
              t_grid_read = adjustl(t_grid(2:8))
              read(t_grid_read,*) var_periods
              if (var_periods > 1) var_sz = var_sz * var_periods
           endif
        endif

        if ((m==start_var) .OR. (size_in_file < max_file_size-var_sz)) then
           size_in_file = size_in_file + var_sz
        else ; exit
        endif

     enddo
     next_var = m

     !query fms_io if there is a filename_appendix (for ensemble runs)
     call get_filename_appendix(filename_appendix)
     if (len_trim(filename_appendix) > 0) then
        length = len_trim(restartname)
        if (restartname(length-2:length) == '.nc') then
           restartname = restartname(1:length-3)//'.'//trim(filename_appendix)//'.nc'
        else
           restartname = restartname(1:length)  //'.'//trim(filename_appendix)
        endif
     endif

     restartpath = trim(directory)//trim(restartname)

     if (num_files < 10) then
        write(suffix,'("_",I1)') num_files
     else
        write(suffix,'("_",I2)') num_files
     endif

     if (num_files > 0) restartpath = trim(restartpath) // trim(suffix)

     do m=start_var,next_var-1
        vars(m-start_var+1) = CS%restart_field(m)%vars
        call query_vardesc(vars(m-start_var+1), t_grid=t_grid, hor_grid=hor_grid,z_grid=z_grid, &
                           caller="save_restart")
        call get_horizontal_grid_coordinates(hor_grid,use_lath,use_lonh,use_latq,use_lonq,pos)
        call get_vertical_grid_coordinates(z_grid,use_layer,use_int)
        call get_time_coordinates(t_grid,use_time,use_period)
     enddo
     !call query_vardesc(vars(1), t_grid=t_grid, hor_grid=hor_grid, caller="save_restart")
     !t_grid = adjustl(t_grid)
     !if (t_grid(1:1) /= 'p') &
     !   call modify_vardesc(vars(1), t_grid='s', caller="save_restart")
     !select case (hor_grid)
     !   case ('q') ; pos = CORNER
     !   case ('h') ; pos = CENTER
     !   case ('u') ; pos = EAST_FACE
     !   case ('v') ; pos = NORTH_FACE
     !   case ('Bu') ; pos = CORNER
     !   case ('T')  ; pos = CENTER
      !  case ('Cu') ; pos = EAST_FACE
     !   case ('Cv') ; pos = NORTH_FACE
     !   case ('1') ; pos = 0
      !  case default ; pos = 0
     !end select
     ! write the axis (dimension) data to the restart file
     if (use_lath) then
        call write_axis_data(CS%fileObj,'lath',G=G)
     endif
     if (use_lonh) then
        call write_axis_data(CS%fileObj,'lonh',G=G)
     endif
     if (use_latq) then
        call write_axis_data(CS%fileObj,'latq',G=G)
     endif
     if (use_lonq) then
        call write_axis_data(CS%fileObj,'lonq',G=G)
     endif
     if (use_layer) then
        call write_axis_data(CS%fileObj,'Layer',GV=GV)
     endif
     if (use_int) then
        call write_axis_data(CS%fileObj,'Interface',GV=GV)
     endif
     if (use_time) then
        call write_axis_data(CS%fileObj,'Time',restart_time_in_days=restart_time)
     endif
     if (use_period) then
        call write_axis_data(CS%fileObj,'Period',t_grid_in=t_grid)
     endif
     !Prepare the checksum of the restart fields to be written to restart files
     call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL)
     do m=start_var,next_var-1
        if (associated(CS%var_ptr3d(m)%p)) then
           check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:))
        elseif (associated(CS%var_ptr2d(m)%p)) then
           check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL))
        elseif (associated(CS%var_ptr4d(m)%p)) then
           check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:))
        elseif (associated(CS%var_ptr1d(m)%p)) then
           check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr1d(m)%p)
        elseif (associated(CS%var_ptr0d(m)%p)) then
           check_val(m-start_var+1,1) = mpp_chksum(CS%var_ptr0d(m)%p,pelist=(/mpp_pe()/))
        endif
     enddo

     !if (CS%parallel_restartfiles) then
     !   call create_file(unit, trim(restartpath), vars, (next_var-start_var), &
     !                    fields, MULTIPLE, G=G, GV=GV, checksums=check_val)
     !else
     !   call create_file(unit, trim(restartpath), vars, (next_var-start_var), &
     !                    fields, SINGLE_FILE, G=G, GV=GV, checksums=check_val)
     !endif
     file_open_success = .false.
     file_open_success = fms2_open_file(CS%fileObj, restartpath,"write", &
                                        G%Domain%mpp_domain,is_restart=.true.)
     if (.not. (file_open_success)) then
        call MOM_error(FATAL,"MOM_restart::save_restart: Failed to open file "//restartpath)
     endif
     
     do m=start_var,next_var-1     
        if (associated(CS%var_ptr3d(m)%p)) then
           !call write_field(unit,fields(m-start_var+1), G%Domain%mpp_domain, &
           !                 CS%var_ptr3d(m)%p, restart_time)
           call fms2_write_data(CS%fileObj,vars(m)%name, CS%var_ptr3d(m)%p)
        elseif (associated(CS%var_ptr2d(m)%p)) then
           !call write_field(unit,fields(m-start_var+1), G%Domain%mpp_domain, &
           !                 CS%var_ptr2d(m)%p, restart_time)
           call fms2_write_data(CS%fileObj,vars(m)%name, CS%var_ptr2d(m)%p)
        elseif (associated(CS%var_ptr4d(m)%p)) then
           !call write_field(unit,fields(m-start_var+1), G%Domain%mpp_domain, &
           !                 CS%var_ptr4d(m)%p, restart_time)
           call fms2_write_data(CS%fileObj,vars(m)%name, CS%var_ptr4d(m)%p)
        elseif (associated(CS%var_ptr1d(m)%p)) then
           !call write_field(unit, fields(m-start_var+1), CS%var_ptr1d(m)%p, &
           !                 restart_time)
           call fms2_write_data(CS%fileObj,vars(m)%name, CS%var_ptr1d(m)%p)
        elseif (associated(CS%var_ptr0d(m)%p)) then
           !call write_field(unit, fields(m-start_var+1), CS%var_ptr0d(m)%p, &
           !                 restart_time)
           call fms2_write_data(CS%fileObj,vars(m)%name, CS%var_ptr0d(m)%p)
        endif
        ! write the checksum
        call fms2_register_variable_attribute(CS%fileObj,vars(m)%name,'checksum',check_val(m,1))      
     enddo

     !call close_file(unit)
     call fms2_close_file(CS%fileObj)
     num_files = num_files+1
  enddo

end subroutine save_restart

!> restore_state reads the model state from previously generated files.  All
!! restart variables are read from the first file in the input filename list
!! in which they are found.
subroutine restore_state(filename, directory, day, G, CS)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files.
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(time_type),       intent(out) :: day       !< The time of the restarted run
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init.

!  This subroutine reads the model state from previously
!  generated files.  All restart variables are read from the first
!  file in the input filename list in which they are found.

  ! Local variables
  character(len=200) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=80) :: fname      ! The name of the current file.
  character(len=8)  :: suffix     ! A suffix (like "_2") that is added to any
                                  ! additional restart files.
  character(len=512) :: mesg      ! A message for warnings.
  character(len=80) :: varname    ! A variable's name.
  integer :: num_file        ! The number of files (restart files and others
                             ! explicitly in filename) that are open.
  integer :: i, n, m, missing_fields, ic
  integer :: isL, ieL, jsL, jeL, is0, js0
  integer :: sizes(7)
  integer :: ndim, nvar, natt, ntime, pos

  integer :: unit(CS%max_fields) ! The mpp unit of all open files.
  character(len=200) :: unit_path(CS%max_fields) ! The file names.
  logical :: unit_is_global(CS%max_fields) ! True if the file is global.

  character(len=8)   :: hor_grid ! Variable grid info.
  real    :: t1, t2 ! Two times.
  real, allocatable :: time_vals(:)
  type(fieldtype), allocatable :: fields(:)
  logical                          :: check_exist, is_there_a_checksum
  integer(LONG_KIND),dimension(3)  :: checksum_file
  integer(kind=8)                  :: checksum_data
  logical :: file_open_success = .false. ! returned by call to fms2_open_file 
  type(MOM_restart_CS) :: fileObj   ! fms2 data structure
  integer :: str_split_index = 1
  integer :: str_end_index = 1
  character(len=200) :: file_path_1,file_path_2
  character(len=80), dimension(:), allocatable :: variable_names(:) ! File variable names
  character(len=64) :: checksum_char
  integer(LONG_KIND) :: checksumh
  integer :: num_checksumh, last, is, k
  logical :: var_exists = .false.

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "restore_state: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)
    
! Get number of restart files and full paths to files 
  if ((LEN_TRIM(filename) == 1) .and. (filename(1:1) == 'F')) then
     num_file = open_restart_units('r', directory, G, CS, &
                     file_paths=unit_path, global_files=unit_is_global)
  else
    num_file = open_restart_units(filename, directory, G, CS, &
                     file_paths=unit_path, global_files=unit_is_global)
  endif

  if (num_file == 0) then
     write(mesg,'("Unable to find any restart files specified by  ",A,"  in directory ",A,".")') &
                  trim(filename), trim(directory)
     call MOM_error(FATAL,"MOM_restart: "//mesg)
  endif

! Get the time from the first file in the list that has one.
  do n=1,num_file
     file_path_1(1:len_trim(unit_path(n))) = trim(unit_path(n))
     str_split_index = INDEX(file_path_1,'.res')
     str_end_index = INDEX(file_path_1,'.nc') 
     if (str_split_index > 1 .and. str_end_index > 1 )then ! remove .res. from the file name since fms read automatically appends it to the file name
         file_path_2 = trim(file_path_1(1:str_split_index-1)// &
                         file_path_1(str_split_index+4:str_end_index-1))//'.nc'                         
     else 
         file_path_2 = trim(file_path_1)
     endif
     file_open_success=fms2_open_file(CS%fileObj, trim(file_path_2),"read",  G%Domain%mpp_domain, is_restart = .true.)
                                                
     if(.not. file_open_success) then 
        write(mesg,'( "ERROR, unable to open restart file  ",A) ') trim(unit_path(n))
        call MOM_error(FATAL,"MOM_restart: "//mesg)
     endif
     call fms2_get_dimension_size(CS%fileObj, "Time", ntime)

     if (ntime < 1) cycle

     allocate(time_vals(ntime))
     call fms2_register_restart_field(CS%fileObj, "Time","double")
     call fms2_read_data(CS%fileObj,"Time", time_vals)
     t1 = time_vals(1)
     deallocate(time_vals)
     day = real_to_time(t1*86400.0)

     call fms2_close_file(CS%fileObj)
     exit
  enddo

  if (n>num_file) call MOM_error(WARNING,"MOM_restart: " // &
                                 "No times found in restart files.")

! Check the remaining files, if any, for different times and issue a warning
! if they differ from the first time.
  if (is_root_pe()) then
     do m = n+1,num_file
        file_path_1(1:len_trim(unit_path(n))) = trim(unit_path(n))
        str_split_index = INDEX(file_path_1,'.res')
        str_end_index = INDEX(file_path_1,'.nc')
        if (str_split_index > 1 .and. str_end_index > 1 )then
           file_path_2 = trim(file_path_1(1:str_split_index-1)// &
                         file_path_1(str_split_index+4:str_end_index-1))//'.nc'                         
        else 
          file_path_2 = trim(file_path_1)
        endif
        file_open_success=fms2_open_file(CS%fileObj, trim(file_path_2),"read", G%Domain%mpp_domain, is_restart = .true.)
                                                
        if(.not. file_open_success) then 
           write(mesg,'( "ERROR, unable to open restart file  ",A) ') trim(unit_path(n))
           call MOM_error(FATAL,"MOM_restart: "//mesg)
        endif
        
        call fms2_get_dimension_size(CS%fileObj, "Time", ntime)

        if (ntime < 1) cycle

        allocate(time_vals(ntime))
        call fms2_register_restart_field(CS%fileObj, "Time","double")
        call fms2_read_data(CS%fileObj,"Time", time_vals)
        t2 = time_vals(1)
        deallocate(time_vals)

        if (t1 /= t2) then
           write(mesg,'("WARNING: Restart file ",I2," has time ",F10.4,"whereas simulation is restarted at ",F10.4," (differing by ",F10.4,").")')&
              m,t1,t2,t1-t2
           call MOM_error(WARNING, "MOM_restart: "//mesg)
        endif
        call fms2_close_file(CS%fileObj)
     enddo
  endif

  if (n>num_file) call MOM_error(WARNING,"MOM_restart: " // &
                            "No times found in restart files.")

! Read each variable from the first file in which it is found.
  do n=1,num_file
     !call get_file_info(unit(n), ndim, nvar, natt, ntime)
     file_path_1(1:len_trim(unit_path(n))) = trim(unit_path(n))
     str_split_index = INDEX(file_path_1,'.res')
     str_end_index = INDEX(file_path_1,'.nc')
     if (str_split_index > 1 .and. str_end_index > 1 )then
        file_path_2 = trim(file_path_1(1:str_split_index-1)// &
                         file_path_1(str_split_index+4:str_end_index-1))//'.nc'                         
     else 
        file_path_2 = trim(file_path_1)
     endif
     file_open_success=fms2_open_file(CS%fileObj, trim(file_path_2),"read", G%Domain%mpp_domain, is_restart = .true.)
                                                
     if (.not. file_open_success) then 
        write(mesg,'( "ERROR, unable to open restart file  ",A )') trim(unit_path(n))
        call MOM_error(FATAL,"MOM_restart: "//mesg)
     endif
     ! register the horizontal axes
     call fms2_register_axis(CS%fileObj,'latq','y')
     call fms2_register_axis(CS%fileObj,'lath','y')
     call fms2_register_axis(CS%fileObj,'lonq','x')
     call fms2_register_axis(CS%fileObj,'lonh','x')

     ! get number of variables in the file
     nvar=fms2_get_num_variables(CS%fileObj)
     ! get the names of the variables in the file
     allocate(variable_names(nvar))
     call fms2_get_variable_names(CS%fileObj, variable_names)
     ! allocate the fields to hold the variable data
     allocate(fields(nvar))
    
     missing_fields = 0

     do m=1,CS%novars
        if (CS%restart_field(m)%initialized) cycle
        call query_vardesc(CS%restart_field(m)%vars, hor_grid=hor_grid, &
                             caller="restore_state")
        select case (hor_grid)
           case ('q') ; pos = CORNER
           case ('h') ; pos = CENTER
           case ('u') ; pos = EAST_FACE
           case ('v') ; pos = NORTH_FACE
           case ('Bu') ; pos = CORNER
           case ('T')  ; pos = CENTER
           case ('Cu') ; pos = EAST_FACE
           case ('Cv') ; pos = NORTH_FACE
           case ('1') ; pos = 0
           case default ; pos = 0
        end select

        call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL)
        do i=1, nvar
           varname = " "
           checksum_char = " "
           var_exists = .false.
           varname(1:len_trim(CS%restart_field(m)%var_name)) = trim(CS%restart_field(m)%var_name)
           ! check if variable is in the restart file
           var_exists = fms2_variable_exists(CS%fileObj, varname)

           if (var_exists) then
              if (.NOT. CS%checksum_required) then
                 is_there_a_checksum = .false. ! Do not need to do data checksumming.
              else
                 check_exist = fms2_attribute_exists(CS%fileObj, varname, "checksum")
                 checksum_file(:) = -1
                 checksum_data = -1
                 is_there_a_checksum = .false.
                 if ( check_exist ) then
                    ! The following checksum conversion proceure is adapted from mpp_get_atts 
                    call fms2_get_variable_attribute(CS%fileObj,varname,"checksum",checksum_char)

                    last = len_trim(checksum_char)
                    is = index(trim(checksum_char),",") ! A value of 0 implies only 1 checksum value
                    ! Scan checksum character array for the ',' delimiter, which indicates that the corresponding variable
                    ! has multiple time levels.
                    checksumh = 0
                    num_checksumh = 1
                    do while ((is > 0) .and. (is < (last-15)))
                       is = is + scan(checksum_char(is:last), "," ) ! move starting pointer after ","
                       num_checksumh = num_checksumh + 1
                    enddo
           
                    is = 1
  
                    do k = 1, num_checksumh
                       read(checksum_char(is:is+15),'(Z16)') checksumh ! Z=hexadecimal integer: Z16 is for 64-bit data types
                       checksum_file(k) = checksumh 
                       is = is + 17 ! Move index past the ',' in checksum_char
                    enddo

                    is_there_a_checksum = .true.
                 endif
              endif
              ! register the restart variable
              call fms2_register_restart_field(CS%fileObj,varname,'real')

              if (associated(CS%var_ptr1d(m)%p))  then
                 call fms2_read_data(CS%fileObj, varname, CS%var_ptr1d(m)%p)

                 if (is_there_a_checksum) checksum_data = mpp_chksum(CS%var_ptr1d(m)%p)
              elseif (associated(CS%var_ptr0d(m)%p)) then ! Read a scalar
                 call fms2_read_data(CS%fileObj,varname,  CS%var_ptr0d(m)%p)
                 if (is_there_a_checksum) checksum_data = mpp_chksum(CS%var_ptr0d(m)%p,pelist=(/mpp_pe()/))
                   
              elseif (associated(CS%var_ptr2d(m)%p)) then  ! Read a 2d array
                 if (pos /= 0) then
                    call fms2_read_data(CS%fileObj, varname, CS%var_ptr2d(m)%p) ! domain-decomposed read
                 !else ! This array is not domain-decomposed.  This variant may be under-tested.
                 !   call read_data(unit_path(n), varname, CS%var_ptr2d(m)%p, &
                 !          no_domain=.true., timelevel=1)
                 endif
                 if (is_there_a_checksum) checksum_data = mpp_chksum(CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL))
              elseif (associated(CS%var_ptr3d(m)%p)) then  ! Read a 3d array.
                 if (pos /= 0) then
                    call fms2_read_data(CS%fileObj, varname,  CS%var_ptr3d(m)%p) ! domain-decomposed read

                 !else ! This array is not domain-decomposed.  This variant may be under-tested.
                 !   call read_data(unit_path(n), varname, CS%var_ptr3d(m)%p, &
                 !                  no_domain=.true., timelevel=1)
                 endif
                 if (is_there_a_checksum) checksum_data = mpp_chksum(CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:))
              elseif (associated(CS%var_ptr4d(m)%p)) then  ! Read a 4d array.
                 if (pos /= 0) then
                     call fms2_read_data(CS%fileObj, varname, CS%var_ptr4d(m)%p) ! domain-decomposed read
                 !else ! This array is not domain-decomposed.  This variant may be under-tested.
                 !    call read_data(unit_path(n), varname, CS%var_ptr4d(m)%p, &
                 !           no_domain=.true., timelevel=1)
                 endif

                 if (is_there_a_checksum) then 
                    checksum_data = mpp_chksum(CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:))
                 endif                 
              else
                 call MOM_error(FATAL, "MOM_restart restore_state: No pointers set for "//trim(varname))
              endif

              if (is_root_pe() .and. is_there_a_checksum .and. (checksum_file(1) /= checksum_data)) then
                 write (mesg,'(a,Z16,a,Z16,a)') "Checksum of input field "// trim(varname)//" ",checksum_data,&
                        " does not match value ", checksum_file(1), &
                        " stored in "//trim(unit_path(n)//"." )
                 call MOM_error(FATAL, "MOM_restart(restore_state): "//trim(mesg) )
              endif

              CS%restart_field(m)%initialized = .true.
               
              exit ! Start search for next restart variable.
           endif
        enddo
        if (i>nvar) missing_fields = missing_fields+1
     enddo

     deallocate(fields)
     deallocate(variable_names)

     call fms2_close_file(CS%fileObj)

     if (missing_fields == 0) exit   
  enddo

! Check whether any mandatory fields have not been found.
  CS%restart = .true.
  do m=1,CS%novars
     if (.not.(CS%restart_field(m)%initialized)) then
        CS%restart = .false.
        if (CS%restart_field(m)%mand_var) then
           call MOM_error(FATAL,"MOM_restart: Unable to find mandatory variable " &
                       //trim(CS%restart_field(m)%var_name)//" in restart files.")
        endif
     endif
  enddo

end subroutine restore_state

!> restart_files_exist determines whether any restart files exist.
function restart_files_exist(filename, directory, G, CS)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files.
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init.
  logical :: restart_files_exist                  !< The function result, which indicates whether
                                                  !! any of the explicitly or automatically named
                                                  !! restart files exist in directory.
  integer :: num_files

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "restart_files_exist: Module must be initialized before it is used.")

  if ((LEN_TRIM(filename) == 1) .and. (filename(1:1) == 'F')) then
    num_files = open_restart_units('r', directory, G, CS)
  else
    num_files = open_restart_units(filename, directory, G, CS)
  endif
  restart_files_exist = (num_files > 0)

end function restart_files_exist

!> determine_is_new_run determines from the value of filename and the existence
!! automatically named restart files in directory whether this would be a new,
!! and as a side effect stores this information in CS.
function determine_is_new_run(filename, directory, G, CS) result(is_new_run)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files.
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init.
  logical :: is_new_run                           !< The function result, which indicates whether
                                                  !! this is a new run, based on the value of
                                                  !! filename and whether restart files exist.

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "determine_is_new_run: Module must be initialized before it is used.")
  if (LEN_TRIM(filename) > 1) then
    CS%new_run = .false.
  elseif (LEN_TRIM(filename) == 0) then
    CS%new_run = .true.
  elseif (filename(1:1) == 'n') then
    CS%new_run = .true.
  elseif (filename(1:1) == 'F') then
    CS%new_run = (open_restart_units('r', directory, G, CS) == 0)
  else
    CS%new_run = .false.
  endif

  CS%new_run_set = .true.
  is_new_run = CS%new_run
end function determine_is_new_run

!> is_new_run returns whether this is going to be a new run based on the
!! information stored in CS by a previous call to determine_is_new_run.
function is_new_run(CS)
  type(MOM_restart_CS),  pointer :: CS !< The control structure returned by a previous
                                       !! call to restart_init.
  logical :: is_new_run                !< The function result, which indicates whether
                                       !! this is a new run, based on the value of
                                       !! filename and whether restart files exist.

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "is_new_run: Module must be initialized before it is used.")
  if (.not.CS%new_run_set) call MOM_error(FATAL, "MOM_restart " // &
      "determine_is_new_run must be called for a restart file before is_new_run.")

  is_new_run = CS%new_run
end function is_new_run

!> open_restart_units determines the number of existing restart files and optionally opens
!! them and returns unit ids, paths and whether the files are global or spatially decomposed.
function open_restart_units(filename, directory, G, CS, units, file_paths, &
                            global_files) result(num_files)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files.
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(ocean_grid_type), intent(in)  :: G         !< The ocean's grid structure
  type(MOM_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to restart_init.
  integer, dimension(:), &
               optional, intent(out) :: units     !< The mpp units of all opened files.
  character(len=*), dimension(:), &
               optional, intent(out) :: file_paths   !< The full paths to open files.
  logical, dimension(:), &
               optional, intent(out) :: global_files !< True if a file is global.

  integer :: num_files  !< The number of files (both automatically named restart
                        !! files and others explicitly in filename) that have been opened.

!  This subroutine reads the model state from previously
!  generated files.  All restart variables are read from the first
!  file in the input filename list in which they are found.

  ! Local variables
  character(len=256) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=256) :: fname     ! The name of the current file.
  character(len=8)   :: suffix    ! A suffix (like "_2") that is added to any
                                  ! additional restart files.
! character(len=256) :: mesg      ! A message for warnings.
  integer :: num_restart     ! The number of restart files that have already
                             ! been opened.
  integer :: start_char      ! The location of the starting character in the
                             ! current file name.
  integer :: n, m, err, length


  logical :: fexists
  character(len=32) :: filename_appendix = '' !fms appendix to filename for ensemble runs
  character(len=80) :: restartname

  if (.not.associated(CS)) call MOM_error(FATAL, "MOM_restart " // &
      "open_restart_units: Module must be initialized before it is used.")

! Get NetCDF ids for all of the restart files.
  num_restart = 0 ; n = 1 ; start_char = 1
  do while (start_char <= len_trim(filename) )
    do m=start_char,len_trim(filename)
      if (filename(m:m) == ' ') exit
    enddo
    fname = filename(start_char:m-1)
    start_char = m
    do while (start_char <= len_trim(filename))
      if (filename(start_char:start_char) == ' ') then
        start_char = start_char + 1
      else
        exit
      endif
    enddo

    if ((fname(1:1)=='r') .and. ( len_trim(fname) == 1)) then
      err = 0
      if (num_restart > 0) err = 1 ! Avoid going through the file list twice.
      do while (err == 0)
        restartname = trim(CS%restartfile)

       !query fms_io if there is a filename_appendix (for ensemble runs)
       call get_filename_appendix(filename_appendix)
       if (len_trim(filename_appendix) > 0) then
         length = len_trim(restartname)
         if (restartname(length-2:length) == '.nc') then
         restartname = restartname(1:length-3)//'.'//trim(filename_appendix)//'.nc'
         else
           restartname = restartname(1:length)  //'.'//trim(filename_appendix)
         endif
        endif
        filepath = trim(directory) // trim(restartname)

        if (num_restart < 10) then
          write(suffix,'("_",I1)') num_restart
        else
          write(suffix,'("_",I2)') num_restart
        endif
        if (num_restart > 0) filepath = trim(filepath) // suffix

        ! if (.not.file_exists(filepath)) &
          filepath = trim(filepath)//".nc"

        num_restart = num_restart + 1
        inquire(file=filepath, exist=fexists)
        if (fexists) then
          if (present(units)) &
            call open_file(units(n), trim(filepath), READONLY_FILE, NETCDF_FILE, &
                           threading = MULTIPLE, fileset = SINGLE_FILE)
          if (present(global_files)) global_files(n) = .true.
        elseif (CS%parallel_restartfiles) then
          ! Look for decomposed files using the I/O Layout.
          fexists = file_exists(filepath, G%Domain)
          if (fexists .and. (present(units))) &
            call open_file(units(n), trim(filepath), READONLY_FILE, NETCDF_FILE, &
                           domain=G%Domain%mpp_domain)
          if (fexists .and. present(global_files)) global_files(n) = .false.
        endif

        if (fexists) then
          if (present(file_paths)) file_paths(n) = filepath
          n = n + 1
          if (is_root_pe() .and. (present(units))) &
            call MOM_error(NOTE, "MOM_restart: MOM run restarted using : "//trim(filepath))
        else
          err = 1 ; exit
        endif
      enddo ! while (err == 0) loop
    else
      filepath = trim(directory)//trim(fname)
      inquire(file=filepath, exist=fexists)
      if (.not. fexists) filepath = trim(filepath)//".nc"

      inquire(file=filepath, exist=fexists)
      if (fexists) then
        if (present(units)) &
          call open_file(units(n), trim(filepath), READONLY_FILE, NETCDF_FILE, &
                       threading = MULTIPLE, fileset = SINGLE_FILE)
        if (present(global_files)) global_files(n) = .true.
        if (present(file_paths)) file_paths(n) = filepath
        n = n + 1
        if (is_root_pe() .and. (present(units))) &
          call MOM_error(NOTE,"MOM_restart: MOM run restarted using : "//trim(filepath))
      else
        if (present(units)) &
          call MOM_error(WARNING,"MOM_restart: Unable to find restart file : "//trim(filepath))
      endif

    endif
  enddo ! while (start_char < strlen(filename)) loop
  num_files = n-1

end function open_restart_units

!> Initialize this module and set up a restart control structure.
subroutine restart_init(param_file, CS, restart_root)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object that is allocated here
  character(len=*), optional, &
                         intent(in) :: restart_root !< A filename root that overrides the value
                                          !! set by RESTARTFILE to enable the use of this module by
                                          !! other components than MOM.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_restart"   ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "restart_init called with an associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "PARALLEL_RESTARTFILES", &
                                CS%parallel_restartfiles, &
                 "If true, each processor writes its own restart file, \n"//&
                 "otherwise a single restart file is generated", &
                 default=.false.)

  if (present(restart_root)) then
    CS%restartfile = restart_root
    call log_param(param_file, mdl, "RESTARTFILE from argument", CS%restartfile)
  else
    call get_param(param_file, mdl, "RESTARTFILE", CS%restartfile, &
                 "The name-root of the restart file.", default="MOM.res")
  endif
  call get_param(param_file, mdl, "LARGE_FILE_SUPPORT", CS%large_file_support, &
                 "If true, use the file-size limits with NetCDF large \n"//&
                 "file support (4Gb), otherwise the limit is 2Gb.", &
                 default=.true.)
  call get_param(param_file, mdl, "MAX_FIELDS", CS%max_fields, &
                 "The maximum number of restart fields that can be used.", &
                 default=100)
  call get_param(param_file, mdl, "RESTART_CHECKSUMS_REQUIRED", CS%checksum_required, &
                 "If true, require the restart checksums to match and error out otherwise. \n"//&
                 "Users may want to avoid this comparison if for example the restarts are  \n"//&
                 "made from a run with a different mask_table than the current run,  \n"//&
                 "in which case the checksums will not match and cause crash.",&
                 default=.true.)

  allocate(CS%restart_field(CS%max_fields))
  allocate(CS%var_ptr0d(CS%max_fields))
  allocate(CS%var_ptr1d(CS%max_fields))
  allocate(CS%var_ptr2d(CS%max_fields))
  allocate(CS%var_ptr3d(CS%max_fields))
  allocate(CS%var_ptr4d(CS%max_fields))

end subroutine restart_init

!> Indicate that all variables have now been registered.
subroutine restart_init_end(CS)
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object

  if (associated(CS)) then
    if (CS%novars == 0) call restart_end(CS)
  endif

end subroutine restart_init_end

!> Deallocate memory associated with a MOM_restart_CS variable.
subroutine restart_end(CS)
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object

  if (associated(CS%restart_field)) deallocate(CS%restart_field)
  if (associated(CS%var_ptr0d)) deallocate(CS%var_ptr0d)
  if (associated(CS%var_ptr1d)) deallocate(CS%var_ptr1d)
  if (associated(CS%var_ptr2d)) deallocate(CS%var_ptr2d)
  if (associated(CS%var_ptr3d)) deallocate(CS%var_ptr3d)
  if (associated(CS%var_ptr4d)) deallocate(CS%var_ptr4d)
  deallocate(CS)

end subroutine restart_end

subroutine restart_error(CS)
  type(MOM_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object

  character(len=16)  :: num  ! String for error messages

  if (CS%novars > CS%max_fields) then
    write(num,'(I0)') CS%novars
    call MOM_error(FATAL,"MOM_restart: Too many fields registered for " // &
           "restart.  Set MAX_FIELDS to be at least " // &
           trim(adjustl(num)) // " in the MOM input file.")
  else
    call MOM_error(FATAL,"MOM_restart: Unspecified fatal error.")
  endif
end subroutine restart_error

!> Return bounds for computing checksums to store in restart files
subroutine get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL)
  type(ocean_grid_type), intent(in)  :: G !< The ocean's grid structure
  integer,               intent(in)  :: pos !< An integer indicating staggering of variable
  integer,               intent(out) :: isL !< i-start for checksum
  integer,               intent(out) :: ieL !< i-end for checksum
  integer,               intent(out) :: jsL !< j-start for checksum
  integer,               intent(out) :: jeL !< j-end for checksum

  ! Regular non-symmetric compute domain
  isL = G%isc-G%isd+1
  ieL = G%iec-G%isd+1
  jsL = G%jsc-G%jsd+1
  jeL = G%jec-G%jsd+1

  ! Expand range east or south for symmetric arrays
  if (G%symmetric) then
    if ((pos == EAST_FACE) .or. (pos == CORNER)) then ! For u-, q-points only
      if (G%idg_offset == 0) isL = isL - 1 ! include western edge in checksums only for western PEs
    endif
    if ((pos == NORTH_FACE) .or. (pos == CORNER)) then ! For v-, q-points only
      if (G%jdg_offset == 0) jsL = jsL - 1 ! include western edge in checksums only for southern PEs
    endif
  endif

end subroutine get_checksum_loop_ranges

subroutine get_horizontal_grid_coordinates(hor_grid,use_lath,use_lonh, &
                                           use_latq,use_lonq, gridPosition)
  character(len=*), intent(in) :: hor_grid !< name of the horizontal grid 
  logical, intent(out) ::  use_lath, use_lonh, use_latq, use_lonq !< flags for lath, lonh,latq, lonq
  integer, intent(out) :: gridPosition !< integer corresponding to the grid position
  use_lath = .false.
  use_lonh = .false.
  use_latq = .false.
  use_lonh = .false.

  select case (hor_grid)
     case ('h') ; use_lath = .true. ; use_lonh = .true.; gridPosition = CENTER
     case ('q') ; use_latq = .true. ; use_lonq = .true.; gridPosition = CORNER
     case ('u') ; use_lath = .true. ; use_lonq = .true.; gridPosition = EAST_FACE
     case ('v') ; use_latq = .true. ; use_lonh = .true.; gridPosition = NORTH_FACE
     case ('T')  ; use_lath = .true. ; use_lonh = .true.; gridPosition = CENTER
     case ('Bu') ; use_latq = .true. ; use_lonq = .true.; gridPosition = CORNER
     case ('Cu') ; use_lath = .true. ; use_lonq = .true.; gridPosition = EAST_FACE
     case ('Cv') ; use_latq = .true. ; use_lonh = .true.; gridPosition = NORTH_FACE
     case ('1') ; gridPosition = 0 
     case default
        call MOM_error(FATAL, "MOM_restart:get_horizontal_grid_coordinates "//&
                        "Unrecognized hor_grid argument "//trim(hor_grid))
  end select
end subroutine get_horizontal_grid_coordinates

subroutine get_vertical_grid_coordinates(z_grid,use_layer,use_int)
  character(len=*), intent(in) :: z_grid
  logical, intent(out) ::  use_layer, use_int

  use_layer=.false.
  use_int=.false.
  select case (trim(z_grid))
     case ('L') ; use_layer = .true.
     case ('i') ; use_int = .true.
     case ('1') ! Do nothing.
     case default
       call MOM_error(FATAL, "MOM_restart: get_vertical_grid_coordinates: "//&
                        " has unrecognized z_grid argument"//trim(z_grid))
  end select
end subroutine get_vertical_grid_coordinates

subroutine get_time_coordinates(t_grid_in, use_time, use_periodic)
  character(len=*), intent(in) :: t_grid_in !< 's', 'a', 'm' for time,
                                            !< 'p' for periodic, '1' for no time axis
  logical, intent(out) :: use_time, use_periodic
  
  ! Local
  character(len=8) :: t_grid, t_grid_read

  use_time = .false.
  use_periodic = .false.

  t_grid = adjustl(t_grid_in)

  select case (t_grid(1:1))
     case ('s', 'a', 'm') ; use_time = .true.
     case ('p')
        if (len_trim(t_grid(2:8)) <= 0) then
            call MOM_error(FATAL,"MOM_restart:get_time_coordinates: "//&
                           "No periodic axis length was specified in "//trim(t_grid))
        endif
        use_periodic = .true.
     case ('1') ! Do nothing.
     case default
        call MOM_error(WARNING, "MOM_restart: get_time_coordinates: "//&
                       "Unrecognized t_grid "//trim(t_grid))
  end select
end subroutine get_time_coordinates

subroutine check_for_restart_axis(fileObject,axisName,axisLength)
   type(FmsNetcdfDomainFile_t), intent(inout) :: fileObject !< file object returned by prior call to fms2_open_file
   character(len=*), intent(in) :: axisName ! name of the restart file axis to register to file
   integer, intent(in) :: axisLength ! length of axis/dimension (only needed for Z and Time)
   ! local  
   logical :: axisExists 
   
   axisExists= .false.
 
   axisExists = fms2_dimension_exists(fileObject,axisName)
   if (.not. (axisExists)) then
      select case (trim(axisName))
         case ('latq'); call fms2_register_axis(fileObject,'latq','y')
         case ('lath'); call fms2_register_axis(fileObject,'lath','y') 
         case ('lonq'); call fms2_register_axis(fileObject,'lonq','x') 
         case ('lonh'); call fms2_register_axis(fileObject,'lath','x')
         case ('Layer'); call fms2_register_axis(fileObject,'Layer',axisLength)
         case ('Interface'); call fms2_register_axis(fileObject,'Interface',axisLength)
         case ('Time'); call fms2_register_axis(fileObject,'Time', axisLength)
         case ('Period'); call fms2_register_axis(fileObject,'Period',axisLength)
      end select
   endif
  
end subroutine check_for_restart_axis

!> Define the axis variable attributes, and write the axis data
!! to the restart file
subroutine write_axis_data(fileObject, axisName, G, dG, GV, timeunit, & 
                           restart_time_in_days, t_grid_in)
  type(FmsNetcdfDomainFile_t), intent(inout) :: fileObject !< file object returned by prior call to fms2_open_file
  character(len=*), intent(in) :: axisName !< Name of the axis
  type(ocean_grid_type), optional, intent(in) :: G !< ocean horizontal grid structure; G or dG
                                                    !! is required if the new file uses any
                                                    !! horizontal grid axes.
  type(dyn_horgrid_type), optional, intent(in) :: dG !< dynamic horizontal grid structure; G or dG
                                                      !! is required if the new file uses any
  type(verticalGrid_type), optional, intent(in) :: GV !< ocean vertical grid structure, which is
                                                     !! required if the new file uses any
                                                     !! vertical grid axes.
  real, optional, intent(in) :: timeunit !< length of the units for time [s]. The
                                          !! default value is 86400.0, for 1 day.
  real, optional, intent(in) :: restart_time_in_days !< restart time in days

  character(len=*), optional, intent(in) :: t_grid_in
  ! local
  type(MOM_domain_type), pointer :: Domain => NULL()
  type(domain1d) :: x_domain, y_domain
  logical :: domain_set = .false.
  real, pointer, dimension(:) :: &
     gridLatT => NULL(), & ! The latitude or longitude of T or B points for
     gridLatB => NULL(), & ! the purpose of labeling the output axes.
     gridLonT => NULL(), & 
     gridLonB => NULL()
  
  integer :: num_axes, num_periods, var_periods, k
  integer :: isg, ieg, jsg, jeg, IsgB, IegB, JsgB, JegB
  character(len=40) :: axis_units,x_axis_units, y_axis_units, time_units
  real, dimension(:), allocatable :: period_val
  character(len=8) :: t_grid, t_grid_read
  character(len=50) :: long_name

  if ((trim(axisName) == "Period") .and. .not.(present(t_grid_in))) then 
      call MOM_error(FATAL,"MOM_restart::write_axis_data: NO argument passed for 't_grid_in', "//&
                     " which is required to determine the periodic time axis.")
  endif

  if (present(G)) then
     domain_set = .true. ; Domain => G%Domain
     gridLatT => G%gridLatT ; gridLatB => G%gridLatB
     gridLonT => G%gridLonT ; gridLonB => G%gridLonB
     x_axis_units = G%x_axis_units ; y_axis_units = G%y_axis_units
     isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
     IsgB = G%IsgB ; IegB = G%IegB ; JsgB = G%JsgB ; JegB = G%JegB
  elseif (present(dG)) then
     domain_set = .true. ; Domain => dG%Domain
     gridLatT => dG%gridLatT ; gridLatB => dG%gridLatB
     gridLonT => dG%gridLonT ; gridLonB => dG%gridLonB
     x_axis_units = dG%x_axis_units ; y_axis_units = dG%y_axis_units
     isg = dG%isg ; ieg = dG%ieg ; jsg = dG%jsg ; jeg = dG%jeg
     IsgB = dG%IsgB ; IegB = dG%IegB ; JsgB = dG%JsgB ; JegB = dG%JegB
  endif
  
  select case (trim(axisName))
     case ('latq')
        call fms2_write_data(fileObject,axisName,gridLatB(JsgB:JegB))
        long_name = 'Latitude'
        axis_units = y_axis_units
     case ('lath')
        call fms2_write_data(fileObject,axisName,gridLatT(jsg:jeg))
        long_name = 'Latitude'
        axis_units = y_axis_units
     case ('lonq')
        call fms2_write_data(fileObject,axisName,gridLonB(IsgB:IegB))
        long_name = 'Longitude'
        axis_units = x_axis_units
     case ('lonh')
        call fms2_write_data(fileObject,axisName,gridLonT(isg:ieg))
        long_name = 'Longitude'
        axis_units = x_axis_units
     case ('Layer')
        if (.not. (present(GV))) then
           call MOM_error(FATAL,"MOM_restart::write_axis_data: No argument passed for 'GV', "//&
                          " which is required to write the Layer axis data.")
        endif
        call fms2_write_data(fileObject,axisName,GV%sLayer(1:GV%ke))
        long_name = 'Layer'
        axis_units = trim(GV%zAxisUnits)
     case ('Interface')
        if (.not. (present(GV))) then
           call MOM_error(FATAL,"MOM_restart::write_axis_data: No argument passed for 'GV', "//&
                          " which is required to write the Interface axis data.")
        endif
        call fms2_write_data(fileObject,axisName,GV%sInterface(1:GV%ke+1))
        long_name="Interface "//trim(GV%zAxisLongName)
        axis_units = trim(GV%zAxisUnits)
     case ('Time') 
        if (present(timeunit)) then
        ! Set appropriate units, depending on the value.
           if (timeunit < 0.0) then
              time_units = "days" ! The default value.
           elseif ((timeunit >= 0.99) .and. (timeunit < 1.01)) then
              time_units = "seconds"
           elseif ((timeunit >= 3599.0) .and. (timeunit < 3601.0)) then
              time_units = "hours"
           elseif ((timeunit >= 86399.0) .and. (timeunit < 86401.0)) then
              time_units = "days"
           elseif ((timeunit >= 3.0e7) .and. (timeunit < 3.2e7)) then
              time_units = "years"
           else
              write(time_units,'(es8.2," s")') timeunit
           endif
           axis_units = time_units
        else
           axis_units = "days"
        endif
        if (present(restart_time_in_days)) then
           call fms2_write_data(fileObject,axisName,restart_time_in_days)
        else
           call fms2_write_data(fileObject,axisName,'1.0')
        endif
        long_name = "Time"
     case ('Period')
        t_grid = adjustl(t_grid_in)
        if (len_trim(t_grid(2:8)) <= 0) then
           call MOM_error(FATAL,"MOM_restart::write_axis_data: No periodic axis"//&
                          "length was specified in "//trim(t_grid)//".")
        endif
        num_periods = 0
        var_periods = -9999999
        t_grid_read = adjustl(t_grid(2:8))
        read(t_grid_read,*) var_periods
        if (var_periods == -9999999) then
           call MOM_error(FATAL, "MOM_restart::write_axis_data: Failed to "//&
                          "read the number of periods from "//trim(t_grid)// ".")
        endif
        if (var_periods < 1) then
           call MOM_error(FATAL, "MOM_restart::write_axis_data: The number of periods "//&
                          trim(t_grid)//" must be positive.")
        endif
        if ((num_periods > 0) .and. (var_periods /= num_periods)) then
           call MOM_error(FATAL, "MOM_restart::write_axis_data: Only one value"//&
                          " of the number of periods can be used for "//trim(t_grid)//".")
        endif

        num_periods = var_periods
        if (num_periods <= 1) then
           call MOM_error(FATAL, "MOM_restart::write_axis_data: num_periods must be at least 1.")
        endif
        ! Define a periodic axis with unit labels.
        allocate(period_val(num_periods))
        do k=1,num_periods
           period_val(k) = real(k)
        enddo
        ! call mpp_write_meta(unit, axis_periodic, name="Period", units="nondimensional", &
        !longname="Periods for cyclical varaiables", cartesian= 't', data=period_val)
        call fms2_write_data(fileObject,axisName,period_val)
        deallocate(period_val)
        axis_units = "nondimensional"
        long_name="Periods for cyclical varaiables"
   end select
   ! write attributes
   call fms2_register_variable_attribute(fileObject,axisName,"units",axis_units)
   call fms2_register_variable_attribute(fileObject,axisName,"longname",long_name)

end subroutine write_axis_data

!subroutine write_variable_axes(var_name, axes, hor_grid,z_grid,t_grid)
!   character(len=*), intent(in) :: var_name !< name of the variable
!   type(axistype), intent(out)  :: axes(4) 
!   character(len=*), intent(in), optional :: hor_grid
!   character(len=*), intent(in), optional :: z_grid
!   character(len=*), intent(in), optional :: t_grid
   ! local
!   type(axistype) :: axis_lath, axis_latq, axis_lonh, axis_lonq
!   type(axistype) :: axis_layer, axis_int, axis_time, axis_periodic
!   integer :: numaxes = 0
!   if (present(hor_grid) then
!      select case (hor_grid)
!         case ('h')  ; num_axes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_lath
!         case ('q')  ; num_axes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_latq
!         case ('u')  ; num_axes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_lath
!         case ('v')  ; num_axes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_latq
!         case ('T')  ; num_axes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_lath
!         case ('Bu') ; num_axes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_latq
!         case ('Cu') ; num_axes = 2 ; axes(1) = axis_lonq ; axes(2) = axis_lath
!!         case ('Cv') ; num_axes = 2 ; axes(1) = axis_lonh ; axes(2) = axis_latq
!         case ('1') ! Do nothing.
!         case default
!            call MOM_error(WARNING, "MOM_restart:get_variable_axes: "//trim(var_name)//&
!                        " has unrecognized hor_grid "//trim(hor_grid))
!      end select
!   endif
!   if (present(z_grid))
!     select case (z_grid)
!      case ('L') ; num_axes = num_axes+1 ; axes(num_axes) = axis_layer
!      case ('i') ; num_axes = num_axes+1 ; axes(num_axes) = axis_int
!      case ('1') ! Do nothing.
!      case default
!        call MOM_error(FATAL, "MOM_restart:get_variable_axes: "//trim(var_name)//&
!                        " has unrecognized z_grid "//trim(z_grid))
!    end select
!    t_grid = adjustl(vars(k)%t_grid)
!    select case (t_grid(1:1))
!      case ('s', 'a', 'm') ; num_axes = num_axes+1 ; axes(num_axes) = axis_time
!      case ('p')           ; num_axes = num_axes+1 ; axes(num_axes) = axis_periodic
!      case ('1') ! Do nothing.
!      case default
!        call MOM_error(WARNING, "MOM_io create_file: "//trim(vars(k)%name)//&
!                        " has unrecognized t_grid "//trim(vars(k)%t_grid))
!    end select
!!    pack = 1

!    if (present(checksums)) then
!       call mpp_write_meta(unit, fields(k), axes(1:numaxes), vars(k)%name, vars(k)%units, &
!           vars(k)%longname, pack = pack, checksum=checksums(k,:))
!    else
!       call mpp_write_meta(unit, fields(k), axes(1:numaxes), vars(k)%name, vars(k)%units, &
!           vars(k)%longname, pack = pack)
!    endif
!enddo
!end subroutine write_variable_axes


end module MOM_restart
