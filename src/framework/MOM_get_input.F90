!> \brief Reads the only Fortran name list needed to boot-strap the model.
!!
!! The name list parameters indicate which directories to use for
!! certain types of input and output, and which files to look in for
!! the full parsable input parameter file(s).
module MOM_get_input

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : open_param_file, param_file_type
use MOM_io, only : file_exists, close_file, slasher, ensembler
use MOM_io, only : open_namelist_file, check_nml_error

implicit none ; private

public get_MOM_input

!> Container for paths and parameter file names.
type, public :: directories
  character(len=240) :: &
    restart_input_dir = ' ',& !< The directory to read restart and input files.
    restart_output_dir = ' ',&!< The directory into which to write restart files.
    output_directory = ' '    !< The directory to use to write the model output.
  character(len=2048) :: &
    input_filename  = ' '     !< A string that indicates the input files or how
                              !! the run segment should be started.
end type directories

contains

!> Get the names of the I/O directories and initialization file.
!! Also calls the subroutine that opens run-time parameter files.
subroutine get_MOM_input(param_file, dirs, check_params, default_input_filename, ensemble_num)
  type(param_file_type), optional, intent(out) :: param_file   !< A structure to parse for run-time parameters.
  type(directories),     optional, intent(out) :: dirs         !< Container for paths and parameter file names.
  logical,               optional, intent(in)  :: check_params !< If present and False will stop error checking for
                                                               !! run-time parameters.
  character(len=*),      optional, intent(in)  :: default_input_filename !< If present, is the value assumed for
                                                               !! input_filename if input_filename is not listed
                                                               !! in the namelist MOM_input_nml.
  integer, optional, intent(in) :: ensemble_num !< The ensemble id of the current member
  ! Local variables
  integer, parameter :: npf = 5 ! Maximum number of parameter files

  character(len=240) :: &
    parameter_filename(npf), & ! List of files containing parameters.
    output_directory,        & ! Directory to use to write the model output.
    restart_input_dir,       & ! Directory for reading restart and input files.
    restart_output_dir         ! Directory into which to write restart files.
  character(len=2048) :: &
    input_filename             ! A string that indicates the input files or how
                               ! the run segment should be started.
  character(len=240) :: output_dir
  integer :: unit, io, ierr, valid_param_files

  namelist /MOM_input_nml/ output_directory, input_filename, parameter_filename, &
                           restart_input_dir, restart_output_dir

  ! Default values in case parameter is not set in file input.nml
  parameter_filename(:) = ' '
  output_directory = ' '
  restart_input_dir = ' '
  restart_output_dir = ' '
  input_filename  = ' '
  if (present(default_input_filename)) input_filename = trim(default_input_filename)

  ! Open namelist
  if (file_exists('input.nml')) then
    unit = open_namelist_file(file='input.nml')
  else
    call MOM_error(FATAL,'Required namelist file input.nml does not exist.')
  endif

  ! Read namelist parameters
  ierr=1 ; do while (ierr /= 0)
    read(unit, nml=MOM_input_nml, iostat=io, end=10)
    ierr = check_nml_error(io, 'MOM_input_nml')
  enddo
10 call close_file(unit)

  ! Store parameters in container
  if (present(dirs)) then
    if (present(ensemble_num)) then
      dirs%output_directory = slasher(ensembler(output_directory,ensemble_num))
      dirs%restart_output_dir = slasher(ensembler(restart_output_dir,ensemble_num))
      dirs%restart_input_dir = slasher(ensembler(restart_input_dir,ensemble_num))
      dirs%input_filename = ensembler(input_filename,ensemble_num)
    else
      dirs%output_directory = slasher(ensembler(output_directory))
      dirs%restart_output_dir = slasher(ensembler(restart_output_dir))
      dirs%restart_input_dir = slasher(ensembler(restart_input_dir))
      dirs%input_filename = ensembler(input_filename)
    endif
  endif

  ! Open run-time parameter file(s)
  if (present(param_file)) then
    output_dir = slasher(ensembler(output_directory))
    valid_param_files = 0
    do io = 1, npf
      if (len_trim(trim(parameter_filename(io))) > 0) then
        if (present(ensemble_num)) then
          call open_param_file(ensembler(parameter_filename(io),ensemble_num), param_file, &
               check_params, doc_file_dir=output_dir)
        else
          call open_param_file(ensembler(parameter_filename(io)), param_file, &
               check_params, doc_file_dir=output_dir)
        endif
        valid_param_files = valid_param_files + 1
      endif
    enddo
    if (valid_param_files == 0) call MOM_error(FATAL, "There must be at "//&
         "least 1 valid entry in input_filename in MOM_input_nml in input.nml.")
  endif

end subroutine get_MOM_input

end module MOM_get_input
