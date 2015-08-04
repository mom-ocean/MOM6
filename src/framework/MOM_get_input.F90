module MOM_get_input
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, April 2013                                     *
!*                                                                     *
!*    The subroutine in this file reads the MOM6 namelist input, which *
!*  indicates which directories to use for certain types of input and  *
!*  output, and where to look for the full parsable input file(s).     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : open_param_file, param_file_type
use MOM_io, only : file_exists, close_file, slasher, ensembler
use MOM_io, only : open_namelist_file, check_nml_error

implicit none ; private

public Get_MOM_Input

! This structure is to simplify communication with the calling code.

type, public :: directories
  character(len=240) :: &
    restart_input_dir = ' ',& ! The directory to read restart and input files.
    restart_output_dir = ' ',&! The directory into which to write restart files.
    output_directory = ' ', & ! The directory to use to write the model output.
    input_filename  = ' '     ! A string that indicates the input files or how
                              ! the run segment should be started.
end type directories

contains

subroutine Get_MOM_Input(param_file, dirs, check_params)
  type(param_file_type), optional, intent(out) :: param_file
  type(directories),     optional, intent(out) :: dirs
  logical,               optional, intent(in)  :: check_params

!    See if the run is to be started from saved conditions, and get  !
!  the names of the I/O directories and initialization file.  This   !
!  subroutine also calls the subroutine that allows run-time changes !
!  in parameters.                                                    !
  integer, parameter :: npf = 5 ! Maximum number of parameter files
  character(len=240) :: &
    parameter_filename(npf) = ' ', & ! List of files containing parameters.
    output_directory = ' ', &   ! Directory to use to write the model output.
    restart_input_dir = ' ', &  ! Directory for reading restart and input files.
    restart_output_dir = ' ', & ! Directory into which to write restart files.
    input_filename  = ' '       ! A string that indicates the input files or how
                                ! the run segment should be started.
  character(len=240) :: output_dir
  integer :: unit, io, ierr, valid_param_files

  namelist /MOM_input_nml/ output_directory, input_filename, parameter_filename, &
                           restart_input_dir, restart_output_dir

  if (file_exists('input.nml')) then
    unit = open_namelist_file(file='input.nml')
  else
    call MOM_error(FATAL,'Required namelist file input.nml does not exist.')
  endif

  ierr=1 ; do while (ierr /= 0)
    read(unit, nml=MOM_input_nml, iostat=io, end=10)
    ierr = check_nml_error(io, 'MOM_input_nml')
  enddo
10 call close_file(unit)

  if (present(dirs)) then
    dirs%output_directory = slasher(ensembler(output_directory))
    dirs%restart_output_dir = slasher(ensembler(restart_output_dir))
    dirs%restart_input_dir = slasher(ensembler(restart_input_dir))
    dirs%input_filename = ensembler(input_filename)
  endif

  if (present(param_file)) then
    output_dir = slasher(ensembler(output_directory))
    valid_param_files = 0
    do io = 1, npf
      if (len_trim(trim(parameter_filename(io))) > 0) then
        call open_param_file(ensembler(parameter_filename(io)), param_file, &
                             check_params, doc_file_dir=output_dir)
        valid_param_files = valid_param_files + 1
      endif
    enddo
    if (valid_param_files == 0) call MOM_error(FATAL, "There must be at "//&
         "least 1 valid entry in input_filename in MOM_input_nml in input.nml.")
  endif

end subroutine Get_MOM_Input

end module MOM_get_input
