program test_MOM_file_parser

use MPI
use MOM_domains, only : MOM_infra_init
use MOM_domains, only : MOM_infra_end
use MOM_file_parser_tests, only : run_file_parser_tests

implicit none

integer, parameter :: comm = MPI_COMM_WORLD
integer, parameter :: root = 0
integer :: rank
logical :: file_exists_on_rank
logical :: input_nml_exists, MOM_input_exists
integer :: io_unit
logical :: is_open, is_file
integer :: rc

! NOTE: Bootstrapping requires external MPI configuration.
!   - FMS initialization requires the presence of input.nml
!   - MOM initialization requires MOM_input (if unspecificed by input.nml)
!   - Any MPI-based I/O prior to MOM and FMS init will MPI initialization
! Thus, we need to do some minimal MPI setup.
call MPI_Init(rc)
call MPI_Comm_rank(comm, rank, rc)

inquire(file='input.nml', exist=file_exists_on_rank)
call MPI_Reduce(file_exists_on_rank, input_nml_exists, 1, MPI_LOGICAL, &
                MPI_LOR, root, comm, rc)

inquire(file='MOM_input', exist=file_exists_on_rank)
call MPI_Reduce(file_exists_on_rank, MOM_input_exists, 1, MPI_LOGICAL, &
                MPI_LOR, root, comm, rc)

if (rank == root) then
  ! Abort if at least one rank sees either input.nml or MOM_input
  if (input_nml_exists) error stop "Remove existing 'input.nml' file."
  if (MOM_input_exists) error stop "Remove existing 'MOM_input' file."

  ! Otherwise, create the (empty) files
  open(newunit=io_unit, file='input.nml', status='replace')
  write(io_unit, '(a)') "&fms2_io_nml /"
  close(io_unit)

  open(newunit=io_unit, file='MOM_input', status='replace')
  close(io_unit)
endif

call MOM_infra_init(comm)

! Run tests
call run_file_parser_tests

! Cleanup
call MOM_infra_end

if (rank == root) then
  open(newunit=io_unit, file='MOM_input')
  close(io_unit, status='delete')

  open(newunit=io_unit, file='input.nml')
  close(io_unit, status='delete')
endif

end program test_MOM_file_parser
