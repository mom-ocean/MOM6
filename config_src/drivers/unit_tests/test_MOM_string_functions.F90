program test_MOM_string_functions

use MOM_string_functions, only : string_functions_unit_tests
use MOM_error_handler, only : set_skip_mpi

call set_skip_mpi(.true.) ! This unit tests is not expecting MPI to be used

if ( string_functions_unit_tests(.true.) ) stop 1

end program test_MOM_string_functions
