program test_MOM_EOS

use MOM_EOS, only           : EOS_unit_tests
use MOM_error_handler, only : set_skip_mpi

call set_skip_mpi(.true.) ! This unit tests is not expecting MPI to be used

if ( EOS_unit_tests(.true.) ) stop 1

end program test_MOM_EOS
