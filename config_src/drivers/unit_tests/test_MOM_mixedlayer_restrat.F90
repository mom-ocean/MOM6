! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

program test_MOM_mixedlayer_restrat

use MOM_mixed_layer_restrat, only : mixedlayer_restrat_unit_tests
use MOM_error_handler, only : set_skip_mpi

call set_skip_mpi(.true.) ! This unit tests is not expecting MPI to be used

if ( mixedlayer_restrat_unit_tests(.true.) ) stop 1

end program test_MOM_mixedlayer_restrat
