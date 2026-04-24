! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

program test_MOM_array_transform

use MOM_array_transform, only : symmetric_sum_unit_tests

if ( symmetric_sum_unit_tests(.true.) ) stop 1

end program test_MOM_array_transform
