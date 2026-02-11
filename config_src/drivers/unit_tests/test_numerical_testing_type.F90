! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

program test_numerical_testing_type

use numerical_testing_type, only : numerical_testing_type_unit_tests

if (numerical_testing_type_unit_tests(.true.)) stop 1

end program test_numerical_testing_type
