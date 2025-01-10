program test_MOM_remapping

use MOM_remapping, only : remapping_unit_tests

if (remapping_unit_tests(.true.)) stop 1

end program test_MOM_remapping
