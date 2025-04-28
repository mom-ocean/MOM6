program test_MOM_remapping

use MOM_remapping, only : remapping_unit_tests

integer :: n !< Number of arguments, or tests
character(len=12) :: cmd_ln_arg !< Command line argument (if any)

n = command_argument_count()

if (n==1) then
  call get_command_argument(1, cmd_ln_arg)
  read(cmd_ln_arg,*) n
else
  n = 3000 ! Fallback value if no argument provided
endif

if (remapping_unit_tests(.true., num_comp_samp=n)) stop 1

end program test_MOM_remapping
