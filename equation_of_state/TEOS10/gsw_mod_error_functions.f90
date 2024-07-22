!==========================================================================
module gsw_mod_error_functions
!==========================================================================

use gsw_mod_kinds

implicit none

logical, public :: gsw_error_check = .true.
logical, public :: gsw_abort_on_error = .true.

real (r8), parameter, public :: gsw_error_limit = 1e10_r8

integer, parameter, private :: nfuncs = 37
integer, parameter, private :: maxlen = 40
character (len=maxlen), dimension(nfuncs), private :: func_list

data func_list / &
		"gsw_ct_from_enthalpy_exact", &
		"gsw_ct_from_enthalpy", &
		"gsw_ct_from_rho", &
		"gsw_deltasa_atlas", &
		"gsw_deltasa_from_sp", &
		"gsw_fdelta", &
		"gsw_frazil_properties", &
		"gsw_frazil_properties_potential", &
		"gsw_frazil_properties_potential_poly", &
		"gsw_geo_strf_dyn_height", &
		"gsw_geo_strf_dyn_height_pc", &
		"gsw_ice_fraction_to_freeze_seawater", &
		"gsw_ipv_vs_fnsquared_ratio", &
		"gsw_melting_ice_into_seawater", &
		"gsw_melting_ice_sa_ct_ratio", &
		"gsw_melting_ice_sa_ct_ratio_poly", &
		"gsw_melting_seaice_into_seawater", &
		"gsw_melting_seaice_sa_ct_ratio", &
		"gsw_melting_seaice_sa_ct_ratio_poly", &
		"gsw_nsquared", &
		"gsw_pressure_freezing_ct", &
		"gsw_rr68_interp_sa_ct", &
		"gsw_saar", &
		"gsw_sa_freezing_from_ct", &
		"gsw_sa_freezing_from_ct_poly", &
		"gsw_sa_freezing_from_t", &
		"gsw_sa_freezing_from_t_poly", &
		"gsw_sa_from_rho", &
		"gsw_sa_from_sp", &
		"gsw_sa_from_sstar", &
		"gsw_seaice_fraction_to_freeze_seawater", &
		"gsw_sp_from_c", &
		"gsw_sp_from_sa", &
		"gsw_sp_from_sstar", &
		"gsw_sstar_from_sa", &
		"gsw_sstar_from_sp", &
		"gsw_turner_rsubrho" /

public :: gsw_error_code
public :: gsw_error_handler

private :: gsw_error_fnum

contains

    elemental function gsw_error_code (err_num, func_name, error_code)

    ! Constructs an error code of the form 9.nabcxyz000000d15
    !
    ! where n   = current error level (1-4)
    !       abc = error code for level #1
    !       xyz = error code for level #2
    !       ...
    ! and level error codes comprise ...
    !       a  = error number for level #1 (0-9)
    !       bc = function number for level #10

    implicit none

    integer, intent(in) :: err_num
    character (*), intent(in) :: func_name
    real (r8), intent(in), optional :: error_code

    integer :: ival, k
    real (r8) :: gsw_error_code, base_code, mult

    if (present(error_code)) then
        k = int(error_code/1.0e14_r8) - 90
	base_code = error_code + 1.0e14_r8
	mult = 10.0_r8**(11-k*3)
    else
        base_code = 9.1e15_r8
	mult = 1.0e11_r8
    end if

    ival = err_num*100 + gsw_error_fnum(func_name)
    gsw_error_code = base_code + ival*mult

    end function gsw_error_code

    !==========================================================================

    elemental function gsw_error_fnum (func_name)

    implicit none

    character (*), intent(in) :: func_name

    integer :: gsw_error_fnum

    integer :: i
    character (len=maxlen) :: fname

    fname = func_name
    do i = 1, nfuncs
        if (fname == func_list(i)) goto 100
    end do
    gsw_error_fnum = 99
    return

100 gsw_error_fnum = i
    return

    end function gsw_error_fnum

    !==========================================================================

    subroutine gsw_error_handler (error_code)

    implicit none

    real (r8), intent(in) :: error_code

    integer (selected_int_kind(14)) :: base_code
    integer :: func_num, ival, i, k
    real (r8) :: gsw_error_code, mult

    character (len=maxlen) :: func_name

    print '(/"Trace for error code: ", es20.13/)', error_code

    base_code = error_code - 9.0e15_r8
    k = int(base_code/1.0e14_r8)
    base_code = base_code/(10**(14-k*3))

    do i = 1, k
    	ival =  mod(base_code,1000)
	func_num = mod(ival,100)
	if (func_num .le. nfuncs) then
	    func_name = func_list(func_num)
	else
	    func_name = "unknown"
	end if
        print '("  Code: ",i1," in function: ",a)', ival/100, func_name
	base_code = base_code/1000
    end do

    if (gsw_abort_on_error) stop

    end subroutine gsw_error_handler

end module gsw_mod_error_functions
