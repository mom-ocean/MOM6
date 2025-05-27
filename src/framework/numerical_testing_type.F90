!> A simple type for keeping track of numerical tests
module numerical_testing_type

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public testing
public testing_type_unit_test

!> Class to assist in unit tests, not to be used outside of Recon1d types
type :: testing
  private
  !> True if any fail has been encountered since this instance of "testing" was created
  logical :: state = .false.
  !> Count of tests checked
  integer :: num_tests_checked = 0
  !> Count of tests failed
  integer :: num_tests_failed = 0
  !> If true, be verbose and write results to stdout. Default True.
  logical :: verbose = .true.
  !> Error channel
  integer, public :: stderr = 0
  !> Standard output channel
  integer, public :: stdout = 6
  !> If true, stop instantly
  logical :: stop_instantly = .false.
  !> If true, ignore fails until ignore_fail=.false.
  logical :: ignore_fail = .false.
  !> Record instances that fail
  integer :: ifailed(100) = 0.
  !> Record label of first instance that failed
  character(len=:), allocatable :: label_first_fail

  contains
    procedure :: test => test           !< Update the testing state
    procedure :: set => set             !< Set attributes
    procedure :: summarize => summarize !< Summarize testing state
    procedure :: real_scalar => real_scalar !< Compare two reals
    procedure :: real_arr => real_arr   !< Compare array of reals
    procedure :: int_arr => int_arr     !< Compare array of integers
end type

contains

!> Update the state with "test"
subroutine test(this, state, label, ignore)
  class(testing),    intent(inout) :: this  !< This testing class
  logical,           intent(in)    :: state !< True to indicate a fail, false otherwise
  character(len=*),  intent(in)    :: label !< Message
  logical, optional, intent(in)    :: ignore !< If present and true, ignore a fail
  ! Local variables
  logical :: ignore_this_fail

  ignore_this_fail = this%ignore_fail
  if (present(ignore)) ignore_this_fail = ignore

  this%num_tests_checked = this%num_tests_checked + 1
  if (state) then
    if (.not. ignore_this_fail) then
      this%state = .true.
      this%num_tests_failed = this%num_tests_failed + 1
      if (this%num_tests_failed<=100) this%ifailed(this%num_tests_failed) = this%num_tests_checked
      if (this%num_tests_failed == 1) this%label_first_fail = label
      write(this%stdout, '(2x,3a)') 'Test "',trim(label),'" FAILED!'
      write(this%stderr, '(2x,3a)') 'Test "',trim(label),'" FAILED!'
    else
      write(this%stdout, '(2x,3a)') 'Test "',trim(label),'" IGNORED!'
      write(this%stderr, '(2x,3a)') 'Test "',trim(label),'" IGNORED!'
    endif
  elseif (this%verbose) then
    write(this%stdout, '(2x,3a)') 'Test "',trim(label),'" passed'
  endif
  if (this%stop_instantly .and. this%state .and. .not. ignore_this_fail) stop 1
end subroutine test

!> Set attributes
subroutine set(this, verbose, stdout, stderr, stop_instantly, ignore_fail)
  class(testing), intent(inout) :: this  !< This testing class
  logical, optional, intent(in) :: verbose !< True or false setting to assign to verbosity
  integer, optional, intent(in) :: stdout !< The stdout channel to use
  integer, optional, intent(in) :: stderr !< The stderr channel to use
  logical, optional, intent(in) :: stop_instantly !< If true, stop immediately on error detection
  logical, optional, intent(in) :: ignore_fail !< If true, ignore fails until this option is set false

  if (present(verbose)) then
    this%verbose = verbose
  endif
  if (present(stdout)) then
    this%stdout = stdout
  endif
  if (present(stderr)) then
    this%stderr = stderr
  endif
  if (present(stop_instantly)) then
    this%stop_instantly = stop_instantly
  endif
  if (present(ignore_fail)) then
    this%ignore_fail = ignore_fail
  endif
end subroutine set

!> Summarize results
logical function summarize(this, label)
  class(testing),  intent(inout) :: this  !< This testing class
  character(len=*),   intent(in) :: label !< Message
  integer :: i

  if (this%state) then
    write(this%stdout,'(a," : ",a,", ",i4," failed of ",i4," tested")') &
         'FAIL', trim(label), this%num_tests_failed, this%num_tests_checked
    write(this%stdout,'(a,100i4)') 'Failed tests:',(this%ifailed(i),i=1,min(100,this%num_tests_failed))
    write(this%stdout,'(a,a)') 'First failed test: ',trim(this%label_first_fail)
    write(this%stderr,'(a,100i4)') 'Failed tests:',(this%ifailed(i),i=1,min(100,this%num_tests_failed))
    write(this%stderr,'(a,a)') 'First failed test: ',trim(this%label_first_fail)
    write(this%stderr,'(a," : ",a)') trim(label),'FAILED'
  else
    write(this%stdout,'(a," : ",a,", all ",i4," tests passed")') &
         'Pass', trim(label), this%num_tests_checked
  endif
  summarize = this%state
end function summarize

!> Compare u_test to u_true, report, and return true if a difference larger than tol is measured
!!
!! If in verbose mode, display results to stdout
!! If a difference is measured, display results to stdout and stderr
subroutine real_scalar(this, u_test, u_true, label, tol, robits, ignore)
  class(testing),  intent(inout) :: this   !< This testing class
  real,               intent(in) :: u_test !< Value to test [A]
  real,               intent(in) :: u_true !< Value to test against (correct answer) [A]
  character(len=*),   intent(in) :: label  !< Message
  real,     optional, intent(in) :: tol    !< The tolerance for differences between u and u_true [A]
  integer,  optional, intent(in) :: robits !< Number of bits of round-off to allow
  logical,  optional, intent(in) :: ignore !< If present and true, ignore a fail
  ! Local variables
  logical :: this_test, ignore_this_fail
  real :: tolerance, err ! Tolerance and error [A]

  tolerance = 0.0
  if (present(tol)) tolerance = tol
  ignore_this_fail = this%ignore_fail
  if (present(ignore)) ignore_this_fail = ignore
  this_test = .false.

  ! Scan for any mismatch between u_test and u_true
  if (present(robits)) tolerance = abs(u_true) * float(robits) * epsilon(err)
  if (abs(u_test - u_true) > tolerance) this_test = .true.

  if (this_test) then
    if (ignore_this_fail) then
      if (this%verbose) then
        write(this%stdout,'(3(a,1p1e24.16,1x),2a)') "Calculated value =",u_test,"Correct value =",u_true, &
               "err =",u_test - u_true, label, " <--- IGNORING"
        write(this%stderr,'(3(a,1p1e24.16,1x),2a)') "Calculated value =",u_test,"Correct value =",u_true, &
               "err =",u_test - u_true, label, " <--- IGNORING"
      endif
      this_test = .false.
    else
      write(this%stdout,'(3(a,1p1e24.16,1x),2a)') "Calculated value =",u_test,"Correct value =",u_true, &
               "err =",u_test - u_true, label, " <--- WRONG"
      write(this%stderr,'(3(a,1p1e24.16,1x),2a)') "Calculated value =",u_test,"Correct value =",u_true, &
               "err =",u_test - u_true, label, " <--- WRONG"
    endif
  elseif (this%verbose) then
    write(this%stdout,'(2(a,1p1e24.16,1x),a)') "Calculated value =",u_test,"Correct value =",u_true,label
  endif

  call this%test( this_test, label, ignore=ignore_this_fail ) ! Updates state and counters in this
end subroutine real_scalar

!> Compare u_test to u_true, report, and return true if a difference larger than tol is measured
!!
!! If in verbose mode, display results to stdout
!! If a difference is measured, display results to stdout and stderr
subroutine real_arr(this, n, u_test, u_true, label, tol, robits, ignore)
  class(testing),  intent(inout) :: this   !< This testing class
  integer,            intent(in) :: n      !< Number of cells in u
  real, dimension(n), intent(in) :: u_test !< Values to test [A]
  real, dimension(n), intent(in) :: u_true !< Values to test against (correct answer) [A]
  character(len=*),   intent(in) :: label  !< Message
  real,     optional, intent(in) :: tol    !< The tolerance for differences between u and u_true [A]
  integer,  optional, intent(in) :: robits !< Number of bits of round-off to allow
  logical,  optional, intent(in) :: ignore !< If present and true, ignore a fail
  ! Local variables
  integer :: k
  logical :: this_test, ignore_this_fail
  real :: tolerance, err ! Tolerance and error [A]

  tolerance = 0.0
  if (present(tol)) tolerance = tol
  ignore_this_fail = this%ignore_fail
  if (present(ignore)) ignore_this_fail = ignore
  this_test = .false.

  ! Scan for any mismatch between u_test and u_true
  do k = 1, n
    if (present(robits)) tolerance = abs(u_true(k)) * float(robits) * epsilon(err)
    if (abs(u_test(k) - u_true(k)) > tolerance) this_test = .true.
  enddo

  ! If either being verbose, or an error was measured then display results
  if (this_test .or. this%verbose) then
    write(this%stdout,'(a4,2a24,1x,a)') 'k','Calculated value','Correct value',label
    if (this_test) write(this%stderr,'(a4,2a24,1x,a)') 'k','Calculated value','Correct value',label
    do k = 1, n
      if (present(robits)) tolerance = abs(u_true(k)) * float(robits) * epsilon(err)
      err = u_test(k) - u_true(k)
      if ( ( abs(err) > tolerance .and. ignore_this_fail ) .or. &
           ( abs(err) > 0. .and. abs(err) <= tolerance ) ) then
        write(this%stdout,'(i4,1p2e24.16,a,1pe24.16,a)') k, u_test(k), u_true(k), &
                         ' err=', err, ' <--- IGNORING'
      elseif (abs(err) > tolerance) then
        write(this%stdout,'(i4,1p2e24.16,a,1pe24.16,a)') k, u_test(k), u_true(k), &
                         ' err=', err, ' <--- WRONG'
        write(this%stderr,'(i4,1p2e24.16,a,1pe24.16,a)') k, u_test(k), u_true(k), &
                         ' err=', err, ' <--- WRONG'
      else
        write(this%stdout,'(i4,1p2e24.16)') k, u_test(k), u_true(k)
      endif
    enddo
  endif

  call this%test( this_test, label, ignore=ignore_this_fail ) ! Updates state and counters in this
end subroutine real_arr

!> Compare i_test to i_true and report and return true if a difference is found
!!
!! If in verbose mode, display results to stdout
!! If a difference is measured, display results to stdout and stderr
subroutine int_arr(this, n, i_test, i_true, label, ignore)
  class(testing),     intent(inout) :: this   !< This testing class
  integer,               intent(in) :: n      !< Number of cells in u
  integer, dimension(n), intent(in) :: i_test !< Values to test [A]
  integer, dimension(n), intent(in) :: i_true !< Values to test against (correct answer) [A]
  character(len=*),      intent(in) :: label  !< Message
  logical,  optional,    intent(in) :: ignore !< If present and true, ignore a fail
  ! Local variables
  integer :: k
  logical :: this_test, ignore_this_fail

  ignore_this_fail = this%ignore_fail
  if (present(ignore)) ignore_this_fail = ignore
  this_test = .false.

  ! Scan for any mismatch between u_test and u_true
  do k = 1, n
    if (i_test(k) /= i_true(k)) this_test = .true.
  enddo

  if (this%verbose) then
    write(this%stdout,'(a14," : calculated =",30i3)') label, i_test
    write(this%stdout,'(14x,"      correct =",30i3)') i_true
    if (this_test) then
      if (ignore_this_fail) then
        write(this%stdout,'(3x,a,8x,"error =",30i3)') 'IGNORE --->', i_test(:) - i_true(:)
      else
        write(this%stdout,'(3x,a,8x,"error =",30i3)') ' FAIL  --->', i_test(:) - i_true(:)
      endif
    endif
  endif

  if (ignore_this_fail) this_test = .false.

  if (this_test) then
    write(this%stderr,'(a14," : calculated =",30i3)') label, i_test
    write(this%stderr,'(14x,"      correct =",30i3)') i_true
    write(this%stderr,'("   FAIL --->        error =",30i3)') i_test(:) - i_true(:)
  endif

  call this%test( this_test, label ) ! Updates state and counters in this
end subroutine int_arr

!> Tests the testing type itself
logical function testing_type_unit_test(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables
  type(testing) :: test ! The instance to be tested
  logical :: tmpflag ! Temporary for return flags

  testing_type_unit_test = .false. ! Assume all is well at the outset
  if (verbose) write(test%stdout,*) "  ===== testing_type: testing_type_unit_test ============"

  call test%set( verbose=verbose ) ! Sets the verbosity flag in test
  call test%set( stderr=0 ) ! Sets stderr
  call test%set( stdout=6 ) ! Sets stdout
  call test%set( stop_instantly=.false. ) ! Sets stop_instantly
  call test%set( ignore_fail=.false. ) ! Sets ignore_fail

  call test%test( .false., "This should pass" )
  if (verbose .and. .not. test%state) then
    write(test%stdout,*) "   => test(F) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  call test%test( .true., "This should fail but be ignored", ignore=.true. )
  if (verbose .and. .not. test%state) then
    write(test%stdout,*) "   => test(T,ignore) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  call test%real_scalar(1., 1., "s == s should pass", robits=0, tol=0.)
  if (verbose .and. .not. test%state) then
    write(test%stdout,*) "   => real(s,s) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  call test%real_scalar(1., 2., "s != t but ignored", ignore=.true.)
  if (verbose .and. .not. test%state) then
    write(test%stdout,*) "   => real(s,t,ignore) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  call test%real_arr(2, (/1.,2./), (/1.,2./), "a == a should pass", robits=0, tol=0.)
  if (verbose .and. .not. test%state) then
    write(test%stdout,*) "   => real(a,a) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  call test%real_arr(2, (/1.,2./), (/3.,4./), "a != b but ignored", ignore=.true.)
  if (verbose .and. .not. test%state) then
    write(test%stdout,*) "   => real(a,b,ignore) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  call test%int_arr(2, (/1,2/), (/1,2/), "i == i should pass")
  if (verbose .and. .not. test%state) then
    write(test%stdout,*) "   => int(a,a) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  call test%int_arr(2, (/1,2/), (/3,4/), "i != j but ignored", ignore=.true.)
  if (verbose .and. .not. test%state) then
    write(test%stdout,*) "   => int(a,b,ignore) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  tmpflag = test%summarize("This summary is for a passing state")
  if (verbose .and. .not. tmpflag) then
    write(test%stdout,*) "   => summarize(F) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  ! This following all fail
  test%state = .false. ! reset
  call test%test( .true., "This should fail" )
  if (verbose .and. test%state) then
    write(test%stdout,*) "   => test(T) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  test%state = .false. ! reset
  call test%real_scalar(1., 2., "s != t should fail")
  if (verbose .and. test%state) then
    write(test%stdout,*) "   => real(s,t) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  test%state = .false. ! reset
  call test%real_arr(2, (/1.,2./), (/3.,4./), "a != b and should fail")
  if (verbose .and. test%state) then
    write(test%stdout,*) "   => real(a,b) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  test%state = .false. ! reset
  call test%int_arr(2, (/1,2/), (/3,4/), "i != j and should fail")
  if (verbose .and. test%state) then
    write(test%stdout,*) "   => int(a,b) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  tmpflag = test%summarize("This summary should have 3 fails")
  if (verbose .and. tmpflag) then
    write(test%stdout,*) "   => summarize(T) passed"
  else; testing_type_unit_test = testing_type_unit_test .or. .true. ; endif

  if (verbose .and. .not. testing_type_unit_test) write(test%stdout,*) "testing_type_unit_test passed"

end function testing_type_unit_test

!> \namespace numerical_testing_type
!!
end module numerical_testing_type
