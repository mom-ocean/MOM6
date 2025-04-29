!> A simple type for keeping track of numerical tests
module numerical_testing_type

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public testing
public numerical_testing_type_unit_tests

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
logical function numerical_testing_type_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables
  type(testing) :: tester ! An instance to record tests
  type(testing) :: test ! The instance used for testing (is mutable)
  logical :: tmpflag ! Temporary for return flags

  numerical_testing_type_unit_tests = .false. ! Assume all is well at the outset
  if (verbose) write(test%stdout,*) "  ===== testing_type: numerical_testing_type_unit_tests ====="
  call tester%set( verbose=verbose ) ! Sets the verbosity flag in tester

  call test%set( verbose=verbose ) ! Sets the verbosity flag in test
  call test%set( stderr=6 ) ! Sets stderr (redirect errors for "test" since they are not real)
  call test%set( stdout=6 ) ! Sets stdout
  call test%set( stop_instantly=.false. ) ! Sets stop_instantly
  call test%set( ignore_fail=.false. ) ! Sets ignore_fail

  ! Check that %summary() reports nothing when %state is unset
  ! (note this has to be confirmed visually since everything is in stdout)
  tmpflag = test%summarize("Summary is for a passing state")
  call tester%test(tmpflag, "test%summarize() with no fails")

  ! Check that %test(.false.,...) leaves %state unchanged
  call test%test( .false., "test(F) should pass" )
  call tester%test(test%state, "test%test(F)")

  ! Check that %test(.true.,...,ignore=.true.) leaves %state unchanged
  call test%test( .true., "test(T) should fail but be ignored", ignore=.true. )
  call tester%test(test%state, "test%test(T,ignore)")

  ! Check that %test(.true.,...) sets %state
  call test%test( .true., "test(T) should fail" )
  call tester%test(.not. test%state, "test%test(T,ignore)")
  test%state = .false. ! reset

  ! Check that %real_scalar(a,a,...) leaves %state unchanged
  call test%real_scalar(1., 1., "real_scalar(s,s) should pass", robits=0, tol=0.)
  call tester%test(test%state, "test%real_scalar(s,s)")

  ! Check that %real_scalar(a,b,...,ignore=.true.) leaves %state unchanged
  call test%real_scalar(1., 2., "real_scalar(s,t) should fail but be ignored", ignore=.true.)
  call tester%test(test%state, "test%real_scalar(s,t,ignore)")

  ! Check that %real_scalar(a,a,...) sets %state
  call test%real_scalar(1., 2., "s != t should fail")
  call tester%test(.not. test%state, "test%real_scalar(s,t)")
  test%state = .false. ! reset

  ! Check that %real_arr(a,a,...) leaves %state unchanged
  call test%real_arr(2, (/1.,2./), (/1.,2./), "real_arr(a,a) should pass", robits=0, tol=0.)
  call tester%test(test%state, "test%real_arr(a,a)")

  ! Check that %real_arr(a,b,...,ignore=.true.) leaves %state unchanged
  call test%real_arr(2, (/1.,2./), (/3.,4./), "real_arr(a,b) should fail but be ignored", ignore=.true.)
  call tester%test(test%state, "test%real_arr(a,b,ignore)")

  ! Check that %real_arr(a,b,...) sets %state
  call test%real_arr(2, (/1.,2./), (/3.,4./), "real(a,b) should fail")
  call tester%test(.not. test%state, "test%real_arr(a,b)")
  test%state = .false. ! reset

  ! Check that %int_arr(a,a,...) leaves %state unchanged
  call test%int_arr(2, (/1,2/), (/1,2/), "int_arr(i,i) should pass")
  call tester%test(test%state, "test%int_arr(i,i)")

  ! Check that %int_arr(a,b,...,ignore=.true.) leaves %state unchanged
  call test%int_arr(2, (/1,2/), (/3,4/), "int_arr(i,j) should fail but be ignored", ignore=.true.)
  call tester%test(test%state, "test%int_arr(i,j,ignore)")

  ! Check that %int_arr(a,b,...) sets %state
  call test%int_arr(2, (/1,2/), (/3,4/), "int(arr(i,j) should fail")
  call tester%test(.not. test%state, "test%int_arr(i,j)")
  test%state = .false. ! reset

  ! Check that %summary() reports nothing when %state is set
  ! (note this has to be confirmed visually since everything is in stdout)
  test%state = .true. ! reset to fail for testing %summary()
  tmpflag = test%summarize("This summary should report 4 fails")
  call tester%test(.not. tmpflag, "test%summarize() with fails")

  numerical_testing_type_unit_tests = tester%summarize("numerical_testing_type_unit_tests")

end function numerical_testing_type_unit_tests

!> \namespace numerical_testing_type
!!
!! numerical_testing_type is a helper class to facilitate implementing
!! tests of a numerical nature.
!! The class helps hide the logic and code associated with handling the
!! results of a test, essentially reducing the multiple lines of `if
!! ... then ... print ... else ... error_mesg ...` into one line.
!!
!! The class is light weight, meaning is does not depend on anything else,
!! allowing to be particularly useful in unit tests and small drivers.
!! However, this means it is up to the user to do something with the results,
!! e.g. `call MOM_error()` appropriately.
!!
!! Each test, e.g. real_scalar(), is expected to pass.
!! If a fail is encountered, it is immediately reported to stderr and stdour,
!! recorded internally, but does not terminate execuation unless
!! `set(stop_instantly=.true.)` was called previously.
!! Most tests take the form of `f(a,b)` where `a` should equal `b`.
!! Only test() takes a single input (boolean) which is expected to
!! be false for the test to pass.
!!
!! summarize() is used to "finalize" the tests.
!! It prints a summary of how many and which tests faield, and returns a logical
!! that is set to .true. if any test failed.
!!
!! Usage by example:
!! \verbatim
!! use numerical_testing_type, only : testing
!! ...
!!
!! !> Runs my unit_tests. Returns .true. if a test fails, .false. otherwise
!! logical function my_unit_tests(verbose)
!!   logical, intent(in) :: verbose !< If true, write results to stdout
!!   ...
!!   type(testing) :: test ! An instance of the numerical_testing_type
!!   ...
!!   call test%set( verbose=.true. ) ! Show intermediate results rather than just the fails
!!   ...
!!
!!   call test%test(flag, 'Flag is not set')               ! Check flag=.false.
!!   call test%real_scalar(a, 1., 'u = 1')                 ! Check a=1
!!   call test%real_arr(3, u, (/1.,2.,3./), 'u = [1,2,3]') ! Check u(:)=[1,2,3]
!!   call test%int_arr(2, iv, (/1,2/), 'iv = [1,2]')       ! Check that iv(:)=[1,2]
!!
!!   my_unit_tests = test%summarize('my_unit_tests') ! Return true if a fail occurs
!! end function my_unit_tests(verbose)
!! \endverbatim

end module numerical_testing_type
