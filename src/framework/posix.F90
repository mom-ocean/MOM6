!> Interface to the libc POSIX API
#include "posix.h"

module posix

use, intrinsic :: iso_c_binding, only : c_char
use, intrinsic :: iso_c_binding, only : c_int
use, intrinsic :: iso_c_binding, only : c_long
use, intrinsic :: iso_c_binding, only : c_null_char
use, intrinsic :: iso_c_binding, only : c_funptr
use, intrinsic :: iso_c_binding, only : c_funloc
use, intrinsic :: iso_c_binding, only : c_f_procpointer

implicit none

!> Container for file metadata from stat
!!
!! NOTE: This is currently just a placeholder containing fields, such as size,
!! uid, mode, etc.  A readable Fortran type may be used in the future.
type, bind(c) :: stat_buf
  private
  character(kind=c_char) :: state(SIZEOF_STAT_BUF)
    !< Byte array containing file metadata
end type stat_buf

!> Container for the jump point buffer created by setjmp().
!!
!! The buffer typically contains the current register values, stack pointers,
!! and any information required to restore the process state.
type, bind(c) :: jmp_buf
  private
  character(kind=c_char) :: state(SIZEOF_JMP_BUF)
    !< Unstructured array of bytes used to store the process state
end type jmp_buf

!> Container for the jump point buffer (with signals) created by sigsetjmp()
!!
!! In addition to the content stored by `jmp_buf`, it also stores signal state.
type, bind(c) :: sigjmp_buf
  private
  character(kind=c_char) :: state(SIZEOF_SIGJMP_BUF)
    !< Unstructured array of bytes used to store the process state
end type sigjmp_buf

! POSIX signals
integer, parameter :: SIGUSR1 = POSIX_SIGUSR1
  !< Signal number for SIGUSR1 (user-defined signal 1)

interface
  !> C interface to POSIX chmod()
  !! Users should use the Fortran-defined chmod() function.
  function chmod_posix(path, mode) result(rc) bind(c, name="chmod")
    ! #include <sys/stat.h>
    ! int chmod(const char *path, mode_t mode);
    import :: c_char, c_int

    character(kind=c_char), dimension(*), intent(in) :: path
      !< Zero-delimited file path
    integer(kind=c_int), value, intent(in) :: mode
      !< File permission to be assigned to file.
    integer(kind=c_int) :: rc
      !< Function return code
  end function chmod_posix

  !> C interface to POSIX mkdir()
  !! Users should use the Fortran-defined mkdir() function.
  function mkdir_posix(path, mode) result(rc) bind(c, name="mkdir")
    ! #include <sys/stat.h>
    ! int mkdir(const char *path, mode_t mode);
    import :: c_char, c_int

    character(kind=c_char), dimension(*), intent(in) :: path
      !< Zero-delimited file path
    integer(kind=c_int), value, intent(in) :: mode
      !< File permission to be assigned to file.
    integer(kind=c_int) :: rc
      !< Function return code
  end function mkdir_posix

  !> C interface to POSIX stat()
  !! Users should use the Fortran-defined stat() function.
  function stat_posix(path, buf) result(rc) bind(c, name="stat")
    import :: c_char, stat_buf, c_int

    character(kind=c_char), dimension(*), intent(in) :: path
      !< Pathname of a POSIX file
    type(stat_buf), intent(inout) :: buf
      !< Information describing the file if it exists
    integer(kind=c_int) :: rc
      !< Function return code
  end function

  !> C interface to POSIX signal()
  !! Users should use the Fortran-defined signal() function.
  function signal_posix(sig, func) result(handle) bind(c, name="signal")
    ! #include <signal.h>
    ! void (*signal(int sig, void (*func)(int)))(int);
    import :: c_int, c_funptr

    integer(kind=c_int), value, intent(in) :: sig
      !< Signal to be configured
    type(c_funptr), value, intent(in) :: func
      !< Function handle to be called when `sig` is raised
    type(c_funptr) :: handle
      !< Prior handle for sig to be replaced by `func`
  end function signal_posix

  !> C interface to POSIX kill()
  !! Users should use the Fortran-defined kill() function.
  function kill_posix(pid, sig) result(rc) bind(c, name="kill")
    ! #include <signal.h>
    ! int kill(pid_t pid, int sig);
    import :: c_int

    integer(kind=c_int), value, intent(in) :: pid
      !< Process ID which is to receive the signal
    integer(kind=c_int), value, intent(in) :: sig
      !< Signal to be sent to the process
    integer(kind=c_int) :: rc
      !< Function return code
  end function kill_posix

  !> C interface to POSIX getpid()
  !! Users should use the Fortran-defined getpid() function.
  function getpid_posix() result(pid) bind(c, name="getpid")
    ! #include <unistd.h>
    ! pid_t getpid(void);
    import :: c_long

    integer(kind=c_long) :: pid
      !< Process ID of the current process.
  end function getpid_posix

  !> C interface to POSIX getppid()
  !! Users should use the Fortran-defined getppid() function.
  function getppid_posix() result(pid) bind(c, name="getppid")
    ! #include <unistd.h>
    ! pid_t getppid(void);
    import :: c_long

    integer(kind=c_long) :: pid
      !< Process ID of the parent process to the current process.
  end function getppid_posix

  !> C interface to POSIX sleep()
  !! Users should use the Fortran-defined sleep() function.
  function sleep_posix(seconds) result(rc) bind(c, name="sleep")
    ! #include <unistd.h>
    ! unsigned int sleep(unsigned int seconds);
    import :: c_int

    integer(kind=c_int), value, intent(in) :: seconds
      !< Number of real-time seconds which the thread should sleep
    integer(kind=c_int) :: rc
      !< Function return code
  end function

  ! NOTE: The C setjmp and sigsetjmp functions *must* be called explicitly by
  ! the Fortran code, rather than through a wrapper Fortran function.
  !
  ! Otherwise, setjmp() will capture the stack inside the wrapper, rather than
  ! the point where setjmp() is called.
  !
  ! Hence, we remove the `_posix` suffix and call these explicitly.
  ! (The integer kind <-> c_int conversion will need to be addressed.)

  ! NOTE: POSIX explicitly says setjmp/sigsetjmp may be either a function or a
  !   macro, and thus bind() may point to a nonexistent function.
  ! e.g. sigsetjmp is a macro to __sigsetjmp in glibc, so we use a macro.

  !> Save the current program execution state to `env`.
  !!
  !! This function creates a snapshot of the process state to `env`, which can
  !! be restored by calling `longjmp`.  When `setjmp` is called, the function
  !! returns 0.  When `longjmp` is later called, the program is restored to the
  !! point where `setjmp` was called, except it now returns a value (rc) as
  !! specified by `longjmp`.
  function setjmp(env) result(rc) bind(c, name=SETJMP_NAME)
    ! #include <setjmp.h>
    ! int setjmp(jmp_buf env);
    import :: jmp_buf, c_int

    type(jmp_buf), intent(in) :: env
      !< Current process state
    integer(kind=c_int) :: rc
      !< Function return code; set to 0 if setjmp() was called, otherwise
      !! specified by the corresponding longjmp() call.
  end function setjmp

  !> Save the current execution and ,optionally, the signal state to `env`.
  !!
  !! This function creates a snapshot of the process state to `env`, which can
  !! be restored by calling `longjmp`.  When `setjmp` is called, the function
  !! returns 0.  When `longjmp` is later called, the program is restored to the
  !! point where `setjmp` was called, except it now returns a value (rc) as
  !! specified by `longjmp`.
  !!
  !! If `savesigs` is set to a nonzero value, then the signal state is included
  !! in the program state.
  function sigsetjmp(env, savesigs) result(rc) bind(c, name=SIGSETJMP_NAME)
    ! #include <setjmp.h>
    ! int sigsetjmp(jmp_buf env, int savesigs);
    import :: sigjmp_buf, c_int

    type(sigjmp_buf), intent(in) :: env
      !< Current process state
    integer(kind=c_int), value, intent(in) :: savesigs
      !< Flag to enable signal state when set to a nonzero value
    integer(kind=c_int) :: rc
      !< Function return code; set to 0 if sigsetjmp() was called, otherwise
      !! specified by the corresponding siglongjmp() call.
  end function sigsetjmp

  !> C interface to POSIX longjmp()
  !! Users should use the Fortran-defined longjmp() function.
  subroutine longjmp_posix(env, val) bind(c, name=LONGJMP_NAME)
    ! #include <setjmp.h>
    ! int longjmp(jmp_buf env, int val);
    import :: jmp_buf, c_int

    type(jmp_buf), intent(in) :: env
      !< Process state to restore
    integer(kind=c_int), value, intent(in) :: val
      !< Return code sent to setjmp()
  end subroutine longjmp_posix

  !> C interface to POSIX siglongjmp()
  !! Users should use the Fortran-defined siglongjmp() function.
  subroutine siglongjmp_posix(env, val) bind(c, name=SIGLONGJMP_NAME)
    ! #include <setjmp.h>
    ! int siglongjmp(jmp_buf env, int val);
    import :: sigjmp_buf, c_int

    type(sigjmp_buf), intent(in) :: env
      !< Process state to restore
    integer(kind=c_int), value, intent(in) :: val
      !< Return code sent to sigsetjmp()
  end subroutine siglongjmp_posix

  ! Note on types:
  ! mode_t:
  !   "According to POSIX, it shall be an integer type."
  ! pid_t:
  !   "According to POSIX, it shall be a signed integer type, and the
  !   implementation shall support one or more programming environments where
  !   the width of pid_t is no greater than the width of the type long.
  ! jmp_buf:
  !   This is a strongly platform-dependent variable, large enough to contain
  !   a complete copy of the process execution state (registers, stack, etc).
  ! sigjmp_buf:
  !   A more comprehensive version of jmp_buf which contains signal state.
end interface

abstract interface
  !> Function interface for signal handlers
  subroutine handler_interface(sig)
    integer, intent(in) :: sig
      !> Input signal to handler
  end subroutine
end interface

contains

!> Change mode of a file
!!
!! This changes the file permission of file `path` to `mode` following POSIX
!! conventions.  If successful, it returns zero.  Otherwise, it returns -1.
function chmod(path, mode) result(rc)
  character(len=*), intent(in) :: path
  integer, intent(in) :: mode
  integer :: rc

  integer(kind=c_int) :: mode_c
  integer(kind=c_int) :: rc_c

  mode_c = int(mode, kind=c_int)
  rc_c = chmod_posix(path//c_null_char, mode_c)
  rc = int(rc_c)
end function chmod

!> Create a file directory
!!
!! This creates a new directory named `path` with permissons set by `mode`.
!! If successful, it returns zero.  Otherwise, it returns -1.
function mkdir(path, mode) result(rc)
  character(len=*), intent(in) :: path
  integer, intent(in) :: mode
  integer :: rc

  integer(kind=c_int) :: mode_c
  integer(kind=c_int) :: rc_c

  mode_c = int(mode, kind=c_int)
  rc_c = mkdir_posix(path//c_null_char, mode_c)
  rc = int(rc_c)
end function mkdir

!> Get file status
!!
!! This obtains information about the named file and writes it to buf.
!! If found, it returns zero.  Otherwise, it returns -1.
function stat(path, buf) result(rc)
  character(len=*), intent(in) :: path
    !< Pathname of file to be inspected
  type(stat_buf), intent(out) :: buf
    !< Buffer containing information about the file if it exists
    ! NOTE: Currently the contents of buf are not readable, but we could move
    ! the contents into a readable Fortran type.
  integer :: rc
    !< Function return code

  integer(kind=c_int) :: rc_c

  rc_c = stat_posix(path//c_null_char, buf)

  rc = int(rc_c)
end function stat

!> Create a signal handler `handle` to be called when `sig` is detected.
!!
!! If successful, the previous handler for `sig` is returned.  Otherwise,
!! SIG_ERR is returned.
function signal(sig, func) result(handle)
  integer, intent(in) :: sig
  procedure(handler_interface) :: func
  procedure(handler_interface), pointer :: handle

  integer(kind=c_int) :: sig_c
  type(c_funptr) :: handle_c

  sig_c = int(sig, kind=c_int)
  handle_c = signal_posix(sig_c, c_funloc(func))
  call c_f_procpointer(handle_c, handle)
end function signal

!> Send signal `sig` to process `pid`.
!!
!! If successful, this function returns 0.  Otherwise, it returns -1.
function kill(pid, sig) result(rc)
  integer, intent(in) :: pid
  integer, intent(in) :: sig
  integer :: rc

  integer(kind=c_int) :: pid_c, sig_c, rc_c

  pid_c = int(pid, kind=c_int)
  sig_c = int(sig, kind=c_int)
  rc_c = kill_posix(pid_c, sig_c)
  rc = int(rc_c)
end function kill

!> Get the ID of the current process.
function getpid() result(pid)
  integer :: pid

  integer(kind=c_long) :: pid_c

  pid_c = getpid_posix()
  pid = int(pid_c)
end function getpid

!> Get the ID of the parent process of the current process.
function getppid() result(pid)
  integer :: pid

  integer(kind=c_long) :: pid_c

  pid_c = getppid_posix()
  pid = int(pid_c)
end function getppid

!> Force the process to a sleep state for `seconds` seconds.
!!
!! The sleep state may be interrupted by a signal.  If it sleeps for the entire
!! duration, then it returns 0.  Otherwise, it returns the number of seconds
!! remaining at the point of interruption.
function sleep(seconds) result(rc)
  ! NOTE: This function may replace an existing compiler `sleep()` extension.
  integer, intent(in) :: seconds
  integer :: rc

  integer(kind=c_int) :: seconds_c
  integer(kind=c_int) :: rc_c

  seconds_c = int(seconds, kind=c_int)
  rc_c = sleep_posix(seconds_c)
  rc = int(rc_c)
end function sleep

!> Restore program to state saved by `env`, and return the value `val`.
!!
!! This "nonlocal goto" alters program execution to the state stored in `env`
!! produced by a prior execution of `setjmp`.  Program execution is moved
!! back to this `setjmp`, except the function will now return `val`.
subroutine longjmp(env, val)
  type(jmp_buf), intent(in) :: env
  integer, intent(in) :: val

  integer(kind=c_int) :: val_c

  val_c = int(val, kind=c_int)
  call longjmp_posix(env, val_c)
end subroutine longjmp

!> Restore program to state saved by `env`, and return the value `val`.
!!
!! This "nonlocal goto" alters program execution to the state stored in `env`
!! produced by a prior execution of `setjmp`.  Program execution is moved back
!! to this `setjmp`, except the function will now return `val`.
!!
!! `siglongjmp` behaves in the same manner as `longjmp`, but also provides
!! predictable handling of the signal state.
subroutine siglongjmp(env, val)
  type(sigjmp_buf), intent(in) :: env
  integer, intent(in) :: val

  integer(kind=c_int) :: val_c

  val_c = int(val, kind=c_int)
  call siglongjmp_posix(env, val_c)
end subroutine siglongjmp


! Symbols in <setjmp.h> may be platform-dependent and may not exist if defined
! as a macro.  The following functions permit compilation when they are
! unavailable, and report a runtime error if used in the program.

!> Placeholder function for a missing or unconfigured setjmp
function setjmp_missing(env) result(rc) bind(c)
  type(jmp_buf), intent(in) :: env
    !< Current process state (unused)
  integer(kind=c_int) :: rc
    !< Function return code (unused)

  print '(a)', 'ERROR: setjmp() is not implemented in this build.'
  print '(a)', 'Recompile with autoconf or -DSETJMP_NAME=\"<symbol name>\".'
  error stop

  ! NOTE: compilers may expect a return value, even if it is unreachable
  read env%state
  rc = -1
end function setjmp_missing

!> Placeholder function for a missing or unconfigured longjmp
subroutine longjmp_missing(env, val) bind(c)
  type(jmp_buf), intent(in) :: env
    !< Current process state (unused)
  integer(kind=c_int), value, intent(in) :: val
    !< Enable signal state flag (unused)

  print '(a)', 'ERROR: longjmp() is not implemented in this build.'
  print '(a)', 'Recompile with autoconf or -DLONGJMP_NAME=\"<symbol name>\".'
  error stop

  read env%state
  read char(val)
end subroutine longjmp_missing

!> Placeholder function for a missing or unconfigured sigsetjmp
function sigsetjmp_missing(env, savesigs) result(rc) bind(c)
  type(sigjmp_buf), intent(in) :: env
    !< Current process state (unused)
  integer(kind=c_int), value, intent(in) :: savesigs
    !< Enable signal state flag (unused)
  integer(kind=c_int) :: rc
    !< Function return code (unused)

  print '(a)', 'ERROR: sigsetjmp() is not implemented in this build.'
  print '(a)', 'Recompile with autoconf or -DSIGSETJMP_NAME=\"<symbol name>\".'
  error stop

  ! NOTE: compilers may expect a return value, even if it is unreachable
  read env%state
  read char(savesigs)
  rc = -1
end function sigsetjmp_missing

!> Placeholder function for a missing or unconfigured siglongjmp
subroutine siglongjmp_missing(env, val) bind(c)
  type(sigjmp_buf), intent(in) :: env
    !< Current process state (unused)
  integer(kind=c_int), value, intent(in) :: val
    !< Enable signal state flag (unused)

  print '(a)', 'ERROR: siglongjmp() is not implemented in this build.'
  print '(a)', 'Recompile with autoconf or -DSIGLONGJMP_NAME=\"<symbol name>\".'
  read env%state
  read char(val)
  error stop
end subroutine siglongjmp_missing

end module posix
