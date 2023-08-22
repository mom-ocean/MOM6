#ifndef MOM6_POSIX_H_
#define MOM6_POSIX_H_

! JMP_BUF_SIZE should be set to sizeof(jmp_buf).
! If unset, then use a typical glibc value (25 long ints)
#ifndef SIZEOF_JMP_BUF
#define SIZEOF_JMP_BUF 200
#endif

! If unset, assume jmp_buf and sigjmp_buf are equivalent (as in glibc).
#ifndef SIZEOF_SIGJMP_BUF
#define SIZEOF_SIGJMP_BUF SIZEOF_JMP_BUF
#endif

! Wrappers to <setjmp.h> are disabled on default.
#ifndef SETJMP_NAME
#define SETJMP_NAME "setjmp_missing"
#endif

#ifndef LONGJMP_NAME
#define LONGJMP_NAME "longjmp_missing"
#endif

#ifndef SIGSETJMP_NAME
#define SIGSETJMP_NAME "sigsetjmp_missing"
#endif

#ifndef SIGLONGJMP_NAME
#define SIGLONGJMP_NAME "siglongjmp_missing"
#endif

! This should be defined by <signal.h>;
! If unset, we use the most common (x86) value
#ifndef POSIX_SIGUSR1
#define POSIX_SIGUSR1 10
#endif

#endif
