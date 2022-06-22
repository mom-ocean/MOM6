#ifndef MOM6_POSIX_H_
#define MOM6_POSIX_H_

! JMP_BUF_SIZE should be set to sizeof(jmp_buf).
! If unset, then use a typical glibc value (25 long ints)
#ifndef JMP_BUF_SIZE
#define JMP_BUF_SIZE 200
#endif

! If unset, assume jmp_buf and sigjmp_buf are equivalent (as in glibc).
#ifndef SIGJMP_BUF_SIZE
#define SIGJMP_BUF_SIZE JMP_BUF_SIZE
#endif

! glibc defines sigsetjmp as __sigsetjmp via macro readable from <setjmp.h>.
! Perhaps autoconf can configure this one...
! TODO: Need a solution here!
#ifndef SIGSETJMP_NAME
#define SIGSETJMP_NAME "__sigsetjmp"
#endif

! This should be defined by /usr/include/signal.h
! If unset, we use the most common (x86) value
#ifndef POSIX_SIGUSR1
#define POSIX_SIGUSR1 10
#endif

#endif
