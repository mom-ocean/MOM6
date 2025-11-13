#ifndef DO_CONCURRENT_COMPAT_H_
#define DO_CONCURRENT_COMPAT_H_

! This macro conditionally applies the locality specifier of a do concurrent
! loop if supported by the compiler.
#ifdef HAVE_FC_DO_CONCURRENT_LOCAL
#define DO_LOCALITY(X) X
#else
#define DO_LOCALITY(X) ;
#endif

#endif
