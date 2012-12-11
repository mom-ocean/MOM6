!***********************************************************************
! This is a header file to define macros for static and dynamic memory *
! allocation.  Define STATIC_MEMORY in MOM_memory.h for static memory *
! allocation.  Otherwise dynamic memory allocation will be assumed.    *
!***********************************************************************

#if (defined(X1) && !defined(NX_HALO))
#  define NX_HALO X1
#endif
#if (defined(Y1) && !defined(NY_HALO))
#  define NY_HALO Y1
#endif

#ifdef STATIC_MEMORY
!   NXMEM and NYMEM are the maximum number of grid points in the
! x- and y-directions on each processsor.
#  define NXMEM (((NXTOT-1)/NXPROC)+1+2*NX_HALO)
#  define NYMEM (((NYTOT-1)/NYPROC)+1+2*NY_HALO)

#  define DEALLOC(x)
#  define ALLOC(x)
#  define ALLOCABLE
#  define PTR_
#  define TO_NULL_
#  define NXMEM_     NXMEM
#  define NYMEM_     NYMEM
#  ifdef SYMMETRIC_MEMORY
#    define NXMEMQ_    0:NXMEM
#    define NYMEMQ_    0:NYMEM
#  else
#    define NXMEMQ_    NXMEM
#    define NYMEMQ_    NYMEM
#  endif
#  define NXMEMQP_    NXMEMQ_
#  define NYMEMQP_    NYMEMQ_
#  define NXMEMQS_    0:NXMEM
#  define NYMEMQS_    0:NYMEM
#  define NKMEM_     NZ
#  define NK_INTERFACE_ NZ+1
#  define C1_        1
#  define C2_        2
#  define C3_        3
#  define SZ1_(u)    NXMEM
#  define SZ2_(u)    NYMEM
#  define SZ3_(u)    NZ
#  define SZI_(G)    NXMEM
#  define SZJ_(G)    NYMEM
#  define SZK_(G)    NZ
#  define SZIQ_(G)   NXMEMQ_
#  define SZJQ_(G)   NYMEMQ_
#  define SZIQS_(G)  0:NXMEM
#  define SZJQS_(G)  0:NYMEM

#else
! Dynamic memory allocation

#  define DEALLOC(x) deallocate(x)
#  define ALLOC(x)   allocate(x)
#  define ALLOCABLE  ,allocatable
#  define PTR_       ,pointer
#  define TO_NULL_   =>NULL()
#  define NXMEM_     :
#  define NYMEM_     :
#  define NXMEMQP_   :
#  define NYMEMQP_   :
#  ifdef SYMMETRIC_MEMORY
#    define NXMEMQ_    0:
#    define NYMEMQ_    0:
#  else
#    define NXMEMQ_    :
#    define NYMEMQ_    :
#  endif
#  define NXMEMQS_    0:
#  define NYMEMQS_    0:
#  define NKMEM_     :
#  define NK_INTERFACE_  :
#  define C1_        :
#  define C2_        :
#  define C3_        :
#  define SZ1_(u)    lbound(u,1):ubound(u,1)
#  define SZ2_(u)    lbound(u,2):ubound(u,2)
#  define SZ3_(u)    lbound(u,3):ubound(u,3)
#  define SZI_(G)    G%isd:G%ied
#  define SZJ_(G)    G%jsd:G%jed
#  define SZK_(G)    G%ks:G%ke
#  define SZIQ_(G)   G%Isdq:G%Iedq
#  define SZJQ_(G)   G%Jsdq:G%Jedq
#  define SZIQS_(G)  G%isd-1:G%ied
#  define SZJQS_(G)  G%jsd-1:G%jed

#endif
