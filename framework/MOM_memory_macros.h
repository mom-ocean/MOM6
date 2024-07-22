!//! \brief Memory macros
!//! \details This is a header file to define macros for static and dynamic memory allocation.
!//! Define STATIC_MEMORY_ in MOM_memory.h for static memory allocation.
!//! Otherwise dynamic memory allocation will be assumed.
!//!
!//! For explanation of symmetric and non-symmetric memory modes see \ref Horizontal_Indexing.
!//! \file MOM_memory_macros.h

#ifdef STATIC_MEMORY_
!/* Static memory allocation section */

!/// Deallocates array x when using dynamic memory mode. Does nothing in static memory mode.
#  define DEALLOC_(x)
!/// Allocates array x when using dynamic memory mode. Does nothing in static memory mode.
#  define ALLOC_(x)
!/// Attaches the ALLOCATABLE attribute to an array in dynamic memory mode. Does nothing in static memory mode.
#  define ALLOCABLE_
!/// Attaches the POINTER attribute to an array in dynamic memory mode. Does nothing in static memory mode.
#  define PTR_
!/// Nullify a pointer in dynamic memory mode. Does nothing in static memory mode.
#  define TO_NULL_

!/* These are the macros that should be used when setting up ALLOCABLE_ or PTR_ (heap) variables. */

!/// Expands to : in dynamic memory mode, or is the i-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at h- or v- points.
#  define NIMEM_ (((NIGLOBAL_-1)/NIPROC_)+1+2*NIHALO_)
!/// Expands to : in dynamic memory mode, or is the j-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at h- or u- points.
#  define NJMEM_ (((NJGLOBAL_-1)/NJPROC_)+1+2*NJHALO_)

#  ifdef SYMMETRIC_MEMORY_
!/// Expands to : or 0: in dynamic memory mode, or is the staggered i-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at q- or u- points.
#    define NIMEMB_   0:NIMEM_
!/// Expands to : or 0: in dynamic memory mode, or is the staggered j-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at q- or v- points.
#    define NJMEMB_   0:NJMEM_
#  else
!/// Expands to : or 0: in dynamic memory mode, or is the staggered i-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at q- or u- points.
#    define NIMEMB_   NIMEM_
!/// Expands to : or 0: in dynamic memory mode, or is the staggered j-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at q- or v- points.
#    define NJMEMB_   NJMEM_
#  endif
!/// Expands to : in dynamic memory mode, or to NIMEMB_ in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at h- or v- points.
#  define NIMEMB_PTR_ NIMEMB_
!/// Expands to : in dynamic memory mode, or to NJMEMB_ in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at h- or u- points.
#  define NJMEMB_PTR_ NJMEMB_
!/// Expands to 0: in dynamic memory mode, or is the staggered i-shape of a tile in static memory mode.
!/// Use for always-symmetric heap (ALLOCABLE_ or PTR_) variables at q- or u- points.
#  define NIMEMB_SYM_ 0:NIMEM_
!/// Expands to 0: in dynamic memory mode, or is the staggered j-shape of a tile in static memory mode.
!/// Use for always-symmetric heap (ALLOCABLE_ or PTR_) variables at q- or v- points.
#  define NJMEMB_SYM_ 0:NJMEM_
!/// Expands to : in dynamic memory mode or is to the number of layers in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) layer variables.
#  define NKMEM_      NK_
!/// Expands to 0: in dynamic memory mode or to 0:NK_ in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) interface variables.
#  define NKMEM0_     0:NK_
!/// Expands to : in dynamic memory mode or to NK_+1 in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) interface variables.
#  define NK_INTERFACE_ NK_+1
!/// Expands to : or 1. UNKNOWN PURPOSE!
#  define C1_         1
!/// Expands to : or 2. UNKNOWN PURPOSE!
#  define C2_         2
!/// Expands to : or 3. UNKNOWN PURPOSE!
#  define C3_         3

!/* These are the macros that should be used for subroutine arguments or for automatically allocated (stack) variables. */

!/// The i-shape of a dummy argument staggered at h- or v-points.
#  define SZI_(G)     NIMEM_
!/// The j-shape of a dummy argument staggered at h- or u-points.
#  define SZJ_(G)     NJMEM_
!/// The k-shape of a layer dummy argument.
#  define SZK_(G)     NK_
!/// The k-shape of an interface dummy argument.
#  define SZK0_(G)    0:NK_
!/// The i-shape of a dummy argument staggered at q- or u-points.
#  define SZIB_(G)    NIMEMB_
!/// The j-shape of a dummy argument staggered at q- or v-points.
#  define SZJB_(G)    NJMEMB_
!/// The i-shape of a symmetric dummy argument staggered at q- or u-points.
#  define SZIBS_(G)   0:NIMEM_
!/// The j-shape of a symmetric dummy argument staggered at q- or v-points.
#  define SZJBS_(G)   0:NJMEM_

#else
!/* Dynamic memory allocation section */

!/// Deallocates array x when using dynamic memory mode. Does nothing in static memory mode.
#  define DEALLOC_(x) deallocate(x)
!/// Allocates array x when using dynamic memory mode. Does nothing in static memory mode.
#  define ALLOC_(x)   allocate(x)
!/// Attaches the ALLOCATABLE attribute to an array in dynamic memory mode. Does nothing in static memory mode.
#  define ALLOCABLE_ ,allocatable
!/// Attaches the POINTER attribute to an array in dynamic memory mode. Does nothing in static memory mode.
#  define PTR_       ,pointer
!/// Nullify a pointer in dynamic memory mode. Does nothing in static memory mode.
#  define TO_NULL_   =>NULL()

!/* These are the macros that should be used when setting up ALLOCABLE_ or PTR_ (heap) variables. */

!/// Expands to : in dynamic memory mode, or is the i-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at h- or v- points.
#  define NIMEM_      :
!/// Expands to : in dynamic memory mode, or is the j-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at h- or u- points.
#  define NJMEM_      :
!/// Expands to : in dynamic memory mode, or to NIMEMB_ in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at h- or v- points.
#  define NIMEMB_PTR_ :
!/// Expands to : in dynamic memory mode, or to NJMEMB_ in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at h- or u- points.
#  define NJMEMB_PTR_ :
#  ifdef SYMMETRIC_MEMORY_
!/// Expands to : or 0: in dynamic memory mode, or is the staggered i-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at q- or u- points.
#    define NIMEMB_   0:
!/// Expands to : or 0: in dynamic memory mode, or is the staggered j-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at q- or v- points.
#    define NJMEMB_   0:
#  else
!/// Expands to : or 0: in dynamic memory mode, or is the staggered i-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at q- or u- points.
#    define NIMEMB_   :
!/// Expands to : or 0: in dynamic memory mode, or is the staggered j-shape of a tile in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) variables at q- or v- points.
#    define NJMEMB_   :
#  endif
!/// Expands to 0: in dynamic memory mode, or is the staggered i-shape of a tile in static memory mode.
!/// Use for always-symmetric heap (ALLOCABLE_ or PTR_) variables at q- or u- points.
#  define NIMEMB_SYM_ 0:
!/// Expands to 0: in dynamic memory mode, or is the staggered j-shape of a tile in static memory mode.
!/// Use for always-symmetric heap (ALLOCABLE_ or PTR_) variables at q- or v- points.
#  define NJMEMB_SYM_ 0:
!/// Expands to : in dynamic memory mode or is to the number of layers in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) layer variables.
#  define NKMEM_      :
!/// Expands to 0: in dynamic memory mode or to 0:NK_ in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) interface variables.
#  define NKMEM0_     0:
!/// Expands to : in dynamic memory mode or to NK_+1 in static memory mode.
!/// Use for heap (ALLOCABLE_ or PTR_) interface variables.
#  define NK_INTERFACE_  :
!/// Expands to : or 1. UNKNOWN PURPOSE!
#  define C1_         :
!/// Expands to : or 2. UNKNOWN PURPOSE!
#  define C2_         :
!/// Expands to : or 3. UNKNOWN PURPOSE!
#  define C3_         :

!/// \todo Explain or remove C1_, C2_ and C3_

!/* These are the macros that should be used for subroutine arguments or for automatically allocated (stack) variables. */

!/// The i-shape of a dummy argument staggered at h- or v-points.
#  define SZI_(G)     G%isd:G%ied
!/// The j-shape of a dummy argument staggered at h- or u-points.
#  define SZJ_(G)     G%jsd:G%jed
!/// The k-shape of a layer dummy argument.
#  define SZK_(G)     G%ke
!/// The k-shape of an interface dummy argument.
#  define SZK0_(G)    0:G%ke
!/// The i-shape of a dummy argument staggered at q- or u-points.
#  define SZIB_(G)    G%IsdB:G%IedB
!/// The j-shape of a dummy argument staggered at q- or v-points.
#  define SZJB_(G)    G%JsdB:G%JedB
!/// The i-shape of a symmetric dummy argument staggered at q- or u-points.
#  define SZIBS_(G)   G%isd-1:G%ied
!/// The j-shape of a symmetric dummy argument staggered at q- or v-points.
#  define SZJBS_(G)   G%jsd-1:G%jed

#endif

!/* These dynamic size macros always give the same results (for now). */

!/// The i-shape of a dynamic dummy argument staggered at h- or v-points.
#define SZDI_(G)  G%isd:G%ied
!/// The i-shape of a dynamic dummy argument staggered at q- or u-points.
#define SZDIB_(G) G%IsdB:G%IedB
!/// The j-shape of a dynamic dummy argument staggered at h- or u-points.
#define SZDJ_(G)  G%jsd:G%jed
!/// The j-shape of a dynamic dummy argument staggered at q- or v-points.
#define SZDJB_(G) G%JsdB:G%JedB
