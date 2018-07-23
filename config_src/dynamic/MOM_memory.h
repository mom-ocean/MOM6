!/// \brief Compile-time memory settings
!/// \details This include file determines the compile-time memory settings.
!/// There are several variants of this file and only one should be in the search path for compilation.
!/// \file MOM_memory.h

!/// The number of thickness grid points in the i-direction of the global domain.
#define NIGLOBAL_ NONSENSE_NIGLOBAL
!/// The number of thickness grid points in the j-direction of the global domain.
#define NJGLOBAL_ NONSENSE_NJGLOBAL
!/// The number of layers in the vertical direction.
#define NK_ NONSENSE_NK

!/// The number of processors in the i-direction.
#define NIPROC_ NONSENSE_NIPROC

!/// The number of processors in the j-direction.
#define NJPROC_ NONSENSE_NJPROC

!/// The maximum permitted number (each) of restart variables, time derivatives, etc.
!/// This is mostly used for the size of pointer arrays, so it should be set generously.
#ifndef MAX_FIELDS_
#define MAX_FIELDS_ 50
#endif

!/// The number of memory halo cells on each side of the computational domain in the i-direction.
#define NIHALO_ 2

!/// The number of memory halo cells on each side of the computational domain in the j-direction.
#define NJHALO_ 2

!/// If SYMMETRIC_MEMORY_() is defined, the velocity point data domain includes every face of the thickness points.
!/// In other words, some arrays are larger than others, depending on where they are on the staggered grid.
#undef  SYMMETRIC_MEMORY_

!/// If STATIC_MEMORY_ is defined, the principle variables have sizes that are statically determined at compile time.
!/// Otherwise the sizes are not determined until run time.
#undef STATIC_MEMORY_

#include <MOM_memory_macros.h>
