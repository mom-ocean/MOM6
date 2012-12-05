!********+*********+*********+*********+*********+*********+*********+*
!*   This include file determines the compile-time settings for the   *
!* Generalized Ocean Layered Dynamics (GOLD) ocean model.             *
!********+*********+*********+*********+*********+*********+*********+*

!  Specify the appropriate dimensionality for the metrics.
#undef  CARTESIAN
                               !    Use a uniform Cartesian grid if CARTESIAN
                               !  is defined, otherwise use a latitude-longitude
                               !  coordinate grid.
#define XMETRIC_J
#define XMETRIC_I
#define YMETRIC_J
#define YMETRIC_I
                               !    Define XMETRIC_J if the x-direction metrics
                               !  vary in the y- (or j^) direction.  Otherwise
                               !  undefine XMETRIC_J.  XMETRIC_I, YMETRIC_J,
                               !  and YMETRIC_I are used similarly.
                               !  CARTESIAN overrides all of these choices.
                               !    For example, on a regular latitude/longitude
                               !  grid, define XMETRIC_J and undefine the rest.
                               !  For a Mercator grid, define only XMETRIC_J and
                               !  YMETRIC_J.  Defining all of these always
                               !  works, but the model will be slower.

!  Specify the numerical domain.
#define NXTOT 360
#define NYTOT 210
                               !    NXTOT and NYTOT are the number of thickness
                               !  grid points in the zonal and meridional
                               !  directions of the physical domain.
#define NZ 63
                               !    The number of layers.

#define STATIC_MEMORY
                               !    If STATIC_MEMORY is defined, the principle
                               !  variables will have sizes that are statically
                               !  determined at compile time.  Otherwise the
                               !  sizes are not determined until run time. The
                               !  STATIC option is substantially faster, but
                               !  does not allow the PE count to be changed at
                               !  run time.
#undef  SYMMETRIC_MEMORY
                               !    If defined, the velocity point data domain
                               !  includes every face of the thickness points.
                               !  In other words, some arrays are larger than
                               !  others, depending on where they are on the 
                               !  staggered grid.

#define NXPROC 10
                               !    NXPROC is the number of processors in the
                               !  x-direction.
#define NYPROC 6
                               !    NYPROC is the number of processors in the
                               !  y-direction.

#define MAX_FIELDS 80
                               !    The maximum permitted number (each) of
                               !  restart variables, time derivatives, etc.
                               !  This is mostly used for the size of pointer
                               !  arrays, so it should be set generously.

#define EPSILON 1.0e-10
                               !    The default minimum layer thickness, in m.
                               !  This is superseded by setting ANGSTROM in the
                               !  GOLD parameter file.

#define NX_HALO 4
#define NY_HALO 4
                               !   NX_HALO and NY_HALO are the sizes of the
                               ! memory halos on each side.
#define BT_HALO 10
                               !   BT_HALO is the size of the memory halos in
                               ! the barotropic solver.

#include <GOLD_memory_macros.h>
#include <GOLD_grid_macros.h>
