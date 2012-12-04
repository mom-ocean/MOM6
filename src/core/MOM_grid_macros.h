! ********+*********+*********+*********+*********+*********+*********+*
! *    This header file defines the metric arrays for MOM to be of     *
! *  appropriate sizes.                                                *
! ********+*********+*********+*********+*********+*********+*********+*

# if defined(XMETRIC_J) && defined(XMETRIC_I)
#  define XMETRIC_I_J
# endif
# if defined(YMETRIC_J) && defined(YMETRIC_I)
#  define YMETRIC_I_J
# endif

# if defined(XMETRIC_I_J) || defined(YMETRIC_I_J) || \
    (defined(XMETRIC_J) && defined(YMETRIC_I)) || \
    (defined(XMETRIC_I) && defined(YMETRIC_J))
#  define XYMETRICS_I_J
# else
#  if defined(XMETRIC_J) || defined(YMETRIC_J)
#   define XYMETRICS_J
#  else
#   if defined(XMETRIC_I) || defined(YMETRIC_I)
#    define XYMETRICS_I
#   endif
!  XMETRIC_I or YMETRIC_I
#  endif
!  XMETRIC_J YMETRIC_J
# endif
!  XMETRIC_I_J or YMETRIC_I_J


# ifdef XMETRIC_I_J
#    define XMETRIC_NX_S NXMEM_
#    define XMETRIC_NY_S NYMEM_
#    define XMETRIC_NX_D isd:ied
#    define XMETRIC_NY_D jsd:jed
#    define XMETRICQ_NX_S NXMEMQP_
#    define XMETRICQ_NY_S NYMEMQP_
#    define XMETRICQ_NX_D Isdq:Iedq
#    define XMETRICQ_NY_D Jsdq:Jedq
  
#    define DXh(i,j)   dxh(i,j)
#    define IDXh(i,j) Idxh(i,j)
#    define DXq(i,j)   dxq(i,j)
#    define IDXq(i,j) Idxq(i,j)
#    define DXu(i,j)   dxu(i,j)
#    define IDXu(i,j) Idxu(i,j)
#    define DXv(i,j)   dxv(i,j)
#    define IDXv(i,j) Idxv(i,j)
# else
#  ifdef XMETRIC_J
#    define XMETRIC_NX_S C1_
#    define XMETRIC_NY_S NYMEM_
#    define XMETRIC_NX_D 1
#    define XMETRIC_NY_D jsd:jed
#    define XMETRICQ_NX_S C1_
#    define XMETRICQ_NY_S NYMEMQP_
#    define XMETRICQ_NX_D 1
#    define XMETRICQ_NY_D Jsdq:Jedq

#    define DXh(i,j)   dxh(1,j)
#    define IDXh(i,j) Idxh(1,j)
#    define DXq(i,j)   dxq(1,j)
#    define IDXq(i,j) Idxq(1,j)
#    define DXu(i,j)   dxu(1,j)
#    define IDXu(i,j) Idxu(1,j)
#    define DXv(i,j)   dxv(1,j)
#    define IDXv(i,j) Idxv(1,j)
#  else
#   ifdef XMETRIC_I
#    define XMETRIC_NX_S NXMEM_
#    define XMETRIC_NY_S C1_
#    define XMETRIC_NX_D isd:ied
#    define XMETRIC_NY_D 1
#    define XMETRICQ_NX_S NXMEMQP_
#    define XMETRICQ_NY_S C1_
#    define XMETRICQ_NX_D Isdq:Iedq
#    define XMETRICQ_NY_D 1

#    define DXh(i,j)   dxh(i,1)
#    define IDXh(i,j) Idxh(i,1)
#    define DXq(i,j)   dxq(i,1)
#    define IDXq(i,j) Idxq(i,1)
#    define DXu(i,j)   dxu(i,1)
#    define IDXu(i,j) Idxu(i,1)
#    define DXv(i,j)   dxv(i,1)
#    define IDXv(i,j) Idxv(i,1)
#   else
#    define XMETRIC_NX_S C1_
#    define XMETRIC_NY_S C1_
#    define XMETRIC_NX_D 1
#    define XMETRIC_NY_D 1
#    define XMETRICQ_NX_S C1_
#    define XMETRICQ_NY_S C1_
#    define XMETRICQ_NX_D 1
#    define XMETRICQ_NY_D 1

#    define DXh(i,j)   dxh(1,1)
#    define IDXh(i,j) Idxh(1,1)
#    define DXq(i,j)   dxq(1,1)
#    define IDXq(i,j) Idxq(1,1)
#    define DXu(i,j)   dxu(1,1)
#    define IDXu(i,j) Idxu(1,1)
#    define DXv(i,j)   dxv(1,1)
#    define IDXv(i,j) Idxv(1,1)
#   endif
!  XMETRIC_I
#  endif
!  XMETRIC_J
# endif
!  XMETRIC_I_J

# ifdef YMETRIC_I_J
#    define YMETRIC_NX_S NXMEM_
#    define YMETRIC_NY_S NYMEM_
#    define YMETRIC_NX_D isd:ied
#    define YMETRIC_NY_D jsd:jed
#    define YMETRICQ_NX_S NXMEMQP_
#    define YMETRICQ_NY_S NYMEMQP_
#    define YMETRICQ_NX_D Isdq:Iedq
#    define YMETRICQ_NY_D Jsdq:Jedq

#    define DYh(i,j)  dyh(i,j)
#    define IDYh(i,j) Idyh(i,j)
#    define DYq(i,j)  dyq(i,j)
#    define IDYq(i,j) Idyq(i,j)
#    define DYu(i,j)  dyu(i,j)
#    define IDYu(i,j) Idyu(i,j)
#    define DYv(i,j)  dyv(i,j)
#    define IDYv(i,j) Idyv(i,j)
# else
#  ifdef YMETRIC_J
#    define YMETRIC_NX_S C1_
#    define YMETRIC_NY_S NYMEM_
#    define YMETRIC_NX_D 1
#    define YMETRIC_NY_D jsd:jed
#    define YMETRICQ_NX_S C1_
#    define YMETRICQ_NY_S NYMEMQP_
#    define YMETRICQ_NX_D 1
#    define YMETRICQ_NY_D Jsdq:Jedq

#    define DYh(i,j)  dyh(1,j)
#    define IDYh(i,j) Idyh(1,j)
#    define DYq(i,j)  dyq(1,j)
#    define IDYq(i,j) Idyq(1,j)
#    define DYu(i,j)  dyu(1,j)
#    define IDYu(i,j) Idyu(1,j)
#    define DYv(i,j)  dyv(1,j)
#    define IDYv(i,j) Idyv(1,j)
#  else
#   ifdef YMETRIC_I
#    define YMETRIC_NX_S NXMEM_
#    define YMETRIC_NY_S C1_
#    define YMETRIC_NX_D isd:ied
#    define YMETRIC_NY_D 1
#    define YMETRICQ_NX_S NXMEMQP_
#    define YMETRICQ_NY_S C1_
#    define YMETRICQ_NX_D Isdq:Iedq
#    define YMETRICQ_NY_D 1

#    define DYh(i,j)  dyh(i,1)
#    define IDYh(i,j) Idyh(i,1)
#    define DYq(i,j)  dyq(i,1)
#    define IDYq(i,j) Idyq(i,1)
#    define DYu(i,j)  dyu(i,1)
#    define IDYu(i,j) Idyu(i,1)
#    define DYv(i,j)  dyv(i,1)
#    define IDYv(i,j) Idyv(i,1)
#   else
#    define YMETRIC_NX_S C1_
#    define YMETRIC_NY_S C1_
#    define YMETRIC_NX_D 1
#    define YMETRIC_NY_D 1
#    define YMETRICQ_NX_S C1_
#    define YMETRICQ_NY_S C1_
#    define YMETRICQ_NX_D 1
#    define YMETRICQ_NY_D 1

#    define DYh(i,j)  dyh(1,1)
#    define IDYh(i,j) Idyh(1,1)
#    define DYq(i,j)  dyq(1,1)
#    define IDYq(i,j) Idyq(1,1)
#    define DYu(i,j)  dyu(1,1)
#    define IDYu(i,j) Idyu(1,1)
#    define DYv(i,j)  dyv(1,1)
#    define IDYv(i,j) Idyv(1,1)
#   endif
!  YMETRIC_I
#  endif
!  YMETRIC_J
# endif
! YMETRIC_I_J

# ifdef XYMETRICS_I_J
#    define XYMETRIC_NX_S NXMEM_
#    define XYMETRIC_NY_S NYMEM_
#    define XYMETRIC_NX_D isd:ied
#    define XYMETRIC_NY_D jsd:jed
#    define XYMETRICQ_NX_S NXMEMQP_
#    define XYMETRICQ_NY_S NYMEMQP_
#    define XYMETRICQ_NX_D Isdq:Iedq
#    define XYMETRICQ_NY_D Jsdq:Jedq

#    define DXDYh(i,j)   dxdyh(i,j)
#    define IDXDYh(i,j) Idxdyh(i,j)
#    define DXDYq(i,j)   dxdyq(i,j)
#    define IDXDYq(i,j) Idxdyq(i,j)
# else
#  ifdef XYMETRICS_J
#    define XYMETRIC_NX_S C1_
#    define XYMETRIC_NY_S NYMEM_
#    define XYMETRIC_NX_D 1
#    define XYMETRIC_NY_D jsd:jed
#    define XYMETRICQ_NX_S C1_
#    define XYMETRICQ_NY_S NYMEMQP_
#    define XYMETRICQ_NX_D 1
#    define XYMETRICQ_NY_D Jsdq:Jedq

#    define DXDYh(i,j)   dxdyh(1,j)
#    define IDXDYh(i,j) Idxdyh(1,j)
#    define DXDYq(i,j)   dxdyq(1,j)
#    define IDXDYq(i,j) Idxdyq(1,j)
#  else
#   ifdef XYMETRICS_I
#    define XYMETRIC_NX_S NXMEM_
#    define XYMETRIC_NY_S C1_
#    define XYMETRIC_NX_D isd:ied
#    define XYMETRIC_NY_D 1
#    define XYMETRICQ_NX_S NXMEMQP_
#    define XYMETRICQ_NY_S C1_
#    define XYMETRICQ_NX_D Isdq:Iedq
#    define XYMETRICQ_NY_D 1

#    define DXDYh(i,j)   dxdyh(i,1)
#    define IDXDYh(i,j) Idxdyh(i,1)
#    define DXDYq(i,j)   dxdyq(i,1)
#    define IDXDYq(i,j) Idxdyq(i,1)
#   else
#    define XYMETRIC_NX_S C1_
#    define XYMETRIC_NY_S C1_
#    define XYMETRIC_NX_D 1
#    define XYMETRIC_NY_D 1
#    define XYMETRICQ_NX_S C1_
#    define XYMETRICQ_NY_S C1_
#    define XYMETRICQ_NX_D 1
#    define XYMETRICQ_NY_D 1

#    define DXDYh(i,j)   dxdyh(1,1)
#    define IDXDYh(i,j) Idxdyh(1,1)
#    define DXDYq(i,j)   dxdyq(1,1)
#    define IDXDYq(i,j) Idxdyq(1,1)
#   endif
!  XYMETRICS_I
#  endif
!  XYMETRICS_J
# endif
!  XYMETRICS_I_J
