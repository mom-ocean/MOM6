!*******+*********+*********+*********+*********+*********+*********+
!                                                                    
!  By Robert Hallberg, March 1999                                    
!                                                                    
! This header file is necessary to define macros which allow the     
! unique metric terms used by MOM_hor_visc.F90 to vary in 0, 1, or 2       
! horizontal dimensions without taking up more memory than is        
! necessary.                                                         
! The names of these 20 macros should be a self-evident description  
! of what they are. These are very similar to the macros defined     
! in MOM_grid_macros.h.                                                      
!*******+*********+*********+*********+*********+*********+*********+

#ifdef CARTESIAN
# define I_SZ C1_
# define J_SZ C1_
# define I_METRIC_SZ 1
# define J_METRIC_SZ 1
# define IQ_SZ C1_
# define JQ_SZ C1_
# define I_METRICQ_SZ 1
# define J_METRICQ_SZ 1
#else
! not CARTESIAN
# ifdef XYMETRICS_I_J
#   define I_SZ NXMEM_
#   define J_SZ NYMEM_
#   define I_METRIC_SZ isd:ied
#   define J_METRIC_SZ jsd:jed
#   define IQ_SZ NXMEMQP_
#   define JQ_SZ NYMEMQP_
#   define I_METRICQ_SZ Isdq:Iedq
#   define J_METRICQ_SZ Jsdq:Jedq
# else
!not XYMETRICS_I_J
#  ifdef XYMETRICS_I
#   define I_SZ NXMEM_
#   define J_SZ C1_
#   define I_METRIC_SZ isd:ied
#   define J_METRIC_SZ 1
#   define IQ_SZ NXMEMQP_
#   define JQ_SZ C1_
#   define I_METRICQ_SZ Isdq:Iedq
#   define J_METRICQ_SZ 1
#  else
! must be XYMETRICS_J 
#   define I_SZ C1_
#   define J_SZ NYMEM_
#   define I_METRIC_SZ 1
#   define J_METRIC_SZ jsd:jed
#   define IQ_SZ C1_
#   define JQ_SZ NYMEMQP_
#   define I_METRICQ_SZ 1
#   define J_METRICQ_SZ Jsdq:Jedq
#  endif
! XYMETRICS_I
# endif
! XYMETRICS_I_J
#endif
! CARTESIAN

#ifdef CARTESIAN
#   define DX2h(i,j)       dx2h(1,1)
#   define DX2q(i,j)       dx2q(1,1)
#   define DY2h(i,j)       dy2h(1,1)
#   define DY2q(i,j)       dy2q(1,1)
#   define DX_DYh(i,j)   dx_dyh(1,1)
#   define DX_DYq(i,j)   dx_dyq(1,1)
#   define DY_DXh(i,j)   dy_dxh(1,1)
#   define DY_DXq(i,j)   dy_dxq(1,1)
#   define IDX2DYu(i,j) Idx2dyu(1,1)
#   define IDX2DYv(i,j) Idx2dyv(1,1)
#   define IDXDY2u(i,j) Idxdy2u(1,1)
#   define IDXDY2v(i,j) Idxdy2v(1,1)

#   define LAPLAC_CONST_xx(i,j) Laplac_Const_xx(1,1)
#   define LAPLAC_CONST_xy(i,j) Laplac_Const_xy(1,1)
#   define BIHARM_CONST_xx(i,j) Biharm_Const_xx(1,1)
#   define BIHARM_CONST_xy(i,j) Biharm_Const_xy(1,1)
#else
! not CARTESIAN
# ifdef XYMETRICS_I_J
#   define DX2h(i,j)       dx2h(i,j)
#   define DX2q(i,j)       dx2q(i,j)
#   define DY2h(i,j)       dy2h(i,j)
#   define DY2q(i,j)       dy2q(i,j)
#   define DX_DYh(i,j)   dx_dyh(i,j)
#   define DX_DYq(i,j)   dx_dyq(i,j)
#   define DY_DXh(i,j)   dy_dxh(i,j)
#   define DY_DXq(i,j)   dy_dxq(i,j)
#   define IDX2DYu(i,j) Idx2dyu(i,j)
#   define IDX2DYv(i,j) Idx2dyv(i,j)
#   define IDXDY2u(i,j) Idxdy2u(i,j)
#   define IDXDY2v(i,j) Idxdy2v(i,j)

#   define LAPLAC_CONST_xx(i,j) Laplac_Const_xx(i,j)
#   define LAPLAC_CONST_xy(i,j) Laplac_Const_xy(i,j)
#   define BIHARM_CONST_xx(i,j) Biharm_Const_xx(i,j)
#   define BIHARM_CONST_xy(i,j) Biharm_Const_xy(i,j)
# else
!not XYMETRICS_I_J
#  ifdef XYMETRICS_I
#   define DX2h(i,j)       dx2h(i,1)
#   define DX2q(i,j)       dx2q(i,1)
#   define DY2h(i,j)       dy2h(i,1)
#   define DY2q(i,j)       dy2q(i,1)
#   define DX_DYh(i,j)   dx_dyh(i,1)
#   define DX_DYq(i,j)   dx_dyq(i,1)
#   define DY_DXh(i,j)   dy_dxh(i,1)
#   define DY_DXq(i,j)   dy_dxq(i,1)
#   define IDX2DYu(i,j) Idx2dyu(i,1)
#   define IDX2DYv(i,j) Idx2dyv(i,1)
#   define IDXDY2u(i,j) Idxdy2u(i,1)
#   define IDXDY2v(i,j) Idxdy2v(i,1)

#   define LAPLAC_CONST_xx(i,j) Laplac_Const_xx(i,1)
#   define LAPLAC_CONST_xy(i,j) Laplac_Const_xy(i,1)
#   define BIHARM_CONST_xx(i,j) Biharm_Const_xx(i,1)
#   define BIHARM_CONST_xy(i,j) Biharm_Const_xy(i,1)
#  else
! must be XYMETRICS_J 
#   define DX2h(i,j)       dx2h(1,j)
#   define DX2q(i,j)       dx2q(1,j)
#   define DY2h(i,j)       dy2h(1,j)
#   define DY2q(i,j)       dy2q(1,j)
#   define DX_DYh(i,j)   dx_dyh(1,j)
#   define DX_DYq(i,j)   dx_dyq(1,j)
#   define DY_DXh(i,j)   dy_dxh(1,j)
#   define DY_DXq(i,j)   dy_dxq(1,j)
#   define IDX2DYu(i,j) Idx2dyu(1,j)
#   define IDX2DYv(i,j) Idx2dyv(1,j)
#   define IDXDY2u(i,j) Idxdy2u(1,j)
#   define IDXDY2v(i,j) Idxdy2v(1,j)

#   define LAPLAC_CONST_xx(i,j) Laplac_Const_xx(1,j)
#   define LAPLAC_CONST_xy(i,j) Laplac_Const_xy(1,j)
#   define BIHARM_CONST_xx(i,j) Biharm_Const_xx(1,j)
#   define BIHARM_CONST_xy(i,j) Biharm_Const_xy(1,j)
#  endif
! XYMETRICS_I
# endif
! XYMETRICS_I_J
#endif
! CARTESIAN
