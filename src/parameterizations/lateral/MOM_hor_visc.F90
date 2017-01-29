module MOM_hor_visc
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, April 1994 - June 2002.                        *
!*                                                                     *
!*    This program contains the subroutine that calculates the         *
!*  effects of horizontal viscosity, including parameterizations of    *
!*  the value of the viscosity itself. horizontal_viscosity calc-      *
!*  ulates the acceleration due to some combination of a biharmonic    *
!*  viscosity and a Laplacian viscosity. Either or both may use a      *
!*  coefficient that depends on the shear and strain of the flow.      *
!*  All metric terms are retained.  The Laplacian is calculated as     *
!*  the divergence of a stress tensor, using the form suggested by     *
!*  Smagorinsky (1993).  The biharmonic is calculated by twice         *
!*  applying the divergence of the stress tensor that is used to       *
!*  calculate the Laplacian, but without the dependence on thickness   *
!*  in the first pass.  This form permits a variable viscosity, and    *
!*  indicates no acceleration for either resting fluid or solid body   *
!*  rotation.                                                          *
!*                                                                     *
!*    set_up_hor_visc calculates and stores the values of a number of  *
!*  metric functions that are used in horizontal_viscosity.  It is     *
!*  called by horizontal_viscosity the first time that the latter is   *
!*  called.                                                            *
!*                                                                     *
!*    The form of the Laplacian viscosity is:                          *
!*                                                                     *
!*    diffu = 1/h * {d/dx[KH*h*sh_xx] + d/dy[KH*h*sh_xy]}              *
!*    diffv = 1/h * {d/dx[KH*h*sh_xy] - d/dy[KH*h*sh_xx]}              *
!*                                                                     *
!*    sh_xx = du/dx - dv/dy         sh_xy = du/dy + dv/dx              *
!*                                                                     *
!*  with appropriate metric terms thrown in.  KH may either be a       *
!*  constant or may vary with the shear, as proposed by Smagorinsky.   *
!*  The form of this term is discussed extensively in Griffies and     *
!*  Hallberg (MWR, 2000), and the implementation here follows that     *
!*  discussion closely.                                                *
!*                                                                     *
!*    Only free slip boundary conditions have been coded, although     *
!*  no slip boundary conditions could be used with the Laplacian       *
!*  viscosity.  For a western boundary, for example, the boundary      *
!*  conditions with the biharmonic operator would be written as:       *
!*    dv/dx = 0, d^3v/dx^3 = 0, u = 0, d^2u/dx^2 = 0 ,                 *
!*  while for a Laplacian operator, they are simply:                   *
!*    dv/dx = 0, u = 0 .                                               *
!*  These boundary conditions are largely dictated by the use of       *
!*  a an Arakawa C-grid and by the varying layer thickness.            *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!* Macros written all in capital letters are defined in MOM_memory.h.  *
!*                                                                     *
!*     A small fragment of the C-grid is shown below:                  *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu, hq, str_xy, sh_xy         *
!*    j+1  > o > o >   At ^:  v, diffv, v0                             *
!*    j    x ^ x ^ x   At >:  u, diffu, u0                             *
!*    j    > o > o >   At o:  h, str_xx, sh_xx                         *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator,         only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_domains,               only : pass_var
use MOM_error_handler,         only : MOM_error, FATAL, WARNING
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types,            only : MEKE_type
use MOM_open_boundary,         only : ocean_OBC_type, OBC_NONE
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_io,                    only : read_data, slasher

implicit none ; private

#include <MOM_memory.h>

public horizontal_viscosity, hor_visc_init, hor_visc_end

type, public :: hor_visc_CS ; private
  logical :: Laplacian       ! Use a Laplacian horizontal viscosity if true.
  logical :: biharmonic      ! Use a biharmonic horizontal viscosity if true.
  logical :: no_slip         ! If true, no slip boundary conditions are used.
                             ! Otherwise free slip boundary conditions are assumed.
                             ! The implementation of the free slip boundary
                             ! conditions on a C-grid is much cleaner than the
                             ! no slip boundary conditions. The use of free slip
                             ! b.c.s is strongly encouraged. The no slip b.c.s
                             ! are not implemented with the biharmonic viscosity.
  logical :: bound_Kh        ! If true, the Laplacian coefficient is locally
                             ! limited to guarantee stability.
  logical :: better_bound_Kh ! If true, use a more careful bounding of the
                             ! Laplacian viscosity to guarantee stability.
  logical :: bound_Ah        ! If true, the biharmonic coefficient is locally
                             ! limited to guarantee stability.
  logical :: better_bound_Ah ! If true, use a more careful bounding of the
                             ! biharmonic viscosity to guarantee stability.
  real    :: bound_coef      ! The nondimensional coefficient of the ratio of
                             ! the viscosity bounds to the theoretical maximum
                             ! for stability without considering other terms.
                             ! The default is 0.8.
  logical :: Smagorinsky_Kh  ! If true, use Smagorinsky nonlinear eddy
                             ! viscosity. KH is the background value.
  logical :: Smagorinsky_Ah  ! If true, use a biharmonic form of Smagorinsky
                             ! nonlinear eddy viscosity. AH is the background.
  logical :: bound_Coriolis  ! If true & SMAGORINSKY_AH is used, the biharmonic
                             ! viscosity is modified to include a term that
                             ! scales quadratically with the velocity shears.
  logical :: use_Kh_bg_2d    ! Read 2d background viscosity from a file.
  real    :: Kh_bg_min       ! The minimum value allowed for Laplacian horizontal
                             ! viscosity. The default is 0.0

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    Kh_bg_xx,        &! The background Laplacian viscosity at h points, in units
                      ! of m2 s-1. The actual viscosity may be the larger of this
                      ! viscosity and the Smagorinsky viscosity.
    Kh_bg_2d,        &! The background Laplacian viscosity at h points, in units
                      ! of m2 s-1. The actual viscosity may be the larger of this
                      ! viscosity and the Smagorinsky viscosity.
    Ah_bg_xx,        &! The background biharmonic viscosity at h points, in units
                      ! of m4 s-1. The actual viscosity may be the larger of this
                      ! viscosity and the Smagorinsky viscosity.
    Kh_Max_xx,       &! The maximum permitted Laplacian viscosity, m2 s-1.
    Ah_Max_xx,       &! The maximum permitted biharmonic viscosity, m4 s-1.
    Biharm_Const2_xx,&! A constant relating the biharmonic viscosity to the
                      ! square of the velocity shear, in m4 s.  This value is
                      ! set to be the magnitude of the Coriolis terms once the
                      ! velocity differences reach a value of order 1/2 MAXVEL.

    reduction_xx      ! The amount by which stresses through h points are reduced
                      ! due to partial barriers. Nondimensional.

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    Kh_bg_xy,        &! The background Laplacian viscosity at q points, in units
                      ! of m2 s-1. The actual viscosity may be the larger of this
                      ! viscosity and the Smagorinsky viscosity.
    Ah_bg_xy,        &! The background biharmonic viscosity at q points, in units
                      ! of m4 s-1. The actual viscosity may be the larger of this
                      ! viscosity and the Smagorinsky viscosity.
    Kh_Max_xy,       &! The maximum permitted Laplacian viscosity, m2 s-1.
    Ah_Max_xy,       &! The maximum permitted biharmonic viscosity, m4 s-1.
    Biharm_Const2_xy,&! A constant relating the biharmonic viscosity to the
                      ! square of the velocity shear, in m4 s.  This value is
                      ! set to be the magnitude of the Coriolis terms once the
                      ! velocity differences reach a value of order 1/2 MAXVEL.
    reduction_xy      ! The amount by which stresses through q points are reduced
                      ! due to partial barriers. Nondimensional.

! The following variables are precalculated combinations of metric terms.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    dx2h, dy2h, &       ! dx^2  and dy^2  at h points, in m2
    dx_dyT, dy_dxT      ! dx/dy and dy/dx at h points, nondim
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    dx2q, dy2q, &       ! dx^2  and dy^2  at q points, in m2
    dx_dyBu, dy_dxBu    ! dx/dy and dy/dx at q points, nondim
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    Idx2dyCu, Idxdy2u   ! 1/(dx^2 dy) and 1/(dx dy^2) at u points, in m-3
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    Idx2dyCv, Idxdy2v   ! 1/(dx^2 dy) and 1/(dx dy^2) at v points, in m-3

! The following variables are precalculated time-invariant combinations of
! parameters and metric terms.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    Laplac_Const_xx, & ! Laplacian  metric-dependent constants (nondim)
    Biharm_Const_xx    ! Biharmonic metric-dependent constants (nondim)
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    Laplac_Const_xy, & ! Laplacian  metric-dependent constants (nondim)
    Biharm_Const_xy    ! Biharmonic metric-dependent constants (nondim)

  type(diag_ctrl), pointer :: diag ! structure to regulate diagnostic timing

  ! diagnostic ids
  integer :: id_diffu     = -1, id_diffv         = -1
  integer :: id_Ah_h      = -1, id_Ah_q          = -1
  integer :: id_Kh_h      = -1, id_Kh_q          = -1
  integer :: id_FrictWork = -1, id_FrictWorkIntz = -1

end type hor_visc_CS

contains

subroutine horizontal_viscosity(u, v, h, diffu, diffv, MEKE, VarMix, G, GV, CS, OBC)
  type(ocean_grid_type),                     intent(in)  :: G
  type(verticalGrid_type),                   intent(in)  :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)  :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)  :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)  :: h
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out) :: diffu
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out) :: diffv
  type(MEKE_type),                           pointer     :: MEKE
  type(VarMix_CS),                           pointer     :: VarMix
  type(hor_visc_CS),                         pointer     :: CS
  type(ocean_OBC_type),              pointer, optional   :: OBC

! Arguments:
!  (in)      u      - zonal velocity (m/s)
!  (in)      v      - meridional velocity (m/s)
!  (in)      h      - layer thickness (m or kg m-2); h units are referred to as H.
!  (out)     diffu  - zonal acceleration due to convergence of
!                     along-coordinate stress tensor (m/s2)
!  (out)     diffv  - meridional acceleration due to convergence of
!                     along-coordinate stress tensor (m/s2)
!  (inout)   MEKE   - pointer to a structure containing fields related to
!                     Mesoscale Eddy Kinetic Energy
!  (in)      VarMix - pointer to a structure with fields that specify the
!                     spatially variable viscosities
!  (in)      G      - ocean grid structure
!  (in)      GV     - The ocean's vertical grid structure.
!  (in)      CS     - control structure returned by a previous call to
!                     hor_visc_init
!  (in)      OBC    - pointer to an open boundary condition type

!  By R. Hallberg, August 1998 - November 1998.
!    This subroutine determines the acceleration due to the
!  horizontal viscosity.  A combination of biharmonic and Laplacian
!  forms can be used.  The coefficient may either be a constant or
!  a shear-dependent form.  The biharmonic is determined by twice
!  taking the divergence of an appropriately defined stress tensor.
!  The Laplacian is determined by doing so once.
!    To work, the following fields must be set outside of the usual
!  is to ie range before this subroutine is called:
!   v[is-2,is-1,ie+1,ie+2], u[is-2,is-1,ie+1,ie+2], and h[is-1,ie+1],
!  with a similarly sized halo in the y-direction.

  real, dimension(SZIB_(G),SZJ_(G)) :: u0 ! Laplacian of u (m-1 s-1)
  real, dimension(SZI_(G),SZJB_(G)) :: v0 ! Laplacian of v (m-1 s-1)

  real, dimension(SZI_(G),SZJ_(G)) :: &
    sh_xx, &      ! horizontal tension (du/dx - dv/dy) (1/sec) including metric terms
    str_xx,&      ! str_xx is the diagonal term in the stress tensor (H m2 s-2)
    bhstr_xx,&    ! A copy of str_xx that only contains the biharmonic contribution (H m2 s-2)
    FrictWorkIntz ! depth integrated energy dissipated by lateral friction (W/m2)

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dvdx, dudy, & ! components in the shearing strain (s-1)
    sh_xy,  &     ! horizontal shearing strain (du/dy + dv/dx) (1/sec) including metric terms
    str_xy, &     ! str_xy is the cross term in the stress tensor (H m2 s-2)
    bhstr_xy      ! A copy of str_xy that only contains the biharmonic contribution (H m2 s-2)

  real, dimension(SZIB_(G),SZJB_(G),SZK_(G)) :: &
    Ah_q, &   ! biharmonic viscosity at corner points (m4/s)
    Kh_q      ! Laplacian viscosity at corner points (m2/s)

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    Ah_h, &          ! biharmonic viscosity at thickness points (m4/s)
    Kh_h, &          ! Laplacian viscosity at thickness points (m2/s)
    FrictWork        ! energy dissipated by lateral friction (W/m2)

  real :: Ah         ! biharmonic viscosity (m4/s)
  real :: Kh         ! Laplacian  viscosity (m2/s)
  real :: AhSm       ! Smagorinsky biharmonic viscosity (m4/s)
  real :: KhSm       ! Smagorinsky Laplacian viscosity  (m2/s)
  real :: Shear_mag  ! magnitude of the shear (1/s)
  real :: huq, hvq   ! temporary variables in units of H^2 (i.e. m2 or kg2 m-4).
  real :: hu, hv     ! thicknesses at velocity points in units of H (i.e. m or kg m-2).
  real :: hq         ! harmonic mean of the harmonic means of the u- & v-
                     ! point thicknesses, in H; guarantees that hq/hu < 4.
  real :: h_neglect  ! thickness so small it can be lost in roundoff and so neglected (H)
  real :: h_neglect3 ! h_neglect^3, in H3
  real :: hrat_min   ! minimum thicknesses at the 4 neighboring
                     ! velocity points divided by the thickness at the stress
                     ! point (h or q point) (nondimensional)
  real :: visc_bound_rem ! fraction of overall viscous bounds that
                         ! remain to be applied (nondim)
  real :: Kh_scale  ! A factor between 0 and 1 by which the horizontal
                    ! Laplacian viscosity is rescaled
  real :: RoScl     ! The scaling function for MEKE source term
  real :: FatH      ! abs(f) at h-point for MEKE source term (s-1)

  logical :: rescale_Kh
  logical :: find_FrictWork
  logical :: apply_OBC = .false.
  logical :: use_MEKE_Ku
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  h_neglect  = GV%H_subroundoff
  h_neglect3 = h_neglect**3

  if (present(OBC)) then ; if (associated(OBC)) then ; if (OBC%OBC_pe) then
    apply_OBC = OBC%Flather_u_BCs_exist_globally .or. OBC%Flather_v_BCs_exist_globally
    apply_OBC = .true.
  endif ; endif ; endif

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "MOM_hor_visc: Module must be initialized before it is used.")

  find_FrictWork = (CS%id_FrictWork > 0)
  if (CS%id_FrictWorkIntz > 0)    find_FrictWork = .true.
  if (associated(MEKE)) then
    if (associated(MEKE%mom_src)) find_FrictWork = .true.
  endif

  rescale_Kh = .false.
  if (associated(VarMix)) then
    rescale_Kh = VarMix%Resoln_scaled_Kh
    if (rescale_Kh .and. &
    (.not.associated(VarMix%Res_fn_h) .or. .not.associated(VarMix%Res_fn_q))) &
      call MOM_error(FATAL, "MOM_hor_visc: VarMix%Res_fn_h and " //&
        "VarMix%Res_fn_q both need to be associated with Resoln_scaled_Kh.")
  endif

  ! Toggle whether to use a Laplacian viscosity derived from MEKE
  use_MEKE_Ku = associated(MEKE%Ku)

!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,nz,CS,G,GV,u,v,is,js,ie,je,h,  &
!$OMP                                  rescale_Kh,VarMix,h_neglect,h_neglect3,        &
!$OMP                                  Kh_h,Ah_h,Kh_q,Ah_q,diffu,apply_OBC,OBC,diffv, &
!$OMP                                  find_FrictWork,FrictWork,use_MEKE_Ku,MEKE)     &
!$OMP                          private(u0, v0, sh_xx, str_xx, visc_bound_rem,         &
!$OMP                                  sh_xy, str_xy, Ah, Kh, AhSm, KhSm,             &
!$OMP                                  bhstr_xx, bhstr_xy,FatH,RoScl, hu, hv,         &
!$OMP                                  Shear_mag, huq, hvq, hq, Kh_scale, hrat_min)
  do k=1,nz

!    This code uses boundary conditions that are consistent with
!  free slip and no normal flow boundary conditions.  The boundary
!  conditions for the western boundary, for example, are:
!    dv/dx = 0,  d^3v/dx^3 = 0,    u = 0,     d^2u/dx^2 = 0 .
!  The overall scheme is second order accurate.
!    All of the metric terms are retained, and the repeated use of
!  the symmetric stress tensor insures that no stress is applied with
!  no flow or solid-body rotation, even with non-constant values of
!  of the biharmonic viscosity.

!  The following are the forms of the horizontal tension and hori-
!  shearing strain advocated by Smagorinsky (1993) and discussed in
!  Griffies and Hallberg (MWR, 2000).
    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      sh_xx(i,j) = (CS%DY_dxT(i,j)*(G%IdyCu(I,j) * u(I,j,k) - &
                                    G%IdyCu(I-1,j) * u(I-1,j,k)) - &
                    CS%DX_dyT(i,j)*(G%IdxCv(i,J) * v(i,J,k) - &
                                    G%IdxCv(i,J-1)*v(i,J-1,k)))
    enddo ; enddo
    ! Components for the shearing strain
    do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
      dvdx(I,J) = CS%DY_dxBu(I,J)*(v(i+1,J,k)*G%IdyCv(i+1,J) - v(i,J,k)*G%IdyCv(i,J))
      dudy(I,J) = CS%DX_dyBu(I,J)*(u(I,j+1,k)*G%IdxCu(I,j+1) - u(I,j,k)*G%IdxCu(I,j))
    enddo ; enddo
    ! Adjust contributions to shearing strain on open boundaries.
    if (apply_OBC) then ; if (OBC%zero_strain .or. OBC%freeslip_strain) then
      do n=1,OBC%number_of_segments
        if (OBC%segment(n)%is_N_or_S) then
          J = OBC%segment(n)%HI%JsdB
          do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
            if (OBC%zero_strain) then
              dvdx(I,J) = 0.
              dudy(I,J) = 0.
            elseif (OBC%freeslip_strain) then
              dudy(I,J) = 0.
            endif
          enddo
        elseif (OBC%segment(n)%is_E_or_W) then
          I = OBC%segment(n)%HI%IsdB
          do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
            if (OBC%zero_strain) then
              dvdx(I,J) = 0.
              dudy(I,J) = 0.
            elseif (OBC%freeslip_strain) then
              dvdx(I,J) = 0.
            endif
          enddo
        endif
      enddo
    endif ; endif
    if (CS%no_slip) then
      do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
        sh_xy(I,J) = (2.0-G%mask2dBu(I,J)) * ( dvdx(I,J) + dudy(I,J) )
      enddo ; enddo
    else
      do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
        sh_xy(I,J) = G%mask2dBu(I,J) * ( dvdx(I,J) + dudy(I,J) )
      enddo ; enddo
    endif

!  Evaluate u0 = x.Div(Grad u) and v0 = y.Div( Grad u)
    if (CS%biharmonic) then
      do j=js-1,Jeq+1 ; do I=Isq-1,Ieq+1
        u0(I,j) = CS%IDXDY2u(I,j)*(CS%DY2h(i+1,j)*sh_xx(i+1,j) - CS%DY2h(i,j)*sh_xx(i,j)) + &
                  CS%IDX2dyCu(I,j)*(CS%DX2q(I,J)*sh_xy(I,J) - CS%DX2q(I,J-1)*sh_xy(I,J-1))
      enddo ; enddo
      do J=Jsq-1,Jeq+1 ; do i=is-1,Ieq+1
        v0(i,J) = CS%IDXDY2v(i,J)*(CS%DY2q(I,J)*sh_xy(I,J) - CS%DY2q(I-1,J)*sh_xy(I-1,J)) - &
                  CS%IDX2dyCv(i,J)*(CS%DX2h(i,j+1)*sh_xx(i,j+1) - CS%DX2h(i,j)*sh_xx(i,j))
      enddo ; enddo
      if (apply_OBC .and. OBC%zero_biharmonic) then
        do n=1,OBC%number_of_segments
          if (OBC%segment(n)%is_N_or_S) then
            J = OBC%segment(n)%HI%JsdB
            do I=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
              v0(i,J) = 0.
            enddo
          elseif (OBC%segment(n)%is_E_or_W) then
            I = OBC%segment(n)%HI%IsdB
            do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
              u0(I,j) = 0.
            enddo
          endif
        enddo
      endif
    endif

    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      if ((CS%Smagorinsky_Kh) .or. (CS%Smagorinsky_Ah)) &
        Shear_mag = sqrt(sh_xx(i,j)*sh_xx(i,j) + &
          0.25*((sh_xy(I-1,J-1)*sh_xy(I-1,J-1) + sh_xy(I,J)*sh_xy(I,J)) + &
                (sh_xy(I-1,J)*sh_xy(I-1,J) + sh_xy(I,J-1)*sh_xy(I,J-1))))
      if (CS%better_bound_Ah .or. CS%better_bound_Kh) then
        hrat_min = min(1.0, &
            0.5*min((h(i+1,j,k) + h(i,j,k)), (h(i,j,k) + h(i-1,j,k)), &
                    (h(i,j+1,k) + h(i,j,k)), (h(i,j-1,k) + h(i,j,k))) / &
            (h(i,j,k) + h_neglect) )
        visc_bound_rem = 1.0
      endif

      if (CS%Laplacian) then
        ! Determine the Laplacian viscosity at h points, using the
        ! largest value from several parameterizations.
        Kh_scale = 1.0 ; if (rescale_Kh) Kh_scale = VarMix%Res_fn_h(i,j)
        if (CS%Smagorinsky_Kh) then
          KhSm = CS%LAPLAC_CONST_xx(i,j) * Shear_mag
          Kh = Kh_scale * MAX(CS%Kh_bg_xx(i,j), KhSm)
          if (CS%bound_Kh .and. .not.CS%better_bound_Kh) &
            Kh = MIN(Kh, CS%Kh_Max_xx(i,j))
        else
          Kh = Kh_scale * CS%Kh_bg_xx(i,j)
        endif
        if (use_MEKE_Ku) then
          Kh = Kh + MEKE%Ku(i,j)
        endif
        if (CS%better_bound_Kh) then
          if (Kh >= hrat_min*CS%Kh_Max_xx(i,j)) then
            visc_bound_rem = 0.0
            Kh = hrat_min*CS%Kh_Max_xx(i,j)
          else
           !visc_bound_rem = 1.0 - abs(Kh) / (hrat_min*CS%Kh_Max_xx(i,j))
            visc_bound_rem = 1.0 - Kh / (hrat_min*CS%Kh_Max_xx(i,j))
          endif
        endif

        if (CS%id_Kh_h>0) Kh_h(i,j,k) = Kh

        str_xx(i,j) = -Kh * sh_xx(i,j)
      else   ! not Laplacian
        str_xx(i,j) = 0.0
      endif ! Laplacian

      if (CS%biharmonic) then
!       Determine the biharmonic viscosity at h points, using the
!       largest value from several parameterizations.
        if (CS%Smagorinsky_Ah) then
          if (CS%bound_Coriolis) then
            AhSm =  Shear_mag * (CS%BIHARM_CONST_xx(i,j) + &
                                 CS%Biharm_Const2_xx(i,j)*Shear_mag)
          else
            AhSm = CS%BIHARM_CONST_xx(i,j) * Shear_mag
          endif
          Ah = MAX(CS%Ah_bg_xx(i,j), AhSm)
          if (CS%bound_Ah .and. .not.CS%better_bound_Ah) &
            Ah = MIN(Ah, CS%Ah_Max_xx(i,j))
        else
          Ah = CS%Ah_bg_xx(i,j)
        endif ! Smagorinsky_Ah
        if (CS%better_bound_Ah) then
          Ah = MIN(Ah, visc_bound_rem*hrat_min*CS%Ah_Max_xx(i,j))
        endif

        if (CS%id_Ah_h>0) Ah_h(i,j,k) = Ah

        str_xx(i,j) = str_xx(i,j) + Ah * &
          (CS%DY_dxT(i,j)*(G%IdyCu(I,j)*u0(I,j) - G%IdyCu(I-1,j)*u0(I-1,j)) - &
           CS%DX_dyT(i,j) *(G%IdxCv(i,J)*v0(i,J) - G%IdxCv(i,J-1)*v0(i,J-1)))

        ! Keep a copy of the biharmonic contribution for backscatter parameterization
        bhstr_xx(i,j) =             Ah * &
          (CS%DY_dxT(i,j)*(G%IdyCu(I,j)*u0(I,j) - G%IdyCu(I-1,j)*u0(I-1,j)) - &
           CS%DX_dyT(i,j) *(G%IdxCv(i,J)*v0(i,J) - G%IdxCv(i,J-1)*v0(i,J-1)))
        bhstr_xx(i,j) = bhstr_xx(i,j) * (h(i,j,k) * CS%reduction_xx(i,j))

      endif  ! biharmonic

      str_xx(i,j) = str_xx(i,j) * (h(i,j,k) * CS%reduction_xx(i,j))
    enddo ; enddo

    if (CS%biharmonic) then
      ! Gradient of Laplacian, for use in bi-harmonic term
      do J=js-1,Jeq ; do I=is-1,Ieq
        dvdx(I,J) = CS%DY_dxBu(I,J)*(v0(i+1,J)*G%IdyCv(i+1,J) - v0(i,J)*G%IdyCv(i,J))
        dudy(I,J) = CS%DX_dyBu(I,J)*(u0(I,j+1)*G%IdxCu(I,j+1) - u0(I,j)*G%IdxCu(I,j))
      enddo ; enddo
      ! Adjust contributions to shearing strain on open boundaries.
      if (apply_OBC) then ; if (OBC%zero_strain .or. OBC%freeslip_strain) then
        do n=1,OBC%number_of_segments
          if (OBC%segment(n)%is_N_or_S) then
            J = OBC%segment(n)%HI%JsdB
            do I=OBC%segment(n)%HI%IsdB,OBC%segment(n)%HI%IedB
              if (OBC%zero_strain) then
                dvdx(I,J) = 0.
                dudy(I,J) = 0.
              elseif (OBC%freeslip_strain) then
                dudy(I,J) = 0.
              endif
            enddo
          elseif (OBC%segment(n)%is_E_or_W) then
            I = OBC%segment(n)%HI%IsdB
            do J=OBC%segment(n)%HI%JsdB,OBC%segment(n)%HI%JedB
              if (OBC%zero_strain) then
                dvdx(I,J) = 0.
                dudy(I,J) = 0.
              elseif (OBC%freeslip_strain) then
                dvdx(I,J) = 0.
              endif
            enddo
          endif
        enddo
      endif ; endif
    endif
    do J=js-1,Jeq ; do I=is-1,Ieq
      if ((CS%Smagorinsky_Kh) .or. (CS%Smagorinsky_Ah)) &
        Shear_mag = sqrt(sh_xy(I,J)*sh_xy(I,J) + &
            0.25*((sh_xx(i,j)*sh_xx(i,j) + sh_xx(i+1,j+1)*sh_xx(i+1,j+1)) + &
                  (sh_xx(i,j+1)*sh_xx(i,j+1) + sh_xx(i+1,j)*sh_xx(i+1,j))))

      huq = (h(i,j,k) + h(i+1,j,k)) * (h(i,j+1,k) + h(i+1,j+1,k))
      hvq = (h(i,j,k) + h(i,j+1,k)) * (h(i+1,j,k) + h(i+1,j+1,k))
      hq = 2.0 * huq * hvq / (h_neglect3 + (huq + hvq) * &
          ((h(i,j,k) + h(i+1,j+1,k)) + (h(i,j+1,k) + h(i+1,j,k))))

      if (CS%better_bound_Ah .or. CS%better_bound_Kh) then
        hrat_min = min(1.0, &
            0.5*min((h(i,j,k) + h(i+1,j,k)), (h(i,j+1,k) + h(i+1,j+1,k)), &
                    (h(i,j,k) + h(i,j+1,k)), (h(i+1,j,k) + h(i+1,j+1,k))) / &
                   (hq + h_neglect) )
        visc_bound_rem = 1.0
      endif

      if (CS%no_slip .and. (G%mask2dBu(I,J) < 0.5)) then
        if ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) + &
            (G%mask2dCv(i,J) + G%mask2dCv(i+1,J)) > 0.0) then
          ! This is a coastal vorticity point, so modify hq and hrat_min.

          hu = 0.5 * (G%mask2dCu(I,j)   * (h(i,j,k) + h(i+1,j,k)) + &
                      G%mask2dCu(I,j+1) * (h(i,j+1,k) + h(i+1,j+1,k)))
          hv = 0.5 * (G%mask2dCv(i,J)   * (h(i,j,k) + h(i,j+1,k)) + &
                      G%mask2dCv(i+1,J) * (h(i+1,j,k) + h(i+1,j+1,k)))
          if ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) * &
              (G%mask2dCv(i,J) + G%mask2dCv(i+1,J)) == 0.0) then
            ! Only one of hu and hv is nonzero, so just add them.
            hq = hu + hv
            hrat_min = 1.0
          else
            ! Both hu and hv are nonzero, so take the harmonic mean.
            hq = 2.0 * (hu * hv) / ((hu + hv) + h_neglect)
            hrat_min = min(1.0, min(hu, hv) / (hq + h_neglect) )
          endif
        endif
      endif

      if (CS%Laplacian) then
        ! Determine the Laplacian viscosity at q points, using the
        ! largest value from several parameterizations.
        Kh_scale = 1.0 ; if (rescale_Kh) Kh_scale = VarMix%Res_fn_q(I,J)
        if (CS%Smagorinsky_Kh) then
          KhSm = CS%LAPLAC_CONST_xy(I,J) * Shear_mag
          Kh = Kh_scale * MAX(CS%Kh_bg_xy(I,J), KhSm)
          if (CS%bound_Kh .and. .not.CS%better_bound_Kh) &
            Kh = MIN(Kh, CS%Kh_Max_xy(I,J))
        else
          Kh = Kh_scale * CS%Kh_bg_xy(I,J)
        endif
        if (use_MEKE_Ku) then
          Kh = Kh + 0.25*( (MEKE%Ku(I,J)+MEKE%Ku(I+1,J+1))    &
                          +(MEKE%Ku(I+1,J)+MEKE%Ku(I,J+1)) )
        endif
        ! Place a floor on the viscosity, if desired.
        Kh = MAX(Kh,CS%Kh_bg_min)

        if (CS%better_bound_Kh) then
          if (Kh >= hrat_min*CS%Kh_Max_xy(I,J)) then
            visc_bound_rem = 0.0
            Kh = hrat_min*CS%Kh_Max_xy(I,J)
          elseif (CS%Kh_Max_xy(I,J)>0.) then
           !visc_bound_rem = 1.0 - abs(Kh) / (hrat_min*CS%Kh_Max_xy(I,J))
            visc_bound_rem = 1.0 - Kh / (hrat_min*CS%Kh_Max_xy(I,J))
          endif
        endif

        if (CS%id_Kh_q>0) Kh_q(I,J,k) = Kh

        str_xy(I,J) = -Kh * sh_xy(I,J)
      else   ! not Laplacian
        str_xy(I,J) = 0.0
      endif ! Laplacian

      if (CS%biharmonic) then
      ! Determine the biharmonic viscosity at q points, using the
      ! largest value from several parameterizations.
        if (CS%Smagorinsky_Ah) then
          if (CS%bound_Coriolis) then
            AhSm =  Shear_mag * (CS%BIHARM_CONST_xy(I,J) + &
                                 CS%Biharm_Const2_xy(I,J)*Shear_mag)
          else
            AhSm = CS%BIHARM_CONST_xy(I,J) * Shear_mag
          endif
          Ah = MAX(CS%Ah_bg_xy(I,J), AhSm)
          if (CS%bound_Ah .and. .not.CS%better_bound_Ah) &
            Ah = MIN(Ah, CS%Ah_Max_xy(I,J))
        else
          Ah = CS%Ah_bg_xy(I,J)
        endif ! Smagorinsky_Ah
        if (CS%better_bound_Ah) then
          Ah = MIN(Ah, visc_bound_rem*hrat_min*CS%Ah_Max_xy(I,J))
        endif

        if (CS%id_Ah_q>0) Ah_q(I,J,k) = Ah

        str_xy(I,J) = str_xy(I,J) + Ah * ( dvdx(I,J) + dudy(I,J) )

        ! Keep a copy of the biharmonic contribution for backscatter parameterization
        bhstr_xy(I,J) = Ah * ( dvdx(I,J) + dudy(I,J) ) * &
                        (hq * G%mask2dBu(I,J) * CS%reduction_xy(I,J))

      endif  ! biharmonic

      if (CS%no_slip) then
        str_xy(I,J) = str_xy(I,J) * (hq * CS%reduction_xy(I,J))
      else
        str_xy(I,J) = str_xy(I,J) * (hq * G%mask2dBu(I,J) * CS%reduction_xy(I,J))
      endif
    enddo ; enddo

    do j=js,je ; do i=isq,ieq
!  Evaluate 1/h x.Div(h Grad u) or the biharmonic equivalent.
      diffu(I,j,k) = ((G%IdyCu(I,j)*(CS%DY2h(i,j) *str_xx(i,j) - &
                                    CS%DY2h(i+1,j)*str_xx(i+1,j)) + &
                       G%IdxCu(I,j)*(CS%DX2q(I,J-1)*str_xy(I,J-1) - &
                                    CS%DX2q(I,J) *str_xy(I,J))) * &
                     G%IareaCu(I,j)) / (0.5*(h(i+1,j,k) + h(i,j,k)) + h_neglect)

    enddo ; enddo
    if (apply_OBC) then
      ! This is not the right boundary condition. If all the masking of tendencies are done
      ! correctly later then eliminating this block should not change answers.
      do n=1,OBC%number_of_segments
        if (OBC%segment(n)%is_E_or_W) then
          I = OBC%segment(n)%HI%IsdB
          do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
            diffu(I,j,k) = 0.
          enddo
        endif
      enddo
    endif

!  Evaluate 1/h y.Div(h Grad u) or the biharmonic equivalent.
    do J=Jsq,Jeq ; do i=is,ie
      diffv(i,J,k) = ((G%IdyCv(i,J)*(CS%DY2q(I-1,J)*str_xy(I-1,J) - &
                                    CS%DY2q(I,J) *str_xy(I,J)) - &
                       G%IdxCv(i,J)*(CS%DX2h(i,j) *str_xx(i,j) - &
                                    CS%DX2h(i,j+1)*str_xx(i,j+1))) * &
                     G%IareaCv(i,J)) / (0.5*(h(i,j+1,k) + h(i,j,k)) + h_neglect)
    enddo ; enddo
    if (apply_OBC) then
      ! This is not the right boundary condition. If all the masking of tendencies are done
      ! correctly later then eliminating this block should not change answers.
      do n=1,OBC%number_of_segments
        if (OBC%segment(n)%is_N_or_S) then
          J = OBC%segment(n)%HI%JsdB
          do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
            diffv(i,J,k) = 0.
          enddo
        endif
      enddo
    endif

    if (find_FrictWork) then ; do j=js,je ; do i=is,ie
    ! Diagnose   str_xx*d_x u - str_yy*d_y v + str_xy*(d_y u + d_x v)
      FrictWork(i,j,k) = GV%H_to_kg_m2 * ( &
              (str_xx(i,j)*(u(I,j,k)-u(I-1,j,k))*G%IdxT(i,j)     &
              -str_xx(i,j)*(v(i,J,k)-v(i,J-1,k))*G%IdyT(i,j))    &
       +0.25*((str_xy(I,J)*(                                     &
                   (u(I,j+1,k)-u(I,j,k))*G%IdyBu(I,J)            &
                  +(v(i+1,J,k)-v(i,J,k))*G%IdxBu(I,J) )          &
              +str_xy(I-1,J-1)*(                                 &
                   (u(I-1,j,k)-u(I-1,j-1,k))*G%IdyBu(I-1,J-1)    &
                  +(v(i,J-1,k)-v(i-1,J-1,k))*G%IdxBu(I-1,J-1) )) &
             +(str_xy(I-1,J)*(                                   &
                   (u(I-1,j+1,k)-u(I-1,j,k))*G%IdyBu(I-1,J)      &
                  +(v(i,J,k)-v(i-1,J,k))*G%IdxBu(I-1,J) )        &
              +str_xy(I,J-1)*(                                   &
                   (u(I,j,k)-u(I,j-1,k))*G%IdyBu(I,J-1)          &
                  +(v(i+1,J-1,k)-v(i,J-1,k))*G%IdxBu(I,J-1) )) ) )
    enddo ; enddo ; endif

    ! Make a similar calculation as for FrictWork above but accumulating into
    ! the vertically integrated MEKE source term, and adjusting for any
    ! energy loss seen as a reduction in the [biharmonic] frictional source term.
    if (find_FrictWork .and. associated(MEKE)) then ; if (associated(MEKE%mom_src)) then
      if (k==1) then
        do j=js,je ; do i=is,ie
          MEKE%mom_src(i,j) = 0.
        enddo ; enddo
      endif
      if (MEKE%backscatter_Ro_c /= 0.) then
        do j=js,je ; do i=is,ie
          FatH = 0.25*( (abs(G%CoriolisBu(I-1,J-1)) + abs(G%CoriolisBu(I,J))) &
                       +(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J-1))) )
          Shear_mag = sqrt(sh_xx(i,j)*sh_xx(i,j) + &
            0.25*((sh_xy(I-1,J-1)*sh_xy(I-1,J-1) + sh_xy(I,J)*sh_xy(I,J)) + &
                  (sh_xy(I-1,J)*sh_xy(I-1,J) + sh_xy(I,J-1)*sh_xy(I,J-1))))
          FatH = FatH ** MEKE%backscatter_Ro_pow ! f^n
          Shear_mag = ( ( Shear_mag ** MEKE%backscatter_Ro_pow ) + 1.e-30 ) &
                      * MEKE%backscatter_Ro_c ! c * D^n
          ! The Rossby number function is g(Ro) = 1/(1+c.Ro^n)
          ! RoScl = 1 - g(Ro)
          RoScl = Shear_mag / ( FatH + Shear_mag ) ! = 1 - f^n/(f^n+c*D^n)
          MEKE%mom_src(i,j) = MEKE%mom_src(i,j) + GV%H_to_kg_m2 * (                   &
                ((str_xx(i,j)-RoScl*bhstr_xx(i,j))*(u(I,j,k)-u(I-1,j,k))*G%IdxT(i,j)  &
                -(str_xx(i,j)-RoScl*bhstr_xx(i,j))*(v(i,J,k)-v(i,J-1,k))*G%IdyT(i,j)) &
         +0.25*(((str_xy(I,J)-RoScl*bhstr_xy(I,J))*(                                  &
                     (u(I,j+1,k)-u(I,j,k))*G%IdyBu(I,J)                               &
                    +(v(i+1,J,k)-v(i,J,k))*G%IdxBu(I,J) )                             &
                +(str_xy(I-1,J-1)-RoScl*bhstr_xy(I-1,J-1))*(                          &
                     (u(I-1,j,k)-u(I-1,j-1,k))*G%IdyBu(I-1,J-1)                       &
                    +(v(i,J-1,k)-v(i-1,J-1,k))*G%IdxBu(I-1,J-1) ))                    &
               +((str_xy(I-1,J)-RoScl*bhstr_xy(I-1,J))*(                              &
                     (u(I-1,j+1,k)-u(I-1,j,k))*G%IdyBu(I-1,J)                         &
                    +(v(i,J,k)-v(i-1,J,k))*G%IdxBu(I-1,J) )                           &
                +(str_xy(I,J-1)-RoScl*bhstr_xy(I,J-1))*(                              &
                     (u(I,j,k)-u(I,j-1,k))*G%IdyBu(I,J-1)                             &
                    +(v(i+1,J-1,k)-v(i,J-1,k))*G%IdxBu(I,J-1) )) ) )
        enddo ; enddo
      else
        do j=js,je ; do i=is,ie
         MEKE%mom_src(i,j) = MEKE%mom_src(i,j) + FrictWork(i,j,k)
        enddo ; enddo
      endif
    endif ; endif

  enddo ! end of k loop

! Offer fields for diagnostic averaging.
  if (CS%id_diffu>0)     call post_data(CS%id_diffu, diffu, CS%diag)
  if (CS%id_diffv>0)     call post_data(CS%id_diffv, diffv, CS%diag)
  if (CS%id_FrictWork>0) call post_data(CS%id_FrictWork, FrictWork, CS%diag)
  if (CS%id_Ah_h>0)      call post_data(CS%id_Ah_h, Ah_h, CS%diag)
  if (CS%id_Ah_q>0)      call post_data(CS%id_Ah_q, Ah_q, CS%diag)
  if (CS%id_Kh_h>0)      call post_data(CS%id_Kh_h, Kh_h, CS%diag)
  if (CS%id_Kh_q>0)      call post_data(CS%id_Kh_q, Kh_q, CS%diag)

  if (CS%id_FrictWorkIntz > 0) then
    do j=js,je
      do i=is,ie ; FrictWorkIntz(i,j) = FrictWork(i,j,1) ; enddo
      do k=2,nz ; do i=is,ie
        FrictWorkIntz(i,j) = FrictWorkIntz(i,j) + FrictWork(i,j,k)
      enddo ; enddo
    enddo
    call post_data(CS%id_FrictWorkIntz, FrictWorkIntz, CS%diag)
  endif


end subroutine horizontal_viscosity


subroutine hor_visc_init(Time, G, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(inout) :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(hor_visc_CS), pointer             :: CS

! This subroutine allocates space for and calculates static variables
! used by this module. The metrics may be 0, 1, or 2-D arrays,
! while fields like the background viscosities are 2-D arrays.
! ALLOC is a macro defined in MOM_memory.h to either allocate
! for dynamic memory, or do nothing when using static memory.
!
! Arguments:
!  (in)      Time       - current model time
!  (in)      G          - ocean grid structure
!  (in)      param_file - structure to parse for model parameter values
!  (in)      diag       - structure to regulate diagnostic output
!  (in/out)  CS         - pointer to the control structure for this module

  real, dimension(SZIB_(G),SZJ_(G)) :: u0u, u0v
  real, dimension(SZI_(G),SZJB_(G)) :: v0u, v0v
                ! u0v is the Laplacian sensitivities to the v velocities
                ! at u points, in m-2, with u0u, v0u, and v0v defined similarly.
  real :: grid_sp_h2       ! Harmonic mean of the squares of the grid
  real :: grid_sp_q2       ! spacings at h and q points (m2)
  real :: Kh_Limit         ! A coefficient (1/s) used, along with the
                           ! grid spacing, to limit Laplacian viscosity.
  real :: fmax             ! maximum absolute value of f at the four
                           ! vorticity points around a thickness point (1/s)
  real :: BoundCorConst    ! constant (s2/m2)
  real :: Ah_Limit         ! coefficient (1/s) used, along with the
                           ! grid spacing, to limit biharmonic viscosity
  real :: Kh               ! Lapacian horizontal viscosity (m2/s)
  real :: Ah               ! biharmonic horizontal viscosity (m4/s)
  real :: Kh_vel_scale     ! this speed (m/s) times grid spacing gives Lap visc
  real :: Ah_vel_scale     ! this speed (m/s) times grid spacing cubed gives bih visc
  real :: Smag_Lap_const   ! nondimensional Laplacian Smagorinsky constant
  real :: Smag_bi_const    ! nondimensional biharmonic Smagorinsky constant
  real :: dt               ! dynamics time step (sec)
  real :: Idt              ! inverse of dt (1/s)
  real :: denom            ! work variable; the denominator of a fraction
  real :: maxvel           ! largest permitted velocity components (m/s)
  real :: bound_Cor_vel    ! grid-scale velocity variations at which value
                           ! the quadratically varying biharmonic viscosity
                           ! balances Coriolis acceleration (m/s)
  logical :: bound_Cor_def ! parameter setting of BOUND_CORIOLIS
  logical :: get_all       ! If true, read and log all parameters, regardless of
                           ! whether they are used, to enable spell-checking of
                           ! valid parameters.
  character(len=64) :: inputdir, filename
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: i, j

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_hor_visc"  ! module name

  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = G%ke
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "hor_visc_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  ! Read parameters and write them to the model log.
  call log_version(param_file, mod, version, "")

  !   It is not clear whether these initialization lines are needed for the
  ! cases where the corresponding parameters are not read.
  CS%bound_Kh = .false. ; CS%better_bound_Kh = .false. ; CS%Smagorinsky_Kh = .false.
  CS%bound_Ah = .false. ; CS%better_bound_Ah = .false. ; CS%Smagorinsky_Ah = .false.
  CS%bound_Coriolis = .false.

  Kh = 0.0 ; Ah = 0.0

  !   If GET_ALL_PARAMS is true, all parameters are read in all cases to enable
  ! parameter spelling checks.
  call get_param(param_file, mod, "GET_ALL_PARAMS", get_all, default=.false.)

  call get_param(param_file, mod, "LAPLACIAN", CS%Laplacian, &
                 "If true, use a Laplacian horizontal viscosity.", &
                 default=.false.)
  if (CS%Laplacian .or. get_all) then
    call get_param(param_file, mod, "KH", Kh,                      &
                 "The background Laplacian horizontal viscosity.", &
                 units = "m2 s-1", default=0.0)
    call get_param(param_file, mod, "KH_BG_MIN", CS%Kh_bg_min, &
                 "The minimum value allowed for Laplacian horizontal viscosity, KH.", &
                 units = "m2 s-1",  default=0.0)
    call get_param(param_file, mod, "KH_VEL_SCALE", Kh_vel_scale, &
                 "The velocity scale which is multiplied by the grid \n"//&
                 "spacing to calculate the Laplacian viscosity. \n"//&
                 "The final viscosity is the largest of this scaled \n"//&
                 "viscosity, the Smagorinsky viscosity and KH.", &
                 units="m s-1", default=0.0)

    call get_param(param_file, mod, "SMAGORINSKY_KH", CS%Smagorinsky_Kh, &
                 "If true, use a Smagorinsky nonlinear eddy viscosity.", &
                 default=.false.)
    if (CS%Smagorinsky_Kh .or. get_all) &
      call get_param(param_file, mod, "SMAG_LAP_CONST", Smag_Lap_const, &
                 "The nondimensional Laplacian Smagorinsky constant, \n"//&
                 "often 0.15.", units="nondim", default=0.0, &
                  fail_if_missing = CS%Smagorinsky_Kh)

    call get_param(param_file, mod, "BOUND_KH", CS%bound_Kh, &
                 "If true, the Laplacian coefficient is locally limited \n"//&
                 "to be stable.", default=.true.)
    call get_param(param_file, mod, "BETTER_BOUND_KH", CS%better_bound_Kh, &
                 "If true, the Laplacian coefficient is locally limited \n"//&
                 "to be stable with a better bounding than just BOUND_KH.", &
                 default=CS%bound_Kh)
  endif

  call get_param(param_file, mod, "BIHARMONIC", CS%biharmonic, &
                 "If true, use a biharmonic horizontal viscosity. \n"//&
                 "BIHARMONIC may be used with LAPLACIAN.", &
                 default=.true.)
  if (CS%biharmonic .or. get_all) then
    call get_param(param_file, mod, "AH", Ah, &
                 "The background biharmonic horizontal viscosity.", &
                 units = "m4 s-1", default=0.0)
    call get_param(param_file, mod, "AH_VEL_SCALE", Ah_vel_scale, &
                 "The velocity scale which is multiplied by the cube of \n"//&
                 "the grid spacing to calculate the biharmonic viscosity. \n"//&
                 "The final viscosity is the largest of this scaled \n"//&
                 "viscosity, the Smagorinsky viscosity and AH.", &
                 units="m s-1", default=0.0)
    call get_param(param_file, mod, "SMAGORINSKY_AH", CS%Smagorinsky_Ah, &
                 "If true, use a biharmonic Smagorinsky nonlinear eddy \n"//&
                 "viscosity.", default=.false.)

    call get_param(param_file, mod, "BOUND_AH", CS%bound_Ah, &
                 "If true, the biharmonic coefficient is locally limited \n"//&
                 "to be stable.", default=.true.)
    call get_param(param_file, mod, "BETTER_BOUND_AH", CS%better_bound_Ah, &
                 "If true, the biharmonic coefficient is locally limited \n"//&
                 "to be stable with a better bounding than just BOUND_AH.", &
                 default=CS%bound_Ah)

    if (CS%Smagorinsky_Ah .or. get_all) then
      call get_param(param_file, mod, "SMAG_BI_CONST",Smag_bi_const, &
                 "The nondimensional biharmonic Smagorinsky constant, \n"//&
                 "typically 0.015 - 0.06.", units="nondim", default=0.0, &
                 fail_if_missing = CS%Smagorinsky_Ah)

      call get_param(param_file, mod, "BOUND_CORIOLIS", bound_Cor_def, default=.false.)
      call get_param(param_file, mod, "BOUND_CORIOLIS_BIHARM", CS%bound_Coriolis, &
                 "If true use a viscosity that increases with the square \n"//&
                 "of the velocity shears, so that the resulting viscous \n"//&
                 "drag is of comparable magnitude to the Coriolis terms \n"//&
                 "when the velocity differences between adjacent grid \n"//&
                 "points is 0.5*BOUND_CORIOLIS_VEL.  The default is the \n"//&
                 "value of BOUND_CORIOLIS (or false).", default=bound_Cor_def)
      if (CS%bound_Coriolis .or. get_all) then
        call get_param(param_file, mod, "MAXVEL", maxvel, default=3.0e8)
        bound_Cor_vel = maxvel
        call get_param(param_file, mod, "BOUND_CORIOLIS_VEL", bound_Cor_vel, &
                 "The velocity scale at which BOUND_CORIOLIS_BIHARM causes \n"//&
                 "the biharmonic drag to have comparable magnitude to the \n"//&
                 "Coriolis acceleration.  The default is set by MAXVEL.", &
                 units="m s-1", default=maxvel)
      endif
    endif
  endif

  if (CS%better_bound_Ah .or. CS%better_bound_Kh .or. get_all) &
    call get_param(param_file, mod, "HORVISC_BOUND_COEF", CS%bound_coef, &
                 "The nondimensional coefficient of the ratio of the \n"//&
                 "viscosity bounds to the theoretical maximum for \n"//&
                 "stability without considering other terms.", units="nondim", &
                 default=0.8)

  call get_param(param_file, mod, "NOSLIP", CS%no_slip, &
                 "If true, no slip boundary conditions are used; otherwise \n"//&
                 "free slip boundary conditions are assumed. The \n"//&
                 "implementation of the free slip BCs on a C-grid is much \n"//&
                 "cleaner than the no slip BCs. The use of free slip BCs \n"//&
                 "is strongly encouraged, and no slip BCs are not used with \n"//&
                 "the biharmonic viscosity.", default=.false.)

  call get_param(param_file, mod, "USE_KH_BG_2D", CS%use_Kh_bg_2d, &
                 "If true, read a file containing 2-d background harmonic  \n"//&
                 "viscosities. The final viscosity is the maximum of the other "//&
                 "terms and this background value.", default=.false.)


  if (CS%bound_Kh .or. CS%bound_Ah .or. CS%better_bound_Kh .or. CS%better_bound_Ah) &
    call get_param(param_file, mod, "DT", dt, &
                 "The (baroclinic) dynamics time step.", units = "s", &
                 fail_if_missing=.true.)

  if (CS%no_slip .and. CS%biharmonic) &
    call MOM_error(FATAL,"ERROR: NOSLIP and BIHARMONIC cannot be defined "// &
                          "at the same time in MOM.")

  if (.not.(CS%Laplacian .or. CS%biharmonic)) call MOM_error(WARNING, &
    "hor_visc_init:  It is usually a very bad idea not to use either "//&
    "LAPLACIAN or BIHARMONIC viscosity.")

  ALLOC_(CS%dx2h(isd:ied,jsd:jed))        ; CS%dx2h(:,:)    = 0.0
  ALLOC_(CS%dy2h(isd:ied,jsd:jed))        ; CS%dy2h(:,:)    = 0.0
  ALLOC_(CS%dx2q(IsdB:IedB,JsdB:JedB))    ; CS%dx2q(:,:)    = 0.0
  ALLOC_(CS%dy2q(IsdB:IedB,JsdB:JedB))    ; CS%dy2q(:,:)    = 0.0
  ALLOC_(CS%dx_dyT(isd:ied,jsd:jed))      ; CS%dx_dyT(:,:)  = 0.0
  ALLOC_(CS%dy_dxT(isd:ied,jsd:jed))      ; CS%dy_dxT(:,:)  = 0.0
  ALLOC_(CS%dx_dyBu(IsdB:IedB,JsdB:JedB)) ; CS%dx_dyBu(:,:) = 0.0
  ALLOC_(CS%dy_dxBu(IsdB:IedB,JsdB:JedB)) ; CS%dy_dxBu(:,:) = 0.0

  if (CS%Laplacian) then
    ALLOC_(CS%Kh_bg_xx(isd:ied,jsd:jed))     ; CS%Kh_bg_xx(:,:) = 0.0
    ALLOC_(CS%Kh_bg_xy(IsdB:IedB,JsdB:JedB)) ; CS%Kh_bg_xy(:,:) = 0.0
    if (CS%bound_Kh .or. CS%better_bound_Kh) then
      ALLOC_(CS%Kh_Max_xx(IsdB:IedB,JsdB:JedB)) ; CS%Kh_Max_xx(:,:) = 0.0
      ALLOC_(CS%Kh_Max_xy(IsdB:IedB,JsdB:JedB)) ; CS%Kh_Max_xy(:,:) = 0.0
    endif
    if (CS%Smagorinsky_Kh) then
      ALLOC_(CS%Laplac_Const_xx(isd:ied,jsd:jed))     ; CS%Laplac_Const_xx(:,:) = 0.0
      ALLOC_(CS%Laplac_Const_xy(IsdB:IedB,JsdB:JedB)) ; CS%Laplac_Const_xy(:,:) = 0.0
    endif
  endif
  ALLOC_(CS%reduction_xx(isd:ied,jsd:jed))     ; CS%reduction_xx(:,:) = 0.0
  ALLOC_(CS%reduction_xy(IsdB:IedB,JsdB:JedB)) ; CS%reduction_xy(:,:) = 0.0

  if (CS%use_Kh_bg_2d) then
    ALLOC_(CS%Kh_bg_2d(isd:ied,jsd:jed))     ; CS%Kh_bg_2d(:,:) = 0.0
    call get_param(param_file, mod, "KH_BG_2D_FILENAME", filename, &
                 'The filename containing a 2d map of "Kh".', &
                 default='KH_background_2d.nc')
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    call read_data(trim(inputdir)//trim(filename), 'Kh', CS%Kh_bg_2d, &
                   domain=G%domain%mpp_domain, timelevel=1)
    call pass_var(CS%Kh_bg_2d, G%domain)
  endif

  if (CS%biharmonic) then
    ALLOC_(CS%Idx2dyCu(IsdB:IedB,jsd:jed)) ; CS%Idx2dyCu(:,:) = 0.0
    ALLOC_(CS%Idx2dyCv(isd:ied,JsdB:JedB)) ; CS%Idx2dyCv(:,:) = 0.0
    ALLOC_(CS%Idxdy2u(IsdB:IedB,jsd:jed))  ; CS%Idxdy2u(:,:)  = 0.0
    ALLOC_(CS%Idxdy2v(isd:ied,JsdB:JedB))  ; CS%Idxdy2v(:,:)  = 0.0

    ALLOC_(CS%Ah_bg_xx(isd:ied,jsd:jed))     ; CS%Ah_bg_xx(:,:) = 0.0
    ALLOC_(CS%Ah_bg_xy(IsdB:IedB,JsdB:JedB)) ; CS%Ah_bg_xy(:,:) = 0.0
    if (CS%bound_Ah .or. CS%better_bound_Ah) then
      ALLOC_(CS%Ah_Max_xx(isd:ied,jsd:jed))     ; CS%Ah_Max_xx(:,:) = 0.0
      ALLOC_(CS%Ah_Max_xy(IsdB:IedB,JsdB:JedB)) ; CS%Ah_Max_xy(:,:) = 0.0
    endif
    if (CS%Smagorinsky_Ah) then
      ALLOC_(CS%Biharm_Const_xx(isd:ied,jsd:jed))     ; CS%Biharm_Const_xx(:,:) = 0.0
      ALLOC_(CS%Biharm_Const_xy(IsdB:IedB,JsdB:JedB)) ; CS%Biharm_Const_xy(:,:) = 0.0
      if (CS%bound_Coriolis) then
        ALLOC_(CS%Biharm_Const2_xx(isd:ied,jsd:jed))     ; CS%Biharm_Const2_xx(:,:) = 0.0
        ALLOC_(CS%Biharm_Const2_xy(IsdB:IedB,JsdB:JedB)) ; CS%Biharm_Const2_xy(:,:) = 0.0
      endif
    endif
  endif

  do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
    CS%DX2q(I,J) = G%dxBu(I,J)*G%dxBu(I,J) ; CS%DY2q(I,J) = G%dyBu(I,J)*G%dyBu(I,J)
    CS%DX_dyBu(I,J) = G%dxBu(I,J)*G%IdyBu(I,J) ; CS%DY_dxBu(I,J) = G%dyBu(I,J)*G%IdxBu(I,J)
  enddo ; enddo
  do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
    CS%DX2h(i,j) = G%dxT(i,j)*G%dxT(i,j) ; CS%DY2h(i,j) = G%dyT(i,j)*G%dyT(i,j)
    CS%DX_dyT(i,j) = G%dxT(i,j)*G%IdyT(i,j) ; CS%DY_dxT(i,j) = G%dyT(i,j)*G%IdxT(i,j)
  enddo ; enddo

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    CS%reduction_xx(i,j) = 1.0
    if ((G%dy_Cu(I,j) > 0.0) .and. (G%dy_Cu(I,j) < G%dyCu(I,j)) .and. &
        (G%dy_Cu(I,j) < G%dyCu(I,j) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dy_Cu(I,j) / G%dyCu(I,j)
    if ((G%dy_Cu(I-1,j) > 0.0) .and. (G%dy_Cu(I-1,j) < G%dyCu(I-1,j)) .and. &
        (G%dy_Cu(I-1,j) < G%dyCu(I-1,j) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dy_Cu(I-1,j) / G%dyCu(I-1,j)
    if ((G%dx_Cv(i,J) > 0.0) .and. (G%dx_Cv(i,J) < G%dxCv(i,J)) .and. &
        (G%dx_Cv(i,J) < G%dxCv(i,J) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dx_Cv(i,J) / G%dxCv(i,J)
    if ((G%dx_Cv(i,J-1) > 0.0) .and. (G%dx_Cv(i,J-1) < G%dxCv(i,J-1)) .and. &
        (G%dx_Cv(i,J-1) < G%dxCv(i,J-1) * CS%reduction_xx(i,j))) &
      CS%reduction_xx(i,j) = G%dx_Cv(i,J-1) / G%dxCv(i,J-1)
  enddo ; enddo

  do J=js-1,Jeq ; do I=is-1,Ieq
    CS%reduction_xy(I,J) = 1.0
    if ((G%dy_Cu(I,j) > 0.0) .and. (G%dy_Cu(I,j) < G%dyCu(I,j)) .and. &
        (G%dy_Cu(I,j) < G%dyCu(I,j) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dy_Cu(I,j) / G%dyCu(I,j)
    if ((G%dy_Cu(I,j+1) > 0.0) .and. (G%dy_Cu(I,j+1) < G%dyCu(I,j+1)) .and. &
        (G%dy_Cu(I,j+1) < G%dyCu(I,j+1) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dy_Cu(I,j+1) / G%dyCu(I,j+1)
    if ((G%dx_Cv(i,J) > 0.0) .and. (G%dx_Cv(i,J) < G%dxCv(i,J)) .and. &
        (G%dx_Cv(i,J) < G%dxCv(i,J) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dx_Cv(i,J) / G%dxCv(i,J)
    if ((G%dx_Cv(i+1,J) > 0.0) .and. (G%dx_Cv(i+1,J) < G%dxCv(i+1,J)) .and. &
        (G%dx_Cv(i+1,J) < G%dxCv(i+1,J) * CS%reduction_xy(I,J))) &
      CS%reduction_xy(I,J) = G%dx_Cv(i+1,J) / G%dxCv(i+1,J)
  enddo ; enddo

  if (CS%Laplacian) then
   ! The 0.3 below was 0.4 in MOM1.10.  The change in hq requires
   ! this to be less than 1/3, rather than 1/2 as before.
    if (CS%bound_Kh .or. CS%bound_Ah) Kh_Limit = 0.3 / (dt*4.0)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      grid_sp_h2 = (2.0*CS%DX2h(i,j)*CS%DY2h(i,j)) / (CS%DX2h(i,j) + CS%DY2h(i,j))
      if (CS%Smagorinsky_Kh) CS%LAPLAC_CONST_xx(i,j) = Smag_Lap_const * grid_sp_h2

      CS%Kh_bg_xx(i,j) = MAX(Kh, Kh_vel_scale * sqrt(grid_sp_h2))

      if (CS%use_Kh_bg_2d) CS%Kh_bg_xx(i,j) = MAX(CS%Kh_bg_2d(i,j), CS%Kh_bg_xx(i,j))

      if (CS%bound_Kh .and. .not.CS%better_bound_Kh) then
        CS%Kh_Max_xx(i,j) = Kh_Limit * grid_sp_h2
        CS%Kh_bg_xx(i,j) = MIN(CS%Kh_bg_xx(i,j), CS%Kh_Max_xx(i,j))
      endif
    enddo ; enddo
    do J=js-1,Jeq ; do I=is-1,Ieq
      grid_sp_q2 = (2.0*CS%DX2q(I,J)*CS%DY2q(I,J)) / (CS%DX2q(I,J) + CS%DY2q(I,J))
      if (CS%Smagorinsky_Kh) CS%LAPLAC_CONST_xy(I,J) = Smag_Lap_const * grid_sp_q2

      CS%Kh_bg_xy(I,J) = MAX(Kh, Kh_vel_scale * sqrt(grid_sp_q2))

      if (CS%use_Kh_bg_2d) CS%Kh_bg_xy(I,J) = MAX(CS%Kh_bg_2d(i,j), CS%Kh_bg_xy(I,J))

      if (CS%bound_Kh .and. .not.CS%better_bound_Kh) then
        CS%Kh_Max_xy(I,J) = Kh_Limit * grid_sp_q2
        CS%Kh_bg_xy(I,J) = MIN(CS%Kh_bg_xy(I,J), CS%Kh_Max_xy(I,J))
      endif
    enddo ; enddo
  endif

  if (CS%biharmonic) then

    do j=js-1,Jeq+1 ; do I=Isq-1,Ieq+1
      CS%IDX2dyCu(I,j) = (G%IdxCu(I,j)*G%IdxCu(I,j)) * G%IdyCu(I,j)
      CS%IDXDY2u(I,j) = G%IdxCu(I,j) * (G%IdyCu(I,j)*G%IdyCu(I,j))
    enddo ; enddo
    do J=Jsq-1,Jeq+1 ; do i=is-1,Ieq+1
      CS%IDX2dyCv(i,J) = (G%IdxCv(i,J)*G%IdxCv(i,J)) * G%IdyCv(i,J)
      CS%IDXDY2v(i,J) = G%IdxCv(i,J) * (G%IdyCv(i,J)*G%IdyCv(i,J))
    enddo ; enddo

    CS%Ah_bg_xy(:,:) = 0.0
   ! The 0.3 below was 0.4 in MOM1.10.  The change in hq requires
   ! this to be less than 1/3, rather than 1/2 as before.
    Ah_Limit = 0.3 / (dt*64.0)
    if (CS%Smagorinsky_Ah .and. CS%bound_Coriolis) &
      BoundCorConst = 1.0 / (5.0*(bound_Cor_vel*bound_Cor_vel))
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      grid_sp_h2 = (2.0*CS%DX2h(i,j)*CS%DY2h(i,j)) / (CS%DX2h(i,j)+CS%DY2h(i,j))

      if (CS%Smagorinsky_Ah) then
        CS%BIHARM_CONST_xx(i,j) = Smag_bi_const * (grid_sp_h2 * grid_sp_h2)
        if (CS%bound_Coriolis) then
          fmax = MAX(abs(G%CoriolisBu(I-1,J-1)), abs(G%CoriolisBu(I,J-1)), &
                     abs(G%CoriolisBu(I-1,J)),   abs(G%CoriolisBu(I,J)))
          CS%Biharm_Const2_xx(i,j) = (grid_sp_h2 * grid_sp_h2 * grid_sp_h2) * &
                                  (fmax * BoundCorConst)
        endif
      endif
      CS%Ah_bg_xx(i,j) = MAX(Ah, Ah_vel_scale * grid_sp_h2 * sqrt(grid_sp_h2))
      if (CS%bound_Ah .and. .not.CS%better_bound_Ah) then
        CS%Ah_Max_xx(i,j) = Ah_Limit * (grid_sp_h2 * grid_sp_h2)
        CS%Ah_bg_xx(i,j) = MIN(CS%Ah_bg_xx(i,j), CS%Ah_Max_xx(i,j))
      endif
    enddo ; enddo
    do J=js-1,Jeq ; do I=is-1,Ieq
      grid_sp_q2 = (2.0*CS%DX2q(I,J)*CS%DY2q(I,J)) / (CS%DX2q(I,J)+CS%DY2q(I,J))
      if (CS%Smagorinsky_Ah) then
        CS%BIHARM_CONST_xy(I,J) = Smag_bi_const * (grid_sp_q2 * grid_sp_q2)
        if (CS%bound_Coriolis) then
          CS%Biharm_Const2_xy(I,J) = (grid_sp_q2 * grid_sp_q2 * grid_sp_q2) * &
                                      (abs(G%CoriolisBu(I,J)) * BoundCorConst)
        endif
      endif
      CS%Ah_bg_xy(I,J) = MAX(Ah, Ah_vel_scale * grid_sp_q2 * sqrt(grid_sp_q2))
      if (CS%bound_Ah .and. .not.CS%better_bound_Ah) then
        CS%Ah_Max_xy(I,J) = Ah_Limit * (grid_sp_q2 * grid_sp_q2)
        CS%Ah_bg_xy(I,J) = MIN(CS%Ah_bg_xy(I,J), CS%Ah_Max_xy(I,J))
      endif
    enddo ; enddo
  endif

  ! The Laplacian bounds should avoid overshoots when CS%bound_coef < 1.
  if (CS%Laplacian .and. CS%better_bound_Kh) then
    Idt = 1.0 / dt
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      denom = max( &
         (CS%DY2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j)) * &
          max(G%IdyCu(I,j)*G%IareaCu(I,j), G%IdyCu(I-1,j)*G%IareaCu(I-1,j)) ), &
         (CS%DX2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1)) * &
          max(G%IdxCv(i,J)*G%IareaCv(i,J), G%IdxCv(i,J-1)*G%IareaCv(i,J-1)) ) )
      CS%Kh_Max_xx(i,j) = 0.0
      if (denom > 0.0) &
        CS%Kh_Max_xx(i,j) = CS%bound_coef * 0.25 * Idt / denom
    enddo ; enddo
    do J=js-1,Jeq ; do I=is-1,Ieq
      denom = max( &
         (CS%DX2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j)) * &
          max(G%IdxCu(I,j)*G%IareaCu(I,j), G%IdxCu(I,j+1)*G%IareaCu(I,j+1)) ), &
         (CS%DY2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J)) * &
          max(G%IdyCv(i,J)*G%IareaCv(i,J), G%IdyCv(i+1,J)*G%IareaCv(i+1,J)) ) )
      CS%Kh_Max_xy(I,J) = 0.0
      if (denom > 0.0) &
        CS%Kh_Max_xy(I,J) = CS%bound_coef * 0.25 * Idt / denom
    enddo ; enddo
  endif

  ! The biharmonic bounds should avoid overshoots when CS%bound_coef < 0.5, but
  ! empirically work for CS%bound_coef <~ 1.0
  if (CS%biharmonic .and. CS%better_bound_Ah) then
    Idt = 1.0 / dt
    do j=js-1,Jeq+1 ; do I=Isq-1,Ieq+1
      u0u(I,j) = CS%IDXDY2u(I,j)*(CS%DY2h(i+1,j)*CS%DY_dxT(i+1,j)*(G%IdyCu(I+1,j) + G%IdyCu(I,j))   + &
                                  CS%DY2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j)) ) + &
                 CS%IDX2dyCu(I,j)*(CS%DX2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j)) + &
                                  CS%DX2q(I,J-1)*CS%DX_dyBu(I,J-1)*(G%IdxCu(I,j) + G%IdxCu(I,j-1)) )

      u0v(I,j) = CS%IDXDY2u(I,j)*(CS%DY2h(i+1,j)*CS%DX_dyT(i+1,j)*(G%IdxCv(i+1,J) + G%IdxCv(i+1,J-1)) + &
                                  CS%DY2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1)) )   + &
                 CS%IDX2dyCu(I,j)*(CS%DX2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J))   + &
                                  CS%DX2q(I,J-1)*CS%DY_dxBu(I,J-1)*(G%IdyCv(i+1,J-1) + G%IdyCv(i,J-1)) )
    enddo ; enddo
    do J=Jsq-1,Jeq+1 ; do i=is-1,Ieq+1
      v0u(i,J) = CS%IDXDY2v(i,J)*(CS%DY2q(I,J) * CS%DX_dyBu(I,J) * (G%IdxCu(I,j+1) + G%IdxCu(I,j))       + &
                                  CS%DY2q(I-1,J)*CS%DX_dyBu(I-1,J)*(G%IdxCu(I-1,j+1) + G%IdxCu(I-1,j)) ) + &
                 CS%IDX2dyCv(i,J)*(CS%DX2h(i,j+1)*CS%DY_dxT(i,j+1)*(G%IdyCu(I,j+1) + G%IdyCu(I-1,j+1))   + &
                                  CS%DX2h(i,j) * CS%DY_dxT(i,j) * (G%IdyCu(I,j) + G%IdyCu(I-1,j)) )

      v0v(i,J) = CS%IDXDY2v(i,J)*(CS%DY2q(I,J) * CS%DY_dxBu(I,J) * (G%IdyCv(i+1,J) + G%IdyCv(i,J))   + &
                                  CS%DY2q(I-1,J)*CS%DY_dxBu(I-1,J)*(G%IdyCv(i,J) + G%IdyCv(i-1,J)) ) + &
                 CS%IDX2dyCv(i,J)*(CS%DX2h(i,j+1)*CS%DX_dyT(i,j+1)*(G%IdxCv(i,J+1) + G%IdxCv(i,J))   + &
                                  CS%DX2h(i,j) * CS%DX_dyT(i,j) * (G%IdxCv(i,J) + G%IdxCv(i,J-1)) )
    enddo ; enddo

    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      denom = max( &
         (CS%DY2h(i,j) * &
          (CS%DY_dxT(i,j)*(G%IdyCu(I,j)*u0u(I,j) + G%IdyCu(I-1,j)*u0u(I-1,j))  + &
           CS%DX_dyT(i,j)*(G%IdxCv(i,J)*v0u(i,J) + G%IdxCv(i,J-1)*v0u(i,J-1))) * &
          max(G%IdyCu(I,j)*G%IareaCu(I,j), G%IdyCu(I-1,j)*G%IareaCu(I-1,j)) ),   &
         (CS%DX2h(i,j) * &
          (CS%DY_dxT(i,j)*(G%IdyCu(I,j)*u0v(I,j) + G%IdyCu(I-1,j)*u0v(I-1,j))  + &
           CS%DX_dyT(i,j)*(G%IdxCv(i,J)*v0v(i,J) + G%IdxCv(i,J-1)*v0v(i,J-1))) * &
          max(G%IdxCv(i,J)*G%IareaCv(i,J), G%IdxCv(i,J-1)*G%IareaCv(i,J-1)) ) )
      CS%Ah_Max_xx(I,J) = 0.0
      if (denom > 0.0) &
        CS%Ah_Max_xx(I,J) = CS%bound_coef * 0.5 * Idt / denom
    enddo ; enddo

    do J=js-1,Jeq ; do I=is-1,Ieq
      denom = max( &
         (CS%DX2q(I,J) * &
          (CS%DX_dyBu(I,J)*(u0u(I,j+1)*G%IdxCu(I,j+1) + u0u(I,j)*G%IdxCu(I,j))  + &
           CS%DY_dxBu(I,J)*(v0u(i+1,J)*G%IdyCv(i+1,J) + v0u(i,J)*G%IdyCv(i,J))) * &
          max(G%IdxCu(I,j)*G%IareaCu(I,j), G%IdxCu(I,j+1)*G%IareaCu(I,j+1)) ),    &
         (CS%DY2q(I,J) * &
          (CS%DX_dyBu(I,J)*(u0v(I,j+1)*G%IdxCu(I,j+1) + u0v(I,j)*G%IdxCu(I,j))  + &
           CS%DY_dxBu(I,J)*(v0v(i+1,J)*G%IdyCv(i+1,J) + v0v(i,J)*G%IdyCv(i,J))) * &
          max(G%IdyCv(i,J)*G%IareaCv(i,J), G%IdyCv(i+1,J)*G%IareaCv(i+1,J)) ) )
      CS%Ah_Max_xy(I,J) = 0.0
      if (denom > 0.0) &
        CS%Ah_Max_xy(I,J) = CS%bound_coef * 0.5 * Idt / denom
    enddo ; enddo
  endif

  ! Register fields for output from this module.

  CS%id_diffu = register_diag_field('ocean_model', 'diffu', diag%axesCuL, Time, &
      'Zonal Acceleration from Horizontal Viscosity', 'meter second-2')

  CS%id_diffv = register_diag_field('ocean_model', 'diffv', diag%axesCvL, Time, &
      'Meridional Acceleration from Horizontal Viscosity', 'meter second-2')

  if (CS%biharmonic) then
    CS%id_Ah_h = register_diag_field('ocean_model', 'Ahh', diag%axesTL, Time,    &
        'Biharmonic Horizontal Viscosity at h Points', 'meter4 second-1',        &
        cmor_field_name='difmxybo', cmor_units='m4 s-1',                        &
        cmor_long_name='Ocean lateral biharmonic viscosity',                     &
        cmor_standard_name='ocean_momentum_xy_biharmonic_diffusivity')

    CS%id_Ah_q = register_diag_field('ocean_model', 'Ahq', diag%axesBL, Time, &
        'Biharmonic Horizontal Viscosity at q Points', 'meter4 second-1')
  endif

  if (CS%Laplacian) then
    CS%id_Kh_h = register_diag_field('ocean_model', 'Khh', diag%axesTL, Time,   &
        'Laplacian Horizontal Viscosity at h Points', 'meter2 second-1',        &
        cmor_field_name='difmxylo', cmor_units='m2 s-1',                        &
        cmor_long_name='Ocean lateral Laplacian viscosity',                     &
        cmor_standard_name='ocean_momentum_xy_laplacian_diffusivity')

    CS%id_Kh_q = register_diag_field('ocean_model', 'Khq', diag%axesBL, Time, &
        'Laplacian Horizontal Viscosity at q Points', 'meter2 second-1')
  endif

  CS%id_FrictWork = register_diag_field('ocean_model','FrictWork',diag%axesTL,Time,&
      'Integral work done by lateral friction terms', 'Watt meter-2')

  CS%id_FrictWorkIntz = register_diag_field('ocean_model','FrictWorkIntz',diag%axesT1,Time,      &
      'Depth integrated work done by lateral friction', 'Watt meter-2',                          &
      cmor_field_name='dispkexyfo', cmor_units='W m-2',                                          &
      cmor_long_name='Depth integrated ocean kinetic energy dissipation due to lateral friction',&
      cmor_standard_name='ocean_kinetic_energy_dissipation_per_unit_area_due_to_xy_friction')

  if (CS%Laplacian .or. get_all) then
  endif

end subroutine hor_visc_init

subroutine hor_visc_end(CS)
! This subroutine deallocates any variables allocated in hor_visc_init.
! Argument:  CS - The control structure returned by a previous call to
!                 hor_visc_init.
  type(hor_visc_CS), pointer :: CS

  DEALLOC_(CS%dx2h) ; DEALLOC_(CS%dx2q) ; DEALLOC_(CS%dy2h) ; DEALLOC_(CS%dy2q)
  DEALLOC_(CS%dx_dyT) ; DEALLOC_(CS%dy_dxT) ; DEALLOC_(CS%dx_dyBu) ; DEALLOC_(CS%dy_dxBu)

  DEALLOC_(CS%reduction_xx) ; DEALLOC_(CS%reduction_xy)

  if (CS%Laplacian) then
    DEALLOC_(CS%Kh_bg_xx) ; DEALLOC_(CS%Kh_bg_xy)
    if (CS%bound_Kh) then
      DEALLOC_(CS%Kh_Max_xx) ; DEALLOC_(CS%Kh_Max_xy)
    endif
    if (CS%Smagorinsky_Kh) then
      DEALLOC_(CS%Laplac_Const_xx) ; DEALLOC_(CS%Laplac_Const_xy)
    endif
  endif

  if (CS%biharmonic) then
    DEALLOC_(CS%Idx2dyCu) ; DEALLOC_(CS%Idx2dyCv)
    DEALLOC_(CS%Idxdy2u) ; DEALLOC_(CS%Idxdy2v)
    DEALLOC_(CS%Ah_bg_xx) ; DEALLOC_(CS%Ah_bg_xy)
    if (CS%bound_Ah) then
      DEALLOC_(CS%Ah_Max_xx) ; DEALLOC_(CS%Ah_Max_xy)
    endif
    if (CS%Smagorinsky_Ah) then
      DEALLOC_(CS%Biharm_Const_xx) ; DEALLOC_(CS%Biharm_Const_xy)
      if (CS%bound_Coriolis) then
        DEALLOC_(CS%Biharm_Const2_xx) ; DEALLOC_(CS%Biharm_Const2_xy)
      endif
    endif
  endif
  deallocate(CS)

end subroutine hor_visc_end

end module MOM_hor_visc
