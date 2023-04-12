! > Calculates Zanna and Bolton 2020 parameterization
module MOM_Zanna_Bolton

use MOM_grid,          only : ocean_grid_type
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_diag_mediator, only : post_data, register_diag_field
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,       only : To_North, To_East
use MOM_domains,       only : pass_var, CORNER
use MOM_coms,          only : reproducing_sum, max_across_PEs, min_across_PEs

implicit none ; private

#include <MOM_memory.h>

public Zanna_Bolton_2020, ZB_2020_init

!> Control structure that contains MEKE parameters and diagnostics handles
type, public :: ZB2020_CS
  ! Parameters
  logical   :: use_ZB2020     !< If true, parameterization works
  real      :: amplitude      !< k_bc = - amplitude * cell_area
  real      :: amp_bottom     !< amplitude in the bottom layer; -1 = use same
  integer   :: ZB_type        !< 0 = Full model, 1 = trace-free part, 2 = only trace part
  integer   :: ZB_cons        !< 0: nonconservative; 1: conservative without interface;
  integer   :: LPF_iter       !< Low-pass filter for Velocity gradient; number of iterations
  integer   :: LPF_order      !< Low-pass filter for Velocity gradient; 1: Laplacian, 2: Bilaplacian
  integer   :: HPF_iter       !< High-pass filter for Velocity gradient; number of iterations
  integer   :: HPF_order      !< High-pass filter for Velocity gradient; 1: Laplacian, 2: Bilaplacian
  integer   :: Stress_iter    !< Low-pass filter for Stress (Momentum Flux); number of iterations
  integer   :: Stress_order   !< Low-pass filter for Stress (Momentum Flux); 1: Laplacian, 2: Bilaplacian
  integer   :: ssd_iter       !< Small-scale dissipation in RHS of momentum eq; -1: off, 0:Laplacian, 4:Laplacian^5
  real      :: ssd_bound_coef !< the viscosity bounds to the theoretical maximum for stability

  real      :: DT            !< The (baroclinic) dynamics time step.
  
  type(diag_ctrl), pointer :: diag => NULL() !< A type that regulates diagnostics output
  !>@{ Diagnostic handles
  integer :: id_ZB2020u = -1, id_ZB2020v = -1, id_KE_ZB2020 = -1, id_kbc = -1
  integer :: id_maskT = -1
  integer :: id_maskq = -1
  integer :: id_S_11f = -1
  integer :: id_S_22f = -1
  integer :: id_S_12f = -1
  !>@}

end type ZB2020_CS

contains

subroutine ZB_2020_init(Time, GV, US, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time       !< The current model time.
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file parser structure.
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(ZB2020_CS),         intent(inout) :: CS         !< ZB2020 control structure.

  ! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_Zanna_Bolton" ! This module's name.

  call log_version(param_file, mdl, version, "")
  
  call get_param(param_file, mdl, "USE_ZB2020", CS%use_ZB2020, &
                 "If true, turns on Zanna-Bolton 2020 parameterization", &
                 default=.true.)

  call get_param(param_file, mdl, "amplitude", CS%amplitude, &
                 "k_bc=-amplitude*cell_area, amplitude=1/24..1", &
                 units="nondim", default=1./24.)
  
  call get_param(param_file, mdl, "amp_bottom", CS%amp_bottom, &
                 "-1=use same amplitude, or specify", &
                 units="nondim", default=-1.)

  if (CS%amp_bottom < -0.5) CS%amp_bottom = CS%amplitude
  
  call get_param(param_file, mdl, "ZB_type", CS%ZB_type, &
                 "Type of parameterization: 0 = full, 1 = trace-free, 2 = trace-only", &
                 default=0)

  call get_param(param_file, mdl, "ZB_cons", CS%ZB_cons, &
                 "0: nonconservative; 1: conservative without interface", &
                 default=1)

  call get_param(param_file, mdl, "LPF_iter", CS%LPF_iter, &
                 "Low-pass filter for Velocity gradient; number of iterations", &
                 default=2)
  
  call get_param(param_file, mdl, "LPF_order", CS%LPF_order, &
                 "Low-pass filter for Velocity gradient; 1: Laplacian, 2: Bilaplacian", &
                 default=2)

  call get_param(param_file, mdl, "HPF_iter", CS%HPF_iter, &
                 "High-pass filter for Velocity gradient; number of iterations", &
                 default=2)
  
  call get_param(param_file, mdl, "HPF_order", CS%HPF_order, &
                 "High-pass filter for Velocity gradient; 1: Laplacian, 2: Bilaplacian", &
                 default=2)

  call get_param(param_file, mdl, "Stress_iter", CS%Stress_iter, &
                 "Low-pass filter for Stress (Momentum Flux); number of iterations", &
                 default=2)
  
  call get_param(param_file, mdl, "Stress_order", CS%Stress_order, &
                 "Low-pass filter for Stress (Momentum Flux); 1: Laplacian, 2: Bilaplacian", &
                 default=2)

  call get_param(param_file, mdl, "ssd_iter", CS%ssd_iter, &
                 "Small-scale dissipation in RHS of momentum eq; -1: off, 0:Laplacian, 4:Laplacian^5", &
                 default=-1)

  call get_param(param_file, mdl, "ssd_bound_coef", CS%ssd_bound_coef, &
                 "The viscosity bounds to the theoretical maximum for stability", units="nondim", &
                 default=0.8)
  
  call get_param(param_file, mdl, "DT", CS%dt, &
                 "The (baroclinic) dynamics time step.", units="s", scale=US%s_to_T, &
                 fail_if_missing=.true.)

  ! Register fields for output from this module.
  CS%diag => diag

  CS%id_ZB2020u = register_diag_field('ocean_model', 'ZB2020u', diag%axesCuL, Time, &
      'Zonal Acceleration from Zanna-Bolton 2020', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_ZB2020v = register_diag_field('ocean_model', 'ZB2020v', diag%axesCvL, Time, &
      'Meridional Acceleration from Zanna-Bolton 2020', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_KE_ZB2020 = register_diag_field('ocean_model', 'KE_ZB2020', diag%axesTL, Time, &
      'Kinetic Energy Source from Horizontal Viscosity', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_kbc = register_diag_field('ocean_model', 'kbc_ZB2020', diag%axesTL, Time, &
      'Kinetic Energy Source from Horizontal Viscosity', &
      'm2', conversion=US%L_to_m**2)

  ! masks and action of filter on test fields
  CS%id_maskT  = register_diag_field('ocean_model', 'maskT', diag%axesTL, Time, &
      'mask of wet T points', '1', conversion=1.)
  
  CS%id_maskq  = register_diag_field('ocean_model', 'maskq', diag%axesBL, Time, &
      'mask of wet q points', '1', conversion=1.)
  
  ! action of filter on momentum flux
  CS%id_S_11f = register_diag_field('ocean_model', 'S_11f', diag%axesTL, Time, &
      '11 momentum flux filtered', 'm2s-2', conversion=US%L_T_to_m_s**2)

  CS%id_S_22f = register_diag_field('ocean_model', 'S_22f', diag%axesTL, Time, &
      '22 momentum flux filtered', 'm2s-2', conversion=US%L_T_to_m_s**2)

  CS%id_S_12f = register_diag_field('ocean_model', 'S_12f', diag%axesBL, Time, &
      '12 momentum flux filtered', 'm2s-2', conversion=US%L_T_to_m_s**2)
  
end subroutine ZB_2020_init

!> Baroclinic parameterization is as follows:
!! eq. 6 in https://laurezanna.github.io/files/Zanna-Bolton-2020.pdf
!! (du/dt, dv/dt) = k_BC * 
!!                  (div(S0) + 1/2 * grad(vort_xy^2 + sh_xy^2 + sh_xx^2))
!! vort_xy = dv/dx - du/dy - relative vorticity
!! sh_xy   = dv/dx + du/dy - shearing deformation (or horizontal shear strain)
!! sh_xx   = du/dx - dv/dy - stretching deformation (or horizontal tension)
!! S0 - 2x2 tensor:
!! S0 = vort_xy * (-sh_xy, sh_xx; sh_xx, sh_xy)
!! Relating k_BC to velocity gradient model,
!! k_BC = - amplitude * cell_area
!! where amplitude = 1/24..1 (approx)
!! 
!! S - is a tensor of full tendency
!! S = (-vort_xy * sh_xy + 1/2 * (vort_xy^2 + sh_xy^2 + sh_xx^2), vort_xy * sh_xx;
!!       vort_xy * sh_xx, vort_xy * sh_xy + 1/2 * (vort_xy^2 + sh_xy^2 + sh_xx^2))
!! So the full parameterization:
!! (du/dt, dv/dt) = k_BC * div(S)
!! In generalized curvilinear orthogonal coordinates (see Griffies 2020,
!! and MOM documentation 
!! https://mom6.readthedocs.io/en/dev-gfdl/api/generated/modules/mom_hor_visc.html#f/mom_hor_visc):
!! du/dx -> dy/dx * delta_i (u / dy)
!! dv/dy -> dx/dy * delta_j (v / dx)
!! dv/dx -> dy/dx * delta_i (v / dy)
!! du/dy -> dx/dy * delta_j (u / dx)
!!
!! vort_xy and sh_xy are in the corner of the cell
!! sh_xx in the center of the cell
!!
!! In order to compute divergence of S, its components must be:
!! S_11, S_22 in center of the cells
!! S_12 (=S_21) in the corner
!!
!! The following interpolations are required:
!! sh_xx center -> corner
!! vort_xy, sh_xy corner -> center
subroutine Zanna_Bolton_2020(u, v, h, fx, fy, G, GV, CS)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  type(ZB2020_CS),               intent(in)  :: CS     !< ZB2020 control structure.

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2].
  
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(out) :: fx     !< Zonal acceleration due to convergence of
                                                       !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(out) :: fy     !< Meridional acceleration due to convergence
                                                       !! of along-coordinate stress tensor [L T-2 ~> m s-2]

  real, dimension(SZI_(G),SZJ_(G)) :: &
    dx_dyT, &          !< Pre-calculated dx/dy at h points [nondim]
    dy_dxT, &          !< Pre-calculated dy/dx at h points [nondim]
    dx2h, &            !< Pre-calculated dx^2 at h points [L2 ~> m2]
    dy2h, &            !< Pre-calculated dy^2 at h points [L2 ~> m2]
    dudx, dvdy, &      ! components in the horizontal tension [T-1 ~> s-1]
    sh_xx, &           ! horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    vort_xy_center, &  ! vort_xy in the center
    sh_xy_center, &    ! sh_xy in the center
    S_11, S_22, &      ! flux tensor in the cell center, multiplied with interface height [m^2/s^2 * h]
    ssd_11, &          ! diagonal part of ssd in cell center
    mask_T             ! mask of wet center points

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dx_dyBu, &         !< Pre-calculated dx/dy at q points [nondim]
    dy_dxBu, &         !< Pre-calculated dy/dx at q points [nondim]
    dx2q, &            !< Pre-calculated dx^2 at q points [L2 ~> m2]
    dy2q, &            !< Pre-calculated dy^2 at q points [L2 ~> m2]
    dvdx, dudy, &      ! components in the shearing strain [T-1 ~> s-1]
    vort_xy, &         ! Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
    sh_xy, &           ! horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
    sh_xx_corner, &    ! sh_xx in the corner
    S_12, &            ! flux tensor in the corner, multiplied with interface height [m^2/s^2 * h]
    ssd_12, &          ! off-diagonal part of ssd in corner
    hq, &              ! harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2]
    mask_q             ! mask of wet corner points

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    kbc_3d,    &          ! k_bc parameter as 3d field [L2 ~> m2]
    mask_T_3d, &          
    S_11_3df,  &                    
    S_22_3df

    real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)) :: &
    mask_q_3d, &          
    S_12_3df      

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    h_u                ! Thickness interpolated to u points [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    h_v                ! Thickness interpolated to v points [H ~> m or kg m-2].

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n
  
  real :: h_neglect  ! thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
  real :: h_neglect3 ! h_neglect^3 [H3 ~> m3 or kg3 m-6]
  real :: h2uq, h2vq ! temporary variables [H2 ~> m2 or kg2 m-4].

  real :: sum_sq     ! squared sum, i.e. 1/2*(vort_xy^2 + sh_xy^2 + sh_xx^2)
  real :: vort_sh    ! multiplication of vort_xt and sh_xy
  real :: amplitude  ! amplitude of ZB parameterization

  real :: k_bc ! free constant in parameterization, k_bc < 0, [k_bc] = m^2  

  ! Line 407 of MOM_hor_visc.F90
  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  h_neglect  = GV%H_subroundoff ! Line 410 on MOM_hor_visc.F90
  h_neglect3 = h_neglect**3

  fx(:,:,:) = 0.
  fy(:,:,:) = 0.
  kbc_3d(:,:,:) = 0.

  ! Calculate metric terms (line 2119 of MOM_hor_visc.F90)
  do J=js-2,Jeq+1 ; do I=is-2,Ieq+1
    dx2q(I,J) = G%dxBu(I,J)*G%dxBu(I,J) ; dy2q(I,J) = G%dyBu(I,J)*G%dyBu(I,J)
    DX_dyBu(I,J) = G%dxBu(I,J)*G%IdyBu(I,J) ; DY_dxBu(I,J) = G%dyBu(I,J)*G%IdxBu(I,J)
  enddo ; enddo

  ! Calculate metric terms (line 2122 of MOM_hor_visc.F90)
  do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
    dx2h(i,j) = G%dxT(i,j)*G%dxT(i,j) ; dy2h(i,j) = G%dyT(i,j)*G%dyT(i,j)
    DX_dyT(i,j) = G%dxT(i,j)*G%IdyT(i,j) ; DY_dxT(i,j) = G%dyT(i,j)*G%IdxT(i,j)
  enddo ; enddo

  do k=1,nz

    if (k==1) amplitude = CS%amplitude
    if (k==2) amplitude = CS%amp_bottom

    sh_xx(:,:) = 0.
    sh_xy(:,:) = 0.
    vort_xy(:,:) = 0.
    S_12(:,:) = 0.
    S_11(:,:) = 0.
    S_22(:,:) = 0.

    ! Calculate horizontal tension (line 590 of MOM_hor_visc.F90)
    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      dudx(i,j) = DY_dxT(i,j)*(G%IdyCu(I,j) * u(I,j,k) - &
                                  G%IdyCu(I-1,j) * u(I-1,j,k))
      dvdy(i,j) = DX_dyT(i,j)*(G%IdxCv(i,J) * v(i,J,k) - &
                                  G%IdxCv(i,J-1) * v(i,J-1,k))
      sh_xx(i,j) = dudx(i,j) - dvdy(i,j) ! center of the cell
    enddo ; enddo

    ! Components for the shearing strain (line 599 of MOM_hor_visc.F90)
    do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
      dvdx(I,J) = DY_dxBu(I,J)*(v(i+1,J,k)*G%IdyCv(i+1,J) - v(i,J,k)*G%IdyCv(i,J))
      dudy(I,J) = DX_dyBu(I,J)*(u(I,j+1,k)*G%IdxCu(I,j+1) - u(I,j,k)*G%IdxCu(I,j))
    enddo ; enddo

    ! Shearing strain with free-slip B.C. (line 751 of MOM_hor_visc.F90)
    ! We use free-slip as cannot guarantee that non-diagonal stress
    ! will accelerate or decelerate currents
    ! Note that as there is no stencil operator, set of indices
    ! is identical to the previous loop, compared to MOM_hor_visc.F90
    do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
      sh_xy(I,J) = G%mask2dBu(I,J) * ( dvdx(I,J) + dudy(I,J) ) ! corner of the cell
    enddo ; enddo

    ! Relative vorticity with free-slip B.C. (line 789 of MOM_hor_visc.F90)
    do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
      vort_xy(I,J) = G%mask2dBu(I,J) * ( dvdx(I,J) - dudy(I,J) ) ! corner of the cell
    enddo ; enddo
    
    call compute_masks(G, GV, h, mask_T, mask_q, k)
    mask_T_3d(:,:,k) = mask_T(:,:)
    mask_q_3d(:,:,k) = mask_q(:,:)

    ! Numerical scheme for ZB2020 requires 
    ! interpolation center <-> corner
    ! This interpolation requires B.C.,
    ! and that is why B.C. for Velocity Gradients should be
    ! well defined
    ! The same B.C. will be used by all filtering operators,
    ! So, it must be applied
    sh_xx(:,:) = sh_xx(:,:) * mask_T(:,:)
    sh_xy(:,:) = sh_xy(:,:) * mask_q(:,:)
    vort_xy(:,:) = vort_xy(:,:) * mask_q(:,:)

    if (CS%ssd_iter > -1) then
      ssd_11(:,:) = sh_xx(:,:) * &
                  CS%ssd_bound_coef * 0.25 / CS%DT * &
                  dx2h(:,:) * dy2h(:,:) / (dx2h(:,:) + dy2h(:,:))
      ssd_12(:,:) = sh_xy(:,:) * &
                  CS%ssd_bound_coef * 0.25 / CS%DT * &
                  dx2q(:,:) * dy2q(:,:) / (dx2q(:,:) + dy2q(:,:))
     
      if (CS%ssd_iter > 0) then
        call filter(G, mask_T, mask_q, -1, CS%ssd_iter, T=ssd_11)
        call filter(G, mask_T, mask_q, -1, CS%ssd_iter, q=ssd_12)
      endif
    endif

    call filter(G, mask_T, mask_q, -CS%HPF_iter, CS%HPF_order, T=sh_xx)
    call filter(G, mask_T, mask_q, +CS%LPF_iter, CS%LPF_order, T=sh_xx)

    call filter(G, mask_T, mask_q, -CS%HPF_iter, CS%HPF_order, q=sh_xy)
    call filter(G, mask_T, mask_q, +CS%LPF_iter, CS%LPF_order, q=sh_xy)

    call filter(G, mask_T, mask_q, -CS%HPF_iter, CS%HPF_order, q=vort_xy)
    call filter(G, mask_T, mask_q, +CS%LPF_iter, CS%LPF_order, q=vort_xy)
    
    ! Corner to center interpolation (line 901 of MOM_hor_visc.F90)
    ! lower index as in loop for sh_xy, but minus 1
    ! upper index is identical 
    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      sh_xy_center(i,j) = 0.25 * ( (sh_xy(I-1,J-1) + sh_xy(I,J)) &
                                 + (sh_xy(I-1,J) + sh_xy(I,J-1)) )
      vort_xy_center(i,j) = 0.25 * ( (vort_xy(I-1,J-1) + vort_xy(I,J)) &
                                   + (vort_xy(I-1,J) + vort_xy(I,J-1)) )
    enddo ; enddo

    ! Center to corner interpolation
    ! lower index as in loop for sh_xx
    ! upper index as in the same loop, but minus 1
    do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
      sh_xx_corner(I,J) = 0.25 * ( (sh_xx(i+1,j+1) + sh_xx(i,j)) &
                                 + (sh_xx(i+1,j) + sh_xx(i,j+1)))
    enddo ; enddo

    ! WITH land mask (line 622 of MOM_hor_visc.F90)
    ! Use of mask eliminates dependence on the 
    ! values on land
    do j=js-2,je+2 ; do I=Isq-1,Ieq+1
      h_u(I,j) = 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i+1,j)*h(i+1,j,k))
    enddo ; enddo
    do J=Jsq-1,Jeq+1 ; do i=is-2,ie+2
      h_v(i,J) = 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i,j+1)*h(i,j+1,k))
    enddo ; enddo

    ! Line 1187 of MOM_hor_visc.F90
    do J=js-1,Jeq ; do I=is-1,Ieq
      h2uq = 4.0 * (h_u(I,j) * h_u(I,j+1))
      h2vq = 4.0 * (h_v(i,J) * h_v(i+1,J))
      hq(I,J) = (2.0 * (h2uq * h2vq)) &
          / (h_neglect3 + (h2uq + h2vq) * ((h_u(I,j) + h_u(I,j+1)) + (h_v(i,J) + h_v(i+1,J))))
    enddo ; enddo

    ! Form S_11 and S_22 tensors
    ! Indices - intersection of loops for
    ! sh_xy_center and sh_xx
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      if (CS%ZB_type == 1) then 
        sum_sq = 0.
      else
        sum_sq = 0.5 * &
        (vort_xy_center(i,j)**2 + sh_xy_center(i,j)**2 + sh_xx(i,j)**2)
      endif
      
      if (CS%ZB_type == 2) then
        vort_sh = 0.
      else
        if (CS%ZB_cons == 1) then
          vort_sh = 0.25 * (                                          &
            (G%areaBu(I-1,J-1) * vort_xy(I-1,J-1) * sh_xy(I-1,J-1)  + &
            G%areaBu(I  ,J  ) * vort_xy(I  ,J  ) * sh_xy(I  ,J  )) + &
            (G%areaBu(I-1,J  ) * vort_xy(I-1,J  ) * sh_xy(I-1,J  )  + &
            G%areaBu(I  ,J-1) * vort_xy(I  ,J-1) * sh_xy(I  ,J-1))   &
            ) * G%IareaT(i,j)
        else if (CS%ZB_cons == 0) then
          vort_sh = vort_xy_center(i,j) * sh_xy_center(i,j)
        endif
      endif
      k_bc = - amplitude * G%areaT(i,j)
      S_11(i,j) = k_bc * (- vort_sh + sum_sq)
      S_22(i,j) = k_bc * (+ vort_sh + sum_sq)

      kbc_3d(i,j,k) = k_bc
    enddo ; enddo

    ! Form S_12 tensor
    ! indices correspond to sh_xx_corner loop
    do J=Jsq-1,Jeq ; do I=Isq-1,Ieq
      if (CS%ZB_type == 2) then
        vort_sh = 0.
      else 
        vort_sh = vort_xy(I,J) * sh_xx_corner(I,J)
      endif
      k_bc = - amplitude * G%areaBu(i,j)
      S_12(I,J) = k_bc * vort_sh
    enddo ; enddo

    call filter(G, mask_T, mask_q, CS%Stress_iter, CS%Stress_order, T=S_11)
    call filter(G, mask_T, mask_q, CS%Stress_iter, CS%Stress_order, T=S_22)
    call filter(G, mask_T, mask_q, CS%Stress_iter, CS%Stress_order, q=S_12)
   
    S_11_3df(:,:,k) = S_11(:,:)
    S_22_3df(:,:,k) = S_22(:,:)
    S_12_3df(:,:,k) = S_12(:,:)

    if (CS%ssd_iter>-1) then
      S_11(:,:) = S_11(:,:) + ssd_11(:,:)
      S_12(:,:) = S_12(:,:) + ssd_12(:,:)
      S_22(:,:) = S_22(:,:) - ssd_11(:,:)
    endif

    ! Weight with interface height (Line 1478 of MOM_hor_visc.F90)
    ! Note that reduction is removed
    do J=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      S_11(i,j) = S_11(i,j) * h(i,j,k)
      S_22(i,j) = S_22(i,j) * h(i,j,k)
    enddo ; enddo

    ! Free slip (Line 1487 of MOM_hor_visc.F90)
    do J=js-1,Jeq ; do I=is-1,Ieq
      S_12(I,J) = S_12(I,J) * (hq(I,J) * G%mask2dBu(I,J))
    enddo ; enddo

    ! Evaluate 1/h x.Div(h S) (Line 1495 of MOM_hor_visc.F90)
    ! Minus occurs because in original file (du/dt) = - div(S),
    ! but here is the discretization of div(S)
    do j=js,je ; do I=Isq,Ieq
      fx(I,j,k) = - ((G%IdyCu(I,j)*(dy2h(i,j)  *S_11(i,j) - &
                                    dy2h(i+1,j)*S_11(i+1,j)) + &
                      G%IdxCu(I,j)*(dx2q(I,J-1)*S_12(I,J-1) - &
                                    dx2q(I,J)  *S_12(I,J))) * &
                      G%IareaCu(I,j)) / (h_u(I,j) + h_neglect)
    enddo ; enddo

    ! Evaluate 1/h y.Div(h S) (Line 1517 of MOM_hor_visc.F90)
    do J=Jsq,Jeq ; do i=is,ie
      fy(i,J,k) = - ((G%IdyCv(i,J)*(dy2q(I-1,J)*S_12(I-1,J) - &
                                    dy2q(I,J)  *S_12(I,J)) + & ! NOTE this plus
                      G%IdxCv(i,J)*(dx2h(i,j)  *S_22(i,j) - &
                                    dx2h(i,j+1)*S_22(i,j+1))) * &
                      G%IareaCv(i,J)) / (h_v(i,J) + h_neglect)
    enddo ; enddo

  enddo ! end of k loop

  if (CS%id_ZB2020u>0)   call post_data(CS%id_ZB2020u, fx, CS%diag)
  if (CS%id_ZB2020v>0)   call post_data(CS%id_ZB2020v, fy, CS%diag)
  if (CS%id_kbc>0)       call post_data(CS%id_kbc, kbc_3d, CS%diag)
  
  if (CS%id_maskT>0)     call post_data(CS%id_maskT, mask_T_3d, CS%diag)
  if (CS%id_maskq>0)     call post_data(CS%id_maskq, mask_q_3d, CS%diag)
  
  if (CS%id_S_11f>0)     call post_data(CS%id_S_11f, S_11_3df, CS%diag)

  if (CS%id_S_22f>0)     call post_data(CS%id_S_22f, S_22_3df, CS%diag)

  if (CS%id_S_12f>0)     call post_data(CS%id_S_12f, S_12_3df, CS%diag)

  call compute_energy_source(u, v, h, fx, fy, G, GV, CS)

end subroutine Zanna_Bolton_2020

! if n_lowpass and n_highpass are positive,
! performs n_lowpass iterations of 
! filter of order 2*n_highpass
! if n_lowpass is negative, returns residual instead
! Input does not require halo
! Output has full halo
! filtering occurs in-place
subroutine filter(G, mask_T, mask_q, n_lowpass, n_highpass, T, q)
  type(ocean_grid_type),              intent(in)              :: G          !< Ocean grid
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)              :: mask_T     !< mask of wet points in T points
  real, dimension(SZIB_(G),SZJB_(G)), intent(in)              :: mask_q     !< mask of wet points in q points
  real, dimension(SZI_(G),SZJ_(G)),   optional, intent(inout) :: T          !< any field at T points
  real, dimension(SZIB_(G),SZJB_(G)), optional, intent(inout) :: q          !< any field at q points
  integer,                            intent(in)              :: n_lowpass  !< number of low-pass iterations 
  integer,                            intent(in)              :: n_highpass !< number of high-pass iterations 

  integer :: i, j
  real, dimension(SZIB_(G),SZJB_(G)) :: q1, q2 ! additional q fields
  real, dimension(SZI_(G),SZJ_(G))   :: T1, T2 ! additional T fields
  real :: max_before, min_before, max_after, min_after    ! for testing
  
  if (n_lowpass==0) then
    return
  endif

  ! Total operator is I - (I-G^n_lowpass)^n_highpass
  if (present(q)) then
    call pass_var(q, G%Domain, position=CORNER, complete=.true.)
    q(:,:) = q(:,:) * mask_q(:,:)
    call min_max(q, min_before, max_before)

    q1(:,:) = q(:,:)
    
    do i=1,n_highpass
      q2(:,:) = q1(:,:)
      ! q2 -> (G^n_lowpass)*q2
      do j=1,ABS(n_lowpass)
        call smooth_Tq(G, mask_T, mask_q, q=q2)
      enddo
      ! q1 -> (I-G^n_lowpass)*q1
      q1(:,:) = q1(:,:) - q2(:,:)
    enddo

    if (n_lowpass>0) then
      ! q -> q - ((I-G^n_lowpass)^n_highpass)*q
      q(:,:) = q(:,:) - q1(:,:)
    else
      ! q -> ((I-G^n_lowpass)^n_highpass)*q
      q(:,:) = q1(:,:)
    endif
    
    !call check_nan(q, 'applying filter at q points')

    if (n_highpass==1 .AND. n_lowpass>0) then
      call min_max(q, min_after, max_after)
      if (max_after > max_before .OR. min_after < min_before) then
        write(*,*) 'filter error: not monotone in q field:', min_before, min_after, max_before, max_after
      endif
    endif
  endif

  if (present(T)) then
    call pass_var(T, G%Domain)
    T(:,:) = T(:,:) * mask_T(:,:)
    call min_max(T, min_before, max_before)

    T1(:,:) = T(:,:)
    
    do i=1,n_highpass
      T2(:,:) = T1(:,:)
      do j=1,ABS(n_lowpass)
        call smooth_Tq(G, mask_T, mask_q, T=T2)
      enddo
      T1(:,:) = T1(:,:) - T2(:,:)
    enddo

    if (n_lowpass>0) then
      T(:,:) = T(:,:) - T1(:,:)
    else
      T(:,:) = T1(:,:)
    endif
    
    !call check_nan(T, 'applying filter at T points')

    if (n_highpass==1 .AND. n_lowpass>0) then
      call min_max(T, min_after, max_after)
      if (max_after > max_before .OR. min_after < min_before) then
        write(*,*) 'filter error: not monotone in T field:', min_before, min_after, max_before, max_after
      endif
    endif
  endif
end subroutine filter

! returns filtered fields in-place and 
! residuals as optional argument
subroutine smooth_Tq(G, mask_T, mask_q, T, q)
  type(ocean_grid_type),              intent(in)              :: G      !< Ocean grid
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)              :: mask_T !< mask of wet points in T points
  real, dimension(SZIB_(G),SZJB_(G)), intent(in)              :: mask_q !< mask of wet points in q points
  real, dimension(SZI_(G),SZJ_(G)),   optional, intent(inout) :: T      !< any field at T points
  real, dimension(SZIB_(G),SZJB_(G)), optional, intent(inout) :: q      !< any field at q points
  
  real, dimension(SZI_(G),SZJ_(G))   :: Tim ! intermediate value of T-field
  real, dimension(SZIB_(G),SZJB_(G)) :: qim ! intermediate value of q-field

  integer :: i, j
  real :: wside   ! weights for side (i+1,j), (i-1,j), (i,j+1), (i,j-1)
  real :: wcorner ! weights for corners (i+1,j+1), (i+1,j-1), (i-1,j-1), (i-1,j+1)
  real :: wcenter ! weight for center point (i,j)
  
  wside = 1. / 8.
  wcorner = 1. / 16.
  wcenter = 1. - (wside*4. + wcorner*4.)

  if (present(q)) then
    call pass_var(q, G%Domain, position=CORNER, complete=.true.)
    qim(:,:) = q(:,:) * mask_q(:,:)
    do J = G%JscB, G%JecB
      do I = G%IscB, G%IecB
        q(I,J) = wcenter * qim(i,j)                &
               + wcorner * (                       &
                  (qim(I-1,J-1)+qim(I+1,J+1))      &
                + (qim(I-1,J+1)+qim(I+1,J-1))      &
               )                                   &
               + wside * (                         &
                  (qim(I-1,J)+qim(I+1,J))          &
                + (qim(I,J-1)+qim(I,J+1))          &
               )
        q(I,J) = q(I,J) * mask_q(I,J)
      enddo
    enddo
    call pass_var(q, G%Domain, position=CORNER, complete=.true.)
  endif

  if (present(T)) then
    call pass_var(T, G%Domain)
    Tim(:,:) = T(:,:) * mask_T(:,:)
    do j = G%jsc, G%jec
      do i = G%isc, G%iec
        T(i,j) = wcenter * Tim(i,j)                &
               + wcorner * (                       &
                  (Tim(i-1,j-1)+Tim(i+1,j+1))      &
                + (Tim(i-1,j+1)+Tim(i+1,j-1))      &
               )                                   &
               + wside * (                         &
                  (Tim(i-1,j)+Tim(i+1,j))          &
                + (Tim(i,j-1)+Tim(i,j+1))          &
               )
        T(i,j) = T(i,j) * mask_T(i,j)
      enddo
    enddo
    call pass_var(T, G%Domain)
  endif

end subroutine smooth_Tq

subroutine min_max(array, min_val, max_val)
  real, dimension(:,:), intent(in) :: array
  real, intent(out)                :: min_val, max_val
  
  min_val = minval(array)
  max_val = maxval(array)
  call min_across_PEs(min_val)
  call max_across_PEs(max_val)
end subroutine

subroutine check_nan(array, str)
  !use mpi
  use MOM_coms_infra, only : PE_here
  real, intent(in) :: array(:,:)
  character(*), intent(in) :: str

  integer :: i,j
  integer :: nx, ny
  logical :: flag = .False.
  integer :: ierr
  character(100) :: out

  nx = size(array,1)
  ny = size(array,2)

  do i=1,nx
    do j=1,ny
      if (isnan(array(i,j))) then
        write(*,'(A,I3,A,I3,A,I3,A,I3,A,I3,A,A)') 'NaN at(',i,',',j,') of (',nx,',',ny,') at PE', PE_here(), ', in ', str 
        flag = .True.
      endif
    enddo
  enddo

end subroutine check_nan

subroutine set_chessnoise(G, T, q)
  type(ocean_grid_type),              intent(in)              :: G      !< Ocean grid
  real, dimension(SZI_(G),SZJ_(G)),   intent(inout) :: T      !< any field at T points
  real, dimension(SZIB_(G),SZJB_(G)), intent(inout) :: q      !< any field at q points

  integer :: i, j, ig, jg

  do j = G%jsd, G%jed
    do i = G%isd, G%ied
      ig = i + G%idg_offset
      jg = j + G%jdg_offset
      T(i,j) = (-1.) ** (ig+jg)
    enddo
  enddo

  do J = G%JsdB, G%JedB
    do I = G%IsdB, G%IedB
      Ig = I + G%idg_offset
      Jg = J + G%jdg_offset
      q(I,J) = (-1.) ** (Ig+Jg)
    enddo
  enddo

end subroutine set_chessnoise

subroutine compute_masks(G, GV, h, mask_T, mask_q, k)
  type(ocean_grid_type),              intent(in)    :: G      !< Ocean grid
  type(verticalGrid_type),            intent(in)    :: GV     !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),   intent(inout) :: mask_T !< mask of wet points in T points
  real, dimension(SZIB_(G),SZJB_(G)), intent(inout) :: mask_q !< mask of wet points in q points
  integer,                            intent(in)    :: k      !< index of vertical layer

  real :: hmin
  integer :: i, j

  hmin = GV%Angstrom_H * 2. ! min thickness beyound which we have boundary

  mask_q(:,:) = 0.
  do J = G%JscB, G%JecB
    do I = G%IscB, G%IecB
      if (h(i+1,j+1,k) < hmin .or. &
          h(i  ,j  ,k) < hmin .or. &
          h(i+1,j  ,k) < hmin .or. &
          h(i  ,j+1,k) < hmin      &
        ) then
          mask_q(I,J) = 0.
        else
          mask_q(I,J) = 1.
        endif
        mask_q(I,J) = mask_q(I,J) * G%mask2dBu(I,J)
    enddo
  enddo
  call pass_var(mask_q, G%Domain, position=CORNER, complete=.true.)

  mask_T(:,:) = 0.
  do j = G%jsc, G%jec
    do i = G%isc, G%iec
      if (h(i,j,k) < hmin) then
        mask_T(i,j) = 0.
      else
        mask_T(i,j) = 1.
      endif
        mask_T(i,j) = mask_T(i,j) * G%mask2dT(i,j)
    enddo
  enddo
  call pass_var(mask_T, G%Domain)

end subroutine compute_masks

! This is copy-paste from MOM_diagnostics.F90, specifically 1125 line
subroutine compute_energy_source(u, v, h, fx, fy, G, GV, CS)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  type(ZB2020_CS),               intent(in)  :: CS     !< ZB2020 control structure.

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2].
  
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in) :: fx     !< Zonal acceleration due to convergence of
                                                      !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in) :: fy     !< Meridional acceleration due to convergence
                                                      !! of along-coordinate stress tensor [L T-2 ~> m s-2]
  
  real :: KE_term(SZI_(G),SZJ_(G),SZK_(GV)) ! A term in the kinetic energy budget
                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]
  real :: tmp(SZI_(G),SZJ_(G),SZK_(GV)) ! temporary array for integration
  real :: KE_u(SZIB_(G),SZJ_(G)) ! The area integral of a KE term in a layer at u-points
                                 ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real :: KE_v(SZI_(G),SZJB_(G)) ! The area integral of a KE term in a layer at v-points
                                 ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]

  type(group_pass_type) :: pass_KE_uv !< A handle used for group halo passes

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k

  real :: uh !< Transport through zonal faces = u*h*dy,
             !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: vh !< Transport through meridional faces = v*h*dx,
             !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: global_integral !< Global integral of the energy effect of ZB2020 [W]

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  
  if (CS%id_KE_ZB2020 > 0) then
    call create_group_pass(pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)
    
    KE_term(:,:,:) = 0.
    tmp(:,:,:) = 0.
    ! Calculate the KE source from Zanna-Bolton2020 [H L2 T-3 ~> m3 s-3].
    do k=1,nz
      KE_u(:,:) = 0.
      KE_v(:,:) = 0.
      do j=js,je ; do I=Isq,Ieq
        uh = u(I,j,k) * 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i+1,j)*h(i+1,j,k)) * &
          G%dyCu(I,j)
        KE_u(I,j) = uh * G%dxCu(I,j) * fx(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        vh = v(i,J,k) * 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i,j+1)*h(i,j+1,k)) * &
          G%dxCv(i,J)
        KE_v(i,J) = vh * G%dyCv(i,J) * fy(i,J,k)
      enddo ; enddo
      call do_group_pass(pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
        ! copy-paste from MOM_spatial_means.F90, line 42
        tmp(i,j,k) = KE_term(i,j,k) * G%areaT(i,j) * G%mask2dT(i,j)
      enddo ; enddo
    enddo

    global_integral = reproducing_sum(tmp)

    !write(*,*) 'Global energy rate of change [W] for ZB2020:', global_integral

    call post_data(CS%id_KE_ZB2020, KE_term, CS%diag)
  endif

end subroutine compute_energy_source

end module MOM_Zanna_Bolton