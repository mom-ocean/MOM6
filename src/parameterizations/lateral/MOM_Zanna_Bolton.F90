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
use MOM_error_handler, only : MOM_error, WARNING

implicit none ; private

#include <MOM_memory.h>

public Zanna_Bolton_2020, ZB_2020_init

!> Control structure for Zanna-Bolton-2020 parameterization.
type, public :: ZB2020_CS ; private
  ! Parameters
  real      :: amplitude      !< The nondimensional scaling factor in ZB model,
                              !! typically 0.1 - 10 [nondim].
  integer   :: ZB_type        !< Select how to compute the trace part of ZB model:
                              !! 0 - both deviatoric and trace components are computed
                              !! 1 - only deviatoric component is computed
                              !! 2 - only trace component is computed
  integer   :: ZB_cons        !< Select a discretization scheme for ZB model
                              !! 0 - non-conservative scheme
                              !! 1 - conservative scheme for deviatoric component
  integer   :: LPF_iter       !< Number of smoothing passes for the Velocity Gradient (VG) components
                              !! in ZB model.
  integer   :: LPF_order      !< The scale selectivity of the smoothing filter
                              !! 1 - Laplacian filter
                              !! 2 - Bilaplacian filter
  integer   :: HPF_iter       !< Number of sharpening passes for the Velocity Gradient (VG) components
                              !! in ZB model.
  integer   :: HPF_order      !< The scale selectivity of the sharpening filter
                              !! 1 - Laplacian filter
                              !! 2 - Bilaplacian filter
  integer   :: Stress_iter    !< Number of smoothing passes for the Stress tensor components
                              !! in ZB model.
  integer   :: Stress_order   !< The scale selectivity of the smoothing filter
                              !! 1 - Laplacian filter
                              !! 2 - Bilaplacian filter
  integer   :: ssd_iter       !< Hyperviscosity parameter. Defines the number of sharpening passes
                              !! in Laplacian viscosity model:
                              !! -1: hyperviscosity is off
                              !!  0: Laplacian viscosity
                              !!  9: (Laplacian)^10 viscosity, ...
  real      :: ssd_bound_coef !< The non-dimensional damping coefficient of the grid harmonic
                              !! by hyperviscous dissipation:
                              !! 0.0: no damping
                              !! 1.0: grid harmonic is removed after a step in time
  real      :: DT             !< The (baroclinic) dynamics time step [T ~> s]

  type(diag_ctrl), pointer :: diag => NULL() !< A type that regulates diagnostics output
  !>@{ Diagnostic handles
  integer :: id_ZB2020u = -1, id_ZB2020v = -1, id_KE_ZB2020 = -1
  integer :: id_maskT = -1
  integer :: id_maskq = -1
  integer :: id_S_11 = -1
  integer :: id_S_22 = -1
  integer :: id_S_12 = -1
  !>@}

end type ZB2020_CS

contains

!> Read parameters and register output fields
!! used in Zanna_Bolton_2020().
subroutine ZB_2020_init(Time, GV, US, param_file, diag, CS, use_ZB2020)
  type(time_type),         intent(in)    :: Time       !< The current model time.
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file parser structure.
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(ZB2020_CS),         intent(inout) :: CS         !< ZB2020 control structure.
  logical,                 intent(out)   :: use_ZB2020 !< If true, turns on ZB scheme.

  ! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_Zanna_Bolton" ! This module's name.

  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "USE_ZB2020", use_ZB2020, &
                 "If true, turns on Zanna-Bolton-2020 (ZB) " //&
                 "subgrid momentum parameterization of mesoscale eddies.", default=.false.)
  if (.not. use_ZB2020) return

  call get_param(param_file, mdl, "ZB_SCALING", CS%amplitude, &
                 "The nondimensional scaling factor in ZB model, " //&
                 "typically 0.1 - 10.", units="nondim", default=0.3)

  call get_param(param_file, mdl, "ZB_TRACE_MODE", CS%ZB_type, &
                 "Select how to compute the trace part of ZB model:\n" //&
                 "\t 0 - both deviatoric and trace components are computed\n" //&
                 "\t 1 - only deviatoric component is computed\n" //&
                 "\t 2 - only trace component is computed", default=0)

  call get_param(param_file, mdl, "ZB_SCHEME", CS%ZB_cons, &
                 "Select a discretization scheme for ZB model:\n" //&
                 "\t 0 - non-conservative scheme\n" //&
                 "\t 1 - conservative scheme for deviatoric component", default=1)

  call get_param(param_file, mdl, "VG_SMOOTH_PASS", CS%LPF_iter, &
                 "Number of smoothing passes for the Velocity Gradient (VG) components " //&
                 "in ZB model.", default=0)

  call get_param(param_file, mdl, "VG_SMOOTH_SEL", CS%LPF_order, &
                 "The scale selectivity of the smoothing filter " //&
                 "for VG components:\n" //&
                 "\t 1 - Laplacian filter\n" //&
                 "\t 2 - Bilaplacian filter, ...", &
                 default=1, do_not_log = CS%LPF_iter==0)

  call get_param(param_file, mdl, "VG_SHARP_PASS", CS%HPF_iter, &
                "Number of sharpening passes for the Velocity Gradient (VG) components " //&
                "in ZB model.", default=0)

  call get_param(param_file, mdl, "VG_SHARP_SEL", CS%HPF_order, &
                "The scale selectivity of the sharpening filter " //&
                "for VG components:\n" //&
                "\t 1 - Laplacian filter\n" //&
                "\t 2 - Bilaplacian filter,...", &
                default=1, do_not_log = CS%HPF_iter==0)

  call get_param(param_file, mdl, "STRESS_SMOOTH_PASS", CS%Stress_iter, &
                 "Number of smoothing passes for the Stress tensor components " //&
                 "in ZB model.", default=0)

  call get_param(param_file, mdl, "STRESS_SMOOTH_SEL", CS%Stress_order, &
                "The scale selectivity of the smoothing filter " //&
                "for the Stress tensor components:\n" //&
                "\t 1 - Laplacian filter\n" //&
                "\t 2 - Bilaplacian filter,...", &
                default=1, do_not_log = CS%Stress_iter==0)

  call get_param(param_file, mdl, "ZB_HYPERVISC", CS%ssd_iter, &
                 "Select an additional hyperviscosity to stabilize the ZB model:\n" //&
                 "\t 0  - off\n" //&
                 "\t 1  - Laplacian viscosity\n" //&
                 "\t 10 - (Laplacian)**10 viscosity, ...", &
                 default=0)
                 ! Convert to the number of sharpening passes
                 ! applied to the Laplacian viscosity model
                 CS%ssd_iter = CS%ssd_iter-1

  call get_param(param_file, mdl, "HYPVISC_GRID_DAMP", CS%ssd_bound_coef, &
                 "The non-dimensional damping coefficient of the grid harmonic " //&
                 "by hyperviscous dissipation:\n" //&
                 "\t 0.0 - no damping\n" //&
                 "\t 1.0 - grid harmonic is removed after a step in time", &
                 units="nondim", default=0.2, do_not_log = CS%ssd_iter==-1)

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

  CS%id_maskT  = register_diag_field('ocean_model', 'maskT', diag%axesTL, Time, &
      'Mask of wet points in T (CENTER) points', '1', conversion=1.)

  CS%id_maskq  = register_diag_field('ocean_model', 'maskq', diag%axesBL, Time, &
      'Mask of wet points in q (CORNER) points', '1', conversion=1.)

  ! action of filter on momentum flux
  CS%id_S_11 = register_diag_field('ocean_model', 'S_11', diag%axesTL, Time, &
      'Diagonal term (11) in the ZB stress tensor', 'm2s-2', conversion=US%L_T_to_m_s**2)

  CS%id_S_22 = register_diag_field('ocean_model', 'S_22', diag%axesTL, Time, &
      'Diagonal term (22) in the ZB stress tensor', 'm2s-2', conversion=US%L_T_to_m_s**2)

  CS%id_S_12 = register_diag_field('ocean_model', 'S_12', diag%axesBL, Time, &
      'Off-diagonal term in the ZB stress tensor', 'm2s-2', conversion=US%L_T_to_m_s**2)

end subroutine ZB_2020_init

!> Baroclinic Zanna-Bolton-2020 parameterization, see
!! eq. 6 in https://laurezanna.github.io/files/Zanna-Bolton-2020.pdf
!! We collect all contributions to a tensor S, with components:
!! (S_11, S_12;
!!  S_12, S_22)
!! Which consists of the deviatoric and trace components, respectively:
!! S =   (-vort_xy * sh_xy, vort_xy * sh_xx;
!!         vort_xy * sh_xx,  vort_xy * sh_xy) +
!! 1/2 * (vort_xy^2 + sh_xy^2 + sh_xx^2, 0;
!!        0, vort_xy^2 + sh_xy^2 + sh_xx^2)
!! Where:
!! vort_xy = dv/dx - du/dy - relative vorticity
!! sh_xy   = dv/dx + du/dy - shearing deformation (or horizontal shear strain)
!! sh_xx   = du/dx - dv/dy - stretching deformation (or horizontal tension)
!! Update of the governing equations:
!! (du/dt, dv/dt) = k_BC * div(S)
!! Where:
!! k_BC = - amplitude * grid_cell_area
!! amplitude = 0.1..10 (approx)

subroutine Zanna_Bolton_2020(u, v, h, fx, fy, G, GV, CS)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  type(ZB2020_CS),               intent(in)  :: CS     !< ZB2020 control structure.

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                                 intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2].

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(out) :: fx     !< Zonal acceleration due to convergence of
                                                       !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(out) :: fy     !< Meridional acceleration due to convergence
                                                       !! of along-coordinate stress tensor [L T-2 ~> m s-2]

  ! Arrays defined in h (CENTER) points
  real, dimension(SZI_(G),SZJ_(G)) :: &
    dx_dyT, &          ! dx/dy at h points [nondim]
    dy_dxT, &          ! dy/dx at h points [nondim]
    dx2h, &            ! dx^2 at h points [L2 ~> m2]
    dy2h, &            ! dy^2 at h points [L2 ~> m2]
    dudx, dvdy, &      ! Components in the horizontal tension [T-1 ~> s-1]
    sh_xx, &           ! Horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    vort_xy_center, &  ! Vorticity interpolated to the center [T-1 ~> s-1]
    sh_xy_center, &    ! Shearing strain interpolated to the center [T-1 ~> s-1]
    S_11, S_22, &      ! Diagonal terms in the ZB stress tensor:
                       ! Above Line 539 [L2 T-2 ~> m2 s-2]
                       ! Below Line 539 it is layer-integrated [H L2 T-2 ~> m3 s-2 or kg s-2]
    ssd_11, &          ! Diagonal component of hyperviscous stress [L2 T-2 ~> m2 s-2]
    ssd_11_coef, &     ! Viscosity coefficient in hyperviscous stress in center points
                       ! [L2 T-1 ~> m2 s-1]
    mask_T             ! Mask of wet points in T (CENTER) points [nondim]

  ! Arrays defined in q (CORNER) points
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dx_dyBu, &         ! dx/dy at q points [nondim]
    dy_dxBu, &         ! dy/dx at q points [nondim]
    dx2q, &            ! dx^2 at q points [L2 ~> m2]
    dy2q, &            ! dy^2 at q points [L2 ~> m2]
    dvdx, dudy, &      ! Components in the shearing strain [T-1 ~> s-1]
    vort_xy, &         ! Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
    sh_xy, &           ! Horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
    sh_xx_corner, &    ! Horizontal tension interpolated to the corner [T-1 ~> s-1]
    S_12, &            ! Off-diagonal term in the ZB stress tensor:
                       ! Above Line 539 [L2 T-2 ~> m2 s-2]
                       ! Below Line 539 it is layer-integrated [H L2 T-2 ~> m3 s-2 or kg s-2]
    ssd_12, &          ! Off-diagonal component of hyperviscous stress [L2 T-2 ~> m2 s-2]
    ssd_12_coef, &     ! Viscosity coefficient in hyperviscous stress in corner points
                       ! [L2 T-1 ~> m2 s-1]
    mask_q             ! Mask of wet points in q (CORNER) points [nondim]

  ! Thickness arrays for computing the horizontal divergence of the stress tensor
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    hq                 ! Thickness in CORNER points [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G))  :: &
    h_u                ! Thickness interpolated to u points [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJB_(G))  :: &
    h_v                ! Thickness interpolated to v points [H ~> m or kg m-2].

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    mask_T_3d, &       ! Mask of wet points in T (CENTER) points [nondim]
    S_11_3d, S_22_3d   ! Diagonal terms in the ZB stress tensor [L2 T-2 ~> m2 s-2]

  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)) :: &
    mask_q_3d, &       ! Mask of wet points in q (CORNER) points [nondim]
    S_12_3d            ! Off-diagonal term in the ZB stress tensor [L2 T-2 ~> m2 s-2]

  real :: h_neglect    ! Thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
  real :: h_neglect3   ! h_neglect^3 [H3 ~> m3 or kg3 m-6]
  real :: h2uq, h2vq   ! Temporary variables [H2 ~> m2 or kg2 m-4].

  real :: sum_sq       ! 1/2*(vort_xy^2 + sh_xy^2 + sh_xx^2) [T-2 ~> s-2]
  real :: vort_sh      ! vort_xy*sh_xy [T-2 ~> s-2]

  real :: k_bc         ! Constant in from of the parameterization [L2 ~> m2]
                       ! Related to the amplitude as follows:
                       ! k_bc = - amplitude * grid_cell_area < 0

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n

  ! Line 407 of MOM_hor_visc.F90
  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  h_neglect  = GV%H_subroundoff ! Line 410 on MOM_hor_visc.F90
  h_neglect3 = h_neglect**3

  fx(:,:,:) = 0.
  fy(:,:,:) = 0.

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

  if (CS%ssd_iter > -1) then
    ssd_11_coef(:,:) = 0.
    ssd_12_coef(:,:) = 0.
    do J=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      ssd_11_coef(i,j) = ((CS%ssd_bound_coef * 0.25) / CS%DT) &
                       * ((dx2h(i,j) * dy2h(i,j)) / (dx2h(i,j) + dy2h(i,j)))
    enddo; enddo

    do J=js-1,Jeq ; do I=is-1,Ieq
      ssd_12_coef(I,J) = ((CS%ssd_bound_coef * 0.25) / CS%DT) &
                       * ((dx2q(I,J) * dy2q(I,J)) / (dx2q(I,J) + dy2q(I,J)))
    enddo; enddo
  endif

  do k=1,nz

    sh_xx(:,:) = 0.
    sh_xy(:,:) = 0.
    vort_xy(:,:) = 0.
    S_12(:,:) = 0.
    S_11(:,:) = 0.
    S_22(:,:) = 0.
    ssd_11(:,:) = 0.
    ssd_12(:,:) = 0.

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
    if (CS%id_maskT>0) then
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        mask_T_3d(i,j,k) = mask_T(i,j)
      enddo; enddo
    endif

    if (CS%id_maskq>0) then
      do J=Jsq-1,Jeq+1 ; do I=Isq-1,Ieq+1
        mask_q_3d(i,j,k) = mask_q(i,j)
      enddo; enddo
    endif

    ! Numerical scheme for ZB2020 requires
    ! interpolation center <-> corner
    ! This interpolation requires B.C.,
    ! and that is why B.C. for Velocity Gradients should be
    ! well defined
    ! The same B.C. will be used by all filtering operators
    do J=Jsq-1,Jeq+2 ; do I=Isq-1,Ieq+2
      sh_xx(i,j) = sh_xx(i,j) * mask_T(i,j)
    enddo ; enddo

    do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
      sh_xy(i,j) = sh_xy(i,j) * mask_q(i,j)
      vort_xy(i,j) = vort_xy(i,j) * mask_q(i,j)
    enddo ; enddo

    if (CS%ssd_iter > -1) then
      do J=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        ssd_11(i,j) = sh_xx(i,j) * ssd_11_coef(i,j)
      enddo; enddo

      do J=js-1,Jeq ; do I=is-1,Ieq
        ssd_12(I,J) = sh_xy(I,J) * ssd_12_coef(I,J)
      enddo; enddo

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
      k_bc = - CS%amplitude * G%areaT(i,j)
      S_11(i,j) = k_bc * (- vort_sh + sum_sq)
      S_22(i,j) = k_bc * (+ vort_sh + sum_sq)
    enddo ; enddo

    ! Form S_12 tensor
    ! indices correspond to sh_xx_corner loop
    do J=Jsq-1,Jeq ; do I=Isq-1,Ieq
      if (CS%ZB_type == 2) then
        vort_sh = 0.
      else
        vort_sh = vort_xy(I,J) * sh_xx_corner(I,J)
      endif
      k_bc = - CS%amplitude * G%areaBu(i,j)
      S_12(I,J) = k_bc * vort_sh
    enddo ; enddo

    call filter(G, mask_T, mask_q, CS%Stress_iter, CS%Stress_order, T=S_11)
    call filter(G, mask_T, mask_q, CS%Stress_iter, CS%Stress_order, T=S_22)
    call filter(G, mask_T, mask_q, CS%Stress_iter, CS%Stress_order, q=S_12)

    if (CS%ssd_iter>-1) then
      do J=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        S_11(i,j) = S_11(i,j) + ssd_11(i,j)
        S_22(i,j) = S_22(i,j) - ssd_11(i,j)
      enddo ; enddo
      do J=js-1,Jeq ; do I=is-1,Ieq
        S_12(I,J) = S_12(I,J) + ssd_12(I,J)
      enddo ; enddo
    endif

    if (CS%id_S_11>0) then
      do J=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        S_11_3d(i,j,k) = S_11(i,j)
      enddo; enddo
    endif

    if (CS%id_S_22>0) then
      do J=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        S_22_3d(i,j,k) = S_22(i,j)
      enddo; enddo
    endif

    if (CS%id_S_12>0) then
      do J=js-1,Jeq ; do I=is-1,Ieq
        S_12_3d(I,J,k) = S_12(I,J)
      enddo; enddo
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

  if (CS%id_maskT>0)     call post_data(CS%id_maskT, mask_T_3d, CS%diag)
  if (CS%id_maskq>0)     call post_data(CS%id_maskq, mask_q_3d, CS%diag)

  if (CS%id_S_11>0)     call post_data(CS%id_S_11, S_11_3d, CS%diag)

  if (CS%id_S_22>0)     call post_data(CS%id_S_22, S_22_3d, CS%diag)

  if (CS%id_S_12>0)     call post_data(CS%id_S_12, S_12_3d, CS%diag)

  call compute_energy_source(u, v, h, fx, fy, G, GV, CS)

end subroutine Zanna_Bolton_2020

!> Filter which is used to smooth velocity gradient tensor
!! or the stress tensor.
!! If n_lowpass and n_highpass are positive,
!! the filter is given by:
!! I - (I-G^n_lowpass)^n_highpass
!! where I is the identity matrix and G is smooth_Tq().
!! It is filter of order 2*n_highpass,
!! where n_lowpass is the number of iterations
!! which defines the filter scale.
!! If n_lowpass is negative, returns residual
!! for the same filter:
!! (I-G^|n_lowpass|)^n_highpass
!! Input does not require halo. Output has full halo.
subroutine filter(G, mask_T, mask_q, n_lowpass, n_highpass, T, q)
  type(ocean_grid_type), intent(in) :: G !< Ocean grid
  integer, intent(in) :: n_lowpass  !< number of low-pass iterations
  integer, intent(in) :: n_highpass !< number of high-pass iterations
  real, dimension(SZI_(G),SZJ_(G)),   &
                          intent(in) :: mask_T !< mask of wet points in T (CENTER) points [nondim]
  real, dimension(SZIB_(G),SZJB_(G)), &
                          intent(in) :: mask_q !< mask of wet points in q (CORNER) points [nondim]
  real, dimension(SZI_(G),SZJ_(G)),   &
             optional, intent(inout) :: T      !< any field at T (CENTER) points [arbitrary]
  real, dimension(SZIB_(G),SZJB_(G)), &
             optional, intent(inout) :: q      !< any field at q (CORNER) points [arbitrary]

  real, dimension(SZIB_(G),SZJB_(G)) :: q1, q2          ! intermediate q-fields [arbitrary]
  real, dimension(SZI_(G),SZJ_(G))   :: T1, T2          ! intermediate T-fields [arbitrary]
  real :: max_before, min_before, max_after, min_after  ! minimum and maximum values of fields
                                                        ! before and after filtering [arbitrary]

  integer :: i_highpass, i_lowpass
  integer :: i, j
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (n_lowpass==0) then
    return
  endif

  ! Total operator is I - (I-G^n_lowpass)^n_highpass
  if (present(q)) then
    call pass_var(q, G%Domain, position=CORNER, complete=.true.)
    do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
      q(I,J) = q(I,J) * mask_q(I,J)
    enddo ; enddo

    if (n_highpass==1 .AND. n_lowpass>0) then
      call min_max(G, min_before, max_before, q=q)
    endif

    do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
      q1(I,J) = q(I,J)
    enddo ; enddo

    ! q1 -> ((I-G^n_lowpass)^n_highpass)*q1
    do i_highpass=1,n_highpass
      do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
        q2(I,J) = q1(I,J)
      enddo ; enddo
      ! q2 -> (G^n_lowpass)*q2
      do i_lowpass=1,ABS(n_lowpass)
        call smooth_Tq(G, mask_T, mask_q, q=q2)
      enddo
      ! q1 -> (I-G^n_lowpass)*q1
      do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
        q1(I,J) = q1(I,J) - q2(I,J)
      enddo ; enddo
    enddo

    if (n_lowpass>0) then
      ! q -> q - ((I-G^n_lowpass)^n_highpass)*q
      do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
        q(I,J) = q(I,J) - q1(I,J)
      enddo ; enddo
    else
      ! q -> ((I-G^n_lowpass)^n_highpass)*q
      do J=Jsq-2,Jeq+2 ; do I=Isq-2,Ieq+2
        q(I,J) = q1(I,J)
      enddo ; enddo
    endif

    if (n_highpass==1 .AND. n_lowpass>0) then
      call min_max(G, min_after, max_after, q=q)
      if (max_after > max_before .OR. min_after < min_before) then
        call MOM_error(WARNING, "MOM_Zanna_Bolton.F90, filter applied in CORNER points "//&
                                "does not preserve [min,max] values. There may be issues with "//&
                                "boundary conditions")
      endif
    endif
  endif

  if (present(T)) then
    call pass_var(T, G%Domain)
    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      T(i,j) = T(i,j) * mask_T(i,j)
    enddo ; enddo

    if (n_highpass==1 .AND. n_lowpass>0) then
      call min_max(G, min_before, max_before, T=T)
    endif

    do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
      T1(i,j) = T(i,j)
    enddo ; enddo

    do i_highpass=1,n_highpass
      do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
        T2(i,j) = T1(i,j)
      enddo ; enddo
      do i_lowpass=1,ABS(n_lowpass)
        call smooth_Tq(G, mask_T, mask_q, T=T2)
      enddo
      do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
        T1(i,j) = T1(i,j) - T2(i,j)
      enddo ; enddo
    enddo

    if (n_lowpass>0) then
      do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
        T(i,j) = T(i,j) - T1(i,j)
      enddo ; enddo
    else
      do j=Jsq-1,Jeq+2 ; do i=Isq-1,Ieq+2
        T(i,j) = T1(i,j)
      enddo ; enddo
    endif

    if (n_highpass==1 .AND. n_lowpass>0) then
      call min_max(G, min_after, max_after, T=T)
      if (max_after > max_before .OR. min_after < min_before) then
        call MOM_error(WARNING, "MOM_Zanna_Bolton.F90, filter applied in CENTER points "//&
                                " does not preserve [min,max] values. There may be issues with "//&
                                " boundary conditions")
      endif
    endif
  endif
end subroutine filter

!> One iteration of 3x3 filter
!! [1 2 1;
!!  2 4 2;
!!  1 2 1]/16
!! removing chess-harmonic.
!! It is used as a buiding block in filter().
!! Zero Dirichlet boundary conditions are applied
!! with mask_T and mask_q.
subroutine smooth_Tq(G, mask_T, mask_q, T, q)
  type(ocean_grid_type), intent(in) :: G !< Ocean grid
  real, dimension(SZI_(G),SZJ_(G)),   &
                          intent(in) :: mask_T !< mask of wet points in T (CENTER) points [nondim]
  real, dimension(SZIB_(G),SZJB_(G)), &
                          intent(in) :: mask_q !< mask of wet points in q (CORNER) points [nondim]
  real, dimension(SZI_(G),SZJ_(G)),   &
             optional, intent(inout) :: T      !< any field at T (CENTER) points [arbitrary]
  real, dimension(SZIB_(G),SZJB_(G)), &
             optional, intent(inout) :: q      !< any field at q (CORNER) points [arbitrary]

  real, dimension(SZI_(G),SZJ_(G))   :: Tim ! intermediate T-field [arbitrary]
  real, dimension(SZIB_(G),SZJB_(G)) :: qim ! intermediate q-field [arbitrary]

  real :: wside   ! weights for side points
                  ! (i+1,j), (i-1,j), (i,j+1), (i,j-1)
                  ! [nondim]
  real :: wcorner ! weights for corner points
                  ! (i+1,j+1), (i+1,j-1), (i-1,j-1), (i-1,j+1)
                  ! [nondim]
  real :: wcenter ! weight for the center point (i,j) [nondim]

  integer :: i, j
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  wside = 1. / 8.
  wcorner = 1. / 16.
  wcenter = 1. - (wside*4. + wcorner*4.)

  if (present(q)) then
    call pass_var(q, G%Domain, position=CORNER, complete=.true.)
    do J = Jsq-1, Jeq+1; do I = Isq-1, Ieq+1
      qim(I,J) = q(I,J) * mask_q(I,J)
    enddo; enddo
    do J = Jsq, Jeq
      do I = Isq, Ieq
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
    do j = js-1, je+1; do i = is-1, ie+1
      Tim(i,j) = T(i,j) * mask_T(i,j)
    enddo; enddo
    do j = js, je
      do i = is, ie
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

!> Returns min and max values of array across all PEs.
!! It is used in filter() to check its monotonicity.
subroutine min_max(G, min_val, max_val, T, q)
  type(ocean_grid_type), intent(in) :: G   !< Ocean grid
  real, dimension(SZI_(G),SZJ_(G)),   &
             optional, intent(inout) :: T  !< any field at T (CENTER) points [arbitrary]
  real, dimension(SZIB_(G),SZJB_(G)), &
             optional, intent(inout) :: q  !< any field at q (CORNER) points [arbitrary]
  real, intent(out) :: min_val, max_val    !< min and max values of array accross PEs [arbitrary]

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (present(q)) then
    min_val = minval(q(Isq:Ieq, Jsq:Jeq))
    max_val = maxval(q(Isq:Ieq, Jsq:Jeq))
  endif

  if (present(T)) then
    min_val = minval(T(is:ie, js:je))
    max_val = maxval(T(is:ie, js:je))
  endif

  call min_across_PEs(min_val)
  call max_across_PEs(max_val)

end subroutine

!> Computes mask of wet points in T (CENTER) and q (CORNER) points.
!! Method: compare layer thicknesses with Angstrom_H.
!! Mask is computed separately for every vertical layer and
!! for every time step.
subroutine compute_masks(G, GV, h, mask_T, mask_q, k)
  type(ocean_grid_type),        intent(in)    :: G      !< Ocean grid
  type(verticalGrid_type),      intent(in)    :: GV     !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),          &
                                intent(inout) :: mask_T !< mask of wet points in T (CENTER) points [nondim]
  real, dimension(SZIB_(G),SZJB_(G)),        &
                                intent(inout) :: mask_q !< mask of wet points in q (CORNER) points [nondim]
  integer,                      intent(in)    :: k      !< index of vertical layer

  real :: hmin       ! Minimum layer thickness
                     ! beyond which we have boundary [H ~> m or kg m-2]
  integer :: i, j

  hmin = GV%Angstrom_H * 2.

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

!> Computes the 3D energy source term for the ZB2020 scheme
!! similarly to MOM_diagnostics.F90, specifically 1125 line.
subroutine compute_energy_source(u, v, h, fx, fy, G, GV, CS)
  type(ocean_grid_type),         intent(in)  :: G    !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV   !< The ocean's vertical grid structure.
  type(ZB2020_CS),               intent(in)  :: CS   !< ZB2020 control structure.

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                                 intent(inout) :: h  !< Layer thicknesses [H ~> m or kg m-2].

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in) :: fx    !< Zonal acceleration due to convergence of
                                                     !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in) :: fy    !< Meridional acceleration due to convergence
                                                     !! of along-coordinate stress tensor [L T-2 ~> m s-2]

  real :: KE_term(SZI_(G),SZJ_(G),SZK_(GV)) ! A term in the kinetic energy budget
                                            ! [H L2 T-3 ~> m3 s-3 or W m-2]
  real :: KE_u(SZIB_(G),SZJ_(G))            ! The area integral of a KE term in a layer at u-points
                                            ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real :: KE_v(SZI_(G),SZJB_(G))            ! The area integral of a KE term in a layer at v-points
                                            ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]

  !real :: tmp(SZI_(G),SZJ_(G),SZK_(GV))    ! temporary array for integration
  !real :: global_integral                  ! Global integral of the energy effect of ZB2020
                                            ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]


  real :: uh                                ! Transport through zonal faces = u*h*dy,
                                            ! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: vh                                ! Transport through meridional faces = v*h*dx,
                                            ! [H L2 T-1 ~> m3 s-1 or kg s-1].

  type(group_pass_type) :: pass_KE_uv       ! A handle used for group halo passes

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (CS%id_KE_ZB2020 > 0) then
    call create_group_pass(pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)

    KE_term(:,:,:) = 0.
    !tmp(:,:,:) = 0.
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
        !tmp(i,j,k) = KE_term(i,j,k) * G%areaT(i,j) * G%mask2dT(i,j)
      enddo ; enddo
    enddo

    !global_integral = reproducing_sum(tmp)

    call post_data(CS%id_KE_ZB2020, KE_term, CS%diag)
  endif

end subroutine compute_energy_source

end module MOM_Zanna_Bolton