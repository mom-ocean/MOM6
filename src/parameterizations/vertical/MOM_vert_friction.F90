!> Implements vertical viscosity (vertvisc)
module MOM_vert_friction

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_domains,       only : pass_var, To_All, Omit_corners
use MOM_domains,       only : pass_vector, Scalar_Pair
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : post_product_u, post_product_sum_u
use MOM_diag_mediator, only : post_product_v, post_product_sum_v
use MOM_diag_mediator, only : diag_ctrl, query_averaging_enabled
use MOM_domains,       only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,       only : To_North, To_East
use MOM_debugging,     only : uvchksum, hchksum
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,  only : mech_forcing, find_ustar
use MOM_get_input,     only : directories
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : MOM_read_data, slasher
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary, only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_PointAccel,    only : write_u_accel, write_v_accel, PointAccel_init
use MOM_PointAccel,    only : PointAccel_CS
use MOM_time_manager,  only : time_type, time_type_to_real, operator(-)
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs, vertvisc_type
use MOM_variables,     only : cont_diag_ptrs, accel_diag_ptrs
use MOM_variables,     only : ocean_internal_state
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_wave_interface, only : wave_parameters_CS
use MOM_set_visc,      only : set_v_at_u, set_u_at_v
use MOM_lateral_mixing_coeffs, only : VarMix_CS

use CVMix_kpp,         only : cvmix_kpp_composite_Gshape

implicit none ; private

#include <MOM_memory.h>

public vertvisc, vertvisc_remnant, vertvisc_coef
public vertvisc_limit_vel, vertvisc_init, vertvisc_end
public updateCFLtruncationValue
public vertFPmix

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure with parameters and memory for the MOM_vert_friction module
type, public :: vertvisc_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  real    :: Hmix            !< The mixed layer thickness [Z ~> m].
  real    :: Hmix_stress     !< The mixed layer thickness over which the wind
                             !! stress is applied with direct_stress [H ~> m or kg m-2].
  real    :: Kvml_invZ2      !< The extra vertical viscosity scale in [H Z T-1 ~> m2 s-1 or Pa s] in a
                             !! surface mixed layer with a characteristic thickness given by Hmix,
                             !! and scaling proportional to (Hmix/z)^2, where z is the distance
                             !! from the surface; this can get very large with thin layers.
  real    :: Kv              !< The interior vertical viscosity [H Z T-1 ~> m2 s-1 or Pa s].
  real    :: Hbbl            !< The static bottom boundary layer thickness [Z ~> m].
  real    :: Hbbl_gl90       !< The static bottom boundary layer thickness used for GL90 [Z ~> m].
  real    :: Kv_extra_bbl    !< An extra vertical viscosity in the bottom boundary layer of thickness
                             !! Hbbl when there is not a bottom drag law in use [H Z T-1 ~> m2 s-1 or Pa s].
  real    :: vonKar          !< The von Karman constant as used for mixed layer viscosity [nondim]

  logical :: use_GL90_in_SSW !< If true, use the GL90 parameterization in stacked shallow water mode (SSW).
                             !! The calculation of the GL90 viscosity coefficient uses the fact that in SSW
                             !! we simply have 1/N^2 = h/g^prime, where g^prime is the reduced gravity.
                             !! This identity does not generalize to non-SSW setups.
  logical :: use_GL90_N2     !< If true, use GL90 vertical viscosity coefficient that is depth-independent;
                             !! this corresponds to a kappa_GM that scales as N^2 with depth.
  real    :: kappa_gl90      !< The scalar diffusivity used in the GL90 vertical viscosity scheme
                             !! [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]
  logical :: read_kappa_gl90 !< If true, read a file containing the spatially varying kappa_gl90
  real    :: alpha_gl90      !< Coefficient used to compute a depth-independent GL90 vertical
                             !! viscosity via Kv_gl90 = alpha_gl90 * f^2. Note that the implied
                             !! Kv_gl90 corresponds to a kappa_gl90 that scales as N^2 with depth.
                             !! [H Z T ~> m2 s or kg s m-1]
  real    :: vel_underflow   !< Velocity components smaller than vel_underflow
                             !! are set to 0 [L T-1 ~> m s-1].
  logical :: CFL_based_trunc !< If true, base truncations on CFL numbers, not
                             !! absolute velocities.
  real    :: CFL_trunc       !< Velocity components will be truncated when they
                             !! are large enough that the corresponding CFL number
                             !! exceeds this value [nondim].
  real    :: CFL_report      !< The value of the CFL number that will cause the
                             !! accelerations to be reported [nondim].  CFL_report
                             !! will often equal CFL_trunc.
  real    :: truncRampTime   !< The time-scale over which to ramp up the value of
                             !! CFL_trunc from CFL_truncS to CFL_truncE [T ~> s]
  real    :: CFL_truncS      !< The start value of CFL_trunc [nondim]
  real    :: CFL_truncE      !< The end/target value of CFL_trunc [nondim]
  logical :: CFLrampingIsActivated = .false. !< True if the ramping has been initialized
  type(time_type) :: rampStartTime !< The time at which the ramping of CFL_trunc starts

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NK_INTERFACE_) :: &
    a_u                !< The u-drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NK_INTERFACE_) :: &
    a_u_gl90           !< The u-drag coefficient associated with GL90 across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_,NKMEM_) :: &
    h_u                !< The effective layer thickness at u-points [H ~> m or kg m-2].
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NK_INTERFACE_) :: &
    a_v                !< The v-drag coefficient across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NK_INTERFACE_) :: &
    a_v_gl90           !< The v-drag coefficient associated with GL90 across an interface [H T-1 ~> m s-1 or Pa s m-1]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_,NKMEM_) :: &
    h_v                !< The effective layer thickness at v-points [H ~> m or kg m-2].
  real, pointer, dimension(:,:) :: a1_shelf_u => NULL() !< The u-momentum coupling coefficient under
                           !! ice shelves [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress under shelves.
  real, pointer, dimension(:,:) :: a1_shelf_v => NULL() !< The v-momentum coupling coefficient under
                           !! ice shelves [H T-1 ~> m s-1 or Pa s m-1]. Retained to determine stress under shelves.

  logical :: split          !< If true, use the split time stepping scheme.
  logical :: bottomdraglaw  !< If true, the  bottom stress is calculated with a
                            !! drag law c_drag*|u|*u. The velocity magnitude
                            !! may be an assumed value or it may be based on the
                            !! actual velocity in the bottommost HBBL, depending
                            !! on whether linear_drag is true.
  logical :: harmonic_visc  !< If true, the harmonic mean thicknesses are used
                            !! to calculate the viscous coupling between layers
                            !! except near the bottom.  Otherwise the arithmetic
                            !! mean thickness is used except near the bottom.
  real    :: harm_BL_val    !< A scale to determine when water is in the boundary
                            !! layers based solely on harmonic mean thicknesses
                            !! for the purpose of determining the extent to which
                            !! the thicknesses used in the viscosities are upwinded [nondim].
  logical :: direct_stress  !< If true, the wind stress is distributed over the topmost Hmix_stress
                            !! of fluid, and an added mixed layer viscosity or a physically based
                            !! boundary layer turbulence parameterization is not needed for stability.
  logical :: dynamic_viscous_ML  !< If true, use the results from a dynamic
                            !! calculation, perhaps based on a bulk Richardson
                            !! number criterion, to determine the mixed layer
                            !! thickness for viscosity.
  logical :: fixed_LOTW_ML  !< If true, use a Law-of-the-wall prescription for the mixed layer
                            !! viscosity within a boundary layer that is the lesser of Hmix and the
                            !! total depth of the ocean in a column.
  logical :: apply_LOTW_floor !< If true, use a Law-of-the-wall prescription to set a lower bound
                            !! on the viscous coupling between layers within the surface boundary
                            !! layer, based the distance of interfaces from the surface.  This only
                            !! acts when there are large changes in the thicknesses of successive
                            !! layers or when the viscosity is set externally and the wind stress
                            !! has subsequently increased.
  integer :: answer_date    !< The vintage of the order of arithmetic and expressions in the viscous
                            !! calculations.  Values below 20190101 recover the answers from the end
                            !! of 2018, while higher values use expressions that do not use an
                            !! arbitrary and hard-coded maximum viscous coupling coefficient between
                            !! layers.  In non-Boussinesq cases, values below 20230601 recover a
                            !! form of the viscosity within  the mixed layer that breaks up the
                            !! magnitude of the wind stress with BULKMIXEDLAYER, DYNAMIC_VISCOUS_ML
                            !! or FIXED_DEPTH_LOTW_ML, but not LOTW_VISCOUS_ML_FLOOR.
  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  integer :: nkml           !< The number of layers in the mixed layer.
  integer, pointer :: ntrunc !< The number of times the velocity has been
                            !! truncated since the last call to write_energy.
  character(len=200) :: u_trunc_file  !< The complete path to a file in which a column of
                            !! u-accelerations are written if velocity truncations occur.
  character(len=200) :: v_trunc_file !< The complete path to a file in which a column of
                            !! v-accelerations are written if velocity truncations occur.
  logical :: StokesMixing   !< If true, do Stokes drift mixing via the Lagrangian current
                            !! (Eulerian plus Stokes drift).  False by default and set
                            !! via STOKES_MIXING_COMBINED.

  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  real, allocatable, dimension(:,:) :: kappa_gl90_2d !< 2D kappa_gl90 at h-points [L2 H Z-1 T-1 ~> m2 s-1 or Pa s]

  !>@{ Diagnostic identifiers
  integer :: id_du_dt_visc = -1, id_dv_dt_visc = -1, id_du_dt_visc_gl90 = -1, id_dv_dt_visc_gl90 = -1
  integer :: id_GLwork = -1
  integer :: id_au_vv = -1, id_av_vv = -1, id_au_gl90_vv = -1, id_av_gl90_vv = -1
  integer :: id_du_dt_str = -1, id_dv_dt_str = -1
  integer :: id_h_u = -1, id_h_v = -1, id_hML_u = -1 , id_hML_v = -1
  integer :: id_Omega_w2x = -1, id_FPtau2s  = -1 , id_FPtau2w = -1
  integer :: id_uE_h  = -1, id_vE_h  = -1
  integer :: id_uStk  = -1, id_vStk  = -1
  integer :: id_uStk0 = -1, id_vStk0 = -1
  integer :: id_uInc_h= -1, id_vInc_h= -1
  integer :: id_taux_bot = -1, id_tauy_bot = -1
  integer :: id_Kv_slow = -1, id_Kv_u = -1, id_Kv_v = -1
  integer :: id_Kv_gl90_u = -1, id_Kv_gl90_v = -1
  ! integer :: id_hf_du_dt_visc    = -1, id_hf_dv_dt_visc    = -1
  integer :: id_h_du_dt_visc    = -1, id_h_dv_dt_visc    = -1
  integer :: id_hf_du_dt_visc_2d = -1, id_hf_dv_dt_visc_2d = -1
  integer :: id_h_du_dt_str    = -1, id_h_dv_dt_str    = -1
  integer :: id_du_dt_str_visc_rem = -1, id_dv_dt_str_visc_rem = -1
  !>@}

  type(PointAccel_CS), pointer :: PointAccel_CSp => NULL() !< A pointer to the control structure
                              !! for recording accelerations leading to velocity truncations

  type(group_pass_type) :: pass_KE_uv !< A handle used for group halo passes
end type vertvisc_CS

contains

!> Add nonlocal stress increments to ui^n and vi^n.
subroutine vertFPmix(ui, vi, uold, vold, hbl_h, h, forces, dt, lpost, Cemp_NL, G, GV, US, CS, OBC, Waves)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: ui     !< Zonal velocity after vertvisc [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vi     !< Meridional velocity after vertvisc [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: uold   !< Old Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: vold   !< Old Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: hbl_h !<  boundary layer depth [H ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h       !< Layer thicknesses [H ~> m or kg m-2]
  type(mech_forcing),      intent(in) :: forces  !< A structure with the driving mechanical forces
  real,                    intent(in) :: dt      !< Time increment [T ~> s]
  real,                    intent(in) :: Cemp_NL !< empirical coefficient of non-local momentum mixing [nondim]
  logical,                 intent(in) :: lpost   !< Compute and make available FPMix diagnostics
  type(unit_scale_type),   intent(in) :: US      !< A dimensional unit scaling type
  type(vertvisc_CS),       pointer    :: CS      !< Vertical viscosity control structure
  type(ocean_OBC_type),    pointer    :: OBC     !< Open boundary condition structure
  type(wave_parameters_CS), &
                   optional, pointer  :: Waves   !< Container for wave/Stokes information

  ! local variables
  real, dimension(SZIB_(G),SZJ_(G))  :: hbl_u   !< boundary layer depth (u-pts) [H ~> m]
  real, dimension(SZI_(G),SZJB_(G))  :: hbl_v   !< boundary layer depth (v-pts) [H ~> m]
  real, dimension(SZIB_(G),SZJ_(G))  :: taux_u  !< kinematic zonal wind stress (u-pts) [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJB_(G))  :: tauy_v  !< kinematic merid wind stress (v-pts) [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(G),SZJ_(G))   :: uS0     !< surface zonal Stokes drift h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G))   :: vS0     !< surface zonal Stokes drift h-pts [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uE_u    !< zonal Eulerian u-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: uE_h    !< zonal Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vE_v    !< merid Eulerian v-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: vE_h    !< merid Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uInc_u  !< zonal Eulerian u-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: uInc_h  !< zonal Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vInc_v  !< merid Eulerian v-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: vInc_h  !< merid Eulerian h-pts [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: uStk    !< zonal Stokes Drift (h-pts) [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)) :: vStk    !< merid Stokes Drift (h-pts) [L T-1 ~> m s-1]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)+1) :: omega_tau2s !< angle stress to shear (h-pts) [rad]
  real, dimension(SZI_(G) ,SZJ_(G),SZK_(GV)+1) :: omega_tau2w !< angle stress to wind  (h-pts) [rad]
  real :: omega_tmp, omega_s2x, omega_tau2x                    !< temporary angle wrt the x axis [rad]
  real :: Irho0        !< Inverse of the mean density rescaled to [Z L-1 R-1 ~> m3 kg-1]
  real :: pi           !< ! The ratio of the circumference of a circle to its diameter [nondim]
  real :: tmp_u, tmp_v !< temporary ocean mask weights on u and v points [nondim]
  real :: fexp         !< temporary exponential function [nondim]
  real :: sigma        !< temporary normalize boundary layer coordinate [nondim]
  real :: Gat1, Gsig, dGdsig !< Shape parameters [nondim]
  real :: du, dv       !< Intermediate velocity differences [L T-1 ~> m s-1]
  real :: depth        !< Cumulative of thicknesses [H ~> m]
  integer :: b, kbld, kp1, k, nz !< band and vertical indices
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq !< horizontal indices

  is = G%isc ; ie = G%iec; js = G%jsc; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke

  pi = 4. * atan2(1.,1.)
  Irho0 = 1.0 / GV%Rho0

  call pass_var(hbl_h , G%Domain, halo=1)

  ! u-points
  do j = js,je
    do I = Isq,Ieq
      taux_u(I,j)  = forces%taux(I,j) * Irho0
      if ( (G%mask2dCu(I,j) > 0.5) ) then
        ! h to u-pts
        tmp_u  = MAX (1.0 ,(G%mask2dT(i,j) + G%mask2dT(i+1,j) ) )
        hbl_u(I,j) = ((G%mask2dT(i,j) * hbl_h(i,j)) + (G%mask2dT(i+1,j) * hbl_h(i+1,j))) / tmp_u
        depth   = 0.
        Gat1  = 0.
        do k=1, nz
          ! cell center
          depth = depth + 0.5*CS%h_u(I,j,k)
          uE_u(I,j,k) = ui(I,j,k) - waves%Us_x(I,j,k)
          if ( depth < hbl_u(I,j) )     then
            sigma = depth / hbl_u(i,j)
            ! cell bottom
            depth = depth + 0.5*CS%h_u(I,j,k)
            call cvmix_kpp_composite_Gshape(sigma,Gat1,Gsig,dGdsig)
            ! nonlocal boundary-layer increment
            uInc_u(I,j,k)  = dt * Cemp_NL * taux_u(I,j) * dGdsig / hbl_u(I,j)
            ui(I,j,k) = ui(I,j,k) + uInc_u(I,j,k)
          else
            uInc_u(I,j,k) = 0.0
          endif
        enddo
      else
        do k=1, nz
          uInc_u(I,j,k) = 0.0
        enddo
      endif
    enddo
  enddo

  ! v-points
  do J = Jsq,Jeq
    do i = is,ie
      tauy_v(i,J)  = forces%tauy(i,J) * Irho0
      if ( (G%mask2dCv(i,J) > 0.5) ) then
        ! h to v-pts
        tmp_v  = max( 1.0 ,(G%mask2dT(i,j) + G%mask2dT(i,j+1)))
        hbl_v(i,J) = (G%mask2dT(i,j) * hbl_h(i,J) + G%mask2dT(i,j+1) * hbl_h(i,j+1)) / tmp_v
        depth = 0.
        Gat1  = 0.
        do k=1, nz
          ! cell center
          depth = depth + 0.5* CS%h_v(i,J,k)
          vE_v(i,J,k) = vi(i,J,k) - waves%Us_y(i,J,k)
          if ( depth < hbl_v(i,J) )    then
            sigma = depth / hbl_v(i,J)
            ! cell bottom
            depth = depth + 0.5* CS%h_v(i,J,k)
            call cvmix_kpp_composite_Gshape(sigma,Gat1,Gsig,dGdsig)
            ! nonlocal boundary-layer increment
            vInc_v(i,J,k) = dt * Cemp_NL * tauy_v(i,J) * dGdsig / hbl_v(i,J)
            vi(i,J,k) = vi(i,J,k) + vInc_v(i,J,k)
          else
            vInc_v(i,J,k)  = 0.0
          endif
        enddo
      else
        do k=1, nz
          vInc_v(i,J,k)  = 0.0
        enddo
      endif
    enddo
  enddo

  ! Compute and store diagnostics, only during the corrector step.
  if (lpost)  then
    call pass_vector(uE_u  ,  vE_v  , G%Domain, To_All)
    call pass_vector(uInc_u, vInc_v , G%Domain, To_All)
    uStk = 0.0
    vStk = 0.0
    uS0  = 0.0
    vS0  = 0.0

    do j = js,je
      do i = is,ie
        if (G%mask2dT(i,j) > 0.5)  then
          ! u to h-pts
          tmp_u  = max( 1.0 ,(G%mask2dCu(i,j) + G%mask2dCu(i-1,j)))
          ! v to h-pts
          tmp_v  = max( 1.0 ,(G%mask2dCv(i,j) + G%mask2dCv(i,j-1)))
          do k = 1,nz
            uE_h(i,j,k)   = (G%mask2dCu(i,j) *   uE_u(i,j,k) + G%mask2dCu(i-1,j) *   uE_u(i-1,j,k)) / tmp_u
            uInc_h(i,j,k) = (G%mask2dCu(i,j) * uInc_u(i,j,k) + G%mask2dCu(i-1,j) * uInc_u(i-1,j,k)) / tmp_u
            vE_h(i,j,k)   = (G%mask2dCv(i,j) *   vE_v(i,j,k) + G%mask2dCv(i,j-1) *   vE_v(i,j-1,k)) / tmp_v
            vInc_h(i,j,k) = (G%mask2dCv(i,j) * vInc_v(i,j,k) + G%mask2dCv(i,j-1) * vInc_v(i,j-1,k)) / tmp_v
          enddo
          ! Wind, Stress and Shear align at surface
          Omega_tau2w(i,j,1) = 0.0
          Omega_tau2s(i,j,1) = 0.0
          do k = 1,nz
            kp1 = min( nz , k+1)
            du = uE_h(i,j,k) - uE_h(i,j,kp1)
            dv = vE_h(i,j,k) - vE_h(i,j,kp1)
            omega_s2x = atan2( dv , du )

            du = du + uInc_h(i,j,k) - uInc_h(i,j,kp1)
            dv = dv + vInc_h(i,j,k) - vInc_h(i,j,kp1)
            omega_tau2x = atan2( dv , du )

            omega_tmp = omega_tau2x - forces%omega_w2x(i,j)
            if ( (omega_tmp  >   pi   ) )  omega_tmp = omega_tmp - 2.*pi
            if ( (omega_tmp  < (0.-pi)) )  omega_tmp = omega_tmp + 2.*pi
            Omega_tau2w(i,j,kp1) = omega_tmp

            omega_tmp = omega_tau2x - omega_s2x
            if ( (omega_tmp  >   pi   ) )  omega_tmp = omega_tmp - 2.*pi
            if ( (omega_tmp  < (0.-pi)) )  omega_tmp = omega_tmp + 2.*pi
            Omega_tau2s(i,j,kp1) = omega_tmp

          enddo
        endif

        ! Stokes drift
        do b=1,waves%NumBands
          uS0(i,j)  = uS0(i,j) + waves%UStk_Hb(i,j,b)    ! or forces%UStkb(i,j,b)
          vS0(i,j)  = vS0(i,j) + waves%VStk_Hb(i,j,b)    ! or forces%VStkb(i,j,b)
        enddo
        depth = 0.0
        do k = 1,nz
          do b  = 1, waves%NumBands
            ! cell center
            fexp = exp(-2. * waves%WaveNum_Cen(b) * (depth+0.5*h(i,j,k)) )
            uStk(i,j,k) = uStk(i,j,k) + waves%UStk_Hb(i,j,b) * fexp
            vStk(i,j,k) = vStk(i,j,k) + waves%VStk_Hb(i,j,b) * fexp
          enddo
          ! cell bottom
          depth = depth + h(i,j,k)
        enddo
      enddo
    enddo

    ! post FPmix diagnostics
    if (CS%id_uE_h    > 0) call post_data(CS%id_uE_h     , uE_h   , CS%diag)
    if (CS%id_vE_h    > 0) call post_data(CS%id_vE_h   , vE_h   , CS%diag)
    if (CS%id_uInc_h  > 0) call post_data(CS%id_uInc_h , uInc_h , CS%diag)
    if (CS%id_vInc_h  > 0) call post_data(CS%id_vInc_h , vInc_h , CS%diag)
    if (CS%id_FPtau2s > 0) call post_data(CS%id_FPtau2s, Omega_tau2s, CS%diag)
    if (CS%id_FPtau2w > 0) call post_data(CS%id_FPtau2w, Omega_tau2w, CS%diag)
    if (CS%id_uStk0   > 0) call post_data(CS%id_uStk0  , uS0 , CS%diag)
    if (CS%id_vStk0   > 0) call post_data(CS%id_vStk0  , vS0    , CS%diag)
    if (CS%id_uStk    > 0) call post_data(CS%id_uStk   , uStk   , CS%diag)
    if (CS%id_vStk    > 0) call post_data(CS%id_vStk   , vStk   , CS%diag)
    if (CS%id_Omega_w2x > 0) call post_data(CS%id_Omega_w2x, forces%omega_w2x, CS%diag)

  endif

end subroutine vertFPmix


!> Expose loop indices to IPO for alias analysis and loop transformation.
function touch_ij(i,j) result(ij)
  integer, intent(in) :: i
    !< Inner loop index
  integer, intent(in) :: j
    !< Outer loop index
  integer:: ij
    !< Trivial operation to prevent removal during optimization

  ij = i * j
end function touch_ij


!> Compute coupling coefficient associated with vertical viscosity parameterization as in Greatbatch and Lamb
!! (1990), hereafter referred to as the GL90 vertical viscosity parameterization. This vertical viscosity scheme
!! redistributes momentum in the vertical, and is the equivalent of the Gent & McWilliams (1990) parameterization,
!! but in a TWA (thickness-weighted averaged) set of equations. The vertical viscosity coefficient nu is computed
!! from kappa_GM via thermal wind balance, and the following relation:
!! nu = kappa_GM * f^2 / N^2.
!! In the following subroutine kappa_GM is assumed either (a) constant or (b) horizontally varying. In both cases,
!! (a) and  (b), one can additionally impose an EBT structure in the vertical for kappa_GM.
!! A third possible formulation of nu is depth-independent:
!! nu = f^2 * alpha
!! The latter formulation would be equivalent to a kappa_GM that varies as N^2 with depth.
!! The vertical viscosity del_z ( nu del_z u) is applied to the momentum equation with stress-free boundary
!! conditions at the top and bottom.
!!
!! In SSW mode, we have 1/N^2 = h/g'. The coupling coefficient is therefore equal to
!! a_cpl_gl90 = nu / h = kappa_GM * f^2 / g'
!! or
!! a_cpl_gl90 = nu / h = f^2 * alpha / h

subroutine find_coupling_coef_gl90(a_cpl_gl90, hvel, do_i, z_i, G, GV, CS, VarMix, work_on_u)
  type(ocean_grid_type),                        intent(in)    :: G   !< Grid structure.
  type(verticalGrid_type),                      intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)),  intent(in)    :: hvel !< Distance between interfaces
                                                                     !! at velocity points [Z ~> m]
  logical, dimension(SZIB_(G),SZJB_(G)),        intent(in)    :: do_i !< If true, determine coupling coefficient
                                                                     !!  for a column
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), intent(in)   :: z_i !< Estimate of interface heights above the
                                                                     !! bottom, normalized by the GL90 bottom
                                                                     !! boundary layer thickness [nondim]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), intent(out) :: a_cpl_gl90 !< Coupling coefficient associated
                                                                     !! with GL90 across interfaces; is not
                                                                     !! included in a_cpl [H T-1 ~> m s-1 or Pa s m-1].
  type(vertvisc_cs),                            intent(in)    :: CS  !< Vertical viscosity control structure
  type(VarMix_CS),                              intent(in)    :: VarMix !< Variable mixing coefficients
  logical,                                      intent(in)    :: work_on_u !< If true, u-points are being calculated,
                                                                     !! otherwise they are v-points.

  ! local variables
  logical :: kdgl90_use_vert_struct  ! use vertical structure for GL90 coefficient
  integer :: i, j, k, is, ie, js, je, nz
  real    :: f2         !< Squared Coriolis parameter at a velocity grid point [T-2 ~> s-2].
  real    :: h_neglect  ! A vertical distance that is so small it is usually lost in roundoff error
                        ! and can be neglected [Z ~> m].
  real    :: botfn      ! A function that is 1 at the bottom and small far from it [nondim]
  real    :: z2         ! The distance from the bottom, normalized by Hbbl_gl90 [nondim]

  if (work_on_u) then
    Is = G%iscB ; Ie = G%iecB
    js = G%jsc ; je = G%jec
  else
    is = G%isc ; ie = G%iec
    Js = G%jscB ; Je = G%jecB
  endif

  nz = GV%ke

  h_neglect = GV%dZ_subroundoff
  kdgl90_use_vert_struct = .false.
  if (VarMix%use_variable_mixing) then
    kdgl90_use_vert_struct = allocated(VarMix%kdgl90_struct)
  endif

  a_cpl_gl90(:,:,:) = 0.

  do K=2,nz
    if (work_on_u) then
      ! compute coupling coefficient at u-points
      do j=js,je ; do I=Is,Ie; if (do_i(I,j)) then
        f2 = 0.25 * (G%CoriolisBu(I,J-1) + G%CoriolisBu(I,J))**2
        if (CS%use_GL90_N2) then
          a_cpl_gl90(I,j,K) = 2. * f2 * CS%alpha_gl90 / (hvel(I,j,k) + hvel(I,j,k-1) + h_neglect)
        else
          if (CS%read_kappa_gl90) then
            a_cpl_gl90(I,j,K) = f2 * 0.5 * (CS%kappa_gl90_2d(i,j) + CS%kappa_gl90_2d(i+1,j)) / GV%g_prime(K)
          else
            a_cpl_gl90(I,j,K) = f2 * CS%kappa_gl90 / GV%g_prime(K)
          endif
          if (kdgl90_use_vert_struct) then
            a_cpl_gl90(I,j,K) = a_cpl_gl90(I,j,K) * 0.5 &
                * (VarMix%kdgl90_struct(i,j,k-1) + VarMix%kdgl90_struct(i+1,j,k-1))
          endif
        endif
        ! botfn determines when a point is within the influence of the GL90 bottom boundary layer,
        ! going from 1 at the bottom to 0 in the interior.
        z2 = z_i(I,j,k)
        botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)
        a_cpl_gl90(I,j,K) = a_cpl_gl90(I,j,K) * (1. - botfn)
      endif; enddo ; enddo
    else
      ! compute viscosities at v-points
      do J=Js,Je ; do i=is,ie ; if (do_i(i,J)) then
        f2 = 0.25 * (G%CoriolisBu(I-1,J) + G%CoriolisBu(I,J))**2
        if (CS%use_GL90_N2) then
          a_cpl_gl90(i,J,K) = 2. * f2 * CS%alpha_gl90 / (hvel(i,J,k) + hvel(i,J,k-1) + h_neglect)
        else
          if (CS%read_kappa_gl90) then
            a_cpl_gl90(i,J,K) = f2 * 0.5 * (CS%kappa_gl90_2d(i,j) + CS%kappa_gl90_2d(i,j+1)) / GV%g_prime(K)
          else
            a_cpl_gl90(i,J,K) = f2 * CS%kappa_gl90 / GV%g_prime(K)
          endif
          if (kdgl90_use_vert_struct) then
            a_cpl_gl90(i,J,K) = a_cpl_gl90(i,J,K) * 0.5 &
                * (VarMix%kdgl90_struct(i,j,k-1) + VarMix%kdgl90_struct(i,j+1,k-1))
          endif
        endif
        ! botfn determines when a point is within the influence of the GL90 bottom boundary layer,
        ! going from 1 at the bottom to 0 in the interior.
        z2 = z_i(i,J,k)
        botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)
        a_cpl_gl90(i,J,K) = a_cpl_gl90(i,J,K) * (1. - botfn)
      endif ; enddo ; enddo
    endif
  enddo
end subroutine find_coupling_coef_gl90


!> Perform a fully implicit vertical diffusion
!! of momentum.  Stress top and bottom boundary conditions are used.
!!
!! This is solving the tridiagonal system
!! \f[ \left(h_k + a_{k + 1/2} + a_{k - 1/2} + r_k\right) u_k^{n+1}
!!     = h_k u_k^n + a_{k + 1/2} u_{k+1}^{n+1} + a_{k - 1/2} u_{k-1}^{n+1} \f]
!! where \f$a_{k + 1/2} = \Delta t \nu_{k + 1/2} / h_{k + 1/2}\f$
!! is the <em>interfacial coupling thickness per time step</em>,
!! encompassing background viscosity as well as contributions from
!! enhanced mixed and bottom layer viscosities.
!! $r_k$ is a Rayleigh drag term due to channel drag.
!! There is an additional stress term on the right-hand side
!! if DIRECT_STRESS is true, applied to the surface layer.
subroutine vertvisc(u, v, h, forces, visc, dt, OBC, ADp, CDp, G, GV, US, CS, &
                    taux_bot, tauy_bot, fpmix, Waves)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  type(mech_forcing),    intent(in)      :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),   intent(inout)   :: visc   !< Viscosities and bottom drag
  real,                  intent(in)      :: dt     !< Time increment [T ~> s]
  type(ocean_OBC_type),  pointer         :: OBC    !< Open boundary condition structure
  type(accel_diag_ptrs), intent(inout)   :: ADp    !< Accelerations in the momentum
                                                   !! equations for diagnostics
  type(cont_diag_ptrs),  intent(inout)   :: CDp    !< Continuity equation terms
  type(vertvisc_CS),     pointer         :: CS     !< Vertical viscosity control structure
  real, dimension(SZIB_(G),SZJ_(G)), &
                   optional, intent(out) :: taux_bot !< Zonal bottom stress from ocean to
                                                     !! rock [R L Z T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), &
                   optional, intent(out) :: tauy_bot !< Meridional bottom stress from ocean to
                                                     !! rock [R L Z T-2 ~> Pa]
  logical,         optional, intent(in)  :: fpmix !< fpmix along Eulerian shear
  type(wave_parameters_CS), &
                   optional, pointer     :: Waves !< Container for wave/Stokes information

  ! Fields from forces used in this subroutine:
  !   taux: Zonal wind stress [R L Z T-2 ~> Pa].
  !   tauy: Meridional wind stress [R L Z T-2 ~> Pa].

  ! Local variables

  real :: b1(SZIB_(G), SZJB_(G))
    ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: c1(SZIB_(G), SZJB_(G), SZK_(GV))
    ! A variable used by the tridiagonal solver [nondim].
  real :: d1(SZIB_(G), SZJB_(G))
    ! d1=1-c1 is used by the tridiagonal solver [nondim].
  real :: Ray(SZIB_(G), SZJB_(G))
    ! Ray is the Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
  real :: b_denom_1              ! The first term in the denominator of b1 [H ~> m or kg m-2].

  real :: Hmix             ! The mixed layer thickness over which stress
                           ! is applied with direct_stress [H ~> m or kg m-2].
  real :: I_Hmix           ! The inverse of Hmix [H-1 ~> m-1 or m2 kg-1].
  real :: Idt              ! The inverse of the time step [T-1 ~> s-1].
  real :: dt_Rho0          ! The time step divided by the mean density [T H Z-1 R-1 ~> s m3 kg-1 or s].
  real :: h_neglect        ! A thickness that is so small it is usually lost
                           ! in roundoff and can be neglected [H ~> m or kg m-2].

  real :: stress           !   The surface stress times the time step, divided
                           ! by the density [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real :: accel_underflow  ! An acceleration magnitude that is so small that values that are less
                           ! than this are diagnosed as 0 [L T-2 ~> m s-2].
  real :: zDS, h_a         ! Temporary thickness variables used with direct_stress [H ~> m or kg m-2]
  real :: hfr              ! Temporary ratio of thicknesses used with direct_stress [nondim]
  real :: surface_stress(SZIB_(G), SZJB_(G))
    ! The same as stress, unless the wind stress is applied as a body force
    ! [H L T-1 ~> m2 s-1 or kg m-1 s-1].
  real, allocatable, dimension(:,:,:) :: KE_term ! A term in the kinetic energy budget
                                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]
  real, allocatable, dimension(:,:,:) :: KE_u ! The area integral of a KE term in a layer at u-points
                                              ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real, allocatable, dimension(:,:,:) :: KE_v ! The area integral of a KE term in a layer at v-points
                                              ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]

  logical :: DoStokesMixing
  logical :: lfpmix

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, n
  is = G%isc ; ie = G%iec; js = G%jsc; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
         "Module must be initialized before it is used.")

  if (.not.CS%initialized) call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
         "Module must be initialized before it is used.")

  if (CS%id_GLwork > 0) then
    allocate(KE_u(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke), source=0.0)
    allocate(KE_v(G%isd:G%ied,G%JsdB:G%JedB,GV%ke), source=0.0)
    allocate(KE_term(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0)
    if (.not.G%symmetric) &
      call create_group_pass(CS%pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)
  endif

  if (CS%direct_stress) then
    Hmix = CS%Hmix_stress
    I_Hmix = 1.0 / Hmix
  endif
  dt_Rho0 = dt / GV%H_to_RZ
  h_neglect = GV%H_subroundoff
  Idt = 1.0 / dt

  accel_underflow = CS%vel_underflow * Idt

  !Check if Stokes mixing allowed if requested (present and associated)
  DoStokesMixing=.false.
  if (CS%StokesMixing) then
    if (present(Waves)) DoStokesMixing = associated(Waves)
    if (.not. DoStokesMixing) &
      call MOM_error(FATAL,"Stokes Mixing called without allocated"//&
                     "Waves Control Structure")
  endif
  lfpmix = .false.
  if ( present(fpmix) ) lfpmix = fpmix

  !   Update the zonal velocity component using a modification of a standard
  ! tridagonal solver.

  ! WGL: Brandon Reichl says the following is obsolete. u(I,j,k) already
  ! includes Stokes.
  ! When mixing down Eulerian current + Stokes drift add before calling solver
  if (DoStokesMixing) then
    do k=1,nz ; do j=G%jsc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
      u(I,j,k) = u(I,j,k) + Waves%Us_x(I,j,k)
    endif ; enddo ; enddo ; enddo
  endif

  if (lfpmix) then
    do k=1,nz ; do j=G%jsc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
      u(I,j,k) = u(I,j,k) - Waves%Us_x(I,j,k)
    endif ; enddo ; enddo ; enddo
  endif

  if (associated(ADp%du_dt_visc)) then
    do k=1,nz ; do j=G%jsc,G%jec ; do I=Isq,Ieq
      ADp%du_dt_visc(I,j,k) = u(I,j,k)
    enddo ; enddo; enddo
  endif

  if (associated(ADp%du_dt_visc_gl90)) then
    do k=1,nz ; do j=G%jsc,G%jec ; do I=Isq,Ieq
      ADp%du_dt_visc_gl90(I,j,k) = u(I,j,k)
    enddo ; enddo ; enddo
  endif

  if (associated(ADp%du_dt_str)) then
    do k=1,nz ; do j=G%jsc,G%jec ; do I=Isq,Ieq
      ADp%du_dt_str(I,j,k) = 0.0
    enddo ; enddo ; enddo
  endif

  !   One option is to have the wind stress applied as a body force
  ! over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
  ! the wind stress is applied as a stress boundary condition.
  if (CS%direct_stress) then
    do j=G%jsc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
      surface_stress(I,j) = 0.0
      zDS = 0.0
      stress = dt_Rho0 * forces%taux(I,j)
      do k=1,nz
        h_a = 0.5 * (h(i,j,k) + h(i+1,j,k)) + h_neglect
        hfr = 1.0 ; if ((zDS+h_a) > Hmix) hfr = (Hmix - zDS) / h_a
        u(I,j,k) = u(I,j,k) + I_Hmix * hfr * stress
        if (associated(ADp%du_dt_str)) ADp%du_dt_str(i,J,k) = (I_Hmix * hfr * stress) * Idt
        zDS = zDS + h_a ; if (zDS >= Hmix) exit
      enddo
    endif ; enddo ; enddo
  else
    do j=G%jsc,G%jec ; do I=Isq,Ieq
      surface_stress(I,j) = dt_Rho0 * (G%mask2dCu(I,j)*forces%taux(I,j))
    enddo ; enddo
  endif

  ! perform forward elimination on the tridiagonal system
  !
  ! denote the diagonal of the system as b_k, the subdiagonal as a_k
  ! and the superdiagonal as c_k. The right-hand side terms are d_k.
  !
  ! ignoring the Rayleigh drag contribution,
  ! we have a_k = -dt * a_u(k)
  !         b_k = h_u(k) + dt * (a_u(k) + a_u(k+1))
  !         c_k = -dt * a_u(k+1)
  !
  ! for forward elimination, we want to:
  ! calculate c'_k = - c_k                / (b_k + a_k c'_(k-1))
  ! and       d'_k = (d_k - a_k d'_(k-1)) / (b_k + a_k c'_(k-1))
  ! where c'_1 = c_1/b_1 and d'_1 = d_1/b_1
  !
  ! This form is mathematically equivalent to Thomas' tridiagonal matrix algorithm, but it
  ! does not suffer from the acute sensitivity to truncation errors of the Thomas algorithm
  ! because it involves no subtraction, as discussed by Schopf & Loughe, MWR, 1995.
  !
  ! b1 is the denominator term 1 / (b_k + a_k c'_(k-1))
  ! b_denom_1 is (b_k + a_k + c_k) - a_k(1 - c'_(k-1))
  !            = (b_k + c_k + c'_(k-1))
  ! this is done so that d1 = b1 * b_denom_1 = 1 - c'_(k-1)
  ! c1(k) is -c'_(k - 1)
  ! and the right-hand-side is destructively updated to be d'_k

  if (allocated(visc%Ray_u)) then
    do j=G%jsc,G%jec ; do I=Isq,Ieq
      Ray(I,j) = visc%Ray_u(I,j,1)
    enddo ; enddo
  else
    do j=G%jsc,G%jec ; do I=Isq,Ieq
      Ray(I,j) = 0.
    enddo ; enddo
  endif

  do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
    b_denom_1 = CS%h_u(I,j,1) + dt * (Ray(I,j) + CS%a_u(I,j,1))
    b1(I,j) = 1.0 / (b_denom_1 + dt*CS%a_u(I,j,2))
    d1(I,j) = b_denom_1 * b1(I,j)
    u(I,j,1) = b1(I,j) * (CS%h_u(I,j,1) * u(I,j,1) + surface_stress(I,j))
  endif ; enddo ; enddo

  if (associated(ADp%du_dt_str)) then
    do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
      ADp%du_dt_str(I,j,1) = b1(I,j) * (CS%h_u(I,j,1) * ADp%du_dt_str(I,j,1) + surface_stress(I,j) * Idt)
    endif ; enddo ; enddo
  endif

  do k=2,nz
    if (allocated(visc%Ray_u)) then
      do j=G%jsc,G%jec ; do I=Isq,Ieq
        Ray(I,j) = visc%Ray_u(I,j,k)
      enddo ; enddo
    endif

    do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
      c1(I,j,k) = dt * CS%a_u(I,j,K) * b1(I,j)
      b_denom_1 = CS%h_u(I,j,k) + dt * (Ray(I,j) + CS%a_u(I,j,K)*d1(I,j))
      b1(I,j) = 1.0 / (b_denom_1 + dt * CS%a_u(I,j,K+1))
      d1(I,j) = b_denom_1 * b1(I,j)
      u(I,j,k) = (CS%h_u(I,j,k) * u(I,j,k) + &
                  dt * CS%a_u(I,j,K) * u(I,j,k-1)) * b1(I,j)
    endif ; enddo ; enddo

    if (associated(ADp%du_dt_str)) then
      do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
        ADp%du_dt_str(I,j,k) = (CS%h_u(I,j,k) * ADp%du_dt_str(I,j,k) &
            + dt * CS%a_u(I,j,K) * ADp%du_dt_str(I,j,k-1)) * b1(I,j)
      endif ; enddo ; enddo
    endif
  enddo

  ! back substitute to solve for the new velocities
  ! u_k = d'_k - c'_k x_(k+1)
  do k=nz-1,1,-1
    do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
      u(I,j,k) = u(I,j,k) + c1(I,j,k+1) * u(I,j,k+1)
    endif ; enddo ; enddo
  enddo

  if (associated(ADp%du_dt_str)) then
    do j=G%isc,G%jec ; do I=Isq,Ieq
      if (abs(ADp%du_dt_str(I,j,nz)) < accel_underflow) &
        ADp%du_dt_str(I,j,nz) = 0.0
    enddo ; enddo

    do k=nz-1,1,-1
      do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
        ADp%du_dt_str(I,j,k) = ADp%du_dt_str(I,j,k) + c1(I,j,k+1) * ADp%du_dt_str(I,j,k+1)

        if (abs(ADp%du_dt_str(I,j,k)) < accel_underflow) &
          ADp%du_dt_str(I,j,k) = 0.0
      endif ; enddo ; enddo
    enddo
  endif

  ! compute vertical velocity tendency that arises from GL90 viscosity;
  ! follow tridiagonal solve method as above; to avoid corrupting u,
  ! use ADp%du_dt_visc_gl90 as a placeholder for updated u (due to GL90) until last do loop
  if ((CS%id_du_dt_visc_gl90 > 0) .or. (CS%id_GLwork > 0)) then
    if (associated(ADp%du_dt_visc_gl90)) then
      do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
        b_denom_1 = CS%h_u(I,j,1)  ! CS%a_u_gl90(I,j,1) is zero
        b1(I,j) = 1.0 / (b_denom_1 + dt*CS%a_u_gl90(I,j,2))
        d1(I,j) = b_denom_1 * b1(I,j)

        ADp%du_dt_visc_gl90(I,j,1) = b1(I,j) * (CS%h_u(I,j,1) * ADp%du_dt_visc_gl90(I,j,1))
      endif ; enddo ; enddo

      do k=2,nz
        do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
          c1(I,j,k) = dt * CS%a_u_gl90(I,j,K) * b1(I,j)
          b_denom_1 = CS%h_u(I,j,k) + dt * (CS%a_u_gl90(I,j,K)*d1(I,j))
          b1(I,j) = 1.0 / (b_denom_1 + dt * CS%a_u_gl90(I,j,K+1))
          d1(I,j) = b_denom_1 * b1(I,j)

          ADp%du_dt_visc_gl90(I,j,k) = (CS%h_u(I,j,k) * ADp%du_dt_visc_gl90(I,j,k) &
              + dt * CS%a_u_gl90(I,j,K) * ADp%du_dt_visc_gl90(I,j,k-1)) * b1(I,j)
        endif ; enddo ; enddo
      enddo

      ! back substitute to solve for new velocities, held by ADp%du_dt_visc_gl90
      do k=nz-1,1,-1
        do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
          ADp%du_dt_visc_gl90(I,j,k) = ADp%du_dt_visc_gl90(I,j,k) &
              + c1(I,j,k+1) * ADp%du_dt_visc_gl90(I,j,k+1)
        endif ; enddo ; enddo
      enddo

      do k=1,nz
        do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
          ! now fill ADp%du_dt_visc_gl90(I,j,k) with actual velocity tendency due to GL90;
          ! note that on RHS: ADp%du_dt_visc(I,j,k) holds the original velocity value u(I,j,k)
          ! and ADp%du_dt_visc_gl90(I,j,k) the updated velocity due to GL90
          ADp%du_dt_visc_gl90(I,j,k) = &
              (ADp%du_dt_visc_gl90(I,j,k) - ADp%du_dt_visc(I,j,k)) * Idt

          if (abs(ADp%du_dt_visc_gl90(I,j,k)) < accel_underflow) then
            ADp%du_dt_visc_gl90(I,j,k) = 0.0
          endif
        endif ; enddo ; enddo
      enddo

      ! to compute energetics, we need to multiply by u*h, where u is original velocity before
      ! velocity update; note that ADp%du_dt_visc(I,j,k) holds the original velocity value u(I,j,k)
      if (CS%id_GLwork > 0) then
        do k=1,nz
          do j=G%isc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
            KE_u(I,j,k) = ADp%du_dt_visc(I,j,k) * CS%h_u(I,j,k) * G%areaCu(I,j) * ADp%du_dt_visc_gl90(I,j,k)
          endif ; enddo ; enddo
        enddo
      endif
    endif
  endif

  if (associated(ADp%du_dt_visc)) then
    do k=1,nz ; do j=G%jsc,G%jec ; do I=Isq,Ieq
      ADp%du_dt_visc(I,j,k) = (u(I,j,k) - ADp%du_dt_visc(I,j,k)) * Idt

      if (abs(ADp%du_dt_visc(I,j,k)) < accel_underflow) &
        ADp%du_dt_visc(I,j,k) = 0.0
    enddo ; enddo ; enddo
  endif

  if (allocated(visc%taux_shelf)) then
    do j=G%jsc,G%jec ; do I=Isq,Ieq
      visc%taux_shelf(I,j) = -GV%H_to_RZ * CS%a1_shelf_u(I,j) * u(I,j,1) ! - u_shelf?
    enddo ; enddo
  endif

  if (present(taux_bot)) then
    do j=G%jsc,G%jec ; do I=Isq,Ieq
      taux_bot(I,j) = GV%H_to_RZ * (u(I,j,nz) * CS%a_u(I,j,nz+1))
    enddo ; enddo

    if (allocated(visc%Ray_u)) then
      do k=1,nz ; do j=G%jsc,G%jec ; do I=Isq,Ieq
        taux_bot(I,j) = taux_bot(I,j) + GV%H_to_RZ * (visc%Ray_u(I,j,k) * u(I,j,k))
      enddo ; enddo ; enddo
    endif
  endif

  ! When mixing down Eulerian current + Stokes drift subtract after calling solver
  if (DoStokesMixing) then
    do k=1,nz ; do j=G%jsc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
      u(I,j,k) = u(I,j,k) - Waves%Us_x(I,j,k)
    endif ; enddo ; enddo ; enddo
  endif

  if (lfpmix) then
    do k=1,nz ; do j=G%jsc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
      u(I,j,k) = u(I,j,k) + Waves%Us_x(I,j,k)
    endif ; enddo ; enddo ; enddo
  endif

  ! == Now work on the meridional velocity component.

  ! When mixing down Eulerian current + Stokes drift add before calling solver
  if (DoStokesMixing) then
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
      v(i,j,k) = v(i,j,k) + Waves%Us_y(i,j,k)
    endif ; enddo ; enddo ; enddo
  endif

  if (lfpmix) then
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
      v(i,j,k) = v(i,j,k) - Waves%Us_y(i,j,k)
    endif ; enddo ; enddo ; enddo
  endif

  if (associated(ADp%dv_dt_visc)) then
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      ADp%dv_dt_visc(i,J,k) = v(i,J,k)
    enddo ; enddo ; enddo
  endif

  if (associated(ADp%dv_dt_visc_gl90)) then
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      ADp%dv_dt_visc_gl90(i,J,k) = v(i,J,k)
    enddo ; enddo ; enddo
  endif

  if (associated(ADp%dv_dt_str)) then
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      ADp%dv_dt_str(i,J,k) = 0.0
    enddo ; enddo ; enddo
  endif

  !   One option is to have the wind stress applied as a body force
  ! over the topmost Hmix fluid.  If DIRECT_STRESS is not defined,
  ! the wind stress is applied as a stress boundary condition.
  if (CS%direct_stress) then
    do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
      surface_stress(i,J) = 0.0
      zDS = 0.0
      stress = dt_Rho0 * forces%tauy(i,J)
      do k=1,nz
        h_a = 0.5 * (h(i,J,k) + h(i,J+1,k)) + h_neglect
        hfr = 1.0 ; if ((zDS+h_a) > Hmix) hfr = (Hmix - zDS) / h_a
        v(i,J,k) = v(i,J,k) + I_Hmix * hfr * stress
        if (associated(ADp%dv_dt_str)) ADp%dv_dt_str(i,J,k) = (I_Hmix * hfr * stress) * Idt
        zDS = zDS + h_a ; if (zDS >= Hmix) exit
      enddo
    endif ; enddo ; enddo
  else
    do J=Jsq,Jeq ; do i=is,ie
      surface_stress(i,J) = dt_Rho0 * (G%mask2dCv(i,J) * forces%tauy(i,J))
    enddo ; enddo
  endif

  if (allocated(visc%Ray_v)) then
    do J=Jsq,Jeq ; do i=is,ie
      Ray(i,J) = visc%Ray_v(i,J,1)
    enddo ; enddo
  else
    do J=Jsq,Jeq ; do i=is,ie
      Ray(i,J) = 0.
    enddo ; enddo
  endif

  do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
    b_denom_1 = CS%h_v(i,J,1) + dt * (Ray(i,J) + CS%a_v(i,J,1))
    b1(i,J) = 1.0 / (b_denom_1 + dt*CS%a_v(i,J,2))
    d1(i,J) = b_denom_1 * b1(i,J)
    v(i,J,1) = b1(i,J) * (CS%h_v(i,J,1) * v(i,J,1) + surface_stress(i,J))
  endif ; enddo ; enddo

  if (associated(ADp%dv_dt_str)) then
    do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
      ADp%dv_dt_str(i,J,1) = b1(i,J) * (CS%h_v(i,J,1) * ADp%dv_dt_str(i,J,1) + surface_stress(i,J) * Idt)
    endif ; enddo ; enddo
  endif

  do k=2,nz
    if (allocated(visc%Ray_v)) then
      do J=Jsq,Jeq ; do i=is,ie
        Ray(i,J) = visc%Ray_v(i,J,k)
      enddo ; enddo
    endif

    do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
      c1(i,J,k) = dt * CS%a_v(i,J,K) * b1(i,J)
      b_denom_1 = CS%h_v(i,J,k) + dt * (Ray(i,J) + CS%a_v(i,J,K)*d1(i,J))
      b1(i,J) = 1.0 / (b_denom_1 + dt * CS%a_v(i,J,K+1))
      d1(i,J) = b_denom_1 * b1(i,J)
      v(i,J,k) = (CS%h_v(i,J,k) * v(i,J,k) + dt * CS%a_v(i,J,K) * v(i,J,k-1)) * b1(i,J)
    endif ; enddo ; enddo

    if (associated(ADp%dv_dt_str)) then
      do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
        ADp%dv_dt_str(i,J,k) = (CS%h_v(i,J,k) * ADp%dv_dt_str(i,J,k) &
            + dt * CS%a_v(i,J,K) * ADp%dv_dt_str(i,J,k-1)) * b1(i,J)
      endif ; enddo ; enddo
    endif
  enddo

  do k=nz-1,1,-1
    do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
      v(i,J,k) = v(i,J,k) + c1(i,J,k+1) * v(i,J,k+1)
    endif ; enddo ; enddo
  enddo

  if (associated(ADp%dv_dt_str)) then
    do J=Jsq,Jeq ; do i=is,ie
      if (abs(ADp%dv_dt_str(i,J,nz)) < accel_underflow) ADp%dv_dt_str(i,J,nz) = 0.0
    enddo ; enddo

    do k=nz-1,1,-1
      do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
        ADp%dv_dt_str(i,J,k) = ADp%dv_dt_str(i,J,k) + c1(i,J,k+1) * ADp%dv_dt_str(i,J,k+1)
        if (abs(ADp%dv_dt_str(i,J,k)) < accel_underflow) ADp%dv_dt_str(i,J,k) = 0.0
      endif ; enddo ; enddo
    enddo
  endif

  ! compute vertical velocity tendency that arises from GL90 viscosity;
  ! follow tridiagonal solve method as above; to avoid corrupting v,
  ! use ADp%dv_dt_visc_gl90 as a placeholder for updated u (due to GL90) until last do loop
  if ((CS%id_dv_dt_visc_gl90 > 0) .or. (CS%id_GLwork > 0)) then
    if (associated(ADp%dv_dt_visc_gl90)) then
      do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
        b_denom_1 = CS%h_v(i,J,1)  ! CS%a_v_gl90(i,J,1) is zero
        b1(i,J) = 1.0 / (b_denom_1 + dt*CS%a_v_gl90(i,J,2))
        d1(i,J) = b_denom_1 * b1(i,J)
        ADp%dv_dt_visc_gl90(I,J,1) = b1(i,J) * (CS%h_v(i,J,1) * ADp%dv_dt_visc_gl90(i,J,1))
      endif ; enddo ; enddo

      do k=2,nz
        do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
          c1(i,J,k) = dt * CS%a_v_gl90(i,J,K) * b1(i,J)
          b_denom_1 = CS%h_v(i,J,k) + dt * (CS%a_v_gl90(i,J,K)*d1(i,J))
          b1(i,J) = 1.0 / (b_denom_1 + dt * CS%a_v_gl90(i,J,K+1))
          d1(i,J) = b_denom_1 * b1(i,J)
          ADp%dv_dt_visc_gl90(i,J,k) = (CS%h_v(i,J,k) * ADp%dv_dt_visc_gl90(i,J,k) + &
                      dt * CS%a_v_gl90(i,J,K) * ADp%dv_dt_visc_gl90(i,J,k-1)) * b1(i,J)
        endif ; enddo ; enddo
      enddo

      ! back substitute to solve for new velocities, held by ADp%dv_dt_visc_gl90
      do k=nz-1,1,-1
        do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
          ADp%dv_dt_visc_gl90(i,J,k) = ADp%dv_dt_visc_gl90(i,J,k) + c1(i,J,k+1) * ADp%dv_dt_visc_gl90(i,J,k+1)
        endif ; enddo ; enddo
      enddo

      do k=1,nz
        do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
          ! now fill ADp%dv_dt_visc_gl90(i,J,k) with actual velocity tendency due to GL90;
          ! note that on RHS: ADp%dv_dt_visc(i,J,k) holds the original velocity value v(i,J,k)
          ! and ADp%dv_dt_visc_gl90(i,J,k) the updated velocity due to GL90
          ADp%dv_dt_visc_gl90(i,J,k) = (ADp%dv_dt_visc_gl90(i,J,k) - ADp%dv_dt_visc(i,J,k))*Idt
          if (abs(ADp%dv_dt_visc_gl90(i,J,k)) < accel_underflow) ADp%dv_dt_visc_gl90(i,J,k) = 0.0
        endif ; enddo ; enddo
      enddo

      ! to compute energetics, we need to multiply by v*h, where u is original velocity before
      ! velocity update; note that ADp%dv_dt_visc(I,j,k) holds the original velocity value v(i,J,k)
      if (CS%id_GLwork > 0) then
        do k=1,nz
          do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
            ! note that on RHS: ADp%dv_dt_visc(I,j,k) holds the original velocity value v(I,j,k)
            KE_v(I,j,k) = ADp%dv_dt_visc(i,J,k) * CS%h_v(i,J,k) * G%areaCv(i,J) * ADp%dv_dt_visc_gl90(i,J,k)
          endif ; enddo ; enddo
        enddo
      endif
    endif
  endif

  if (associated(ADp%dv_dt_visc)) then
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      ADp%dv_dt_visc(i,J,k) = (v(i,J,k) - ADp%dv_dt_visc(i,J,k))*Idt
      if (abs(ADp%dv_dt_visc(i,J,k)) < accel_underflow) ADp%dv_dt_visc(i,J,k) = 0.0
    enddo ; enddo ; enddo
  endif

  if (allocated(visc%tauy_shelf)) then
    do J=Jsq,Jeq ; do i=is,ie
      visc%tauy_shelf(i,J) = -GV%H_to_RZ * CS%a1_shelf_v(i,J) * v(i,J,1) ! - v_shelf?
    enddo ; enddo
  endif

  if (present(tauy_bot)) then
    do J=Jsq,Jeq ; do i=is,ie
      tauy_bot(i,J) = GV%H_to_RZ * (v(i,J,nz) * CS%a_v(i,J,nz+1))
    enddo; enddo

    if (allocated(visc%Ray_v)) then
      do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
        tauy_bot(i,J) = tauy_bot(i,J) + GV%H_to_RZ * (visc%Ray_v(i,J,k)*v(i,J,k))
      enddo ; enddo ; enddo
    endif
  endif

  ! When mixing down Eulerian current + Stokes drift subtract after calling solver
  if (DoStokesMixing) then
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
      v(i,J,k) = v(i,J,k) - Waves%Us_y(i,J,k)
    endif ; enddo ; enddo ; enddo
  endif

  if (lfpmix) then
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
      v(i,J,k) = v(i,J,k) + Waves%Us_y(i,J,k)
    endif ; enddo ; enddo ; enddo
  endif

  ! Calculate the KE source from GL90 vertical viscosity [H L2 T-3 ~> m3 s-3].
  ! We do the KE-rate calculation here (rather than in MOM_diagnostics) to ensure
  ! a sign-definite term. MOM_diagnostics does not have access to the velocities
  ! and thicknesses used in the vertical solver, but rather uses a time-mean
  ! barotropic transport [uv]h.
  if (CS%id_GLwork > 0) then
    if (.not.G%symmetric) &
      call do_group_pass(CS%pass_KE_uv, G%domain)
    do k=1,nz
      do j=js,je ; do i=is,ie
        KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
            * (KE_u(I,j,k) + KE_u(I-1,j,k) + KE_v(i,J,k) + KE_v(i,J-1,k))
      enddo ; enddo
    enddo
    call post_data(CS%id_GLwork, KE_term, CS%diag)
  endif

  call vertvisc_limit_vel(u, v, h, ADp, CDp, forces, visc, dt, G, GV, US, CS)

  ! Here the velocities associated with open boundary conditions are applied.
  if (associated(OBC)) then
    do n=1,OBC%number_of_segments
      if (OBC%segment(n)%specified) then
        if (OBC%segment(n)%is_N_or_S) then
          J = OBC%segment(n)%HI%JsdB
          do k=1,nz ; do i=OBC%segment(n)%HI%isd,OBC%segment(n)%HI%ied
            v(i,J,k) = OBC%segment(n)%normal_vel(i,J,k)
          enddo ; enddo
        elseif (OBC%segment(n)%is_E_or_W) then
          I = OBC%segment(n)%HI%IsdB
          do k=1,nz ; do j=OBC%segment(n)%HI%jsd,OBC%segment(n)%HI%jed
            u(I,j,k) = OBC%segment(n)%normal_vel(I,j,k)
          enddo ; enddo
        endif
      endif
    enddo
  endif

  ! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_du_dt_visc > 0) &
      call post_data(CS%id_du_dt_visc, ADp%du_dt_visc, CS%diag)
    if (CS%id_du_dt_visc_gl90 > 0) &
      call post_data(CS%id_du_dt_visc_gl90, ADp%du_dt_visc_gl90, CS%diag)
    if (CS%id_dv_dt_visc > 0) &
      call post_data(CS%id_dv_dt_visc, ADp%dv_dt_visc, CS%diag)
    if (CS%id_dv_dt_visc_gl90 > 0) &
      call post_data(CS%id_dv_dt_visc_gl90, ADp%dv_dt_visc_gl90, CS%diag)
    if (present(taux_bot) .and. (CS%id_taux_bot > 0)) &
      call post_data(CS%id_taux_bot, taux_bot, CS%diag)
    if (present(tauy_bot) .and. (CS%id_tauy_bot > 0)) &
      call post_data(CS%id_tauy_bot, tauy_bot, CS%diag)
    if (CS%id_du_dt_str > 0) &
      call post_data(CS%id_du_dt_str, ADp%du_dt_str, CS%diag)
    if (CS%id_dv_dt_str > 0) &
      call post_data(CS%id_dv_dt_str, ADp%dv_dt_str, CS%diag)

    if (associated(ADp%du_dt_visc) .and. associated(ADp%du_dt_visc)) then
      ! Diagnostics of the fractional thicknesses times momentum budget terms
      ! 3D diagnostics of hf_du(dv)_dt_visc are commented because there is no clarity on proper remapping grid option.
      ! The code is retained for debugging purposes in the future.
      !if (CS%id_hf_du_dt_visc > 0) &
      !  call post_product_u(CS%id_hf_du_dt_visc, ADp%du_dt_visc, ADp%diag_hfrac_u, G, nz, CS%diag)
      !if (CS%id_hf_dv_dt_visc > 0) &
      !  call post_product_v(CS%id_hf_dv_dt_visc, ADp%dv_dt_visc, ADp%diag_hfrac_v, G, nz, CS%diag)

      ! Diagnostics for thickness-weighted vertically averaged viscous accelerations
      if (CS%id_hf_du_dt_visc_2d > 0) &
        call post_product_sum_u(CS%id_hf_du_dt_visc_2d, ADp%du_dt_visc, ADp%diag_hfrac_u, G, nz, CS%diag)
      if (CS%id_hf_dv_dt_visc_2d > 0) &
        call post_product_sum_v(CS%id_hf_dv_dt_visc_2d, ADp%dv_dt_visc, ADp%diag_hfrac_v, G, nz, CS%diag)

      ! Diagnostics for thickness x viscous accelerations
      if (CS%id_h_du_dt_visc > 0) call post_product_u(CS%id_h_du_dt_visc, ADp%du_dt_visc, ADp%diag_hu, G, nz, CS%diag)
      if (CS%id_h_dv_dt_visc > 0) call post_product_v(CS%id_h_dv_dt_visc, ADp%dv_dt_visc, ADp%diag_hv, G, nz, CS%diag)
    endif

    if (associated(ADp%du_dt_str) .and.  associated(ADp%dv_dt_str)) then
      ! Diagnostics for thickness x wind stress accelerations
      if (CS%id_h_du_dt_str > 0) call post_product_u(CS%id_h_du_dt_str, ADp%du_dt_str, ADp%diag_hu, G, nz, CS%diag)
      if (CS%id_h_dv_dt_str > 0) call post_product_v(CS%id_h_dv_dt_str, ADp%dv_dt_str, ADp%diag_hv, G, nz, CS%diag)

      ! Diagnostics for wind stress accelerations multiplied by visc_rem_[uv],
      if (CS%id_du_dt_str_visc_rem > 0) &
        call post_product_u(CS%id_du_dt_str_visc_rem, ADp%du_dt_str, ADp%visc_rem_u, G, nz, CS%diag)
      if (CS%id_dv_dt_str_visc_rem > 0) &
        call post_product_v(CS%id_dv_dt_str_visc_rem, ADp%dv_dt_str, ADp%visc_rem_v, G, nz, CS%diag)
    endif
  endif

end subroutine vertvisc


!> Calculate the fraction of momentum originally in a layer that remains in the water column
!! after a time-step of viscosity, equivalently the fraction of a time-step's worth of
!! barotropic acceleration that a layer experiences after viscosity is applied.
subroutine vertvisc_remnant(visc, visc_rem_u, visc_rem_v, dt, G, GV, US, CS)
  type(ocean_grid_type), intent(in)   :: G    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV   !< Ocean vertical grid structure
  type(vertvisc_type),   intent(in)   :: visc !< Viscosities and bottom drag
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                         intent(inout) :: visc_rem_u !< Fraction of a time-step's worth of a
                                              !! barotropic acceleration that a layer experiences after
                                              !! viscosity is applied in the zonal direction [nondim]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                         intent(inout) :: visc_rem_v !< Fraction of a time-step's worth of a
                                              !! barotropic acceleration that a layer experiences after
                                              !! viscosity is applied in the meridional direction [nondim]
  real,                  intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  type(vertvisc_CS),     pointer       :: CS  !< Vertical viscosity control structure

  ! Local variables

  real :: b1(SZIB_(G),SZJB_(G))
    ! A variable used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].
  real :: c1(SZIB_(G),SZJB_(G),SZK_(GV))
    ! A variable used by the tridiagonal solver [nondim].
  real :: d1(SZIB_(G),SZJB_(G))
    ! d1=1-c1 is used by the tridiagonal solver [nondim].
  real :: Ray(SZIB_(G),SZJB_(G))
    ! Ray is the Rayleigh-drag velocity [H T-1 ~> m s-1 or Pa s m-1]
  real :: b_denom_1   ! The first term in the denominator of b1 [H ~> m or kg m-2].

  integer :: i, j, k, is, ie, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke

  if (.not.associated(CS)) call MOM_error(FATAL,"MOM_vert_friction(visc): "// &
         "Module must be initialized before it is used.")

  if (.not.CS%initialized) call MOM_error(FATAL,"MOM_vert_friction(remant): "// &
         "Module must be initialized before it is used.")

  ! Find the zonal viscous remnant using a modification of a standard tridagonal solver.
  if (allocated(visc%Ray_u)) then
    do j=G%jsc,G%jec ; do I=Isq,Ieq
      Ray(I,j) = visc%Ray_u(I,j,1)
    enddo ; enddo
  else
    do j=G%jsc,G%jec ; do I=Isq,Ieq
      Ray(I,j) = 0.
    enddo ; enddo
  endif

  do j=G%jsc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
    b_denom_1 = CS%h_u(I,j,1) + dt * (Ray(I,j) + CS%a_u(I,j,1))
    b1(I,j) = 1.0 / (b_denom_1 + dt * CS%a_u(I,j,2))
    d1(I,j) = b_denom_1 * b1(I,j)
    visc_rem_u(I,j,1) = b1(I,j) * CS%h_u(I,j,1)
  endif ; enddo ; enddo

  do k=2,nz
    if (allocated(visc%Ray_u)) then
      do j=G%jsc,G%jec ; do I=Isq,Ieq
        Ray(I,j) = visc%Ray_u(I,j,k)
      enddo ; enddo
    endif

    do j=G%jsc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
      c1(I,j,k) = dt * CS%a_u(I,j,K)*b1(I,j)
      b_denom_1 = CS%h_u(I,j,k) + dt * (Ray(I,j) + CS%a_u(I,j,K) * d1(I,j))
      b1(I,j) = 1.0 / (b_denom_1 + dt * CS%a_u(I,j,K+1))
      d1(I,j) = b_denom_1 * b1(I,j)
      visc_rem_u(I,j,k) = (CS%h_u(I,j,k) + dt * CS%a_u(I,j,K) * visc_rem_u(I,j,k-1)) * b1(I,j)
    endif ; enddo ; enddo
  enddo

  do k=nz-1,1,-1
    do j=G%jsc,G%jec ; do I=Isq,Ieq ; if (G%mask2dCu(I,j) > 0.) then
      visc_rem_u(I,j,k) = visc_rem_u(I,j,k) + c1(I,j,k+1) * visc_rem_u(I,j,k+1)
    endif ; enddo ; enddo
  enddo

  ! Now find the meridional viscous remnant using the robust tridiagonal solver.
  if (allocated(visc%Ray_v)) then
    do J=Jsq,Jeq ; do i=is,ie
      Ray(i,J) = visc%Ray_v(i,J,1)
    enddo ; enddo
  else
    do J=Jsq,Jeq ; do i=is,ie
      Ray(i,J) = 0.
    enddo ; enddo
  endif

  do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
    b_denom_1 = CS%h_v(i,J,1) + dt * (Ray(i,J) + CS%a_v(i,J,1))
    b1(i,J) = 1.0 / (b_denom_1 + dt*CS%a_v(i,J,2))
    d1(i,J) = b_denom_1 * b1(i,J)
    visc_rem_v(i,J,1) = b1(i,J) * CS%h_v(i,J,1)
  endif ; enddo ; enddo

  do k=2,nz
    if (allocated(visc%Ray_v)) then
      do J=Jsq,Jeq ; do i=is,ie
        Ray(i,J) = visc%Ray_v(i,J,k)
      enddo ; enddo
    endif

    do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
      c1(i,J,k) = dt * CS%a_v(i,J,K) * b1(i,J)
      b_denom_1 = CS%h_v(i,J,k) + dt * (Ray(i,J) + CS%a_v(i,J,K) * d1(i,J))
      b1(i,J) = 1.0 / (b_denom_1 + dt * CS%a_v(i,J,K+1))
      d1(i,J) = b_denom_1 * b1(i,J)
      visc_rem_v(i,J,k) = (CS%h_v(i,J,k) + dt * CS%a_v(i,J,K) * visc_rem_v(i,J,k-1)) * b1(i,J)
    endif ; enddo ; enddo
  enddo

  do k=nz-1,1,-1
    do J=Jsq,Jeq ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.) then
      visc_rem_v(i,J,k) = visc_rem_v(i,J,k) + c1(i,J,k+1) * visc_rem_v(i,J,k+1)
    endif ; enddo ; enddo ! i and k loops
  enddo

  if (CS%debug) then
    call uvchksum("visc_rem_[uv]", visc_rem_u, visc_rem_v, G%HI, haloshift=0, &
                  scalar_pair=.true.)
  endif
end subroutine vertvisc_remnant


!> Calculate the coupling coefficients (CS%a_u, CS%a_v, CS%a_u_gl90, CS%a_v_gl90)
!! and effective layer thicknesses (CS%h_u and CS%h_v) for later use in the
!! applying the implicit vertical viscosity via vertvisc().
subroutine vertvisc_coef(u, v, h, dz, forces, visc, tv, dt, G, GV, US, CS, OBC, VarMix)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dz     !< Vertical distance across layers [Z ~> m]
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(in)    :: visc   !< Viscosities and bottom drag
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any available
                                                   !! thermodynamic fields.
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(vertvisc_CS),       intent(inout) :: CS     !< Vertical viscosity control structure
  type(ocean_OBC_type),    pointer       :: OBC    !< Open boundary condition structure
  type(VarMix_CS),         intent(in) :: VarMix !< Variable mixing coefficients
  ! Field from forces used in this subroutine:
  !   ustar: the friction velocity [Z T-1 ~> m s-1], used here as the mixing
  !     velocity in the mixed layer if NKML > 1 in a bulk mixed layer.

  ! Local variables

  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)) :: &
    hvel, &     ! hvel is the thickness used at a velocity grid point [H ~> m or kg m-2].
    dz_harm, &  ! Harmonic mean of the vertical distances around a velocity grid point,
                ! given by 2*(h+ * h-)/(h+ + h-) [Z ~> m].
    dz_vel, &   ! The vertical distance between interfaces used at a velocity grid point [Z ~> m].
    hvel_shelf, & ! The equivalent of hvel under shelves [H ~> m or kg m-2].
    dz_vel_shelf ! The equivalent of dz_vel under shelves [Z ~> m].
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    h_harm, &   ! Harmonic mean of the thicknesses around a velocity grid point,
                ! given by 2*(h+ * h-)/(h+ + h-) [H ~> m or kg m-2].
    h_arith, &  ! The arithmetic mean thickness [H ~> m or kg m-2].
    h_delta, &  ! The lateral difference of thickness [H ~> m or kg m-2].
    dz_arith    ! The arithmetic mean of the vertical distances around a velocity grid point [Z ~> m]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1) :: &
    z_i, &      ! An estimate of each interface's height above the bottom,
                ! normalized by the bottom boundary layer thickness [nondim]
    z_i_gl90, & ! An estimate of each interface's height above the bottom,
                ! normalized by the GL90 bottom boundary layer thickness [nondim]
    a_cpl, &    ! The drag coefficients across interfaces [H T-1 ~> m s-1 or Pa s m-1].  a_cpl times
                ! the velocity difference gives the stress across an interface.
    a_cpl_gl90, & ! The drag coefficients across interfaces associated with GL90 [H T-1 ~> m s-1 or Pa s m-1].
                ! a_cpl_gl90 times the velocity difference gives the GL90 stress across an interface.
                ! a_cpl_gl90 is part of a_cpl.
    a_shelf     ! The drag coefficients across interfaces in water columns under
                ! ice shelves [H T-1 ~> m s-1 or Pa s m-1].
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    kv_bbl, &     ! The bottom boundary layer viscosity [H Z T-1 ~> m2 s-1 or Pa s].
    bbl_thick, &  ! The bottom boundary layer thickness [Z ~> m].
    I_Hbbl, &     ! The inverse of the bottom boundary layer thickness [Z-1 ~> m-1].
    I_Hbbl_gl90, &! The inverse of the bottom boundary layer thickness used for the GL90 scheme
                  ! [Z-1 ~> m-1].
    I_HTbl, &     ! The inverse of the top boundary layer thickness [Z-1 ~> m-1].
    Ztop_min, &   ! The deeper of the two adjacent surface heights [Z ~> m].
    Dmin, &       ! The shallower of the two adjacent bottom depths [Z ~> m].
    zh, &         ! An estimate of the interface's distance from the bottom
                  ! based on harmonic mean thicknesses [Z ~> m].
    h_ml          ! The mixed layer depth [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    Ustar_2d    ! The wind friction velocity, calculated using the Boussinesq reference density or
                ! the time-evolving surface density in non-Boussinesq mode [Z T-1 ~> m s-1]
  real, allocatable, dimension(:,:) :: hML_u ! Diagnostic of the mixed layer depth at u points [Z ~> m].
  real, allocatable, dimension(:,:) :: hML_v ! Diagnostic of the mixed layer depth at v points [Z ~> m].
  real, allocatable, dimension(:,:,:) :: Kv_u ! Total vertical viscosity at u-points in
                                              ! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, allocatable, dimension(:,:,:) :: Kv_v ! Total vertical viscosity at v-points in
                                              ! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, allocatable, dimension(:,:,:) :: Kv_gl90_u ! GL90 vertical viscosity at u-points in
                                              ! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real, allocatable, dimension(:,:,:) :: Kv_gl90_v ! GL90 vertical viscosity at v-points in
                                              ! thickness-based units [H2 T-1 ~> m2 s-1 or kg2 m-4 s-1].
  real :: zcol(SZI_(G), SZJ_(G)) ! The height of an interface at h-points [Z ~> m].
  real :: botfn   ! A function which goes from 1 at the bottom to 0 much more
                  ! than Hbbl into the interior [nondim].
  real :: topfn   ! A function which goes from 1 at the top to 0 much more
                  ! than Htbl into the interior [nondim].
  real :: z2      ! The distance from the bottom, normalized by Hbbl [nondim]
  real :: z2_wt   ! A nondimensional (0-1) weight used when calculating z2 [nondim].
  real :: z_clear ! The clearance of an interface above the surrounding topography [Z ~> m].
  real :: a_cpl_max  ! The maximum drag coefficient across interfaces, set so that it will be
                     ! representable as a 32-bit float in MKS units [H T-1 ~> m s-1 or Pa s m-1]
  real :: h_neglect  ! A thickness that is so small it is usually lost
                     ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: dz_neglect ! A vertical distance that is so small it is usually lost
                     ! in roundoff and can be neglected [Z ~> m].

  real :: I_valBL ! The inverse of a scaling factor determining when water is
                  ! still within the boundary layer, as determined by the sum
                  ! of the harmonic mean thicknesses [nondim].
  logical :: do_i(SZIB_(G), SZJB_(G))
    ! Land mask
  logical :: do_i_shelf(SZIB_(G), SZJB_(G))
    ! Land mask with fractional shelf
  logical :: do_any_shelf
  integer, dimension(SZIB_(G), SZJB_(G)) :: &
    zi_dir   !  A trinary logical array indicating which thicknesses to use for
             !  finding z_clear.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, ij
  integer :: is_N_OBC, is_S_OBC, Is_E_OBC, Is_W_OBC, ie_N_OBC, ie_S_OBC, Ie_E_OBC, Ie_W_OBC
  integer :: js_N_OBC, js_S_OBC, Js_E_OBC, Js_W_OBC, je_N_OBC, je_S_OBC, Je_E_OBC, Je_W_OBC

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke

  if (.not.CS%initialized) call MOM_error(FATAL,"MOM_vert_friction(coef): "// &
         "Module must be initialized before it is used.")

  h_neglect = GV%H_subroundoff
  dz_neglect = GV%dZ_subroundoff
  a_cpl_max = 1.0e37 * GV%m_to_H * US%T_to_s
  I_Hbbl(:,:) = 1. / (CS%Hbbl + dz_neglect)
  if (CS%use_GL90_in_SSW) then
    I_Hbbl_gl90(:,:) = 1. / (CS%Hbbl_gl90 + dz_neglect)
  endif
  I_valBL = 0.0 ; if (CS%harm_BL_val > 0.0) I_valBL = 1.0 / CS%harm_BL_val

  if (CS%id_Kv_u > 0) allocate(Kv_u(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke), source=0.0)

  if (CS%id_Kv_v > 0) allocate(Kv_v(G%isd:G%ied,G%JsdB:G%JedB,GV%ke), source=0.0)

  if (CS%id_Kv_gl90_u > 0) allocate(Kv_gl90_u(G%IsdB:G%IedB,G%jsd:G%jed,GV%ke), source=0.0)

  if (CS%id_Kv_gl90_v > 0) allocate(Kv_gl90_v(G%isd:G%ied,G%JsdB:G%JedB,GV%ke), source=0.0)

  if (CS%debug .or. (CS%id_hML_u > 0)) allocate(hML_u(G%IsdB:G%IedB,G%jsd:G%jed), source=0.0)
  if (CS%debug .or. (CS%id_hML_v > 0)) allocate(hML_v(G%isd:G%ied,G%JsdB:G%JedB), source=0.0)

  if ((allocated(visc%taux_shelf) .or. associated(forces%frac_shelf_u)) .and. &
      .not.associated(CS%a1_shelf_u)) then
    allocate(CS%a1_shelf_u(G%IsdB:G%IedB,G%jsd:G%jed), source=0.0)
  endif
  if ((allocated(visc%tauy_shelf) .or. associated(forces%frac_shelf_v)) .and. &
      .not.associated(CS%a1_shelf_v)) then
    allocate(CS%a1_shelf_v(G%isd:G%ied,G%JsdB:G%JedB), source=0.0)
  endif

  if (associated(OBC)) then
    ! Set the ranges that contain various orientations of OBCs on this PE.
    is_N_OBC = max(is, OBC%is_v_N_obc) ; ie_N_OBC = min(ie, OBC%ie_v_N_obc)
    is_S_OBC = max(is, OBC%is_v_S_obc) ; ie_S_OBC = min(ie, OBC%ie_v_S_obc)
    Js_N_OBC = max(Jsq, OBC%Js_v_N_obc) ; Je_N_OBC = min(Jeq, OBC%Je_v_N_obc)
    Js_S_OBC = max(Jsq, OBC%Js_v_S_obc) ; Je_S_OBC = min(Jeq, OBC%Je_v_S_obc)
    Is_E_OBC = max(Isq, OBC%Is_u_E_obc) ; Ie_E_OBC = min(Ieq, OBC%Ie_u_E_obc)
    Is_W_OBC = max(Isq, OBC%Is_u_W_obc) ; Ie_W_OBC = min(Ieq, OBC%Ie_u_W_obc)
    js_E_OBC = max(js, OBC%Js_u_E_obc) ; je_E_OBC = min(je, OBC%je_u_E_obc)
    js_W_OBC = max(js, OBC%Js_u_W_obc) ; je_W_OBC = min(je, OBC%je_u_W_obc)
  endif

  call find_ustar(forces, tv, Ustar_2d, G, GV, US, halo=1)

  ! First do u-points

  ! Force IPO optimizations (e.g. Intel)
  ij = touch_ij(i,j)

  do j=js,je ; do I=Isq,Ieq
    do_i(I,j) = G%mask2dCu(I,j) > 0.
  enddo ; enddo

  if (CS%bottomdraglaw) then
    do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
      kv_bbl(I,j) = visc%Kv_bbl_u(I,j)
      bbl_thick(I,j) = visc%bbl_thick_u(I,j) + dz_neglect
      I_Hbbl(I,j) = 1. / bbl_thick(I,j)
    endif ; enddo ; enddo
  endif

  do j=js,je ; do I=Isq,Ieq
    Dmin(I,j) = min(G%bathyT(i,j), G%bathyT(i+1,j))
    zi_dir(I,j) = 0
  enddo ; enddo

  ! Project thickness outward across OBCs using a zero-gradient condition.
  if (associated(OBC)) then
    if (OBC%u_E_OBCs_on_PE) then
      do j=js_E_OBC,je_E_OBC ; do I=Is_E_OBC,Ie_E_OBC
        if (do_i(I,j) .and. OBC%segnum_u(I,j) > 0) then
          Dmin(I,j) = G%bathyT(i,j)
          zi_dir(I,j) = -1
        endif
      enddo ; enddo
    endif

    if (OBC%u_W_OBCs_on_PE) then
      do j=js_W_OBC,je_W_OBC ; do I=Is_W_OBC,Is_W_OBC
        if (do_i(I,j) .and. OBC%segnum_u(I,j) < 0) then
          Dmin(I,j) = G%bathyT(i+1,j)
          zi_dir(I,j) = 1
        endif
      enddo ; enddo
    endif
  endif

  ! The following block calculates the thicknesses at velocity grid points for
  ! the vertical viscosity (hvel and dz_vel).  Near the bottom an upwind biased
  ! thickness is used to control the effect of spurious Montgomery potential
  ! gradients at the bottom where nearly massless layers layers ride over the
  ! topography.

  do j=js,je ; do I=Isq,Ieq
    z_i(I,j,nz+1) = 0.
  enddo ; enddo

  if (.not. CS%harmonic_visc) then
    do j=js,je ; do I=Isq,Ieq
      zh(I,j) = 0.
    enddo ; enddo

    do j=js,je ; do I=Isq,Ieq+1
      zcol(i,j) = -G%bathyT(i,j)
    enddo ; enddo
  endif

  if (CS%use_GL90_in_SSW) then
    do j=js,je ; do I=Isq,Ieq
      z_i_gl90(I,j,nz+1) = 0.
    enddo ; enddo
  endif

  do k=nz,1,-1
    do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
      h_harm(I,j) = 2. * h(i,j,k) * h(i+1,j,k) / (h(i,j,k) + h(i+1,j,k) + h_neglect)
      h_arith(I,j) = 0.5 * (h(i+1,j,k) + h(i,j,k))
      h_delta(I,j) = h(i+1,j,k) - h(i,j,k)
      dz_harm(I,j,k) = 2. * dz(i,j,k) * dz(i+1,j,k) / (dz(i,j,k) + dz(i+1,j,k) + dz_neglect)
      dz_arith(I,j) = 0.5 * (dz(i+1,j,k) + dz(i,j,k))
    endif ; enddo ; enddo

    ! Project thickness outward across OBCs using a zero-gradient condition.
    if (associated(OBC)) then
      if (OBC%u_E_OBCs_on_PE) then
        do j=js_E_OBC,je_E_OBC ; do I=Is_E_OBC,Ie_E_OBC
          if (do_i(I,j) .and. OBC%segnum_u(I,j) > 0) then
            h_harm(I,j) = h(i,j,k)
            h_arith(I,j) = h(i,j,k)
            h_delta(I,j) = 0.
            dz_harm(I,j,k) = dz(i,j,k)
            dz_arith(I,j) = dz(i,j,k)
          endif
        enddo ; enddo
      endif

      if (OBC%u_W_OBCs_on_PE) then
        do j=js_W_OBC,je_W_OBC ; do I=Is_W_OBC,Ie_W_OBC
          if (do_i(I,j) .and. OBC%segnum_u(I,j) < 0) then
            h_harm(I,j) = h(i+1,j,k)
            h_arith(I,j) = h(i+1,j,k)
            h_delta(I,j) = 0.
            dz_harm(I,j,k) = dz(i+1,j,k)
            dz_arith(I,j) = dz(i+1,j,k)
          endif
        enddo ; enddo
      endif
    endif

    if (CS%harmonic_visc) then
      ! The following block calculates the thicknesses at velocity grid points
      ! for the vertical viscosity (hvel and dz_vel).  Near the bottom an
      ! upwind biased thickness is used to control the effect of spurious
      ! Montgomery potential gradients at the bottom where nearly massless
      ! layers ride over the topography.

      do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
        hvel(I,j,k) = h_harm(I,j)
        dz_vel(I,j,k) = dz_harm(I,j,k)

        if (u(I,j,k) * h_delta(I,j) < 0) then
          z2 = z_i(I,j,k+1)
          botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

          hvel(I,j,k) = (1. - botfn) * h_harm(I,j) + botfn * h_arith(I,j)
          dz_vel(I,j,k) = (1. - botfn) * dz_harm(I,j,k) + botfn * dz_arith(I,j)
        endif

        z_i(I,j,k) =  z_i(I,j,k+1) + dz_harm(I,j,k) * I_Hbbl(I,j)
      endif ; enddo ; enddo
    else
      do j=js,je ; do I=Isq,Ieq+1
        zcol(i,j) = zcol(i,j) + dz(i,j,k)
      enddo ; enddo

      do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
        zh(I,j) = zh(I,j) + dz_harm(I,j,k)

        z_clear = max(zcol(i,j),zcol(i+1,j)) + Dmin(I,j)
        if (zi_dir(I,j) < 0) z_clear = zcol(i,j) + Dmin(I,j)
        if (zi_dir(I,j) > 0) z_clear = zcol(i+1,j) + Dmin(I,j)

        z_i(I,j,k) = max(zh(I,j), z_clear) * I_Hbbl(I,j)

        hvel(I,j,k) = h_arith(I,j)
        dz_vel(I,j,k) = dz_arith(I,j)

        if (u(I,j,k) * h_delta(I,j) > 0.) then
          if (zh(I,j) * I_Hbbl(I,j) < CS%harm_BL_val) then
            hvel(I,j,k) = h_harm(I,j)
            dz_vel(I,j,k) = dz_harm(I,j,k)
          else
            z2_wt = 1.
            if (zh(I,j) * I_Hbbl(I,j) < 2. * CS%harm_BL_val) &
              z2_wt = max(0., min(1., zh(I,j) * I_Hbbl(I,j) * I_valBL - 1.))

            z2 = z2_wt * (max(zh(I,j), z_clear) * I_Hbbl(I,j))
            botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

            hvel(I,j,k) = (1. - botfn) * h_arith(I,j) + botfn * h_harm(I,j)
            dz_vel(I,j,k) = (1. - botfn) * dz_arith(I,j) + botfn * dz_harm(I,j,k)
          endif
        endif
      endif ; enddo ; enddo
    endif

    if (CS%use_GL90_in_SSW) then
      ! The following block calculates the normalized height above the GL90 BBL
      ! (z_i_gl90), using a harmonic mean between layer thicknesses. For the
      ! GL90 BBL we use simply a constant (Hbbl_gl90). The purpose is that the
      ! GL90 coupling coefficient is zeroed out within Hbbl_gl90, to ensure
      ! that no momentum gets fluxed into vanished layers. The scheme is not
      ! sensitive to the exact value of Hbbl_gl90, as long as it is in a
      ! reasonable range (~1-20 m): large enough to capture vanished layers
      ! over topography, small enough to not contaminate the interior.

      do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
        z_i_gl90(I,j,k) = z_i_gl90(I,j,k+1) + dz_harm(I,j,k) * I_Hbbl_gl90(I,j)
      endif ; enddo ; enddo
    endif
  enddo

  call find_coupling_coef(a_cpl, dz_vel, do_i, dz_harm, bbl_thick, kv_bbl, z_i, &
      h_ml, dt, G, GV, US, CS, visc, Ustar_2d, tv, work_on_u=.true., OBC=OBC)

  if (allocated(hML_u)) then
    do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
      hML_u(I,j) = h_ml(I,j)
    endif ; enddo ; enddo
  endif

  if (CS%use_GL90_in_SSW) then
    call find_coupling_coef_gl90(a_cpl_gl90, dz_vel, do_i, z_i_gl90, G, GV, &
        CS, VarMix, work_on_u=.true.)
  endif

  do_any_shelf = .false.
  if (associated(forces%frac_shelf_u)) then
    do j=js,je ; do I=Isq,Ieq
      CS%a1_shelf_u(I,j) = 0.
      do_i_shelf(I,j) = do_i(I,j) .and. forces%frac_shelf_u(I,j) > 0.
    enddo ; enddo
    do_any_shelf = any(do_i_shelf)

    if (do_any_shelf) then
      if (.not. CS%harmonic_visc) then
        do j=js,je ; do I=Isq,Ieq ; if (do_i_shelf(I,j)) then
          zh(I,j) = 0.
          Ztop_min(I,j) = min(zcol(i,j), zcol(i+1,j))
          I_HTbl(I,j) = 1. / (visc%tbl_thick_shelf_u(I,j) + dz_neglect)
        endif ; enddo ; enddo
      endif

      do k=1,nz
        if (CS%harmonic_visc) then
          do j=js,je ; do I=Isq,Ieq
            hvel_shelf(I,j,k) = hvel(I,j,k)
            dz_vel_shelf(I,j,k) = dz_vel(I,j,k)
          enddo ; enddo
        else
          ! Find upwind-biased thickness near the surface.
          ! (Perhaps this needs to be done more carefully, via find_eta.)
          do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
            h_harm(I,j) = 2. * h(i,j,k) * h(i+1,j,k) &
                / (h(i,j,k) + h(i+1,j,k) + h_neglect)
            h_arith(I,j) = 0.5 * (h(i+1,j,k) + h(i,j,k))
            h_delta(I,j) = h(i+1,j,k) - h(i,j,k)
            dz_arith(I,j) = 0.5 * (dz(i+1,j,k) + dz(i,j,k))
          endif ; enddo ; enddo

          if (associated(OBC)) then
            if (OBC%u_E_OBCs_on_PE) then
              do j=js_E_OBC,je_E_OBC ; do I=Is_E_OBC,Ie_E_OBC
                if (do_i(I,j) .and. OBC%segnum_u(I,j) > 0) then
                  h_harm(I,j) = h(i,j,k)
                  h_arith(I,j) = h(i,j,k)
                  h_delta(I,j) = 0.
                  dz_arith(I,j) = dz(i,j,k)
                endif
              enddo ; enddo
            endif

            if (OBC%u_W_OBCs_on_PE) then
              do j=js_W_OBC,je_W_OBC ; do I=Is_W_OBC,Ie_W_OBC
                if (do_i(I,j) .and. OBC%segnum_u(I,j) < 0) then
                  h_harm(I,j) = h(i+1,j,k)
                  h_arith(I,j) = h(i+1,j,k)
                  h_delta(I,j) = 0.
                  dz_arith(I,j) = dz(i+1,j,k)
                endif
              enddo ; enddo
            endif
          endif

          do j=js,je ; do i=Isq,Ieq+1
            zcol(i,j) = zcol(i,j) - dz(i,j,k)
          enddo ; enddo

          do j=js,je ; do I=Isq,Ieq ; if (do_i_shelf(I,j)) then
            zh(I,j) = zh(I,j) + dz_harm(I,j,k)

            hvel_shelf(I,j,k) = hvel(I,j,k)
            dz_vel_shelf(I,j,k) = dz_vel(I,j,k)

            if (u(I,j,k) * h_delta(I,j) > 0) then
              if (zh(I,j) * I_HTbl(I,j) < CS%harm_BL_val) then
                hvel_shelf(I,j,k) = min(hvel(I,j,k), h_harm(I,j))
                dz_vel_shelf(I,j,k) = min(dz_vel(I,j,k), dz_harm(I,j,k))
              else
                z2_wt = 1.
                if (zh(I,j) * I_HTbl(I,j) < 2. * CS%harm_BL_val) &
                  z2_wt = max(0., min(1., zh(I,j) * I_HTbl(I,j) * I_valBL - 1.))

                z2 = z2_wt * (max(zh(I,j), Ztop_min(I,j) - min(zcol(i,j),zcol(i+1,j))) * I_HTbl(I,j))
                topfn = 1. / (1. + 0.09 * z2**6)

                h_arith(I,j) = 0.5 * (h(i+1,j,k) + h(i,j,k))
                dz_arith(I,j) = 0.5 * (dz(i+1,j,k) + dz(i,j,k))

                hvel_shelf(I,j,k) = min(hvel(I,j,k), (1. - topfn) * h_arith(I,j) + topfn * h_harm(I,j))
                dz_vel_shelf(I,j,k) = min(dz_vel(I,j,k), (1. - topfn) * dz_arith(I,j) + topfn * dz_harm(I,j,k))
              endif
            endif
          endif ; enddo ; enddo
        endif
      enddo

      call find_coupling_coef(a_shelf, dz_vel_shelf, do_i_shelf, dz_harm, &
          bbl_thick, kv_bbl, z_i, h_ml, dt, G, GV, US, CS, visc, Ustar_2d, &
          tv, work_on_u=.true., OBC=OBC, shelf=.true.)

      do j=js,je ; do I=Isq,Ieq ; if (do_i_shelf(I,j)) then
        CS%a1_shelf_u(I,j) = a_shelf(I,j,1)
      endif ; enddo ; enddo
    endif
  endif

  if (do_any_shelf) then
    if (CS%use_GL90_in_SSW) then
      do K=1,nz+1
        do j=js,je ; do I=Isq,Ieq
          if (do_i_shelf(I,j)) then
            CS%a_u(I,j,K) = min(a_cpl_max, (forces%frac_shelf_u(I,j) * a_shelf(I,j,K) + &
                                           (1. - forces%frac_shelf_u(I,j)) * a_cpl(I,j,K)) + a_cpl_gl90(I,j,K))

            ! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
            ! CS%a_u(I,j,K) = min(a_cpl_max, forces%frac_shelf_u(I,j) * max(a_shelf(I,j,K), a_cpl(I,j,K)) + &
            !                                (1. - forces%frac_shelf_u(I,j)) * a_cpl(I,j,K))
            CS%a_u_gl90(I,j,K) = min(a_cpl_max, a_cpl_gl90(I,j,K))
          elseif (do_i(I,j)) then
            CS%a_u(I,j,K) = min(a_cpl_max, a_cpl(I,j,K) + a_cpl_gl90(I,j,K))
            CS%a_u_gl90(I,j,K) = min(a_cpl_max, a_cpl_gl90(I,j,K))
          endif
        enddo ; enddo
      enddo
    else
      do K=1,nz+1
        do j=js,je ; do I=Isq,Ieq
          if (do_i_shelf(I,j)) then
            CS%a_u(I,j,K) = min(a_cpl_max, (forces%frac_shelf_u(I,j) * a_shelf(I,j,K) + &
                                           (1. - forces%frac_shelf_u(I,j)) * a_cpl(I,j,K)))

            ! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
            ! CS%a_u(I,j,K) = min(a_cpl_max, forces%frac_shelf_u(I,j) * max(a_shelf(I,j,K), a_cpl(I,j,K)) + &
            !                                (1. - forces%frac_shelf_u(I,j)) * a_cpl(I,j,K))
          elseif (do_i(I,j)) then
            CS%a_u(I,j,K) = min(a_cpl_max, a_cpl(I,j,K))
          endif
        enddo ; enddo
      enddo
    endif

    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        if (do_i_shelf(I,j)) then
          ! Should we instead take the inverse of the average of the inverses?
          CS%h_u(I,j,k) = forces%frac_shelf_u(I,j) * hvel_shelf(I,j,k) &
              + (1. - forces%frac_shelf_u(I,j)) * hvel(I,j,k) + h_neglect
        elseif (do_i(I,j)) then
          CS%h_u(I,j,k) = hvel(I,j,k) + h_neglect
        endif
      enddo ; enddo
    enddo
  else
    if (CS%use_GL90_in_SSW) then
      do K=1,nz+1
        do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
          a_cpl(I,j,K) = a_cpl(I,j,K) + a_cpl_gl90(I,j,K)
        endif; enddo ; enddo
      enddo

      do K=1,nz+1
        do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
          CS%a_u_gl90(I,j,K) = min(a_cpl_max, a_cpl_gl90(I,j,K))
        endif; enddo ; enddo
      enddo
    endif

    do K=1,nz+1
      do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
        CS%a_u(I,j,K) = min(a_cpl_max, a_cpl(I,j,K))
      endif; enddo ; enddo
    enddo

    do k=1,nz
      do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
        CS%h_u(I,j,k) = hvel(I,j,k) + h_neglect
      endif; enddo ; enddo
    enddo
  endif

  ! Diagnose total Kv at u-points
  if (CS%id_Kv_u > 0) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
        Kv_u(I,j,k) = 0.5 * (CS%a_u(I,j,K) + CS%a_u(I,j,K+1)) * CS%h_u(I,j,k)
      endif ; enddo ; enddo
    enddo
  endif

  ! Diagnose GL90 Kv at u-points
  if (CS%id_Kv_gl90_u > 0) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq ; if (do_i(I,j)) then
        Kv_gl90_u(I,j,k) = 0.5 * (CS%a_u_gl90(I,j,K) + CS%a_u_gl90(I,j,K+1)) * CS%h_u(I,j,k)
      endif ; enddo ; enddo
    enddo
  endif

  ! Now work on v-points.

  ! Force IPO optimizations (e.g. Intel)
  ij = touch_ij(i,j)

  do J=Jsq,Jeq ; do i=is,ie
    do_i(i,J) = G%mask2dCv(i,J) > 0.
  enddo ; enddo

  if (CS%bottomdraglaw) then
    do J=Jsq,Jeq ; do i=is,ie ; if(do_i(i,J)) then
      kv_bbl(i,J) = visc%Kv_bbl_v(i,J)
      bbl_thick(i,J) = visc%bbl_thick_v(i,J) + dz_neglect
      I_Hbbl(i,J) = 1. / bbl_thick(i,J)
    endif ; enddo ; enddo
  endif

  do J=Jsq,Jeq ; do i=is,ie
    Dmin(i,J) = min(G%bathyT(i,j), G%bathyT(i,j+1))
    zi_dir(i,J) = 0
  enddo ; enddo

  ! Project thickness outward across OBCs using a zero-gradient condition.
  if (associated(OBC)) then
    if (OBC%v_N_OBCs_on_PE) then
      do J=Js_N_OBC,Je_N_OBC ; do i=is_N_OBC,ie_N_OBC
        if (do_i(i,J) .and. OBC%segnum_v(i,J) > 0) then
          Dmin(I,J) = G%bathyT(i,j)
          zi_dir(I,J) = -1
        endif
      enddo ; enddo
    endif

    if (OBC%v_S_OBCs_on_PE) then
      do J=Js_S_OBC,Je_S_OBC ; do i=is_S_OBC,ie_S_OBC
        if (do_i(i,J) .and. OBC%segnum_v(i,J) < 0) then
          Dmin(i,J) = G%bathyT(i,j+1)
          zi_dir(i,J) = 1
        endif
      enddo ; enddo
    endif
  endif

  do J=Jsq,Jeq ; do i=is,ie
    z_i(i,J,nz+1) = 0.
  enddo ; enddo

  if (.not. CS%harmonic_visc) then
    do J=Jsq,Jeq ; do i=is,ie
      zh(i,J) = 0.
    enddo ; enddo

    do J=Jsq,Jeq+1 ; do i=is,ie
      zcol(i,j) = -G%bathyT(i,j)
    enddo ; enddo
  endif

  if (CS%use_GL90_in_SSW) then
    do j=Jsq,Jeq ; do i=is,ie
      z_i_gl90(i,J,nz+1) = 0.
    enddo ; enddo
  endif

  do k=nz,1,-1
    do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
      h_harm(i,J) = 2. * h(i,j,k) * h(i,j+1,k) / (h(i,j,k) + h(i,j+1,k) + h_neglect)
      h_arith(i,J) = 0.5 * (h(i,j+1,k) + h(i,j,k))
      h_delta(i,J) = h(i,j+1,k) - h(i,j,k)
      dz_harm(i,J,k) = 2. * dz(i,j,k) * dz(i,j+1,k) / (dz(i,j,k) + dz(i,j+1,k) + dz_neglect)
      dz_arith(i,J) = 0.5 * (dz(i,j+1,k) + dz(i,j,k))
    endif ; enddo ; enddo

    ! Project thickness outward across OBCs using a zero-gradient condition.
    if (associated(OBC)) then
      if (OBC%v_N_OBCs_on_PE) then
        do J=Js_N_OBC,Je_N_OBC ; do i=is_N_OBC,ie_N_OBC
          if (do_i(i,J) .and. OBC%segnum_v(i,J) > 0) then
            h_harm(i,J) = h(i,j,k)
            h_arith(i,J) = h(i,j,k)
            h_delta(i,J) = 0.
            dz_harm(i,J,k) = dz(i,j,k)
            dz_arith(i,J) = dz(i,j,k)
          endif
        enddo ; enddo
      endif

      if (OBC%v_S_OBCs_on_PE) then
        do J=Js_S_OBC,Je_S_OBC ; do i=is_S_OBC,ie_S_OBC
          if (do_i(i,J) .and. OBC%segnum_v(i,J) < 0) then
            h_harm(i,J) = h(i,j+1,k)
            h_arith(i,J) = h(i,j+1,k)
            h_delta(i,J) = 0.
            dz_harm(i,J,k) = dz(i,j+1,k)
            dz_arith(i,J) = dz(i,j+1,k)
          endif
        enddo ; enddo
      endif
    endif

    if (CS%harmonic_visc) then
      ! The following block calculates the thicknesses at velocity grid points
      ! for the vertical viscosity (hvel and dz_vel).  Near the bottom an
      ! upwind biased thickness is used to control the effect of spurious
      ! Montgomery potential gradients at the bottom where nearly massless
      ! layers ride over the topography.

      do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
        hvel(i,J,k) = h_harm(i,J)
        dz_vel(i,J,k) = dz_harm(i,J,k)

        if (v(i,J,k) * h_delta(i,J) < 0) then
          z2 = z_i(i,J,k+1)
          botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

          hvel(i,J,k) = (1. - botfn) * h_harm(i,J) + botfn * h_arith(i,J)
          dz_vel(i,J,k) = (1. - botfn) * dz_harm(i,J,k) + botfn * dz_arith(i,J)
        endif

        z_i(i,J,k) = z_i(i,J,k+1) + dz_harm(i,J,k)*I_Hbbl(i,J)
      endif ; enddo ; enddo
    else ! Not harmonic_visc
      do J=Jsq,Jeq+1 ; do i=is,ie
        zcol(i,j) = zcol(i,j) + dz(i,j,k)
      enddo ; enddo

      do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
        zh(i,J) = zh(i,J) + dz_harm(i,J,k)

        z_clear = max(zcol(i,j), zcol(i,j+1)) + Dmin(i,J)
        if (zi_dir(i,J) < 0) z_clear = zcol(i,j) + Dmin(i,J)
        if (zi_dir(i,J) > 0) z_clear = zcol(i,j+1) + Dmin(i,J)

        z_i(i,J,k) = max(zh(i,J), z_clear) * I_Hbbl(i,J)

        hvel(i,J,k) = h_arith(i,J)
        dz_vel(i,J,k) = dz_arith(i,J)

        if (v(i,J,k) * h_delta(i,J) > 0) then
          if (zh(i,J) * I_Hbbl(i,J) < CS%harm_BL_val) then
            hvel(i,J,k) = h_harm(i,J)
            dz_vel(i,J,k) = dz_harm(i,J,k)
          else
            z2_wt = 1.
            if (zh(i,J) * I_Hbbl(i,J) < 2. * CS%harm_BL_val) &
              z2_wt = max(0., min(1., zh(i,J) * I_Hbbl(i,J) * I_valBL - 1.))

            ! TODO: should z_clear be used here?
            z2 = z2_wt * (max(zh(i,J), max(zcol(i,j), zcol(i,j+1)) + Dmin(i,J)) * I_Hbbl(i,J))
            botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

            hvel(i,J,k) = (1. - botfn) * h_arith(i,J) + botfn * h_harm(i,J)
            dz_vel(i,J,k) = (1. - botfn) * dz_arith(i,J) + botfn * dz_harm(i,J,k)
          endif
        endif
      endif ; enddo ; enddo
    endif

    if (CS%use_GL90_in_SSW) then
      ! The following block calculates the normalized height above the GL90 BBL
      ! (z_i_gl90), using a harmonic mean between layer thicknesses. For the
      ! GL90 BBL we use simply a constant (Hbbl_gl90). The purpose is that the
      ! GL90 coupling coefficient is zeroed out within Hbbl_gl90, to ensure
      ! that no momentum gets fluxed into vanished layers. The scheme is not
      ! sensitive to the exact value of Hbbl_gl90, as long as it is in a
      ! reasonable range (~1-20 m): large enough to capture vanished layers
      ! over topography, small enough to not contaminate the interior.

      do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
        z_i_gl90(i,J,k) = z_i_gl90(i,J,k+1) + dz_harm(i,J,k) * I_Hbbl_gl90(i,J)
      endif ; enddo ; enddo
    endif
  enddo

  call find_coupling_coef(a_cpl, dz_vel, do_i, dz_harm, bbl_thick, kv_bbl, z_i, &
      h_ml, dt, G, GV, US, CS, visc, Ustar_2d, tv, work_on_u=.false., OBC=OBC)

  if ( allocated(hML_v)) then
    do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
      hML_v(i,J) = h_ml(i,J)
    endif ; enddo ; enddo
  endif

  if (CS%use_GL90_in_SSW) then
    call find_coupling_coef_gl90(a_cpl_gl90, dz_vel, do_i, z_i_gl90, G, GV, &
        CS, VarMix, work_on_u=.false.)
  endif

  do_any_shelf = .false.
  if (associated(forces%frac_shelf_v)) then
    do J=Jsq,Jeq ; do i=is,ie
      CS%a1_shelf_v(i,J) = 0.
      do_i_shelf(i,J) = do_i(i,J) .and. forces%frac_shelf_v(i,J) > 0.
    enddo ; enddo
    do_any_shelf = any(do_i_shelf)

    if (do_any_shelf) then
      ! Initialize non-harmonic depths
      if (.not. CS%harmonic_visc) then
        do J=Jsq,Jeq ; do i=is,ie ; if (do_i_shelf(i,J)) then
          zh(i,J) = 0.
          Ztop_min(i,J) = min(zcol(i,j), zcol(i,j+1))
          I_HTbl(i,J) = 1. / (visc%tbl_thick_shelf_v(i,J) + dz_neglect)
        endif ; enddo ; enddo
      endif

      do k=1,nz
        if (CS%harmonic_visc) then
          do J=Jsq,Jeq ; do i=is,ie
            hvel_shelf(i,J,k) = hvel(i,J,k)
            dz_vel_shelf(i,J,k) = dz_vel(i,J,k)
          enddo ; enddo
        else
          ! Find upwind-biased thickness near the surface.
          ! Perhaps this needs to be done more carefully, via find_eta.
          do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
            h_harm(i,J) = 2. * h(i,j,k) * h(i,j+1,k) &
                / (h(i,j,k) + h(i,j+1,k) + h_neglect)
            h_arith(i,J) = 0.5 * (h(i,j+1,k) + h(i,j,k))
            h_delta(i,J) = h(i,j+1,k) - h(i,j,k)
            dz_arith(i,J) = 0.5 * (dz(i,j+1,k) + dz(i,j,k))
          endif ; enddo ; enddo

          ! Project thickness outward across OBCs using a zero-gradient condition.
          if (associated(OBC)) then
            if (OBC%v_N_OBCs_on_PE) then
              do J=Js_N_OBC,Je_N_OBC ; do i=is_N_OBC,ie_N_OBC
                if (do_i(i,J) .and. OBC%segnum_v(i,J) > 0) then
                  h_harm(i,J) = h(i,j,k)
                  h_arith(i,J) = h(i,j,k)
                  h_delta(i,J) = 0.
                  dz_arith(i,J) = dz(i,j,k)
                endif
              enddo ; enddo
            endif

            if (OBC%v_S_OBCs_on_PE) then
              do J=Js_S_OBC,Je_S_OBC ; do i=is_S_OBC,ie_S_OBC
                if (do_i(i,J) .and. OBC%segnum_v(i,J) < 0) then
                  h_harm(i,J) = h(i,j+1,k)
                  h_arith(i,J) = h(i,j+1,k)
                  h_delta(i,J) = 0.
                  dz_arith(i,J) = dz(i,j+1,k)
                endif
              enddo ; enddo
            endif
          endif

          do J=Jsq,Jeq+1 ; do i=is,ie
            zcol(i,j) = zcol(i,j) - dz(i,j,k)
          enddo ; enddo

          do J=Jsq,Jeq ; do i=is,je ; if (do_i_shelf(i,J)) then
            zh(i,J) = zh(i,J) + dz_harm(i,J,k)

            hvel_shelf(i,J,k) = hvel(i,J,k)
            dz_vel_shelf(i,J,k) = dz_vel(i,J,k)

            if (v(i,J,k) * h_delta(i,J) > 0.) then
              if (zh(i,J) * I_HTbl(i,J) < CS%harm_BL_val) then
                hvel_shelf(i,J,k) = min(hvel(i,J,k), h_harm(i,J))
                dz_vel_shelf(i,J,k) = min(dz_vel(i,J,k), dz_harm(i,J,k))
              else
                z2_wt = 1.
                if (zh(i,J) * I_HTbl(i,J) < 2. * CS%harm_BL_val) &
                  z2_wt = max(0., min(1., zh(i,J) * I_HTbl(i,J) * I_valBL - 1.))

                z2 = z2_wt * (max(zh(i,J), Ztop_min(i,J) - min(zcol(i,j), zcol(i,j+1))) * I_HTbl(i,J))
                topfn = 1. / (1. + 0.09 * z2**6)

                h_arith(i,J) = 0.5 * (h(i,j+1,k) + h(i,j,k))
                dz_arith(i,J) = 0.5 * (dz(i,j+1,k) + dz(i,j,k))

                hvel_shelf(i,J,k) = min(hvel(i,J,k), (1. - topfn) * h_arith(i,J) + topfn * h_harm(i,J))
                dz_vel_shelf(i,J,k) = min(dz_vel(i,J,k), (1. - topfn) * dz_arith(i,J) + topfn * dz_harm(i,J,k))
              endif
            endif
          endif ; enddo ; enddo
        endif
      enddo

      call find_coupling_coef(a_shelf, dz_vel_shelf, do_i_shelf, dz_harm, &
          bbl_thick, kv_bbl, z_i, h_ml, dt, G, GV, US, CS, visc, Ustar_2d, &
          tv, work_on_u=.false., OBC=OBC, shelf=.true.)

      do J=Jsq,Jeq ; do i=is,ie ; if (do_i_shelf(i,J)) then
        CS%a1_shelf_v(i,J) = a_shelf(i,J,1)
      endif ; enddo ; enddo
    endif
  endif

  if (do_any_shelf) then
    if (CS%use_GL90_in_SSW) then
      do K=1,nz+1
        do J=Jsq,Jeq ; do i=is,ie
          if (do_i_shelf(i,J)) then
            CS%a_v(i,J,K) = min(a_cpl_max, (forces%frac_shelf_v(i,J) * a_shelf(i,J,k) + &
                                           (1. - forces%frac_shelf_v(i,J)) * a_cpl(i,J,K)) + a_cpl_gl90(i,J,K))
            ! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
            ! CS%a_v(i,J,K) = min(a_cpl_max, forces%frac_shelf_v(i,J) * max(a_shelf(i,J,K), a_cpl(i,J,K)) + &
            !                                (1. - forces%frac_shelf_v(i,J)) * a_cpl(i,J,K))
            CS%a_v_gl90(i,J,K) = min(a_cpl_max, a_cpl_gl90(i,J,K))
          elseif (do_i(i,J)) then
            CS%a_v(i,J,K) = min(a_cpl_max, a_cpl(i,J,K) + a_cpl_gl90(i,J,K))
            CS%a_v_gl90(i,J,K) = min(a_cpl_max, a_cpl_gl90(i,J,K))
          endif
        enddo ; enddo
      enddo
    else
      do K=1,nz+1
        do J=Jsq,Jeq ; do i=is,ie
          if (do_i_shelf(i,J)) then
            CS%a_v(i,J,K) = min(a_cpl_max, (forces%frac_shelf_v(i,J) * a_shelf(i,J,k) + &
                                           (1. - forces%frac_shelf_v(i,J)) * a_cpl(i,J,K)))
            ! This is Alistair's suggestion, but it destabilizes the model. I do not know why. RWH
            ! CS%a_v(i,J,K) = min(a_cpl_max, forces%frac_shelf_v(i,J) * max(a_shelf(i,J,K), a_cpl(i,J,K)) + &
            !                                (1. - forces%frac_shelf_v(i,J)) * a_cpl(i,J,K))
          elseif (do_i(i,J)) then
            CS%a_v(i,J,K) = min(a_cpl_max, a_cpl(i,J,K))
          endif
        enddo ; enddo
      enddo
    endif

    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie
        if (do_i_shelf(i,J)) then
          ! Should we instead take the inverse of the average of the inverses?
          CS%h_v(i,J,k) = forces%frac_shelf_v(i,J)  * hvel_shelf(i,J,k) + &
                     (1. - forces%frac_shelf_v(i,J)) * hvel(i,J,k) + h_neglect
        elseif (do_i(i,J)) then
          CS%h_v(i,J,k) = hvel(i,J,k) + h_neglect
        endif
      enddo ; enddo
    enddo
  else
    if (CS%use_GL90_in_SSW) then
      do K=1,nz+1
        do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
          a_cpl(i,J,K) = a_cpl(i,J,K) + a_cpl_gl90(i,J,K)
        endif ; enddo ; enddo
      enddo

      do K=1,nz+1
        do J=Jsq,Jeq; do i=is,ie ; if (do_i(i,J)) then
          CS%a_v_gl90(i,J,K) = min(a_cpl_max, a_cpl_gl90(i,J,K))
        endif ; enddo ; enddo
      enddo
    endif

    do K=1,nz+1
      do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
        CS%a_v(i,J,K) = min(a_cpl_max, a_cpl(i,J,K))
      endif ; enddo ; enddo
    enddo

    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
        CS%h_v(i,J,k) = hvel(i,J,k) + h_neglect
      endif; enddo ; enddo
    enddo
  endif

  ! Diagnose total Kv at v-points
  if (CS%id_Kv_v > 0) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
        Kv_v(i,J,k) = 0.5 * (CS%a_v(i,J,K)+CS%a_v(i,J,K+1)) * CS%h_v(i,J,k)
      endif ; enddo ; enddo
    enddo
  endif

  ! Diagnose GL90 Kv at v-points
  if (CS%id_Kv_gl90_v > 0) then
    do k=1,nz
      do J=Jsq,Jeq ; do i=is,ie ; if (do_i(i,J)) then
        Kv_gl90_v(i,J,k) = 0.5 * (CS%a_v_gl90(i,J,K)+CS%a_v_gl90(i,J,K+1)) * CS%h_v(i,J,k)
      endif ; enddo ; enddo
    enddo
  endif

  if (CS%debug) then
    call uvchksum("vertvisc_coef h_[uv]", CS%h_u, CS%h_v, G%HI, haloshift=0, &
                  unscale=GV%H_to_m, scalar_pair=.true.)
    call uvchksum("vertvisc_coef a_[uv]", CS%a_u, CS%a_v, G%HI, haloshift=0, &
                  unscale=GV%H_to_m*US%s_to_T, scalar_pair=.true.)
    if (allocated(hML_u) .and. allocated(hML_v)) &
      call uvchksum("vertvisc_coef hML_[uv]", hML_u, hML_v, G%HI, &
                    haloshift=0, unscale=US%Z_to_m, scalar_pair=.true.)
  endif

! Offer diagnostic fields for averaging.
  if (query_averaging_enabled(CS%diag)) then
    if (associated(visc%Kv_slow) .and. (CS%id_Kv_slow > 0)) &
        call post_data(CS%id_Kv_slow, visc%Kv_slow, CS%diag)
    if (CS%id_Kv_u > 0) call post_data(CS%id_Kv_u, Kv_u, CS%diag)
    if (CS%id_Kv_v > 0) call post_data(CS%id_Kv_v, Kv_v, CS%diag)
    if (CS%id_Kv_gl90_u > 0) call post_data(CS%id_Kv_gl90_u, Kv_gl90_u, CS%diag)
    if (CS%id_Kv_gl90_v > 0) call post_data(CS%id_Kv_gl90_v, Kv_gl90_v, CS%diag)
    if (CS%id_au_vv > 0) call post_data(CS%id_au_vv, CS%a_u, CS%diag)
    if (CS%id_av_vv > 0) call post_data(CS%id_av_vv, CS%a_v, CS%diag)
    if (CS%id_au_gl90_vv > 0) call post_data(CS%id_au_gl90_vv, CS%a_u_gl90, CS%diag)
    if (CS%id_av_gl90_vv > 0) call post_data(CS%id_av_gl90_vv, CS%a_v_gl90, CS%diag)
    if (CS%id_h_u > 0) call post_data(CS%id_h_u, CS%h_u, CS%diag)
    if (CS%id_h_v > 0) call post_data(CS%id_h_v, CS%h_v, CS%diag)
    if (CS%id_hML_u > 0) call post_data(CS%id_hML_u, hML_u, CS%diag)
    if (CS%id_hML_v > 0) call post_data(CS%id_hML_v, hML_v, CS%diag)
  endif

  if (allocated(hML_u)) deallocate(hML_u)
  if (allocated(hML_v)) deallocate(hML_v)

end subroutine vertvisc_coef


!> Calculate the 'coupling coefficient' (a_cpl) at the interfaces.
!! If BOTTOMDRAGLAW is defined, the minimum of Hbbl and half the adjacent
!! layer thicknesses are used to calculate a_cpl near the bottom.
subroutine find_coupling_coef(a_cpl, hvel, do_i, h_harm, bbl_thick, kv_bbl, z_i, h_ml, &
                              dt, G, GV, US, CS, visc, Ustar_2d, tv, work_on_u, OBC, shelf)
  type(ocean_grid_type),     intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),   intent(in)  :: GV !< Ocean vertical grid structure
  type(unit_scale_type),     intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), &
                             intent(out) :: a_cpl !< Coupling coefficient across interfaces [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)), &
                             intent(in)  :: hvel !< Distance between interfaces at velocity points [Z ~> m]
  logical, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)  :: do_i !< If true, determine coupling coefficient for a column
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)), &
                             intent(in)  :: h_harm !< Harmonic mean of thicknesses around a velocity
                                                   !! grid point [Z ~> m]
  real, dimension(SZIB_(G),SZJB_(G)), intent(in)  :: bbl_thick !< Bottom boundary layer thickness [Z ~> m]
  real, dimension(SZIB_(G),SZJB_(G)), intent(in)  :: kv_bbl !< Bottom boundary layer viscosity, exclusive of
                                                   !! any depth-dependent contributions from
                                                   !! visc%Kv_shear [H Z T-1 ~> m2 s-1 or Pa s]
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), &
                             intent(in)  :: z_i  !< Estimate of interface heights above the bottom,
                                                 !! normalized by the bottom boundary layer thickness [nondim]
  real, dimension(SZIB_(G),SZJB_(G)), intent(out) :: h_ml !< Mixed layer depth [Z ~> m]
  real,                      intent(in)  :: dt   !< Time increment [T ~> s]
  type(vertvisc_CS),         intent(in)  :: CS   !< Vertical viscosity control structure
  type(vertvisc_type),       intent(in)  :: visc !< Structure containing viscosities and bottom drag
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)  :: Ustar_2d !< The wind friction velocity, calculated using
                                                 !! the Boussinesq reference density or the
                                                 !! time-evolving surface density in non-Boussinesq
                                                 !! mode [Z T-1 ~> m s-1]
  type(thermo_var_ptrs),     intent(in)  :: tv   !< A structure containing pointers to any available
                                                 !! thermodynamic fields.
  logical,                   intent(in)  :: work_on_u !< If true, u-points are being calculated,
                                                  !! otherwise they are v-points
  type(ocean_OBC_type),      pointer     :: OBC   !< Open boundary condition structure
  logical,         optional, intent(in)  :: shelf !< If present and true, use a surface boundary
                                                  !! condition appropriate for an ice shelf.

  ! Local variables

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    u_star, &   ! ustar at a velocity point [Z T-1 ~> m s-1]
    tau_mag, &  ! The magnitude of the wind stress at a velocity point including gustiness [H Z T-2 ~> m2 s-2 or Pa]
    absf, &     ! The average of the neighboring absolute values of f [T-1 ~> s-1].
    rho_av1, &  ! The harmonic mean surface layer density at velocity points [R ~> kg m-3]
    z_t, &      ! The distance from the top, sometimes normalized
                ! by Hmix, [Z ~> m] or [nondim].
    kv_TBL, &   ! The viscosity in a top boundary layer under ice [H Z T-1 ~> m2 s-1 or Pa s]
    tbl_thick, &! The thickness of the top boundary layer [Z ~> m]
    Kv_add, &   ! A viscosity to add [H Z T-1 ~> m2 s-1 or Pa s]
    Kv_tot      ! The total viscosity at an interface [H Z T-1 ~> m2 s-1 or Pa s]
  integer, dimension(SZIB_(G),SZJB_(G)) :: &
    nk_in_ml      ! The index of the deepest interface in the mixed layer.
  real :: h_shear ! The distance over which shears occur [Z ~> m].
  real :: dhc     ! The distance between the center of adjacent layers [Z ~> m].
  real :: visc_ml ! The mixed layer viscosity [H Z T-1 ~> m2 s-1 or Pa s].
  real :: I_Hmix  ! The inverse of the mixed layer thickness [Z-1 ~> m-1].
  real :: a_ml    ! The layer coupling coefficient across an interface in
                  ! the mixed layer [H T-1 ~> m s-1 or Pa s m-1].
  real :: a_floor ! A lower bound on the layer coupling coefficient across an interface in
                  ! the mixed layer [H T-1 ~> m s-1 or Pa s m-1].
  real :: I_amax  ! The inverse of the maximum coupling coefficient [T H-1 ~> s m-1 or s m2 kg-1].
  real :: temp1   ! A temporary variable [Z2 ~> m2]
  real :: ustar2_denom ! A temporary variable in the surface boundary layer turbulence
                  ! calculations [H Z-1 T-1 ~> s-1 or kg m-3 s-1]
  real :: h_neglect ! A vertical distance that is so small it is usually lost
                  ! in roundoff and can be neglected [Z ~> m].
  real :: z2      ! A copy of z_i [nondim]
  real :: botfn   ! A function that is 1 at the bottom and small far from it [nondim]
  real :: topfn   ! A function that is 1 at the top and small far from it [nondim]
  real :: kv_top  ! A viscosity associated with the top boundary layer [H Z T-1 ~> m2 s-1 or Pa s]
  logical :: do_shelf, do_OBCs, can_exit
  integer :: i, j, k
  integer :: is, ie, js, je
  integer :: nz, max_nk
  integer :: is_N_OBC, is_S_OBC, Is_E_OBC, Is_W_OBC, ie_N_OBC, ie_S_OBC, Ie_E_OBC, Ie_W_OBC
  integer :: js_N_OBC, js_S_OBC, Js_E_OBC, Js_W_OBC, je_N_OBC, je_S_OBC, Je_E_OBC, Je_W_OBC

  if (work_on_u) then
    Is = G%IscB ; Ie = G%IecB
    js = G%jsc ; je = G%jec
  else
    is = G%isc ; ie = G%iec
    Js = G%JscB ; Je = G%JecB
  endif
  nz = GV%ke

  h_neglect = GV%dZ_subroundoff

  if (CS%answer_date < 20190101) then
    !   The maximum coupling coefficient was originally introduced to avoid
    ! truncation error problems in the tridiagonal solver. Effectively, the 1e-10
    ! sets the maximum coupling coefficient increment to 1e10 m per timestep.
    I_amax = (1.0e-10*GV%H_to_m) * dt
  else
    I_amax = 0.0
  endif

  do_shelf = .false. ; if (present(shelf)) do_shelf = shelf

  do_OBCs = .false.
  if (associated(OBC)) then
    if (work_on_u) then
      do_OBCS = OBC%u_E_OBCs_on_PE .or. OBC%u_W_OBCs_on_PE
      Is_E_OBC = max(G%IscB, OBC%Is_u_E_obc) ; Ie_E_OBC = min(G%IecB, OBC%Ie_u_E_obc)
      Is_W_OBC = max(G%IscB, OBC%Is_u_W_obc) ; Ie_W_OBC = min(G%IecB, OBC%Ie_u_W_obc)
      js_E_OBC = max(G%jsc, OBC%js_u_E_obc) ; je_E_OBC = min(G%jec, OBC%je_u_E_obc)
      js_W_OBC = max(G%jsc, OBC%js_u_W_obc) ; je_W_OBC = min(G%jec, OBC%je_u_W_obc)
    else
      do_OBCS = OBC%v_N_OBCs_on_PE .or. OBC%v_S_OBCs_on_PE
      is_N_OBC = max(G%isc, OBC%is_v_N_obc) ; ie_N_OBC = min(G%iec, OBC%ie_v_N_obc)
      is_S_OBC = max(G%isc, OBC%is_v_S_obc) ; ie_S_OBC = min(G%iec, OBC%ie_v_S_obc)
      Js_N_OBC = max(G%JscB, OBC%Js_v_N_obc) ; Je_N_OBC = min(G%JecB, OBC%Je_v_N_obc)
      Js_S_OBC = max(G%JscB, OBC%Js_v_S_obc) ; Je_S_OBC = min(G%JecB, OBC%Je_v_S_obc)
    endif
  endif

  a_cpl(:,:,:) = 0.0
  h_ml(:,:) = 0.

  if (CS%Kvml_invZ2 > 0. .and. .not. do_shelf) then
    I_Hmix = 1. / (CS%Hmix + h_neglect)
    do j=js,je ; do i=is,ie
      z_t(i,j) = h_neglect * I_Hmix
    enddo ; enddo
  endif

  do K=2,nz
    do j=js,je ; do i=is,ie
      Kv_tot(i,j) = CS%Kv
    enddo ; enddo

    if (CS%Kvml_invZ2 > 0. .and. .not. do_shelf) then
      ! This is an older (vintage ~1997) way to prevent wind stresses from driving very
      ! large flows in nearly massless near-surface layers when there is not a physically-
      ! based surface boundary layer parameterization.  It does not have a plausible
      ! physical basis, and probably should not be used.
      do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
        z_t(i,j) = z_t(i,j) + h_harm(i,j,k-1) * I_Hmix
        Kv_tot(i,j) = CS%Kv + CS%Kvml_invZ2 / ((z_t(i,j)*z_t(i,j)) *  &
                 (1. + 0.09 * z_t(i,j) * z_t(i,j) * z_t(i,j) * z_t(i,j) * z_t(i,j) * z_t(i,j)))
      endif ; enddo ; enddo
    endif

    if (associated(visc%Kv_shear)) then
      ! Add in viscosities that are determined by physical processes that are handled in
      ! other modules, and which do not respond immediately to the changing layer thicknesses.
      ! These processes may include shear-driven mixing or contributions from some boundary
      ! layer turbulence schemes.  Other viscosity contributions that respond to the evolving
      ! layer thicknesses or the surface wind stresses are added later.
      if (work_on_u) then
        ! FIXME: Uppercase i?
        do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
          Kv_add(i,j) = 0.5*(visc%Kv_shear(i,j,k) + visc%Kv_shear(i+1,j,k))
        endif ; enddo ; enddo

        if (do_OBCs) then
          if (OBC%u_E_OBCs_on_PE) then
            do j=js_E_OBC,je_E_OBC ; do I=Is_E_OBC,Ie_E_OBC
              if (do_i(I,j) .and. OBC%segnum_u(I,j) > 0) then
                Kv_add(i,j) = visc%Kv_shear(i,j,k)
              endif
            enddo ; enddo
          endif

          if (OBC%u_W_OBCs_on_PE) then
            do j=js_W_OBC,je_W_OBC ; do I=Is_W_OBC,Ie_W_OBC
              if (do_i(I,j) .and. OBC%segnum_u(I,j) < 0) then
                Kv_add(i,j) = visc%Kv_shear(i+1,j,k)
              endif
            enddo ; enddo
          endif
        endif

        do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
          Kv_tot(i,j) = Kv_tot(i,j) + Kv_add(i,j)
        endif ; enddo ; enddo
      else
        ! FIXME: Uppercase j?
        do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
          Kv_add(i,j) = 0.5*(visc%Kv_shear(i,j,k) + visc%Kv_shear(i,j+1,k))
        endif ; enddo ; enddo

        if (do_OBCs) then
          if (OBC%v_N_OBCs_on_PE) then
            do J=Js_N_OBC,Je_N_OBC ; do i=is_N_OBC,ie_N_OBC
              if (do_i(i,J) .and. OBC%segnum_v(i,J) > 0) then
                Kv_add(i,j) = visc%Kv_shear(i,j,k)
              endif
            enddo ; enddo
          endif

          if (OBC%v_S_OBCs_on_PE) then
            do J=Js_S_OBC,Je_S_OBC ; do i=is_S_OBC,ie_S_OBC
              if (do_i(i,J) .and. OBC%segnum_v(i,J) < 0) then
                Kv_add(i,j) = visc%Kv_shear(i,j+1,k)
              endif
            enddo ; enddo
          endif
        endif

        do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
          Kv_tot(i,j) = Kv_tot(i,j) + Kv_add(i,j)
        endif ; enddo ; enddo
      endif
    endif

    if (associated(visc%Kv_shear_Bu)) then
      ! This is similar to what was done above, but for contributions coming from the corner
      ! (vorticity) points.  Because OBCs run through the faces and corners there is no need
      ! to further modify these viscosities here to take OBCs into account.
      if (work_on_u) then
        do J=Js,Je ; do I=Is,Ie ; If (do_i(i,j)) then
          Kv_tot(I,J) = Kv_tot(I,J) + 0.5 * (visc%Kv_shear_Bu(I,J-1,k) + visc%Kv_shear_Bu(I,J,k))
        endif ; enddo ; enddo
      else
        do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
          Kv_tot(i,j) = Kv_tot(i,j) + 0.5 * (visc%Kv_shear_Bu(I-1,J,k) + visc%Kv_shear_Bu(I,J,k))
        endif ; enddo ; enddo
      endif
    endif

    ! Set the viscous coupling coefficients, excluding surface mixed layer contributions
    ! for now, but including viscous bottom drag, working up from the bottom.
    if (CS%bottomdraglaw) then
      do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
        !    botfn determines when a point is within the influence of the bottom
        !  boundary layer, going from 1 at the bottom to 0 in the interior.
        z2 = z_i(i,j,k)
        botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

        Kv_tot(i,j) = Kv_tot(i,j) + (kv_bbl(i,j) - CS%Kv)*botfn
        dhc = 0.5 * (hvel(i,j,k) + hvel(i,j,k-1))
        if (dhc > bbl_thick(i,j)) then
          h_shear = ((1. - botfn) * dhc + botfn*bbl_thick(i,j)) + h_neglect
        else
          h_shear = dhc + h_neglect
        endif

        ! Calculate the coupling coefficients from the viscosities.
        a_cpl(i,j,K) = Kv_tot(i,j) / (h_shear + (I_amax * Kv_tot(i,j)))
      endif ; enddo ; enddo ! i & k loops
    elseif (abs(CS%Kv_extra_bbl) > 0.0) then
      ! There is a simple enhancement of the near-bottom viscosities, but no
      ! adjustment of the viscous coupling length scales to give a particular
      ! bottom stress.

      do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
        !    botfn determines when a point is within the influence of the bottom
        !  boundary layer, going from 1 at the bottom to 0 in the interior.
        z2 = z_i(i,j,k)
        botfn = 1. / (1. + 0.09 * z2 * z2 * z2 * z2 * z2 * z2)

        Kv_tot(i,j) = Kv_tot(i,j) + CS%Kv_extra_bbl*botfn
        h_shear = 0.5 * (hvel(i,j,k) + hvel(i,j,k-1) + h_neglect)

        ! Calculate the coupling coefficients from the viscosities.
        a_cpl(i,j,K) = Kv_tot(i,j) / (h_shear + I_amax * Kv_tot(i,j))
      endif ; enddo ; enddo ! i & k loops
    else
      ! Any near-bottom viscous enhancements were already incorporated into
      ! Kv_tot, and there is no adjustment of the viscous coupling length
      ! scales to give a particular bottom stress.

      do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
        h_shear = 0.5*(hvel(i,j,k) + hvel(i,j,k-1) + h_neglect)
        ! Calculate the coupling coefficients from the viscosities.
        a_cpl(i,j,K) = Kv_tot(i,j) / (h_shear + I_amax*Kv_tot(i,j))
      endif ; enddo ; enddo ! i & k loops
    endif
  enddo

  ! Assign the bottom coupling coefficients
  if (CS%bottomdraglaw) then
    do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
      dhc = hvel(i,j,nz)*0.5
      a_cpl(i,j,nz+1) = kv_bbl(i,j) / ((min(dhc, bbl_thick(i,j)) + h_neglect) + I_amax*kv_bbl(i,j))
    endif ; enddo ; enddo
  elseif (abs(CS%Kv_extra_bbl) > 0.0) then
    do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
      a_cpl(i,j,nz+1) = (CS%Kv + CS%Kv_extra_bbl) &
          / ((0.5 * hvel(i,j,nz) + h_neglect) + I_amax * (CS%Kv + CS%Kv_extra_bbl))
    endif ; enddo ; enddo
  else
    do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
      a_cpl(i,j,nz+1) = CS%Kv / ((0.5 * hvel(i,j,nz) + h_neglect) + I_amax * CS%Kv)
    endif ; enddo ; enddo
  endif

  ! Add surface intensified viscous coupling, either as a no-slip boundary condition under a
  ! rigid ice-shelf, or due to wind-stress driven surface boundary layer mixing that has not
  ! already been added via visc%Kv_shear.
  if (do_shelf) then
    ! Set the coefficients to include the no-slip surface stress.
    do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
      if (work_on_u) then
        kv_TBL(i,j) = visc%Kv_tbl_shelf_u(I,j)
        tbl_thick(i,j) = visc%tbl_thick_shelf_u(I,j) + h_neglect
      else
        kv_TBL(i,j) = visc%Kv_tbl_shelf_v(i,J)
        tbl_thick(i,j) = visc%tbl_thick_shelf_v(i,J) + h_neglect
      endif
      z_t(i,j) = 0.0

      ! If a_cpl(i,1) were not already 0, it would be added here.
      if (0.5*hvel(i,j,1) > tbl_thick(i,j)) then
        a_cpl(i,j,1) = kv_TBL(i,j) / (tbl_thick(i,j) + I_amax * kv_TBL(i,j))
      else
        a_cpl(i,j,1) = kv_TBL(i,j) &
            / ((0.5 * hvel(i,j,1) + h_neglect) + I_amax * kv_TBL(i,j))
      endif
    endif ; enddo ; enddo

    do K=2,nz
      do j=js,je ; do i=is,ie ;  if (do_i(i,j)) then
        z_t(i,j) = z_t(i,j) + hvel(i,j,k-1) / tbl_thick(i,j)
        topfn = 1. / (1. + 0.09 * z_t(i,j)**6)

        dhc = 0.5 * (hvel(i,j,k) + hvel(i,j,k-1))
        if (dhc > tbl_thick(i,j)) then
          h_shear = ((1. - topfn) * dhc + topfn * tbl_thick(i,j)) + h_neglect
        else
          h_shear = dhc + h_neglect
        endif

        kv_top = topfn * kv_TBL(i,j)
        a_cpl(i,j,K) = a_cpl(i,j,K) + kv_top / (h_shear + I_amax * kv_top)
      endif ; enddo ; enddo
    enddo
  elseif (CS%dynamic_viscous_ML .or. (GV%nkml>0) .or. CS%fixed_LOTW_ML .or. CS%apply_LOTW_floor) then

    ! Find the friction velocity and the absolute value of the Coriolis parameter at this point.
    u_star(:,:) = 0.0  ! Zero out the friction velocity on land points.
    tau_mag(:,:) = 0.0  ! Zero out the friction velocity on land points.

    if (allocated(tv%SpV_avg)) then
      rho_av1(:,:) = 0.0

      if (work_on_u) then
        do j=js,je ; do I=is,ie ; if (do_i(i,j)) then
          u_star(I,j) = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i+1,j))
          rho_av1(I,j) = 2. / (tv%SpV_avg(i,j,1) + tv%SpV_avg(i+1,j,1))
          absf(I,j) = 0.5 * (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))
        endif ; enddo ; enddo

        if (do_OBCs) then
          if (OBC%u_E_OBCs_on_PE) then
            do j=js_E_OBC,je_E_OBC ; do I=Is_E_OBC,Ie_E_OBC
              if (do_i(I,j) .and. OBC%segnum_u(I,j) > 0) then
                u_star(I,j) = Ustar_2d(i,j)
                rho_av1(I,j) = 1. / tv%SpV_avg(i,j,1)
              endif
            enddo ; enddo
          endif

          if (OBC%u_W_OBCs_on_PE) then
            do j=js_W_OBC,je_W_OBC ; do I=Is_W_OBC,Ie_W_OBC
              if (do_i(I,j) .and. OBC%segnum_u(I,j) < 0) then
                u_star(I,j) = Ustar_2d(i+1,j)
                rho_av1(I,j) = 1. / tv%SpV_avg(i+1,j,1)
              endif
            enddo ; enddo
          endif
        endif
      else
        do J=Js,Je ; do i=is,ie ; if (do_i(i,J)) then
          u_star(i,J) = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i,j+1))
          rho_av1(i,J) = 2. / (tv%SpV_avg(i,j,1) + tv%SpV_avg(i,j+1,1))
          absf(i,J) = 0.5 * (abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
        endif ; enddo ; enddo

        if (do_OBCs) then
          if (OBC%v_N_OBCs_on_PE) then
            do J=Js_N_OBC,Je_N_OBC ; do i=is_N_OBC,ie_N_obc
              if (do_i(i,J) .and. OBC%segnum_v(i,J) > 0) then
                u_star(i,J) = Ustar_2d(i,j)
                rho_av1(i,J) = 1. / tv%SpV_avg(i,j,1)
              endif
            enddo ; enddo
          endif

          if (OBC%v_S_OBCs_on_PE) then
            do J=Js_S_OBC,Je_S_OBC ; do i=is_S_OBC,ie_S_obc
              if (do_i(i,J) .and. OBC%segnum_v(i,J) < 0) then
                u_star(i,J) = Ustar_2d(i,j+1)
                rho_av1(i,J) = 1. / tv%SpV_avg(i,j+1,1)
              endif
            enddo ; enddo
          endif
        endif
      endif

      do J=Js,Je ; do I=is,ie
        tau_mag(I,J) = GV%RZ_to_H * rho_av1(i,j) * u_star(I,J)**2
      enddo ; enddo
    else ! (.not.allocated(tv%SpV_avg))
      if (work_on_u) then
        do j=js,je ; do I=is,ie ; if (do_i(I,j)) then
          u_star(I,j) = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i+1,j))
          absf(I,j) = 0.5 * (abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))
        endif ; enddo ; enddo

        if (do_OBCs) then
          if (OBC%u_E_OBCs_on_PE) then
            do j=js_E_OBC,je_E_OBC ; do I=Is_E_OBC,Ie_E_OBC
              if (do_i(I,j) .and. OBC%segnum_u(I,j) > 0) then
                u_star(I,j) = Ustar_2d(i,j)
              endif
            enddo ; enddo
          endif

          if (OBC%u_W_OBCs_on_PE) then
            do j=js_W_OBC,je_W_OBC ; do I=Is_W_OBC,Ie_W_OBC
              if (do_i(I,j) .and. OBC%segnum_u(I,j) < 0) then
                u_star(I,j) = Ustar_2d(i+1,j)
              endif
            enddo ; enddo
          endif
        endif
      else
        do J=Js,Je ; do i=is,ie ; if (do_i(i,J)) then
          u_star(i,J) = 0.5 * (Ustar_2d(i,j) + Ustar_2d(i,j+1))
          absf(i,J) = 0.5 * (abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
        endif ; enddo ; enddo

        if (do_OBCs) then
          if (OBC%v_N_OBCs_on_PE) then
            do J=Js_N_OBC,Je_N_OBC ; do i=is_N_OBC,ie_N_OBC
              if (do_i(i,J) .and. OBC%segnum_v(i,J) > 0) then
                u_star(i,J) = Ustar_2d(i,j)
              endif
            enddo ; enddo
          endif

          if (OBC%v_S_OBCs_on_PE) then
            do J=Js_S_OBC,Je_S_OBC ; do i=is_S_OBC,ie_S_OBC
              if (do_i(i,J) .and. OBC%segnum_v(i,J) < 0) then
                u_star(i,J) = Ustar_2d(i,j+1)
              endif
            enddo ; enddo
          endif
        endif
      endif
      do J=Js,Je ; do I=is,ie
        tau_mag(I,J) = GV%Z_to_H*u_star(I,J)**2
      enddo ; enddo
    endif

    ! Determine the thickness of the surface ocean boundary layer and its extent in index space.
    nk_in_ml(:,:) = 0
    if (CS%dynamic_viscous_ML) then
      ! The fractional number of layers that are within the viscous boundary layer were
      ! previously stored in visc%nkml_visc_[uv].
      h_ml(:,:) = h_neglect
      max_nk = 0
      if (work_on_u) then
        do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
          nk_in_ml(I,j) = ceiling(visc%nkml_visc_u(I,j))
          max_nk = max(max_nk, nk_in_ml(I,j))
        endif ; enddo ; enddo

        do k=1,max_nk
          do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
            if (k <= visc%nkml_visc_u(I,j)) then ! This layer is all in the ML.
              h_ml(i,j) = h_ml(i,j) + hvel(i,j,k)
            elseif (k < visc%nkml_visc_u(I,j) + 1.) then ! Part of this layer is in the ML.
              h_ml(i,j) = h_ml(i,j) + ((visc%nkml_visc_u(I,j) + 1.) - k) * hvel(i,j,k)
            endif
          endif ; enddo ; enddo
        enddo
      else
        do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
          nk_in_ml(i,j) = ceiling(visc%nkml_visc_v(i,J))
          max_nk = max(max_nk, nk_in_ml(i,j))
        endif ; enddo ; enddo

        do k=1,max_nk
          do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
            if (k <= visc%nkml_visc_v(i,J)) then ! This layer is all in the ML.
              h_ml(i,j) = h_ml(i,j) + hvel(i,j,k)
            elseif (k < visc%nkml_visc_v(i,J) + 1.) then ! Part of this layer is in the ML.
              h_ml(i,j) = h_ml(i,j) + ((visc%nkml_visc_v(i,J) + 1.) - k) * hvel(i,j,k)
            endif
          endif ; enddo ; enddo
        enddo
      endif
    elseif (GV%nkml>0) then
      ! This is a simple application of a refined-bulk mixed layer with GV%nkml sublayers.
      max_nk = GV%nkml
      do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
        nk_in_ml(i,j) = GV%nkml
      endif ; enddo ; enddo

      h_ml(:,:) = h_neglect

      do k=1,GV%nkml
        do j=js,je ; do i=is,ie ; if (do_i(i,j)) then
          h_ml(i,j) = h_ml(i,j) + hvel(i,j,k)
        endif ; enddo ; enddo
      enddo
    elseif (CS%fixed_LOTW_ML .or. CS%apply_LOTW_floor) then
      ! Determine which interfaces are within CS%Hmix of the surface, and set the viscous
      ! boundary layer thickness to the the smaller of CS%Hmix and the depth of the ocean.
      h_ml(:,:) = 0.0
      do k=1,nz
        can_exit = .true.
        do j=js,je ; do i=is,ie ; if (do_i(i,j) .and. (h_ml(i,j) < CS%Hmix)) then
          nk_in_ml(i,j) = k

          if (h_ml(i,j) + hvel(i,j,k) < CS%Hmix) then
            h_ml(i,j) = h_ml(i,j) + hvel(i,j,k)
            can_exit = .false.  ! Part of the next deeper layer is also in the mixed layer.
          else
            h_ml(i,j) = CS%Hmix
          endif
        endif ; enddo ; enddo

        if (can_exit) exit  ! All remaining layers in this row are below the mixed layer depth.
      enddo

      max_nk = 0
      do j=js,je ; do i=is,ie
        max_nk = max(max_nk, nk_in_ml(i,j))
      enddo ; enddo
    endif

    ! Avoid working on land or on columns where the viscous coupling could not be increased.
    do j=js,je ; do i=is,ie ; if ((u_star(i,j)<=0.0) .or. (.not.do_i(i,j))) then
      nk_in_ml(i,j) = 0
    endif ; enddo ; enddo

    ! Set the viscous coupling at the interfaces as the larger of what was previously
    ! set and the contributions from the surface boundary layer.
    z_t(:,:) = 0.0
    if (CS%apply_LOTW_floor .and. &
        (CS%dynamic_viscous_ML .or. (GV%nkml>0) .or. CS%fixed_LOTW_ML)) then
      do K=2,max_nk
        do j=js,je ; do i=is,ie ; if (k <= nk_in_ml(i,j)) then
          z_t(i,j) = z_t(i,j) + hvel(i,j,k-1)

          !   The viscosity in visc_ml is set to go to 0 at the mixed layer top and bottom
          ! (in a log-layer) and be further limited by rotation to give the natural Ekman length.
          temp1 = (z_t(i,j)*h_ml(i,j) - z_t(i,j)*z_t(i,j))
          if (GV%Boussinesq) then
            ustar2_denom = (CS%vonKar * GV%Z_to_H * u_star(i,j)**2) &
                / (absf(i,j) * temp1 + (h_ml(i,j) + h_neglect) * u_star(i,j))
          else
            ustar2_denom = (CS%vonKar * tau_mag(i,j)) &
                / (absf(i,j) * temp1 + (h_ml(i,j) + h_neglect) * u_star(i,j))
          endif

          visc_ml = temp1 * ustar2_denom
          ! Set the viscous coupling based on the model's vertical resolution.  The omission of
          ! the I_amax factor here is consistent with answer dates above 20190101.
          a_ml = visc_ml / (0.25 * (hvel(i,j,k) + hvel(i,j,k-1) + h_neglect))

          ! As a floor on the viscous coupling, assume that the length scale in the denominator can
          ! not be larger than the distance from the surface, consistent with a logarithmic velocity
          ! profile.  This is consistent with visc_ml, but cancels out common factors of z_t.
          a_floor = (h_ml(i,j) - z_t(i,j)) * ustar2_denom

          ! Choose the largest estimate of a_cpl.
          a_cpl(i,j,K) = max(a_cpl(i,j,K), a_ml, a_floor)
          ! An option could be added to change this to: a_cpl(i,K) = max(a_cpl(i,K) + a_ml, a_floor)
        endif ; enddo ; enddo
      enddo
    elseif (CS%apply_LOTW_floor) then
      do K=2,max_nk
        do j=js,je ; do i=is,ie ; if (k <= nk_in_ml(i,j)) then
          z_t(i,j) = z_t(i,j) + hvel(i,j,k-1)

          temp1 = (z_t(i,j)*h_ml(i,j) - z_t(i,j) * z_t(i,j))
          if (GV%Boussinesq) then
            ustar2_denom = (CS%vonKar * GV%Z_to_H * u_star(i,j)**2) &
                / (absf(i,j) * temp1 + (h_ml(i,j) + h_neglect) * u_star(i,j))
          else
            ustar2_denom = (CS%vonKar * tau_mag(i,j)) &
                / (absf(i,j) * temp1 + (h_ml(i,j) + h_neglect) * u_star(i,j))
          endif

          ! As a floor on the viscous coupling, assume that the length scale in the denominator can not
          ! be larger than the distance from the surface, consistent with a logarithmic velocity profile.
          a_cpl(i,j,K) = max(a_cpl(i,j,K), (h_ml(i,j) - z_t(i,j)) * ustar2_denom)
        endif ; enddo ; enddo
      enddo
    else
      do K=2,max_nk
        do j=js,je ; do i=is,ie ; if (k <= nk_in_ml(i,j)) then
          z_t(i,j) = z_t(i,j) + hvel(i,j,k-1)

          temp1 = (z_t(i,j) * h_ml(i,j) - z_t(i,j) * z_t(i,j))
          !   This viscosity is set to go to 0 at the mixed layer top and bottom (in a log-layer)
          ! and be further limited by rotation to give the natural Ekman length.
          ! The following expressions are mathematically equivalent.
          if (GV%Boussinesq .or. (CS%answer_date < 20230601)) then
            visc_ml = u_star(i,j) * CS%vonKar * (GV%Z_to_H * temp1 * u_star(i,j)) &
                / (absf(i,j) * temp1 + (h_ml(i,j)+h_neglect) * u_star(i,j))
          else
            visc_ml = CS%vonKar * (temp1 * tau_mag(i,j)) &
                / (absf(i,j) * temp1 + (h_ml(i,j) + h_neglect) * u_star(i,j))
          endif
          a_ml = visc_ml / (0.25 * (hvel(i,j,k) + hvel(i,j,k-1) + h_neglect) + 0.5 * I_amax * visc_ml)

          ! Choose the largest estimate of a_cpl, but these could be changed to be additive.
          a_cpl(i,j,K) = max(a_cpl(i,j,K), a_ml)
          ! An option could be added to change this to: a_cpl(i,K) = a_cpl(i,K) + a_ml
        endif ; enddo ; enddo
      enddo
    endif
  endif
end subroutine find_coupling_coef

!> Velocity components which exceed a threshold for physically reasonable values
!! are truncated. Optionally, any column with excessive velocities may be sent
!! to a diagnostic reporting subroutine.
subroutine vertvisc_limit_vel(u, v, h, ADp, CDp, forces, visc, dt, G, GV, US, CS)
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u      !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v      !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  type(accel_diag_ptrs),   intent(in)    :: ADp    !< Acceleration diagnostic pointers
  type(cont_diag_ptrs),    intent(in)    :: CDp    !< Continuity diagnostic pointers
  type(mech_forcing),      intent(in)    :: forces !< A structure with the driving mechanical forces
  type(vertvisc_type),     intent(in)    :: visc   !< Viscosities and bottom drag
  real,                    intent(in)    :: dt     !< Time increment [T ~> s]
  type(vertvisc_CS),       pointer       :: CS     !< Vertical viscosity control structure

  ! Local variables
  real :: CFL              ! The local CFL number [nondim]
  real :: H_report         ! A thickness below which not to report truncations [H ~> m or kg m-2]
  real :: vel_report(SZIB_(G),SZJB_(G))   ! The velocity to report [L T-1 ~> m s-1]
  real :: u_old(SZIB_(G),SZJ_(G),SZK_(GV)) ! The previous u-velocity [L T-1 ~> m s-1]
  real :: v_old(SZI_(G),SZJB_(G),SZK_(GV)) ! The previous v-velocity [L T-1 ~> m s-1]
  logical :: trunc_any, dowrite(SZIB_(G),SZJB_(G))
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  H_report = 6.0 * GV%Angstrom_H

  if (len_trim(CS%u_trunc_file) > 0) then
    !$OMP parallel do default(shared) private(trunc_any,CFL)
    do j=js,je
      trunc_any = .false.
      do I=Isq,Ieq ; dowrite(I,j) = .false. ; enddo
      do I=Isq,Ieq ; vel_report(i,j) = 3.0e8*US%m_s_to_L_T ; enddo ! Speed of light default.
      do k=1,nz ; do I=Isq,Ieq
        if (abs(u(I,j,k)) < CS%vel_underflow) u(I,j,k) = 0.0
        if (u(I,j,k) < 0.0) then
          CFL = (-u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
        else
          CFL = (u(I,j,k) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
        endif
        if (CFL > CS%CFL_trunc) trunc_any = .true.
        if (CFL > CS%CFL_report) then
          dowrite(I,j) = .true.
          vel_report(I,j) = MIN(vel_report(I,j), abs(u(I,j,k)))
        endif
      enddo ; enddo

      do I=Isq,Ieq ; if (dowrite(I,j)) then
        u_old(I,j,:) = u(I,j,:)
      endif ; enddo

      if (trunc_any) then
        do k=1,nz ; do I=Isq,Ieq
          if ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i+1,j) < -CS%CFL_trunc) then
            u(I,j,k) = (-0.9*CS%CFL_trunc) * (G%areaT(i+1,j) / (dt * G%dy_Cu(I,j)))
            if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
          elseif ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i,j) > CS%CFL_trunc) then
            u(I,j,k) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dy_Cu(I,j)))
            if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
          endif
        enddo ; enddo
      endif
    enddo ! j-loop
  else  ! Do not report accelerations leading to large velocities.
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=js,je ; do I=Isq,Ieq
      if (abs(u(I,j,k)) < CS%vel_underflow) then ; u(I,j,k) = 0.0
      elseif ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i+1,j) < -CS%CFL_trunc) then
        u(I,j,k) = (-0.9*CS%CFL_trunc) * (G%areaT(i+1,j) / (dt * G%dy_Cu(I,j)))
        if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
      elseif ((u(I,j,k) * (dt * G%dy_Cu(I,j))) * G%IareaT(i,j) > CS%CFL_trunc) then
        u(I,j,k) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dy_Cu(I,j)))
        if (h(i,j,k) + h(i+1,j,k) > H_report) CS%ntrunc = CS%ntrunc + 1
      endif
    enddo ; enddo ; enddo
  endif

  if (len_trim(CS%u_trunc_file) > 0) then
    do j=js,je ; do I=Isq,Ieq ; if (dowrite(I,j)) then
      ! Call a diagnostic reporting subroutines are called if unphysically large values are found.
      call write_u_accel(I, j, u_old, h, ADp, CDp, dt, G, GV, US, CS%PointAccel_CSp, &
                         vel_report(I,j), forces%taux(I,j), a=CS%a_u, hv=CS%h_u)
    endif ; enddo ; enddo
  endif

  if (len_trim(CS%v_trunc_file) > 0) then
    !$OMP parallel do default(shared) private(trunc_any,CFL)
    do J=Jsq,Jeq
      trunc_any = .false.
      do i=is,ie ; dowrite(i,J) = .false. ; enddo
      do i=is,ie ; vel_report(i,J) = 3.0e8*US%m_s_to_L_T ; enddo ! Speed of light default.
      do k=1,nz ; do i=is,ie
        if (abs(v(i,J,k)) < CS%vel_underflow) v(i,J,k) = 0.0
        if (v(i,J,k) < 0.0) then
          CFL = (-v(i,J,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
        else
          CFL = (v(i,J,k) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
        endif
        if (CFL > CS%CFL_trunc) trunc_any = .true.
        if (CFL > CS%CFL_report) then
          dowrite(i,J) = .true.
          vel_report(i,J) = MIN(vel_report(i,J), abs(v(i,J,k)))
        endif
      enddo ; enddo

      do i=is,ie ; if (dowrite(i,J)) then
        v_old(i,J,:) = v(i,J,:)
      endif ; enddo

      if (trunc_any) then
        do k=1,nz ; do i=is,ie
          if ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j+1) < -CS%CFL_trunc) then
            v(i,J,k) = (-0.9*CS%CFL_trunc) * (G%areaT(i,j+1) / (dt * G%dx_Cv(i,J)))
            if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
          elseif ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j) > CS%CFL_trunc) then
            v(i,J,k) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dx_Cv(i,J)))
            if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
          endif
        enddo ; enddo
      endif
    enddo ! J-loop
  else  ! Do not report accelerations leading to large velocities.
    !$OMP parallel do default(shared)
    do k=1,nz ; do J=Jsq,Jeq ; do i=is,ie
      if (abs(v(i,J,k)) < CS%vel_underflow) then ; v(i,J,k) = 0.0
      elseif ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j+1) < -CS%CFL_trunc) then
        v(i,J,k) = (-0.9*CS%CFL_trunc) * (G%areaT(i,j+1) / (dt * G%dx_Cv(i,J)))
        if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
      elseif ((v(i,J,k) * (dt * G%dx_Cv(i,J))) * G%IareaT(i,j) > CS%CFL_trunc) then
        v(i,J,k) = (0.9*CS%CFL_trunc) * (G%areaT(i,j) / (dt * G%dx_Cv(i,J)))
        if (h(i,j,k) + h(i,j+1,k) > H_report) CS%ntrunc = CS%ntrunc + 1
      endif
    enddo ; enddo ; enddo
  endif

  if (len_trim(CS%v_trunc_file) > 0) then
    do J=Jsq,Jeq ; do i=is,ie ; if (dowrite(i,J)) then
      ! Call a diagnostic reporting subroutines are called if unphysically large values are found.
      call write_v_accel(i, J, v_old, h, ADp, CDp, dt, G, GV, US, CS%PointAccel_CSp, &
                         vel_report(i,J), forces%tauy(i,J), a=CS%a_v, hv=CS%h_v)
    endif ; enddo ; enddo
  endif

end subroutine vertvisc_limit_vel

!> Initialize the vertical friction module
subroutine vertvisc_init(MIS, Time, G, GV, US, param_file, diag, ADp, dirs, &
                          ntrunc, CS, fpmix)
  type(ocean_internal_state), &
                   target, intent(in)    :: MIS    !< The "MOM Internal State", a set of pointers
                                                   !! to the fields and accelerations that make
                                                   !! up the ocean's physical state
  type(time_type), target, intent(in)    :: Time   !< Current model time
  type(ocean_grid_type),   intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< File to parse for parameters
  type(diag_ctrl), target, intent(inout) :: diag   !< Diagnostic control structure
  type(accel_diag_ptrs),   intent(inout) :: ADp    !< Acceleration diagnostic pointers
  type(directories),       intent(in)    :: dirs   !< Relevant directory paths
  integer, target,         intent(inout) :: ntrunc !< Number of velocity truncations
  type(vertvisc_CS),       pointer       :: CS     !< Vertical viscosity control structure
  logical, optional,       intent(in)    :: fpmix  !< Nonlocal momentum mixing

  ! Local variables

  real :: Kv_BBL  ! A viscosity in the bottom boundary layer with a simple scheme [H Z T-1 ~> m2 s-1 or Pa s]
  real :: Kv_back_z  ! A background kinematic viscosity [Z2 T-1 ~> m2 s-1]
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  logical :: lfpmix
  character(len=200) :: kappa_gl90_file, inputdir, kdgl90_varname
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_vert_friction" ! This module's name.
  character(len=40)  :: thickness_units
  real :: Kv_mks ! KVML in MKS [m2 s-1]

  if (associated(CS)) then
    call MOM_error(WARNING, "vertvisc_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%initialized = .true.

  if (GV%Boussinesq) then; thickness_units = "m"
  else; thickness_units = "kg m-2"; endif

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%diag => diag ; CS%ntrunc => ntrunc ; ntrunc = 0

  lfpmix = .false.
  if (present(fpmix)) lfpmix = fpmix

! Default, read and log parameters
  call log_version(param_file, mdl, version, "", log_to_all=.true., debugging=.true.)
  call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
  call get_param(param_file, mdl, "VERT_FRICTION_ANSWER_DATE", CS%answer_date, &
                 "The vintage of the order of arithmetic and expressions in the viscous "//&
                 "calculations.  Values below 20190101 recover the answers from the end of 2018, "//&
                 "while higher values use expressions that do not use an arbitrary hard-coded "//&
                 "maximum viscous coupling coefficient between layers.  Values below 20230601 "//&
                 "recover a form of the viscosity within the mixed layer that breaks up the "//&
                 "magnitude of the wind stress in some non-Boussinesq cases.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
  if (.not.GV%Boussinesq) CS%answer_date = max(CS%answer_date, 20230701)

  call get_param(param_file, mdl, "BOTTOMDRAGLAW", CS%bottomdraglaw, &
                 "If true, the bottom stress is calculated with a drag "//&
                 "law of the form c_drag*|u|*u. The velocity magnitude "//&
                 "may be an assumed value or it may be based on the "//&
                 "actual velocity in the bottommost HBBL, depending on "//&
                 "LINEAR_DRAG.", default=.true.)
  call get_param(param_file, mdl, "DIRECT_STRESS", CS%direct_stress, &
                 "If true, the wind stress is distributed over the topmost HMIX_STRESS of fluid "//&
                 "(like in HYCOM), and an added mixed layer viscosity or a physically based "//&
                 "boundary layer turbulence parameterization is not needed for stability.", &
                 default=.false.)
  call get_param(param_file, mdl, "DYNAMIC_VISCOUS_ML", CS%dynamic_viscous_ML, &
                 "If true, use a bulk Richardson number criterion to "//&
                 "determine the mixed layer thickness for viscosity.", &
                 default=.false.)
  call get_param(param_file, mdl, "FIXED_DEPTH_LOTW_ML", CS%fixed_LOTW_ML, &
                 "If true, use a Law-of-the-wall prescription for the mixed layer viscosity "//&
                 "within a boundary layer that is the lesser of HMIX_FIXED and the total "//&
                 "depth of the ocean in a column.", default=.false.)
  call get_param(param_file, mdl, "LOTW_VISCOUS_ML_FLOOR", CS%apply_LOTW_floor, &
                 "If true, use a Law-of-the-wall prescription to set a lower bound on the "//&
                 "viscous coupling between layers within the surface boundary layer, based "//&
                 "the distance of interfaces from the surface.  This only acts when there "//&
                 "are large changes in the thicknesses of successive layers or when the "//&
                 "viscosity is set externally and the wind stress has subsequently increased.", &
                 default=.false.)
  call get_param(param_file, mdl, 'VON_KARMAN_CONST', CS%vonKar, &
                 'The value the von Karman constant as used for mixed layer viscosity.', &
                 units='nondim', default=0.41)
  call get_param(param_file, mdl, "U_TRUNC_FILE", CS%u_trunc_file, &
                 "The absolute path to a file into which the accelerations "//&
                 "leading to zonal velocity truncations are written. "//&
                 "Undefine this for efficiency if this diagnostic is not needed.", &
                 default=" ", debuggingParam=.true.)
  call get_param(param_file, mdl, "V_TRUNC_FILE", CS%v_trunc_file, &
                 "The absolute path to a file into which the accelerations "//&
                 "leading to meridional velocity truncations are written. "//&
                 "Undefine this for efficiency if this diagnostic is not needed.", &
                 default=" ", debuggingParam=.true.)
  call get_param(param_file, mdl, "HARMONIC_VISC", CS%harmonic_visc, &
                 "If true, use the harmonic mean thicknesses for "//&
                 "calculating the vertical viscosity.", default=.false.)
  call get_param(param_file, mdl, "HARMONIC_BL_SCALE", CS%harm_BL_val, &
                 "A scale to determine when water is in the boundary "//&
                 "layers based solely on harmonic mean thicknesses for "//&
                 "the purpose of determining the extent to which the "//&
                 "thicknesses used in the viscosities are upwinded.", &
                 default=0.0, units="nondim")
  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)

  if (GV%nkml < 1) then
    call get_param(param_file, mdl, "HMIX_FIXED", CS%Hmix, &
                 "The prescribed depth over which the near-surface viscosity and "//&
                 "diffusivity are elevated when the bulk mixed layer is not used.", &
                 units="m", scale=US%m_to_Z, fail_if_missing=.true.)
  endif
  if (CS%direct_stress) then
    if (GV%nkml < 1) then
      call get_param(param_file, mdl, "HMIX_STRESS", CS%Hmix_stress, &
                 "The depth over which the wind stress is applied if DIRECT_STRESS is true.", &
                 units="m", default=US%Z_to_m*CS%Hmix, scale=GV%m_to_H)
    else
      call get_param(param_file, mdl, "HMIX_STRESS", CS%Hmix_stress, &
                 "The depth over which the wind stress is applied if DIRECT_STRESS is true.", &
                 units="m", fail_if_missing=.true., scale=GV%m_to_H)
    endif
    if (CS%Hmix_stress <= 0.0) call MOM_error(FATAL, "vertvisc_init: " // &
       "HMIX_STRESS must be set to a positive value if DIRECT_STRESS is true.")
  endif
  call get_param(param_file, mdl, "KV", Kv_back_z, &
                 "The background kinematic viscosity in the interior. "//&
                 "The molecular value, ~1e-6 m2 s-1, may be used.", &
                 units="m2 s-1", fail_if_missing=.true., scale=US%m2_s_to_Z2_T)
  ! Convert input kinematic viscosity to dynamic viscosity when non-Boussinesq.
  CS%Kv = (US%Z2_T_to_m2_s*GV%m2_s_to_HZ_T) * Kv_back_z

  call get_param(param_file, mdl, "USE_GL90_IN_SSW", CS%use_GL90_in_SSW, &
                 "If true, use simpler method to calculate 1/N^2 in GL90 vertical "// &
                 "viscosity coefficient. This method is valid in stacked shallow water mode.", &
                 default=.false.)
  call get_param(param_file, mdl, "KD_GL90", CS%kappa_gl90, &
                 "The scalar diffusivity used in GL90 vertical viscosity scheme.", &
                 units="m2 s-1", default=0.0, scale=US%m_to_L*US%Z_to_L*GV%m_to_H*US%T_to_s, &
                 do_not_log=.not.CS%use_GL90_in_SSW)
  call get_param(param_file, mdl, "READ_KD_GL90", CS%read_kappa_gl90, &
                 "If true, read a file (given by KD_GL90_FILE) containing the "//&
                 "spatially varying diffusivity KD_GL90 used in the GL90 scheme.", default=.false., &
                 do_not_log=.not.CS%use_GL90_in_SSW)
  if (CS%read_kappa_gl90) then
    if (CS%kappa_gl90 > 0) then
        call MOM_error(FATAL, "MOM_vert_friction.F90, vertvisc_init: KD_GL90 > 0 "// &
              "is not compatible with READ_KD_GL90 = .TRUE. ")
    endif
    call get_param(param_file, mdl, "INPUTDIR", inputdir, &
                 "The directory in which all input files are found.", &
                 default=".", do_not_log=.true.)
    inputdir = slasher(inputdir)
    call get_param(param_file, mdl, "KD_GL90_FILE", kappa_gl90_file, &
                 "The file containing the spatially varying diffusivity used in the "// &
                 "GL90 scheme.", default="kd_gl90.nc", do_not_log=.not.CS%use_GL90_in_SSW)
    call get_param(param_file, mdl, "KD_GL90_VARIABLE", kdgl90_varname, &
                 "The name of the GL90 diffusivity variable to read "//&
                 "from KD_GL90_FILE.", default="kd_gl90", do_not_log=.not.CS%use_GL90_in_SSW)
    kappa_gl90_file = trim(inputdir) // trim(kappa_gl90_file)

    allocate(CS%kappa_gl90_2d(G%isd:G%ied, G%jsd:G%jed), source=0.0)
    call MOM_read_data(kappa_gl90_file, kdgl90_varname, CS%kappa_gl90_2d(:,:), G%domain, &
                       scale=US%m_to_L*US%Z_to_L*GV%m_to_H*US%T_to_s)
    call pass_var(CS%kappa_gl90_2d, G%domain)
  endif
  call get_param(param_file, mdl, "USE_GL90_N2", CS%use_GL90_N2, &
                 "If true, use GL90 vertical viscosity coefficient that is depth-independent; "// &
                 "this corresponds to a kappa_GM that scales as N^2 with depth.", &
                 default=.false., do_not_log=.not.CS%use_GL90_in_SSW)
  if (CS%use_GL90_N2) then
    if (.not. CS%use_GL90_in_SSW) call MOM_error(FATAL, &
           "MOM_vert_friction.F90, vertvisc_init: "//&
           "When USE_GL90_N2=True, USE_GL90_in_SSW must also be True.")
    if (CS%kappa_gl90 > 0) then
        call MOM_error(FATAL, "MOM_vert_friction.F90, vertvisc_init: KD_GL90 > 0 "// &
              "is not compatible with USE_GL90_N2 = .TRUE. ")
    endif
    if (CS%read_kappa_gl90) call MOM_error(FATAL, &
           "MOM_vert_friction.F90, vertvisc_init: "//&
           "READ_KD_GL90 = .TRUE. is not compatible with USE_GL90_N2 = .TRUE.")
    call get_param(param_file, mdl, "alpha_GL90", CS%alpha_gl90, &
                   "Coefficient used to compute a depth-independent GL90 vertical "//&
                   "viscosity via Kv_GL90 = alpha_GL90 * f2. Is only used "// &
                   "if USE_GL90_N2 is true. Note that the implied Kv_GL90 "// &
                   "corresponds to a KD_GL90 that scales as N^2 with depth.", &
                   units="m2 s", default=0.0, scale=GV%m_to_H*US%m_to_Z*US%s_to_T, &
                   do_not_log=.not.CS%use_GL90_in_SSW)
  endif
  call get_param(param_file, mdl, "HBBL_GL90", CS%Hbbl_gl90, &
                 "The thickness of the GL90 bottom boundary layer, "//&
                 "which defines the range over which the GL90 coupling "//&
                 "coefficient is zeroed out, in order to avoid fluxing "//&
                 "momentum into vanished layers over steep topography.", &
                 units="m", default=5.0, scale=US%m_to_Z, do_not_log=.not.CS%use_GL90_in_SSW)

  CS%Kvml_invZ2 = 0.0
  if (GV%nkml < 1) then
    call get_param(param_file, mdl, "KVML", Kv_mks, &
                 "The scale for an extra kinematic viscosity in the mixed layer", &
                 units="m2 s-1", default=-1.0, do_not_log=.true.)
    if (Kv_mks >= 0.0) then
      call MOM_error(WARNING, "KVML is a deprecated parameter. Use KV_ML_INVZ2 instead.")
    else
      Kv_mks = 0.0
    endif
    call get_param(param_file, mdl, "KV_ML_INVZ2", CS%Kvml_invZ2, &
                 "An extra kinematic viscosity in a mixed layer of thickness HMIX_FIXED, "//&
                 "with the actual viscosity scaling as 1/(z*HMIX_FIXED)^2, where z is the "//&
                 "distance from the surface, to allow for finite wind stresses to be "//&
                 "transmitted through infinitesimally thin surface layers.  This is an "//&
                 "older option for numerical convenience without a strong physical basis, "//&
                 "and its use is now discouraged.", &
                 units="m2 s-1", default=Kv_mks, scale=GV%m2_s_to_HZ_T)
  endif

  if (.not.CS%bottomdraglaw) then
    call get_param(param_file, mdl, "KV_EXTRA_BBL", CS%Kv_extra_bbl, &
                 "An extra kinematic viscosity in the benthic boundary layer. "//&
                 "KV_EXTRA_BBL is not used if BOTTOMDRAGLAW is true.", &
                 units="m2 s-1", default=0.0, scale=GV%m2_s_to_HZ_T, do_not_log=.true.)
    if (CS%Kv_extra_bbl == 0.0) then
      call get_param(param_file, mdl, "KVBBL", Kv_BBL, &
                 "An extra kinematic viscosity in the benthic boundary layer. "//&
                 "KV_EXTRA_BBL is not used if BOTTOMDRAGLAW is true.", &
                 units="m2 s-1", default=US%Z2_T_to_m2_s*Kv_back_z, scale=GV%m2_s_to_HZ_T, &
                 do_not_log=.true.)
      if (abs(Kv_BBL - CS%Kv) > 1.0e-15*abs(CS%Kv)) then
        call MOM_error(WARNING, "KVBBL is a deprecated parameter. Use KV_EXTRA_BBL instead.")
        CS%Kv_extra_bbl = Kv_BBL - CS%Kv
      endif
    endif
    call log_param(param_file, mdl, "KV_EXTRA_BBL", CS%Kv_extra_bbl, &
                 "An extra kinematic viscosity in the benthic boundary layer. "//&
                 "KV_EXTRA_BBL is not used if BOTTOMDRAGLAW is true.", &
                 units="m2 s-1", default=0.0, unscale=GV%HZ_T_to_m2_s)
  endif
  call get_param(param_file, mdl, "HBBL", CS%Hbbl, &
                 "The thickness of a bottom boundary layer with a viscosity increased by "//&
                 "KV_EXTRA_BBL if BOTTOMDRAGLAW is not defined, or the thickness over which "//&
                 "near-bottom velocities are averaged for the drag law if BOTTOMDRAGLAW is "//&
                 "defined but LINEAR_DRAG is not.", &
                 units="m", fail_if_missing=.true., scale=US%m_to_Z)
  call get_param(param_file, mdl, "CFL_TRUNCATE", CS%CFL_trunc, &
                 "The value of the CFL number that will cause velocity "//&
                 "components to be truncated; instability can occur past 0.5.", &
                 units="nondim", default=0.5)
  call get_param(param_file, mdl, "CFL_REPORT", CS%CFL_report, &
                 "The value of the CFL number that causes accelerations "//&
                 "to be reported; the default is CFL_TRUNCATE.", &
                 units="nondim", default=CS%CFL_trunc)
  call get_param(param_file, mdl, "CFL_TRUNCATE_RAMP_TIME", CS%truncRampTime, &
                 "The time over which the CFL truncation value is ramped "//&
                 "up at the beginning of the run.", &
                 units="s", default=0., scale=US%s_to_T)
  CS%CFL_truncE = CS%CFL_trunc
  call get_param(param_file, mdl, "CFL_TRUNCATE_START", CS%CFL_truncS, &
                 "The start value of the truncation CFL number used when "//&
                 "ramping up CFL_TRUNC.", &
                 units="nondim", default=0.)
  call get_param(param_file, mdl, "STOKES_MIXING_COMBINED", CS%StokesMixing, &
                 "Flag to use Stokes drift Mixing via the Lagrangian "//&
                 " current (Eulerian plus Stokes drift). "//&
                 " Still needs work and testing, so not recommended for use.",&
                 default=.false.)
  !BGR 04/04/2018{
  ! StokesMixing is required for MOM6 for some Langmuir mixing parameterization.
  !   The code used here has not been developed for vanishing layers or in
  !   conjunction with any bottom friction.  Therefore, the following line is
  !   added so this functionality cannot be used without user intervention in
  !   the code.  This will prevent general use of this functionality until proper
  !   care is given to the previously mentioned issues.  Comment out the following
  !   MOM_error to use, but do so at your own risk and with these points in mind.
  !}
  if (CS%StokesMixing) then
    call MOM_error(FATAL, "Stokes mixing requires user intervention in the code.\n"//&
                          "  Model now exiting.  See MOM_vert_friction.F90 for \n"//&
                          "  details (search 'BGR 04/04/2018' to locate comment).")
  endif
  call get_param(param_file, mdl, "VEL_UNDERFLOW", CS%vel_underflow, &
                 "A negligibly small velocity magnitude below which velocity "//&
                 "components are set to 0.  A reasonable value might be "//&
                 "1e-30 m/s, which is less than an Angstrom divided by "//&
                 "the age of the universe.", units="m s-1", default=0.0, scale=US%m_s_to_L_T)

  ALLOC_(CS%a_u(IsdB:IedB,jsd:jed,nz+1)) ; CS%a_u(:,:,:) = 0.0
  ALLOC_(CS%a_u_gl90(IsdB:IedB,jsd:jed,nz+1)) ; CS%a_u_gl90(:,:,:) = 0.0
  ALLOC_(CS%h_u(IsdB:IedB,jsd:jed,nz))   ; CS%h_u(:,:,:) = 0.0
  ALLOC_(CS%a_v(isd:ied,JsdB:JedB,nz+1)) ; CS%a_v(:,:,:) = 0.0
  ALLOC_(CS%a_v_gl90(isd:ied,JsdB:JedB,nz+1)) ; CS%a_v_gl90(:,:,:) = 0.0
  ALLOC_(CS%h_v(isd:ied,JsdB:JedB,nz))   ; CS%h_v(:,:,:) = 0.0

  CS%id_Kv_slow = register_diag_field('ocean_model', 'Kv_slow', diag%axesTi, Time, &
      'Slow varying vertical viscosity', 'm2 s-1', conversion=GV%HZ_T_to_m2_s)

  CS%id_Kv_u = register_diag_field('ocean_model', 'Kv_u', diag%axesCuL, Time, &
      'Total vertical viscosity at u-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_Kv_v = register_diag_field('ocean_model', 'Kv_v', diag%axesCvL, Time, &
      'Total vertical viscosity at v-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_Kv_gl90_u = register_diag_field('ocean_model', 'Kv_gl90_u', diag%axesCuL, Time, &
      'GL90 vertical viscosity at u-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_Kv_gl90_v = register_diag_field('ocean_model', 'Kv_gl90_v', diag%axesCvL, Time, &
      'GL90 vertical viscosity at v-points', 'm2 s-1', conversion=GV%H_to_m**2*US%s_to_T)

  CS%id_au_vv = register_diag_field('ocean_model', 'au_visc', diag%axesCui, Time, &
      'Zonal Viscous Vertical Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_av_vv = register_diag_field('ocean_model', 'av_visc', diag%axesCvi, Time, &
      'Meridional Viscous Vertical Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_au_gl90_vv = register_diag_field('ocean_model', 'au_gl90_visc', diag%axesCui, Time, &
      'Zonal Viscous Vertical GL90 Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_av_gl90_vv = register_diag_field('ocean_model', 'av_gl90_visc', diag%axesCvi, Time, &
      'Meridional Viscous Vertical GL90 Coupling Coefficient', 'm s-1', conversion=GV%H_to_m*US%s_to_T)

  CS%id_h_u = register_diag_field('ocean_model', 'Hu_visc', diag%axesCuL, Time, &
      'Thickness at Zonal Velocity Points for Viscosity', &
      thickness_units, conversion=GV%H_to_MKS)
      ! Alternately, to always give this variable in 'm' use the following line instead:
      ! 'm', conversion=GV%H_to_m)

  CS%id_h_v = register_diag_field('ocean_model', 'Hv_visc', diag%axesCvL, Time, &
      'Thickness at Meridional Velocity Points for Viscosity', &
      thickness_units, conversion=GV%H_to_MKS)

  CS%id_hML_u = register_diag_field('ocean_model', 'HMLu_visc', diag%axesCu1, Time, &
      'Mixed Layer Thickness at Zonal Velocity Points for Viscosity', &
      thickness_units, conversion=US%Z_to_m)

  CS%id_hML_v = register_diag_field('ocean_model', 'HMLv_visc', diag%axesCv1, Time, &
      'Mixed Layer Thickness at Meridional Velocity Points for Viscosity', &
      thickness_units, conversion=US%Z_to_m)

 if (lfpmix) then
  CS%id_uE_h = register_diag_field('ocean_model', 'uE_h' , CS%diag%axesTL, &
      Time, 'x-zonal Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vE_h = register_diag_field('ocean_model', 'vE_h' , CS%diag%axesTL, &
      Time, 'y-merid Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_uInc_h = register_diag_field('ocean_model','uInc_h',CS%diag%axesTL, &
      Time, 'x-zonal Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vInc_h = register_diag_field('ocean_model','vInc_h',CS%diag%axesTL, &
      Time, 'x-zonal Eulerian' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_uStk = register_diag_field('ocean_model', 'uStk' , CS%diag%axesTL, &
      Time, 'x-FP du increment' , 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vStk = register_diag_field('ocean_model', 'vStk' , CS%diag%axesTL, &
      Time, 'y-FP dv increment' , 'm s-1', conversion=US%L_T_to_m_s)

  CS%id_FPtau2s = register_diag_field('ocean_model','Omega_tau2s',CS%diag%axesTi, &
      Time, 'Stress direction from shear','radians')
  CS%id_FPtau2w = register_diag_field('ocean_model','Omega_tau2w',CS%diag%axesTi, &
      Time, 'Stress direction from wind','radians')
  CS%id_uStk0 = register_diag_field('ocean_model', 'uStk0' , diag%axesT1, &
      Time, 'Zonal Surface Stokes', 'm s-1', conversion=US%L_T_to_m_s)
  CS%id_vStk0 = register_diag_field('ocean_model', 'vStk0' , diag%axesT1, &
      Time, 'Merid Surface Stokes', 'm s-1', conversion=US%L_T_to_m_s)
  endif

  CS%id_du_dt_visc = register_diag_field('ocean_model', 'du_dt_visc', diag%axesCuL, Time, &
      'Zonal Acceleration from Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_du_dt_visc > 0) call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  CS%id_dv_dt_visc = register_diag_field('ocean_model', 'dv_dt_visc', diag%axesCvL, Time, &
      'Meridional Acceleration from Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_dv_dt_visc > 0) call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  CS%id_GLwork = register_diag_field('ocean_model', 'GLwork', diag%axesTL, Time, &
      'Sign-definite Kinetic Energy Source from GL90 Vertical Viscosity', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_du_dt_visc_gl90 = register_diag_field('ocean_model', 'du_dt_visc_gl90', diag%axesCuL, Time, &
      'Zonal Acceleration from GL90 Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_du_dt_visc_gl90 > 0) .or. (CS%id_GLwork > 0)) then
    call safe_alloc_ptr(ADp%du_dt_visc_gl90,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  endif
  CS%id_dv_dt_visc_gl90 = register_diag_field('ocean_model', 'dv_dt_visc_gl90', diag%axesCvL, Time, &
      'Meridional Acceleration from GL90 Vertical Viscosity', 'm s-2', conversion=US%L_T2_to_m_s2)
  if ((CS%id_dv_dt_visc_gl90 > 0) .or. (CS%id_GLwork > 0)) then
    call safe_alloc_ptr(ADp%dv_dt_visc_gl90,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  endif
  CS%id_du_dt_str = register_diag_field('ocean_model', 'du_dt_str', diag%axesCuL, Time, &
      'Zonal Acceleration from Surface Wind Stresses', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_du_dt_str > 0) call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
  CS%id_dv_dt_str = register_diag_field('ocean_model', 'dv_dt_str', diag%axesCvL, Time, &
      'Meridional Acceleration from Surface Wind Stresses', 'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_dv_dt_str > 0) call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)

  CS%id_taux_bot = register_diag_field('ocean_model', 'taux_bot', diag%axesCu1, &
      Time, 'Zonal Bottom Stress from Ocean to Earth', &
      'Pa', conversion=US%RZ_to_kg_m2*US%L_T2_to_m_s2)
  CS%id_tauy_bot = register_diag_field('ocean_model', 'tauy_bot', diag%axesCv1, &
      Time, 'Meridional Bottom Stress from Ocean to Earth', &
      'Pa', conversion=US%RZ_to_kg_m2*US%L_T2_to_m_s2)

  !CS%id_hf_du_dt_visc = register_diag_field('ocean_model', 'hf_du_dt_visc', diag%axesCuL, Time, &
  !    'Fractional Thickness-weighted Zonal Acceleration from Vertical Viscosity', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_du_dt_visc > 0) then
  !  call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
  !  call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  !endif

  !CS%id_hf_dv_dt_visc = register_diag_field('ocean_model', 'hf_dv_dt_visc', diag%axesCvL, Time, &
  !    'Fractional Thickness-weighted Meridional Acceleration from Vertical Viscosity', &
  !    'm s-2', v_extensive=.true., conversion=US%L_T2_to_m_s2)
  !if (CS%id_hf_dv_dt_visc > 0) then
  !  call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  !  call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)
  !endif

  CS%id_hf_du_dt_visc_2d = register_diag_field('ocean_model', 'hf_du_dt_visc_2d', diag%axesCu1, Time, &
      'Depth-sum Fractional Thickness-weighted Zonal Acceleration from Vertical Viscosity', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_du_dt_visc_2d > 0) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%diag_hfrac_u,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_hf_dv_dt_visc_2d = register_diag_field('ocean_model', 'hf_dv_dt_visc_2d', diag%axesCv1, Time, &
      'Depth-sum Fractional Thickness-weighted Meridional Acceleration from Vertical Viscosity', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_hf_dv_dt_visc_2d > 0) then
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%diag_hfrac_v,isd,ied,JsdB,JedB,nz)
  endif

  CS%id_h_du_dt_visc = register_diag_field('ocean_model', 'h_du_dt_visc', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_du_dt_visc > 0) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%diag_hu,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_h_dv_dt_visc = register_diag_field('ocean_model', 'h_dv_dt_visc', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Acceleration from Horizontal Viscosity', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_dv_dt_visc > 0) then
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%diag_hv,isd,ied,JsdB,JedB,nz)
  endif

  CS%id_h_du_dt_str = register_diag_field('ocean_model', 'h_du_dt_str', diag%axesCuL, Time, &
      'Thickness Multiplied Zonal Acceleration from Surface Wind Stresses', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_du_dt_str > 0) then
    call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%diag_hu,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_h_dv_dt_str = register_diag_field('ocean_model', 'h_dv_dt_str', diag%axesCvL, Time, &
      'Thickness Multiplied Meridional Acceleration from Surface Wind Stresses', &
      'm2 s-2', conversion=GV%H_to_m*US%L_T2_to_m_s2)
  if (CS%id_h_dv_dt_str > 0) then
    call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%diag_hv,isd,ied,JsdB,JedB,nz)
  endif

  CS%id_du_dt_str_visc_rem = register_diag_field('ocean_model', 'du_dt_str_visc_rem', diag%axesCuL, Time, &
      'Zonal Acceleration from Surface Wind Stresses multiplied by viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_du_dt_str_visc_rem > 0) then
    call safe_alloc_ptr(ADp%du_dt_str,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%visc_rem_u,IsdB,IedB,jsd,jed,nz)
  endif

  CS%id_dv_dt_str_visc_rem = register_diag_field('ocean_model', 'dv_dt_str_visc_rem', diag%axesCvL, Time, &
      'Meridional Acceleration from Surface Wind Stresses multiplied by viscous remnant', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  if (CS%id_dv_dt_str_visc_rem > 0) then
    call safe_alloc_ptr(ADp%dv_dt_str,isd,ied,JsdB,JedB,nz)
    call safe_alloc_ptr(ADp%visc_rem_v,isd,ied,JsdB,JedB,nz)
  endif

  if ((len_trim(CS%u_trunc_file) > 0) .or. (len_trim(CS%v_trunc_file) > 0)) &
    call PointAccel_init(MIS, Time, G, param_file, diag, dirs, CS%PointAccel_CSp)

end subroutine vertvisc_init

!> Update the CFL truncation value as a function of time.
!! If called with the optional argument activate=.true., record the
!! value of Time as the beginning of the ramp period.
subroutine updateCFLtruncationValue(Time, CS, US, activate)
  type(time_type), target, intent(in)    :: Time     !< Current model time
  type(vertvisc_CS),       pointer       :: CS       !< Vertical viscosity control structure
  type(unit_scale_type),   intent(in)    :: US       !< A dimensional unit scaling type
  logical, optional,       intent(in)    :: activate !< Specify whether to record the value of
                                                     !! Time as the beginning of the ramp period

  ! Local variables
  real :: deltaTime ! The time since CS%rampStartTime [T ~> s], which may be negative.
  real :: wghtA     ! The relative weight of the final value [nondim]
  character(len=12) :: msg

  if (CS%truncRampTime==0.) return ! This indicates to ramping is turned off

  ! We use the optional argument to indicate this Time should be recorded as the
  ! beginning of the ramp-up period.
  if (present(activate)) then
    if (activate) then
      CS%rampStartTime = Time ! Record the current time
      CS%CFLrampingIsActivated = .true.
    endif
  endif
  if (.not.CS%CFLrampingIsActivated) return
  deltaTime = max( 0., US%s_to_T*time_type_to_real( Time - CS%rampStartTime ) )
  if (deltaTime >= CS%truncRampTime) then
    CS%CFL_trunc = CS%CFL_truncE
    CS%truncRampTime = 0. ! This turns off ramping after this call
  else
    wghtA = min( 1., deltaTime / CS%truncRampTime ) ! Linear profile in time
    !wghtA = wghtA*wghtA ! Convert linear profile to parabolic profile in time
    !wghtA = wghtA*wghtA*(3. - 2.*wghtA) ! Convert linear profile to cosine profile
    wghtA = 1. - ( (1. - wghtA)**2 ) ! Convert linear profile to inverted parabolic profile
    CS%CFL_trunc = CS%CFL_truncS + wghtA * ( CS%CFL_truncE - CS%CFL_truncS )
  endif
  write(msg(1:12),'(es12.3)') CS%CFL_trunc
  call MOM_error(NOTE, "MOM_vert_friction: updateCFLtruncationValue set CFL"// &
                       " limit to "//trim(msg))
end subroutine updateCFLtruncationValue

!> Clean up and deallocate the vertical friction module
subroutine vertvisc_end(CS)
  type(vertvisc_CS), intent(inout) :: CS  !< Vertical viscosity control structure that
                                          !! will be deallocated in this subroutine.

  if ((len_trim(CS%u_trunc_file) > 0) .or. (len_trim(CS%v_trunc_file) > 0)) &
    deallocate(CS%PointAccel_CSp)

  DEALLOC_(CS%a_u) ; DEALLOC_(CS%h_u)
  DEALLOC_(CS%a_v) ; DEALLOC_(CS%h_v)
  if (associated(CS%a1_shelf_u)) deallocate(CS%a1_shelf_u)
  if (associated(CS%a1_shelf_v)) deallocate(CS%a1_shelf_v)
  if (allocated(CS%kappa_gl90_2d)) deallocate(CS%kappa_gl90_2d)
end subroutine vertvisc_end

!> \namespace mom_vert_friction
!! \author Robert Hallberg
!! \date April 1994 - October 2006
!!
!!  The vertical diffusion of momentum is fully implicit.  This is
!!  necessary to allow for vanishingly small layers.  The coupling
!!  is based on the distance between the centers of adjacent layers,
!!  except where a layer is close to the bottom compared with a
!!  bottom boundary layer thickness when a bottom drag law is used.
!!  A stress top b.c. and a no slip bottom  b.c. are used.  There
!!  is no limit on the time step for vertvisc.
!!
!!  Near the bottom, the horizontal thickness interpolation scheme
!!  changes to an upwind biased estimate to control the effect of
!!  spurious Montgomery potential gradients at the bottom where
!!  nearly massless layers layers ride over the topography.  Within a
!!  few boundary layer depths of the bottom, the harmonic mean
!!  thickness (i.e. (2 h+ h-) / (h+ + h-) ) is used if the velocity
!!  is from the thinner side and the arithmetic mean thickness
!!  (i.e. (h+ + h-)/2) is used if the velocity is from the thicker
!!  side.  Both of these thickness estimates are second order
!!  accurate.  Above this the arithmetic mean thickness is used.
!!
!!  In addition, vertvisc truncates any velocity component that exceeds a
!!  maximum CFL number to a fraction of this value.  This basically keeps
!!  instabilities spatially localized.  The number of times the velocity is
!!  truncated is reported each time the energies are saved, and if
!!  exceeds CS%Maxtrunc the model will stop itself and change the time
!!  to a large value.  This has proven very useful in (1) diagnosing
!!  model failures and (2) letting the model settle down to a
!!  meaningful integration from a poorly specified initial condition.
!!
!!  The same code is used for the two velocity components, by
!!  indirectly referencing the velocities and defining a handful of
!!  direction-specific defined variables.
!!
!!  Macros written all in capital letters are defined in MOM_memory.h.
!!
!!     A small fragment of the grid is shown below:
!! \verbatim
!!    j+1  x ^ x ^ x   At x:  q
!!    j+1  > o > o >   At ^:  v, frhatv, tauy
!!    j    x ^ x ^ x   At >:  u, frhatu, taux
!!    j    > o > o >   At o:  h
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!! \endverbatim
!!
!!  The boundaries always run through q grid points (x).
end module MOM_vert_friction
