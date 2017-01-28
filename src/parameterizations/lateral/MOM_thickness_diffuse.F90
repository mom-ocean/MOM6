!> Thickness diffusion (or Gent McWilliams)
module MOM_thickness_diffuse

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,             only : hchksum, uchksum, vchksum
use MOM_diag_mediator,         only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator,         only : register_diag_field, safe_alloc_ptr, time_type
use MOM_diag_mediator,         only : diag_update_remap_grids
use MOM_error_handler,         only : MOM_error, FATAL, WARNING
use MOM_EOS,                   only : calculate_density, calculate_density_derivs
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_interface_heights,     only : find_eta
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types,            only : MEKE_type
use MOM_variables,             only : thermo_var_ptrs, cont_diag_ptrs
use MOM_verticalGrid,          only : verticalGrid_type

implicit none ; private

public thickness_diffuse, thickness_diffuse_init, thickness_diffuse_end
public vert_fill_TS

#include <MOM_memory.h>

!> Control structure for thickness diffusion
type, public :: thickness_diffuse_CS ; private
  real    :: Khth                !< Background interface depth diffusivity (m2 s-1)
  real    :: Khth_Slope_Cff      !< Slope dependence coefficient of Khth (m2 s-1)
  real    :: max_Khth_CFL        !< Maximum value of the diffusive CFL for thickness diffusion
  real    :: Khth_Min            !< Minimum value of Khth (m2 s-1)
  real    :: Khth_Max            !< Maximum value of Khth (m2 s-1), or 0 for no max
  real    :: slope_max           !< Slopes steeper than slope_max are limited in some way.
  real    :: kappa_smooth        !< Vertical diffusivity used to interpolate more
                                 !! sensible values of T & S into thin layers.
  logical :: thickness_diffuse   !< If true, interfaces heights are diffused.
  logical :: detangle_interfaces !< If true, add 3-d structured interface height
                                 !! diffusivities to horizontally smooth jagged layers.
  real    :: detangle_time       !< If detangle_interfaces is true, this is the
                                 !! timescale over which maximally jagged grid-scale
                                 !! thickness variations are suppressed.  This must be
                                 !! longer than DT, or 0 (the default) to use DT.
  integer :: nkml                !< number of layers within mixed layer
  logical :: debug               !< write verbose checksums for debugging purposes
  type(diag_ctrl), pointer :: diag ! structure used to regulate timing of diagnostics
  real, pointer :: GMwork(:,:)       => NULL()  !< Work by thickness diffusivity (W m-2)
  real, pointer :: diagSlopeX(:,:,:) => NULL()  !< Diagnostic: zonal neutral slope (nondim)
  real, pointer :: diagSlopeY(:,:,:) => NULL()  !< Diagnostic: zonal neutral slope (nondim)

  !>@{
  !! Diagnostic identifier
  integer :: id_uhGM    = -1, id_vhGM    = -1, id_GMwork = -1
  integer :: id_KH_u    = -1, id_KH_v    = -1, id_KH_t   = -1
  integer :: id_KH_u1   = -1, id_KH_v1   = -1, id_KH_t1  = -1
  integer :: id_slope_x = -1, id_slope_y = -1
  !>@}
 ! integer :: id_sfn_slope_x = -1, id_sfn_slope_y = -1, id_sfn_x = -1, id_sfn_y = -1
end type thickness_diffuse_CS

contains

!> Calculates thickness diffusion coefficients and applies thickness diffusion to layer
!! thicknesses, h. Diffusivities are limited to ensure stability.
!! Also returns along-layer mass fluxes used in the continuity equation.
subroutine thickness_diffuse(h, uhtr, vhtr, tv, dt, G, GV, MEKE, VarMix, CDp, CS)
  type(ocean_grid_type),                     intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(inout) :: h      !< Layer thickness (m or kg/m2)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: uhtr   !< Accumulated zonal mass flux (m2 H)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: vhtr   !< Accumulated meridional mass flux (m2 H)
  type(thermo_var_ptrs),                     intent(in)    :: tv     !< Thermodynamics structure
  real,                                      intent(in)    :: dt     !< Time increment (s)
  type(MEKE_type),                           intent(inout) :: MEKE   !< MEKE control structure
  type(VarMix_CS),                           pointer       :: VarMix !< Variable mixing coefficients
  type(cont_diag_ptrs),                      intent(inout) :: CDp    !< Diagnostics for the continuity equation
  type(thickness_diffuse_CS),                pointer       :: CS     !< Control structure for thickness diffusion
  ! Local variables
  real :: e(SZI_(G), SZJ_(G), SZK_(G)+1) ! heights of interfaces, relative to mean
                                         ! sea level,in H units, positive up.
  real :: uhD(SZIB_(G), SZJ_(G), SZK_(G)) ! uhD & vhD are the diffusive u*h &
  real :: vhD(SZI_(G), SZJB_(G), SZK_(G)) ! v*h fluxes (m2 H s-1)

  real, dimension(SZIB_(G), SZJ_(G), SZK_(G)+1) :: &
    KH_u, &       ! interface height diffusivities in u-columns (m2 s-1)
    int_slope_u   ! A nondimensional ratio from 0 to 1 that gives the relative
                  ! weighting of the interface slopes to that calculated also
                  ! using density gradients at u points.  The physically correct
                  ! slopes occur at 0, while 1 is used for numerical closures.
  real, dimension(SZI_(G), SZJB_(G), SZK_(G)+1) :: &
    KH_v, &       ! interface height diffusivities in v-columns (m2 s-1)
    int_slope_v   ! A nondimensional ratio from 0 to 1 that gives the relative
                  ! weighting of the interface slopes to that calculated also
                  ! using density gradients at v points.  The physically correct
                  ! slopes occur at 0, while 1 is used for numerical closures.
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: &
    KH_t          ! diagnosed diffusivity at tracer points (m^2/s)

  real, dimension(SZIB_(G), SZJ_(G)) :: &
    KH_u_CFL      ! The maximum stable interface height diffusivity at u grid points (m2 s-1)
  real, dimension(SZI_(G), SZJB_(G)) :: &
    KH_v_CFL      ! The maximum stable interface height diffusivity at v grid points (m2 s-1)
  real :: Khth_Loc_u(SZIB_(G), SZJ_(G))
  real :: Khth_Loc(SZIB_(G), SZJB_(G))  ! locally calculated thickness diffusivity (m2/s)
  real :: H_to_m, m_to_H   ! Local copies of unit conversion factors.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in H.
  logical :: use_VarMix, Resoln_scaled, use_stored_slopes, khth_use_ebt_struct
  integer :: i, j, k, is, ie, js, je, nz
  logical :: MEKE_not_null
  real :: hu(SZI_(G), SZJ_(G))       ! u-thickness (H)
  real :: hv(SZI_(G), SZJ_(G))       ! v-thickness (H)
  real :: KH_u_lay(SZI_(G), SZJ_(G)) ! layer ave thickness diffusivities (m2/sec)
  real :: KH_v_lay(SZI_(G), SZJ_(G)) ! layer ave thickness diffusivities (m2/sec)

  if (.not. ASSOCIATED(CS)) call MOM_error(FATAL, "MOM_thickness_diffuse:"// &
         "Module must be initialized before it is used.")
  MEKE_not_null = (LOC(MEKE) .NE. 0)

  if ((.not.CS%thickness_diffuse) .or. &
       .not.( CS%Khth > 0.0 .or. associated(VarMix) .or. MEKE_not_null ) ) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  h_neglect = GV%H_subroundoff
  H_to_m = GV%H_to_m ; m_to_H = GV%m_to_H

  if (MEKE_not_null) then
    if (ASSOCIATED(MEKE%GM_src)) then
      do j=js,je ; do i=is,ie ; MEKE%GM_src(i,j) = 0. ; enddo ; enddo
    endif
  endif

  use_VarMix = .false. ; Resoln_scaled = .false. ; use_stored_slopes = .false.
  khth_use_ebt_struct = .false.
  if (Associated(VarMix)) then
    use_VarMix = VarMix%use_variable_mixing
    Resoln_scaled = VarMix%Resoln_scaled_KhTh
    use_stored_slopes = VarMix%use_stored_slopes
    khth_use_ebt_struct = VarMix%khth_use_ebt_struct
  endif

!$OMP parallel do default(none) shared(is,ie,js,je,KH_u_CFL,dt,G,CS)
  do j=js,je ; do I=is-1,ie
    KH_u_CFL(I,j) = (0.25*CS%max_Khth_CFL) /  &
      (dt*(G%IdxCu(I,j)*G%IdxCu(I,j) + G%IdyCu(I,j)*G%IdyCu(I,j)))
  enddo ; enddo
!$OMP parallel do default(none) shared(is,ie,js,je,KH_v_CFL,dt,G,CS)
  do j=js-1,je ; do I=is,ie
    KH_v_CFL(i,J) = (0.25*CS%max_Khth_CFL) / &
      (dt*(G%IdxCv(i,J)*G%IdxCv(i,J) + G%IdyCv(i,J)*G%IdyCv(i,J)))
  enddo ; enddo

  call find_eta(h, tv, GV%g_Earth, G, GV, e, halo_size=1)

  ! Set the diffusivities.
!$OMP parallel default(none) shared(is,ie,js,je,Khth_Loc_u,CS,use_VarMix,VarMix,    &
!$OMP                               MEKE_not_null,MEKE,Resoln_scaled,KH_u,          &
!$OMP                               KH_u_CFL,nz,Khth_Loc,KH_v,KH_v_CFL,int_slope_u, &
!$OMP                               int_slope_v,khth_use_ebt_struct)
!$OMP do
  do j=js,je; do I=is-1,ie
    Khth_Loc_u(I,j) = CS%Khth
  enddo ; enddo

  if (use_VarMix) then
!$OMP do
    do j=js,je ; do I=is-1,ie
      Khth_Loc_u(I,j) = Khth_Loc_u(I,j) + CS%KHTH_Slope_Cff*VarMix%L2u(I,j)*VarMix%SN_u(I,j)
    enddo ; enddo
  endif

  if (MEKE_not_null) then ; if (associated(MEKE%Kh)) then
!$OMP do
    do j=js,je ; do I=is-1,ie
      Khth_Loc_u(I,j) = Khth_Loc_u(I,j) + MEKE%KhTh_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i+1,j))
    enddo ; enddo
  endif ; endif

  if (Resoln_scaled) then
!$OMP do
    do j=js,je; do I=is-1,ie
      Khth_Loc_u(I,j) = Khth_Loc_u(I,j) * VarMix%Res_fn_u(I,j)
    enddo ; enddo
  endif

  if (CS%Khth_Max > 0) then
!$OMP do
    do j=js,je; do I=is-1,ie
      Khth_Loc_u(I,j) = max(CS%Khth_min, min(Khth_Loc_u(I,j),CS%Khth_Max))
    enddo ; enddo
  else
!$OMP do
    do j=js,je; do I=is-1,ie
      Khth_Loc_u(I,j) = max(CS%Khth_min, Khth_Loc_u(I,j))
    enddo ; enddo
  endif
!$OMP do
  do j=js,je; do I=is-1,ie
    KH_u(I,j,1) = min(KH_u_CFL(I,j), Khth_Loc_u(I,j))
  enddo ; enddo

  if (khth_use_ebt_struct) then
!$OMP do
    do K=2,nz+1 ; do j=js,je ; do I=is-1,ie
      KH_u(I,j,K) = KH_u(I,j,1) * 0.5 * ( VarMix%ebt_struct(i,j,k-1) + VarMix%ebt_struct(i+1,j,k-1) )
    enddo ; enddo ; enddo
  else
!$OMP do
    do K=2,nz+1 ; do j=js,je ; do I=is-1,ie
      KH_u(I,j,K) = KH_u(I,j,1)
    enddo ; enddo ; enddo
  endif

!$OMP do
  do J=js-1,je ; do i=is,ie
    Khth_Loc(i,j) = CS%Khth
  enddo ; enddo

  if (use_VarMix) then
!$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_Loc(i,j) = Khth_Loc(i,j) + CS%KHTH_Slope_Cff*VarMix%L2v(i,J)*VarMix%SN_v(i,J)
    enddo ; enddo
  endif
  if (MEKE_not_null) then ; if (associated(MEKE%Kh)) then
!$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_Loc(i,j) = Khth_Loc(i,j) + MEKE%KhTh_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i,j+1))
    enddo ; enddo
  endif ; endif

  if (Resoln_scaled) then
!$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_Loc(i,j) = Khth_Loc(i,j) * VarMix%Res_fn_v(i,J)
    enddo ; enddo
  endif

  if (CS%Khth_Max > 0) then
!$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_Loc(i,j) = max(CS%Khth_min, min(Khth_Loc(i,j),CS%Khth_Max))
    enddo ; enddo
  else
!$OMP do
    do J=js-1,je ; do i=is,ie
      Khth_Loc(i,j) = max(CS%Khth_min, Khth_Loc(i,j))
    enddo ; enddo
  endif

  if (CS%max_Khth_CFL > 0.0) then
!$OMP do
    do J=js-1,je ; do i=is,ie
      KH_v(i,J,1) = min(KH_v_CFL(i,J), Khth_Loc(i,j))
    enddo ; enddo
  endif
  if (khth_use_ebt_struct) then
!$OMP do
    do K=2,nz+1 ; do J=js-1,je ; do i=is,ie
      KH_v(i,J,K) = KH_v(i,J,1) * 0.5 * ( VarMix%ebt_struct(i,j,k-1) + VarMix%ebt_struct(i,j+1,k-1) )
    enddo ; enddo ; enddo
  else
!$OMP do
    do K=2,nz+1 ; do J=js-1,je ; do i=is,ie
      KH_v(i,J,K) = KH_v(i,J,1)
    enddo ; enddo ; enddo
  endif
!$OMP do
  do K=1,nz+1 ; do j=js,je ; do I=is-1,ie ; int_slope_u(I,j,K) = 0.0 ; enddo ; enddo ; enddo
!$OMP do
  do K=1,nz+1 ; do J=js-1,je ; do i=is,ie ; int_slope_v(i,J,K) = 0.0 ; enddo ; enddo ; enddo
!$OMP end parallel

  if (CS%detangle_interfaces) then
    call add_detangling_Kh(h, e, Kh_u, Kh_v, KH_u_CFL, KH_v_CFL, tv, dt, G, GV, &
                           CS, int_slope_u, int_slope_v)
  endif

  if (CS%debug) then
    call uchksum(Kh_u(:,:,:),"Kh_u",G%HI,haloshift=0)
    call vchksum(Kh_v(:,:,:),"Kh_v",G%HI,haloshift=0)
    call uchksum(int_slope_u(:,:,:),"int_slope_u",G%HI,haloshift=0)
    call vchksum(int_slope_v(:,:,:),"int_slope_v",G%HI,haloshift=0)
    call hchksum(h(:,:,:)*H_to_m,"thickness_diffuse_1 h",G%HI,haloshift=1)
    call hchksum(e(:,:,:),"thickness_diffuse_1 e",G%HI,haloshift=1)
    if (use_stored_slopes) then
      call uchksum(VarMix%slope_x(:,:,:),"VarMix%slope_x",G%HI,haloshift=0)
      call vchksum(VarMix%slope_y(:,:,:),"VarMix%slope_y",G%HI,haloshift=0)
    endif
    if (associated(tv%eqn_of_state)) then
      call hchksum(tv%T(:,:,:),"thickness_diffuse T",G%HI,haloshift=1)
      call hchksum(tv%S(:,:,:),"thickness_diffuse S",G%HI,haloshift=1)
    endif
  endif

  ! Calculate uhD, vhD from h, e, KH_u, KH_v, tv%T/S
  if (use_stored_slopes) then
    call thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, dt, G, GV, MEKE, CS, &
                                int_slope_u, int_slope_v, VarMix%slope_x, VarMix%slope_y)
  else
    call thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, dt, G, GV, MEKE, CS, &
                                int_slope_u, int_slope_v)
  endif
!$OMP parallel do default(none) shared(is,ie,js,je,nz,uhtr,uhD,dt,vhtr,CDp,vhD,h,G,GV)
  do k=1,nz
    do j=js,je ; do I=is-1,ie
      uhtr(I,j,k) = uhtr(I,j,k) + uhD(I,j,k)*dt
      if (ASSOCIATED(CDp%uhGM)) CDp%uhGM(I,j,k) = uhD(I,j,k)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      vhtr(i,J,k) = vhtr(i,J,k) + vhD(i,J,k)*dt
      if (ASSOCIATED(CDp%vhGM)) CDp%vhGM(i,J,k) = vhD(i,J,k)
    enddo ; enddo
    do j=js,je ; do i=is,ie
      h(i,j,k) = h(i,j,k) - dt * G%IareaT(i,j) * &
          ((uhD(I,j,k) - uhD(I-1,j,k)) + (vhD(i,J,k) - vhD(i,J-1,k)))
      if (h(i,j,k) < GV%Angstrom) h(i,j,k) = GV%Angstrom
    enddo ; enddo
  enddo

  ! Whenever thickness changes let the diag manager know, target grids
  ! for vertical remapping may need to be regenerated.
  ! This needs to happen after the H update and before the next post_data.
  call diag_update_remap_grids(CS%diag)

  if (MEKE_not_null .AND. ASSOCIATED(VarMix)) then
    if (ASSOCIATED(MEKE%Rd_dx_h) .and. ASSOCIATED(VarMix%Rd_dx_h)) then
!$OMP parallel do default(none) shared(is,ie,js,je,MEKE,VarMix)
      do j=js,je ; do i=is,ie
        MEKE%Rd_dx_h(i,j) = VarMix%Rd_dx_h(i,j)
      enddo ; enddo
    endif
  endif

  if (CS%debug) then
    call uchksum(uhD(:,:,:)*H_to_m,"thickness_diffuse uhD",G%HI,haloshift=0)
    call vchksum(vhD(:,:,:)*H_to_m,"thickness_diffuse vhD",G%HI,haloshift=0)
    call uchksum(uhtr(:,:,:)*H_to_m,"thickness_diffuse uhtr",G%HI,haloshift=0)
    call vchksum(vhtr(:,:,:)*H_to_m,"thickness_diffuse vhtr",G%HI,haloshift=0)
    call hchksum(h(:,:,:)*H_to_m,"thickness_diffuse h",G%HI,haloshift=0)
  endif

  ! offer diagnostic fields for averaging
  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_uhGM > 0)   call post_data(CS%id_uhGM, CDp%uhGM, CS%diag)
    if (CS%id_vhGM > 0)   call post_data(CS%id_vhGM, CDp%vhGM, CS%diag)
    if (CS%id_GMwork > 0) call post_data(CS%id_GMwork, CS%GMwork, CS%diag)
    if (CS%id_KH_u > 0)   call post_data(CS%id_KH_u, KH_u, CS%diag)
    if (CS%id_KH_v > 0)   call post_data(CS%id_KH_v, KH_v, CS%diag)
    if (CS%id_KH_u1 > 0)  call post_data(CS%id_KH_u1, KH_u(:,:,1), CS%diag)
    if (CS%id_KH_v1 > 0)  call post_data(CS%id_KH_v1, KH_v(:,:,1), CS%diag)

    ! Diagnose diffusivity at T-cell point.  Do simple average, rather than
    ! thickness-weighted average, in order that KH_t is depth-independent
    ! in the case where KH_u and KH_v are depth independent.  Otherwise,
    ! if use thickess weighted average, the variations of thickness with
    ! depth will place a spurious depth dependence to the diagnosed KH_t.
    if (CS%id_KH_t > 0 .or. CS%id_KH_t1 > 0) then
      do k=1,nz
        ! thicknesses across u and v faces, converted to 0/1 mask;
        ! layer average of the interface diffusivities KH_u and KH_v
        do j=js,je ; do I=is-1,ie
          hu(I,j)       = 2.0*h(i,j,k)*h(i+1,j,k)/(h(i,j,k)+h(i+1,j,k)+h_neglect)
          if(hu(I,j) /= 0.0) hu(I,j) = 1.0
          KH_u_lay(I,j) = 0.5*(KH_u(I,j,k)+KH_u(I,j,k+1))
        enddo ; enddo
        do J=js-1,je ; do i=is,ie
          hv(i,J)       = 2.0*h(i,j,k)*h(i,j+1,k)/(h(i,j,k)+h(i,j+1,k)+h_neglect)
          if(hv(i,J) /= 0.0) hv(i,J) = 1.0
          KH_v_lay(i,J) = 0.5*(KH_v(i,J,k)+KH_v(i,J,k+1))
        enddo ; enddo
        ! diagnose diffusivity at T-point
        do j=js,je ; do i=is,ie
          KH_t(i,j,k) = ((hu(I-1,j)*KH_u_lay(i-1,j)+hu(I,j)*KH_u_lay(I,j))  &
                        +(hv(i,J-1)*KH_v_lay(i,J-1)+hv(i,J)*KH_v_lay(i,J))) &
                       / (hu(I-1,j)+hu(I,j)+hv(i,J-1)+hv(i,J)+h_neglect)
        enddo ; enddo
      enddo
      if(CS%id_KH_t  > 0) call post_data(CS%id_KH_t,  KH_t,        CS%diag)
      if(CS%id_KH_t1 > 0) call post_data(CS%id_KH_t1, KH_t(:,:,1), CS%diag)
    endif

  endif

end subroutine thickness_diffuse

!> Calculates parameterized layer transports for use in the continuity equation.
!! Fluxes are limited to give positive definite thicknesses.
!! Called by thickness_diffuse().
subroutine thickness_diffuse_full(h, e, Kh_u, Kh_v, tv, uhD, vhD, dt, G, GV, MEKE, &
                                  CS, int_slope_u, int_slope_v, slope_x, slope_y)
  type(ocean_grid_type),                       intent(in)  :: G      !< Ocean grid structure
  type(verticalGrid_type),                     intent(in)  :: GV     !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),    intent(in)  :: h      !< Layer thickness (m or kg/m2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1),  intent(in)  :: e      !< Interface positions (m)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), intent(in)  :: Kh_u   !< Thickness diffusivity on interfaces at u points (m2/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), intent(in)  :: Kh_v   !< Thickness diffusivity on interfaces at v points (m2/s)
  type(thermo_var_ptrs),                       intent(in)  :: tv     !< Thermodynamics structure
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)),   intent(out) :: uhD    !< Zonal mass fluxes (m3/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)),   intent(out) :: vhD    !< Meridional mass fluxes (m3/s)
  real,                                        intent(in)  :: dt     !< Time increment (s)
  type(MEKE_type),                             intent(inout) :: MEKE !< MEKE control structue
  type(thickness_diffuse_CS),                  pointer     :: CS     !< Control structure for thickness diffusion
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), optional, intent(in)  :: int_slope_u !< Ratio that determine how much of
                                                                     !! the isopycnal slopes are taken directly from the
                                                                     !! interface slopes without consideration of density gradients.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), optional, intent(in)  :: int_slope_v !< Ratio that determine how much of
                                                                     !! the isopycnal slopes are taken directly from the
                                                                     !! interface slopes without consideration of density gradients.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), optional, intent(in)  :: slope_x
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), optional, intent(in)  :: slope_y
  ! Local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: &
    T, &          ! The temperature (or density) in C, with the values in
                  ! in massless layers filled vertically by diffusion.
    S, &          ! The filled salinity, in PSU, with the values in
                  ! in massless layers filled vertically by diffusion.
    Rho, &        ! Density itself, when a nonlinear equation of state is
                  ! not in use.
    h_avail, &    ! The mass available for diffusion out of each face, divided
                  ! by dt, in m3 s-1.
    h_frac        ! The fraction of the mass in the column above the bottom
                  ! interface of a layer that is within a layer, ND. 0<h_frac<=1
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)+1) :: &
    pres, &       ! The pressure at an interface, in Pa.
    h_avail_rsum  ! The running sum of h_avail above an interface, in m3 s-1.
  real, dimension(SZIB_(G)) :: &
    drho_dT_u, &  ! The derivatives of density with temperature and
    drho_dS_u     ! salinity at u points, in kg m-3 K-1 and kg m-3 psu-1.
  real, dimension(SZI_(G)) :: &
    drho_dT_v, &  ! The derivatives of density with temperature and
    drho_dS_v     ! salinity at v points, in kg m-3 K-1 and kg m-3 psu-1.
  real :: uhtot(SZIB_(G), SZJ_(G))  ! The vertical sum of uhD, in m3 s-1.
  real :: vhtot(SZI_(G), SZJB_(G))  ! The vertical sum of vhD, in m3 s-1.
  real, dimension(SZIB_(G)) :: &
    T_u, S_u, &   ! Temperature, salinity, and pressure on the interface at
    pres_u        ! the u-point in the horizontal.
  real, dimension(SZI_(G)) :: &
    T_v, S_v, &   ! Temperature, salinity, and pressure on the interface at
    pres_v        ! the v-point in the horizontal.
  real :: Work_u(SZIB_(G), SZJ_(G)) ! The work being done by the thickness
  real :: Work_v(SZI_(G), SZJB_(G)) ! diffusion integrated over a cell, in W.
  real :: Work_h                    ! The work averaged over an h-cell in W m-2.
  real :: I4dt                      ! 1 / 4 dt
  real :: drdiA, drdiB  ! Along layer zonal- and meridional- potential density
  real :: drdjA, drdjB  ! gradients in the layers above (A) and below(B) the
                        ! interface times the grid spacing, in kg m-3.
  real :: drdkL, drdkR  ! Vertical density differences across an interface,
                        ! in kg m-3.
  real :: hg2A, hg2B, hg2L, hg2R
  real :: haA, haB, haL, haR
  real :: dzaL, dzaR
  real :: wtA, wtB, wtL, wtR
  real :: drdx, drdy, drdz  ! Zonal, meridional, and vertical density gradients,
                            ! in units of kg m-4.
  real :: Sfn_est       ! Two preliminary estimates (before limiting) of the
  real :: Sfn_unlim     ! overturning streamfunction, both in m3 s-1.
  real :: Sfn           ! The overturning streamfunction, in m3 s-1.
  real :: Sfn_safe      ! The streamfunction that goes linearly back to 0 at the
                        ! top.  This is a good thing to use when the slope is
                        ! so large as to be meaningless.
  real :: Slope         ! The slope of density surfaces, calculated in a way
                        ! that is always between -1 and 1.
  real :: mag_grad2     ! The squared magnitude of the 3-d density gradient, in kg2 m-8.
  real :: slope2_Ratio  ! The ratio of the slope squared to slope_max squared.
  real :: I_slope_max2  ! The inverse of slope_max squared, nondimensional.
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in H.
  real :: h_neglect2    ! h_neglect^2, in H2.
  real :: dz_neglect    ! A thickness in m that is so small it is usually lost
                        ! in roundoff and can be neglected, in m.
  real :: G_scale       ! The gravitational accerlation times the conversion
                        ! factor from thickness to m, in m s-2 or m4 s-2 kg-1.
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  logical :: find_work  ! If true, find the change in energy due to the fluxes.
  integer :: nk_linear  ! The number of layers over which the streamfunction
                        ! goes to 0.
  real :: H_to_m, m_to_H   ! Local copies of unit conversion factors.

! Diagnostics that should be eliminated altogether later...
 ! real, dimension(SZIB_(G), SZJ_(G), SZK_(G)+1) :: sfn_x, sfn_slope_x
 ! real, dimension(SZI_(G), SZJB_(G), SZK_(G)+1) :: sfn_y, sfn_slope_y
  logical :: MEKE_not_null, present_int_slope_u, present_int_slope_v
  logical :: present_slope_x, present_slope_y, calc_derivatives
  integer :: is, ie, js, je, nz, IsdB
  integer :: i, j, k
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke ; IsdB = G%IsdB

  H_to_m = GV%H_to_m ; m_to_H = GV%m_to_H
  I4dt = 0.25 / dt
  I_slope_max2 = 1.0 / (CS%slope_max**2)
  G_scale = GV%g_Earth * H_to_m
  h_neglect = GV%H_subroundoff ; h_neglect2 = h_neglect**2
  dz_neglect = GV%H_subroundoff*H_to_m

  use_EOS = associated(tv%eqn_of_state)
  MEKE_not_null = (LOC(MEKE) .NE. 0)
  present_int_slope_u = PRESENT(int_slope_u)
  present_int_slope_v = PRESENT(int_slope_v)
  present_slope_x = PRESENT(slope_x)
  present_slope_y = PRESENT(slope_y)

  nk_linear = max(GV%nkml, 1)

  find_work = .false.
  if (MEKE_not_null) find_work = ASSOCIATED(MEKE%GM_src)
  find_work = (ASSOCIATED(CS%GMwork) .or. find_work)

  if (use_EOS) then
    call vert_fill_TS(h, tv%T, tv%S, CS%kappa_smooth, dt, T, S, G, GV, 1)
  endif

!$OMP parallel default(none) shared(is,ie,js,je,h_avail_rsum,pres,h_avail,I4dt, &
!$OMP                               G,GV,h,h_frac,nz,uhtot,Work_u,vhtot,Work_v )
  ! Find the maximum and minimum permitted streamfunction.
!$OMP do
  do j=js-1,je+1 ; do i=is-1,ie+1
    h_avail_rsum(i,j,1) = 0.0
    pres(i,j,1) = 0.0  ! ### This should be atmospheric pressure.

    h_avail(i,j,1) = max(I4dt*G%areaT(i,j)*(h(i,j,1)-GV%Angstrom),0.0)
    h_avail_rsum(i,j,2) = h_avail(i,j,1)
    h_frac(i,j,1) = 1.0
    pres(i,j,2) = pres(i,j,1) + GV%H_to_Pa*h(i,j,1)
  enddo ; enddo
!$OMP do
  do j=js-1,je+1
    do k=2,nz ; do i=is-1,ie+1
      h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-GV%Angstrom),0.0)
      h_avail_rsum(i,j,k+1) = h_avail_rsum(i,j,k) + h_avail(i,j,k)
      h_frac(i,j,k) = 0.0 ; if (h_avail(i,j,k) > 0.0) &
        h_frac(i,j,k) = h_avail(i,j,k) / h_avail_rsum(i,j,k+1)
      pres(i,j,K+1) = pres(i,j,K) + GV%H_to_Pa*h(i,j,k)
    enddo ; enddo
  enddo
!$OMP do
  do j=js,je ; do I=is-1,ie
    uhtot(I,j) = 0.0 ; Work_u(I,j) = 0.0
    ! sfn_x(I,j,1) = 0.0 ; sfn_slope_x(I,j,1) = 0.0
    ! sfn_x(I,j,nz+1) = 0.0 ; sfn_slope_x(I,j,nz+1) = 0.0
  enddo ; enddo
!$OMP do
  do J=js-1,je ; do i=is,ie
    vhtot(i,J) = 0.0 ; Work_v(i,J) = 0.0
    ! sfn_y(i,J,1) = 0.0 ; sfn_slope_y(i,J,1) = 0.0
    ! sfn_y(i,J,nz+1) = 0.0 ; sfn_slope_y(i,J,nz+1) = 0.0
  enddo ; enddo
!$OMP end parallel

!$OMP parallel do default(none) shared(nz,is,ie,js,je,find_work,use_EOS,G,GV,pres,T,S, &
!$OMP                                  nk_linear,IsdB,tv,h,h_neglect,e,dz_neglect,  &
!$OMP                                  I_slope_max2,h_neglect2,present_int_slope_u, &
!$OMP                                  int_slope_u,KH_u,uhtot,h_frac,h_avail_rsum,  &
!$OMP                                  uhD,h_avail,G_scale,work_u,CS,slope_x,       &
!$OMP                                  slope_y,present_slope_x,H_to_m,m_to_H) &
!$OMP                          private(drdiA,drdiB,drdkL,drdkR,pres_u,T_u,S_u,      &
!$OMP                                  drho_dT_u,drho_dS_u,hg2A,hg2B,hg2L,hg2R,haA, &
!$OMP                                  haB,haL,haR,dzaL,dzaR,wtA,wtB,wtL,wtR,drdz,  &
!$OMP                                  drdx,mag_grad2,Slope,slope2_Ratio,Sfn_unlim, &
!$OMP                                  Sfn_safe,Sfn_est,Sfn,calc_derivatives)
  do j=js,je ; do K=nz,2,-1
    if (find_work .and. .not.(use_EOS)) then
      drdiA = 0.0 ; drdiB = 0.0
!       drdkL = GV%g_prime(k) ; drdkR = GV%g_prime(k)
      drdkL = GV%Rlay(k)-GV%Rlay(k-1) ; drdkR = GV%Rlay(k)-GV%Rlay(k-1)
    endif

    calc_derivatives = use_EOS .and. (k >= nk_linear) .and. &
                (find_work .or. .not. present_slope_x)

    ! Calculate the zonal fluxes and gradients.
    if (calc_derivatives) then
      do I=is-1,ie
        pres_u(I) = 0.5*(pres(i,j,K) + pres(i+1,j,K))
        T_u(I) = 0.25*((T(i,j,k) + T(i+1,j,k)) + (T(i,j,k-1) + T(i+1,j,k-1)))
        S_u(I) = 0.25*((S(i,j,k) + S(i+1,j,k)) + (S(i,j,k-1) + S(i+1,j,k-1)))
      enddo
      call calculate_density_derivs(T_u, S_u, pres_u, drho_dT_u, &
                   drho_dS_u, (is-IsdB+1)-1, ie-is+2, tv%eqn_of_state)
    endif

    do I=is-1,ie
      if (calc_derivatives) then
        ! Estimate the horizontal density gradients along layers.
        drdiA = drho_dT_u(I) * (T(i+1,j,k-1)-T(i,j,k-1)) + &
                drho_dS_u(I) * (S(i+1,j,k-1)-S(i,j,k-1))
        drdiB = drho_dT_u(I) * (T(i+1,j,k)-T(i,j,k)) + &
                drho_dS_u(I) * (S(i+1,j,k)-S(i,j,k))

        ! Estimate the vertical density gradients times the grid spacing.
        drdkL = (drho_dT_u(I) * (T(i,j,k)-T(i,j,k-1)) + &
                 drho_dS_u(I) * (S(i,j,k)-S(i,j,k-1)))
        drdkR = (drho_dT_u(I) * (T(i+1,j,k)-T(i+1,j,k-1)) + &
                 drho_dS_u(I) * (S(i+1,j,k)-S(i+1,j,k-1)))
      endif

      if (k > nk_linear) then
        if (use_EOS) then
          if (present_slope_x) then
            Slope = slope_x(I,j,k)
            slope2_Ratio = Slope**2 * I_slope_max2
          else
            hg2A = h(i,j,k-1)*h(i+1,j,k-1) + h_neglect2
            hg2B = h(i,j,k)*h(i+1,j,k) + h_neglect2
            hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
            hg2R = h(i+1,j,k-1)*h(i+1,j,k) + h_neglect2
            haA = 0.5*(h(i,j,k-1) + h(i+1,j,k-1))
            haB = 0.5*(h(i,j,k) + h(i+1,j,k)) + h_neglect
            haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
            haR = 0.5*(h(i+1,j,k-1) + h(i+1,j,k)) + h_neglect
            if (GV%Boussinesq) then
              dzaL = haL * H_to_m ; dzaR = haR * H_to_m
            else
              dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
              dzaR = 0.5*(e(i+1,j,K-1) - e(i+1,j,K+1)) + dz_neglect
            endif
            ! Use the harmonic mean thicknesses to weight the horizontal gradients.
            ! These unnormalized weights have been rearranged to minimize divisions.
            wtA = hg2A*haB ; wtB = hg2B*haA
            wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

            drdz = (wtL * drdkL + wtR * drdkR) / (dzaL*wtL + dzaR*wtR)
            ! The expression for drdz above is mathematically equivalent to:
            !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
            !          ((hg2L/haL) + (hg2R/haR))
            ! This is the gradient of density along geopotentials.
            drdx = ((wtA * drdiA + wtB * drdiB) / (wtA + wtB) - &
                    drdz * (e(i,j,K)-e(i+1,j,K))) * G%IdxCu(I,j)

            ! This estimate of slope is accurate for small slopes, but bounded
            ! to be between -1 and 1.
            mag_grad2 = drdx**2 + drdz**2
            if (mag_grad2 > 0.0) then
              Slope = drdx / sqrt(mag_grad2)
              slope2_Ratio = Slope**2 * I_slope_max2
            else ! Just in case mag_grad2 = 0 ever.
              Slope = 0.0
              slope2_Ratio = 1.0e20  ! Force the use of the safe streamfunction.
            endif
          endif

          ! Adjust real slope by weights that bias towards slope of interfaces
          ! that ignore density gradients along layers.
          if (present_int_slope_u) then
            Slope = (1.0 - int_slope_u(I,j,K)) * Slope + &
                    int_slope_u(I,j,K) * ((e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j))
            slope2_Ratio = (1.0 - int_slope_u(I,j,K)) * slope2_Ratio
          endif
          if (CS%id_slope_x > 0) CS%diagSlopeX(I,j,k) = Slope

          ! Estimate the streamfunction at each interface.
          Sfn_unlim = -((KH_u(I,j,K)*G%dy_Cu(I,j))*Slope) * m_to_H
          if (uhtot(I,j) <= 0.0) then
            ! The transport that must balance the transport below is positive.
            Sfn_safe = uhtot(I,j) * (1.0 - h_frac(i,j,k))
          else !  (uhtot(I,j) > 0.0)
            Sfn_safe = uhtot(I,j) * (1.0 - h_frac(i+1,j,k))
          endif

          ! Avoid moving dense water upslope from below the level of
          ! the bottom on the receiving side.
          if (Sfn_unlim > 0.0) then ! The flow below this interface is positive.
            if (e(i,j,K) < e(i+1,j,nz+1)) then
              Sfn_unlim = 0.0 ! This is not uhtot, because it may compensate for
                              ! deeper flow in very unusual cases.
            elseif (e(i+1,j,nz+1) > e(i,j,K+1)) then
              ! Scale the transport with the fraction of the donor layer above
              ! the bottom on the receiving side.
              Sfn_unlim = Sfn_unlim * ((e(i,j,K) - e(i+1,j,nz+1)) / &
                                       ((e(i,j,K) - e(i,j,K+1)) + dz_neglect))
            endif
          else
            if (e(i+1,j,K) < e(i,j,nz+1)) then ; Sfn_unlim = 0.0
            elseif (e(i,j,nz+1) > e(i+1,j,K+1)) then
              Sfn_unlim = Sfn_unlim * ((e(i+1,j,K) - e(i,j,nz+1)) / &
                                     ((e(i+1,j,K) - e(i+1,j,K+1)) + dz_neglect))
            endif
          endif

          Sfn_est = (Sfn_unlim + slope2_Ratio*Sfn_safe) / (1.0 + slope2_Ratio)
        else  ! With .not.use_EOS, the layers are constant density.
          if (present_slope_x) then
            Slope = slope_x(I,j,k)
          else
            Slope = ((e(i,j,K)-e(i+1,j,K))*G%IdxCu(I,j)) * m_to_H
          endif
          Sfn_est = (KH_u(I,j,K)*G%dy_Cu(I,j)) * Slope
                  !  ((e(i,j,K)-e(i+1,j,K))*G%IdxCu(I,j))) * m_to_H
          if (CS%id_slope_x > 0) CS%diagSlopeX(I,j,k) = Slope
        endif

        ! Make sure that there is enough mass above to allow the streamfunction
        ! to satisfy the boundary condition of 0 at the surface.
        Sfn = min(max(Sfn_est, -h_avail_rsum(i,j,K)), h_avail_rsum(i+1,j,K))

        ! The actual transport is limited by the mass available in the two
        ! neighboring grid cells.
        uhD(I,j,k) = max(min((Sfn - uhtot(I,j)), h_avail(i,j,k)), &
                         -h_avail(i+1,j,k))

 !       sfn_x(I,j,K) = max(min(Sfn, uhtot(I,j)+h_avail(i,j,k)), &
 !                          uhtot(I,j)-h_avail(i+1,j,K))
 !       sfn_slope_x(I,j,K) = max(uhtot(I,j)-h_avail(i+1,j,k), &
 !                                min(uhtot(I,j)+h_avail(i,j,k), &
 !             min(h_avail_rsum(i+1,j,K), max(-h_avail_rsum(i,j,K), &
 !             (KH_u(I,j,K)*G%dy_Cu(I,j)) * ((e(i,j,K)-e(i+1,j,K))*G%IdxCu(I,j)) )) ))
      else ! k <= nk_linear
        ! Balance the deeper flow with a return flow uniformly distributed
        ! though the remaining near-surface layers.  This is the same as
        ! using Sfn_safe above.  There is no need to apply the limiters in
        ! this case.
        if (uhtot(I,j) <= 0.0) then
          uhD(I,j,k) = -uhtot(I,j) * h_frac(i,j,k)
        else !  (uhtot(I,j) > 0.0)
          uhD(I,j,k) = -uhtot(I,j) * h_frac(i+1,j,k)
        endif

 !       sfn_x(I,j,K) = sfn_x(I,j,K+1) + uhD(I,j,k)
 !       if (sfn_slope_x(I,j,K+1) <= 0.0) then
 !         sfn_slope_x(I,j,K) = sfn_slope_x(I,j,K+1) * (1.0 - h_frac(i,j,k))
 !       else
 !         sfn_slope_x(I,j,K) = sfn_slope_x(I,j,K+1) * (1.0 - h_frac(i+1,j,k))
 !       endif
      endif

      uhtot(I,j) = uhtot(I,j) + uhD(I,j,k)

      if (find_work) then
        !   This is the energy tendency based on the original profiles, and does
        ! not include any nonlinear terms due to a finite time step (which would
        ! involve interactions between the fluxes through the different faces.
        !   A second order centered estimate is used for the density transfered
        ! between water columns.

        Work_u(I,j) = Work_u(I,j) + G_scale * &
          ( uhtot(I,j) * (drdkR * e(i+1,j,K) - drdkL * e(i,j,K)) - &
            (uhD(I,j,K) * drdiB) * 0.25 * &
            ((e(i,j,K) + e(i,j,K+1)) + (e(i+1,j,K) + e(i+1,j,K+1))) )
      endif

    enddo
  enddo ; enddo ! end of j-loop

    ! Calculate the meridional fluxes and gradients.
!$OMP parallel do default(none) shared(nz,is,ie,js,je,find_work,use_EOS,G,GV,pres,T,S, &
!$OMP                                  nk_linear,IsdB,tv,h,h_neglect,e,dz_neglect,  &
!$OMP                                  I_slope_max2,h_neglect2,present_int_slope_v, &
!$OMP                                  int_slope_v,KH_v,vhtot,h_frac,h_avail_rsum,  &
!$OMP                                  vhD,h_avail,G_scale,Work_v,CS,slope_y,       &
!$OMP                                  present_slope_y,present_slope_x,m_to_H,H_to_m) &
!$OMP                          private(drdjA,drdjB,drdkL,drdkR,pres_v,T_v,S_v,      &
!$OMP                                  drho_dT_v,drho_dS_v,hg2A,hg2B,hg2L,hg2R,haA, &
!$OMP                                  haB,haL,haR,dzaL,dzaR,wtA,wtB,wtL,wtR,drdz,  &
!$OMP                                  drdy,mag_grad2,Slope,slope2_Ratio,Sfn_unlim, &
!$OMP                                  Sfn_safe,Sfn_est,Sfn,calc_derivatives)
  do J=js-1,je ; do K=nz,2,-1
    if (find_work .and. .not.(use_EOS)) then
      drdjA = 0.0 ; drdjB = 0.0
!       drdkL = GV%g_prime(k) ; drdkR = GV%g_prime(k)
      drdkL = GV%Rlay(k)-GV%Rlay(k-1) ; drdkR = GV%Rlay(k)-GV%Rlay(k-1)
    endif

    calc_derivatives = use_EOS .and. (k >= nk_linear) .and. &
                (find_work .or. .not. present_slope_x)

    if (calc_derivatives) then
      do i=is,ie
        pres_v(i) = 0.5*(pres(i,j,K) + pres(i,j+1,K))
        T_v(i) = 0.25*((T(i,j,k) + T(i,j+1,k)) + (T(i,j,k-1) + T(i,j+1,k-1)))
        S_v(i) = 0.25*((S(i,j,k) + S(i,j+1,k)) + (S(i,j,k-1) + S(i,j+1,k-1)))
      enddo
      call calculate_density_derivs(T_v, S_v, pres_v, drho_dT_v, &
                   drho_dS_v, is, ie-is+1, tv%eqn_of_state)
    endif
    do i=is,ie
      if (calc_derivatives) then
        ! Estimate the horizontal density gradients along layers.
        drdjA = drho_dT_v(i) * (T(i,j+1,k-1)-T(i,j,k-1)) + &
                drho_dS_v(i) * (S(i,j+1,k-1)-S(i,j,k-1))
        drdjB = drho_dT_v(i) * (T(i,j+1,k)-T(i,j,k)) + &
                drho_dS_v(i) * (S(i,j+1,k)-S(i,j,k))

        ! Estimate the vertical density gradients times the grid spacing.
        drdkL = (drho_dT_v(i) * (T(i,j,k)-T(i,j,k-1)) + &
                 drho_dS_v(i) * (S(i,j,k)-S(i,j,k-1)))
        drdkR = (drho_dT_v(i) * (T(i,j+1,k)-T(i,j+1,k-1)) + &
                 drho_dS_v(i) * (S(i,j+1,k)-S(i,j+1,k-1)))
      endif

      if (k > nk_linear) then
        if (use_EOS) then
          if (present_slope_y) then
            Slope = slope_y(i,J,k)
            slope2_Ratio = Slope**2 * I_slope_max2
          else
            hg2A = h(i,j,k-1)*h(i,j+1,k-1) + h_neglect2
            hg2B = h(i,j,k)*h(i,j+1,k) + h_neglect2
            hg2L = h(i,j,k-1)*h(i,j,k) + h_neglect2
            hg2R = h(i,j+1,k-1)*h(i,j+1,k) + h_neglect2
            haA = 0.5*(h(i,j,k-1) + h(i,j+1,k-1)) + h_neglect
            haB = 0.5*(h(i,j,k) + h(i,j+1,k)) + h_neglect
            haL = 0.5*(h(i,j,k-1) + h(i,j,k)) + h_neglect
            haR = 0.5*(h(i,j+1,k-1) + h(i,j+1,k)) + h_neglect
            if (GV%Boussinesq) then
              dzaL = haL * H_to_m ; dzaR = haR * H_to_m
            else
              dzaL = 0.5*(e(i,j,K-1) - e(i,j,K+1)) + dz_neglect
              dzaR = 0.5*(e(i,j+1,K-1) - e(i,j+1,K+1)) + dz_neglect
            endif
            ! Use the harmonic mean thicknesses to weight the horizontal gradients.
            ! These unnormalized weights have been rearranged to minimize divisions.
            wtA = hg2A*haB ; wtB = hg2B*haA
            wtL = hg2L*(haR*dzaR) ; wtR = hg2R*(haL*dzaL)

            drdz = (wtL * drdkL + wtR * drdkR) / (dzaL*wtL + dzaR*wtR)
            ! The expression for drdz above is mathematically equivalent to:
            !   drdz = ((hg2L/haL) * drdkL/dzaL + (hg2R/haR) * drdkR/dzaR) / &
            !          ((hg2L/haL) + (hg2R/haR))
            ! This is the gradient of density along geopotentials.
            drdy = ((wtA * drdjA + wtB * drdjB) / (wtA + wtB) - &
                    drdz * (e(i,j,K)-e(i,j+1,K))) * G%IdyCv(i,J)

            ! This estimate of slope is accurate for small slopes, but bounded
            ! to be between -1 and 1.
            mag_grad2 = drdy**2 + drdz**2
            if (mag_grad2 > 0.0) then
              Slope = drdy / sqrt(mag_grad2)
              slope2_Ratio = Slope**2 * I_slope_max2
            else ! Just in case mag_grad2 = 0 ever.
              Slope = 0.0
              slope2_Ratio = 1.0e20  ! Force the use of the safe streamfunction.
            endif
          endif

          ! Adjust real slope by weights that bias towards slope of interfaces
          ! that ignore density gradients along layers.
          if (present_int_slope_v) then
            Slope = (1.0 - int_slope_v(i,J,K)) * Slope + &
                    int_slope_v(i,J,K) * ((e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J))
            slope2_Ratio = (1.0 - int_slope_v(i,J,K)) * slope2_Ratio
          endif
          if (CS%id_slope_y > 0) CS%diagSlopeY(I,j,k) = Slope

          ! Estimate the streamfunction at each interface.
          Sfn_unlim = -((KH_v(i,J,K)*G%dx_Cv(i,J))*Slope) * m_to_H
          if (vhtot(i,J) <= 0.0) then
            Sfn_safe = vhtot(i,J) * (1.0 - h_frac(i,j,k))
          else !  (vhtot(I,j) > 0.0)
            Sfn_safe = vhtot(i,J) * (1.0 - h_frac(i,j+1,k))
          endif

          ! Avoid moving dense water upslope from below the level of
          ! the bottom on the receiving side.
          if (Sfn_unlim > 0.0) then ! The flow below this interface is positive.
            if (e(i,j,K) < e(i,j+1,nz+1)) then
              Sfn_unlim = 0.0 ! This is not vhtot, because it may compensate for
                              ! deeper flow in very unusual cases.
            elseif (e(i,j+1,nz+1) > e(i,j,K+1)) then
              ! Scale the transport with the fraction of the donor layer above
              ! the bottom on the receiving side.
              Sfn_unlim = Sfn_unlim * ((e(i,j,K) - e(i,j+1,nz+1)) / &
                                       ((e(i,j,K) - e(i,j,K+1)) + dz_neglect))
            endif
          else
            if (e(i,j+1,K) < e(i,j,nz+1)) then ; Sfn_unlim = 0.0
            elseif (e(i,j,nz+1) > e(i,j+1,K+1)) then
              Sfn_unlim = Sfn_unlim * ((e(i,j+1,K) - e(i,j,nz+1)) / &
                                     ((e(i,j+1,K) - e(i,j+1,K+1)) + dz_neglect))
            endif
          endif

          ! Estimate the streamfunction at each interface.
          Sfn_est = (Sfn_unlim + slope2_Ratio*Sfn_safe) / (1.0 + slope2_Ratio)
        else      ! With .not.use_EOS, the layers are constant density.
          if (present_slope_y) then
            Slope = slope_y(i,J,k)
          else
            Slope = ((e(i,j,K)-e(i,j+1,K))*G%IdyCv(i,J)) * m_to_H
          endif
          Sfn_est = (KH_v(i,J,K)*G%dx_Cv(i,J)) * Slope
                  !  ((e(i,j,K)-e(i,j+1,K))*G%IdyCv(i,J))) * m_to_H
          if (CS%id_slope_y > 0) CS%diagSlopeY(I,j,k) = Slope
        endif

        ! Make sure that there is enough mass above to allow the streamfunction
        ! to satisfy the boundary condition of 0 at the surface.
        Sfn = min(max(Sfn_est, -h_avail_rsum(i,j,K)), h_avail_rsum(i,j+1,K))

        ! The actual transport is limited by the mass available in the two
        ! neighboring grid cells.
        vhD(i,J,k) = max(min((Sfn - vhtot(i,J)), h_avail(i,j,k)), &
                         -h_avail(i,j+1,k))

  !      sfn_y(i,J,K) = max(min(Sfn, vhtot(i,J)+h_avail(i,j,k)), &
  !                         vhtot(i,J)-h_avail(i,j+1,k))
  !      sfn_slope_y(i,J,K) = max(vhtot(i,J)-h_avail(i,j+1,k), &
  !                               min(vhtot(i,J)+h_avail(i,j,k), &
  !            min(h_avail_rsum(i,j+1,K), max(-h_avail_rsum(i,j,K), &
  !            (KH_v(i,J,K)*G%dx_Cv(i,J)) * ((e(i,j,K)-e(i,j+1,K))*G%IdyCv(i,J)) )) ))
      else  ! k <= nk_linear
        ! Balance the deeper flow with a return flow uniformly distributed
        ! though the remaining near-surface layers.
        if (vhtot(i,J) <= 0.0) then
          vhD(i,J,k) = -vhtot(i,J) * h_frac(i,j,k)
        else !  (vhtot(i,J) > 0.0)
          vhD(i,J,k) = -vhtot(i,J) * h_frac(i,j+1,k)
        endif

  !      sfn_y(i,J,K) = sfn_y(i,J,K+1) + vhD(i,J,k)
  !      if (sfn_slope_y(i,J,K+1) <= 0.0) then
  !        sfn_slope_y(i,J,K) = sfn_slope_y(i,J,K+1) * (1.0 - h_frac(i,j,k))
  !      else
  !        sfn_slope_y(i,J,K) = sfn_slope_y(i,J,K+1) * (1.0 - h_frac(i,j+1,k))
  !      endif
      endif

      vhtot(i,J) = vhtot(i,J)  + vhD(i,J,k)

      if (find_work) then
        !   This is the energy tendency based on the original profiles, and does
        ! not include any nonlinear terms due to a finite time step (which would
        ! involve interactions between the fluxes through the different faces.
        !   A second order centered estimate is used for the density transfered
        ! between water columns.

        Work_v(i,J) = Work_v(i,J) + G_scale * &
          ( vhtot(i,J) * (drdkR * e(i,j+1,K) - drdkL * e(i,j,K)) - &
           (vhD(i,J,K) * drdjB) * 0.25 * &
           ((e(i,j,K) + e(i,j,K+1)) + (e(i,j+1,K) + e(i,j+1,K+1))) )
      endif
    enddo
  enddo ; enddo! j-loop

  ! In layer 1, enforce the boundary conditions that Sfn(z=0) = 0.0
  if (.not.find_work .or. .not.(use_EOS)) then
    do j=js,je ; do I=is-1,ie ; uhD(I,j,1) = -uhtot(I,j) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; vhD(i,J,1) = -vhtot(i,J) ; enddo ; enddo
  else
!$OMP parallel do default(none) shared(is,ie,js,je,pres,T,S,IsdB,tv,uhD,uhtot, &
!$OMP                                  Work_u,G_scale,use_EOS,e)   &
!$OMP                          private(pres_u,T_u,S_u,drho_dT_u,drho_dS_u,drdiB)
    do j=js,je
      if (use_EOS) then
        do I=is-1,ie
          pres_u(I) = 0.5*(pres(i,j,1) + pres(i+1,j,1))
          T_u(I) = 0.5*(T(i,j,1) + T(i+1,j,1))
          S_u(I) = 0.5*(S(i,j,1) + S(i+1,j,1))
        enddo
        call calculate_density_derivs(T_u, S_u, pres_u, drho_dT_u, &
                   drho_dS_u, (is-IsdB+1)-1, ie-is+2, tv%eqn_of_state)
      endif
      do I=is-1,ie
        uhD(I,j,1) = -uhtot(I,j)

        if (use_EOS) then
          drdiB = drho_dT_u(I) * (T(i+1,j,1)-T(i,j,1)) + &
                  drho_dS_u(I) * (S(i+1,j,1)-S(i,j,1))
        endif
        Work_u(I,j) = Work_u(I,j) + G_scale * ( (uhD(I,j,1) * drdiB) * 0.25 * &
            ((e(i,j,1) + e(i,j,2)) + (e(i+1,j,1) + e(i+1,j,2))) )

      enddo
    enddo

    do J=js-1,je
      if (use_EOS) then
        do i=is,ie
          pres_v(i) = 0.5*(pres(i,j,1) + pres(i,j+1,1))
          T_v(i) = 0.5*(T(i,j,1) + T(i,j+1,1))
          S_v(i) = 0.5*(S(i,j,1) + S(i,j+1,1))
        enddo
        call calculate_density_derivs(T_v, S_v, pres_v, drho_dT_v, &
                   drho_dS_v, is, ie-is+1, tv%eqn_of_state)
      endif
      do i=is,ie
        vhD(i,J,1) = -vhtot(i,J)

        if (use_EOS) then
          drdjB = drho_dT_v(i) * (T(i,j+1,1)-T(i,j,1)) + &
                  drho_dS_v(i) * (S(i,j+1,1)-S(i,j,1))
        endif
        Work_v(i,J) = Work_v(i,J) - G_scale * ( (vhD(i,J,1) * drdjB) * 0.25 * &
            ((e(i,j,1) + e(i,j,2)) + (e(i,j+1,1) + e(i,j+1,2))) )
      enddo
    enddo
  endif

  if (find_work) then ; do j=js,je ; do i=is,ie
    ! Note that the units of Work_v and Work_u are W, while Work_h is W m-2.
    Work_h = 0.5 * G%IareaT(i,j) * &
      ((Work_u(I-1,j) + Work_u(I,j)) + (Work_v(i,J-1) + Work_v(i,J)))
    if (ASSOCIATED(CS%GMwork)) CS%GMwork(i,j) = Work_h
    if (MEKE_not_null) then ; if (ASSOCIATED(MEKE%GM_src)) then
      MEKE%GM_src(i,j) = MEKE%GM_src(i,j) + Work_h
    endif ; endif
  enddo ; enddo ; endif

  if (CS%id_slope_x > 0) call post_data(CS%id_slope_x, CS%diagSlopeX, CS%diag)
  if (CS%id_slope_y > 0) call post_data(CS%id_slope_y, CS%diagSlopeY, CS%diag)
 ! if (CS%id_sfn_x > 0) call post_data(CS%id_sfn_x, sfn_x, CS%diag)
 ! if (CS%id_sfn_y > 0) call post_data(CS%id_sfn_y, sfn_y, CS%diag)
 ! if (CS%id_sfn_slope_x > 0) call post_data(CS%id_sfn_slope_x, sfn_slope_x, CS%diag)
 ! if (CS%id_sfn_slope_y > 0) call post_data(CS%id_sfn_slope_y, sfn_slope_y, CS%diag)

end subroutine thickness_diffuse_full

!> Modifies thickness diffusivities to untangle layer structures
subroutine add_detangling_Kh(h, e, Kh_u, Kh_v, KH_u_CFL, KH_v_CFL, tv, dt, G, GV, CS, &
                             int_slope_u, int_slope_v)
  type(ocean_grid_type),                       intent(in)    :: G    !< Ocean grid structure
  type(verticalGrid_type),                     intent(in)    :: GV   !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),    intent(in)    :: h    !< Layer thickness (H)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1),  intent(in)    :: e    !< Interface positions (m)
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: Kh_u !< Thickness diffusivity on interfaces at u points (m2/s)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), intent(inout) :: Kh_v !< Thickness diffusivity on interfaces at u points (m2/s)
  real, dimension(SZIB_(G),SZJ_(G)),           intent(in)    :: Kh_u_CFL !< Maximum stable thickness diffusivity at u points (m2/s)
  real, dimension(SZI_(G),SZJB_(G)),           intent(in)    :: Kh_v_CFL !< Maximum stable thickness diffusivity at v points (m2/s)
  type(thermo_var_ptrs),                       intent(in)    :: tv   !< Thermodynamics structure
  real,                                        intent(in)    :: dt   !< Time increment (s)
  type(thickness_diffuse_CS),                  pointer       :: CS   !< Control structure for thickness diffusion
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)+1), intent(inout) :: int_slope_u !< Ratio that determine how much of
                                                                     !! the isopycnal slopes are taken directly from the
                                                                     !! interface slopes without consideration of density gradients.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)+1), intent(inout) :: int_slope_v !< Ratio that determine how much of
                                                                     !! the isopycnal slopes are taken directly from the
                                                                     !! interface slopes without consideration of density gradients.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    de_top     ! The distances between the top of a layer and the top of the
               ! region where the detangling is applied, in H.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: &
    Kh_lay_u   ! The tentative interface height diffusivity for each layer at
               ! u points, in m2 s-1.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: &
    Kh_lay_v   ! The tentative interface height diffusivity for each layer at
               ! v points, in m2 s-1.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    de_bot     ! The distances from the bottom of the region where the
               ! detangling is applied, in H.
  real :: h1, h2    ! The thinner and thicker surrounding thicknesses, in H,
                    ! with the thinner modified near the boundaries to mask out
                    ! thickness variations due to topography, etc.
  real :: jag_Rat   ! The nondimensional jaggedness ratio for a layer, going
                    ! from 0 (smooth) to 1 (jagged).  This is the difference
                    ! between the arithmetic and harmonic mean thicknesses
                    ! normalized by the arithmetic mean thickness.
  real :: Kh_scale  ! A ratio by which Kh_u_CFL is scaled for maximally jagged
                    ! layers, nondim.
  real :: Kh_det    ! The detangling diffusivity, in m2 s-1.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in H.

  real :: I_sl      ! The absolute value of the larger in magnitude of the slopes
                    ! above and below.
  real :: Rsl       ! The ratio of the smaller magnitude slope to the larger
                    ! magnitude one, ND. 0 <= Rsl <1.
  real :: IRsl      ! The (limited) inverse of Rsl, ND. 1 < IRsl <= 1e9.
  real :: dH        ! The thickness gradient divided by the damping timescale
                    ! and the ratio of the face length to the adjacent cell
                    ! areas for comparability with the diffusivities, in m2 s-1.
  real :: adH       ! The absolute value of dH, in m2 s-1.
  real :: sign      ! 1 or -1, with the same sign as the layer thickness gradient.
  real :: sl_K      ! The sign-corrected slope of the interface above, ND.
  real :: sl_Kp1    ! The sign-corrected slope of the interface below, ND.
  real :: I_sl_K    ! The (limited) inverse of sl_K, ND.
  real :: I_sl_Kp1  ! The (limited) inverse of sl_Kp1, ND.
  real :: I_4t      ! A quarter of inverse of the damping timescale, in s-1.
  real :: Fn_R      ! A function of Rsl, such that Rsl < Fn_R < 1.
  real :: denom, I_denom ! A denominator and its inverse, various units.
  real :: Kh_min    ! A local floor on the diffusivity, in m2 s-1.
  real :: Kh_max    ! A local ceiling on the diffusivity, in m2 s-1.
  real :: wt1, wt2  ! Nondimensional weights.
  !   Variables used only in testing code.
  ! real, dimension(SZK_(G)) :: uh_here
  ! real, dimension(SZK_(G)+1) :: Sfn
  real :: dKh       ! An increment in the diffusivity, in m2 s-1.

  real, dimension(SZIB_(G),SZK_(G)+1) :: &
    Kh_bg, &        ! The background (floor) value of Kh, in m2 s-1.
    Kh, &           ! The tentative value of Kh, in m2 s-1.
    Kh_detangle, &  ! The detangling diffusivity that could be used, in m2 s-1.
    Kh_min_max_p, & ! The smallest ceiling that can be placed on Kh(I,K)
                    ! based on the value of Kh(I,K+1), in m2 s-1.
    Kh_min_max_m, & ! The smallest ceiling that can be placed on Kh(I,K)
                    ! based on the value of Kh(I,K-1), in m2 s-1.
    ! The following are variables that define the relationships between
    ! successive values of Kh.
    ! Search for Kh that satisfy...
    !    Kh(I,K) >= Kh_min_m(I,K)*Kh(I,K-1) + Kh0_min_m(I,K)
    !    Kh(I,K) >= Kh_min_p(I,K)*Kh(I,K+1) + Kh0_min_p(I,K)
    !    Kh(I,K) <= Kh_max_m(I,K)*Kh(I,K-1) + Kh0_max_m(I,K)
    !    Kh(I,K) <= Kh_max_p(I,K)*Kh(I,K+1) + Kh0_max_p(I,K)
    Kh_min_m , &   ! See above, ND.
    Kh0_min_m , &  ! See above, in m2 s-1.
    Kh_max_m , &   ! See above, ND.
    Kh0_max_m, &   ! See above, in m2 s-1.
    Kh_min_p , &   ! See above, ND.
    Kh0_min_p , &  ! See above, in m2 s-1.
    Kh_max_p , &   ! See above, ND.
    Kh0_max_p      ! See above, in m2 s-1.
  real, dimension(SZIB_(G)) :: &
    Kh_max_max  ! The maximum diffusivity permitted in a column.
  logical, dimension(SZIB_(G)) :: &
    do_i        ! If true, work on a column.
  integer :: i, j, k, n, ish, jsh, is, ie, js, je, nz, k_top
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  k_top = GV%nk_rho_varies + 1
  h_neglect = GV%H_subroundoff
  !   The 0.5 is because we are not using uniform weightings, but are
  ! distributing the diffusivities more effectively (with wt1 & wt2), but this
  ! means that the additions to a single interface can be up to twice as large.
  Kh_scale = 0.5
  if (CS%detangle_time > dt) Kh_scale = 0.5 * dt / CS%detangle_time

  do j=js-1,je+1 ; do i=is-1,ie+1
    de_top(i,j,k_top) = 0.0 ; de_bot(i,j) = 0.0
  enddo ; enddo
  do k=k_top+1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
    de_top(i,j,k) = de_top(i,j,k-1) + h(i,j,k-1)
  enddo ; enddo ; enddo

  do j=js,je ; do I=is-1,ie
    Kh_lay_u(I,j,nz) = 0.0 ; Kh_lay_u(I,j,k_top) = 0.0
  enddo ; enddo
  do J=js-1,je ; do i=is,ie
    Kh_lay_v(i,J,nz) = 0.0 ; Kh_lay_v(i,J,k_top) = 0.0
  enddo ; enddo

  do k=nz-1,k_top+1,-1
    ! Find the diffusivities associated with each layer.
    do j=js-1,je+1 ; do i=is-1,ie+1
      de_bot(i,j) = de_bot(i,j) + h(i,j,k+1)
    enddo ; enddo

    do j=js,je ; do I=is-1,ie ; if (G%mask2dCu(I,j) > 0.0) then
      if (h(i,j,k) > h(i+1,j,k)) then
        h2 = h(i,j,k)
        h1 = max( h(i+1,j,k), h2 - min(de_bot(i+1,j), de_top(i+1,j,k)) )
      else
        h2 = h(i+1,j,k)
        h1 = max( h(i,j,k), h2 - min(de_bot(i,j), de_top(i,j,k)) )
      endif
      jag_Rat = (h2 - h1)**2 / (h2 + h1 + h_neglect)**2
      Kh_lay_u(I,j,k) = (Kh_scale * Kh_u_CFL(I,j)) * jag_Rat**2
    endif ; enddo ; enddo

    do J=js-1,je ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.0) then
      if (h(i,j,k) > h(i,j+1,k)) then
        h2 = h(i,j,k)
        h1 = max( h(i,j+1,k), h2 - min(de_bot(i,j+1), de_top(i,j+1,k)) )
      else
        h2 = h(i,j+1,k)
        h1 = max( h(i,j,k), h2 - min(de_bot(i,j), de_top(i,j,k)) )
      endif
      jag_Rat = (h2 - h1)**2 / (h2 + h1 + h_neglect)**2
      Kh_lay_v(i,J,k) = (Kh_scale * Kh_v_CFL(i,J)) * jag_Rat**2
    endif ; enddo ; enddo
  enddo

  ! Limit the diffusivities

  I_4t = Kh_scale / (4.0*dt)

  do n=1,2
    if (n==1) then ; jsh = js ; ish = is-1
    else ; jsh = js-1 ; ish = is ; endif

    do j=jsh,je

      ! First, populate the diffusivities
      if (n==1) then ! This is a u-column.
        do i=ish,ie
          do_i(I) = (G%mask2dCu(I,j) > 0.0)
          Kh_max_max(I) = Kh_u_CFL(I,j)
        enddo
        do K=1,nz+1 ; do i=ish,ie
          Kh_bg(I,K) = Kh_u(I,j,K) ; Kh(I,K) = Kh_bg(I,K)
          Kh_min_max_p(I,K) = Kh_bg(I,K) ; Kh_min_max_m(I,K) = Kh_bg(I,K)
          Kh_detangle(I,K) = 0.0
        enddo ; enddo
      else ! This is a v-column.
        do i=ish,ie
          do_i(i) = (G%mask2dCv(i,J) > 0.0) ; Kh_max_max(I) = Kh_v_CFL(i,J)
        enddo
        do K=1,nz+1 ; do i=ish,ie
          Kh_bg(I,K) = Kh_v(I,j,K) ; Kh(I,K) = Kh_bg(I,K)
          Kh_min_max_p(I,K) = Kh_bg(I,K) ; Kh_min_max_m(I,K) = Kh_bg(I,K)
          Kh_detangle(I,K) = 0.0
        enddo ; enddo
      endif

      ! Determine the limits on the diffusivities.
      do k=k_top,nz ; do i=ish,ie ; if (do_i(i)) then
        if (n==1) then ! This is a u-column.
          dH = 0.0
          denom = ((G%IareaT(i+1,j) + G%IareaT(i,j))*G%dy_Cu(I,j))
          !   This expression uses differences in e in place of h for better
          ! consistency with the slopes.
          if (denom > 0.0) &
            dH = I_4t * ((e(i+1,j,K) - e(i+1,j,K+1)) - &
                         (e(i,j,K) - e(i,j,K+1))) / denom
           ! dH = I_4t * (h(i+1,j,k) - h(i,j,k)) / denom

          adH = abs(dH)
          sign = 1.0 ; if (dH < 0) sign = -1.0
          sl_K = sign * (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
          sl_Kp1 = sign * (e(i+1,j,K+1)-e(i,j,K+1)) * G%IdxCu(I,j)

          ! Add the incremental diffusivites to the surrounding interfaces.
          ! Adding more to the more steeply sloping layers (as below) makes
          ! the diffusivities more than twice as effective.
          denom = (sl_K**2 + sl_Kp1**2)
          wt1 = 0.5 ; wt2 = 0.5
          if (denom > 0.0) then
            wt1 = sl_K**2 / denom ; wt2 = sl_Kp1**2 / denom
          endif
          Kh_detangle(I,K) = Kh_detangle(I,K) + wt1*Kh_lay_u(I,j,k)
          Kh_detangle(I,K+1) = Kh_detangle(I,K+1) + wt2*Kh_lay_u(I,j,k)
        else ! This is a v-column.
          dH = 0.0
          denom = ((G%IareaT(i,j+1) + G%IareaT(i,j))*G%dx_Cv(I,j))
          if (denom > 0.0) &
            dH = I_4t * ((e(i,j+1,K) - e(i,j+1,K+1)) - &
                         (e(i,j,K) - e(i,j,K+1))) / denom
           ! dH = I_4t * (h(i,j+1,k) - h(i,j,k)) / denom

          adH = abs(dH)
          sign = 1.0 ; if (dH < 0) sign = -1.0
          sl_K = sign * (e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J)
          sl_Kp1 = sign * (e(i,j+1,K+1)-e(i,j,K+1)) * G%IdyCv(i,J)

          ! Add the incremental diffusviites to the surrounding interfaces.
          ! Adding more to the more steeply sloping layers (as below) makes
          ! the diffusivities more than twice as effective.
          denom = (sl_K**2 + sl_Kp1**2)
          wt1 = 0.5 ; wt2 = 0.5
          if (denom > 0.0) then
            wt1 = sl_K**2 / denom ; wt2 = sl_Kp1**2 / denom
          endif
          Kh_detangle(I,K) = Kh_detangle(I,K) + wt1*Kh_lay_v(i,J,k)
          Kh_detangle(I,K+1) = Kh_detangle(I,K+1) + wt2*Kh_lay_v(i,J,k)
        endif

        if (adH == 0.0) then
          Kh_min_m(I,K+1) = 1.0 ; Kh0_min_m(I,K+1) = 0.0
          Kh_max_m(I,K+1) = 1.0 ; Kh0_max_m(I,K+1) = 0.0
          Kh_min_p(I,K) = 1.0 ; Kh0_min_p(I,K) = 0.0
          Kh_max_p(I,K) = 1.0 ; Kh0_max_p(I,K) = 0.0
        elseif (adH > 0.0) then
          if (sl_K <= sl_Kp1) then
            ! This case should only arise from nonlinearities in the equation of state.
            ! Treat it as though dedx(K) = dedx(K+1) & dH = 0.
            Kh_min_m(I,K+1) = 1.0 ; Kh0_min_m(I,K+1) = 0.0
            Kh_max_m(I,K+1) = 1.0 ; Kh0_max_m(I,K+1) = 0.0
            Kh_min_p(I,K) = 1.0 ; Kh0_min_p(I,K) = 0.0
            Kh_max_p(I,K) = 1.0 ; Kh0_max_p(I,K) = 0.0
          elseif (sl_K <= 0.0) then   ! Both slopes are opposite to dH
            I_sl = -1.0 / sl_Kp1
            Rsl = -sl_K * I_sl                            ! 0 <= Rsl < 1
            IRsl = 1e9 ; if (Rsl > 1e-9) IRsl = 1.0/Rsl   ! 1 < IRsl <= 1e9

            Fn_R = Rsl
            if (Kh_max_max(I) > 0) &
              Fn_R = min(sqrt(Rsl), Rsl + (adH * I_sl) / Kh_max_max(I))

            Kh_min_m(I,K+1) = Fn_R ; Kh0_min_m(I,K+1) = 0.0
            Kh_max_m(I,K+1) = Rsl ; Kh0_max_m(I,K+1) = adH * I_sl
            Kh_min_p(I,K) = IRsl ; Kh0_min_p(I,K) = -adH * (I_sl*IRsl)
            Kh_max_p(I,K) = 1.0/(Fn_R + 1.0e-30) ; Kh0_max_p(I,K) = 0.0
          elseif (sl_Kp1 < 0.0) then  ! Opposite (nonzero) signs of slopes.
            I_sl_K = 1e18 ; if (sl_K > 1e-18) I_sl_K = 1.0 / sl_K
            I_sl_Kp1 = 1e18 ; if (-sl_Kp1 > 1e-18) I_sl_Kp1 = -1.0 / sl_Kp1

            Kh_min_m(I,K+1) = 0.0 ; Kh0_min_m(I,K+1) = 0.0
            Kh_max_m(I,K+1) = - sl_K*I_sl_Kp1 ; Kh0_max_m(I,K+1) = adH*I_sl_Kp1
            Kh_min_p(I,K) = 0.0 ; Kh0_min_p(I,K) = 0.0
            Kh_max_p(I,K) = sl_Kp1*I_sl_K ; Kh0_max_p(I,K) = adH*I_sl_K

            ! This limit does not use the slope weighting so that potentially
            ! sharp gradients in diffusivities are not forced to occur.
            Kh_max = adH / (sl_K - sl_Kp1)
            Kh_min_max_p(I,K) = max(Kh_min_max_p(I,K), Kh_max)
            Kh_min_max_m(I,K+1) = max(Kh_min_max_m(I,K+1), Kh_max)
          else ! Both slopes are of the same sign as dH.
            I_sl = 1.0 / sl_K
            Rsl = sl_Kp1 * I_sl                           ! 0 <= Rsl < 1
            IRsl = 1e9 ; if (Rsl > 1e-9) IRsl = 1.0/Rsl   ! 1 < IRsl <= 1e9

            ! Rsl <= Fn_R <= 1
            Fn_R = Rsl
            if (Kh_max_max(I) > 0) &
              Fn_R = min(sqrt(Rsl), Rsl + (adH * I_sl) / Kh_max_max(I))

            Kh_min_m(I,K+1) = IRsl ; Kh0_min_m(I,K+1) = -adH * (I_sl*IRsl)
            Kh_max_m(I,K+1) = 1.0/(Fn_R + 1.0e-30) ; Kh0_max_m(I,K+1) = 0.0
            Kh_min_p(I,K) = Fn_R ; Kh0_min_p(I,K) = 0.0
            Kh_max_p(I,K) = Rsl ; Kh0_max_p(I,K) = adH * I_sl
          endif
        endif
      endif ; enddo ; enddo ! I-loop & k-loop

      do k=k_top,nz+1,nz+1-k_top ; do i=ish,ie ; if (do_i(i)) then
        ! The diffusivities at k_top and nz+1 are both fixed.
        Kh_min_m(I,k) = 0.0 ; Kh0_min_m(I,k) = 0.0
        Kh_max_m(I,k) = 0.0 ; Kh0_max_m(I,k) = 0.0
        Kh_min_p(I,k) = 0.0 ; Kh0_min_p(I,k) = 0.0
        Kh_max_p(I,k) = 0.0 ; Kh0_max_p(I,k) = 0.0
        Kh_min_max_p(I,K) = Kh_bg(I,K)
        Kh_min_max_m(I,K) = Kh_bg(I,K)
      endif ; enddo ; enddo ! I-loop and k_top/nz+1 loop

      ! Search for Kh that satisfy...
      !    Kh(I,K) >= Kh_min_m(I,K)*Kh(I,K-1) + Kh0_min_m(I,K)
      !    Kh(I,K) >= Kh_min_p(I,K)*Kh(I,K+1) + Kh0_min_p(I,K)
      !    Kh(I,K) <= Kh_max_m(I,K)*Kh(I,K-1) + Kh0_max_m(I,K)
      !    Kh(I,K) <= Kh_max_p(I,K)*Kh(I,K+1) + Kh0_max_p(I,K)

      ! Increase the diffusivies to satisfy the min constraints.
      ! All non-zero min constraints on one diffusivity are max constraints on another.
      do K=k_top+1,nz ; do i=ish,ie ; if (do_i(i)) then
        Kh(I,K) = max(Kh_bg(I,K), Kh_detangle(I,K), &
                      min(Kh_min_m(I,K)*Kh(I,K-1) + Kh0_min_m(I,K), Kh(I,K-1)))

        if (Kh0_max_m(I,K) > Kh_bg(I,K)) Kh(I,K) = min(Kh(I,K), Kh0_max_m(I,K))
        if (Kh0_max_p(I,K) > Kh_bg(I,K)) Kh(I,K) = min(Kh(I,K), Kh0_max_p(I,K))
      endif ; enddo ; enddo ! I-loop & k-loop
      ! This is still true... do i=ish,ie ; Kh(I,nz+1) = Kh_bg(I,nz+1) ; enddo
      do K=nz,k_top+1,-1 ; do i=ish,ie ; if (do_i(i)) then
        Kh(I,k) = max(Kh(I,K), min(Kh_min_p(I,K)*Kh(I,K+1) + Kh0_min_p(I,K), Kh(I,K+1)))

        Kh_max = max(Kh_min_max_p(I,K), Kh_max_p(I,K)*Kh(I,K+1) + Kh0_max_p(I,K))
        Kh(I,k) = min(Kh(I,k), Kh_max)
      endif ; enddo ; enddo ! I-loop & k-loop
      !  All non-zero min constraints on one diffusivity are max constraints on
      ! another layer, so the min constraints can now be discounted.

      ! Decrease the diffusivities to satisfy the max constraints.
        do K=k_top+1,nz ; do i=ish,ie ; if (do_i(i)) then
          Kh_max = max(Kh_min_max_m(I,K), Kh_max_m(I,K)*Kh(I,K-1) + Kh0_max_m(I,K))
          if (Kh(I,k) > Kh_max) Kh(I,k) = Kh_Max
        endif ; enddo ; enddo  ! i- and K-loops

      ! This code tests the solutions...
!     do i=ish,ie
!       Sfn(:) = 0.0 ; uh_here(:) = 0.0
!       do K=k_top,nz
!         if ((Kh(i,K) > Kh_bg(i,K)) .or. (Kh(i,K+1) > Kh_bg(i,K+1))) then
!           if (n==1) then ! u-point.
!             if ((h(i+1,j,k) - h(i,j,k)) * &
!                 ((e(i+1,j,K)-e(i+1,j,K+1)) - (e(i,j,K)-e(i,j,K+1))) > 0.0) then
!               Sfn(K) = -Kh(i,K) * (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
!               Sfn(K+1) = -Kh(i,K+1) * (e(i+1,j,K+1)-e(i,j,K+1)) * G%IdxCu(I,j)
!               uh_here(k) = (Sfn(K) - Sfn(K+1))*G%dy_Cu(I,j)
!               if (abs(uh_here(k))*min(G%IareaT(i,j), G%IareaT(i+1,j)) > &
!                   (1e-10*GV%m_to_H)) then
!                 if (uh_here(k) * (h(i+1,j,k) - h(i,j,k)) > 0.0) then
!                   call MOM_error(WARNING, &
!                          "Corrective u-transport is up the thickness gradient.", .true.)
!                 endif
!                 if (((h(i,j,k) - 4.0*dt*G%IareaT(i,j)*uh_here(k)) - &
!                      (h(i+1,j,k) + 4.0*dt*G%IareaT(i+1,j)*uh_here(k))) * &
!                     (h(i,j,k) - h(i+1,j,k)) < 0.0) then
!                   call MOM_error(WARNING, &
!                          "Corrective u-transport is too large.", .true.)
!                 endif
!               endif
!             endif
!           else ! v-point
!             if ((h(i,j+1,k) - h(i,j,k)) * &
!                 ((e(i,j+1,K)-e(i,j+1,K+1)) - (e(i,j,K)-e(i,j,K+1))) > 0.0) then
!               Sfn(K) = -Kh(i,K) * (e(i,j+1,K)-e(i,j,K)) * G%IdyCv(i,J)
!               Sfn(K+1) = -Kh(i,K+1) * (e(i,j+1,K+1)-e(i,j,K+1)) * G%IdyCv(i,J)
!               uh_here(k) = (Sfn(K) - Sfn(K+1))*G%dx_Cv(i,J)
!               if (abs(uh_here(K))*min(G%IareaT(i,j), G%IareaT(i,j+1)) > &
!                   (1e-10*GV%m_to_H)) then
!                 if (uh_here(K) * (h(i,j+1,k) - h(i,j,k)) > 0.0) then
!                   call MOM_error(WARNING, &
!                          "Corrective v-transport is up the thickness gradient.", .true.)
!                 endif
!                 if (((h(i,j,k) - 4.0*dt*G%IareaT(i,j)*uh_here(K)) - &
!                      (h(i,j+1,k) + 4.0*dt*G%IareaT(i,j+1)*uh_here(K))) * &
!                     (h(i,j,k) - h(i,j+1,k)) < 0.0) then
!                   call MOM_error(WARNING, &
!                          "Corrective v-transport is too large.", .true.)
!                 endif
!               endif
!             endif
!           endif ! u- or v- selection.
!          !  de_dx(I,K) = (e(i+1,j,K)-e(i,j,K)) * G%IdxCu(I,j)
!         endif
!       enddo
!     enddo

      if (n==1) then ! This is a u-column.
        do K=k_top+1,nz ; do i=ish,ie
          if (Kh(I,K) > Kh_u(I,j,K)) then
            dKh = (Kh(I,K) - Kh_u(I,j,K))
            int_slope_u(I,j,K) = dKh / Kh(I,K)
            Kh_u(I,j,K) = Kh(I,K)
          endif
        enddo ; enddo
      else ! This is a v-column.
        do K=k_top+1,nz ; do i=ish,ie
          if (Kh(i,K) > Kh_v(i,J,K)) then
            dKh = Kh(i,K) - Kh_v(i,J,K)
            int_slope_v(i,J,K) = dKh / Kh(i,K)
            Kh_v(i,J,K) = Kh(i,K)
          endif
        enddo ; enddo
      endif

    enddo ! j-loop
  enddo  ! n-loop over u- and v- directions.

end subroutine add_detangling_Kh

!> Fills tracer values in massless layers with sensible values by diffusing
!! vertically with a (small) constant diffusivity.
subroutine vert_fill_TS(h, T_in, S_in, kappa, dt, T_f, S_f, G, GV, halo_here)
  type(ocean_grid_type),                    intent(in)  :: G     !< Ocean grid structure
  type(verticalGrid_type),                  intent(in)  :: GV    !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: h     !< Layer thickness (m or kg/m2)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: T_in  !< Input temperature (C)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: S_in  !< Input salinity (ppt)
  real,                                     intent(in)  :: kappa !< Constant diffusivity to use (m2/s)
  real,                                     intent(in)  :: dt    !< Time increment (s)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: T_f   !< Filled temperature (C)
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out) :: S_f   !< Filled salinity (ppt)
  integer,                        optional, intent(in)  :: halo_here !< Number of halo points to work on,
                                                                 !! 0 by default
  ! Local variables
  real :: ent(SZI_(G),SZK_(G)+1)   ! The diffusive entrainment (kappa*dt)/dz
                                   ! between layers in a timestep in m or kg m-2.
  real :: b1(SZI_(G)), d1(SZI_(G)) ! b1, c1, and d1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))      ! tridiagonal solver.
  real :: kap_dt_x2                ! The product of 2*kappa*dt, converted to
                                   ! the same units as h, in m2 or kg2 m-4.
  real :: h0                       ! A negligible thickness, in m or kg m-2, to
                                   ! allow for zero thicknesses.
  real :: h_neglect                ! A thickness that is so small it is usually
                                   ! lost in roundoff and can be neglected
                                   ! (m for Bouss and kg/m^2 for non-Bouss).
                                   ! 0 < h_neglect << h0.
  real :: h_tr                     ! h_tr is h at tracer points with a tiny thickness
                                   ! added to ensure positive definiteness
                                   ! (m for Bouss, kg/m^2 for non-Bouss)
  integer :: i, j, k, is, ie, js, je, nz, halo

  halo=0 ; if (present(halo_here)) halo = max(halo_here,0)

  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  nz = G%ke
  h_neglect = GV%H_subroundoff
  kap_dt_x2 = (2.0*kappa*dt)*GV%m_to_H**2
  h0 = 1.0e-16*sqrt(kappa*dt)*GV%m_to_H

  if (kap_dt_x2 <= 0.0) then
!$OMP parallel do default(none) shared(is,ie,js,je,nz,T_f,T_in,S_f,S_in)
    do k=1,nz ; do j=js,je ; do i=is,ie
      T_f(i,j,k) = T_in(i,j,k) ; S_f(i,j,k) = S_in(i,j,k)
    enddo ; enddo ; enddo
  else
!$OMP parallel do default(none) private(ent,b1,d1,c1,h_tr)   &
!$OMP             shared(is,ie,js,je,nz,kap_dt_x2,h,h0,h_neglect,T_f,S_f,T_in,S_in)
    do j=js,je
      do i=is,ie
        ent(i,2) = kap_dt_x2 / ((h(i,j,1)+h(i,j,2)) + h0)
        h_tr = h(i,j,1) + h_neglect
        b1(i) = 1.0 / (h_tr + ent(i,2))
        d1(i) = b1(i) * h(i,j,1)
        T_f(i,j,1) = (b1(i)*h_tr)*T_in(i,j,1)
        S_f(i,j,1) = (b1(i)*h_tr)*S_in(i,j,1)
      enddo
      do k=2,nz-1 ; do i=is,ie
        ent(i,K+1) = kap_dt_x2 / ((h(i,j,k)+h(i,j,k+1)) + h0)
        h_tr = h(i,j,k) + h_neglect
        c1(i,k) = ent(i,K) * b1(i)
        b1(i) = 1.0 / ((h_tr + d1(i)*ent(i,K)) + ent(i,K+1))
        d1(i) = b1(i) * (h_tr + d1(i)*ent(i,K))
        T_f(i,j,k) = b1(i) * (h_tr*T_in(i,j,k) + ent(i,K)*T_f(i,j,k-1))
        S_f(i,j,k) = b1(i) * (h_tr*S_in(i,j,k) + ent(i,K)*S_f(i,j,k-1))
      enddo ; enddo
      do i=is,ie
        c1(i,nz) = ent(i,nz) * b1(i)
        h_tr = h(i,j,nz) + h_neglect
        b1(i) = 1.0 / (h_tr + d1(i)*ent(i,nz))
        T_f(i,j,nz) = b1(i) * (h_tr*T_in(i,j,nz) + ent(i,nz)*T_f(i,j,nz-1))
        S_f(i,j,nz) = b1(i) * (h_tr*S_in(i,j,nz) + ent(i,nz)*S_f(i,j,nz-1))
      enddo
      do k=nz-1,1,-1 ; do i=is,ie
        T_f(i,j,k) = T_f(i,j,k) + c1(i,k+1)*T_f(i,j,k+1)
        S_f(i,j,k) = S_f(i,j,k) + c1(i,k+1)*S_f(i,j,k+1)
      enddo ; enddo
    enddo
  endif

end subroutine vert_fill_TS

!> Initialize the thickness diffusion module/strucutre
subroutine thickness_diffuse_init(Time, G, GV, param_file, diag, CDp, CS)
  type(time_type),         intent(in) :: Time    !< Current model time
  type(ocean_grid_type),   intent(in) :: G       !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV      !< Vertical grid structure
  type(param_file_type),   intent(in) :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(cont_diag_ptrs),    intent(inout) :: CDp  !< Continuity equation diagnostics
  type(thickness_diffuse_CS), pointer    :: CS   !< Control structure for thickness diffusion

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_thickness_diffuse" ! This module's name.
  character(len=48)  :: flux_units

  if (associated(CS)) then
    call MOM_error(WARNING, &
      "Thickness_diffuse_init called with an associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "THICKNESSDIFFUSE", CS%thickness_diffuse, &
                 "If true, interface heights are diffused with a \n"//&
                 "coefficient of KHTH.", default=.false.)
  call get_param(param_file, mod, "KHTH", CS%Khth, &
                 "The background horizontal thickness diffusivity.", &
                 units = "m2 s-1", default=0.0)
  call get_param(param_file, mod, "KHTH_SLOPE_CFF", CS%KHTH_Slope_Cff, &
                 "The nondimensional coefficient in the Visbeck formula \n"//&
                 "for the interface depth diffusivity", units="nondim", &
                 default=0.0)
  call get_param(param_file, mod, "KHTH_MIN", CS%KHTH_Min, &
                 "The minimum horizontal thickness diffusivity.", &
                 units = "m2 s-1", default=0.0)
  call get_param(param_file, mod, "KHTH_MAX", CS%KHTH_Max, &
                 "The maximum horizontal thickness diffusivity.", &
                 units = "m2 s-1", default=0.0)
  call get_param(param_file, mod, "KHTH_MAX_CFL", CS%max_Khth_CFL, &
                 "The maximum value of the local diffusive CFL ratio that \n"//&
                 "is permitted for the thickness diffusivity. 1.0 is the \n"//&
                 "marginally unstable value in a pure layered model, but \n"//&
                 "much smaller numbers (e.g. 0.1) seem to work better for \n"//&
                 "ALE-based models.", units = "nondimensional", default=0.8)
  if (CS%max_Khth_CFL < 0.0) CS%max_Khth_CFL = 0.0
  call get_param(param_file, mod, "DETANGLE_INTERFACES", CS%detangle_interfaces, &
                 "If defined add 3-d structured enhanced interface height \n"//&
                 "diffusivities to horizonally smooth jagged layers.", &
                 default=.false.)
  CS%detangle_time = 0.0
  if (CS%detangle_interfaces) &
    call get_param(param_file, mod, "DETANGLE_TIMESCALE", CS%detangle_time, &
                 "A timescale over which maximally jagged grid-scale \n"//&
                 "thickness variations are suppressed.  This must be \n"//&
                 "longer than DT, or 0 to use DT.", units = "s", default=0.0)
  call get_param(param_file, mod, "KHTH_SLOPE_MAX", CS%slope_max, &
                 "A slope beyond which the calculated isopycnal slope is \n"//&
                 "not reliable and is scaled away.", units="nondim", default=0.01)
  call get_param(param_file, mod, "KD_SMOOTH", CS%kappa_smooth, &
                 "A diapycnal diffusivity that is used to interpolate \n"//&
                 "more sensible values of T & S into thin layers.", &
                 default=1.0e-6)
  call get_param(param_file, mod, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)


  if (GV%Boussinesq) then ; flux_units = "meter3 second-1"
  else ; flux_units = "kilogram second-1" ; endif

  CS%id_uhGM = register_diag_field('ocean_model', 'uhGM', diag%axesCuL, Time, &
           'Time Mean Diffusive Zonal Thickness Flux', flux_units, &
           y_cell_method='sum', v_extensive=.true.)
  if (CS%id_uhGM > 0) call safe_alloc_ptr(CDp%uhGM,G%IsdB,G%IedB,G%jsd,G%jed,G%ke)
  CS%id_vhGM = register_diag_field('ocean_model', 'vhGM', diag%axesCvL, Time, &
           'Time Mean Diffusive Meridional Thickness Flux', flux_units, &
           x_cell_method='sum', v_extensive=.true.)
  if (CS%id_vhGM > 0) call safe_alloc_ptr(CDp%vhGM,G%isd,G%ied,G%JsdB,G%JedB,G%ke)

  CS%id_GMwork = register_diag_field('ocean_model', 'GMwork', diag%axesT1, Time,                     &
   'Integrated Tendency of Ocean Mesoscale Eddy KE from Parameterized Eddy Advection',               &
   'Watt meter-2', cmor_field_name='tnkebto',                                                        &
   cmor_long_name='Integrated Tendency of Ocean Mesoscale Eddy KE from Parameterized Eddy Advection',&
   cmor_units='W m-2',                                                                               &
   cmor_standard_name='tendency_of_ocean_eddy_kinetic_energy_content_due_to_parameterized_eddy_advection')
  if (CS%id_GMwork > 0) call safe_alloc_ptr(CS%GMwork,G%isd,G%ied,G%jsd,G%jed)

  CS%id_KH_u = register_diag_field('ocean_model', 'KHTH_u', diag%axesCui, Time, &
           'Parameterized mesoscale eddy advection diffusivity at U-point', 'meter second-2')
  CS%id_KH_v = register_diag_field('ocean_model', 'KHTH_v', diag%axesCvi, Time, &
           'Parameterized mesoscale eddy advection diffusivity at V-point', 'meter second-2')
  CS%id_KH_t = register_diag_field('ocean_model', 'KHTH_t', diag%axesTL, Time,               &
       'Ocean Tracer Diffusivity due to Parameterized Mesoscale Advection', 'meter second-2',&
   cmor_field_name='diftrblo',                                                               &
   cmor_long_name='Ocean Tracer Diffusivity due to Parameterized Mesoscale Advection',       &
   cmor_units='m2 s-1',                                                                      &
   cmor_standard_name='ocean_tracer_diffusivity_due_to_parameterized_mesoscale_advection')

  CS%id_KH_u1 = register_diag_field('ocean_model', 'KHTH_u1', diag%axesCu1, Time,         &
           'Parameterized mesoscale eddy advection diffusivity at U-points (2-D)', 'meter second-2')
  CS%id_KH_v1 = register_diag_field('ocean_model', 'KHTH_v1', diag%axesCv1, Time,         &
           'Parameterized mesoscale eddy advection diffusivity at V-points (2-D)', 'meter second-2')
  CS%id_KH_t1 = register_diag_field('ocean_model', 'KHTH_t1', diag%axesT1, Time,&
           'Parameterized mesoscale eddy advection diffusivity at T-points (2-D)', 'meter second-2')

  CS%id_slope_x =  register_diag_field('ocean_model', 'neutral_slope_x', diag%axesCui, Time, &
           'Zonal slope of neutral surface', 'nondim')
  if (CS%id_slope_x > 0) call safe_alloc_ptr(CS%diagSlopeX,G%IsdB,G%IedB,G%jsd,G%jed,G%ke+1)
  CS%id_slope_y =  register_diag_field('ocean_model', 'neutral_slope_y', diag%axesCvi, Time, &
           'Meridional slope of neutral surface', 'nondim')
  if (CS%id_slope_y > 0) call safe_alloc_ptr(CS%diagSlopeY,G%isd,G%ied,G%JsdB,G%JedB,G%ke+1)
 ! CS%id_sfn_x =  register_diag_field('ocean_model', 'sfn_x', diag%axesCui, Time, &
 !          'Parameterized Zonal Overturning Streamfunction', 'meter3 second-1')
 ! CS%id_sfn_y =  register_diag_field('ocean_model', 'sfn_y', diag%axesCvi, Time, &
 !          'Parameterized Meridional Overturning Streamfunction', 'meter3 second-1')
 ! CS%id_sfn_slope_x =  register_diag_field('ocean_model', 'sfn_sl_x', diag%axesCui, Time, &
 !          'Parameterized Zonal Overturning Streamfunction from Interface Slopes', 'meter3 second-1')
 ! CS%id_sfn_slope_y =  register_diag_field('ocean_model', 'sfn_sl_y', diag%axesCvi, Time, &
 !          'Parameterized Meridional Overturning Streamfunction from Interface Slopes', 'meter3 second-1')

end subroutine thickness_diffuse_init

!> Deallocate the thickness diffusion control structure
subroutine thickness_diffuse_end(CS)
  type(thickness_diffuse_CS), pointer :: CS   !< Control structure for thickness diffusion
  if(associated(CS)) deallocate(CS)
end subroutine thickness_diffuse_end

!> \namespace mom_thickness_diffuse
!!
!! Thickness diffusion is implemented via along-layer mass fluxes
!!
!! \f[
!! h^\dagger \leftarrow h^n - \nabla \cdot ( \vec{uh}^* )
!! \f]
!!
!! where the mass fluxes are cast as the difference in stream-functions proportional to the isoneutral slope
!!
!! \f[
!! \vec{uh}^* = \delta_k \vec{\psi} = \delta_k \left( \frac{g\kappa_h}{\rho_o} \frac{\nabla \rho}{N^2} \right)
!! \f]
!!
!! \todo Check signs of GM stream function in documentation.
!!
!! Thickness diffusivities are calculated independently at u- and v-points using the following expression
!!
!! \f[
!! \kappa_h = \left( \kappa_o + \alpha_{s} L_{s}^2 S N + \alpha_{M} \kappa_{M} \right) r(\Delta x,L_d)
!! \f]
!!
!! where \f$ S \f$ is the isoneutral slope magnitude, \f$ N \f$ is the square root of Brunt-Vaisala frequency,
!! \f$\kappa_{M}\f$ is the diffusivity calculated by the MEKE parameterization (mom_meke module) and \f$ r(\Delta x,L_d) \f$ is
!! a function of the local resolution (ratio of grid-spacing, \f$\Delta x\f$, to deformation radius, \f$L_d\f$).
!! The length \f$L_s\f$ is provided by the mom_lateral_mixing_coeffs module (enabled with
!! <code>USE_VARIABLE_MIXING=True</code>.
!!
!! The result of the above expression is subsequently bounded by minimum and maximum values, including an upper
!! diffusivity consistent with numerical stability (\f$ \kappa_{cfl} \f$ is calculated internally).
!!
!! \f[
!! \kappa_h \leftarrow \min{\left( \kappa_{max}, \kappa_{cfl}, \max{\left( \kappa_{min}, \kappa_h \right)} \right)} f(c_g,z)
!! \f]
!!
!! where \f$f(c_g,z)\f$ is a vertical structure function.
!! \f$f(c_g,z)\f$ is calculated in module mom_lateral_mixing_coeffs.
!! If <code>KHTH_USE_EBT_STRUCT=True</code> then \f$f(c_g,z)\f$ is set to look like the equivalent barotropic modal velocity structure.
!! Otherwise \f$f(c_g,z)=1\f$ and the diffusivity is independent of depth.
!!
!! In order to calculate meaningful slopes in vanished layers, temporary copies of the thermodynamic variables
!! are passed through a vertical smoother, function vert_fill_ts():
!! \f{eqnarray*}{
!! \left[ 1 + \Delta t \kappa_{smth} \frac{\partial^2}{\partial_z^2} \right] \theta & \leftarrow & \theta \\
!! \left[ 1 + \Delta t \kappa_{smth} \frac{\partial^2}{\partial_z^2} \right] s & \leftarrow & s
!! \f}
!!
!! \subsection section_khth_module_parameters Module mom_thickness_diffuse parameters
!!
!! | Symbol                | Module parameter |
!! | ------                | --------------- |
!! | -                     | <code>THICKNESSDIFFUSE</code> |
!! | \f$ \kappa_o \f$      | <code>KHTH</code> |
!! | \f$ \alpha_{s} \f$    | <code>KHTH_SLOPE_CFF</code> |
!! | \f$ \kappa_{min} \f$  | <code>KHTH_MIN</code> |
!! | \f$ \kappa_{max} \f$  | <code>KHTH_MAX</code> |
!! | -                     | <code>KHTH_MAX_CFL</code> |
!! | \f$ \kappa_{smth} \f$ | <code>KD_SMOOTH</code> |
!! | \f$ \alpha_{M} \f$    | <code>MEKE_KHTH_FAC</code> (from mom_meke module) |
!! | -                     | <code>KHTH_USE_EBT_STRUCT</code> (from mom_lateral_mixing_coeffs module) |

end module MOM_thickness_diffuse
