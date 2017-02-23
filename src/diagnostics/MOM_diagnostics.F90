module MOM_diagnostics

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
!*  By Robert Hallberg, February 2001                                  *
!*                                                                     *
!*    This subroutine calculates any requested diagnostic quantities   *
!*  that are not calculated in the various subroutines.  Diagnostic    *
!*  quantities are requested by allocating them memory.                *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, bathyT                                *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_coms,              only : reproducing_sum
use MOM_diag_mediator,     only : post_data, post_data_1d_k
use MOM_diag_mediator,     only : register_diag_field, register_scalar_field
use MOM_diag_mediator,     only : diag_ctrl, time_type, safe_alloc_ptr
use MOM_domains,           only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,           only : To_North, To_East
use MOM_EOS,               only : calculate_density, int_density_dz
use MOM_error_handler,     only : MOM_error, FATAL, WARNING
use MOM_file_parser,       only : get_param, log_version, param_file_type
use MOM_forcing_type,      only : forcing
use MOM_grid,              only : ocean_grid_type
use MOM_interface_heights, only : find_eta
use MOM_spatial_means,     only : global_area_mean, global_layer_mean, global_volume_mean
use MOM_variables,         only : thermo_var_ptrs, ocean_internal_state, p3d
use MOM_variables,         only : accel_diag_ptrs, cont_diag_ptrs
use MOM_verticalGrid,      only : verticalGrid_type
use MOM_wave_speed,        only : wave_speed, wave_speed_CS, wave_speed_init

implicit none ; private

#include <MOM_memory.h>

public calculate_diagnostic_fields
public register_time_deriv
public find_eta
public MOM_diagnostics_init
public MOM_diagnostics_end


type, public :: diagnostics_CS ; private
  real :: mono_N2_column_fraction = 0. !< The lower fraction of water column over which N2 is limited as
                                       !! monotonic for the purposes of calculating the equivalent barotropic wave speed.
  real :: mono_N2_depth = -1.          !< The depth below which N2 is limited as monotonic for the purposes of
                                       !! calculating the equivalent barotropic wave speed. (m)

  type(diag_ctrl), pointer :: diag ! structure to regulate diagnostics timing

  ! following arrays store diagnostics calculated here and unavailable outside.

  ! following fields have nz+1 levels.
  real, pointer, dimension(:,:,:) :: &
    e => NULL(), &   ! interface height (metre)
    e_D => NULL()    ! interface height above bottom (metre)

  ! following fields have nz layers.
  real, pointer, dimension(:,:,:) :: &
    du_dt => NULL(), &    ! net i-acceleration in m/s2
    dv_dt => NULL(), &    ! net j-acceleration in m/s2
    dh_dt => NULL(), &    ! thickness rate of change in (m/s) or kg/(m2*s)

    h_Rlay => NULL(),    & ! layer thicknesses in layered potential density
                           ! coordinates, in m (Bouss) or kg/m2 (non-Bouss)
    uh_Rlay   => NULL(), & ! zonal and meridional transports in layered
    vh_Rlay   => NULL(), & ! potential rho coordinates: m3/s(Bouss) kg/s(non-Bouss)
    uhGM_Rlay => NULL(), & ! zonal and meridional Gent-McWilliams transports in layered
    vhGM_Rlay => NULL(), & ! potential density coordinates, m3/s (Bouss) kg/s(non-Bouss)
    p_ebt     => NULL()    ! Equivalent barotropic modal structure

  ! following fields are 2-D.
  real, pointer, dimension(:,:) :: &
    cg1 => NULL(),       & ! first baroclinic gravity wave speed, in m s-1
    Rd1 => NULL(),       & ! first baroclinic deformation radius, in m
    cfl_cg1 => NULL(),   & ! CFL for first baroclinic gravity wave speed, nondim
    cfl_cg1_x => NULL(), & ! i-component of CFL for first baroclinic gravity wave speed, nondim
    cfl_cg1_y => NULL()    ! j-component of CFL for first baroclinic gravity wave speed, nondim

  ! arrays to hold diagnostics in the layer-integrated energy budget.
  ! all except KE have units of m3 s-3 (when Boussinesq).
  real, pointer, dimension(:,:,:) :: &
    KE        => NULL(), &  ! KE per unit mass, in m2 s-2
    dKE_dt    => NULL(), &  ! time derivative of the layer KE
    PE_to_KE  => NULL(), &  ! potential energy to KE term
    KE_CorAdv => NULL(), &  ! KE source from the combined Coriolis and
                            ! advection terms.  The Coriolis source should be
                            ! zero, but is not due to truncation errors.  There
                            ! should be near-cancellation of the global integral
                            ! of this spurious Coriolis source.
    KE_adv     => NULL(),&  ! KE source from along-layer advection
    KE_visc    => NULL(),&  ! KE source from vertical viscosity
    KE_horvisc => NULL(),&  ! KE source from horizontal viscosity
    KE_dia     => NULL(),&  ! KE source from diapycnal diffusion
    diag_tmp3d => NULL()    ! 3D re-usable arrays for diagnostics

  ! diagnostic IDs
  integer :: id_e              = -1, id_e_D            = -1
  integer :: id_du_dt          = -1, id_dv_dt          = -1
  integer :: id_col_ht         = -1, id_dh_dt          = -1
  integer :: id_KE             = -1, id_dKEdt          = -1
  integer :: id_PE_to_KE       = -1, id_KE_Coradv      = -1
  integer :: id_KE_adv         = -1, id_KE_visc        = -1
  integer :: id_KE_horvisc     = -1, id_KE_dia         = -1
  integer :: id_uh_Rlay        = -1, id_vh_Rlay        = -1
  integer :: id_uhGM_Rlay      = -1, id_vhGM_Rlay      = -1
  integer :: id_h_Rlay         = -1, id_Rd1            = -1
  integer :: id_Rml            = -1, id_Rcv            = -1
  integer :: id_cg1            = -1, id_cfl_cg1        = -1
  integer :: id_cfl_cg1_x      = -1, id_cfl_cg1_y      = -1
  integer :: id_cg_ebt         = -1, id_Rd_ebt         = -1
  integer :: id_p_ebt          = -1
  integer :: id_temp_int       = -1, id_salt_int       = -1
  integer :: id_mass_wt        = -1, id_col_mass       = -1
  integer :: id_masscello      = -1, id_masso          = -1
  integer :: id_thetaoga       = -1, id_soga           = -1
  integer :: id_sosga          = -1, id_tosga          = -1
  integer :: id_temp_layer_ave = -1, id_salt_layer_ave = -1
  integer :: id_pbo            = -1
  integer :: id_thkcello       = -1, id_rhoinsitu      = -1
  integer :: id_rhopot0        = -1, id_rhopot2        = -1

  type(wave_speed_CS), pointer :: wave_speed_CSp => NULL()

  ! pointers used in calculation of time derivatives
  type(p3d) :: var_ptr(MAX_FIELDS_)
  type(p3d) :: deriv(MAX_FIELDS_)
  type(p3d) :: prev_val(MAX_FIELDS_)
  integer   :: nlay(MAX_FIELDS_)
  integer   :: num_time_deriv = 0

  ! for group halo pass
  type(group_pass_type) :: pass_KE_uv

end type diagnostics_CS

contains

subroutine calculate_diagnostic_fields(u, v, h, uh, vh, tv, ADp, CDp, fluxes, &
                                       dt, G, GV, CS, eta_bt)
  type(ocean_grid_type),                     intent(inout) :: G
  type(verticalGrid_type),                   intent(in)    :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: uh
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: vh
  type(thermo_var_ptrs),                     intent(in)    :: tv
  type(accel_diag_ptrs),                     intent(in)    :: ADp
  type(cont_diag_ptrs),                      intent(in)    :: CDp
  type(forcing),                             intent(in)    :: fluxes
  real,                                      intent(in)    :: dt
  type(diagnostics_CS),                      intent(inout) :: CS
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)    :: eta_bt

! Diagnostics not more naturally calculated elsewhere are computed here.

! Arguments:
!  (in)      u   - zonal velocity component (m/s)
!  (in)      v   - meridional velocity component (m/s)
!  (in)      h   - layer thickness,  meter(Bouss)  kg/m^2(non-Bouss)
!  (in)      uh  - transport through zonal faces = u*h*dy, m3/s(Bouss) kg/s(non-Bouss)
!  (in)      vh  - transport through meridional faces = v*h*dx, m3/s(Bouss) kg/s(non-Bouss)
!  (in)      tv  - structure pointing to various thermodynamic variables.
!  (in)      ADp - structure with pointers to accelerations in momentum equation
!  (in)      CDp - structure with pointers to terms in continuity equation
!  (in)      dt  - time difference in s since the last call to this subroutine
!  (in)      G   - ocean grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS  - control structure returned by a previous call to diagnostics_init
!  (in,opt)  eta_bt - An optional barotropic variable that gives the "correct"
!                     free surface height (Boussinesq) or total water column
!                     mass per unit area (non-Boussinesq).  This is used to
!                     dilate the layer thicknesses when calculating interface
!                     heights, in m or kg m-2.

  integer i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb

  ! coordinate variable potential density, in kg m-3.
  real :: Rcv(SZI_(G),SZJ_(G),SZK_(G))

  ! tmp array for surface properties
  real :: surface_field(SZI_(G),SZJ_(G))
  real :: pressure_1d(SZI_(G)) ! Temporary array for pressure when calling EOS
  real :: wt, wt_p

  ! squared Coriolis parameter at to h-points (1/s2)
  real :: f2_h

  ! magnitude of the gradient of f (1/(m*s))
  real :: mag_beta

  ! frequency squared used to avoid division by 0 (1/s2)
  ! value is roughly (pi / (the age of the universe) )^2.
  real, parameter :: absurdly_small_freq2 = 1e-34

  integer :: k_list

  real, dimension(SZK_(G)) :: temp_layer_ave, salt_layer_ave
  real :: thetaoga, soga, masso, tosga, sosga

  is  = G%isc  ; ie   = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq  = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  nz  = G%ke   ; nkmb = GV%nk_rho_varies

  ! smg: is the following robust to ALE? It seems a bit opaque.
  ! If the model is NOT in isopycnal mode then nkmb=0. But we need all the
  ! following diagnostics to treat all layers as variable density, so we set
  ! nkmb = nz, on the expectation that loops nkmb+1,nz will not iterate.
  ! This behavior is ANSI F77 but some compiler options can force at least
  ! one iteration that would break the following one-line workaround!
  if (nkmb==0) nkmb = nz

  if (loc(CS)==0) call MOM_error(FATAL, &
         "calculate_diagnostic_fields: Module must be initialized before used.")

  call calculate_derivs(dt, G, CS)

  if (ASSOCIATED(CS%e)) then
    call find_eta(h, tv, GV%g_Earth, G, GV, CS%e, eta_bt)
    if (CS%id_e > 0) call post_data(CS%id_e, CS%e, CS%diag)
  endif

  if (ASSOCIATED(CS%e_D)) then
    if (ASSOCIATED(CS%e)) then
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        CS%e_D(i,j,k) = CS%e(i,j,k) + G%bathyT(i,j)
      enddo ; enddo ; enddo
    else
      call find_eta(h, tv, GV%g_Earth, G, GV, CS%e_D, eta_bt)
      do k=1,nz+1 ; do j=js,je ; do i=is,ie
        CS%e_D(i,j,k) = CS%e_D(i,j,k) + G%bathyT(i,j)
      enddo ; enddo ; enddo
    endif

    if (CS%id_e_D > 0) call post_data(CS%id_e_D, CS%e_D, CS%diag)
  endif

  ! mass per area of grid cell (for Bouss, use Rho0)
  if (CS%id_masscello > 0) then
    do k=1,nz; do j=js,je ; do i=is,ie
       CS%diag_tmp3d(i,j,k) = GV%H_to_kg_m2*h(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_masscello, CS%diag_tmp3d, CS%diag)
  endif

  ! mass of liquid ocean (for Bouss, use Rho0)
  if (CS%id_masso > 0) then
    do k=1,nz; do j=js,je ; do i=is,ie
       CS%diag_tmp3d(i,j,k) = GV%H_to_kg_m2*h(i,j,k)*G%areaT(i,j)
    enddo ; enddo ; enddo
    masso = (reproducing_sum(sum(CS%diag_tmp3d,3)))
    call post_data(CS%id_masso, masso, CS%diag)
  endif

  ! diagnose thickness of grid cells (meter)
  if (CS%id_thkcello > 0) then

    ! thkcello = h for Boussinesq
    if (GV%Boussinesq) then
      call post_data(CS%id_thkcello, GV%H_to_m*h, CS%diag)

    ! thkcello = dp/(rho*g) for non-Boussinesq
    else
      do j=js,je

        ! pressure loading at top of surface layer (Pa)
        if(ASSOCIATED(fluxes%p_surf)) then
          do i=is,ie
            pressure_1d(i) = fluxes%p_surf(i,j)
          enddo
        else
          do i=is,ie
            pressure_1d(i) = 0.0
          enddo
        endif

        do k=1,nz
          ! pressure for EOS at the layer center (Pa)
          do i=is,ie
            pressure_1d(i) = pressure_1d(i) + 0.5*(GV%g_Earth*GV%H_to_kg_m2)*h(i,j,k)
          enddo
          ! store in-situ density (kg/m3) in diag_tmp3d
          call calculate_density(tv%T(:,j,k),tv%S(:,j,k), pressure_1d, &
                                 CS%diag_tmp3d(:,j,k), is, ie-is+1, tv%eqn_of_state)
          ! cell thickness = dz = dp/(g*rho) (meter); store in diag_tmp3d
          do i=is,ie
            CS%diag_tmp3d(i,j,k) = (GV%H_to_kg_m2*h(i,j,k))/CS%diag_tmp3d(i,j,k)
          enddo
          ! pressure for EOS at the bottom interface (Pa)
          do i=is,ie
            pressure_1d(i) = pressure_1d(i) + 0.5*(GV%g_Earth*GV%H_to_kg_m2)*h(i,j,k)
          enddo
        enddo ! k

      enddo ! j
      call post_data(CS%id_thkcello, CS%diag_tmp3d, CS%diag)
    endif
  endif

  ! volume mean potential temperature
  if (CS%id_thetaoga>0) then
    thetaoga = global_volume_mean(tv%T, h, G, GV)
    call post_data(CS%id_thetaoga, thetaoga, CS%diag)
  endif

  ! area mean SST
  if (CS%id_tosga > 0) then
    do j=js,je ; do i=is,ie
       surface_field(i,j) = tv%T(i,j,1)
    enddo ; enddo
    tosga = global_area_mean(surface_field, G)
    call post_data(CS%id_tosga, tosga, CS%diag)
  endif

  ! volume mean salinity
  if (CS%id_soga>0) then
    soga = global_volume_mean(tv%S, h, G, GV)
    call post_data(CS%id_soga, soga, CS%diag)
  endif

  ! area mean SSS
  if (CS%id_sosga > 0) then
    do j=js,je ; do i=is,ie
       surface_field(i,j) = tv%S(i,j,1)
    enddo ; enddo
    sosga = global_area_mean(surface_field, G)
    call post_data(CS%id_sosga, sosga, CS%diag)
  endif

  ! layer mean potential temperature
  if (CS%id_temp_layer_ave>0) then
    temp_layer_ave = global_layer_mean(tv%T, h, G, GV)
    call post_data_1d_k(CS%id_temp_layer_ave, temp_layer_ave, CS%diag)
  endif

  ! layer mean salinity
  if (CS%id_salt_layer_ave>0) then
    salt_layer_ave = global_layer_mean(tv%S, h, G, GV)
    call post_data_1d_k(CS%id_salt_layer_ave, salt_layer_ave, CS%diag)
  endif

  call calculate_vertical_integrals(h, tv, fluxes, G, GV, CS)

  if ((CS%id_Rml > 0) .or. (CS%id_Rcv > 0) .or. ASSOCIATED(CS%h_Rlay) .or. &
      ASSOCIATED(CS%uh_Rlay) .or. ASSOCIATED(CS%vh_Rlay) .or. &
      ASSOCIATED(CS%uhGM_Rlay) .or. ASSOCIATED(CS%vhGM_Rlay)) then

    if (associated(tv%eqn_of_state)) then
      pressure_1d(:) = tv%P_Ref
!$OMP parallel do default(none) shared(tv,Rcv,is,ie,js,je,nz,pressure_1d)
      do k=1,nz ; do j=js-1,je+1
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), pressure_1d, &
                               Rcv(:,j,k), is-1, ie-is+3, tv%eqn_of_state)
      enddo ; enddo
    else ! Rcv should not be used much in this case, so fill in sensible values.
      do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
        Rcv(i,j,k) = GV%Rlay(k)
      enddo ; enddo ; enddo
    endif
    if (CS%id_Rml > 0) call post_data(CS%id_Rml, Rcv, CS%diag)
    if (CS%id_Rcv > 0) call post_data(CS%id_Rcv, Rcv, CS%diag)

    if (ASSOCIATED(CS%h_Rlay)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(is,ie,js,je,nz,nkmb,CS,Rcv,h,GV) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do i=is,ie
          CS%h_Rlay(i,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          CS%h_Rlay(i,j,k) = h(i,j,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, Rcv(i,j,k), k_list, nz, wt, wt_p)
          CS%h_Rlay(i,j,k_list)   = CS%h_Rlay(i,j,k_list)   + h(i,j,k)*wt
          CS%h_Rlay(i,j,k_list+1) = CS%h_Rlay(i,j,k_list+1) + h(i,j,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_h_Rlay > 0) call post_data(CS%id_h_Rlay, CS%h_Rlay, CS%diag)
    endif

    if (ASSOCIATED(CS%uh_Rlay)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(Isq,Ieq,js,je,nz,nkmb,Rcv,CS,GV,uh) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do I=Isq,Ieq
          CS%uh_Rlay(I,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do I=Isq,Ieq
          CS%uh_Rlay(I,j,k) = uh(I,j,k)
        enddo ; enddo
        k_list = nz/2
        do k=1,nkmb ; do I=Isq,Ieq
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i+1,j,k)), k_list, nz, wt, wt_p)
          CS%uh_Rlay(I,j,k_list)   = CS%uh_Rlay(I,j,k_list)   + uh(I,j,k)*wt
          CS%uh_Rlay(I,j,k_list+1) = CS%uh_Rlay(I,j,k_list+1) + uh(I,j,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_uh_Rlay > 0) call post_data(CS%id_uh_Rlay, CS%uh_Rlay, CS%diag)
    endif

    if (ASSOCIATED(CS%vh_Rlay)) then
      k_list = nz/2
!$OMP parallel do default(none)  shared(Jsq,Jeq,is,ie,nz,nkmb,Rcv,CS,GV,vh) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do J=Jsq,Jeq
        do k=1,nkmb ; do i=is,ie
          CS%vh_Rlay(i,J,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          CS%vh_Rlay(i,J,k) = vh(i,J,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i,j+1,k)), k_list, nz, wt, wt_p)
          CS%vh_Rlay(i,J,k_list)   = CS%vh_Rlay(i,J,k_list)   + vh(i,J,k)*wt
          CS%vh_Rlay(i,J,k_list+1) = CS%vh_Rlay(i,J,k_list+1) + vh(i,J,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_vh_Rlay > 0) call post_data(CS%id_vh_Rlay, CS%vh_Rlay, CS%diag)
    endif

    if (ASSOCIATED(CS%uhGM_Rlay) .and. ASSOCIATED(CDp%uhGM)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(Isq,Ieq,js,je,nz,nkmb,Rcv,CDP,CS,GV) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do j=js,je
        do k=1,nkmb ; do I=Isq,Ieq
          CS%uhGM_Rlay(I,j,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do I=Isq,Ieq
          CS%uhGM_Rlay(I,j,k) = CDp%uhGM(I,j,k)
        enddo ; enddo
        do k=1,nkmb ; do I=Isq,Ieq
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i+1,j,k)), k_list, nz, wt, wt_p)
          CS%uhGM_Rlay(I,j,k_list)   = CS%uhGM_Rlay(I,j,k_list)   + CDp%uhGM(I,j,k)*wt
          CS%uhGM_Rlay(I,j,k_list+1) = CS%uhGM_Rlay(I,j,k_list+1) + CDp%uhGM(I,j,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_uh_Rlay > 0) call post_data(CS%id_uhGM_Rlay, CS%uhGM_Rlay, CS%diag)
    endif

    if (ASSOCIATED(CS%vhGM_Rlay) .and. ASSOCIATED(CDp%vhGM)) then
      k_list = nz/2
!$OMP parallel do default(none) shared(is,ie,Jsq,Jeq,nz,nkmb,CS,CDp,Rcv,GV) &
!$OMP                          private(wt,wt_p) firstprivate(k_list)
      do J=Jsq,Jeq
        do k=1,nkmb ; do i=is,ie
          CS%vhGM_Rlay(i,J,k) = 0.0
        enddo ; enddo
        do k=nkmb+1,nz ; do i=is,ie
          CS%vhGM_Rlay(i,J,k) = CDp%vhGM(i,J,k)
        enddo ; enddo
        do k=1,nkmb ; do i=is,ie
          call find_weights(GV%Rlay, 0.5*(Rcv(i,j,k)+Rcv(i,j+1,k)), k_list, nz, wt, wt_p)
          CS%vhGM_Rlay(i,J,k_list)   = CS%vhGM_Rlay(i,J,k_list)   + CDp%vhGM(i,J,k)*wt
          CS%vhGM_Rlay(i,J,k_list+1) = CS%vhGM_Rlay(i,J,k_list+1) + CDp%vhGM(i,J,k)*wt_p
        enddo ; enddo
      enddo

      if (CS%id_vhGM_Rlay > 0) call post_data(CS%id_vhGM_Rlay, CS%vhGM_Rlay, CS%diag)
    endif
  endif

  if (associated(tv%eqn_of_state)) then
    if (CS%id_rhopot0 > 0) then
      pressure_1d(:) = 0.
!$OMP parallel do default(none) shared(tv,Rcv,is,ie,js,je,nz,pressure_1d)
      do k=1,nz ; do j=js,je
        call calculate_density(tv%T(:,j,k),tv%S(:,j,k),pressure_1d, &
                               Rcv(:,j,k),is,ie-is+1, tv%eqn_of_state)
      enddo ; enddo
      if (CS%id_rhopot0 > 0) call post_data(CS%id_rhopot0, Rcv, CS%diag)
    endif
    if (CS%id_rhopot2 > 0) then
      pressure_1d(:) = 2.E7 ! 2000 dbars
!$OMP parallel do default(none) shared(tv,Rcv,is,ie,js,je,nz,pressure_1d)
      do k=1,nz ; do j=js,je
        call calculate_density(tv%T(:,j,k),tv%S(:,j,k),pressure_1d, &
                               Rcv(:,j,k),is,ie-is+1, tv%eqn_of_state)
      enddo ; enddo
      if (CS%id_rhopot2 > 0) call post_data(CS%id_rhopot2, Rcv, CS%diag)
    endif
    if (CS%id_rhoinsitu > 0) then
!$OMP parallel do default(none) shared(tv,Rcv,is,ie,js,je,nz,pressure_1d,h,GV)
      do j=js,je
        pressure_1d(:) = 0. ! Start at p=0 Pa at surface
        do k=1,nz
          pressure_1d(:) =  pressure_1d(:) + 0.5 * h(:,j,k) * GV%H_to_Pa ! Pressure in middle of layer k
          call calculate_density(tv%T(:,j,k),tv%S(:,j,k),pressure_1d, &
                                 Rcv(:,j,k),is,ie-is+1, tv%eqn_of_state)
          pressure_1d(:) =  pressure_1d(:) + 0.5 * h(:,j,k) * GV%H_to_Pa ! Pressure at bottom of layer k
        enddo
      enddo
      if (CS%id_rhoinsitu > 0) call post_data(CS%id_rhoinsitu, Rcv, CS%diag)
    endif
  endif

  if ((CS%id_cg1>0) .or. (CS%id_Rd1>0) .or. (CS%id_cfl_cg1>0) .or. &
      (CS%id_cfl_cg1_x>0) .or. (CS%id_cfl_cg1_y>0)) then
    call wave_speed(h, tv, G, GV, CS%cg1, CS%wave_speed_CSp)
    if (CS%id_cg1>0) call post_data(CS%id_cg1, CS%cg1, CS%diag)
    if (CS%id_Rd1>0) then
!$OMP parallel do default(none) shared(is,ie,js,je,G,CS) &
!$OMP                          private(f2_h,mag_beta)
      do j=js,je ; do i=is,ie
        ! Blend the equatorial deformation radius with the standard one.
        f2_h = absurdly_small_freq2 + 0.25 * &
            ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
             (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2))
        mag_beta = sqrt(0.5 * ( &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2) ))
        CS%Rd1(i,j) = CS%cg1(i,j) / sqrt(f2_h + CS%cg1(i,j) * mag_beta)

      enddo ; enddo
      call post_data(CS%id_Rd1, CS%Rd1, CS%diag)
    endif
    if (CS%id_cfl_cg1>0) then
      do j=js,je ; do i=is,ie
        CS%cfl_cg1(i,j) = (dt*CS%cg1(i,j)) * (G%IdxT(i,j) + G%IdyT(i,j))
      enddo ; enddo
      call post_data(CS%id_cfl_cg1, CS%cfl_cg1, CS%diag)
    endif
    if (CS%id_cfl_cg1_x>0) then
      do j=js,je ; do i=is,ie
        CS%cfl_cg1_x(i,j) = (dt*CS%cg1(i,j)) * G%IdxT(i,j)
      enddo ; enddo
      call post_data(CS%id_cfl_cg1_x, CS%cfl_cg1_x, CS%diag)
    endif
    if (CS%id_cfl_cg1_y>0) then
      do j=js,je ; do i=is,ie
        CS%cfl_cg1_y(i,j) = (dt*CS%cg1(i,j)) * G%IdyT(i,j)
      enddo ; enddo
      call post_data(CS%id_cfl_cg1_y, CS%cfl_cg1_y, CS%diag)
    endif
  endif
  if ((CS%id_cg_ebt>0) .or. (CS%id_Rd_ebt>0) .or. (CS%id_p_ebt>0)) then
    if (CS%id_p_ebt>0) then
      call wave_speed(h, tv, G, GV, CS%cg1, CS%wave_speed_CSp, use_ebt_mode=.true., &
                      mono_N2_column_fraction=CS%mono_N2_column_fraction, &
                      mono_N2_depth=CS%mono_N2_depth, modal_structure=CS%p_ebt)
      call post_data(CS%id_p_ebt, CS%p_ebt, CS%diag)
    else
      call wave_speed(h, tv, G, GV, CS%cg1, CS%wave_speed_CSp, use_ebt_mode=.true., &
                      mono_N2_column_fraction=CS%mono_N2_column_fraction, &
                      mono_N2_depth=CS%mono_N2_depth)
    endif
    if (CS%id_cg_ebt>0) call post_data(CS%id_cg_ebt, CS%cg1, CS%diag)
    if (CS%id_Rd_ebt>0) then
!$OMP parallel do default(none) shared(is,ie,js,je,G,CS) &
!$OMP                          private(f2_h,mag_beta)
      do j=js,je ; do i=is,ie
        ! Blend the equatorial deformation radius with the standard one.
        f2_h = absurdly_small_freq2 + 0.25 * &
            ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
             (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2))
        mag_beta = sqrt(0.5 * ( &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I-1,J)) * G%IdxCv(i,J))**2 + &
             ((G%CoriolisBu(I,J-1)-G%CoriolisBu(I-1,J-1)) * G%IdxCv(i,J-1))**2) + &
            (((G%CoriolisBu(I,J)-G%CoriolisBu(I,J-1)) * G%IdyCu(I,j))**2 + &
             ((G%CoriolisBu(I-1,J)-G%CoriolisBu(I-1,J-1)) * G%IdyCu(I-1,j))**2) ))
        CS%Rd1(i,j) = CS%cg1(i,j) / sqrt(f2_h + CS%cg1(i,j) * mag_beta)

      enddo ; enddo
      call post_data(CS%id_Rd_ebt, CS%Rd1, CS%diag)
    endif
  endif

  if (dt > 0.0) then
    if (CS%id_du_dt>0) call post_data(CS%id_du_dt, CS%du_dt, CS%diag)

    if (CS%id_dv_dt>0) call post_data(CS%id_dv_dt, CS%dv_dt, CS%diag)

    if (CS%id_dh_dt>0) call post_data(CS%id_dh_dt, CS%dh_dt, CS%diag)

    call calculate_energy_diagnostics(u, v, h, uh, vh, ADp, CDp, G, CS)
  endif

end subroutine calculate_diagnostic_fields


subroutine find_weights(Rlist, R_in, k, nz, wt, wt_p)
  real,     intent(in)    :: Rlist(:), R_in
  integer,  intent(inout) :: k
  integer,  intent(in)    :: nz
  real,     intent(out)   :: wt, wt_p

  ! This subroutine finds location of R_in in an increasing ordered
  ! list, Rlist, returning as k the element such that
  !  Rlist(k) <= R_in < Rlist(k+1), and where wt and wt_p are the linear
  ! weights that should be assigned to elements k and k+1.

  integer :: k_upper, k_lower, k_new, inc

  ! First, bracket the desired point.
  if ((k < 1) .or. (k > nz)) k = nz/2

  k_upper = k ; k_lower = k ; inc = 1
  if (R_in < Rlist(k)) then
    do
      k_lower = max(k_lower-inc, 1)
      if ((k_lower == 1) .or. (R_in >= Rlist(k_lower))) exit
      k_upper = k_lower
      inc = inc*2
    end do
  else
    do
      k_upper = min(k_upper+inc, nz)
      if ((k_upper == nz) .or. (R_in < Rlist(k_upper))) exit
      k_lower = k_upper
      inc = inc*2
    end do
  endif

  if ((k_lower == 1) .and. (R_in <= Rlist(k_lower))) then
    k = 1 ; wt = 1.0 ; wt_p = 0.0
  else if ((k_upper == nz) .and. (R_in >= Rlist(k_upper))) then
    k = nz-1 ; wt = 0.0 ; wt_p = 1.0
  else
    do
      if (k_upper <= k_lower+1) exit
      k_new = (k_upper + k_lower) / 2
      if (R_in < Rlist(k_new)) then
        k_upper = k_new
      else
        k_lower = k_new
      endif
    end do

!   Uncomment this as a code check
!    if ((R_in < Rlist(k_lower)) .or. (R_in >= Rlist(k_upper)) .or. (k_upper-k_lower /= 1)) &
!      write (*,*) "Error: ",R_in," is not between R(",k_lower,") = ", &
!        Rlist(k_lower)," and R(",k_upper,") = ",Rlist(k_upper),"."
    k = k_lower
    wt = (Rlist(k_upper) - R_in) / (Rlist(k_upper) - Rlist(k_lower))
    wt_p = 1.0 - wt

  endif

end subroutine find_weights

subroutine calculate_vertical_integrals(h, tv, fluxes, G, GV, CS)
  type(ocean_grid_type),                    intent(inout) :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  type(thermo_var_ptrs),                    intent(in)    :: tv
  type(forcing),                            intent(in)    :: fluxes
  type(diagnostics_CS),                     intent(inout) :: CS

! Subroutine calculates vertical integrals of several tracers, along
! with the mass-weight of these tracers, the total column mass, and the
! carefully calculated column height.

! Arguments:
!  (in)      h  - layer thickness: metre (Bouss) or kg/ m2 (non-Bouss)
!  (in)      tv - structure pointing to thermodynamic variables
!  (in)      fluxes - a structure containing the surface fluxes.
!  (in)      G  - ocean grid structure
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - control structure returned by a previous call to diagnostics_init

  real, dimension(SZI_(G), SZJ_(G)) :: &
    z_top, &  ! Height of the top of a layer or the ocean, in m.
    z_bot, &  ! Height of the bottom of a layer (for id_mass) or the
              ! (positive) depth of the ocean (for id_col_ht), in m.
    mass, &   ! integrated mass of the water column, in kg m-2.  For
              ! non-Boussinesq models this is rho*dz. For Boussiensq
              ! models, this is either the integral of in-situ density
              ! (rho*dz for col_mass) or reference dens (Rho_0*dz for mass_wt).
    btm_pres,&! The pressure at the ocean bottom, or CMIP variable 'pbo'.
              ! This is the column mass multiplied by gravity plus the pressure
              ! at the ocean surface.
    dpress, &    ! Change in hydrostatic pressure across a layer, in Pa.
    tr_int    ! vertical integral of a tracer times density,
              ! (Rho_0 in a Boussinesq model) in TR kg m-2.
  real    :: IG_Earth  ! Inverse of gravitational acceleration, in s2 m-1.

  integer :: i, j, k, is, ie, js, je, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (CS%id_mass_wt > 0) then
    do j=js,je ; do i=is,ie ; mass(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      mass(i,j) = mass(i,j) + GV%H_to_kg_m2*h(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_mass_wt, mass, CS%diag)
  endif

  if (CS%id_temp_int > 0) then
    do j=js,je ; do i=is,ie ; tr_int(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      tr_int(i,j) = tr_int(i,j) + (GV%H_to_kg_m2*h(i,j,k))*tv%T(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_temp_int, tr_int, CS%diag)
  endif

  if (CS%id_salt_int > 0) then
    do j=js,je ; do i=is,ie ; tr_int(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js,je ; do i=is,ie
      tr_int(i,j) = tr_int(i,j) + (GV%H_to_kg_m2*h(i,j,k))*tv%S(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_salt_int, tr_int, CS%diag)
  endif

  if (CS%id_col_ht > 0) then
    call find_eta(h, tv, GV%g_Earth, G, GV, z_top)
    do j=js,je ; do i=is,ie
      z_bot(i,j) = z_top(i,j) + G%bathyT(i,j)
    enddo ; enddo
    call post_data(CS%id_col_ht, z_bot, CS%diag)
  endif

  if (CS%id_col_mass > 0 .or. CS%id_pbo > 0) then
    do j=js,je ; do i=is,ie ; mass(i,j) = 0.0 ; enddo ; enddo
    if (GV%Boussinesq) then
      if (associated(tv%eqn_of_state)) then
        IG_Earth = 1.0 / GV%g_Earth
!       do j=js,je ; do i=is,ie ; z_bot(i,j) = -P_SURF(i,j)/GV%H_to_Pa ; enddo ; enddo
        do j=js,je ; do i=is,ie ; z_bot(i,j) = 0.0 ; enddo ; enddo
        do k=1,nz
          do j=js,je ; do i=is,ie
            z_top(i,j) = z_bot(i,j)
            z_bot(i,j) = z_top(i,j) - GV%H_to_m*h(i,j,k)
          enddo ; enddo
          call int_density_dz(tv%T(:,:,k), tv%S(:,:,k), &
                              z_top, z_bot, 0.0, GV%H_to_kg_m2, GV%g_Earth, &
                              G%HI, G%HI, tv%eqn_of_state, dpress)
          do j=js,je ; do i=is,ie
            mass(i,j) = mass(i,j) + dpress(i,j) * IG_Earth
          enddo ; enddo
        enddo
      else
        do k=1,nz ; do j=js,je ; do i=is,ie
          mass(i,j) = mass(i,j) + (GV%H_to_m*GV%Rlay(k))*h(i,j,k)
        enddo ; enddo ; enddo
      endif
    else
      do k=1,nz ; do j=js,je ; do i=is,ie
        mass(i,j) = mass(i,j) + GV%H_to_kg_m2*h(i,j,k)
      enddo ; enddo ; enddo
    endif
    if (CS%id_col_mass > 0) then
      call post_data(CS%id_col_mass, mass, CS%diag)
    endif
    if (CS%id_pbo > 0) then
      do j=js,je ; do i=is,ie ; btm_pres(i,j) = 0.0 ; enddo ; enddo
      ! 'pbo' is defined as the sea water pressure at the sea floor
      !     pbo = (mass * g) + pso
      ! where pso is the sea water pressure at sea water surface
      ! note that pso is equivalent to fluxes%p_surf
      do j=js,je ; do i=is,ie
        btm_pres(i,j) = mass(i,j) * GV%g_Earth
        if (ASSOCIATED(fluxes%p_surf)) then
          btm_pres(i,j) = btm_pres(i,j) + fluxes%p_surf(i,j)
        endif
      enddo ; enddo
      call post_data(CS%id_pbo, btm_pres, CS%diag)
    endif
  endif

end subroutine calculate_vertical_integrals


subroutine calculate_energy_diagnostics(u, v, h, uh, vh, ADp, CDp, G, CS)
  type(ocean_grid_type),                     intent(inout) :: G
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)    :: uh
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)    :: vh
  type(accel_diag_ptrs),                     intent(in)    :: ADp
  type(cont_diag_ptrs),                      intent(in)    :: CDp
  type(diagnostics_CS),                      intent(inout) :: CS

! This subroutine calculates terms in the mechanical energy budget.

! Arguments:
!  (in)      u   - zonal velocity component (m/s)
!  (in)      v   - meridional velocity componnent (m/s)
!  (in)      h   - layer thickness: metre(Bouss) of kg/m2(non-Bouss)
!  (in)      uh  - transport through zonal faces=u*h*dy: m3/s (Bouss) kg/s(non-Bouss)
!  (in)      vh  - transport through merid faces=v*h*dx: m3/s (Bouss) kg/s(non-Bouss)
!  (in)      ADp - structure pointing to accelerations in momentum equation
!  (in)      CDp - structure pointing to terms in continuity equations
!  (in)      G   - ocean grid structure
!  (in)      CS  - control structure returned by a previous call to diagnostics_init

  real :: KE_u(SZIB_(G),SZJ_(G))
  real :: KE_v(SZI_(G),SZJB_(G))
  real :: KE_h(SZI_(G),SZJ_(G))

  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  do j=js-1,je ; do i=is-1,ie
    KE_u(I,j) = 0.0 ; KE_v(i,J) = 0.0
  enddo ; enddo

  if (ASSOCIATED(CS%KE)) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      CS%KE(i,j,k) = ((u(I,j,k)*u(I,j,k) + u(I-1,j,k)*u(I-1,j,k)) + &
          (v(i,J,k)*v(i,J,k) + v(i,J-1,k)*v(i,J-1,k)))*0.25
      ! DELETE THE FOLLOWING...  Make this 0 to test the momentum balance,
      ! or a huge number to test the continuity balance.
      ! CS%KE(i,j,k) *= 1e20
    enddo ; enddo ; enddo
    if (CS%id_KE > 0) call post_data(CS%id_KE, CS%KE, CS%diag)
  endif

  if(.not.G%symmetric) then
    if(ASSOCIATED(CS%dKE_dt) .OR. ASSOCIATED(CS%PE_to_KE) .OR. ASSOCIATED(CS%KE_CorAdv) .OR. &
       ASSOCIATED(CS%KE_adv) .OR. ASSOCIATED(CS%KE_visc)  .OR. ASSOCIATED(CS%KE_horvisc).OR. &
       ASSOCIATED(CS%KE_dia) ) then
        call create_group_pass(CS%pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)
    endif
  endif

  if (ASSOCIATED(CS%dKE_dt)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*CS%du_dt(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*CS%dv_dt(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = CS%KE(i,j,k)*CS%dh_dt(i,j,k)
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%dKE_dt(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_dKEdt > 0) call post_data(CS%id_dKEdt, CS%dKE_dt, CS%diag)
  endif

  if (ASSOCIATED(CS%PE_to_KE)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%PFu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%PFv(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%PE_to_KE(i,j,k) = 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_PE_to_KE > 0) call post_data(CS%id_PE_to_KE, CS%PE_to_KE, CS%diag)
  endif

  if (ASSOCIATED(CS%KE_CorAdv)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%CAu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%CAv(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = -CS%KE(i,j,k) * G%IareaT(i,j) * &
            (uh(I,j,k) - uh(I-1,j,k) + vh(i,J,k) - vh(i,J-1,k))
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_CorAdv(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_Coradv > 0) call post_data(CS%id_KE_Coradv, CS%KE_Coradv, CS%diag)
  endif

  if (ASSOCIATED(CS%KE_adv)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%gradKEu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%gradKEv(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = -CS%KE(i,j,k) * G%IareaT(i,j) * &
            (uh(I,j,k) - uh(I-1,j,k) + vh(i,J,k) - vh(i,J-1,k))
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_adv(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_adv > 0) call post_data(CS%id_KE_adv, CS%KE_adv, CS%diag)
  endif

  if (ASSOCIATED(CS%KE_visc)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%du_dt_visc(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%dv_dt_visc(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_visc(i,j,k) = 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_visc > 0) call post_data(CS%id_KE_visc, CS%KE_visc, CS%diag)
  endif

  if (ASSOCIATED(CS%KE_horvisc)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%diffu(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%diffv(i,J,k)
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_horvisc(i,j,k) = 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_horvisc > 0) call post_data(CS%id_KE_horvisc, CS%KE_horvisc, CS%diag)
  endif

  if (ASSOCIATED(CS%KE_dia)) then
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        KE_u(I,j) = uh(I,j,k)*G%dxCu(I,j)*ADp%du_dt_dia(I,j,k)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        KE_v(i,J) = vh(i,J,k)*G%dyCv(i,J)*ADp%dv_dt_dia(i,J,k)
      enddo ; enddo
      do j=js,je ; do i=is,ie
        KE_h(i,j) = CS%KE(i,j,k) * &
            (CDp%diapyc_vel(i,j,k) - CDp%diapyc_vel(i,j,k+1))
      enddo ; enddo
      if (.not.G%symmetric) &
         call do_group_pass(CS%pass_KE_uv, G%domain)
      do j=js,je ; do i=is,ie
        CS%KE_dia(i,j,k) = KE_h(i,j) + 0.5 * G%IareaT(i,j) * &
            (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      enddo ; enddo
    enddo
    if (CS%id_KE_dia > 0) call post_data(CS%id_KE_dia, CS%KE_dia, CS%diag)
  endif

end subroutine calculate_energy_diagnostics


subroutine register_time_deriv(f_ptr, deriv_ptr, CS)
  real, dimension(:,:,:), target :: f_ptr, deriv_ptr
  type(diagnostics_CS),  pointer :: CS

! This subroutine registers fields to calculate a diagnostic time derivative.
! Arguments:
!  (target)  f_ptr     - field whose derivative is taken
!  (in)      deriv_ptr - field in which the calculated time derivatives placed
!  (in)      num_lay   - number of layers in this field
!  (in)      CS        - control structure returned by previous call to diagnostics_init

  integer :: m

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "register_time_deriv: Module must be initialized before it is used.")

  if (CS%num_time_deriv >= MAX_FIELDS_) then
    call MOM_error(WARNING,"MOM_diagnostics:  Attempted to register more than " // &
                   "MAX_FIELDS_ diagnostic time derivatives via register_time_deriv.")
    return
  endif

  m = CS%num_time_deriv+1 ; CS%num_time_deriv = m

  CS%nlay(m) = size(f_ptr(:,:,:),3)
  CS%deriv(m)%p => deriv_ptr
  allocate(CS%prev_val(m)%p(size(f_ptr(:,:,:),1), size(f_ptr(:,:,:),2), CS%nlay(m)) )

  CS%var_ptr(m)%p => f_ptr
  CS%prev_val(m)%p(:,:,:) = f_ptr(:,:,:)

end subroutine register_time_deriv


subroutine calculate_derivs(dt, G, CS)
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(inout) :: G
  type(diagnostics_CS),  intent(inout) :: CS

! This subroutine calculates all registered time derivatives.
! Arguments:
!  (in)      dt - time interval in s over which differences occur
!  (in)      G  - ocean grid structure.
!  (in)      CS - control structure returned by previous call to diagnostics_init

  integer i, j, k, m
  real Idt

  if (dt > 0.0) then ; Idt = 1.0/dt
  else ; return ; endif

  do m=1,CS%num_time_deriv
    do k=1,CS%nlay(m) ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      CS%deriv(m)%p(i,j,k) = (CS%var_ptr(m)%p(i,j,k) - CS%prev_val(m)%p(i,j,k)) * Idt
      CS%prev_val(m)%p(i,j,k) = CS%var_ptr(m)%p(i,j,k)
    enddo ; enddo ; enddo
  enddo

end subroutine calculate_derivs

subroutine MOM_diagnostics_init(MIS, ADp, CDp, Time, G, GV, param_file, diag, CS)
  type(ocean_internal_state), intent(in)    :: MIS
  type(accel_diag_ptrs),      intent(inout) :: ADp
  type(cont_diag_ptrs),       intent(inout) :: CDp
  type(time_type),            intent(in)    :: Time
  type(ocean_grid_type),      intent(in)    :: G
  type(verticalGrid_type),    intent(in)    :: GV
  type(param_file_type),      intent(in)    :: param_file
  type(diag_ctrl), target,    intent(inout) :: diag
  type(diagnostics_CS),       pointer       :: CS

! Arguments
!  (in)     MIS    - For "MOM Internal State" a set of pointers to the fields and
!                    accelerations that make up the ocean's internal physical
!                    state.
!  (inout)  ADp    - structure with pointers to momentum equation terms
!  (inout)  CDp    - structure with pointers to continuity equation terms
!  (in)     Time   - current model time
!  (in)     G      - ocean grid structure
!  (in)     GV     - The ocean's vertical grid structure.
!  (in) param_file - structure indicating the open file to parse for
!                     model parameter values
!  (in)     diag   - structure to regulate diagnostic output
!  (in/out) CS     - pointer set to point to control structure for this module

! This include declares and sets the variable "version".
#include "version_variable.h"

  character(len=40)  :: mod = "MOM_diagnostics" ! This module's name.
  real :: omega, f2_min
  character(len=48) :: thickness_units, flux_units
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz, nkml, nkbl
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j

  is   = G%isc  ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec
  Isq  = G%IscB ; Ieq  = G%IecB ; Jsq  = G%JscB ; Jeq  = G%JecB
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "MOM_diagnostics_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "DIAG_EBT_MONO_N2_COLUMN_FRACTION", CS%mono_N2_column_fraction, &
                 "The lower fraction of water column over which N2 is limited as monotonic\n"// &
                 "for the purposes of calculating the equivalent barotropic wave speed.", &
                 units='nondim', default=0.)
  call get_param(param_file, mod, "DIAG_EBT_MONO_N2_DEPTH", CS%mono_N2_depth, &
                 "The depth below which N2 is limited as monotonic for the\n"// &
                 "purposes of calculating the equivalent barotropic wave speed.", &
                 units='m', default=-1.)

  if (GV%Boussinesq) then
    thickness_units = "meter" ; flux_units = "meter3 second-1"
  else
    thickness_units = "kilogram meter-2" ; flux_units = "kilogram second-1"
  endif

  CS%id_temp_layer_ave = register_diag_field('ocean_model', 'temp_layer_ave', diag%axesZL, Time, &
      'Layer Average Ocean Temperature', 'Celsius')

  CS%id_salt_layer_ave = register_diag_field('ocean_model', 'salt_layer_ave', diag%axesZL, Time, &
      'Layer Average Ocean Salinity', 'ppt')

  CS%id_masscello = register_diag_field('ocean_model', 'masscello', diag%axesTL,&
      Time, 'Mass per unit area of liquid ocean grid cell', 'kg m-2',           &
      standard_name='sea_water_mass_per_unit_area', v_extensive=.true.)

  CS%id_masso = register_scalar_field('ocean_model', 'masso', Time,  &
      diag, 'Mass of liquid ocean', 'kg', standard_name='sea_water_mass')

  CS%id_thkcello = register_diag_field('ocean_model', 'thkcello', diag%axesTL, Time, &
      long_name = 'Cell Thickness', standard_name='cell_thickness', units='m', v_extensive=.true.)

  if (((CS%id_masscello>0) .or. (CS%id_masso>0) .or. (CS%id_thkcello>0.and..not.GV%Boussinesq)) &
      .and. .not.ASSOCIATED(CS%diag_tmp3d)) then
    call safe_alloc_ptr(CS%diag_tmp3d,isd,ied,jsd,jed,nz)
  endif

  CS%id_thetaoga = register_scalar_field('ocean_model', 'thetaoga',    &
      Time, diag, 'Global Mean Ocean Potential Temperature', 'Celsius',&
      standard_name='sea_water_potential_temperature')

  CS%id_soga = register_scalar_field('ocean_model', 'soga', &
      Time, diag, 'Global Mean Ocean Salinity', 'ppt',      &
      standard_name='sea_water_salinity')

  CS%id_tosga = register_scalar_field('ocean_model', 'sst_global', Time, diag,&
      long_name='Global Area Average Sea Surface Temperature',                &
      units='degC', standard_name='sea_surface_temperature',                  &
      cmor_field_name='tosga', cmor_standard_name='sea_surface_temperature',  &
      cmor_units='degC', cmor_long_name='Sea Surface Temperature')

  CS%id_sosga = register_scalar_field('ocean_model', 'sss_global', Time, diag,&
      long_name='Global Area Average Sea Surface Salinity',                   &
      units='ppt', standard_name='sea_surface_salinity',                      &
      cmor_field_name='sosga', cmor_standard_name='sea_surface_salinity',     &
      cmor_units='ppt', cmor_long_name='Sea Surface Salinity')

  CS%id_e = register_diag_field('ocean_model', 'e', diag%axesTi, Time, &
      'Interface Height Relative to Mean Sea Level', 'meter')
  if (CS%id_e>0) call safe_alloc_ptr(CS%e,isd,ied,jsd,jed,nz+1)

  CS%id_e_D = register_diag_field('ocean_model', 'e_D', diag%axesTi, Time, &
      'Interface Height above the Seafloor', 'meter')
  if (CS%id_e_D>0) call safe_alloc_ptr(CS%e_D,isd,ied,jsd,jed,nz+1)

  CS%id_Rml = register_diag_field('ocean_model', 'Rml', diag%axesTL, Time, &
      'Mixed Layer Coordinate Potential Density', 'kg meter-3')

  CS%id_Rcv = register_diag_field('ocean_model', 'Rho_cv', diag%axesTL, Time, &
      'Coordinate Potential Density', 'kg meter-3')

  CS%id_rhopot0 = register_diag_field('ocean_model', 'rhopot0', diag%axesTL, Time, &
      'Potential density referenced to surface', 'kg meter-3')
  CS%id_rhopot2 = register_diag_field('ocean_model', 'rhopot2', diag%axesTL, Time, &
      'Potential density referenced to 2000 dbar', 'kg meter-3')
  CS%id_rhoinsitu = register_diag_field('ocean_model', 'rhoinsitu', diag%axesTL, Time, &
      'In situ density', 'kg meter-3')

  CS%id_du_dt = register_diag_field('ocean_model', 'dudt', diag%axesCuL, Time, &
      'Zonal Acceleration', 'meter second-2')
  if ((CS%id_du_dt>0) .and. .not.ASSOCIATED(CS%du_dt)) then
    call safe_alloc_ptr(CS%du_dt,IsdB,IedB,jsd,jed,nz)
    call register_time_deriv(MIS%u, CS%du_dt, CS)
  endif

  CS%id_dv_dt = register_diag_field('ocean_model', 'dvdt', diag%axesCvL, Time, &
      'Meridional Acceleration', 'meter second-2')
  if ((CS%id_dv_dt>0) .and. .not.ASSOCIATED(CS%dv_dt)) then
    call safe_alloc_ptr(CS%dv_dt,isd,ied,JsdB,JedB,nz)
    call register_time_deriv(MIS%v, CS%dv_dt, CS)
  endif

  CS%id_dh_dt = register_diag_field('ocean_model', 'dhdt', diag%axesTL, Time, &
      'Thickness tendency', trim(thickness_units)//" second-1")
  if ((CS%id_dh_dt>0) .and. .not.ASSOCIATED(CS%dh_dt)) then
    call safe_alloc_ptr(CS%dh_dt,isd,ied,jsd,jed,nz)
    call register_time_deriv(MIS%h, CS%dh_dt, CS)
  endif

  ! layer thickness variables
  !if (GV%nk_rho_varies > 0) then
    CS%id_h_Rlay = register_diag_field('ocean_model', 'h_rho', diag%axesTL, Time, &
        'Layer thicknesses in pure potential density coordinates', thickness_units)
    if (CS%id_h_Rlay>0) call safe_alloc_ptr(CS%h_Rlay,isd,ied,jsd,jed,nz)

    CS%id_uh_Rlay = register_diag_field('ocean_model', 'uh_rho', diag%axesCuL, Time, &
        'Zonal volume transport in pure potential density coordinates', flux_units)
    if (CS%id_uh_Rlay>0) call safe_alloc_ptr(CS%uh_Rlay,IsdB,IedB,jsd,jed,nz)

    CS%id_vh_Rlay = register_diag_field('ocean_model', 'vh_rho', diag%axesCvL, Time, &
        'Meridional volume transport in pure potential density coordinates', flux_units)
    if (CS%id_vh_Rlay>0) call safe_alloc_ptr(CS%vh_Rlay,isd,ied,JsdB,JedB,nz)

    CS%id_uhGM_Rlay = register_diag_field('ocean_model', 'uhGM_rho', diag%axesCuL, Time, &
        'Zonal volume transport due to interface height diffusion in pure potential &
        &density coordinates', flux_units)
    if (CS%id_uhGM_Rlay>0) call safe_alloc_ptr(CS%uhGM_Rlay,IsdB,IedB,jsd,jed,nz)

    CS%id_vhGM_Rlay = register_diag_field('ocean_model', 'vhGM_rho', diag%axesCvL, Time, &
        'Meridional volume transport due to interface height diffusion in pure &
        &potential density coordinates', flux_units)
    if (CS%id_vhGM_Rlay>0) call safe_alloc_ptr(CS%vhGM_Rlay,isd,ied,JsdB,JedB,nz)
  !endif


  ! terms in the kinetic energy budget
  CS%id_KE = register_diag_field('ocean_model', 'KE', diag%axesTL, Time, &
      'Layer kinetic energy per unit mass', 'meter2 second-2')
  if (CS%id_KE>0) call safe_alloc_ptr(CS%KE,isd,ied,jsd,jed,nz)

  CS%id_dKEdt = register_diag_field('ocean_model', 'dKE_dt', diag%axesTL, Time, &
      'Kinetic Energy Tendency of Layer', 'meter3 second-3')
  if (CS%id_dKEdt>0) call safe_alloc_ptr(CS%dKE_dt,isd,ied,jsd,jed,nz)

  CS%id_PE_to_KE = register_diag_field('ocean_model', 'PE_to_KE', diag%axesTL, Time, &
      'Potential to Kinetic Energy Conversion of Layer', 'meter3 second-3')
  if (CS%id_PE_to_KE>0) call safe_alloc_ptr(CS%PE_to_KE,isd,ied,jsd,jed,nz)

  CS%id_KE_Coradv = register_diag_field('ocean_model', 'KE_Coradv', diag%axesTL, Time, &
      'Kinetic Energy Source from Coriolis and Advection', 'meter3 second-3')
  if (CS%id_KE_Coradv>0) call safe_alloc_ptr(CS%KE_Coradv,isd,ied,jsd,jed,nz)

  CS%id_KE_adv = register_diag_field('ocean_model', 'KE_adv', diag%axesTL, Time, &
      'Kinetic Energy Source from Advection', 'meter3 second-3')
  if (CS%id_KE_adv>0) call safe_alloc_ptr(CS%KE_adv,isd,ied,jsd,jed,nz)

  CS%id_KE_visc = register_diag_field('ocean_model', 'KE_visc', diag%axesTL, Time, &
      'Kinetic Energy Source from Vertical Viscosity and Stresses', 'meter3 second-3')
  if (CS%id_KE_visc>0) call safe_alloc_ptr(CS%KE_visc,isd,ied,jsd,jed,nz)

  CS%id_KE_horvisc = register_diag_field('ocean_model', 'KE_horvisc', diag%axesTL, Time, &
      'Kinetic Energy Source from Horizontal Viscosity', 'meter3 second-3')
  if (CS%id_KE_horvisc>0) call safe_alloc_ptr(CS%KE_horvisc,isd,ied,jsd,jed,nz)

  CS%id_KE_dia = register_diag_field('ocean_model', 'KE_dia', diag%axesTL, Time, &
      'Kinetic Energy Source from Diapycnal Diffusion', 'meter3 second-3')
  if (CS%id_KE_dia>0) call safe_alloc_ptr(CS%KE_dia,isd,ied,jsd,jed,nz)


  ! gravity wave CFLs
  CS%id_cg1 = register_diag_field('ocean_model', 'cg1', diag%axesT1, Time, &
      'First baroclinic gravity wave speed', 'meter second-1')
  CS%id_Rd1 = register_diag_field('ocean_model', 'Rd1', diag%axesT1, Time, &
      'First baroclinic deformation radius', 'meter')
  CS%id_cfl_cg1 = register_diag_field('ocean_model', 'CFL_cg1', diag%axesT1, Time, &
      'CFL of first baroclinic gravity wave = dt*cg1*(1/dx+1/dy)', 'nondim')
  CS%id_cfl_cg1_x = register_diag_field('ocean_model', 'CFL_cg1_x', diag%axesT1, Time, &
      'i-component of CFL of first baroclinic gravity wave = dt*cg1*/dx', 'nondim')
  CS%id_cfl_cg1_y = register_diag_field('ocean_model', 'CFL_cg1_y', diag%axesT1, Time, &
      'j-component of CFL of first baroclinic gravity wave = dt*cg1*/dy', 'nondim')
  CS%id_cg_ebt = register_diag_field('ocean_model', 'cg_ebt', diag%axesT1, Time, &
      'Equivalent barotropic gravity wave speed', 'meter second-1')
  CS%id_Rd_ebt = register_diag_field('ocean_model', 'Rd_ebt', diag%axesT1, Time, &
      'Equivalent barotropic deformation radius', 'meter')
  CS%id_p_ebt = register_diag_field('ocean_model', 'p_ebt', diag%axesTL, Time, &
      'Equivalent barotropic modal strcuture', 'nondim')

  if ((CS%id_cg1>0) .or. (CS%id_Rd1>0) .or. (CS%id_cfl_cg1>0) .or. &
      (CS%id_cfl_cg1_x>0) .or. (CS%id_cfl_cg1_y>0) .or. &
      (CS%id_cg_ebt>0) .or. (CS%id_Rd_ebt>0) .or. (CS%id_p_ebt>0)) then
    call wave_speed_init(CS%wave_speed_CSp)
    call safe_alloc_ptr(CS%cg1,isd,ied,jsd,jed)
    if (CS%id_Rd1>0)       call safe_alloc_ptr(CS%Rd1,isd,ied,jsd,jed)
    if (CS%id_Rd_ebt>0)    call safe_alloc_ptr(CS%Rd1,isd,ied,jsd,jed)
    if (CS%id_cfl_cg1>0)   call safe_alloc_ptr(CS%cfl_cg1,isd,ied,jsd,jed)
    if (CS%id_cfl_cg1_x>0) call safe_alloc_ptr(CS%cfl_cg1_x,isd,ied,jsd,jed)
    if (CS%id_cfl_cg1_y>0) call safe_alloc_ptr(CS%cfl_cg1_y,isd,ied,jsd,jed)
    if (CS%id_p_ebt>0) call safe_alloc_ptr(CS%p_ebt,isd,ied,jsd,jed,nz)
  endif

  CS%id_mass_wt = register_diag_field('ocean_model', 'mass_wt', diag%axesT1, Time,  &
      'The column mass for calculating mass-weighted average properties', 'kg m-2')

  CS%id_temp_int = register_diag_field('ocean_model', 'temp_int', diag%axesT1, Time,                &
      'Density weighted column integrated potential temperature', 'degC kg m-2',                    &
      cmor_field_name='opottempmint',                                                               &
      cmor_long_name='integral_wrt_depth_of_product_of_sea_water_density_and_potential_temperature',&
      cmor_units='degC kg m-2', cmor_standard_name='Depth integrated density times potential temperature')

  CS%id_salt_int = register_diag_field('ocean_model', 'salt_int', diag%axesT1, Time,   &
      'Density weighted column integrated salinity', 'ppt kg m-2',                     &
      cmor_field_name='somint',                                                        &
      cmor_long_name='integral_wrt_depth_of_product_of_sea_water_density_and_salinity',&
      cmor_units='ppt kg m-2', cmor_standard_name='Depth integrated density times salinity')

  CS%id_col_mass = register_diag_field('ocean_model', 'col_mass', diag%axesT1, Time, &
      'The column integrated in situ density', 'kg m-2')

  CS%id_col_ht = register_diag_field('ocean_model', 'col_height', diag%axesT1, Time, &
      'The height of the water column', 'm')
  CS%id_pbo = register_diag_field('ocean_model', 'pbo', diag%axesT1, Time, &
      long_name='Sea Water Pressure at Sea Floor', standard_name='sea_water_pressure_at_sea_floor', &
      units='Pa')

  call set_dependent_diagnostics(MIS, ADp, CDp, G, CS)

end subroutine MOM_diagnostics_init


subroutine set_dependent_diagnostics(MIS, ADp, CDp, G, CS)
  type(ocean_internal_state), intent(in)    :: MIS
  type(accel_diag_ptrs),      intent(inout) :: ADp
  type(cont_diag_ptrs),       intent(inout) :: CDp
  type(ocean_grid_type),      intent(in)    :: G
  type(diagnostics_CS),       pointer       :: CS

! This subroutine sets up diagnostics upon which other diagnostics depend.
! Arguments:
!  (in)      MIS - For "MOM Internal State" a set of pointers to the fields and
!                  accelerations making up ocean internal physical state.
!  (inout)   ADp - structure pointing to accelerations in momentum equation
!  (inout)   CDp - structure pointing to terms in continuity equation
!  (in)      G   - ocean grid structure
!  (in)      CS -  pointer to the control structure for this module

  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, nz
  isd  = G%isd  ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed ; nz = G%ke
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (ASSOCIATED(CS%dKE_dt) .or. ASSOCIATED(CS%PE_to_KE) .or. &
      ASSOCIATED(CS%KE_CorAdv) .or. ASSOCIATED(CS%KE_adv) .or. &
      ASSOCIATED(CS%KE_visc) .or. ASSOCIATED(CS%KE_horvisc) .or. &
      ASSOCIATED(CS%KE_dia)) &
    call safe_alloc_ptr(CS%KE,isd,ied,jsd,jed,nz)

  if (ASSOCIATED(CS%dKE_dt)) then
    if (.not.ASSOCIATED(CS%du_dt)) then
      call safe_alloc_ptr(CS%du_dt,IsdB,IedB,jsd,jed,nz)
      call register_time_deriv(MIS%u, CS%du_dt, CS)
    endif
    if (.not.ASSOCIATED(CS%dv_dt)) then
      call safe_alloc_ptr(CS%dv_dt,isd,ied,JsdB,JedB,nz)
      call register_time_deriv(MIS%v, CS%dv_dt, CS)
    endif
    if (.not.ASSOCIATED(CS%dh_dt)) then
      call safe_alloc_ptr(CS%dh_dt,isd,ied,jsd,jed,nz)
      call register_time_deriv(MIS%h, CS%dh_dt, CS)
    endif
  endif

  if (ASSOCIATED(CS%KE_adv)) then
    call safe_alloc_ptr(ADp%gradKEu,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%gradKEv,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%KE_visc)) then
    call safe_alloc_ptr(ADp%du_dt_visc,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_visc,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%KE_dia)) then
    call safe_alloc_ptr(ADp%du_dt_dia,IsdB,IedB,jsd,jed,nz)
    call safe_alloc_ptr(ADp%dv_dt_dia,isd,ied,JsdB,JedB,nz)
  endif

  if (ASSOCIATED(CS%uhGM_Rlay)) call safe_alloc_ptr(CDp%uhGM,IsdB,IedB,jsd,jed,nz)
  if (ASSOCIATED(CS%vhGM_Rlay)) call safe_alloc_ptr(CDp%vhGM,isd,ied,JsdB,JedB,nz)

end subroutine set_dependent_diagnostics


subroutine MOM_diagnostics_end(CS, ADp)
  type(diagnostics_CS),   pointer       :: CS
  type(accel_diag_ptrs),  intent(inout) :: ADp
  integer :: m

  if (ASSOCIATED(CS%e))          deallocate(CS%e)
  if (ASSOCIATED(CS%e_D))        deallocate(CS%e_D)
  if (ASSOCIATED(CS%KE))         deallocate(CS%KE)
  if (ASSOCIATED(CS%dKE_dt))     deallocate(CS%dKE_dt)
  if (ASSOCIATED(CS%PE_to_KE))   deallocate(CS%PE_to_KE)
  if (ASSOCIATED(CS%KE_Coradv))  deallocate(CS%KE_Coradv)
  if (ASSOCIATED(CS%KE_adv))     deallocate(CS%KE_adv)
  if (ASSOCIATED(CS%KE_visc))    deallocate(CS%KE_visc)
  if (ASSOCIATED(CS%KE_horvisc)) deallocate(CS%KE_horvisc)
  if (ASSOCIATED(CS%KE_dia))     deallocate(CS%KE_dia)
  if (ASSOCIATED(CS%dv_dt))      deallocate(CS%dv_dt)
  if (ASSOCIATED(CS%dh_dt))      deallocate(CS%dh_dt)
  if (ASSOCIATED(CS%du_dt))      deallocate(CS%du_dt)
  if (ASSOCIATED(CS%h_Rlay))     deallocate(CS%h_Rlay)
  if (ASSOCIATED(CS%uh_Rlay))    deallocate(CS%uh_Rlay)
  if (ASSOCIATED(CS%vh_Rlay))    deallocate(CS%vh_Rlay)
  if (ASSOCIATED(CS%uhGM_Rlay))  deallocate(CS%uhGM_Rlay)
  if (ASSOCIATED(CS%vhGM_Rlay))  deallocate(CS%vhGM_Rlay)
  if (ASSOCIATED(CS%diag_tmp3d)) deallocate(CS%diag_tmp3d)

  if (ASSOCIATED(ADp%gradKEu))    deallocate(ADp%gradKEu)
  if (ASSOCIATED(ADp%gradKEu))    deallocate(ADp%gradKEu)
  if (ASSOCIATED(ADp%du_dt_visc)) deallocate(ADp%du_dt_visc)
  if (ASSOCIATED(ADp%dv_dt_visc)) deallocate(ADp%dv_dt_visc)
  if (ASSOCIATED(ADp%du_dt_dia))  deallocate(ADp%du_dt_dia)
  if (ASSOCIATED(ADp%dv_dt_dia))  deallocate(ADp%dv_dt_dia)
  if (ASSOCIATED(ADp%du_other))   deallocate(ADp%du_other)
  if (ASSOCIATED(ADp%dv_other))   deallocate(ADp%dv_other)

  do m=1,CS%num_time_deriv ; deallocate(CS%prev_val(m)%p) ; enddo

  deallocate(CS)

end subroutine MOM_diagnostics_end

end module MOM_diagnostics
