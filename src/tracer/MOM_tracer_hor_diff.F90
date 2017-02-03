!> Main routine for lateral (along surface or neutral) diffusion of tracers
module MOM_tracer_hor_diff

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,             only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,         only : post_data, diag_ctrl
use MOM_diag_mediator,         only : register_diag_field, safe_alloc_ptr, time_type
use MOM_domains,               only : sum_across_PEs, max_across_PEs
use MOM_domains,               only : create_group_pass, do_group_pass, group_pass_type
use MOM_domains,               only : pass_vector
use MOM_debugging,             only : hchksum, uchksum, vchksum
use MOM_EOS,                   only : calculate_density
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_error_handler,         only : MOM_set_verbosity, callTree_showQuery
use MOM_error_handler,         only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_grid,                  only : ocean_grid_type
use MOM_lateral_mixing_coeffs, only : VarMix_CS
use MOM_MEKE_types,            only : MEKE_type
use MOM_neutral_diffusion,     only : neutral_diffusion_init, neutral_diffusion_end
use MOM_neutral_diffusion,     only : neutral_diffusion_CS
use MOM_neutral_diffusion,     only : neutral_diffusion_calc_coeffs, neutral_diffusion
use MOM_tracer_registry,       only : tracer_registry_type, tracer_type, MOM_tracer_chksum
use MOM_variables,             only : thermo_var_ptrs
use MOM_verticalGrid,          only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public tracer_hordiff, tracer_hor_diff_init, tracer_hor_diff_end

type, public :: tracer_hor_diff_CS ; private
  real    :: dt             ! The baroclinic dynamics time step, in s.
  real    :: KhTr           ! The along-isopycnal tracer diffusivity in m2/s.
  real    :: KhTr_Slope_Cff ! The non-dimensional coefficient in KhTr formula
  real    :: KhTr_min       ! Minimum along-isopycnal tracer diffusivity in m2/s.
  real    :: KhTr_max       ! Maximum along-isopycnal tracer diffusivity in m2/s.
  real    :: KhTr_passivity_coeff ! Passivity coefficient that scales Rd/dx (default = 0)
                                  ! where passivity is the ratio between along-isopycnal
                                  ! tracer mixing and thickness mixing
  real    :: KhTr_passivity_min   ! Passivity minimum (default = 1/2)
  real    :: ML_KhTR_scale  ! With Diffuse_ML_interior, the ratio of the truly
                            ! horizontal diffusivity in the mixed layer to the
                            ! epipycnal diffusivity.  Nondim.
  logical :: Diffuse_ML_interior  ! If true, diffuse along isopycnals between
                                  ! the mixed layer and the interior.
  logical :: check_diffusive_CFL  ! If true, automatically iterate the diffusion
                                  ! to ensure that the diffusive equivalent of the CFL
                                  ! limit is not violated.
  logical :: use_neutral_diffusion ! If true, use the neutral_diffusion module from within
                                   ! tracer_hor_diff.

  type(neutral_diffusion_CS), pointer :: neutral_diffusion_CSp => NULL() ! Control structure for neutral diffusion.
  type(diag_ctrl), pointer :: diag ! structure to regulate timing of diagnostic output.
  logical :: debug                 ! If true, write verbose checksums for debugging purposes.
  logical :: show_call_tree        ! Display the call tree while running. Set by VERBOSITY level.
  logical :: first_call = .true.
  integer :: id_KhTr_u  = -1
  integer :: id_KhTr_v  = -1
  integer :: id_KhTr_h  = -1
  integer :: id_CFL     = -1
  integer :: id_khdt_x  = -1
  integer :: id_khdt_y  = -1

  type(group_pass_type) :: pass_t !For group halo pass, used in both
                                  !tracer_hordiff and tracer_epipycnal_ML_diff
end type tracer_hor_diff_CS

type p2d
  real, dimension(:,:), pointer :: p => NULL()
end type p2d
type p2di
  integer, dimension(:,:), pointer :: p => NULL()
end type p2di

integer :: id_clock_diffuse, id_clock_epimix, id_clock_pass, id_clock_sync

contains

!> Compute along-coordinate diffusion of all tracers
!! using the diffusivity in CS%KhTr, or using space-dependent diffusivity.
!! Multiple iterations are used (if necessary) so that there is no limit
!! on the acceptable time increment.
subroutine tracer_hordiff(h, dt, MEKE, VarMix, G, GV, CS, Reg, tv, do_online_flag, read_khdt_x, read_khdt_y)
  type(ocean_grid_type),                 intent(inout) :: G       !< Grid type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h       !< Layer thickness (m or kg m-2)
  real,                                  intent(in)    :: dt      !< time step (seconds)
  type(MEKE_type),                       pointer       :: MEKE    !< MEKE type
  type(VarMix_CS),                       pointer       :: VarMix  !< Variable mixing type
  type(verticalGrid_type),               intent(in)    :: GV      !< ocean vertical grid structure
  type(tracer_hor_diff_CS),              pointer       :: CS      !< module control structure
  type(tracer_registry_type),            intent(inout) :: Reg     !< registered tracers
  type(thermo_var_ptrs),                 intent(in)    :: tv      !< A structure containing pointers to any available
                                                                  !! thermodynamic fields, including potential temp and
                                                                  !! salinity or mixed layer density. Absent fields have
                                                                  !! NULL ptrs, and these may (probably will) point to
                                                                  !! some of the same arrays as Tr does.  tv is required
                                                                  !! for epipycnal mixing between mixed layer and the interior.
  ! Optional inputs for offline tracer transport
  logical,                           optional             :: do_online_flag
  real, dimension(SZIB_(G),SZJ_(G)), optional, intent(in) :: read_khdt_x
  real, dimension(SZI_(G),SZJB_(G)), optional, intent(in) :: read_khdt_y


  real, dimension(SZI_(G),SZJ_(G)) :: &
    Ihdxdy, &     ! The inverse of the volume or mass of fluid in a layer in a
                  ! grid cell, in m-3 or kg-1.
    Kh_h, &       ! The tracer diffusivity averaged to tracer points, in m2 s-1.
    CFL, &        ! A diffusive CFL number for each cell, nondim.
    dTr           ! The change in a tracer's concentration, in units of
                  ! concentration.

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    khdt_x, &     ! The value of Khtr*dt times the open face width divided by
                  ! the distance between adjacent tracer points, in m2.
    Coef_x, &     ! The coefficients relating zonal tracer differences
                  ! to time-integrated fluxes, in m3 or kg.
    Kh_u          ! Tracer mixing coefficient at u-points, in m2 s-1.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    khdt_y, &     ! The value of Khtr*dt times the open face width divided by
                  ! the distance between adjacent tracer points, in m2.
    Coef_y, &     ! The coefficients relating meridional tracer differences
                  ! to time-integrated fluxes, in m3 or kg.
    Kh_v          ! Tracer mixing coefficient at u-points, in m2 s-1.

  real :: max_CFL ! The global maximum of the diffusive CFL number.
  logical :: use_VarMix, Resoln_scaled, do_online
  integer :: i, j, k, m, is, ie, js, je, nz, ntr, itt, num_itts
  real :: I_numitts  ! The inverse of the number of iterations, num_itts.
  real :: scale      ! The fraction of khdt_x or khdt_y that is applied in this
                     ! layer for this iteration, nondim.
  real :: Idt        ! The inverse of the time step, in s-1.
  real :: h_neglect  ! A thickness that is so small it is usually lost
                     ! in roundoff and can be neglected, in m.
  real :: Kh_loc     ! The local value of Kh, in m2 s-1.
  real :: Res_Fn     ! The local value of the resolution function, nondim.
  real :: Rd_dx      ! The local value of deformation radius over grid-spacing, nondim.
  real :: normalize  ! normalization used for diagnostic Kh_h; diffusivity averaged to h-points.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do_online = .true.
  if (present(do_online_flag)) do_online = do_online_flag

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_tracer_hor_diff: "// &
       "register_tracer must be called before tracer_hordiff.")
  if (LOC(Reg)==0) call MOM_error(FATAL, "MOM_tracer_hor_diff: "// &
       "register_tracer must be called before tracer_hordiff.")
  if ((Reg%ntr==0) .or. ((CS%KhTr <= 0.0) .and. .not.associated(VarMix)) ) return

  if (CS%show_call_tree) call callTree_enter("tracer_hordiff(), MOM_tracer_hor_diff.F90")

  call cpu_clock_begin(id_clock_diffuse)

  ntr = Reg%ntr
  Idt = 1.0/dt
  h_neglect = GV%H_subroundoff

  if (CS%Diffuse_ML_interior .and. CS%first_call) then ; if (is_root_pe()) then
    do m=1,ntr ; if (associated(Reg%Tr(m)%df_x) .or. associated(Reg%Tr(m)%df_y)) &
      call MOM_error(WARNING, "tracer_hordiff: Tracer "//trim(Reg%Tr(m)%name)// &
          " has associated 3-d diffusive flux diagnostics.  These are not "//&
          "valid when DIFFUSE_ML_TO_INTERIOR is defined. Use 2-d tracer "//&
          "diffusion diagnostics instead to get accurate total fluxes.")
    enddo
  endif ; endif
  CS%first_call = .false.

  if (CS%debug) call MOM_tracer_chksum("Before tracer diffusion ", Reg%Tr, ntr, G)

  use_VarMix = .false. ; Resoln_scaled = .false.
  if (Associated(VarMix)) then
    use_VarMix = VarMix%use_variable_mixing
    Resoln_scaled = VarMix%Resoln_scaled_KhTr
  endif

  call cpu_clock_begin(id_clock_pass)
  do m=1,ntr
    call create_group_pass(CS%pass_t, Reg%Tr(m)%t(:,:,:), G%Domain)
  enddo
  call cpu_clock_end(id_clock_pass)

  if (CS%show_call_tree) call callTree_waypoint("Calculating diffusivity (tracer_hordiff)")

  if (do_online) then
      if (use_VarMix) then
    !$OMP parallel default(none) shared(is,ie,js,je,CS,VarMix,MEKE,Resoln_scaled, &
    !$OMP                               Kh_u,Kh_v,khdt_x,dt,G,khdt_y)                        &
    !$OMP                       private(Kh_loc,Rd_dx)
    !$OMP do
        do j=js,je ; do I=is-1,ie
          Kh_loc = CS%KhTr + CS%KhTr_Slope_Cff*VarMix%L2u(I,j)*VarMix%SN_u(I,j)
          if (associated(MEKE%Kh)) &
            Kh_Loc = Kh_Loc + MEKE%KhTr_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i+1,j))
          if (CS%KhTr_max > 0.) Kh_loc = min(Kh_loc, CS%KhTr_max)
          if (Resoln_scaled) &
            Kh_Loc = Kh_Loc * 0.5*(VarMix%Res_fn_h(i,j) + VarMix%Res_fn_h(i+1,j))
          Kh_u(I,j) = max(Kh_loc, CS%KhTr_min)
          if (CS%KhTr_passivity_coeff>0.) then ! Apply passivity
            Rd_dx=0.5*( VarMix%Rd_dx_h(i,j)+VarMix%Rd_dx_h(i+1,j) ) ! Rd/dx at u-points
            Kh_loc=Kh_u(I,j)*max( CS%KhTr_passivity_min, CS%KhTr_passivity_coeff*Rd_dx )
            if (CS%KhTr_max > 0.) Kh_loc = min(Kh_loc, CS%KhTr_max) ! Re-apply max
            Kh_u(I,j) = max(Kh_loc, CS%KhTr_min) ! Re-apply min
          endif
        enddo ; enddo
    !$OMP do
        do J=js-1,je ;  do i=is,ie
          Kh_loc = CS%KhTr + CS%KhTr_Slope_Cff*VarMix%L2v(i,J)*VarMix%SN_v(i,J)
          if (associated(MEKE%Kh)) &
            Kh_Loc = Kh_Loc + MEKE%KhTr_fac*sqrt(MEKE%Kh(i,j)*MEKE%Kh(i,j+1))
          if (CS%KhTr_max > 0.) Kh_loc = min(Kh_loc, CS%KhTr_max)
          if (Resoln_scaled) &
            Kh_Loc = Kh_Loc * 0.5*(VarMix%Res_fn_h(i,j) + VarMix%Res_fn_h(i,j+1))
          Kh_v(i,J) = max(Kh_loc, CS%KhTr_min)
          if (CS%KhTr_passivity_coeff>0.) then ! Apply passivity
            Rd_dx=0.5*( VarMix%Rd_dx_h(i,j)+VarMix%Rd_dx_h(i,j+1) ) ! Rd/dx at v-points
            Kh_loc=Kh_v(I,j)*max( CS%KhTr_passivity_min, CS%KhTr_passivity_coeff*Rd_dx )
            if (CS%KhTr_max > 0.) Kh_loc = min(Kh_loc, CS%KhTr_max) ! Re-apply max
            Kh_v(i,J) = max(Kh_loc, CS%KhTr_min) ! Re-apply min
          endif
        enddo ; enddo

    !$OMP do
        do j=js,je ; do I=is-1,ie
          khdt_x(I,j) = dt*(Kh_u(I,j)*(G%dy_Cu(I,j)*G%IdxCu(I,j)))
        enddo ; enddo
    !$OMP do
        do J=js-1,je ; do i=is,ie
          khdt_y(i,J) = dt*(Kh_v(i,J)*(G%dx_Cv(i,J)*G%IdyCv(i,J)))
        enddo ; enddo
    !$OMP end parallel
      elseif (Resoln_scaled) then
    !$OMP parallel default(none) shared(is,ie,js,je,VarMix,Kh_u,Kh_v,khdt_x,khdt_y,CS,dt,G) &
    !$OMP                       private(Res_fn)
    !$OMP do
        do j=js,je ; do I=is-1,ie
          Res_fn = 0.5 * (VarMix%Res_fn_h(i,j) + VarMix%Res_fn_h(i+1,j))
          Kh_u(I,j) = max(CS%KhTr * Res_fn, CS%KhTr_min)
          khdt_x(I,j) = dt*(CS%KhTr*(G%dy_Cu(I,j)*G%IdxCu(I,j))) * Res_fn
        enddo ; enddo
    !$OMP do
        do J=js-1,je ;  do i=is,ie
          Res_fn = 0.5*(VarMix%Res_fn_h(i,j) + VarMix%Res_fn_h(i,j+1))
          Kh_v(i,J) = max(CS%KhTr * Res_fn, CS%KhTr_min)
          khdt_y(i,J) = dt*(CS%KhTr*(G%dx_Cv(i,J)*G%IdyCv(i,J))) * Res_fn
        enddo ; enddo
    !$OMP end parallel
      else
    !$OMP parallel default(none) shared(is,ie,js,je,Kh_u,Kh_v,khdt_x,khdt_y,CS,G,dt)
        if (CS%id_KhTr_u > 0) then
    !$OMP do
          do j=js,je ; do I=is-1,ie
            Kh_u(I,j) = CS%KhTr
            khdt_x(I,j) = dt*(CS%KhTr*(G%dy_Cu(I,j)*G%IdxCu(I,j)))
          enddo ; enddo
        else
    !$OMP do
          do j=js,je ; do I=is-1,ie
            khdt_x(I,j) = dt*(CS%KhTr*(G%dy_Cu(I,j)*G%IdxCu(I,j)))
          enddo ; enddo
        endif
        if (CS%id_KhTr_v > 0) then
    !$OMP do
          do J=js-1,je ;  do i=is,ie
            Kh_v(i,J) = CS%KhTr
            khdt_y(i,J) = dt*(CS%KhTr*(G%dx_Cv(i,J)*G%IdyCv(i,J)))
          enddo ; enddo
        else
    !$OMP do
          do J=js-1,je ;  do i=is,ie
            khdt_y(i,J) = dt*(CS%KhTr*(G%dx_Cv(i,J)*G%IdyCv(i,J)))
          enddo ; enddo
        endif
    !$OMP end parallel
      endif ! VarMix
  else ! .not. do_online
    khdt_x = read_khdt_x
    khdt_y = read_khdt_y
    call pass_vector(khdt_x,khdt_y,G%Domain)
  endif ! do_online



  if (CS%check_diffusive_CFL) then
    if (CS%show_call_tree) call callTree_waypoint("Checking diffusive CFL (tracer_hordiff)")
    max_CFL = 0.0
    do j=js,je ; do i=is,ie
      CFL(i,j) = 2.0*((khdt_x(I-1,j) + khdt_x(I,j)) + &
                      (khdt_y(i,J-1) + khdt_y(i,J))) * G%IareaT(i,j)
      if (max_CFL < CFL(i,j)) max_CFL = CFL(i,j)
    enddo ; enddo
    call cpu_clock_begin(id_clock_sync)
    call max_across_PEs(max_CFL)
    call cpu_clock_end(id_clock_sync)
    num_itts = max(1,ceiling(max_CFL))
    I_numitts = 1.0 ; if (num_itts > 1) I_numitts = 1.0 / (real(num_itts))
    if(CS%id_CFL > 0) call post_data(CS%id_CFL, CFL, CS%diag, mask=G%mask2dT)
  else
    num_itts = 1 ; I_numitts = 1.0
  endif

  do m=1,ntr
    if (associated(Reg%Tr(m)%df_x)) then
      do k=1,nz ; do j=js,je ; do I=is-1,ie
        Reg%Tr(m)%df_x(I,j,k) = 0.0
      enddo ; enddo ; enddo
    endif
    if (associated(Reg%Tr(m)%df_y)) then
      do k=1,nz ; do J=js-1,je ; do i=is,ie
        Reg%Tr(m)%df_y(i,J,k) = 0.0
      enddo ; enddo ; enddo
    endif
    if (associated(Reg%Tr(m)%df2d_x)) then
      do j=js,je ; do I=is-1,ie ; Reg%Tr(m)%df2d_x(I,j) = 0.0 ; enddo ; enddo
    endif
    if (associated(Reg%Tr(m)%df2d_y)) then
      do J=js-1,je ; do i=is,ie ; Reg%Tr(m)%df2d_y(i,J) = 0.0 ; enddo ; enddo
    endif
  enddo

  if (CS%use_neutral_diffusion) then

    if (CS%show_call_tree) call callTree_waypoint("Calling neutral diffusion coeffs (tracer_hordiff)")
    call cpu_clock_begin(id_clock_pass)
    call do_group_pass(CS%pass_t, G%Domain)
    call cpu_clock_end(id_clock_pass)
    ! We are assuming that neutral surfaces do not evolve (much) as a result of multiple
    ! lateral diffusion iterations. Otherwise the call to neutral_diffusion_calc_coeffs()
    ! would be inside the itt-loop. -AJA

    call neutral_diffusion_calc_coeffs(G, GV, h, tv%T, tv%S, tv%eqn_of_state, &
                                       CS%neutral_diffusion_CSp)
    do J=js-1,je ; do i=is,ie
      Coef_y(i,J) = I_numitts * khdt_y(i,J)
    enddo ; enddo
    do j=js,je
      do I=is-1,ie
        Coef_x(I,j) = I_numitts * khdt_x(I,j)
      enddo
    enddo

    do itt=1,num_itts
      if (CS%show_call_tree) call callTree_waypoint("Calling neutral diffusion (tracer_hordiff)",itt)
      if (itt>1) then ! Update halos for subsequent iterations
        call cpu_clock_begin(id_clock_pass)
        call do_group_pass(CS%pass_t, G%Domain)
        call cpu_clock_end(id_clock_pass)
      endif
      do m=1,ntr ! for each tracer
        call neutral_diffusion(G, GV,  h, Coef_x, Coef_y, Reg%Tr(m)%t, m, I_numitts*dt, &
                               Reg%Tr(m)%name, CS%neutral_diffusion_CSp)
      enddo ! m
    enddo ! itt

  else    ! following if not using neutral diffusion, but instead along-surface diffusion

    if (CS%show_call_tree) call callTree_waypoint("Calculating horizontal diffusion (tracer_hordiff)")
    do itt=1,num_itts
      call cpu_clock_begin(id_clock_pass)
      call do_group_pass(CS%pass_t, G%Domain)
      call cpu_clock_end(id_clock_pass)
!$OMP parallel do default(none) shared(is,ie,js,je,nz,I_numitts,CS,G,GV,khdt_y,h, &
!$OMP                                    h_neglect,khdt_x,ntr,Idt,Reg)           &
!$OMP                            private(scale,Coef_y,Coef_x,Ihdxdy,dTr)
      do k=1,nz
        scale = I_numitts
        if (CS%Diffuse_ML_interior) then
          if (k<=GV%nkml) then
            if (CS%ML_KhTr_scale <= 0.0) cycle
            scale = I_numitts * CS%ML_KhTr_scale
          endif
          if ((k>GV%nkml) .and. (k<=GV%nk_rho_varies)) cycle
        endif

        do J=js-1,je ; do i=is,ie
          Coef_y(i,J) = ((scale * khdt_y(i,J))*2.0*(h(i,j,k)*h(i,j+1,k))) / &
                                                   (h(i,j,k)+h(i,j+1,k)+h_neglect)
        enddo ; enddo

        do j=js,je
          do I=is-1,ie
            Coef_x(I,j) = ((scale * khdt_x(I,j))*2.0*(h(i,j,k)*h(i+1,j,k))) / &
                                                     (h(i,j,k)+h(i+1,j,k)+h_neglect)
          enddo

          do i=is,ie
            Ihdxdy(i,j) = G%IareaT(i,j) / (h(i,j,k)+h_neglect)
          enddo
        enddo

        do m=1,ntr
          do j=js,je ; do i=is,ie
            dTr(i,j) = Ihdxdy(i,j) * &
              ((Coef_x(I-1,j) * (Reg%Tr(m)%t(i-1,j,k) - Reg%Tr(m)%t(i,j,k)) - &
                Coef_x(I,j) * (Reg%Tr(m)%t(i,j,k) - Reg%Tr(m)%t(i+1,j,k))) + &
               (Coef_y(i,J-1) * (Reg%Tr(m)%t(i,j-1,k) - Reg%Tr(m)%t(i,j,k)) - &
                Coef_y(i,J) * (Reg%Tr(m)%t(i,j,k) - Reg%Tr(m)%t(i,j+1,k))))
          enddo ; enddo
          if (associated(Reg%Tr(m)%df_x)) then ; do j=js,je ; do I=G%IscB,G%IecB
            Reg%Tr(m)%df_x(I,j,k) = Reg%Tr(m)%df_x(I,j,k) + Coef_x(I,j) * &
                                (Reg%Tr(m)%t(i,j,k) - Reg%Tr(m)%t(i+1,j,k))*Idt
          enddo ; enddo ; endif
          if (associated(Reg%Tr(m)%df_y)) then ; do J=G%JscB,G%JecB ; do i=is,ie
            Reg%Tr(m)%df_y(i,J,k) = Reg%Tr(m)%df_y(i,J,k) + Coef_y(i,J) * &
                                (Reg%Tr(m)%t(i,j,k) - Reg%Tr(m)%t(i,j+1,k))*Idt
          enddo ; enddo ; endif
          if (associated(Reg%Tr(m)%df2d_x)) then ; do j=js,je ; do I=G%IscB,G%IecB
            Reg%Tr(m)%df2d_x(I,j) = Reg%Tr(m)%df2d_x(I,j) + Coef_x(I,j) * &
                                (Reg%Tr(m)%t(i,j,k) - Reg%Tr(m)%t(i+1,j,k))*Idt
          enddo ; enddo ; endif
          if (associated(Reg%Tr(m)%df2d_y)) then ; do J=G%JscB,G%JecB ; do i=is,ie
            Reg%Tr(m)%df2d_y(i,J) = Reg%Tr(m)%df2d_y(i,J) + Coef_y(i,J) * &
                                (Reg%Tr(m)%t(i,j,k) - Reg%Tr(m)%t(i,j+1,k))*Idt
          enddo ; enddo ; endif
          do j=js,je ; do i=is,ie
            Reg%Tr(m)%t(i,j,k) = Reg%Tr(m)%t(i,j,k) + dTr(i,j)
          enddo ; enddo
        enddo

      enddo ! End of k loop.

    enddo ! End of "while" loop.

  endif   ! endif for CS%use_neutral_diffusion
  call cpu_clock_end(id_clock_diffuse)


  if (CS%Diffuse_ML_interior) then
    if (CS%show_call_tree) call callTree_waypoint("Calling epipycnal_ML_diff (tracer_hordiff)")
    if (CS%debug) call MOM_tracer_chksum("Before epipycnal diff ", Reg%Tr, ntr, G)

    call cpu_clock_begin(id_clock_epimix)
    call tracer_epipycnal_ML_diff(h, dt, Reg%Tr, ntr, khdt_x, khdt_y, G, GV, &
                                  CS, tv, num_itts)
    call cpu_clock_end(id_clock_epimix)
  endif

  if (CS%debug) call MOM_tracer_chksum("After tracer diffusion ", Reg%Tr, ntr, G)

  ! post diagnostics for 2d tracer diffusivity
  if (CS%id_KhTr_u > 0) then
    do j=js,je ; do I=is-1,ie
      Kh_u(I,j) = G%mask2dCu(I,j)*Kh_u(I,j)
    enddo ; enddo
    call post_data(CS%id_KhTr_u, Kh_u, CS%diag, mask=G%mask2dCu)
  endif
  if (CS%id_KhTr_v > 0) then
    do J=js-1,je ; do i=is,ie
      Kh_v(i,J) = G%mask2dCv(i,J)*Kh_v(i,J)
    enddo ; enddo
    call post_data(CS%id_KhTr_v, Kh_v, CS%diag, mask=G%mask2dCv)
  endif
  if (CS%id_KhTr_h > 0) then
    Kh_h(:,:) = 0.0
    do j=js,je ; do I=is-1,ie
      Kh_u(I,j) = G%mask2dCu(I,j)*Kh_u(I,j)
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      Kh_v(i,J) = G%mask2dCv(i,J)*Kh_v(i,J)
    enddo ; enddo
    do j=js,je ; do i=is,ie
      normalize = 1.0 / ((G%mask2dCu(I-1,j)+G%mask2dCu(I,j)) + &
                  (G%mask2dCv(i,J-1)+G%mask2dCv(i,J)) + GV%H_subroundoff)
      Kh_h(i,j) = normalize*G%mask2dT(i,j)*((Kh_u(I-1,j)+Kh_u(I,j)) + &
                                            (Kh_v(i,J-1)+Kh_v(i,J)))
    enddo ; enddo
    call post_data(CS%id_KhTr_h, Kh_h, CS%diag, mask=G%mask2dT)
  endif


  if (CS%debug) then
    call uchksum(khdt_x,"After tracer diffusion khdt_x", G%HI, haloshift=2)
    call vchksum(khdt_y,"After tracer diffusion khdt_y", G%HI, haloshift=2)
    call uchksum(Coef_x,"After tracer diffusion Coef_x", G%HI, haloshift=2)
    call vchksum(Coef_y,"After tracer diffusion Coef_y", G%HI, haloshift=2)
  endif

  if (CS%id_khdt_x > 0) call post_data(CS%id_khdt_x, khdt_x, CS%diag)
  if (CS%id_khdt_y > 0) call post_data(CS%id_khdt_y, khdt_y, CS%diag)

  if (CS%show_call_tree) call callTree_leave("tracer_hordiff()")

end subroutine tracer_hordiff

!> This subroutine does epipycnal diffusion of all tracers between the mixed
!! and buffer layers and the interior, using the diffusivity in CS%KhTr.
!! Multiple iterations are used (if necessary) so that there is no limit on the
!! acceptable time increment.
subroutine tracer_epipycnal_ML_diff(h, dt, Tr, ntr, khdt_epi_x, khdt_epi_y, G, &
                                    GV, CS, tv, num_itts)
  type(ocean_grid_type),                    intent(inout) :: G          !< ocean grid structure
  type(verticalGrid_type),                  intent(in)    :: GV         !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h          !< layer thickness (m or kg m-2)
  real,                                     intent(in)    :: dt         !< time step
  type(tracer_type),                        intent(inout) :: Tr(:)      !< tracer array
  integer,                                  intent(in)    :: ntr        !< number of tracers
  real, dimension(SZIB_(G),SZJ_(G)),        intent(in)    :: khdt_epi_x !< needs a comment
  real, dimension(SZI_(G),SZJB_(G)),        intent(in)    :: khdt_epi_y !< needs a comment
  type(tracer_hor_diff_CS),                 intent(inout) :: CS         !< module control structure
  type(thermo_var_ptrs),                    intent(in)    :: tv         !< thermodynamic structure
  integer,                                  intent(in)    :: num_itts   !< number of iterations (usually=1)


  real, dimension(SZI_(G), SZJ_(G)) :: &
    Rml_max  ! The maximum coordinate density within the mixed layer, in kg m-3.
  real, dimension(SZI_(G), SZJ_(G), max(1,GV%nk_rho_varies)) :: &
    rho_coord ! The coordinate density that is used to mix along, in kg m-3.

  ! The naming mnemnonic is a=above,b=below,L=Left,R=Right,u=u-point,v=v-point.
  ! These are 1-D arrays of pointers to 2-d arrays to minimize memory usage.
  type(p2d), dimension(SZJ_(G)) :: &
    deep_wt_Lu, deep_wt_Ru, &  ! The relative weighting of the deeper of a pair, ND.
    hP_Lu, hP_Ru       ! The total thickness on each side for each pair, in m or kg m-2.

  type(p2d), dimension(SZJB_(G)) :: &
    deep_wt_Lv, deep_wt_Rv, & ! The relative weighting of the deeper of a pair, ND.
    hP_Lv, hP_Rv       ! The total thickness on each side for each pair, in m or kg m-2.

  type(p2di), dimension(SZJ_(G)) :: &
    k0b_Lu, k0a_Lu, &  ! The original k-indices of the layers that participate
    k0b_Ru, k0a_Ru     ! in each pair of mixing at u-faces.
  type(p2di), dimension(SZJB_(G)) :: &
    k0b_Lv, k0a_Lv, &  ! The original k-indices of the layers that participate
    k0b_Rv, k0a_Rv     ! in each pair of mixing at v-faces.

  real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: &
    tr_flux_conv  ! The flux convergence of tracers, in TR m3 or TR kg.
  real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: Tr_flux_3d, Tr_adj_vert_L, Tr_adj_vert_R

  real, dimension(SZI_(G), SZK_(G), SZJ_(G)) :: &
    rho_srt, & ! The density of each layer of the sorted columns, in kg m-3.
    h_srt      ! The thickness of each layer of the sorted columns, in m or kg m-2.
  integer, dimension(SZI_(G), SZK_(G), SZJ_(G)) :: &
    k0_srt     ! The original k-index that each layer of the sorted column
               ! corresponds to.

  real, dimension(SZK_(G)) :: &
    h_demand_L, & ! The thickness in the left (_L) or right (_R) column that
    h_demand_R, & ! is demanded to match the thickness in the counterpart, in H.
    h_used_L, &   ! The summed thickness from the left or right columns that
    h_used_R, &   ! have actually been used, in m or kg m-2 (H).
    h_supply_frac_L, &  ! The fraction of the demanded thickness that can
    h_supply_frac_R     ! actually be supplied from a layer.
  integer, dimension(SZK_(G)) :: &
    kbs_Lp, &   ! The sorted indicies of the Left and Right columns for
    kbs_Rp      ! each pairing.

  integer, dimension(SZI_(G), SZJ_(G))  :: &
    num_srt, &   ! The number of layers that are sorted in each column.
    k_end_srt, & ! The maximum index in each column that might need to be
                 ! sorted, based on neighboring values of max_kRho
    max_kRho     ! The index of the layer whose target density is just denser
                 ! than the densest part of the mixed layer.
  integer, dimension(SZJ_(G))           :: &
    max_srt      ! The maximum value of num_srt in a k-row.
  integer, dimension(SZIB_(G), SZJ_(G)) :: &
    nPu          ! The number of epipycnal pairings at each u-point.
  integer, dimension(SZI_(G), SZJB_(G)) :: &
    nPv          ! The number of epipycnal pairings at each v-point.
  real :: h_exclude    ! A thickness that layers must attain to be considered
                       ! for inclusion in mixing, in m.
  real :: Idt        ! The inverse of the time step, in s-1.
  real :: I_maxitt   ! The inverse of the maximum number of iterations.
  real :: rho_pair, rho_a, rho_b  ! Temporary densities, in kg m-3.
  real :: Tr_min_face  ! The minimum and maximum tracer concentrations
  real :: Tr_max_face  ! associated with a pairing, in conc.
  real :: Tr_La, Tr_Lb ! The 4 tracer concentrations that might be
  real :: Tr_Ra, Tr_Rb ! associated with a pairing, in conc.
  real :: Tr_av_L    ! The average tracer concentrations on the left and right
  real :: Tr_av_R    ! sides of a pairing, in conc.
  real :: Tr_flux    ! The tracer flux from left to right in a pair, in conc m3.
  real :: Tr_adj_vert  ! A downward vertical adjustment to Tr_flux between the
                     ! two cells that make up one side of the pairing, in conc m3.
  real :: h_L, h_R   ! Thicknesses to the left and right, in m or kg m-2 (H).
  real :: wt_a, wt_b ! Fractional weights of layers above and below, ND.
  real :: vol        ! A cell volume or mass, in m3 or kg (H m2).
  logical, dimension(SZK_(G)) :: &
    left_set, &  ! If true, the left or right point determines the density of
    right_set    ! of the trio.  If densities are exactly equal, both are true.
  real :: tmp
  real :: p_ref_cv(SZI_(G))

  integer :: k_max, k_min, k_test, itmp
  integer :: i, j, k, k2, m, is, ie, js, je, nz, nkmb
  integer :: isd, ied, jsd, jed, IsdB, IedB, k_size
  integer :: kL, kR, kLa, kLb, kRa, kRb, nP, itt, ns, max_itt
  integer :: PEmax_kRho
  integer :: isv, iev, jsv, jev ! The valid range of the indices.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB
  Idt = 1.0/dt
  nkmb = GV%nk_rho_varies

  if (num_itts <= 1) then
    max_itt = 1 ; I_maxitt = 1.0
  else
    max_itt = num_itts ; I_maxitt = 1.0 / (real(max_itt))
  endif

  do i=is-2,ie+2 ; p_ref_cv(i) = tv%P_Ref ; enddo

  call cpu_clock_begin(id_clock_pass)
  call do_group_pass(CS%pass_t, G%Domain)
  call cpu_clock_end(id_clock_pass)
  ! Determine which layers the mixed- and buffer-layers map into...
!$OMP parallel do default(none) shared(nkmb,is,ie,js,je,tv,p_ref_cv,rho_coord)
  do k=1,nkmb
    do j=js-2,je+2
       call calculate_density(tv%T(:,j,k),tv%S(:,j,k), p_ref_cv, &
                           rho_coord(:,j,k), is-2, ie-is+5, tv%eqn_of_state)
    enddo
 enddo

  do j=js-2,je+2 ; do i=is-2,ie+2
    Rml_max(i,j) = rho_coord(i,j,1)
    num_srt(i,j) = 0 ; max_kRho(i,j) = 0
  enddo ; enddo
  do k=2,nkmb ; do j=js-2,je+2 ; do i=is-2,ie+2
    if (Rml_max(i,j) < rho_coord(i,j,k)) Rml_max(i,j) = rho_coord(i,j,k)
  enddo ; enddo ; enddo
  !   Use bracketing and bisection to find the k-level that the densest of the
  ! mixed and buffer layer corresponds to, such that:
  !     GV%Rlay(max_kRho-1) < Rml_max <= GV%Rlay(max_kRho)
!$OMP parallel do default(none) shared(is,ie,js,je,nz,nkmb,G,GV,Rml_max,max_kRho) &
!$OMP                          private(k_min,k_max,k_test)
  do j=js-2,je+2 ; do i=is-2,ie+2 ; if (G%mask2dT(i,j) > 0.5) then
    if (Rml_max(i,j) > GV%Rlay(nz)) then ; max_kRho(i,j) = nz+1
    elseif (Rml_max(i,j) <= GV%Rlay(nkmb+1)) then ; max_kRho(i,j) = nkmb+1
    else
      k_min = nkmb+2 ; k_max = nz
      do
        k_test = (k_min + k_max) / 2
        if (Rml_max(i,j) <= GV%Rlay(k_test-1)) then ; k_max = k_test-1
        elseif (GV%Rlay(k_test) < Rml_max(i,j)) then ; k_min = k_test+1
        else ; max_kRho(i,j) = k_test ; exit ; endif

        if (k_min == k_max) then ; max_kRho(i,j) = k_max ; exit ; endif
      enddo
    endif
  endif ; enddo ; enddo

  PEmax_kRho = 0
  do j=js-1,je+1 ; do i=is-1,ie+1
    k_end_srt(i,j) = max(max_kRho(i,j), max_kRho(i-1,j), max_kRho(i+1,j), &
                         max_kRho(i,j-1), max_kRho(i,j+1))
    if (PEmax_kRho < k_end_srt(i,j)) PEmax_kRho = k_end_srt(i,j)
  enddo ; enddo
  if (PEmax_kRho > nz) PEmax_kRho = nz ! PEmax_kRho could have been nz+1.

  h_exclude = 10.0*(GV%Angstrom + GV%H_subroundoff)
!$OMP parallel default(none) shared(is,ie,js,je,nkmb,G,GV,h,h_exclude,num_srt,k0_srt, &
!$OMP                               rho_srt,h_srt,PEmax_kRho,k_end_srt,rho_coord,max_srt) &
!$OMP                       private(ns,tmp,itmp)
!$OMP do
  do j=js-1,je+1
    do k=1,nkmb ; do i=is-1,ie+1 ; if (G%mask2dT(i,j) > 0.5) then
      if (h(i,j,k) > h_exclude) then
        num_srt(i,j) = num_srt(i,j) + 1 ; ns = num_srt(i,j)
        k0_srt(i,ns,j) = k
        rho_srt(i,ns,j) = rho_coord(i,j,k)
        h_srt(i,ns,j) = h(i,j,k)
      endif
    endif ; enddo ; enddo
    do k=nkmb+1,PEmax_kRho ; do i=is-1,ie+1 ; if (G%mask2dT(i,j) > 0.5) then
      if ((k<=k_end_srt(i,j)) .and. (h(i,j,k) > h_exclude)) then
        num_srt(i,j) = num_srt(i,j) + 1 ; ns = num_srt(i,j)
        k0_srt(i,ns,j) = k
        rho_srt(i,ns,j) = GV%Rlay(k)
        h_srt(i,ns,j) = h(i,j,k)
      endif
    endif ; enddo ; enddo
  enddo
  ! Sort each column by increasing density.  This should already be close,
  ! and the size of the arrays are small, so straight insertion is used.
!$OMP do
   do j=js-1,je+1; do i=is-1,ie+1
    do k=2,num_srt(i,j) ; if (rho_srt(i,k,j) < rho_srt(i,k-1,j)) then
      ! The last segment needs to be shuffled earlier in the list.
      do k2 = k,2,-1 ; if (rho_srt(i,k2,j) >= rho_srt(i,k2-1,j)) exit
        itmp = k0_srt(i,k2-1,j) ; k0_srt(i,k2-1,j) = k0_srt(i,k2,j) ; k0_srt(i,k2,j) = itmp
        tmp = rho_srt(i,k2-1,j) ; rho_srt(i,k2-1,j) = rho_srt(i,k2,j) ; rho_srt(i,k2,j) = tmp
        tmp = h_srt(i,k2-1,j) ; h_srt(i,k2-1,j) = h_srt(i,k2,j) ; h_srt(i,k2,j) = tmp
      enddo
    endif ; enddo
  enddo; enddo
!$OMP do
  do j=js-1,je+1
    max_srt(j) = 0
    do i=is-1,ie+1 ; max_srt(j) = max(max_srt(j), num_srt(i,j)) ; enddo
  enddo
!$OMP end parallel

  do j=js,je
    k_size = max(2*max_srt(j),1)
    allocate(deep_wt_Lu(j)%p(IsdB:IedB,k_size))
    allocate(deep_wt_Ru(j)%p(IsdB:IedB,k_size))
    allocate(hP_Lu(j)%p(IsdB:IedB,k_size))
    allocate(hP_Ru(j)%p(IsdB:IedB,k_size))
    allocate(k0a_Lu(j)%p(IsdB:IedB,k_size))
    allocate(k0a_Ru(j)%p(IsdB:IedB,k_size))
    allocate(k0b_Lu(j)%p(IsdB:IedB,k_size))
    allocate(k0b_Ru(j)%p(IsdB:IedB,k_size))
  enddo

!$OMP parallel do default(none) shared(is,ie,js,je,G,num_srt,rho_srt,k0b_Lu,k0_srt, &
!$OMP                                  k0b_Ru,k0a_Lu,k0a_Ru,deep_wt_Lu,deep_wt_Ru,  &
!$OMP                                  h_srt,nkmb,nPu,hP_Lu,hP_Ru)                  &
!$OMP                          private(h_demand_L,h_used_L,h_demand_R,h_used_R,     &
!$OMP                                  kR,kL,nP,rho_pair,kbs_Lp,kbs_Rp,rho_a,rho_b, &
!$OMP                                  wt_b,left_set,right_set,h_supply_frac_R,     &
!$OMP                                  h_supply_frac_L)
  do j=js,je ; do I=is-1,ie ; if (G%mask2dCu(I,j) > 0.5) then
    ! Set up the pairings for fluxes through the zonal faces.

    do k=1,num_srt(i,j)   ; h_demand_L(k) = 0.0 ; h_used_L(k) = 0.0 ; enddo
    do k=1,num_srt(i+1,j) ; h_demand_R(k) = 0.0 ; h_used_R(k) = 0.0 ; enddo

    ! First merge the left and right lists into a single, sorted list.

    !   Discard any layers that are lighter than the lightest in the other
    ! column.  They can only participate in mixing as the lighter part of a
    ! pair of points.
    if (rho_srt(i,1,j) < rho_srt(i+1,1,j)) then
      kR = 1
      do kL=2,num_srt(i,j) ; if (rho_srt(i,kL,j) >= rho_srt(i+1,1,j)) exit ; enddo
    elseif (rho_srt(i+1,1,j) < rho_srt(i,1,j)) then
      kL = 1
      do kR=2,num_srt(i+1,j) ; if (rho_srt(i+1,kR,j) >= rho_srt(i,1,j)) exit ; enddo
    else
      kL = 1 ; kR = 1
    endif
    nP = 0
    do ! Loop to accumulate pairs of columns.
      if ((kL > num_srt(i,j)) .or. (kR > num_srt(i+1,j))) exit

      if (rho_srt(i,kL,j) > rho_srt(i+1,kR,j)) then
      ! The right point is lighter and defines the density for this trio.
        nP = nP+1 ; k = nP
        rho_pair = rho_srt(i+1,kR,j)

        k0b_Lu(j)%p(I,k) = k0_srt(i,kL,j) ; k0b_Ru(j)%p(I,k) = k0_srt(i+1,kR,j)
        k0a_Lu(j)%p(I,k) = k0_srt(i,kL-1,j) ; k0a_Ru(j)%p(I,k) = k0b_Ru(j)%p(I,k)
        kbs_Lp(k) = kL ; kbs_Rp(k) = kR

        rho_a = rho_srt(i,kL-1,j) ; rho_b = rho_srt(i,kL,j)
        wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
          wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
        deep_wt_Lu(j)%p(I,k) = wt_b ; deep_wt_Ru(j)%p(I,k) = 1.0

        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i+1,kR,j) * wt_b
        h_demand_L(kL-1) = h_demand_L(kL-1) + 0.5*h_srt(i+1,kR,j) * (1.0-wt_b)

        kR = kR+1 ; left_set(k) = .false. ; right_set(k) = .true.
      elseif (rho_srt(i,kL,j) < rho_srt(i+1,kR,j)) then
      ! The left point is lighter and defines the density for this trio.
        nP = nP+1 ; k = nP
        rho_pair = rho_srt(i,kL,j)
        k0b_Lu(j)%p(I,k) = k0_srt(i,kL,j) ; k0b_Ru(j)%p(I,k) = k0_srt(i+1,kR,j)
        k0a_Lu(j)%p(I,k) = k0b_Lu(j)%p(I,k) ; k0a_Ru(j)%p(I,k) = k0_srt(i+1,kR-1,j)

        kbs_Lp(k) = kL ; kbs_Rp(k) = kR

        rho_a = rho_srt(i+1,kR-1,j) ; rho_b = rho_srt(i+1,kR,j)
        wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
          wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
        deep_wt_Lu(j)%p(I,k) = 1.0 ; deep_wt_Ru(j)%p(I,k) = wt_b

        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j) * wt_b
        h_demand_R(kR-1) = h_demand_R(kR-1) + 0.5*h_srt(i,kL,j) * (1.0-wt_b)

        kL = kL+1 ; left_set(k) = .true. ; right_set(k) = .false.
      elseif ((k0_srt(i,kL,j) <= nkmb) .or. (k0_srt(i+1,kR,j) <= nkmb)) then
        ! The densities are exactly equal and one layer is above the interior.
        nP = nP+1 ; k = nP
        k0b_Lu(j)%p(I,k) = k0_srt(i,kL,j) ; k0b_Ru(j)%p(I,k) = k0_srt(i+1,kR,j)
        k0a_Lu(j)%p(I,k) = k0b_Lu(j)%p(I,k) ; k0a_Ru(j)%p(I,k) = k0b_Ru(j)%p(I,k)
        kbs_Lp(k) = kL ; kbs_Rp(k) = kR
        deep_wt_Lu(j)%p(I,k) = 1.0 ; deep_wt_Ru(j)%p(I,k) = 1.0

        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i+1,kR,j)
        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)

        kL = kL+1 ; kR = kR+1 ; left_set(k) = .true. ; right_set(k) = .true.
      else ! The densities are exactly equal and in the interior.
        ! Mixing in this case has already occurred, so accumulate the thickness
        ! demanded for that mixing and skip onward.
        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i+1,kR,j)
        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)

        kL = kL+1 ; kR = kR+1
      endif
    enddo ! Loop to accumulate pairs of columns.
    nPu(I,j) = nP ! This is the number of active pairings.

    ! Determine what fraction of the thickness "demand" can be supplied.
    do k=1,num_srt(i+1,j)
      h_supply_frac_R(k) = 1.0
      if (h_demand_R(k) > 0.5*h_srt(i+1,k,j)) &
        h_supply_frac_R(k) = 0.5*h_srt(i+1,k,j) / h_demand_R(k)
    enddo
    do k=1,num_srt(i,j)
      h_supply_frac_L(k) = 1.0
      if (h_demand_L(k) > 0.5*h_srt(i,k,j)) &
        h_supply_frac_L(k) = 0.5*h_srt(i,k,j) / h_demand_L(k)
    enddo

    !  Distribute the "exported" thicknesses proportionately.
    do k=1,nPu(I,j)
      kL = kbs_Lp(k) ; kR = kbs_Rp(k)
      hP_Lu(j)%p(I,k) = 0.0 ; hP_Ru(j)%p(I,k) = 0.0
      if (left_set(k)) then ! Add the contributing thicknesses on the right.
        if (deep_wt_Ru(j)%p(I,k) < 1.0) then
          hP_Ru(j)%p(I,k) = 0.5*h_srt(i,kL,j) * min(h_supply_frac_R(kR), h_supply_frac_R(kR-1))
          wt_b = deep_wt_Ru(j)%p(I,k)
          h_used_R(kR-1) = h_used_R(kR-1) + (1.0 - wt_b)*hP_Ru(j)%p(I,k)
          h_used_R(kR) = h_used_R(kR) + wt_b*hP_Ru(j)%p(I,k)
        else
          hP_Ru(j)%p(I,k) = 0.5*h_srt(i,kL,j) * h_supply_frac_R(kR)
          h_used_R(kR) = h_used_R(kR) + hP_Ru(j)%p(I,k)
        endif
      endif
      if (right_set(k)) then ! Add the contributing thicknesses on the left.
        if (deep_wt_Lu(j)%p(I,k) < 1.0) then
          hP_Lu(j)%p(I,k) = 0.5*h_srt(i+1,kR,j) * min(h_supply_frac_L(kL), h_supply_frac_L(kL-1))
          wt_b = deep_wt_Lu(j)%p(I,k)
          h_used_L(kL-1) = h_used_L(kL-1) + (1.0 - wt_b)*hP_Lu(j)%p(I,k)
          h_used_L(kL) = h_used_L(kL) + wt_b*hP_Lu(j)%p(I,k)
        else
          hP_Lu(j)%p(I,k) = 0.5*h_srt(i+1,kR,j) * h_supply_frac_L(kL)
          h_used_L(kL) = h_used_L(kL) + hP_Lu(j)%p(I,k)
        endif
      endif
    enddo

    !   The left-over thickness (at least half the layer thickness) is now
    ! added to the thicknesses of the importing columns.
    do k=1,nPu(I,j)
      if (left_set(k)) hP_Lu(j)%p(I,k) = hP_Lu(j)%p(I,k) + &
                           (h_srt(i,kbs_Lp(k),j) - h_used_L(kbs_Lp(k)))
      if (right_set(k)) hP_Ru(j)%p(I,k) = hP_Ru(j)%p(I,k) + &
                            (h_srt(i+1,kbs_Rp(k),j) - h_used_R(kbs_Rp(k)))
    enddo

  endif ; enddo ; enddo ! i- & j- loops over zonal faces.

  do J=js-1,je
    k_size = max(max_srt(j)+max_srt(j+1),1)
    allocate(deep_wt_Lv(J)%p(isd:ied,k_size))
    allocate(deep_wt_Rv(J)%p(isd:ied,k_size))
    allocate(hP_Lv(J)%p(isd:ied,k_size))
    allocate(hP_Rv(J)%p(isd:ied,k_size))
    allocate(k0a_Lv(J)%p(isd:ied,k_size))
    allocate(k0a_Rv(J)%p(isd:ied,k_size))
    allocate(k0b_Lv(J)%p(isd:ied,k_size))
    allocate(k0b_Rv(J)%p(isd:ied,k_size))
  enddo

!$OMP parallel do default(none) shared(is,ie,js,je,G,num_srt,rho_srt,k0b_Lv,k0b_Rv, &
!$OMP                                  k0_srt,k0a_Lv,k0a_Rv,deep_wt_Lv,deep_wt_Rv,  &
!$OMP                                  h_srt,nkmb,nPv,hP_Lv,hP_Rv)                  &
!$OMP                          private(h_demand_L,h_used_L,h_demand_R,h_used_R,     &
!$OMP                                  kR,kL,nP,rho_pair,kbs_Lp,kbs_Rp,rho_a,rho_b, &
!$OMP                                  wt_b,left_set,right_set,h_supply_frac_R,     &
!$OMP                                  h_supply_frac_L)
  do J=js-1,je ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.5) then
    ! Set up the pairings for fluxes through the meridional faces.

    do k=1,num_srt(i,j)   ; h_demand_L(k) = 0.0 ; h_used_L(k) = 0.0 ; enddo
    do k=1,num_srt(i,j+1) ; h_demand_R(k) = 0.0 ; h_used_R(k) = 0.0 ; enddo

    ! First merge the left and right lists into a single, sorted list.

    !   Discard any layers that are lighter than the lightest in the other
    ! column.  They can only participate in mixing as the lighter part of a
    ! pair of points.
    if (rho_srt(i,1,j) < rho_srt(i,1,j+1)) then
      kR = 1
      do kL=2,num_srt(i,j) ; if (rho_srt(i,kL,j) >= rho_srt(i,1,j+1)) exit ; enddo
    elseif (rho_srt(i,1,j+1) < rho_srt(i,1,j)) then
      kL = 1
      do kR=2,num_srt(i,j+1) ; if (rho_srt(i,kR,j+1) >= rho_srt(i,1,j)) exit ; enddo
    else
      kL = 1 ; kR = 1
    endif
    nP = 0
    do ! Loop to accumulate pairs of columns.
      if ((kL > num_srt(i,j)) .or. (kR > num_srt(i,j+1))) exit

      if (rho_srt(i,kL,j) > rho_srt(i,kR,j+1)) then
      ! The right point is lighter and defines the density for this trio.
        nP = nP+1 ; k = nP
        rho_pair = rho_srt(i,kR,j+1)

        k0b_Lv(J)%p(i,k) = k0_srt(i,kL,j)   ; k0b_Rv(J)%p(i,k) = k0_srt(i,kR,j+1)
        k0a_Lv(J)%p(i,k) = k0_srt(i,kL-1,j) ; k0a_Rv(J)%p(i,k) = k0b_Rv(J)%p(i,k)
        kbs_Lp(k) = kL ; kbs_Rp(k) = kR

        rho_a = rho_srt(i,kL-1,j) ; rho_b = rho_srt(i,kL,j)
        wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
          wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
        deep_wt_Lv(J)%p(i,k) = wt_b ; deep_wt_Rv(J)%p(i,k) = 1.0

        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i,kR,j+1) * wt_b
        h_demand_L(kL-1) = h_demand_L(kL-1) + 0.5*h_srt(i,kR,j+1) * (1.0-wt_b)

        kR = kR+1 ; left_set(k) = .false. ; right_set(k) = .true.
      elseif (rho_srt(i,kL,j) < rho_srt(i,kR,j+1)) then
      ! The left point is lighter and defines the density for this trio.
        nP = nP+1 ; k = nP
        rho_pair = rho_srt(i,kL,j)
        k0b_Lv(J)%p(i,k) = k0_srt(i,kL,j) ; k0b_Rv(J)%p(i,k) = k0_srt(i,kR,j+1)
        k0a_Lv(J)%p(i,k) = k0b_Lv(J)%p(i,k) ; k0a_Rv(J)%p(i,k) = k0_srt(i,kR-1,j+1)

        kbs_Lp(k) = kL ; kbs_Rp(k) = kR

        rho_a = rho_srt(i,kR-1,j+1) ; rho_b = rho_srt(i,kR,j+1)
        wt_b = 1.0 ; if (abs(rho_a - rho_b) > abs(rho_pair - rho_a)) &
          wt_b = (rho_pair - rho_a) / (rho_b - rho_a)
        deep_wt_Lv(J)%p(i,k) = 1.0 ; deep_wt_Rv(J)%p(i,k) = wt_b

        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j) * wt_b
        h_demand_R(kR-1) = h_demand_R(kR-1) + 0.5*h_srt(i,kL,j) * (1.0-wt_b)

        kL = kL+1 ; left_set(k) = .true. ; right_set(k) = .false.
      elseif ((k0_srt(i,kL,j) <= nkmb) .or. (k0_srt(i,kR,j+1) <= nkmb)) then
        ! The densities are exactly equal and one layer is above the interior.
        nP = nP+1 ; k = nP
        k0b_Lv(J)%p(i,k) = k0_srt(i,kL,j) ; k0b_Rv(J)%p(i,k) = k0_srt(i,kR,j+1)
        k0a_Lv(J)%p(i,k) = k0b_Lv(J)%p(i,k)  ; k0a_Rv(J)%p(i,k) = k0b_Rv(J)%p(i,k)
        kbs_Lp(k) = kL ; kbs_Rp(k) = kR
        deep_wt_Lv(J)%p(i,k) = 1.0 ; deep_wt_Rv(J)%p(i,k) = 1.0

        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i,kR,j+1)
        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)

        kL = kL+1 ; kR = kR+1 ; left_set(k) = .true. ; right_set(k) = .true.
      else ! The densities are exactly equal and in the interior.
        ! Mixing in this case has already occurred, so accumulate the thickness
        ! demanded for that mixing and skip onward.
        h_demand_L(kL) = h_demand_L(kL) + 0.5*h_srt(i,kR,j+1)
        h_demand_R(kR) = h_demand_R(kR) + 0.5*h_srt(i,kL,j)

        kL = kL+1 ; kR = kR+1
      endif
    enddo ! Loop to accumulate pairs of columns.
    nPv(i,J) = nP ! This is the number of active pairings.

    ! Determine what fraction of the thickness "demand" can be supplied.
    do k=1,num_srt(i,j+1)
      h_supply_frac_R(k) = 1.0
      if (h_demand_R(k) > 0.5*h_srt(i,k,j+1)) &
        h_supply_frac_R(k) = 0.5*h_srt(i,k,j+1) / h_demand_R(k)
    enddo
    do k=1,num_srt(i,j)
      h_supply_frac_L(k) = 1.0
      if (h_demand_L(k) > 0.5*h_srt(i,k,j)) &
        h_supply_frac_L(k) = 0.5*h_srt(i,k,j) / h_demand_L(k)
    enddo

    !  Distribute the "exported" thicknesses proportionately.
    do k=1,nPv(i,J)
      kL = kbs_Lp(k) ; kR = kbs_Rp(k)
      hP_Lv(J)%p(i,k) = 0.0 ; hP_Rv(J)%p(i,k) = 0.0
      if (left_set(k)) then ! Add the contributing thicknesses on the right.
        if (deep_wt_Rv(J)%p(i,k) < 1.0) then
          hP_Rv(J)%p(i,k) = 0.5*h_srt(i,kL,j) * min(h_supply_frac_R(kR), h_supply_frac_R(kR-1))
          wt_b = deep_wt_Rv(J)%p(i,k)
          h_used_R(kR-1) = h_used_R(kR-1) + (1.0 - wt_b) * hP_Rv(J)%p(i,k)
          h_used_R(kR) = h_used_R(kR) + wt_b * hP_Rv(J)%p(i,k)
        else
          hP_Rv(J)%p(i,k) = 0.5*h_srt(i,kL,j) * h_supply_frac_R(kR)
          h_used_R(kR) = h_used_R(kR) + hP_Rv(J)%p(i,k)
        endif
      endif
      if (right_set(k)) then ! Add the contributing thicknesses on the left.
        if (deep_wt_Lv(J)%p(i,k) < 1.0) then
          hP_Lv(J)%p(i,k) = 0.5*h_srt(i,kR,j+1) * min(h_supply_frac_L(kL), h_supply_frac_L(kL-1))
          wt_b = deep_wt_Lv(J)%p(i,k)
          h_used_L(kL-1) = h_used_L(kL-1) + (1.0 - wt_b) * hP_Lv(J)%p(i,k)
          h_used_L(kL) = h_used_L(kL) + wt_b * hP_Lv(J)%p(i,k)
        else
          hP_Lv(J)%p(i,k) = 0.5*h_srt(i,kR,j+1) * h_supply_frac_L(kL)
          h_used_L(kL) = h_used_L(kL) + hP_Lv(J)%p(i,k)
        endif
      endif
    enddo

    !   The left-over thickness (at least half the layer thickness) is now
    ! added to the thicknesses of the importing columns.
    do k=1,nPv(i,J)
      if (left_set(k)) hP_Lv(J)%p(i,k) = hP_Lv(J)%p(i,k) + &
                            (h_srt(i,kbs_Lp(k),j) - h_used_L(kbs_Lp(k)))
      if (right_set(k)) hP_Rv(J)%p(i,k) = hP_Rv(J)%p(i,k) + &
                             (h_srt(i,kbs_Rp(k),j+1) - h_used_R(kbs_Rp(k)))
    enddo


  endif ; enddo ; enddo ! i- & j- loops over meridional faces.

! The tracer-specific calculations start here.

  ! Zero out tracer tendencies.
  do k=1,PEmax_kRho ; do j=js-1,je+1 ; do i=is-1,ie+1
    tr_flux_conv(i,j,k) = 0.0
  enddo ; enddo ; enddo

  do itt=1,max_itt

    if (itt > 1) then ! The halos have already been filled if itt==1.
      call cpu_clock_begin(id_clock_pass)
      call do_group_pass(CS%pass_t, G%Domain)
      call cpu_clock_end(id_clock_pass)
    endif

    do m=1,ntr
!$OMP parallel do default(none) shared(is,ie,js,je,G,Tr,nkmb,nPu,m,max_kRho,nz,h,h_exclude, &
!$OMP                                  k0b_Lu,k0b_Ru,deep_wt_Lu,k0a_Lu,deep_wt_Ru,k0a_Ru,   &
!$OMP                                  hP_Lu,hP_Ru,I_maxitt,khdt_epi_x,tr_flux_conv,Idt) &
!$OMP                          private(Tr_min_face,Tr_max_face,kLa,kLb,kRa,kRb,Tr_La, &
!$OMP                                     Tr_Lb,Tr_Ra,Tr_Rb,Tr_av_L,wt_b,Tr_av_R,h_L,h_R, &
!$OMP                                     Tr_flux,Tr_adj_vert,wt_a,vol)
      do j=js,je ; do I=is-1,ie ; if (G%mask2dCu(I,j) > 0.5) then
        ! Determine the fluxes through the zonal faces.

        ! Find the acceptable range of tracer concentration around this face.
        if (nPu(I,j) >= 1) then
          Tr_min_face = min(Tr(m)%t(i,j,1), Tr(m)%t(i+1,j,1))
          Tr_max_face = max(Tr(m)%t(i,j,1), Tr(m)%t(i+1,j,1))
          do k=2,nkmb
            Tr_min_face = min(Tr_min_face, Tr(m)%t(i,j,k), Tr(m)%t(i+1,j,k))
            Tr_max_face = max(Tr_max_face, Tr(m)%t(i,j,k), Tr(m)%t(i+1,j,k))
          enddo

          ! Include the next two layers denser than the densest buffer layer.
          kLa = nkmb+1 ; if (max_kRho(i,j) < nz+1) kLa = max_kRho(i,j)
          kLb = kLa ; if (max_kRho(i,j) < nz) kLb = max_kRho(i,j)+1
          kRa = nkmb+1 ; if (max_kRho(i+1,j) < nz+1) kRa = max_kRho(i+1,j)
          kRb = kRa ; if (max_kRho(i+1,j) < nz) kRb = max_kRho(i+1,j)+1
          Tr_La = Tr_min_face ; Tr_Lb = Tr_La ; Tr_Ra = Tr_La ; Tr_Rb = Tr_La
          if (h(i,j,kLa) > h_exclude) Tr_La = Tr(m)%t(i,j,kLa)
          if (h(i,j,kLb) > h_exclude) Tr_La = Tr(m)%t(i,j,kLb)
          if (h(i+1,j,kRa) > h_exclude) Tr_Ra = Tr(m)%t(i+1,j,kRa)
          if (h(i+1,j,kRb) > h_exclude) Tr_Rb = Tr(m)%t(i+1,j,kRb)
          Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)

          ! Include all points in diffusive pairings at this face.
          do k=1,nPu(I,j)
            Tr_Lb = Tr(m)%t(i,j,k0b_Lu(j)%p(I,k))
            Tr_Rb = Tr(m)%t(i+1,j,k0b_Ru(j)%p(I,k))
            Tr_La = Tr_Lb ; Tr_Ra = Tr_Rb
            if (deep_wt_Lu(j)%p(I,k) < 1.0) Tr_La = Tr(m)%t(i,j,k0a_Lu(j)%p(I,k))
            if (deep_wt_Ru(j)%p(I,k) < 1.0) Tr_Ra = Tr(m)%t(i+1,j,k0a_Ru(j)%p(I,k))
            Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
            Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          enddo
        endif

        do k=1,nPu(I,j)
          kLb = k0b_Lu(j)%p(I,k) ; Tr_Lb = Tr(m)%t(i,j,kLb) ; Tr_av_L = Tr_Lb
          if (deep_wt_Lu(j)%p(I,k) < 1.0) then
            kLa = k0a_Lu(j)%p(I,k) ; Tr_La = Tr(m)%t(i,j,kLa)
            wt_b = deep_wt_Lu(j)%p(I,k)
            Tr_av_L = wt_b*Tr_Lb + (1.0-wt_b)*Tr_La
          endif

          kRb = k0b_Ru(j)%p(I,k) ; Tr_Rb = Tr(m)%t(i+1,j,kRb) ; Tr_av_R = Tr_Rb
          if (deep_wt_Ru(j)%p(I,k) < 1.0) then
            kRa = k0a_Ru(j)%p(I,k) ; Tr_Ra = Tr(m)%t(i+1,j,kRa)
            wt_b = deep_wt_Ru(j)%p(I,k)
            Tr_av_R = wt_b*Tr_Rb + (1.0-wt_b)*Tr_Ra
          endif

          h_L = hP_Lu(j)%p(I,k) ; h_R = hP_Ru(j)%p(I,k)
          Tr_flux = I_maxitt * khdt_epi_x(I,j) * (Tr_av_L - Tr_av_R) * &
            ((2.0 * h_L * h_R) / (h_L + h_R))


          if (deep_wt_Lu(j)%p(I,k) >= 1.0) then
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - Tr_flux
          else
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Lu(j)%p(I,k) ; wt_a = 1.0 - wt_b
            vol = hP_Lu(j)%p(I,k) * G%areaT(i,j)

            !   Ensure that the tracer flux does not drive the tracer values
            ! outside of the range Tr_min_face <= Tr <= Tr_max_face, or if it
            ! does that the concentration in both contributing peices exceed
            ! this range equally. With downgradient fluxes and the initial tracer
            ! concentrations determining the valid range, the latter condition
            ! only enters for large values of the effective diffusive CFL number.
            if (Tr_flux > 0.0) then
              if (Tr_La < Tr_Lb) then ; if (vol*(Tr_La-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = -wt_a * min(Tr_flux - vol * (Tr_La-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Lb - Tr_La))
              else ; if (vol*(Tr_Lb-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = wt_b * min(Tr_flux - vol * (Tr_Lb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_La - Tr_Lb))
              endif
            elseif (Tr_flux < 0.0) then
              if (Tr_La > Tr_Lb) then ; if (vol * (Tr_max_face-Tr_La) < -Tr_flux) &
                Tr_adj_vert = wt_a * min(-Tr_flux - vol * (Tr_max_face-Tr_La), &
                                         (vol*wt_b) * (Tr_La - Tr_Lb))
              else ; if (vol*(Tr_max_face-Tr_Lb) < -Tr_flux) &
                Tr_adj_vert = -wt_b * min(-Tr_flux - vol * (Tr_max_face-Tr_Lb), &
                                          (vol*wt_a)*(Tr_Lb - Tr_La))
              endif
            endif

            tr_flux_conv(i,j,kLa) = tr_flux_conv(i,j,kLa) - (wt_a*Tr_flux + Tr_adj_vert)
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - (wt_b*Tr_flux - Tr_adj_vert)
          endif

          if (deep_wt_Ru(j)%p(I,k) >= 1.0) then
            tr_flux_conv(i+1,j,kRb) = tr_flux_conv(i+1,j,kRb) + Tr_flux
          else
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Ru(j)%p(I,k) ; wt_a = 1.0 - wt_b
            vol = hP_Ru(j)%p(I,k) * G%areaT(i+1,j)

            !   Ensure that the tracer flux does not drive the tracer values
            ! outside of the range Tr_min_face <= Tr <= Tr_max_face, or if it
            ! does that the concentration in both contributing peices exceed
            ! this range equally. With downgradient fluxes and the initial tracer
            ! concentrations determining the valid range, the latter condition
            ! only enters for large values of the effective diffusive CFL number.
            if (Tr_flux < 0.0) then
              if (Tr_Ra < Tr_Rb) then ; if (vol * (Tr_Ra-Tr_min_face) < -Tr_flux) &
                Tr_adj_vert = -wt_a * min(-Tr_flux - vol * (Tr_Ra-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Rb - Tr_Ra))
              else ; if (vol*(Tr_Rb-Tr_min_face) < (-Tr_flux)) &
                Tr_adj_vert = wt_b * min(-Tr_flux - vol * (Tr_Rb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_Ra - Tr_Rb))
              endif
            elseif (Tr_flux > 0.0) then
              if (Tr_Ra > Tr_Rb) then ; if (vol * (Tr_max_face-Tr_Ra) < Tr_flux) &
                Tr_adj_vert = wt_a * min(Tr_flux - vol * (Tr_max_face-Tr_Ra), &
                                         (vol*wt_b) * (Tr_Ra - Tr_Rb))
              else ; if (vol*(Tr_max_face-Tr_Rb) < Tr_flux) &
                Tr_adj_vert = -wt_b * min(Tr_flux - vol * (Tr_max_face-Tr_Rb), &
                                          (vol*wt_a)*(Tr_Rb - Tr_Ra))
              endif
            endif

            tr_flux_conv(i+1,j,kRa) = tr_flux_conv(i+1,j,kRa) + &
                                            (wt_a*Tr_flux - Tr_adj_vert)
            tr_flux_conv(i+1,j,kRb) = tr_flux_conv(i+1,j,kRb) + &
                                            (wt_b*Tr_flux + Tr_adj_vert)
          endif
          if (associated(Tr(m)%df2d_x)) &
            Tr(m)%df2d_x(I,j) = Tr(m)%df2d_x(I,j) + Tr_flux * Idt
        enddo ! Loop over pairings at faces.
      endif ; enddo ; enddo ! i- & j- loops over zonal faces.

!$OMP parallel do default(none) shared(is,ie,js,je,G,Tr,nkmb,nPv,m,max_kRho,nz,h,h_exclude, &
!$OMP                                  k0b_Lv,k0b_Rv,deep_wt_Lv,k0a_Lv,deep_wt_Rv,k0a_Rv,   &
!$OMP                                  hP_Lv,hP_Rv,I_maxitt,khdt_epi_y,Tr_flux_3d,          &
!$OMP                                  Tr_adj_vert_L,Tr_adj_vert_R,Idt)                     &
!$OMP                          private(Tr_min_face,Tr_max_face,kLa,kLb,kRa,kRb,             &
!$OMP                                  Tr_La,Tr_Lb,Tr_Ra,Tr_Rb,Tr_av_L,wt_b,Tr_av_R,        &
!$OMP                                  h_L,h_R,Tr_flux,Tr_adj_vert,wt_a,vol)
      do J=js-1,je ; do i=is,ie ; if (G%mask2dCv(i,J) > 0.5) then
        ! Determine the fluxes through the meridional faces.

        ! Find the acceptable range of tracer concentration around this face.
        if (nPv(i,J) >= 1) then
          Tr_min_face = min(Tr(m)%t(i,j,1), Tr(m)%t(i,j+1,1))
          Tr_max_face = max(Tr(m)%t(i,j,1), Tr(m)%t(i,j+1,1))
          do k=2,nkmb
            Tr_min_face = min(Tr_min_face, Tr(m)%t(i,j,k), Tr(m)%t(i,j+1,k))
            Tr_max_face = max(Tr_max_face, Tr(m)%t(i,j,k), Tr(m)%t(i,j+1,k))
          enddo

          ! Include the next two layers denser than the densest buffer layer.
          kLa = nkmb+1 ; if (max_kRho(i,j) < nz+1) kLa = max_kRho(i,j)
          kLb = kLa ; if (max_kRho(i,j) < nz) kLb = max_kRho(i,j)+1
          kRa = nkmb+1 ; if (max_kRho(i,j+1) < nz+1) kRa = max_kRho(i,j+1)
          kRb = kRa ; if (max_kRho(i,j+1) < nz) kRb = max_kRho(i,j+1)+1
          Tr_La = Tr_min_face ; Tr_Lb = Tr_La ; Tr_Ra = Tr_La ; Tr_Rb = Tr_La
          if (h(i,j,kLa) > h_exclude) Tr_La = Tr(m)%t(i,j,kLa)
          if (h(i,j,kLb) > h_exclude) Tr_La = Tr(m)%t(i,j,kLb)
          if (h(i,j+1,kRa) > h_exclude) Tr_Ra = Tr(m)%t(i,j+1,kRa)
          if (h(i,j+1,kRb) > h_exclude) Tr_Rb = Tr(m)%t(i,j+1,kRb)
          Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)

          ! Include all points in diffusive pairings at this face.
          do k=1,nPv(i,J)
            Tr_Lb = Tr(m)%t(i,j,k0b_Lv(J)%p(i,k)) ; Tr_Rb = Tr(m)%t(i,j+1,k0b_Rv(J)%p(i,k))
            Tr_La = Tr_Lb ; Tr_Ra = Tr_Rb
            if (deep_wt_Lv(J)%p(i,k) < 1.0) Tr_La = Tr(m)%t(i,j,k0a_Lv(J)%p(i,k))
            if (deep_wt_Rv(J)%p(i,k) < 1.0) Tr_Ra = Tr(m)%t(i,j+1,k0a_Rv(J)%p(i,k))
            Tr_min_face = min(Tr_min_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
            Tr_max_face = max(Tr_max_face, Tr_La, Tr_Lb, Tr_Ra, Tr_Rb)
          enddo
        endif

        do k=1,nPv(i,J)
          kLb = k0b_Lv(J)%p(i,k) ; Tr_Lb = Tr(m)%t(i,j,kLb) ; Tr_av_L = Tr_Lb
          if (deep_wt_Lv(J)%p(i,k) < 1.0) then
            kLa = k0a_Lv(J)%p(i,k) ; Tr_La = Tr(m)%t(i,j,kLa)
            wt_b = deep_wt_Lv(J)%p(i,k)
            Tr_av_L = wt_b * Tr_Lb + (1.0-wt_b) * Tr_La
          endif

          kRb = k0b_Rv(J)%p(i,k) ; Tr_Rb = Tr(m)%t(i,j+1,kRb) ; Tr_av_R = Tr_Rb
          if (deep_wt_Rv(J)%p(i,k) < 1.0) then
            kRa = k0a_Rv(J)%p(i,k) ; Tr_Ra = Tr(m)%t(i,j+1,kRa)
            wt_b = deep_wt_Rv(J)%p(i,k)
            Tr_av_R = wt_b * Tr_Rb + (1.0-wt_b) * Tr_Ra
          endif

          h_L = hP_Lv(J)%p(i,k) ; h_R = hP_Rv(J)%p(i,k)
          Tr_flux = I_maxitt * ((2.0 * h_L * h_R) / (h_L + h_R)) * &
                    khdt_epi_y(i,J) * (Tr_av_L - Tr_av_R)
          Tr_flux_3d(i,j,k) = Tr_flux

          if (deep_wt_Lv(J)%p(i,k) < 1.0) then
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Lv(J)%p(i,k) ; wt_a = 1.0 - wt_b
            vol = hP_Lv(J)%p(i,k) * G%areaT(i,j)

            !   Ensure that the tracer flux does not drive the tracer values
            ! outside of the range Tr_min_face <= Tr <= Tr_max_face.
            if (Tr_flux > 0.0) then
              if (Tr_La < Tr_Lb) then ; if (vol * (Tr_La-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = -wt_a * min(Tr_flux - vol * (Tr_La-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Lb - Tr_La))
              else ; if (vol*(Tr_Lb-Tr_min_face) < Tr_flux) &
                Tr_adj_vert = wt_b * min(Tr_flux - vol * (Tr_Lb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_La - Tr_Lb))
              endif
            elseif (Tr_flux < 0.0) then
              if (Tr_La > Tr_Lb) then ; if (vol * (Tr_max_face-Tr_La) < -Tr_flux) &
                Tr_adj_vert = wt_a * min(-Tr_flux - vol * (Tr_max_face-Tr_La), &
                                         (vol*wt_b) * (Tr_La - Tr_Lb))
              else ; if (vol*(Tr_max_face-Tr_Lb) < -Tr_flux) &
                Tr_adj_vert = -wt_b * min(-Tr_flux - vol * (Tr_max_face-Tr_Lb), &
                                          (vol*wt_a)*(Tr_Lb - Tr_La))
              endif
            endif
            Tr_adj_vert_L(i,j,k) = Tr_adj_vert
          endif

          if (deep_wt_Rv(J)%p(i,k) < 1.0) then
            Tr_adj_vert = 0.0
            wt_b = deep_wt_Rv(J)%p(i,k) ; wt_a = 1.0 - wt_b
            vol = hP_Rv(J)%p(i,k) * G%areaT(i,j+1)

            !   Ensure that the tracer flux does not drive the tracer values
            ! outside of the range Tr_min_face <= Tr <= Tr_max_face.
            if (Tr_flux < 0.0) then
              if (Tr_Ra < Tr_Rb) then ; if (vol * (Tr_Ra-Tr_min_face) < -Tr_flux) &
                Tr_adj_vert = -wt_a * min(-Tr_flux - vol * (Tr_Ra-Tr_min_face), &
                                          (vol*wt_b) * (Tr_Rb - Tr_Ra))
              else ; if (vol*(Tr_Rb-Tr_min_face) < (-Tr_flux)) &
                Tr_adj_vert = wt_b * min(-Tr_flux - vol * (Tr_Rb-Tr_min_face), &
                                         (vol*wt_a) * (Tr_Ra - Tr_Rb))
              endif
            elseif (Tr_flux > 0.0) then
              if (Tr_Ra > Tr_Rb) then ; if (vol * (Tr_max_face-Tr_Ra) < Tr_flux) &
                Tr_adj_vert = wt_a * min(Tr_flux - vol * (Tr_max_face-Tr_Ra), &
                                         (vol*wt_b) * (Tr_Ra - Tr_Rb))
              else ; if (vol*(Tr_max_face-Tr_Rb) < Tr_flux) &
                Tr_adj_vert = -wt_b * min(Tr_flux - vol * (Tr_max_face-Tr_Rb), &
                                          (vol*wt_a)*(Tr_Rb - Tr_Ra))
              endif
            endif
            Tr_adj_vert_R(i,j,k) = Tr_adj_vert
          endif
          if (associated(Tr(m)%df2d_y)) &
            Tr(m)%df2d_y(i,J) = Tr(m)%df2d_y(i,J) + Tr_flux * Idt
        enddo ! Loop over pairings at faces.
      endif ; enddo ; enddo ! i- & j- loops over meridional faces.
!$OMP parallel do default(none) shared(is,ie,js,je,G,nPv,k0b_Lv,k0b_Rv,deep_wt_Lv,  &
!$OMP                                  tr_flux_conv,Tr_flux_3d,k0a_Lv,Tr_adj_vert_L,&
!$OMP                                  deep_wt_Rv,k0a_Rv,Tr_adj_vert_R) &
!$OMP                          private(kLa,kLb,kRa,kRb,wt_b,wt_a)
      do i=is,ie ; do J=js-1,je ; if (G%mask2dCv(i,J) > 0.5) then
        do k=1,nPv(i,J)
          kLb = k0b_Lv(J)%p(i,k); kRb = k0b_Rv(J)%p(i,k)
          if (deep_wt_Lv(J)%p(i,k) >= 1.0) then
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - Tr_flux_3d(i,j,k)
          else
            kLa = k0a_Lv(J)%p(i,k)
            wt_b = deep_wt_Lv(J)%p(i,k) ; wt_a = 1.0 - wt_b
            tr_flux_conv(i,j,kLa) = tr_flux_conv(i,j,kLa) - (wt_a*Tr_flux_3d(i,j,k) + Tr_adj_vert_L(i,j,k))
            tr_flux_conv(i,j,kLb) = tr_flux_conv(i,j,kLb) - (wt_b*Tr_flux_3d(i,j,k) - Tr_adj_vert_L(i,j,k))
          endif
          if (deep_wt_Rv(J)%p(i,k) >= 1.0) then
            tr_flux_conv(i,j+1,kRb) = tr_flux_conv(i,j+1,kRb) + tr_flux_3d(i,j,k)
          else
            kRa = k0a_Rv(J)%p(i,k)
            wt_b = deep_wt_Rv(J)%p(i,k) ; wt_a = 1.0 - wt_b
            tr_flux_conv(i,j+1,kRa) = tr_flux_conv(i,j+1,kRa) + &
                                            (wt_a*Tr_flux_3d(i,j,k) - Tr_adj_vert_R(i,j,k))
            tr_flux_conv(i,j+1,kRb) = tr_flux_conv(i,j+1,kRb) + &
                                            (wt_b*Tr_flux_3d(i,j,k) + Tr_adj_vert_R(i,j,k))
          endif
        enddo
      endif ; enddo ; enddo
!$OMP parallel do default(none) shared(PEmax_kRho,is,ie,js,je,G,h,Tr,tr_flux_conv,m)
      do k=1,PEmax_kRho ; do j=js,je ; do i=is,ie
        if ((G%mask2dT(i,j) > 0.5) .and. (h(i,j,k) > 0.0)) then
          Tr(m)%t(i,j,k) = Tr(m)%t(i,j,k) + tr_flux_conv(i,j,k) / &
                                            (h(i,j,k)*G%areaT(i,j))
          tr_flux_conv(i,j,k) = 0.0
        endif
      enddo ; enddo ; enddo

    enddo ! Loop over tracers
  enddo ! Loop over iterations

  do j=js,je
    deallocate(deep_wt_Lu(j)%p) ; deallocate(deep_wt_Ru(j)%p)
    deallocate(Hp_Lu(j)%p)  ; deallocate(Hp_Ru(j)%p)
    deallocate(k0a_Lu(j)%p) ; deallocate(k0a_Ru(j)%p)
    deallocate(k0b_Lu(j)%p) ; deallocate(k0b_Ru(j)%p)
  enddo

  do J=js-1,je
    deallocate(deep_wt_Lv(J)%p) ; deallocate(deep_wt_Rv(J)%p)
    deallocate(Hp_Lv(J)%p)  ; deallocate(Hp_Rv(J)%p)
    deallocate(k0a_Lv(J)%p) ; deallocate(k0a_Rv(J)%p)
    deallocate(k0b_Lv(J)%p) ; deallocate(k0b_Rv(J)%p)
  enddo

end subroutine tracer_epipycnal_ML_diff


!> Initialize lateral tracer diffusion module
subroutine tracer_hor_diff_init(Time, G, param_file, diag, CS, CSnd)
  type(time_type), target,    intent(in)    :: Time       !< current model time
  type(ocean_grid_type),      intent(in)    :: G          !< ocean grid structure
  type(diag_ctrl), target,    intent(inout) :: diag       !< diagnostic control
  type(param_file_type),      intent(in)    :: param_file !< parameter file
  type(tracer_hor_diff_CS),   pointer       :: CS         !< horz diffusion control structure
  type(neutral_diffusion_CS), pointer       :: CSnd       !< pointer to neutral diffusion CS

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_tracer_hor_diff" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (associated(CS)) then
    call MOM_error(WARNING, "tracer_hor_diff_init called with associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag
  CS%show_call_tree = callTree_showQuery()

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "KHTR", CS%KhTr, &
                 "The background along-isopycnal tracer diffusivity.", &
                 units="m2 s-1", default=0.0)
  call get_param(param_file, mod, "KHTR_SLOPE_CFF", CS%KhTr_Slope_Cff, &
                 "The scaling coefficient for along-isopycnal tracer \n"//&
                 "diffusivity using a shear-based (Visbeck-like) \n"//&
                 "parameterization.  A non-zero value enables this param.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mod, "KHTR_MIN", CS%KhTr_Min, &
                 "The minimum along-isopycnal tracer diffusivity.", &
                 units="m2 s-1", default=0.0)
  call get_param(param_file, mod, "KHTR_MAX", CS%KhTr_Max, &
                 "The maximum along-isopycnal tracer diffusivity.", &
                 units="m2 s-1", default=0.0)
  call get_param(param_file, mod, "KHTR_PASSIVITY_COEFF", CS%KhTr_passivity_coeff, &
               "The coefficient that scales deformation radius over \n"//&
               "grid-spacing in passivity, where passiviity is the ratio \n"//&
               "between along isopycnal mxiing of tracers to thickness mixing. \n"//&
               "A non-zero value enables this parameterization.", &
               units="nondim", default=0.0)
  call get_param(param_file, mod, "KHTR_PASSIVITY_MIN", CS%KhTr_passivity_min, &
               "The minimum passivity which is the ratio between \n"//&
               "along isopycnal mxiing of tracers to thickness mixing. \n", &
               units="nondim", default=0.5)
  call get_param(param_file, mod, "DT", CS%dt, fail_if_missing=.true., &
          desc="The (baroclinic) dynamics time step.", units="s")


  call get_param(param_file, mod, "DIFFUSE_ML_TO_INTERIOR", CS%Diffuse_ML_interior, &
                 "If true, enable epipycnal mixing between the surface \n"//&
                 "boundary layer and the interior.", default=.false.)
  call get_param(param_file, mod, "CHECK_DIFFUSIVE_CFL", CS%check_diffusive_CFL, &
                 "If true, use enough iterations the diffusion to ensure \n"//&
                 "that the diffusive equivalent of the CFL limit is not \n"//&
                 "violated.  If false, always use 1 iteration.", default=.false.)
  CS%ML_KhTR_scale = 1.0
  if (CS%Diffuse_ML_interior) then
    call get_param(param_file, mod, "ML_KHTR_SCALE", CS%ML_KhTR_scale, &
                 "With Diffuse_ML_interior, the ratio of the truly \n"//&
                 "horizontal diffusivity in the mixed layer to the \n"//&
                 "epipycnal diffusivity.  The valid range is 0 to 1.", &
                 units="nondim", default=1.0)
  endif

  CS%use_neutral_diffusion = neutral_diffusion_init(Time, G, param_file, diag, CS%neutral_diffusion_CSp)
  CSnd => CS%neutral_diffusion_CSp
  if (CS%use_neutral_diffusion .and. CS%Diffuse_ML_interior) call MOM_error(FATAL, "MOM_tracer_hor_diff: "// &
       "USE_NEUTRAL_DIFFUSION and DIFFUSE_ML_TO_INTERIOR are mutually exclusive!")

  call get_param(param_file, mod, "DEBUG", CS%debug, default=.false.)

  id_clock_diffuse = cpu_clock_id('(Ocean diffuse tracer)',          grain=CLOCK_MODULE)
  id_clock_epimix  = cpu_clock_id('(Ocean epipycnal diffuse tracer)',grain=CLOCK_MODULE)
  id_clock_pass    = cpu_clock_id('(Ocean tracer halo updates)',     grain=CLOCK_ROUTINE)
  id_clock_sync    = cpu_clock_id('(Ocean tracer global synch)',     grain=CLOCK_ROUTINE)

  CS%id_KhTr_u = -1
  CS%id_KhTr_v = -1
  CS%id_KhTr_h = -1
  CS%id_CFL    = -1

  CS%id_KhTr_u = register_diag_field('ocean_model', 'KHTR_u', diag%axesCu1, Time, &
     'Epipycnal tracer diffusivity at zonal faces of tracer cell', 'meter2 second-1')
  CS%id_KhTr_v = register_diag_field('ocean_model', 'KHTR_v', diag%axesCv1, Time, &
     'Epipycnal tracer diffusivity at meridional faces of tracer cell', 'meter2 second-1')
  CS%id_KhTr_h = register_diag_field('ocean_model', 'KHTR_h', diag%axesT1, Time,&
     'Epipycnal tracer diffusivity at tracer cell center', 'meter2 second-1',   &
     cmor_field_name='diftrelo', cmor_units='m2 sec-1',                         &
     cmor_standard_name= 'ocean_tracer_epineutral_laplacian_diffusivity',       &
     cmor_long_name = 'Ocean Tracer Epineutral Laplacian Diffusivity')

  CS%id_khdt_x = register_diag_field('ocean_model', 'KHDT_x', diag%axesCu1, Time, &
     'Epipycnal tracer diffusivity operator at zonal faces of tracer cell', 'meter2')
  CS%id_khdt_y = register_diag_field('ocean_model', 'KHDT_y', diag%axesCv1, Time, &
     'Epipycnal tracer diffusivity operator at meridional faces of tracer cell', 'meter2')
  if (CS%check_diffusive_CFL) then
    CS%id_CFL = register_diag_field('ocean_model', 'CFL_lateral_diff', diag%axesT1, Time,&
       'Grid CFL number for lateral/neutral tracer diffusion', 'dimensionless')
  endif


end subroutine tracer_hor_diff_init

subroutine tracer_hor_diff_end(CS)
  type(tracer_hor_diff_CS), pointer :: CS

  call neutral_diffusion_end(CS%neutral_diffusion_CSp)
  if (associated(CS)) deallocate(CS)

end subroutine tracer_hor_diff_end


!> \namespace mom_tracer_hor_diff
!!
!! \section section_intro Introduction to the module
!!
!!    This module contains subroutines that handle horizontal
!!  diffusion (i.e., isoneutral or along layer) of tracers.
!!
!!    Each of the tracers are subject to Fickian along-coordinate
!!  diffusion if Khtr is defined and positive.  The tracer diffusion
!!  can use a suitable number of iterations to guarantee stability
!!  with an arbitrarily large time step.
!!
!!  \section section_gridlayout MOM grid layout
!!
!!  A small fragment of the grid is shown below:
!!
!! \verbatim
!!    j+1  x ^ x ^ x
!!
!!    j+1  > o > o >
!!
!!    j    x ^ x ^ x
!!
!!    j    > o > o >
!!
!!    j-1  x ^ x ^ x
!!
!!        i-1  i  i+1
!!
!!           i  i+1
!!
!! \endverbatim
!!
!!  Fields at each point
!!  * x =  q, CoriolisBu
!!  * ^ =  v, PFv, CAv, vh, diffv, tauy, vbt, vhtr
!!  * > =  u, PFu, CAu, uh, diffu, taux, ubt, uhtr
!!  * o =  h, bathyT, eta, T, S, tr
!!
!!  The boundaries always run through q grid points (x).
!!

end module MOM_tracer_hor_diff
