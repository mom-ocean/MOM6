!> Shear-dependent mixing following Jackson et al. 2008.
module MOM_kappa_shear

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_debugging, only : hchksum, Bchksum
use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public Calculate_kappa_shear, Calc_kappa_shear_vertex, kappa_shear_init
public kappa_shear_is_used, kappa_shear_at_vertex

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure holds the parameters that regulate shear mixing
type, public :: Kappa_shear_CS ; private
  real    :: RiNo_crit       !< The critical shear Richardson number for
                             !! shear-entrainment [nondim]. The theoretical value is 0.25.
                             !! The values found by Jackson et al. are 0.25-0.35.
  real    :: Shearmix_rate   !< A nondimensional rate scale for shear-driven
                             !! entrainment [nondim].  The value given by Jackson et al.
                             !! is 0.085-0.089.
  real    :: FRi_curvature   !<   A constant giving the curvature of the function
                             !! of the Richardson number that relates shear to
                             !! sources in the kappa equation [nondim].
                             !! The values found by Jackson et al. are -0.97 - -0.89.
  real    :: C_N             !<   The coefficient for the decay of TKE due to
                             !! stratification (i.e. proportional to N*tke) [nondim].
                             !! The values found by Jackson et al. are 0.24-0.28.
  real    :: C_S             !<   The coefficient for the decay of TKE due to
                             !! shear (i.e. proportional to |S|*tke) [nondim].
                             !! The values found by Jackson et al. are 0.14-0.12.
  real    :: lambda          !<   The coefficient for the buoyancy length scale
                             !! in the kappa equation [nondim].
                             !! The values found by Jackson et al. are 0.82-0.81.
  real    :: lambda2_N_S     !<   The square of the ratio of the coefficients of
                             !! the buoyancy and shear scales in the diffusivity
                             !! equation, 0 to eliminate the shear scale [nondim].
  real    :: TKE_bg          !<   The background level of TKE [Z2 T-2 ~> m2 s-2].
  real    :: kappa_0         !<   The background diapycnal diffusivity [Z2 T-1 ~> m2 s-1].
  real    :: kappa_trunc     !< Diffusivities smaller than this are rounded to 0 [Z2 T-1 ~> m2 s-1].
  real    :: kappa_tol_err   !<   The fractional error in kappa that is tolerated [nondim].
  real    :: Prandtl_turb    !< Prandtl number used to convert Kd_shear into viscosity [nondim].
  integer :: nkml            !<   The number of layers in the mixed layer, as
                             !! treated in this routine.  If the pieces of the
                             !! mixed layer are not to be treated collectively,
                             !! nkml is set to 1.
  integer :: max_RiNo_it     !< The maximum number of iterations that may be used
                             !! to estimate the instantaneous shear-driven mixing.
  integer :: max_KS_it       !< The maximum number of iterations that may be used
                             !! to estimate the time-averaged diffusivity.
  logical :: dKdQ_iteration_bug !< If true. use an older, dimensionally inconsistent estimate of
                             !! the derivative of diffusivity with energy in the Newton's method
                             !! iteration.  The bug causes undercorrections when dz > 1m.
  logical :: KS_at_vertex    !< If true, do the calculations of the shear-driven mixing
                             !! at the cell vertices (i.e., the vorticity points).
  logical :: eliminate_massless !< If true, massless layers are merged with neighboring
                             !! massive layers in this calculation.
                             !  I can think of no good reason why this should be false. - RWH
  real    :: vel_underflow   !< Velocity components smaller than vel_underflow
                             !! are set to 0 [L T-1 ~> m s-1].
  real    :: kappa_src_max_chg !< The maximum permitted increase in the kappa source within an
                             !! iteration relative to the local source [nondim].  This must be
                             !! greater than 1.  The lower limit for the permitted fractional
                             !! decrease is (1 - 0.5/kappa_src_max_chg).  These limits could
                             !! perhaps be made dynamic with an improved iterative solver.
  logical :: all_layer_TKE_bug !< If true, report back the latest estimate of TKE instead of the
                             !! time average TKE when there is mass in all layers.  Otherwise always
                             !! report the time-averaged TKE, as is currently done when there
                             !! are some massless layers.
!  logical :: layer_stagger = .false. ! If true, do the calculations centered at
                             !  layers, rather than the interfaces.
  logical :: debug = .false. !< If true, write verbose debugging messages.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
  !>@{ Diagnostic IDs
  integer :: id_Kd_shear = -1, id_TKE = -1, id_ILd2 = -1, id_dz_Int = -1
  !>@}
end type Kappa_shear_CS

! integer :: id_clock_project, id_clock_KQ, id_clock_avg, id_clock_setup

contains

!> Subroutine for calculating shear-driven diffusivity and TKE in tracer columns
subroutine Calculate_kappa_shear(u_in, v_in, h, tv, p_surf, kappa_io, tke_io, &
                                 kv_io, dt, G, GV, US, CS, initialize_all)
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: u_in   !< Initial zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: v_in   !< Initial meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent fields
                                                   !! have NULL ptrs.
  real, dimension(:,:),    pointer       :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa] (or NULL).
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(inout) :: kappa_io !< The diapycnal diffusivity at each interface
                                                   !! (not layer!) [Z2 T-1 ~> m2 s-1].  Initially this is the
                                                   !! value from the previous timestep, which may
                                                   !! accelerate the iteration toward convergence.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out) :: tke_io   !< The turbulent kinetic energy per unit mass at
                                                   !! each interface (not layer!) [Z2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(inout) :: kv_io  !< The vertical viscosity at each interface
                                                   !! (not layer!) [Z2 T-1 ~> m2 s-1]. This discards any
                                                   !! previous value (i.e. it is intent out) and
                                                   !! simply sets Kv = Prandtl * Kd_shear
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  type(Kappa_shear_CS),    pointer       :: CS     !< The control structure returned by a previous
                                                   !! call to kappa_shear_init.
  logical,       optional, intent(in)    :: initialize_all !< If present and false, the previous
                                                   !! value of kappa is used to start the iterations

  ! Local variables
  real, dimension(SZI_(G),SZK_(GV)) :: &
    h_2d, &             ! A 2-D version of h, but converted to [Z ~> m].
    u_2d, v_2d, &       ! 2-D versions of u_in and v_in, converted to [L T-1 ~> m s-1].
    T_2d, S_2d, rho_2d  ! 2-D versions of T [degC], S [ppt], and rho [R ~> kg m-3].
  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    kappa_2d, & ! 2-D version of kappa_io [Z2 T-1 ~> m2 s-1].
    tke_2d      ! 2-D version tke_io [Z2 T-2 ~> m2 s-2].
  real, dimension(SZK_(GV)) :: &
    Idz, &      ! The inverse of the distance between TKE points [Z-1 ~> m-1].
    dz, &       ! The layer thickness [Z ~> m].
    u0xdz, &    ! The initial zonal velocity times dz [Z L T-1 ~> m2 s-1].
    v0xdz, &    ! The initial meridional velocity times dz [Z L T-1 ~> m2 s-1].
    T0xdz, &    ! The initial temperature times dz [degC Z ~> degC m].
    S0xdz       ! The initial salinity times dz [ppt Z ~> ppt m].
  real, dimension(SZK_(GV)+1) :: &
    kappa, &    ! The shear-driven diapycnal diffusivity at an interface [Z2 T-1 ~> m2 s-1].
    tke, &      ! The Turbulent Kinetic Energy per unit mass at an interface [Z2 T-2 ~> m2 s-2].
    kappa_avg, & ! The time-weighted average of kappa [Z2 T-1 ~> m2 s-1].
    tke_avg     ! The time-weighted average of TKE [Z2 T-2 ~> m2 s-2].
  real :: f2   ! The squared Coriolis parameter of each column [T-2 ~> s-2].
  real :: surface_pres  ! The top surface pressure [R L2 T-2 ~> Pa].

  real :: dz_in_lay     !   The running sum of the thickness in a layer [Z ~> m].
  real :: k0dt          ! The background diffusivity times the timestep [Z2 ~> m2].
  real :: dz_massless   ! A layer thickness that is considered massless [Z ~> m].
  logical :: use_temperature  !  If true, temperature and salinity have been
                        ! allocated and are being used as state variables.
  logical :: new_kappa = .true. ! If true, ignore the value of kappa from the
                        ! last call to this subroutine.

  integer, dimension(SZK_(GV)+1) :: kc ! The index map between the original
                        ! interfaces and the interfaces with massless layers
                        ! merged into nearby massive layers.
  real, dimension(SZK_(GV)+1) :: kf ! The fractional weight of interface kc+1 for
                        ! interpolating back to the original index space [nondim].
  integer :: is, ie, js, je, i, j, k, nz, nzc

  is = G%isc ; ie = G%iec; js = G%jsc ; je = G%jec ; nz = GV%ke

  use_temperature = .false. ; if (associated(tv%T)) use_temperature = .true.
  new_kappa = .true. ; if (present(initialize_all)) new_kappa = initialize_all

  k0dt = dt*CS%kappa_0
  dz_massless = 0.1*sqrt(k0dt)

  !$OMP parallel do default(private) shared(js,je,is,ie,nz,h,u_in,v_in,use_temperature,new_kappa, &
  !$OMP                                tv,G,GV,US,CS,kappa_io,dz_massless,k0dt,p_surf,dt,tke_io,kv_io)
  do j=js,je
    do k=1,nz ; do i=is,ie
      h_2d(i,k) = h(i,j,k)*GV%H_to_Z
      u_2d(i,k) = u_in(i,j,k) ; v_2d(i,k) = v_in(i,j,k)
    enddo ; enddo
    if (use_temperature) then ; do k=1,nz ; do i=is,ie
      T_2d(i,k) = tv%T(i,j,k) ; S_2d(i,k) = tv%S(i,j,k)
    enddo ; enddo ; else ; do k=1,nz ; do i=is,ie
      rho_2d(i,k) = GV%Rlay(k) ! Could be tv%Rho(i,j,k) ?
    enddo ; enddo ; endif
    if (.not.new_kappa) then ; do K=1,nz+1 ; do i=is,ie
      kappa_2d(i,K) = kappa_io(i,j,K)
    enddo ; enddo ; endif

!---------------------------------------
! Work on each column.
!---------------------------------------
    do i=is,ie ; if (G%mask2dT(i,j) > 0.5) then
    ! call cpu_clock_begin(id_clock_setup)
      ! Store a transposed version of the initial arrays.
      ! Any elimination of massless layers would occur here.
      if (CS%eliminate_massless) then
        nzc = 1
        do k=1,nz
          ! Zero out the thicknesses of all layers, even if they are unused.
          dz(k) = 0.0 ; u0xdz(k) = 0.0 ; v0xdz(k) = 0.0
          T0xdz(k) = 0.0 ; S0xdz(k) = 0.0

          ! Add a new layer if this one has mass.
!          if ((dz(nzc) > 0.0) .and. (h_2d(i,k) > dz_massless)) nzc = nzc+1
          if ((k>CS%nkml) .and. (dz(nzc) > 0.0) .and. &
              (h_2d(i,k) > dz_massless)) nzc = nzc+1

          ! Only merge clusters of massless layers.
!         if ((dz(nzc) > dz_massless) .or. &
!             ((dz(nzc) > 0.0) .and. (h_2d(i,k) > dz_massless))) nzc = nzc+1

          kc(k) = nzc
          dz(nzc) = dz(nzc) + h_2d(i,k)
          u0xdz(nzc) = u0xdz(nzc) + u_2d(i,k)*h_2d(i,k)
          v0xdz(nzc) = v0xdz(nzc) + v_2d(i,k)*h_2d(i,k)
          if (use_temperature) then
            T0xdz(nzc) = T0xdz(nzc) + T_2d(i,k)*h_2d(i,k)
            S0xdz(nzc) = S0xdz(nzc) + S_2d(i,k)*h_2d(i,k)
          else
            T0xdz(nzc) = T0xdz(nzc) + rho_2d(i,k)*h_2d(i,k)
            S0xdz(nzc) = S0xdz(nzc) + rho_2d(i,k)*h_2d(i,k)
          endif
        enddo
        kc(nz+1) = nzc+1

        ! Set up Idz as the inverse of layer thicknesses.
        do k=1,nzc ; Idz(k) = 1.0 / dz(k) ; enddo

        !   Now determine kf, the fractional weight of interface kc when
        ! interpolating between interfaces kc and kc+1.
        kf(1) = 0.0 ; dz_in_lay = h_2d(i,1)
        do k=2,nz
          if (kc(k) > kc(k-1)) then
            kf(k) = 0.0 ; dz_in_lay = h_2d(i,k)
          else
            kf(k) = dz_in_lay*Idz(kc(k)) ; dz_in_lay = dz_in_lay + h_2d(i,k)
          endif
        enddo
        kf(nz+1) = 0.0
      else
        do k=1,nz
          dz(k) = h_2d(i,k)
          u0xdz(k) = u_2d(i,k)*dz(k) ; v0xdz(k) = v_2d(i,k)*dz(k)
        enddo
        if (use_temperature) then
          do k=1,nz
            T0xdz(k) = T_2d(i,k)*dz(k) ; S0xdz(k) = S_2d(i,k)*dz(k)
          enddo
        else
          do k=1,nz
            T0xdz(k) = rho_2d(i,k)*dz(k) ; S0xdz(k) = rho_2d(i,k)*dz(k)
          enddo
        endif
        nzc = nz
        do k=1,nzc+1 ; kc(k) = k ; kf(k) = 0.0 ; enddo
      endif
      f2 = 0.25 * ((G%CoriolisBu(I,j)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                   (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J)**2))
      surface_pres = 0.0 ; if (associated(p_surf)) surface_pres = p_surf(i,j)

    ! ----------------------------------------------------    I_Ld2_1d, dz_Int_1d

    ! Set the initial guess for kappa, here defined at interfaces.
    ! ----------------------------------------------------
      if (new_kappa) then
        do K=1,nzc+1 ; kappa(K) = US%m2_s_to_Z2_T*1.0 ; enddo
      else
        do K=1,nzc+1 ; kappa(K) = kappa_2d(i,K) ; enddo
      endif

      call kappa_shear_column(kappa, tke, dt, nzc, f2, surface_pres, &
                              dz, u0xdz, v0xdz, T0xdz, S0xdz, kappa_avg, &
                              tke_avg, tv, CS, GV, US)

    ! call cpu_clock_begin(id_clock_setup)
    ! Extrapolate from the vertically reduced grid back to the original layers.
      if (nz == nzc) then
        do K=1,nz+1
          kappa_2d(i,K) = kappa_avg(K)
          if (CS%all_layer_TKE_bug) then
            tke_2d(i,K) = tke(K)
          else
            tke_2d(i,K) = tke_avg(K)
          endif
        enddo
      else
        do K=1,nz+1
          if (kf(K) == 0.0) then
            kappa_2d(i,K) = kappa_avg(kc(K))
            tke_2d(i,K) = tke_avg(kc(K))
          else
            kappa_2d(i,K) = (1.0-kf(K)) * kappa_avg(kc(K)) + &
                             kf(K) * kappa_avg(kc(K)+1)
            tke_2d(i,K) = (1.0-kf(K)) * tke_avg(kc(K)) + &
                           kf(K) * tke_avg(kc(K)+1)
          endif
        enddo
      endif
    ! call cpu_clock_end(id_clock_setup)
    else  ! Land points, still inside the i-loop.
      do K=1,nz+1
        kappa_2d(i,K) = 0.0 ; tke_2d(i,K) = 0.0
      enddo
    endif ; enddo ! i-loop

    do K=1,nz+1 ; do i=is,ie
      kappa_io(i,j,K) = G%mask2dT(i,j) * kappa_2d(i,K)
      tke_io(i,j,K) = G%mask2dT(i,j) * tke_2d(i,K)
      kv_io(i,j,K) = ( G%mask2dT(i,j) * kappa_2d(i,K) ) * CS%Prandtl_turb
    enddo ; enddo

  enddo ! end of j-loop

  if (CS%debug) then
    call hchksum(kappa_io, "kappa", G%HI, scale=US%Z2_T_to_m2_s)
    call hchksum(tke_io, "tke", G%HI, scale=US%Z_to_m**2*US%s_to_T**2)
  endif

  if (CS%id_Kd_shear > 0) call post_data(CS%id_Kd_shear, kappa_io, CS%diag)
  if (CS%id_TKE > 0) call post_data(CS%id_TKE, tke_io, CS%diag)

end subroutine Calculate_kappa_shear


!> Subroutine for calculating shear-driven diffusivity and TKE in corner columns
subroutine Calc_kappa_shear_vertex(u_in, v_in, h, T_in, S_in, tv, p_surf, kappa_io, tke_io, &
                                   kv_io, dt, G, GV, US, CS, initialize_all)
  type(ocean_grid_type),   intent(in)    :: G      !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),    intent(in)   :: US     !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: u_in   !< Initial zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   &
                           intent(in)    :: v_in   !< Initial meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: h      !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: T_in   !< Layer potential temperatures [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   &
                           intent(in)    :: S_in   !< Layer salinities in ppt.
  type(thermo_var_ptrs),   intent(in)    :: tv     !< A structure containing pointers to any
                                                   !! available thermodynamic fields. Absent fields
                                                   !! have NULL ptrs.
  real, dimension(:,:),    pointer       :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa]
                                                   !! (or NULL).
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(out)   :: kappa_io !< The diapycnal diffusivity at each interface
                                                   !! (not layer!) [Z2 T-1 ~> m2 s-1].
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), &
                           intent(out)   :: tke_io !< The turbulent kinetic energy per unit mass at
                                                   !! each interface (not layer!) [Z2 T-2 ~> m2 s-2].
  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)+1), &
                           intent(inout) :: kv_io  !< The vertical viscosity at each interface [Z2 T-1 ~> m2 s-1].
                                                   !! The previous value is used to initialize kappa
                                                   !! in the vertex columes as Kappa = Kv/Prandtl
                                                   !! to accelerate the iteration toward covergence.
  real,                    intent(in)    :: dt     !< Time increment [T ~> s].
  type(Kappa_shear_CS),    pointer       :: CS     !< The control structure returned by a previous
                                                   !! call to kappa_shear_init.
  logical,       optional, intent(in)    :: initialize_all !< If present and false, the previous
                                                   !! value of kappa is used to start the iterations

  ! Local variables
  real, dimension(SZIB_(G),SZK_(GV)) :: &
    h_2d, &             ! A 2-D version of h, but converted to [Z ~> m].
    u_2d, v_2d, &       ! 2-D versions of u_in and v_in, converted to [L T-1 ~> m s-1].
    T_2d, S_2d, rho_2d  ! 2-D versions of T [degC], S [ppt], and rho [R ~> kg m-3].
  real, dimension(SZIB_(G),SZK_(GV)+1,2) :: &
    kappa_2d    ! Quasi 2-D versions of kappa_io [Z2 T-1 ~> m2 s-1].
  real, dimension(SZIB_(G),SZK_(GV)+1) :: &
    tke_2d      ! 2-D version tke_io [Z2 T-2 ~> m2 s-2].
  real, dimension(SZK_(GV)) :: &
    Idz, &      ! The inverse of the distance between TKE points [Z-1 ~> m-1].
    dz, &       ! The layer thickness [Z ~> m].
    u0xdz, &    ! The initial zonal velocity times dz [L Z T-1 ~> m2 s-1].
    v0xdz, &    ! The initial meridional velocity times dz [L Z T-1 ~> m2 s-1].
    T0xdz, &    ! The initial temperature times dz [degC Z ~> degC m].
    S0xdz       ! The initial salinity times dz [ppt Z ~> ppt m].
  real, dimension(SZK_(GV)+1) :: &
    kappa, &    ! The shear-driven diapycnal diffusivity at an interface [Z2 T-1 ~> m2 s-1].
    tke, &      ! The Turbulent Kinetic Energy per unit mass at an interface [Z2 T-2 ~> m2 s-2].
    kappa_avg, & ! The time-weighted average of kappa [Z2 T-1 ~> m2 s-1].
    tke_avg     ! The time-weighted average of TKE [Z2 T-2 ~> m2 s-2].
  real :: f2   ! The squared Coriolis parameter of each column [T-2 ~> s-2].
  real :: surface_pres  ! The top surface pressure [R L2 T-2 ~> Pa].

  real :: dz_in_lay     !   The running sum of the thickness in a layer [Z ~> m].
  real :: k0dt          ! The background diffusivity times the timestep [Z2 ~> m2].
  real :: dz_massless   ! A layer thickness that is considered massless [Z ~> m].
  real :: I_hwt         ! The inverse of the masked thickness weights [H-1 ~> m-1 or m2 kg-1].
  real :: I_Prandtl     ! The inverse of the turbulent Prandtl number [nondim].
  logical :: use_temperature  !  If true, temperature and salinity have been
                        ! allocated and are being used as state variables.
  logical :: new_kappa = .true. ! If true, ignore the value of kappa from the
                        ! last call to this subroutine.
  logical :: do_i       ! If true, work on this column.

  integer, dimension(SZK_(GV)+1) :: kc ! The index map between the original
                        ! interfaces and the interfaces with massless layers
                        ! merged into nearby massive layers.
  real, dimension(SZK_(GV)+1) :: kf ! The fractional weight of interface kc+1 for
                        ! interpolating back to the original index space [nondim].
  integer :: IsB, IeB, JsB, JeB, i, j, k, nz, nzc, J2, J2m1

  ! Diagnostics that should be deleted?
  isB = G%isc-1 ; ieB = G%iecB ; jsB = G%jsc-1 ; jeB = G%jecB ; nz = GV%ke

  use_temperature = .false. ; if (associated(tv%T)) use_temperature = .true.
  new_kappa = .true. ; if (present(initialize_all)) new_kappa = initialize_all

  k0dt =  dt*CS%kappa_0
  dz_massless = 0.1*sqrt(k0dt)
  I_Prandtl = 0.0 ; if (CS%Prandtl_turb > 0.0) I_Prandtl = 1.0 / CS%Prandtl_turb

  !$OMP parallel do default(private) shared(jsB,jeB,isB,ieB,nz,h,u_in,v_in,use_temperature,new_kappa, &
  !$OMP                                tv,G,GV,US,CS,kappa_io,dz_massless,k0dt,p_surf,dt,tke_io,kv_io,I_Prandtl)
  do J=JsB,JeB
    J2 = mod(J,2)+1 ; J2m1 = 3-J2 ! = mod(J-1,2)+1

    ! Interpolate the various quantities to the corners, using masks.
    do k=1,nz ; do I=IsB,IeB
      u_2d(I,k) = (u_in(I,j,k)   * (G%mask2dCu(I,j)   * (h(i,j,k)   + h(i+1,j,k))) + &
                   u_in(I,j+1,k) * (G%mask2dCu(I,j+1) * (h(i,j+1,k) + h(i+1,j+1,k))) ) / &
                  ((G%mask2dCu(I,j)   * (h(i,j,k)   + h(i+1,j,k)) + &
                    G%mask2dCu(I,j+1) * (h(i,j+1,k) + h(i+1,j+1,k))) + GV%H_subroundoff)
      v_2d(I,k) = (v_in(i,J,k)   * (G%mask2dCv(i,J)   * (h(i,j,k)   + h(i,j+1,k))) + &
                   v_in(i+1,J,k) * (G%mask2dCv(i+1,J) * (h(i+1,j,k) + h(i+1,j+1,k))) ) / &
                  ((G%mask2dCv(i,J)   * (h(i,j,k)   + h(i,j+1,k)) + &
                    G%mask2dCv(i+1,J) * (h(i+1,j,k) + h(i+1,j+1,k))) + GV%H_subroundoff)
      I_hwt = 1.0 / (((G%mask2dT(i,j) * h(i,j,k) + G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) + &
                      (G%mask2dT(i+1,j) * h(i+1,j,k) + G%mask2dT(i,j+1) * h(i,j+1,k))) + &
                     GV%H_subroundoff)
      if (use_temperature) then
        T_2d(I,k) = ( ((G%mask2dT(i,j) * h(i,j,k)) * T_in(i,j,k) + &
                       (G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) * T_in(i+1,j+1,k)) + &
                      ((G%mask2dT(i+1,j) * h(i+1,j,k)) * T_in(i+1,j,k) + &
                       (G%mask2dT(i,j+1) * h(i,j+1,k)) * T_in(i,j+1,k)) ) * I_hwt
        S_2d(I,k) = ( ((G%mask2dT(i,j) * h(i,j,k)) * S_in(i,j,k) + &
                       (G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) * S_in(i+1,j+1,k)) + &
                      ((G%mask2dT(i+1,j) * h(i+1,j,k)) * S_in(i+1,j,k) + &
                       (G%mask2dT(i,j+1) * h(i,j+1,k)) * S_in(i,j+1,k)) ) * I_hwt
      endif
      h_2d(I,k) = GV%H_to_Z * ((G%mask2dT(i,j) * h(i,j,k) + G%mask2dT(i+1,j+1) * h(i+1,j+1,k)) + &
                               (G%mask2dT(i+1,j) * h(i+1,j,k) + G%mask2dT(i,j+1) * h(i,j+1,k)) ) / &
                              ((G%mask2dT(i,j) + G%mask2dT(i+1,j+1)) + &
                               (G%mask2dT(i+1,j) + G%mask2dT(i,j+1)) + 1.0e-36 )
!      h_2d(I,k) = 0.25*((h(i,j,k) + h(i+1,j+1,k)) + (h(i+1,j,k) + h(i,j+1,k)))*GV%H_to_Z
!      h_2d(I,k) = ((h(i,j,k)**2 + h(i+1,j+1,k)**2) + &
!                   (h(i+1,j,k)**2 + h(i,j+1,k)**2))*GV%H_to_Z * I_hwt
    enddo ; enddo
    if (.not.use_temperature) then ; do k=1,nz ; do I=IsB,IeB
      rho_2d(I,k) = GV%Rlay(k)
    enddo ; enddo ; endif
    if (.not.new_kappa) then ; do K=1,nz+1 ; do I=IsB,IeB
      kappa_2d(I,K,J2) = kv_io(I,J,K) * I_Prandtl
    enddo ; enddo ; endif

!---------------------------------------
! Work on each column.
!---------------------------------------
    do I=IsB,IeB ; if ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) + &
                       (G%mask2dCv(i,J) + G%mask2dCv(i+1,J)) > 0.0) then
    ! call cpu_clock_begin(Id_clock_setup)
      ! Store a transposed version of the initial arrays.
      ! Any elimination of massless layers would occur here.
      if (CS%eliminate_massless) then
        nzc = 1
        do k=1,nz
          ! Zero out the thicknesses of all layers, even if they are unused.
          dz(k) = 0.0 ; u0xdz(k) = 0.0 ; v0xdz(k) = 0.0
          T0xdz(k) = 0.0 ; S0xdz(k) = 0.0

          ! Add a new layer if this one has mass.
!          if ((dz(nzc) > 0.0) .and. (h_2d(I,k) > dz_massless)) nzc = nzc+1
          if ((k>CS%nkml) .and. (dz(nzc) > 0.0) .and. &
              (h_2d(I,k) > dz_massless)) nzc = nzc+1

          ! Only merge clusters of massless layers.
!         if ((dz(nzc) > dz_massless) .or. &
!             ((dz(nzc) > 0.0) .and. (h_2d(I,k) > dz_massless))) nzc = nzc+1

          kc(k) = nzc
          dz(nzc) = dz(nzc) + h_2d(I,k)
          u0xdz(nzc) = u0xdz(nzc) + u_2d(I,k)*h_2d(I,k)
          v0xdz(nzc) = v0xdz(nzc) + v_2d(I,k)*h_2d(I,k)
          if (use_temperature) then
            T0xdz(nzc) = T0xdz(nzc) + T_2d(I,k)*h_2d(I,k)
            S0xdz(nzc) = S0xdz(nzc) + S_2d(I,k)*h_2d(I,k)
          else
            T0xdz(nzc) = T0xdz(nzc) + rho_2d(I,k)*h_2d(I,k)
            S0xdz(nzc) = S0xdz(nzc) + rho_2d(I,k)*h_2d(I,k)
          endif
        enddo
        kc(nz+1) = nzc+1

        ! Set up Idz as the inverse of layer thicknesses.
        do k=1,nzc ; Idz(k) = 1.0 / dz(k) ; enddo

        !   Now determine kf, the fractional weight of interface kc when
        ! interpolating between interfaces kc and kc+1.
        kf(1) = 0.0 ; dz_in_lay = h_2d(I,1)
        do k=2,nz
          if (kc(k) > kc(k-1)) then
            kf(k) = 0.0 ; dz_in_lay = h_2d(I,k)
          else
            kf(k) = dz_in_lay*Idz(kc(k)) ; dz_in_lay = dz_in_lay + h_2d(I,k)
          endif
        enddo
        kf(nz+1) = 0.0
      else
        do k=1,nz
          dz(k) = h_2d(I,k)
          u0xdz(k) = u_2d(I,k)*dz(k) ; v0xdz(k) = v_2d(I,k)*dz(k)
        enddo
        if (use_temperature) then
          do k=1,nz
            T0xdz(k) = T_2d(I,k)*dz(k) ; S0xdz(k) = S_2d(I,k)*dz(k)
          enddo
        else
          do k=1,nz
            T0xdz(k) = rho_2d(I,k)*dz(k) ; S0xdz(k) = rho_2d(I,k)*dz(k)
          enddo
        endif
        nzc = nz
        do k=1,nzc+1 ; kc(k) = k ; kf(k) = 0.0 ; enddo
      endif
      f2 = G%CoriolisBu(I,J)**2
      surface_pres = 0.0 ; if (associated(p_surf)) &
        surface_pres = 0.25 * ((p_surf(i,j) + p_surf(i+1,j+1)) + &
                               (p_surf(i+1,j) + p_surf(i,j+1)))

    ! ----------------------------------------------------
    ! Set the initial guess for kappa, here defined at interfaces.
    ! ----------------------------------------------------
      if (new_kappa) then
        do K=1,nzc+1 ; kappa(K) = US%m2_s_to_Z2_T*1.0 ; enddo
      else
        do K=1,nzc+1 ; kappa(K) = kappa_2d(I,K,J2) ; enddo
      endif

      call kappa_shear_column(kappa, tke, dt, nzc, f2, surface_pres, &
                              dz, u0xdz, v0xdz, T0xdz, S0xdz, kappa_avg, &
                              tke_avg, tv, CS, GV, US)
    ! call cpu_clock_begin(Id_clock_setup)
    ! Extrapolate from the vertically reduced grid back to the original layers.
      if (nz == nzc) then
        do K=1,nz+1
          kappa_2d(I,K,J2) = kappa_avg(K)
          if (CS%all_layer_TKE_bug) then
            tke_2d(i,K) = tke(K)
          else
            tke_2d(i,K) = tke_avg(K)
          endif
        enddo
      else
        do K=1,nz+1
          if (kf(K) == 0.0) then
            kappa_2d(I,K,J2) = kappa_avg(kc(K))
            tke_2d(I,K) = tke_avg(kc(K))
          else
            kappa_2d(I,K,J2) = (1.0-kf(K)) * kappa_avg(kc(K)) + kf(K) * kappa_avg(kc(K)+1)
            tke_2d(I,K) = (1.0-kf(K)) * tke_avg(kc(K)) + kf(K) * tke_avg(kc(K)+1)
          endif
        enddo
      endif
    ! call cpu_clock_end(Id_clock_setup)
    else  ! Land points, still inside the i-loop.
      do K=1,nz+1
        kappa_2d(I,K,J2) = 0.0 ; tke_2d(I,K) = 0.0
      enddo
    endif ; enddo ! i-loop

    do K=1,nz+1 ; do I=IsB,IeB
      tke_io(I,J,K) = G%mask2dBu(I,J) * tke_2d(I,K)
      kv_io(I,J,K) = ( G%mask2dBu(I,J) * kappa_2d(I,K,J2) ) * CS%Prandtl_turb
    enddo ; enddo
    if (J>=G%jsc) then ; do K=1,nz+1 ; do i=G%isc,G%iec
      ! Set the diffusivities in tracer columns from the values at vertices.
      kappa_io(i,j,K) = G%mask2dT(i,j) * 0.25 * &
                        ((kappa_2d(I-1,K,J2m1) + kappa_2d(I,K,J2)) + &
                         (kappa_2d(I-1,K,J2)   + kappa_2d(I,K,J2m1)))
    enddo ; enddo ; endif

  enddo ! end of J-loop

  if (CS%debug) then
    call hchksum(kappa_io, "kappa", G%HI, scale=US%Z2_T_to_m2_s)
    call Bchksum(tke_io, "tke", G%HI, scale=US%Z_to_m**2*US%s_to_T**2)
  endif

  if (CS%id_Kd_shear > 0) call post_data(CS%id_Kd_shear, kappa_io, CS%diag)
  if (CS%id_TKE > 0) call post_data(CS%id_TKE, tke_io, CS%diag)

end subroutine Calc_kappa_shear_vertex


!> This subroutine calculates shear-driven diffusivity and TKE in a single column
subroutine kappa_shear_column(kappa, tke, dt, nzc, f2, surface_pres, &
                              dz, u0xdz, v0xdz, T0xdz, S0xdz, kappa_avg, &
                              tke_avg, tv, CS, GV, US, I_Ld2_1d, dz_Int_1d)
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)+1), &
                     intent(inout) :: kappa !< The time-weighted average of kappa [Z2 T-1 ~> m2 s-1].
  real, dimension(SZK_(GV)+1), &
                     intent(out)   :: tke  !< The Turbulent Kinetic Energy per unit mass at
                                           !! an interface [Z2 T-2 ~> m2 s-2].
  integer,           intent(in)    :: nzc  !< The number of active layers in the column.
  real,              intent(in)    :: f2   !< The square of the Coriolis parameter [T-2 ~> s-2].
  real,              intent(in)    :: surface_pres  !< The surface pressure [R L2 T-2 ~> Pa].
  real, dimension(SZK_(GV)), &
                     intent(in)    :: dz   !< The layer thickness [Z ~> m].
  real, dimension(SZK_(GV)), &
                     intent(in)    :: u0xdz !< The initial zonal velocity times dz [Z L T-1 ~> m2 s-1].
  real, dimension(SZK_(GV)), &
                     intent(in)    :: v0xdz !< The initial meridional velocity times dz [Z L T-1 ~> m2 s-1].
  real, dimension(SZK_(GV)), &
                     intent(in)    :: T0xdz !< The initial temperature times dz [degC Z ~> degC m].
  real, dimension(SZK_(GV)), &
                     intent(in)    :: S0xdz !< The initial salinity times dz [ppt Z ~> ppt m].
  real, dimension(SZK_(GV)+1), &
                     intent(out)   :: kappa_avg !< The time-weighted average of kappa [Z2 T-1 ~> m2 s-1].
  real, dimension(SZK_(GV)+1), &
                     intent(out)   :: tke_avg  !< The time-weighted average of TKE [Z2 T-2 ~> m2 s-2].
  real,                    intent(in)    :: dt !< Time increment [T ~> s].
  type(thermo_var_ptrs),   intent(in)    :: tv !< A structure containing pointers to any
                                               !! available thermodynamic fields. Absent fields
                                               !! have NULL ptrs.
  type(Kappa_shear_CS),    pointer       :: CS !< The control structure returned by a previous
                                               !! call to kappa_shear_init.
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  real,  dimension(SZK_(GV)+1), &
           optional, intent(out)   :: I_Ld2_1d !< The inverse of the squared mixing length [Z-2 ~> m-2].
  real,  dimension(SZK_(GV)+1), &
           optional, intent(out)   :: dz_Int_1d !< The extent of a finite-volume space surrounding an interface,
                                               !! as used in calculating kappa and TKE [Z ~> m].

  real, dimension(nzc) :: &
    u, &        ! The zonal velocity after a timestep of mixing [L T-1 ~> m s-1].
    v, &        ! The meridional velocity after a timestep of mixing [L T-1 ~> m s-1].
    Idz, &      ! The inverse of the distance between TKE points [Z-1 ~> m-1].
    T, &        ! The potential temperature after a timestep of mixing [degC].
    Sal, &      ! The salinity after a timestep of mixing [ppt].
    u_test, v_test, & ! Temporary velocities [L T-1 ~> m s-1].
    T_test, S_test ! Temporary temperatures [degC] and salinities [ppt].

  real, dimension(nzc+1) :: &
    N2, &       ! The squared buoyancy frequency at an interface [T-2 ~> s-2].
    dz_Int, &   ! The extent of a finite-volume space surrounding an interface,
                ! as used in calculating kappa and TKE [Z ~> m].
    I_dz_int, & ! The inverse of the distance between velocity & density points
                ! above and below an interface [Z-1 ~> m-1].  This is used to
                ! calculate N2, shear, and fluxes, and it might differ from
                ! 1/dz_Int, as they have different uses.
    S2, &       ! The squared shear at an interface [T-2 ~> s-2].
    a1, &       ! a1 is the coupling between adjacent interfaces in the TKE,
                ! velocity, and density equations [Z s-1 ~> m s-1] or [Z ~> m]
    c1, &       ! c1 is used in the tridiagonal (and similar) solvers.
    k_src, &    ! The shear-dependent source term in the kappa equation [T-1 ~> s-1].
    kappa_src, & ! The shear-dependent source term in the kappa equation [T-1 ~> s-1].
    kappa_out, & ! The kappa that results from the kappa equation [Z2 T-1 ~> m2 s-1].
    kappa_mid, & ! The average of the initial and predictor estimates of kappa [Z2 T-1 ~> m2 s-1].
    tke_pred, & ! The value of TKE from a predictor step [Z2 T-2 ~> m2 s-2].
    kappa_pred, & ! The value of kappa from a predictor step [Z2 T-1 ~> m2 s-1].
    pressure, & ! The pressure at an interface [R L2 T-2 ~> Pa].
    T_int, &    ! The temperature interpolated to an interface [degC].
    Sal_int, &  ! The salinity interpolated to an interface [ppt].
    dbuoy_dT, & ! The partial derivatives of buoyancy with changes in temperature
    dbuoy_dS, & ! and salinity, [Z T-2 degC-1 ~> m s-2 degC-1] and [Z T-2 ppt-1 ~> m s-2 ppt-1].
    I_L2_bdry, &   ! The inverse of the square of twice the harmonic mean
                   ! distance to the top and bottom boundaries [Z-2 ~> m-2].
    K_Q, &         ! Diffusivity divided by TKE [T ~> s].
    K_Q_tmp, &     ! A temporary copy of diffusivity divided by TKE [T ~> s].
    local_src_avg, & ! The time-integral of the local source [nondim].
    tol_min, & ! Minimum tolerated ksrc for the corrector step [T-1 ~> s-1].
    tol_max, & ! Maximum tolerated ksrc for the corrector step [T-1 ~> s-1].
    tol_chg, & ! The tolerated kappa change integrated over a timestep [nondim].
    dist_from_top, &  ! The distance from the top surface [Z ~> m].
    local_src     ! The sum of all sources of kappa, including kappa_src and
                  ! sources from the elliptic term [T-1 ~> s-1].

  real :: dist_from_bot ! The distance from the bottom surface [Z ~> m].
  real :: b1            ! The inverse of the pivot in the tridiagonal equations.
  real :: bd1           ! A term in the denominator of b1.
  real :: d1            ! 1 - c1 in the tridiagonal equations.
  real :: gR0           ! A conversion factor from Z to pressure, given by Rho_0 times g
                        ! [R L2 T-2 Z-1 ~> kg m-2 s-2].
  real :: g_R0          ! g_R0 is a rescaled version of g/Rho [Z R-1 T-2 ~> m4 kg-1 s-2].
  real :: Norm          ! A factor that normalizes two weights to 1 [Z-2 ~> m-2].
  real :: tol_dksrc     ! Tolerance for the change in the kappa source within an iteration
                        ! relative to the local source [nondim].  This must be greater than 1.
  real :: tol2          ! The tolerance for the change in the kappa source within an iteration
                        ! relative to the average local source over previous iterations [nondim].
  real :: tol_dksrc_low ! The tolerance for the fractional decrease in ksrc
                        ! within an iteration [nondim].  0 < tol_dksrc_low < 1.
  real :: Ri_crit       !   The critical shear Richardson number for shear-
                        ! driven mixing [nondim]. The theoretical value is 0.25.
  real :: dt_rem        !   The remaining time to advance the solution [T ~> s].
  real :: dt_now        !   The time step used in the current iteration [T ~> s].
  real :: dt_wt         !   The fractional weight of the current iteration [nondim].
  real :: dt_test       !   A time-step that is being tested for whether it
                        ! gives acceptably small changes in k_src [T ~> s].
  real :: Idtt          !   Idtt = 1 / dt_test [T-1 ~> s-1].
  real :: dt_inc        !   An increment to dt_test that is being tested [T ~> s].

  real :: k0dt          ! The background diffusivity times the timestep [Z2 ~> m2].
  logical :: valid_dt   ! If true, all levels so far exhibit acceptably small changes in k_src.
  logical :: use_temperature  !  If true, temperature and salinity have been
                        ! allocated and are being used as state variables.
  integer :: ks_kappa, ke_kappa  ! The k-range with nonzero kappas.
  integer :: dt_halvings   ! The number of times that the time-step is halved
                           ! in seeking an acceptable timestep.  If none is
                           ! found, dt_rem*0.5^dt_halvings is used.
  integer :: dt_refinements ! The number of 2-fold refinements that will be used
                           ! to estimate the maximum permitted time step.  I.e.,
                           ! the resolution is 1/2^dt_refinements.
  integer :: k, itt, itt_dt

  ! This calculation of N2 is for debugging only.
  ! real, dimension(SZK_(GV)+1) :: &
  !   N2_debug, & ! A version of N2 for debugging [T-2 ~> s-2]

  Ri_crit = CS%Rino_crit
  gR0 = GV%Rho0 * GV%g_Earth
  g_R0 = (US%L_to_Z**2 * GV%g_Earth) / (GV%Rho0)
  k0dt = dt*CS%kappa_0

  tol_dksrc = CS%kappa_src_max_chg
  if (tol_dksrc == 10.0) then
    ! This is equivalent to the expression below, but avoids changes at roundoff for the default value.
    tol_dksrc_low = 0.95
  else
    tol_dksrc_low = (tol_dksrc - 0.5)/tol_dksrc
  endif
  tol2 = 2.0*CS%kappa_tol_err
  dt_refinements = 5 ! Selected so that 1/2^dt_refinements < 1-tol_dksrc_low
  use_temperature = .false. ; if (associated(tv%T)) use_temperature = .true.


  ! Set up Idz as the inverse of layer thicknesses.
  do k=1,nzc ; Idz(k) = 1.0 / dz(k) ; enddo
  !   Set up I_dz_int as the inverse of the distance between
  ! adjacent layer centers.
  I_dz_int(1) = 2.0 / dz(1)
  dist_from_top(1) = 0.0
  do K=2,nzc
    I_dz_int(K) = 2.0 / (dz(k-1) + dz(k))
    dist_from_top(K) = dist_from_top(K-1) + dz(k-1)
  enddo
  I_dz_int(nzc+1) = 2.0 / dz(nzc)

  !   Determine the velocities and thicknesses after eliminating massless
  ! layers and applying a time-step of background diffusion.
  if (nzc > 1) then
    a1(2) = k0dt*I_dz_int(2)
    b1 = 1.0 / (dz(1) + a1(2))
    u(1) = b1 * u0xdz(1) ; v(1) = b1 * v0xdz(1)
    T(1) = b1 * T0xdz(1) ; Sal(1) = b1 * S0xdz(1)
    c1(2) = a1(2) * b1 ; d1 = dz(1) * b1 ! = 1 - c1
    do k=2,nzc-1
      bd1 = dz(k) + d1*a1(k)
      a1(k+1) = k0dt*I_dz_int(k+1)
      b1 = 1.0 / (bd1 + a1(k+1))
      u(k) = b1 * (u0xdz(k) + a1(k)*u(k-1))
      v(k) = b1 * (v0xdz(k) + a1(k)*v(k-1))
      T(k) = b1 * (T0xdz(k) + a1(k)*T(k-1))
      Sal(k) = b1 * (S0xdz(k) + a1(k)*Sal(k-1))
      c1(k+1) = a1(k+1) * b1 ; d1 = bd1 * b1 ! d1 = 1 - c1
    enddo
    ! rho or T and S have insulating boundary conditions, u & v use no-slip
    ! bottom boundary conditions (if kappa0 > 0).
    ! For no-slip bottom boundary conditions
    b1 = 1.0 / ((dz(nzc) + d1*a1(nzc)) + k0dt*I_dz_int(nzc+1))
    u(nzc) = b1 * (u0xdz(nzc) + a1(nzc)*u(nzc-1))
    v(nzc) = b1 * (v0xdz(nzc) + a1(nzc)*v(nzc-1))
    ! For insulating boundary conditions
    b1 = 1.0 / (dz(nzc) + d1*a1(nzc))
    T(nzc) = b1 * (T0xdz(nzc) + a1(nzc)*T(nzc-1))
    Sal(nzc) = b1 * (S0xdz(nzc) + a1(nzc)*Sal(nzc-1))
    do k=nzc-1,1,-1
      u(k) = u(k) + c1(k+1)*u(k+1) ; v(k) = v(k) + c1(k+1)*v(k+1)
      T(k) = T(k) + c1(k+1)*T(k+1) ; Sal(k) = Sal(k) + c1(k+1)*Sal(k+1)
    enddo
  else
    ! This is correct, but probably unnecessary.
    b1 = 1.0 / (dz(1) + k0dt*I_dz_int(2))
    u(1) = b1 * u0xdz(1) ; v(1) = b1 * v0xdz(1)
    b1 = 1.0 / dz(1)
    T(1) = b1 * T0xdz(1) ; Sal(1) = b1 * S0xdz(1)
  endif

  ! This uses half the harmonic mean of thicknesses to provide two estimates
  ! of the boundary between cells, and the inverse of the harmonic mean to
  ! weight the two estimates.  The net effect is that interfaces around thin
  ! layers have thin cells, and the total thickness adds up properly.
  ! The top- and bottom- interfaces have zero thickness, consistent with
  ! adding additional zero thickness layers.
  dz_Int(1) = 0.0 ; dz_Int(2) = dz(1)
  do K=2,nzc-1
    Norm = 1.0 / (dz(k)*(dz(k-1)+dz(k+1)) + 2.0*dz(k-1)*dz(k+1))
    dz_Int(K) = dz_Int(K) + dz(k) * ( ((dz(k)+dz(k+1)) * dz(k-1)) * Norm)
    dz_Int(K+1) = dz(k) * ( ((dz(k-1)+dz(k)) * dz(k+1)) * Norm)
  enddo
  dz_Int(nzc) = dz_Int(nzc) + dz(nzc) ; dz_Int(nzc+1) = 0.0

  dist_from_bot = 0.0
  do K=nzc,2,-1
    dist_from_bot = dist_from_bot + dz(k)
    I_L2_bdry(K) = (dist_from_top(K) + dist_from_bot)**2 / &
                   (dist_from_top(K) * dist_from_bot)**2
  enddo

  ! Calculate thermodynamic coefficients and an initial estimate of N2.
  if (use_temperature) then
    pressure(1) = surface_pres
    do K=2,nzc
      pressure(K) = pressure(K-1) + gR0*dz(k-1)
      T_int(K) = 0.5*(T(k-1) + T(k))
      Sal_int(K) = 0.5*(Sal(k-1) + Sal(k))
    enddo
    call calculate_density_derivs(T_int, Sal_int, pressure, dbuoy_dT, dbuoy_dS, &
                                  tv%eqn_of_state, (/2,nzc/), scale=-g_R0 )
  else
    do K=1,nzc+1 ; dbuoy_dT(K) = -g_R0 ; dbuoy_dS(K) = 0.0 ; enddo
  endif

  ! N2_debug(1) = 0.0 ; N2_debug(nzc+1) = 0.0
  ! do K=2,nzc
  !   N2_debug(K) = max((dbuoy_dT(K) * (T0xdz(k-1)*Idz(k-1) - T0xdz(k)*Idz(k)) + &
  !                      dbuoy_dS(K) * (S0xdz(k-1)*Idz(k-1) - S0xdz(k)*Idz(k))) * &
  !                      I_dz_int(K), 0.0)
  ! enddo

  ! This call just calculates N2 and S2.
  call calculate_projected_state(kappa, u, v, T, Sal, 0.0, nzc, dz, I_dz_int, &
                                 dbuoy_dT, dbuoy_dS, u, v, T, Sal, GV, US, &
                                 N2=N2, S2=S2, vel_underflow=CS%vel_underflow)
! ----------------------------------------------------
! Iterate
! ----------------------------------------------------
  dt_rem = dt
  do K=1,nzc+1
    K_Q(K) = 0.0
    kappa_avg(K) = 0.0 ; tke_avg(K) = 0.0
    local_src_avg(K) = 0.0
    ! Use the grid spacings to scale errors in the source.
    if ( dz_Int(K) > 0.0 ) &
      local_src_avg(K) = 0.1 * k0dt * I_dz_int(K) / dz_Int(K)
  enddo

! call cpu_clock_end(id_clock_setup)

! do itt=1,CS%max_RiNo_it
  do itt=1,CS%max_KS_it

! ----------------------------------------------------
! Calculate new values of u, v, rho, N^2 and S.
! ----------------------------------------------------

  ! call cpu_clock_begin(id_clock_KQ)
    call find_kappa_tke(N2, S2, kappa, Idz, dz_Int, I_L2_bdry, f2, &
                        nzc, CS, GV, US, K_Q, tke, kappa_out, kappa_src, local_src)
  ! call cpu_clock_end(id_clock_KQ)

  ! call cpu_clock_begin(id_clock_avg)
    ! Determine the range of non-zero values of kappa_out.
    ks_kappa = GV%ke+1 ; ke_kappa = 0
    do K=2,nzc ; if (kappa_out(K) > 0.0) then
      ks_kappa = K ; exit
    endif ; enddo
    do k=nzc,ks_kappa,-1 ; if (kappa_out(K) > 0.0) then
      ke_kappa = K ; exit
    endif ; enddo
    if (ke_kappa == nzc) kappa_out(nzc+1) = 0.0
  ! call cpu_clock_end(id_clock_avg)

    ! Determine how long to use this value of kappa (dt_now).

  ! call cpu_clock_begin(id_clock_project)
    if ((ke_kappa < ks_kappa) .or. (itt==CS%max_RiNo_it)) then
      dt_now = dt_rem
    else
      ! Limit dt_now so that |k_src(k)-kappa_src(k)| < tol * local_src(k)
      dt_test = dt_rem
      do K=2,nzc
        tol_max(K) = kappa_src(K) + tol_dksrc * local_src(K)
        tol_min(K) = kappa_src(K) - tol_dksrc_low * local_src(K)
        tol_chg(K) = tol2 * local_src_avg(K)
      enddo

      do itt_dt=1,(CS%max_KS_it+1-itt)/2
        !   The maximum number of times that the time-step is halved in
        ! seeking an acceptable timestep is reduced with each iteration,
        ! so that as the maximum number of iterations is approached, the
        ! whole remaining timestep is used.  Typically, an acceptable
        ! timestep is found long before the minimum is reached, so the
        ! value of max_KS_it may be unimportant, especially if it is large
        ! enough.
        call calculate_projected_state(kappa_out, u, v, T, Sal, 0.5*dt_test, nzc, dz, I_dz_int, &
                                       dbuoy_dT, dbuoy_dS, u_test, v_test, T_test, S_test, &
                                       GV, US, N2, S2, ks_int=ks_kappa, ke_int=ke_kappa, &
                                       vel_underflow=CS%vel_underflow)
        valid_dt = .true.
        Idtt = 1.0 / dt_test
        do K=max(ks_kappa-1,2),min(ke_kappa+1,nzc)
          if (N2(K) < Ri_crit * S2(K)) then ! Equivalent to Ri < Ri_crit.
            K_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
                       ((Ri_crit*S2(K) - N2(K)) / (Ri_crit*S2(K) + CS%FRi_curvature*N2(K)))
            if ((K_src(K) > max(tol_max(K), kappa_src(K) + Idtt*tol_chg(K))) .or. &
                (K_src(K) < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K)))) then
              valid_dt = .false. ; exit
            endif
          else
            if (0.0 < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K))) then
              valid_dt = .false. ; k_src(K) = 0.0 ; exit
            endif
          endif
        enddo

        if (valid_dt) exit
        dt_test = 0.5*dt_test
      enddo
      if ((dt_test < dt_rem) .and. valid_dt) then
        dt_inc = 0.5*dt_test
        do itt_dt=1,dt_refinements
          call calculate_projected_state(kappa_out, u, v, T, Sal, 0.5*(dt_test+dt_inc), &
                   nzc, dz, I_dz_int, dbuoy_dT, dbuoy_dS, u_test, v_test, T_test, S_test, &
                   GV, US, N2, S2, ks_int=ks_kappa, ke_int=ke_kappa, vel_underflow=CS%vel_underflow)
          valid_dt = .true.
          Idtt = 1.0 / (dt_test+dt_inc)
          do K=max(ks_kappa-1,2),min(ke_kappa+1,nzc)
            if (N2(K) < Ri_crit * S2(K)) then ! Equivalent to Ri < Ri_crit.
              K_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
                         ((Ri_crit*S2(K) - N2(K)) / (Ri_crit*S2(K) + CS%FRi_curvature*N2(K)))
              if ((K_src(K) > max(tol_max(K), kappa_src(K) + Idtt*tol_chg(K))) .or. &
                  (K_src(K) < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K)))) then
                valid_dt = .false. ; exit
              endif
            else
              if (0.0 < min(tol_min(K), kappa_src(K) - Idtt*tol_chg(K))) then
                valid_dt = .false. ; k_src(K) = 0.0 ; exit
              endif
            endif
          enddo

          if (valid_dt) dt_test = dt_test + dt_inc
          dt_inc = 0.5*dt_inc
        enddo
      else
        dt_inc = 0.0
      endif

      dt_now = min(dt_test*(1.0+CS%kappa_tol_err)+dt_inc, dt_rem)
      do K=2,nzc
        local_src_avg(K) = local_src_avg(K) + dt_now * local_src(K)
      enddo
    endif  ! Are all the values of kappa_out 0?
  ! call cpu_clock_end(id_clock_project)

    ! The state has already been projected forward. Now find new values of kappa.

    if (ke_kappa < ks_kappa) then
      ! There is no mixing now, and will not be again.
    ! call cpu_clock_begin(id_clock_avg)
      dt_wt = dt_rem / dt ; dt_rem = 0.0
      do K=1,nzc+1
        kappa_mid(K) = 0.0
        ! This would be here but does nothing.
        ! kappa_avg(K) = kappa_avg(K) + kappa_mid(K)*dt_wt
        tke_avg(K) = tke_avg(K) + dt_wt*tke(K)
      enddo
    ! call cpu_clock_end(id_clock_avg)
    else
    ! call cpu_clock_begin(id_clock_project)
      call calculate_projected_state(kappa_out, u, v, T, Sal, dt_now, nzc, dz, I_dz_int, &
                                     dbuoy_dT, dbuoy_dS, u_test, v_test, T_test, S_test, &
                                     GV, US, N2=N2, S2=S2, ks_int=ks_kappa, ke_int=ke_kappa, &
                                     vel_underflow=CS%vel_underflow)
    ! call cpu_clock_end(id_clock_project)

    ! call cpu_clock_begin(id_clock_KQ)
      do K=1,nzc+1 ; K_Q_tmp(K) = K_Q(K) ; enddo
      call find_kappa_tke(N2, S2, kappa_out, Idz, dz_Int, I_L2_bdry, f2, &
                          nzc, CS, GV, US, K_Q_tmp, tke_pred, kappa_pred)
    ! call cpu_clock_end(id_clock_KQ)

      ks_kappa = GV%ke+1 ; ke_kappa = 0
      do K=1,nzc+1
        kappa_mid(K) = 0.5*(kappa_out(K) + kappa_pred(K))
        if ((kappa_mid(K) > 0.0) .and. (K<ks_kappa)) ks_kappa = K
        if (kappa_mid(K) > 0.0) ke_kappa = K
      enddo

    ! call cpu_clock_begin(id_clock_project)
      call calculate_projected_state(kappa_mid, u, v, T, Sal, dt_now, nzc, dz, I_dz_int, &
                                     dbuoy_dT, dbuoy_dS, u_test, v_test, T_test, S_test, &
                                     GV, US, N2=N2, S2=S2, ks_int=ks_kappa, ke_int=ke_kappa, &
                                     vel_underflow=CS%vel_underflow)
    ! call cpu_clock_end(id_clock_project)

    ! call cpu_clock_begin(id_clock_KQ)
      call find_kappa_tke(N2, S2, kappa_out, Idz, dz_Int, I_L2_bdry, f2, &
                          nzc, CS, GV, US, K_Q, tke_pred, kappa_pred)
    ! call cpu_clock_end(id_clock_KQ)

    ! call cpu_clock_begin(id_clock_avg)
      dt_wt = dt_now / dt ; dt_rem = dt_rem - dt_now
      do K=1,nzc+1
        kappa_mid(K) = 0.5*(kappa_out(K) + kappa_pred(K))
        kappa_avg(K) = kappa_avg(K) + kappa_mid(K)*dt_wt
        tke_avg(K) = tke_avg(K) + dt_wt*0.5*(tke_pred(K) + tke(K))
        kappa(K) = kappa_pred(K) ! First guess for the next iteration.
      enddo
    ! call cpu_clock_end(id_clock_avg)
    endif

    if (dt_rem > 0.0) then
      ! Update the values of u, v, T, Sal, N2, and S2 for the next iteration.
    ! call cpu_clock_begin(id_clock_project)
      call calculate_projected_state(kappa_mid, u, v, T, Sal, dt_now, nzc, &
                                     dz, I_dz_int, dbuoy_dT, dbuoy_dS, u, v, T, Sal, &
                                     GV, US, N2, S2, vel_underflow=CS%vel_underflow)
    ! call cpu_clock_end(id_clock_project)
    endif

    if (dt_rem <= 0.0) exit

  enddo ! end itt loop

  if (present(I_Ld2_1d)) then
    do K=1,GV%ke+1 ; I_Ld2_1d(K) = 0.0 ; enddo
    do K=2,nzc ; if (TKE(K) > 0.0) &
      I_Ld2_1d(K) = I_L2_bdry(K) + (N2(K) / CS%lambda**2 + f2) / TKE(K)
    enddo
  endif
  if (present(dz_Int_1d)) then
    do K=1,nzc+1 ; dz_Int_1d(K) = dz_Int(K) ; enddo
    do K=nzc+2,GV%ke ; dz_Int_1d(K) = 0.0 ; enddo
  endif

end subroutine kappa_shear_column

!>   This subroutine calculates the velocities, temperature and salinity that
!! the water column will have after mixing for dt with diffusivities kappa.  It
!! may also calculate the projected buoyancy frequency and shear.
subroutine calculate_projected_state(kappa, u0, v0, T0, S0, dt, nz, &
                                     dz, I_dz_int, dbuoy_dT, dbuoy_dS, &
                                     u, v, T, Sal, GV, US, N2, S2, ks_int, ke_int, vel_underflow)
  integer,               intent(in)    :: nz  !< The number of layers (after eliminating massless
                                              !! layers?).
  real, dimension(nz+1), intent(in)    :: kappa !< The diapycnal diffusivity at interfaces,
                                              !! [Z2 T-1 ~> m2 s-1].
  real, dimension(nz),   intent(in)    :: u0  !< The initial zonal velocity [L T-1 ~> m s-1].
  real, dimension(nz),   intent(in)    :: v0  !< The initial meridional velocity [L T-1 ~> m s-1].
  real, dimension(nz),   intent(in)    :: T0  !< The initial temperature [degC].
  real, dimension(nz),   intent(in)    :: S0  !< The initial salinity [ppt].
  real, dimension(nz),   intent(in)    :: dz  !< The grid spacing of layers [Z ~> m].
  real, dimension(nz+1), intent(in)    :: I_dz_int !< The inverse of the layer's thicknesses
                                              !! [Z-1 ~> m-1].
  real, dimension(nz+1), intent(in)    :: dbuoy_dT !< The partial derivative of buoyancy with
                                              !! temperature [Z T-2 degC-1 ~> m s-2 degC-1].
  real, dimension(nz+1), intent(in)    :: dbuoy_dS !< The partial derivative of buoyancy with
                                              !! salinity [Z T-2 ppt-1 ~> m s-2 ppt-1].
  real,                  intent(in)    :: dt  !< The time step [T ~> s].
  real, dimension(nz),   intent(inout) :: u   !< The zonal velocity after dt [L T-1 ~> m s-1].
  real, dimension(nz),   intent(inout) :: v   !< The meridional velocity after dt [L T-1 ~> m s-1].
  real, dimension(nz),   intent(inout) :: T   !< The temperature after dt [degC].
  real, dimension(nz),   intent(inout) :: Sal !< The salinity after dt [ppt].
  type(verticalGrid_type), intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(nz+1), optional, &
                         intent(inout) :: N2  !< The buoyancy frequency squared at interfaces [T-2 ~> s-2].
  real, dimension(nz+1), optional, &
                         intent(inout) :: S2  !< The squared shear at interfaces [T-2 ~> s-2].
  integer, optional,     intent(in)    :: ks_int !< The topmost k-index with a non-zero diffusivity.
  integer, optional,     intent(in)    :: ke_int !< The bottommost k-index with a non-zero
                                              !! diffusivity.
  real,    optional,     intent(in)    :: vel_underflow !< If present and true, any velocities that
                                              !! are smaller in magnitude than this value are
                                              !! set to 0 [L T-1 ~> m s-1].

  ! Local variables
  real, dimension(nz+1) :: c1
  real :: L2_to_Z2       ! A conversion factor from horizontal length units to vertical depth
                         ! units squared [Z2 s2 T-2 m-2 ~> 1].
  real :: underflow_vel  ! Velocities smaller in magnitude than underflow_vel are set to 0 [L T-1 ~> m s-1].
  real :: a_a, a_b, b1, d1, bd1, b1nz_0
  integer :: k, ks, ke

  ks = 1 ; ke = nz
  if (present(ks_int)) ks = max(ks_int-1,1)
  if (present(ke_int)) ke = min(ke_int,nz)
  underflow_vel = 0.0 ; if (present(vel_underflow)) underflow_vel = vel_underflow

  if (ks > ke) return

  if (dt > 0.0) then
    a_b = dt*(kappa(ks+1)*I_dz_int(ks+1))
    b1 = 1.0 / (dz(ks) + a_b)
    c1(ks+1) = a_b * b1 ; d1 = dz(ks) * b1 ! = 1 - c1

    u(ks) = (b1 * dz(ks))*u0(ks) ; v(ks) = (b1 * dz(ks))*v0(ks)
    T(ks) = (b1 * dz(ks))*T0(ks) ; Sal(ks) = (b1 * dz(ks))*S0(ks)
    do K=ks+1,ke-1
      a_a = a_b
      a_b = dt*(kappa(K+1)*I_dz_int(K+1))
      bd1 = dz(k) + d1*a_a
      b1 = 1.0 / (bd1 + a_b)
      c1(K+1) = a_b * b1 ; d1 = bd1 * b1 ! d1 = 1 - c1

      u(k) = b1 * (dz(k)*u0(k) + a_a*u(k-1))
      v(k) = b1 * (dz(k)*v0(k) + a_a*v(k-1))
      T(k) = b1 * (dz(k)*T0(k) + a_a*T(k-1))
      Sal(k) = b1 * (dz(k)*S0(k) + a_a*Sal(k-1))
    enddo
    !   T and S have insulating boundary conditions, u & v use no-slip
    ! bottom boundary conditions at the solid bottom.

    ! For insulating boundary conditions or mixing simply stopping, use...
    a_a = a_b
    b1 = 1.0 / (dz(ke) + d1*a_a)
    T(ke) = b1 * (dz(ke)*T0(ke) + a_a*T(ke-1))
    Sal(ke) = b1 * (dz(ke)*S0(ke) + a_a*Sal(ke-1))

    !   There is no distinction between the effective boundary conditions for
    ! tracers and velocities if the mixing is separated from the bottom, but if
    ! the mixing goes all the way to the bottom, use no-slip BCs for velocities.
    if (ke == nz) then
      a_b = dt*(kappa(nz+1)*I_dz_int(nz+1))
      b1nz_0 = 1.0 / ((dz(nz) + d1*a_a) + a_b)
    else
      b1nz_0 = b1
    endif
    u(ke) = b1nz_0 * (dz(ke)*u0(ke) + a_a*u(ke-1))
    v(ke) = b1nz_0 * (dz(ke)*v0(ke) + a_a*v(ke-1))
    if (abs(u(ke)) < underflow_vel) u(ke) = 0.0
    if (abs(v(ke)) < underflow_vel) v(ke) = 0.0

    do k=ke-1,ks,-1
      u(k) = u(k) + c1(k+1)*u(k+1)
      v(k) = v(k) + c1(k+1)*v(k+1)
      if (abs(u(k)) < underflow_vel) u(k) = 0.0
      if (abs(v(k)) < underflow_vel) v(k) = 0.0
      T(k) = T(k) + c1(k+1)*T(k+1)
      Sal(k) = Sal(k) + c1(k+1)*Sal(k+1)
    enddo
  else ! dt <= 0.0
    do k=1,nz
      u(k) = u0(k) ; v(k) = v0(k) ; T(k) = T0(k) ; Sal(k) = S0(k)
      if (abs(u(k)) < underflow_vel) u(k) = 0.0
      if (abs(v(k)) < underflow_vel) v(k) = 0.0
    enddo
  endif

  if (present(S2)) then
    ! L2_to_Z2 = US%m_to_Z**2 * US%T_to_s**2
    L2_to_Z2 = US%L_to_Z**2
    S2(1) = 0.0 ; S2(nz+1) = 0.0
    if (ks > 1) &
      S2(ks) = ((u(ks)-u0(ks-1))**2 + (v(ks)-v0(ks-1))**2) * (L2_to_Z2*I_dz_int(ks)**2)
    do K=ks+1,ke
      S2(K) = ((u(k)-u(k-1))**2 + (v(k)-v(k-1))**2) * (L2_to_Z2*I_dz_int(K)**2)
    enddo
    if (ke<nz) &
      S2(ke+1) = ((u0(ke+1)-u(ke))**2 + (v0(ke+1)-v(ke))**2) * (L2_to_Z2*I_dz_int(ke+1)**2)
  endif

  if (present(N2)) then
    N2(1) = 0.0 ; N2(nz+1) = 0.0
    if (ks > 1) &
      N2(ks) = max(0.0, I_dz_int(ks) * &
        (dbuoy_dT(ks) * (T0(ks-1)-T(ks)) + dbuoy_dS(ks) * (S0(ks-1)-Sal(ks))))
    do K=ks+1,ke
      N2(K) = max(0.0, I_dz_int(K) * &
        (dbuoy_dT(K) * (T(k-1)-T(k)) + dbuoy_dS(K) * (Sal(k-1)-Sal(k))))
    enddo
    if (ke<nz) &
      N2(ke+1) = max(0.0, I_dz_int(ke+1) * &
        (dbuoy_dT(ke+1) * (T(ke)-T0(ke+1)) + dbuoy_dS(ke+1) * (Sal(ke)-S0(ke+1))))
  endif

end subroutine calculate_projected_state

!> This subroutine calculates new, consistent estimates of TKE and kappa.
subroutine find_kappa_tke(N2, S2, kappa_in, Idz, dz_Int, I_L2_bdry, f2, &
                          nz, CS, GV, US, K_Q, tke, kappa, kappa_src, local_src)
  integer,               intent(in)    :: nz  !< The number of layers to work on.
  real, dimension(nz+1), intent(in)    :: N2  !< The buoyancy frequency squared at interfaces [T-2 ~> s-2].
  real, dimension(nz+1), intent(in)    :: S2  !< The squared shear at interfaces [T-2 ~> s-2].
  real, dimension(nz+1), intent(in)    :: kappa_in  !< The initial guess at the diffusivity
                                              !! [Z2 T-1 ~> m2 s-1].
  real, dimension(nz+1), intent(in)    :: dz_Int !< The thicknesses associated with interfaces
                                              !! [Z-1 ~> m-1].
  real, dimension(nz+1), intent(in)    :: I_L2_bdry !< The inverse of the squared distance to
                                              !! boundaries [Z-2 ~> m-2].
  real, dimension(nz),   intent(in)    :: Idz !< The inverse grid spacing of layers [Z-1 ~> m-1].
  real,                  intent(in)    :: f2  !< The squared Coriolis parameter [T-2 ~> s-2].
  type(Kappa_shear_CS),  pointer       :: CS  !< A pointer to this module's control structure.
  type(verticalGrid_type), intent(in)  :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type), intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(nz+1), intent(inout) :: K_Q !< The shear-driven diapycnal diffusivity divided by
                                              !! the turbulent kinetic energy per unit mass at
                                              !! interfaces [T ~> s].
  real, dimension(nz+1), intent(out)   :: tke !< The turbulent kinetic energy per unit mass at
                                              !! interfaces [Z2 T-2 ~> m2 s-2].
  real, dimension(nz+1), intent(out)   :: kappa  !< The diapycnal diffusivity at interfaces
                                              !! [Z2 T-1 ~> m2 s-1].
  real, dimension(nz+1), optional, &
                         intent(out)   :: kappa_src !< The source term for kappa [T-1 ~> s-1].
  real, dimension(nz+1), optional, &
                         intent(out)   :: local_src !< The sum of all local sources for kappa,
                                              !! [T-1 ~> s-1].
!   This subroutine calculates new, consistent estimates of TKE and kappa.

  ! Local variables
  real, dimension(nz) :: &
    aQ, &       ! aQ is the coupling between adjacent interfaces in the TKE equations [Z T-1 ~> m s-1].
    dQdz        ! Half the partial derivative of TKE with depth [Z T-2 ~> m s-2].
  real, dimension(nz+1) :: &
    dK, &         ! The change in kappa [Z2 T-1 ~> m2 s-1].
    dQ, &         ! The change in TKE [Z2 T-2 ~> m2 s-2].
    cQ, cK, &     ! cQ and cK are the upward influences in the tridiagonal and
                  ! hexadiagonal solvers for the TKE and kappa equations [nondim].
    I_Ld2, &      ! 1/Ld^2, where Ld is the effective decay length scale for kappa [Z-2 ~> m-2].
    TKE_decay, &  ! The local TKE decay rate [T-1 ~> s-1].
    k_src, &      ! The source term in the kappa equation [T-1 ~> s-1].
    dQmdK, &      ! With Newton's method the change in dQ(k-1) due to dK(k) [T ~> s].
    dKdQ, &       ! With Newton's method the change in dK(k) due to dQ(k) [T-1 ~> s-1].
    e1            ! The fractional change in a layer TKE due to a change in the
                  ! TKE of the layer above when all the kappas below are 0.
                  ! e1 is nondimensional, and 0 < e1 < 1.
  real :: tke_src       ! The net source of TKE due to mixing against the shear
                        ! and stratification [Z2 T-3 ~> m2 s-3].  (For convenience,
                        ! a term involving the non-dissipation of q0 is also
                        ! included here.)
  real :: bQ            ! The inverse of the pivot in the tridiagonal equations [T Z-1 ~> s m-1].
  real :: bK            ! The inverse of the pivot in the tridiagonal equations [Z-1 ~> m-1].
  real :: bQd1          ! A term in the denominator of bQ [Z T-1 ~> m s-1].
  real :: bKd1          ! A term in the denominator of bK [Z ~> m].
  real :: cQcomp, cKcomp ! 1 - cQ or 1 - cK in the tridiagonal equations.
  real :: c_s2          !   The coefficient for the decay of TKE due to
                        ! shear (i.e. proportional to |S|*tke), nondimensional.
  real :: c_n2          !   The coefficient for the decay of TKE due to
                        ! stratification (i.e. proportional to N*tke) [nondim].
  real :: Ri_crit       !   The critical shear Richardson number for shear-
                        ! driven mixing. The theoretical value is 0.25.
  real :: q0            !   The background level of TKE [Z2 T-2 ~> m2 s-2].
  real :: Ilambda2      ! 1.0 / CS%lambda**2 [nondim]
  real :: TKE_min       !   The minimum value of shear-driven TKE that can be
                        ! solved for [Z2 T-2 ~> m2 s-2].
  real :: kappa0        ! The background diapycnal diffusivity [Z2 T-1 ~> m2 s-1].
  real :: kappa_trunc   ! Diffusivities smaller than this are rounded to 0 [Z2 T-1 ~> m2 s-1].

  real :: eden1, eden2, I_eden, ome  ! Variables used in calculating e1.
  real :: diffusive_src ! The diffusive source in the kappa equation [Z T-1 ~> m s-1].
  real :: chg_by_k0     ! The value of k_src that leads to an increase of
                        ! kappa_0 if only the diffusive term is a sink [T-1 ~> s-1].

  real :: kappa_mean    ! A mean value of kappa [Z2 T-1 ~> m2 s-1].
  real :: Newton_test   ! The value of relative error that will cause the next
                        ! iteration to use Newton's method.
  ! Temporary variables used in the Newton's method iterations.
  real :: decay_term_k  ! The decay term in the diffusivity equation
  real :: decay_term_Q  ! The decay term in the TKE equation - proportional to [T-1 ~> s-1]
  real :: I_Q           ! The inverse of TKE [T2 Z-2 ~> s2 m-2]
  real :: kap_src
  real :: v1            ! A temporary variable proportional to [T-1 ~> s-1]
  real :: v2
  real :: tol_err        ! The tolerance for max_err that determines when to
                         ! stop iterating.
  real :: Newton_err     ! The tolerance for max_err that determines when to
                         ! start using Newton's method.  Empirically, an initial
                         ! value of about 0.2 seems to be most efficient.
  real, parameter :: roundoff = 1.0e-16 ! A negligible fractional change in TKE.
                         ! This could be larger but performance gains are small.

  logical :: tke_noflux_bottom_BC = .false. ! Specify the boundary conditions
  logical :: tke_noflux_top_BC = .false.    ! that are applied to the TKE eqns.
  logical :: do_Newton    ! If .true., use Newton's method for the next iteration.
  logical :: abort_Newton ! If .true., an Newton's method has encountered a 0
                          ! pivot, and should not have been used.
  logical :: was_Newton   ! The value of do_Newton before checking convergence.
  logical :: within_tolerance ! If .true., all points are within tolerance to
                          ! enable this subroutine to return.
  integer :: ks_src, ke_src ! The range indices that have nonzero k_src.
  integer :: ks_kappa, ke_kappa, ke_tke   ! The ranges of k-indices that are or
  integer :: ks_kappa_prev, ke_kappa_prev ! were being worked on.
  integer :: itt, k, k2

  ! These variables are used only for debugging.
  logical, parameter :: debug_soln = .false.
  real :: K_err_lin, Q_err_lin, TKE_src_norm
  real, dimension(nz+1) :: &
    I_Ld2_debug, & ! A separate version of I_Ld2 for debugging [Z-2 ~> m-2].
    kappa_prev, & ! The value of kappa at the start of the current iteration [Z2 T-1 ~> m2 s-1].
    TKE_prev   ! The value of TKE at the start of the current iteration [Z2 T-2 ~> m2 s-2].

  c_N2 = CS%C_N**2 ; c_S2 = CS%C_S**2
  q0 = CS%TKE_bg ; kappa0 = CS%kappa_0
  TKE_min = max(CS%TKE_bg, 1.0E-20*US%m_to_Z**2*US%T_to_s**2)
  Ri_crit = CS%Rino_crit
  Ilambda2 = 1.0 / CS%lambda**2
  kappa_trunc = CS%kappa_trunc
  do_Newton = .false. ; abort_Newton = .false.
  tol_err = CS%kappa_tol_err
  Newton_err = 0.2     ! This initial value may be automatically reduced later.

  ks_kappa = 2 ; ke_kappa = nz ; ks_kappa_prev = 2 ; ke_kappa_prev = nz

  ke_src = 0 ; ks_src = nz+1
  do K=2,nz
    if (N2(K) < Ri_crit * S2(K)) then ! Equivalent to Ri < Ri_crit.
!       Ri = N2(K) / S2(K)
!       k_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
!                  ((Ri_crit - Ri) / (Ri_crit + CS%FRi_curvature*Ri))
      K_src(K) = (2.0 * CS%Shearmix_rate * sqrt(S2(K))) * &
                 ((Ri_crit*S2(K) - N2(K)) / (Ri_crit*S2(K) + CS%FRi_curvature*N2(K)))
      ke_src = K
      if (ks_src > k) ks_src = K
    else
      k_src(K) = 0.0
    endif
  enddo

  ! If there is no source anywhere, return kappa(K) = 0.
  if (ks_src > ke_src) then
    do K=1,nz+1
      kappa(K) = 0.0 ; K_Q(K) = 0.0 ; tke(K) = TKE_min
    enddo
    if (present(kappa_src)) then ; do K=1,nz+1 ; kappa_src(K) = 0.0 ; enddo ; endif
    if (present(local_src)) then ; do K=1,nz+1 ; local_src(K) = 0.0 ; enddo ; endif
    return
  endif

  do K=1,nz+1
    kappa(K) = kappa_in(K)
!     TKE_decay(K) = c_n*sqrt(N2(K)) + c_s*sqrt(S2(K)) ! The expression in JHL.
    TKE_decay(K) = sqrt(c_n2*N2(K) + c_s2*S2(K))
    if ((kappa(K) > 0.0) .and. (K_Q(K) > 0.0)) then
      TKE(K) = kappa(K) / K_Q(K) ! Perhaps take the max with TKE_min
    else
      TKE(K) = TKE_min
    endif
  enddo
  ! Apply boundary conditions to kappa.
  kappa(1) = 0.0 ; kappa(nz+1) = 0.0

  ! Calculate the term (e1) that allows changes in TKE to be calculated quickly
  ! below the deepest nonzero value of kappa.  If kappa = 0, below interface
  ! k-1, the final changes in TKE are related by dQ(K+1) = e1(K+1)*dQ(K).
  eden2 = kappa0 * Idz(nz)
  if (tke_noflux_bottom_BC) then
    eden1 = dz_Int(nz+1)*TKE_decay(nz+1)
    I_eden = 1.0 / (eden2 + eden1)
    e1(nz+1) = eden2 * I_eden ; ome = eden1 * I_eden
  else
    e1(nz+1) = 0.0 ; ome = 1.0
  endif
  do k=nz,2,-1
    eden1 = dz_Int(K)*TKE_decay(K) + ome * eden2
    eden2 = kappa0 * Idz(k-1)
    I_eden = 1.0 / (eden2 + eden1)
    e1(K) = eden2 * I_eden ; ome = eden1 * I_eden ! = 1-e1
  enddo
  e1(1) = 0.0


  ! Iterate here to convergence to within some tolerance of order tol_err.
  do itt=1,CS%max_RiNo_it

  ! ----------------------------------------------------
  ! Calculate TKE
  ! ----------------------------------------------------

    if (debug_soln) then ; do K=1,nz+1 ; kappa_prev(K) = kappa(K) ; TKE_prev(K) = TKE(K) ; enddo ; endif

    if (.not.do_Newton) then
      !   Use separate steps of the TKE and kappa equations, that are
      ! explicit in the nonlinear source terms, implicit in a linearized
      ! version of the nonlinear sink terms, and implicit in the linear
      ! terms.

      ke_tke = max(ke_kappa,ke_kappa_prev)+1
      ! aQ is the coupling between adjacent interfaces [Z T-1 ~> m s-1].
      do k=1,min(ke_tke,nz)
        aQ(k) = (0.5*(kappa(K)+kappa(K+1)) + kappa0) * Idz(k)
      enddo
      dQ(1) = -TKE(1)
      if (tke_noflux_top_BC) then
        tke_src = kappa0*S2(1) + q0 * TKE_decay(1) ! Uses that kappa(1) = 0
        bQd1 = dz_Int(1) * TKE_decay(1)
        bQ = 1.0 / (bQd1 +  aQ(1))
        tke(1) = bQ * (dz_Int(1)*tke_src)
        cQ(2) = aQ(1) * bQ ; cQcomp = bQd1 * bQ ! = 1 - cQ
      else
        tke(1) = q0 ; cQ(2) = 0.0 ; cQcomp = 1.0
      endif
      do K=2,ke_tke-1
        dQ(K) = -TKE(K)
        tke_src = (kappa(K) + kappa0)*S2(K) + q0*TKE_decay(K)
        bQd1 = dz_Int(K)*(TKE_decay(K) + N2(K)*K_Q(K)) + cQcomp*aQ(k-1)
        bQ = 1.0 / (bQd1 + aQ(k))
        tke(K) = bQ * (dz_Int(K)*tke_src + aQ(k-1)*tke(K-1))
        cQ(K+1) = aQ(k) * bQ ; cQcomp = bQd1 * bQ ! = 1 - cQ
      enddo
      if ((ke_tke == nz+1) .and. .not.(tke_noflux_bottom_BC)) then
        tke(nz+1) = TKE_min
        dQ(nz+1) = 0.0
      else
        k = ke_tke
        tke_src = kappa0*S2(K) + q0*TKE_decay(K) ! Uses that kappa(ke_tke) = 0
        if (K == nz+1) then
          dQ(K) = -TKE(K)
          bQ = 1.0 / (dz_Int(K)*TKE_decay(K) + cQcomp*aQ(k-1))
          tke(K) = max(TKE_min, bQ * (dz_Int(K)*tke_src + aQ(k-1)*tke(K-1)))
          dQ(K) = tke(K) + dQ(K)
        else
          bQ = 1.0 / ((dz_Int(K)*TKE_decay(K) + cQcomp*aQ(k-1)) + aQ(k))
          cQ(K+1) = aQ(k) * bQ
          ! Account for all changes deeper in the water column.
          dQ(K) = -TKE(K)
          tke(K) = max((bQ * (dz_Int(K)*tke_src + aQ(k-1)*tke(K-1)) + &
                        cQ(K+1)*(tke(K+1) - e1(K+1)*tke(K))) / (1.0 - cQ(K+1)*e1(K+1)), TKE_min)
          dQ(K) = tke(K) + dQ(K)

          ! Adjust TKE deeper in the water column in case ke_tke increases.
          ! This might not be strictly necessary?
          do K=ke_tke+1,nz+1
            dQ(K) = e1(K)*dQ(K-1)
            tke(K) = max(tke(K) + dQ(K), TKE_min)
            if (abs(dQ(K)) < roundoff*tke(K)) exit
          enddo
          do K2=K+1,nz
            if (dQ(K2) == 0.0) exit
            dQ(K2) = 0.0
          enddo
        endif
      endif
      do K=ke_tke-1,1,-1
        tke(K) = max(tke(K) + cQ(K+1)*tke(K+1), TKE_min)
        dQ(K) = tke(K) + dQ(K)
      enddo

  ! ----------------------------------------------------
  ! Calculate kappa, here defined at interfaces.
  ! ----------------------------------------------------

      ke_kappa_prev = ke_kappa ; ks_kappa_prev = ks_kappa

      dK(1) = 0.0 ! kappa takes boundary values of 0.
      cK(2) = 0.0 ; cKcomp = 1.0
      if (itt == 1) then ; dO K=2,nz
        I_Ld2(K) = (N2(K)*Ilambda2 + f2) / tke(K) + I_L2_bdry(K)
      enddo ; endif
      do K=2,nz
        dK(K) = -kappa(K)
        if (itt>1) &
          I_Ld2(K) = (N2(K)*Ilambda2 + f2) / tke(K) + I_L2_bdry(K)
        bKd1 = dz_Int(K)*I_Ld2(K) + cKcomp*Idz(k-1)
        bK = 1.0 / (bKd1 + Idz(k))

        kappa(K) = bK * (Idz(k-1)*kappa(K-1) + dz_Int(K) * K_src(K))
        cK(K+1) = Idz(k) * bK ; cKcomp = bKd1 * bK ! = 1 - cK(K+1)

        ! Neglect values that are smaller than kappa_trunc.
        if (kappa(K) < cKcomp*kappa_trunc) then
          kappa(K) = 0.0
          if (K > ke_src) then ; ke_kappa = k-1 ; K_Q(K) = 0.0 ; exit ; endif
        elseif (kappa(K) < 2.0*cKcomp*kappa_trunc) then
          kappa(K) = 2.0 * (kappa(K) - cKcomp*kappa_trunc)
        endif
      enddo
      K_Q(ke_kappa) = kappa(ke_kappa) / tke(ke_kappa)
      dK(ke_kappa) = dK(ke_kappa) + kappa(ke_kappa)
      do K=ke_kappa+2,ke_kappa_prev
        dK(K) = -kappa(K) ; kappa(K) = 0.0 ; K_Q(K) = 0.0
      enddo
      do K=ke_kappa-1,2,-1
        kappa(K) = kappa(K) + cK(K+1)*kappa(K+1)
        ! Neglect values that are smaller than kappa_trunc.
        if (kappa(K) <= kappa_trunc) then
          kappa(K) = 0.0
          if (K < ks_src) then ; ks_kappa = k+1 ; K_Q(K) = 0.0 ; exit ; endif
        elseif (kappa(K) < 2.0*kappa_trunc) then
          kappa(K) = 2.0 * (kappa(K) - kappa_trunc)
        endif

        dK(K) = dK(K) + kappa(K)
        K_Q(K) = kappa(K) / tke(K)
      enddo
      do K=ks_kappa_prev,ks_kappa-2 ; kappa(K) = 0.0 ; K_Q(K) = 0.0 ; enddo

    else ! do_Newton is .true.
!   Once the solutions are close enough, use a Newton's method solver of the
!  whole system to accelerate convergence.
      ks_kappa_prev = ks_kappa ; ke_kappa_prev = ke_kappa ; ke_kappa = nz
      ks_kappa = 2
      dK(1) = 0.0 ; cK(2) = 0.0 ; cKcomp = 1.0 ; dKdQ(1) = 0.0
      aQ(1) = (0.5*(kappa(1)+kappa(2))+kappa0) * Idz(1)
      dQdz(1) = 0.5*(TKE(1) - TKE(2))*Idz(1)
      if (tke_noflux_top_BC) then
        tke_src = dz_Int(1) * (kappa0*S2(1) - (TKE(1) - q0)*TKE_decay(1)) - &
                  aQ(1) * (TKE(1) - TKE(2))

        bQ = 1.0 / (aQ(1) + dz_Int(1)*TKE_decay(1))
        cQ(2) = aQ(1) * bQ
        cQcomp = (dz_Int(1)*TKE_decay(1)) * bQ ! = 1 - cQ(2)
        dQmdK(2) = -dQdz(1) * bQ
        dQ(1) = bQ * tke_src
      else
        dQ(1) = 0.0 ; cQ(2) = 0.0 ; cQcomp = 1.0 ; dQmdK(2) = 0.0
      endif
      do K=2,nz
        I_Q = 1.0 / TKE(K)
        I_Ld2(K) = (N2(K)*Ilambda2 + f2) * I_Q + I_L2_bdry(K)

        kap_src = dz_Int(K) * (K_src(K) - I_Ld2(K)*kappa(K)) + &
                            Idz(k-1)*(kappa(K-1)-kappa(K)) - Idz(k)*(kappa(K)-kappa(K+1))

        ! Ensure that the pivot is always positive, and that 0 <= cK <= 1.
        ! Otherwise do not use Newton's method.
        decay_term_k = -Idz(k-1)*dQmdK(K)*dKdQ(K-1) + dz_Int(K)*I_Ld2(K)
        if (decay_term_k < 0.0) then ; abort_Newton = .true. ; exit ; endif
        bK = 1.0 / (Idz(k) + Idz(k-1)*cKcomp + decay_term_k)

        cK(K+1) = bK * Idz(k)
        cKcomp = bK * (Idz(k-1)*cKcomp + decay_term_k) ! = 1-cK(K+1)
        if (CS%dKdQ_iteration_bug) then
          dKdQ(K) = bK * (Idz(k-1)*dKdQ(K-1)*cQ(K) + &
                      US%m_to_Z*(N2(K)*Ilambda2 + f2) * I_Q**2 * kappa(K) )
        else
          dKdQ(K) = bK * (Idz(k-1)*dKdQ(K-1)*cQ(K) + &
                      dz_Int(K)*(N2(K)*Ilambda2 + f2) * I_Q**2 * kappa(K) )
        endif
        dK(K) = bK * (kap_src + Idz(k-1)*dK(K-1) + Idz(k-1)*dKdQ(K-1)*dQ(K-1))

        ! Truncate away negligibly small values of kappa.
        if (dK(K) <= cKcomp*(kappa_trunc - kappa(K))) then
          dK(K) = -cKcomp*kappa(K)
!         if (K > ke_src) then ; ke_kappa = k-1 ; K_Q(K) = 0.0 ; exit ; endif
        elseif (dK(K) < cKcomp*(2.0*kappa_trunc - kappa(K))) then
          dK(K) = 2.0 * dK(K) - cKcomp*(2.0*kappa_trunc - kappa(K))
        endif

        ! Solve for dQ(K)...
        aQ(k) = (0.5*(kappa(K)+kappa(K+1))+kappa0) * Idz(k)
        dQdz(k) = 0.5*(TKE(K) - TKE(K+1))*Idz(k)
        tke_src = dz_Int(K) * (((kappa(K) + kappa0)*S2(K) - kappa(k)*N2(K)) - &
                                  (TKE(k) - q0)*TKE_decay(k)) - &
                  (aQ(k) * (TKE(K) - TKE(K+1)) - aQ(k-1) * (TKE(K-1) - TKE(K)))
        v1 = aQ(k-1) + dQdz(k-1)*dKdQ(K-1)
        v2 = (v1*dQmdK(K) + dQdz(k-1)*cK(K)) + &
             ((dQdz(k-1) - dQdz(k)) + dz_Int(K)*(S2(K) - N2(K)))

        ! Ensure that the pivot is always positive, and that 0 <= cQ <= 1.
        ! Otherwise do not use Newton's method.
        decay_term_Q = dz_Int(K)*TKE_decay(K) - dQdz(k-1)*dKdQ(K-1)*cQ(K) - v2*dKdQ(K)
        if (decay_term_Q < 0.0) then ; abort_Newton = .true. ; exit ; endif
        bQ = 1.0 / (aQ(k) + (cQcomp*aQ(k-1) + decay_term_Q))

        cQ(K+1) = aQ(k) * bQ
        cQcomp = (cQcomp*aQ(k-1) + decay_term_Q) * bQ
        dQmdK(K+1) = (v2 * cK(K+1) - dQdz(k)) * bQ

        ! Ensure that TKE+dQ will not drop below 0.5*TKE.
        dQ(K) = max(bQ * ((v1 * dQ(K-1) + dQdz(k-1)*dK(k-1)) + &
                          (v2 * dK(K) + tke_src)), cQcomp*(-0.5*TKE(K)))

        ! Check whether the next layer will be affected by any nonzero kappas.
        if ((itt > 1) .and. (K > ke_src) .and. (dK(K) == 0.0) .and. &
            ((kappa(K) + kappa(K+1)) == 0.0)) then
        ! Could also do  .and. (bQ*abs(tke_src) < roundoff*TKE(K)) then
          ke_kappa = k-1 ; exit
        endif
      enddo
      if ((ke_kappa == nz) .and. (.not. abort_Newton)) then
        dK(nz+1) = 0.0 ; dKdQ(nz+1) = 0.0
        if (tke_noflux_bottom_BC) then
          K = nz+1
          tke_src = dz_Int(K) * (kappa0*S2(K) - (TKE(K) - q0)*TKE_decay(K)) + &
                    aQ(k-1) * (TKE(K-1) - TKE(K))

          v1 = aQ(k-1) + dQdz(k-1)*dKdQ(K-1)
          decay_term_Q = max(0.0, dz_Int(K)*TKE_decay(K) - dQdz(k-1)*dKdQ(K-1)*cQ(K))
          if (decay_term_Q < 0.0) then
            abort_Newton = .true.
          else
            bQ = 1.0 / (aQ(k) + (cQcomp*aQ(k-1) + decay_term_Q))
          ! Ensure that TKE+dQ will not drop below 0.5*TKE.
            dQ(K) = max(bQ * ((v1 * dQ(K-1) + dQdz(k-1)*dK(K-1)) + tke_src), -0.5*TKE(K))
            TKE(K) = max(TKE(K) + dQ(K), TKE_min)
          endif
        else
          dQ(nz+1) = 0.0
        endif
      elseif (.not. abort_Newton) then
        ! Alter the first-guess determination of dQ(K).
        dQ(ke_kappa+1) = dQ(ke_kappa+1) / (1.0 - cQ(ke_kappa+2)*e1(ke_kappa+2))
        TKE(ke_kappa+1) = max(TKE(ke_kappa+1) + dQ(ke_kappa+1), TKE_min)
        do k=ke_kappa+2,nz+1
          if (debug_soln .and. (K < nz+1)) then
          ! Ignore this source?
            aQ(k) = (0.5*(kappa(K)+kappa(K+1))+kappa0) * Idz(k)
        !    tke_src_norm = (dz_Int(K) * (kappa0*S2(K) - (TKE(K)-q0)*TKE_decay(K)) - &
        !                   (aQ(k) * (TKE(K) - TKE(K+1)) - aQ(k-1) * (TKE(K-1) - TKE(K))) ) / &
        !                   (aQ(k) + (aQ(k-1) + dz_Int(K)*TKE_decay(K)))
          endif
          dK(K) = 0.0
        ! Ensure that TKE+dQ will not drop below 0.5*TKE.
          dQ(K) = max(e1(K)*dQ(K-1),-0.5*TKE(K))
          TKE(K) = max(TKE(K) + dQ(K), TKE_min)
          if (abs(dQ(K)) < roundoff*TKE(K)) exit
        enddo
        if (debug_soln) then ; do K2=K+1,nz+1 ; dQ(K2) = 0.0 ; dK(K2) = 0.0 ; enddo ; endif
      endif
      if (.not. abort_Newton) then
        do K=ke_kappa,2,-1
          ! Ensure that TKE+dQ will not drop below 0.5*TKE.
          dQ(K) = max(dQ(K) + (cQ(K+1)*dQ(K+1) + dQmdK(K+1) * dK(K+1)), -0.5*TKE(K))
          TKE(K) = max(TKE(K) + dQ(K), TKE_min)
          dK(K) = dK(K) + (cK(K+1)*dK(K+1) + dKdQ(K) * dQ(K))
          ! Truncate away negligibly small values of kappa.
          if (dK(K) <= kappa_trunc - kappa(K)) then
            dK(K) = -kappa(K)
            kappa(K) = 0.0
            if ((K < ks_src) .and. (K+1 > ks_kappa)) ks_kappa = K+1
          elseif (dK(K) < 2.0*kappa_trunc - kappa(K)) then
            dK(K) =  2.0*dK(K) - (2.0*kappa_trunc - kappa(K))
            kappa(K) = max(kappa(K) + dK(K), 0.0) ! The max is for paranoia.
            if (K<=ks_kappa) ks_kappa = 2
          else
            kappa(K) = kappa(K) + dK(K)
            if (K<=ks_kappa) ks_kappa = 2
          endif
        enddo
        dQ(1) = max(dQ(1) + cQ(2)*dQ(2) + dQmdK(2) * dK(2), TKE_min - TKE(1))
        TKE(1) = max(TKE(1) + dQ(1), TKE_min)
        dK(1) = 0.0
      endif

      ! Check these solutions for consistency.
      !  The unit conversions here have not been carefully tested.
      if (debug_soln) then ; do K=2,nz
        ! In these equations, K_err_lin and Q_err_lin should be at round-off levels
        ! compared with the dominant terms, perhaps, dz_Int*I_Ld2*kappa and
        ! dz_Int*TKE_decay*TKE.  The exception is where, either 1) the decay term has been
        ! been increased to ensure a positive pivot, or 2) negative TKEs have been
        ! truncated, or 3) small or negative kappas have been rounded toward 0.
        I_Q = 1.0 / TKE(K)
        I_Ld2_debug(K) = (N2(K)*Ilambda2 + f2) * I_Q + I_L2_bdry(K)

        kap_src = dz_Int(K) * (K_src(K) - I_Ld2(K)*kappa_prev(K)) + &
                            (Idz(k-1)*(kappa_prev(k-1)-kappa_prev(k)) - &
                             Idz(k)*(kappa_prev(k)-kappa_prev(k+1)))
        K_err_lin = -Idz(k-1)*(dK(K-1)-dK(K)) + Idz(k)*(dK(K)-dK(K+1)) + &
                     dz_Int(K)*I_Ld2_debug(K)*dK(K) - kap_src - &
                     dz_Int(K)*(N2(K)*Ilambda2 + f2)*I_Q**2*kappa_prev(K) * dQ(K)

        tke_src = dz_Int(K) * ((kappa_prev(K) + kappa0)*S2(K) - &
                     kappa_prev(K)*N2(K) - (TKE_prev(K) - q0)*TKE_decay(K)) - &
                  (aQ(k) * (TKE_prev(K) - TKE_prev(K+1)) -  aQ(k-1) * (TKE_prev(K-1) - TKE_prev(K)))
        Q_err_lin = tke_src + (aQ(k-1) * (dQ(K-1)-dQ(K)) - aQ(k) * (dQ(k)-dQ(k+1))) - &
                    0.5*(TKE_prev(K)-TKE_prev(K+1))*Idz(k)  * (dK(K) + dK(K+1)) - &
                    0.5*(TKE_prev(K)-TKE_prev(K-1))*Idz(k-1)* (dK(K-1) + dK(K)) + &
                    dz_Int(K) * (dK(K) * (S2(K) - N2(K)) - dQ(K)*TKE_decay(K))
      enddo ; endif

    endif  ! End of the Newton's method solver.

    ! Test kappa for convergence...
    if ((tol_err < Newton_err) .and. (.not.abort_Newton)) then
      !  A lower tolerance is used to switch to Newton's method than to switch back.
      Newton_test = Newton_err ; if (do_Newton) Newton_test = 2.0*Newton_err
      was_Newton = do_Newton
      within_tolerance = .true. ; do_Newton = .true.
      do K=min(ks_kappa,ks_kappa_prev),max(ke_kappa,ke_kappa_prev)
        kappa_mean = kappa0 + (kappa(K) - 0.5*dK(K))
        if (abs(dK(K)) > Newton_test * kappa_mean) then
          if (do_Newton) abort_Newton = .true.
          within_tolerance = .false. ; do_Newton = .false. ; exit
        elseif (abs(dK(K)) > tol_err * kappa_mean) then
          within_tolerance = .false. ; if (.not.do_Newton) exit
        endif
        if (abs(dQ(K)) > Newton_test*(tke(K) - 0.5*dQ(K))) then
          if (do_Newton) abort_Newton = .true.
          do_Newton = .false. ; if (.not.within_tolerance) exit
        endif
      enddo

    else  ! Newton's method will not be used again, so no need to check.
      within_tolerance = .true.
      do K=min(ks_kappa,ks_kappa_prev),max(ke_kappa,ke_kappa_prev)
        if (abs(dK(K)) > tol_err * (kappa0 + (kappa(K) - 0.5*dK(K)))) then
          within_tolerance = .false. ;  exit
        endif
      enddo
    endif

    if (abort_Newton) then
      do_Newton = .false. ; abort_Newton = .false.
      ! We went to Newton too quickly last time, so restrict the tolerance.
      Newton_err = 0.5*Newton_err
      ke_kappa_prev = nz
      do K=2,nz ; K_Q(K) = kappa(K) / max(TKE(K), TKE_min) ; enddo
    endif

    if (within_tolerance) exit

  enddo

  if (do_Newton) then  ! K_Q needs to be calculated.
    do K=1,ks_kappa-1 ;  K_Q(K) = 0.0 ; enddo
    do K=ks_kappa,ke_kappa ; K_Q(K) = kappa(K) / TKE(K) ; enddo
    do K=ke_kappa+1,nz+1 ; K_Q(K) = 0.0 ; enddo
  endif

  if (present(local_src)) then
    local_src(1) = 0.0 ; local_src(nz+1) = 0.0
    do K=2,nz
      diffusive_src = Idz(k-1)*(kappa(K-1)-kappa(K)) + Idz(k)*(kappa(K+1)-kappa(K))
      chg_by_k0 = kappa0 * ((Idz(k-1)+Idz(k)) / dz_Int(K) + I_Ld2(K))
      if (diffusive_src <= 0.0) then
        local_src(K) = K_src(K) + chg_by_k0
      else
        local_src(K) = (K_src(K) + chg_by_k0) + diffusive_src / dz_Int(K)
      endif
    enddo
  endif
  if (present(kappa_src)) then
    kappa_src(1) = 0.0 ; kappa_src(nz+1) = 0.0
    do K=2,nz
      kappa_src(K) = K_src(K)
    enddo
  endif

end subroutine find_kappa_tke

!> This subroutineinitializesthe parameters that regulate shear-driven mixing
function kappa_shear_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic
                                                 !! output.
  type(Kappa_shear_CS),    pointer       :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  logical :: kappa_shear_init !< True if module is to be used, False otherwise

  ! Local variables
  logical :: merge_mixedlayer
  logical :: just_read ! If true, this module is not used, so only read the parameters.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_kappa_shear"  ! This module's name.
  real :: kappa_0_unscaled  ! The value of kappa_0 in MKS units [m2 s-1]
  real :: KD_normal ! The KD of the main model, read here only as a parameter
                    ! for setting the default of KD_SMOOTH in MKS units [m2 s-1]

  if (associated(CS)) then
    call MOM_error(WARNING, "kappa_shear_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  !   The Jackson-Hallberg-Legg shear mixing parameterization uses the following
  ! 6 nondimensional coefficients.  That paper gives 3 best fit parameter sets.
  !    Ri_Crit  Rate    FRi_Curv  K_buoy  TKE_N  TKE_Shear
  ! p1: 0.25    0.089    -0.97     0.82    0.24    0.14
  ! p2: 0.30    0.085    -0.94     0.86    0.26    0.13
  ! p3: 0.35    0.088    -0.89     0.81    0.28    0.12
  !   Future research will reveal how these should be modified to take
  ! subgridscale inhomogeneity into account.

! Set default, read and log parameters
  call get_param(param_file, mdl, "USE_JACKSON_PARAM", kappa_shear_init, default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, &
    "Parameterization of shear-driven turbulence following Jackson, Hallberg and Legg, JPO 2008", &
    log_to_all=.true., debugging=kappa_shear_init, all_default=.not.kappa_shear_init)
  call get_param(param_file, mdl, "USE_JACKSON_PARAM", kappa_shear_init, &
                 "If true, use the Jackson-Hallberg-Legg (JPO 2008) "//&
                 "shear mixing parameterization.", default=.false.)
  just_read = .not.kappa_shear_init
  call get_param(param_file, mdl, "VERTEX_SHEAR", CS%KS_at_vertex, &
                 "If true, do the calculations of the shear-driven mixing "//&
                 "at the cell vertices (i.e., the vorticity points).", &
                 default=.false., do_not_log=just_read)
  call get_param(param_file, mdl, "RINO_CRIT", CS%RiNo_crit, &
                 "The critical Richardson number for shear mixing.", &
                 units="nondim", default=0.25, do_not_log=just_read)
  call get_param(param_file, mdl, "SHEARMIX_RATE", CS%Shearmix_rate, &
                 "A nondimensional rate scale for shear-driven entrainment. "//&
                 "Jackson et al find values in the range of 0.085-0.089.", &
                 units="nondim", default=0.089, do_not_log=just_read)
  call get_param(param_file, mdl, "MAX_RINO_IT", CS%max_RiNo_it, &
                 "The maximum number of iterations that may be used to "//&
                 "estimate the Richardson number driven mixing.", &
                 units="nondim", default=50, do_not_log=just_read)
  call get_param(param_file, mdl, "KD", KD_normal, default=0.0, do_not_log=.true.)
  call get_param(param_file, mdl, "KD_KAPPA_SHEAR_0", CS%kappa_0, &
                 "The background diffusivity that is used to smooth the "//&
                 "density and shear profiles before solving for the "//&
                 "diffusivities.  The default is the greater of KD and 1e-7 m2 s-1.", &
                 units="m2 s-1", default=max(KD_normal, 1.0e-7), scale=US%m2_s_to_Z2_T, &
                 unscaled=kappa_0_unscaled, do_not_log=just_read)
  call get_param(param_file, mdl, "KD_TRUNC_KAPPA_SHEAR", CS%kappa_trunc, &
                 "The value of shear-driven diffusivity that is considered negligible "//&
                 "and is rounded down to 0. The default is 1% of KD_KAPPA_SHEAR_0.", &
                 units="m2 s-1", default=0.01*kappa_0_unscaled, scale=US%m2_s_to_Z2_T, do_not_log=just_read)
  call get_param(param_file, mdl, "FRI_CURVATURE", CS%FRi_curvature, &
                 "The nondimensional curvature of the function of the "//&
                 "Richardson number in the kappa source term in the "//&
                 "Jackson et al. scheme.", units="nondim", default=-0.97, do_not_log=just_read)
  call get_param(param_file, mdl, "TKE_N_DECAY_CONST", CS%C_N, &
                 "The coefficient for the decay of TKE due to "//&
                 "stratification (i.e. proportional to N*tke). "//&
                 "The values found by Jackson et al. are 0.24-0.28.", &
                 units="nondim", default=0.24, do_not_log=just_read)
!  call get_param(param_file, mdl, "LAYER_KAPPA_STAGGER", CS%layer_stagger, &
!                 default=.false., do_not_log=just_read)
  call get_param(param_file, mdl, "TKE_SHEAR_DECAY_CONST", CS%C_S, &
                 "The coefficient for the decay of TKE due to shear (i.e. "//&
                 "proportional to |S|*tke). The values found by Jackson "//&
                 "et al. are 0.14-0.12.", units="nondim", default=0.14, do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_BUOY_SCALE_COEF", CS%lambda, &
                 "The coefficient for the buoyancy length scale in the "//&
                 "kappa equation.  The values found by Jackson et al. are "//&
                 "in the range of 0.81-0.86.", units="nondim", default=0.82, do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_N_OVER_S_SCALE_COEF2", CS%lambda2_N_S, &
                 "The square of the ratio of the coefficients of the "//&
                 "buoyancy and shear scales in the diffusivity equation, "//&
                 "Set this to 0 (the default) to eliminate the shear scale. "//&
                 "This is only used if USE_JACKSON_PARAM is true.", &
                 units="nondim", default=0.0, do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_SHEAR_TOL_ERR", CS%kappa_tol_err, &
                 "The fractional error in kappa that is tolerated. "//&
                 "Iteration stops when changes between subsequent "//&
                 "iterations are smaller than this everywhere in a "//&
                 "column.  The peak diffusivities usually converge most "//&
                 "rapidly, and have much smaller errors than this.", &
                 units="nondim", default=0.1, do_not_log=just_read)
  call get_param(param_file, mdl, "TKE_BACKGROUND", CS%TKE_bg, &
                 "A background level of TKE used in the first iteration "//&
                 "of the kappa equation.  TKE_BACKGROUND could be 0.", &
                 units="m2 s-2", default=0.0, scale=US%m_to_Z**2*US%T_to_s**2)
  call get_param(param_file, mdl, "KAPPA_SHEAR_ELIM_MASSLESS", CS%eliminate_massless, &
                 "If true, massless layers are merged with neighboring "//&
                 "massive layers in this calculation.  The default is "//&
                 "true and I can think of no good reason why it should "//&
                 "be false. This is only used if USE_JACKSON_PARAM is true.", &
                 default=.true., do_not_log=just_read)
  call get_param(param_file, mdl, "MAX_KAPPA_SHEAR_IT", CS%max_KS_it, &
                 "The maximum number of iterations that may be used to "//&
                 "estimate the time-averaged diffusivity.", units="nondim", &
                 default=13, do_not_log=just_read)
  call get_param(param_file, mdl, "PRANDTL_TURB", CS%Prandtl_turb, &
                 "The turbulent Prandtl number applied to shear instability.", &
                 units="nondim", default=1.0, do_not_log=.true.)
  call get_param(param_file, mdl, "VEL_UNDERFLOW", CS%vel_underflow, &
                 "A negligibly small velocity magnitude below which velocity components are set "//&
                 "to 0.  A reasonable value might be 1e-30 m/s, which is less than an "//&
                 "Angstrom divided by the age of the universe.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T, do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_SHEAR_MAX_KAP_SRC_CHG", CS%kappa_src_max_chg, &
                 "The maximum permitted increase in the kappa source within an iteration relative "//&
                 "to the local source; this must be greater than 1.  The lower limit for the "//&
                 "permitted fractional decrease is (1 - 0.5/kappa_src_max_chg).  These limits "//&
                 "could perhaps be made dynamic with an improved iterative solver.", &
                 default=10.0, units="nondim", do_not_log=just_read)

  call get_param(param_file, mdl, "DEBUG_KAPPA_SHEAR", CS%debug, &
                 "If true, write debugging data for the kappa-shear code. \n"//&
                 "Caution: this option is _very_ verbose and should only "//&
                 "be used in single-column mode!", &
                 default=.false., debuggingParam=.true., do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_SHEAR_ITER_BUG", CS%dKdQ_iteration_bug, &
                 "If true, use an older, dimensionally inconsistent estimate of the "//&
                 "derivative of diffusivity with energy in the Newton's method iteration.  "//&
                 "The bug causes undercorrections when dz > 1 m.", default=.false., do_not_log=just_read)
  call get_param(param_file, mdl, "KAPPA_SHEAR_ALL_LAYER_TKE_BUG", CS%all_layer_TKE_bug, &
                 "If true, report back the latest estimate of TKE instead of the time average "//&
                 "TKE when there is mass in all layers.  Otherwise always report the time "//&
                 "averaged TKE, as is currently done when there are some massless layers.", &
                 default=.false., do_not_log=just_read)
!    id_clock_KQ = cpu_clock_id('Ocean KS kappa_shear', grain=CLOCK_ROUTINE)
!    id_clock_avg = cpu_clock_id('Ocean KS avg', grain=CLOCK_ROUTINE)
!    id_clock_project = cpu_clock_id('Ocean KS project', grain=CLOCK_ROUTINE)
!    id_clock_setup = cpu_clock_id('Ocean KS setup', grain=CLOCK_ROUTINE)

  CS%nkml = 1
  if (GV%nkml>0) then
    call get_param(param_file, mdl, "KAPPA_SHEAR_MERGE_ML",merge_mixedlayer, &
                 "If true, combine the mixed layers together before solving the "//&
                 "kappa-shear equations.", default=.true., do_not_log=just_read)
    if (merge_mixedlayer) CS%nkml = GV%nkml
  endif

! Forego remainder of initialization if not using this scheme
  if (.not. kappa_shear_init) return

  CS%diag => diag

  CS%id_Kd_shear = register_diag_field('ocean_model','Kd_shear', diag%axesTi, Time, &
      'Shear-driven Diapycnal Diffusivity', 'm2 s-1', conversion=US%Z2_T_to_m2_s)
  CS%id_TKE = register_diag_field('ocean_model','TKE_shear', diag%axesTi, Time, &
      'Shear-driven Turbulent Kinetic Energy', 'm2 s-2', conversion=US%Z_to_m**2*US%s_to_T**2)

end function kappa_shear_init

!> This function indicates to other modules whether the Jackson et al shear mixing
!! parameterization will be used without needing to duplicate the log entry.
logical function kappa_shear_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
! Reads the parameter "USE_JACKSON_PARAM" and returns state.
  character(len=40)  :: mdl = "MOM_kappa_shear"  ! This module's name.

  call get_param(param_file, mdl, "USE_JACKSON_PARAM", kappa_shear_is_used, &
                 default=.false., do_not_log=.true.)
end function kappa_shear_is_used

!> This function indicates to other modules whether the Jackson et al shear mixing
!! parameterization will be used without needing to duplicate the log entry.
logical function kappa_shear_at_vertex(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
! Reads the parameter "USE_JACKSON_PARAM" and returns state.
  character(len=40)  :: mdl = "MOM_kappa_shear"  ! This module's name.

  logical :: do_kappa_shear

  call get_param(param_file, mdl, "USE_JACKSON_PARAM", do_kappa_shear, &
                 default=.false., do_not_log=.true.)
  kappa_shear_at_vertex = .false.
  if (do_Kappa_Shear) &
    call get_param(param_file, mdl, "VERTEX_SHEAR", kappa_shear_at_vertex, &
                 "If true, do the calculations of the shear-driven mixing "//&
                 "at the cell vertices (i.e., the vorticity points).", &
                 default=.false., do_not_log=.true.)

end function kappa_shear_at_vertex

!> \namespace mom_kappa_shear
!!
!! By Laura Jackson and Robert Hallberg, 2006-2008
!!
!!   This file contains the subroutines that determine the diapycnal
!! diffusivity driven by resolved shears, as specified by the
!! parameterizations described in Jackson and Hallberg (JPO, 2008).
!!
!!   The technique by which the 6 equations (for kappa, TKE, u, v, T,
!! and S) are solved simultaneously has been dramatically revised
!! from the previous version. The previous version was not converging
!! in some cases, especially near the surface mixed layer, while the
!! revised version does.  The revised version solves for kappa and
!! TKE with shear and stratification fixed, then marches the density
!! and velocities forward with an adaptive (and aggressive) time step
!! in a predictor-corrector-corrector emulation of a trapezoidal
!! scheme.  Run-time-settable parameters determine the tolerence to
!! which the kappa and TKE equations are solved and the minimum time
!! step that can be taken.

end module MOM_kappa_shear
