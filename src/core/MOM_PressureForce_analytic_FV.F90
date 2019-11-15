!> Analytically integrated finite volume pressure gradient
module MOM_PressureForce_AFV

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field
use MOM_diag_mediator, only : safe_alloc_ptr, diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_PressureForce_Mont, only : set_pbce_Bouss, set_pbce_nonBouss
use MOM_tidal_forcing, only : calc_tidal_forcing, tidal_forcing_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs
use MOM_EOS, only : int_density_dz, int_specific_vol_dp
use MOM_EOS, only : int_density_dz_generic_plm, int_density_dz_generic_ppm
use MOM_EOS, only : int_spec_vol_dp_generic_plm
use MOM_EOS, only : int_density_dz_generic, int_spec_vol_dp_generic
use MOM_ALE, only : pressure_gradient_plm, pressure_gradient_ppm, ALE_CS

implicit none ; private

#include <MOM_memory.h>

public PressureForce_AFV, PressureForce_AFV_init, PressureForce_AFV_end
public PressureForce_AFV_Bouss, PressureForce_AFV_nonBouss

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Finite volume pressure gradient control structure
type, public :: PressureForce_AFV_CS ; private
  logical :: tides          !< If true, apply tidal momentum forcing.
  real    :: Rho0           !< The density used in the Boussinesq
                            !! approximation [kg m-3].
  real    :: GFS_scale      !< A scaling of the surface pressure gradients to
                            !! allow the use of a reduced gravity model [nondim].
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                            !! timing of diagnostic output.
  logical :: useMassWghtInterp !< Use mass weighting in T/S interpolation
  logical :: boundary_extrap !< Indicate whether high-order boundary
                            !! extrapolation should be used within boundary cells

  logical :: reconstruct    !< If true, polynomial profiles of T & S will be
                            !! reconstructed and used in the integrals for the
                            !! finite volume pressure gradient calculation.
                            !! The default depends on whether regridding is being used.

  integer :: Recon_Scheme   !< Order of the polynomial of the reconstruction of T & S
                            !! for the finite volume pressure gradient calculation.
                            !! By the default (1) is for a piecewise linear method

  integer :: id_e_tidal = -1 !< Diagnostic identifier
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL() !< Tides control structure
end type PressureForce_AFV_CS

contains

!> Thin interface between the model and the Boussinesq and non-Boussinesq
!! pressure force routines.
subroutine PressureForce_AFV(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, p_atm, pbce, eta)
  type(ocean_grid_type),                     intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !< Vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(inout) :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out)   :: PFu !< Zonal acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out)   :: PFv !< Meridional acceleration [L T-2 ~> m s-2]
  type(PressureForce_AFV_CS),                pointer       :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                              pointer       :: ALE_CSp !< ALE control structure
  real, dimension(:,:),                      optional, pointer :: p_atm !< The pressure at the ice-ocean
                                                           !! or atmosphere-ocean interface [Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  optional, intent(out) :: pbce !< The baroclinic pressure
                                                           !! anomaly in each layer due to eta anomalies
                                                           !! [m2 s-2 H-1 ~> m s-2 or m4 s-2 kg-1].
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< The bottom mass used to
                                                           !! calculate PFu and PFv [H ~> m or kg m-2], with any tidal
                                                           !! contributions or compressibility compensation.

  if (GV%Boussinesq) then
    call PressureForce_AFV_bouss(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, p_atm, pbce, eta)
  else
    call PressureForce_AFV_nonbouss(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, p_atm, pbce, eta)
  endif

end subroutine PressureForce_AFV

!> \brief Non-Boussinesq analytically-integrated finite volume form of pressure gradient
!!
!! Determines the acceleration due to hydrostatic pressure forces, using
!! the analytic finite volume form of the Pressure gradient, and does not
!! make the Boussinesq approximation.
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_AFV_nonBouss(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, p_atm, pbce, eta)
  type(ocean_grid_type),                     intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),                     intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)  :: h   !< Layer thickness [H ~> kg/m2]
  type(thermo_var_ptrs),                     intent(in)  :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out) :: PFu !< Zonal acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out) :: PFv !< Meridional acceleration [L T-2 ~> m s-2]
  type(PressureForce_AFV_CS),                pointer     :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                              pointer     :: ALE_CSp !< ALE control structure
  real, dimension(:,:),                      optional, pointer :: p_atm !< The pressure at the ice-ocean
                                                           !! or atmosphere-ocean interface [Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  optional, intent(out) :: pbce !< The baroclinic pressure
                                                           !! anomaly in each layer due to eta anomalies
                                                           !! [m2 s-2 H-1 ~> m s-2 or m4 s-2 kg-1].
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< The bottom mass used to
                                                           !! calculate PFu and PFv [H ~> m or kg m-2], with any tidal
                                                           !! contributions or compressibility compensation.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: p ! Interface pressure [Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [degC].
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [ppt].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    S_t, &      ! Top and bottom edge values for linear reconstructions
    S_b, &      ! of salinity within each layer [ppt].
    T_t, &      ! Top and bottom edge values for linear reconstructions
    T_b         ! of temperature within each layer [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G))  :: &
    dza, &      ! The change in geopotential anomaly between the top and bottom
                ! of a layer [m2 s-2].
    intp_dza    ! The vertical integral in depth of the pressure anomaly less
                ! the pressure anomaly at the top of the layer [Pa m2 s-2].
  real, dimension(SZI_(G),SZJ_(G))  :: &
    dp, &       ! The (positive) change in pressure across a layer [Pa].
    SSH, &      ! The sea surface height anomaly, in depth units [Z ~> m].
    e_tidal, &  ! The bottom geopotential anomaly due to tidal forces from
                ! astronomical sources and self-attraction and loading [Z ~> m].
    dM, &       ! The barotropic adjustment to the Montgomery potential to
                ! account for a reduced gravity model [L2 T-2 ~> m2 s-2].
    za          ! The geopotential anomaly (i.e. g*e + alpha_0*pressure) at the
                ! interface atop a layer [m2 s-2].

  real, dimension(SZI_(G)) :: Rho_cv_BL !  The coordinate potential density in the deepest variable
                ! density near-surface layer [R ~> kg m-3].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    intx_za     ! The zonal integral of the geopotential anomaly along the
                ! interface below a layer, divided by the grid spacing [m2 s-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: &
    intx_dza    ! The change in intx_za through a layer [m2 s-2].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    inty_za     ! The meridional integral of the geopotential anomaly along the
                ! interface below a layer, divided by the grid spacing [m2 s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: &
    inty_dza    ! The change in inty_za through a layer [m2 s-2].
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density, [Pa] (usually 2e7 Pa = 2000 dbar).

  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [Pa].
  real :: g_Earth_z          ! A scaled version of g_Earth [m2 Z-1 s-2 ~> m s-2].
  real :: I_gEarth           ! The inverse of g_Earth_z [s2 Z m-2 ~> s2 m-1]
  real :: alpha_anom         ! The in-situ specific volume, averaged over a
                             ! layer, less alpha_ref [m3 kg-1].
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_ALE         ! If true, use an ALE pressure reconstruction.
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.

  real :: alpha_ref ! A reference specific volume [m3 kg-1], that is used
                    ! to reduce the impact of truncation errors.
  real :: rho_in_situ(SZI_(G)) ! The in situ density [kg m-3].
  real :: Pa_to_H   ! A factor to convert from Pa to the thicknesss units (H).
!  real :: oneatm = 101325.0  ! 1 atm in [Pa] = [kg m-1 s-2]
  real, parameter :: C1_6 = 1.0/6.0
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: i, j, k

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  nkmb=GV%nk_rho_varies
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.associated(CS)) call MOM_error(FATAL, &
       "MOM_PressureForce_AFV_nonBouss: Module must be initialized before it is used.")

  use_p_atm = .false.
  if (present(p_atm)) then ; if (associated(p_atm)) use_p_atm = .true. ; endif
  use_EOS = associated(tv%eqn_of_state)
  use_ALE = .false.
  if (associated(ALE_CSp)) use_ALE = CS%reconstruct .and. use_EOS

  dp_neglect = GV%H_to_Pa * GV%H_subroundoff
  alpha_ref = 1.0/CS%Rho0
  g_Earth_z = US%L_T_to_m_s**2 * GV%g_Earth
  I_gEarth = 1.0 / g_Earth_z

  if (use_p_atm) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      p(i,j,1) = p_atm(i,j)
    enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      p(i,j,1) = 0.0 ! or oneatm
    enddo ; enddo
  endif
  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1 ; do k=2,nz+1 ; do i=Isq,Ieq+1
    p(i,j,K) = p(i,j,K-1) + GV%H_to_Pa * h(i,j,k-1)
  enddo ; enddo ; enddo

  if (use_EOS) then
  !   With a bulk mixed layer, replace the T & S of any layers that are
  ! lighter than the the buffer layer with the properties of the buffer
  ! layer.  These layers will be massless anyway, and it avoids any
  ! formal calculations with hydrostatically unstable profiles.
    if (nkmb>0) then
      tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
      tv_tmp%eqn_of_state => tv%eqn_of_state
      do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
      !$OMP parallel do default(shared) private(Rho_cv_BL)
      do j=Jsq,Jeq+1
        do k=1,nkmb ; do i=Isq,Ieq+1
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        enddo ; enddo
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, &
                        Rho_cv_BL(:), Isq, Ieq-Isq+2, tv%eqn_of_state, scale=US%kg_m3_to_R)
        do k=nkmb+1,nz ; do i=Isq,Ieq+1
          if (GV%Rlay(k) < Rho_cv_BL(i)) then
            tv_tmp%T(i,j,k) = tv%T(i,j,nkmb) ; tv_tmp%S(i,j,k) = tv%S(i,j,nkmb)
          else
            tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
          endif
        enddo ; enddo
      enddo
    else
      tv_tmp%T => tv%T ; tv_tmp%S => tv%S
      tv_tmp%eqn_of_state => tv%eqn_of_state
    endif
  endif

  ! If regridding is activated, do a linear reconstruction of salinity
  ! and temperature across each layer. The subscripts 't' and 'b' refer
  ! to top and bottom values within each layer (these are the only degrees
  ! of freedeom needed to know the linear profile).
  if ( use_ALE ) then
    if ( CS%Recon_Scheme == 1 ) then
      call pressure_gradient_plm(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
    elseif ( CS%Recon_Scheme == 2) then
      call pressure_gradient_ppm(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
    endif
  endif

  !$OMP parallel do default(shared) private(alpha_anom,dp)
  do k=1,nz
    ! Calculate 4 integrals through the layer that are required in the
    ! subsequent calculation.
    if (use_EOS) then
      if ( use_ALE ) then
        if ( CS%Recon_Scheme == 1 ) then
          call int_spec_vol_dp_generic_plm( T_t(:,:,k), T_b(:,:,k), &
                    S_t(:,:,k), S_b(:,:,k), p(:,:,K), p(:,:,K+1), &
                    alpha_ref, dp_neglect, p(:,:,nz+1), G%HI, &
                    tv%eqn_of_state, dza(:,:,k), intp_dza(:,:,k), &
                    intx_dza(:,:,k), inty_dza(:,:,k), &
                    useMassWghtInterp = CS%useMassWghtInterp)
          i=k
        elseif ( CS%Recon_Scheme == 2 ) then
          call MOM_error(FATAL, "PressureForce_AFV_nonBouss: "//&
                         "int_spec_vol_dp_generic_ppm does not exist yet.")
        !  call int_spec_vol_dp_generic_ppm ( tv%T(:,:,k), T_t(:,:,k), T_b(:,:,k), &
        !            tv%S(:,:,k), S_t(:,:,k), S_b(:,:,k), p(:,:,K), p(:,:,K+1), &
        !            alpha_ref, G%HI, tv%eqn_of_state, dza(:,:,k), intp_dza(:,:,k), &
        !            intx_dza(:,:,k), inty_dza(:,:,k))
        endif
      else
        call int_specific_vol_dp(tv_tmp%T(:,:,k), tv_tmp%S(:,:,k), p(:,:,K), &
                               p(:,:,K+1), alpha_ref, G%HI, tv%eqn_of_state, &
                               dza(:,:,k), intp_dza(:,:,k), intx_dza(:,:,k), &
                               inty_dza(:,:,k), bathyP=p(:,:,nz+1), dP_tiny=dp_neglect, &
                               useMassWghtInterp = CS%useMassWghtInterp)
      endif
    else
      alpha_anom = 1.0/(US%R_to_kg_m3*GV%Rlay(k)) - alpha_ref
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dp(i,j) = GV%H_to_Pa * h(i,j,k)
        dza(i,j,k) = alpha_anom * dp(i,j)
        intp_dza(i,j,k) = 0.5 * alpha_anom * dp(i,j)**2
      enddo ; enddo
      do j=js,je ; do I=Isq,Ieq
        intx_dza(i,j,k) = 0.5 * alpha_anom * (dp(i,j)+dp(i+1,j))
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        inty_dza(i,j,k) = 0.5 * alpha_anom * (dp(i,j)+dp(i,j+1))
      enddo ; enddo
    endif
  enddo

  !   The bottom geopotential anomaly is calculated first so that the increments
  ! to the geopotential anomalies can be reused.  Alternately, the surface
  ! geopotential could be calculated directly with separate calls to
  ! int_specific_vol_dp with alpha_ref=0, and the anomalies used going
  ! downward, which would relieve the need for dza, intp_dza, intx_dza, and
  ! inty_dza to be 3-D arrays.

  ! Sum vertically to determine the surface geopotential anomaly.
  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1
    do i=Isq,Ieq+1
      za(i,j) = alpha_ref*p(i,j,nz+1) - g_Earth_z*G%bathyT(i,j)
    enddo
    do k=nz,1,-1 ; do i=Isq,Ieq+1
      za(i,j) = za(i,j) + dza(i,j,k)
    enddo ; enddo
  enddo

  if (CS%tides) then
    ! Find and add the tidal geopotential anomaly.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      SSH(i,j) = (za(i,j) - alpha_ref*p(i,j,1)) * I_gEarth
    enddo ; enddo
    call calc_tidal_forcing(CS%Time, SSH, e_tidal, G, CS%tides_CSp, m_to_Z=US%m_to_Z)
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      za(i,j) = za(i,j) - g_Earth_z * e_tidal(i,j)
    enddo ; enddo
  endif

  if (CS%GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    if (use_EOS) then
      !$OMP parallel do default(shared) private(rho_in_situ)
      do j=Jsq,Jeq+1
        call calculate_density(tv_tmp%T(:,j,1), tv_tmp%S(:,j,1), p(:,j,1), &
                               rho_in_situ, Isq, Ieq-Isq+2, tv%eqn_of_state)

        do i=Isq,Ieq+1
          dM(i,j) = (CS%GFS_scale - 1.0) * US%m_s_to_L_T**2 * &
            (p(i,j,1)*(1.0/rho_in_situ(i) - alpha_ref) + za(i,j))
        enddo
      enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dM(i,j) = (CS%GFS_scale - 1.0) * US%m_s_to_L_T**2 * &
          (p(i,j,1)*(1.0/(US%R_to_kg_m3*GV%Rlay(1)) - alpha_ref) + za(i,j))
      enddo ; enddo
    endif
!  else
!    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1 ; dM(i,j) = 0.0 ; enddo ; enddo
  endif

  !   This order of integrating upward and then downward again is necessary with
  ! a nonlinear equation of state, so that the surface geopotentials will go
  ! linearly between the values at thickness points, but the bottom
  ! geopotentials will not now be linear at the sub-grid-scale.  Doing this
  ! ensures no motion with flat isopycnals, even with a nonlinear equation of state.
  !$OMP parallel do default(shared)
  do j=js,je ; do I=Isq,Ieq
    intx_za(I,j) = 0.5*(za(i,j) + za(i+1,j))
  enddo ; enddo
  !$OMP parallel do default(shared)
  do J=Jsq,Jeq ; do i=is,ie
    inty_za(i,J) = 0.5*(za(i,j) + za(i,j+1))
  enddo ; enddo
  do k=1,nz
    ! These expressions for the acceleration have been carefully checked in
    ! a set of idealized cases, and should be bug-free.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      dp(i,j) = GV%H_to_Pa*h(i,j,k)
      za(i,j) = za(i,j) - dza(i,j,k)
    enddo ; enddo
    !$OMP parallel do default(shared)
    do j=js,je ; do I=Isq,Ieq
      intx_za(I,j) = intx_za(I,j) - intx_dza(I,j,k)
      PFu(I,j,k) = (((za(i,j)*dp(i,j) + intp_dza(i,j,k)) - &
                   (za(i+1,j)*dp(i+1,j) + intp_dza(i+1,j,k))) + &
                   ((dp(i+1,j) - dp(i,j)) * intx_za(I,j) - &
                   (p(i+1,j,K) - p(i,j,K)) * intx_dza(I,j,k))) * &
                   (US%m_s_to_L_T**2 * 2.0*G%IdxCu(I,j) / &
                   ((dp(i,j) + dp(i+1,j)) + dp_neglect))
    enddo ; enddo
    !$OMP parallel do default(shared)
    do J=Jsq,Jeq ; do i=is,ie
      inty_za(i,J) = inty_za(i,J) - inty_dza(i,J,k)
      PFv(i,J,k) = (((za(i,j)*dp(i,j) + intp_dza(i,j,k)) - &
                   (za(i,j+1)*dp(i,j+1) + intp_dza(i,j+1,k))) + &
                   ((dp(i,j+1) - dp(i,j)) * inty_za(i,J) - &
                   (p(i,j+1,K) - p(i,j,K)) * inty_dza(i,J,k))) * &
                   (US%m_s_to_L_T**2 * 2.0*G%IdyCv(i,J) / &
                   ((dp(i,j) + dp(i,j+1)) + dp_neglect))
    enddo ; enddo

    if (CS%GFS_scale < 1.0) then
      ! Adjust the Montgomery potential to make this a reduced gravity model.
      !$OMP parallel do default(shared)
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = PFu(I,j,k) - (dM(i+1,j) - dM(i,j)) * G%IdxCu(I,j)
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = PFv(i,J,k) - (dM(i,j+1) - dM(i,j)) * G%IdyCv(i,J)
      enddo ; enddo
    endif
  enddo

  if (present(pbce)) then
    call set_pbce_nonBouss(p, tv_tmp, G, GV, US, CS%GFS_scale, pbce)
  endif

  if (present(eta)) then
    Pa_to_H = 1.0 / GV%H_to_Pa
    if (use_p_atm) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = (p(i,j,nz+1) - p_atm(i,j))*Pa_to_H ! eta has the same units as h.
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = p(i,j,nz+1)*Pa_to_H ! eta has the same units as h.
      enddo ; enddo
    endif
  endif

  if (CS%id_e_tidal>0) call post_data(CS%id_e_tidal, e_tidal, CS%diag)

end subroutine PressureForce_AFV_nonBouss

!> \brief Boussinesq analytically-integrated finite volume form of pressure gradient
!!
!! Determines the acceleration due to hydrostatic pressure forces, using
!! the finite volume form of the terms and analytic integrals in depth.
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_AFV_Bouss(h, tv, PFu, PFv, G, GV, US, CS, ALE_CSp, p_atm, pbce, eta)
  type(ocean_grid_type),                     intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV  !< Vertical grid structure
  type(unit_scale_type),                     intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)  :: h   !< Layer thickness [H ~> m]
  type(thermo_var_ptrs),                     intent(in)  :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out) :: PFu !< Zonal acceleration [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out) :: PFv !< Meridional acceleration [L T-2 ~> m s-2]
  type(PressureForce_AFV_CS),                pointer     :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                              pointer     :: ALE_CSp !< ALE control structure
  real, dimension(:,:),                      optional, pointer :: p_atm !< The pressure at the ice-ocean
                                                         !! or atmosphere-ocean interface [Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  optional, intent(out) :: pbce !< The baroclinic pressure
                                                         !! anomaly in each layer due to eta anomalies
                                                         !! [m2 s-2 H-1 ~> m s-2 or m4 s-2 kg-1].
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< The bottom mass used to
                                                         !! calculate PFu and PFv [H ~> m or kg m-2], with any
                                                         !! tidal contributions or compressibility compensation.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: e ! Interface height in depth units [Z ~> m].
  real, dimension(SZI_(G),SZJ_(G))  :: &
    e_tidal, &  ! The bottom geopotential anomaly due to tidal forces from
                ! astronomical sources and self-attraction and loading [Z ~> m].
    dM          ! The barotropic adjustment to the Montgomery potential to
                ! account for a reduced gravity model [L2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G)) :: &
    Rho_cv_BL   !   The coordinate potential density in the deepest variable
                ! density near-surface layer [R ~> kg m-3].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    dz_geo, &   ! The change in geopotential thickness through a layer times some dimensional
                ! rescaling factors [kg m-1 R-1 s-2 ~> m2 s-2].
    pa, &       ! The pressure anomaly (i.e. pressure + g*RHO_0*e) at the
                ! the interface atop a layer [Pa].
    dpa, &      ! The change in pressure anomaly between the top and bottom
                ! of a layer [Pa].
    intz_dpa    ! The vertical integral in depth of the pressure anomaly less the
                ! pressure anomaly at the top of the layer [H Pa ~> m Pa or kg m-2 Pa].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    intx_pa, &  ! The zonal integral of the pressure anomaly along the interface
                ! atop a layer, divided by the grid spacing [Pa].
    intx_dpa    ! The change in intx_pa through a layer [Pa].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    inty_pa, &  ! The meridional integral of the pressure anomaly along the
                ! interface atop a layer, divided by the grid spacing [Pa].
    inty_dpa    ! The change in inty_pa through a layer [Pa].

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [degC].
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [ppt].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    S_t, S_b, T_t, T_b ! Top and bottom edge values for linear reconstructions
                       ! of salinity and temperature within each layer.
  real :: rho_in_situ(SZI_(G)) ! The in situ density [R ~> kg m-3].
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density, [Pa] (usually 2e7 Pa = 2000 dbar).
  real :: p0(SZI_(G))        ! An array of zeros to use for pressure [Pa].
  real :: h_neglect          ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [H ~> m].
  real :: g_Earth_mks_z      ! A scaled version of g_Earth [m2 Z-1 s-2 ~> m s-2].
  real :: g_Earth_z_geo      ! Another scaled version of g_Earth [R m5 kg-1 Z-1 s-2 ~> m s-2].
  real :: I_Rho0             ! 1/Rho0 times unit scaling factors [L2 m kg-1 s2 T-2 ~> m3 kg-1].
  real :: G_Rho0             ! G_Earth / Rho0 in [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1].
  real :: Rho_ref            ! The reference density [R ~> kg m-3].
  real :: Rho_ref_mks        ! The reference density in mks units [kg m-3].
  real :: dz_neglect         ! A minimal thickness [Z ~> m], like e.
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_ALE         ! If true, use an ALE pressure reconstruction.
  logical :: use_EOS         ! If true, density is calculated from T & S using an equation of state.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.

  real, parameter :: C1_6 = 1.0/6.0
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: i, j, k

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  nkmb=GV%nk_rho_varies
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.associated(CS)) call MOM_error(FATAL, &
       "MOM_PressureForce_AFV_Bouss: Module must be initialized before it is used.")

  use_p_atm = .false.
  if (present(p_atm)) then ; if (associated(p_atm)) use_p_atm = .true. ; endif
  use_EOS = associated(tv%eqn_of_state)
  do i=Isq,Ieq+1 ; p0(i) = 0.0 ; enddo
  use_ALE = .false.
  if (associated(ALE_CSp)) use_ALE = CS%reconstruct .and. use_EOS

  h_neglect = GV%H_subroundoff
  dz_neglect = GV%H_subroundoff * GV%H_to_Z
  I_Rho0 = US%m_s_to_L_T**2 / (US%R_to_kg_m3*GV%Rho0)
  g_Earth_mks_z = US%L_T_to_m_s**2 * GV%g_Earth
  g_Earth_z_geo = US%R_to_kg_m3*US%L_T_to_m_s**2 * GV%g_Earth
  G_Rho0 = GV%g_Earth / GV%Rho0
  rho_ref_mks = CS%Rho0
  rho_ref = rho_ref_mks*US%kg_m3_to_R

  if (CS%tides) then
    !   Determine the surface height anomaly for calculating self attraction
    ! and loading.  This should really be based on bottom pressure anomalies,
    ! but that is not yet implemented, and the current form is correct for
    ! barotropic tides.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        e(i,j,1) = -G%bathyT(i,j)
      enddo
      do k=1,nz ; do i=Isq,Ieq+1
        e(i,j,1) = e(i,j,1) + h(i,j,k)*GV%H_to_Z
      enddo ; enddo
    enddo
    call calc_tidal_forcing(CS%Time, e(:,:,1), e_tidal, G, CS%tides_CSp, m_to_Z=US%m_to_Z)
  endif

!    Here layer interface heights, e, are calculated.
  if (CS%tides) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,nz+1) = -(G%bathyT(i,j) + e_tidal(i,j))
    enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,nz+1) = -G%bathyT(i,j)
    enddo ; enddo
  endif
  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1; do k=nz,1,-1 ; do i=Isq,Ieq+1
    e(i,j,K) = e(i,j,K+1) + h(i,j,k)*GV%H_to_Z
  enddo ; enddo ; enddo

  if (use_EOS) then
! With a bulk mixed layer, replace the T & S of any layers that are
! lighter than the the buffer layer with the properties of the buffer
! layer.  These layers will be massless anyway, and it avoids any
! formal calculations with hydrostatically unstable profiles.

    if (nkmb>0) then
      tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
      tv_tmp%eqn_of_state => tv%eqn_of_state

      do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
     !$OMP parallel do default(shared) private(Rho_cv_BL)
      do j=Jsq,Jeq+1
        do k=1,nkmb ; do i=Isq,Ieq+1
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        enddo ; enddo
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, &
                        Rho_cv_BL(:), Isq, Ieq-Isq+2, tv%eqn_of_state, scale=US%kg_m3_to_R)

        do k=nkmb+1,nz ; do i=Isq,Ieq+1
          if (GV%Rlay(k) < Rho_cv_BL(i)) then
            tv_tmp%T(i,j,k) = tv%T(i,j,nkmb) ; tv_tmp%S(i,j,k) = tv%S(i,j,nkmb)
          else
            tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
          endif
        enddo ; enddo
      enddo
    else
      tv_tmp%T => tv%T ; tv_tmp%S => tv%S
      tv_tmp%eqn_of_state => tv%eqn_of_state
    endif
  endif

  if (CS%GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    if (use_EOS) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1
        if (use_p_atm) then
          call calculate_density(tv_tmp%T(:,j,1), tv_tmp%S(:,j,1), p_atm(:,j), rho_in_situ, &
                                 Isq, Ieq-Isq+2, tv%eqn_of_state, scale=US%kg_m3_to_R)
        else
          call calculate_density(tv_tmp%T(:,j,1), tv_tmp%S(:,j,1), p0, rho_in_situ, &
                                 Isq, Ieq-Isq+2, tv%eqn_of_state, scale=US%kg_m3_to_R)
        endif
        do i=Isq,Ieq+1
          dM(i,j) = (CS%GFS_scale - 1.0) * (G_Rho0 * rho_in_situ(i)) * e(i,j,1)
        enddo
      enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dM(i,j) = (CS%GFS_scale - 1.0) * (G_Rho0 * GV%Rlay(1)) * e(i,j,1)
      enddo ; enddo
    endif
  endif
  ! I have checked that rho_0 drops out and that the 1-layer case is right. RWH.

  ! If regridding is activated, do a linear reconstruction of salinity
  ! and temperature across each layer. The subscripts 't' and 'b' refer
  ! to top and bottom values within each layer (these are the only degrees
  ! of freedeom needed to know the linear profile).
  if ( use_ALE ) then
    if ( CS%Recon_Scheme == 1 ) then
      call pressure_gradient_plm(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
    elseif ( CS%Recon_Scheme == 2 ) then
      call pressure_gradient_ppm(ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h, CS%boundary_extrap)
    endif
  endif

  ! Set the surface boundary conditions on pressure anomaly and its horizontal
  ! integrals, assuming that the surface pressure anomaly varies linearly
  ! in x and y.
  if (use_p_atm) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pa(i,j) = (rho_ref*g_Earth_z_geo)*e(i,j,1) + p_atm(i,j)
    enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pa(i,j) = (rho_ref*g_Earth_z_geo)*e(i,j,1)
    enddo ; enddo
  endif
  !$OMP parallel do default(shared)
  do j=js,je ; do I=Isq,Ieq
    intx_pa(I,j) = 0.5*(pa(i,j) + pa(i+1,j))
  enddo ; enddo
  !$OMP parallel do default(shared)
  do J=Jsq,Jeq ; do i=is,ie
    inty_pa(i,J) = 0.5*(pa(i,j) + pa(i,j+1))
  enddo ; enddo

  do k=1,nz
    ! Calculate 4 integrals through the layer that are required in the
    ! subsequent calculation.

    if (use_EOS) then
      ! The following routine computes the integrals that are needed to
      ! calculate the pressure gradient force. Linear profiles for T and S are
      ! assumed when regridding is activated. Otherwise, the previous version
      ! is used, whereby densities within each layer are constant no matter
      ! where the layers are located.
      if ( use_ALE ) then
        if ( CS%Recon_Scheme == 1 ) then
          call int_density_dz_generic_plm( T_t(:,:,k), T_b(:,:,k), &
                    S_t(:,:,k), S_b(:,:,k), e(:,:,K), e(:,:,K+1), &
                    rho_ref_mks, CS%Rho0, g_Earth_mks_z, &
                    dz_neglect, G%bathyT, G%HI, G%HI, &
                    tv%eqn_of_state, dpa, intz_dpa, intx_dpa, inty_dpa, &
                    useMassWghtInterp = CS%useMassWghtInterp)
        elseif ( CS%Recon_Scheme == 2 ) then
          call int_density_dz_generic_ppm( tv%T(:,:,k), T_t(:,:,k), T_b(:,:,k), &
                    tv%S(:,:,k), S_t(:,:,k), S_b(:,:,k), e(:,:,K), e(:,:,K+1), &
                    rho_ref_mks, CS%Rho0, g_Earth_mks_z, &
                    G%HI, G%HI, tv%eqn_of_state, dpa, intz_dpa, &
                    intx_dpa, inty_dpa)
        endif
      else
        call int_density_dz(tv_tmp%T(:,:,k), tv_tmp%S(:,:,k), e(:,:,K), e(:,:,K+1), &
                  rho_ref_mks, CS%Rho0, g_Earth_mks_z, G%HI, G%HI, tv%eqn_of_state, &
                  dpa, intz_dpa, intx_dpa, inty_dpa, &
                  G%bathyT, dz_neglect, CS%useMassWghtInterp)
      endif
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        intz_dpa(i,j) = intz_dpa(i,j)*GV%Z_to_H
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dz_geo(i,j) = g_Earth_z_geo * GV%H_to_Z*h(i,j,k)
        dpa(i,j) = (GV%Rlay(k) - rho_ref) * dz_geo(i,j)
        intz_dpa(i,j) = 0.5*(GV%Rlay(k) - rho_ref) * dz_geo(i,j)*h(i,j,k)
      enddo ; enddo
      !$OMP parallel do default(shared)
      do j=js,je ; do I=Isq,Ieq
        intx_dpa(I,j) = 0.5*(GV%Rlay(k) - rho_ref) * (dz_geo(i,j) + dz_geo(i+1,j))
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=Jsq,Jeq ; do i=is,ie
        inty_dpa(i,J) = 0.5*(GV%Rlay(k) - rho_ref) * (dz_geo(i,j) + dz_geo(i,j+1))
      enddo ; enddo
    endif

    ! Compute pressure gradient in x direction
    !$OMP parallel do default(shared)
    do j=js,je ; do I=Isq,Ieq
      PFu(I,j,k) = (((pa(i,j)*h(i,j,k) + intz_dpa(i,j)) - &
                   (pa(i+1,j)*h(i+1,j,k) + intz_dpa(i+1,j))) + &
                   ((h(i+1,j,k) - h(i,j,k)) * intx_pa(I,j) - &
                   (e(i+1,j,K+1) - e(i,j,K+1)) * intx_dpa(I,j) * GV%Z_to_H)) * &
                   ((2.0*I_Rho0*G%IdxCu(I,j)) / &
                   ((h(i,j,k) + h(i+1,j,k)) + h_neglect))
      intx_pa(I,j) = intx_pa(I,j) + intx_dpa(I,j)
    enddo ; enddo
    ! Compute pressure gradient in y direction
    !$OMP parallel do default(shared)
    do J=Jsq,Jeq ; do i=is,ie
      PFv(i,J,k) = (((pa(i,j)*h(i,j,k) + intz_dpa(i,j)) - &
                   (pa(i,j+1)*h(i,j+1,k) + intz_dpa(i,j+1))) + &
                   ((h(i,j+1,k) - h(i,j,k)) * inty_pa(i,J) - &
                   (e(i,j+1,K+1) - e(i,j,K+1)) * inty_dpa(i,J) * GV%Z_to_H)) * &
                   ((2.0*I_Rho0*G%IdyCv(i,J)) / &
                   ((h(i,j,k) + h(i,j+1,k)) + h_neglect))
      inty_pa(i,J) = inty_pa(i,J) + inty_dpa(i,J)
    enddo ; enddo
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pa(i,j) = pa(i,j) + dpa(i,j)
    enddo ; enddo
  enddo

  if (CS%GFS_scale < 1.0) then
    do k=1,nz
      !$OMP parallel do default(shared)
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = PFu(I,j,k) - (dM(i+1,j) - dM(i,j)) * G%IdxCu(I,j)
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = PFv(i,J,k) - (dM(i,j+1) - dM(i,j)) * G%IdyCv(i,J)
      enddo ; enddo
    enddo
  endif

  if (present(pbce)) then
    call set_pbce_Bouss(e, tv_tmp, G, GV, US, CS%Rho0, CS%GFS_scale, pbce)
  endif

  if (present(eta)) then
    if (CS%tides) then
    ! eta is the sea surface height relative to a time-invariant geoid, for
    ! comparison with what is used for eta in btstep.  See how e was calculated
    ! about 200 lines above.
    !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = e(i,j,1)*GV%Z_to_H + e_tidal(i,j)*GV%Z_to_H
      enddo ; enddo
    else
    !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = e(i,j,1)*GV%Z_to_H
      enddo ; enddo
    endif
  endif

  if (CS%id_e_tidal>0) call post_data(CS%id_e_tidal, e_tidal, CS%diag)

end subroutine PressureForce_AFV_Bouss

!> Initializes the finite volume pressure gradient control structure
subroutine PressureForce_AFV_init(Time, G, GV, US, param_file, diag, CS, tides_CSp)
  type(time_type), target,    intent(in)    :: Time !< Current model time
  type(ocean_grid_type),      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),      intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target,    intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_AFV_CS), pointer       :: CS !< Finite volume PGF control structure
  type(tidal_forcing_CS), optional, pointer :: tides_CSp !< Tides control structure
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl  ! This module's name.
  logical :: use_ALE

  if (associated(CS)) then
    call MOM_error(WARNING, "PressureForce_init called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag ; CS%Time => Time
  if (present(tides_CSp)) then
    if (associated(tides_CSp)) CS%tides_CSp => tides_CSp
  endif

  mdl = "MOM_PressureForce_AFV"
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mdl, "TIDES", CS%tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  call get_param(param_file, "MOM", "USE_REGRIDDING", use_ALE, &
                 "If True, use the ALE algorithm (regridding/remapping). "//&
                 "If False, use the layered isopycnal algorithm.", default=.false. )
  call get_param(param_file, mdl, "MASS_WEIGHT_IN_PRESSURE_GRADIENT", CS%useMassWghtInterp, &
                 "If true, use mass weighting when interpolating T/S for "//&
                 "integrals near the bathymetry in AFV pressure gradient "//&
                 "calculations.", default=.false.)
  call get_param(param_file, mdl, "RECONSTRUCT_FOR_PRESSURE", CS%reconstruct, &
                 "If True, use vertical reconstruction of T & S within "//&
                 "the integrals of the FV pressure gradient calculation. "//&
                 "If False, use the constant-by-layer algorithm. "//&
                 "The default is set by USE_REGRIDDING.", &
                 default=use_ALE )
  call get_param(param_file, mdl, "PRESSURE_RECONSTRUCTION_SCHEME", CS%Recon_Scheme, &
                 "Order of vertical reconstruction of T/S to use in the "//&
                 "integrals within the FV pressure gradient calculation.\n"//&
                 " 0: PCM or no reconstruction.\n"//&
                 " 1: PLM reconstruction.\n"//&
                 " 2: PPM reconstruction.", default=1)
  call get_param(param_file, mdl, "BOUNDARY_EXTRAPOLATION_PRESSURE", CS%boundary_extrap, &
                 "If true, the reconstruction of T & S for pressure in "//&
                 "boundary cells is extrapolated, rather than using PCM "//&
                 "in these cells. If true, the same order polynomial is "//&
                 "used as is used for the interior cells.", default=.true.)

  if (CS%tides) then
    CS%id_e_tidal = register_diag_field('ocean_model', 'e_tidal', diag%axesT1, &
        Time, 'Tidal Forcing Astronomical and SAL Height Anomaly', 'meter', conversion=US%Z_to_m)
  endif

  CS%GFS_scale = 1.0
  if (GV%g_prime(1) /= GV%g_Earth) CS%GFS_scale = GV%g_prime(1) / GV%g_Earth

  call log_param(param_file, mdl, "GFS / G_EARTH", CS%GFS_scale)

end subroutine PressureForce_AFV_init

!> Deallocates the finite volume pressure gradient control structure
subroutine PressureForce_AFV_end(CS)
  type(PressureForce_AFV_CS), pointer :: CS !< Finite volume pressure control structure that
                                            !! will be deallocated in this subroutine.
  if (associated(CS)) deallocate(CS)
end subroutine PressureForce_AFV_end

!> \namespace mom_pressureforce_afv
!!
!! Provides the Boussinesq and non-Boussinesq forms of horizontal accelerations
!! due to pressure gradients using a 2nd-order analytically vertically integrated
!! finite volume form, as described by Adcroft et al., 2008.
!!
!! This form eliminates the thermobaric instabilities that had been a problem with
!! previous forms of the pressure gradient force calculation, as described by
!! Hallberg, 2005.
!!
!! Adcroft, A., R. Hallberg, and M. Harrison, 2008: A finite volume discretization
!! of the pressure gradient force using analytic integration. Ocean Modelling, 22,
!! 106-113. http://doi.org/10.1016/j.ocemod.2008.02.001
!!
!! Hallberg, 2005: A thermobaric instability of Lagrangian vertical coordinate
!! ocean models. Ocean Modelling, 8, 279-300.
!! http://dx.doi.org/10.1016/j.ocemod.2004.01.001

end module MOM_PressureForce_AFV
