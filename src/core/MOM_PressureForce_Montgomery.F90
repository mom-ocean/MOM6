!> Provides the Montgomery potential form of pressure gradient
module MOM_PressureForce_Mont

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_density_integrals, only : int_specific_vol_dp
use MOM_diag_mediator, only : post_data, register_diag_field
use MOM_diag_mediator, only : safe_alloc_ptr, diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_tidal_forcing, only : calc_tidal_forcing, tidal_forcing_CS
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs
use MOM_EOS, only : query_compressible

implicit none ; private

#include <MOM_memory.h>

public PressureForce_Mont_Bouss, PressureForce_Mont_nonBouss, Set_pbce_Bouss
public Set_pbce_nonBouss, PressureForce_Mont_init, PressureForce_Mont_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for the Montgomery potential form of pressure gradient
type, public :: PressureForce_Mont_CS ; private
  logical :: tides          !< If true, apply tidal momentum forcing.
  real    :: Rho0           !< The density used in the Boussinesq
                            !! approximation [R ~> kg m-3].
  real    :: GFS_scale      !< Ratio between gravity applied to top interface and the
                            !! gravitational acceleration of the planet [nondim].
                            !! Usually this ratio is 1.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate
                            !! the timing of diagnostic output.
  real, pointer :: PFu_bc(:,:,:) => NULL() !< Zonal accelerations due to pressure gradients
                            !! deriving from density gradients within layers [L T-2 ~> m s-2].
  real, pointer :: PFv_bc(:,:,:) => NULL() !< Meridional accelerations due to pressure gradients
                            !! deriving from density gradients within layers [L T-2 ~> m s-2].
  !>@{ Diagnostic IDs
  integer :: id_PFu_bc = -1, id_PFv_bc = -1, id_e_tidal = -1
  !>@}
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL() !< The tidal forcing control structure
end type PressureForce_Mont_CS

contains

!> \brief Non-Boussinesq Montgomery-potential form of pressure gradient
!!
!! Determines the acceleration due to pressure forces in a
!! non-Boussinesq fluid using the compressibility compensated (if appropriate)
!! Montgomery-potential form described in Hallberg (Ocean Mod., 2005).
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_Mont_nonBouss(h, tv, PFu, PFv, G, GV, US, CS, p_atm, pbce, eta)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness, [H ~> kg m-2].
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: PFu !< Zonal acceleration due to pressure gradients
                                                                 !! (equal to -dM/dx) [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: PFv !< Meridional acceleration due to pressure gradients
                                                                 !! (equal to -dM/dy) [L T-2 ~> m s-2].
  type(PressureForce_Mont_CS),                pointer     :: CS  !< Control structure for Montgomery potential PGF
  real, dimension(:,:),             optional, pointer     :: p_atm !< The pressure at the ice-ocean or
                                                                 !! atmosphere-ocean [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    optional, intent(out) :: pbce !< The baroclinic pressure anomaly in
                                                                 !! each layer due to free surface height anomalies,
                                                                 !! [L2 T-2 H-1 ~> m s-2 or m4 kg-1 s-2].
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(out) :: eta !< Free surface height [H ~> kg m-1].

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    M, &          ! The Montgomery potential, M = (p/rho + gz)  [L2 T-2 ~> m2 s-2].
    alpha_star, & ! Compression adjusted specific volume [R-1 ~> m3 kg-1].
    dz_geo        !   The change in geopotential across a layer [L2 T-2 ~> m2 s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: p ! Interface pressure [R L2 T-2 ~> Pa].
                ! p may be adjusted (with a nonlinear equation of state) so that
                ! its derivative compensates for the adiabatic compressibility
                ! in seawater, but p will still be close to the pressure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [degC].
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [ppt].

  real, dimension(SZI_(G)) :: Rho_cv_BL  ! The coordinate potential density in the
                  ! deepest variable density near-surface layer [R ~> kg m-3].

  real, dimension(SZI_(G),SZJ_(G)) :: &
    dM, &         !   A barotropic correction to the Montgomery potentials to enable the use
                  ! of a reduced gravity form of the equations [L2 T-2 ~> m2 s-2].
    dp_star, &    ! Layer thickness after compensation for compressibility [R L2 T-2 ~> Pa].
    SSH, &        ! The sea surface height anomaly, in depth units [Z ~> m].
    e_tidal, &    !   Bottom geopotential anomaly due to tidal forces from
                  ! astronomical sources and self-attraction and loading [Z ~> m].
    geopot_bot    !   Bottom geopotential relative to time-mean sea level,
                  ! including any tidal contributions [L2 T-2 ~> m2 s-2].
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density [R L2 T-2 ~> Pa] (usually 2e7 Pa = 2000 dbar).
  real :: rho_in_situ(SZI_(G)) !In-situ density of a layer [R ~> kg m-3].
  real :: PFu_bc, PFv_bc     ! The pressure gradient force due to along-layer
                             ! compensated density gradients [L T-2 ~> m s-2]
  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [R L2 T-2 ~> Pa].
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_EOS         ! If true, density is calculated from T & S using an equation of state.
  logical :: is_split        ! A flag indicating whether the pressure gradient terms are to be
                             ! split into barotropic and baroclinic pieces.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.

  real :: I_gEarth           ! The inverse of g_Earth [T2 Z L-2 ~> s2 m-1]
!  real :: dalpha
  real :: Pa_to_H     ! A factor to convert from R L2 T-2 to the thickness units (H).
  real :: alpha_Lay(SZK_(GV)) ! The specific volume of each layer [R-1 ~> m3 kg-1].
  real :: dalpha_int(SZK_(GV)+1) ! The change in specific volume across each
                             ! interface [R-1 ~> m3 kg-1].
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: i, j, k
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  nkmb=GV%nk_rho_varies
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)

  use_p_atm = .false.
  if (present(p_atm)) then ; if (associated(p_atm)) use_p_atm = .true. ; endif
  is_split = .false. ; if (present(pbce)) is_split = .true.
  use_EOS = associated(tv%eqn_of_state)

  if (.not.associated(CS)) call MOM_error(FATAL, &
      "MOM_PressureForce_Mont: Module must be initialized before it is used.")
  if (use_EOS) then
    if (query_compressible(tv%eqn_of_state)) call MOM_error(FATAL, &
      "PressureForce_Mont_nonBouss: The Montgomery form of the pressure force "//&
      "can no longer be used with a compressible EOS. Use #define ANALYTIC_FV_PGF.")
  endif

  I_gEarth = 1.0 / GV%g_Earth
  dp_neglect = GV%g_Earth * GV%H_to_RZ * GV%H_subroundoff
  do k=1,nz ; alpha_Lay(k) = 1.0 / (GV%Rlay(k)) ; enddo
  do k=2,nz ; dalpha_int(K) = alpha_Lay(k-1) - alpha_Lay(k) ; enddo

  if (use_p_atm) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1 ; p(i,j,1) = p_atm(i,j) ; enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1 ; p(i,j,1) = 0.0 ; enddo ; enddo
  endif
  !$OMP parallel do default(shared)
  do j=Jsq,Jeq+1 ; do k=1,nz ; do i=Isq,Ieq+1
    p(i,j,K+1) = p(i,j,K) + GV%g_Earth * GV%H_to_RZ * h(i,j,k)
  enddo ; enddo ; enddo

  if (present(eta)) then
    Pa_to_H = 1.0 / (GV%g_Earth * GV%H_to_RZ)
    if (use_p_atm) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = (p(i,j,nz+1) - p_atm(i,j)) * Pa_to_H ! eta has the same units as h.
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = p(i,j,nz+1) * Pa_to_H ! eta has the same units as h.
      enddo ; enddo
    endif
  endif

  if (CS%tides) then
    !   Determine the sea surface height anomalies, to enable the calculation
    ! of self-attraction and loading.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      SSH(i,j) = -G%bathyT(i,j)
    enddo ; enddo
    if (use_EOS) then
      !$OMP parallel do default(shared)
      do k=1,nz
        call int_specific_vol_dp(tv%T(:,:,k), tv%S(:,:,k), p(:,:,k), p(:,:,k+1), &
                                 0.0, G%HI, tv%eqn_of_state, US, dz_geo(:,:,k), halo_size=1)
      enddo
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do k=1,nz ; do i=Isq,Ieq+1
        SSH(i,j) = SSH(i,j) + I_gEarth * dz_geo(i,j,k)
      enddo ; enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1 ; do k=1,nz ; do i=Isq,Ieq+1
        SSH(i,j) = SSH(i,j) + GV%H_to_RZ * h(i,j,k) * alpha_Lay(k)
      enddo ; enddo ; enddo
    endif

    call calc_tidal_forcing(CS%Time, SSH, e_tidal, G, CS%tides_CSp, m_to_Z=US%m_to_Z)
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      geopot_bot(i,j) = -GV%g_Earth*(e_tidal(i,j) + G%bathyT(i,j))
    enddo ; enddo
  else
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      geopot_bot(i,j) = -GV%g_Earth*G%bathyT(i,j)
    enddo ; enddo
  endif

  if (use_EOS) then
    !   Calculate in-situ specific volumes (alpha_star).

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
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, Rho_cv_BL(:), &
                               tv%eqn_of_state, EOSdom)
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
      do i=Isq,Ieq+1 ; p_ref(i) = 0 ; enddo
    endif
    !$OMP parallel do default(shared) private(rho_in_situ)
    do k=1,nz ; do j=Jsq,Jeq+1
      call calculate_density(tv_tmp%T(:,j,k), tv_tmp%S(:,j,k), p_ref, rho_in_situ, &
                             tv%eqn_of_state, EOSdom)
      do i=Isq,Ieq+1 ; alpha_star(i,j,k) = 1.0 / rho_in_situ(i) ; enddo
    enddo ; enddo
  endif                                               ! use_EOS

  if (use_EOS) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        M(i,j,nz) = geopot_bot(i,j) + p(i,j,nz+1) * alpha_star(i,j,nz)
      enddo
      do k=nz-1,1,-1 ; do i=Isq,Ieq+1
        M(i,j,k) = M(i,j,k+1) + p(i,j,K+1) * (alpha_star(i,j,k) - alpha_star(i,j,k+1))
      enddo ; enddo
    enddo
  else ! not use_EOS
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        M(i,j,nz) = geopot_bot(i,j) + p(i,j,nz+1) * alpha_Lay(nz)
      enddo
      do k=nz-1,1,-1 ; do i=Isq,Ieq+1
        M(i,j,k) = M(i,j,k+1) + p(i,j,K+1) * dalpha_int(K+1)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (CS%GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      dM(i,j) = (CS%GFS_scale - 1.0) * M(i,j,1)
    enddo ; enddo
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      M(i,j,k) = M(i,j,k) + dM(i,j)
    enddo ; enddo ; enddo

    !   Could instead do the following, to avoid taking small differences
    ! of large numbers...
!   do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
!     M(i,j,1) = CS%GFS_scale * M(i,j,1)
!   enddo ; enddo
!   if (use_EOS) then
!     do k=2,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
!       M(i,j,k) = M(i,j,k-1) - p(i,j,K) * (alpha_star(i,j,k-1) - alpha_star(i,j,k))
!     enddo ; enddo ; enddo
!   else ! not use_EOS
!     do k=2,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
!        M(i,j,k) = M(i,j,k-1) - p(i,j,K) * dalpha_int(K)
!     enddo ; enddo ; enddo
!   endif ! use_EOS

  endif

  ! Note that ddM/dPb = alpha_star(i,j,1)
  if (present(pbce)) then
    call Set_pbce_nonBouss(p, tv_tmp, G, GV, US, CS%GFS_scale, pbce, alpha_star)
  endif

!    Calculate the pressure force. On a Cartesian grid,
!      PFu = - dM/dx   and  PFv = - dM/dy.
  if (use_EOS) then
    !$OMP parallel do default(shared) private(dp_star,PFu_bc,PFv_bc)
    do k=1,nz
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dp_star(i,j) = (p(i,j,K+1) - p(i,j,K)) + dp_neglect
      enddo ; enddo
      do j=js,je ; do I=Isq,Ieq
        ! PFu_bc = p* grad alpha*
        PFu_bc = (alpha_star(i+1,j,k) - alpha_star(i,j,k)) * (G%IdxCu(I,j) * &
            ((dp_star(i,j)*dp_star(i+1,j) + (p(i,j,K)*dp_star(i+1,j) + p(i+1,j,K)*dp_star(i,j))) / &
             (dp_star(i,j) + dp_star(i+1,j))))
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * G%IdxCu(I,j) + PFu_bc
        if (associated(CS%PFu_bc)) CS%PFu_bc(i,j,k) = PFu_bc
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv_bc = (alpha_star(i,j+1,k) - alpha_star(i,j,k)) * (G%IdyCv(i,J) * &
            ((dp_star(i,j)*dp_star(i,j+1) + (p(i,j,K)*dp_star(i,j+1) + p(i,j+1,K)*dp_star(i,j))) / &
             (dp_star(i,j) + dp_star(i,j+1))))
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * G%IdyCv(i,J) + PFv_bc
        if (associated(CS%PFv_bc)) CS%PFv_bc(i,j,k) = PFv_bc
      enddo ; enddo
    enddo ! k-loop
  else ! .not. use_EOS
    !$OMP parallel do default(shared)
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * G%IdxCu(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * G%IdyCv(i,J)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (CS%id_PFu_bc>0) call post_data(CS%id_PFu_bc, CS%PFu_bc, CS%diag)
  if (CS%id_PFv_bc>0) call post_data(CS%id_PFv_bc, CS%PFv_bc, CS%diag)
  if (CS%id_e_tidal>0) call post_data(CS%id_e_tidal, e_tidal, CS%diag)

end subroutine PressureForce_Mont_nonBouss

!> \brief Boussinesq Montgomery-potential form of pressure gradient
!!
!! Determines the acceleration due to pressure forces.
!!
!! To work, the following fields must be set outside of the usual (is:ie,js:je)
!! range before this subroutine is called:
!!   h(isB:ie+1,jsB:je+1), T(isB:ie+1,jsB:je+1), and S(isB:ie+1,jsB:je+1).
subroutine PressureForce_Mont_Bouss(h, tv, PFu, PFv, G, GV, US, CS, p_atm, pbce, eta)
  type(ocean_grid_type),                      intent(in)  :: G   !< Ocean grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV  !< Vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)  :: h   !< Layer thickness [H ~> m].
  type(thermo_var_ptrs),                      intent(in)  :: tv  !< Thermodynamic variables.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(out) :: PFu !< Zonal acceleration due to pressure gradients
                                                                 !! (equal to -dM/dx) [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(out) :: PFv !< Meridional acceleration due to pressure gradients
                                                                 !! (equal to -dM/dy) [L T-2 ~> m s2].
  type(PressureForce_Mont_CS),                pointer     :: CS  !< Control structure for Montgomery potential PGF
  real, dimension(:,:),                      optional, pointer     :: p_atm !< The pressure at the ice-ocean or
                                                                !! atmosphere-ocean [R L2 T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), optional, intent(out) :: pbce !< The baroclinic pressure anomaly in
                                                                !! each layer due to free surface height anomalies
                                                                !! [L2 T-2 H-1 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< Free surface height [H ~> m].
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    M, &        ! The Montgomery potential, M = (p/rho + gz) [L2 T-2 ~> m2 s-2].
    rho_star    ! In-situ density divided by the derivative with depth of the
                ! corrected e times (G_Earth/Rho0) [m2 Z-1 s-2 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: e ! Interface height in m.
                ! e may be adjusted (with a nonlinear equation of state) so that
                ! its derivative compensates for the adiabatic compressibility
                ! in seawater, but e will still be close to the interface depth.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [degC].
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties [ppt].

  real :: Rho_cv_BL(SZI_(G)) !   The coordinate potential density in
                ! the deepest variable density near-surface layer [R ~> kg m-3].
  real :: h_star(SZI_(G),SZJ_(G)) ! Layer thickness after compensation
                             ! for compressibility [Z ~> m].
  real :: e_tidal(SZI_(G),SZJ_(G)) ! Bottom geopotential anomaly due to tidal
                             ! forces from astronomical sources and self-
                             ! attraction and loading, in depth units [Z ~> m].
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density [R L2 T-2 ~> Pa] (usually 2e7 Pa = 2000 dbar).
  real :: I_Rho0             ! 1/Rho0 [R-1 ~> m3 kg-1].
  real :: G_Rho0             ! G_Earth / Rho0 [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1].
  real :: PFu_bc, PFv_bc     ! The pressure gradient force due to along-layer
                             ! compensated density gradients [L T-2 ~> m s-2]
  real :: h_neglect          ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [Z ~> m].
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_EOS         ! If true, density is calculated from T & S using
                             ! an equation of state.
  logical :: is_split        ! A flag indicating whether the pressure
                             ! gradient terms are to be split into
                             ! barotropic and baroclinic pieces.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: i, j, k

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  nkmb=GV%nk_rho_varies
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)

  use_p_atm = .false.
  if (present(p_atm)) then ; if (associated(p_atm)) use_p_atm = .true. ; endif
  is_split = .false. ; if (present(pbce)) is_split = .true.
  use_EOS = associated(tv%eqn_of_state)

  if (.not.associated(CS)) call MOM_error(FATAL, &
       "MOM_PressureForce_Mont: Module must be initialized before it is used.")
  if (use_EOS) then
    if (query_compressible(tv%eqn_of_state)) call MOM_error(FATAL, &
      "PressureForce_Mont_Bouss: The Montgomery form of the pressure force "//&
      "can no longer be used with a compressible EOS. Use #define ANALYTIC_FV_PGF.")
  endif

  h_neglect = GV%H_subroundoff * GV%H_to_Z
  I_Rho0 = 1.0/CS%Rho0
  G_Rho0 = GV%g_Earth / GV%Rho0

  if (CS%tides) then
    !   Determine the surface height anomaly for calculating self attraction
    ! and loading.  This should really be based on bottom pressure anomalies,
    ! but that is not yet implemented, and the current form is correct for
    ! barotropic tides.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1 ; e(i,j,1) = -G%bathyT(i,j) ; enddo
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
  do j=Jsq,Jeq+1 ; do k=nz,1,-1 ; do i=Isq,Ieq+1
    e(i,j,K) = e(i,j,K+1) + h(i,j,k)*GV%H_to_Z
  enddo ; enddo ; enddo

  if (use_EOS) then
!   Calculate in-situ densities (rho_star).

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
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, Rho_cv_BL(:), &
                               tv%eqn_of_state, EOSdom)

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
      do i=Isq,Ieq+1 ; p_ref(i) = 0.0 ; enddo
    endif

    ! This no longer includes any pressure dependency, since this routine
    ! will come down with a fatal error if there is any compressibility.
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=Jsq,Jeq+1
      call calculate_density(tv_tmp%T(:,j,k), tv_tmp%S(:,j,k), p_ref, rho_star(:,j,k), &
                             tv%eqn_of_state, EOSdom)
      do i=Isq,Ieq+1 ; rho_star(i,j,k) = G_Rho0*rho_star(i,j,k) ; enddo
    enddo ; enddo
  endif                                               ! use_EOS

!    Here the layer Montgomery potentials, M, are calculated.
  if (use_EOS) then
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        M(i,j,1) = CS%GFS_scale * (rho_star(i,j,1) * e(i,j,1))
        if (use_p_atm) M(i,j,1) = M(i,j,1) + p_atm(i,j) * I_Rho0
      enddo
      do k=2,nz ; do i=Isq,Ieq+1
        M(i,j,k) = M(i,j,k-1) + (rho_star(i,j,k) - rho_star(i,j,k-1)) * e(i,j,K)
      enddo ; enddo
    enddo
  else ! not use_EOS
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        M(i,j,1) = GV%g_prime(1) * e(i,j,1)
        if (use_p_atm) M(i,j,1) = M(i,j,1) + p_atm(i,j) * I_Rho0
      enddo
      do k=2,nz ; do i=Isq,Ieq+1
        M(i,j,k) = M(i,j,k-1) + GV%g_prime(K) * e(i,j,K)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (present(pbce)) then
    call Set_pbce_Bouss(e, tv_tmp, G, GV, US, CS%Rho0, CS%GFS_scale, pbce, rho_star)
  endif

!    Calculate the pressure force. On a Cartesian grid,
!      PFu = - dM/dx   and  PFv = - dM/dy.
  if (use_EOS) then
    !$OMP parallel do default(shared) private(h_star,PFu_bc,PFv_bc)
    do k=1,nz
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        h_star(i,j) = (e(i,j,K) - e(i,j,K+1)) + h_neglect
      enddo ; enddo
      do j=js,je ; do I=Isq,Ieq
        PFu_bc = -1.0*(rho_star(i+1,j,k) - rho_star(i,j,k)) * (G%IdxCu(I,j) * &
          ((h_star(i,j) * h_star(i+1,j) - (e(i,j,K) * h_star(i+1,j) + &
          e(i+1,j,K) * h_star(i,j))) / (h_star(i,j) + h_star(i+1,j))))
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * G%IdxCu(I,j) + PFu_bc
        if (associated(CS%PFu_bc)) CS%PFu_bc(i,j,k) = PFu_bc
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv_bc = -1.0*(rho_star(i,j+1,k) - rho_star(i,j,k)) * (G%IdyCv(i,J) * &
          ((h_star(i,j) * h_star(i,j+1) - (e(i,j,K) * h_star(i,j+1) + &
          e(i,j+1,K) * h_star(i,j))) / (h_star(i,j) + h_star(i,j+1))))
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * G%IdyCv(i,J) + PFv_bc
        if (associated(CS%PFv_bc)) CS%PFv_bc(i,j,k) = PFv_bc
      enddo ; enddo
    enddo ! k-loop
  else ! .not. use_EOS
    !$OMP parallel do default(shared)
    do k=1,nz
      do j=js,je ; do I=Isq,Ieq
        PFu(I,j,k) = -(M(i+1,j,k) - M(i,j,k)) * G%IdxCu(I,j)
      enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie
        PFv(i,J,k) = -(M(i,j+1,k) - M(i,j,k)) * G%IdyCv(i,J)
      enddo ; enddo
    enddo
  endif ! use_EOS

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

  if (CS%id_PFu_bc>0) call post_data(CS%id_PFu_bc, CS%PFu_bc, CS%diag)
  if (CS%id_PFv_bc>0) call post_data(CS%id_PFv_bc, CS%PFv_bc, CS%diag)
  if (CS%id_e_tidal>0) call post_data(CS%id_e_tidal, e_tidal, CS%diag)

end subroutine PressureForce_Mont_Bouss

!> Determines the partial derivative of the acceleration due
!! to pressure forces with the free surface height.
subroutine Set_pbce_Bouss(e, tv, G, GV, US, Rho0, GFS_scale, pbce, rho_star)
  type(ocean_grid_type),                intent(in)  :: G    !< Ocean grid structure
  type(verticalGrid_type),              intent(in)  :: GV   !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: e !< Interface height [Z ~> m].
  type(thermo_var_ptrs),                intent(in)  :: tv   !< Thermodynamic variables
  type(unit_scale_type),                intent(in)  :: US   !< A dimensional unit scaling type
  real,                                 intent(in)  :: Rho0 !< The "Boussinesq" ocean density [R ~> kg m-3].
  real,                                 intent(in)  :: GFS_scale !< Ratio between gravity applied to top
                                                            !! interface and the gravitational acceleration of
                                                            !! the planet [nondim]. Usually this ratio is 1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                        intent(out) :: pbce !< The baroclinic pressure anomaly in each layer due
                                                            !! to free surface height anomalies
                                                            !! [L2 T-2 H-1 ~> m s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              optional, intent(in)  :: rho_star !< The layer densities (maybe compressibility
                                                            !! compensated), times g/rho_0 [L2 Z-1 T-2 ~> m s-2].

  ! Local variables
  real :: Ihtot(SZI_(G))     ! The inverse of the sum of the layer thicknesses [H-1 ~> m-1 or m2 kg-1].
  real :: press(SZI_(G))     ! Interface pressure [R L2 T-2 ~> Pa].
  real :: T_int(SZI_(G))     ! Interface temperature [degC].
  real :: S_int(SZI_(G))     ! Interface salinity [ppt].
  real :: dR_dT(SZI_(G))     ! Partial derivative of density with temperature [R degC-1 ~> kg m-3 degC-1].
  real :: dR_dS(SZI_(G))     ! Partial derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1].
  real :: rho_in_situ(SZI_(G)) ! In-situ density at the top of a layer [R ~> kg m-3].
  real :: G_Rho0             ! A scaled version of g_Earth / Rho0 [L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]
  real :: Rho0xG             ! g_Earth * Rho0 [kg s-2 m-1 Z-1 ~> kg s-2 m-2]
  logical :: use_EOS         ! If true, density is calculated from T & S using
                             ! an equation of state.
  real :: z_neglect          ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [Z ~> m].
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: Isq, Ieq, Jsq, Jeq, nz, i, j, k

  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke
  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)

  Rho0xG = Rho0 * GV%g_Earth
  G_Rho0 = GV%g_Earth / GV%Rho0
  use_EOS = associated(tv%eqn_of_state)
  z_neglect = GV%H_subroundoff*GV%H_to_Z

  if (use_EOS) then
    if (present(rho_star)) then
     !$OMP parallel do default(shared) private(Ihtot)
      do j=Jsq,Jeq+1
        do i=Isq,Ieq+1
          Ihtot(i) = GV%H_to_Z / ((e(i,j,1)-e(i,j,nz+1)) + z_neglect)
          pbce(i,j,1) = GFS_scale * rho_star(i,j,1) * GV%H_to_Z
        enddo
        do k=2,nz ; do i=Isq,Ieq+1
          pbce(i,j,k) = pbce(i,j,k-1) + (rho_star(i,j,k)-rho_star(i,j,k-1)) * &
                        ((e(i,j,K) - e(i,j,nz+1)) * Ihtot(i))
        enddo ; enddo
      enddo ! end of j loop
    else
      !$OMP parallel do default(shared) private(Ihtot,press,rho_in_situ,T_int,S_int,dR_dT,dR_dS)
      do j=Jsq,Jeq+1
        do i=Isq,Ieq+1
          Ihtot(i) = GV%H_to_Z / ((e(i,j,1)-e(i,j,nz+1)) + z_neglect)
          press(i) = -Rho0xG*e(i,j,1)
        enddo
        call calculate_density(tv%T(:,j,1), tv%S(:,j,1), press, rho_in_situ, &
                               tv%eqn_of_state, EOSdom)
        do i=Isq,Ieq+1
          pbce(i,j,1) = G_Rho0*(GFS_scale * rho_in_situ(i)) * GV%H_to_Z
        enddo
        do k=2,nz
          do i=Isq,Ieq+1
            press(i) = -Rho0xG*e(i,j,K)
            T_int(i) = 0.5*(tv%T(i,j,k-1)+tv%T(i,j,k))
            S_int(i) = 0.5*(tv%S(i,j,k-1)+tv%S(i,j,k))
          enddo
          call calculate_density_derivs(T_int, S_int, press, dR_dT, dR_dS, &
                                        tv%eqn_of_state, EOSdom)
          do i=Isq,Ieq+1
            pbce(i,j,k) = pbce(i,j,k-1) + G_Rho0 * &
               ((e(i,j,K) - e(i,j,nz+1)) * Ihtot(i)) * &
               (dR_dT(i)*(tv%T(i,j,k)-tv%T(i,j,k-1)) + &
                dR_dS(i)*(tv%S(i,j,k)-tv%S(i,j,k-1)))
          enddo
        enddo
      enddo ! end of j loop
    endif
  else ! not use_EOS
    !$OMP parallel do default(shared) private(Ihtot)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        Ihtot(i) = 1.0 / ((e(i,j,1)-e(i,j,nz+1)) + z_neglect)
        pbce(i,j,1) = GV%g_prime(1) * GV%H_to_Z
      enddo
      do k=2,nz ; do i=Isq,Ieq+1
        pbce(i,j,k) = pbce(i,j,k-1) + &
                      (GV%g_prime(K)*GV%H_to_Z) * ((e(i,j,K) - e(i,j,nz+1)) * Ihtot(i))
      enddo ; enddo
    enddo ! end of j loop
  endif ! use_EOS

end subroutine Set_pbce_Bouss

!> Determines the partial derivative of the acceleration due
!! to pressure forces with the column mass.
subroutine Set_pbce_nonBouss(p, tv, G, GV, US, GFS_scale, pbce, alpha_star)
  type(ocean_grid_type),                intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type),              intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: p !< Interface pressures [R L2 T-2 ~> Pa].
  type(thermo_var_ptrs),                intent(in)  :: tv !< Thermodynamic variables
  type(unit_scale_type),                intent(in)  :: US !< A dimensional unit scaling type
  real,                                 intent(in)  :: GFS_scale !< Ratio between gravity applied to top
                                                          !! interface and the gravitational acceleration of
                                                          !! the planet [nondim]. Usually this ratio is 1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: pbce !< The baroclinic pressure anomaly in each
                                                          !! layer due to free surface height anomalies
                                                          !! [L2 H-1 T-2 ~> m4 kg-1 s-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), optional, intent(in) :: alpha_star !< The layer specific volumes
                                                          !! (maybe compressibility compensated) [R-1 ~> m3 kg-1].
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    dpbce, &      !   A barotropic correction to the pbce to enable the use of
                  ! a reduced gravity form of the equations [L2 H-1 T-2 ~> m4 kg-1 s-2].
    C_htot        ! dP_dH divided by the total ocean pressure [H-1 ~> m2 kg-1].
  real :: T_int(SZI_(G))     ! Interface temperature [degC].
  real :: S_int(SZI_(G))     ! Interface salinity [ppt].
  real :: dR_dT(SZI_(G))     ! Partial derivative of density with temperature [R degC-1 ~> kg m-3 degC-1].
  real :: dR_dS(SZI_(G))     ! Partial derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1].
  real :: rho_in_situ(SZI_(G)) ! In-situ density at an interface [R ~> kg m-3].
  real :: alpha_Lay(SZK_(GV)) ! The specific volume of each layer [R-1 ~> m3 kg-1].
  real :: dalpha_int(SZK_(GV)+1) ! The change in specific volume across each interface [R-1 ~> m3 kg-1].
  real :: dP_dH              ! A factor that converts from thickness to pressure times other dimensional
                             ! conversion factors [R L2 T-2 H-1 ~> Pa m2 kg-1].
  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected [R L2 T-2 ~> Pa].
  logical :: use_EOS         ! If true, density is calculated from T & S using an equation of state.
  integer, dimension(2) :: EOSdom ! The computational domain for the equation of state
  integer :: Isq, Ieq, Jsq, Jeq, nz, i, j, k

  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB ; nz = GV%ke
  EOSdom(1) = Isq - (G%isd-1) ;  EOSdom(2) = G%iec+1 - (G%isd-1)

  use_EOS = associated(tv%eqn_of_state)

  dP_dH = GV%g_Earth * GV%H_to_RZ
  dp_neglect = GV%g_Earth * GV%H_to_RZ * GV%H_subroundoff

  if (use_EOS) then
    if (present(alpha_star)) then
      !$OMP parallel do default(shared)
      do j=Jsq,Jeq+1
        do i=Isq,Ieq+1
          C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
          pbce(i,j,nz) = dP_dH * alpha_star(i,j,nz)
        enddo
        do k=nz-1,1,-1 ; do i=Isq,Ieq+1
          pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,K+1)-p(i,j,1)) * C_htot(i,j)) * &
              (alpha_star(i,j,k) - alpha_star(i,j,k+1))
        enddo ; enddo
      enddo
    else
      !$OMP parallel do default(shared) private(T_int,S_int,dR_dT,dR_dS,rho_in_situ)
      do j=Jsq,Jeq+1
        call calculate_density(tv%T(:,j,nz), tv%S(:,j,nz), p(:,j,nz+1), rho_in_situ, &
                               tv%eqn_of_state, EOSdom)
        do i=Isq,Ieq+1
          C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
          pbce(i,j,nz) = dP_dH / (rho_in_situ(i))
        enddo
        do k=nz-1,1,-1
          do i=Isq,Ieq+1
            T_int(i) = 0.5*(tv%T(i,j,k)+tv%T(i,j,k+1))
            S_int(i) = 0.5*(tv%S(i,j,k)+tv%S(i,j,k+1))
          enddo
          call calculate_density(T_int, S_int, p(:,j,k+1), rho_in_situ, tv%eqn_of_state, EOSdom)
          call calculate_density_derivs(T_int, S_int, p(:,j,k+1), dR_dT, dR_dS, &
                                        tv%eqn_of_state, EOSdom)
          do i=Isq,Ieq+1
            pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,K+1)-p(i,j,1))*C_htot(i,j)) *  &
                ((dR_dT(i)*(tv%T(i,j,k+1)-tv%T(i,j,k)) + &
                  dR_dS(i)*(tv%S(i,j,k+1)-tv%S(i,j,k))) / (rho_in_situ(i)**2))
          enddo
        enddo
      enddo
    endif
  else ! not use_EOS

    do k=1,nz ; alpha_Lay(k) = 1.0 / (GV%Rlay(k)) ; enddo
    do k=2,nz ; dalpha_int(K) = alpha_Lay(k-1) - alpha_Lay(k) ; enddo

    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        C_htot(i,j) = dP_dH / ((p(i,j,nz+1)-p(i,j,1)) + dp_neglect)
        pbce(i,j,nz) = dP_dH * alpha_Lay(nz)
      enddo
      do k=nz-1,1,-1 ; do i=Isq,Ieq+1
        pbce(i,j,k) = pbce(i,j,k+1) + ((p(i,j,K+1)-p(i,j,1))*C_htot(i,j)) * dalpha_int(K+1)
      enddo ; enddo
    enddo
  endif ! use_EOS

  if (GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    !$OMP parallel do default(shared)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      dpbce(i,j) = (GFS_scale - 1.0) * pbce(i,j,1)
    enddo ; enddo
    !$OMP parallel do default(shared)
    do k=1,nz ; do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      pbce(i,j,k) = pbce(i,j,k) + dpbce(i,j)
    enddo ; enddo ; enddo
  endif

end subroutine Set_pbce_nonBouss

!> Initialize the Montgomery-potential form of PGF control structure
subroutine PressureForce_Mont_init(Time, G, GV, US, param_file, diag, CS, tides_CSp)
  type(time_type), target, intent(in)    :: Time !< Current model time
  type(ocean_grid_type),   intent(in)    :: G  !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target, intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_Mont_CS),  pointer  :: CS !< Montgomery PGF control structure
  type(tidal_forcing_CS), optional, pointer :: tides_CSp !< Tides control structure

  ! Local variables
  logical :: use_temperature, use_EOS
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl   ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "PressureForce_init called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag ; CS%Time => Time
  if (present(tides_CSp)) then
    if (associated(tides_CSp)) CS%tides_CSp => tides_CSp
  endif

  mdl = "MOM_PressureForce_Mont"
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to "//&
                 "calculate accelerations and the mass for conservation "//&
                 "properties, or with BOUSSINSEQ false to convert some "//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0, scale=US%R_to_kg_m3)
  call get_param(param_file, mdl, "TIDES", CS%tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  call get_param(param_file, mdl, "USE_EOS", use_EOS, default=.true., &
                 do_not_log=.true.) ! Input for diagnostic use only.

  if (use_EOS) then
    CS%id_PFu_bc = register_diag_field('ocean_model', 'PFu_bc', diag%axesCuL, Time, &
         'Density Gradient Zonal Pressure Force Accel.', "meter second-2", conversion=US%L_T2_to_m_s2)
    CS%id_PFv_bc = register_diag_field('ocean_model', 'PFv_bc', diag%axesCvL, Time, &
         'Density Gradient Meridional Pressure Force Accel.', "meter second-2", conversion=US%L_T2_to_m_s2)
    if (CS%id_PFu_bc > 0) then
      call safe_alloc_ptr(CS%PFu_bc,G%IsdB,G%IedB,G%jsd,G%jed,GV%ke)
      CS%PFu_bc(:,:,:) = 0.0
    endif
    if (CS%id_PFv_bc > 0) then
      call safe_alloc_ptr(CS%PFv_bc,G%isd,G%ied,G%JsdB,G%JedB,GV%ke)
      CS%PFv_bc(:,:,:) = 0.0
    endif
  endif

  if (CS%tides) then
    CS%id_e_tidal = register_diag_field('ocean_model', 'e_tidal', diag%axesT1, &
        Time, 'Tidal Forcing Astronomical and SAL Height Anomaly', 'meter', conversion=US%Z_to_m)
  endif

  CS%GFS_scale = 1.0
  if (GV%g_prime(1) /= GV%g_Earth) CS%GFS_scale = GV%g_prime(1) / GV%g_Earth

  call log_param(param_file, mdl, "GFS / G_EARTH", CS%GFS_scale)

end subroutine PressureForce_Mont_init

!> Deallocates the Montgomery-potential form of PGF control structure
subroutine PressureForce_Mont_end(CS)
  type(PressureForce_Mont_CS), pointer :: CS  !< Control structure for Montgomery potential PGF
  if (associated(CS)) deallocate(CS)
end subroutine PressureForce_Mont_end

!>\namespace mom_pressureforce_mont
!!
!! Provides the Boussunesq and non-Boussinesq forms of the horizontal
!! accelerations due to pressure gradients using the Montgomery potential. A
!! second-order accurate, centered scheme is used. If a split time stepping
!! scheme is used, the vertical decomposition into barotropic and baroclinic
!! contributions described by Hallberg (J Comp Phys 1997) is used.  With a
!! nonlinear equation of state, compressibility is added along the lines proposed
!! by Sun et al. (JPO 1999), but with compressibility coefficients based on a fit
!! to a user-provided reference profile.

end module MOM_PressureForce_Mont
