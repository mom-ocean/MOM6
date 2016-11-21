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
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs
use MOM_EOS, only : int_density_dz, int_specific_vol_dp
use MOM_EOS, only : int_density_dz_generic_plm, int_density_dz_generic_ppm
use MOM_EOS, only : int_density_dz_generic_plm_analytic
use MOM_ALE, only : pressure_gradient_plm, pressure_gradient_ppm
use MOM_ALE, only : usePressureReconstruction, pressureReconstructionScheme
use MOM_ALE, only : ALE_CS
use regrid_defs, only: PRESSURE_RECONSTRUCTION_PLM, PRESSURE_RECONSTRUCTION_PPM

implicit none ; private

#include <MOM_memory.h>

public PressureForce_AFV, PressureForce_AFV_init, PressureForce_AFV_end
public PressureForce_AFV_Bouss, PressureForce_AFV_nonBouss

!> Finite volume pressure gradient control structure
type, public :: PressureForce_AFV_CS ; private
  logical :: tides          !< If true, apply tidal momentum forcing.
  real    :: Rho0           !< The density used in the Boussinesq
                            !! approximation, in kg m-3.
  real    :: GFS_scale      !< A scaling of the surface pressure gradients to
                            !! allow the use of a reduced gravity model.
  type(time_type), pointer :: Time !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                            !! timing of diagnostic output.
  logical :: useMassWghtInterp !< Use mass weighting in T/S interpolation
  integer :: id_e_tidal = -1 !< Diagnostic identifier
  type(tidal_forcing_CS), pointer :: tides_CSp => NULL() !< Tides control structure
end type PressureForce_AFV_CS

contains

!> Thin interface between the model and the Boussinesq and non-Boussinesq
!! pressure force routines.
subroutine PressureForce_AFV(h, tv, PFu, PFv, G, GV, CS, ALE_CSp, p_atm, pbce, eta)
  type(ocean_grid_type),                     intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h   !< Layer thickness (m or kg/m2)
  type(thermo_var_ptrs),                     intent(inout) :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out)   :: PFu !< Zonal acceleration (m/s2)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out)   :: PFv !< Meridional acceleration (m/s2)
  type(PressureForce_AFV_CS),                pointer       :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                              pointer       :: ALE_CSp !< ALE control structure
  real, dimension(:,:),                      optional, pointer :: p_atm !< The pressure at the ice-ocean
                                                           !! or atmosphere-ocean interface in Pa.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  optional, intent(out) :: pbce !< The baroclinic pressure
                                                           !! anomaly in each layer due to eta anomalies,
                                                           !! in m2 s-2 H-1.
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< The bottom mass used to
                                                           !! calculate PFu and PFv, in H, with any tidal
                                                           !! contributions or compressibility compensation.

  if (GV%Boussinesq) then
    call PressureForce_AFV_bouss(h, tv, PFu, PFv, G, GV, CS, ALE_CSp, p_atm, pbce, eta)
  else
    call PressureForce_AFV_nonbouss(h, tv, PFu, PFv, G, GV, CS, p_atm, pbce, eta)
  endif

end subroutine PressureForce_AFV

!> \brief Non-Boussinesq analytically-integrated finite volume form of pressure gradient
!!
!! Determines the acceleration due to hydrostatic pressure forces, using
!! the analytic finite volume form of the Pressure gradient, and does not
!! make the Boussinesq approximation.
!!
!! To work, the following fields must be set outside of the usual
!! ie to ie, je to je range before this subroutine is called:
!!  h[ie+1] and h[je+1] and (if tv%eqn_of_state is set) T[ie+1], S[ie+1],
!!  T[je+1], and S[je+1].
subroutine PressureForce_AFV_nonBouss(h, tv, PFu, PFv, G, GV, CS, p_atm, pbce, eta)
  type(ocean_grid_type),                     intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV  !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h   !< Layer thickness (kg/m2)
  type(thermo_var_ptrs),                     intent(in)    :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out)   :: PFu !< Zonal acceleration (m/s2)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out)   :: PFv !< Meridional acceleration (m/s2)
  type(PressureForce_AFV_CS),                pointer       :: CS  !< Finite volume PGF control structure
  real, dimension(:,:),                      optional, pointer :: p_atm !< The pressure at the ice-ocean
                                                           !! or atmosphere-ocean interface in Pa.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  optional, intent(out) :: pbce !< The baroclinic pressure
                                                           !! anomaly in each layer due to eta anomalies,
                                                           !! in m2 s-2 H-1.
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< The bottom mass used to
                                                           !! calculate PFu and PFv, in H, with any tidal
                                                           !! contributions or compressibility compensation.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: p ! Interface pressure in Pa.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties, in C.
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties, in psu.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G))  :: &
    dza, &      ! The change in geopotential anomaly between the top and bottom
                ! of a layer, in m2 s-2.
    intp_dza    ! The vertical integral in depth of the pressure anomaly less
                ! the pressure anomaly at the top of the layer, in Pa m2 s-2.
  real, dimension(SZI_(G),SZJ_(G))  :: &
    dp, &       ! The (positive) change in pressure across a layer, in Pa.
    SSH, &      ! The sea surface height anomaly, in m.
    e_tidal, &  ! The bottom geopotential anomaly due to tidal forces from
                ! astronomical sources and self-attraction and loading, in m.
    dM, &       ! The barotropic adjustment to the Montgomery potential to
                ! account for a reduced gravity model, in m2 s-2.
    za          ! The geopotential anomaly (i.e. g*e + alpha_0*pressure) at the
                ! interface atop a layer, in m2 s-2.
  real, dimension(SZDI_(G%Block(1)),SZDJ_(G%Block(1))) :: & ! on block indices
    dp_bk, &    ! The (positive) change in pressure across a layer, in Pa.
    za_bk       ! The geopotential anomaly (i.e. g*e + alpha_0*pressure) at the
                ! interface atop a layer, in m2 s-2.

  real, dimension(SZI_(G)) :: Rho_cv_BL !  The coordinate potential density in the deepest variable
                ! density near-surface layer, in kg m-3.
  real, dimension(SZDIB_(G%Block(1)),SZDJ_(G%Block(1))) :: & ! on block indices
    intx_za_bk ! The zonal integral of the geopotential anomaly along the
               ! interface below a layer, divided by the grid spacing, m2 s-2.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)) :: &
    intx_dza    ! The change in intx_za through a layer, in m2 s-2.
  real, dimension(SZDI_(G%Block(1)),SZDJB_(G%Block(1))) :: & ! on block indices
    inty_za_bk ! The meridional integral of the geopotential anomaly along the
               ! interface below a layer, divided by the grid spacing, m2 s-2.
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)) :: &
    inty_dza    ! The change in inty_za through a layer, in m2 s-2.
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density, in Pa (usually 2e7 Pa = 2000 dbar).

  real :: dp_neglect         ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected, in Pa.
  real :: alpha_anom         ! The in-situ specific volume, averaged over a
                             ! layer, less alpha_ref, in m3 kg-1.
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.

  real :: alpha_ref ! A reference specific volume, in m3 kg-1, that is used
                    ! to reduce the impact of truncation errors.
  real :: rho_in_situ(SZI_(G)) ! The in situ density, in kg m-3.
  real :: Pa_to_H   ! A factor to convert from Pa to the thicknesss units (H).
!  real :: oneatm = 101325.0  ! 1 atm in Pa (kg/ms2)
  real :: I_gEarth
  real, parameter :: C1_6 = 1.0/6.0
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: is_bk, ie_bk, js_bk, je_bk, Isq_bk, Ieq_bk, Jsq_bk, Jeq_bk
  integer :: i, j, k, n, ib, jb, ioff_bk, joff_bk

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  nkmb=GV%nk_rho_varies
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  use_p_atm = .false.
  if (present(p_atm)) then ; if (associated(p_atm)) use_p_atm = .true. ; endif
  use_EOS = associated(tv%eqn_of_state)

  if (.not.associated(CS)) call MOM_error(FATAL, &
       "MOM_PressureForce: Module must be initialized before it is used.")

  dp_neglect = GV%H_to_Pa * GV%H_subroundoff
  alpha_ref = 1.0/CS%Rho0

!$OMP parallel default(none) shared(Isq,Ieq,Jsq,Jeq,nz,use_p_atm,p,p_atm,GV,h)
  if (use_p_atm) then
!$OMP do
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      p(i,j,1) = p_atm(i,j)
    enddo ; enddo
  else
!$OMP do
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      p(i,j,1) = 0.0 ! or oneatm
    enddo ; enddo
  endif
!$OMP do
  do j=Jsq,Jeq+1 ; do k=2,nz+1 ; do i=Isq,Ieq+1
    p(i,j,K) = p(i,j,K-1) + GV%H_to_Pa * h(i,j,k-1)
  enddo ; enddo ; enddo
!$OMP end parallel

  I_gEarth = 1.0 / GV%g_Earth

  if (use_EOS) then
  !   With a bulk mixed layer, replace the T & S of any layers that are
  ! lighter than the the buffer layer with the properties of the buffer
  ! layer.  These layers will be massless anyway, and it avoids any
  ! formal calculations with hydrostatically unstable profiles.
    if (nkmb>0) then
      tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
      tv_tmp%eqn_of_state => tv%eqn_of_state
      do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,nz,nkmb,tv_tmp,tv,p_ref,GV) &
!$OMP                          private(Rho_cv_BL)
      do j=Jsq,Jeq+1
        do k=1,nkmb ; do i=Isq,Ieq+1
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        enddo ; enddo
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, &
                        Rho_cv_BL(:), Isq, Ieq-Isq+2, tv%eqn_of_state)
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

!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,nz,is,ie,js,je,tv_tmp,alpha_ref, &
!$OMP                                  p,h,G,GV,tv,dza,intp_dza,intx_dza,inty_dza,use_EOS) &
!$OMP                          private(alpha_anom,dp)
  do k=1,nz
    ! Calculate 4 integrals through the layer that are required in the
    ! subsequent calculation.
    if (use_EOS) then
      call int_specific_vol_dp(tv_tmp%T(:,:,k), tv_tmp%S(:,:,k), p(:,:,K), &
                               p(:,:,K+1), alpha_ref, G%HI, tv%eqn_of_state, &
                               dza(:,:,k), intp_dza(:,:,k), intx_dza(:,:,k), &
                               inty_dza(:,:,k))
    else
      alpha_anom = 1.0/GV%Rlay(k) - alpha_ref
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
!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,nz,za,alpha_ref,p,G,GV,dza)
  do j=Jsq,Jeq+1
    do i=Isq,Ieq+1
      za(i,j) = alpha_ref*p(i,j,nz+1) - GV%g_Earth*G%bathyT(i,j)
    enddo
    do k=nz,1,-1 ; do i=Isq,Ieq+1
    za(i,j) = za(i,j) + dza(i,j,k)
    enddo ; enddo
  enddo

  if (CS%tides) then
    ! Find and add the tidal geopotential anomaly.
!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,SSH,za,alpha_ref,p,I_gEarth)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      SSH(i,j) = (za(i,j) - alpha_ref*p(i,j,1)) * I_gEarth
    enddo ; enddo
    call calc_tidal_forcing(CS%Time, SSH, e_tidal, G, CS%tides_CSp)
!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,za,G,GV,e_tidal)
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      za(i,j) = za(i,j) - GV%g_Earth*e_tidal(i,j)
    enddo ; enddo
  endif

  if (CS%GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    if (use_EOS) then
!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,tv_tmp,p,tv,dM,CS,alpha_ref,za) &
!$OMP                          private(rho_in_situ)
      do j=Jsq,Jeq+1
        call calculate_density(tv_tmp%T(:,j,1), tv_tmp%S(:,j,1), p(:,j,1), &
                               rho_in_situ, Isq, Ieq-Isq+2, tv%eqn_of_state)

        do i=Isq,Ieq+1
          dM(i,j) = (CS%GFS_scale - 1.0) * &
            (p(i,j,1)*(1.0/rho_in_situ(i) - alpha_ref) + za(i,j))
        enddo
      enddo
    else
!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,dM,CS,p,GV,alpha_ref,za)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dM(i,j) = (CS%GFS_scale - 1.0) * &
          (p(i,j,1)*(1.0/GV%Rlay(1) - alpha_ref) + za(i,j))
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
!$OMP parallel do default(none) shared(nz,za,G,GV,dza,intx_dza,h,PFu, &
!$OMP                                  intp_dza,p,dp_neglect,inty_dza,PFv,CS,dM) &
!$OMP                          private(is_bk,ie_bk,js_bk,je_bk,Isq_bk,Ieq_bk,Jsq_bk, &
!$OMP                                  Jeq_bk,ioff_bk,joff_bk,i,j,za_bk,intx_za_bk,  &
!$OMP                                  inty_za_bk,dp_bk)
  do n = 1, G%nblocks
    is_bk=G%block(n)%isc      ; ie_bk=G%block(n)%iec
    js_bk=G%block(n)%jsc      ; je_bk=G%block(n)%jec
    Isq_bk=G%block(n)%IscB    ; Ieq_bk=G%block(n)%IecB
    Jsq_bk=G%block(n)%JscB    ; Jeq_bk=G%block(n)%JecB
    ioff_bk = G%Block(n)%idg_offset - G%HI%idg_offset
    joff_bk = G%Block(n)%jdg_offset - G%HI%jdg_offset
    do jb=Jsq_bk,Jeq_bk+1 ; do ib=Isq_bk,Ieq_bk+1
      i = ib+ioff_bk ; j = jb+joff_bk
      za_bk(ib,jb) = za(i,j)
    enddo ; enddo
    do jb=js_bk,je_bk ; do Ib=Isq_bk,Ieq_bk
      I = Ib+ioff_bk ; j = jb+joff_bk
      intx_za_bk(Ib,jb) = 0.5*(za_bk(ib,jb) + za_bk(ib+1,jb))
    enddo ; enddo
    do Jb=Jsq_bk,Jeq_bk ; do ib=is_bk,ie_bk
      i = ib+ioff_bk ; J = Jb+joff_bk
      inty_za_bk(ib,Jb) = 0.5*(za_bk(ib,jb) + za_bk(ib,jb+1))
    enddo ; enddo
    do k=1,nz
      ! These expressions for the acceleration have been carefully checked in
      ! a set of idealized cases, and should be bug-free.
      do jb=Jsq_bk,Jeq_bk+1 ; do ib=Isq_bk,Ieq_bk+1
        i = ib+ioff_bk ; j = jb+joff_bk
        dp_bk(ib,jb) = GV%H_to_Pa*h(i,j,k)
        za_bk(ib,jb) = za_bk(ib,jb) - dza(i,j,k)
      enddo ; enddo
      do jb=js_bk,je_bk ; do Ib=Isq_bk,Ieq_bk
        I = Ib+ioff_bk ; j = jb+joff_bk
        intx_za_bk(Ib,jb) = intx_za_bk(Ib,jb) - intx_dza(I,j,k)
        PFu(I,j,k) = (((za_bk(ib,jb)*dp_bk(ib,jb) + intp_dza(i,j,k)) - &
                     (za_bk(ib+1,jb)*dp_bk(ib+1,jb) + intp_dza(i+1,j,k))) + &
                     ((dp_bk(ib+1,jb) - dp_bk(ib,jb)) * intx_za_bk(Ib,jb) - &
                     (p(i+1,j,K) - p(i,j,K)) * intx_dza(I,j,k))) * &
                     (2.0*G%IdxCu(I,j) / ((dp_bk(ib,jb) + dp_bk(ib+1,jb)) + &
                     dp_neglect))
      enddo ; enddo
      do Jb=Jsq_bk,Jeq_bk ; do ib=is_bk,ie_bk
        i = ib+ioff_bk ; J = Jb+joff_bk
        inty_za_bk(ib,Jb) = inty_za_bk(ib,Jb) - inty_dza(i,J,k)
        PFv(i,J,k) = (((za_bk(ib,jb)*dp_bk(ib,jb) + intp_dza(i,j,k)) - &
                     (za_bk(ib,jb+1)*dp_bk(ib,jb+1) + intp_dza(i,j+1,k))) + &
                     ((dp_bk(ib,jb+1) - dp_bk(ib,jb)) * inty_za_bk(ib,Jb) - &
                     (p(i,j+1,K) - p(i,j,K)) * inty_dza(i,J,k))) * &
                     (2.0*G%IdyCv(i,J) / ((dp_bk(ib,jb) + dp_bk(ib,jb+1)) + &
                     dp_neglect))
      enddo ; enddo

      if (CS%GFS_scale < 1.0) then
        ! Adjust the Montgomery potential to make this a reduced gravity model.
        do j=js_bk+joff_bk,je_bk+joff_bk ; do I=Isq_bk+ioff_bk,Ieq_bk+ioff_bk
          PFu(I,j,k) = PFu(I,j,k) - (dM(i+1,j) - dM(i,j)) * G%IdxCu(I,j)
        enddo ; enddo
        do J=Jsq_bk+joff_bk,Jeq_bk+joff_bk ; do i=is_bk+ioff_bk,ie_bk+ioff_bk
          PFv(i,J,k) = PFv(i,J,k) - (dM(i,j+1) - dM(i,j)) * G%IdyCv(i,J)
        enddo ; enddo
      endif
    enddo
  enddo

  if (present(pbce)) then
    call set_pbce_nonBouss(p, tv_tmp, G, GV, GV%g_Earth, CS%GFS_scale, pbce)
  endif

  if (present(eta)) then
    Pa_to_H = 1.0 / GV%H_to_Pa
    if (use_p_atm) then
!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,nz,eta,p,p_atm,Pa_to_H)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = (p(i,j,nz+1) - p_atm(i,j))*Pa_to_H ! eta has the same units as h.
      enddo ; enddo
    else
!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,nz,eta,p,Pa_to_H)
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
!! To work, the following fields must be set outside of the usual
!! ie to ie, je to je range before this subroutine is called:
!!  h[ie+1] and h[je+1] and (if tv%eqn_of_state is set) T[ie+1], S[ie+1],
!!  T[je+1], and S[je+1].
subroutine PressureForce_AFV_Bouss(h, tv, PFu, PFv, G, GV, CS, ALE_CSp, p_atm, pbce, eta)
  type(ocean_grid_type),                     intent(in)  :: G   !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV  !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)  :: h   !< Layer thickness (kg/m2)
  type(thermo_var_ptrs),                     intent(in)  :: tv  !< Thermodynamic variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(out) :: PFu !< Zonal acceleration (m/s2)
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(out) :: PFv !< Meridional acceleration (m/s2)
  type(PressureForce_AFV_CS),                pointer     :: CS  !< Finite volume PGF control structure
  type(ALE_CS),                              pointer     :: ALE_CSp !< ALE control structure
  real, dimension(:,:),                      optional, pointer :: p_atm !< The pressure at the ice-ocean
                                                         !! or atmosphere-ocean interface in Pa.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  optional, intent(out) :: pbce !< The baroclinic pressure
                                                         !! anomaly in each layer due to eta anomalies,
                                                         !! in m2 s-2 H-1.
  real, dimension(SZI_(G),SZJ_(G)),          optional, intent(out) :: eta !< The bottom mass used to
                                                         !! calculate PFu and PFv, in H, with any tidal
                                                         !! contributions or compressibility compensation.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: e ! Interface height in m.
  real, dimension(SZI_(G),SZJ_(G))  :: &
    e_tidal, &  ! The bottom geopotential anomaly due to tidal forces from
                ! astronomical sources and self-attraction and loading, in m.
    dM          ! The barotropic adjustment to the Montgomery potential to
                ! account for a reduced gravity model, in m2 s-2.
  real, dimension(SZI_(G)) :: &
    Rho_cv_BL   !   The coordinate potential density in the deepest variable
                ! density near-surface layer, in kg m-3.
  real, dimension(SZDI_(G%Block(1)),SZDJ_(G%Block(1))) :: &  ! on block indices
    dz_bk, &     ! The change in geopotential thickness through a layer, m2 s-2.
    pa_bk, &     ! The pressure anomaly (i.e. pressure + g*RHO_0*e) at the
                 ! the interface atop a layer, in Pa.
    dpa_bk, &    ! The change in pressure anomaly between the top and bottom
                 ! of a layer, in Pa.
    intz_dpa_bk  ! The vertical integral in depth of the pressure anomaly less
                 ! the pressure anomaly at the top of the layer, in H Pa (m Pa).
  real, dimension(SZDIB_(G%Block(1)),SZDJ_(G%Block(1))) :: & ! on block indices
    intx_pa_bk, & ! The zonal integral of the pressure anomaly along the interface
                  ! atop a layer, divided by the grid spacing, in Pa.
    intx_dpa_bk   ! The change in intx_pa through a layer, in Pa.
  real, dimension(SZDI_(G%Block(1)),SZDJB_(G%Block(1))) :: & ! on block indices
    inty_pa_bk, & ! The meridional integral of the pressure anomaly along the
                  ! interface atop a layer, divided by the grid spacing, in Pa.
    inty_dpa_bk   ! The change in inty_pa through a layer, in Pa.

  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), target :: &
    T_tmp, &    ! Temporary array of temperatures where layers that are lighter
                ! than the mixed layer have the mixed layer's properties, in C.
    S_tmp       ! Temporary array of salinities where layers that are lighter
                ! than the mixed layer have the mixed layer's properties, in psu.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    S_t, S_b, T_t, T_b ! Top and bottom edge values for linear reconstructions
                       ! of salinity and temperature within each layer.
  real :: rho_in_situ(SZI_(G)) ! The in situ density, in kg m-3.
  real :: p_ref(SZI_(G))     !   The pressure used to calculate the coordinate
                             ! density, in Pa (usually 2e7 Pa = 2000 dbar).
  real :: p0(SZI_(G)) ! An array of zeros to use for pressure in Pa.
  real :: h_neglect          ! A thickness that is so small it is usually lost
                             ! in roundoff and can be neglected, in m.
  real :: I_Rho0             ! 1/Rho0.
  real :: G_Rho0             ! G_Earth / Rho0 in m4 s-2 kg-1.
  real :: Rho_ref            ! The reference density in kg m-3.
  logical :: use_p_atm       ! If true, use the atmospheric pressure.
  logical :: use_ALE         ! If true, use an ALE pressure reconstruction.
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  type(thermo_var_ptrs) :: tv_tmp! A structure of temporary T & S.

  real, parameter :: C1_6 = 1.0/6.0
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz, nkmb
  integer :: is_bk, ie_bk, js_bk, je_bk, Isq_bk, Ieq_bk, Jsq_bk, Jeq_bk
  integer :: ioff_bk, joff_bk
  integer :: i, j, k, n, ib, jb
  integer :: PRScheme

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  nkmb=GV%nk_rho_varies
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.associated(CS)) call MOM_error(FATAL, &
       "MOM_PressureForce: Module must be initialized before it is used.")

  use_p_atm = .false.
  if (present(p_atm)) then ; if (associated(p_atm)) use_p_atm = .true. ; endif
  use_EOS = associated(tv%eqn_of_state)
  do i=Isq,Ieq+1 ; p0(i) = 0.0 ; enddo
  use_ALE = .false.
  if (associated(ALE_CSp)) use_ALE = usePressureReconstruction(ALE_CSp) .and. use_EOS

  PRScheme = pressureReconstructionScheme(ALE_CSp)
  h_neglect = GV%H_subroundoff
  I_Rho0 = 1.0/GV%Rho0
  G_Rho0 = GV%g_Earth/GV%Rho0
  rho_ref = CS%Rho0

  if (CS%tides) then
    !   Determine the surface height anomaly for calculating self attraction
    ! and loading.  This should really be based on bottom pressure anomalies,
    ! but that is not yet implemented, and the current form is correct for
    ! barotropic tides.
!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,nz,e,G,GV,h)
    do j=Jsq,Jeq+1
      do i=Isq,Ieq+1
        e(i,j,1) = -1.0*G%bathyT(i,j)
      enddo
      do k=1,nz ; do i=Isq,Ieq+1
        e(i,j,1) = e(i,j,1) + h(i,j,k)*GV%H_to_m
      enddo ; enddo
    enddo
    call calc_tidal_forcing(CS%Time, e(:,:,1), e_tidal, G, CS%tides_CSp)
  endif

!    Here layer interface heights, e, are calculated.
!$OMP parallel default(none) shared(Isq,Ieq,Jsq,Jeq,nz,e,G,GV,h,CS,e_tidal)
  if (CS%tides) then
!$OMP do
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,nz+1) = -1.0*G%bathyT(i,j) - e_tidal(i,j)
    enddo ; enddo
  else
!$OM do
    do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
      e(i,j,nz+1) = -1.0*G%bathyT(i,j)
    enddo ; enddo
  endif
!$OMP do
  do j=Jsq,Jeq+1; do k=nz,1,-1 ; do i=Isq,Ieq+1
    e(i,j,K) = e(i,j,K+1) + h(i,j,k)*GV%H_to_m
  enddo ; enddo ; enddo
!$OMP end parallel


  if (use_EOS) then
! With a bulk mixed layer, replace the T & S of any layers that are
! lighter than the the buffer layer with the properties of the buffer
! layer.  These layers will be massless anyway, and it avoids any
! formal calculations with hydrostatically unstable profiles.

    if (nkmb>0) then
      tv_tmp%T => T_tmp ; tv_tmp%S => S_tmp
      tv_tmp%eqn_of_state => tv%eqn_of_state

      do i=Isq,Ieq+1 ; p_ref(i) = tv%P_Ref ; enddo
!$OMP parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,nkmb,nz,GV,tv_tmp,tv,p_ref) &
!$OMP                          private(Rho_cv_BL)
      do j=Jsq,Jeq+1
        do k=1,nkmb ; do i=Isq,Ieq+1
          tv_tmp%T(i,j,k) = tv%T(i,j,k) ; tv_tmp%S(i,j,k) = tv%S(i,j,k)
        enddo ; enddo
        call calculate_density(tv%T(:,j,nkmb), tv%S(:,j,nkmb), p_ref, &
                        Rho_cv_BL(:), Isq, Ieq-Isq+2, tv%eqn_of_state)

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

!$OMP parallel default(none) shared(Jsq,Jeq,Isq,Ieq,tv_tmp,p_atm,rho_in_situ,tv, &
!$OMP                               p0,dM,CS,G_Rho0,e,use_p_atm,use_EOS,GV,    &
!$OMP                               rho_ref,js,je,is,ie)
  if (CS%GFS_scale < 1.0) then
    ! Adjust the Montgomery potential to make this a reduced gravity model.
    if (use_EOS) then
!$OMP do
      do j=Jsq,Jeq+1
        if (use_p_atm) then
          call calculate_density(tv_tmp%T(:,j,1), tv_tmp%S(:,j,1), p_atm(:,j), &
                                 rho_in_situ, Isq, Ieq-Isq+2, tv%eqn_of_state)
        else
          call calculate_density(tv_tmp%T(:,j,1), tv_tmp%S(:,j,1), p0, &
                                 rho_in_situ, Isq, Ieq-Isq+2, tv%eqn_of_state)
        endif
        do i=Isq,Ieq+1
          dM(i,j) = (CS%GFS_scale - 1.0) * (G_Rho0 * rho_in_situ(i)) * e(i,j,1)
        enddo
      enddo
    else
!$OMP do
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        dM(i,j) = (CS%GFS_scale - 1.0) * (G_Rho0 * GV%Rlay(1)) * e(i,j,1)
      enddo ; enddo
    endif
  endif
!$OMP end parallel

! Have checked that rho_0 drops out and that the 1-layer case is right. RWH.

  ! If regridding is activated, do a linear reconstruction of salinity
  ! and temperature across each layer. The subscripts 't' and 'b' refer
  ! to top and bottom values within each layer (these are the only degrees
  ! of freedeom needed to know the linear profile).
  if ( use_ALE ) then
    if ( PRScheme == PRESSURE_RECONSTRUCTION_PLM ) then
      call pressure_gradient_plm (ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h);
    elseif ( PRScheme == PRESSURE_RECONSTRUCTION_PPM ) then
      call pressure_gradient_ppm (ALE_CSp, S_t, S_b, T_t, T_b, G, GV, tv, h);
    endif
  endif

!$OMP parallel do default(none) shared(use_p_atm,rho_ref,G,GV,e,     &
!$OMP                                  p_atm,nz,use_EOS,use_ALE,PRScheme,T_t,T_b,S_t, &
!$OMP                                  S_b,CS,tv,tv_tmp,h,PFu,I_Rho0,h_neglect,PFv,dM)&
!$OMP                          private(is_bk,ie_bk,js_bk,je_bk,Isq_bk,Ieq_bk,Jsq_bk,  &
!$OMP                                  Jeq_bk,ioff_bk,joff_bk,pa_bk,  &
!$OMP                                  intx_pa_bk,inty_pa_bk,dpa_bk,intz_dpa_bk,      &
!$OMP                                  intx_dpa_bk,inty_dpa_bk,dz_bk,i,j)
  do n = 1, G%nblocks
    is_bk=G%Block(n)%isc      ; ie_bk=G%Block(n)%iec
    js_bk=G%Block(n)%jsc      ; je_bk=G%Block(n)%jec
    Isq_bk=G%Block(n)%IscB    ; Ieq_bk=G%Block(n)%IecB
    Jsq_bk=G%Block(n)%JscB    ; Jeq_bk=G%Block(n)%JecB
    ioff_bk = G%Block(n)%idg_offset - G%HI%idg_offset
    joff_bk = G%Block(n)%jdg_offset - G%HI%jdg_offset

    ! Set the surface boundary conditions on pressure anomaly and its horizontal
    ! integrals, assuming that the surface pressure anomaly varies linearly
    ! in x and y.
    if (use_p_atm) then
      do jb=Jsq_bk,Jeq_bk+1 ; do ib=Isq_bk,Ieq_bk+1
        i = ib+ioff_bk ; j = jb+joff_bk
        pa_bk(ib,jb) = (rho_ref*GV%g_Earth)*e(i,j,1) + p_atm(i,j)
      enddo ; enddo
    else
      do jb=Jsq_bk,Jeq_bk+1 ; do ib=Isq_bk,Ieq_bk+1
        i = ib+ioff_bk ; j = jb+joff_bk
        pa_bk(ib,jb) = (rho_ref*GV%g_Earth)*e(i,j,1)
      enddo ; enddo
    endif
    do jb=js_bk,je_bk ; do Ib=Isq_bk,Ieq_bk
      intx_pa_bk(Ib,jb) = 0.5*(pa_bk(ib,jb) + pa_bk(ib+1,jb))
    enddo ; enddo
    do Jb=Jsq_bk,Jeq_bk ; do ib=is_bk,ie_bk
      inty_pa_bk(ib,Jb) = 0.5*(pa_bk(ib,jb) + pa_bk(ib,jb+1))
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
          if ( PRScheme == PRESSURE_RECONSTRUCTION_PLM ) then
            call int_density_dz_generic_plm ( T_t(:,:,k), T_b(:,:,k), &
                      S_t(:,:,k), S_b(:,:,k), e(:,:,K), e(:,:,K+1), &
                      rho_ref, CS%Rho0, GV%g_Earth,    &
                      GV%H_subroundoff, G%bathyT, G%HI, G%Block(n), &
                      tv%eqn_of_state, dpa_bk, intz_dpa_bk, intx_dpa_bk, inty_dpa_bk, &
                      useMassWghtInterp = CS%useMassWghtInterp)
          elseif ( PRScheme == PRESSURE_RECONSTRUCTION_PPM ) then
            call int_density_dz_generic_ppm ( tv%T(:,:,k), T_t(:,:,k), T_b(:,:,k), &
                      tv%S(:,:,k), S_t(:,:,k), S_b(:,:,k), e(:,:,K), e(:,:,K+1), &
                      rho_ref, CS%Rho0, GV%g_Earth, &
                      G%HI, G%Block(n), tv%eqn_of_state, dpa_bk, intz_dpa_bk,    &
                      intx_dpa_bk, inty_dpa_bk)
          endif
        else
          call int_density_dz(tv_tmp%T(:,:,k), tv_tmp%S(:,:,k), &
                    e(:,:,K), e(:,:,K+1),             &
                    rho_ref, CS%Rho0, GV%g_Earth, G%HI, G%Block(n), tv%eqn_of_state, &
                    dpa_bk, intz_dpa_bk, intx_dpa_bk, inty_dpa_bk )
        endif
        intz_dpa_bk(:,:) = intz_dpa_bk(:,:)*GV%m_to_H
      else
        do jb=Jsq_bk,Jeq_bk+1 ; do ib=Isq_bk,Ieq_bk+1
          i = ib+ioff_bk ; j = jb+joff_bk
          dz_bk(ib,jb) = GV%g_Earth*GV%H_to_m*h(i,j,k)
          dpa_bk(ib,jb) = (GV%Rlay(k) - rho_ref)*dz_bk(ib,jb)
          intz_dpa_bk(ib,jb) = 0.5*(GV%Rlay(k) - rho_ref)*dz_bk(ib,jb)*h(i,j,k)
        enddo ; enddo
        do jb=js_bk,je_bk ; do Ib=Isq_bk,Ieq_bk
          intx_dpa_bk(Ib,jb) = 0.5*(GV%Rlay(k) - rho_ref) * (dz_bk(ib,jb)+dz_bk(ib+1,jb))
        enddo ; enddo
        do Jb=Jsq_bk,Jeq_bk ; do ib=is_bk,ie_bk
          inty_dpa_bk(ib,Jb) = 0.5*(GV%Rlay(k) - rho_ref) * (dz_bk(ib,jb)+dz_bk(ib,jb+1))
        enddo ; enddo
      endif

      ! Compute pressure gradient in x direction
      do jb=js_bk,je_bk ; do Ib=Isq_bk,Ieq_bk
        I = Ib+ioff_bk ; j = jb+joff_bk
        PFu(I,j,k) = (((pa_bk(ib,jb)*h(i,j,k) + intz_dpa_bk(ib,jb)) - &
                     (pa_bk(ib+1,jb)*h(i+1,j,k) + intz_dpa_bk(ib+1,jb))) + &
                     ((h(i+1,j,k) - h(i,j,k)) * intx_pa_bk(Ib,jb) - &
                     (e(i+1,j,K+1) - e(i,j,K+1)) * intx_dpa_bk(Ib,jb) * GV%m_to_H)) * &
                     ((2.0*I_Rho0*G%IdxCu(I,j)) / &
                     ((h(i,j,k) + h(i+1,j,k)) + h_neglect))
        intx_pa_bk(Ib,jb) = intx_pa_bk(Ib,jb) + intx_dpa_bk(Ib,jb)
      enddo ; enddo
      ! Compute pressure gradient in y direction
      do Jb=Jsq_bk,Jeq_bk ; do ib=is_bk,ie_bk
        i = ib+ioff_bk ; J = Jb+joff_bk
        PFv(i,J,k) = (((pa_bk(ib,jb)*h(i,j,k) + intz_dpa_bk(ib,jb)) - &
                     (pa_bk(ib,jb+1)*h(i,j+1,k) + intz_dpa_bk(ib,jb+1))) + &
                     ((h(i,j+1,k) - h(i,j,k)) * inty_pa_bk(ib,Jb) - &
                     (e(i,j+1,K+1) - e(i,j,K+1)) * inty_dpa_bk(ib,Jb) * GV%m_to_H)) * &
                     ((2.0*I_Rho0*G%IdyCv(i,J)) / &
                     ((h(i,j,k) + h(i,j+1,k)) + h_neglect))
        inty_pa_bk(ib,Jb) = inty_pa_bk(ib,Jb) + inty_dpa_bk(ib,Jb)
      enddo ; enddo
      do jb=Jsq_bk,Jeq_bk+1 ; do ib=Isq_bk,Ieq_bk+1
        pa_bk(ib,jb) = pa_bk(ib,jb) + dpa_bk(ib,jb)
      enddo ; enddo
    enddo

    if (CS%GFS_scale < 1.0) then
      do k=1,nz
        do j=js_bk+joff_bk,je_bk+joff_bk ; do I=Isq_bk+ioff_bk,Ieq_bk+ioff_bk
          PFu(I,j,k) = PFu(I,j,k) - (dM(i+1,j) - dM(i,j)) * G%IdxCu(I,j)
        enddo ; enddo
        do J=Jsq_bk+joff_bk,Jeq_bk+joff_bk ; do i=is_bk+ioff_bk,ie_bk+ioff_bk
          PFv(i,J,k) = PFv(i,J,k) - (dM(i,j+1) - dM(i,j)) * G%IdyCv(i,J)
        enddo ; enddo
      enddo
    endif
  enddo

  if (present(pbce)) then
    call set_pbce_Bouss(e, tv_tmp, G, GV, GV%g_Earth, CS%Rho0, CS%GFS_scale, pbce)
  endif

  if (present(eta)) then
    if (CS%tides) then
    ! eta is the sea surface height relative to a time-invariant geoid, for
    ! comparison with what is used for eta in btstep.  See how e was calculated
    ! about 200 lines above.
!$OM parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,eta,e,e_tidal)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = e(i,j,1)*GV%m_to_H + e_tidal(i,j)*GV%m_to_H
      enddo ; enddo
    else
!$OM parallel do default(none) shared(Isq,Ieq,Jsq,Jeq,eta,e)
      do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
        eta(i,j) = e(i,j,1)*GV%m_to_H
      enddo ; enddo
    endif
  endif

  if (CS%id_e_tidal>0) call post_data(CS%id_e_tidal, e_tidal, CS%diag)

end subroutine PressureForce_AFV_Bouss

!> Initializes the finite volume pressure gradient control structure
subroutine PressureForce_AFV_init(Time, G, GV, param_file, diag, CS, tides_CSp)
  type(time_type), target,    intent(in)    :: Time !< Current model time
  type(ocean_grid_type),      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV !< Vertical grid structure
  type(param_file_type),      intent(in)    :: param_file !< Parameter file handles
  type(diag_ctrl), target,    intent(inout) :: diag !< Diagnostics control structure
  type(PressureForce_AFV_CS), pointer       :: CS !< Finite volume PGF control structure
  type(tidal_forcing_CS), optional, pointer :: tides_CSp !< Tides control structure
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "PressureForce_init called with an associated "// &
                            "control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag ; CS%Time => Time
  if (present(tides_CSp)) then
    if (associated(tides_CSp)) CS%tides_CSp => tides_CSp
  endif

  mod = "MOM_PressureForce_AFV"
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "RHO_0", CS%Rho0, &
                 "The mean ocean density used with BOUSSINESQ true to \n"//&
                 "calculate accelerations and the mass for conservation \n"//&
                 "properties, or with BOUSSINSEQ false to convert some \n"//&
                 "parameters from vertical units of m to kg m-2.", &
                 units="kg m-3", default=1035.0)
  call get_param(param_file, mod, "TIDES", CS%tides, &
                 "If true, apply tidal momentum forcing.", default=.false.)
  call get_param(param_file, mod, "MASS_WEIGHT_IN_PRESSURE_GRADIENT", CS%useMassWghtInterp, &
                 "If true, use mass weighting when interpolation T/S for\n"//&
                 "top/bottom integrals in AFV pressure gradient calculation.", default=.false.)

  if (CS%tides) then
    CS%id_e_tidal = register_diag_field('ocean_model', 'e_tidal', diag%axesT1, &
        Time, 'Tidal Forcing Astronomical and SAL Height Anomaly', 'meter')
  endif

  CS%GFS_scale = 1.0
  if (GV%g_prime(1) /= GV%g_Earth) CS%GFS_scale = GV%g_prime(1) / GV%g_Earth

  call log_param(param_file, mod, "GFS / G_EARTH", CS%GFS_scale)

end subroutine PressureForce_AFV_init

!> Deallocates the finite volume pressure gradient control structure
subroutine PressureForce_AFV_end(CS)
  type(PressureForce_AFV_CS), pointer :: CS
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
