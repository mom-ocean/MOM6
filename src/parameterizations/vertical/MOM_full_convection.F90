!> Does full convective adjustment of unstable regions via a strong diffusivity.
module MOM_full_convection

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density_derivs, EOS_domain

implicit none ; private

#include <MOM_memory.h>

public full_convection

contains

!> Calculate new temperatures and salinities that have been subject to full convective mixing.
subroutine full_convection(G, GV, US, h, tv, T_adj, S_adj, p_surf, Kddt_smooth, &
                           Kddt_convect, halo)
  type(ocean_grid_type),   intent(in)    :: G     !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)    :: tv    !< A structure pointing to various
                                                  !! thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(out)   :: T_adj !< Adjusted potential temperature [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(out)   :: S_adj !< Adjusted salinity [ppt].
  real, dimension(:,:),    pointer       :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa] (or NULL).
  real,                    intent(in)    :: Kddt_smooth  !< A smoothing vertical
                                                  !! diffusivity times a timestep [H2 ~> m2 or kg2 m-4].
  real,          optional, intent(in)    :: Kddt_convect !< A large convecting vertical
                                                  !! diffusivity times a timestep [H2 ~> m2 or kg2 m-4].
  integer,       optional, intent(in)    :: halo  !< Halo width over which to compute

  ! Local variables
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    dRho_dT, &  ! The derivative of density with temperature [R degC-1 ~> kg m-3 degC-1]
    dRho_dS     ! The derivative of density with salinity [R ppt-1 ~> kg m-3 ppt-1].
  real :: h_neglect, h0 ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
! logical :: use_EOS    ! If true, density is calculated from T & S using an equation of state.
  real, dimension(SZI_(G),SZK0_(G)) :: &
    Te_a, & ! A partially updated temperature estimate including the influnce from
            ! mixing with layers above rescaled by a factor of d_a [degC].
            ! This array is discreted on tracer cells, but contains an extra
            ! layer at the top for algorithmic convenience.
    Se_a    ! A partially updated salinity estimate including the influnce from
            ! mixing with layers above rescaled by a factor of d_a [ppt].
            ! This array is discreted on tracer cells, but contains an extra
            ! layer at the top for algorithmic convenience.
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    Te_b, & ! A partially updated temperature estimate including the influnce from
            ! mixing with layers below rescaled by a factor of d_b [degC].
            ! This array is discreted on tracer cells, but contains an extra
            ! layer at the bottom for algorithmic convenience.
    Se_b    ! A partially updated salinity estimate including the influnce from
            ! mixing with layers below rescaled by a factor of d_b [ppt].
            ! This array is discreted on tracer cells, but contains an extra
            ! layer at the bottom for algorithmic convenience.
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    c_a, &  ! The fractional influence of the properties of the layer below
            ! in the final properies with a downward-first solver, nondim.
    d_a, &  ! The fractional influence of the properties of the layer in question
            ! and layers above in the final properies with a downward-first solver, nondim.
            ! d_a = 1.0 - c_a
    c_b, &  ! The fractional influence of the properties of the layer above
            ! in the final properies with a upward-first solver, nondim.
    d_b     ! The fractional influence of the properties of the layer in question
            ! and layers below in the final properies with a upward-first solver, nondim.
            ! d_b = 1.0 - c_b
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    mix     !< The amount of mixing across the interface between layers [H ~> m or kg m-2].
  real :: mix_len  ! The length-scale of mixing, when it is active [H ~> m or kg m-2]
  real :: h_b, h_a ! The thicknessses of the layers above and below an interface [H ~> m or kg m-2]
  real :: b_b, b_a ! Inverse pivots used by the tridiagonal solver [H-1 ~> m-1 or m2 kg-1].

  real :: kap_dt_x2 ! The product of 2*kappa*dt [H2 ~> m2 or kg2 m-4].

  logical, dimension(SZI_(G)) :: do_i ! Do more work on this column.
  logical, dimension(SZI_(G)) :: last_down ! The last setup pass was downward.
  integer, dimension(SZI_(G)) :: change_ct ! The number of interfaces where the
                         ! mixing has changed this iteration.
  integer :: changed_col ! The number of colums whose mixing changed.
  integer :: i, j, k, is, ie, js, je, nz, itt

  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  else
    is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  endif
  nz = G%ke

  if (.not.associated(tv%eqn_of_state)) return

  h_neglect = GV%H_subroundoff
  kap_dt_x2 = 0.0
  if (present(Kddt_convect)) kap_dt_x2 = 2.0*Kddt_convect
  mix_len = (1.0e20 * nz) * (G%max_depth * GV%Z_to_H)
  h0 = 1.0e-16*sqrt(Kddt_smooth) + h_neglect

  do j=js,je
    mix(:,:) = 0.0 ; d_b(:,:) = 1.0
    ! These would be Te_b(:,:) = tv%T(:,j,:), etc., but the values are not used
    Te_b(:,:) = 0.0 ; Se_b(:,:) = 0.0

    call smoothed_dRdT_dRdS(h, tv, Kddt_smooth, dRho_dT, dRho_dS, G, GV, US, j, p_surf, halo)

    do i=is,ie
      do_i(i) = (G%mask2dT(i,j) > 0.0)

      d_a(i,1) = 1.0
      last_down(i) = .true. ! This is set for debuggers.
      ! These are extra values are used for convenience in the stability test
      Te_a(i,0) = 0.0 ; Se_a(i,0) = 0.0
    enddo

    do itt=1,nz ! At least 2 interfaces will change with each full pass, or the
                ! iterations stop, so the maximum count of nz is very conservative.

      do i=is,ie ; change_ct(i) = 0 ; enddo
      ! Move down the water column, finding unstable interfaces, and building up the
      ! temporary arrays for the tridiagonal solver.
      do K=2,nz ; do i=is,ie ; if (do_i(i)) then

        h_a = h(i,j,k-1) + h_neglect ; h_b = h(i,j,k) + h_neglect
        if (mix(i,K) <= 0.0) then
          if (is_unstable(dRho_dT(i,K), dRho_dS(i,K), h_a, h_b, mix(i,K-1), mix(i,K+1), &
                          tv%T(i,j,k-1), tv%T(i,j,k), tv%S(i,j,k-1), tv%S(i,j,k), &
                          Te_a(i,k-2), Te_b(i,k+1), Se_a(i,k-2), Se_b(i,k+1), &
                          d_a(i,K-1), d_b(i,K+1))) then
            mix(i,K) = mix_len
            if (kap_dt_x2 > 0.0) mix(i,K) = kap_dt_x2 / ((h(i,j,k-1)+h(i,j,k)) + h0)
            change_ct(i) = change_ct(i) + 1
          endif
        endif

        b_a = 1.0 / ((h_a + d_a(i,K-1)*mix(i,K-1)) + mix(i,K))
        if (mix(i,K) <= 0.0) then
          c_a(i,K) = 0.0 ; d_a(i,K) = 1.0
        else
          d_a(i,K) = b_a * (h_a + d_a(i,K-1)*mix(i,K-1)) ! = 1.0-c_a(i,K)
          c_a(i,K) = 1.0 ; if (d_a(i,K) > epsilon(b_a)) c_a(i,K) = b_a * mix(i,K)
        endif

        if (K>2) then
          Te_a(i,k-1) = b_a * (h_a*tv%T(i,j,k-1) + mix(i,K-1)*Te_a(i,k-2))
          Se_a(i,k-1) = b_a * (h_a*tv%S(i,j,k-1) + mix(i,K-1)*Se_a(i,k-2))
        else
          Te_a(i,k-1) = b_a * (h_a*tv%T(i,j,k-1))
          Se_a(i,k-1) = b_a * (h_a*tv%S(i,j,k-1))
        endif
      endif ; enddo ; enddo

      ! Determine which columns might have further instabilities.
      changed_col = 0
      do i=is,ie ; if (do_i(i)) then
        if (change_ct(i) == 0) then
          last_down(i) = .true. ; do_i(i) = .false.
        else
          changed_col = changed_col + 1 ; change_ct(i) = 0
        endif
      endif ; enddo
      if (changed_col == 0) exit ! No more columns are unstable.

      ! This is the same as above, but with the direction reversed (bottom to top)
      do K=nz,2,-1 ; do i=is,ie ; if (do_i(i)) then

        h_a = h(i,j,k-1) + h_neglect ; h_b = h(i,j,k) + h_neglect
        if (mix(i,K) <= 0.0) then
          if (is_unstable(dRho_dT(i,K), dRho_dS(i,K), h_a, h_b, mix(i,K-1), mix(i,K+1), &
                          tv%T(i,j,k-1), tv%T(i,j,k), tv%S(i,j,k-1), tv%S(i,j,k), &
                          Te_a(i,k-2), Te_b(i,k+1), Se_a(i,k-2), Se_b(i,k+1), &
                          d_a(i,K-1), d_b(i,K+1))) then
            mix(i,K) = mix_len
            if (kap_dt_x2 > 0.0) mix(i,K) = kap_dt_x2 / ((h(i,j,k-1)+h(i,j,k)) + h0)
            change_ct(i) = change_ct(i) + 1
          endif
        endif

        b_b = 1.0 / ((h_b + d_b(i,K+1)*mix(i,K+1)) + mix(i,K))
        if (mix(i,K) <= 0.0) then
          c_b(i,K) = 0.0 ; d_b(i,K) = 1.0
        else
          d_b(i,K) = b_b * (h_b + d_b(i,K+1)*mix(i,K+1)) ! = 1.0-c_b(i,K)
          c_b(i,K) = 1.0 ; if (d_b(i,K) > epsilon(b_b)) c_b(i,K) = b_b * mix(i,K)
        endif

        if (k<nz) then
          Te_b(i,k) = b_b * (h_b*tv%T(i,j,k) + mix(i,K+1)*Te_b(i,k+1))
          Se_b(i,k) = b_b * (h_b*tv%S(i,j,k) + mix(i,K+1)*Se_b(i,k+1))
        else
          Te_b(i,k) = b_b * (h_b*tv%T(i,j,k))
          Se_b(i,k) = b_b * (h_b*tv%S(i,j,k))
        endif
      endif ; enddo ; enddo

      ! Determine which columns might have further instabilities.
      changed_col = 0
      do i=is,ie ; if (do_i(i)) then
        if (change_ct(i) == 0) then
          last_down(i) = .false. ; do_i(i) = .false.
        else
          changed_col = changed_col + 1 ; change_ct(i) = 0
        endif
      endif ; enddo
      if (changed_col == 0) exit ! No more columns are unstable.

    enddo  ! End of iterations, all columns are now stable.

    ! Do the final return pass on the columns where the penultimate pass was downward.
    do i=is,ie ; do_i(i) = ((G%mask2dT(i,j) > 0.0) .and. last_down(i)) ; enddo
    do i=is,ie ; if (do_i(i)) then
      h_a = h(i,j,nz) + h_neglect
      b_a = 1.0 / (h_a + d_a(i,nz)*mix(i,nz))
      T_adj(i,j,nz) = b_a * (h_a*tv%T(i,j,nz) + mix(i,nz)*Te_a(i,nz-1))
      S_adj(i,j,nz) = b_a * (h_a*tv%S(i,j,nz) + mix(i,nz)*Se_a(i,nz-1))
    endif ; enddo
    do k=nz-1,1,-1 ; do i=is,ie ; if (do_i(i)) then
      T_adj(i,j,k) = Te_a(i,k) + c_a(i,K+1)*T_adj(i,j,k+1)
      S_adj(i,j,k) = Se_a(i,k) + c_a(i,K+1)*S_adj(i,j,k+1)
    endif ; enddo ; enddo

    do i=is,ie ; if (do_i(i)) then
      k = 1 ! A hook for debugging.
    endif ; enddo

    ! Do the final return pass on the columns where the penultimate pass was upward.
    ! Also do a simple copy of T & S values on land points.
    do i=is,ie
      do_i(i) = ((G%mask2dT(i,j) > 0.0) .and. .not.last_down(i))
      if (do_i(i)) then
        h_b = h(i,j,1) + h_neglect
        b_b = 1.0 / (h_b + d_b(i,2)*mix(i,2))
        T_adj(i,j,1) = b_b * (h_b*tv%T(i,j,1) + mix(i,2)*Te_b(i,2))
        S_adj(i,j,1) = b_b * (h_b*tv%S(i,j,1) + mix(i,2)*Se_b(i,2))
      elseif (G%mask2dT(i,j) <= 0.0) then
        T_adj(i,j,1) = tv%T(i,j,1) ; S_adj(i,j,1) = tv%S(i,j,1)
      endif
    enddo
    do k=2,nz ; do i=is,ie
      if (do_i(i)) then
        T_adj(i,j,k) = Te_b(i,k) + c_b(i,K)*T_adj(i,j,k-1)
        S_adj(i,j,k) = Se_b(i,k) + c_b(i,K)*S_adj(i,j,k-1)
      elseif (G%mask2dT(i,j) <= 0.0) then
        T_adj(i,j,k) = tv%T(i,j,k) ; S_adj(i,j,k) = tv%S(i,j,k)
      endif
    enddo ; enddo

    do i=is,ie ; if (do_i(i)) then
      k = 1 ! A hook for debugging.
    endif ; enddo

  enddo ! j-loop

  k = 1 ! A hook for debugging.

  ! The following set of expressions for the final values are derived from the the partial
  ! updates for the estimated temperatures and salinities around an interface, then directly
  ! solving for the final temperatures and salinities.  They are here for later reference
  ! and to document an intermediate step in the stability calculation.
    ! hp_a = (h_a + d_a(i,K-1)*mix(i,K-1))
    ! hp_b = (h_b + d_b(i,K+1)*mix(i,K+1))
    ! b2_c = 1.0 / (hp_a*hp_b + (hp_a + hp_b) * mix(i,K))
    ! Th_a = h_a*tv%T(i,j,k-1) + mix(i,K-1)*Te_a(i,k-2)
    ! Th_b = h_b*tv%T(i,j,k)   + mix(i,K+1)*Te_b(i,k+1)
    ! T_fin(i,k)   = ( (hp_a + mix(i,K)) * Th_b  + Th_a * mix(i,K) ) * b2_c
    ! T_fin(i,k-1) = ( (hp_b + mix(i,K)) * Th_a  + Th_b * mix(i,K) ) * b2_c
    ! Sh_a = h_a*tv%S(i,j,k-1) + mix(i,K-1)*Se_a(i,k-2)
    ! Sh_b = h_b*tv%S(i,j,k)   + mix(i,K+1)*Se_b(i,k+1)
    ! S_fin(i,k)   = ( (hp_a + mix(i,K)) * Sh_b  + Sh_a * mix(i,K) ) * b2_c
    ! S_fin(i,k-1) = ( (hp_b + mix(i,K)) * Sh_a  + Sh_b * mix(i,K) ) * b2_c

end subroutine full_convection

!> This function returns True if the profiles around the given interface will be
!! statically unstable after mixing above below.  The arguments are the ocean state
!! above and below, including partial calculations from a tridiagonal solver.
function is_unstable(dRho_dT, dRho_dS, h_a, h_b, mix_A, mix_B, T_a, T_b, S_a, S_b, &
                     Te_aa, Te_bb, Se_aa, Se_bb, d_A, d_B)
  real, intent(in) :: dRho_dT !< The derivative of in situ density with temperature [R degC-1 ~> kg m-3 degC-1]
  real, intent(in) :: dRho_dS !< The derivative of in situ density with salinity [R ppt-1 ~> kg m-3 ppt-1]
  real, intent(in) :: h_a     !< The thickness of the layer above [H ~> m or kg m-2]
  real, intent(in) :: h_b     !< The thickness of the layer below [H ~> m or kg m-2]
  real, intent(in) :: mix_A   !< The time integrated mixing rate of the interface above [H ~> m or kg m-2]
  real, intent(in) :: mix_B   !< The time integrated mixing rate of the interface below [H ~> m or kg m-2]
  real, intent(in) :: T_a     !< The initial temperature of the layer above [degC]
  real, intent(in) :: T_b     !< The initial temperature of the layer below [degC]
  real, intent(in) :: S_a     !< The initial salinity of the layer below [ppt]
  real, intent(in) :: S_b     !< The initial salinity of the layer below [ppt]
  real, intent(in) :: Te_aa   !< The estimated temperature two layers above rescaled by d_A [degC]
  real, intent(in) :: Te_bb   !< The estimated temperature two layers below rescaled by d_B [degC]
  real, intent(in) :: Se_aa   !< The estimated salinity two layers above rescaled by d_A [ppt]
  real, intent(in) :: Se_bb   !< The estimated salinity two layers below rescaled by d_B [ppt]
  real, intent(in) :: d_A     !< The rescaling dependency across the interface above, nondim.
  real, intent(in) :: d_B     !< The rescaling dependency across the interface below, nondim.
  logical :: is_unstable !< The return value, true if the profile is statically unstable
                         !! around the interface in question.

  ! These expressions for the local stability are long, but they have been carefully
  ! grouped for accuracy even when the mixing rates are huge or tiny, and common
  ! positive definite factors that would appear in the final expression for the
  ! locally referenced potential density difference across an interface have been omitted.
  is_unstable = (dRho_dT * ((h_a * h_b * (T_b - T_a) + &
                             mix_A*mix_B * (d_A*Te_bb - d_B*Te_aa)) + &
                            (h_a*mix_B * (Te_bb - d_B*T_a) + &
                             h_b*mix_A * (d_A*T_b - Te_aa)) ) + &
                 dRho_dS * ((h_a * h_b * (S_b - S_a) + &
                             mix_A*mix_B * (d_A*Se_bb - d_B*Se_aa)) + &
                            (h_a*mix_B * (Se_bb - d_B*S_a) + &
                             h_b*mix_A * (d_A*S_b - Se_aa)) ) < 0.0)
end function is_unstable

!> Returns the partial derivatives of locally referenced potential density with
!! temperature and salinity after the properties have been smoothed with a small
!! constant diffusivity.
subroutine smoothed_dRdT_dRdS(h, tv, Kddt, dR_dT, dR_dS, G, GV, US, j, p_surf, halo)
  type(ocean_grid_type),   intent(in)  :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)  :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in)  :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)  :: tv   !< A structure pointing to various
                                               !! thermodynamic variables
  real,                    intent(in)  :: Kddt !< A diffusivity times a time increment [H2 ~> m2 or kg2 m-4].
  real, dimension(SZI_(G),SZK_(G)+1), &
                           intent(out) :: dR_dT !< Derivative of locally referenced
                                               !! potential density with temperature [R degC-1 ~> kg m-3 degC-1]
  real, dimension(SZI_(G),SZK_(G)+1), &
                           intent(out) :: dR_dS !< Derivative of locally referenced
                                               !! potential density with salinity [R degC-1 ~> kg m-3 ppt-1]
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  integer,                 intent(in)  :: j    !< The j-point to work on.
  real, dimension(:,:),    pointer     :: p_surf !< The pressure at the ocean surface [R L2 T-2 ~> Pa].
  integer,       optional, intent(in)  :: halo !< Halo width over which to compute

  ! Local variables
  real :: mix(SZI_(G),SZK_(G)+1)   ! The diffusive mixing length (kappa*dt)/dz
                                   ! between layers within in a timestep [H ~> m or kg m-2].
  real :: b1(SZI_(G)), d1(SZI_(G)) ! b1, c1, and d1 are variables used by the
  real :: c1(SZI_(G),SZK_(G))      ! tridiagonal solver.
  real :: T_f(SZI_(G),SZK_(G))     ! Filtered temperatures [degC]
  real :: S_f(SZI_(G),SZK_(G))     ! Filtered salinities [ppt]
  real :: pres(SZI_(G))            ! Interface pressures [R L2 T-2 ~> Pa].
  real :: T_EOS(SZI_(G))           ! Filtered and vertically averaged temperatures [degC]
  real :: S_EOS(SZI_(G))           ! Filtered and vertically averaged salinities [ppt]
  real :: kap_dt_x2                ! The product of 2*kappa*dt [H2 ~> m2 or kg2 m-4].
  real :: h_neglect, h0            ! Negligible thicknesses to allow for zero thicknesses,
                                   ! [H ~> m or kg m-2].
  real :: h_tr                     ! The thickness at tracer points, plus h_neglect [H ~> m or kg m-2].
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, k, is, ie, nz

  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo
  else
    is = G%isc ; ie = G%iec
  endif
  nz = G%ke

  h_neglect = GV%H_subroundoff
  kap_dt_x2 = 2.0*Kddt

  if (Kddt <= 0.0) then
    do k=1,nz ; do i=is,ie
      T_f(i,k) = tv%T(i,j,k) ; S_f(i,k) = tv%S(i,j,k)
    enddo ; enddo
  else
    h0 = 1.0e-16*sqrt(Kddt) + h_neglect
    do i=is,ie
      mix(i,2) = kap_dt_x2 / ((h(i,j,1)+h(i,j,2)) + h0)

      h_tr = h(i,j,1) + h_neglect
      b1(i) = 1.0 / (h_tr + mix(i,2))
      d1(i) = b1(i) * h(i,j,1)
      T_f(i,1) = (b1(i)*h_tr)*tv%T(i,j,1)
      S_f(i,1) = (b1(i)*h_tr)*tv%S(i,j,1)
    enddo
    do k=2,nz-1 ; do i=is,ie
      mix(i,K+1) = kap_dt_x2 / ((h(i,j,k)+h(i,j,k+1)) + h0)

      c1(i,k) = mix(i,K) * b1(i)
      h_tr = h(i,j,k) + h_neglect
      b1(i) = 1.0 / ((h_tr + d1(i)*mix(i,K)) + mix(i,K+1))
      d1(i) = b1(i) * (h_tr + d1(i)*mix(i,K))
      T_f(i,k) = b1(i) * (h_tr*tv%T(i,j,k) + mix(i,K)*T_f(i,k-1))
      S_f(i,k) = b1(i) * (h_tr*tv%S(i,j,k) + mix(i,K)*S_f(i,k-1))
    enddo ; enddo
    do i=is,ie
      c1(i,nz) = mix(i,nz) * b1(i)
      h_tr = h(i,j,nz) + h_neglect
      b1(i) = 1.0 / (h_tr + d1(i)*mix(i,nz))
      T_f(i,nz) = b1(i) * (h_tr*tv%T(i,j,nz) + mix(i,nz)*T_f(i,nz-1))
      S_f(i,nz) = b1(i) * (h_tr*tv%S(i,j,nz) + mix(i,nz)*S_f(i,nz-1))
    enddo
    do k=nz-1,1,-1 ; do i=is,ie
      T_f(i,k) = T_f(i,k) + c1(i,k+1)*T_f(i,k+1)
      S_f(i,k) = S_f(i,k) + c1(i,k+1)*S_f(i,k+1)
    enddo ; enddo
  endif

  if (associated(p_surf)) then
    do i=is,ie ; pres(i) = p_surf(i,j) ; enddo
  else
    do i=is,ie ; pres(i) = 0.0 ; enddo
  endif
  EOSdom(:) = EOS_domain(G%HI)
  call calculate_density_derivs(T_f(:,1), S_f(:,1), pres, dR_dT(:,1), dR_dS(:,1), tv%eqn_of_state, EOSdom)
  do i=is,ie ; pres(i) = pres(i) + h(i,j,1)*(GV%H_to_RZ*GV%g_Earth) ; enddo
  do K=2,nz
    do i=is,ie
      T_EOS(i) = 0.5*(T_f(i,k-1) + T_f(i,k))
      S_EOS(i) = 0.5*(S_f(i,k-1) + S_f(i,k))
    enddo
    call calculate_density_derivs(T_EOS, S_EOS, pres, dR_dT(:,K), dR_dS(:,K), tv%eqn_of_state, EOSdom)
    do i=is,ie ; pres(i) = pres(i) + h(i,j,k)*(GV%H_to_RZ*GV%g_Earth) ; enddo
  enddo
  call calculate_density_derivs(T_f(:,nz), S_f(:,nz), pres, dR_dT(:,nz+1), dR_dS(:,nz+1), &
                                tv%eqn_of_state, EOSdom)

end subroutine smoothed_dRdT_dRdS

end module MOM_full_convection
