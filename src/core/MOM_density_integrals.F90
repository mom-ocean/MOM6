!> Provides integrals of density
module MOM_density_integrals

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS,              only : EOS_type
use MOM_EOS,              only : EOS_quadrature
use MOM_EOS,              only : analytic_int_density_dz
use MOM_EOS,              only : analytic_int_specific_vol_dp
use MOM_EOS,              only : calculate_density
use MOM_EOS,              only : calculate_spec_vol
use MOM_EOS,              only : calculate_specific_vol_derivs
use MOM_error_handler,    only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_file_parser,      only : get_param, log_version, param_file_type
use MOM_hor_index,        only : hor_index_type
use MOM_string_functions, only : uppercase
use MOM_variables,        only : thermo_var_ptrs
use MOM_unit_scaling,     only : unit_scale_type
use MOM_verticalGrid,     only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public int_density_dz
public int_density_dz_generic_pcm
public int_density_dz_generic_plm
public int_density_dz_generic_ppm
public int_specific_vol_dp
public int_spec_vol_dp_generic_pcm
public int_spec_vol_dp_generic_plm
public find_depth_of_pressure_in_cell

contains

!> Calls the appropriate subroutine to calculate analytical and nearly-analytical
!! integrals in z across layers of pressure anomalies, which are
!! required for calculating the finite-volume form pressure accelerations in a
!! Boussinesq model.
subroutine int_density_dz(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, EOS, US, dpa, &
                          intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp, Z_0p)
  type(hor_index_type), intent(in)  :: HI  !< Ocean horizontal index structures for the arrays
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S   !< Salinity [ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: z_t !< Height at the top of the layer in depth units [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: z_b !< Height at the bottom of the layer [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3] or [kg m-3], that is
                                           !! subtracted out to reduce the magnitude of each of the
                                           !! integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3] or [kg m-3], that is used
                                           !! to calculate the pressure (as p~=-z*rho_0*G_e)
                                           !! used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration
                                           !! [L2 Z-1 T-2 ~> m s-2] or [m2 Z-1 s-2 ~> m s-2]
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                      intent(inout) :: dpa !< The change in the pressure anomaly
                                           !! across the layer [R L2 T-2 ~> Pa] or [Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
            optional, intent(inout) :: intz_dpa !< The integral through the thickness of the
                                           !! layer of the pressure anomaly relative to the
                                           !! anomaly at the top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
            optional, intent(inout) :: intx_dpa !< The integral in x of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJB_(HI)), &
            optional, intent(inout) :: inty_dpa !< The integral in y of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the y grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  real,       optional, intent(in)  :: dz_neglect !< A minuscule thickness change [Z ~> m]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting to
                                           !! interpolate T/S for top and bottom integrals.
  real,       optional, intent(in)  :: Z_0p !< The height at which the pressure is 0 [Z ~> m]

  if (EOS_quadrature(EOS)) then
    call int_density_dz_generic_pcm(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, EOS, US, dpa, &
                                    intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp, Z_0p=Z_0p)
  else
    call analytic_int_density_dz(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, EOS, dpa, &
                                 intz_dpa, intx_dpa, inty_dpa, bathyT, dz_neglect, useMassWghtInterp, Z_0p=Z_0p)
  endif

end subroutine int_density_dz


!> Calculates (by numerical quadrature) integrals of pressure anomalies across layers, which
!! are required for calculating the finite-volume form pressure accelerations in a Boussinesq model.
subroutine int_density_dz_generic_pcm(T, S, z_t, z_b, rho_ref, rho_0, G_e, HI, &
                                      EOS, US, dpa, intz_dpa, intx_dpa, inty_dpa, bathyT, &
                                      dz_neglect, useMassWghtInterp, use_inaccurate_form, Z_0p)
  type(hor_index_type), intent(in)  :: HI  !< Horizontal index type for input variables.
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T  !< Potential temperature of the layer [degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S  !< Salinity of the layer [ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: z_t !< Height at the top of the layer in depth units [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: z_b !< Height at the bottom of the layer [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3] or [kg m-3], that is
                                          !! subtracted out to reduce the magnitude
                                          !! of each of the integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3] or [kg m-3], that is used
                                          !! to calculate the pressure (as p~=-z*rho_0*G_e)
                                          !! used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration
                                          !! [L2 Z-1 T-2 ~> m s-2] or [m2 Z-1 s-2 ~> m s-2]
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                      intent(inout) :: dpa !< The change in the pressure anomaly
                                          !! across the layer [R L2 T-2 ~> Pa] or [Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
            optional, intent(inout) :: intz_dpa !< The integral through the thickness of the
                                          !! layer of the pressure anomaly relative to the
                                          !! anomaly at the top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
            optional, intent(inout) :: intx_dpa !< The integral in x of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJB_(HI)), &
            optional, intent(inout) :: inty_dpa !< The integral in y of the difference between
                                          !! the pressure anomaly at the top and bottom of the
                                          !! layer divided by the y grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  real,       optional, intent(in)  :: dz_neglect !< A minuscule thickness change [Z ~> m]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting to
                                          !! interpolate T/S for top and bottom integrals.
  logical,    optional, intent(in)  :: use_inaccurate_form !< If true, uses an inaccurate form of
                                          !! density anomalies, as was used prior to March 2018.
  real,       optional, intent(in)  :: Z_0p !< The height at which the pressure is 0 [Z ~> m]

  ! Local variables
  real :: T5(5), S5(5) ! Temperatures and salinities at five quadrature points [degC] and [ppt]
  real :: p5(5)      ! Pressures at five quadrature points, never rescaled from Pa [Pa]
  real :: r5(5)      ! Densities at five quadrature points [R ~> kg m-3] or [kg m-3]
  real :: rho_anom   ! The depth averaged density anomaly [R ~> kg m-3] or [kg m-3]
  real :: w_left, w_right ! Left and right weights [nondim]
  real, parameter :: C1_90 = 1.0/90.0  ! Rational constants.
  real :: GxRho      ! The gravitational acceleration times density and unit conversion factors [Pa Z-1 ~> kg m-2 s-2]
  real :: I_Rho      ! The inverse of the Boussinesq density [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: rho_scale  ! A scaling factor for densities from kg m-3 to R [R m3 kg-1 ~> 1]
  real :: rho_ref_mks ! The reference density in MKS units, never rescaled from kg m-3 [kg m-3]
  real :: dz         ! The layer thickness [Z ~> m]
  real :: z0pres     ! The height at which the pressure is zero [Z ~> m]
  real :: hWght      ! A pressure-thickness below topography [Z ~> m]
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [Z ~> m]
  real :: iDenom     ! The inverse of the denominator in the weights [Z-2 ~> m-2]
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim]
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim]
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim]
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim]
  real :: intz(5)    ! The gravitational acceleration times the integrals of density
                     ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa] or [Pa]
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  logical :: use_rho_ref ! Pass rho_ref to the equation of state for more accurate calculation
                         ! of density anomalies.
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j, m, n

  ! These array bounds work for the indexing convention of the input arrays, but
  ! on the computational domain defined for the output arrays.
  Isq = HI%IscB ; Ieq = HI%IecB
  Jsq = HI%JscB ; Jeq = HI%JecB
  is = HI%isc ; ie = HI%iec
  js = HI%jsc ; je = HI%jec

  rho_scale = US%kg_m3_to_R
  GxRho = US%RL2_T2_to_Pa * G_e * rho_0
  rho_ref_mks = rho_ref * US%R_to_kg_m3
  I_Rho = 1.0 / rho_0
  z0pres = 0.0 ; if (present(Z_0p)) z0pres = Z_0p
  use_rho_ref = .true.
  if (present(use_inaccurate_form)) then
    if (use_inaccurate_form) use_rho_ref = .not. use_inaccurate_form
  endif

  do_massWeight = .false.
  if (present(useMassWghtInterp)) then ; if (useMassWghtInterp) then
    do_massWeight = .true.
    if (.not.present(bathyT)) call MOM_error(FATAL, "int_density_dz_generic: "//&
        "bathyT must be present if useMassWghtInterp is present and true.")
    if (.not.present(dz_neglect)) call MOM_error(FATAL, "int_density_dz_generic: "//&
        "dz_neglect must be present if useMassWghtInterp is present and true.")
  endif ; endif

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    dz = z_t(i,j) - z_b(i,j)
    do n=1,5
      T5(n) = T(i,j) ; S5(n) = S(i,j)
      p5(n) = -GxRho*((z_t(i,j) - z0pres) - 0.25*real(n-1)*dz)
    enddo
    if (use_rho_ref) then
      if (rho_scale /= 1.0) then
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
      else
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks)
      endif
      ! Use Boole's rule to estimate the pressure anomaly change.
      rho_anom = C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3))
    else
      if (rho_scale /= 1.0) then
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, scale=rho_scale)
      else
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS)
      endif
      ! Use Boole's rule to estimate the pressure anomaly change.
      rho_anom = C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)) - rho_ref
    endif

    dpa(i,j) = G_e*dz*rho_anom
    ! Use a Boole's-rule-like fifth-order accurate estimate of the double integral of
    ! the pressure anomaly.
    if (present(intz_dpa)) intz_dpa(i,j) = 0.5*G_e*dz**2 * &
          (rho_anom - C1_90*(16.0*(r5(4)-r5(2)) + 7.0*(r5(5)-r5(1))) )
  enddo ; enddo

  if (present(intx_dpa)) then ; do j=js,je ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., -bathyT(i,j)-z_t(i+1,j), -bathyT(i+1,j)-z_t(i,j))
    if (hWght > 0.) then
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i+1,j) - z_b(i+1,j)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)
    do m=2,4
      ! T, S, and z are interpolated in the horizontal.  The z interpolation
      ! is linear, but for T and S it may be thickness weighted.
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR
      dz = wt_L*(z_t(i,j) - z_b(i,j)) + wt_R*(z_t(i+1,j) - z_b(i+1,j))
      T5(1) = wtT_L*T(i,j) + wtT_R*T(i+1,j)
      S5(1) = wtT_L*S(i,j) + wtT_R*S(i+1,j)
      p5(1) = -GxRho*((wt_L*z_t(i,j) + wt_R*z_t(i+1,j)) - z0pres)
      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1) ; p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo
      if (use_rho_ref) then
        if (rho_scale /= 1.0) then
          call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
        else
          call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks)
        endif
        ! Use Boole's rule to estimate the pressure anomaly change.
        intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)))
      else
        if (rho_scale /= 1.0) then
          call calculate_density(T5, S5, p5, r5, 1, 5, EOS, scale=rho_scale)
        else
          call calculate_density(T5, S5, p5, r5, 1, 5, EOS)
        endif
        ! Use Boole's rule to estimate the pressure anomaly change.
        intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)) - rho_ref )
      endif

    enddo
    ! Use Boole's rule to integrate the bottom pressure anomaly values in x.
    intx_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                           12.0*intz(3))
  enddo ; enddo ; endif

  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=is,ie
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., -bathyT(i,j)-z_t(i,j+1), -bathyT(i,j+1)-z_t(i,j))
    if (hWght > 0.) then
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i,j+1) - z_b(i,j+1)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)
    do m=2,4
      ! T, S, and z are interpolated in the horizontal.  The z interpolation
      ! is linear, but for T and S it may be thickness weighted.
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR
      dz = wt_L*(z_t(i,j) - z_b(i,j)) + wt_R*(z_t(i,j+1) - z_b(i,j+1))
      T5(1) = wtT_L*T(i,j) + wtT_R*T(i,j+1)
      S5(1) = wtT_L*S(i,j) + wtT_R*S(i,j+1)
      p5(1) = -GxRho*((wt_L*z_t(i,j) + wt_R*z_t(i,j+1)) - z0pres)
      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1)
        p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo
      if (use_rho_ref) then
        if (rho_scale /= 1.0) then
          call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
        else
          call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks)
        endif
        ! Use Boole's rule to estimate the pressure anomaly change.
        intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)))
      else
        if (rho_scale /= 1.0) then
          call calculate_density(T5, S5, p5, r5, 1, 5, EOS, scale=rho_scale)
        else
          call calculate_density(T5, S5, p5, r5, 1, 5, EOS)
        endif
        ! Use Boole's rule to estimate the pressure anomaly change.
        intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)) - rho_ref )
      endif

    enddo
    ! Use Boole's rule to integrate the values.
    inty_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                                     12.0*intz(3))
  enddo ; enddo ; endif
end subroutine int_density_dz_generic_pcm


!> Compute pressure gradient force integrals by quadrature for the case where
!! T and S are linear profiles.
subroutine int_density_dz_generic_plm(k, tv, T_t, T_b, S_t, S_b, e, rho_ref, &
                                      rho_0, G_e, dz_subroundoff, bathyT, HI, GV, EOS, US, dpa, &
                                      intz_dpa, intx_dpa, inty_dpa, useMassWghtInterp, &
                                      use_inaccurate_form, Z_0p)
  integer,              intent(in)  :: k   !< Layer index to calculate integrals for
  type(hor_index_type), intent(in)  :: HI  !< Ocean horizontal index structures for the input arrays
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(thermo_var_ptrs), intent(in) :: tv  !< Thermodynamic variables
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: T_t !< Potential temperature at the cell top [degC]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: T_b !< Potential temperature at the cell bottom [degC]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: S_t !< Salinity at the cell top [ppt]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: S_b !< Salinity at the cell bottom [ppt]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)+1), &
                        intent(in)  :: e   !< Height of interfaces [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3] or [kg m-3], that is subtracted
                                           !! out to reduce the magnitude of each of the integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3] or [kg m-3], that is used to calculate
                                           !! the pressure (as p~=-z*rho_0*G_e) used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,                 intent(in)  :: dz_subroundoff !< A minuscule thickness change [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dpa !< The change in the pressure anomaly across the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intz_dpa !< The integral through the thickness of the layer of
                                           !! the pressure anomaly relative to the anomaly at the
                                           !! top of the layer [R L2 Z T-2 ~> Pa Z]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dpa !< The integral in x of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dpa !< The integral in y of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the y grid spacing [R L2 T-2 ~> Pa]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting to
                                           !! interpolate T/S for top and bottom integrals.
  logical,    optional, intent(in)  :: use_inaccurate_form !< If true, uses an inaccurate form of
                                           !! density anomalies, as was used prior to March 2018.
  real,       optional, intent(in)  :: Z_0p !< The height at which the pressure is 0 [Z ~> m]

! This subroutine calculates (by numerical quadrature) integrals of
! pressure anomalies across layers, which are required for calculating the
! finite-volume form pressure accelerations in a Boussinesq model.  The one
! potentially dodgy assumption here is that rho_0 is used both in the denominator
! of the accelerations, and in the pressure used to calculated density (the
! latter being -z*rho_0*G_e).  These two uses could be separated if need be.
!
! It is assumed that the salinity and temperature profiles are linear in the
! vertical. The top and bottom values within each layer are provided and
! a linear interpolation is used to compute intermediate values.

  ! Local variables
  real :: T5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Temperatures along a line of subgrid locations [degC]
  real :: S5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Salinities along a line of subgrid locations [ppt]
  real :: T25((5*HI%iscB+1):(5*(HI%iecB+2))) ! SGS temperature variance along a line of subgrid locations [degC2]
  real :: TS5((5*HI%iscB+1):(5*(HI%iecB+2))) ! SGS temp-salt covariance along a line of subgrid locations [degC ppt]
  real :: S25((5*HI%iscB+1):(5*(HI%iecB+2))) ! SGS salinity variance along a line of subgrid locations [ppt2]
  real :: p5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Pressures along a line of subgrid locations, never
                                             ! rescaled from Pa [Pa]
  real :: r5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Densities anomalies along a line of subgrid
                                             ! locations [R ~> kg m-3] or [kg m-3]
  real :: u5((5*HI%iscB+1):(5*(HI%iecB+2)))  ! Densities anomalies along a line of subgrid locations
                                             ! (used for inaccurate form) [R ~> kg m-3] or [kg m-3]
  real :: T15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Temperatures at an array of subgrid locations [degC]
  real :: S15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Salinities at an array of subgrid locations [ppt]
  real :: T215((15*HI%iscB+1):(15*(HI%iecB+1))) ! SGS temperature variance along a line of subgrid locations [degC2]
  real :: TS15((15*HI%iscB+1):(15*(HI%iecB+1))) ! SGS temp-salt covariance along a line of subgrid locations [degC ppt]
  real :: S215((15*HI%iscB+1):(15*(HI%iecB+1))) ! SGS salinity variance along a line of subgrid locations [ppt2]
  real :: p15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Pressures at an array of subgrid locations [Pa]
  real :: r15((15*HI%iscB+1):(15*(HI%iecB+1))) ! Densities at an array of subgrid locations
                                                 ! [R ~> kg m-3] or [kg m-3]
  real :: wt_t(5), wt_b(5)          ! Top and bottom weights [nondim]
  real :: rho_anom                  ! A density anomaly [R ~> kg m-3] or [kg m-3]
  real :: w_left, w_right           ! Left and right weights [nondim]
  real :: intz(5)    ! The gravitational acceleration times the integrals of density
                     ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa] or [Pa]
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant [nondim]
  real :: GxRho      ! The gravitational acceleration times density and unit conversion factors [Pa Z-1 ~> kg m-2 s-2]
  real :: I_Rho      ! The inverse of the Boussinesq density [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: rho_scale  ! A scaling factor for densities from kg m-3 to R [R m3 kg-1 ~> 1]
  real :: rho_ref_mks ! The reference density in MKS units, never rescaled from kg m-3 [kg m-3]
  real :: dz(HI%iscB:HI%iecB+1)   ! Layer thicknesses at tracer points [Z ~> m]
  real :: dz_x(5,HI%iscB:HI%iecB) ! Layer thicknesses along an x-line of subrid locations [Z ~> m]
  real :: dz_y(5,HI%isc:HI%iec)   ! Layer thicknesses along a y-line of subrid locations [Z ~> m]
  real :: massWeightToggle          ! A non-dimensional toggle factor (0 or 1) [nondim]
  real :: Ttl, Tbl, Ttr, Tbr        ! Temperatures at the velocity cell corners [degC]
  real :: Stl, Sbl, Str, Sbr        ! Salinities at the velocity cell corners [ppt]
  real :: z0pres                    ! The height at which the pressure is zero [Z ~> m]
  real :: hWght                     ! A topographically limited thicknes weight [Z ~> m]
  real :: hL, hR                    ! Thicknesses to the left and right [Z ~> m]
  real :: iDenom                    ! The denominator of the thickness weight expressions [Z-2 ~> m-2]
  logical :: use_stanley_eos ! True is SGS variance fields exist in tv.
  logical :: use_rho_ref ! Pass rho_ref to the equation of state for more accurate calculation
                         ! of density anomalies.
  logical :: use_varT, use_varS, use_covarTS ! Logicals for SGS variances fields
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n
  integer :: pos

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB

  rho_scale = US%kg_m3_to_R
  GxRho = US%RL2_T2_to_Pa * G_e * rho_0
  rho_ref_mks = rho_ref * US%R_to_kg_m3
  I_Rho = 1.0 / rho_0
  z0pres = 0.0 ; if (present(Z_0p)) z0pres = Z_0p
  massWeightToggle = 0.
  if (present(useMassWghtInterp)) then
    if (useMassWghtInterp) massWeightToggle = 1.
  endif
  use_rho_ref = .true.
  if (present(use_inaccurate_form)) then
    if (use_inaccurate_form) use_rho_ref = .not. use_inaccurate_form
  endif

  use_varT = associated(tv%varT)
  use_covarTS = associated(tv%covarTS)
  use_varS = associated(tv%varS)
  use_stanley_eos = use_varT .or. use_covarTS .or. use_varS
  T25(:) = 0.
  TS5(:) = 0.
  S25(:) = 0.
  T215(:) = 0.
  TS15(:) = 0.
  S215(:) = 0.

  do n = 1, 5
    wt_t(n) = 0.25 * real(5-n)
    wt_b(n) = 1.0 - wt_t(n)
  enddo

  ! 1. Compute vertical integrals
  do j=Jsq,Jeq+1
    do i = Isq,Ieq+1
      dz(i) = e(i,j,K) - e(i,j,K+1)
      do n=1,5
        p5(i*5+n) = -GxRho*((e(i,j,K) - z0pres) - 0.25*real(n-1)*dz(i))
        ! Salinity and temperature points are linearly interpolated
        S5(i*5+n) = wt_t(n) * S_t(i,j,k) + wt_b(n) * S_b(i,j,k)
        T5(i*5+n) = wt_t(n) * T_t(i,j,k) + wt_b(n) * T_b(i,j,k)
      enddo
      if (use_varT) T25(i*5+1:i*5+5) = tv%varT(i,j,k)
      if (use_covarTS) TS5(i*5+1:i*5+5) = tv%covarTS(i,j,k)
      if (use_varS) S25(i*5+1:i*5+5) = tv%varS(i,j,k)
    enddo
    if (use_Stanley_eos) then
      if (rho_scale /= 1.0) then
        call calculate_density(T5, S5, p5, T25, TS5, S25, r5, 1, (ieq-isq+2)*5, EOS, &
                               rho_ref=rho_ref_mks, scale=rho_scale)
      else
        call calculate_density(T5, S5, p5, T25, TS5, S25, r5, 1, (ieq-isq+2)*5, EOS, &
                               rho_ref=rho_ref_mks)
      endif
    else
      if (use_rho_ref) then
        if (rho_scale /= 1.0) then
          call calculate_density(T5, S5, p5, r5, 1, (ieq-isq+2)*5, EOS, rho_ref=rho_ref_mks, &
                                 scale=rho_scale)
        else
          call calculate_density(T5, S5, p5, r5, 1, (ieq-isq+2)*5, EOS, rho_ref=rho_ref_mks)
        endif
      else
        if (rho_scale /= 1.0) then
          call calculate_density(T5, S5, p5, r5, 1, (ieq-isq+2)*5, EOS, scale=rho_scale)
        else
          call calculate_density(T5, S5, p5, r5, 1, (ieq-isq+2)*5, EOS)
        endif
        u5(:) = r5(:) - rho_ref
      endif
    endif

    if (use_rho_ref) then
      do i=Isq,Ieq+1
        ! Use Boole's rule to estimate the pressure anomaly change.
        rho_anom = C1_90*(7.0*(r5(i*5+1)+r5(i*5+5)) + 32.0*(r5(i*5+2)+r5(i*5+4)) + 12.0*r5(i*5+3))
        dpa(i,j) = G_e*dz(i)*rho_anom
        if (present(intz_dpa)) then
        ! Use a Boole's-rule-like fifth-order accurate estimate of
        ! the double integral of the pressure anomaly.
          intz_dpa(i,j) = 0.5*G_e*dz(i)**2 * &
                  (rho_anom - C1_90*(16.0*(r5(i*5+4)-r5(i*5+2)) + 7.0*(r5(i*5+5)-r5(i*5+1))) )
        endif
      enddo
    else
      do i=Isq,Ieq+1
        ! Use Boole's rule to estimate the pressure anomaly change.
        rho_anom = C1_90*(7.0*(r5(i*5+1)+r5(i*5+5)) + 32.0*(r5(i*5+2)+r5(i*5+4)) + 12.0*r5(i*5+3)) &
                   - rho_ref
        dpa(i,j) = G_e*dz(i)*rho_anom
        if (present(intz_dpa)) then
        ! Use a Boole's-rule-like fifth-order accurate estimate of
        ! the double integral of the pressure anomaly.
          intz_dpa(i,j) = 0.5*G_e*dz(i)**2 * &
                  (rho_anom - C1_90*(16.0*(u5(i*5+4)-u5(i*5+2)) + 7.0*(u5(i*5+5)-u5(i*5+1))) )
        endif
      enddo
    endif
  enddo ! end loops on j

  ! 2. Compute horizontal integrals in the x direction
  if (present(intx_dpa)) then ; do j=HI%jsc,HI%jec
    do I=Isq,Ieq
      ! Corner values of T and S
      ! hWght is the distance measure by which the cell is violation of
      ! hydrostatic consistency. For large hWght we bias the interpolation
      ! of T,S along the top and bottom integrals, almost like thickness
      ! weighting.
      ! Note: To work in terrain following coordinates we could offset
      ! this distance by the layer thickness to replicate other models.
      hWght = massWeightToggle * &
              max(0., -bathyT(i,j)-e(i+1,j,K), -bathyT(i+1,j)-e(i,j,K))
      if (hWght > 0.) then
        hL = (e(i,j,K) - e(i,j,K+1)) + dz_subroundoff
        hR = (e(i+1,j,K) - e(i+1,j,K+1)) + dz_subroundoff
        hWght = hWght * ( (hL-hR)/(hL+hR) )**2
        iDenom = 1./( hWght*(hR + hL) + hL*hR )
        Ttl = ( (hWght*hR)*T_t(i+1,j,k) + (hWght*hL + hR*hL)*T_t(i,j,k) ) * iDenom
        Ttr = ( (hWght*hL)*T_t(i,j,k) + (hWght*hR + hR*hL)*T_t(i+1,j,k) ) * iDenom
        Tbl = ( (hWght*hR)*T_b(i+1,j,k) + (hWght*hL + hR*hL)*T_b(i,j,k) ) * iDenom
        Tbr = ( (hWght*hL)*T_b(i,j,k) + (hWght*hR + hR*hL)*T_b(i+1,j,k) ) * iDenom
        Stl = ( (hWght*hR)*S_t(i+1,j,k) + (hWght*hL + hR*hL)*S_t(i,j,k) ) * iDenom
        Str = ( (hWght*hL)*S_t(i,j,k) + (hWght*hR + hR*hL)*S_t(i+1,j,k) ) * iDenom
        Sbl = ( (hWght*hR)*S_b(i+1,j,k) + (hWght*hL + hR*hL)*S_b(i,j,k) ) * iDenom
        Sbr = ( (hWght*hL)*S_b(i,j,k) + (hWght*hR + hR*hL)*S_b(i+1,j,k) ) * iDenom
      else
        Ttl = T_t(i,j,k); Tbl = T_b(i,j,k); Ttr = T_t(i+1,j,k); Tbr = T_b(i+1,j,k)
        Stl = S_t(i,j,k); Sbl = S_b(i,j,k); Str = S_t(i+1,j,k); Sbr = S_b(i+1,j,k)
      endif

      do m=2,4
        w_left = wt_t(m) ; w_right = wt_b(m)
        dz_x(m,i) = w_left*(e(i,j,K) - e(i,j,K+1)) + w_right*(e(i+1,j,K) - e(i+1,j,K+1))

        ! Salinity and temperature points are linearly interpolated in
        ! the horizontal. The subscript (1) refers to the top value in
        ! the vertical profile while subscript (5) refers to the bottom
        ! value in the vertical profile.
        pos = i*15+(m-2)*5
        T15(pos+1) = w_left*Ttl + w_right*Ttr
        T15(pos+5) = w_left*Tbl + w_right*Tbr

        S15(pos+1) = w_left*Stl + w_right*Str
        S15(pos+5) = w_left*Sbl + w_right*Sbr

        p15(pos+1) = -GxRho*((w_left*e(i,j,K) + w_right*e(i+1,j,K)) - z0pres)

        ! Pressure
        do n=2,5
          p15(pos+n) = p15(pos+n-1) + GxRho*0.25*dz_x(m,i)
        enddo

        ! Salinity and temperature (linear interpolation in the vertical)
        do n=2,4
          S15(pos+n) = wt_t(n) * S15(pos+1) + wt_b(n) * S15(pos+5)
          T15(pos+n) = wt_t(n) * T15(pos+1) + wt_b(n) * T15(pos+5)
        enddo
        if (use_varT) T215(pos+1:pos+5) = w_left*tv%varT(i,j,k) + w_right*tv%varT(i+1,j,k)
        if (use_covarTS) TS15(pos+1:pos+5) = w_left*tv%covarTS(i,j,k) + w_right*tv%covarTS(i+1,j,k)
        if (use_varS) S215(pos+1:pos+5) = w_left*tv%varS(i,j,k) + w_right*tv%varS(i+1,j,k)
      enddo
    enddo

    if (use_stanley_eos) then
      if (rho_scale /= 1.0) then
        call calculate_density(T15, S15, p15, T215, TS15, S215, r15, 1, 15*(ieq-isq+1), EOS, &
                               rho_ref=rho_ref_mks, scale=rho_scale)
      else
        call calculate_density(T15, S15, p15, T215, TS15, S215, r15, 1, 15*(ieq-isq+1), EOS, &
                               rho_ref=rho_ref_mks)
      endif
    else
      if (use_rho_ref) then
        if (rho_scale /= 1.0) then
          call calculate_density(T15, S15, p15, r15, 1, 15*(ieq-isq+1), EOS, rho_ref=rho_ref_mks, &
                                 scale=rho_scale)
        else
          call calculate_density(T15, S15, p15, r15, 1, 15*(ieq-isq+1), EOS, rho_ref=rho_ref_mks)
        endif
      else
        if (rho_scale /= 1.0) then
          call calculate_density(T15, S15, p15, r15, 1, 15*(ieq-isq+1), EOS, scale=rho_scale)
        else
          call calculate_density(T15, S15, p15, r15, 1, 15*(ieq-isq+1), EOS)
        endif
      endif
    endif

    do I=Isq,Ieq
      intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)

      ! Use Boole's rule to estimate the pressure anomaly change.
      if (use_rho_ref) then
        do m = 2,4
          pos = i*15+(m-2)*5
          intz(m) = G_e*dz_x(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + 32.0*(r15(pos+2)+r15(pos+4)) + &
                            12.0*r15(pos+3)) )
        enddo
      else
        do m = 2,4
          pos = i*15+(m-2)*5
          intz(m) = G_e*dz_x(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + 32.0*(r15(pos+2)+r15(pos+4)) + &
                            12.0*r15(pos+3)) - rho_ref )
        enddo
      endif
      ! Use Boole's rule to integrate the bottom pressure anomaly values in x.
      intx_dpa(I,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    enddo
  enddo ; endif

  ! 3. Compute horizontal integrals in the y direction
  if (present(inty_dpa)) then ; do J=Jsq,Jeq
    do i=HI%isc,HI%iec
    ! Corner values of T and S
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, almost like thickness
    ! weighting.
    ! Note: To work in terrain following coordinates we could offset
    ! this distance by the layer thickness to replicate other models.
      hWght = massWeightToggle * &
              max(0., -bathyT(i,j)-e(i,j+1,K), -bathyT(i,j+1)-e(i,j,K))
      if (hWght > 0.) then
        hL = (e(i,j,K) - e(i,j,K+1)) + dz_subroundoff
        hR = (e(i,j+1,K) - e(i,j+1,K+1)) + dz_subroundoff
        hWght = hWght * ( (hL-hR)/(hL+hR) )**2
        iDenom = 1./( hWght*(hR + hL) + hL*hR )
        Ttl = ( (hWght*hR)*T_t(i,j+1,k) + (hWght*hL + hR*hL)*T_t(i,j,k) ) * iDenom
        Ttr = ( (hWght*hL)*T_t(i,j,k) + (hWght*hR + hR*hL)*T_t(i,j+1,k) ) * iDenom
        Tbl = ( (hWght*hR)*T_b(i,j+1,k) + (hWght*hL + hR*hL)*T_b(i,j,k) ) * iDenom
        Tbr = ( (hWght*hL)*T_b(i,j,k) + (hWght*hR + hR*hL)*T_b(i,j+1,k) ) * iDenom
        Stl = ( (hWght*hR)*S_t(i,j+1,k) + (hWght*hL + hR*hL)*S_t(i,j,k) ) * iDenom
        Str = ( (hWght*hL)*S_t(i,j,k) + (hWght*hR + hR*hL)*S_t(i,j+1,k) ) * iDenom
        Sbl = ( (hWght*hR)*S_b(i,j+1,k) + (hWght*hL + hR*hL)*S_b(i,j,k) ) * iDenom
        Sbr = ( (hWght*hL)*S_b(i,j,k) + (hWght*hR + hR*hL)*S_b(i,j+1,k) ) * iDenom
      else
        Ttl = T_t(i,j,k); Tbl = T_b(i,j,k); Ttr = T_t(i,j+1,k); Tbr = T_b(i,j+1,k)
        Stl = S_t(i,j,k); Sbl = S_b(i,j,k); Str = S_t(i,j+1,k); Sbr = S_b(i,j+1,k)
      endif

      do m=2,4
        w_left = wt_t(m) ; w_right = wt_b(m)
        dz_y(m,i) = w_left*(e(i,j,K) - e(i,j,K+1)) + w_right*(e(i,j+1,K) - e(i,j+1,K+1))

        ! Salinity and temperature points are linearly interpolated in
        ! the horizontal. The subscript (1) refers to the top value in
        ! the vertical profile while subscript (5) refers to the bottom
        ! value in the vertical profile.
        pos = i*15+(m-2)*5
        T15(pos+1) = w_left*Ttl + w_right*Ttr
        T15(pos+5) = w_left*Tbl + w_right*Tbr

        S15(pos+1) = w_left*Stl + w_right*Str
        S15(pos+5) = w_left*Sbl + w_right*Sbr

        p15(pos+1) = -GxRho*((w_left*e(i,j,K) + w_right*e(i,j+1,K)) - z0pres)

        ! Pressure
        do n=2,5
          p15(pos+n) = p15(pos+n-1) + GxRho*0.25*dz_y(m,i)
        enddo

        ! Salinity and temperature (linear interpolation in the vertical)
        do n=2,4
          S15(pos+n) = wt_t(n) * S15(pos+1) + wt_b(n) * S15(pos+5)
          T15(pos+n) = wt_t(n) * T15(pos+1) + wt_b(n) * T15(pos+5)
        enddo
        if (use_varT) T215(pos+1:pos+5) = w_left*tv%varT(i,j,k) + w_right*tv%varT(i,j+1,k)
        if (use_covarTS) TS15(pos+1:pos+5) = w_left*tv%covarTS(i,j,k) + w_right*tv%covarTS(i,j+1,k)
        if (use_varS) S215(pos+1:pos+5) = w_left*tv%varS(i,j,k) + w_right*tv%varS(i,j+1,k)
      enddo
    enddo

    if (use_stanley_eos) then
      if (rho_scale /= 1.0) then
        call calculate_density(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                               T215(15*HI%isc+1:), TS15(15*HI%isc+1:), S215(15*HI%isc+1:), &
                               r15(15*HI%isc+1:), 1, 15*(HI%iec-HI%isc+1), EOS, &
                               rho_ref=rho_ref_mks, scale=rho_scale)
      else
        call calculate_density(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                               T215(15*HI%isc+1:), TS15(15*HI%isc+1:), S215(15*HI%isc+1:), &
                               r15(15*HI%isc+1:), 1, 15*(HI%iec-HI%isc+1), EOS, rho_ref=rho_ref_mks)
      endif
    else
      if (use_rho_ref) then
        if (rho_scale /= 1.0) then
          call calculate_density(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                                 r15(15*HI%isc+1:), 1, 15*(HI%iec-HI%isc+1), EOS, &
                                 rho_ref=rho_ref_mks, scale=rho_scale)
        else
          call calculate_density(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                                 r15(15*HI%isc+1:), 1, 15*(HI%iec-HI%isc+1), EOS, rho_ref=rho_ref_mks)
        endif
      else
        if (rho_scale /= 1.0) then
          call calculate_density(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                                 r15(15*HI%isc+1:), 1, 15*(HI%iec-HI%isc+1), EOS, &
                                 scale=rho_scale)
        else
          call calculate_density(T15(15*HI%isc+1:), S15(15*HI%isc+1:), p15(15*HI%isc+1:), &
                                 r15(15*HI%isc+1:), 1, 15*(HI%iec-HI%isc+1), EOS)
        endif
      endif
    endif

    do i=HI%isc,HI%iec
      intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)

      ! Use Boole's rule to estimate the pressure anomaly change.
      if (use_rho_ref) then
        do m = 2,4
          pos = i*15+(m-2)*5
          intz(m) = G_e*dz_y(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + &
                                           32.0*(r15(pos+2)+r15(pos+4)) + &
                                           12.0*r15(pos+3)) )
        enddo
      else
        do m = 2,4
          pos = i*15+(m-2)*5
          intz(m) = G_e*dz_y(m,i)*( C1_90*(7.0*(r15(pos+1)+r15(pos+5)) + &
                                           32.0*(r15(pos+2)+r15(pos+4)) + &
                                           12.0*r15(pos+3)) - rho_ref )
        enddo
      endif
      ! Use Boole's rule to integrate the values.
      inty_dpa(i,J) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    enddo
  enddo ; endif

end subroutine int_density_dz_generic_plm


!> Compute pressure gradient force integrals for layer "k" and the case where T and S
!! are parabolic profiles
subroutine int_density_dz_generic_ppm(k, tv, T_t, T_b, S_t, S_b, e, &
                                      rho_ref, rho_0, G_e, dz_subroundoff, bathyT, HI, GV, EOS, US, &
                                      dpa, intz_dpa, intx_dpa, inty_dpa, useMassWghtInterp, Z_0p)
  integer,              intent(in)  :: k   !< Layer index to calculate integrals for
  type(hor_index_type), intent(in)  :: HI  !< Ocean horizontal index structures for the input arrays
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(thermo_var_ptrs), intent(in) :: tv  !< Thermodynamic variables
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: T_t !< Potential temperature at the cell top [degC]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: T_b !< Potential temperature at the cell bottom [degC]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: S_t !< Salinity at the cell top [ppt]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)), &
                        intent(in)  :: S_b !< Salinity at the cell bottom [ppt]
  real, dimension(SZI_(HI),SZJ_(HI),SZK_(GV)+1), &
                        intent(in)  :: e   !< Height of interfaces [Z ~> m]
  real,                 intent(in)  :: rho_ref !< A mean density [R ~> kg m-3] or [kg m-3], that is
                                           !! subtracted out to reduce the magnitude of each of the integrals.
  real,                 intent(in)  :: rho_0 !< A density [R ~> kg m-3] or [kg m-3], that is used to calculate
                                           !! the pressure (as p~=-z*rho_0*G_e) used in the equation of state.
  real,                 intent(in)  :: G_e !< The Earth's gravitational acceleration [m s-2]
  real,                 intent(in)  :: dz_subroundoff !< A minuscule thickness change [Z ~> m]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: bathyT !< The depth of the bathymetry [Z ~> m]
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dpa !< The change in the pressure anomaly across the layer [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intz_dpa !< The integral through the thickness of the layer of
                                           !! the pressure anomaly relative to the anomaly at the
                                           !! top of the layer [R L2 Z T-2 ~> Pa m]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dpa !< The integral in x of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the x grid spacing [R L2 T-2 ~> Pa]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dpa !< The integral in y of the difference between the
                                           !! pressure anomaly at the top and bottom of the layer
                                           !! divided by the y grid spacing [R L2 T-2 ~> Pa]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting to
                                           !! interpolate T/S for top and bottom integrals.
  real,       optional, intent(in)  :: Z_0p !< The height at which the pressure is 0 [Z ~> m]

! This subroutine calculates (by numerical quadrature) integrals of
! pressure anomalies across layers, which are required for calculating the
! finite-volume form pressure accelerations in a Boussinesq model.  The one
! potentially dodgy assumption here is that rho_0 is used both in the denominator
! of the accelerations, and in the pressure used to calculated density (the
! latter being -z*rho_0*G_e).  These two uses could be separated if need be.
!
! It is assumed that the salinity and temperature profiles are parabolic in the
! vertical. The top and bottom values within each layer are provided and
! a parabolic interpolation is used to compute intermediate values.

  ! Local variables
  real :: T5(5) ! Temperatures along a line of subgrid locations [degC]
  real :: S5(5) ! Salinities along a line of subgrid locations [ppt]
  real :: T25(5) ! SGS temperature variance along a line of subgrid locations [degC2]
  real :: TS5(5) ! SGS temperature-salinity covariance along a line of subgrid locations [degC ppt]
  real :: S25(5) ! SGS salinity variance along a line of subgrid locations [ppt2]
  real :: p5(5) ! Pressures at five quadrature points, never rescaled from Pa [Pa]
  real :: r5(5) ! Density anomalies from rho_ref at quadrature points [R ~> kg m-3] or [kg m-3]
  real :: wt_t(5), wt_b(5) ! Top and bottom weights [nondim]
  real :: rho_anom ! The integrated density anomaly [R ~> kg m-3] or [kg m-3]
  real :: w_left, w_right  ! Left and right weights [nondim]
  real :: intz(5) ! The gravitational acceleration times the integrals of density
                  ! with height at the 5 sub-column locations [R L2 T-2 ~> Pa] or [Pa]
  real, parameter :: C1_90 = 1.0/90.0  ! Rational constants.
  real :: GxRho ! The gravitational acceleration times density and unit conversion factors [Pa Z-1 ~> kg m-2 s-2]
  real :: I_Rho ! The inverse of the Boussinesq density [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: rho_scale ! A scaling factor for densities from kg m-3 to R [R m3 kg-1 ~> 1]
  real :: rho_ref_mks ! The reference density in MKS units, never rescaled from kg m-3 [kg m-3]
  real :: dz ! Layer thicknesses at tracer points [Z ~> m]
  real :: massWeightToggle ! A non-dimensional toggle factor (0 or 1) [nondim]
  real :: Ttl, Tbl, Tml, Ttr, Tbr, Tmr ! Temperatures at the velocity cell corners [degC]
  real :: Stl, Sbl, Sml, Str, Sbr, Smr ! Salinities at the velocity cell corners [ppt]
  real :: s6 ! PPM curvature coefficient for S [ppt]
  real :: t6 ! PPM curvature coefficient for T [degC]
  real :: T_top, T_mn, T_bot ! Left edge, cell mean and right edge values used in PPM reconstructions of T
  real :: S_top, S_mn, S_bot ! Left edge, cell mean and right edge values used in PPM reconstructions of S
  real :: z0pres ! The height at which the pressure is zero [Z ~> m]
  real :: hWght  ! A topographically limited thicknes weight [Z ~> m]
  real :: hL, hR ! Thicknesses to the left and right [Z ~> m]
  real :: iDenom ! The denominator of the thickness weight expressions [Z-2 ~> m-2]
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n
  logical :: use_PPM ! If false, assume zero curvature in reconstruction, i.e. PLM
  logical :: use_stanley_eos ! True is SGS variance fields exist in tv.
  logical :: use_varT, use_varS, use_covarTS

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB

  rho_scale = US%kg_m3_to_R
  GxRho = US%RL2_T2_to_Pa * G_e * rho_0
  rho_ref_mks = rho_ref * US%R_to_kg_m3
  I_Rho = 1.0 / rho_0
  z0pres = 0.0 ; if (present(Z_0p)) z0pres = Z_0p
  massWeightToggle = 0.
  if (present(useMassWghtInterp)) then
    if (useMassWghtInterp) massWeightToggle = 1.
  endif

  ! In event PPM calculation is bypassed with use_PPM=False
  s6 = 0.
  t6 = 0.
  use_PPM = .true. ! This is a place-holder to allow later re-use of this function

  use_varT = associated(tv%varT)
  use_covarTS = associated(tv%covarTS)
  use_varS = associated(tv%varS)
  use_stanley_eos = use_varT .or. use_covarTS .or. use_varS
  T25(:) = 0.
  TS5(:) = 0.
  S25(:) = 0.

  do n = 1, 5
    wt_t(n) = 0.25 * real(5-n)
    wt_b(n) = 1.0 - wt_t(n)
  enddo

  ! 1. Compute vertical integrals
  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    if (use_PPM) then
      ! Curvature coefficient of the parabolas
      s6 = 3.0 * ( 2.0*tv%S(i,j,k) - ( S_t(i,j,k) + S_b(i,j,k) ) )
      t6 = 3.0 * ( 2.0*tv%T(i,j,k) - ( T_t(i,j,k) + T_b(i,j,k) ) )
    endif
    dz = e(i,j,K) - e(i,j,K+1)
    do n=1,5
      p5(n) = -GxRho*((e(i,j,K) - z0pres) - 0.25*real(n-1)*dz)
      ! Salinity and temperature points are reconstructed with PPM
      S5(n) = wt_t(n) * S_t(i,j,k) + wt_b(n) * ( S_b(i,j,k) + s6 * wt_t(n) )
      T5(n) = wt_t(n) * T_t(i,j,k) + wt_b(n) * ( T_b(i,j,k) + t6 * wt_t(n) )
    enddo
    if (use_stanley_eos) then
      if (use_varT) T25(:) = tv%varT(i,j,k)
      if (use_covarTS) TS5(:) = tv%covarTS(i,j,k)
      if (use_varS) S25(:) = tv%varS(i,j,k)
      call calculate_density(T5, S5, p5, T25, TS5, S25, r5, &
                             1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
    else
      call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
    endif

    ! Use Boole's rule to estimate the pressure anomaly change.
    rho_anom = C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3))
    dpa(i,j) = G_e*dz*rho_anom
    if (present(intz_dpa)) then
      ! Use a Boole's-rule-like fifth-order accurate estimate of
      ! the double integral of the pressure anomaly.
      intz_dpa(i,j) = 0.5*G_e*dz**2 * &
                      (rho_anom - C1_90*(16.0*(r5(4)-r5(2)) + 7.0*(r5(5)-r5(1))) )
    endif
  enddo ; enddo ! end loops on j and i

  ! 2. Compute horizontal integrals in the x direction
  if (present(intx_dpa)) then ; do j=HI%jsc,HI%jec ; do I=Isq,Ieq
    ! Corner values of T and S
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, almost like thickness
    ! weighting.
    ! Note: To work in terrain following coordinates we could offset
    ! this distance by the layer thickness to replicate other models.
    hWght = massWeightToggle * &
            max(0., -bathyT(i,j)-e(i+1,j,K), -bathyT(i+1,j)-e(i,j,K))
    if (hWght > 0.) then
      hL = (e(i,j,K) - e(i,j,K+1)) + dz_subroundoff
      hR = (e(i+1,j,K) - e(i+1,j,K+1)) + dz_subroundoff
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1./( hWght*(hR + hL) + hL*hR )
      Ttl = ( (hWght*hR)*T_t(i+1,j,k) + (hWght*hL + hR*hL)*T_t(i,j,k) ) * iDenom
      Tbl = ( (hWght*hR)*T_b(i+1,j,k) + (hWght*hL + hR*hL)*T_b(i,j,k) ) * iDenom
      Tml = ( (hWght*hR)*tv%T(i+1,j,k)+ (hWght*hL + hR*hL)*tv%T(i,j,k) ) * iDenom
      Ttr = ( (hWght*hL)*T_t(i,j,k) + (hWght*hR + hR*hL)*T_t(i+1,j,k) ) * iDenom
      Tbr = ( (hWght*hL)*T_b(i,j,k) + (hWght*hR + hR*hL)*T_b(i+1,j,k) ) * iDenom
      Tmr = ( (hWght*hL)*tv%T(i,j,k) + (hWght*hR + hR*hL)*tv%T(i+1,j,k) ) * iDenom
      Stl = ( (hWght*hR)*S_t(i+1,j,k) + (hWght*hL + hR*hL)*S_t(i,j,k) ) * iDenom
      Sbl = ( (hWght*hR)*S_b(i+1,j,k) + (hWght*hL + hR*hL)*S_b(i,j,k) ) * iDenom
      Sml = ( (hWght*hR)*tv%S(i+1,j,k) + (hWght*hL + hR*hL)*tv%S(i,j,k) ) * iDenom
      Str = ( (hWght*hL)*S_t(i,j,k) + (hWght*hR + hR*hL)*S_t(i+1,j,k) ) * iDenom
      Sbr = ( (hWght*hL)*S_b(i,j,k) + (hWght*hR + hR*hL)*S_b(i+1,j,k) ) * iDenom
      Smr = ( (hWght*hL)*tv%S(i,j,k) + (hWght*hR + hR*hL)*tv%S(i+1,j,k) ) * iDenom
    else
      Ttl = T_t(i,j,k); Tbl = T_b(i,j,k); Ttr = T_t(i+1,j,k); Tbr = T_b(i+1,j,k)
      Tml = tv%T(i,j,k); Tmr = tv%T(i+1,j,k)
      Stl = S_t(i,j,k); Sbl = S_b(i,j,k); Str = S_t(i+1,j,k); Sbr = S_b(i+1,j,k)
      Sml = tv%S(i,j,k); Smr = tv%S(i+1,j,k)
    endif

    do m=2,4
      w_left = wt_t(m) ; w_right = wt_b(m)

      ! Salinity and temperature points are linearly interpolated in
      ! the horizontal. The subscript (1) refers to the top value in
      ! the vertical profile while subscript (5) refers to the bottom
      ! value in the vertical profile.
      T_top = w_left*Ttl + w_right*Ttr
      T_mn = w_left*Tml + w_right*Tmr
      T_bot = w_left*Tbl + w_right*Tbr

      S_top = w_left*Stl + w_right*Str
      S_mn = w_left*Sml + w_right*Smr
      S_bot = w_left*Sbl + w_right*Sbr

      ! Pressure
      dz = w_left*(e(i,j,K) - e(i,j,K+1)) + w_right*(e(i+1,j,K) - e(i+1,j,K+1))
      p5(1) = -GxRho*((w_left*e(i,j,K) + w_right*e(i+1,j,K)) - z0pres)
      do n=2,5
        p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo

      ! Parabolic reconstructions in the vertical for T and S
      if (use_PPM) then
        ! Coefficients of the parabolas
        s6 = 3.0 * ( 2.0*S_mn - ( S_top + S_bot ) )
        t6 = 3.0 * ( 2.0*T_mn - ( T_top + T_bot ) )
      endif
      do n=1,5
        S5(n) = wt_t(n) * S_top + wt_b(n) * ( S_bot + s6 * wt_t(n) )
        T5(n) = wt_t(n) * T_top + wt_b(n) * ( T_bot + t6 * wt_t(n) )
      enddo

      if (use_stanley_eos) then
        if (use_varT) T25(:) = w_left*tv%varT(i,j,k) + w_right*tv%varT(i+1,j,k)
        if (use_covarTS) TS5(:) = w_left*tv%covarTS(i,j,k) + w_right*tv%covarTS(i+1,j,k)
        if (use_varS) S25(:) = w_left*tv%varS(i,j,k) + w_right*tv%varS(i+1,j,k)
        call calculate_density(T5, S5, p5, T25, TS5, S25, r5, &
                               1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
      else
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
      endif

      ! Use Boole's rule to estimate the pressure anomaly change.
      intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)) )
    enddo ! m
    intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)

    ! Use Boole's rule to integrate the bottom pressure anomaly values in x.
    intx_dpa(I,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + 12.0*intz(3))

  enddo ; enddo ; endif

  ! 3. Compute horizontal integrals in the y direction
  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=HI%isc,HI%iec
    ! Corner values of T and S
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, almost like thickness
    ! weighting.
    ! Note: To work in terrain following coordinates we could offset
    ! this distance by the layer thickness to replicate other models.
    hWght = massWeightToggle * &
            max(0., -bathyT(i,j)-e(i,j+1,K), -bathyT(i,j+1)-e(i,j,K))
    if (hWght > 0.) then
      hL = (e(i,j,K) - e(i,j,K+1)) + dz_subroundoff
      hR = (e(i,j+1,K) - e(i,j+1,K+1)) + dz_subroundoff
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1./( hWght*(hR + hL) + hL*hR )
      Ttl = ( (hWght*hR)*T_t(i,j+1,k) + (hWght*hL + hR*hL)*T_t(i,j,k) ) * iDenom
      Tbl = ( (hWght*hR)*T_b(i,j+1,k) + (hWght*hL + hR*hL)*T_b(i,j,k) ) * iDenom
      Tml = ( (hWght*hR)*tv%T(i,j+1,k)+ (hWght*hL + hR*hL)*tv%T(i,j,k) ) * iDenom
      Ttr = ( (hWght*hL)*T_t(i,j,k) + (hWght*hR + hR*hL)*T_t(i,j+1,k) ) * iDenom
      Tbr = ( (hWght*hL)*T_b(i,j,k) + (hWght*hR + hR*hL)*T_b(i,j+1,k) ) * iDenom
      Tmr = ( (hWght*hL)*tv%T(i,j,k) + (hWght*hR + hR*hL)*tv%T(i,j+1,k) ) * iDenom
      Stl = ( (hWght*hR)*S_t(i,j+1,k) + (hWght*hL + hR*hL)*S_t(i,j,k) ) * iDenom
      Sbl = ( (hWght*hR)*S_b(i,j+1,k) + (hWght*hL + hR*hL)*S_b(i,j,k) ) * iDenom
      Sml = ( (hWght*hR)*tv%S(i,j+1,k)+ (hWght*hL + hR*hL)*tv%S(i,j,k) ) * iDenom
      Str = ( (hWght*hL)*S_t(i,j,k) + (hWght*hR + hR*hL)*S_t(i,j+1,k) ) * iDenom
      Sbr = ( (hWght*hL)*S_b(i,j,k) + (hWght*hR + hR*hL)*S_b(i,j+1,k) ) * iDenom
      Smr = ( (hWght*hL)*tv%S(i,j,k) + (hWght*hR + hR*hL)*tv%S(i,j+1,k) ) * iDenom
    else
      Ttl = T_t(i,j,k); Tbl = T_b(i,j,k); Ttr = T_t(i,j+1,k); Tbr = T_b(i,j+1,k)
      Tml = tv%T(i,j,k); Tmr = tv%T(i,j+1,k)
      Stl = S_t(i,j,k); Sbl = S_b(i,j,k); Str = S_t(i,j+1,k); Sbr = S_b(i,j+1,k)
      Sml = tv%S(i,j,k); Smr = tv%S(i,j+1,k)
    endif

    do m=2,4
      w_left = wt_t(m) ; w_right = wt_b(m)

      ! Salinity and temperature points are linearly interpolated in
      ! the horizontal. The subscript (1) refers to the top value in
      ! the vertical profile while subscript (5) refers to the bottom
      ! value in the vertical profile.
      T_top = w_left*Ttl + w_right*Ttr
      T_mn = w_left*Tml + w_right*Tmr
      T_bot = w_left*Tbl + w_right*Tbr

      S_top = w_left*Stl + w_right*Str
      S_mn = w_left*Sml + w_right*Smr
      S_bot = w_left*Sbl + w_right*Sbr

      ! Pressure
      dz = w_left*(e(i,j,K) - e(i,j,K+1)) + w_right*(e(i,j+1,K) - e(i,j+1,K+1))
      p5(1) = -GxRho*((w_left*e(i,j,K) + w_right*e(i,j+1,K)) - z0pres)
      do n=2,5
        p5(n) = p5(n-1) + GxRho*0.25*dz
      enddo

      ! Parabolic reconstructions in the vertical for T and S
      if (use_PPM) then
        ! Coefficients of the parabolas
        s6 = 3.0 * ( 2.0*S_mn - ( S_top + S_bot ) )
        t6 = 3.0 * ( 2.0*T_mn - ( T_top + T_bot ) )
      endif
      do n=1,5
        S5(n) = wt_t(n) * S_top + wt_b(n) * ( S_bot + s6 * wt_t(n) )
        T5(n) = wt_t(n) * T_top + wt_b(n) * ( T_bot + t6 * wt_t(n) )
      enddo

      if (use_stanley_eos) then
        if (use_varT) T25(:) = w_left*tv%varT(i,j,k) + w_right*tv%varT(i,j+1,k)
        if (use_covarTS) TS5(:) = w_left*tv%covarTS(i,j,k) + w_right*tv%covarTS(i,j+1,k)
        if (use_varS) S25(:) = w_left*tv%varS(i,j,k) + w_right*tv%varS(i,j+1,k)
        call calculate_density(T5, S5, p5, T25, TS5, S25, r5, &
                               1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
      else
        call calculate_density(T5, S5, p5, r5, 1, 5, EOS, rho_ref=rho_ref_mks, scale=rho_scale)
      endif

      ! Use Boole's rule to estimate the pressure anomaly change.
      intz(m) = G_e*dz*( C1_90*(7.0*(r5(1)+r5(5)) + 32.0*(r5(2)+r5(4)) + 12.0*r5(3)) )
    enddo ! m
    intz(1) = dpa(i,j) ; intz(5) = dpa(i,j+1)

    ! Use Boole's rule to integrate the bottom pressure anomaly values in y.
    inty_dpa(i,J) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + 12.0*intz(3))

  enddo ; enddo ; endif

end subroutine int_density_dz_generic_ppm

!> Calls the appropriate subroutine to calculate analytical and nearly-analytical
!! integrals in pressure across layers of geopotential anomalies, which are
!! required for calculating the finite-volume form pressure accelerations in a
!! non-Boussinesq model.  There are essentially no free assumptions, apart from the
!! use of Boole's rule to do the horizontal integrals, and from a truncation in the
!! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.
subroutine int_specific_vol_dp(T, S, p_t, p_b, alpha_ref, HI, EOS, US, &
                               dza, intp_dza, intx_dza, inty_dza, halo_size, &
                               bathyP, dP_tiny, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI  !< The horizontal index structure
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T   !< Potential temperature referenced to the surface [degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S   !< Salinity [ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_t !< Pressure at the top of the layer [R L2 T-2 ~> Pa] or [Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_b !< Pressure at the bottom of the layer [R L2 T-2 ~> Pa] or [Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but this reduces the effects of roundoff.
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dza !< The change in the geopotential anomaly across
                            !! the layer [L2 T-2 ~> m2 s-2] or [m2 s-2]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of the
                            !! geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2] or [Pa m2 s-2]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dza !< The integral in x of the difference between the
                            !! geopotential anomaly at the top and bottom of the layer divided by
                            !! the x grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dza !< The integral in y of the difference between the
                            !! geopotential anomaly at the top and bottom of the layer divided by
                            !! the y grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate dza.
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: bathyP  !< The pressure at the bathymetry [R L2 T-2 ~> Pa] or [Pa]
  real,       optional, intent(in)  :: dP_tiny !< A minuscule pressure change with
                            !! the same units as p_t [R L2 T-2 ~> Pa] or [Pa]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting
                            !! to interpolate T/S for top and bottom integrals.

  if (EOS_quadrature(EOS)) then
    call int_spec_vol_dp_generic_pcm(T, S, p_t, p_b, alpha_ref, HI, EOS, US, &
                                     dza, intp_dza, intx_dza, inty_dza, halo_size, &
                                     bathyP, dP_tiny, useMassWghtInterp)
  else
    call analytic_int_specific_vol_dp(T, S, p_t, p_b, alpha_ref, HI, EOS, &
                                      dza, intp_dza, intx_dza, inty_dza, halo_size, &
                                      bathyP, dP_tiny, useMassWghtInterp)
  endif

end subroutine int_specific_vol_dp


!>   This subroutine calculates integrals of specific volume anomalies in
!! pressure across layers, which are required for calculating the finite-volume
!! form pressure accelerations in a non-Boussinesq model.  There are essentially
!! no free assumptions, apart from the use of Boole's rule quadrature to do the integrals.
subroutine int_spec_vol_dp_generic_pcm(T, S, p_t, p_b, alpha_ref, HI, EOS, US, dza, &
                                       intp_dza, intx_dza, inty_dza, halo_size, &
                                       bathyP, dP_neglect, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI !< A horizontal index type structure.
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T  !< Potential temperature of the layer [degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S  !< Salinity of the layer [ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_t !< Pressure atop the layer [R L2 T-2 ~> Pa] or [Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_b !< Pressure below the layer [R L2 T-2 ~> Pa] or [Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but alpha_ref alters the effects of roundoff, and
                            !! answers do change.
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dza !< The change in the geopotential anomaly
                            !! across the layer [L2 T-2 ~> m2 s-2] or [m2 s-2]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of
                            !! the geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2] or [Pa m2 s-2]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dza  !< The integral in x of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the x grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dza  !< The integral in y of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the y grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate dza.
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(in)  :: bathyP !< The pressure at the bathymetry [R L2 T-2 ~> Pa] or [Pa]
  real,       optional, intent(in)  :: dP_neglect !< A minuscule pressure change with
                                             !! the same units as p_t [R L2 T-2 ~> Pa] or [Pa]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting
                            !! to interpolate T/S for top and bottom integrals.

!   This subroutine calculates analytical and nearly-analytical integrals in
! pressure across layers of geopotential anomalies, which are required for
! calculating the finite-volume form pressure accelerations in a non-Boussinesq
! model.  There are essentially no free assumptions, apart from the use of
! Boole's rule to do the horizontal integrals, and from a truncation in the
! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.

  ! Local variables
  real :: T5(5)      ! Temperatures at five quadrature points [degC]
  real :: S5(5)      ! Salinities at five quadrature points [ppt]
  real :: p5(5)      ! Pressures at five quadrature points, scaled back to Pa if necessary [Pa]
  real :: a5(5)      ! Specific volumes at five quadrature points [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: alpha_anom ! The depth averaged specific density anomaly [R-1 ~> m3 kg-1]
  real :: dp         ! The pressure change through a layer [R L2 T-2 ~> Pa]
  real :: hWght      ! A pressure-thickness below topography [R L2 T-2 ~> Pa]
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [R L2 T-2 ~> Pa]
  real :: alpha_ref_mks ! The reference specific volume in MKS units, never rescaled from m3 kg-1 [m3 kg-1]
  real :: iDenom     ! The inverse of the denominator in the weights [T4 R-2 L-4 ~> Pa-2]
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim]
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim]
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim]
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim]
  real :: intp(5)    ! The integrals of specific volume with pressure at the
                     ! 5 sub-column locations [L2 T-2 ~> m2 s-2]
  real :: RL2_T2_to_Pa  ! A unit conversion factor from the rescaled units of pressure to Pa [Pa T2 R-1 L-2 ~> 1]
  real :: SV_scale   ! A multiplicative factor by which to scale specific
                     ! volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant.
  integer :: Isq, Ieq, Jsq, Jeq, ish, ieh, jsh, jeh, i, j, m, n, halo

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB
  halo = 0 ; if (present(halo_size)) halo = MAX(halo_size,0)
  ish = HI%isc-halo ; ieh = HI%iec+halo ; jsh = HI%jsc-halo ; jeh = HI%jec+halo
  if (present(intx_dza)) then ; ish = MIN(Isq,ish) ; ieh = MAX(Ieq+1,ieh); endif
  if (present(inty_dza)) then ; jsh = MIN(Jsq,jsh) ; jeh = MAX(Jeq+1,jeh); endif

  SV_scale = US%R_to_kg_m3
  RL2_T2_to_Pa = US%RL2_T2_to_Pa
  alpha_ref_mks = alpha_ref * US%kg_m3_to_R

  do_massWeight = .false.
  if (present(useMassWghtInterp)) then ; if (useMassWghtInterp) then
    do_massWeight = .true.
    if (.not.present(bathyP)) call MOM_error(FATAL, "int_spec_vol_dp_generic: "//&
        "bathyP must be present if useMassWghtInterp is present and true.")
    if (.not.present(dP_neglect)) call MOM_error(FATAL, "int_spec_vol_dp_generic: "//&
        "dP_neglect must be present if useMassWghtInterp is present and true.")
  endif ; endif

  do j=jsh,jeh ; do i=ish,ieh
    dp = p_b(i,j) - p_t(i,j)
    do n=1,5
      T5(n) = T(i,j) ; S5(n) = S(i,j)
      p5(n) = RL2_T2_to_Pa * (p_b(i,j) - 0.25*real(n-1)*dp)
    enddo

    if (SV_scale /= 1.0) then
      call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks, scale=SV_scale)
    else
      call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks)
    endif

    ! Use Boole's rule to estimate the interface height anomaly change.
    alpha_anom = C1_90*(7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4)) + 12.0*a5(3))
    dza(i,j) = dp*alpha_anom
    ! Use a Boole's-rule-like fifth-order accurate estimate of the double integral of
    ! the interface height anomaly.
    if (present(intp_dza)) intp_dza(i,j) = 0.5*dp**2 * &
          (alpha_anom - C1_90*(16.0*(a5(4)-a5(2)) + 7.0*(a5(5)-a5(1))) )
  enddo ; enddo

  if (present(intx_dza)) then ; do j=HI%jsc,HI%jec ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., bathyP(i,j)-p_t(i+1,j), bathyP(i+1,j)-p_t(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i+1,j) - p_t(i+1,j)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i+1,j)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness weighted.
      p5(1) = RL2_T2_to_Pa * (wt_L*p_b(i,j) + wt_R*p_b(i+1,j))
      dp = wt_L*(p_b(i,j) - p_t(i,j)) + wt_R*(p_b(i+1,j) - p_t(i+1,j))
      T5(1) = wtT_L*T(i,j) + wtT_R*T(i+1,j)
      S5(1) = wtT_L*S(i,j) + wtT_R*S(i+1,j)

      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1) ; p5(n) = p5(n-1) - RL2_T2_to_Pa * 0.25*dp
      enddo
      if (SV_scale /= 1.0) then
        call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks, scale=SV_scale)
      else
        call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks)
      endif

    ! Use Boole's rule to estimate the interface height anomaly change.
      intp(m) = dp*( C1_90*(7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4)) + &
                                12.0*a5(3)))
    enddo
    ! Use Boole's rule to integrate the interface height anomaly values in x.
    intx_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

  if (present(inty_dza)) then ; do J=Jsq,Jeq ; do i=HI%isc,HI%iec
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., bathyP(i,j)-p_t(i,j+1), bathyP(i,j+1)-p_t(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i,j+1) - p_t(i,j+1)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i,j+1)
    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness weighted.
      p5(1) = RL2_T2_to_Pa * (wt_L*p_b(i,j) + wt_R*p_b(i,j+1))
      dp = wt_L*(p_b(i,j) - p_t(i,j)) + wt_R*(p_b(i,j+1) - p_t(i,j+1))
      T5(1) = wtT_L*T(i,j) + wtT_R*T(i,j+1)
      S5(1) = wtT_L*S(i,j) + wtT_R*S(i,j+1)
      do n=2,5
        T5(n) = T5(1) ; S5(n) = S5(1) ; p5(n) = RL2_T2_to_Pa * (p5(n-1) - 0.25*dp)
      enddo
      if (SV_scale /= 1.0) then
        call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks, scale=SV_scale)
      else
        call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks)
      endif

    ! Use Boole's rule to estimate the interface height anomaly change.
      intp(m) = dp*( C1_90*(7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4)) + &
                                12.0*a5(3)))
    enddo
    ! Use Boole's rule to integrate the interface height anomaly values in y.
    inty_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

end subroutine int_spec_vol_dp_generic_pcm

!>   This subroutine calculates integrals of specific volume anomalies in
!! pressure across layers, which are required for calculating the finite-volume
!! form pressure accelerations in a non-Boussinesq model.  There are essentially
!! no free assumptions, apart from the use of Boole's rule quadrature to do the integrals.
subroutine int_spec_vol_dp_generic_plm(T_t, T_b, S_t, S_b, p_t, p_b, alpha_ref, &
                             dP_neglect, bathyP, HI, EOS, US, dza, &
                             intp_dza, intx_dza, inty_dza, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI !< A horizontal index type structure.
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T_t  !< Potential temperature at the top of the layer [degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: T_b  !< Potential temperature at the bottom of the layer [degC]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S_t  !< Salinity at the top the layer [ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: S_b  !< Salinity at the bottom the layer [ppt]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_t !< Pressure atop the layer [R L2 T-2 ~> Pa] or [Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: p_b !< Pressure below the layer [R L2 T-2 ~> Pa] or [Pa]
  real,                 intent(in)  :: alpha_ref !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1]
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but alpha_ref alters the effects of roundoff, and
                            !! answers do change.
  real,                 intent(in)  :: dP_neglect !<!< A miniscule pressure change with
                                             !! the same units as p_t [R L2 T-2 ~> Pa] or [Pa]
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(in)  :: bathyP !< The pressure at the bathymetry [R L2 T-2 ~> Pa] or [Pa]
  type(EOS_type),       pointer     :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US !< A dimensional unit scaling type
  real, dimension(SZI_(HI),SZJ_(HI)), &
                        intent(inout) :: dza !< The change in the geopotential anomaly
                            !! across the layer [L2 T-2 ~> m2 s-2]
  real, dimension(SZI_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intp_dza !< The integral in pressure through the layer of
                            !! the geopotential anomaly relative to the anomaly at the bottom of the
                            !! layer [R L4 T-4 ~> Pa m2 s-2] or [Pa m2 s-2]
  real, dimension(SZIB_(HI),SZJ_(HI)), &
              optional, intent(inout) :: intx_dza  !< The integral in x of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the x grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  real, dimension(SZI_(HI),SZJB_(HI)), &
              optional, intent(inout) :: inty_dza  !< The integral in y of the difference between
                            !! the geopotential anomaly at the top and bottom of the layer divided
                            !! by the y grid spacing [L2 T-2 ~> m2 s-2] or [m2 s-2]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting
                            !! to interpolate T/S for top and bottom integrals.

!   This subroutine calculates analytical and nearly-analytical integrals in
! pressure across layers of geopotential anomalies, which are required for
! calculating the finite-volume form pressure accelerations in a non-Boussinesq
! model.  There are essentially no free assumptions, apart from the use of
! Boole's rule to do the horizontal integrals, and from a truncation in the
! series for log(1-eps/1+eps) that assumes that |eps| < 0.34.

  real :: T5(5)      ! Temperatures at five quadrature points [degC]
  real :: S5(5)      ! Salinities at five quadrature points [ppt]
  real :: p5(5)      ! Pressures at five quadrature points, scaled back to Pa as necessary [Pa]
  real :: a5(5)      ! Specific volumes at five quadrature points [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: T15(15)    ! Temperatures at fifteen interior quadrature points [degC]
  real :: S15(15)    ! Salinities at fifteen interior quadrature points [ppt]
  real :: p15(15)    ! Pressures at fifteen quadrature points, scaled back to Pa as necessary [Pa]
  real :: a15(15)    ! Specific volumes at fifteen quadrature points [R-1 ~> m3 kg-1] or [m3 kg-1]
  real :: wt_t(5), wt_b(5) ! Weights of top and bottom values at quadrature points [nondim]
  real :: T_top, T_bot, S_top, S_bot, P_top, P_bot

  real :: alpha_anom ! The depth averaged specific density anomaly [m3 kg-1]
  real :: dp         ! The pressure change through a layer [R L2 T-2 ~> Pa]
  real :: dp_90(2:4) ! The pressure change through a layer divided by 90 [R L2 T-2 ~> Pa]
  real :: hWght      ! A pressure-thickness below topography [R L2 T-2 ~> Pa]
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [R L2 T-2 ~> Pa]
  real :: alpha_ref_mks ! The reference specific volume in MKS units, never rescaled from m3 kg-1 [m3 kg-1]
  real :: iDenom     ! The inverse of the denominator in the weights [T4 R-2 L-4 ~> Pa-2]
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim]
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim]
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim]
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim]
  real :: intp(5)    ! The integrals of specific volume with pressure at the
                     ! 5 sub-column locations [L2 T-2 ~> m2 s-2]
  real :: RL2_T2_to_Pa  ! A unit conversion factor from the rescaled units of pressure to Pa [Pa T2 R-1 L-2 ~> 1]
  real :: SV_scale   ! A multiplicative factor by which to scale specific
                     ! volume from m3 kg-1 to the desired units [kg m-3 R-1 ~> 1]
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant.
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  integer :: Isq, Ieq, Jsq, Jeq, i, j, m, n, pos

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB

  do_massWeight = .false.
  if (present(useMassWghtInterp)) do_massWeight = useMassWghtInterp

  SV_scale = US%R_to_kg_m3
  RL2_T2_to_Pa = US%RL2_T2_to_Pa
  alpha_ref_mks = alpha_ref * US%kg_m3_to_R

  do n = 1, 5 ! Note that these are reversed from int_density_dz.
    wt_t(n) = 0.25 * real(n-1)
    wt_b(n) = 1.0 - wt_t(n)
  enddo

  ! 1. Compute vertical integrals
  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    dp = p_b(i,j) - p_t(i,j)
    do n=1,5 ! T, S and p are linearly interpolated in the vertical.
      p5(n) = RL2_T2_to_Pa * (wt_t(n) * p_t(i,j) + wt_b(n) * p_b(i,j))
      S5(n) = wt_t(n) * S_t(i,j) + wt_b(n) * S_b(i,j)
      T5(n) = wt_t(n) * T_t(i,j) + wt_b(n) * T_b(i,j)
    enddo
    if (SV_scale /= 1.0) then
      call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks, scale=SV_scale)
    else
      call calculate_spec_vol(T5, S5, p5, a5, 1, 5, EOS, alpha_ref_mks)
    endif

    ! Use Boole's rule to estimate the interface height anomaly change.
    alpha_anom = C1_90*((7.0*(a5(1)+a5(5)) + 32.0*(a5(2)+a5(4))) + 12.0*a5(3))
    dza(i,j) = dp*alpha_anom
    ! Use a Boole's-rule-like fifth-order accurate estimate of the double integral of
    ! the interface height anomaly.
    if (present(intp_dza)) intp_dza(i,j) = 0.5*dp**2 * &
          (alpha_anom - C1_90*(16.0*(a5(4)-a5(2)) + 7.0*(a5(5)-a5(1))) )
  enddo ; enddo

  ! 2. Compute horizontal integrals in the x direction
  if (present(intx_dza)) then ; do j=HI%jsc,HI%jec ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, almost like thickness
    ! weighting. Note: To work in terrain following coordinates we could
    ! offset this distance by the layer thickness to replicate other models.
    hWght = 0.0
    if (do_massWeight) &
      hWght =  max(0., bathyP(i,j)-p_t(i+1,j), bathyP(i+1,j)-p_t(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i+1,j) - p_t(i+1,j)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness weighted.
      P_top = wt_L*p_t(i,j) + wt_R*p_t(i+1,j)
      P_bot = wt_L*p_b(i,j) + wt_R*p_b(i+1,j)
      T_top = wtT_L*T_t(i,j) + wtT_R*T_t(i+1,j)
      T_bot = wtT_L*T_b(i,j) + wtT_R*T_b(i+1,j)
      S_top = wtT_L*S_t(i,j) + wtT_R*S_t(i+1,j)
      S_bot = wtT_L*S_b(i,j) + wtT_R*S_b(i+1,j)
      dp_90(m) = C1_90*(P_bot - P_top)

      ! Salinity, temperature and pressure with linear interpolation in the vertical.
      pos = (m-2)*5
      do n=1,5
        p15(pos+n) = RL2_T2_to_Pa * (wt_t(n) * P_top + wt_b(n) * P_bot)
        S15(pos+n) = wt_t(n) * S_top + wt_b(n) * S_bot
        T15(pos+n) = wt_t(n) * T_top + wt_b(n) * T_bot
      enddo
    enddo

    if (SV_scale /= 1.0) then
      call calculate_spec_vol(T15, S15, p15, a15, 1, 15, EOS, alpha_ref_mks, scale=SV_scale)
    else
      call calculate_spec_vol(T15, S15, p15, a15, 1, 15, EOS, alpha_ref_mks)
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i+1,j)
    do m=2,4
      ! Use Boole's rule to estimate the interface height anomaly change.
      ! The integrals at the ends of the segment are already known.
      pos = (m-2)*5
      intp(m) = dp_90(m)*((7.0*(a15(pos+1)+a15(pos+5)) + &
                          32.0*(a15(pos+2)+a15(pos+4))) + 12.0*a15(pos+3))
    enddo
    ! Use Boole's rule to integrate the interface height anomaly values in x.
    intx_dza(I,j) = C1_90*((7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4))) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

  ! 3. Compute horizontal integrals in the y direction
  if (present(inty_dza)) then ; do J=Jsq,Jeq ; do i=HI%isc,HI%iec
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation
    ! of T,S along the top and bottom integrals, like thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., bathyP(i,j)-p_t(i,j+1), bathyP(i,j+1)-p_t(i,j))
    if (hWght > 0.) then
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i,j+1) - p_t(i,j+1)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom
    else
      hWt_LL = 1.0 ; hWt_LR = 0.0 ; hWt_RR = 1.0 ; hWt_RL = 0.0
    endif

    do m=2,4
      wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
      wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

      ! T, S, and p are interpolated in the horizontal.  The p interpolation
      ! is linear, but for T and S it may be thickness weighted.
      P_top = wt_L*p_t(i,j) + wt_R*p_t(i,j+1)
      P_bot = wt_L*p_b(i,j) + wt_R*p_b(i,j+1)
      T_top = wtT_L*T_t(i,j) + wtT_R*T_t(i,j+1)
      T_bot = wtT_L*T_b(i,j) + wtT_R*T_b(i,j+1)
      S_top = wtT_L*S_t(i,j) + wtT_R*S_t(i,j+1)
      S_bot = wtT_L*S_b(i,j) + wtT_R*S_b(i,j+1)
      dp_90(m) = C1_90*(P_bot - P_top)

      ! Salinity, temperature and pressure with linear interpolation in the vertical.
      pos = (m-2)*5
      do n=1,5
        p15(pos+n) = RL2_T2_to_Pa * (wt_t(n) * P_top + wt_b(n) * P_bot)
        S15(pos+n) = wt_t(n) * S_top + wt_b(n) * S_bot
        T15(pos+n) = wt_t(n) * T_top + wt_b(n) * T_bot
      enddo
    enddo

    if (SV_scale /= 1.0) then
      call calculate_spec_vol(T15, S15, p15, a15, 1, 15, EOS, alpha_ref_mks, scale=SV_scale)
    else
      call calculate_spec_vol(T15, S15, p15, a15, 1, 15, EOS, alpha_ref_mks)
    endif

    intp(1) = dza(i,j) ; intp(5) = dza(i,j+1)
    do m=2,4
      ! Use Boole's rule to estimate the interface height anomaly change.
      ! The integrals at the ends of the segment are already known.
      pos = (m-2)*5
      intp(m) = dp_90(m) * ((7.0*(a15(pos+1)+a15(pos+5)) + &
                            32.0*(a15(pos+2)+a15(pos+4))) + 12.0*a15(pos+3))
    enddo
    ! Use Boole's rule to integrate the interface height anomaly values in x.
    inty_dza(i,J) = C1_90*((7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4))) + &
                           12.0*intp(3))
  enddo ; enddo ; endif

end subroutine int_spec_vol_dp_generic_plm


!> Find the depth at which the reconstructed pressure matches P_tgt
subroutine find_depth_of_pressure_in_cell(T_t, T_b, S_t, S_b, z_t, z_b, P_t, P_tgt, &
                       rho_ref, G_e, EOS, US, P_b, z_out, z_tol)
  real,                  intent(in)  :: T_t !< Potential temperature at the cell top [degC]
  real,                  intent(in)  :: T_b !< Potential temperature at the cell bottom [degC]
  real,                  intent(in)  :: S_t !< Salinity at the cell top [ppt]
  real,                  intent(in)  :: S_b !< Salinity at the cell bottom [ppt]
  real,                  intent(in)  :: z_t !< Absolute height of top of cell [Z ~> m]   (Boussinesq ????)
  real,                  intent(in)  :: z_b !< Absolute height of bottom of cell [Z ~> m]
  real,                  intent(in)  :: P_t !< Anomalous pressure of top of cell, relative
                                            !! to g*rho_ref*z_t [R L2 T-2 ~> Pa]
  real,                  intent(in)  :: P_tgt !< Target pressure at height z_out, relative
                                            !! to g*rho_ref*z_out [R L2 T-2 ~> Pa]
  real,                  intent(in)  :: rho_ref !< Reference density with which calculation
                                            !! are anomalous to [R ~> kg m-3]
  real,                  intent(in)  :: G_e !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  type(EOS_type),        pointer     :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in)  :: US !< A dimensional unit scaling type
  real,                  intent(out) :: P_b !< Pressure at the bottom of the cell [R L2 T-2 ~> Pa]
  real,                  intent(out) :: z_out !< Absolute depth at which anomalous pressure = p_tgt [Z ~> m]
  real, optional,        intent(in)  :: z_tol !< The tolerance in finding z_out [Z ~> m]

  ! Local variables
  real :: dp    ! Pressure thickness of the layer [R L2 T-2 ~> Pa]
  real :: F_guess, F_l, F_r  ! Fractional positions [nondim]
  real :: GxRho ! The product of the gravitational acceleration and reference density [R L2 Z-1 T-2 ~> Pa m-1]
  real :: Pa, Pa_left, Pa_right, Pa_tol ! Pressure anomalies, P = integral of g*(rho-rho_ref) dz [R L2 T-2 ~> Pa]
  character(len=240) :: msg

  GxRho = G_e * rho_ref

  ! Anomalous pressure difference across whole cell
  dp = frac_dp_at_pos(T_t, T_b, S_t, S_b, z_t, z_b, rho_ref, G_e, 1.0, EOS)

  P_b = P_t + dp ! Anomalous pressure at bottom of cell

  if (P_tgt <= P_t ) then
    z_out = z_t
    return
  endif

  if (P_tgt >= P_b) then
    z_out = z_b
    return
  endif

  F_l = 0.
  Pa_left = P_t - P_tgt ! Pa_left < 0
  F_r = 1.
  Pa_right = P_b - P_tgt ! Pa_right > 0
  Pa_tol = GxRho * 1.0e-5*US%m_to_Z
  if (present(z_tol)) Pa_tol = GxRho * z_tol

  F_guess = F_l - Pa_left / (Pa_right - Pa_left) * (F_r - F_l)
  Pa = Pa_right - Pa_left ! To get into iterative loop
  do while ( abs(Pa) > Pa_tol )

    z_out = z_t + ( z_b - z_t ) * F_guess
    Pa = frac_dp_at_pos(T_t, T_b, S_t, S_b, z_t, z_b, rho_ref, G_e, F_guess, EOS) - ( P_tgt - P_t )

    if (Pa<Pa_left) then
      write(msg,*) Pa_left,Pa,Pa_right,P_t-P_tgt,P_b-P_tgt
      call MOM_error(FATAL, 'find_depth_of_pressure_in_cell out of bounds negative: /n'//msg)
    elseif (Pa<0.) then
      Pa_left = Pa
      F_l = F_guess
    elseif (Pa>Pa_right) then
      write(msg,*) Pa_left,Pa,Pa_right,P_t-P_tgt,P_b-P_tgt
      call MOM_error(FATAL, 'find_depth_of_pressure_in_cell out of bounds positive: /n'//msg)
    elseif (Pa>0.) then
      Pa_right = Pa
      F_r = F_guess
    else ! Pa == 0
      return
    endif
    F_guess = F_l - Pa_left / (Pa_right - Pa_left) * (F_r - F_l)

  enddo

end subroutine find_depth_of_pressure_in_cell


!> Returns change in anomalous pressure change from top to non-dimensional
!! position pos between z_t and z_b
real function frac_dp_at_pos(T_t, T_b, S_t, S_b, z_t, z_b, rho_ref, G_e, pos, EOS)
  real,           intent(in)  :: T_t !< Potential temperature at the cell top [degC]
  real,           intent(in)  :: T_b !< Potential temperature at the cell bottom [degC]
  real,           intent(in)  :: S_t !< Salinity at the cell top [ppt]
  real,           intent(in)  :: S_b !< Salinity at the cell bottom [ppt]
  real,           intent(in)  :: z_t !< The geometric height at the top of the layer [Z ~> m]
  real,           intent(in)  :: z_b !< The geometric height at the bottom of the layer [Z ~> m]
  real,           intent(in)  :: rho_ref !< A mean density [R ~> kg m-3], that is subtracted out to
                                     !! reduce the magnitude of each of the integrals.
  real,           intent(in)  :: G_e !< The Earth's gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real,           intent(in)  :: pos !< The fractional vertical position, 0 to 1 [nondim]
  type(EOS_type), pointer     :: EOS !< Equation of state structure
  real                        :: fract_dp_at_pos !< The change in pressure from the layer top to
                                     !! fractional position pos [R L2 T-2 ~> Pa]
  ! Local variables
  real, parameter :: C1_90 = 1.0/90.0  ! A rational constant [nondim]
  real :: dz                 ! Distance from the layer top [Z ~> m]
  real :: top_weight, bottom_weight ! Fractional weights at quadrature points [nondim]
  real :: rho_ave            ! Average density [R ~> kg m-3]
  real, dimension(5) :: T5   ! Temperatures at quadrature points [degC]
  real, dimension(5) :: S5   ! Salinities at quadrature points [ppt]
  real, dimension(5) :: p5   ! Pressures at quadrature points [R L2 T-2 ~> Pa]
  real, dimension(5) :: rho5 ! Densities at quadrature points [R ~> kg m-3]
  integer :: n

  do n=1,5
    ! Evaluate density at five quadrature points
    bottom_weight = 0.25*real(n-1) * pos
    top_weight = 1.0 - bottom_weight
    ! Salinity and temperature points are linearly interpolated
    S5(n) = top_weight * S_t + bottom_weight * S_b
    T5(n) = top_weight * T_t + bottom_weight * T_b
    p5(n) = ( top_weight * z_t + bottom_weight * z_b ) * ( G_e * rho_ref )
  enddo
  call calculate_density(T5, S5, p5, rho5, EOS)
  rho5(:) = rho5(:) !- rho_ref ! Work with anomalies relative to rho_ref

  ! Use Boole's rule to estimate the average density
  rho_ave = C1_90*(7.0*(rho5(1)+rho5(5)) + 32.0*(rho5(2)+rho5(4)) + 12.0*rho5(3))

  dz = ( z_t - z_b ) * pos
  frac_dp_at_pos = G_e * dz * rho_ave
end function frac_dp_at_pos

end module MOM_density_integrals

!> \namespace mom_density_integrals
!!
