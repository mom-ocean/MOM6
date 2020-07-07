!> A simple linear equation of state for sea water with constant coefficients
module MOM_EOS_linear

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_hor_index, only : hor_index_type

implicit none ; private

#include <MOM_memory.h>

public calculate_compress_linear, calculate_density_linear, calculate_spec_vol_linear
public calculate_density_derivs_linear, calculate_density_derivs_scalar_linear
public calculate_specvol_derivs_linear
public calculate_density_scalar_linear, calculate_density_array_linear
public calculate_density_second_derivs_linear
public int_density_dz_linear, int_spec_vol_dp_linear

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Compute the density of sea water (in kg/m^3), or its anomaly from a reference density,
!! using a simple linear equation of state from salinity (in psu), potential temperature (in deg C)
!! and pressure [Pa].
interface calculate_density_linear
  module procedure calculate_density_scalar_linear, calculate_density_array_linear
end interface calculate_density_linear

!> Compute the specific volume of sea water (in m^3/kg), or its anomaly from a reference value,
!! using a simple linear equation of state from salinity (in psu), potential temperature (in deg C)
!! and pressure [Pa].
interface calculate_spec_vol_linear
  module procedure calculate_spec_vol_scalar_linear, calculate_spec_vol_array_linear
end interface calculate_spec_vol_linear

!> For a given thermodynamic state, return the derivatives of density with temperature and
!! salinity using the simple linear equation of state
interface calculate_density_derivs_linear
  module procedure calculate_density_derivs_scalar_linear, calculate_density_derivs_array_linear
end interface calculate_density_derivs_linear

!> For a given thermodynamic state, return the second derivatives of density with various
!! combinations of temperature, salinity, and pressure.  Note that with a simple linear
!! equation of state these second derivatives are all 0.
interface calculate_density_second_derivs_linear
  module procedure calculate_density_second_derivs_scalar_linear, calculate_density_second_derivs_array_linear
end interface calculate_density_second_derivs_linear

contains

!> This subroutine computes the density of sea water with a trivial
!! linear equation of state (in [kg m-3]) from salinity (sal [PSU]),
!! potential temperature (T [degC]), and pressure [Pa].
subroutine calculate_density_scalar_linear(T, S, pressure, rho, &
                                           Rho_T0_S0, dRho_dT, dRho_dS, rho_ref)
  real,           intent(in)  :: T        !< Potential temperature relative to the surface [degC].
  real,           intent(in)  :: S        !< Salinity [PSU].
  real,           intent(in)  :: pressure !< pressure [Pa].
  real,           intent(out) :: rho      !< In situ density [kg m-3].
  real,           intent(in)  :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3].
  real,           intent(in)  :: dRho_dT  !< The derivatives of density with temperature
                                          !! [kg m-3 degC-1].
  real,           intent(in)  :: dRho_dS  !< The derivatives of density with salinity
                                          !! in [kg m-3 ppt-1].
  real, optional, intent(in)  :: rho_ref  !< A reference density [kg m-3].

  if (present(rho_ref)) then
    rho = (Rho_T0_S0 - rho_ref) + (dRho_dT*T + dRho_dS*S)
  else
    rho = Rho_T0_S0 + dRho_dT*T + dRho_dS*S
  endif

end subroutine calculate_density_scalar_linear

!> This subroutine computes the density of sea water with a trivial
!! linear equation of state (in kg/m^3) from salinity (sal in psu),
!! potential temperature (T [degC]), and pressure [Pa].
subroutine calculate_density_array_linear(T, S, pressure, rho, start, npts, &
                                          Rho_T0_S0, dRho_dT, dRho_dS, rho_ref)
  real, dimension(:), intent(in)  :: T        !< potential temperature relative to the surface [degC].
  real, dimension(:), intent(in)  :: S        !< salinity [PSU].
  real, dimension(:), intent(in)  :: pressure !< pressure [Pa].
  real, dimension(:), intent(out) :: rho      !< in situ density [kg m-3].
  integer,            intent(in)  :: start    !< the starting point in the arrays.
  integer,            intent(in)  :: npts     !< the number of values to calculate.
  real,               intent(in)  :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3].
  real,               intent(in)  :: dRho_dT  !< The derivatives of density with temperature
                                              !! [kg m-3 degC-1].
  real,               intent(in)  :: dRho_dS  !< The derivatives of density with salinity
                                              !! in [kg m-3 ppt-1].
  real,     optional, intent(in)  :: rho_ref  !< A reference density [kg m-3].
  ! Local variables
  integer :: j

  if (present(rho_ref)) then ; do j=start,start+npts-1
    rho(j) = (Rho_T0_S0 - rho_ref) + (dRho_dT*T(j) + dRho_dS*S(j))
  enddo ; else ; do j=start,start+npts-1
    rho(j) = Rho_T0_S0 + dRho_dT*T(j) + dRho_dS*S(j)
  enddo ; endif

end subroutine calculate_density_array_linear

!> This subroutine computes the in situ specific volume of sea water (specvol in
!! [m3 kg-1]) from salinity (S [PSU]), potential temperature (T [degC])
!! and pressure [Pa], using a trivial linear equation of state for density.
!! If spv_ref is present, specvol is an anomaly from spv_ref.
subroutine calculate_spec_vol_scalar_linear(T, S, pressure, specvol, &
                                            Rho_T0_S0, dRho_dT, dRho_dS, spv_ref)
  real,    intent(in)  :: T        !< potential temperature relative to the surface
                                   !! [degC].
  real,    intent(in)  :: S        !< Salinity [PSU].
  real,    intent(in)  :: pressure !< Pressure [Pa].
  real,    intent(out) :: specvol  !< In situ specific volume [m3 kg-1].
  real,    intent(in)  :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3].
  real,    intent(in)  :: dRho_dT  !< The derivatives of density with temperature [kg m-3 degC-1].
  real,    intent(in)  :: dRho_dS  !< The derivatives of density with salinity [kg m-3 ppt-1].
  real, optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1].
  ! Local variables
  integer :: j

  if (present(spv_ref)) then
    specvol = ((1.0 - Rho_T0_S0*spv_ref) + spv_ref*(dRho_dT*T + dRho_dS*S)) / &
             ( Rho_T0_S0 + (dRho_dT*T + dRho_dS*S))
  else
    specvol = 1.0 / ( Rho_T0_S0 + (dRho_dT*T + dRho_dS*S))
  endif

end subroutine calculate_spec_vol_scalar_linear

!> This subroutine computes the in situ specific volume of sea water (specvol in
!! [m3 kg-1]) from salinity (S [PSU]), potential temperature (T [degC])
!! and pressure [Pa], using a trivial linear equation of state for density.
!! If spv_ref is present, specvol is an anomaly from spv_ref.
subroutine calculate_spec_vol_array_linear(T, S, pressure, specvol, start, npts, &
                                           Rho_T0_S0, dRho_dT, dRho_dS, spv_ref)
  real, dimension(:), intent(in)  :: T        !< potential temperature relative to the surface
                                              !! [degC].
  real, dimension(:), intent(in)  :: S        !< Salinity [PSU].
  real, dimension(:), intent(in)  :: pressure !< Pressure [Pa].
  real, dimension(:), intent(out) :: specvol  !< in situ specific volume [m3 kg-1].
  integer,            intent(in)  :: start    !< the starting point in the arrays.
  integer,            intent(in)  :: npts     !< the number of values to calculate.
  real,               intent(in)  :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3].
  real,               intent(in)  :: dRho_dT  !< The derivatives of density with temperature [kg m-3 degC-1].
  real,               intent(in)  :: dRho_dS  !< The derivatives of density with salinity [kg m-3 ppt-1].
  real,     optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1].
  ! Local variables
  integer :: j

  if (present(spv_ref)) then ; do j=start,start+npts-1
    specvol(j) = ((1.0 - Rho_T0_S0*spv_ref) + spv_ref*(dRho_dT*T(j) + dRho_dS*S(j))) / &
                 ( Rho_T0_S0 + (dRho_dT*T(j) + dRho_dS*S(j)))
  enddo ; else ; do j=start,start+npts-1
    specvol(j) = 1.0 / ( Rho_T0_S0 + (dRho_dT*T(j) + dRho_dS*S(j)))
  enddo ; endif

end subroutine calculate_spec_vol_array_linear

!> This subroutine calculates the partial derivatives of density    *
!! with potential temperature and salinity.
subroutine calculate_density_derivs_array_linear(T, S, pressure, drho_dT_out, &
                       drho_dS_out, Rho_T0_S0, dRho_dT, dRho_dS, start, npts)
  real,    intent(in),  dimension(:) :: T           !< Potential temperature relative to the surface
                                                    !! [degC].
  real,    intent(in),  dimension(:) :: S           !< Salinity [PSU].
  real,    intent(in),  dimension(:) :: pressure    !< Pressure [Pa].
  real,    intent(out), dimension(:) :: drho_dT_out !< The partial derivative of density with
                                                    !! potential temperature [kg m-3 degC-1].
  real,    intent(out), dimension(:) :: drho_dS_out !< The partial derivative of density with
                                                    !! salinity [kg m-3 ppt-1].
  real,    intent(in)                :: Rho_T0_S0   !< The density at T=0, S=0 [kg m-3].
  real,    intent(in)                :: dRho_dT     !< The derivative of density with temperature [kg m-3 degC-1].
  real,    intent(in)                :: dRho_dS     !< The derivative of density with salinity [kg m-3 ppt-1].
  integer, intent(in)                :: start       !< The starting point in the arrays.
  integer, intent(in)                :: npts        !< The number of values to calculate.
  ! Local variables
  integer :: j

  do j=start,start+npts-1
    drho_dT_out(j) = dRho_dT
    drho_dS_out(j) = dRho_dS
  enddo

end subroutine calculate_density_derivs_array_linear

!> This subroutine calculates the partial derivatives of density    *
!! with potential temperature and salinity for a single point.
subroutine calculate_density_derivs_scalar_linear(T, S, pressure, drho_dT_out, &
                       drho_dS_out, Rho_T0_S0, dRho_dT, dRho_dS)
  real,    intent(in)  :: T           !< Potential temperature relative to the surface
                                      !! [degC].
  real,    intent(in)  :: S           !< Salinity [PSU].
  real,    intent(in)  :: pressure    !< pressure [Pa].
  real,    intent(out) :: drho_dT_out !< The partial derivative of density with
                                      !! potential temperature [kg m-3 degC-1].
  real,    intent(out) :: drho_dS_out !< The partial derivative of density with
                                      !! salinity [kg m-3 ppt-1].
  real,    intent(in)  :: Rho_T0_S0   !< The density at T=0, S=0 [kg m-3].
  real,    intent(in)  :: dRho_dT     !< The derivatives of density with temperature [kg m-3 degC-1].
  real,    intent(in)  :: dRho_dS     !< The derivatives of density with salinity [kg m-3 ppt-1].
  drho_dT_out = dRho_dT
  drho_dS_out = dRho_dS

end subroutine calculate_density_derivs_scalar_linear

!> This subroutine calculates the five, partial second derivatives of density w.r.t.
!! potential temperature and salinity and pressure which for a linear equation of state should all be 0.
subroutine calculate_density_second_derivs_scalar_linear(T, S, pressure, drho_dS_dS, drho_dS_dT, &
                                                         drho_dT_dT, drho_dS_dP, drho_dT_dP)
  real, intent(in)  :: T           !< Potential temperature relative to the surface [degC].
  real, intent(in)  :: S           !< Salinity [PSU].
  real, intent(in)  :: pressure    !< pressure [Pa].
  real, intent(out) :: drho_dS_dS  !< The second derivative of density with
                                   !! salinity [kg m-3 PSU-2].
  real, intent(out) :: drho_dS_dT  !< The second derivative of density with
                                   !! temperature and salinity [kg m-3 ppt-1 degC-1].
  real, intent(out) :: drho_dT_dT  !< The second derivative of density with
                                   !! temperature [kg m-3 degC-2].
  real, intent(out) :: drho_dS_dP  !< The second derivative of density with
                                   !! salinity and pressure [kg m-3 PSU-1 Pa-1].
  real, intent(out) :: drho_dT_dP  !< The second derivative of density with
                                   !! temperature and pressure [kg m-3 degC-1 Pa-1].

  drho_dS_dS = 0.
  drho_dS_dT = 0.
  drho_dT_dT = 0.
  drho_dS_dP = 0.
  drho_dT_dP = 0.

end subroutine calculate_density_second_derivs_scalar_linear

!> This subroutine calculates the five, partial second derivatives of density w.r.t.
!! potential temperature and salinity and pressure which for a linear equation of state should all be 0.
subroutine calculate_density_second_derivs_array_linear(T, S,pressure,  drho_dS_dS, drho_dS_dT, drho_dT_dT,&
                        drho_dS_dP, drho_dT_dP, start, npts)
  real, dimension(:), intent(in)  :: T           !< Potential temperature relative to the surface [degC].
  real, dimension(:), intent(in)  :: S           !< Salinity [PSU].
  real, dimension(:), intent(in)  :: pressure    !< pressure [Pa].
  real, dimension(:), intent(out) :: drho_dS_dS  !< The second derivative of density with
                                                 !! salinity [kg m-3 PSU-2].
  real, dimension(:), intent(out) :: drho_dS_dT  !< The second derivative of density with
                                                 !! temperature and salinity [kg m-3 ppt-1 degC-1].
  real, dimension(:), intent(out) :: drho_dT_dT  !< The second derivative of density with
                                                 !! temperature [kg m-3 degC-2].
  real, dimension(:), intent(out) :: drho_dS_dP  !< The second derivative of density with
                                                 !! salinity and pressure [kg m-3 PSU-1 Pa-1].
  real, dimension(:), intent(out) :: drho_dT_dP  !< The second derivative of density with
                                                 !! temperature and pressure [kg m-3 degC-1 Pa-1].
  integer, intent(in)  :: start       !< The starting point in the arrays.
  integer, intent(in)  :: npts        !< The number of values to calculate.
  ! Local variables
  integer :: j
  do j=start,start+npts-1
    drho_dS_dS(j) = 0.
    drho_dS_dT(j) = 0.
    drho_dT_dT(j) = 0.
    drho_dS_dP(j) = 0.
    drho_dT_dP(j) = 0.
  enddo

end subroutine calculate_density_second_derivs_array_linear

!> Calculate the derivatives of specific volume with temperature and salinity
subroutine calculate_specvol_derivs_linear(T, S, pressure, dSV_dT, dSV_dS, &
                             start, npts, Rho_T0_S0, dRho_dT, dRho_dS)
  real,    intent(in),  dimension(:) :: T         !< Potential temperature relative to the surface
                                                  !! [degC].
  real,    intent(in),  dimension(:) :: S         !< Salinity [PSU].
  real,    intent(in),  dimension(:) :: pressure  !< pressure [Pa].
  real,    intent(out), dimension(:) :: dSV_dS    !< The partial derivative of specific volume with
                                                  !! salinity [m3 kg-1 PSU-1].
  real,    intent(out), dimension(:) :: dSV_dT    !< The partial derivative of specific volume with
                                                  !! potential temperature [m3 kg-1 degC-1].
  integer, intent(in)                :: start     !< The starting point in the arrays.
  integer, intent(in)                :: npts      !< The number of values to calculate.
  real,    intent(in)                :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3].
  real,    intent(in)                :: dRho_dT   !< The derivative of density with
                                                  !! temperature, [kg m-3 degC-1].
  real,    intent(in)                :: dRho_dS   !< The derivative of density with
                                                  !! salinity [kg m-3 ppt-1].
  ! Local variables
  real :: I_rho2
  integer :: j

  do j=start,start+npts-1
    ! Sv = 1.0 / (Rho_T0_S0 + dRho_dT*T(j) + dRho_dS*S(j))
    I_rho2 = 1.0 / (Rho_T0_S0 + (dRho_dT*T(j) + dRho_dS*S(j)))**2
    dSV_dT(j) = -dRho_dT * I_rho2
    dSV_dS(j) = -dRho_dS * I_rho2
  enddo

end subroutine calculate_specvol_derivs_linear

!> This subroutine computes the in situ density of sea water (rho)
!! and the compressibility (drho/dp == C_sound^-2) at the given
!! salinity, potential temperature, and pressure.
subroutine calculate_compress_linear(T, S, pressure, rho, drho_dp, start, npts,&
                                     Rho_T0_S0, dRho_dT, dRho_dS)
  real,    intent(in),  dimension(:) :: T         !< Potential temperature relative to the surface
                                                  !! [degC].
  real,    intent(in),  dimension(:) :: S         !< Salinity [PSU].
  real,    intent(in),  dimension(:) :: pressure  !< pressure [Pa].
  real,    intent(out), dimension(:) :: rho       !< In situ density [kg m-3].
  real,    intent(out), dimension(:) :: drho_dp   !< The partial derivative of density with pressure
                                                  !! (also the inverse of the square of sound speed)
                                                  !! [s2 m-2].
  integer, intent(in)                :: start     !< The starting point in the arrays.
  integer, intent(in)                :: npts      !< The number of values to calculate.
  real,    intent(in)                :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3].
  real,    intent(in)                :: dRho_dT   !< The derivative of density with
                                                  !! temperature [kg m-3 degC-1].
  real,    intent(in)                :: dRho_dS   !< The derivative of density with
                                                  !! salinity [kg m-3 ppt-1].
  !  Local variables
  integer :: j

  do j=start,start+npts-1
    rho(j) = Rho_T0_S0 + dRho_dT*T(j) + dRho_dS*S(j)
    drho_dp(j) = 0.0
  enddo
end subroutine calculate_compress_linear

!>   This subroutine calculates analytical and nearly-analytical integrals of
!! pressure anomalies across layers, which are required for calculating the
!! finite-volume form pressure accelerations in a Boussinesq model.
subroutine int_density_dz_linear(T, S, z_t, z_b, rho_ref, rho_0_pres, G_e, HI, &
                 Rho_T0_S0, dRho_dT, dRho_dS, dpa, intz_dpa, intx_dpa, inty_dpa, &
                 bathyT, dz_neglect, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI        !< The horizontal index type for the arrays.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: T         !< Potential temperature relative to the surface
                                                 !! [degC].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: S         !< Salinity [PSU].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_t       !< Height at the top of the layer in depth units [Z ~> m].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(in)  :: z_b       !< Height at the top of the layer [Z ~> m].
  real,                 intent(in)  :: rho_ref   !< A mean density [R ~> kg m-3] or [kg m-3], that
                                                 !! is subtracted out to reduce the magnitude of
                                                 !! each of the integrals.
  real,                 intent(in)  :: rho_0_pres !< A density [R ~> kg m-3], used to calculate
                                                 !! the pressure (as p~=-z*rho_0_pres*G_e) used in
                                                 !! the equation of state. rho_0_pres is not used.
  real,                 intent(in)  :: G_e       !< The Earth's gravitational acceleration
                                                 !! [L2 Z-1 T-2 ~> m s-2] or [m2 Z-1 s-2 ~> m s-2].
  real,                 intent(in)  :: Rho_T0_S0 !< The density at T=0, S=0 [R ~> kg m-3] or [kg m-3].
  real,                 intent(in)  :: dRho_dT   !< The derivative of density with temperature,
                                                 !! [R degC-1 ~> kg m-3 degC-1] or [kg m-3 degC-1].
  real,                 intent(in)  :: dRho_dS   !< The derivative of density with salinity,
                                                 !! in [R ppt-1 ~> kg m-3 ppt-1] or [kg m-3 ppt-1].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(out) :: dpa       !< The change in the pressure anomaly across the
                                                 !! layer [R L2 T-2 ~> Pa] or [Pa].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(out) :: intz_dpa  !< The integral through the thickness of the layer
                                                 !! of the pressure anomaly relative to the anomaly
                                                 !! at the top of the layer [R L2 Z T-2 ~> Pa Z] or [Pa Z].
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed),  &
              optional, intent(out) :: intx_dpa  !< The integral in x of the difference between the
                                                 !! pressure anomaly at the top and bottom of the
                                                 !! layer divided by the x grid spacing [R L2 T-2 ~> Pa] or [Pa].
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB),  &
              optional, intent(out) :: inty_dpa  !< The integral in y of the difference between the
                                                 !! pressure anomaly at the top and bottom of the
                                                 !! layer divided by the y grid spacing [R L2 T-2 ~> Pa] or [Pa].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyT    !< The depth of the bathymetry [Z ~> m].
  real,       optional, intent(in)  :: dz_neglect !< A miniscule thickness change [Z ~> m].
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting to
                                                 !! interpolate T/S for top and bottom integrals.

  ! Local variables
  real :: rho_anom      ! The density anomaly from rho_ref [R ~> kg m-3].
  real :: raL, raR      ! rho_anom to the left and right [R ~> kg m-3].
  real :: dz, dzL, dzR  ! Layer thicknesses [Z ~> m].
  real :: hWght      ! A pressure-thickness below topography [Z ~> m].
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [Z ~> m].
  real :: iDenom     ! The inverse of the denominator in the weights [Z-2 ~> m-2].
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim].
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim].
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim].
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim].
  real :: intz(5)    ! The integrals of density with height at the
                     ! 5 sub-column locations [R L2 T-2 ~> Pa] or [Pa].
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  real, parameter :: C1_6 = 1.0/6.0, C1_90 = 1.0/90.0  ! Rational constants.
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, i, j, m

  ! These array bounds work for the indexing convention of the input arrays, but
  ! on the computational domain defined for the output arrays.
  Isq = HI%IscB ; Ieq = HI%IecB
  Jsq = HI%JscB ; Jeq = HI%JecB
  is = HI%isc ; ie = HI%iec
  js = HI%jsc ; je = HI%jec

  do_massWeight = .false.
  if (present(useMassWghtInterp)) then ; if (useMassWghtInterp) then
    do_massWeight = .true.
  ! if (.not.present(bathyT)) call MOM_error(FATAL, "int_density_dz_generic: "//&
  !     "bathyT must be present if useMassWghtInterp is present and true.")
  ! if (.not.present(dz_neglect)) call MOM_error(FATAL, "int_density_dz_generic: "//&
  !     "dz_neglect must be present if useMassWghtInterp is present and true.")
  endif ; endif

  do j=Jsq,Jeq+1 ; do i=Isq,Ieq+1
    dz = z_t(i,j) - z_b(i,j)
    rho_anom = (Rho_T0_S0 - rho_ref) + dRho_dT*T(i,j) + dRho_dS*S(i,j)
    dpa(i,j) = G_e*rho_anom*dz
    if (present(intz_dpa)) intz_dpa(i,j) = 0.5*G_e*rho_anom*dz**2
  enddo ; enddo

  if (present(intx_dpa)) then ; do j=js,je ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., -bathyT(i,j)-z_t(i+1,j), -bathyT(i+1,j)-z_t(i,j))

    if (hWght <= 0.0) then
      dzL = z_t(i,j) - z_b(i,j) ; dzR = z_t(i+1,j) - z_b(i+1,j)
      raL = (Rho_T0_S0 - rho_ref) + (dRho_dT*T(i,j) + dRho_dS*S(i,j))
      raR = (Rho_T0_S0 - rho_ref) + (dRho_dT*T(i+1,j) + dRho_dS*S(i+1,j))

      intx_dpa(i,j) = G_e*C1_6 * (dzL*(2.0*raL + raR) + dzR*(2.0*raR + raL))
    else
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i+1,j) - z_b(i+1,j)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom

      intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)
      do m=2,4
        wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
        wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

        dz = wt_L*(z_t(i,j) - z_b(i,j)) + wt_R*(z_t(i+1,j) - z_b(i+1,j))
        rho_anom = (Rho_T0_S0 - rho_ref) + &
                   (dRho_dT * (wtT_L*T(i,j) + wtT_R*T(i+1,j)) + &
                    dRho_dS * (wtT_L*S(i,j) + wtT_R*S(i+1,j)))
        intz(m) = G_e*rho_anom*dz
      enddo
      ! Use Bode's rule to integrate the values.
      intx_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    endif
  enddo ; enddo ; endif

  if (present(inty_dpa)) then ; do J=Jsq,Jeq ; do i=is,ie
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., -bathyT(i,j)-z_t(i,j+1), -bathyT(i,j+1)-z_t(i,j))

    if (hWght <= 0.0) then
      dzL = z_t(i,j) - z_b(i,j) ; dzR = z_t(i,j+1) - z_b(i,j+1)
      raL = (Rho_T0_S0 - rho_ref) + (dRho_dT*T(i,j) + dRho_dS*S(i,j))
      raR = (Rho_T0_S0 - rho_ref) + (dRho_dT*T(i,j+1) + dRho_dS*S(i,j+1))

      inty_dpa(i,j) = G_e*C1_6 * (dzL*(2.0*raL + raR) + dzR*(2.0*raR + raL))
    else
      hL = (z_t(i,j) - z_b(i,j)) + dz_neglect
      hR = (z_t(i,j+1) - z_b(i,j+1)) + dz_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom

      intz(1) = dpa(i,j) ; intz(5) = dpa(i+1,j)
      do m=2,4
        wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
        wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

        dz = wt_L*(z_t(i,j) - z_b(i,j)) + wt_R*(z_t(i,j+1) - z_b(i,j+1))
        rho_anom = (Rho_T0_S0 - rho_ref) + &
                   (dRho_dT * (wtT_L*T(i,j) + wtT_R*T(i,j+1)) + &
                    dRho_dS * (wtT_L*S(i,j) + wtT_R*S(i,j+1)))
        intz(m) = G_e*rho_anom*dz
      enddo
      ! Use Bode's rule to integrate the values.
      inty_dpa(i,j) = C1_90*(7.0*(intz(1)+intz(5)) + 32.0*(intz(2)+intz(4)) + &
                             12.0*intz(3))
    endif

  enddo ; enddo ; endif
end subroutine int_density_dz_linear

!> Calculates analytical and nearly-analytical integrals in
!! pressure across layers of geopotential anomalies, which are required for
!! calculating the finite-volume form pressure accelerations in a non-Boussinesq
!! model.  Specific volume is assumed to vary linearly between adjacent points.
subroutine int_spec_vol_dp_linear(T, S, p_t, p_b, alpha_ref, HI, Rho_T0_S0, &
               dRho_dT, dRho_dS, dza, intp_dza, intx_dza, inty_dza, halo_size, &
               bathyP, dP_neglect, useMassWghtInterp)
  type(hor_index_type), intent(in)  :: HI        !< The ocean's horizontal index type.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed),  &
                        intent(in)  :: T         !< Potential temperature relative to the surface
                                                 !! [degC].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed),  &
                        intent(in)  :: S         !< Salinity [PSU].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed),  &
                        intent(in)  :: p_t       !< Pressure at the top of the layer [R L2 T-2 ~> Pa] or [Pa].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed),  &
                        intent(in)  :: p_b       !< Pressure at the top of the layer [R L2 T-2 ~> Pa] or [Pa].
  real,                 intent(in)  :: alpha_ref   !< A mean specific volume that is subtracted out
                            !! to reduce the magnitude of each of the integrals [R-1 ~> m3 kg-1].
                            !! The calculation is mathematically identical with different values of
                            !! alpha_ref, but this reduces the effects of roundoff.
  real,                 intent(in)  :: Rho_T0_S0 !< The density at T=0, S=0 [R ~> kg m-3] or [kg m-3].
  real,                 intent(in)  :: dRho_dT   !< The derivative of density with temperature
                                                 !! [R degC-1 ~> kg m-3 degC-1] or [kg m-3 degC-1].
  real,                 intent(in)  :: dRho_dS   !< The derivative of density with salinity,
                                                 !! in [R ppt-1 ~> kg m-3 ppt-1] or [kg m-3 ppt-1].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                        intent(out) :: dza       !< The change in the geopotential anomaly across
                                                 !! the layer [L2 T-2 ~> m2 s-2] or [m2 s-2].
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(out) :: intp_dza  !< The integral in pressure through the layer of the
                                                 !! geopotential anomaly relative to the anomaly at the
                                                 !! bottom of the layer [R L4 T-4 ~> Pa m2 s-2] or [Pa m2 s-2].
  real, dimension(HI%IsdB:HI%IedB,HI%jsd:HI%jed), &
              optional, intent(out) :: intx_dza  !< The integral in x of the difference between the
                                                 !! geopotential anomaly at the top and bottom of
                                                 !! the layer divided by the x grid spacing
                                                 !! [L2 T-2 ~> m2 s-2] or [m2 s-2].
  real, dimension(HI%isd:HI%ied,HI%JsdB:HI%JedB), &
              optional, intent(out) :: inty_dza  !< The integral in y of the difference between the
                                                 !! geopotential anomaly at the top and bottom of
                                                 !! the layer divided by the y grid spacing
                                                 !! [L2 T-2 ~> m2 s-2] or [m2 s-2].
  integer,    optional, intent(in)  :: halo_size !< The width of halo points on which to calculate dza.
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
              optional, intent(in)  :: bathyP    !< The pressure at the bathymetry [R L2 T-2 ~> Pa] or [Pa]
  real,       optional, intent(in)  :: dP_neglect !< A miniscule pressure change with
                                                 !! the same units as p_t [R L2 T-2 ~> Pa] or [Pa]
  logical,    optional, intent(in)  :: useMassWghtInterp !< If true, uses mass weighting
                            !! to interpolate T/S for top and bottom integrals.
  ! Local variables
  real :: dRho_TS       ! The density anomaly due to T and S [R ~> kg m-3] or [kg m-3].
  real :: alpha_anom    ! The specific volume anomaly from 1/rho_ref [R-1 ~> m3 kg-1] or [m3 kg-1].
  real :: aaL, aaR      ! The specific volume anomaly to the left and right [R-1 ~> m3 kg-1] or [m3 kg-1].
  real :: dp, dpL, dpR  ! Layer pressure thicknesses [R L2 T-2 ~> Pa] or [Pa].
  real :: hWght      ! A pressure-thickness below topography [R L2 T-2 ~> Pa] or [Pa].
  real :: hL, hR     ! Pressure-thicknesses of the columns to the left and right [R L2 T-2 ~> Pa] or [Pa].
  real :: iDenom     ! The inverse of the denominator in the weights [T4 R-2 L-2 ~> Pa-2] or [Pa-2].
  real :: hWt_LL, hWt_LR ! hWt_LA is the weighted influence of A on the left column [nondim].
  real :: hWt_RL, hWt_RR ! hWt_RA is the weighted influence of A on the right column [nondim].
  real :: wt_L, wt_R ! The linear weights of the left and right columns [nondim].
  real :: wtT_L, wtT_R ! The weights for tracers from the left and right columns [nondim].
  real :: intp(5)    ! The integrals of specific volume with pressure at the
                     ! 5 sub-column locations [L2 T-2 ~> m2 s-2] or [m2 s-2].
  logical :: do_massWeight ! Indicates whether to do mass weighting.
  real, parameter :: C1_6 = 1.0/6.0, C1_90 = 1.0/90.0  ! Rational constants.
  integer :: Isq, Ieq, Jsq, Jeq, ish, ieh, jsh, jeh, i, j, m, halo

  Isq = HI%IscB ; Ieq = HI%IecB ; Jsq = HI%JscB ; Jeq = HI%JecB
  halo = 0 ; if (present(halo_size)) halo = MAX(halo_size,0)
  ish = HI%isc-halo ; ieh = HI%iec+halo ; jsh = HI%jsc-halo ; jeh = HI%jec+halo
  if (present(intx_dza)) then ; ish = MIN(Isq,ish) ; ieh = MAX(Ieq+1,ieh); endif
  if (present(inty_dza)) then ; jsh = MIN(Jsq,jsh) ; jeh = MAX(Jeq+1,jeh); endif

  do_massWeight = .false.
  if (present(useMassWghtInterp)) then ; if (useMassWghtInterp) then
    do_massWeight = .true.
!    if (.not.present(bathyP)) call MOM_error(FATAL, "int_spec_vol_dp_generic: "//&
!        "bathyP must be present if useMassWghtInterp is present and true.")
!    if (.not.present(dP_neglect)) call MOM_error(FATAL, "int_spec_vol_dp_generic: "//&
!        "dP_neglect must be present if useMassWghtInterp is present and true.")
  endif ; endif

  do j=jsh,jeh ; do i=ish,ieh
    dp = p_b(i,j) - p_t(i,j)
    dRho_TS = dRho_dT*T(i,j) + dRho_dS*S(i,j)
    ! alpha_anom = 1.0/(Rho_T0_S0  + dRho_TS)) - alpha_ref
    alpha_anom = ((1.0-Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)
    dza(i,j) = alpha_anom*dp
    if (present(intp_dza)) intp_dza(i,j) = 0.5*alpha_anom*dp**2
  enddo ; enddo

  if (present(intx_dza)) then ; do j=HI%jsc,HI%jec ; do I=Isq,Ieq
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., bathyP(i,j)-p_t(i+1,j), bathyP(i+1,j)-p_t(i,j))

    if (hWght <= 0.0) then
      dpL = p_b(i,j) - p_t(i,j) ; dpR = p_b(i+1,j) - p_t(i+1,j)
      dRho_TS = dRho_dT*T(i,j) + dRho_dS*S(i,j)
      aaL = ((1.0 - Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)
      dRho_TS = dRho_dT*T(i+1,j) + dRho_dS*S(i+1,j)
      aaR = ((1.0 - Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)

      intx_dza(i,j) = C1_6 * (2.0*(dpL*aaL + dpR*aaR) + (dpL*aaR + dpR*aaL))
    else
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i+1,j) - p_t(i+1,j)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom

      intp(1) = dza(i,j) ; intp(5) = dza(i+1,j)
      do m=2,4
        wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
        wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

        ! T, S, and p are interpolated in the horizontal.  The p interpolation
        ! is linear, but for T and S it may be thickness wekghted.
        dp = wt_L*(p_b(i,j) - p_t(i,j)) + wt_R*(p_b(i+1,j) - p_t(i+1,j))

        dRho_TS = dRho_dT*(wtT_L*T(i,j) + wtT_R*T(i+1,j)) + &
                  dRho_dS*(wtT_L*S(i,j) + wtT_R*S(i+1,j))
        ! alpha_anom = 1.0/(Rho_T0_S0  + dRho_TS)) - alpha_ref
        alpha_anom = ((1.0-Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)
        intp(m) = alpha_anom*dp
      enddo
      ! Use Bode's rule to integrate the interface height anomaly values in y.
      intx_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                             12.0*intp(3))
    endif
  enddo ; enddo ; endif

  if (present(inty_dza)) then ; do J=Jsq,Jeq ; do i=HI%isc,HI%iec
    ! hWght is the distance measure by which the cell is violation of
    ! hydrostatic consistency. For large hWght we bias the interpolation of
    ! T & S along the top and bottom integrals, akin to thickness weighting.
    hWght = 0.0
    if (do_massWeight) &
      hWght = max(0., bathyP(i,j)-p_t(i,j+1), bathyP(i,j+1)-p_t(i,j))

    if (hWght <= 0.0) then
      dpL = p_b(i,j) - p_t(i,j) ; dpR = p_b(i,j+1) - p_t(i,j+1)
      dRho_TS = dRho_dT*T(i,j) + dRho_dS*S(i,j)
      aaL = ((1.0 - Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)
      dRho_TS = dRho_dT*T(i,j+1) + dRho_dS*S(i,j+1)
      aaR = ((1.0 - Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)

      inty_dza(i,j) = C1_6 * (2.0*(dpL*aaL + dpR*aaR) + (dpL*aaR + dpR*aaL))
    else
      hL = (p_b(i,j) - p_t(i,j)) + dP_neglect
      hR = (p_b(i,j+1) - p_t(i,j+1)) + dP_neglect
      hWght = hWght * ( (hL-hR)/(hL+hR) )**2
      iDenom = 1.0 / ( hWght*(hR + hL) + hL*hR )
      hWt_LL = (hWght*hL + hR*hL) * iDenom ; hWt_LR = (hWght*hR) * iDenom
      hWt_RR = (hWght*hR + hR*hL) * iDenom ; hWt_RL = (hWght*hL) * iDenom

      intp(1) = dza(i,j) ; intp(5) = dza(i,j+1)
      do m=2,4
        wt_L = 0.25*real(5-m) ; wt_R = 1.0-wt_L
        wtT_L = wt_L*hWt_LL + wt_R*hWt_RL ; wtT_R = wt_L*hWt_LR + wt_R*hWt_RR

        ! T, S, and p are interpolated in the horizontal.  The p interpolation
        ! is linear, but for T and S it may be thickness wekghted.
        dp = wt_L*(p_b(i,j) - p_t(i,j)) + wt_R*(p_b(i,j+1) - p_t(i,j+1))

        dRho_TS = dRho_dT*(wtT_L*T(i,j) + wtT_R*T(i,j+1)) + &
                  dRho_dS*(wtT_L*S(i,j) + wtT_R*S(i,j+1))
        ! alpha_anom = 1.0/(Rho_T0_S0  + dRho_TS)) - alpha_ref
        alpha_anom = ((1.0-Rho_T0_S0*alpha_ref) - dRho_TS*alpha_ref) / (Rho_T0_S0 + dRho_TS)
        intp(m) = alpha_anom*dp
      enddo
      ! Use Bode's rule to integrate the interface height anomaly values in y.
      inty_dza(i,j) = C1_90*(7.0*(intp(1)+intp(5)) + 32.0*(intp(2)+intp(4)) + &
                             12.0*intp(3))
    endif
  enddo ; enddo ; endif
end subroutine int_spec_vol_dp_linear

end module MOM_EOS_linear
