!> A generic type for equations of state
module MOM_EOS_base_type

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public EOS_base

!> The base class for implementations of the equation of state
type, abstract :: EOS_base

contains

  ! The following functions/subroutines are deferred and must be provided specifically by each EOS

  !> Deferred implementation of the in-situ density as an elemental function [kg m-3]
  procedure(i_density_elem), deferred :: density_elem
  !> Deferred implementation of the in-situ density anomaly as an elemental function [kg m-3]
  procedure(i_density_anomaly_elem), deferred :: density_anomaly_elem
  !> Deferred implementation of the in-situ specific volume as an elemental function [m3 kg-1]
  procedure(i_spec_vol_elem), deferred :: spec_vol_elem
  !> Deferred implementation of the in-situ specific volume anomaly as an elemental function [m3 kg-1]
  procedure(i_spec_vol_anomaly_elem), deferred :: spec_vol_anomaly_elem
  !> Deferred implementation of the calculation of derivatives of density
  procedure(i_calculate_density_derivs_elem), deferred :: calculate_density_derivs_elem
  !> Deferred implementation of the calculation of second derivatives of density
  procedure(i_calculate_density_second_derivs_elem), deferred :: calculate_density_second_derivs_elem
  !> Deferred implementation of the calculation of derivatives of specific volume
  procedure(i_calculate_specvol_derivs_elem), deferred :: calculate_specvol_derivs_elem
  !> Deferred implementation of the calculation of compressibility
  procedure(i_calculate_compress_elem), deferred :: calculate_compress_elem
  !> Deferred implementation of the range query function
  procedure(i_EOS_fit_range), deferred :: EOS_fit_range

  ! The following functions/subroutines are shared across all EOS and provided by this module
  !> Returns the in-situ density or density anomaly [kg m-3]
  procedure :: density_fn => a_density_fn
  !> Returns the in-situ specific volume or specific volume anomaly [m3 kg-1]
  procedure :: spec_vol_fn => a_spec_vol_fn
  !> Calculates the in-situ density or density anomaly for scalar inputs [m3 kg-1]
  procedure :: calculate_density_scalar => a_calculate_density_scalar
  !> Calculates the in-situ density or density anomaly for array inputs [m3 kg-1]
  procedure :: calculate_density_array => a_calculate_density_array
  !> Calculates the in-situ specific volume or specific volume anomaly for scalar inputs [m3 kg-1]
  procedure :: calculate_spec_vol_scalar => a_calculate_spec_vol_scalar
  !> Calculates the in-situ specific volume or specific volume anomaly for array inputs [m3 kg-1]
  procedure :: calculate_spec_vol_array => a_calculate_spec_vol_array
  !> Calculates the derivatives of density for scalar inputs
  procedure :: calculate_density_derivs_scalar => a_calculate_density_derivs_scalar
  !> Calculates the derivatives of density for array inputs
  procedure :: calculate_density_derivs_array => a_calculate_density_derivs_array
  !> Calculates the second derivatives of density for scalar inputs
  procedure :: calculate_density_second_derivs_scalar => a_calculate_density_second_derivs_scalar
  !> Calculates the second derivatives of density for array inputs
  procedure :: calculate_density_second_derivs_array => a_calculate_density_second_derivs_array
  !> Calculates the derivatives of specific volume for array inputs
  procedure :: calculate_specvol_derivs_array => a_calculate_specvol_derivs_array
  !> Calculates the compressibility for array inputs
  procedure :: calculate_compress_array => a_calculate_compress_array

end type EOS_base

interface

  !> In situ density [kg m-3]
  !!
  !! This is an elemental function that can be applied to any combination of scalar and array inputs.
  real elemental function i_density_elem(this, T, S, pressure)
    import :: EOS_base
    class(EOS_base), intent(in) :: this     !< This EOS
    real,            intent(in) :: T        !< Potential temperature relative to the surface [degC]
    real,            intent(in) :: S        !< Salinity [PSU]
    real,            intent(in) :: pressure !< Pressure [Pa]

  end function i_density_elem

  !> In situ density anomaly [kg m-3]
  !!
  !! This is an elemental function that can be applied to any combination of scalar and array inputs.
  real elemental function i_density_anomaly_elem(this, T, S, pressure, rho_ref)
    import :: EOS_base
    class(EOS_base), intent(in) :: this     !< This EOS
    real,            intent(in) :: T        !< Potential temperature relative to the surface [degC]
    real,            intent(in) :: S        !< Salinity [PSU]
    real,            intent(in) :: pressure !< Pressure [Pa]
    real,            intent(in) :: rho_ref  !< A reference density [kg m-3]

  end function i_density_anomaly_elem

  !> In situ specific volume [m3 kg-1]
  !!
  !! This is an elemental function that can be applied to any combination of scalar and array inputs.
  real elemental function i_spec_vol_elem(this, T, S, pressure)
    import :: EOS_base
    class(EOS_base), intent(in) :: this     !< This EOS
    real,            intent(in) :: T        !< Potential temperature relative to the surface [degC]
    real,            intent(in) :: S        !< Salinity [PSU]
    real,            intent(in) :: pressure !< Pressure [Pa]

  end function i_spec_vol_elem

  !> In situ specific volume anomaly [m3 kg-1]
  !!
  !! This is an elemental function that can be applied to any combination of scalar and array inputs.
  real elemental function i_spec_vol_anomaly_elem(this, T, S, pressure, spv_ref)
    import :: EOS_base
    class(EOS_base), intent(in) :: this     !< This EOS
    real,            intent(in) :: T        !< Potential temperature relative to the surface [degC]
    real,            intent(in) :: S        !< Salinity [PSU]
    real,            intent(in) :: pressure !< Pressure [Pa]
    real,            intent(in) :: spv_ref  !< A reference specific volume [m3 kg-1]

  end function i_spec_vol_anomaly_elem

  !> Calculate the partial derivatives of density with potential temperature and salinity
  elemental subroutine i_calculate_density_derivs_elem(this, T, S, pressure, drho_dT, drho_dS)
    import :: EOS_base
    class(EOS_base), intent(in) :: this      !< This EOS
    real,            intent(in)  :: T        !< Potential temperature relative to the surface [degC]
    real,            intent(in)  :: S        !< Salinity [PSU]
    real,            intent(in)  :: pressure !< Pressure [Pa]
    real,            intent(out) :: drho_dT  !< The partial derivative of density with potential
                                             !! temperature [kg m-3 degC-1]
    real,            intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                             !! in [kg m-3 PSU-1]

  end subroutine i_calculate_density_derivs_elem

  !> Calculate the partial derivatives of specific volume with temperature and salinity
  elemental subroutine i_calculate_specvol_derivs_elem(this, T, S, pressure, dSV_dT, dSV_dS)
    import :: EOS_base
    class(EOS_base), intent(in)    :: this     !< This EOS
    real,            intent(in)    :: T        !< Potential temperature [degC]
    real,            intent(in)    :: S        !< Salinity [PSU]
    real,            intent(in)    :: pressure !< Pressure [Pa]
    real,            intent(inout) :: dSV_dT   !< The partial derivative of specific volume with
                                               !! potential temperature [m3 kg-1 degC-1]
    real,            intent(inout) :: dSV_dS   !< The partial derivative of specific volume with
                                               !! salinity [m3 kg-1 PSU-1]

  end subroutine i_calculate_specvol_derivs_elem

  !> Calculate second derivatives of density with respect to temperature, salinity, and pressure
  elemental subroutine i_calculate_density_second_derivs_elem(this, T, S, pressure, &
                          drho_ds_ds, drho_ds_dt, drho_dt_dt, drho_ds_dp, drho_dt_dp)
    import :: EOS_base
    class(EOS_base), intent(in)    :: this     !< This EOS
    real,            intent(in)    :: T !< Potential temperature referenced to 0 dbar [degC]
    real,            intent(in)    :: S !< Salinity [PSU]
    real,            intent(in)    :: pressure !< Pressure [Pa]
    real,            intent(inout) :: drho_ds_ds !< Partial derivative of beta with respect
                                                 !! to S [kg m-3 PSU-2]
    real,            intent(inout) :: drho_ds_dt !< Partial derivative of beta with respect
                                                 !! to T [kg m-3 PSU-1 degC-1]
    real,            intent(inout) :: drho_dt_dt !< Partial derivative of alpha with respect
                                                 !! to T [kg m-3 degC-2]
    real,            intent(inout) :: drho_ds_dp !< Partial derivative of beta with respect
                                                 !! to pressure [kg m-3 PSU-1 Pa-1] = [s2 m-2 PSU-1]
    real,            intent(inout) :: drho_dt_dp !< Partial derivative of alpha with respect
                                                 !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]

  end subroutine i_calculate_density_second_derivs_elem

  !> Compute the in situ density of sea water (rho) and the compressibility (drho/dp == C_sound^-2)
  !! at the given salinity, potential temperature and pressure
  elemental subroutine i_calculate_compress_elem(this, T, S, pressure, rho, drho_dp)
    import :: EOS_base
    class(EOS_base), intent(in)  :: this     !< This EOS
    real,            intent(in)  :: T        !< Potential temperature relative to the surface [degC]
    real,            intent(in)  :: S        !< Salinity [PSU]
    real,            intent(in)  :: pressure !< Pressure [Pa]
    real,            intent(out) :: rho      !< In situ density [kg m-3]
    real,            intent(out) :: drho_dp  !< The partial derivative of density with pressure (or
                                             !! the inverse of the square of sound speed) [s2 m-2]

  end subroutine i_calculate_compress_elem

  !> Return the range of temperatures, salinities and pressures for which the equations of state has been
  !! fitted or is valid. Care should be taken when applying this equation of state outside of its fit range.
  subroutine i_EOS_fit_range(this, T_min, T_max, S_min, S_max, p_min, p_max)
    import :: EOS_base
    class(EOS_base), intent(in) :: this     !< This EOS
    real, optional, intent(out) :: T_min !< The minimum potential temperature over which this EoS is fitted [degC]
    real, optional, intent(out) :: T_max !< The maximum potential temperature over which this EoS is fitted [degC]
    real, optional, intent(out) :: S_min !< The minimum practical salinity over which this EoS is fitted [PSU]
    real, optional, intent(out) :: S_max !< The maximum practical salinity over which this EoS is fitted [PSU]
    real, optional, intent(out) :: p_min !< The minimum pressure over which this EoS is fitted [Pa]
    real, optional, intent(out) :: p_max !< The maximum pressure over which this EoS is fitted [Pa]

  end subroutine i_EOS_fit_range

end interface

contains

  !> In situ density [kg m-3]
  real function a_density_fn(this, T, S, pressure, rho_ref)
    class(EOS_base), intent(in) :: this     !< This EOS
    real,           intent(in) :: T        !< Potential temperature relative to the surface [degC]
    real,           intent(in) :: S        !< Salinity [PSU]
    real,           intent(in) :: pressure !< Pressure [Pa]
    real, optional, intent(in) :: rho_ref  !< A reference density [kg m-3]

    if (present(rho_ref)) then
      a_density_fn = this%density_anomaly_elem(T, S, pressure, rho_ref)
    else
      a_density_fn = this%density_elem(T, S, pressure)
    endif

  end function a_density_fn

  !> Calculate the in-situ density for scalar inputs and outputs.
  subroutine a_calculate_density_scalar(this, T, S, pressure, rho, rho_ref)
    class(EOS_base), intent(in) :: this     !< This EOS
    real,           intent(in)  :: T        !< Potential temperature relative to the surface [degC]
    real,           intent(in)  :: S        !< Salinity [PSU]
    real,           intent(in)  :: pressure !< Pressure [Pa]
    real,           intent(out) :: rho      !< In situ density [kg m-3]
    real, optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]

    if (present(rho_ref)) then
      rho = this%density_anomaly_elem(T, S, pressure, rho_ref)
    else
      rho = this%density_elem(T, S, pressure)
    endif

  end subroutine a_calculate_density_scalar

  !> Calculate the in-situ density for 1D arraya inputs and outputs.
  subroutine a_calculate_density_array(this, T, S, pressure, rho, start, npts, rho_ref)
    class(EOS_base), intent(in) :: this     !< This EOS
    real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
    real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
    real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
    real, dimension(:), intent(out) :: rho      !< In situ density [kg m-3]
    integer,            intent(in)  :: start    !< The starting index for calculations
    integer,            intent(in)  :: npts     !< The number of values to calculate
    real,     optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]

    ! Local variables
    integer :: js, je

    js = start
    je = start+npts-1

    if (present(rho_ref)) then
      rho(js:je) = this%density_anomaly_elem(T(js:je), S(js:je), pressure(js:je), rho_ref)
    else
      rho(js:je) = this%density_elem(T(js:je), S(js:je), pressure(js:je))
    endif

  end subroutine a_calculate_density_array

  !> In situ specific volume [m3 kg-1]
  real function a_spec_vol_fn(this, T, S, pressure, spv_ref)
    class(EOS_base), intent(in) :: this     !< This EOS
    real,           intent(in) :: T        !< Potential temperature relative to the surface [degC]
    real,           intent(in) :: S        !< Salinity [PSU]
    real,           intent(in) :: pressure !< Pressure [Pa]
    real, optional, intent(in) :: spv_ref  !< A reference specific volume [m3 kg-1]

    if (present(spv_ref)) then
      a_spec_vol_fn = this%spec_vol_anomaly_elem(T, S, pressure, spv_ref)
    else
      a_spec_vol_fn = this%spec_vol_elem(T, S, pressure)
    endif

  end function a_spec_vol_fn

  !> Calculate the in-situ specific volume for scalar inputs and outputs.
  subroutine a_calculate_spec_vol_scalar(this, T, S, pressure, specvol, spv_ref)
    class(EOS_base), intent(in) :: this     !< This EOS
    real,           intent(in)  :: T        !< Potential temperature relative to the surface [degC]
    real,           intent(in)  :: S        !< Salinity [PSU]
    real,           intent(in)  :: pressure !< Pressure [Pa]
    real,           intent(out) :: specvol  !< In situ specific volume [m3 kg-1]
    real, optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1]

    if (present(spv_ref)) then
      specvol = this%spec_vol_anomaly_elem(T, S, pressure, spv_ref)
    else
      specvol = this%spec_vol_elem(T, S, pressure)
    endif

  end subroutine a_calculate_spec_vol_scalar

  !> Calculate the in-situ specific volume for 1D array inputs and outputs.
  subroutine a_calculate_spec_vol_array(this, T, S, pressure, specvol, start, npts, spv_ref)
    class(EOS_base), intent(in) :: this     !< This EOS
    real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
    real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
    real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
    real, dimension(:), intent(out) :: specvol  !< In situ specific volume [m3 kg-1]
    integer,            intent(in)  :: start    !< The starting index for calculations
    integer,            intent(in)  :: npts     !< The number of values to calculate
    real,     optional, intent(in)  :: spv_ref  !< A reference specific volume [m3 kg-1]

    ! Local variables
    integer :: js, je

    js = start
    je = start+npts-1

    if (present(spv_ref)) then
      specvol(js:je) = this%spec_vol_anomaly_elem(T(js:je), S(js:je), pressure(js:je), spv_ref)
    else
      specvol(js:je) = this%spec_vol_elem(T(js:je), S(js:je), pressure(js:je) )
    endif

  end subroutine a_calculate_spec_vol_array

  !> Calculate the derivatives of density with respect to temperature, salinity and pressure
  !! for scalar inputs
  subroutine a_calculate_density_derivs_scalar(this, T, S, P, drho_dT, drho_dS)
    class(EOS_base), intent(in) :: this !< This EOS
    real, intent(in)  :: T       !< Potential temperature referenced to 0 dbar
    real, intent(in)  :: S       !< Salinity [PSU]
    real, intent(in)  :: P       !< Pressure [Pa]
    real, intent(out) :: drho_dT !< The partial derivative of density with potential
                                 !! temperature [kg m-3 degC-1]
    real, intent(out) :: drho_dS !< The partial derivative of density with salinity,
                                 !! in [kg m-3 PSU-1]

    call this%calculate_density_derivs_elem(T, S, P, drho_dt, drho_ds)

  end subroutine a_calculate_density_derivs_scalar

  !> Calculate the derivatives of density with respect to temperature, salinity and pressure
  !! for array inputs
  subroutine a_calculate_density_derivs_array(this, T, S, pressure, drho_dT, drho_dS, start, npts)
    class(EOS_base),    intent(in)  :: this     !< This EOS
    real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
    real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
    real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
    real, dimension(:), intent(out) :: drho_dT  !< The partial derivative of density with potential
                                                !! temperature [kg m-3 degC-1]
    real, dimension(:), intent(out) :: drho_dS  !< The partial derivative of density with salinity,
                                                !! in [kg m-3 PSU-1]
    integer,            intent(in)  :: start    !< The starting index for calculations
    integer,            intent(in)  :: npts     !< The number of values to calculate

    ! Local variables
    integer :: js, je

    js = start
    je = start+npts-1

    call this%calculate_density_derivs_elem(T(js:je), S(js:je), pressure(js:je), drho_dt(js:je), drho_ds(js:je))

  end subroutine a_calculate_density_derivs_array

  !> Calculate the second derivatives of density with respect to temperature, salinity and pressure
  !! for scalar inputs
  subroutine a_calculate_density_second_derivs_scalar(this, T, S, pressure, &
                     drho_ds_ds, drho_ds_dt, drho_dt_dt, drho_ds_dp, drho_dt_dp)
    class(EOS_base), intent(in)  :: this       !< This EOS
    real,            intent(in)  :: T          !< Potential temperature referenced to 0 dbar
    real,            intent(in)  :: S          !< Salinity [PSU]
    real,            intent(in)  :: pressure   !< Pressure [Pa]
    real,            intent(out) :: drho_ds_ds !< Partial derivative of beta with respect
                                               !! to S [kg m-3 PSU-2]
    real,            intent(out) :: drho_ds_dt !< Partial derivative of beta with respect
                                               !! to T [kg m-3 PSU-1 degC-1]
    real,            intent(out) :: drho_dt_dt !< Partial derivative of alpha with respect
                                               !! to T [kg m-3 degC-2]
    real,            intent(out) :: drho_ds_dp !< Partial derivative of beta with respect
                                               !! to pressure [kg m-3 PSU-1 Pa-1] = [s2 m-2 PSU-1]
    real,            intent(out) :: drho_dt_dp !< Partial derivative of alpha with respect
                                               !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]

    call this%calculate_density_second_derivs_elem(T, S, pressure, &
                      drho_ds_ds, drho_ds_dt, drho_dt_dt, drho_ds_dp, drho_dt_dp)

  end subroutine a_calculate_density_second_derivs_scalar

  !> Calculate the second derivatives of density with respect to temperature, salinity and pressure
  !! for array inputs
  subroutine a_calculate_density_second_derivs_array(this, T, S, pressure, &
                     drho_ds_ds, drho_ds_dt, drho_dt_dt, drho_ds_dp, drho_dt_dp, start, npts)
    class(EOS_base),    intent(in)  :: this       !< This EOS
    real, dimension(:), intent(in)  :: T          !< Potential temperature referenced to 0 dbar
    real, dimension(:), intent(in)  :: S          !< Salinity [PSU]
    real, dimension(:), intent(in)  :: pressure   !< Pressure [Pa]
    real, dimension(:), intent(out) :: drho_ds_ds !< Partial derivative of beta with respect
                                                  !! to S [kg m-3 PSU-2]
    real, dimension(:), intent(out) :: drho_ds_dt !< Partial derivative of beta with respect
                                                  !! to T [kg m-3 PSU-1 degC-1]
    real, dimension(:), intent(out) :: drho_dt_dt !< Partial derivative of alpha with respect
                                                  !! to T [kg m-3 degC-2]
    real, dimension(:), intent(out) :: drho_ds_dp !< Partial derivative of beta with respect
                                                  !! to pressure [kg m-3 PSU-1 Pa-1] = [s2 m-2 PSU-1]
    real, dimension(:), intent(out) :: drho_dt_dp !< Partial derivative of alpha with respect
                                                  !! to pressure [kg m-3 degC-1 Pa-1] = [s2 m-2 degC-1]
    integer,            intent(in)  :: start      !< The starting index for calculations
    integer,            intent(in)  :: npts       !< The number of values to calculate

    ! Local variables
    integer :: js, je

    js = start
    je = start+npts-1

    call this%calculate_density_second_derivs_elem(T(js:je), S(js:je), pressure(js:je), &
                              drho_ds_ds(js:je), drho_ds_dt(js:je), drho_dt_dt(js:je), &
                              drho_ds_dp(js:je), drho_dt_dp(js:je))

  end subroutine a_calculate_density_second_derivs_array

  !> Calculate the partial derivatives of specific volume with temperature and salinity
  !! for array inputs
  subroutine a_calculate_specvol_derivs_array(this, T, S, pressure, dSV_dT, dSV_dS, start, npts)
    class(EOS_base),    intent(in)    :: this     !< This EOS
    real, dimension(:), intent(in)    :: T        !< Potential temperature [degC]
    real, dimension(:), intent(in)    :: S        !< Salinity [PSU]
    real, dimension(:), intent(in)    :: pressure !< Pressure [Pa]
    real, dimension(:), intent(inout) :: dSV_dT   !< The partial derivative of specific volume with
                                                  !! potential temperature [m3 kg-1 degC-1]
    real, dimension(:), intent(inout) :: dSV_dS   !< The partial derivative of specific volume with
                                                  !! salinity [m3 kg-1 PSU-1]
    integer,            intent(in)    :: start    !< The starting index for calculations
    integer,            intent(in)    :: npts     !< The number of values to calculate

    ! Local variables
    integer :: js, je

    js = start
    je = start+npts-1

    call this%calculate_specvol_derivs_elem(T(js:je), S(js:je), pressure(js:je), &
                                            dSV_dT(js:je), dSV_dS(js:je))

  end subroutine a_calculate_specvol_derivs_array

  !> Compute the in situ density of sea water (rho) and the compressibility (drho/dp == C_sound^-2)
  !! at the given salinity, potential temperature and pressure for array inputs
  subroutine a_calculate_compress_array(this, T, S, pressure, rho, drho_dp, start, npts)
    class(EOS_base),    intent(in)  :: this     !< This EOS
    real, dimension(:), intent(in)  :: T        !< Potential temperature relative to the surface [degC]
    real, dimension(:), intent(in)  :: S        !< Salinity [PSU]
    real, dimension(:), intent(in)  :: pressure !< Pressure [Pa]
    real, dimension(:), intent(out) :: rho      !< In situ density [kg m-3]
    real, dimension(:), intent(out) :: drho_dp  !< The partial derivative of density with pressure (or
                                                !! the inverse of the square of sound speed) [s2 m-2]
    integer,            intent(in)  :: start    !< The starting index for calculations
    integer,            intent(in)  :: npts     !< The number of values to calculate

    ! Local variables
    integer :: js, je

    js = start
    je = start+npts-1

    call this%calculate_compress_elem(T(js:je), S(js:je), pressure(js:je), &
                                      rho(js:je), drho_dp(js:je))

  end subroutine a_calculate_compress_array

!> \namespace mom_eos_base_type
!!
!! \section section_EOS_base_type Generic EOS type
!!

end module MOM_EOS_base_type
