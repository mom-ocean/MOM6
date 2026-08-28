! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module provides a common datatype to provide forcing for MARBL tracers
!! regardless of driver
module MARBL_forcing_mod

!! This module exists to house code used by multiple drivers in config_src/
!! for passing forcing fields to MARBL
!! (This comment can go in the wiki on the NCAR fork?)

use MOM_diag_mediator,        only : safe_alloc_ptr, diag_ctrl
use MOM_time_manager,         only : time_type
use MOM_error_handler,        only : MOM_error, WARNING, FATAL
use MOM_file_parser,          only : get_param, log_param, param_file_type
use MOM_grid,                 only : ocean_grid_type
use MOM_unit_scaling,         only : unit_scale_type
use MOM_interpolate,          only : external_field, init_external_field, time_interp_external
use MOM_io,                   only : slasher
use marbl_constants_mod,      only : molw_Fe
use MOM_forcing_type,         only : forcing

implicit none ; private

#include <MOM_memory.h>

public :: MARBL_forcing_init
public :: convert_driver_fields_to_forcings

!> Control structure for this module
type, public :: marbl_forcing_CS ; private
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.

  real    :: dust_ratio_thres               !< coarse/fine dust ratio threshold [1]
  real    :: fe_bioavail_frac_offset        !< offset for iron bioavailability fraction [1]
  real    :: atm_fe_to_bc_ratio             !< atmospheric iron to black carbon ratio [1]
  real    :: seaice_fe_to_bc_ratio          !< sea-ice iron to black carbon ratio [1]
  real    :: iron_frac_in_atm_fine_dust     !< Fraction of fine dust from the atmosphere that is iron [1]
  real    :: iron_frac_in_atm_coarse_dust   !< Fraction of coarse dust from the atmosphere that is iron [1]
  real    :: iron_frac_in_seaice_dust       !< Fraction of dust from the sea ice that is iron [1]
  real    :: atm_co2_const                  !< atmospheric CO2 (if specifying a constant value) [ppm]
  real    :: atm_alt_co2_const              !< alternate atmospheric CO2 for _ALT_CO2 tracers
                                            !! (if specifying a constant value) [ppm]

  logical :: use_MARBL_tracers    !< most functions can return immediately
                                  !! MARBL tracers are turned off
  integer :: atm_co2_iopt         !< Integer version of atm_co2_opt, which determines source of atm_co2
  integer :: atm_alt_co2_iopt     !< Integer version of atm_alt_co2_opt, which determines source of atm_alt_co2

end type marbl_forcing_CS

! Module parameters
integer, parameter :: atm_co2_constant_iopt = 0     !< module parameter denoting atm_co2_opt = 'constant'
integer, parameter :: atm_co2_prognostic_iopt = 1   !< module parameter denoting atm_co2_opt = 'diagnostic'
integer, parameter :: atm_co2_diagnostic_iopt = 2   !< module parameter denoting atm_co2_opt = 'prognostic'

contains

  subroutine MARBL_forcing_init(G, US, param_file, diag, day, inputdir, use_MARBL_tracers, CS)
    type(ocean_grid_type),           intent(in)    :: G           !< The ocean's grid structure
    type(unit_scale_type),           intent(in)    :: US          !< A dimensional unit scaling type
    type(param_file_type),           intent(in)    :: param_file  !< A structure to parse for run-time parameters
    type(diag_ctrl), target,         intent(in)    :: diag        !< Structure used to regulate diagnostic output.
    type(time_type), target,         intent(in)    :: day         !< Time of the start of the run.
    character(len=*),                intent(in)    :: inputdir    !< Directory containing input files
    logical,                         intent(in)    :: use_MARBL_tracers   !< Is MARBL tracer package active?
    type(marbl_forcing_CS), pointer, intent(inout) :: CS          !< A pointer that is set to point to control
                                                                  !! structure for MARBL forcing

    character(len=40)  :: mdl = "MARBL_forcing_mod"  ! This module's name.
    character(len=15)  :: atm_co2_opt
    character(len=200) :: err_message

    if (associated(CS)) then
      call MOM_error(WARNING, "marbl_forcing_init called with an associated control structure.")
      return
    endif

    allocate(CS)
    CS%diag => diag

    CS%use_MARBL_tracers = .true.
    if (.not. use_MARBL_tracers) then
      CS%use_MARBL_tracers = .false.
      return
    endif

    call get_param(param_file, mdl, "DUST_RATIO_THRES", CS%dust_ratio_thres, &
        "coarse/fine dust ratio threshold", units="1", default=90.)
    call get_param(param_file, mdl, "ATM_FE_TO_BC_RATIO", CS%atm_fe_to_bc_ratio, &
        "atmospheric iron to black carbon ratio", units="1", default=1.33)
    call get_param(param_file, mdl, "SEAICE_FE_TO_BC_RATIO", CS%seaice_fe_to_bc_ratio, &
        "sea-ice iron to black carbon ratio", units="1", default=1.33)
    call get_param(param_file, mdl, "IRON_FRAC_IN_ATM_FINE_DUST", CS%iron_frac_in_atm_fine_dust, &
        "Fraction of fine dust from the atmosphere that is iron", units="1", default=0.035)
    call get_param(param_file, mdl, "IRON_FRAC_IN_ATM_COARSE_DUST", CS%iron_frac_in_atm_coarse_dust, &
        "Fraction of coarse dust from the atmosphere that is iron", units="1", default=0.035)
    call get_param(param_file, mdl, "IRON_FRAC_IN_SEAICE_DUST", CS%iron_frac_in_seaice_dust, &
        "Fraction of dust from sea ice that is iron", units="1", default=0.035)
    call get_param(param_file, mdl, "ATM_CO2_OPT", atm_co2_opt, &
        "Source of atmospheric CO2 [constant, diagnostic, or prognostic]", &
        default="constant")
    select case (trim(atm_co2_opt))
      case("prognostic")
        CS%atm_co2_iopt = atm_co2_prognostic_iopt
      case("diagnostic")
        CS%atm_co2_iopt = atm_co2_diagnostic_iopt
      case("constant")
        CS%atm_co2_iopt = atm_co2_constant_iopt
      case DEFAULT
        write(err_message, "(3A)") "'", trim(atm_co2_opt), "' is not a valid ATM_CO2_OPT value"
        call MOM_error(FATAL, err_message)
    end select
    if (CS%atm_co2_iopt == atm_co2_constant_iopt) then
      call get_param(param_file, mdl, "ATM_CO2_CONST", CS%atm_co2_const, &
          "Value to send to MARBL as xco2", &
          default=284.317, units="ppm")
    endif
    call get_param(param_file, mdl, "ATM_ALT_CO2_OPT", atm_co2_opt, &
        "Source of alternate atmospheric CO2 [constant, diagnostic, or prognostic]", &
        default="constant")
    select case (trim(atm_co2_opt))
      case("prognostic")
        CS%atm_alt_co2_iopt = atm_co2_prognostic_iopt
      case("diagnostic")
        CS%atm_alt_co2_iopt = atm_co2_diagnostic_iopt
      case("constant")
        CS%atm_alt_co2_iopt = atm_co2_constant_iopt
      case DEFAULT
        write(err_message, "(3A)") "'", trim(atm_co2_opt), "' is not a valid ATM_ALT_CO2_OPT value"
        call MOM_error(FATAL, err_message)
    end select
    if (CS%atm_alt_co2_iopt == atm_co2_constant_iopt) then
      call get_param(param_file, mdl, "ATM_ALT_CO2_CONST", CS%atm_alt_co2_const, &
          "Value to send to MARBL as xco2_alt_co2", &
          default=284.317, units="ppm")
    endif

  end subroutine MARBL_forcing_init

  ! Note: ice fraction and u10_sqr are handled in mom_surface_forcing because of CFCs
  subroutine convert_driver_fields_to_forcings(atm_fine_dust_flux, atm_coarse_dust_flux, &
                                               seaice_dust_flux, atm_bc_flux, seaice_bc_flux, &
                                               nhx_dep, noy_dep, atm_co2_prog, atm_co2_diag, &
                                               afracr, swnet_afracr, ifrac_n, &
                                               swpen_ifrac_n, Time, G, US, i0, j0, fluxes, CS)

    real, dimension(:,:),   pointer, intent(in)    :: atm_fine_dust_flux   !< atmosphere fine dust flux from IOB
                                                                           !! [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: atm_coarse_dust_flux !< atmosphere coarse dust flux from IOB
                                                                           !! [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: seaice_dust_flux     !< sea ice dust flux from IOB [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: atm_bc_flux          !< atmosphere black carbon flux from IOB
                                                                           !! [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: seaice_bc_flux       !< sea ice black carbon flux from IOB
                                                                           !! [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: nhx_dep              !< NHx flux from atmosphere [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: noy_dep              !< NOy flux from atmosphere [kg m-2 s-1]
    real, dimension(:,:),   pointer, intent(in)    :: atm_co2_prog         !< Prognostic atmospheric CO2 concentration
                                                                           !! [ppm]
    real, dimension(:,:),   pointer, intent(in)    :: atm_co2_diag         !< Diagnostic atmospheric CO2 concentration
                                                                           !! [ppm]
    real, dimension(:,:),   pointer, intent(in)    :: afracr               !< open ocean fraction [1]
    real, dimension(:,:),   pointer, intent(in)    :: swnet_afracr         !< shortwave flux * open ocean fraction
                                                                           !! [W m-2]
    real, dimension(:,:,:), pointer, intent(in)    :: ifrac_n              !< per-category ice fraction [1]
    real, dimension(:,:,:), pointer, intent(in)    :: swpen_ifrac_n        !< per-category shortwave flux * ice fraction
                                                                           !! [W m-2]
    type(time_type),                 intent(in)    :: Time                 !< The time of the fluxes, used for
                                                                           !! interpolating the salinity to the
                                                                           !! right time, when it is being
                                                                           !! restored.
    type(ocean_grid_type),           intent(in)    :: G                    !< The ocean's grid structure
    type(unit_scale_type),           intent(in)    :: US                   !< A dimensional unit scaling type
    integer,                         intent(in)    :: i0                   !< i index offset
    integer,                         intent(in)    :: j0                   !< j index offset
    type(forcing),                   intent(inout) :: fluxes               !< MARBL-specific forcing fields
    type(marbl_forcing_CS), pointer, intent(inout) :: CS                   !< A pointer that is set to point to
                                                                           !! control structure for MARBL forcing

    integer :: i, j, is, ie, js, je, m
    real :: atm_fe_bioavail_frac     !< Fraction of iron from the atmosphere available for biological uptake [1]
    real :: dust_ratio               !< Ratio of coarse to fine dust from the atmosphere [1]
    real :: seaice_fe_bioavail_frac  !< Fraction of iron from sea ice available for biological uptake [1]
    ! Note: following two conversion factors are used to both convert from km m-2 s-1 -> mmol m-2 s-1
    !!      AND cast in MOM6's unique dimensional consistency scaling system [conc Z T-1]
    real :: iron_flux_conversion     !< Factor to convert iron flux from kg m-2 s-1 -> mmol m-3 (m s-1)
                                     !! [s m2 kg-1 conc Z T-1 ~> mmol kg-1]
    real :: ndep_conversion          !< Factor to convert nitrogen deposition from kg m-2 s-1 -> mmol m-3 (m s-1)
                                     !! [s m2 kg-1 conc Z T-1 ~> mmol kg-1]

    if (.not. CS%use_MARBL_tracers) return

    is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
    ndep_conversion = (1.e6/14.) * (US%m_to_Z * US%T_to_s)
    iron_flux_conversion = (1.e6 / molw_Fe) * (US%m_to_Z * US%T_to_s)

    do j=js,je ; do i=is,ie
      ! Components of dust flux
      fluxes%atm_fine_dust_flux(i,j) = (G%mask2dT(i,j) * US%kg_m2s_to_RZ_T) * atm_fine_dust_flux(i-i0,j-j0)
      fluxes%atm_coarse_dust_flux(i,j) = (G%mask2dT(i,j) * US%kg_m2s_to_RZ_T) * atm_coarse_dust_flux(i-i0,j-j0)
      fluxes%seaice_dust_flux(i,j) = (G%mask2dT(i,j) * US%kg_m2s_to_RZ_T) * seaice_dust_flux(i-i0,j-j0)

      ! Components of black carbon flux
      fluxes%atm_bc_flux(i,j) = (G%mask2dT(i,j) * US%kg_m2s_to_RZ_T) * atm_bc_flux(i-i0,j-j0)
      fluxes%seaice_bc_flux(i,j) = (G%mask2dT(i,j) * US%kg_m2s_to_RZ_T) * seaice_bc_flux(i-i0,j-j0)

      ! Nitrogen Deposition
      fluxes%nhx_dep(i,j) = (G%mask2dT(i,j) * ndep_conversion) * nhx_dep(i-i0,j-j0)
      fluxes%noy_dep(i,j) = (G%mask2dT(i,j) * ndep_conversion) * noy_dep(i-i0,j-j0)
    enddo ; enddo

    ! Atmospheric CO2
    select case (CS%atm_co2_iopt)
      case (atm_co2_prognostic_iopt)
        if (associated(atm_co2_prog)) then
          do j=js,je ; do i=is,ie
            fluxes%atm_co2(i,j) = G%mask2dT(i,j) * atm_co2_prog(i-i0,j-j0)
          enddo ; enddo
        else
          call MOM_error(FATAL, &
              "ATM_CO2_OPT = 'prognostic' but atmosphere is not providing this field")
        endif
      case (atm_co2_diagnostic_iopt)
        if (associated(atm_co2_diag)) then
          do j=js,je ; do i=is,ie
            fluxes%atm_co2(i,j) = G%mask2dT(i,j) * atm_co2_diag(i-i0,j-j0)
          enddo ; enddo
        else
          call MOM_error(FATAL, &
              "ATM_CO2_OPT = 'diagnostic' but atmosphere is not providing this field")
        endif
      case (atm_co2_constant_iopt)
        do j=js,je ; do i=is,ie
          fluxes%atm_co2(i,j) = G%mask2dT(i,j) * CS%atm_co2_const
        enddo ; enddo
    end select

    ! Alternate Atmospheric CO2
    select case (CS%atm_alt_co2_iopt)
      case (atm_co2_prognostic_iopt)
        if (associated(atm_co2_prog)) then
          do j=js,je ; do i=is,ie
            fluxes%atm_alt_co2(i,j) = G%mask2dT(i,j) * atm_co2_prog(i-i0,j-j0)
          enddo ; enddo
        else
          call MOM_error(FATAL, &
              "ATM_ALT_CO2_OPT = 'prognostic' but atmosphere is not providing this field")
        endif
      case (atm_co2_diagnostic_iopt)
        if (associated(atm_co2_diag)) then
          do j=js,je ; do i=is,ie
            fluxes%atm_alt_co2(i,j) = G%mask2dT(i,j) * atm_co2_diag(i-i0,j-j0)
          enddo ; enddo
        else
          call MOM_error(FATAL, &
              "ATM_ALT_CO2_OPT = 'diagnostic' but atmosphere is not providing this field")
        endif
      case (atm_co2_constant_iopt)
        do j=js,je ; do i=is,ie
          fluxes%atm_alt_co2(i,j) = G%mask2dT(i,j) * CS%atm_co2_const
        enddo ; enddo
    end select

    ! Dust flux
    if (associated(atm_fine_dust_flux)) then
      do j=js,je ; do i=is,ie
        fluxes%dust_flux(i,j) = (US%kg_m2s_to_RZ_T * G%mask2dT(i,j)) * &
            ((atm_fine_dust_flux(i-i0,j-j0) + atm_coarse_dust_flux(i-i0,j-j0)) + &
            seaice_dust_flux(i-i0,j-j0))
      enddo ; enddo
    endif

    if (associated(atm_bc_flux)) then
      do j=js,je ; do i=is,ie
        ! TODO: abort if atm_fine_dust_flux and atm_coarse_dust_flux are not associated?
        ! Contribution of atmospheric dust to iron flux
        atm_fe_bioavail_frac = 0.005
        if ((atm_coarse_dust_flux(i-i0,j-j0) > 0.) .and. (atm_fine_dust_flux(i-i0,j-j0)) > 0.) then
          dust_ratio = max(atm_coarse_dust_flux(i-i0,j-j0) / atm_fine_dust_flux(i-i0,j-j0), 9.903)
        else
          dust_ratio = 9.903
        endif
        dust_ratio = dust_ratio - 5.5
        if (dust_ratio < CS%dust_ratio_thres) &
          atm_fe_bioavail_frac = dust_ratio**(-0.9) - 0.0134

        ! Contribution of atmospheric dust to iron flux
        fluxes%iron_flux(i,j) = (atm_fe_bioavail_frac * &
            (CS%iron_frac_in_atm_fine_dust * atm_fine_dust_flux(i-i0,j-j0) + &
            CS%iron_frac_in_atm_coarse_dust * atm_coarse_dust_flux(i-i0,j-j0)))

        ! Contribution of atmospheric black carbon to iron flux
        fluxes%iron_flux(i,j) = fluxes%iron_flux(i,j) + (atm_bc_flux(i-i0,j-j0) * &
            (atm_fe_bioavail_frac * CS%atm_fe_to_bc_ratio))

        seaice_fe_bioavail_frac = atm_fe_bioavail_frac
        ! Contribution of seaice dust to iron flux
        fluxes%iron_flux(i,j) = fluxes%iron_flux(i,j) + (seaice_fe_bioavail_frac * &
            (CS%iron_frac_in_seaice_dust * seaice_dust_flux(i-i0,j-j0)))

        ! Contribution of seaice black carbon to iron flux
        fluxes%iron_flux(i,j) = fluxes%iron_flux(i,j) + (seaice_bc_flux(i-i0,j-j0) * &
            (seaice_fe_bioavail_frac * CS%seaice_fe_to_bc_ratio))

        ! Unit conversion (kg m-2 s-1 -> conc Z T-1)
        fluxes%iron_flux(i,j) = (G%mask2dT(i,j) * iron_flux_conversion) * fluxes%iron_flux(i,j)

      enddo ; enddo
    endif

      ! Per ice-category forcings
      ! If the cap receives per-category fields, memory should be allocated in fluxes
    if (associated(ifrac_n)) then
      do j=js,je ; do i=is,ie
        fluxes%fracr_cat(i,j,1) = min(1., afracr(i-i0,j-j0))
        fluxes%qsw_cat(i,j,1) = swnet_afracr(i-i0,j-j0)
        do m=1,size(ifrac_n, 3)
          fluxes%fracr_cat(i,j,m+1) = min(1., ifrac_n(i-i0,j-j0,m))
          fluxes%qsw_cat(i,j,m+1)   = swpen_ifrac_n(i-i0,j-j0,m)
        enddo
        where (fluxes%fracr_cat(i,j,:) > 0.)
          fluxes%qsw_cat(i,j,:) = fluxes%qsw_cat(i,j,:) / fluxes%fracr_cat(i,j,:)
        elsewhere
          fluxes%fracr_cat(i,j,:) = 0.
          fluxes%qsw_cat(i,j,:) = 0.
        endwhere
        fluxes%fracr_cat(i,j,:) = G%mask2dT(i,j) * fluxes%fracr_cat(i,j,:)
        fluxes%qsw_cat(i,j,:)   = (US%W_m2_to_QRZ_T * G%mask2dT(i,j)) * fluxes%qsw_cat(i,j,:)
      enddo ; enddo
    endif

  end subroutine convert_driver_fields_to_forcings

end module MARBL_forcing_mod
