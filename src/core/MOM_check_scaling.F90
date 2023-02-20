!> This module is used to check the dimensional scaling factors used by the MOM6 ocean model
module MOM_check_scaling

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler,        only : MOM_error, MOM_mesg, FATAL, WARNING, assert, MOM_get_verbosity
use MOM_unique_scales,        only : check_scaling_uniqueness, scales_to_powers
use MOM_unit_scaling,         only : unit_scale_type
use MOM_verticalGrid,         only : verticalGrid_type

implicit none ; private

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

public check_MOM6_scaling_factors

contains

!> Evaluate whether the dimensional scaling factors provide unique tests for all of the combinations
!! of dimensions that are used in MOM6 (or perhaps widely used), and if they are not unique, explore
!! whether another combination of scaling factors can be found that is unique or has less common
!! cases with coinciding scaling.
subroutine check_MOM6_scaling_factors(GV, US)
  type(verticalGrid_type), pointer    :: GV         !< The container for vertical grid data
  type(unit_scale_type),   intent(in) :: US         !< A dimensional unit scaling type

  ! Local variables
  integer, parameter :: ndims = 8 ! The number of rescalable dimensional factors.
  real,    dimension(ndims) :: scales ! An array of scaling factors for each of the basic units.
  integer, dimension(ndims) :: scale_pow2 ! The powers of 2 that give each element of scales.
  character(len=2), dimension(ndims) :: key
  integer, allocatable :: weights(:)
  character(len=80), allocatable :: descriptions(:)
  integer :: n, ns, max_pow

  ! If no scaling is being done, simply return.
  if ((US%Z_to_m == 1.) .and. (GV%H_to_MKS == 1.) .and. (US%L_to_m == 1.) .and. &
      (US%T_to_s == 1.) .and. (US%R_to_kg_m3 == 1.) .and. (US%Q_to_J_kg == 1.) .and. &
      (US%C_to_degC == 1.) .and. (US%S_to_ppt == 1.)) return

  ! Set the names and scaling factors of the dimensions being rescaled.
  key(:) = ["Z", "H", "L", "T", "R", "Q", "C", "S"]
  scales(:) = (/ US%Z_to_m, GV%H_to_MKS, US%L_to_m, US%T_to_s, US%R_to_kg_m3, US%Q_to_J_kg, &
                 US%C_to_degC, US%S_to_ppt/)
  call scales_to_powers(scales, scale_pow2)
  max_pow = 40 ! 60

  ! The first call is just to find out how many elements are in the list of scaling combinations.
  call compose_dimension_list(ns, descriptions, weights)

  allocate(descriptions(ns))
  do n=1,ns ; descriptions(n) = "" ; enddo
  allocate(weights(ns), source=0)
  ! This call records all the list of powers, the descriptions, and their weights.
  call compose_dimension_list(ns, descriptions, weights)

  call check_scaling_uniqueness("MOM6", descriptions, weights, key, scale_pow2, max_pow)

  deallocate(weights)
  deallocate(descriptions)

end subroutine check_MOM6_scaling_factors


!> This routine composes a list of the commonly used dimensional scaling factors in the MOM6
!! code, along with weights reflecting the frequency of their occurrence in the MOM6 code or
!! other considerations of how likely the variables are be used.
subroutine compose_dimension_list(ns, des, wts)
  integer,                       intent(out)   :: ns     !< The running sum of valid descriptions
  character(len=*), allocatable, intent(inout) :: des(:) !< The unit descriptions that have been converted
  integer,          allocatable, intent(inout) :: wts(:) !< A list of the integer weights for each scaling factor,
                                                         !! perhaps the number of times it occurs in the MOM6 code.

  ns = 0
  ! Accumulate a list of units used in MOM6, in approximate descending order of frequency of occurrence in
  ! doxygen comments (i.e., arguments and elements in types), excluding the code in the user, ice_shelf and
  ! framework directories and the passive tracer packages.
  call add_scaling(ns, des, wts, "[H ~> m or kg m-2]", 716)   ! Layer thicknesses
  call add_scaling(ns, des, wts, "[L T-1 ~> m s-1]", 264)     ! Horizontal velocities
  call add_scaling(ns, des, wts, "[Z ~> m]", 244)             ! Depths and vertical distance
  call add_scaling(ns, des, wts, "[T ~> s]", 154)             ! Time intervals
  call add_scaling(ns, des, wts, "[S ~> ppt]", 135)           ! Salinities
  call add_scaling(ns, des, wts, "[C ~> degC]", 135)          ! Temperatures
  call add_scaling(ns, des, wts, "[R L2 T-2 ~> Pa]", 133)     ! Dynamic pressure
  ! call add_scaling(ns, des, wts, "[R L2 T-2 ~> J m-3]")     ! Energy density
  call add_scaling(ns, des, wts, "[Z2 T-1 ~> m2 s-1]", 132)   ! Vertical viscosities and diffusivities
  call add_scaling(ns, des, wts, "[R ~> kg m-3]", 122)        ! Densities

  call add_scaling(ns, des, wts, "[H L2 T-1 ~> m3 s-1 or kg s-1]", 97) ! Volume or mass transports
  call add_scaling(ns, des, wts, "[H L2 ~> m3 or kg]", 91)    ! Cell volumes or masses
  call add_scaling(ns, des, wts, "[L T-2 ~> m s-2]", 82)      ! Horizontal accelerations
  call add_scaling(ns, des, wts, "[T-1 ~> s-1]", 67)          ! Rates
  call add_scaling(ns, des, wts, "[Z T-1 ~> m s-1]", 56)      ! Friction velocities and viscous coupling
  call add_scaling(ns, des, wts, "[Q R Z T-1 ~> W m-2]", 42)  ! Vertical heat fluxes
  call add_scaling(ns, des, wts, "[L2 T-1 ~> m2 s-1]", 45)    ! Horizontal viscosity or diffusivity
  call add_scaling(ns, des, wts, "[L2 T-2 ~> m2 s-2]", 37)    ! Resolved kinetic energy per unit mass
  call add_scaling(ns, des, wts, "[L ~> m]", 35)              ! Horizontal distances
  call add_scaling(ns, des, wts, "[T-2 ~> s-2]", 33)          ! Squared shears and buoyancy frequency

  call add_scaling(ns, des, wts, "[R Z L T-2 ~> Pa]", 33)     ! Wind stresses
  call add_scaling(ns, des, wts, "[H L ~> m2 or kg m-1]", 32) ! Lateral cell face areas
  call add_scaling(ns, des, wts, "[L2 ~> m2]", 31)            ! Horizontal areas
  call add_scaling(ns, des, wts, "[R C-1 ~> kg m-3 degC-1]", 26) ! Thermal expansion coefficients
  call add_scaling(ns, des, wts, "[L2 Z-1 T-2 ~> m s-2]", 26) ! Gravitational acceleration
  call add_scaling(ns, des, wts, "[R S-1 ~> kg m-3 ppt-1]", 23) ! Haline contraction coefficients
  call add_scaling(ns, des, wts, "[R Z3 T-3 ~> W m-2]", 23)   ! Integrated turbulent kinetic energy sources
  call add_scaling(ns, des, wts, "[R Z T-1 ~> kg m-2 s-1]", 19) ! Vertical mass fluxes
  call add_scaling(ns, des, wts, "[C H ~> degC m or degC kg m-2]", 17) ! Heat content
  call add_scaling(ns, des, wts, "[H-1 ~> m-1 or m2 kg-1]", 17) ! Inverse cell thicknesses

  call add_scaling(ns, des, wts, "[Z-1 ~> m-1]", 14)          ! Inverse vertical distances
  call add_scaling(ns, des, wts, "[R-1 ~> m3 kg-1]", 14)      ! Specific volumes
  call add_scaling(ns, des, wts, "[Z L-1 ~> nondim]", 12)     ! Slopes
  call add_scaling(ns, des, wts, "[L-1 ~> m-1]", 12)          ! Inverse horizontal distances
  call add_scaling(ns, des, wts, "[L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1]", 12) ! pbce or gtot
  call add_scaling(ns, des, wts, "[R Z ~> kg m-2]", 11)       ! Layer or column mass loads
  call add_scaling(ns, des, wts, "[Z L2 T-2 ~> m3 s-2]", 11)  ! Integrated energy per unit mass
  call add_scaling(ns, des, wts, "[R Z3 T-2 ~> J m-2]", 11)   ! Integrated turbulent kinetic energy density
  call add_scaling(ns, des, wts, "[H T-1 ~> m s-1 or kg m-2 s-1]", 9) ! Vertical thickness fluxes
  call add_scaling(ns, des, wts, "[L-1 T-1 ~> m-1 s-1]", 9)   ! Laplacian of velocity

  call add_scaling(ns, des, wts, "[Z3 T-3 ~> m3 s-3]", 9)     ! Integrated turbulent kinetic energy sources
  call add_scaling(ns, des, wts, "[S H ~> ppt m or ppt kg m-2]", 8) ! Depth integrated salinity
  call add_scaling(ns, des, wts, "[Z2 T-2 ~> m2 s-2]", 8)     ! Turbulent kinetic energy
  call add_scaling(ns, des, wts, "[R L2 Z T-2 ~> Pa m]", 7)   ! Vertically integrated pressure anomalies
  call add_scaling(ns, des, wts, "[Z2 T-1 / Z3 T-3 = T2 Z-1 ~> s2 m-1]", 7) ! (TKE_to_Kd)
  call add_scaling(ns, des, wts, "[L4 T-1 ~> m4 s-1]", 7)     ! Biharmonic viscosity
  call add_scaling(ns, des, wts, "[L3 ~> m3]", 7)             ! Metric dependent constants for viscosity
  call add_scaling(ns, des, wts, "[L2 T-3 ~> m2 s-3]", 7)     ! Buoyancy flux or MEKE sources [L2 T-3 ~> W kg-1]
  call add_scaling(ns, des, wts, "[H2 ~> m2 or kg2 m-4]", 7)  ! Squared layer thicknesses
  call add_scaling(ns, des, wts, "[C H T-1 ~> degC m s-1 or degC kg m-2 s-1]", 7) ! vertical heat fluxes

  call add_scaling(ns, des, wts, "[L-2 ~> m-2]", 6)           ! Inverse areas
  call add_scaling(ns, des, wts, "[R Z L2 T-3 ~> W m-2]", 6)  ! Energy sources, including for MEKE
  call add_scaling(ns, des, wts, "[Z2 T-3 ~> m2 s-3]", 5)     ! Certain buoyancy fluxes
  call add_scaling(ns, des, wts, "[Z2 ~> m2]", 5)             ! Squared vertical distances
  call add_scaling(ns, des, wts, "[S H T-1 ~> ppt m s-1 or ppt kg m-2 s-1]", 5) ! vertical salinity fluxes
  call add_scaling(ns, des, wts, "[R-1 C-1 ~> m3 kg-1 degC-1]", 5) ! Specific volume temperature gradient
  call add_scaling(ns, des, wts, "[R-1 S-1 ~> m3 kg-1 ppt-1]", 4) ! Specific volume salnity gradient
  call add_scaling(ns, des, wts, "[Q R Z ~> J m-2]", 4)       ! time-integrated frazil heat flux
  call add_scaling(ns, des, wts, "[Z C-1 ~> m degC-1]", 4)    ! Inverse temperature gradients
  call add_scaling(ns, des, wts, "[Z S-1 ~> m ppt-1]", 4)     ! Inverse salinity gradients

  call add_scaling(ns, des, wts, "[R Z3 T-2 H-1 ~> J m-3 or J kg-1]", 4) ! Partial derivatives of energy
  call add_scaling(ns, des, wts, "[R Z3 T-2 S-1 ~> J m-2 ppt-1]", 4) ! Sensitity of energy change to salinity
  call add_scaling(ns, des, wts, "[R Z3 T-2 C-1 ~> J m-2 degC-1]", 4) ! Sensitity of energy change to temperature
  call add_scaling(ns, des, wts, "[R L4 T-4 ~> Pa m2 s-2]", 4) ! Integral in geopotential of pressure
  call add_scaling(ns, des, wts, "[Q ~> J kg-1]", 4)          ! Latent heats
  call add_scaling(ns, des, wts, "[Q C-1 ~> J kg-1 degC-1]", 4) ! Heat capacity
  call add_scaling(ns, des, wts, "[L-3 ~> m-3]", 4)           ! Metric dependent constants for viscosity
  call add_scaling(ns, des, wts, "[L2 Z-2 T-2 ~> s-2]", 4)    ! Buoyancy frequency in some params.
  call add_scaling(ns, des, wts, "[H R ~> kg m-2 or kg2 m-5]", 4) ! Layer-integrated density
  call add_scaling(ns, des, wts, "[H L T-1 ~> m2 s-1 or kg m-1 s-1]", 4) ! Layer integrated velocities

  call add_scaling(ns, des, wts, "[H T2 L-1 ~> s2 or kg s2 m-3]", 4) ! BT_cont_type face curvature fit
  call add_scaling(ns, des, wts, "[H L-1 ~> nondim or kg m-3]", 4)   ! BT_cont_type face curvature fit
  call add_scaling(ns, des, wts, "[C2 ~> degC2]", 4)          ! Squared temperature anomalies
  call add_scaling(ns, des, wts, "[S2 ~> ppt2]", 3)           ! Squared salinity anomalies
  call add_scaling(ns, des, wts, "[C S ~> degC ppt]", 3)      ! Covariance of temperature and salinity anomalies
  call add_scaling(ns, des, wts, "[S R Z ~> gSalt m-2]", 3)   ! Total ocean column salt
  call add_scaling(ns, des, wts, "[C R Z ~> degC kg m-2]", 3) ! Total ocean column temperature
  call add_scaling(ns, des, wts, "[Pa T2 R-1 L-2 ~> 1]", 3)   ! Pressure conversions
  call add_scaling(ns, des, wts, "[Z H-1 ~> nondim or m3 kg-1]", 3) ! Thickness to height conversion
  call add_scaling(ns, des, wts, "[R Z2 T-2 ~> J m-3]", 3)    ! Potential energy height derivatives

  call add_scaling(ns, des, wts, "[H-2 ~> m-2 or m4 kg-2]", 3) ! Mixed layer local work variables
  call add_scaling(ns, des, wts, "[C S-1 ~> degC ppt-1]", 2)  ! T / S gauge transformation
  call add_scaling(ns, des, wts, "[R S-2 ~> kg m-3 ppt-2]", 2) ! Second derivative of density
  call add_scaling(ns, des, wts, "[R C-2 ~> kg m-3 degC-2]", 2) ! Second derivative of density
  call add_scaling(ns, des, wts, "[R S-1 C-1 ~> kg m-3 ppt-1 degC-1]", 2) ! Second derivative of density
  call add_scaling(ns, des, wts, "[T2 S-1 L-2 ~> kg m-3 ppt-1 Pa-1]", 2)  ! Second derivative of density
  call add_scaling(ns, des, wts, "[T2 C-1 L-2 ~> kg m-3 degC-1 Pa-1]", 2) ! Second derivative of density
  call add_scaling(ns, des, wts, "[T2 L-2 ~> s2 m-2]", 2)     ! Inverse velocities squared
  call add_scaling(ns, des, wts, "[R Z2 T-3 ~> W m-3]", 2)    ! Kinetic energy dissipation rates
  call add_scaling(ns, des, wts, "[R H-1 ~> kg m-4 or m-1]", 2) ! Vertical density gradients

  call add_scaling(ns, des, wts, "[L4 ~> m4]", 2)             ! Metric dependent constants for viscosity
  call add_scaling(ns, des, wts, "[Z L T-1 ~> m2 s-1]", 2)    ! Layer integrated velocities
  call add_scaling(ns, des, wts, "[C Z ~> degC m]", 2)        ! Depth integrated temperature
  call add_scaling(ns, des, wts, "[S Z ~> ppt m]", 1)         ! Layer integrated salinity
  call add_scaling(ns, des, wts, "[T L4 ~> s m4]", 2)         ! Biharmonic metric dependent constant
  call add_scaling(ns, des, wts, "[L6 ~> m6]", 2)             ! Biharmonic Leith metric dependent constant
  call add_scaling(ns, des, wts, "[L4 Z-1 T-1 ~> m3 s-1]", 2) ! Rigidity of ice
  call add_scaling(ns, des, wts, "[L4 Z-2 T-1 ~> m2 s-1]", 1) ! Ice rigidity term
  call add_scaling(ns, des, wts, "[R-1 Z-1 ~> m2 kg-1]", 1)   ! Inverse of column mass
  call add_scaling(ns, des, wts, "[Z-2 ~> m-2]", 1)           ! Inverse of denominator in some weighted averages

  call add_scaling(ns, des, wts, "[R Z2 T-1 ~> J s m-3]", 1)  ! River mixing term
  call add_scaling(ns, des, wts, "[R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1]", 1) ! Thickness to pressure conversion
  call add_scaling(ns, des, wts, "[Z T2 R-1 L-2 ~> m Pa-1]", 1) ! Atmospheric pressure SSH correction
  call add_scaling(ns, des, wts, "[T Z ~> s m] ", 1)          ! Time integrated SSH
  call add_scaling(ns, des, wts, "[Z-1 T-1 ~> m-1 s-1]", 1)   ! barotropic PV
  call add_scaling(ns, des, wts, "[L2 T ~> m2 s]", 1)         ! Greatbatch & Lamb 90 coefficient
  call add_scaling(ns, des, wts, "[Z L2 T-1 ~> m3 s-1]", 1)   ! Overturning (GM) streamfunction
  call add_scaling(ns, des, wts, "[kg H-1 L-2 ~> kg m-3 or 1]", 1) ! Diagnostic conversions to mass
  call add_scaling(ns, des, wts, "[S-1 ~> ppt-1]", 1)         ! Unscaling salinity
  call add_scaling(ns, des, wts, "[C-1 ~> degC-1]", 1)        ! Unscaling temperature

  call add_scaling(ns, des, wts, "[R Z H-1 ~> kg m-3 or 1] ", 1) ! A unit conversion factor
  call add_scaling(ns, des, wts, "[H R-1 Z-1 ~> m3 kg-2 or 1]", 1) ! A unit conversion factor
  call add_scaling(ns, des, wts, "[H Z-1 ~> 1 or kg m-3]", 1) ! A unit conversion factor
  call add_scaling(ns, des, wts, "[m T s-1 L-1 ~> 1]", 1)     ! A unit conversion factor

end subroutine compose_dimension_list

!> Augment the count the valid unit descriptions, and add the provided description and its weight
!! to the end of the list if that list is allocated.
subroutine add_scaling(ns, descs, wts, scaling, weight)
  integer,                       intent(inout) :: ns       !< The running sum of valid descriptions.
  character(len=*), allocatable, intent(inout) :: descs(:) !< The unit descriptions that have been converted
  integer,          allocatable, intent(inout) :: wts(:)   !< A list of the integers for each scaling
  character(len=*),              intent(in)    :: scaling  !< The unit description that will be converted
  integer,             optional, intent(in)    :: weight   !< An optional weight or occurrence count
                                                           !! for this unit description, 1 by default.

  integer :: iend

  iend = index(scaling, "~>")
  if (iend <= 1) then
    call MOM_mesg("No scaling indicator ~> found for "//trim(scaling))
  else
    ! Count and perhaps store this description and its weight.
    ns = ns + 1
    if (allocated(descs)) descs(ns) = scaling
    if (allocated(wts)) then
      wts(ns) = 1 ; if (present(weight)) wts(ns) = weight
    endif
  endif

end subroutine add_scaling

end module MOM_check_scaling
