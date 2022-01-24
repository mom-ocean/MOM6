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
  integer, parameter :: ndims = 6 ! The number of rescalable dimensional factors.
  real,    dimension(ndims) :: scales ! An array of scaling factors for each of the basic units.
  integer, dimension(ndims) :: scale_pow2 ! The powers of 2 that give each element of scales.
  character(len=2), dimension(ndims) :: key
  ! character(len=128) :: mesg, msg_frag
  integer, allocatable :: weights(:)
  character(len=80), allocatable :: descriptions(:)
  ! logical :: verbose, very_verbose
  integer :: n, ns, max_pow

  ! Set the names and scaling factors of the dimensions being rescaled.
  key(:) = ["Z", "H", "L", "T", "R", "Q"]
  scales(:) = (/ US%Z_to_m, GV%H_to_MKS, US%L_to_m, US%T_to_s, US%R_to_kg_m3, US%Q_to_J_kg /)
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
  ! Accumulate a list of units used in MOM6, in approximate descending order of frequency of occurrence.
  call add_scaling(ns, des, wts, "[H ~> m or kg m-2]", 1239)  ! Layer thicknesses
  call add_scaling(ns, des, wts, "[Z ~> m]", 660)             ! Depths and vertical distance
  call add_scaling(ns, des, wts, "[L T-1 ~> m s-1]", 506)     ! Horizontal velocities
  call add_scaling(ns, des, wts, "[R ~> kg m-3]", 356)        ! Densities
  call add_scaling(ns, des, wts, "[T-1 ~> s-1]", 247)         ! Rates
  call add_scaling(ns, des, wts, "[T ~> s]", 237)             ! Time intervals
  call add_scaling(ns, des, wts, "[R L2 T-2 ~> Pa]", 231)     ! Dynamic pressure
  ! call add_scaling(ns, des, wts, "[R L2 T-2 ~> J m-3]")     ! Energy density
  call add_scaling(ns, des, wts, "[Z2 T-1 ~> m2 s-1]", 181)   ! Vertical viscosities and diffusivities
  call add_scaling(ns, des, wts, "[H L2 ~> m3 or kg]", 174)   ! Cell volumes or masses
  call add_scaling(ns, des, wts, "[H L2 T-1 ~> m3 s-1 or kg s-1]", 163) ! Volume or mass transports
  call add_scaling(ns, des, wts, "[L T-2 ~> m s-2]", 136)     ! Horizontal accelerations
  call add_scaling(ns, des, wts, "[L ~> m]", 107)             ! Horizontal distances
  call add_scaling(ns, des, wts, "[Z T-1 ~> m s-1]", 104)     ! Friction velocities and viscous coupling
  call add_scaling(ns, des, wts, "[H-1 ~> m-1 or m2 kg-1]", 89) ! Inverse cell thicknesses
  call add_scaling(ns, des, wts, "[L2 T-2 ~> m2 s-2]", 88)    ! Resolved kinetic energy per unit mass
  call add_scaling(ns, des, wts, "[R Z3 T-2 ~> J m-2]", 85)   ! Integrated turbulent kinetic energy density
  call add_scaling(ns, des, wts, "[L2 T-1 ~> m2 s-1]", 78)    ! Horizontal viscosity or diffusivity
  call add_scaling(ns, des, wts, "[T-2 ~> s-2]", 69)          ! Squared shears and buoyancy frequency
  call add_scaling(ns, des, wts, "[H L ~> m2 or kg m-1]", 68) ! Lateral cell face areas
  call add_scaling(ns, des, wts, "[L2 ~> m2]", 67)            ! Horizontal areas

  call add_scaling(ns, des, wts, "[R-1 ~> m3 kg-1]", 61)      ! Specific volumes
  call add_scaling(ns, des, wts, "[Q R Z T-1 ~> W m-2]", 62)  ! Vertical heat fluxes
  call add_scaling(ns, des, wts, "[Z-1 ~> m-1]", 60)          ! Inverse vertical distances
  call add_scaling(ns, des, wts, "[L2 Z-1 T-2 ~> m s-2]", 57) ! Gravitational acceleration
  call add_scaling(ns, des, wts, "[R Z T-1 ~> kg m-2 s-1]", 52) ! Vertical mass fluxes
  call add_scaling(ns, des, wts, "[H T-1 ~> m s-1 or kg m-2 s-1]", 51) ! Vertical thickness fluxes
  call add_scaling(ns, des, wts, "[R Z3 T-3 ~> W m-2]", 45)   ! Integrated turbulent kinetic energy sources
  call add_scaling(ns, des, wts, "[R Z ~> kg m-2]", 42)       ! Layer or column mass loads
  call add_scaling(ns, des, wts, "[Z3 T-3 ~> m3 s-3]", 33)    ! Integrated turbulent kinetic energy sources
  call add_scaling(ns, des, wts, "[H2 ~> m2 or kg2 m-4]", 35) ! Squared layer thicknesses
  call add_scaling(ns, des, wts, "[Z2 T-2 ~> m2 s-2]", 33)    ! Turbulent kinetic energy
  call add_scaling(ns, des, wts, "[L-1 ~> m-1]", 32)          ! Inverse horizontal distances
  call add_scaling(ns, des, wts, "[R L Z T-2 ~> Pa]", 27)     ! Wind stresses
  call add_scaling(ns, des, wts, "[T2 L-2 ~> s2 m-2]", 33)    ! Inverse velocities squared
  call add_scaling(ns, des, wts, "[R Z L2 T-2 ~> J m-2]", 25) ! Integrated energy
  ! call add_scaling(ns, des, wts, "[R L2 Z T-2 ~> Pa m]")    ! Depth integral of pressures (25)
  call add_scaling(ns, des, wts, "[Z L2 T-2 ~> m3 s-2]", 25)  ! Integrated energy
  call add_scaling(ns, des, wts, "[H R ~> kg m-2 or kg2 m-5]", 24) ! Layer-integrated density
  call add_scaling(ns, des, wts, "[L2 T-2 H-1 ~> m s-2 or m4 s-2 kg-1]", 20) ! pbce or gtot
  call add_scaling(ns, des, wts, "[L-1 T-1 ~> m-1 s-1]", 19)  ! Laplacian of velocity

  call add_scaling(ns, des, wts, "[L4 T-1 ~> m4 s-1]", 18)    ! Biharmonic viscosity
  call add_scaling(ns, des, wts, "[Z L T-1 ~> m2 s-1]", 17)   ! Layer integrated velocities
  call add_scaling(ns, des, wts, "[Z L-1 ~> nondim]", 15)     ! Slopes
  call add_scaling(ns, des, wts, "[Z L2 ~> m3]", 14)          ! Diagnostic volumes
  call add_scaling(ns, des, wts, "[H L T-1 ~> m2 s-1 or kg m-1 s-1]", 12) ! Layer integrated velocities
  call add_scaling(ns, des, wts, "[L2 T-3 ~> m2 s-3]", 14)    ! Buoyancy flux or MEKE sources [L2 T-3 ~> W kg-1]
  call add_scaling(ns, des, wts, "[Z2 ~> m2]", 12)            ! Squared vertical distances
  call add_scaling(ns, des, wts, "[R Z L2 T-1 ~> kg s-1]", 12) ! Mass fluxes
  call add_scaling(ns, des, wts, "[L-2 ~> m-2]", 12)          ! Inverse areas
  call add_scaling(ns, des, wts, "[L2 Z-1 T-2 R-1 ~> m4 s-2 kg-1]", 11) ! Gravitational acceleration over density
  call add_scaling(ns, des, wts, "[Z T-2 ~> m s-2]", 10)      ! Buoyancy differences or their derivatives
  ! Could also add [Z T-2 degC-1 ~> m s-2 degC-1] or [Z T-2 ppt-1 ~> m s-2 ppt-1]
  call add_scaling(ns, des, wts, "[R Z L2 T-3 ~> W m-2]", 10) ! Energy sources, including for MEKE
  call add_scaling(ns, des, wts, "[L3 ~> m3]", 10)            ! Metric dependent constants for viscosity
  call add_scaling(ns, des, wts, "[Z-2 ~> m-2]", 9)           ! Inverse of denominator in some weighted averages
  call add_scaling(ns, des, wts, "[H-2 ~> m-2 or m4 kg-2]", 9) ! Mixed layer local work variables
  call add_scaling(ns, des, wts, "[Z L2 T-1 ~> m3 s-1]", 9)   ! Overturning (GM) streamfunction
  call add_scaling(ns, des, wts, "[L2 Z-2 T-2 ~> s-2]", 9)    ! Buoyancy frequency in some params.
  call add_scaling(ns, des, wts, "[Q R Z ~> J m-2]", 8)       ! time-integrated frazil heat flux
  call add_scaling(ns, des, wts, "[Z2 T-1 / Z3 T-3 = T2 Z-1 ~> s2 m-1]", 7) ! (TKE_to_Kd)
  call add_scaling(ns, des, wts, "[Q degC-1 ~> J kg-1 degC-1]", 7) ! Heat capacity

  call add_scaling(ns, des, wts, "[R Z2 T-2 ~> J m-3]", 6)    ! Potential energy height derivatives
  call add_scaling(ns, des, wts, "[R Z3 T-2 H-1 ~> J m-3 or J kg-1]", 7) ! Partial derivatives of energy
  call add_scaling(ns, des, wts, "[R L2 T-2 Z-1 ~> Pa m-1]", 7) ! Converts depth to pressure
  call add_scaling(ns, des, wts, "[L4 Z-1 T-1 ~> m3 s-1]", 7) ! Rigidity of ice
  call add_scaling(ns, des, wts, "[H L2 T-3 ~> m3 s-3]", 9)   ! Kinetic energy diagnostics
  call add_scaling(ns, des, wts, "[H-1 T-1 ~> m-1 s-1 or m2 kg-1 s-1]", 6) ! Layer potential vorticity
  call add_scaling(ns, des, wts, "[R Z2 T-3 ~> W m-3]", 3)    ! Kinetic energy dissipation rates
  call add_scaling(ns, des, wts, "[Z2 L-2 ~> 1]", 1)          !  Slopes squared
  call add_scaling(ns, des, wts, "[Z H-1 ~> nondim or m3 kg-1]", 6) ! Thickness to height conversion
  call add_scaling(ns, des, wts, "[Pa T2 R-1 L-2 ~> 1]", 6)   ! Pressure conversion factor
    ! Could also add [m T2 R-1 L-2 ~> m Pa-1]
    ! Could also add [degC T2 R-1 L-2 ~> degC Pa-1]
  call add_scaling(ns, des, wts, "[R H-1 ~> kg m-4 or m-1]", 5) ! Vertical density gradients
  call add_scaling(ns, des, wts, "[R L4 T-4 ~> Pa m2 s-2]", 4) ! Integral in geopotential of pressure
  call add_scaling(ns, des, wts, "[L Z-1 ~> nondim]", 4)      ! Inverse slopes
  call add_scaling(ns, des, wts, "[L-3 ~> m-3]", 4)           ! Metric dependent constants for viscosity
  call add_scaling(ns, des, wts, "[H T2 L-1 ~> s2 or kg s2 m-3]", 4) ! BT_cont_type face curvature fit
  call add_scaling(ns, des, wts, "[H L-1 ~> nondim or kg m-3]", 4)   ! BT_cont_type face curvature fit
  call add_scaling(ns, des, wts, "[kg H-1 L-2 ~> kg m-3 or 1]", 20)  ! Diagnostic conversions to mass
    ! Could also add [m3 H-1 L-2 ~> 1 or m3 kg-1]
  call add_scaling(ns, des, wts, "[Z T-2 R-1 ~> m4 s-2 kg-1]", 9) ! Gravitational acceleration over density
  call add_scaling(ns, des, wts, "[R Z L4 T-3 ~> kg m2 s-3]", 9) ! MEKE fluxes
  call add_scaling(ns, des, wts, "[R L2 T-2 H-1 ~> Pa m-1 or Pa m2 kg-1]", 3) ! Thickness to pressure conversion

  call add_scaling(ns, des, wts, "[R-1 Z-1 ~> m2 kg-1]", 3)   ! Inverse of column mass
  call add_scaling(ns, des, wts, "[L4 ~> m4]", 3)             ! Metric dependent constants for viscosity
  call add_scaling(ns, des, wts, "[T-1 Z-1 ~> s-1 m-1]", 2)   ! Barotropic PV, for some options
  call add_scaling(ns, des, wts, "[R Z2 T-1 ~> J s m-3]", 2)  ! River mixing term [R Z2 T-1 ~> Pa s]
  call add_scaling(ns, des, wts, "[degC Q-1 ~> kg degC J-1]", 2) ! Inverse heat capacity
     ! Could add call add_scaling(ns, des, wts, "[Q-1 ~> kg J-1]", 1) ! Inverse heat content
  call add_scaling(ns, des, wts, "[L4 Z-2 T-1 ~> m2 s-1]", 2) ! Ice rigidity term
  call add_scaling(ns, des, wts, "[R Z-1 ~> kg m-4]", 3)      ! Vertical density gradient
  call add_scaling(ns, des, wts, "[R Z L2 ~> kg]", 3)         ! Depth and time integrated mass fluxes
  call add_scaling(ns, des, wts, "[R L2 T-3 ~> W m-2]", 3)    ! Depth integrated friction work
  call add_scaling(ns, des, wts, "[ppt2 R-2 ~> ppt2 m6 kg-2]", 3) ! T / S gauge transformation
  call add_scaling(ns, des, wts, "[R L-1 ~> kg m-4]", 2)      ! Horizontal density gradient
     ! Could add call add_scaling(ns, des, wts, "[H Z ~> m2 or kg m-1]", 2)  ! Temporary variables
  call add_scaling(ns, des, wts, "[Z3 R2 T-2 H-2 ~> kg2 m-5 s-2 or m s-2]", 2) ! Heating to PE change
  call add_scaling(ns, des, wts, "[R2 L2 Z2 T-4 ~> Pa2]", 2)  ! Squared wind stresses
  call add_scaling(ns, des, wts, "[L-2 T-2 ~> m-2 s-2]", 2)   ! Squared Laplacian of velocity
  call add_scaling(ns, des, wts, "[T H Z-1 ~> s or s kg m-3]", 2) ! Time step times thickness conversion
  call add_scaling(ns, des, wts, "[T H Z-1 R-1 ~> s m3 kg-1 or s]", 2) ! Time step over density with conversion
  call add_scaling(ns, des, wts, "[H-3 ~> m-3 or m6 kg-3]", 1) ! A local term in ePBL
  call add_scaling(ns, des, wts, "[H-4 ~> m-4 or m8 kg-4]", 1) ! A local term in ePBL
  call add_scaling(ns, des, wts, "[H T Z-2 ~> s m-1 or kg s m-4]", 1) ! A local term in ePBL

  call add_scaling(ns, des, wts, "[H3 ~> m3 or kg3 m-6]", 1)  ! Thickness cubed in a denominator
  call add_scaling(ns, des, wts, "[H2 T-2 ~> m2 s-2 or kg2 m-4 s-2]", 1) ! Thickness times f squared
  call add_scaling(ns, des, wts, "[H T2 R-1 Z-2 ~> m Pa-1 or s2 m-1]", 1) ! Pressure to thickness conversion
  call add_scaling(ns, des, wts, "[L2 Z-2 ~> nondim]", 1)     ! Inverse slope squared
  call add_scaling(ns, des, wts, "[H R L2 T-2 ~> m Pa]", 1)   ! Integral in thickness of pressure
  call add_scaling(ns, des, wts, "[R T2 Z-1 ~> kg s2 m-4]", 1) ! Density divided by gravitational acceleration

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
