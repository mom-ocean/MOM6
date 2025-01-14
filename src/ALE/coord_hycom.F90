!> Regrid columns for the HyCOM coordinate
module coord_hycom

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_remapping,     only : remapping_CS, remapping_core_h
use MOM_EOS,           only : EOS_type, calculate_density
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid, regridding_set_ppolys
use regrid_interp,     only : DEGREE_MAX

implicit none ; private

!> Control structure containing required parameters for the HyCOM coordinate
type, public :: hycom_CS ; private

  !> Number of layers/levels in generated grid
  integer :: nk

  !> Nominal near-surface resolution [Z ~> m]
  real, allocatable, dimension(:) :: coordinateResolution

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  !> Maximum depths of interfaces [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: max_interface_depths

  !> Maximum thicknesses of layers [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: max_layer_thickness

  !> If true, an interface only moves if it improves the density fit
  logical :: only_improves = .false.

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS
end type hycom_CS

public init_coord_hycom, set_hycom_params, build_hycom1_column, end_coord_hycom

contains

!> Initialise a hycom_CS with pointers to parameters
subroutine init_coord_hycom(CS, nk, coordinateResolution, target_density, interp_CS)
  type(hycom_CS),       pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk !< Number of layers in generated grid
  real, dimension(nk),  intent(in) :: coordinateResolution !< Nominal near-surface resolution [Z ~> m]
  real, dimension(nk+1),intent(in) :: target_density !< Interface target densities [R ~> kg m-3]
  type(interp_CS_type), intent(in) :: interp_CS !< Controls for interpolation

  if (associated(CS)) call MOM_error(FATAL, "init_coord_hycom: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk))
  allocate(CS%target_density(nk+1))

  CS%nk                      = nk
  CS%coordinateResolution(:) = coordinateResolution(:)
  CS%target_density(:)       = target_density(:)
  CS%interp_CS               = interp_CS

end subroutine init_coord_hycom

!> This subroutine deallocates memory in the control structure for the coord_hycom module
subroutine end_coord_hycom(CS)
  type(hycom_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS%target_density)
  if (allocated(CS%max_interface_depths)) deallocate(CS%max_interface_depths)
  if (allocated(CS%max_layer_thickness)) deallocate(CS%max_layer_thickness)
  deallocate(CS)
end subroutine end_coord_hycom

!> This subroutine can be used to set the parameters for the coord_hycom module
subroutine set_hycom_params(CS, max_interface_depths, max_layer_thickness, only_improves, interp_CS)
  type(hycom_CS),                 pointer    :: CS !< Coordinate control structure
  real, dimension(:),   optional, intent(in) :: max_interface_depths !< Maximum depths of interfaces [H ~> m or kg m-2]
  real, dimension(:),   optional, intent(in) :: max_layer_thickness  !< Maximum thicknesses of layers [H ~> m or kg m-2]
  logical, optional, intent(in) :: only_improves !< If true, an interface only moves if it improves the density fit
  type(interp_CS_type), optional, intent(in) :: interp_CS !< Controls for interpolation

  if (.not. associated(CS)) call MOM_error(FATAL, "set_hycom_params: CS not associated")

  if (present(max_interface_depths)) then
    if (size(max_interface_depths) /= CS%nk+1) &
      call MOM_error(FATAL, "set_hycom_params: max_interface_depths inconsistent size")
    allocate(CS%max_interface_depths(CS%nk+1))
    CS%max_interface_depths(:) = max_interface_depths(:)
  endif

  if (present(max_layer_thickness)) then
    if (size(max_layer_thickness) /= CS%nk) &
      call MOM_error(FATAL, "set_hycom_params: max_layer_thickness inconsistent size")
    allocate(CS%max_layer_thickness(CS%nk))
    CS%max_layer_thickness(:) = max_layer_thickness(:)
  endif

  if (present(only_improves)) CS%only_improves = only_improves

  if (present(interp_CS)) CS%interp_CS = interp_CS
end subroutine set_hycom_params

!> Build a HyCOM coordinate column
subroutine build_hycom1_column(CS, remapCS, eqn_of_state, nz, depth, h, T, S, p_col, &
                               z_col, z_col_new, zScale, h_neglect, h_neglect_edge)
  type(hycom_CS),        intent(in)    :: CS    !< Coordinate control structure
  type(remapping_CS),    intent(in)    :: remapCS !< Remapping parameters and options
  type(EOS_type),        intent(in)    :: eqn_of_state !< Equation of state structure
  integer,               intent(in)    :: nz    !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive [H ~> m or kg m-2])
  real, dimension(nz),   intent(in)    :: T     !< Temperature of column [C ~> degC]
  real, dimension(nz),   intent(in)    :: S     !< Salinity of column [S ~> ppt]
  real, dimension(nz),   intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz),   intent(in)    :: p_col !< Layer pressure [R L2 T-2 ~> Pa]
  real, dimension(nz+1), intent(in)    :: z_col !< Interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: z_col_new !< Absolute positions of interfaces [H ~> m or kg m-2]
  real, optional,        intent(in)    :: zScale !< Scaling factor from the input coordinate thicknesses in [Z ~> m]
                                                !! to desired units for zInterface, perhaps GV%Z_to_H in which
                                                !! case this has units of [H Z-1 ~> nondim or kg m-3]
  real,                  intent(in)    :: h_neglect !< A negligibly small width for the purpose of
                                                !! cell reconstruction [H ~> m or kg m-2]
  real,        optional, intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose of
                                                !! edge value calculation [H ~> m or kg m-2]

  ! Local variables
  integer   :: k
  real, dimension(nz)      :: rho_col   ! Layer densities in a column [R ~> kg m-3]
  real, dimension(CS%nk)   :: h_col_new ! New layer thicknesses [H ~> m or kg m-2]
  real, dimension(CS%nk)   :: r_col_new ! New layer densities [R ~> kg m-3]
  real, dimension(CS%nk)   :: T_col_new ! New layer temperatures [C ~> degC]
  real, dimension(CS%nk)   :: S_col_new ! New layer salinities [S ~> ppt]
  real, dimension(CS%nk)   :: p_col_new ! New layer pressure [R L2 T-2 ~> Pa]
  real, dimension(CS%nk+1) :: RiA_ini   ! Initial nk+1 interface density anomaly w.r.t. the
                                        ! interface target densities [R ~> kg m-3]
  real, dimension(CS%nk+1) :: RiA_new   ! New interface density anomaly w.r.t. the
                                        ! interface target densities [R ~> kg m-3]
  real :: z_1, z_nz  ! mid point of 1st and last layers [H ~> m or kg m-2]
  real :: z_scale    ! A scaling factor from the input thicknesses to the target thicknesses,
                     ! perhaps 1 or a factor in [H Z-1 ~> 1 or kg m-3]
  real :: stretching ! z* stretching, converts z* to z [nondim].
  real :: nominal_z ! Nominal depth of interface when using z* [H ~> m or kg m-2]
  logical :: maximum_depths_set ! If true, the maximum depths of interface have been set.
  logical :: maximum_h_set      ! If true, the maximum layer thicknesses have been set.

  maximum_depths_set = allocated(CS%max_interface_depths)
  maximum_h_set = allocated(CS%max_layer_thickness)

  z_scale = 1.0 ; if (present(zScale)) z_scale = zScale

  if (CS%only_improves .and. nz == CS%nk) then
    call build_hycom1_target_anomaly(CS, remapCS, eqn_of_state, CS%nk, depth, &
        h, T, S, p_col, rho_col, RiA_ini, h_neglect, h_neglect_edge)
  else
    ! Work bottom recording potential density
    call calculate_density(T, S, p_col, rho_col, eqn_of_state)
    ! This ensures the potential density profile is monotonic
    ! although not necessarily single valued.
    do k = nz-1, 1, -1
      rho_col(k) = min( rho_col(k), rho_col(k+1) )
    enddo
  endif

  ! Interpolates for the target interface position with the rho_col profile
  ! Based on global density profile, interpolate to generate a new grid
  call build_and_interpolate_grid(CS%interp_CS, rho_col, nz, h(:), z_col, &
           CS%target_density, CS%nk, h_col_new, z_col_new, h_neglect, h_neglect_edge)
  if (CS%only_improves .and. nz == CS%nk) then
    ! Only move an interface if it improves the density fit
    z_1 = 0.5 * ( z_col(1) + z_col(2) )
    z_nz  = 0.5 * ( z_col(nz) + z_col(nz+1) )
    do k = 1,CS%nk
      p_col_new(k) = p_col(1) + ( 0.5 * ( z_col_new(K) + z_col_new(K+1) ) - z_1 ) / ( z_nz - z_1 ) * &
          ( p_col(nz) - p_col(1) )
    enddo
    ! Remap from original h and T,S to get T,S_col_new
    call remapping_core_h(remapCS, nz, h(:), T, CS%nk, h_col_new, T_col_new)
    call remapping_core_h(remapCS, nz, h(:), S, CS%nk, h_col_new, S_col_new)
    call build_hycom1_target_anomaly(CS, remapCS, eqn_of_state, CS%nk, depth, &
        h_col_new, T_col_new, S_col_new, p_col_new, r_col_new, RiA_new, h_neglect, h_neglect_edge)
    do k= 2,CS%nk
      if     ( abs(RiA_ini(K)) <= abs(RiA_new(K)) .and. z_col(K) > z_col_new(K-1) .and. &
               z_col(K) < z_col_new(K+1)) then
        z_col_new(K) = z_col(K)
      endif
    enddo
  endif !only_improves

  ! Sweep down the interfaces and make sure that the interface is at least
  ! as deep as a nominal target z* grid
  nominal_z = 0.
  stretching = z_col(nz+1) / depth ! Stretches z* to z
  do k = 2, CS%nk+1
    nominal_z = nominal_z + (z_scale * CS%coordinateResolution(k-1)) * stretching
    z_col_new(k) = max( z_col_new(k), nominal_z )
    z_col_new(k) = min( z_col_new(k), z_col(nz+1) )
  enddo

  if (maximum_depths_set .and. maximum_h_set) then ; do k=2,CS%nk
    ! The loop bounds are 2 & nz so the top and bottom interfaces do not move.
    ! Recall that z_col_new is positive downward.
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K), &
                       z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; elseif (maximum_depths_set) then ; do K=2,CS%nk
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K))
  enddo ; elseif (maximum_h_set) then ; do k=2,CS%nk
    z_col_new(K) = min(z_col_new(K), z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; endif
end subroutine build_hycom1_column

!> Calculate interface density anomaly w.r.t. the target.
subroutine build_hycom1_target_anomaly(CS, remapCS, eqn_of_state, nz, depth, h, T, S, p_col, &
                                       R, RiAnom, h_neglect, h_neglect_edge)
  type(hycom_CS),        intent(in)  :: CS     !< Coordinate control structure
  type(remapping_CS),    intent(in)  :: remapCS !< Remapping parameters and options
  type(EOS_type),        intent(in)  :: eqn_of_state !< Equation of state structure
  integer,               intent(in)  :: nz     !< Number of levels
  real,                  intent(in)  :: depth  !< Depth of ocean bottom (positive [H ~> m or kg m-2])
  real, dimension(nz),   intent(in)  :: T      !< Temperature of column [C ~> degC]
  real, dimension(nz),   intent(in)  :: S      !< Salinity of column [S ~> ppt]
  real, dimension(nz),   intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz),   intent(in)  :: p_col  !< Layer pressure [R L2 T-2 ~> Pa]
  real, dimension(nz),   intent(out) :: R      !< Layer density [R ~> kg m-3]
  real, dimension(nz+1), intent(out) :: RiAnom !< The interface density anomaly
                                               !! w.r.t. the interface target
                                               !! densities [R ~> kg m-3]
  real,                  intent(in)  :: h_neglect !< A negligibly small width for the purpose of
                                               !! cell reconstruction [H ~> m or kg m-2]
  real,        optional, intent(in)  :: h_neglect_edge !< A negligibly small width for the purpose of
                                                !! edge value calculation [H ~> m or kg m-2]
  ! Local variables
  integer   :: degree,k
  real, dimension(nz)   :: rho_col ! Layer densities in a column [R ~> kg m-3]
  real, dimension(nz,2) :: ppoly_E ! Polynomial edge values [R ~> kg m-3]
  real, dimension(nz,2) :: ppoly_S ! Polynomial edge slopes [R H-1]
  real, dimension(nz,DEGREE_MAX+1) :: ppoly_C ! Polynomial interpolant coeficients on the local 0-1 grid [R ~> kg m-3]

  ! Work bottom recording potential density
  call calculate_density(T, S, p_col, rho_col, eqn_of_state)
  ! This ensures the potential density profile is monotonic
  ! although not necessarily single valued.
  do k = nz-1, 1, -1
    rho_col(k) = min( rho_col(k), rho_col(k+1) )
  enddo

  call regridding_set_ppolys(CS%interp_CS, rho_col, nz, h, ppoly_E, ppoly_S, ppoly_C, &
                             degree, h_neglect, h_neglect_edge)

  R(1) = rho_col(1)
  RiAnom(1) = ppoly_E(1,1) - CS%target_density(1)
  do k= 2,nz
    R(k) = rho_col(k)
    if (ppoly_E(k-1,2) > CS%target_density(k)) then
      RiAnom(k) = ppoly_E(k-1,2) - CS%target_density(k)  !interface is heavier than target
    elseif (ppoly_E(k,1) < CS%target_density(k)) then
      RiAnom(k) = ppoly_E(k,1)   - CS%target_density(k)  !interface is lighter than target
    else
      RiAnom(k) = 0.0  !interface spans the target
    endif
  enddo
  RiAnom(nz+1) = ppoly_E(nz,2) - CS%target_density(nz+1)

end subroutine build_hycom1_target_anomaly

end module coord_hycom
