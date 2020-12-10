!> Regrid columns for the continuous isopycnal (rho) coordinate
module coord_rho

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_remapping,     only : remapping_CS, remapping_core_h
use MOM_EOS,           only : EOS_type, calculate_density
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid, DEGREE_MAX

implicit none ; private

!> Control structure containing required parameters for the rho coordinate
type, public :: rho_CS ; private

  !> Number of layers
  integer :: nk

  !> Minimum thickness allowed for layers, often in [H ~> m or kg m-2]
  real :: min_thickness = 0.

  !> Reference pressure for density calculations [R L2 T-2 ~> Pa]
  real :: ref_pressure

  !> If true, integrate for interface positions from the top downward.
  !! If false, integrate from the bottom upward, as does the rest of the model.
  logical :: integrate_downward_for_e = .false.

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS
end type rho_CS

public init_coord_rho, set_rho_params, build_rho_column, old_inflate_layers_1d, end_coord_rho

contains

!> Initialise a rho_CS with pointers to parameters
subroutine init_coord_rho(CS, nk, ref_pressure, target_density, interp_CS)
  type(rho_CS),         pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk !< Number of layers in the grid
  real,                 intent(in) :: ref_pressure !< Coordinate reference pressure [R L2 T-2 ~> Pa]
  real, dimension(:),   intent(in) :: target_density !< Nominal density of interfaces [R ~> kg m-3]
  type(interp_CS_type), intent(in) :: interp_CS !< Controls for interpolation

  if (associated(CS)) call MOM_error(FATAL, "init_coord_rho: CS already associated!")
  allocate(CS)
  allocate(CS%target_density(nk+1))

  CS%nk                = nk
  CS%ref_pressure      = ref_pressure
  CS%target_density(:) = target_density(:)
  CS%interp_CS         = interp_CS

end subroutine init_coord_rho

!> This subroutine deallocates memory in the control structure for the coord_rho module
subroutine end_coord_rho(CS)
  type(rho_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%target_density)
  deallocate(CS)
end subroutine end_coord_rho

!> This subroutine can be used to set the parameters for the coord_rho module
subroutine set_rho_params(CS, min_thickness, integrate_downward_for_e, interp_CS)
  type(rho_CS),      pointer    :: CS !< Coordinate control structure
  real,    optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]
  logical, optional, intent(in) :: integrate_downward_for_e !< If true, integrate for interface
                                      !! positions from the top downward.  If false, integrate
                                      !! from the bottom upward, as does the rest of the model.

  type(interp_CS_type), optional, intent(in) :: interp_CS !< Controls for interpolation

  if (.not. associated(CS)) call MOM_error(FATAL, "set_rho_params: CS not associated")

  if (present(min_thickness)) CS%min_thickness = min_thickness
  if (present(integrate_downward_for_e)) CS%integrate_downward_for_e = integrate_downward_for_e
  if (present(interp_CS)) CS%interp_CS = interp_CS
end subroutine set_rho_params

!> Build a rho coordinate column
!!
!! 1. Density profiles are calculated on the source grid.
!! 2. Positions of target densities (for interfaces) are found by interpolation.
subroutine build_rho_column(CS, nz, depth, h, T, S, eqn_of_state, z_interface, &
                            z_rigid_top, eta_orig, h_neglect, h_neglect_edge)
  type(rho_CS),        intent(in)    :: CS !< coord_rho control structure
  integer,             intent(in)    :: nz !< Number of levels on source grid (i.e. length of  h, T, S)
  real,                intent(in)    :: depth !< Depth of ocean bottom (positive downward) [H ~> m or kg m-2]
  real, dimension(nz), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz), intent(in)    :: T  !< Temperature for source column [degC]
  real, dimension(nz), intent(in)    :: S  !< Salinity for source column [ppt]
  type(EOS_type),      pointer       :: eqn_of_state !< Equation of state structure
  real, dimension(CS%nk+1), &
       intent(inout) :: z_interface !< Absolute positions of interfaces
  real, optional,           intent(in)    :: z_rigid_top !< The height of a rigid top (positive upward in the same
  !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real, optional,           intent(in)    :: eta_orig !< The actual original height of the top in the same
                                                   !! units as depth) [Z ~> m] or [H ~> m or kg m-2]
  real,      optional, intent(in)    :: h_neglect !< A negligibly small width for the purpose
                                             !! of cell reconstructions [H ~> m or kg m-2]
  real,      optional, intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose
                                             !! of edge value calculations [H ~> m or kg m-2]

  ! Local variables
  integer :: k, count_nonzero_layers
  integer, dimension(nz) :: mapping
  real, dimension(nz) :: pres     ! Pressures used to calculate density [R L2 T-2 ~> Pa]
  real, dimension(nz) :: h_nv     ! Thicknesses of non-vanishing layers [H ~> m or kg m-2]
  real, dimension(nz) :: densities ! Layer density [R ~> kg m-3]
  real, dimension(nz+1) :: xTmp   ! Temporary positions [H ~> m or kg m-2]
  real, dimension(CS%nk) :: h_new ! New thicknesses [H ~> m or kg m-2]
  real, dimension(CS%nk+1) :: x1  ! Interface heights [H ~> m or kg m-2]
  real :: z0_top, eta ! Thicknesses or heights [Z ~> m] or [H ~> m or kg m-2]

  ! Construct source column with vanished layers removed (stored in h_nv)
  call copy_finite_thicknesses(nz, h, CS%min_thickness, count_nonzero_layers, h_nv, mapping)

  z0_top = 0.
  eta=0.0
  if (present(z_rigid_top)) then
    z0_top = z_rigid_top
    eta=z0_top
    if (present(eta_orig)) then
       eta=eta_orig
    endif
  endif


  if (count_nonzero_layers > 1) then
    xTmp(1) = 0.0
    do k = 1,count_nonzero_layers
      xTmp(k+1) = xTmp(k) + h_nv(k)
    enddo

    ! Compute densities on source column
    pres(:) = CS%ref_pressure
    call calculate_density(T, S, pres, densities, eqn_of_state)
    do k = 1,count_nonzero_layers
      densities(k) = densities(mapping(k))
    enddo

    ! Based on source column density profile, interpolate to generate a new grid
    call build_and_interpolate_grid(CS%interp_CS, densities, count_nonzero_layers, &
                                    h_nv, xTmp, CS%target_density, CS%nk, h_new, &
                                    x1, h_neglect, h_neglect_edge)

    ! Inflate vanished layers
    call old_inflate_layers_1d(CS%min_thickness, CS%nk, h_new)

    ! Comment: The following adjustment of h_new, and re-calculation of h_new via x1 needs to be removed
    x1(1) = 0.0 ; do k = 1,CS%nk ; x1(k+1) = x1(k) + h_new(k) ; enddo
    do k = 1,CS%nk
      h_new(k) = x1(k+1) - x1(k)
    enddo

  else ! count_nonzero_layers <= 1
    if (nz == CS%nk) then
      h_new(:) = h(:) ! This keeps old behavior
    else
      h_new(:) = 0.
      h_new(1) = h(1)
    endif
  endif

  ! Return interface positions
  if (CS%integrate_downward_for_e) then
    ! Remapping is defined integrating from zero
    z_interface(1) = 0.
    do k = 1,CS%nk
      z_interface(k+1) = z_interface(k) - h_new(k)
    enddo
  else
    ! The rest of the model defines grids integrating up from the bottom
    z_interface(CS%nk+1) = -depth
    do k = CS%nk,1,-1
      z_interface(k) = z_interface(k+1) + h_new(k)
    enddo
  endif

end subroutine build_rho_column

!### build_rho_column_iteratively is never used or called.

!> Iteratively build a rho coordinate column
!!
!! The algorithm operates as follows within each column:
!!
!! 1. Given T & S within each layer, the layer densities are computed.
!! 2. Based on these layer densities, a global density profile is reconstructed
!!    (this profile is monotonically increasing and may be discontinuous)
!! 3. The new grid interfaces are determined based on the target interface
!!    densities.
!! 4. T & S are remapped onto the new grid.
!! 5. Return to step 1 until convergence or until the maximum number of
!!    iterations is reached, whichever comes first.
subroutine build_rho_column_iteratively(CS, remapCS, nz, depth, h, T, S, eqn_of_state, &
                                        zInterface, h_neglect, h_neglect_edge, dev_tol)
  type(rho_CS),          intent(in)    :: CS !< Regridding control structure
  type(remapping_CS),    intent(in)    :: remapCS !< Remapping parameters and options
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom [Z ~> m]
  real, dimension(nz),   intent(in)    :: h  !< Layer thicknesses in Z coordinates [Z ~> m]
  real, dimension(nz),   intent(in)    :: T  !< T for column [degC]
  real, dimension(nz),   intent(in)    :: S  !< S for column [ppt]
  type(EOS_type),        pointer       :: eqn_of_state !< Equation of state structure
  real, dimension(nz+1), intent(inout) :: zInterface !< Absolute positions of interfaces
  real,        optional, intent(in)    :: h_neglect !< A negligibly small width for the
                                             !! purpose of cell reconstructions
                                             !! in the same units as h [Z ~> m]
  real,        optional, intent(in)    :: h_neglect_edge !< A negligibly small width
                                             !! for the purpose of edge value calculations
                                             !! in the same units as h [Z ~> m]
  real,        optional, intent(in)    :: dev_tol !< The tolerance for the deviation between
                                             !! successive grids for determining when the
                                             !! iterative solver has converged [Z ~> m]

  ! Local variables
  real, dimension(nz+1) :: x0, x1, xTmp ! Temporary interface heights [Z ~> m]
  real, dimension(nz) :: pres       ! The pressure used in the equation of state [R L2 T-2 ~> Pa].
  real, dimension(nz) :: densities  ! Layer densities [R ~> kg m-3]
  real, dimension(nz) :: T_tmp, S_tmp ! A temporary profile of temperature [degC] and salinity [ppt].
  real, dimension(nz) :: Tmp        ! A temporary variable holding a remapped variable.
  real, dimension(nz) :: h0, h1, hTmp ! Temporary thicknesses [Z ~> m]
  real :: deviation            ! When iterating to determine the final grid, this is the
                               ! deviation between two successive grids [Z ~> m].
  real :: deviation_tol        ! Deviation tolerance between succesive grids in
                               ! regridding iterations [Z ~> m]
  real :: threshold            ! The minimum thickness for a layer to be considered to exist [Z ~> m]
  integer, dimension(nz) :: mapping ! The indices of the massive layers in the initial column.
  integer :: k, m, count_nonzero_layers

  !  Maximum number of regridding iterations
  integer, parameter :: NB_REGRIDDING_ITERATIONS = 1

  threshold = CS%min_thickness
  pres(:) = CS%ref_pressure
  T_tmp(:) = T(:)
  S_tmp(:) = S(:)
  h0(:) = h(:)

  ! Start iterations to build grid
  m = 1
  deviation_tol = 1.0e-15*depth ; if (present(dev_tol)) deviation_tol = dev_tol

  do m=1,NB_REGRIDDING_ITERATIONS

    ! Construct column with vanished layers removed
    call copy_finite_thicknesses(nz, h0, threshold, count_nonzero_layers, hTmp, mapping)
    if ( count_nonzero_layers <= 1 ) then
      h1(:) = h0(:)
      exit  ! stop iterations here
    endif

    xTmp(1) = 0.0
    do k = 1,count_nonzero_layers
      xTmp(k+1) = xTmp(k) + hTmp(k)
    enddo

    ! Compute densities within current water column
    call calculate_density( T_tmp, S_tmp, pres, densities, eqn_of_state)

    do k = 1,count_nonzero_layers
      densities(k) = densities(mapping(k))
    enddo

    ! One regridding iteration
    ! Based on global density profile, interpolate to generate a new grid
    call build_and_interpolate_grid(CS%interp_CS, densities, count_nonzero_layers, &
         hTmp, xTmp, CS%target_density, nz, h1, x1, h_neglect, h_neglect_edge)

    call old_inflate_layers_1d( CS%min_thickness, nz, h1 )
    x1(1) = 0.0 ; do k = 1,nz ; x1(k+1) = x1(k) + h1(k) ; enddo

    ! Remap T and S from previous grid to new grid
    do k = 1,nz
      h1(k) = x1(k+1) - x1(k)
    enddo

    call remapping_core_h(remapCS, nz, h0, S, nz, h1, Tmp, h_neglect, h_neglect_edge)
    S_tmp(:) = Tmp(:)

    call remapping_core_h(remapCS, nz, h0, T, nz, h1, Tmp, h_neglect, h_neglect_edge)
    T_tmp(:) = Tmp(:)

    ! Compute the deviation between two successive grids
    deviation = 0.0
    x0(1) = 0.0
    x1(1) = 0.0
    do k = 2,nz
      x0(k) = x0(k-1) + h0(k-1)
      x1(k) = x1(k-1) + h1(k-1)
      deviation = deviation + (x0(k)-x1(k))**2
    enddo
    deviation = sqrt( deviation / (nz-1) )

    if ( deviation <= deviation_tol ) exit

    ! Copy final grid onto start grid for next iteration
    h0(:) = h1(:)
  enddo ! end regridding iterations

  if (CS%integrate_downward_for_e) then
    zInterface(1) = 0.
    do k = 1,nz
      zInterface(k+1) = zInterface(k) - h1(k)
      ! Adjust interface position to accommodate inflating layers
      ! without disturbing the interface above
    enddo
  else
    ! The rest of the model defines grids integrating up from the bottom
    zInterface(nz+1) = -depth
    do k = nz,1,-1
      zInterface(k) = zInterface(k+1) + h1(k)
      ! Adjust interface position to accommodate inflating layers
      ! without disturbing the interface above
    enddo
  endif

end subroutine build_rho_column_iteratively

!> Copy column thicknesses with vanished layers removed
subroutine copy_finite_thicknesses(nk, h_in, thresh, nout, h_out, mapping)
  integer,                intent(in)  :: nk      !< Number of layer for h_in, T_in, S_in
  real, dimension(nk),    intent(in)  :: h_in    !< Thickness of input column [H ~> m or kg m-2] or [Z ~> m]
  real,                   intent(in)  :: thresh  !< Thickness threshold defining vanished
                                                 !! layers [H ~> m or kg m-2] or [Z ~> m]
  integer,                intent(out) :: nout    !< Number of non-vanished layers
  real, dimension(nk),    intent(out) :: h_out   !< Thickness of output column [H ~> m or kg m-2] or [Z ~> m]
  integer, dimension(nk), intent(out) :: mapping !< Index of k-out corresponding to k-in
  ! Local variables
  integer :: k, k_thickest
  real :: thickness_in_vanished ! Summed thicknesses in discarded layers [H ~> m or kg m-2] or [Z ~> m]
  real :: thickest_h_out        ! Thickness of the thickest layer [H ~> m or kg m-2] or [Z ~> m]

  ! Build up new grid
  nout = 0
  thickness_in_vanished = 0.0
  thickest_h_out = h_in(1)
  k_thickest = 1
  do k = 1, nk
    mapping(k) = nout ! Note k>=nout always
    h_out(k) = 0.  ! Make sure h_out is set everywhere
    if (h_in(k) > thresh) then
      ! For non-vanished layers
      nout = nout + 1
      mapping(nout) = k
      h_out(nout) = h_in(k)
      if (h_out(nout) > thickest_h_out) then
        thickest_h_out = h_out(nout)
        k_thickest = nout
      endif
    else
      ! Add up mass in vanished layers
      thickness_in_vanished = thickness_in_vanished + h_in(k)
    endif
  enddo

  ! No finite layers
  if (nout <= 1) return

  ! Adjust for any lost volume in vanished layers
  h_out(k_thickest) = h_out(k_thickest) + thickness_in_vanished

end subroutine copy_finite_thicknesses

!------------------------------------------------------------------------------
!> Inflate vanished layers to finite (nonzero) width
subroutine old_inflate_layers_1d( min_thickness, nk, h )

  ! Argument
  real,               intent(in)    :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]
  integer,            intent(in)    :: nk  !< Number of layers in the grid
  real, dimension(:), intent(inout) :: h   !< Layer thicknesses [H ~> m or kg m-2]

  ! Local variable
  integer   :: k
  integer   :: k_found
  integer   :: count_nonzero_layers
  real      :: delta
  real      :: correction
  real      :: maxThickness

  ! Count number of nonzero layers
  count_nonzero_layers = 0
  do k = 1,nk
    if ( h(k) > min_thickness ) then
      count_nonzero_layers = count_nonzero_layers + 1
    endif
  enddo

  ! If all layer thicknesses are greater than the threshold, exit routine
  if ( count_nonzero_layers == nk ) return

  ! If all thicknesses are zero, inflate them all and exit
  if ( count_nonzero_layers == 0 ) then
    do k = 1,nk
      h(k) = min_thickness
    enddo
    return
  endif

  ! Inflate zero layers
  correction = 0.0
  do k = 1,nk
    if ( h(k) <= min_thickness ) then
      delta = min_thickness - h(k)
      correction = correction + delta
      h(k) = h(k) + delta
    endif
  enddo

  ! Modify thicknesses of nonzero layers to ensure volume conservation
  maxThickness = h(1)
  k_found = 1
  do k = 1,nk
    if ( h(k) > maxThickness ) then
      maxThickness = h(k)
      k_found = k
    endif
  enddo

  h(k_found) = h(k_found) - correction

end subroutine old_inflate_layers_1d

end module coord_rho
