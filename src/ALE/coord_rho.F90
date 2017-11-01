!> Regrid columns for the continuous isopycnal (rho) coordinate
module coord_rho

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_remapping,     only : remapping_CS, remapping_core_h
use MOM_EOS,           only : EOS_type, calculate_density
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid, DEGREE_MAX

implicit none ; private

!> Control structure containing required parameters for the rho coordinate
type, public :: rho_CS
  private

  !> Number of layers
  integer :: nk

  !> Minimum thickness allowed for layers
  real :: min_thickness = 0.

  !> Reference pressure for density calculations
  real :: ref_pressure

  !> If true, integrate for interface positions from the top downward.
  !! If false, integrate from the bottom upward, as does the rest of the model.
  logical :: integrate_downward_for_e = .false.

  !> Nominal density of interfaces
  real, allocatable, dimension(:) :: target_density

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS
end type rho_CS

!> Maximum number of regridding iterations
integer, parameter :: NB_REGRIDDING_ITERATIONS = 1
!> Deviation tolerance between succesive grids in regridding iterations
real, parameter    :: DEVIATION_TOLERANCE = 1e-10
! This CPP macro embeds some safety checks

public init_coord_rho, set_rho_params, build_rho_column, old_inflate_layers_1d, end_coord_rho

contains

!> Initialise a rho_CS with pointers to parameters
subroutine init_coord_rho(CS, nk, ref_pressure, target_density, interp_CS)
  type(rho_CS),         pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk
  real,                 intent(in) :: ref_pressure
  real, dimension(:),   intent(in) :: target_density
  type(interp_CS_type), intent(in) :: interp_CS

  if (associated(CS)) call MOM_error(FATAL, "init_coord_rho: CS already associated!")
  allocate(CS)
  allocate(CS%target_density(nk+1))

  CS%nk                = nk
  CS%ref_pressure      = ref_pressure
  CS%target_density(:) = target_density(:)
  CS%interp_CS         = interp_CS
end subroutine init_coord_rho

subroutine end_coord_rho(CS)
  type(rho_CS), pointer :: CS

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%target_density)
  deallocate(CS)
end subroutine end_coord_rho

subroutine set_rho_params(CS, min_thickness, integrate_downward_for_e, interp_CS)
  type(rho_CS),                   pointer    :: CS
  real,                 optional, intent(in) :: min_thickness
  logical,              optional, intent(in) :: integrate_downward_for_e
  type(interp_CS_type), optional, intent(in) :: interp_CS

  if (.not. associated(CS)) call MOM_error(FATAL, "set_rho_params: CS not associated")

  if (present(min_thickness)) CS%min_thickness = min_thickness
  if (present(integrate_downward_for_e)) CS%integrate_downward_for_e = integrate_downward_for_e
  if (present(interp_CS)) CS%interp_CS = interp_CS
end subroutine set_rho_params

!> Build a rho coordinate column
!!
!! 1. Density profiles are calculated on the source grid.
!! 2. Positions of target densities (for interfaces) are found by interpolation.
subroutine build_rho_column(CS, nz, depth, h, T, S, eqn_of_state, z_interface)
  type(rho_CS),             intent(in)    :: CS !< coord_rho control structure
  integer,                  intent(in)    :: nz !< Number of levels on source grid (i.e. length of  h, T, S)
  real,                     intent(in)    :: depth !< Depth of ocean bottom (positive in m)
  real, dimension(nz),      intent(in)    :: h  !< Layer thicknesses, in m
  real, dimension(nz),      intent(in)    :: T !< T for source column
  real, dimension(nz),      intent(in)    :: S !< S for source column
  type(EOS_type),           pointer       :: eqn_of_state !< Equation of state structure
  real, dimension(CS%nk+1), intent(inout) :: z_interface !< Absolute positions of interfaces
  ! Local variables
  integer :: k, count_nonzero_layers
  integer, dimension(nz) :: mapping
  real, dimension(nz) :: p, densities, h_nv
  real, dimension(nz+1) :: xTmp
  real, dimension(CS%nk) :: h_new ! New thicknesses
  real, dimension(CS%nk+1) :: x1

  ! Construct source column with vanished layers removed (stored in h_nv)
  call copy_finite_thicknesses(nz, h, CS%min_thickness, count_nonzero_layers, h_nv, mapping)

  if (count_nonzero_layers > 1) then
    xTmp(1) = 0.0
    do k = 1,count_nonzero_layers
      xTmp(k+1) = xTmp(k) + h_nv(k)
    end do

    ! Compute densities on source column
    p(:) = CS%ref_pressure
    call calculate_density(T, S, p, densities, 1, nz, eqn_of_state)
    do k = 1,count_nonzero_layers
      densities(k) = densities(mapping(k))
    end do

    ! Based on source column density profile, interpolate to generate a new grid
    call build_and_interpolate_grid(CS%interp_CS, densities, count_nonzero_layers, &
                                    h_nv, xTmp, CS%target_density, CS%nk, h_new, x1)

    ! Inflate vanished layers
    call old_inflate_layers_1d(CS%min_thickness, CS%nk, h_new)

    ! Comment: The following adjustment of h_new, and re-calculation of h_new via x1 needs to be removed
    x1(1) = 0.0 ; do k = 1,CS%nk ; x1(k+1) = x1(k) + h_new(k) ; end do
    do k = 1,CS%nk
      h_new(k) = x1(k+1) - x1(k)
    end do

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

subroutine build_rho_column_iteratively(CS, remapCS, nz, depth, h, T, S, eqn_of_state, zInterface)
  !< Iteratively uild a rho coordinate column
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

  type(rho_CS),          intent(in)    :: CS !< Regridding control structure
  type(remapping_CS),    intent(in)    :: remapCS !< Remapping parameters and options
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive in m)
  real, dimension(nz),   intent(in)    :: h  !< Layer thicknesses, in m
  real, dimension(nz),   intent(in)    :: T, S !< T and S for column
  type(EOS_type),        pointer       :: eqn_of_state !< Equation of state structure
  real, dimension(nz+1), intent(inout) :: zInterface !< Absolute positions of interfaces

  ! Local variables
  integer   :: k, m
  integer   :: count_nonzero_layers
  real      :: deviation            ! When iterating to determine the final
                                    ! grid, this is the deviation between two
                                    ! successive grids.
  real      :: threshold
  real, dimension(nz) :: p, densities, T_tmp, S_tmp, Tmp
  integer, dimension(nz) :: mapping
  real, dimension(nz) :: h0, h1, hTmp
  real, dimension(nz+1) :: x0, x1, xTmp

  threshold = CS%min_thickness
  p(:) = CS%ref_pressure
  T_tmp(:) = T(:)
  S_tmp(:) = S(:)
  h0(:) = h(:)

  ! Start iterations to build grid
  m = 1
  deviation = 1e10
  do while ( ( m <= NB_REGRIDDING_ITERATIONS ) .and. &
             ( deviation > DEVIATION_TOLERANCE ) )

    ! Construct column with vanished layers removed
    call copy_finite_thicknesses(nz, h0, threshold, count_nonzero_layers, hTmp, mapping)
    if ( count_nonzero_layers <= 1 ) then
      h1(:) = h0(:)
      exit  ! stop iterations here
    end if

    xTmp(1) = 0.0
    do k = 1,count_nonzero_layers
      xTmp(k+1) = xTmp(k) + hTmp(k)
    end do

    ! Compute densities within current water column
    call calculate_density( T_tmp, S_tmp, p, densities,&
                             1, nz, eqn_of_state )

    do k = 1,count_nonzero_layers
      densities(k) = densities(mapping(k))
    end do

    ! One regridding iteration
    ! Based on global density profile, interpolate to generate a new grid
    call build_and_interpolate_grid(CS%interp_CS, densities, count_nonzero_layers, &
         hTmp, xTmp, CS%target_density, nz, h1, x1)

    call old_inflate_layers_1d( CS%min_thickness, nz, h1 )
    x1(1) = 0.0 ; do k = 1,nz ; x1(k+1) = x1(k) + h1(k) ; end do

    ! Remap T and S from previous grid to new grid
    do k = 1,nz
      h1(k) = x1(k+1) - x1(k)
    end do

    call remapping_core_h(remapCS, nz, h0, S, nz, h1, Tmp)
    S_tmp(:) = Tmp(:)

    call remapping_core_h(remapCS, nz, h0, T, nz, h1, Tmp)
    T_tmp(:) = Tmp(:)

    ! Compute the deviation between two successive grids
    deviation = 0.0
    x0(1) = 0.0
    x1(1) = 0.0
    do k = 2,nz
      x0(k) = x0(k-1) + h0(k-1)
      x1(k) = x1(k-1) + h1(k-1)
      deviation = deviation + (x0(k)-x1(k))**2
    end do
    deviation = sqrt( deviation / (nz-1) )

    m = m + 1

    ! Copy final grid onto start grid for next iteration
    h0(:) = h1(:)

  end do ! end regridding iterations

  if (CS%integrate_downward_for_e) then
    zInterface(1) = 0.
    do k = 1,nz
      zInterface(k+1) = zInterface(k) - h1(k)
      ! Adjust interface position to accomodate inflating layers
      ! without disturbing the interface above
    enddo
  else
    ! The rest of the model defines grids integrating up from the bottom
    zInterface(nz+1) = -depth
    do k = nz,1,-1
      zInterface(k) = zInterface(k+1) + h1(k)
      ! Adjust interface position to accomodate inflating layers
      ! without disturbing the interface above
    enddo
  endif

end subroutine build_rho_column_iteratively

!> Copy column thicknesses with vanished layers removed
subroutine copy_finite_thicknesses(nk, h_in, threshold, nout, h_out, mapping)
  integer,                intent(in)  :: nk !< Number of layer for h_in, T_in, S_in
  real, dimension(nk),    intent(in)  :: h_in !< Thickness of input column
  real,                   intent(in)  :: threshold !< Thickness threshold defining vanished layers
  integer,                intent(out) :: nout !< Number of non-vanished layers
  real, dimension(nk),    intent(out) :: h_out !< Thickness of output column
  integer, dimension(nk), intent(out) :: mapping !< Index of k-out corresponding to k-in
  ! Local variables
  integer :: k, k_thickest
  real :: thickness_in_vanished, thickest_h_out

  ! Build up new grid
  nout = 0
  thickness_in_vanished = 0.0
  thickest_h_out = h_in(1)
  k_thickest = 1
  do k = 1, nk
    mapping(k) = nout ! Note k>=nout always
    h_out(k) = 0.  ! Make sure h_out is set everywhere
    if (h_in(k) > threshold) then
      ! For non-vanished layers
      nout = nout + 1
      mapping(nout) = k
      h_out(nout) = h_in(k)
      if (h_out(nout) > thickest_h_out) then
        thickest_h_out = h_out(nout)
        k_thickest = nout
      end if
    else
      ! Add up mass in vanished layers
      thickness_in_vanished = thickness_in_vanished + h_in(k)
    end if
  end do

  ! No finite layers
  if (nout <= 1) return

  ! Adjust for any lost volume in vanished layers
  h_out(k_thickest) = h_out(k_thickest) + thickness_in_vanished

end subroutine copy_finite_thicknesses

!------------------------------------------------------------------------------
! Inflate vanished layers to finite (nonzero) width
!------------------------------------------------------------------------------
subroutine old_inflate_layers_1d( minThickness, N, h )

  ! Argument
  real,                intent(in) :: minThickness
  integer,             intent(in) :: N
  real,                intent(inout) :: h(:)

  ! Local variable
  integer   :: k
  integer   :: k_found
  integer   :: count_nonzero_layers
  real      :: delta
  real      :: correction
  real      :: maxThickness

  ! Count number of nonzero layers
  count_nonzero_layers = 0
  do k = 1,N
    if ( h(k) > minThickness ) then
      count_nonzero_layers = count_nonzero_layers + 1
    end if
  end do

  ! If all layer thicknesses are greater than the threshold, exit routine
  if ( count_nonzero_layers == N ) return

  ! If all thicknesses are zero, inflate them all and exit
  if ( count_nonzero_layers == 0 ) then
    do k = 1,N
      h(k) = minThickness
    end do
    return
  end if

  ! Inflate zero layers
  correction = 0.0
  do k = 1,N
    if ( h(k) <= minThickness ) then
      delta = minThickness - h(k)
      correction = correction + delta
      h(k) = h(k) + delta
    end if
  end do

  ! Modify thicknesses of nonzero layers to ensure volume conservation
  maxThickness = h(1)
  k_found = 1
  do k = 1,N
    if ( h(k) > maxThickness ) then
      maxThickness = h(k)
      k_found = k
    end if
  end do

  h(k_found) = h(k_found) - correction

end subroutine old_inflate_layers_1d

end module coord_rho
