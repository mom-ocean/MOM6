!> Regrid columns for the continuous isopycnal (rho) coordinate
module coord_rho

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
  real :: min_thickness

  !> Reference pressure for density calculations
  real :: ref_pressure

  !> If true, integrate for interface positions from the top downward.
  !! If false, integrate from the bottom upward, as does the rest of the model.
  logical :: integrate_downward_for_e

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

subroutine build_rho_column(CS, remapCS, nz, depth, h, T, S, eqn_of_state, zInterface)
  !< Build a rho coordinate column
  !!
  !! The algorithn operates as follows within each column:
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
  integer   :: map_index
  integer   :: k_found
  integer   :: count_nonzero_layers
  real      :: deviation            ! When iterating to determine the final
                                    ! grid, this is the deviation between two
                                    ! successive grids.
  real      :: threshold
  real      :: max_thickness
  real      :: correction
  real, dimension(nz) :: p, densities, T_tmp, S_tmp, Tmp
  integer, dimension(nz) :: mapping
  real :: dh
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

    ! Count number of nonzero layers within current water column
    count_nonzero_layers = 0
    do k = 1,nz
      if ( h0(k) > threshold ) then
        count_nonzero_layers = count_nonzero_layers + 1
      end if
    end do

    ! If there is at most one nonzero layer, stop here (no regridding)
    if ( count_nonzero_layers <= 1 ) then
      h1(:) = h0(:)
      exit  ! stop iterations here
    end if

    ! Build new grid containing only nonzero layers
    map_index = 1
    correction = 0.0
    do k = 1,nz
      if ( h0(k) > threshold ) then
        mapping(map_index) = k
        hTmp(map_index) = h0(k)
        map_index = map_index + 1
      else
        correction = correction + h0(k)
      end if
    end do

    max_thickness = hTmp(1)
    k_found = 1
    do k = 1,count_nonzero_layers
      if ( hTmp(k) > max_thickness ) then
        max_thickness = hTmp(k)
        k_found = k
      end if
    end do

    hTmp(k_found) = hTmp(k_found) + correction

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

end subroutine build_rho_column

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
