! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A generic type for vertical 1D reconstructions
module Recon1d_type

use numerical_testing_type, only : testing

implicit none ; private

public Recon1d
public testing

!> The base class for implementations of 1D reconstructions
type, abstract :: Recon1d

  integer :: n = 0 !< Number of cells in column
  real, allocatable, dimension(:) :: u_mean !< Cell mean [A]
  real :: h_neglect = 0. !< A negligibly small width used in cell reconstructions in the same units as h [H]
  real :: x_tolerance = 1. * epsilon(1.) !< Solver tolerance for x in element (0,1) [nondim]
  logical :: check = .false. !< If true, enable some consistency checking

  logical :: debug = .false. !< If true, dump info as calculations are made (do not enable)
contains

  ! The following functions/subroutines are deferred and must be provided specifically by each scheme

  !> Deferred implementation of initialization
  procedure(i_init), deferred :: init
  !> Deferred implementation of reconstruction function
  procedure(i_reconstruct), deferred :: reconstruct
  !> Deferred implementation of the average over an interval
  procedure(i_average), deferred :: average
  !> Deferred implementation of evaluating the reconstruction at a point
  procedure(i_f), deferred :: f
  !> Deferred implementation of the derivative of the reconstruction at a point
  procedure(i_dfdx), deferred :: dfdx
  !> Deferred implementation of check_reconstruction
  !!
  !! Returns True if a check fails. Returns False if all checks pass.
  !! Checks are about internal, or inferred, state for arbitrary inputs.
  !! Checks should cover all the expected properties of a reconstruction.
  procedure(i_check_reconstruction), deferred :: check_reconstruction
  !> Deferred implementation of unit tests for the reconstruction
  !!
  !! Returns True if a test fails. Returns False if all tests pass.
  !! Tests in unit_tests() are usually checks against known (e.g. analytic) solutions.
  procedure(i_unit_tests), deferred :: unit_tests
  !> Deferred implementation of deallocation
  procedure(i_destroy), deferred :: destroy

  ! The following functions/subroutines are shared across all reconstructions and provided by this module
  ! unless replaced for the purpose of optimization

  !> Solves for x such that f(x)=t
  procedure :: x => x
  !> Remaps the column to subgrid h_sub
  procedure :: remap_to_sub_grid => remap_to_sub_grid
  !> Set debugging
  procedure :: set_debug => a_set_debug

  ! The following functions usually point to the same implementation as above but
  ! for derived secondary children these allow invocation of the parent class function.

  !> Second interface to init(), used to reach the primary class if derived from a primary implementation
  procedure(i_init_parent), deferred :: init_parent
  !> Second interface to reconstruct(), used to reach the primary class if derived from a primary implementation
  procedure(i_reconstruct_parent), deferred :: reconstruct_parent

end type Recon1d

interface

  !> Initialize a 1D reconstruction for n cells
  subroutine i_init(this, n, h_neglect, check)
    import :: Recon1d
    class(Recon1d),    intent(out) :: this !< This reconstruction
    integer,           intent(in)  :: n    !< Number of cells in this column
    real, optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
    logical, optional, intent(in)  :: check !< If true, enable some consistency checking
  end subroutine i_init

  !> Calculate a 1D reconstructions based on h(:) and u(:)
  subroutine i_reconstruct(this, h, u)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this !< This reconstruction
    real,           intent(in)    :: h(*) !< Grid spacing (thickness), typically in [H]
    real,           intent(in)    :: u(*) !< Cell mean values [A]
  end subroutine i_reconstruct

  !> Average between xa and xb for cell k of a 1D reconstruction [A]
  !!
  !! It is assumed that 0<=xa<=1, 0<=xb<=1, and xa<=xb
  real function i_average(this, k, xa, xb)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    integer,        intent(in) :: k    !< Cell number
    real,           intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
    real,           intent(in) :: xb   !< End of averaging interval on element (0 to 1)
  end function i_average

  !> Point-wise value of reconstruction [A]
  !!
  !! The function is only valid for 0 <= x <= 1. x is effectively clipped to this range.
  real function i_f(this, k, x)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    integer,        intent(in) :: k    !< Cell number
    real,           intent(in) :: x    !< Non-dimensional position within element [nondim]
  end function i_f

  !> Point-wise value of derivative reconstruction [A]
  !!
  !! The function is only valid for 0 <= x <= 1. x is effectively clipped to this range.
  real function i_dfdx(this, k, x)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    integer,        intent(in) :: k    !< Cell number
    real,           intent(in) :: x    !< Non-dimensional position within element [nondim]
  end function i_dfdx

  !> Point-wise solver for x: f(x)=t [nondim]
  !!
  !! The function solves for the non-dimensional position x within the cell where
  !! the reconstruction f(x)=t. The solver returns x=0 or x=1 if the target, t,
  !! is outside of the cell.
  real function i_x(this, k, t)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    integer,        intent(in) :: k    !< Cell number
    real,           intent(in) :: t    !< Value to solve for [A]
  end function i_x

  !> Returns true if some inconsistency is detected, false otherwise
  !!
  !! The nature of "consistency" is defined by the implementations
  !! and might be no-ops.
  logical function i_check_reconstruction(this, h, u)
    import :: Recon1d
    class(Recon1d), intent(in) :: this !< This reconstruction
    real,           intent(in) :: h(*) !< Grid spacing (thickness), typically in [H]
    real,           intent(in) :: u(*) !< Cell mean values [A]
  end function i_check_reconstruction

  !> Deallocate a 1D reconstruction
  subroutine i_destroy(this)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this !< This reconstruction
  end subroutine i_destroy

  !> Second interface to init(), or to parent init()
  subroutine i_init_parent(this, n, h_neglect, check)
    import :: Recon1d
    class(Recon1d), intent(out) :: this !< This reconstruction
    integer,        intent(in)  :: n    !< Number of cells in this column
    real, optional, intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
    logical, optional, intent(in)  :: check !< If true, enable some consistency checking
  end subroutine i_init_parent

  !> Second interface to reconstruct(), or to parent reconstruct()
  subroutine i_reconstruct_parent(this, h, u)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this !< This reconstruction
    real,           intent(in)    :: h(*) !< Grid spacing (thickness), typically in [H]
    real,           intent(in)    :: u(*) !< Cell mean values [A]
  end subroutine i_reconstruct_parent

  !> Runs reconstruction unit tests and returns True for any fails, False otherwise
  !!
  !! Assumes single process/thread context
  logical function i_unit_tests(this, verbose, stdout, stderr)
    import :: Recon1d
    class(Recon1d), intent(inout) :: this    !< This reconstruction
    logical,        intent(in)    :: verbose !< True, if verbose
    integer,        intent(in)    :: stdout  !< I/O channel for stdout
    integer,        intent(in)    :: stderr  !< I/O channel for stderr
  end function i_unit_tests

end interface

contains

!> Solve for x such that f(x)=t
!!
!! This solver uses bounded Newton-Raphson method with a fixed
!! number of iterations
real function x(this, k, t)
  class(Recon1d), intent(in) :: this !< This reconstruction
  integer,        intent(in) :: k    !< Cell number
  real,           intent(in) :: t    !< Value to solve for [A]
  real :: xl, xr, xo ! Left/right bounds and guess [nondim]
  real :: fl, fr ! Left right values [A]
  real :: slp ! Difference across cell or derivative wrt nondim x [A]
  real :: f_at_x ! Value at current x [A]
  integer :: iter

  x = 0.5 ! Fall back for special conditions
  fl = this%f(k, 0.)
  fr = this%f(k, 1.)
  slp = fr - fl
  if ( ( fl - t ) * ( t - fr ) > 0. ) then
    ! t is inside the range fl..fr
    xl = 0.
    xr = 1.
    xo = ( t - this%f(k, 0.) ) / slp ! First guess by regula falsi
    f_at_x = this%f(k, xo)
    do iter = 1,10
      slp = this%dfdx(k, xo)
      x = xo - ( f_at_x - t ) / slp ! Newton-Raphson step
      if ( x < xl ) x = 0.5 * ( xl + xo ) ! Replace with bi-section
      if ( x > xr ) x = 0.5 * ( xr + xo ) ! Replace with bi-section
      f_at_x = this%f(k, x)
      if ( abs(f_at_x - t) <= 0. .or. abs(x - xo) < this%x_tolerance ) return
      if ( f_at_x < t ) xl = x ! Replace left bound
      if ( f_at_x > t ) xr = x ! Replace right bound
      xo = x
    enddo
  elseif ( abs(slp) > 0. ) then
    slp = sign(1., slp)
    ! if t>u_mean & slp=1 then x=1
    ! if t<u_mean & slp=1 then x=0
    ! if t>u_mean & slp=-1 then x=0
    ! if t<u_mean & slp=-1 then x=1
    x = 0.5 + slp * sign(0.5, t - this%u_mean(k))
  else
    ! slp=0 so estimate "direction" from neighbors
    slp = this%f(min(k+1,this%n), 0.) - this%f(max(k-1,1), 1.)
    if ( abs(slp) > 0. ) slp = sign(1., slp)
    ! if t>u_mean & slp=1 then x=1
    ! if t<u_mean & slp=1 then x=0
    ! if t>u_mean & slp=-1 then x=0
    ! if t<u_mean & slp=-1 then x=1
    ! if t=u_mean then x=0.5
    ! if slp=0 then x=0.5
    if ( abs(t - this%u_mean(k)) > 0. ) x = 0.5 + slp * sign(0.5, t - this%u_mean(k))
  endif
end function x

!> Remaps the column to subgrid h_sub
!!
!! It is assumed that h_sub is a perfect sub-grid of h0, meaning each h0 cell
!! can be constructed by joining a contiguous set of h_sub cells. The integer
!! indices isrc_start, isrc_end, isub_src provide this mapping, and are
!! calculated in MOM_remapping
subroutine remap_to_sub_grid(this, h0, u0, n1, h_sub, &
                                   isrc_start, isrc_end, isrc_max, isub_src, &
                                   u_sub, uh_sub, u02_err)
  class(Recon1d), intent(in) :: this !< 1-D reconstruction type
  real,    intent(in)  :: h0(*)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: u0(*)  !< Source grid widths (size n0) [H]
  integer, intent(in)  :: n1      !< Number of cells in target grid
  real,    intent(in)  :: h_sub(*) !< Overlapping sub-cell thicknesses, h_sub [H]
  integer, intent(in)  :: isrc_start(*) !< Index of first sub-cell within each source cell
  integer, intent(in)  :: isrc_end(*) !< Index of last sub-cell within each source cell
  integer, intent(in)  :: isrc_max(*) !< Index of thickest sub-cell within each source cell
  integer, intent(in)  :: isub_src(*) !< Index of source cell for each sub-cell
  real,    intent(out) :: u_sub(*) !< Sub-cell cell averages (size n1) [A]
  real,    intent(out) :: uh_sub(*) !< Sub-cell cell integrals (size n1) [A H]
  real,    intent(out) :: u02_err !< Integrated reconstruction error estimates [A H]
  ! Local variables
  integer :: i_sub ! Index of sub-cell
  integer :: i0 ! Index into h0(1:n0), source column
  integer :: i_max ! Used to record which sub-cell is the largest contribution of a source cell
  real :: dh_max ! Used to record which sub-cell is the largest contribution of a source cell [H]
  real :: xa, xb ! Non-dimensional position within a source cell (0..1) [nondim]
  real :: dh ! The width of the sub-cell [H]
  real :: duh ! The total amount of accumulated stuff (u*h) [A H]
  real :: dh0_eff ! Running sum of source cell thickness [H]
  integer :: i0_last_thick_cell, n0
! real :: u0_min(this%n), u0_max(this%n) ! Min/max of u0 for each source cell [A]
! real :: ul,ur ! Left/right edge values [A]

  n0 = this%n

  i0_last_thick_cell = 0
  do i0 = 1, n0
!   ul = this%f(i0, 0.)
!   ur = this%f(i0, 1.)
!   u0_min(i0) = min(ul, ur)
!   u0_max(i0) = max(ul, ur)
    if (h0(i0)>0.) i0_last_thick_cell = i0
  enddo

  ! Loop over each sub-cell to calculate average/integral values within each sub-cell.
  ! Uses: h_sub, isub_src, h0_eff
  ! Sets: u_sub, uh_sub
  xa = 0.
  dh0_eff = 0.
  u02_err = 0.
  do i_sub = 1, n0+n1

    ! Sub-cell thickness from loop above
    dh = h_sub(i_sub)

    ! Source cell
    i0 = isub_src(i_sub)

    ! Evaluate average and integral for sub-cell i_sub.
    ! Integral is over distance dh but expressed in terms of non-dimensional
    ! positions with source cell from xa to xb  (0 <= xa <= xb <= 1).
    dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
    if (h0(i0)>0.) then
      xb = dh0_eff / h0(i0) ! This expression yields xa <= xb <= 1.0
      xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
      u_sub(i_sub) = this%average( i0, xa, xb )
    else ! Vanished cell
      xb = 1.
      u_sub(i_sub) = u0(i0)
    endif
!   u_sub(i_sub) = max( u_sub(i_sub), u0_min(i0) )
!   u_sub(i_sub) = min( u_sub(i_sub), u0_max(i0) )
    uh_sub(i_sub) = dh * u_sub(i_sub)

    if (isub_src(i_sub+1) /= i0) then
      ! If the next sub-cell is in a different source cell, reset the position counters
      dh0_eff = 0.
      xa = 0.
    else
      xa = xb ! Next integral will start at end of last
    endif

  enddo
  i_sub = n0+n1+1
  ! Sub-cell thickness from loop above
  dh = h_sub(i_sub)
  ! Source cell
  i0 = isub_src(i_sub)

  ! Evaluate average and integral for sub-cell i_sub.
  ! Integral is over distance dh but expressed in terms of non-dimensional
  ! positions with source cell from xa to xb  (0 <= xa <= xb <= 1).
  dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
  if (h0(i0)>0.) then
    xb = dh0_eff / h0(i0) ! This expression yields xa <= xb <= 1.0
    xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
    u_sub(i_sub) = this%average( i0, xa, xb )
  else ! Vanished cell
    xb = 1.
    u_sub(i_sub) = u0(i0)
  endif
! u_sub(i_sub) = max( u_sub(i_sub), u0_min(i0) )
! u_sub(i_sub) = min( u_sub(i_sub), u0_max(i0) )
  uh_sub(i_sub) = dh * u_sub(i_sub)

  ! Loop over each source cell substituting the integral/average for the thickest sub-cell (within
  ! the source cell) with the residual of the source cell integral minus the other sub-cell integrals
  ! aka a genius algorithm for accurate conservation when remapping from Robert Hallberg (\@Hallberg-NOAA).
  ! Uses: i0_last_thick_cell, isrc_max, h_sub, isrc_start, isrc_end, uh_sub, u0, h0
  ! Updates: uh_sub
  do i0 = 1, i0_last_thick_cell
    i_max = isrc_max(i0)
    dh_max = h_sub(i_max)
    if (dh_max > 0.) then
      ! duh will be the sum of sub-cell integrals within the source cell except for the thickest sub-cell.
      duh = 0.
      do i_sub = isrc_start(i0), isrc_end(i0)
        if (i_sub /= i_max) duh = duh + uh_sub(i_sub)
      enddo
      uh_sub(i_max) = u0(i0)*h0(i0) - duh
      u02_err = u02_err + max( abs(uh_sub(i_max)), abs(u0(i0)*h0(i0)), abs(duh) )
    endif
  enddo

  ! This should not generally be used
  if (this%check) then
    if ( this%check_reconstruction(h0, u0) ) stop 910 ! A debugger is required to understand why this failed
  endif

end subroutine remap_to_sub_grid

!> Turns on debugging
subroutine a_set_debug(this)
  class(Recon1d), intent(inout) :: this !< 1-D reconstruction type

  this%debug = .true.

end subroutine a_set_debug

!> \namespace recon1d_type
!!
!! \section section_recon1d_type Generic vertical reconstruction type
!!
!! A class to describe generic reconstruction in 1-D. This module has no implementations
!! but defines the interfaces for members that implement a reconstruction.
!!
!! e.g. a chain of derived reconstructions might look like
!!   Recon1d_type <- Recond1d_XYZ <- Recon1d_XYZ_v2
!! where
!!   Recon1d_type      - defines the interfaces (this module)
!!   Recon1d_XYZ       - extends Recon1d_type, implements the XYZ reconstruction in reconstruct(),
!!                       and reconstruc_parent() -> reconstruct() of the same Recon1d_XYZ module
!!   Recon1d_XYZ_v2    - implements a slight variant of Recon1d_XYZ via reconstruct()
!!                       but reconstruc_parent() is not redefined so that it still is defined by Recon1d_XYZ
!!
!! The schemes that use this structure are described in \ref Vertical_Reconstruction
end module Recon1d_type
