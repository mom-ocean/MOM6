!> Monotonized Piecewise Linear Method 1D reconstruction
!!
!! This implementation of PLM follows White and Adcroft, 2008 \cite white2008.
!! The PLM slopes are first limited following Colella and Woodward, 1984, but are then
!! further limited to ensure the edge values moving across cell boundaries are monotone.
!! The first and last cells are always limited to PCM.
!!
!! This differs from recon1d_mplm_wa_poly in the internally not polynomial representations
!! are referred to.
module Recon1d_MPLM_WA

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_PLM_CW, only : PLM_CW, testing

implicit none ; private

public MPLM_WA, testing

!> Limited Monotonic PLM reconstruction following White and Adcroft, 2008
!!
!! The source for the methods ultimately used by this class are:
!! - init()                 -> recon1d_plm_cw.init()
!! - reconstruct()             *locally defined
!! - average()              -> recon1d_plm_cw.average()
!! - f()                    -> recon1d_plm_cw.f()
!! - dfdx()                 -> recon1d_plm_cw.dfdx()
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()              -> recon1d_plm_cw.destroy()
!! - remap_to_sub_grid()    -> recon1d_plm_cw           -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> recon1d_plm_cw.init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (PLM_CW) :: MPLM_WA

contains
  !> Implementation of the MPLM_WA reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of check reconstruction for the MPLM_WA reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the MPLM_WA reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type MPLM_WA

contains

!> Calculate a 1D PLM reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(MPLM_WA), intent(inout) :: this !< This reconstruction
  real,           intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real :: slp(this%n) ! The PLM slopes (difference across cell) [A]
  real :: mslp(this%n) ! The monotonized PLM slopes [A]
  integer :: k, n
  real :: u_tmp, u_min, u_max ! Working values of cells [A]

  n = this%n

  ! Loop over all cells
  do k = 1, n
    this%u_mean(k) = u(k)
  enddo

  ! Loop on interior cells
  do k = 2, n-1
    slp(k) = PLM_slope_wa(h(k-1), h(k), h(k+1), this%h_neglect, u(k-1), u(k), u(k+1))
  enddo ! end loop on interior cells

  ! Boundary cells use PCM. Extrapolation is handled after monotonization.
  slp(1) = 0.
  slp(n) = 0.

  ! This loop adjusts the slope so that edge values are monotonic.
  do k = 2, n-1
    mslp(k) = PLM_monotonized_slope( u(k-1), u(k), u(k+1), slp(k-1), slp(k), slp(k+1) )
  enddo ! end loop on interior cells
  mslp(1) = 0.
  mslp(n) = 0.

  ! Store edge values
  this%ul(1) = u(1)
  this%ur(1) = u(1)
  do k = 2, n-1
    u_tmp = u(k-1) + 0.5 * mslp(k-1) ! Right edge value of cell k-1
    u_min = min( u(k), u_tmp )
    u_max = max( u(k), u_tmp )
    u_tmp = u(k) - 0.5 * mslp(k) ! Left edge value of cell k
    this%ul(k) = max( min( u_tmp, u_max), u_min ) ! Bounded to handle roundoff
    u_tmp = u(k+1) - 0.5 * mslp(k-1) ! Left edge value of cell k+1
    u_min = min( u(k), u_tmp )
    u_max = max( u(k), u_tmp )
    u_tmp = u(k) + 0.5 * mslp(k) ! Right edge value of cell k
    this%ur(k) = max( min( u_tmp, u_max), u_min ) ! Bounded to handle roundoff
  enddo
  this%ul(n) = u(n)
  this%ur(n) = u(n)

end subroutine reconstruct

!> Returns a limited PLM slope following White and Adcroft, 2008, in the same arbitrary
!! units [A] as the input values.
!! Note that this is not the same as the Colella and Woodward method.
real elemental pure function PLM_slope_wa(h_l, h_c, h_r, h_neglect, u_l, u_c, u_r)
  real, intent(in) :: h_l !< Thickness of left cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_c !< Thickness of center cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_r !< Thickness of right cell in arbitrary grid thickness units [H]
  real, intent(in) :: h_neglect !< A negligible thickness [H]
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  real, intent(in) :: u_r !< Value of right cell in arbitrary units [A]
  ! Local variables
  real :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
                                    ! differences across the cell [A]
  real :: u_min, u_max ! Minimum and maximum value across cell [A]

  ! Side differences
  sigma_r = u_r - u_c
  sigma_l = u_c - u_l

  ! Quasi-second order difference
  sigma_c = 2.0 * ( u_r - u_l ) * ( h_c / ( h_l + 2.0*h_c + h_r + h_neglect) )

  ! Limit slope so that reconstructions are bounded by neighbors
  u_min = min( u_l, u_c, u_r )
  u_max = max( u_l, u_c, u_r )
  if ( (sigma_l * sigma_r) > 0.0 ) then
    ! This limits the slope so that the edge values are bounded by the
    ! two cell averages spanning the edge.
    PLM_slope_wa = sign( min( abs(sigma_c), 2.*min( u_c - u_min, u_max - u_c ) ), sigma_c )
  else
    ! Extrema in the mean values require a PCM reconstruction avoid generating
    ! larger extreme values.
    PLM_slope_wa = 0.0
  endif

end function PLM_slope_wa

!> Returns a limited PLM slope following Colella and Woodward 1984, in the same
!! arbitrary units as the input values [A].
real elemental pure function PLM_monotonized_slope(u_l, u_c, u_r, s_l, s_c, s_r)
  real, intent(in) :: u_l !< Value of left cell in arbitrary units [A]
  real, intent(in) :: u_c !< Value of center cell in arbitrary units [A]
  real, intent(in) :: u_r !< Value of right cell in arbitrary units [A]
  real, intent(in) :: s_l !< PLM slope of left cell [A]
  real, intent(in) :: s_c !< PLM slope of center cell [A]
  real, intent(in) :: s_r !< PLM slope of right cell [A]
  ! Local variables
  real :: neighbor_edge ! Edge value of nieghbor cell [A]
  real :: this_edge ! Edge value of this cell [A]
  real :: slp ! Magnitude of PLM central slope [A]

  ! Comparison are made assuming +ve slopes
  slp = abs(s_c)

  ! Check that left edge is between right edge of cell to the left and this cell mean
  neighbor_edge = u_l + 0.5 * s_l
  this_edge = u_c - 0.5 * s_c
  if ( ( this_edge - neighbor_edge ) * ( u_c - this_edge ) < 0. ) then
    ! Using the midpoint works because the neighbor is similarly adjusted
    this_edge = 0.5 * ( this_edge + neighbor_edge )
    slp = min( slp, abs( this_edge - u_c ) * 2. )
  endif

  ! Check that right edge is between left edge of cell to the right and this cell mean
  neighbor_edge = u_r - 0.5 * s_r
  this_edge = u_c + 0.5 * s_c
  if ( ( this_edge - u_c ) * ( neighbor_edge - this_edge ) < 0. ) then
    ! Using the midpoint works because the neighbor is similarly adjusted
    this_edge = 0.5 * ( this_edge + neighbor_edge )
    slp = min( slp, abs( this_edge - u_c ) * 2. )
  endif

  PLM_monotonized_slope = sign( slp, s_c )

end function PLM_monotonized_slope

!> Checks the MPLM_WA reconstruction for consistency
logical function check_reconstruction(this, h, u)
  class(MPLM_WA), intent(in) :: this !< This reconstruction
  real,           intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real,           intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k

  check_reconstruction = .false.

  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0. ) check_reconstruction = .true.
  enddo

  ! Check the cell reconstruction is monotonic within each cell (it should be as a straight line)
  do k = 1, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ur(k) - this%u_mean(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! This next test fails abysmally!
  ! Using intel/2023.2.0 on gaea, MOM_remapping:test_recon_consistency iter=6
  !    um~0.581492556923472  ul~0.402083491713151 ur~0.749082615698503
  ! Check the cell is a straight line (to within machine precision)
! do k = 1, this%n
!   if ( abs(2. * this%u_mean(k) - ( this%ul(k) + this%ur(k) )) > epsilon(this%u_mean(1)) * &
!        max(abs(2. * this%u_mean(k)), abs(this%ul(k)), abs(this%ur(k))) ) check_reconstruction = .true.
! enddo

  ! Check bounding of right edges, w.r.t. the cell means
  do K = 1, this%n-1
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%u_mean(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check bounding of left edges, w.r.t. the cell means
  do K = 2, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%u_mean(k-1) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check order of u, ur, ul
  do K = 1, this%n-1
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%ul(k+1) - this%ur(k) ) < 0. ) check_reconstruction = .true.
  enddo

  ! Check bounding of left edges, w.r.t. this cell mean and the previous cell right edge
  do K = 2, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%ur(k-1) ) < 0. ) check_reconstruction = .true.
  enddo

end function check_reconstruction

!> Runs PLM reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(MPLM_WA), intent(inout) :: this    !< This reconstruction
  logical,        intent(in)    :: verbose !< True, if verbose
  integer,        intent(in)    :: stdout  !< I/O channel for stdout
  integer,        intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables
  real, allocatable :: ul(:), ur(:), um(:) ! test values [A]
  real, allocatable :: ull(:), urr(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k

  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  call this%init(3)
  call test%test( this%n /= 3, 'Setting number of levels')
  allocate( um(3), ul(3), ur(3), ull(3), urr(3) )

  call this%reconstruct( (/2.,2.,2./), (/1.,3.,5./) )
  call test%real_arr(3, this%u_mean, (/1.,3.,5./), 'Setting cell values')

  do k = 1, 3
    ul(k) = this%f(k, 0.)
    um(k) = this%f(k, 0.5)
    ur(k) = this%f(k, 1.)
  enddo
  call test%real_arr(3, ul, (/1.,2.,5./), 'Evaluation on left edge')
  call test%real_arr(3, um, (/1.,3.,5./), 'Evaluation in center')
  call test%real_arr(3, ur, (/1.,4.,5./), 'Evaluation on right edge')

  do k = 1, 3
    ul(k) = this%dfdx(k, 0.)
    um(k) = this%dfdx(k, 0.5)
    ur(k) = this%dfdx(k, 1.)
  enddo
  call test%real_arr(3, ul, (/0.,2.,0./), 'dfdx on left edge')
  call test%real_arr(3, um, (/0.,2.,0./), 'dfdx in center')
  call test%real_arr(3, ur, (/0.,2.,0./), 'dfdx on right edge')

  do k = 1, 3
    um(k) = this%average(k, 0.5, 0.75) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(3, um, (/1.,3.25,5./), 'Return interval average')

  unit_tests = test%summarize('MPLM_WA:unit_tests')

end function unit_tests

!> \namespace recon1d_mplm_wa
!!

end module Recon1d_MPLM_WA
