!> Provides column-wise vertical remapping functions
module MOM_remapping

! This file is part of MOM6. See LICENSE.md for the license.
! Original module written by Laurent White, 2008.06.09

use MOM_error_handler, only : MOM_error, FATAL
use MOM_string_functions, only : uppercase
use numerical_testing_type, only : testing
use regrid_edge_values, only : edge_values_explicit_h4, edge_values_implicit_h4
use regrid_edge_values, only : edge_values_explicit_h4cw
use regrid_edge_values, only : edge_values_implicit_h4, edge_values_implicit_h6
use regrid_edge_values, only : edge_slopes_implicit_h3, edge_slopes_implicit_h5
use PCM_functions, only : PCM_reconstruction
use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation
use PPM_functions, only : PPM_monotonicity
use PQM_functions, only : PQM_reconstruction, PQM_boundary_extrapolation_v1
use MOM_hybgen_remap, only : hybgen_plm_coefs, hybgen_ppm_coefs, hybgen_weno_coefs

use Recon1d_type, only : Recon1d
use Recon1d_PCM, only : PCM
use Recon1d_PLM_CW, only : PLM_CW
use Recon1d_PLM_hybgen, only : PLM_hybgen
use Recon1d_PLM_CWK, only : PLM_CWK
use Recon1d_MPLM_CWK, only : MPLM_CWK
use Recon1d_EMPLM_CWK, only : EMPLM_CWK
use Recon1d_MPLM_WA, only : MPLM_WA
use Recon1d_EMPLM_WA, only : EMPLM_WA
use Recon1d_MPLM_WA_poly, only : MPLM_WA_poly
use Recon1d_EMPLM_WA_poly, only : EMPLM_WA_poly
use Recon1d_PPM_CW, only : PPM_CW
use Recon1d_PPM_hybgen, only : PPM_hybgen
use Recon1d_PPM_CWK, only : PPM_CWK
use Recon1d_EPPM_CWK, only : EPPM_CWK
use Recon1d_PPM_H4_2019, only : PPM_H4_2019
use Recon1d_PPM_H4_2018, only : PPM_H4_2018

implicit none ; private

!> Container for remapping parameters
type, public :: remapping_CS ; private
  !> Determines which reconstruction to use
  integer :: remapping_scheme = -911
  !> Degree of polynomial reconstruction
  integer :: degree = 0
  !> If true, extrapolate boundaries
  logical :: boundary_extrapolation = .true.
  !> If true, reconstructions are checked for consistency.
  logical :: check_reconstruction = .false.
  !> If true, the result of remapping are checked for conservation and bounds.
  logical :: check_remapping = .false.
  !> If true, the intermediate values used in remapping are forced to be bounded.
  logical :: force_bounds_in_subcell = .false.
  !> If true, impose bounds on the remapping from sub-cells to target grid
  logical :: force_bounds_in_target = .true.
  !> If true, impose bounds on the remapping from non-vanished sub-cells to target grid
  logical :: better_force_bounds_in_target = .false.
  !> If true, calculate and use an offset when summing sub-cells to the target grid
  logical :: offset_tgt_summation = .false.
  !> The vintage of the expressions to use for remapping. Values below 20190101 result
  !! in the use of older, less accurate expressions.
  integer :: answer_date = 99991231
  !> If true, use the OM4 version of the remapping algorithm that makes poor assumptions
  !! about the reconstructions in top and bottom layers of the source grid
  logical :: om4_remap_via_sub_cells = .false.

  !> A negligibly small width for the purpose of cell reconstructions in the same units
  !! as the h0 argument to remapping_core_h [H]
  real :: h_neglect
  !> A negligibly small width for the purpose of edge value calculations in the same units
  !! as the h0 argument to remapping_core_h [H]
  real :: h_neglect_edge

  !> If true, do some debugging as operations proceed
  logical :: debug = .false.

  !> The instance of the actual equation of state
  class(Recon1d), pointer :: reconstruction => Null()
end type

! The following routines are visible to the outside world
public remapping_core_h, remapping_core_w
public initialize_remapping, end_remapping, remapping_set_param, extract_member_remapping_CS
public remapping_unit_tests, build_reconstructions_1d, average_value_ppoly
public interpolate_column, reintegrate_column, dzFromH1H2

! The following are private parameter constants
integer, parameter  :: REMAPPING_PCM        = 0 !< O(h^1) remapping scheme
integer, parameter  :: REMAPPING_PLM        = 2 !< O(h^2) remapping scheme
integer, parameter  :: REMAPPING_PLM_HYBGEN = 3 !< O(h^2) remapping scheme
integer, parameter  :: REMAPPING_PPM_CW     =10 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PPM_H4     = 4 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PPM_IH4    = 5 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PPM_HYBGEN = 6 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_WENO_HYBGEN= 7 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PQM_IH4IH3 = 8 !< O(h^4) remapping scheme
integer, parameter  :: REMAPPING_PQM_IH6IH5 = 9 !< O(h^5) remapping scheme
integer, parameter  :: REMAPPING_VIA_CLASS  =99 !< Scheme is controlled by Recon1d class

integer, parameter  :: INTEGRATION_PCM = 0  !< Piecewise Constant Method
integer, parameter  :: INTEGRATION_PLM = 1  !< Piecewise Linear Method
integer, parameter  :: INTEGRATION_PPM = 3  !< Piecewise Parabolic Method
integer, parameter  :: INTEGRATION_PQM = 5  !< Piecewise Quartic Method

!> Documentation for external callers
character(len=360), public :: remappingSchemesDoc = &
                 "PCM         (1st-order accurate)\n"//&
                 "PLM         (2nd-order accurate)\n"//&
                 "PLM_HYBGEN  (2nd-order accurate)\n"//&
                 "PPM_H4      (3rd-order accurate)\n"//&
                 "PPM_IH4     (3rd-order accurate)\n"//&
                 "PPM_HYBGEN  (3rd-order accurate)\n"//&
                 "WENO_HYBGEN (3rd-order accurate)\n"//&
                 "PQM_IH4IH3  (4th-order accurate)\n"//&
                 "PQM_IH6IH5  (5th-order accurate)\n"
character(len=3), public :: remappingDefaultScheme = "PLM" !< Default remapping method

contains

!> Set parameters within remapping object
subroutine remapping_set_param(CS, remapping_scheme, boundary_extrapolation,  &
               check_reconstruction, check_remapping, force_bounds_in_subcell, &
               force_bounds_in_target, better_force_bounds_in_target, offset_tgt_summation, &
               om4_remap_via_sub_cells, answers_2018, answer_date, nk, &
               h_neglect, h_neglect_edge)
  type(remapping_CS),         intent(inout) :: CS !< Remapping control structure
  character(len=*), optional, intent(in)    :: remapping_scheme !< Remapping scheme to use
  logical, optional,          intent(in)    :: boundary_extrapolation !< Indicate to extrapolate in boundary cells
  logical, optional,          intent(in)    :: check_reconstruction !< Indicate to check reconstructions
  logical, optional,          intent(in)    :: check_remapping !< Indicate to check results of remapping
  logical, optional,          intent(in)    :: force_bounds_in_subcell !< Force subcells values to be bounded
  logical, optional,          intent(in)    :: force_bounds_in_target !< Force target values to be bounded
  logical, optional,          intent(in)    :: better_force_bounds_in_target !< Force target values to be bounded
  logical, optional,          intent(in)    :: offset_tgt_summation !< Use an offset when summing sub-cells
  logical, optional,          intent(in)    :: om4_remap_via_sub_cells !< If true, use OM4 remapping algorithm
  logical, optional,          intent(in)    :: answers_2018 !< If true use older, less accurate expressions.
  integer, optional,          intent(in)    :: answer_date  !< The vintage of the expressions to use
  real,    optional,          intent(in)    :: h_neglect !< A negligibly small width for the purpose of cell
                                                         !! reconstructions in the same units as the h0 argument
                                                         !! to remapping_core_h [H]
  real,    optional,          intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose of edge
                                                         !! value calculations in the same units as as the h0
                                                         !! argument to remapping_core_h [H]
  integer, optional,          intent(in)    :: nk !< Number of levels to initialize reconstruction class with

  if (present(remapping_scheme)) then
    call setReconstructionType( remapping_scheme, CS )
    if (index(trim(remapping_scheme),'C_')>0) then
      if (present(nk)) then
        call CS%reconstruction%init(nk, h_neglect=h_neglect)
      else
        call MOM_error( FATAL, 'MOM_remapping, remapping_set_param: '//&
           'Using the Recon1d class for remapping requires nk to be passed' )
      endif
    endif
  endif
  if (present(boundary_extrapolation)) then
    CS%boundary_extrapolation = boundary_extrapolation
  endif
  if (present(check_reconstruction)) then
    CS%check_reconstruction = check_reconstruction
  endif
  if (present(check_remapping)) then
    CS%check_remapping = check_remapping
  endif
  if (present(force_bounds_in_subcell)) then
    CS%force_bounds_in_subcell = force_bounds_in_subcell
  endif
  if (present(force_bounds_in_target)) then
    CS%force_bounds_in_target = force_bounds_in_target
  endif
  if (present(better_force_bounds_in_target)) then
    CS%better_force_bounds_in_target = better_force_bounds_in_target
  endif
  if (present(offset_tgt_summation)) then
    CS%offset_tgt_summation = offset_tgt_summation
  endif
  if (present(om4_remap_via_sub_cells)) then
    CS%om4_remap_via_sub_cells = om4_remap_via_sub_cells
  endif
  if (present(answers_2018)) then
    if (answers_2018) then
      CS%answer_date = 20181231
    else
      CS%answer_date = 20190101
    endif
  endif
  if (present(answer_date)) then
    CS%answer_date = answer_date
  endif
  if (present(h_neglect)) then
    CS%h_neglect = h_neglect
  endif
  if (present(h_neglect_edge)) then
    CS%h_neglect_edge = h_neglect_edge
  endif

end subroutine remapping_set_param

subroutine extract_member_remapping_CS(CS, remapping_scheme, degree, boundary_extrapolation, check_reconstruction, &
                                       check_remapping, force_bounds_in_subcell, force_bounds_in_target, &
                                       better_force_bounds_in_target, offset_tgt_summation)
  type(remapping_CS), intent(in) :: CS !< Control structure for remapping module
  integer, optional, intent(out) :: remapping_scheme        !< Determines which reconstruction scheme to use
  integer, optional, intent(out) :: degree                  !< Degree of polynomial reconstruction
  logical, optional, intent(out) :: boundary_extrapolation  !< If true, extrapolate boundaries
  logical, optional, intent(out) :: check_reconstruction    !< If true, reconstructions are checked for consistency.
  logical, optional, intent(out) :: check_remapping         !< If true, the result of remapping are checked
                                                            !!  for conservation and bounds.
  logical, optional, intent(out) :: force_bounds_in_subcell !< If true, the intermediate values used in
                                                            !! remapping are forced to be bounded.
  logical, optional, intent(out) :: force_bounds_in_target  !< Force target values to be bounded
  logical, optional, intent(out) :: better_force_bounds_in_target  !< Force target values to be bounded
  logical, optional, intent(out) :: offset_tgt_summation    !< Use an offset when summing sub-cells

  if (present(remapping_scheme)) remapping_scheme = CS%remapping_scheme
  if (present(degree)) degree = CS%degree
  if (present(boundary_extrapolation)) boundary_extrapolation = CS%boundary_extrapolation
  if (present(check_reconstruction)) check_reconstruction = CS%check_reconstruction
  if (present(check_remapping)) check_remapping = CS%check_remapping
  if (present(force_bounds_in_subcell)) force_bounds_in_subcell = CS%force_bounds_in_subcell
  if (present(force_bounds_in_target)) force_bounds_in_target = CS%force_bounds_in_target
  if (present(better_force_bounds_in_target)) better_force_bounds_in_target = CS%better_force_bounds_in_target
  if (present(offset_tgt_summation)) offset_tgt_summation = CS%offset_tgt_summation

end subroutine extract_member_remapping_CS

!> Remaps column of values u0 on grid h0 to grid h1 assuming the top edge is aligned and using the OM4
!! reconstruction methods
!!
!! \todo Remove h_neglect argument by moving into remapping_CS
!! \todo Remove PCM_cell argument by adding new method in Recon1D class
subroutine remapping_core_h(CS, n0, h0, u0, n1, h1, u1, net_err, PCM_cell)
  type(remapping_CS),  intent(in)  :: CS !< Remapping control structure
  integer,             intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0), intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0), intent(in)  :: u0 !< Cell averages on source grid [A]
  integer,             intent(in)  :: n1 !< Number of cells on target grid
  real, dimension(n1), intent(in)  :: h1 !< Cell widths on target grid [H]
  real, dimension(n1), intent(out) :: u1 !< Cell averages on target grid [A]
  real, optional,      intent(out) :: net_err !< Error in total column [A H]
  logical, dimension(n0), optional, intent(in) :: PCM_cell !< If present, use PCM remapping for
                                         !! cells in the source grid where this is true.
  ! Local variables
  real, dimension(n0+n1+1) :: h_sub ! Width of each each sub-cell [H]
  real, dimension(n0+n1+1) :: uh_sub ! Integral of u*h over each sub-cell [A H]
  real, dimension(n0+n1+1) :: u_sub ! Average of u over each sub-cell [A]
  integer, dimension(n0+n1+1) :: isub_src ! Index of source cell for each sub-cell
  integer, dimension(n0) :: isrc_start ! Index of first sub-cell within each source cell
  integer, dimension(n0) :: isrc_end ! Index of last sub-cell within each source cell
  integer, dimension(n0) :: isrc_max ! Index of thickest sub-cell within each source cell
  real, dimension(n0) :: h0_eff ! Effective thickness of source cells [H]
  integer, dimension(n1) :: itgt_start ! Index of first sub-cell within each target cell
  integer, dimension(n1) :: itgt_end ! Index of last sub-cell within each target cell
  ! For error checking/debugging
  real :: u02_err ! Integrated reconstruction error estimates [H A]
  real, dimension(n0,2)           :: ppoly_r_E     ! Edge value of polynomial [A]
  real, dimension(n0,2)           :: ppoly_r_S     ! Edge slope of polynomial [A H-1]
  real, dimension(n0,CS%degree+1) :: ppoly_r_coefs ! Coefficients of polynomial reconstructions [A]
  real :: uh_err       ! A bound on the error in the sum of u*h, as estimated by the remapping code [A H]
  integer :: iMethod   ! An integer indicating the integration method used

  ! Calculate sub-layer thicknesses and indices connecting sub-layers to source and target grids
  ! Sets: h_sub, h0_eff, isrc_start, isrc_end, isrc_max, isub_src, itgt_start, itgt_end
  call intersect_src_tgt_grids(n0, h0, n1, h1, h_sub, h0_eff, &
                               isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src)

  if (CS%remapping_scheme == REMAPPING_VIA_CLASS) then

!   if (CS%debug) call CS%reconstruction%set_debug() ! Sets an internal flag

    call CS%reconstruction%reconstruct(h0, u0)

    ! Adjust h_sub so that the Hallberg conservation trick works properly
!   call adjust_h_sub( n0, h0, n1, isrc_start, isrc_end, isrc_max, h_sub )

    ! Loop over each sub-cell to calculate average/integral values within each sub-cell.
    ! Uses: h_sub, isrc_start, isrc_end, isrc_max, isub_src
    ! Sets: u_sub, uh_sub
    call CS%reconstruction%remap_to_sub_grid(h0, u0, n1, h_sub, &
                                             isrc_start, isrc_end, isrc_max, isub_src, &
                                             u_sub, uh_sub, u02_err)

    ! Loop over each target cell summing the integrals from sub-cells within the target cell.
    ! Uses: itgt_start, itgt_end, h1, h_sub, uh_sub, u_sub
    ! Sets: u1, uh_err
    call remap_sub_to_tgt_grid(n0, n1, h1, h_sub, u_sub, uh_sub, itgt_start, itgt_end, &
                               CS%force_bounds_in_target, CS%offset_tgt_summation, &
                               CS%better_force_bounds_in_target, u1, uh_err)

    ! Include the error remapping from source to sub-cells in the estimate of total remapping error
    uh_err = uh_err + u02_err

  else ! Uses the OM4-era reconstruction functions

    call build_reconstructions_1d(CS, n0, h0, u0, ppoly_r_coefs, ppoly_r_E, ppoly_r_S, iMethod, &
                                  CS%h_neglect, CS%h_neglect_edge, PCM_cell, debug=CS%debug)

    if (CS%check_reconstruction) call check_reconstructions_1d(n0, h0, u0, CS%degree, &
                                   CS%boundary_extrapolation, ppoly_r_coefs, ppoly_r_E)

    ! Loop over each sub-cell to calculate average/integral values within each sub-cell.
    ! Uses: h_sub, h0_eff, isub_src
    ! Sets: u_sub, uh_sub
    if (CS%om4_remap_via_sub_cells) then ! Uses the version from OM4 with a bug at the bottom

      call remap_src_to_sub_grid_om4(n0, h0, u0, ppoly_r_E, ppoly_r_coefs, n1, h_sub, &
                                     h0_eff, isrc_start, isrc_end, isrc_max, isub_src, &
                                     iMethod, CS%force_bounds_in_subcell, u_sub, uh_sub, u02_err)

    else ! i.e. if (CS%om4_remap_via_sub_cells == .false.)

      call remap_src_to_sub_grid(n0, h0, u0, ppoly_r_E, ppoly_r_coefs, n1, h_sub, &
                                 isrc_start, isrc_end, isrc_max, isub_src, &
                                 iMethod, CS%force_bounds_in_subcell, u_sub, uh_sub, u02_err)

    endif

    ! Loop over each target cell summing the integrals from sub-cells within the target cell.
    ! Uses: itgt_start, itgt_end, h1, h_sub, uh_sub, u_sub
    ! Sets: u1, uh_err
    call remap_sub_to_tgt_grid_om4(n0, n1, h1, h_sub, u_sub, uh_sub, itgt_start, itgt_end, &
                               CS%force_bounds_in_target, u1, uh_err)
    ! Include the error remapping from source to sub-cells in the estimate of total remapping error
    uh_err = uh_err + u02_err

    if (CS%check_remapping) call check_remapped_values(n0, h0, u0, ppoly_r_E, CS%degree, ppoly_r_coefs, &
                                                       n1, h1, u1, iMethod, uh_err, "remapping_core_h")

  endif

 if (present(net_err)) net_err = uh_err

end subroutine remapping_core_h

!> Remaps column of values u0 on grid h0 to implied grid h1
!! where the interfaces of h1 differ from those of h0 by dx.
subroutine remapping_core_w( CS, n0, h0, u0, n1, dx, u1)
  type(remapping_CS),    intent(in)  :: CS !< Remapping control structure
  integer,               intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),   intent(in)  :: u0 !< Cell averages on source grid [A]
  integer,               intent(in)  :: n1 !< Number of cells on target grid
  real, dimension(n1+1), intent(in)  :: dx !< Cell widths on target grid [H]
  real, dimension(n1),   intent(out) :: u1 !< Cell averages on target grid [A]

  ! Local variables
  real, dimension(n0+n1+1) :: h_sub ! Width of each each sub-cell [H]
  real, dimension(n0+n1+1) :: uh_sub ! Integral of u*h over each sub-cell [A H]
  real, dimension(n0+n1+1) :: u_sub ! Average of u over each sub-cell [A]
  integer, dimension(n0+n1+1) :: isub_src ! Index of source cell for each sub-cell
  integer, dimension(n0) :: isrc_start ! Index of first sub-cell within each source cell
  integer, dimension(n0) :: isrc_end ! Index of last sub-cell within each source cell
  integer, dimension(n0) :: isrc_max ! Index of thickest sub-cell within each source cell
  real, dimension(n0) :: h0_eff ! Effective thickness of source cells [H]
  integer, dimension(n1) :: itgt_start ! Index of first sub-cell within each target cell
  integer, dimension(n1) :: itgt_end ! Index of last sub-cell within each target cell
  ! For error checking/debugging
  real :: u02_err ! Integrated reconstruction error estimates [H A]
  real, dimension(n0,2)           :: ppoly_r_E     ! Edge value of polynomial [A]
  real, dimension(n0,2)           :: ppoly_r_S     ! Edge slope of polynomial [A H-1]
  real, dimension(n0,CS%degree+1) :: ppoly_r_coefs ! Coefficients of polynomial reconstructions [A]
  real, dimension(n1) :: h1 !< Cell widths on target grid [H]
  real :: uh_err       ! A bound on the error in the sum of u*h, as estimated by the remapping code [A H]
  integer :: iMethod   ! An integer indicating the integration method used
  integer :: k

  call build_reconstructions_1d( CS, n0, h0, u0, ppoly_r_coefs, ppoly_r_E, ppoly_r_S, iMethod,&
                                 CS%h_neglect, CS%h_neglect_edge )

  if (CS%check_reconstruction) call check_reconstructions_1d(n0, h0, u0, CS%degree, &
                                   CS%boundary_extrapolation, ppoly_r_coefs, ppoly_r_E)

  ! This is a temporary step prior to switching to remapping_core_h()
  do k = 1, n1
    if (k<=n0) then
      h1(k) = max( 0., h0(k) + ( dx(k+1) - dx(k) ) )
    else
      h1(k) = max( 0., dx(k+1) - dx(k) )
    endif
  enddo

  ! Calculate sub-layer thicknesses and indices connecting sub-layers to source and target grids
  call intersect_src_tgt_grids( n0, h0, n1, h1, h_sub, h0_eff, &
                                isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )

  ! Loop over each sub-cell to calculate average/integral values within each sub-cell.
  ! Uses: h_sub, h0_eff, isub_src
  ! Sets: u_sub, uh_sub
  call remap_src_to_sub_grid_om4(n0, h0, u0, ppoly_r_E, ppoly_r_coefs, n1, h_sub, &
                             h0_eff, isrc_start, isrc_end, isrc_max, isub_src, &
                             iMethod, CS%force_bounds_in_subcell, u_sub, uh_sub, u02_err)

  ! Loop over each target cell summing the integrals from sub-cells within the target cell.
  ! Uses: itgt_start, itgt_end, h1, h_sub, uh_sub, u_sub
  ! Sets: u1, uh_err
  call remap_sub_to_tgt_grid_om4(n0, n1, h1, h_sub, u_sub, uh_sub, itgt_start, itgt_end, &
                             CS%force_bounds_in_target, u1, uh_err)

  ! Include the error remapping from source to sub-cells in the estimate of total remapping error
  uh_err = uh_err + u02_err

  if (CS%check_remapping) call check_remapped_values(n0, h0, u0, ppoly_r_E, CS%degree, ppoly_r_coefs, &
                                                     n1, h1, u1, iMethod, uh_err, "remapping_core_w")

end subroutine remapping_core_w

!> Creates polynomial reconstructions of u0 on the source grid h0.
subroutine build_reconstructions_1d( CS, n0, h0, u0, ppoly_r_coefs, &
                                     ppoly_r_E, ppoly_r_S, iMethod, h_neglect, &
                                     h_neglect_edge, PCM_cell, debug )
  type(remapping_CS),    intent(in)  :: CS !< Remapping control structure
  integer,               intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),   intent(in)  :: u0 !< Cell averages on source grid [A]
  real, dimension(n0,CS%degree+1), &
                         intent(out) :: ppoly_r_coefs !< Coefficients of polynomial [A]
  real, dimension(n0,2), intent(out) :: ppoly_r_E !< Edge value of polynomial [A]
  real, dimension(n0,2), intent(out) :: ppoly_r_S !< Edge slope of polynomial [A H-1]
  integer,               intent(out) :: iMethod !< Integration method
  real,                  intent(in)  :: h_neglect !< A negligibly small width for the
                                         !! purpose of cell reconstructions
                                         !! in the same units as h0 [H]
  real, optional,        intent(in)  :: h_neglect_edge !< A negligibly small width for the purpose
                                         !! of edge value calculations in the same units as h0 [H].
                                         !! The default is h_neglect.
  logical, optional,     intent(in)  :: PCM_cell(n0) !< If present, use PCM remapping for
                                         !! cells from the source grid where this is true.
  logical, optional,     intent(in) :: debug !< If true, enable debugging

  ! Local variables
  real :: h_neg_edge  ! A negligibly small width for the purpose of edge value
                      ! calculations in the same units as h0 [H]
  integer :: local_remapping_scheme
  integer :: k, n
  logical :: deb ! Do debugging

  deb=.false.; if (present(debug)) deb=debug

  h_neg_edge = h_neglect ; if (present(h_neglect_edge)) h_neg_edge = h_neglect_edge

  ! Reset polynomial
  ppoly_r_E(:,:) = 0.0
  ppoly_r_S(:,:) = 0.0
  ppoly_r_coefs(:,:) = 0.0
  iMethod = -999

  local_remapping_scheme = CS%remapping_scheme
  if (n0<=1) then
    local_remapping_scheme = REMAPPING_PCM
  elseif (n0<=3) then
    local_remapping_scheme = min( local_remapping_scheme, REMAPPING_PLM )
  elseif (n0<=4 .and. local_remapping_scheme /= REMAPPING_PPM_CW ) then
    local_remapping_scheme = min( local_remapping_scheme, REMAPPING_PPM_H4 )
  endif
  select case ( local_remapping_scheme )
    case ( REMAPPING_PCM )
      call PCM_reconstruction( n0, u0, ppoly_r_E, ppoly_r_coefs)
      iMethod = INTEGRATION_PCM
    case ( REMAPPING_PLM )
      call PLM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect )
      if ( CS%boundary_extrapolation ) then
        call PLM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect)
      endif
      iMethod = INTEGRATION_PLM
    case ( REMAPPING_PLM_HYBGEN )
      call hybgen_PLM_coefs(u0, h0, ppoly_r_coefs(:,2), n0, 1, h_neglect)
      do k=1,n0
        ppoly_r_E(k,1) = u0(k) - 0.5 * ppoly_r_coefs(k,2) ! Left edge value of cell k
        ppoly_r_E(k,2) = u0(k) + 0.5 * ppoly_r_coefs(k,2) ! Right edge value of cell k
        ppoly_r_coefs(k,1) = ppoly_r_E(k,1)
      enddo
      if ( CS%boundary_extrapolation ) &
        call PLM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect )
      iMethod = INTEGRATION_PLM
    case ( REMAPPING_PPM_CW )
      ! identical to REMAPPING_PPM_HYBGEN
      call edge_values_explicit_h4cw( n0, h0, u0, ppoly_r_E, h_neg_edge )
      call PPM_monotonicity(   n0,     u0, ppoly_r_E )
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect, answer_date=CS%answer_date )
      if ( CS%boundary_extrapolation ) then
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect )
      endif
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PPM_H4 )
      call edge_values_explicit_h4( n0, h0, u0, ppoly_r_E, h_neg_edge, answer_date=CS%answer_date )
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect, answer_date=CS%answer_date )
      if ( CS%boundary_extrapolation ) then
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect )
      endif
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PPM_IH4 )
      call edge_values_implicit_h4( n0, h0, u0, ppoly_r_E, h_neg_edge, answer_date=CS%answer_date )
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect, answer_date=CS%answer_date )
      if ( CS%boundary_extrapolation ) then
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect )
      endif
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PPM_HYBGEN )
      call hybgen_PPM_coefs(u0, h0, ppoly_r_E, n0, 1, h_neglect)
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect, answer_date=99991231 )
      if ( CS%boundary_extrapolation ) &
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect )
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_WENO_HYBGEN )
      call hybgen_weno_coefs(u0, h0, ppoly_r_E, n0, 1, h_neglect)
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect, answer_date=99991231 )
      if ( CS%boundary_extrapolation ) &
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect )
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PQM_IH4IH3 )
      call edge_values_implicit_h4( n0, h0, u0, ppoly_r_E, h_neg_edge, answer_date=CS%answer_date )
      call edge_slopes_implicit_h3( n0, h0, u0, ppoly_r_S, h_neglect, answer_date=CS%answer_date )
      call PQM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefs, h_neglect, &
                               answer_date=CS%answer_date )
      if ( CS%boundary_extrapolation ) then
        call PQM_boundary_extrapolation_v1( n0, h0, u0, ppoly_r_E, ppoly_r_S, &
                                            ppoly_r_coefs, h_neglect )
      endif
      iMethod = INTEGRATION_PQM
    case ( REMAPPING_PQM_IH6IH5 )
      call edge_values_implicit_h6( n0, h0, u0, ppoly_r_E, h_neg_edge, answer_date=CS%answer_date )
      call edge_slopes_implicit_h5( n0, h0, u0, ppoly_r_S, h_neglect, answer_date=CS%answer_date )
      call PQM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefs, h_neglect, &
                               answer_date=CS%answer_date )
      if ( CS%boundary_extrapolation ) then
        call PQM_boundary_extrapolation_v1( n0, h0, u0, ppoly_r_E, ppoly_r_S, &
                                            ppoly_r_coefs, h_neglect )
      endif
      iMethod = INTEGRATION_PQM
    case ( REMAPPING_VIA_CLASS )
      call MOM_error( FATAL, 'MOM_remapping, build_reconstructions_1d: '//&
           'Should not reach this point if using Recon1d class for remapping' )
    case default
      call MOM_error( FATAL, 'MOM_remapping, build_reconstructions_1d: '//&
           'The selected remapping method is invalid' )
  end select

  if (present(PCM_cell)) then
    ! Change the coefficients to those for the piecewise constant method in indicated cells.
    do k=1,n0 ; if (PCM_cell(k)) then
      ppoly_r_coefs(k,1) = u0(k)
      ppoly_r_E(k,1:2) = u0(k)
      ppoly_r_S(k,1:2) = 0.0
      do n=2,CS%degree+1 ; ppoly_r_coefs(k,n) = 0.0 ; enddo
    endif ; enddo
  endif

end subroutine build_reconstructions_1d

!> Checks that edge values and reconstructions satisfy bounds
subroutine check_reconstructions_1d(n0, h0, u0, deg, boundary_extrapolation, &
                                    ppoly_r_coefs, ppoly_r_E)
  integer,                  intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),      intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),      intent(in)  :: u0 !< Cell averages on source grid [A]
  integer,                  intent(in)  :: deg !< Degree of polynomial reconstruction
  logical,                  intent(in)  :: boundary_extrapolation !< Extrapolate at boundaries if true
  real, dimension(n0,deg+1),intent(in)  :: ppoly_r_coefs !< Coefficients of polynomial [A]
  real, dimension(n0,2),    intent(in)  :: ppoly_r_E !< Edge value of polynomial [A]
  ! Local variables
  integer :: i0, n
  real :: u_l, u_c, u_r ! Cell averages [A]
  real :: u_min, u_max  ! Cell extrema [A]
  logical :: problem_detected

  problem_detected = .false.
  do i0 = 1, n0
    u_l = u0(max(1,i0-1))
    u_c = u0(i0)
    u_r = u0(min(n0,i0+1))
    if (i0 > 1 .or. .not. boundary_extrapolation) then
      u_min = min(u_l, u_c)
      u_max = max(u_l, u_c)
      if (ppoly_r_E(i0,1) < u_min) then
        write(0,'(a,i4,5(1x,a,1pe24.16))') 'Left edge undershoot at',i0,'u(i0-1)=',u_l,'u(i0)=',u_c, &
                                           'edge=',ppoly_r_E(i0,1),'err=',ppoly_r_E(i0,1)-u_min
        problem_detected = .true.
      endif
      if (ppoly_r_E(i0,1) > u_max) then
        write(0,'(a,i4,5(1x,a,1pe24.16))') 'Left edge overshoot at',i0,'u(i0-1)=',u_l,'u(i0)=',u_c, &
                                           'edge=',ppoly_r_E(i0,1),'err=',ppoly_r_E(i0,1)-u_max
        problem_detected = .true.
      endif
    endif
    if (i0 < n0 .or. .not. boundary_extrapolation) then
      u_min = min(u_c, u_r)
      u_max = max(u_c, u_r)
      if (ppoly_r_E(i0,2) < u_min) then
        write(0,'(a,i4,5(1x,a,1pe24.16))') 'Right edge undershoot at',i0,'u(i0)=',u_c,'u(i0+1)=',u_r, &
                                           'edge=',ppoly_r_E(i0,2),'err=',ppoly_r_E(i0,2)-u_min
        problem_detected = .true.
      endif
      if (ppoly_r_E(i0,2) > u_max) then
        write(0,'(a,i4,5(1x,a,1pe24.16))') 'Right edge overshoot at',i0,'u(i0)=',u_c,'u(i0+1)=',u_r, &
                                           'edge=',ppoly_r_E(i0,2),'err=',ppoly_r_E(i0,2)-u_max
        problem_detected = .true.
      endif
    endif
    if (i0 > 1) then
      if ( (u_c-u_l)*(ppoly_r_E(i0,1)-ppoly_r_E(i0-1,2)) < 0.) then
        write(0,'(a,i4,5(1x,a,1pe24.16))') 'Non-monotonic edges at',i0,'u(i0-1)=',u_l,'u(i0)=',u_c, &
                                           'right edge=',ppoly_r_E(i0-1,2),'left edge=',ppoly_r_E(i0,1)
        write(0,'(5(a,1pe24.16,1x))') 'u(i0)-u(i0-1)',u_c-u_l,'edge diff=',ppoly_r_E(i0,1)-ppoly_r_E(i0-1,2)
        problem_detected = .true.
      endif
    endif
    if (problem_detected) then
      write(0,'(a,1p9e24.16)') 'Polynomial coeffs:',ppoly_r_coefs(i0,:)
      write(0,'(3(a,1pe24.16,1x))') 'u_l=',u_l,'u_c=',u_c,'u_r=',u_r
      write(0,'(a4,10a24)') 'i0','h0(i0)','u0(i0)','left edge','right edge','Polynomial coefficients'
      do n = 1, n0
        write(0,'(i4,1p10e24.16)') n,h0(n),u0(n),ppoly_r_E(n,1),ppoly_r_E(n,2),ppoly_r_coefs(n,:)
      enddo
      call MOM_error(FATAL, 'MOM_remapping, check_reconstructions_1d: '// &
                   'Edge values or polynomial coefficients were inconsistent!')
    endif
  enddo

end subroutine check_reconstructions_1d

!> Returns the intersection of source and targets grids along with and auxiliary lists or indices.
!!
!! For source grid with thicknesses h0(1:n0) and target grid with thicknesses  h1(1:n1) the intersection
!! or "subgrid" has thicknesses h_sub(1:n0+n1+1).
!! h0 and h1 must have the same units. h_sub will return with the same units as h0 and h1.
!!
!! Notes on the algorithm:
!! Internally, grids are defined by the interfaces (although we describe grids via thicknesses for accuracy).
!! The intersection or union of two grids is thus defined by the super set of both lists of interfaces.
!! Because both source and target grids can contain vanished cells, we do not eliminate repeated
!! interfaces from the union.
!! That is, the total number of interfaces of the sub-cells is equal to the total numer of interfaces of
!! the source grid (n0+1) plus the total number of interfaces of the target grid (n1+1), i.e. n0+n1+2.
!! Whenever target and source interfaces align, then the retention of identical interfaces leads to a
!! vanished subcell.
!! The remapping uses a common point of reference to the left (top) so there is always a vanished subcell
!! at the left (top).
!! If the total column thicknesses are the same, then the right (bottom) interfaces are also aligned and
!! so the last subcell will also be vanished.
subroutine intersect_src_tgt_grids( n0, h0, n1, h1, h_sub, h0_eff, &
                                    isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )
  integer, intent(in)  :: n0      !< Number of cells in source grid
  real,    intent(in)  :: h0(n0)  !< Source grid widths (size n0) [H]
  integer, intent(in)  :: n1      !< Number of cells in target grid
  real,    intent(in)  :: h1(n1)  !< Target grid widths (size n1) [H]
  real,    intent(out) :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  real,    intent(out) :: h0_eff(n0) !< Effective thickness of source cells [H]
  integer, intent(out) :: isrc_start(n0) !< Index of first sub-cell within each source cell
  integer, intent(out) :: isrc_end(n0) !< Index of last sub-cell within each source cell
  integer, intent(out) :: isrc_max(n0) !< Index of thickest sub-cell within each source cell
  integer, intent(out) :: itgt_start(n1) !< Index of first sub-cell within each target cell
  integer, intent(out) :: itgt_end(n1) !< Index of last sub-cell within each target cell
  integer, intent(out) :: isub_src(n0+n1+1) !< Index of source cell for each sub-cell
  ! Local variables
  integer :: i_sub ! Index of sub-cell
  integer :: i0 ! Index into h0(1:n0), source column
  integer :: i1 ! Index into h1(1:n1), target column
  integer :: i_start0 ! Used to record which sub-cells map to source cells
  integer :: i_start1 ! Used to record which sub-cells map to target cells
  integer :: i_max ! Used to record which sub-cell is the largest contribution of a source cell
  real :: dh_max ! Used to record which sub-cell is the largest contribution of a source cell [H]
  real :: h0_supply, h1_supply ! The amount of width available for constructing sub-cells [H]
  real :: dh ! The width of the sub-cell [H]
  real :: dh0_eff ! Running sum of source cell thickness [H]
  ! For error checking/debugging
  logical :: src_has_volume !< True if h0 has not been consumed
  logical :: tgt_has_volume !< True if h1 has not been consumed

  ! Initialize algorithm
  h0_supply = h0(1)
  h1_supply = h1(1)
  src_has_volume = .true.
  tgt_has_volume = .true.
  i0 = 1 ; i1 = 1
  i_start0 = 1 ; i_start1 = 1
  i_max = 1
  dh_max = 0.
  dh0_eff = 0.

  ! First sub-cell is always vanished
  h_sub(1) = 0.
  isrc_start(1) = 1
  isrc_end(1) = 1
  isrc_max(1) = 1
  isub_src(1) = 1

  ! Loop over each sub-cell to calculate intersections with source and target grids
  do i_sub = 2, n0+n1+1

    ! This is the width of the sub-cell, determined by which ever column has the least
    ! supply available to consume.
    dh = min(h0_supply, h1_supply)

    ! This is the running sum of the source cell thickness. After summing over each
    ! sub-cell, the sum of sub-cell thickness might differ from the original source
    ! cell thickness due to round off.
    dh0_eff = dh0_eff + min(dh, h0_supply)

    ! Record the source index (i0) that this sub-cell integral belongs to. This
    ! is needed to index the reconstruction coefficients for the source cell
    ! used in the integrals of the sub-cell width.
    isub_src(i_sub) = i0
    h_sub(i_sub) = dh

    ! For recording the largest sub-cell within a source cell.
    if (dh >= dh_max) then
      i_max = i_sub
      dh_max = dh
    endif

    ! Which ever column (source or target) has the least width left to consume determined
    ! the width, dh, of sub-cell i_sub in the expression for dh above.
    if (h0_supply <= h1_supply .and. src_has_volume) then
      ! h0_supply is smaller than h1_supply) so we consume h0_supply and increment the
      ! source cell index.
      h1_supply = h1_supply - dh ! Although this is a difference the result will
                                 ! be non-negative because of the conditional.
      ! Record the sub-cell start/end index that span the source cell i0.
      isrc_start(i0) = i_start0
      isrc_end(i0) = i_sub
      i_start0 = i_sub + 1
      ! Record the sub-cell that is the largest fraction of the source cell.
      isrc_max(i0) = i_max
      i_max = i_sub + 1
      dh_max = 0.
      ! Record the source cell thickness found by summing the sub-cell thicknesses.
      h0_eff(i0) = dh0_eff
      ! Move the source index.
      if (i0 < n0) then
        i0 = i0 + 1
        h0_supply = h0(i0)
        dh0_eff = 0.
      else
        h0_supply = 0.
        src_has_volume = .false.
      endif
    elseif (h0_supply >= h1_supply .and. tgt_has_volume) then
      ! h1_supply is smaller than h0_supply) so we consume h1_supply and increment the
      ! target cell index.
      h0_supply = h0_supply - dh ! Although this is a difference the result will
                                 ! be non-negative because of the conditional.
      ! Record the sub-cell start/end index that span the target cell i1.
      itgt_start(i1) = i_start1
      itgt_end(i1) = i_sub
      i_start1 = i_sub + 1
      ! Move the target index.
      if (i1 < n1) then
        i1 = i1 + 1
        h1_supply = h1(i1)
      else
        h1_supply = 0.
        tgt_has_volume = .false.
      endif
    elseif (src_has_volume) then
      ! We ran out of target volume but still have source cells to consume
      h_sub(i_sub) = h0_supply
      ! Record the sub-cell start/end index that span the source cell i0.
      isrc_start(i0) = i_start0
      isrc_end(i0) = i_sub
      i_start0 = i_sub + 1
      ! Record the sub-cell that is the largest fraction of the source cell.
      isrc_max(i0) = i_max
      i_max = i_sub + 1
      dh_max = 0.
      ! Record the source cell thickness found by summing the sub-cell thicknesses.
      h0_eff(i0) = dh0_eff
      if (i0 < n0) then
        i0 = i0 + 1
        h0_supply = h0(i0)
        dh0_eff = 0.
      else
        h0_supply = 0.
        src_has_volume = .false.
      endif
    elseif (tgt_has_volume) then
      ! We ran out of source volume but still have target cells to consume
      h_sub(i_sub) = h1_supply
      ! Record the sub-cell start/end index that span the target cell i1.
      itgt_start(i1) = i_start1
      itgt_end(i1) = i_sub
      i_start1 = i_sub + 1
      ! Move the target index.
      if (i1 < n1) then
        i1 = i1 + 1
        h1_supply = h1(i1)
      else
        h1_supply = 0.
        tgt_has_volume = .false.
      endif
    else
      stop 'intersect_src_tgt_grids: THIS SHOULD NEVER HAPPEN!'
    endif

  enddo

end subroutine intersect_src_tgt_grids

!> Adjust h_sub to ensure accurate conservation
!!
!! Loop over each source cell substituting the thickest sub-cell (within the source cell) with the
!! residual of the source cell thickness minus the sum of other sub-cells
!! aka a genius algorithm for accurate conservation when remapping from Robert Hallberg (\@Hallberg-NOAA).
!subroutine adjust_h_sub( n0, h0, n1, isrc_start, isrc_end, isrc_max, h_sub )
!  integer, intent(in)    :: n0      !< Number of cells in source grid
!  real,    intent(in)    :: h0(n0)  !< Source grid widths (size n0) [H]
!  integer, intent(in)    :: n1      !< Number of cells in target grid
!  integer, intent(in)    :: isrc_start(n0) !< Index of first sub-cell within each source cell
!  integer, intent(in)    :: isrc_end(n0) !< Index of last sub-cell within each source cell
!  integer, intent(in)    :: isrc_max(n0) !< Index of thickest sub-cell within each source cell
!  real,    intent(inout) :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
!  ! Local variables
!  integer :: i_sub ! Index of sub-cell
!  integer :: i0 ! Index into h0(1:n0), source column
!  integer :: i_max ! Used to record which sub-cell is the largest contribution of a source cell
!  real :: dh_max ! Used to record which sub-cell is the largest contribution of a source cell [H]
!  real :: dh ! The width of the sub-cell [H]
!  integer :: i0_last_thick_cell ! Last h0 cell with finite thickness
!
!  i0_last_thick_cell = 0
!  do i0 = 1, n0
!    if (h0(i0)>0.) i0_last_thick_cell = i0
!  enddo
!
!  do i0 = 1, i0_last_thick_cell
!    i_max = isrc_max(i0)
!    dh_max = h_sub(i_max)
!    if (dh_max > 0.) then
!      ! dh will be the sum of sub-cell thicknesses within the source cell except for the thickest sub-cell.
!      dh = 0.
!      do i_sub = isrc_start(i0), isrc_end(i0)
!        if (i_sub /= i_max) dh = dh + h_sub(i_sub)
!      enddo
!      h_sub(i_max) = h0(i0) - dh
!    endif
!  enddo
!
!end subroutine adjust_h_sub

!> Remaps column of n0 values u0 on grid h0 to subgrid h_sub
!!
!! This includes an error for the scenario where the source grid is much thicker than
!! the target grid and extrapolation is needed.
subroutine remap_src_to_sub_grid_om4(n0, h0, u0, ppoly0_E, ppoly0_coefs, n1, h_sub, &
                                 h0_eff, isrc_start, isrc_end, isrc_max, isub_src, &
                                 method, force_bounds_in_subcell, u_sub, uh_sub, u02_err)
  integer, intent(in)  :: n0      !< Number of cells in source grid
  real,    intent(in)  :: h0(n0)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: u0(n0)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: ppoly0_E(n0,2)    !< Edge value of polynomial [A]
  real,    intent(in)  :: ppoly0_coefs(:,:) !< Coefficients of polynomial [A]
  integer, intent(in)  :: n1      !< Number of cells in target grid
  real,    intent(in)  :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  real,    intent(in)  :: h0_eff(n0) !< Effective thickness of source cells [H]
  integer, intent(in)  :: isrc_start(n0) !< Index of first sub-cell within each source cell
  integer, intent(in)  :: isrc_end(n0) !< Index of last sub-cell within each source cell
  integer, intent(in)  :: isrc_max(n0) !< Index of thickest sub-cell within each source cell
  integer, intent(in)  :: isub_src(n0+n1+1) !< Index of source cell for each sub-cell
  integer, intent(in)  :: method  !< Remapping scheme to use
  logical, intent(in)  :: force_bounds_in_subcell !< Force sub-cell values to be bounded
  real,    intent(out) :: u_sub(n0+n1+1) !< Sub-cell cell averages (size n1) [A]
  real,    intent(out) :: uh_sub(n0+n1+1) !< Sub-cell cell integrals (size n1) [A H]
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
  real :: u0_min(n0), u0_max(n0) !< Min/max of u0 for each source cell [A]
  ! For error checking/debugging
  logical, parameter :: adjust_thickest_subcell = .true. ! To fix round-off conservation issues
  integer :: i0_last_thick_cell
  real :: u_orig              ! The original value of the reconstruction in a cell [A]

  i0_last_thick_cell = 0
  do i0 = 1, n0
    u0_min(i0) = min(ppoly0_E(i0,1), ppoly0_E(i0,2))
    u0_max(i0) = max(ppoly0_E(i0,1), ppoly0_E(i0,2))
    if (h0(i0)>0.) i0_last_thick_cell = i0
  enddo

  ! Loop over each sub-cell to calculate average/integral values within each sub-cell.
  ! Uses: h_sub, isub_src, h0_eff
  ! Sets: u_sub, uh_sub
  xa = 0.
  dh0_eff = 0.
  uh_sub(1) = 0.
  u_sub(1) = ppoly0_E(1,1)
  u02_err = 0.
  do i_sub = 2, n0+n1

    ! Sub-cell thickness from loop above
    dh = h_sub(i_sub)

    ! Source cell
    i0 = isub_src(i_sub)

    ! Evaluate average and integral for sub-cell i_sub.
    ! Integral is over distance dh but expressed in terms of non-dimensional
    ! positions with source cell from xa to xb  (0 <= xa <= xb <= 1).
    dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
    if (h0_eff(i0)>0.) then
      xb = dh0_eff / h0_eff(i0) ! This expression yields xa <= xb <= 1.0
      xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
      u_sub(i_sub) = average_value_ppoly( n0, u0, ppoly0_E, ppoly0_coefs, method, i0, xa, xb)
    else ! Vanished cell
      xb = 1.
      u_sub(i_sub) = u0(i0)
    endif
    if (force_bounds_in_subcell) then
      ! These next two lines should not be needed but when using PQM we found roundoff
      ! can lead to overshoots. These lines sweep issues under the rug which need to be
      ! properly .. later. -AJA
      u_orig = u_sub(i_sub)
      u_sub(i_sub) = max( u_sub(i_sub), u0_min(i0) )
      u_sub(i_sub) = min( u_sub(i_sub), u0_max(i0) )
      u02_err = u02_err + dh*abs( u_sub(i_sub) - u_orig )
    endif
    uh_sub(i_sub) = dh * u_sub(i_sub)

    if (isub_src(i_sub+1) /= i0) then
      ! If the next sub-cell is in a different source cell, reset the position counters
      dh0_eff = 0.
      xa = 0.
    else
      xa = xb ! Next integral will start at end of last
    endif

  enddo
  u_sub(n0+n1+1) = ppoly0_E(n0,2)                   ! This value is only needed when total target column
  uh_sub(n0+n1+1) = ppoly0_E(n0,2) * h_sub(n0+n1+1) ! is wider than the source column

  if (adjust_thickest_subcell) then
    ! Loop over each source cell substituting the integral/average for the thickest sub-cell (within
    ! the source cell) with the residual of the source cell integral minus the other sub-cell integrals
    ! aka a genius algorithm for accurate conservation when remapping from Robert Hallberg (\@Hallberg-NOAA).
    ! Uses: i0_last_thick_cell, isrc_max, h_sub, isrc_start, isrc_end, uh_sub, u0, h0
    ! Updates: uh_sub, u_sub
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
  endif

end subroutine remap_src_to_sub_grid_om4

!> Remaps column of n0 values u0 on grid h0 to subgrid h_sub
subroutine remap_src_to_sub_grid(n0, h0, u0, ppoly0_E, ppoly0_coefs, n1, h_sub, &
                                 isrc_start, isrc_end, isrc_max, isub_src, &
                                 method, force_bounds_in_subcell, u_sub, uh_sub, u02_err)
  integer, intent(in)  :: n0      !< Number of cells in source grid
  real,    intent(in)  :: h0(n0)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: u0(n0)  !< Source grid widths (size n0) [H]
  real,    intent(in)  :: ppoly0_E(n0,2)    !< Edge value of polynomial [A]
  real,    intent(in)  :: ppoly0_coefs(:,:) !< Coefficients of polynomial [A]
  integer, intent(in)  :: n1      !< Number of cells in target grid
  real,    intent(in)  :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  integer, intent(in)  :: isrc_start(n0) !< Index of first sub-cell within each source cell
  integer, intent(in)  :: isrc_end(n0) !< Index of last sub-cell within each source cell
  integer, intent(in)  :: isrc_max(n0) !< Index of thickest sub-cell within each source cell
  integer, intent(in)  :: isub_src(n0+n1+1) !< Index of source cell for each sub-cell
  integer, intent(in)  :: method  !< Remapping scheme to use
  logical, intent(in)  :: force_bounds_in_subcell !< Force sub-cell values to be bounded
  real,    intent(out) :: u_sub(n0+n1+1) !< Sub-cell cell averages (size n1) [A]
  real,    intent(out) :: uh_sub(n0+n1+1) !< Sub-cell cell integrals (size n1) [A H]
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
  real :: u0_min(n0), u0_max(n0) ! Min/max of u0 for each source cell [A]
  ! For error checking/debugging
  logical, parameter :: adjust_thickest_subcell = .true. ! To fix round-off conservation issues
  integer :: i0_last_thick_cell
  real :: u_orig              ! The original value of the reconstruction in a cell [A]

  i0_last_thick_cell = 0
  do i0 = 1, n0
    u0_min(i0) = min(ppoly0_E(i0,1), ppoly0_E(i0,2))
    u0_max(i0) = max(ppoly0_E(i0,1), ppoly0_E(i0,2))
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
      u_sub(i_sub) = average_value_ppoly( n0, u0, ppoly0_E, ppoly0_coefs, method, i0, xa, xb)
    else ! Vanished cell
      xb = 1.
      u_sub(i_sub) = u0(i0)
    endif
    if (force_bounds_in_subcell) then
      ! These next two lines should not be needed but when using PQM we found roundoff
      ! can lead to overshoots. These lines sweep issues under the rug which need to be
      ! properly .. later. -AJA
      u_orig = u_sub(i_sub)
      u_sub(i_sub) = max( u_sub(i_sub), u0_min(i0) )
      u_sub(i_sub) = min( u_sub(i_sub), u0_max(i0) )
      u02_err = u02_err + dh*abs( u_sub(i_sub) - u_orig )
    endif
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
    u_sub(i_sub) = average_value_ppoly( n0, u0, ppoly0_E, ppoly0_coefs, method, i0, xa, xb)
  else ! Vanished cell
    xb = 1.
    u_sub(i_sub) = u0(i0)
  endif
  if (force_bounds_in_subcell) then
    ! These next two lines should not be needed but when using PQM we found roundoff
    ! can lead to overshoots. These lines sweep issues under the rug which need to be
    ! properly .. later. -AJA
    u_orig = u_sub(i_sub)
    u_sub(i_sub) = max( u_sub(i_sub), u0_min(i0) )
    u_sub(i_sub) = min( u_sub(i_sub), u0_max(i0) )
    u02_err = u02_err + dh*abs( u_sub(i_sub) - u_orig )
  endif
  uh_sub(i_sub) = dh * u_sub(i_sub)

  if (adjust_thickest_subcell) then
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
  endif

end subroutine remap_src_to_sub_grid

!> Remaps column of n0+n1+1 values usub on sub-grid h_sub to targets on grid h1
!! using the OM4-era algorithm
subroutine remap_sub_to_tgt_grid_om4(n0, n1, h1, h_sub, u_sub, uh_sub, &
                                 itgt_start, itgt_end, force_bounds_in_target, u1, uh_err)
  integer, intent(in)  :: n0     !< Number of cells in source grid
  integer, intent(in)  :: n1     !< Number of cells in target grid
  real,    intent(in)  :: h1(n1) !< Target grid widths (size n1) [H]
  real,    intent(in)  :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  real,    intent(in)  :: u_sub(n0+n1+1) !< Sub-cell cell averages (size n1) [A]
  real,    intent(in)  :: uh_sub(n0+n1+1) !< Sub-cell cell integrals (size n1) [A H]
  integer, intent(in)  :: itgt_start(n1) !< Index of first sub-cell within each target cell
  integer, intent(in)  :: itgt_end(n1) !< Index of last sub-cell within each target cell
  logical, intent(in)  :: force_bounds_in_target !< Force sub-cell values to be bounded
  real,    intent(out) :: u1(n1) !< Target cell averages (size n1) [A]
  real,    intent(out) :: uh_err !< Estimate of bound on error in sum of u*h [A H]
  ! Local variables
  integer :: i1 ! tgt loop index
  integer :: i_sub ! index to sub-layer
  real :: dh ! The width of the sub-cell [H]
  real :: duh ! The total amount of accumulated stuff (u*h)  [A H]
  real :: u1min, u1max ! Minimum and maximum values of reconstructions [A]
  real :: u_orig ! The original value of the reconstruction in a cell prior to bounding [A]

  u1min = 0. ! Not necessary, but avoids an overzealous compiler ...
  u1max = 0. ! ... warning about uninitialized variables

  ! Loop over each target cell summing the integrals from sub-cells within the target cell.
  ! Uses: itgt_start, itgt_end, h_sub, uh_sub, u_sub
  ! Sets: u1, uh_err
  uh_err = 0.
  do i1 = 1, n1
    if (h1(i1) > 0.) then
      duh = 0. ; dh = 0.
      i_sub = itgt_start(i1)
      if (force_bounds_in_target) then
        u1min = u_sub(i_sub)
        u1max = u_sub(i_sub)
      endif
      do i_sub = itgt_start(i1), itgt_end(i1)
        if (force_bounds_in_target) then
          u1min = min(u1min, u_sub(i_sub))
          u1max = max(u1max, u_sub(i_sub))
        endif
        dh = dh + h_sub(i_sub)
        duh = duh + uh_sub(i_sub)
        ! This accumulates the contribution to the error bound for the sum of u*h
        uh_err = uh_err + max(abs(duh),abs(uh_sub(i_sub)))*epsilon(duh)
      enddo
      u1(i1) = duh / dh
      ! This is the contribution from the division to the error bound for the sum of u*h
      uh_err = uh_err + abs(duh)*epsilon(duh)
      if (force_bounds_in_target) then
        u_orig = u1(i1)
        u1(i1) = max(u1min, min(u1max, u1(i1)))
        ! Adjusting to be bounded contributes to the error for the sum of u*h
        uh_err = uh_err + dh*abs( u1(i1)-u_orig )
      endif
    else
      u1(i1) = u_sub(itgt_start(i1))
    endif
  enddo

end subroutine remap_sub_to_tgt_grid_om4

!> Remaps column of n0+n1+1 values usub on sub-grid h_sub to targets on grid h1
subroutine remap_sub_to_tgt_grid(n0, n1, h1, h_sub, u_sub, uh_sub, &
                                 itgt_start, itgt_end, force_bounds_in_target, &
                                 better_force_bounds_in_target, offset_summation, u1, uh_err)
  integer, intent(in)  :: n0     !< Number of cells in source grid
  integer, intent(in)  :: n1     !< Number of cells in target grid
  real,    intent(in)  :: h1(n1) !< Target grid widths (size n1) [H]
  real,    intent(in)  :: h_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  real,    intent(in)  :: u_sub(n0+n1+1) !< Sub-cell cell averages (size n1) [A]
  real,    intent(in)  :: uh_sub(n0+n1+1) !< Sub-cell cell integrals (size n1) [A H]
  integer, intent(in)  :: itgt_start(n1) !< Index of first sub-cell within each target cell
  integer, intent(in)  :: itgt_end(n1) !< Index of last sub-cell within each target cell
  logical, intent(in)  :: force_bounds_in_target !< Force sub-cell values to be bounded
  logical, intent(in)  :: better_force_bounds_in_target !< Force sub-cell values to be bounded
  logical, intent(in)  :: offset_summation !< Offset values in summation for accuracy
  real,    intent(out) :: u1(n1) !< Target cell averages (size n1) [A]
  real,    intent(out) :: uh_err !< Estimate of bound on error in sum of u*h [A H]
  ! Local variables
  integer :: i1 ! tgt loop index
  integer :: i_sub ! index to sub-layer
  real :: dh ! The width of the sub-cell [H]
  real :: duh ! The total amount of accumulated stuff (u*h)  [A H]
  real :: u1min, u1max ! Minimum and maximum values of reconstructions [A]
  real :: u_orig ! The original value of the reconstruction in a cell prior to bounding [A]
  real :: u_ref ! A value to offest the summation to gain accuracy [A]
  real :: h_max ! Thickest cell encountered [H]

  u1min = 0. ! Not necessary, but avoids an overzealous compiler ...
  u1max = 0. ! ... warning about uninitialized variables
  u_ref = 0. ! An offset of 0. should do no harm
  h_max = 0.

  ! Loop over each target cell summing the integrals from sub-cells within the target cell.
  ! Uses: itgt_start, itgt_end, h_sub, uh_sub, u_sub
  ! Sets: u1, uh_err
  uh_err = 0.
  do i1 = 1, n1
    if (h1(i1) > 0.) then
      duh = 0. ; dh = 0.
      i_sub = itgt_start(i1)
      if (force_bounds_in_target) then
        u1min = u_sub(i_sub)
        u1max = u_sub(i_sub)
      endif
      if (offset_summation) then
        u_ref = 0. ! An offset of 0. should do no harm
        h_max = 0.
        do i_sub = itgt_start(i1), itgt_end(i1)
          if (h_sub(i_sub) > h_max) then
            u_ref = u_sub(i_sub)
            h_max = h_sub(i_sub)
          endif
        enddo
      endif
      do i_sub = itgt_start(i1), itgt_end(i1)
        if (force_bounds_in_target .or. better_force_bounds_in_target .and. h_sub(i_sub)>0.) then
          u1min = min(u1min, u_sub(i_sub))
          u1max = max(u1max, u_sub(i_sub))
        endif
        dh = dh + h_sub(i_sub)
        ! Ideally u_ref would be already be substracted in uh_sub
        duh = duh + ( uh_sub(i_sub) - h_sub(i_sub) * u_ref )
        ! This accumulates the contribution to the error bound for the sum of u*h
        uh_err = uh_err + max(abs(duh),abs(uh_sub(i_sub)))*epsilon(duh)
      enddo
      u1(i1) = duh / dh + u_ref
      ! This is the contribution from the division to the error bound for the sum of u*h
      uh_err = uh_err + abs(duh)*epsilon(duh)
      if (force_bounds_in_target) then
        u_orig = u1(i1)
        u1(i1) = max(u1min, min(u1max, u1(i1)))
        ! Adjusting to be bounded contributes to the error for the sum of u*h
        uh_err = uh_err + dh*abs( u1(i1)-u_orig )
      endif
    else
      u1(i1) = u_sub(itgt_start(i1))
    endif
  enddo

end subroutine remap_sub_to_tgt_grid

!> Linearly interpolate interface data, u_src, from grid h_src to a grid h_dest
subroutine interpolate_column(nsrc, h_src, u_src, ndest, h_dest, u_dest, mask_edges)
  integer,                  intent(in)    :: nsrc   !< Number of source cells
  real, dimension(nsrc),    intent(in)    :: h_src  !< Thickness of source cells [H]
  real, dimension(nsrc+1),  intent(in)    :: u_src  !< Values at source cell interfaces [A]
  integer,                  intent(in)    :: ndest  !< Number of destination cells
  real, dimension(ndest),   intent(in)    :: h_dest !< Thickness of destination cells [H]
  real, dimension(ndest+1), intent(inout) :: u_dest !< Interpolated value at destination cell interfaces [A]
  logical,                  intent(in)    :: mask_edges !< If true, mask the values outside of massless
                                                    !! layers at the top and bottom of the column.

  ! Local variables
  real :: x_dest            ! Relative position of target interface [H]
  real :: dh                ! Source cell thickness [H]
  real :: frac_pos(ndest+1) ! Fractional position of the destination interface
                            ! within the source layer [nondim], 0 <= frac_pos <= 1.
  integer :: k_src(ndest+1) ! Source grid layer index of destination interface, 1 <= k_src <= ndest.
  integer :: ks, k_dest     ! Index of cell in src and dest columns

  ! The following forces the "do while" loop to do one cycle that will set u1, u2, dh.
  ks = 0
  dh = 0.
  x_dest = 0.

  ! Find the layer index and fractional position of the interfaces of the target
  ! grid on the source grid.
  do k_dest=1,ndest+1
    do while (dh<=x_dest .and. ks<nsrc)
      ! Move positions pointers forward until the interval 0 .. dh spans x_dest.
      x_dest = x_dest - dh
      ks = ks + 1
      dh = h_src(ks) ! Thickness of layer ks
    enddo
    k_src(k_dest) = ks

    if (dh>0.) then
      frac_pos(k_dest) = max(0., min(1., x_dest / dh)) ! Weight of u2
    else  ! For a vanished source layer we need to do something reasonable...
      frac_pos(k_dest) = 0.5
    endif

    if (k_dest <= ndest) then
      x_dest = x_dest + h_dest(k_dest) ! Position of interface k_dest+1
    endif
  enddo

  do k_dest=1,ndest+1
    ! Linear interpolation between surrounding edge values.
    ks = k_src(k_dest)
    u_dest(k_dest) = (1.0 - frac_pos(k_dest)) * u_src(ks) + frac_pos(k_dest) * u_src(ks+1)
  enddo

  if (mask_edges) then
    ! Mask vanished layers at the surface which would be under an ice-shelf.
    ! When the layer k_dest is vanished and all layers above are also vanished,
    ! the k_dest interface value should be missing.
    do k_dest=1,ndest
      if (h_dest(k_dest) > 0.) exit
      u_dest(k_dest) = 0.0
    enddo

    ! Mask interfaces below vanished layers at the bottom
    do k_dest=ndest,1,-1
      if (h_dest(k_dest) > 0.) exit
      u_dest(k_dest+1) = 0.0
    enddo
  endif

end subroutine interpolate_column

!> Conservatively calculate integrated data, uh_dest, on grid h_dest, from layer-integrated data, uh_src, on grid h_src
subroutine reintegrate_column(nsrc, h_src, uh_src, ndest, h_dest, uh_dest)
  integer,                intent(in)    :: nsrc    !< Number of source cells
  real, dimension(nsrc),  intent(in)    :: h_src   !< Thickness of source cells [H]
  real, dimension(nsrc),  intent(in)    :: uh_src  !< Values at source cell interfaces [A H]
  integer,                intent(in)    :: ndest   !< Number of destination cells
  real, dimension(ndest), intent(in)    :: h_dest  !< Thickness of destination cells [H]
  real, dimension(ndest), intent(inout) :: uh_dest !< Interpolated value at destination cell interfaces [A H]

  ! Local variables
  real :: h_src_rem, h_dest_rem, dh ! Incremental thicknesses [H]
  real :: uh_src_rem, duh  ! Incremental amounts of stuff [A H]
  integer :: k_src, k_dest ! Index of cell in src and dest columns
  logical :: src_ran_out

  uh_dest(:) = 0.0

  k_src = 0
  k_dest = 0
  h_dest_rem = 0.
  h_src_rem = 0.
  uh_src_rem = 0.
  src_ran_out = .false.

  do while(.true.)
    if (h_src_rem==0. .and. k_src<nsrc) then
      ! Supply is empty so move to the next source cell
      k_src = k_src + 1
      h_src_rem = h_src(k_src)
      uh_src_rem = uh_src(k_src)
      if (h_src_rem==0.) cycle
    endif
    if (h_dest_rem==0. .and. k_dest<ndest) then
      ! Sink has no capacity so move to the next destination cell
      k_dest = k_dest + 1
      h_dest_rem = h_dest(k_dest)
      uh_dest(k_dest) = 0.
      if (h_dest_rem==0.) cycle
    endif
    if (k_src==nsrc .and. h_src_rem==0.) then
      if (src_ran_out) exit ! This is the second time implying there is no more src
      src_ran_out = .true.
      cycle
    endif
    duh = 0.
    if (h_src_rem<h_dest_rem) then
      ! The source cell is fully within the destination cell
      dh = h_src_rem
      if (dh>0.) duh = uh_src_rem
      h_src_rem = 0.
      uh_src_rem = 0.
      h_dest_rem = max(0., h_dest_rem - dh)
    elseif (h_src_rem>h_dest_rem) then
      ! Only part of the source cell can be used up
      dh = h_dest_rem
      duh = (dh / h_src_rem) * uh_src_rem
      h_src_rem = max(0., h_src_rem - dh)
      uh_src_rem = uh_src_rem - duh
      h_dest_rem = 0.
    else ! h_src_rem==h_dest_rem
      ! The source cell exactly fits the destination cell
      duh = uh_src_rem
      h_src_rem = 0.
      uh_src_rem = 0.
      h_dest_rem = 0.
    endif
    uh_dest(k_dest) = uh_dest(k_dest) + duh
    if (k_dest==ndest .and. (k_src==nsrc .or. h_dest_rem==0.)) exit
  enddo

end subroutine reintegrate_column

!> Returns the average value of a reconstruction within a single source cell, i0,
!! between the non-dimensional positions xa and xb (xa<=xb) with dimensional
!! separation dh.
real function average_value_ppoly( n0, u0, ppoly0_E, ppoly0_coefs, method, i0, xa, xb)
  integer,       intent(in)    :: n0     !< Number of cells in source grid
  real,          intent(in)    :: u0(n0) !< Cell means [A]
  real,          intent(in)    :: ppoly0_E(n0,2)    !< Edge value of polynomial [A]
  real,          intent(in)    :: ppoly0_coefs(:,:) !< Coefficients of polynomial [A]
  integer,       intent(in)    :: method !< Remapping scheme to use
  integer,       intent(in)    :: i0     !< Source cell index
  real,          intent(in)    :: xa     !< Non-dimensional start position within source cell [nondim]
  real,          intent(in)    :: xb     !< Non-dimensional end position within source cell [nondim]
  ! Local variables
  real :: u_ave                 ! The average value of the polynomial over the specified range [A]
  real :: xapxb                 ! A sum of fracional positions [nondim]
  real :: mx, Ya, Yb, my        ! Various fractional positions [nondim]
  real :: xa_2, xb_2            ! Squared fractional positions [nondim]
  real :: xa2pxb2,  xa2b2ab, Ya2b2ab  ! Sums of squared fractional positions [nondim]
  real :: a_L, a_R, u_c, a_c    ! Values of the polynomial at various locations [A]
  real, parameter :: r_3 = 1.0/3.0 ! Used in evaluation of integrated polynomials [nondim]

  u_ave = 0. ! Avoids warnings about "potentially unset values"; u_ave is always calculated for legitimate schemes
  if (xb > xa) then
    select case ( method )
      case ( INTEGRATION_PCM )
        u_ave = u0(i0)
      case ( INTEGRATION_PLM )
        u_ave = (                                           &
            ppoly0_coefs(i0,1)                       &
          + ppoly0_coefs(i0,2) * 0.5 * ( xb + xa ) )
      case ( INTEGRATION_PPM )
        mx = 0.5 * ( xa + xb )
        a_L = ppoly0_E(i0, 1)
        a_R = ppoly0_E(i0, 2)
        u_c = u0(i0)
        a_c = 0.5 * ( ( u_c - a_L ) + ( u_c - a_R ) ) ! a_6 / 6
        if (mx<0.5) then
          ! This integration of the PPM reconstruction is expressed in distances from the left edge
          xa2b2ab = (xa*xa+xb*xb)+xa*xb
          u_ave = a_L + ( ( a_R - a_L ) * mx &
                          + a_c * ( 3. * ( xb + xa ) - 2.*xa2b2ab ) )
        else
          ! This integration of the PPM reconstruction is expressed in distances from the right edge
          Ya = 1. - xa
          Yb = 1. - xb
          my = 0.5 * ( Ya + Yb )
          Ya2b2ab = (Ya*Ya+Yb*Yb)+Ya*Yb
          u_ave = a_R  + ( ( a_L - a_R ) * my &
                           + a_c * ( 3. * ( Yb + Ya ) - 2.*Ya2b2ab ) )
        endif
      case ( INTEGRATION_PQM )
        xa_2 = xa*xa
        xb_2 = xb*xb
        xa2pxb2 = xa_2 + xb_2
        xapxb = xa + xb
        u_ave = (                                                                               &
              ppoly0_coefs(i0,1)                                                         &
          + ( ppoly0_coefs(i0,2) * 0.5 * ( xapxb )                                       &
          + ( ppoly0_coefs(i0,3) * r_3 * ( xa2pxb2 + xa*xb )                             &
          + ( ppoly0_coefs(i0,4) * 0.25* ( xa2pxb2 * xapxb )                             &
          +   ppoly0_coefs(i0,5) * 0.2 * ( ( xb*xb_2 + xa*xa_2 ) * xapxb + xa_2*xb_2 ) ) ) ) )
      case default
        call MOM_error( FATAL,'The selected integration method is invalid' )
    end select
  else ! dh == 0.
    select case ( method )
      case ( INTEGRATION_PCM )
        u_ave =        ppoly0_coefs(i0,1)
      case ( INTEGRATION_PLM )
       !u_ave =        ppoly0_coefs(i0,1)   &
       !      + xa *   ppoly0_coefs(i0,2)
        a_L = ppoly0_E(i0, 1)
        a_R = ppoly0_E(i0, 2)
        Ya = 1. - xa
        if (xa < 0.5) then
          u_ave = a_L + xa * ( a_R - a_L )
        else
          u_ave = a_R + Ya * ( a_L - a_R )
        endif
      case ( INTEGRATION_PPM )
       !u_ave =        ppoly0_coefs(i0,1)   &
       !      + xa * ( ppoly0_coefs(i0,2)   &
       !      + xa *   ppoly0_coefs(i0,3) )
        a_L = ppoly0_E(i0, 1)
        a_R = ppoly0_E(i0, 2)
        u_c = u0(i0)
        a_c = 3. * ( ( u_c - a_L ) + ( u_c - a_R ) ) ! a_6
        Ya = 1. - xa
        if (xa < 0.5) then
          u_ave = a_L + xa * ( ( a_R - a_L ) + a_c * Ya )
        else
          u_ave = a_R + Ya * ( ( a_L - a_R ) + a_c * xa )
        endif
      case ( INTEGRATION_PQM )
        u_ave =        ppoly0_coefs(i0,1)   &
              + xa * ( ppoly0_coefs(i0,2)   &
              + xa * ( ppoly0_coefs(i0,3)   &
              + xa * ( ppoly0_coefs(i0,4)   &
              + xa *   ppoly0_coefs(i0,5) ) ) )
      case default
        u_ave = 0.
        call MOM_error( FATAL,'The selected integration method is invalid' )
    end select
  endif
  average_value_ppoly = u_ave

end function average_value_ppoly

!> This subroutine checks for sufficient consistence in the extrema and total amounts on the old
!! and new grids.
subroutine check_remapped_values(n0, h0, u0, ppoly_r_E, deg, ppoly_r_coefs, &
                                 n1, h1, u1, iMethod, uh_err, caller)
  integer,               intent(in) :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in) :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),   intent(in) :: u0 !< Cell averages on source grid [A]
  real, dimension(n0,2), intent(in) :: ppoly_r_E  !< Edge values of polynomial fits [A]
  integer,               intent(in) :: deg !< Degree of the piecewise polynomial reconstrution
  real, dimension(n0,deg+1), intent(in) :: ppoly_r_coefs !< Coefficients of the piecewise
                                          !! polynomial reconstructions [A]
  integer,               intent(in) :: n1 !< Number of cells on target grid
  real, dimension(n1),   intent(in) :: h1 !< Cell widths on target grid [H]
  real, dimension(n1),   intent(in) :: u1 !< Cell averages on target grid [A]
  integer,               intent(in) :: iMethod !< An integer indicating the integration method used
  real,                  intent(in) :: uh_err  !< A bound on the error in the sum of u*h as
                                               !! estimated by the remapping code [H A]
  character(len=*),      intent(in) :: caller  !< The name of the calling routine.

  ! Local variables
  real :: h0tot, h0err ! Sum of source cell widths and round-off error in this sum [H]
  real :: h1tot, h1err ! Sum of target cell widths and round-off error in this sum [H]
  real :: u0tot, u0err ! Integrated values on the source grid and round-off error in this sum [H A]
  real :: u1tot, u1err ! Integrated values on the target grid and round-off error in this sum [H A]
  real :: u0min, u0max, u1min, u1max ! Extrema of values on the two grids [A]
  integer :: k

  ! Check errors and bounds
  call measure_input_bounds( n0, h0, u0, ppoly_r_E, h0tot, h0err, u0tot, u0err, u0min, u0max )
  call measure_output_bounds( n1, h1, u1, h1tot, h1err, u1tot, u1err, u1min, u1max )

  if (iMethod<5) return ! We except PQM until we've debugged it

  if ( (abs(u1tot-u0tot)>(u0err+u1err)+uh_err .and. abs(h1tot-h0tot)<h0err+h1err) &
      .or. (u1min<u0min .or. u1max>u0max) ) then
    write(0,*) 'iMethod = ',iMethod
    write(0,*) 'H: h0tot=',h0tot,'h1tot=',h1tot,'dh=',h1tot-h0tot,'h0err=',h0err,'h1err=',h1err
    if (abs(h1tot-h0tot)>h0err+h1err) &
      write(0,*) 'H non-conservation difference=',h1tot-h0tot,'allowed err=',h0err+h1err,' <-----!'
    write(0,*) 'UH: u0tot=',u0tot,'u1tot=',u1tot,'duh=',u1tot-u0tot,'u0err=',u0err,'u1err=',u1err,'uh_err=',uh_err
    if (abs(u1tot-u0tot)>(u0err+u1err)+uh_err) &
      write(0,*) 'U non-conservation difference=',u1tot-u0tot,'allowed err=',u0err+u1err+uh_err,' <-----!'
    write(0,*) 'U: u0min=',u0min,'u1min=',u1min
    if (u1min<u0min) write(0,*) 'U minimum overshoot=',u1min-u0min,' <-----!'
    write(0,*) 'U: u0max=',u0max,'u1max=',u1max
    if (u1max>u0max) write(0,*) 'U maximum overshoot=',u1max-u0max,' <-----!'
    write(0,'(a3,6a24)') 'k','h0','left edge','u0','right edge','h1','u1'
    do k = 1, max(n0,n1)
      if (k<=min(n0,n1)) then
        write(0,'(i3,1p6e24.16)') k,h0(k),ppoly_r_E(k,1),u0(k),ppoly_r_E(k,2),h1(k),u1(k)
      elseif (k>n0) then
        write(0,'(i3,96x,1p2e24.16)') k,h1(k),u1(k)
      else
        write(0,'(i3,1p4e24.16)') k,h0(k),ppoly_r_E(k,1),u0(k),ppoly_r_E(k,2)
      endif
    enddo
    write(0,'(a3,2a24)') 'k','u0','Polynomial coefficients'
    do k = 1, n0
      write(0,'(i3,1p6e24.16)') k,u0(k),ppoly_r_coefs(k,:)
    enddo
    call MOM_error( FATAL, 'MOM_remapping, '//trim(caller)//': '//&
           'Remapping result is inconsistent!' )
  endif

end subroutine check_remapped_values

!> Measure totals and bounds on source grid
subroutine measure_input_bounds( n0, h0, u0, edge_values, h0tot, h0err, u0tot, u0err, u0min, u0max )
  integer,               intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),   intent(in)  :: u0 !< Cell averages on source grid [A]
  real, dimension(n0,2), intent(in)  :: edge_values !< Cell edge values on source grid [A]
  real,                  intent(out) :: h0tot !< Sum of cell widths [H]
  real,                  intent(out) :: h0err !< Magnitude of round-off error in h0tot [H]
  real,                  intent(out) :: u0tot !< Sum of cell widths times values [H A]
  real,                  intent(out) :: u0err !< Magnitude of round-off error in u0tot [H A]
  real,                  intent(out) :: u0min !< Minimum value in reconstructions of u0 [A]
  real,                  intent(out) :: u0max !< Maximum value in reconstructions of u0 [A]
  ! Local variables
  real :: eps  ! The smallest representable fraction of a number [nondim]
  integer :: k

  eps = epsilon(h0(1))
  h0tot = h0(1)
  h0err = 0.
  u0tot = h0(1) * u0(1)
  u0err = 0.
  u0min = min( edge_values(1,1), edge_values(1,2) )
  u0max = max( edge_values(1,1), edge_values(1,2) )
  do k = 2, n0
    h0tot = h0tot + h0(k)
    h0err = h0err + eps * max(h0tot, h0(k))
    u0tot = u0tot + h0(k) * u0(k)
    u0err = u0err + eps * max(abs(u0tot), abs(h0(k) * u0(k)))
    u0min = min( u0min, edge_values(k,1), edge_values(k,2) )
    u0max = max( u0max, edge_values(k,1), edge_values(k,2) )
  enddo

end subroutine measure_input_bounds

!> Measure totals and bounds on destination grid
subroutine measure_output_bounds( n1, h1, u1, h1tot, h1err, u1tot, u1err, u1min, u1max )
  integer,               intent(in)  :: n1 !< Number of cells on destination grid
  real, dimension(n1),   intent(in)  :: h1 !< Cell widths on destination grid [H]
  real, dimension(n1),   intent(in)  :: u1 !< Cell averages on destination grid [A]
  real,                  intent(out) :: h1tot !< Sum of cell widths [H]
  real,                  intent(out) :: h1err !< Magnitude of round-off error in h1tot [H]
  real,                  intent(out) :: u1tot !< Sum of cell widths times values [H A]
  real,                  intent(out) :: u1err !< Magnitude of round-off error in u1tot [H A]
  real,                  intent(out) :: u1min !< Minimum value in reconstructions of u1 [A]
  real,                  intent(out) :: u1max !< Maximum value in reconstructions of u1 [A]
  ! Local variables
  real :: eps  ! The smallest representable fraction of a number [nondim]
  integer :: k

  eps = epsilon(h1(1))
  h1tot = h1(1)
  h1err = 0.
  u1tot = h1(1) * u1(1)
  u1err = 0.
  u1min = u1(1)
  u1max = u1(1)
  do k = 2, n1
    h1tot = h1tot + h1(k)
    h1err = h1err + eps * max(h1tot, h1(k))
    u1tot = u1tot + h1(k) * u1(k)
    u1err = u1err + eps * max(abs(u1tot), abs(h1(k) * u1(k)))
    u1min = min(u1min, u1(k))
    u1max = max(u1max, u1(k))
  enddo

end subroutine measure_output_bounds

!> Calculates the change in interface positions based on h1 and h2
subroutine dzFromH1H2( n1, h1, n2, h2, dx )
  integer,            intent(in)  :: n1 !< Number of cells on source grid
  real, dimension(:), intent(in)  :: h1 !< Cell widths of source grid (size n1) [H]
  integer,            intent(in)  :: n2 !< Number of cells on target grid
  real, dimension(:), intent(in)  :: h2 !< Cell widths of target grid (size n2) [H]
  real, dimension(:), intent(out) :: dx !< Change in interface position (size n2+1) [H]
  ! Local variables
  integer :: k
  real :: x1, x2 ! Interface positions [H]

  x1 = 0.
  x2 = 0.
  dx(1) = 0.
  do K = 1, max(n1,n2)
    if (k <= n1) x1 = x1 + h1(k) ! Interface k+1, right of source cell k
    if (k <= n2) then
      x2 = x2 + h2(k) ! Interface k+1, right of target cell k
      dx(K+1) = x2 - x1 ! Change of interface k+1, target - source
    endif
  enddo

end subroutine dzFromH1H2

!> Constructor for remapping control structure
subroutine initialize_remapping( CS, remapping_scheme, boundary_extrapolation, &
                check_reconstruction, check_remapping, force_bounds_in_subcell, &
                force_bounds_in_target, better_force_bounds_in_target, offset_tgt_summation, &
                om4_remap_via_sub_cells, answers_2018, answer_date, nk, &
                h_neglect, h_neglect_edge)
  ! Arguments
  type(remapping_CS), intent(inout) :: CS !< Remapping control structure
  character(len=*),   intent(in)    :: remapping_scheme !< Remapping scheme to use
  logical, optional,  intent(in)    :: boundary_extrapolation !< Indicate to extrapolate in boundary cells
  logical, optional,  intent(in)    :: check_reconstruction !< Indicate to check reconstructions
  logical, optional,  intent(in)    :: check_remapping !< Indicate to check results of remapping
  logical, optional,  intent(in)    :: force_bounds_in_subcell !< Force subcells values to be bounded
  logical, optional,  intent(in)    :: force_bounds_in_target !< Force target values to be bounded
  logical, optional,  intent(in)    :: better_force_bounds_in_target !< Force target values to be bounded
  logical, optional,  intent(in)    :: offset_tgt_summation !< Use an offset when summing sub-cells
  logical, optional,  intent(in)    :: om4_remap_via_sub_cells !< If true, use OM4 remapping algorithm
  logical, optional,  intent(in)    :: answers_2018 !< If true use older, less accurate expressions.
  integer, optional,  intent(in)    :: answer_date  !< The vintage of the expressions to use
  real,    optional,  intent(in)    :: h_neglect !< A negligibly small width for the purpose of cell
                                                 !! reconstructions in the same units as h0 [H]
  real,    optional,  intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose of edge
                                                      !! value calculations in the same units as h0 [H].
  integer, optional,  intent(in)    :: nk !< Number of levels to initialize reconstruction class with

  ! Note that remapping_scheme is mandatory for initialize_remapping()
  call remapping_set_param(CS, &
               remapping_scheme=remapping_scheme, &
               boundary_extrapolation=boundary_extrapolation,  &
               check_reconstruction=check_reconstruction, &
               check_remapping=check_remapping, &
               force_bounds_in_subcell=force_bounds_in_subcell, &
               om4_remap_via_sub_cells=om4_remap_via_sub_cells, &
               force_bounds_in_target=force_bounds_in_target, &
               better_force_bounds_in_target=better_force_bounds_in_target, &
               offset_tgt_summation=offset_tgt_summation, &
               answers_2018=answers_2018, &
               answer_date=answer_date, &
               nk=nk, &
               h_neglect=h_neglect, &
               h_neglect_edge=h_neglect_edge)

end subroutine initialize_remapping

!> Changes the method of reconstruction
!! Use this routine to parse a string parameter specifying the reconstruction
!! and re-allocates work arrays appropriately. It is called from
!! initialize_remapping but can be called from an external module too.
subroutine setReconstructionType(string,CS)
  character(len=*),   intent(in)    :: string !< String to parse for method
  type(remapping_CS), intent(inout) :: CS !< Remapping control structure
  ! Local variables
  integer :: degree
  degree = -99
  if (associated(CS%reconstruction)) then
    ! We have a choice of being careless and allowing easy re-use (e.g. when testing)
    CS%remapping_scheme = -911
    call CS%reconstruction%destroy()
    deallocate( CS%reconstruction )
    ! or being careful and make sure we've properly clean up...
    !  call MOM_error(FATAL, "setReconstructionType: "//&
    !   "Recon1d type is already associated when initializing.")
  endif
  select case ( uppercase(trim(string)) )
    case ("PCM")
      CS%remapping_scheme = REMAPPING_PCM
      degree = 0
    case ("PLM")
      CS%remapping_scheme = REMAPPING_PLM
      degree = 1
    case ("PLM_HYBGEN")
      CS%remapping_scheme = REMAPPING_PLM_HYBGEN
      degree = 1
    case ("PPM_CW")
      CS%remapping_scheme = REMAPPING_PPM_CW
      degree = 2
    case ("PPM_H4")
      CS%remapping_scheme = REMAPPING_PPM_H4
      degree = 2
    case ("PPM_IH4")
      CS%remapping_scheme = REMAPPING_PPM_IH4
      degree = 2
    case ("PPM_HYBGEN")
      CS%remapping_scheme = REMAPPING_PPM_HYBGEN
      degree = 2
    case ("WENO_HYBGEN")
      CS%remapping_scheme = REMAPPING_WENO_HYBGEN
      degree = 2
    case ("PQM_IH4IH3")
      CS%remapping_scheme = REMAPPING_PQM_IH4IH3
      degree = 4
    case ("PQM_IH6IH5")
      CS%remapping_scheme = REMAPPING_PQM_IH6IH5
      degree = 4
    case ("C_PCM")
      allocate( PCM :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_PLM_CW")
      allocate( PLM_CW :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_PLM_HYBGEN")
      allocate( PLM_hybgen :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_MPLM_WA")
      allocate( MPLM_WA :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_EMPLM_WA")
      allocate( EMPLM_WA :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_MPLM_WA_POLY")
      allocate( MPLM_WA_poly :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_EMPLM_WA_POLY")
      allocate( EMPLM_WA_poly :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_PLM_CWK")
      allocate( PLM_CWK :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_MPLM_CWK")
      allocate( MPLM_CWK :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_EMPLM_CWK")
      allocate( EMPLM_CWK :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_PPM_CW")
      allocate( PPM_CW :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_PPM_HYBGEN")
      allocate( PPM_hybgen :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_PPM_CWK")
      allocate( PPM_CWK :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_EPPM_CWK")
      allocate( EPPM_CWK :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_PPM_H4_2019")
      allocate( PPM_H4_2019 :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case ("C_PPM_H4_2018")
      allocate( PPM_H4_2018 :: CS%reconstruction )
      CS%remapping_scheme = REMAPPING_VIA_CLASS
    case default
      call MOM_error(FATAL, "setReconstructionType: "//&
       "Unrecognized choice for REMAPPING_SCHEME ("//trim(string)//").")
  end select

  CS%degree = degree

end subroutine setReconstructionType

!> Destrcutor for remapping control structure
subroutine end_remapping(CS)
  type(remapping_CS), intent(inout) :: CS !< Remapping control structure

  CS%degree = 0

end subroutine end_remapping

!> Test if interpolate_column() produces the wrong answer
subroutine test_interp(test, msg, nsrc, h_src, u_src, ndest, h_dest, u_true)
  type(testing),         intent(inout) :: test   !< Unit testing convenience functions
  character(len=*),         intent(in) :: msg    !< Message to label test
  integer,                  intent(in) :: nsrc   !< Number of source cells
  real, dimension(nsrc),    intent(in) :: h_src  !< Thickness of source cells [H]
  real, dimension(nsrc+1),  intent(in) :: u_src  !< Values at source cell interfaces [A]
  integer,                  intent(in) :: ndest  !< Number of destination cells
  real, dimension(ndest),   intent(in) :: h_dest !< Thickness of destination cells [H]
  real, dimension(ndest+1), intent(in) :: u_true !< Correct value at destination cell interfaces [A]
  ! Local variables
  real, dimension(ndest+1) :: u_dest ! Interpolated value at destination cell interfaces [A]

  ! Interpolate from src to dest
  call interpolate_column(nsrc, h_src, u_src, ndest, h_dest, u_dest, .true.)
  call test%real_arr(ndest, u_dest, u_true, msg)
end subroutine test_interp

!> Test if reintegrate_column() produces the wrong answer
subroutine test_reintegrate(test, msg, nsrc, h_src, uh_src, ndest, h_dest, uh_true)
  type(testing),       intent(inout) :: test    !< Unit testing convenience functions
  character(len=*),       intent(in) :: msg     !< Message to label test
  integer,                intent(in) :: nsrc    !< Number of source cells
  real, dimension(nsrc),  intent(in) :: h_src   !< Thickness of source cells [H]
  real, dimension(nsrc),  intent(in) :: uh_src  !< Values of source cell stuff [A H]
  integer,                intent(in) :: ndest   !< Number of destination cells
  real, dimension(ndest), intent(in) :: h_dest  !< Thickness of destination cells [H]
  real, dimension(ndest), intent(in) :: uh_true !< Correct value of destination cell stuff [A H]
  ! Local variables
  real, dimension(ndest) :: uh_dest ! Reintegrated value on destination cells [A H]

  ! Interpolate from src to dest
  call reintegrate_column(nsrc, h_src, uh_src, ndest, h_dest, uh_dest)
  call test%real_arr(ndest, uh_dest, uh_true, msg)

end subroutine test_reintegrate

!> Test class-based remapping for internal consistency on random data
subroutine test_recon_consistency(test, scheme, n0, niter, h_neglect)
  type(testing),      intent(inout) :: test    !< Unit testing convenience functions
  character(len=*),   intent(in)    :: scheme  !< Name of scheme to use
  integer,            intent(in)    :: n0      !< Number of source cells
  integer,            intent(in)    :: niter   !< Number of randomized columns to try
  real,               intent(in)    :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  ! Local
  type(remapping_CS) :: remapCS !< Remapping control structure
  real :: h0(n0) ! Source grid [H but really nondim]
  real :: u0(n0) ! Source values [A]
  logical :: error ! Indicates a divergence
  integer :: iter ! Loop counter
  integer :: seed_size ! Number of integers used by seed
  integer, allocatable :: seed(:) ! Random number seed
  character(len=8) :: label ! Generated label

  call initialize_remapping(remapCS, scheme, nk=n0, h_neglect=h_neglect, &
                            force_bounds_in_subcell=.false. )

  call random_seed(size=seed_size)
  allocate( seed(seed_Size) )
  seed(:) = 102030405 ! Repeatable sequences
  call random_seed(put=seed)

  error = .false.
  do iter = 1, niter
    call random_number( h0 ) ! In range 0-1
    h0(:) = max(0., h0(:) - 0.05) ! Make 5% of values equal to zero
    call random_number( u0 ) ! In range 0-1

    call remapCS%reconstruction%reconstruct(h0, u0)
    if ( remapCS%reconstruction%check_reconstruction(h0, u0) ) then
      if ( .not. error ) then ! Only dump first error
        print *,'iter=',iter
        print *,'h0',h0
        print *,'u0',u0
      endif
      error = .true.
    endif

  enddo

  write(label(1:8),'(i8)') niter
  call test%test( error, trim(adjustl(label))//' consistency tests of '//scheme )

  call remapCS%reconstruction%destroy()

end subroutine test_recon_consistency

!> Test that remapping a uniform field remains uniform
subroutine test_preserve_uniform(test, scheme, n0, niter, h_neglect)
  type(testing),      intent(inout) :: test    !< Unit testing convenience functions
  character(len=*),   intent(in)    :: scheme  !< Name of scheme to use
  integer,            intent(in)    :: n0      !< Number of source cells
  integer,            intent(in)    :: niter   !< Number of randomized columns to try
  real,               intent(in)    :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  ! Local
  type(remapping_CS) :: remapCS !< Remapping control structure
  real :: h0(n0), h1(n0) ! Source and target grids [H but really nondim]
  real :: u0(n0), u1(n0) ! Source and target values [A]
  logical :: error ! Indicates a divergence
  integer :: iter ! Loop counter
  integer :: seed_size ! Number of integers used by seed
  integer, allocatable :: seed(:) ! Random number seed
  character(len=8) :: label ! Generated label

  call initialize_remapping(remapCS, scheme, nk=n0, h_neglect=h_neglect, &
                            force_bounds_in_subcell=.true., &
                            force_bounds_in_target=.true., &
                            better_force_bounds_in_target=.true., &
                            offset_tgt_summation=.false., &
                            om4_remap_via_sub_cells=.false.)

  call random_seed(size=seed_size)
  allocate( seed(seed_Size) )
  seed(:) = 102030405 ! Repeatable sequences
  call random_seed(put=seed)

  error = .false.
  do iter = 1, niter
    call random_number( h0 ) ! In range 0-1
    h0(:) = max(0., h0(:) - 0.05) ! Make 5% of values equal to zero
    call random_number( h1 ) ! In range 0-1
    h1(:) = max(0., h1(:) - 0.05) ! Make 5% of values equal to zero
    call random_number( u0(1) ) ! In range 0-1
    u0(:) = u0(1) ! Make u0 uniform

    call remapping_core_h( remapCS, n0, h0, u0, n0, h1, u1 )
    if ( maxval( abs( u1(:) - u0(1) ) ) > 0. ) then
      if ( .not. error ) then ! Only dump first error
        print *,'iter=',iter
        print *,'u0(1)',u0(1)
        print *,'u1',u1
        print *,'u1-u0(1)',u1 - u0(1)
      endif
      error = .true.
    endif

  enddo

  write(label(1:8),'(i8)') niter
  call test%test( error, trim(adjustl(label))//' uniformity tests of '//scheme )

end subroutine test_preserve_uniform

!> Test that remapping to the same grid preserves answers
!!
!! Notes:
!! 1) this test is currently imperfect since occasionally we see round-off
!!    implying that     ( A * B ) / A != B
!! 2) this test does not work for vanished layers
subroutine test_unchanged_grid(test, scheme, n0, niter, h_neglect)
  type(testing),      intent(inout) :: test    !< Unit testing convenience functions
  character(len=*),   intent(in)    :: scheme  !< Name of scheme to use
  integer,            intent(in)    :: n0      !< Number of source cells
  integer,            intent(in)    :: niter   !< Number of randomized columns to try
  real,               intent(in)    :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  ! Local
  type(remapping_CS) :: remapCS !< Remapping control structure
  real :: h0(n0), h1(n0) ! Source and target grids [H but really nondim]
  real :: u0(n0), u1(n0) ! Source and target values [A]
  logical :: error ! Indicates a divergence
  integer :: iter ! Loop counter
  character(len=8) :: label ! Generated label

  call initialize_remapping(remapCS, scheme, nk=n0, h_neglect=h_neglect, &
                            force_bounds_in_subcell=.true., &
                            force_bounds_in_target=.false., &
                            better_force_bounds_in_target=.true., &
                            offset_tgt_summation=.true., &
                            om4_remap_via_sub_cells=.false.)

  error = .false.
  do iter = 1, niter
    call random_number( h0 ) ! In range 0-1
    h0(:) = max(0., h0(:) - 0.00) ! Note we do NOT test with vanished layers
    h1(:) = h0(:) ! Exact copy
    call random_number( u0 ) ! In range 0-1

    call remapping_core_h( remapCS, n0, h0, u0, n0, h1, u1 )
    if ( maxval( abs( u1(:) - u0(:) ) ) > epsilon(h0(1)) * maxval( abs( u0 ) ) ) then
      if ( .not. error ) then ! Only dump first error
        print *,'iter=',iter
        print *,'h0',h0
        print *,'u0',u0
        print *,'u1',u1
        print *,'u1-u0',u1 - u0
      endif
      error = .true.
    endif

  enddo

  write(label(1:8),'(i8)') niter
  call test%test( error, trim(adjustl(label))//' unchanged grid tests of '//scheme )

  call remapCS%reconstruction%destroy()

end subroutine test_unchanged_grid

!> Test class-based remapping bitwise reproduces original implementation
subroutine compare_two_schemes(test, CS1, CS2, n0, n1, niter, msg)
  type(testing),      intent(inout) :: test  !< Unit testing convenience functions
  type(remapping_CS), intent(inout) :: CS1   !< Remapping control structure configured for
                                             !! original implementation
  type(remapping_CS), intent(inout) :: CS2   !< Remapping control structure configured for
                                             !! class-based implementation
  integer,            intent(in)    :: n0    !< Number of source cells
  integer,            intent(in)    :: n1    !< Number of destination cells
  integer,            intent(in)    :: niter !< Number of randomized columns to try
  character(len=*),   intent(in)    :: msg   !< Message to label test
  ! Local
  real :: h0(n0), h1(n1) ! Source and target grids [H but really nondim]
  real :: u0(n0), u1(n1), u2(n1)  ! Source and two target values [A]
  logical :: error ! Indicates a divergence
  integer :: iter ! Loop counter
  integer :: seed_size ! Number of integers used by seed
  integer, allocatable :: seed(:) ! Random number seed
  character(len=8) :: label ! Generated label

  call random_seed(size=seed_size)
  allocate( seed(seed_Size) )
  seed(:) = 102030405 ! Repeatable sequences
  call random_seed(put=seed)

  error = .false.
  do iter = 1, niter
    call random_number( h0 ) ! In range 0-1
    h0(:) = max(0., h0(:) - 0.00) ! Make 5% of values equal to zero
    h0(:) = h0(:) / sum( h0 ) ! Approximately normalize to total depth of 1
    call random_number(h1) ! In range 0-1
    h1(:) = max(0., h1(:) - 0.00) ! Make 5% of values equal to zero
    h1(:) = h1(:) / sum( h1 ) ! Approximately normalize to total depth of 1
    call random_number( u0 ) ! In range 0-1

    call remapping_core_h( CS1, n0, h0, u0, n1, h1, u1 )
    call remapping_core_h( CS2, n0, h0, u0, n1, h1, u2 )
    error = sum( abs( u2(:) - u1(:) ) ) > 0.
    if (error) then
      print *,'iter=',iter
      print *,'h1',h1
      print *,'h0',h0
      print *,'u0',u0
      print *,'u1',u1
      print *,'u2',u2
      print *,'e',u2-u1
    ! CS1%debug = .true.
    ! call remapping_core_h( CS1, n0, h0, u0, n1, h1, u1 )
    ! CS2%debug = .true.
    ! call remapping_core_h( CS2, n0, h0, u0, n1, h1, u2 )
      exit
    endif
  enddo

  write(label(1:8),'(i8)') niter
  call test%test( error, trim(adjustl(label))//' comparisons of '//msg )

end subroutine compare_two_schemes

!> Runs unit tests on remapping functions.
!! Should only be called from a single/root thread
!! Returns True if a test fails, otherwise False
logical function remapping_unit_tests(verbose, num_comp_samp)
  logical, intent(in) :: verbose !< If true, write results to stdout
  integer, optional, intent(in) :: num_comp_samp !< If present, number of samples to
                                 !! try comparing class-based cade against OM4 code
  ! Local variables
  integer :: n0, n1, n2
  real, allocatable :: h0(:), h1(:), h2(:) ! Thicknesses for test columns [H]
  real, allocatable :: u0(:), u1(:), u2(:) ! Values for test profiles [A]
  real, allocatable :: dx1(:) ! Change in interface position [H]
  type(remapping_CS) :: CS, CS2 !< Remapping control structures
  real, allocatable, dimension(:,:) :: ppoly0_E     ! Edge values of polynomials [A]
  real, allocatable, dimension(:,:) :: ppoly0_S     ! Edge slopes of polynomials [A H-1]
  real, allocatable, dimension(:,:) :: ppoly0_coefs ! Coefficients of polynomials [A]
  real, allocatable, dimension(:) :: h_sub, h0_eff ! Subgrid and effective source thicknesses [H]
  real, allocatable, dimension(:) :: u_sub, uh_sub ! Subgrid values and totals [A, A H]
  real :: u02_err ! Error in remaping [A]
  integer, allocatable, dimension(:) :: isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src ! Indices
  integer :: answer_date  ! The vintage of the expressions to test
  real :: err                         ! Errors in the remapped thicknesses [H] or values [A]
  real :: h_neglect, h_neglect_edge   ! Tiny thicknesses used in remapping [H]
  integer :: seed_size ! Number of integers used by seed
  integer, allocatable :: seed(:) ! Random number seed
  type(testing) :: test ! Unit testing convenience functions
  integer :: om4 ! Loop parameter, 0 or 1
  integer :: ntests ! Number of iterations when brute force testing
  character(len=4) :: om4_tag ! Generated label
  type(PCM) :: PCM
  type(PLM_CW) :: PLM_CW
  type(PLM_hybgen) :: PLM_hybgen
  type(MPLM_WA) :: MPLM_WA
  type(EMPLM_WA) :: EMPLM_WA
  type(MPLM_WA_poly) :: MPLM_WA_poly
  type(EMPLM_WA_poly) :: EMPLM_WA_poly
  type(PLM_CWK) :: PLM_CWK
  type(MPLM_CWK) :: MPLM_CWK
  type(EMPLM_CWK) :: EMPLM_CWK
  type(PPM_H4_2019) :: PPM_H4_2019
  type(PPM_H4_2018) :: PPM_H4_2018
  type(PPM_CW) :: PPM_CW
  type(PPM_hybgen) :: PPM_hybgen
  type(PPM_CWK) :: PPM_CWK
  type(EPPM_CWK) :: EPPM_CWK

  call test%set( verbose=verbose ) ! Sets the verbosity flag in test
! call test%set( stop_instantly=.true. ) ! While debugging

  answer_date = 20190101 ! 20181231
  h_neglect = 1.0e-30
  h_neglect_edge = h_neglect ; if (answer_date < 20190101) h_neglect_edge = 1.0e-10

  if (verbose) write(test%stdout,*) '  ===== MOM_remapping: remapping_unit_tests ================='

  if (verbose) write(test%stdout,*) '  - - - - - 1st generation tests - - - - -'

  call initialize_remapping(CS, 'PPM_H4', answer_date=answer_date, &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)

  ! Profile 0: 4 layers of thickness 0.75 and total depth 3, with du/dz=8
  n0 = 4
  allocate( h0(n0), u0(n0) )
  h0 = (/0.75, 0.75, 0.75, 0.75/)
  u0 = (/9., 3., -3., -9./)

  ! Profile 1: 3 layers of thickness 1.0 and total depth 3
  n1 = 3
  allocate( h1(n1), u1(n1), dx1(n1+1) )
  h1 = (/1.0, 1.0, 1.0/)

  ! Profile 2: 6 layers of thickness 0.5 and total depth 3
  n2 = 6
  allocate( h2(n2), u2(n2) )
  h2 = (/0.5, 0.5, 0.5, 0.5, 0.5, 0.5/)

  ! Mapping u1 from h1 to h2
  call dzFromH1H2( n0, h0, n1, h1, dx1 )
  call remapping_core_w( CS, n0, h0, u0, n1, dx1, u1 )
  call test%real_arr(3, u1, (/8.,0.,-8./), 'remapping_core_w() PPM_H4')

  allocate(ppoly0_E(n0,2), ppoly0_S(n0,2), ppoly0_coefs(n0,CS%degree+1))
  ppoly0_E(:,:) = 0.0
  ppoly0_S(:,:) = 0.0
  ppoly0_coefs(:,:) = 0.0

  call initialize_remapping(CS, 'PPM_H4', force_bounds_in_subcell=.false., answer_date=answer_date)

  call remapping_core_h( CS, n0, h0, u0, n2, h2, u2, net_err=err )
  call test%real_arr(6, u2, (/10.,6.,2.,-2.,-6.,-10./), 'remapping_core_h() 2')

  call remapping_core_h( CS, n0, h0, u0, 6, (/.125,.125,.125,.125,.125,.125/), u2, net_err=err )
  call test%real_arr(6, u2, (/11.5,10.5,9.5,8.5,7.5,6.5/), 'remapping_core_h() 3')

  call remapping_core_h( CS, n0, h0, u0, 3, (/2.25,1.5,1./), u2, net_err=err )
  call test%real_arr(3, u2, (/3.,-10.5,-12./), 'remapping_core_h() 4')

  deallocate(h0, u0, h1, u1, h2, u2, ppoly0_E, ppoly0_S, ppoly0_coefs)
  call end_remapping(CS)

  ! ===============================================
  ! This section tests the reconstruction functions
  ! ===============================================
  if (verbose) write(test%stdout,*) '  - - - - - reconstruction tests - - - - -'

  allocate( ppoly0_coefs(5,6), ppoly0_E(5,2), ppoly0_S(5,2), u2(2) )

  call PCM_reconstruction(3, (/1.,2.,4./), &
                          ppoly0_E(1:3,:), ppoly0_coefs(1:3,:) )
  call test%real_arr(3, ppoly0_E(:,1), (/1.,2.,4./), 'PCM: left edges')
  call test%real_arr(3, ppoly0_E(:,2), (/1.,2.,4./), 'PCM: right edges')
  call test%real_arr(3, ppoly0_coefs(:,1), (/1.,2.,4./), 'PCM: P0')

  call PLM_reconstruction(3, (/1.,1.,1./), (/1.,3.,5./), &
                          ppoly0_E(1:3,:), ppoly0_coefs(1:3,:), h_neglect )
  call test%real_arr(3, ppoly0_E(:,1), (/1.,2.,5./), 'Unlim PLM: left edges')
  call test%real_arr(3, ppoly0_E(:,2), (/1.,4.,5./), 'Unlim PLM: right edges')
  call test%real_arr(3, ppoly0_coefs(:,1), (/1.,2.,5./), 'Unlim PLM: P0')
  call test%real_arr(3, ppoly0_coefs(:,2), (/0.,2.,0./), 'Unlim PLM: P1')

  call PLM_reconstruction(3, (/1.,1.,1./), (/1.,2.,7./), &
                          ppoly0_E(1:3,:), ppoly0_coefs(1:3,:), h_neglect )
  call test%real_arr(3, ppoly0_E(:,1), (/1.,1.,7./), 'Left lim PLM: left edges')
  call test%real_arr(3, ppoly0_E(:,2), (/1.,3.,7./), 'Left lim PLM: right edges')
  call test%real_arr(3, ppoly0_coefs(:,1), (/1.,1.,7./), 'Left lim PLM: P0')
  call test%real_arr(3, ppoly0_coefs(:,2), (/0.,2.,0./), 'Left lim PLM: P1')

  call PLM_reconstruction(3, (/1.,1.,1./), (/1.,6.,7./), &
                          ppoly0_E(1:3,:), ppoly0_coefs(1:3,:), h_neglect )
  call test%real_arr(3, ppoly0_E(:,1), (/1.,5.,7./), 'Right lim PLM: left edges')
  call test%real_arr(3, ppoly0_E(:,2), (/1.,7.,7./), 'Right lim PLM: right edges')
  call test%real_arr(3, ppoly0_coefs(:,1), (/1.,5.,7./), 'Right lim PLM: P0')
  call test%real_arr(3, ppoly0_coefs(:,2), (/0.,2.,0./), 'Right lim PLM: P1')

  call PLM_reconstruction(3, (/1.,2.,3./), (/1.,4.,9./), &
                          ppoly0_E(1:3,:), ppoly0_coefs(1:3,:), h_neglect )
  call test%real_arr(3, ppoly0_E(:,1), (/1.,2.,9./), 'Non-uniform line PLM: left edges')
  call test%real_arr(3, ppoly0_E(:,2), (/1.,6.,9./), 'Non-uniform line PLM: right edges')
  call test%real_arr(3, ppoly0_coefs(:,1), (/1.,2.,9./), 'Non-uniform line PLM: P0')
  call test%real_arr(3, ppoly0_coefs(:,2), (/0.,4.,0./), 'Non-uniform line PLM: P1')

  call edge_values_explicit_h4(5, (/1.,1.,1.,1.,1./), (/1.,3.,5.,7.,9./), &
                               ppoly0_E, h_neglect=1e-10, answer_date=answer_date )
  ! The next two tests currently fail due to roundoff, but pass when given a reasonable tolerance.
  call test%real_arr(5, ppoly0_E(:,1), (/0.,2.,4.,6.,8./), 'Line H4: left edges', tol=8.0e-15)
  call test%real_arr(5, ppoly0_E(:,2), (/2.,4.,6.,8.,10./), 'Line H4: right edges', tol=1.0e-14)

  ppoly0_E(:,1) = (/0.,2.,4.,6.,8./)
  ppoly0_E(:,2) = (/2.,4.,6.,8.,10./)
  call PPM_reconstruction(5, (/1.,1.,1.,1.,1./), (/1.,3.,5.,7.,9./), ppoly0_E(1:5,:), &
                              ppoly0_coefs(1:5,:), h_neglect, answer_date=answer_date )
  call test%real_arr(5, ppoly0_coefs(:,1), (/1.,2.,4.,6.,9./), 'Line PPM: P0')
  call test%real_arr(5, ppoly0_coefs(:,2), (/0.,2.,2.,2.,0./), 'Line PPM: P1')
  call test%real_arr(5, ppoly0_coefs(:,3), (/0.,0.,0.,0.,0./), 'Line PPM: P2')

  call edge_values_explicit_h4( 5, (/1.,1.,1.,1.,1./), (/1.,1.,7.,19.,37./), ppoly0_E, &
                                h_neglect=1e-10, answer_date=answer_date )
  ! The next two tests are now passing when answer_date >= 20190101, but otherwise only work to roundoff.
  call test%real_arr(5, ppoly0_E(:,1), (/3.,0.,3.,12.,27./), 'Parabola H4: left edges', tol=2.7e-14)
  call test%real_arr(5, ppoly0_E(:,2), (/0.,3.,12.,27.,48./), 'Parabola H4: right edges', tol=4.8e-14)
  ppoly0_E(:,1) = (/0.,0.,3.,12.,27./)
  ppoly0_E(:,2) = (/0.,3.,12.,27.,48./)
  call PPM_reconstruction(5, (/1.,1.,1.,1.,1./), (/0.,1.,7.,19.,37./), ppoly0_E(1:5,:), &
                          ppoly0_coefs(1:5,:), h_neglect, answer_date=answer_date )
  call test%real_arr(5, ppoly0_E(:,1), (/0.,0.,3.,12.,37./), 'Parabola PPM: left edges')
  call test%real_arr(5, ppoly0_E(:,2), (/0.,3.,12.,27.,37./), 'Parabola PPM: right edges')
  call test%real_arr(5, ppoly0_coefs(:,1), (/0.,0.,3.,12.,37./), 'Parabola PPM: P0')
  call test%real_arr(5, ppoly0_coefs(:,2), (/0.,0.,6.,12.,0./), 'Parabola PPM: P1')
  call test%real_arr(5, ppoly0_coefs(:,3), (/0.,3.,3.,3.,0./), 'Parabola PPM: P2')

  ppoly0_E(:,1) = (/0.,0.,6.,10.,15./)
  ppoly0_E(:,2) = (/0.,6.,12.,17.,15./)
  call PPM_reconstruction(5, (/1.,1.,1.,1.,1./), (/0.,5.,7.,16.,15./), ppoly0_E(1:5,:), &
                          ppoly0_coefs(1:5,:), h_neglect, answer_date=answer_date )
  call test%real_arr(5, ppoly0_E(:,1), (/0.,3.,6.,16.,15./), 'Limits PPM: left edges')
  call test%real_arr(5, ppoly0_E(:,2), (/0.,6.,9.,16.,15./), 'Limits PPM: right edges')
  call test%real_arr(5, ppoly0_coefs(:,1), (/0.,3.,6.,16.,15./), 'Limits PPM: P0')
  call test%real_arr(5, ppoly0_coefs(:,2), (/0.,6.,0.,0.,0./), 'Limits PPM: P1')
  call test%real_arr(5, ppoly0_coefs(:,3), (/0.,-3.,3.,0.,0./), 'Limits PPM: P2')

  deallocate(ppoly0_E, ppoly0_S, ppoly0_coefs, u2)

  ! ==============================================================
  ! This section tests the components of remapping_core_h()
  ! ==============================================================

  if (verbose) write(test%stdout,*) '  - - - - - remapping algororithm tests - - - - -'

  ! Test 1: n0=2, n1=3  Maps uniform grids with one extra target layer and no implicitly-vanished interior sub-layers
  ! h_src =     |       3        |        3       |
  ! h_tgt =     |     2    |     2     |    2     |
  ! h_sub =     |0|   2    |  1  |  1  |    2   |0|
  ! isrc_start  |1               |  4             |
  ! isrc_end    |             3  |          5     |
  ! isrc_max    |     2          |          5     |
  ! itgt_start  |1         |  3        |    5     |
  ! itgt_end    |     2    |        4  |         6|
  ! isub_src    |1|   1    |  1  |  2  |    2   |2|
  allocate( h_sub(6), h0_eff(2), isrc_start(2), isrc_end(2), isrc_max(2), itgt_start(3), itgt_end(3), isub_src(6) )
  call intersect_src_tgt_grids( 2, (/3., 3./), &  ! n0, h0
                                3, (/2., 2., 2./), &  ! n1, h1
                                h_sub, h0_eff, &
                                isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )
  if (verbose) write(test%stdout,*) "intersect_src_tgt_grids test 1: n0=2, n1=3"
  if (verbose) write(test%stdout,*) "  h_src =     |     3     |     3     |"
  if (verbose) write(test%stdout,*) "  h_tgt =     |   2   |   2   |   2   |"
  call test%real_arr(6, h_sub, (/0.,2.,1.,1.,2.,0./), 'h_sub')
  call test%real_arr(2, h0_eff, (/3.,3./), 'h0_eff')
  call test%int_arr(2, isrc_start, (/1,4/), 'isrc_start')
  call test%int_arr(2, isrc_end, (/3,5/), 'isrc_end')
  call test%int_arr(2, isrc_max, (/2,5/), 'isrc_max')
  call test%int_arr(3, itgt_start, (/1,3,5/), 'itgt_start')
  call test%int_arr(3, itgt_end, (/2,4,6/), 'itgt_end')
  call test%int_arr(6, isub_src, (/1,1,1,2,2,2/), 'isub_src')
  deallocate( h_sub, h0_eff, isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )

  ! Test 2: n0=3, n1=2  Reverses "test 1" with more source than target layers
  ! h_src =     |    2    |     2     |    2    |
  ! h_tgt =     |      3        |        3      |
  ! h_sub =     |0|  2    |  1  |  1  |    2  |0|
  ! isrc_start  |1        |  3        |    5    |
  ! isrc_end    |    2    |        4  |    5    |
  ! isrc_max    |    2    |        4  |    5    |
  ! itgt_start  |1              |  4            |
  ! itgt_end    |            3  |              6|
  ! isub_src    |1|  1    |  2  |  2  |    3  |3|
  allocate( h_sub(6), h0_eff(3), isrc_start(3), isrc_end(3), isrc_max(3), itgt_start(2), itgt_end(2), isub_src(6) )
  call intersect_src_tgt_grids( 3, (/2., 2., 2./), &  ! n0, h0
                                2, (/3., 3./), &  ! n1, h1
                                h_sub, h0_eff, &
                                isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )
  if (verbose) write(test%stdout,*) "intersect_src_tgt_grids test 2: n0=3, n1=2"
  if (verbose) write(test%stdout,*) "  h_src =     |   2   |   2   |   2   |"
  if (verbose) write(test%stdout,*) "  h_tgt =     |     3     |     3     |"
  call test%real_arr(6, h_sub, (/0.,2.,1.,1.,2.,0./), 'h_sub')
  call test%real_arr(3, h0_eff, (/2.,2.,2./), 'h0_eff')
  call test%int_arr(3, isrc_start, (/1,3,5/), 'isrc_start')
  call test%int_arr(3, isrc_end, (/2,4,5/), 'isrc_end')
  call test%int_arr(3, isrc_max, (/2,4,5/), 'isrc_max')
  call test%int_arr(2, itgt_start, (/1,4/), 'itgt_start')
  call test%int_arr(2, itgt_end, (/3,6/), 'itgt_end')
  call test%int_arr(6, isub_src, (/1,1,2,2,3,3/), 'isub_src')
  deallocate( h_sub, h0_eff, isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )

  ! Test 3: n0=2, n1=3  With aligned interfaces that lead to implicitly-vanished interior sub-layers
  n0 = 2 ; n1 = 3
  allocate( h0_eff(n0), isrc_start(n0), isrc_end(n0), isrc_max(n0), h0(n0), u0(n0) )
  allocate( itgt_start(n1), itgt_end(n1), h1(n1), u1(n1) )
  allocate( h_sub(n0+n1+1), isub_src(n0+n1+1) )
  u0 =         (/     2.    ,          5.         /)
  h0 =         (/     2.    ,          4.         /)
  h1 =         (/     2.    ,     2.   ,    2.    /)
  ! h_src =     |<-   2   ->|<-        4        ->|
  ! h_tgt =     |<-   2   ->|<-   2  ->|<-  2   ->|
  ! h_sub =    |0|<-  2  ->|0|<-  2  ->|<-  2  ->|0|
  ! isrc_start  |1          |3                    |
  ! isrc_end    |     2     |               5     |
  ! isrc_max    |     2     |               5     |
  ! itgt_start  |1          |     4    |    5     |
  ! itgt_end    |            3|   4    |         6|
  ! isub_src    |1|   1     |2|   2    |    2   |2|
  call intersect_src_tgt_grids( n0, h0, n1, h1, h_sub, h0_eff, &
                                isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )
  if (verbose) write(test%stdout,*) "intersect_src_tgt_grids test 3: n0=2, n1=3"
  if (verbose) write(test%stdout,*) "  h_src =     |   2   |       4       |"
  if (verbose) write(test%stdout,*) "  h_tgt =     |   2   |   2   |   2   |"
  call test%real_arr(6, h_sub, (/0.,2.,0.,2.,2.,0./), 'h_sub')
  call test%real_arr(2, h0_eff, (/2.,4./), 'h0_eff')
  call test%int_arr(2, isrc_start, (/1,3/), 'isrc_start')
  call test%int_arr(2, isrc_end, (/2,5/), 'isrc_end')
  call test%int_arr(2, isrc_max, (/2,5/), 'isrc_max')
  call test%int_arr(3, itgt_start, (/1,4,5/), 'itgt_start')
  call test%int_arr(3, itgt_end, (/3,4,6/), 'itgt_end')
  call test%int_arr(6, isub_src, (/1,1,2,2,2,2/), 'isub_src')
  allocate(ppoly0_coefs(n0,2), ppoly0_E(n0,2), ppoly0_S(n0,2))
  ! h_src =     |<-   2   ->|<-        4        ->|
  ! h_sub =    |0|<-  2  ->|0|<-  2  ->|<-  2  ->|0|
  ! u_src =     |     2     |          5          |
  !  edge =     |1         3|3                   7|
  ! u_sub =    |1|    2    |3|    4    |    6    |7|
  call PLM_reconstruction(n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect )
  call PLM_boundary_extrapolation(n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect)
  allocate(u_sub(n0+n1+1), uh_sub(n0+n1+1))
  call remap_src_to_sub_grid_om4(n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                             n1, h_sub, h0_eff, isrc_start, isrc_end, isrc_max, isub_src, &
                             INTEGRATION_PLM, .false., u_sub, uh_sub, u02_err)
  call test%real_arr(6, u_sub, (/1.,2.,3.,4.,6.,7./), 'u_sub om4')
  call remap_src_to_sub_grid(n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                             n1, h_sub, isrc_start, isrc_end, isrc_max, isub_src, &
                             INTEGRATION_PLM, .false., u_sub, uh_sub, u02_err)
  call test%real_arr(6, u_sub, (/1.,2.,3.,4.,6.,7./), 'u_sub')
  ! h_sub =    |0|<-  2  ->|0|<-  2  ->|<-  2  ->|0|
  ! u_sub =    |1|    2    |3|    4    |    6    |7|
  ! h_tgt =     |<-   2   ->|<-   2  ->|<-  2   ->|
  ! u_tgt =     |     2     |     4    |    6     |
  call remap_sub_to_tgt_grid(n0, n1, h1, h_sub, u_sub, uh_sub, itgt_start, itgt_end, &
                             .false., .false., .false., u1, u02_err)
  call test%real_arr(3, u1, (/2.,4.,6./), 'u1')
  call remap_sub_to_tgt_grid(n0, n1, h1, h_sub, u_sub, uh_sub, itgt_start, itgt_end, &
                             .true., .false., .false., u1, u02_err)
  call test%real_arr(3, u1, (/2.,4.,6./), 'u1.b')
  deallocate( ppoly0_coefs, ppoly0_E, ppoly0_S, u_sub, uh_sub, h0, u0, h1, u1)
  deallocate( h_sub, h0_eff, isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )

  ! Test 4: n0=2, n1=3  Incomplete target column, sum(h_tgt)<sum(h_src), useful for diagnostics
  n0 = 2 ; n1 = 3
  allocate( h0_eff(n0), isrc_start(n0), isrc_end(n0), isrc_max(n0), h0(n0), u0(n0) )
  allocate( itgt_start(n1), itgt_end(n1), h1(n1), u1(n1) )
  allocate( h_sub(n0+n1+1), isub_src(n0+n1+1) )
  u0 =         (/     2.    ,           5.         /)
  h0 =         (/     2.    ,           4.         /)
  h1 =         (/     2.    ,     2.   ,   1. /)
  ! h_src =     |<-   2   ->|<-         4         ->|
  ! h_tgt =     |<-   2   ->|<-   2   ->|< 1 >|
  ! h_sub =    |0|<-  2  ->|0|<-  2   ->|< 1 >|< 1 >|
  ! isrc_start  |1          |3                      |
  ! isrc_end    |     2     |                    6  |
  ! isrc_max    |     2     |     4                 |
  ! itgt_start  |1          |     4     |  5  |
  ! itgt_end    |          3|     4     |  5  |
  ! isub_src   |1|    1    |2|    2     |  2  |  2  |
  call intersect_src_tgt_grids( n0, h0, n1, h1, h_sub, h0_eff, &
                                isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )
  if (verbose) write(test%stdout,*) "intersect_src_tgt_grids test 4: n0=2, n1=3"
  if (verbose) write(test%stdout,*) "  h_src =     |   2   |       4       |"
  if (verbose) write(test%stdout,*) "  h_tgt =     |   2   |   2   | 1 |"
  call test%real_arr(6, h_sub, (/0.,2.,0.,2.,1.,1./), 'h_sub')
  call test%real_arr(2, h0_eff, (/2.,3./), 'h0_eff')
  call test%int_arr(2, isrc_start, (/1,3/), 'isrc_start')
  call test%int_arr(2, isrc_end, (/2,6/), 'isrc_end')
  call test%int_arr(2, isrc_max, (/2,4/), 'isrc_max')
  call test%int_arr(3, itgt_start, (/1,4,5/), 'itgt_start')
  call test%int_arr(3, itgt_end, (/3,4,5/), 'itgt_end')
  call test%int_arr(6, isub_src, (/1,1,2,2,2,2/), 'isub_src')
  allocate(ppoly0_coefs(n0,2), ppoly0_E(n0,2), ppoly0_S(n0,2))
  ! h_src =     |<-   2   ->|<-       4         ->|
  ! h_sub =     |0|<- 2   ->|0|<- 2 ->|<-1->|<-1->|
  ! u_src =     |     2     |         5           |
  !  edge =     |1         3|3                   7|
  ! u_sub =     |1|   2     |3|   4   | 5.5 | 6.5 |
  call PLM_reconstruction(n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect )
  call PLM_boundary_extrapolation(n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect)
  allocate(u_sub(n0+n1+1), uh_sub(n0+n1+1))
  call remap_src_to_sub_grid(2, (/2.,4./), (/2.,5./), ppoly0_E, ppoly0_coefs, &
                             3, h_sub, isrc_start, isrc_end, isrc_max, isub_src, &
                             INTEGRATION_PLM, .false., u_sub, uh_sub, u02_err)
  call test%real_arr(6, u_sub, (/1.,2.,3.,4.,5.5,6.5/), 'u_sub')
  deallocate( ppoly0_coefs, ppoly0_E, ppoly0_S, u_sub, uh_sub, h0, u0, h1, u1)
  deallocate( h_sub, h0_eff, isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )

  ! Test 5: n0=3, n1=2  Target column exceeds source column, sum(h_tgt)>sum(h_src), useful for diagnostics
  n0 = 3 ; n1 = 2
  allocate( h0_eff(n0), isrc_start(n0), isrc_end(n0), isrc_max(n0), h0(n0), u0(n0) )
  allocate( itgt_start(n1), itgt_end(n1), h1(n1), u1(n1) )
  allocate( h_sub(n0+n1+1), isub_src(n0+n1+1) )
  u0 =         (/     2.    ,     4.   ,  5.5 /)
  h0 =         (/     2.    ,     2.   ,   1. /)
  h1 =         (/     2.    ,           4.         /)
  ! h_src =     |<-   2   ->|<-   2   ->|< 1 >|
  ! h_tgt =     |<-   2   ->|<-         4         ->|
  ! h_sub =    |0|<-  2  ->|0|<-  2   ->|< 1 >|< 1 >|
  ! isrc_start  |1          |3          |  5  |
  ! isrc_end    |     2     |     4     |  5  |
  ! isrc_max    |     2     |     4     |  5  |
  ! itgt_start  |1          |     4                 |
  ! itgt_end    |            3|                6  |
  ! isub_src    |1|   1     |2|  2    |  3  |  3  |
  call intersect_src_tgt_grids( n0, h0, n1, h1, h_sub, h0_eff, &
                                isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )
  if (verbose) write(test%stdout,*) "intersect_src_tgt_grids test 5: n0=3, n1=2"
  if (verbose) write(test%stdout,*) "  h_src =     |   2   |   2   | 1 |"
  if (verbose) write(test%stdout,*) "  h_tgt =     |   2   |       4       |"
  call test%real_arr(6, h_sub, (/0.,2.,0.,2.,1.,1./), 'h_sub')
  call test%real_arr(3, h0_eff, (/2.,2.,1./), 'h0_eff')
  call test%int_arr(3, isrc_start, (/1,3,5/), 'isrc_start')
  call test%int_arr(3, isrc_end, (/2,4,5/), 'isrc_end')
  call test%int_arr(3, isrc_max, (/2,4,5/), 'isrc_max')
  call test%int_arr(2, itgt_start, (/1,4/), 'itgt_start')
  call test%int_arr(2, itgt_end, (/3,6/), 'itgt_end')
  call test%int_arr(6, isub_src, (/1,1,2,2,3,3/), 'isub_src')
  allocate(ppoly0_coefs(n0,2), ppoly0_E(n0,2), ppoly0_S(n0,2))
  ! h_src =     |<-   2   ->|<-   2   ->|< 1 >|
  ! h_sub =    |0|<-  2  ->|0|<-  2   ->|< 1 >|< 1 >|
  ! u_src =     |     2     |     4     | 5.5 |
  !  edge =     |1         3|3         5|5   6|
  ! u_sub =    |1|   2     |3|    4     | 5.5 |  6  |
  call PLM_reconstruction(n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect )
  call PLM_boundary_extrapolation(n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect)
  allocate(u_sub(n0+n1+1), uh_sub(n0+n1+1))
  call remap_src_to_sub_grid_om4(n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                             n1, h_sub, h0_eff, isrc_start, isrc_end, isrc_max, isub_src, &
                             INTEGRATION_PLM, .false., u_sub, uh_sub, u02_err)
  call test%real_arr(6, u_sub, (/1.,2.,3.,4.,5.5,6./), 'u_sub om4')
  call remap_src_to_sub_grid(n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                             n1, h_sub, isrc_start, isrc_end, isrc_max, isub_src, &
                             INTEGRATION_PLM, .false., u_sub, uh_sub, u02_err)
  call test%real_arr(6, u_sub, (/1.,2.,3.,4.,5.5,6./), 'u_sub')
  ! h_sub =    |0|<-  2  ->|0|<-  2   ->|< 1 >|< 1 >|
  ! u_sub =    |1|   2     |3|    4     | 5.5 |  6  |
  ! h_tgt =     |<-   2   ->|<-         4         ->|
  ! u_tgt =     |     2     |          4 7/8        |
  call remap_sub_to_tgt_grid(n0, n1, h1, h_sub, u_sub, uh_sub, itgt_start, itgt_end, &
                             .false., .false., .false., u1, u02_err)
  call test%real_arr(2, u1, (/2.,4.875/), 'u1')
  call remap_sub_to_tgt_grid(n0, n1, h1, h_sub, u_sub, uh_sub, itgt_start, itgt_end, &
                             .true., .false., .false., u1, u02_err)
  call test%real_arr(2, u1, (/2.,4.875/), 'u1.b')
  deallocate( ppoly0_coefs, ppoly0_E, ppoly0_S, u_sub, uh_sub, h0, u0, h1, u1)
  deallocate( h_sub, h0_eff, isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )

  ! Test 6: n0=3, n1=5  Source and targets with vanished layers
  n0 = 3 ; n1 = 5
  allocate( h0_eff(n0), isrc_start(n0), isrc_end(n0), isrc_max(n0), h0(n0), u0(n0) )
  allocate( itgt_start(n1), itgt_end(n1), h1(n1), u1(n1) )
  allocate( h_sub(n0+n1+1), isub_src(n0+n1+1) )
  u0 =         (/        2.       ,3.,        4.        /)
  h0 =         (/        2.       ,0.,        2.        /)
  h1 =         (/   1.  ,0.,  1.  ,0.,        2.        /)
  ! h_src =     |<-      2      ->|0|<-       2       ->|
  ! h_tgt =     |<- 1 ->|0|<- 1 ->|0|<-       2       ->|
  ! h_sub =    |0|< 1 ->|0|< 1 >|0|0|0|<-     2      ->|0|
  ! isrc_start  |1                |5|6                  |
  ! isrc_end    |            4    |5|         8         |
  ! isrc_max    |            4    |5|         8         |
  ! itgt_start  |1      |3|  4    |7|         8         |
  ! itgt_end    |   2   |3|      6|7|                  9|
  ! isub_src   |1|  1   |1|  1  |2|3|3|       3        |3|
  call intersect_src_tgt_grids( n0, h0, n1, h1, h_sub, h0_eff, &
                                isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )
  if (verbose) write(test%stdout,*) "intersect_src_tgt_grids test 6: n0=3, n1=5"
  if (verbose) write(test%stdout,*) "  h_src =     |    2    |0|    2    |"
  if (verbose) write(test%stdout,*) "  h_tgt =     | 1 |0| 1 |0|    2    |"
  call test%real_arr(9, h_sub, (/0.,1.,0.,1.,0.,0.,0.,2.,0./), 'h_sub')
  call test%real_arr(3, h0_eff, (/2.,0.,2./), 'h0_eff')
  call test%int_arr(3, isrc_start, (/1,5,6/), 'isrc_start')
  call test%int_arr(3, isrc_end, (/4,5,8/), 'isrc_end')
  call test%int_arr(3, isrc_max, (/4,5,8/), 'isrc_max')
  call test%int_arr(5, itgt_start, (/1,3,4,7,8/), 'itgt_start')
  call test%int_arr(5, itgt_end, (/2,3,6,7,9/), 'itgt_end')
  call test%int_arr(9, isub_src, (/1,1,1,1,2,3,3,3,3/), 'isub_src')
  allocate(ppoly0_coefs(n0,2), ppoly0_E(n0,2), ppoly0_S(n0,2))
  ! h_src =     |<-      2      ->|0|<-       2       ->|
  ! h_sub =    |0|< 1 ->|0|< 1 >|0|0|0|<-     2      ->|0|
  ! u_src =     |        2        |3|         4         |
  !  edge =     |1               3|3|3                 5|
  ! u_sub =    |1| 1.5  |2| 2.5 |3|3|3|       4        |5|
  call PLM_reconstruction(n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect )
  call PLM_boundary_extrapolation(n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect)
  allocate(u_sub(n0+n1+1), uh_sub(n0+n1+1))
  call remap_src_to_sub_grid_om4(n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                             n1, h_sub, h0_eff, isrc_start, isrc_end, isrc_max, isub_src, &
                             INTEGRATION_PLM, .false., u_sub, uh_sub, u02_err)
  call test%real_arr(9, u_sub, (/1.,1.5,2.,2.5,3.,3.,3.,4.,5./), 'u_sub om4')
  call remap_src_to_sub_grid(n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                             n1, h_sub, isrc_start, isrc_end, isrc_max, isub_src, &
                             INTEGRATION_PLM, .false., u_sub, uh_sub, u02_err)
  call test%real_arr(9, u_sub, (/1.,1.5,2.,2.5,3.,3.,3.,4.,5./), 'u_sub')
  ! h_sub =    |0|< 1 ->|0|< 1 >|0|0|0|<-     2      ->|0|
  ! u_sub =    |1| 1.5  |2| 2.5 |3|3|3|       4        |5|
  ! h_tgt =     |<- 1 ->|0|<- 1 ->|0|<-       2       ->|
  ! u_tgt =     |  1.5  |2|  2.5  |3|         4         |
  call remap_sub_to_tgt_grid(n0, n1, h1, h_sub, u_sub, uh_sub, itgt_start, itgt_end, &
                             .false., .false., .false., u1, u02_err)
  call test%real_arr(5, u1, (/1.5,2.,2.5,3.,4./), 'u1')
  call remap_sub_to_tgt_grid(n0, n1, h1, h_sub, u_sub, uh_sub, itgt_start, itgt_end, &
                             .true., .false., .false., u1, u02_err)
  call test%real_arr(5, u1, (/1.5,2.,2.5,3.,4./), 'u1.b')
  deallocate( ppoly0_coefs, ppoly0_E, ppoly0_S, u_sub, uh_sub, h0, u0, h1, u1)
  deallocate( h_sub, h0_eff, isrc_start, isrc_end, isrc_max, itgt_start, itgt_end, isub_src )

  ! ============================================================
  ! This section tests remapping_core_h()
  ! ============================================================
  if (verbose) write(test%stdout,*) '- - - - - - - - - - remapping_core_h() tests  - - - - - - - - -'

  allocate(u2(2))

  call initialize_remapping(CS, 'PLM', force_bounds_in_subcell=.false., answer_date=answer_date)

  ! Remapping to just the two interior layers yields the same values as u_src(2:3)
  call remapping_core_h(CS, 4, (/0.,1.,1.,0./), (/5.,4.,2.,1./), 2, (/1.,1./), u2)
  call test%real_arr(2, u2, (/4.,2./), 'PLM: remapped  h=0110->h=11 om4')
  call remapping_core_h(CS, 4, (/0.,1.,1.,0./), (/5.,4.,2.,1./), 2, (/1.,1./), u2)
  call test%real_arr(2, u2, (/4.,2./), 'PLM: remapped  h=0110->h=11')

  ! Remapping to two layers that are deeper. For the bottom layer of thickness 4,
  ! the first 1/4 has average 2, the remaining 3/4 has the bottom edge value or 1
  ! yield ing and average or 1.25
  call remapping_core_h(CS, 4, (/0.,1.,1.,0./), (/5.,4.,2.,1./), 2, (/1.,4./), u2)
  call test%real_arr(2, u2, (/4.,1.25/), 'PLM: remapped  h=0110->h=14 om4')
  call remapping_core_h(CS, 4, (/0.,1.,1.,0./), (/5.,4.,2.,1./), 2, (/1.,4./), u2)
  call test%real_arr(2, u2, (/4.,1.25/), 'PLM: remapped  h=0110->h=14')

  ! Remapping to two layers with lowest layer not reach the bottom.
  ! Here, the bottom layer samples top half of source yeilding 2.5.
  ! Note: OM4 used the value as if the target layer was the same thickness as source.
  call remapping_set_param(CS, om4_remap_via_sub_cells=.true.)
  call remapping_core_h(CS, 4, (/0.,4.,4.,0./), (/5.,4.,2.,1./), 2, (/4.,2./), u2)
  call test%real_arr(2, u2, (/4.,2./), 'PLM: remapped  h=0440->h=42 om4 (with known bug)')
  call remapping_set_param(CS, om4_remap_via_sub_cells=.false.)
  call remapping_core_h(CS, 4, (/0.,4.,4.,0./), (/5.,4.,2.,1./), 2, (/4.,2./), u2)
  call test%real_arr(2, u2, (/4.,2.5/), 'PLM: remapped  h=0440->h=42')

  ! Remapping to two layers with no layers sampling the bottom source layer
  ! The first layer samples the top half of u1, yielding 4.5
  ! The second layer samples the next quarter of u1, yielding 3.75
  call remapping_set_param(CS, om4_remap_via_sub_cells=.true.)
  call remapping_core_h(CS, 4, (/0.,5.,5.,0./), (/5.,4.,2.,1./), 2, (/2.,2./), u2)
  call test%real_arr(2, u2, (/4.5,3.5/), 'PLM: remapped  h=0880->h=21 om4 (with known bug)')
  call remapping_set_param(CS, om4_remap_via_sub_cells=.false.)
  call remapping_core_h(CS, 4, (/0.,4.,4.,0./), (/5.,4.,2.,1./), 2, (/2.,1./), u2)
  call test%real_arr(2, u2, (/4.5,3.75/), 'PLM: remapped  h=0440->h=21')

  deallocate(u2)

  ! Profile 0: 8 layers, 1x top/2x bottom vanished, and the rest with thickness 1.0, total depth 5, u(z) = 1 + z
  n0 = 8
  allocate( h0(n0), u0(n0) )
  h0 = (/0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0/)
  u0 = (/1.0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.0, 6.0/)
  allocate( u1(8) )

  call initialize_remapping(CS, 'PLM', answer_date=99990101, h_neglect=1.e-17, h_neglect_edge=1.e-2)

  do om4 = 0, 1
    if ( om4 == 0 ) then
      CS%om4_remap_via_sub_cells = .false.
      om4_tag(:) = '    '
    else
      CS%om4_remap_via_sub_cells = .true.
      om4_tag(:) = ' om4'
    endif

    ! Unchanged grid
    call remapping_core_h( CS, n0, h0, u0, 8, [0.,1.,1.,1.,1.,1.,0.,0.], u1)
    call test%real_arr(8, u1, (/1.0,1.5,2.5,3.5,4.5,5.5,6.0,6.0/), 'PLM: remapped  h=01111100->h=01111100'//om4_tag)

    ! Removing vanished layers (unchanged values for non-vanished layers, layer centers 0.5, 1.5, 2.5, 3.5, 4.5)
    call remapping_core_h( CS, n0, h0, u0, 5, [1.,1.,1.,1.,1.], u1)
    call test%real_arr(5, u1, (/1.5,2.5,3.5,4.5,5.5/), 'PLM: remapped  h=01111100->h=11111'//om4_tag)

    ! Remapping to variable thickness layers (layer centers 0.25, 1.0, 2.25, 4.0)
    call remapping_core_h( CS, n0, h0, u0, 4, [0.5,1.,1.5,2.], u1)
    call test%real_arr(4, u1, (/1.25,2.,3.25,5./), 'PLM: remapped  h=01111100->h=h1t2'//om4_tag)

    ! Remapping to variable thickness + vanished layers (layer centers 0.25, 1.0, 1.5, 2.25, 4.0)
    call remapping_core_h( CS, n0, h0, u0, 6, [0.5,1.,0.,1.5,2.,0.], u1)
    call test%real_arr(6, u1, (/1.25,2.,2.5,3.25,5.,6./), 'PLM: remapped  h=01111100->h=h10t20'//om4_tag)

    ! Remapping to deeper water column (layer centers 0.75, 2.25, 3., 5., 8.)
    call remapping_core_h( CS, n0, h0, u0, 5, [1.5,1.5,0.,4.,2.], u1)
    call test%real_arr(5, u1, (/1.75,3.25,4.,5.5,6./), 'PLM: remapped  h=01111100->h=tt02'//om4_tag)

    ! Remapping to slightly shorter water column (layer centers 0.5, 1.5, 2.5,, 3.5, 4.25)
    call remapping_core_h( CS, n0, h0, u0, 5, [1.,1.,1.,1.,0.5], u1)
    if ( om4 == 0 ) then
      call test%real_arr(5, u1, (/1.5,2.5,3.5,4.5,5.25/), 'PLM: remapped  h=01111100->h=1111h')
    else
      call test%real_arr(5, u1, (/1.5,2.5,3.5,4.5,5.5/), 'PLM: remapped  h=01111100->h=1111h om4 (known bug)')
    endif

    ! Remapping to much shorter water column (layer centers 0.25, 0.5, 1.)
    call remapping_core_h( CS, n0, h0, u0, 3, [0.5,0.,1.], u1)
    if ( om4 == 0 ) then
      call test%real_arr(3, u1, (/1.25,1.5,2./), 'PLM: remapped  h=01111100->h=h01')
    else
      call test%real_arr(3, u1, (/1.25,1.5,1.875/), 'PLM: remapped  h=01111100->h=h01 om4 (known bug)')
    endif

  enddo ! om4

  call end_remapping(CS)
  deallocate( h0, u0, u1 )

  ! ============================================================
  ! This section tests interpolation and reintegration functions
  ! ============================================================
  if (verbose) write(test%stdout,*) '- - - - - - - - - - interpolation tests  - - - - - - - - -'

  call test_interp(test, 'Identity: 3 layer', &
                     3, (/1.,2.,3./), (/1.,2.,3.,4./), &
                     3, (/1.,2.,3./), (/1.,2.,3.,4./) )

  call test_interp(test, 'A: 3 layer to 2', &
                     3, (/1.,1.,1./), (/1.,2.,3.,4./), &
                     2, (/1.5,1.5/), (/1.,2.5,4./) )

  call test_interp(test, 'B: 2 layer to 3', &
                     2, (/1.5,1.5/), (/1.,4.,7./), &
                     3, (/1.,1.,1./), (/1.,3.,5.,7./) )

  call test_interp(test, 'C: 3 layer (vanished middle) to 2', &
                     3, (/1.,0.,2./), (/1.,2.,2.,3./), &
                     2, (/1.,2./), (/1.,2.,3./) )

  call test_interp(test, 'D: 3 layer (deep) to 3', &
                     3, (/1.,2.,3./), (/1.,2.,4.,7./), &
                     2, (/2.,2./), (/1.,3.,5./) )

  call test_interp(test, 'E: 3 layer to 3 (deep)', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     3, (/2.,3.,4./), (/1.,3.,6.,8./) )

  call test_interp(test, 'F: 3 layer to 4 with vanished top/botton', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     4, (/0.,2.,5.,0./), (/0.,1.,3.,8.,0./) )

  call test_interp(test, 'Fs: 3 layer to 4 with vanished top/botton (shallow)', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     4, (/0.,2.,4.,0./), (/0.,1.,3.,7.,0./) )

  call test_interp(test, 'Fd: 3 layer to 4 with vanished top/botton (deep)', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     4, (/0.,2.,6.,0./), (/0.,1.,3.,8.,0./) )

  if (verbose) write(test%stdout,*) '  - - - - - reintegration tests - - - - -'

  call test_reintegrate(test, 'Identity: 3 layer', &
                     3, (/1.,2.,3./), (/-5.,2.,1./), &
                     3, (/1.,2.,3./), (/-5.,2.,1./) )

  call test_reintegrate(test, 'A: 3 layer to 2', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     2, (/3.,3./), (/-4.,2./) )

  call test_reintegrate(test, 'A: 3 layer to 2 (deep)', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     2, (/3.,4./), (/-4.,2./) )

  call test_reintegrate(test, 'A: 3 layer to 2 (shallow)', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     2, (/3.,2./), (/-4.,1.5/) )

  call test_reintegrate(test, 'B: 3 layer to 4 with vanished top/bottom', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     4, (/0.,3.,3.,0./), (/0.,-4.,2.,0./) )

  call test_reintegrate(test, 'C: 3 layer to 4 with vanished top//middle/bottom', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     5, (/0.,3.,0.,3.,0./), (/0.,-4.,0.,2.,0./) )

  call test_reintegrate(test, 'D: 3 layer to 3 (vanished)', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     3, (/0.,0.,0./), (/0.,0.,0./) )

  call test_reintegrate(test, 'D: 3 layer (vanished) to 3', &
                     3, (/0.,0.,0./), (/-5.,2.,1./), &
                     3, (/2.,2.,2./), (/0., 0., 0./) )

  call test_reintegrate(test, 'D: 3 layer (vanished) to 3 (vanished)', &
                     3, (/0.,0.,0./), (/-5.,2.,1./), &
                     3, (/0.,0.,0./), (/0.,0.,0./) )

  call test_reintegrate(test, 'D: 3 layer (vanished) to 3 (vanished)', &
                     3, (/0.,0.,0./), (/0.,0.,0./), &
                     3, (/0.,0.,0./), (/0.,0.,0./) )

  if (verbose) write(test%stdout,*) '- - - - - - - - - - Recon1d PCM tests  - - - - - - - - -'
  call test%test( PCM%unit_tests(verbose, test%stdout, test%stderr), 'PCM unit test')
  call test%test( MPLM_WA%unit_tests(verbose, test%stdout, test%stderr), 'MPLM_WA unit test')
  call test%test( EMPLM_WA%unit_tests(verbose, test%stdout, test%stderr), 'EMPLM_WA unit test')
  call test%test( MPLM_WA_poly%unit_tests(verbose, test%stdout, test%stderr), 'MPLM_WA_poly unit test')
  call test%test( EMPLM_WA_poly%unit_tests(verbose, test%stdout, test%stderr), 'EMPLM_WA_poly unit test')
  call test%test( PLM_hybgen%unit_tests(verbose, test%stdout, test%stderr), 'PLM_hybgen unit test')
  call test%test( PLM_CW%unit_tests(verbose, test%stdout, test%stderr), 'PLM_CW unit test')
  call test%test( PLM_CWK%unit_tests(verbose, test%stdout, test%stderr), 'PLM_CWK unit test')
  call test%test( MPLM_CWK%unit_tests(verbose, test%stdout, test%stderr), 'MPLM_CWK unit test')
  call test%test( EMPLM_CWK%unit_tests(verbose, test%stdout, test%stderr), 'EMPLM_CWK unit test')
  call test%test( PPM_H4_2019%unit_tests(verbose, test%stdout, test%stderr), 'PPM_H4_2019 unit test')
  call test%test( PPM_H4_2018%unit_tests(verbose, test%stdout, test%stderr), 'PPM_H4_2018 unit test')
  call test%test( PPM_hybgen%unit_tests(verbose, test%stdout, test%stderr), 'PPM_hybgen unit test')
  call test%test( PPM_CW%unit_tests(verbose, test%stdout, test%stderr), 'PPM_CW unit test')
  call test%test( PPM_CWK%unit_tests(verbose, test%stdout, test%stderr), 'PPM_CWK unit test')
  call test%test( EPPM_CWK%unit_tests(verbose, test%stdout, test%stderr), 'EPPM_CWK unit test')

  ! Randomized, brute force tests
  ntests = 3000
  if (present(num_comp_samp)) ntests = num_comp_samp

  call random_seed(size=seed_size)
  allocate( seed(seed_Size) )
  seed(:) = 102030405
  call random_seed(put=seed)

  n0 = 9

  ! Internal consistency
  call test_recon_consistency(test, 'C_PCM', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_PLM_CW', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_PLM_HYBGEN', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_MPLM_WA', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_EMPLM_WA', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_MPLM_WA_poly', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_EMPLM_WA_poly', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_PLM_CWK', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_MPLM_CWK', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_EMPLM_CWK', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_PPM_H4_2018', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_PPM_H4_2019', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_PPM_HYBGEN', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_PPM_CW', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_PPM_CWK', n0, ntests, h_neglect)
  call test_recon_consistency(test, 'C_EPPM_CWK', n0, ntests, h_neglect)

  call test_preserve_uniform(test, 'PCM', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_PCM', n0, ntests, h_neglect)
! call test_preserve_uniform(test, 'PLM', n0, ntests, h_neglect) ! Fails
! call test_preserve_uniform(test, 'PLM_HYBGEN', n0, ntests, h_neglect) ! Fails
! call test_preserve_uniform(test, 'PPM_H4', n0, ntests, h_neglect) ! Fails
! call test_preserve_uniform(test, 'PPM_IH4', n0, ntests, h_neglect) ! Fails
! call test_preserve_uniform(test, 'PPM_HYBGEN', n0, ntests, h_neglect) ! Fails
! call test_preserve_uniform(test, 'PPM_CW', n0, ntests, h_neglect) ! Fails
! call test_preserve_uniform(test, 'WENO_HYBGEN', n0, ntests, h_neglect) ! Fails
! call test_preserve_uniform(test, 'PQM_IH4IH3', n0, ntests, h_neglect) ! Fails
  call test_preserve_uniform(test, 'C_PLM_CW', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_PLM_HYBGEN', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_MPLM_WA', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_EMPLM_WA', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_MPLM_WA_poly', n0, ntests, h_neglect) ! Surprised this passes -AJA
! call test_preserve_uniform(test, 'C_EMPLM_WA_poly', n0, ntests, h_neglect) ! This is known to fail
  call test_preserve_uniform(test, 'C_PLM_CWK', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_MPLM_CWK', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_EMPLM_CWK', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_PPM_H4_2019', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_PPM_H4_2018', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_PPM_HYBGEN', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_PPM_CW', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_PPM_CWK', n0, ntests, h_neglect)
  call test_preserve_uniform(test, 'C_EPPM_CWK', n0, ntests, h_neglect)

  call test_unchanged_grid(test, 'C_PCM', n0, ntests, h_neglect)
  call test_unchanged_grid(test, 'C_PLM_CW', n0, ntests, h_neglect)
  call test_unchanged_grid(test, 'C_PLM_HYBGEN', n0, ntests, h_neglect)
  call test_unchanged_grid(test, 'C_PLM_CWK', n0, ntests, h_neglect)
  call test_unchanged_grid(test, 'C_MPLM_CWK', n0, ntests, h_neglect)
  call test_unchanged_grid(test, 'C_EMPLM_CWK', n0, ntests, h_neglect)
  call test_unchanged_grid(test, 'C_PPM_HYBGEN', n0, ntests, h_neglect)
  call test_unchanged_grid(test, 'C_PPM_CW', n0, ntests, h_neglect)
  call test_unchanged_grid(test, 'C_PPM_CWK', n0, ntests, h_neglect)
  call test_unchanged_grid(test, 'C_EPPM_CWK', n0, ntests, h_neglect)

  ! Check that remapping to the exact same grid leaves values unchanged
  allocate( h0(8), u0(8) )
  h0 = (/0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0/)
  u0 = (/1.0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.0, 6.0/)
  allocate( u1(8) )
  call initialize_remapping(CS, 'C_PLM_CW', nk=8)
  call remapping_core_h( CS, 8, h0, u0, 8, [0.,1.,1.,1.,1.,1.,0.,0.], u1 )
  call test%real_arr(8, u1, u0, 'remapping_core to unchanged grid with class')

  call end_remapping(CS)
  deallocate( h0, u0, u1 )

  ! Brute force test that we have bitwise identical answers with the new classes
  n0 = 7
  n1 = 4

  ! PPM_CW and PPM_HYBGEN are identical, but are different options in build_reconstructions_1d()
  call initialize_remapping(CS, 'PPM_CW', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.false., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'PPM_HYBGEN', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.false., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PPM_CW <-> PPM_HYBGEN')

  ! PPM_CW <-> PPM_HYBGEN, as above but with OM4 subcells
  call initialize_remapping(CS, 'PPM_CW', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.true., force_bounds_in_subcell=.false., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'PPM_HYBGEN', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.true., force_bounds_in_subcell=.false., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PPM_CW <-> PPM_HYBGEN OM4')

  ! PPM_CW <-> PPM_HYBGEN, as above but with extrapolation
  call initialize_remapping(CS, 'PPM_CW', answer_date=99990101, boundary_extrapolation=.true., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.false., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'PPM_HYBGEN', answer_date=99990101, boundary_extrapolation=.true., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.false., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PPM_CW <-> PPM_HYBGEN Extrap')

  ! PPM_CW <-> PPM_HYBGEN, as above but with OM4 subcells and subcell bounds
  call initialize_remapping(CS, 'PPM_CW', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.true., force_bounds_in_subcell=.true., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'PPM_HYBGEN', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.true., force_bounds_in_subcell=.true., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PPM_CW <-> PPM_HYBGEN')

  ! PCM <-> C_PCM
  call initialize_remapping(CS, 'PCM', answer_date=99990101, om4_remap_via_sub_cells=.false., &
                            force_bounds_in_subcell=.false.)
  call initialize_remapping(CS2, 'C_PCM', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PCM <-> C_PCM')

  ! PLM <-> C_MPLM_WA_POLY
  call initialize_remapping(CS, 'PLM', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.false., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'C_MPLM_WA_POLY', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PLM <-> C_MPLM_WA_poly')

  ! PLM (with subcell bounds) <-> C_MPLM_WA_POLY
  call initialize_remapping(CS, 'PLM', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.true., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'C_MPLM_WA_POLY', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PLM bounded <-> C_MPLM_WA_poly')

  ! PLM + extrapolation <-> C_EMPLM_WA_POLY
  call initialize_remapping(CS, 'PLM', answer_date=99990101, boundary_extrapolation=.true., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.false., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'C_EMPLM_WA_POLY', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PLM <-> C_EMPLM_WA_poly')

  ! PLM + extrapolation (with subcell bounds) <-> C_EMPLM_WA_POLY
  call initialize_remapping(CS, 'PLM', answer_date=99990101, boundary_extrapolation=.true., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.true., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'C_EMPLM_WA_POLY', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PLM bounded <-> C_EMPLM_WA_poly')

  ! PPM_H4 (2018 answers) <-> C_PPM_H4_2018
  call initialize_remapping(CS, 'PPM_H4', answer_date=20180101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.false., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'C_PPM_H4_2018', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PPM_H4 2018 <-> C_PPM_H4_2018')

  ! PPM_H4 (2018 answers with subcell bounds) <-> C_PPM_H4_2018
  call initialize_remapping(CS, 'PPM_H4', answer_date=20180101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.true., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'C_PPM_H4_2018', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PPM_H4 2018 bounded <-> C_PPM_H4_2018')

  ! PPM_H4 (latest answers) <-> C_PPM_H4_2019
  call initialize_remapping(CS, 'PPM_H4', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.false., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'C_PPM_H4_2019', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PPM_H4 <-> C_PPM_H4_2019')

  ! PPM_H4 (latest answers with subcell bounds) <-> C_PPM_H4_2019
  call initialize_remapping(CS, 'PPM_H4', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.true., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'C_PPM_H4_2019', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PPM_H4 bounded <-> C_PPM_H4_2019')

  ! PLM_HYBGEN (latest answers with subcell bounds) <-> C_PLM_hybgen
  call initialize_remapping(CS, 'PLM_HYBGEN', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.true., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'C_PLM_hybgen', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PLM_HYBGEN bounded <-> C_PLM_hygen')

  ! PPM_HYBGEN (latest answers with subcell bounds) <-> C_PPM_hybgen
  call initialize_remapping(CS, 'PPM_HYBGEN', answer_date=99990101, boundary_extrapolation=.false., &
                            om4_remap_via_sub_cells=.false., force_bounds_in_subcell=.true., &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
  call initialize_remapping(CS2, 'C_PPM_HYBGEN', nk=n0, h_neglect=h_neglect)
  call compare_two_schemes(test, CS, CS2, n0, n1, ntests, 'PPM_HYBGEN bounded <-> C_PPM_HYGEN')

  call end_remapping(CS)
  call end_remapping(CS2)

  remapping_unit_tests = test%summarize('remapping_unit_tests')

end function remapping_unit_tests

end module MOM_remapping
