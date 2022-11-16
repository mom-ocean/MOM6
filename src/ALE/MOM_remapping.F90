!> Provides column-wise vertical remapping functions
module MOM_remapping

! This file is part of MOM6. See LICENSE.md for the license.
! Original module written by Laurent White, 2008.06.09

use MOM_error_handler, only : MOM_error, FATAL
use MOM_io,            only : stdout, stderr
use MOM_string_functions, only : uppercase
use regrid_edge_values, only : edge_values_explicit_h4, edge_values_implicit_h4
use regrid_edge_values, only : edge_values_implicit_h4, edge_values_implicit_h6
use regrid_edge_values, only : edge_slopes_implicit_h3, edge_slopes_implicit_h5
use remapping_attic,    only : remapping_attic_unit_tests
use PCM_functions, only : PCM_reconstruction
use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation
use PQM_functions, only : PQM_reconstruction, PQM_boundary_extrapolation_v1
use MOM_hybgen_remap, only : hybgen_plm_coefs, hybgen_ppm_coefs, hybgen_weno_coefs

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
  !> The vintage of the expressions to use for remapping. Values below 20190101 result
  !! in the use of older, less accurate expressions.
  integer :: answer_date = 99991231
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
integer, parameter  :: REMAPPING_PPM_H4     = 4 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PPM_IH4    = 5 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PPM_HYBGEN = 6 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_WENO_HYBGEN= 7 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PQM_IH4IH3 = 8 !< O(h^4) remapping scheme
integer, parameter  :: REMAPPING_PQM_IH6IH5 = 9 !< O(h^5) remapping scheme

integer, parameter  :: INTEGRATION_PCM = 0  !< Piecewise Constant Method
integer, parameter  :: INTEGRATION_PLM = 1  !< Piecewise Linear Method
integer, parameter  :: INTEGRATION_PPM = 3  !< Piecewise Parabolic Method
integer, parameter  :: INTEGRATION_PQM = 5  !< Piecewise Quartic Method

character(len=40)  :: mdl = "MOM_remapping" !< This module's name.

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
               check_reconstruction, check_remapping, force_bounds_in_subcell, answers_2018, answer_date)
  type(remapping_CS),         intent(inout) :: CS !< Remapping control structure
  character(len=*), optional, intent(in)    :: remapping_scheme !< Remapping scheme to use
  logical, optional,          intent(in)    :: boundary_extrapolation !< Indicate to extrapolate in boundary cells
  logical, optional,          intent(in)    :: check_reconstruction !< Indicate to check reconstructions
  logical, optional,          intent(in)    :: check_remapping !< Indicate to check results of remapping
  logical, optional,          intent(in)    :: force_bounds_in_subcell !< Force subcells values to be bounded
  logical, optional,          intent(in)    :: answers_2018 !< If true use older, less accurate expressions.
  integer, optional,          intent(in)    :: answer_date  !< The vintage of the expressions to use

  if (present(remapping_scheme)) then
    call setReconstructionType( remapping_scheme, CS )
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

end subroutine remapping_set_param

subroutine extract_member_remapping_CS(CS, remapping_scheme, degree, boundary_extrapolation, check_reconstruction, &
                                       check_remapping, force_bounds_in_subcell)
  type(remapping_CS), intent(in) :: CS !< Control structure for remapping module
  integer, optional, intent(out) :: remapping_scheme        !< Determines which reconstruction scheme to use
  integer, optional, intent(out) :: degree                  !< Degree of polynomial reconstruction
  logical, optional, intent(out) :: boundary_extrapolation  !< If true, extrapolate boundaries
  logical, optional, intent(out) :: check_reconstruction    !< If true, reconstructions are checked for consistency.
  logical, optional, intent(out) :: check_remapping         !< If true, the result of remapping are checked
                                                            !!  for conservation and bounds.
  logical, optional, intent(out) :: force_bounds_in_subcell !< If true, the intermediate values used in
                                                            !! remapping are forced to be bounded.

  if (present(remapping_scheme)) remapping_scheme = CS%remapping_scheme
  if (present(degree)) degree = CS%degree
  if (present(boundary_extrapolation)) boundary_extrapolation = CS%boundary_extrapolation
  if (present(check_reconstruction)) check_reconstruction = CS%check_reconstruction
  if (present(check_remapping)) check_remapping = CS%check_remapping
  if (present(force_bounds_in_subcell)) force_bounds_in_subcell = CS%force_bounds_in_subcell

end subroutine extract_member_remapping_CS

!> Calculate edge coordinate x from cell width h
subroutine buildGridFromH(nz, h, x)
  integer,               intent(in)    :: nz !< Number of cells
  real, dimension(nz),   intent(in)    :: h  !< Cell widths [H]
  real, dimension(nz+1), intent(inout) :: x  !< Edge coordinates starting at x(1)=0 [H]
  ! Local variables
  integer :: k

  x(1) = 0.0
  do k = 1,nz
    x(k+1) = x(k) + h(k)
  enddo

end subroutine buildGridFromH

!> Remaps column of values u0 on grid h0 to grid h1 assuming the top edge is aligned.
subroutine remapping_core_h(CS, n0, h0, u0, n1, h1, u1, h_neglect, h_neglect_edge, PCM_cell)
  type(remapping_CS),  intent(in)  :: CS !< Remapping control structure
  integer,             intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0), intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0), intent(in)  :: u0 !< Cell averages on source grid [A]
  integer,             intent(in)  :: n1 !< Number of cells on target grid
  real, dimension(n1), intent(in)  :: h1 !< Cell widths on target grid [H]
  real, dimension(n1), intent(out) :: u1 !< Cell averages on target grid [A]
  real, optional,      intent(in)  :: h_neglect !< A negligibly small width for the
                                         !! purpose of cell reconstructions
                                         !! in the same units as h0 [H]
  real, optional,      intent(in)  :: h_neglect_edge !< A negligibly small width
                                         !! for the purpose of edge value
                                         !! calculations in the same units as h0 [H]
  logical, dimension(n0), optional, intent(in) :: PCM_cell !< If present, use PCM remapping for
                                         !! cells in the source grid where this is true.

  ! Local variables
  real, dimension(n0,2)           :: ppoly_r_E     ! Edge value of polynomial [A]
  real, dimension(n0,2)           :: ppoly_r_S     ! Edge slope of polynomial [A H-1]
  real, dimension(n0,CS%degree+1) :: ppoly_r_coefs ! Coefficients of polynomial reconstructions [A]
  real :: uh_err       ! A bound on the error in the sum of u*h, as estimated by the remapping code [A H]
  real :: hNeglect, hNeglect_edge ! Negligibly small cell widths in the same units as h0 [H]
  integer :: iMethod   ! An integer indicating the integration method used
  integer :: k

  hNeglect = 1.0e-30 ; if (present(h_neglect)) hNeglect = h_neglect
  hNeglect_edge = 1.0e-10 ; if (present(h_neglect_edge)) hNeglect_edge = h_neglect_edge

  call build_reconstructions_1d( CS, n0, h0, u0, ppoly_r_coefs, ppoly_r_E, ppoly_r_S, iMethod, &
                               hNeglect, hNeglect_edge, PCM_cell )

  if (CS%check_reconstruction) call check_reconstructions_1d(n0, h0, u0, CS%degree, &
                                 CS%boundary_extrapolation, ppoly_r_coefs, ppoly_r_E, ppoly_r_S)

  call remap_via_sub_cells( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, n1, h1, iMethod, &
                          CS%force_bounds_in_subcell, u1, uh_err )

  if (CS%check_remapping) call check_remapped_values(n0, h0, u0, ppoly_r_E, CS%degree, ppoly_r_coefs, &
                                                     n1, h1, u1, iMethod, uh_err, "remapping_core_h")

end subroutine remapping_core_h

!> Remaps column of values u0 on grid h0 to implied grid h1
!! where the interfaces of h1 differ from those of h0 by dx.
subroutine remapping_core_w( CS, n0, h0, u0, n1, dx, u1, h_neglect, h_neglect_edge )
  type(remapping_CS),    intent(in)  :: CS !< Remapping control structure
  integer,               intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),   intent(in)  :: u0 !< Cell averages on source grid [A]
  integer,               intent(in)  :: n1 !< Number of cells on target grid
  real, dimension(n1+1), intent(in)  :: dx !< Cell widths on target grid [H]
  real, dimension(n1),   intent(out) :: u1 !< Cell averages on target grid [A]
  real, optional,        intent(in)  :: h_neglect !< A negligibly small width for the
                                           !! purpose of cell reconstructions
                                           !! in the same units as h0 [H].
  real, optional,        intent(in)  :: h_neglect_edge !< A negligibly small width
                                           !! for the purpose of edge value
                                           !! calculations in the same units as h0 [H].
  ! Local variables
  real, dimension(n0,2)           :: ppoly_r_E     ! Edge value of polynomial [A]
  real, dimension(n0,2)           :: ppoly_r_S     ! Edge slope of polynomial [A H-1]
  real, dimension(n0,CS%degree+1) :: ppoly_r_coefs ! Coefficients of polynomial reconstructions [A]
  real, dimension(n1) :: h1 !< Cell widths on target grid [H]
  real :: uh_err       ! A bound on the error in the sum of u*h, as estimated by the remapping code [A H]
  real :: hNeglect, hNeglect_edge  ! Negligibly small thicknesses [H]
  integer :: iMethod   ! An integer indicating the integration method used
  integer :: k

  hNeglect = 1.0e-30 ; if (present(h_neglect)) hNeglect = h_neglect
  hNeglect_edge = 1.0e-10 ; if (present(h_neglect_edge)) hNeglect_edge = h_neglect_edge

  call build_reconstructions_1d( CS, n0, h0, u0, ppoly_r_coefs, ppoly_r_E, ppoly_r_S, iMethod,&
                                  hNeglect, hNeglect_edge )

  if (CS%check_reconstruction) call check_reconstructions_1d(n0, h0, u0, CS%degree, &
                                   CS%boundary_extrapolation, ppoly_r_coefs, ppoly_r_E, ppoly_r_S)

  ! This is a temporary step prior to switching to remapping_core_h()
  do k = 1, n1
    if (k<=n0) then
      h1(k) = max( 0., h0(k) + ( dx(k+1) - dx(k) ) )
    else
      h1(k) = max( 0., dx(k+1) - dx(k) )
    endif
  enddo
  call remap_via_sub_cells( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, n1, h1, iMethod, &
                            CS%force_bounds_in_subcell, u1, uh_err )
! call remapByDeltaZ( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, n1, dx, iMethod, u1, hNeglect )
! call remapByProjection( n0, h0, u0, CS%ppoly_r, n1, h1, iMethod, u1, hNeglect )

  if (CS%check_remapping) call check_remapped_values(n0, h0, u0, ppoly_r_E, CS%degree, ppoly_r_coefs, &
                                                     n1, h1, u1, iMethod, uh_err, "remapping_core_w")

end subroutine remapping_core_w

!> Creates polynomial reconstructions of u0 on the source grid h0.
subroutine build_reconstructions_1d( CS, n0, h0, u0, ppoly_r_coefs, &
                                     ppoly_r_E, ppoly_r_S, iMethod, h_neglect, &
                                     h_neglect_edge, PCM_cell )
  type(remapping_CS),    intent(in)  :: CS !< Remapping control structure
  integer,               intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),   intent(in)  :: u0 !< Cell averages on source grid [A]
  real, dimension(n0,CS%degree+1), &
                         intent(out) :: ppoly_r_coefs !< Coefficients of polynomial [A]
  real, dimension(n0,2), intent(out) :: ppoly_r_E !< Edge value of polynomial [A]
  real, dimension(n0,2), intent(out) :: ppoly_r_S !< Edge slope of polynomial [A H-1]
  integer,               intent(out) :: iMethod !< Integration method
  real, optional,        intent(in)  :: h_neglect !< A negligibly small width for the
                                         !! purpose of cell reconstructions
                                         !! in the same units as h0 [H]
  real, optional,        intent(in)  :: h_neglect_edge !< A negligibly small width
                                         !! for the purpose of edge value
                                         !! calculations in the same units as h0 [H]
  logical, optional,     intent(in)  :: PCM_cell(n0) !< If present, use PCM remapping for
                                         !! cells from the source grid where this is true.

  ! Local variables
  integer :: local_remapping_scheme
  integer :: k, n

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
  elseif (n0<=4) then
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
    case ( REMAPPING_PPM_H4 )
      call edge_values_explicit_h4( n0, h0, u0, ppoly_r_E, h_neglect_edge, answer_date=CS%answer_date )
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect, answer_date=CS%answer_date )
      if ( CS%boundary_extrapolation ) then
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefs, h_neglect )
      endif
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PPM_IH4 )
      call edge_values_implicit_h4( n0, h0, u0, ppoly_r_E, h_neglect_edge, answer_date=CS%answer_date )
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
      call edge_values_implicit_h4( n0, h0, u0, ppoly_r_E, h_neglect_edge, answer_date=CS%answer_date )
      call edge_slopes_implicit_h3( n0, h0, u0, ppoly_r_S, h_neglect, answer_date=CS%answer_date )
      call PQM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefs, h_neglect, &
                               answer_date=CS%answer_date )
      if ( CS%boundary_extrapolation ) then
        call PQM_boundary_extrapolation_v1( n0, h0, u0, ppoly_r_E, ppoly_r_S, &
                                            ppoly_r_coefs, h_neglect )
      endif
      iMethod = INTEGRATION_PQM
    case ( REMAPPING_PQM_IH6IH5 )
      call edge_values_implicit_h6( n0, h0, u0, ppoly_r_E, h_neglect_edge, answer_date=CS%answer_date )
      call edge_slopes_implicit_h5( n0, h0, u0, ppoly_r_S, h_neglect, answer_date=CS%answer_date )
      call PQM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefs, h_neglect, &
                               answer_date=CS%answer_date )
      if ( CS%boundary_extrapolation ) then
        call PQM_boundary_extrapolation_v1( n0, h0, u0, ppoly_r_E, ppoly_r_S, &
                                            ppoly_r_coefs, h_neglect )
      endif
      iMethod = INTEGRATION_PQM
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
                                    ppoly_r_coefs, ppoly_r_E, ppoly_r_S)
  integer,                  intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),      intent(in)  :: h0 !< Cell widths on source grid [H]
  real, dimension(n0),      intent(in)  :: u0 !< Cell averages on source grid [A]
  integer,                  intent(in)  :: deg !< Degree of polynomial reconstruction
  logical,                  intent(in)  :: boundary_extrapolation !< Extrapolate at boundaries if true
  real, dimension(n0,deg+1),intent(out) :: ppoly_r_coefs !< Coefficients of polynomial [A]
  real, dimension(n0,2),    intent(out) :: ppoly_r_E !< Edge value of polynomial [A]
  real, dimension(n0,2),    intent(out) :: ppoly_r_S !< Edge slope of polynomial [A H-1]
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

!> Remaps column of n0 values u0 on grid h0 to grid h1 with n1 cells by calculating
!! the n0+n1+1 sub-integrals of the intersection of h0 and h1, and the summing the
!! appropriate integrals into the h1*u1 values.  h0 and h1 must have the same units.
subroutine remap_via_sub_cells( n0, h0, u0, ppoly0_E, ppoly0_coefs, n1, h1, method, &
                                force_bounds_in_subcell, u1, uh_err, ah_sub, aisub_src, aiss, aise )
  integer,           intent(in)    :: n0      !< Number of cells in source grid
  real,              intent(in)    :: h0(n0)  !< Source grid widths (size n0) [H]
  real,              intent(in)    :: u0(n0)  !< Source cell averages (size n0) [A]
  real,              intent(in)    :: ppoly0_E(n0,2)    !< Edge value of polynomial [A]
  real,              intent(in)    :: ppoly0_coefs(:,:) !< Coefficients of polynomial [A]
  integer,           intent(in)    :: n1      !< Number of cells in target grid
  real,              intent(in)    :: h1(n1)  !< Target grid widths (size n1) [H]
  integer,           intent(in)    :: method  !< Remapping scheme to use
  logical,           intent(in)    :: force_bounds_in_subcell !< Force sub-cell values to be bounded
  real,              intent(out)   :: u1(n1)  !< Target cell averages (size n1) [A]
  real,              intent(out)   :: uh_err  !< Estimate of bound on error in sum of u*h [A H]
  real, optional,    intent(out)   :: ah_sub(n0+n1+1) !< Overlapping sub-cell thicknesses, h_sub [H]
  integer, optional, intent(out)   :: aisub_src(n0+n1+1) !< i_sub_src
  integer, optional, intent(out)   :: aiss(n0) !< isrc_start
  integer, optional, intent(out)   :: aise(n0) !< isrc_ens
  ! Local variables
  integer :: i_sub ! Index of sub-cell
  integer :: i0 ! Index into h0(1:n0), source column
  integer :: i1 ! Index into h1(1:n1), target column
  integer :: i_start0 ! Used to record which sub-cells map to source cells
  integer :: i_start1 ! Used to record which sub-cells map to target cells
  integer :: i_max ! Used to record which sub-cell is the largest contribution of a source cell
  real :: dh_max ! Used to record which sub-cell is the largest contribution of a source cell [H]
  real, dimension(n0+n1+1) :: h_sub ! Width of each each sub-cell [H]
  real, dimension(n0+n1+1) :: uh_sub ! Integral of u*h over each sub-cell [A H]
  real, dimension(n0+n1+1) :: u_sub ! Average of u over each sub-cell [A]
  integer, dimension(n0+n1+1) :: isub_src ! Index of source cell for each sub-cell
  integer, dimension(n0) :: isrc_start ! Index of first sub-cell within each source cell
  integer, dimension(n0) :: isrc_end ! Index of last sub-cell within each source cell
  integer, dimension(n0) :: isrc_max ! Index of thickest sub-cell within each source cell
  real, dimension(n0) :: h0_eff ! Effective thickness of source cells [H]
  real, dimension(n0) :: u0_min ! Minimum value of reconstructions in source cell [A]
  real, dimension(n0) :: u0_max ! Minimum value of reconstructions in source cell [A]
  integer, dimension(n1) :: itgt_start ! Index of first sub-cell within each target cell
  integer, dimension(n1) :: itgt_end ! Index of last sub-cell within each target cell
  real :: xa, xb ! Non-dimensional position within a source cell (0..1) [nondim]
  real :: h0_supply, h1_supply ! The amount of width available for constructing sub-cells [H]
  real :: dh ! The width of the sub-cell [H]
  real :: duh ! The total amount of accumulated stuff (u*h)  [A H]
  real :: dh0_eff ! Running sum of source cell thickness [H]
  ! For error checking/debugging
  logical, parameter :: force_bounds_in_target = .true. ! To fix round-off issues
  logical, parameter :: adjust_thickest_subcell = .true. ! To fix round-off conservation issues
  logical, parameter :: debug_bounds = .false. ! For debugging overshoots etc.
  integer :: k, i0_last_thick_cell
  real :: h0tot, h1tot, h2tot ! Summed thicknesses used for debugging [H]
  real :: h0err, h1err, h2err ! Estimates of round-off errors used for debugging [H]
  real :: u02_err, u0err, u1err, u2err  ! Integrated reconstruction error estimates [H A]
  real :: u0tot, u1tot, u2tot ! Integrated reconstruction values [H A]
  real :: u_orig              ! The original value of the reconstruction in a cell [A]
  real :: u0min, u0max, u1min, u1max, u2min, u2max ! Minimum and maximum values of reconstructions [A]
  logical :: src_has_volume !< True if h0 has not been consumed
  logical :: tgt_has_volume !< True if h1 has not been consumed

  i0_last_thick_cell = 0
  do i0 = 1, n0
    u0_min(i0) = min(ppoly0_E(i0,1), ppoly0_E(i0,2))
    u0_max(i0) = max(ppoly0_E(i0,1), ppoly0_E(i0,2))
    if (h0(i0)>0.) i0_last_thick_cell = i0
  enddo

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
      stop 'remap_via_sub_cells: THIS SHOULD NEVER HAPPEN!'
    endif

  enddo

  ! Loop over each sub-cell to calculate average/integral values within each sub-cell.
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
    if (debug_bounds) then
      if (method<5 .and.(u_sub(i_sub)<u0_min(i0) .or. u_sub(i_sub)>u0_max(i0))) then
        write(0,*) 'Sub cell average is out of bounds',i_sub,'method=',method
        write(0,*) 'xa,xb: ',xa,xb
        write(0,*) 'Edge values: ',ppoly0_E(i0,:),'mean',u0(i0)
        write(0,*) 'a_c: ',(u0(i0)-ppoly0_E(i0,1))+(u0(i0)-ppoly0_E(i0,2))
        write(0,*) 'Polynomial coeffs: ',ppoly0_coefs(i0,:)
        write(0,*) 'Bounds min=',u0_min(i0),'max=',u0_max(i0)
        write(0,*) 'Average: ',u_sub(i_sub),'rel to min=',u_sub(i_sub)-u0_min(i0),'rel to max=',u_sub(i_sub)-u0_max(i0)
        call MOM_error( FATAL, 'MOM_remapping, remap_via_sub_cells: '//&
             'Sub-cell average is out of bounds!' )
      endif
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
    ! aka a genius algorithm for accurate conservation when remapping from Robert Hallberg (@Hallberg-NOAA).
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

  ! Loop over each target cell summing the integrals from sub-cells within the target cell.
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

  ! Check errors and bounds
  if (debug_bounds) then
    call measure_input_bounds( n0, h0, u0, ppoly0_E, h0tot, h0err, u0tot, u0err, u0min, u0max )
    call measure_output_bounds( n1, h1, u1, h1tot, h1err, u1tot, u1err, u1min, u1max )
    call measure_output_bounds( n0+n1+1, h_sub, u_sub, h2tot, h2err, u2tot, u2err, u2min, u2max )
    if (method<5) then ! We except PQM until we've debugged it
    if (     (abs(u1tot-u0tot)>(u0err+u1err)+uh_err+u02_err .and. abs(h1tot-h0tot)<h0err+h1err) &
        .or. (abs(u2tot-u0tot)>u0err+u2err+u02_err .and. abs(h2tot-h0tot)<h0err+h2err) &
        .or. (u1min<u0min .or. u1max>u0max) ) then
      write(0,*) 'method = ',method
      write(0,*) 'Source to sub-cells:'
      write(0,*) 'H: h0tot=',h0tot,'h2tot=',h2tot,'dh=',h2tot-h0tot,'h0err=',h0err,'h2err=',h2err
      if (abs(h2tot-h0tot)>h0err+h2err) &
        write(0,*) 'H non-conservation difference=',h2tot-h0tot,'allowed err=',h0err+h2err,' <-----!'
      write(0,*) 'UH: u0tot=',u0tot,'u2tot=',u2tot,'duh=',u2tot-u0tot,'u0err=',u0err,'u2err=',u2err,&
                 'adjustment err=',u02_err
      if (abs(u2tot-u0tot)>u0err+u2err) &
        write(0,*) 'U non-conservation difference=',u2tot-u0tot,'allowed err=',u0err+u2err,' <-----!'
      write(0,*) 'Sub-cells to target:'
      write(0,*) 'H: h2tot=',h2tot,'h1tot=',h1tot,'dh=',h1tot-h2tot,'h2err=',h2err,'h1err=',h1err
      if (abs(h1tot-h2tot)>h2err+h1err) &
        write(0,*) 'H non-conservation difference=',h1tot-h2tot,'allowed err=',h2err+h1err,' <-----!'
      write(0,*) 'UH: u2tot=',u2tot,'u1tot=',u1tot,'duh=',u1tot-u2tot,'u2err=',u2err,'u1err=',u1err,'uh_err=',uh_err
      if (abs(u1tot-u2tot)>u2err+u1err) &
        write(0,*) 'U non-conservation difference=',u1tot-u2tot,'allowed err=',u2err+u1err,' <-----!'
      write(0,*) 'Source to target:'
      write(0,*) 'H: h0tot=',h0tot,'h1tot=',h1tot,'dh=',h1tot-h0tot,'h0err=',h0err,'h1err=',h1err
      if (abs(h1tot-h0tot)>h0err+h1err) &
        write(0,*) 'H non-conservation difference=',h1tot-h0tot,'allowed err=',h0err+h1err,' <-----!'
      write(0,*) 'UH: u0tot=',u0tot,'u1tot=',u1tot,'duh=',u1tot-u0tot,'u0err=',u0err,'u1err=',u1err,'uh_err=',uh_err
      if (abs(u1tot-u0tot)>(u0err+u1err)+uh_err) &
        write(0,*) 'U non-conservation difference=',u1tot-u0tot,'allowed err=',u0err+u1err+uh_err,' <-----!'
      write(0,*) 'U: u0min=',u0min,'u1min=',u1min,'u2min=',u2min
      if (u1min<u0min) write(0,*) 'U minimum overshoot=',u1min-u0min,' <-----!'
      if (u2min<u0min) write(0,*) 'U2 minimum overshoot=',u2min-u0min,' <-----!'
      write(0,*) 'U: u0max=',u0max,'u1max=',u1max,'u2max=',u2max
      if (u1max>u0max) write(0,*) 'U maximum overshoot=',u1max-u0max,' <-----!'
      if (u2max>u0max) write(0,*) 'U2 maximum overshoot=',u2max-u0max,' <-----!'
      write(0,'(a3,6a24,2a3)') 'k','h0','left edge','u0','right edge','h1','u1','is','ie'
      do k = 1, max(n0,n1)
        if (k<=min(n0,n1)) then
          write(0,'(i3,1p6e24.16,2i3)') k,h0(k),ppoly0_E(k,1),u0(k),ppoly0_E(k,2),h1(k),u1(k),itgt_start(k),itgt_end(k)
        elseif (k>n0) then
          write(0,'(i3,96x,1p2e24.16,2i3)') k,h1(k),u1(k),itgt_start(k),itgt_end(k)
        else
          write(0,'(i3,1p4e24.16)') k,h0(k),ppoly0_E(k,1),u0(k),ppoly0_E(k,2)
        endif
      enddo
      write(0,'(a3,2a24)') 'k','u0','Polynomial coefficients'
      do k = 1, n0
        write(0,'(i3,1p6e24.16)') k,u0(k),ppoly0_coefs(k,:)
      enddo
      write(0,'(a3,3a24,a3,2a24)') 'k','Sub-cell h','Sub-cell u','Sub-cell hu','i0','xa','xb'
      xa = 0.
      dh0_eff = 0.
      do k = 1, n0+n1+1
        dh = h_sub(k)
        i0 = isub_src(k)
        dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
        xb = dh0_eff / h0_eff(i0) ! This expression yields xa <= xb <= 1.0
        xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
        write(0,'(i3,1p3e24.16,i3,1p2e24.16)') k,h_sub(k),u_sub(k),uh_sub(k),i0,xa,xb
        if (k<=n0+n1) then
          if (isub_src(k+1) /= i0) then
            dh0_eff = 0.; xa = 0.
          else
            xa = xb
          endif
        endif
      enddo
      call MOM_error( FATAL, 'MOM_remapping, remap_via_sub_cells: '//&
             'Remapping result is inconsistent!' )
    endif
    endif ! method<5
  endif ! debug_bounds

  ! Include the error remapping from source to sub-cells in the estimate of total remapping error
  uh_err = uh_err + u02_err

  if (present(ah_sub)) ah_sub(1:n0+n1+1) = h_sub(1:n0+n1+1)
  if (present(aisub_src)) aisub_src(1:n0+n1+1) = isub_src(1:n0+n1+1)
  if (present(aiss)) aiss(1:n0) = isrc_start(1:n0)
  if (present(aise)) aise(1:n0) = isrc_end(1:n0)

end subroutine remap_via_sub_cells

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
  real,          intent(in)    :: u0(:)  !< Cell means [A]
  real,          intent(in)    :: ppoly0_E(:,:)     !< Edge value of polynomial [A]
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
                check_reconstruction, check_remapping, force_bounds_in_subcell, answers_2018, answer_date)
  ! Arguments
  type(remapping_CS), intent(inout) :: CS !< Remapping control structure
  character(len=*),   intent(in)    :: remapping_scheme !< Remapping scheme to use
  logical, optional,  intent(in)    :: boundary_extrapolation !< Indicate to extrapolate in boundary cells
  logical, optional,  intent(in)    :: check_reconstruction !< Indicate to check reconstructions
  logical, optional,  intent(in)    :: check_remapping !< Indicate to check results of remapping
  logical, optional,  intent(in)    :: force_bounds_in_subcell !< Force subcells values to be bounded
  logical, optional,  intent(in)    :: answers_2018 !< If true use older, less accurate expressions.
  integer, optional,  intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Note that remapping_scheme is mandatory for initialize_remapping()
  call remapping_set_param(CS, remapping_scheme=remapping_scheme, boundary_extrapolation=boundary_extrapolation,  &
               check_reconstruction=check_reconstruction, check_remapping=check_remapping, &
               force_bounds_in_subcell=force_bounds_in_subcell, answers_2018=answers_2018, answer_date=answer_date)

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

!> Runs unit tests on remapping functions.
!! Should only be called from a single/root thread
!! Returns True if a test fails, otherwise False
logical function remapping_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables
  integer, parameter :: n0 = 4, n1 = 3, n2 = 6
  real :: h0(n0), x0(n0+1), u0(n0)  ! Thicknesses [H], interface heights [H] and values [A] for profile 0
  real :: h1(n1), x1(n1+1), u1(n1)  ! Thicknesses [H], interface heights [H] and values [A] for profile 1
  real :: dx1(n1+1)                 ! Interface height changes for profile 1 [H]
  real :: h2(n2), x2(n2+1), u2(n2)  ! Thicknesses [H], interface heights [H] and values [A] for profile 2
  data u0 /9., 3., -3., -9./   ! Linear profile, 4 at surface to -4 at bottom [A]
  data h0 /4*0.75/ ! 4 uniform layers with total depth of 3 [H]
  data h1 /3*1./   ! 3 uniform layers with total depth of 3 [H]
  data h2 /6*0.5/  ! 6 uniform layers with total depth of 3 [H]
  type(remapping_CS) :: CS !< Remapping control structure
  real, allocatable, dimension(:,:) :: ppoly0_E     ! Edge values of polynomials [A]
  real, allocatable, dimension(:,:) :: ppoly0_S     ! Edge slopes of polynomials [A H-1]
  real, allocatable, dimension(:,:) :: ppoly0_coefs ! Coefficients of polynomials [A]
  integer :: answer_date  ! The vintage of the expressions to test
  integer :: i
  real, parameter :: hNeglect_dflt = 1.0e-30 ! A thickness [H ~> m or kg m-2] that can be
                                      ! added to thicknesses in a denominator without
                                      ! changing the numerical result, except where
                                      ! a division by zero would otherwise occur.
  real :: err                         ! Errors in the remapped thicknesses [H] or values [A]
  real :: h_neglect, h_neglect_edge   ! Tiny thicknesses used in remapping [H]
  logical :: thisTest, v, fail

  v = verbose
  answer_date = 20190101 ! 20181231
  h_neglect = hNeglect_dflt
  h_neglect_edge = hNeglect_dflt ; if (answer_date < 20190101) h_neglect_edge = 1.0e-10

  write(stdout,*) '==== MOM_remapping: remapping_unit_tests ================='
  remapping_unit_tests = .false. ! Normally return false

  thisTest = .false.
  call buildGridFromH(n0, h0, x0)
  do i=1,n0+1
    err=x0(i)-0.75*real(i-1)
    if (abs(err)>real(i-1)*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(stdout,*) 'remapping_unit_tests: Failed buildGridFromH() 1'
  remapping_unit_tests = remapping_unit_tests .or. thisTest
  call buildGridFromH(n1, h1, x1)
  do i=1,n1+1
    err=x1(i)-real(i-1)
    if (abs(err)>real(i-1)*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(stdout,*) 'remapping_unit_tests: Failed buildGridFromH() 2'
  remapping_unit_tests = remapping_unit_tests .or. thisTest

  thisTest = .false.
  call initialize_remapping(CS, 'PPM_H4', answer_date=answer_date)
  if (verbose) write(stdout,*) 'h0 (test data)'
  if (verbose) call dumpGrid(n0,h0,x0,u0)

  call dzFromH1H2( n0, h0, n1, h1, dx1 )
  call remapping_core_w( CS, n0, h0, u0, n1, dx1, u1, h_neglect, h_neglect_edge)
  do i=1,n1
    err=u1(i)-8.*(0.5*real(1+n1)-real(i))
    if (abs(err)>real(n1-1)*epsilon(err)) thisTest = .true.
  enddo
  if (verbose) write(stdout,*) 'h1 (by projection)'
  if (verbose) call dumpGrid(n1,h1,x1,u1)
  if (thisTest) write(stdout,*) 'remapping_unit_tests: Failed remapping_core_w()'
  remapping_unit_tests = remapping_unit_tests .or. thisTest

  thisTest = .false.
  allocate(ppoly0_E(n0,2))
  allocate(ppoly0_S(n0,2))
  allocate(ppoly0_coefs(n0,CS%degree+1))

  ppoly0_E(:,:) = 0.0
  ppoly0_S(:,:) = 0.0
  ppoly0_coefs(:,:) = 0.0

  call edge_values_explicit_h4( n0, h0, u0, ppoly0_E, h_neglect=1e-10, answer_date=answer_date )
  call PPM_reconstruction( n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=answer_date )
  call PPM_boundary_extrapolation( n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect )

  thisTest = .false.
  call buildGridFromH(n2, h2, x2)

  if (verbose) write(stdout,*) 'Via sub-cells'
  thisTest = .false.
  call remap_via_sub_cells( n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                            n2, h2, INTEGRATION_PPM, .false., u2, err )
  if (verbose) call dumpGrid(n2,h2,x2,u2)

  do i=1,n2
    err=u2(i)-8./2.*(0.5*real(1+n2)-real(i))
    if (abs(err)>2.*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(stdout,*) 'remapping_unit_tests: Failed remap_via_sub_cells() 2'
  remapping_unit_tests = remapping_unit_tests .or. thisTest

  call remap_via_sub_cells( n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                            6, (/.125,.125,.125,.125,.125,.125/), INTEGRATION_PPM, .false., u2, err )
  if (verbose) call dumpGrid(6,h2,x2,u2)

  call remap_via_sub_cells( n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                            3, (/2.25,1.5,1./), INTEGRATION_PPM, .false., u2, err )
  if (verbose) call dumpGrid(3,h2,x2,u2)

  if (.not. remapping_unit_tests) write(stdout,*) 'Pass'

  write(stdout,*) '===== MOM_remapping: new remapping_unit_tests =================='

  deallocate(ppoly0_E, ppoly0_S, ppoly0_coefs)
  allocate(ppoly0_coefs(5,6))
  allocate(ppoly0_E(5,2))
  allocate(ppoly0_S(5,2))

  call PCM_reconstruction(3, (/1.,2.,4./), ppoly0_E(1:3,:), &
                          ppoly0_coefs(1:3,:) )
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_E(:,1), (/1.,2.,4./), 'PCM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_E(:,2), (/1.,2.,4./), 'PCM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_coefs(:,1), (/1.,2.,4./), 'PCM: P0')

  call PLM_reconstruction(3, (/1.,1.,1./), (/1.,3.,5./), ppoly0_E(1:3,:), &
                          ppoly0_coefs(1:3,:), h_neglect )
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_E(:,1), (/1.,2.,5./), 'Unlim PLM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_E(:,2), (/1.,4.,5./), 'Unlim PLM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_coefs(:,1), (/1.,2.,5./), 'Unlim PLM: P0')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_coefs(:,2), (/0.,2.,0./), 'Unlim PLM: P1')

  call PLM_reconstruction(3, (/1.,1.,1./), (/1.,2.,7./), ppoly0_E(1:3,:), &
                          ppoly0_coefs(1:3,:), h_neglect )
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_E(:,1), (/1.,1.,7./), 'Left lim PLM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_E(:,2), (/1.,3.,7./), 'Left lim PLM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_coefs(:,1), (/1.,1.,7./), 'Left lim PLM: P0')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_coefs(:,2), (/0.,2.,0./), 'Left lim PLM: P1')

  call PLM_reconstruction(3, (/1.,1.,1./), (/1.,6.,7./), ppoly0_E(1:3,:), &
                          ppoly0_coefs(1:3,:), h_neglect )
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_E(:,1), (/1.,5.,7./), 'Right lim PLM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_E(:,2), (/1.,7.,7./), 'Right lim PLM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_coefs(:,1), (/1.,5.,7./), 'Right lim PLM: P0')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_coefs(:,2), (/0.,2.,0./), 'Right lim PLM: P1')

  call PLM_reconstruction(3, (/1.,2.,3./), (/1.,4.,9./), ppoly0_E(1:3,:), &
                          ppoly0_coefs(1:3,:), h_neglect )
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_E(:,1), (/1.,2.,9./), 'Non-uniform line PLM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_E(:,2), (/1.,6.,9./), 'Non-uniform line PLM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_coefs(:,1), (/1.,2.,9./), 'Non-uniform line PLM: P0')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 3, ppoly0_coefs(:,2), (/0.,4.,0./), 'Non-uniform line PLM: P1')

  call edge_values_explicit_h4( 5, (/1.,1.,1.,1.,1./), (/1.,3.,5.,7.,9./), ppoly0_E, &
                                h_neglect=1e-10, answer_date=answer_date )
  ! The next two tests currently fail due to roundoff, but pass when given a reasonable tolerance.
  thisTest = test_answer(v, 5, ppoly0_E(:,1), (/0.,2.,4.,6.,8./), 'Line H4: left edges', tol=8.0e-15)
  remapping_unit_tests = remapping_unit_tests .or. thisTest
  thisTest = test_answer(v, 5, ppoly0_E(:,2), (/2.,4.,6.,8.,10./), 'Line H4: right edges', tol=1.0e-14)
  remapping_unit_tests = remapping_unit_tests .or. thisTest
  ppoly0_E(:,1) = (/0.,2.,4.,6.,8./)
  ppoly0_E(:,2) = (/2.,4.,6.,8.,10./)
  call PPM_reconstruction(5, (/1.,1.,1.,1.,1./), (/1.,3.,5.,7.,9./), ppoly0_E(1:5,:), &
                              ppoly0_coefs(1:5,:), h_neglect, answer_date=answer_date )
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_coefs(:,1), (/1.,2.,4.,6.,9./), 'Line PPM: P0')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_coefs(:,2), (/0.,2.,2.,2.,0./), 'Line PPM: P1')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_coefs(:,3), (/0.,0.,0.,0.,0./), 'Line PPM: P2')

  call edge_values_explicit_h4( 5, (/1.,1.,1.,1.,1./), (/1.,1.,7.,19.,37./), ppoly0_E, &
                                h_neglect=1e-10, answer_date=answer_date )
  ! The next two tests are now passing when answer_date >= 20190101, but otherwise only work to roundoff.
  thisTest = test_answer(v, 5, ppoly0_E(:,1), (/3.,0.,3.,12.,27./), 'Parabola H4: left edges', tol=2.7e-14)
  remapping_unit_tests = remapping_unit_tests .or. thisTest
  thisTest = test_answer(v, 5, ppoly0_E(:,2), (/0.,3.,12.,27.,48./), 'Parabola H4: right edges', tol=4.8e-14)
  remapping_unit_tests = remapping_unit_tests .or. thisTest
  ppoly0_E(:,1) = (/0.,0.,3.,12.,27./)
  ppoly0_E(:,2) = (/0.,3.,12.,27.,48./)
  call PPM_reconstruction(5, (/1.,1.,1.,1.,1./), (/0.,1.,7.,19.,37./), ppoly0_E(1:5,:), &
                          ppoly0_coefs(1:5,:), h_neglect, answer_date=answer_date )
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_E(:,1), (/0.,0.,3.,12.,37./), 'Parabola PPM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_E(:,2), (/0.,3.,12.,27.,37./), 'Parabola PPM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_coefs(:,1), (/0.,0.,3.,12.,37./), 'Parabola PPM: P0')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_coefs(:,2), (/0.,0.,6.,12.,0./), 'Parabola PPM: P1')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_coefs(:,3), (/0.,3.,3.,3.,0./), 'Parabola PPM: P2')

  ppoly0_E(:,1) = (/0.,0.,6.,10.,15./)
  ppoly0_E(:,2) = (/0.,6.,12.,17.,15./)
  call PPM_reconstruction(5, (/1.,1.,1.,1.,1./), (/0.,5.,7.,16.,15./), ppoly0_E(1:5,:), &
                          ppoly0_coefs(1:5,:), h_neglect, answer_date=answer_date )
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_E(:,1), (/0.,3.,6.,16.,15./), 'Limits PPM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_E(:,2), (/0.,6.,9.,16.,15./), 'Limits PPM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_coefs(:,1), (/0.,3.,6.,16.,15./), 'Limits PPM: P0')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_coefs(:,2), (/0.,6.,0.,0.,0./), 'Limits PPM: P1')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 5, ppoly0_coefs(:,3), (/0.,-3.,3.,0.,0./), 'Limits PPM: P2')

  call PLM_reconstruction(4, (/0.,1.,1.,0./), (/5.,4.,2.,1./), ppoly0_E(1:4,:), &
                          ppoly0_coefs(1:4,:), h_neglect )
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 4, ppoly0_E(1:4,1), (/5.,5.,3.,1./), 'PPM: left edges h=0110')
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 4, ppoly0_E(1:4,2), (/5.,3.,1.,1./), 'PPM: right edges h=0110')
  call remap_via_sub_cells( 4, (/0.,1.,1.,0./), (/5.,4.,2.,1./), ppoly0_E(1:4,:), &
                            ppoly0_coefs(1:4,:), &
                            2, (/1.,1./), INTEGRATION_PLM, .false., u2, err )
  remapping_unit_tests = remapping_unit_tests .or. &
    test_answer(v, 2, u2, (/4.,2./), 'PLM: remapped  h=0110->h=11')

  deallocate(ppoly0_E, ppoly0_S, ppoly0_coefs)

  ! This line carries out tests on some older remapping schemes.
  remapping_unit_tests = remapping_unit_tests .or. remapping_attic_unit_tests(verbose)

  if (.not. remapping_unit_tests) write(stdout,*) 'Pass'

  write(stdout,*) '=== MOM_remapping: interpolation and reintegration unit tests ==='
  if (verbose) write(stdout,*) '- - - - - - - - - - interpolation tests  - - - - - - - - -'

  fail = test_interp(verbose, 'Identity: 3 layer', &
                     3, (/1.,2.,3./), (/1.,2.,3.,4./), &
                     3, (/1.,2.,3./), (/1.,2.,3.,4./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_interp(verbose, 'A: 3 layer to 2', &
                     3, (/1.,1.,1./), (/1.,2.,3.,4./), &
                     2, (/1.5,1.5/), (/1.,2.5,4./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_interp(verbose, 'B: 2 layer to 3', &
                     2, (/1.5,1.5/), (/1.,4.,7./), &
                     3, (/1.,1.,1./), (/1.,3.,5.,7./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_interp(verbose, 'C: 3 layer (vanished middle) to 2', &
                     3, (/1.,0.,2./), (/1.,2.,2.,3./), &
                     2, (/1.,2./), (/1.,2.,3./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_interp(verbose, 'D: 3 layer (deep) to 3', &
                     3, (/1.,2.,3./), (/1.,2.,4.,7./), &
                     2, (/2.,2./), (/1.,3.,5./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_interp(verbose, 'E: 3 layer to 3 (deep)', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     3, (/2.,3.,4./), (/1.,3.,6.,8./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_interp(verbose, 'F: 3 layer to 4 with vanished top/botton', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     4, (/0.,2.,5.,0./), (/0.,1.,3.,8.,0./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_interp(verbose, 'Fs: 3 layer to 4 with vanished top/botton (shallow)', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     4, (/0.,2.,4.,0./), (/0.,1.,3.,7.,0./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_interp(verbose, 'Fd: 3 layer to 4 with vanished top/botton (deep)', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     4, (/0.,2.,6.,0./), (/0.,1.,3.,8.,0./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  if (verbose) write(stdout,*) '- - - - - - - - - - reintegration tests  - - - - - - - - -'

  fail = test_reintegrate(verbose, 'Identity: 3 layer', &
                     3, (/1.,2.,3./), (/-5.,2.,1./), &
                     3, (/1.,2.,3./), (/-5.,2.,1./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_reintegrate(verbose, 'A: 3 layer to 2', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     2, (/3.,3./), (/-4.,2./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_reintegrate(verbose, 'A: 3 layer to 2 (deep)', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     2, (/3.,4./), (/-4.,2./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_reintegrate(verbose, 'A: 3 layer to 2 (shallow)', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     2, (/3.,2./), (/-4.,1.5/) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_reintegrate(verbose, 'B: 3 layer to 4 with vanished top/bottom', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     4, (/0.,3.,3.,0./), (/0.,-4.,2.,0./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_reintegrate(verbose, 'C: 3 layer to 4 with vanished top//middle/bottom', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     5, (/0.,3.,0.,3.,0./), (/0.,-4.,0.,2.,0./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_reintegrate(verbose, 'D: 3 layer to 3 (vanished)', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     3, (/0.,0.,0./), (/0.,0.,0./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_reintegrate(verbose, 'D: 3 layer (vanished) to 3', &
                     3, (/0.,0.,0./), (/-5.,2.,1./), &
                     3, (/2.,2.,2./), (/0., 0., 0./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_reintegrate(verbose, 'D: 3 layer (vanished) to 3 (vanished)', &
                     3, (/0.,0.,0./), (/-5.,2.,1./), &
                     3, (/0.,0.,0./), (/0., 0., 0./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  fail = test_reintegrate(verbose, 'D: 3 layer (vanished) to 3 (vanished)', &
                     3, (/0.,0.,0./), (/0.,0.,0./), &
                     3, (/0.,0.,0./), (/0., 0., 0./) )
  remapping_unit_tests = remapping_unit_tests .or. fail

  if (.not. remapping_unit_tests) write(stdout,*) 'Pass'

end function remapping_unit_tests

!> Returns true if any cell of u and u_true are not identical. Returns false otherwise.
logical function test_answer(verbose, n, u, u_true, label, tol)
  logical,            intent(in) :: verbose !< If true, write results to stdout
  integer,            intent(in) :: n      !< Number of cells in u
  real, dimension(n), intent(in) :: u      !< Values to test [A]
  real, dimension(n), intent(in) :: u_true !< Values to test against (correct answer) [A]
  character(len=*),   intent(in) :: label  !< Message
  real,     optional, intent(in) :: tol    !< The tolerance for differences between u and u_true [A]
  ! Local variables
  real :: tolerance ! The tolerance for differences between u and u_true [A]
  integer :: k

  tolerance = 0.0 ; if (present(tol)) tolerance = tol
  test_answer = .false.
  do k = 1, n
    if (abs(u(k) - u_true(k)) > tolerance) test_answer = .true.
  enddo
  if (test_answer .or. verbose) then
    write(stdout,'(a4,2a24,1x,a)') 'k','Calculated value','Correct value',label
    do k = 1, n
      if (abs(u(k) - u_true(k)) > tolerance) then
        write(stdout,'(i4,1p2e24.16,a,1pe24.16,a)') k,u(k),u_true(k),' err=',u(k)-u_true(k),' < wrong'
        write(stderr,'(i4,1p2e24.16,a,1pe24.16,a)') k,u(k),u_true(k),' err=',u(k)-u_true(k),' < wrong'
      else
        write(stdout,'(i4,1p2e24.16)') k,u(k),u_true(k)
      endif
    enddo
  endif

end function test_answer

!> Returns true if a test of interpolate_column() produces the wrong answer
logical function test_interp(verbose, msg, nsrc, h_src, u_src, ndest, h_dest, u_true)
  logical,                  intent(in) :: verbose !< If true, write results to stdout
  character(len=*),         intent(in) :: msg   !< Message to label test
  integer,                  intent(in) :: nsrc  !< Number of source cells
  real, dimension(nsrc),    intent(in) :: h_src !< Thickness of source cells [H]
  real, dimension(nsrc+1),  intent(in) :: u_src !< Values at source cell interfaces [A]
  integer,                  intent(in) :: ndest !< Number of destination cells
  real, dimension(ndest),   intent(in) :: h_dest !< Thickness of destination cells [H]
  real, dimension(ndest+1), intent(in) :: u_true !< Correct value at destination cell interfaces [A]
  ! Local variables
  real, dimension(ndest+1) :: u_dest ! Interpolated value at destination cell interfaces [A]
  integer :: k
  real :: error

  ! Interpolate from src to dest
  call interpolate_column(nsrc, h_src, u_src, ndest, h_dest, u_dest, .true.)

  test_interp = .false.
  do k=1,ndest+1
    if (u_dest(k)/=u_true(k)) test_interp = .true.
  enddo
  if (verbose .or. test_interp) then
    write(stdout,'(2a)') ' Test: ',msg
    write(stdout,'(a3,3(a24))') 'k','u_result','u_true','error'
    do k=1,ndest+1
      error = u_dest(k)-u_true(k)
      if (error==0.) then
        write(stdout,'(i3,3(1pe24.16))') k,u_dest(k),u_true(k),u_dest(k)-u_true(k)
      else
        write(stdout,'(i3,3(1pe24.16),1x,a)') k,u_dest(k),u_true(k),u_dest(k)-u_true(k),'<--- WRONG!'
        write(stderr,'(i3,3(1pe24.16),1x,a)') k,u_dest(k),u_true(k),u_dest(k)-u_true(k),'<--- WRONG!'
      endif
    enddo
  endif
end function test_interp

!> Returns true if a test of reintegrate_column() produces the wrong answer
logical function test_reintegrate(verbose, msg, nsrc, h_src, uh_src, ndest, h_dest, uh_true)
  logical,                intent(in) :: verbose !< If true, write results to stdout
  character(len=*),       intent(in) :: msg   !< Message to label test
  integer,                intent(in) :: nsrc  !< Number of source cells
  real, dimension(nsrc),  intent(in) :: h_src !< Thickness of source cells [H]
  real, dimension(nsrc),  intent(in) :: uh_src !< Values of source cell stuff [A H]
  integer,                intent(in) :: ndest  !< Number of destination cells
  real, dimension(ndest), intent(in) :: h_dest !< Thickness of destination cells [H]
  real, dimension(ndest), intent(in) :: uh_true !< Correct value of destination cell stuff [A H]
  ! Local variables
  real, dimension(ndest) :: uh_dest ! Reintegrated value on destination cells [A H]
  integer :: k
  real :: error

  ! Interpolate from src to dest
  call reintegrate_column(nsrc, h_src, uh_src, ndest, h_dest, uh_dest)

  test_reintegrate = .false.
  do k=1,ndest
    if (uh_dest(k)/=uh_true(k)) test_reintegrate = .true.
  enddo
  if (verbose .or. test_reintegrate) then
    write(stdout,'(2a)') ' Test: ',msg
    write(stdout,'(a3,3(a24))') 'k','uh_result','uh_true','error'
    do k=1,ndest
      error = uh_dest(k)-uh_true(k)
      if (error==0.) then
        write(stdout,'(i3,3(1pe24.16))') k,uh_dest(k),uh_true(k),uh_dest(k)-uh_true(k)
      else
        write(stdout,'(i3,3(1pe24.16),1x,a)') k,uh_dest(k),uh_true(k),uh_dest(k)-uh_true(k),'<--- WRONG!'
        write(stderr,'(i3,3(1pe24.16),1x,a)') k,uh_dest(k),uh_true(k),uh_dest(k)-uh_true(k),'<--- WRONG!'
      endif
    enddo
  endif
end function test_reintegrate

!> Convenience function for printing grid to screen
subroutine dumpGrid(n,h,x,u)
  integer, intent(in) :: n !< Number of cells
  real, dimension(:), intent(in) :: h !< Cell thickness [H]
  real, dimension(:), intent(in) :: x !< Interface delta [H]
  real, dimension(:), intent(in) :: u !< Cell average values [A]
  integer :: i
  write(stdout,'("i=",20i10)') (i,i=1,n+1)
  write(stdout,'("x=",20es10.2)') (x(i),i=1,n+1)
  write(stdout,'("i=",5x,20i10)') (i,i=1,n)
  write(stdout,'("h=",5x,20es10.2)') (h(i),i=1,n)
  write(stdout,'("u=",5x,20es10.2)') (u(i),i=1,n)
end subroutine dumpGrid

end module MOM_remapping
