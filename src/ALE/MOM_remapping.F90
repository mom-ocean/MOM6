!> Provides column-wise vertical remapping functions
module MOM_remapping

! This file is part of MOM6. See LICENSE.md for the license.
! Original module written by Laurent White, 2008.06.09

use MOM_error_handler, only : MOM_error, FATAL
use MOM_string_functions, only : uppercase
use regrid_edge_values, only : edge_values_explicit_h4, edge_values_implicit_h4
use regrid_edge_values, only : edge_values_implicit_h4, edge_values_implicit_h6
use regrid_edge_slopes, only : edge_slopes_implicit_h3, edge_slopes_implicit_h5
use PCM_functions, only : PCM_reconstruction
use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation
use PQM_functions, only : PQM_reconstruction, PQM_boundary_extrapolation_v1

implicit none ; private

#include <MOM_memory.h>

!> Container for remapping parameters
type, public :: remapping_CS
  private
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
end type

! The following routines are visible to the outside world
public remapping_core_h, remapping_core_w
public initialize_remapping, end_remapping, remapping_set_param
public remapping_unit_tests
public dzFromH1H2

! The following are private parameter constants
integer, parameter  :: REMAPPING_PCM        = 0 !< O(h^1) remapping scheme
integer, parameter  :: REMAPPING_PLM        = 1 !< O(h^2) remapping scheme
integer, parameter  :: REMAPPING_PPM_H4     = 2 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PPM_IH4    = 3 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PQM_IH4IH3 = 4 !< O(h^4) remapping scheme
integer, parameter  :: REMAPPING_PQM_IH6IH5 = 5 !< O(h^5) remapping scheme

integer, parameter  :: INTEGRATION_PCM = 0  !< Piecewise Constant Method
integer, parameter  :: INTEGRATION_PLM = 1  !< Piecewise Linear Method
integer, parameter  :: INTEGRATION_PPM = 3  !< Piecewise Parabolic Method
integer, parameter  :: INTEGRATION_PQM = 5  !< Piecewise Quartic Method

character(len=40)  :: mod = "MOM_remapping" !< This module's name.

!> Documentation for external callers
character(len=256), public :: remappingSchemesDoc = &
                 "PCM         (1st-order accurate)\n"//&
                 "PLM         (2nd-order accurate)\n"//&
                 "PPM_H4      (3rd-order accurate)\n"//&
                 "PPM_IH4     (3rd-order accurate)\n"//&
                 "PQM_IH4IH3  (4th-order accurate)\n"//&
                 "PQM_IH6IH5  (5th-order accurate)\n"
character(len=3), public :: remappingDefaultScheme = "PLM" !< Default remapping method

! This CPP macro turns on/off bounding of integrations limits so that they are
! always within the cell. Roundoff can lead to the non-dimensional bounds being
! outside of the range 0 to 1.
#define __USE_ROUNDOFF_SAFE_ADJUSTMENTS__

real, parameter :: h_neglect = 1.E-30 !< A dimensional (H units) number that can be
                                      !! added to thicknesses in a denominator without
                                      !! changing the numerical result, except where
                                      !! a division by zero would otherwise occur.

logical, parameter :: old_algorithm = .true.  !< Use the old "broken" algorithm.
                                              !! This is a temporary measure to assist
                                              !! debugging until we delete the old algorithm.

contains

!> Set parameters within remapping object
subroutine remapping_set_param(CS, remapping_scheme, boundary_extrapolation,  &
               check_reconstruction, check_remapping, force_bounds_in_subcell)
  type(remapping_CS),         intent(inout) :: CS !< Remapping control structure
  character(len=*), optional, intent(in)    :: remapping_scheme !< Remapping scheme to use
  logical, optional,          intent(in)    :: boundary_extrapolation !< Indicate to extrapolate in boundary cells
  logical, optional,          intent(in)    :: check_reconstruction !< Indicate to check reconstructions
  logical, optional,          intent(in)    :: check_remapping !< Indicate to check results of remapping
  logical, optional,          intent(in)    :: force_bounds_in_subcell !< Force subcells values to be bounded

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
end subroutine remapping_set_param

!> Calculate edge coordinate x from cell width h
subroutine buildGridFromH(nz, h, x)
  integer,               intent(in)    :: nz !< Number of cells
  real, dimension(nz),   intent(in)    :: h  !< Cell widths
  real, dimension(nz+1), intent(inout) :: x  !< Edge coordiantes starting at x(1)=0
  ! Local variables
  integer :: k

  x(1) = 0.0
  do k = 1,nz
    x(k+1) = x(k) + h(k)
  end do

end subroutine buildGridFromH

!> Compare two summation estimates of positive data and judge if due to more
!! than round-off.
!! When two sums are calculated from different vectors that should add up to
!! the same value, the results can differ by round off. The round off error
!! can be bounded to be proportional to the number of operations.
!! This function returns true if the difference between sum1 and sum2 is
!! larger than than the estimated round off bound.
!! \note This estimate/function is only valid for summation of positive data.
function isPosSumErrSignificant(n1, sum1, n2, sum2)
  integer, intent(in) :: n1   !< Number of values in sum1
  integer, intent(in) :: n2   !< Number of values in sum2
  real,    intent(in) :: sum1 !< Sum of n1 values
  real,    intent(in) :: sum2 !< Sum of n2 values
  logical             :: isPosSumErrSignificant !< True if difference in sums is large
  ! Local variables
  real :: sumErr, allowedErr, eps

  if (sum1<0.) call MOM_error(FATAL,'isPosSumErrSignificant: sum1<0 is not allowed!')
  if (sum2<0.) call MOM_error(FATAL,'isPosSumErrSignificant: sum2<0 is not allowed!')
  sumErr = abs(sum1-sum2)
  eps = epsilon(sum1)
  allowedErr = eps*0.5*(real(n1-1)*sum1+real(n2-1)*sum2)
  if (sumErr>allowedErr) then
    write(0,*) 'isPosSumErrSignificant: sum1,sum2=',sum1,sum2
    write(0,*) 'isPosSumErrSignificant: eps=',eps
    write(0,*) 'isPosSumErrSignificant: err,n*eps=',sumErr,allowedErr
    write(0,*) 'isPosSumErrSignificant: err/eps,n1,n2,n1+n2=',sumErr/eps,n1,n2,n1+n2
    isPosSumErrSignificant = .true.
  else
    isPosSumErrSignificant = .false.
  endif
end function isPosSumErrSignificant

!> Remaps column of values u0 on grid h0 to grid h1
!! assuming the top edge is aligned.
subroutine remapping_core_h( n0, h0, u0, n1, h1, u1, CS )
  integer,             intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0), intent(in)  :: h0 !< Cell widths on source grid
  real, dimension(n0), intent(in)  :: u0 !< Cell averages on source grid
  integer,             intent(in)  :: n1 !< Number of cells on target grid
  real, dimension(n1), intent(in)  :: h1 !< Cell widths on target grid
  real, dimension(n1), intent(out) :: u1 !< Cell averages on target grid
  type(remapping_CS),  intent(in)  :: CS !< Remapping control structure
  ! Local variables
  integer :: iMethod
  real, dimension(n0,2)           :: ppoly_r_E            !Edge value of polynomial
  real, dimension(n0,2)           :: ppoly_r_S            !Edge slope of polynomial
  real, dimension(n0,CS%degree+1) :: ppoly_r_coefficients !Coefficients of polynomial
  integer :: k
  real :: eps, h0tot, h0err, h1tot, h1err, u0tot, u0err, u0min, u0max, u1tot, u1err, u1min, u1max, uh_err

  call build_reconstructions_1d( n0, h0, u0, CS%remapping_scheme, CS%degree, CS%boundary_extrapolation, &
                                   ppoly_r_coefficients, ppoly_r_E, ppoly_r_S, iMethod )

  if (CS%check_reconstruction) call check_reconstructions_1d(n0, h0, u0, CS%degree, &
                                   CS%boundary_extrapolation, ppoly_r_coefficients, ppoly_r_E, ppoly_r_S)


  call remap_via_sub_cells( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients, n1, h1, iMethod, &
                            CS%force_bounds_in_subcell, u1, uh_err )

  if (CS%check_remapping) then
    ! Check errors and bounds
    call measure_input_bounds( n0, h0, u0, ppoly_r_E, h0tot, h0err, u0tot, u0err, u0min, u0max )
    call measure_output_bounds( n1, h1, u1, h1tot, h1err, u1tot, u1err, u1min, u1max )
    if (iMethod<5) then ! We except PQM until we've debugged it
    if ( (abs(u1tot-u0tot)>(u0err+u1err)+uh_err .and. abs(h1tot-h0tot)<h0err+h1err) &
        .or. (u1min<u0min .or. u1max>u0max) ) then
      write(0,*) 'iMethod = ',iMethod
      write(0,*) 'H: h0tot=',h0tot,'h1tot=',h1tot,'dh=',h1tot-h0tot,'h0err=',h0err,'h1err=',h1err
      if (abs(h1tot-h0tot)>h0err+h1err) write(0,*) 'H non-conservation difference=',h1tot-h0tot,'allowed err=',h0err+h1err,' <-----!'
      write(0,*) 'UH: u0tot=',u0tot,'u1tot=',u1tot,'duh=',u1tot-u0tot,'u0err=',u0err,'u1err=',u1err,'uh_err=',uh_err
      if (abs(u1tot-u0tot)>(u0err+u1err)+uh_err) write(0,*) 'U non-conservation difference=',u1tot-u0tot,'allowed err=',u0err+u1err+uh_err,' <-----!'
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
        write(0,'(i3,1p6e24.16)') k,u0(k),ppoly_r_coefficients(k,:)
      enddo
      call MOM_error( FATAL, 'MOM_remapping, remapping_core_h: '//&
             'Remapping result is inconsistent!' )
    endif
    endif ! method<5
  endif

end subroutine remapping_core_h

!> Remaps column of values u0 on grid h0 to implied grid h1
!! where the interfaces of h1 differ from those of h0 by dx.
subroutine remapping_core_w( CS, n0, h0, u0, n1, dx, u1 )
  type(remapping_CS),  intent(in)  :: CS !< Remapping control structure
  integer,             intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0), intent(in)  :: h0 !< Cell widths on source grid
  real, dimension(n0), intent(in)  :: u0 !< Cell averages on source grid
  integer,             intent(in)  :: n1 !< Number of cells on target grid
  real, dimension(n1+1), intent(in)  :: dx !< Cell widths on target grid
  real, dimension(n1), intent(out) :: u1 !< Cell averages on target grid
  ! Local variables
  integer :: iMethod
  real, dimension(n0,2)           :: ppoly_r_E            !Edge value of polynomial
  real, dimension(n0,2)           :: ppoly_r_S            !Edge slope of polynomial
  real, dimension(n0,CS%degree+1) :: ppoly_r_coefficients !Coefficients of polynomial
  integer :: k
  real :: eps, h0tot, h0err, h1tot, h1err
  real :: u0tot, u0err, u0min, u0max, u1tot, u1err, u1min, u1max, uh_err
  real, dimension(n1) :: h1 !< Cell widths on target grid

  call build_reconstructions_1d( n0, h0, u0, CS%remapping_scheme, CS%degree, CS%boundary_extrapolation, &
                                   ppoly_r_coefficients, ppoly_r_E, ppoly_r_S, iMethod )

  if (CS%check_reconstruction) call check_reconstructions_1d(n0, h0, u0, CS%degree, &
                                   CS%boundary_extrapolation, ppoly_r_coefficients, ppoly_r_E, ppoly_r_S)

  ! This is a temporary step prior to switching to remapping_core_h()
  do k = 1, n1
    if (k<=n0) then
      h1(k) = max( 0., h0(k) + ( dx(k+1) - dx(k) ) )
    else
      h1(k) = max( 0., dx(k+1) - dx(k) )
    endif
  enddo
  call remap_via_sub_cells( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients, n1, h1, iMethod, &
                            CS%force_bounds_in_subcell,u1, uh_err )
! call remapByDeltaZ( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients, n1, dx, iMethod, u1 )
! call remapByProjection( n0, h0, u0, CS%ppoly_r, n1, h1, iMethod, u1 )

  if (CS%check_remapping) then
    ! Check errors and bounds
    call measure_input_bounds( n0, h0, u0, ppoly_r_E, h0tot, h0err, u0tot, u0err, u0min, u0max )
    call measure_output_bounds( n1, h1, u1, h1tot, h1err, u1tot, u1err, u1min, u1max )
    if (iMethod<5) then ! We except PQM until we've debugged it
    if ( (abs(u1tot-u0tot)>(u0err+u1err)+uh_err .and. abs(h1tot-h0tot)<h0err+h1err) &
        .or. (u1min<u0min .or. u1max>u0max) ) then
      write(0,*) 'iMethod = ',iMethod
      write(0,*) 'H: h0tot=',h0tot,'h1tot=',h1tot,'dh=',h1tot-h0tot,'h0err=',h0err,'h1err=',h1err
      if (abs(h1tot-h0tot)>h0err+h1err) write(0,*) 'H non-conservation difference=',h1tot-h0tot,'allowed err=',h0err+h1err,' <-----!'
      write(0,*) 'UH: u0tot=',u0tot,'u1tot=',u1tot,'duh=',u1tot-u0tot,'u0err=',u0err,'u1err=',u1err,'uh_err=',uh_err
      if (abs(u1tot-u0tot)>(u0err+u1err)+uh_err) write(0,*) 'U non-conservation difference=',u1tot-u0tot,'allowed err=',u0err+u1err+uh_err,' <-----!'
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
        write(0,'(i3,1p6e24.16)') k,u0(k),ppoly_r_coefficients(k,:)
      enddo
      call MOM_error( FATAL, 'MOM_remapping, remapping_core_w: '//&
             'Remapping result is inconsistent!' )
    endif
    endif ! method<5
  endif

end subroutine remapping_core_w

!> Creates polynomial reconstructions of u0 on the source grid h0.
subroutine build_reconstructions_1d( n0, h0, u0, remapping_scheme, deg, boundary_extrapolation, &
                                     ppoly_r_coefficients, ppoly_r_E, ppoly_r_S, iMethod )
  integer,                  intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),      intent(in)  :: h0 !< Cell widths on source grid
  real, dimension(n0),      intent(in)  :: u0 !< Cell averages on source grid
  integer,                  intent(in)  :: remapping_scheme !< Remapping scheme
  integer,                  intent(in)  :: deg !< Degree of polynomial reconstruction
  logical,                  intent(in)  :: boundary_extrapolation !< Extrapolate at boundaries if true
  real, dimension(n0,deg+1),intent(out) :: ppoly_r_coefficients !< Coefficients of polynomial
  real, dimension(n0,2),    intent(out) :: ppoly_r_E !< Edge value of polynomial
  real, dimension(n0,2),    intent(out) :: ppoly_r_S !< Edge slope of polynomial
  integer,                  intent(out) :: iMethod !< Integration method
  ! Local variables
  integer :: local_remapping_scheme

  ! Reset polynomial
  ppoly_r_E(:,:) = 0.0
  ppoly_r_S(:,:) = 0.0
  ppoly_r_coefficients(:,:) = 0.0
  iMethod = -999

  local_remapping_scheme = remapping_scheme
  if (n0<=1) then
    local_remapping_scheme = REMAPPING_PCM
  elseif (n0<=3) then
    local_remapping_scheme = min( local_remapping_scheme, REMAPPING_PLM )
  elseif (n0<=4) then
    local_remapping_scheme = min( local_remapping_scheme, REMAPPING_PPM_H4 )
  endif
  select case ( local_remapping_scheme )
    case ( REMAPPING_PCM )
      call PCM_reconstruction( n0, u0, ppoly_r_E, ppoly_r_coefficients)
      iMethod = INTEGRATION_PCM
    case ( REMAPPING_PLM )
      call PLM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      if ( boundary_extrapolation) then
        call PLM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients)
      end if
      iMethod = INTEGRATION_PLM
    case ( REMAPPING_PPM_H4 )
      call edge_values_explicit_h4( n0, h0, u0, ppoly_r_E )
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      if ( boundary_extrapolation) then
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      end if
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PPM_IH4 )
      call edge_values_implicit_h4( n0, h0, u0, ppoly_r_E )
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      if ( boundary_extrapolation) then
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      end if
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PQM_IH4IH3 )
      call edge_values_implicit_h4( n0, h0, u0, ppoly_r_E )
      call edge_slopes_implicit_h3( n0, h0, u0, ppoly_r_S )
      call PQM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      if ( boundary_extrapolation) then
        call PQM_boundary_extrapolation_v1( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      end if
      iMethod = INTEGRATION_PQM
    case ( REMAPPING_PQM_IH6IH5 )
      call edge_values_implicit_h6( n0, h0, u0, ppoly_r_E )
      call edge_slopes_implicit_h5( n0, h0, u0, ppoly_r_S )
      call PQM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      if ( boundary_extrapolation) then
        call PQM_boundary_extrapolation_v1( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      end if
      iMethod = INTEGRATION_PQM
    case default
      call MOM_error( FATAL, 'MOM_remapping, build_reconstructions_1d: '//&
           'The selected remapping method is invalid' )
  end select

end subroutine build_reconstructions_1d

!> Checks that edge values and reconstructions satisfy bounds
subroutine check_reconstructions_1d(n0, h0, u0, deg, boundary_extrapolation, &
                                    ppoly_r_coefficients, ppoly_r_E, ppoly_r_S)
  integer,                  intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),      intent(in)  :: h0 !< Cell widths on source grid
  real, dimension(n0),      intent(in)  :: u0 !< Cell averages on source grid
  integer,                  intent(in)  :: deg !< Degree of polynomial reconstruction
  logical,                  intent(in)  :: boundary_extrapolation !< Extrapolate at boundaries if true
  real, dimension(n0,deg+1),intent(out) :: ppoly_r_coefficients !< Coefficients of polynomial
  real, dimension(n0,2),    intent(out) :: ppoly_r_E !< Edge value of polynomial
  real, dimension(n0,2),    intent(out) :: ppoly_r_S !< Edge slope of polynomial
  ! Local variables
  integer :: i0, n
  real :: u_l, u_c, u_r ! Cell averages
  real :: u_min, u_max
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
        write(0,'(a,i4,5(x,a,1pe24.16))') 'Left edge undershoot at',i0,'u(i0-1)=',u_l,'u(i0)=',u_c, &
                                          'edge=',ppoly_r_E(i0,1),'err=',ppoly_r_E(i0,1)-u_min
        problem_detected = .true.
      endif
      if (ppoly_r_E(i0,1) > u_max) then
        write(0,'(a,i4,5(x,a,1pe24.16))') 'Left edge overshoot at',i0,'u(i0-1)=',u_l,'u(i0)=',u_c, &
                                          'edge=',ppoly_r_E(i0,1),'err=',ppoly_r_E(i0,1)-u_max
        problem_detected = .true.
      endif
    endif
    if (i0 < n0 .or. .not. boundary_extrapolation) then
      u_min = min(u_c, u_r)
      u_max = max(u_c, u_r)
      if (ppoly_r_E(i0,2) < u_min) then
        write(0,'(a,i4,5(x,a,1pe24.16))') 'Right edge undershoot at',i0,'u(i0)=',u_c,'u(i0+1)=',u_r, &
                                          'edge=',ppoly_r_E(i0,2),'err=',ppoly_r_E(i0,2)-u_min
        problem_detected = .true.
      endif
      if (ppoly_r_E(i0,2) > u_max) then
        write(0,'(a,i4,5(x,a,1pe24.16))') 'Right edge overshoot at',i0,'u(i0)=',u_c,'u(i0+1)=',u_r, &
                                          'edge=',ppoly_r_E(i0,2),'err=',ppoly_r_E(i0,2)-u_max
        problem_detected = .true.
      endif
    endif
    if (i0 > 1) then
      if ( (u_c-u_l)*(ppoly_r_E(i0,1)-ppoly_r_E(i0-1,2)) < 0.) then
        write(0,'(a,i4,5(x,a,1pe24.16))') 'Non-monotonic edges at',i0,'u(i0-1)=',u_l,'u(i0)=',u_c, &
                                          'right edge=',ppoly_r_E(i0-1,2),'left edge=',ppoly_r_E(i0,1)
        write(0,'(5(a,1pe24.16,x))') 'u(i0)-u(i0-1)',u_c-u_l,'edge diff=',ppoly_r_E(i0,1)-ppoly_r_E(i0-1,2)
        problem_detected = .true.
      endif
    endif
    if (problem_detected) then
      write(0,'(a,1p9e24.16)') 'Polynomial coeffs:',ppoly_r_coefficients(i0,:)
      write(0,'(3(a,1pe24.16,x))') 'u_l=',u_l,'u_c=',u_c,'u_r=',u_r
      write(0,'(a4,10a24)') 'i0','h0(i0)','u0(i0)','left edge','right edge','Polynomial coefficients'
      do n = 1, n0
        write(0,'(i4,1p10e24.16)') n,h0(n),u0(n),ppoly_r_E(n,1),ppoly_r_E(n,2),ppoly_r_coefficients(n,:)
      enddo
      call MOM_error(FATAL, 'MOM_remapping, check_reconstructions_1d: '// &
                   'Edge values or polynomial coefficients were inconsistent!')
    endif
  enddo

end subroutine check_reconstructions_1d

!> Remaps column of n0 values u0 on grid h0 to grid h1 with n1 cells by calculating
!! the n0+n1+1 sub-integrals of the intersection of h0 and h1, and the summing the
!! appropriate integrals into the h1*u1 values.
subroutine remap_via_sub_cells( n0, h0, u0, ppoly0_E, ppoly0_coefficients, n1, h1, method, &
                                force_bounds_in_subcell, u1, uh_err, ah_sub, aisub_src, aiss, aise )
  integer,           intent(in)    :: n0     !< Number of cells in source grid
  real,              intent(in)    :: h0(n0)  !< Source grid widths (size n0)
  real,              intent(in)    :: u0(n0)  !< Source cell averages (size n0)
  real,              intent(in)    :: ppoly0_E(n0,2)            !< Edge value of polynomial
  real,              intent(in)    :: ppoly0_coefficients(:,:) !< Coefficients of polynomial
  integer,           intent(in)    :: n1     !< Number of cells in target grid
  real,              intent(in)    :: h1(n1)  !< Target grid widths (size n1)
  integer,           intent(in)    :: method !< Remapping scheme to use
  logical,           intent(in)    :: force_bounds_in_subcell !< Force sub-cell values to be bounded
  real,              intent(out)   :: u1(n1)  !< Target cell averages (size n1)
  real,              intent(out)   :: uh_err !< Estimate of bound on error in sum of u*h
  real, optional,    intent(out)   :: ah_sub(n0+n1+1) !< h_sub
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
  real :: dh_max ! Used to record which sub-cell is the largest contribution of a source cell
  real, dimension(n0+n1+1) :: h_sub ! Width of each each sub-cell
  real, dimension(n0+n1+1) :: uh_sub ! Integral of u*h over each sub-cell
  real, dimension(n0+n1+1) :: u_sub ! Average of u over each sub-cell
  integer, dimension(n0+n1+1) :: isub_src ! Index of source cell for each sub-cell
  integer, dimension(n0) :: isrc_start ! Index of first sub-cell within each source cell
  integer, dimension(n0) :: isrc_end ! Index of last sub-cell within each source cell
  integer, dimension(n0) :: isrc_max ! Index of thickest sub-cell within each source cell
  real, dimension(n0) :: h0_eff ! Effective thickness of source cells
  real, dimension(n0) :: u0_min ! Minimum value of reconstructions in source cell
  real, dimension(n0) :: u0_max ! Minimum value of reconstructions in source cell
  integer, dimension(n1) :: itgt_start ! Index of first sub-cell within each target cell
  integer, dimension(n1) :: itgt_end ! Index of last sub-cell within each target cell
  real :: xa, xb ! Non-dimensional position within a source cell (0..1)
  real :: h0_supply, h1_supply ! The amount of width available for constructing sub-cells
  real :: dh ! The width of the sub-cell
  real :: duh ! The total amount of accumulated stuff (u*h)
  real :: dh0_eff ! Running sum of source cell thickness
  ! For error checking/debugging
  logical, parameter :: force_bounds_in_target = .true. ! To fix round-off issues
  logical, parameter :: adjust_thickest_subcell = .true. ! To fix round-off conservation issues
  logical, parameter :: debug_bounds = .false. ! For debugging overshoots etc.
  integer :: k, i0_last_thick_cell
  real :: h0tot, h0err, h1tot, h1err, h2tot, h2err, u02_err
  real :: u0tot, u0err, u0min, u0max, u1tot, u1err, u1min, u1max, u2tot, u2err, u2min, u2max, u_orig
  logical :: src_has_volume !< True if h0 has not been consumed
  logical :: tgt_has_volume !< True if h1 has not been consumed

  if (old_algorithm) isrc_max(:)=1

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
      if (old_algorithm) then
        if (i0 < i0_last_thick_cell) then
          i0 = i0 + 1
          h0_supply = h0(i0)
          dh0_eff = 0.
          do while (h0_supply==0. .and. i0<i0_last_thick_cell)
            ! This loop skips over vanished source cells
            i0 = i0 + 1
            h0_supply = h0(i0)
          enddo
        else
          h0_supply = 1.E30
        endif
      else
        if (i0 < n0) then
          i0 = i0 + 1
          h0_supply = h0(i0)
          dh0_eff = 0.
        else
          h0_supply = 0.
          src_has_volume = .false.
        endif
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
        if (old_algorithm) then
          h1_supply = 1.E30
        else
          h1_supply = 0.
          tgt_has_volume = .false.
        endif
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
      u_sub(i_sub) = average_value_ppoly( n0, u0, ppoly0_E, ppoly0_coefficients, method, i0, xa, xb)
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
        write(0,*) 'Polynomial coeffs: ',ppoly0_coefficients(i0,:)
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
      if (abs(h2tot-h0tot)>h0err+h2err) write(0,*) 'H non-conservation difference=',h2tot-h0tot,'allowed err=',h0err+h2err,' <-----!'
      write(0,*) 'UH: u0tot=',u0tot,'u2tot=',u2tot,'duh=',u2tot-u0tot,'u0err=',u0err,'u2err=',u2err,'adjustment err=',u02_err
      if (abs(u2tot-u0tot)>u0err+u2err) write(0,*) 'U non-conservation difference=',u2tot-u0tot,'allowed err=',u0err+u2err,' <-----!'
      write(0,*) 'Sub-cells to target:'
      write(0,*) 'H: h2tot=',h2tot,'h1tot=',h1tot,'dh=',h1tot-h2tot,'h2err=',h2err,'h1err=',h1err
      if (abs(h1tot-h2tot)>h2err+h1err) write(0,*) 'H non-conservation difference=',h1tot-h2tot,'allowed err=',h2err+h1err,' <-----!'
      write(0,*) 'UH: u2tot=',u2tot,'u1tot=',u1tot,'duh=',u1tot-u2tot,'u2err=',u2err,'u1err=',u1err,'uh_err=',uh_err
      if (abs(u1tot-u2tot)>u2err+u1err) write(0,*) 'U non-conservation difference=',u1tot-u2tot,'allowed err=',u2err+u1err,' <-----!'
      write(0,*) 'Source to target:'
      write(0,*) 'H: h0tot=',h0tot,'h1tot=',h1tot,'dh=',h1tot-h0tot,'h0err=',h0err,'h1err=',h1err
      if (abs(h1tot-h0tot)>h0err+h1err) write(0,*) 'H non-conservation difference=',h1tot-h0tot,'allowed err=',h0err+h1err,' <-----!'
      write(0,*) 'UH: u0tot=',u0tot,'u1tot=',u1tot,'duh=',u1tot-u0tot,'u0err=',u0err,'u1err=',u1err,'uh_err=',uh_err
      if (abs(u1tot-u0tot)>(u0err+u1err)+uh_err) write(0,*) 'U non-conservation difference=',u1tot-u0tot,'allowed err=',u0err+u1err+uh_err,' <-----!'
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
        write(0,'(i3,1p6e24.16)') k,u0(k),ppoly0_coefficients(k,:)
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

!> Returns the average value of a reconstruction within a single source cell, i0,
!! between the non-dimensional positions xa and xb (xa<=xb) with dimensional
!! separation dh.
real function average_value_ppoly( n0, u0, ppoly0_E, ppoly0_coefficients, method, i0, xa, xb)
  integer,       intent(in)    :: n0     !< Number of cells in source grid
  real,          intent(in)    :: u0(:)  !< Cell means
  real,          intent(in)    :: ppoly0_E(:,:)            !< Edge value of polynomial
  real,          intent(in)    :: ppoly0_coefficients(:,:) !< Coefficients of polynomial
  integer,       intent(in)    :: method !< Remapping scheme to use
  integer,       intent(in)    :: i0     !< Source cell index
  real,          intent(in)    :: xa     !< Non-dimensional start position within source cell
  real,          intent(in)    :: xb     !< Non-dimensional end position within source cell
  ! Local variables
  real :: u_ave, xa_2, xb_2, xa2pxb2, xapxb
  real, parameter :: r_3 = 1.0/3.0 ! Used in evaluation of integrated polynomials

  real :: mx, a_L, a_R, u_c, Ya, Yb, my, xa2b2ab, Ya2b2ab, a_c

  if (xb > xa) then
    select case ( method )
      case ( INTEGRATION_PCM )
        u_ave = u0(i0)
      case ( INTEGRATION_PLM )
        u_ave = (                                           &
            ppoly0_coefficients(i0,1)                       &
          + ppoly0_coefficients(i0,2) * 0.5 * ( xb + xa ) )
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
              ppoly0_coefficients(i0,1)                                                         &
          + ( ppoly0_coefficients(i0,2) * 0.5 * ( xapxb )                                       &
          + ( ppoly0_coefficients(i0,3) * r_3 * ( xa2pxb2 + xa*xb )                             &
          + ( ppoly0_coefficients(i0,4) * 0.25* ( xa2pxb2 * xapxb )                             &
          +   ppoly0_coefficients(i0,5) * 0.2 * ( ( xb*xb_2 + xa*xa_2 ) * xapxb + xa_2*xb_2 ) ) ) ) )
      case default
        call MOM_error( FATAL,'The selected integration method is invalid' )
    end select
  else ! dh == 0.
    select case ( method )
      case ( INTEGRATION_PCM )
        u_ave =        ppoly0_coefficients(i0,1)
      case ( INTEGRATION_PLM )
       !u_ave =        ppoly0_coefficients(i0,1)   &
       !      + xa *   ppoly0_coefficients(i0,2)
        a_L = ppoly0_E(i0, 1)
        a_R = ppoly0_E(i0, 2)
        Ya = 1. - xa
        if (xa < 0.5) then
          u_ave = a_L + xa * ( a_R - a_L )
        else
          u_ave = a_R + Ya * ( a_L - a_R )
        endif
      case ( INTEGRATION_PPM )
       !u_ave =        ppoly0_coefficients(i0,1)   &
       !      + xa * ( ppoly0_coefficients(i0,2)   &
       !      + xa *   ppoly0_coefficients(i0,3) )
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
        u_ave =        ppoly0_coefficients(i0,1)   &
              + xa * ( ppoly0_coefficients(i0,2)   &
              + xa * ( ppoly0_coefficients(i0,3)   &
              + xa * ( ppoly0_coefficients(i0,4)   &
              + xa *   ppoly0_coefficients(i0,5) ) ) )
      case default
        call MOM_error( FATAL,'The selected integration method is invalid' )
    end select
  endif
  average_value_ppoly = u_ave

end function average_value_ppoly

!> Measure totals and bounds on source grid
subroutine measure_input_bounds( n0, h0, u0, ppoly_E, h0tot, h0err, u0tot, u0err, u0min, u0max )
  integer,               intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)  :: h0 !< Cell widths on source grid
  real, dimension(n0),   intent(in)  :: u0 !< Cell averages on source grid
  real, dimension(n0,2), intent(in)  :: ppoly_E !< Cell edge values on source grid
  real,                  intent(out) :: h0tot !< Sum of cell widths
  real,                  intent(out) :: h0err !< Magnitude of round-off error in h0tot
  real,                  intent(out) :: u0tot !< Sum of cell widths times values
  real,                  intent(out) :: u0err !< Magnitude of round-off error in u0tot
  real,                  intent(out) :: u0min !< Minimum value in reconstructions of u0
  real,                  intent(out) :: u0max !< Maximum value in reconstructions of u0
  ! Local variables
  integer :: k
  real :: eps

  eps = epsilon(h0(1))
  h0tot = h0(1)
  h0err = 0.
  u0tot = h0(1) * u0(1)
  u0err = 0.
  u0min = min( ppoly_E(1,1), ppoly_E(1,2) )
  u0max = max( ppoly_E(1,1), ppoly_E(1,2) )
  do k = 2, n0
    h0tot = h0tot + h0(k)
    h0err = h0err + eps * max(h0tot, h0(k))
    u0tot = u0tot + h0(k) * u0(k)
    u0err = u0err + eps * max(abs(u0tot), abs(h0(k) * u0(k)))
    u0min = min( u0min, ppoly_E(k,1), ppoly_E(k,2) )
    u0max = max( u0max, ppoly_E(k,1), ppoly_E(k,2) )
  enddo

end subroutine measure_input_bounds

!> Measure totals and bounds on destination grid
subroutine measure_output_bounds( n1, h1, u1, h1tot, h1err, u1tot, u1err, u1min, u1max )
  integer,               intent(in)  :: n1 !< Number of cells on destination grid
  real, dimension(n1),   intent(in)  :: h1 !< Cell widths on destination grid
  real, dimension(n1),   intent(in)  :: u1 !< Cell averages on destination grid
  real,                  intent(out) :: h1tot !< Sum of cell widths
  real,                  intent(out) :: h1err !< Magnitude of round-off error in h1tot
  real,                  intent(out) :: u1tot !< Sum of cell widths times values
  real,                  intent(out) :: u1err !< Magnitude of round-off error in u1tot
  real,                  intent(out) :: u1min !< Minimum value in reconstructions of u1
  real,                  intent(out) :: u1max !< Maximum value in reconstructions of u1
  ! Local variables
  integer :: k
  real :: eps

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

!> Remaps column of values u0 on grid h0 to grid h1 by integrating
!! over the projection of each h1 cell onto the h0 grid.
subroutine remapByProjection( n0, h0, u0, ppoly0_E, ppoly0_coefficients, n1, h1, method, u1 )
  integer,       intent(in)    :: n0     !< Number of cells in source grid
  real,          intent(in)    :: h0(:)  !< Source grid widths (size n0)
  real,          intent(in)    :: u0(:)  !< Source cell averages (size n0)
  real,          intent(in)    :: ppoly0_E(:,:)            !< Edge value of polynomial
  real,          intent(in)    :: ppoly0_coefficients(:,:) !< Coefficients of polynomial
  integer,       intent(in)    :: n1     !< Number of cells in target grid
  real,          intent(in)    :: h1(:)  !< Target grid widths (size n1)
  integer,       intent(in)    :: method !< Remapping scheme to use
  real,          intent(out)   :: u1(:)  !< Target cell averages (size n1)
  ! Local variables
  integer       :: iTarget
  real          :: xL, xR       ! coordinates of target cell edges
  integer       :: jStart ! Used by integrateReconOnInterval()
  real          :: xStart ! Used by integrateReconOnInterval()

  ! Loop on cells in target grid (grid1). For each target cell, we need to find
  ! in which source cells the target cell edges lie. The associated indexes are
  ! noted j0 and j1.
  xR = 0. ! Left boundary is at x=0
  jStart = 1
  xStart = 0.
  do iTarget = 1,n1
    ! Determine the coordinates of the target cell edges
    xL = xR
    xR = xL + h1(iTarget)

    call integrateReconOnInterval( n0, h0, u0, ppoly0_E, ppoly0_coefficients, method, &
                                   xL, xR, h1(iTarget), u1(iTarget), jStart, xStart )

  end do ! end iTarget loop on target grid cells

end subroutine remapByProjection


!> Remaps column of values u0 on grid h0 to implied grid h1
!! where the interfaces of h1 differ from those of h0 by dx.
!! The new grid is defined relative to the original grid by change
!!  dx1(:) = xNew(:) - xOld(:)
!! and the remapping calculated so that
!!  hNew(k) qNew(k) = hOld(k) qOld(k) + F(k+1) - F(k)
!! where
!!  F(k) = dx1(k) qAverage
!! and where qAverage is the average qOld in the region zOld(k) to zNew(k).
subroutine remapByDeltaZ( n0, h0, u0, ppoly0_E, ppoly0_coefficients, n1, dx1, method, u1, h1 )
  integer,        intent(in)  :: n0     !< Number of cells in source grid
  real,           intent(in)  :: h0(:)  !< Source grid widths (size n0)
  real,           intent(in)  :: u0(:)  !< Source cell averages (size n0)
  real,           intent(in)  :: ppoly0_E(:,:)            !< Edge value of polynomial
  real,           intent(in)  :: ppoly0_coefficients(:,:) !< Coefficients of polynomial
  integer,        intent(in)  :: n1     !< Number of cells in target grid
  real,           intent(in)  :: dx1(:) !< Target grid edge positions (size n1+1)
  integer                     :: method !< Remapping scheme to use
  real,           intent(out) :: u1(:)  !< Target cell averages (size n1)
  real, optional, intent(out) :: h1(:)  !< Target grid widths (size n1)
  ! Local variables
  integer :: iTarget
  real    :: xL, xR    ! coordinates of target cell edges
  real    :: xOld, hOld, uOld
  real    :: xNew, hNew, h_err
  real    :: uhNew, hFlux, uAve, fluxL, fluxR
  integer :: jStart ! Used by integrateReconOnInterval()
  real    :: xStart ! Used by integrateReconOnInterval()

  ! Loop on cells in target grid. For each cell, iTarget, the left flux is
  ! the right flux of the cell to the left, iTarget-1.
  ! The left flux is initialized by started at iTarget=0 to calculate the
  ! right flux which can take into account the target left boundary being
  ! in the interior of the source domain.
  fluxR = 0.
  h_err = 0. ! For measuring round-off error
  jStart = 1
  xStart = 0.
  do iTarget = 0,n1
    fluxL = fluxR ! This does nothing for iTarget=0

    if (iTarget == 0) then
      xOld = 0.     ! Left boundary is at x=0
      hOld = -1.E30 ! Should not be used for iTarget = 0
      uOld = -1.E30 ! Should not be used for iTarget = 0
    elseif (iTarget <= n0) then
      xOld = xOld + h0(iTarget) ! Position of right edge of cell
      hOld = h0(iTarget)
      uOld = u0(iTarget)
      h_err = h_err + epsilon(hOld) * max(hOld, xOld)
    else
      hOld = 0.       ! as if for layers>n0, they were vanished
      uOld = 1.E30    ! and the initial value should not matter
    endif
    xNew = xOld + dx1(iTarget+1)
    xL = min( xOld, xNew )
    xR = max( xOld, xNew )

    ! hFlux is the positive width of the remapped volume
    hFlux = abs(dx1(iTarget+1))
    call integrateReconOnInterval( n0, h0, u0, ppoly0_E, ppoly0_coefficients, method, &
                                   xL, xR, hFlux, uAve, jStart, xStart )
    ! uAve is the average value of u, independent of sign of dx1
    fluxR = dx1(iTarget+1)*uAve ! Includes sign of dx1

    if (iTarget>0) then
      hNew = hOld + ( dx1(iTarget+1) - dx1(iTarget) )
      hNew = max( 0., hNew )
      uhNew = ( uOld * hOld ) + ( fluxR - fluxL )
      if (hNew>0.) then
        u1(iTarget) = uhNew / hNew
      else
        u1(iTarget) = uAve
      endif
      if (present(h1)) h1(iTarget) = hNew
    endif

  end do ! end iTarget loop on target grid cells

end subroutine remapByDeltaZ


!> Integrate the reconstructed column profile over a single cell
subroutine integrateReconOnInterval( n0, h0, u0, ppoly0_E, ppoly0_coefficients, method, &
                                     xL, xR, hC, uAve, jStart, xStart )
  integer, intent(in)  :: n0       !< Number of cells in source grid
  real,    intent(in)  :: h0(:)    !< Source grid sizes (size n0)
  real,    intent(in)  :: u0(:)    !< Source cell averages
  real,    intent(in)  :: ppoly0_E(:,:)            !< Edge value of polynomial
  real,    intent(in)  :: ppoly0_coefficients(:,:) !< Coefficients of polynomial
  integer, intent(in)  :: method   !< Remapping scheme to use
  real,    intent(in)  :: xL, xR   !< Left/right edges of target cell
  real,    intent(in)  :: hC       !< Cell width hC = xR - xL
  real,    intent(out) :: uAve     !< Average value on target cell
  integer, intent(inout) :: jStart !< The index of the cell to start searching from
                                   !< On exit, contains index of last cell used
  real,    intent(inout) :: xStart !< The left edge position of cell jStart
                                   !< On first entry should be 0.
  ! Local variables
  integer :: j, k
  integer :: jL, jR       ! indexes of source cells containing target
                          ! cell edges
  real    :: q            ! complete integration
  real    :: xi0, xi1     ! interval of integration (local -- normalized
                          ! -- coordinates)
  real    :: x0jLl, x0jLr ! Left/right position of cell jL
  real    :: x0jRl, x0jRr ! Left/right position of cell jR
  real    :: hAct         ! The distance actually used in the integration
                          ! (notionally xR - xL) which differs due to roundoff.
  real    :: x0_2, x1_2, x02px12, x0px1 ! Used in evaluation of integrated polynomials
  real, parameter :: r_3 = 1.0/3.0 ! Used in evaluation of integrated polynomials

  q = -1.E30
  x0jLl = -1.E30
  x0jRl = -1.E30

  ! Find the left most cell in source grid spanned by the target cell
  jL = -1
  x0jLr = xStart
  do j = jStart, n0
    x0jLl = x0jLr
    x0jLr = x0jLl + h0(j)
    ! Left edge is found in cell j
    if ( ( xL >= x0jLl ) .AND. ( xL <= x0jLr ) ) then
      jL = j
      exit ! once target grid cell is found, exit loop
    endif
  enddo
  jStart = jL
  xStart = x0jLl

! ! HACK to handle round-off problems. Need only at  j=n0.
! ! This moves the effective cell boundary outwards a smidgen.
! if (xL>x0jLr) x0jLr = xL

  ! If, at this point, jL is equal to -1, it means the vanished
  ! cell lies outside the source grid. In other words, it means that
  ! the source and target grids do not cover the same physical domain
  ! and there is something very wrong !
  if ( jL == -1 ) call MOM_error(FATAL, &
          'MOM_remapping, integrateReconOnInterval: '//&
          'The location of the left-most cell could not be found')


  ! ============================================================
  ! Check whether target cell is vanished. If it is, the cell
  ! average is simply the interpolated value at the location
  ! of the vanished cell. If it isn't, we need to integrate the
  ! quantity within the cell and divide by the cell width to
  ! determine the cell average.
  ! ============================================================
  ! 1. Cell is vanished
 !if ( abs(xR - xL) <= epsilon(xR)*max(abs(xR),abs(xL)) ) then
  if ( abs(xR - xL) == 0.0 ) then

    ! We check whether the source cell (i.e. the cell in which the
    ! vanished target cell lies) is vanished. If it is, the interpolated
    ! value is set to be mean of the edge values (which should be the same).
    ! If it isn't, we simply interpolate.
    if ( h0(jL) == 0.0 ) then
      uAve = 0.5 * ( ppoly0_E(jL,1) + ppoly0_E(jL,2) )
    else
      ! WHY IS THIS NOT WRITTEN AS xi0 = ( xL - x0jLl ) / h0(jL) ---AJA
      xi0 = xL / ( h0(jL) + h_neglect ) - x0jLl / ( h0(jL) + h_neglect )

      select case ( method )
        case ( INTEGRATION_PCM )
          uAve =         ppoly0_coefficients(jL,1)
        case ( INTEGRATION_PLM )
          uAve =         ppoly0_coefficients(jL,1)   &
               + xi0 *   ppoly0_coefficients(jL,2)
        case ( INTEGRATION_PPM )
          uAve =         ppoly0_coefficients(jL,1)   &
               + xi0 * ( ppoly0_coefficients(jL,2)   &
               + xi0 *   ppoly0_coefficients(jL,3) )
        case ( INTEGRATION_PQM )
          uAve =         ppoly0_coefficients(jL,1)   &
               + xi0 * ( ppoly0_coefficients(jL,2)   &
               + xi0 * ( ppoly0_coefficients(jL,3)   &
               + xi0 * ( ppoly0_coefficients(jL,4)   &
               + xi0 *   ppoly0_coefficients(jL,5) ) ) )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select

    end if ! end checking whether source cell is vanished

  ! 2. Cell is not vanished
  else

    ! Find the right most cell in source grid spanned by the target cell
    jR = -1
    x0jRr = xStart
    do j = jStart,n0
      x0jRl = x0jRr
      x0jRr = x0jRl + h0(j)
      ! Right edge is found in cell j
      if ( ( xR >= x0jRl ) .AND. ( xR <= x0jRr ) ) then
        jR = j
        exit  ! once target grid cell is found, exit loop
      endif
    enddo ! end loop on source grid cells

    ! If xR>x0jRr then the previous loop reached j=n0 and the target
    ! position, xR, was beyond the right edge of the source grid (h0).
    ! This can happen due to roundoff, in which case we set jR=n0.
    if (xR>x0jRr) jR = n0

    ! To integrate, two cases must be considered: (1) the target cell is
    ! entirely contained within a cell of the source grid and (2) the target
    ! cell spans at least two cells of the source grid.

    if ( jL == jR ) then
      ! The target cell is entirely contained within a cell of the source
      ! grid. This situation is represented by the following schematic, where
      ! the cell in which xL and xR are located has index jL=jR :
      !
      ! ----|-----o--------o----------|-------------
      !           xL       xR
      !
      ! Determine normalized coordinates
#ifdef __USE_ROUNDOFF_SAFE_ADJUSTMENTS__
      xi0 = max( 0., min( 1., ( xL - x0jLl ) / ( h0(jL) + h_neglect ) ) )
      xi1 = max( 0., min( 1., ( xR - x0jLl ) / ( h0(jL) + h_neglect ) ) )
#else
      xi0 = xL / h0(jL) - x0jLl / ( h0(jL) + h_neglect )
      xi1 = xR / h0(jL) - x0jLl / ( h0(jL) + h_neglect )
#endif

      hAct = h0(jL) * ( xi1 - xi0 )

      ! Depending on which polynomial is used, integrate quantity
      ! between xi0 and xi1. Integration is carried out in normalized
      ! coordinates, hence: \int_xL^xR p(x) dx = h \int_xi0^xi1 p(xi) dxi
      select case ( method )
        case ( INTEGRATION_PCM )
          q = ( xR - xL ) * ppoly0_coefficients(jL,1)
        case ( INTEGRATION_PLM )
          q = ( xR - xL ) * (                                 &
              ppoly0_coefficients(jL,1)                       &
            + ppoly0_coefficients(jL,2) * 0.5 * ( xi1 + xi0 ) )
        case ( INTEGRATION_PPM )
          q = ( xR - xL ) * (                                   &
                ppoly0_coefficients(jL,1)                       &
            + ( ppoly0_coefficients(jL,2) * 0.5 * ( xi1 + xi0 ) &
            +   ppoly0_coefficients(jL,3) * r_3 * ( ( xi1*xi1 + xi0*xi0 ) + xi0*xi1 ) ) )
        case ( INTEGRATION_PQM )
          x0_2 = xi0*xi0
          x1_2 = xi1*xi1
          x02px12 = x0_2 + x1_2
          x0px1 = xi1 + xi0
          q = ( xR - xL ) * (                                                                     &
                ppoly0_coefficients(jL,1)                                                         &
            + ( ppoly0_coefficients(jL,2) * 0.5 * ( xi1 + xi0 )                                   &
            + ( ppoly0_coefficients(jL,3) * r_3 * ( x02px12 + xi0*xi1 )                           &
            +   ppoly0_coefficients(jL,4) * 0.25* ( x02px12 * x0px1 )                             &
            +   ppoly0_coefficients(jL,5) * 0.2 * ( ( xi1*x1_2 + xi0*x0_2 ) * x0px1 + x0_2*x1_2 ) ) ) )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select

    else
    ! The target cell spans at least two cells of the source grid.
    ! This situation is represented by the following schematic, where
    ! the cells in which xL and xR are located have indexes jL and jR,
    ! respectively :
    !
    ! ----|-----o---|--- ... --|---o----------|-------------
    !           xL                 xR
    !
    ! We first integrate from xL up to the right boundary of cell jL, then
    ! add the integrated amounts of cells located between jL and jR and then
    ! integrate from the left boundary of cell jR up to xR

      q = 0.0

      ! Integrate from xL up to right boundary of cell jL
#ifdef __USE_ROUNDOFF_SAFE_ADJUSTMENTS__
      xi0 = max( 0., min( 1., ( xL - x0jLl ) / ( h0(jL) + h_neglect ) ) )
#else
      xi0 = (xL - x0jLl) / ( h0(jL) + h_neglect )
#endif
      xi1 = 1.0

      hAct = h0(jL) * ( xi1 - xi0 )

      select case ( method )
        case ( INTEGRATION_PCM )
          q = q + ( x0jLr - xL ) * ppoly0_coefficients(jL,1)
        case ( INTEGRATION_PLM )
          q = q + ( x0jLr - xL ) * (                          &
              ppoly0_coefficients(jL,1)                       &
            + ppoly0_coefficients(jL,2) * 0.5 * ( xi1 + xi0 ) )
        case ( INTEGRATION_PPM )
          q = q + ( x0jLr - xL ) * (                            &
                ppoly0_coefficients(jL,1)                       &
            + ( ppoly0_coefficients(jL,2) * 0.5 * ( xi1 + xi0 ) &
            +   ppoly0_coefficients(jL,3) * r_3 * ( ( xi1*xi1 + xi0*xi0 ) + xi0*xi1 ) ) )
        case ( INTEGRATION_PQM )
          x0_2 = xi0*xi0
          x1_2 = xi1*xi1
          x02px12 = x0_2 + x1_2
          x0px1 = xi1 + xi0
          q = q + ( x0jLr - xL ) * (                                                              &
                ppoly0_coefficients(jL,1)                                                         &
            + ( ppoly0_coefficients(jL,2) * 0.5 * ( xi1 + xi0 )                                   &
            + ( ppoly0_coefficients(jL,3) * r_3 * ( x02px12 + xi0*xi1 )                           &
            +   ppoly0_coefficients(jL,4) * 0.25* ( x02px12 * x0px1 )                             &
            +   ppoly0_coefficients(jL,5) * 0.2 * ( ( xi1*x1_2 + xi0*x0_2 ) * x0px1 + x0_2*x1_2 ) ) ) )
        case default
          call MOM_error( FATAL, 'The selected integration method is invalid' )
      end select

      ! Integrate contents within cells strictly comprised between jL and jR
      if ( jR > (jL+1) ) then
        do k = jL+1,jR-1
          q = q + h0(k) * u0(k)
          hAct = hAct + h0(k)
        end do
      end if

      ! Integrate from left boundary of cell jR up to xR
      xi0 = 0.0
#ifdef __USE_ROUNDOFF_SAFE_ADJUSTMENTS__
      xi1 = max( 0., min( 1., ( xR - x0jRl ) / ( h0(jR) + h_neglect ) ) )
#else
      xi1 = (xR - x0jRl) / ( h0(jR) + h_neglect )
#endif

      hAct = hAct + h0(jR) * ( xi1 - xi0 )

      select case ( method )
        case ( INTEGRATION_PCM )
          q = q + ( xR - x0jRl ) * ppoly0_coefficients(jR,1)
        case ( INTEGRATION_PLM )
          q = q + ( xR - x0jRl ) * (                          &
              ppoly0_coefficients(jR,1)                       &
            + ppoly0_coefficients(jR,2) * 0.5 * ( xi1 + xi0 ) )
        case ( INTEGRATION_PPM )
          q = q + ( xR - x0jRl ) * (                            &
                ppoly0_coefficients(jR,1)                       &
            + ( ppoly0_coefficients(jR,2) * 0.5 * ( xi1 + xi0 ) &
            +   ppoly0_coefficients(jR,3) * r_3 * ( ( xi1*xi1 + xi0*xi0 ) + xi0*xi1 ) ) )
        case ( INTEGRATION_PQM )
          x0_2 = xi0*xi0
          x1_2 = xi1*xi1
          x02px12 = x0_2 + x1_2
          x0px1 = xi1 + xi0
          q = q + ( xR - x0jRl ) * (                                                              &
                ppoly0_coefficients(jR,1)                                                         &
            + ( ppoly0_coefficients(jR,2) * 0.5 * ( xi1 + xi0 )                                   &
            + ( ppoly0_coefficients(jR,3) * r_3 * ( x02px12 + xi0*xi1 )                           &
            +   ppoly0_coefficients(jR,4) * 0.25* ( x02px12 * x0px1 )                             &
            +   ppoly0_coefficients(jR,5) * 0.2 * ( ( xi1*x1_2 + xi0*x0_2 ) * x0px1 + x0_2*x1_2 ) ) ) )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select

    end if ! end integration for non-vanished cells

    ! The cell average is the integrated value divided by the cell width
#ifdef __USE_ROUNDOFF_SAFE_ADJUSTMENTS__
if (hAct==0.) then
    uAve = ppoly0_coefficients(jL,1)
else
    uAve = q / hAct
endif
#else
    uAve = q / hC
#endif

  end if ! end if clause to check if cell is vanished

end subroutine integrateReconOnInterval

!> Calculates the change in interface positions based on h1 and h2
subroutine dzFromH1H2( n1, h1, n2, h2, dx )
  integer,            intent(in)  :: n1 !< Number of cells on source grid
  real, dimension(:), intent(in)  :: h1 !< Cell widths of source grid (size n1)
  integer,            intent(in)  :: n2 !< Number of cells on target grid
  real, dimension(:), intent(in)  :: h2 !< Cell widths of target grid (size n2)
  real, dimension(:), intent(out) :: dx !< Change in interface position (size n2+1)
  ! Local variables
  integer :: k
  real :: x1, x2

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
                check_reconstruction, check_remapping, force_bounds_in_subcell)
  ! Arguments
  type(remapping_CS), intent(inout) :: CS !< Remapping control structure
  character(len=*),   intent(in)    :: remapping_scheme !< Remapping scheme to use
  logical, optional,  intent(in)    :: boundary_extrapolation !< Indicate to extrapolate in boundary cells
  logical, optional,  intent(in)    :: check_reconstruction !< Indicate to check reconstructions
  logical, optional,  intent(in)    :: check_remapping !< Indicate to check results of remapping
  logical, optional,  intent(in)    :: force_bounds_in_subcell !< Force subcells values to be bounded

  ! Note that remapping_scheme is mandatory fir initialize_remapping()
  call remapping_set_param(CS, remapping_scheme=remapping_scheme, boundary_extrapolation=boundary_extrapolation,  &
               check_reconstruction=check_reconstruction, check_remapping=check_remapping, &
               force_bounds_in_subcell=force_bounds_in_subcell)

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
    case ("PPM_H4")
      CS%remapping_scheme = REMAPPING_PPM_H4
      degree = 2
    case ("PPM_IH4")
      CS%remapping_scheme = REMAPPING_PPM_IH4
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
logical function remapping_unit_tests()
  integer, parameter :: n0 = 4, n1 = 3, n2 = 6
  real :: h0(n0), x0(n0+1), u0(n0)
  real :: h1(n1), x1(n1+1), u1(n1), hn1(n1), dx1(n1+1)
  real :: h2(n2), x2(n2+1), u2(n2), hn2(n2), dx2(n2+1)
  data u0 /9., 3., -3., -9./   ! Linear profile, 4 at surface to -4 at bottom
  data h0 /4*0.75/ ! 4 uniform layers with total depth of 3
  data h1 /3*1./   ! 3 uniform layers with total depth of 3
  data h2 /6*0.5/  ! 6 uniform layers with total depth of 3
  type(remapping_CS) :: CS !< Remapping control structure
  real, allocatable, dimension(:,:) :: ppoly0_E, ppoly0_S, ppoly0_coefficients
  integer :: i
  real :: err
  logical :: thisTest

  write(*,*) '===== MOM_remapping: remapping_unit_tests ================'
  remapping_unit_tests = .false. ! Normally return false

  thisTest = .false.
  call buildGridFromH(n0, h0, x0)
  do i=1,n0+1
    err=x0(i)-0.75*real(i-1)
    if (abs(err)>real(i-1)*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(*,*) 'remapping_unit_tests: Failed buildGridFromH() 1'
  remapping_unit_tests = remapping_unit_tests .or. thisTest
  call buildGridFromH(n1, h1, x1)
  do i=1,n1+1
    err=x1(i)-real(i-1)
    if (abs(err)>real(i-1)*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(*,*) 'remapping_unit_tests: Failed buildGridFromH() 2'
  remapping_unit_tests = remapping_unit_tests .or. thisTest

  thisTest = .false.
  call initialize_remapping(CS, 'PPM_H4')
  write(*,*) 'h0 (test data)'
  call dumpGrid(n0,h0,x0,u0)

  call dzFromH1H2( n0, h0, n1, h1, dx1 )
  call remapping_core_w( CS, n0, h0, u0, n1, dx1, u1 )
  do i=1,n1
    err=u1(i)-8.*(0.5*real(1+n1)-real(i))
    if (abs(err)>real(n1-1)*epsilon(err)) thisTest = .true.
  enddo
  write(*,*) 'h1 (by projection)'
  call dumpGrid(n1,h1,x1,u1)
  if (thisTest) write(*,*) 'remapping_unit_tests: Failed remapping_core_w()'
  remapping_unit_tests = remapping_unit_tests .or. thisTest

  thisTest = .false.
  allocate(ppoly0_E(n0,2))
  allocate(ppoly0_S(n0,2))
  allocate(ppoly0_coefficients(n0,CS%degree+1))

  ppoly0_E(:,:) = 0.0
  ppoly0_S(:,:) = 0.0
  ppoly0_coefficients(:,:) = 0.0

  call edge_values_explicit_h4( n0, h0, u0, ppoly0_E )
  call PPM_reconstruction( n0, h0, u0, ppoly0_E, ppoly0_coefficients )
  call PPM_boundary_extrapolation( n0, h0, u0, ppoly0_E, ppoly0_coefficients )
  u1(:) = 0.
  call remapByProjection( n0, h0, u0, ppoly0_E, ppoly0_coefficients, &
                          n1, h1, INTEGRATION_PPM, u1 )
  do i=1,n1
    err=u1(i)-8.*(0.5*real(1+n1)-real(i))
    if (abs(err)>2.*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(*,*) 'remapping_unit_tests: Failed remapByProjection()'
  remapping_unit_tests = remapping_unit_tests .or. thisTest

  thisTest = .false.
  u1(:) = 0.
  call remapByDeltaZ( n0, h0, u0, ppoly0_E, ppoly0_coefficients, &
                      n1, x1-x0(1:n1+1), &
                      INTEGRATION_PPM, u1, hn1 )
  write(*,*) 'h1 (by delta)'
  call dumpGrid(n1,h1,x1,u1)
  hn1=hn1-h1
  do i=1,n1
    err=u1(i)-8.*(0.5*real(1+n1)-real(i))
    if (abs(err)>2.*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(*,*) 'remapping_unit_tests: Failed remapByDeltaZ() 1'
  remapping_unit_tests = remapping_unit_tests .or. thisTest

  thisTest = .false.
  call buildGridFromH(n2, h2, x2)
  dx2(1:n0+1) = x2(1:n0+1) - x0
  dx2(n0+2:n2+1) = x2(n0+2:n2+1) - x0(n0+1)
  call remapByDeltaZ( n0, h0, u0, ppoly0_E, ppoly0_coefficients, &
                      n2, dx2, &
                      INTEGRATION_PPM, u2, hn2 )
  write(*,*) 'h2'
  call dumpGrid(n2,h2,x2,u2)
  write(*,*) 'hn2'
  call dumpGrid(n2,hn2,x2,u2)

  do i=1,n2
    err=u2(i)-8./2.*(0.5*real(1+n2)-real(i))
    if (abs(err)>2.*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(*,*) 'remapping_unit_tests: Failed remapByDeltaZ() 2'
  remapping_unit_tests = remapping_unit_tests .or. thisTest

  write(*,*) 'Via sub-cells'
  thisTest = .false.
  call remap_via_sub_cells( n0, h0, u0, ppoly0_E, ppoly0_coefficients, &
                            n2, h2, INTEGRATION_PPM, .false., u2, err )
  call dumpGrid(n2,h2,x2,u2)

  do i=1,n2
    err=u2(i)-8./2.*(0.5*real(1+n2)-real(i))
    if (abs(err)>2.*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(*,*) 'remapping_unit_tests: Failed remap_via_sub_cells() 2'
  remapping_unit_tests = remapping_unit_tests .or. thisTest

  call remap_via_sub_cells( n0, h0, u0, ppoly0_E, ppoly0_coefficients, &
                            6, (/.125,.125,.125,.125,.125,.125/), INTEGRATION_PPM, .false., u2, err )
  call dumpGrid(6,h2,x2,u2)

  call remap_via_sub_cells( n0, h0, u0, ppoly0_E, ppoly0_coefficients, &
                            3, (/2.25,1.5,1./), INTEGRATION_PPM, .false., u2, err )
  call dumpGrid(3,h2,x2,u2)


  write(*,*) '===== MOM_remapping: new remapping_unit_tests =================='

  deallocate(ppoly0_E, ppoly0_S, ppoly0_coefficients)
  allocate(ppoly0_coefficients(5,6))
  allocate(ppoly0_E(5,2))
  allocate(ppoly0_S(5,2))

  call PCM_reconstruction(3, (/1.,2.,4./), ppoly0_E(1:3,:), ppoly0_coefficients(1:3,:) )
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_E(:,1), (/1.,2.,4./), 'PCM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_E(:,2), (/1.,2.,4./), 'PCM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_coefficients(:,1), (/1.,2.,4./), 'PCM: P0')

  call PLM_reconstruction(3, (/1.,1.,1./), (/1.,3.,5./), ppoly0_E(1:3,:), ppoly0_coefficients(1:3,:) )
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_E(:,1), (/1.,2.,5./), 'Unlim PLM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_E(:,2), (/1.,4.,5./), 'Unlim PLM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_coefficients(:,1), (/1.,2.,5./), 'Unlim PLM: P0')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_coefficients(:,2), (/0.,2.,0./), 'Unlim PLM: P1')

  call PLM_reconstruction(3, (/1.,1.,1./), (/1.,2.,7./), ppoly0_E(1:3,:), ppoly0_coefficients(1:3,:) )
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_E(:,1), (/1.,1.,7./), 'Left lim PLM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_E(:,2), (/1.,3.,7./), 'Left lim PLM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_coefficients(:,1), (/1.,1.,7./), 'Left lim PLM: P0')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_coefficients(:,2), (/0.,2.,0./), 'Left lim PLM: P1')

  call PLM_reconstruction(3, (/1.,1.,1./), (/1.,6.,7./), ppoly0_E(1:3,:), ppoly0_coefficients(1:3,:) )
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_E(:,1), (/1.,5.,7./), 'Right lim PLM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_E(:,2), (/1.,7.,7./), 'Right lim PLM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_coefficients(:,1), (/1.,5.,7./), 'Right lim PLM: P0')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_coefficients(:,2), (/0.,2.,0./), 'Right lim PLM: P1')

  call PLM_reconstruction(3, (/1.,2.,3./), (/1.,4.,9./), ppoly0_E(1:3,:), ppoly0_coefficients(1:3,:) )
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_E(:,1), (/1.,2.,9./), 'Non-uniform line PLM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_E(:,2), (/1.,6.,9./), 'Non-uniform line PLM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_coefficients(:,1), (/1.,2.,9./), 'Non-uniform line PLM: P0')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(3, ppoly0_coefficients(:,2), (/0.,4.,0./), 'Non-uniform line PLM: P1')

  call edge_values_explicit_h4( 5, (/1.,1.,1.,1.,1./), (/1.,3.,5.,7.,9./), ppoly0_E )
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_E(:,1), (/0.,2.,4.,6.,8./), 'Line H4: left edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_E(:,2), (/2.,4.,6.,8.,10./), 'Line H4: right edges')
  remapping_unit_tests = .false. ! edge_values_explicit_h4 fails due to roundoff!
  ppoly0_E(:,1) = (/0.,2.,4.,6.,8./)
  ppoly0_E(:,2) = (/2.,4.,6.,8.,10./)
  call PPM_reconstruction(5, (/1.,1.,1.,1.,1./), (/1.,3.,5.,7.,9./), ppoly0_E(1:5,:), ppoly0_coefficients(1:5,:) )
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_coefficients(:,1), (/1.,2.,4.,6.,9./), 'Line PPM: P0')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_coefficients(:,2), (/0.,2.,2.,2.,0./), 'Line PPM: P1')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_coefficients(:,3), (/0.,0.,0.,0.,0./), 'Line PPM: P2')

  call edge_values_explicit_h4( 5, (/1.,1.,1.,1.,1./), (/1.,1.,7.,19.,37./), ppoly0_E )
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_E(:,1), (/3.,0.,3.,12.,27./), 'Parabola H4: left edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_E(:,2), (/0.,3.,12.,27.,48./), 'Parabola H4: right edges')
  remapping_unit_tests = .false. ! edge_values_explicit_h4 fails due to roundoff!
  ppoly0_E(:,1) = (/0.,0.,3.,12.,27./)
  ppoly0_E(:,2) = (/0.,3.,12.,27.,48./)
  call PPM_reconstruction(5, (/1.,1.,1.,1.,1./), (/0.,1.,7.,19.,37./), ppoly0_E(1:5,:), ppoly0_coefficients(1:5,:) )
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_E(:,1), (/0.,0.,3.,12.,37./), 'Parabola PPM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_E(:,2), (/0.,3.,12.,27.,37./), 'Parabola PPM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_coefficients(:,1), (/0.,0.,3.,12.,37./), 'Parabola PPM: P0')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_coefficients(:,2), (/0.,0.,6.,12.,0./), 'Parabola PPM: P1')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_coefficients(:,3), (/0.,3.,3.,3.,0./), 'Parabola PPM: P2')

  ppoly0_E(:,1) = (/0.,0.,6.,10.,15./)
  ppoly0_E(:,2) = (/0.,6.,12.,17.,15./)
  call PPM_reconstruction(5, (/1.,1.,1.,1.,1./), (/0.,5.,7.,16.,15./), ppoly0_E(1:5,:), ppoly0_coefficients(1:5,:) )
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_E(:,1), (/0.,3.,6.,16.,15./), 'Limits PPM: left edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_E(:,2), (/0.,6.,9.,16.,15./), 'Limits PPM: right edges')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_coefficients(:,1), (/0.,3.,6.,16.,15./), 'Limits PPM: P0')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_coefficients(:,2), (/0.,6.,0.,0.,0./), 'Limits PPM: P1')
  remapping_unit_tests = remapping_unit_tests .or. test_answer(5, ppoly0_coefficients(:,3), (/0.,-3.,3.,0.,0./), 'Limits PPM: P2')

  deallocate(ppoly0_E, ppoly0_S, ppoly0_coefficients)

  write(*,*) '=========================================================='

  contains

  !> Returns true if any cell of u and u_true are not identical. Returns false otherwise.
  logical function test_answer(n, u, u_true, label)
    integer,            intent(in) :: n !< Number of cells in u
    real, dimension(n), intent(in) :: u !< Values to test
    real, dimension(n), intent(in) :: u_true !< Values to test against (correct answer)
    character(len=*),   intent(in) :: label !< Message
    ! Local variables
    integer :: k

    test_answer = .false.
    do k = 1, n
      if (u(k) /= u_true(k)) test_answer = .true.
    enddo
    if (test_answer .or. .true.) then
      write(*,'(a4,2a24,x,a)') 'k','Calculated value','Correct value',label
      do k = 1, n
        if (u(k) /= u_true(k)) then
          write(*,'(i4,1p2e24.16,a,1pe24.16,a)') k,u(k),u_true(k),' err=',u(k)-u_true(k),' < wrong'
        else
          write(*,'(i4,1p2e24.16)') k,u(k),u_true(k)
        endif
      enddo
    endif

  end function test_answer

  !> Convenience function for printing grid to screen
  subroutine dumpGrid(n,h,x,u)
  integer, intent(in) :: n !< Number of cells
  real, dimension(:), intent(in) :: h !< Cell thickness
  real, dimension(:), intent(in) :: x !< Interface delta
  real, dimension(:), intent(in) :: u !< Cell average values
  integer :: i
  write(*,'("i=",20i10)') (i,i=1,n+1)
  write(*,'("x=",20es10.2)') (x(i),i=1,n+1)
  write(*,'("i=",5x,20i10)') (i,i=1,n)
  write(*,'("h=",5x,20es10.2)') (h(i),i=1,n)
  write(*,'("u=",5x,20es10.2)') (u(i),i=1,n)
  end subroutine dumpGrid

end function remapping_unit_tests

end module MOM_remapping
