!> Retains older versions of column-wise vertical remapping functions that are
!! no longer used in MOM6, but may be useful later for documenting the development
!! of the schemes that are used in MOM6.
module remapping_attic

! This file is part of MOM6. See LICENSE.md for the license.
! Original module written by Laurent White, 2008.06.09

use MOM_error_handler,  only : MOM_error, FATAL
use MOM_io,             only : stdout
use PPM_functions,      only : PPM_reconstruction, PPM_boundary_extrapolation
use regrid_edge_values, only : edge_values_explicit_h4

implicit none ; private

! The following routines are visible to the outside world
public remapping_attic_unit_tests, remapByProjection, remapByDeltaZ
public isPosSumErrSignificant

! The following are private parameter constants
integer, parameter  :: INTEGRATION_PCM = 0  !< Piecewise Constant Method
integer, parameter  :: INTEGRATION_PLM = 1  !< Piecewise Linear Method
integer, parameter  :: INTEGRATION_PPM = 3  !< Piecewise Parabolic Method
integer, parameter  :: INTEGRATION_PQM = 5  !< Piecewise Quartic Method

! This CPP macro turns on/off bounding of integrations limits so that they are
! always within the cell. Roundoff can lead to the non-dimensional bounds being
! outside of the range 0 to 1.
#define __USE_ROUNDOFF_SAFE_ADJUSTMENTS__

contains

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
  real,    intent(in) :: sum1 !< Sum of n1 values in arbitrary units [A]
  real,    intent(in) :: sum2 !< Sum of n2 values [A]
  logical             :: isPosSumErrSignificant !< True if difference in sums is large
  ! Local variables
  real :: sumErr      ! The absolutde difference in the sums [A]
  real :: allowedErr  ! The tolerance for the integrated reconstruction [A]
  real :: eps         ! A tiny fractional error [nondim]

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

!> Remaps column of values u0 on grid h0 to grid h1 by integrating
!! over the projection of each h1 cell onto the h0 grid.
subroutine remapByProjection( n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                              n1, h1, method, u1, h_neglect )
  integer,       intent(in)    :: n0     !< Number of cells in source grid
  real,          intent(in)    :: h0(:)  !< Source grid widths (size n0) in thickness units [H]
  real,          intent(in)    :: u0(:)  !< Source cell averages (size n0) in arbitrary units [A]
  real,          intent(in)    :: ppoly0_E(:,:)     !< Edge value of polynomial [A]
  real,          intent(in)    :: ppoly0_coefs(:,:) !< Coefficients of polynomial [A]
  integer,       intent(in)    :: n1     !< Number of cells in target grid
  real,          intent(in)    :: h1(:)  !< Target grid widths (size n1) [H]
  integer,       intent(in)    :: method !< Remapping scheme to use
  real,          intent(out)   :: u1(:)  !< Target cell averages (size n1) [A]
  real,          intent(in)    :: h_neglect !< A negligibly small width for the
                                           !! purpose of cell reconstructions
                                           !! in the same units as h [H].
  ! Local variables
  integer       :: iTarget
  real          :: xL, xR       ! coordinates of target cell edges [H]
  integer       :: jStart ! Used by integrateReconOnInterval()
  real          :: xStart ! Used by integrateReconOnInterval() [H]

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

    call integrateReconOnInterval( n0, h0, u0, ppoly0_E, ppoly0_coefs, method, &
                                   xL, xR, h1(iTarget), u1(iTarget), jStart, xStart, h_neglect )

  enddo ! end iTarget loop on target grid cells

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
subroutine remapByDeltaZ( n0, h0, u0, ppoly0_E, ppoly0_coefs, n1, dx1, &
                          method, u1, h1, h_neglect )
  integer,              intent(in)  :: n0     !< Number of cells in source grid
  real, dimension(:),   intent(in)  :: h0     !< Source grid sizes (size n0) in thickness units [H]
  real, dimension(:),   intent(in)  :: u0     !< Source cell averages (size n0) in arbitrary units [A]
  real, dimension(:,:), intent(in)  :: ppoly0_E !< Edge value of polynomial [A]
  real, dimension(:,:), intent(in)  :: ppoly0_coefs !< Coefficients of polynomial [A]
  integer,              intent(in)  :: n1     !< Number of cells in target grid
  real, dimension(:),   intent(in)  :: dx1    !< Target grid edge positions (size n1+1) [H]
  integer,              intent(in)  :: method !< Remapping scheme to use
  real, dimension(:),   intent(out) :: u1     !< Target cell averages (size n1) [A]
  real, dimension(:), &
              optional, intent(out) :: h1     !< Target grid widths (size n1) [H]
  real,                 intent(in)  :: h_neglect !< A negligibly small width for the
                                           !! purpose of cell reconstructions
                                           !! in the same units as h [H].
  ! Local variables
  integer :: iTarget
  real    :: xL, xR     ! Coordinates of target cell edges [H]
  real    :: xOld, xNew ! Edge positions on the old and new grids [H]
  real    :: hOld, hNew ! Cell thicknesses on the old and new grids [H]
  real    :: uOld       ! A source cell average of u [A]
  real    :: h_err      ! An estimate of the error in the reconstructed thicknesses [H]
  real    :: uhNew      ! Cell integrated u on the new grid [A H]
  real    :: hFlux      ! Width of the remapped volume [H]
  real    :: uAve       ! Target cell average of u [A]
  real    :: fluxL, fluxR ! Fluxes of u through the two cell faces [A H]
  integer :: jStart ! Used by integrateReconOnInterval()
  real    :: xStart ! Used by integrateReconOnInterval() [H]

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
    call integrateReconOnInterval( n0, h0, u0, ppoly0_E, ppoly0_coefs, method, &
                                   xL, xR, hFlux, uAve, jStart, xStart, h_neglect )
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

  enddo ! end iTarget loop on target grid cells

end subroutine remapByDeltaZ

!> Integrate the reconstructed column profile over a single cell
subroutine integrateReconOnInterval( n0, h0, u0, ppoly0_E, ppoly0_coefs, method, &
                                     xL, xR, hC, uAve, jStart, xStart, h_neglect )
  integer,              intent(in)    :: n0     !< Number of cells in source grid
  real, dimension(:),   intent(in)    :: h0     !< Source grid sizes (size n0) in thickness units [H]
  real, dimension(:),   intent(in)    :: u0     !< Source cell averages in arbitrary units [A]
  real, dimension(:,:), intent(in)    :: ppoly0_E !< Edge value of polynomial [A]
  real, dimension(:,:), intent(in)    :: ppoly0_coefs !< Coefficients of polynomial [A]
  integer,              intent(in)    :: method !< Remapping scheme to use
  real,                 intent(in)    :: xL     !< Left edges of target cell [H]
  real,                 intent(in)    :: xR     !< Right edges of target cell [H]
  real,                 intent(in)    :: hC     !< Cell width hC = xR - xL [H]
  real,                 intent(out)   :: uAve   !< Average value on target cell [A]
  integer,              intent(inout) :: jStart !< The index of the cell to start searching from
                                   !< On exit, contains index of last cell used
  real,                 intent(inout) :: xStart !< The left edge position of cell jStart [H]
                                   !< On first entry should be 0.
  real,                 intent(in)    :: h_neglect !< A negligibly small width for the
                                          !! purpose of cell reconstructions
                                          !! in the same units as h [H]
  ! Local variables
  integer :: j, k
  integer :: jL, jR       ! indexes of source cells containing target cell edges
  real    :: q            ! complete integration [A H]
  real    :: xi0, xi1     ! interval of integration (local -- normalized -- coordinates) [nondim]
  real    :: x0jLl, x0jLr ! Left/right position of cell jL [H]
  real    :: x0jRl, x0jRr ! Left/right position of cell jR [H]
  real    :: hAct         ! The distance actually used in the integration
                          ! (notionally xR - xL) which differs due to roundoff [H].
  real    :: x0_2, x1_2   ! Squares of normalized positions used to evaluate polynomials [nondim]
  real    :: x0px1, x02px12 ! Sums of normalized positions and their squares [nondim]
  real, parameter :: r_3 = 1.0/3.0 ! Used in evaluation of integrated polynomials [nondim]

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
          uAve =         ppoly0_coefs(jL,1)
        case ( INTEGRATION_PLM )
          uAve =         ppoly0_coefs(jL,1)   &
               + xi0 *   ppoly0_coefs(jL,2)
        case ( INTEGRATION_PPM )
          uAve =         ppoly0_coefs(jL,1)   &
               + xi0 * ( ppoly0_coefs(jL,2)   &
               + xi0 *   ppoly0_coefs(jL,3) )
        case ( INTEGRATION_PQM )
          uAve =         ppoly0_coefs(jL,1)   &
               + xi0 * ( ppoly0_coefs(jL,2)   &
               + xi0 * ( ppoly0_coefs(jL,3)   &
               + xi0 * ( ppoly0_coefs(jL,4)   &
               + xi0 *   ppoly0_coefs(jL,5) ) ) )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select

    endif ! end checking whether source cell is vanished

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
          q = ( xR - xL ) * ppoly0_coefs(jL,1)
        case ( INTEGRATION_PLM )
          q = ( xR - xL ) * (                            &
              ppoly0_coefs(jL,1)                         &
            + ppoly0_coefs(jL,2) * 0.5 * ( xi1 + xi0 ) )
        case ( INTEGRATION_PPM )
          q = ( xR - xL ) * (                            &
                ppoly0_coefs(jL,1)                       &
            + ( ppoly0_coefs(jL,2) * 0.5 * ( xi1 + xi0 ) &
            +   ppoly0_coefs(jL,3) * r_3 * ( ( xi1*xi1 + xi0*xi0 ) + xi0*xi1 ) ) )
        case ( INTEGRATION_PQM )
          x0_2 = xi0*xi0
          x1_2 = xi1*xi1
          x02px12 = x0_2 + x1_2
          x0px1 = xi1 + xi0
          q = ( xR - xL ) * (                                    &
                ppoly0_coefs(jL,1)                               &
            + ( ppoly0_coefs(jL,2) * 0.5 * ( xi1 + xi0 )         &
            + ( ppoly0_coefs(jL,3) * r_3 * ( x02px12 + xi0*xi1 ) &
            +   ppoly0_coefs(jL,4) * 0.25* ( x02px12 * x0px1 )   &
            +   ppoly0_coefs(jL,5) * 0.2 * ( ( xi1*x1_2 + xi0*x0_2 ) * x0px1 + x0_2*x1_2 ) ) ) )
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
          q = q + ( x0jLr - xL ) * ppoly0_coefs(jL,1)
        case ( INTEGRATION_PLM )
          q = q + ( x0jLr - xL ) * (                     &
              ppoly0_coefs(jL,1)                         &
            + ppoly0_coefs(jL,2) * 0.5 * ( xi1 + xi0 ) )
        case ( INTEGRATION_PPM )
          q = q + ( x0jLr - xL ) * (                     &
                ppoly0_coefs(jL,1)                       &
            + ( ppoly0_coefs(jL,2) * 0.5 * ( xi1 + xi0 ) &
            +   ppoly0_coefs(jL,3) * r_3 * ( ( xi1*xi1 + xi0*xi0 ) + xi0*xi1 ) ) )
        case ( INTEGRATION_PQM )
          x0_2 = xi0*xi0
          x1_2 = xi1*xi1
          x02px12 = x0_2 + x1_2
          x0px1 = xi1 + xi0
          q = q + ( x0jLr - xL ) * (                             &
                ppoly0_coefs(jL,1)                               &
            + ( ppoly0_coefs(jL,2) * 0.5 * ( xi1 + xi0 )         &
            + ( ppoly0_coefs(jL,3) * r_3 * ( x02px12 + xi0*xi1 ) &
            +   ppoly0_coefs(jL,4) * 0.25* ( x02px12 * x0px1 )   &
            +   ppoly0_coefs(jL,5) * 0.2 * ( ( xi1*x1_2 + xi0*x0_2 ) * x0px1 + x0_2*x1_2 ) ) ) )
        case default
          call MOM_error( FATAL, 'The selected integration method is invalid' )
      end select

      ! Integrate contents within cells strictly comprised between jL and jR
      if ( jR > (jL+1) ) then
        do k = jL+1,jR-1
          q = q + h0(k) * u0(k)
          hAct = hAct + h0(k)
        enddo
      endif

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
          q = q + ( xR - x0jRl ) * ppoly0_coefs(jR,1)
        case ( INTEGRATION_PLM )
          q = q + ( xR - x0jRl ) * (                     &
              ppoly0_coefs(jR,1)                         &
            + ppoly0_coefs(jR,2) * 0.5 * ( xi1 + xi0 ) )
        case ( INTEGRATION_PPM )
          q = q + ( xR - x0jRl ) * (                     &
                ppoly0_coefs(jR,1)                       &
            + ( ppoly0_coefs(jR,2) * 0.5 * ( xi1 + xi0 ) &
            +   ppoly0_coefs(jR,3) * r_3 * ( ( xi1*xi1 + xi0*xi0 ) + xi0*xi1 ) ) )
        case ( INTEGRATION_PQM )
          x0_2 = xi0*xi0
          x1_2 = xi1*xi1
          x02px12 = x0_2 + x1_2
          x0px1 = xi1 + xi0
          q = q + ( xR - x0jRl ) * (                             &
                ppoly0_coefs(jR,1)                               &
            + ( ppoly0_coefs(jR,2) * 0.5 * ( xi1 + xi0 )         &
            + ( ppoly0_coefs(jR,3) * r_3 * ( x02px12 + xi0*xi1 ) &
            +   ppoly0_coefs(jR,4) * 0.25* ( x02px12 * x0px1 )   &
            +   ppoly0_coefs(jR,5) * 0.2 * ( ( xi1*x1_2 + xi0*x0_2 ) * x0px1 + x0_2*x1_2 ) ) ) )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select

    endif ! end integration for non-vanished cells

    ! The cell average is the integrated value divided by the cell width
#ifdef __USE_ROUNDOFF_SAFE_ADJUSTMENTS__
if (hAct==0.) then
    uAve = ppoly0_coefs(jL,1)
else
    uAve = q / hAct
endif
#else
    uAve = q / hC
#endif

  endif ! endif clause to check if cell is vanished

end subroutine integrateReconOnInterval

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

!> Runs unit tests on archaic remapping functions.
!! Should only be called from a single/root thread
!! Returns True if a test fails, otherwise False
logical function remapping_attic_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables
  integer, parameter :: n0 = 4, n1 = 3, n2 = 6
  real :: h0(n0), x0(n0+1) ! Test cell widths and edge coordinates [H]
  real :: u0(n0)           ! Test values for remapping in arbitrary units [A]
  real :: h1(n1), x1(n1+1) ! Test cell widths and edge coordinates [H]
  real :: u1(n1)           ! Test values for remapping [A]
  real :: h2(n2), x2(n2+1) ! Test cell widths and edge coordinates [H]
  real :: u2(n2)           ! Test values for remapping [A]
  real :: hn1(n1), hn2(n2) ! Updated grid thicknesses [H]
  real :: dx1(n1+1), dx2(n2+1) ! Differences in interface positions [H]
  data u0 /9., 3., -3., -9./   ! Linear profile, 4 at surface to -4 at bottom
  data h0 /4*0.75/ ! 4 uniform layers with total depth of 3
  data h1 /3*1./   ! 3 uniform layers with total depth of 3
  data h2 /6*0.5/  ! 6 uniform layers with total depth of 3
  real, allocatable, dimension(:,:) :: ppoly0_E, ppoly0_S ! Polynomial edge values [A]
  real, allocatable, dimension(:,:) :: ppoly0_coefs ! Polynomial reconstruction coefficients [A]
  integer :: answer_date  ! The vintage of the expressions to test
  integer :: i, degree
  real :: err  ! Difference between a remapped value and its expected value [A]
  real :: h_neglect, h_neglect_edge ! Negligible thicknesses used in remapping [H]
  logical :: thisTest, v

  v = verbose
  answer_date = 20190101 ! 20181231
  h_neglect = 1.0E-30
  h_neglect_edge = h_neglect ; if (answer_date < 20190101) h_neglect_edge = 1.0e-10

  write(stdout,*) '==== remapping_attic: remapping_attic_unit_tests ================='
  remapping_attic_unit_tests = .false. ! Normally return false

  call buildGridFromH(n0, h0, x0)
  call buildGridFromH(n1, h1, x1)

  thisTest = .false.
  degree = 2
  if (verbose) write(stdout,*) 'h0 (test data)'
  if (verbose) call dumpGrid(n0,h0,x0,u0)

  call dzFromH1H2( n0, h0, n1, h1, dx1 )

  thisTest = .false.
  allocate(ppoly0_E(n0,2))
  allocate(ppoly0_S(n0,2))
  allocate(ppoly0_coefs(n0,degree+1))

  ppoly0_E(:,:) = 0.0
  ppoly0_S(:,:) = 0.0
  ppoly0_coefs(:,:) = 0.0

  call edge_values_explicit_h4( n0, h0, u0, ppoly0_E, h_neglect=1e-10, answer_date=answer_date )
  call PPM_reconstruction( n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect, answer_date=answer_date )
  call PPM_boundary_extrapolation( n0, h0, u0, ppoly0_E, ppoly0_coefs, h_neglect )
  u1(:) = 0.
  call remapByProjection( n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                          n1, h1, INTEGRATION_PPM, u1, h_neglect )
  do i=1,n1
    err = u1(i)-8.*(0.5*real(1+n1)-real(i))
    if (abs(err)>2.*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(stdout,*) 'remapping_attic_unit_tests: Failed remapByProjection()'
  remapping_attic_unit_tests = remapping_attic_unit_tests .or. thisTest

  thisTest = .false.
  u1(:) = 0.
  call remapByDeltaZ( n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                      n1, x1-x0(1:n1+1), &
                      INTEGRATION_PPM, u1, hn1, h_neglect )
  if (verbose) write(stdout,*) 'h1 (by delta)'
  if (verbose) call dumpGrid(n1,h1,x1,u1)
  hn1 = hn1-h1
  do i=1,n1
    err = u1(i)-8.*(0.5*real(1+n1)-real(i))
    if (abs(err)>2.*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(stdout,*) 'remapping_attic_unit_tests: Failed remapByDeltaZ() 1'
  remapping_attic_unit_tests = remapping_attic_unit_tests .or. thisTest

  thisTest = .false.
  call buildGridFromH(n2, h2, x2)
  dx2(1:n0+1) = x2(1:n0+1) - x0
  dx2(n0+2:n2+1) = x2(n0+2:n2+1) - x0(n0+1)
  call remapByDeltaZ( n0, h0, u0, ppoly0_E, ppoly0_coefs, &
                      n2, dx2, &
                      INTEGRATION_PPM, u2, hn2, h_neglect )
  if (verbose) write(stdout,*) 'h2'
  if (verbose) call dumpGrid(n2,h2,x2,u2)
  if (verbose) write(stdout,*) 'hn2'
  if (verbose) call dumpGrid(n2,hn2,x2,u2)

  do i=1,n2
    err = u2(i)-8./2.*(0.5*real(1+n2)-real(i))
    if (abs(err)>2.*epsilon(err)) thisTest = .true.
  enddo
  if (thisTest) write(stdout,*) 'remapping_attic_unit_tests: Failed remapByDeltaZ() 2'
  remapping_attic_unit_tests = remapping_attic_unit_tests .or. thisTest

  if (.not. remapping_attic_unit_tests) write(stdout,*) 'Pass'

end function remapping_attic_unit_tests

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

end module remapping_attic
