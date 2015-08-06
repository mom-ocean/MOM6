module MOM_remapping
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.09
! L. White
!
! This module contains the main remapping routines. 
!
!==============================================================================
use MOM_error_handler, only : MOM_error, FATAL
use MOM_string_functions, only : uppercase
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use polynomial_functions, only : evaluation_polynomial, integration_polynomial
use regrid_edge_values, only : edge_values_explicit_h4, edge_values_implicit_h4
use regrid_edge_values, only : edge_values_implicit_h4, edge_values_implicit_h6
use regrid_edge_slopes, only : edge_slopes_implicit_h3, edge_slopes_implicit_h5
use PCM_functions, only : PCM_reconstruction
use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation
use PQM_functions, only : PQM_reconstruction, PQM_boundary_extrapolation_v1

implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Container for private (module-wise) variables and parameters
! -----------------------------------------------------------------------------
type, public :: remapping_CS
  private
  ! Parameters
  integer :: nk = 0                    ! Number of layers/levels in vertical
  integer :: remapping_scheme = -911   ! Determines which reconstruction to use
  integer :: degree=0                  ! Degree of polynomial reconstruction
  logical :: boundary_extrapolation = .true.  ! If true, extrapolate boundaries
end type

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public remapping_main, remapping_core
public initialize_remapping, end_remapping
public remapEnableBoundaryExtrapolation, remapDisableBoundaryExtrapolation
public setReconstructionType
public remappingUnitTests
public dzFromH1H2

! -----------------------------------------------------------------------------
! The following are private parameter constants
! -----------------------------------------------------------------------------
! List of remapping schemes
integer, parameter  :: REMAPPING_PCM        = 0 ! O(h^1)
integer, parameter  :: REMAPPING_PLM        = 1 ! O(h^2)
integer, parameter  :: REMAPPING_PPM_H4     = 2 ! O(h^3)
integer, parameter  :: REMAPPING_PPM_IH4    = 3 ! O(h^3)
integer, parameter  :: REMAPPING_PQM_IH4IH3 = 4 ! O(h^4)
integer, parameter  :: REMAPPING_PQM_IH6IH5 = 5 ! O(h^5)

! These control what routine to use for the remapping integration
integer, parameter  :: INTEGRATION_PCM = 0  ! scope: global
integer, parameter  :: INTEGRATION_PLM = 1  ! scope: global
integer, parameter  :: INTEGRATION_PPM = 3  ! scope: global
integer, parameter  :: INTEGRATION_PQM = 5  ! scope: global

character(len=40)  :: mod = "MOM_remapping" ! This module's name.

! Documentation for external callers
character(len=256), public :: remappingSchemesDoc = &
                 "PCM         (1st-order accurate)\n"//&
                 "PLM         (2nd-order accurate)\n"//&
                 "PPM_H4      (3rd-order accurate)\n"//&
                 "PPM_IH4     (3rd-order accurate)\n"//&
                 "PQM_IH4IH3  (4th-order accurate)\n"//&
                 "PQM_IH6IH5  (5th-order accurate)\n"
character(len=3), public :: remappingDefaultScheme = "PLM"

! This CPP macro embeds some safety checks
#define __DO_SAFETY_CHECKS__

! This CPP macro turns on/off bounding of integrations limits so that they are
! always within the cell. Roundoff can lead to the non-dimensional bounds being
! outside of the range 0 to 1.
#define __USE_ROUNDOFF_SAFE_ADJUSTMENTS__

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! General remapping routine 
!------------------------------------------------------------------------------
subroutine remapping_main( CS, G, h, dxInterface, tv, u, v )
!------------------------------------------------------------------------------
! This routine takes care of remapping all variable between the old and the
! new grids. When velocity components need to be remapped, thicknesses at
! velocity points are taken to be arithmetic averages of tracer thicknesses.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(remapping_CS),                               intent(in)    :: CS
  type(ocean_grid_type),                            intent(in)    :: G
  real, dimension(NIMEM_,NJMEM_,NKMEM_),            intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_),     intent(in)    :: dxInterface
  type(thermo_var_ptrs),                            intent(inout) :: tv       
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), optional, intent(inout) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), optional, intent(inout) :: v
  
  ! Local variables
  integer               :: i, j, k
  integer               :: nz
  real, dimension(G%ke+1) :: dx
  real, dimension(G%ke) :: h1, u_column

  nz = G%ke

  ! Remap tracer
!$OMP parallel default(none) shared(G,h,dxInterface,CS,nz,tv,u,v) &
!$OMP                       private(h1,dx,u_column)
  if (associated(tv%S)) then ! Assume T and S are either both associated or both not
!$OMP do
    do j = G%jsc,G%jec
      do i = G%isc,G%iec
        if (G%mask2dT(i,j)>0.) then
          ! Build the start and final grids
          h1(:) = h(i,j,:)
          dx(:) = dxInterface(i,j,:)
          call remapping_core(CS, nz, h1, tv%S(i,j,:), nz, dx, u_column)
          tv%S(i,j,:) = u_column(:)
          call remapping_core(CS, nz, h1, tv%T(i,j,:), nz, dx, u_column)
          tv%T(i,j,:) = u_column(:)
        endif
      enddo
    enddo
  endif
  
  ! Remap u velocity component
  if ( present(u) ) then
!$OMP do
    do j = G%jsc,G%jec
      do i = G%iscB,G%iecB
        if (G%mask2dCu(i,j)>0.) then
          ! Build the start and final grids
          h1(:) = 0.5 * ( h(i,j,:) + h(i+1,j,:) )
          dx(:) = 0.5 * ( dxInterface(i,j,:) + dxInterface(i+1,j,:) )
          call remapping_core(CS, nz, h1, u(i,j,:), nz, dx, u_column)
          u(i,j,:) = u_column(:)
        endif
      enddo
    enddo
  endif
  
  ! Remap v velocity component
  if ( present(v) ) then
!$OMP do
    do j = G%jscB,G%jecB
      do i = G%isc,G%iec
        if (G%mask2dCv(i,j)>0.) then
          ! Build the start and final grids
          h1(:) = 0.5 * ( h(i,j,:) + h(i,j+1,:) )
          dx(:) = 0.5 * ( dxInterface(i,j,:) + dxInterface(i,j+1,:) )
          call remapping_core(CS, nz, h1, v(i,j,:), nz, dx, u_column)
          v(i,j,:) = u_column(:)
        endif
      enddo
    enddo
  endif
!$OMP end parallel

end subroutine remapping_main


!------------------------------------------------------------------------------
! Build a grid from h
!------------------------------------------------------------------------------
subroutine buildGridFromH(nz, h, x)
!------------------------------------------------------------------------------
! This routine calculates the coordinates x by integrating h.
!------------------------------------------------------------------------------

  ! Arguments
  integer,               intent(in)    :: nz
  real, dimension(nz),   intent(in)    :: h
  real, dimension(nz+1), intent(inout) :: x

  integer :: k

  ! Build start grid
  x(1) = 0.0
  do k = 1,nz
    x(k+1) = x(k) + h(k)
  end do

end subroutine buildGridFromH


!------------------------------------------------------------------------------
! Check that two grids are consistent
!------------------------------------------------------------------------------
subroutine checkConsistantCoords(ns, xs, nf, xf, strict, msg)
!------------------------------------------------------------------------------
! Checks that xs and xf are consistent to within roundoff.
! If strict=False, the end points of xs and xf are allowed to differ by
! numerical roundoff due to the nature of summation to obtain xs, xf.
! If strict=True, the end points must be identical.
!------------------------------------------------------------------------------

  ! Arguments
  integer, intent(in) :: ns, nf
  real,    intent(in) :: xs(ns+1), xf(nf+1)
  logical, intent(in) :: strict
  character(len=*), intent(in) :: msg

  ! Local variables
  integer :: k
  real    :: sumHs, sumHf

  sumHs = xs(ns+1)-xs(1)
  sumHf = xf(nf+1)-xf(1)

  if (strict) then
    if (sumHf /= sumHs) call &
        MOM_error(FATAL,'MOM_remapping, checkConsistantCoords: '//&
                        'Total thickness of two grids are not exactly equal.'//&
                        ' Called from '//trim(msg) )
  else ! not strict
    if (isPosSumErrSignificant(ns, sumHs, nf, sumhf)) then
      write(0,*) 'Start/final/start-final grid'
      do k = 1,max(ns,nf)+1
        if (k<=min(ns+1,nf+1)) then
          write(0,'(i4,3es12.3)') k,xs(k),xf(k),xs(k)-xf(k)
        elseif (k>ns+1) then
          write(0,'(i4,12x,1es12.3)') k,xf(k)
        else
          write(0,'(i4,1es12.3)') k,xs(k)
        endif
      enddo
      call MOM_error(FATAL,'MOM_remapping, checkConsistantCoords: '//&
              'Total thickness of two grids do not match to within round-off.'//&
              ' Called from '//trim(msg) )
    endif
  endif

end subroutine checkConsistantCoords


!------------------------------------------------------------------------------
! Compare two summation estimates of positive data and judge if due to more
! than round-off
!------------------------------------------------------------------------------
function isPosSumErrSignificant(n1, sum1, n2, sum2)
!------------------------------------------------------------------------------
! When two sums are calculated from different vectors that should add up to
! the same value, the results can differ by round off. The round off error
! can be bounded to be proportional to the number of operations.
! This function returns true if the difference between sum1 and sum2 is
! larger than than the estimated round off bound.
! NOTE: This estimate/function is only valid for summation of positive data.
!------------------------------------------------------------------------------
  ! Arguments
  integer, intent(in) :: n1, n2
  real,    intent(in) :: sum1, sum2
  logical             :: isPosSumErrSignificant
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


!------------------------------------------------------------------------------
! Sum the product of two arrays
!------------------------------------------------------------------------------
subroutine sumHtimesQ(nz, h, q, sumHQ, sumErr)
!------------------------------------------------------------------------------
! This routine calculates the sum of h(:)*q(:), and optionally returns a
! bound on the roundoff error in the sum.
!------------------------------------------------------------------------------

  ! Arguments
  integer,             intent(in)  :: nz
  real, dimension(nz), intent(in)  :: h, q
  real,                intent(out) :: sumHQ
  real, optional,      intent(out) :: sumErr

  integer :: k
  real :: hq, eps

  if (present(sumErr)) then
    ! Calculate the sum and estimate errors
    eps = epsilon(q(1))
    sumErr=0.
    sumHQ = 0.
    do k = 1,nz
      hq = h(k)*q(k)
      sumHQ = sumHQ + hq
      if (k>1) sumErr = sumErr + eps*max(abs(sumHQ),abs(hq))
    end do
  else
    ! Calculate the sum
    sumHQ = 0.
    do k = 1,nz
      sumHQ = sumHQ + h(k)*q(k)
    end do
  endif

end subroutine sumHtimesQ


!------------------------------------------------------------------------------
! Compare two summation estimates of signed data and judge if due to more
! than round-off
!------------------------------------------------------------------------------
function isSignedSumErrSignificant(n1, maxTerm1, sum1, n2, maxTerm2, sum2)
!------------------------------------------------------------------------------
! When two sums are calculated from different vectors that should add up to
! the same value, the results can differ by round off. The round off error
! can be bounded to be proportional to the number of operations.
! This function returns true if the difference between sum1 and sum2 is
! larger than than the estimated round off bound.
!------------------------------------------------------------------------------
  ! Arguments
  integer, intent(in) :: n1, n2
  real,    intent(in) :: maxTerm1, sum1, maxTerm2, sum2
  logical             :: isSignedSumErrSignificant
  ! Local variables
  real :: sumErr, allowedErr, eps

  sumErr = abs(sum1-sum2)
  eps = epsilon(sumErr)
  allowedErr = eps*0.5*( real(n1-1)*max(abs(maxTerm1),abs(sum1)) &
                       + real(n2-1)*max(abs(maxTerm2),abs(sum2)) )
  if (sumErr>allowedErr) then
    write(0,*) 'isSignedSumErrSignificant: maxTerm1,maxTerm2=',maxTerm1,maxTerm2
    write(0,*) 'isSignedSumErrSignificant: sum1,sum2=',sum1,sum2
    write(0,*) 'isSignedSumErrSignificant: eps=',eps
    write(0,*) 'isSignedSumErrSignificant: err,n*eps*maxTerm=',sumErr,allowedErr
    isSignedSumErrSignificant = .true.
  else
    isSignedSumErrSignificant = .false.
  endif
end function isSignedSumErrSignificant


!------------------------------------------------------------------------------
! Remapping core routine
!------------------------------------------------------------------------------
subroutine remapping_core( CS, n0, h0, u0, n1, dx, u1 )
!------------------------------------------------------------------------------
! This routine is basic in that it simply takes two grids and remaps the
! field known on the first grid onto the second grid, following the rules
! stored in the structure CS.
!------------------------------------------------------------------------------

  ! Arguments
  type(remapping_CS), intent(in)    :: CS
  integer,            intent(in)    :: n0 ! Number of cells on source grid
  real, dimension(:), intent(in)    :: h0 ! cell widths on source grid
  real, dimension(:), intent(in)    :: u0 ! cell averages on source grid
  integer,            intent(in)    :: n1 ! Number of cells on target grid
  real, dimension(:), intent(in)    :: dx ! Change in interface positions
  real, dimension(:), intent(out)   :: u1 ! cell averages on target grid

  ! Local variables
  integer :: iMethod
  real, dimension(CS%nk,2)           :: ppoly_r_E            !Edge value of polynomial
  real, dimension(CS%nk,2)           :: ppoly_r_S            !Edge slope of polynomial
  real, dimension(CS%nk,CS%degree+1) :: ppoly_r_coefficients !Coefficients of polynomial
  integer :: remapping_scheme

#ifdef __DO_SAFETY_CHECKS__
  integer :: k
  real :: hTmp, totalH0, totalHf, eps
  real :: err0, totalHU0, err2, totalHU2
  real :: z0, z1

  if (dx(1) /= 0.) call MOM_error( FATAL, 'MOM_remapping, remapping_core: '//&
             'Non-zero surface flux!' ) ! This is technically allowed but in avoided practice 
  totalH0 = 0.
  do k=1, n0
    totalH0 = totalH0 + h0(k)
  enddo
  totalHf = 0.
  do k=1, n1
    if (k <= n0) then
      hTmp = h0(k) + ( dx(k+1) - dx(k) )
      if (hTmp < 0.) then
        write(0,*) 'k,h0(k),hTmp,dx(k+1),dx(k)=',k,h0(k),hTmp,dx(k+1),dx(k)
        call MOM_error( FATAL, 'MOM_remapping, remapping_core: '//&
             'negative h implied by fluxes' )
      endif
    else
      hTmp = ( dx(k+1) - dx(k) )
      if (hTmp < 0.) then
        write(0,*) 'k,hNew,dx(+1),dx(0)=',k,dx(k+1),dx(k)
        call MOM_error( FATAL, 'MOM_remapping, remapping_core: '//&
             'negative h implied by fluxes' )
      endif
    endif
    totalHf = totalHf + hTmp
  end do
  eps = epsilon(hTmp)*totalH0
  if (abs(totalHf-totalH0) > 0.5*real(n0+n1-1)*eps) then
    write(0,*) 'H0,Hf=',totalH0,totalHf,totalHf-totalH0,eps
    call MOM_error( FATAL, 'MOM_remapping, remapping_core: '//&
         'Total thicknesses of h0 and h2 differ by more than roundoff' )
  endif
#endif

  iMethod = -999
  
  ! Reset polynomial
  ppoly_r_E(:,:) = 0.0
  ppoly_r_S(:,:) = 0.0
  ppoly_r_coefficients(:,:) = 0.0

  remapping_scheme = CS%remapping_scheme
  if (n0<=1) then
    remapping_scheme = REMAPPING_PCM
  elseif (n0<=3) then
    remapping_scheme = min( remapping_scheme, REMAPPING_PLM )
  elseif (n0<=4) then
    remapping_scheme = min( remapping_scheme, REMAPPING_PPM_H4 )
  endif
  select case ( remapping_scheme )
    case ( REMAPPING_PCM )
      call PCM_reconstruction( n0, u0, ppoly_r_E, ppoly_r_coefficients)
      iMethod = INTEGRATION_PCM
    case ( REMAPPING_PLM )
      call PLM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      if ( CS%boundary_extrapolation) then
        call PLM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients)
      end if    
      iMethod = INTEGRATION_PLM
    case ( REMAPPING_PPM_H4 )
      call edge_values_explicit_h4( n0, h0, u0, ppoly_r_E )
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      if ( CS%boundary_extrapolation) then
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      end if
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PPM_IH4 )
      call edge_values_implicit_h4( n0, h0, u0, ppoly_r_E )
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      if ( CS%boundary_extrapolation) then
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      end if    
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PQM_IH4IH3 )
      call edge_values_implicit_h4( n0, h0, u0, ppoly_r_E )
      call edge_slopes_implicit_h3( n0, h0, u0, ppoly_r_S )
      call PQM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      if ( CS%boundary_extrapolation) then
        call PQM_boundary_extrapolation_v1( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      end if    
      iMethod = INTEGRATION_PQM
    case ( REMAPPING_PQM_IH6IH5 )
      call edge_values_implicit_h6( n0, h0, u0, ppoly_r_E )
      call edge_slopes_implicit_h5( n0, h0, u0, ppoly_r_S )
      call PQM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      if ( CS%boundary_extrapolation) then
        call PQM_boundary_extrapolation_v1( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      end if    
      iMethod = INTEGRATION_PQM
    case default
      call MOM_error( FATAL, 'MOM_remapping, remapping_core: '//&
           'The selected remapping method is invalid' )
  end select
#ifdef __DO_SAFETY_CHECKS__
  do k = 1, n0
    if (ppoly_r_coefficients(k,2) /= ppoly_r_coefficients(k,2)) then
      write(0,*) 'NaN in PPOLY at k=',k,' extrap=',CS%boundary_extrapolation
      write(0,*) 'h0=',h0(k),' u0(k)',u0(k)
      write(0,*) 'ppoly_r_E=',ppoly_r_E(k,:)
      write(0,*) 'ppoly_r_coefficients(k)=',ppoly_r_coefficients(k,:)
    endif
  enddo
#endif

  call remapByDeltaZ( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients, n1, dx, iMethod, u1 )
! call remapByProjection( n0, h0, u0, CS%ppoly_r, n1, h1, iMethod, u1 )

#ifdef __DO_SAFETY_CHECKS__
  totalHU0 = 0.
  err0 = 0.
  do k = 1, n0
    hTmp = h0(k) * u0(k)
    totalHU0 = totalHU0 + hTmp
    err0 = err0 + epsilon(err0)*max(err0,abs(hTmp))
  enddo
  totalHU2 = 0. ; err2 = 0. ; z0 =0. ; z1 = 0.
  do k = 1, n1
    if (k <= n0) then
      hTmp = h0(k) + ( dx(k+1) - dx(k) )
    else
      hTmp = ( dx(k+1) - dx(k) )
    endif
    hTmp = hTmp * u1(k)
    totalHU2 = totalHU2 + hTmp
    err2 = err2 + epsilon(err2)*max(err2,abs(hTmp))
    if (u1(k) /= u1(k)) then
      write(0,*) 'NaN detected at k=',k
      write(0,*) 'h0(k)=',h0(k),' u0(k)=',u0(k)
      write(0,*) 'z0(k)=',z0,' z0(k+1)=',z0+h0(k)
      write(0,*) 'w(k)=',dx(k),' w(k+1)=',dx(k+1)
      write(0,*) 'z1(k)=',z1,' z1(k+1)=',z1+(h0(k)+(dx(k+1)-dx(k)))
      write(0,*) 'z0(k)+w(k)=',z0-dx(k),' z1(k+1)+w(k+1)=',(z0+h0(k))-dx(k+1)
      write(0,*) 'h1(k)=',h0(k) + ( dx(k+1) - dx(k) )
      hTmp = 0.;err0 = 0.
      do iMethod = 1, n1
        hTmp = hTmp + h0(iMethod)
        err0 = err0 + (h0(iMethod)+(dx(iMethod+1)-dx(iMethod)))
        write(0,*) iMethod,hTmp,h0(iMethod),u0(iMethod),err0,h0(iMethod)+(dx(iMethod+1)-dx(iMethod)),u1(iMethod),dx(iMethod)
      enddo
      write(0,*) dx(n1+1)
      call MOM_error( FATAL, 'MOM_remapping, remapping_core: '//&
         'NaN detected!' )
    endif
    if (k<=n0) then; hTmp = h0(k); else; hTmp = 0.; endif
    z0 = z0 + hTmp ; z1 = z1 + ( hTmp + ( dx(k+1) - dx(k) ) )
  enddo
  ! Maximum error based on guess at maximum roundoff
  if (abs(totalHU2-totalHU0) > (err0+err2)*max(real(n0), real(n1)) .and. (err0+err2)/=0.) then
    ! Maximum relative error
    if (abs(totalHU2-totalHU0) / totalHU2 > 1e-09) then
      ! Maximum absolute error
      if (abs(totalHU2-totalHU0) > 1e-18) then
        write(0,*) 'h0=',h0
        write(0,*) 'hf=',h0(1:n1)+dx(2:n1+1)-dx(1:n1)
        write(0,*) 'u0=',u0
        write(0,*) 'u1=',u1
        write(0,*) 'total HU0,HUf,f-0=',totalHU0,totalHU2,totalHU2-totalHU0
        write(0,*) 'err0,errF=',err0,err2
        call MOM_error( FATAL, 'MOM_remapping, remapping_core: '//&
             'Total stuff on h0 and hF differ by more than maximum errors' )
      endif
    endif
  endif
#endif

end subroutine remapping_core


! -----------------------------------------------------------------------------
! remapByProjection (integration of reconstructed profile)
! -----------------------------------------------------------------------------
subroutine remapByProjection( n0, h0, u0, ppoly0_E, ppoly0_coefficients, n1, h1, method, u1 )
  ! Arguments
  integer,       intent(in)    :: n0     ! number of cells in source grid
  real,          intent(in)    :: h0(:)  ! source grid widths (size n0)
  real,          intent(in)    :: u0(:)  ! source cell averages (size n0)
  real,          intent(in)    :: ppoly0_E(:,:)            !Edge value of polynomial
  real,          intent(in)    :: ppoly0_coefficients(:,:) !Coefficients of polynomial 
  integer,       intent(in)    :: n1     ! number of cells in target grid
  real,          intent(in)    :: h1(:)  ! target grid widths (size n1)
  integer,       intent(in)    :: method ! remapping scheme to use
  real,          intent(out)   :: u1(:)  ! target cell averages (size n1)
  
  ! Local variables
  integer       :: iTarget
  real          :: xL, xR       ! coordinates of target cell edges  

  ! Loop on cells in target grid (grid1). For each target cell, we need to find
  ! in which source cells the target cell edges lie. The associated indexes are 
  ! noted j0 and j1.
  xR = 0. ! Left boundary is at x=0
  do iTarget = 1,n1
    ! Determine the coordinates of the target cell edges
    xL = xR
    xR = xL + h1(iTarget)

    call integrateReconOnInterval( n0, h0, u0, ppoly0_E, ppoly0_coefficients, method, &
                                   xL, xR, h1(iTarget), u1(iTarget) )
    
  end do ! end iTarget loop on target grid cells

end subroutine remapByProjection


! -----------------------------------------------------------------------------
! Remap using change in interface positions
! -----------------------------------------------------------------------------
subroutine remapByDeltaZ( n0, h0, u0, ppoly0_E, ppoly0_coefficients, n1, dx1, method, u1, h1 )
! The new grid is defined relative to the original grid by change
!  dx1(:) = xNew(:) - xOld(:)
! and the remapping calculated so that
!  hNew(k) qNew(k) = hOld(k) qOld(k) + F(k+1) - F(k)
! where
!  F(k) = dx1(k) qAverage
! and where qAverage is the average qOld in the region zOld(k) to zNew(k).
  ! Arguments
  integer,        intent(in)  :: n0     ! number of cells in source grid
  real,           intent(in)  :: h0(:)  ! source grid widths (size n0)
  real,           intent(in)  :: u0(:)  ! source cell averages (size n0)
  real,           intent(in)  :: ppoly0_E(:,:)            !Edge value of polynomial
  real,           intent(in)  :: ppoly0_coefficients(:,:) !Coefficients of polynomial 
  integer,        intent(in)  :: n1     ! number of cells in target grid
  real,           intent(in)  :: dx1(:) ! target grid edge positions (size n1+1)
  integer                     :: method ! remapping scheme to use
  real,           intent(out) :: u1(:)  ! target cell averages (size n1)
  real, optional, intent(out) :: h1(:)  ! target grid widths (size n1)
  
  ! Local variables
  integer :: iTarget
  real    :: xL, xR    ! coordinates of target cell edges  
  real    :: xOld, hOld, uOld
  real    :: xNew, hNew
  real    :: uhNew, hFlux, uAve, fluxL, fluxR
integer :: k
#ifdef __DO_SAFETY_CHECKS__
  real    :: h0Total

  h0Total = 0.
  do iTarget = 1, n0
    h0Total = h0Total + h0(iTarget)
  enddo
#endif

  ! Loop on cells in target grid. For each cell, iTarget, the left flux is
  ! the right flux of the cell to the left, iTarget-1.
  ! The left flux is initialized by started at iTarget=0 to calculate the
  ! right flux which can take into account the target left boundary being
  ! in the interior of the source domain.
  fluxR = 0.
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
    else
      hOld = 0.       ! as if for layers>n0, they were vanished
      uOld = 1.E30    ! and the initial value should not matter
    endif
    xNew = xOld + dx1(iTarget+1)
    xL = min( xOld, xNew )
    xR = max( xOld, xNew )

#ifdef __DO_SAFETY_CHECKS__
    if (xL < 0.) then
      write(0,*) 'h0=',h0
      write(0,*) 'dx1=',dx1
      write(0,*) 'xOld,xNew,xL,xR,i=',xOld,xNew,xL,xR,iTarget
      call MOM_error(FATAL,'MOM_remapping, remapByDeltaZ: xL too negative')
    endif
    if (xR > h0Total) then
      write(0,*) 'h0=',h0
      write(0,*) 'dx1=',dx1
      write(0,*) 'xOld,xNew,xL,xR,i=',xOld,xNew,xL,xR,iTarget
      call MOM_error(FATAL,'MOM_remapping, remapByDeltaZ: xR too positive')
    endif
#endif

    ! hFlux is the positive width of the remapped volume
    hFlux = abs(dx1(iTarget+1))
    call integrateReconOnInterval( n0, h0, u0, ppoly0_E, ppoly0_coefficients, method, &
                                   xL, xR, hFlux, uAve )
    ! uAve is the average value of u, independent of sign of dx1
    fluxR = dx1(iTarget+1)*uAve ! Includes sign of dx1

#ifdef XXX__DO_SAFETY_CHECKS__
    ! Valid for CFL<1
    if (dx1(iTarget+1)<h0(iTarget+1) .and. dx1(iTarget+1)>-hOld) then
      if (uAve<min(u0(iTarget),u0(iTarget+1))) then
        write(0,*) 'u,u(k),u(k+1)=',uAve,u0(iTarget),u0(iTarget+1)
        call MOM_error(FATAL,'MOM_remapping, remapByDeltaZ: undershoot in U')
      endif
      if (uAve>max(u0(iTarget),u0(iTarget+1))) then
        write(0,*) 'u,u(k),u(k+1)=',uAve,u0(iTarget),u0(iTarget+1)
        call MOM_error(FATAL,'MOM_remapping, remapByDeltaZ: overshoot in U')
      endif
    endif
#endif

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


! -----------------------------------------------------------------------------
! integrate the reconstructed profile over a single cell
! -----------------------------------------------------------------------------
subroutine integrateReconOnInterval( n0, h0, u0, ppoly0_E, ppoly0_coefficients, method, &
                                     xL, xR, hC, uAve )
  ! Arguments
  integer, intent(in)  :: n0       ! number of cells in source grid
  real,    intent(in)  :: h0(:)    ! source grid sizes (size n0)
  real,    intent(in)  :: u0(:)       ! source cell averages
  real,    intent(in)  :: ppoly0_E(:,:)            !Edge value of polynomial
  real,    intent(in)  :: ppoly0_coefficients(:,:) !Coefficients of polynomial 
  integer, intent(in)  :: method   ! remapping scheme to use
  real,    intent(in)  :: xL, xR   ! left/right edges of target cell
  real,    intent(in)  :: hC       ! cell width hC = xR - xL
  real,    intent(out) :: uAve     ! average value on target cell
  
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

#ifdef __DO_SAFETY_CHECKS__
  real    :: h0Total

  h0Total = 0.
  do k = 1, n0
    h0Total = h0Total + h0(k)
  enddo
! if (xL < 0.) call MOM_error(FATAL, &
!         'MOM_remapping, integrateReconOnInterval: '//&
!         'The target cell starts beyond the left edge of the source grid')
! if (xR < 0.) call MOM_error(FATAL, &
!         'MOM_remapping, integrateReconOnInterval: '//&
!         'The target cell ends beyond the left edge of the source grid')
! if (xL > h0Total) call MOM_error(FATAL, &
!         'MOM_remapping, integrateReconOnInterval: '//&
!         'The target cell starts beyond the right edge of the source grid')
! if (xR > h0Total) call MOM_error(FATAL, &
!         'MOM_remapping, integrateReconOnInterval: '//&
!         'The target cell ends beyond the right edge of the source grid')
#endif

  q = -1.E30
  x0jLl = -1.E30
  x0jRl = -1.E30

  ! Find the left most cell in source grid spanned by the target cell
  jL = -1
  x0jLr = 0.
  do j = 1,n0
    x0jLl = x0jLr
    x0jLr = x0jLl + h0(j)
    ! Left edge is found in cell j
    if ( ( xL >= x0jLl ) .AND. ( xL <= x0jLr ) ) then
      jL = j
      exit ! once target grid cell is found, exit loop
    endif
  enddo

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
  if ( abs(xR - xL) == 0.0 ) then
    
    ! We check whether the source cell (i.e. the cell in which the
    ! vanished target cell lies) is vanished. If it is, the interpolated 
    ! value is set to be mean of the edge values (which should be the same).
    ! If it isn't, we simply interpolate.
    if ( h0(jL) == 0.0 ) then
      uAve = 0.5 * ( ppoly0_E(jL,1) + ppoly0_E(jL,2) )
    else
      ! WHY IS THIS NOT WRITTEN AS xi0 = ( xL - x0jLl ) / h0(jL) ---AJA
      xi0 = xL / h0(jL) - x0jLl / h0(jL)
  
      select case ( method )
        case ( INTEGRATION_PCM )   
          uAve = ppoly0_coefficients(jL,1)
        case ( INTEGRATION_PLM )  
          uAve = evaluation_polynomial( ppoly0_coefficients(jL,:), 2, xi0 )
        case ( INTEGRATION_PPM )
          uAve = evaluation_polynomial( ppoly0_coefficients(jL,:), 3, xi0 )
        case ( INTEGRATION_PQM )
          uAve = evaluation_polynomial( ppoly0_coefficients(jL,:), 5, xi0 )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select   
      
    end if ! end checking whether source cell is vanished
  
  ! 2. Cell is not vanished
  else

    ! Find the right most cell in source grid spanned by the target cell
    jR = -1
    x0jRr = 0.
    do j = 1,n0
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

#ifdef __DO_SAFETY_CHECKS__
    if ( jR == -1 ) call MOM_error(FATAL, &
          'MOM_remapping, integrateReconOnInterval: '//&
          'The location of the right-most cell could not be found')
#endif

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
      xi0 = max( 0., min( 1., ( xL - x0jLl ) / h0(jL) ) )
      xi1 = max( 0., min( 1., ( xR - x0jLl ) / h0(jL) ) )
#else
      xi0 = xL / h0(jL) - x0jLl / h0(jL)
      xi1 = xR / h0(jL) - x0jLl / h0(jL)
#endif

      hAct = h0(jL) * ( xi1 - xi0 )

      ! Depending on which polynomial is used, integrate quantity
      ! between xi0 and xi1. Integration is carried out in normalized
      ! coordinates, hence: \int_xL^xR p(x) dx = h \int_xi0^xi1 p(xi) dxi
      select case ( method )
        case ( INTEGRATION_PCM )     
          q = ppoly0_coefficients(jL,1) * ( xR - xL )
        case ( INTEGRATION_PLM )    
          q = h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0_coefficients(jL,:), 1 )
        case ( INTEGRATION_PPM )
          q = h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0_coefficients(jL,:), 2 )
        case ( INTEGRATION_PQM )
          q = h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0_coefficients(jL,:), 4 )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select     
#ifdef __DO_SAFETY_CHECKS__
      if (q /= q) then
        write(0,*) 'Nan at jL==jR: jL=',jL,' jR=',jR
        write(0,*) 'xL=',XL,' xR=',xR
        write(0,*) 'xi0=',xi0,' xi1=',xi1
        stop 'Nan during __DO_SAFETY_CHECKS__'
      endif
#endif

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
      xi0 = max( 0., min( 1., ( xL - x0jLl ) / h0(jL) ) )
#else
      xi0 = (xL - x0jLl) / h0(jL)
#endif
      xi1 = 1.0

      hAct = h0(jL) * ( xi1 - xi0 )

      select case ( method )
        case ( INTEGRATION_PCM )     
          q = q + ppoly0_coefficients(jL,1) * ( x0jLr - xL )
        case ( INTEGRATION_PLM )    
          q = q + h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0_coefficients(jL,:), 1 )
        case ( INTEGRATION_PPM )
          q = q + h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0_coefficients(jL,:), 2 )
        case ( INTEGRATION_PQM )
          q = q + h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0_coefficients(jL,:), 4 )
        case default
          call MOM_error( FATAL, 'The selected integration method is invalid' )
      end select     
#ifdef __DO_SAFETY_CHECKS__
      if (q /= q) then
        write(0,*) 'Nan on left segment: jL=',jL,' jR=',jR
        write(0,*) 'xL=',XL,' xR=',xR
        write(0,*) 'xi0=',xi0,' xi1=',xi1
        stop 'Nan during __DO_SAFETY_CHECKS__'
      endif
#endif
  
      ! Integrate contents within cells strictly comprised between jL and jR
      if ( jR > (jL+1) ) then
        do k = jL+1,jR-1
          q = q + h0(k) * u0(k)
          hAct = hAct + h0(k)
        end do
      end if
#ifdef __DO_SAFETY_CHECKS__
      if (q /= q) then
        write(0,*) 'Nan on middle segment: jL=',jL,' jR=',jR
        write(0,*) 'xL=',XL,' xR=',xR
        stop 'Nan during __DO_SAFETY_CHECKS__'
      endif
#endif

      ! Integrate from left boundary of cell jR up to xR
      xi0 = 0.0
#ifdef __USE_ROUNDOFF_SAFE_ADJUSTMENTS__
      xi1 = max( 0., min( 1., ( xR - x0jRl ) / h0(jR) ) )
#else
      xi1 = (xR - x0jRl) / h0(jR)
#endif

      hAct = hAct + h0(jR) * ( xi1 - xi0 )
    
      select case ( method )
        case ( INTEGRATION_PCM )     
          q = q + ppoly0_coefficients(jR,1) * ( xR - x0jRl )
        case ( INTEGRATION_PLM )    
          q = q + h0(jR) * &
              integration_polynomial( xi0, xi1, ppoly0_coefficients(jR,:), 1 )
        case ( INTEGRATION_PPM )
          q = q + h0(jR) * &
              integration_polynomial( xi0, xi1, ppoly0_coefficients(jR,:), 2 )
        case ( INTEGRATION_PQM )
          q = q + h0(jR) * &
              integration_polynomial( xi0, xi1, ppoly0_coefficients(jR,:), 4 )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select     
#ifdef __DO_SAFETY_CHECKS__
      if (q /= q) then
        write(0,*) 'Nan on right segment: jL=',jL,' jR=',jR
        write(0,*) 'h0(jR)=',h0(jR)
        write(0,*) 'xL=',xL,' xR=',xR
        write(0,*) 'xi0=',xi0,' xi1=',xi1
        write(0,*) 'ppoly0_coeff(jR)=',ppoly0_coefficients(jR,:)
        stop 'Nan during __DO_SAFETY_CHECKS__'
      endif
#endif

    end if ! end integration for non-vanished cells 

    ! The cell average is the integrated value divided by the cell width
#ifdef __USE_ROUNDOFF_SAFE_ADJUSTMENTS__
    uAve = q / hAct
#else
    uAve = q / hC
#endif
  
  end if ! end if clause to check if cell is vanished
    
end subroutine integrateReconOnInterval


!------------------------------------------------------------------------------
! dzFromH1H2
!------------------------------------------------------------------------------
subroutine dzFromH1H2( n1, h1, n2, h2, dx )
! ------------------------------------------------------------------------------
! Calculates the change in interface positions based on h1 and h2
! ------------------------------------------------------------------------------
  
  ! Arguments
  integer,            intent(in)  :: n1 ! Number of cells on source grid
  real, dimension(:), intent(in)  :: h1 ! cell widths of source grid (size n1)
  integer,            intent(in)  :: n2 ! Number of cells on target grid
  real, dimension(:), intent(in)  :: h2 ! cell widths of target grid (size n2)
  real, dimension(:), intent(out) :: dx ! Change in interface position (size n2+1)
    
  ! Local variables
  integer :: k
  real :: x1, x2

  x1 = 0.
  x2 = 0.
  dx(1) = 0.
  do k = 1, max(n1,n2)
    if (k <= n1) x1 = x1 + h1(k) ! Interface k+1, right of source cell k
    if (k <= n2) then
      x2 = x2 + h2(k) ! Interface k+1, right of target cell k
      dx(k+1) = x2 - x1 ! Change of interface k+1, target - source
    endif
  enddo
#ifdef __DO_SAFETY_CHECKS__
  if (abs(x2-x1) > 0.5*(real(n1-1)*x1+real(n2-1)*x2)*epsilon(x1)) then
    write(0,'(a4,3a12)') 'k','h1','h2','dh'
    do k = 1,max(n1,n2)
      if (k<=min(n1,n2)) then
        write(0,'(i4,3es12.3)') k,h1(k),h2(k),dx(k)
      elseif (k>n1) then
        write(0,'(i4,12x,2es12.3)') k,h2(k),dx(k)
      else
        write(0,'(i4,es12.3)') k,h1(k)
      endif
    enddo
    write(0,'(i4,24x,es12.3)') n2+1,dx(n2+1)
    write(0,*) 'x1,x2,x2-x1',x1,x2,x2-x1
    call MOM_error(FATAL,'MOM_remapping, dzFromH1H2: Bottom has moved!')
  endif
#endif

end subroutine dzFromH1H2


!------------------------------------------------------------------------------
! Constructor for remapping
!------------------------------------------------------------------------------
subroutine initialize_remapping( nk, remappingScheme, CS)
  ! Arguments
  integer, intent(in)                  :: nk
  character(len=*),      intent(in)    :: remappingScheme
  type(remapping_CS),    intent(inout) :: CS
  
  CS%nk = nk

  call setReconstructionType( remappingScheme, CS )

end subroutine initialize_remapping


!------------------------------------------------------------------------------
! Set the type of reconstruction
! Use this routine to parse a string parameter specifying the reconstruction 
! and re-allocates work arrays appropriately. It is called from
! initialize_remapping but can be called from an external module too.
!------------------------------------------------------------------------------
subroutine setReconstructionType(string,CS)
  ! Arguments
  character(len=*),   intent(in)    :: string
  type(remapping_CS), intent(inout) :: CS
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

!------------------------------------------------------------------------------
! Function to enable extrapolation in boundary cells
!------------------------------------------------------------------------------
subroutine remapEnableBoundaryExtrapolation(CS)
! Use this to enable extrapolation at boundaries
  type(remapping_CS), intent(inout) :: CS
  CS%boundary_extrapolation = .true.
end subroutine remapEnableBoundaryExtrapolation

!------------------------------------------------------------------------------
! Function to disable extrapolation in boundary cells
!------------------------------------------------------------------------------
subroutine remapDisableBoundaryExtrapolation(CS)
! Use this to disable extrapolation at boundaries
  type(remapping_CS), intent(inout) :: CS
  CS%boundary_extrapolation = .false.
end subroutine remapDisableBoundaryExtrapolation

!------------------------------------------------------------------------------
! Memory deallocation for remapping
!------------------------------------------------------------------------------
subroutine end_remapping(CS)
  ! Arguments
  type(remapping_CS), intent(inout) :: CS

  CS%degree = 0

end subroutine end_remapping

logical function remappingUnitTests()
  ! Should only be called from a single/root thread
  ! Returns True if a test fails, otherwise False
  integer, parameter :: n0 = 4, n1 = 3, n2 = 6
  real :: h0(n0), x0(n0+1), u0(n0)
  real :: h1(n1), x1(n1+1), u1(n1), hn1(n1), dx1(n1+1)
  real :: h2(n2), x2(n2+1), u2(n2), hn2(n2), dx2(n2+1)
  data u0 /3., 1., -1., -3./   ! Linear profile, 4 at surface to -4 at bottom
  data h0 /4*0.75/ ! 4 uniform layers with total depth of 3
  data h1 /3*1./   ! 3 uniform layers with total depth of 3
  data h2 /6*0.5/  ! 6 uniform layers with total depth of 3
  type(remapping_CS) :: CS 
  real, allocatable, dimension(:,:) :: ppoly0_E, ppoly0_S, ppoly0_coefficients
  integer :: i
  real :: err

  write(*,*) '===== MOM_remapping: remappingUnitTests =================='
  remappingUnitTests = .false. ! Normally return false

  call buildGridFromH(n0, h0, x0)
  do i=1,n0+1
    err=x0(i)-0.75*real(i-1)
    if (abs(err)>real(i-1)*epsilon(err)) remappingUnitTests = .true.
  enddo
  call buildGridFromH(n1, h1, x1)
  do i=1,n1+1
    err=x1(i)-real(i-1)
    if (abs(err)>real(i-1)*epsilon(err)) remappingUnitTests = .true.
  enddo

  call initialize_remapping(n0, 'PPM_H4', CS)
  write(*,*) 'h0 (test data)'
  call dumpGrid(n0,h0,x0,u0)

  call dzFromH1H2( n0, h0, n1, h1, dx1 )
  call remapping_core( CS, n0, h0, u0, n1, dx1, u1 )
  do i=1,n1
    err=u1(i)-8./3.*(0.5*real(1+n1)-real(i))
    if (abs(err)>epsilon(err)) remappingUnitTests = .true.
  enddo
  write(*,*) 'h1 (by projection)'
  call dumpGrid(n1,h1,x1,u1)

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
    err=u1(i)-8./3.*(0.5*real(1+n1)-real(i))
    if (abs(err)>2.*epsilon(err)) remappingUnitTests = .true.
  enddo

  u1(:) = 0.
  call remapByDeltaZ( n0, h0, u0, ppoly0_E, ppoly0_coefficients, &
                      n1, x1-x0(1:n1+1), &
                      INTEGRATION_PPM, u1, hn1 )
  write(*,*) 'h1 (by delta)'
  call dumpGrid(n1,h1,x1,u1)
  hn1=hn1-h1
  do i=1,n1
    err=u1(i)-8./3.*(0.5*real(1+n1)-real(i))
    if (abs(err)>2.*epsilon(err)) remappingUnitTests = .true.
  enddo

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
    err=u2(i)-8./6.*(0.5*real(1+n2)-real(i))
    if (abs(err)>2.*epsilon(err)) remappingUnitTests = .true.
  enddo

  deallocate(ppoly0_E, ppoly0_S, ppoly0_coefficients)

  write(*,*) '=========================================================='

  contains

  subroutine dumpGrid(n,h,x,u)
  integer, intent(in) :: n
  real, dimension(:), intent(in) :: h,x,u
  integer :: i
  write(*,'("i=",20i10)') (i,i=1,n+1)
  write(*,'("x=",20es10.2)') (x(i),i=1,n+1)
  write(*,'("i=",5x,20i10)') (i,i=1,n)
  write(*,'("h=",5x,20es10.2)') (h(i),i=1,n)
  write(*,'("u=",5x,20es10.2)') (u(i),i=1,n)
  end subroutine dumpGrid

end function remappingUnitTests

end module MOM_remapping
