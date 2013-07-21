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
use regrid_grid1d_class, only : grid1D_t, grid1Dconstruct, grid1Ddestroy
use regrid_ppoly_class, only : ppoly_t, ppoly_init, ppoly_destroy
use polynomial_functions, only : evaluation_polynomial, integration_polynomial
use regrid_edge_values, only : edgeValueArrays
use regrid_edge_values, only : edge_values_explicit_h4, edge_values_implicit_h4
use regrid_edge_values, only : edge_values_implicit_h4, edge_values_implicit_h6
use regrid_edge_values, only : triDiagEdgeWorkAllocate, triDiagEdgeWorkDeallocate
use regrid_edge_slopes, only : edgeSlopeArrays
use regrid_edge_slopes, only : edge_slopes_implicit_h3, edge_slopes_implicit_h5
use regrid_edge_slopes, only : triDiagSlopeWorkAllocate, triDiagSlopeWorkDeallocate
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
  ! Work arrays
  type(grid1D_t)                  :: grid_start ! starting grid
  type(grid1D_t)                  :: grid_final ! final grid
  type(ppoly_t)                   :: ppoly_r    ! reconstruction ppoly
  real, dimension(:), allocatable :: u_column   ! generic variable
  type(edgeValueArrays)           :: edgeValueWrk ! Work space for edge values
  type(edgeSlopeArrays)           :: edgeSlopeWrk ! Work space for edge slopes
  ! Parameters
  integer :: nk = 0                    ! Number of layers/levels in vertical
  integer :: remapping_scheme = -911   ! Determines which reconstruction to use
  integer :: degree                    ! Degree of polynomical reconstruction
  logical :: boundary_extrapolation = .true.  ! If true, extrapolate boundaries
end type

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public remapping_main, remapping_core
public initialize_remapping, end_remapping
public rempaEnableBoundaryExtrapolation, remapDisableBoundaryExtrapolation
public setReconstructionType

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
#define __DO_SAFTEY_CHECKS__

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! General remapping routine 
!------------------------------------------------------------------------------
subroutine remapping_main( CS, G, h, dzInterface, h_new, tv, u, v )
!------------------------------------------------------------------------------
! This routine takes care of remapping all variable between the old and the
! new grids. When velocity components need to be remapped, thicknesses at
! velocity points are taken to be arithmetic averages of tracer thicknesses.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(remapping_CS),                               intent(inout) :: CS
  type(ocean_grid_type),                            intent(in)    :: G
  real, dimension(NIMEM_,NJMEM_,NKMEM_),            intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_),     intent(in)    :: dzInterface
  real, dimension(NIMEM_,NJMEM_,NKMEM_),            intent(in)    :: h_new
  type(thermo_var_ptrs),                            intent(inout) :: tv       
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), optional, intent(inout) :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), optional, intent(inout) :: v
  
  ! Local variables
  integer               :: i, j, k
  integer               :: nz
  real                  :: val, new_val
  integer               :: problem

  nz = G%ke
#ifdef __DO_SAFTEY_CHECKS__
    if (nz>CS%grid_start%nb_cells) call MOM_error(FATAL,'nz>nk_start')
    if (nz>CS%grid_final%nb_cells) call MOM_error(FATAL,'nz>nk_final')
#endif

  ! Remap tracer
  do j = G%jsc,G%jec
    do i = G%isc,G%iec
    
      ! Build the start and final grids
      CS%grid_start%h(:) = h(i,j,:)
      CS%grid_final%h(:) = h_new(i,j,:)
      call buildConsistentGrids(nz, CS%grid_start%h, CS%grid_final%h, CS%grid_start%x, CS%grid_final%x)
      
      do k = 1,nz
        CS%grid_start%h(k) = CS%grid_start%x(k+1) - CS%grid_start%x(k)
        CS%grid_final%h(k) = CS%grid_final%x(k+1) - CS%grid_final%x(k)
      end do
      
      call remapping_core(CS, CS%grid_start, tv%S(i,j,:), CS%grid_final, CS%u_column)
      
      tv%S(i,j,:) = CS%u_column(:)
      
      call remapping_core(CS, CS%grid_start, tv%T(i,j,:), CS%grid_final, CS%u_column)
     
      tv%T(i,j,:) = CS%u_column(:)

    end do
  end do
  
  ! Remap u velocity component
  if ( present(u) ) then
  do j = G%jsc,G%jec
    do i = G%iscB,G%iecB
    
      ! Build the start and final grids
      CS%grid_start%h(:) = 0.5 * ( h(i,j,:) + h(i+1,j,:) )
      CS%grid_final%h(:) = 0.5 * ( h_new(i,j,:) + h_new(i+1,j,:) )
      call buildConsistentGrids(nz, CS%grid_start%h, CS%grid_final%h, CS%grid_start%x, CS%grid_final%x)
      
      do k = 1,nz
        CS%grid_start%h(k) = CS%grid_start%x(k+1) - CS%grid_start%x(k)
        CS%grid_final%h(k) = CS%grid_final%x(k+1) - CS%grid_final%x(k)
      end do
  
      call remapping_core(CS, CS%grid_start, u(i,j,:), CS%grid_final, CS%u_column)
     
      u(i,j,:) = CS%u_column(:)
      
    end do
  end do
  end if
  
  ! Remap v velocity component
  if ( present(v) ) then
  do j = G%jscB,G%jecB
    do i = G%isc,G%iec

      ! Build the start and final grids
      CS%grid_start%h(:) = 0.5 * ( h(i,j,:) + h(i,j+1,:) )
      CS%grid_final%h(:) = 0.5 * ( h_new(i,j,:) + h_new(i,j+1,:) )
      call buildConsistentGrids(nz, CS%grid_start%h, CS%grid_final%h, CS%grid_start%x, CS%grid_final%x)

      do k = 1,nz
        CS%grid_start%h(k) = CS%grid_start%x(k+1) - CS%grid_start%x(k)
        CS%grid_final%h(k) = CS%grid_final%x(k+1) - CS%grid_final%x(k)
      end do

      call remapping_core(CS, CS%grid_start, v(i,j,:), CS%grid_final, CS%u_column)
     
      v(i,j,:) = CS%u_column(:)
      
    end do
  end do
  end if

end subroutine remapping_main


!------------------------------------------------------------------------------
! Build a final grid 
!------------------------------------------------------------------------------
subroutine buildConsistentGrids(nz, hs, hf, xs, xf)
!------------------------------------------------------------------------------
! This routine calculates the coordinates xs and xf consistently from
! hs and hf so that the edges of the domain line up.
! If suM(hs) and sum(hf) differ significantly and error is generated.
!------------------------------------------------------------------------------

  ! Arguments
  integer,               intent(in)    :: nz
  real, dimension(nz),   intent(in)    :: hs, hf
  real, dimension(nz+1), intent(inout) :: xs, xf

  integer :: k
  real    :: sumH1, sumH2, nonDimPos

  ! Build start grid
  call buildGridFromH(nz, hs, xs)
  sumH1 = xs(nz+1)

  ! Initial guess at final grid
  call buildGridFromH(nz, hf, xf)
  sumH2 = xf(nz+1)

#ifdef __DO_SAFTEY_CHECKS__
    call checkGridConsistentcies(nz, xs, nz, xf, strict=.false.)
    if (abs(sumH1-sumH2)>0.5*real(nz)*epsilon(sumH2)*(sumH1+sumH2)) then
      write(0,*) 'Start/final/start-final grid'
      do k = 1,nz+1
        write(0,'(i4,3es12.3)') k,xs(k),xf(k),xs(k)-xf(k)
      enddo
      write(0,*) 'eps,H*eps',epsilon(sumH2),0.5*epsilon(sumH2)*(sumH1+sumH2)
      call MOM_error(FATAL,'MOM_remapping, buildConsistentGrids: '//&
                     'Final and start grids do not match.')
    endif
#endif

! call makeGridsConsistent(nz, xs, nz, hf, xf)

end subroutine buildConsistentGrids


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
subroutine checkGridConsistentcies(ns, xs, nf, xf, strict)
!------------------------------------------------------------------------------
! Checks that xs and xf are consistent to within roundoff.
! If strict=False, the end points of xs and xf are allowed to differ by
! numerical roundoff due to the nature of summation to obtain xs, xf.
! If strict=True, the edn points must be identical.
!------------------------------------------------------------------------------

  ! Arguments
  integer, intent(in) :: ns, nf
  real,    intent(in) :: xs(ns+1), xf(nf+1)
  logical, intent(in) :: strict

  ! Local variables
  integer :: k
  real    :: sumHs, sumHf

  sumHs = xs(ns+1)-xs(1)
  sumHf = xf(nf+1)-xf(1)

  if (strict) then
    if (sumHf /= sumHs) call &
        MOM_error(FATAL,'MOM_remapping, checkGridConsistentcies: '//&
                        'Total thickness of two grids are not exactly equal..')
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
      call MOM_error(FATAL,'MOM_remapping, checkGridConsistentcies: '//&
              'Total thickness of two grids do not match to within round-off.')
    endif
  endif

end subroutine checkGridConsistentcies


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
! NOTE: This estimate/function is only valid for summation of postive data.
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
! Check that data remapped between two grids are conserved
!------------------------------------------------------------------------------
subroutine checkGridConservation(ns, hs, xs, us, nf, hf, xf, uf)
!------------------------------------------------------------------------------
! Checks that the sum of hs*us and hf*uf match. Also checks that the
! analgous sums in terms of sx and xf are consistant.
!------------------------------------------------------------------------------

  ! Arguments
  integer, intent(in) :: ns, nf
  real,    intent(in) :: hs(ns), xs(ns+1), us(ns)
  real,    intent(in) :: hf(nf), xf(nf+1), uf(nf)

  ! Local variables
  integer :: k
  real    :: sumHUs, errHUs, sumXUs, errXUs
  real    :: sumHUf, errHUf, sumXUf, errXUf

  call sumHtimesQ(ns, hs, us, sumHUs, errHUs)
  call sumHtimesQ(ns, xs(2:ns+1)-xs(1:ns), us, sumXUs, errXUs)
  call sumHtimesQ(nf, hf, uf, sumHUf, errHUf)
  call sumHtimesQ(ns, xf(2:nf+1)-xf(1:nf), uf, sumXUf, errXUf)
  if (abs(sumHUs-sumXUs)>errHUs+errXUs) then
    write(0,'("ns=",i4)') ns
    do k = 1,ns+1
      write(0,'(i4,"xs=",es12.3)') k,xs(k)
      if (k<=ns) write(0,'(i4,"hs,us=",2es12.3)') k,hs(k),us(k)
    enddo
    write(0,'("sumHUs,sumXUs=",2es12.3)') sumHUs,sumXUs
    write(0,'("err,errHUs,errXUs=",3es12.3)') abs(sumHUs-sumXUs),errHUs,errXUs
    call MOM_error(FATAL,'MOM_remapping, checkGridConservation: '//&
       'Total amount of stuff on start grid differs by more than round-off.')
  endif
  if (abs(sumHUf-sumXUf)>errHUf+errXUf) then
    write(0,'("nf=",i4)') nf
    do k = 1,nf+1
      write(0,'(i4,"xf=",es12.3)') k,xf(k)
      if (k<=nf) write(0,'(i4,"hf,uf=",2es12.3)') k,hf(k),uf(k)
    enddo
    write(0,'("sumHUf,sumXUf=",2es12.3)') sumHUf,sumXUf
    write(0,'("err,errHUf,errXUf=",3es12.3)') abs(sumHUf-sumXUf),errHUf,errXUf
    call MOM_error(FATAL,'MOM_remapping, checkGridConservation: '//&
       'Total amount of stuff on final grid differs by more than round-off.')
  endif
#ifdef DISABLE_CONSEEVATION_CHECK_BECAUSE_IT_FAILS______
  if (abs(sumHUf-sumHUs)>errHUf+errHUs) then
    write(0,'("ns,nf=",2i4)') ns,nf
    do k = 1,max(ns,nf)+1
      if (k<=min(ns+1,nf+1)) then
        write(0,'(i4,"xs,xf=",2es12.3)') k,xs(k),xf(k)
      elseif (k>ns+1) then
        write(0,'(i4,"   xf=",12x,es12.3)') k,xf(k)
      else
        write(0,'(i4,"xs   =",es12.3)') k,xs(k)
      endif
      if (k<=min(ns,nf)) then
        write(0,'(i4,"hs,us,hf,uf=",4es12.3)') k,hs(k),us(k),hf(k),uf(k)
      elseif (k>ns .and. k<=nf) then
        write(0,'(i4,"      hf,uf=",24x,2es12.3)') k,hf(k),uf(k)
      elseif (k>nf .and. k<=ns) then
        write(0,'(i4,"hs,us      =",2es12.3)') k,hs(k),us(k)
      endif
    enddo
    write(0,'("sumHUf,sumHUs=",2es12.3)') sumHUf,sumHUs
    write(0,'("err,errHUf,errHUs=",3es12.3)') abs(sumHUf-sumHUs),errHUf,errHUs
    call MOM_error(FATAL,'MOM_remapping, checkGridConservation: '//&
       'Total amount of stuff on two grids differs by more than round-off.')
  endif
#endif

end subroutine checkGridConservation


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
  integer :: n
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
! Make a second grid consistent with the first
!------------------------------------------------------------------------------
subroutine makeGridsConsistent(ns, xs, nf, hf, xf)
!------------------------------------------------------------------------------
! Adjusts xf so that the end points exactly match those of xs.
! It is best to have called checkGridConsistentcies with strict=false
! to ensure that the grids are already close and only differ due to
! round-off.
!------------------------------------------------------------------------------

  ! Arguments
  integer, intent(in)    :: ns, nf
  real,    intent(in)    :: xs(ns+1), hf(nf)
  real,    intent(inout) :: xf(nf+1)

  ! Local variables
  integer :: n, k
  real    :: nonDimPos, sumHs, sumHf

#ifdef __DO_SAFTEY_CHECKS__
    if (xf(1) /= xs(1)) call &
         MOM_error(FATAL,'MOM_remapping, makeGridsConsistent: '//&
                         'Starting point of two grids do not match.')
#endif

  ! Adjust new grid so that end-points match those of the start grid.
  sumHs = xs(ns+1)
  sumHf = xf(nf+1)
  if (sumHf/=sumHs) then
    xf(nf+1) = sumHs
    do k = nf,1,-1
      nonDimPos = xf(k) / sumHf ! Position of xf interface within column
      ! When nonDimPos -> 1, adjust xf towards a bottom-up integration
      ! When nonDomPos -> 0, keep xf at the original top-down integration
      xf(k) = (1.-nonDimPos) * xf(k) + nonDimPos * (xf(k+1) - hf(k))
    end do
  endif

#ifdef __DO_SAFTEY_CHECKS__
    call checkGridConsistentcies(ns, xs, nf, xf, strict=.true.)
#endif

end subroutine makeGridsConsistent


!------------------------------------------------------------------------------
! Remapping core routine
!------------------------------------------------------------------------------
subroutine remapping_core( CS, grid0, u0, grid1, u1 )
!------------------------------------------------------------------------------
! This routine is basic in that it simply takes two grids and remaps the
! field known on the first grid onto the second grid, following the rules
! stored in the structure CS.
!------------------------------------------------------------------------------

  ! Arguments
  type(remapping_CS), intent(inout)   :: CS
  type(grid1D_t), intent(in)          :: grid0
  real, dimension(:), intent(in)      :: u0
  type(grid1D_t), intent(in)          :: grid1
  real, dimension(:), intent(inout)   :: u1

  ! Local variables
  integer :: n0, n1

  n0 = grid0%nb_cells
  n1 = grid1%nb_cells
  
  ! Reset polynomial
  CS%ppoly_r%E(:,:) = 0.0
  CS%ppoly_r%S(:,:) = 0.0
  CS%ppoly_r%coefficients(:,:) = 0.0

  select case ( CS%remapping_scheme )
    case ( REMAPPING_PCM )
      call PCM_reconstruction( grid0, u0, CS%ppoly_r )
      call remapping_integration( n0, grid0%h, grid0%x, u0, CS%ppoly_r, &
                                  n1, grid1%h, grid1%x, u1, &
                                  INTEGRATION_PCM )
    case ( REMAPPING_PLM )
      call PLM_reconstruction( grid0, u0, CS%ppoly_r )
      if ( CS%boundary_extrapolation) then
        call PLM_boundary_extrapolation( grid0, u0, CS%ppoly_r )
      end if    
      call remapping_integration( n0, grid0%h, grid0%x, u0, CS%ppoly_r, &
                                  n1, grid1%h, grid1%x, u1, &
                                  INTEGRATION_PLM )
    case ( REMAPPING_PPM_H4 )
      call edge_values_explicit_h4( grid0, u0, CS%ppoly_r%E )
      call PPM_reconstruction( grid0, u0, CS%ppoly_r )
      if ( CS%boundary_extrapolation) then
        call PPM_boundary_extrapolation( grid0, u0, CS%ppoly_r )
      end if    
      call remapping_integration( n0, grid0%h, grid0%x, u0, CS%ppoly_r, &
                                  n1, grid1%h, grid1%x, u1, &
                                  INTEGRATION_PPM )
    case ( REMAPPING_PPM_IH4 )
      call edge_values_implicit_h4( grid0, CS%edgeValueWrk, u0, CS%ppoly_r%E )
      call PPM_reconstruction( grid0, u0, CS%ppoly_r )
      if ( CS%boundary_extrapolation) then
        call PPM_boundary_extrapolation( grid0, u0, CS%ppoly_r )
      end if    
      call remapping_integration( n0, grid0%h, grid0%x, u0, CS%ppoly_r, &
                                  n1, grid1%h, grid1%x, u1, &
                                  INTEGRATION_PPM )
    case ( REMAPPING_PQM_IH4IH3 )
      call edge_values_implicit_h4( grid0, CS%edgeValueWrk, u0, CS%ppoly_r%E )
      call edge_slopes_implicit_h3( grid0, CS%edgeSlopeWrk, u0, CS%ppoly_r%S )
      call PQM_reconstruction( grid0, u0, CS%ppoly_r )
      if ( CS%boundary_extrapolation) then
        call PQM_boundary_extrapolation_v1( grid0, u0, CS%ppoly_r )
      end if    
      call remapping_integration( n0, grid0%h, grid0%x, u0, CS%ppoly_r, &
                                  n1, grid1%h, grid1%x, u1, &
                                  INTEGRATION_PQM )
    case ( REMAPPING_PQM_IH6IH5 )
      call edge_values_implicit_h6( grid0, CS%edgeValueWrk, u0, CS%ppoly_r%E )
      call edge_slopes_implicit_h5( grid0, CS%edgeSlopeWrk, u0, CS%ppoly_r%S )
      call PQM_reconstruction( grid0, u0, CS%ppoly_r )
      if ( CS%boundary_extrapolation) then
        call PQM_boundary_extrapolation_v1( grid0, u0, CS%ppoly_r )
      end if    
      call remapping_integration( n0, grid0%h, grid0%x, u0, CS%ppoly_r, &
                                  n1, grid1%h, grid1%x, u1, &
                                  INTEGRATION_PQM )
    case default
      call MOM_error( FATAL, 'The selected remapping method is invalid' )
  end select

#ifdef __DO_SAFTEY_CHECKS__
    call checkGridConservation(grid0%nb_cells, grid0%h, grid0%x, u0, &
                               grid1%nb_cells, grid1%h, grid1%x, u1)
#endif

end subroutine remapping_core


! -----------------------------------------------------------------------------
! remapping_integration (integration of reconstructed profile)
! -----------------------------------------------------------------------------
subroutine remapping_integration( n0, h0, x0, u0, ppoly0, n1, h1, x1, u1, method )
  ! Arguments
  integer,       intent(in)    :: n0       ! number of cells in source grid
  real,          intent(in)    :: h0(n0)   ! source grid widths
  real,          intent(in)    :: x0(n0+1) ! source grid edge positions
  real,          intent(in)    :: u0(n0)   ! source cell averages
  type(ppoly_t), intent(in)    :: ppoly0   ! source piecewise polynomial
  integer,       intent(in)    :: n1       ! number of cells in target grid
  real,          intent(in)    :: h1(n1)   ! target grid widths
  real,          intent(in)    :: x1(n1+1) ! target grid edge positions
  real,          intent(inout) :: u1(n1)   ! target cell averages
  integer                      :: method   ! remapping scheme to use
  
  ! Local variables
  integer       :: iTarget, j, k
  real          :: xL, xR       ! coordinates of target cell edges  

  ! Loop on cells in target grid (grid1). For each target cell, we need to find
  ! in which source cells the target cell edges lie. The associated indexes are 
  ! noted j0 and j1.
  do iTarget = 1,n1
    ! Determine the coordinates of the target cell edges
    xL = x1(iTarget)
    xR = x1(iTarget+1)

!   if (x1<=(1.+epsilon(x1))*x0(n0+1)) then
!     if (x1>x0(n0+1)) then   ! HACK ALERT !!!!!! -----AJA
!       x1 = min( x0(n0+1), x1) ! Bound target grid to be within source grid
!     endif
!   endif

    call integrateReconOnCell( n0, x0, h0, u0, ppoly0, method, &
                               xL, xR, h1(iTarget), u1(iTarget) )
    
  end do ! end iTarget loop on target grid cells

end subroutine remapping_integration


! -----------------------------------------------------------------------------
! integrate the reconstructed profile over a single cell
! -----------------------------------------------------------------------------
subroutine integrateReconOnCell( n0, x0, h0, u0, ppoly0, method, &
                                 xL, xR, hC, uAve )
  ! Arguments
  integer,            intent(in)    :: n0       ! number of cells in source grid
  real,               intent(in)    :: x0(n0+1) ! source grid edges
  real,               intent(in)    :: h0(n0)   ! source grid sizes
  real, dimension(:), intent(in)    :: u0       ! source cell averages
  type(ppoly_t),      intent(in)    :: ppoly0   ! source piecewise polynomial
  integer,            intent(in)    :: method   ! remapping scheme to use
  real,               intent(in)    :: xL, xR   ! left/right edges of target cell
  real,               intent(in)    :: hC       ! cell width hC = xR - xL
  real,               intent(inout) :: uAve     ! average value on target cell
  
  ! Local variables
  integer       :: j, k
  integer       :: jL, jR       ! indexes of source cells containing target 
                                ! cell edges
  real          :: q0, q1       ! partially integrated quantities in source 
                                ! cells j0 and j1
  real          :: q            ! complete integration
  real          :: a, b         ! interval of integration (global coordinates)
  real          :: xi0, xi1     ! interval of integration (local -- normalized 
                                ! -- coordinates)

#ifdef __DO_SAFTEY_CHECKS__
! if (xL < x0(1)) call MOM_error(FATAL, &
!         'MOM_remapping, integrateReconOnCell: '//&
!         'The target cell starts beyond the left edge of the source grid')
! if (xR < x0(1)) call MOM_error(FATAL, &
!         'MOM_remapping, integrateReconOnCell: '//&
!         'The target cell ends beyond the left edge of the source grid')
! if (xL > x0(n0+1)) call MOM_error(FATAL, &
!         'MOM_remapping, integrateReconOnCell: '//&
!         'The target cell starts beyond the right edge of the source grid')
! if (xR > x0(n0+1)) call MOM_error(FATAL, &
!         'MOM_remapping, integrateReconOnCell: '//&
!         'The target cell ends beyond the right edge of the source grid')
#endif

  ! Find the left most cell in source grid spanned by the target cell
  jL = -1
  do j = 1,n0
    ! Left edge is found in cell j
    if ( ( xL >= x0(j) ) .AND. ( xL <= x0(j+1) ) ) then
      jL = j
      exit ! once target grid cell is found, exit loop
    endif
  enddo

  ! If, at this point, jL is equal to -1, it means the vanished
  ! cell lies outside the source grid. In other words, it means that
  ! the source and target grids do not cover the same physical domain
  ! and there is something very wrong !
  if ( jL == -1 ) call MOM_error(FATAL, &
          'MOM_remapping, integrateReconOnCell: '//&
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
      uAve = 0.5 * ( ppoly0%E(jL,1) + ppoly0%E(jL,2) )
    else
      ! WHY IS THIS NOT WRITTEN AS xi0 = ( xL - x0(jL) ) / h0(jL) ---AJA
      xi0 = xL / h0(jL) - x0(jL) / h0(jL)
  
      select case ( method )
        case ( INTEGRATION_PCM )   
          uAve = ppoly0%coefficients(jL,1)
        case ( INTEGRATION_PLM )  
          uAve = evaluation_polynomial( ppoly0%coefficients(jL,:), 2, xi0 )
        case ( INTEGRATION_PPM )
          uAve = evaluation_polynomial( ppoly0%coefficients(jL,:), 3, xi0 )
        case ( INTEGRATION_PQM )
          uAve = evaluation_polynomial( ppoly0%coefficients(jL,:), 5, xi0 )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select   
      
    end if ! end checking whether source cell is vanished
  
  ! 2. Cell is not vanished
  else

    ! Find the right most cell in source grid spanned by the target cell
    jR = -1
    do j = 1,n0
      ! Right edge is found in cell j
      if ( ( xR >= x0(j) ) .AND. ( xR <= x0(j+1) ) ) then
        jR = j
        exit  ! once target grid cell is found, exit loop
      endif
    enddo ! end loop on source grid cells

    ! HACK to avoid roundoff problems  THIS NEEDS TO BE REMOVED ---AJA
    if (xR>x0(n0+1)) jR = n0

#ifdef __DO_SAFTEY_CHECKS__
    if ( jR == -1 ) call MOM_error(FATAL, &
          'MOM_remapping, integrateReconOnCell: '//&
          'The location of the right-most cell could not be found')
#endif

    ! To integrate, two cases must be considered: (1) the target cell is
    ! entirely comtained within a cell of the source grid and (2) the target
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
      ! WHY IS THIS NOT WRITTEN AS xi0 = ( xL - x0(jL) ) / h0(jL) ---AJA
      ! WHY IS THIS NOT WRITTEN AS xi0 = ( xR - x0(jL) ) / h0(jL) ---AJA
      xi0 = xL / h0(jL) - x0(jL) / h0(jL)
      xi1 = xR / h0(jL) - x0(jL) / h0(jL)

      ! Depending on which polynomial is used, integrate quantity
      ! between xi0 and xi1. Integration is carried out in normalized
      ! coordinates, hence: \int_xL^xR p(x) dx = h \int_xi0^xi1 p(xi) dxi
      select case ( method )
        case ( INTEGRATION_PCM )     
          q = ppoly0%coefficients(jL,1) * ( xR - xL )
        case ( INTEGRATION_PLM )    
          q = h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(jL,:), 1 )
        case ( INTEGRATION_PPM )
          q = h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(jL,:), 2 )
        case ( INTEGRATION_PQM )
          q = h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(jL,:), 4 )
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
      !xi0 = xL / h0(jL) - x0(jL) / h0(jL)
      xi0 = (xL - x0(jL)) / h0(jL)
      xi1 = 1.0
      select case ( method )
        case ( INTEGRATION_PCM )     
          q = q + ppoly0%coefficients(jL,1) * ( x0(jL+1) - xL )
        case ( INTEGRATION_PLM )    
          q = q + h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(jL,:), 1 )
        case ( INTEGRATION_PPM )
          q = q + h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(jL,:), 2 )
        case ( INTEGRATION_PQM )
          q = q + h0(jL) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(jL,:), 4 )
        case default
          call MOM_error( FATAL, 'The selected integration method is invalid' )
      end select     
  
      ! Integrate contents within cells strictly comprised between jL and jR
      if ( jR > (jL+1) ) then
        do k = jL+1,jR-1
          q = q + h0(k) * u0(k)
        end do
      end if

      ! Integrate from left boundary of cell jR up to xR
      xi0 = 0.0
      !xi1 = xR / h0(jR) - x0(jR) / h0(jR)
      xi1 = (xR - x0(jR)) / h0(jR)
    
      select case ( method )
        case ( INTEGRATION_PCM )     
          q = q + ppoly0%coefficients(jR,1) * ( xR - x0(jR) )
        case ( INTEGRATION_PLM )    
          q = q + h0(jR) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(jR,:), 1 )
        case ( INTEGRATION_PPM )
          q = q + h0(jR) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(jR,:), 2 )
        case ( INTEGRATION_PQM )
          q = q + h0(jR) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(jR,:), 4 )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select     
      
    end if ! end integration for non-vanished cells 
    
    ! The cell average is the integrated value divided by the cell width
    uAve = q / hC
  
  end if ! end if clause to check if cell is vanished
    
end subroutine integrateReconOnCell


!------------------------------------------------------------------------------
! Constructor for remapping
!------------------------------------------------------------------------------
subroutine initialize_remapping( nk, remappingScheme, CS)
  ! Arguments
  integer, intent(in)                  :: nk
  character(len=40),     intent(in)    :: remappingScheme
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
  character(len=*), intent(in) :: string
  type(remapping_CS), intent(inout) :: CS
  ! Local variables
  integer :: degree
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

  if (allocated(CS%u_column) .and. degree/=CS%degree) then
    ! If the degree has changed then deallocate to force a re-allocation
    call end_remapping(CS)
  endif
  CS%degree = degree
  if (.not. allocated(CS%u_column)) then
    call allocate_remapping( CS )
  endif
  
end subroutine setReconstructionType

!------------------------------------------------------------------------------
! Functino to enable extraplation in boundary cells
!------------------------------------------------------------------------------
subroutine rempaEnableBoundaryExtrapolation(CS)
! Use this to enable extrapolation at boundaries
  type(remapping_CS), intent(inout) :: CS
  CS%boundary_extrapolation = .true.
end subroutine rempaEnableBoundaryExtrapolation

!------------------------------------------------------------------------------
! Functino to disable extraplation in boundary cells
!------------------------------------------------------------------------------
subroutine remapDisableBoundaryExtrapolation(CS)
! Use this to disable extrapolation at boundaries
  type(remapping_CS), intent(inout) :: CS
  CS%boundary_extrapolation = .false.
end subroutine remapDisableBoundaryExtrapolation

!------------------------------------------------------------------------------
! Memory allocation for remapping
!------------------------------------------------------------------------------
subroutine allocate_remapping( CS )
  ! Arguments
  type(remapping_CS),    intent(inout) :: CS
  
  call grid1Dconstruct( CS%grid_start, CS%nk )
  call grid1Dconstruct( CS%grid_final, CS%nk )
  call ppoly_init( CS%ppoly_r, CS%nk, CS%degree )
  allocate( CS%u_column(CS%nk) ); CS%u_column = 0.0
  call triDiagEdgeWorkAllocate( CS%nk, CS%edgeValueWrk )
  call triDiagSlopeWorkAllocate( CS%nk, CS%edgeSlopeWrk )

end subroutine allocate_remapping


!------------------------------------------------------------------------------
! Memory deallocation for remapping
!------------------------------------------------------------------------------
subroutine end_remapping(CS)
  ! Arguments
  type(remapping_CS), intent(inout) :: CS

  ! Deallocate memory for grid
  call grid1Ddestroy( CS%grid_start )
  call grid1Ddestroy( CS%grid_final )
  call ppoly_destroy( CS%ppoly_r )
  deallocate( CS%u_column )
  call triDiagEdgeWorkDeallocate( CS%edgeValueWrk )
  call triDiagSlopeWorkDeallocate( CS%edgeSlopeWrk )

end subroutine end_remapping

end module MOM_remapping
