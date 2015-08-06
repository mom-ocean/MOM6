module MOM_regridding
!==============================================================================
!
! This file is part of MOM.
!
! Date of creation: 2008.06.09
! L. White
!
! This module contains the main regridding routines. 
!
! Regridding must be understood as comprising two steps:
! (1) Interpolation and creation of a new grid based on target interface 
!     densities (or any other criterion).
! (2) Remapping of quantities between old grid and new grid.
!
!==============================================================================
use MOM_error_handler, only : MOM_error, FATAL
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use MOM_EOS,           only : calculate_density
use MOM_string_functions,only : uppercase

use regrid_edge_values, only : edge_values_explicit_h2, edge_values_explicit_h4
use regrid_edge_values, only : edge_values_implicit_h4, edge_values_implicit_h6
use regrid_edge_slopes, only : edge_slopes_implicit_h3, edge_slopes_implicit_h5

use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation
use PQM_functions, only : PQM_reconstruction, PQM_boundary_extrapolation_v1

use P1M_functions, only : P1M_interpolation, P1M_boundary_extrapolation
use P3M_functions, only : P3M_interpolation, P3M_boundary_extrapolation
use MOM_remapping, only : remapping_core
use MOM_remapping, only : remapping_CS
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA
use regrid_consts, only : REGRIDDING_ARBITRARY

implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Private (module-wise) variables
! -----------------------------------------------------------------------------

type, public :: regridding_CS
  private
                                                    
  ! This array is set by function setCoordinateResolution
  ! It contains the "resolution" or delta coordinate of the target
  ! coorindate. It has the units of the target coordiante, e.g.
  ! meters for z*, non-dimensional for sigma, etc.
  real, dimension(:), allocatable :: coordinateResolution
  ! This array is set by function setcoordinateInterfaces
  ! This array is nominal coordinate of interfaces and is the
  ! running sum of coordinateResolution. i.e.
  !  coordinateInterfaces(k) = coordinateResolution(k+1)-coordinateResolution(k)
  ! It is only used in "rho" mode.
  real, dimension(:), allocatable :: coordinateInterfaces

  integer   :: nk ! Number of layers/levels

  integer   :: degree_i=4 !Degree of interpolation polynomial

  ! Indicates which grid to use in the vertical (z*, sigma, target interface
  ! densities)
  integer   :: regridding_scheme

  ! The following parameter is only relevant when used with the target
  ! interface densities regridding scheme. It indicates which interpolation
  ! to use to determine the grid.
  integer   :: interpolation_scheme

  ! Indicate whether high-order boundary extrapolation should be used within
  ! boundary cells
  logical   :: boundary_extrapolation

  ! Minimum thickness allowed when building the new grid through regridding
  real      :: min_thickness

end type

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public initialize_regridding
public end_regridding
public regridding_main 
public check_grid_integrity
public setRegriddingBoundaryExtrapolation
public setRegriddingMinimumThickness
public uniformResolution
public setCoordinateResolution
public setCoordinateInterfaces
public getCoordinateResolution
public getCoordinateInterfaces
public getCoordinateUnits
public getCoordinateShortName
public getStaticThickness
public buildGridZstarColumn

public DEFAULT_COORDINATE_MODE
character(len=158), parameter, public :: regriddingCoordinateModeDoc = &
                 " LAYER - Isopycnal or stacked shallow water layers\n"//&
                 " Z*    - stetched geopotential z*\n"//&
                 " SIGMA - terrain following coordinates\n"//&
                 " RHO   - continuous isopycnal\n"
character(len=338), parameter, public :: regriddingInterpSchemeDoc = &
                 " P1M_H2     (2nd-order accurate)\n"//&
                 " P1M_H4     (2nd-order accurate)\n"//&
                 " P1M_IH4    (2nd-order accurate)\n"//&
                 " PLM        (2nd-order accurate)\n"//&
                 " PPM_H4     (3rd-order accurate)\n"//&
                 " PPM_IH4    (3rd-order accurate)\n"//&
                 " P3M_IH4IH3 (4th-order accurate)\n"//&
                 " P3M_IH6IH5 (4th-order accurate)\n"//&
                 " PQM_IH4IH3 (4th-order accurate)\n"//&
                 " PQM_IH6IH5 (5th-order accurate)"
character(len=6), parameter, public :: regriddingDefaultInterpScheme = "P1M_H2"
logical, parameter, public :: regriddingDefaultBoundaryExtrapolation = .false.
real, parameter, public :: regriddingDefaultMinThickness = 1.e-3

! -----------------------------------------------------------------------------
! The following are private constants
! -----------------------------------------------------------------------------
! List of interpolation schemes
integer, parameter :: INTERPOLATION_P1M_H2     = 0 ! O(h^2)
integer, parameter :: INTERPOLATION_P1M_H4     = 1 ! O(h^2)
integer, parameter :: INTERPOLATION_P1M_IH4    = 2 ! O(h^2)
integer, parameter :: INTERPOLATION_PLM        = 3 ! O(h^2)
integer, parameter :: INTERPOLATION_PPM_H4     = 4 ! O(h^3)
integer, parameter :: INTERPOLATION_PPM_IH4    = 5 ! O(h^3)
integer, parameter :: INTERPOLATION_P3M_IH4IH3 = 6 ! O(h^4)
integer, parameter :: INTERPOLATION_P3M_IH6IH5 = 7 ! O(h^4)
integer, parameter :: INTERPOLATION_PQM_IH4IH3 = 8 ! O(h^4)
integer, parameter :: INTERPOLATION_PQM_IH6IH5 = 9 ! O(h^5)

! List of interpolant degrees
integer, parameter :: DEGREE_1 = 1, DEGREE_2 = 2, DEGREE_3 = 3, DEGREE_4 = 4

! Maximum number of regridding iterations
integer, parameter :: NB_REGRIDDING_ITERATIONS = 1
! Deviation tolerance between succesive grids in regridding iterations
real, parameter    :: DEVIATION_TOLERANCE = 1e-10
! Maximum number of Newton-Raphson iterations. Newton-Raphson iterations are
! used to build the new grid by finding the coordinates associated with
! target densities and interpolations of degree larger than 1.
integer, parameter :: NR_ITERATIONS = 8
! Tolerance for Newton-Raphson iterations (stop when increment falls below this)
real, parameter    :: NR_TOLERANCE = 1e-12
! When the N-R algorithm produces an estimate that lies outside [0,1], the
! estimate is set to be equal to the boundary location, 0 or 1, plus or minus
! an offset, respectively, when the derivative is zero at the boundary.
real, parameter    :: NR_OFFSET = 1e-6

! This CPP macro embeds some safety checks
#undef __DO_SAFETY_CHECKS__

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! Initialization of regridding options
!------------------------------------------------------------------------------
subroutine initialize_regridding( nk, coordMode, interpScheme, CS )
!------------------------------------------------------------------------------
! This routine is typically called (from initialize_MOM in file MOM.F90)
! before the main time integration loop to initialize the regridding stuff.
! We read the MOM_input file to register the values of different
! regridding  parameters.
!------------------------------------------------------------------------------
  
  ! Arguments
  integer,               intent(in)    :: nk
  character(len=*),      intent(in)    :: coordMode
  character(len=*),      intent(in)    :: interpScheme
  type(regridding_CS),   intent(inout) :: CS

  ! Local variables

  CS%nk = nk

  CS%regridding_scheme = coordinateMode(coordMode)

  select case ( uppercase(trim(interpScheme)) )
    case ("P1M_H2");     CS%interpolation_scheme = INTERPOLATION_P1M_H2
    case ("P1M_H4");     CS%interpolation_scheme = INTERPOLATION_P1M_H4
    case ("P1M_IH2");    CS%interpolation_scheme = INTERPOLATION_P1M_IH4
    case ("PLM");        CS%interpolation_scheme = INTERPOLATION_PLM
    case ("PPM_H4");     CS%interpolation_scheme = INTERPOLATION_PPM_H4
    case ("PPM_IH4");    CS%interpolation_scheme = INTERPOLATION_PPM_IH4
    case ("P3M_IH4IH3"); CS%interpolation_scheme = INTERPOLATION_P3M_IH4IH3
    case ("P3M_IH6IH5"); CS%interpolation_scheme = INTERPOLATION_P3M_IH6IH5
    case ("PQM_IH4IH3"); CS%interpolation_scheme = INTERPOLATION_PQM_IH4IH3
    case ("PQM_IH6IH5"); CS%interpolation_scheme = INTERPOLATION_PQM_IH6IH5
    case default ; call MOM_error(FATAL, "read_regridding_options: "//&
     "Unrecognized choice for INTERPOLATION_SCHEME ("//trim(interpScheme)//").")
  end select

  CS%boundary_extrapolation = regriddingDefaultBoundaryExtrapolation

  CS%min_thickness = regriddingDefaultMinThickness

  call allocate_regridding( CS )
  
end subroutine initialize_regridding



!------------------------------------------------------------------------------
! End of regridding (memory deallocation)
!------------------------------------------------------------------------------
subroutine end_regridding(CS)
!------------------------------------------------------------------------------
! This routine is typically called (from MOM_end in file MOM.F90)
! after the main time integration loop to deallocate the regridding stuff.
!------------------------------------------------------------------------------
  type(regridding_CS), intent(inout) :: CS
  
  ! Deallocate memory used for the regridding
  call regridding_memory_deallocation( CS )

end subroutine end_regridding


!------------------------------------------------------------------------------
! Dispatching regridding routine: regridding & remapping
!------------------------------------------------------------------------------
subroutine regridding_main( remapCS, CS, G, h, tv, dzInterface )
!------------------------------------------------------------------------------
! This routine takes care of (1) building a new grid and (2) remapping between
! the old grid and the new grid. The creation of the new grid can be based
! on z coordinates, target interface densities, sigma coordinates or any
! arbitrary coordinate system.
!   The MOM6 interface positions are always calculated from the bottom up by
! accumulating the layer thicknessea starting at z=-G%bathyT.  z increases
! upwards (decreasing k-index).
!   The new grid is defined by the change in position of those interfaces in z
!       dzInterface = zNew - zOld.
!   Thus, if the regridding inflates the top layer, hNew(1) > hOld(1), then the
! second interface moves downward, zNew(2) < zOld(2), and dzInterface(2) < 0.
!       hNew(k) = hOld(k) - dzInterface(k+1) + dzInterface(k)
! IMPORTANT NOTE:
!   This is the converse of the sign convention used in the remapping code!
!------------------------------------------------------------------------------
  
  ! Arguments
  type(remapping_CS),                      intent(in)    :: remapCS ! Remapping parameters and options
  type(regridding_CS),                     intent(in)    :: CS     ! Regridding parameters and options
  type(ocean_grid_type),                   intent(in)    :: G      ! Ocean grid informations
  real, dimension(NIMEM_,NJMEM_, NKMEM_),  intent(inout) :: h      ! Current 3D grid obtained after the last time step
  type(thermo_var_ptrs),                   intent(inout) :: tv     ! Thermodynamical variables (T, S, ...)  
  real, dimension(NIMEM_,NJMEM_, NK_INTERFACE_), intent(inout) :: dzInterface  ! The change in position of each interface

  real :: trickGnuCompiler

  ! Local variables
  
  select case ( CS%regridding_scheme )

    case ( REGRIDDING_ZSTAR )
      call buildGridZstar( CS, G, h, dzInterface )

    case ( REGRIDDING_SIGMA )
      call buildGridSigma( CS, G, h, dzInterface )
    
    case ( REGRIDDING_RHO )  
      call convective_adjustment(G, h, tv)
      call buildGridRho( G, h, tv, dzInterface, remapCS, CS )

    case ( REGRIDDING_ARBITRARY )
      call build_grid_arbitrary( G, h, dzInterface, trickGnuCompiler, CS )

    case default
      call MOM_error(FATAL,'MOM_regridding, regridding_main: '//&
                     'Unknown regridding scheme selected!')
      
  end select ! type of grid 
  
#ifdef __DO_SAFETY_CHECKS__
  call checkGridsMatch(G, h, dzInterface)
#endif

end subroutine regridding_main


!------------------------------------------------------------------------------
! Check that the total thickness of two grids match
!------------------------------------------------------------------------------
subroutine checkGridsMatch( G, h, dzInterface )
!------------------------------------------------------------------------------
! This routine calculates the total thickness of 
!------------------------------------------------------------------------------
  
  ! Arguments
  type(ocean_grid_type),                        intent(in) :: G
  real, dimension(NIMEM_,NJMEM_,NKMEM_),        intent(in) :: h
  real, dimension(NIMEM_,NJMEM_,NK_INTERFACE_), intent(in) :: dzInterface
  
  ! Local variables
  integer :: i, j, k
  integer :: nz
  real    :: totalHold, totalHnewF, eps, hNewF, zOld, zNewF

  nz = G%ke
  eps =1. ; eps = epsilon(eps)
!$OMP parallel do default(none) shared(G,nz,h,dzInterface,eps) &
!$OMP                          private(totalHold,zOld,totalHnewF,zNewF,hNewF)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      ! Total thickness of grid h
      totalHold = 0.
      do k = 1,nz
        totalHold = totalHold + h(i,j,k)
      enddo

      ! Integrate upwards for the interfaces consistent with the rest of MOM6
      zOld = - G%bathyT(i,j)
      totalHnewF = 0.
      do k = nz,1,-1
        zOld = zOld + h(i,j,k) ! Old interface position
        zNewF = zOld + dzInterface(i,j,k) ! New interface position based on dzInterface
       !hNewF = ( h(i,j,k) - dzInterface(i,j,k+1) ) + dzInterface(i,j,k)
        hNewF = h(i,j,k) + ( dzInterface(i,j,k) - dzInterface(i,j,k+1) )
        if (hNewF<0.) then
          write(0,*) 'k,h,hnew=',k,h(i,j,k),hNewF
          write(0,*) 'dzI(k+1),dzI(k)=',dzInterface(i,j,k+1),dzInterface(i,j,k)
          call MOM_error( FATAL, 'MOM_regridding, checkGridsMatch: '//&
            'Flux form led to negative layer thickness')
        endif
        totalHnewF = totalHnewF + hNewF

      enddo

      ! Conservation by implied hNewF
      if (abs(totalHnewF-totalHold)>real(nz-1)*0.5*(totalHold+totalHnewF)*eps) then
        write(0,*) 'i,j,nz=',i,j,nz
        do k = 1,nz
          write(0,*) 'k,h,hnew=',k,h(i,j,k),h(i,j,k)+(dzInterface(i,j,k)-dzInterface(i,j,k+1))
        enddo
        write(0,*) 'Hold,Hnew,Hnew-Hold=',totalHold,totalHnewF,totalHnewF-totalHold
        write(0,*) 'eps,(n)/2*eps*H=',eps,real(nz-1)*0.5*(totalHold+totalHnewF)*eps
        call MOM_error( FATAL, 'MOM_regridding, checkGridsMatch: '//&
          'Flux form did NOT conserve total thickness to within roundoff')
      endif

      if (dzInterface(i,j,1) /= 0.) call MOM_error( FATAL, &
          'MOM_regridding, checkGridsMatch: '//&
          'Non-zero dzInterface at surface!')
      if (dzInterface(i,j,nz+1) /= 0.) call MOM_error( FATAL, &
          'MOM_regridding, checkGridsMatch: '//&
          'Non-zero dzInterface at bottom!')
    enddo
  enddo

end subroutine checkGridsMatch

!------------------------------------------------------------------------------
! Build uniform z*-ccordinate grid with partial steps
!------------------------------------------------------------------------------
subroutine buildGridZstar( CS, G, h, dzInterface )
!------------------------------------------------------------------------------
! This routine builds a grid where the distribution of levels is based on a
! z* coordinate system with partial steps (Adcroft and Campin, 2004).
! The module parameter coordinateResolution(:) determines the nominal
! resolution, dz*(:), in the absence of topography. z* is defined by
!   z* = (z-eta)/(H+eta)*H  s.t. z*=0 when z=eta and z*=-H when z=-H .
!------------------------------------------------------------------------------
  
  ! Arguments
  type(regridding_CS),                           intent(in)    :: CS
  type(ocean_grid_type),                         intent(in)    :: G
  real, dimension(NIMEM_,NJMEM_,NKMEM_),         intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_, NK_INTERFACE_), intent(inout) :: dzInterface
  
  ! Local variables
  integer :: i, j, k
  integer :: nz
  real    :: nominalDepth, totalThickness, eta, stretching, dh
  real, dimension(SZK_(G)+1) :: zOld, zNew
  real :: minThickness

  nz = G%ke
  
!$OMP parallel do default(none) shared(G,dzInterface,CS,nz,h)                    &
!$OMP                          private(nominalDepth,totalThickness,minThickness, &
!$OMP                                  eta,stretching,zNew,dh,zOld)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      if (G%mask2dT(i,j)==0.) then
        dzInterface(i,j,:) = 0.
        cycle
      endif

      ! Local depth (G%bathyT is positive)
      nominalDepth = G%bathyT(i,j)*G%m_to_H

      ! Determine water column thickness
      totalThickness = 0.0
      do k = 1,nz
        totalThickness = totalThickness + h(i,j,k)
      end do

      call buildGridZStarColumn(CS, nz, nominalDepth, totalThickness, zNew)

      zOld(nz+1) = - nominalDepth
      do k = nz,1,-1
        zOld(k) = zOld(k+1) + h(i,j,k)
      enddo

      ! Define regridding in terms of a movement of interfaces
      dzInterface(i,j,1) = 0.
      do k = 2,nz
        dzInterface(i,j,k) = zNew(k) - zOld(k)
      enddo
      dzInterface(i,j,nz+1) = 0.

#ifdef __DO_SAFETY_CHECKS__
      dh=max(nominalDepth,totalThickness)
      if (abs(zNew(1)-zOld(1))>(nz-1)*0.5*epsilon(eta)*dh) then
        write(0,*) 'min_thickness=',minThickness
        write(0,*) 'eta=',eta,'nominalDepth=',nominalDepth,'totalThickness=',totalThickness
        write(0,*) 'dzInterface(1) = ',dzInterface(i,j,1),epsilon(eta),nz
        do k=1,nz+1
          write(0,*) k,zOld(k),zNew(k)
        enddo
        write(0,*) 'stretching=',stretching
        do k=1,nz
          write(0,*) k,h(i,j,k),zNew(k)-zNew(k+1),stretching * CS%coordinateResolution(k),CS%coordinateResolution(k)
        enddo
        call MOM_error( FATAL, &
               'MOM_regridding, buildGridZstar: top surface has moved!!!' )
      endif
#endif

    end do
  end do

end subroutine buildGridZstar

subroutine buildGridZstarColumn( CS, nz, depth, totalThickness, zInterface)

  ! Arguments
  type(regridding_CS), intent(in)    :: CS
  integer, intent(in) :: nz
  real, intent(in) :: depth
  real, intent(in) :: totalThickness
  real, dimension(nz+1), intent(inout) :: zInterface

  real :: eta, stretching, dh
  real :: minThickness
  integer :: k

  minThickness = min( CS%min_thickness, totalThickness/float(nz) )

  ! Position of free-surface
  eta = totalThickness - depth

  ! z* = (z-eta) / stretching   where stretching = (H+eta)/H
  ! z = eta + stretching * z*
  stretching = totalThickness / depth

  ! Integrate down from the top for a notional new grid, ignoring topography
  zInterface(1) = eta
  do k = 1,nz
    dh = stretching * CS%coordinateResolution(k) ! Notional grid spacing
    zInterface(k+1) = zInterface(k) - dh
  enddo

  ! Integrating up from the bottom adjusting interface position to accomodate
  ! inflating layers without disturbing the interface above
  zInterface(nz+1) = -depth
  do k = nz,1,-1
    if ( zInterface(k) < (zInterface(k+1) + minThickness) ) then
      zInterface(k) = zInterface(k+1) + minThickness
    endif
  enddo

end subroutine buildGridZstarColumn


!------------------------------------------------------------------------------
! Build sigma grid
!------------------------------------------------------------------------------
subroutine buildGridSigma( CS, G, h, dzInterface )
!------------------------------------------------------------------------------
! This routine builds a grid based on terrain-following coordinates.
! The module parameter coordinateResolution(:) determines the resolution in
! sigma coordinate, dSigma(:). sigma-coordinates are defined by
!   sigma = (eta-z)/(H+eta)  s.t. sigma=0 at z=eta and sigma=1 at z=-H .
!------------------------------------------------------------------------------
  
  ! Arguments
  type(regridding_CS),                           intent(in)    :: CS
  type(ocean_grid_type),                         intent(in)    :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_),        intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_, NK_INTERFACE_), intent(inout) :: dzInterface
  
  ! Local variables
  integer :: i, j, k
  integer :: nz
  real    :: nominalDepth, totalThickness, dh
  real, dimension(SZK_(G)+1) :: zOld, zNew

  nz = G%ke
  
  do i = G%isc-1,G%iec+1
    do j = G%jsc-1,G%jec+1
      
      ! Determine water column height
      totalThickness = 0.0
      do k = 1,nz
        totalThickness = totalThickness + h(i,j,k)
      end do
          
      ! The rest of the model defines grids integrating up from the bottom
      nominalDepth = G%bathyT(i,j)*G%m_to_H
      zOld(nz+1) = - nominalDepth
      zNew(nz+1) = - nominalDepth
      do k = nz,1,-1
        zNew(k) = zNew(k+1) + ( totalThickness * CS%coordinateResolution(k) )
        ! Adjust interface position to accomodate inflating layers
        ! without disturbing the interface above
        if ( zNew(k) < (zNew(k+1) + CS%min_thickness) ) then
          zNew(k) = zNew(k+1) + CS%min_thickness
        endif
        zOld(k) = zOld(k+1) + h(i,j,k)
      enddo

      ! Define regridding in terms of a movement of interfaces
      dzInterface(i,j,1) = 0.
      do k = 2,nz
        dzInterface(i,j,k) = zNew(k) - zOld(k)
      enddo
      dzInterface(i,j,nz+1) = 0.

#ifdef __DO_SAFETY_CHECKS__
      dh=max(nominalDepth,totalThickness)
      if (abs(zNew(1)-zOld(1))>(nz-1)*0.5*epsilon(dh)*dh) then
        write(0,*) 'min_thickness=',CS%min_thickness
        write(0,*) 'nominalDepth=',nominalDepth,'totalThickness=',totalThickness
        write(0,*) 'dzInterface(1) = ',dzInterface(i,j,1),epsilon(dh),nz
        do k=1,nz+1
          write(0,*) k,zOld(k),zNew(k)
        enddo
        do k=1,nz
          write(0,*) k,h(i,j,k),zNew(k)-zNew(k+1),totalThickness*CS%coordinateResolution(k),CS%coordinateResolution(k)
        enddo
        call MOM_error( FATAL, &
               'MOM_regridding, buildGridSigma: top surface has moved!!!' )
      endif
      dzInterface(i,j,1) = 0.
      dzInterface(i,j,nz+1) = 0.
#endif

    end do
  end do
  
end subroutine buildGridSigma


!------------------------------------------------------------------------------
! Build grid based on target interface densities
!------------------------------------------------------------------------------
subroutine buildGridRho( G, h, tv, dzInterface, remapCS, CS )
!------------------------------------------------------------------------------
! This routine builds a new grid based on a given set of target interface
! densities (these target densities are computed by taking the mean value
! of given layer densities). The algorithn operates as follows within each
! column:
! 1. Given T & S within each layer, the layer densities are computed.
! 2. Based on these layer densities, a global density profile is reconstructed
!    (this profile is monotonically increasing and may be discontinuous)
! 3. The new grid interfaces are determined based on the target interface
!    densities.
! 4. T & S are remapped onto the new grid.
! 5. Return to step 1 until convergence or until the maximum number of
!    iterations is reached, whichever comes first.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(ocean_grid_type),                         intent(in)    :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_),        intent(in)    :: h
  type(thermo_var_ptrs),                         intent(in)    :: tv     
  real, dimension(NIMEM_,NJMEM_, NK_INTERFACE_), intent(inout) :: dzInterface
  type(remapping_CS),                            intent(in)    :: remapCS
  type(regridding_CS),                           intent(in)    :: CS
  
  ! Local variables
  integer   :: i, j, k, m
  integer   :: map_index
  integer   :: nz
  integer   :: k_found
  integer   :: count_nonzero_layers
  real      :: deviation            ! When iterating to determine the final 
                                    ! grid, this is the deviation between two
                                    ! successive grids.
  real      :: threshold
  real      :: max_thickness
  real      :: correction
  real, dimension(CS%nk,2) :: ppoly_i_E            !Edge value of polynomial
  real, dimension(CS%nk,2) :: ppoly_i_S            !Edge slope of polynomial
  real, dimension(CS%nk,CS%degree_i+1) :: ppoly_i_coefficients !Coefficients of polynomial
  real, dimension(SZK_(G)) :: p_column, densities, T_column, S_column, Tmp_column
  integer, dimension(SZK_(G)) :: mapping
  real    :: nominalDepth, totalThickness, dh
  real, dimension(SZK_(G)+1) :: zOld, zNew
  real, dimension(SZK_(G)) :: h0, h1, hTmp
  real, dimension(SZK_(G)+1) :: x0, x1, xTmp, dx
                                    
  nz = G%ke
  threshold = CS%min_thickness
  p_column(:) = 0.

  ! Build grid based on target interface densities
  do i = G%isc-1,G%iec+1
    do j = G%jsc-1,G%jec+1
    
      ! Copy T and S onto new variables so as to not alter the original values
      ! of T and S (these are remapped at the end of the regridding iterations
      ! once the final grid has been determined).
      T_column = tv%T(i,j,:)
      S_column = tv%S(i,j,:)
        
      ! Copy original grid
      h0(1:nz) = h(i,j,1:nz)
      x0(1) = 0.0
      do k = 1,nz
        x0(k+1) = x0(k) + h0(k)
      end do
    
      ! Start iterations to build grid
      m = 1
      deviation = 1e10
      do while ( ( m .le. NB_REGRIDDING_ITERATIONS ) .and. &
                 ( deviation .gt. DEVIATION_TOLERANCE ) )

        ! Count number of nonzero layers within current water column
        count_nonzero_layers = 0
        do k = 1,nz
          if ( h0(k) .gt. threshold ) then
            count_nonzero_layers = count_nonzero_layers + 1
          end if
        end do

        ! If there is at most one nonzero layer, stop here (no regridding)
        if ( count_nonzero_layers .le. 1 ) then 
          h1(1:nz) = h0(1:nz)
          exit  ! stop iterations here
        end if

        ! Build new grid containing only nonzero layers
        map_index = 1
        correction = 0.0
        do k = 1,nz
          if ( h0(k) .gt. threshold ) then
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
          if ( hTmp(k) .gt. max_thickness ) then
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
        call calculate_density( T_column, S_column, p_column, densities,&
                                 1, nz, tv%eqn_of_state )

        do k = 1,count_nonzero_layers
          densities(k) = densities(mapping(k))
        end do

        ! One regridding iteration
        call regridding_iteration( densities, CS%coordinateInterfaces, CS,&
                                   count_nonzero_layers, hTmp, xTmp, &
                                   ppoly_i_E, ppoly_i_S, ppoly_i_coefficients, &
                                   nz, h1, x1 )
        
        ! Remap T and S from previous grid to new grid
        do k = 1,nz
          h0(k) = x0(k+1) - x0(k)
          h1(k) = x1(k+1) - x1(k)
        end do
        dx(:) = x1(:) - x0(:)
        dx(1) = 0.
        dx(nz+1) = 0.
        
        call remapping_core(remapCS, nz, h0, S_column, nz, dx, Tmp_column)
        S_column(:) = Tmp_column(:)
        
        call remapping_core(remapCS, nz, h0, T_column, nz, dx, Tmp_column)
        T_column(:) = Tmp_column(:)

        ! Compute the deviation between two successive grids
        deviation = 0.0
        do k = 2,nz
          deviation = deviation + (x0(k)-x1(k))**2
        end do
        deviation = sqrt( deviation / (nz-1) )

    
        m = m + 1
        
        ! Copy final grid onto start grid for next iteration
        x0(:) = x1(:)
        h0(:) = h1(:)

      end do ! end regridding iterations               

      ! Local depth (G%bathyT is positive)
      nominalDepth = G%bathyT(i,j)*G%m_to_H

      ! The rest of the model defines grids integrating up from the bottom
      totalThickness = 0.0
      zOld(nz+1) = - nominalDepth
      zNew(nz+1) = - nominalDepth
      do k = nz,1,-1
        totalThickness = totalThickness + h(i,j,k)
        zNew(k) = zNew(k+1) + h1(k)
        ! Adjust interface position to accomodate inflating layers
        ! without disturbing the interface above
  !     if ( zNew(k) < (zNew(k+1) + CS%min_thickness) ) then
  !       zNew(k) = zNew(k+1) + CS%min_thickness
  !     endif
        zOld(k) = zOld(k+1) + h(i,j,k)
      enddo

      ! Define regridding in terms of a movement of interfaces
      dzInterface(i,j,1) = 0.
      do k = 2,nz
        dzInterface(i,j,k) = zNew(k) - zOld(k)
#ifdef __DO_SAFETY_CHECKS__
        if (zNew(k) > zOld(1)) then
          write(0,*) 'zOld=',zOld
          write(0,*) 'zNew=',zNew
          call MOM_error( FATAL, 'MOM_regridding, buildGridRho: '//&
               'interior interface above surface!' )
        endif
        if (zNew(k) > zNew(k-1)) then
          write(0,*) 'zOld=',zOld
          write(0,*) 'zNew=',zNew
          call MOM_error( FATAL, 'MOM_regridding, buildGridRho: '//&
               'interior interfaces cross!' )
        endif
#endif
      enddo
      dzInterface(i,j,nz+1) = 0.

#ifdef __DO_SAFETY_CHECKS__
      dh=max(nominalDepth,totalThickness)
      if (abs(zNew(1)-zOld(1))>(nz-1)*0.5*epsilon(dh)*dh) then
        write(0,*) 'min_thickness=',CS%min_thickness
        write(0,*) 'nominalDepth=',nominalDepth,'totalThickness=',totalThickness
        write(0,*) 'zNew(1)-zOld(1) = ',zNew(1)-zOld(1),epsilon(dh),nz
        do k=1,nz+1
          write(0,*) k,zOld(k),zNew(k)
        enddo
        do k=1,nz
          write(0,*) k,h(i,j,k),zNew(k)-zNew(k+1)
        enddo
        call MOM_error( FATAL, &
               'MOM_regridding, buildGridRho: top surface has moved!!!' )
      endif
#endif

    end do  ! end loop on j 
  end do  ! end loop on i

end subroutine buildGridRho


!------------------------------------------------------------------------------
! Build arbitrary grid
!------------------------------------------------------------------------------
subroutine build_grid_arbitrary( G, h, dzInterface, h_new, CS )
!------------------------------------------------------------------------------
! This routine builds a grid based on arbitrary rules
!------------------------------------------------------------------------------
  
  ! Arguments
  type(ocean_grid_type),                         intent(in)    :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_),        intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_, NK_INTERFACE_), intent(inout) :: dzInterface
  real,                                          intent(inout) :: h_new
  type(regridding_CS),                           intent(in)    :: CS
  
  ! Local variables
  integer   :: i, j, k
  integer   :: nz
  real      :: z_inter(SZK_(G)+1)
  real      :: total_height
  real      :: delta_h
  real      :: max_depth
  real      :: min_thickness
  real      :: eta              ! local elevation
  real      :: local_depth
  real      :: x1, y1, x2, y2
  real      :: x, t

  nz = G%ke
  
  max_depth = G%max_depth
  min_thickness = CS%min_thickness

  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      ! Local depth
      local_depth = G%bathyT(i,j)*G%m_to_H
      
      ! Determine water column height
      total_height = 0.0
      do k = 1,nz
        total_height = total_height + h(i,j,k)
      end do
          
      eta = total_height - local_depth
      
      ! Compute new thicknesses based on stretched water column
      delta_h = (max_depth + eta) / nz
      
      ! Define interfaces
      z_inter(1) = eta
      do k = 1,nz
        z_inter(k+1) = z_inter(k) - delta_h
      end do

      ! Refine grid in the middle
      do k = 1,nz+1
        x1 = 0.35; y1 = 0.45; x2 = 0.65; y2 = 0.55
    
        x = - ( z_inter(k) - eta ) / max_depth
      
        if ( x .le. x1 ) then
          t = y1*x/x1
        else if ( (x .gt. x1 ) .and. ( x .lt. x2 )) then
          t = y1 + (y2-y1) * (x-x1) / (x2-x1)
        else  
          t = y2 + (1.0-y2) * (x-x2) / (1.0-x2)
        end if

        z_inter(k) = -t * max_depth + eta
      
      end do
    
      ! Modify interface heights to account for topography
      z_inter(nz+1) = - local_depth

      ! Modify interface heights to avoid layers of zero thicknesses
      do k = nz,1,-1
        if ( z_inter(k) .LT. (z_inter(k+1) + min_thickness) ) then
          z_inter(k) = z_inter(k+1) + min_thickness
        end if
      end do

      ! Chnage in interface position
      x = 0. ! Left boundary at x=0
      dzInterface(i,j,1) = 0.
      do k = 2,nz
        x = x + h(i,j,k)
        dzInterface(i,j,k) = z_inter(k) - x
      end do
      dzInterface(i,j,nz+1) = 0.
    
    end do
  end do
  
stop 'OOOOOOPS' ! For some reason the gnu compiler will not let me delete this
                ! routine????
  
end subroutine build_grid_arbitrary


!------------------------------------------------------------------------------
! Regridding iterations
!------------------------------------------------------------------------------
subroutine regridding_iteration( densities, target_values, CS, &
                                 n0, h0, x0, ppoly0_E, ppoly0_S, ppoly0_coefficients, n1, h1, x1 )
! ------------------------------------------------------------------------------
! This routine performs one single iteration of the regridding based on target
! interface densities. Given the set of target values and cell densities, this
! routine builds an interpolated profile and determines the location of the new
! grid. It may happen that, given a high-order interpolator, the number of
! available layers is insufficient (e.g., there are two available layers for
! a third-order PPM ih4 scheme). In these cases, we resort to the simplest 
! continuous linear scheme (P1M h2).
! ------------------------------------------------------------------------------
  
  ! Arguments
  real, dimension(:),  intent(in)    :: densities ! Actual cell densities
  real, dimension(:),  intent(in)    :: target_values ! Target interface densities
  integer,             intent(in)    :: n0 ! Number of cells on source grid
  real, dimension(:),  intent(in)    :: h0 ! cell widths on source grid
  real, dimension(:),  intent(in)    :: x0 ! interface positions on source grid
  real, dimension(:,:),intent(inout) :: ppoly0_E            !Edge value of polynomial
  real, dimension(:,:),intent(inout) :: ppoly0_S            !Edge slope of polynomial
  real, dimension(:,:),intent(inout) :: ppoly0_coefficients !Coefficients of polynomial 
  integer,             intent(in)    :: n1 ! Number of cells on target grid
  real, dimension(:),  intent(out)   :: h1 ! cell widths on target grid
  real, dimension(:),  intent(out)   :: x1 ! interface positions on target grid
  type(regridding_CS), intent(in)    :: CS   ! Parameters used for regridding

  ! Local variables
  integer :: degree, k

  ! Reset piecewise polynomials
  ppoly0_E = 0.0
  ppoly0_S = 0.0
  ppoly0_coefficients = 0.0
  
  ! 1. Compute the interpolated profile of the density field and build grid
  select case ( CS%interpolation_scheme )
  
    case ( INTERPOLATION_P1M_H2 )
      degree = DEGREE_1
      call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
      call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
      if ( CS%boundary_extrapolation) then
        call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
      end if    
    
    case ( INTERPOLATION_P1M_H4 )
      degree = DEGREE_1
      if ( n0 .ge. 4 ) then
        call edge_values_explicit_h4( n0, h0, densities, ppoly0_E )
      else
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
      end if
      call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
      if ( CS%boundary_extrapolation) then
        call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
      end if    
    
    case ( INTERPOLATION_P1M_IH4 )
      degree = DEGREE_1
      if ( n0 .ge. 4 ) then
        call edge_values_implicit_h4( n0, h0, densities, ppoly0_E )
      else
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
      end if
      call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
      if ( CS%boundary_extrapolation) then
        call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
      end if    
    
    case ( INTERPOLATION_PLM )   
      degree = DEGREE_1
      call PLM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
      if ( CS%boundary_extrapolation) then
        call PLM_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
      end if    
    
    case ( INTERPOLATION_PPM_H4 )
      if ( n0 .ge. 4 ) then
        degree = DEGREE_2
        call edge_values_explicit_h4( n0, h0, densities, ppoly0_E )
        call PPM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call PPM_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        end if
      end if
    
    case ( INTERPOLATION_PPM_IH4 )

      if ( n0 .ge. 4 ) then
        degree = DEGREE_2
        call edge_values_implicit_h4( n0, h0, densities, ppoly0_E )
        call PPM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call PPM_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        end if  
      end if
    
    case ( INTERPOLATION_P3M_IH4IH3 )
      
      if ( n0 .ge. 4 ) then
        degree = DEGREE_3
        call edge_values_implicit_h4( n0, h0, densities, ppoly0_E )
        call edge_slopes_implicit_h3( n0, h0, densities, ppoly0_S )
        call P3M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_S, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call P3M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_S, ppoly0_coefficients )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        end if  
      end if
      
    case ( INTERPOLATION_P3M_IH6IH5 )
      if ( n0 .ge. 6 ) then
        degree = DEGREE_3
        call edge_values_implicit_h6( n0, h0, densities, ppoly0_E )
        call edge_slopes_implicit_h5( n0, h0, densities, ppoly0_S )
        call P3M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_S, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call P3M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_S, ppoly0_coefficients )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        end if  
      end if
    
    case ( INTERPOLATION_PQM_IH4IH3 )
    
      if ( n0 .ge. 4 ) then
        degree = DEGREE_4
        call edge_values_implicit_h4( n0, h0, densities, ppoly0_E )
        call edge_slopes_implicit_h3( n0, h0, densities, ppoly0_S )
        call PQM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_S, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call PQM_boundary_extrapolation_v1( n0, h0, densities, ppoly0_E, ppoly0_S, ppoly0_coefficients )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        end if  
      end if
    
    case ( INTERPOLATION_PQM_IH6IH5 )
      if ( n0 .ge. 6 ) then
        degree = DEGREE_4
        call edge_values_implicit_h6( n0, h0, densities, ppoly0_E )
        call edge_slopes_implicit_h5( n0, h0, densities, ppoly0_S )
        call PQM_reconstruction( n0, h0, densities, ppoly0_E, ppoly0_S, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call PQM_boundary_extrapolation_v1( n0, h0, densities, ppoly0_E, ppoly0_S, ppoly0_coefficients )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2( n0, h0, densities, ppoly0_E )
        call P1M_interpolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        if ( CS%boundary_extrapolation) then
          call P1M_boundary_extrapolation( n0, h0, densities, ppoly0_E, ppoly0_coefficients )
        end if  
      end if
    
  end select
  
  ! Based on global density profile, interpolate new grid and 
  ! inflate vanished layers    
  call interpolate_grid( n0, h0, x0, ppoly0_E, ppoly0_coefficients, target_values, degree, n1, h1, x1 )
  call inflate_vanished_layers( CS%min_thickness, n1, h1 )
  x1(1) = 0.0
  do k = 1,n1
    x1(k+1) = x1(k) + h1(k)
  end do
    
end subroutine regridding_iteration


!------------------------------------------------------------------------------
! Given target values (e.g., density), build new grid based on polynomial 
!------------------------------------------------------------------------------
subroutine interpolate_grid( n0, h0, x0, ppoly0_E, ppoly0_coefficients, target_values, degree, n1, h1, x1 )
! ------------------------------------------------------------------------------
! Given the grid 'grid0' and the piecewise polynomial interpolant 
! 'ppoly0' (possibly discontinuous), the coordinates of the new grid 'grid1' 
! are determined by finding the corresponding target interface densities.
! ------------------------------------------------------------------------------
  
  ! Arguments
  integer,            intent(in)    :: n0
  real, dimension(:), intent(in)    :: h0
  real, dimension(:), intent(in)    :: x0
  real, dimension(:,:), intent(in)  :: ppoly0_E            !Edge value of polynomial
  real, dimension(:,:), intent(in)  :: ppoly0_coefficients !Coefficients of polynomial  
  real, dimension(:), intent(in)    :: target_values
  integer,            intent(in)    :: degree
  integer,            intent(in)    :: n1
  real, dimension(:), intent(inout) :: h1
  real, dimension(:), intent(inout) :: x1
    
  ! Local variables
  integer        :: k   ! loop index
  real           :: t   ! current interface target density
  
  ! Make sure boundary coordinates of new grid coincide with boundary 
  ! coordinates of previous grid
  x1(1) = x0(1)
  x1(n1+1) = x0(n0+1)
      
  ! Find coordinates for interior target values
  do k = 2,n1
    t = target_values(k)
    x1(k) = get_polynomial_coordinate ( n0, h0, ppoly0_E, ppoly0_coefficients, t, degree )
    h1(k-1) = x1(k) - x1(k-1)
  end do
  h1(n1) = x1(n1+1) - x1(n1)
    
end subroutine interpolate_grid


!------------------------------------------------------------------------------
! Given target value, find corresponding coordinate for given polynomial
!------------------------------------------------------------------------------
real function get_polynomial_coordinate ( N, h, ppoly_E, ppoly_coefficients, target_value, degree )
! ------------------------------------------------------------------------------
! Here, 'ppoly' is assumed to be a piecewise discontinuous polynomial of degree
! 'degree' throughout the domain defined by 'grid'. A target value is given 
! and we need to determine the corresponding grid coordinate to define the 
! new grid. 
!
! If the target value is out of range, the grid coordinate is simply set to
! be equal to one of the boundary coordinates, which defines vanished layers
! near the boundaries.
!
! IT IS ASSUMED THAT THE PIECEWISE POLYNOMIAL IS MONOTONICALLY INCREASING.
! IF THIS IS NOT THE CASE, THE NEW GRID MAY BE ILL-DEFINED.
! 
! It is assumed that the number of cells defining 'grid' and 'ppoly' are the
! same.
! ------------------------------------------------------------------------------

  ! Arguments
  integer,        intent(in)  :: N
  real,           intent(in)  :: h(:) ! size n
  real,           intent(in)  :: ppoly_E(:,:)            !Edge value of polynomial
  real,           intent(in)  :: ppoly_coefficients(:,:) !Coefficients of polynomial
  real,           intent(in)  :: target_value
  integer,        intent(in)  :: degree

  ! Local variables
  integer            :: i, k            ! loop indices
  integer            :: k_found         ! index of target cell
  integer            :: iter
  real               :: x_l, x_r        ! end coordinates of target cell
  real               :: xi0             ! normalized target coordinate
  real, dimension(5) :: a               ! polynomial coefficients
  real               :: numerator
  real               :: denominator
  real               :: delta           ! Newton-Raphson increment
  real               :: x               ! global target coordinate
  real               :: eps                 ! offset used to get away from
                                        ! boundaries
  real               :: g               ! gradient during N-R iterations

  eps = NR_OFFSET
  
  k_found = -1
  x_l = -1.E30

  ! If the target value is outside the range of all values, we
  ! force the target coordinate to be equal to the lowest or
  ! largest value, depending on which bound is overtaken
  if ( target_value .LE. ppoly_E(1,1) ) then
    x = 0. ! Left boundary is at x=0
    get_polynomial_coordinate = x
    return  ! return because there is no need to look further
  end if
  
  ! Since discontinuous edge values are allowed, we check whether the target
  ! value lies between two discontinuous edge values at interior interfaces
    x = 0. ! Left boundary is at x=0
  do k = 2,N
    x = x + h(k-1) ! Position of interface k
    if ( ( target_value .GE. ppoly_E(k-1,2) ) .AND. &
       ( target_value .LE. ppoly_E(k,1) ) ) then
       get_polynomial_coordinate = x
       return   ! return because there is no need to look further
       exit
    end if   
  end do

  ! If the target value is outside the range of all values, we
  ! force the target coordinate to be equal to the lowest or
  ! largest value, depending on which bound is overtaken
  x = x + h(N) ! Position of right boundary
  if ( target_value .GE. ppoly_E(N,2) ) then
    get_polynomial_coordinate = x
    return  ! return because there is no need to look further
  end if
  
  ! At this point, we know that the target value is bounded and does not
  ! lie between discontinuous, monotonic edge values. Therefore,
  ! there is a unique solution. We loop on all cells and find which one
  ! contains the target value. The variable k_found holds the index value
  ! of the cell where the taregt value lies.
  x_r = 0. ! Left boundary is at x=0
  do k = 1,N
    x_l = x_r
    x_r = x_l + h(k)
    if ( ( target_value .GT. ppoly_E(k,1) ) .AND. &
       ( target_value .LT. ppoly_E(k,2) ) ) then
       k_found = k
       exit
    end if   
  end do

  ! At this point, 'k_found' should be strictly positive. If not, this is
  ! a major failure because it means we could not find any target cell
  ! despite the fact that the target value lies between the extremes. It
  ! means there is a major problem with the interpolant. This needs to be
  ! reported.
  if ( k_found .EQ. -1 ) then
      write(*,*) target_value, ppoly_E(1,1), ppoly_E(N,2)
      write(*,*) 'Could not find target coordinate in ' //&
                 '"get_polynomial_coordinate". This is caused by an '//&
                 'inconsistent interpolant (perhaps not monotonically '//&
                 'increasing)'
      call MOM_error( FATAL, 'Aborting execution' )
  end if

  ! Reset all polynomial coefficients to 0 and copy those pertaining to 
  ! the found cell
  a(:) = 0.0
  do i = 1,degree+1
    a(i) = ppoly_coefficients(k_found,i)
  end do
    
  ! Guess value to start Newton-Raphson iterations (middle of cell) 
  xi0 = 0.5
  iter = 1
  delta = 1e10

  ! Newton-Raphson iterations
  do 
    if ( ( iter .GT. NR_ITERATIONS ) .OR. &
         ( abs(delta) .LT. NR_TOLERANCE ) ) then
      exit    
    end if
    
    numerator = a(1) + a(2)*xi0 + a(3)*xi0*xi0 + a(4)*xi0*xi0*xi0 + &
                a(5)*xi0*xi0*xi0*xi0 - target_value
                
    denominator = a(2) + 2*a(3)*xi0 + 3*a(4)*xi0*xi0 + 4*a(5)*xi0*xi0*xi0
  
    delta = - ( numerator ) / &
              ( denominator )

    xi0 = xi0 + delta

    ! Check whether new estimate is out of bounds. If the new estimate is
    ! indeed out of bounds, we manually set it to be equal to the overtaken 
    ! bound with a small offset towards the interior when the gradient of
    ! the function at the boundary is zero (in which case, the Newton-Raphson
    ! algorithm does not converge).
    if ( xi0 .LT. 0.0 ) then
      xi0 = 0.0
      g = a(2)
      if ( g .EQ. 0.0 ) xi0 = xi0 + eps
    end if
    
    if ( xi0 .GT. 1.0 ) then
      xi0 = 1.0
      g = a(2) + 2*a(3) + 3*a(4) + 4*a(5)
      if ( g .EQ. 0.0 ) xi0 = xi0 - eps
    end if

    iter = iter + 1

  end do ! end Newton-Raphson iterations

  x = x_l + xi0 * h(k_found)
    
  get_polynomial_coordinate = x

end function get_polynomial_coordinate


!------------------------------------------------------------------------------
! Check grid integrity
!------------------------------------------------------------------------------
subroutine check_grid_integrity( CS, G, h )
!------------------------------------------------------------------------------
! This routine is called when initializing the regridding options. The 
! objective is to make sure all layers are at least as thick as the minimum
! thickness allowed for regridding purposes (this parameter is set in the
! MOM_input file or defaulted to 1.0e-3). When layers are too thin, they
! are inflated up to the minmum thickness.
!------------------------------------------------------------------------------

  ! Arguments
  type(regridding_CS),                    intent(in)    :: CS
  type(ocean_grid_type),                  intent(in)    :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout) :: h

  ! Local variables
  integer :: i, j, k
  real    :: hTmp(G%ke)

  do i = G%isc-1,G%iec+1
    do j = G%jsc-1,G%jec+1
    
      ! Build grid for current column
      do k = 1,G%ke
        hTmp(k) = h(i,j,k)
      end do

      call inflate_vanished_layers( CS%min_thickness, G%ke, hTmp )

      ! Save modified grid
      do k = 1,G%ke
        h(i,j,k) = hTmp(k)
      end do
    
    end do
  end do

end subroutine check_grid_integrity


!------------------------------------------------------------------------------
! Inflate vanished layers to finite (nonzero) width
!------------------------------------------------------------------------------
subroutine inflate_vanished_layers( minThickness, N, h )

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
    if ( h(k) .GT. minThickness ) then
      count_nonzero_layers = count_nonzero_layers + 1
    end if
  end do

  ! If all layer thicknesses are greater than the threshold, exit routine
  if ( count_nonzero_layers .eq. N ) return

  ! If all thicknesses are zero, inflate them all and exit
  if ( count_nonzero_layers .eq. 0 ) then  
    do k = 1,N
      h(k) = minThickness
    end do
    return
  end if    
  
  ! Inflate zero layers
  correction = 0.0
  do k = 1,N
    if ( h(k) .le. minThickness ) then
      delta = minThickness - h(k)
      correction = correction + delta
      h(k) = h(k) + delta
    end if  
  end do
  
  ! Modify thicknesses of nonzero layers to ensure volume conservation
  maxThickness = h(1)
  k_found = 1
  do k = 1,N
    if ( h(k) .gt. maxThickness ) then
      maxThickness = h(k)
      k_found = k
    end if  
  end do
  
  h(k_found) = h(k_found) - correction
  
end subroutine inflate_vanished_layers


!------------------------------------------------------------------------------
! Convective adjustment by swapping layers
!------------------------------------------------------------------------------
subroutine convective_adjustment(G, h, tv)
!------------------------------------------------------------------------------
! Check each water column to see whether it is stratified. If not, sort the
! layers by successive swappings of water masses (bubble sort algorithm)
!------------------------------------------------------------------------------

  ! Arguments
  type(ocean_grid_type), intent(in)                  :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout) :: h
  type(thermo_var_ptrs), intent(inout)               :: tv     
  
  ! Local variables
  integer   :: i, j, k
  real      :: T0, T1       ! temperatures
  real      :: S0, S1       ! salinities
  real      :: r0, r1       ! densities
  real      :: h0, h1
  logical   :: stratified
  real, dimension(G%kE) :: p_column, densities

  p_column(:) = 0.
  
  ! Loop on columns 
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1
        
      ! Compute densities within current water column
      call calculate_density( tv%T(i,j,:), tv%S(i,j,:), p_column, &
                              densities, 1, G%ke, tv%eqn_of_state )
     
      ! Repeat restratification until complete  
      do

        stratified = .true.
        do k = 1,G%ke-1
          ! Gather information of current and next cells
          T0 = tv%T(i,j,k)
          T1 = tv%T(i,j,k+1)
          S0 = tv%S(i,j,k)
          S1 = tv%S(i,j,k+1)
          r0 = densities(k)
          r1 = densities(k+1)
          h0 = h(i,j,k)
          h1 = h(i,j,k+1)
          ! If the density of the current cell is larger than the density
          ! below it, we swap the cells and recalculate the densitiies
          ! within the swapped cells    
          if ( r0 .gt. r1 ) then
            tv%T(i,j,k)   = T1
            tv%T(i,j,k+1) = T0
            tv%S(i,j,k)   = S1
            tv%S(i,j,k+1) = S0
            h(i,j,k)      = h1
            h(i,j,k+1)    = h0
            ! Recompute densities at levels k and k+1
            call calculate_density( tv%T(i,j,k), tv%S(i,j,k), &
                                     p_column(k), &
                                     densities(k), tv%eqn_of_state )
            call calculate_density( tv%T(i,j,k+1), tv%S(i,j,k+1), &
                                     p_column(k+1), &
                                     densities(k+1), tv%eqn_of_state )
            stratified = .false.
          end if
        end do  ! k 
    
        if ( stratified ) exit        

      end do    

    end do  ! i
  end do  ! j   

end subroutine convective_adjustment


!------------------------------------------------------------------------------
! Return uniform resolution vector based on coordiante mode
!------------------------------------------------------------------------------
function uniformResolution(nk,coordMode,maxDepth,rhoLight,rhoHeavy)
!------------------------------------------------------------------------------
! Calculate a vector of uniform resolution in the units of the coordinate
!------------------------------------------------------------------------------
  ! Arguments
  integer,          intent(in) :: nk
  character(len=*), intent(in) :: coordMode
  real,             intent(in) :: maxDepth, rhoLight, rhoHeavy
  real                         :: uniformResolution(nk)

  ! Local variables
  integer :: scheme
  
  scheme = coordinateMode(coordMode)
  select case ( scheme )

    case ( REGRIDDING_ZSTAR )
      uniformResolution(:) = maxDepth / real(nk)

    case ( REGRIDDING_RHO )  
      uniformResolution(:) = (rhoHeavy - rhoLight) / real(nk)

    case ( REGRIDDING_SIGMA )
      uniformResolution(:) = 1. / real(nk)

    case default
      call MOM_error(FATAL, "MOM_regridding, uniformResolution: "//&
       "Unrecognized choice for coordinate mode ("//trim(coordMode)//").")

  end select ! type of grid 

end function uniformResolution


!------------------------------------------------------------------------------
! Set the fixed resolution data
!------------------------------------------------------------------------------
subroutine setCoordinateResolution( dz, CS )
  real, dimension(:),  intent(in)    :: dz
  type(regridding_CS), intent(inout) :: CS

  if (size(dz)/=CS%nk) call MOM_error( FATAL, &
      'setCoordinateResolution: inconsistent number of levels' )

  CS%coordinateResolution(:) = dz(:)
  
end subroutine setCoordinateResolution

!-----------------------------------------------------------------------------
! Set the running sum of coordinateResolution
!-----------------------------------------------------------------------------
subroutine setCoordinateInterfaces( G, CS )
  type(ocean_grid_type),  intent(in) :: G      ! Ocean grid informations
  type(regridding_CS), intent(inout) :: CS
  integer :: k, nz

  nz = CS%nk
  CS%coordinateInterfaces(1)    = G%Rlay(1)+0.5*(G%Rlay(1)-G%Rlay(2))
  CS%coordinateInterfaces(nz+1) = G%Rlay(nz)+0.5*(G%Rlay(nz)-G%Rlay(nz-1))
  do k = 2,nz
    CS%coordinateInterfaces(k) = CS%coordinateInterfaces(k-1) + CS%coordinateResolution(k)
  end do

end subroutine setCoordinateInterfaces

!------------------------------------------------------------------------------
! Query the fixed resolution data
!------------------------------------------------------------------------------
function getCoordinateResolution( CS )
  type(regridding_CS), intent(in) :: CS
  real, dimension(CS%nk)          :: getCoordinateResolution

  getCoordinateResolution(:) = CS%coordinateResolution(:)
  
end function getCoordinateResolution

!------------------------------------------------------------------------------
! Query the target coordinate interfaces positions
!------------------------------------------------------------------------------
function getCoordinateInterfaces( CS )
  type(regridding_CS), intent(in) :: CS
  real, dimension(CS%nk+1)        :: getCoordinateInterfaces

  integer :: k

  getCoordinateInterfaces(1) = 0.
  do k = 1, CS%nk
    getCoordinateInterfaces(k+1) = getCoordinateInterfaces(k) &
                                  -CS%coordinateResolution(k)
  enddo
  ! The following line has an "abs()" to allow ferret users to reference
  ! data by index. It is a temporary work around...  :(  -AJA
  getCoordinateInterfaces(:) = abs( getCoordinateInterfaces(:) )

end function getCoordinateInterfaces

!------------------------------------------------------------------------------
! Query the target coordinate units
!------------------------------------------------------------------------------
function getCoordinateUnits( CS )
  type(regridding_CS), intent(in) :: CS
  character(len=20)               :: getCoordinateUnits

  select case ( CS%regridding_scheme )
    case ( REGRIDDING_ZSTAR )
      getCoordinateUnits = 'meter'
    case ( REGRIDDING_SIGMA )
      getCoordinateUnits = 'fraction'
    case ( REGRIDDING_RHO )  
      getCoordinateUnits = 'kg/m3'
    case ( REGRIDDING_ARBITRARY )
      getCoordinateUnits = 'unknown'
    case default
      call MOM_error(FATAL,'MOM_regridding, getCoordinateUnits: '//&
                     'Unknown regridding scheme selected!')
  end select ! type of grid 

end function getCoordinateUnits

!------------------------------------------------------------------------------
! Query the short name of the coordinate
!------------------------------------------------------------------------------
function getCoordinateShortName( CS )
  type(regridding_CS), intent(in) :: CS
  character(len=20)               :: getCoordinateShortName

  select case ( CS%regridding_scheme )
    case ( REGRIDDING_ZSTAR )
      !getCoordinateShortName = 'z*'
      ! The following line is a temporary work around...  :(  -AJA
      getCoordinateShortName = 'pseaduo-depth, -z*'
    case ( REGRIDDING_SIGMA )
      getCoordinateShortName = 'sigma'
    case ( REGRIDDING_RHO )  
      getCoordinateShortName = 'rho'
    case ( REGRIDDING_ARBITRARY )
      getCoordinateShortName = 'coordinate'
    case default
      call MOM_error(FATAL,'MOM_regridding, getCoordinateShortName: '//&
                     'Unknown regridding scheme selected!')
  end select ! type of grid 

end function getCoordinateShortName

!------------------------------------------------------------------------------
! Control the extrapolation of boundary data
!------------------------------------------------------------------------------
subroutine setRegriddingBoundaryExtrapolation( onOff, CS )
  logical,             intent(in)    :: onOff
  type(regridding_CS), intent(inout) :: CS

  CS%boundary_extrapolation = onOff
  
end subroutine setRegriddingBoundaryExtrapolation

!------------------------------------------------------------------------------
! Control the minimum thickness permitted in regridding
!------------------------------------------------------------------------------
subroutine setRegriddingMinimumThickness( minThickness, CS )
  real   ,             intent(in)    :: minThickness
  type(regridding_CS), intent(inout) :: CS

  CS%min_thickness = minThickness
  
end subroutine setRegriddingMinimumThickness

!------------------------------------------------------------------------------
! Return coordinate-derived thicknesses for fixed coordinate systems
!------------------------------------------------------------------------------
function getStaticThickness( CS, SSH, depth )
  type(regridding_CS), intent(in) :: CS
  real,                intent(in) :: SSH
  real,                intent(in) :: depth
  real, dimension(CS%nk)          :: getStaticThickness
  ! Local
  integer :: k
  real :: z, dz

  select case ( CS%regridding_scheme )
    case ( REGRIDDING_ZSTAR )
      if (depth>0.) then
        z = ssh
        do k = 1, CS%nk
          dz = CS%coordinateResolution(k) * ( 1. + ssh/depth ) ! Nominal dz*
          dz = max(dz, 0.)                                     ! Avoid negative incase ssh=-depth
          dz = min(dz, depth - z)                              ! Clip if below topography
          z = z + dz                                           ! Bottom of layer
          getStaticThickness(k) = dz
        enddo
      else
        getStaticThickness(:) = 0. ! On land ...
      endif
    case ( REGRIDDING_SIGMA )
      getStaticThickness(:) = CS%coordinateResolution(:) * ( depth + ssh )
    case ( REGRIDDING_RHO )  
      getStaticThickness(:) = 0. ! Not applicable
    case ( REGRIDDING_ARBITRARY )
      getStaticThickness(:) = 0.  ! Not applicable
    case default
      call MOM_error(FATAL,'MOM_regridding, getStaticThickness: '//&
                     'Unknown regridding scheme selected!')
  end select ! type of grid 
  
end function getStaticThickness

!------------------------------------------------------------------------------
! Allocate memory for regridding
!------------------------------------------------------------------------------
subroutine allocate_regridding( CS )
!------------------------------------------------------------------------------
! In this routine, we allocate the memory needed to carry out regridding
! steps in the course of the simulation. 

! For example, to compute implicit edge-value estimates, a tridiagonal system
! must be solved. We allocate the needed memory at the beginning of the
! simulation because the number of layers never changes.
!------------------------------------------------------------------------------

  ! Arguments
  type(regridding_CS), intent(inout) :: CS

  ! Local variables

  ! Target values
  allocate( CS%coordinateInterfaces(CS%nk+1) )

  ! Target resolution (for fixed coordinates)
  allocate( CS%coordinateResolution(CS%nk) ); CS%coordinateResolution(:) = -1.E30

end subroutine allocate_regridding


!------------------------------------------------------------------------------
! Deallocate memory for regridding
!------------------------------------------------------------------------------
subroutine regridding_memory_deallocation( CS )
!------------------------------------------------------------------------------
! In this routine, we reclaim the memory that was allocated for the regridding. 
!------------------------------------------------------------------------------
  
  type(regridding_CS), intent(inout) :: CS
  
  ! Target values
  deallocate( CS%coordinateInterfaces )
  deallocate( CS%coordinateResolution )

end subroutine regridding_memory_deallocation

end module MOM_regridding
