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

logical, parameter :: doSafetyChecks = .false.

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! General remapping routine 
!------------------------------------------------------------------------------
subroutine remapping_main( CS, G, h, h_new, tv, u, v )
!------------------------------------------------------------------------------
! This routine takes care of remapping all variable between the old and the
! new grids. When velocity components need to be remapped, thicknesses at
! velocity points are taken to be arithmetic averages of tracer thicknesses.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(remapping_CS),                     intent(inout) :: CS
  type(ocean_grid_type),                  intent(in)    :: G
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(in)    :: h_new
  type(thermo_var_ptrs),                  intent(inout) :: tv       
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout), optional  :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout), optional  :: v
  
  ! Local variables
  integer               :: i, j, k
  integer               :: nz
  real                  :: val, new_val
  integer               :: problem

  nz = G%ke
  if (doSafetyChecks) then
    if (nz>CS%grid_start%nb_cells) call MOM_error(FATAL,'nz>nk_start')
    if (nz>CS%grid_final%nb_cells) call MOM_error(FATAL,'nz>nk_final')
  endif

  ! Remap tracer
  do j = G%jsc,G%jec
    do i = G%isc,G%iec
    
      ! Build the start and final grids
      CS%grid_start%h(:) = h(i,j,:)
      CS%grid_final%h(:) = h_new(i,j,:)
      call buildFinalGrid(nz, CS%grid_start%h, CS%grid_final%h, CS%grid_start%x, CS%grid_final%x)
      
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
      call buildFinalGrid(nz, CS%grid_start%h, CS%grid_final%h, CS%grid_start%x, CS%grid_final%x)
      
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
      call buildFinalGrid(nz, CS%grid_start%h, CS%grid_final%h, CS%grid_start%x, CS%grid_final%x)

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
subroutine buildFinalGrid(nz, hs, hf, xs, xf)
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
  xs(1) = 0.0
  do k = 1,nz
    xs(k+1) = xs(k) + hs(k)
  end do
  sumH1 = xs(nz+1)

  ! Initial guess at final grid
  xf(1) = 0.0
  do k = 1,nz
    xf(k+1) = xf(k) + hf(k)
  end do
  sumH2 = xf(nz+1)

  if (doSafetyChecks) then
    if (abs(sumH1-sumH2)>0.5*real(nz)*epsilon(sumH2)*(sumH1+sumH2)) then
      write(0,*) 'Start/final/start-final grid'
      do k = 1,nz+1
        write(0,'(i4,3es12.3)') k,xs(k),xf(k),xs(k)-xf(k)
      enddo
      write(0,*) 'eps,H*eps',epsilon(sumH2),0.5*epsilon(sumH2)*(sumH1+sumH2)
      call MOM_error(FATAL,'MOM_remapping, buildFinalGrid: '//&
                     'Final and start grids do not match.')
    endif
  endif

  ! Adjust new grid to match start grid to handle round-off differences
  ! This ensures that the final grid matches the start grid at top and bottom
! if (sumH1/=sumH2) then
!   xf(nz+1) = xs(nz+1)
!   do k = nz,1,-1
!     nonDimPos = xf(k) / sumH2
!     xf(k) = (1.-nonDimPos) * xf(k) + nonDimPos  * (xf(k+1) - hf(k))
!   end do
! endif

end subroutine buildFinalGrid


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
  
  ! Reset polynomial
  CS%ppoly_r%E(:,:) = 0.0
  CS%ppoly_r%S(:,:) = 0.0
  CS%ppoly_r%coefficients(:,:) = 0.0

  select case ( CS%remapping_scheme )
    case ( REMAPPING_PCM )
      call PCM_reconstruction( grid0, u0, CS%ppoly_r )
      call remapping_integration( grid0, u0, CS%ppoly_r, grid1, u1, &
                                   INTEGRATION_PCM )
    case ( REMAPPING_PLM )
      call PLM_reconstruction( grid0, u0, CS%ppoly_r )
      if ( CS%boundary_extrapolation) then
        call PLM_boundary_extrapolation( grid0, u0, CS%ppoly_r )
      end if    
      call remapping_integration( grid0, u0, CS%ppoly_r, grid1, u1, &
                                   INTEGRATION_PLM )
    case ( REMAPPING_PPM_H4 )
      call edge_values_explicit_h4( grid0, u0, CS%ppoly_r%E )
      call PPM_reconstruction( grid0, u0, CS%ppoly_r )
      if ( CS%boundary_extrapolation) then
        call PPM_boundary_extrapolation( grid0, u0, CS%ppoly_r )
      end if    
      call remapping_integration( grid0, u0, CS%ppoly_r, grid1, u1, &
                                   INTEGRATION_PPM )
    case ( REMAPPING_PPM_IH4 )
      call edge_values_implicit_h4( grid0, CS%edgeValueWrk, u0, CS%ppoly_r%E )
      call PPM_reconstruction( grid0, u0, CS%ppoly_r )
      if ( CS%boundary_extrapolation) then
        call PPM_boundary_extrapolation( grid0, u0, CS%ppoly_r )
      end if    
      call remapping_integration( grid0, u0, CS%ppoly_r, grid1, u1, &
                                   INTEGRATION_PPM )
    case ( REMAPPING_PQM_IH4IH3 )
      call edge_values_implicit_h4( grid0, CS%edgeValueWrk, u0, CS%ppoly_r%E )
      call edge_slopes_implicit_h3( grid0, CS%edgeSlopeWrk, u0, CS%ppoly_r%S )
      call PQM_reconstruction( grid0, u0, CS%ppoly_r )
      if ( CS%boundary_extrapolation) then
        call PQM_boundary_extrapolation_v1( grid0, u0, CS%ppoly_r )
      end if    
      call remapping_integration( grid0, u0, CS%ppoly_r, grid1, u1, &
                                   INTEGRATION_PQM )
    case ( REMAPPING_PQM_IH6IH5 )
      call edge_values_implicit_h6( grid0, CS%edgeValueWrk, u0, CS%ppoly_r%E )
      call edge_slopes_implicit_h5( grid0, CS%edgeSlopeWrk, u0, CS%ppoly_r%S )
      call PQM_reconstruction( grid0, u0, CS%ppoly_r )
      if ( CS%boundary_extrapolation) then
        call PQM_boundary_extrapolation_v1( grid0, u0, CS%ppoly_r )
      end if    
      call remapping_integration( grid0, u0, CS%ppoly_r, grid1, u1, &
                                   INTEGRATION_PQM )
    case default
      call MOM_error( FATAL, 'The selected remapping method is invalid' )
  end select

end subroutine remapping_core


! -----------------------------------------------------------------------------
! remapping_integration (integration of reconstructed profile)
! -----------------------------------------------------------------------------
subroutine remapping_integration( grid0, u0, ppoly0, grid1, u1, method )
  ! Arguments
  type(grid1D_t), intent(in)        :: grid0    ! source grid
  real, dimension(:), intent(in)    :: u0       ! source cell averages
  type(ppoly_t), intent(in)         :: ppoly0   ! source piecewise polynomial
  type(grid1D_t), intent(in)        :: grid1    ! target grid
  real, dimension(:), intent(inout) :: u1       ! target cell averages
  integer                           :: method   ! remapping scheme to use
  
  ! Local variables
  integer       :: iTarget, j, k
  integer       :: n0, n1       ! nb of cells in grid0 and grid1, respectively
  integer       :: j0, j1       ! indexes of source cells containing target 
                                ! cell edges
  real          :: x0, x1       ! coordinates of target cell edges  
  real          :: q0, q1       ! partially integrated quantities in source 
                                ! cells j0 and j1
  real          :: q            ! complete integration
  real          :: a, b         ! interval of integration (global coordinates)
  real          :: xi0, xi1     ! interval of integration (local -- normalized 
                                ! -- coordinates)

  ! A priori, both grids contains the same number of cells but, who knows...
  n0 = grid0%nb_cells
  n1 = grid1%nb_cells

  ! Loop on cells in target grid (grid1). For each target cell, we need to find
  ! in which source cells the target cell edges lie. The associated indexes are 
  ! noted j0 and j1.
  do iTarget = 1,grid1%nb_cells
    ! Determine the coordinates of the target cell edges
    x0 = grid1%x(iTarget)
    x1 = grid1%x(iTarget+1)

    ! ============================================================
    ! Check whether target cell is vanished. If it is, the cell
    ! average is simply the interpolated value at the location
    ! of the vanished cell. If it isn't, we need to integrate the
    ! quantity within the cell and divide by the cell width to
    ! determine the cell average.
    ! ============================================================
    ! 1. Cell is vanished
    if ( abs(x0 - x1) .EQ. 0.0 ) then
    
      j0 = -1
    
      do j = 1,grid0%nb_cells
        ! Left edge is found in cell j
        if ( ( x0 .GE. grid0%x(j) ) .AND. ( x0 .LE. grid0%x(j+1) ) ) then
          j0 = j
          exit ! once target grid cell is found, exit loop
        end if
      end do

      ! If, at this point, j0 is equal to -1, it means the vanished
      ! cell lies outside the source grid. In other words, it means that
      ! the source and target grids do not cover the same physical domain
      ! and there is something very wrong !
      if ( j0 .EQ. -1 ) then
        write(0,*) iTarget
        call MOM_error(FATAL, 'The location of the vanished cell could '//&
                              'not be found in "remapping_integration"' )
      end if

      ! We check whether the source cell (i.e. the cell in which the
      ! vanished target cell lies) is vanished. If it is, the interpolated 
      ! value is set to be mean of the edge values (which should be the same).
      ! If it isn't, we simply interpolate.
      if ( grid0%h(j0) .EQ. 0.0 ) then
        u1(iTarget) = 0.5 * ( ppoly0%E(j0,1) + ppoly0%E(j0,2) )
      else
        xi0 = x0 / grid0%h(j0) - grid0%x(j0) / grid0%h(j0)
    
        select case ( method )
          case ( INTEGRATION_PCM )   
            u1(iTarget) = ppoly0%coefficients(j0,1)
          case ( INTEGRATION_PLM )  
            u1(iTarget) = evaluation_polynomial( ppoly0%coefficients(j0,:), 2, xi0 )
          case ( INTEGRATION_PPM )
            u1(iTarget) = evaluation_polynomial( ppoly0%coefficients(j0,:), 3, xi0 )
          case ( INTEGRATION_PQM )
            u1(iTarget) = evaluation_polynomial( ppoly0%coefficients(j0,:), 5, xi0 )
          case default
            call MOM_error( FATAL,'The selected integration method is invalid' )
        end select   
        
      end if ! end checking whether source cell is vanished
    
    ! 2. Cell is not vanished
    else

    ! Find the cells in source grid containing the target cell edges
    j0 = -1
    j1 = -1
    
    do j = 1,grid0%nb_cells
    
        ! Left edge is found in cell j
        if ( ( x0 .GE. grid0%x(j) ) .AND. ( x0 .LE. grid0%x(j+1) ) ) then
            j0 = j
            exit   ! once target grid cell is found, exit loop
        end if
        
    end do  

    do j = 1,grid0%nb_cells
        
        ! Right edge is found in cell j
        if ( ( x1 .GE. grid0%x(j) ) .AND. ( x1 .LE. grid0%x(j+1) ) ) then
            j1 = j
            exit  ! once target grid cell is found, exit loop
        end if
        
    end do ! end loop on source grid cells

    ! Here, we make sure that the boundary edges of boundary cells
    ! coincide
    if ( iTarget .EQ. 1 ) then
      j0 = 1
    end if  
    
    if ( iTarget .EQ. grid1%nb_cells ) then
      j1 = grid1%nb_cells
    end if  

    ! To integrate, two cases must be considered: (1) the target cell is
    ! entirely comtained within a cell of the source grid and (2) the target
    ! cell spans at least two cells of the source grid.

    if ( j0 .EQ. j1 ) then
    ! The target cell is entirely contained within a cell of the source
    ! grid. This situation is represented by the following schematic, where
    ! the cell in which x0 and x1 are located has index j0=j1 :
    ! 
    ! ----|-----o--------o----------|-------------
    !           x0       x1
    !
      ! Determine normalized coordinates
      xi0 = x0 / grid0%h(j0) - grid0%x(j0) / grid0%h(j0)
      xi1 = x1 / grid0%h(j0) - grid0%x(j0) / grid0%h(j0)

      ! Depending on which polynomial is used, integrate quantity
      ! between x0 and xi1. Integration is carried out in normalized
      ! coordinates, hence: \int_x0^x1 p(x) dx = h \int_xi0^xi1 p(xi) dxi
      select case ( method )
        case ( INTEGRATION_PCM )     
          q = ppoly0%coefficients(j0,1) * ( x1 - x0 )
        case ( INTEGRATION_PLM )    
          q = grid0%h(j0) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(j0,:), 1 )
        case ( INTEGRATION_PPM )
          q = grid0%h(j0) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(j0,:), 2 )
        case ( INTEGRATION_PQM )
          q = grid0%h(j0) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(j0,:), 4 )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select     
    
    else
    ! The target cell spans at least two cells of the source grid.
    ! This situation is represented by the following schematic, where
    ! the cells in which x0 and x1 are located have indexes j0 and j1,
    ! respectively :
    ! 
    ! ----|-----o---|--- ... --|---o----------|-------------
    !           x0                 x1
    !
    ! We first integrate from x0 up to the right boundary of cell j0, then
    ! add the integrated amounts of cells located between j0 and j1 and then
    ! integrate from the left boundary of cell j1 up to x1

      q = 0.0

      ! Integrate from x0 up to right boundary of cell j0
      !xi0 = x0 / grid0%h(j0) - grid0%x(j0) / grid0%h(j0)
      xi0 = (x0 - grid0%x(j0)) / grid0%h(j0)
      xi1 = 1.0
      select case ( method )
        case ( INTEGRATION_PCM )     
          q = q + ppoly0%coefficients(j0,1) * ( grid0%x(j0+1) - x0 )
        case ( INTEGRATION_PLM )    
          q = q + grid0%h(j0) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(j0,:), 1 )
        case ( INTEGRATION_PPM )
          q = q + grid0%h(j0) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(j0,:), 2 )
        case ( INTEGRATION_PQM )
          q = q + grid0%h(j0) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(j0,:), 4 )
        case default
          call MOM_error( FATAL, 'The selected integration method is invalid' )
      end select     
    
      ! Integrate contents within cells strictly comprised between j0 and j1
      if ( j1 .GT. (j0+1) ) then
        do k = j0+1,j1-1
          q = q + grid0%h(k) * u0(k)
        end do
      end if

      ! Integrate from left boundary of cell j1 up to x1
      xi0 = 0.0
      !xi1 = x1 / grid0%h(j1) - grid0%x(j1) / grid0%h(j1)
      xi1 = (x1 - grid0%x(j1)) / grid0%h(j1)
      
      select case ( method )
        case ( INTEGRATION_PCM )     
          q = q + ppoly0%coefficients(j1,1) * ( x1 - grid0%x(j1) )
        case ( INTEGRATION_PLM )    
          q = q + grid0%h(j1) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(j1,:), 1 )
        case ( INTEGRATION_PPM )
          q = q + grid0%h(j1) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(j1,:), 2 )
        case ( INTEGRATION_PQM )
          q = q + grid0%h(j1) * &
              integration_polynomial( xi0, xi1, ppoly0%coefficients(j1,:), 4 )
        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select     
      
    end if ! end integration for non-vanished cells 
    
    ! The cell average is the integrated value divided by the cell width
    u1(iTarget) = q / grid1%h(iTarget)
    
    end if ! end if clause to check if cell is vanished
    
  end do ! end iTarget loop on target grid cells

end subroutine remapping_integration


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
