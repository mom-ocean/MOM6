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
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use regrid_grid1d_class ! see 'regrid_grid1d_class.F90'
use regrid_ppoly_class  ! see 'regrid_ppoly.F90'
use regrid_polynomial   ! see 'regrid_polynomial.F90'
use regrid_edge_values  ! see 'regrid_edge_values.F90'
use regrid_edge_slopes  ! see 'regrid_edge_slopes.F90'
use regrid_pcm          ! see 'regrid_pcm.F90'
use regrid_plm          ! see 'regrid_plm.F90'
use regrid_ppm          ! see 'regrid_ppm.F90'
use regrid_pqm          ! see 'regrid_pqm.F90'
use regrid_p1m          ! see 'regrid_p1m.F90'
use regrid_p3m          ! see 'regrid_p3m.F90'
use regrid_defs         ! see 'regrid_defs.F90' (contains types and parameters)
            
implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Private (module-wise) variables
! -----------------------------------------------------------------------------
type(grid1d_t)                  :: grid_start;      ! starting grid
type(grid1d_t)                  :: grid_final;      ! final grid
type(ppoly_t)                   :: ppoly_r;         ! reconstruction ppoly
real, dimension(:), allocatable :: u_column;        ! generic variable used for
                                                    ! remapping

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public remapping_main
public remapping_core
public remapping_memory_allocation
public remapping_memory_deallocation

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! General remapping routine 
!------------------------------------------------------------------------------
subroutine remapping_main ( G, regridding_opts, h, h_new, tv, u, v )
!------------------------------------------------------------------------------
! This routine takes care of remapping all variable between the old and the
! new grids. When velocity components need to be remapped, thicknesses at
! velocity points are taken to be arithmetic averages of tracer thicknesses.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(ocean_grid_type), intent(in)                       :: G
  type(regridding_opts_t), intent(in)                     :: regridding_opts
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)       :: h
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout)    :: h_new
  type(thermo_var_ptrs), intent(inout)                    :: tv       
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout), optional  :: u
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout), optional  :: v
  
  ! Local variables
  integer               :: i, j, k
  integer               :: nz
  real                  :: val, new_val
  integer               :: problem

  nz = G%ke

  ! Remap tracer
  do j = G%jsc,G%jec
    do i = G%isc,G%iec
    
      ! Build initial grid
      grid_start%h(:) = h(i,j,:)
      
      grid_start%x(1) = 0.0
      do k = 1,nz
        grid_start%x(k+1) = grid_start%x(k) + grid_start%h(k)
      end do

      ! Build new grid
      grid_final%h(:) = h_new(i,j,:)
      
      grid_final%x(1) = 0.0
      do k = 1,nz
        grid_final%x(k+1) = grid_final%x(k) + grid_final%h(k)
      end do
      
      do k = 1,nz
        grid_start%h(k) = grid_start%x(k+1) - grid_start%x(k)
        grid_final%h(k) = grid_final%x(k+1) - grid_final%x(k)
      end do
      
      call remapping_core ( grid_start, tv%S(i,j,:), grid_final, u_column, &
                            ppoly_r, regridding_opts )
      
      tv%S(i,j,:) = u_column(:)
      
      call remapping_core ( grid_start, tv%T(i,j,:), grid_final, u_column, &
                            ppoly_r, regridding_opts )
     
      tv%T(i,j,:) = u_column(:)

    end do
  end do
  
  ! Remap u velocity component
  if ( present(u) ) then
  do j = G%jsc,G%jec
    do i = G%iscB,G%iecB
    
      ! Build initial grid
      grid_start%h(:) = 0.5 * ( h(i,j,:) + h(i+1,j,:) )

      grid_start%x(1) = 0.0
      do k = 1,nz
        grid_start%x(k+1) = grid_start%x(k) + grid_start%h(k)
      end do
      
      ! Build final grid
      grid_final%h(:) = 0.5 * ( h_new(i,j,:) + h_new(i+1,j,:) )
      
      grid_final%x(1) = 0.0
      do k = 1,nz
        grid_final%x(k+1) = grid_final%x(k) + grid_final%h(k)
      end do
      
      do k = 1,nz
        grid_start%h(k) = grid_start%x(k+1) - grid_start%x(k)
        grid_final%h(k) = grid_final%x(k+1) - grid_final%x(k)
      end do
  
      call remapping_core ( grid_start, u(i,j,:), grid_final, u_column, &
                            ppoly_r, regridding_opts )
     
      u(i,j,:) = u_column(:)
      
    end do
  end do
  end if
  
  ! Remap v velocity component
  if ( present(v) ) then
  do j = G%jscB,G%jecB
    do i = G%isc,G%iec

      ! Build initial grid
      grid_start%h(:) = 0.5 * ( h(i,j,:) + h(i,j+1,:) )
      grid_start%x(1) = 0.0
      do k = 1,nz
        grid_start%x(k+1) = grid_start%x(k) + grid_start%h(k)
      end do
      
      ! Build final grid
      grid_final%h(:) = 0.5 * ( h_new(i,j,:) + h_new(i,j+1,:) )
      grid_final%x(1) = 0.0
      do k = 1,nz
        grid_final%x(k+1) = grid_final%x(k) + grid_final%h(k)
      end do

      do k = 1,nz
        grid_start%h(k) = grid_start%x(k+1) - grid_start%x(k)
        grid_final%h(k) = grid_final%x(k+1) - grid_final%x(k)
      end do

      call remapping_core ( grid_start, v(i,j,:), grid_final, u_column, &
                            ppoly_r, regridding_opts )
     
      v(i,j,:) = u_column(:)
      
    end do
  end do
  end if

end subroutine remapping_main


!------------------------------------------------------------------------------
! Remapping core routine
!------------------------------------------------------------------------------
subroutine remapping_core ( grid0, u0, grid1, u1, ppoly, regridding_opts )
!------------------------------------------------------------------------------
! This routine is basic in that it simply takes two grids and remaps the
! field known on the first grid onto the second grid, following the rules
! stored in the structure regridding_opts.
!------------------------------------------------------------------------------

  ! Arguments
  type(grid1d_t), intent(in)          :: grid0
  real, dimension(:), intent(in)      :: u0
  type(ppoly_t), intent(inout)        :: ppoly
  type(grid1d_t), intent(in)          :: grid1
  real, dimension(:), intent(inout)   :: u1
  type(regridding_opts_t), intent(in) :: regridding_opts
  
  ! Reset polynomial
  ppoly%E(:,:) = 0.0
  ppoly%S(:,:) = 0.0
  ppoly%coefficients(:,:) = 0.0

  select case ( regridding_opts%remapping_scheme )

    case ( REMAPPING_PCM )
      
      call pcm_reconstruction ( grid0, ppoly, u0 )
      call remapping_integration ( grid0, u0, ppoly, grid1, u1, &
                                   INTEGRATION_PCM )
    
    case ( REMAPPING_PLM )
      
      call plm_reconstruction ( grid0, ppoly, u0 )
      if ( regridding_opts%boundary_extrapolation) then
        call plm_boundary_extrapolation ( grid0, ppoly, u0 )
      end if    
      call remapping_integration ( grid0, u0, ppoly, grid1, u1, &
                                   INTEGRATION_PLM )
    
    case ( REMAPPING_PPM_H4 )
      
      call edge_values_explicit_h4 ( grid0, u0, ppoly%E )
      call ppm_reconstruction ( grid0, ppoly, u0 )
      if ( regridding_opts%boundary_extrapolation) then
        call ppm_boundary_extrapolation ( grid0, ppoly, u0 )
      end if    
      call remapping_integration ( grid0, u0, ppoly, grid1, u1, &
                                   INTEGRATION_PPM )
    
    case ( REMAPPING_PPM_IH4 )
    
      call edge_values_implicit_h4 ( grid0, u0, ppoly%E )
      call ppm_reconstruction ( grid0, ppoly, u0 )
      if ( regridding_opts%boundary_extrapolation) then
        call ppm_boundary_extrapolation ( grid0, ppoly, u0 )
      end if    
      call remapping_integration ( grid0, u0, ppoly, grid1, u1, &
                                   INTEGRATION_PPM )
      
    case ( REMAPPING_PQM_IH4IH3 )
      
      call edge_values_implicit_h4 ( grid0, u0, ppoly%E )
      call edge_slopes_implicit_h3 ( grid0, u0, ppoly%S )
      call pqm_reconstruction ( grid0, ppoly, u0 )
      if ( regridding_opts%boundary_extrapolation) then
        call pqm_boundary_extrapolation_v1 ( grid0, ppoly, u0 )
      end if    
      call remapping_integration ( grid0, u0, ppoly, grid1, u1, &
                                   INTEGRATION_PQM )
    
    case ( REMAPPING_PQM_IH6IH5 )
      
      call edge_values_implicit_h6 ( grid0, u0, ppoly%E )
      call edge_slopes_implicit_h5 ( grid0, u0, ppoly%S )
      call pqm_reconstruction ( grid0, ppoly, u0 )
      if ( regridding_opts%boundary_extrapolation) then
        call pqm_boundary_extrapolation_v1 ( grid0, ppoly, u0 )
      end if    
      call remapping_integration ( grid0, u0, ppoly, grid1, u1, &
                                   INTEGRATION_PQM )

    case default

      call MOM_error ( FATAL, 'The selected remapping method is invalid' )
      
  end select

end subroutine remapping_core


! -----------------------------------------------------------------------------
! remapping_integration (integration of reconstructed profile)
! -----------------------------------------------------------------------------
subroutine remapping_integration ( grid0, u0, ppoly0, grid1, u1, method )
  ! Arguments
  type(grid1d_t), intent(in)        :: grid0;   ! source grid
  real, dimension(:), intent(in)    :: u0;      ! source cell averages
  type(ppoly_t), intent(in)         :: ppoly0;  ! source piecewise polynomial
  type(grid1d_t), intent(in)        :: grid1;   ! target grid
  real, dimension(:), intent(inout) :: u1;      ! target cell averages
  integer                           :: method;  ! remapping scheme to use
  
  ! Local variables
  integer       :: i, j, k
  integer       :: n0, n1;      ! nb of cells in grid0 and grid1, respectively
  integer       :: j0, j1;      ! indexes of source cells containing target 
                                ! cell edges
  real          :: x0, x1;      ! coordinates of target cell edges  
  real          :: q0, q1;      ! partially integrated quantities in source 
                                ! cells j0 and j1
  real          :: q;           ! complete integration
  real          :: a, b;        ! interval of integration (global coordinates)
  real          :: xi0, xi1;    ! interval of integration (local -- normalized 
                                ! -- coordinates)

  ! A priori, both grids contains the same number of cells but, who knows...
  n0 = grid0%nb_cells
  n1 = grid1%nb_cells

  ! Loop on cells in target grid (grid1). For each target cell, we need to find
  ! in which source cells the target cell edges lie. The associated indexes are 
  ! noted j0 and j1.
  do i = 1,grid1%nb_cells
    ! Determine the coordinates of the target cell edges
    x0 = grid1%x(i)
    x1 = grid1%x(i+1)

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
          exit; ! once target grid cell is found, exit loop
        end if
      end do

      ! If, at this point, j0 is equal to -1, it means the vanished
      ! cell lies outside the source grid. In other words, it means that
      ! the source and target grids do not cover the same physical domain
      ! and there is something very wrong !
      if ( j0 .EQ. -1 ) then
        write(*,*) i 
        call MOM_error ( FATAL, 'The location of the vanished cell could &
                                  not be found in "remapping_integration"' )
      end if

      ! We check whether the source cell (i.e. the cell in which the
      ! vanished target cell lies) is vanished. If it is, the interpolated 
      ! value is set to be mean of the edge values (which should be the same).
      ! If it isn't, we simply interpolate.
      if ( grid0%h(j0) .EQ. 0.0 ) then
        u1(i) = 0.5 * ( ppoly0%E(j0,1) + ppoly0%E(j0,2) )
      else
        xi0 = x0 / grid0%h(j0) - grid0%x(j0) / grid0%h(j0)
    
        select case ( method )
        
          case ( INTEGRATION_PCM )   
            u1(i) = ppoly0%coefficients(j0,1)
            
          case ( INTEGRATION_PLM )  
            
            u1(i) = evaluation_polynomial ( ppoly0%coefficients(j0,:), 2, xi0 )
        
          case ( INTEGRATION_PPM )
            
            u1(i) = evaluation_polynomial ( ppoly0%coefficients(j0,:), 3, xi0 )
        
          case ( INTEGRATION_PQM )
            
            u1(i) = evaluation_polynomial ( ppoly0%coefficients(j0,:), 5, xi0 )
            
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
            exit;   ! once target grid cell is found, exit loop
        end if
        
    end do  

    do j = 1,grid0%nb_cells
        
        ! Right edge is found in cell j
        if ( ( x1 .GE. grid0%x(j) ) .AND. ( x1 .LE. grid0%x(j+1) ) ) then
            j1 = j
            exit;   ! once target grid cell is found, exit loop
        end if
        
    end do ! end loop on source grid cells

    ! Here, we make sure that the boundary edges of boundary cells
    ! coincide
    if ( i .EQ. 1 ) then
      j0 = 1
    end if  
    
    if ( i .EQ. grid1%nb_cells ) then
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
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 1 )
        
        case ( INTEGRATION_PPM )
          q = grid0%h(j0) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 2 )
    
        case ( INTEGRATION_PQM )
          q = grid0%h(j0) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 4 )

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
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 1 )
    
        case ( INTEGRATION_PPM )
          q = q + grid0%h(j0) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 2 )
    
        case ( INTEGRATION_PQM )
          q = q + grid0%h(j0) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 4 )

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
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j1,:), 1 )
    
        case ( INTEGRATION_PPM )
          q = q + grid0%h(j1) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j1,:), 2 )
    
        case ( INTEGRATION_PQM )
          q = q + grid0%h(j1) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j1,:), 4 )

        case default
          call MOM_error( FATAL,'The selected integration method is invalid' )
      end select     
      
    end if ! end integration for non-vanished cells 
    
    ! The cell average is the integrated value divided by the cell width
    u1(i) = q / grid1%h(i)
    
    end if ! end if clause to check if cell is vanished
    
  end do ! end i loop on target grid cells

end subroutine remapping_integration


!------------------------------------------------------------------------------
! Memory allocation for remapping
!------------------------------------------------------------------------------
subroutine remapping_memory_allocation ( G, regridding_opts )

  ! Arguments
  type(ocean_grid_type), intent(in)   :: G
  type(regridding_opts_t), intent(in) :: regridding_opts
  
  ! Local variables
  integer   :: nz
  integer   :: degree;      ! Degree of polynomials used for the reconstruction 
  
  nz = G%ke
  
  ! Allocate memory for grids
  call grid1d_init ( grid_start, nz )
  call grid1d_init ( grid_final, nz )
  
  ! Piecewise polynomials used for remapping
  select case ( regridding_opts%remapping_scheme )
    case ( REMAPPING_PCM )   
      degree = 0
    case ( REMAPPING_PLM )   
      degree = 1
    case ( REMAPPING_PPM_H4, REMAPPING_PPM_IH4 )  
      degree = 2
    case ( REMAPPING_PQM_IH4IH3, REMAPPING_PQM_IH6IH5 )
      degree = 4
    case default
      call MOM_error ( FATAL, 'The selected remapping scheme is invalid' )
  end select
  
  call ppoly_init ( ppoly_r, nz, degree )
  
  allocate ( u_column(nz) ); u_column = 0.0

end subroutine remapping_memory_allocation


!------------------------------------------------------------------------------
! Memory deallocation for remapping
!------------------------------------------------------------------------------
subroutine remapping_memory_deallocation ( )

  ! Deallocate memory for grid
  call grid1d_destroy ( grid_start )
  call grid1d_destroy ( grid_final )
  
  ! Piecewise polynomials
  call ppoly_destroy ( ppoly_r )

  deallocate ( u_column )

end subroutine remapping_memory_deallocation


end module MOM_remapping
