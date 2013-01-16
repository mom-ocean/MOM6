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
use MOM_file_parser,   only : get_param, param_file_type
use MOM_EOS,           only : calculate_density

use regrid_grid1d_class ! see 'regrid_grid.F90'
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
use MOM_remapping      ! see 'MOM_remapping.F90'
use regrid_defs         ! see 'regrid_defs.F90' (contains types and parameters)

implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Private (module-wise) variables
! -----------------------------------------------------------------------------

! Generic grid used for various purposes throughout the code (the same 
! grid is used to avoid having dynamical memory allocation)
type(grid1d_t)                  :: grid_generic

! Generic linear piecewise polynomial used for various purposes throughout 
! the code (the same ppoly is used to avoid having dynamical memory allocation)
type(ppoly_t)                   :: ppoly_linear

! Generic parabolic piecewise polynomial used for various purposes 
! throughout the code (the same ppoly is used to avoid having dynamical 
! memory allocation)
type(ppoly_t)                   :: ppoly_parab

type(grid1d_t)                  :: grid_start      ! starting grid
type(grid1d_t)                  :: grid_trans      ! transition/iterated grid
type(grid1d_t)                  :: grid_final      ! final grid
type(ppoly_t)                   :: ppoly_i         ! interpolation ppoly
type(ppoly_t)                   :: ppoly_r         ! reconstruction ppoly
                                                    
! These variables are private to this module and memory is allocated at
! the beginning of the simulation. Note that 'u_column' is a generic variable
! used for remapping. It is not necessarily associated with the zonal velocity
! component.
real, dimension(:), allocatable :: target_values
real, dimension(:), allocatable :: densities
real, dimension(:), allocatable :: mapping
real, dimension(:), allocatable :: u_column
real, dimension(:), allocatable :: T_column
real, dimension(:), allocatable :: S_column
real, dimension(:), allocatable :: p_column

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public initialize_regridding
public end_regridding
public regridding_main 
public pressure_gradient_plm
public pressure_gradient_ppm

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! Initialization of regridding
!------------------------------------------------------------------------------
subroutine initialize_regridding ( param_file, regridding_opts, G, h, h_aux, &
                                   u, v, tv )
!------------------------------------------------------------------------------
! This routine is typically called (from initialize_MOM in file MOM.F90)
! before the main time integration loop to initialize the regridding stuff.
! We read the MOM_input file to register the values of different
! regridding/remapping parameters.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(param_file_type), intent(in)                      :: param_file
  type(regridding_opts_t), intent(inout)                 :: regridding_opts
  type(ocean_grid_type), intent(in)                      :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_,C2_), intent(inout) :: h
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout)     :: h_aux
  real, dimension(NIMEMB_,NJMEM_, NKMEM_), intent(inout)    :: u
  real, dimension(NIMEM_,NJMEMB_, NKMEM_), intent(inout)    :: v
  type(thermo_var_ptrs), intent(inout)                   :: tv

  ! Local variables
  logical       :: verbose
  integer       :: k, m
  integer       :: i, j

  verbose = .false.
   
  ! Read regridding/remapping parameters from MOM_input file
  call initialize_regridding_options ( param_file, regridding_opts )
  
  ! If regridding must be used, the following steps are performed
  if ( regridding_opts%use_regridding ) then 
    
    ! Memory allocation for regridding
    call regridding_memory_allocation ( G, regridding_opts )

    ! Check grid integrity with respect to minimum allowed thickness
    do m = 1,size(h,4)
      call check_grid_integrity ( G, h(:,:,:,m), regridding_opts )
    end do  

    if ( verbose ) then
        i = 20
        j = 4
        write(*,*) 'Before regridding/remapping', i, j
        do k = 1,G%ke
          write(*,*) h(i,j,k,1), tv%T(i,j,k)
        end do
        i = 22
        j = 4
        write(*,*) 'Before regridding/remapping', i, j
        do k = 1,G%ke
          write(*,*) h(i,j,k,1), tv%T(i,j,k)
        end do
    end if

    ! Perform one regridding/remapping step -- This should NOT modify
    ! neither the initial grid nor the initial cell averages. This
    ! step is therefore not strictly necessary but is included for historical
    ! reasons when I needed to check whether the combination 'initial 
    ! conditions - regridding/remapping' was consistently implemented.
    call regridding_main ( G, h(:,:,:,1), h_aux, u, v, tv, regridding_opts )
    
    if ( verbose ) then
        i = 20
        j = 4
        write(*,*) 'After regridding/remapping', i, j
        do k = 1,G%ke
          write(*,*) h(i,j,k,1), tv%T(i,j,k)
        end do
        i = 22
        j = 4
        write(*,*) 'After regridding/remapping', i, j
        do k = 1,G%ke
          write(*,*) h(i,j,k,1), tv%T(i,j,k)
        end do
    end if
  
  end if ! use_regridding

end subroutine initialize_regridding


!------------------------------------------------------------------------------
! Initialization of regridding options
!------------------------------------------------------------------------------
subroutine initialize_regridding_options ( param_file, regridding_opts )
!------------------------------------------------------------------------------
! Read the regridding/remapping parameters in the MOM_input file and
! update the structure that is passed as argument all over the place.
!------------------------------------------------------------------------------

  ! Arguments
  type(param_file_type), intent(in)        :: param_file
  type(regridding_opts_t), intent(inout)   :: regridding_opts
  ! Local variables
  character(len=40)  :: mod = "MOM_regridding" ! This module's name.

  ! This sets whether we want to use regridding or not. By default, 
  ! regridding is NOT used but this can be overridden in the input
  ! file
  call get_param(param_file, mod, "USE_REGRIDDING", &
                 regridding_opts%use_regridding , &
                 "If True, use the ALE algorithm (regridding/remapping).\n"//&
                 "If False, use the layered isopycnal algorithm.", default=.false. )

  ! When not using ALE, we always use the layer-wise constant integration of T/S
  ! in the FV pressure gradient integrals.
  regridding_opts%reconstructForPressure = .false.

  if (regridding_opts%use_regridding) then
  ! The following options are only relevant when the property 'use_regridding'
  ! is true.

  ! --- TYPE OF VERTICAL GRID ---
  ! This sets which kind of grid we want to use in the vertical. If none
  ! is specified, target interface densities are used to build the grid
  call get_param(param_file, mod, "REGRIDDING_SCHEME", &
                 regridding_opts%regridding_scheme, &
                 "Type of grid to build in the vertical."//&
                 "Choose among the following possibilities (must\n"//&
                 "be an integer !):\n"//&
                 " 0: z*.\n"//&
                 " 1: target interface densities.\n"//&
                 " 2: sigma.\n", fail_if_missing=.true.)

  ! --- REMAPPING SCHEME ---
  ! This sets which remapping scheme we want to use to remap all variables
  ! betwenn grids. If none is specified, PLM is used for remapping.
  call get_param(param_file, mod, "REMAPPING_SCHEME", &
                 regridding_opts%remapping_scheme, &
                 "This sets the remapping scheme to use to\n"//&
                 "remap all variables between successive grids.\n"//&
                 "It can be one of the following schemes (must be\n"//&
                 "an integer !):\n"//&
                 "0: PCM         (1st-order accurate)\n"//&
                 "1: PLM         (2nd-order accurate)\n"//&
                 "2: PPM_H4         (3rd-order accurate)\n"//&
                 "3: PPM_IH4     (3rd-order accurate)\n"//&
                 "4: PQM_IH4IH3     (4th-order accurate)\n"//&
                 "5: PQM_IH6IH5     (5th-order accurate)\n", &
                 default=REMAPPING_PLM)
  
  ! --- INTERPOLATION SCHEME ---
  ! This sets which interpolation scheme we want to use to define the new
  ! grid when regridding is based upon target interface densities. If none
  ! is specified, the p1m h2 interpolation scheme is used.
  call get_param(param_file, mod, "INTERPOLATION_SCHEME", &
                 regridding_opts%interpolation_scheme, &
                 "This sets the interpolation scheme to use to\n"//&
                 "determine the new grid. These parameters are\n"//&
                 "only relevant when REGRIDDING_SCHEME is set to\n"//&
                 "1. Otherwise, it is not used.\n"//&
                 "It can be one of the following schemes (must be\n"//&
                 "an integer !):\n"//&
                 "0: P1M_H2         (2nd-order accurate)\n"//&
                 "1: P1M_H4         (2nd-order accurate)\n"//&
                 "2: P1M_IH4     (2nd-order accurate)\n"//&
                 "3: PLM         (2nd-order accurate)\n"//&
                 "4: PPM_H4        (3rd-order accurate)\n"//&
                 "5: PPM_IH4        (3rd-order accurate)\n"//&
                 "6: P3M_IH4IH3     (4th-order accurate)\n"//&
                 "7: P3M_IH6IH5     (4th-order accurate)\n"//&
                 "8: PQM_IH4IH3    (4th-order accurate)\n"//&
                 "9: PQM_IH6IH5    (5th-order accurate)\n", &
                 default=INTERPOLATION_P1M_H2)
  
  ! --- BOUNDARY EXTRAPOLATION --
  ! This sets whether high-order (rather than PCM) reconstruction schemes
  ! should be used within boundary cells
  call get_param(param_file, mod, "BOUNDARY_EXTRAPOLATION", &
                 regridding_opts%boundary_extrapolation, &
                 "When defined, a proper high-order reconstruction\n"//&
                 "scheme is used within boundary cells rather\n"//&
                 "than PCM. E.g., if PPM is used for remapping, a\n"//&
                 "PPM reconstruction will also be used within\n"//&
                 "boundary cells.", default=.false.)

  ! --- PRESSURE GRADIENT CALCULATION ---
  call get_param(param_file, mod, "RECONSTRUCT_FOR_PRESSURE", &
                 regridding_opts%reconstructForPressure , &
                 "If True, use vertical reconstruction of T/S within\n"//&
                 "the integrals of teh FV pressure gradient calculation.\n"//&
                 "If False, use the constant-by-layer algorithm.\n"//&
                 "By default, this is True when using ALE and False otherwise.", &
                 default=.true. )

  call get_param(param_file, mod, "PRESSURE_RECONSTRUCTION_SCHEME", &
                 regridding_opts%pressureReconstructionScheme, &
                 "Type of vertical reconstruction of T/S to use in integrals\n"//&
                 "within the FV pressure gradient calculation."//&
                 " 1: PLM reconstruction.\n"//&
                 " 2: PPM reconstruction.", default=PRESSURE_RECONSTRUCTION_PLM)

  ! --- MINIMUM THICKNESS ---
  call get_param(param_file, mod, "MIN_THICKNESS", &
                 regridding_opts%min_thickness, &
                 "When regridding, this is the minimum layer\n"//&
                 "thickness allowed.", default=1.e-3 )

  endif
  
end subroutine initialize_regridding_options



!------------------------------------------------------------------------------
! End of regridding (memory deallocation)
!------------------------------------------------------------------------------
subroutine end_regridding ( regridding_opts )
!------------------------------------------------------------------------------
! This routine is typically called (from MOM_end in file MOM.F90)
! after the main time integration loop to deallocate the regridding stuff.
!------------------------------------------------------------------------------
  type(regridding_opts_t), intent(in)   :: regridding_opts
  
  ! Deallocate memory used for the regridding
  if ( regridding_opts%use_regridding ) then
       call regridding_memory_deallocation ( regridding_opts )
  end if

end subroutine end_regridding


!------------------------------------------------------------------------------
! Dispatching regridding routine: regridding & remapping
!------------------------------------------------------------------------------
subroutine regridding_main ( G, h, h_new, u, v, tv, regridding_opts )
!------------------------------------------------------------------------------
! This routine takes care of (1) building a new grid and (2) remapping between
! the old grid and the new grid. The creation of the new grid can be based
! on z coordinates, target interface densities, sigma coordinates or any
! arbitrary coordinate system.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(ocean_grid_type), intent(in)                    :: &
  G;     ! Ocean grid informations
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout)   :: &
  h;     ! Current 3D grid obtained after the last time step
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout)   :: &
  h_new; ! The new 3D grid obtained via regridding
  real, dimension(NIMEMB_,NJMEM_, NKMEM_), intent(inout)  :: &
  u;     ! Zonal velocity field
  real, dimension(NIMEM_,NJMEMB_, NKMEM_), intent(inout)  :: &
  v;     ! Meridional velocity field
  type(thermo_var_ptrs), intent(inout)                 :: &
  tv;    ! Thermodynamical variables (T, S, ...)  
  type(regridding_opts_t), intent(in)                  :: &
  regridding_opts; ! Regridding parameters and options

  ! Local variables
  integer   :: i, j, k
  
  ! Build new grid. The new grid is stored in h_new. The old grid is h.
  ! Both are needed for the subsequent remapping of variables.
  select case ( regridding_opts%regridding_scheme )

    case ( REGRIDDING_Z )
      call build_grid_uniform ( G, h, h_new, regridding_opts )

    case ( REGRIDDING_RHO )  
      call convective_adjustment ( G, h, tv )
      call build_grid_target_densities ( G, h, h_new, tv, regridding_opts )

    case ( REGRIDDING_SIGMA )
      call build_grid_sigma ( G, h, h_new )
    
    case ( REGRIDDING_ARBITRARY )
      call build_grid_arbitrary ( G, h, h_new, regridding_opts )

  end select ! type of grid 
  
  ! Remap all variables from old grid h onto new grid h_new
  call remapping_main ( G, regridding_opts, h, h_new, tv, u, v )
  
  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
  h(:,:,:) = h_new(:,:,:)

end subroutine regridding_main


!------------------------------------------------------------------------------
! Use plm reconstruction for pressure gradient (determine edge values)
!------------------------------------------------------------------------------
subroutine pressure_gradient_plm ( S_t, S_b, T_t, T_b, G, tv, h )
!------------------------------------------------------------------------------
! By using a PLM (limited piecewise linear method) reconstruction, this 
! routine determines the edge values for the salinity and temperature 
! within each layer. These edge values are returned and are used to compute 
! the pressure gradient (by computing the densities).
!------------------------------------------------------------------------------

  ! Arguments
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: &
  S_t, S_b; ! Salinity at the top and bottom edges of each layer
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: &
  T_t, T_b; ! Temperature at the top and bottom edges of each layer
  type(ocean_grid_type), intent(in)                 :: G
  type(thermo_var_ptrs), intent(in)                 :: tv
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: &
  h;        ! Three-dimensional ocean grid

  ! Local variables
  integer           :: i, j, k

  ! NOTE: the variables 'grid_generic' and 'ppoly_linear' are declared at
  ! the module level. Memory is allocated once at the beginning of the run
  ! in 'regridding_memory_allocation'.

  ! Determine reconstruction within each column
  do i = G%isc,G%iec+1
    do j = G%jsc,G%jec+1
     
      ! Build current grid
      grid_generic%h(:) = h(i,j,:)
      grid_generic%x(1) = 0.0
      do k = 1,G%ke
        grid_generic%x(k+1) = grid_generic%x(k) + grid_generic%h(k)
      end do
      
      ! Reconstruct salinity profile    
      ppoly_linear%E = 0.0
      ppoly_linear%coefficients = 0.0
      call plm_reconstruction ( grid_generic, ppoly_linear, tv%S(i,j,:) )
      call plm_boundary_extrapolation ( grid_generic, ppoly_linear, &
                                        tv%S(i,j,:) )
      
      do k = 1,G%ke
        S_t(i,j,k) = ppoly_linear%E(k,1)
        S_b(i,j,k) = ppoly_linear%E(k,2)
      end do
      
      ! Reconstruct temperature profile 
      ppoly_linear%E = 0.0
      ppoly_linear%coefficients = 0.0
      call plm_reconstruction ( grid_generic, ppoly_linear, tv%T(i,j,:) )
      call plm_boundary_extrapolation ( grid_generic, ppoly_linear, &
                                        tv%T(i,j,:) )
      
      do k = 1,G%ke
        T_t(i,j,k) = ppoly_linear%E(k,1)
        T_b(i,j,k) = ppoly_linear%E(k,2)
      end do
      
    end do
  end do

end subroutine pressure_gradient_plm


!------------------------------------------------------------------------------
! Use ppm reconstruction for pressure gradient (determine edge values)
!------------------------------------------------------------------------------
subroutine pressure_gradient_ppm ( S_t, S_b, T_t, T_b, G, tv, h )
!------------------------------------------------------------------------------
! By using a PPM (limited piecewise linear method) reconstruction, this 
! routine determines the edge values for the salinity and temperature 
! within each layer. These edge values are returned and are used to compute 
! the pressure gradient (by computing the densities).
!------------------------------------------------------------------------------

  ! Arguments
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: &
  S_t, S_b; ! Salinity at the top and bottom edges of each layer
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: &
  T_t, T_b; ! Temperature at the top and bottom edges of each layer
  type(ocean_grid_type), intent(in)                 :: G
  type(thermo_var_ptrs), intent(in)                 :: tv
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: &
  h;        ! Three-dimensional ocean grid

  ! Local variables
  integer           :: i, j, k

  ! NOTE: the variables 'grid_generic' and 'ppoly_parab' are declared at
  ! the module level. Memory is allocated once at the beginning of the run
  ! in 'regridding_memory_allocation'.

  ! Determine reconstruction within each column
  do i = G%isc,G%iec+1
    do j = G%jsc,G%jec+1
     
      ! Build current grid
      grid_generic%h(:) = h(i,j,:)
      grid_generic%x(1) = 0.0
      do k = 1,G%ke
        grid_generic%x(k+1) = grid_generic%x(k) + grid_generic%h(k)
      end do
      
      ! Reconstruct salinity profile    
      ppoly_parab%E = 0.0
      ppoly_parab%coefficients = 0.0
      call edge_values_implicit_h4 ( grid_generic, tv%S(i,j,:), ppoly_parab%E )
      call ppm_reconstruction ( grid_generic, ppoly_parab, tv%S(i,j,:) )
      call ppm_boundary_extrapolation ( grid_generic, ppoly_parab, &
                                        tv%S(i,j,:) )
      
      do k = 1,G%ke
        S_t(i,j,k) = ppoly_parab%E(k,1)
        S_b(i,j,k) = ppoly_parab%E(k,2)
      end do
      
      ! Reconstruct temperature profile 
      ppoly_parab%E = 0.0
      ppoly_parab%coefficients = 0.0
      call edge_values_implicit_h4 ( grid_generic, tv%T(i,j,:), ppoly_parab%E )
      call ppm_reconstruction ( grid_generic, ppoly_parab, tv%T(i,j,:) )
      call ppm_boundary_extrapolation ( grid_generic, ppoly_parab, &
                                        tv%T(i,j,:) )
      
      do k = 1,G%ke
        T_t(i,j,k) = ppoly_parab%E(k,1)
        T_b(i,j,k) = ppoly_parab%E(k,2)
      end do
      
    end do
  end do

end subroutine pressure_gradient_ppm


!------------------------------------------------------------------------------
! Build uniform z*-ccordinate grid with partial steps
!------------------------------------------------------------------------------
subroutine build_grid_uniform ( G, h, h_new, regridding_opts )
!------------------------------------------------------------------------------
! This routine builds a grid where the distribution of levels is based on a
! z* coordinate system with partial steps. Within each water column, a uniform
! layer thickness is determined based on the local free-surface elevation and
! the global maximum depth. Layers are then distributed vertically, modulated by
! topography.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(ocean_grid_type), intent(in)                  :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout) :: h_new
  type(regridding_opts_t), intent(in)                :: regridding_opts
  
  ! Local variables
  integer   :: i, j, k
  integer   :: nz
  real      :: z_inter(SZK_(G)+1)
  real      :: total_height
  real      :: delta_h
  real      :: max_depth
  real      :: min_thickness
  real      :: eta;             ! local elevation
  real      :: local_depth

  nz = G%ke
  
  max_depth = G%max_depth
  min_thickness = regridding_opts%min_thickness

  do j = G%jsc,G%jec+1
    do i = G%isc,G%iec+1

      ! Local depth
      local_depth = G%bathyT(i,j)
      
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
    
      ! Modify interface heights to account for topography
      z_inter(nz+1) = - local_depth

      ! Modify interface heights to avoid layers of zero thicknesses
      do k = nz,1,-1
        if ( z_inter(k) .LT. (z_inter(k+1) + min_thickness) ) then
          z_inter(k) = z_inter(k+1) + min_thickness
        end if
      end do
    
      ! Define thicknesses in terms of interface heights
      do k = 1,nz
        h_new(i,j,k) = z_inter(k) - z_inter(k+1)
      end do    

    end do
  end do

end subroutine build_grid_uniform


!------------------------------------------------------------------------------
! Build sigma grid
!------------------------------------------------------------------------------
subroutine build_grid_sigma ( G, h, h_new )
!------------------------------------------------------------------------------
! This routine builds a grid based on terrain-following coordinates.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(ocean_grid_type), intent(in)                  :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout) :: h_new
  
  ! Local variables
  integer   :: i, j, k
  integer   :: nz
  real      :: total_height
  real      :: delta_h

  nz = G%ke
  
  do i = G%isc,G%iec+1
    do j = G%jsc,G%jec+1
      
      ! Determine water column height
      total_height = 0.0
      do k = 1,nz
        total_height = total_height + h(i,j,k)
      end do
          
      ! Compute new thicknesses based on stretched water column
      delta_h = total_height / nz
      
      ! Define thicknesses in terms of interface heights
      do k = 1,nz
        h_new(i,j,k) = delta_h
      end do    
      
    end do
  end do
  
end subroutine build_grid_sigma


!------------------------------------------------------------------------------
! Build arbitrary grid
!------------------------------------------------------------------------------
subroutine build_grid_arbitrary ( G, h, h_new, regridding_opts )
!------------------------------------------------------------------------------
! This routine builds a grid based on arbitrary rules
!------------------------------------------------------------------------------
  
  ! Arguments
  type(ocean_grid_type), intent(in)                  :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout) :: h_new
  type(regridding_opts_t), intent(in)                :: regridding_opts
  
  ! Local variables
  integer   :: i, j, k
  integer   :: nz
  real      :: z_inter(SZK_(G)+1)
  real      :: total_height
  real      :: delta_h
  real      :: max_depth
  real      :: min_thickness
  real      :: eta;             ! local elevation
  real      :: local_depth
  real      :: x1, y1, x2, y2
  real      :: x, t

  nz = G%ke
  
  max_depth = G%max_depth
  min_thickness = regridding_opts%min_thickness

  do j = G%jsc,G%jec+1
    do i = G%isc,G%iec+1

      ! Local depth
      local_depth = G%bathyT(i,j)
      
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
    
      ! Define thicknesses in terms of interface heights
      do k = 1,nz
        h_new(i,j,k) = z_inter(k) - z_inter(k+1)
      end do    

    end do
  end do
  
  
end subroutine build_grid_arbitrary


!------------------------------------------------------------------------------
! Build grid based on target interface densities
!------------------------------------------------------------------------------
subroutine build_grid_target_densities ( G, h, h_new, tv, regridding_opts )
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
  type(ocean_grid_type), intent(in)                  :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(in)    :: h
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout) :: h_new
  type(thermo_var_ptrs), intent(in)                  :: tv     
  type(regridding_opts_t), intent(in)                :: regridding_opts
  
  ! Local variables
  integer   :: i, j, k, m
  integer   :: map_index
  integer   :: nz
  integer   :: k_found
  integer   :: count_nonzero_layers
  real      :: deviation;           ! When iterating to determine the final 
                                    ! grid, this is the deviation between two
                                    ! successive grids.
  real      :: threshold
  real      :: max_thickness
  real      :: correction
  
  real      :: x1, y1, x2, y2, x
  real      :: t
                                    
  nz = G%ke
  threshold = regridding_opts%min_thickness
  
  ! Prescribe target values
  target_values(1)    = G%Rlay(1)+0.5*(G%Rlay(1)-G%Rlay(2))
  target_values(nz+1) = G%Rlay(nz)+0.5*(G%Rlay(nz)-G%Rlay(nz-1))

  do k = 2,nz
    target_values(k) = target_values(k-1) + ( G%Rlay(nz) - G%Rlay(1) ) / (nz-1)
  end do
  
  ! Build grid based on target interface densities
  do i = G%isc,G%iec+1
    do j = G%jsc,G%jec+1
    
      ! Copy T and S onto new variables so as to not alter the original values
      ! of T and S (these are remapped at the end of the regridding iterations
      ! once the final grid has been determined).
      T_column = tv%T(i,j,:)
      S_column = tv%S(i,j,:)
        
      ! Copy original grid
      grid_start%h(1:nz) = h(i,j,1:nz)
      grid_start%x(1) = 0.0
      do k = 1,nz
        grid_start%x(k+1) = grid_start%x(k) + grid_start%h(k)
      end do
    
      ! Start iterations to build grid
      m = 1
      deviation = 1e10
      do while ( ( m .le. NB_REGRIDDING_ITERATIONS ) .and. &
                 ( deviation .gt. DEVIATION_TOLERANCE ) )

        ! Count number of nonzero layers within current water column
        count_nonzero_layers = 0
        do k = 1,nz
          if ( grid_start%h(k) .gt. threshold ) then
            count_nonzero_layers = count_nonzero_layers + 1
          end if
        end do

        ! If there is at most one nonzero layer, stop here (no regridding)
        if ( count_nonzero_layers .le. 1 ) then 
          grid_final%h(1:nz) = grid_start%h(1:nz)
          exit  ! stop iterations here
        end if

        ! Limit number of usable cells
        grid_trans%nb_cells = count_nonzero_layers

        ! Build new grid containing only nonzero layers
        map_index = 1
        correction = 0.0
        do k = 1,nz
          if ( grid_start%h(k) .gt. threshold ) then
            mapping(map_index) = k
            grid_trans%h(map_index) = grid_start%h(k)
            map_index = map_index + 1
          else
            correction = correction + grid_start%h(k)
          end if
        end do

        max_thickness = grid_trans%h(1)
        k_found = 1
        do k = 1,grid_trans%nb_cells
          if ( grid_trans%h(k) .gt. max_thickness ) then
            max_thickness = grid_trans%h(k)
            k_found = k
          end if  
        end do

        grid_trans%h(k_found) = grid_trans%h(k_found) + correction

        grid_trans%x(1) = 0.0
        do k = 1,grid_trans%nb_cells
          grid_trans%x(k+1) = grid_trans%x(k) + grid_trans%h(k)
        end do


        ! Compute densities within current water column
        call calculate_density ( T_column, S_column, p_column, densities,&
                                 1, nz, tv%eqn_of_state )

        do k = 1,count_nonzero_layers
          densities(k) = densities(mapping(k))
        end do

        ! One regridding iteration
        call regridding_iteration ( densities, target_values, regridding_opts,&
                                    grid_trans, ppoly_i, grid_final )
        
        ! Remap T and S from previous grid to new grid
        do k = 1,nz
          grid_start%h(k) = grid_start%x(k+1) - grid_start%x(k)
          grid_final%h(k) = grid_final%x(k+1) - grid_final%x(k)
        end do
        
        call remapping_core ( grid_start, S_column, grid_final,& 
                              S_column, ppoly_r, regridding_opts )
        
        call remapping_core ( grid_start, T_column, grid_final,& 
                              T_column, ppoly_r, regridding_opts )

        ! Compute the deviation between two successive grids
        deviation = 0.0
        do k = 2,nz
          deviation = deviation + (grid_start%x(k)-grid_final%x(k))**2
        end do
        deviation = sqrt ( deviation / (nz-1) )

    
        m = m + 1
        
        ! Copy final grid onto start grid for next iteration
        grid_start%x = grid_final%x
        grid_start%h = grid_final%h

      end do ! end regridding iterations               

      ! The new grid is that obtained after the iterations
      do k = 1,nz
        h_new(i,j,k) = grid_final%h(k)
      end do
        
    end do; ! end loop on j 
  end do; ! end loop on i
      
end subroutine build_grid_target_densities


!------------------------------------------------------------------------------
! Regridding iterations
!------------------------------------------------------------------------------
subroutine regridding_iteration ( densities, target_values, regridding_opts, &
                                  grid0, ppoly0, grid1 )
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
  real, dimension(:), intent(in)      :: &
  densities;        ! Actual cell densities
  real, dimension(:), intent(in)      :: &
  target_values;    ! Target interface densities
  type(grid1d_t), intent(in)          :: &
  grid0             ! The grid on which cell densities are known                
  type(ppoly_t), intent(inout)        :: &
  ppoly0;           ! Piecewise polynomial for density interpolation
  type(grid1d_t), intent(inout)       :: &
  grid1;            ! The new grid based on target interface densities          
  type(regridding_opts_t), intent(in) :: &
  regridding_opts;  ! Parameters used for regridding

  ! Local variables
  integer        :: degree

  ! Reset piecewise polynomials
  ppoly0%E = 0.0
  ppoly0%S = 0.0
  ppoly0%coefficients = 0.0
  
  ! 1. Compute the interpolated profile of the density field and build grid
  select case ( regridding_opts%interpolation_scheme )
  
    case ( INTERPOLATION_P1M_H2 )
      degree = DEGREE_1
      call edge_values_explicit_h2 ( grid0, densities, ppoly0%E )
      call p1m_interpolation ( grid0, ppoly0, densities )
      if ( regridding_opts%boundary_extrapolation) then
        call p1m_boundary_extrapolation ( grid0, ppoly0, densities )
      end if    
    
    case ( INTERPOLATION_P1M_H4 )
      degree = DEGREE_1
      if ( grid0%nb_cells .ge. 4 ) then
        call edge_values_explicit_h4 ( grid0, densities, ppoly0%E )
      else
        call edge_values_explicit_h2 ( grid0, densities, ppoly0%E )
      end if
      call p1m_interpolation ( grid0, ppoly0, densities )
      if ( regridding_opts%boundary_extrapolation) then
        call p1m_boundary_extrapolation ( grid0, ppoly0, densities )
      end if    
    
    case ( INTERPOLATION_P1M_IH4 )
      degree = DEGREE_1
      if ( grid0%nb_cells .ge. 4 ) then
        call edge_values_implicit_h4 ( grid0, densities, ppoly0%E )
      else
        call edge_values_explicit_h2 ( grid0, densities, ppoly0%E )
      end if
      call p1m_interpolation ( grid0, ppoly0, densities )
      if ( regridding_opts%boundary_extrapolation) then
        call p1m_boundary_extrapolation ( grid0, ppoly0, densities )
      end if    
    
    case ( INTERPOLATION_PLM )   
      degree = DEGREE_1
      call plm_reconstruction ( grid0, ppoly0, densities )
      if ( regridding_opts%boundary_extrapolation) then
        call plm_boundary_extrapolation ( grid0, ppoly0, densities )
      end if    
    
    case ( INTERPOLATION_PPM_H4 )
      if ( grid0%nb_cells .ge. 4 ) then
        degree = DEGREE_2
        call edge_values_explicit_h4 ( grid0, densities, ppoly0%E )
        call ppm_reconstruction ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call ppm_boundary_extrapolation ( grid0, ppoly0, densities )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2 ( grid0, densities, ppoly0%E )
        call p1m_interpolation ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call p1m_boundary_extrapolation ( grid0, ppoly0, densities )
        end if
      end if
    
    case ( INTERPOLATION_PPM_IH4 )

      if ( grid0%nb_cells .ge. 4 ) then
        degree = DEGREE_2
        call edge_values_implicit_h4 ( grid0, densities, ppoly0%E )
        call ppm_reconstruction ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call ppm_boundary_extrapolation ( grid0, ppoly0, densities )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2 ( grid0, densities, ppoly0%E )
        call p1m_interpolation ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call p1m_boundary_extrapolation ( grid0, ppoly0, densities )
        end if  
      end if
    
    case ( INTERPOLATION_P3M_IH4IH3 )
      
      if ( grid0%nb_cells .ge. 4 ) then
        degree = DEGREE_3
        call edge_values_implicit_h4 ( grid0, densities, ppoly0%E )
        call edge_slopes_implicit_h3 ( grid0, densities, ppoly0%S )
        call p3m_interpolation ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call p3m_boundary_extrapolation ( grid0, ppoly0, densities )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2 ( grid0, densities, ppoly0%E )
        call p1m_interpolation ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call p1m_boundary_extrapolation ( grid0, ppoly0, densities )
        end if  
      end if
      
    case ( INTERPOLATION_P3M_IH6IH5 )
      if ( grid0%nb_cells .ge. 6 ) then
        degree = DEGREE_3
        call edge_values_implicit_h6 ( grid0, densities, ppoly0%E )
        call edge_slopes_implicit_h5 ( grid0, densities, ppoly0%S )
        call p3m_interpolation ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call p3m_boundary_extrapolation ( grid0, ppoly0, densities )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2 ( grid0, densities, ppoly0%E )
        call p1m_interpolation ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call p1m_boundary_extrapolation ( grid0, ppoly0, densities )
        end if  
      end if
    
    case ( INTERPOLATION_PQM_IH4IH3 )
    
      if ( grid0%nb_cells .ge. 4 ) then
        degree = DEGREE_4
        call edge_values_implicit_h4 ( grid0, densities, ppoly0%E )
        call edge_slopes_implicit_h3 ( grid0, densities, ppoly0%S )
        call pqm_reconstruction ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call pqm_boundary_extrapolation_v1 ( grid0, ppoly0, densities )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2 ( grid0, densities, ppoly0%E )
        call p1m_interpolation ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call p1m_boundary_extrapolation ( grid0, ppoly0, densities )
        end if  
      end if
    
    case ( INTERPOLATION_PQM_IH6IH5 )
      if ( grid0%nb_cells .ge. 6 ) then
        degree = DEGREE_4
        call edge_values_implicit_h6 ( grid0, densities, ppoly0%E )
        call edge_slopes_implicit_h5 ( grid0, densities, ppoly0%S )
        call pqm_reconstruction ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call pqm_boundary_extrapolation_v1 ( grid0, ppoly0, densities )
        end if  
      else
        degree = DEGREE_1
        call edge_values_explicit_h2 ( grid0, densities, ppoly0%E )
        call p1m_interpolation ( grid0, ppoly0, densities )
        if ( regridding_opts%boundary_extrapolation) then
          call p1m_boundary_extrapolation ( grid0, ppoly0, densities )
        end if  
      end if
    
  end select
  
  ! Based on global density profile, interpolate new grid and 
  ! inflate vanished layers    
  call interpolate_grid ( grid0, ppoly0, grid1, target_values, degree )
  call inflate_vanished_layers ( grid1, regridding_opts )
    
end subroutine regridding_iteration


!------------------------------------------------------------------------------
! Given target values (e.g., density), build new grid based on polynomial 
!------------------------------------------------------------------------------
subroutine interpolate_grid ( grid0, ppoly0, grid1, target_values, degree )
! ------------------------------------------------------------------------------
! Given the grid 'grid0' and the piecewise polynomial interpolant 
! 'ppoly0' (possibly discontinuous), the coordinates of the new grid 'grid1' 
! are determined by finding the corresponding target interface densities.
! ------------------------------------------------------------------------------
  
  ! Arguments
  type(grid1d_t), intent(in)      :: grid0              
  type(ppoly_t), intent(in)       :: ppoly0
  type(grid1d_t), intent(inout)   :: grid1              
  real, dimension(:), intent(in)  :: target_values
  integer, intent(in)             :: degree
    
  ! Local variables
  integer        :: k;  ! loop index
  integer        :: m
  integer        :: N;  ! number of grid cells
  integer        :: nz; ! number of layers
  real           :: t;  ! current interface target density
  integer        :: count_nonzero_layers
  real           :: delta_h
  
  N = grid0%nb_cells
  nz = grid1%nb_cells

  ! Make sure boundary coordinates of new grid coincide with boundary 
  ! coordinates of previous grid
  grid1%x(1) = grid0%x(1)
  grid1%x(nz+1) = grid0%x(N+1)
      
  ! Find coordinates for interior target values
  do k = 2,nz
    t = target_values(k)
    grid1%x(k) = get_polynomial_coordinate ( grid0, ppoly0, t, degree )
  end do
    
  ! Determine cell widths
  do k = 1,nz
    grid1%h(k) = grid1%x(k+1) - grid1%x(k)
  end do

end subroutine interpolate_grid


!------------------------------------------------------------------------------
! Given target value, find corresponding coordinate for given polynomial
!------------------------------------------------------------------------------
real function get_polynomial_coordinate ( grid, ppoly, target_value, degree )
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
  type(grid1d_t), intent(in)  :: grid               
  type(ppoly_t),  intent(in)  :: ppoly
  real,           intent(in)  :: target_value
  integer,        intent(in)  :: degree

  ! Local variables
  integer            :: i, k;           ! loop indices
  integer            :: N
  integer            :: k_found;        ! index of target cell
  integer            :: iter
  real               :: x_l, x_r;       ! end coordinates of target cell
  real               :: xi0;            ! normalized target coordinate
  real, dimension(5) :: a;              ! polynomial coefficients
  real               :: numerator
  real               :: denominator
  real               :: delta;          ! Newton-Raphson increment
  real               :: x;              ! global target coordinate
  real               :: eps;                ! offset used to get away from
                                        ! boundaries
  real               :: g;              ! gradient during N-R iterations

  N = grid%nb_cells
  eps = NR_OFFSET
  
  k_found = -1

  ! If the target value is outside the range of all values, we
  ! force the target coordinate to be equal to the lowest or
  ! largest value, depending on which bound is overtaken
  if ( target_value .LE. ppoly%E(1,1) ) then
    x = grid%x(1)
    get_polynomial_coordinate = x
    return; ! return because there is no need to look further
  end if
  
  if ( target_value .GE. ppoly%E(N,2) ) then
    x = grid%x(N+1)
    get_polynomial_coordinate = x
    return; ! return because there is no need to look further
  end if
  
  ! Since discontinuous edge values are allowed, we check whether the target
  ! value lies between two discontinuous edge values at interior interfaces
  do k = 2,N
    if ( ( target_value .GE. ppoly%E(k-1,2) ) .AND. &
       ( target_value .LE. ppoly%E(k,1) ) ) then
       x = grid%x(k)
       get_polynomial_coordinate = x
       return;  ! return because there is no need to look further
       exit
    end if   
  end do

  ! At this point, we know that the target value is bounded and does not
  ! lie between discontinuous, monotonic edge values. Therefore,
  ! there is a unique solution. We loop on all cells and find which one
  ! contains the target value. The variable k_found holds the index value
  ! of the cell where the taregt value lies.
  do k = 1,N
    if ( ( target_value .GT. ppoly%E(k,1) ) .AND. &
       ( target_value .LT. ppoly%E(k,2) ) ) then
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
      write(*,*) target_value, ppoly%E(1,1), ppoly%E(N,2)
      write(*,*) 'Could not find target coordinate in ' //&
                 '"get_polynomial_coordinate". This is caused by an '//&
                 'inconsistent interpolant (perhaps not monotonically '//&
                 'increasing)'
      call MOM_error ( FATAL, 'Aborting execution' )
  end if

  ! Reset all polynomial coefficients to 0 and copy those pertaining to 
  ! the found cell
  a(:) = 0.0
  do i = 1,degree+1
    a(i) = ppoly%coefficients(k_found,i)
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

  x_l = grid%x(k_found)
  x_r = grid%x(k_found+1)
  x = x_l + xi0 * grid%h(k_found)
    
  get_polynomial_coordinate = x

end function get_polynomial_coordinate


!------------------------------------------------------------------------------
! Check grid integrity
!------------------------------------------------------------------------------
subroutine check_grid_integrity ( G, h, regridding_opts )
!------------------------------------------------------------------------------
! This routine is called when initializing the regridding options. The 
! objective is to make sure all layers are at least as thick as the minimum
! thickness allowed for regridding purposes (this parameter is set in the
! MOM_input file or defaulted to 1.0e-3). When layers are too thin, they
! are inflated up to the minmum thickness.
!------------------------------------------------------------------------------

  ! Arguments
  type(ocean_grid_type), intent(in)                    :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout)   :: h
  type(regridding_opts_t), intent(in)                  :: regridding_opts

  ! Local variables
  integer           :: i, j, k
  type(grid1d_t)    :: grid

  ! Initialize grid 
  call grid1d_init ( grid, G%ke )

  do i = G%isc,G%iec+1
    do j = G%jsc,G%jec+1
    
      ! Build grid for current column
      do k = 1,G%ke
        grid%h(k) = h(i,j,k)
      end do

      grid%x(1) = 0.0
      do k = 1,G%ke
        grid%x(k+1) = grid%x(k) + grid%h(k)
      end do
      
      call inflate_vanished_layers ( grid, regridding_opts )

      ! Save modified grid
      do k = 1,G%ke
        h(i,j,k) = grid%h(k)
      end do
    
    end do
  end do

  call grid1d_destroy ( grid )

end subroutine check_grid_integrity


!------------------------------------------------------------------------------
! Inflate vanished layers to finite (nonzero) width
!------------------------------------------------------------------------------
subroutine inflate_vanished_layers ( grid, regridding_opts )

  ! Argument
  type(grid1d_t), intent(inout)       :: grid
  type(regridding_opts_t), intent(in) :: regridding_opts
    
  ! Local variable
  integer   :: N
  integer   :: k
  integer   :: k_found
  integer   :: count_nonzero_layers
  real      :: delta
  real      :: correction
  real      :: min_thickness
  real      :: max_thickness

  N = grid%nb_cells
  min_thickness = regridding_opts%min_thickness
  
  ! Count number of nonzero layers
  count_nonzero_layers = 0
  do k = 1,N
    if ( grid%h(k) .GT. min_thickness ) then
      count_nonzero_layers = count_nonzero_layers + 1
    end if
  end do

  ! If all layer thicknesses are greater than the threshold, exit routine
  if ( count_nonzero_layers .eq. N ) return

  ! If all thicknesses are zero, inflate them all and exit
  if ( count_nonzero_layers .eq. 0 ) then  
    do k = 1,N
      grid%h(k) = min_thickness
    end do
    return
  end if    
  
  ! Inflate zero layers
  correction = 0.0
  do k = 1,N
    if ( grid%h(k) .le. min_thickness ) then
      delta = min_thickness - grid%h(k)
      correction = correction + delta
      grid%h(k) = grid%h(k) + delta
    end if  
  end do
  
  ! Modify thicknesses of nonzero layers to ensure volume conservation
  max_thickness = grid%h(1)
  k_found = 1
  do k = 1,grid%nb_cells
    if ( grid%h(k) .gt. max_thickness ) then
      max_thickness = grid%h(k)
      k_found = k
    end if  
  end do
  
  grid%h(k_found) = grid%h(k_found) - correction
  
  ! Redefine grid coordinates according to new layer thicknesses
  grid%x(1) = 0.0
  do k = 1,N
    grid%x(k+1) = grid%x(k) + grid%h(k)
  end do    
  
end subroutine inflate_vanished_layers


!------------------------------------------------------------------------------
! Convective adjustment by swapping layers
!------------------------------------------------------------------------------
subroutine convective_adjustment ( G, h, tv )
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
  real      :: T0, T1;      ! temperatures
  real      :: S0, S1;      ! salinities
  real      :: r0, r1;      ! densities
  real      :: h0, h1
  logical   :: stratified
  
  ! Loop on columns 
  do j = G%jsc,G%jec+1
    do i = G%isc,G%iec+1
        
      ! Compute densities within current water column
      call calculate_density ( tv%T(i,j,:), tv%S(i,j,:), p_column, &
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
            call calculate_density ( tv%T(i,j,k), tv%S(i,j,k), &
                                     p_column(k), &
                                     densities(k), 1, 1, tv%eqn_of_state )
            call calculate_density ( tv%T(i,j,k+1), tv%S(i,j,k+1), &
                                     p_column(k+1), &
                                     densities(k+1), 1, 1, tv%eqn_of_state )
            stratified = .false.
          end if
        end do; ! k 
    
        if ( stratified ) exit        

      end do    

    end do; ! i
  end do; ! j   

end subroutine convective_adjustment


!------------------------------------------------------------------------------
! Allocate memory for regridding
!------------------------------------------------------------------------------
subroutine regridding_memory_allocation ( G, regridding_opts )
!------------------------------------------------------------------------------
! In this routine, we allocate the memory needed to carry out regridding
! steps in the course of the simulation. 

! For example, to compute implicit edge-value estimates, a tridiagonal system
! must be solved. We allocate the needed memory at the beginning of the
! simulation because the number of layers never changes.
!------------------------------------------------------------------------------

  ! Arguments
  type(ocean_grid_type), intent(in)   :: G
  type(regridding_opts_t), intent(in) :: regridding_opts

  ! Local variables
  integer   :: nz
  integer   :: degree;      ! Degree of polynomials used for the reconstruction 
  
  nz = G%ke

  ! Allocate memory for the tridiagonal system
  call tridiagonal_system_0_memory_allocation ( nz )
  call tridiagonal_system_1_memory_allocation ( nz )

  ! Target values
  allocate ( target_values(nz+1) )

  ! Allocate memory for grids
  call grid1d_init ( grid_generic, nz )
  call grid1d_init ( grid_start, nz )
  call grid1d_init ( grid_trans, nz )
  call grid1d_init ( grid_final, nz )

  ! Piecewise polynomials used for remapping
  call ppoly_init ( ppoly_r, nz, 4 )
  
  ! Piecewise polynomials used for interpolation
  call ppoly_init ( ppoly_i, nz, 4 )
  
  ! Generic linear piecewise polynomial
  call ppoly_init ( ppoly_linear, nz, 1 )
  
  ! Generic parabolic piecewise polynomial
  call ppoly_init ( ppoly_parab, nz, 2 )
  
  ! Memory allocation for one column
  allocate ( mapping(nz) ); mapping = 0.0
  allocate ( u_column(nz) ); u_column = 0.0
  allocate ( T_column(nz) ); T_column = 0.0
  allocate ( S_column(nz) ); S_column = 0.0
  allocate ( p_column(nz) ); p_column = 0.0
  allocate ( densities(nz) ); densities = 0.0

  ! Remapping
  call remapping_memory_allocation ( G, regridding_opts )

end subroutine regridding_memory_allocation


!------------------------------------------------------------------------------
! Deallocate memory for regridding
!------------------------------------------------------------------------------
subroutine regridding_memory_deallocation ( regridding_opts )
!------------------------------------------------------------------------------
! In this routine, we reclaim the memory that was allocated for the regridding. 
!------------------------------------------------------------------------------
  
  type(regridding_opts_t), intent(in) :: regridding_opts
  
  ! Reclaim memory for the tridiagonal system
  call tridiagonal_system_0_memory_deallocation ( )
  call tridiagonal_system_1_memory_deallocation ( )
  
  ! Target values
  deallocate ( target_values )

  ! Deallocate memory for grid
  call grid1d_destroy ( grid_generic )
  call grid1d_destroy ( grid_start )
  call grid1d_destroy ( grid_trans )
  call grid1d_destroy ( grid_final )
  
  ! Piecewise polynomials
  call ppoly_destroy ( ppoly_i )
  call ppoly_destroy ( ppoly_r )
  call ppoly_destroy ( ppoly_linear )
  call ppoly_destroy ( ppoly_parab )

  deallocate ( mapping )
  deallocate ( u_column )
  deallocate ( T_column )
  deallocate ( S_column )
  deallocate ( p_column )
  deallocate ( densities )

  ! Remapping
  call remapping_memory_deallocation ( )

end subroutine regridding_memory_deallocation

end module MOM_regridding
