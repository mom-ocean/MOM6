module MOM_ALE
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
use MOM_file_parser,   only : get_param, param_file_type, uppercase
use MOM_EOS,           only : calculate_density

use regrid_grid1d_class, only : grid1D_t, grid1Dconstruct, grid1Ddestroy
use regrid_ppoly_class, only : ppoly_t, ppoly_init, ppoly_destroy
use regrid_edge_values, only : edgeValueArrays
use regrid_edge_values, only : edge_values_implicit_h4
use regrid_edge_values, only : triDiagEdgeWorkAllocate, triDiagEdgeWorkDeallocate
use regrid_edge_slopes, only : edgeSlopeArrays
use regrid_edge_slopes, only : triDiagSlopeWorkAllocate, triDiagSlopeWorkDeallocate

use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation

use P1M_functions, only : P1M_interpolation, P1M_boundary_extrapolation
use P3M_functions, only : P3M_interpolation, P3M_boundary_extrapolation
use MOM_regridding, only : initialize_regridding, regridding_main , end_regridding
use MOM_regridding, only : regridding_CS
use MOM_remapping, only : initialize_remapping, allocate_remapping, remapping_main, end_remapping
use MOM_remapping, only : remapping_CS
use regrid_defs, only : PRESSURE_RECONSTRUCTION_PLM, PRESSURE_RECONSTRUCTION_PPM
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE


implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Private (module-wise) variables
! -----------------------------------------------------------------------------

type, public :: ALE_CS
  private
  ! Generic grid used for various purposes throughout the code (the same 
  ! grid is used to avoid having dynamical memory allocation)
  type(grid1D_t)                  :: grid_generic

  ! Generic linear piecewise polynomial used for various purposes throughout 
  ! the code (the same ppoly is used to avoid having dynamical memory allocation)
  type(ppoly_t)                   :: ppoly_linear

  ! Generic parabolic piecewise polynomial used for various purposes 
  ! throughout the code (the same ppoly is used to avoid having dynamical 
  ! memory allocation)
  type(ppoly_t)                   :: ppoly_parab

  type(grid1D_t)                  :: grid_start      ! starting grid
  type(grid1D_t)                  :: grid_trans      ! transition/iterated grid
  type(grid1D_t)                  :: grid_final      ! final grid
  type(ppoly_t)                   :: ppoly_i         ! interpolation ppoly
  type(ppoly_t)                   :: ppoly_r         ! reconstruction ppoly
                                                    
  ! These variables are private to this module and memory is allocated at
  ! the beginning of the simulation. Note that 'u_column' is a generic variable
  ! used for remapping. It is not necessarily associated with the zonal velocity
  ! component.
  real, dimension(:), allocatable :: target_values
  real, dimension(:), allocatable :: densities
  integer, dimension(:), allocatable :: mapping
  real, dimension(:), allocatable :: u_column
  real, dimension(:), allocatable :: T_column
  real, dimension(:), allocatable :: S_column
  real, dimension(:), allocatable :: p_column

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

  ! Indicates whether integrals for FV pressure gradient calculation will
  ! use reconstruction of T/S.
  ! By default, it is true if regridding has been initialized, otherwise false.
  logical   :: reconstructForPressure = .false.

  ! The form of the reconstruction of T/S for FV pressure gradient calculation.
  ! By default, it is =1 (PLM)
  integer   :: pressureReconstructionScheme

  type(regridding_CS) :: regridCS ! Regridding parameters and work arrays
  type(remapping_CS) :: remapCS ! Remapping parameters and work arrays
  type(edgeValueArrays) :: edgeValueWrk ! Work space for edge values
  type(edgeSlopeArrays) :: edgeSlopeWrk ! Work space for edge slopes
end type

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public initialize_ALE
public end_ALE
public ALE_main 
public pressure_gradient_plm
public pressure_gradient_ppm
public usePressureReconstruction
public pressureReconstructionScheme

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

! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! Initialization of regridding
!------------------------------------------------------------------------------
subroutine initialize_ALE( param_file, G, h, h_aux, &
                                  u, v, tv, CS )
!------------------------------------------------------------------------------
! This routine is typically called (from initialize_MOM in file MOM.F90)
! before the main time integration loop to initialize the regridding stuff.
! We read the MOM_input file to register the values of different
! regridding/remapping parameters.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(param_file_type), intent(in)                      :: param_file
  type(ocean_grid_type), intent(in)                      :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_,C2_), intent(inout) :: h
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout)     :: h_aux
  real, dimension(NIMEMB_,NJMEM_, NKMEM_), intent(inout)    :: u
  real, dimension(NIMEM_,NJMEMB_, NKMEM_), intent(inout)    :: v
  type(thermo_var_ptrs), intent(inout)                   :: tv
  type(ALE_CS), intent(inout)                     :: CS

  ! Local variables
  logical       :: verbose
  integer       :: k, m
  integer       :: i, j

  verbose = .false.
 
  ! Memory allocation for regridding
  call ALE_memory_allocation( G, CS )

  call read_ALE_options( param_file, CS )

  call initialize_regridding( param_file, G, CS%regridCS )

  call initialize_remapping( param_file, G%ke, CS%remapCS )
  call allocate_remapping( CS%remapCS )

  ! Check grid integrity with respect to minimum allowed thickness
  do m = 1,size(h,4)
    call check_grid_integrity( G, h(:,:,:,m), CS )
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
  call regridding_main( CS%remapCS, CS%regridCS, G, h(:,:,:,1), u, v, tv, h_aux )
  call remapping_main( CS%remapCS, G, h(:,:,:,1), h_aux, tv, u, v )
  h(:,:,:,1) = h_aux(:,:,:)
  
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

end subroutine initialize_ALE


!------------------------------------------------------------------------------
! Initialization of regridding options
!------------------------------------------------------------------------------
subroutine read_ALE_options( param_file, CS )
!------------------------------------------------------------------------------
! Read the regridding/remapping parameters in the MOM_input file and
! update the structure that is passed as argument all over the place.
!------------------------------------------------------------------------------

  ! Arguments
  type(param_file_type), intent(in)        :: param_file
  type(ALE_CS), intent(inout)   :: CS
  ! Local variables
  character(len=40)  :: mod = "MOM_ALE" ! This module's name.
  character(len=40)  :: string ! Temporary string

  ! --- TYPE OF VERTICAL GRID ---
  ! This sets which kind of grid we want to use in the vertical. If none
  ! is specified, target interface densities are used to build the grid
  call get_param(param_file, mod, "REGRIDDING_COORDINATE_MODE", string, &
                 "Coordinate mode for vertical regridding.\n"//&
                 "Choose among the following possibilities:\n"//&
                 " LAYER - Isopycnal or stacked shallow water layers\n"//&
                 " Z*    - stetched geopotential z*\n"//&
                 " SIGMA - terrain following coordinates\n"//&
                 " RHO   - continuous isopycnal\n",&
                 default=DEFAULT_COORDINATE_MODE, fail_if_missing=.true.)
  CS%regridding_scheme = coordinateMode(string)

  ! --- INTERPOLATION SCHEME ---
  ! This sets which interpolation scheme we want to use to define the new
  ! grid when regridding is based upon target interface densities. If none
  ! is specified, the p1m h2 interpolation scheme is used.
  call get_param(param_file, mod, "INTERPOLATION_SCHEME", string, &
                 "This sets the interpolation scheme to use to\n"//&
                 "determine the new grid. These parameters are\n"//&
                 "only relevant when REGRIDDING_COORDINATE_MODE is\n"//&
                 "set to a function of state. Otherwise, it is not\n"//&
                 "used. It can be one of the following schemes:\n"//&
                 " P1M_H2     (2nd-order accurate)\n"//&
                 " P1M_H4     (2nd-order accurate)\n"//&
                 " P1M_IH4    (2nd-order accurate)\n"//&
                 " PLM        (2nd-order accurate)\n"//&
                 " PPM_H4     (3rd-order accurate)\n"//&
                 " PPM_IH4    (3rd-order accurate)\n"//&
                 " P3M_IH4IH3 (4th-order accurate)\n"//&
                 " P3M_IH6IH5 (4th-order accurate)\n"//&
                 " PQM_IH4IH3 (4th-order accurate)\n"//&
                 " PQM_IH6IH5 (5th-order accurate)", &
                 default="P1M_H2")
  select case ( uppercase(trim(string)) )
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
    case default ; call MOM_error(FATAL, "read_ALE_options: "//&
       "Unrecognized choice for INTERPOLATION_SCHEME ("//trim(string)//").")
  end select

  ! --- BOUNDARY EXTRAPOLATION --
  ! This sets whether high-order (rather than PCM) reconstruction schemes
  ! should be used within boundary cells
  call get_param(param_file, mod, "BOUNDARY_EXTRAPOLATION", &
                 CS%boundary_extrapolation, &
                 "When defined, a proper high-order reconstruction\n"//&
                 "scheme is used within boundary cells rather\n"//&
                 "than PCM. E.g., if PPM is used for remapping, a\n"//&
                 "PPM reconstruction will also be used within\n"//&
                 "boundary cells.", default=.false.)

  ! --- PRESSURE GRADIENT CALCULATION ---
  call get_param(param_file, mod, "RECONSTRUCT_FOR_PRESSURE", &
                 CS%reconstructForPressure , &
                 "If True, use vertical reconstruction of T/S within\n"//&
                 "the integrals of teh FV pressure gradient calculation.\n"//&
                 "If False, use the constant-by-layer algorithm.\n"//&
                 "By default, this is True when using ALE and False otherwise.", &
                 default=.true. )

  call get_param(param_file, mod, "PRESSURE_RECONSTRUCTION_SCHEME", &
                 CS%pressureReconstructionScheme, &
                 "Type of vertical reconstruction of T/S to use in integrals\n"//&
                 "within the FV pressure gradient calculation."//&
                 " 1: PLM reconstruction.\n"//&
                 " 2: PPM reconstruction.", default=PRESSURE_RECONSTRUCTION_PLM)

  ! --- MINIMUM THICKNESS ---
  call get_param(param_file, mod, "MIN_THICKNESS", &
                 CS%min_thickness, &
                 "When regridding, this is the minimum layer\n"//&
                 "thickness allowed.", default=1.e-3 )
  
end subroutine read_ALE_options



!------------------------------------------------------------------------------
! End of regridding (memory deallocation)
!------------------------------------------------------------------------------
subroutine end_ALE(CS)
!------------------------------------------------------------------------------
! This routine is typically called (from MOM_end in file MOM.F90)
! after the main time integration loop to deallocate the regridding stuff.
!------------------------------------------------------------------------------
  type(ALE_CS), intent(inout) :: CS
  
  ! Deallocate memory used for the regridding
  call end_remapping( CS%remapCS )
  call end_regridding( CS%regridCS )
  call ALE_memory_deallocation( CS )

end subroutine end_ALE


!------------------------------------------------------------------------------
! Dispatching regridding routine: regridding & remapping
!------------------------------------------------------------------------------
subroutine ALE_main( G, h, h_new, u, v, tv, CS )
!------------------------------------------------------------------------------
! This routine takes care of (1) building a new grid and (2) remapping between
! the old grid and the new grid. The creation of the new grid can be based
! on z coordinates, target interface densities, sigma coordinates or any
! arbitrary coordinate system.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(ocean_grid_type), intent(in)                    :: &
  G      ! Ocean grid informations
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout)   :: &
  h      ! Current 3D grid obtained after the last time step
  real, dimension(NIMEM_,NJMEM_, NKMEM_), intent(inout)   :: &
  h_new  ! The new 3D grid obtained via regridding
  real, dimension(NIMEMB_,NJMEM_, NKMEM_), intent(inout)  :: &
  u      ! Zonal velocity field
  real, dimension(NIMEM_,NJMEMB_, NKMEM_), intent(inout)  :: &
  v      ! Meridional velocity field
  type(thermo_var_ptrs), intent(inout)                 :: &
  tv     ! Thermodynamical variables (T, S, ...)  
  type(ALE_CS), intent(inout) :: CS ! Regridding parameters and options

  ! Local variables
  
  ! Build new grid. The new grid is stored in h_new. The old grid is h.
  ! Both are needed for the subsequent remapping of variables.
  call regridding_main( CS%remapCS, CS%regridCS, G, h, u, v, tv, h_new )
  
  ! Remap all variables from old grid h onto new grid h_new
  call remapping_main( CS%remapCS, G, h, h_new, tv, u, v )
  
  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
  h(:,:,:) = h_new(:,:,:)

end subroutine ALE_main


!------------------------------------------------------------------------------
! Use plm reconstruction for pressure gradient (determine edge values)
!------------------------------------------------------------------------------
subroutine pressure_gradient_plm( CS, S_t, S_b, T_t, T_b, G, tv, h )
!------------------------------------------------------------------------------
! By using a PLM (limited piecewise linear method) reconstruction, this 
! routine determines the edge values for the salinity and temperature 
! within each layer. These edge values are returned and are used to compute 
! the pressure gradient (by computing the densities).
!------------------------------------------------------------------------------

  ! Arguments
  type(ALE_CS), intent(inout) :: CS ! Regridding parameters and options
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: &
  S_t, S_b  ! Salinity at the top and bottom edges of each layer
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: &
  T_t, T_b  ! Temperature at the top and bottom edges of each layer
  type(ocean_grid_type), intent(in)                 :: G
  type(thermo_var_ptrs), intent(in)                 :: tv
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: &
  h         ! Three-dimensional ocean grid

  ! Local variables
  integer           :: i, j, k

  ! NOTE: the variables 'CS%grid_generic' and 'CS%ppoly_linear' are declared at
  ! the module level. Memory is allocated once at the beginning of the run
  ! in 'ALE_memory_allocation'.

  ! Determine reconstruction within each column
  do i = G%isc,G%iec+1
    do j = G%jsc,G%jec+1
     
      ! Build current grid
      CS%grid_generic%h(:) = h(i,j,:)
      CS%grid_generic%x(1) = 0.0
      do k = 1,G%ke
        CS%grid_generic%x(k+1) = CS%grid_generic%x(k) + CS%grid_generic%h(k)
      end do
      
      ! Reconstruct salinity profile    
      CS%ppoly_linear%E = 0.0
      CS%ppoly_linear%coefficients = 0.0
      call PLM_reconstruction( CS%grid_generic, tv%S(i,j,:), CS%ppoly_linear )
      call PLM_boundary_extrapolation( CS%grid_generic, tv%S(i,j,:), CS%ppoly_linear )
      
      do k = 1,G%ke
        S_t(i,j,k) = CS%ppoly_linear%E(k,1)
        S_b(i,j,k) = CS%ppoly_linear%E(k,2)
      end do
      
      ! Reconstruct temperature profile 
      CS%ppoly_linear%E = 0.0
      CS%ppoly_linear%coefficients = 0.0
      call PLM_reconstruction( CS%grid_generic, tv%T(i,j,:), CS%ppoly_linear )
      call PLM_boundary_extrapolation( CS%grid_generic, tv%T(i,j,:), CS%ppoly_linear )
      
      do k = 1,G%ke
        T_t(i,j,k) = CS%ppoly_linear%E(k,1)
        T_b(i,j,k) = CS%ppoly_linear%E(k,2)
      end do
      
    end do
  end do

end subroutine pressure_gradient_plm


!------------------------------------------------------------------------------
! Use ppm reconstruction for pressure gradient (determine edge values)
!------------------------------------------------------------------------------
subroutine pressure_gradient_ppm( CS, S_t, S_b, T_t, T_b, G, tv, h )
!------------------------------------------------------------------------------
! By using a PPM (limited piecewise linear method) reconstruction, this 
! routine determines the edge values for the salinity and temperature 
! within each layer. These edge values are returned and are used to compute 
! the pressure gradient (by computing the densities).
!------------------------------------------------------------------------------

  ! Arguments
  type(ALE_CS), intent(inout) :: CS
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: &
  S_t, S_b  ! Salinity at the top and bottom edges of each layer
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(inout) :: &
  T_t, T_b  ! Temperature at the top and bottom edges of each layer
  type(ocean_grid_type), intent(in)                 :: G
  type(thermo_var_ptrs), intent(in)                 :: tv
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: &
  h         ! Three-dimensional ocean grid

  ! Local variables
  integer           :: i, j, k

  ! NOTE: the variables 'CS%grid_generic' and 'CS%ppoly_parab' are declared at
  ! the module level. Memory is allocated once at the beginning of the run
  ! in 'ALE_memory_allocation'.

  ! Determine reconstruction within each column
  do i = G%isc,G%iec+1
    do j = G%jsc,G%jec+1
     
      ! Build current grid
      CS%grid_generic%h(:) = h(i,j,:)
      CS%grid_generic%x(1) = 0.0
      do k = 1,G%ke
        CS%grid_generic%x(k+1) = CS%grid_generic%x(k) + CS%grid_generic%h(k)
      end do
      
      ! Reconstruct salinity profile    
      CS%ppoly_parab%E = 0.0
      CS%ppoly_parab%coefficients = 0.0
      call edge_values_implicit_h4( CS%grid_generic, CS%edgeValueWrk, tv%S(i,j,:), CS%ppoly_parab%E )
      call PPM_reconstruction( CS%grid_generic, tv%S(i,j,:), CS%ppoly_parab )
      call PPM_boundary_extrapolation( CS%grid_generic, tv%S(i,j,:), CS%ppoly_parab )
      
      do k = 1,G%ke
        S_t(i,j,k) = CS%ppoly_parab%E(k,1)
        S_b(i,j,k) = CS%ppoly_parab%E(k,2)
      end do
      
      ! Reconstruct temperature profile 
      CS%ppoly_parab%E = 0.0
      CS%ppoly_parab%coefficients = 0.0
      call edge_values_implicit_h4( CS%grid_generic, CS%edgeValueWrk, tv%T(i,j,:), CS%ppoly_parab%E )
      call PPM_reconstruction( CS%grid_generic, tv%T(i,j,:), CS%ppoly_parab )
      call PPM_boundary_extrapolation( CS%grid_generic, tv%T(i,j,:), CS%ppoly_parab )
      
      do k = 1,G%ke
        T_t(i,j,k) = CS%ppoly_parab%E(k,1)
        T_b(i,j,k) = CS%ppoly_parab%E(k,2)
      end do
      
    end do
  end do

end subroutine pressure_gradient_ppm


!------------------------------------------------------------------------------
! Check grid integrity
!------------------------------------------------------------------------------
subroutine check_grid_integrity( G, h, CS )
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
  type(ALE_CS), intent(in)                  :: CS

  ! Local variables
  integer           :: i, j, k
  type(grid1D_t)    :: grid

  ! Initialize grid 
  call grid1Dconstruct( grid, G%ke )

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
      
      call inflate_vanished_layers( grid, CS )

      ! Save modified grid
      do k = 1,G%ke
        h(i,j,k) = grid%h(k)
      end do
    
    end do
  end do

  call grid1Ddestroy( grid )

end subroutine check_grid_integrity


!------------------------------------------------------------------------------
! Inflate vanished layers to finite (nonzero) width
!------------------------------------------------------------------------------
subroutine inflate_vanished_layers( grid, CS )

  ! Argument
  type(grid1D_t), intent(inout)       :: grid
  type(ALE_CS), intent(in) :: CS
    
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
  min_thickness = CS%min_thickness
  
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
subroutine convective_adjustment( CS, G, h, tv )
!------------------------------------------------------------------------------
! Check each water column to see whether it is stratified. If not, sort the
! layers by successive swappings of water masses (bubble sort algorithm)
!------------------------------------------------------------------------------

  ! Arguments
  type(ALE_CS), intent(inout) :: CS
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
  
  ! Loop on columns 
  do j = G%jsc,G%jec+1
    do i = G%isc,G%iec+1
        
      ! Compute densities within current water column
      call calculate_density( tv%T(i,j,:), tv%S(i,j,:), CS%p_column, &
                              CS%densities, 1, G%ke, tv%eqn_of_state )
     
      ! Repeat restratification until complete  
      do

        stratified = .true.
        do k = 1,G%ke-1
          ! Gather information of current and next cells
          T0 = tv%T(i,j,k)
          T1 = tv%T(i,j,k+1)
          S0 = tv%S(i,j,k)
          S1 = tv%S(i,j,k+1)
          r0 = CS%densities(k)
          r1 = CS%densities(k+1)
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
                                     CS%p_column(k), &
                                     CS%densities(k), 1, 1, tv%eqn_of_state )
            call calculate_density( tv%T(i,j,k+1), tv%S(i,j,k+1), &
                                     CS%p_column(k+1), &
                                     CS%densities(k+1), 1, 1, tv%eqn_of_state )
            stratified = .false.
          end if
        end do  ! k 
    
        if ( stratified ) exit        

      end do    

    end do  ! i
  end do  ! j   

end subroutine convective_adjustment


!------------------------------------------------------------------------------
! Allocate memory for regridding
!------------------------------------------------------------------------------
subroutine ALE_memory_allocation( G, CS )
!------------------------------------------------------------------------------
! In this routine, we allocate the memory needed to carry out regridding
! steps in the course of the simulation. 

! For example, to compute implicit edge-value estimates, a tridiagonal system
! must be solved. We allocate the needed memory at the beginning of the
! simulation because the number of layers never changes.
!------------------------------------------------------------------------------

  ! Arguments
  type(ocean_grid_type), intent(in)   :: G
  type(ALE_CS), intent(inout) :: CS

  ! Local variables
  integer   :: nz
  integer   :: degree       ! Degree of polynomials used for the reconstruction 
  
  nz = G%ke

  ! Allocate memory for the tridiagonal system
  call triDiagEdgeWorkAllocate( nz, CS%edgeValueWrk )
  call triDiagSlopeWorkAllocate( nz, CS%edgeSlopeWrk )

  ! Target values
  allocate( CS%target_values(nz+1) )

  ! Allocate memory for grids
  call grid1Dconstruct( CS%grid_generic, nz )
  call grid1Dconstruct( CS%grid_start, nz )
  call grid1Dconstruct( CS%grid_trans, nz )
  call grid1Dconstruct( CS%grid_final, nz )

  ! Piecewise polynomials used for remapping
  call ppoly_init( CS%ppoly_r, nz, 4 )
  
  ! Piecewise polynomials used for interpolation
  call ppoly_init( CS%ppoly_i, nz, 4 )
  
  ! Generic linear piecewise polynomial
  call ppoly_init( CS%ppoly_linear, nz, 1 )
  
  ! Generic parabolic piecewise polynomial
  call ppoly_init( CS%ppoly_parab, nz, 2 )
  
  ! Memory allocation for one column
  allocate( CS%mapping(nz) ); CS%mapping = 0
  allocate( CS%u_column(nz) ); CS%u_column = 0.0
  allocate( CS%T_column(nz) ); CS%T_column = 0.0
  allocate( CS%S_column(nz) ); CS%S_column = 0.0
  allocate( CS%p_column(nz) ); CS%p_column = 0.0
  allocate( CS%densities(nz) ); CS%densities = 0.0

end subroutine ALE_memory_allocation


!------------------------------------------------------------------------------
! Deallocate memory for regridding
!------------------------------------------------------------------------------
subroutine ALE_memory_deallocation( CS )
!------------------------------------------------------------------------------
! In this routine, we reclaim the memory that was allocated for the regridding. 
!------------------------------------------------------------------------------
  
  type(ALE_CS), intent(inout) :: CS
  
  ! Reclaim memory for the tridiagonal system
  call triDiagEdgeWorkDeallocate( CS%edgeValueWrk )
  call triDiagSlopeWorkDeallocate( CS%edgeSlopeWrk )
  
  ! Target values
  deallocate( CS%target_values )

  ! Deallocate memory for grid
  call grid1Ddestroy( CS%grid_generic )
  call grid1Ddestroy( CS%grid_start )
  call grid1Ddestroy( CS%grid_trans )
  call grid1Ddestroy( CS%grid_final )
  
  ! Piecewise polynomials
  call ppoly_destroy( CS%ppoly_i )
  call ppoly_destroy( CS%ppoly_r )
  call ppoly_destroy( CS%ppoly_linear )
  call ppoly_destroy( CS%ppoly_parab )

  deallocate( CS%mapping )
  deallocate( CS%u_column )
  deallocate( CS%T_column )
  deallocate( CS%S_column )
  deallocate( CS%p_column )
  deallocate( CS%densities )

end subroutine ALE_memory_deallocation

logical function usePressureReconstruction(CS)
  type(ALE_CS), intent(in) :: CS
  usePressureReconstruction=CS%reconstructForPressure
end function usePressureReconstruction

integer function pressureReconstructionScheme(CS)
  type(ALE_CS), intent(in) :: CS
  pressureReconstructionScheme=CS%pressureReconstructionScheme
end function pressureReconstructionScheme

end module MOM_ALE
