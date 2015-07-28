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
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use MOM_file_parser,   only : get_param, param_file_type, log_param
use MOM_io,            only : file_exists, field_exists, MOM_read_data
use MOM_io,            only : vardesc, fieldtype, SINGLE_FILE
use MOM_io,            only : create_file, write_field, close_file, slasher
use MOM_EOS,           only : calculate_density
use MOM_string_functions, only : uppercase, extractWord
use MOM_verticalGrid,  only : verticalGrid_type
use regrid_edge_values, only : edge_values_implicit_h4
use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation
use P1M_functions, only : P1M_interpolation, P1M_boundary_extrapolation
use P3M_functions, only : P3M_interpolation, P3M_boundary_extrapolation
use MOM_regridding, only : initialize_regridding, regridding_main , end_regridding
use MOM_regridding, only : uniformResolution
use MOM_regridding, only : check_grid_integrity, setCoordinateResolution
use MOM_regridding, only : setcoordinateinterfaces
use MOM_regridding, only : regriddingCoordinateModeDoc, DEFAULT_COORDINATE_MODE
use MOM_regridding, only : regriddingInterpSchemeDoc, regriddingDefaultInterpScheme
use MOM_regridding, only : setRegriddingBoundaryExtrapolation
use MOM_regridding, only : regriddingDefaultBoundaryExtrapolation
use MOM_regridding, only : setRegriddingMinimumThickness, regriddingDefaultMinThickness
use MOM_regridding, only : regridding_CS
use MOM_regridding, only : getCoordinateInterfaces, getCoordinateResolution
use MOM_regridding, only : getCoordinateUnits, getCoordinateShortName
use MOM_regridding, only : getStaticThickness
use MOM_remapping, only : initialize_remapping, remapping_main, end_remapping
use MOM_remapping, only : remappingSchemesDoc, remappingDefaultScheme
use MOM_remapping, only : remapDisableBoundaryExtrapolation, remapEnableBoundaryExtrapolation
use MOM_remapping, only : remapping_CS
use regrid_defs, only : PRESSURE_RECONSTRUCTION_PLM
!use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : coordinateUnits, coordinateMode
use regrid_consts, only : REGRIDDING_ZSTAR


implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Private (module-wise) variables
! -----------------------------------------------------------------------------

type, public :: ALE_CS
  private

  ! Indicate whether high-order boundary extrapolation should be used within
  ! boundary cells
  logical :: boundary_extrapolation_for_pressure

  ! Indicates whether integrals for FV pressure gradient calculation will
  ! use reconstruction of T/S.
  ! By default, it is true if regridding has been initialized, otherwise false.
  logical :: reconstructForPressure = .false.

  ! The form of the reconstruction of T/S for FV pressure gradient calculation.
  ! By default, it is =1 (PLM)
  integer :: pressureReconstructionScheme

  type(regridding_CS) :: regridCS ! Regridding parameters and work arrays
  type(remapping_CS) :: remapCS ! Remapping parameters and work arrays

  ! Used only for queries, not directly by this module
  integer :: nk

  integer :: degree_linear=1  ! Degree of linear piecewise polynomial
  integer :: degree_parab=2   ! Degree of parabolic piecewise polynomial


  ! Work space for communicating between regridding and remapping
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NK_INTERFACE_) :: dzRegrid

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
public adjustGridForIntegrity
public ALE_initRegridding
public ALE_getCoordinate
public ALE_getCoordinateUnits
public ALE_writeCoordinateFile
public ALE_updateVerticalGridType
public ALE_initThicknessToCoord

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
subroutine initialize_ALE( param_file, G, CS )
!------------------------------------------------------------------------------
! This routine is typically called (from initialize_MOM in file MOM.F90)
! before the main time integration loop to initialize the regridding stuff.
! We read the MOM_input file to register the values of different
! regridding/remapping parameters.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(param_file_type), intent(in)                      :: param_file
  type(ocean_grid_type), intent(in)                      :: G
  type(ALE_CS), pointer                                  :: CS

  ! Local variables
  real, dimension(:), allocatable :: dz
  character(len=40)  :: mod = "MOM_ALE" ! This module's name.
  character(len=80) :: string ! Temporary strings

  if (associated(CS)) then
    call MOM_error(WARNING, "initialize_ALE called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  ! Memory allocation for regridding
  call ALE_memory_allocation( G, CS )

  ! --- BOUNDARY EXTRAPOLATION --
  ! This sets whether high-order (rather than PCM) reconstruction schemes
  ! should be used within boundary cells
  call get_param(param_file, mod, "BOUNDARY_EXTRAPOLATION_PRESSURE", &
                 CS%boundary_extrapolation_for_pressure, &
                 "When defined, the reconstruction is extrapolated\n"//&
                 "within boundary cells rather than assume PCM for the.\n"//&
                 "calculation of pressure. e.g. if PPM is used, a\n"//&
                 "PPM reconstruction will also be used within\n"//&
                 "boundary cells.", default=.true.)

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

  ! Initialize and configure regridding
  allocate( dz(G%ke) )
  call ALE_initRegridding( G, param_file, mod, CS%regridCS, dz )
  deallocate( dz )

  ! Initialize and configure remapping
  call get_param(param_file, mod, "REMAPPING_SCHEME", string, &
                 "This sets the reconstruction scheme used\n"//&
                 "for vertical remapping for all variables.\n"//&
                 "It can be one of the following schemes:\n"//&
                 trim(remappingSchemesDoc), default=remappingDefaultScheme)
  call initialize_remapping( G%ke, string, CS%remapCS )
  call remapDisableBoundaryExtrapolation( CS%remapCS )

  ! Keep a record of values for subsequent queries
  CS%nk = G%ke

end subroutine initialize_ALE


!------------------------------------------------------------------------------
! Crudely adjust (initial) grid for integrity
!------------------------------------------------------------------------------
subroutine adjustGridForIntegrity( CS, G, h )
!------------------------------------------------------------------------------
! This routine is typically called (from initialize_MOM in file MOM.F90)
! before the main time integration loop to initialize the regridding stuff.
! We read the MOM_input file to register the values of different
! regridding/remapping parameters.
!------------------------------------------------------------------------------
  
  ! Arguments
  type(ALE_CS), pointer                                  :: CS
  type(ocean_grid_type), intent(in)                      :: G
  real, dimension(NIMEM_,NJMEM_, NKMEM_),  intent(inout) :: h

  call check_grid_integrity( CS%regridCS, G, h(:,:,:) )

end subroutine adjustGridForIntegrity


!------------------------------------------------------------------------------
! End of regridding (memory deallocation)
!------------------------------------------------------------------------------
subroutine end_ALE(CS)
!------------------------------------------------------------------------------
! This routine is typically called (from MOM_end in file MOM.F90)
! after the main time integration loop to deallocate the regridding stuff.
!------------------------------------------------------------------------------
  type(ALE_CS), pointer :: CS
  
  ! Deallocate memory used for the regridding
  call end_remapping( CS%remapCS )
  call end_regridding( CS%regridCS )
  call ALE_memory_deallocation( CS )

  deallocate(CS)

end subroutine end_ALE


!------------------------------------------------------------------------------
! Dispatching regridding routine: regridding & remapping
!------------------------------------------------------------------------------
subroutine ALE_main( G, h, u, v, tv, CS )
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
  real, dimension(NIMEMB_,NJMEM_, NKMEM_), intent(inout)  :: &
  u      ! Zonal velocity field
  real, dimension(NIMEM_,NJMEMB_, NKMEM_), intent(inout)  :: &
  v      ! Meridional velocity field
  type(thermo_var_ptrs), intent(inout)                 :: &
  tv     ! Thermodynamical variables (T, S, ...)  
  type(ALE_CS), intent(inout) :: CS ! Regridding parameters and options

  ! Local variables
  integer :: nk, i, j, k, isd, ied, jsd, jed
  
  ! Build new grid. The new grid is stored in h_new. The old grid is h.
  ! Both are needed for the subsequent remapping of variables.
  call regridding_main( CS%remapCS, CS%regridCS, G, h, tv, CS%dzRegrid )
  
  ! Remap all variables from old grid h onto new grid h_new
  call remapping_main( CS%remapCS, G, h, -CS%dzRegrid, tv, u, v )
  
  ! Override old grid with new one. The new grid 'h_new' is built in
  ! one of the 'build_...' routines above.
  nk = G%ke; isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,nk,h,CS)
  do k = 1,nk
    do j = jsd,jed ; do i = isd,ied   
      h(i,j,k) = h(i,j,k) + ( CS%dzRegrid(i,j,k) - CS%dzRegrid(i,j,k+1) )
    enddo ; enddo
  enddo

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
  integer :: i, j, k
  real    :: hTmp(G%ke)
  real    :: tmp(G%ke)
  real, dimension(CS%nk,2) :: &
      ppoly_linear_E            !Edge value of polynomial
  real, dimension(CS%nk,CS%degree_linear+1) :: &
      ppoly_linear_coefficients !Coefficients of polynomial


  ! NOTE: the variables 'CS%grid_generic' and 'CS%ppoly_linear' are declared at
  ! the module level. Memory is allocated once at the beginning of the run
  ! in 'ALE_memory_allocation'.

  ! Determine reconstruction within each column
!$OMP parallel do default(none) shared(G,h,tv,CS,S_t,S_b,T_t,T_b)                     &
!$OMP                          private(hTmp,ppoly_linear_E,ppoly_linear_coefficients,tmp)
  do j = G%jsc,G%jec+1
    do i = G%isc,G%iec+1
      ! Build current grid
      hTmp(:) = h(i,j,:)*G%H_to_m
      tmp(:) = tv%S(i,j,:)
      ! Reconstruct salinity profile    
      ppoly_linear_E = 0.0
      ppoly_linear_coefficients = 0.0
      call PLM_reconstruction( G%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PLM_boundary_extrapolation( G%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      
      do k = 1,G%ke
        S_t(i,j,k) = ppoly_linear_E(k,1)
        S_b(i,j,k) = ppoly_linear_E(k,2)
      end do
      
      ! Reconstruct temperature profile 
      ppoly_linear_E = 0.0
      ppoly_linear_coefficients = 0.0
      tmp(:) = tv%T(i,j,:)
      call PLM_reconstruction( G%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PLM_boundary_extrapolation( G%ke, hTmp, tmp, ppoly_linear_E, ppoly_linear_coefficients )
      
      do k = 1,G%ke
        T_t(i,j,k) = ppoly_linear_E(k,1)
        T_b(i,j,k) = ppoly_linear_E(k,2)
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
  integer :: i, j, k
  real    :: hTmp(G%ke)
  real    :: tmp(G%ke)
  real, dimension(CS%nk,2) :: &
      ppoly_parab_E            !Edge value of polynomial
  real, dimension(CS%nk,CS%degree_parab+1) :: &
      ppoly_parab_coefficients !Coefficients of polynomial


  ! NOTE: the variables 'CS%grid_generic' and 'CS%ppoly_parab' are declared at
  ! the module level. Memory is allocated once at the beginning of the run
  ! in 'ALE_memory_allocation'.

  ! Determine reconstruction within each column
!$OMP parallel do default(none) shared(G,h,tv,CS,S_t,S_b,T_t,T_b) &
!$OMP                          private(hTmp,tmp,ppoly_parab_E,ppoly_parab_coefficients)
  do j = G%jsc,G%jec+1
    do i = G%isc,G%iec+1
     
      ! Build current grid
      hTmp(:) = h(i,j,:) * G%H_to_m
      tmp(:) = tv%S(i,j,:)
      
      ! Reconstruct salinity profile    
      ppoly_parab_E = 0.0
      ppoly_parab_coefficients = 0.0
      call edge_values_implicit_h4( G%ke, hTmp, tmp, ppoly_parab_E )
      call PPM_reconstruction( G%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PPM_boundary_extrapolation( G%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      
      do k = 1,G%ke
        S_t(i,j,k) = ppoly_parab_E(k,1)
        S_b(i,j,k) = ppoly_parab_E(k,2)
      end do
      
      ! Reconstruct temperature profile 
      ppoly_parab_E = 0.0
      ppoly_parab_coefficients = 0.0
      tmp(:) = tv%T(i,j,:)
      call edge_values_implicit_h4( G%ke, hTmp, tmp, ppoly_parab_E )
      call PPM_reconstruction( G%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      if (CS%boundary_extrapolation_for_pressure) call &
        PPM_boundary_extrapolation( G%ke, hTmp, tmp, ppoly_parab_E, ppoly_parab_coefficients )
      
      do k = 1,G%ke
        T_t(i,j,k) = ppoly_parab_E(k,1)
        T_b(i,j,k) = ppoly_parab_E(k,2)
      end do
      
    end do
  end do

end subroutine pressure_gradient_ppm


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
  
  nz = G%ke

  ! Work space
  ALLOC_(CS%dzRegrid(G%isd:G%ied,G%jsd:G%jed,nz+1)); CS%dzRegrid(:,:,:) = 0.
  
end subroutine ALE_memory_allocation


!------------------------------------------------------------------------------
! Deallocate memory for regridding
!------------------------------------------------------------------------------
subroutine ALE_memory_deallocation( CS )
!------------------------------------------------------------------------------
! In this routine, we reclaim the memory that was allocated for the regridding. 
!------------------------------------------------------------------------------
  
  type(ALE_CS), intent(inout) :: CS
  
  ! Work space
  DEALLOC_(CS%dzRegrid)

end subroutine ALE_memory_deallocation

logical function usePressureReconstruction(CS)
  type(ALE_CS), pointer :: CS
  if (associated(CS)) then
    usePressureReconstruction=CS%reconstructForPressure
  else
    usePressureReconstruction=.false.
  endif
end function usePressureReconstruction

integer function pressureReconstructionScheme(CS)
  type(ALE_CS), pointer :: CS
  if (associated(CS)) then
    pressureReconstructionScheme=CS%pressureReconstructionScheme
  else
    pressureReconstructionScheme=-1
  endif
end function pressureReconstructionScheme

subroutine ALE_initRegridding( G, param_file, mod, regridCS, dz )
  ! Arguments
  type(ocean_grid_type), intent(in)  :: G
  type(param_file_type), intent(in)  :: param_file
  character(len=*),      intent(in)  :: mod ! Name of calling module
  type(regridding_CS),   intent(out) :: regridCS ! Regridding parameters and work arrays
  real, dimension(:),    intent(out) :: dz ! Resolution (thickness) in units of coordinate
  ! Local variables
  character(len=80) :: string, varName ! Temporary strings
  character(len=40) :: coordMode, interpScheme, coordUnits ! Temporary strings
  character(len=200) :: inputdir, fileName
  character(len=320) :: message ! Temporary strings
  integer :: ke
  logical :: tmpLogical
  real :: tmpReal

  ke = size(dz) ! Number of levels in resolution vector

  call get_param(param_file, mod, "REGRIDDING_COORDINATE_MODE", coordMode, &
                 "Coordinate mode for vertical regridding.\n"//&
                 "Choose among the following possibilities:\n"//&
                 trim(regriddingCoordinateModeDoc),&
                 default=DEFAULT_COORDINATE_MODE, fail_if_missing=.true.)
  call get_param(param_file, mod, "REGRIDDING_COORDINATE_UNITS", coordUnits, &
                 "Units of the regridding coordinuate.",&
                 default=coordinateUnits(coordMode))

  call get_param(param_file, mod, "INTERPOLATION_SCHEME", interpScheme, &
                 "This sets the interpolation scheme to use to\n"//&
                 "determine the new grid. These parameters are\n"//&
                 "only relevant when REGRIDDING_COORDINATE_MODE is\n"//&
                 "set to a function of state. Otherwise, it is not\n"//&
                 "used. It can be one of the following schemes:\n"//&
                 trim(regriddingInterpSchemeDoc),&
                 default=regriddingDefaultInterpScheme)
  call initialize_regridding( G%ke, coordMode, interpScheme, regridCS )

  call get_param(param_file, mod, "ALE_COORDINATE_CONFIG", string, &
                 "Determines how to specify the coordinate\n"//&
                 "resolution. Valid options are:\n"//&
                 " PARAM       - use the vector-parameter ALE_RESOLUTION\n"//&
                 " UNIFORM     - uniformly distributed\n"//&
                 " FILE:string - read from a file. The string specifies\n"//&
                 "               the filename and variable name, separated\n"//&
                 "               by a comma or space, e.g. FILE:lev.nc,Z",&
                 default='UNIFORM')
  message = "The distribution of vertical resolution for the target\n"//&
            "grid used for Eulerian-like coordinates. For example,\n"//&
            "in z-coordinate mode, the parameter is a list of level\n"//&
            "thicknesses (in m). In sigma-coordinate mode, the list\n"//&
            "is of non-dimensional fractions of the water column."
  select case ( trim(string) )
    case ("UNIFORM")
      dz(:) = uniformResolution(G%ke, coordMode, G%max_depth, &
                 G%Rlay(1)+0.5*(G%Rlay(1)-G%Rlay(2)), &
                 G%Rlay(G%ke)+0.5*(G%Rlay(G%ke)-G%Rlay(G%ke-1)) )
      call log_param(param_file, mod, "!ALE_RESOLUTION", dz, &
                   trim(message), units=trim(coordUnits))
    case ("PARAM")
      call get_param(param_file, mod, "ALE_RESOLUTION", dz, &
                   trim(message), units=trim(coordUnits), fail_if_missing=.true.)
    case default 
      if (index(trim(string),'FILE:')==1) then
        call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
        inputdir = slasher(inputdir)

        if (string(6:6)=='.' .or. string(6:6)=='/') then
          ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
          fileName = trim( extractWord(trim(string(6:80)), 1) )
        else
          ! Otherwise assume we should look for the file in INPUTDIR
          fileName = trim(inputdir) // trim( extractWord(trim(string(6:80)), 1) )
        endif
        if (.not. file_exists(fileName)) call MOM_error(FATAL,"ALE_initRegridding: "// &
          "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

        varName = trim( extractWord(trim(string(6:80)), 2) )
        if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,"ALE_initRegridding: "// &
          "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
        if (len_trim(varName)==0) then
          if (field_exists(fileName,'dz')) then; varName = 'dz'
          elseif (field_exists(fileName,'dsigma')) then; varName = 'dsigma'
          elseif (field_exists(fileName,'ztest')) then; varName = 'ztest'
          endif
        endif
        if (len_trim(varName)==0) call MOM_error(FATAL,"ALE_initRegridding: "// &
          "Coordinate variable not specified and none could be guessed.")
        call MOM_read_data(trim(fileName), trim(varName), dz)
        call log_param(param_file, mod, "!ALE_RESOLUTION", dz, &
                   trim(message), units=coordinateUnits(coordMode))
      else
        call MOM_error(FATAL,"ALE_initRegridding: "// &
          "Unrecognized coordinate configuraiton"//trim(string))
      endif
  end select
  if (coordinateMode(coordMode) == REGRIDDING_ZSTAR) then
    ! Adjust target grid to be consistent with G%max_depth
    ! This is a work around to the from_Z initialization...  ???
    tmpReal = sum( dz(:) )
    if (tmpReal < G%max_depth) then
      dz(ke) = dz(ke) + ( G%max_depth - tmpReal )
    elseif (tmpReal > G%max_depth) then
      if ( dz(ke) + ( G%max_depth - tmpReal ) > 0. ) then
        dz(ke) = dz(ke) + ( G%max_depth - tmpReal )
      else
        call MOM_error(FATAL,"ALE_initRegridding: "// &
          "MAaIMUMX_DEPTH was too shallow to adjust bottom layer of DZ!"//trim(string))
      endif
    endif
  endif
  call setCoordinateResolution( dz, regridCS )
  call setCoordinateInterfaces( G, regridCS )

  call get_param(param_file, mod, "MIN_THICKNESS", tmpReal, &
                 "When regridding, this is the minimum layer\n"//&
                 "thickness allowed.", units="m",&
                 default=regriddingDefaultMinThickness )
  call setRegriddingMinimumThickness( tmpReal, regridCS )

  call get_param(param_file, mod, "BOUNDARY_EXTRAPOLATION", tmpLogical, &
                 "When defined, a proper high-order reconstruction\n"//&
                 "scheme is used within boundary cells rather\n"//&
                 "than PCM. E.g., if PPM is used for remapping, a\n"//&
                 "PPM reconstruction will also be used within\n"//&
                 "boundary cells.", default=regriddingDefaultBoundaryExtrapolation)
  call setRegriddingBoundaryExtrapolation( tmpLogical, regridCS )

end subroutine ALE_initRegridding

!------------------------------------------------------------------------------
! Query the target coordinate interfaces positions
!------------------------------------------------------------------------------
function ALE_getCoordinate( CS )
  type(ALE_CS), pointer    :: CS
  real, dimension(CS%nk+1) :: ALE_getCoordinate

  ALE_getCoordinate(:) = getCoordinateInterfaces( CS%regridCS )

end function ALE_getCoordinate

!------------------------------------------------------------------------------
! Query the target coordinate units
!------------------------------------------------------------------------------
function ALE_getCoordinateUnits( CS )
  type(ALE_CS), pointer    :: CS
  character(len=20)               :: ALE_getCoordinateUnits

  ALE_getCoordinateUnits = getCoordinateUnits( CS%regridCS )

end function ALE_getCoordinateUnits

!------------------------------------------------------------------------------
! Update the vertical grid type with ALE information
!------------------------------------------------------------------------------
subroutine ALE_updateVerticalGridType( CS, GV )
  type(ALE_CS),            pointer :: CS
  type(verticalGrid_type), pointer :: GV
!   This subroutine sets information in the verticalGrid_type to be
! consistent with the sue of ALE mode
  integer :: nk

  nk = GV%ke
  GV%sInterface(1:nk+1) = getCoordinateInterfaces( CS%regridCS )
  GV%sLayer(1:nk) = 0.5*( GV%sInterface(1:nk) + GV%sInterface(2:nk+1) )
  GV%zAxisUnits = getCoordinateUnits( CS%regridCS )
  GV%zAxisLongName = getCoordinateShortName( CS%regridCS )
  GV%direction = -1 ! Because of ferret in z* mode. Need method to set
                    ! as function of coordinae mode.

end subroutine ALE_updateVerticalGridType

!------------------------------------------------------------------------------
! Write the vertical coordinate information into a file
!------------------------------------------------------------------------------
subroutine ALE_writeCoordinateFile( CS, G, directory )
  type(ALE_CS),          pointer       :: CS
  type(ocean_grid_type), intent(inout) :: G
  character(len=*),      intent(in)    :: directory
!   This subroutine writes out a file containing any available data related
! to the vertical grid used by the MOM ocean model when in ALE mode.
  character(len=120) :: filepath
  type(vardesc) :: vars(2)
  type(fieldtype) :: fields(2)
  integer :: unit
  real :: ds(G%ke), dsi(G%ke+1)

  filepath = trim(directory) // trim("Vertical_coordinate")
  ds(:) = getCoordinateResolution( CS%regridCS )
  dsi(1) = 0.5*ds(1)
  dsi(2:G%ke) = 0.5*( ds(1:G%ke-1) + ds(2:G%ke) )
  dsi(G%ke+1) = 0.5*ds(G%ke)

  vars(1) = vardesc('ds','Layer Coordinate Thickness','1','L','1', &
                    getCoordinateUnits( CS%regridCS ) )
  vars(2) = vardesc('ds_interface','Layer Center Coordinate Separation','1','i','1', &
                    getCoordinateUnits( CS%regridCS ) )

  call create_file(unit, trim(filepath), vars, 2, G, fields, SINGLE_FILE)
  call write_field(unit, fields(1), ds)
  call write_field(unit, fields(2), dsi)
  call close_file(unit)

end subroutine ALE_writeCoordinateFile

!------------------------------------------------------------------------------
! Set h to coordinate values for fixed coordinate systems
!------------------------------------------------------------------------------
subroutine ALE_initThicknessToCoord( CS, G, h )
  ! Arguments
  type(ALE_CS), intent(inout) :: CS ! Regridding parameters and options
  type(ocean_grid_type), intent(in)                  :: G
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(out) :: h ! Three-dimensional ocean grid
  ! Local variables
  integer :: i, j, k

  do j = G%jsd,G%jed ; do i = G%isd,G%ied
    h(i,j,:) = getStaticThickness( CS%regridCS, 0., G%bathyT(i,j) )
  enddo; enddo

end subroutine ALE_initThicknessToCoord

end module MOM_ALE
