!> Generates vertical grids as part of the ALE algorithm
module MOM_regridding

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : param_file_type, get_param, log_param
use MOM_io,            only : file_exists, field_exists, field_size, MOM_read_data
use MOM_io,            only : vardesc, var_desc, fieldtype, SINGLE_FILE
use MOM_io,            only : create_file, write_field, close_file, slasher
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS,           only : EOS_type, calculate_density, calculate_density_derivs
use MOM_EOS,           only : calculate_compress
use MOM_string_functions,only : uppercase, extractWord, extract_integer, extract_real

use regrid_edge_values, only : edge_values_explicit_h2, edge_values_explicit_h4
use regrid_edge_values, only : edge_values_implicit_h4, edge_values_implicit_h6
use regrid_edge_slopes, only : edge_slopes_implicit_h3, edge_slopes_implicit_h5

use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation
use PQM_functions, only : PQM_reconstruction, PQM_boundary_extrapolation_v1

use P1M_functions, only : P1M_interpolation, P1M_boundary_extrapolation
use P3M_functions, only : P3M_interpolation, P3M_boundary_extrapolation
use MOM_remapping, only : remapping_core_h
use MOM_remapping, only : remapping_CS
use regrid_consts, only : state_dependent, coordinateUnits
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA
use regrid_consts, only : REGRIDDING_ARBITRARY, REGRIDDING_SIGMA_SHELF_ZSTAR
use regrid_consts, only : REGRIDDING_HYCOM1, REGRIDDING_SLIGHT

use netcdf ! Used by check_grid_def()

implicit none ; private

#include <MOM_memory.h>

!> Regridding control structure
type, public :: regridding_CS
  private

  !> This array is set by function setCoordinateResolution()
  !! It contains the "resolution" or delta coordinate of the target
  !! coorindate. It has the units of the target coordiante, e.g.
  !! meters for z*, non-dimensional for sigma, etc.
  real, dimension(:), allocatable :: coordinateResolution

  !> This array is set by function set_target_densities()
  !! This array is the nominal coordinate of interfaces and is the
  !! running sum of coordinateResolution. i.e.
  !!  target_density(k+1) = coordinateResolution(k) + coordinateResolution(k)
  !! It is only used in "rho" mode.
  real, dimension(:), allocatable :: target_density
  logical :: target_density_set = .false.

  !> This array is set by function set_regrid_max_depths()
  !! It specifies the maximum depth that every interface is allowed to take, in H.
  real, dimension(:), allocatable :: max_interface_depths

  !> This array is set by function set_regrid_max_thickness()
  !! It specifies the maximum depth that every interface is allowed to take, in H.
  real, dimension(:), allocatable :: max_layer_thickness

  integer :: nk !< Number of layers/levels

  integer :: degree_i=4 !< Degree of interpolation polynomial

  !> Indicates which grid to use in the vertical (z*, sigma, target interface
  !! densities)
  integer :: regridding_scheme

  !> The following parameter is only relevant when used with the target
  !! interface densities regridding scheme. It indicates which interpolation
  !! to use to determine the grid.
  integer :: interpolation_scheme

  ! Indicate whether high-order boundary extrapolation should be used within
  !! boundary cells
  logical :: boundary_extrapolation

  !> Minimum thickness allowed when building the new grid through regridding
  real :: min_thickness

  !> Reference pressure for potential density calculations (Pa)
  real :: ref_pressure = 2.e7

  ! The following 4 parameters were introduced for use with the SLight coordinate:
  !> Depth over which to average to determine the mixed layer potential density (m)
  real :: Rho_ML_avg_depth = 1.0

  !> Number of layers to offset the mixed layer density to find resolved stratification (nondim)
  real :: nlay_ml_offset = 2.0

  !> The number of fixed-thickess layers at the top of the model
  integer :: nz_fixed_surface = 2

  !> The fixed resolution in the topmost SLight_nkml_min layers (m)
  real :: dz_ml_min = 1.0

  !> Weight given to old coordinate when blending between new and old grids (nondim)
  !! Used only below depth_of_time_filter_shallow, with a cubic variation
  !! from zero to full effect between depth_of_time_filter_shallow and
  !! depth_of_time_filter_deep.
  real :: old_grid_weight = 0.

  !> Depth above which no time-filtering of grid is applied (H units)
  real :: depth_of_time_filter_shallow = 0.

  !> Depth below which time-filtering of grid is applied at full effect (H units)
  real :: depth_of_time_filter_deep = 0.

  !> Fraction (between 0 and 1) of compressibility to add to potential density
  !! profiles when interpolating for target grid positions. (nondim)
  real :: compressibility_fraction = 0.

  !> If true, detect regions with much weaker stratification in the coordinate
  !! than based on in-situ density, and use a stretched coordinate there.
  logical :: fix_haloclines = .false.

  !> A length scale over which to filter T & S when looking for spuriously
  !! unstable water mass profiles, in m.
  real :: halocline_filter_length = 2.0

  !> A value of the stratification ratio that defines a problematic halocline region.
  real :: halocline_strat_tol = 0.25

  !> If true, each interface is given a maximum depth based on a rescaling of
  !! the indexing of coordinateResolution.
  logical :: set_maximum_depths = .false.

  !> A scaling factor (> 1) of the rate at which the coordinateResolution list
  !! is traversed to set the minimum depth of interfaces.
  real :: max_depth_index_scale = 2.0

  !> If true, integrate for interface positions from the top downward.
  !! If false, integrate from the bottom upward, as does the rest of the model.
  logical :: integrate_downward_for_e = .true.

end type

! The following routines are visible to the outside world
public initialize_regridding, end_regridding, regridding_main
public inflate_vanished_layers_old, check_remapping_grid, check_grid_column
public adjust_interface_motion
public set_regrid_params, get_regrid_size
public uniformResolution, setCoordinateResolution
public build_zstar_column, build_rho_column, build_sigma_column
public set_target_densities_from_GV, set_target_densities
public set_regrid_max_depths, set_regrid_max_thickness
public getCoordinateResolution, getCoordinateInterfaces
public getCoordinateUnits, getCoordinateShortName, getStaticThickness
public DEFAULT_COORDINATE_MODE

!> Documentation for coordinate options
character(len=322), parameter, public :: regriddingCoordinateModeDoc = &
                 " LAYER - Isopycnal or stacked shallow water layers\n"//&
                 " ZSTAR, Z* - stetched geopotential z*\n"//&
                 " SIGMA_SHELF_ZSTAR - stetched geopotential z* ignoring shelf\n"//&
                 " SIGMA - terrain following coordinates\n"//&
                 " RHO   - continuous isopycnal\n"//&
                 " HYCOM1 - HyCOM-like hybrid coordinate\n"//&
                 " SLIGHT - stretched coordinates above continuous isopycnal"

! Documentation for regridding interpolation schemes
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

! The following are private constants

! List of interpolation schemes
integer, parameter :: INTERPOLATION_P1M_H2     = 0 !< O(h^2)
integer, parameter :: INTERPOLATION_P1M_H4     = 1 !< O(h^2)
integer, parameter :: INTERPOLATION_P1M_IH4    = 2 !< O(h^2)
integer, parameter :: INTERPOLATION_PLM        = 3 !< O(h^2)
integer, parameter :: INTERPOLATION_PPM_H4     = 4 !< O(h^3)
integer, parameter :: INTERPOLATION_PPM_IH4    = 5 !< O(h^3)
integer, parameter :: INTERPOLATION_P3M_IH4IH3 = 6 !< O(h^4)
integer, parameter :: INTERPOLATION_P3M_IH6IH5 = 7 !< O(h^4)
integer, parameter :: INTERPOLATION_PQM_IH4IH3 = 8 !< O(h^4)
integer, parameter :: INTERPOLATION_PQM_IH6IH5 = 9 !< O(h^5)

!> List of interpolant degrees
integer, parameter :: DEGREE_1 = 1, DEGREE_2 = 2, DEGREE_3 = 3, DEGREE_4 = 4

!> Maximum number of regridding iterations
integer, parameter :: NB_REGRIDDING_ITERATIONS = 1
!> Deviation tolerance between succesive grids in regridding iterations
real, parameter    :: DEVIATION_TOLERANCE = 1e-10
!> Maximum number of Newton-Raphson iterations. Newton-Raphson iterations are
!! used to build the new grid by finding the coordinates associated with
!! target densities and interpolations of degree larger than 1.
integer, parameter :: NR_ITERATIONS = 8
!> Tolerance for Newton-Raphson iterations (stop when increment falls below this)
real, parameter    :: NR_TOLERANCE = 1e-12
!> When the N-R algorithm produces an estimate that lies outside [0,1], the
!! estimate is set to be equal to the boundary location, 0 or 1, plus or minus
!! an offset, respectively, when the derivative is zero at the boundary.
real, parameter    :: NR_OFFSET = 1e-6

! This CPP macro embeds some safety checks
#undef __DO_SAFETY_CHECKS__

contains

!> Initialization and configures a regridding control structure based on customizable run-time parameters
subroutine initialize_regridding(CS, GV, max_depth, param_file, mod, coord_mode, param_prefix, param_suffix)
  type(regridding_CS),        intent(inout) :: CS !< Regridding control structure
  type(verticalGrid_type),    intent(in)    :: GV         !< Ocean vertical grid structure
  real,                       intent(in)    :: max_depth  !< The maximum depth of the ocean, in m.
  type(param_file_type),      intent(in)    :: param_file !< Parameter file
  character(len=*),           intent(in)    :: mod        !< Name of calling module.
  character(len=*),           intent(in)    :: coord_mode !< Coordinate mode
  character(len=*),           intent(in)    :: param_prefix !< String to prefix to parameter names.
                                                            !! If empty, causes main model parameters to be used.
  character(len=*),           intent(in)    :: param_suffix !< String to append to parameter names.
  ! Local variables
  integer :: ke ! Number of levels
  character(len=80)  :: string, string2, varName ! Temporary strings
  character(len=40)  :: coord_units, param_name, coord_res_param ! Temporary strings
  character(len=200) :: inputdir, fileName
  character(len=320) :: message ! Temporary strings
  character(len=12) :: expected_units ! Temporary strings
  logical :: tmpLogical, fix_haloclines, set_max, do_sum, main_parameters
  logical :: coord_is_state_dependent, ierr
  real :: filt_len, strat_tol, index_scale, tmpReal
  real :: dz_fixed_sfc, Rho_avg_depth, nlay_sfc_int
  integer :: nz_fixed_sfc, k, nzf(4)
  real, dimension(:), allocatable :: dz     ! Resolution (thickness) in units of coordinate
  real, dimension(:), allocatable :: h_max  ! Maximum layer thicknesses, in m.
  real, dimension(:), allocatable :: dz_max ! Thicknesses used to find maximum interface depths, in m.
  real, dimension(:), allocatable :: z_max  ! Maximum interface depths, in m.
  real, dimension(:), allocatable :: rho_target ! Target density used in HYBRID mode
  ! Thicknesses that give level centers corresponding to table 2 of WOA09
  real, dimension(40) :: woa09_dz = (/ 5.,  10.,  10.,  15.,  22.5, 25., 25.,  25.,  &
                                      37.5, 50.,  50.,  75., 100., 100., 100., 100., &
                                     100., 100., 100., 100., 100., 100., 100., 175., &
                                     250., 375., 500., 500., 500., 500., 500., 500., &
                                     500., 500., 500., 500., 500., 500., 500., 500. /)

  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  main_parameters=.false.
  if (len_trim(param_prefix)==0) main_parameters=.true.
  if (main_parameters .and. len_trim(param_suffix)>0) call MOM_error(FATAL,trim(mod)//', initialize_regridding: '// &
              'Suffix provided without prefix for parameter names!')

  CS%nk = 0
  CS%regridding_scheme = coordinateMode(coord_mode)
  coord_is_state_dependent = state_dependent(coord_mode)

  if (main_parameters) then
    ! Read coordinate units parameter (main model = REGRIDDING_COORDINATE_UNITS)
    call get_param(param_file, mod, "REGRIDDING_COORDINATE_UNITS", coord_units, &
                 "Units of the regridding coordinuate.",&
                 default=coordinateUnits(coord_mode))
  else
    coord_units=coordinateUnits(coord_mode)
  endif

  if (coord_is_state_dependent) then
    if (main_parameters) then
      param_name = "INTERPOLATION_SCHEME"
      string2 = regriddingDefaultInterpScheme
    else
      param_name = trim(param_prefix)//"_INTERP_SCHEME_"//trim(param_suffix)
      string2 = 'PPM_H4' ! Default for diagnostics
    endif
    call get_param(param_file, mod, "INTERPOLATION_SCHEME", string, &
                 "This sets the interpolation scheme to use to\n"//&
                 "determine the new grid. These parameters are\n"//&
                 "only relevant when REGRIDDING_COORDINATE_MODE is\n"//&
                 "set to a function of state. Otherwise, it is not\n"//&
                 "used. It can be one of the following schemes:\n"//&
                 trim(regriddingInterpSchemeDoc), default=trim(string2))
    call set_regrid_params(CS, interp_scheme=string)
  else
    CS%interpolation_scheme = -1 ! Cause error if ever used
  endif

  if (main_parameters .and. coord_is_state_dependent) then
    call get_param(param_file, mod, "REGRID_COMPRESSIBILITY_FRACTION", tmpReal, &
                 "When interpolating potential density profiles we can add\n"//&
                 "some artificial compressibility solely to make homogenous\n"//&
                 "regions appear stratified.", default=0.)
    call set_regrid_params(CS, compress_fraction=tmpReal)
  endif

  ! Read coordinate configuration parameter (main model = ALE_COORDINATE_CONFIG)
  if (main_parameters) then
    param_name = "ALE_COORDINATE_CONFIG"
    coord_res_param = "ALE_RESOLUTION"
    string2 = 'UNIFORM'
  else
    param_name = trim(param_prefix)//"_DEF_"//trim(param_suffix)
    coord_res_param = trim(param_prefix)//"_RES_"//trim(param_suffix)
    string2 = 'UNIFORM'
    if (max_depth>3000.) string2='WOA09' ! For convenience
  endif
  call get_param(param_file, mod, param_name, string, &
                 "Determines how to specify the coordinate\n"//&
                 "resolution. Valid options are:\n"//&
                 " PARAM       - use the vector-parameter "//trim(coord_res_param)//"\n"//&
                 " UNIFORM[:N] - uniformly distributed\n"//&
                 " FILE:string - read from a file. The string specifies\n"//&
                 "               the filename and variable name, separated\n"//&
                 "               by a comma or space, e.g. FILE:lev.nc,dz\n"//&
                 "               or FILE:lev.nc,interfaces=zw\n"//&
                 " WOA09[:N]   - the WOA09 vertical grid (approximately)\n"//&
                 " FNC1:string - FNC1:dz_min,H_total,power,precision\n"//&
                 " HYBRID:string - read from a file. The string specifies\n"//&
                 "               the filename and two variable names, separated\n"//&
                 "               by a comma or space, for sigma-2 and dz. e.g.\n"//&
                 "               HYBRID:vgrid.nc,sigma2,dz",&
                 default=trim(string2))
  message = "The distribution of vertical resolution for the target\n"//&
            "grid used for Eulerian-like coordinates. For example,\n"//&
            "in z-coordinate mode, the parameter is a list of level\n"//&
            "thicknesses (in m). In sigma-coordinate mode, the list\n"//&
            "is of non-dimensional fractions of the water column."
  if (index(trim(string),'UNIFORM')==1) then
    if (len_trim(string)==7) then
      ke = GV%ke ! Use model nk by default
      tmpReal = max_depth
    elseif (index(trim(string),'UNIFORM:')==1 .and. len_trim(string)>8) then
      ! Format is "UNIFORM:N" or "UNIFORM:N,dz"
      ke = extract_integer(string(9:len_trim(string)),'',1)
      tmpReal = extract_real(string(9:len_trim(string)),',',2,missing_value=max_depth)
    else
      call MOM_error(FATAL,trim(mod)//', initialize_regridding: '// &
          'Unable to interpret "'//trim(string)//'".')
    endif
    allocate(dz(ke))
    if (ke==1) then
      dz(:) = uniformResolution(ke, coord_mode, tmpReal, GV%Rlay(1), GV%Rlay(1))
    else
      dz(:) = uniformResolution(ke, coord_mode, tmpReal, &
                   GV%Rlay(1)+0.5*(GV%Rlay(1)-GV%Rlay(2)), &
                   GV%Rlay(ke)+0.5*(GV%Rlay(ke)-GV%Rlay(ke-1)) )
    endif
    if (main_parameters) call log_param(param_file, mod, "!"//coord_res_param, dz, &
                   trim(message), units=trim(coord_units))
  elseif (trim(string)=='PARAM') then
    ! Read coordinate resolution (main model = ALE_RESOLUTION)
    ke = GV%ke ! Use model nk by default
    allocate(dz(ke))
    call get_param(param_file, mod, coord_res_param, dz, &
                   trim(message), units=trim(coord_units), fail_if_missing=.true.)
  elseif (index(trim(string),'FILE:')==1) then
    ! FILE:filename,var_name is assumed to be reading level thickness variables
    ! FILE:filename,interfaces=var_name reads positions
    if (string(6:6)=='.' .or. string(6:6)=='/') then
      ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
      fileName = trim( extractWord(trim(string(6:80)), 1) )
    else
      ! Otherwise assume we should look for the file in INPUTDIR
      fileName = trim(inputdir) // trim( extractWord(trim(string(6:80)), 1) )
    endif
    if (.not. file_exists(fileName)) call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
            "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

    varName = trim( extractWord(trim(string(6:)), 2) )
    if (len_trim(varName)==0) then
      if (field_exists(fileName,'dz')) then; varName = 'dz'
      elseif (field_exists(fileName,'dsigma')) then; varName = 'dsigma'
      elseif (field_exists(fileName,'ztest')) then; varName = 'ztest'
      else ;  call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
                    "Coordinate variable not specified and none could be guessed.")
      endif
    endif
    ! This check fails when the variable is a dimension variable! -AJA
   !if (.not. field_exists(fileName,trim(varName))) call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
   !             "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
    if (CS%regridding_scheme == REGRIDDING_SIGMA) then
      expected_units = 'nondim'
    elseif (CS%regridding_scheme == REGRIDDING_RHO) then
      expected_units = 'kg m-3'
    else
      expected_units = 'meters'
    endif
    if (index(trim(varName),'interfaces=')==1) then
      varName=trim(varName(12:))
      call check_grid_def(filename, varName, expected_units, message, ierr)
      if (ierr) call MOM_error(FATAL,trim(mod)//", initialize_regridding: "//&
                  "Unsupported format in grid definition '"//trim(filename)//"'. Error message "//trim(message))
      call field_size(trim(fileName), trim(varName), nzf)
      ke = nzf(1)-1
      allocate(dz(ke))
      allocate(z_max(ke+1))
      call MOM_read_data(trim(fileName), trim(varName), z_max)
      dz(:) = abs(z_max(1:ke) - z_max(2:ke+1))
      deallocate(z_max)
    else
      ! Assume reading resolution
      call field_size(trim(fileName), trim(varName), nzf)
      ke = nzf(1)
      allocate(dz(ke))
      call MOM_read_data(trim(fileName), trim(varName), dz)
    endif
    if (main_parameters .and. ke/=GV%ke) then
      call MOM_error(FATAL,trim(mod)//', initialize_regridding: '// &
                 'Mismatch in number of model levels and "'//trim(string)//'".')
    endif
    if (main_parameters) call log_param(param_file, mod, "!"//coord_res_param, dz, &
               trim(message), units=coordinateUnits(coord_mode))
  elseif (index(trim(string),'FNC1:')==1) then
    ke = GV%ke; allocate(dz(ke))
    call dz_function1( trim(string(6:)), dz )
    if (main_parameters) call log_param(param_file, mod, "!"//coord_res_param, dz, &
               trim(message), units=coordinateUnits(coord_mode))
  elseif (index(trim(string),'HYBRID:')==1) then
    ke = GV%ke; allocate(dz(ke))
    ! The following assumes the FILE: syntax of above but without "FILE:" in the string
    allocate(rho_target(ke+1))
    fileName = trim( extractWord(trim(string(8:)), 1) )
    if (fileName(1:1)/='.' .and. filename(1:1)/='/') fileName = trim(inputdir) // trim( fileName )
    if (.not. file_exists(fileName)) call MOM_error(FATAL,trim(mod)//", initialize_regridding: HYBRID "// &
      "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")
    varName = trim( extractWord(trim(string(8:)), 2) )
    if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,trim(mod)//", initialize_regridding: HYBRID "// &
      "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
    call MOM_read_data(trim(fileName), trim(varName), rho_target)
    varName = trim( extractWord(trim(string(8:)), 3) )
    if (varName(1:5) == 'FNC1:') then ! Use FNC1 to calculate dz
      call dz_function1( trim(string((index(trim(string),'FNC1:')+5):)), dz )
    else ! Read dz from file
      if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,trim(mod)//", initialize_regridding: HYBRID "// &
        "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
      call MOM_read_data(trim(fileName), trim(varName), dz)
    endif
    if (main_parameters) then
      call log_param(param_file, mod, "!"//coord_res_param, dz, &
               trim(message), units=coordinateUnits(coord_mode))
      call log_param(param_file, mod, "!TARGET_DENSITIES", rho_target, &
               'HYBRID target densities for itnerfaces', units=coordinateUnits(coord_mode))
    endif
  elseif (index(trim(string),'WOA09')==1) then
    if (len_trim(string)==5) then
      tmpReal = 0. ; ke = 0
      do while (tmpReal<max_depth)
        ke = ke + 1
        tmpReal = tmpReal + woa09_dz(ke)
      enddo
    elseif (index(trim(string),'WOA09:')==1) then
      if (len_trim(string)==6) call MOM_error(FATAL,trim(mod)//', initialize_regridding: '// &
                 'Expected string of form "WOA09:N" but got "'//trim(string)//'".')
      ke = extract_integer(string(7:len_trim(string)),'',1)
    endif
    if (ke>40 .or. ke<1) call MOM_error(FATAL,trim(mod)//', initialize_regridding: '// &
                 'For "WOA05:N" N must 0<N<41 but got "'//trim(string)//'".')
    allocate(dz(ke))
    dz(1:ke) = woa09_dz(1:ke)
  else
    call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
      "Unrecognized coordinate configuraiton"//trim(string))
  endif

  if (main_parameters) then
    ! This is a work around to apparently needed to work with the from_Z initialization...  ???
    if (coordinateMode(coord_mode) == REGRIDDING_ZSTAR .or. &
        coordinateMode(coord_mode) == REGRIDDING_HYCOM1 .or. &
        coordinateMode(coord_mode) == REGRIDDING_SLIGHT) then
      ! Adjust target grid to be consistent with max_depth
      tmpReal = sum( dz(:) )
      if (tmpReal < max_depth) then
        dz(ke) = dz(ke) + ( max_depth - tmpReal )
      elseif (tmpReal > max_depth) then
        if ( dz(ke) + ( max_depth - tmpReal ) > 0. ) then
          dz(ke) = dz(ke) + ( max_depth - tmpReal )
        else
          call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
            "MAXIMUM_DEPTH was too shallow to adjust bottom layer of DZ!"//trim(string))
        endif
      endif
    endif
  endif

  CS%nk=ke
  ! Target resolution (for fixed coordinates)
  allocate( CS%coordinateResolution(CS%nk) ); CS%coordinateResolution(:) = -1.E30
  if (state_dependent(CS%regridding_scheme)) then
    ! Target values
    allocate( CS%target_density(CS%nk+1) ); CS%target_density(:) = -1.E30
  endif

  call setCoordinateResolution(dz, CS)
  if (allocated(rho_target)) then
    call set_target_densities(CS, rho_target)
    deallocate(rho_target)
  endif
  ! \todo This line looks like it would overwrite the target densities set just above?
  if (coordinateMode(coord_mode) == REGRIDDING_RHO) call set_target_densities_from_GV(GV, CS)

  if (main_parameters) then
    call get_param(param_file, mod, "MIN_THICKNESS", tmpReal, &
                 "When regridding, this is the minimum layer\n"//&
                 "thickness allowed.", units="m",&
                 default=regriddingDefaultMinThickness )
    call set_regrid_params(CS, min_thickness=tmpReal)
  else
    call set_regrid_params(CS, min_thickness=0.)
  endif

  if (main_parameters .and. coord_is_state_dependent) then
    call get_param(param_file, mod, "BOUNDARY_EXTRAPOLATION", tmpLogical, &
                 "When defined, a proper high-order reconstruction\n"//&
                 "scheme is used within boundary cells rather\n"//&
                 "than PCM. E.g., if PPM is used for remapping, a\n"//&
                 "PPM reconstruction will also be used within\n"//&
                 "boundary cells.", default=regriddingDefaultBoundaryExtrapolation)
    call set_regrid_params(CS, boundary_extrapolation=tmpLogical)
  else
    call set_regrid_params(CS, boundary_extrapolation=.false.)
  endif

  if (coordinateMode(coord_mode) == REGRIDDING_SLIGHT) then
    ! Set SLight-specific regridding parameters.
    call get_param(param_file, mod, "SLIGHT_DZ_SURFACE", dz_fixed_sfc, &
                 "The nominal thickness of fixed thickness near-surface\n"//&
                 "layers with the SLight coordinate.", units="m", default=1.0)
    call get_param(param_file, mod, "SLIGHT_NZ_SURFACE_FIXED", nz_fixed_sfc, &
                 "The number of fixed-depth surface layers with the SLight\n"//&
                 "coordinate.", units="nondimensional", default=2)
    call get_param(param_file, mod, "SLIGHT_SURFACE_AVG_DEPTH", Rho_avg_depth, &
                 "The thickness of the surface region over which to average\n"//&
                 "when calculating the density to use to define the interior\n"//&
                 "with the SLight coordinate.", units="m", default=1.0)
    call get_param(param_file, mod, "SLIGHT_NLAY_TO_INTERIOR", nlay_sfc_int, &
                 "The number of layers to offset the surface density when\n"//&
                 "defining where the interior ocean starts with SLight.", &
                 units="nondimensional", default=2.0)
    call get_param(param_file, mod, "SLIGHT_FIX_HALOCLINES", fix_haloclines, &
                 "If true, identify regions above the reference pressure\n"//&
                 "where the reference pressure systematically underestimates\n"//&
                 "the stratification and use this in the definition of the\n"//&
                 "interior with the SLight coordinate.", default=.false.)

    call set_regrid_params(CS, dz_min_surface=dz_fixed_sfc, &
                nz_fixed_surface=nz_fixed_sfc, Rho_ML_avg_depth=Rho_avg_depth, &
                nlay_ML_to_interior=nlay_sfc_int, fix_haloclines=fix_haloclines)
    if (fix_haloclines) then
      ! Set additional parameters related to SLIGHT_FIX_HALOCLINES.
      call get_param(param_file, mod, "HALOCLINE_FILTER_LENGTH", filt_len, &
                 "A length scale over which to smooth the temperature and\n"//&
                 "salinity before identifying erroneously unstable haloclines.", &
                 units="m", default=2.0)
      call get_param(param_file, mod, "HALOCLINE_STRAT_TOL", strat_tol, &
                 "A tolerance for the ratio of the stratification of the\n"//&
                 "apparent coordinate stratification to the actual value\n"//&
                 "that is used to identify erroneously unstable haloclines.\n"//&
                 "This ratio is 1 when they are equal, and sensible values \n"//&
                 "are between 0 and 0.5.", units="nondimensional", default=0.2)
      call set_regrid_params(CS, halocline_filt_len=filt_len, &
                             halocline_strat_tol=strat_tol)
    endif

  endif

  if (main_parameters .and. coord_is_state_dependent) then
    call get_param(param_file, mod, "MAXIMUM_INT_DEPTH_CONFIG", string, &
                 "Determines how to specify the maximum interface depths.\n"//&
                 "Valid options are:\n"//&
                 " NONE        - there are no maximum interface depths\n"//&
                 " PARAM       - use the vector-parameter MAXIMUM_INTERFACE_DEPTHS\n"//&
                 " FILE:string - read from a file. The string specifies\n"//&
                 "               the filename and variable name, separated\n"//&
                 "               by a comma or space, e.g. FILE:lev.nc,Z\n"//&
                 " FNC1:string - FNC1:dz_min,H_total,power,precision",&
                 default='NONE')
    message = "The list of maximum depths for each interface."
    allocate(z_max(ke+1))
    allocate(dz_max(ke))
    if ( trim(string) == "NONE") then
      ! Do nothing.
    elseif ( trim(string) ==  "PARAM") then
      call get_param(param_file, mod, "MAXIMUM_INTERFACE_DEPTHS", z_max, &
                   trim(message), units="m", fail_if_missing=.true.)
      call set_regrid_max_depths(CS, z_max, GV%m_to_H)
    elseif (index(trim(string),'FILE:')==1) then
      if (string(6:6)=='.' .or. string(6:6)=='/') then
        ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
        fileName = trim( extractWord(trim(string(6:80)), 1) )
      else
        ! Otherwise assume we should look for the file in INPUTDIR
        fileName = trim(inputdir) // trim( extractWord(trim(string(6:80)), 1) )
      endif
      if (.not. file_exists(fileName)) call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
        "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

      do_sum = .false.
      varName = trim( extractWord(trim(string(6:)), 2) )
      if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
        "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
      if (len_trim(varName)==0) then
        if (field_exists(fileName,'z_max')) then; varName = 'z_max'
        elseif (field_exists(fileName,'dz')) then; varName = 'dz' ; do_sum = .true.
        elseif (field_exists(fileName,'dz_max')) then; varName = 'dz_max' ; do_sum = .true.
        else ; call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
          "MAXIMUM_INT_DEPTHS variable not specified and none could be guessed.")
        endif
      endif
      if (do_sum) then
        call MOM_read_data(trim(fileName), trim(varName), dz_max)
        z_max(1) = 0.0 ; do K=1,ke ; z_max(K+1) = z_max(K) + dz_max(k) ; enddo
      else
        call MOM_read_data(trim(fileName), trim(varName), z_max)
      endif
      call log_param(param_file, mod, "!MAXIMUM_INT_DEPTHS", z_max, &
                 trim(message), units=coordinateUnits(coord_mode))
      call set_regrid_max_depths(CS, z_max, GV%m_to_H)
    elseif (index(trim(string),'FNC1:')==1) then
      call dz_function1( trim(string(6:)), dz_max )
      if ((coordinateMode(coord_mode) == REGRIDDING_SLIGHT) .and. &
          (dz_fixed_sfc > 0.0)) then
        do k=1,nz_fixed_sfc ; dz_max(k) = dz_fixed_sfc ; enddo
      endif
      z_max(1) = 0.0 ; do K=1,ke ; z_max(K+1) = z_max(K) + dz_max(K) ; enddo
      call log_param(param_file, mod, "!MAXIMUM_INT_DEPTHS", z_max, &
                 trim(message), units=coordinateUnits(coord_mode))
      call set_regrid_max_depths(CS, z_max, GV%m_to_H)
    else
      call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
        "Unrecognized MAXIMUM_INT_DEPTH_CONFIG "//trim(string))
    endif
    deallocate(z_max)
    deallocate(dz_max)

    ! Optionally specify maximum thicknesses for each layer, enforced by moving
    ! the interface below a layer downward.
    call get_param(param_file, mod, "MAX_LAYER_THICKNESS_CONFIG", string, &
                   "Determines how to specify the maximum layer thicknesses.\n"//&
                   "Valid options are:\n"//&
                   " NONE        - there are no maximum layer thicknesses\n"//&
                   " PARAM       - use the vector-parameter MAX_LAYER_THICKNESS\n"//&
                   " FILE:string - read from a file. The string specifies\n"//&
                   "               the filename and variable name, separated\n"//&
                   "               by a comma or space, e.g. FILE:lev.nc,Z\n"//&
                   " FNC1:string - FNC1:dz_min,H_total,power,precision",&
                   default='NONE')
    message = "The list of maximum thickness for each layer."
    allocate(h_max(ke))
    if ( trim(string) == "NONE") then
      ! Do nothing.
    elseif ( trim(string) ==  "PARAM") then
      call get_param(param_file, mod, "MAX_LAYER_THICKNESS", h_max, &
                   trim(message), units="m", fail_if_missing=.true.)
      call set_regrid_max_thickness(CS, h_max, GV%m_to_H)
    elseif (index(trim(string),'FILE:')==1) then
      if (string(6:6)=='.' .or. string(6:6)=='/') then
        ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
        fileName = trim( extractWord(trim(string(6:80)), 1) )
      else
        ! Otherwise assume we should look for the file in INPUTDIR
        fileName = trim(inputdir) // trim( extractWord(trim(string(6:80)), 1) )
      endif
      if (.not. file_exists(fileName)) call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
        "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

      varName = trim( extractWord(trim(string(6:)), 2) )
      if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
        "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
      if (len_trim(varName)==0) then
        if (field_exists(fileName,'h_max')) then; varName = 'h_max'
        elseif (field_exists(fileName,'dz_max')) then; varName = 'dz_max'
        else ; call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
          "MAXIMUM_INT_DEPTHS variable not specified and none could be guessed.")
        endif
      endif
      call MOM_read_data(trim(fileName), trim(varName), h_max)
      call log_param(param_file, mod, "!MAX_LAYER_THICKNESS", h_max, &
                 trim(message), units=coordinateUnits(coord_mode))
      call set_regrid_max_thickness(CS, h_max, GV%m_to_H)
    elseif (index(trim(string),'FNC1:')==1) then
      call dz_function1( trim(string(6:)), h_max )
      call log_param(param_file, mod, "!MAX_LAYER_THICKNESS", h_max, &
                 trim(message), units=coordinateUnits(coord_mode))
      call set_regrid_max_thickness(CS, h_max, GV%m_to_H)
    else
      call MOM_error(FATAL,trim(mod)//", initialize_regridding: "// &
        "Unrecognized MAX_LAYER_THICKNESS_CONFIG "//trim(string))
    endif
    deallocate(h_max)
  endif

  deallocate(dz)
end subroutine initialize_regridding

!> Do some basic checks on the vertical grid definition file, variable
subroutine check_grid_def(filename, varname, expected_units, msg, ierr)
  character(len=*), intent(in)    :: filename !< File name
  character(len=*), intent(in)    :: varname !< Variable name
  character(len=*), intent(in)    :: expected_units !< Expected units of variable
  character(len=*), intent(inout) :: msg !< Message to use for errors
  logical,          intent(out)   :: ierr !< True if an error occurs
  ! Local variables
  character (len=200) :: units, long_name
  integer :: ncid, status, intid, vid

  ierr = .false.
  status = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid);
  if (status /= NF90_NOERR) then
    ierr = .true.
    msg = 'File not found: '//trim(filename)
    return
  endif

  status = NF90_INQ_VARID(ncid, trim(varname), vid)
  if (status /= NF90_NOERR) then
    ierr = .true.
    msg = 'Var not found: '//trim(varname)
    return
  endif

  status = NF90_GET_ATT(ncid, vid, "units", units)
  if (status /= NF90_NOERR) then
    ierr = .true.
    msg = 'Attribute not found: units'
    return
  endif

  if (trim(units) /= trim(expected_units)) then
    if (trim(expected_units) == "meters") then
      if (trim(units) /= "m") then
        ierr = .true.
      endif
    else
      ierr = .true.
    endif
  endif

  if (ierr) then
    msg = 'Units incorrect: '//trim(units)//' /= '//trim(expected_units)
  endif

end subroutine check_grid_def

!> Deallocation of regridding memory
subroutine end_regridding(CS)
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure

  deallocate( CS%coordinateResolution )
  if (allocated(CS%target_density)) deallocate( CS%target_density )
  if (allocated(CS%max_interface_depths) ) deallocate( CS%max_interface_depths )
  if (allocated(CS%max_layer_thickness) ) deallocate( CS%max_layer_thickness )

end subroutine end_regridding

!> Numeric value of interpolation_scheme corresponding to scheme name
integer function interpolation_scheme(interp_scheme)
  character(len=*), intent(in) :: interp_scheme !< Name of interpolation scheme

  select case ( uppercase(trim(interp_scheme)) )
    case ("P1M_H2");     interpolation_scheme = INTERPOLATION_P1M_H2
    case ("P1M_H4");     interpolation_scheme = INTERPOLATION_P1M_H4
    case ("P1M_IH2");    interpolation_scheme = INTERPOLATION_P1M_IH4
    case ("PLM");        interpolation_scheme = INTERPOLATION_PLM
    case ("PPM_H4");     interpolation_scheme = INTERPOLATION_PPM_H4
    case ("PPM_IH4");    interpolation_scheme = INTERPOLATION_PPM_IH4
    case ("P3M_IH4IH3"); interpolation_scheme = INTERPOLATION_P3M_IH4IH3
    case ("P3M_IH6IH5"); interpolation_scheme = INTERPOLATION_P3M_IH6IH5
    case ("PQM_IH4IH3"); interpolation_scheme = INTERPOLATION_PQM_IH4IH3
    case ("PQM_IH6IH5"); interpolation_scheme = INTERPOLATION_PQM_IH6IH5
    case default ; call MOM_error(FATAL, "MOM_regridding: "//&
     "Unrecognized choice for INTERPOLATION_SCHEME ("//trim(interp_scheme)//").")
  end select

end function interpolation_scheme

!------------------------------------------------------------------------------
! Dispatching regridding routine: regridding & remapping
!------------------------------------------------------------------------------
subroutine regridding_main( remapCS, CS, G, GV, h, tv, h_new, dzInterface, frac_shelf_h)
!------------------------------------------------------------------------------
! This routine takes care of (1) building a new grid and (2) remapping between
! the old grid and the new grid. The creation of the new grid can be based
! on z coordinates, target interface densities, sigma coordinates or any
! arbitrary coordinate system.
!   The MOM6 interface positions are always calculated from the bottom up by
! accumulating the layer thicknesses starting at z=-G%bathyT.  z increases
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
  type(remapping_CS),                         intent(in)    :: remapCS !< Remapping parameters and options
  type(regridding_CS),                        intent(in)    :: CS     !< Regridding control structure
  type(ocean_grid_type),                      intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV     !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)), intent(inout) :: h      !< Current 3D grid obtained after the last time step
  type(thermo_var_ptrs),                      intent(inout) :: tv     !< Thermodynamical variables (T, S, ...)
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)), intent(inout) :: h_new  !< New 3D grid consistent with target coordinate
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)+1), intent(inout) :: dzInterface !< The change in position of each interface
  real, dimension(:,:),                   optional, pointer :: frac_shelf_h !< Fractional ice shelf coverage
  ! Local variables
  real :: trickGnuCompiler
  logical :: use_ice_shelf

  use_ice_shelf = .false.
  if (present(frac_shelf_h)) then
    if (associated(frac_shelf_h)) use_ice_shelf = .true.
  endif

  select case ( CS%regridding_scheme )

    case ( REGRIDDING_ZSTAR )
      if (use_ice_shelf) then
        call build_zstar_grid( CS, G, GV, h, dzInterface, frac_shelf_h )
      else
        call build_zstar_grid( CS, G, GV, h, dzInterface )
      endif
      call calc_h_new_by_dz(G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_SIGMA_SHELF_ZSTAR)
      call build_zstar_grid( CS, G, GV, h, dzInterface )
      call calc_h_new_by_dz(G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_SIGMA )
      call build_sigma_grid( CS, G, GV, h, dzInterface )
      call calc_h_new_by_dz(G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_RHO )
      call convective_adjustment(G, GV, h, tv)
      call build_rho_grid( G, GV, h, tv, dzInterface, remapCS, CS )
      call calc_h_new_by_dz(G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_ARBITRARY )
      call build_grid_arbitrary( G, GV, h, dzInterface, trickGnuCompiler, CS )
      call calc_h_new_by_dz(G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_HYCOM1 )
      call build_grid_HyCOM1( G, GV, h, tv, dzInterface, remapCS, CS )
      call calc_h_new_by_dz(G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_SLIGHT )
      call build_grid_SLight( G, GV, h, tv, dzInterface, remapCS, CS )
      call calc_h_new_by_dz(G, GV, h, dzInterface, h_new)

    case default
      call MOM_error(FATAL,'MOM_regridding, regridding_main: '//&
                     'Unknown regridding scheme selected!')

  end select ! type of grid

#ifdef __DO_SAFETY_CHECKS__
  call check_remapping_grid(G, GV, h, dzInterface,'in regridding_main')
#endif

end subroutine regridding_main

!> Calculates h_new from h + delta_k dzInterface
subroutine calc_h_new_by_dz(G, GV, h, dzInterface, h_new)
  type(ocean_grid_type),                      intent(in)    :: G !< Grid structure
  type(verticalGrid_type),                    intent(in)    :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h !< Old layer thicknesses (m)
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in)    :: dzInterface !< Change in interface positions (m)
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(inout) :: h_new !< New layer thicknesses (m)
  ! Local variables
  integer :: i, j, k

!$OMP parallel do default(none) shared(G,GV,h,dzInterface,h_new)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1
      if (G%mask2dT(i,j)>0.) then
        do k=1,GV%ke
          h_new(i,j,k) = max( 0., h(i,j,k) + ( dzInterface(i,j,k) - dzInterface(i,j,k+1) ) )
        enddo
      else
        h_new(i,j,:) = h(i,j,:)
      endif
    enddo
  enddo

end subroutine calc_h_new_by_dz

!> Check that the total thickness of two grids match
subroutine check_remapping_grid( G, GV, h, dzInterface, msg )
  type(ocean_grid_type),                       intent(in) :: G !< Grid structure
  type(verticalGrid_type),                     intent(in) :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in) :: h !< Layer thicknesses (m)
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: dzInterface !< Change in interface positions (m)
  character(len=*),                            intent(in) :: msg !< Message to append to errors
  ! Local variables
  integer :: i, j

!$OMP parallel do default(none) shared(G,GV,h,dzInterface,msg)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1
      if (G%mask2dT(i,j)>0.) call check_grid_column( GV%ke, G%bathyT(i,j), h(i,j,:), dzInterface(i,j,:), msg )
    enddo
  enddo

end subroutine check_remapping_grid

!> Check that the total thickness of new and old grids are consistent
subroutine check_grid_column( nk, depth, h, dzInterface, msg )
  integer,               intent(in) :: nk !< Number of cells
  real,                  intent(in) :: depth !< Depth of bottom (m)
  real, dimension(nk),   intent(in) :: h !< Cell thicknesses (m)
  real, dimension(nk+1), intent(in) :: dzInterface !< Change in interface positions (m)
  character(len=*),      intent(in) :: msg !< Message to append to errors
  ! Local variables
  integer :: k
  real    :: eps, total_h_old, total_h_new, h_new, z_old, z_new

  eps =1. ; eps = epsilon(eps)

  ! Total thickness of grid h
  total_h_old = 0.
  do k = 1,nk
    total_h_old = total_h_old + h(k)
  enddo

  ! Integrate upwards for the interfaces consistent with the rest of MOM6
  z_old = - depth
  if (depth == 0.) z_old = - total_h_old
  total_h_new = 0.
  do k = nk,1,-1
    z_old = z_old + h(k) ! Old interface position above layer k
    z_new = z_old + dzInterface(k) ! New interface position based on dzInterface
    h_new = h(k) + ( dzInterface(k) - dzInterface(k+1) ) ! New thickness
    if (h_new<0.) then
      write(0,*) 'k,h,hnew=',k,h(k),h_new
      write(0,*) 'dzI(k+1),dzI(k)=',dzInterface(k+1),dzInterface(k)
      call MOM_error( FATAL, 'MOM_regridding, check_grid_column: '//&
        'Negative layer thickness implied by re-gridding, '//trim(msg))
    endif
    total_h_new = total_h_new + h_new

  enddo

  ! Conservation by implied h_new
  if (abs(total_h_new-total_h_old)>real(nk-1)*0.5*(total_h_old+total_h_new)*eps) then
    write(0,*) 'nk=',nk
    do k = 1,nk
      write(0,*) 'k,h,hnew=',k,h(k),h(k)+(dzInterface(k)-dzInterface(k+1))
    enddo
    write(0,*) 'Hold,Hnew,Hnew-Hold=',total_h_old,total_h_new,total_h_new-total_h_old
    write(0,*) 'eps,(n)/2*eps*H=',eps,real(nk-1)*0.5*(total_h_old+total_h_new)*eps
    call MOM_error( FATAL, 'MOM_regridding, check_grid_column: '//&
      'Re-gridding did NOT conserve total thickness to within roundoff '//trim(msg))
  endif

  ! Check that the top and bottom are intentionally moving
  if (dzInterface(1) /= 0.) call MOM_error( FATAL, &
    'MOM_regridding, check_grid_column: Non-zero dzInterface at surface! '//trim(msg))
  if (dzInterface(nk+1) /= 0.) call MOM_error( FATAL, &
    'MOM_regridding, check_grid_column: Non-zero dzInterface at bottom! '//trim(msg))

end subroutine check_grid_column

!> Returns the change in interface position motion after filtering and
!! assuming the top and bottom interfaces do not move.  The filtering is
!! a function of depth, and is applied as the integrated average filtering
!! over the trajectory of the interface.  By design, this code can not give
!! tangled interfaces provided that z_old and z_new are not already tangled.
subroutine filtered_grid_motion( CS, nk, z_old, z_new, dz_g )
  type(regridding_CS),                           intent(in)    :: CS !< Regridding control structure
  integer,               intent(in)    :: nk !< Number of cells
  real, dimension(nk+1), intent(in)    :: z_old !< Old grid position (m)
  real, dimension(nk+1), intent(in)    :: z_new !< New grid position (m)
  real, dimension(nk+1), intent(inout) :: dz_g !< Change in interface positions (m)
  ! Local variables
  real :: sgn  ! The sign convention for downward.
  real :: dz_tgt, zr1
  real :: Aq, Bq, dz0, z0, F0
  real :: zs, zd, dzwt, Idzwt
  real :: wtd, Iwtd
  real :: Int_zs, Int_zd, dInt_zs_zd
! For debugging:
  real, dimension(nk+1) :: z_act
!  real, dimension(nk+1) :: ddz_g_s, ddz_g_d
  logical :: debug = .false.
  integer :: k

  if ((z_old(nk+1) - z_old(1)) * (z_new(nk+1) - z_new(1)) < 0.0) then
    call MOM_error(FATAL, "filtered_grid_motion: z_old and z_new use different sign conventions.")
  elseif ((z_old(nk+1) - z_old(1)) * (z_new(nk+1) - z_new(1)) == 0.0) then
    ! This is a massless column, so do nothing and return.
    do k=1,nk+1 ; dz_g(k) = 0.0 ; enddo ; return
  elseif ((z_old(nk+1) - z_old(1)) + (z_new(nk+1) - z_new(1)) > 0.0) then
    sgn = 1.0
  else
    sgn = -1.0
  endif

  if (debug) then
    do k=2,nk+1
      if (sgn*(z_new(k)-z_new(k-1)) < -5e-16*(abs(z_new(k))+abs(z_new(k-1))) ) &
        call MOM_error(FATAL, "filtered_grid_motion: z_new is tangled.")
      if (sgn*(z_old(k)-z_old(k-1)) < -5e-16*(abs(z_old(k))+abs(z_old(k-1))) ) &
        call MOM_error(FATAL, "filtered_grid_motion: z_old is tangled.")
    enddo
    ! ddz_g_s(:) = 0.0 ; ddz_g_d(:) = 0.0
  endif

  zs = CS%depth_of_time_filter_shallow
  zd = CS%depth_of_time_filter_deep
  wtd = 1.0 - CS%old_grid_weight
  Iwtd = 1.0 / wtd

  dzwt = (zd - zs)
  Idzwt = 0.0 ; if (abs(zd - zs) > 0.0) Idzwt = 1.0 / (zd - zs)
  dInt_zs_zd = 0.5*(1.0 + Iwtd) * (zd - zs)
  Aq = 0.5*(Iwtd - 1.0)

  dz_g(1) = 0.0
  do k = 2,nk
    ! zr1 is positive and increases with depth, and dz_tgt is positive downward.
    dz_tgt = sgn*(z_new(k) - z_old(k))
    zr1 = sgn*(z_old(k) - z_old(1))

    !   First, handle the two simple and common cases that do not pass through
    ! the adjustment rate transition zone.
    if ((zr1 > zd) .and. (zr1 + wtd * dz_tgt > zd)) then
      dz_g(k) = sgn * wtd * dz_tgt
    elseif ((zr1 < zs) .and. (zr1 + dz_tgt < zs)) then
      dz_g(k) = sgn * dz_tgt
    else
      ! Find the new value by inverting the equation
      !   integral(0 to dz_new) Iwt(z) dz = dz_tgt
      ! This is trivial where Iwt is a constant, and agrees with the two limits above.

      ! Take test values at the transition points to figure out which segment
      ! the new value will be found in.
      if (zr1 >= zd) then
        Int_zd = Iwtd*(zd - zr1)
        Int_zs = Int_zd - dInt_zs_zd
      elseif (zr1 <= zs) then
        Int_zs = (zs - zr1)
        Int_zd = dInt_zs_zd + (zs - zr1)
      else
!        Int_zd = (zd - zr1) * (Iwtd + 0.5*(1.0 - Iwtd) * (zd - zr1) / (zd - zs))
        Int_zd = (zd - zr1) * (Iwtd*(0.5*(zd+zr1) - zs) + 0.5*(zd - zr1)) * Idzwt
        Int_zs = (zs - zr1) * (0.5*Iwtd * ((zr1 - zs)) + (zd - 0.5*(zr1+zs))) * Idzwt
        ! It has been verified that  Int_zs = Int_zd - dInt_zs_zd to within roundoff.
      endif

      if (dz_tgt >= Int_zd) then ! The new location is in the deep, slow region.
        dz_g(k) = sgn * ((zd-zr1) + wtd*(dz_tgt - Int_zd))
      elseif (dz_tgt <= Int_zs) then ! The new location is in the shallow region.
        dz_g(k) = sgn * ((zs-zr1) + (dz_tgt - Int_zs))
      else  ! We need to solve a quadratic equation for z_new.
        ! For accuracy, do the integral from the starting depth or the nearest
        ! edge of the transition region.  The results with each choice are
        ! mathematically equivalent, but differ in roundoff, and this choice
        ! should minimize the likelihood of inadvertently overlapping interfaces.
        if (zr1 <= zs) then ; dz0 = zs-zr1 ; z0 = zs ; F0 = dz_tgt - Int_zs
        elseif (zr1 >= zd) then ; dz0 = zd-zr1 ; z0 = zd ; F0 = dz_tgt - Int_zd
        else ; dz0 = 0.0 ; z0 = zr1 ; F0 = dz_tgt ; endif

        Bq = (dzwt + 2.0*Aq*(z0-zs))
        ! Solve the quadratic: Aq*(zn-z0)**2 + Bq*(zn-z0) - F0*dzwt = 0
        ! Note that b>=0, and the two terms in the standard form cancel for the right root.
        dz_g(k) = sgn * (dz0 + 2.0*F0*dzwt / (Bq + sqrt(Bq**2 + 4.0*Aq*F0*dzwt) ))

!       if (debug) then
!         dz0 = zs-zr1 ; z0 = zs ; F0 = dz_tgt - Int_zs ; Bq = (dzwt + 2.0*Aq*(z0-zs))
!         ddz_g_s(k) = sgn * (dz0 + 2.0*F0*dzwt / (Bq + sqrt(Bq**2 + 4.0*Aq*F0*dzwt) )) - dz_g(k)
!         dz0 = zd-zr1 ; z0 = zd ; F0 = dz_tgt - Int_zd ; Bq = (dzwt + 2.0*Aq*(z0-zs))
!         ddz_g_d(k) = sgn * (dz0 + 2.0*F0*dzwt / (Bq + sqrt(Bq**2 + 4.0*Aq*F0*dzwt) )) - dz_g(k)
!
!         if (abs(ddz_g_s(k)) > 1e-12*(abs(dz_g(k)) + abs(dz_g(k)+ddz_g_s(k)))) &
!           call MOM_error(WARNING, "filtered_grid_motion: Expect z_output to be tangled (sc).")
!         if (abs(ddz_g_d(k) - ddz_g_s(k)) > 1e-12*(abs(dz_g(k)+ddz_g_d(k)) + abs(dz_g(k)+ddz_g_s(k)))) &
!           call MOM_error(WARNING, "filtered_grid_motion: Expect z_output to be tangled.")
!       endif
      endif

    endif
  enddo
  dz_g(nk+1) = 0.0

  if (debug) then
    do k=1,nk+1 ; z_act(k) = z_old(k) + dz_g(k) ; enddo
    do k=2,nk+1
      if (sgn*((z_act(k))-z_act(k-1)) < -1e-15*(abs(z_act(k))+abs(z_act(k-1))) ) &
        call MOM_error(FATAL, "filtered_grid_motion: z_output is tangled.")
    enddo
  endif

end subroutine filtered_grid_motion

!> Builds a z*-ccordinate grid with partial steps (Adcroft and Campin, 2004).
!! z* is defined as
!!   z* = (z-eta)/(H+eta)*H  s.t. z*=0 when z=eta and z*=-H when z=-H .
subroutine build_zstar_grid( CS, G, GV, h, dzInterface, frac_shelf_h)

  ! Arguments
  type(regridding_CS),                          intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                        intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h  !< Layer thicknesses, in H
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)+1), intent(inout) :: dzInterface !< The change in interface depth in H.
  real, dimension(:,:),               optional, pointer       :: frac_shelf_h !< Fractional ice shelf coverage.
  ! Local variables
  integer :: i, j, k
  integer :: nz
  real    :: nominalDepth, totalThickness, dh
  real, dimension(SZK_(GV)+1) :: zOld, zNew
  real :: minThickness
  logical :: ice_shelf

  nz = GV%ke
  minThickness = CS%min_thickness
  ice_shelf = .false.
  if (present(frac_shelf_h)) then
    if (associated(frac_shelf_h)) ice_shelf = .true.
  endif

!$OMP parallel do default(none) shared(G,GV,dzInterface,CS,nz,h,frac_shelf_h, &
!$OMP                                  ice_shelf,minThickness) &
!$OMP                          private(nominalDepth,totalThickness, &
!$OMP                                  zNew,dh,zOld)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      if (G%mask2dT(i,j)==0.) then
        dzInterface(i,j,:) = 0.
        cycle
      endif

      ! Local depth (G%bathyT is positive)
      nominalDepth = G%bathyT(i,j)*GV%m_to_H

      ! Determine water column thickness
      totalThickness = 0.0
      do k = 1,nz
        totalThickness = totalThickness + h(i,j,k)
      end do

      zOld(nz+1) = - nominalDepth
      do k = nz,1,-1
        zOld(k) = zOld(k+1) + h(i,j,k)
      enddo

      if (ice_shelf) then
        if (frac_shelf_h(i,j) > 0.) then ! under ice shelf
           call build_zstar_column(CS, nz, nominalDepth, totalThickness, zNew, &
                                z_rigid_top = totalThickness-nominalDepth, &
                                eta_orig = zOld(1))
        else
           call build_zstar_column(CS, nz, nominalDepth, totalThickness, zNew)
        endif
      else
        call build_zstar_column(CS, nz, nominalDepth, totalThickness, zNew)
      endif

      ! Calculate the final change in grid position after blending new and old grids
      call filtered_grid_motion( CS, nz, zOld, zNew, dzInterface(i,j,:) )

#ifdef __DO_SAFETY_CHECKS__
      dh=max(nominalDepth,totalThickness)
      if (abs(zNew(1)-zOld(1))>(nz-1)*0.5*epsilon(dh)*dh) then
        write(0,*) 'min_thickness=',minThickness
        write(0,*) 'nominalDepth=',nominalDepth,'totalThickness=',totalThickness
        write(0,*) 'dzInterface(1) = ',dzInterface(i,j,1),epsilon(dh),nz
        do k=1,nz+1
          write(0,*) k,zOld(k),zNew(k)
        enddo
        do k=1,nz
          write(0,*) k,h(i,j,k),zNew(k)-zNew(k+1),CS%coordinateResolution(k)
        enddo
        call MOM_error( FATAL, &
               'MOM_regridding, build_zstar_grid(): top surface has moved!!!' )
      endif
#endif

      call adjust_interface_motion( nz, CS%min_thickness, h(i,j,:), dzInterface(i,j,:) )

    end do
  end do

end subroutine build_zstar_grid

!> Builds a z* coordinate with a minimum thickness
subroutine build_zstar_column(CS, nz, depth, total_thickness, zInterface, z_rigid_top, eta_orig)
  type(regridding_CS),   intent(in)    :: CS !< Regridding control structure
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive in m)
  real,                  intent(in)    :: total_thickness !< Column thickness (positive in m)
  real, dimension(nz+1), intent(inout) :: zInterface !< Absolute positions of interfaces
  real, optional,        intent(in)    :: z_rigid_top !< The height of a rigid top (negative in m)
  real, optional,        intent(in)    :: eta_orig !< The actual original height of the top (m)
  ! Local variables
  real :: eta, stretching, dh, min_thickness, z0_top, z_star
  integer :: k
  logical :: new_zstar_def

  new_zstar_def = .false.
  min_thickness = min( CS%min_thickness, total_thickness/real(nz) )
  z0_top = 0.
  if (present(z_rigid_top)) then
    z0_top = z_rigid_top
    new_zstar_def = .true.
  endif

  ! Position of free-surface (or the rigid top, for which eta ~ z0_top)
  eta = total_thickness - depth
  if (present(eta_orig)) eta = eta_orig

  ! Conventional z* coordinate:
  !   z* = (z-eta) / stretching   where stretching = (H+eta)/H
  !   z = eta + stretching * z*
  ! The above gives z*(z=eta) = 0, z*(z=-H) = -H.
  ! With a rigid top boundary at eta = z0_top then
  !   z* = z0 + (z-eta) / stretching   where stretching = (H+eta)/(H+z0)
  !   z = eta + stretching * (z*-z0) * stretching
  stretching = total_thickness / ( depth + z0_top )

  if (new_zstar_def) then
    ! z_star is the notional z* coordinate in absence of upper/lower topography
    z_star = 0. ! z*=0 at the free-surface
    zInterface(1) = eta ! The actual position of the top of the column
    do k = 2,nz
      z_star = z_star - CS%coordinateResolution(k-1)
      ! This ensures that z is below a rigid upper surface (ice shelf bottom)
      zInterface(k) = min( eta + stretching * ( z_star - z0_top ), z0_top )
      ! This ensures that the layer in inflated
      zInterface(k) = min( zInterface(k), zInterface(k-1) - min_thickness )
      ! This ensures that z is above or at the topography
      zInterface(k) = max( zInterface(k), -depth + real(nz+1-k) * min_thickness )
    enddo
    zInterface(nz+1) = -depth

  else
    ! Integrate down from the top for a notional new grid, ignoring topography
    ! The starting position is offset by z0_top which, if z0_top<0, will place
    ! interfaces above the rigid boundary.
    zInterface(1) = eta
    do k = 1,nz
      dh = stretching * CS%coordinateResolution(k) ! Notional grid spacing
      zInterface(k+1) = zInterface(k) - dh
    enddo

    ! Integrating up from the bottom adjusting interface position to accommodate
    ! inflating layers without disturbing the interface above
    zInterface(nz+1) = -depth
    do k = nz,1,-1
      if ( zInterface(k) < (zInterface(k+1) + min_thickness) ) then
        zInterface(k) = zInterface(k+1) + min_thickness
      endif
    enddo
  endif

end subroutine build_zstar_column


!------------------------------------------------------------------------------
! Build sigma grid
!------------------------------------------------------------------------------
subroutine build_sigma_grid( CS, G, GV, h, dzInterface )
!------------------------------------------------------------------------------
! This routine builds a grid based on terrain-following coordinates.
! The module parameter coordinateResolution(:) determines the resolution in
! sigma coordinate, dSigma(:). sigma-coordinates are defined by
!   sigma = (eta-z)/(H+eta)  s.t. sigma=0 at z=eta and sigma=1 at z=-H .
!------------------------------------------------------------------------------

  ! Arguments
  type(regridding_CS),                          intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                        intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h  !< Layer thicknesses, in H
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)+1), intent(inout) :: dzInterface !< The change in interface depth in H.

  ! Local variables
  integer :: i, j, k
  integer :: nz
  real    :: nominalDepth, totalThickness, dh
  real, dimension(SZK_(GV)+1) :: zOld, zNew

  nz = GV%ke

  do i = G%isc-1,G%iec+1
    do j = G%jsc-1,G%jec+1

      ! The rest of the model defines grids integrating up from the bottom
      nominalDepth = G%bathyT(i,j)*GV%m_to_H

      ! Determine water column height
      totalThickness = 0.0
      do k = 1,nz
        totalThickness = totalThickness + h(i,j,k)
      end do

      call build_sigma_column(CS, nz, nominalDepth, totalThickness, zNew)

      ! Calculate the final change in grid position after blending new and old grids
      zOld(nz+1) =  -nominalDepth
      do k = nz,1,-1
        zOld(k) = zOld(k+1) + h(i, j, k)
      end do

      call filtered_grid_motion( CS, nz, zOld, zNew, dzInterface(i,j,:) )

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
               'MOM_regridding, build_sigma_grid: top surface has moved!!!' )
      endif
      dzInterface(i,j,1) = 0.
      dzInterface(i,j,nz+1) = 0.
#endif

    end do
  end do

end subroutine build_sigma_grid

subroutine build_sigma_column(CS, nz, depth, totalThickness, zInterface)
  type(regridding_CS),   intent(in)    :: CS !< Regridding control structure
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive in m)
  real,                  intent(in)    :: totalThickness !< Column thickness (positive in m)
  real, dimension(nz+1), intent(inout) :: zInterface !< Absolute positions of interfaces

  ! Local variables
  integer :: k

  zInterface(nz+1) = -depth
  do k = nz,1,-1
    zInterface(k) = zInterface(k+1) + (totalThickness * CS%coordinateResolution(k))
    ! Adjust interface position to accomodate inflating layers
    ! without disturbing the interface above
    if (zInterface(k) < (zInterface(k+1) + CS%min_thickness)) then
      zInterface(k) = zInterface(k+1) + CS%min_thickness
    endif
  enddo

end subroutine build_sigma_column

!------------------------------------------------------------------------------
! Build grid based on target interface densities
!------------------------------------------------------------------------------
subroutine build_rho_grid( G, GV, h, tv, dzInterface, remapCS, CS )
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
  type(ocean_grid_type),                        intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h  !< Layer thicknesses, in H
  type(thermo_var_ptrs),                        intent(in)    :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)+1), intent(inout) :: dzInterface !< The change in interface depth in H
  type(remapping_CS),                           intent(in)    :: remapCS !< The remapping control structure
  type(regridding_CS),                          intent(in)    :: CS !< Regridding control structure

  ! Local variables
  integer :: nz
  integer :: i, j, k
  real    :: nominalDepth, totalThickness
  real, dimension(SZK_(GV)+1) :: zOld, zNew

  nz = GV%ke

  if (.not.CS%target_density_set) call MOM_error(FATAL, "build_rho_grid: "//&
        "Target densities must be set before build_rho_grid is called.")

  ! Build grid based on target interface densities
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      ! Local depth (G%bathyT is positive)
      nominalDepth = G%bathyT(i,j)*GV%m_to_H

      call build_rho_column(CS, remapCS, nz, nominalDepth, h(i, j, :)*GV%H_to_m, &
                            tv%T(i, j, :), tv%S(i, j, :), tv%eqn_of_state, zNew)

      if (CS%integrate_downward_for_e) then
        zOld(1) = 0.
        do k = 1,nz
          zOld(k+1) = zOld(k) - h(i,j,k)
        enddo
      else
        ! The rest of the model defines grids integrating up from the bottom
        zOld(nz+1) = - nominalDepth
        do k = nz,1,-1
          zOld(k) = zOld(k+1) + h(i,j,k)
        enddo
      endif

      ! Calculate the final change in grid position after blending new and old grids
      call filtered_grid_motion( CS, nz, zOld, zNew, dzInterface(i,j,:) )

#ifdef __DO_SAFETY_CHECKS__
      do k = 2,nz
        if (zNew(k) > zOld(1)) then
          write(0,*) 'zOld=',zOld
          write(0,*) 'zNew=',zNew
          call MOM_error( FATAL, 'MOM_regridding, build_rho_grid: '//&
               'interior interface above surface!' )
        endif
        if (zNew(k) > zNew(k-1)) then
          write(0,*) 'zOld=',zOld
          write(0,*) 'zNew=',zNew
          call MOM_error( FATAL, 'MOM_regridding, build_rho_grid: '//&
               'interior interfaces cross!' )
        endif
      enddo

      totalThickness = 0.0
      do k = 1,nz
        totalThickness = totalThickness + h(i,j,k)
      enddo

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
               'MOM_regridding, build_rho_grid: top surface has moved!!!' )
      endif
#endif

    end do  ! end loop on i
  end do  ! end loop on j

end subroutine build_rho_grid


subroutine build_rho_column(CS, remapCS, nz, depth, h, T, S, eqn_of_state, zInterface)
! The algorithn operates as follows within each
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

  type(regridding_CS),   intent(in)    :: CS !< Regridding control structure
  type(remapping_CS),    intent(in)    :: remapCS !< Remapping parameters and options
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive in m)
  real, dimension(nz),   intent(in)    :: h  !< Layer thicknesses, in m
  real, dimension(nz),   intent(in)    :: T, S !< T and S for column
  type(EOS_type),        pointer       :: eqn_of_state !< Equation of state structure
  real, dimension(nz+1), intent(inout) :: zInterface !< Absolute positions of interfaces

  ! Local variables
  integer   :: k, m
  integer   :: map_index
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
  integer   :: ppoly_degree         ! The actual degree of the polynomials.
  real, dimension(nz) :: p, densities, T_tmp, S_tmp, Tmp
  integer, dimension(nz) :: mapping
  real :: dh
  real, dimension(nz) :: h0, h1, hTmp
  real, dimension(nz+1) :: x0, x1, xTmp

  threshold = CS%min_thickness
  p(:) = CS%ref_pressure
  T_tmp(:) = T(:)
  S_tmp(:) = S(:)
  h0(:) = h(:)

  ! Start iterations to build grid
  m = 1
  deviation = 1e10
  do while ( ( m <= NB_REGRIDDING_ITERATIONS ) .and. &
             ( deviation > DEVIATION_TOLERANCE ) )

    ! Count number of nonzero layers within current water column
    count_nonzero_layers = 0
    do k = 1,nz
      if ( h0(k) > threshold ) then
        count_nonzero_layers = count_nonzero_layers + 1
      end if
    end do

    ! If there is at most one nonzero layer, stop here (no regridding)
    if ( count_nonzero_layers <= 1 ) then
      h1(:) = h0(:)
      exit  ! stop iterations here
    end if

    ! Build new grid containing only nonzero layers
    map_index = 1
    correction = 0.0
    do k = 1,nz
      if ( h0(k) > threshold ) then
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
      if ( hTmp(k) > max_thickness ) then
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
    call calculate_density( T_tmp, S_tmp, p, densities,&
                             1, nz, eqn_of_state )

    do k = 1,count_nonzero_layers
      densities(k) = densities(mapping(k))
    end do

    ! One regridding iteration
    call regridding_set_ppolys(densities, CS, count_nonzero_layers, hTmp, &
                               ppoly_i_E, ppoly_i_S, ppoly_i_coefficients, ppoly_degree)
    ! Based on global density profile, interpolate to generate a new grid
    call interpolate_grid(count_nonzero_layers, hTmp, xTmp, ppoly_i_E, ppoly_i_coefficients, &
                          CS%target_density, ppoly_degree, nz, h1, x1 )

    call old_inflate_layers_1d( CS%min_thickness, nz, h1 )
    x1(1) = 0.0 ; do k = 1,nz ; x1(k+1) = x1(k) + h1(k) ; end do

    ! Remap T and S from previous grid to new grid
    do k = 1,nz
      h1(k) = x1(k+1) - x1(k)
    end do

    call remapping_core_h(nz, h0, S, nz, h1, Tmp, remapCS)
    S_tmp(:) = Tmp(:)

    call remapping_core_h(nz, h0, T, nz, h1, Tmp, remapCS)
    T_tmp(:) = Tmp(:)

    ! Compute the deviation between two successive grids
    deviation = 0.0
    x0(1) = 0.0
    x1(1) = 0.0
    do k = 2,nz
      x0(k) = x0(k-1) + h0(k-1)
      x1(k) = x1(k-1) + h1(k-1)
      deviation = deviation + (x0(k)-x1(k))**2
    end do
    deviation = sqrt( deviation / (nz-1) )

    m = m + 1

    ! Copy final grid onto start grid for next iteration
    h0(:) = h1(:)

  end do ! end regridding iterations

  if (CS%integrate_downward_for_e) then
    zInterface(1) = 0.
    do k = 1,nz
      zInterface(k+1) = zInterface(k) - h1(k)
      ! Adjust interface position to accomodate inflating layers
      ! without disturbing the interface above
    enddo
  else
    ! The rest of the model defines grids integrating up from the bottom
    zInterface(nz+1) = -depth
    do k = nz,1,-1
      zInterface(k) = zInterface(k+1) + h1(k)
      ! Adjust interface position to accomodate inflating layers
      ! without disturbing the interface above
    enddo
  endif

end subroutine build_rho_column


!> Builds a simple HyCOM-like grid with the deepest location of potential
!! density interpolated from the column profile and a clipping of depth for
!! each interface to a fixed z* or p* grid.  This should probably be (optionally?)
!! changed to find the nearest location of the target density.
!! \remark { Based on Bleck, 2002: An oceanice general circulation model framed in
!! hybrid isopycnic-Cartesian coordinates, Ocean Modelling 37, 55-88.
!! http://dx.doi.org/10.1016/S1463-5003(01)00012-9 }
subroutine build_grid_HyCOM1( G, GV, h, tv, dzInterface, remapCS, CS )
  type(ocean_grid_type),                       intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),                     intent(in)    :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h  !< Existing model thickness, in H units
  type(thermo_var_ptrs),                       intent(in)    :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: dzInterface !< Changes in interface position
  type(remapping_CS),                          intent(in)    :: remapCS !< Remapping control structure
  type(regridding_CS),                         intent(in)    :: CS !< Regridding control structure

  ! Local variables
  real, dimension(SZK_(GV)+1) :: z_col, z_col_new ! Interface positions relative to the surface in H units (m or kg m-2)
  real, dimension(SZK_(GV)+1) :: dz_col  ! The realized change in z_col in H units (m or kg m-2)
  real, dimension(SZK_(GV))   :: p_col   ! Layer pressure in Pa
  integer   :: i, j, k, nz
  real :: depth

  nz = GV%ke

  if (.not.CS%target_density_set) call MOM_error(FATAL, "build_grid_HyCOM1 : "//&
        "Target densities must be set before build_grid_HyCOM1 is called.")

  ! Build grid based on target interface densities
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
    if (G%mask2dT(i,j)>0.) then

      depth = G%bathyT(i,j) * GV%m_to_H

      z_col(1) = 0. ! Work downward rather than bottom up
      do K = 1, nz
        z_col(K+1) = z_col(K) + h(i,j,k) ! Work in units of h (m or Pa)
        p_col(k) = CS%ref_pressure + CS%compressibility_fraction * &
             ( 0.5 * ( z_col(K) + z_col(K+1) ) * GV%H_to_Pa - CS%ref_pressure )
      enddo

      call build_hycom1_column(CS, remapCS, tv%eqn_of_state, nz, depth, &
                               h(i, j, :), tv%T(i, j, :), tv%S(i, j, :), p_col, z_col, z_col_new)

      ! Calculate the final change in grid position after blending new and old grids
      call filtered_grid_motion( CS, nz, z_col, z_col_new, dz_col )
      do K=1,nz+1 ; dzInterface(i,j,K) = -dz_col(K) ; enddo

      ! This adjusts things robust to round-off errors
      call adjust_interface_motion( nz, CS%min_thickness, h(i,j,:), dzInterface(i,j,:) )

    else ! on land
      dzInterface(i,j,:) = 0.
    endif ! mask2dT
  enddo; enddo ! i,j

end subroutine build_grid_HyCOM1

subroutine build_hycom1_column(CS, remapCS, eqn_of_state, nz, depth, h, T, S, p_col, z_col, z_col_new)
  type(regridding_CS),   intent(in)    :: CS !< Regridding control structure
  type(remapping_CS),    intent(in)    :: remapCS !< Remapping parameters and options
  type(EOS_type),        pointer       :: eqn_of_state !< Equation of state structure
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive in H)
  real, dimension(nz),   intent(in)    :: T, S !< T and S for column
  real, dimension(nz),   intent(in)    :: h  !< Layer thicknesses, in m
  real, dimension(nz),   intent(in)    :: p_col !< Layer pressure in Pa
  real, dimension(nz+1), intent(in)    :: z_col ! Interface positions relative to the surface in H units (m or kg m-2)
  real, dimension(nz+1), intent(inout) :: z_col_new !< Absolute positions of interfaces

  ! Local variables
  integer   :: k
  real, dimension(nz) :: rho_col, h_col_new ! Layer quantities
  real, dimension(CS%nk,2) :: ppoly_i_E ! Edge value of polynomial
  real, dimension(CS%nk,2) :: ppoly_i_S ! Edge slope of polynomial
  real, dimension(CS%nk,CS%degree_i+1) :: ppoly_i_coefficients ! Coefficients of polynomial
  real :: stretching ! z* stretching, converts z* to z.
  real :: nominal_z ! Nominal depth of interface is using z* (m or Pa)
  real :: hNew
  logical :: maximum_depths_set ! If true, the maximum depths of interface have been set.
  logical :: maximum_h_set      ! If true, the maximum layer thicknesses have been set.
  integer :: ppoly_degree

  maximum_depths_set = allocated(CS%max_interface_depths)
  maximum_h_set = allocated(CS%max_layer_thickness)

  ! Work bottom recording potential density
  call calculate_density(T, S, p_col, rho_col, 1, nz, eqn_of_state)
  ! This ensures the potential density profile is monotonic
  ! although not necessarily single valued.
  do k = nz-1, 1, -1
    rho_col(k) = min( rho_col(k), rho_col(k+1) )
  enddo

  ! Interpolates for the target interface position with the rho_col profile
  call regridding_set_ppolys(rho_col, CS, nz, h(:), ppoly_i_E, ppoly_i_S, &
                             ppoly_i_coefficients, ppoly_degree)
  ! Based on global density profile, interpolate to generate a new grid
  call interpolate_grid(nz, h(:), z_col, ppoly_i_E, ppoly_i_coefficients, &
                        CS%target_density, ppoly_degree, nz, h_col_new, z_col_new)

  ! Sweep down the interfaces and make sure that the interface is at least
  ! as deep as a nominal target z* grid
  nominal_z = 0.
  stretching = z_col(nz+1) / depth ! Stretches z* to z
  do k = 2, nz+1
    nominal_z = nominal_z + CS%coordinateResolution(k-1) * stretching
    z_col_new(k) = max( z_col_new(k), nominal_z )
    z_col_new(k) = min( z_col_new(k), z_col(nz+1) )
  enddo

  if (maximum_depths_set .and. maximum_h_set) then ; do k=2,nz
    ! The loop bounds are 2 & nz so the top and bottom interfaces do not move.
    ! Recall that z_col_new is positive downward.
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K), &
                       z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; elseif (maximum_depths_set) then ; do K=2,nz
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K))
  enddo ; elseif (maximum_h_set) then ; do k=2,nz
    z_col_new(K) = min(z_col_new(K), z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; endif

end subroutine build_hycom1_column


!> Builds a grid that tracks density interfaces for water that is denser than
!! the surface density plus an increment of some number of layers, and uses all
!! lighter layers uniformly above this location.  Note that this amounts to
!! interpolating to find the depth of an arbitrary (non-integer) interface index
!! which should make the results vary smoothly in space to the extent that the
!! surface density and interior stratification vary smoothly in space.  Over
!! shallow topography, this will tend to give a uniform sigma-like coordinate.
!! For sufficiently shallow water, a minimum grid spacing is used to avoid
!! certain instabilities.
subroutine build_grid_SLight( G, GV, h, tv, dzInterface, remapCS, CS )
  type(ocean_grid_type),                       intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),                     intent(in)    :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h  !< Existing model thickness, in H units
  type(thermo_var_ptrs),                       intent(in)    :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: dzInterface !< Changes in interface position
  type(remapping_CS),                          intent(in)    :: remapCS !< Remapping control structure
  type(regridding_CS),                         intent(in)    :: CS !< Regridding control structure

  real, dimension(SZK_(GV)+1) :: z_col, z_col_new ! Interface positions relative to the surface in H units (m or kg m-2)
  real, dimension(SZK_(GV)+1) :: dz_col  ! The realized change in z_col in H units (m or kg m-2)
  real, dimension(SZK_(GV))   :: p_col   ! Layer pressure in Pa
  real :: depth
  integer :: i, j, k, nz

  nz = GV%ke

  if (.not.CS%target_density_set) call MOM_error(FATAL, "build_grid_SLight : "//&
        "Target densities must be set before build_grid_SLight is called.")

  ! Build grid based on target interface densities
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
    if (G%mask2dT(i,j)>0.) then

      depth = G%bathyT(i,j) * GV%m_to_H
      z_col(1) = 0. ! Work downward rather than bottom up
      do K=1,nz
        z_col(K+1) = z_col(K) + h(i, j, k) ! Work in units of h (m or Pa)
        p_col(k) = CS%ref_pressure + CS%compressibility_fraction * &
                    ( 0.5 * ( z_col(K) + z_col(K+1) ) * GV%H_to_Pa - CS%ref_pressure )
      enddo

      call build_slight_column(CS, remapCS, tv%eqn_of_state, GV%H_to_Pa, GV%m_to_H, &
                          GV%H_subroundoff, nz, depth, &
                          h(i, j, :), tv%T(i, j, :), tv%S(i, j, :), p_col, z_col, z_col_new)

      ! Calculate the final change in grid position after blending new and old grids
      call filtered_grid_motion( CS, nz, z_col, z_col_new, dz_col )
      do K=1,nz+1 ; dzInterface(i,j,K) = -dz_col(K) ; enddo
#ifdef __DO_SAFETY_CHECKS__
      if (dzInterface(i,j,1) /= 0.) stop 'build_grid_SLight: Surface moved?!'
      if (dzInterface(i,j,nz+1) /= 0.) stop 'build_grid_SLight: Bottom moved?!'
#endif

      ! This adjusts things robust to round-off errors
      call adjust_interface_motion( nz, CS%min_thickness, h(i,j,:), dzInterface(i,j,:) )

    else ! on land
      dzInterface(i,j,:) = 0.
    endif ! mask2dT
  enddo; enddo ! i,j

end subroutine build_grid_SLight

subroutine build_slight_column(CS, remapCS, eqn_of_state, H_to_Pa, m_to_H, H_subroundoff, &
                                nz, depth, h_col, T_col, S_col, p_col, z_col, z_col_new)
  type(regridding_CS),   intent(in)    :: CS !< Regridding control structure
  type(remapping_CS),    intent(in)    :: remapCS !< Remapping parameters and options
  type(EOS_type),        pointer       :: eqn_of_state !< Equation of state structure
  real,                  intent(in)    :: H_to_Pa !< GV%H_to_Pa
  real,                  intent(in)    :: m_to_H  !< GV%m_to_H
  real,                  intent(in)    :: H_subroundoff !< GV%H_subroundoff
  integer,               intent(in)    :: nz !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive in m)
  real, dimension(nz),   intent(in)    :: T_col, S_col !< T and S for column
  real, dimension(nz),   intent(in)    :: h_col !< Layer thicknesses, in m
  real, dimension(nz),   intent(in)    :: p_col !< Layer quantities
  real, dimension(nz+1), intent(in)    :: z_col !< Interface positions relative to the surface in H units (m or kg m-2)
  real, dimension(nz+1), intent(inout) :: z_col_new !< Absolute positions of interfaces

  ! Local variables
  real, dimension(nz) :: rho_col ! Layer quantities
  real, dimension(nz) :: T_f, S_f  ! Filtered ayer quantities
  logical, dimension(nz+1) :: reliable  ! If true, this interface is in a reliable position.
  real, dimension(nz+1) :: T_int, S_int ! Temperature and salinity interpolated to interfaces.
  real, dimension(nz+1) :: rho_tmp, drho_dp, p_IS, p_R
  real, dimension(nz+1) :: drhoIS_dT, drhoIS_dS
  real, dimension(nz+1) :: drhoR_dT, drhoR_dS
  real, dimension(nz+1) :: strat_rat
  real :: H_to_cPa
  real :: drIS, drR, Fn_now, I_HStol, Fn_zero_val
  real :: z_int_unst
  real :: dz      ! A uniform layer thickness in very shallow water, in H.
  real :: dz_ur   ! The total thickness of an unstable region, in H.
  real :: wgt, cowgt  ! A weight and its complement, nondim.
  real :: rho_ml_av ! The average potential density in a near-surface region, in kg m-3.
  real :: H_ml_av ! A thickness to try to use in taking the near-surface average, in H.
  real :: rho_x_z ! A cumulative integral of a density, in kg m-3 H.
  real :: z_wt    ! The thickness actually used in taking the near-surface average, in H.
  real :: k_interior  ! The (real) value of k where the interior grid starts.
  real :: k_int2      ! The (real) value of k where the interior grid starts.
  real :: z_interior  ! The depth where the interior grid starts, in H.
  real :: z_ml_fix    ! The depth at which the fixed-thickness near-surface layers end, in H.
  real :: dz_dk       ! The thickness of layers between the fixed-thickness
                      ! near-surface layars and the interior, in H.
  real :: Lfilt       ! A filtering lengthscale, in H.
  logical :: maximum_depths_set ! If true, the maximum depths of interface have been set.
  logical :: maximum_h_set      ! If true, the maximum layer thicknesses have been set.
  real :: k2_used, k2here, dz_sum, z_max
  integer :: k2
  real :: h_tr, b_denom_1, b1, d1 ! Temporary variables used by the tridiagonal solver.
  real, dimension(nz) :: c1  ! Temporary variables used by the tridiagonal solver.
  integer :: kur1, kur2  ! The indicies at the top and bottom of an unreliable region.
  integer :: kur_ss      ! The index to start with in the search for the next unstable region.
  integer :: i, j, k, nkml

  maximum_depths_set = allocated(CS%max_interface_depths)
  maximum_h_set = allocated(CS%max_layer_thickness)

  if (z_col(nz+1) - z_col(1) < nz*CS%min_thickness) then
    ! This is a nearly massless total depth, so distribute the water evenly.
    dz = (z_col(nz+1) - z_col(1)) / real(nz)
    do K=2,nz ; z_col_new(K) = z_col(1) + dz*real(K-1) ; enddo
  else
    call calculate_density(T_col, S_col, p_col, rho_col, 1, nz, &
                           eqn_of_state)

    ! Find the locations of the target potential densities, flagging
    ! locations in apparently unstable regions as not reliable.
    call rho_interfaces_col(rho_col, h_col, z_col, CS%target_density, nz, &
                            z_col_new, CS, reliable, debug=.true.)

    ! Ensure that the interfaces are at least CS%min_thickness apart.
    if (CS%min_thickness > 0.0) then
      ! Move down interfaces below overly thin layers.
      do K=2,nz ; if (z_col_new(K) < z_col_new(K-1) + CS%min_thickness) then
        z_col_new(K) = z_col_new(K-1) + CS%min_thickness
      endif ; enddo
      ! Now move up any interfaces that are too close to the bottom.
      do K=nz,2,-1 ; if (z_col_new(K) > z_col_new(K+1) - CS%min_thickness) then
        z_col_new(K) = z_col_new(K+1) - CS%min_thickness
      else
        exit ! No more interfaces can be too close to the bottom.
      endif ; enddo
    endif

    ! Fix up the unreliable regions.
    kur_ss = 2 ! reliable(1) and reliable(nz+1) must always be true.
    do
      ! Search for the uppermost unreliable interface postion.
      kur1 = nz+2
      do K=kur_ss,nz ; if (.not.reliable(K)) then
        kur1 = K ; exit
      endif ; enddo
      if (kur1 > nz) exit ! Everything is now reliable.

      kur2 = kur1-1 ! For error checking.
      do K=kur1+1,nz+1 ; if (reliable(K)) then
        kur2 = K-1 ; kur_ss = K ; exit
      endif ; enddo
      if (kur2 < kur1) call MOM_error(FATAL, "Bad unreliable range.")

      dz_ur = z_col_new(kur2+1) - z_col_new(kur1-1)
  !        drho = CS%target_density(kur2+1) - CS%target_density(kur1-1)
      ! Perhaps reset the wgt and cowgt depending on how bad the old interface
      ! locations were.
      wgt = 1.0 ; cowgt = 0.0 ! = 1.0-wgt
      do K=kur1,kur2
        z_col_new(K) = cowgt*z_col_new(K) + &
              wgt * (z_col_new(kur1-1) + dz_ur*(K - (kur1-1)) / ((kur2 - kur1) + 2))
      enddo
    enddo

    ! Determine which interfaces are in the s-space region and the depth extent
    ! of this region.
    z_wt = 0.0 ; rho_x_z = 0.0
    H_ml_av = m_to_H*CS%Rho_ml_avg_depth
    do k=1,nz
      if (z_wt + h_col(k) >= H_ml_av) then
        rho_x_z = rho_x_z + rho_col(k) * (H_ml_av - z_wt)
        z_wt = H_ml_av
        exit
      else
        rho_x_z =  rho_x_z + rho_col(k) * h_col(k)
        z_wt = z_wt + h_col(k)
      endif
    enddo
    if (z_wt > 0.0) rho_ml_av = rho_x_z / z_wt

    nkml = CS%nz_fixed_surface
    ! Find the interface that matches rho_ml_av.
    if (rho_ml_av <= CS%target_density(nkml)) then
      k_interior = CS%nlay_ml_offset + real(nkml)
    elseif (rho_ml_av > CS%target_density(nz+1)) then
      k_interior = real(nz+1)
    else ; do K=nkml,nz
      if ((rho_ml_av >= CS%target_density(K)) .and. &
          (rho_ml_av <  CS%target_density(K+1))) then
        k_interior = (CS%nlay_ml_offset + K) + &
                (rho_ml_av - CS%target_density(K)) / &
                (CS%target_density(K+1) - CS%target_density(K))
        exit
      endif
    enddo ; endif
    if (k_interior > real(nz+1)) k_interior = real(nz+1)

    ! Linearly interpolate to find z_interior.  This could be made more sophisticated.
    K = int(ceiling(k_interior))
    z_interior = (K-k_interior)*z_col_new(K-1) + (1.0+(k_interior-K))*z_col_new(K)

    if (CS%fix_haloclines) then
  !       ! Identify regions above the reference pressure where the chosen
  !       ! potential density significantly underestimates the actual
  !       ! stratification, and use these to find a second estimate of
  !       ! z_int_unst and k_interior.

      if (CS%halocline_filter_length > 0.0) then
        Lfilt = CS%halocline_filter_length*m_to_H

        ! Filter the temperature and salnity with a fixed lengthscale.
        h_tr = h_col(1) + H_subroundoff
        b1 = 1.0 / (h_tr + Lfilt) ; d1 = h_tr * b1
        T_f(1) = (b1*h_tr)*T_col(1) ;  S_f(1) = (b1*h_tr)*S_col(1)
        do k=2,nz
          c1(k) = Lfilt * b1
          h_tr = h_col(k) + H_subroundoff ; b_denom_1 = h_tr + d1*Lfilt
          b1 = 1.0 / (b_denom_1 + Lfilt) ; d1 = b_denom_1 * b1
          T_f(k) = b1 * (h_tr*T_col(k) + Lfilt*T_f(k-1))
          S_f(k) = b1 * (h_tr*S_col(k) + Lfilt*S_f(k-1))
        enddo
        do k=nz-1,1,-1
          T_f(k) = T_f(k) + c1(k+1)*T_f(k+1) ; S_f(k) = S_f(k) + c1(k+1)*S_f(k+1)
        enddo
      else
        do k=1,nz ; T_f(k) = T_col(k) ; S_f(k) = S_col(k) ; enddo
      endif

      T_int(1) = T_f(1) ; S_int(1) = S_f(1)
      do K=2,nz
        T_int(K) = 0.5*(T_f(k-1) + T_f(k)) ; S_int(K) = 0.5*(S_f(k-1) + S_f(k))
        p_IS(K) = z_col(K) * H_to_Pa
        p_R(K) = CS%ref_pressure + CS%compressibility_fraction * ( p_IS(K) - CS%ref_pressure )
      enddo
      T_int(nz+1) = T_f(nz) ; S_int(nz+1) = S_f(nz)
      p_IS(nz+1) = z_col(nz+1) * H_to_Pa
      call calculate_density_derivs(T_int, S_int, p_IS, drhoIS_dT, drhoIS_dS, 2, nz-1, &
                                    eqn_of_state)
      call calculate_density_derivs(T_int, S_int, p_R, drhoR_dT, drhoR_dS, 2, nz-1, &
                                    eqn_of_state)
      if (CS%compressibility_fraction > 0.0) then
        call calculate_compress(T_int, S_int, p_R, rho_tmp, drho_dp, 2, nz-1, &
                                      eqn_of_state)
      else
        do K=2,nz ; drho_dp(K) = 0.0 ; enddo
      endif

      H_to_cPa = CS%compressibility_fraction*H_to_Pa
      strat_rat(1) = 1.0
      do K=2,nz
        drIS = drhoIS_dT(K) * (T_f(k) - T_f(k-1)) + &
               drhoIS_dS(K) * (S_f(k) - S_f(k-1))
        drR = (drhoR_dT(K) * (T_f(k) - T_f(k-1)) + &
               drhoR_dS(K) * (S_f(k) - S_f(k-1))) + &
              drho_dp(K) * (H_to_cPa*0.5*(h_col(k) + h_col(k-1)))

        if (drIS <= 0.0) then
          strat_rat(K) = 2.0 ! Maybe do this? => ; if (drR < 0.0) strat_rat(K) = -2.0
        else
          strat_rat(K) = 2.0*max(drR,0.0) / (drIS + abs(drR))
        endif
      enddo
      strat_rat(nz+1) = 1.0

      z_int_unst = 0.0 ; Fn_now = 0.0
      Fn_zero_val = min(2.0*CS%halocline_strat_tol, &
                        0.5*(1.0 + CS%halocline_strat_tol))
      if (CS%halocline_strat_tol > 0.0) then
        ! Use Adcroft's reciprocal rule.
        I_HStol = 0.0 ; if (Fn_zero_val - CS%halocline_strat_tol > 0.0) &
          I_HStol = 1.0 / (Fn_zero_val - CS%halocline_strat_tol)
        do k=nz,1,-1 ; if (CS%ref_pressure > p_IS(k+1)) then
          z_int_unst = z_int_unst + Fn_now * h_col(k)
          if (strat_rat(K) <= Fn_zero_val) then
            if (strat_rat(K) <= CS%halocline_strat_tol) then ; Fn_now = 1.0
            else
              Fn_now = max(Fn_now, (Fn_zero_val - strat_rat(K)) * I_HStol)
            endif
          endif
        endif ; enddo
      else
        do k=nz,1,-1 ; if (CS%ref_pressure > p_IS(k+1)) then
          z_int_unst = z_int_unst + Fn_now * h_col(k)
          if (strat_rat(K) <= CS%halocline_strat_tol) Fn_now = 1.0
        endif ; enddo
      endif

      if (z_interior < z_int_unst) then
        ! Find a second estimate of the extent of the s-coordinate region.
        kur1 = max(int(ceiling(k_interior)),2)
        if (z_col_new(kur1-1) < z_interior) then
          k_int2 = kur1
          do K = kur1,nz+1 ; if (z_col_new(K) >= z_int_unst) then
            ! This is linear interpolation again.
            if (z_col_new(K-1) >= z_int_unst) &
              call MOM_error(FATAL,"build_grid_SLight, bad halocline structure.")
            k_int2 = real(K-1) + (z_int_unst - z_col_new(K-1)) / &
                                     (z_col_new(K) - z_col_new(K-1))
            exit
          endif ; enddo
          if (z_col_new(nz+1) < z_int_unst) then
            ! This should be unnecessary.
            z_int_unst = z_col_new(nz+1) ; k_int2 = real(nz+1)
          endif

          ! Now take the larger values.
          if (k_int2 > k_interior) then
            k_interior = k_int2 ; z_interior = z_int_unst
          endif
        endif
      endif
    endif  ! fix_haloclines

    z_col_new(1) = 0.0
    do K=2,nkml+1
      z_col_new(K) = min((K-1)*CS%dz_ml_min, &
                         z_col_new(nz+1) - CS%min_thickness*(nz+1-K))
    enddo
    z_ml_fix = z_col_new(nkml+1)
    if (z_interior > z_ml_fix) then
      dz_dk = (z_interior - z_ml_fix) / (k_interior - (nkml+1))
      do K=nkml+2,int(floor(k_interior))
        z_col_new(K) = z_ml_fix + dz_dk * (K - (nkml+1))
      enddo
    else ! The fixed-thickness z-region penetrates into the interior.
      do K=nkml+2,nz
        if (z_col_new(K) <= z_col_new(CS%nz_fixed_surface+1)) then
          z_col_new(K) = z_col_new(CS%nz_fixed_surface+1)
        else ; exit ; endif
      enddo
    endif

    if (maximum_depths_set .and. maximum_h_set) then ; do k=2,nz
      ! The loop bounds are 2 & nz so the top and bottom interfaces do not move.
      ! Recall that z_col_new is positive downward.
      z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K), &
                            z_col_new(K-1) + CS%max_layer_thickness(k-1))
    enddo ; elseif (maximum_depths_set) then ; do K=2,nz
      z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K))
    enddo ; elseif (maximum_h_set) then ; do k=2,nz
      z_col_new(K) = min(z_col_new(K), z_col_new(K-1) + CS%max_layer_thickness(k-1))
    enddo ; endif

  endif ! Total thickness exceeds nz*CS%min_thickness.

end subroutine build_slight_column

!> Finds the new interface locations in a column of water that match the
!! prescribed target densities.
subroutine rho_interfaces_col(rho_col, h_col, z_col, rho_tgt, nz, z_col_new, &
                              CS, reliable, debug)
  integer,               intent(in)    :: nz      !< Number of layers
  real, dimension(nz),   intent(in)    :: rho_col !< Initial layer reference densities.
  real, dimension(nz),   intent(in)    :: h_col   !< Initial layer thicknesses.
  real, dimension(nz+1), intent(in)    :: z_col   !< Initial interface heights.
  real, dimension(nz+1), intent(in)    :: rho_tgt !< Interface target densities.
  real, dimension(nz+1), intent(inout) :: z_col_new !< New interface heights.
  type(regridding_CS),   intent(in)    :: CS      !< Regridding control structure
  logical, dimension(nz+1), intent(inout) :: reliable !< If true, the interface positions
                                                  !! are well defined from a stable region.
  logical, optional,     intent(in) :: debug      !< If present and true, do debugging checks.

  real, dimension(nz+1) :: ru_max_int ! The maximum and minimum densities in
  real, dimension(nz+1) :: ru_min_int ! an unstable region around an interface.
  real, dimension(nz)   :: ru_max_lay ! The maximum and minimum densities in
  real, dimension(nz)   :: ru_min_lay ! an unstable region containing a layer.
  real, dimension(nz,2) :: ppoly_i_E ! Edge value of polynomial
  real, dimension(nz,2) :: ppoly_i_S ! Edge slope of polynomial
  real, dimension(nz,CS%degree_i+1) :: ppoly_i_coefficients ! Coefficients of polynomial
  logical, dimension(nz)   :: unstable_lay ! If true, this layer is in an unstable region.
  logical, dimension(nz+1) :: unstable_int ! If true, this interface is in an unstable region.
  real :: rt  ! The current target density, in kg m-3.
  real :: zf  ! The fractional z-position within a layer of the target density.
  real :: rfn
  real :: a(5) ! Coefficients of a local polynomial minus the target density.
  real :: zf1, zf2, rfn1, rfn2
  real :: drfn_dzf, sgn, delta_zf, zf_prev
  real :: tol
  logical :: k_found ! If true, the position has been found.
  integer :: k_layer ! The index of the stable layer containing an interface.
  integer :: ppoly_degree
  integer :: k, k1, k1_min, itt, max_itt, m

  real :: z_sgn  ! 1 or -1, depending on whether z increases with increasing K.
  logical :: debugging

  debugging = .false. ; if (present(debug)) debugging = debug
  max_itt = NR_ITERATIONS
  tol = NR_TOLERANCE

  z_sgn = 1.0 ; if ( z_col(1) > z_col(nz+1) ) z_sgn = -1.0
  if (debugging) then
    do K=1,nz
      if (abs((z_col(K+1) - z_col(K)) - z_sgn*h_col(k)) > &
          1.0e-14*(abs(z_col(K+1)) + abs(z_col(K)) + abs(h_col(k))) ) &
        call MOM_error(FATAL, "rho_interfaces_col: Inconsistent z_col and h_col")
    enddo
  endif

  if ( z_col(1) == z_col(nz+1) ) then
    ! This is a massless column!
    do K=1,nz+1 ; z_col_new(K) = z_col(1) ; reliable(K) = .true. ; enddo
    return
  endif

  ! This sets up the piecewise polynomials based on the rho_col profile.
  call regridding_set_ppolys(rho_col, CS, nz, h_col, ppoly_i_E, ppoly_i_S, &
                             ppoly_i_coefficients, ppoly_degree)

  ! Determine the density ranges of unstably stratified segments.
  ! Interfaces that start out in an unstably stratified segment can
  ! only escape if they are outside of the bounds of that segment, and no
  ! interfaces are ever mapped into an unstable segment.
  unstable_int(1) = .false.
  ru_max_int(1) = ppoly_i_E(1,1)

  unstable_lay(1) = (ppoly_i_E(1,1) > ppoly_i_E(1,2))
  ru_max_lay(1) = max(ppoly_i_E(1,1), ppoly_i_E(1,2))

  do K=2,nz
    unstable_int(K) = (ppoly_i_E(k-1,2) > ppoly_i_E(k,1))
    ru_max_int(K) = max(ppoly_i_E(k-1,2), ppoly_i_E(k,1))
    ru_min_int(K) = min(ppoly_i_E(k-1,2), ppoly_i_E(k,1))
    if (unstable_int(K) .and. unstable_lay(k-1)) &
      ru_max_int(K) = max(ru_max_lay(k-1), ru_max_int(K))

    unstable_lay(k) = (ppoly_i_E(k,1) > ppoly_i_E(k,2))
    ru_max_lay(k) = max(ppoly_i_E(k,1), ppoly_i_E(k,2))
    ru_min_lay(k) = min(ppoly_i_E(k,1), ppoly_i_E(k,2))
    if (unstable_lay(k) .and. unstable_int(K)) &
      ru_max_lay(k) = max(ru_max_int(K), ru_max_lay(k))
  enddo
  unstable_int(nz+1) = .false.
  ru_min_int(nz+1) = ppoly_i_E(nz,2)

  do K=nz,1,-1
    if (unstable_lay(k) .and. unstable_int(K+1)) &
      ru_min_lay(k) = min(ru_min_int(K+1), ru_min_lay(k))

    if (unstable_int(K) .and. unstable_lay(k)) &
      ru_min_int(K) = min(ru_min_lay(k), ru_min_int(K))
  enddo

  z_col_new(1) = z_col(1) ; reliable(1) = .true.
  k1_min = 1
  do K=2,nz ! Find the locations of the various target densities for the interfaces.
    rt = rho_tgt(K)
    k_layer = -1
    k_found = .false.

    ! Many light layers are found at the top, so start there.
    if (rt <= ppoly_i_E(k1_min,1)) then
      z_col_new(K) = z_col(k1_min)
      k_found = .true.
      ! Do not change k1_min for the next layer.
    elseif (k1_min == nz+1) then
      z_col_new(K) = z_col(nz+1)
    else
      ! Start with the previous location and search outward.
      if (unstable_int(K) .and. (rt >= ru_min_int(K)) .and. (rt <= ru_max_int(K))) then
        ! This interface started in an unstable region and should not move due to remapping.
        z_col_new(K) = z_col(K) ; reliable(K) = .false.
        k1_min = K ; k_found = .true.
      elseif ((rt >= ppoly_i_E(k-1,2)) .and. (rt <= ppoly_i_E(k,1))) then
        ! This interface is already in the right place and does not move.
        z_col_new(K) = z_col(K) ; reliable(K) = .true.
        k1_min = K ; k_found = .true.
      elseif (rt < ppoly_i_E(k-1,2)) then   ! Search upward
        do k1=K-1,k1_min,-1
          ! Check whether rt is in layer k.
          if ((rt < ppoly_i_E(k1,2)) .and. (rt > ppoly_i_E(k1,1))) then
            ! rt is in layer k.
            k_layer = k1
            k1_min = k1 ; k_found = .true. ; exit
          elseif (unstable_lay(k1) .and. (rt >= ru_min_lay(k1)) .and. (rt <= ru_max_lay(K1))) then
            ! rt would be found at unstable layer that it can not penetrate.
            !   It is possible that this can never happen?
            z_col_new(K) = z_col(K1+1) ; reliable(K) = .false.
            k1_min = k1 ; k_found = .true. ; exit
          endif
          ! Check whether rt is at interface K.
          if (k1 > 1) then ; if ((rt <= ppoly_i_E(k1,1)) .and. (rt >= ppoly_i_E(k1-1,2))) then
            ! rt is at interface K1
            z_col_new(K) = z_col(K1) ; reliable(K) = .true.
            k1_min = k1 ; k_found = .true. ; exit
          elseif (unstable_int(K1) .and. (rt >= ru_min_int(k1)) .and. (rt <= ru_max_int(K1))) then
            ! rt would be found at an unstable interface that it can not pass.
            !   It is possible that this can never happen?
            z_col_new(K) = z_col(K1) ; reliable(K) = .false.
            k1_min = k1 ; k_found = .true. ; exit
          endif ; endif
        enddo

        if (.not.k_found) then
          ! This should not happen unless k1_min = 1.
          if (k1_min < 2) then
            z_col_new(K) = z_col(k1_min)
          else
            z_col_new(K) = z_col(k1_min)
          endif
        endif

      else  ! Search downward
        do k1=K,nz
          if ((rt < ppoly_i_E(k1,2)) .and. (rt > ppoly_i_E(k1,1))) then
            ! rt is in layer k.
            k_layer = k1
            k1_min = k1 ; k_found = .true. ; exit
          elseif (unstable_lay(k1) .and. (rt >= ru_min_lay(k1)) .and. (rt <= ru_max_lay(K1))) then
            ! rt would be found at unstable layer that it can not penetrate.
            !   It is possible that this can never happen?
            z_col_new(K) = z_col(K1)
            reliable(K) = .false.
            k1_min = k1 ; k_found = .true. ; exit
          endif
          if (k1 < nz) then ; if ((rt <= ppoly_i_E(k1+1,1)) .and. (rt >= ppoly_i_E(k1,2))) then
            ! rt is at interface K1+1

            z_col_new(K) = z_col(K1+1) ; reliable(K) = .true.
            k1_min = k1+1 ; k_found = .true. ; exit
          elseif (unstable_int(K1+1) .and. (rt >= ru_min_int(k1+1)) .and. (rt <= ru_max_int(K1+1))) then
            ! rt would be found at an unstable interface that it can not pass.
            !   It is possible that this can never happen?
            z_col_new(K) = z_col(K1+1)
            reliable(K) = .false.
            k1_min = k1+1 ; k_found = .true. ; exit
          endif ; endif
        enddo
        if (.not.k_found) then
          z_col_new(K) = z_col(nz+1)
          if (rt >= ppoly_i_E(nz,2)) then
            reliable(K) = .true.
          else
            reliable(K) = .false.
          endif
        endif
      endif

      if (k_layer > 0) then  ! The new location is inside of layer k_layer.
        ! Note that this is coded assuming that this layer is stably stratified.
        if (.not.(ppoly_i_E(k1,2) > ppoly_i_E(k1,1))) call MOM_error(FATAL, &
          "build_grid_SLight: Erroneously searching for an interface in an unstratified layer.") !### COMMENT OUT LATER?

        ! Use the false position method to find the location (degree <= 1) or the first guess.
        zf = (rt - ppoly_i_E(k1,1)) / (ppoly_i_E(k1,2) - ppoly_i_E(k1,1))

        if (ppoly_degree > 1) then ! Iterate to find the solution.
          a(:) = 0.0 ; a(1) = ppoly_i_coefficients(k_layer,1) - rt
          do m=2,ppoly_degree+1 ; a(m) = ppoly_i_coefficients(k_layer,m) ; enddo
          ! Bracket the root.
          zf1 = 0.0 ; rfn1 = a(1)
          zf2 = 1.0 ; rfn2 =  a(1) + (a(2) + (a(3) + (a(4) + a(5))))
          if (rfn1 * rfn2 > 0.0) call MOM_error(FATAL, "build_grid_SLight: Bad bracketing.") !### COMMENT OUT LATER?

          do itt=1,max_itt
            rfn = a(1) + zf*(a(2) + zf*(a(3) + zf*(a(4) + zf*a(5))))
            ! Reset one of the ends of the bracket.
            if (rfn * rfn1 > 0.0) then
              zf1 = zf ; rfn1 = rfn
            else
              zf2 = zf ; rfn2 = rfn
            endif
            if (rfn1 == rfn2) exit

            drfn_dzf = (a(2) + zf*(2.0*a(3) + zf*(3.0*a(4) + zf*4.0*a(5))))
            sgn = 1.0 ; if (drfn_dzf < 0.0) sgn = -1.0

            if ((sgn*(zf - rfn) >= zf1 * abs(drfn_dzf)) .and. &
                (sgn*(zf - rfn) <= zf2 * abs(drfn_dzf))) then
              delta_zf = -rfn / drfn_dzf
              zf = zf + delta_zf
            else ! Newton's method goes out of bounds, so use a false position method estimate
              zf_prev = zf
              zf = ( rfn2 * zf1 - rfn1 * zf2 ) / (rfn2 - rfn1)
              delta_zf = zf - zf_prev
            endif

            if (abs(delta_zf) < tol) exit
          enddo
        endif
        z_col_new(K) = z_col(k_layer) + zf * z_sgn * h_col(k_layer)
        reliable(K) = .true.
      endif

    endif

  enddo
  z_col_new(nz+1) = z_col(nz+1) ; reliable(nz+1) = .true.

end subroutine rho_interfaces_col

!> Adjust dz_Interface to ensure non-negative future thicknesses
subroutine adjust_interface_motion( nk, min_thickness, h_old, dz_int )
  integer,               intent(in)    :: nk !< Number of layers
  real,                  intent(in)    :: min_thickness !< Minium allowed thickness of h (H units)
  real, dimension(nk),   intent(in)    :: h_old !< Minium allowed thickness of h (H units)
  real, dimension(nk+1), intent(inout) :: dz_int !< Minium allowed thickness of h (H units)
  ! Local variables
  integer :: k
  real :: h_new, eps, h_total, h_err

  eps = 1. ; eps = epsilon(eps)

  h_total = 0. ; h_err = 0.
  do k = 1, nk
    h_total = h_total + h_old(k)
    h_err = h_err + max( h_old(k), abs(dz_int(k)), abs(dz_int(k+1)) )*eps
    h_new = h_old(k) + ( dz_int(k) - dz_int(k+1) )
    if (h_new < -3.0*h_err) then
      write(0,*) 'h<0 at k=',k,'h_old=',h_old(k), &
        'wup=',dz_int(k),'wdn=',dz_int(k+1),'dw_dz=',dz_int(k) - dz_int(k+1), &
        'h_new=',h_new,'h_err=',h_err
      call MOM_error( FATAL, 'MOM_regridding: adjust_interface_motion() - '//&
                     'implied h<0 is larger than roundoff!')
    endif
  enddo
  do k = nk,2,-1
    h_new = h_old(k) + ( dz_int(k) - dz_int(k+1) )
    if (h_new<min_thickness) dz_int(k) = ( dz_int(k+1) - h_old(k) ) + min_thickness ! Implies next h_new = min_thickness
    h_new = h_old(k) + ( dz_int(k) - dz_int(k+1) )
    if (h_new<0.) dz_int(k) = ( 1. - eps ) * ( dz_int(k+1) - h_old(k) ) ! Backup in case min_thickness==0
    h_new = h_old(k) + ( dz_int(k) - dz_int(k+1) )
    if (h_new<0.) then
      write(0,*) 'h<0 at k=',k,'h_old=',h_old(k), &
        'wup=',dz_int(k),'wdn=',dz_int(k+1),'dw_dz=',dz_int(k) - dz_int(k+1), &
        'h_new=',h_new
      stop 'Still did not work!'
      call MOM_error( FATAL, 'MOM_regridding: adjust_interface_motion() - '//&
                     'Repeated adjustment for roundoff h<0 failed!')
    endif
  enddo
 !if (dz_int(1)/=0.) stop 'MOM_regridding: adjust_interface_motion() surface moved'

end subroutine adjust_interface_motion

!------------------------------------------------------------------------------
! Build arbitrary grid
!------------------------------------------------------------------------------
subroutine build_grid_arbitrary( G, GV, h, dzInterface, h_new, CS )
!------------------------------------------------------------------------------
! This routine builds a grid based on arbitrary rules
!------------------------------------------------------------------------------

  ! Arguments
  type(ocean_grid_type),                        intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                      intent(in)    :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)),   intent(in)    :: h  !< Original ayer thicknesses, in H
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)+1), intent(inout) :: dzInterface !< The change in interface depth in H
  real,                                         intent(inout) :: h_new !< New layer thicknesses, in H
  type(regridding_CS),                          intent(in)    :: CS !< Regridding control structure

  ! Local variables
  integer   :: i, j, k
  integer   :: nz
  real      :: z_inter(SZK_(GV)+1)
  real      :: total_height
  real      :: delta_h
  real      :: max_depth
  real      :: min_thickness
  real      :: eta              ! local elevation
  real      :: local_depth
  real      :: x1, y1, x2, y2
  real      :: x, t

  nz = GV%ke

  max_depth = G%max_depth
  min_thickness = CS%min_thickness

  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      ! Local depth
      local_depth = G%bathyT(i,j)*GV%m_to_H

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

        if ( x <= x1 ) then
          t = y1*x/x1
        else if ( (x > x1 ) .and. ( x < x2 )) then
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
        if ( z_inter(k) < (z_inter(k+1) + min_thickness) ) then
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


!> Given the set of target values and cell densities, this routine
!! builds an interpolated profile for the densities within each grid cell.
!! It may happen that, given a high-order interpolator, the number of
!! available layers is insufficient (e.g., there are two available layers for
!! a third-order PPM ih4 scheme). In these cases, we resort to the simplest
!! continuous linear scheme (P1M h2).
subroutine regridding_set_ppolys( densities, CS, n0, h0, ppoly0_E, ppoly0_S, &
                                  ppoly0_coefficients, degree)

  real, dimension(:),  intent(in)    :: densities !< Actual cell densities
  integer,             intent(in)    :: n0 !< Number of cells on source grid
  real, dimension(:),  intent(in)    :: h0 !< cell widths on source grid
  real, dimension(:,:),intent(inout) :: ppoly0_E            !< Edge value of polynomial
  real, dimension(:,:),intent(inout) :: ppoly0_S            !< Edge slope of polynomial
  real, dimension(:,:),intent(inout) :: ppoly0_coefficients !< Coefficients of polynomial
  integer,             intent(inout) :: degree  !< The degree of the polynomials
  type(regridding_CS), intent(in)    :: CS !< Parameters used for regridding

  ! Reset piecewise polynomials
  ppoly0_E(:,:) = 0.0
  ppoly0_S(:,:) = 0.0
  ppoly0_coefficients(:,:) = 0.0

  ! Compute the interpolated profile of the density field and build grid
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
      if ( n0 >= 4 ) then
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
      if ( n0 >= 4 ) then
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
      if ( n0 >= 4 ) then
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

      if ( n0 >= 4 ) then
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

      if ( n0 >= 4 ) then
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
      if ( n0 >= 6 ) then
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

      if ( n0 >= 4 ) then
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
      if ( n0 >= 6 ) then
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

end subroutine regridding_set_ppolys


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
    x1(k) = get_polynomial_coordinate ( n0, h0, x0, ppoly0_E, ppoly0_coefficients, t, degree )
    h1(k-1) = x1(k) - x1(k-1)
  end do
  h1(n1) = x1(n1+1) - x1(n1)

end subroutine interpolate_grid


!------------------------------------------------------------------------------
! Given target value, find corresponding coordinate for given polynomial
!------------------------------------------------------------------------------
function get_polynomial_coordinate ( N, h, x_g, ppoly_E, ppoly_coefficients, &
                                     target_value, degree ) result ( x_tgt )
! ------------------------------------------------------------------------------
! Here, 'ppoly' is assumed to be a piecewise discontinuous polynomial of degree
! 'degree' throughout the domain defined by 'grid'. A target value is given
! and we need to determine the corresponding grid coordinate to define the
! new grid.
!
! If the target value is out of range, the grid coordinate is simply set to
! be equal to one of the boundary coordinates, which results in vanished layers
! near the boundaries.
!
! IT IS ASSUMED THAT THE PIECEWISE POLYNOMIAL IS MONOTONICALLY INCREASING.
! IF THIS IS NOT THE CASE, THE NEW GRID MAY BE ILL-DEFINED.
!
! It is assumed that the number of cells defining 'grid' and 'ppoly' are the
! same.
! ------------------------------------------------------------------------------

  ! Arguments
  integer,              intent(in) :: N     ! The number of grid cells
  real, dimension(:),   intent(in) :: h     ! Grid cell thicknesses    (size N)
  real, dimension(:),   intent(in) :: x_g   ! Grid interface locations (size N+1)
  real, dimension(:,:), intent(in) :: ppoly_E  !Edge value of polynomial
  real, dimension(:,:), intent(in) :: ppoly_coefficients !Coefficients of polynomial
  real,                 intent(in) :: target_value
  integer,              intent(in) :: degree ! The degree of the polynomials

  real :: x_tgt      !< The position of x_g at which target_value is found.

  ! Local variables
  integer            :: i, k            ! loop indices
  integer            :: k_found         ! index of target cell
  integer            :: iter
  real               :: xi0             ! normalized target coordinate
  real, dimension(5) :: a               ! polynomial coefficients
  real               :: numerator
  real               :: denominator
  real               :: delta           ! Newton-Raphson increment
  real               :: x               ! global target coordinate
  real               :: eps                 ! offset used to get away from
                                        ! boundaries
  real               :: grad            ! gradient during N-R iterations

  eps = NR_OFFSET

  k_found = -1

  ! If the target value is outside the range of all values, we
  ! force the target coordinate to be equal to the lowest or
  ! largest value, depending on which bound is overtaken
  if ( target_value <= ppoly_E(1,1) ) then
    x_tgt = x_g(1)
    return  ! return because there is no need to look further
  end if

  ! Since discontinuous edge values are allowed, we check whether the target
  ! value lies between two discontinuous edge values at interior interfaces
  do k = 2,N
    if ( ( target_value >= ppoly_E(k-1,2) ) .AND. &
      ( target_value <= ppoly_E(k,1) ) ) then
      x_tgt = x_g(k)
      return   ! return because there is no need to look further
      exit
    end if
  end do

  ! If the target value is outside the range of all values, we
  ! force the target coordinate to be equal to the lowest or
  ! largest value, depending on which bound is overtaken
  if ( target_value >= ppoly_E(N,2) ) then
    x_tgt = x_g(N+1)
    return  ! return because there is no need to look further
  end if

  ! At this point, we know that the target value is bounded and does not
  ! lie between discontinuous, monotonic edge values. Therefore,
  ! there is a unique solution. We loop on all cells and find which one
  ! contains the target value. The variable k_found holds the index value
  ! of the cell where the taregt value lies.
  do k = 1,N
    if ( ( target_value > ppoly_E(k,1) ) .AND. &
         ( target_value < ppoly_E(k,2) ) ) then
      k_found = k
      exit
    end if
  end do

  ! At this point, 'k_found' should be strictly positive. If not, this is
  ! a major failure because it means we could not find any target cell
  ! despite the fact that the target value lies between the extremes. It
  ! means there is a major problem with the interpolant. This needs to be
  ! reported.
  if ( k_found == -1 ) then
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
    if ( ( iter > NR_ITERATIONS ) .OR. &
         ( abs(delta) < NR_TOLERANCE ) ) then
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
    if ( xi0 < 0.0 ) then
      xi0 = 0.0
      grad = a(2)
      if ( grad == 0.0 ) xi0 = xi0 + eps
    end if

    if ( xi0 > 1.0 ) then
      xi0 = 1.0
      grad = a(2) + 2*a(3) + 3*a(4) + 4*a(5)
      if ( grad == 0.0 ) xi0 = xi0 - eps
    end if

    iter = iter + 1

  end do ! end Newton-Raphson iterations

  x_tgt = x_g(k_found) + xi0 * h(k_found)

end function get_polynomial_coordinate


!------------------------------------------------------------------------------
! Check grid integrity
!------------------------------------------------------------------------------
subroutine inflate_vanished_layers_old( CS, G, GV, h )
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
  type(verticalGrid_type),                intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)), intent(inout) :: h

  ! Local variables
  integer :: i, j, k
  real    :: hTmp(GV%ke)

  do i = G%isc-1,G%iec+1
    do j = G%jsc-1,G%jec+1

      ! Build grid for current column
      do k = 1,GV%ke
        hTmp(k) = h(i,j,k)
      end do

      call old_inflate_layers_1d( CS%min_thickness, GV%ke, hTmp )

      ! Save modified grid
      do k = 1,GV%ke
        h(i,j,k) = hTmp(k)
      end do

    end do
  end do

end subroutine inflate_vanished_layers_old

!------------------------------------------------------------------------------
! Inflate vanished layers to finite (nonzero) width
!------------------------------------------------------------------------------
subroutine old_inflate_layers_1d( minThickness, N, h )

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
    if ( h(k) > minThickness ) then
      count_nonzero_layers = count_nonzero_layers + 1
    end if
  end do

  ! If all layer thicknesses are greater than the threshold, exit routine
  if ( count_nonzero_layers == N ) return

  ! If all thicknesses are zero, inflate them all and exit
  if ( count_nonzero_layers == 0 ) then
    do k = 1,N
      h(k) = minThickness
    end do
    return
  end if

  ! Inflate zero layers
  correction = 0.0
  do k = 1,N
    if ( h(k) <= minThickness ) then
      delta = minThickness - h(k)
      correction = correction + delta
      h(k) = h(k) + delta
    end if
  end do

  ! Modify thicknesses of nonzero layers to ensure volume conservation
  maxThickness = h(1)
  k_found = 1
  do k = 1,N
    if ( h(k) > maxThickness ) then
      maxThickness = h(k)
      k_found = k
    end if
  end do

  h(k_found) = h(k_found) - correction

end subroutine old_inflate_layers_1d


!------------------------------------------------------------------------------
! Convective adjustment by swapping layers
!------------------------------------------------------------------------------
subroutine convective_adjustment(G, GV, h, tv)
!------------------------------------------------------------------------------
! Check each water column to see whether it is stratified. If not, sort the
! layers by successive swappings of water masses (bubble sort algorithm)
!------------------------------------------------------------------------------

  ! Arguments
  type(ocean_grid_type), intent(in)                  :: G
  type(verticalGrid_type), intent(in)                :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h
  type(thermo_var_ptrs), intent(inout)               :: tv

  ! Local variables
  integer   :: i, j, k
  real      :: T0, T1       ! temperatures
  real      :: S0, S1       ! salinities
  real      :: r0, r1       ! densities
  real      :: h0, h1
  logical   :: stratified
  real, dimension(GV%ke) :: p_col, densities

  p_col(:) = 0.

  ! Loop on columns
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1

    ! Compute densities within current water column
    call calculate_density( tv%T(i,j,:), tv%S(i,j,:), p_col, &
                            densities, 1, GV%ke, tv%eqn_of_state )

    ! Repeat restratification until complete
    do
      stratified = .true.
      do k = 1,GV%ke-1
        ! Gather information of current and next cells
        T0 = tv%T(i,j,k)  ; T1 = tv%T(i,j,k+1)
        S0 = tv%S(i,j,k)  ; S1 = tv%S(i,j,k+1)
        r0 = densities(k) ; r1 = densities(k+1)
        h0 = h(i,j,k) ; h1 = h(i,j,k+1)
        ! If the density of the current cell is larger than the density
        ! below it, we swap the cells and recalculate the densitiies
        ! within the swapped cells
        if ( r0 > r1 ) then
          tv%T(i,j,k) = T1 ; tv%T(i,j,k+1) = T0
          tv%S(i,j,k) = S1 ; tv%S(i,j,k+1) = S0
          h(i,j,k)    = h1 ; h(i,j,k+1)    = h0
          ! Recompute densities at levels k and k+1
          call calculate_density( tv%T(i,j,k), tv%S(i,j,k), p_col(k), &
                                  densities(k), tv%eqn_of_state )
          call calculate_density( tv%T(i,j,k+1), tv%S(i,j,k+1), p_col(k+1), &
                                  densities(k+1), tv%eqn_of_state )
          stratified = .false.
        end if
      enddo  ! k

      if ( stratified ) exit
    enddo

  enddo ; enddo  ! i & j

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

    case ( REGRIDDING_ZSTAR, REGRIDDING_HYCOM1, REGRIDDING_SLIGHT, REGRIDDING_SIGMA_SHELF_ZSTAR )
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

!> Set target densities based on the old Rlay variable
subroutine set_target_densities_from_GV( GV, CS )
  type(verticalGrid_type), intent(in)    :: GV !< Ocean vertical grid structure
  type(regridding_CS),     intent(inout) :: CS !< Regridding control structure
  ! Local variables
  integer :: k, nz

  nz = CS%nk
  CS%target_density(1)    = GV%Rlay(1)+0.5*(GV%Rlay(1)-GV%Rlay(2))
  CS%target_density(nz+1) = GV%Rlay(nz)+0.5*(GV%Rlay(nz)-GV%Rlay(nz-1))
  do k = 2,nz
    CS%target_density(k) = CS%target_density(k-1) + CS%coordinateResolution(k)
  end do
  CS%target_density_set = .true.

end subroutine set_target_densities_from_GV

!> Set target densities based on vector of interface values
subroutine set_target_densities( CS, rho_int )
  type(regridding_CS),      intent(inout) :: CS !< Regridding control structure
  real, dimension(CS%nk+1), intent(in)    :: rho_int !< Interface densities

  CS%target_density(:) = rho_int(:)
  CS%target_density_set = .true.

end subroutine set_target_densities

!> Set maximum interface depths based on a vector of input values.
subroutine set_regrid_max_depths( CS, max_depths, units_to_H )
  type(regridding_CS),      intent(inout) :: CS !< Regridding control structure
  real, dimension(CS%nk+1), intent(in)    :: max_depths !< Maximum interface depths, in arbitrary units
  real, optional,           intent(in)    :: units_to_H !< A conversion factor for max_depths into H units
  ! Local variables
  real :: val_to_H
  integer :: K

  if (.not.allocated(CS%max_interface_depths)) allocate(CS%max_interface_depths(1:CS%nk+1))

  val_to_H = 1.0 ; if (present( units_to_H)) val_to_H = units_to_H
  if (max_depths(CS%nk+1) < max_depths(1)) val_to_H = -1.0*val_to_H

  ! Check for sign reversals in the depths.
  if (max_depths(CS%nk+1) < max_depths(1)) then
    do K=1,CS%nk ; if (max_depths(K+1) > max_depths(K)) &
      call MOM_error(FATAL, "Unordered list of maximum depths sent to set_regrid_max_depths!")
    enddo
  else
    do K=1,CS%nk ; if (max_depths(K+1) < max_depths(K)) &
      call MOM_error(FATAL, "Unordered list of maximum depths sent to set_regrid_max_depths.")
    enddo
  endif

  do K=1,CS%nk+1
    CS%max_interface_depths(K) = val_to_H * max_depths(K)
  enddo

end subroutine set_regrid_max_depths

!> Set maximum layer thicknesses based on a vector of input values.
subroutine set_regrid_max_thickness( CS, max_h, units_to_H )
  type(regridding_CS),      intent(inout) :: CS !< Regridding control structure
  real, dimension(CS%nk+1), intent(in)    :: max_h !< Maximum interface depths, in arbitrary units
  real, optional,           intent(in)    :: units_to_H !< A conversion factor for max_h into H units
  ! Local variables
  real :: val_to_H
  integer :: K

  if (.not.allocated(CS%max_layer_thickness)) allocate(CS%max_layer_thickness(1:CS%nk))

  val_to_H = 1.0 ; if (present( units_to_H)) val_to_H = units_to_H

  do k=1,CS%nk
    CS%max_layer_thickness(k) = val_to_H * max_h(k)
  enddo

end subroutine set_regrid_max_thickness


!------------------------------------------------------------------------------
! Query the fixed resolution data
!------------------------------------------------------------------------------
function getCoordinateResolution( CS )
  type(regridding_CS), intent(in) :: CS
  real, dimension(CS%nk)          :: getCoordinateResolution

  getCoordinateResolution(:) = CS%coordinateResolution(:)

end function getCoordinateResolution

!> Query the target coordinate interface positions
function getCoordinateInterfaces( CS )
  type(regridding_CS), intent(in) :: CS                      !< Regridding control structure
  real, dimension(CS%nk+1)        :: getCoordinateInterfaces !< Interface positions in target coordinate

  integer :: k

  ! When using a coordinate with target densities, we need to get the actual
  ! densities, rather than computing the interfaces based on resolution
  if (CS%regridding_scheme == REGRIDDING_RHO) then
    if (.not. CS%target_density_set) &
      call MOM_error(FATAL, 'MOM_regridding, getCoordinateInterfaces: '//&
                            'target densities not set!')

    getCoordinateInterfaces(:) = CS%target_density(:)
  else
    getCoordinateInterfaces(1) = 0.
    do k = 1, CS%nk
      getCoordinateInterfaces(k+1) = getCoordinateInterfaces(k) &
                                    -CS%coordinateResolution(k)
    enddo
    ! The following line has an "abs()" to allow ferret users to reference
    ! data by index. It is a temporary work around...  :(  -AJA
    getCoordinateInterfaces(:) = abs( getCoordinateInterfaces(:) )
  end if

end function getCoordinateInterfaces

!------------------------------------------------------------------------------
! Query the target coordinate units
!------------------------------------------------------------------------------
function getCoordinateUnits( CS )
  type(regridding_CS), intent(in) :: CS
  character(len=20)               :: getCoordinateUnits

  select case ( CS%regridding_scheme )
    case ( REGRIDDING_ZSTAR, REGRIDDING_HYCOM1, REGRIDDING_SLIGHT )
      getCoordinateUnits = 'meter'
    case ( REGRIDDING_SIGMA_SHELF_ZSTAR )
      getCoordinateUnits = 'meter/fraction'
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
      getCoordinateShortName = 'pseudo-depth, -z*'
    case ( REGRIDDING_SIGMA_SHELF_ZSTAR )
      getCoordinateShortName = 'pseudo-depth, -z*/sigma'
    case ( REGRIDDING_SIGMA )
      getCoordinateShortName = 'sigma'
    case ( REGRIDDING_RHO )
      getCoordinateShortName = 'rho'
    case ( REGRIDDING_ARBITRARY )
      getCoordinateShortName = 'coordinate'
    case ( REGRIDDING_HYCOM1 )
      getCoordinateShortName = 'z-rho'
    case ( REGRIDDING_SLIGHT )
      getCoordinateShortName = 's-rho'
    case default
      call MOM_error(FATAL,'MOM_regridding, getCoordinateShortName: '//&
                     'Unknown regridding scheme selected!')
  end select ! type of grid

end function getCoordinateShortName

!> Can be used to set any of the parameters for MOM_regridding.
subroutine set_regrid_params( CS, boundary_extrapolation, min_thickness, old_grid_weight, &
             interp_scheme, depth_of_time_filter_shallow, depth_of_time_filter_deep, &
             compress_fraction, dz_min_surface, nz_fixed_surface, Rho_ML_avg_depth, &
             nlay_ML_to_interior, fix_haloclines, halocline_filt_len, &
             halocline_strat_tol, integrate_downward_for_e)
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure
  logical, optional, intent(in) :: boundary_extrapolation !< Extrapolate in boundary cells
  real,    optional, intent(in) :: min_thickness !< Minimum thickness allowed when building the new grid (m)
  real,    optional, intent(in) :: old_grid_weight !< Weight given to old coordinate when time-filtering grid
  character(len=*), optional, intent(in) :: interp_scheme !< Interpolation method for state-dependent coordinates
  real,    optional, intent(in) :: depth_of_time_filter_shallow !< Depth to start cubic (H units)
  real,    optional, intent(in) :: depth_of_time_filter_deep !< Depth to end cubic (H units)
  real,    optional, intent(in) :: compress_fraction !< Fraction of compressibility to add to potential density
  real,    optional, intent(in) :: dz_min_surface !< The fixed resolution in the topmost SLight_nkml_min layers (m)
  integer, optional, intent(in) :: nz_fixed_surface !< The number of fixed-thickess layers at the top of the model
  real,    optional, intent(in) :: Rho_ml_avg_depth !< Averaging depth over which to determine mixed layer potential density (m)
  real,    optional, intent(in) :: nlay_ML_to_interior !< Number of layers to offset the mixed layer density to find resolved stratification (nondim)
  logical, optional, intent(in) :: fix_haloclines !< Detect regions with much weaker stratification in the coordinate
  real,    optional, intent(in) :: halocline_filt_len !< Length scale over which to filter T & S when looking for spuriously unstable water mass profiles (m)
  real,    optional, intent(in) :: halocline_strat_tol !< Value of the stratification ratio that defines a problematic halocline region.
  logical, optional, intent(in) :: integrate_downward_for_e !< If true, integrate for interface positions downward from the top.

  if (present(boundary_extrapolation)) CS%boundary_extrapolation = boundary_extrapolation
  if (present(min_thickness)) CS%min_thickness = min_Thickness
  if (present(old_grid_weight)) then
    if (old_grid_weight<0. .or. old_grid_weight>1.) &
      call MOM_error(FATAL,'MOM_regridding, set_regrid_params: Weight is out side the range 0..1!')
    CS%old_grid_weight = old_grid_weight
  endif
  if (present(interp_scheme)) CS%interpolation_scheme = interpolation_scheme(interp_scheme)
  if (present(depth_of_time_filter_shallow)) CS%depth_of_time_filter_shallow = depth_of_time_filter_shallow
  if (present(depth_of_time_filter_deep)) CS%depth_of_time_filter_deep = depth_of_time_filter_deep
  if (present(depth_of_time_filter_shallow) .or. present(depth_of_time_filter_deep)) then
    if (CS%depth_of_time_filter_deep<CS%depth_of_time_filter_shallow) call MOM_error(FATAL,'MOM_regridding, '//&
                     'set_regrid_params: depth_of_time_filter_deep<depth_of_time_filter_shallow!')
  endif
  if (present(compress_fraction)) CS%compressibility_fraction = compress_fraction
  if (present(dz_min_surface)) CS%dz_ml_min = dz_min_surface
  if (present(nz_fixed_surface)) CS%nz_fixed_surface = nz_fixed_surface
  if (present(Rho_ML_avg_depth)) CS%Rho_ML_avg_depth = Rho_ML_avg_depth
  if (present(nlay_ML_to_interior)) CS%nlay_ML_offset = nlay_ML_to_interior
  if (present(fix_haloclines)) CS%fix_haloclines = fix_haloclines
  if (present(halocline_filt_len)) CS%halocline_filter_length = halocline_filt_len
  if (present(halocline_strat_tol)) then
    if (halocline_strat_tol > 1.0) call MOM_error(FATAL, "set_regrid_params: "//&
        "HALOCLINE_STRAT_TOL must not exceed 1.0.")
    CS%halocline_strat_tol = halocline_strat_tol
  endif
  if (present(integrate_downward_for_e)) CS%integrate_downward_for_e = integrate_downward_for_e

end subroutine set_regrid_params

!> Returns the number of levels/layers in the regridding control structure
integer function get_regrid_size(CS)
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure

  get_regrid_size = CS%nk

end function get_regrid_size

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
    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA_SHELF_ZSTAR, REGRIDDING_HYCOM1, REGRIDDING_SLIGHT )
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

!> Parses a string and generates a dz(:) profile that goes like k**power.
subroutine dz_function1( string, dz )
  character(len=*),   intent(in)    :: string !< String with list of parameters in form
                                              !! dz_min, H_total, power, precision
  real, dimension(:), intent(inout) :: dz     !< Profile of nominal thicknesses
  ! Local variables
  integer :: nk, k
  real    :: dz_min, power, prec, H_total

  nk = size(dz) ! Number of cells
  prec = -1024.
  read( string, *) dz_min, H_total, power, prec
  if (prec == -1024.) call MOM_error(FATAL,"dz_function1: "// &
          "Problem reading FNC1: string  ="//trim(string))
  ! Create profile of ( dz - dz_min )
  do k = 1, nk
    dz(k) = (real(k-1)/real(nk-1))**power
  enddo
  dz(:) = ( H_total - real(nk) * dz_min ) * ( dz(:) / sum(dz) ) ! Rescale to so total is H_total
  dz(:) = anint( dz(:) / prec ) * prec ! Rounds to precision prec
  dz(:) = ( H_total - real(nk) * dz_min ) * ( dz(:) / sum(dz) ) ! Rescale to so total is H_total
  dz(:) = anint( dz(:) / prec ) * prec ! Rounds to precision prec
  dz(nk) = dz(nk) + ( H_total - sum( dz(:) + dz_min ) ) ! Adjust bottommost layer
  dz(:) = anint( dz(:) / prec ) * prec ! Rounds to precision prec
  dz(:) = dz(:) + dz_min ! Finally add in the constant dz_min

end subroutine dz_function1

!> \namespace mom_regridding
!!
!! A vertical grid is defined solely by the cell thicknesses, \f$h\f$.
!! Most calculations in this module start with the coordinate at the bottom
!! of the column set to -depth, and use a increasing value of coordinate with
!! decreasing k. This is consistent with the rest of MOM6 that uses position,
!! \f$z\f$ which is a negative quantity for most of the ocean.
!!
!! A change in grid is define through a change in position of the interfaces:
!! \f[
!! z^n_{k+1/2} = z^{n-1}_{k+1/2} + \Delta z_{k+1/2}
!! \f]
!! with the positive upward coordinate convention
!! \f[
!! z_{k-1/2} = z_{k+1/2} + h_k
!! \f]
!! so that
!! \f[
!! h^n_k = h^{n-1}_k + ( \Delta z_{k-1/2} - \Delta z_{k+1/2} )
!! \f]
!!
!! Original date of creation: 2008.06.09 by L. White

end module MOM_regridding
