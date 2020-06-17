!> Generates vertical grids as part of the ALE algorithm
module MOM_regridding

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser,   only : param_file_type, get_param, log_param
use MOM_io,            only : file_exists, field_exists, field_size, MOM_read_data
use MOM_io,            only : slasher
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS,           only : EOS_type, calculate_density
use MOM_string_functions, only : uppercase, extractWord, extract_integer, extract_real

use MOM_remapping, only : remapping_CS
use regrid_consts, only : state_dependent, coordinateUnits
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA
use regrid_consts, only : REGRIDDING_ARBITRARY, REGRIDDING_SIGMA_SHELF_ZSTAR
use regrid_consts, only : REGRIDDING_HYCOM1, REGRIDDING_SLIGHT, REGRIDDING_ADAPTIVE
use regrid_interp, only : interp_CS_type, set_interp_scheme, set_interp_extrap

use coord_zlike,  only : init_coord_zlike, zlike_CS, set_zlike_params, build_zstar_column, end_coord_zlike
use coord_sigma,  only : init_coord_sigma, sigma_CS, set_sigma_params, build_sigma_column, end_coord_sigma
use coord_rho,    only : init_coord_rho, rho_CS, set_rho_params, build_rho_column, end_coord_rho
use coord_rho,    only : old_inflate_layers_1d
use coord_hycom,  only : init_coord_hycom, hycom_CS, set_hycom_params, build_hycom1_column, end_coord_hycom
use coord_slight, only : init_coord_slight, slight_CS, set_slight_params, build_slight_column, end_coord_slight
use coord_adapt,  only : init_coord_adapt, adapt_CS, set_adapt_params, build_adapt_column, end_coord_adapt

use netcdf ! Used by check_grid_def()

implicit none ; private

#include <MOM_memory.h>

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Regridding control structure
type, public :: regridding_CS ; private

  !> This array is set by function setCoordinateResolution()
  !! It contains the "resolution" or delta coordinate of the target
  !! coordinate.  It has the units of the target coordinate, e.g.
  !! [Z ~> m] for z*, non-dimensional for sigma, etc.
  real, dimension(:), allocatable :: coordinateResolution

  !> This is a scaling factor that restores coordinateResolution to values in
  !! the natural units for output.
  real :: coord_scale = 1.0

  !> This array is set by function set_target_densities()
  !! This array is the nominal coordinate of interfaces and is the
  !! running sum of coordinateResolution, in [R ~> kg m-3]. i.e.
  !!  target_density(k+1) = coordinateResolution(k) + coordinateResolution(k)
  !! It is only used in "rho", "SLight" or "Hycom" mode.
  real, dimension(:), allocatable :: target_density

  !> A flag to indicate that the target_density arrays has been filled with data.
  logical :: target_density_set = .false.

  !> This array is set by function set_regrid_max_depths()
  !! It specifies the maximum depth that every interface is allowed to take [H ~> m or kg m-2].
  real, dimension(:), allocatable :: max_interface_depths

  !> This array is set by function set_regrid_max_thickness()
  !! It specifies the maximum depth that every interface is allowed to take [H ~> m or kg m-2].
  real, dimension(:), allocatable :: max_layer_thickness

  integer :: nk !< Number of layers/levels in generated grid

  !> Indicates which grid to use in the vertical (z*, sigma, target interface
  !! densities)
  integer :: regridding_scheme

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS

  !> Minimum thickness allowed when building the new grid through regridding [H ~> m or kg m-2].
  real :: min_thickness

  !> Reference pressure for potential density calculations [R L2 T-2 ~> Pa]
  real :: ref_pressure = 2.e7

  !> Weight given to old coordinate when blending between new and old grids [nondim]
  !! Used only below depth_of_time_filter_shallow, with a cubic variation
  !! from zero to full effect between depth_of_time_filter_shallow and
  !! depth_of_time_filter_deep.
  real :: old_grid_weight = 0.

  !> Depth above which no time-filtering of grid is applied [H ~> m or kg m-2]
  real :: depth_of_time_filter_shallow = 0.

  !> Depth below which time-filtering of grid is applied at full effect [H ~> m or kg m-2]
  real :: depth_of_time_filter_deep = 0.

  !> Fraction (between 0 and 1) of compressibility to add to potential density
  !! profiles when interpolating for target grid positions. [nondim]
  real :: compressibility_fraction = 0.

  !> If true, each interface is given a maximum depth based on a rescaling of
  !! the indexing of coordinateResolution.
  logical :: set_maximum_depths = .false.

  !> A scaling factor (> 1) of the rate at which the coordinateResolution list
  !! is traversed to set the minimum depth of interfaces.
  real :: max_depth_index_scale = 2.0

  !> If true, integrate for interface positions from the top downward.
  !! If false, integrate from the bottom upward, as does the rest of the model.
  logical :: integrate_downward_for_e = .true.

  !> If true, use the order of arithmetic and expressions that recover the remapping answers from 2018.
  !! If false, use more robust forms of the same remapping expressions.
  logical :: remap_answers_2018 = .true.

  type(zlike_CS),  pointer :: zlike_CS  => null() !< Control structure for z-like coordinate generator
  type(sigma_CS),  pointer :: sigma_CS  => null() !< Control structure for sigma coordinate generator
  type(rho_CS),    pointer :: rho_CS    => null() !< Control structure for rho coordinate generator
  type(hycom_CS),  pointer :: hycom_CS  => null() !< Control structure for hybrid coordinate generator
  type(slight_CS), pointer :: slight_CS => null() !< Control structure for Slight-coordinate generator
  type(adapt_CS),  pointer :: adapt_CS  => null() !< Control structure for adaptive coordinate generator

end type

! The following routines are visible to the outside world
public initialize_regridding, end_regridding, regridding_main
public inflate_vanished_layers_old, check_remapping_grid, check_grid_column
public set_regrid_params, get_regrid_size
public uniformResolution, setCoordinateResolution
public build_rho_column
public set_target_densities_from_GV, set_target_densities
public set_regrid_max_depths, set_regrid_max_thickness
public getCoordinateResolution, getCoordinateInterfaces
public getCoordinateUnits, getCoordinateShortName, getStaticThickness
public DEFAULT_COORDINATE_MODE
public get_zlike_CS, get_sigma_CS, get_rho_CS

!> Documentation for coordinate options
character(len=*), parameter, public :: regriddingCoordinateModeDoc = &
                 " LAYER - Isopycnal or stacked shallow water layers\n"//&
                 " ZSTAR, Z* - stretched geopotential z*\n"//&
                 " SIGMA_SHELF_ZSTAR - stretched geopotential z* ignoring shelf\n"//&
                 " SIGMA - terrain following coordinates\n"//&
                 " RHO   - continuous isopycnal\n"//&
                 " HYCOM1 - HyCOM-like hybrid coordinate\n"//&
                 " SLIGHT - stretched coordinates above continuous isopycnal\n"//&
                 " ADAPTIVE - optimize for smooth neutral density surfaces"

!> Documentation for regridding interpolation schemes
character(len=*), parameter, public :: regriddingInterpSchemeDoc = &
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

!> Default interpolation scheme
character(len=*), parameter, public :: regriddingDefaultInterpScheme = "P1M_H2"
!> Default mode for boundary extrapolation
logical, parameter, public :: regriddingDefaultBoundaryExtrapolation = .false.
!> Default minimum thickness for some coordinate generation modes
real, parameter, public :: regriddingDefaultMinThickness = 1.e-3

#undef __DO_SAFETY_CHECKS__

contains

!> Initialization and configures a regridding control structure based on customizable run-time parameters
subroutine initialize_regridding(CS, GV, US, max_depth, param_file, mdl, coord_mode, param_prefix, param_suffix)
  type(regridding_CS),        intent(inout) :: CS  !< Regridding control structure
  type(verticalGrid_type),    intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),      intent(in)    :: US  !< A dimensional unit scaling type
  real,                       intent(in)    :: max_depth  !< The maximum depth of the ocean [Z ~> m].
  type(param_file_type),      intent(in)    :: param_file !< Parameter file
  character(len=*),           intent(in)    :: mdl        !< Name of calling module.
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
  logical :: default_2018_answers, remap_answers_2018
  real :: filt_len, strat_tol, index_scale, tmpReal, P_Ref
  real :: maximum_depth ! The maximum depth of the ocean [m] (not in Z).
  real :: dz_fixed_sfc, Rho_avg_depth, nlay_sfc_int
  real :: adaptTimeRatio, adaptZoom, adaptZoomCoeff, adaptBuoyCoeff, adaptAlpha
  real :: adaptDrho0 ! Reference density difference for stratification-dependent diffusion. [R ~> kg m-3]
  integer :: nz_fixed_sfc, k, nzf(4)
  real, dimension(:), allocatable :: dz     ! Resolution (thickness) in units of coordinate, which may be [m]
                                            ! or [Z ~> m] or [H ~> m or kg m-2] or [R ~> kg m-3] or other units.
  real, dimension(:), allocatable :: h_max  ! Maximum layer thicknesses [H ~> m or kg m-2]
  real, dimension(:), allocatable :: z_max  ! Maximum interface depths [H ~> m or kg m-2] or other
                                            ! units depending on the coordinate
  real, dimension(:), allocatable :: dz_max ! Thicknesses used to find maximum interface depths
                                            ! [H ~> m or kg m-2] or other units
  real, dimension(:), allocatable :: rho_target ! Target density used in HYBRID mode [kg m-3]
  ! Thicknesses [m] that give level centers corresponding to table 2 of WOA09
  real, dimension(40) :: woa09_dz = (/ 5.,  10.,  10.,  15.,  22.5, 25., 25.,  25.,  &
                                      37.5, 50.,  50.,  75., 100., 100., 100., 100., &
                                     100., 100., 100., 100., 100., 100., 100., 175., &
                                     250., 375., 500., 500., 500., 500., 500., 500., &
                                     500., 500., 500., 500., 500., 500., 500., 500. /)

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)

  main_parameters=.false.
  if (len_trim(param_prefix)==0) main_parameters=.true.
  if (main_parameters .and. len_trim(param_suffix)>0) call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
              'Suffix provided without prefix for parameter names!')

  CS%nk = 0
  CS%regridding_scheme = coordinateMode(coord_mode)
  coord_is_state_dependent = state_dependent(coord_mode)
  maximum_depth = US%Z_to_m*max_depth

  if (main_parameters) then
    ! Read coordinate units parameter (main model = REGRIDDING_COORDINATE_UNITS)
    call get_param(param_file, mdl, "REGRIDDING_COORDINATE_UNITS", coord_units, &
                 "Units of the regridding coordinate.", default=coordinateUnits(coord_mode))
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
    call get_param(param_file, mdl, "INTERPOLATION_SCHEME", string, &
                 "This sets the interpolation scheme to use to "//&
                 "determine the new grid. These parameters are "//&
                 "only relevant when REGRIDDING_COORDINATE_MODE is "//&
                 "set to a function of state. Otherwise, it is not "//&
                 "used. It can be one of the following schemes: "//&
                 trim(regriddingInterpSchemeDoc), default=trim(string2))
    call set_regrid_params(CS, interp_scheme=string)

    call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false.)
    call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", remap_answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
    call set_regrid_params(CS, remap_answers_2018=remap_answers_2018)
  endif

  if (main_parameters .and. coord_is_state_dependent) then
    call get_param(param_file, mdl, "BOUNDARY_EXTRAPOLATION", tmpLogical, &
                 "When defined, a proper high-order reconstruction "//&
                 "scheme is used within boundary cells rather "//&
                 "than PCM. E.g., if PPM is used for remapping, a "//&
                 "PPM reconstruction will also be used within "//&
                 "boundary cells.", default=regriddingDefaultBoundaryExtrapolation)
    call set_regrid_params(CS, boundary_extrapolation=tmpLogical)
  else
    call set_regrid_params(CS, boundary_extrapolation=.false.)
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
    if (maximum_depth>3000.) string2='WOA09' ! For convenience
  endif
  call get_param(param_file, mdl, param_name, string, &
                 "Determines how to specify the coordinate "//&
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
      tmpReal = maximum_depth
    elseif (index(trim(string),'UNIFORM:')==1 .and. len_trim(string)>8) then
      ! Format is "UNIFORM:N" or "UNIFORM:N,dz"
      ke = extract_integer(string(9:len_trim(string)),'',1)
      tmpReal = extract_real(string(9:len_trim(string)),',',2,missing_value=maximum_depth)
    else
      call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
          'Unable to interpret "'//trim(string)//'".')
    endif
    allocate(dz(ke))
    dz(:) = uniformResolution(ke, coord_mode, tmpReal, &
                              US%R_to_kg_m3*(GV%Rlay(1) + 0.5*(GV%Rlay(1)-GV%Rlay(min(2,ke)))), &
                              US%R_to_kg_m3*(GV%Rlay(ke) + 0.5*(GV%Rlay(ke)-GV%Rlay(max(ke-1,1)))) )
    if (main_parameters) call log_param(param_file, mdl, "!"//coord_res_param, dz, &
                   trim(message), units=trim(coord_units))
  elseif (trim(string)=='PARAM') then
    ! Read coordinate resolution (main model = ALE_RESOLUTION)
    ke = GV%ke ! Use model nk by default
    allocate(dz(ke))
    call get_param(param_file, mdl, coord_res_param, dz, &
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
    if (.not. file_exists(fileName)) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
            "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

    varName = trim( extractWord(trim(string(6:)), 2) )
    if (len_trim(varName)==0) then
      if (field_exists(fileName,'dz')) then; varName = 'dz'
      elseif (field_exists(fileName,'dsigma')) then; varName = 'dsigma'
      elseif (field_exists(fileName,'ztest')) then; varName = 'ztest'
      else ;  call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
                    "Coordinate variable not specified and none could be guessed.")
      endif
    endif
    ! This check fails when the variable is a dimension variable! -AJA
   !if (.not. field_exists(fileName,trim(varName))) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
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
      if (ierr) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "//&
                  "Unsupported format in grid definition '"//trim(filename)//"'. Error message "//trim(message))
      call field_size(trim(fileName), trim(varName), nzf)
      ke = nzf(1)-1
      if (CS%regridding_scheme == REGRIDDING_RHO) then
        allocate(rho_target(ke+1))
        call MOM_read_data(trim(fileName), trim(varName), rho_target)
      else
        allocate(dz(ke))
        allocate(z_max(ke+1))
        call MOM_read_data(trim(fileName), trim(varName), z_max)
        dz(:) = abs(z_max(1:ke) - z_max(2:ke+1))
        deallocate(z_max)
      endif
    else
      ! Assume reading resolution
      call field_size(trim(fileName), trim(varName), nzf)
      ke = nzf(1)
      allocate(dz(ke))
      call MOM_read_data(trim(fileName), trim(varName), dz)
    endif
    if (main_parameters .and. ke/=GV%ke) then
      call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
                 'Mismatch in number of model levels and "'//trim(string)//'".')
    endif
    if (main_parameters) call log_param(param_file, mdl, "!"//coord_res_param, dz, &
               trim(message), units=coordinateUnits(coord_mode))
  elseif (index(trim(string),'FNC1:')==1) then
    ke = GV%ke; allocate(dz(ke))
    call dz_function1( trim(string(6:)), dz )
    if (main_parameters) call log_param(param_file, mdl, "!"//coord_res_param, dz, &
               trim(message), units=coordinateUnits(coord_mode))
  elseif (index(trim(string),'RFNC1:')==1) then
    ! Function used for set target interface densities
    ke = rho_function1( trim(string(7:)), rho_target )
  elseif (index(trim(string),'HYBRID:')==1) then
    ke = GV%ke; allocate(dz(ke))
    ! The following assumes the FILE: syntax of above but without "FILE:" in the string
    allocate(rho_target(ke+1))
    fileName = trim( extractWord(trim(string(8:)), 1) )
    if (fileName(1:1)/='.' .and. filename(1:1)/='/') fileName = trim(inputdir) // trim( fileName )
    if (.not. file_exists(fileName)) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: HYBRID "// &
      "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")
    varName = trim( extractWord(trim(string(8:)), 2) )
    if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: HYBRID "// &
      "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
    call MOM_read_data(trim(fileName), trim(varName), rho_target)
    varName = trim( extractWord(trim(string(8:)), 3) )
    if (varName(1:5) == 'FNC1:') then ! Use FNC1 to calculate dz
      call dz_function1( trim(string((index(trim(string),'FNC1:')+5):)), dz )
    else ! Read dz from file
      if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: HYBRID "// &
        "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
      call MOM_read_data(trim(fileName), trim(varName), dz)
    endif
    if (main_parameters) then
      call log_param(param_file, mdl, "!"//coord_res_param, dz, &
               trim(message), units=coordinateUnits(coord_mode))
      call log_param(param_file, mdl, "!TARGET_DENSITIES", rho_target, &
               'HYBRID target densities for interfaces', units=coordinateUnits(coord_mode))
    endif
  elseif (index(trim(string),'WOA09')==1) then
    if (len_trim(string)==5) then
      tmpReal = 0. ; ke = 0
      do while (tmpReal<maximum_depth)
        ke = ke + 1
        tmpReal = tmpReal + woa09_dz(ke)
      enddo
    elseif (index(trim(string),'WOA09:')==1) then
      if (len_trim(string)==6) call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
                 'Expected string of form "WOA09:N" but got "'//trim(string)//'".')
      ke = extract_integer(string(7:len_trim(string)),'',1)
    endif
    if (ke>40 .or. ke<1) call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
                 'For "WOA05:N" N must 0<N<41 but got "'//trim(string)//'".')
    allocate(dz(ke))
    dz(1:ke) = woa09_dz(1:ke)
  else
    call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
      "Unrecognized coordinate configuration"//trim(string))
  endif

  if (main_parameters) then
    ! This is a work around to apparently needed to work with the from_Z initialization...  ???
    if (coordinateMode(coord_mode) == REGRIDDING_ZSTAR .or. &
        coordinateMode(coord_mode) == REGRIDDING_HYCOM1 .or. &
        coordinateMode(coord_mode) == REGRIDDING_SLIGHT .or. &
        coordinateMode(coord_mode) == REGRIDDING_ADAPTIVE) then
      ! Adjust target grid to be consistent with maximum_depth
      tmpReal = sum( dz(:) )
      if (tmpReal < maximum_depth) then
        dz(ke) = dz(ke) + ( maximum_depth - tmpReal )
      elseif (tmpReal > maximum_depth) then
        if ( dz(ke) + ( maximum_depth - tmpReal ) > 0. ) then
          dz(ke) = dz(ke) + ( maximum_depth - tmpReal )
        else
          call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
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
    allocate( CS%target_density(CS%nk+1) ); CS%target_density(:) = -1.E30*US%kg_m3_to_R
  endif

  if (allocated(dz)) then
    if (coordinateMode(coord_mode) == REGRIDDING_SIGMA) then
      call setCoordinateResolution(dz, CS, scale=1.0)
    elseif (coordinateMode(coord_mode) == REGRIDDING_RHO) then
      call setCoordinateResolution(dz, CS, scale=US%kg_m3_to_R)
      CS%coord_scale = US%R_to_kg_m3
    elseif (coordinateMode(coord_mode) == REGRIDDING_ADAPTIVE) then
      call setCoordinateResolution(dz, CS, scale=GV%m_to_H)
      CS%coord_scale = GV%H_to_m
    else
      call setCoordinateResolution(dz, CS, scale=US%m_to_Z)
      CS%coord_scale = US%Z_to_m
    endif
  endif

  if (allocated(rho_target)) then
    call set_target_densities(CS, US%kg_m3_to_R*rho_target)
    deallocate(rho_target)

  ! \todo This line looks like it would overwrite the target densities set just above?
  elseif (coordinateMode(coord_mode) == REGRIDDING_RHO) then
    call set_target_densities_from_GV(GV, US, CS)
    call log_param(param_file, mdl, "!TARGET_DENSITIES", US%R_to_kg_m3*CS%target_density(:), &
             'RHO target densities for interfaces', units=coordinateUnits(coord_mode))
  endif

  ! initialise coordinate-specific control structure
  call initCoord(CS, GV, US, coord_mode)

  if (main_parameters .and. coord_is_state_dependent) then
    call get_param(param_file, mdl, "P_REF", P_Ref, &
                 "The pressure that is used for calculating the coordinate "//&
                 "density.  (1 Pa = 1e4 dbar, so 2e7 is commonly used.) "//&
                 "This is only used if USE_EOS and ENABLE_THERMODYNAMICS are true.", &
                 units="Pa", default=2.0e7, scale=US%kg_m3_to_R*US%m_s_to_L_T**2)
    call get_param(param_file, mdl, "REGRID_COMPRESSIBILITY_FRACTION", tmpReal, &
                 "When interpolating potential density profiles we can add "//&
                 "some artificial compressibility solely to make homogeneous "//&
                 "regions appear stratified.", units="nondim", default=0.)
    call set_regrid_params(CS, compress_fraction=tmpReal, ref_pressure=P_Ref)
  endif

  if (main_parameters) then
    call get_param(param_file, mdl, "MIN_THICKNESS", tmpReal, &
                 "When regridding, this is the minimum layer "//&
                 "thickness allowed.", units="m", scale=GV%m_to_H, &
                 default=regriddingDefaultMinThickness )
    call set_regrid_params(CS, min_thickness=tmpReal)
  else
    call set_regrid_params(CS, min_thickness=0.)
  endif

  if (coordinateMode(coord_mode) == REGRIDDING_SLIGHT) then
    ! Set SLight-specific regridding parameters.
    call get_param(param_file, mdl, "SLIGHT_DZ_SURFACE", dz_fixed_sfc, &
                 "The nominal thickness of fixed thickness near-surface "//&
                 "layers with the SLight coordinate.", units="m", default=1.0, scale=GV%m_to_H)
    call get_param(param_file, mdl, "SLIGHT_NZ_SURFACE_FIXED", nz_fixed_sfc, &
                 "The number of fixed-depth surface layers with the SLight "//&
                 "coordinate.", units="nondimensional", default=2)
    call get_param(param_file, mdl, "SLIGHT_SURFACE_AVG_DEPTH", Rho_avg_depth, &
                 "The thickness of the surface region over which to average "//&
                 "when calculating the density to use to define the interior "//&
                 "with the SLight coordinate.", units="m", default=1.0, scale=GV%m_to_H)
    call get_param(param_file, mdl, "SLIGHT_NLAY_TO_INTERIOR", nlay_sfc_int, &
                 "The number of layers to offset the surface density when "//&
                 "defining where the interior ocean starts with SLight.", &
                 units="nondimensional", default=2.0)
    call get_param(param_file, mdl, "SLIGHT_FIX_HALOCLINES", fix_haloclines, &
                 "If true, identify regions above the reference pressure "//&
                 "where the reference pressure systematically underestimates "//&
                 "the stratification and use this in the definition of the "//&
                 "interior with the SLight coordinate.", default=.false.)

    call set_regrid_params(CS, dz_min_surface=dz_fixed_sfc, &
                nz_fixed_surface=nz_fixed_sfc, Rho_ML_avg_depth=Rho_avg_depth, &
                nlay_ML_to_interior=nlay_sfc_int, fix_haloclines=fix_haloclines)
    if (fix_haloclines) then
      ! Set additional parameters related to SLIGHT_FIX_HALOCLINES.
      call get_param(param_file, mdl, "HALOCLINE_FILTER_LENGTH", filt_len, &
                 "A length scale over which to smooth the temperature and "//&
                 "salinity before identifying erroneously unstable haloclines.", &
                 units="m", default=2.0, scale=GV%m_to_H)
      call get_param(param_file, mdl, "HALOCLINE_STRAT_TOL", strat_tol, &
                 "A tolerance for the ratio of the stratification of the "//&
                 "apparent coordinate stratification to the actual value "//&
                 "that is used to identify erroneously unstable haloclines. "//&
                 "This ratio is 1 when they are equal, and sensible values "//&
                 "are between 0 and 0.5.", units="nondimensional", default=0.2)
      call set_regrid_params(CS, halocline_filt_len=filt_len, &
                             halocline_strat_tol=strat_tol)
    endif

  endif

  if (coordinateMode(coord_mode) == REGRIDDING_ADAPTIVE) then
    call get_param(param_file, mdl, "ADAPT_TIME_RATIO", adaptTimeRatio, &
                 "Ratio of ALE timestep to grid timescale.", units="nondim", default=1.0e-1)
    call get_param(param_file, mdl, "ADAPT_ZOOM_DEPTH", adaptZoom, &
                 "Depth of near-surface zooming region.", units="m", default=200.0, scale=GV%m_to_H)
    call get_param(param_file, mdl, "ADAPT_ZOOM_COEFF", adaptZoomCoeff, &
                 "Coefficient of near-surface zooming diffusivity.", units="nondim", default=0.2)
    call get_param(param_file, mdl, "ADAPT_BUOY_COEFF", adaptBuoyCoeff, &
                 "Coefficient of buoyancy diffusivity.", units="nondim", default=0.8)
    call get_param(param_file, mdl, "ADAPT_ALPHA", adaptAlpha, &
                 "Scaling on optimization tendency.", units="nondim", default=1.0)
    call get_param(param_file, mdl, "ADAPT_DO_MIN_DEPTH", tmpLogical, &
                 "If true, make a HyCOM-like mixed layer by preventing interfaces "//&
                 "from being shallower than the depths specified by the regridding coordinate.", &
                 default=.false.)
    call get_param(param_file, mdl, "ADAPT_DRHO0", adaptDrho0, &
                 "Reference density difference for stratification-dependent diffusion.", &
                 units="kg m-3", default=0.5, scale=US%kg_m3_to_R)

    call set_regrid_params(CS, adaptTimeRatio=adaptTimeRatio, adaptZoom=adaptZoom, &
         adaptZoomCoeff=adaptZoomCoeff, adaptBuoyCoeff=adaptBuoyCoeff, adaptAlpha=adaptAlpha, &
         adaptDoMin=tmpLogical, adaptDrho0=adaptDrho0)
  endif

  if (main_parameters .and. coord_is_state_dependent) then
    call get_param(param_file, mdl, "MAXIMUM_INT_DEPTH_CONFIG", string, &
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
      call get_param(param_file, mdl, "MAXIMUM_INTERFACE_DEPTHS", z_max, &
                   trim(message), units="m", scale=GV%m_to_H, fail_if_missing=.true.)
      call set_regrid_max_depths(CS, z_max)
    elseif (index(trim(string),'FILE:')==1) then
      if (string(6:6)=='.' .or. string(6:6)=='/') then
        ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
        fileName = trim( extractWord(trim(string(6:80)), 1) )
      else
        ! Otherwise assume we should look for the file in INPUTDIR
        fileName = trim(inputdir) // trim( extractWord(trim(string(6:80)), 1) )
      endif
      if (.not. file_exists(fileName)) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
        "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

      do_sum = .false.
      varName = trim( extractWord(trim(string(6:)), 2) )
      if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
        "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
      if (len_trim(varName)==0) then
        if (field_exists(fileName,'z_max')) then; varName = 'z_max'
        elseif (field_exists(fileName,'dz')) then; varName = 'dz' ; do_sum = .true.
        elseif (field_exists(fileName,'dz_max')) then; varName = 'dz_max' ; do_sum = .true.
        else ; call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
          "MAXIMUM_INT_DEPTHS variable not specified and none could be guessed.")
        endif
      endif
      if (do_sum) then
        call MOM_read_data(trim(fileName), trim(varName), dz_max)
        z_max(1) = 0.0 ; do K=1,ke ; z_max(K+1) = z_max(K) + dz_max(k) ; enddo
      else
        call MOM_read_data(trim(fileName), trim(varName), z_max)
      endif
      call log_param(param_file, mdl, "!MAXIMUM_INT_DEPTHS", z_max, &
                 trim(message), units=coordinateUnits(coord_mode))
      call set_regrid_max_depths(CS, z_max, GV%m_to_H)
    elseif (index(trim(string),'FNC1:')==1) then
      call dz_function1( trim(string(6:)), dz_max )
      if ((coordinateMode(coord_mode) == REGRIDDING_SLIGHT) .and. &
          (dz_fixed_sfc > 0.0)) then
        do k=1,nz_fixed_sfc ; dz_max(k) = dz_fixed_sfc ; enddo
      endif
      z_max(1) = 0.0 ; do K=1,ke ; z_max(K+1) = z_max(K) + dz_max(K) ; enddo
      call log_param(param_file, mdl, "!MAXIMUM_INT_DEPTHS", z_max, &
                 trim(message), units=coordinateUnits(coord_mode))
      call set_regrid_max_depths(CS, z_max, GV%m_to_H)
    else
      call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
        "Unrecognized MAXIMUM_INT_DEPTH_CONFIG "//trim(string))
    endif
    deallocate(z_max)
    deallocate(dz_max)

    ! Optionally specify maximum thicknesses for each layer, enforced by moving
    ! the interface below a layer downward.
    call get_param(param_file, mdl, "MAX_LAYER_THICKNESS_CONFIG", string, &
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
      call get_param(param_file, mdl, "MAX_LAYER_THICKNESS", h_max, &
                   trim(message), units="m", fail_if_missing=.true., scale=GV%m_to_H)
      call set_regrid_max_thickness(CS, h_max)
    elseif (index(trim(string),'FILE:')==1) then
      if (string(6:6)=='.' .or. string(6:6)=='/') then
        ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
        fileName = trim( extractWord(trim(string(6:80)), 1) )
      else
        ! Otherwise assume we should look for the file in INPUTDIR
        fileName = trim(inputdir) // trim( extractWord(trim(string(6:80)), 1) )
      endif
      if (.not. file_exists(fileName)) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
        "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(string)//")")

      varName = trim( extractWord(trim(string(6:)), 2) )
      if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
        "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(string)//")")
      if (len_trim(varName)==0) then
        if (field_exists(fileName,'h_max')) then; varName = 'h_max'
        elseif (field_exists(fileName,'dz_max')) then; varName = 'dz_max'
        else ; call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
          "MAXIMUM_INT_DEPTHS variable not specified and none could be guessed.")
        endif
      endif
      call MOM_read_data(trim(fileName), trim(varName), h_max)
      call log_param(param_file, mdl, "!MAX_LAYER_THICKNESS", h_max, &
                 trim(message), units=coordinateUnits(coord_mode))
      call set_regrid_max_thickness(CS, h_max, GV%m_to_H)
    elseif (index(trim(string),'FNC1:')==1) then
      call dz_function1( trim(string(6:)), h_max )
      call log_param(param_file, mdl, "!MAX_LAYER_THICKNESS", h_max, &
                 trim(message), units=coordinateUnits(coord_mode))
      call set_regrid_max_thickness(CS, h_max, GV%m_to_H)
    else
      call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
        "Unrecognized MAX_LAYER_THICKNESS_CONFIG "//trim(string))
    endif
    deallocate(h_max)
  endif

  if (allocated(dz)) deallocate(dz)
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
  integer :: i

  ierr = .false.
  status = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid)
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
  ! NF90_GET_ATT can return attributes with null characters, which TRIM will not truncate.
  ! This loop replaces any null characters with a space so that the following check between
  ! the read units and the expected units will pass
  do i=1,LEN_TRIM(units)
    if (units(i:i) == CHAR(0)) units(i:i) = " "
  enddo

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

  if (associated(CS%zlike_CS))  call end_coord_zlike(CS%zlike_CS)
  if (associated(CS%sigma_CS))  call end_coord_sigma(CS%sigma_CS)
  if (associated(CS%rho_CS))    call end_coord_rho(CS%rho_CS)
  if (associated(CS%hycom_CS))  call end_coord_hycom(CS%hycom_CS)
  if (associated(CS%slight_CS)) call end_coord_slight(CS%slight_CS)
  if (associated(CS%adapt_CS))  call end_coord_adapt(CS%adapt_CS)

  deallocate( CS%coordinateResolution )
  if (allocated(CS%target_density)) deallocate( CS%target_density )
  if (allocated(CS%max_interface_depths) ) deallocate( CS%max_interface_depths )
  if (allocated(CS%max_layer_thickness) ) deallocate( CS%max_layer_thickness )

end subroutine end_regridding

!------------------------------------------------------------------------------
!> Dispatching regridding routine for orchestrating regridding & remapping
subroutine regridding_main( remapCS, CS, G, GV, h, tv, h_new, dzInterface, frac_shelf_h, conv_adjust)
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
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)), intent(inout) :: h      !< Current 3D grid obtained after
                                                                      !! the last time step
  type(thermo_var_ptrs),                      intent(inout) :: tv     !< Thermodynamical variables (T, S, ...)
  real, dimension(SZI_(G),SZJ_(G), CS%nk),    intent(inout) :: h_new  !< New 3D grid consistent with target coordinate
  real, dimension(SZI_(G),SZJ_(G), CS%nk+1),  intent(inout) :: dzInterface !< The change in position of each interface
  real, dimension(:,:),             optional, pointer       :: frac_shelf_h !< Fractional ice shelf coverage
  logical,                          optional, intent(in   ) :: conv_adjust !< If true, do convective adjustment
  ! Local variables
  real :: trickGnuCompiler
  logical :: use_ice_shelf
  logical :: do_convective_adjustment

  do_convective_adjustment = .true.
  if (present(conv_adjust)) do_convective_adjustment = conv_adjust

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
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_SIGMA_SHELF_ZSTAR)
      call build_zstar_grid( CS, G, GV, h, dzInterface )
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_SIGMA )
      call build_sigma_grid( CS, G, GV, h, dzInterface )
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_RHO )
      if (do_convective_adjustment) call convective_adjustment(G, GV, h, tv)
      call build_rho_grid( G, GV, G%US, h, tv, dzInterface, remapCS, CS )
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_ARBITRARY )
      call build_grid_arbitrary( G, GV, h, dzInterface, trickGnuCompiler, CS )
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_HYCOM1 )
      call build_grid_HyCOM1( G, GV, G%US, h, tv, h_new, dzInterface, CS )

    case ( REGRIDDING_SLIGHT )
      call build_grid_SLight( G, GV, G%US, h, tv, dzInterface, CS )
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_ADAPTIVE )
      call build_grid_adaptive(G, GV, G%US, h, tv, dzInterface, remapCS, CS)
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)

    case default
      call MOM_error(FATAL,'MOM_regridding, regridding_main: '//&
                     'Unknown regridding scheme selected!')

  end select ! type of grid

#ifdef __DO_SAFETY_CHECKS__
  call check_remapping_grid(G, GV, h, dzInterface,'in regridding_main')
#endif

end subroutine regridding_main

!> Calculates h_new from h + delta_k dzInterface
subroutine calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                     intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Old layer thicknesses (arbitrary units)
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),  intent(in)    :: dzInterface !< Change in interface positions (same as h)
  real, dimension(SZI_(G),SZJ_(G),CS%nk),    intent(inout) :: h_new !< New layer thicknesses (same as h)
  ! Local variables
  integer :: i, j, k, nki

  nki = min(CS%nk, GV%ke)

  !$OMP parallel do default(shared)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1
      if (G%mask2dT(i,j)>0.) then
        do k=1,nki
          h_new(i,j,k) = max( 0., h(i,j,k) + ( dzInterface(i,j,k) - dzInterface(i,j,k+1) ) )
        enddo
        if (CS%nk > GV%ke) then
          do k=nki+1, CS%nk
            h_new(i,j,k) = max( 0., dzInterface(i,j,k) - dzInterface(i,j,k+1) )
          enddo
        endif
      else
        h_new(i,j,1:nki) = h(i,j,1:nki)
        if (CS%nk > GV%ke) h_new(i,j,nki+1:CS%nk) = 0.
        ! On land points, why are we keeping the original h rather than setting to zero? -AJA
      endif
    enddo
  enddo

end subroutine calc_h_new_by_dz

!> Check that the total thickness of two grids match
subroutine check_remapping_grid( G, GV, h, dzInterface, msg )
  type(ocean_grid_type),                       intent(in) :: G   !< Grid structure
  type(verticalGrid_type),                     intent(in) :: GV  !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(in) :: dzInterface !< Change in interface positions
                                                                 !! [H ~> m or kg m-2]
  character(len=*),                            intent(in) :: msg !< Message to append to errors
  ! Local variables
  integer :: i, j

  !$OMP parallel do default(shared)
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
    if (G%mask2dT(i,j)>0.) call check_grid_column( GV%ke, GV%Z_to_H*G%bathyT(i,j), h(i,j,:), dzInterface(i,j,:), msg )
  enddo ; enddo

end subroutine check_remapping_grid

!> Check that the total thickness of new and old grids are consistent
subroutine check_grid_column( nk, depth, h, dzInterface, msg )
  integer,               intent(in) :: nk !< Number of cells
  real,                  intent(in) :: depth !< Depth of bottom [Z ~> m] or arbitrary units
  real, dimension(nk),   intent(in) :: h  !< Cell thicknesses [Z ~> m] or arbitrary units
  real, dimension(nk+1), intent(in) :: dzInterface !< Change in interface positions (same units as h)
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
  type(regridding_CS),      intent(in)    :: CS !< Regridding control structure
  integer,                  intent(in)    :: nk !< Number of cells in source grid
  real, dimension(nk+1),    intent(in)    :: z_old !< Old grid position [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(in)    :: z_new !< New grid position [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: dz_g !< Change in interface positions [H ~> m or kg m-2]
  ! Local variables
  real :: sgn  ! The sign convention for downward.
  real :: dz_tgt, zr1, z_old_k
  real :: Aq, Bq, dz0, z0, F0
  real :: zs, zd, dzwt, Idzwt
  real :: wtd, Iwtd
  real :: Int_zs, Int_zd, dInt_zs_zd
! For debugging:
  real, dimension(nk+1) :: z_act
!  real, dimension(nk+1) :: ddz_g_s, ddz_g_d
  logical :: debug = .false.
  integer :: k

  if ((z_old(nk+1) - z_old(1)) * (z_new(CS%nk+1) - z_new(1)) < 0.0) then
    call MOM_error(FATAL, "filtered_grid_motion: z_old and z_new use different sign conventions.")
  elseif ((z_old(nk+1) - z_old(1)) * (z_new(CS%nk+1) - z_new(1)) == 0.0) then
    ! This is a massless column, so do nothing and return.
    do k=1,CS%nk+1 ; dz_g(k) = 0.0 ; enddo ; return
  elseif ((z_old(nk+1) - z_old(1)) + (z_new(CS%nk+1) - z_new(1)) > 0.0) then
    sgn = 1.0
  else
    sgn = -1.0
  endif

  if (debug) then
    do k=2,CS%nk+1
      if (sgn*(z_new(k)-z_new(k-1)) < -5e-16*(abs(z_new(k))+abs(z_new(k-1))) ) &
        call MOM_error(FATAL, "filtered_grid_motion: z_new is tangled.")
    enddo
    do k=2,nk+1
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
  z_old_k = z_old(1)
  do k = 2,CS%nk+1
    if (k<=nk+1) z_old_k = z_old(k) ! This allows for virtual z_old interface at bottom of the model
    ! zr1 is positive and increases with depth, and dz_tgt is positive downward.
    dz_tgt = sgn*(z_new(k) - z_old_k)
    zr1 = sgn*(z_old_k - z_old(1))

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
 !dz_g(CS%nk+1) = 0.0

  if (debug) then
    z_old_k = z_old(1)
    do k=1,CS%nk+1
      if (k<=nk+1) z_old_k = z_old(k) ! This allows for virtual z_old interface at bottom of the model
      z_act(k) = z_old_k + dz_g(k)
    enddo
    do k=2,CS%nk+1
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
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                     intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G), CS%nk+1), intent(inout) :: dzInterface !< The change in interface depth
                                                                 !! [H ~> m or kg m-2].
  real, dimension(:,:),            optional, pointer       :: frac_shelf_h !< Fractional ice shelf coverage [nondim].
  ! Local variables
  real    :: nominalDepth, totalThickness, dh  ! Depths and thicknesses [H ~> m or kg m-2]
  real, dimension(SZK_(GV)+1) :: zOld, zNew    ! Coordinate interface heights [H ~> m or kg m-2]
  integer :: i, j, k, nz
  logical :: ice_shelf

  nz = GV%ke
  ice_shelf = .false.
  if (present(frac_shelf_h)) then
    if (associated(frac_shelf_h)) ice_shelf = .true.
  endif

!$OMP parallel do default(none) shared(G,GV,dzInterface,CS,nz,h,frac_shelf_h,ice_shelf) &
!$OMP                          private(nominalDepth,totalThickness,zNew,dh,zOld)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      if (G%mask2dT(i,j)==0.) then
        dzInterface(i,j,:) = 0.
        cycle
      endif

      ! Local depth (G%bathyT is positive)
      nominalDepth = G%bathyT(i,j)*GV%Z_to_H

      ! Determine water column thickness
      totalThickness = 0.0
      do k = 1,nz
        totalThickness = totalThickness + h(i,j,k)
      enddo

      zOld(nz+1) = - nominalDepth
      do k = nz,1,-1
        zOld(k) = zOld(k+1) + h(i,j,k)
      enddo

      if (ice_shelf) then
        if (frac_shelf_h(i,j) > 0.) then ! under ice shelf
          call build_zstar_column(CS%zlike_CS, nominalDepth, totalThickness, zNew, &
                                z_rigid_top = totalThickness-nominalDepth, &
                                eta_orig=zOld(1), zScale=GV%Z_to_H)
        else
          call build_zstar_column(CS%zlike_CS, nominalDepth, totalThickness, &
                                zNew, zScale=GV%Z_to_H)
        endif
      else
        call build_zstar_column(CS%zlike_CS, nominalDepth, totalThickness, &
                                zNew, zScale=GV%Z_to_H)
      endif

      ! Calculate the final change in grid position after blending new and old grids
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
          write(0,*) k,h(i,j,k),zNew(k)-zNew(k+1),CS%coordinateResolution(k)
        enddo
        call MOM_error( FATAL, &
               'MOM_regridding, build_zstar_grid(): top surface has moved!!!' )
      endif
#endif

      call adjust_interface_motion( CS, nz, h(i,j,:), dzInterface(i,j,:) )

    enddo
  enddo

end subroutine build_zstar_grid

!------------------------------------------------------------------------------
! Build sigma grid
!> This routine builds a grid based on terrain-following coordinates.
subroutine build_sigma_grid( CS, G, GV, h, dzInterface )
!------------------------------------------------------------------------------
! This routine builds a grid based on terrain-following coordinates.
! The module parameter coordinateResolution(:) determines the resolution in
! sigma coordinate, dSigma(:). sigma-coordinates are defined by
!   sigma = (eta-z)/(H+eta)  s.t. sigma=0 at z=eta and sigma=1 at z=-H .
!------------------------------------------------------------------------------

  ! Arguments
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                     intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G), CS%nk+1), intent(inout) :: dzInterface !< The change in interface depth
                                                                 !! [H ~> m or kg m-2]

  ! Local variables
  integer :: i, j, k
  integer :: nz
  real    :: nominalDepth, totalThickness, dh
  real, dimension(SZK_(GV)+1) :: zOld, zNew

  nz = GV%ke

  do i = G%isc-1,G%iec+1
    do j = G%jsc-1,G%jec+1

      if (G%mask2dT(i,j)==0.) then
        dzInterface(i,j,:) = 0.
        cycle
      endif

      ! The rest of the model defines grids integrating up from the bottom
      nominalDepth = G%bathyT(i,j)*GV%Z_to_H

      ! Determine water column height
      totalThickness = 0.0
      do k = 1,nz
        totalThickness = totalThickness + h(i,j,k)
      enddo

      call build_sigma_column(CS%sigma_CS, nominalDepth, totalThickness, zNew)

      ! Calculate the final change in grid position after blending new and old grids
      zOld(nz+1) =  -nominalDepth
      do k = nz,1,-1
        zOld(k) = zOld(k+1) + h(i, j, k)
      enddo

      call filtered_grid_motion( CS, nz, zOld, zNew, dzInterface(i,j,:) )

#ifdef __DO_SAFETY_CHECKS__
      dh=max(nominalDepth,totalThickness)
      if (abs(zNew(1)-zOld(1))>(CS%nk-1)*0.5*epsilon(dh)*dh) then
        write(0,*) 'min_thickness=',CS%min_thickness
        write(0,*) 'nominalDepth=',nominalDepth,'totalThickness=',totalThickness
        write(0,*) 'dzInterface(1) = ',dzInterface(i,j,1),epsilon(dh),nz,CS%nk
        do k=1,nz+1
          write(0,*) k,zOld(k),zNew(k)
        enddo
        do k=1,CS%nk
          write(0,*) k,h(i,j,k),zNew(k)-zNew(k+1),totalThickness*CS%coordinateResolution(k),CS%coordinateResolution(k)
        enddo
        call MOM_error( FATAL, &
               'MOM_regridding, build_sigma_grid: top surface has moved!!!' )
      endif
      dzInterface(i,j,1) = 0.
      dzInterface(i,j,CS%nk+1) = 0.
#endif

    enddo
  enddo

end subroutine build_sigma_grid

!------------------------------------------------------------------------------
! Build grid based on target interface densities
!------------------------------------------------------------------------------
!> This routine builds a new grid based on a given set of target interface densities.
subroutine build_rho_grid( G, GV, US, h, tv, dzInterface, remapCS, CS )
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
  type(unit_scale_type),                        intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                        intent(in)    :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)+1), intent(inout) :: dzInterface !< The change in interface depth
                                                                    !! [H ~> m or kg m-2]
  type(remapping_CS),                           intent(in)    :: remapCS !< The remapping control structure
  type(regridding_CS),                          intent(in)    :: CS !< Regridding control structure

  ! Local variables
  integer :: nz
  integer :: i, j, k
  real    :: nominalDepth   ! Depth of the bottom of the ocean, positive downward [H ~> m or kg m-2]
  real, dimension(SZK_(GV)+1) :: zOld, zNew ! Old and new interface heights [H ~> m or kg m-2]
  real :: h_neglect, h_neglect_edge ! Negligible thicknesses [H ~> m or kg m-2]
#ifdef __DO_SAFETY_CHECKS__
  real    :: totalThickness
  real    :: dh
#endif

  if (.not.CS%remap_answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  nz = GV%ke

  if (.not.CS%target_density_set) call MOM_error(FATAL, "build_rho_grid: "//&
        "Target densities must be set before build_rho_grid is called.")

  ! Build grid based on target interface densities
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      if (G%mask2dT(i,j)==0.) then
        dzInterface(i,j,:) = 0.
        cycle
      endif


      ! Local depth (G%bathyT is positive)
      nominalDepth = G%bathyT(i,j)*GV%Z_to_H

      call build_rho_column(CS%rho_CS, nz, nominalDepth, h(i, j, :), &
                            tv%T(i, j, :), tv%S(i, j, :), tv%eqn_of_state, zNew, &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)

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

    enddo  ! end loop on i
  enddo  ! end loop on j

end subroutine build_rho_grid

!> Builds a simple HyCOM-like grid with the deepest location of potential
!! density interpolated from the column profile and a clipping of depth for
!! each interface to a fixed z* or p* grid.  This should probably be (optionally?)
!! changed to find the nearest location of the target density.
!! \remark { Based on Bleck, 2002: An oceanice general circulation model framed in
!! hybrid isopycnic-Cartesian coordinates, Ocean Modelling 37, 55-88.
!! http://dx.doi.org/10.1016/S1463-5003(01)00012-9 }
subroutine build_grid_HyCOM1( G, GV, US, h, tv, h_new, dzInterface, CS )
  type(ocean_grid_type),                     intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Existing model thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Thermodynamics structure
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk),    intent(inout) :: h_new !< New layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),  intent(inout) :: dzInterface !< Changes in interface position

  ! Local variables
  real, dimension(SZK_(GV)+1) :: z_col ! Source interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(CS%nk+1) :: z_col_new ! New interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(SZK_(GV)+1) :: dz_col  ! The realized change in z_col [H ~> m or kg m-2]
  real, dimension(SZK_(GV))   :: p_col   ! Layer center pressure [R L2 T-2 ~> Pa]
  real :: ref_pres  ! The reference pressure [R L2 T-2 ~> Pa]
  integer   :: i, j, k, nki
  real :: depth
  real :: h_neglect, h_neglect_edge

  if (.not.CS%remap_answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  if (.not.CS%target_density_set) call MOM_error(FATAL, "build_grid_HyCOM1 : "//&
        "Target densities must be set before build_grid_HyCOM1 is called.")

  nki = min(GV%ke, CS%nk)

  ! Build grid based on target interface densities
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
    if (G%mask2dT(i,j)>0.) then

      depth = G%bathyT(i,j) * GV%Z_to_H

      z_col(1) = 0. ! Work downward rather than bottom up
      do K = 1, GV%ke
        z_col(K+1) = z_col(K) + h(i,j,k)
        p_col(k) = tv%P_Ref + CS%compressibility_fraction * &
             ( 0.5 * ( z_col(K) + z_col(K+1) ) * (GV%H_to_RZ*GV%g_Earth) - tv%P_Ref )
      enddo

      call build_hycom1_column(CS%hycom_CS, tv%eqn_of_state, GV%ke, depth, &
                               h(i,j,:), tv%T(i,j,:), tv%S(i,j,:), p_col, &
                               z_col, z_col_new, zScale=GV%Z_to_H, &
                               h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)

      ! Calculate the final change in grid position after blending new and old grids
      call filtered_grid_motion( CS, GV%ke, z_col, z_col_new, dz_col )

      ! This adjusts things robust to round-off errors
      dz_col(:) = -dz_col(:)
      call adjust_interface_motion( CS, GV%ke, h(i,j,:), dz_col(:) )

      dzInterface(i,j,1:nki+1) = dz_col(1:nki+1)
      if (nki<CS%nk) dzInterface(i,j,nki+2:CS%nk+1) = 0.

    else ! on land
      dzInterface(i,j,:) = 0.
    endif ! mask2dT
  enddo ; enddo ! i,j

  call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)

end subroutine build_grid_HyCOM1

!> This subroutine builds an adaptive grid that follows density surfaces where
!! possible, subject to constraints on the smoothness of interface heights.
subroutine build_grid_adaptive(G, GV, US, h, tv, dzInterface, remapCS, CS)
  type(ocean_grid_type),                       intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                       intent(in)    :: tv   !< A structure pointing to various
                                                                     !! thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: dzInterface !< The change in interface depth
                                                                     !! [H ~> m or kg m-2]
  type(remapping_CS),                          intent(in)    :: remapCS !< The remapping control structure
  type(regridding_CS),                         intent(in)    :: CS   !< Regridding control structure

  ! local variables
  integer :: i, j, k, nz ! indices and dimension lengths
  ! temperature, salinity and pressure on interfaces
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: tInt, sInt
  ! current interface positions and after tendency term is applied
  ! positive downward
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: zInt ! Interface depths [H ~> m or kg m-2]
  real, dimension(SZK_(GV)+1) :: zNext  ! New interface depths [H ~> m or kg m-2]

  nz = GV%ke

  ! position surface at z = 0.
  zInt(:,:,1) = 0.

  ! work on interior interfaces
  do K = 2, nz ; do j = G%jsc-2,G%jec+2 ; do i = G%isc-2,G%iec+2
    tInt(i,j,K) = 0.5 * (tv%T(i,j,k-1) + tv%T(i,j,k))
    sInt(i,j,K) = 0.5 * (tv%S(i,j,k-1) + tv%S(i,j,k))
    zInt(i,j,K) = zInt(i,j,K-1) + h(i,j,k-1) ! zInt in [H]
  enddo ; enddo ; enddo

  ! top and bottom temp/salt interfaces are just the layer
  ! average values
  tInt(:,:,1) = tv%T(:,:,1) ; tInt(:,:,nz+1) = tv%T(:,:,nz)
  sInt(:,:,1) = tv%S(:,:,1) ; sInt(:,:,nz+1) = tv%S(:,:,nz)

  ! set the bottom interface depth
  zInt(:,:,nz+1)  = zInt(:,:,nz) + h(:,:,nz)

  ! calculate horizontal density derivatives (alpha/beta)
  ! between cells in a 5-point stencil, columnwise
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
    if (G%mask2dT(i,j) < 0.5) then
      dzInterface(i,j,:) = 0. ! land point, don't move interfaces, and skip
      cycle
    endif

    call build_adapt_column(CS%adapt_CS, G, GV, US, tv, i, j, zInt, tInt, sInt, h, zNext)

    call filtered_grid_motion(CS, nz, zInt(i,j,:), zNext, dzInterface(i,j,:))
    ! convert from depth to z
    do K = 1, nz+1 ; dzInterface(i,j,K) = -dzInterface(i,j,K) ; enddo
    call adjust_interface_motion(CS, nz, h(i,j,:), dzInterface(i,j,:))
  enddo ; enddo
end subroutine build_grid_adaptive

!> Builds a grid that tracks density interfaces for water that is denser than
!! the surface density plus an increment of some number of layers, and uses all
!! lighter layers uniformly above this location.  Note that this amounts to
!! interpolating to find the depth of an arbitrary (non-integer) interface index
!! which should make the results vary smoothly in space to the extent that the
!! surface density and interior stratification vary smoothly in space.  Over
!! shallow topography, this will tend to give a uniform sigma-like coordinate.
!! For sufficiently shallow water, a minimum grid spacing is used to avoid
!! certain instabilities.
subroutine build_grid_SLight(G, GV, US, h, tv, dzInterface, CS)
  type(ocean_grid_type),                       intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),                     intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h  !< Existing model thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                       intent(in)    :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), intent(inout) :: dzInterface !< Changes in interface position
  type(regridding_CS),                         intent(in)    :: CS !< Regridding control structure

  real, dimension(SZK_(GV)+1) :: z_col   ! Interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(SZK_(GV)+1) :: z_col_new ! Interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(SZK_(GV)+1) :: dz_col  ! The realized change in z_col [H ~> m or kg m-2]
  real, dimension(SZK_(GV))   :: p_col   ! Layer center pressure [R L2 T-2 ~> Pa]
  real :: depth
  integer :: i, j, k, nz
  real :: h_neglect, h_neglect_edge

  if (.not.CS%remap_answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  nz = GV%ke

  if (.not.CS%target_density_set) call MOM_error(FATAL, "build_grid_SLight : "//&
        "Target densities must be set before build_grid_SLight is called.")

  ! Build grid based on target interface densities
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
    if (G%mask2dT(i,j)>0.) then

      depth = G%bathyT(i,j) * GV%Z_to_H
      z_col(1) = 0. ! Work downward rather than bottom up
      do K=1,nz
        z_col(K+1) = z_col(K) + h(i,j,k)
        p_col(k) = tv%P_Ref + CS%compressibility_fraction * &
                    ( 0.5 * ( z_col(K) + z_col(K+1) ) * (GV%H_to_RZ*GV%g_Earth) - tv%P_Ref )
      enddo

      call build_slight_column(CS%slight_CS, tv%eqn_of_state, GV%H_to_RZ*GV%g_Earth, &
                          GV%H_subroundoff, nz, depth, h(i, j, :), &
                          tv%T(i, j, :), tv%S(i, j, :), p_col, z_col, z_col_new, &
                          h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)

      ! Calculate the final change in grid position after blending new and old grids
      call filtered_grid_motion( CS, nz, z_col, z_col_new, dz_col )
      do K=1,nz+1 ; dzInterface(i,j,K) = -dz_col(K) ; enddo
#ifdef __DO_SAFETY_CHECKS__
      if (dzInterface(i,j,1) /= 0.) stop 'build_grid_SLight: Surface moved?!'
      if (dzInterface(i,j,nz+1) /= 0.) stop 'build_grid_SLight: Bottom moved?!'
#endif

      ! This adjusts things robust to round-off errors
      call adjust_interface_motion( CS, nz, h(i,j,:), dzInterface(i,j,:) )

    else ! on land
      dzInterface(i,j,:) = 0.
    endif ! mask2dT
  enddo ; enddo ! i,j

end subroutine build_grid_SLight

!> Adjust dz_Interface to ensure non-negative future thicknesses
subroutine adjust_interface_motion( CS, nk, h_old, dz_int )
  type(regridding_CS),      intent(in)    :: CS !< Regridding control structure
  integer,                  intent(in)    :: nk !< Number of layers in h_old
  real, dimension(nk),      intent(in)    :: h_old !< Minium allowed thickness of h [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: dz_int !< Minium allowed thickness of h [H ~> m or kg m-2]
  ! Local variables
  integer :: k
  real :: h_new, eps, h_total, h_err

  eps = 1. ; eps = epsilon(eps)

  h_total = 0. ; h_err = 0.
  do k = 1, min(CS%nk,nk)
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
  if (CS%nk>nk) then
    do k = nk+1, CS%nk
      h_err = h_err + max( abs(dz_int(k)), abs(dz_int(k+1)) )*eps
      h_new = ( dz_int(k) - dz_int(k+1) )
      if (h_new < -3.0*h_err) then
        write(0,*) 'h<0 at k=',k,'h_old was empty',&
          'wup=',dz_int(k),'wdn=',dz_int(k+1),'dw_dz=',dz_int(k) - dz_int(k+1), &
          'h_new=',h_new,'h_err=',h_err
        call MOM_error( FATAL, 'MOM_regridding: adjust_interface_motion() - '//&
                       'implied h<0 is larger than roundoff!')
      endif
    enddo
  endif
  do k = min(CS%nk,nk),2,-1
    h_new = h_old(k) + ( dz_int(k) - dz_int(k+1) )
    if (h_new<CS%min_thickness) &
      dz_int(k) = ( dz_int(k+1) - h_old(k) ) + CS%min_thickness ! Implies next h_new = min_thickness
    h_new = h_old(k) + ( dz_int(k) - dz_int(k+1) )
    if (h_new<0.) &
      dz_int(k) = ( 1. - eps ) * ( dz_int(k+1) - h_old(k) ) ! Backup in case min_thickness==0
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
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)),   intent(in)    :: h  !< Original layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)+1), intent(inout) :: dzInterface !< The change in interface
                                                                    !! depth [H ~> m or kg m-2]
  real,                                         intent(inout) :: h_new !< New layer thicknesses [H ~> m or kg m-2]
  type(regridding_CS),                          intent(in)    :: CS !< Regridding control structure

  ! Local variables
  integer   :: i, j, k
  integer   :: nz
  real      :: z_inter(SZK_(GV)+1)
  real      :: total_height
  real      :: delta_h
  real      :: max_depth
  real      :: eta              ! local elevation
  real      :: local_depth
  real      :: x1, y1, x2, y2
  real      :: x, t

  nz = GV%ke
  max_depth = G%max_depth*GV%Z_to_H

  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      ! Local depth
      local_depth = G%bathyT(i,j)*GV%Z_to_H

      ! Determine water column height
      total_height = 0.0
      do k = 1,nz
        total_height = total_height + h(i,j,k)
      enddo

      eta = total_height - local_depth

      ! Compute new thicknesses based on stretched water column
      delta_h = (max_depth + eta) / nz

      ! Define interfaces
      z_inter(1) = eta
      do k = 1,nz
        z_inter(k+1) = z_inter(k) - delta_h
      enddo

      ! Refine grid in the middle
      do k = 1,nz+1
        x1 = 0.35; y1 = 0.45; x2 = 0.65; y2 = 0.55

        x = - ( z_inter(k) - eta ) / max_depth

        if ( x <= x1 ) then
          t = y1*x/x1
        elseif ( (x > x1 ) .and. ( x < x2 )) then
          t = y1 + (y2-y1) * (x-x1) / (x2-x1)
        else
          t = y2 + (1.0-y2) * (x-x2) / (1.0-x2)
        endif

        z_inter(k) = -t * max_depth + eta

      enddo

      ! Modify interface heights to account for topography
      z_inter(nz+1) = - local_depth

      ! Modify interface heights to avoid layers of zero thicknesses
      do k = nz,1,-1
        if ( z_inter(k) < (z_inter(k+1) + CS%min_thickness) ) then
          z_inter(k) = z_inter(k+1) + CS%min_thickness
        endif
      enddo

      ! Chnage in interface position
      x = 0. ! Left boundary at x=0
      dzInterface(i,j,1) = 0.
      do k = 2,nz
        x = x + h(i,j,k)
        dzInterface(i,j,k) = z_inter(k) - x
      enddo
      dzInterface(i,j,nz+1) = 0.

    enddo
  enddo

stop 'OOOOOOPS' ! For some reason the gnu compiler will not let me delete this
                ! routine????

end subroutine build_grid_arbitrary



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
  type(regridding_CS),                    intent(in)    :: CS   !< Regridding control structure
  type(ocean_grid_type),                  intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(GV)), intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2]

  ! Local variables
  integer :: i, j, k
  real    :: hTmp(GV%ke)

  do i = G%isc-1,G%iec+1
    do j = G%jsc-1,G%jec+1

      ! Build grid for current column
      do k = 1,GV%ke
        hTmp(k) = h(i,j,k)
      enddo

      call old_inflate_layers_1d( CS%min_thickness, GV%ke, hTmp )

      ! Save modified grid
      do k = 1,GV%ke
        h(i,j,k) = hTmp(k)
      enddo

    enddo
  enddo

end subroutine inflate_vanished_layers_old

!------------------------------------------------------------------------------
!> Achieve convective adjustment by swapping layers
subroutine convective_adjustment(G, GV, h, tv)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(inout) :: tv   !< A structure pointing to various thermodynamic variables
!------------------------------------------------------------------------------
! Check each water column to see whether it is stratified. If not, sort the
! layers by successive swappings of water masses (bubble sort algorithm)
!------------------------------------------------------------------------------

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
    call calculate_density( tv%T(i,j,:), tv%S(i,j,:), p_col, densities, tv%eqn_of_state)

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
          call calculate_density( tv%T(i,j,k), tv%S(i,j,k), p_col(k), densities(k), tv%eqn_of_state)
          call calculate_density( tv%T(i,j,k+1), tv%S(i,j,k+1), p_col(k+1), &
                                  densities(k+1), tv%eqn_of_state )
          stratified = .false.
        endif
      enddo  ! k

      if ( stratified ) exit
    enddo

  enddo ; enddo  ! i & j

end subroutine convective_adjustment


!------------------------------------------------------------------------------
!> Return a uniform resolution vector in the units of the coordinate
function uniformResolution(nk,coordMode,maxDepth,rhoLight,rhoHeavy)
!------------------------------------------------------------------------------
! Calculate a vector of uniform resolution in the units of the coordinate
!------------------------------------------------------------------------------
  ! Arguments
  integer,          intent(in) :: nk !< Number of cells in source grid
  character(len=*), intent(in) :: coordMode !< A string indicating the coordinate mode.
                                            !! See the documenttion for regrid_consts
                                            !! for the recognized values.
  real,             intent(in) :: maxDepth  !< The range of the grid values in some modes
  real,             intent(in) :: rhoLight  !< The minimum value of the grid in RHO mode
  real,             intent(in) :: rhoHeavy  !< The maximum value of the grid in RHO mode

  real                         :: uniformResolution(nk) !< The returned uniform resolution grid.

  ! Local variables
  integer :: scheme

  scheme = coordinateMode(coordMode)
  select case ( scheme )

    case ( REGRIDDING_ZSTAR, REGRIDDING_HYCOM1, REGRIDDING_SLIGHT, REGRIDDING_SIGMA_SHELF_ZSTAR, &
           REGRIDDING_ADAPTIVE )
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

!> Initialize the coordinate resolutions by calling the appropriate initialization
!! routine for the specified coordinate mode.
subroutine initCoord(CS, GV, US, coord_mode)
  type(regridding_CS),     intent(inout) :: CS !< Regridding control structure
  character(len=*),        intent(in)    :: coord_mode !< A string indicating the coordinate mode.
                                               !! See the documentation for regrid_consts
                                               !! for the recognized values.
  type(verticalGrid_type), intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type

  select case (coordinateMode(coord_mode))
  case (REGRIDDING_ZSTAR)
    call init_coord_zlike(CS%zlike_CS, CS%nk, CS%coordinateResolution)
  case (REGRIDDING_SIGMA_SHELF_ZSTAR)
    call init_coord_zlike(CS%zlike_CS, CS%nk, CS%coordinateResolution)
  case (REGRIDDING_SIGMA)
    call init_coord_sigma(CS%sigma_CS, CS%nk, CS%coordinateResolution)
  case (REGRIDDING_RHO)
    call init_coord_rho(CS%rho_CS, CS%nk, CS%ref_pressure, CS%target_density, CS%interp_CS)
  case (REGRIDDING_HYCOM1)
    call init_coord_hycom(CS%hycom_CS, CS%nk, CS%coordinateResolution, CS%target_density, &
                          CS%interp_CS)
  case (REGRIDDING_SLIGHT)
    call init_coord_slight(CS%slight_CS, CS%nk, CS%ref_pressure, CS%target_density, &
                           CS%interp_CS, GV%m_to_H)
  case (REGRIDDING_ADAPTIVE)
    call init_coord_adapt(CS%adapt_CS, CS%nk, CS%coordinateResolution, GV%m_to_H, US%kg_m3_to_R)
  end select
end subroutine initCoord

!------------------------------------------------------------------------------
!> Set the fixed resolution data
subroutine setCoordinateResolution( dz, CS, scale )
  real, dimension(:),  intent(in)    :: dz !< A vector of vertical grid spacings
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure
  real,      optional, intent(in)    :: scale !< A scaling factor converting dz to coordRes

  if (size(dz)/=CS%nk) call MOM_error( FATAL, &
      'setCoordinateResolution: inconsistent number of levels' )

  if (present(scale)) then
    CS%coordinateResolution(:) = scale*dz(:)
  else
    CS%coordinateResolution(:) = dz(:)
  endif

end subroutine setCoordinateResolution

!> Set target densities based on the old Rlay variable
subroutine set_target_densities_from_GV( GV, US, CS )
  type(verticalGrid_type), intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(regridding_CS),     intent(inout) :: CS !< Regridding control structure
  ! Local variables
  integer :: k, nz

  nz = CS%nk
  CS%target_density(1)    = (GV%Rlay(1) + 0.5*(GV%Rlay(1)-GV%Rlay(2)))
  CS%target_density(nz+1) = (GV%Rlay(nz) + 0.5*(GV%Rlay(nz)-GV%Rlay(nz-1)))
  do k = 2,nz
    CS%target_density(k) = CS%target_density(k-1) + CS%coordinateResolution(k)
  enddo
  CS%target_density_set = .true.

end subroutine set_target_densities_from_GV

!> Set target densities based on vector of interface values
subroutine set_target_densities( CS, rho_int )
  type(regridding_CS),      intent(inout) :: CS !< Regridding control structure
  real, dimension(CS%nk+1), intent(in)    :: rho_int !< Interface densities [R ~> kg m-3]

  if (size(CS%target_density)/=size(rho_int)) then
    call MOM_error(FATAL, "set_target_densities inconsistent args!")
  endif

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

  val_to_H = 1.0 ; if (present(units_to_H)) val_to_H = units_to_H
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

  ! set max depths for coordinate
  select case (CS%regridding_scheme)
  case (REGRIDDING_HYCOM1)
    call set_hycom_params(CS%hycom_CS, max_interface_depths=CS%max_interface_depths)
  case (REGRIDDING_SLIGHT)
    call set_slight_params(CS%slight_CS, max_interface_depths=CS%max_interface_depths)
  end select
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

  ! set max thickness for coordinate
  select case (CS%regridding_scheme)
  case (REGRIDDING_HYCOM1)
    call set_hycom_params(CS%hycom_CS, max_layer_thickness=CS%max_layer_thickness)
  case (REGRIDDING_SLIGHT)
    call set_slight_params(CS%slight_CS, max_layer_thickness=CS%max_layer_thickness)
  end select
end subroutine set_regrid_max_thickness


!------------------------------------------------------------------------------
!> Query the fixed resolution data
function getCoordinateResolution( CS, undo_scaling )
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  logical,   optional, intent(in) :: undo_scaling !< If present and true, undo any internal
                                        !! rescaling of the resolution data.
  real, dimension(CS%nk)          :: getCoordinateResolution

  logical :: unscale
  unscale = .false. ; if (present(undo_scaling)) unscale = undo_scaling

  if (unscale) then
    getCoordinateResolution(:) = CS%coord_scale * CS%coordinateResolution(:)
  else
    getCoordinateResolution(:) = CS%coordinateResolution(:)
  endif

end function getCoordinateResolution

!> Query the target coordinate interface positions
function getCoordinateInterfaces( CS, undo_scaling )
  type(regridding_CS), intent(in) :: CS                      !< Regridding control structure
  logical,   optional, intent(in) :: undo_scaling            !< If present and true, undo any internal
                                                             !! rescaling of the resolution data.
  real, dimension(CS%nk+1)        :: getCoordinateInterfaces !< Interface positions in target coordinate

  integer :: k
  logical :: unscale
  unscale = .false. ; if (present(undo_scaling)) unscale = undo_scaling

  ! When using a coordinate with target densities, we need to get the actual
  ! densities, rather than computing the interfaces based on resolution
  if (CS%regridding_scheme == REGRIDDING_RHO) then
    if (.not. CS%target_density_set) &
      call MOM_error(FATAL, 'MOM_regridding, getCoordinateInterfaces: '//&
                            'target densities not set!')

    if (unscale) then
      getCoordinateInterfaces(:) = CS%coord_scale * CS%target_density(:)
    else
      getCoordinateInterfaces(:) = CS%target_density(:)
    endif
  else
    if (unscale) then
      getCoordinateInterfaces(1) = 0.
      do k = 1, CS%nk
        getCoordinateInterfaces(K+1) = getCoordinateInterfaces(K) - &
                                       CS%coord_scale * CS%coordinateResolution(k)
      enddo
    else
      getCoordinateInterfaces(1) = 0.
      do k = 1, CS%nk
        getCoordinateInterfaces(K+1) = getCoordinateInterfaces(K) - &
                                       CS%coordinateResolution(k)
      enddo
    endif
    ! The following line has an "abs()" to allow ferret users to reference
    ! data by index. It is a temporary work around...  :(  -AJA
    getCoordinateInterfaces(:) = abs( getCoordinateInterfaces(:) )
  endif

end function getCoordinateInterfaces

!------------------------------------------------------------------------------
!> Query the target coordinate units
function getCoordinateUnits( CS )
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  character(len=20)               :: getCoordinateUnits

  select case ( CS%regridding_scheme )
    case ( REGRIDDING_ZSTAR, REGRIDDING_HYCOM1, REGRIDDING_SLIGHT, REGRIDDING_ADAPTIVE )
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
!> Query the short name of the coordinate
function getCoordinateShortName( CS )
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
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
    case ( REGRIDDING_ADAPTIVE )
      getCoordinateShortName = 'adaptive'
    case default
      call MOM_error(FATAL,'MOM_regridding, getCoordinateShortName: '//&
                     'Unknown regridding scheme selected!')
  end select ! type of grid

end function getCoordinateShortName

!> Can be used to set any of the parameters for MOM_regridding.
subroutine set_regrid_params( CS, boundary_extrapolation, min_thickness, old_grid_weight, &
             interp_scheme, depth_of_time_filter_shallow, depth_of_time_filter_deep, &
             compress_fraction, ref_pressure, dz_min_surface, nz_fixed_surface, Rho_ML_avg_depth, &
             nlay_ML_to_interior, fix_haloclines, halocline_filt_len, &
             halocline_strat_tol, integrate_downward_for_e, remap_answers_2018, &
             adaptTimeRatio, adaptZoom, adaptZoomCoeff, adaptBuoyCoeff, adaptAlpha, adaptDoMin, adaptDrho0)
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure
  logical, optional, intent(in) :: boundary_extrapolation !< Extrapolate in boundary cells
  real,    optional, intent(in) :: min_thickness    !< Minimum thickness allowed when building the
                                                    !! new grid [H ~> m or kg m-2]
  real,    optional, intent(in) :: old_grid_weight  !< Weight given to old coordinate when time-filtering grid
  character(len=*), optional, intent(in) :: interp_scheme !< Interpolation method for state-dependent coordinates
  real,    optional, intent(in) :: depth_of_time_filter_shallow !< Depth to start cubic [H ~> m or kg m-2]
  real,    optional, intent(in) :: depth_of_time_filter_deep !< Depth to end cubic [H ~> m or kg m-2]
  real,    optional, intent(in) :: compress_fraction !< Fraction of compressibility to add to potential density [nondim]
  real,    optional, intent(in) :: ref_pressure     !< The reference pressure for density-dependent
                                                    !! coordinates [R L2 T-2 ~> Pa]
  real,    optional, intent(in) :: dz_min_surface   !< The fixed resolution in the topmost
                                                    !! SLight_nkml_min layers [H ~> m or kg m-2]
  integer, optional, intent(in) :: nz_fixed_surface !< The number of fixed-thickness layers at the top of the model
  real,    optional, intent(in) :: Rho_ml_avg_depth !< Averaging depth over which to determine mixed layer potential
                                                    !! density [H ~> m or kg m-2]
  real,    optional, intent(in) :: nlay_ML_to_interior !< Number of layers to offset the mixed layer density to find
                                                    !! resolved stratification [nondim]
  logical, optional, intent(in) :: fix_haloclines   !< Detect regions with much weaker stratification in the coordinate
  real,    optional, intent(in) :: halocline_filt_len !< Length scale over which to filter T & S when looking for
                                                    !! spuriously unstable water mass profiles [H ~> m or kg m-2]
  real,    optional, intent(in) :: halocline_strat_tol !< Value of the stratification ratio that defines a problematic
                                                    !! halocline region.
  logical, optional, intent(in) :: integrate_downward_for_e !< If true, integrate for interface positions downward
                                                    !! from the top.
  logical, optional, intent(in) :: remap_answers_2018 !< If true, use the order of arithmetic and expressions
                                                    !! that recover the remapping answers from 2018.  Otherwise
                                                    !! use more robust but mathematically equivalent expressions.
  real,    optional, intent(in) :: adaptTimeRatio   !< Ratio of the ALE timestep to the grid timescale [nondim].
  real,    optional, intent(in) :: adaptZoom        !< Depth of near-surface zooming region [H ~> m or kg m-2].
  real,    optional, intent(in) :: adaptZoomCoeff   !< Coefficient of near-surface zooming diffusivity [nondim].
  real,    optional, intent(in) :: adaptBuoyCoeff   !< Coefficient of buoyancy diffusivity [nondim].
  real,    optional, intent(in) :: adaptAlpha       !< Scaling factor on optimization tendency [nondim].
  logical, optional, intent(in) :: adaptDoMin       !< If true, make a HyCOM-like mixed layer by
                                                    !! preventing interfaces from being shallower than
                                                    !! the depths specified by the regridding coordinate.
  real,    optional, intent(in) :: adaptDrho0       !< Reference density difference for stratification-dependent
                                                    !! diffusion. [R ~> kg m-3]

  if (present(interp_scheme)) call set_interp_scheme(CS%interp_CS, interp_scheme)
  if (present(boundary_extrapolation)) call set_interp_extrap(CS%interp_CS, boundary_extrapolation)

  if (present(old_grid_weight)) then
    if (old_grid_weight<0. .or. old_grid_weight>1.) &
      call MOM_error(FATAL,'MOM_regridding, set_regrid_params: Weight is out side the range 0..1!')
    CS%old_grid_weight = old_grid_weight
  endif
  if (present(depth_of_time_filter_shallow)) CS%depth_of_time_filter_shallow = depth_of_time_filter_shallow
  if (present(depth_of_time_filter_deep)) CS%depth_of_time_filter_deep = depth_of_time_filter_deep
  if (present(depth_of_time_filter_shallow) .or. present(depth_of_time_filter_deep)) then
    if (CS%depth_of_time_filter_deep<CS%depth_of_time_filter_shallow) call MOM_error(FATAL,'MOM_regridding, '//&
                     'set_regrid_params: depth_of_time_filter_deep<depth_of_time_filter_shallow!')
  endif

  if (present(min_thickness)) CS%min_thickness = min_thickness
  if (present(compress_fraction)) CS%compressibility_fraction = compress_fraction
  if (present(ref_pressure)) CS%ref_pressure = ref_pressure
  if (present(integrate_downward_for_e)) CS%integrate_downward_for_e = integrate_downward_for_e
  if (present(remap_answers_2018)) CS%remap_answers_2018 = remap_answers_2018

  select case (CS%regridding_scheme)
  case (REGRIDDING_ZSTAR)
    if (present(min_thickness)) call set_zlike_params(CS%zlike_CS, min_thickness=min_thickness)
  case (REGRIDDING_SIGMA_SHELF_ZSTAR)
    if (present(min_thickness)) call set_zlike_params(CS%zlike_CS, min_thickness=min_thickness)
  case (REGRIDDING_SIGMA)
    if (present(min_thickness)) call set_sigma_params(CS%sigma_CS, min_thickness=min_thickness)
  case (REGRIDDING_RHO)
    if (present(min_thickness)) call set_rho_params(CS%rho_CS, min_thickness=min_thickness)
    if (present(integrate_downward_for_e)) &
      call set_rho_params(CS%rho_CS, integrate_downward_for_e=integrate_downward_for_e)
    if (associated(CS%rho_CS) .and. (present(interp_scheme) .or. present(boundary_extrapolation))) &
      call set_rho_params(CS%rho_CS, interp_CS=CS%interp_CS)
  case (REGRIDDING_HYCOM1)
    if (associated(CS%hycom_CS) .and. (present(interp_scheme) .or. present(boundary_extrapolation))) &
      call set_hycom_params(CS%hycom_CS, interp_CS=CS%interp_CS)
  case (REGRIDDING_SLIGHT)
    if (present(min_thickness))       call set_slight_params(CS%slight_CS, min_thickness=min_thickness)
    if (present(dz_min_surface))      call set_slight_params(CS%slight_CS, dz_ml_min=dz_min_surface)
    if (present(nz_fixed_surface))    call set_slight_params(CS%slight_CS, nz_fixed_surface=nz_fixed_surface)
    if (present(Rho_ML_avg_depth))    call set_slight_params(CS%slight_CS, Rho_ML_avg_depth=Rho_ML_avg_depth)
    if (present(nlay_ML_to_interior)) call set_slight_params(CS%slight_CS, nlay_ML_offset=nlay_ML_to_interior)
    if (present(fix_haloclines))      call set_slight_params(CS%slight_CS, fix_haloclines=fix_haloclines)
    if (present(halocline_filt_len))  call set_slight_params(CS%slight_CS, halocline_filter_length=halocline_filt_len)
    if (present(halocline_strat_tol)) call set_slight_params(CS%slight_CS, halocline_strat_tol=halocline_strat_tol)
    if (present(compress_fraction))   call set_slight_params(CS%slight_CS, compressibility_fraction=compress_fraction)
    if (associated(CS%slight_CS) .and. (present(interp_scheme) .or. present(boundary_extrapolation))) &
      call set_slight_params(CS%slight_CS, interp_CS=CS%interp_CS)
  case (REGRIDDING_ADAPTIVE)
    if (present(adaptTimeRatio)) call set_adapt_params(CS%adapt_CS, adaptTimeRatio=adaptTimeRatio)
    if (present(adaptZoom))      call set_adapt_params(CS%adapt_CS, adaptZoom=adaptZoom)
    if (present(adaptZoomCoeff)) call set_adapt_params(CS%adapt_CS, adaptZoomCoeff=adaptZoomCoeff)
    if (present(adaptBuoyCoeff)) call set_adapt_params(CS%adapt_CS, adaptBuoyCoeff=adaptBuoyCoeff)
    if (present(adaptAlpha))     call set_adapt_params(CS%adapt_CS, adaptAlpha=adaptAlpha)
    if (present(adaptDoMin))     call set_adapt_params(CS%adapt_CS, adaptDoMin=adaptDoMin)
    if (present(adaptDrho0))     call set_adapt_params(CS%adapt_CS, adaptDrho0=adaptDrho0)
  end select

end subroutine set_regrid_params

!> Returns the number of levels/layers in the regridding control structure
integer function get_regrid_size(CS)
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure

  get_regrid_size = CS%nk

end function get_regrid_size

!> This returns a copy of the zlike_CS stored in the regridding control structure.
function get_zlike_CS(CS)
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  type(zlike_CS) :: get_zlike_CS

  get_zlike_CS = CS%zlike_CS
end function get_zlike_CS

!> This returns a copy of the sigma_CS stored in the regridding control structure.
function get_sigma_CS(CS)
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  type(sigma_CS) :: get_sigma_CS

  get_sigma_CS = CS%sigma_CS
end function get_sigma_CS

!> This returns a copy of the rho_CS stored in the regridding control structure.
function get_rho_CS(CS)
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  type(rho_CS) :: get_rho_CS

  get_rho_CS = CS%rho_CS
end function get_rho_CS

!------------------------------------------------------------------------------
!> Return coordinate-derived thicknesses for fixed coordinate systems
function getStaticThickness( CS, SSH, depth )
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  real,                intent(in) :: SSH   !< The sea surface height, in the same units as depth
  real,                intent(in) :: depth !< The maximum depth of the grid, often [Z ~> m]
  real, dimension(CS%nk)          :: getStaticThickness !< The returned thicknesses in the units of depth
  ! Local
  integer :: k
  real :: z, dz

  select case ( CS%regridding_scheme )
    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA_SHELF_ZSTAR, REGRIDDING_HYCOM1, REGRIDDING_SLIGHT, REGRIDDING_ADAPTIVE )
      if (depth>0.) then
        z = ssh
        do k = 1, CS%nk
          dz = CS%coordinateResolution(k) * ( 1. + ssh/depth ) ! Nominal dz*
          dz = max(dz, 0.)              ! Avoid negative incase ssh=-depth
          dz = min(dz, depth - z)       ! Clip if below topography
          z = z + dz                    ! Bottom of layer
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

!> Parses a string and generates a rho_target(:) profile with refined resolution downward
!! and returns the number of levels
integer function rho_function1( string, rho_target )
  character(len=*),   intent(in)    :: string !< String with list of parameters in form
                                              !! dz_min, H_total, power, precision
  real, dimension(:), allocatable, intent(inout) :: rho_target !< Profile of interface densities [kg m-3]
  ! Local variables
  integer :: nki, k, nk
  real    :: ddx, dx, rho_1, rho_2, rho_3, drho, rho_4, drho_min

  read( string, *) nk, rho_1, rho_2, rho_3, drho, rho_4, drho_min
  allocate(rho_target(nk+1))
  nki = nk + 1 - 4 ! Number of interfaces minus 4 specified values
  rho_target(1) = rho_1
  rho_target(2) = rho_2
  dx = 0.
  do k = 0, nki
    ddx = max( drho_min, real(nki-k)/real(nki*nki) )
    dx = dx + ddx
    rho_target(3+k) = rho_3 + (2. * drho) * dx
  enddo
  rho_target(nki+4) = rho_4

  rho_function1 = nk

end function rho_function1

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
