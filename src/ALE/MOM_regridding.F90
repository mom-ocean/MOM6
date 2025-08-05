!> Generates vertical grids as part of the ALE algorithm
module MOM_regridding

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING, assert
use MOM_file_parser,   only : param_file_type, get_param, log_param
use MOM_io,            only : file_exists, field_exists, field_size, MOM_read_data
use MOM_io,            only : vardesc, var_desc, SINGLE_FILE
use MOM_io,            only : MOM_netCDF_file, MOM_field
use MOM_io,            only : create_MOM_file, MOM_write_field
use MOM_io,            only : verify_variable_units, slasher
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
use regrid_consts, only : REGRIDDING_HYCOM1, REGRIDDING_HYBGEN, REGRIDDING_ADAPTIVE
use regrid_interp, only : interp_CS_type, set_interp_scheme, set_interp_extrap, set_interp_answer_date

use coord_zlike,  only : init_coord_zlike, zlike_CS, set_zlike_params, build_zstar_column, end_coord_zlike
use coord_sigma,  only : init_coord_sigma, sigma_CS, set_sigma_params, build_sigma_column, end_coord_sigma
use coord_rho,    only : init_coord_rho, rho_CS, set_rho_params, build_rho_column, end_coord_rho
use coord_rho,    only : old_inflate_layers_1d
use coord_hycom,  only : init_coord_hycom, hycom_CS, set_hycom_params, build_hycom1_column, end_coord_hycom
use coord_adapt,  only : init_coord_adapt, adapt_CS, set_adapt_params, build_adapt_column, end_coord_adapt
use MOM_hybgen_regrid, only : hybgen_regrid, hybgen_regrid_CS, init_hybgen_regrid, end_hybgen_regrid
use MOM_hybgen_regrid, only : write_Hybgen_coord_file

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
  !! [Z ~> m] for z*, [nondim] for sigma, etc.
  real, dimension(:), allocatable :: coordinateResolution

  !> This is a scaling factor that restores coordinateResolution to values in
  !! the natural units for output, perhaps [nondim]
  real :: coord_scale = 1.0

  !> This array is set by function set_target_densities()
  !! This array is the nominal coordinate of interfaces and is the
  !! running sum of coordinateResolution, in [R ~> kg m-3]. i.e.
  !!  target_density(k+1) = coordinateResolution(k) + coordinateResolution(k)
  !! It is only used in "rho" or "Hycom" mode.
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
  !! profiles when interpolating for target grid positions [nondim]
  real :: compressibility_fraction = 0.

  !> If true, each interface is given a maximum depth based on a rescaling of
  !! the indexing of coordinateResolution.
  logical :: set_maximum_depths = .false.

  !> If true, integrate for interface positions from the top downward.
  !! If false, integrate from the bottom upward, as does the rest of the model.
  logical :: integrate_downward_for_e = .true.

  !> The vintage of the order of arithmetic and expressions to use for remapping.
  !! Values below 20190101 recover the remapping answers from 2018.
  !! Higher values use more robust forms of the same remapping expressions.
  integer :: remap_answer_date = 99991231

  logical :: use_hybgen_unmix = .false.  !< If true, use the hybgen unmixing code before remapping

  type(zlike_CS),  pointer :: zlike_CS  => null() !< Control structure for z-like coordinate generator
  type(sigma_CS),  pointer :: sigma_CS  => null() !< Control structure for sigma coordinate generator
  type(rho_CS),    pointer :: rho_CS    => null() !< Control structure for rho coordinate generator
  type(hycom_CS),  pointer :: hycom_CS  => null() !< Control structure for hybrid coordinate generator
  type(adapt_CS),  pointer :: adapt_CS  => null() !< Control structure for adaptive coordinate generator
  type(hybgen_regrid_CS), pointer :: hybgen_CS => NULL() !< Control structure for hybgen regridding

end type

! The following routines are visible to the outside world
public initialize_regridding, end_regridding, regridding_main
public regridding_preadjust_reqs, convective_adjustment
public inflate_vanished_layers_old, check_grid_column
public set_regrid_params, get_regrid_size, write_regrid_file
public uniformResolution, setCoordinateResolution
public set_target_densities_from_GV, set_target_densities
public set_regrid_max_depths, set_regrid_max_thickness
public getCoordinateResolution, getCoordinateInterfaces
public getCoordinateUnits, getCoordinateShortName, getStaticThickness
public DEFAULT_COORDINATE_MODE
public set_h_neglect, set_dz_neglect
public get_zlike_CS, get_sigma_CS, get_rho_CS

!> Documentation for coordinate options
character(len=*), parameter, public :: regriddingCoordinateModeDoc = &
                 " LAYER - Isopycnal or stacked shallow water layers\n"//&
                 " ZSTAR, Z* - stretched geopotential z*\n"//&
                 " SIGMA_SHELF_ZSTAR - stretched geopotential z* ignoring shelf\n"//&
                 " SIGMA - terrain following coordinates\n"//&
                 " RHO   - continuous isopycnal\n"//&
                 " HYCOM1 - HyCOM-like hybrid coordinate\n"//&
                 " HYBGEN - Hybrid coordinate from the Hycom hybgen code\n"//&
                 " ADAPTIVE - optimize for smooth neutral density surfaces"

!> Documentation for regridding interpolation schemes
character(len=*), parameter, public :: regriddingInterpSchemeDoc = &
                 " P1M_H2     (2nd-order accurate)\n"//&
                 " P1M_H4     (2nd-order accurate)\n"//&
                 " P1M_IH4    (2nd-order accurate)\n"//&
                 " PLM        (2nd-order accurate)\n"//&
                 " PPM_CW     (3rd-order accurate)\n"//&
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
!> Default minimum thickness for some coordinate generation modes [m]
real, parameter, public :: regriddingDefaultMinThickness = 1.e-3

!> Maximum length of parameters
integer, parameter :: MAX_PARAM_LENGTH = 120

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
  character(len=40)  :: coord_units, coord_res_param ! Temporary strings
  character(len=MAX_PARAM_LENGTH) :: param_name
  character(len=200) :: inputdir, fileName, longString
  character(len=320) :: message ! Temporary strings
  character(len=12) :: expected_units, alt_units ! Temporary strings
  logical :: tmpLogical, do_sum, main_parameters
  logical :: coord_is_state_dependent, ierr
  integer :: default_answer_date  ! The default setting for the various ANSWER_DATE flags.
  integer :: remap_answer_date    ! The vintage of the remapping expressions to use.
  integer :: regrid_answer_date   ! The vintage of the regridding expressions to use.
  real :: tmpReal  ! A temporary variable used in setting other variables [various]
  real :: P_Ref    ! The coordinate variable reference pression [R L2 T-2 ~> Pa]
  real :: maximum_depth ! The maximum depth of the ocean [m] (not in Z).
  real :: dz_extra      ! The thickness of an added layer to append to the woa09_dz profile when
                        ! maximum_depth is large [m] (not in Z).
  real :: adaptTimeRatio, adaptZoomCoeff ! Temporary variables for input parameters [nondim]
  real :: adaptBuoyCoeff, adaptAlpha     ! Temporary variables for input parameters [nondim]
  real :: adaptZoom  ! The thickness of the near-surface zooming region with the adaptive coordinate [H ~> m or kg m-2]
  real :: adaptDrho0 ! Reference density difference for stratification-dependent diffusion. [R ~> kg m-3]
  integer :: k, nzf(4)
  real, dimension(:), allocatable :: dz     ! Resolution (thickness) in units of coordinate, which may be [m]
                                            ! or [Z ~> m] or [H ~> m or kg m-2] or [R ~> kg m-3] or other units.
  real, dimension(:), allocatable :: h_max  ! Maximum layer thicknesses [H ~> m or kg m-2]
  real, dimension(:), allocatable :: z_max  ! Maximum interface depths [H ~> m or kg m-2] or other
                                            ! units depending on the coordinate
  real, dimension(:), allocatable :: dz_max ! Thicknesses used to find maximum interface depths
                                            ! [H ~> m or kg m-2] or other units
  real, dimension(:), allocatable :: rho_target ! Target density used in HYBRID mode [kg m-3]
  ! Thicknesses [m] that give level centers approximately corresponding to table 2 of WOA09
  ! These are approximate because the WOA09 depths are not smoothly spaced. Levels
  ! 1, 4, 5, 9, 12, 24, and 36 are 2.5, 2.5, 1.25 12.5, 37.5 and 62.5 m deeper than WOA09
  ! but all others are identical.
  real, dimension(40) :: woa09_dz_approx = (/ 5.,  10.,  10.,  15.,  22.5, 25.,  25.,  25.,  &
                                             37.5, 50.,  50.,  75., 100., 100., 100., 100., &
                                            100., 100., 100., 100., 100., 100., 100., 175., &
                                            250., 375., 500., 500., 500., 500., 500., 500., &
                                            500., 500., 500., 500., 500., 500., 500., 500. /)
  ! These are the actual spacings [m] between WOA09 depths which, if used for layer thickness, places
  ! the interfaces at the WOA09 depths.
  real, dimension(39) :: woa09_dzi = (/ 10.,  10.,  10.,  20.,  25.,  25.,  25.,  25.,  &
                                        50.,  50.,  50., 100., 100., 100., 100., 100., &
                                       100., 100., 100., 100., 100., 100., 100., 250., &
                                       250., 500., 500., 500., 500., 500., 500., 500., &
                                       500., 500., 500., 500., 500., 500., 500. /)
  ! These are the spacings [m] between WOA23 depths from table 3 of
  ! https://www.ncei.noaa.gov/data/oceans/woa/WOA13/DOC/woa13documentation.pdf
  real, dimension(136) :: woa23_dzi = (/ 5.,   5.,   5.,   5.,   5.,   5.,   5.,   5.,   5.,   5., &
                                         5.,   5.,   5.,   5.,   5.,   5.,   5.,   5.,   5.,   5., &
                                        25.,  25.,  25.,  25.,  25.,  25.,  25.,  25.,  25.,  25., &
                                        25.,  25.,  25.,  25.,  25.,  25.,  50.,  50.,  50.,  50., &
                                        50.,  50.,  50.,  50.,  50.,  50.,  50.,  50.,  50.,  50., &
                                        50.,  50.,  50.,  50.,  50.,  50.,  50.,  50.,  50.,  50., &
                                        50.,  50.,  50.,  50.,  50.,  50., 100., 100., 100., 100., &
                                       100., 100., 100., 100., 100., 100., 100., 100., 100., 100., &
                                       100., 100., 100., 100., 100., 100., 100., 100., 100., 100., &
                                       100., 100., 100., 100., 100., 100., 100., 100., 100., 100., &
                                       100., 100., 100., 100., 100., 100., 100., 100., 100., 100., &
                                       100., 100., 100., 100., 100., 100., 100., 100., 100., 100., &
                                       100., 100., 100., 100., 100., 100., 100., 100., 100., 100., &
                                       100., 100., 100., 100., 100., 100. /)

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
      param_name = create_coord_param(param_prefix, "INTERP_SCHEME", param_suffix)
      string2 = 'PPM_H4' ! Default for diagnostics
    endif
    call get_param(param_file, mdl, param_name, string, &
                 "This sets the interpolation scheme to use to "//&
                 "determine the new grid. These parameters are "//&
                 "only relevant when REGRIDDING_COORDINATE_MODE is "//&
                 "set to a function of state. Otherwise, it is not "//&
                 "used. It can be one of the following schemes: \n"//&
                 trim(regriddingInterpSchemeDoc), default=trim(string2))
    call set_regrid_params(CS, interp_scheme=string)

    call get_param(param_file, mdl, "DEFAULT_ANSWER_DATE", default_answer_date, &
                 "This sets the default value for the various _ANSWER_DATE parameters.", &
                 default=99991231)
    call get_param(param_file, mdl, "REMAPPING_ANSWER_DATE", remap_answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for remapping.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
    if (.not.GV%Boussinesq) remap_answer_date = max(remap_answer_date, 20230701)
    call set_regrid_params(CS, remap_answer_date=remap_answer_date)
    call get_param(param_file, mdl, "REGRIDDING_ANSWER_DATE", regrid_answer_date, &
                 "The vintage of the expressions and order of arithmetic to use for regridding.  "//&
                 "Values below 20190101 result in the use of older, less accurate expressions "//&
                 "that were in use at the end of 2018.  Higher values result in the use of more "//&
                 "robust and accurate forms of mathematically equivalent expressions.", &
                 default=default_answer_date, do_not_log=.not.GV%Boussinesq)
    if (.not.GV%Boussinesq) regrid_answer_date = max(regrid_answer_date, 20230701)
    call set_regrid_params(CS, regrid_answer_date=regrid_answer_date)
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
    param_name = create_coord_param(param_prefix, "DEF", param_suffix)
    coord_res_param = create_coord_param(param_prefix, "RES", param_suffix)
    string2 = 'UNIFORM'
    if ((maximum_depth>3000.) .and. (maximum_depth<9250.)) string2='WOA09' ! For convenience
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
                 " WOA09INT[:N] - layers spanned by the WOA09 depths\n"//&
                 " WOA23INT[:N] - layers spanned by the WOA23 depths\n"//&
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
      expected_units = 'nondim' ; alt_units = expected_units
    elseif (CS%regridding_scheme == REGRIDDING_RHO) then
      expected_units = 'kg m-3' ; alt_units = expected_units
    else
      expected_units = 'meters' ; alt_units = 'm'
    endif
    if (index(trim(varName),'interfaces=')==1) then
      varName=trim(varName(12:))
      call verify_variable_units(filename, varName, expected_units, message, ierr, alt_units)
      if (ierr) call MOM_error(FATAL, trim(mdl)//", initialize_regridding: "//&
                  "Unsupported format in grid definition '"//trim(filename)//"'. Error message "//trim(message))
      call field_size(trim(fileName), trim(varName), nzf)
      ke = nzf(1)-1
      if (ke < 1) call MOM_error(FATAL, trim(mdl)//" initialize_regridding via Var "//trim(varName)//&
                                 "in FILE "//trim(filename)//" requires at least 2 target interface values.")
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
    if (main_parameters .and. (ke/=GV%ke)) then
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
  elseif (index(trim(string),'WOA09INT')==1) then
    if (len_trim(string)==8) then ! string=='WOA09INT'
      tmpReal = 0. ; ke = 0 ; dz_extra = 0.
      do while (tmpReal<maximum_depth)
        ke = ke + 1
        if (ke > size(woa09_dzi)) then
          dz_extra = maximum_depth - tmpReal
          exit
        endif
        tmpReal = tmpReal + woa09_dzi(ke)
      enddo
    elseif (index(trim(string),'WOA09INT:')==1) then ! string starts with 'WOA09INT:'
      if (len_trim(string)==9) call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
                 'Expected string of form "WOA09INT:N" but got "'//trim(string)//'".')
      ke = extract_integer(string(10:len_trim(string)),'',1)
      if (ke>39 .or. ke<1) call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
                   'For "WOA05INT:N" N must 0<N<40 but got "'//trim(string)//'".')
    endif
    allocate(dz(ke))
    do k=1,min(ke, size(woa09_dzi))
      dz(k) = woa09_dzi(k)
    enddo
    if (ke > size(woa09_dzi)) dz(ke) = dz_extra
  elseif (index(trim(string),'WOA23INT')==1) then
    if (len_trim(string)==8) then ! string=='WOA23INT'
      tmpReal = 0. ; ke = 0 ; dz_extra = 0.
      do while (tmpReal<maximum_depth)
        ke = ke + 1
        if (ke > size(woa23_dzi)) then
          dz_extra = maximum_depth - tmpReal
          exit
        endif
        tmpReal = tmpReal + woa23_dzi(ke)
      enddo
    elseif (index(trim(string),'WOA23INT:')==1) then ! string starts with 'WOA23INT:'
      if (len_trim(string)==9) call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
                 'Expected string of form "WOA23INT:N" but got "'//trim(string)//'".')
      ke = extract_integer(string(10:len_trim(string)),'',1)
      if (ke>39 .or. ke<1) call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
                   'For "WOA05INT:N" N must 0<N<40 but got "'//trim(string)//'".')
    endif
    allocate(dz(ke))
    do k=1,min(ke, size(woa23_dzi))
      dz(k) = woa23_dzi(k)
    enddo
    if (ke > size(woa23_dzi)) dz(ke) = dz_extra
  elseif (index(trim(string),'WOA09')==1) then
    if (len_trim(string)==5) then ! string=='WOA09'
      tmpReal = 0. ; ke = 0 ; dz_extra = 0.
      do while (tmpReal<maximum_depth)
        ke = ke + 1
        if (ke > size(woa09_dz_approx)) then
          dz_extra = maximum_depth - tmpReal
          exit
        endif
        tmpReal = tmpReal + woa09_dz_approx(ke)
      enddo
    elseif (index(trim(string),'WOA09:')==1) then ! string starts with 'WOA09:'
      if (len_trim(string)==6) call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
                 'Expected string of form "WOA09:N" but got "'//trim(string)//'".')
      ke = extract_integer(string(7:len_trim(string)),'',1)
      if (ke>40 .or. ke<1) call MOM_error(FATAL,trim(mdl)//', initialize_regridding: '// &
                   'For "WOA05:N" N must 0<N<41 but got "'//trim(string)//'".')
    endif
    allocate(dz(ke))
    do k=1,min(ke, size(woa09_dz_approx))
      dz(k) = woa09_dz_approx(k)
    enddo
    if (ke > size(woa09_dz_approx)) dz(ke) = dz_extra
  else
    call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
      "Unrecognized coordinate configuration"//trim(string))
  endif

  if (main_parameters) then
    ! This is a work around to apparently needed to work with the from_Z initialization...  ???
    if (coordinateMode(coord_mode) == REGRIDDING_ZSTAR .or. &
        coordinateMode(coord_mode) == REGRIDDING_HYCOM1 .or. &
        coordinateMode(coord_mode) == REGRIDDING_HYBGEN .or. &
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
  allocate( CS%coordinateResolution(CS%nk), source=-1.E30 )
  if (state_dependent(CS%regridding_scheme)) then
    ! Target values
    allocate( CS%target_density(CS%nk+1), source=-1.E30*US%kg_m3_to_R )
  endif

  if (allocated(dz)) then
    if (coordinateMode(coord_mode) == REGRIDDING_SIGMA) then
      call setCoordinateResolution(dz, CS, scale=1.0)
    elseif (coordinateMode(coord_mode) == REGRIDDING_RHO) then
      call setCoordinateResolution(dz, CS, scale=US%kg_m3_to_R)
    elseif (coordinateMode(coord_mode) == REGRIDDING_ADAPTIVE) then
      call setCoordinateResolution(dz, CS, scale=GV%m_to_H)
      CS%coord_scale = GV%H_to_m
    else
      call setCoordinateResolution(dz, CS, scale=US%m_to_Z)
      CS%coord_scale = US%Z_to_m
    endif
  endif

  ! set coord_scale for RHO regridding independent of allocation status of dz
  if (coordinateMode(coord_mode) == REGRIDDING_RHO) then
    CS%coord_scale = US%R_to_kg_m3
  endif

  ! ensure CS%ref_pressure is rescaled properly
  CS%ref_pressure = US%Pa_to_RL2_T2 * CS%ref_pressure

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
  call initCoord(CS, GV, US, coord_mode, param_file)

  if (coord_is_state_dependent) then
    if (main_parameters) then
      call get_param(param_file, mdl, create_coord_param(param_prefix, "P_REF", param_suffix), P_Ref, &
                   "The pressure that is used for calculating the coordinate "//&
                   "density.  (1 Pa = 1e4 dbar, so 2e7 is commonly used.) "//&
                   "This is only used if USE_EOS and ENABLE_THERMODYNAMICS are true.", &
                   units="Pa", default=2.0e7, scale=US%Pa_to_RL2_T2)
    else
      call get_param(param_file, mdl, create_coord_param(param_prefix, "P_REF", param_suffix), P_Ref, &
                   "The pressure that is used for calculating the diagnostic coordinate "//&
                   "density.  (1 Pa = 1e4 dbar, so 2e7 is commonly used.) "//&
                   "This is only used for the RHO coordinate.", &
                   units="Pa", default=2.0e7, scale=US%Pa_to_RL2_T2)
    endif
    call get_param(param_file, mdl, create_coord_param(param_prefix, "REGRID_COMPRESSIBILITY_FRACTION", param_suffix), &
                 tmpReal, &
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

  if (main_parameters .and. coordinateMode(coord_mode) == REGRIDDING_HYCOM1) then
    call get_param(param_file, mdl, "HYCOM1_ONLY_IMPROVES", tmpLogical, &
              "When regridding, an interface is only moved if this improves the fit to the target density.", &
              default=.false.)
    call set_hycom_params(CS%hycom_CS, only_improves=tmpLogical)
  endif

  CS%use_hybgen_unmix = .false.
  if (coordinateMode(coord_mode) == REGRIDDING_HYBGEN) then
    call get_param(param_file, mdl, "USE_HYBGEN_UNMIX", CS%use_hybgen_unmix, &
              "If true, use hybgen unmixing code before regridding.", &
              default=.false.)
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
    call get_param(param_file, mdl, "MAX_LAYER_THICKNESS_CONFIG", longString, &
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
    if ( trim(longString) == "NONE") then
      ! Do nothing.
    elseif ( trim(longString) ==  "PARAM") then
      call get_param(param_file, mdl, "MAX_LAYER_THICKNESS", h_max, &
                   trim(message), units="m", fail_if_missing=.true., scale=GV%m_to_H)
      call set_regrid_max_thickness(CS, h_max)
    elseif (index(trim(longString),'FILE:')==1) then
      if (longString(6:6)=='.' .or. longString(6:6)=='/') then
        ! If we specified "FILE:./xyz" or "FILE:/xyz" then we have a relative or absolute path
        fileName = trim( extractWord(trim(longString(6:200)), 1) )
      else
        ! Otherwise assume we should look for the file in INPUTDIR
        fileName = trim(inputdir) // trim( extractWord(trim(longString(6:200)), 1) )
      endif
      if (.not. file_exists(fileName)) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
        "Specified file not found: Looking for '"//trim(fileName)//"' ("//trim(longString)//")")

      varName = trim( extractWord(trim(longString(6:)), 2) )
      if (.not. field_exists(fileName,varName)) call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
        "Specified field not found: Looking for '"//trim(varName)//"' ("//trim(longString)//")")
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
    elseif (index(trim(longString),'FNC1:')==1) then
      call dz_function1( trim(longString(6:)), h_max )
      call log_param(param_file, mdl, "!MAX_LAYER_THICKNESS", h_max, &
                 trim(message), units=coordinateUnits(coord_mode))
      call set_regrid_max_thickness(CS, h_max, GV%m_to_H)
    else
      call MOM_error(FATAL,trim(mdl)//", initialize_regridding: "// &
        "Unrecognized MAX_LAYER_THICKNESS_CONFIG "//trim(longString))
    endif
    deallocate(h_max)
  endif

  if (allocated(dz)) deallocate(dz)
end subroutine initialize_regridding



!> Deallocation of regridding memory
subroutine end_regridding(CS)
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure

  if (associated(CS%zlike_CS))  call end_coord_zlike(CS%zlike_CS)
  if (associated(CS%sigma_CS))  call end_coord_sigma(CS%sigma_CS)
  if (associated(CS%rho_CS))    call end_coord_rho(CS%rho_CS)
  if (associated(CS%hycom_CS))  call end_coord_hycom(CS%hycom_CS)
  if (associated(CS%adapt_CS))  call end_coord_adapt(CS%adapt_CS)
  if (associated(CS%hybgen_CS)) call end_hybgen_regrid(CS%hybgen_CS)

  deallocate( CS%coordinateResolution )
  if (allocated(CS%target_density)) deallocate( CS%target_density )
  if (allocated(CS%max_interface_depths) ) deallocate( CS%max_interface_depths )
  if (allocated(CS%max_layer_thickness) ) deallocate( CS%max_layer_thickness )

end subroutine end_regridding

!------------------------------------------------------------------------------
!> Dispatching regridding routine for orchestrating regridding & remapping
subroutine regridding_main( remapCS, CS, G, GV, US, h, tv, h_new, dzInterface, &
                            frac_shelf_h, PCM_cell)
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
  type(unit_scale_type),                      intent(in)    :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h      !< Current 3D grid obtained after
                                                                      !! the last time step [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)    :: tv     !< Thermodynamical variables (T, S, ...)
  real, dimension(SZI_(G),SZJ_(G),CS%nk),     intent(inout) :: h_new  !< New 3D grid consistent with target
                                                                      !! coordinate [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),   intent(inout) :: dzInterface !< The change in position of each
                                                                      !! interface [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in   ) :: frac_shelf_h !< Fractional ice shelf coverage [nomdim]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                    optional, intent(out  ) :: PCM_cell !< Use PCM remapping in cells where true

  ! Local variables
  real :: nom_depth_H(SZI_(G),SZJ_(G))  !< The nominal ocean depth at each point in thickness units [H ~> m or kg m-2]
  real :: tot_h(SZI_(G),SZJ_(G))  !< The total thickness of the water column [H ~> m or kg m-2]
  real :: tot_dz(SZI_(G),SZJ_(G)) !< The total distance between the top and bottom of the water column [Z ~> m]
  real :: Z_to_H  ! A conversion factor used by some routines to convert coordinate
                  ! parameters to depth units [H Z-1 ~> nondim or kg m-3]
  character(len=128) :: mesg    ! A string for error messages
  integer :: i, j, k

  if (present(PCM_cell)) PCM_cell(:,:,:) = .false.

  Z_to_H = US%Z_to_m * GV%m_to_H  ! Often this is equivalent to GV%Z_to_H.

  if ((allocated(tv%SpV_avg)) .and. (tv%valid_SpV_halo < 1)) then
    if (tv%valid_SpV_halo < 0) then
      mesg = "invalid values of SpV_avg."
    else
      mesg = "insufficiently large SpV_avg halos of width 0 but 1 is needed."
    endif
    call MOM_error(FATAL, "Regridding_main called in fully non-Boussinesq mode with "//trim(mesg))
  endif

  if (allocated(tv%SpV_avg)) then  ! This is the fully non-Boussinesq case
    do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      tot_h(i,j) = 0.0 ; tot_dz(i,j) = 0.0
    enddo ; enddo
    do k=1,GV%ke ; do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      tot_h(i,j) = tot_h(i,j) + h(i,j,k)
      tot_dz(i,j) = tot_dz(i,j) + GV%H_to_RZ * tv%SpV_avg(i,j,k) * h(i,j,k)
    enddo ; enddo ; enddo
    do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      if ((tot_dz(i,j) > 0.0) .and. (G%bathyT(i,j)+G%Z_ref > 0.0)) then
        nom_depth_H(i,j) = (G%bathyT(i,j)+G%Z_ref) * (tot_h(i,j) / tot_dz(i,j))
      else
        nom_depth_H(i,j) = 0.0
      endif
    enddo ; enddo
  else
    do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      nom_depth_H(i,j) = max((G%bathyT(i,j)+G%Z_ref) * Z_to_H, 0.0)
    enddo ; enddo
  endif

  select case ( CS%regridding_scheme )

    case ( REGRIDDING_ZSTAR )
      call build_zstar_grid( CS, G, GV, h, nom_depth_H, dzInterface, frac_shelf_h, zScale=Z_to_H )
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)
    case ( REGRIDDING_SIGMA_SHELF_ZSTAR)
      call build_zstar_grid( CS, G, GV, h, nom_depth_H, dzInterface, zScale=Z_to_H )
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)
    case ( REGRIDDING_SIGMA )
      call build_sigma_grid( CS, G, GV, h, nom_depth_H, dzInterface )
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)
    case ( REGRIDDING_RHO )
      call build_rho_grid( G, GV, G%US, h, nom_depth_H, tv, dzInterface, remapCS, CS, frac_shelf_h )
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)
    case ( REGRIDDING_HYCOM1 )
      call build_grid_HyCOM1( G, GV, G%US, h, nom_depth_H, tv, h_new, dzInterface, remapCS, CS, &
                              frac_shelf_h, zScale=Z_to_H )
    case ( REGRIDDING_HYBGEN )
      call hybgen_regrid(G, GV, G%US, h, nom_depth_H, tv, CS%hybgen_CS, dzInterface, PCM_cell)
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)
    case ( REGRIDDING_ADAPTIVE )
      call build_grid_adaptive(G, GV, G%US, h, nom_depth_H, tv, dzInterface, remapCS, CS)
      call calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)

    case ( REGRIDDING_ARBITRARY )
      call MOM_error(FATAL,'MOM_regridding, regridding_main: '//&
                     'Regridding mode "ARB" is not implemented.')
    case default
      call MOM_error(FATAL,'MOM_regridding, regridding_main: '//&
                     'Unknown regridding scheme selected!')

  end select ! type of grid

#ifdef __DO_SAFETY_CHECKS__
  if (CS%nk == GV%ke) then
    do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1 ; if (G%mask2dT(i,j)>0.) then
      call check_grid_column( GV%ke, h(i,j,:), dzInterface(i,j,:), 'in regridding_main')
    endif ; enddo ; enddo
  endif
#endif
  do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (G%mask2dT(i,j) > 0.) then
    if (minval(h(i,j,:)) < 0.0) then
      write(0,*) 'regridding_main check_grid: i,j=', i, j, 'h_new(i,j,:)=', h_new(i,j,:)
      call MOM_error(FATAL, "regridding_main: negative thickness encountered.")
    endif
  endif ; enddo ; enddo

end subroutine regridding_main

!------------------------------------------------------------------------------
!> This routine returns flags indicating which pre-remapping state adjustments
!! are needed depending on the coordinate mode in use.
subroutine regridding_preadjust_reqs(CS, do_conv_adj, do_hybgen_unmix, hybgen_CS)

  ! Arguments
  type(regridding_CS), intent(in)  :: CS          !< Regridding control structure
  logical,             intent(out) :: do_conv_adj !< Convective adjustment should be done
  logical,             intent(out) :: do_hybgen_unmix !< Hybgen unmixing should be done
  type(hybgen_regrid_CS), pointer, &
             optional, intent(out) :: hybgen_CS   !< Control structure for hybgen regridding for sharing parameters.


  do_conv_adj = .false. ; do_hybgen_unmix = .false.
  select case ( CS%regridding_scheme )

    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA_SHELF_ZSTAR, REGRIDDING_SIGMA, REGRIDDING_ARBITRARY, &
           REGRIDDING_HYCOM1, REGRIDDING_ADAPTIVE )
      do_conv_adj = .false. ; do_hybgen_unmix = .false.
    case ( REGRIDDING_RHO )
      do_conv_adj = .true. ; do_hybgen_unmix = .false.
    case ( REGRIDDING_HYBGEN )
      do_conv_adj = .false. ; do_hybgen_unmix = CS%use_hybgen_unmix
    case default
      call MOM_error(FATAL,'MOM_regridding, regridding_preadjust_reqs: '//&
                     'Unknown regridding scheme selected!')
  end select ! type of grid

  if (present(hybgen_CS) .and. do_hybgen_unmix) hybgen_CS => CS%hybgen_CS

end subroutine regridding_preadjust_reqs


!> Calculates h_new from h + delta_k dzInterface
subroutine calc_h_new_by_dz(CS, G, GV, h, dzInterface, h_new)
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                     intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Old layer thicknesses [H ~> m or kg m-2]
                                                                 !! or other units
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),  intent(in)    :: dzInterface !< Change in interface positions
                                                                 !! in the same units as h [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk),    intent(inout) :: h_new !< New layer thicknesses in the same
                                                                 !! units as h [H ~> m or kg m-2]
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


!> Check that the total thickness of new and old grids are consistent
subroutine check_grid_column( nk, h, dzInterface, msg )
  integer,               intent(in) :: nk !< Number of cells
  real, dimension(nk),   intent(in) :: h  !< Cell thicknesses [Z ~> m] or arbitrary units
  real, dimension(nk+1), intent(in) :: dzInterface !< Change in interface positions (same units as h), often [Z ~> m]
  character(len=*),      intent(in) :: msg !< Message to append to errors
  ! Local variables
  integer :: k
  real :: eps          ! A tiny relative thickness [nondim]
  real :: total_h_old  ! The total thickness in the old column, in [Z ~> m] or arbitrary units
  real :: total_h_new  ! The total thickness in the updated column, in [Z ~> m] or arbitrary units
  real :: h_new        ! A thickness in the updated column, in [Z ~> m] or arbitrary units

  eps =1. ; eps = epsilon(eps)

  ! Total thickness of grid h
  total_h_old = 0.
  do k = 1,nk
    total_h_old = total_h_old + h(k)
  enddo

  total_h_new = 0.
  do k = nk,1,-1
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
  real, dimension(CS%nk+1), intent(in)    :: z_new !< New grid position before filtering [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: dz_g  !< Change in interface positions including
                                                   !! the effects of filtering [H ~> m or kg m-2]
  ! Local variables
  real :: sgn     ! The sign convention for downward [nondim].
  real :: dz_tgt  ! The target grid movement of the unfiltered grid [H ~> m or kg m-2]
  real :: zr1     ! The old grid position of an interface relative to the surface [H ~> m or kg m-2]
  real :: z_old_k ! The corrected position of the old grid [H ~> m or kg m-2]
  real :: Aq      ! A temporary variable related to the grid weights [nondim]
  real :: Bq      ! A temporary variable used in the linear term in the quadratic expression for the
                  ! filtered grid movement [H ~> m or kg m-2]
  real :: z0, dz0 ! Together these give the position of an interface relative to a reference hieght
                  ! that may be adjusted for numerical accuracy in a solver [H ~> m or kg m-2]
  real :: F0      ! An estimated grid movement  [H ~> m or kg m-2]
  real :: zs      ! The depth at which the shallow filtering timescale applies [H ~> m or kg m-2]
  real :: zd      ! The depth at which the deep filtering timescale applies [H ~> m or kg m-2]
  real :: dzwt    ! The depth range over which the transition in the filtering timescale occurs [H ~> m or kg m-2]
  real :: Idzwt   ! The Adcroft reciprocal of dzwt [H-1 ~> m-1 or m2 kg-1]
  real :: wtd     ! The weight given to the new grid when time filtering [nondim]
  real :: Iwtd    ! The inverse of wtd [nondim]
  real :: Int_zs  ! A depth integral of the weights in [H ~> m or kg m-2]
  real :: Int_zd  ! A depth integral of the weights in [H ~> m or kg m-2]
  real :: dInt_zs_zd ! The depth integral of the weights between the deep and shallow depths in [H ~> m or kg m-2]
! For debugging:
  real, dimension(nk+1) :: z_act ! The final grid positions after the filtered movement [H ~> m or kg m-2]
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

!> Builds a z*-coordinate grid with partial steps (Adcroft and Campin, 2004).
!! z* is defined as
!!   z* = (z-eta)/(H+eta)*H  s.t. z*=0 when z=eta and z*=-H when z=-H .
subroutine build_zstar_grid( CS, G, GV, h, nom_depth_H, dzInterface, frac_shelf_h, zScale)

  ! Arguments
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                     intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                 !! relative to mean sea level or another locally
                                                                 !! valid reference height, converted to thickness
                                                                 !! units [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),  intent(inout) :: dzInterface !< The change in interface depth
                                                                 !! [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)), optional,intent(in)    :: frac_shelf_h !< Fractional
                                                                 !! ice shelf coverage [nondim].
  real,                            optional, intent(in)    :: zScale !< Scaling factor from the target coordinate
                                                                 !! resolution in Z to desired units for zInterface,
                                                                 !! usually Z_to_H in which case it is in
                                                                 !! units of [H Z-1 ~> nondim or kg m-3]
  ! Local variables
  real   :: nominalDepth, minThickness, totalThickness  ! Depths and thicknesses [H ~> m or kg m-2]
#ifdef __DO_SAFETY_CHECKS__
  real :: dh    ! The larger of the total column thickness or bathymetric depth [H ~> m or kg m-2]
#endif
  real, dimension(SZK_(GV)+1) :: zOld    ! Previous coordinate interface heights [H ~> m or kg m-2]
  real, dimension(CS%nk+1)    :: zNew    ! New coordinate interface heights [H ~> m or kg m-2]
  integer :: i, j, k, nz
  logical :: ice_shelf

  nz = GV%ke
  minThickness = CS%min_thickness
  ice_shelf = present(frac_shelf_h)

  !$OMP parallel do default(none) shared(G,GV,dzInterface,CS,nz,h,frac_shelf_h, &
  !$OMP                                  ice_shelf,minThickness,zScale,nom_depth_H) &
  !$OMP                          private(nominalDepth,totalThickness, &
#ifdef __DO_SAFETY_CHECKS__
  !$OMP                                  dh, &
#endif
  !$OMP                                  zNew,zOld)
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      if (G%mask2dT(i,j)==0.) then
        dzInterface(i,j,:) = 0.
        cycle
      endif

      ! Local depth (positive downward)
      nominalDepth = nom_depth_H(i,j)

      ! Determine water column thickness
      totalThickness = 0.0
      do k = 1,nz
        totalThickness = totalThickness + h(i,j,k)
      enddo

      ! if (GV%Boussinesq) then
      zOld(nz+1) = - nominalDepth
      do k = nz,1,-1
        zOld(k) = zOld(k+1) + h(i,j,k)
      enddo
      ! else   ! Work downward?
      ! endif

      if (ice_shelf) then
        if (frac_shelf_h(i,j) > 0.) then ! under ice shelf
          call build_zstar_column(CS%zlike_CS, nominalDepth, totalThickness, zNew, &
                                z_rigid_top=totalThickness-nominalDepth, &
                                eta_orig=zOld(1), zScale=zScale)
        else
          call build_zstar_column(CS%zlike_CS, nominalDepth, totalThickness, &
                                zNew, zScale=zScale)
        endif
      else
        call build_zstar_column(CS%zlike_CS, nominalDepth, totalThickness, &
                                zNew, zScale=zScale)
      endif

      ! Calculate the final change in grid position after blending new and old grids
      call filtered_grid_motion( CS, nz, zOld, zNew, dzInterface(i,j,:) )

#ifdef __DO_SAFETY_CHECKS__
      dh = max(nominalDepth,totalThickness)
      if (abs(zNew(1)-zOld(1)) > (nz-1)*0.5*epsilon(dh)*dh) then
        write(0,*) 'min_thickness=',CS%min_thickness
        write(0,*) 'nominalDepth=',nominalDepth,'totalThickness=',totalThickness
        write(0,*) 'dzInterface(1) = ', dzInterface(i,j,1), epsilon(dh), nz, CS%nk
        do k=1,min(nz,CS%nk)+1
          write(0,*) k,zOld(k),zNew(k)
        enddo
        do k=min(nz,CS%nk)+2,CS%nk+1
          write(0,*) k,zOld(nz+1),zNew(k)
        enddo
        do k=1,min(nz,CS%nk)
          write(0,*) k,h(i,j,k),zNew(k)-zNew(k+1),CS%coordinateResolution(k)
        enddo
        do k=min(nz,CS%nk)+1,CS%nk
          write(0,*) k, 0.0, zNew(k)-zNew(k+1), CS%coordinateResolution(k)
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
subroutine build_sigma_grid( CS, G, GV, h, nom_depth_H, dzInterface )
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
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                 !! relative to mean sea level or another locally
                                                                 !! valid reference height, converted to thickness
                                                                 !! units [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),  intent(inout) :: dzInterface !< The change in interface depth
                                                                 !! [H ~> m or kg m-2]

  ! Local variables
  real :: nominalDepth    ! The nominal depth of the sea-floor in thickness units [H ~> m or kg m-2]
  real :: totalThickness  ! The total thickness of the water column [H ~> m or kg m-2]
#ifdef __DO_SAFETY_CHECKS__
  real :: dh  ! The larger of the total column thickness or bathymetric depth [H ~> m or kg m-2]
#endif
  real, dimension(SZK_(GV)+1) :: zOld    ! Previous coordinate interface heights [H ~> m or kg m-2]
  real, dimension(CS%nk+1)    :: zNew    ! New coordinate interface heights [H ~> m or kg m-2]
  integer :: i, j, k, nz

  nz = GV%ke

  do i = G%isc-1,G%iec+1
    do j = G%jsc-1,G%jec+1

      if (G%mask2dT(i,j)==0.) then
        dzInterface(i,j,:) = 0.
        cycle
      endif

      ! Determine water column height
      totalThickness = 0.0
      do k = 1,nz
        totalThickness = totalThickness + h(i,j,k)
      enddo

      ! In sigma coordinates, the bathymetric depth is only used as an arbitrary offset that
      ! cancels out when determining coordinate motion, so referencing the column postions to
      ! the surface is perfectly acceptable, but for preservation of previous answers the
      ! referencing is done relative to the bottom when in Boussinesq or semi-Boussinesq mode.
      if (GV%Boussinesq .or. GV%semi_Boussinesq) then
        nominalDepth = nom_depth_H(i,j)
      else
        nominalDepth = totalThickness
      endif

      call build_sigma_column(CS%sigma_CS, nominalDepth, totalThickness, zNew)

      ! Calculate the final change in grid position after blending new and old grids
      zOld(nz+1) =  -nominalDepth
      do k = nz,1,-1
        zOld(k) = zOld(k+1) + h(i,j,k)
      enddo

      call filtered_grid_motion( CS, nz, zOld, zNew, dzInterface(i,j,:) )

#ifdef __DO_SAFETY_CHECKS__
      dh = max(nominalDepth,totalThickness)
      if (abs(zNew(1)-zOld(1)) > (CS%nk-1)*0.5*epsilon(dh)*dh) then
        write(0,*) 'min_thickness=',CS%min_thickness
        write(0,*) 'nominalDepth=',nominalDepth,'totalThickness=',totalThickness
        write(0,*) 'dzInterface(1) = ',dzInterface(i,j,1),epsilon(dh),nz,CS%nk
        do k=1,min(nz,CS%nk)+1
          write(0,*) k,zOld(k),zNew(k)
        enddo
        do k=min(nz,CS%nk)+2,CS%nk+1
          write(0,*) k,zOld(nz+1),zNew(k)
        enddo
        do k=1,min(nz,CS%nk)
          write(0,*) k,h(i,j,k),zNew(k)-zNew(k+1),totalThickness*CS%coordinateResolution(k),CS%coordinateResolution(k)
        enddo
        do k=min(nz,CS%nk)+1,CS%nk
          write(0,*) k,0.0,zNew(k)-zNew(k+1),totalThickness*CS%coordinateResolution(k),CS%coordinateResolution(k)
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
subroutine build_rho_grid( G, GV, US, h, nom_depth_H, tv, dzInterface, remapCS, CS, frac_shelf_h )
!------------------------------------------------------------------------------
! This routine builds a new grid based on a given set of target interface
! densities (these target densities are computed by taking the mean value
! of given layer densities). The algorithm operates as follows within each
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
  type(regridding_CS),                        intent(in)    :: CS !< Regridding control structure
  type(ocean_grid_type),                      intent(in)    :: G  !< Ocean grid structure
  type(verticalGrid_type),                    intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),                      intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),           intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                  !! relative to mean sea level or another locally
                                                                  !! valid reference height, converted to thickness
                                                                  !! units [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)    :: tv !< Thermodynamics structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),   intent(inout) :: dzInterface !< The change in interface depth
                                                                  !! [H ~> m or kg m-2]
  type(remapping_CS),                         intent(in)    :: remapCS !< The remapping control structure
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)    :: frac_shelf_h  !< Fractional ice
                                                                  !! shelf coverage [nondim]
  ! Local variables
  integer :: nz  ! The number of layers in the input grid
  integer :: i, j, k
  real    :: nominalDepth   ! Depth of the bottom of the ocean, positive downward [H ~> m or kg m-2]
  real, dimension(SZK_(GV)+1) :: zOld    ! Previous coordinate interface heights [H ~> m or kg m-2]
  real, dimension(CS%nk+1)    :: zNew    ! New coordinate interface heights [H ~> m or kg m-2]
  real :: h_neglect, h_neglect_edge ! Negligible thicknesses [H ~> m or kg m-2]
  real    :: totalThickness ! Total thicknesses [H ~> m or kg m-2]
#ifdef __DO_SAFETY_CHECKS__
  real    :: dh    ! The larger of the total column thickness or bathymetric depth [H ~> m or kg m-2]
#endif
  logical :: ice_shelf

  h_neglect = set_h_neglect(GV, CS%remap_answer_date, h_neglect_edge)

  nz = GV%ke
  ice_shelf = present(frac_shelf_h)

  if (.not.CS%target_density_set) call MOM_error(FATAL, "build_rho_grid: "//&
        "Target densities must be set before build_rho_grid is called.")

  ! Build grid based on target interface densities
  do j = G%jsc-1,G%jec+1
    do i = G%isc-1,G%iec+1

      if (G%mask2dT(i,j)==0.) then
        dzInterface(i,j,:) = 0.
        cycle
      endif

      ! Determine total water column thickness
      totalThickness = 0.0
      do k=1,nz
        totalThickness = totalThickness + h(i,j,k)
      enddo

      ! In rho coordinates, the bathymetric depth is only used as an arbitrary offset that
      ! cancels out when determining coordinate motion, so referencing the column postions to
      ! the surface is perfectly acceptable, but for preservation of previous answers the
      ! referencing is done relative to the bottom when in Boussinesq or semi-Boussinesq mode.
      if (GV%Boussinesq .or. GV%semi_Boussinesq) then
        nominalDepth = nom_depth_H(i,j)
      else
        nominalDepth = totalThickness
      endif

      ! Determine absolute interface positions
      zOld(nz+1) = - nominalDepth
      do k = nz,1,-1
        zOld(k) = zOld(k+1) + h(i,j,k)
      enddo

      if (ice_shelf) then
         call build_rho_column(CS%rho_CS, nz, nominalDepth, h(i,j,:), &
              tv%T(i,j,:), tv%S(i,j,:), tv%eqn_of_state, zNew, &
              z_rigid_top=totalThickness - nominalDepth, eta_orig = zOld(1), &
              h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
      else
         call build_rho_column(CS%rho_CS, nz, nominalDepth, h(i,j,:), &
                            tv%T(i,j,:), tv%S(i,j,:), tv%eqn_of_state, zNew, &
                            h_neglect=h_neglect, h_neglect_edge=h_neglect_edge)
      endif

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
      do k=2,CS%nk
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

      dh = max(nominalDepth, totalThickness)
      if (abs(zNew(1)-zOld(1)) > (nz-1)*0.5*epsilon(dh)*dh) then
        write(0,*) 'min_thickness=',CS%min_thickness
        write(0,*) 'nominalDepth=',nominalDepth,'totalThickness=',totalThickness
        write(0,*) 'zNew(1)-zOld(1) = ',zNew(1)-zOld(1),epsilon(dh),nz
        do k=1,min(nz,CS%nk)+1
          write(0,*) k,zOld(k),zNew(k)
        enddo
        do k=min(nz,CS%nk)+2,CS%nk+1
          write(0,*) k,zOld(nz+1),zNew(k)
        enddo
        do k=1,nz
          write(0,*) k,h(i,j,k),zNew(k)-zNew(k+1)
        enddo
        do k=min(nz,CS%nk)+1,CS%nk
          write(0,*) k, 0.0, zNew(k)-zNew(k+1), CS%coordinateResolution(k)
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
!! \remark { Based on Bleck, 2002: An ocean-ice general circulation model framed in
!! hybrid isopycnic-Cartesian coordinates, Ocean Modelling 37, 55-88.
!! http://dx.doi.org/10.1016/S1463-5003(01)00012-9 }
subroutine build_grid_HyCOM1( G, GV, US, h, nom_depth_H, tv, h_new, dzInterface, remapCS, CS, frac_shelf_h, zScale )
  type(ocean_grid_type),                     intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),                     intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Existing model thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),          intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                 !! relative to mean sea level or another locally
                                                                 !! valid reference height, converted to thickness
                                                                 !! units [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Thermodynamics structure
  type(remapping_CS),                        intent(in)    :: remapCS !< The remapping control structure
  type(regridding_CS),                       intent(in)    :: CS !< Regridding control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk),    intent(inout) :: h_new !< New layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),  intent(inout) :: dzInterface !< Changes in interface position
                                                                 !! in thickness units [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)   :: frac_shelf_h !< Fractional ice shelf
                                                                 !! coverage [nondim]
  real,                            optional, intent(in)    :: zScale !< Scaling factor from the target coordinate
                                                                 !! resolution in Z to desired units for zInterface,
                                                                 !! usually Z_to_H in which case it is in
                                                                 !! units of [H Z-1 ~> nondim or kg m-3]

  ! Local variables
  real, dimension(SZK_(GV)+1) :: z_col  ! Source interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(SZK_(GV))   :: p_col  ! Layer center pressure in the input column [R L2 T-2 ~> Pa]
  real, dimension(CS%nk+1) :: z_col_new ! New interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(CS%nk+1) :: dz_col    ! The realized change in z_col [H ~> m or kg m-2]
  real :: nominalDepth    ! The nominal depth of the seafloor in thickness units [H ~> m or kg m-2]
  real :: h_neglect, h_neglect_edge ! Negligible thicknesses used for remapping [H ~> m or kg m-2]
  real :: z_top_col       ! The nominal height of the sea surface or ice-ocean interface
                          ! in thickness units [H ~> m or kg m-2]
  real :: totalThickness  ! The total thickness of the water column [H ~> m or kg m-2]
  logical :: ice_shelf
  integer :: i, j, k, nki

  h_neglect = set_h_neglect(GV, CS%remap_answer_date, h_neglect_edge)

  if (.not.CS%target_density_set) call MOM_error(FATAL, "build_grid_HyCOM1 : "//&
        "Target densities must be set before build_grid_HyCOM1 is called.")

  nki = min(GV%ke, CS%nk)
  ice_shelf = present(frac_shelf_h)

  ! Build grid based on target interface densities
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1
    if (G%mask2dT(i,j)>0.) then

      nominalDepth = nom_depth_H(i,j)

      if (ice_shelf) then
        totalThickness = 0.0
        do k=1,GV%ke
          totalThickness = totalThickness + h(i,j,k)
        enddo
        z_top_col = max(nominalDepth-totalThickness,0.0)
      else
        z_top_col = 0.0
      endif

      z_col(1) = z_top_col ! Work downward rather than bottom up
      do K = 1, GV%ke
        z_col(K+1) = z_col(K) + h(i,j,k)
        p_col(k) = tv%P_Ref + CS%compressibility_fraction * &
             ( 0.5 * ( z_col(K) + z_col(K+1) ) * (GV%H_to_RZ*GV%g_Earth) - tv%P_Ref )
      enddo

      call build_hycom1_column(CS%hycom_CS, remapCS, tv%eqn_of_state, GV%ke, nominalDepth, &
           h(i,j,:), tv%T(i,j,:), tv%S(i,j,:), p_col, &
           z_col, z_col_new, zScale=zScale, &
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
subroutine build_grid_adaptive(G, GV, US, h, nom_depth_H, tv, dzInterface, remapCS, CS)
  type(ocean_grid_type),                       intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),                       intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),   intent(in)    :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),            intent(in)    :: nom_depth_H !< The bathymetric depth of this column
                                                                     !! relative to mean sea level or another locally
                                                                     !! valid reference height, converted to thickness
                                                                     !! units [H ~> m or kg m-2]
  type(thermo_var_ptrs),                       intent(in)    :: tv   !< A structure pointing to various
                                                                     !! thermodynamic variables
  type(regridding_CS),                         intent(in)    :: CS   !< Regridding control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1),    intent(inout) :: dzInterface !< The change in interface depth
                                                                     !! [H ~> m or kg m-2]
  type(remapping_CS),                          intent(in)    :: remapCS !< The remapping control structure

  ! local variables
  integer :: i, j, k, nz ! indices and dimension lengths
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: tInt ! Temperature on interfaces [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: sInt ! Salinity on interfaces [S ~> ppt]
  ! current interface positions and after tendency term is applied
  ! positive downward
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: zInt ! Interface depths [H ~> m or kg m-2]
  real, dimension(SZK_(GV)+1) :: zNext  ! New interface depths [H ~> m or kg m-2]

  nz = GV%ke

  call assert((GV%ke == CS%nk), "build_grid_adaptive is only written to work "//&
                                "with the same number of input and target layers.")

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

    call build_adapt_column(CS%adapt_CS, G, GV, US, tv, i, j, zInt, tInt, sInt, h, nom_depth_H, zNext)

    call filtered_grid_motion(CS, nz, zInt(i,j,:), zNext, dzInterface(i,j,:))
    ! convert from depth to z
    do K = 1, nz+1 ; dzInterface(i,j,K) = -dzInterface(i,j,K) ; enddo
    call adjust_interface_motion(CS, nz, h(i,j,:), dzInterface(i,j,:))
  enddo ; enddo
end subroutine build_grid_adaptive

!> Adjust dz_Interface to ensure non-negative future thicknesses
subroutine adjust_interface_motion( CS, nk, h_old, dz_int )
  type(regridding_CS),      intent(in)    :: CS !< Regridding control structure
  integer,                  intent(in)    :: nk !< Number of layers in h_old
  real, dimension(nk),      intent(in)    :: h_old  !< Layer thicknesses on the old grid [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: dz_int !< Interface movements, adjusted to keep the thicknesses
                                                    !! thicker than their minimum value [H ~> m or kg m-2]
  ! Local variables
  real :: h_new   ! A layer thickness on the new grid [H ~> m or kg m-2]
  real :: eps     ! A tiny relative thickness [nondim]
  real :: h_total ! The total thickness of the old grid [H ~> m or kg m-2]
  real :: h_err   ! An error tolerance that use used to flag unacceptably large negative layer thicknesses
                  ! that can not be explained by roundoff errors [H ~> m or kg m-2]
  integer :: k

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
!> make sure all layers are at least as thick as the minimum thickness allowed
!! for regridding purposes by inflating thin layers.  This breaks mass conservation
!! and adds mass to the model when there are excessively thin layers.
subroutine inflate_vanished_layers_old( CS, G, GV, h )
!------------------------------------------------------------------------------
! This routine is called when initializing the regridding options. The
! objective is to make sure all layers are at least as thick as the minimum
! thickness allowed for regridding purposes (this parameter is set in the
! MOM_input file or defaulted to 1.0e-3). When layers are too thin, they
! are inflated up to the minimum thickness.
!------------------------------------------------------------------------------

  ! Arguments
  type(regridding_CS),                       intent(in)    :: CS   !< Regridding control structure
  type(ocean_grid_type),                     intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2]

  ! Local variables
  integer :: i, j, k
  real    :: hTmp(GV%ke) ! A copy of a 1-d column of h [H ~> m or kg m-2]

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
  real      :: T0, T1       ! temperatures of two layers [C ~> degC]
  real      :: S0, S1       ! salinities of two layers [S ~> ppt]
  real      :: r0, r1       ! densities of two layers [R ~> kg m-3]
  real      :: h0, h1       ! Layer thicknesses  [H ~> m or kg m-2]
  real, dimension(GV%ke) :: p_col  ! A column of zero pressures [R L2 T-2 ~> Pa]
  real, dimension(GV%ke) :: densities ! Densities in the column [R ~> kg m-3]
  logical   :: stratified
  integer   :: i, j, k

  !### Doing convective adjustment based on potential densities with zero pressure seems
  !    questionable, although it does avoid ambiguous sorting. -RWH
  p_col(:) = 0.

  ! Loop on columns
  do j = G%jsc-1,G%jec+1 ; do i = G%isc-1,G%iec+1

    ! Compute densities within current water column
    call calculate_density(tv%T(i,j,:), tv%S(i,j,:), p_col, densities, tv%eqn_of_state)

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
          call calculate_density(tv%T(i,j,k), tv%S(i,j,k), p_col(k), densities(k), tv%eqn_of_state)
          call calculate_density(tv%T(i,j,k+1), tv%S(i,j,k+1), p_col(k+1), &
                                  densities(k+1), tv%eqn_of_state )
          ! Because p_col is has uniform values, these calculate_density calls are equivalent to
          ! densities(k) = r1 ; densities(k+1) = r0
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
                                            !! See the documentation for regrid_consts
                                            !! for the recognized values.
  real,             intent(in) :: maxDepth  !< The range of the grid values in some modes, in coordinate
                                            !! dependent units that might be [m] or [kg m-3] or [nondim]
                                            !! or something else.
  real,             intent(in) :: rhoLight  !< The minimum value of the grid in RHO mode [kg m-3]
  real,             intent(in) :: rhoHeavy  !< The maximum value of the grid in RHO mode [kg m-3]

  real                         :: uniformResolution(nk) !< The returned uniform resolution grid, in
                                            !! coordinate dependent units that might be [m] or
                                            !! [kg m-3] or [nondim] or something else.

  ! Local variables
  integer :: scheme

  scheme = coordinateMode(coordMode)
  select case ( scheme )

    case ( REGRIDDING_ZSTAR, REGRIDDING_HYCOM1, REGRIDDING_HYBGEN, &
           REGRIDDING_SIGMA_SHELF_ZSTAR, REGRIDDING_ADAPTIVE )
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
subroutine initCoord(CS, GV, US, coord_mode, param_file)
  type(regridding_CS),     intent(inout) :: CS !< Regridding control structure
  character(len=*),        intent(in)    :: coord_mode !< A string indicating the coordinate mode.
                                               !! See the documentation for regrid_consts
                                               !! for the recognized values.
  type(verticalGrid_type), intent(in)    :: GV !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file

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
  case (REGRIDDING_HYBGEN)
    call init_hybgen_regrid(CS%hybgen_CS, GV, US, param_file)
  case (REGRIDDING_ADAPTIVE)
    call init_coord_adapt(CS%adapt_CS, CS%nk, CS%coordinateResolution, GV%m_to_H, US%kg_m3_to_R)
  end select
end subroutine initCoord

!------------------------------------------------------------------------------
!> Set the fixed resolution data
subroutine setCoordinateResolution( dz, CS, scale )
  real, dimension(:),  intent(in)    :: dz !< A vector of vertical grid spacings, in arbitrary coordinate
                                           !! dependent units, such as [m] for a z-coordinate or [kg m-3]
                                           !! for a density coordinate.
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure
  real,      optional, intent(in)    :: scale !< A scaling factor converting dz to the internal represetation
                                           !! of coordRes, in various units that depend on the coordinate,
                                           !! such as [Z m-1 ~> 1] for a z-coordinate or [R m3 kg-1 ~> 1] for
                                           !! a density coordinate.

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
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(regridding_CS),     intent(inout) :: CS !< Regridding control structure
  ! Local variables
  integer :: k, nz

  nz = CS%nk
  if (nz == 1) then ! Set a broad range of bounds.  Regridding may not be meaningful in this case.
    CS%target_density(1)    = 0.0
    CS%target_density(2)    = 2.0*GV%Rlay(1)
  else
    CS%target_density(1)    = (GV%Rlay(1) + 0.5*(GV%Rlay(1)-GV%Rlay(2)))
    CS%target_density(nz+1) = (GV%Rlay(nz) + 0.5*(GV%Rlay(nz)-GV%Rlay(nz-1)))
    do k=2,nz
      CS%target_density(k)  = CS%target_density(k-1) + CS%coordinateResolution(k)
    enddo
  endif
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
  real, dimension(CS%nk+1), intent(in)    :: max_depths !< Maximum interface depths, in arbitrary units, often [m]
  real, optional,           intent(in)    :: units_to_H !< A conversion factor for max_depths into H units,
                                                        !! often in [H m-1 ~> 1 or kg m-3]
  ! Local variables
  real :: val_to_H ! A conversion factor from the units for max_depths into H units, often [H m-1 ~> 1 or kg m-3]
                   ! if units_to_H is present, or [nondim] if it is absent.
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
  end select
end subroutine set_regrid_max_depths

!> Set maximum layer thicknesses based on a vector of input values.
subroutine set_regrid_max_thickness( CS, max_h, units_to_H )
  type(regridding_CS),      intent(inout) :: CS !< Regridding control structure
  real, dimension(CS%nk+1), intent(in)    :: max_h !< Maximum layer thicknesses, in arbitrary units, often [m]
  real, optional,           intent(in)    :: units_to_H !< A conversion factor for max_h into H units,
                                                        !! often [H m-1 ~> 1 or kg m-3]
  ! Local variables
  real :: val_to_H ! A conversion factor from the units for max_h into H units, often [H m-1 ~> 1 or kg m-3]
                   ! if units_to_H is present, or [nondim] if it is absent.
  integer :: k

  if (.not.allocated(CS%max_layer_thickness)) allocate(CS%max_layer_thickness(1:CS%nk))

  val_to_H = 1.0 ; if (present( units_to_H)) val_to_H = units_to_H

  do k=1,CS%nk
    CS%max_layer_thickness(k) = val_to_H * max_h(k)
  enddo

  ! set max thickness for coordinate
  select case (CS%regridding_scheme)
  case (REGRIDDING_HYCOM1)
    call set_hycom_params(CS%hycom_CS, max_layer_thickness=CS%max_layer_thickness)
  end select
end subroutine set_regrid_max_thickness


!> Write the vertical coordinate information into a file.
!! This subroutine writes out a file containing any available data related
!! to the vertical grid used by the MOM ocean model when in ALE mode.
subroutine write_regrid_file( CS, GV, filepath )
  type(regridding_CS),     intent(in) :: CS        !< Regridding control structure
  type(verticalGrid_type), intent(in) :: GV        !< ocean vertical grid structure
  character(len=*),        intent(in) :: filepath  !< The full path to the file to write

  type(vardesc)      :: vars(2)
  type(MOM_field)    :: fields(2)
  type(MOM_netCDF_file) :: IO_handle ! The I/O handle of the fileset
  real               :: ds(GV%ke), dsi(GV%ke+1)  ! The labeling layer and interface coordinates for output
                                                 ! in axes in files, in coordinate-dependent units that can
                                                 ! be obtained from getCoordinateUnits [various]

  if (CS%regridding_scheme == REGRIDDING_HYBGEN) then
    call write_Hybgen_coord_file(GV, CS%hybgen_CS, filepath)
    return
  endif

  ds(:)       = CS%coord_scale * CS%coordinateResolution(:)
  dsi(1)      = 0.5*ds(1)
  dsi(2:GV%ke) = 0.5*( ds(1:GV%ke-1) + ds(2:GV%ke) )
  dsi(GV%ke+1) = 0.5*ds(GV%ke)

  vars(1) = var_desc('ds', getCoordinateUnits( CS ), &
                     'Layer Coordinate Thickness', '1', 'L', '1')
  vars(2) = var_desc('ds_interface', getCoordinateUnits( CS ), &
                     'Layer Center Coordinate Separation', '1', 'i', '1')

  call create_MOM_file(IO_handle, trim(filepath), vars, 2, fields, &
      SINGLE_FILE, GV=GV)
  call MOM_write_field(IO_handle, fields(1), ds)
  call MOM_write_field(IO_handle, fields(2), dsi)
  call IO_handle%close()

end subroutine write_regrid_file

!> Set appropriate values for the negligible thicknesses used for remapping based on an answer date.
function set_h_neglect(GV, remap_answer_date, h_neglect_edge) result(h_neglect)
  type(verticalGrid_type), intent(in)  :: GV   !< Ocean vertical grid structure
  integer,                 intent(in)  :: remap_answer_date !< The vintage of the expressions to use
                                               !! for remapping.  Values below 20190101 recover the
                                               !! remapping answers from 2018. Higher values use more
                                               !! robust forms of the same remapping algorithms.
  real,                    intent(out) :: h_neglect_edge !< A negligibly small thickness used in
                                               !! remapping edge value calculations [H ~> m or kg m-2]
  real                                 :: h_neglect !< A negligibly small thickness used in
                                               !! remapping cell reconstructions [H ~> m or kg m-2]

  if (remap_answer_date >= 20190101) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif
end function set_h_neglect

!> Set appropriate values for the negligible vertical layer extents used for remapping based on an answer date.
function set_dz_neglect(GV, US, remap_answer_date, dz_neglect_edge) result(dz_neglect)
  type(verticalGrid_type), intent(in)  :: GV   !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  integer,                 intent(in)  :: remap_answer_date !< The vintage of the expressions to use
                                               !! for remapping.  Values below 20190101 recover the
                                               !! remapping answers from 2018. Higher values use more
                                               !! robust forms of the same remapping algorithms.
  real,                    intent(out) :: dz_neglect_edge !< A negligibly small vertical layer extent
                                               !! used in remapping edge value calculations [Z ~> m]
  real                                 :: dz_neglect !< A negligibly small vertical layer extent
                                               !! used in remapping cell reconstructions [Z ~> m]

  if (remap_answer_date >= 20190101) then
    dz_neglect = GV%dZ_subroundoff ; dz_neglect_edge = GV%dZ_subroundoff
  elseif (GV%Boussinesq) then
    dz_neglect = US%m_to_Z*1.0e-30 ; dz_neglect_edge = US%m_to_Z*1.0e-10
  else
    dz_neglect = GV%kg_m2_to_H * (GV%H_to_m*US%m_to_Z) * 1.0e-30
    dz_neglect_edge = GV%kg_m2_to_H * (GV%H_to_m*US%m_to_Z) * 1.0e-10
  endif
end function set_dz_neglect

!------------------------------------------------------------------------------
!> Query the fixed resolution data
function getCoordinateResolution( CS, undo_scaling )
  type(regridding_CS), intent(in) :: CS !< Regridding control structure
  logical,   optional, intent(in) :: undo_scaling !< If present and true, undo any internal
                                        !! rescaling of the resolution data.
  real, dimension(CS%nk)          :: getCoordinateResolution !< The resolution or delta of the target coordinate,
                                                             !! in units that depend on the coordinate [various]

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
  real, dimension(CS%nk+1)        :: getCoordinateInterfaces !< Interface positions in target coordinate,
                                                             !! in units that depend on the coordinate [various]

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
    case ( REGRIDDING_ZSTAR, REGRIDDING_HYCOM1, REGRIDDING_HYBGEN, &
           REGRIDDING_ADAPTIVE )
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
    case ( REGRIDDING_HYBGEN )
      getCoordinateShortName = 'hybrid'
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
             compress_fraction, ref_pressure, &
             integrate_downward_for_e, remap_answers_2018, remap_answer_date, regrid_answer_date, &
             adaptTimeRatio, adaptZoom, adaptZoomCoeff, adaptBuoyCoeff, adaptAlpha, adaptDoMin, adaptDrho0)
  type(regridding_CS), intent(inout) :: CS !< Regridding control structure
  logical, optional, intent(in) :: boundary_extrapolation !< Extrapolate in boundary cells
  real,    optional, intent(in) :: min_thickness    !< Minimum thickness allowed when building the
                                                    !! new grid [H ~> m or kg m-2]
  real,    optional, intent(in) :: old_grid_weight  !< Weight given to old coordinate when time-filtering grid [nondim]
  character(len=*), optional, intent(in) :: interp_scheme !< Interpolation method for state-dependent coordinates
  real,    optional, intent(in) :: depth_of_time_filter_shallow !< Depth to start cubic [H ~> m or kg m-2]
  real,    optional, intent(in) :: depth_of_time_filter_deep !< Depth to end cubic [H ~> m or kg m-2]
  real,    optional, intent(in) :: compress_fraction !< Fraction of compressibility to add to potential density [nondim]
  real,    optional, intent(in) :: ref_pressure     !< The reference pressure for density-dependent
                                                    !! coordinates [R L2 T-2 ~> Pa]
  logical, optional, intent(in) :: integrate_downward_for_e !< If true, integrate for interface positions downward
                                                    !! from the top.
  logical, optional, intent(in) :: remap_answers_2018 !< If true, use the order of arithmetic and expressions
                                                    !! that recover the remapping answers from 2018.  Otherwise
                                                    !! use more robust but mathematically equivalent expressions.
  integer, optional, intent(in) :: remap_answer_date !< The vintage of the expressions to use for remapping
  integer, optional, intent(in) :: regrid_answer_date !< The vintage of the expressions to use for regridding
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
  if (present(regrid_answer_date)) call set_interp_answer_date(CS%interp_CS, regrid_answer_date)

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
  if (present(remap_answers_2018)) then
    if (remap_answers_2018) then
      CS%remap_answer_date = 20181231
    else
      CS%remap_answer_date = 20190101
    endif
  endif
  if (present(remap_answer_date)) CS%remap_answer_date = remap_answer_date

  select case (CS%regridding_scheme)
  case (REGRIDDING_ZSTAR)
    if (present(min_thickness)) call set_zlike_params(CS%zlike_CS, min_thickness=min_thickness)
  case (REGRIDDING_SIGMA_SHELF_ZSTAR)
    if (present(min_thickness)) call set_zlike_params(CS%zlike_CS, min_thickness=min_thickness)
  case (REGRIDDING_SIGMA)
    if (present(min_thickness)) call set_sigma_params(CS%sigma_CS, min_thickness=min_thickness)
  case (REGRIDDING_RHO)
    if (present(min_thickness)) call set_rho_params(CS%rho_CS, min_thickness=min_thickness)
    if (present(ref_pressure)) call set_rho_params(CS%rho_CS, ref_pressure=ref_pressure)
    if (present(integrate_downward_for_e)) &
      call set_rho_params(CS%rho_CS, integrate_downward_for_e=integrate_downward_for_e)
    if (associated(CS%rho_CS) .and. (present(interp_scheme) .or. present(boundary_extrapolation))) &
      call set_rho_params(CS%rho_CS, interp_CS=CS%interp_CS)
  case (REGRIDDING_HYCOM1)
    if (associated(CS%hycom_CS) .and. (present(interp_scheme) .or. present(boundary_extrapolation))) &
      call set_hycom_params(CS%hycom_CS, interp_CS=CS%interp_CS)
  case (REGRIDDING_HYBGEN)
    ! Do nothing for now.
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
  real,                intent(in) :: SSH   !< The sea surface height, in the same units as depth, often [Z ~> m]
  real,                intent(in) :: depth !< The maximum depth of the grid, often [Z ~> m]
  real, dimension(CS%nk)          :: getStaticThickness !< The returned thicknesses in the units of
                                           !! depth, often [Z ~> m]
  ! Local
  integer :: k
  real :: z, dz  ! Vertical positions and grid spacing [Z ~> m]

  select case ( CS%regridding_scheme )
    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA_SHELF_ZSTAR, REGRIDDING_HYCOM1, REGRIDDING_HYBGEN, &
           REGRIDDING_ADAPTIVE )
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
  real, dimension(:), intent(inout) :: dz     !< Profile of nominal thicknesses [m] or other units
  ! Local variables
  integer :: nk, k
  real    :: dz_min  ! minimum grid spacing [m] or other units
  real    :: power   ! A power to raise the relative position in index space [nondim]
  real    :: prec    ! The precision with which positions are returned [m] or other units
  real    :: H_total ! The sum of the nominal thicknesses [m] or other units

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

!> Construct the name of a parameter for a specific coordinate based on param_prefix and param_suffix. For the main,
!! prognostic coordinate this will simply return the parameter name (e.g. P_REF)
function create_coord_param(param_prefix, param_name, param_suffix) result(coord_param)
  character(len=*) :: param_name   !< The base name of the parameter (e.g. the one used for the main coordinate)
  character(len=*) :: param_prefix !< String to prefix to parameter names.
  character(len=*) :: param_suffix !< String to append to parameter names.
  character(len=MAX_PARAM_LENGTH) :: coord_param  !< Parameter name prepended by param_prefix
                                                  !! and appended with param_suffix
  integer :: out_length

  if (len_trim(param_prefix) + len_trim(param_suffix) == 0) then
    coord_param = param_name
  else
    ! Note the +2 is because of two underscores
    out_length = len_trim(param_name)+len_trim(param_prefix)+len_trim(param_suffix)+2
    if (out_length > MAX_PARAM_LENGTH) then
      call MOM_error(FATAL,"Coordinate parameter is too long; increase MAX_PARAM_LENGTH")
    endif
    coord_param = TRIM(param_prefix)//"_"//TRIM(param_name)//"_"//TRIM(param_suffix)
  endif

end function create_coord_param

!> Parses a string and generates a rho_target(:) profile with refined resolution downward
!! and returns the number of levels
integer function rho_function1( string, rho_target )
  character(len=*),   intent(in)    :: string !< String with list of parameters in form
                                              !! dz_min, H_total, power, precision
  real, dimension(:), allocatable, intent(inout) :: rho_target !< Profile of interface densities [kg m-3]
  ! Local variables
  integer :: nki, k, nk
  real    :: dx   ! Fractional distance from interface nki [nondim]
  real    :: ddx  ! Change in dx between interfaces [nondim]
  real    :: rho_1, rho_2 ! Density of the top two layers in a profile [kg m-3]
  real    :: rho_3    ! Density in the third layer, below which the density increase linearly
                      ! in subsequent layers [kg m-3]
  real    :: drho     ! Change in density over the linear region [kg m-3]
  real    :: rho_4    ! The densest density in this profile [kg m-3], which might be very large.
  real    :: drho_min ! A minimal fractional density difference [nondim]?

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
