!> Code that initializes fixed aspects of the model grid, such as horizontal
!! grid metrics, topography and Coriolis, and can be shared between components.
module MOM_shared_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms, only : max_across_PEs, reproducing_sum
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE, AGRID
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, log_param, param_file_type, log_version
use MOM_io, only : close_file, create_file, file_type, fieldtype, file_exists, field_size
use MOM_io, only : MOM_read_data, MOM_read_vector, read_variable, stdout
use MOM_io, only : open_file_to_read, close_file_to_read, SINGLE_FILE, MULTIPLE
use MOM_io, only : slasher, vardesc, MOM_write_field, var_desc
use MOM_string_functions, only : uppercase
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

public MOM_shared_init_init
public MOM_initialize_rotation, MOM_calculate_grad_Coriolis
public initialize_topography_from_file, apply_topography_edits_from_file
public initialize_topography_named, limit_topography, diagnoseMaximumDepth
public set_rotation_planetary, set_rotation_beta_plane, initialize_grid_rotation_angle
public reset_face_lengths_named, reset_face_lengths_file, reset_face_lengths_list
public read_face_length_list, set_velocity_depth_max, set_velocity_depth_min
public compute_global_grid_integrals, write_ocean_geometry_file

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

! -----------------------------------------------------------------------------
!> MOM_shared_init_init just writes the code version.
subroutine MOM_shared_init_init(PF)
  type(param_file_type),   intent(in)    :: PF   !< A structure indicating the open file
                                                 !! to parse for model parameter values.

  character(len=40)  :: mdl = "MOM_shared_initialization" ! This module's name.

! This include declares and sets the variable "version".
#include "version_variable.h"
  call log_version(PF, mdl, version, &
   "Sharable code to initialize time-invariant fields, like bathymetry and Coriolis parameters.")

end subroutine MOM_shared_init_init
! -----------------------------------------------------------------------------

!> MOM_initialize_rotation makes the appropriate call to set up the Coriolis parameter.
subroutine MOM_initialize_rotation(f, G, PF, US)
  type(dyn_horgrid_type),                       intent(in)  :: G  !< The dynamic horizontal grid type
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), intent(out) :: f  !< The Coriolis parameter [T-1 ~> s-1]
  type(param_file_type),                        intent(in)  :: PF !< Parameter file structure
  type(unit_scale_type),              optional, intent(in)  :: US !< A dimensional unit scaling type

!   This subroutine makes the appropriate call to set up the Coriolis parameter.
! This is a separate subroutine so that it can be made public and shared with
! the ice-sheet code or other components.
! Set up the Coriolis parameter, f, either analytically or from file.
  character(len=40)  :: mdl = "MOM_initialize_rotation" ! This subroutine's name.
  character(len=200) :: config

  call callTree_enter(trim(mdl)//"(), MOM_shared_initialization.F90")
  call get_param(PF, mdl, "ROTATION", config, &
                 "This specifies how the Coriolis parameter is specified: \n"//&
                 " \t 2omegasinlat - Use twice the planetary rotation rate \n"//&
                 " \t\t times the sine of latitude.\n"//&
                 " \t betaplane - Use a beta-plane or f-plane.\n"//&
                 " \t USER - call a user modified routine.", &
                 default="2omegasinlat")
  select case (trim(config))
    case ("2omegasinlat"); call set_rotation_planetary(f, G, PF, US)
    case ("beta"); call set_rotation_beta_plane(f, G, PF, US)
    case ("betaplane"); call set_rotation_beta_plane(f, G, PF, US)
   !case ("nonrotating") ! Note from AJA: Missing case?
    case default ; call MOM_error(FATAL,"MOM_initialize: "// &
      "Unrecognized rotation setup "//trim(config))
  end select
  call callTree_leave(trim(mdl)//'()')
end subroutine MOM_initialize_rotation

!> Calculates the components of grad f (Coriolis parameter)
subroutine MOM_calculate_grad_Coriolis(dF_dx, dF_dy, G, US)
  type(dyn_horgrid_type),             intent(inout) :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out)   :: dF_dx !< x-component of grad f [T-1 L-1 ~> s-1 m-1]
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out)   :: dF_dy !< y-component of grad f [T-1 L-1 ~> s-1 m-1]
  type(unit_scale_type),    optional, intent(in)    :: US !< A dimensional unit scaling type
  ! Local variables
  integer :: i,j
  real :: m_to_L  ! A unit conversion factor [L m-1 ~> nondim]
  real :: f1, f2

  m_to_L = 1.0 ; if (present(US)) m_to_L = US%m_to_L

  if ((LBOUND(G%CoriolisBu,1) > G%isc-1) .or. &
      (LBOUND(G%CoriolisBu,2) > G%isc-1)) then
    ! The gradient of the Coriolis parameter can not be calculated with this grid.
    dF_dx(:,:) = 0.0 ; dF_dy(:,:) = 0.0
    return
  endif

  do j=G%jsc, G%jec ; do i=G%isc, G%iec
    f1 = 0.5*( G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1) )
    f2 = 0.5*( G%CoriolisBu(I-1,J) + G%CoriolisBu(I-1,J-1) )
    dF_dx(i,j) = G%IdxT(i,j) * ( f1 - f2 )
    f1 = 0.5*( G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J) )
    f2 = 0.5*( G%CoriolisBu(I,J-1) + G%CoriolisBu(I-1,J-1) )
    dF_dy(i,j) = G%IdyT(i,j) * ( f1 - f2 )
  enddo ; enddo
  call pass_vector(dF_dx, dF_dy, G%Domain, stagger=AGRID)

end subroutine MOM_calculate_grad_Coriolis

!> Return the global maximum ocean bottom depth in the same units as the input depth.
function diagnoseMaximumDepth(D, G)
  type(dyn_horgrid_type),  intent(in) :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(in) :: D !< Ocean bottom depth in m or Z
  real :: diagnoseMaximumDepth             !< The global maximum ocean bottom depth in m or Z
  ! Local variables
  integer :: i,j
  diagnoseMaximumDepth = D(G%isc,G%jsc)
  do j=G%jsc, G%jec ; do i=G%isc, G%iec
    diagnoseMaximumDepth = max(diagnoseMaximumDepth,D(i,j))
  enddo ; enddo
  call max_across_PEs(diagnoseMaximumDepth)
end function diagnoseMaximumDepth


!> Read gridded depths from file
subroutine initialize_topography_from_file(D, G, param_file, US)
  type(dyn_horgrid_type),           intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: D !< Ocean bottom depth in m or Z if US is present
  type(param_file_type),            intent(in)  :: param_file !< Parameter file structure
  type(unit_scale_type),  optional, intent(in)  :: US !< A dimensional unit scaling type
  ! Local variables
  real :: m_to_Z  ! A dimensional rescaling factor.
  character(len=200) :: filename, topo_file, inputdir ! Strings for file/path
  character(len=200) :: topo_varname                  ! Variable name in file
  character(len=40)  :: mdl = "initialize_topography_from_file" ! This subroutine's name.

  call callTree_enter(trim(mdl)//"(), MOM_shared_initialization.F90")

  m_to_Z = 1.0 ; if (present(US)) m_to_Z = US%m_to_Z

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "TOPO_FILE", topo_file, &
                 "The file from which the bathymetry is read.", &
                 default="topog.nc")
  call get_param(param_file, mdl, "TOPO_VARNAME", topo_varname, &
                 "The name of the bathymetry variable in TOPO_FILE.", &
                 default="depth")

  filename = trim(inputdir)//trim(topo_file)
  call log_param(param_file, mdl, "INPUTDIR/TOPO_FILE", filename)

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topography_from_file: Unable to open "//trim(filename))

  D(:,:) = -9.e30*m_to_Z ! Initializing to a very large negative depth (tall mountains) everywhere
                         ! before reading from a file should do nothing. However, in the instance of
                         ! masked-out PEs, halo regions are not updated when a processor does not
                         ! exist. We need to ensure the depth in masked-out PEs appears to be that
                         ! of land so this line does that in the halo regions. For non-masked PEs
                         ! the halo region is filled properly with a later pass_var().
  call MOM_read_data(filename, trim(topo_varname), D, G%Domain, scale=m_to_Z)

  call apply_topography_edits_from_file(D, G, param_file, US)

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_topography_from_file

!> Applies a list of topography overrides read from a netcdf file
subroutine apply_topography_edits_from_file(D, G, param_file, US)
  type(dyn_horgrid_type),           intent(in)    :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(inout) :: D !< Ocean bottom depth in m or Z if US is present
  type(param_file_type),            intent(in)    :: param_file !< Parameter file structure
  type(unit_scale_type),  optional, intent(in)    :: US !< A dimensional unit scaling type

  ! Local variables
  real :: m_to_Z  ! A dimensional rescaling factor.
  real, dimension(:), allocatable :: new_depth ! The new values of the depths [m]
  integer, dimension(:), allocatable :: ig, jg ! The global indicies of the points to modify
  character(len=200) :: topo_edits_file, inputdir ! Strings for file/path
  character(len=40)  :: mdl = "apply_topography_edits_from_file" ! This subroutine's name.
  integer :: i, j, n, ncid, n_edits, i_file, j_file, ndims, sizes(8)
  logical :: found

  call callTree_enter(trim(mdl)//"(), MOM_shared_initialization.F90")

  m_to_Z = 1.0 ; if (present(US)) m_to_Z = US%m_to_Z

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "TOPO_EDITS_FILE", topo_edits_file, &
                 "The file from which to read a list of i,j,z topography overrides.", &
                 default="")

  if (len_trim(topo_edits_file)==0) return

  topo_edits_file = trim(inputdir)//trim(topo_edits_file)
  if (is_root_PE()) then
    if (.not.file_exists(topo_edits_file, G%Domain)) &
      call MOM_error(FATAL, trim(mdl)//': Unable to find file '//trim(topo_edits_file))
    call open_file_to_read(topo_edits_file, ncid)
  else
    ncid = -1
  endif

  ! Read and check the values of ni and nj in the file for consistency with this configuration.
  call read_variable(topo_edits_file, 'ni', i_file, ncid_in=ncid)
  call read_variable(topo_edits_file, 'nj', j_file, ncid_in=ncid)
  if (i_file /= G%ieg) call MOM_error(FATAL, trim(mdl)//': Incompatible i-dimension of grid in '//&
                                      trim(topo_edits_file))
  if (j_file /= G%jeg) call MOM_error(FATAL, trim(mdl)//': Incompatible j-dimension of grid in '//&
                                      trim(topo_edits_file))

  ! Get nEdits
  call field_size(topo_edits_file, 'zEdit', sizes, ndims=ndims, ncid_in=ncid)
  if (ndims /= 1) call MOM_error(FATAL, "The variable zEdit has an "//&
            "unexpected number of dimensions in "//trim(topo_edits_file) )
  n_edits = sizes(1)
  allocate(ig(n_edits))
  allocate(jg(n_edits))
  allocate(new_depth(n_edits))

  ! Read iEdit, jEdit and zEdit
  call read_variable(topo_edits_file, 'iEdit', ig, ncid_in=ncid)
  call read_variable(topo_edits_file, 'jEdit', jg, ncid_in=ncid)
  call read_variable(topo_edits_file, 'zEdit', new_depth, ncid_in=ncid)
  call close_file_to_read(ncid, topo_edits_file)

  do n = 1, n_edits
    i = ig(n) - G%isd_global + 2 ! +1 for python indexing and +1 for ig-isd_global+1
    j = jg(n) - G%jsd_global + 2
    if (i>=G%isc .and. i<=G%iec .and. j>=G%jsc .and. j<=G%jec) then
      if (new_depth(n)/=0.) then
        write(stdout,'(a,3i5,f8.2,a,f8.2,2i4)') &
          'Ocean topography edit: ', n, ig(n), jg(n), D(i,j)/m_to_Z, '->', abs(new_depth(n)), i, j
        D(i,j) = abs(m_to_Z*new_depth(n)) ! Allows for height-file edits (i.e. converts negatives)
      else
        call MOM_error(FATAL, trim(mdl)//': A zero depth edit would change the land mask and '//&
          "is not allowed in"//trim(topo_edits_file))
      endif
    endif
  enddo

  deallocate( ig, jg, new_depth )

  call callTree_leave(trim(mdl)//'()')
end subroutine apply_topography_edits_from_file

!> initialize the bathymetry based on one of several named idealized configurations
subroutine initialize_topography_named(D, G, param_file, topog_config, max_depth, US)
  type(dyn_horgrid_type),           intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: D !< Ocean bottom depth in m or Z if US is present
  type(param_file_type),            intent(in)  :: param_file !< Parameter file structure
  character(len=*),                 intent(in)  :: topog_config !< The name of an idealized
                                                              !! topographic configuration
  real,                             intent(in)  :: max_depth  !< Maximum depth of model in the units of D
  type(unit_scale_type),  optional, intent(in)  :: US !< A dimensional unit scaling type

  ! This subroutine places the bottom depth in m into D(:,:), shaped according to the named config.

  ! Local variables
  real :: m_to_Z               ! A dimensional rescaling factor.
  real :: min_depth            ! The minimum depth [Z ~> m].
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum  basin depth MAXIMUM_DEPTH.
  real :: expdecay             ! A decay scale of associated with the sloping boundaries [m].
  real :: Dedge                ! The depth [Z ~> m], at the basin edge
! real :: south_lat, west_lon, len_lon, len_lat, Rad_earth
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  character(len=40)  :: mdl = "initialize_topography_named" ! This subroutine's name.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call callTree_enter(trim(mdl)//"(), MOM_shared_initialization.F90")
  call MOM_mesg("  MOM_shared_initialization.F90, initialize_topography_named: "//&
                 "TOPO_CONFIG = "//trim(topog_config), 5)

  m_to_Z = 1.0 ; if (present(US)) m_to_Z = US%m_to_Z

  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0, scale=m_to_Z)
  if (max_depth<=0.) call MOM_error(FATAL,"initialize_topography_named: "// &
      "MAXIMUM_DEPTH has a non-sensical value! Was it set?")

  if (trim(topog_config) /= "flat") then
    call get_param(param_file, mdl, "EDGE_DEPTH", Dedge, &
                   "The depth at the edge of one of the named topographies.", &
                   units="m", default=100.0, scale=m_to_Z)
!   call get_param(param_file, mdl, "SOUTHLAT", south_lat, &
!                  "The southern latitude of the domain.", units="degrees", &
!                  fail_if_missing=.true.)
!   call get_param(param_file, mdl, "LENLAT", len_lat, &
!                  "The latitudinal length of the domain.", units="degrees", &
!                  fail_if_missing=.true.)
!   call get_param(param_file, mdl, "WESTLON", west_lon, &
!                  "The western longitude of the domain.", units="degrees", &
!                  default=0.0)
!   call get_param(param_file, mdl, "LENLON", len_lon, &
!                  "The longitudinal length of the domain.", units="degrees", &
!                  fail_if_missing=.true.)
!   call get_param(param_file, mdl, "RAD_EARTH", Rad_Earth, &
!                  "The radius of the Earth.", units="m", default=6.378e6)
    call get_param(param_file, mdl, "TOPOG_SLOPE_SCALE", expdecay, &
                   "The exponential decay scale used in defining some of "//&
                   "the named topographies.", units="m", default=400000.0)
  endif


  PI = 4.0*atan(1.0)

  if (trim(topog_config) == "flat") then
    do i=is,ie ; do j=js,je ; D(i,j) = max_depth ; enddo ; enddo
  elseif (trim(topog_config) == "spoon") then
    D0 = (max_depth - Dedge) / &
             ((1.0 - exp(-0.5*G%len_lat*G%Rad_earth*PI/(180.0 *expdecay))) * &
              (1.0 - exp(-0.5*G%len_lat*G%Rad_earth*PI/(180.0 *expdecay))))
    do i=is,ie ; do j=js,je
  !  This sets a bowl shaped (sort of) bottom topography, with a       !
  !  maximum depth of max_depth.                                   !
      D(i,j) =  Dedge + D0 * &
             (sin(PI * (G%geoLonT(i,j) - (G%west_lon)) / G%len_lon) * &
           (1.0 - exp((G%geoLatT(i,j) - (G%south_lat+G%len_lat))*G%Rad_earth*PI / &
                      (180.0*expdecay)) ))
    enddo ; enddo
  elseif (trim(topog_config) == "bowl") then
    D0 = (max_depth - Dedge) / &
             ((1.0 - exp(-0.5*G%len_lat*G%Rad_earth*PI/(180.0 *expdecay))) * &
              (1.0 - exp(-0.5*G%len_lat*G%Rad_earth*PI/(180.0 *expdecay))))

  !  This sets a bowl shaped (sort of) bottom topography, with a
  !  maximum depth of max_depth.
    do i=is,ie ; do j=js,je
      D(i,j) =  Dedge + D0 * &
             (sin(PI * (G%geoLonT(i,j) - G%west_lon) / G%len_lon) * &
             ((1.0 - exp(-(G%geoLatT(i,j) - G%south_lat)*G%Rad_Earth*PI/ &
                          (180.0*expdecay))) * &
             (1.0 - exp((G%geoLatT(i,j) - (G%south_lat+G%len_lat))* &
                         G%Rad_Earth*PI/(180.0*expdecay)))))
    enddo ; enddo
  elseif (trim(topog_config) == "halfpipe") then
    D0 = max_depth - Dedge
    do i=is,ie ; do j=js,je
      D(i,j) =  Dedge + D0 * ABS(sin(PI*(G%geoLatT(i,j) - G%south_lat)/G%len_lat))
    enddo ; enddo
  else
    call MOM_error(FATAL,"initialize_topography_named: "// &
      "Unrecognized topography name "//trim(topog_config))
  endif

  ! This is here just for safety.  Hopefully it doesn't do anything.
  do i=is,ie ; do j=js,je
    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
  enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine initialize_topography_named
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> limit_topography ensures that  min_depth < D(x,y) < max_depth
subroutine limit_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type), intent(in)    :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                          intent(inout) :: D !< Ocean bottom depth in m or Z if US is present
  type(param_file_type),  intent(in)    :: param_file !< Parameter file structure
  real,                   intent(in)    :: max_depth  !< Maximum depth of model in the units of D
  type(unit_scale_type), optional, intent(in) :: US   !< A dimensional unit scaling type

  ! Local variables
  real :: m_to_Z  ! A dimensional rescaling factor.
  integer :: i, j
  character(len=40)  :: mdl = "limit_topography" ! This subroutine's name.
  real :: min_depth, mask_depth

  call callTree_enter(trim(mdl)//"(), MOM_shared_initialization.F90")

  m_to_Z = 1.0 ; if (present(US)) m_to_Z = US%m_to_Z

  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "If MASKING_DEPTH is unspecified, then anything shallower than "//&
                 "MINIMUM_DEPTH is assumed to be land and all fluxes are masked out. "//&
                 "If MASKING_DEPTH is specified, then all depths shallower than "//&
                 "MINIMUM_DEPTH but deeper than MASKING_DEPTH are rounded to MINIMUM_DEPTH.", &
                 units="m", default=0.0, scale=m_to_Z)
  call get_param(param_file, mdl, "MASKING_DEPTH", mask_depth, &
                 "The depth below which to mask the ocean as land.", &
                 units="m", default=-9999.0, scale=m_to_Z, do_not_log=.true.)

! Make sure that min_depth < D(x,y) < max_depth
  if (mask_depth < -9990.*m_to_Z) then
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      D(i,j) = min( max( D(i,j), 0.5*min_depth ), max_depth )
    enddo ; enddo
  else
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      if (D(i,j)>0.) then
        D(i,j) = min( max( D(i,j), min_depth ), max_depth )
      else
        D(i,j) = 0.
      endif
    enddo ; enddo
  endif

  call callTree_leave(trim(mdl)//'()')
end subroutine limit_topography
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> This subroutine sets up the Coriolis parameter for a sphere
subroutine set_rotation_planetary(f, G, param_file, US)
  type(dyn_horgrid_type), intent(in)  :: G  !< The dynamic horizontal grid
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(out) :: f  !< Coriolis parameter (vertical component) [T-1 ~> s-1]
  type(param_file_type),  intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(unit_scale_type), optional, intent(in) :: US !< A dimensional unit scaling type

! This subroutine sets up the Coriolis parameter for a sphere
  character(len=30) :: mdl = "set_rotation_planetary" ! This subroutine's name.
  integer :: I, J
  real    :: PI
  real    :: omega  ! The planetary rotation rate [T-1 ~> s-1]
  real    :: T_to_s ! A time unit conversion factor

  call callTree_enter(trim(mdl)//"(), MOM_shared_initialization.F90")

  T_to_s = 1.0 ; if (present(US)) T_to_s = US%T_to_s

  call get_param(param_file, "set_rotation_planetary", "OMEGA", omega, &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5, scale=T_to_s)
  PI = 4.0*atan(1.0)

  do I=G%IsdB,G%IedB ; do J=G%JsdB,G%JedB
    f(I,J) = ( 2.0 * omega ) * sin( ( PI * G%geoLatBu(I,J) ) / 180.)
  enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine set_rotation_planetary
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> This subroutine sets up the Coriolis parameter for a beta-plane or f-plane
subroutine set_rotation_beta_plane(f, G, param_file, US)
  type(dyn_horgrid_type), intent(in)  :: G  !< The dynamic horizontal grid
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(out) :: f  !< Coriolis parameter (vertical component) [T-1 ~> s-1]
  type(param_file_type),  intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(unit_scale_type), optional, intent(in) :: US !< A dimensional unit scaling type

! This subroutine sets up the Coriolis parameter for a beta-plane
  integer :: I, J
  real    :: f_0    ! The reference value of the Coriolis parameter [T-1 ~> s-1]
  real    :: beta   ! The meridional gradient of the Coriolis parameter [T-1 m-1 ~> s-1 m-1]
  real    :: y_scl, Rad_Earth
  real    :: T_to_s ! A time unit conversion factor
  real    :: PI
  character(len=40)  :: mdl = "set_rotation_beta_plane" ! This subroutine's name.
  character(len=200) :: axis_units

  call callTree_enter(trim(mdl)//"(), MOM_shared_initialization.F90")

  T_to_s = 1.0 ; if (present(US)) T_to_s = US%T_to_s

  call get_param(param_file, mdl, "F_0", f_0, &
                 "The reference value of the Coriolis parameter with the "//&
                 "betaplane option.", units="s-1", default=0.0, scale=T_to_s)
  call get_param(param_file, mdl, "BETA", beta, &
                 "The northward gradient of the Coriolis parameter with "//&
                 "the betaplane option.", units="m-1 s-1", default=0.0, scale=T_to_s)
  call get_param(param_file, mdl, "AXIS_UNITS", axis_units, default="degrees")

  PI = 4.0*atan(1.0)
  select case (axis_units(1:1))
    case ("d")
      call get_param(param_file, mdl, "RAD_EARTH", Rad_Earth, &
                   "The radius of the Earth.", units="m", default=6.378e6)
      y_scl = Rad_Earth/PI
    case ("k"); y_scl = 1.E3
    case ("m"); y_scl = 1.
    case ("c"); y_scl = 1.E-2
    case default ; call MOM_error(FATAL, &
      " set_rotation_beta_plane: unknown AXIS_UNITS = "//trim(axis_units))
  end select

  do I=G%IsdB,G%IedB ; do J=G%JsdB,G%JedB
    f(I,J) = f_0 + beta * ( G%geoLatBu(I,J) * y_scl )
  enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine set_rotation_beta_plane

!> initialize_grid_rotation_angle initializes the arrays with the sine and
!!   cosine of the angle between logical north on the grid and true north.
subroutine initialize_grid_rotation_angle(G, PF)
  type(dyn_horgrid_type), intent(inout) :: G   !< The dynamic horizontal grid
  type(param_file_type),  intent(in)    :: PF  !< A structure indicating the open file
                                               !! to parse for model parameter values.

  real    :: angle, lon_scale
  real    :: len_lon    ! The periodic range of longitudes, usually 360 degrees.
  real    :: pi_720deg  ! One quarter the conversion factor from degrees to radians.
  real    :: lonB(2,2)  ! The longitude of a point, shifted to have about the same value.
  character(len=40)  :: mdl = "initialize_grid_rotation_angle" ! This subroutine's name.
  logical :: use_bugs
  integer :: i, j, m, n

  call get_param(PF, mdl, "GRID_ROTATION_ANGLE_BUGS", use_bugs, &
                 "If true, use an older algorithm to calculate the sine and "//&
                 "cosines needed rotate between grid-oriented directions and "//&
                 "true north and east.  Differences arise at the tripolar fold.", &
                 default=.false.)

  if (use_bugs) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      lon_scale    = cos((G%geoLatBu(I-1,J-1) + G%geoLatBu(I,J-1  ) + &
                          G%geoLatBu(I-1,J) + G%geoLatBu(I,J)) * atan(1.0)/180)
      angle        = atan2((G%geoLonBu(I-1,J) + G%geoLonBu(I,J) - &
                            G%geoLonBu(I-1,J-1) - G%geoLonBu(I,J-1))*lon_scale, &
                            G%geoLatBu(I-1,J) + G%geoLatBu(I,J) - &
                            G%geoLatBu(I-1,J-1) - G%geoLatBu(I,J-1) )
      G%sin_rot(i,j) = sin(angle) ! angle is the clockwise angle from lat/lon to ocean
      G%cos_rot(i,j) = cos(angle) ! grid (e.g. angle of ocean "north" from true north)
    enddo ; enddo

    ! This is not right at a tripolar or cubed-sphere fold.
    call pass_var(G%cos_rot, G%Domain)
    call pass_var(G%sin_rot, G%Domain)
  else
    pi_720deg = atan(1.0) / 180.0
    len_lon = 360.0 ; if (G%len_lon > 0.0) len_lon = G%len_lon
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      do n=1,2 ; do m=1,2
        lonB(m,n) = modulo_around_point(G%geoLonBu(I+m-2,J+n-2), G%geoLonT(i,j), len_lon)
      enddo ; enddo
      lon_scale = cos(pi_720deg*((G%geoLatBu(I-1,J-1) + G%geoLatBu(I,J)) + &
                                 (G%geoLatBu(I,J-1) + G%geoLatBu(I-1,J)) ) )
      angle = atan2(lon_scale*((lonB(1,2) - lonB(2,1)) + (lonB(2,2) - lonB(1,1))), &
                    (G%geoLatBu(I-1,J) - G%geoLatBu(I,J-1)) + &
                    (G%geoLatBu(I,J) - G%geoLatBu(I-1,J-1)) )
      G%sin_rot(i,j) = sin(angle) ! angle is the clockwise angle from lat/lon to ocean
      G%cos_rot(i,j) = cos(angle) ! grid (e.g. angle of ocean "north" from true north)
    enddo ; enddo

    call pass_vector(G%cos_rot, G%sin_rot, G%Domain, stagger=AGRID)
  endif

end subroutine initialize_grid_rotation_angle

! -----------------------------------------------------------------------------
!> Return the modulo value of x in an interval [xc-(Lx/2) xc+(Lx/2)]
!! If Lx<=0, then it returns x without applying modulo arithmetic.
function modulo_around_point(x, xc, Lx) result(x_mod)
  real, intent(in) :: x  !< Value to which to apply modulo arithmetic
  real, intent(in) :: xc !< Center of modulo range
  real, intent(in) :: Lx !< Modulo range width
  real :: x_mod          !< x shifted by an integer multiple of Lx to be close to xc.

  if (Lx > 0.0) then
    x_mod = modulo(x - (xc - 0.5*Lx), Lx) + (xc - 0.5*Lx)
  else
    x_mod = x
  endif
end function modulo_around_point

! -----------------------------------------------------------------------------
!>   This subroutine sets the open face lengths at selected points to restrict
!! passages to their observed widths based on a named set of sizes.
subroutine reset_face_lengths_named(G, param_file, name, US)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters
  character(len=*),       intent(in)    :: name !< The name for the set of face lengths. Only "global_1deg"
                                                !! is currently implemented.
  type(unit_scale_type), optional, intent(in) :: US !< A dimensional unit scaling type

  ! Local variables
  character(len=256) :: mesg    ! Message for error messages.
  real :: m_to_L  ! A unit conversion factor [L m-1 ~> nondim]
  real :: L_to_m  ! A unit conversion factor [m L-1 ~> nondim]
  real    :: dx_2 = -1.0, dy_2 = -1.0
  real    :: pi_180
  integer :: option = -1
  integer :: i, j, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  pi_180 = (4.0*atan(1.0))/180.0

  select case ( trim(name) )
    case ("global_1deg")    ; option = 1 ; dx_2 = 0.5*1.0
    case default ; call MOM_error(FATAL, "reset_face_lengths_named: "//&
      "Unrecognized channel configuration name "//trim(name))
  end select

  m_to_L = 1.0 ; if (present(US)) m_to_L = US%m_to_L
  L_to_m = 1.0 ; if (present(US)) L_to_m = US%L_to_m

  if (option==1) then ! 1-degree settings.
    do j=jsd,jed ; do I=IsdB,IedB  ! Change any u-face lengths within this loop.
      dy_2 = dx_2 * G%dyCu(I,j)*G%IdxCu(I,j) * cos(pi_180 * G%geoLatCu(I,j))

      if ((abs(G%geoLatCu(I,j)-35.5) < dy_2) .and. (G%geoLonCu(I,j) < -4.5) .and. &
          (G%geoLonCu(I,j) > -6.5)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*12000.0*m_to_L   ! Gibraltar

      if ((abs(G%geoLatCu(I,j)-12.5) < dy_2) .and. (abs(G%geoLonCu(I,j)-43.0) < dx_2)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*10000.0*m_to_L   ! Red Sea

      if ((abs(G%geoLatCu(I,j)-40.5) < dy_2) .and. (abs(G%geoLonCu(I,j)-26.0) < dx_2)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*5000.0*m_to_L   ! Dardanelles

      if ((abs(G%geoLatCu(I,j)-41.5) < dy_2) .and. (abs(G%geoLonCu(I,j)+220.0) < dx_2)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*35000.0*m_to_L   ! Tsugaru strait at 140.0e

      if ((abs(G%geoLatCu(I,j)-45.5) < dy_2) .and. (abs(G%geoLonCu(I,j)+217.5) < 0.9)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*15000.0*m_to_L   ! Betw Hokkaido and Sakhalin at 217&218 = 142e

      ! Greater care needs to be taken in the tripolar region.
      if ((abs(G%geoLatCu(I,j)-80.84) < 0.2) .and. (abs(G%geoLonCu(I,j)+64.9) < 0.8)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*38000.0*m_to_L   ! Smith Sound in Canadian Arch - tripolar region

    enddo ; enddo

    do J=JsdB,JedB ; do i=isd,ied  ! Change any v-face lengths within this loop.
      dy_2 = dx_2 * G%dyCv(i,J)*G%IdxCv(i,J) * cos(pi_180 * G%geoLatCv(i,J))
      if ((abs(G%geoLatCv(i,J)-41.0) < dy_2) .and. (abs(G%geoLonCv(i,J)-28.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*2500.0*m_to_L   ! Bosporus - should be 1000.0 m wide.

      if ((abs(G%geoLatCv(i,J)-13.0) < dy_2) .and. (abs(G%geoLonCv(i,J)-42.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*10000.0*m_to_L   ! Red Sea

      if ((abs(G%geoLatCv(i,J)+2.8) < 0.8) .and. (abs(G%geoLonCv(i,J)+241.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*40000.0*m_to_L   ! Makassar Straits at 241.5 W = 118.5 E

      if ((abs(G%geoLatCv(i,J)-0.56) < 0.5) .and. (abs(G%geoLonCv(i,J)+240.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*80000.0*m_to_L   ! entry to Makassar Straits at 240.5 W = 119.5 E

      if ((abs(G%geoLatCv(i,J)-0.19) < 0.5) .and. (abs(G%geoLonCv(i,J)+230.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*25000.0*m_to_L   ! Channel betw N Guinea and Halmahara 230.5 W = 129.5 E

      if ((abs(G%geoLatCv(i,J)-0.19) < 0.5) .and. (abs(G%geoLonCv(i,J)+229.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*25000.0*m_to_L   ! Channel betw N Guinea and Halmahara 229.5 W = 130.5 E

      if ((abs(G%geoLatCv(i,J)-0.0) < 0.25) .and. (abs(G%geoLonCv(i,J)+228.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*25000.0*m_to_L   ! Channel betw N Guinea and Halmahara 228.5 W = 131.5 E

      if ((abs(G%geoLatCv(i,J)+8.5) < 0.5) .and. (abs(G%geoLonCv(i,J)+244.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*20000.0*m_to_L   ! Lombok Straits at 244.5 W = 115.5 E

      if ((abs(G%geoLatCv(i,J)+8.5) < 0.5) .and. (abs(G%geoLonCv(i,J)+235.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*20000.0*m_to_L   ! Timor Straits at 235.5 W = 124.5 E

      if ((abs(G%geoLatCv(i,J)-52.5) < dy_2) .and. (abs(G%geoLonCv(i,J)+218.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*2500.0*m_to_L    ! Russia and Sakhalin Straits at 218.5 W = 141.5 E

      ! Greater care needs to be taken in the tripolar region.
      if ((abs(G%geoLatCv(i,J)-76.8) < 0.06) .and. (abs(G%geoLonCv(i,J)+88.7) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*8400.0*m_to_L    ! Jones Sound in Canadian Arch - tripolar region

    enddo ; enddo
  endif

  ! These checks apply regardless of the chosen option.

  do j=jsd,jed ; do I=IsdB,IedB
    if (L_to_m*G%dy_Cu(I,j) > L_to_m*G%dyCu(I,j)) then
      write(mesg,'("dy_Cu of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4," at lon/lat of ", ES11.4, ES11.4)') &
                   L_to_m*G%dy_Cu(I,j), L_to_m*G%dyCu(I,j), L_to_m*G%dy_Cu(I,j)-L_to_m*G%dyCu(I,j), &
                   G%geoLonCu(I,j), G%geoLatCu(I,j)
      call MOM_error(FATAL,"reset_face_lengths_named "//mesg)
    endif
    G%areaCu(I,j) = G%dxCu(I,j) * G%dy_Cu(I,j)
    G%IareaCu(I,j) = 0.0
    if (G%areaCu(I,j) > 0.0) G%IareaCu(I,j) = G%mask2dCu(I,j) / (G%areaCu(I,j))
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (L_to_m*G%dx_Cv(i,J) > L_to_m*G%dxCv(i,J)) then
      write(mesg,'("dx_Cv of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4, " at lon/lat of ", ES11.4, ES11.4)') &
                   L_to_m*G%dx_Cv(i,J), L_to_m*G%dxCv(i,J), L_to_m*G%dx_Cv(i,J)-L_to_m*G%dxCv(i,J), &
                   G%geoLonCv(i,J), G%geoLatCv(i,J)

      call MOM_error(FATAL,"reset_face_lengths_named "//mesg)
    endif
    G%areaCv(i,J) = G%dyCv(i,J) * G%dx_Cv(i,J)
    G%IareaCv(i,J) = 0.0
    if (G%areaCv(i,J) > 0.0) G%IareaCv(i,J) = G%mask2dCv(i,J) / (G%areaCv(i,J))
  enddo ; enddo

end subroutine reset_face_lengths_named
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> This subroutine sets the open face lengths at selected points to restrict
!! passages to their observed widths from a arrays read from a file.
subroutine reset_face_lengths_file(G, param_file, US)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(unit_scale_type), optional, intent(in) :: US !< A dimensional unit scaling type

  ! Local variables
  character(len=40)  :: mdl = "reset_face_lengths_file" ! This subroutine's name.
  character(len=256) :: mesg    ! Message for error messages.
  character(len=200) :: filename, chan_file, inputdir ! Strings for file/path
  real :: m_to_L  ! A unit conversion factor [L m-1 ~> nondim]
  real :: L_to_m  ! A unit conversion factor [m L-1 ~> nondim]
  integer :: i, j, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  ! These checks apply regardless of the chosen option.

  call callTree_enter(trim(mdl)//"(), MOM_shared_initialization.F90")
  m_to_L = 1.0 ; if (present(US)) m_to_L = US%m_to_L
  L_to_m = 1.0 ; if (present(US)) L_to_m = US%L_to_m

  call get_param(param_file, mdl, "CHANNEL_WIDTH_FILE", chan_file, &
                 "The file from which the list of narrowed channels is read.", &
                 default="ocean_geometry.nc")
  call get_param(param_file,  mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(chan_file)
  call log_param(param_file, mdl, "INPUTDIR/CHANNEL_WIDTH_FILE", filename)

  if (is_root_pe()) then ; if (.not.file_exists(filename)) &
    call MOM_error(FATAL," reset_face_lengths_file: Unable to open "//&
                           trim(filename))
  endif

  call MOM_read_vector(filename, "dyCuo", "dxCvo", G%dy_Cu, G%dx_Cv, G%Domain, scale=m_to_L)
  call pass_vector(G%dy_Cu, G%dx_Cv, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)

  do j=jsd,jed ; do I=IsdB,IedB
    if (L_to_m*G%dy_Cu(I,j) > L_to_m*G%dyCu(I,j)) then
      write(mesg,'("dy_Cu of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4," at lon/lat of ", ES11.4, ES11.4)') &
                   L_to_m*G%dy_Cu(I,j), L_to_m*G%dyCu(I,j), L_to_m*G%dy_Cu(I,j)-L_to_m*G%dyCu(I,j), &
                   G%geoLonCu(I,j), G%geoLatCu(I,j)
      call MOM_error(FATAL,"reset_face_lengths_file "//mesg)
    endif
    G%areaCu(I,j) = G%dxCu(I,j) * G%dy_Cu(I,j)
    G%IareaCu(I,j) = 0.0
    if (G%areaCu(I,j) > 0.0) G%IareaCu(I,j) = G%mask2dCu(I,j) / (G%areaCu(I,j))
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (L_to_m*G%dx_Cv(i,J) > L_to_m*G%dxCv(i,J)) then
      write(mesg,'("dx_Cv of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4, " at lon/lat of ", ES11.4, ES11.4)') &
                   L_to_m*G%dx_Cv(i,J), L_to_m*G%dxCv(i,J), L_to_m*G%dx_Cv(i,J)-L_to_m*G%dxCv(i,J), &
                   G%geoLonCv(i,J), G%geoLatCv(i,J)

      call MOM_error(FATAL,"reset_face_lengths_file "//mesg)
    endif
    G%areaCv(i,J) = G%dyCv(i,J) * G%dx_Cv(i,J)
    G%IareaCv(i,J) = 0.0
    if (G%areaCv(i,J) > 0.0) G%IareaCv(i,J) = G%mask2dCv(i,J) / (G%areaCv(i,J))
  enddo ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine reset_face_lengths_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> This subroutine sets the open face lengths at selected points to restrict
!! passages to their observed widths from a list read from a file.
subroutine reset_face_lengths_list(G, param_file, US)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(unit_scale_type), optional, intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
  character(len=120), pointer, dimension(:) :: lines => NULL()
  character(len=120) :: line
  character(len=200) :: filename, chan_file, inputdir, mesg ! Strings for file/path
  character(len=40)  :: mdl = "reset_face_lengths_list" ! This subroutine's name.
  real, allocatable, dimension(:,:) :: &
    u_lat, u_lon, v_lat, v_lon ! The latitude and longitude ranges of faces [degrees]
  real, allocatable, dimension(:) :: &
    u_width, v_width      ! The open width of faces [m]
  integer, allocatable, dimension(:) :: &
    u_line_no, v_line_no, &  ! The line numbers in lines of u- and v-face lines
    u_line_used, v_line_used ! The number of times each u- and v-line is used.
  real    :: m_to_L       ! A unit conversion factor [L m-1 ~> nondim]
  real    :: L_to_m       ! A unit conversion factor [m L-1 ~> nondim]
  real    :: lat, lon     ! The latitude and longitude of a point.
  real    :: len_lon      ! The periodic range of longitudes, usually 360 degrees.
  real    :: len_lat      ! The range of latitudes, usually 180 degrees.
  real    :: lon_p, lon_m ! The longitude of a point shifted by 360 degrees.
  logical :: check_360    ! If true, check for longitudes that are shifted by
                          ! +/- 360 degrees from the specified range of values.
  logical :: found_u, found_v
  logical :: unit_in_use
  logical :: fatal_unused_lengths
  integer :: unused
  integer :: ios, iounit, isu, isv
  integer :: last, num_lines, nl_read, ln, npt, u_pt, v_pt
  integer :: i, j, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call callTree_enter(trim(mdl)//"(), MOM_shared_initialization.F90")
  m_to_L = 1.0 ; if (present(US)) m_to_L = US%m_to_L
  L_to_m = 1.0 ; if (present(US)) L_to_m = US%L_to_m

  call get_param(param_file, mdl, "CHANNEL_LIST_FILE", chan_file, &
                 "The file from which the list of narrowed channels is read.", &
                 default="MOM_channel_list")
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(chan_file)
  call log_param(param_file, mdl, "INPUTDIR/CHANNEL_LIST_FILE", filename)
  call get_param(param_file, mdl, "CHANNEL_LIST_360_LON_CHECK", check_360, &
                 "If true, the channel configuration list works for any "//&
                 "longitudes in the range of -360 to 360.", default=.true.)
  call get_param(param_file, mdl, "FATAL_UNUSED_CHANNEL_WIDTHS", fatal_unused_lengths, &
                 "If true, trigger a fatal error if there are any channel widths in "//&
                 "CHANNEL_LIST_FILE that do not cause any open face widths to change.", &
                 default=.false.)

  if (is_root_pe()) then
    ! Open the input file.
    if (.not.file_exists(filename)) call MOM_error(FATAL, &
        " reset_face_lengths_list: Unable to open "//trim(filename))

    ! Find an unused unit number.
    do iounit=10,512
      INQUIRE(iounit,OPENED=unit_in_use) ; if (.not.unit_in_use) exit
    enddo
    if (iounit >= 512) call MOM_error(FATAL, &
        "reset_face_lengths_list: No unused file unit could be found.")

    ! Open the parameter file.
    open(iounit, file=trim(filename), access='SEQUENTIAL', &
         form='FORMATTED', action='READ', position='REWIND', iostat=ios)
    if (ios /= 0) call MOM_error(FATAL, &
            "reset_face_lengths_list: Error opening "//trim(filename))

    ! Count the number of u_width and v_width entries.
    call read_face_length_list(iounit, filename, num_lines, lines)
  endif

  len_lon = 360.0 ; if (G%len_lon > 0.0) len_lon = G%len_lon
  len_lat = 180.0 ; if (G%len_lat > 0.0) len_lat = G%len_lat
  ! Broadcast the number of lines and allocate the required space.
  call broadcast(num_lines, root_PE())
  u_pt = 0 ; v_pt = 0
  if (num_lines > 0) then
    allocate(lines(num_lines))

    allocate(u_lat(2,num_lines)) ; u_lat(:,:) = -1e34
    allocate(u_lon(2,num_lines)) ; u_lon(:,:) = -1e34
    allocate(u_width(num_lines)) ; u_width(:) = -1e34
    allocate(u_line_used(num_lines)) ; u_line_used(:) = 0
    allocate(u_line_no(num_lines)) ; u_line_no(:) = 0

    allocate(v_lat(2,num_lines)) ; v_lat(:,:) = -1e34
    allocate(v_lon(2,num_lines)) ; v_lon(:,:) = -1e34
    allocate(v_width(num_lines)) ; v_width(:) = -1e34
    allocate(v_line_used(num_lines)) ; v_line_used(:) = 0
    allocate(v_line_no(num_lines)) ; v_line_no(:) = 0

    ! Actually read the lines.
    if (is_root_pe()) then
      call read_face_length_list(iounit, filename, nl_read, lines)
      if (nl_read /= num_lines) &
        call MOM_error(FATAL, 'reset_face_lengths_list : Found different '// &
                  'number of valid lines on second reading of '//trim(filename))
      close(iounit) ; iounit = -1
    endif

    ! Broadcast the lines.
    call broadcast(lines, 120, root_PE())

    ! Populate the u_width, etc., data.
    do ln=1,num_lines
      line = lines(ln)
      ! Detect keywords
      found_u = .false.; found_v = .false.
      isu = index(uppercase(line), "U_WIDTH" ); if (isu > 0) found_u = .true.
      isv = index(uppercase(line), "V_WIDTH" ); if (isv > 0) found_v = .true.

      ! Store and check the relevant values.
      if (found_u) then
        u_pt = u_pt + 1
        read(line(isu+8:),*) u_lon(1:2,u_pt), u_lat(1:2,u_pt), u_width(u_pt)
        u_line_no(u_pt) = ln
        if (is_root_PE()) then
          if (check_360) then
            if ((abs(u_lon(1,u_pt)) > len_lon) .or. (abs(u_lon(2,u_pt)) > len_lon)) &
              call MOM_error(WARNING, "reset_face_lengths_list : Out-of-bounds "//&
                 "u-longitude found when reading line "//trim(line)//" from file "//&
                 trim(filename))
            if ((abs(u_lat(1,u_pt)) > len_lat) .or. (abs(u_lat(2,u_pt)) > len_lat)) &
              call MOM_error(WARNING, "reset_face_lengths_list : Out-of-bounds "//&
                 "u-latitude found when reading line "//trim(line)//" from file "//&
                 trim(filename))
          endif
          if (u_lat(1,u_pt) > u_lat(2,u_pt)) &
            call MOM_error(WARNING, "reset_face_lengths_list : Out-of-order "//&
               "u-face latitudes found when reading line "//trim(line)//" from file "//&
               trim(filename))
          if (u_lon(1,u_pt) > u_lon(2,u_pt)) &
            call MOM_error(WARNING, "reset_face_lengths_list : Out-of-order "//&
               "u-face longitudes found when reading line "//trim(line)//" from file "//&
               trim(filename))
          if (u_width(u_pt) < 0.0) &
            call MOM_error(WARNING, "reset_face_lengths_list : Negative "//&
               "u-width found when reading line "//trim(line)//" from file "//&
               trim(filename))
        endif
      elseif (found_v) then
        v_pt = v_pt + 1
        read(line(isv+8:),*) v_lon(1:2,v_pt), v_lat(1:2,v_pt), v_width(v_pt)
        v_line_no(v_pt) = ln
        if (is_root_PE()) then
          if (check_360) then
            if ((abs(v_lon(1,v_pt)) > len_lon) .or. (abs(v_lon(2,v_pt)) > len_lon)) &
              call MOM_error(WARNING, "reset_face_lengths_list : Out-of-bounds "//&
                 "v-longitude found when reading line "//trim(line)//" from file "//&
                 trim(filename))
            if ((abs(v_lat(1,v_pt)) > len_lat) .or. (abs(v_lat(2,v_pt)) > len_lat)) &
              call MOM_error(WARNING, "reset_face_lengths_list : Out-of-bounds "//&
                 "v-latitude found when reading line "//trim(line)//" from file "//&
                 trim(filename))
          endif
          if (v_lat(1,v_pt) > v_lat(2,v_pt)) &
            call MOM_error(WARNING, "reset_face_lengths_list : Out-of-order "//&
               "v-face latitudes found when reading line "//trim(line)//" from file "//&
               trim(filename))
          if (v_lon(1,v_pt) > v_lon(2,v_pt)) &
            call MOM_error(WARNING, "reset_face_lengths_list : Out-of-order "//&
               "v-face longitudes found when reading line "//trim(line)//" from file "//&
               trim(filename))
          if (v_width(v_pt) < 0.0) &
            call MOM_error(WARNING, "reset_face_lengths_list : Negative "//&
               "v-width found when reading line "//trim(line)//" from file "//&
               trim(filename))
        endif
      endif
    enddo

  endif

  do j=jsd,jed ; do I=IsdB,IedB
    lat = G%geoLatCu(I,j) ; lon = G%geoLonCu(I,j)
    if (check_360) then ; lon_p = lon+len_lon ; lon_m = lon-len_lon
    else ; lon_p = lon ; lon_m = lon ; endif

    do npt=1,u_pt
      if (((lat >= u_lat(1,npt)) .and. (lat <= u_lat(2,npt))) .and. &
          (((lon >= u_lon(1,npt)) .and. (lon <= u_lon(2,npt))) .or. &
           ((lon_p >= u_lon(1,npt)) .and. (lon_p <= u_lon(2,npt))) .or. &
           ((lon_m >= u_lon(1,npt)) .and. (lon_m <= u_lon(2,npt)))) ) then

        G%dy_Cu(I,j) = G%mask2dCu(I,j) * m_to_L*min(L_to_m*G%dyCu(I,j), max(u_width(npt), 0.0))
        if (j>=G%jsc .and. j<=G%jec .and. I>=G%isc .and. I<=G%iec) then ! Limit messages/checking to compute domain
          if ( G%mask2dCu(I,j) == 0.0 )  then
            write(stdout,'(A,2F8.2,A,4F8.2,A)') "read_face_lengths_list : G%mask2dCu=0 at ",lat,lon," (",&
                u_lat(1,npt), u_lat(2,npt), u_lon(1,npt), u_lon(2,npt),") so grid metric is unmodified."
          else
            u_line_used(npt) = u_line_used(npt) + 1
            write(stdout,'(A,2F8.2,A,4F8.2,A5,F9.2,A1)') &
                  "read_face_lengths_list : Modifying dy_Cu gridpoint at ",lat,lon," (",&
                  u_lat(1,npt), u_lat(2,npt), u_lon(1,npt), u_lon(2,npt),") to ",L_to_m*G%dy_Cu(I,j),"m"
          endif
        endif
      endif
    enddo

    G%areaCu(I,j) = G%dxCu(I,j) * G%dy_Cu(I,j)
    G%IareaCu(I,j) = 0.0
    if (G%areaCu(I,j) > 0.0) G%IareaCu(I,j) = G%mask2dCu(I,j) / (G%areaCu(I,j))
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    lat = G%geoLatCv(i,J) ; lon = G%geoLonCv(i,J)
    if (check_360) then ; lon_p = lon+len_lon ; lon_m = lon-len_lon
    else ; lon_p = lon ; lon_m = lon ; endif

    do npt=1,v_pt
      if (((lat >= v_lat(1,npt)) .and. (lat <= v_lat(2,npt))) .and. &
          (((lon >= v_lon(1,npt)) .and. (lon <= v_lon(2,npt))) .or. &
           ((lon_p >= v_lon(1,npt)) .and. (lon_p <= v_lon(2,npt))) .or. &
           ((lon_m >= v_lon(1,npt)) .and. (lon_m <= v_lon(2,npt)))) ) then
        G%dx_Cv(i,J) = G%mask2dCv(i,J) * m_to_L*min(L_to_m*G%dxCv(i,J), max(v_width(npt), 0.0))
        if (i>=G%isc .and. i<=G%iec .and. J>=G%jsc .and. J<=G%jec) then ! Limit messages/checking to compute domain
          if ( G%mask2dCv(i,J) == 0.0 )  then
            write(stdout,'(A,2F8.2,A,4F8.2,A)') "read_face_lengths_list : G%mask2dCv=0 at ",lat,lon," (",&
                  v_lat(1,npt), v_lat(2,npt), v_lon(1,npt), v_lon(2,npt),") so grid metric is unmodified."
          else
            v_line_used(npt) = v_line_used(npt) + 1
            write(stdout,'(A,2F8.2,A,4F8.2,A5,F9.2,A1)') &
                  "read_face_lengths_list : Modifying dx_Cv gridpoint at ",lat,lon," (",&
                  v_lat(1,npt), v_lat(2,npt), v_lon(1,npt), v_lon(2,npt),") to ",L_to_m*G%dx_Cv(I,j),"m"
          endif
        endif
      endif
    enddo

    G%areaCv(i,J) = G%dyCv(i,J) * G%dx_Cv(i,J)
    G%IareaCv(i,J) = 0.0
    if (G%areaCv(i,J) > 0.0) G%IareaCv(i,J) = G%mask2dCv(i,J) / (G%areaCv(i,J))
  enddo ; enddo

  ! Verify that all channel widths have been used
  unused = 0
  if (u_pt > 0) call sum_across_PEs(u_line_used, u_pt)
  if (v_pt > 0) call sum_across_PEs(v_line_used, v_pt)
  if (is_root_PE()) then
    unused = 0
    do npt=1,u_pt ; if (u_line_used(npt) == 0) then
      call MOM_error(WARNING, "reset_face_lengths_list unused u-face line: "//&
                     trim(lines(u_line_no(npt))) )
      unused = unused + 1
    endif ; enddo
    do npt=1,v_pt ; if (v_line_used(npt) == 0) then
      call MOM_error(WARNING, "reset_face_lengths_list unused v-face line: "//&
                     trim(lines(v_line_no(npt))) )
      unused = unused + 1
    endif ; enddo
    if (fatal_unused_lengths .and. (unused > 0)) call MOM_error(FATAL, &
      "reset_face_lengths_list causing MOM6 abort due to unused face length lines.")
  endif

  if (num_lines > 0) then
    deallocate(lines)
    deallocate(u_line_used, v_line_used, u_line_no, v_line_no)
    deallocate(u_lat) ; deallocate(u_lon) ; deallocate(u_width)
    deallocate(v_lat) ; deallocate(v_lon) ; deallocate(v_width)
  endif

  call callTree_leave(trim(mdl)//'()')
end subroutine reset_face_lengths_list
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!>   This subroutine reads and counts the non-blank lines in the face length list file, after removing comments.
subroutine read_face_length_list(iounit, filename, num_lines, lines)
  integer,                          intent(in)  :: iounit    !< An open I/O unit number for the file
  character(len=*),                 intent(in)  :: filename  !< The name of the face-length file to read
  integer,                          intent(out) :: num_lines !< The number of non-blank lines in the file
  character(len=120), dimension(:), pointer     :: lines  !< The non-blank lines, after removing comments

  !   This subroutine reads and counts the non-blank lines in the face length
  ! list file, after removing comments.
  character(len=120) :: line, line_up
  logical :: found_u, found_v
  integer :: isu, isv, icom, verbose
  integer :: last

  num_lines = 0

  if (iounit <= 0) return
  rewind(iounit)
  do while(.true.)
    read(iounit, '(a)', end=8, err=9) line
    last = len_trim(line)
    ! Eliminate either F90 or C comments from the line.
    icom = index(line(:last), "!") ; if (icom > 0) last = icom-1
    icom = index(line(:last), "/*") ; if (icom > 0) last = icom-1
    if (last < 1) cycle

    ! Detect keywords
    line_up = uppercase(line)
    found_u = .false.; found_v = .false.
    isu = index(line_up(:last), "U_WIDTH" ); if (isu > 0) found_u = .true.
    isv = index(line_up(:last), "V_WIDTH" ); if (isv > 0) found_v = .true.

    if (found_u .and. found_v) call MOM_error(FATAL, &
      "read_face_length_list : both U_WIDTH and V_WIDTH found when "//&
      "reading the line "//trim(line(:last))//" in file "//trim(filename))
    if (found_u .or. found_v) then
      num_lines = num_lines + 1
      if (associated(lines)) then
        lines(num_lines) = line(1:last)
      endif
    endif
  enddo ! while (.true.)

8 continue
  return

9 call MOM_error(FATAL, "read_face_length_list : "//&
                  "Error while reading file "//trim(filename))

end subroutine read_face_length_list
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set the bathymetry at velocity points to be the maximum of the depths at the
!! neighoring tracer points.
subroutine set_velocity_depth_max(G)
  type(dyn_horgrid_type), intent(inout) :: G   !< The dynamic horizontal grid
  ! This subroutine sets the 4 bottom depths at velocity points to be the
  ! maximum of the adjacent depths.
  integer :: i, j

  do I=G%isd,G%ied-1 ; do j=G%jsd,G%jed
    G%Dblock_u(I,j) = G%mask2dCu(I,j) * max(G%bathyT(i,j), G%bathyT(i+1,j))
    G%Dopen_u(I,j) = G%Dblock_u(I,j)
  enddo ; enddo
  do i=G%isd,G%ied ; do J=G%jsd,G%jed-1
    G%Dblock_v(I,J) = G%mask2dCv(i,J) * max(G%bathyT(i,j), G%bathyT(i,j+1))
    G%Dopen_v(I,J) = G%Dblock_v(I,J)
  enddo ; enddo
end subroutine set_velocity_depth_max
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Set the bathymetry at velocity points to be the minimum of the depths at the
!! neighoring tracer points.
subroutine set_velocity_depth_min(G)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid
  ! This subroutine sets the 4 bottom depths at velocity points to be the
  ! minimum of the adjacent depths.
  integer :: i, j

  do I=G%isd,G%ied-1 ; do j=G%jsd,G%jed
    G%Dblock_u(I,j) = G%mask2dCu(I,j) * min(G%bathyT(i,j), G%bathyT(i+1,j))
    G%Dopen_u(I,j) = G%Dblock_u(I,j)
  enddo ; enddo
  do i=G%isd,G%ied ; do J=G%jsd,G%jed-1
    G%Dblock_v(I,J) = G%mask2dCv(i,J) * min(G%bathyT(i,j), G%bathyT(i,j+1))
    G%Dopen_v(I,J) = G%Dblock_v(I,J)
  enddo ; enddo
end subroutine set_velocity_depth_min
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Pre-compute global integrals of grid quantities (like masked ocean area) for
!! later use in reporting diagnostics
subroutine compute_global_grid_integrals(G, US)
  type(dyn_horgrid_type),          intent(inout) :: G  !< The dynamic horizontal grid
  type(unit_scale_type), optional, intent(in)    :: US !< A dimensional unit scaling type

  ! Local variables
  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: tmpForSumming
  real :: area_scale  ! A scaling factor for area into MKS units
  integer :: i,j

  area_scale = 1.0 ; if (present(US)) area_scale = US%L_to_m**2

  tmpForSumming(:,:) = 0.
  G%areaT_global = 0.0 ; G%IareaT_global = 0.0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    tmpForSumming(i,j) = area_scale*G%areaT(i,j) * G%mask2dT(i,j)
  enddo ; enddo
  G%areaT_global = reproducing_sum(tmpForSumming)

  if (G%areaT_global == 0.0) &
    call MOM_error(FATAL, "compute_global_grid_integrals: "//&
                    "zero ocean area (check topography?)")

  G%IareaT_global = 1.0 / (G%areaT_global)
end subroutine compute_global_grid_integrals
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Write out a file describing the topography, Coriolis parameter, grid locations
!! and various other fixed fields from the grid.
subroutine write_ocean_geometry_file(G, param_file, directory, geom_file, US)
  type(dyn_horgrid_type),       intent(inout) :: G         !< The dynamic horizontal grid
  type(param_file_type),        intent(in)    :: param_file !< Parameter file structure
  character(len=*),             intent(in)    :: directory !< The directory into which to place the geometry file.
  character(len=*),   optional, intent(in)    :: geom_file !< If present, the name of the geometry file
                                                           !! (otherwise the file is "ocean_geometry")
  type(unit_scale_type), optional, intent(in) :: US        !< A dimensional unit scaling type

  ! Local variables.
  character(len=240) :: filepath
  character(len=40)  :: mdl = "write_ocean_geometry_file"
  integer, parameter :: nFlds=23
  type(vardesc) :: vars(nFlds)
  type(fieldtype) :: fields(nFlds)
  real :: Z_to_m_scale ! A unit conversion factor from Z to m
  real :: s_to_T_scale ! A unit conversion factor from T-1 to s-1
  real :: L_to_m_scale ! A unit conversion factor from L to m
  type(file_type) :: IO_handle ! The I/O handle of the fileset
  integer :: file_threading
  integer :: nFlds_used
  logical :: multiple_files

  call callTree_enter('write_ocean_geometry_file()')

  Z_to_m_scale = 1.0 ; if (present(US)) Z_to_m_scale = US%Z_to_m
  s_to_T_scale = 1.0 ; if (present(US)) s_to_T_scale = US%s_to_T
  L_to_m_scale = 1.0 ; if (present(US)) L_to_m_scale = US%L_to_m

  !   var_desc populates a type defined in MOM_io.F90.  The arguments, in order, are:
  ! (1) the variable name for the NetCDF file
  ! (2) the units of the variable when output
  ! (3) the variable's long name
  ! (4) a character indicating the  horizontal grid, which may be '1' (column),
  !     'h', 'q', 'u', or 'v', for the corresponding C-grid variable
  ! (5) a character indicating the vertical grid, which may be 'L' (layer),
  !     'i' (interface), or '1' (no vertical location)
  ! (6) a character indicating the time levels of the field, which may be
  !    's' (snap-shot), 'p' (periodic), or '1' (no time variation)
  vars(1) = var_desc("geolatb","degree","latitude at corner (Bu) points",'q','1','1')
  vars(2) = var_desc("geolonb","degree","longitude at corner (Bu) points",'q','1','1')
  vars(3) = var_desc("geolat","degree", "latitude at tracer (T) points", 'h','1','1')
  vars(4) = var_desc("geolon","degree","longitude at tracer (T) points",'h','1','1')
  vars(5) = var_desc("D","meter","Basin Depth",'h','1','1')
  vars(6) = var_desc("f","s-1","Coriolis Parameter",'q','1','1')
  vars(7) = var_desc("dxCv","m","Zonal grid spacing at v points",'v','1','1')
  vars(8) = var_desc("dyCu","m","Meridional grid spacing at u points",'u','1','1')
  vars(9) = var_desc("dxCu","m","Zonal grid spacing at u points",'u','1','1')
  vars(10)= var_desc("dyCv","m","Meridional grid spacing at v points",'v','1','1')
  vars(11)= var_desc("dxT","m","Zonal grid spacing at h points",'h','1','1')
  vars(12)= var_desc("dyT","m","Meridional grid spacing at h points",'h','1','1')
  vars(13)= var_desc("dxBu","m","Zonal grid spacing at q points",'q','1','1')
  vars(14)= var_desc("dyBu","m","Meridional grid spacing at q points",'q','1','1')
  vars(15)= var_desc("Ah","m2","Area of h cells",'h','1','1')
  vars(16)= var_desc("Aq","m2","Area of q cells",'q','1','1')

  vars(17)= var_desc("dxCvo","m","Open zonal grid spacing at v points",'v','1','1')
  vars(18)= var_desc("dyCuo","m","Open meridional grid spacing at u points",'u','1','1')
  vars(19)= var_desc("wet", "nondim", "land or ocean?", 'h','1','1')

  vars(20) = var_desc("Dblock_u","m","Blocked depth at u points",'u','1','1')
  vars(21) = var_desc("Dopen_u","m","Open depth at u points",'u','1','1')
  vars(22) = var_desc("Dblock_v","m","Blocked depth at v points",'v','1','1')
  vars(23) = var_desc("Dopen_v","m","Open depth at v points",'v','1','1')


  nFlds_used = 19 ; if (G%bathymetry_at_vel) nFlds_used = 23

  if (present(geom_file)) then
    filepath = trim(directory) // trim(geom_file)
  else
    filepath = trim(directory) // "ocean_geometry"
  endif

  call get_param(param_file, mdl, "PARALLEL_RESTARTFILES", multiple_files, &
                 "If true, each processor writes its own restart file, "//&
                 "otherwise a single restart file is generated", &
                 default=.false.)
  file_threading = SINGLE_FILE
  if (multiple_files) file_threading = MULTIPLE

  call create_file(IO_handle, trim(filepath), vars, nFlds_used, fields, file_threading, dG=G)

  call MOM_write_field(IO_handle, fields(1), G%Domain, G%geoLatBu)
  call MOM_write_field(IO_handle, fields(2), G%Domain, G%geoLonBu)
  call MOM_write_field(IO_handle, fields(3), G%Domain, G%geoLatT)
  call MOM_write_field(IO_handle, fields(4), G%Domain, G%geoLonT)

  call MOM_write_field(IO_handle, fields(5), G%Domain, G%bathyT, scale=Z_to_m_scale)
  call MOM_write_field(IO_handle, fields(6), G%Domain, G%CoriolisBu, scale=s_to_T_scale)

  call MOM_write_field(IO_handle, fields(7),  G%Domain, G%dxCv, scale=L_to_m_scale)
  call MOM_write_field(IO_handle, fields(8),  G%Domain, G%dyCu, scale=L_to_m_scale)
  call MOM_write_field(IO_handle, fields(9),  G%Domain, G%dxCu, scale=L_to_m_scale)
  call MOM_write_field(IO_handle, fields(10), G%Domain, G%dyCv, scale=L_to_m_scale)
  call MOM_write_field(IO_handle, fields(11), G%Domain, G%dxT, scale=L_to_m_scale)
  call MOM_write_field(IO_handle, fields(12), G%Domain, G%dyT, scale=L_to_m_scale)
  call MOM_write_field(IO_handle, fields(13), G%Domain, G%dxBu, scale=L_to_m_scale)
  call MOM_write_field(IO_handle, fields(14), G%Domain, G%dyBu, scale=L_to_m_scale)

  call MOM_write_field(IO_handle, fields(15), G%Domain, G%areaT, scale=L_to_m_scale**2)
  call MOM_write_field(IO_handle, fields(16), G%Domain, G%areaBu, scale=L_to_m_scale**2)

  call MOM_write_field(IO_handle, fields(17), G%Domain, G%dx_Cv, scale=L_to_m_scale)
  call MOM_write_field(IO_handle, fields(18), G%Domain, G%dy_Cu, scale=L_to_m_scale)
  call MOM_write_field(IO_handle, fields(19), G%Domain, G%mask2dT)

  if (G%bathymetry_at_vel) then
    call MOM_write_field(IO_handle, fields(20), G%Domain, G%Dblock_u, scale=Z_to_m_scale)
    call MOM_write_field(IO_handle, fields(21), G%Domain, G%Dopen_u, scale=Z_to_m_scale)
    call MOM_write_field(IO_handle, fields(22), G%Domain, G%Dblock_v, scale=Z_to_m_scale)
    call MOM_write_field(IO_handle, fields(23), G%Domain, G%Dopen_v, scale=Z_to_m_scale)
  endif

  call close_file(IO_handle)

  call callTree_leave('write_ocean_geometry_file()')
end subroutine write_ocean_geometry_file

end module MOM_shared_initialization
