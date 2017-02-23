!> Code that initializes fixed aspects of the model grid, such as horizontal
!! grid metrics, topography and Coriolis, and can be shared between components.
module MOM_shared_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms, only : max_across_PEs
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE, AGRID
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, log_param, param_file_type, log_version
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : read_data, SINGLE_FILE, MULTIPLE
use MOM_io, only : slasher, vardesc, write_field, var_desc
use MOM_string_functions, only : uppercase

use netcdf

implicit none ; private

public MOM_shared_init_init
public MOM_initialize_rotation, MOM_calculate_grad_Coriolis
public initialize_topography_from_file, apply_topography_edits_from_file
public initialize_topography_named, limit_topography, diagnoseMaximumDepth
public set_rotation_planetary, set_rotation_beta_plane, initialize_grid_rotation_angle
public reset_face_lengths_named, reset_face_lengths_file, reset_face_lengths_list
public read_face_length_list, set_velocity_depth_max, set_velocity_depth_min
public compute_global_grid_integrals, write_ocean_geometry_file

contains

! -----------------------------------------------------------------------------
!> MOM_shared_init_init just writes the code version.
subroutine MOM_shared_init_init(PF)
  type(param_file_type),   intent(in)    :: PF   !< A structure indicating the open file
                                                 !! to parse for model parameter values.

  character(len=40)  :: mod = "MOM_shared_initialization" ! This module's name.

! This include declares and sets the variable "version".
#include "version_variable.h"
  call log_version(PF, mod, version, &
   "Sharable code to initialize time-invariant fields, like bathymetry and Coriolis parameters.")

end subroutine MOM_shared_init_init
! -----------------------------------------------------------------------------

!> MOM_initialize_rotation makes the appropriate call to set up the Coriolis parameter.
subroutine MOM_initialize_rotation(f, G, PF)
  type(dyn_horgrid_type),                       intent(in)  :: G  !< The dynamic horizontal grid type
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), intent(out) :: f  !< The Coriolis parameter in s-1
  type(param_file_type),                        intent(in)  :: PF !< Parameter file structure

!   This subroutine makes the appropriate call to set up the Coriolis parameter.
! This is a separate subroutine so that it can be made public and shared with
! the ice-sheet code or other components.
! Set up the Coriolis parameter, f, either analytically or from file.
  character(len=40)  :: mod = "MOM_initialize_rotation" ! This subroutine's name.
  character(len=200) :: config

  call callTree_enter(trim(mod)//"(), MOM_shared_initialization.F90")
  call get_param(PF, mod, "ROTATION", config, &
                 "This specifies how the Coriolis parameter is specified: \n"//&
                 " \t 2omegasinlat - Use twice the planetary rotation rate \n"//&
                 " \t\t times the sine of latitude.\n"//&
                 " \t betaplane - Use a beta-plane or f-plane. \n"//&
                 " \t USER - call a user modified routine.", &
                 default="2omegasinlat")
  select case (trim(config))
    case ("2omegasinlat"); call set_rotation_planetary(f, G, PF)
    case ("beta"); call set_rotation_beta_plane(f, G, PF)
    case ("betaplane"); call set_rotation_beta_plane(f, G, PF)
   !case ("nonrotating") ! Note from AJA: Missing case?
    case default ; call MOM_error(FATAL,"MOM_initialize: "// &
      "Unrecognized rotation setup "//trim(config))
  end select
  call callTree_leave(trim(mod)//'()')
end subroutine MOM_initialize_rotation

!> Calculates the components of grad f (Coriolis parameter)
subroutine MOM_calculate_grad_Coriolis(dF_dx, dF_dy, G)
  type(dyn_horgrid_type),             intent(inout) :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out)   :: dF_dx !< x-component of grad f
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out)   :: dF_dy !< y-component of grad f
  ! Local variables
  integer :: i,j
  real :: f1, f2

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

!> Return the global maximum ocean bottom depth in m.
function diagnoseMaximumDepth(D,G)
  type(dyn_horgrid_type),  intent(in) :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(in) :: D !< Ocean bottom depth in m
  real :: diagnoseMaximumDepth             !< The global maximum ocean bottom depth in m
  ! Local variables
  integer :: i,j
  diagnoseMaximumDepth=D(G%isc,G%jsc)
  do j=G%jsc, G%jec
    do i=G%isc, G%iec
      diagnoseMaximumDepth=max(diagnoseMaximumDepth,D(i,j))
    enddo
  enddo
  call max_across_PEs(diagnoseMaximumDepth)
end function diagnoseMaximumDepth


!> Read gridded depths from file
subroutine initialize_topography_from_file(D, G, param_file)
  type(dyn_horgrid_type),           intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),            intent(in)  :: param_file !< Parameter file structure
  ! Local variables
  character(len=200) :: filename, topo_file, inputdir ! Strings for file/path
  character(len=200) :: topo_varname                  ! Variable name in file
  character(len=40)  :: mod = "initialize_topography_from_file" ! This subroutine's name.

  call callTree_enter(trim(mod)//"(), MOM_shared_initialization.F90")

  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mod, "TOPO_FILE", topo_file, &
                 "The file from which the bathymetry is read.", &
                 default="topog.nc")
  call get_param(param_file, mod, "TOPO_VARNAME", topo_varname, &
                 "The name of the bathymetry variable in TOPO_FILE.", &
                 default="depth")

  filename = trim(inputdir)//trim(topo_file)
  call log_param(param_file, mod, "INPUTDIR/TOPO_FILE", filename)

  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topography_from_file: Unable to open "//trim(filename))

  D(:,:) = -9.E30 ! Initializing to a very large negative depth (tall mountains)
                  ! everywhere before reading from a file should do nothing.
                  ! However, in the instance of masked-out PEs, halo regions
                  ! are not updated when a processor does not exist. We need to
                  ! ensure the depth in masked-out PEs appears to be that of land
                  ! so this line does that in the halo regions. For non-masked PEs
                  ! the halo region is filled properly with a later pass_var().
  call read_data(filename,trim(topo_varname),D,domain=G%Domain%mpp_domain)

  call apply_topography_edits_from_file(D, G, param_file)

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_topography_from_file

!> Applies a list of topography overrides read from a netcdf file
subroutine apply_topography_edits_from_file(D, G, param_file)
  type(dyn_horgrid_type),           intent(in)    :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(inout) :: D !< Ocean bottom depth in m
  type(param_file_type),            intent(in)    :: param_file !< Parameter file structure

  ! Local variables
  character(len=200) :: topo_edits_file, inputdir ! Strings for file/path
  character(len=40)  :: mod = "apply_topography_edits_from_file" ! This subroutine's name.
  integer :: n_edits, n, ashape(5), i, j, ncid, id, ncstatus, iid, jid, zid
  integer, dimension(:), allocatable :: ig, jg
  real, dimension(:), allocatable :: new_depth

  call callTree_enter(trim(mod)//"(), MOM_shared_initialization.F90")

  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mod, "TOPO_EDITS_FILE", topo_edits_file, &
                 "The file from which to read a list of i,j,z topography overrides.", &
                 default="")

  if (len_trim(topo_edits_file)==0) return

  topo_edits_file = trim(inputdir)//trim(topo_edits_file)
  if (.not.file_exists(topo_edits_file, G%Domain)) call MOM_error(FATAL, &
     'initialize_topography_from_file: Unable to open '//trim(topo_edits_file))

  ncstatus = nf90_open(trim(topo_edits_file), NF90_NOWRITE, ncid)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                                'Failed to open '//trim(topo_edits_file))

  ! Get nEdits
  ncstatus = nf90_inq_dimid(ncid, 'nEdits', id)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                                'Failed to inq_dimid nEdits for '//trim(topo_edits_file))
  ncstatus = nf90_inquire_dimension(ncid, id, len=n_edits)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                                'Failed to inquire_dimension nEdits for '//trim(topo_edits_file))

  ! Read ni
  ncstatus = nf90_inq_varid(ncid, 'ni', id)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                                'Failed to inq_varid ni for '//trim(topo_edits_file))
  ncstatus = nf90_get_var(ncid, id, i)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                              'Failed to get_var ni for '//trim(topo_edits_file))
  if (i /= G%ieg) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                              'Incompatible i-dimension of grid in '//trim(topo_edits_file))

  ! Read nj
  ncstatus = nf90_inq_varid(ncid, 'nj', id)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                                'Failed to inq_varid nj for '//trim(topo_edits_file))
  ncstatus = nf90_get_var(ncid, id, j)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                              'Failed to get_var nj for '//trim(topo_edits_file))
  if (j /= G%jeg) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                              'Incompatible j-dimension of grid in '//trim(topo_edits_file))

  ! Read iEdit
  ncstatus = nf90_inq_varid(ncid, 'iEdit', id)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                                'Failed to inq_varid iEdit for '//trim(topo_edits_file))
  allocate(ig(n_edits))
  ncstatus = nf90_get_var(ncid, id, ig)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                              'Failed to get_var iEdit for '//trim(topo_edits_file))

  ! Read jEdit
  ncstatus = nf90_inq_varid(ncid, 'jEdit', id)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                                'Failed to inq_varid jEdit for '//trim(topo_edits_file))
  allocate(jg(n_edits))
  ncstatus = nf90_get_var(ncid, id, jg)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                              'Failed to get_var jEdit for '//trim(topo_edits_file))

  ! Read zEdit
  ncstatus = nf90_inq_varid(ncid, 'zEdit', id)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                                'Failed to inq_varid zEdit for '//trim(topo_edits_file))
  allocate(new_depth(n_edits))
  ncstatus = nf90_get_var(ncid, id, new_depth)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                              'Failed to get_var zEdit for '//trim(topo_edits_file))

  ! Close file
  ncstatus = nf90_close(ncid)
  if (ncstatus /= NF90_NOERR) call MOM_error(FATAL, 'apply_topography_edits_from_file: '//&
                                'Failed to close '//trim(topo_edits_file))

  do n = 1, n_edits
    i = ig(n) - G%isd_global + 2 ! +1 for python indexing and +1 for ig-isd_global+1
    j = jg(n) - G%jsd_global + 2
    if (i>=G%isc .and. i<=G%iec .and. j>=G%jsc .and. j<=G%jec) then
      if (new_depth(n)/=0.) then
        write(*,'(a,3i5,f8.2,a,f8.2,2i4)') 'Ocean topography edit: ',n,ig(n),jg(n),D(i,j),'->',abs(new_depth(n)),i,j
        D(i,j) = abs(new_depth(n)) ! Allows for height-file edits (i.e. converts negatives)
      else
        call MOM_error(FATAL, ' apply_topography_edits_from_file: '//&
          "A zero depth edit would change the land mask and is not allowed in"//trim(topo_edits_file))
      endif
    endif
  enddo

  deallocate( ig, jg, new_depth )

  call callTree_leave(trim(mod)//'()')
end subroutine apply_topography_edits_from_file

!> initialize the bathymetry based on one of several named idealized configurations
subroutine initialize_topography_named(D, G, param_file, topog_config, max_depth)
  type(dyn_horgrid_type),           intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),            intent(in)  :: param_file !< Parameter file structure
  character(len=*),                 intent(in)  :: topog_config !< The name of an idealized
                                                              !! topographic configuration
  real,                             intent(in)  :: max_depth  !< Maximum depth of model in m

! Arguments: D          - the bottom depth in m. Intent out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      topog_config - The name of an idealized topographic configuration.
!  (in)      max_depth  - The maximum depth in m.

! This subroutine places the bottom depth in m into D(:,:), shaped in a spoon
  real :: min_depth            ! The minimum depth in m.
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: D0                   ! A constant to make the maximum     !
                               ! basin depth MAXIMUM_DEPTH.         !
  real :: expdecay             ! A decay scale of associated with   !
                               ! the sloping boundaries, in m.      !
  real :: Dedge                ! The depth in m at the basin edge.  !
! real :: south_lat, west_lon, len_lon, len_lat, Rad_earth
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  character(len=40)  :: mod = "initialize_topography_named" ! This subroutine's name.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call callTree_enter(trim(mod)//"(), MOM_shared_initialization.F90")
  call MOM_mesg("  MOM_shared_initialization.F90, initialize_topography_named: "//&
                 "TOPO_CONFIG = "//trim(topog_config), 5)

  call get_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)
  if (max_depth<=0.) call MOM_error(FATAL,"initialize_topography_named: "// &
      "MAXIMUM_DEPTH has a non-sensical value! Was it set?")

  if (trim(topog_config) /= "flat") then
    call get_param(param_file, mod, "EDGE_DEPTH", Dedge, &
                   "The depth at the edge of one of the named topographies.", &
                   units="m", default=100.0)
!   call get_param(param_file, mod, "SOUTHLAT", south_lat, &
!                  "The southern latitude of the domain.", units="degrees", &
!                  fail_if_missing=.true.)
!   call get_param(param_file, mod, "LENLAT", len_lat, &
!                  "The latitudinal length of the domain.", units="degrees", &
!                  fail_if_missing=.true.)
!   call get_param(param_file, mod, "WESTLON", west_lon, &
!                  "The western longitude of the domain.", units="degrees", &
!                  default=0.0)
!   call get_param(param_file, mod, "LENLON", len_lon, &
!                  "The longitudinal length of the domain.", units="degrees", &
!                  fail_if_missing=.true.)
!   call get_param(param_file, mod, "RAD_EARTH", Rad_Earth, &
!                  "The radius of the Earth.", units="m", default=6.378e6)
    call get_param(param_file, mod, "TOPOG_SLOPE_SCALE", expdecay, &
                   "The exponential decay scale used in defining some of \n"//&
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

  call callTree_leave(trim(mod)//'()')
end subroutine initialize_topography_named
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> limit_topography ensures that  min_depth < D(x,y) < max_depth
subroutine limit_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type), intent(in)    :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                          intent(inout) :: D !< Ocean bottom depth in m
  type(param_file_type),  intent(in)    :: param_file !< Parameter file structure
  real,                   intent(in)    :: max_depth  !< Maximum depth of model in m
! Arguments: D          - the bottom depth in m. Intent in/out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      max_depth  - The maximum depth in m.

! This subroutine ensures that    min_depth < D(x,y) < max_depth
  integer :: i, j
  character(len=40)  :: mod = "limit_topography" ! This subroutine's name.
  real :: min_depth, mask_depth

  call callTree_enter(trim(mod)//"(), MOM_shared_initialization.F90")

  call get_param(param_file, mod, "MINIMUM_DEPTH", min_depth, &
                 "If MASKING_DEPTH is unspecified, then anything shallower than\n"//&
                 "MINIMUM_DEPTH is assumed to be land and all fluxes are masked out.\n"//&
                 "If MASKING_DEPTH is specified, then all depths shallower than\n"//&
                 "MINIMUM_DEPTH but deeper than MASKING_DEPTH are rounded to MINIMUM_DEPTH.", &
                 units="m", default=0.0)
  call get_param(param_file, mod, "MASKING_DEPTH", mask_depth, &
                 "The depth below which to mask the ocean as land.", units="m", &
                 default=-9999.0, do_not_log=.true.)

! Make sure that min_depth < D(x,y) < max_depth
  if (mask_depth<-9990.) then
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

  call callTree_leave(trim(mod)//'()')
end subroutine limit_topography
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_rotation_planetary(f, G, param_file)
  type(dyn_horgrid_type),                        intent(in)  :: G
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), intent(out) :: f
  type(param_file_type),                        intent(in)  :: param_file
! Arguments: f          - Coriolis parameter (vertical component) in s^-1
!     (in)   G          - grid type
!     (in)   param_file - parameter file type

! This subroutine sets up the Coriolis parameter for a sphere
  character(len=30) :: mod = "set_rotation_planetary" ! This subroutine's name.
  integer :: I, J
  real    :: PI, omega

  call callTree_enter(trim(mod)//"(), MOM_shared_initialization.F90")

  call get_param(param_file, "set_rotation_planetary", "OMEGA", omega, &
                 "The rotation rate of the earth.", units="s-1", &
                 default=7.2921e-5)
  PI = 4.0*atan(1.0)

  do I=G%IsdB,G%IedB ; do J=G%JsdB,G%JedB
    f(I,J) = ( 2.0 * omega ) * sin( ( PI * G%geoLatBu(I,J) ) / 180.)
  enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_rotation_planetary
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_rotation_beta_plane(f, G, param_file)
  type(dyn_horgrid_type),                        intent(in)  :: G
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), intent(out) :: f
  type(param_file_type),                        intent(in)  :: param_file
! Arguments: f          - Coriolis parameter (vertical component) in s^-1
!     (in)   G          - grid type
!     (in)   param_file - parameter file type

! This subroutine sets up the Coriolis parameter for a beta-plane
  integer :: I, J
  real    :: f_0, beta, y_scl, Rad_Earth, PI
  character(len=40)  :: mod = "set_rotation_beta_plane" ! This subroutine's name.
  character(len=200) :: axis_units

  call callTree_enter(trim(mod)//"(), MOM_shared_initialization.F90")

  call get_param(param_file, mod, "F_0", f_0, &
                 "The reference value of the Coriolis parameter with the \n"//&
                 "betaplane option.", units="s-1", default=0.0)
  call get_param(param_file, mod, "BETA", beta, &
                 "The northward gradient of the Coriolis parameter with \n"//&
                 "the betaplane option.", units="m-1 s-1", default=0.0)
  call get_param(param_file, mod, "AXIS_UNITS", axis_units, default="degrees")

  PI = 4.0*atan(1.0)
  select case (axis_units(1:1))
    case ("d")
      call get_param(param_file, mod, "RAD_EARTH", Rad_Earth, &
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

  call callTree_leave(trim(mod)//'()')
end subroutine set_rotation_beta_plane

!> initialize_grid_rotation_angle initializes the arrays with the sine and
!!   cosine of the angle between logical north on the grid and true north.
subroutine initialize_grid_rotation_angle(G, PF)
  type(dyn_horgrid_type), intent(inout) :: G   !< The model's horizontal grid structure.
  type(param_file_type),  intent(in)    :: PF  !< A structure indicating the open file
                                               !! to parse for model parameter values.

  real    :: angle, lon_scale
  integer :: i, j

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

  ! ### THIS DOESN'T SEEM RIGHT AT A CUBED-SPHERE FOLD -RWH
  call pass_var(G%cos_rot, G%Domain)
  call pass_var(G%sin_rot, G%Domain)

end subroutine initialize_grid_rotation_angle

! -----------------------------------------------------------------------------
subroutine reset_face_lengths_named(G, param_file, name)
  type(dyn_horgrid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  character(len=*),      intent(in)    :: name
!   This subroutine sets the open face lengths at selected points to restrict
! passages to their observed widths.

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      name - The name for the set of face lengths.
  character(len=256) :: mesg    ! Message for error messages.
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

  if (option==1) then ! 1-degree settings.
    do j=jsd,jed ; do I=IsdB,IedB  ! Change any u-face lengths within this loop.
      dy_2 = dx_2 * G%dyCu(I,j)*G%IdxCu(I,j) * cos(pi_180 * G%geoLatCu(I,j))

      if ((abs(G%geoLatCu(I,j)-35.5) < dy_2) .and. (G%geoLonCu(I,j) < -4.5) .and. &
          (G%geoLonCu(I,j) > -6.5)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*12000.0   ! Gibraltar

      if ((abs(G%geoLatCu(I,j)-12.5) < dy_2) .and. (abs(G%geoLonCu(I,j)-43.0) < dx_2)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*10000.0   ! Red Sea

      if ((abs(G%geoLatCu(I,j)-40.5) < dy_2) .and. (abs(G%geoLonCu(I,j)-26.0) < dx_2)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*5000.0   ! Dardanelles

      if ((abs(G%geoLatCu(I,j)-41.5) < dy_2) .and. (abs(G%geoLonCu(I,j)+220.0) < dx_2)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*35000.0   ! Tsugaru strait at 140.0e

      if ((abs(G%geoLatCu(I,j)-45.5) < dy_2) .and. (abs(G%geoLonCu(I,j)+217.5) < 0.9)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*15000.0   ! Betw Hokkaido and Sakhalin at 217&218 = 142e


      ! Greater care needs to be taken in the tripolar region.
      if ((abs(G%geoLatCu(I,j)-80.84) < 0.2) .and. (abs(G%geoLonCu(I,j)+64.9) < 0.8)) &
        G%dy_Cu(I,j) = G%mask2dCu(I,j)*38000.0   ! Smith Sound in Canadian Arch - tripolar region

    enddo ; enddo

    do J=JsdB,JedB ; do i=isd,ied  ! Change any v-face lengths within this loop.
      dy_2 = dx_2 * G%dyCv(i,J)*G%IdxCv(i,J) * cos(pi_180 * G%geoLatCv(i,J))
      if ((abs(G%geoLatCv(i,J)-41.0) < dy_2) .and. (abs(G%geoLonCv(i,J)-28.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*2500.0   ! Bosporus - should be 1000.0 m wide.

      if ((abs(G%geoLatCv(i,J)-13.0) < dy_2) .and. (abs(G%geoLonCv(i,J)-42.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*10000.0   ! Red Sea

      if ((abs(G%geoLatCv(i,J)+2.8) < 0.8) .and. (abs(G%geoLonCv(i,J)+241.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*40000.0   ! Makassar Straits at 241.5 W = 118.5 E

      if ((abs(G%geoLatCv(i,J)-0.56) < 0.5) .and. (abs(G%geoLonCv(i,J)+240.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*80000.0   ! entry to Makassar Straits at 240.5 W = 119.5 E

      if ((abs(G%geoLatCv(i,J)-0.19) < 0.5) .and. (abs(G%geoLonCv(i,J)+230.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*25000.0   ! Channel betw N Guinea and Halmahara 230.5 W = 129.5 E

      if ((abs(G%geoLatCv(i,J)-0.19) < 0.5) .and. (abs(G%geoLonCv(i,J)+229.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*25000.0   ! Channel betw N Guinea and Halmahara 229.5 W = 130.5 E

      if ((abs(G%geoLatCv(i,J)-0.0) < 0.25) .and. (abs(G%geoLonCv(i,J)+228.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*25000.0   ! Channel betw N Guinea and Halmahara 228.5 W = 131.5 E

      if ((abs(G%geoLatCv(i,J)+8.5) < 0.5) .and. (abs(G%geoLonCv(i,J)+244.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*20000.0   ! Lombok Straits at 244.5 W = 115.5 E

      if ((abs(G%geoLatCv(i,J)+8.5) < 0.5) .and. (abs(G%geoLonCv(i,J)+235.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*20000.0   ! Timor Straits at 235.5 W = 124.5 E

      if ((abs(G%geoLatCv(i,J)-52.5) < dy_2) .and. (abs(G%geoLonCv(i,J)+218.5) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*2500.0    ! Russia and Sakhalin Straits at 218.5 W = 141.5 E

      ! Greater care needs to be taken in the tripolar region.
      if ((abs(G%geoLatCv(i,J)-76.8) < 0.06) .and. (abs(G%geoLonCv(i,J)+88.7) < dx_2)) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J)*8400.0    ! Jones Sound in Canadian Arch - tripolar region

    enddo ; enddo
  endif

  ! These checks apply regardless of the chosen option.

  do j=jsd,jed ; do I=IsdB,IedB
    if (G%dy_Cu(I,j) > G%dyCu(I,j)) then
      write(mesg,'("dy_Cu of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4," at lon/lat of ", ES11.4, ES11.4)') &
                   G%dy_Cu(I,j), G%dyCu(I,j), G%dy_Cu(I,j)-G%dyCu(I,j), &
                   G%geoLonCu(I,j), G%geoLatCu(I,j)
      call MOM_error(FATAL,"reset_face_lengths_named "//mesg)
    endif
    G%areaCu(I,j) = G%dxCu(I,j)*G%dy_Cu(I,j)
    G%IareaCu(I,j) = 0.0
    if (G%areaCu(I,j) > 0.0) G%IareaCu(I,j) = G%mask2dCu(I,j) / G%areaCu(I,j)
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (G%dx_Cv(i,J) > G%dxCv(i,J)) then
      write(mesg,'("dx_Cv of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4, " at lon/lat of ", ES11.4, ES11.4)') &
                   G%dx_Cv(i,J), G%dxCv(i,J), G%dx_Cv(i,J)-G%dxCv(i,J), &
                   G%geoLonCv(i,J), G%geoLatCv(i,J)

      call MOM_error(FATAL,"reset_face_lengths_named "//mesg)
    endif
    G%areaCv(i,J) = G%dyCv(i,J)*G%dx_Cv(i,J)
    G%IareaCv(i,J) = 0.0
    if (G%areaCv(i,J) > 0.0) G%IareaCv(i,J) = G%mask2dCv(i,J) / G%areaCv(i,J)
  enddo ; enddo

end subroutine reset_face_lengths_named
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine reset_face_lengths_file(G, param_file)
  type(dyn_horgrid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
!   This subroutine sets the open face lengths at selected points to restrict
! passages to their observed widths.

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
  character(len=40)  :: mod = "reset_face_lengths_file" ! This subroutine's name.
  character(len=256) :: mesg    ! Message for error messages.
  character(len=200) :: filename, chan_file, inputdir ! Strings for file/path
  integer :: i, j, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  ! These checks apply regardless of the chosen option.

  call callTree_enter(trim(mod)//"(), MOM_shared_initialization.F90")

  call get_param(param_file, mod, "CHANNEL_WIDTH_FILE", chan_file, &
                 "The file from which the list of narrowed channels is read.", &
                 default="ocean_geometry.nc")
  call get_param(param_file,  mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(chan_file)
  call log_param(param_file, mod, "INPUTDIR/CHANNEL_WIDTH_FILE", filename)

  if (is_root_pe()) then ; if (.not.file_exists(filename)) &
    call MOM_error(FATAL," reset_face_lengths_file: Unable to open "//&
                           trim(filename))
  endif

  call read_data(filename,"dyCuo",G%dy_Cu,domain=G%Domain%mpp_domain)
  call read_data(filename,"dxCvo",G%dx_Cv,domain=G%Domain%mpp_domain)
  call pass_vector(G%dy_Cu, G%dx_Cv, G%Domain, To_All+SCALAR_PAIR, CGRID_NE)

  do j=jsd,jed ; do I=IsdB,IedB
    if (G%dy_Cu(I,j) > G%dyCu(I,j)) then
      write(mesg,'("dy_Cu of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4," at lon/lat of ", ES11.4, ES11.4)') &
                   G%dy_Cu(I,j), G%dyCu(I,j), G%dy_Cu(I,j)-G%dyCu(I,j), &
                   G%geoLonCu(I,j), G%geoLatCu(I,j)
      call MOM_error(FATAL,"reset_face_lengths_file "//mesg)
    endif
    G%areaCu(I,j) = G%dxCu(I,j)*G%dy_Cu(I,j)
    G%IareaCu(I,j) = 0.0
    if (G%areaCu(I,j) > 0.0) G%IareaCu(I,j) = G%mask2dCu(I,j) / G%areaCu(I,j)
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (G%dx_Cv(i,J) > G%dxCv(i,J)) then
      write(mesg,'("dx_Cv of ",ES11.4," exceeds unrestricted width of ",ES11.4,&
                   &" by ",ES11.4, " at lon/lat of ", ES11.4, ES11.4)') &
                   G%dx_Cv(i,J), G%dxCv(i,J), G%dx_Cv(i,J)-G%dxCv(i,J), &
                   G%geoLonCv(i,J), G%geoLatCv(i,J)

      call MOM_error(FATAL,"reset_face_lengths_file "//mesg)
    endif
    G%areaCv(i,J) = G%dyCv(i,J)*G%dx_Cv(i,J)
    G%IareaCv(i,J) = 0.0
    if (G%areaCv(i,J) > 0.0) G%IareaCv(i,J) = G%mask2dCv(i,J) / G%areaCv(i,J)
  enddo ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine reset_face_lengths_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine reset_face_lengths_list(G, param_file)
  type(dyn_horgrid_type), intent(inout) :: G
  type(param_file_type),  intent(in)    :: param_file
!   This subroutine sets the open face lengths at selected points to restrict
! passages to their observed widths.

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
  character(len=120), pointer, dimension(:) :: lines => NULL()
  character(len=120) :: line
  character(len=200) :: filename, chan_file, inputdir ! Strings for file/path
  character(len=40)  :: mod = "reset_face_lengths_list" ! This subroutine's name.
  real, pointer, dimension(:,:) :: &
    u_lat => NULL(), u_lon => NULL(), v_lat => NULL(), v_lon => NULL()
  real, pointer, dimension(:) :: &
    u_width => NULL(), v_width => NULL()
  real    :: lat, lon     ! The latitude and longitude of a point.
  real    :: lon_p, lon_m ! The longitude of a point shifted by 360 degrees.
  logical :: check_360    ! If true, check for longitudes that are shifted by
                          ! +/- 360 degrees from the specified range of values.
  logical :: found_u, found_v
  logical :: unit_in_use
  integer :: ios, iounit, isu, isv
  integer :: last, num_lines, nl_read, ln, npt, u_pt, v_pt
  integer :: i, j, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call callTree_enter(trim(mod)//"(), MOM_shared_initialization.F90")

  call get_param(param_file, mod, "CHANNEL_LIST_FILE", chan_file, &
                 "The file from which the list of narrowed channels is read.", &
                 default="MOM_channel_list")
  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(chan_file)
  call log_param(param_file, mod, "INPUTDIR/CHANNEL_LIST_FILE", filename)
  call get_param(param_file, mod, "CHANNEL_LIST_360_LON_CHECK", check_360, &
                 "If true, the channel configuration list works for any \n"//&
                 "longitudes in the range of -360 to 360.", default=.true.)

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

  ! Broadcast the number of lines and allocate the required space.
  call broadcast(num_lines, root_PE())
  u_pt = 0 ; v_pt = 0
  if (num_lines > 0) then
    allocate (lines(num_lines))
    if (num_lines > 0) then
      allocate(u_lat(2,num_lines)) ; u_lat(:,:) = -1e34
      allocate(u_lon(2,num_lines)) ; u_lon(:,:) = -1e34
      allocate(u_width(num_lines)) ; u_width(:) = -1e34

      allocate(v_lat(2,num_lines)) ; v_lat(:,:) = -1e34
      allocate(v_lon(2,num_lines)) ; v_lon(:,:) = -1e34
      allocate(v_width(num_lines)) ; v_width(:) = -1e34
    endif

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
        if (is_root_PE()) then
          if (check_360) then
            if ((abs(u_lon(1,u_pt)) > 360.0) .or. (abs(u_lon(2,u_pt)) > 360.0)) &
              call MOM_error(WARNING, "reset_face_lengths_list : Out-of-bounds "//&
                 "u-longitude found when reading line "//trim(line)//" from file "//&
                 trim(filename))
            if ((abs(u_lat(1,u_pt)) > 180.0) .or. (abs(u_lat(2,u_pt)) > 180.0)) &
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
        if (is_root_PE()) then
          if (check_360) then
            if ((abs(v_lon(1,v_pt)) > 360.0) .or. (abs(v_lon(2,v_pt)) > 360.0)) &
              call MOM_error(WARNING, "reset_face_lengths_list : Out-of-bounds "//&
                 "v-longitude found when reading line "//trim(line)//" from file "//&
                 trim(filename))
            if ((abs(v_lat(1,v_pt)) > 180.0) .or. (abs(v_lat(2,v_pt)) > 180.0)) &
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

    deallocate(lines)
  endif

  do j=jsd,jed ; do I=IsdB,IedB
    lat = G%geoLatCu(I,j) ; lon = G%geoLonCu(I,j)
    if (check_360) then ; lon_p = lon+360.0 ; lon_m = lon-360.0
    else ; lon_p = lon ; lon_m = lon ; endif

    do npt=1,u_pt
      if (((lat >= u_lat(1,npt)) .and. (lat <= u_lat(2,npt))) .and. &
          (((lon >= u_lon(1,npt)) .and. (lon <= u_lon(2,npt))) .or. &
           ((lon_p >= u_lon(1,npt)) .and. (lon_p <= u_lon(2,npt))) .or. &
           ((lon_m >= u_lon(1,npt)) .and. (lon_m <= u_lon(2,npt)))) ) &

      G%dy_Cu(I,j) = G%mask2dCu(I,j) * min(G%dyCu(I,j), max(u_width(npt), 0.0))
    enddo

    G%areaCu(I,j) = G%dxCu(I,j)*G%dy_Cu(I,j)
    G%IareaCu(I,j) = 0.0
    if (G%areaCu(I,j) > 0.0) G%IareaCu(I,j) = G%mask2dCu(I,j) / G%areaCu(I,j)
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    lat = G%geoLatCv(i,J) ; lon = G%geoLonCv(i,J)
    if (check_360) then ; lon_p = lon+360.0 ; lon_m = lon-360.0
    else ; lon_p = lon ; lon_m = lon ; endif

    do npt=1,v_pt
      if (((lat >= v_lat(1,npt)) .and. (lat <= v_lat(2,npt))) .and. &
          (((lon >= v_lon(1,npt)) .and. (lon <= v_lon(2,npt))) .or. &
           ((lon_p >= v_lon(1,npt)) .and. (lon_p <= v_lon(2,npt))) .or. &
           ((lon_m >= v_lon(1,npt)) .and. (lon_m <= v_lon(2,npt)))) ) &
        G%dx_Cv(i,J) = G%mask2dCv(i,J) * min(G%dxCv(i,J), max(v_width(npt), 0.0))
    enddo

    G%areaCv(i,J) = G%dyCv(i,J)*G%dx_Cv(i,J)
    G%IareaCv(i,J) = 0.0
    if (G%areaCv(i,J) > 0.0) G%IareaCv(i,J) = G%mask2dCv(i,J) / G%areaCv(i,J)
  enddo ; enddo

  if (num_lines > 0) then
    deallocate(u_lat) ; deallocate(u_lon) ; deallocate(u_width)
    deallocate(v_lat) ; deallocate(v_lon) ; deallocate(v_width)
  endif

  call callTree_leave(trim(mod)//'()')
end subroutine reset_face_lengths_list
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine read_face_length_list(iounit, filename, num_lines, lines)
  integer,                          intent(in)  :: iounit
  character(len=*),                 intent(in)  :: filename
  integer,                          intent(out) :: num_lines
  character(len=120), dimension(:), pointer     :: lines

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
subroutine compute_global_grid_integrals(G)
  type(dyn_horgrid_type), intent(inout) :: G  !< The dynamic horizontal grid
  ! Subroutine to pre-compute global integrals of grid quantities for
  ! later use in reporting diagnostics
  integer :: i,j

  G%areaT_global = 0.0 ; G%IareaT_global = 0.0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    G%areaT_global = G%areaT_global + ( G%areaT(i,j) * G%mask2dT(i,j) )
  enddo ; enddo
  call sum_across_PEs( G%areaT_global )

  if (G%areaT_global == 0.0) &
    call MOM_error(FATAL, "compute_global_grid_integrals: "//&
                    "zero ocean area (check topography?)")

  G%IareaT_global = 1. / G%areaT_global
end subroutine compute_global_grid_integrals
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> Write out a file describing the topography, Coriolis parameter, grid locations
!! and various other fixed fields from the grid.
subroutine write_ocean_geometry_file(G, param_file, directory, geom_file)
  type(dyn_horgrid_type),     intent(inout) :: G         !< The dynamic horizontal grid
  type(param_file_type),      intent(in)    :: param_file !< Parameter file structure
  character(len=*),           intent(in)    :: directory !< The directory into which to place the geometry file.
  character(len=*), optional, intent(in)    :: geom_file !< If present, the name of the geometry file
                                                         !! (otherwise the file is "ocean_geometry")
!   This subroutine writes out a file containing all of the ocean geometry
! and grid data uses by the MOM ocean model.
! Arguments: G - The ocean's grid structure.  Effectively intent in.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      directory - The directory into which to place the file.
  character(len=240) :: filepath
  character(len=40)  :: mod = "write_ocean_geometry_file"
  integer, parameter :: nFlds=23
  type(vardesc) :: vars(nFlds)
  type(fieldtype) :: fields(nFlds)
  integer :: unit
  integer :: file_threading
  integer :: nFlds_used
  integer :: i, j, is, ie, js, je, Isq, Ieq, Jsq, Jeq
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  logical :: multiple_files
  real, dimension(G%isd :G%ied ,G%jsd :G%jed ) :: out_h
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB) :: out_q
  real, dimension(G%IsdB:G%IedB,G%jsd :G%jed ) :: out_u
  real, dimension(G%isd :G%ied ,G%JsdB:G%JedB) :: out_v

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

!   vardesc is a structure defined in MOM_io.F90.  The elements of
! this structure, in order, are:
! (1) the variable name for the NetCDF file
! (2) the variable's long name
! (3) a character indicating the  horizontal grid, which may be '1' (column),
!     'h', 'q', 'u', or 'v', for the corresponding C-grid variable
! (4) a character indicating the vertical grid, which may be 'L' (layer),
!     'i' (interface), or '1' (no vertical location)
! (5) a character indicating the time levels of the field, which may be
!    's' (snap-shot), 'p' (periodic), or '1' (no time variation)
! (6) the variable's units
  vars(1) = var_desc("geolatb","degree","latitude at corner (Bu) points",'q','1','1')
  vars(2) = var_desc("geolonb","degree","longitude at corner (Bu) points",'q','1','1')
  vars(3) = var_desc("geolat","degree", "latitude at tracer (T) points", 'h','1','1')
  vars(4) = var_desc("geolon","degree","longitude at tracer (T) points",'h','1','1')
  vars(5) = var_desc("D","meter","Basin Depth",'h','1','1')
  vars(6) = var_desc("f","second-1","Coriolis Parameter",'q','1','1')
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
  vars(19)= var_desc("wet", "none", "land or ocean?", 'h','1','1')

  vars(20) = var_desc("Dblock_u","meter","Blocked depth at u points",'u','1','1')
  vars(21) = var_desc("Dopen_u","meter","Open depth at u points",'u','1','1')
  vars(22) = var_desc("Dblock_v","meter","Blocked depth at v points",'v','1','1')
  vars(23) = var_desc("Dopen_v","meter","Open depth at v points",'v','1','1')

  nFlds_used = 19 ; if (G%bathymetry_at_vel) nFlds_used = 23

  if (present(geom_file)) then
    filepath = trim(directory) // trim(geom_file)
  else
    filepath = trim(directory) // "ocean_geometry"
  endif

  out_h(:,:) = 0.0
  out_u(:,:) = 0.0
  out_v(:,:) = 0.0
  out_q(:,:) = 0.0

  call get_param(param_file, mod, "PARALLEL_RESTARTFILES", multiple_files, &
                 "If true, each processor writes its own restart file, \n"//&
                 "otherwise a single restart file is generated", &
                 default=.false.)
  file_threading = SINGLE_FILE
  if (multiple_files) file_threading = MULTIPLE

  call create_file(unit, trim(filepath), vars, nFlds_used, fields, &
                   file_threading, dG=G)

  do J=Jsq,Jeq; do I=Isq,Ieq; out_q(I,J) = G%geoLatBu(I,J); enddo; enddo
  call write_field(unit, fields(1), G%Domain%mpp_domain, out_q)
  do J=Jsq,Jeq; do I=Isq,Ieq; out_q(I,J) = G%geoLonBu(I,J); enddo; enddo
  call write_field(unit, fields(2), G%Domain%mpp_domain, out_q)
  call write_field(unit, fields(3), G%Domain%mpp_domain, G%geoLatT)
  call write_field(unit, fields(4), G%Domain%mpp_domain, G%geoLonT)

  call write_field(unit, fields(5), G%Domain%mpp_domain, G%bathyT)
  call write_field(unit, fields(6), G%Domain%mpp_domain, G%CoriolisBu)

  !   I think that all of these copies are holdovers from a much earlier
  ! ancestor code in which many of the metrics were macros that could have
  ! had reduced dimensions, and that they are no longer needed in MOM6. -RWH
  do J=Jsq,Jeq ; do i=is,ie ; out_v(i,J) = G%dxCv(i,J) ; enddo ; enddo
  call write_field(unit, fields(7), G%Domain%mpp_domain, out_v)
  do j=js,je ; do I=Isq,Ieq ; out_u(I,j) = G%dyCu(I,j) ; enddo ; enddo
  call write_field(unit, fields(8), G%Domain%mpp_domain, out_u)

  do j=js,je ; do I=Isq,Ieq ; out_u(I,j) = G%dxCu(I,j) ; enddo ; enddo
  call write_field(unit, fields(9), G%Domain%mpp_domain, out_u)
  do J=Jsq,Jeq ; do i=is,ie ; out_v(i,J) = G%dyCv(i,J) ; enddo ; enddo
  call write_field(unit, fields(10), G%Domain%mpp_domain, out_v)

  do j=js,je ; do i=is,ie ; out_h(i,j) = G%dxT(i,j); enddo; enddo
  call write_field(unit, fields(11), G%Domain%mpp_domain, out_h)
  do j=js,je ; do i=is,ie ; out_h(i,j) = G%dyT(i,j) ; enddo ; enddo
  call write_field(unit, fields(12), G%Domain%mpp_domain, out_h)

  do J=Jsq,Jeq ; do I=Isq,Ieq ; out_q(i,J) = G%dxBu(I,J) ; enddo ; enddo
  call write_field(unit, fields(13), G%Domain%mpp_domain, out_q)
  do J=Jsq,Jeq ; do I=Isq,Ieq ; out_q(I,J) = G%dyBu(I,J) ; enddo ; enddo
  call write_field(unit, fields(14), G%Domain%mpp_domain, out_q)

  do j=js,je ; do i=is,ie ; out_h(i,j) = G%areaT(i,j) ; enddo ; enddo
  call write_field(unit, fields(15), G%Domain%mpp_domain, out_h)
  do J=Jsq,Jeq ; do I=Isq,Ieq ; out_q(I,J) = G%areaBu(I,J) ; enddo ; enddo
  call write_field(unit, fields(16), G%Domain%mpp_domain, out_q)

  call write_field(unit, fields(17), G%Domain%mpp_domain, G%dx_Cv)
  call write_field(unit, fields(18), G%Domain%mpp_domain, G%dy_Cu)
  call write_field(unit, fields(19), G%Domain%mpp_domain, G%mask2dT)

  if (G%bathymetry_at_vel) then
    call write_field(unit, fields(20), G%Domain%mpp_domain, G%Dblock_u)
    call write_field(unit, fields(21), G%Domain%mpp_domain, G%Dopen_u)
    call write_field(unit, fields(22), G%Domain%mpp_domain, G%Dblock_v)
    call write_field(unit, fields(23), G%Domain%mpp_domain, G%Dopen_v)
  endif

  call close_file(unit)

end subroutine write_ocean_geometry_file

end module MOM_shared_initialization
