!> Initializes fixed aspects of the model, such as horizontal grid metrics,
!! topography and Coriolis.
module MOM_fixed_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums, only : hchksum, qchksum, uchksum, vchksum, chksum
use MOM_coms, only : max_across_PEs
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE, AGRID
use MOM_EOS, only : calculate_density, EOS_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_file_parser, only : log_version
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type, isPointInCell
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE, MULTIPLE
use MOM_io, only : slasher, vardesc, write_field, var_desc
use MOM_io, only : EAST_FACE, NORTH_FACE
use MOM_grid_initialize, only : initialize_masks, set_grid_metrics
use MOM_string_functions, only : uppercase
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type, setVerticalGridAxes
use user_initialization, only : user_set_coord, user_initialize_topography
use DOME_initialization, only : DOME_initialize_topography
use ISOMIP_initialization, only : ISOMIP_initialize_topography
use benchmark_initialization, only : benchmark_initialize_topography
use DOME2d_initialization, only : DOME2d_initialize_topography
use sloshing_initialization, only : sloshing_initialize_topography
use seamount_initialization, only : seamount_initialize_topography
use Phillips_initialization, only : Phillips_initialize_topography

use mpp_domains_mod, only  : mpp_global_field, mpp_get_compute_domain
use mpp_mod, only          : mpp_broadcast,mpp_root_pe,mpp_sync,mpp_sync_self

use horiz_interp_mod, only : horiz_interp_new, horiz_interp,horiz_interp_type
use horiz_interp_mod, only : horiz_interp_init, horiz_interp_del

use netcdf

implicit none ; private

#include <MOM_memory.h>

public MOM_initialize_fixed, MOM_initialize_rotation, MOM_initialize_topography

character(len=40) :: mod = "MOM_fixed_initialization" ! This module's name.

contains

! -----------------------------------------------------------------------------
subroutine MOM_initialize_fixed(G, GV, PF, dirs, tv)
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(inout) :: GV   !< Ocean vertical grid structure
  type(param_file_type),   intent(in)    :: PF   !< A structure indicating the open file
                                                 !! to parse for model parameter values.
  type(directories),       intent(in)    :: dirs !< A structure containing relevant paths.
  type(thermo_var_ptrs),   intent(inout) :: tv   !< TO BE DELETED -aja
  ! Local
  character(len=200) :: filename   ! The name of an input file.
  character(len=200) :: filename2  ! The name of an input files.
  character(len = 200) :: inputdir ! The directory where NetCDF input files are.
  character(len=200) :: config
  integer :: write_geom
  logical :: debug, new_sim
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  type(EOS_type), pointer :: eos => NULL() ! TO BE DELETED -aja

  if (associated(tv%eqn_of_state)) eos => tv%eqn_of_state ! TO BE DELETED -aja

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call callTree_enter("MOM_initialize_fixed(), MOM_fixed_initialization.F90")
  call log_version(PF, mod, version)
  call get_param(PF, mod, "DEBUG", debug, default=.false.)

  call get_param(PF, mod, "INPUTDIR", inputdir, &
         "The directory in which input files are found.", default=".")
  inputdir = slasher(inputdir)

! Set up the parameters of the physical domain (i.e. the grid), G
  call set_grid_metrics(G, PF)

! Set up the bottom depth, G%bathyT either analytically or from file
  call MOM_initialize_topography(G%bathyT, G%max_depth, G, PF)

! ====================================================================
!    Initialize fields that are time invariant - metrics, topography,
!  masks, vertical coordinate, Coriolis parameter.
! ====================================================================

! Set-up the layer densities, GV%Rlay, and reduced gravities, GV%g_prime.
  call get_param(PF, mod, "COORD_CONFIG", config, &
                 "This specifies how layers are to be defined: \n"//&
                 " \t file - read coordinate information from the file \n"//&
                 " \t\t specified by (COORD_FILE).\n"//&
                 " \t linear - linear based on interfaces not layers \n"//&
                 " \t layer_ref - linear based on layer densities \n"//&
                 " \t ts_ref - use reference temperature and salinity \n"//&
                 " \t ts_range - use range of temperature and salinity \n"//&
                 " \t\t (T_REF and S_REF) to determine surface density \n"//&
                 " \t\t and GINT calculate internal densities. \n"//&
                 " \t gprime - use reference density (RHO_0) for surface \n"//&
                 " \t\t density and GINT calculate internal densities. \n"//&
                 " \t ts_profile - use temperature and salinity profiles \n"//&
                 " \t\t (read from COORD_FILE) to set layer densities. \n"//&
                 " \t USER - call a user modified routine.", &
                 fail_if_missing=.true.)
  select case ( trim(config) )
    case ("gprime")
      call set_coord_from_gprime(GV%Rlay, GV%g_prime, G, GV, PF)
    case ("layer_ref")
      call set_coord_from_layer_density(GV%Rlay, GV%g_prime, G, GV, PF)
    case ("linear")
      call set_coord_linear(GV%Rlay, GV%g_prime, G, GV, PF)
    case ("ts_ref")
      call set_coord_from_ts_ref(GV%Rlay, GV%g_prime, G, GV, PF, eos, tv%P_Ref)
    case ("ts_profile")
      call set_coord_from_TS_profile(GV%Rlay, GV%g_prime, G, GV, PF, eos, tv%P_Ref)
    case ("ts_range")
      call set_coord_from_TS_range(GV%Rlay, GV%g_prime, G, GV, PF, eos, tv%P_Ref)
    case ("file")
      call set_coord_from_file(GV%Rlay, GV%g_prime, G, GV, PF)
    case ("USER")
      call user_set_coord(GV%Rlay, GV%g_prime, G, PF, eos)
    case ("none")
    case default ; call MOM_error(FATAL,"MOM_initialize_fixed: "// &
      "Unrecognized coordinate setup"//trim(config))
  end select
  if (debug) call chksum(GV%Rlay, "MOM_initialize_fixed: Rlay ", 1, nz)
  if (debug) call chksum(GV%g_prime, "MOM_initialize_fixed: g_prime ", 1, nz)
  call setVerticalGridAxes( GV%Rlay, GV )


!    This call sets seamasks that prohibit flow over any point with  !
!  a bottom that is shallower than min_depth from PF.                !
  call initialize_masks(G, PF)
  if (debug) then
    call hchksum(G%bathyT, 'MOM_initialize_fixed: depth ', G, haloshift=1)
    call hchksum(G%mask2dT, 'MOM_initialize_fixed: mask2dT ', G)
    call uchksum(G%mask2dCu, 'MOM_initialize_fixed: mask2dCu ', G)
    call vchksum(G%mask2dCv, 'MOM_initialize_fixed: mask2dCv ', G)
    call qchksum(G%mask2dBu, 'MOM_initialize_fixed: mask2dBu ', G)
  endif

! Modulate geometric scales according to geography.
  call get_param(PF, mod, "CHANNEL_CONFIG", config, &
                 "A parameter that determines which set of channels are \n"//&
                 "restricted to specific  widths.  Options are:\n"//&
                 " \t none - All channels have the grid width.\n"//&
                 " \t global_1deg - Sets 16 specific channels appropriate \n"//&
                 " \t\t for a 1-degree model, as used in CM2G.\n"//&
                 " \t list - Read the channel locations and widths from a \n"//&
                 " \t\t text file, like MOM_channel_list in the MOM_SIS \n"//&
                 " \t\t test case.\n"//&
                 " \t file - Read open face widths everywhere from a \n"//&
                 " \t\t NetCDF file on the model grid.", &
                 default="none")
  select case ( trim(config) )
    case ("none")
    case ("list") ; call reset_face_lengths_list(G, PF)
    case ("file") ; call reset_face_lengths_file(G, PF)
    case ("global_1deg") ; call reset_face_lengths_named(G, PF, trim(config))
    case default ; call MOM_error(FATAL, "MOM_initialize_fixed: "// &
      "Unrecognized channel configuration "//trim(config))
  end select

!   This call sets the topography at velocity points.
  if (G%bathymetry_at_vel) then
    call get_param(PF, mod, "VELOCITY_DEPTH_CONFIG", config, &
                   "A string that determines how the topography is set at \n"//&
                   "velocity points. This may be 'min' or 'max'.", &
                   default="max")
    select case ( trim(config) )
      case ("max") ; call set_velocity_depth_max(G)
      case ("min") ; call set_velocity_depth_min(G)
      case default ; call MOM_error(FATAL, "MOM_initialize_fixed: "// &
        "Unrecognized velocity depth configuration "//trim(config))
    end select
  endif

!    Calculate the value of the Coriolis parameter at the latitude   !
!  of the q grid points, in s-1.
  call MOM_initialize_rotation(G%CoriolisBu, G, PF)
!   Calculate the components of grad f (beta)
  call MOM_calculate_grad_Coriolis(G%dF_dx, G%dF_dy, G)
  if (debug) then
    call qchksum(G%CoriolisBu, "MOM_initialize_fixed: f ", G)
    call hchksum(G%dF_dx, "MOM_initialize_fixed: dF_dx ", G)
    call hchksum(G%dF_dy, "MOM_initialize_fixed: dF_dy ", G)
  endif

! Compute global integrals of grid values for later use in scalar diagnostics !
  call compute_global_grid_integrals(G)

! Write out all of the grid data used by this run.
  call get_param(PF, mod, "WRITE_GEOM", write_geom, &
                 "If =0, never write the geometry and vertical grid files.\n"//&
                 "If =1, write the geometry and vertical grid files only for\n"//&
                 "a new simulation. If =2, always write the geometry and\n"//&
                 "vertical grid files. Other values are invalid.", default=1)
  if (write_geom<0 .or. write_geom>2) call MOM_error(FATAL,"MOM_initialize_fixed: "//&
         "WRITE_GEOM must be equal to 0, 1 or 2.")
  new_sim = .false.
  if ((dirs%input_filename(1:1)=='n') .and. (LEN_TRIM(dirs%input_filename)==1)) new_sim = .true.
  if ((write_geom==1 .and. new_sim) .or. write_geom==2) then
    call write_ocean_geometry_file(G, PF, dirs%output_directory)
    call write_vertgrid_file(GV, G, PF, dirs%output_directory)
  endif

  call callTree_leave('MOM_initialize_fixed()')

end subroutine MOM_initialize_fixed
! -----------------------------------------------------------------------------

! The set_coord routines deal with initializing aspects of the vertical grid.
! -----------------------------------------------------------------------------
subroutine set_coord_from_gprime(Rlay, g_prime, G, GV, param_file)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(verticalGrid_type), intent(in)  :: GV
  type(param_file_type), intent(in)  :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: g_int   ! Reduced gravities across the internal interfaces, in m s-2.
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  character(len=40)  :: mod = "set_coord_from_gprime" ! This subroutine's name.
  integer :: k, nz
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

  call get_param(param_file, mod, "GFS" , g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call get_param(param_file, mod, "GINT", g_int, &
                 "The reduced gravity across internal interfaces.", &
                 units="m s-2", fail_if_missing=.true.)

  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = g_int ; enddo
  Rlay(1) = GV%Rho0
  do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(GV%Rho0/G%g_Earth) ; enddo

  call callTree_leave(trim(mod)//'()')

end subroutine set_coord_from_gprime
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_layer_density(Rlay, g_prime, G, GV, param_file)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(verticalGrid_type), intent(in)  :: GV
  type(param_file_type), intent(in)  :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  real :: Rlay_Ref! The surface layer's target density, in kg m-3.
  real :: RLay_range ! The range of densities, in kg m-3.
  character(len=40)  :: mod = "set_coord_from_layer_density" ! This subroutine's name.
  integer :: k, nz
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call get_param(param_file, mod, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for layer 1.", &
                 units="kg m-3", default=GV%Rho0)
  call get_param(param_file, mod, "DENSITY_RANGE", Rlay_range, &
                 "The range of reference potential densities in the layers.", &
                 units="kg m-3", default=2.0)

  g_prime(1) = g_fs
  Rlay(1) = Rlay_Ref
  do k=2,nz
     Rlay(k) = Rlay(k-1) + RLay_range/(real(nz-1))
  enddo
!    These statements set the interface reduced gravities.           !
  do k=2,nz
     g_prime(k) = (G%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1))
  enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_from_layer_density
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_TS_ref(Rlay, g_prime, G, GV, param_file, eqn_of_state, &
                                 P_Ref)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(verticalGrid_type), intent(in)  :: GV
  type(param_file_type), intent(in)  :: param_file
  type(EOS_type),        pointer     :: eqn_of_state
  real,                  intent(in)  :: P_Ref
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer selecting the equation of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: T_ref   ! Reference temperature
  real :: S_ref   ! Reference salinity
  real :: g_int   ! Reduced gravities across the internal interfaces, in m s-2.
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  character(len=40)  :: mod = "set_coord_from_TS_ref" ! This subroutine's name.
  integer :: k, nz
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

  call get_param(param_file, mod, "T_REF", T_Ref, &
                 "The initial temperature of the lightest layer.", units="degC", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "S_REF", S_Ref, &
                 "The initial salinities.", units="PSU", default=35.0)
  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call get_param(param_file, mod, "GINT", g_int, &
                 "The reduced gravity across internal interfaces.", &
                 units="m s-2", fail_if_missing=.true.)
                                      !
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = g_int ; enddo

!    The uppermost layer's density is set here.  Subsequent layers'  !
!  densities are determined from this value and the g values.        !
!        T0 = 28.228 ; S0 = 34.5848 ; Pref = P_Ref
  call calculate_density(T_ref, S_ref, P_ref, Rlay(1), eqn_of_state)

!    These statements set the layer densities.                       !
  do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(GV%Rho0/G%g_Earth) ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_from_TS_ref
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_TS_profile(Rlay, g_prime, G, GV, param_file, &
                                     eqn_of_state, P_Ref)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(verticalGrid_type), intent(in)  :: GV
  type(param_file_type), intent(in)  :: param_file
  type(EOS_type),        pointer     :: eqn_of_state
  real,                  intent(in)  :: P_Ref
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer that selects equation of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real, dimension(SZK_(G)) :: T0, S0,  Pref
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  integer :: k, nz
  character(len=40)  :: mod = "set_coord_from_TS_profile" ! This subroutine's name.
  character(len=200) :: filename, coord_file, inputdir ! Strings for file/path
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call get_param(param_file, mod, "COORD_FILE", coord_file, &
                 "The file from which the coordinate temperatures and \n"//&
                 "salnities are read.", fail_if_missing=.true.)

  call get_param(param_file,  mod, "INPUTDIR", inputdir, default=".")
  filename = trim(slasher(inputdir))//trim(coord_file)
  call log_param(param_file, mod, "INPUTDIR/COORD_FILE", filename)

  call read_data(filename,"PTEMP",T0(:),domain=G%Domain%mpp_domain)
  call read_data(filename,"SALT",S0(:),domain=G%Domain%mpp_domain)

  if (.not.file_exists(filename)) call MOM_error(FATAL, &
      " set_coord_from_TS_profile: Unable to open " //trim(filename))
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=1,nz ; Pref(k) = P_ref ; enddo
  call calculate_density(T0, S0, Pref, Rlay, 1,nz,eqn_of_state)
  do k=2,nz; g_prime(k) = (G%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1)); enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_from_TS_profile
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_TS_range(Rlay, g_prime, G, GV, param_file, &
                                     eqn_of_state, P_Ref)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(verticalGrid_type), intent(in)  :: GV
  type(param_file_type), intent(in)  :: param_file
  type(EOS_type),        pointer     :: eqn_of_state
  real,                  intent(in)  :: P_Ref
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      eqn_of_state - integer that selects equation of state
!  (in)      P_Ref - The coordinate-density reference pressure in Pa.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real, dimension(SZK_(G)) :: T0, S0,  Pref
  real :: S_Ref, S_Light, S_Dense ! Salnity range parameters in PSU.
  real :: T_Ref, T_Light, T_Dense ! Temperature range parameters in dec C.
  real :: res_rat ! The ratio of density space resolution in the denser part
                  ! of the range to that in the lighter part of the range.
                  ! Setting this greater than 1 increases the resolution for
                  ! the denser water.
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  real :: a1, frac_dense, k_frac
  integer :: k, nz, k_light
  character(len=40)  :: mod = "set_coord_from_TS_range" ! This subroutine's name.
  character(len=200) :: filename, coord_file, inputdir ! Strings for file/path
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

  call get_param(param_file, mod, "T_REF", T_Ref, &
                 "The default initial temperatures.", units="degC", default=10.0)
  call get_param(param_file, mod, "TS_RANGE_T_LIGHT", T_Light, &
                 "The initial temperature of the lightest layer when \n"//&
                 "COORD_CONFIG is set to ts_range.", units="degC", default=T_Ref)
  call get_param(param_file, mod, "TS_RANGE_T_DENSE", T_Dense, &
                 "The initial temperature of the densest layer when \n"//&
                 "COORD_CONFIG is set to ts_range.", units="degC", default=T_Ref)

  call get_param(param_file, mod, "S_REF", S_Ref, &
                 "The default initial salinities.", units="PSU", default=35.0)
  call get_param(param_file, mod, "TS_RANGE_S_LIGHT", S_Light, &
                 "The initial lightest salinities when COORD_CONFIG \n"//&
                 "is set to ts_range.", default = S_Ref, units="PSU")
  call get_param(param_file, mod, "TS_RANGE_S_DENSE", S_Dense, &
                 "The initial densest salinities when COORD_CONFIG \n"//&
                 "is set to ts_range.", default = S_Ref, units="PSU")
 
  call get_param(param_file, mod, "TS_RANGE_RESOLN_RATIO", res_rat, &
                 "The ratio of density space resolution in the densest \n"//&
                 "part of the range to that in the lightest part of the \n"//&
                 "range when COORD_CONFIG is set to ts_range. Values \n"//&
                 "greater than 1 increase the resolution of the denser water.",&
                 default=1.0, units="nondim")

  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)

  k_light = GV%nk_rho_varies + 1

  ! Set T0(k) to range from T_LIGHT to T_DENSE, and simliarly for S0(k).
  T0(k_light) = T_light ; S0(k_light) = S_light
  a1 = 2.0 * res_rat / (1.0 + res_rat)
  do k=k_light+1,nz
    k_frac = real(k-k_light)/real(nz-k_light)
    frac_dense = a1 * k_frac + (1.0 - a1) * k_frac**2
    T0(k) = frac_dense * (T_Dense - T_Light) + T_Light
    S0(k) = frac_dense * (S_Dense - S_Light) + S_Light
  enddo

  g_prime(1) = g_fs
  do k=1,nz ; Pref(k) = P_ref ; enddo
  call calculate_density(T0, S0, Pref, Rlay, k_light,nz-k_light+1,eqn_of_state)
  ! Extrapolate target densities for the variable density mixed and buffer layers.
  do k=k_light-1,1,-1
    Rlay(k) = 2.0*Rlay(k+1) -  Rlay(k+2)
  enddo
  do k=2,nz; g_prime(k) = (G%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1)); enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_from_TS_range
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_from_file(Rlay, g_prime, G, GV, param_file)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(verticalGrid_type), intent(in)  :: GV
  type(param_file_type), intent(in)  :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  !
! reduced gravities (g).                                             !
  real :: g_fs    ! Reduced gravity across the free surface, in m s-2.
  integer :: k, nz
  character(len=40)  :: mod = "set_coord_from_file" ! This subroutine's name.
  character(len=40)  :: coord_var
  character(len=200) :: filename,coord_file,inputdir ! Strings for file/path
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)
  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mod, "COORD_FILE", coord_file, &
                 "The file from which the coordinate densities are read.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mod, "COORD_VAR", coord_var, &
                 "The variable in COORD_FILE that is to be used for the \n"//&
                 "coordinate densities.", default="Layer")
  filename = trim(inputdir)//trim(coord_file)
  call log_param(param_file, mod, "INPUTDIR/COORD_FILE", filename)
  if (.not.file_exists(filename)) call MOM_error(FATAL, &
      " set_coord_from_file: Unable to open "//trim(filename))

  call read_axis_data(filename, coord_var, Rlay)
  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = (G%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1)) ; enddo
  do k=1,nz ; if (g_prime(k) <= 0.0) then
    call MOM_error(FATAL, "MOM_initialization set_coord_from_file: "//&
       "Zero or negative g_primes read from variable "//"Layer"//" in file "//&
       trim(filename))
  endif ; enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_from_file
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine set_coord_linear(Rlay, g_prime, G, GV, param_file)
  real, dimension(:),    intent(out) :: Rlay, g_prime
  type(ocean_grid_type), intent(in)  :: G
  type(verticalGrid_type), intent(in)  :: GV
  type(param_file_type), intent(in)  :: param_file
! Arguments: Rlay - the layers' target coordinate values (potential density).
!  (out)     g_prime - the reduced gravity across the interfaces, in m s-2.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This subroutine sets the layer densities (Rlay) and the interface  
! reduced gravities (g) according to a linear profile starting at a
! reference surface layer density and spanning a range of densities 
! to the bottom defined by the parameter RLAY_RANGE
! (defaulting to 2.0 if not defined)
  character(len=40)  :: mod = "set_coord_linear" ! This subroutine
  real :: Rlay_ref, Rlay_range, g_fs
  integer :: k, nz
  nz = G%ke

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

  call get_param(param_file, mod, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for the surface \n"// &
                 "interface.", units="kg m-3", default=GV%Rho0)
  call get_param(param_file, mod, "DENSITY_RANGE", Rlay_range, &
                 "The range of reference potential densities across \n"// &
                 "all interfaces.", units="kg m-3", default=2.0)
  call get_param(param_file, mod, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=G%g_Earth)

  ! This following sets the target layer densities such that a the
  ! surface interface has density Rlay_ref and the bottom
  ! is Rlay_range larger
  do k=1,nz
     Rlay(k) = Rlay_Ref + RLay_range*((real(k)-0.5)/real(nz))
  enddo
  ! These statements set the interface reduced gravities.
  g_prime(1) = g_fs
  do k=2,nz 
     g_prime(k) = (G%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1))
  enddo

  call callTree_leave(trim(mod)//'()')
end subroutine set_coord_linear

! The remainder of the files deals with initializing the horizontal grid.
! -----------------------------------------------------------------------------

subroutine MOM_initialize_rotation(f, G, PF)
  type(ocean_grid_type),                        intent(in)  :: G
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), intent(out) :: f
  type(param_file_type),                        intent(in)  :: PF
! Arguments: f  - the Coriolis parameter in s-1. Intent out.
!  (in)      G  - The ocean's grid structure.
!  (in)      PF - A structure indicating the open file to parse for
!                         model parameter values.

!   This subroutine makes the appropriate call to set up the Coriolis parameter.
! This is a separate subroutine so that it can be made public and shared with
! the ice-sheet code or other components.
! Set up the Coriolis parameter, f, either analytically or from file.
  character(len=40)  :: mod = "MOM_initialize_rotation" ! This subroutine's name.
  character(len=200) :: config

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")
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
  real, dimension(NIMEM_,NJMEM_), intent(out)   :: dF_dx !< x-component of grad f
  real, dimension(NIMEM_,NJMEM_), intent(out)   :: dF_dy !< y-component of grad f
  type(ocean_grid_type),          intent(inout) :: G !< Grid type
  ! Local variables
  integer :: i,j
  real :: f1, f2
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

subroutine MOM_initialize_topography(D, max_depth, G, PF)
  real, dimension(NIMEM_,NJMEM_), intent(out) :: D
  real,                           intent(out) :: max_depth
  type(ocean_grid_type),          intent(in)  :: G
  type(param_file_type),          intent(in)  :: PF
! Arguments: D  - the bottom depth in m. Intent out.
!  (in)      G  - The ocean's grid structure.
!  (in)      PF - A structure indicating the open file to parse for
!                         model parameter values.

!  This subroutine makes the appropriate call to set up the bottom depth.
!  This is a separate subroutine so that it can be made public and shared with
!  the ice-sheet code or other components.
! Set up the bottom depth, G%bathyT either analytically or from file
  character(len=40)  :: mod = "MOM_initialize_topography" ! This subroutine's name.
  character(len=200) :: config

  call get_param(PF, mod, "TOPO_CONFIG", config, &
                 "This specifies how bathymetry is specified: \n"//&
                 " \t file - read bathymetric information from the file \n"//&
                 " \t\t specified by (TOPO_FILE).\n"//&
                 " \t flat - flat bottom set to MAXIMUM_DEPTH. \n"//&
                 " \t bowl - an analytically specified bowl-shaped basin \n"//&
                 " \t\t ranging between MAXIMUM_DEPTH and MINIMUM_DEPTH. \n"//&
                 " \t spoon - a similar shape to 'bowl', but with an vertical \n"//&
                 " \t\t wall at the southern face. \n"//&
                 " \t halfpipe - a zonally uniform channel with a half-sine \n"//&
                 " \t\t profile in the meridional direction. \n"//&
                 " \t benchmark - use the benchmark test case topography. \n"//&
                 " \t DOME - use a slope and channel configuration for the \n"//&
                 " \t\t DOME sill-overflow test case. \n"//&
                 " \t ISOMIP - use a slope and channel configuration for the \n"//&
                 " \t\t ISOMIP test case. \n"//&
                 " \t DOME2D - use a shelf and slope configuration for the \n"//&
                 " \t\t DOME2D gravity current/overflow test case. \n"//&
                 " \t seamount - Gaussian bump for spontaneous motion test case.\n"//&
                 " \t Phillips - ACC-like idealized topography used in the Phillips config.\n"//&
                 " \t USER - call a user modified routine.", &
                 fail_if_missing=.true.)
  max_depth = -1.e9; call read_param(PF, "MAXIMUM_DEPTH", max_depth)
  select case ( trim(config) )
    case ("file");      call initialize_topography_from_file(D, G, PF)
    case ("flat");      call initialize_topography_named(D, G, PF, config, max_depth)
    case ("spoon");     call initialize_topography_named(D, G, PF, config, max_depth)
    case ("bowl");      call initialize_topography_named(D, G, PF, config, max_depth)
    case ("halfpipe");  call initialize_topography_named(D, G, PF, config, max_depth)
    case ("DOME");      call DOME_initialize_topography(D, G, PF, max_depth)
    case ("ISOMIP");      call ISOMIP_initialize_topography(D, G, PF, max_depth)
    case ("benchmark"); call benchmark_initialize_topography(D, G, PF, max_depth)
    case ("DOME2D");    call DOME2d_initialize_topography(D, G, PF, max_depth)
    case ("sloshing");  call sloshing_initialize_topography(D, G, PF, max_depth)
    case ("seamount");  call seamount_initialize_topography(D, G, PF, max_depth)
    case ("Phillips");  call Phillips_initialize_topography(D, G, PF)
    case ("USER");      call user_initialize_topography(D, G, PF)
    case default ;      call MOM_error(FATAL,"MOM_initialize_topography: "// &
      "Unrecognized topography setup '"//trim(config)//"'")
  end select
  if (max_depth>0.) then
    call log_param(PF, mod, "MAXIMUM_DEPTH", max_depth, &
                   "The maximum depth of the ocean.", units="m")
  else
    max_depth = diagnoseMaximumDepth(D,G)
    call log_param(PF, mod, "!MAXIMUM_DEPTH", max_depth, &
                   "The (diagnosed) maximum depth of the ocean.", units="m")
  endif
  if (trim(config) .ne. "DOME") then
    call limit_topography(D, G, PF, max_depth)
  endif
  
end subroutine MOM_initialize_topography

! -----------------------------------------------------------------------------
function diagnoseMaximumDepth(D,G)
  real, dimension(NIMEM_,NJMEM_), intent(in) :: D
  type(ocean_grid_type),          intent(in) :: G
  real :: diagnoseMaximumDepth
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
! -----------------------------------------------------------------------------

!> Read gridded depths from file
subroutine initialize_topography_from_file(D, G, param_file)
  real, dimension(NIMEM_,NJMEM_), intent(out) :: D !< Bottom depth (positive, in m)
  type(ocean_grid_type),          intent(in)  :: G !< Grid structure
  type(param_file_type),          intent(in)  :: param_file !< Parameter file structure
  ! Local variables
  character(len=200) :: filename, topo_file, inputdir ! Strings for file/path
  character(len=200) :: topo_varname                  ! Variable name in file
  character(len=40)  :: mod = "initialize_topography_from_file" ! This subroutine's name.

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

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
  real, dimension(NIMEM_,NJMEM_), intent(inout) :: D !< Bottom depth (positive, in m)
  type(ocean_grid_type),          intent(in)    :: G !< Grid structure
  type(param_file_type),          intent(in)    :: param_file !< Parameter file structure
  ! Local variables
  character(len=200) :: topo_edits_file, inputdir ! Strings for file/path
  character(len=40)  :: mod = "apply_topography_edits_from_file" ! This subroutine's name.
  integer :: n_edits, n, ashape(5), i, j, ncid, id, ncstatus, iid, jid, zid
  integer, dimension(:), allocatable :: ig, jg
  real, dimension(:), allocatable :: new_depth

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

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

! -----------------------------------------------------------------------------
subroutine initialize_topography_named(D, G, param_file, topog_config, max_depth)
  real, dimension(NIMEM_,NJMEM_), intent(out) :: D
  type(ocean_grid_type),          intent(in)  :: G
  type(param_file_type),          intent(in)  :: param_file
  character(len=*),               intent(in)  :: topog_config
  real,                           intent(in)  :: max_depth
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

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")
  call MOM_mesg("  MOM_fixed_initialization.F90, initialize_topography_named: "//&
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
subroutine limit_topography(D, G, param_file, max_depth)
  real, dimension(NIMEM_,NJMEM_), intent(inout) :: D
  type(ocean_grid_type),          intent(in)    :: G
  type(param_file_type),          intent(in)    :: param_file
  real,                           intent(in)    :: max_depth
! Arguments: D          - the bottom depth in m. Intent in/out.
!  (in)      G          - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      max_depth  - The maximum depth in m.

! This subroutine ensures that    min_depth < D(x,y) < max_depth
  integer :: i, j
  character(len=40)  :: mod = "limit_topography" ! This subroutine's name.
  real :: min_depth, mask_depth

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

  call get_param(param_file, "MOM_grid_init initialize_masks", "MINIMUM_DEPTH", min_depth, &
                 "If MASKING_DEPTH is unspecified, then anything shallower than\n"//&
                 "MINIMUM_DEPTH is assumed to be land and all fluxes are masked out.\n"//&
                 "If MASKING_DEPTH is specified, then all depths shallower than\n"//&
                 "MINIMUM_DEPTH but depper than MASKING_DEPTH are rounded to MINIMUM_DEPTH.", &
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
  type(ocean_grid_type),                        intent(in)  :: G
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), intent(out) :: f
  type(param_file_type),                        intent(in)  :: param_file
! Arguments: f          - Coriolis parameter (vertical component) in s^-1
!     (in)   G          - grid type
!     (in)   param_file - parameter file type

! This subroutine sets up the Coriolis parameter for a sphere
  character(len=30) :: mod = "set_rotation_planetary" ! This subroutine's name.
  integer :: I, J
  real    :: PI, omega

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

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
  type(ocean_grid_type),                        intent(in)  :: G
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

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

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
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine reset_face_lengths_named(G, param_file, name)
  type(ocean_grid_type), intent(inout) :: G
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

      if ((abs(G%geoLatCu(i,j)-40.5) < dy_2) .and. (abs(G%geoLonCu(i,j)-26.0) < dx_2)) &
        G%dy_Cu(i,j) = G%mask2dCu(i,j)*5000.0   ! Dardanelles

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
  type(ocean_grid_type), intent(inout) :: G
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

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

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
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
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

  call callTree_enter(trim(mod)//"(), MOM_fixed_initialization.F90")

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
subroutine set_velocity_depth_max(G)
  type(ocean_grid_type), intent(inout) :: G
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
subroutine compute_global_grid_integrals(G)
  type(ocean_grid_type), intent(inout) :: G
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
subroutine set_velocity_depth_min(G)
  type(ocean_grid_type), intent(inout) :: G
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
subroutine write_ocean_geometry_file(G, param_file, directory)
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  character(len=*),      intent(in)    :: directory
!   This subroutine writes out a file containing all of the ocean geometry
! and grid data uses by the MOM ocean model.
! Arguments: G - The ocean's grid structure.  Effectively intent in.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      directory - The directory into which to place the file.
  character(len=120) :: filepath
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
  real :: out_h(SZI_(G),SZJ_(G))
  real :: out_u(SZIB_(G),SZJ_(G))
  real :: out_v(SZI_(G),SZJB_(G))
  real :: out_q(SZIB_(G),SZJB_(G))
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

  filepath = trim(directory) // "ocean_geometry"

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

  call create_file(unit, trim(filepath), vars, nFlds_used, G, fields, file_threading)

  do J=Jsq,Jeq; do I=Isq,Ieq; out_q(I,J) = G%geoLatBu(I,J); enddo; enddo
  call write_field(unit, fields(1), G%Domain%mpp_domain, out_q)
  do J=Jsq,Jeq; do I=Isq,Ieq; out_q(I,J) = G%geoLonBu(I,J); enddo; enddo
  call write_field(unit, fields(2), G%Domain%mpp_domain, out_q)
  call write_field(unit, fields(3), G%Domain%mpp_domain, G%geoLatT)
  call write_field(unit, fields(4), G%Domain%mpp_domain, G%geoLonT)

  call write_field(unit, fields(5), G%Domain%mpp_domain, G%bathyT)
  call write_field(unit, fields(6), G%Domain%mpp_domain, G%CoriolisBu)

  do J=Jsq,Jeq; do i=is,ie; out_v(i,J) = G%dxCv(i,J); enddo; enddo
  call write_field(unit, fields(7), G%Domain%mpp_domain, out_v)
  do j=js,je; do I=Isq,Ieq; out_u(I,j) = G%dyCu(I,j); enddo; enddo
  call write_field(unit, fields(8), G%Domain%mpp_domain, out_u)

  do J=Jsq,Jeq; do i=is,ie; out_u(i,J) = G%dxCu(i,J); enddo; enddo
  call write_field(unit, fields(9), G%Domain%mpp_domain, out_u)
  do j=js,je; do I=Isq,Ieq; out_v(I,j) = G%dyCv(I,j); enddo; enddo
  call write_field(unit, fields(10), G%Domain%mpp_domain, out_v)

  do J=Jsq,Jeq; do i=is,ie; out_h(i,J) = G%dxT(i,J); enddo; enddo
  call write_field(unit, fields(11), G%Domain%mpp_domain, out_h)
  do j=js,je; do I=Isq,Ieq; out_h(I,j) = G%dyT(I,j); enddo; enddo
  call write_field(unit, fields(12), G%Domain%mpp_domain, out_h)

  do J=Jsq,Jeq; do i=is,ie; out_q(i,J) = G%dxBu(i,J); enddo; enddo
  call write_field(unit, fields(13), G%Domain%mpp_domain, out_q)
  do j=js,je; do I=Isq,Ieq; out_q(I,j) = G%dyBu(I,j); enddo; enddo
  call write_field(unit, fields(14), G%Domain%mpp_domain, out_q)

  do j=js,je; do i=is,ie; out_h(i,j) = G%areaT(i,j); enddo; enddo
  call write_field(unit, fields(15), G%Domain%mpp_domain, out_h)
  do j=js,je; do i=is,ie; out_q(i,j) = G%areaBu(i,j); enddo; enddo
  call write_field(unit, fields(16), G%Domain%mpp_domain, out_q)

!  do J=Jsq,Jeq; do i=is,ie; out_v(i,J) = G%dx_Cv(i,J); enddo; enddo
  call write_field(unit, fields(17), G%Domain%mpp_domain, G%dx_Cv)
!  do j=js,je; do I=Isq,Ieq; out_u(I,j) = G%dy_Cu(I,j); enddo; enddo
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
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
subroutine write_vertgrid_file(GV, G, param_file, directory)
  type(verticalGrid_type), intent(in)  :: GV
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  character(len=*),      intent(in)    :: directory
!   This subroutine writes out a file containing any available data related
! to the vertical grid used by the MOM ocean model.
! Arguments: Gv - The container for the vertical grid data.
!  (in)      G - The ocean's grid structure.  Effectively intent in.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      directory - The directory into which to place the file.
  character(len=120) :: filepath
  type(vardesc) :: vars(2)
  type(fieldtype) :: fields(2)
  integer :: unit

  filepath = trim(directory) // trim("Vertical_coordinate")

  vars(1) = var_desc("R","kilogram meter-3","Target Potential Density",'1','L','1')
  vars(2) = var_desc("g","meter second-2","Reduced gravity",'1','L','1')

  call create_file(unit, trim(filepath), vars, 2, G, fields, SINGLE_FILE, GV=GV)

  call write_field(unit, fields(1), GV%Rlay)
  call write_field(unit, fields(2), GV%g_prime)

  call close_file(unit)

end subroutine write_vertgrid_file
! -----------------------------------------------------------------------------

end module MOM_fixed_initialization
