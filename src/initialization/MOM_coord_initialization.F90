!> Initializes fixed aspects of the related to its vertical coordinate.
module MOM_coord_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging, only : chksum
use MOM_EOS, only : calculate_density, EOS_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_file_parser, only : log_version
use MOM_io, only : close_file, create_file, fieldtype, file_exists
use MOM_io, only : open_file, MOM_read_data, read_axis_data, SINGLE_FILE, MULTIPLE
use MOM_io, only : slasher, vardesc, write_field, var_desc
use MOM_string_functions, only : uppercase
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type, setVerticalGridAxes
use user_initialization, only : user_set_coord
use BFB_initialization, only : BFB_set_coord

use netcdf

implicit none ; private

public MOM_initialize_coord

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

character(len=40) :: mdl = "MOM_coord_initialization" !< This module's name.

contains

!> MOM_initialize_coord sets up time-invariant quantities related to MOM6's
!!   vertical coordinate.
subroutine MOM_initialize_coord(GV, US, PF, write_geom, output_dir, tv, max_depth)
  type(verticalGrid_type), intent(inout) :: GV         !< Ocean vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: PF         !< A structure indicating the open file
                                                       !! to parse for model parameter values.
  logical,                 intent(in)    :: write_geom !< If true, write grid geometry files.
  character(len=*),        intent(in)    :: output_dir !< The directory into which to write files.
  type(thermo_var_ptrs),   intent(inout) :: tv         !< The thermodynamic variable structure.
  real,                    intent(in)    :: max_depth  !< The ocean's maximum depth [Z ~> m].
  ! Local
  character(len=200) :: config
  logical :: debug
! This include declares and sets the variable "version".
#include "version_variable.h"
  integer :: nz
  type(EOS_type), pointer :: eos => NULL()

  if (associated(tv%eqn_of_state)) eos => tv%eqn_of_state

  nz = GV%ke

  call callTree_enter("MOM_initialize_coord(), MOM_coord_initialization.F90")
  call log_version(PF, mdl, version, "")
  call get_param(PF, mdl, "DEBUG", debug, default=.false.)

! Set-up the layer densities, GV%Rlay, and reduced gravities, GV%g_prime.
  call get_param(PF, mdl, "COORD_CONFIG", config, &
                 "This specifies how layers are to be defined: \n"//&
                 " \t ALE or none - used to avoid defining layers in ALE mode \n"//&
                 " \t file - read coordinate information from the file \n"//&
                 " \t\t specified by (COORD_FILE).\n"//&
                 " \t BFB - Custom coords for buoyancy-forced basin case \n"//&
                 " \t\t based on SST_S, T_BOT and DRHO_DT.\n"//&
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
                 default="none")
  select case ( trim(config) )
    case ("gprime")
      call set_coord_from_gprime(GV%Rlay, GV%g_prime, GV, US, PF)
    case ("layer_ref")
      call set_coord_from_layer_density(GV%Rlay, GV%g_prime, GV, US, PF)
    case ("linear")
      call set_coord_linear(GV%Rlay, GV%g_prime, GV, US, PF)
    case ("ts_ref")
      call set_coord_from_TS_ref(GV%Rlay, GV%g_prime, GV, US, PF, eos, tv%P_Ref)
    case ("ts_profile")
      call set_coord_from_TS_profile(GV%Rlay, GV%g_prime, GV, US, PF, eos, tv%P_Ref)
    case ("ts_range")
      call set_coord_from_TS_range(GV%Rlay, GV%g_prime, GV, US, PF, eos, tv%P_Ref)
    case ("file")
      call set_coord_from_file(GV%Rlay, GV%g_prime, GV, US, PF)
    case ("USER")
      call user_set_coord(GV%Rlay, GV%g_prime, GV, US, PF, eos)
    case ("BFB")
      call BFB_set_coord(GV%Rlay, GV%g_prime, GV, US, PF, eos)
    case ("none", "ALE")
      call set_coord_to_none(GV%Rlay, GV%g_prime, GV, US, PF)
    case default ; call MOM_error(FATAL,"MOM_initialize_coord: "// &
      "Unrecognized coordinate setup"//trim(config))
  end select
  ! There are nz+1 values of g_prime because it is an interface field, but the value at the bottom
  ! should not matter.  This is here just to avoid having an uninitialized value in some output.
  GV%g_prime(nz+1) = 10.0*GV%g_Earth

  if (debug) call chksum(US%R_to_kg_m3*GV%Rlay(:), "MOM_initialize_coord: Rlay ", 1, nz)
  if (debug) call chksum(US%m_to_Z*US%L_to_m**2*US%s_to_T**2*GV%g_prime(:), "MOM_initialize_coord: g_prime ", 1, nz)
  call setVerticalGridAxes( GV%Rlay, GV, scale=US%R_to_kg_m3 )

! Copy the maximum depth across from the input argument
  GV%max_depth = max_depth

! Write out all of the grid data used by this run.
  if (write_geom) call write_vertgrid_file(GV, US, PF, output_dir)

  call callTree_leave('MOM_initialize_coord()')

end subroutine MOM_initialize_coord

!  The set_coord routines deal with initializing aspects of the vertical grid.

!> Sets the layer densities (Rlay) and the interface reduced gravities (g).
subroutine set_coord_from_gprime(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV         !< The ocean's vertical grid structure.
  real, dimension(GV%ke),   intent(out) :: Rlay       !< The layers' target coordinate values
                                                      !! (potential density) [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime    !< The reduced gravity across the interfaces
                                                      !! [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US         !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  ! Local variables
  real :: g_int   ! Reduced gravities across the internal interfaces [L2 Z-1 T-2 ~> m s-2].
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  character(len=40)  :: mdl = "set_coord_from_gprime" ! This subroutine's name.
  integer :: k, nz
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "GFS" , g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "GINT", g_int, &
                 "The reduced gravity across internal interfaces.", &
                 units="m s-2", fail_if_missing=.true., scale=US%m_s_to_L_T**2*US%Z_to_m)

  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = g_int ; enddo
  Rlay(1) = GV%Rho0
  do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(GV%Rho0/GV%g_Earth) ; enddo

  call callTree_leave(trim(mdl)//'()')

end subroutine set_coord_from_gprime

!> Sets the layer densities (Rlay) and the interface reduced gravities (g).
subroutine set_coord_from_layer_density(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV         !< The ocean's vertical grid structure.
  real, dimension(GV%ke),   intent(out) :: Rlay       !< The layers' target coordinate values
                                                      !! (potential density) [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime    !< The reduced gravity across the interfaces
                                                      !! [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US         !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters

  ! Local variables
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  real :: Rlay_Ref! The surface layer's target density [R ~> kg m-3].
  real :: RLay_range ! The range of densities [R ~> kg m-3].
  character(len=40)  :: mdl = "set_coord_from_layer_density" ! This subroutine's name.
  integer :: k, nz
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for layer 1.", &
                 units="kg m-3", default=US%R_to_kg_m3*GV%Rho0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "DENSITY_RANGE", Rlay_range, &
                 "The range of reference potential densities in the layers.", &
                 units="kg m-3", default=2.0, scale=US%kg_m3_to_R)

  Rlay(1) = Rlay_Ref
  do k=2,nz
    Rlay(k) = Rlay(k-1) + RLay_range/(real(nz-1))
  enddo
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=2,nz
    g_prime(k) = (GV%g_Earth/(GV%Rho0)) * (Rlay(k) - Rlay(k-1))
  enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine set_coord_from_layer_density

!> Sets the layer densities (Rlay) and the interface reduced gravities (g) from a profile of g'.
subroutine set_coord_from_TS_ref(Rlay, g_prime, GV, US, param_file, eqn_of_state, P_Ref)
  type(verticalGrid_type),  intent(in)  :: GV         !< The ocean's vertical grid structure.
  real, dimension(GV%ke),   intent(out) :: Rlay       !< The layers' target coordinate values
                                                      !! (potential density) [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime    !< The reduced gravity across the interfaces
                                                      !! [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US         !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(EOS_type),           pointer     :: eqn_of_state !< Equation of state structure
  real,                     intent(in)  :: P_Ref      !< The coordinate-density reference pressure
                                                      !! [R L2 T-2 ~> Pa].

  ! Local variables
  real :: T_ref   ! Reference temperature
  real :: S_ref   ! Reference salinity
  real :: g_int   ! Reduced gravities across the internal interfaces [L2 Z-1 T-2 ~> m s-2].
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  character(len=40)  :: mdl = "set_coord_from_TS_ref" ! This subroutine's name.
  integer :: k, nz
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "T_REF", T_Ref, &
                 "The initial temperature of the lightest layer.", units="degC", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "S_REF", S_Ref, &
                 "The initial salinities.", units="PSU", default=35.0)
  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "GINT", g_int, &
                 "The reduced gravity across internal interfaces.", &
                 units="m s-2", fail_if_missing=.true., scale=US%m_s_to_L_T**2*US%Z_to_m)

!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = g_int ; enddo

!    The uppermost layer's density is set here.  Subsequent layers'  !
!  densities are determined from this value and the g values.        !
!        T0 = 28.228 ; S0 = 34.5848 ; Pref = P_Ref
  call calculate_density(T_ref, S_ref, P_ref, Rlay(1), eqn_of_state)

!    These statements set the layer densities.                       !
  do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(GV%Rho0/GV%g_Earth) ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine set_coord_from_TS_ref

!> Sets the layer densities (Rlay) and the interface reduced gravities (g) from a T-S profile.
subroutine set_coord_from_TS_profile(Rlay, g_prime, GV, US, param_file, eqn_of_state, P_Ref)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(EOS_type),           pointer     :: eqn_of_state !< Equation of state structure
  real,                     intent(in)  :: P_Ref   !< The coordinate-density reference pressure
                                                   !! [R L2 T-2 ~> Pa].

  ! Local variables
  real, dimension(GV%ke) :: T0, S0,  Pref
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  integer :: k, nz
  character(len=40)  :: mdl = "set_coord_from_TS_profile" ! This subroutine's name.
  character(len=200) :: filename, coord_file, inputdir ! Strings for file/path
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "COORD_FILE", coord_file, &
                 "The file from which the coordinate temperatures and "//&
                 "salinities are read.", fail_if_missing=.true.)

  call get_param(param_file,  mdl, "INPUTDIR", inputdir, default=".")
  filename = trim(slasher(inputdir))//trim(coord_file)
  call log_param(param_file, mdl, "INPUTDIR/COORD_FILE", filename)

  call MOM_read_data(filename,"PTEMP",T0(:))
  call MOM_read_data(filename,"SALT",S0(:))

  if (.not.file_exists(filename)) call MOM_error(FATAL, &
      " set_coord_from_TS_profile: Unable to open " //trim(filename))
!    These statements set the interface reduced gravities.           !
  g_prime(1) = g_fs
  do k=1,nz ; Pref(k) = P_Ref ; enddo
  call calculate_density(T0, S0, Pref, Rlay, eqn_of_state, (/1,nz/) )
  do k=2,nz; g_prime(k) = (GV%g_Earth/(GV%Rho0)) * (Rlay(k) - Rlay(k-1)) ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine set_coord_from_TS_profile

!> Sets the layer densities (Rlay) and the interface reduced gravities (g) from a linear T-S profile.
subroutine set_coord_from_TS_range(Rlay, g_prime, GV, US, param_file, eqn_of_state, P_Ref)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  type(EOS_type),           pointer     :: eqn_of_state !< Equation of state structure
  real,                     intent(in)  :: P_Ref   !< The coordinate-density reference pressure
                                                   !! [R L2 T-2 ~> Pa].

  ! Local variables
  real, dimension(GV%ke) :: T0, S0,  Pref
  real :: S_Ref, S_Light, S_Dense ! Salinity range parameters [ppt].
  real :: T_Ref, T_Light, T_Dense ! Temperature range parameters [decC].
  real :: res_rat ! The ratio of density space resolution in the denser part
                  ! of the range to that in the lighter part of the range.
                  ! Setting this greater than 1 increases the resolution for
                  ! the denser water.
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  real :: a1, frac_dense, k_frac
  integer :: k, nz, k_light
  character(len=40)  :: mdl = "set_coord_from_TS_range" ! This subroutine's name.
  character(len=200) :: filename, coord_file, inputdir ! Strings for file/path
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "T_REF", T_Ref, &
                 "The default initial temperatures.", units="degC", default=10.0)
  call get_param(param_file, mdl, "TS_RANGE_T_LIGHT", T_Light, &
                 "The initial temperature of the lightest layer when "//&
                 "COORD_CONFIG is set to ts_range.", units="degC", default=T_Ref)
  call get_param(param_file, mdl, "TS_RANGE_T_DENSE", T_Dense, &
                 "The initial temperature of the densest layer when "//&
                 "COORD_CONFIG is set to ts_range.", units="degC", default=T_Ref)

  call get_param(param_file, mdl, "S_REF", S_Ref, &
                 "The default initial salinities.", units="PSU", default=35.0)
  call get_param(param_file, mdl, "TS_RANGE_S_LIGHT", S_Light, &
                 "The initial lightest salinities when COORD_CONFIG "//&
                 "is set to ts_range.", default = S_Ref, units="PSU")
  call get_param(param_file, mdl, "TS_RANGE_S_DENSE", S_Dense, &
                 "The initial densest salinities when COORD_CONFIG "//&
                 "is set to ts_range.", default = S_Ref, units="PSU")

  call get_param(param_file, mdl, "TS_RANGE_RESOLN_RATIO", res_rat, &
                 "The ratio of density space resolution in the densest "//&
                 "part of the range to that in the lightest part of the "//&
                 "range when COORD_CONFIG is set to ts_range. Values "//&
                 "greater than 1 increase the resolution of the denser water.",&
                 default=1.0, units="nondim")

  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)

  if ((GV%nk_rho_varies > 0) .and. (nz < GV%nk_rho_varies+2)) &
    call MOM_error(FATAL, "set_coord_from_TS_range requires that NZ >= NKML+NKBL+2.")

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
  do k=1,nz ; Pref(k) = P_Ref ; enddo
  call calculate_density(T0, S0, Pref, Rlay, eqn_of_state, (/k_light,nz/) )
  ! Extrapolate target densities for the variable density mixed and buffer layers.
  do k=k_light-1,1,-1
    Rlay(k) = 2.0*Rlay(k+1) - Rlay(k+2)
  enddo
  do k=2,nz ; g_prime(k) = (GV%g_Earth/GV%Rho0) * (Rlay(k) - Rlay(k-1)) ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine set_coord_from_TS_range

! Sets the layer densities (Rlay) and the interface reduced gravities (g) from data in file.
subroutine set_coord_from_file(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters

  ! Local variables
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  integer :: k, nz
  character(len=40)  :: mdl = "set_coord_from_file" ! This subroutine's name.
  character(len=40)  :: coord_var
  character(len=200) :: filename,coord_file,inputdir ! Strings for file/path
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)
  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "COORD_FILE", coord_file, &
                 "The file from which the coordinate densities are read.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "COORD_VAR", coord_var, &
                 "The variable in COORD_FILE that is to be used for the "//&
                 "coordinate densities.", default="Layer")
  filename = trim(inputdir)//trim(coord_file)
  call log_param(param_file, mdl, "INPUTDIR/COORD_FILE", filename)
  if (.not.file_exists(filename)) call MOM_error(FATAL, &
      " set_coord_from_file: Unable to open "//trim(filename))

  call read_axis_data(filename, coord_var, Rlay)
  do k=1,nz ; Rlay(k) = US%kg_m3_to_R*Rlay(k) ; enddo
  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = (GV%g_Earth/(GV%Rho0)) * (Rlay(k) - Rlay(k-1)) ; enddo
  do k=1,nz ; if (g_prime(k) <= 0.0) then
    call MOM_error(FATAL, "MOM_initialization set_coord_from_file: "//&
       "Zero or negative g_primes read from variable "//"Layer"//" in file "//&
       trim(filename))
  endif ; enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine set_coord_from_file

!> Sets the layer densities (Rlay) and the interface
!! reduced gravities (g) according to a linear profile starting at a
!! reference surface layer density and spanning a range of densities
!! to the bottom defined by the parameter RLAY_RANGE
!! (defaulting to 2.0 if not defined)
subroutine set_coord_linear(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters

  ! Local variables
  character(len=40)  :: mdl = "set_coord_linear" ! This subroutine
  real :: Rlay_ref, Rlay_range, g_fs
  integer :: k, nz
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "LIGHTEST_DENSITY", Rlay_Ref, &
                 "The reference potential density used for the surface interface.", &
                 units="kg m-3", default=US%R_to_kg_m3*GV%Rho0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "DENSITY_RANGE", Rlay_range, &
                 "The range of reference potential densities across all interfaces.", &
                 units="kg m-3", default=2.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "GFS", g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)

  ! This following sets the target layer densities such that a the
  ! surface interface has density Rlay_ref and the bottom
  ! is Rlay_range larger
  do k=1,nz
    Rlay(k) = Rlay_Ref + RLay_range*((real(k)-0.5)/real(nz))
  enddo
  ! These statements set the interface reduced gravities.
  g_prime(1) = g_fs
  do k=2,nz
    g_prime(k) = (GV%g_Earth/(GV%Rho0)) * (Rlay(k) - Rlay(k-1))
  enddo

  call callTree_leave(trim(mdl)//'()')
end subroutine set_coord_linear

!> Sets Rlay to Rho0 and g_prime to zero except for the free surface.
!! This is for use only in ALE mode where Rlay should not be used and g_prime(1) alone
!! might be used.
subroutine set_coord_to_none(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure to parse for run-time parameters
  ! Local variables
  real :: g_fs    ! Reduced gravity across the free surface [L2 Z-1 T-2 ~> m s-2].
  character(len=40)  :: mdl = "set_coord_to_none" ! This subroutine's name.
  integer :: k, nz
  nz = GV%ke

  call callTree_enter(trim(mdl)//"(), MOM_coord_initialization.F90")

  call get_param(param_file, mdl, "GFS" , g_fs, &
                 "The reduced gravity at the free surface.", units="m s-2", &
                 default=GV%g_Earth*US%L_T_to_m_s**2*US%m_to_Z, scale=US%m_s_to_L_T**2*US%Z_to_m)

  g_prime(1) = g_fs
  do k=2,nz ; g_prime(k) = 0. ; enddo
  Rlay(1) = GV%Rho0
  do k=2,nz ; Rlay(k) = Rlay(k-1) + g_prime(k)*(GV%Rho0/GV%g_Earth) ; enddo

  call callTree_leave(trim(mdl)//'()')

end subroutine set_coord_to_none

!> Writes out a file containing any available data related
!! to the vertical grid used by the MOM ocean model.
subroutine write_vertgrid_file(GV, US, param_file, directory)
  type(verticalGrid_type), intent(in)  :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)  :: US         !< A dimensional unit scaling type
  type(param_file_type), intent(in)    :: param_file !< A structure to parse for run-time parameters
  character(len=*),      intent(in)    :: directory  !< The directory into which to place the file.
  ! Local variables
  character(len=240) :: filepath
  type(vardesc) :: vars(2)
  type(fieldtype) :: fields(2)
  integer :: unit

  filepath = trim(directory) // trim("Vertical_coordinate")

  vars(1) = var_desc("R","kilogram meter-3","Target Potential Density",'1','L','1')
  vars(2) = var_desc("g","meter second-2","Reduced gravity",'1','L','1')

  call create_file(unit, trim(filepath), vars, 2, fields, SINGLE_FILE, GV=GV)

  call write_field(unit, fields(1), US%R_to_kg_m3*GV%Rlay(:))
  call write_field(unit, fields(2), US%L_T_to_m_s**2*US%m_to_Z*GV%g_prime(:))

  call close_file(unit)

end subroutine write_vertgrid_file

end module MOM_coord_initialization
