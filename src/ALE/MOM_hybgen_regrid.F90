!> This module contains the hybgen regridding routines from HYCOM, with minor
!! modifications to follow the MOM6 coding conventions
module MOM_hybgen_regrid

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS,              only : EOS_type, calculate_density
use MOM_error_handler,    only : MOM_mesg, MOM_error, FATAL, WARNING, assert
use MOM_file_parser,      only : get_param, param_file_type, log_param
use MOM_io,               only : close_file, create_file, file_type, fieldtype, file_exists
use MOM_io,               only : MOM_read_data, MOM_write_field, vardesc, var_desc, SINGLE_FILE
use MOM_string_functions, only : slasher
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : ocean_grid_type, thermo_var_ptrs
use MOM_verticalGrid,     only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

!> Control structure containing required parameters for the hybgen coordinate generator
type, public :: hybgen_regrid_CS ; private

  real :: min_thickness !< Minimum thickness allowed for layers [H ~> m or kg m-2]

  integer :: nk         !< Number of layers on the target grid

  !> Reference pressure for density calculations [R L2 T-2 ~> Pa]
  real :: ref_pressure

  !> Hybgen uses PCM if layer is within hybiso of target density [R ~> kg m-3]
  real :: hybiso
  !> Number of sigma levels used by HYBGEN
  integer :: nsigma

  real :: dp00i    !< Deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
  real :: qhybrlx  !< Fractional relaxation within a regridding step [nondim]

  real, allocatable, dimension(:) ::  &
    dp0k, & !< minimum deep    z-layer separation [H ~> m or kg m-2]
    ds0k    !< minimum shallow z-layer separation [H ~> m or kg m-2]

  real :: coord_scale = 1.0     !< A scaling factor to restores the depth coordinates to values in m
  real :: Rho_coord_scale = 1.0 !< A scaling factor to restores the denesity coordinates to values in kg m-3

  real :: dpns  !< depth to start terrain following [H ~> m or kg m-2]
  real :: dsns  !< depth to stop terrain following [H ~> m or kg m-2]
  real :: min_dilate !< The minimum amount of dilation that is permitted when converting target
                     !! coordinates from z to z* [nondim].  This limit applies when wetting occurs.
  real :: max_dilate !< The maximum amount of dilation that is permitted when converting target
                     !! coordinates from z to z* [nondim].  This limit applies when drying occurs.

  real :: thkbot !< Thickness of a bottom boundary layer, within which hybgen does
                 !! something different. [H ~> m or kg m-2]

  !> Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  real :: topiso_const
  ! real, dimension(:,:), allocatable :: topiso

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:) :: target_density

  real :: onem       !< Nominally one m in thickness units [H ~> m or kg m-2]

end type hybgen_regrid_CS


public hybgen_regrid, init_hybgen_regrid, end_hybgen_regrid
public hybgen_column_init, get_hybgen_regrid_params, write_Hybgen_coord_file

contains

!> Initialise a hybgen_regrid_CS control structure and store its parameters
subroutine init_hybgen_regrid(CS, GV, US, param_file)
  type(hybgen_regrid_CS),  pointer    :: CS  !< Unassociated pointer to hold the control structure
  type(verticalGrid_type), intent(in) :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in) :: param_file !< Parameter file

  character(len=40) :: mdl = "MOM_hybgen_regrid" ! This module's name.
  real    :: hybrlx  ! The number of remappings over which to move toward the target coordinate [timesteps]
  character(len=40)  :: dp0_coord_var, ds0_coord_var, rho_coord_var
  character(len=200) :: filename, coord_file, inputdir ! Strings for file/path
  logical :: use_coord_file
  integer :: k

  if (associated(CS)) call MOM_error(FATAL, "init_hybgen_regrid: CS already associated!")
  allocate(CS)

  CS%nk = GV%ke

  allocate(CS%target_density(CS%nk))
  allocate(CS%dp0k(CS%nk), source=0.0) ! minimum deep z-layer separation
  allocate(CS%ds0k(CS%nk), source=0.0) ! minimum shallow z-layer separation

  do k=1,CS%nk ;  CS%target_density(k) = GV%Rlay(k) ; enddo

  call get_param(param_file, mdl, "P_REF", CS%ref_pressure, &
                 "The pressure that is used for calculating the coordinate "//&
                 "density.  (1 Pa = 1e4 dbar, so 2e7 is commonly used.) "//&
                 "This is only used if USE_EOS and ENABLE_THERMODYNAMICS are true.", &
                 units="Pa", default=2.0e7, scale=US%kg_m3_to_R*US%m_s_to_L_T**2)

  call get_param(param_file, mdl, "HYBGEN_MIN_THICKNESS", CS%min_thickness, &
                 "The minimum layer thickness allowed when regridding with Hybgen.",  &
                 units="m", default=0.0, scale=GV%m_to_H )

  call get_param(param_file, mdl, "HYBGEN_N_SIGMA", CS%nsigma, &
                 "The number of sigma-coordinate (terrain-following) layers with Hybgen regridding.", &
                 default=0)
  call get_param(param_file, mdl, "HYBGEN_COORD_FILE", coord_file, &
                 "The file from which the Hybgen profile is read, or blank to use a list of "//&
                 "real input parameters from the MOM_input file.", default="")

  use_coord_file = (len_trim(coord_file) > 0)
  call get_param(param_file, mdl, "HYBGEN_DEEP_DZ_PR0FILE", CS%dp0k, &
                 "The layerwise list of deep z-level minimum thicknesses for Hybgen (dp0k in Hycom).", &
                 units="m", default=0.0, scale=GV%m_to_H, do_not_log=use_coord_file)
  call get_param(param_file, mdl, "HYBGEN_SHALLOW_DZ_PR0FILE", CS%ds0k, &
                 "The layerwise list of shallow z-level minimum thicknesses for Hybgen (ds0k in Hycom).", &
                 units="m", default=0.0, scale=GV%m_to_H, do_not_log=use_coord_file)

  if (use_coord_file) then
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)
    filename = trim(inputdir)//trim(coord_file)
    call log_param(param_file, mdl, "INPUTDIR/HYBGEN_COORD_FILE", filename)
    if (.not.file_exists(filename)) call MOM_error(FATAL, &
        " set_coord_from_file: Unable to open "//trim(filename))

    call get_param(param_file, mdl, "HYBGEN_DEEP_DZ_VAR", dp0_coord_var, &
                 "The variable in HYBGEN_COORD_FILE that is to be used for the "//&
                 "deep z-level minimum thicknesses for Hybgen (dp0k in Hycom).", &
                 default="dp0")
    call get_param(param_file, mdl, "HYBGEN_SHALLOW_DZ_VAR", ds0_coord_var, &
                 "The variable in HYBGEN_COORD_FILE that is to be used for the "//&
                 "shallow z-level minimum thicknesses for Hybgen (ds0k in Hycom).", &
                 default="ds0")
    call get_param(param_file, mdl, "HYBGEN_TGT_DENSITY_VAR", rho_coord_var, &
                 "The variable in HYBGEN_COORD_FILE that is to be used for the Hybgen "//&
                 "target layer densities, or blank to reuse the values in GV%Rlay.", &
                 default="")

    call MOM_read_data(filename, dp0_coord_var, CS%dp0k, scale=GV%m_to_H)

    call MOM_read_data(filename, ds0_coord_var, CS%ds0k, scale=GV%m_to_H)

    if (len_trim(rho_coord_var) > 0) &
      call MOM_read_data(filename, rho_coord_var, CS%target_density, scale=US%kg_m3_to_R)
  endif

  call get_param(param_file, mdl, "HYBGEN_ISOPYCNAL_DZ_MIN", CS%dp00i, &
                 "The Hybgen deep isopycnal spacing minimum thickness (dp00i in Hycom)", &
                 units="m", default=0.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "HYBGEN_MIN_ISO_DEPTH", CS%topiso_const, &
                 "The Hybgen shallowest depth for isopycnal layers (isotop in Hycom)", &
                 units="m", default=0.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "HYBGEN_RELAX_PERIOD", hybrlx, &
                 "The Hybgen coordinate relaxation period in timesteps, or 1 to move to "//&
                 "the new target coordinates in a single step.  This must be >= 1.", &
                 units="timesteps", default=1.0)
  if (hybrlx < 1.0) call MOM_error(FATAL, "init_hybgen_regrid: HYBGEN_RELAX_PERIOD must be at least 1.")
  CS%qhybrlx = 1.0 / hybrlx
  call get_param(param_file, mdl, "HYBGEN_BBL_THICKNESS", CS%thkbot, &
                 "A bottom boundary layer thickness within which Hybgen is able to move "//&
                 "overlying layers upward to match a target density.", &
                 units="m", default=0.0, scale=GV%m_to_H)
  call get_param(param_file, mdl, "HYBGEN_REMAP_DENSITY_MATCH", CS%hybiso, &
                 "A tolerance between the layer densities and their target, within which "//&
                 "Hybgen determines that remapping uses PCM for a layer.", &
                 units="kg m-3", default=0.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "HYBGEN_REMAP_MIN_ZSTAR_DILATE", CS%min_dilate, &
                 "The maximum amount of dilation that is permitted when converting target "//&
                 "coordinates from z to z* [nondim].  This limit applies when drying occurs.", &
                 default=0.5)
  call get_param(param_file, mdl, "HYBGEN_REMAP_MAX_ZSTAR_DILATE", CS%max_dilate, &
                 "The maximum amount of dilation that is permitted when converting target "//&
                 "coordinates from z to z* [nondim].  This limit applies when drying occurs.", &
                 default=2.0)

  CS%onem = 1.0 * GV%m_to_H

  do k=1,CS%nk ; CS%dp0k(k) = max(CS%dp0k(k), CS%min_thickness) ; enddo
  CS%dp00i = max(CS%dp00i, CS%min_thickness)

  ! Determine the depth range over which to use a sigma (terrain-following) coordinate.
  ! --- terrain following starts at depth dpns and ends at depth dsns
  if (CS%nsigma == 0) then
    CS%dpns = CS%dp0k(1)
    CS%dsns = 0.0
  else
    CS%dpns = 0.0
    CS%dsns = 0.0
    do k=1,CS%nsigma
      CS%dpns = CS%dpns + CS%dp0k(k)
      CS%dsns = CS%dsns + CS%ds0k(k)
    enddo !k
  endif !nsigma

  CS%coord_scale = GV%H_to_m
  CS%Rho_coord_scale = US%R_to_kg_m3

end subroutine init_hybgen_regrid

!> Writes out a file containing any available data related
!! to the vertical grid used by the MOM ocean model.
subroutine write_Hybgen_coord_file(GV, CS, filepath)
  type(verticalGrid_type), intent(in)  :: GV        !< The ocean's vertical grid structure
  type(hybgen_regrid_CS),  intent(in)  :: CS        !< Control structure for this module
  character(len=*),        intent(in)  :: filepath  !< The full path to the file to write
  ! Local variables
  type(vardesc) :: vars(3)
  type(fieldtype) :: fields(3)
  type(file_type) :: IO_handle ! The I/O handle of the fileset

  vars(1) = var_desc("dp0", "meter", "Deep z-level minimum thicknesses for Hybgen", '1', 'L', '1')
  vars(2) = var_desc("ds0", "meter", "Shallow z-level minimum thicknesses for Hybgen", '1', 'L', '1')
  vars(3) = var_desc("Rho_tgt", "kg m-3", "Target coordinate potential densities for Hybgen", '1', 'L', '1')
  call create_file(IO_handle, trim(filepath), vars, 3, fields, SINGLE_FILE, GV=GV)

  call MOM_write_field(IO_handle, fields(1), CS%dp0k, scale=CS%coord_scale)
  call MOM_write_field(IO_handle, fields(2), CS%ds0k, scale=CS%coord_scale)
  call MOM_write_field(IO_handle, fields(3), CS%target_density, scale=CS%Rho_coord_scale)

  call close_file(IO_handle)

end subroutine write_Hybgen_coord_file

!> This subroutine deallocates memory in the control structure for the hybgen module
subroutine end_hybgen_regrid(CS)
  type(hybgen_regrid_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return

  deallocate(CS%target_density)
  deallocate(CS%dp0k, CS%ds0k)
  deallocate(CS)
end subroutine end_hybgen_regrid

!> This subroutine can be used to retrieve the parameters for the hybgen regrid module
subroutine get_hybgen_regrid_params(CS, nk, ref_pressure, hybiso, nsigma, dp00i, qhybrlx, &
                                    dp0k, ds0k, dpns, dsns, min_dilate, max_dilate, &
                                    thkbot, topiso_const, target_density)
  type(hybgen_regrid_CS),  pointer    :: CS !< Coordinate regridding control structure
  integer, optional, intent(out) :: nk  !< Number of layers on the target grid
  real,    optional, intent(out) :: ref_pressure !< Reference pressure for density calculations [R L2 T-2 ~> Pa]
  real,    optional, intent(out) :: hybiso  !< Hybgen uses PCM if layer is within hybiso of target density [R ~> kg m-3]
  integer, optional, intent(out) :: nsigma  !< Number of sigma levels used by HYBGEN
  real,    optional, intent(out) :: dp00i   !< Deep isopycnal spacing minimum thickness (m)
  real,    optional, intent(out) :: qhybrlx !< Fractional relaxation amount per timestep, 0 < qyhbrlx <= 1 [nondim]
  real,    optional, intent(out) :: dp0k(:) !< minimum deep    z-layer separation [H ~> m or kg m-2]
  real,    optional, intent(out) :: ds0k(:) !< minimum shallow z-layer separation [H ~> m or kg m-2]
  real,    optional, intent(out) :: dpns    !< depth to start terrain following [H ~> m or kg m-2]
  real,    optional, intent(out) :: dsns    !< depth to stop terrain following [H ~> m or kg m-2]
  real,    optional, intent(out) :: min_dilate !< The minimum amount of dilation that is permitted when
                                            !! converting target coordinates from z to z* [nondim].
                                            !! This limit applies when wetting occurs.
  real,    optional, intent(out) :: max_dilate !< The maximum amount of dilation that is permitted when
                                            !! converting target coordinates from z to z* [nondim].
                                            !! This limit applies when drying occurs.
  real,    optional, intent(out) :: thkbot  !< Thickness of a bottom boundary layer, within which
                                            !! hybgen does something different. [H ~> m or kg m-2]
  real,    optional, intent(out) :: topiso_const !< Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  ! real, dimension(:,:), allocatable :: topiso
  real,    optional, intent(out) :: target_density(:) !< Nominal density of interfaces [R ~> kg m-3]

  if (.not. associated(CS)) call MOM_error(FATAL, "get_hybgen_params: CS not associated")

  if (present(nk))      nk = CS%nk
  if (present(ref_pressure)) ref_pressure = CS%ref_pressure
  if (present(hybiso))  hybiso = CS%hybiso
  if (present(nsigma))  nsigma = CS%nsigma
  if (present(dp00i))   dp00i = CS%dp00i
  if (present(qhybrlx)) qhybrlx = CS%qhybrlx
  if (present(dp0k)) then
    if (size(dp0k) < CS%nk) call MOM_error(FATAL, "get_hybgen_regrid_params: "//&
                                    "The dp0k argument is not allocated with enough space.")
    dp0k(1:CS%nk) = CS%dp0k(1:CS%nk)
  endif
  if (present(ds0k)) then
    if (size(ds0k) < CS%nk) call MOM_error(FATAL, "get_hybgen_regrid_params: "//&
                                    "The ds0k argument is not allocated with enough space.")
    ds0k(1:CS%nk) = CS%ds0k(1:CS%nk)
  endif
  if (present(dpns))    dpns = CS%dpns
  if (present(dsns))    dsns = CS%dsns
  if (present(min_dilate)) min_dilate = CS%min_dilate
  if (present(max_dilate)) max_dilate = CS%max_dilate
  if (present(thkbot))  thkbot = CS%thkbot
  if (present(topiso_const)) topiso_const = CS%topiso_const
  if (present(target_density)) then
    if (size(target_density) < CS%nk) call MOM_error(FATAL, "get_hybgen_regrid_params: "//&
                                    "The target_density argument is not allocated with enough space.")
    target_density(1:CS%nk) = CS%target_density(1:CS%nk)
  endif

end subroutine get_hybgen_regrid_params


!> Modify the input grid to give a new vertical grid based on the HYCOM hybgen code.
subroutine hybgen_regrid(G, GV, US, dp, tv, CS, dzInterface, PCM_cell)
  type(ocean_grid_type),   intent(in)    :: G   !< Ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV  !< Ocean vertical grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: dp  !< Source grid layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),   intent(in)    :: tv  !< Thermodynamics structure
  type(hybgen_regrid_CS),  intent(in)    :: CS  !< hybgen control structure
  real, dimension(SZI_(G),SZJ_(G),CS%nk+1), &
                           intent(inout) :: dzInterface !< The change in height of each interface,
                                                !! using a sign convention opposite to the change
                                                !! in pressure [H ~> m or kg m-2]
  logical, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(inout) :: PCM_cell !< If true, PCM remapping should be used in a cell.
                                                !! This is effectively intent out, but values in wide
                                                !! halo regions and land points are reused.

  ! --- -------------------------------------
  ! --- hybrid grid generator from HYCOM
  ! --- -------------------------------------

  ! These notes on the parameters for the hybrid grid generator are inhereted from the
  ! Hycom source code for these algorithms.
  !
  ! From blkdat.input (units may have changed from m to pressure):
  !
  ! --- 'nsigma' = number of sigma  levels
  ! --- 'dp0k  ' = layer k deep    z-level spacing minimum thickness (m)
  ! ---              k=1,nk
  ! --- 'ds0k  ' = layer k shallow z-level spacing minimum thickness (m)
  ! ---              k=1,nsigma
  ! --- 'dp00i'  = deep isopycnal spacing minimum thickness (m)
  ! --- 'isotop' = shallowest depth for isopycnal layers     (m)
  !                now in topiso(:,:)
  ! --- 'sigma ' = isopycnal layer target densities (sigma units)
  ! ---            now in theta(:,:,1:nk)
  !
  ! --- the above specifies a vertical coord. that is isopycnal or:
  ! ---  near surface z in    deep water, based on dp0k
  ! ---  near surface z in shallow water, based on ds0k and nsigma
  ! ---   terrain-following between them, based on ds0k and nsigma
  !
  ! --- terrain following starts at depth dpns=sum(dp0k(k),k=1,nsigma) and
  ! --- ends at depth dsns=sum(ds0k(k),k=1,nsigma), and the depth of the
  ! --- k-th layer interface varies linearly with total depth between
  ! --- these two reference depths, i.e. a z-sigma-z fixed coordinate.
  !
  ! --- near the surface (i.e. shallower than isotop), layers are always
  ! --- fixed depth (z or sigma).
  ! --  layer 1 is always fixed, so isotop=0.0 is not realizable.
  ! --- near surface layers can also be forced to be fixed depth
  ! --- by setting target densities (sigma(k)) very small.
  !
  ! --- away from the surface, the minimum layer thickness is dp00i.
  !
  ! --- for fixed depth targets to be:
  ! ---  z-only set nsigma=0,
  ! ---  sigma-z (shallow-deep) use a very small ds0k(:),
  ! ---  sigma-only set nsigma=nk, dp0k large, and ds0k small.

  ! These arrays work with the input column
  real :: p_col(GV%ke)      ! A column of reference pressures [R L2 T-2 ~> Pa]
  real :: temp_in(GV%ke)    ! A column of input potential temperatures [degC]
  real :: saln_in(GV%ke)    ! A column of input layer salinities [ppt]
  real :: Rcv_in(GV%ke)     ! An input column of coordinate potential density [R ~> kg m-3]
  real :: dp_in(GV%ke)      ! The input column of layer thicknesses [H ~> m or kg m-2]
  logical :: PCM_lay(GV%ke) ! If true for a layer, use PCM remapping for that layer

  ! These arrays are on the target grid.
  real :: Rcv_tgt(CS%nk)    ! Target potential density [R ~> kg m-3]
  real :: Rcv(CS%nk)        ! Initial values of coordinate potential density on the target grid [R ~> kg m-3]
  real :: h_col(CS%nk)      ! A column of layer thicknesses [H ~> m or kg m-2]
  real :: dz_int(CS%nk+1)   ! The change in interface height due to remapping [H ~> m or kg m-2]
  real :: Rcv_integral      ! Integrated coordinate potential density in a layer [R H ~> kg m-2 or kg2 m-5]

  real :: qhrlx(CS%nk+1)    ! Fractional relaxation within a timestep (between 0 and 1) [nondim]
  real :: dp0ij(CS%nk)      ! minimum layer thickness [H ~> m or kg m-2]
  real :: dp0cum(CS%nk+1)   ! minimum interface depth [H ~> m or kg m-2]

  real :: h_tot             ! Total thickness of the water column [H ~> m or kg m-2]
  real :: nominalDepth      ! Depth of ocean bottom (positive downward) [H ~> m or kg m-2]
  real :: dilate            ! A factor by which to dilate the target positions from z to z* [nondim]
  integer :: fixlay         ! Deepest fixed coordinate layer
  integer, dimension(0:CS%nk) :: k_end ! The index of the deepest source layer that contributes to
                            ! each target layer, in the unusual case where the the input grid is
                            ! larger than the new grid.  This situation only occurs during certain
                            ! types of initialization or when generating output diagnostics.
  integer :: i, j, k, nk, m, k2, nk_in

  nk = CS%nk

  p_col(:) = CS%ref_pressure

  do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1 ; if (G%mask2dT(i,j)>0.) then

    ! Store one-dimensional arrays of thicknesses for the 'old'  vertical grid before regridding
    h_tot = 0.0
    do K=1,GV%ke
      temp_in(k) = tv%T(i,j,k)
      saln_in(k) = tv%S(i,j,k)
      dp_in(k) = dp(i,j,k)
      h_tot = h_tot + dp_in(k)
    enddo

    ! This sets the input column's coordinate potential density from T and S.
    call calculate_density(temp_in, saln_in, p_col, Rcv_in, tv%eqn_of_state)

    ! Set the initial properties on the new grid from the old grid.
    nk_in = GV%ke
    if (GV%ke > CS%nk) then ; do k=GV%ke,CS%nk+1,-1
      ! Remove any excess massless layers from the bottom of the input column.
      if (dp_in(k) > 0.0) exit
      nk_in = k-1
    enddo ; endif

    if (CS%nk >= nk_in) then
      ! Simply copy over the common layers.  This is the usual case.
      do k=1,min(CS%nk,GV%ke)
        h_col(k) = dp_in(k)
        Rcv(k) = Rcv_in(k)
      enddo
      if (CS%nk > GV%ke) then
        ! Pad out the input column with additional massless layers with the bottom properties.
        ! This case only occurs during initialization or perhaps when writing diagnostics.
        do k=GV%ke+1,CS%nk
          Rcv(k) = Rcv_in(GV%ke)
          h_col(k) = 0.0
        enddo
      endif
    else ! (CS%nk < nk_in)
      ! The input column has more data than the output.  For now, combine layers to
      ! make them the same size, but there may be better approaches that should be taken.
      ! This case only occurs during initialization or perhaps when writing diagnostics.
      ! This case was not handled by the original Hycom code in hybgen.F90.
      do k=0,CS%nk ; k_end(k) = (k * nk_in) / CS%nk ; enddo
      do k=1,CS%nk
        h_col(k) = 0.0 ; Rcv_integral = 0.0
        do k2=k_end(k-1) + 1,k_end(k)
          h_col(k) = h_col(k) + dp_in(k2)
          Rcv_integral = Rcv_integral + dp_in(k2)*Rcv_in(k2)
        enddo
        if (h_col(k) > GV%H_subroundoff) then
          ! Take the volume-weighted average properties.
          Rcv(k) = Rcv_integral / h_col(k)
        else ! Take the properties of the topmost source layer that contributes.
          Rcv(k) = Rcv_in(k_end(k-1) + 1)
        endif
      enddo
    endif

    ! Set the target densities for the new layers.
    do k=1,CS%nk
      ! Rcv_tgt(k) = theta(i,j,k)  ! If a 3-d target density were set up in theta, use that here.
      Rcv_tgt(k) = CS%target_density(k)  ! MOM6 does not yet support 3-d target densities.
    enddo

    ! The following block of code is used to trigger z* stretching of the targets heights.
    nominalDepth = (G%bathyT(i,j) + G%Z_ref)*GV%Z_to_H
    if (h_tot <= CS%min_dilate*nominalDepth) then
      dilate = CS%min_dilate
    elseif (h_tot >= CS%max_dilate*nominalDepth) then
      dilate = CS%max_dilate
    else
      dilate = h_tot / nominalDepth
    endif

    ! Convert the regridding parameters into specific constraints for this column.
    call hybgen_column_init(nk, CS%nsigma, CS%dp0k, CS%ds0k, CS%dp00i, &
                            CS%topiso_const, CS%qhybrlx, CS%dpns, CS%dsns, h_tot, dilate, &
                            h_col, fixlay, qhrlx, dp0ij, dp0cum)

    ! Determine whether to require the use of PCM remapping from each source layer.
    do k=1,GV%ke
      if (CS%hybiso > 0.0) then
        ! --- thin or isopycnal source layers are remapped with PCM.
        PCM_lay(k) = (k > fixlay) .and. (abs(Rcv(k) - Rcv_tgt(k)) < CS%hybiso)
      else ! hybiso==0.0, so purely isopycnal layers use PCM
        PCM_lay(k) = .false.
      endif ! hybiso
    enddo !k

    ! Determine the new layer thicknesses.
    call hybgen_column_regrid(CS, nk, CS%thkbot, CS%onem, &
                              1.0e-11*US%kg_m3_to_R, Rcv_tgt, fixlay, qhrlx, dp0ij, &
                              dp0cum, Rcv, h_col, dz_int)

    ! Store the output from hybgenaij_regrid in 3-d arrays.
    if (present(PCM_cell)) then ; do k=1,GV%ke
      PCM_cell(i,j,k) = PCM_lay(k)
    enddo ; endif

    do K=1,nk+1
      ! Note that dzInterface uses the opposite sign convention from the change in p.
      dzInterface(i,j,K) = -dz_int(K)
    enddo

  else
    if (present(PCM_cell)) then ; do k=1,GV%ke
      PCM_cell(i,j,k) = .false.
    enddo ; endif
    do k=1,CS%nk+1 ; dzInterface(i,j,k) = 0.0 ; enddo
  endif ; enddo ; enddo !i & j.

end subroutine hybgen_regrid

!> Initialize some of the variables that are used for regridding or unmixing, including the
!! stretched contraits on where the new interfaces can be.
subroutine hybgen_column_init(nk, nsigma, dp0k, ds0k, dp00i, topiso_i_j, &
                          qhybrlx, dpns, dsns, h_tot, dilate, h_col, &
                          fixlay, qhrlx, dp0ij, dp0cum)
  integer, intent(in)    :: nk           !< The number of layers in the new grid
  integer, intent(in)    :: nsigma       !< The number of sigma  levels
  real,    intent(in)    :: dp0k(nk)     !< Layer deep z-level spacing minimum thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: ds0k(nsigma) !< Layer shallow z-level spacing minimum thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: dp00i        !< Deep isopycnal spacing minimum thickness [H ~> m or kg m-2]
  real,    intent(in)    :: topiso_i_j   !< Shallowest depth for isopycnal layers [H ~> m or kg m-2]
  real,    intent(in)    :: qhybrlx      !< Fractional relaxation amount per timestep, 0 < qyhbrlx <= 1 [nondim]
  real,    intent(in)    :: h_tot        !< The sum of the initial layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: dilate       !< A factor by which to dilate the target positions
                                         !! from z to z* [nondim]
  real,    intent(in)    :: h_col(nk)    !< Initial layer thicknesses [H ~> m or kg m-2]
  real,    intent(in)    :: dpns         !< Vertical sum of dp0k [H ~> m or kg m-2]
  real,    intent(in)    :: dsns         !< Vertical sum of ds0k [H ~> m or kg m-2]
  integer, intent(out)   :: fixlay       !< Deepest fixed coordinate layer
  real,    intent(out)   :: qhrlx(nk+1)  !< Fractional relaxation within a timestep (between 0 and 1) [nondim]
  real,    intent(out)   :: dp0ij(nk)    !< minimum layer thickness [H ~> m or kg m-2]
  real,    intent(out)   :: dp0cum(nk+1) !< minimum interface depth [H ~> m or kg m-2]

  ! --- --------------------------------------------------------------
  ! --- hybrid grid generator, single column - initialization.
  ! --- --------------------------------------------------------------

  ! Local variables
  character(len=256) :: mesg  ! A string for output messages
  real :: hybrlx  ! The relaxation rate in the hybrid region [timestep-1]?
  real :: qdep    ! Total water column thickness as a fraction of dp0k (vs ds0k) [nondim]
  real :: q       ! A portion of the thickness that contributes to the new cell [H ~> m or kg m-2]
  real :: p_int(nk+1)  ! Interface depths [H ~> m or kg m-2]
  integer :: k, fixall

  ! --- dpns = sum(dp0k(k),k=1,nsigma)
  ! --- dsns = sum(ds0k(k),k=1,nsigma)
  ! --- terrain following starts (on the deep side) at depth dpns and ends (on the
  ! --- shallow side) at depth dsns and the depth of the k-th layer interface varies
  ! --- linearly with total depth between these two reference depths.
  if ((h_tot >= dilate * dpns) .or. (dpns <= dsns)) then
    qdep = 1.0  ! Not terrain following - this column is too thick or terrain following is disabled.
  elseif (h_tot <= dilate * dsns) then
    qdep = 0.0  ! Not terrain following - this column is too thin
  else
    qdep = (h_tot - dilate * dsns) / (dilate * (dpns - dsns))
  endif

  if (qdep < 1.0) then
    ! Terrain following or shallow fixed coordinates, qhrlx=1 and ignore dp00
    p_int( 1) = 0.0
    dp0cum(1) = 0.0
    qhrlx( 1) = 1.0
    dp0ij( 1) = dilate * (qdep*dp0k(1) + (1.0-qdep)*ds0k(1))

    dp0cum(2) = dp0cum(1) + dp0ij(1)
    qhrlx( 2) = 1.0
    p_int( 2) = p_int(1) + h_col(1)
    do k=2,nk
      qhrlx( k+1) = 1.0
      dp0ij( k)   = dilate * (qdep*dp0k(k) + (1.0-qdep)*ds0k(k))
      dp0cum(k+1) = dp0cum(k) + dp0ij(k)
      p_int( k+1) = p_int(k) + h_col(k)
    enddo !k
  else
    ! Not terrain following
    p_int( 1) = 0.0
    dp0cum(1) = 0.0
    qhrlx( 1) = 1.0 !no relaxation in top layer
    dp0ij( 1) = dilate * dp0k(1)

    dp0cum(2) = dp0cum(1) + dp0ij(1)
    qhrlx( 2) = 1.0 !no relaxation in top layer
    p_int( 2) = p_int(1) + h_col(1)
    do k=2,nk
      if ((dp0k(k) <= dp00i) .or. (dilate * dp0k(k) >= p_int(k) - dp0cum(k))) then
        ! This layer is in fixed surface coordinates.
        dp0ij(k) = dp0k(k)
        qhrlx(k+1) = 1.0
      else
        q = dp0k(k) * (dilate * dp0k(k) / ( p_int(k) - dp0cum(k)) ) ! A fraction between 0 and 1 of dp0 to use here.
        if (dp00i >= q) then
          ! This layer is much deeper than the fixed surface coordinates.
          dp0ij(k) = dp00i
          qhrlx(k+1) = qhybrlx
        else
          ! This layer spans the margins of the fixed surface coordinates.
          ! In this case dp00i < q < dp0k.
          dp0ij(k) = dilate * q
          qhrlx(k+1) = qhybrlx * (dp0k(k) - dp00i) / &
                       ((dp0k(k) - q) + (q - dp00i)*qhybrlx) ! 1 at dp0k, qhybrlx at dp00i
        endif

        ! The old equivalent code is:
        ! hybrlx = 1.0 / qhybrlx
        ! q = max( dp00i, dp0k(k) * (dp0k(k) / max(dp0k( k), p_int(k) - dp0cum(k)) ) )
        ! qts = 1.0 - (q-dp00i) / (dp0k(k) - dp00i)  !0 at q = dp0k, 1 at q=dp00i
        ! qhrlx( k+1) = 1.0 / (1.0 + qts*(hybrlx-1.0))  !1 at dp0k, qhybrlx at dp00i
      endif
      dp0cum(k+1) = dp0cum(k) + dp0ij(k)
      p_int(k+1) = p_int(k) + h_col(k)
    enddo !k
  endif !qdep<1:else

  ! Identify the current fixed coordinate layers
  fixlay = 1  !layer 1 always fixed
  do k=2,nk
    if (dp0cum(k) >= dilate * topiso_i_j) then
      exit  !layers k to nk might be isopycnal
    endif
    ! Top of layer is above topiso, i.e. always fixed coordinate layer
    qhrlx(k+1) = 1.0  !no relaxation in fixed layers
    fixlay     = fixlay+1
  enddo !k

  fixall = fixlay
  do k=fixall+1,nk
    if (p_int(k+1) > dp0cum(k+1) + 0.1*dp0ij(k)) then
      if ( (fixlay > fixall) .and. (p_int(k) > dp0cum(k)) ) then
        ! --- The previous layer should remain fixed.
        fixlay = fixlay-1
      endif
      exit  !layers k to nk might be isopycnal
    endif
    ! Sometimes fixed coordinate layer
    qhrlx(k) = 1.0  !no relaxation in fixed layers
    fixlay   = fixlay+1
  enddo !k

end subroutine hybgen_column_init

!> The cushion function from Bleck & Benjamin, 1992, which returns a smoothly varying
!! but limited value that goes between dp0 and delp
real function cushn(delp, dp0)
  real, intent(in) :: delp  ! A thickness change [H ~> m or kg m-2]
  real, intent(in) :: dp0   ! A non-negative reference thickness [H ~> m or kg m-2]

  real :: qq  ! A limited ratio of delp/dp0 [nondim]

  ! These are the nondimensional parameters that define the cushion function.
  real, parameter :: qqmn=-4.0, qqmx=2.0  ! shifted range for cushn [nondim]
! real, parameter :: qqmn=-2.0, qqmx=4.0  ! traditional range for cushn [nondim]
! real, parameter :: qqmn=-4.0, qqmx=6.0  ! somewhat wider range for cushn [nondim]
  ! These are derivative nondimensional parameters.
  ! real, parameter :: cusha = qqmn**2 * (qqmx-1.0) / (qqmx-qqmn)**2
  ! real, parameter :: I_qqmn = 1.0 / qqmn
  real, parameter :: qq_scale = (qqmx-1.0) / (qqmx-qqmn)**2
  real, parameter :: I_qqmx = 1.0 / qqmx

  ! --- if delp >= qqmx*dp0 >>  dp0, cushn returns delp.
  ! --- if delp <= qqmn*dp0 << -dp0, cushn returns dp0.

  ! This is the original version from Hycom.
  ! qq = max(qqmn, min(qqmx, delp/dp0))
  ! cushn = dp0 * (1.0 + cusha * (1.0-I_qqmn*qq)**2) * max(1.0, delp/(dp0*qqmx))

  ! This is mathematically equivalent, has one fewer divide, and works as intended even if dp0 = 0.
  if (delp >= qqmx*dp0) then
    cushn = delp
  elseif (delp < qqmn*dp0) then
    cushn = max(dp0, delp * I_qqmx)
  else
    cushn = max(dp0, delp * I_qqmx) * (1.0 + qq_scale * ((delp / dp0) - qqmn)**2)
  endif

end function cushn

!> Create a new grid for a column of water using the Hybgen algorithm.
subroutine hybgen_column_regrid(CS, nk, thkbot, onem, epsil, Rcv_tgt, &
                                fixlay, qhrlx, dp0ij, dp0cum, Rcv, h_in, dp_int)
  type(hybgen_regrid_CS), intent(in)    :: CS  !< hybgen regridding control structure
  integer, intent(in)    :: nk            !< number of layers
  real,    intent(in)    :: thkbot        !< thickness of bottom boundary layer [H ~> m or kg m-2]
  real,    intent(in)    :: onem          !< one m in pressure units [H ~> m or kg m-2]
  real,    intent(in)    :: epsil         !< small nonzero density to prevent division by zero [R ~> kg m-3]
  real,    intent(in)    :: Rcv_tgt(nk)   !< Target potential density [R ~> kg m-3]
  integer, intent(in)    :: fixlay        !< deepest fixed coordinate layer
  real,    intent(in)    :: qhrlx( nk+1)  !< relaxation coefficient per timestep [nondim]
  real,    intent(in)    :: dp0ij( nk)    !< minimum layer thickness [H ~> m or kg m-2]
  real,    intent(in)    :: dp0cum(nk+1)  !< minimum interface depth [H ~> m or kg m-2]
  real,    intent(in)    :: Rcv(nk)       !< Coordinate potential density [R ~> kg m-3]
  real,    intent(in)    :: h_in(nk)      !< Layer thicknesses [H ~> m or kg m-2]
  real,    intent(out)   :: dp_int(nk+1)  !< The change in interface positions [H ~> m or kg m-2]

  ! --- ------------------------------------------------------
  ! --- hybrid grid generator, single column - regrid.
  ! --- ------------------------------------------------------

  ! Local variables
  real :: p_new  ! A new interface position [H ~> m or kg m-2]
  real :: pres_in(nk+1) ! layer interface positions [H ~> m or kg m-2]
  real :: p_int(nk+1)   ! layer interface positions [H ~> m or kg m-2]
  real :: h_col(nk)     ! Updated layer thicknesses [H ~> m or kg m-2]
  real :: q_frac ! A fraction of a layer to entrain [nondim]
  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2]
  real :: h_hat3 ! Thickness movement upward across the interface between layers k-2 and k-3 [H ~> m or kg m-2]
  real :: h_hat2 ! Thickness movement upward across the interface between layers k-1 and k-2 [H ~> m or kg m-2]
  real :: h_hat  ! Thickness movement upward across the interface between layers k and k-1 [H ~> m or kg m-2]
  real :: h_hat0 ! A first guess at thickness movement upward across the interface
                 ! between layers k and k-1 [H ~> m or kg m-2]
  real :: dh_cor ! Thickness changes [H ~> m or kg m-2]
  real :: h1_tgt ! A target thickness for the top layer [H ~> m or kg m-2]
  real :: tenm   ! ten m  in pressure units [H ~> m or kg m-2]
  real :: onemm  ! one mm in pressure units [H ~> m or kg m-2]
  logical :: trap_errors
  integer :: k
  character(len=256) :: mesg  ! A string for output messages

  ! This line needs to be consistent with the parameters set in cushn().
  real, parameter :: qqmn=-4.0, qqmx=2.0  ! shifted range for cushn
! real, parameter :: qqmn=-2.0, qqmx=4.0  ! traditional range for cushn
! real, parameter :: qqmn=-4.0, qqmx=6.0  ! somewhat wider range for cushn

  !### These hard-coded parameters should be changed to run-time variables.
  tenm = 10.0*onem
  onemm = 0.001*onem

  trap_errors = .true.

  do K=1,nk+1 ; dp_int(K) = 0.0 ; enddo

  p_int(1) = 0.0
  do k=1,nk
    h_col(k) = max(h_in(k), 0.0)
    p_int(K+1) = p_int(K) + h_col(k)
  enddo
  h_min = min( CS%min_thickness, p_int(nk+1)/real(CS%nk) )

  if (trap_errors) then
    do K=1,nk+1 ; pres_in(K) = p_int(K) ; enddo
  endif

  ! Try to restore isopycnic conditions by moving layer interfaces
  ! qhrlx(k) are relaxation amounts per timestep.

  ! Maintain prescribed thickness in layer k <= fixlay
  ! There may be massless layers at the bottom, so work upwards.
  do k=min(nk-1,fixlay),1,-1
    p_new = min(dp0cum(k+1), p_int(nk+1) - (nk-k)*h_min) ! This could be positive or negative.
    dh_cor = p_new - p_int(K+1)
    if (k<fixlay) dh_cor = min(dh_cor, h_col(k+1) - h_min)
    h_col(k)    = h_col(k)    + dh_cor
    h_col(k+1)  = h_col(k+1)  - dh_cor
    dp_int(K+1) = dp_int(K+1) + dh_cor
    p_int(K+1) = p_new
  enddo

  ! Eliminate negative thicknesses below the fixed layers, entraining from below as necessary.
  do k=fixlay+1,nk-1
    if (h_col(k) >= h_min) exit  ! usually get here quickly
    dh_cor = h_min - h_col(k)  ! This is positive.
    h_col(k)    = h_min ! = h_col(k) + dh_cor
    h_col(k+1)  = h_col(k+1)  - dh_cor
    dp_int(k+1) = dp_int(k+1) + dh_cor
    p_int(k+1) = p_int(fixlay+1)
  enddo
  if (h_col(nk) < h_min) then  ! This should be uncommon, and should only arise at the level of roundoff.
    do k=nk,2,-1
      if (h_col(k) >= h_min) exit
      dh_cor = h_col(k) - h_min ! dh_cor is negative.
      h_col(k-1) = h_col(k-1) + dh_cor
      h_col(k)   = h_min ! = h_col(k)  - dh_cor
      dp_int(k)  = dp_int(k) + dh_cor
      p_int(k)   = p_int(k) + dh_cor
    enddo
  endif

  ! Remap the non-fixed layers.

  ! In the Hycom version, this loop was fused the loop correcting water that is
  ! too light, and it ran down the water column, but if there are a set of layers
  ! that are very dense, that structure can lead to all of the water being remapped
  ! into a single thick layer.  Splitting the loops and running the loop upwards
  ! (as is done here avoids that catastrophic problem for layers that are far from
  ! their targets.  However, this code is still prone to a thin-thick-thin null mode.
  do k=nk,fixlay+2,-1
    !  This is how the Hycom code would do this loop: do k=fixlay+1,nk ; if (k>fixlay+1) then

    if ((Rcv(k) > Rcv_tgt(k) + epsil)) then
      ! Water in layer k is too dense, so try to dilute with water from layer k-1
      ! Do not move interface if k = fixlay + 1

      if ((Rcv(k-1) >= Rcv_tgt(k-1)) .or. &
          (p_int(k) <= dp0cum(k) + onem) .or. &
          (h_col(k) <= h_col(k-1))) then
        ! If layer k-1 is too light, there is a conflict in the direction the
        ! inteface between them should move, so thicken the thinner of the two.

        if ((Rcv_tgt(k) - Rcv(k-1)) <= epsil) then
          ! layer k-1 is far too dense, take the entire layer
          ! If this code is working downward and this branch is repeated in a series
          ! of successive layers, it can accumulate into a very thick homogenous layers.
          h_hat0 = 0.0  ! This line was not in the Hycom version of hybgen.F90.
          h_hat = dp0ij(k-1) - h_col(k-1)
        else
          ! Entrain enough from the layer above to bring layer k to its target density.
          q_frac = (Rcv_tgt(k) - Rcv(k)) / (Rcv_tgt(k) - Rcv(k-1))    ! -1 <= q_frac < 0
          h_hat0 = q_frac*h_col(k)  ! -h_col(k-1) <= h_hat0 < 0
          if (k == fixlay+2) then
            ! Treat layer k-1 as fixed.
            h_hat = max(h_hat0, dp0ij(k-1) - h_col(k-1))
          else
            ! Maintain the minimum thickess of layer k-1.
            h_hat = cushn(h_hat0 + h_col(k-1), dp0ij(k-1)) - h_col(k-1)
          endif !fixlay+2:else
        endif
        ! h_hat is usually negative, so this check may be unnecessary if the values of
        ! dp0ij are limited to not be below the seafloor?
        h_hat = min(h_hat, p_int(nk+1) - p_int(k))

        ! If isopycnic conditions cannot be achieved because of a blocking
        ! layer (thinner than its minimum thickness) in the interior ocean,
        ! move interface k-1 (and k-2 if necessary) upward
        ! Only work on layers that are sufficiently far from the fixed near-surface layers.
        if ((h_hat >= 0.0) .and. (k > fixlay+2) .and. (p_int(k-1) > dp0cum(k-1) + tenm)) then

          ! Only act if interface k-1 is near the bottom or layer k-2 could donate water.
          if ( (p_int(nk+1) - p_int(k-1) < thkbot) .or. &
               (h_col(k-2) > qqmx*dp0ij(k-2)) ) then
            ! Determine how much water layer k-2 could supply without becoming too thin.
            if (k == fixlay+3) then
              ! Treat layer k-2 as fixed.
              h_hat2 = max(h_hat0 - h_hat, dp0ij(k-2) - h_col(k-2))
            else
              ! Maintain minimum thickess of layer k-2.
              h_hat2 = cushn(h_col(k-2) + (h_hat0 - h_hat), dp0ij(k-2)) - h_col(k-2)
            endif !fixlay+3:else

            if (h_hat2 < -onemm) then
              dh_cor = qhrlx(k-1) * max(h_hat2, -h_hat - h_col(k-1))
              h_col(k-2)  = h_col(k-2) + dh_cor
              h_col(k-1)  = h_col(k-1) - dh_cor
              dp_int(k-1) = dp_int(k-1) + dh_cor
              p_int(k-1)  = p_int(k-1) + dh_cor
              ! Recalculate how much layer k-1 could donate to layer k.
              h_hat = cushn(h_hat0 + h_col(k-1), dp0ij(k-1)) - h_col(k-1)
            elseif (k <= fixlay+3) then
              ! Do nothing.
            elseif (p_int(k-2) > dp0cum(k-2) + tenm .and. &
                   (p_int(nk+1) - p_int(k-2) < thkbot .or. &
                    h_col(k-3) > qqmx*dp0ij(k-3))) then

              ! Determine how much water layer k-3 could supply without becoming too thin.
              if (k == fixlay+4) then
                ! Treat layer k-3 as fixed.
                h_hat3 = max(h_hat0 - h_hat, dp0ij(k-3) - h_col(k-3))
              else
                ! Maintain minimum thickess of layer k-3.
                h_hat3 = cushn(h_col(k-3) + (h_hat0 - h_hat), dp0ij(k-3)) - h_col(k-3)
              endif !fixlay+4:else
              if (h_hat3 < -onemm) then
                ! Water is moved from layer k-3 to k-2, but do not dilute layer k-2 too much.
                dh_cor = qhrlx(k-2) * max(h_hat3, -h_col(k-2))
                h_col(k-3) = h_col(k-3) + dh_cor
                h_col(k-2) = h_col(k-2) - dh_cor
                dp_int(k-2) = dp_int(k-2) + dh_cor
                p_int(k-2) = p_int(k-2) + dh_cor

                ! Now layer k-2 might be able donate to layer k-1.
                h_hat2 = cushn(h_col(k-2) + (h_hat0 - h_hat), dp0ij(k-2)) - h_col(k-2)
                if (h_hat2 < -onemm) then
                  dh_cor = qhrlx(k-1) * (max(h_hat2, -h_hat - h_col(k-1)) )
                  h_col(k-2) = h_col(k-2) + dh_cor
                  h_col(k-1) = h_col(k-1) - dh_cor
                  dp_int(k-1) = dp_int(k-1) + dh_cor
                  p_int(k-1) = p_int(k-1) + dh_cor
                  ! Recalculate how much layer k-1 could donate to layer k.
                  h_hat = cushn(h_hat0 + h_col(k-1), dp0ij(k-1)) - h_col(k-1)
                endif !h_hat2
              endif !h_hat3
            endif !h_hat2:blocking
          endif ! Layer k-2 could move.
        endif ! blocking, i.e., h_hat >= 0, and far enough from the fixed layers to permit a change.

        if (h_hat < 0.0) then
          ! entrain layer k-1 water into layer k, move interface up.
          dh_cor = qhrlx(k) * h_hat
          h_col(k-1) = h_col(k-1) + dh_cor
          h_col(k)   = h_col(k)   - dh_cor
          dp_int(k) = dp_int(k) + dh_cor
          p_int(k) = p_int(k) + dh_cor
        endif !entrain

      endif  !too-dense adjustment
    endif

  ! In the original Hycom version, there is not a break between these two loops.
  enddo

  do k=fixlay+1,nk
    if (Rcv(k) < Rcv_tgt(k) - epsil) then   ! layer too light
      ! Water in layer k is too light, so try to dilute with water from layer k+1.
      ! Entrainment is not possible if layer k touches the bottom.
      if (p_int(k+1) < p_int(nk+1)) then  ! k<nk
        if ((Rcv(k+1) <= Rcv_tgt(k+1)) .or. &
            (p_int(k+1) <= dp0cum(k+1) + onem) .or. &
            (h_col(k) < h_col(k+1))) then
          ! If layer k+1 is too dense, there is a conflict in the direction the
          ! inteface between them should move, so thicken the thinner of the two.

          if     ((Rcv(k+1) - Rcv_tgt(k)) <= epsil) then
            ! layer k+1 is far too light, so take the entire layer
            ! Because this code is working downward, this flux does not accumulate across
            ! successive layers.
            h_hat = h_col(k+1)
          else
            q_frac = (Rcv_tgt(k) - Rcv(k)) / (Rcv(k+1) - Rcv_tgt(k)) ! 0 < q_frac <= 1
            h_hat = q_frac*h_col(k)
          endif

          ! If layer k+1, or layer k+2, does not touch the bottom, maintain minimum
          ! thicknesses of layers k and k+1 as much as possible. otherwise, permit
          ! layers to collapse to zero thickness at the bottom.
          if (p_int(min(k+3,nk+1)) < p_int(nk+1)) then
            if (p_int(nk+1) - p_int(k) > dp0ij(k) + dp0ij(k+1)) then
              h_hat = h_col(k+1) - cushn(h_col(k+1) - h_hat, dp0ij(k+1))
            endif
            ! Try to bring layer layer k up to its minimum thickness.
            h_hat = max(h_hat, dp0ij(k) - h_col(k))
            ! Do not drive layer k+1 below its minimum thickness or take more than half of it.
            h_hat = min(h_hat, max(0.5*h_col(k+1), h_col(k+1) - dp0ij(k+1)) )
          else
            ! Layers that touch the bottom can lose their entire contents.
            h_hat = min(h_col(k+1), h_hat)
          endif !p.k+2<p.nk+1

          if (h_hat > 0.0) then
            ! Entrain layer k+1 water into layer k.
            dh_cor = qhrlx(k+1) * h_hat
            h_col(k)   = h_col(k)   + dh_cor
            h_col(k+1) = h_col(k+1) - dh_cor
            dp_int(k+1) = dp_int(k+1) + dh_cor
            p_int(k+1) = p_int(k+1) + dh_cor
          endif !entrain

        endif !too-light adjustment
      endif !above bottom
    endif !too light

    ! If layer above is still too thin, move interface down.
    dh_cor = min(qhrlx(k-1) * min(dp0ij(k-1) - h_col(k-1), p_int(nk+1) - p_int(k)), h_col(k))
    if (dh_cor > 0.0) then
      h_col(k-1) = h_col(k-1) + dh_cor
      h_col(k)   = h_col(k)   - dh_cor
      dp_int(k) = dp_int(k) + dh_cor
      p_int(k) = p_int(k) + dh_cor
    endif

  enddo !k  Hybrid vertical coordinate relocation moving interface downward

  if (trap_errors) then
    ! Verify that everything is consistent.
    do k=1,nk
      if (abs((h_col(k) - h_in(k)) + (dp_int(K) - dp_int(K+1))) > 1.0e-13*max(p_int(nk+1), onem)) then
        write(mesg, '("k ",i4," h ",es13.4," h_in ",es13.4, " dp ",2es13.4," err ",es13.4)') &
              k, h_col(k), h_in(k), dp_int(K), dp_int(K+1), (h_col(k) - h_in(k)) + (dp_int(K) - dp_int(K+1))
        call MOM_error(FATAL, "Mismatched thickness changes in hybgen_regrid: "//trim(mesg))
      endif
      if (h_col(k) < 0.0) then ! Could instead do: -1.0e-15*max(p_int(nk+1), onem)) then
        write(mesg, '("k ",i4," h ",es13.4," h_in ",es13.4, " dp ",2es13.4, " fixlay ",i4)') &
              k, h_col(k), h_in(k), dp_int(K), dp_int(K+1), fixlay
        call MOM_error(FATAL, "Significantly negative final thickness in hybgen_regrid: "//trim(mesg))
      endif
    enddo
    do K=1,nk+1
      if (abs(dp_int(K) - (p_int(K) - pres_in(K))) > 1.0e-13*max(p_int(nk+1), onem)) then
        call MOM_error(FATAL, "Mismatched interface height changes in hybgen_regrid.")
      endif
    enddo
  endif

end subroutine hybgen_column_regrid

end module MOM_hybgen_regrid

! This code was translated in 2022 from the HYCOM hybgen code, which was primarily developed
! between 2000 and 2015, with some minor subsequent changes and bug fixes.
