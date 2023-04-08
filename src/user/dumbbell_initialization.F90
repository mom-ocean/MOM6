!> Configures the model for the idealized dumbbell test case.
module dumbbell_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains, only : sum_across_PEs
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_interface_heights, only : dz_to_thickness, dz_to_thickness_simple
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA, REGRIDDING_HYCOM1
use MOM_ALE_sponge,    only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge

implicit none ; private

#include <MOM_memory.h>

character(len=40) :: mdl = "dumbbell_initialization" !< This module's name.

public dumbbell_initialize_topography
public dumbbell_initialize_thickness
public dumbbell_initialize_temperature_salinity
public dumbbell_initialize_sponges

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Initialization of topography.
subroutine dumbbell_initialize_topography( D, G, param_file, max_depth )
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth [Z ~> m]

  ! Local variables
  real    :: x, y   ! Fractional x- and y- positions [nondim]
  real    :: dblen  ! Lateral length scale for dumbbell [km] or [m]
  real    :: dbfrac ! Meridional fraction for narrow part of dumbbell [nondim]
  logical :: dbrotate ! If true, rotate this configuration
  integer :: i, j

  call get_param(param_file, mdl, "DUMBBELL_LEN", dblen, &
                'Lateral Length scale for dumbbell.', &
                 units=G%x_ax_unit_short, default=600., do_not_log=.false.)
  call get_param(param_file, mdl, "DUMBBELL_FRACTION", dbfrac, &
                'Meridional fraction for narrow part of dumbbell.', &
                 units='nondim', default=0.5, do_not_log=.false.)
  call get_param(param_file, mdl, "DUMBBELL_ROTATION", dbrotate, &
                'Logical for rotation of dumbbell domain.', &
                 default=.false., do_not_log=.false.)

  if (G%x_axis_units(1:1) == 'm') then
    dblen = dblen*1.e3
  endif

  if (dbrotate) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      ! Compute normalized zonal coordinates (x,y=0 at center of domain)
      x = ( G%geoLonT(i,j) ) / G%len_lon
      y = ( G%geoLatT(i,j)  ) / dblen
      D(i,j) = G%max_depth
      if ((y>=-0.25 .and. y<=0.25) .and. (x <= -0.5*dbfrac .or. x >= 0.5*dbfrac)) then
        D(i,j) = 0.0
      endif
    enddo ; enddo
  else
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      ! Compute normalized zonal coordinates (x,y=0 at center of domain)
      x = ( G%geoLonT(i,j) ) / dblen
      y = ( G%geoLatT(i,j)  ) / G%len_lat
      D(i,j) = G%max_depth
      if ((x>=-0.25 .and. x<=0.25) .and. (y <= -0.5*dbfrac .or. y >= 0.5*dbfrac)) then
        D(i,j) = 0.0
      endif
    enddo ; enddo
  endif

end subroutine dumbbell_initialize_topography

!> Initializes the layer thicknesses to be uniform in the dumbbell test case
subroutine dumbbell_initialize_thickness ( h, depth_tot, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot   !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

  real :: e0(SZK_(GV)+1)  ! The resting interface heights [Z ~> m], usually
                          ! negative because it is positive upward.
  real :: eta1D(SZK_(GV)+1) ! Interface height relative to the sea surface
                          ! positive upward [Z ~> m].
  real :: min_thickness   ! The minimum layer thicknesses [Z ~> m].
  real :: S_ref           ! A default value for salinities [S ~> ppt].
  real :: S_surf          ! The surface salinity [S ~> ppt]
  real :: S_range         ! The range of salinities in this test case [S ~> ppt]
  real :: S_light, S_dense ! The lightest and densest salinities in the sponges [S ~> ppt].
  real :: eta_IC_quanta   ! The granularity of quantization of initial interface heights [Z-1 ~> m-1].
  real :: x               ! Along-channel position in the axis units [m] or [km] or [deg]
  logical :: dbrotate     ! If true, rotate the domain.
  logical :: use_ALE      ! True if ALE is being used, False if in layered mode

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=20) :: verticalCoordinate
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("dumbbell_initialization.F90, dumbbell_initialize_thickness: setting thickness")

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl,"MIN_THICKNESS", min_thickness, &
                'Minimum thickness for layer', &
                 units='m', default=1.0e-3, scale=US%m_to_Z, do_not_log=just_read)
  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl, "USE_REGRIDDING", use_ALE, default=.false., do_not_log=.true.)
  if (.not. use_ALE) verticalCoordinate = "LAYER"

  ! WARNING: this routine specifies the interface heights so that the last layer
  !          is vanished, even at maximum depth. In order to have a uniform
  !          layer distribution, use this line of code within the loop:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz)
  !          To obtain a thickness distribution where the last layer is
  !          vanished and the other thicknesses uniformly distributed, use:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz-1)
  !do k=1,nz+1
  !  e0(k) = -G%max_depth * real(k-1) / real(nz)
  !enddo

  select case ( coordinateMode(verticalCoordinate) )
  case ( REGRIDDING_LAYER) ! Initial thicknesses for isopycnal coordinates
    call get_param(param_file, mdl, "DUMBBELL_ROTATION", dbrotate, &
                'Logical for rotation of dumbbell domain.', &
                 default=.false., do_not_log=just_read)
    do j=js,je
      do i=is,ie
        ! Compute normalized zonal coordinates (x,y=0 at center of domain)
        if (dbrotate) then
          ! This is really y in the rotated case
          x = G%geoLatT(i,j)
        else
          x = G%geoLonT(i,j)
        endif
        eta1D(1) = 0.0
        eta1D(nz+1) = -depth_tot(i,j)
        if (x<0.0) then
          do k=nz,2, -1
            eta1D(k) =  eta1D(k+1) + min_thickness
          enddo
        else
          do k=2,nz
            eta1D(k) =  eta1D(k-1) - min_thickness
          enddo
        endif
        do k=1,nz
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        enddo
      enddo
    enddo

  case ( REGRIDDING_RHO, REGRIDDING_HYCOM1) ! Initial thicknesses for isopycnal coordinates
    call get_param(param_file, mdl, "INITIAL_SSS", S_surf, &
                   units='1e-3', default=34., scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(param_file, mdl, "INITIAL_S_RANGE", S_range, &
                   units='1e-3', default=2., scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(param_file, mdl, "S_REF", S_ref, &
                   units='1e-3', default=35.0, scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(param_file, mdl, "TS_RANGE_S_LIGHT", S_light, &
                   units='1e-3', default=US%S_to_ppt*S_Ref, scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(param_file, mdl, "TS_RANGE_S_DENSE", S_dense, &
                   units='1e-3', default=US%S_to_ppt*S_Ref, scale=US%ppt_to_S, do_not_log=.true.)
    call get_param(param_file, mdl, "INTERFACE_IC_QUANTA", eta_IC_quanta, &
                   "The granularity of initial interface height values "//&
                   "per meter, to avoid sensivity to order-of-arithmetic changes.", &
                   default=2048.0, units="m-1", scale=US%Z_to_m, do_not_log=just_read)
    if (just_read) return ! All run-time parameters have been read, so return.

    do K=1,nz+1
      ! Salinity of layer k is S_light + (k-1)/(nz-1) * (S_dense - S_light)
      ! Salinity of interface K is S_light + (K-3/2)/(nz-1) * (S_dense - S_light)
      ! Salinity at depth z should be S(z) = S_surf - S_range * z/max_depth
      ! Equating: S_surf - S_range * z/max_depth = S_light + (K-3/2)/(nz-1) * (S_dense - S_light)
      ! Equating: - S_range * z/max_depth = S_light - S_surf + (K-3/2)/(nz-1) * (S_dense - S_light)
      ! Equating: z/max_depth = - ( S_light - S_surf + (K-3/2)/(nz-1) * (S_dense - S_light) ) / S_range
      e0(K) = - G%max_depth * ( ( S_light  - S_surf ) + ( S_dense - S_light ) * &
                ( (real(K)-1.5) / real(nz-1) ) ) / S_range
      ! Force round numbers ... the above expression has irrational factors ...
      if (eta_IC_quanta > 0.0) &
        e0(K) = nint(eta_IC_quanta*e0(K)) / eta_IC_quanta
      e0(K) = min(real(1-K)*GV%Angstrom_Z, e0(K)) ! Bound by surface
      e0(K) = max(-G%max_depth, e0(K)) ! Bound by bottom
    enddo
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -depth_tot(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_Z)) then
          eta1D(k) = eta1D(k+1) + GV%Angstrom_Z
          h(i,j,k) = GV%Angstrom_Z
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  case ( REGRIDDING_ZSTAR )                       ! Initial thicknesses for z coordinates
    if (just_read) return ! All run-time parameters have been read, so return.
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -depth_tot(i,j)
      do k=nz,1,-1
        eta1D(k) = -G%max_depth * real(k-1) / real(nz)
        if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
          eta1D(k) = eta1D(k+1) + min_thickness
          h(i,j,k) = min_thickness
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

  case ( REGRIDDING_SIGMA )             ! Initial thicknesses for sigma coordinates
    if (just_read) return ! All run-time parameters have been read, so return.
    do j=js,je ; do i=is,ie
      h(i,j,:) = depth_tot(i,j) / real(nz)
    enddo ; enddo

end select

end subroutine dumbbell_initialize_thickness

!> Initial values for temperature and salinity for the dumbbell test case
subroutine dumbbell_initialize_temperature_salinity ( T, S, h, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T !< Potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S !< Salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h !< Layer thickness [Z ~> m]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  logical,                                   intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing h.

  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz
  real    :: S_surf    ! The surface salinity [S ~> ppt]
  real    :: S_range   ! The range of salinities in this test case [S ~> ppt]
  real    :: T_surf    ! The surface temperature [C ~> degC]
  real    :: x         ! The fractional position in the domain [nondim]
  real    :: dblen     ! The size of the dumbbell test case [km] or [m]
  logical :: dbrotate  ! If true, rotate the domain.
  logical :: use_ALE   ! If false, use layer mode.
  character(len=20) :: verticalCoordinate, density_profile

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  ! layer mode
  call get_param(param_file, mdl, "USE_REGRIDDING", use_ALE, default=.false., do_not_log=.true.)
  if (.not. use_ALE) call MOM_error(FATAL,  "dumbbell_initialize_temperature_salinity: "//&
               "Please use 'fit' for 'TS_CONFIG' in the LAYER mode.")

  call get_param(param_file, mdl, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl, "INITIAL_DENSITY_PROFILE", density_profile, &
                 'Initial profile shape. Valid values are "linear", "parabolic" '// &
                 'and "exponential".', default='linear', do_not_log=just_read)
  call get_param(param_file, mdl, "DUMBBELL_T_SURF", T_surf, &
                 'Initial surface temperature in the DUMBBELL configuration', &
                 units='degC', default=20., scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "DUMBBELL_SREF", S_surf, &
                 'DUMBBELL REFERENCE SALINITY', &
                 units='1e-3', default=34., scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "DUMBBELL_S_RANGE", S_range, &
                 'DUMBBELL salinity range (right-left)', &
                 units='1e-3', default=2., scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "DUMBBELL_LEN", dblen, &
                 'Lateral Length scale for dumbbell ', &
                 units=G%x_ax_unit_short, default=600., do_not_log=just_read)
  call get_param(param_file, mdl, "DUMBBELL_ROTATION", dbrotate, &
                'Logical for rotation of dumbbell domain.', &
                 default=.false., do_not_log=just_read)

  if (G%x_axis_units(1:1) == 'm') then
    dblen = dblen*1.e3
  endif

  do j=G%jsc,G%jec
    do i=G%isc,G%iec
    ! Compute normalized zonal coordinates (x,y=0 at center of domain)
      if (dbrotate) then
        ! This is really y in the rotated case
        x = ( G%geoLatT(i,j) ) / dblen
      else
        x = ( G%geoLonT(i,j) ) / dblen
      endif
      do k=1,nz
        T(i,j,k) = T_surf
      enddo
      if (x>=0. ) then
        do k=1,nz
          S(i,j,k) = S_surf + 0.5*S_range
        enddo
      endif
      if (x<0. ) then
        do k=1,nz
          S(i,j,k) = S_surf - 0.5*S_range
        enddo
      endif

    enddo
  enddo

end subroutine dumbbell_initialize_temperature_salinity

!> Initialize the restoring sponges for the dumbbell test case
subroutine dumbbell_initialize_sponges(G, GV, US, tv, h_in, depth_tot, param_file, use_ALE, CSp, ACSp)
  type(ocean_grid_type),   intent(in) :: G !< Horizontal grid control structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid control structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv !< Thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h_in !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: param_file !< Parameter file structure
  logical,                 intent(in) :: use_ALE !< ALE flag
  type(sponge_CS),         pointer    :: CSp !< Layered sponge control structure pointer
  type(ALE_sponge_CS),     pointer    :: ACSp !< ALE sponge control structure pointer

  real :: sponge_time_scale  ! The damping time scale [T ~> s]

  real, dimension(SZI_(G),SZJ_(G)) :: Idamp ! inverse damping timescale [T-1 ~> s-1]
  real :: dz(SZI_(G),SZJ_(G),SZK_(GV)) ! Sponge thicknesses in height units [Z ~> m]
  real :: h(SZI_(G),SZJ_(G),SZK_(GV))  ! Sponge thicknesses [H ~> m or kg m-2]
  real :: S(SZI_(G),SZJ_(G),SZK_(GV))  ! Sponge salinities [S ~> ppt]
  real :: T(SZI_(G),SZJ_(G),SZK_(GV))  ! Sponge tempertures [C ~> degC], used only to convert thicknesses
                                       ! in non-Boussinesq mode
  real, dimension(SZK_(GV)+1) :: eta1D ! Interface positions for ALE sponge [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: eta ! A temporary array for interface heights [Z ~> m].

  integer :: i, j, k, nz
  real :: x              ! The fractional position in the domain [nondim]
  real :: dblen          ! The size of the dumbbell test case [km] or [m]
  real :: min_thickness  ! The minimum layer thickness [Z ~> m]
  real :: S_ref, S_range ! A reference salinity and the range of salinities in this test case [S ~> ppt]
  real :: T_surf         ! The surface temperature [C ~> degC]
  logical :: dbrotate    ! If true, rotate the domain.

  call get_param(param_file, mdl,"DUMBBELL_LEN",dblen, &
                'Lateral Length scale for dumbbell ', &
                 units='km', default=600., do_not_log=.true.)
  call get_param(param_file, mdl, "DUMBBELL_ROTATION", dbrotate, &
                'Logical for rotation of dumbbell domain.', &
                 default=.false., do_not_log=.true.)

  if (G%x_axis_units(1:1) == 'm') then
    dblen = dblen*1.e3
  endif

  nz = GV%ke

  call get_param(param_file, mdl, "DUMBBELL_SPONGE_TIME_SCALE", sponge_time_scale, &
                 "The time scale in the reservoir for restoring. If zero, the sponge is disabled.", &
                 units="s", default=0., scale=US%s_to_T)
  call get_param(param_file, mdl, "DUMBBELL_T_SURF", T_surf, &
                 'Initial surface temperature in the DUMBBELL configuration', &
                 units='degC', default=20., scale=US%degC_to_C, do_not_log=.true.)
  call get_param(param_file, mdl, "DUMBBELL_SREF", S_ref, &
                 'DUMBBELL REFERENCE SALINITY', &
                 units='1e-3', default=34., scale=US%ppt_to_S, do_not_log=.true.)
  call get_param(param_file, mdl, "DUMBBELL_S_RANGE", S_range, &
                 'DUMBBELL salinity range (right-left)', &
                 units='1e-3', default=2., scale=US%ppt_to_S, do_not_log=.true.)
  call get_param(param_file, mdl,"MIN_THICKNESS", min_thickness, &
                'Minimum thickness for layer', &
                 units='m', default=1.0e-3, scale=US%m_to_Z, do_not_log=.true.)

  ! no active sponges
  if (sponge_time_scale <= 0.) return

  ! everywhere is initially unsponged
  Idamp(:,:) = 0.0

  do j = G%jsc, G%jec
    do i = G%isc,G%iec
      if (G%mask2dT(i,j) > 0.) then
        ! nondimensional x position
        if (dbrotate) then
          ! This is really y in the rotated case
          x = ( G%geoLatT(i,j) ) / dblen
        else
          x = ( G%geoLonT(i,j) ) / dblen
        endif
        if (x > 0.25 .or. x < -0.25) then
          ! scale restoring by depth into sponge
          Idamp(i,j) = 1. / sponge_time_scale
        endif
      endif
    enddo
  enddo

  if (use_ALE) then
    ! construct a uniform grid for the sponge
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      eta1D(nz+1) =  depth_tot(i,j)
      do k=nz,1,-1
        eta1D(k) = -G%max_depth * real(k-1) / real(nz)
        if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
          eta1D(k) = eta1D(k+1) + min_thickness
          dz(i,j,k) = min_thickness
        else
          dz(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo ; enddo

    ! construct temperature and salinity for the sponge
    ! start with initial condition
    S(:,:,:) = 0.0
    T(:,:,:) = T_surf

    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      ! Compute normalized zonal coordinates (x,y=0 at center of domain)
      if (dbrotate) then
        ! This is really y in the rotated case
        x = ( G%geoLatT(i,j) ) / dblen
      else
        x = ( G%geoLonT(i,j) ) / dblen
      endif
      if (x>=0.25 ) then
        do k=1,nz
          S(i,j,k) = S_ref + 0.5*S_range
        enddo
      endif
      if (x<=-0.25 ) then
        do k=1,nz
          S(i,j,k) = S_ref - 0.5*S_range
        enddo
      endif
    enddo ; enddo

    ! Convert thicknesses from height units to thickness units
    if (associated(tv%eqn_of_state)) then
      call dz_to_thickness(dz, T, S, tv%eqn_of_state, h, G, GV, US)
    else
      call dz_to_thickness_simple(dz, h, G, GV, US, layer_mode=.true.)
    endif

    ! Store damping rates and the grid on which the T/S sponge data will reside
    call initialize_ALE_sponge(Idamp, G, GV, param_file, ACSp, h, nz)

    if (associated(tv%S)) call set_up_ALE_sponge_field(S, G, GV, tv%S, ACSp, 'salt', &
                          sp_long_name='salinity', sp_unit='g kg-1 s-1')
  else
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      eta(i,j,1) = 0.0
      do k=2,nz
        eta(i,j,k) = eta(i,j,k-1) - GV%H_to_Z * h_in(i,j,k-1)
      enddo
      eta(i,j,nz+1) = -depth_tot(i,j)
      do k=1,nz
        S(i,j,k)= tv%S(i,j,k)
      enddo
    enddo ; enddo

    !  This call sets up the damping rates and interface heights.
    !  This sets the inverse damping timescale fields in the sponges.    !
    call initialize_sponge(Idamp, eta, G, param_file, CSp, GV)

    !  The remaining calls to set_up_sponge_field can be in any order. !
    if ( associated(tv%S) ) call set_up_sponge_field(S, tv%S, G, GV, nz, CSp)
  endif

end subroutine dumbbell_initialize_sponges

end module dumbbell_initialization
