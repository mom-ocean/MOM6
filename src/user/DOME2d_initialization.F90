!> Initialization of the 2D DOME experiment with density water initialized on a coastal shelf.
module DOME2d_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_ALE_sponge, only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

implicit none ; private

#include <MOM_memory.h>

! Public functions
public DOME2d_initialize_topography
public DOME2d_initialize_thickness
public DOME2d_initialize_temperature_salinity
public DOME2d_initialize_sponges

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

character(len=40) :: mdl = "DOME2D_initialization" !< This module's name.

contains

!> Initialize topography with a shelf and slope in a 2D domain
subroutine DOME2d_initialize_topography( D, G, param_file, max_depth )
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth in the units of depth_max
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth in arbitrary units

  ! Local variables
  integer :: i, j
  real    :: x, bay_depth, l1, l2
  real    :: dome2d_width_bay, dome2d_width_bottom, dome2d_depth_bay
  ! This include declares and sets the variable "version".
# include "version_variable.h"

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DOME2D_SHELF_WIDTH", dome2d_width_bay, &
                 'Width of shelf, as fraction of domain, in 2d DOME configuration.', &
                 units='nondim',default=0.1)
  call get_param(param_file, mdl, "DOME2D_BASIN_WIDTH", dome2d_width_bottom, &
                 'Width of deep ocean basin, as fraction of domain, in 2d DOME configuration.', &
                 units='nondim',default=0.3)
  call get_param(param_file, mdl, "DOME2D_SHELF_DEPTH", dome2d_depth_bay, &
                 'Depth of shelf, as fraction of basin depth, in 2d DOME configuration.', &
                 units='nondim',default=0.2)

  ! location where downslope starts
  l1 = dome2d_width_bay

  ! location where downslope reaches maximum depth
  l2 = 1.0 - dome2d_width_bottom

  bay_depth = dome2d_depth_bay

  do j=G%jsc,G%jec ; do i=G%isc,G%iec

    ! Compute normalized zonal coordinate
    x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon

    if ( x <= l1 ) then
      D(i,j) = bay_depth * max_depth
    elseif (( x > l1 ) .and. ( x < l2 )) then
      D(i,j) = bay_depth * max_depth + (1.0-bay_depth) * max_depth * &
               ( x - l1 ) / (l2 - l1)
    else
      D(i,j) = max_depth
    endif

  enddo ; enddo

end subroutine DOME2d_initialize_topography

!> Initialize thicknesses according to coordinate mode
subroutine DOME2d_initialize_thickness ( h, G, GV, US, param_file, just_read_params )
  type(ocean_grid_type),   intent(in)  :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in)  :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2].
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.

  ! Local variables
  real :: e0(SZK_(GV))     ! The resting interface heights, in depth units [Z ~> m], usually
                           ! negative because it is positive upward.
  real :: eta1D(SZK_(GV)+1)! Interface height relative to the sea surface
                           ! positive upward, in depth units [Z ~> m].
  integer :: i, j, k, is, ie, js, je, nz
  real    :: x
  real    :: delta_h
  real    :: min_thickness
  real    :: dome2d_width_bay, dome2d_width_bottom, dome2d_depth_bay
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) &
    call MOM_mesg("MOM_initialization.F90, DOME2d_initialize_thickness: setting thickness")

  call get_param(param_file, mdl,"MIN_THICKNESS",min_thickness, &
                 default=1.e-3, units="m", do_not_log=.true., scale=US%m_to_Z)
  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=.true.)
  call get_param(param_file, mdl, "DOME2D_SHELF_WIDTH", dome2d_width_bay, &
                 default=0.1, do_not_log=.true.)
  call get_param(param_file, mdl, "DOME2D_BASIN_WIDTH", dome2d_width_bottom, &
                 default=0.3, do_not_log=.true.)
  call get_param(param_file, mdl, "DOME2D_SHELF_DEPTH", dome2d_depth_bay, &
                 default=0.2, do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  ! WARNING: this routine specifies the interface heights so that the last layer
  !          is vanished, even at maximum depth. In order to have a uniform
  !          layer distribution, use this line of code within the loop:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz)
  !          To obtain a thickness distribution where the last layer is
  !          vanished and the other thicknesses uniformly distributed, use:
  !          e0(k) = -G%max_depth * real(k-1) / real(nz-1)
  do k=1,nz
    e0(k) = -G%max_depth * real(k-1) / real(nz)
  enddo

  select case ( coordinateMode(verticalCoordinate) )

    case ( REGRIDDING_LAYER, REGRIDDING_RHO )

      do j=js,je ; do i=is,ie
        eta1D(nz+1) = -G%bathyT(i,j)
        do k=nz,1,-1
          eta1D(k) = e0(k)
          if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_Z)) then
            eta1D(k) = eta1D(k+1) + GV%Angstrom_Z
            h(i,j,k) = GV%Angstrom_H
          else
            h(i,j,k) = GV%Z_to_H * (eta1D(k) - eta1D(k+1))
          endif
        enddo

        x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon
        if ( x <= dome2d_width_bay ) then
          h(i,j,1:nz-1) = GV%Angstrom_H
          h(i,j,nz) = GV%Z_to_H * dome2d_depth_bay * G%max_depth - (nz-1) * GV%Angstrom_H
        endif

      enddo ; enddo

 !  case ( IC_RHO_C )
 !
 !    do j=js,je ; do i=is,ie
 !       eta1D(nz+1) = -G%bathyT(i,j)
 !       do k=nz,1,-1
 !         eta1D(k) = e0(k)
 !         if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
 !           eta1D(k) = eta1D(k+1) + min_thickness
 !           h(i,j,k) = GV%Z_to_H * min_thickness
 !         else
 !           h(i,j,k) = GV%Z_to_H * (eta1D(k) - eta1D(k+1))
 !         endif
 !       enddo
 !
 !       x = G%geoLonT(i,j) / G%len_lon
 !       if ( x <= dome2d_width_bay ) then
 !         h(i,j,1:nz-1) = GV%Z_to_H * min_thickness
 !         h(i,j,nz) = GV%Z_to_H * (dome2d_depth_bay * G%max_depth - (nz-1) * min_thickness)
 !       endif
 !
 !    enddo ; enddo

    case ( REGRIDDING_ZSTAR )

      do j=js,je ; do i=is,ie
        eta1D(nz+1) = -G%bathyT(i,j)
        do k=nz,1,-1
          eta1D(k) = e0(k)
          if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
            eta1D(k) = eta1D(k+1) + min_thickness
            h(i,j,k) = GV%Z_to_H * min_thickness
          else
            h(i,j,k) = GV%Z_to_H * (eta1D(k) - eta1D(k+1))
          endif
        enddo
      enddo ; enddo

    case ( REGRIDDING_SIGMA )
      do j=js,je ; do i=is,ie
        h(i,j,:) = GV%Z_to_H*G%bathyT(i,j) / nz
      enddo ; enddo

    case default
      call MOM_error(FATAL,"dome2d_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

end subroutine DOME2d_initialize_thickness


!> Initialize temperature and salinity in the 2d DOME configuration
subroutine DOME2d_initialize_temperature_salinity ( T, S, h, G, GV, param_file, &
                     eqn_of_state, just_read_params)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T !< Potential temperature [degC]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S !< Salinity [ppt]
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h !< Layer thickness [H ~> m or kg m-2]
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  type(EOS_type),                            pointer     :: eqn_of_state !< Equation of state structure
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing T & S.

  ! Local variables
  integer   :: i, j, k, is, ie, js, je, nz
  real      :: x
  integer   :: index_bay_z
  real      :: delta_S, delta_T
  real      :: S_ref, T_ref;        ! Reference salinity and temperature within surface layer
  real      :: S_range, T_range;    ! Range of salinities and temperatures over the vertical
  real      :: xi0, xi1
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40) :: verticalCoordinate
  real    :: dome2d_width_bay, dome2d_width_bottom, dome2d_depth_bay

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=.true.)
  call get_param(param_file, mdl, "DOME2D_SHELF_WIDTH", dome2d_width_bay, &
                 default=0.1, do_not_log=.true.)
  call get_param(param_file, mdl, "DOME2D_BASIN_WIDTH", dome2d_width_bottom, &
                 default=0.3, do_not_log=.true.)
  call get_param(param_file, mdl, "DOME2D_SHELF_DEPTH", dome2d_depth_bay, &
                 default=0.2, do_not_log=.true.)
  call get_param(param_file, mdl, "S_REF", S_ref, 'Reference salinity', &
                 default=35.0, units='1e-3', do_not_log=just_read)
  call get_param(param_file, mdl,"T_REF",T_ref,'Reference temperature', units='degC', &
                fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"S_RANGE",S_range,'Initial salinity range', &
                units='1e-3', default=2.0, do_not_log=just_read)
  call get_param(param_file, mdl,"T_RANGE",T_range,'Initial temperature range', &
                units='1e-3', default=0.0, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  T(:,:,:) = 0.0
  S(:,:,:) = 0.0

  ! Linear salinity profile

  select case ( coordinateMode(verticalCoordinate) )

    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA )

      do j=js,je ; do i=is,ie
        xi0 = 0.0
        do k = 1,nz
          xi1 = xi0 + (GV%H_to_Z * h(i,j,k)) / G%max_depth
          S(i,j,k) = 34.0 + 0.5 * S_range * (xi0 + xi1)
          xi0 = xi1
        enddo
      enddo ; enddo

    case ( REGRIDDING_RHO )

      do j=js,je ; do i=is,ie
        xi0 = 0.0
        do k = 1,nz
          xi1 = xi0 + (GV%H_to_Z * h(i,j,k)) / G%max_depth
          S(i,j,k) = 34.0 + 0.5 * S_range * (xi0 + xi1)
          xi0 = xi1
        enddo
        x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon
        if ( x <= dome2d_width_bay ) then
          S(i,j,nz) = 34.0 + S_range
        endif
      enddo ; enddo

    case ( REGRIDDING_LAYER )

      delta_S = S_range / ( G%ke - 1.0 )
      S(:,:,1) = S_ref
      do k = 2,G%ke
        S(:,:,k) = S(:,:,k-1) + delta_S
      enddo

    case default
      call MOM_error(FATAL,"dome2d_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

  ! Modify salinity and temperature when z coordinates are used
  if ( coordinateMode(verticalCoordinate) == REGRIDDING_ZSTAR ) then
    index_bay_z = Nint ( dome2d_depth_bay * G%ke )
    do j = G%jsc,G%jec ; do i = G%isc,G%iec
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon
      if ( x <= dome2d_width_bay ) then
        S(i,j,1:index_bay_z) = S_ref + S_range; ! Use for z coordinates
        T(i,j,1:index_bay_z) = 1.0;             ! Use for z coordinates
      endif
    enddo ; enddo ! i and j loops
  endif ! Z initial conditions

  ! Modify salinity and temperature when sigma coordinates are used
  if ( coordinateMode(verticalCoordinate) == REGRIDDING_SIGMA ) then
    do i = G%isc,G%iec ; do j = G%jsc,G%jec
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon
      if ( x <= dome2d_width_bay ) then
        S(i,j,1:G%ke) = S_ref + S_range;    ! Use for sigma coordinates
        T(i,j,1:G%ke) = 1.0;                ! Use for sigma coordinates
      endif
    enddo ; enddo
  endif

  ! Modify temperature when rho coordinates are used
  T(G%isc:G%iec,G%jsc:G%jec,1:G%ke) = 0.0
  if (( coordinateMode(verticalCoordinate) == REGRIDDING_RHO ) .or. &
      ( coordinateMode(verticalCoordinate) == REGRIDDING_LAYER )) then
    do i = G%isc,G%iec ; do j = G%jsc,G%jec
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon
      if ( x <= dome2d_width_bay ) then
        T(i,j,G%ke) = 1.0
      endif
    enddo ; enddo
  endif

end subroutine DOME2d_initialize_temperature_salinity

!> Set up sponges in 2d DOME configuration
subroutine DOME2d_initialize_sponges(G, GV, US, tv, param_file, use_ALE, CSp, ACSp)
  type(ocean_grid_type),   intent(in) :: G  !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv !< Thermodynamics structure
  type(param_file_type),   intent(in) :: param_file !< Parameter file structure
  logical,                 intent(in) :: use_ALE !< If true, indicates model is in ALE mode
  type(sponge_CS),         pointer    :: CSp !< Layer-mode sponge structure
  type(ALE_sponge_CS),     pointer    :: ACSp !< ALE-mode sponge structure
  ! Local variables
  real :: T(SZI_(G),SZJ_(G),SZK_(G))   ! A temporary array for temp [degC]
  real :: S(SZI_(G),SZJ_(G),SZK_(G))   ! A temporary array for salt [ppt]
  real :: h(SZI_(G),SZJ_(G),SZK_(G))   ! A temporary array for thickness [H ~> m or kg m-2].
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for thickness [Z ~> m]
  real :: Idamp(SZI_(G),SZJ_(G))       ! The inverse damping rate [T-1 ~> s-1].
  real :: S_ref, T_ref                 ! Reference salinity and temerature within surface layer
  real :: S_range, T_range             ! Range of salinities and temperatures over the vertical
  real :: e0(SZK_(G)+1)             ! The resting interface heights [Z ~> m],
                                    ! usually negative because it is positive upward.
  real :: eta1D(SZK_(G)+1)          ! Interface height relative to the sea surface
                                    ! positive upward [Z ~> m].
  real :: d_eta(SZK_(G))            ! The layer thickness in a column [Z ~> m].
  real :: dome2d_width_bay, dome2d_width_bottom, dome2d_depth_bay
  real :: dome2d_west_sponge_time_scale, dome2d_east_sponge_time_scale ! Sponge timescales [T ~> s]
  real :: dome2d_west_sponge_width, dome2d_east_sponge_width
  real :: dummy1, x, z
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call get_param(param_file, mdl, "DOME2D_WEST_SPONGE_TIME_SCALE", dome2d_west_sponge_time_scale, &
                 'The time-scale on the west edge of the domain for restoring T/S '//&
                 'in the sponge. If zero, the western sponge is disabled', &
                 units='s', default=0., scale=US%s_to_T)
  call get_param(param_file, mdl, "DOME2D_EAST_SPONGE_TIME_SCALE", dome2d_east_sponge_time_scale, &
                 'The time-scale on the east edge of the domain for restoring T/S '//&
                 'in the sponge. If zero, the eastern sponge is disabled', &
                 units='s', default=0., scale=US%s_to_T)
  call get_param(param_file, mdl, "DOME2D_WEST_SPONGE_WIDTH", dome2d_west_sponge_width, &
                 'The fraction of the domain in which the western sponge for restoring T/S '//&
                 'is active.', &
                 units='nondim', default=0.1)
  call get_param(param_file, mdl, "DOME2D_EAST_SPONGE_WIDTH", dome2d_east_sponge_width, &
                 'The fraction of the domain in which the eastern sponge for restoring T/S '//&
                 'is active.', &
                 units='nondim', default=0.1)

  ! Return if sponges are not in use
  if (dome2d_west_sponge_time_scale <= 0. .and. dome2d_east_sponge_time_scale <= 0.) return

  if (associated(CSp)) call MOM_error(FATAL, &
     "DOME2d_initialize_sponges called with an associated control structure.")
  if (associated(ACSp)) call MOM_error(FATAL, &
     "DOME2d_initialize_sponges called with an associated ALE-sponge control structure.")

  call get_param(param_file, mdl, "DOME2D_SHELF_WIDTH", dome2d_width_bay, &
                 default=0.1, do_not_log=.true.)
  call get_param(param_file, mdl, "DOME2D_BASIN_WIDTH", dome2d_width_bottom, &
                 default=0.3, do_not_log=.true.)
  call get_param(param_file, mdl, "DOME2D_SHELF_DEPTH", dome2d_depth_bay, &
                 default=0.2, do_not_log=.true.)
  call get_param(param_file, mdl, "S_REF", S_ref, default=35.0)
  call get_param(param_file, mdl, "T_REF", T_ref)
  call get_param(param_file, mdl, "S_RANGE", S_range, default=2.0)
  call get_param(param_file, mdl, "T_RANGE", T_range, default=0.0)


  ! Set the inverse damping rate as a function of position
  Idamp(:,:) = 0.0
  do j=js,je ; do i=is,ie
    if (G%mask2dT(i,j) > 0.) then ! Only set damping rate for wet points
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon ! Non-dimensional position within domain (0,1)
      if ( dome2d_west_sponge_time_scale > 0. .and. x < dome2d_west_sponge_width ) then
        ! Within half the shelf width from the left edge
        dummy1 = 1. - x / dome2d_west_sponge_width
        Idamp(i,j) = 1./dome2d_west_sponge_time_scale * max(0., min(1., dummy1))
      elseif ( dome2d_east_sponge_time_scale > 0. .and. x > ( 1. - dome2d_east_sponge_width ) ) then
        ! Within a quarter of the basin width from the right
        dummy1 = 1. - ( 1. - x ) / dome2d_east_sponge_width
        Idamp(i,j) = 1./dome2d_east_sponge_time_scale * max(0., min(1., dummy1))
      else
        Idamp(i,j) = 0.
      endif
    else
      Idamp(i,j) = 0.
    endif
  enddo ; enddo


  if (use_ALE) then

    ! Construct a grid (somewhat arbitrarily) to describe  the sponge T/S on
    do k=1,nz
     e0(k) = -G%max_depth * ( real(k-1) / real(nz) )
    enddo
    e0(nz+1) = -G%max_depth
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_Z)) then
          eta1D(k) = eta1D(k+1) + GV%Angstrom_Z
          h(i,j,k) = GV%Angstrom_H
        else
          h(i,j,k) = GV%Z_to_H * (eta1D(k) - eta1D(k+1))
        endif
      enddo
    enddo ; enddo
    ! Store the grid on which the T/S sponge data will reside
    call initialize_ALE_sponge(Idamp, G, param_file, ACSp, h, nz)

    ! Construct temperature and salinity on the arbitrary grid
    T(:,:,:) = 0.0 ; S(:,:,:) = 0.0
    do j=js,je ; do i=is,ie
      z = -G%bathyT(i,j)
      do k = nz,1,-1
        z = z + 0.5 * GV%H_to_Z * h(i,j,k) ! Position of the center of layer k
        S(i,j,k) = 34.0 - 1.0 * (z / (G%max_depth))
        if ( ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon < dome2d_west_sponge_width ) &
          S(i,j,k) = S_ref + S_range
        z = z + 0.5 *  GV%H_to_Z * h(i,j,k) ! Position of the interface k
      enddo
    enddo ; enddo

    if ( associated(tv%T) ) then
      call set_up_ALE_sponge_field(T, G, tv%T, ACSp)
    endif
    if ( associated(tv%S) ) then
      call set_up_ALE_sponge_field(S, G, tv%S, ACSp)
    endif

  else

    ! Construct interface heights to restore toward
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(K) = -G%max_depth * real(k-1) / real(nz)
        if (eta1D(K) < (eta1D(K+1) + GV%Angstrom_Z)) then
          eta1D(K) = eta1D(K+1) + GV%Angstrom_Z
          d_eta(k) = GV%Angstrom_Z
        else
          d_eta(k) = (eta1D(K) - eta1D(K+1))
        endif
      enddo

      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon
      if ( x <= dome2d_width_bay ) then
        do k=1,nz-1 ; d_eta(k) = GV%Angstrom_Z ; enddo
        d_eta(nz) = dome2d_depth_bay * G%max_depth - (nz-1) * GV%Angstrom_Z
      endif

      eta(i,j,nz+1) = -G%bathyT(i,j)
      do K=nz,1,-1
        eta(i,j,K) = eta(i,j,K+1) + d_eta(k)
      enddo
    enddo ; enddo
    call initialize_sponge(Idamp, eta, G, param_file, CSp, GV)

  endif

end subroutine DOME2d_initialize_sponges

end module DOME2d_initialization
