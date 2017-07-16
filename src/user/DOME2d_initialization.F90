module DOME2d_initialization

use MOM_ALE_sponge, only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : close_file, fieldtype, file_exists
use MOM_io, only : open_file, read_data, read_axis_data, SINGLE_FILE
use MOM_io, only : write_field, slasher, vardesc
use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
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

character(len=40) :: mod = "DOEM2D_initialization" !< This module's name.

contains

!> Initialize topography with a shelf and slope in a 2D domain
subroutine DOME2d_initialize_topography ( D, G, param_file, max_depth )
  ! Arguments
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m
  ! Local variables
  integer :: i, j
  real    :: x, bay_depth, l1, l2
  real    :: dome2d_width_bay, dome2d_width_bottom, dome2d_depth_bay

  call get_param(param_file, mod, "DOME2D_SHELF_WIDTH", dome2d_width_bay, &
                 'Width of shelf, as fraction of domain, in 2d DOME configuration.', &
                 units='nondim',default=0.1)
  call get_param(param_file, mod, "DOME2D_BASIN_WIDTH", dome2d_width_bottom, &
                 'Width of deep ocean basin, as fraction of domain, in 2d DOME configuration.', &
                 units='nondim',default=0.3)
  call get_param(param_file, mod, "DOME2D_SHELF_DEPTH", dome2d_depth_bay, &
                 'Depth of shelf, as fraction of basin depth, in 2d DOME configuration.', &
                 units='nondim',default=0.2)

  ! location where downslope starts
  l1 = dome2d_width_bay

  ! location where downslope reaches maximum depth
  l2 = 1.0 - dome2d_width_bottom

  bay_depth = dome2d_depth_bay

  do i=G%isc,G%iec
    do j=G%jsc,G%jec

      ! Compute normalized zonal coordinate
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon;

      if ( x .le. l1 ) then
        D(i,j) = bay_depth * max_depth
      else if (( x .gt. l1 ) .and. ( x .lt. l2 )) then
        D(i,j) = bay_depth * max_depth + (1.0-bay_depth) * max_depth * &
                 ( x - l1 ) / (l2 - l1)
      else
        D(i,j) = max_depth
      end if

    enddo
  enddo
end subroutine DOME2d_initialize_topography

!> Initialize thicknesses according to coordinate mode
subroutine DOME2d_initialize_thickness ( h, G, GV, param_file )
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: h !< Layer thicknesses
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  ! Local variables
  real :: e0(SZK_(G))     ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  integer :: i, j, k, is, ie, js, je, nz
  real    :: x
  real    :: delta_h
  real    :: min_thickness
  character(len=40) :: verticalCoordinate
  real    :: dome2d_width_bay, dome2d_width_bottom, dome2d_depth_bay

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call MOM_mesg("MOM_initialization.F90, DOME2d_initialize_thickness: setting thickness")

  call get_param(param_file,mod,"MIN_THICKNESS",min_thickness, &
                 default=1.e-3, do_not_log=.true.)
  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=.true.)
  call get_param(param_file, mod, "DOME2D_SHELF_WIDTH", dome2d_width_bay, &
                 default=0.1, do_not_log=.true.)
  call get_param(param_file, mod, "DOME2D_BASIN_WIDTH", dome2d_width_bottom, &
                 default=0.3, do_not_log=.true.)
  call get_param(param_file, mod, "DOME2D_SHELF_DEPTH", dome2d_depth_bay, &
                 default=0.2, do_not_log=.true.)

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
        eta1D(nz+1) = -1.0*G%bathyT(i,j)
        do k=nz,1,-1
          eta1D(k) = e0(k)
          if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_z)) then
            eta1D(k) = eta1D(k+1) + GV%Angstrom_z
            h(i,j,k) = GV%Angstrom_z
          else
            h(i,j,k) = eta1D(k) - eta1D(k+1)
          endif
        enddo

         x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon;
         if ( x .le. dome2d_width_bay ) then
           h(i,j,1:nz-1) = GV%Angstrom;
           h(i,j,nz) = dome2d_depth_bay * G%max_depth - (nz-1) * GV%Angstrom;
         end if

      end do ; end do

 !  case ( IC_RHO_C )
 !
 !    do j=js,je ; do i=is,ie
 !        eta1D(nz+1) = -1.0*G%bathyT(i,j)
 !        do k=nz,1,-1
 !          eta1D(k) = e0(k)
 !          if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
 !            eta1D(k) = eta1D(k+1) + min_thickness
 !            h(i,j,k) = min_thickness
 !          else
 !            h(i,j,k) = eta1D(k) - eta1D(k+1)
 !          endif
 !       enddo
 !
 !       x = G%geoLonT(i,j) / G%len_lon;
 !       if ( x .le. dome2d_width_bay ) then
 !         h(i,j,1:nz-1) = min_thickness;
 !         h(i,j,nz) = dome2d_depth_bay * G%max_depth - (nz-1) * min_thickness;
 !       end if
 !
 !    enddo ; enddo

    case ( REGRIDDING_ZSTAR )

      do j=js,je ; do i=is,ie
          eta1D(nz+1) = -1.0*G%bathyT(i,j)
          do k=nz,1,-1
            eta1D(k) = e0(k)
            if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
              eta1D(k) = eta1D(k+1) + min_thickness
              h(i,j,k) = min_thickness
            else
              h(i,j,k) = eta1D(k) - eta1D(k+1)
            endif
         enddo
      enddo ; enddo

    case ( REGRIDDING_SIGMA )
      do j=js,je ; do i=is,ie
        delta_h = G%bathyT(i,j) / nz;
        h(i,j,:) = delta_h;
      end do ; end do

    case default
      call MOM_error(FATAL,"dome2d_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

end subroutine DOME2d_initialize_thickness


!> Initialize temperature and salinity in the 2d DOME configuration
subroutine DOME2d_initialize_temperature_salinity ( T, S, h, G, param_file, eqn_of_state)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T !< Potential temperature (degC)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S !< Salinity (ppt)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h !< Layer thickness (m or Pa)
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  type(EOS_type),                            pointer     :: eqn_of_state !< Equation of state structure
  ! Local variables
  integer   :: i, j, k, is, ie, js, je, nz
  real      :: x;
  integer   :: index_bay_z;
  real      :: delta_S, delta_T;
  real      :: S_ref, T_ref;        ! Reference salinity and temperature within surface layer
  real      :: S_range, T_range;    ! Range of salinities and temperatures over the vertical
  real      :: xi0, xi1;
  character(len=40) :: verticalCoordinate
  real    :: dome2d_width_bay, dome2d_width_bottom, dome2d_depth_bay

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  call get_param(param_file,mod,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=.true.)
  call get_param(param_file, mod, "DOME2D_SHELF_WIDTH", dome2d_width_bay, &
                 default=0.1, do_not_log=.true.)
  call get_param(param_file, mod, "DOME2D_BASIN_WIDTH", dome2d_width_bottom, &
                 default=0.3, do_not_log=.true.)
  call get_param(param_file, mod, "DOME2D_SHELF_DEPTH", dome2d_depth_bay, &
                 default=0.2, do_not_log=.true.)
  call get_param(param_file,mod,"S_REF",S_ref,'Reference salinity',units='1e-3',fail_if_missing=.true.)
  call get_param(param_file,mod,"T_REF",T_ref,'Refernce temperature',units='C',fail_if_missing=.true.)
  call get_param(param_file,mod,"S_RANGE",S_range,'Initial salinity range',units='1e-3',default=2.0)
  call get_param(param_file,mod,"T_RANGE",T_range,'Initial temperature range',units='1e-3',default=0.0)

  T(:,:,:) = 0.0
  S(:,:,:) = 0.0

  ! Linear salinity profile

  select case ( coordinateMode(verticalCoordinate) )

    case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA )

      do j=js,je ; do i=is,ie
        xi0 = 0.0;
        do k = 1,nz
          xi1 = xi0 + h(i,j,k) / G%max_depth;
          S(i,j,k) = 34.0 + 0.5 * S_range * (xi0 + xi1);
          xi0 = xi1;
        enddo
      enddo ; enddo

    case ( REGRIDDING_RHO )

      do j=js,je ; do i=is,ie
        xi0 = 0.0;
        do k = 1,nz
          xi1 = xi0 + h(i,j,k) / G%max_depth;
          S(i,j,k) = 34.0 + 0.5 * S_range * (xi0 + xi1);
          xi0 = xi1;
        enddo
        x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon;
        if ( x .le. dome2d_width_bay ) then
          S(i,j,nz) = 34.0 + S_range;
        endif
      enddo ; enddo

    case ( REGRIDDING_LAYER )

      delta_S = S_range / ( G%ke - 1.0 );
      S(:,:,1) = S_ref;
      do k = 2,G%ke
        S(:,:,k) = S(:,:,k-1) + delta_S;
      enddo

    case default
      call MOM_error(FATAL,"dome2d_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

  ! Modify salinity and temperature when z coordinates are used
  if ( coordinateMode(verticalCoordinate) .eq. REGRIDDING_ZSTAR ) then
    index_bay_z = Nint ( dome2d_depth_bay * G%ke );
    do j = G%jsc,G%jec ; do i = G%isc,G%iec
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon;
      if ( x .le. dome2d_width_bay ) then
        S(i,j,1:index_bay_z) = S_ref + S_range; ! Use for z coordinates
        T(i,j,1:index_bay_z) = 1.0;             ! Use for z coordinates
      endif
    enddo ; enddo ! i and j loops
  endif ! Z initial conditions

  ! Modify salinity and temperature when sigma coordinates are used
  if ( coordinateMode(verticalCoordinate) .eq. REGRIDDING_SIGMA ) then
    do i = G%isc,G%iec ; do j = G%jsc,G%jec
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon;
      if ( x .le. dome2d_width_bay ) then
        S(i,j,1:G%ke) = S_ref + S_range;    ! Use for sigma coordinates
        T(i,j,1:G%ke) = 1.0;                ! Use for sigma coordinates
      endif
    enddo ; enddo
  endif

  ! Modify temperature when rho coordinates are used
  T(G%isc:G%iec,G%jsc:G%jec,1:G%ke) = 0.0;
  if (( coordinateMode(verticalCoordinate) .eq. REGRIDDING_RHO ) .or. ( coordinateMode(verticalCoordinate) .eq. REGRIDDING_LAYER )) then
    do i = G%isc,G%iec ; do j = G%jsc,G%jec
      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon;
      if ( x .le. dome2d_width_bay ) then
        T(i,j,G%ke) = 1.0;
      end if
    end do ; end do
  end if

end subroutine DOME2d_initialize_temperature_salinity

!> Set up sponges in 2d DOME configuration
subroutine DOME2d_initialize_sponges(G, GV, tv, param_file, use_ALE, CSp, ACSp)
  type(ocean_grid_type),   intent(in) :: G !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  type(thermo_var_ptrs),   intent(in) :: tv !< Thermodynamics structure
  type(param_file_type),   intent(in) :: param_file !< Parameter file structure
  logical,                 intent(in) :: use_ALE !< If true, indicates model is in ALE mode
  type(sponge_CS),         pointer    :: CSp !< Layer-mode sponge structure
  type(ALE_sponge_CS),     pointer    :: ACSp !< ALE-mode sponge structure
  ! Local variables
  real :: T(SZI_(G),SZJ_(G),SZK_(G))   ! A temporary array for temp
  real :: S(SZI_(G),SZJ_(G),SZK_(G))   ! A temporary array for salt
  real :: RHO(SZI_(G),SZJ_(G),SZK_(G)) ! A temporary array for RHO
  real :: h(SZI_(G),SZJ_(G),SZK_(G))   ! A temporary array for thickness
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for thickness
  real :: Idamp(SZI_(G),SZJ_(G))       ! The inverse damping rate, in s-1.
  real :: S_ref, T_ref                 ! Reference salinity and temerature within surface layer
  real :: S_range, T_range             ! Range of salinities and temperatures over the vertical
  real :: e0(SZK_(G)+1)             ! The resting interface heights, in m, usually !
                                    ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)          ! Interface height relative to the sea surface !
                                    ! positive upward, in m.
  real :: dome2d_width_bay, dome2d_width_bottom, dome2d_depth_bay
  real :: dome2d_west_sponge_time_scale, dome2d_east_sponge_time_scale
  real :: dome2d_west_sponge_width, dome2d_east_sponge_width
  real :: dummy1, x, z
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call get_param(param_file, mod, "DOME2D_WEST_SPONGE_TIME_SCALE", dome2d_west_sponge_time_scale, &
                 'The time-scale on the west edge of the domain for restoring T/S\n' //&
                 'in the sponge. If zero, the western sponge is disabled', &
                 units='s', default=0.)
  call get_param(param_file, mod, "DOME2D_EAST_SPONGE_TIME_SCALE", dome2d_east_sponge_time_scale, &
                 'The time-scale on the east edge of the domain for restoring T/S\n' //&
                 'in the sponge. If zero, the eastern sponge is disabled', &
                 units='s', default=0.)
  call get_param(param_file, mod, "DOME2D_WEST_SPONGE_WIDTH", dome2d_west_sponge_width, &
                 'The fraction of the domain in which the western sponge for restoring T/S\n' //&
                 'is active.', &
                 units='nondim', default=0.1)
  call get_param(param_file, mod, "DOME2D_EAST_SPONGE_WIDTH", dome2d_east_sponge_width, &
                 'The fraction of the domain in which the eastern sponge for restoring T/S\n' //&
                 'is active.', &
                 units='nondim', default=0.1)

  ! Return if sponges are not in use
  if (dome2d_west_sponge_time_scale <= 0. .and. dome2d_east_sponge_time_scale <= 0.) return

  if (associated(CSp)) call MOM_error(FATAL, &
     "DOME2d_initialize_sponges called with an associated control structure.")
  if (associated(ACSp)) call MOM_error(FATAL, &
     "DOME2d_initialize_sponges called with an associated ALE-sponge control structure.")

  call get_param(param_file, mod, "DOME2D_SHELF_WIDTH", dome2d_width_bay, &
                 default=0.1, do_not_log=.true.)
  call get_param(param_file, mod, "DOME2D_BASIN_WIDTH", dome2d_width_bottom, &
                 default=0.3, do_not_log=.true.)
  call get_param(param_file, mod, "DOME2D_SHELF_DEPTH", dome2d_depth_bay, &
                 default=0.2, do_not_log=.true.)
  call get_param(param_file,mod,"S_REF",S_ref)
  call get_param(param_file,mod,"T_REF",T_ref)
  call get_param(param_file,mod,"S_RANGE",S_range,default=2.0)
  call get_param(param_file,mod,"T_RANGE",T_range,default=0.0)


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
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) = e0(k)
        if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_z)) then
          eta1D(k) = eta1D(k+1) + GV%Angstrom_z
          h(i,j,k) = GV%Angstrom_z
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo
    enddo;  enddo
    ! Store the grid on which the T/S sponge data will reside
    call initialize_ALE_sponge(Idamp, h, nz, G, param_file, ACSp)

    ! Construct temperature and salinity on the arbitrary grid
    T(:,:,:) = 0.0 ; S(:,:,:) = 0.0
    do j=js,je ; do i=is,ie
      z = -G%bathyT(i,j)
      do k = nz,1,-1
        z = z + 0.5 * h(i,j,k) ! Position of the center of layer k
        S(i,j,k) = 34.0 - 1.0 * (z/G%max_depth)
        if ( ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon < dome2d_west_sponge_width ) S(i,j,k) = S_ref + S_range
        z = z + 0.5 * h(i,j,k) ! Position of the interface k
      enddo
    enddo ; enddo

    if ( associated(tv%T) ) then
      call set_up_ALE_sponge_field(T,G,tv%T,ACSp)
    endif
    if ( associated(tv%S) ) then
      call set_up_ALE_sponge_field(S,G,tv%S,ACSp)
    endif

  else

    ! Construct thicknesses to restore to
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) = -G%max_depth * real(k-1) / real(nz)
        if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_z)) then
          eta1D(k) = eta1D(k+1) + GV%Angstrom_z
          h(i,j,k) = GV%Angstrom_z
        else
          h(i,j,k) = eta1D(k) - eta1D(k+1)
        endif
      enddo

      x = ( G%geoLonT(i,j) - G%west_lon ) / G%len_lon;
      if ( x .le. dome2d_width_bay ) then
        h(i,j,1:nz-1) = GV%Angstrom;
        h(i,j,nz) = dome2d_depth_bay * G%max_depth - (nz-1) * GV%Angstrom;
      end if

      eta(i,j,nz+1) = -G%bathyT(i,j)
      do K=nz,1,-1
        eta(i,j,K) = eta(i,j,K+1) + h(i,j,k)
      enddo
    enddo ; enddo
    call initialize_sponge(Idamp, eta, G, param_file, CSp)

  endif

end subroutine DOME2d_initialize_sponges

end module DOME2d_initialization
