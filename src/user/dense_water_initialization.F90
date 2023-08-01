!> Initialization routines for the dense water formation
!! and overflow experiment.
module dense_water_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_ALE_sponge,    only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_dyn_horgrid,   only : dyn_horgrid_type
use MOM_EOS,           only : EOS_type
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser,   only : get_param, param_file_type
use MOM_interface_heights, only : dz_to_thickness, dz_to_thickness_simple
use MOM_grid,          only : ocean_grid_type
use MOM_sponge,        only : sponge_CS
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public dense_water_initialize_topography
public dense_water_initialize_TS
public dense_water_initialize_sponges

character(len=40) :: mdl = "dense_water_initialization" !< Module name

real, parameter :: default_sill  = 0.2  !< Default depth of the sill [nondim]
real, parameter :: default_shelf = 0.4  !< Default depth of the shelf [nondim]
real, parameter :: default_mld   = 0.25 !< Default depth of the mixed layer [nondim]

contains

!> Initialize the topography field for the dense water experiment
subroutine dense_water_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth [Z ~> m]

  ! Local variables
  real, dimension(5) :: domain_params ! nondimensional widths of all domain sections [nondim]
  real :: sill_frac     ! Depth of the sill separating downslope from upslope, as a fraction of
                        ! the basin depth [nondim]
  real :: shelf_frac    ! Depth of the shelf region accumulating dense water for overflow,
                        ! as a fraction the basin depth [nondim]
  real :: x             ! Horizontal position normalized by the domain width [nondim]
  integer :: i, j

  call get_param(param_file, mdl, "DENSE_WATER_DOMAIN_PARAMS", domain_params, &
       "Fractional widths of all the domain sections for the dense water experiment.\n"//&
       "As a 5-element vector:\n"//&
       "  - open ocean, the section at maximum depth\n"//&
       "  - downslope, the downward overflow slope\n"//&
       "  - sill separating downslope from upslope\n"//&
       "  - upslope, the upward slope accumulating dense water\n"//&
       "  - the shelf in the dense formation region.", &
       units="nondim", fail_if_missing=.true.)
  call get_param(param_file, mdl, "DENSE_WATER_SILL_DEPTH", sill_frac, &
       "Depth of the sill separating downslope from upslope, as fraction of basin depth.", &
       units="nondim", default=default_sill)
  call get_param(param_file, mdl, "DENSE_WATER_SHELF_DEPTH", shelf_frac, &
       "Depth of the shelf region accumulating dense water for overflow, as fraction of basin depth.", &
       units="nondim", default=default_shelf)

  do i = 2, 5
    ! turn widths into positions
    domain_params(i) = domain_params(i-1) + domain_params(i)
  enddo

  do j = G%jsc,G%jec
    do i = G%isc,G%iec
      ! compute normalised zonal coordinate
      x = (G%geoLonT(i,j) - G%west_lon) / G%len_lon

      if (x <= domain_params(1)) then
        ! open ocean region
        D(i,j) = max_depth
      elseif (x <= domain_params(2)) then
        ! downslope region, linear
        D(i,j) = max_depth - (1.0 - sill_frac) * max_depth * &
             (x - domain_params(1)) / (domain_params(2) - domain_params(1))
      elseif (x <= domain_params(3)) then
        ! sill region
        D(i,j) = sill_frac * max_depth
      elseif (x <= domain_params(4)) then
        ! upslope region
        D(i,j) = sill_frac * max_depth + (shelf_frac - sill_frac) * max_depth * &
             (x - domain_params(3)) / (domain_params(4) - domain_params(3))
      else
        ! shelf region
        D(i,j) = shelf_frac * max_depth
      endif
    enddo
  enddo

end subroutine dense_water_initialize_topography

!> Initialize the temperature and salinity for the dense water experiment
subroutine dense_water_initialize_TS(G, GV, US, param_file, T, S, h, just_read)
  type(ocean_grid_type),                     intent(in)  :: G !< Horizontal grid control structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid control structure
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T !< Output temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S !< Output salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h !< Layer thicknesses [Z ~> m]
  logical,                                   intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing T & S.
  ! Local variables
  real :: mld             ! The initial mixed layer depth as a fraction of the maximum depth [nondim]
  real :: S_ref, S_range  ! The reference salinity and its range in the initial conditions [S ~> ppt]
  real :: T_ref           ! The reference temperature [C ~> degC]
  real :: zi, zmid        ! Depths from the surface nondimensionalized by the maximum depth [nondim]
  integer :: i, j, k, nz

  nz = GV%ke

  call get_param(param_file, mdl, "DENSE_WATER_MLD", mld, &
       "Depth of unstratified mixed layer as a fraction of the water column.", &
       units="nondim", default=default_mld, do_not_log=just_read)
  call get_param(param_file, mdl, "S_REF", S_ref, 'Reference salinity', &
                 default=35.0, units='1e-3', scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl,"T_REF", T_ref, 'Reference temperature', &
                units='degC', scale=US%degC_to_C, fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl,"S_RANGE", S_range, 'Initial salinity range', &
                units='1e-3', default=2.0, scale=US%ppt_to_S, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  ! uniform temperature everywhere
  T(:,:,:) = T_ref

  do j = G%jsc,G%jec
    do i = G%isc,G%iec
      zi = 0.
      do k = 1,nz
        ! nondimensional middle of layer
        zmid = zi + 0.5 * h(i,j,k) / G%max_depth

        if (zmid < mld) then
          ! use reference salinity in the mixed layer
          S(i,j,k) = S_ref
        else
          ! linear between bottom of mixed layer and bottom
          S(i,j,k) = S_ref + S_range * (zmid - mld) / (1.0 - mld)
        endif

        zi = zi + h(i,j,k) / G%max_depth
      enddo
    enddo
  enddo
end subroutine dense_water_initialize_TS

!> Initialize the restoring sponges for the dense water experiment
subroutine dense_water_initialize_sponges(G, GV, US, tv, depth_tot, param_file, use_ALE, CSp, ACSp)
  type(ocean_grid_type),   intent(in) :: G !< Horizontal grid control structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid control structure
  type(unit_scale_type),   intent(in) :: US !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv !< Thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: param_file !< Parameter file structure
  logical,                 intent(in) :: use_ALE !< ALE flag
  type(sponge_CS),         pointer    :: CSp !< Layered sponge control structure pointer
  type(ALE_sponge_CS),     pointer    :: ACSp !< ALE sponge control structure pointer

  ! Local variables
  real :: west_sponge_time_scale, east_sponge_time_scale ! Sponge timescales [T ~> s]
  real :: west_sponge_width ! The fraction of the domain in which the western (outflow) sponge is active [nondim]
  real :: east_sponge_width ! The fraction of the domain in which the eastern (outflow) sponge is active [nondim]

  real, dimension(SZI_(G),SZJ_(G)) :: Idamp ! inverse damping timescale [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: dz ! sponge layer thicknesses in height units [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h  ! sponge layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: T  ! sponge temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: S  ! sponge salinity [S ~> ppt]
  real, dimension(SZK_(GV)+1) :: e0, eta1D ! interface positions for ALE sponge [Z ~> m]
  real :: x         ! Horizontal position normalized by the domain width [nondim]
  real :: zi, zmid  ! Depths from the surface nondimensionalized by the maximum depth [nondim]
  real :: dist      ! Distance from the edge of a sponge normalized by the width of that sponge [nondim]
  real :: mld       ! The initial mixed layer depth as a fraction of the maximum depth [nondim]
  real :: S_ref, S_range  ! The reference salinity and its range in the initial conditions [S ~> ppt]
  real :: S_dense   ! The salinity of the dense water being formed on the shelf [S ~> ppt]
  real :: T_ref     ! The reference temperature [C ~> degC]
  real :: sill_frac ! Fractional depths of the sill, relative to the maximum depth [nondim]
  integer :: i, j, k, nz

  nz = GV%ke

  call get_param(param_file, mdl, "DENSE_WATER_WEST_SPONGE_TIME_SCALE", west_sponge_time_scale, &
                 "The time scale on the west (outflow) of the domain for restoring. "//&
                 "If zero, the sponge is disabled.", units="s", default=0., scale=US%s_to_T)
  call get_param(param_file, mdl, "DENSE_WATER_WEST_SPONGE_WIDTH", west_sponge_width, &
                 "The fraction of the domain in which the western (outflow) sponge is active.", &
                 units="nondim", default=0.1)
  call get_param(param_file, mdl, "DENSE_WATER_EAST_SPONGE_TIME_SCALE", east_sponge_time_scale, &
                 "The time scale on the east (outflow) of the domain for restoring. "//&
                 "If zero, the sponge is disabled.", units="s", default=0., scale=US%s_to_T)
  call get_param(param_file, mdl, "DENSE_WATER_EAST_SPONGE_WIDTH", east_sponge_width, &
                 "The fraction of the domain in which the eastern (outflow) sponge is active.", &
                 units="nondim", default=0.1)
  call get_param(param_file, mdl, "DENSE_WATER_EAST_SPONGE_SALT", S_dense, &
                 "Salt anomaly of the dense water being formed in the overflow region.", &
                 units="1e-3", default=4.0, scale=US%ppt_to_S)

  call get_param(param_file, mdl, "DENSE_WATER_MLD", mld, &
                 units="nondim", default=default_mld, do_not_log=.true.)
  call get_param(param_file, mdl, "DENSE_WATER_SILL_DEPTH", sill_frac, &
                 units="nondim", default=default_sill, do_not_log=.true.)

  call get_param(param_file, mdl, "S_REF", S_ref, &
                 units='1e-3', default=35.0, scale=US%ppt_to_S, do_not_log=.true.)
  call get_param(param_file, mdl, "S_RANGE", S_range, &
                 units='1e-3', default=2.0, scale=US%ppt_to_S, do_not_log=.true.)
  call get_param(param_file, mdl, "T_REF", T_ref, &
                 units='degC', scale=US%degC_to_C, fail_if_missing=.true., do_not_log=.true.)

  ! no active sponges
  if (west_sponge_time_scale <= 0. .and. east_sponge_time_scale <= 0.) return

  ! everywhere is initially unsponged
  Idamp(:,:) = 0.0

  do j = G%jsc, G%jec
    do i = G%isc,G%iec
      if (G%mask2dT(i,j) > 0.) then
        ! nondimensional x position
        x = (G%geoLonT(i,j) - G%west_lon) / G%len_lon

        if (west_sponge_time_scale > 0. .and. x < west_sponge_width) then
          dist = 1. - x / west_sponge_width
          ! scale restoring by depth into sponge
          Idamp(i,j) = 1. / west_sponge_time_scale * max(0., min(1., dist))
        elseif (east_sponge_time_scale > 0. .and. x > (1. - east_sponge_width)) then
          dist = 1. - (1. - x) / east_sponge_width
          Idamp(i,j) = 1. / east_sponge_time_scale * max(0., min(1., dist))
        endif
      endif
    enddo
  enddo

  if (use_ALE) then
    ! construct a uniform grid for the sponge
    do k = 1,nz
      e0(k) = -G%max_depth * (real(k - 1) / real(nz))
    enddo
    e0(nz+1) = -G%max_depth

    do j = G%jsc,G%jec
      do i = G%isc,G%iec
        eta1D(nz+1) = -depth_tot(i,j)
        do k = nz,1,-1
          eta1D(k) = e0(k)

          if (eta1D(k) < (eta1D(k+1) + GV%Angstrom_Z)) then
            ! is this layer vanished?
            eta1D(k) = eta1D(k+1) + GV%Angstrom_Z
            dz(i,j,k) = GV%Angstrom_Z
          else
            dz(i,j,k) = eta1D(k) - eta1D(k+1)
          endif
        enddo
      enddo
    enddo

    ! construct temperature and salinity for the sponge
    ! start with initial condition
    T(:,:,:) = T_ref
    S(:,:,:) = S_ref

    do j = G%jsc,G%jec
      do i = G%isc,G%iec
        zi = 0.
        x = (G%geoLonT(i,j) - G%west_lon) / G%len_lon
        do k = 1,nz
          ! nondimensional middle of layer
          zmid = zi + 0.5 * dz(i,j,k) / G%max_depth

          if (x > (1. - east_sponge_width)) then
            !if (zmid >= 0.9 * sill_frac) &
              S(i,j,k) = S_ref + S_dense
          else
            ! linear between bottom of mixed layer and bottom
            if (zmid >= mld) &
              S(i,j,k) = S_ref + S_range * (zmid - mld) / (1.0 - mld)
          endif

          zi = zi + dz(i,j,k) / G%max_depth
        enddo
      enddo
    enddo

    ! Convert thicknesses from height units to thickness units
    if (associated(tv%eqn_of_state)) then
      call dz_to_thickness(dz, T, S, tv%eqn_of_state, h, G, GV, US)
    else
      call dz_to_thickness_simple(dz, h, G, GV, US, layer_mode=.true.)
    endif

    ! This call sets up the damping rates and interface heights in the sponges.
    call initialize_ALE_sponge(Idamp, G, GV, param_file, ACSp, h, nz)

    if ( associated(tv%T) ) call set_up_ALE_sponge_field(T, G, GV, tv%T, ACSp, 'temp', &
        sp_long_name='temperature', sp_unit='degC s-1')
    if ( associated(tv%S) ) call set_up_ALE_sponge_field(S, G, GV, tv%S, ACSp, 'salt', &
        sp_long_name='salinity', sp_unit='g kg-1 s-1')
  else
    call MOM_error(FATAL, "dense_water_initialize_sponges: trying to use non ALE sponge")
  endif
end subroutine dense_water_initialize_sponges

end module dense_water_initialization

!> \namespace dense_water_initialization
!!
!! This experiment consists of a shelf accumulating dense water, which spills
!! over an upward slope and a sill, before flowing down a slope into an open
!! ocean region. It's intended as a test of one of the motivating situations
!! for the adaptive coordinate.
!!
!! The nondimensional widths of the 5 regions are controlled by the
!! <code>DENSE_WATER_DOMAIN_PARAMS</code>, and the heights of the sill and shelf
!! as a fraction of the total domain depth are controlled by
!! <code>DENSE_WATER_SILL_DEPTH</code> and <code>DENSE_WATER_SHELF_DEPTH</code>.
!!
!! The density in the domain is governed by a linear equation of state, and
!! is set up with a mixed layer of non-dimensional depth <code>DENSE_WATER_MLD</code>
!! below which there is a linear stratification from <code>S_REF</code>, increasing
!! by <code>S_RANGE</code> to the bottom.
!!
!! To force the experiment, there are sponges on the inflow and outflow of the
!! domain. The inflow sponge has a salinity anomaly of
!! <code>DENSE_WATER_EAST_SPONGE_SALT</code> through the entire depth. The outflow
!! sponge simply restores to the initial condition. Both sponges have controllable
!! widths and restoring timescales.
