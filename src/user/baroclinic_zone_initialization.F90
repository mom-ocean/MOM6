!> Initial conditions for an idealized baroclinic zone
module baroclinic_zone_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_file_parser,   only : openParameterBlock, closeParameterBlock
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>
#include "version_variable.h"

! Private (module-wise) parameters
character(len=40) :: mdl = "baroclinic_zone_initialization" !< This module's name.

public baroclinic_zone_init_temperature_salinity

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Reads the parameters unique to this module
subroutine bcz_params(G, GV, US, param_file, S_ref, dSdz, delta_S, dSdx, T_ref, dTdz, &
                      delta_T, dTdx, L_zone, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G          !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US    !< A dimensional unit scaling type
  type(param_file_type),   intent(in)  :: param_file !< Parameter file handle
  real,                    intent(out) :: S_ref      !< Reference salinity [ppt]
  real,                    intent(out) :: dSdz       !< Salinity stratification [ppt Z-1 ~> ppt m-1]
  real,                    intent(out) :: delta_S    !< Salinity difference across baroclinic zone [ppt]
  real,                    intent(out) :: dSdx       !< Linear salinity gradient [ppt m-1]
  real,                    intent(out) :: T_ref      !< Reference temperature [degC]
  real,                    intent(out) :: dTdz       !< Temperature stratification [degC Z-1 ~> degC m-1]
  real,                    intent(out) :: delta_T    !< Temperature difference across baroclinic zone [degC]
  real,                    intent(out) :: dTdx       !< Linear temperature gradient in [degC G%x_axis_units-1]
  real,                    intent(out) :: L_zone     !< Width of baroclinic zone in [G%x_axis_units]
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                     !! only read parameters without changing h.

  logical :: just_read    ! If true, just read parameters but set nothing.

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) &
    call log_version(param_file, mdl, version, 'Initialization of an analytic baroclinic zone')
  call openParameterBlock(param_file,'BCZIC')
  call get_param(param_file, mdl, "S_REF", S_ref, 'Reference salinity', units='ppt', &
                 default=35., do_not_log=just_read)
  call get_param(param_file, mdl, "DSDZ", dSdz, 'Salinity stratification', &
                 units='ppt/m', default=0.0, scale=US%Z_to_m, do_not_log=just_read)
  call get_param(param_file, mdl,"DELTA_S",delta_S,'Salinity difference across baroclinic zone', &
                 units='ppt', default=0.0, do_not_log=just_read)
  call get_param(param_file, mdl,"DSDX",dSdx,'Meridional salinity difference', &
                 units='ppt/'//trim(G%x_axis_units), default=0.0, do_not_log=just_read)
  call get_param(param_file, mdl,"T_REF",T_ref,'Reference temperature',units='C', &
                 default=10., do_not_log=just_read)
  call get_param(param_file, mdl, "DTDZ", dTdz, 'Temperature stratification', &
                 units='C/m', default=0.0, scale=US%Z_to_m, do_not_log=just_read)
  call get_param(param_file, mdl,"DELTA_T",delta_T,'Temperature difference across baroclinic zone', &
                 units='C', default=0.0, do_not_log=just_read)
  call get_param(param_file, mdl,"DTDX",dTdx,'Meridional temperature difference', &
                 units='C/'//trim(G%x_axis_units), default=0.0, do_not_log=just_read)
  call get_param(param_file, mdl,"L_ZONE",L_zone,'Width of baroclinic zone', &
                 units=G%x_axis_units, default=0.5*G%len_lat, do_not_log=just_read)
  call closeParameterBlock(param_file)

end subroutine bcz_params

!> Initialization of temperature and salinity with the baroclinic zone initial conditions
subroutine baroclinic_zone_init_temperature_salinity(T, S, h, depth_tot, G, GV, US, param_file, &
                                                     just_read_params)
  type(ocean_grid_type),                     intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< Potential temperature [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< Salinity [ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h  !< The model thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),                     intent(in)  :: param_file  !< Parameter file handle
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing T & S.

  integer   :: i, j, k, is, ie, js, je, nz
  real      :: T_ref, dTdz, dTdx, delta_T ! Parameters describing temperature distribution
  real      :: S_ref, dSdz, dSdx, delta_S ! Parameters describing salinity distribution
  real      :: L_zone ! Width of baroclinic zone in [G%axis_units]
  real      :: zc, zi ! Depths in depth units [Z ~> m]
  real      :: x, xd, xs, y, yd, fn
  real      :: PI                   ! 3.1415926... calculated as 4*atan(1)
  logical :: just_read    ! If true, just read parameters but set nothing.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call bcz_params(G, GV, US, param_file, S_ref, dSdz, delta_S, dSdx, T_ref, dTdz, &
                  delta_T, dTdx, L_zone, just_read_params)

  if (just_read) return ! All run-time parameters have been read, so return.

  T(:,:,:) = 0.
  S(:,:,:) = 0.
  PI = 4.*atan(1.)

  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    zi = -depth_tot(i,j)
    x = G%geoLonT(i,j) - (G%west_lon + 0.5*G%len_lon) ! Relative to center of domain
    xd = x / G%len_lon ! -1/2 < xd 1/2
    y = G%geoLatT(i,j) - (G%south_lat + 0.5*G%len_lat) ! Relative to center of domain
    yd = y / G%len_lat ! -1/2 < yd 1/2
    if (L_zone/=0.) then
      xs = min(1., max(-1., x/L_zone)) ! -1 < ys < 1
      fn = sin((0.5*PI)*xs)
    else
      xs = sign(1., x) ! +/- 1
      fn = xs
    endif
    do k = nz, 1, -1
      zc = zi + 0.5*h(i,j,k)*GV%H_to_Z ! Position of middle of cell
      zi = zi + h(i,j,k)*GV%H_to_Z    ! Top interface position
      T(i,j,k) = T_ref + dTdz * zc  & ! Linear temperature stratification
                 + dTdx * x         & ! Linear gradient
                 + delta_T * fn       ! Smooth fn of width L_zone
      S(i,j,k) = S_ref + dSdz * zc  & ! Linear temperature stratification
                 + dSdx * x         & ! Linear gradient
                 + delta_S * fn       ! Smooth fn of width L_zone
    enddo
  enddo ; enddo

end subroutine baroclinic_zone_init_temperature_salinity

!> \namespace baroclinic_zone_initialization
!!
!! \section section_baroclinic_zone Description of the baroclinic zone initial conditions
!!
!! yada yada yada

end module baroclinic_zone_initialization
