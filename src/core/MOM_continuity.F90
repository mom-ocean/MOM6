!> Solve the layer continuity equation.
module MOM_continuity

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_continuity_PPM, only : continuity_PPM, zonal_mass_flux, meridional_mass_flux
use MOM_continuity_PPM, only : continuity_PPM_stencil, continuity_PPM_init
use MOM_continuity_PPM, only : continuity_PPM_CS
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_string_functions, only : uppercase
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : BT_cont_type, porous_barrier_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public continuity, continuity_fluxes, continuity_adjust_vel, continuity_init, continuity_stencil

!> Control structure for mom_continuity
type, public :: continuity_CS ; private
  integer :: continuity_scheme !< Selects the discretization for the continuity solver.
                               !! Valid values are:
                               !! - PPM - A directionally split piecewise parabolic reconstruction solver.
                               !! The default, PPM, seems most appropriate for use with our current
                               !! time-splitting strategies.
  type(continuity_PPM_CS) :: PPM  !< Control structure for mom_continuity_ppm
end type continuity_CS

integer, parameter :: PPM_SCHEME = 1 !< Enumerated constant to select PPM
character(len=20), parameter :: PPM_STRING = "PPM" !< String to select PPM

!> Finds the thickness fluxes from the continuity solver or their vertical sum without
!! actually updating the layer thicknesses.
interface continuity_fluxes
  module procedure continuity_3d_fluxes, continuity_2d_fluxes
end interface continuity_fluxes

contains

!> Time steps the layer thicknesses, using a monotonically or positive-definite limited, directionally
!! split PPM scheme, based on Lin (1994).
subroutine continuity(u, v, hin, h, uh, vh, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, &
                      visc_rem_u, visc_rem_v, u_cor, v_cor, BT_cont)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: hin !< Initial layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(inout) :: h   !< Final layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: uh  !< Volume flux through zonal faces =
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out)   :: vh  !< Volume flux through meridional faces =
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_CS),     intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(in)    :: uhbt !< The vertically summed volume flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(in)    :: vhbt !< The vertically summed volume flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u !< Both the fraction of
          !! zonal momentum that remains after a time-step of viscosity, and the fraction of a time-step's
          !! worth of a barotropic acceleration that a layer experiences after viscosity is applied [nondim].
          !! Non-dimensional between 0 (at the bottom) and 1 (far above the bottom).  When this column is
          !! under an ice shelf, this can also go to 0 at the top due to the no-slip boundary condition there.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_v !< Both the fraction of
          !! meridional momentum that remains after a time-step of viscosity, and the fraction of a time-step's
          !! worth of a barotropic acceleration that a layer experiences after viscosity is applied [nondim].
          !! Non-dimensional between 0 (at the bottom) and 1 (far above the bottom).  When this column is
          !! under an ice shelf, this can also go to 0 at the top due to the no-slip boundary condition there.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: u_cor !< The zonal velocities that
          !! give uhbt as the depth-integrated transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(out)   :: v_cor !< The meridional velocities that
          !! give vhbt as the depth-integrated transport [L T-1 ~> m s-1].
  type(BT_cont_type), &
                 optional, pointer       :: BT_cont !< A structure with elements
          !! that describe the effective open face areas as a function of barotropic flow.

  if (present(visc_rem_u) .neqv. present(visc_rem_v)) call MOM_error(FATAL, &
      "MOM_continuity: Either both visc_rem_u and visc_rem_v or neither"// &
       " one must be present in call to continuity.")
  if (present(u_cor) .neqv. present(v_cor)) call MOM_error(FATAL, &
      "MOM_continuity: Either both u_cor and v_cor or neither"// &
       " one must be present in call to continuity.")

  if (CS%continuity_scheme == PPM_SCHEME) then
    call continuity_PPM(u, v, hin, h, uh, vh, dt, G, GV, US, CS%PPM, OBC, pbv, uhbt, vhbt, &
                        visc_rem_u, visc_rem_v, u_cor, v_cor, BT_cont=BT_cont)
  else
    call MOM_error(FATAL, "continuity: Unrecognized value of continuity_scheme")
  endif

end subroutine continuity

!> Finds the thickness fluxes from the continuity solver without actually updating the
!! layer thicknesses.  Because the fluxes in the two directions are calculated based on the
!! input thicknesses, which are not updated between the direcitons, the fluxes returned here
!! are not the same as those that would be returned by a call to continuity.
subroutine continuity_3d_fluxes(u, v, h, uh, vh, dt, G, GV, US, CS, OBC, pbv)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: uh  !< Thickness fluxes through zonal faces,
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out)   :: vh  !< Thickness fluxes through meridional faces,
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_CS),     intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  if (CS%continuity_scheme == PPM_SCHEME) then
    call zonal_mass_flux(u, h, uh, dt, G, GV, US, CS%PPM, OBC, pbv%por_face_areaU)

    call meridional_mass_flux(v, h, vh, dt, G, GV, US, CS%PPM, OBC, pbv%por_face_areaV)
  else
    call MOM_error(FATAL, "continuity: Unrecognized value of continuity_scheme")
  endif

end subroutine continuity_3d_fluxes

!> Find the vertical sum of the thickness fluxes from the continuity solver without actually
!! updating the layer thicknesses.  Because the fluxes in the two directions are calculated
!! based on the input thicknesses, which are not updated between the directions, the fluxes
!! returned here are not the same as those that would be returned by a call to continuity.
subroutine continuity_2d_fluxes(u, v, h, uhbt, vhbt, dt, G, GV, US, CS, OBC, pbv)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(out)   :: uhbt !< Vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(out)   :: vhbt !< Vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_CS),     intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uh  ! Thickness fluxes through zonal faces,
                                                    ! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vh  ! Thickness fluxes through meridional faces,
                                                    ! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  integer :: i, j, k

  uh(:,:,:) = 0.0
  vh(:,:,:) = 0.0

  if (CS%continuity_scheme == PPM_SCHEME) then
    call zonal_mass_flux(u, h, uh, dt, G, GV, US, CS%PPM, OBC, pbv%por_face_areaU)

    call meridional_mass_flux(v, h, vh, dt, G, GV, US, CS%PPM, OBC, pbv%por_face_areaV)
  else
    call MOM_error(FATAL, "continuity: Unrecognized value of continuity_scheme")
  endif

  uhbt(:,:) = 0.0
  vhbt(:,:) = 0.0

  do k=1,GV%ke ; do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
    uhbt(I,j) = uhbt(I,j) + uh(I,j,k)
  enddo ; enddo ; enddo

  do k=1,GV%ke ; do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
    vhbt(I,j) = vhbt(I,j) + vh(I,j,k)
  enddo ; enddo ; enddo

end subroutine continuity_2d_fluxes


!> Correct the velocities to give the specified depth-integrated transports by applying a
!! barotropic acceleration (subject to viscous drag) to the velocities.
subroutine continuity_adjust_vel(u, v, h, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, visc_rem_u, visc_rem_v)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u   !< Zonal velocity, which will be adjusted to
                                                !! give uhbt as the depth-integrated
                                                !! transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v   !< Meridional velocity, which will be adjusted
                                                !! to give vhbt as the depth-integrated
                                                !! transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_CS),     intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: uhbt !< The vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in)    :: vhbt !< The vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u !< Both the fraction of the zonal momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_v !< Both the fraction of the meridional momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u_in  !< Input zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: v_in  !< Input meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uh  !< Volume flux through zonal faces =
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vh  !< Volume flux through meridional faces =
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].

  ! It might not be necessary to separate the input velocity array from the adjusted velocities,
  ! but it seems safer to do so, even if it might be less efficient.
  u_in(:,:,:) = u(:,:,:)
  v_in(:,:,:) = v(:,:,:)

  if (CS%continuity_scheme == PPM_SCHEME) then
    call zonal_mass_flux(u_in, h, uh, dt, G, GV, US, CS%PPM, OBC, pbv%por_face_areaU, &
                         uhbt=uhbt, visc_rem_u=visc_rem_u, u_cor=u)

    call meridional_mass_flux(v_in, h, vh, dt, G, GV, US, CS%PPM, OBC, pbv%por_face_areaV, &
                              vhbt=vhbt, visc_rem_v=visc_rem_v, v_cor=v)
  else
    call MOM_error(FATAL, "continuity: Unrecognized value of continuity_scheme")
  endif

end subroutine continuity_adjust_vel

!> Initializes continuity_cs
subroutine continuity_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time       !< Current model time.
  type(ocean_grid_type),   intent(in)    :: G          !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handles.
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics control structure.
  type(continuity_CS),     intent(inout) :: CS         !< Control structure for mom_continuity.

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_continuity" ! This module's name.
  character(len=20)  :: tmpstr

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "CONTINUITY_SCHEME", tmpstr, &
                 "CONTINUITY_SCHEME selects the discretization for the "//&
                 "continuity solver. The only valid value currently is: \n"//&
                 "\t PPM - use a positive-definite (or monotonic) \n"//&
                 "\t       piecewise parabolic reconstruction solver.", &
                 default=PPM_STRING)

  tmpstr = uppercase(tmpstr) ; CS%continuity_scheme = 0
  select case (trim(tmpstr))
    case (PPM_STRING) ; CS%continuity_scheme = PPM_SCHEME
    case default
      call MOM_mesg('continuity_init: CONTINUITY_SCHEME ="'//trim(tmpstr)//'"', 0)
      call MOM_mesg("continuity_init: The only valid value is currently "// &
                     trim(PPM_STRING), 0)
      call MOM_error(FATAL, "continuity_init: Unrecognized setting "// &
            "#define CONTINUITY_SCHEME "//trim(tmpstr)//" found in input file.")
  end select

  if (CS%continuity_scheme == PPM_SCHEME) then
    call continuity_PPM_init(Time, G, GV, US, param_file, diag, CS%PPM)
  endif

end subroutine continuity_init


!> continuity_stencil returns the continuity solver stencil size
function continuity_stencil(CS) result(stencil)
  type(continuity_CS), intent(in) :: CS !< Module's control structure.
  integer ::  stencil !< The continuity solver stencil size with the current settings.

  stencil = 1

  if (CS%continuity_scheme == PPM_SCHEME) then
    stencil = continuity_PPM_stencil(CS%PPM)
  endif
end function continuity_stencil

end module MOM_continuity
