!> Solve the layer continuity equation.
module MOM_continuity

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_continuity_PPM, only : continuity_PPM, continuity_PPM_init
use MOM_continuity_PPM, only : continuity_PPM_stencil
use MOM_continuity_PPM, only : continuity_PPM_end, continuity_PPM_CS
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_string_functions, only : uppercase
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : BT_cont_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public continuity, continuity_init, continuity_end, continuity_stencil

!> Control structure for mom_continuity
type, public :: continuity_CS ; private
  integer :: continuity_scheme !< Selects the discretization for the continuity solver.
                               !! Valid values are:
                               !! - PPM - A directionally split piecewise parabolic reconstruction solver.
                               !! The default, PPM, seems most appropriate for use with our current
                               !! time-splitting strategies.
  type(continuity_PPM_CS), pointer :: PPM_CSp => NULL() !< Control structure for mom_continuity_ppm
end type continuity_CS

integer, parameter :: PPM_SCHEME = 1 !< Enumerated constant to select PPM
character(len=20), parameter :: PPM_STRING = "PPM" !< String to select PPM

contains

!> Time steps the layer thicknesses, using a monotonically limited, directionally split PPM scheme,
!! based on Lin (1994).
subroutine continuity(u, v, hin, h, uh, vh, dt, G, GV, US, CS, uhbt, vhbt, OBC, &
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
  type(continuity_CS),     pointer       :: CS  !< Control structure for mom_continuity.
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(in)    :: uhbt !< The vertically summed volume
                                                !! flux through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(in)    :: vhbt !< The vertically summed volume
                                                !! flux through meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(ocean_OBC_type), &
                 optional, pointer       :: OBC !< Open boundaries control structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u !< Both the fraction of
          !! zonal momentum that remains after a time-step of viscosity, and the fraction of a time-step's
          !! worth of a barotropic acceleration that a layer experiences after viscosity is applied.
          !! Non-dimensional between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_v !< Both the fraction of
          !! meridional momentum that remains after a time-step of viscosity, and the fraction of a time-step's
          !! worth of a barotropic acceleration that a layer experiences after viscosity is applied.
          !! Non-dimensional between 0 (at the bottom) and 1 (far above the bottom).
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
    call continuity_PPM(u, v, hin, h, uh, vh, dt, G, GV, US, CS%PPM_CSp, uhbt, vhbt, OBC, &
                        visc_rem_u, visc_rem_v, u_cor, v_cor, BT_cont=BT_cont)
  else
    call MOM_error(FATAL, "continuity: Unrecognized value of continuity_scheme")
  endif

end subroutine continuity

!> Initializes continuity_cs
subroutine continuity_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time       !< Current model time.
  type(ocean_grid_type),   intent(in)    :: G          !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handles.
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics control structure.
  type(continuity_CS),     pointer       :: CS         !< Control structure for mom_continuity.

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_continuity" ! This module's name.
  character(len=20)  :: tmpstr

  if (associated(CS)) then
    call MOM_error(WARNING, "continuity_init called with associated control structure.")
    return
  endif
  allocate(CS)

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
    call continuity_PPM_init(Time, G, GV, US, param_file, diag, CS%PPM_CSp)
  endif

end subroutine continuity_init


!> continuity_stencil returns the continuity solver stencil size
function continuity_stencil(CS) result(stencil)
  type(continuity_CS), pointer       :: CS  !< Module's control structure.
  integer ::  stencil !< The continuity solver stencil size with the current settings.

  stencil = 1

  if (CS%continuity_scheme == PPM_SCHEME) then
    stencil = continuity_PPM_stencil(CS%PPM_CSp)
  endif

end function continuity_stencil

!> Destructor for continuity_cs.
subroutine continuity_end(CS)
  type(continuity_CS), intent(inout) :: CS !< Control structure for mom_continuity.

  if (CS%continuity_scheme == PPM_SCHEME) then
    call continuity_PPM_end(CS%PPM_CSp)
  endif
end subroutine continuity_end

end module MOM_continuity
