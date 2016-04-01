module MOM_continuity
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg and Alistair Adcroft, September 2006.           *
!*                                                                     *
!*    This file contains the driver routine which selects which        *
!*  continuity solver will be called, based on run-time input.         *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, vh                                    *
!*    j    x ^ x ^ x   At >:  u, uh                                    *
!*    j    > o > o >   At o:  h, hin                                   *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_continuity_PPM, only : continuity_PPM, continuity_PPM_init
use MOM_continuity_PPM, only : continuity_PPM_end, continuity_PPM_CS
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_string_functions, only : uppercase
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : ocean_OBC_type, BT_cont_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public continuity, continuity_init, continuity_end

integer :: id_clock_pass, id_clock_vertvisc

type, public :: continuity_CS ; private
  integer :: continuity_scheme ! CONTINUITY_SCHEME selects the discretization
                               ! for the continuity solver. Valid values are:
                               !  PPM - A directionally split peicewise
                               !        parabolic reconstruction solver.
                               ! The default, PPM, seems most appropriate for
                               ! use with our current time-splitting strategies.
  type(continuity_PPM_CS), pointer :: PPM_CSp => NULL()
end type continuity_CS

integer, parameter :: PPM_SCHEME = 1
character(len=20), parameter :: PPM_STRING = "PPM"

contains

subroutine continuity(u, v, hin, h, uh, vh, dt, G, GV, CS, uhbt, vhbt, OBC, &
                      visc_rem_u, visc_rem_v, u_cor, v_cor, &
                      uhbt_aux, vhbt_aux, u_cor_aux, v_cor_aux, BT_cont)
  real, intent(in),  dimension(NIMEMB_,NJMEM_,NKMEM_) :: u
  real, intent(in),  dimension(NIMEM_,NJMEMB_,NKMEM_) :: v
  real, intent(in),  dimension(NIMEM_,NJMEM_,NKMEM_)  :: hin
  real, intent(inout), dimension(NIMEM_,NJMEM_,NKMEM_)  :: h
  real, intent(out), dimension(NIMEMB_,NJMEM_,NKMEM_) :: uh
  real, intent(out), dimension(NIMEM_,NJMEMB_,NKMEM_) :: vh
  real, intent(in)                                    :: dt
  type(ocean_grid_type), intent(inout)                :: G
  type(verticalGrid_type), intent(in)                 :: GV
  type(continuity_CS), pointer                        :: CS
  real, intent(in), optional, dimension(NIMEMB_,NJMEM_) :: uhbt
  real, intent(in), optional, dimension(NIMEM_,NJMEMB_) :: vhbt
  type(ocean_OBC_type), pointer, optional             :: OBC
  real, intent(in), optional, dimension(NIMEMB_,NJMEM_,NKMEM_) :: visc_rem_u
  real, intent(in), optional, dimension(NIMEM_,NJMEMB_,NKMEM_) :: visc_rem_v
  real, intent(out), optional, dimension(NIMEMB_,NJMEM_,NKMEM_) :: u_cor
  real, intent(out), optional, dimension(NIMEM_,NJMEMB_,NKMEM_) :: v_cor
  real, intent(in), optional, dimension(NIMEMB_,NJMEM_) :: uhbt_aux
  real, intent(in), optional, dimension(NIMEM_,NJMEMB_) :: vhbt_aux
  real, intent(inout), optional, dimension(NIMEMB_,NJMEM_,NKMEM_) :: u_cor_aux
  real, intent(inout), optional, dimension(NIMEM_,NJMEMB_,NKMEM_) :: v_cor_aux
  type(BT_cont_type),                  pointer,     optional :: BT_cont
!    This subroutine time steps the layer thicknesses, using a monotonically
!  limit, directionally split PPM scheme, based on Lin (1994).

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      hin - Initial layer thickness, in m.
!  (out)     h - Final layer thickness, in m.
!  (out)     uh - Volume flux through zonal faces = u*h*dy, m3 s-1.
!  (out)     vh - Volume flux through meridional faces = v*h*dx,
!                  in m3 s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_init.
!  (in, opt) uhbt - The summed volume flux through zonal faces, m3 s-1.
!  (in, opt) vhbt - The summed volume flux through meridional faces, m3 s-1.
!  (in, opt) OBC - This open boundary condition type specifies whether, where,
!                  and what open boundary conditions are used.
!  (in, opt) visc_rem_u - Both the fraction of the momentum originally in a
!  (in, opt) visc_rem_v - layer that remains after a time-step of viscosity,
!                         and the fraction of a time-step's worth of a
!                         barotropic acceleration that a layer experiences
!                         after viscosity is applied, in the zonal (_u) and
!                         meridional (_v) directions.  Nondimensional between
!                         0 (at the bottom) and 1 (far above the bottom).
!  (out, opt) u_cor - The zonal velocities that give uhbt as the depth-
!                     integrated transport, in m s-1.
!  (out, opt) v_cor - The meridional velocities that give vhbt as the
!                     depth-integrated transport, in m s-1.
!  (in, opt) uhbt_aux - A second set of summed volume fluxes through zonal
!  (in, opt) vhbt_aux - and meridional faces, both in m3 s-1.
!  (out, opt) u_cor_aux - The zonal and meridional velocities that give uhbt_aux
!  (out, opt) v_cor_aux - and vhbt_aux as the depth-integrated transports,
!                         both in m s-1.
!  (out, opt) BT_cont - A structure with elements that describe the effective
!                       open face areas as a function of barotropic flow.
  if (present(visc_rem_u) .neqv. present(visc_rem_v)) call MOM_error(FATAL, &
      "MOM_continuity: Either both visc_rem_u and visc_rem_v or neither"// &
       " one must be present in call to continuity.")
  if (present(u_cor) .neqv. present(v_cor)) call MOM_error(FATAL, &
      "MOM_continuity: Either both u_cor and v_cor or neither"// &
       " one must be present in call to continuity.")
  if (present(uhbt_aux) .neqv. present(vhbt_aux)) call MOM_error(FATAL, &
      "MOM_continuity: Either both uhbt_aux and uhbt_aux or neither"// &
       " one must be present in call to continuity.")
  if (present(u_cor_aux) .neqv. present(v_cor_aux)) call MOM_error(FATAL, &
      "MOM_continuity: Either both u_cor_aux and v_cor_aux or neither"// &
       " one must be present in call to continuity.")
  if (present(u_cor_aux) .neqv. present(uhbt_aux)) call MOM_error(FATAL, &
      "MOM_continuity: u_cor_aux can only be calculated if uhbt_aux is"// &
      " provided, and uhbt_aux has no other purpose.  Include both arguments"//&
      " or neither.")

  if (CS%continuity_scheme == PPM_SCHEME) then
    call continuity_PPM(u, v, hin, h, uh, vh, dt, G, GV, CS%PPM_CSp, uhbt, vhbt, OBC, &
                        visc_rem_u, visc_rem_v, u_cor, v_cor, &
                        uhbt_aux, vhbt_aux, u_cor_aux, v_cor_aux, BT_cont)
  else
    call MOM_error(FATAL, "continuity: Unrecognized value of continuity_scheme")
  endif

end subroutine continuity

subroutine continuity_init(Time, G, GV, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(verticalGrid_type), intent(in)    :: GV
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(continuity_CS),     pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_continuity" ! This module's name.
  character(len=20)  :: tmpstr

  if (associated(CS)) then
    call MOM_error(WARNING, "continuity_init called with associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "CONTINUITY_SCHEME", tmpstr, &
                 "CONTINUITY_SCHEME selects the discretization for the \n"//&
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
    call continuity_PPM_init(Time, G, GV, param_file, diag, CS%PPM_CSp)
  endif

end subroutine continuity_init

subroutine continuity_end(CS)
  type(continuity_CS),     pointer       :: CS

  if (CS%continuity_scheme == PPM_SCHEME) then
    call continuity_PPM_end(CS%PPM_CSp)
  endif

  deallocate(CS)

end subroutine continuity_end

end module MOM_continuity
