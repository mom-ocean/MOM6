module MOM_wave_structure
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
!
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Ben Mater, September, 2015                                      *
!*                                                                     *
!*    The subroutine in this module calculates the vertical structure  *
!*    functions of the first baroclinic mode internal wave speed.      *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, vh, vav                               *
!*    j    x ^ x ^ x   At >:  u, uh, uav                               *
!*    j    > o > o >   At o:  h                                        *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

!use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
!use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
!use MOM_error_handler, only : MOM_error, FATAL, WARNING
!use MOM_file_parser, only : log_version, param_file_type
use MOM_grid, only : ocean_grid_type

implicit none ; private

#include <MOM_memory.h>

public wave_structure, wave_structure_init

type, public :: wave_structure_CS ; private
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                                   ! timing of diagnostic output. 
end type wave_structure_CS

contains

subroutine wave_structure(cg1, G, CS, full_halos, Igu_map, Igl_map, wmode, umode)
  real, dimension(NIMEM_,NJMEM_),        intent(in)  :: cg1
  type(ocean_grid_type),                 intent(in)  :: G
  type(wave_speed_CS), optional,         pointer     :: CS
  logical,             optional,         intent(in)  :: full_halos
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)  :: Igl_map
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)  :: Igu_map
!    This subroutine determines the first mode internal wave velocity structure.
! Arguments: 
!  (in)      cg1 - The first mode internal gravity wave speed, in m s-1.
!  (in)      Igu_map - The inverse of the reduced gravity across an interface times
!                      the thickness above it for all columns, in  s2 m-2.
!  (in)      Igl_map - The inverse of the reduced gravity across an interface times
!                      the thickness below it for all columns, in  s2 m-2.
!  (in)      G - The ocean's grid structure.
!  (out)     wmode - Vertical structure of vertical velocity (normalized), in m s-1.
!  (out)     umode - Vertical structure of horizontal velocity (normalized), in m s-1.

! This subroutine solves for the eigen vector [vertical structure, e(k)] associated with 
! the first baroclinic mode speed [i.e., smallest eigen value (lam = 1/c^2)] of the 
! system d2e/dz2 = -(N2/cn2)e, or (A-lam*I)e = 0, where A = -(1/N2)(d2/dz2), lam = 1/c^2,
! and I is the identity matrix. 2nd order discretization in the vertical lets this system
! be represented as 
!
!   -Igu(k)*e(k-1) + (Igu(k)+Igl(k)-lam)*e(k) - Igl(k)*e(k+1) = 0.0
!
! with rigid lid boundary conditions e(1) = e(nz+1) = 0.0 giving
!
!   (Igu(2)+Igl(2)-lam)*e(2) - Igl(2)*e(3) = 0.0
!   -Igu(nz)*e(nz-1) + (Igu(nz)+Igl(nz)-lam)*e(nz) = 0.0
!
! where, upon noting N2 = reduced gravity/layer thickness, we get
!    Igl(k) = 1.0/(gprime(k)*H(k)) ; Igu(k) = 1.0/(gprime(k)*H(k-1))
!
! The eigen value for this system is approximated using "wave_speed." This subroutine uses 
! these eigen values (mode speeds) to estimate the corresponding eigen vectors (velocity
! structure) using the "inverse iteration with shift" method. The algorithm is 
! 
!   Pick a starting vector reasonably close to mode structure and with unit magnitude, b_guess
!   For n=1,2,3,...
!     Solve (A-lam*I)e = e_guess for e
!     Set e_guess=e/|e| and repeat, with each iteration refining the estimate of e

  real, dimension(SZK_(G))   :: e_guess
                  ! initial guess at eigen vector with unit magnitude
  !integer, parameter :: max_itt = 10
  !real :: lam_it(max_itt), det_it(max_itt), ddet_it(max_itt)
  !integer :: kc
  integer :: i, j, k, k2, itt, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (present(CS)) then
    if (.not. associated(CS)) call MOM_error(FATAL, "MOM_wave_structure: "// &
           "Module must be initialized before it is used.")
  endif
  if (present(full_halos)) then ; if (full_halos) then
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif ; endif

  do j=js,je
    do i=is,ie
      lam = 1/cg1(i,j)**2
      z = 
      e_guess = cos()
      e_guess = e_guess/sqrt(sum(e_guess**2))
      do itt=1,max_itt
                
      
    endif ; enddo ! i-loop
  enddo ! j-loop

end subroutine wave_speed

subroutine wave_structure_init(Time, G, param_file, diag, CS)
  type(time_type),             intent(in)    :: Time
  type(ocean_grid_type),       intent(in)    :: G
  type(param_file_type),       intent(in)    :: param_file
  type(diag_ctrl), target,     intent(inout) :: diag
  type(wave_structure_CS),     pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                  for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_wave_structure"  ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "wave_structure_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag

  ! Write all relevant parameters to the model log.
  call log_version(param_file, mod, version, "")

end subroutine wave_structure_init

end module MOM_wave_structure
