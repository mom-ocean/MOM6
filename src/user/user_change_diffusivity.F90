module user_change_diffusivity
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

use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs, vertvisc_type, p3d
use MOM_EOS, only : calculate_density

implicit none ; private

#include <MOM_memory.h>

public user_change_diff, user_change_diff_init
public user_change_diff_end

type, public :: user_change_diff_CS ; private
  real :: Kd_add        !   The scale of a diffusivity that is added everywhere
                        ! without any filtering or scaling, in m2 s-1.
  real :: lat_range(4)  !   4 values that define the latitude range over which
                        ! a diffusivity scaled by Kd_add is added, in deg.
  real :: rho_range(4)  !   4 values that define the coordinate potential
                        ! density range over which a diffusivity scaled by
                        ! Kd_add is added, in kg m-3.
  logical :: use_abs_lat  ! If true, use the absolute value of latitude when
                          ! setting lat_range.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
end type user_change_diff_CS

contains

!> This subroutine provides an interface for a user to use to modify the
!! main code to alter the diffusivities as needed.  The specific example
!! implemented here augments the diffusivity for a specified range of latitude
!! and coordinate potential density.
subroutine user_change_diff(h, tv, G, CS, Kd, Kd_int, T_f, S_f, Kd_int_add)
  type(ocean_grid_type),                    intent(in)    :: G   !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h   !< Layer thickness, in m or kg m-2.
  type(thermo_var_ptrs),                    intent(in)    :: tv  !< A structure containing pointers
                                                                 !! to any available thermodynamic
                                                                 !! fields. Absent fields have NULL ptrs.
  type(user_change_diff_CS),                pointer       :: CS  !< This module's control structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   optional, intent(inout) :: Kd !< The diapycnal diffusivity of
                                                                  !! each layer in m2 s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), optional, intent(inout) :: Kd_int !< The diapycnal diffusivity
                                                                  !! at each interface in m2 s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   optional, intent(in)    :: T_f !< Temperature with massless
                                                                  !! layers filled in vertically.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   optional, intent(in)    :: S_f !< Salinity with massless
                                                                  !! layers filled in vertically.
  real, dimension(:,:,:),                     optional, pointer       :: Kd_int_add !< The diapycnal
                                                                  !! diffusivity that is being added at
                                                                  !! each interface in m2 s-1.

  real :: Rcv(SZI_(G),SZK_(G)) ! The coordinate density in layers in kg m-3.
  real :: p_ref(SZI_(G))       ! An array of tv%P_Ref pressures.
  real :: rho_fn      ! The density dependence of the input function, 0-1, ND.
  real :: lat_fn      ! The latitude dependence of the input function, 0-1, ND.
  logical :: use_EOS  ! If true, density is calculated from T & S using an
                      ! equation of state.
  logical :: store_Kd_add  ! Save the added diffusivity as a diagnostic if true.
  integer :: i, j, k, is, ie, js, je, nz
  integer :: isd, ied, jsd, jed

  real :: kappa_fill  ! diffusivity used to fill massless layers
  real :: dt_fill     ! timestep used to fill massless layers
  character(len=200) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) call MOM_error(FATAL,"user_set_diffusivity: "//&
         "Module must be initialized before it is used.")

  use_EOS = associated(tv%eqn_of_state)
  if (.not.use_EOS) return
  store_Kd_add = .false.
  if (present(Kd_int_add)) store_Kd_add = associated(Kd_int_add)

  if (.not.range_OK(CS%lat_range)) then
    write(mesg, '(4(1pe15.6))') CS%lat_range(1:4)
    call MOM_error(FATAL, "user_set_diffusivity: bad latitude range: \n  "//&
                    trim(mesg))
  endif
  if (.not.range_OK(CS%rho_range)) then
    write(mesg, '(4(1pe15.6))') CS%rho_range(1:4)
    call MOM_error(FATAL, "user_set_diffusivity: bad density range: \n  "//&
                    trim(mesg))
  endif

  if (store_Kd_add) Kd_int_add(:,:,:) = 0.0

  do i=is,ie ; p_ref(i) = tv%P_Ref ; enddo
  do j=js,je
    if (present(T_f) .and. present(S_f)) then
      do k=1,nz
        call calculate_density(T_f(:,j,k),S_f(:,j,k),p_ref,Rcv(:,k),&
                               is,ie-is+1,tv%eqn_of_state)
      enddo
    else
      do k=1,nz
        call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p_ref,Rcv(:,k),&
                               is,ie-is+1,tv%eqn_of_state)
      enddo
    endif

    if (present(Kd)) then
      do k=1,nz ; do i=is,ie
        if (CS%use_abs_lat) then
          lat_fn = val_weights(abs(G%geoLatT(i,j)), CS%lat_range)
        else
          lat_fn = val_weights(G%geoLatT(i,j), CS%lat_range)
        endif
        rho_fn = val_weights(Rcv(i,k), CS%rho_range)
        if (rho_fn * lat_fn > 0.0) &
          Kd(i,j,k) = Kd(i,j,k) + CS%Kd_add * rho_fn * lat_fn
      enddo ; enddo
    endif
    if (present(Kd_int)) then
      do K=2,nz ; do i=is,ie
        if (CS%use_abs_lat) then
          lat_fn = val_weights(abs(G%geoLatT(i,j)), CS%lat_range)
        else
          lat_fn = val_weights(G%geoLatT(i,j), CS%lat_range)
        endif
        !   rho_int = 0.5*(Rcv(i,k-1) + Rcv(i,k))
        rho_fn = val_weights( 0.5*(Rcv(i,k-1) + Rcv(i,k)), CS%rho_range)
        if (rho_fn * lat_fn > 0.0) then
          Kd_int(i,j,K) = Kd_int(i,j,K) + CS%Kd_add * rho_fn * lat_fn
          if (store_Kd_add) Kd_int_add(i,j,K) = CS%Kd_add * rho_fn * lat_fn
        endif
      enddo ; enddo
    endif
  enddo

end subroutine user_change_diff

!> This subroutine checks whether the 4 values of range are in ascending order.
function range_OK(range) result(OK)
  real, dimension(4), intent(in) :: range  !< Four values to check.
  logical                        :: OK     !< Return value.


  OK = ((range(1) <= range(2)) .and. (range(2) <= range(3)) .and. &
        (range(3) <= range(4)))

end function range_OK

!> This subroutine returns a value that goes smoothly from 0 to 1, stays
!! at 1, and then goes smoothly back to 0 at the four values of range.  The
!! transitions are cubic, and have zero first derivatives where the curves
!! hit 0 and 1.  The values in range must be in ascending order, as can be
!! checked by calling range_OK.
function val_weights(val, range) result(ans)
  real,               intent(in) :: val    !< Value for which we need an answer.
  real, dimension(4), intent(in) :: range  !< Range over which the answer is non-zero.
  real                           :: ans    !< Return value.

  real :: x   ! A nondimensional number between 0 and 1.

  ans = 0.0
  if ((val > range(1)) .and. (val < range(4))) then
    if (val < range(2)) then
      ! x goes from 0 to 1; ans goes from 0 to 1, with 0 derivatives at the ends.
      x = (val - range(1)) / (range(2) - range(1))
      ans = x**2 * (3.0 - 2.0 * x)
    elseif (val > range(3)) then
      ! x goes from 0 to 1; ans goes from 0 to 1, with 0 derivatives at the ends.
      x = (range(4) - val) / (range(4) - range(3))
      ans = x**2 * (3.0 - 2.0 * x)
    else
      ans = 1.0
    endif
  endif

end function val_weights

!> Set up the module control structure.
subroutine user_change_diff_init(Time, G, param_file, diag, CS)
  type(time_type),           intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),     intent(in)    :: G          !< The ocean's grid structure.
  type(param_file_type),     intent(in)    :: param_file !< A structure indicating the
                                                         !! open file to parse for
                                                         !! model parameter values.
  type(diag_ctrl), target,   intent(inout) :: diag       !< A structure that is used to
                                                         !! regulate diagnostic output.
  type(user_change_diff_CS), pointer       :: CS         !< A pointer that is set to
                                                         !! point to the control
                                                         !! structure for this module.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "user_set_diffusivity"  ! This module's name.
  character(len=200) :: mesg
  integer :: i, j, is, ie, js, je

  if (associated(CS)) then
    call MOM_error(WARNING, "diabatic_entrain_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "USER_KD_ADD", CS%Kd_add, &
                 "A user-specified additional diffusivity over a range of \n"//&
                 "latitude and density.", units="m2 s-1", default=0.0)
  if (CS%Kd_add /= 0.0) then
    call get_param(param_file, mod, "USER_KD_ADD_LAT_RANGE", CS%lat_range(:), &
                 "Four successive values that define a range of latitudes \n"//&
                 "over which the user-specified extra diffusivity is \n"//&
                 "applied.  The four values specify the latitudes at \n"//&
                 "which the extra diffusivity starts to increase from 0, \n"//&
                 "hits its full value, starts to decrease again, and is \n"//&
                 "back to 0.", units="degree", default=-1.0e9)
    call get_param(param_file, mod, "USER_KD_ADD_RHO_RANGE", CS%rho_range(:), &
                 "Four successive values that define a range of potential \n"//&
                 "densities over which the user-given extra diffusivity \n"//&
                 "is applied.  The four values specify the density at \n"//&
                 "which the extra diffusivity starts to increase from 0, \n"//&
                 "hits its full value, starts to decrease again, and is \n"//&
                 "back to 0.", units="kg m-3", default=-1.0e9)
    call get_param(param_file, mod, "USER_KD_ADD_USE_ABS_LAT", CS%use_abs_lat, &
                 "If true, use the absolute value of latitude when \n"//&
                 "checking whether a point fits into range of latitudes.", &
                 default=.false.)
  endif

 if (.not.range_OK(CS%lat_range)) then
    write(mesg, '(4(1pe15.6))') CS%lat_range(1:4)
    call MOM_error(FATAL, "user_set_diffusivity: bad latitude range: \n  "//&
                    trim(mesg))
  endif
  if (.not.range_OK(CS%rho_range)) then
    write(mesg, '(4(1pe15.6))') CS%rho_range(1:4)
    call MOM_error(FATAL, "user_set_diffusivity: bad density range: \n  "//&
                    trim(mesg))
  endif

end subroutine user_change_diff_init

!> Clean up the module control structure.
subroutine user_change_diff_end(CS)
  type(user_change_diff_CS), pointer :: CS         !< A pointer that is set to
                                                   !! point to the control
                                                   !! structure for this module.

  if (associated(CS)) deallocate(CS)

end subroutine user_change_diff_end

!> \class user_change_diffusivity
!!
!!  By Robert Hallberg, May 2012
!!
!!    This file contains a subroutine that increments the diapycnal
!!  diffusivity in a specified band of latitudes and densities.
!!
!!     A small fragment of the grid is shown below:
!!
!!    j+1  x ^ x ^ x   At x:  q
!!    j+1  > o > o >   At ^:  v
!!    j    x ^ x ^ x   At >:  u
!!    j    > o > o >   At o:  h, T, S, Kd, etc.
!!    j-1  x ^ x ^ x
!!        i-1  i  i+1  At x & ^:
!!           i  i+1    At > & o:
!!
!!  The boundaries always run through q grid points (x).

end module user_change_diffusivity
