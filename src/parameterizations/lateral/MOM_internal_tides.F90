module MOM_internal_tides
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
!*  By Robert Hallberg, January 2013                                   *
!*                                                                     *
!*    This program contains the subroutines that use the ray-tracing   *
!*  equations to propagate the internal tide energy density.           *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, tauy                                  *
!*    j    x ^ x ^ x   At >:  u, taux                                  *
!*    j    > o > o >   At o:  h, fluxes.                               *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_axis_init
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_diag_mediator, only : axesType, defineAxes
use MOM_domains, only : AGRID, To_South, To_West, To_All
use MOM_domains, only : create_group_pass, do_group_pass
use MOM_domains, only : group_pass_type, start_group_pass, complete_group_pass
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : vardesc
use MOM_restart, only : register_restart_field, MOM_restart_CS
use MOM_time_manager, only : time_type, operator(+), operator(/), operator(-)
use MOM_time_manager, only : get_time, get_date, set_time, set_date
use MOM_time_manager, only : time_type_to_real
use MOM_variables, only : surface
!   Forcing is a structure containing pointers to the forcing fields
! which may be used to drive MOM.  All fluxes are positive downward.
!   Surface is a structure containing pointers to various fields that
! may be used describe the surface state of MOM.

implicit none ; private

#include <MOM_memory.h>

public propagate_int_tide, register_int_tide_restarts
public internal_tides_init, internal_tides_end

type, public :: int_tide_CS ; private
  logical :: do_int_tides    ! If true, use the internal tide code.
  integer :: nFreq = 0
  integer :: nMode = 1
  integer :: nAngle = 24
  integer :: energized_angle = -1
  logical :: upwind_1st      ! If true, use a first-order upwind scheme.
  logical :: simple_2nd      ! If true, use a simple second order (arithmetic
                             ! mean) interpolation of the edge values instead
                             ! of the higher order interpolation.
  logical :: vol_CFL         ! If true, use the ratio of the open face lengths
                             ! to the tracer cell areas when estimating CFL
                             ! numbers.  Without aggress_adjust, the default is
                             ! false; it is always true with.

  real :: decay_rate    ! A constant rate at which internal tide energy is
                        ! lost to the interior ocean internal wave field.
  real :: cdrag         ! The bottom drag coefficient for MEKE (non-dim).
  logical :: apply_drag ! If true, apply a quadratic bottom drag as a sink.
  real, allocatable, dimension(:,:,:,:,:) :: &
    En                  ! The internal wave energy density as a function of
                        ! (i,j,angle,frequency,mode)
  real, allocatable, dimension(:) :: &
    frequency           ! The frequency of each band.

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                        ! timing of diagnostic output.
  integer :: id_tot_En = -1, id_itide_drag = -1
  integer, allocatable, dimension(:,:) :: id_En_mode, id_En_ang_mode
end type int_tide_CS

type :: loop_bounds_type ; private
  integer :: ish, ieh, jsh, jeh
end type loop_bounds_type

contains

subroutine propagate_int_tide(cg1, En_in, vel_btTide, dt, G, CS)
  real, dimension(NIMEM_,NJMEM_), intent(in) :: cg1, En_in, vel_btTide
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(inout) :: G
  type(int_tide_CS), pointer       :: CS
! This subroutine calls any of the other subroutines in this file
! that are needed to specify the current surface forcing fields.
!
! Arguments: cg1 - The first mode internal gravity wave speed, in m s-1.
!  (in)      En_in - The energy input to the internal waves, in W m-2.
!  (in)      dt - Length of time over which these fluxes will be applied, in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure returned by a previous
!                 call to int_tide_init.
  real, dimension(SZI_(G),SZJ_(G),2) :: &
    test
  real, dimension(SZI_(G),SZJ_(G),CS%nMode) :: &
    c1
  real, dimension(SZI_(G),SZJB_(G)) :: &
    flux_heat_y, &
    flux_prec_y
  real, dimension(SZI_(G),SZJ_(G)) :: &
    tot_En, drag_scale
  real :: frac_per_sector, f2, I_rho0, I_D_here

  integer :: a, m, fr, i, j, is, ie, js, je, isd, ied, jsd, jed, nAngle
  type(group_pass_type), save :: pass_test, pass_En

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nAngle = CS%NAngle
  I_rho0 = 1.0 / G%Rho0

  ! Set the wave speeds for the modes, using that cg(n) ~ cg(1)/n.  This is
  ! wrong, of course, but it works reasonably in some cases.
  do m=1,CS%nMode ; do j=jsd,jed ; do i=isd,ied
    c1(i,j,m) = cg1(i,j) / real(m)
  enddo ; enddo ; enddo

  ! Add the forcing.
  if (CS%energized_angle <= 0) then
    frac_per_sector = 1.0 / real(CS%nAngle * CS%nMode * CS%nFreq)
    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
      f2 = 0.25*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                 (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J)**2))
      if (CS%frequency(fr)**2 > f2) &
        CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) + dt*frac_per_sector*En_in(i,j)
    enddo ; enddo ; enddo ; enddo ; enddo
  elseif (CS%energized_angle <= CS%nAngle) then
    frac_per_sector = 1.0 / real(CS%nMode * CS%nFreq)
    a = CS%energized_angle
    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do j=js,je ; do i=is,ie
      f2 = 0.25*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                 (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J)**2))
      if (CS%frequency(fr)**2 > f2) &
        CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) + dt*frac_per_sector*En_in(i,j)
    enddo ; enddo ; enddo ; enddo

    !### Delete this later.
    do m=1,CS%nMode ; do j=jsd,jed ; do i=isd,ied
      c1(i,j,m) = 1.0 / real(m)
    enddo ; enddo ; enddo
  else
    call MOM_error(WARNING, "Internal tide energy is being put into a angular "//&
                            "band that does not exist.")
  endif

  ! Pass a test vector to check for grid rotation in the halo updates.
  do j=jsd,jed ; do i=isd,ied ; test(i,j,1) = 1.0 ; test(i,j,2) = 0.0 ; enddo ; enddo
  do m=1,CS%nMode ; do fr=1,CS%nFreq
    call create_group_pass(pass_En, CS%En(:,:,:,fr,m), G%domain)
  enddo; enddo
  call create_group_pass(pass_test, test(:,:,1), test(:,:,2), G%domain, stagger=AGRID)
  call start_group_pass(pass_test, G%domain)
  
  ! Apply half the refraction.
  do m=1,CS%nMode ; do fr=1,CS%nFreq
    call refract(CS%En(:,:,:,fr,m), c1(:,:,m), CS%frequency(fr), 0.5*dt, G, CS%nAngle)
  enddo ; enddo
  call do_group_pass(pass_En, G%domain)

  call complete_group_pass(pass_test, G%domain)

  ! Rotate points in the halos as necessary.
  call correct_halo_rotation(CS%En, test, G, CS%nAngle)

  ! Propagate the waves.
  do m=1,CS%NMode ; do fr=1,CS%Nfreq
    call propagate(CS%En(:,:,:,fr,m), c1(:,:,m), CS%frequency(fr), dt, G, CS, CS%NAngle)
  enddo ; enddo
  
  ! Apply the other half of the refraction.
  do m=1,CS%NMode ; do fr=1,CS%Nfreq
    call refract(CS%En(:,:,:,fr,m), c1(:,:,m), CS%frequency(fr), 0.5*dt, G, CS%NAngle)
  enddo ; enddo

  if (CS%apply_drag .or. (CS%id_tot_En > 0)) then
    tot_En(:,:) = 0.0
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=js,je ; do i=is,ie
        tot_En(i,j) = tot_En(i,j) + CS%En(i,j,a,fr,m)
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  ! Extract the energy for mixing.
  if (CS%apply_drag) then
    do j=jsd,jed ; do i=isd,ied
      I_D_here = 1.0 / max(G%bathyT(i,j), 1.0)
      drag_scale(i,j) = CS%decay_rate + CS%cdrag * &
        sqrt(max(0.0, vel_btTide(i,j)**2 + tot_En(i,j) * I_rho0 * I_D_here)) * I_D_here
    enddo ; enddo
  else
    do j=jsd,jed ; do i=isd,ied
      drag_scale(i,j) = CS%decay_rate
    enddo ; enddo
  endif

  do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
    CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) / (1.0 + dt*drag_scale(i,j))
  enddo ; enddo ; enddo ; enddo ; enddo

  if (query_averaging_enabled(CS%diag)) then
    if (CS%id_tot_En > 0) call post_data(CS%id_tot_En, tot_En, CS%diag)

    if (CS%id_itide_drag > 0) call post_data(CS%id_itide_drag, drag_scale, CS%diag)

    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; if (CS%id_En_mode(fr,m) > 0) then
      tot_En(:,:) = 0.0
      do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
        tot_En(i,j) = tot_En(i,j) + CS%En(i,j,a,fr,m)
      enddo ; enddo ; enddo
      call post_data(CS%id_En_mode(fr,m), tot_En, CS%diag)
    endif ; enddo ; enddo
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; if (CS%id_En_mode(fr,m) > 0) then
      call post_data(CS%id_En_ang_mode(fr,m), CS%En(:,:,:,fr,m) , CS%diag)
    endif ; enddo ; enddo
  endif

end subroutine propagate_int_tide

subroutine refract(En, cg, freq, dt, G, NAngle)
  type(ocean_grid_type),  intent(in)    :: G
  integer,                intent(in)    :: NAngle
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), intent(inout) :: En
  real, dimension(G%isd:G%ied,G%jsd:G%jed),        intent(in)    :: cg
  real,                   intent(in)    :: freq
  real,                   intent(in)    :: dt
!  This subroutine does refraction on the internal waves at a single frequency.

! Arguments: En - the internal gravity wave energy density as a function of space
!                 and angular resolution, in J m-2 radian-1.
  integer, parameter :: stensil = 2
  real, dimension(SZI_(G),1-stensil:NAngle+stensil) :: &
    En2d
  real, dimension(1-stensil:NAngle+stensil) :: &
    cos_angle, sin_angle
  real, dimension(SZI_(G)) :: &
    Dk_Dt_Kmag, Dl_Dt_Kmag
  real, dimension(SZI_(G),0:nAngle) :: &
    Flux_E
  real :: f2   ! The squared Coriolis parameter, in s-2.
  real :: df2_dy, df2_dx  ! The x- and y- gradients of the squared Coriolis parameter, in s-2 m-1.
  real :: dlnCg_dx  ! The x-gradient of the wave speed divided by itself in m-1.
  real :: dlnCg_dy  ! The y-gradient of the wave speed divided by itself in m-1.
  real :: Angle_size, dt_Angle_size, angle
  real :: Ifreq, Kmag2, I_Kmag, CFL_ang
  real, parameter :: cg_subRO = 1e-100
  integer :: is, ie, js, je, asd, aed, na
  integer :: i, j, a

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; na = size(En,3)
  asd = 1-stensil ; aed = NAngle+stensil

  Ifreq = 1.0 / freq

  Angle_size = (8.0*atan(1.0)) / (real(NAngle))
  dt_Angle_size = dt / Angle_size

  do A=0,na
    angle = (real(A) - 0.5) * Angle_size
    cos_angle(A) = cos(angle) ; sin_angle(A) = sin(angle)
  enddo

!### There should also be refraction due to cg.grad(grid_orientation).

  do j=js,je
  ! Copy En into angle space with halos.
    do a=1,na ; do i=is,ie
      En2d(i,a) = En(i,j,a)
    enddo ; enddo
    do a=asd,0 ; do i=is,ie
      En2d(i,a) = En2d(i,a+NAngle)
      En2d(i,NAngle+stensil+a) = En2d(i,stensil+a)
    enddo ; enddo

  ! Do the refraction.
    do i=is,ie
      f2 = 0.25*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                 (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J)**2))
      df2_dx = 0.5*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I,J-1)**2) - &
                    (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I-1,J-1)**2)) * &
               G%IdxT(i,j)
      dlnCg_dx = 0.5*( G%IdxCu(I,j) * (cg(i+1,j) - cg(i,j)) / &
                       (cg(i+1,j) + cg(i,j) + cg_subRO) + &
                       G%IdxCu(I-1,j) * (cg(i,j) - cg(i-1,j)) / &
                       (cg(i,j) + cg(i-1,j) + cg_subRO) )
      df2_dy = 0.5*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J)**2) - &
                    (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J-1)**2)) * &
               G%IdyT(i,j)
      dlnCg_dy = 0.5*( G%IdyCv(i,J) * (cg(i,j+1) - cg(i,j)) / &
                       (cg(i,j+1) + cg(i,j) + cg_subRO) + &
                       G%IdyCv(i,J-1) * (cg(i,j) - cg(i,j-1)) / &
                       (cg(i,j) + cg(i,j-1) + cg_subRO) )
      Kmag2 = (freq**2 - f2) / (cg(i,j)**2 + cg_subRO**2)
      if (Kmag2 > 0.0) then
        I_Kmag = 1.0 / sqrt(Kmag2)
        Dk_Dt_Kmag(i) = -Ifreq * (df2_dx + (freq**2 - f2) * dlnCg_dx) * I_Kmag
        Dl_Dt_Kmag(i) = -Ifreq * (df2_dy + (freq**2 - f2) * dlnCg_dy) * I_Kmag
      else
        Dk_Dt_Kmag(i) = 0.0
        Dl_Dt_Kmag(i) = 0.0
      endif
    enddo
    ! Determine the energy fluxes in angular orientation space.
    do A=0,na ; do i=is,ie
      CFL_ang = (sin_angle(A) * Dk_Dt_Kmag(i) - cos_angle(A) * Dl_Dt_Kmag(i)) * &
                dt_Angle_size
      ! This is upwind for now.  Replace it with monotonic PPM later.
      if (abs(CFL_ang) > 1.0) then 
        call MOM_error(WARNING, "refract: CFL exceeds 1.", .true.)
        if (CFL_ang > 0.0) then ; CFL_ang = 1.0 ; else ; CFL_ang = -1.0 ; endif
      endif
      if (CFL_ang > 0.0) then
        Flux_E(i,A) = CFL_ang * En2d(i,A)
      else
        Flux_E(i,A) = CFL_ang * En2d(i,A+1)
      endif
    enddo ; enddo
  
  ! Update and copy back to En.
    do a=1,na ; do i=is,ie
      En(i,j,a) = En2d(i,a) + (Flux_E(i,A-1) - Flux_E(i,A))
    enddo ; enddo
  end do ! j-loop
 
end subroutine refract

subroutine propagate(En, cg, freq, dt, G, CS, NAngle)
  type(ocean_grid_type),  intent(in)    :: G
  integer,                intent(in)    :: NAngle
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), intent(inout) :: En
  real, dimension(G%isd:G%ied,G%jsd:G%jed),        intent(in)    :: cg
  real,                   intent(in)    :: freq
  real,                   intent(in)    :: dt
  type(int_tide_CS),      pointer       :: CS
!  This subroutine does refraction on the internal waves at a single frequency.

! Arguments: En - the internal gravity wave energy density as a function of space
!                 and angular resolution, in J m-2 radian-1.
  integer, parameter :: stensil = 2
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    speed_x  ! The magnitude of the group velocity at the Cu points, in m s-1.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    speed_y  ! The magnitude of the group velocity at the Cv points, in m s-1.
  real, dimension(0:NAngle) :: &
    cos_angle, sin_angle
  real, dimension(NAngle) :: &
    Cgx_av, Cgy_av, dCgx, dCgy
  real :: f2   ! The squared Coriolis parameter, in s-2.
  real :: Angle_size, I_Angle_size, angle
  real :: Ifreq, freq2
  real, parameter :: cg_subRO = 1e-100
  type(loop_bounds_type) :: LB
  integer :: is, ie, js, je, asd, aed, na
  integer :: ish, ieh, jsh, jeh
  integer :: i, j, a

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; na = size(En,3)
  asd = 1-stensil ; aed = NAngle+stensil

  Ifreq = 1.0 / freq
  freq2 = freq**2

  Angle_size = (8.0*atan(1.0)) / real(NAngle)
  I_Angle_size = 1.0 / Angle_size

  ! These could be in the control structure, as they do not vary.
  do A=0,na
    ! These are the angles at the cell edges...
    angle = (real(A) - 0.5) * Angle_size
    cos_angle(A) = cos(angle) ; sin_angle(A) = sin(angle)
  enddo

  do a=1,na
    Cgx_av(a) = (sin_angle(A) - sin_angle(A-1)) * I_Angle_size
    Cgy_av(a) = -(cos_angle(A) - cos_angle(A-1)) * I_Angle_size
    dCgx(a) = sqrt(0.5 + 0.5*(sin_angle(A)*cos_angle(A) - &
                              sin_angle(A-1)*cos_angle(A-1)) * I_Angle_size - &
                   Cgx_av(a)**2) 
    dCgy(a) = sqrt(0.5 - 0.5*(sin_angle(A)*cos_angle(A) - &
                              sin_angle(A-1)*cos_angle(A-1)) * I_Angle_size - &
                   Cgy_av(a)**2)
  enddo

  jsh = js-2 ; jeh = je+2 ; ish = is ; ieh = ie

  do j=jsh,jeh ; do I=ish-1,ieh
    f2 = 0.5*(G%CoriolisBu(I,J)**2 + G%CoriolisBu(I,J-1)**2)
    speed_x(I,j) = 0.5*(cg(i,j) + cg(i+1,j)) * G%mask2dCu(I,j) * &
                   sqrt(max(freq2 - f2, 0.0)) * Ifreq
  enddo ; enddo
  do J=jsh-1,jeh ; do i=ish,ieh
    f2 = 0.5*(G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J)**2)
    speed_y(I,j) = 0.5*(cg(i,j) + cg(i,j+1)) * G%mask2dCv(i,J) * &
                   sqrt(max(freq2 - f2, 0.0)) * Ifreq
  enddo ; enddo


  do a=1,na
    ! Apply the propagation in alternating directions.
    LB%jsh = js-2 ; LB%jeh = je+2 ; LB%ish = is ; LB%ieh = ie
    call propagate_x(En(:,:,a), speed_x, Cgx_av(a), dCgx(a), dt, G, CS, LB)

    LB%jsh = js ; LB%jeh = je ; LB%ish = is ; LB%ieh = ie
    call propagate_y(En(:,:,a), speed_y, Cgy_av(a), dCgy(a), dt, G, CS, LB)
  end do ! a-loop

end subroutine propagate
    
subroutine propagate_x(En, speed_x, Cgx_av, dCgx, dt, G, CS, LB)
  type(ocean_grid_type),  intent(in) :: G
  real, dimension(G%isd:G%ied,G%jsd:G%jed),   intent(inout) :: En
  real, dimension(G%IsdB:G%IedB,G%jsd:G%jed), intent(in)    :: speed_x
  real,                   intent(in)    :: Cgx_av, dCgx
  real,                   intent(in)    :: dt
  type(int_tide_CS),      pointer       :: CS
  type(loop_bounds_type), intent(in)    :: LB
! Arguments: En - The energy density integrated over an angular band, in W m-2,
!                 intent in/out.
!  (in)      speed_x - The magnitude of the group velocity at the Cu
!                      points, in m s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_PPM_init.
!  (in)      LB - A structure with the active energy loop bounds.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    EnL, EnR    ! Left and right face energy densities, in J m-2.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    flux_x      ! The internal wave energy flux, in J s-1.
  real, dimension(SZIB_(G)) :: &
    cg_p, cg_m, flux1, flux2
  integer :: i, j, k, ish, ieh, jsh, jeh

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh
  ! This sets EnL and EnR.
  if (CS%upwind_1st) then
    do j=jsh,jeh ; do i=ish-1,ieh+1
      EnL(i,j) = En(i,j) ; EnR(i,j) = En(i,j)
    enddo ; enddo
  else
    call PPM_reconstruction_x(En, EnL, EnR, G, LB, simple_2nd=CS%simple_2nd)
  endif

  do j=jsh,jeh
    !   This is done twice with different speeds to mimic the angular spread
    ! of the wave orientations within a sector.
    do I=ish-1,ieh
      cg_p(I) = speed_x(I,j) * (Cgx_av + dCgx)
      cg_m(I) = speed_x(I,j) * (Cgx_av - dCgx)
    enddo
    call zonal_flux_En(cg_p, En(:,j), EnL(:,j), EnR(:,j), flux1, &
                       dt, G, j, ish, ieh, CS%vol_CFL)
    call zonal_flux_En(cg_m, En(:,j), EnL(:,j), EnR(:,j), flux2, &
                       dt, G, j, ish, ieh, CS%vol_CFL)
    do I=ish-1,ieh ; flux_x(I,j) = flux1(I) + flux2(I) ; enddo
  enddo
  do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
    En(i,j) = En(i,j) - dt* G%IareaT(i,j) * (flux_x(I,j) - flux_x(I-1,j))
  enddo ; enddo

end subroutine propagate_x

    
subroutine propagate_y(En, speed_y, Cgy_av, dCgy, dt, G, CS, LB)
  type(ocean_grid_type),  intent(in) :: G
  real, dimension(G%isd:G%ied,G%jsd:G%jed),   intent(inout) :: En
  real, dimension(G%isd:G%ied,G%JsdB:G%JedB), intent(in)    :: speed_y
  real,                   intent(in)    :: Cgy_av, dCgy
  real,                   intent(in)    :: dt
  type(int_tide_CS),      pointer       :: CS
  type(loop_bounds_type), intent(in)    :: LB
! Arguments: En - The energy density integrated over an angular band, in W m-2,
!                 intent in/out.
!  (in)      speed_y - The magnitude of the group velocity at the Cv
!                      points, in m s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 continuity_PPM_init.
!  (in)      LB - A structure with the active energy loop bounds.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    EnL, EnR    ! Left and right face energy densities, in J m-2.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    flux_y      ! The internal wave energy flux, in J s-1.
  real, dimension(SZI_(G)) :: &
    cg_p, cg_m, flux1, flux2
  integer :: i, j, k, ish, ieh, jsh, jeh

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh
  ! This sets EnL and EnR.
  if (CS%upwind_1st) then
    do j=jsh-1,jeh+1 ; do i=ish,ieh
      EnL(i,j) = En(i,j) ; EnR(i,j) = En(i,j)
    enddo ; enddo
  else
    call PPM_reconstruction_y(En, EnL, EnR, G, LB, simple_2nd=CS%simple_2nd)
  endif

  do J=jsh-1,jeh
    !   This is done twice with different speeds to mimic the angular spread
    ! of the wave orientations within a sector.
    do i=ish,ieh
      cg_p(i) = speed_y(i,J) * (Cgy_av + dCgy)
      cg_m(i) = speed_y(i,J) * (Cgy_av - dCgy)
    enddo
    call merid_flux_En(cg_p, En, EnL, EnR, flux1, &
                       dt, G, J, ish, ieh, CS%vol_CFL)
    call merid_flux_En(cg_m, En, EnL, EnR, flux2, &
                       dt, G, j, ish, ieh, CS%vol_CFL)
    do i=ish,ieh ; flux_y(i,J) = flux1(i) + flux2(i) ; enddo
  enddo
  do j=jsh,jeh ; do i=ish,ieh
    En(i,j) = En(i,j) - dt* G%IareaT(i,j) * (flux_y(i,J) - flux_y(i,J-1))
  enddo ; enddo

end subroutine propagate_y

subroutine zonal_flux_En(u, h, hL, hR, uh, dt, G, j, &
                         ish, ieh, vol_CFL)
  real, dimension(NIMEMB_),    intent(in)    :: u
  real, dimension(NIMEM_),     intent(in)    :: h, hL, hR
  real, dimension(NIMEMB_),    intent(inout) :: uh
  real,                        intent(in)    :: dt
  type(ocean_grid_type),       intent(in)    :: G
  integer,                     intent(in)    :: j, ish, ieh
  logical,                     intent(in)    :: vol_CFL
!   This subroutines evaluates the zonal mass or volume fluxes in a layer.
!
! Arguments: u - Zonal velocity, in m s-1.
!  (in)      h - Energy density used to calculate the fluxes, in W m-2.
!  (in)      hL, hR - Left- and right- Energy densities in the reconstruction, in W m-2.
!  (out)     uh - The zonal energy transport, in J s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      j, ish, ieh - The index range to work on.
!  (in)      vol_CFL - If true, rescale the ratio of face areas to the cell
!                      areas when estimating the CFL number.
  real :: CFL  ! The CFL number based on the local velocity and grid spacing, ND.
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
  integer :: i

  do I=ish-1,ieh
    ! Set new values of uh and duhdu.
    if (u(I) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      else ; CFL = u(I) * dt * G%IdxT(i,j) ; endif
      curv_3 = (hL(i) + hR(i)) - 2.0*h(i)
      uh(I) = G%dy_Cu(I,j) * u(I) * &
          (hR(i) + CFL * (0.5*(hL(i) - hR(i)) + curv_3*(CFL - 1.5)))
    elseif (u(I) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I) * dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else ; CFL = -u(I) * dt * G%IdxT(i+1,j) ; endif
      curv_3 = (hL(i+1) + hR(i+1)) - 2.0*h(i+1)
      uh(I) = G%dy_Cu(I,j) * u(I) * &
          (hL(i+1) + CFL * (0.5*(hR(i+1)-hL(i+1)) + curv_3*(CFL - 1.5)))
    else
      uh(I) = 0.0
    endif
  enddo

end subroutine zonal_flux_En

subroutine merid_flux_En(v, h, hL, hR, vh, dt, G, J, &
                         ish, ieh, vol_CFL)
  real, dimension(NIMEM_),        intent(in)    :: v
  real, dimension(NIMEM_,NJMEM_), intent(in)    :: h, hL, hR
  real, dimension(NIMEM_),        intent(inout) :: vh
  real,                           intent(in)    :: dt
  type(ocean_grid_type),          intent(in)    :: G
  integer,                        intent(in)    :: J, ish, ieh
  logical,                        intent(in)    :: vol_CFL
!   This subroutines evaluates the meridional mass or volume fluxes in a layer.
!
! Arguments: v - Meridional velocity, in m s-1.
!  (in)      h - Energy density used to calculate the fluxes, in W m-2.
!  (in)      hL, hR - Left- and right- Energy densities in the reconstruction, in W m-2.
!  (out)     vh - The meridional energy transport, in J s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      J, ish, ieh - The index range to work on.
!  (in)      vol_CFL - If true, rescale the ratio of face areas to the cell
!                       areas when estimating the CFL number.
  real :: CFL ! The CFL number based on the local velocity and grid spacing, ND.
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
  integer :: i

  do i=ish,ieh
    if (v(i) > 0.0) then
      if (vol_CFL) then ; CFL = (v(i) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
      else ; CFL = v(i) * dt * G%IdyT(i,j) ; endif
      curv_3 = hL(i,j) + hR(i,j) - 2.0*h(i,j)
      vh(i) = G%dx_Cv(i,J) * v(i) * ( hR(i,j) + CFL * &
          (0.5*(hL(i,j) - hR(i,j)) + curv_3*(CFL - 1.5)) )
    elseif (v(i) < 0.0) then
      if (vol_CFL) then ; CFL = (-v(i) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
      else ; CFL = -v(i) * dt * G%IdyT(i,j+1) ; endif
      curv_3 = hL(i,j+1) + hR(i,j+1) - 2.0*h(i,j+1)
      vh(i) = G%dx_Cv(i,J) * v(i) * ( hL(i,j+1) + CFL * &
          (0.5*(hR(i,j+1)-hL(i,j+1)) + curv_3*(CFL - 1.5)) )
    else
      vh(i) = 0.0
    endif
  enddo

end subroutine merid_flux_En


subroutine correct_halo_rotation(En, test, G, NAngle)
  real, dimension(:,:,:,:,:),         intent(inout) :: En
  real, dimension(NIMEM_,NJMEM_,C2_), intent(in)    :: test
  type(ocean_grid_type),              intent(in)    :: G
  integer,                            intent(in)    :: NAngle
  !   This subroutine rotates points in the halos where required to accomodate
  ! changes in grid orientation, such as at the tripolar fold.

  real, dimension(G%isd:G%ied,NAngle) :: En2d
  integer, dimension(G%isd:G%ied) :: a_shift
  integer :: i_first, i_last, a_new
  integer :: a, i, j, isd, ied, jsd, jed, m, fr
  character(len=80) :: mesg
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed
    i_first = ied+1 ; i_last = isd-1
    do i=isd,ied
      a_shift(i) = 0
      if (test(i,j,1) /= 1.0) then
        if (i<i_first) i_first = i
        if (i>i_last) i_last = i

        if (test(i,j,1) == -1.0) then ; a_shift(i) = nAngle/2
        elseif (test(i,j,2) == 1.0) then ; a_shift(i) = -nAngle/4
        elseif (test(i,j,2) == -1.0) then ; a_shift(i) = nAngle/4
        else
          write(mesg,'("Unrecognized rotation test vector ",2ES9.2," at ",F7.2," E, ",&
                       &F7.2," N; i,j=",2i4)') &
                test(i,j,1), test(i,j,2), G%GeoLonT(i,j), G%GeoLatT(i,j), i, j
          call MOM_error(FATAL, mesg)
        endif
      endif
    enddo

    if (i_first <= i_last) then
      ! At least one point in this row needs to be rotated.
      do m=1,size(En,5) ; do fr=1,size(En,4)
        do a=1,nAngle ; do i=i_first,i_last ; if (a_shift(i) /= 0) then
          a_new = a + a_shift(i)
          if (a_new < 1) a_new = a_new + nAngle
          if (a_new > nAngle) a_new = a_new - nAngle
          En2d(i,a_new) = En(i,j,a,fr,m)
        endif ; enddo ; enddo
        do a=1,nAngle ; do i=i_first,i_last ; if (a_shift(i) /= 0) then
          En(i,j,a,fr,m) = En2d(i,a)
        endif ; enddo ; enddo
      enddo ; enddo
    endif 
  enddo

end subroutine correct_halo_rotation


subroutine PPM_reconstruction_x(h_in, h_l, h_r, G, LB, simple_2nd)
  real, dimension(NIMEM_,NJMEM_), intent(in)  :: h_in
  real, dimension(NIMEM_,NJMEM_), intent(out) :: h_l, h_r
  type(ocean_grid_type),          intent(in)  :: G
  type(loop_bounds_type),         intent(in)  :: LB
  logical, optional,              intent(in)  :: simple_2nd
! This subroutine calculates left/right edge values for PPM reconstruction.
! Arguments: h_in    - Energy density in a sector (2D)
!  (out)     h_l,h_r - left/right edge value of reconstruction (2D)
!  (in)      G - The ocean's grid structure.
!  (in)      LB - A structure with the active loop bounds.
!  (in, opt) simple_2nd - If true, use the arithmetic mean energy densities as
!                         default edge values for a simple 2nd order scheme.

! Local variables with useful mnemonic names.
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_ip1, h_im1
  real :: dMx, dMn
  logical :: use_CW84, use_2nd
  character(len=256) :: mesg
  integer :: i, j, isl, iel, jsl, jel, stensil

  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  isl = LB%ish-1 ; iel = LB%ieh+1 ; jsl = LB%jsh ; jel = LB%jeh

  ! This is the stensil of the reconstruction, not the scheme overall.
  stensil = 2 ; if (use_2nd) stensil = 1

  if ((isl-stensil < G%isd) .or. (iel+stensil > G%ied)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_x called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               stensil + max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl < G%jsd) .or. (jel > G%jed)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_x called with a ", &
               & "y-halo that needs to be increased by ",i2,".")') &
               max(G%jsd-jsl,jel-G%jed)
    call MOM_error(FATAL,mesg)
  endif

  if (use_2nd) then
    do j=jsl,jel ; do i=isl,iel
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_im1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) )
    enddo ; enddo
  else
    do j=jsl,jel ; do i=isl-1,iel+1
      if ((G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j)) == 0.0) then
        slp(i,j) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j) = 0.5 * (h_in(i+1,j) - h_in(i-1,j))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i+1,j), h_in(i-1,j), h_in(i,j)) - h_in(i,j)
        dMn = h_in(i,j) - min(h_in(i+1,j), h_in(i-1,j), h_in(i,j))
        slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j))
      endif
    enddo; enddo

    do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
    ! h_im1 = h_in(i-1,j,k) ; if (G%mask2dT(i-1,j) < 0.5) h_im1 = h_in(i,j)
    ! h_ip1 = h_in(i+1,j,k) ; if (G%mask2dT(i+1,j) < 0.5) h_ip1 = h_in(i,j)
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_l(i,j) = 0.5*( h_im1 + h_in(i,j) ) + oneSixth*( slp(i-1,j) - slp(i,j) )
      h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i+1,j) )
    enddo; enddo
  endif

  call PPM_limit_pos(h_in, h_l, h_r, 0.0, isl, iel, jsl, jel)

end subroutine PPM_reconstruction_x

subroutine PPM_reconstruction_y(h_in, h_l, h_r, G, LB, simple_2nd)
  real, dimension(NIMEM_,NJMEM_), intent(in)  :: h_in
  real, dimension(NIMEM_,NJMEM_), intent(out) :: h_l, h_r
  type(ocean_grid_type),          intent(in)  :: G
  type(loop_bounds_type),         intent(in)  :: LB
  logical, optional,              intent(in)  :: simple_2nd
! This subroutine calculates left/right edge valus for PPM reconstruction.
! Arguments: h_in    - Energy density in a sector (2D)
!  (out)     h_l,h_r - left/right edge value of reconstruction (2D)
!  (in)      G - The ocean's grid structure.
!  (in)      LB - A structure with the active loop bounds.
!  (in, opt) simple_2nd - If true, use the arithmetic mean energy densities as
!                         default edge values for a simple 2nd order scheme.

! Local variables with useful mnemonic names.
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_jp1, h_jm1
  real :: dMx, dMn
  logical :: use_2nd
  character(len=256) :: mesg
  integer :: i, j, isl, iel, jsl, jel, stensil

  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  isl = LB%ish ; iel = LB%ieh ; jsl = LB%jsh-1 ; jel = LB%jeh+1

  ! This is the stensil of the reconstruction, not the scheme overall.
  stensil = 2 ; if (use_2nd) stensil = 1

  if ((isl < G%isd) .or. (iel > G%ied)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_y called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl-stensil < G%jsd) .or. (jel+stensil > G%jed)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_y called with a ", &
                 & "y-halo that needs to be increased by ",i2,".")') &
                 stensil + max(G%jsd-jsl,jel-G%jed)
    call MOM_error(FATAL,mesg)
  endif

  if (use_2nd) then
    do j=jsl,jel ; do i=isl,iel
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) )
    enddo ; enddo
  else
    do j=jsl-1,jel+1 ; do i=isl,iel
      if ((G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1)) == 0.0) then
        slp(i,j) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j) = 0.5 * (h_in(i,j+1) - h_in(i,j-1))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i,j+1), h_in(i,j-1), h_in(i,j)) - h_in(i,j)
        dMn = h_in(i,j) - min(h_in(i,j+1), h_in(i,j-1), h_in(i,j))
        slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1))
      endif
    enddo ; enddo

    do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) ) + oneSixth*( slp(i,j-1) - slp(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i,j+1) )
    enddo ; enddo
  endif

  call PPM_limit_pos(h_in, h_l, h_r, 0.0, isl, iel, jsl, jel)

end subroutine PPM_reconstruction_y


subroutine PPM_limit_pos(h_in, h_L, h_R, h_min, iis, iie, jis, jie)
  real, dimension(NIMEM_,NJMEM_), intent(in)     :: h_in
  real, dimension(NIMEM_,NJMEM_), intent(inout)  :: h_L, h_R
  real,                           intent(in)     :: h_min
  integer,                        intent(in)     :: iis, iie, jis, jie
! This subroutine limits the left/right edge values of the PPM reconstruction
! to give a reconstruction that is positive-definite.  Here this is
! reinterpreted as giving a constant thickness if the mean thickness is less
! than h_min, with a minimum of h_min otherwise.
! Arguments: h_in    - thickness of layer (2D)
!  (inout)   h_L     - left edge value (2D)
!  (inout)   h_R     - right edge value (2D)
!  (in)      h_min   - The minimum thickness that can be obtained by a
!                      concave parabolic fit.
!  (in)      iis, iie, jis, jie - Index range for computation.

! Local variables
  real    :: curv, dh, scale
  character(len=256) :: mesg
  integer :: i,j

  do j=jis,jie ; do i=iis,iie
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    curv = 3.0*(h_L(i,j) + h_R(i,j) - 2.0*h_in(i,j))
    if (curv > 0.0) then ! Only minima are limited.
      dh = h_R(i,j) - h_L(i,j)
      if (abs(dh) < curv) then ! The parabola's minimum is within the cell.
        if (h_in(i,j) <= h_min) then
          h_L(i,j) = h_in(i,j) ; h_R(i,j) = h_in(i,j)
        elseif (12.0*curv*(h_in(i,j) - h_min) < (curv**2 + 3.0*dh**2)) then
          ! The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
          ! be limited in this case.  0 < scale < 1.
          scale = 12.0*curv*(h_in(i,j) - h_min) / (curv**2 + 3.0*dh**2)
          h_L(i,j) = h_in(i,j) + scale*(h_L(i,j) - h_in(i,j))
          h_R(i,j) = h_in(i,j) + scale*(h_R(i,j) - h_in(i,j))
        endif
      endif
    endif
  enddo ; enddo

end subroutine PPM_limit_pos


subroutine register_int_tide_restarts(G, param_file, CS, restart_CS)
  type(ocean_grid_type), intent(in) :: G
  type(param_file_type), intent(in) :: param_file
  type(int_tide_CS),     pointer :: CS
  type(MOM_restart_CS),  pointer :: restart_CS
! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.
!  (in)      restart_CS - A pointer to the restart control structure.
! This subroutine is used to allocate and register any fields in this module
! that should be written to or read from the restart file.
  logical :: use_int_tides, use_temperature
  character (len=8) :: period_str
  type(vardesc) :: vd
  integer :: num_freq, num_angle , num_mode
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, a
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "register_int_tide_restarts called "//&
                             "with an associated control structure.")
    return
  endif

  use_int_tides = .false.
  call read_param(param_file, "INTERNAL_TIDES", use_int_tides)
  if (.not.use_int_tides) return

  use_temperature = .true.
  call read_param(param_file, "ENABLE_THERMODYNAMICS", use_temperature)
  if (.not.use_temperature) call MOM_error(FATAL, &
    "register_int_tide_restarts: internal_tides only works with "//&
    "ENABLE_THERMODYNAMICS defined.")

  num_freq = 0 ; num_angle = 24 ; num_mode = 1
  call read_param(param_file, "INTERNAL_TIDE_FREQS", num_freq)
  call read_param(param_file, "INTERNAL_TIDE_ANGLES", num_angle)
  call read_param(param_file, "INTERNAL_TIDE_MODES", num_mode)

  if (.not.((num_freq > 0) .and. (num_angle > 0) .and. (num_mode > 0))) return

  ! The internal tide code will be used.
  allocate(CS)
  CS%do_int_tides = .true.
  CS%nFreq = num_freq ; CS%nAngle = num_angle ; CS%nMode = num_mode

  allocate(CS%En(isd:ied, jsd:jed, num_angle, num_freq, num_mode))
  CS%En(:,:,:,:,:) = 0.0

  allocate(CS%frequency(num_freq))
  do a=1,num_freq  ! freq = 2 pi a / 1 day
    CS%frequency(a) = 8.0*atan(1.0) * (real(a)) / 86400.0
  enddo

! if (CS%do_integrated) then
!   vd = vardesc("Ctrl_heat","Control Integrative Heating",'h','1','s',"W m-2")
!   call register_restart_field(CS%heat_0, vd, .false., restart_CS)
! endif

end subroutine register_int_tide_restarts


subroutine internal_tides_init(Time, G, param_file, diag, CS)
  type(time_type),           intent(in) :: Time
  type(ocean_grid_type),     intent(in) :: G
  type(param_file_type),     intent(in) :: param_file
  type(diag_ctrl), target,   intent(in) :: diag
  type(int_tide_CS),     pointer    :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  real    :: angle_size
  real, allocatable :: angles(:)
  logical :: use_int_tides, use_temperature
  integer :: num_angle, num_freq, num_mode, m, fr
  integer :: isd, ied, jsd, jed, a, id_ang
  type(axesType) :: axes_ang
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_internal_tides" ! This module's name.
  character(len=16), dimension(8) :: freq_name
  character(len=40)  :: var_name
  character(len=160) :: var_descript

  ! These should have already been called.
  ! call read_param(param_file, "INTERNAL_TIDES", controlled)

  if (associated(CS)) then
    use_int_tides = CS%do_int_tides
  else
!   use_int_tides = .false.

! ###This is here just until the restart code is ready.
      isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
      use_int_tides = .false.
      call read_param(param_file, "INTERNAL_TIDES", use_int_tides)
      if (.not.use_int_tides) return

      use_temperature = .true.
      call read_param(param_file, "ENABLE_THERMODYNAMICS", use_temperature)
      if (.not.use_temperature) call MOM_error(FATAL, &
        "register_int_tide_restarts: internal_tides only works with "//&
        "ENABLE_THERMODYNAMICS defined.")

      num_freq = 1 ; num_angle = 24 ; num_mode = 1
      call read_param(param_file, "INTERNAL_TIDE_FREQS", num_freq)
      call read_param(param_file, "INTERNAL_TIDE_ANGLES", num_angle)
      call read_param(param_file, "INTERNAL_TIDE_MODES", num_mode)

      if (.not.((num_freq > 0) .and. (num_angle > 0) .and. (num_mode > 0))) return

      ! The internal tide code will be used.
      allocate(CS)
      CS%do_int_tides = .true. ;  use_int_tides = .true.
      CS%nFreq = num_freq ; CS%nAngle = num_angle ; CS%nMode = num_mode

      allocate(CS%En(isd:ied, jsd:jed, num_angle, num_freq, num_mode))
      CS%En(:,:,:,:,:) = 0.0

      allocate(CS%frequency(num_freq))
      do a=1,num_freq  ! freq = 2 pi a / 1 day
        CS%frequency(a) = 8.0*atan(1.0) * (real(a)) / 86400.0
      enddo

  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "INTERNAL_TIDE_FREQS", num_freq, &
                 "The number of distinct internal tide frequency bands \n"//&
                 "that will be calculated.", default=1)
  call get_param(param_file, mod, "INTERNAL_TIDE_MODES", num_mode, &
                 "The number of distinct internal tide modes \n"//&
                 "that will be calculated.", default=1)
  call get_param(param_file, mod, "INTERNAL_TIDE_ANGLES", num_angle, &
                 "The number of angular resolution bands for the internal \n"//&
                 "tide calculations.", default=24)

  if (use_int_tides) then
    if ((num_freq <= 0) .and. (num_mode <= 0) .and. (num_angle <= 0)) then
      call MOM_error(WARNING, "Internal tides were enabled, but the number "//&
             "of requested frequencies, modes and angles were not all positive.")
      return
    endif
  else
    if ((num_freq > 0) .and. (num_mode > 0) .and. (num_angle > 0)) then
      call MOM_error(WARNING, "Internal tides were not enabled, even though "//&
             "the number of requested frequencies, modes and angles were all "//&
             "positive.")
      return
    endif
  endif
  
  if (CS%NFreq /= num_freq) call MOM_error(FATAL, "Internal_tides_init: "//&
      "Inconsistent number of frequencies.")
  if (CS%NAngle /= num_angle) call MOM_error(FATAL, "Internal_tides_init: "//&
      "Inconsistent number of angles.")
  if (CS%NMode /= num_mode) call MOM_error(FATAL, "Internal_tides_init: "//&
      "Inconsistent number of angles.")
  if (4*(num_angle/4) /= num_angle) call MOM_error(FATAL, &
    "Internal_tides_init: INTERNAL_TIDE_ANGLES must be a multiple of 4.")

  CS%diag => diag
  
  call get_param(param_file, mod, "INTERNAL_TIDE_DECAY_RATE", CS%decay_rate, &
                 "The rate at which internal tide energy is lost to the \n"//&
                 "interior ocean internal wave field.", units="s-1", default=0.0)
  call get_param(param_file, mod, "INTERNAL_TIDE_VOLUME_BASED_CFL", CS%vol_CFL, &
                 "If true, use the ratio of the open face lengths to the \n"//&
                 "tracer cell areas when estimating CFL numbers in the \n"//&
                 "internal tide code.", default=.false.)
  call get_param(param_file, mod, "INTERNAL_TIDE_SIMPLE_2ND_PPM", CS%simple_2nd, &
                 "If true, CONTINUITY_PPM uses a simple 2nd order \n"//&
                 "(arithmetic mean) interpolation of the edge values. \n"//&
                 "This may give better PV conservation propterties. While \n"//&
                 "it formally reduces the accuracy of the continuity \n"//&
                 "solver itself in the strongly advective limit, it does \n"//&
                 "not reduce the overall order of accuracy of the dynamic \n"//&
                 "core.", default=.false.)
  call get_param(param_file, mod, "INTERNAL_TIDE_UPWIND_1ST", CS%upwind_1st, &
                 "If true, the internal tide ray-tracing advection uses \n"//&
                 "1st-order upwind advection.  This scheme is highly \n"//&
                 "continuity solver.  This scheme is highly \n"//&
                 "diffusive but may be useful for debugging.", default=.false.)
  call get_param(param_file, mod, "INTERNAL_TIDE_QUAD_DRAG", CS%apply_drag, &
                 "If true, the internal tide ray-tracing advection uses \n"//&
                 "a quadratic bottom drag term as a sink.", default=.true.)
  call get_param(param_file, mod, "CDRAG", CS%cdrag, &
                 "CDRAG is the drag coefficient relating the magnitude of \n"//&
                 "the velocity field to the bottom stress.", units="nondim", &
                 default=0.003)
  call get_param(param_file, mod, "INTERNAL_TIDE_ENERGIZED_ANGLE", CS%energized_angle, &
                 "If positive, only one angular band of the internal tides \n"//&
                 "gets all of the energy.  (This is for debugging.)", default=-1)

  CS%id_tot_En = register_diag_field('ocean_model', 'ITide_tot_En', diag%axesT1, &
                 Time, 'Internal tide total energy density', 'J m-2')
  CS%id_itide_drag = register_diag_field('ocean_model', 'ITide_drag', diag%axesT1, &
                 Time, 'Interior and bottom drag internal tide decay timescale', 's-1')

  allocate(CS%id_En_mode(CS%nFreq,CS%nMode)) ; CS%id_En_mode(:,:) = -1
  allocate(CS%id_En_ang_mode(CS%nFreq,CS%nMode)) ; CS%id_En_ang_mode(:,:) = -1

  allocate(angles(CS%NAngle)) ; angles(:) = 0.0
  Angle_size = (8.0*atan(1.0)) / (real(num_angle))
  do a=1,num_angle ; angles(a) = (real(a) - 1) * Angle_size ; enddo

  id_ang = diag_axis_init("angle", angles, "Radians", "N", "Angular Orienation of Fluxes")
  call defineAxes(diag, (/ diag%axesT1%handles(1), diag%axesT1%handles(2), id_ang /), axes_ang)

  do fr=1,CS%nFreq ; write(freq_name(fr), '("K",i1)') fr ; enddo
  do m=1,CS%nMode ; do fr=1,CS%nFreq
    write(var_name, '("Itide_en_K",i1,"_M",i1)') fr, m
    write(var_descript, '("Internal tide energy density in frequency K",i1," mode ",i1)') fr, m
    CS%id_En_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'J m-2')
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    write(var_name, '("Itide_en_ang_K",i1,"_M",i1)') fr, m
    write(var_descript, '("Internal tide angular energy density in frequency K",i1," mode ",i1)') fr, m
    CS%id_En_ang_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 axes_ang, Time, var_descript, 'J m-2 band-1')
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
  enddo ; enddo

end subroutine internal_tides_init


subroutine internal_tides_end(CS)
  type(int_tide_CS), pointer :: CS
! Arguments:  CS - A pointer to the control structure returned by a previous
!                  call to internal_tides_init, it will be deallocated here.

  if (associated(CS)) then
    if (allocated(CS%En)) deallocate(CS%En)
    if (allocated(CS%frequency)) deallocate(CS%frequency)
    if (allocated(CS%id_En_mode)) deallocate(CS%id_En_mode)
    deallocate(CS)
  endif
  CS => NULL()

end subroutine internal_tides_end

end module MOM_internal_tides
