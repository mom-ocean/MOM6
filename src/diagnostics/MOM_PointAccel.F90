module MOM_PointAccel
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

!***********************************************************************
!*                                                                     *
!*     The two subroutines in this file write out all of the terms     *
!*  in the u- or v-momentum balance at a given point.  Usually         *
!*  these subroutines are called after the velocities exceed some      *
!*  threshold, in order to determine which term is culpable.           *
!*  often this is done for debugging purposes.                         *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h  *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q, CoriolisBu                            *
!*    j+1  > o > o >   At ^:  v, PFv, CAv, vh, diffv, vbt, vhtr        *
!*    j    x ^ x ^ x   At >:  u, PFu, CAu, uh, diffu, ubt, uhtr        *
!*    j    > o > o >   At o:  h, bathyT, tr, T, S                      *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator, only : diag_ctrl
use MOM_domains, only : pe_here
use MOM_error_handler, only : MOM_error, NOTE
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : open_file
use MOM_io, only : APPEND_FILE, ASCII_FILE, MULTIPLE, SINGLE_FILE
use MOM_time_manager, only : time_type, get_time, get_date, set_date, operator(-)
use MOM_variables, only : ocean_internal_state, accel_diag_ptrs, cont_diag_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public write_u_accel, write_v_accel, PointAccel_init

type, public :: PointAccel_CS ; private
  character(len=200) :: u_trunc_file ! The complete path to files in which a
  character(len=200) :: v_trunc_file ! column's worth of accelerations are
                                     ! written if velocity truncations occur.
  integer :: u_file, v_file ! The unit numbers for opened u- or v- truncation
                            ! files, or -1 if they have not yet been opened.
  integer :: cols_written   ! The number of columns whose output has been
                            ! written by this PE during the current run.
  integer :: max_writes     ! The maximum number of times any PE can write out
                            ! a column's worth of accelerations during a run.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag ! A pointer to a structure of shareable
                            ! ocean diagnostic fields.
! The following are pointers to many of the state variables and accelerations
! that are used to step the physical model forward.  They all use the same
! names as the variables they point to in MOM.F90
  real, pointer, dimension(:,:,:) :: &
    u_av => NULL(), v_av => NULL(), & ! Time average velocities in m s-1.
    u_prev => NULL(), v_prev => NULL(), & ! Previous velocities in m s-1.
    T => NULL(), S => NULL(), &     ! Temperature and salinity in C and psu.
    pbce => NULL(), &               ! pbce times eta gives the baroclinic
                                    ! pressure anomaly in each layer due to
                                    ! free surface height anomalies.
                                    ! pbce has units of m s-2.
    u_accel_bt => NULL(), &         ! Barotropic acclerations in m s-2.
    v_accel_bt => NULL()

end type PointAccel_CS

contains

subroutine write_u_accel(I, j, um, hin, ADp, CDp, dt, G, GV, CS, &
                         maxvel, minvel, str, a, hv)
  integer,                                intent(in) :: I, j
  type(ocean_grid_type),                  intent(in) :: G
  type(verticalGrid_type),                intent(in) :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in) :: um
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: hin
  type(accel_diag_ptrs),                  intent(in) :: ADp
  type(cont_diag_ptrs),                   intent(in) :: CDp
  real,                                   intent(in) :: dt
  type(PointAccel_CS),                    pointer    :: CS
  real,                                   intent(in) :: maxvel, minvel
  real, optional,                         intent(in) :: str
  real, dimension(SZIB_(G),SZK_(G)), optional, intent(in) :: a, hv
! This subroutine writes to an output file all of the accelerations
! that have been applied to a column of zonal velocities over the
! previous timestep.  This subroutine is called from vertvisc.

! Arguments: I - The zonal index of the column to be documented.
!  (in)      j - The meridional index of the column to be documented.
!  (in)      um - The new zonal velocity, in m s-1.
!  (in)      hin - The layer thickness, in m.
!  (in)      ADp - A structure pointing to the various accelerations in
!                  the momentum equations.
!  (in)      CDp - A structure with pointers to various terms in the continuity
!                  equations.
!  (in)      dt - The model's dynamics time step.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 PointAccel_init.
!  (in)      str - The surface wind stress integrated over a time
!                  step, in m2 s-1.
!  (in)      a - The layer coupling coefficients from vertvisc, m.
!  (in)      hv - The layer thicknesses at velocity grid points, from
!                 vertvisc, in m.

  real    :: f_eff, CFL
  real    :: Angstrom
  real    :: truncvel, du
  real    :: Inorm(SZK_(G))
  real    :: e(SZK_(G)+1)
  integer :: yr, mo, day, hr, minute, sec, yearday
  integer :: k, ks, ke
  integer :: nz
  logical :: do_k(SZK_(G)+1)
  logical :: prev_avail
  integer :: file

  Angstrom = GV%Angstrom + GV%H_subroundoff

!  if (.not.associated(CS)) return
  nz = G%ke
  if (CS%cols_written < CS%max_writes) then
    CS%cols_written = CS%cols_written + 1

    ks = 1 ; ke = nz
    do_k(:) = .false.

  ! Open up the file for output if this is the first call.
    if (CS%u_file < 0) then
      if (len_trim(CS%u_trunc_file) < 1) return
      call open_file(CS%u_file, trim(CS%u_trunc_file), action=APPEND_FILE, &
                     form=ASCII_FILE, threading=MULTIPLE, fileset=SINGLE_FILE)
      if (CS%u_file < 0) then
        call MOM_error(NOTE, 'Unable to open file '//trim(CS%u_trunc_file)//'.')
        return
      endif
    endif
    file = CS%u_file

    prev_avail = (associated(CS%u_prev) .and. associated(CS%v_prev))

  ! Determine which layers to write out accelerations for.
    do k=1,nz
      if (((max(CS%u_av(I,j,k),um(I,j,k)) >= maxvel) .or. &
           (min(CS%u_av(I,j,k),um(I,j,k)) <= minvel)) .and. &
          ((hin(i,j,k) + hin(i+1,j,k)) > 3.0*Angstrom)) exit
    enddo
    ks = k
    do k=nz,1,-1
      if (((max(CS%u_av(I,j,k), um(I,j,k)) >= maxvel) .or. &
           (min(CS%u_av(I,j,k), um(I,j,k)) <= minvel)) .and. &
          ((hin(i,j,k) + hin(i+1,j,k)) > 3.0*Angstrom)) exit
    enddo
    ke = k
    if (ke < ks) then
      ks = 1; ke = nz; write(file,'("U: Unable to set ks & ke.")')
    endif

    call get_date(CS%Time, yr, mo, day, hr, minute, sec)
    call get_time((CS%Time - set_date(yr, 1, 1, 0, 0, 0)), sec, yearday)
    write (file,'(/,"--------------------------")')
    write (file,'(/,"Time ",i5,i4,F6.2," U-velocity violation at ",I4,": ",2(I3), &
        & " (",F7.2," E "F7.2," N) Layers ",I3," to ",I3,". dt = ",1PG10.4)') &
        yr, yearday, (REAL(sec)/3600.0), pe_here(), I, j, &
        G%geoLonCu(I,j), G%geoLatCu(I,j), ks, ke, dt

    if (ks <= GV%nk_rho_varies) ks = 1
    do k=ks,ke
      if ((hin(i,j,k) + hin(i+1,j,k)) > 3.0*Angstrom) do_k(k) = .true.
    enddo

    write(file,'(/,"Layers:",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(I10," ",$)') (k); enddo
    write(file,'(/,"u(m):  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (um(I,j,k)); enddo
    if (prev_avail) then
      write(file,'(/,"u(mp): ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (CS%u_prev(I,j,k)); enddo
    endif
    write(file,'(/,"u(3):  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (CS%u_av(I,j,k)); enddo

    write(file,'(/,"CFL u: ",$)')
    do k=ks,ke ; if (do_k(k)) then
      CFL = abs(um(I,j,k)) * dt * G%dy_Cu(I,j)
      if (um(I,j,k) < 0.0) then ; CFL = CFL * G%IareaT(i+1,j)
      else ; CFL = CFL * G%IareaT(i,j) ; endif
      write(file,'(ES10.3," ",$)') CFL
    endif ; enddo
    write(file,'(/,"CFL0 u:",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                    abs(um(I,j,k)) * dt * G%IdxCu(I,j) ; enddo

    if (prev_avail) then
      write(file,'(/,"du:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      ((um(I,j,k)-CS%u_prev(I,j,k))); enddo
    endif
    write(file,'(/,"CAu:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (dt*ADp%CAu(I,j,k)); enddo
    write(file,'(/,"PFu:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (dt*ADp%PFu(I,j,k)); enddo
    write(file,'(/,"diffu: ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (dt*ADp%diffu(I,j,k)); enddo

    if (ASSOCIATED(ADp%gradKEu)) then
      write(file,'(/,"KEu:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      (dt*ADp%gradKEu(I,j,k)); enddo
    endif
    if (ASSOCIATED(ADp%rv_x_v)) then
      write(file,'(/,"Coru:  ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
          dt*(ADp%CAu(I,j,k)-ADp%rv_x_v(I,j,k)); enddo
    endif
    if (ASSOCIATED(ADp%du_dt_visc)) then
      write(file,'(/,"ubv:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
          (um(I,j,k)-dt*ADp%du_dt_visc(I,j,k)); enddo
      write(file,'(/,"duv:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      (dt*ADp%du_dt_visc(I,j,k)); enddo
    endif
    if (ASSOCIATED(ADp%du_other)) then
      write(file,'(/,"du_other: ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      (ADp%du_other(I,j,k)); enddo
    endif
    if (present(a)) then
      write(file,'(/,"a:     ",$)')
      do k=ks,ke+1 ; if (do_k(k)) write(file,'(ES10.3," ",$)') a(I,k); enddo
    endif
    if (present(hv)) then
      write(file,'(/,"hvel:  ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') hv(I,k); enddo
    endif
    write(file,'(/,"Stress:  ",ES10.3)') str

    if (ASSOCIATED(CS%u_accel_bt)) then
      write(file,'("dubt:  ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      (dt*CS%u_accel_bt(I,j,k)) ; enddo
      write(file,'(/)')
    endif

    write(file,'(/,"h--:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (hin(i,j-1,k)); enddo
    write(file,'(/,"h+-:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (hin(i+1,j-1,k)); enddo
    write(file,'(/,"h-0:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (hin(i,j,k)); enddo
    write(file,'(/,"h+0:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (hin(i+1,j,k)); enddo
    write(file,'(/,"h-+:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (hin(i,j+1,k)); enddo
    write(file,'(/,"h++:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (hin(i+1,j+1,k)); enddo


    e(nz+1) = -G%bathyT(i,j)
    do k=nz,1,-1 ; e(K) = e(K+1) + hin(i,j,k) ; enddo
    write(file,'(/,"e-:    ",$)')
    write(file,'(ES10.3," ",$)') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ",$)') e(K); enddo

    e(nz+1) = -G%bathyT(i+1,j)
    do k=nz,1,-1 ; e(K) = e(K+1) + hin(i+1,j,k) ; enddo
    write(file,'(/,"e+:    ",$)')
    write(file,'(ES10.3," ",$)') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ",$)') e(K) ; enddo
    if (ASSOCIATED(CS%T)) then
      write(file,'(/,"T-:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%T(i,j,k); enddo
      write(file,'(/,"T+:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%T(i+1,j,k); enddo
    endif
    if (ASSOCIATED(CS%S)) then
      write(file,'(/,"S-:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%S(i,j,k); enddo
      write(file,'(/,"S+:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%S(i+1,j,k); enddo
    endif

    if (prev_avail) then
      write(file,'(/,"v--:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (CS%v_prev(i,J-1,k)); enddo
      write(file,'(/,"v-+:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (CS%v_prev(i,J,k)); enddo
      write(file,'(/,"v+-:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (CS%v_prev(i+1,J-1,k)); enddo
      write(file,'(/,"v++:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (CS%v_prev(i+1,J,k)); enddo
    endif

    write(file,'(/,"vh--:  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                    (CDp%vh(i,J-1,k)*G%IdxCv(i,J-1)); enddo
    write(file,'(/," vhC--:",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                        (0.5*CS%v_av(i,j-1,k)*(hin(i,j-1,k) + hin(i,j,k))); enddo
    if (prev_avail) then
      write(file,'(/," vhCp--:",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                          (0.5*CS%v_prev(i,j-1,k)*(hin(i,j-1,k) + hin(i,j,k))); enddo
    endif

    write(file,'(/,"vh-+:  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                    (CDp%vh(i,J,k)*G%IdxCv(i,J)); enddo
    write(file,'(/," vhC-+:",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                        (0.5*CS%v_av(i,J,k)*(hin(i,j,k) + hin(i,j+1,k))); enddo
    if (prev_avail) then
      write(file,'(/," vhCp-+:",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                          (0.5*CS%v_prev(i,J,k)*(hin(i,j,k) + hin(i,j+1,k))); enddo
    endif

    write(file,'(/,"vh+-:  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      (CDp%vh(i+1,J-1,k)*G%IdxCv(i+1,J-1)); enddo
    write(file,'(/," vhC+-:",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                    (0.5*CS%v_av(i+1,J-1,k)*(hin(i+1,j-1,k) + hin(i+1,j,k))); enddo
    if (prev_avail) then
      write(file,'(/," vhCp+-:",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                      (0.5*CS%v_prev(i+1,J-1,k)*(hin(i+1,j-1,k) + hin(i+1,j,k))); enddo
    endif

    write(file,'(/,"vh++:  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                          (CDp%vh(i+1,J,k)*G%IdxCv(i+1,J)); enddo
    write(file,'(/," vhC++:",$)')
         do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                     (0.5*CS%v_av(i+1,J,k)*(hin(i+1,j,k) + hin(i+1,j+1,k))); enddo
    if (prev_avail) then
      write(file,'(/," vhCp++:",$)')
           do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                       (0.5*CS%v_av(i+1,J,k)*(hin(i+1,j,k) + hin(i+1,j+1,k))); enddo
    endif

    write(file,'(/,"D:     ",2(ES10.3))') G%bathyT(i,j),G%bathyT(i+1,j)

  !  From here on, the normalized accelerations are written.
    if (prev_avail) then
      do k=ks,ke
        du = um(I,j,k)-CS%u_prev(I,j,k)
        if (abs(du) < 1.0e-6) du = 1.0e-6
        Inorm(k) = 1.0 / du
      enddo

      write(file,'(2/,"Norm:  ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') (1.0/Inorm(k)); enddo

      write(file,'(/,"du:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                        ((um(I,j,k)-CS%u_prev(I,j,k))*Inorm(k)); enddo

      write(file,'(/,"CAu:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                      (dt*ADp%CAu(I,j,k)*Inorm(k)); enddo

      write(file,'(/,"PFu:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                      (dt*ADp%PFu(I,j,k)*Inorm(k)); enddo

      write(file,'(/,"diffu: ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                      (dt*ADp%diffu(I,j,k)*Inorm(k)); enddo

      if (ASSOCIATED(ADp%gradKEu)) then
        write(file,'(/,"KEu:   ",$)')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                        (dt*ADp%gradKEu(I,j,k)*Inorm(k)); enddo
      endif
      if (ASSOCIATED(ADp%rv_x_v)) then
        write(file,'(/,"Coru:  ",$)')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
            dt*(ADp%CAu(I,j,k)-ADp%rv_x_v(I,j,k))*Inorm(k); enddo
      endif
      if (ASSOCIATED(ADp%du_dt_visc)) then
        write(file,'(/,"duv:   ",$)')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
            (dt*ADp%du_dt_visc(I,j,k))*Inorm(k); enddo
      endif
      if (ASSOCIATED(ADp%du_other)) then
        write(file,'(/,"du_other: ",$)')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
            (ADp%du_other(I,j,k))*Inorm(k); enddo
      endif
      if (ASSOCIATED(CS%u_accel_bt)) then
        write(file,'(/,"dubt:  ",$)')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                        (dt*CS%u_accel_bt(I,j,k)*Inorm(k)) ; enddo
      endif
    endif

    write(file,'(2/)')

    call flush(file)
  endif

end subroutine write_u_accel


subroutine write_v_accel(i, J, vm, hin, ADp, CDp, dt, G, GV, CS, &
                         maxvel, minvel, str, a, hv)
  integer,                                intent(in) :: i, J
  type(ocean_grid_type),                  intent(in) :: G
  type(verticalGrid_type),                intent(in) :: GV
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in) :: vm
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in) :: hin
  type(accel_diag_ptrs),                  intent(in) :: ADp
  type(cont_diag_ptrs),                   intent(in) :: CDp
  real,                                   intent(in) :: dt
  type(PointAccel_CS),                    pointer    :: CS
  real,                                   intent(in) :: maxvel, minvel
  real, optional,                         intent(in) :: str
  real, dimension(SZI_(G),SZK_(G)), optional, intent(in) :: a, hv

! This subroutine writes to an output file all of the accelerations
! that have been applied to a column of meridional velocities over
! the previous timestep.  This subroutine is called from vertvisc.

! Arguments: i - The zonal index of the column to be documented.
!  (in)      J - The meridional index of the column to be documented.
!  (in)      vm - The new meridional velocity, in m s-1.
!  (in)      hin - The layer thickness, in m.
!  (in)      ADp - A structure pointing to the various accelerations in
!                  the momentum equations.
!  (in)      CDp - A structure with pointers to various terms in the continuity
!                  equations.
!  (in)      dt - The model's dynamics time step.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 PointAccel_init.
!  (in)      str - The surface wind stress integrated over a time
!                  step, in m2 s-1.
!  (in)      a - The layer coupling coefficients from vertvisc, m.
!  (in)      hv - The layer thicknesses at velocity grid points, from
!                 vertvisc, in m.

  real    :: f_eff, CFL
  real    :: Angstrom
  real    :: truncvel, dv
  real    :: Inorm(SZK_(G))
  real    :: e(SZK_(G)+1)
  integer :: yr, mo, day, hr, minute, sec, yearday
  integer :: k, ks, ke
  integer :: nz
  logical :: do_k(SZK_(G)+1)
  logical :: prev_avail
  integer :: file

  Angstrom = GV%Angstrom + GV%H_subroundoff

!  if (.not.associated(CS)) return
  nz = G%ke
  if (CS%cols_written < CS%max_writes) then
    CS%cols_written = CS%cols_written + 1

    ks = 1 ; ke = nz
    do_k(:) = .false.

  ! Open up the file for output if this is the first call.
    if (CS%v_file < 0) then
      if (len_trim(CS%v_trunc_file) < 1) return
      call open_file(CS%v_file, trim(CS%v_trunc_file), action=APPEND_FILE, &
                     form=ASCII_FILE, threading=MULTIPLE, fileset=SINGLE_FILE)
      if (CS%v_file < 0) then
        call MOM_error(NOTE, 'Unable to open file '//trim(CS%v_trunc_file)//'.')
        return
      endif
    endif
    file = CS%v_file

    prev_avail = (associated(CS%u_prev) .and. associated(CS%v_prev))

    do k=1,nz
      if (((max(CS%v_av(i,J,k), vm(i,J,k)) >= maxvel) .or. &
           (min(CS%v_av(i,J,k), vm(i,J,k)) <= minvel)) .and. &
          ((hin(i,j,k) + hin(i,j+1,k)) > 3.0*Angstrom)) exit
    enddo
    ks = k
    do k=nz,1,-1
      if (((max(CS%v_av(i,J,k), vm(i,J,k)) >= maxvel) .or. &
           (min(CS%v_av(i,J,k), vm(i,J,k)) <= minvel)) .and. &
          ((hin(i,j,k) + hin(i,j+1,k)) > 3.0*Angstrom)) exit
    enddo
    ke = k
    if (ke < ks) then
      ks = 1; ke = nz; write(file,'("V: Unable to set ks & ke.")')
    endif

    call get_date(CS%Time, yr, mo, day, hr, minute, sec)
    call get_time((CS%Time - set_date(yr, 1, 1, 0, 0, 0)), sec, yearday)
    write (file,'(/,"--------------------------")')
    write (file,'(/,"Time ",i5,i4,F6.2," V-velocity violation at ",I4,": ",2(I3), &
        & " (",F7.2," E ",F7.2," N) Layers ",I3," to ",I3,". dt = ",1PG10.4)') &
        yr, yearday, (REAL(sec)/3600.0), pe_here(), i, J, &
        G%geoLonCv(i,J), G%geoLatCv(i,J), ks, ke, dt

    if (ks <= GV%nk_rho_varies) ks = 1
    do k=ks,ke
      if ((hin(i,j,k) + hin(i,j+1,k)) > 3.0*Angstrom) do_k(k) = .true.
    enddo

    write(file,'(/,"Layers:",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(I10," ",$)') (k); enddo
    write(file,'(/,"v(m):  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (vm(i,J,k)); enddo

    if (prev_avail) then
      write(file,'(/,"v(mp): ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (CS%v_prev(i,J,k)); enddo
    endif

    write(file,'(/,"v(3):  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (CS%v_av(i,J,k)); enddo
    write(file,'(/,"CFL v: ",$)')
    do k=ks,ke ; if (do_k(k)) then
      CFL = abs(vm(i,J,k)) * dt * G%dx_Cv(i,J)
      if (vm(i,J,k) < 0.0) then ; CFL = CFL * G%IareaT(i,j+1)
      else ; CFL = CFL * G%IareaT(i,j) ; endif
      write(file,'(ES10.3," ",$)') CFL
    endif ; enddo
    write(file,'(/,"CFL0 v:",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                    abs(vm(i,J,k)) * dt * G%IdyCv(i,J) ; enddo

    if (prev_avail) then
      write(file,'(/,"dv:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      ((vm(i,J,k)-CS%v_prev(i,J,k))); enddo
    endif

    write(file,'(/,"CAv:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (dt*ADp%CAv(i,J,k)); enddo

    write(file,'(/,"PFv:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (dt*ADp%PFv(i,J,k)); enddo

    write(file,'(/,"diffv: ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') (dt*ADp%diffv(i,J,k)); enddo

    if (ASSOCIATED(ADp%gradKEv)) then
      write(file,'(/,"KEv:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      (dt*ADp%gradKEv(i,J,k)); enddo
    endif
    if (ASSOCIATED(ADp%rv_x_u)) then
      write(file,'(/,"Corv:  ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                 dt*(ADp%CAv(i,J,k)-ADp%rv_x_u(i,J,k)); enddo
    endif
    if (ASSOCIATED(ADp%dv_dt_visc)) then
      write(file,'(/,"vbv:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
          (vm(i,J,k)-dt*ADp%dv_dt_visc(i,J,k)); enddo

      write(file,'(/,"dvv:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      (dt*ADp%dv_dt_visc(i,J,k)); enddo
    endif
    if (ASSOCIATED(ADp%dv_other)) then
      write(file,'(/,"dv_other: ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      (ADp%dv_other(i,J,k)); enddo
    endif
    if (present(a)) then
      write(file,'(/,"a:     ",$)')
      do k=ks,ke+1 ; if (do_k(k)) write(file,'(ES10.3," ",$)') a(i,k); enddo
    endif
    if (present(hv)) then
      write(file,'(/,"hvel:  ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') hv(i,k); enddo
    endif
    write(file,'(/,"Stress:  ",ES10.3)') str

    if (ASSOCIATED(CS%v_accel_bt)) then
      write(file,'("dvbt:  ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                      (dt*CS%v_accel_bt(i,J,k)) ; enddo
      write(file,'(/)')
    endif

    write(file,'("h--:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') hin(i-1,j,k); enddo
    write(file,'(/,"h0-:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') hin(i,j,k); enddo
    write(file,'(/,"h+-:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') hin(i+1,j,k); enddo
    write(file,'(/,"h-+:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') hin(i-1,j+1,k); enddo
    write(file,'(/,"h0+:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') hin(i,j+1,k); enddo
    write(file,'(/,"h++:   ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') hin(i+1,j+1,k); enddo

    e(nz+1) = -G%bathyT(i,j)
    do k=nz,1,-1 ; e(K) = e(K+1) + hin(i,j,k); enddo
    write(file,'(/,"e-:    ",$)')
    write(file,'(ES10.3," ",$)') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ",$)') e(K); enddo

    e(nz+1) = -G%bathyT(i,j+1)
    do k=nz,1,-1 ; e(K) = e(K+1) + hin(i,j+1,k) ; enddo
    write(file,'(/,"e+:    ",$)')
    write(file,'(ES10.3," ",$)') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ",$)') e(K); enddo
    if (ASSOCIATED(CS%T)) then
      write(file,'(/,"T-:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%T(i,j,k); enddo
      write(file,'(/,"T+:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%T(i,j+1,k); enddo
    endif
    if (ASSOCIATED(CS%S)) then
      write(file,'(/,"S-:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%S(i,j,k); enddo
      write(file,'(/,"S+:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%S(i,j+1,k); enddo
    endif

    if (prev_avail) then
      write(file,'(/,"u--:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%u_prev(I-1,j,k); enddo
      write(file,'(/,"u-+:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%u_prev(I-1,j+1,k); enddo
      write(file,'(/,"u+-:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%u_prev(I,j,k); enddo
      write(file,'(/,"u++:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') CS%u_prev(I,j+1,k); enddo
    endif

    write(file,'(/,"uh--:  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                    (CDp%uh(I-1,j,k)*G%IdyCu(I-1,j)); enddo
    write(file,'(/," uhC--: ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
            (CS%u_av(I-1,j,k) * 0.5*(hin(i-1,j,k) + hin(i,j,k))); enddo
    if (prev_avail) then
      write(file,'(/," uhCp--:",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                          (0.5*CS%u_prev(I-1,j,k)*(hin(i-1,j,k) + hin(i,j,k))); enddo
    endif

    write(file,'(/,"uh-+:  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                    (CDp%uh(I-1,j+1,k)*G%IdyCu(I-1,j+1)); enddo
    write(file,'(/," uhC-+: ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
            (CS%u_av(I-1,j+1,k) * 0.5*(hin(i-1,j+1,k) + hin(i,j+1,k))); enddo
    if (prev_avail) then
      write(file,'(/," uhCp-+:",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                    (0.5*CS%u_prev(I-1,j+1,k)*(hin(i-1,j+1,k) + hin(i,j+1,k))); enddo
    endif

    write(file,'(/,"uh+-:  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                    (CDp%uh(I,j,k)*G%IdyCu(I,j)); enddo
    write(file,'(/," uhC+-: ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
            (CS%u_av(I,j,k) * 0.5*(hin(i,j,k) + hin(i+1,j,k))); enddo
    if (prev_avail) then
      write(file,'(/," uhCp+-:",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                            (0.5*CS%u_prev(I,j,k)*(hin(i,j,k) + hin(i+1,j,k))); enddo
    endif

    write(file,'(/,"uh++:  ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                                    (CDp%uh(I,j+1,k)*G%IdyCu(I,j+1)); enddo
    write(file,'(/," uhC++: ",$)')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
            (CS%u_av(I,j+1,k) * 0.5*(hin(i,j+1,k) + hin(i+1,j+1,k))); enddo
    if (prev_avail) then
      write(file,'(/," uhCp++:",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ",$)') &
                      (0.5*CS%u_prev(I,j+1,k)*(hin(i,j+1,k) + hin(i+1,j+1,k))); enddo
    endif

    write(file,'(/,"D:     ",2(ES10.3))') G%bathyT(i,j),G%bathyT(i,j+1)

  !  From here on, the normalized accelerations are written.
    if (prev_avail) then
      do k=ks,ke
        dv = vm(i,J,k)-CS%v_prev(i,J,k)
        if (abs(dv) < 1.0e-6) dv = 1.0e-6
        Inorm(k) = 1.0 / dv
      enddo

      write(file,'(2/,"Norm:  ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') (1.0/Inorm(k)); enddo
      write(file,'(/,"dv:    ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                      ((vm(i,J,k)-CS%v_prev(i,J,k))*Inorm(k)); enddo
      write(file,'(/,"CAv:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                      (dt*ADp%CAv(i,J,k)*Inorm(k)); enddo
      write(file,'(/,"PFv:   ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                      (dt*ADp%PFv(i,J,k)*Inorm(k)); enddo
      write(file,'(/,"diffv: ",$)')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                      (dt*ADp%diffv(i,J,k)*Inorm(k)); enddo

      if (ASSOCIATED(ADp%gradKEu)) then
        write(file,'(/,"KEv:   ",$)')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                      (dt*ADp%gradKEv(i,J,k)*Inorm(k)); enddo
      endif
      if (ASSOCIATED(ADp%rv_x_u)) then
        write(file,'(/,"Corv:  ",$)')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
            dt*(ADp%CAv(i,J,k)-ADp%rv_x_u(i,J,k))*Inorm(k); enddo
      endif
      if (ASSOCIATED(ADp%dv_dt_visc)) then
        write(file,'(/,"dvv:   ",$)')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
            (dt*ADp%dv_dt_visc(i,J,k)*Inorm(k)); enddo
      endif
      if (ASSOCIATED(ADp%dv_other)) then
        write(file,'(/,"dv_other: ",$)')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
            (ADp%dv_other(i,J,k)*Inorm(k)); enddo
      endif
      if (ASSOCIATED(CS%v_accel_bt)) then
        write(file,'(/,"dvbt:  ",$)')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ",$)') &
                                        (dt*CS%v_accel_bt(i,J,k)*Inorm(k)) ; enddo
      endif
    endif

    write(file,'(2/)')

    call flush(file)
  endif

end subroutine write_v_accel

subroutine PointAccel_init(MIS, Time, G, param_file, diag, dirs, CS)
  type(ocean_internal_state), target, intent(in) :: MIS
  type(time_type), target, intent(in) :: Time
  type(ocean_grid_type),   intent(in) :: G
  type(param_file_type),   intent(in) :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(directories),       intent(in) :: dirs
  type(PointAccel_CS),     pointer    :: CS
! Arguments: MIS - For "MOM Internal State" a set of pointers to the fields and
!                  accelerations that make up the ocean's physical state.
!  (in)      Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in)      dirs - A structure containing several relevant directory paths.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_PointAccel" ! This module's name.

  if (associated(CS)) return
  allocate(CS)

  CS%diag => diag ; CS%Time => Time

  CS%T => MIS%T ; CS%S => MIS%S ; CS%pbce => MIS%pbce
  CS%u_accel_bt => MIS%u_accel_bt ; CS%v_accel_bt => MIS%v_accel_bt
  CS%u_prev => MIS%u_prev ; CS%v_prev => MIS%v_prev
  CS%u_av => MIS%u_av; if (.not.associated(MIS%u_av)) CS%u_av => MIS%u(:,:,:)
  CS%v_av => MIS%v_av; if (.not.associated(MIS%v_av)) CS%v_av => MIS%v(:,:,:)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "U_TRUNC_FILE", CS%u_trunc_file, &
                 "The absolute path to the file where the accelerations \n"//&
                 "leading to zonal velocity truncations are written. \n"//&
                 "Leave this empty for efficiency if this diagnostic is \n"//&
                 "not needed.", default="")
  call get_param(param_file, mod, "V_TRUNC_FILE", CS%v_trunc_file, &
                 "The absolute path to the file where the accelerations \n"//&
                 "leading to meridional velocity truncations are written. \n"//&
                 "Leave this empty for efficiency if this diagnostic is \n"//&
                 "not needed.", default="")
  call get_param(param_file, mod, "MAX_TRUNC_FILE_SIZE_PER_PE", CS%max_writes, &
                 "The maximum number of colums of truncations that any PE \n"//&
                 "will write out during a run.", default=50)

  if (len_trim(dirs%output_directory) > 0) then
    if (len_trim(CS%u_trunc_file) > 0) &
      CS%u_trunc_file = trim(dirs%output_directory)//trim(CS%u_trunc_file)
    if (len_trim(CS%v_trunc_file) > 0) &
      CS%v_trunc_file = trim(dirs%output_directory)//trim(CS%v_trunc_file)
    call log_param(param_file, mod, "output_dir/U_TRUNC_FILE", CS%u_trunc_file)
    call log_param(param_file, mod, "output_dir/V_TRUNC_FILE", CS%v_trunc_file)
  endif
  CS%u_file = -1 ; CS%v_file = -1 ; CS%cols_written = 0

end subroutine PointAccel_init
end module MOM_PointAccel
