!> Debug accelerations at a given point
!!
!!    The two subroutines in this file write out all of the terms
!! in the u- or v-momentum balance at a given point.  Usually
!! these subroutines are called after the velocities exceed some
!! threshold, in order to determine which term is culpable.
!! often this is done for debugging purposes.
module MOM_PointAccel

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : diag_ctrl
use MOM_domains, only : pe_here
use MOM_error_handler, only : MOM_error, NOTE
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : open_ASCII_file, APPEND_FILE, MULTIPLE, SINGLE_FILE
use MOM_time_manager, only : time_type, get_time, get_date, set_date, operator(-)
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : ocean_internal_state, accel_diag_ptrs, cont_diag_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public write_u_accel, write_v_accel, PointAccel_init

!> The control structure for the MOM_PointAccel module
type, public :: PointAccel_CS ; private
  character(len=200) :: u_trunc_file !< The complete path to the file in which a column's worth of
                                     !! u-accelerations are written if u-velocity truncations occur.
  character(len=200) :: v_trunc_file !< The complete path to the file in which a column's worth of
                                     !! v-accelerations are written if v-velocity truncations occur.
  integer :: u_file         !< The unit number for an opened u-truncation files, or -1 if it has not yet been opened.
  integer :: v_file         !< The unit number for an opened v-truncation files, or -1 if it has not yet been opened.
  integer :: cols_written   !< The number of columns whose output has been
                            !! written by this PE during the current run.
  integer :: max_writes     !< The maximum number of times any PE can write out
                            !! a column's worth of accelerations during a run.
  logical :: full_column    !< If true, write out the accelerations in all massive layers,
                            !! otherwise just document the ones with large velocities.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
! The following are pointers to many of the state variables and accelerations
! that are used to step the physical model forward.  They all use the same
! names as the variables they point to in MOM.F90
  real, pointer, dimension(:,:,:) :: &
    u_av => NULL(), &       !< Time average u-velocity [L T-1 ~> m s-1]
    v_av => NULL(), &       !< Time average velocity [L T-1 ~> m s-1]
    u_prev => NULL(), &     !< Previous u-velocity [L T-1 ~> m s-1]
    v_prev => NULL(), &     !< Previous v-velocity [L T-1 ~> m s-1]
    T => NULL(), &          !< Temperature [C ~> degC]
    S => NULL(), &          !< Salinity [S ~> ppt]
    u_accel_bt => NULL(), & !< Barotropic u-accelerations [L T-2 ~> m s-2]
    v_accel_bt => NULL()    !< Barotropic v-accelerations [L T-2 ~> m s-2]
end type PointAccel_CS

contains

!> This subroutine writes to an output file all of the accelerations
!! that have been applied to a column of zonal velocities over the
!! previous timestep.  This subroutine is called from vertvisc.
subroutine write_u_accel(I, j, um, hin, ADp, CDp, dt, G, GV, US, CS, vel_rpt, str, a, hv)
  integer,                     intent(in) :: I   !< The zonal index of the column to be documented.
  integer,                     intent(in) :: j   !< The meridional index of the column to be documented.
  type(ocean_grid_type),       intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),     intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),       intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                               intent(in) :: um  !< The new zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                               intent(in) :: hin !< The layer thickness [H ~> m or kg m-2].
  type(accel_diag_ptrs),       intent(in) :: ADp !< A structure pointing to the various
                                                 !! accelerations in the momentum equations.
  type(cont_diag_ptrs),        intent(in) :: CDp !<  A structure with pointers to various terms
                                                 !! in the continuity equations.
  real,                        intent(in) :: dt  !< The ocean dynamics time step [T ~> s].
  type(PointAccel_CS),         pointer    :: CS  !< The control structure returned by a previous
                                                 !! call to PointAccel_init.
  real,                        intent(in) :: vel_rpt !< The velocity magnitude that triggers a report [L T-1 ~> m s-1]
  real, optional,              intent(in) :: str !< The surface wind stress [R L Z T-2 ~> Pa]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1), &
                     optional, intent(in) :: a   !< The layer coupling coefficients from vertvisc
                                                 !! [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                     optional, intent(in) :: hv  !< The layer thicknesses at velocity grid points,
                                                 !! from vertvisc [H ~> m or kg m-2].

  ! Local variables
  real    :: CFL              ! The local velocity-based CFL number [nondim]
  real    :: Angstrom         ! A negligibly small thickness [H ~> m or kg m-2]
  real    :: du               ! A velocity change [L T-1 ~> m s-1]
  real    :: Inorm(SZK_(GV))  ! The inverse of the normalized velocity change [L T-1 ~> m s-1]
  real    :: e(SZK_(GV)+1)    ! Simple estimates of interface heights based on the sum of thicknesses [m]
  real    :: h_scale          ! A scaling factor for thicknesses [m H-1 ~> 1 or m3 kg-1]
  real    :: vel_scale        ! A scaling factor for velocities [m T s-1 L-1 ~> 1]
  real    :: uh_scale         ! A scaling factor for transport per unit length [m2 T s-1 L-1 H-1 ~> 1 or m3 kg-1]
  real    :: temp_scale       ! A scaling factor for temperatures [degC C-1 ~> 1]
  real    :: saln_scale       ! A scaling factor for salinities [ppt S-1 ~> 1]
  integer :: yr, mo, day, hr, minute, sec, yearday
  integer :: k, ks, ke
  integer :: nz
  logical :: do_k(SZK_(GV)+1)
  logical :: prev_avail
  integer :: file

  Angstrom = GV%Angstrom_H + GV%H_subroundoff
  h_scale = GV%H_to_m ; vel_scale = US%L_T_to_m_s ; uh_scale = GV%H_to_m*US%L_T_to_m_s
  temp_scale = US%C_to_degC ; saln_scale = US%S_to_ppt

!  if (.not.associated(CS)) return
  nz = GV%ke
  if (CS%cols_written < CS%max_writes) then
    CS%cols_written = CS%cols_written + 1

    ks = 1 ; ke = nz
    do_k(:) = .false.

  ! Open up the file for output if this is the first call.
    if (CS%u_file < 0) then
      if (len_trim(CS%u_trunc_file) < 1) return
      call open_ASCII_file(CS%u_file, trim(CS%u_trunc_file), action=APPEND_FILE, &
                           threading=MULTIPLE, fileset=SINGLE_FILE)
      if (CS%u_file < 0) then
        call MOM_error(NOTE, 'Unable to open file '//trim(CS%u_trunc_file)//'.')
        return
      endif
    endif
    file = CS%u_file

    prev_avail = (associated(CS%u_prev) .and. associated(CS%v_prev))

  ! Determine which layers to write out accelerations for.
    do k=1,nz
      if (((max(CS%u_av(I,j,k),um(I,j,k)) >= vel_rpt) .or. &
           (min(CS%u_av(I,j,k),um(I,j,k)) <= -vel_rpt)) .and. &
          ((hin(i,j,k) + hin(i+1,j,k)) > 3.0*Angstrom)) exit
    enddo
    ks = k
    do k=nz,1,-1
      if (((max(CS%u_av(I,j,k), um(I,j,k)) >= vel_rpt) .or. &
           (min(CS%u_av(I,j,k), um(I,j,k)) <= -vel_rpt)) .and. &
          ((hin(i,j,k) + hin(i+1,j,k)) > 3.0*Angstrom)) exit
    enddo
    ke = k
    if (ke < ks) then
      ks = 1; ke = nz; write(file,'("U: Unable to set ks & ke.")')
    endif
    if (CS%full_column) then
      ks = 1 ; ke = nz
    endif

    call get_date(CS%Time, yr, mo, day, hr, minute, sec)
    call get_time((CS%Time - set_date(yr, 1, 1, 0, 0, 0)), sec, yearday)
    write (file,'(/,"--------------------------")')
    write (file,'(/,"Time ",i5,i4,F6.2," U-velocity violation at ",I4,": ",2(I3), &
        & " (",F7.2," E ",F7.2," N) Layers ",I3," to ",I3,". dt = ",1PG10.4)') &
        yr, yearday, (REAL(sec)/3600.0), pe_here(), I, j, &
        G%geoLonCu(I,j), G%geoLatCu(I,j), ks, ke, US%T_to_s*dt

    if (ks <= GV%nk_rho_varies) ks = 1
    do k=ks,ke
      if ((hin(i,j,k) + hin(i+1,j,k)) > 3.0*Angstrom) do_k(k) = .true.
    enddo

    write(file,'(/,"Layers:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(I10," ")', advance='no') (k) ; enddo
    write(file,'(/,"u(m):  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*um(I,j,k)) ; enddo
    if (prev_avail) then
      write(file,'(/,"u(mp): ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%u_prev(I,j,k)) ; enddo
    endif
    write(file,'(/,"u(3):  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%u_av(I,j,k)) ; enddo

    write(file,'(/,"CFL u: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) then
      CFL = abs(um(I,j,k)) * dt * G%dy_Cu(I,j)
      if (um(I,j,k) < 0.0) then ; CFL = CFL * G%IareaT(i+1,j)
      else ; CFL = CFL * G%IareaT(i,j) ; endif
      write(file,'(ES10.3," ")', advance='no') CFL
    endif ; enddo
    write(file,'(/,"CFL0 u:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    abs(um(I,j,k)) * dt * G%IdxCu(I,j) ; enddo

    if (prev_avail) then
      write(file,'(/,"du:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*(um(I,j,k)-CS%u_prev(I,j,k))) ; enddo
    endif
    write(file,'(/,"CAu:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%CAu(I,j,k)) ; enddo
    write(file,'(/,"PFu:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%PFu(I,j,k)) ; enddo
    write(file,'(/,"diffu: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%diffu(I,j,k)) ; enddo

    if (associated(ADp%gradKEu)) then
      write(file,'(/,"KEu:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*ADp%gradKEu(I,j,k)) ; enddo
    endif
    if (associated(ADp%rv_x_v)) then
      write(file,'(/,"Coru:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
          vel_scale*dt*(ADp%CAu(I,j,k)-ADp%rv_x_v(I,j,k)) ; enddo
    endif
    if (associated(ADp%du_dt_visc)) then
      write(file,'(/,"ubv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
          vel_scale*(um(I,j,k) - dt*ADp%du_dt_visc(I,j,k)) ; enddo
      write(file,'(/,"duv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*ADp%du_dt_visc(I,j,k)) ; enddo
    endif
    if (associated(ADp%du_other)) then
      write(file,'(/,"du_other: ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*ADp%du_other(I,j,k)) ; enddo
    endif
    if (present(a)) then
      write(file,'(/,"a:     ",ES10.3," ")', advance='no') h_scale*a(I,j,ks)*dt
      do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') (h_scale*a(I,j,K)*dt) ; enddo
    endif
    if (present(hv)) then
      write(file,'(/,"hvel:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hv(I,j,k) ; enddo
    endif
    if (present(str)) then
      write(file,'(/,"Stress:  ",ES10.3)', advance='no') vel_scale*US%Z_to_m * (str*dt / GV%Rho0)
    endif

    if (associated(CS%u_accel_bt)) then
      write(file,'(/,"dubt:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*CS%u_accel_bt(I,j,k)) ; enddo
    endif
    write(file,'(/)')

    write(file,'(/,"h--:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i,j-1,k)) ; enddo
    write(file,'(/,"h+-:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i+1,j-1,k)) ; enddo
    write(file,'(/,"h-0:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i,j,k)) ; enddo
    write(file,'(/,"h+0:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i+1,j,k)) ; enddo
    write(file,'(/,"h-+:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i,j+1,k)) ; enddo
    write(file,'(/,"h++:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (h_scale*hin(i+1,j+1,k)) ; enddo


    e(nz+1) = -US%Z_to_m*(G%bathyT(i,j) + G%Z_ref)
    do k=nz,1,-1 ; e(K) = e(K+1) + h_scale*hin(i,j,k) ; enddo
    write(file,'(/,"e-:    ",ES10.3," ")', advance='no') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') e(K) ; enddo

    e(nz+1) = -US%Z_to_m*(G%bathyT(i+1,j) + G%Z_ref)
    do k=nz,1,-1 ; e(K) = e(K+1) + h_scale*hin(i+1,j,k) ; enddo
    write(file,'(/,"e+:    ",ES10.3," ")', advance='no') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') e(K) ; enddo
    if (associated(CS%T)) then
      write(file,'(/,"T-:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') temp_scale*CS%T(i,j,k) ; enddo
      write(file,'(/,"T+:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') temp_scale*CS%T(i+1,j,k) ; enddo
    endif
    if (associated(CS%S)) then
      write(file,'(/,"S-:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') saln_scale*CS%S(i,j,k) ; enddo
      write(file,'(/,"S+:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') saln_scale*CS%S(i+1,j,k) ; enddo
    endif

    if (prev_avail) then
      write(file,'(/,"v--:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_prev(i,J-1,k)) ; enddo
      write(file,'(/,"v-+:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_prev(i,J,k)) ; enddo
      write(file,'(/,"v+-:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_prev(i+1,J-1,k)) ; enddo
      write(file,'(/,"v++:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_prev(i+1,J,k)) ; enddo
    endif

    write(file,'(/,"vh--:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%vh(i,J-1,k)*G%IdxCv(i,J-1)) ; enddo
    write(file,'(/," vhC--:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                        (0.5*CS%v_av(i,j-1,k)*uh_scale*(hin(i,j-1,k) + hin(i,j,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," vhCp--:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                          (0.5*CS%v_prev(i,j-1,k)*uh_scale*(hin(i,j-1,k) + hin(i,j,k))) ; enddo
    endif

    write(file,'(/,"vh-+:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%vh(i,J,k)*G%IdxCv(i,J)) ; enddo
    write(file,'(/," vhC-+:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                        (0.5*CS%v_av(i,J,k)*uh_scale*(hin(i,j,k) + hin(i,j+1,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," vhCp-+:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                          (0.5*CS%v_prev(i,J,k)*uh_scale*(hin(i,j,k) + hin(i,j+1,k))) ; enddo
    endif

    write(file,'(/,"vh+-:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (uh_scale*CDp%vh(i+1,J-1,k)*G%IdxCv(i+1,J-1)) ; enddo
    write(file,'(/," vhC+-:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                    (0.5*CS%v_av(i+1,J-1,k)*uh_scale*(hin(i+1,j-1,k) + hin(i+1,j,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," vhCp+-:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                      (0.5*CS%v_prev(i+1,J-1,k)*uh_scale*(hin(i+1,j-1,k) + hin(i+1,j,k))) ; enddo
    endif

    write(file,'(/,"vh++:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                          (uh_scale*CDp%vh(i+1,J,k)*G%IdxCv(i+1,J)) ; enddo
    write(file,'(/," vhC++:")', advance='no')
         do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                     (0.5*CS%v_av(i+1,J,k)*uh_scale*(hin(i+1,j,k) + hin(i+1,j+1,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," vhCp++:")', advance='no')
           do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                       (0.5*CS%v_av(i+1,J,k)*uh_scale*(hin(i+1,j,k) + hin(i+1,j+1,k))) ; enddo
    endif

    write(file,'(/,"D:     ",2(ES10.3))') US%Z_to_m*(G%bathyT(i,j) + G%Z_ref), US%Z_to_m*(G%bathyT(i+1,j) + G%Z_ref)

  !  From here on, the normalized accelerations are written.
    if (prev_avail) then
      do k=ks,ke
        du = um(I,j,k) - CS%u_prev(I,j,k)
        if (abs(du) < 1.0e-6*US%m_s_to_L_T) du = 1.0e-6*US%m_s_to_L_T
        Inorm(k) = 1.0 / du
      enddo

      write(file,'(2/,"Norm:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') (vel_scale / Inorm(k)) ; enddo

      write(file,'(/,"du:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      ((um(I,j,k)-CS%u_prev(I,j,k)) * Inorm(k)) ; enddo

      write(file,'(/,"CAu:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%CAu(I,j,k) * Inorm(k)) ; enddo

      write(file,'(/,"PFu:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%PFu(I,j,k) * Inorm(k)) ; enddo

      write(file,'(/,"diffu: ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%diffu(I,j,k) * Inorm(k)) ; enddo

      if (associated(ADp%gradKEu)) then
        write(file,'(/,"KEu:   ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*ADp%gradKEu(I,j,k) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%rv_x_v)) then
        write(file,'(/,"Coru:  ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*(ADp%CAu(I,j,k)-ADp%rv_x_v(I,j,k)) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%du_dt_visc)) then
        write(file,'(/,"duv:   ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*ADp%du_dt_visc(I,j,k) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%du_other)) then
        write(file,'(/,"du_other: ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (ADp%du_other(I,j,k) * Inorm(k)) ; enddo
      endif
      if (associated(CS%u_accel_bt)) then
        write(file,'(/,"dubt:  ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*CS%u_accel_bt(I,j,k) * Inorm(k)) ; enddo
      endif
    endif

    write(file,'(2/)')

    flush(file)
  endif

end subroutine write_u_accel

!> This subroutine writes to an output file all of the accelerations
!! that have been applied to a column of meridional velocities over
!! the previous timestep.  This subroutine is called from vertvisc.
subroutine write_v_accel(i, J, vm, hin, ADp, CDp, dt, G, GV, US, CS, vel_rpt, str, a, hv)
  integer,                     intent(in) :: i   !< The zonal index of the column to be documented.
  integer,                     intent(in) :: J   !< The meridional index of the column to be documented.
  type(ocean_grid_type),       intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),     intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),       intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                               intent(in) :: vm  !< The new meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                               intent(in) :: hin !< The layer thickness [H ~> m or kg m-2].
  type(accel_diag_ptrs),       intent(in) :: ADp !< A structure pointing to the various
                                                 !! accelerations in the momentum equations.
  type(cont_diag_ptrs),        intent(in) :: CDp !< A structure with pointers to various terms in
                                                 !! the continuity equations.
  real,                        intent(in) :: dt  !< The ocean dynamics time step [T ~> s].
  type(PointAccel_CS),         pointer    :: CS  !< The control structure returned by a previous
                                                 !! call to PointAccel_init.
  real,                        intent(in) :: vel_rpt !< The velocity magnitude that triggers a report [L T-1 ~> m s-1]
  real, optional,              intent(in) :: str !< The surface wind stress [R L Z T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1), &
                     optional, intent(in) :: a   !< The layer coupling coefficients from vertvisc
                                                 !! [H T-1 ~> m s-1 or Pa s m-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                     optional, intent(in) :: hv  !< The layer thicknesses at velocity grid points,
                                                 !! from vertvisc [H ~> m or kg m-2].

  ! Local variables
  real    :: CFL              ! The local velocity-based CFL number [nondim]
  real    :: Angstrom         ! A negligibly small thickness [H ~> m or kg m-2]
  real    :: dv               ! A velocity change [L T-1 ~> m s-1]
  real    :: Inorm(SZK_(GV))  ! The inverse of the normalized velocity change [L T-1 ~> m s-1]
  real    :: e(SZK_(GV)+1)    ! Simple estimates of interface heights based on the sum of thicknesses [m]
  real    :: h_scale          ! A scaling factor for thicknesses [m H-1 ~> 1 or m3 kg-1]
  real    :: vel_scale        ! A scaling factor for velocities [m T s-1 L-1 ~> 1]
  real    :: uh_scale         ! A scaling factor for transport per unit length [m2 T s-1 L-1 H-1 ~> 1 or m3 kg-1]
  real    :: temp_scale       ! A scaling factor for temperatures [degC C-1 ~> 1]
  real    :: saln_scale       ! A scaling factor for salinities [ppt S-1 ~> 1]
  integer :: yr, mo, day, hr, minute, sec, yearday
  integer :: k, ks, ke
  integer :: nz
  logical :: do_k(SZK_(GV)+1)
  logical :: prev_avail
  integer :: file

  Angstrom = GV%Angstrom_H + GV%H_subroundoff
  h_scale = GV%H_to_m ; vel_scale = US%L_T_to_m_s ; uh_scale = GV%H_to_m*US%L_T_to_m_s
  temp_scale = US%C_to_degC ; saln_scale = US%S_to_ppt

!  if (.not.associated(CS)) return
  nz = GV%ke
  if (CS%cols_written < CS%max_writes) then
    CS%cols_written = CS%cols_written + 1

    ks = 1 ; ke = nz
    do_k(:) = .false.

  ! Open up the file for output if this is the first call.
    if (CS%v_file < 0) then
      if (len_trim(CS%v_trunc_file) < 1) return
      call open_ASCII_file(CS%v_file, trim(CS%v_trunc_file), action=APPEND_FILE, &
                           threading=MULTIPLE, fileset=SINGLE_FILE)
      if (CS%v_file < 0) then
        call MOM_error(NOTE, 'Unable to open file '//trim(CS%v_trunc_file)//'.')
        return
      endif
    endif
    file = CS%v_file

    prev_avail = (associated(CS%u_prev) .and. associated(CS%v_prev))

    do k=1,nz
      if (((max(CS%v_av(i,J,k), vm(i,J,k)) >= vel_rpt) .or. &
           (min(CS%v_av(i,J,k), vm(i,J,k)) <= -vel_rpt)) .and. &
          ((hin(i,j,k) + hin(i,j+1,k)) > 3.0*Angstrom)) exit
    enddo
    ks = k
    do k=nz,1,-1
      if (((max(CS%v_av(i,J,k), vm(i,J,k)) >= vel_rpt) .or. &
           (min(CS%v_av(i,J,k), vm(i,J,k)) <= -vel_rpt)) .and. &
          ((hin(i,j,k) + hin(i,j+1,k)) > 3.0*Angstrom)) exit
    enddo
    ke = k
    if (ke < ks) then
      ks = 1; ke = nz; write(file,'("V: Unable to set ks & ke.")')
    endif
    if (CS%full_column) then
      ks = 1 ; ke = nz
    endif

    call get_date(CS%Time, yr, mo, day, hr, minute, sec)
    call get_time((CS%Time - set_date(yr, 1, 1, 0, 0, 0)), sec, yearday)
    write (file,'(/,"--------------------------")')
    write (file,'(/,"Time ",i5,i4,F6.2," V-velocity violation at ",I4,": ",2(I3), &
        & " (",F7.2," E ",F7.2," N) Layers ",I3," to ",I3,". dt = ",1PG10.4)') &
        yr, yearday, (REAL(sec)/3600.0), pe_here(), i, J, &
        G%geoLonCv(i,J), G%geoLatCv(i,J), ks, ke, US%T_to_s*dt

    if (ks <= GV%nk_rho_varies) ks = 1
    do k=ks,ke
      if ((hin(i,j,k) + hin(i,j+1,k)) > 3.0*Angstrom) do_k(k) = .true.
    enddo

    write(file,'(/,"Layers:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(I10," ")', advance='no') (k) ; enddo
    write(file,'(/,"v(m):  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*vm(i,J,k)) ; enddo

    if (prev_avail) then
      write(file,'(/,"v(mp): ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_prev(i,J,k)) ; enddo
    endif

    write(file,'(/,"v(3):  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*CS%v_av(i,J,k)) ; enddo
    write(file,'(/,"CFL v: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) then
      CFL = abs(vm(i,J,k)) * dt * G%dx_Cv(i,J)
      if (vm(i,J,k) < 0.0) then ; CFL = CFL * G%IareaT(i,j+1)
      else ; CFL = CFL * G%IareaT(i,j) ; endif
      write(file,'(ES10.3," ")', advance='no') CFL
    endif ; enddo
    write(file,'(/,"CFL0 v:")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    abs(vm(i,J,k)) * dt * G%IdyCv(i,J) ; enddo

    if (prev_avail) then
      write(file,'(/,"dv:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*(vm(i,J,k)-CS%v_prev(i,J,k))) ; enddo
    endif

    write(file,'(/,"CAv:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%CAv(i,J,k)) ; enddo

    write(file,'(/,"PFv:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%PFv(i,J,k)) ; enddo

    write(file,'(/,"diffv: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') (vel_scale*dt*ADp%diffv(i,J,k)) ; enddo

    if (associated(ADp%gradKEv)) then
      write(file,'(/,"KEv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*ADp%gradKEv(i,J,k)) ; enddo
    endif
    if (associated(ADp%rv_x_u)) then
      write(file,'(/,"Corv:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                 vel_scale*dt*(ADp%CAv(i,J,k)-ADp%rv_x_u(i,J,k)) ; enddo
    endif
    if (associated(ADp%dv_dt_visc)) then
      write(file,'(/,"vbv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
          vel_scale*(vm(i,J,k) - dt*ADp%dv_dt_visc(i,J,k)) ; enddo

      write(file,'(/,"dvv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*ADp%dv_dt_visc(i,J,k)) ; enddo
    endif
    if (associated(ADp%dv_other)) then
      write(file,'(/,"dv_other: ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*ADp%dv_other(i,J,k)) ; enddo
    endif
    if (present(a)) then
      write(file,'(/,"a:     ",ES10.3," ")', advance='no') h_scale*a(i,J,ks)*dt
      do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') (h_scale*a(i,J,K)*dt) ; enddo
    endif
    if (present(hv)) then
      write(file,'(/,"hvel:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hv(i,J,k) ; enddo
    endif
    if (present(str)) then
      write(file,'(/,"Stress:  ",ES10.3)', advance='no') vel_scale*US%Z_to_m * (str*dt / GV%Rho0)
    endif

    if (associated(CS%v_accel_bt)) then
      write(file,'("dvbt:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                      (vel_scale*dt*CS%v_accel_bt(i,J,k)) ; enddo
    endif
    write(file,'(/)')

    write(file,'("h--:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i-1,j,k) ; enddo
    write(file,'(/,"h0-:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i,j,k) ; enddo
    write(file,'(/,"h+-:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i+1,j,k) ; enddo
    write(file,'(/,"h-+:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i-1,j+1,k) ; enddo
    write(file,'(/,"h0+:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i,j+1,k) ; enddo
    write(file,'(/,"h++:   ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') h_scale*hin(i+1,j+1,k) ; enddo

    e(nz+1) = -US%Z_to_m*(G%bathyT(i,j) + G%Z_ref)
    do k=nz,1,-1 ; e(K) = e(K+1) + h_scale*hin(i,j,k) ; enddo
    write(file,'(/,"e-:    ",ES10.3," ")', advance='no') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') e(K) ; enddo

    e(nz+1) = -US%Z_to_m*(G%bathyT(i,j+1) + G%Z_ref)
    do k=nz,1,-1 ; e(K) = e(K+1) + h_scale*hin(i,j+1,k) ; enddo
    write(file,'(/,"e+:    ",ES10.3," ")', advance='no') e(ks)
    do K=ks+1,ke+1 ; if (do_k(k-1)) write(file,'(ES10.3," ")', advance='no') e(K) ; enddo
    if (associated(CS%T)) then
      write(file,'(/,"T-:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') temp_scale*CS%T(i,j,k) ; enddo
      write(file,'(/,"T+:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') temp_scale*CS%T(i,j+1,k) ; enddo
    endif
    if (associated(CS%S)) then
      write(file,'(/,"S-:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') saln_scale*CS%S(i,j,k) ; enddo
      write(file,'(/,"S+:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') saln_scale*CS%S(i,j+1,k) ; enddo
    endif

    if (prev_avail) then
      write(file,'(/,"u--:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') vel_scale*CS%u_prev(I-1,j,k) ; enddo
      write(file,'(/,"u-+:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') vel_scale*CS%u_prev(I-1,j+1,k) ; enddo
      write(file,'(/,"u+-:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') vel_scale*CS%u_prev(I,j,k) ; enddo
      write(file,'(/,"u++:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') vel_scale*CS%u_prev(I,j+1,k) ; enddo
    endif

    write(file,'(/,"uh--:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%uh(I-1,j,k)*G%IdyCu(I-1,j)) ; enddo
    write(file,'(/," uhC--: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_av(I-1,j,k) * uh_scale*0.5*(hin(i-1,j,k) + hin(i,j,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," uhCp--:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_prev(I-1,j,k) * uh_scale*0.5*(hin(i-1,j,k) + hin(i,j,k))) ; enddo
    endif

    write(file,'(/,"uh-+:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%uh(I-1,j+1,k)*G%IdyCu(I-1,j+1)) ; enddo
    write(file,'(/," uhC-+: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_av(I-1,j+1,k) * uh_scale*0.5*(hin(i-1,j+1,k) + hin(i,j+1,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," uhCp-+:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_prev(I-1,j+1,k) * uh_scale*0.5*(hin(i-1,j+1,k) + hin(i,j+1,k))) ; enddo
    endif

    write(file,'(/,"uh+-:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%uh(I,j,k)*G%IdyCu(I,j)) ; enddo
    write(file,'(/," uhC+-: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_av(I,j,k) * uh_scale*0.5*(hin(i,j,k) + hin(i+1,j,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," uhCp+-:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_prev(I,j,k) * uh_scale*0.5*(hin(i,j,k) + hin(i+1,j,k))) ; enddo
    endif

    write(file,'(/,"uh++:  ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
                                    (uh_scale*CDp%uh(I,j+1,k)*G%IdyCu(I,j+1)) ; enddo
    write(file,'(/," uhC++: ")', advance='no')
    do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_av(I,j+1,k) * uh_scale*0.5*(hin(i,j+1,k) + hin(i+1,j+1,k))) ; enddo
    if (prev_avail) then
      write(file,'(/," uhCp++:")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(ES10.3," ")', advance='no') &
            (CS%u_prev(I,j+1,k) * uh_scale*0.5*(hin(i,j+1,k) + hin(i+1,j+1,k))) ; enddo
    endif

    write(file,'(/,"D:     ",2(ES10.3))') US%Z_to_m*(G%bathyT(i,j) + G%Z_ref), US%Z_to_m*(G%bathyT(i,j+1) + G%Z_ref)

  !  From here on, the normalized accelerations are written.
    if (prev_avail) then
      do k=ks,ke
        dv = vm(i,J,k) - CS%v_prev(i,J,k)
        if (abs(dv) < 1.0e-6*US%m_s_to_L_T) dv = 1.0e-6*US%m_s_to_L_T
        Inorm(k) = 1.0 / dv
      enddo

      write(file,'(2/,"Norm:  ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') (vel_scale / Inorm(k)) ; enddo
      write(file,'(/,"dv:    ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      ((vm(i,J,k)-CS%v_prev(i,J,k)) * Inorm(k)) ; enddo
      write(file,'(/,"CAv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%CAv(i,J,k) * Inorm(k)) ; enddo
      write(file,'(/,"PFv:   ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%PFv(i,J,k) * Inorm(k)) ; enddo
      write(file,'(/,"diffv: ")', advance='no')
      do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                      (dt*ADp%diffv(i,J,k) * Inorm(k)) ; enddo

      if (associated(ADp%gradKEu)) then
        write(file,'(/,"KEv:   ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*ADp%gradKEv(i,J,k) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%rv_x_u)) then
        write(file,'(/,"Corv:  ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*(ADp%CAv(i,J,k)-ADp%rv_x_u(i,J,k)) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%dv_dt_visc)) then
        write(file,'(/,"dvv:   ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*ADp%dv_dt_visc(i,J,k) * Inorm(k)) ; enddo
      endif
      if (associated(ADp%dv_other)) then
        write(file,'(/,"dv_other: ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (ADp%dv_other(i,J,k) * Inorm(k)) ; enddo
      endif
      if (associated(CS%v_accel_bt)) then
        write(file,'(/,"dvbt:  ")', advance='no')
        do k=ks,ke ; if (do_k(k)) write(file,'(F10.6," ")', advance='no') &
                                        (dt*CS%v_accel_bt(i,J,k) * Inorm(k)) ; enddo
      endif
    endif

    write(file,'(2/)')

    flush(file)
  endif

end subroutine write_v_accel

!> This subroutine initializes the parameters regulating how truncations are logged.
subroutine PointAccel_init(MIS, Time, G, param_file, diag, dirs, CS)
  type(ocean_internal_state), &
                        target, intent(in)    :: MIS  !< For "MOM Internal State" a set of pointers
                                                      !! to the fields and accelerations that make
                                                      !! up the ocean's physical state.
  type(time_type),      target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),        intent(in)    :: G    !< The ocean's grid structure.
  type(param_file_type),        intent(in)    :: param_file !< A structure to parse for run-time
                                                      !! parameters.
  type(diag_ctrl),      target, intent(inout) :: diag !< A structure that is used to regulate
                                                      !! diagnostic output.
  type(directories),            intent(in)    :: dirs !< A structure containing several relevant
                                                      !! directory paths.
  type(PointAccel_CS),          pointer       :: CS   !< A pointer that is set to point to the
                                                      !! control structure for this module.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_PointAccel" ! This module's name.

  if (associated(CS)) return
  allocate(CS)

  CS%diag => diag ; CS%Time => Time

  CS%T => MIS%T ; CS%S => MIS%S
  CS%u_accel_bt => MIS%u_accel_bt ; CS%v_accel_bt => MIS%v_accel_bt
  CS%u_prev => MIS%u_prev ; CS%v_prev => MIS%v_prev
  CS%u_av => MIS%u_av; if (.not.associated(MIS%u_av)) CS%u_av => MIS%u(:,:,:)
  CS%v_av => MIS%v_av; if (.not.associated(MIS%v_av)) CS%v_av => MIS%v(:,:,:)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "", debugging=.true.)
  call get_param(param_file, mdl, "U_TRUNC_FILE", CS%u_trunc_file, &
                 "The absolute path to the file where the accelerations "//&
                 "leading to zonal velocity truncations are written. \n"//&
                 "Leave this empty for efficiency if this diagnostic is "//&
                 "not needed.", default="", debuggingParam=.true.)
  call get_param(param_file, mdl, "V_TRUNC_FILE", CS%v_trunc_file, &
                 "The absolute path to the file where the accelerations "//&
                 "leading to meridional velocity truncations are written. \n"//&
                 "Leave this empty for efficiency if this diagnostic is "//&
                 "not needed.", default="", debuggingParam=.true.)
  call get_param(param_file, mdl, "MAX_TRUNC_FILE_SIZE_PER_PE", CS%max_writes, &
                 "The maximum number of columns of truncations that any PE "//&
                 "will write out during a run.", default=50, debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_FULL_COLUMN", CS%full_column, &
                 "If true, write out the accelerations in all massive layers; otherwise "//&
                 "just document the ones with large velocities.", &
                 default=.false., debuggingParam=.true.)

  if (len_trim(dirs%output_directory) > 0) then
    if (len_trim(CS%u_trunc_file) > 0) &
      CS%u_trunc_file = trim(dirs%output_directory)//trim(CS%u_trunc_file)
    if (len_trim(CS%v_trunc_file) > 0) &
      CS%v_trunc_file = trim(dirs%output_directory)//trim(CS%v_trunc_file)
    call log_param(param_file, mdl, "output_dir/U_TRUNC_FILE", CS%u_trunc_file, debuggingParam=.true.)
    call log_param(param_file, mdl, "output_dir/V_TRUNC_FILE", CS%v_trunc_file, debuggingParam=.true.)
  endif
  CS%u_file = -1 ; CS%v_file = -1 ; CS%cols_written = 0

end subroutine PointAccel_init

end module MOM_PointAccel
