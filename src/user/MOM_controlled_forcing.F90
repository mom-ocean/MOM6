module MOM_controlled_forcing
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
use MOM_diag_mediator, only : post_data, query_averaging_enabled
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_domains, only : pass_var, pass_vector, AGRID, To_South, To_West, To_All
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_io, only : vardesc, var_desc
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

public apply_ctrl_forcing, register_ctrl_forcing_restarts
public controlled_forcing_init, controlled_forcing_end

type, public :: ctrl_forcing_CS ; private
  logical :: use_temperature ! If true, temperature and salinity are used as
                             ! state variables.
  logical :: do_integrated   ! If true, use time-integrated anomalies to control
                             ! the surface state.
  integer :: num_cycle       ! The number of elements in the forcing cycle.
  real    :: heat_int_rate  ! The rate at which heating anomalies accumulate, in s-1.
  real    :: prec_int_rate  ! The rate at which precipitation anomalies accumulate, in s-1.
  real    :: heat_cyc_rate  ! The rate at which cyclical heating anomaliess
                            ! accumulate, in s-1.
  real    :: prec_cyc_rate  ! The rate at which cyclical precipitation anomaliess
                            ! accumulate, in s-1.
  real    :: Len2           ! The square of the length scale over which the anomalies
                            ! are smoothed via a Laplacian filter, in m2.
  real    :: lam_heat       ! A constant of proportionality between SST anomalies
                            ! and heat fluxes, in W m-2 K-1.
  real    :: lam_prec       ! A constant of proportionality between SSS anomalies
                            ! (normalised by mean SSS) and precipitation, in kg m-2.
  real    :: lam_cyc_heat   ! A constant of proportionality between cyclical SST
                            ! anomalies and corrective heat fluxes, in W m-2 K-1.
  real    :: lam_cyc_prec   ! A constant of proportionality between cyclical SSS
                            ! anomalies (normalised by mean SSS) and corrective
                            ! precipitation, in kg m-2.

  real, pointer, dimension(:) :: &
    avg_time => NULL()
  real, pointer, dimension(:,:) :: &
    heat_0 => NULL(), &
    precip_0 => NULL()
  real, pointer, dimension(:,:,:) :: &
    heat_cyc => NULL(), &
    precip_cyc => NULL(), &
    avg_SST_anom => NULL(), &
    avg_SSS_anom => NULL(), &
    avg_SSS => NULL()
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer :: id_heat_0 = -1 ! See if these are neede later...
end type ctrl_forcing_CS

contains

!> This subroutine calls any of the other subroutines in this file
!! that are needed to specify the current surface forcing fields.
subroutine apply_ctrl_forcing(SST_anom, SSS_anom, SSS_mean, virt_heat, virt_precip, &
                              day_start, dt, G, CS)
  type(ocean_grid_type), intent(inout) :: G                    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: SST_anom  !< The sea surface temperature
                                                               !! anomalies, in deg C.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: SSS_anom  !< The sea surface salinity
                                                               !! anomlies, in g kg-1.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: SSS_mean  !< The mean sea surface
                                                               !! salinity, in g kg-1.
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: virt_heat !< Virtual (corrective) heat
                                                               !! fluxes that are augmented
                                                               !! in this subroutine, in W m-2.
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: virt_precip !< Virtual (corrective)
                                                               !! precipitation fluxes that
                                                               !! are augmented in this
                                                               !! subroutine, in kg m-2 s-1.
  type(time_type),       intent(in)    :: day_start      !< Start time of the fluxes.
  real,                  intent(in)    :: dt             !< Length of time over which these
                                                         !! fluxes will be applied, in s.
  type(ctrl_forcing_CS), pointer       :: CS             !< A pointer to the control structure
                                                         !! returned by a previous call to
                                                         !! ctrl_forcing_init.
!
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    flux_heat_x, &
    flux_prec_x
  real, dimension(SZI_(G),SZJB_(G)) :: &
    flux_heat_y, &
    flux_prec_y
  type(time_type) :: day_end
  real    :: coef    ! A heat-flux coefficient with units of m2.
  real    :: mr_st, mr_end, mr_mid, mr_prev, mr_next
  real    :: dt_wt, dt_heat_rate, dt_prec_rate
  real    :: dt1_heat_rate, dt1_prec_rate, dt2_heat_rate, dt2_prec_rate
  real    :: wt_per1, wt_st, wt_end, wt_mid
  integer :: m_st, m_end, m_mid, m_u1, m_u2, m_u3
  integer :: yr, mon, day, hr, min, sec
  integer :: i, j, is, ie, js, je

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) return
  if ((CS%num_cycle <= 0) .and. (.not.CS%do_integrated)) return

  day_end = day_start + set_time(floor(dt+0.5))

  do j=js,je ; do i=is,ie
    virt_heat(i,j) = 0.0 ; virt_precip(i,j) = 0.0
  enddo ; enddo

  if (CS%do_integrated) then
    dt_heat_rate = dt * CS%heat_int_rate
    dt_prec_rate = dt * CS%prec_int_rate
    call pass_var(CS%heat_0, G%Domain, complete=.false.)
    call pass_var(CS%precip_0, G%Domain)

    do j=js,je ; do I=is-1,ie
      coef = CS%Len2 * (G%dy_Cu(I,j)*G%IdxCu(I,j))
      flux_heat_x(I,j) = coef * (CS%heat_0(i,j) - CS%heat_0(i+1,j))
      flux_prec_x(I,j) = coef * (CS%precip_0(i,j) - CS%precip_0(i+1,j))
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      coef = CS%Len2 * (G%dx_Cv(i,J)*G%IdyCv(i,J))
      flux_heat_y(i,J) = coef * (CS%heat_0(i,j) - CS%heat_0(i,j+1))
      flux_prec_y(i,J) = coef * (CS%precip_0(i,j) - CS%precip_0(i,j+1))
    enddo ; enddo
    do j=js,je ; do i=is,ie
      CS%heat_0(i,j) = CS%heat_0(i,j) + dt_heat_rate * ( &
         -CS%lam_heat*G%mask2dT(i,j)*SST_anom(i,j) + &
        (G%IareaT(i,j) * ((flux_heat_x(I-1,j) - flux_heat_x(I,j)) + &
                          (flux_heat_y(i,J-1) - flux_heat_y(i,J))) ) )

      CS%precip_0(i,j) = CS%precip_0(i,j) + dt_prec_rate * ( &
         CS%lam_prec * G%mask2dT(i,j)*(SSS_anom(i,j) / SSS_mean(i,j)) + &
        (G%IareaT(i,j) * ((flux_prec_x(I-1,j) - flux_prec_x(I,j)) + &
                          (flux_prec_y(i,J-1) - flux_prec_y(i,J))) ) )

      virt_heat(i,j) = virt_heat(i,j) + CS%heat_0(i,j)
      virt_precip(i,j) = virt_precip(i,j) + CS%precip_0(i,j)
    enddo ; enddo
  endif

  if (CS%num_cycle > 0) then
    ! Determine the current period, with values that run from 0 to CS%num_cycle.
    call get_date(day_start, yr, mon, day, hr, min, sec)
    mr_st = CS%num_cycle * (time_type_to_real(day_start - set_date(yr, 1, 1)) / &
                   time_type_to_real(set_date(yr+1, 1, 1) - set_date(yr, 1, 1)))

    call get_date(day_end, yr, mon, day, hr, min, sec)
    mr_end = CS%num_cycle * (time_type_to_real(day_end - set_date(yr, 1, 1)) / &
                   time_type_to_real(set_date(yr+1, 1, 1) - set_date(yr, 1, 1)))

    ! The Chapeau functions are centered at whole integer values that are nominally
    ! the end of the month to enable simple conversion from the fractional-years times
    ! CS%num_cycle.

    ! The month-average temperatures have as an index the month number.

    m_end = periodic_int(real(ceiling(mr_end)), CS%num_cycle)
    m_mid = periodic_int(real(ceiling(mr_st)), CS%num_cycle)
    m_st = periodic_int(mr_st, CS%num_cycle)

    mr_st = periodic_real(mr_st, CS%num_cycle)
    mr_end = periodic_real(mr_end, CS%num_cycle)
      !  mr_mid = periodic_real(ceiling(mr_st), CS%num_cycle)
    mr_prev = periodic_real(real(floor(mr_st)), CS%num_cycle)
    mr_next = periodic_real(real(m_end), CS%num_cycle)
    if (m_mid == m_end) then ; mr_mid = mr_end ! There is only one cell.
    else ; mr_mid = periodic_real(real(m_mid), CS%num_cycle) ; endif

    ! There may be two cells that run from mr_st to mr_mid and mr_mid to mr_end.

    ! The values of m for weights are all calculated relative to mr_prev, so
    ! check whether mr_mid, etc., need to be shifted by CS%num_cycle, so that these
    ! values satisfiy  mr_prev <= mr_st < mr_mid <= mr_end <= mr_next.
    if (mr_st < mr_prev) mr_prev = mr_prev - CS%num_cycle
    if (mr_mid < mr_st) mr_mid = mr_mid + CS%num_cycle
    if (mr_end < mr_st) mr_end = mr_end + CS%num_cycle
    if (mr_next < mr_prev) mr_next = mr_next + CS%num_cycle

    !### These might be removed later - they are to check the coding.
    if ((mr_mid < mr_st) .or. (mr_mid > mr_prev + 1.)) call MOM_error(FATAL, &
          "apply ctrl_forcing: m_mid interpolation out of bounds; fix the code.")
    if ((mr_end < mr_st) .or. (mr_end > mr_prev + 2.)) call MOM_error(FATAL, &
          "apply ctrl_forcing: m_end interpolation out of bounds; fix the code.")
    if (mr_end > mr_next) call MOM_error(FATAL, &
          "apply ctrl_forcing: mr_next interpolation out of bounds; fix the code.")

    wt_per1 = 1.0
    if (mr_mid < mr_end) wt_per1 = (mr_mid - mr_st) / (mr_end - mr_st)

    ! Find the 3 Chapeau-function weights, bearing in mind that m_end may be m_mid.
    wt_st = wt_per1 * (1. + (mr_prev - 0.5*(mr_st + mr_mid)))
    wt_end = (1.0-wt_per1) * (1. + (0.5*(mr_end + mr_mid) - mr_next))
    wt_mid = 1.0 - (wt_st + wt_end)
    if ((wt_st < 0.0) .or. (wt_end < 0.0) .or. (wt_mid < 0.0)) &
      call MOM_error(FATAL, "apply_ctrl_forcing: Negative m weights")
    if ((wt_st > 1.0) .or. (wt_end > 1.0) .or. (wt_mid > 1.0)) &
      call MOM_error(FATAL, "apply_ctrl_forcing: Excessive m weights")

    ! Add to vert_heat and vert_precip.
    do j=js,je ; do i=is,ie
      virt_heat(i,j) = virt_heat(i,j) + (wt_st * CS%heat_cyc(i,j,m_st) + &
                        (wt_mid * CS%heat_cyc(i,j,m_mid) + &
                         wt_end * CS%heat_cyc(i,j,m_end)))
      virt_precip(i,j) = virt_precip(i,j) + (wt_st * CS%precip_cyc(i,j,m_st) + &
                        (wt_mid * CS%precip_cyc(i,j,m_mid) + &
                         wt_end * CS%precip_cyc(i,j,m_end)))
    enddo ; enddo

    ! If different from the last period, take the average and determine the
    ! chapeau weighting

    ! The Chapeau functions are centered at whole integer values that are nominally
    ! the end of the month to enable simple conversion from the fractional-years times
    ! CS%num_cycle.

    ! The month-average temperatures have as an index the month number, so the averages
    ! apply to indicies m_end and m_mid.

    if (CS%avg_time(m_end) <= 0.0) then ! zero out the averages.
      CS%avg_time(m_end) = 0.0
      do j=js,je ; do i=is,ie
        CS%avg_SST_anom(i,j,m_end) = 0.0
        CS%avg_SSS_anom(i,j,m_end) = 0.0 ; CS%avg_SSS(i,j,m_end) = 0.0
      enddo ; enddo
    endif
    if (CS%avg_time(m_mid) <= 0.0) then ! zero out the averages.
      CS%avg_time(m_mid) = 0.0
      do j=js,je ; do i=is,ie
        CS%avg_SST_anom(i,j,m_mid) = 0.0
        CS%avg_SSS_anom(i,j,m_mid) = 0.0 ; CS%avg_SSS(i,j,m_mid) = 0.0
      enddo ; enddo
    endif

    ! Accumulate the average anomalies for this period.
    dt_wt = wt_per1 * dt
    CS%avg_time(m_mid) = CS%avg_time(m_mid) + dt_wt
    do j=js,je ; do i=is,ie
      CS%avg_SST_anom(i,j,m_mid) = CS%avg_SST_anom(i,j,m_mid) + &
                                   dt_wt * G%mask2dT(i,j) * SST_anom(i,j)
      CS%avg_SSS_anom(i,j,m_mid) = CS%avg_SSS_anom(i,j,m_mid) + &
                                   dt_wt * G%mask2dT(i,j) * SSS_anom(i,j)
      CS%avg_SSS(i,j,m_mid) = CS%avg_SSS(i,j,m_mid) + dt_wt * SSS_mean(i,j)
    enddo ; enddo
    if (wt_per1 < 1.0) then
      dt_wt = (1.0-wt_per1) * dt
      CS%avg_time(m_end) = CS%avg_time(m_end) + dt_wt
      do j=js,je ; do i=is,ie
        CS%avg_SST_anom(i,j,m_end) = CS%avg_SST_anom(i,j,m_end) + &
                                     dt_wt * G%mask2dT(i,j) * SST_anom(i,j)
        CS%avg_SSS_anom(i,j,m_end) = CS%avg_SSS_anom(i,j,m_end) + &
                                     dt_wt * G%mask2dT(i,j) * SSS_anom(i,j)
        CS%avg_SSS(i,j,m_end) = CS%avg_SSS(i,j,m_end) + dt_wt * SSS_mean(i,j)
      enddo ; enddo
    endif

    ! Update the Chapeau magnitudes for 4 cycles ago.
    m_u1 = periodic_int(m_st - 4.0, CS%num_cycle)
    m_u2 = periodic_int(m_st - 3.0, CS%num_cycle)
    m_u3 = periodic_int(m_st - 2.0, CS%num_cycle)

    if (CS%avg_time(m_u1) > 0.0) then
      do j=js,je ; do i=is,ie
        CS%avg_SST_anom(i,j,m_u1) = CS%avg_SST_anom(i,j,m_u1) / CS%avg_time(m_u1)
        CS%avg_SSS_anom(i,j,m_u1) = CS%avg_SSS_anom(i,j,m_u1) / CS%avg_time(m_u1)
        CS%avg_SSS(i,j,m_u1) = CS%avg_SSS(i,j,m_u1) / CS%avg_time(m_u1)
      enddo ; enddo
      CS%avg_time(m_u1) = -1.0
    endif
    if (CS%avg_time(m_u2) > 0.0) then
      do j=js,je ; do i=is,ie
        CS%avg_SST_anom(i,j,m_u2) = CS%avg_SST_anom(i,j,m_u2) / CS%avg_time(m_u2)
        CS%avg_SSS_anom(i,j,m_u2) = CS%avg_SSS_anom(i,j,m_u2) / CS%avg_time(m_u2)
        CS%avg_SSS(i,j,m_u2) = CS%avg_SSS(i,j,m_u2) / CS%avg_time(m_u2)
      enddo ; enddo
      CS%avg_time(m_u2) = -1.0
    endif
    if (CS%avg_time(m_u3) > 0.0) then
      do j=js,je ; do i=is,ie
        CS%avg_SST_anom(i,j,m_u3) = CS%avg_SST_anom(i,j,m_u3) / CS%avg_time(m_u3)
        CS%avg_SSS_anom(i,j,m_u3) = CS%avg_SSS_anom(i,j,m_u3) / CS%avg_time(m_u3)
        CS%avg_SSS(i,j,m_u3) = CS%avg_SSS(i,j,m_u3) / CS%avg_time(m_u3)
      enddo ; enddo
      CS%avg_time(m_u3) = -1.0
    endif

    dt1_heat_rate = wt_per1 * dt * CS%heat_cyc_rate
    dt1_prec_rate = wt_per1 * dt * CS%prec_cyc_rate
    dt2_heat_rate = (1.0-wt_per1) * dt * CS%heat_cyc_rate
    dt2_prec_rate = (1.0-wt_per1) * dt * CS%prec_cyc_rate

    if (wt_per1 < 1.0) then
      call pass_var(CS%heat_cyc(:,:,m_u2), G%Domain, complete=.false.)
      call pass_var(CS%precip_cyc(:,:,m_u2), G%Domain, complete=.false.)
    endif
    call pass_var(CS%heat_cyc(:,:,m_u1), G%Domain, complete=.false.)
    call pass_var(CS%precip_cyc(:,:,m_u1), G%Domain)

    if ((CS%avg_time(m_u1) == -1.0) .and. (CS%avg_time(m_u2) == -1.0)) then
      do j=js,je ; do I=is-1,ie
        coef = CS%Len2 * (G%dy_Cu(I,j)*G%IdxCu(I,j))
        flux_heat_x(I,j) = coef * (CS%heat_cyc(i,j,m_u1) - CS%heat_cyc(i+1,j,m_u1))
        flux_prec_x(I,j) = coef * (CS%precip_cyc(i,j,m_u1) - CS%precip_cyc(i+1,j,m_u1))
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        coef = CS%Len2 * (G%dx_Cv(i,J)*G%IdyCv(i,J))
        flux_heat_y(i,J) = coef * (CS%heat_cyc(i,j,m_u1) - CS%heat_cyc(i,j+1,m_u1))
        flux_prec_y(i,J) = coef * (CS%precip_cyc(i,j,m_u1) - CS%precip_cyc(i,j+1,m_u1))
      enddo ; enddo
      do j=js,je ; do i=is,ie
        CS%heat_cyc(i,j,m_u1) = CS%heat_cyc(i,j,m_u1) + dt1_heat_rate * ( &
           -CS%lam_cyc_heat*(CS%avg_SST_anom(i,j,m_u2) - CS%avg_SST_anom(i,j,m_u1)) + &
          (G%IareaT(i,j) * ((flux_heat_x(I-1,j) - flux_heat_x(I,j)) + &
                            (flux_heat_y(i,J-1) - flux_heat_y(i,J))) ) )

        CS%precip_cyc(i,j,m_u1) = CS%precip_cyc(i,j,m_u1) + dt1_prec_rate * ( &
          CS%lam_cyc_prec * (CS%avg_SSS_anom(i,j,m_u2) - CS%avg_SSS_anom(i,j,m_u1)) / &
                            (0.5*(CS%avg_SSS(i,j,m_u2) + CS%avg_SSS(i,j,m_u1))) + &
          (G%IareaT(i,j) * ((flux_prec_x(I-1,j) - flux_prec_x(I,j)) + &
                            (flux_prec_y(i,J-1) - flux_prec_y(i,J))) ) )
      enddo ; enddo
    endif

    if ((wt_per1 < 1.0) .and. (CS%avg_time(m_u1) == -1.0) .and. (CS%avg_time(m_u2) == -1.0))  then
      do j=js,je ; do I=is-1,ie
        coef = CS%Len2 * (G%dy_Cu(I,j)*G%IdxCu(I,j))
        flux_heat_x(I,j) = coef * (CS%heat_cyc(i,j,m_u2) - CS%heat_cyc(i+1,j,m_u2))
        flux_prec_x(I,j) = coef * (CS%precip_cyc(i,j,m_u2) - CS%precip_cyc(i+1,j,m_u2))
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        coef = CS%Len2 * (G%dx_Cv(i,J)*G%IdyCv(i,J))
        flux_heat_y(i,J) = coef * (CS%heat_cyc(i,j,m_u2) - CS%heat_cyc(i,j+1,m_u2))
        flux_prec_y(i,J) = coef * (CS%precip_cyc(i,j,m_u2) - CS%precip_cyc(i,j+1,m_u2))
      enddo ; enddo
      do j=js,je ; do i=is,ie
        CS%heat_cyc(i,j,m_u2) = CS%heat_cyc(i,j,m_u2) + dt1_heat_rate * ( &
         -CS%lam_cyc_heat*(CS%avg_SST_anom(i,j,m_u3) - CS%avg_SST_anom(i,j,m_u2)) + &
          (G%IareaT(i,j) * ((flux_heat_x(I-1,j) - flux_heat_x(I,j)) + &
                            (flux_heat_y(i,J-1) - flux_heat_y(i,J))) ) )

        CS%precip_cyc(i,j,m_u2) = CS%precip_cyc(i,j,m_u2) + dt1_prec_rate * ( &
          CS%lam_cyc_prec * (CS%avg_SSS_anom(i,j,m_u3) - CS%avg_SSS_anom(i,j,m_u2)) / &
                             (0.5*(CS%avg_SSS(i,j,m_u3) + CS%avg_SSS(i,j,m_u2))) + &
          (G%IareaT(i,j) * ((flux_prec_x(I-1,j) - flux_prec_x(I,j)) + &
                            (flux_prec_y(i,J-1) - flux_prec_y(i,J))) ) )
      enddo ; enddo
    endif

  endif ! (CS%num_cycle > 0)

end subroutine apply_ctrl_forcing

!> This function maps rval into an integer in the range from 1 to num_period.
function periodic_int(rval, num_period) result (m)
  real,    intent(in) :: rval       !< Input for mapping.
  integer, intent(in) :: num_period !< Maximum output.
  integer             :: m          !< Return value.

  m = floor(rval)
  if (m <= 0) then
    m = m + num_period * (1 + (abs(m) / num_period))
  elseif (m > num_period) then
    m = m - num_period * ((m-1) / num_period)
  endif
end function

!> This function shifts rval by an integer multiple of num_period so that
!! 0 <= val_out < num_period.
function periodic_real(rval, num_period) result(val_out)
  real,    intent(in) :: rval       !< Input to be shifted into valid range.
  integer, intent(in) :: num_period !< Maximum valid value.
  real                :: val_out    !< Return value.
  integer :: nshft

  if (rval < 0) then ; nshft = floor(abs(rval) / num_period) + 1
  elseif (rval < num_period) then ; nshft = 0
  else ; nshft = -1*floor(rval / num_period) ; endif

  val_out = rval + nshft * num_period
end function


!> This subroutine is used to allocate and register any fields in this module
!! that should be written to or read from the restart file.
subroutine register_ctrl_forcing_restarts(G, param_file, CS, restart_CS)
  type(ocean_grid_type), intent(in) :: G          !< The ocean's grid structure.
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.
  type(ctrl_forcing_CS), pointer :: CS            !< A pointer that is set to point to the
                                                  !! control structure for this module.
  type(MOM_restart_CS),  pointer :: restart_CS    !< A pointer to the restart control structure.

  logical :: controlled, use_temperature
  character (len=8) :: period_str
  type(vardesc) :: vd
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(WARNING, "register_ctrl_forcing_restarts called "//&
                             "with an associated control structure.")
    return
  endif

  controlled = .false.
  call read_param(param_file, "CONTROLLED_FORCING", controlled)
  if (.not.controlled) return

  use_temperature = .true.
  call read_param(param_file, "ENABLE_THERMODYNAMICS", use_temperature)
  if (.not.use_temperature) call MOM_error(FATAL, &
    "register_ctrl_forcing_restarts: CONTROLLED_FORCING only works with "//&
    "ENABLE_THERMODYNAMICS defined.")

  allocate(CS)

  CS%do_integrated = .true. ; CS%num_cycle = 0
  call read_param(param_file, "CTRL_FORCE_INTEGRATED", CS%do_integrated)
  call read_param(param_file, "CTRL_FORCE_NUM_CYCLE", CS%num_cycle)

  if (CS%do_integrated) then
    call safe_alloc_ptr(CS%heat_0,isd,ied,jsd,jed) ; CS%heat_0(:,:) = 0.0
    call safe_alloc_ptr(CS%precip_0,isd,ied,jsd,jed) ; CS%precip_0(:,:) = 0.0
    vd = var_desc("Ctrl_heat","W m-2","Control Integrative Heating",z_grid='1')
    call register_restart_field(CS%heat_0, vd, .false., restart_CS)
    vd = var_desc("Ctrl_precip","kg m-2 s-1","Control Integrative Precipitation",z_grid='1')
    call register_restart_field(CS%precip_0, vd, .false., restart_CS)
  endif

  if (CS%num_cycle > 0) then
    write (period_str, '(i8)') CS%num_cycle
    period_str = trim('p ')//trim(adjustl(period_str))
    call safe_alloc_ptr(CS%heat_cyc,isd,ied,jsd,jed,CS%num_cycle) ; CS%heat_cyc(:,:,:) = 0.0
    call safe_alloc_ptr(CS%precip_cyc,isd,ied,jsd,jed,CS%num_cycle) ; CS%precip_cyc(:,:,:) = 0.0
    vd = var_desc("Ctrl_heat_cycle", "W m-2","Cyclical Control Heating",&
                  z_grid='1', t_grid=period_str)
    call register_restart_field(CS%heat_cyc, vd, .false., restart_CS)
    vd = var_desc("Ctrl_precip_cycle","kg m-2 s-1","Cyclical Control Precipitation", &
                  z_grid='1', t_grid=period_str)
    call register_restart_field(CS%precip_cyc, vd, .false., restart_CS)

    call safe_alloc_ptr(CS%avg_time,CS%num_cycle) ; CS%avg_time(:) = 0.0
    vd = var_desc("avg_time","sec","Cyclical accumulated averaging time", &
                  '1',z_grid='1',t_grid=period_str)
    call register_restart_field(CS%avg_time, vd, .false., restart_CS)

    call safe_alloc_ptr(CS%avg_SST_anom,isd,ied,jsd,jed,CS%num_cycle) ; CS%avg_SST_anom(:,:,:) = 0.0
    call safe_alloc_ptr(CS%avg_SSS_anom,isd,ied,jsd,jed,CS%num_cycle) ; CS%avg_SSS_anom(:,:,:) = 0.0
    vd = var_desc("avg_SST_anom","deg C","Cyclical average SST Anomaly", &
                  z_grid='1',t_grid=period_str)
    call register_restart_field(CS%avg_SST_anom, vd, .false., restart_CS)
    vd = var_desc("avg_SSS_anom","g kg-1","Cyclical average SSS Anomaly", &
                  z_grid='1',t_grid=period_str)
    call register_restart_field(CS%avg_SSS_anom, vd, .false., restart_CS)
  endif

end subroutine register_ctrl_forcing_restarts

!> Set up this modules control structure.
subroutine controlled_forcing_init(Time, G, param_file, diag, CS)
  type(time_type),           intent(in) :: Time       !< The current model time.
  type(ocean_grid_type),     intent(in) :: G          !< The ocean's grid structure.
  type(param_file_type),     intent(in) :: param_file !< A structure indicating the
                                                      !! open file to parse for model
                                                      !! parameter values.
  type(diag_ctrl), target,   intent(in) :: diag       !< A structure that is used to regulate
                                                      !! diagnostic output.
  type(ctrl_forcing_CS),     pointer    :: CS         !< A pointer that is set to point to the
                                                      !! control structure for this module.
  real :: smooth_len
  logical :: do_integrated
  integer :: num_cycle
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_controlled_forcing" ! This module's name.

  ! These should have already been called.
  ! call read_param(param_file, "CTRL_FORCE_INTEGRATED", CS%do_integrated)
  ! call read_param(param_file, "CTRL_FORCE_NUM_CYCLE", CS%num_cycle)

  if (associated(CS)) then
    do_integrated = CS%do_integrated ; num_cycle = CS%num_cycle
  else
    do_integrated = .false. ; num_cycle = 0
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call log_param(param_file, mod, "CTRL_FORCE_INTEGRATED", do_integrated, &
                 "If true, use a PI controller to determine the surface \n"//&
                 "forcing that is consistent with the observed mean properties.", &
                 default=.false.)
  call log_param(param_file, mod, "CTRL_FORCE_NUM_CYCLE", num_cycle, &
                 "The number of cycles per year in the controlled forcing, \n"//&
                 "or 0 for no cyclic forcing.", default=0)

  if (.not.associated(CS)) return

  CS%diag => diag

  call get_param(param_file, mod, "CTRL_FORCE_HEAT_INT_RATE", CS%heat_int_rate, &
                 "The integrated rate at which heat flux anomalies are \n"//&
                 "accumulated.", units="s-1", default=0.0)
  call get_param(param_file, mod, "CTRL_FORCE_PREC_INT_RATE", CS%prec_int_rate, &
                 "The integrated rate at which precipitation anomalies \n"//&
                 "are accumulated.", units="s-1", default=0.0)
  call get_param(param_file, mod, "CTRL_FORCE_HEAT_CYC_RATE", CS%heat_cyc_rate, &
                 "The integrated rate at which cyclical heat flux \n"//&
                 "anomalies are accumulated.", units="s-1", default=0.0)
  call get_param(param_file, mod, "CTRL_FORCE_PREC_CYC_RATE", CS%prec_cyc_rate, &
                 "The integrated rate at which cyclical precipitation \n"//&
                 "anomalies are accumulated.", units="s-1", default=0.0)
  call get_param(param_file, mod, "CTRL_FORCE_SMOOTH_LENGTH", smooth_len, &
                 "The length scales over which controlled forcing \n"//&
                 "anomalies are smoothed.", units="m", default=0.0)
  call get_param(param_file, mod, "CTRL_FORCE_LAMDA_HEAT", CS%lam_heat, &
                 "A constant of proportionality between SST anomalies \n"//&
                 "and controlling heat fluxes", "W m-2 K-1", default=0.0)
  call get_param(param_file, mod, "CTRL_FORCE_LAMDA_PREC", CS%lam_prec, &
                 "A constant of proportionality between SSS anomalies \n"//&
                 "(normalised by mean SSS) and controlling precipitation.", &
                 "kg m-2", default=0.0)
  call get_param(param_file, mod, "CTRL_FORCE_LAMDA_CYC_HEAT", CS%lam_cyc_heat, &
                 "A constant of proportionality between SST anomalies \n"//&
                 "and cyclical controlling heat fluxes", "W m-2 K-1", default=0.0)
  call get_param(param_file, mod, "CTRL_FORCE_LAMDA_CYC_PREC", CS%lam_cyc_prec, &
                 "A constant of proportionality between SSS anomalies \n"//&
                 "(normalised by mean SSS) and cyclical controlling \n"//&
                 "precipitation.", "kg m-2", default=0.0)

  CS%Len2 = smooth_len**2

! ### REPLACE THIS WITH ANY DIAGNOSTICS FROM THIS MODULE.
!  CS%id_taux = register_diag_field('ocean_model', 'taux', diag%axesu1, Time, &
!        'Zonal Wind Stress', 'Pascal')

end subroutine controlled_forcing_init

!> Clean up this modules control structure.
subroutine controlled_forcing_end(CS)
  type(ctrl_forcing_CS),    pointer :: CS !< A pointer to the control structure
                                          !! returned by a previous call to
                                          !! controlled_forcing_init, it will be
                                          !! deallocated here.

  if (associated(CS)) then
    if (associated(CS%heat_0))       deallocate(CS%heat_0)
    if (associated(CS%precip_0))     deallocate(CS%precip_0)
    if (associated(CS%heat_cyc))     deallocate(CS%heat_cyc)
    if (associated(CS%precip_cyc))   deallocate(CS%precip_cyc)
    if (associated(CS%avg_SST_anom)) deallocate(CS%avg_SST_anom)
    if (associated(CS%avg_SSS_anom)) deallocate(CS%avg_SSS_anom)
    if (associated(CS%avg_SSS))      deallocate(CS%avg_SSS)

    deallocate(CS)
  endif
  CS => NULL()

end subroutine controlled_forcing_end

!> \class MOM_controlled_forcing
!!                                                                     *
!!  By Robert Hallberg, July 2011                                      *
!!                                                                     *
!!    This program contains the subroutines that use control-theory    *
!!  to adjust the surface heat flux and precipitation, based on the    *
!!  time-mean or periodically (seasonally) varying anomalies from the  *
!!  observed state.  The techniques behind this are described in       *
!!  Hallberg and Adcroft (2011, in prep.).                             *
!!                                                                     *
!!  Macros written all in capital letters are defined in MOM_memory.h. *
!!                                                                     *
!!     A small fragment of the grid is shown below:                    *
!!                                                                     *
!!    j+1  x ^ x ^ x   At x:  q                                        *
!!    j+1  > o > o >   At ^:  v, tauy                                  *
!!    j    x ^ x ^ x   At >:  u, taux                                  *
!!    j    > o > o >   At o:  h, fluxes.                               *
!!    j-1  x ^ x ^ x                                                   *
!!        i-1  i  i+1  At x & ^:                                       *
!!           i  i+1    At > & o:                                       *
!!                                                                     *
!!  The boundaries always run through q grid points (x).               *
end module MOM_controlled_forcing
