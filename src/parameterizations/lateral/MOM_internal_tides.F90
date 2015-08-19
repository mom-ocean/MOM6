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

use MOM_spatial_means, only : global_area_mean ! BDM
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_axis_init
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_diag_mediator, only : axesType, defineAxes
use MOM_domains, only : AGRID, To_South, To_West, To_All
use MOM_domains, only : create_group_pass, do_group_pass, pass_var
use MOM_domains, only : group_pass_type, start_group_pass, complete_group_pass
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : slasher, vardesc
use MOM_restart, only : register_restart_field, MOM_restart_CS, restart_init, save_restart
use MOM_time_manager, only : time_type, operator(+), operator(/), operator(-)
use MOM_time_manager, only : get_time, get_date, set_time, set_date 
use MOM_time_manager, only : time_type_to_real
use MOM_variables, only : surface
use fms_mod, only : read_data
!   Forcing is a structure containing pointers to the forcing fields
! which may be used to drive MOM.  All fluxes are positive downward.
!   Surface is a structure containing pointers to various fields that
! may be used describe the surface state of MOM.

!use, intrinsic :: IEEE_ARITHMETIC

implicit none ; private

#include <MOM_memory.h>

public propagate_int_tide, register_int_tide_restarts
public internal_tides_init, internal_tides_end
public sum_itidal_lowmode_loss

type, public :: int_tide_CS ; private
  logical :: do_int_tides    ! If true, use the internal tide code.
  integer :: nFreq = 0
  integer :: nMode = 1
  integer :: nAngle = 24
  integer :: energized_angle = -1
  logical :: corner_adv      ! If true, use a corner advection rather than PPM.
  logical :: upwind_1st      ! If true, use a first-order upwind scheme.
  logical :: simple_2nd      ! If true, use a simple second order (arithmetic
                             ! mean) interpolation of the edge values instead
                             ! of the higher order interpolation.
  logical :: vol_CFL         ! If true, use the ratio of the open face lengths
                             ! to the tracer cell areas when estimating CFL
                             ! numbers.  Without aggress_adjust, the default is
                             ! false; it is always true with.
  logical :: use_PPMang      ! If true, use PPM for advection of energy in
			     ! angular space. (BDM)

  real, allocatable, dimension(:,:) :: refl_angle  
                        ! local coastline/ridge/shelf angles read from file (BDM)
                        ! (could be in G control structure)
  real, allocatable, dimension(:,:) :: refl_pref  
                        ! partial reflection coeff for each ``coast cell" (BDM)
                        ! (could be in G control structure)
  logical, allocatable, dimension(:,:) :: refl_pref_logical  
                        ! true if reflecting cell with partial reflection (BDM)
                        ! (could be in G control structure)
  logical, allocatable, dimension(:,:) :: refl_dbl  
                        ! identifies reflection cells where double reflection 
                        ! is possible (i.e. ridge cells) (BDM)
                        ! (could be in G control structure)   
  real, allocatable, dimension(:,:) :: TKE_itidal_loss_fixed  
                        ! fixed part of the energy lost due to small-scale drag
                        ! [kg m-2] here; will be multiplied by N and En to get 
                        ! into [W m-2] (BDM)
  real, allocatable, dimension(:,:,:,:,:) :: TKE_itidal_loss  
                        ! energy lost due to small-scale drag [W m-2] (BDM)
  real :: q_itides      ! fraction of local dissipation (nondimensional) (BDM)
  real :: En_sum        ! global sum of energy for use in debugging (BDM)    
  type(time_type),pointer    :: Time  
                        ! The current model time (BDM) 
  character(len=200) :: inputdir 
                        ! directory to look for coastline angle file (BDM)
  type(MOM_restart_CS), pointer :: restart_CSp => NULL()
                        ! for restart (BDM)
  character(len=200) :: restart_dir 
                        ! directory to write to for restart files (BDM)
  real :: decay_rate    ! A constant rate at which internal tide energy is
                        ! lost to the interior ocean internal wave field.
  real :: cdrag         ! The bottom drag coefficient for MEKE (non-dim).
  logical :: apply_drag ! If true, apply a quadratic bottom drag as a sink.
  real, dimension(:,:,:,:,:), pointer :: &
    En                  ! The internal wave energy density as a function of
                        ! (i,j,angle,frequency,mode)
  real, dimension(:,:,:), pointer :: &
    En_restart          ! The internal wave energy density as a function of
                        ! (i,j,angle); temporary for restart (BDM)
  real, allocatable, dimension(:) :: &
    frequency           ! The frequency of each band.

  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                        ! timing of diagnostic output.
  integer :: id_tot_En = -1, id_itide_drag = -1
  integer :: id_refl_pref = -1, id_refl_ang = -1, id_land_mask = -1 !(BDM)
  integer :: id_dx_Cv = -1, id_dy_Cu = -1 !(BDM)
  integer :: id_TKE_itidal_input = -1 !(BDM)
  integer, allocatable, dimension(:,:) :: id_En_mode, id_En_ang_mode
  integer, allocatable, dimension(:,:) :: id_TKE_loss_mode !(BDM)
  integer, allocatable, dimension(:,:) :: id_TKE_loss_ang_mode !(BDM)
  
end type int_tide_CS

type :: loop_bounds_type ; private
  integer :: ish, ieh, jsh, jeh
end type loop_bounds_type

contains


subroutine propagate_int_tide(cg1, TKE_itidal_input, vel_btTide, N2_bot, dt, G, CS)
  real, dimension(NIMEM_,NJMEM_), intent(in) :: cg1, TKE_itidal_input
  real, dimension(NIMEM_,NJMEM_), intent(in) :: vel_btTide, N2_bot
  real,                  intent(in)    :: dt
  type(ocean_grid_type), intent(inout) :: G
  type(int_tide_CS), pointer       :: CS
  ! This subroutine calls any of the other subroutines in this file
  ! that are needed to specify the current surface forcing fields.
  !
  ! Arguments: cg1 - The first mode internal gravity wave speed, in m s-1.
  !  (in)      TKE_itidal_input - The energy input to the internal waves, in W m-2.
  !  (in)      vel_btTide - Barotropic velocity read from file, in m s-1
  !  (in)      N2_bot - Squared near-bottom buoyancy frequency, in s-2
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
  integer :: isd_g, jsd_g         ! start indices on data domain but referenced 
                                  ! to global indexing (for debuggin-BDM)
  integer :: id_g, jd_g           ! global (decomp-invar) indices
                                  ! (for debuggin-BDM)
  type(group_pass_type), save :: pass_test, pass_En

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nAngle = CS%NAngle
  I_rho0 = 1.0 / G%Rho0
  
  isd_g = G%isd_global ; jsd_g = G%jsd_global ! for debugging (BDM)

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
        CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) + & 
                            dt*frac_per_sector*(1-CS%q_itides)*TKE_itidal_input(i,j)
    enddo ; enddo ; enddo ; enddo ; enddo
  elseif (CS%energized_angle <= CS%nAngle) then
    frac_per_sector = 1.0 / real(CS%nMode * CS%nFreq)
    a = CS%energized_angle
    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do j=js,je ; do i=is,ie
      f2 = 0.25*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                 (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J)**2))
      if (CS%frequency(fr)**2 > f2) &
        CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) + & 
                            dt*frac_per_sector**(1-CS%q_itides)*TKE_itidal_input(i,j)
    enddo ; enddo ; enddo ; enddo

    !### Delete this later.
    !do m=1,CS%nMode ; do j=jsd,jed ; do i=isd,ied
    !  c1(i,j,m) = 1.0 / real(m)
    !enddo ; enddo ; enddo
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
    call refract(CS%En(:,:,:,fr,m), c1(:,:,m), CS%frequency(fr), 0.5*dt, G, CS%nAngle, CS%use_PPMang)
  enddo ; enddo
  call do_group_pass(pass_En, G%domain)

  call complete_group_pass(pass_test, G%domain)

  ! Rotate points in the halos as necessary.
  call correct_halo_rotation(CS%En, test, G, CS%nAngle)
  
  ! Propagate the waves.
  do m=1,CS%NMode ; do fr=1,CS%Nfreq
    call propagate(CS%En(:,:,:,fr,m), c1(:,:,m), CS%frequency(fr), dt, G, CS, CS%NAngle)
  enddo ; enddo
  
  ! Test if energy has passed coast for debugging only; delete later(BDM)
  !do j=js,je
  !  jd_g = jsd_g + j - 1
  !  do i=is,ie 
  !    id_g = isd_g + i - 1
  !    if (id_g == 106 .and. jd_g == 55 ) then
        !print *, 'After propagation:'
        !print *, 'En_O  =', CS%En(i,j,:,1,1), 'refl_angle=', CS%refl_angle(i,j) 
        !print *, 'En_W =', CS%En(i-1,j,:,1,1), 'refl_angle=', CS%refl_angle(i-1,j)
        !print *, 'En_NW =', CS%En(i-1,j+1,:,1,1), 'refl_angle=', CS%refl_angle(i-1,j+1)
        !print *, 'En_N =', CS%En(i,j+1,:,1,1), 'refl_angle=', CS%refl_angle(i,j+1)
        !print *, 'En_NE =', CS%En(i+1,j+1,:,1,1), 'refl_angle=', CS%refl_angle(i+1,j+1)
        !print *, 'En_E =', CS%En(i+1,j,:,1,1), 'refl_angle=', CS%refl_angle(i+1,j)        
 !     endif
 !   enddo
 ! enddo

  ! Apply the other half of the refraction.
  do m=1,CS%NMode ; do fr=1,CS%Nfreq
    call refract(CS%En(:,:,:,fr,m), c1(:,:,m), CS%frequency(fr), 0.5*dt, G, CS%NAngle, CS%use_PPMang)
  enddo ; enddo

  if (CS%apply_drag .or. (CS%id_tot_En > 0)) then
    tot_En(:,:) = 0.0
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
      do j=js,je ; do i=is,ie
        tot_En(i,j) = tot_En(i,j) + CS%En(i,j,a,fr,m)
      enddo ; enddo
    enddo ; enddo ; enddo
  endif

  ! Extract the energy for mixing (original).
  !if (CS%apply_drag) then
  !  do j=jsd,jed ; do i=isd,ied
  !    I_D_here = 1.0 / max(G%bathyT(i,j), 1.0)
  !    drag_scale(i,j) = CS%decay_rate + CS%cdrag * &
  !      sqrt(max(0.0, vel_btTide(i,j)**2 + tot_En(i,j) * I_rho0 * I_D_here)) * I_D_here
  !  enddo ; enddo
  !else
  !  do j=jsd,jed ; do i=isd,ied
  !    drag_scale(i,j) = CS%decay_rate
  !  enddo ; enddo
  !endif  
  !do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
  !  CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) / (1.0 + dt*drag_scale(i,j))
  !enddo ; enddo ; enddo ; enddo ; enddo
  
  ! Extract the energy for mixing due to scattering -BDM.
  if (CS%apply_drag) then
    call itidal_lowmode_loss(G, CS, N2_bot, CS%En, CS%TKE_itidal_loss_fixed, &
    CS%TKE_itidal_loss, dt)
  endif
  
  ! Check for energy conservation on computational domain (BDM)
  call sum_En(G,CS,CS%En(:,:,:,1,1),'prop_int_tide')

  if (query_averaging_enabled(CS%diag)) then
    ! Output two-dimensional diagnostistics
    if (CS%id_tot_En > 0) call post_data(CS%id_tot_En, tot_En, CS%diag)
    if (CS%id_itide_drag > 0) call post_data(CS%id_itide_drag, drag_scale, CS%diag)
    if (CS%id_refl_ang > 0) call post_data(CS%id_refl_ang, CS%refl_angle, CS%diag) !(BDM)
    if (CS%id_refl_pref > 0) call post_data(CS%id_refl_pref, CS%refl_pref, CS%diag) !(BDM)
    if (CS%id_dx_Cv > 0) call post_data(CS%id_dx_Cv, G%dx_Cv, CS%diag) !(BDM)
    if (CS%id_dy_Cu > 0) call post_data(CS%id_dy_Cu, G%dy_Cu, CS%diag) !(BDM)
    if (CS%id_TKE_itidal_input > 0) call post_data(CS%id_TKE_itidal_input, &
      TKE_itidal_input, CS%diag) !(BDM)
    if (CS%id_land_mask > 0) call post_data(CS%id_land_mask, G%mask2dT, CS%diag) !(BDM)
    
    ! Output 2-D energy density (summed over angles) for each freq and mode
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; if (CS%id_En_mode(fr,m) > 0) then
      tot_En(:,:) = 0.0
      do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
        tot_En(i,j) = tot_En(i,j) + CS%En(i,j,a,fr,m)
      enddo ; enddo ; enddo
      call post_data(CS%id_En_mode(fr,m), tot_En, CS%diag)
    endif ; enddo ; enddo    
    ! Output 3-D (i,j,a) energy density for each freq and mode
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; if (CS%id_En_ang_mode(fr,m) > 0) then
      call post_data(CS%id_En_ang_mode(fr,m), CS%En(:,:,:,fr,m) , CS%diag)
    endif ; enddo ; enddo
    
    ! Output 2-D energy loss (summed over angles) for each freq and mode
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; if (CS%id_TKE_loss_mode(fr,m) > 0) then
      tot_En(:,:) = 0.0
      do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
        tot_En(i,j) = tot_En(i,j) + CS%TKE_itidal_loss(i,j,a,fr,m)
      enddo ; enddo ; enddo
      call post_data(CS%id_TKE_loss_mode(fr,m), tot_En, CS%diag)
    endif ; enddo ; enddo
    ! Output 3-D (i,j,a) energy loss for each freq and mode
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; if (CS%id_TKE_loss_ang_mode(fr,m) > 0) then
      call post_data(CS%id_TKE_loss_ang_mode(fr,m), CS%TKE_itidal_loss(:,:,:,fr,m) , CS%diag)
    endif ; enddo ; enddo
    
  endif
     
end subroutine propagate_int_tide

subroutine sum_En(G, CS, En, label)
  type(ocean_grid_type),  intent(in)    :: G
  type(int_tide_CS), pointer            :: CS
  !real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle,CS%nFreq,CS%nMode), intent(in) :: En
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle), intent(in) :: En
  character(len=*), intent(in) :: label
     
  ! This subroutine checks for energy conservation on computational domain  
  integer :: m,fr,a
  real :: En_sum, tmpForSumming, En_sum_diff, En_sum_pdiff
  integer :: seconds
  real :: Isecs_per_day = 1.0 / 86400.0
  real :: days
    
  call get_time(CS%Time, seconds)
  days = real(seconds) * Isecs_per_day
  
  !do m=1,CS%NMode
  !  do fr=1,CS%Nfreq
      En_sum = 0.0;
      tmpForSumming = 0.0
      do a=1,CS%nAngle
        !tmpForSumming = global_area_mean(En(:,:,a,fr,m),G)*G%areaT_global
        tmpForSumming = global_area_mean(En(:,:,a),G)*G%areaT_global
        En_sum = En_sum + tmpForSumming
      enddo
      En_sum_diff = En_sum - CS%En_sum
      if (CS%En_sum .ne. 0.0) then
        En_sum_pdiff= (En_sum_diff/CS%En_sum)*100.0
      else
        En_sum_pdiff= 0.0;
      endif
      CS%En_sum = En_sum
      if (is_root_pe()) then
        print *, label,':','days =', days
        !print *, 'mode(',m,'), freq(',fr,'): 
        print *, 'En_sum=', En_sum
        print *, 'En_sum_diff=',    En_sum_diff
        print *, 'Percent change=', En_sum_pdiff, '%'
        !if (abs(En_sum_pdiff) > 1.0) then ; stop ; endif        
      endif
  !  enddo
  !enddo  
end subroutine sum_En

subroutine itidal_lowmode_loss(G, CS, Nb2, En, TKE_loss_fixed, TKE_loss, dt)
  type(ocean_grid_type),  intent(in)    :: G
  type(int_tide_CS), pointer            :: CS
  real, dimension(G%isd:G%ied,G%jsd:G%jed), intent(in) :: Nb2
  real, dimension(G%isd:G%ied,G%jsd:G%jed), intent(in) :: TKE_loss_fixed
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle,CS%nFreq,CS%nMode), intent(inout) :: En
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle,CS%nFreq,CS%nMode), intent(out)   :: TKE_loss
  real :: dt
     
  ! This subroutine calculates the energy lost from the propagating internal tide due to
  ! scattering over small-scale roughness along the lines of Jayne & St. Laurent (2001).
  !
  ! Arguments: 
  !  (in)      Nb2 - near-bottom stratification, in s-2.
  !  (inout)   En - energy density of the internal waves, in J m-2.
  !  (in)      TKE_loss_fixed - fixed part of energy loss, in kg m-2 (formerly "coef_1")
  !  (out)     TKE_loss - energy loss rate, in W m-2
  !  (in)      dt - time increment, in s 
    
  integer :: j,i,m,fr,a
  real :: Nb, modal_vel_bot, TKE_loss_coef
  
  do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      Nb = 0.0
      if (Nb2(i,j) > 0.0) then
        Nb = G%mask2dT(i,j)*sqrt(Nb2(i,j))
      endif
      ! Calculate bottom velocity for mode; 
      ! using a placeholder here - could be done in diabatic driver as is N2_bot.
      modal_vel_bot = 0.0
      if (En(i,j,a,fr,m) > 0.0) then
        modal_vel_bot = 0.001*sqrt(En(i,j,a,fr,m))
      endif
      ! Calculate TKE loss rate; units of [J m-2 = kg s-2] here
      TKE_loss_coef = TKE_loss_fixed(i,j) * modal_vel_bot**2
      ! Calculate TKE loss rate; units of [W m-2] here.
      TKE_loss(i,j,a,fr,m) = TKE_loss_coef * Nb
      ! Update energy remaining in original mode (this is an explicit calc for now)
      En(i,j,a,fr,m) = En(i,j,a,fr,m) - TKE_loss(i,j,a,fr,m)*dt
    enddo ; enddo
  enddo ; enddo ; enddo
end subroutine itidal_lowmode_loss


subroutine sum_itidal_lowmode_loss(i,j,CS,G,TKE_loss_sum)
  integer, intent(in)                :: i,j
  type(ocean_grid_type),  intent(in) :: G
  type(int_tide_CS), pointer         :: CS
  real, intent(out)                  :: TKE_loss_sum  
  ! This subroutine sums the energy lost from the propagating internal across all angles,
  ! modes, and frequencies for a given location.
  !
  ! Arguments: 
  !  (out)      TKE_loss_sum - total energy loss rate, in W m-2.  
  integer :: m, fr, a
  
  TKE_loss_sum = 0.0
  do a = 1,CS%nAngle ; do fr = 1,CS%nFreq ; do m = 1,CS%nMode
    TKE_loss_sum = TKE_loss_sum + CS%TKE_itidal_loss(i,j,a,fr,m)
  enddo ; enddo ; enddo  
  
end subroutine sum_itidal_lowmode_loss


subroutine refract(En, cg, freq, dt, G, NAngle, use_PPMang)
  type(ocean_grid_type),  intent(in)    :: G
  integer,                intent(in)    :: NAngle
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), intent(inout) :: En
  real, dimension(G%isd:G%ied,G%jsd:G%jed),        intent(in)    :: cg
  real,                   intent(in)    :: freq
  real,                   intent(in)    :: dt
  logical,  	            intent(in)    :: use_PPMang !BDM
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
  real, dimension(SZI_(G),SZJ_(G),1-stensil:NAngle+stensil) :: &
    CFL_ang ! MOVED TO HERE BY BDM
  real :: f2   ! The squared Coriolis parameter, in s-2.
  real :: favg ! The average Coriolis parameter at a point, in s-1. Added by BDM.
  real :: df2_dy, df2_dx  ! The x- and y- gradients of the squared Coriolis parameter, in s-2 m-1.
  real :: df_dy, df_dx  ! The x- and y- gradients of the Coriolis parameter, in s-1 m-1. Added by BDM.
  real :: dlnCg_dx  ! The x-gradient of the wave speed divided by itself in m-1.
  real :: dlnCg_dy  ! The y-gradient of the wave speed divided by itself in m-1.
  real :: Angle_size, dt_Angle_size, angle
  real :: Ifreq, Kmag2, I_Kmag !CFL_ang REMOVED BY BDM
  real, parameter :: cg_subRO = 1e-100
  integer :: is, ie, js, je, asd, aed, na
  integer :: i, j, a

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; na = size(En,3)
  asd = 1-stensil ; aed = NAngle+stensil

  Ifreq = 1.0 / freq

  Angle_size = (8.0*atan(1.0)) / (real(NAngle))
  dt_Angle_size = dt / Angle_size

  !do A=0,na
  do A=asd,aed ! BDM
    angle = (real(A) - 0.5) * Angle_size
    cos_angle(A) = cos(angle) ; sin_angle(A) = sin(angle)
  enddo

  !### There should also be refraction due to cg.grad(grid_orientation).
  CFL_ang(:,:,:) = 0.0; ! INITIALIZING AS ZEROS (ADDED BY BDM)
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
      favg = 0.25*((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1)) + &
                 (G%CoriolisBu(I,J-1) + G%CoriolisBu(I-1,J))) ! added by BDM
      df2_dx = 0.5*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I,J-1)**2) - &
                    (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I-1,J-1)**2)) * &
               G%IdxT(i,j)
      df_dx = 0.5*((G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1)) - &
                    (G%CoriolisBu(I-1,J) + G%CoriolisBu(I-1,J-1))) * &
               G%IdxT(i,j) ! added by BDM
      dlnCg_dx = 0.5*( G%IdxCu(I,j) * (cg(i+1,j) - cg(i,j)) / &
                       (0.5*(cg(i+1,j) + cg(i,j)) + cg_subRO) + &
                       G%IdxCu(I-1,j) * (cg(i,j) - cg(i-1,j)) / &
                       (0.5*(cg(i,j) + cg(i-1,j)) + cg_subRO) )
      df2_dy = 0.5*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J)**2) - &
                    (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J-1)**2)) * &
               G%IdyT(i,j)
      df_dy = 0.5*((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J)) - &
                    (G%CoriolisBu(I,J-1) + G%CoriolisBu(I-1,J-1))) * &
               G%IdyT(i,j) ! added by BDM
      dlnCg_dy = 0.5*( G%IdyCv(i,J) * (cg(i,j+1) - cg(i,j)) / &
                       (0.5*(cg(i,j+1) + cg(i,j)) + cg_subRO) + &
                       G%IdyCv(i,J-1) * (cg(i,j) - cg(i,j-1)) / &
                       (0.5*(cg(i,j) + cg(i,j-1)) + cg_subRO) )
      Kmag2 = (freq**2 - f2) / (cg(i,j)**2 + cg_subRO**2)
      if (Kmag2 > 0.0) then
        I_Kmag = 1.0 / sqrt(Kmag2)
        !Dk_Dt_Kmag(i) = -Ifreq * (df2_dx + (freq**2 - f2) * dlnCg_dx) * I_Kmag
        !Dl_Dt_Kmag(i) = -Ifreq * (df2_dy + (freq**2 - f2) * dlnCg_dy) * I_Kmag
        Dk_Dt_Kmag(i) = -Ifreq * (favg*df_dx + (freq**2 - f2) * dlnCg_dx) * I_Kmag ! updated by BDM
        Dl_Dt_Kmag(i) = -Ifreq * (favg*df_dy + (freq**2 - f2) * dlnCg_dy) * I_Kmag ! updated by BDM
      else
        Dk_Dt_Kmag(i) = 0.0
        Dl_Dt_Kmag(i) = 0.0
      endif
    enddo

    ! Determine the energy fluxes in angular orientation space.
    !do A=0,na ; do i=is,ie
    do A=asd,aed ; do i=is,ie ! BDM
      CFL_ang(i,j,A) = (cos_angle(A) * Dl_Dt_Kmag(i) - sin_angle(A) * Dk_Dt_Kmag(i)) * &
                dt_Angle_size ! corrected by BDM
      if (abs(CFL_ang(i,j,A)) > 1.0) then 
        call MOM_error(WARNING, "refract: CFL exceeds 1.", .true.)
        if (CFL_ang(i,j,A) > 0.0) then ; CFL_ang(i,j,A) = 1.0 ; else ; CFL_ang(i,j,A) = -1.0 ; endif
      endif
    enddo; enddo
    
    !ADVECT IN ANGULAR SPACE: SIMPLE UPWIND AS BEFORE - BDM
    if(.not.use_PPMang) then
      do  A=0,na ; do i=is,ie
        if (CFL_ang(i,j,A) > 0.0) then
          Flux_E(i,A) = CFL_ang(i,j,A) * En2d(i,A)
        else
          Flux_E(i,A) = CFL_ang(i,j,A) * En2d(i,A+1)
        endif
      enddo; enddo
    else    
    !ADVECT IN ANGULAR SPACE: NEW PPM ADVECTION IN ANGULAR SPACE - BDM
      do i=is,ie
        call PPM_angular_advect(En2d(i,:),CFL_ang(i,j,:),Flux_E(i,:),NAngle,dt,stensil)
      enddo
    endif
      
  ! Update and copy back to En.
    do a=1,na ; do i=is,ie
      En(i,j,a) = En2d(i,a) + (Flux_E(i,A-1) - Flux_E(i,A))
    enddo ; enddo
  enddo ! j-loop
end subroutine refract


subroutine PPM_angular_advect(En2d, CFL_ang, Flux_En, NAngle, dt, halo_ang)
  integer,                     intent(in)    :: NAngle
  real,                        intent(in)    :: dt
  integer,                     intent(in)    :: halo_ang
  real, dimension(1-halo_ang:NAngle+halo_ang),   intent(in)    :: En2d
  real, dimension(1-halo_ang:NAngle+halo_ang),   intent(in)    :: CFL_ang
  real, dimension(0:NAngle),   intent(out)   :: Flux_En
  
  !   This subroutine calculates the 1-d flux for advection in angular space 
  ! using a monotonic piecewise parabolic scheme. Should be within i and j spatial
  ! loops - BDM
  real :: flux
  real :: u_ang
  real :: Angle_size
  real :: I_Angle_size
  real :: I_dt
  integer :: a
  real :: aR, aL, dMx, dMn, Ep, Ec, Em, dA, mA, a6

  I_dt = 1 / dt
  Angle_size = (8.0*atan(1.0)) / (real(NAngle))
  I_Angle_size = 1 / Angle_size
  Flux_En(:) = 0
  
  do A=0,NAngle
    u_ang = CFL_ang(A)*Angle_size*I_dt  
    if (u_ang >= 0.0) then
      ! Implementation of PPM-H3
      Ep = En2d(a+1)*I_Angle_size !MEAN ANGULAR ENERGY DENSITY FOR WEDGE (Jm-2/rad)  
      Ec = En2d(a)  *I_Angle_size !MEAN ANGULAR ENERGY DENSITY FOR WEDGE (Jm-2/rad) 
      Em = En2d(a-1)*I_Angle_size !MEAN ANGULAR ENERGY DENSITY FOR WEDGE (Jm-2/rad)
      aL = ( 5.*Ec + ( 2.*Em - Ep ) )/6. ! H3 estimate
      aL = max( min(Ec,Em), aL) ; aL = min( max(Ec,Em), aL) ! Bound
      aR = ( 5.*Ec + ( 2.*Ep - Em ) )/6. ! H3 estimate
      aR = max( min(Ec,Ep), aR) ; aR = min( max(Ec,Ep), aR) ! Bound
      dA = aR - aL ; mA = 0.5*( aR + aL )
      if ((Ep-Ec)*(Ec-Em) <= 0.) then
        aL = Ec ; aR = Ec ! PCM for local extremum
      elseif ( dA*(Ec-mA) > (dA*dA)/6. ) then
        aL = 3.*Ec - 2.*aR !?
      elseif ( dA*(Ec-mA) < - (dA*dA)/6. ) then
        aR = 3.*Ec - 2.*aL !?
      endif
      a6 = 6.*Ec - 3. * (aR + aL) ! Curvature
      ! CALCULATE FLUX RATE (Jm-2/s)
      flux = u_ang*( aR + 0.5 * CFL_ang(A) * ( ( aL - aR ) + a6 * ( 1. - 2./3. * CFL_ang(A) ) ) )
      !flux = u_ang*( aR - 0.5 * CFL_ang(A) * ( ( aR - aL ) - a6 * ( 1. - 2./3. * CFL_ang(A) ) ) )
      ! CALCULATE AMOUNT FLUXED (Jm-2)
      Flux_En(A) = dt * flux
      !Flux_En(A) = (dt * I_Angle_size) * flux
    else
      ! Implementation of PPM-H3
      Ep = En2d(a+2)*I_Angle_size !MEAN ANGULAR ENERGY DENSITY FOR WEDGE (Jm-2/rad)  
      Ec = En2d(a+1)*I_Angle_size !MEAN ANGULAR ENERGY DENSITY FOR WEDGE (Jm-2/rad) 
      Em = En2d(a)  *I_Angle_size !MEAN ANGULAR ENERGY DENSITY FOR WEDGE (Jm-2/rad)
      aL = ( 5.*Ec + ( 2.*Em - Ep ) )/6. ! H3 estimate
      aL = max( min(Ec,Em), aL) ; aL = min( max(Ec,Em), aL) ! Bound
      aR = ( 5.*Ec + ( 2.*Ep - Em ) )/6. ! H3 estimate
      aR = max( min(Ec,Ep), aR) ; aR = min( max(Ec,Ep), aR) ! Bound
      dA = aR - aL ; mA = 0.5*( aR + aL )
      if ((Ep-Ec)*(Ec-Em) <= 0.) then
        aL = Ec ; aR = Ec ! PCM for local extremum
      elseif ( dA*(Ec-mA) > (dA*dA)/6. ) then
        aL = 3.*Ec - 2.*aR
      elseif ( dA*(Ec-mA) < - (dA*dA)/6. ) then
        aR = 3.*Ec - 2.*aL
      endif
      a6 = 6.*Ec - 3. * (aR + aL) ! Curvature
      ! CALCULATE FLUX RATE (Jm-2/s)
      flux = u_ang*( aR + 0.5 * CFL_ang(A) * ( ( aL - aR ) + a6 * ( 1. - 2./3. * CFL_ang(A) ) ) )
      !flux = u_ang*( aL + 0.5 * CFL_ang(A) * ( ( aR - aL ) + a6 * ( 1. - 2./3. * CFL_ang(A) ) ) )
      ! CALCULATE AMOUNT FLUXED (Jm-2)
      Flux_En(A) = dt * flux
      !Flux_En(A) = (dt * I_Angle_size) * flux
    endif
  enddo
end subroutine PPM_angular_advect


subroutine propagate(En, cg, freq, dt, G, CS, NAngle)
  type(ocean_grid_type),  intent(inout) :: G
  integer,                intent(in)    :: NAngle
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), intent(inout) :: En
  real, dimension(G%isd:G%ied,G%jsd:G%jed),        intent(in)    :: cg
  real,                   intent(in)    :: freq
  real,                   intent(in)    :: dt
  type(int_tide_CS),      pointer       :: CS
  !  This subroutine does refraction on the internal waves at a single frequency.

  ! Arguments: En - the internal gravity wave energy density as a function of space
  !                 and angular resolution, in J m-2 radian-1.
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB) :: &
    speed  ! The magnitude of the group velocity at the q points for corner adv, in m s-1.
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
  
  ! Define loop bounds: Need extensions on j-loop so propagate_y 
  ! (done after propagate_x) will have updated values in the halo 
  ! for correct PPM reconstruction. Use if no teleporting and 
  ! no pass_var between propagate_x and propagate_y (BDM). 
  !jsh = js-3 ; jeh = je+3 ; ish = is ; ieh = ie  
  
  ! Define loop bounds: Need 1-pt extensions on loops because 
  ! teleporting eats up a halo point. Use if teleporting.
  ! Also requires pass_var before propagate_y (BDM).
  jsh = js-1 ; jeh = je+1 ; ish = is-1 ; ieh = ie+1
    
  Angle_size = (8.0*atan(1.0)) / real(NAngle)
  I_Angle_size = 1.0 / Angle_size
  
  if (CS%corner_adv) then
    ! IMPLEMENT CORNER ADVECTION IN HORIZONTAL--------------------
    ! FIND AVERAGE GROUP VELOCITY (SPEED) AT CELL CORNERS (BDM)
    ! Fix indexing here later
    speed(:,:) = 0;
    do J=jsh-1,jeh ; do I=ish-1,ieh
      f2 = G%CoriolisBu(I,J)**2
      speed(I,J) = 0.25*(cg(i,j) + cg(i+1,j) + cg(i+1,j+1) + cg(i,j+1)) * &
                     sqrt(max(freq2 - f2, 0.0)) * Ifreq
    enddo ; enddo  
    do a=1,na
      ! Apply the propagation WITH CORNER ADVECTION/FINITE VOLUME APPROACH (BDM).
      LB%jsh = js ; LB%jeh = je ; LB%ish = is ; LB%ieh = ie
      call propagate_corner_spread(En(:,:,a), a, NAngle, speed, dt, G, CS, LB)
    end do ! a-loop
  else 
    ! IMPLEMENT PPM ADVECTION IN HORIZONTAL-----------------------
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

    do j=jsh,jeh ; do I=ish-1,ieh
      f2 = 0.5*(G%CoriolisBu(I,J)**2 + G%CoriolisBu(I,J-1)**2)
      speed_x(I,j) = 0.5*(cg(i,j) + cg(i+1,j)) * G%mask2dCu(I,j) * &
                     sqrt(max(freq2 - f2, 0.0)) * Ifreq
    enddo ; enddo
    do J=jsh-1,jeh ; do i=ish,ieh
      f2 = 0.5*(G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J)**2)
      speed_y(i,J) = 0.5*(cg(i,j) + cg(i,j+1)) * G%mask2dCv(i,J) * &
                     sqrt(max(freq2 - f2, 0.0)) * Ifreq
    enddo ; enddo
    
    ! Apply propagation in x-direction (reflection included)
    LB%jsh = jsh ; LB%jeh = jeh ; LB%ish = ish ; LB%ieh = ieh
    call propagate_x(En(:,:,:), speed_x, Cgx_av(:), dCgx(:), dt, G, CS%nAngle, CS, LB)
    
    ! Check for energy conservation on computational domain (BDM)
    !call sum_En(G,CS,En(:,:,:),'post-propagate_x')
    
    ! Update halos
    call pass_var(En(:,:,:),G%domain)

    ! Apply propagation in y-direction (reflection included)
    ! LB%jsh = js ; LB%jeh = je ; LB%ish = is ; LB%ieh = ie ! Use if no teleport
    LB%jsh = jsh ; LB%jeh = jeh ; LB%ish = ish ; LB%ieh = ieh
    call propagate_y(En(:,:,:), speed_y, Cgy_av(:), dCgy(:), dt, G, CS%nAngle, CS, LB)
    
    ! Check for energy conservation on computational domain (BDM)
    !call sum_En(G,CS,En(:,:,:),'post-propagate_y')
 
  endif
end subroutine propagate


subroutine propagate_corner_spread(En, energized_wedge, NAngle, speed, dt, G, CS, LB)
  type(ocean_grid_type),  intent(in) :: G
  real, dimension(G%isd:G%ied,G%jsd:G%jed),   intent(inout) :: En
  real, dimension(G%IsdB:G%IedB,G%Jsd:G%Jed), intent(in)    :: speed
  integer,                intent(in)    :: energized_wedge
  integer,                intent(in)    :: NAngle
  real,                   intent(in)    :: dt
  type(int_tide_CS),      pointer       :: CS
  type(loop_bounds_type), intent(in)    :: LB
  ! Arguments: En - The energy density integrated over an angular band, in W m-2,
  !                 intent in/out.
  !  (in)      energized_wedge - index of current ray direction
  !  (in)      speed - The magnitude of the group velocity at the cell corner
  !                      points, in m s-1.
  !  (in)      dt - Time increment in s.
  !  (in)      G - The ocean's grid structure.
  !  (in)      CS - The control structure returned by a previous call to
  !                 continuity_PPM_init.
  !  (in)      LB - A structure with the active energy loop bounds.

  integer :: i, j, k, ish, ieh, jsh, jeh, m
  real :: TwoPi, Angle_size
  real :: energized_angle ! angle through center of current wedge
  real :: theta ! angle at edge of wedge
  real :: Nsubrays     ! number of sub-rays for averaging; 
                       ! count includes the two rays that bound the current wedge, 
                       ! i.e. those at -dtheta/2 and +dtheta/2 from energized angle
  real :: I_Nsubwedges ! inverse of number of sub-wedges 
  real :: cos_thetaDT, sin_thetaDT ! cos(theta)*dt, sin(theta)*dt
  real :: xNE,xNW,xSW,xSE,yNE,yNW,ySW,ySE ! corner point coordinates of advected fluid parcel
  real :: CFL_xNE,CFL_xNW,CFL_xSW,CFL_xSE,CFL_yNE,CFL_yNW,CFL_ySW,CFL_ySE,CFL_max
  real :: xN,xS,xE,xW,yN,yS,yE,yW ! intersection point coordinates of parcel edges and grid
  real :: xCrn,yCrn ! grid point contained within advected fluid parcel
  real :: xg,yg ! grid point of interest
  real :: slopeN,slopeW,slopeS,slopeE, bN,bW,bS,bE ! parameters defining parcel sides
  real :: aNE,aN,aNW,aW,aSW,aS,aSE,aE,aC ! sub-areas of advected parcel
  real :: a_total ! total area of advected parcel
  real :: a1,a2,a3,a4 ! areas used in calculating polygon areas (sub-areas) of advected parcel
  real, dimension(G%IsdB:G%IedB,G%Jsd:G%Jed) :: x,y ! coordinates of cell corners
  real, dimension(G%IsdB:G%IedB,G%Jsd:G%Jed) :: Idx,Idy ! inverse of dx,dy at cell corners
  real, dimension(G%IsdB:G%IedB,G%Jsd:G%Jed) :: dx,dy ! dx,dy at cell corners
  real, dimension(2) :: E_new ! energy in cell after advection for subray; set size here to 
                              ! define Nsubrays - this should be made an input option later!

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh
  TwoPi = (8.0*atan(1.0))
  Nsubrays = real(size(E_new))
  I_Nsubwedges = 1./(Nsubrays - 1)
  
  Angle_size = TwoPi / real(NAngle)
  energized_angle = Angle_size * real(energized_wedge - 1) ! for a=1 aligned with x-axis
  !energized_angle = Angle_size * real(energized_wedge - 1) + 2.0*Angle_size ! (match 297)
  !energized_angle = Angle_size * real(energized_wedge - 1) + 0.5*Angle_size ! (match 297)
  x = G%geoLonBu
  y = G%geoLatBu
  Idx = G%IdxBu; dx = G%dxBu
  Idy = G%IdyBu; dy = G%dyBu
 
  do j=jsh,jeh; do i=ish,ieh 
    do m=1,int(Nsubrays)
      theta = energized_angle - 0.5*Angle_size + real(m - 1)*Angle_size*I_Nsubwedges
      if (theta < 0.0) then
        theta = theta + TwoPi
      elseif (theta > TwoPi) then
        theta = theta - TwoPi
      endif
      cos_thetaDT = cos(theta)*dt
      sin_thetaDT = sin(theta)*dt

      ! corner point coordinates of advected fluid parcel ----------
      xg = x(I,J); yg = y(I,J)
      xNE = xg - speed(I,J)*cos_thetaDT
      yNE = yg - speed(I,J)*sin_thetaDT
      CFL_xNE = (xg-xNE)*Idx(I,J)
      CFL_yNE = (yg-yNE)*Idy(I,J)

      xg = x(I-1,J); yg = y(I-1,J)
      xNW = xg - speed(I-1,J)*cos_thetaDT
      yNW = yg - speed(I-1,J)*sin_thetaDT
      CFL_xNW = (xg-xNW)*Idx(I-1,J)
      CFL_yNW = (yg-yNW)*Idy(I-1,J)

      xg = x(I-1,J-1); yg = y(I-1,J-1)
      xSW = xg - speed(I-1,J-1)*cos_thetaDT
      ySW = yg - speed(I-1,J-1)*sin_thetaDT
      CFL_xSW = (xg-xSW)*Idx(I-1,J-1)
      CFL_ySW = (yg-ySW)*Idy(I-1,J-1)

      xg = x(I,J-1); yg = y(I,J-1)
      xSE = xg - speed(I,J-1)*cos_thetaDT
      ySE = yg - speed(I,J-1)*sin_thetaDT
      CFL_xSE = (xg-xSE)*Idx(I,J-1)
      CFL_ySE = (yg-ySE)*Idy(I,J-1)
      
      CFL_max = max(abs(CFL_xNE),abs(CFL_xNW),abs(CFL_xSW), &
                    abs(CFL_xSE),abs(CFL_yNE),abs(CFL_yNW), & 
                    abs(CFL_ySW),abs(CFL_ySE))
      if (CFL_max > 1.0) then 
        call MOM_error(WARNING, "propagate_corner_spread: CFL exceeds 1.", .true.)
      endif

      ! intersection point coordinates of parcel edges and cell edges ---
      if (0.0 <= theta .and. theta < 0.25*TwoPi) then
          xN = x(I-1,J-1)
          yW = y(I-1,J-1)
      elseif (0.25*TwoPi <= theta .and. theta < 0.5*TwoPi) then
          xN = x(I,J-1)
          yW = y(I,J-1)
      elseif (0.5*TwoPi <= theta .and. theta < 0.75*TwoPi) then 
          xN = x(I,J)
          yW = y(I,J)
      elseif (0.75*TwoPi <= theta .and. theta <= 1.00*TwoPi) then
          xN = x(I-1,J)
          yW = y(I-1,J)
      endif
      xS = xN
      yE = yW

      ! north intersection
      slopeN = (yNE - yNW)/(xNE - xNW)
      bN = -slopeN*xNE + yNE
      yN = slopeN*xN + bN
      ! west intersection
      if (xNW == xSW) then
        xW = xNW
      else
        slopeW = (yNW - ySW)/(xNW - xSW)
        bW = -slopeW*xNW + yNW
        xW = (yW - bW)/slopeW
      endif
      ! south intersection
      slopeS = (ySW - ySE)/(xSW - xSE)
      bS = -slopeS*xSW + ySW
      yS = slopeS*xS + bS
      ! east intersection
      if (xNE == xSE) then
        xE = xNE
      else
        slopeE = (ySE - yNE)/(xSE - xNE)
        bE = -slopeE*xSE + ySE
        xE = (yE - bE)/slopeE
      endif

      ! areas --------------------------------------------
      aNE = 0.0; aN = 0.0; aNW = 0.0; ! initialize areas
      aW = 0.0; aSW = 0.0; aS = 0.0; ! initialize areas 
      aSE = 0.0; aE = 0.0; aC = 0.0; ! initialize areas
      if (0.0 <= theta .and. theta < 0.25*TwoPi) then
          xCrn = x(I-1,J-1); yCrn = y(I-1,J-1);
          ! west area
          a1 = (yN - yCrn)*(0.5*(xN + xCrn))
          a2 = (yCrn - yW)*(0.5*(xCrn + xW))
          a3 = (yW - yNW)*(0.5*(xW + xNW))
          a4 = (yNW - yN)*(0.5*(xNW + xN))
          aW = a1 + a2 + a3 + a4
          ! southwest area
          a1 = (yCrn - yS)*(0.5*(xCrn + xS))
          a2 = (yS - ySW)*(0.5*(xS + xSW))
          a3 = (ySW - yW)*(0.5*(xSW + xW))
          a4 = (yW - yCrn)*(0.5*(xW + xCrn))
          aSW = a1 + a2 + a3 + a4
          ! south area
          a1 = (yE - ySE)*(0.5*(xE + xSE))
          a2 = (ySE - yS)*(0.5*(xSE + xS))
          a3 = (yS - yCrn)*(0.5*(xS + xCrn))
          a4 = (yCrn - yE)*(0.5*(xCrn + xE))
          aS = a1 + a2 + a3 + a4
          ! area within cell
          a1 = (yNE - yE)*(0.5*(xNE + xE))
          a2 = (yE - yCrn)*(0.5*(xE + xCrn))
          a3 = (yCrn - yN)*(0.5*(xCrn + xN))
          a4 = (yN - yNE)*(0.5*(xN + xNE))
          aC = a1 + a2 + a3 + a4
      elseif (0.25*TwoPi <= theta .and. theta < 0.5*TwoPi) then
          xCrn = x(I,J-1); yCrn = y(I,J-1);
          ! south area
          a1 = (yCrn - yS)*(0.5*(xCrn + xS))
          a2 = (yS - ySW)*(0.5*(xS + xSW))
          a3 = (ySW - yW)*(0.5*(xSW + xW))
          a4 = (yW - yCrn)*(0.5*(xW + xCrn))
          aS = a1 + a2 + a3 + a4        
          ! southeast area
          a1 = (yE - ySE)*(0.5*(xE + xSE))
          a2 = (ySE - yS)*(0.5*(xSE + xS))
          a3 = (yS - yCrn)*(0.5*(xS + xCrn))
          a4 = (yCrn - yE)*(0.5*(xCrn + xE))
          aSE = a1 + a2 + a3 + a4
          ! east area
          a1 = (yNE - yE)*(0.5*(xNE + xE))
          a2 = (yE - yCrn)*(0.5*(xE + xCrn))
          a3 = (yCrn - yN)*(0.5*(xCrn + xN))
          a4 = (yN - yNE)*(0.5*(xN + xNE))
          aE = a1 + a2 + a3 + a4
          ! area within cell
          a1 = (yN - yCrn)*(0.5*(xN + xCrn))
          a2 = (yCrn - yW)*(0.5*(xCrn + xW))
          a3 = (yW - yNW)*(0.5*(xW + xNW))
          a4 = (yNW - yN)*(0.5*(xNW + xN))
          aC = a1 + a2 + a3 + a4
      elseif (0.5*TwoPi <= theta .and. theta < 0.75*TwoPi) then
          xCrn = x(I,J); yCrn = y(I,J);
          ! east area
          a1 = (yE - ySE)*(0.5*(xE + xSE))
          a2 = (ySE - yS)*(0.5*(xSE + xS))
          a3 = (yS - yCrn)*(0.5*(xS + xCrn))
          a4 = (yCrn - yE)*(0.5*(xCrn + xE))
          aE = a1 + a2 + a3 + a4
          ! northeast area
          a1 = (yNE - yE)*(0.5*(xNE + xE))
          a2 = (yE - yCrn)*(0.5*(xE + xCrn))
          a3 = (yCrn - yN)*(0.5*(xCrn + xN))
          a4 = (yN - yNE)*(0.5*(xN + xNE))
          aNE = a1 + a2 + a3 + a4
          ! north area
          a1 = (yN - yCrn)*(0.5*(xN + xCrn))
          a2 = (yCrn - yW)*(0.5*(xCrn + xW))
          a3 = (yW - yNW)*(0.5*(xW + xNW))
          a4 = (yNW - yN)*(0.5*(xNW + xN))
          aN = a1 + a2 + a3 + a4
          ! area within cell
          a1 = (yCrn - yS)*(0.5*(xCrn + xS))
          a2 = (yS - ySW)*(0.5*(xS + xSW))
          a3 = (ySW - yW)*(0.5*(xSW + xW))
          a4 = (yW - yCrn)*(0.5*(xW + xCrn))
          aC = a1 + a2 + a3 + a4
      elseif (0.75*TwoPi <= theta .and. theta <= 1.00*TwoPi) then
          xCrn = x(I-1,J); yCrn = y(I-1,J);
          ! north area
          a1 = (yNE - yE)*(0.5*(xNE + xE))
          a2 = (yE - yCrn)*(0.5*(xE + xCrn))
          a3 = (yCrn - yN)*(0.5*(xCrn + xN))
          a4 = (yN - yNE)*(0.5*(xN + xNE))
          aN = a1 + a2 + a3 + a4
          ! northwest area
          a1 = (yN - yCrn)*(0.5*(xN + xCrn))
          a2 = (yCrn - yW)*(0.5*(xCrn + xW))
          a3 = (yW - yNW)*(0.5*(xW + xNW))
          a4 = (yNW - yN)*(0.5*(xNW + xN))
          aNW = a1 + a2 + a3 + a4;
          ! west area
          a1 = (yCrn - yS)*(0.5*(xCrn + xS))
          a2 = (yS - ySW)*(0.5*(xS + xSW))
          a3 = (ySW - yW)*(0.5*(xSW + xW))
          a4 = (yW - yCrn)*(0.5*(xW + xCrn))
          aW = a1 + a2 + a3 + a4
          ! area within cell
          a1 = (yE - ySE)*(0.5*(xE + xSE))
          a2 = (ySE - yS)*(0.5*(xSE + xS))
          a3 = (yS - yCrn)*(0.5*(xS + xCrn))
          a4 = (yCrn - yE)*(0.5*(xCrn + xE))
          aC = a1 + a2 + a3 + a4
      endif
     
      ! energy weighting ----------------------------------------
      a_total = aNE + aN + aNW + aW + aSW + aS + aSE + aE + aC
      E_new(m) = (aNE*En(i+1,j+1) + aN*En(i,j+1) + aNW*En(i-1,j+1) + &
                  aW*En(i-1,j) + aSW*En(i-1,j-1) + aS*En(i,j-1) + &
                  aSE*En(i+1,j-1) + aE*En(i+1,j) + aC*En(i,j)) / (dx(i,j)*dy(i,j))
    enddo ! m-loop
    ! update energy in cell
    En(i,j) = sum(E_new)/Nsubrays
  enddo; enddo
end subroutine propagate_corner_spread
    
subroutine propagate_x(En, speed_x, Cgx_av, dCgx, dt, G, Nangle, CS, LB)
  type(ocean_grid_type),  intent(in) :: G
  integer,                intent(in) :: NAngle
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle),   intent(inout) :: En
  real, dimension(G%IsdB:G%IedB,G%jsd:G%jed), intent(in)    :: speed_x
  real, dimension(Nangle),                    intent(in)    :: Cgx_av, dCgx
  real,                   intent(in)    :: dt
  type(int_tide_CS),      pointer       :: CS
  type(loop_bounds_type), intent(in)    :: LB
  ! Arguments: En - The energy density integrated over an angular band, in J m-2,
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
  !real, dimension(SZI_(G),SZJB_(G),Nangle) :: En_m, En_p 
  real, dimension(SZI_(G),SZJB_(G),Nangle) :: &
    Fdt_m, Fdt_p! Left and right energy fluxes, in J 
  integer :: i, j, k, ish, ieh, jsh, jeh, a
  
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh
  do a=1,Nangle
  ! This sets EnL and EnR.
    if (CS%upwind_1st) then
      do j=jsh,jeh ; do i=ish-1,ieh+1
        EnL(i,j) = En(i,j,a) ; EnR(i,j) = En(i,j,a)
      enddo ; enddo
    else
      call PPM_reconstruction_x(En(:,:,a), EnL, EnR, G, LB, simple_2nd=CS%simple_2nd)
    endif

    do j=jsh,jeh
      ! This is done once with single speed (GARDEN SPRINKLER EFFECT POSSIBLE)
      do I=ish-1,ieh
        cg_p(I) = speed_x(I,j) * (Cgx_av(a))
      enddo
      call zonal_flux_En(cg_p, En(:,j,a), EnL(:,j), EnR(:,j), flux1, &
                         dt, G, j, ish, ieh, CS%vol_CFL)
      do I=ish-1,ieh ; flux_x(I,j) = flux1(I); enddo      
    enddo
  
    do j=jsh,jeh ; do i=ish,ieh
      Fdt_m(i,j,a) = dt*flux_x(I-1,j) ! left face influx  (J)
      Fdt_p(i,j,a) = -dt*flux_x(I,j)  ! right face influx (J)
    enddo ; enddo
  
    ! test with old (take out later)
    !do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
    !  En(i,j,a) = En(i,j,a) - dt* G%IareaT(i,j) * (flux_x(I,j) - flux_x(I-1,j))
    !enddo ; enddo       
    
  enddo ! a-loop
  
  ! Only reflect newly arrived energy; existing energy in incident wedge
  ! is not reflected and will eventually propagate out of cell (BDM)
  ! (only reflects if En > 0)  
  call reflect(Fdt_m(:,:,:), Nangle, CS, G, LB)
  !call teleport(Fdt_m(:,:,:), Nangle, CS, G, LB)
  call reflect(Fdt_p(:,:,:), Nangle, CS, G, LB)
  !call teleport(Fdt_p(:,:,:), Nangle, CS, G, LB)
 
  ! Update reflected energy (Jm-2)
  do j=jsh,jeh ; do i=ish,ieh
    En(i,j,:) = En(i,j,:) + G%IareaT(i,j)*(Fdt_m(i,j,:) + Fdt_p(i,j,:))
  enddo ; enddo
  
end subroutine propagate_x
   
subroutine propagate_y(En, speed_y, Cgy_av, dCgy, dt, G, Nangle, CS, LB)
  type(ocean_grid_type),  intent(in) :: G
  integer,                intent(in) :: NAngle
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle),   intent(inout) :: En
  real, dimension(G%isd:G%ied,G%JsdB:G%JedB), intent(in)    :: speed_y
  real, dimension(Nangle),                    intent(in)    :: Cgy_av, dCgy
  real,                   intent(in)    :: dt
  type(int_tide_CS),      pointer       :: CS
  type(loop_bounds_type), intent(in)    :: LB
  ! Arguments: En - The energy density integrated over an angular band, in J m-2,
  !                 intent in/out.
  !  (in)      speed_y - The magnitude of the group velocity at the Cv
  !                      points, in m s-1.
  !  (in)      dt - Time increment in s.
  !  (in)      G - The ocean's grid structure.
  !  (in)      CS - The control structure returned by a previous call to
  !                 continuity_PPM_init.
  !  (in)      LB - A structure with the active energy loop bounds.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    EnL, EnR    ! South and north face energy densities, in J m-2.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    flux_y      ! The internal wave energy flux, in J s-1.
  real, dimension(SZI_(G)) :: &
    cg_p, cg_m, flux1, flux2
  !real, dimension(SZI_(G),SZJB_(G),Nangle) :: En_m, En_p
  real, dimension(SZI_(G),SZJB_(G),Nangle) :: &
    Fdt_m, Fdt_p! South and north energy fluxes, in J  
  integer :: i, j, k, ish, ieh, jsh, jeh, a

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh
  do a=1,Nangle
    ! This sets EnL and EnR.
    if (CS%upwind_1st) then
      do j=jsh-1,jeh+1 ; do i=ish,ieh
        EnL(i,j) = En(i,j,a) ; EnR(i,j) = En(i,j,a)
      enddo ; enddo
    else
      call PPM_reconstruction_y(En(:,:,a), EnL, EnR, G, LB, simple_2nd=CS%simple_2nd)
    endif

    do J=jsh-1,jeh
      !   This is done once with single speed (GARDEN SPRINKLER EFFECT POSSIBLE)
      do i=ish,ieh
        cg_p(i) = speed_y(i,J) * (Cgy_av(a))
      enddo
      call merid_flux_En(cg_p, En(:,:,a), EnL(:,:), EnR(:,:), flux1, &
                         dt, G, J, ish, ieh, CS%vol_CFL)
      do i=ish,ieh ; flux_y(i,J) = flux1(i); enddo
    enddo
    
    do j=jsh,jeh ; do i=ish,ieh
      Fdt_m(i,j,a) = dt*flux_y(i,J-1) ! south face influx (J)
      Fdt_p(i,j,a) = -dt*flux_y(i,J)  ! north face influx (J)
    enddo ; enddo   
    
    ! test with old (take out later)
    !do j=jsh,jeh ; do i=ish,ieh
    ! En(i,j,a) = En(i,j,a) - dt* G%IareaT(i,j) * (flux_y(i,J) - flux_y(i,J-1))
    !enddo ; enddo
         
  enddo ! a-loop
  
  ! Only reflect newly arrived energy; existing energy in incident wedge
  ! is not reflected and will eventually propagate out of cell (BDM)
  ! (only reflects if En > 0)    
  call reflect(Fdt_m(:,:,:), Nangle, CS, G, LB)
  !call teleport(Fdt_m(:,:,:), Nangle, CS, G, LB)
  call reflect(Fdt_p(:,:,:), Nangle, CS, G, LB)
  !call teleport(Fdt_p(:,:,:), Nangle, CS, G, LB)
  
  ! Update reflected energy (Jm-2)
  do j=jsh,jeh ; do i=ish,ieh
    En(i,j,:) = En(i,j,:) + G%IareaT(i,j)*(Fdt_m(i,j,:) + Fdt_p(i,j,:))
  enddo ; enddo

end subroutine propagate_y


subroutine zonal_flux_En(u, h, hL, hR, uh, dt, G, j, ish, ieh, vol_CFL)
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
  !  (in)      h - Energy density used to calculate the fluxes, in J m-2.
  !  (in)      hL, hR - Left- and right- Energy densities in the reconstruction, in J m-2.
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


subroutine merid_flux_En(v, h, hL, hR, vh, dt, G, J, ish, ieh, vol_CFL)
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
  !  (in)      h - Energy density used to calculate the fluxes, in J m-2.
  !  (in)      hL, hR - Left- and right- Energy densities in the reconstruction, in J m-2.
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


subroutine reflect(En, NAngle, CS, G, LB)
  type(ocean_grid_type),  intent(in)    :: G
  integer,                intent(in)    :: NAngle
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), intent(inout) :: En  
  type(int_tide_CS),      pointer       :: CS
  type(loop_bounds_type), intent(in)    :: LB

  !  This subroutine does reflection of the internal waves at a single frequency.
  
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: angle_c  
                                           ! angle of boudary wrt equator
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: part_refl  
                                           ! fraction of wave energy reflected
                                           ! values should collocate with angle_c
  logical, dimension(G%isd:G%ied,G%jsd:G%jed) :: ridge  
                                           ! tags of cells with double reflection

  real    :: TwoPi                         ! 2*pi
  real    :: Angle_size                    ! size of beam wedge (rad)
  real    :: angle_wall                    ! angle of coast/ridge/shelf wrt equator
  real, dimension(1:NAngle) :: angle_i     ! angle of incident ray wrt equator
  real    :: angle_r                       ! angle of reflected ray wrt equator
  real, dimension(1:Nangle) :: En_reflected  
  integer :: i, j, a, a_r, na  
  !integer :: isd, ied, jsd, jed   ! start and end local indices on data domain
  !                                ! (values include halos)
  integer :: isc, iec, jsc, jec   ! start and end local indices on PE
                                  ! (values exclude halos)
  integer :: ish, ieh, jsh, jeh   ! start and end local indices on data domain
                                  ! leaving out outdated halo points (march in)
  integer :: isd_g, jsd_g         ! start indices on data domain but referenced 
                                  ! to global indexing
  integer :: id_g, jd_g           ! global (decomp-invar) indices
  
  !isd = G%isd  ; ied = G%ied  ; jsd = G%jsd  ; jed = G%jed
  isc = G%isc  ; iec = G%iec  ; jsc = G%jsc  ; jec = G%jec
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh
  isd_g = G%isd_global ;        jsd_g = G%jsd_global
  
  TwoPi = 8.0*atan(1.0);
  Angle_size = TwoPi / (real(NAngle))
  
  do a=1,NAngle
    ! These are the angles at the cell centers  
    ! (should do this elsewhere since doesn't change with time)
    angle_i(a) = Angle_size * real(a - 1) ! for a=1 aligned with x-axis
  enddo
  
  angle_c   = CS%refl_angle
  part_refl = CS%refl_pref
  ridge     = CS%refl_dbl
  En_reflected(:) = 0.0

  !do j=jsc-1,jec+1
  do j=jsh,jeh
    jd_g = jsd_g + j - 1
    !do i=isc-1,iec+1 
    do i=ish,ieh
      id_g = isd_g + i - 1
      ! redistribute energy in angular space if ray will hit boundary 
      ! i.e., if energy is in a reflecting cell
      if (.not. isnan(angle_c(i,j))) then
        do a=1,NAngle
          if (En(i,j,a) > 0.0) then
            ! if ray is incident, keep specified boundary angle
            if (sin(angle_i(a) - angle_c(i,j)) >= 0.0) then
              angle_wall = angle_c(i,j)
            ! if ray is not incident but in ridge cell, use complementary angle
            elseif (ridge(i,j)) then
              angle_wall = angle_c(i,j) + 0.5*TwoPi
              if (angle_wall > TwoPi) then
                angle_wall = angle_wall - TwoPi*floor(abs(angle_wall)/TwoPi)
              elseif (angle_wall < 0.0) then
                angle_wall = angle_wall + TwoPi*ceiling(abs(angle_wall)/TwoPi)
              endif  
            ! if ray is not incident and not in a ridge cell, keep specified angle
            else
              angle_wall = angle_c(i,j)   
            endif 
            ! do reflection
            if (sin(angle_i(a) - angle_wall) >= 0.0) then
              angle_r = 2.0*angle_wall - angle_i(a)
              if (angle_r > TwoPi) then
                angle_r = angle_r - TwoPi*floor(abs(angle_r)/TwoPi)
              elseif (angle_r < 0.0) then
                angle_r = angle_r + TwoPi*ceiling(abs(angle_r)/TwoPi)
              endif
              a_r = nint(angle_r/Angle_size) + 1
              do while (a_r > Nangle) ; a_r = a_r - Nangle ; enddo
              if (a .ne. a_r) then
                En_reflected(a_r) = part_refl(i,j)*En(i,j,a)
                En(i,j,a)   = (1.0-part_refl(i,j))*En(i,j,a)
              endif
            endif
          endif
        enddo ! a-loop 
        En(i,j,:) = En(i,j,:) + En_reflected(:) 
        En_reflected(:) = 0.0
      endif
    enddo ! i-loop
  enddo ! j-loop
 
  ! Check to make sure no energy gets onto land (only run for debugging)
  do j=jsc,jec
    jd_g = jsd_g + j - 1
    do i=isc,iec
      id_g = isd_g + i - 1
      do a=1,NAngle 
        if (En(i,j,a) > 0.001 .and. G%mask2dT(i,j) == 0) then
          print *, 'En=', En(i,j,a), 'a=', a, 'ig_g=',id_g, 'jg_g=',jd_g
          !stop 'Energy detected out of bounds!' 
        endif
      enddo ! a-loop
    enddo ! i-loop
  enddo ! j-loop

end subroutine reflect

subroutine teleport(En, NAngle, CS, G, LB)
  type(ocean_grid_type),  intent(in)    :: G
  integer,                intent(in)    :: NAngle
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), intent(inout) :: En
  type(int_tide_CS),      pointer       :: CS
  type(loop_bounds_type), intent(in)    :: LB
    
  ! This subroutine moves energy across lines of partial reflection to prevent
  ! reflection of energy that is supposed to get across.

  real, dimension(G%isd:G%ied,G%jsd:G%jed)    :: angle_c  
                                              ! angle of boudary wrt equator
  real, dimension(G%isd:G%ied,G%jsd:G%jed)    :: part_refl  
                                              ! fraction of wave energy reflected
                                              ! values should collocate with angle_c
  logical, dimension(G%isd:G%ied,G%jsd:G%jed) :: pref_cell 
                                              ! flag for partial reflection
  logical, dimension(G%isd:G%ied,G%jsd:G%jed) :: ridge  
                                           ! tags of cells with double reflection
  real                        :: TwoPi      ! size of beam wedge (rad)
  real                        :: Angle_size ! size of beam wedge (rad)
  real, dimension(1:NAngle)   :: angle_i    ! angle of incident ray wrt equator
  real, dimension(1:NAngle)   :: cos_angle, sin_angle
  real                        :: En_tele    ! energy to be "teleported"
  integer :: i, j, a
  !integer :: isd, ied, jsd, jed    ! start and end local indices on data domain
  !                                 ! (values include halos)
  !integer :: isc, iec, jsc, jec    ! start and end local indices on PE
  !                                 ! (values exclude halos)
  integer :: ish, ieh, jsh, jeh     ! start and end local indices on data domain
                                    ! leaving out outdated halo points (march in)
  integer :: isd_g, jsd_g           ! start indices on data domain but referenced 
                                    ! to global indexing
  integer :: id_g, jd_g             ! global (decomp-invar) indices
  integer :: jos, ios               ! offsets
  real    :: cos_normal, sin_normal, angle_wall
                                    ! cos/sin of cross-ridge normal, ridge angle
  
  !isd = G%isd  ; ied = G%ied  ; jsd = G%jsd  ; jed = G%jed
  !isc = G%isc  ; iec = G%iec  ; jsc = G%jsc  ; jec = G%jec
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh
  isd_g = G%isd_global ;        jsd_g = G%jsd_global ! for debugging (BDM)
  
  TwoPi = 8.0*atan(1.0)
  Angle_size = TwoPi / (real(NAngle))
  
  do a=1,Nangle
    ! These are the angles at the cell centers  
    ! (should do this elsewhere since doesn't change with time)
    angle_i(a) = Angle_size * real(a - 1) ! for a=1 aligned with x-axis
    cos_angle(a) = cos(angle_i(a)) ; sin_angle(a) = sin(angle_i(a))
  enddo
  
  angle_c   = CS%refl_angle
  part_refl = CS%refl_pref
  pref_cell = CS%refl_pref_logical
  ridge     = CS%refl_dbl

  do j=jsh,jeh
    jd_g = jsd_g + j - 1
    do i=ish,ieh
      id_g = isd_g + i - 1
      if (pref_cell(i,j)) then
        do a=1,Nangle 
          if (En(i,j,a) > 0) then
            ! if ray is incident, keep specified boundary angle
            if (sin(angle_i(a) - angle_c(i,j)) >= 0.0) then
              angle_wall = angle_c(i,j)
            ! if ray is not incident but in ridge cell, use complementary angle
            elseif (ridge(i,j)) then
              angle_wall = angle_c(i,j) + 0.5*TwoPi
            ! if ray is not incident and not in a ridge cell, keep specified angle          
            else
              angle_wall = angle_c(i,j)            
            endif 
            ! teleport if incident
            if (sin(angle_i(a) - angle_wall) >= 0.0) then
              En_tele = En(i,j,a)
              cos_normal = cos(angle_wall + 0.25*TwoPi)
              sin_normal = sin(angle_wall + 0.25*TwoPi)
              ! find preferred zonal offset based on shelf/ridge angle
              ios = int(sign(1.,cos_normal))
              ! find preferred meridional offset based on shelf/ridge angle
              jos = int(sign(1.,sin_normal))
              ! find receptive ocean cell in direction of offset
              if (.not. pref_cell(i+ios,j+jos)) then
                En(i,j,a) = En(i,j,a) - En_tele
                En(i+ios,j+jos,a) = En(i+ios,j+jos,a) + En_tele
              else
                call MOM_error(WARNING, "teleport: no receptive ocean cell", .true.)
                print *, 'idg=',id_g,'jd_g=',jd_g,'a=',a
                stop
              endif
            endif ! incidence check
          endif ! energy check
        enddo ! a-loop
      endif ! pref check
    enddo ! i-loop
  enddo ! j-loop
   
end subroutine teleport

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
  type(ocean_grid_type), intent(inout) :: G
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
  logical :: use_int_tides
  type(vardesc) :: vd
  integer :: num_freq, num_angle , num_mode, period_1
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
  
  allocate(CS)
  
  num_angle = 24
  call read_param(param_file, "INTERNAL_TIDE_ANGLES", num_angle)
  allocate(CS%En_restart(isd:ied, jsd:jed, num_angle))
  CS%En_restart(:,:,:) = 0.0
  
  vd = vardesc("En_restart", &
    "The internal wave energy density as a function of (i,j,angle,frequency,mode)", &
    'h','1','1',"J m-2")
  call register_restart_field(CS%En_restart, vd, .false., restart_CS)  
  
  !--------------------check----------------------------------------------
  if (is_root_pe()) then
    print *,'register_int_tide_restarts: CS and CS%En_restart allocated!'
    print *,'register_int_tide_restarts: CS%En_restart registered!'
    print *,'register_int_tide_restarts: done!'
  endif
  !-----------------------------------------------------------------------
  
end subroutine register_int_tide_restarts


subroutine internal_tides_init(Time, G, param_file, diag, CS)
  !type(time_type),           intent(in) :: Time
  type(time_type), target,   intent(in) :: Time
  type(ocean_grid_type),     intent(inout) :: G
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
  real, allocatable, dimension(:,:) :: h2 ! (BDM)
  real, allocatable :: ridge_temp(:,:) ! intermediate array (BDM)
  logical :: use_int_tides, use_temperature
  integer :: num_angle, num_freq, num_mode, m, fr, period_1
  real :: kappa_itides, kappa_h2_factor ! (BDM)
  integer :: isd, ied, jsd, jed, a, id_ang, i, j
  type(axesType) :: axes_ang 
  ! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_internal_tides" ! This module's name.
  character(len=16), dimension(8) :: freq_name
  character(len=40)  :: var_name
  character(len=160) :: var_descript
  character(len=200) :: filename ! (BDM)
  character(len=200) :: refl_angle_file, land_mask_file, refl_pref_file ! (BDM)
  character(len=200) :: refl_dbl_file ! (BDM)
  character(len=200) :: dy_Cu_file, dx_Cv_file ! (BDM)
  character(len=200) :: h2_file ! (BDM)
  
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  use_int_tides = .false.
  call read_param(param_file, "INTERNAL_TIDES", use_int_tides)
  CS%do_int_tides = use_int_tides
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
  CS%nFreq = num_freq ; CS%nAngle = num_angle ; CS%nMode = num_mode

  allocate(CS%En(isd:ied, jsd:jed, num_angle, num_freq, num_mode))
  CS%En(:,:,:,:,:) = 0.0    
  CS%En(:,:,:,1,1) = CS%En_restart(:,:,:) ! added here as work-around (BDM)
  
  !print *, 'register_int_tide_init: sum(En_restart)=', sum(CS%En_restart)
      
  !--------------------check----------------------------------------------
  if (is_root_pe()) then
    print *,'internal_tides_init: CS%En_restart added to CS%En!' !BDM
  endif
  !-----------------------------------------------------------------------

  allocate(CS%frequency(num_freq))
  call read_param(param_file, "FIRST_MODE_PERIOD", period_1); ! ADDED BDM
  do a=1,num_freq  
    CS%frequency(a) = (8.0*atan(1.0) * (real(a)) / period_1) ! ADDED BDM
  enddo

  ! Read all relevant parameters and write them to the model log.
  CS%Time => Time ! direct a pointer to the current model time target (BDM)
  
  call get_param(param_file, mod, "INPUTDIR", CS%inputdir, default=".")
  CS%inputdir = slasher(CS%inputdir) ! (BDM)
  
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
  !call get_param(param_file, mod, "MIN_ZBOT_ITIDES", CS%min_zbot_itides, &
  !              "Turn off internal tidal dissipation when the total \n"//&
  !              "ocean depth is less than this value.", units="m", default=0.0) !BDM
  call get_param(param_file, mod, "INTERNAL_TIDE_VOLUME_BASED_CFL", CS%vol_CFL, &
                 "If true, use the ratio of the open face lengths to the \n"//&
                 "tracer cell areas when estimating CFL numbers in the \n"//&
                 "internal tide code.", default=.false.)
  call get_param(param_file, mod, "INTERNAL_TIDE_CORNER_ADVECT", CS%corner_adv, &
                 "If true, internal tide ray-tracing advection uses a \n"//&
                 " corner-advection scheme rather than PPM.\n", default=.false.)
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
  call get_param(param_file, mod, "USE_PPM_ANGULAR", CS%use_PPMang, &
                 "If true, use PPM for advection of energy in angular  \n"//&
                 "space.", default=.false.)
  call get_param(param_file, mod, "GAMMA_ITIDES", CS%q_itides, &
                 "The fraction of the internal tidal energy that is \n"//&
                 "dissipated locally with INT_TIDE_DISSIPATION.  \n"//&
                 "THIS NAME COULD BE BETTER.", &
                 units="nondim", default=0.3333) !(BDM)
                 
  ! Compute the fixed part of the bottom drag loss from baroclinic modes (BDM)
  allocate(h2(isd:ied,jsd:jed)) ; h2(:,:) = 0.0
  allocate(CS%TKE_itidal_loss_fixed(isd:ied,jsd:jed)) 
  CS%TKE_itidal_loss_fixed = 0.0
  allocate(CS%TKE_itidal_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode))
  CS%TKE_itidal_loss(:,:,:,:,:) = 0.0
  call get_param(param_file, mod, "KAPPA_ITIDES", kappa_itides, &
               "A topographic wavenumber used with INT_TIDE_DISSIPATION. \n"//&
               "The default is 2pi/10 km, as in St.Laurent et al. 2002.", &
               units="m-1", default=8.e-4*atan(1.0))
  call get_param(param_file, mod, "KAPPA_H2_FACTOR", kappa_h2_factor, &
               "A scaling factor for the roughness amplitude with n"//&
               "INT_TIDE_DISSIPATION.",  units="nondim", default=1.0)
  call get_param(param_file, mod, "H2_FILE", h2_file, &
               "The path to the file containing the sub-grid-scale \n"//&
               "topographic roughness amplitude with INT_TIDE_DISSIPATION.", &
               fail_if_missing=.true.)
  filename = trim(CS%inputdir) // trim(h2_file)
  call log_param(param_file, mod, "INPUTDIR/H2_FILE", filename)
  call read_data(filename, 'h2', h2, domain=G%domain%mpp_domain, timelevel=1)               
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! Restrict rms topo to 10 percent of column depth.
    h2(i,j) = min(0.01*G%bathyT(i,j)**2, h2(i,j))
    ! Compute the fixed part; units are [kg m-2] here; 
    ! will be multiplied by N and En to get into [W m-2]
    CS%TKE_itidal_loss_fixed(i,j) = 0.5*kappa_h2_factor*G%Rho0*&
         kappa_itides * h2(i,j)
  enddo; enddo
  
  ! Read in prescribed coast/ridge/shelf angles from file (BDM)
  call get_param(param_file, mod, "REFL_ANGLE_FILE", refl_angle_file, &
               "The path to the file containing the local angle of \n"//&
               "the coastline/ridge/shelf with respect to the equator.", &
               fail_if_missing=.false.)
  filename = trim(CS%inputdir) // trim(refl_angle_file)
  call log_param(param_file, mod, "INPUTDIR/REFL_ANGLE_FILE", filename)
  allocate(CS%refl_angle(isd:ied,jsd:jed)) ; CS%refl_angle(:,:) = 0.0
  call read_data(filename, 'refl_angle', CS%refl_angle, &
                 domain=G%domain%mpp_domain, timelevel=1)
  call pass_var(CS%refl_angle,G%domain)
  
  ! Read in prescribed partial reflection coefficients from file (BDM)
  call get_param(param_file, mod, "REFL_PREF_FILE", refl_pref_file, &
               "The path to the file containing the reflection coefficients.", &
               fail_if_missing=.false.)
  filename = trim(CS%inputdir) // trim(refl_pref_file)
  call log_param(param_file, mod, "INPUTDIR/REFL_PREF_FILE", filename)
  allocate(CS%refl_pref(isd:ied,jsd:jed)) ; CS%refl_pref(:,:) = 1.0
  call read_data(filename, 'refl_pref', CS%refl_pref, &
                 domain=G%domain%mpp_domain, timelevel=1)
  !CS%refl_pref = CS%refl_pref*1 ! adjust partial reflection
  call pass_var(CS%refl_pref,G%domain)
  
  ! Tag reflection cells with partial reflection (done here for speed)
  allocate(CS%refl_pref_logical(isd:ied,jsd:jed)) ; CS%refl_pref_logical(:,:) = .false.
  do j=jsd,jed
    do i=isd,ied
      ! flag cells with partial reflection
      if (.not. isnan(CS%refl_angle(i,j)) .and. &
        CS%refl_pref(i,j) < 1.0 .and. CS%refl_pref(i,j) > 0.0) then
        CS%refl_pref_logical(i,j) = .true.
      endif
    enddo
  enddo
  
  ! Read in double-reflective (ridge) tags from file (BDM)
  call get_param(param_file, mod, "REFL_DBL_FILE", refl_dbl_file, &
               "The path to the file containing the double-reflective ridge tags.", &
               fail_if_missing=.false.)
  filename = trim(CS%inputdir) // trim(refl_dbl_file)
  call log_param(param_file, mod, "INPUTDIR/REFL_DBL_FILE", filename)
  allocate(ridge_temp(isd:ied,jsd:jed)) ; ridge_temp(:,:) = 0.0
  call read_data(filename, 'refl_dbl', ridge_temp, &
                 domain=G%domain%mpp_domain, timelevel=1)
  call pass_var(ridge_temp,G%domain)               
  allocate(CS%refl_dbl(isd:ied,jsd:jed)) ; CS%refl_dbl(:,:) = .false.
  do i=isd,ied; do j=jsd,jed 
    if (ridge_temp(i,j) == 1) then; CS%refl_dbl(i,j) = .true.
    else ; CS%refl_dbl(i,j) = .false. ; endif
  enddo; enddo                 
  
  ! Read in prescribed land mask from file. 
  ! This should be done in MOM_initialize_topography subroutine 
  ! defined in MOM_fixed_initialization.F90 (BDM)   
  !call get_param(param_file, mod, "LAND_MASK_FILE", land_mask_file, &
  !             "The path to the file containing the land mask.", &
  !             fail_if_missing=.false.) 
  !filename = trim(CS%inputdir) // trim(land_mask_file)
  !call log_param(param_file, mod, "INPUTDIR/LAND_MASK_FILE", filename)
  !G%mask2dCu(:,:) = 1 ; G%mask2dCv(:,:) = 1 ; G%mask2dT(:,:)  = 1
  !call read_data(filename, 'land_mask', G%mask2dCu, &
  !               domain=G%domain%mpp_domain, timelevel=1)
  !call read_data(filename, 'land_mask', G%mask2dCv, &
  !               domain=G%domain%mpp_domain, timelevel=1)
  !call read_data(filename, 'land_mask', G%mask2dT, &
  !               domain=G%domain%mpp_domain, timelevel=1)
  !call pass_var(G%mask2dCu,G%domain)
  !call pass_var(G%mask2dCv,G%domain)
  !call pass_var(G%mask2dT,G%domain)

  ! Read in prescribed partial east face blockages from file (BDM)
  !call get_param(param_file, mod, "dy_Cu_FILE", dy_Cu_file, &
  !             "The path to the file containing the east face blockages.", &
  !             fail_if_missing=.false.) 
  !filename = trim(CS%inputdir) // trim(dy_Cu_file)
  !call log_param(param_file, mod, "INPUTDIR/dy_Cu_FILE", filename)
  !G%dy_Cu(:,:) = 0.0
  !call read_data(filename, 'dy_Cu', G%dy_Cu, &
  !               domain=G%domain%mpp_domain, timelevel=1)  
  !call pass_var(G%dy_Cu,G%domain)
                 
  ! Read in prescribed partial north face blockages from file (BDM)  
  !call get_param(param_file, mod, "dx_Cv_FILE", dx_Cv_file, &
  !             "The path to the file containing the north face blockages.", &
  !             fail_if_missing=.false.) 
  !filename = trim(CS%inputdir) // trim(dx_Cv_file)
  !call log_param(param_file, mod, "INPUTDIR/dx_Cv_FILE", filename)
  !G%dx_Cv(:,:) = 0.0
  !call read_data(filename, 'dx_Cv', G%dx_Cv, &
  !               domain=G%domain%mpp_domain, timelevel=1)  
  !call pass_var(G%dx_Cv,G%domain)                                      
  
  CS%id_refl_ang = register_diag_field('ocean_model', 'refl_angle', diag%axesT1, &
                 Time, 'Local angle of coastline/ridge/shelf with respect to equator', 'rad') ! (BDM)
  CS%id_refl_pref = register_diag_field('ocean_model', 'refl_pref', diag%axesT1, &
                 Time, 'Partial reflection coefficients', '') ! (BDM)
  CS%id_dx_Cv = register_diag_field('ocean_model', 'dx_Cv', diag%axesT1, &
                 Time, 'North face unblocked width', 'm') ! (BDM)
  CS%id_dy_Cu = register_diag_field('ocean_model', 'dy_Cu', diag%axesT1, &
                 Time, 'East face unblocked width', 'm') ! (BDM)
  CS%id_land_mask = register_diag_field('ocean_model', 'land_mask', diag%axesT1, &
                 Time, 'Land mask', 'logical') !(BDM)
  CS%id_tot_En = register_diag_field('ocean_model', 'ITide_tot_En', diag%axesT1, &
                 Time, 'Internal tide total energy density', 'J m-2') ! (BDM)
  CS%id_itide_drag = register_diag_field('ocean_model', 'ITide_drag', diag%axesT1, &
                 Time, 'Interior and bottom drag internal tide decay timescale', 's-1')
  CS%id_TKE_itidal_input = register_diag_field('ocean_model', 'TKE_itidal_input', diag%axesT1, &
                 Time, 'Conversion from barotropic to baroclinic tide, \n'//&
                 'a fraction of which goes into rays', 'W m-2') ! (BDM)

  allocate(CS%id_En_mode(CS%nFreq,CS%nMode)) ; CS%id_En_mode(:,:) = -1
  allocate(CS%id_En_ang_mode(CS%nFreq,CS%nMode)) ; CS%id_En_ang_mode(:,:) = -1
  allocate(CS%id_TKE_loss_mode(CS%nFreq,CS%nMode)) ; CS%id_TKE_loss_mode(:,:) = -1
  allocate(CS%id_TKE_loss_ang_mode(CS%nFreq,CS%nMode)) ; CS%id_TKE_loss_ang_mode(:,:) = -1

  allocate(angles(CS%NAngle)) ; angles(:) = 0.0
  Angle_size = (8.0*atan(1.0)) / (real(num_angle))
  do a=1,num_angle ; angles(a) = (real(a) - 1) * Angle_size ; enddo

  id_ang = diag_axis_init("angle", angles, "Radians", "N", "Angular Orienation of Fluxes")
  call defineAxes(diag, (/ diag%axesT1%handles(1), diag%axesT1%handles(2), id_ang /), axes_ang)

  do fr=1,CS%nFreq ; write(freq_name(fr), '("freq",i1)') fr ; enddo
  do m=1,CS%nMode ; do fr=1,CS%nFreq  
    ! Register 2-D energy density (summed over angles) for each freq and mode
    write(var_name, '("Itide_En_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy density in frequency ",i1," mode ",i1)') fr, m
    CS%id_En_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'J m-2')
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    
    ! Register 3-D (i,j,a) energy density for each freq and mode
    write(var_name, '("Itide_En_ang_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide angular energy density in frequency ",i1," mode ",i1)') fr, m
    CS%id_En_ang_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 axes_ang, Time, var_descript, 'J m-2 band-1')
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    
    ! Register 2-D energy loss (summed over angles) for each freq and mode
    write(var_name, '("Itide_TKE_loss_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy loss from frequency ",i1," mode ",i1)') fr, m
    CS%id_TKE_loss_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'W m-2')
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    
    ! Register 3-D (i,j,a) energy loss for each freq and mode
    write(var_name, '("Itide_TKE_loss_ang_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy loss from frequency ",i1," mode ",i1)') fr, m
    CS%id_TKE_loss_ang_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 axes_ang, Time, var_descript, 'W m-2 band-1')
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    
  enddo ; enddo
  
  !--------------------check----------------------------------------------
  if (is_root_pe()) then
    print *,'internal_tides_init: done!' !BDM
  endif
  !-----------------------------------------------------------------------
  
end subroutine internal_tides_init


subroutine internal_tides_end(CS)
  type(int_tide_CS),            pointer :: CS
  ! Arguments:  CS - A pointer to the control structure returned by a previous
  !                  call to internal_tides_init, it will be deallocated here.
  
  print *, 'sum En(:,:,:,1,1) = ', sum(CS%En(:,:,:,1,1)) ! BDM
  
  CS%En_restart(:,:,:) = CS%En(:,:,:,1,1)  ! added here as work-around (BDM) 
  !call save_restart(CS%restart_dir, CS%Time, G, CS%restart_CSp) !(BDM) 
   
  if (associated(CS)) then
    if (associated(CS%En)) deallocate(CS%En)
    if (allocated(CS%frequency)) deallocate(CS%frequency)
    if (allocated(CS%id_En_mode)) deallocate(CS%id_En_mode)
    deallocate(CS)
  endif
  CS => NULL()
  
  !--------------------check----------------------------------------------
  if (is_root_pe()) then
    print *,'internal_tides_end: CS%En added to CS%En_restart' !BDM
    print *,'internal_tides_end: done!' !BDM
  endif
  !-----------------------------------------------------------------------
  
end subroutine internal_tides_end


end module MOM_internal_tides
