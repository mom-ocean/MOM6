!> Subroutines that use the ray-tracing equations to propagate the internal tide energy density.
!!
!! \author Benjamin Mater & Robert Hallberg, 2015
module MOM_internal_tides

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : is_NaN
use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_axis_init
use MOM_diag_mediator, only : disable_averaging, enable_averages
use MOM_diag_mediator, only : register_diag_field, diag_ctrl, safe_alloc_ptr
use MOM_diag_mediator, only : axes_grp, define_axes_group
use MOM_domains, only       : AGRID, To_South, To_West, To_All
use MOM_domains, only       : create_group_pass, do_group_pass, pass_var
use MOM_domains, only       : group_pass_type, start_group_pass, complete_group_pass
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser, only   : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only          : ocean_grid_type
use MOM_io, only            : slasher, vardesc, MOM_read_data, file_exists
use MOM_restart, only       : register_restart_field, MOM_restart_CS, restart_init, save_restart
use MOM_spatial_means, only : global_area_mean
use MOM_time_manager, only  : time_type, time_type_to_real, operator(+), operator(/), operator(-)
use MOM_unit_scaling, only  : unit_scale_type
use MOM_variables, only     : surface, thermo_var_ptrs
use MOM_verticalGrid, only  : verticalGrid_type
use MOM_wave_structure, only: wave_structure_init, wave_structure, wave_structure_CS

!use, intrinsic :: IEEE_ARITHMETIC

implicit none ; private

#include <MOM_memory.h>

public propagate_int_tide !, register_int_tide_restarts
public internal_tides_init, internal_tides_end
public get_lowmode_loss

!> This control structure has parameters for the MOM_internal_tides module
type, public :: int_tide_CS ; private
  logical :: do_int_tides    !< If true, use the internal tide code.
  integer :: nFreq = 0       !< The number of internal tide frequency bands
  integer :: nMode = 1       !< The number of internal tide vertical modes
  integer :: nAngle = 24     !< The number of internal tide angular orientations
  integer :: energized_angle = -1 !< If positive, only this angular band is energized for debugging purposes
  logical :: corner_adv      !< If true, use a corner advection rather than PPM.
  logical :: upwind_1st      !< If true, use a first-order upwind scheme.
  logical :: simple_2nd      !< If true, use a simple second order (arithmetic mean) interpolation
                             !! of the edge values instead of the higher order interpolation.
  logical :: vol_CFL         !< If true, use the ratio of the open face lengths to the tracer cell
                             !! areas when estimating CFL numbers.  Without aggress_adjust,
                             !! the default is false; it is always true with aggress_adjust.
  logical :: use_PPMang      !< If true, use PPM for advection of energy in angular space.

  real, allocatable, dimension(:,:) :: refl_angle
                        !< local coastline/ridge/shelf angles read from file
                        ! (could be in G control structure)
  real :: nullangle = -999.9 !< placeholder value in cells with no reflection
  real, allocatable, dimension(:,:) :: refl_pref
                        !< partial reflection coeff for each "coast cell"
                        ! (could be in G control structure)
  logical, allocatable, dimension(:,:) :: refl_pref_logical
                        !< true if reflecting cell with partial reflection
                        ! (could be in G control structure)
  logical, allocatable, dimension(:,:) :: refl_dbl
                        !< identifies reflection cells where double reflection
                        !! is possible (i.e. ridge cells)
                        ! (could be in G control structure)
  real, allocatable, dimension(:,:,:,:) :: cp
                        !< horizontal phase speed [L T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_leak_loss
                        !< energy lost due to misc background processes [R Z3 T-3 ~> W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_quad_loss
                        !< energy lost due to quadratic bottom drag [R Z3 T-3 ~> W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_Froude_loss
                        !< energy lost due to wave breaking [R Z3 T-3 ~> W m-2]
  real, allocatable, dimension(:,:) :: TKE_itidal_loss_fixed
                        !< Fixed part of the energy lost due to small-scale drag [R L-2 Z3 ~> kg m-2] here;
                        !! This will be multiplied by N and the squared near-bottom velocity to get
                        !! the energy losses in [R Z3 T-3 ~> W m-2]
  real, allocatable, dimension(:,:,:,:,:) :: TKE_itidal_loss
                        !< energy lost due to small-scale wave drag [R Z3 T-3 ~> W m-2]
  real, allocatable, dimension(:,:) :: tot_leak_loss !< Energy loss rates due to misc bakground processes,
                        !! summed over angle, frequency and mode [R Z3 T-3 ~> W m-2]
  real, allocatable, dimension(:,:) :: tot_quad_loss !< Energy loss rates due to quadratic bottom drag,
                        !! summed over angle, frequency and mode [R Z3 T-3 ~> W m-2]
  real, allocatable, dimension(:,:) :: tot_itidal_loss !< Energy loss rates due to small-scale drag,
                        !! summed over angle, frequency and mode [R Z3 T-3 ~> W m-2]
  real, allocatable, dimension(:,:) :: tot_Froude_loss !< Energy loss rates due to wave breaking,
                        !! summed over angle, frequency and mode [R Z3 T-3 ~> W m-2]
  real, allocatable, dimension(:,:) :: tot_allprocesses_loss !< Energy loss rates due to all processes,
                        !! summed over angle, frequency and mode [R Z3 T-3 ~> W m-2]
  real :: q_itides      !< fraction of local dissipation [nondim]
  real :: En_sum        !< global sum of energy for use in debugging [R Z3 T-2 ~> J m-2]
  type(time_type), pointer :: Time => NULL() !< A pointer to the model's clock.
  character(len=200) :: inputdir !< directory to look for coastline angle file
  real :: decay_rate    !< A constant rate at which internal tide energy is
                        !! lost to the interior ocean internal wave field.
  real :: cdrag         !< The bottom drag coefficient [nondim].
  logical :: apply_background_drag
                        !< If true, apply a drag due to background processes as a sink.
  logical :: apply_bottom_drag
                        !< If true, apply a quadratic bottom drag as a sink.
  logical :: apply_wave_drag
                        !< If true, apply scattering due to small-scale roughness as a sink.
  logical :: apply_Froude_drag
                        !< If true, apply wave breaking as a sink.
  real, dimension(:,:,:,:,:), pointer :: En => NULL()
                        !< The internal wave energy density as a function of (i,j,angle,frequency,mode)
                        !! integrated within an angular and frequency band [R Z3 T-2 ~> J m-2]
  real, dimension(:,:,:), pointer :: En_restart => NULL()
                        !< The internal wave energy density as a function of (i,j,angle); temporary for restart
  real, allocatable, dimension(:) :: frequency  !< The frequency of each band [T-1 ~> s-1].

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                        !! timing of diagnostic output.
  type(wave_structure_CS), pointer :: wave_structure_CSp => NULL()
                        !< A pointer to the wave_structure module control structure

  !>@{ Diag handles
  ! Diag handles relevant to all modes, frequencies, and angles
  integer :: id_tot_En = -1, id_TKE_itidal_input = -1, id_itide_drag = -1
  integer :: id_refl_pref = -1, id_refl_ang = -1, id_land_mask = -1
  integer :: id_dx_Cv = -1, id_dy_Cu = -1
  ! Diag handles considering: sums over all modes, frequencies, and angles
  integer :: id_tot_leak_loss = -1, id_tot_quad_loss = -1, id_tot_itidal_loss = -1
  integer :: id_tot_Froude_loss = -1, id_tot_allprocesses_loss = -1
  ! Diag handles considering: all modes & freqs; summed over angles
  integer, allocatable, dimension(:,:) :: &
             id_En_mode, &
             id_itidal_loss_mode, &
             id_allprocesses_loss_mode, &
             id_Ub_mode, &
             id_cp_mode
  ! Diag handles considering: all modes, freqs, and angles
  integer, allocatable, dimension(:,:) :: &
             id_En_ang_mode, &
             id_itidal_loss_ang_mode
  !>@}

end type int_tide_CS

!> A structure with the active energy loop bounds.
type :: loop_bounds_type ; private
  !>@{ The active loop bounds
  integer :: ish, ieh, jsh, jeh
  !>@}
end type loop_bounds_type

contains

!> Calls subroutines in this file that are needed to refract, propagate,
!! and dissipate energy density of the internal tide.
subroutine propagate_int_tide(h, tv, cn, TKE_itidal_input, vel_btTide, Nb, dt, &
                              G, GV, US, CS)
  type(ocean_grid_type),            intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),          intent(in)    :: GV !< The ocean's vertical grid structure.
  type(unit_scale_type),            intent(in)    :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                    intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),            intent(in)    :: tv !< Pointer to thermodynamic variables
                                                        !! (needed for wave structure).
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: TKE_itidal_input !< The energy input to the
                                                        !! internal waves [R Z3 T-3 ~> W m-2].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: vel_btTide !< Barotropic velocity read
                                                        !! from file [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: Nb !< Near-bottom buoyancy frequency [T-1 ~> s-1].
  real,                             intent(in)    :: dt !< Length of time over which to advance
                                                        !! the internal tides [T ~> s].
  type(int_tide_CS),                pointer       :: CS !< The control structure returned by a
                                                        !! previous call to int_tide_init.
  real, dimension(SZI_(G),SZJ_(G),CS%nMode), &
                                    intent(in)    :: cn !< The internal wave speeds of each
                                                        !! mode [L T-1 ~> m s-1].
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),2) :: &
    test
  real, dimension(SZI_(G),SZJ_(G),CS%nFreq,CS%nMode) :: &
    tot_En_mode, & ! energy summed over angles only [R Z3 T-2 ~> J m-2]
    Ub, &          ! near-bottom horizontal velocity of wave (modal) [L T-1 ~> m s-1]
    Umax           ! Maximum horizontal velocity of wave (modal) [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    flux_heat_y, &
    flux_prec_y
  real, dimension(SZI_(G),SZJ_(G)) :: &
    tot_En, &      ! energy summed over angles, modes, frequencies [R Z3 T-2 ~> J m-2]
    tot_leak_loss, tot_quad_loss, tot_itidal_loss, tot_Froude_loss, tot_allprocesses_loss, &
                   ! energy loss rates summed over angle, freq, and mode [R Z3 T-3 ~> W m-2]
    drag_scale, &  ! bottom drag scale [T-1 ~> s-1]
    itidal_loss_mode, allprocesses_loss_mode
                   ! energy loss rates for a given mode and frequency (summed over angles) [R Z3 T-3 ~> W m-2]
  real :: frac_per_sector, f2, Kmag2
  real :: I_D_here ! The inverse of the local depth [Z-1 ~> m-1]
  real :: I_rho0  ! The inverse fo the Boussinesq density [R-1 ~> m3 kg-1]
  real :: freq2 ! The frequency squared [T-2 ~> s-2]
  real :: c_phase ! The phase speed [L T-1 ~> m s-1]
  real :: loss_rate  ! An energy loss rate [T-1 ~> s-1]
  real :: Fr2_max
  real :: cn_subRO        ! A tiny wave speed to prevent division by zero [L T-1 ~> m s-1]
  real :: En_new, En_check                           ! Energies for debugging [R Z3 T-2 ~> J m-2]
  real :: En_initial, Delta_E_check                  ! Energies for debugging [R Z3 T-2 ~> J m-2]
  real :: TKE_Froude_loss_check, TKE_Froude_loss_tot ! Energy losses for debugging [R Z3 T-3 ~> W m-2]
  character(len=160) :: mesg  ! The text of an error message
  integer :: a, m, fr, i, j, is, ie, js, je, isd, ied, jsd, jed, nAngle, nzm
  integer :: id_g, jd_g         ! global (decomp-invar) indices (for debugging)
  type(group_pass_type), save :: pass_test, pass_En
  type(time_type) :: time_end
  logical:: avg_enabled

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nAngle = CS%NAngle
  I_rho0 = 1.0 / GV%Rho0
  cn_subRO = 1e-100*US%m_s_to_L_T  ! The hard-coded value here might need to increase.

  ! Set the wave speeds for the modes, using cg(n) ~ cg(1)/n.**********************
  ! This is wrong, of course, but it works reasonably in some cases.
  ! Uncomment if wave_speed is not used to calculate the true values (BDM).
  !do m=1,CS%nMode ; do j=jsd,jed ; do i=isd,ied
  !  cn(i,j,m) = cn(i,j,1) / real(m)
  !enddo ; enddo ; enddo

  ! Add the forcing.***************************************************************
  if (CS%energized_angle <= 0) then
    frac_per_sector = 1.0 / real(CS%nAngle * CS%nMode * CS%nFreq)
    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
      f2 = 0.25*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                 (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2))
      if (CS%frequency(fr)**2 > f2) &
        CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) + dt*frac_per_sector*(1.0-CS%q_itides) * &
                            TKE_itidal_input(i,j)
    enddo ; enddo ; enddo ; enddo ; enddo
  elseif (CS%energized_angle <= CS%nAngle) then
    frac_per_sector = 1.0 / real(CS%nMode * CS%nFreq)
    a = CS%energized_angle
    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do j=js,je ; do i=is,ie
      f2 = 0.25*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                 (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2))
      if (CS%frequency(fr)**2 > f2) &
        CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) + dt*frac_per_sector*(1.0-CS%q_itides) * &
                            TKE_itidal_input(i,j)
    enddo ; enddo ; enddo ; enddo
  else
    call MOM_error(WARNING, "Internal tide energy is being put into a angular "//&
                            "band that does not exist.")
  endif

  ! Pass a test vector to check for grid rotation in the halo updates.
  do j=jsd,jed ; do i=isd,ied ; test(i,j,1) = 1.0 ; test(i,j,2) = 0.0 ; enddo ; enddo
  do m=1,CS%nMode ; do fr=1,CS%nFreq
    call create_group_pass(pass_En, CS%En(:,:,:,fr,m), G%domain)
  enddo ; enddo
  call create_group_pass(pass_test, test(:,:,1), test(:,:,2), G%domain, stagger=AGRID)
  call start_group_pass(pass_test, G%domain)

  ! Apply half the refraction.
  do m=1,CS%nMode ; do fr=1,CS%nFreq
    call refract(CS%En(:,:,:,fr,m), cn(:,:,m), CS%frequency(fr), 0.5*dt, &
                 G, US, CS%nAngle, CS%use_PPMang)
  enddo ; enddo

  ! Check for En<0 - for debugging, delete later
  do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
    do j=js,je ; do i=is,ie
      if (CS%En(i,j,a,fr,m)<0.0) then
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
        write(mesg,*) 'After first refraction: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                      'En=',CS%En(i,j,a,fr,m)
        call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg))
        CS%En(i,j,a,fr,m) = 0.0
!        call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
      endif
    enddo ; enddo
  enddo ; enddo ; enddo

  call do_group_pass(pass_En, G%domain)

  call complete_group_pass(pass_test, G%domain)

  ! Rotate points in the halos as necessary.
  call correct_halo_rotation(CS%En, test, G, CS%nAngle)

  ! Propagate the waves.
  do m=1,CS%NMode ; do fr=1,CS%Nfreq
    call propagate(CS%En(:,:,:,fr,m), cn(:,:,m), CS%frequency(fr), dt, &
                   G, US, CS, CS%NAngle)
  enddo ; enddo

  ! Check for En<0 - for debugging, delete later
  do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
    do j=js,je ; do i=is,ie
      if (CS%En(i,j,a,fr,m)<0.0) then
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset
        if (abs(CS%En(i,j,a,fr,m))>1.0) then ! only print if large
          write(mesg,*)  'After propagation: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                         'En=', CS%En(i,j,a,fr,m)
          call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg))
!          call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
        endif
        CS%En(i,j,a,fr,m) = 0.0
      endif
    enddo ; enddo
  enddo ; enddo ; enddo

  ! Apply the other half of the refraction.
  do m=1,CS%NMode ; do fr=1,CS%Nfreq
    call refract(CS%En(:,:,:,fr,m), cn(:,:,m), CS%frequency(fr), 0.5*dt, &
                 G, US, CS%NAngle, CS%use_PPMang)
  enddo ; enddo

  ! Check for En<0 - for debugging, delete later
  do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
    do j=js,je ; do i=is,ie
      if (CS%En(i,j,a,fr,m)<0.0) then
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
        write(mesg,*) 'After second refraction: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                      'En=',CS%En(i,j,a,fr,m)
        call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg))
        CS%En(i,j,a,fr,m) = 0.0
!        call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
      endif
    enddo ; enddo
  enddo ; enddo ; enddo

  ! Apply various dissipation mechanisms.
  if (CS%apply_background_drag .or. CS%apply_bottom_drag &
      .or. CS%apply_wave_drag .or. CS%apply_Froude_drag &
      .or. (CS%id_tot_En > 0)) then
    tot_En(:,:) = 0.0
    tot_En_mode(:,:,:,:) = 0.0
    do m=1,CS%NMode ; do fr=1,CS%Nfreq
      do j=jsd,jed ; do i=isd,ied ; do a=1,CS%nAngle
        tot_En(i,j) = tot_En(i,j) + CS%En(i,j,a,fr,m)
        tot_En_mode(i,j,fr,m) = tot_En_mode(i,j,fr,m) + CS%En(i,j,a,fr,m)
      enddo ; enddo ; enddo
    enddo ; enddo
  endif

  ! Extract the energy for mixing due to misc. processes (background leakage)------
  if (CS%apply_background_drag) then
    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
    ! Calculate loss rate and apply loss over the time step ; apply the same drag timescale
    ! to each En component (technically not correct; fix later)
    CS%TKE_leak_loss(i,j,a,fr,m)  = CS%En(i,j,a,fr,m) * CS%decay_rate ! loss rate [R Z3 T-3 ~> W m-2]
    CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) / (1.0 + dt * CS%decay_rate) ! implicit update
    enddo ; enddo ; enddo ; enddo ; enddo
  endif
  ! Check for En<0 - for debugging, delete later
  do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
    do j=js,je ; do i=is,ie
      if (CS%En(i,j,a,fr,m)<0.0) then
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
        write(mesg,*) 'After leak loss: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                      'En=',CS%En(i,j,a,fr,m)
        call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg), all_print=.true.)
        CS%En(i,j,a,fr,m) = 0.0
!       call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
      endif
    enddo ; enddo
  enddo ; enddo ; enddo

  ! Extract the energy for mixing due to bottom drag-------------------------------
  if (CS%apply_bottom_drag) then
    do j=jsd,jed ; do i=isd,ied
      ! Note the 1 m dimensional scale here.  Should this be a parameter?
      I_D_here = 1.0 / (max(G%bathyT(i,j), 1.0*US%m_to_Z))
      drag_scale(i,j) = CS%cdrag * sqrt(max(0.0, US%L_to_Z**2*vel_btTide(i,j)**2 + &
                        tot_En(i,j) * I_rho0 * I_D_here)) * I_D_here
    enddo ; enddo
    do m=1,CS%nMode ; do fr=1,CS%nFreq ; do a=1,CS%nAngle ; do j=jsd,jed ; do i=isd,ied
      ! Calculate loss rate and apply loss over the time step ; apply the same drag timescale
      ! to each En component (technically not correct; fix later)
      CS%TKE_quad_loss(i,j,a,fr,m)  = CS%En(i,j,a,fr,m) * drag_scale(i,j) ! loss rate
      CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m) / (1.0 + dt * drag_scale(i,j)) ! implicit update
    enddo ; enddo ; enddo ; enddo ; enddo
  endif
  ! Check for En<0 - for debugging, delete later
  do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
    do j=js,je ; do i=is,ie
      if (CS%En(i,j,a,fr,m)<0.0) then
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
        write(mesg,*) 'After bottom loss: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                      'En=',CS%En(i,j,a,fr,m)
        call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg), all_print=.true.)
        CS%En(i,j,a,fr,m) = 0.0
!       call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
        !stop
      endif
    enddo ; enddo
  enddo ; enddo ; enddo

  ! Extract the energy for mixing due to scattering (wave-drag)--------------------
  ! still need to allow a portion of the extracted energy to go to higher modes.
  ! First, find velocity profiles
  if (CS%apply_wave_drag .or. CS%apply_Froude_drag) then
    do m=1,CS%NMode ; do fr=1,CS%Nfreq
      ! Calculate modal structure for given mode and frequency
      call wave_structure(h, tv, G, GV, US, cn(:,:,m), m, CS%frequency(fr), &
                          CS%wave_structure_CSp, tot_En_mode(:,:,fr,m), full_halos=.true.)
      ! Pick out near-bottom and max horizontal baroclinic velocity values at each point
      do j=jsd,jed ; do i=isd,ied
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
        nzm = CS%wave_structure_CSp%num_intfaces(i,j)
        Ub(i,j,fr,m) = CS%wave_structure_CSp%Uavg_profile(i,j,nzm)
        Umax(i,j,fr,m) = maxval(CS%wave_structure_CSp%Uavg_profile(i,j,1:nzm))
      enddo ; enddo ! i-loop, j-loop
    enddo ; enddo ! fr-loop, m-loop
  endif ! apply_wave or _Froude_drag (Ub or Umax needed)
  ! Finally, apply loss
  if (CS%apply_wave_drag) then
    ! Calculate loss rate and apply loss over the time step
    call itidal_lowmode_loss(G, US, CS, Nb, Ub, CS%En, CS%TKE_itidal_loss_fixed, &
                             CS%TKE_itidal_loss, dt, full_halos=.false.)
  endif
  ! Check for En<0 - for debugging, delete later
  do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
    do j=js,je ; do i=is,ie
      if (CS%En(i,j,a,fr,m)<0.0) then
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
        write(mesg,*) 'After wave drag loss: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                      'En=',CS%En(i,j,a,fr,m)
        call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg), all_print=.true.)
        CS%En(i,j,a,fr,m) = 0.0
!       call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
      endif
    enddo ; enddo
  enddo ; enddo ; enddo

  ! Extract the energy for mixing due to wave breaking-----------------------------
  if (CS%apply_Froude_drag) then
    ! Pick out maximum baroclinic velocity values; calculate Fr=max(u)/cg
    do m=1,CS%NMode ; do fr=1,CS%Nfreq
      freq2 = CS%frequency(fr)**2
      do j=jsd,jed ; do i=isd,ied
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset ! for debugging
        ! Calculate horizontal phase velocity magnitudes
        f2 = 0.25*((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                   (G%CoriolisBu(I-1,J)**2 + G%CoriolisBu(I,J-1)**2))
        Kmag2 = (freq2 - f2) / (cn(i,j,m)**2 + cn_subRO**2)
        c_phase = 0.0
        if (Kmag2 > 0.0) then
          c_phase = sqrt(freq2/Kmag2)
          nzm = CS%wave_structure_CSp%num_intfaces(i,j)
          Fr2_max = (Umax(i,j,fr,m) / c_phase)**2
          ! Dissipate energy if Fr>1; done here with an arbitrary time scale
          if (Fr2_max > 1.0) then
            En_initial = sum(CS%En(i,j,:,fr,m)) ! for debugging
            ! Calculate effective decay rate [T-1 ~> s-1] if breaking occurs over a time step
            loss_rate = (1.0 - Fr2_max) / (Fr2_max * dt)
            do a=1,CS%nAngle
              ! Determine effective dissipation rate (Wm-2)
              CS%TKE_Froude_loss(i,j,a,fr,m) = CS%En(i,j,a,fr,m) * abs(loss_rate)
              ! Update energy
              En_new = CS%En(i,j,a,fr,m)/Fr2_max ! for debugging
              En_check = CS%En(i,j,a,fr,m) - CS%TKE_Froude_loss(i,j,a,fr,m)*dt ! for debugging
              ! Re-scale (reduce) energy due to breaking
              CS%En(i,j,a,fr,m) = CS%En(i,j,a,fr,m)/Fr2_max
              ! Check (for debugging only)
              if (abs(En_new - En_check) > 1e-10*US%kg_m3_to_R*US%m_to_Z**3*US%T_to_s**2) then
                call MOM_error(WARNING, "MOM_internal_tides: something is wrong with Fr-breaking.", &
                               all_print=.true.)
                write(mesg,*) "En_new=", En_new , "En_check=", En_check
                call MOM_error(WARNING, "MOM_internal_tides: "//trim(mesg), all_print=.true.)
              endif
            enddo
            ! Check (for debugging)
            Delta_E_check = En_initial - sum(CS%En(i,j,:,fr,m))
            TKE_Froude_loss_check = abs(Delta_E_check)/dt
            TKE_Froude_loss_tot = sum(CS%TKE_Froude_loss(i,j,:,fr,m))
            if (abs(TKE_Froude_loss_check - TKE_Froude_loss_tot) > 1e-10) then
              call MOM_error(WARNING, "MOM_internal_tides: something is wrong with Fr energy update.", &
                             all_print=.true.)
              write(mesg,*) "TKE_Froude_loss_check=", TKE_Froude_loss_check, &
                            "TKE_Froude_loss_tot=", TKE_Froude_loss_tot
              call MOM_error(WARNING, "MOM_internal_tides: "//trim(mesg), all_print=.true.)
            endif
          endif ! Fr2>1
        endif ! Kmag2>0
        CS%cp(i,j,fr,m) = c_phase
      enddo ; enddo
    enddo ; enddo
  endif
  ! Check for En<0 - for debugging, delete later
  do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle
    do j=js,je ; do i=is,ie
      if (CS%En(i,j,a,fr,m)<0.0) then
        id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset
        write(mesg,*) 'After Froude loss: En<0.0 at ig=', id_g, ', jg=', jd_g, &
                      'En=',CS%En(i,j,a,fr,m)
        call MOM_error(WARNING, "propagate_int_tide: "//trim(mesg), all_print=.true.)
        CS%En(i,j,a,fr,m) = 0.0
!       call MOM_error(FATAL, "propagate_int_tide: stopped due to negative energy.")
        !stop
      endif
    enddo ; enddo
  enddo ; enddo ; enddo

  ! Check for energy conservation on computational domain.*************************
  do m=1,CS%NMode ; do fr=1,CS%Nfreq
    call sum_En(G,CS,CS%En(:,:,:,fr,m),'prop_int_tide')
  enddo ; enddo

  ! Output diagnostics.************************************************************
  avg_enabled = query_averaging_enabled(CS%diag, time_end=time_end)
  call enable_averages(dt, time_end, CS%diag)

  if (query_averaging_enabled(CS%diag)) then
    ! Output two-dimensional diagnostistics
    if (CS%id_tot_En > 0)     call post_data(CS%id_tot_En, tot_En, CS%diag)
    if (CS%id_itide_drag > 0) call post_data(CS%id_itide_drag, drag_scale, CS%diag)
    if (CS%id_TKE_itidal_input > 0) call post_data(CS%id_TKE_itidal_input, &
                                                   TKE_itidal_input, CS%diag)

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

    ! Output 2-D energy loss (summed over angles, freq, modes)
    tot_leak_loss(:,:)   = 0.0
    tot_quad_loss(:,:)   = 0.0
    tot_itidal_loss(:,:) = 0.0
    tot_Froude_loss(:,:) = 0.0
    tot_allprocesses_loss(:,:) = 0.0
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
      tot_leak_loss(i,j)   = tot_leak_loss(i,j)   + CS%TKE_leak_loss(i,j,a,fr,m)
      tot_quad_loss(i,j)   = tot_quad_loss(i,j)   + CS%TKE_quad_loss(i,j,a,fr,m)
      tot_itidal_loss(i,j) = tot_itidal_loss(i,j) + CS%TKE_itidal_loss(i,j,a,fr,m)
      tot_Froude_loss(i,j) = tot_Froude_loss(i,j) + CS%TKE_Froude_loss(i,j,a,fr,m)
    enddo ; enddo ; enddo ; enddo ; enddo
    do j=js,je ; do i=is,ie
      tot_allprocesses_loss(i,j) = tot_leak_loss(i,j) + tot_quad_loss(i,j) + &
                          tot_itidal_loss(i,j) + tot_Froude_loss(i,j)
    enddo ; enddo
    CS%tot_leak_loss         = tot_leak_loss
    CS%tot_quad_loss         = tot_quad_loss
    CS%tot_itidal_loss       = tot_itidal_loss
    CS%tot_Froude_loss       = tot_Froude_loss
    CS%tot_allprocesses_loss = tot_allprocesses_loss
    if (CS%id_tot_leak_loss > 0) then
      call post_data(CS%id_tot_leak_loss, tot_leak_loss, CS%diag)
    endif
    if (CS%id_tot_quad_loss > 0) then
      call post_data(CS%id_tot_quad_loss, tot_quad_loss, CS%diag)
    endif
    if (CS%id_tot_itidal_loss > 0) then
      call post_data(CS%id_tot_itidal_loss, tot_itidal_loss, CS%diag)
    endif
    if (CS%id_tot_Froude_loss > 0) then
      call post_data(CS%id_tot_Froude_loss, tot_Froude_loss, CS%diag)
    endif
    if (CS%id_tot_allprocesses_loss > 0) then
      call post_data(CS%id_tot_allprocesses_loss, tot_allprocesses_loss, CS%diag)
    endif

    ! Output 2-D energy loss (summed over angles) for each freq and mode
    do m=1,CS%NMode ; do fr=1,CS%Nfreq
    if (CS%id_itidal_loss_mode(fr,m) > 0 .or. CS%id_allprocesses_loss_mode(fr,m) > 0) then
      itidal_loss_mode(:,:) = 0.0 ! wave-drag processes (could do others as well)
      allprocesses_loss_mode(:,:) = 0.0 ! all processes summed together
      do a=1,CS%nAngle ; do j=js,je ; do i=is,ie
        itidal_loss_mode(i,j)       = itidal_loss_mode(i,j) + CS%TKE_itidal_loss(i,j,a,fr,m)
        allprocesses_loss_mode(i,j) = allprocesses_loss_mode(i,j) + &
                                 CS%TKE_leak_loss(i,j,a,fr,m) + CS%TKE_quad_loss(i,j,a,fr,m) + &
                                 CS%TKE_itidal_loss(i,j,a,fr,m) + CS%TKE_Froude_loss(i,j,a,fr,m)
      enddo ; enddo ; enddo
      call post_data(CS%id_itidal_loss_mode(fr,m), itidal_loss_mode, CS%diag)
      call post_data(CS%id_allprocesses_loss_mode(fr,m), allprocesses_loss_mode, CS%diag)
    endif ; enddo ; enddo

    ! Output 3-D (i,j,a) energy loss for each freq and mode
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; if (CS%id_itidal_loss_ang_mode(fr,m) > 0) then
      call post_data(CS%id_itidal_loss_ang_mode(fr,m), CS%TKE_itidal_loss(:,:,:,fr,m) , CS%diag)
    endif ; enddo ; enddo

    ! Output 2-D period-averaged horizontal near-bottom mode velocity for each freq and mode
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; if (CS%id_Ub_mode(fr,m) > 0) then
      call post_data(CS%id_Ub_mode(fr,m), Ub(:,:,fr,m), CS%diag)
    endif ; enddo ; enddo

    ! Output 2-D horizontal phase velocity for each freq and mode
    do m=1,CS%NMode ; do fr=1,CS%Nfreq ; if (CS%id_cp_mode(fr,m) > 0) then
      call post_data(CS%id_cp_mode(fr,m), CS%cp(:,:,fr,m), CS%diag)
    endif ; enddo ; enddo

  endif

  call disable_averaging(CS%diag)

end subroutine propagate_int_tide

!> Checks for energy conservation on computational domain
subroutine sum_En(G, CS, En, label)
  type(ocean_grid_type),  intent(in) :: G  !< The ocean's grid structure.
  type(int_tide_CS),      pointer    :: CS !< The control structure returned by a
                                           !! previous call to int_tide_init.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle), &
                          intent(in) :: En !< The energy density of the internal tides [R Z3 T-2 ~> J m-2].
  character(len=*),       intent(in) :: label !< A label to use in error messages
  ! Local variables
  real :: En_sum   ! The total energy [R Z3 T-2 ~> J m-2]
  real :: tmpForSumming
  integer :: m,fr,a
  ! real :: En_sum_diff, En_sum_pdiff
  ! character(len=160) :: mesg  ! The text of an error message
  ! real :: days

  En_sum = 0.0
  tmpForSumming = 0.0
  do a=1,CS%nAngle
    tmpForSumming = global_area_mean(En(:,:,a),G)*G%areaT_global
    En_sum = En_sum + tmpForSumming
  enddo
  CS%En_sum = En_sum
  !En_sum_diff = En_sum - CS%En_sum
  !if (CS%En_sum /= 0.0) then
  !  En_sum_pdiff= (En_sum_diff/CS%En_sum)*100.0
  !else
  !  En_sum_pdiff= 0.0
  !endif
  !! Print to screen
  !if (is_root_pe()) then
  !  days = time_type_to_real(CS%Time) / 86400.0
  !  write(mesg,*) trim(label)//': days =', days, ', En_sum=', En_sum, &
  !                ', En_sum_diff=', En_sum_diff, ', Percent change=', En_sum_pdiff, '%'
  !  call MOM_mesg(mesg)
  !if (is_root_pe() .and. (abs(En_sum_pdiff) > 1.0)) &
  !  call MOM_error(FATAL, "Run stopped due to excessive internal tide energy change.")
  !endif

end subroutine sum_En

!> Calculates the energy lost from the propagating internal tide due to
!! scattering over small-scale roughness along the lines of Jayne & St. Laurent (2001).
subroutine itidal_lowmode_loss(G, US, CS, Nb, Ub, En, TKE_loss_fixed, TKE_loss, dt, full_halos)
  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure.
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  type(int_tide_CS),         pointer       :: CS !< The control structure returned by a
                                                 !! previous call to int_tide_init.
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                             intent(in)    :: Nb !< Near-bottom stratification [T-1 ~> s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%nFreq,CS%nMode), &
                             intent(inout) :: Ub !< RMS (over one period) near-bottom horizontal
                                                 !! mode velocity [L T-1 ~> m s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                             intent(in) :: TKE_loss_fixed !< Fixed part of energy loss [R L-2 Z3 ~> kg m-2]
                                                 !! (rho*kappa*h^2).
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle,CS%nFreq,CS%nMode), &
                             intent(inout) :: En !< Energy density of the internal waves [R Z3 T-2 ~> J m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,CS%NAngle,CS%nFreq,CS%nMode), &
                             intent(out)   :: TKE_loss    !< Energy loss rate [R Z3 T-3 ~> W m-2]
                                                 !! (q*rho*kappa*h^2*N*U^2).
  real,                      intent(in)    :: dt !< Time increment [T ~> s].
  logical,optional,          intent(in)    :: full_halos  !< If true, do the calculation over the
                                                 !! entirecomputational domain.
  ! Local variables
  integer :: j,i,m,fr,a, is, ie, js, je
  real    :: En_tot          ! energy for a given mode, frequency, and point summed over angles [R Z3 T-2 ~> J m-2]
  real    :: TKE_loss_tot    ! dissipation for a given mode, frequency, and point summed over angles [R Z3 T-3 ~> W m-2]
  real    :: TKE_sum_check   ! temporary for check summing
  real    :: frac_per_sector ! fraction of energy in each wedge
  real    :: q_itides        ! fraction of energy actually lost to mixing (remainder, 1-q, is
                             ! assumed to stay in propagating mode for now - BDM)
  real    :: loss_rate       ! approximate loss rate for implicit calc [T-1 ~> s-1]
  real    :: En_negl         ! negilibly small number to prevent division by zero

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  q_itides = CS%q_itides
  En_negl = 1e-30*US%kg_m3_to_R*US%m_to_Z**3*US%T_to_s**2

  if (present(full_halos)) then ; if (full_halos) then
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif ; endif

  do j=js,je ; do i=is,ie ; do m=1,CS%nMode ; do fr=1,CS%nFreq

    ! Sum energy across angles
    En_tot = 0.0
    do a=1,CS%nAngle
      En_tot = En_tot + En(i,j,a,fr,m)
    enddo

    ! Calculate TKE loss rate; units of [R Z3 T-3 ~> W m-2] here.
    TKE_loss_tot = q_itides * TKE_loss_fixed(i,j) * Nb(i,j) * Ub(i,j,fr,m)**2

    ! Update energy remaining (this is a pseudo implicit calc)
    ! (E(t+1)-E(t))/dt = -TKE_loss(E(t+1)/E(t)), which goes to zero as E(t+1) goes to zero
    if (En_tot > 0.0) then
      do a=1,CS%nAngle
        frac_per_sector = En(i,j,a,fr,m)/En_tot
        TKE_loss(i,j,a,fr,m) = frac_per_sector*TKE_loss_tot           ! Wm-2
        loss_rate = TKE_loss(i,j,a,fr,m) / (En(i,j,a,fr,m) + En_negl) ! [T-1 ~> s-1]
        En(i,j,a,fr,m) = En(i,j,a,fr,m) / (1.0 + dt*loss_rate)
      enddo
    else
      ! no loss if no energy
      TKE_loss(i,j,:,fr,m) = 0.0
    endif

    ! Update energy remaining (this is the old explicit calc)
    !if (En_tot > 0.0) then
    !  do a=1,CS%nAngle
    !    frac_per_sector = En(i,j,a,fr,m)/En_tot
    !    TKE_loss(i,j,a,fr,m) = frac_per_sector*TKE_loss_tot
    !    if (TKE_loss(i,j,a,fr,m)*dt <= En(i,j,a,fr,m))then
    !      En(i,j,a,fr,m) = En(i,j,a,fr,m) - TKE_loss(i,j,a,fr,m)*dt
    !    else
    !      call MOM_error(WARNING, "itidal_lowmode_loss: energy loss greater than avalable, "// &
    !                        " setting En to zero.", all_print=.true.)
    !      En(i,j,a,fr,m) = 0.0
    !    endif
    !  enddo
    !else
    !  ! no loss if no energy
    !  TKE_loss(i,j,:,fr,m) = 0.0
    !endif

  enddo ; enddo ; enddo ; enddo

end subroutine itidal_lowmode_loss

!> This subroutine extracts the energy lost from the propagating internal which has
!> been summed across all angles, frequencies, and modes for a given mechanism and location.
!!
!> It can be called from another module to get values from this module's (private) CS.
subroutine get_lowmode_loss(i,j,G,CS,mechanism,TKE_loss_sum)
  integer,               intent(in)  :: i   !< The i-index of the value to be reported.
  integer,               intent(in)  :: j   !< The j-index of the value to be reported.
  type(ocean_grid_type), intent(in)  :: G   !< The ocean's grid structure
  type(int_tide_CS),     pointer     :: CS  !< The control structure returned by a
                                            !! previous call to int_tide_init.
  character(len=*),      intent(in)  :: mechanism    !< The named mechanism of loss to return
  real,                  intent(out) :: TKE_loss_sum !< Total energy loss rate due to specified
                                                     !! mechanism [R Z3 T-3 ~> W m-2].

  if (mechanism == 'LeakDrag') TKE_loss_sum = CS%tot_leak_loss(i,j)   ! not used for mixing yet
  if (mechanism == 'QuadDrag') TKE_loss_sum = CS%tot_quad_loss(i,j)   ! not used for mixing yet
  if (mechanism == 'WaveDrag') TKE_loss_sum = CS%tot_itidal_loss(i,j) ! currently used for mixing
  if (mechanism == 'Froude')   TKE_loss_sum = CS%tot_Froude_loss(i,j) ! not used for mixing yet

end subroutine get_lowmode_loss

!> Implements refraction on the internal waves at a single frequency.
subroutine refract(En, cn, freq, dt, G, US, NAngle, use_PPMang)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure.
  integer,               intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                         intent(inout) :: En   !< The internal gravity wave energy density as a
                                               !! function of space and angular resolution,
                                               !! [R Z3 T-2 ~> J m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed),        &
                         intent(in)    :: cn   !< Baroclinic mode speed [L T-1 ~> m s-1].
  real,                  intent(in)    :: freq !< Wave frequency [T-1 ~> s-1].
  real,                  intent(in)    :: dt   !< Time step [T ~> s].
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  logical,               intent(in)    :: use_PPMang !< If true, use PPM for advection rather
                                               !! than upwind.
  ! Local variables
  integer, parameter :: stencil = 2
  real, dimension(SZI_(G),1-stencil:NAngle+stencil) :: &
    En2d
  real, dimension(1-stencil:NAngle+stencil) :: &
    cos_angle, sin_angle
  real, dimension(SZI_(G)) :: &
    Dk_Dt_Kmag, Dl_Dt_Kmag
  real, dimension(SZI_(G),0:nAngle) :: &
    Flux_E
  real, dimension(SZI_(G),SZJ_(G),1-stencil:NAngle+stencil) :: &
    CFL_ang
  real :: f2              ! The squared Coriolis parameter [T-2 ~> s-2].
  real :: favg            ! The average Coriolis parameter at a point [T-1 ~> s-1].
  real :: df_dy, df_dx    ! The x- and y- gradients of the Coriolis parameter [T-1 L-1 ~> s-1 m-1].
  real :: dlnCn_dx        ! The x-gradient of the wave speed divided by itself [L-1 ~> m-1].
  real :: dlnCn_dy        ! The y-gradient of the wave speed divided by itself [L-1 ~> m-1].
  real :: Angle_size, dt_Angle_size, angle
  real :: Ifreq, Kmag2, I_Kmag
  real :: cn_subRO        ! A tiny wave speed to prevent division by zero [L T-1 ~> m s-1]
  integer :: is, ie, js, je, asd, aed, na
  integer :: i, j, a

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; na = size(En,3)
  asd = 1-stencil ; aed = NAngle+stencil

  Ifreq = 1.0 / freq
  cn_subRO = 1e-100*US%m_s_to_L_T  ! The hard-coded value here might need to increase.
  Angle_size = (8.0*atan(1.0)) / (real(NAngle))
  dt_Angle_size = dt / Angle_size

  do A=asd,aed
    angle = (real(A) - 0.5) * Angle_size
    cos_angle(A) = cos(angle) ; sin_angle(A) = sin(angle)
  enddo

  !### There should also be refraction due to cn.grad(grid_orientation).
  CFL_ang(:,:,:) = 0.0
  do j=js,je
  ! Copy En into angle space with halos.
    do a=1,na ; do i=is,ie
      En2d(i,a) = En(i,j,a)
    enddo ; enddo
    do a=asd,0 ; do i=is,ie
      En2d(i,a) = En2d(i,a+NAngle)
      En2d(i,NAngle+stencil+a) = En2d(i,stencil+a)
    enddo ; enddo

  ! Do the refraction.
    do i=is,ie
      f2 = 0.25* ((G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J-1)**2) + &
                 (G%CoriolisBu(I,J-1)**2 + G%CoriolisBu(I-1,J)**2))
      favg = 0.25*((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J-1)) + &
                   (G%CoriolisBu(I,J-1) + G%CoriolisBu(I-1,J)))
      df_dx = 0.5*((G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1)) - &
                    (G%CoriolisBu(I-1,J) + G%CoriolisBu(I-1,J-1))) * G%IdxT(i,j)
      dlnCn_dx = 0.5*( G%IdxCu(I,j) * (cn(i+1,j) - cn(i,j)) / &
                       (0.5*(cn(i+1,j) + cn(i,j)) + cn_subRO) + &
                       G%IdxCu(I-1,j) * (cn(i,j) - cn(i-1,j)) / &
                       (0.5*(cn(i,j) + cn(i-1,j)) + cn_subRO) )
      df_dy = 0.5*((G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J)) - &
                   (G%CoriolisBu(I,J-1) + G%CoriolisBu(I-1,J-1))) * G%IdyT(i,j)
      dlnCn_dy = 0.5*( G%IdyCv(i,J) * (cn(i,j+1) - cn(i,j)) / &
                       (0.5*(cn(i,j+1) + cn(i,j)) + cn_subRO) + &
                       G%IdyCv(i,J-1) * (cn(i,j) - cn(i,j-1)) / &
                       (0.5*(cn(i,j) + cn(i,j-1)) + cn_subRO) )
      Kmag2 = (freq**2 - f2) / (cn(i,j)**2 + cn_subRO**2)
      if (Kmag2 > 0.0) then
        I_Kmag = 1.0 / sqrt(Kmag2)
        Dk_Dt_Kmag(i) = -Ifreq * (favg*df_dx + (freq**2 - f2) * dlnCn_dx) * I_Kmag
        Dl_Dt_Kmag(i) = -Ifreq * (favg*df_dy + (freq**2 - f2) * dlnCn_dy) * I_Kmag
      else
        Dk_Dt_Kmag(i) = 0.0
        Dl_Dt_Kmag(i) = 0.0
      endif
    enddo

    ! Determine the energy fluxes in angular orientation space.
    do A=asd,aed ; do i=is,ie
      CFL_ang(i,j,A) = (cos_angle(A) * Dl_Dt_Kmag(i) - sin_angle(A) * Dk_Dt_Kmag(i)) * dt_Angle_size
      if (abs(CFL_ang(i,j,A)) > 1.0) then
        call MOM_error(WARNING, "refract: CFL exceeds 1.", .true.)
        if (CFL_ang(i,j,A) > 0.0) then ; CFL_ang(i,j,A) = 1.0 ; else ; CFL_ang(i,j,A) = -1.0 ; endif
      endif
    enddo ; enddo

    ! Advect in angular space
    if (.not.use_PPMang) then
      ! Use simple upwind
      do  A=0,na ; do i=is,ie
        if (CFL_ang(i,j,A) > 0.0) then
          Flux_E(i,A) = CFL_ang(i,j,A) * En2d(i,A)
        else
          Flux_E(i,A) = CFL_ang(i,j,A) * En2d(i,A+1)
        endif
      enddo ; enddo
    else
      ! Use PPM
      do i=is,ie
        call PPM_angular_advect(En2d(i,:),CFL_ang(i,j,:),Flux_E(i,:),NAngle,dt,stencil)
      enddo
    endif

  ! Update and copy back to En.
    do a=1,na ; do i=is,ie
      !if (En2d(i,a)+(Flux_E(i,A-1)-Flux_E(i,A)) < 0.0) then ! for debugging
      !  call MOM_error(FATAL, "refract: OutFlux>Available")
      !endif
      En(i,j,a) = En2d(i,a) + (Flux_E(i,A-1) - Flux_E(i,A))
    enddo ; enddo
  enddo ! j-loop
end subroutine refract

!> This subroutine calculates the 1-d flux for advection in angular space using a monotonic
!! piecewise parabolic scheme. This needs to be called from within i and j spatial loops.
subroutine PPM_angular_advect(En2d, CFL_ang, Flux_En, NAngle, dt, halo_ang)
  integer,                   intent(in)    :: NAngle  !< The number of wave orientations in the
                                                      !! discretized wave energy spectrum.
  real,                      intent(in)    :: dt      !< Time increment [T ~> s].
  integer,                   intent(in)    :: halo_ang !< The halo size in angular space
  real, dimension(1-halo_ang:NAngle+halo_ang),   &
                             intent(in)    :: En2d    !< The internal gravity wave energy density as a
                                                      !! function of angular resolution [R Z3 T-2 ~> J m-2].
  real, dimension(1-halo_ang:NAngle+halo_ang),   &
                             intent(in)    :: CFL_ang !< The CFL number of the energy advection across angles
  real, dimension(0:NAngle), intent(out)   :: Flux_En !< The time integrated internal wave energy flux
                                                      !! across angles  [R Z3 T-2 ~> J m-2].
  ! Local variables
  real :: flux
  real :: u_ang        ! Angular propagation speed [Rad T-1 ~> Rad s-1]
  real :: Angle_size   ! The size of each orientation wedge in radians [Rad]
  real :: I_Angle_size ! The inverse of the the orientation wedges [Rad-1]
  real :: I_dt         ! The inverse of the timestep [T-1 ~> s-1]
  real :: aR, aL       ! Left and right edge estimates of energy density [R Z3 T-2 rad-1 ~> J m-2 rad-1]
  real :: dMx, dMn
  real :: Ep, Ec, Em   ! Mean angular energy density for three successive wedges in angular
                       ! orientation [R Z3 T-2 rad-1 ~> J m-2 rad-1]
  real :: dA, curv_3   ! Difference and curvature of energy density [R Z3 T-2 rad-1 ~> J m-2 rad-1]
  real, parameter :: oneSixth = 1.0/6.0  ! One sixth [nondim]
  integer :: a

  I_dt = 1 / dt
  Angle_size = (8.0*atan(1.0)) / (real(NAngle))
  I_Angle_size = 1 / Angle_size
  Flux_En(:) = 0

  do A=0,NAngle
    u_ang = CFL_ang(A)*Angle_size*I_dt
    if (u_ang >= 0.0) then
      ! Implementation of PPM-H3
      ! Convert wedge-integrated energy density into angular energy densities for three successive
      ! wedges around the source wedge for this flux [R Z3 T-2 rad-1 ~> J m-2 rad-1].
      Ep = En2d(a+1)*I_Angle_size
      Ec = En2d(a)  *I_Angle_size
      Em = En2d(a-1)*I_Angle_size
      ! Calculate and bound edge values of energy density.
      aL = ( 5.*Ec + ( 2.*Em - Ep ) ) * oneSixth ! H3 estimate
      aL = max( min(Ec,Em), aL) ; aL = min( max(Ec,Em), aL) ! Bound
      aR = ( 5.*Ec + ( 2.*Ep - Em ) ) * oneSixth ! H3 estimate
      aR = max( min(Ec,Ep), aR) ; aR = min( max(Ec,Ep), aR) ! Bound
      dA = aR - aL
      if ((Ep-Ec)*(Ec-Em) <= 0.) then
        aL = Ec ; aR = Ec    ! use PCM for local extremum
      elseif ( 3.0*dA*(2.*Ec - (aR + aL)) > (dA*dA) ) then
        aL = 3.*Ec - 2.*aR   ! Flatten the profile to move the extremum to the left edge
      elseif ( 3.0*dA*(2.*Ec - (aR + aL)) < - (dA*dA) ) then
        aR = 3.*Ec - 2.*aL   ! Flatten the profile to move the extremum to the right edge
      endif
      curv_3 = (aR + aL) - 2.0*Ec ! Curvature
      ! Calculate angular flux rate [R Z3 T-3 ~> W m-2]
      flux = u_ang*( aR + CFL_ang(A) * ( 0.5*(aL - aR) + curv_3 * (CFL_ang(A) - 1.5) ) )
      ! Calculate amount of energy fluxed between wedges [R Z3 T-2 ~> J m-2]
      Flux_En(A) = dt * flux
      !Flux_En(A) = (dt * I_Angle_size) * flux
    else
      ! Implementation of PPM-H3
      ! Convert wedge-integrated energy density into angular energy densities for three successive
      ! wedges around the source wedge for this flux [R Z3 T-2 rad-1 ~> J m-2 rad-1].
      Ep = En2d(a+2)*I_Angle_size
      Ec = En2d(a+1)*I_Angle_size
      Em = En2d(a)  *I_Angle_size
      ! Calculate and bound edge values of energy density.
      aL = ( 5.*Ec + ( 2.*Em - Ep ) ) * oneSixth ! H3 estimate
      aL = max( min(Ec,Em), aL) ; aL = min( max(Ec,Em), aL) ! Bound
      aR = ( 5.*Ec + ( 2.*Ep - Em ) ) * oneSixth ! H3 estimate
      aR = max( min(Ec,Ep), aR) ; aR = min( max(Ec,Ep), aR) ! Bound
      dA = aR - aL
      if ((Ep-Ec)*(Ec-Em) <= 0.) then
        aL = Ec ; aR = Ec    ! use PCM for local extremum
      elseif ( 3.0*dA*(2.*Ec - (aR + aL)) > (dA*dA) ) then
        aL = 3.*Ec - 2.*aR   ! Flatten the profile to move the extremum to the left edge
      elseif ( 3.0*dA*(2.*Ec - (aR + aL)) < - (dA*dA) ) then
        aR = 3.*Ec - 2.*aL   ! Flatten the profile to move the extremum to the right edge
      endif
      curv_3 = (aR + aL) - 2.0*Ec ! Curvature
      ! Calculate angular flux rate [R Z3 T-3 ~> W m-2]
      ! Note that CFL_ang is negative here, so it looks odd compared with equivalent expressions.
      flux = u_ang*( aL - CFL_ang(A) * ( 0.5*(aR - aL) + curv_3 * (-CFL_ang(A) - 1.5) ) )
      ! Calculate amount of energy fluxed between wedges [R Z3 T-2 ~> J m-2]
      Flux_En(A) = dt * flux
      !Flux_En(A) = (dt * I_Angle_size) * flux
    endif
  enddo
end subroutine PPM_angular_advect

!> Propagates internal waves at a single frequency.
subroutine propagate(En, cn, freq, dt, G, US, CS, NAngle)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure.
  integer,               intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                         intent(inout) :: En   !< The internal gravity wave energy density as a
                                               !! function of space and angular resolution,
                                               !! [R Z3 T-2 ~> J m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed),        &
                         intent(in)    :: cn   !< Baroclinic mode speed [L T-1 ~> m s-1].
  real,                  intent(in)    :: freq !< Wave frequency [T-1 ~> s-1].
  real,                  intent(in)    :: dt   !< Time step [T ~> s].
  type(unit_scale_type), intent(in)    :: US   !< A dimensional unit scaling type
  type(int_tide_CS),     pointer       :: CS   !< The control structure returned by a
                                               !! previous call to int_tide_init.
  ! Local variables
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB) :: &
    speed  ! The magnitude of the group velocity at the q points for corner adv [L T-1 ~> m s-1].
  integer, parameter :: stencil = 2
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    speed_x  ! The magnitude of the group velocity at the Cu points [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    speed_y  ! The magnitude of the group velocity at the Cv points [L T-1 ~> m s-1].
  real, dimension(0:NAngle) :: &
    cos_angle, sin_angle
  real, dimension(NAngle) :: &
    Cgx_av, Cgy_av, dCgx, dCgy
  real :: f2   ! The squared Coriolis parameter [T-2 ~> s-2].
  real :: Angle_size, I_Angle_size, angle
  real :: Ifreq ! The inverse of the frequency [T ~> s]
  real :: freq2 ! The frequency squared [T-2 ~> s-2]
  type(loop_bounds_type) :: LB
  integer :: is, ie, js, je, asd, aed, na
  integer :: ish, ieh, jsh, jeh
  integer :: i, j, a

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; na = size(En,3)
  asd = 1-stencil ; aed = NAngle+stencil

  Ifreq = 1.0 / freq
  freq2 = freq**2

  ! Define loop bounds: Need extensions on j-loop so propagate_y
  ! (done after propagate_x) will have updated values in the halo
  ! for correct PPM reconstruction. Use if no teleporting and
  ! no pass_var between propagate_x and propagate_y.
  !jsh = js-3 ; jeh = je+3 ; ish = is ; ieh = ie

  ! Define loop bounds: Need 1-pt extensions on loops because
  ! teleporting eats up a halo point. Use if teleporting.
  ! Also requires pass_var before propagate_y.
  jsh = js-1 ; jeh = je+1 ; ish = is-1 ; ieh = ie+1

  Angle_size = (8.0*atan(1.0)) / real(NAngle)
  I_Angle_size = 1.0 / Angle_size

  if (CS%corner_adv) then
    ! IMPLEMENT CORNER ADVECTION IN HORIZONTAL--------------------
    ! FIND AVERAGE GROUP VELOCITY (SPEED) AT CELL CORNERS
    ! NOTE: THIS HAS NOT BE ADAPTED FOR REFLECTION YET (BDM)!!
    ! Fix indexing here later
    speed(:,:) = 0.0
    do J=jsh-1,jeh ; do I=ish-1,ieh
      f2 = G%CoriolisBu(I,J)**2
      speed(I,J) = 0.25*(cn(i,j) + cn(i+1,j) + cn(i+1,j+1) + cn(i,j+1)) * &
                     sqrt(max(freq2 - f2, 0.0)) * Ifreq
    enddo ; enddo
    do a=1,na
      ! Apply the propagation WITH CORNER ADVECTION/FINITE VOLUME APPROACH.
      LB%jsh = js ; LB%jeh = je ; LB%ish = is ; LB%ieh = ie
      call propagate_corner_spread(En(:,:,a), a, NAngle, speed, dt, G, CS, LB)
    enddo ! a-loop
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
      f2 = 0.5 * (G%CoriolisBu(I,J)**2 + G%CoriolisBu(I,J-1)**2)
      speed_x(I,j) = 0.5*(cn(i,j) + cn(i+1,j)) * G%mask2dCu(I,j) * &
                     sqrt(max(freq2 - f2, 0.0)) * Ifreq
    enddo ; enddo
    do J=jsh-1,jeh ; do i=ish,ieh
      f2 = 0.5 * (G%CoriolisBu(I,J)**2 + G%CoriolisBu(I-1,J)**2)
      speed_y(i,J) = 0.5*(cn(i,j) + cn(i,j+1)) * G%mask2dCv(i,J) * &
                     sqrt(max(freq2 - f2, 0.0)) * Ifreq
    enddo ; enddo

    ! Apply propagation in x-direction (reflection included)
    LB%jsh = jsh ; LB%jeh = jeh ; LB%ish = ish ; LB%ieh = ieh
    call propagate_x(En, speed_x, Cgx_av, dCgx, dt, G, US, CS%nAngle, CS, LB)

    ! Check for energy conservation on computational domain (for debugging)
    !call sum_En(G, CS, En, 'post-propagate_x')

    ! Update halos
    call pass_var(En, G%domain)

    ! Apply propagation in y-direction (reflection included)
    ! LB%jsh = js ; LB%jeh = je ; LB%ish = is ; LB%ieh = ie ! Use if no teleport
    LB%jsh = jsh ; LB%jeh = jeh ; LB%ish = ish ; LB%ieh = ieh
    call propagate_y(En, speed_y, Cgy_av, dCgy, dt, G, US, CS%nAngle, CS, LB)

    ! Check for energy conservation on computational domain (for debugging)
    !call sum_En(G, CS, En, 'post-propagate_y')
  endif

end subroutine propagate

!> This subroutine does first-order corner advection. It was written with the hopes
!! of smoothing out the garden sprinkler effect, but is too numerically diffusive to
!! be of much use as of yet. It is not yet compatible with reflection schemes (BDM).
subroutine propagate_corner_spread(En, energized_wedge, NAngle, speed, dt, G, CS, LB)
  type(ocean_grid_type),  intent(in)    :: G     !< The ocean's grid structure.
  real, dimension(G%isd:G%ied,G%jsd:G%jed),   &
                          intent(inout) :: En    !< The energy density integrated over an angular
                                                 !! band [R Z3 T-2 ~> J m-2].
  real, dimension(G%IsdB:G%IedB,G%Jsd:G%Jed), &
                          intent(in)    :: speed !< The magnitude of the group velocity at the cell
                                                 !! corner points [L T-1 ~> m s-1].
  integer,                intent(in)    :: energized_wedge !< Index of current ray direction.
  integer,                intent(in)    :: NAngle !< The number of wave orientations in the
                                                 !! discretized wave energy spectrum.
  real,                   intent(in)    :: dt    !< Time increment [T ~> s].
  type(int_tide_CS),      pointer       :: CS    !< The control structure returned by a previous
                                                 !! call to continuity_PPM_init.
  type(loop_bounds_type), intent(in)    :: LB    !< A structure with the active energy loop bounds.
  ! Local variables
  integer :: i, j, k, ish, ieh, jsh, jeh, m
  real :: TwoPi, Angle_size
  real :: energized_angle ! angle through center of current wedge
  real :: theta ! angle at edge of wedge
  real :: Nsubrays     ! number of sub-rays for averaging
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
  real, dimension(2) :: E_new ! Energy in cell after advection for subray [R Z3 T-2 ~> J m-2]; set size
                              ! here to define Nsubrays - this should be made an input option later!

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh
  TwoPi = (8.0*atan(1.0))
  Nsubrays = real(size(E_new))
  I_Nsubwedges = 1./(Nsubrays - 1)

  Angle_size = TwoPi / real(NAngle)
  energized_angle = Angle_size * real(energized_wedge - 1) ! for a=1 aligned with x-axis
  !energized_angle = Angle_size * real(energized_wedge - 1) + 2.0*Angle_size !
  !energized_angle = Angle_size * real(energized_wedge - 1) + 0.5*Angle_size !
  do J=jsh-1,jeh ; do I=ish-1,ieh
    ! This will only work for a Cartesian grid for which G%geoLonBu is in the same units has dx.
    ! This needs to be extensively revised to work for a general grid.
    x(I,J) = G%US%m_to_L*G%geoLonBu(I,J)
    y(I,J) = G%US%m_to_L*G%geoLatBu(I,J)
    Idx(I,J) = G%IdxBu(I,J) ; dx(I,J) = G%dxBu(I,J)
    Idy(I,J) = G%IdyBu(I,J) ; dy(I,J) = G%dyBu(I,J)
  enddo ; enddo

  do j=jsh,jeh ; do i=ish,ieh
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
          xCrn = x(I-1,J-1); yCrn = y(I-1,J-1)
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
          xCrn = x(I,J-1); yCrn = y(I,J-1)
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
          xCrn = x(I,J); yCrn = y(I,J)
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
          xCrn = x(I-1,J); yCrn = y(I-1,J)
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
          aNW = a1 + a2 + a3 + a4
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
  enddo ; enddo
end subroutine propagate_corner_spread

!> Propagates the internal wave energy in the logical x-direction.
subroutine propagate_x(En, speed_x, Cgx_av, dCgx, dt, G, US, Nangle, CS, LB)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  integer,                 intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle),   &
                           intent(inout) :: En !< The energy density integrated over an angular
                                               !! band [R Z3 T-2 ~> J m-2].
  real, dimension(G%IsdB:G%IedB,G%jsd:G%jed),        &
                           intent(in)    :: speed_x !< The magnitude of the group velocity at the
                                               !! Cu points [L T-1 ~> m s-1].
  real, dimension(Nangle), intent(in)    :: Cgx_av !< The average x-projection in each angular band.
  real, dimension(Nangle), intent(in)    :: dCgx !< The difference in x-projections between the
                                               !! edges of each angular band.
  real,                    intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(int_tide_CS),       pointer       :: CS !< The control structure returned by a previous call
                                               !! to continuity_PPM_init.
  type(loop_bounds_type),  intent(in)    :: LB !< A structure with the active energy loop bounds.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    EnL, EnR    ! Left and right face energy densities [R Z3 T-2 ~> J m-2].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    flux_x      ! The internal wave energy flux [J T-1 ~> J s-1].
  real, dimension(SZIB_(G)) :: &
    cg_p, cg_m, flux1, flux2
  !real, dimension(SZI_(G),SZJB_(G),Nangle) :: En_m, En_p
  real, dimension(SZI_(G),SZJB_(G),Nangle) :: &
    Fdt_m, Fdt_p! Left and right energy fluxes [J]
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
                         dt, G, US, j, ish, ieh, CS%vol_CFL)
      do I=ish-1,ieh ; flux_x(I,j) = flux1(I); enddo
    enddo

    do j=jsh,jeh ; do i=ish,ieh
      Fdt_m(i,j,a) = dt*flux_x(I-1,j) ! left face influx  (J)
      Fdt_p(i,j,a) = -dt*flux_x(I,j)  ! right face influx (J)
    enddo ; enddo

  enddo ! a-loop

  ! Only reflect newly arrived energy; existing energy in incident wedge is not reflected
  ! and will eventually propagate out of cell. (This code only reflects if En > 0.)
  call reflect(Fdt_m, Nangle, CS, G, LB)
  call teleport(Fdt_m, Nangle, CS, G, LB)
  call reflect(Fdt_p, Nangle, CS, G, LB)
  call teleport(Fdt_p, Nangle, CS, G, LB)

  ! Update reflected energy [R Z3 T-2 ~> J m-2]
  do a=1,Nangle ; do j=jsh,jeh ; do i=ish,ieh
    !  if ((En(i,j,a) + G%IareaT(i,j)*(Fdt_m(i,j,a) + Fdt_p(i,j,a))) < 0.0) & ! for debugging
    !    call MOM_error(FATAL, "propagate_x: OutFlux>Available")
    En(i,j,a) = En(i,j,a) + G%IareaT(i,j)*(Fdt_m(i,j,a) + Fdt_p(i,j,a))
  enddo ; enddo ; enddo

end subroutine propagate_x

!> Propagates the internal wave energy in the logical y-direction.
subroutine propagate_y(En, speed_y, Cgy_av, dCgy, dt, G, US, Nangle, CS, LB)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  integer,                 intent(in)    :: NAngle !< The number of wave orientations in the
                                               !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,Nangle), &
                           intent(inout) :: En !< The energy density integrated over an angular
                                               !! band [R Z3 T-2 ~> J m-2].
  real, dimension(G%isd:G%ied,G%JsdB:G%JedB),      &
                           intent(in)    :: speed_y !< The magnitude of the group velocity at the
                                               !! Cv points [L T-1 ~> m s-1].
  real, dimension(Nangle), intent(in)    :: Cgy_av !< The average y-projection in each angular band.
  real, dimension(Nangle), intent(in)    :: dCgy !< The difference in y-projections between the
                                               !! edges of each angular band.
  real,                    intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US !< A dimensional unit scaling type
  type(int_tide_CS),       pointer       :: CS !< The control structure returned by a previous call
                                               !! to continuity_PPM_init.
  type(loop_bounds_type),  intent(in)    :: LB !< A structure with the active energy loop bounds.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    EnL, EnR    ! South and north face energy densities [R Z3 T-2 ~> J m-2].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    flux_y      ! The internal wave energy flux [J T-1 ~> J s-1].
  real, dimension(SZI_(G)) :: &
    cg_p, cg_m, flux1, flux2
  !real, dimension(SZI_(G),SZJB_(G),Nangle) :: En_m, En_p
  real, dimension(SZI_(G),SZJB_(G),Nangle) :: &
    Fdt_m, Fdt_p! South and north energy fluxes [J]
  character(len=160) :: mesg  ! The text of an error message
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
                         dt, G, US, J, ish, ieh, CS%vol_CFL)
      do i=ish,ieh ; flux_y(i,J) = flux1(i); enddo
    enddo

    do j=jsh,jeh ; do i=ish,ieh
      Fdt_m(i,j,a) = dt*flux_y(i,J-1) ! south face influx (J)
      Fdt_p(i,j,a) = -dt*flux_y(i,J)  ! north face influx (J)
      !if ((En(i,j,a) + G%IareaT(i,j)*(Fdt_m(i,j,a) + Fdt_p(i,j,a))) < 0.0) then ! for debugging
      !  call MOM_error(WARNING, "propagate_y: OutFlux>Available prior to reflection", .true.)
      !  write(mesg,*) "flux_y_south=",flux_y(i,J-1),"flux_y_north=",flux_y(i,J),"En=",En(i,j,a), &
      !           "cn_south=", speed_y(i,J-1) * (Cgy_av(a)), "cn_north=", speed_y(i,J) * (Cgy_av(a))
      !  call MOM_error(WARNING, mesg, .true.)
      !endif
    enddo ; enddo

  enddo ! a-loop

  ! Only reflect newly arrived energy; existing energy in incident wedge is not reflected
  ! and will eventually propagate out of cell. (This code only reflects if En > 0.)
  call reflect(Fdt_m, Nangle, CS, G, LB)
  call teleport(Fdt_m, Nangle, CS, G, LB)
  call reflect(Fdt_p, Nangle, CS, G, LB)
  call teleport(Fdt_p, Nangle, CS, G, LB)

  ! Update reflected energy [R Z3 T-2 ~> J m-2]
  do a=1,Nangle ; do j=jsh,jeh ; do i=ish,ieh
    !  if ((En(i,j,a) + G%IareaT(i,j)*(Fdt_m(i,j,a) + Fdt_p(i,j,a))) < 0.0) & ! for debugging
    !    call MOM_error(FATAL, "propagate_y: OutFlux>Available", .true.)
    En(i,j,a) = En(i,j,a) + G%IareaT(i,j)*(Fdt_m(i,j,a) + Fdt_p(i,j,a))
  enddo ; enddo ; enddo

end subroutine propagate_y

!> Evaluates the zonal mass or volume fluxes in a layer.
subroutine zonal_flux_En(u, h, hL, hR, uh, dt, G, US, j, ish, ieh, vol_CFL)
  type(ocean_grid_type),     intent(in)    :: G  !< The ocean's grid structure.
  real, dimension(SZIB_(G)), intent(in)    :: u  !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G)),  intent(in)    :: h  !< Energy density used to calculate the fluxes
                                                 !! [R Z3 T-2 ~> J m-2].
  real, dimension(SZI_(G)),  intent(in)    :: hL !< Left- Energy densities in the reconstruction
                                                 !! [R Z3 T-2 ~> J m-2].
  real, dimension(SZI_(G)),  intent(in)    :: hR !< Right- Energy densities in the reconstruction
                                                 !! [R Z3 T-2 ~> J m-2].
  real, dimension(SZIB_(G)), intent(inout) :: uh !< The zonal energy transport [R Z3 L2 T-3 ~> J s-1].
  real,                      intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),     intent(in)    :: US !< A dimensional unit scaling type
  integer,                   intent(in)    :: j  !< The j-index to work on.
  integer,                   intent(in)    :: ish !< The start i-index range to work on.
  integer,                   intent(in)    :: ieh !< The end i-index range to work on.
  logical,                   intent(in)    :: vol_CFL !< If true, rescale the ratio of face areas to
                                                 !! the cell areas when estimating the CFL number.
  ! Local variables
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim].
  real :: curv_3 ! A measure of the energy density curvature over a grid length [R Z3 T-2 ~> J m-2]
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

!> Evaluates the meridional mass or volume fluxes in a layer.
subroutine merid_flux_En(v, h, hL, hR, vh, dt, G, US, J, ish, ieh, vol_CFL)
  type(ocean_grid_type),            intent(in)    :: G  !< The ocean's grid structure.
  real, dimension(SZI_(G)),         intent(in)    :: v  !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: h  !< Energy density used to calculate the
                                                        !! fluxes [R Z3 T-2 ~> J m-2].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: hL !< Left- Energy densities in the
                                                        !! reconstruction [R Z3 T-2 ~> J m-2].
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: hR !< Right- Energy densities in the
                                                        !! reconstruction [R Z3 T-2 ~> J m-2].
  real, dimension(SZI_(G)),         intent(inout) :: vh !< The meridional energy transport [R Z3 L2 T-3 ~> J s-1].
  real,                             intent(in)    :: dt !< Time increment [T ~> s].
  type(unit_scale_type),            intent(in)    :: US !< A dimensional unit scaling type
  integer,                          intent(in)    :: J  !< The j-index to work on.
  integer,                          intent(in)    :: ish !< The start i-index range to work on.
  integer,                          intent(in)    :: ieh !< The end i-index range to work on.
  logical,                          intent(in)    :: vol_CFL !< If true, rescale the ratio of face
                                                        !! areas to the cell areas when estimating
                                                        !! the CFL number.
  ! Local variables
  real :: CFL ! The CFL number based on the local velocity and grid spacing [nondim].
  real :: curv_3 ! A measure of the energy density curvature over a grid length [R Z3 T-2 ~> J m-2]
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

!> Reflection of the internal waves at a single frequency.
subroutine reflect(En, NAngle, CS, G, LB)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure
  integer,                intent(in)    :: NAngle !< The number of wave orientations in the
                                              !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                          intent(inout) :: En !< The internal gravity wave energy density as a
                                              !! function of space and angular resolution
                                              !! [R Z3 T-2 ~> J m-2].
  type(int_tide_CS),      pointer       :: CS !< The control structure returned by a
                                              !! previous call to int_tide_init.
  type(loop_bounds_type), intent(in)    :: LB !< A structure with the active energy loop bounds.
  ! Local variables
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: angle_c
                                           ! angle of boundary wrt equator [rad]
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: part_refl
                                           ! fraction of wave energy reflected
                                           ! values should collocate with angle_c [nondim]
  logical, dimension(G%isd:G%ied,G%jsd:G%jed) :: ridge
                                           ! tags of cells with double reflection

  real    :: TwoPi                         ! 2*pi = 6.2831853... [nondim]
  real    :: Angle_size                    ! size of beam wedge [rad]
  real    :: angle_wall                    ! angle of coast/ridge/shelf wrt equator [rad]
  real, dimension(1:NAngle) :: angle_i     ! angle of incident ray wrt equator [rad]
  real    :: angle_r                       ! angle of reflected ray wrt equator [rad]
  real, dimension(1:Nangle) :: En_reflected
  integer :: i, j, a, a_r, na
  !integer :: isd, ied, jsd, jed   ! start and end local indices on data domain
  !                                ! (values include halos)
  integer :: isc, iec, jsc, jec   ! start and end local indices on PE
                                  ! (values exclude halos)
  integer :: ish, ieh, jsh, jeh   ! start and end local indices on data domain
                                  ! leaving out outdated halo points (march in)

  !isd = G%isd  ; ied = G%ied  ; jsd = G%jsd  ; jed = G%jed
  isc = G%isc  ; iec = G%iec  ; jsc = G%jsc  ; jec = G%jec
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh

  TwoPi = 8.0*atan(1.0)
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

  do j=jsh,jeh ; do i=ish,ieh
    ! redistribute energy in angular space if ray will hit boundary
    ! i.e., if energy is in a reflecting cell
    if (angle_c(i,j) /= CS%nullangle) then
      do a=1,NAngle ; if (En(i,j,a) > 0.0) then
        if (sin(angle_i(a) - angle_c(i,j)) >= 0.0) then
          ! if ray is incident, keep specified boundary angle
          angle_wall = angle_c(i,j)
        elseif (ridge(i,j)) then
         ! if ray is not incident but in ridge cell, use complementary angle
         angle_wall = angle_c(i,j) + 0.5*TwoPi
          if (angle_wall > TwoPi) then
            angle_wall = angle_wall - TwoPi*floor(abs(angle_wall)/TwoPi)
          elseif (angle_wall < 0.0) then
            angle_wall = angle_wall + TwoPi*ceiling(abs(angle_wall)/TwoPi)
          endif
        else
          ! if ray is not incident and not in a ridge cell, keep specified angle
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
          if (a /= a_r) then
            En_reflected(a_r) = part_refl(i,j)*En(i,j,a)
            En(i,j,a)   = (1.0-part_refl(i,j))*En(i,j,a)
          endif
        endif
      endif ; enddo ! a-loop
      do a=1,NAngle
        En(i,j,a) = En(i,j,a) + En_reflected(a)
        En_reflected(a) = 0.0
      enddo ! a-loop
    endif
  enddo ; enddo ! i- and j-loops

  ! Check to make sure no energy gets onto land (only run for debugging)
  ! do a=1,NAngle ; do j=jsc,jec ; do i=isc,iec
  !   if (En(i,j,a) > 0.001 .and. G%mask2dT(i,j) == 0) then
  !     write (mesg,*) 'En=', En(i,j,a), 'a=', a, 'ig_g=',i+G%idg_offset, 'jg_g=',j+G%jdg_offset
  !     call MOM_error(FATAL, "reflect: Energy detected out of bounds: "//trim(mesg), .true.)
  !   endif
  ! enddo ; enddo ; enddo

end subroutine reflect

!> Moves energy across lines of partial reflection to prevent
!! reflection of energy that is supposed to get across.
subroutine teleport(En, NAngle, CS, G, LB)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure.
  integer,                intent(in)    :: NAngle !< The number of wave orientations in the
                                              !! discretized wave energy spectrum.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,NAngle), &
                          intent(inout) :: En !< The internal gravity wave energy density as a
                                              !! function of space and angular resolution
                                              !! [R Z3 T-2 ~> J m-2].
  type(int_tide_CS),      pointer       :: CS !< The control structure returned by a
                                              !! previous call to int_tide_init.
  type(loop_bounds_type), intent(in)    :: LB !< A structure with the active energy loop bounds.
  ! Local variables
  real, dimension(G%isd:G%ied,G%jsd:G%jed)    :: angle_c
                                              ! angle of boundary wrt equator [rad]
  real, dimension(G%isd:G%ied,G%jsd:G%jed)    :: part_refl
                                              ! fraction of wave energy reflected
                                              ! values should collocate with angle_c [nondim]
  logical, dimension(G%isd:G%ied,G%jsd:G%jed) :: pref_cell
                                              ! flag for partial reflection
  logical, dimension(G%isd:G%ied,G%jsd:G%jed) :: ridge
                                              ! tags of cells with double reflection
  real                        :: TwoPi      ! 2*pi = 6.2831853... [nondim]
  real                        :: Angle_size ! size of beam wedge [rad]
  real, dimension(1:NAngle)   :: angle_i    ! angle of incident ray wrt equator [rad]
  real, dimension(1:NAngle)   :: cos_angle, sin_angle
  real                        :: En_tele    ! energy to be "teleported" [R Z3 T-2 ~> J m-2]
  character(len=160) :: mesg  ! The text of an error message
  integer :: i, j, a
  !integer :: isd, ied, jsd, jed    ! start and end local indices on data domain
  !                                 ! (values include halos)
  !integer :: isc, iec, jsc, jec    ! start and end local indices on PE
  !                                 ! (values exclude halos)
  integer :: ish, ieh, jsh, jeh     ! start and end local indices on data domain
                                    ! leaving out outdated halo points (march in)
  integer :: id_g, jd_g             ! global (decomp-invar) indices
  integer :: jos, ios               ! offsets
  real    :: cos_normal, sin_normal, angle_wall
                                    ! cos/sin of cross-ridge normal, ridge angle

  !isd = G%isd  ; ied = G%ied  ; jsd = G%jsd  ; jed = G%jed
  !isc = G%isc  ; iec = G%iec  ; jsc = G%jsc  ; jec = G%jec
  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh

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
    do i=ish,ieh
      id_g = i + G%idg_offset ; jd_g = j + G%jdg_offset
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
                write(mesg,*) 'idg=',id_g,'jd_g=',jd_g,'a=',a
                call MOM_error(FATAL, "teleport: no receptive ocean cell at "//trim(mesg), .true.)
              endif
            endif ! incidence check
          endif ! energy check
        enddo ! a-loop
      endif ! pref check
    enddo ! i-loop
  enddo ! j-loop

end subroutine teleport

!> Rotates points in the halos where required to accommodate
!! changes in grid orientation, such as at the tripolar fold.
subroutine correct_halo_rotation(En, test, G, NAngle)
  type(ocean_grid_type),      intent(in)    :: G    !< The ocean's grid structure
  real, dimension(:,:,:,:,:), intent(inout) :: En   !< The internal gravity wave energy density as a
                                       !! function of space, angular orientation, frequency,
                                       !! and vertical mode [R Z3 T-2 ~> J m-2].
  real, dimension(SZI_(G),SZJ_(G),2), &
                              intent(in)    :: test !< An x-unit vector that has been passed through
                                       !! the halo updates, to enable the rotation of the
                                       !! wave energies in the halo region to be corrected.
  integer,                    intent(in)    :: NAngle !< The number of wave orientations in the
                                                      !! discretized wave energy spectrum.
  ! Local variables
  real, dimension(G%isd:G%ied,NAngle) :: En2d
  integer, dimension(G%isd:G%ied) :: a_shift
  integer :: i_first, i_last, a_new
  integer :: a, i, j, isd, ied, jsd, jed, m, fr
  character(len=160) :: mesg  ! The text of an error message
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

!> Calculates left/right edge values for PPM reconstruction in x-direction.
subroutine PPM_reconstruction_x(h_in, h_l, h_r, G, LB, simple_2nd)
  type(ocean_grid_type),            intent(in)  :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Energy density in a sector (2D).
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l  !< Left edge value of reconstruction (2D).
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r  !< Right edge value of reconstruction (2D).
  type(loop_bounds_type),           intent(in)  :: LB   !< A structure with the active loop bounds.
  logical, optional,                intent(in)  :: simple_2nd !< If true, use the arithmetic mean
                                                        !! energy densities as default edge values
                                                        !! for a simple 2nd order scheme.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_ip1, h_im1
  real :: dMx, dMn
  logical :: use_CW84, use_2nd
  character(len=256) :: mesg  ! The text of an error message
  integer :: i, j, isl, iel, jsl, jel, stencil

  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  isl = LB%ish-1 ; iel = LB%ieh+1 ; jsl = LB%jsh ; jel = LB%jeh

  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (use_2nd) stencil = 1

  if ((isl-stencil < G%isd) .or. (iel+stencil > G%ied)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_x called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               stencil + max(G%isd-isl,iel-G%ied)
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
    enddo ; enddo

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
    enddo ; enddo
  endif

  call PPM_limit_pos(h_in, h_l, h_r, 0.0, G, isl, iel, jsl, jel)
end subroutine PPM_reconstruction_x

!> Calculates left/right edge valus for PPM reconstruction in y-direction.
subroutine PPM_reconstruction_y(h_in, h_l, h_r, G, LB, simple_2nd)
  type(ocean_grid_type),            intent(in)  :: G    !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Energy density in a sector (2D).
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l  !< Left edge value of reconstruction (2D).
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r  !< Right edge value of reconstruction (2D).
  type(loop_bounds_type),           intent(in)  :: LB   !< A structure with the active loop bounds.
  logical, optional,                intent(in)  :: simple_2nd !< If true, use the arithmetic mean
                                                        !! energy densities as default edge values
                                                        !! for a simple 2nd order scheme.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_jp1, h_jm1
  real :: dMx, dMn
  logical :: use_2nd
  character(len=256) :: mesg  ! The text of an error message
  integer :: i, j, isl, iel, jsl, jel, stencil

  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  isl = LB%ish ; iel = LB%ieh ; jsl = LB%jsh-1 ; jel = LB%jeh+1

  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (use_2nd) stencil = 1

  if ((isl < G%isd) .or. (iel > G%ied)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_y called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               max(G%isd-isl,iel-G%ied)
    call MOM_error(FATAL,mesg)
  endif
  if ((jsl-stencil < G%jsd) .or. (jel+stencil > G%jed)) then
    write(mesg,'("In MOM_internal_tides, PPM_reconstruction_y called with a ", &
                 & "y-halo that needs to be increased by ",i2,".")') &
                 stencil + max(G%jsd-jsl,jel-G%jed)
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

  call PPM_limit_pos(h_in, h_l, h_r, 0.0, G, isl, iel, jsl, jel)
end subroutine PPM_reconstruction_y

!> Limits the left/right edge values of the PPM reconstruction
!! to give a reconstruction that is positive-definite.  Here this is
!! reinterpreted as giving a constant thickness if the mean thickness is less
!! than h_min, with a minimum of h_min otherwise.
subroutine PPM_limit_pos(h_in, h_L, h_R, h_min, G, iis, iie, jis, jie)
  type(ocean_grid_type),            intent(in)     :: G     !< The ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)), intent(in)     :: h_in  !< Thickness of layer (2D).
  real, dimension(SZI_(G),SZJ_(G)), intent(inout)  :: h_L   !< Left edge value (2D).
  real, dimension(SZI_(G),SZJ_(G)), intent(inout)  :: h_R   !< Right edge value (2D).
  real,                             intent(in)     :: h_min !< The minimum thickness that can be
                                                            !! obtained by a concave parabolic fit.
  integer,                          intent(in)     :: iis   !< Start i-index for computations
  integer,                          intent(in)     :: iie   !< End i-index for computations
  integer,                          intent(in)     :: jis   !< Start j-index for computations
  integer,                          intent(in)     :: jie   !< End j-index for computations
  ! Local variables
  real    :: curv, dh, scale
  character(len=256) :: mesg  ! The text of an error message
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

! subroutine register_int_tide_restarts(G, param_file, CS, restart_CS)
!   type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
!   type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
!   type(int_tide_CS),     pointer       :: CS  !< The control structure returned by a
!                                               !! previous call to int_tide_init.
!   type(MOM_restart_CS),  pointer :: restart_CS !<  A pointer to the restart control structure.

!   ! This subroutine is not currently in use!!

!   ! This subroutine is used to allocate and register any fields in this module
!   ! that should be written to or read from the restart file.
!   logical :: use_int_tides
!   type(vardesc) :: vd
!   integer :: num_freq, num_angle , num_mode, period_1
!   integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, a
!   isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
!   IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

!   if (associated(CS)) then
!     call MOM_error(WARNING, "register_int_tide_restarts called "//&
!                              "with an associated control structure.")
!     return
!   endif

!   use_int_tides = .false.
!   call read_param(param_file, "INTERNAL_TIDES", use_int_tides)
!   if (.not.use_int_tides) return

!   allocate(CS)

!   num_angle = 24
!   call read_param(param_file, "INTERNAL_TIDE_ANGLES", num_angle)
!   allocate(CS%En_restart(isd:ied, jsd:jed, num_angle))
!   CS%En_restart(:,:,:) = 0.0

!   vd = vardesc("En_restart", &
!     "The internal wave energy density as a function of (i,j,angle,frequency,mode)", &
!     'h','1','1',"J m-2")
!   call register_restart_field(CS%En_restart, vd, .false., restart_CS)

! end subroutine register_int_tide_restarts

!> This subroutine initializes the internal tides module.
subroutine internal_tides_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target,   intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),     intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type),   intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),     intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),     intent(in)    :: param_file !< A structure to parse for run-time
                                                   !! parameters.
  type(diag_ctrl), target,   intent(in)    :: diag !< A structure that is used to regulate
                                                   !! diagnostic output.
  type(int_tide_CS),pointer                :: CS   !< A pointer that is set to point to the control
                                                   !! structure for this module.
  ! Local variables
  real                              :: Angle_size ! size of wedges, rad
  real, allocatable                 :: angles(:)  ! orientations of wedge centers, rad
  real, allocatable, dimension(:,:) :: h2         ! topographic roughness scale, m^2
  real                              :: kappa_itides, kappa_h2_factor
                                                  ! characteristic topographic wave number
                                                  ! and a scaling factor
  real, allocatable                 :: ridge_temp(:,:)
                                                  ! array for temporary storage of flags
                                                  ! of cells with double-reflecting ridges
  logical                           :: use_int_tides, use_temperature
  real    :: period_1  ! The period of the gravest modeled mode [T ~> s]
  integer :: num_angle, num_freq, num_mode, m, fr
  integer :: isd, ied, jsd, jed, a, id_ang, i, j
  type(axes_grp) :: axes_ang
  ! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_internal_tides" ! This module's name.
  character(len=16), dimension(8) :: freq_name
  character(len=40)  :: var_name
  character(len=160) :: var_descript
  character(len=200) :: filename
  character(len=200) :: refl_angle_file, land_mask_file
  character(len=200) :: refl_pref_file, refl_dbl_file
  character(len=200) :: dy_Cu_file, dx_Cv_file
  character(len=200) :: h2_file

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
    call MOM_error(WARNING, "internal_tides_init called "//&
                             "with an associated control structure.")
    return
  else
    allocate(CS)
  endif

  use_int_tides = .false.
  call read_param(param_file, "INTERNAL_TIDES", use_int_tides)
  CS%do_int_tides = use_int_tides
  if (.not.use_int_tides) return

  use_temperature = .true.
  call read_param(param_file, "ENABLE_THERMODYNAMICS", use_temperature)
  if (.not.use_temperature) call MOM_error(FATAL, &
    "register_int_tide_restarts: internal_tides only works with "//&
    "ENABLE_THERMODYNAMICS defined.")

  ! Set number of frequencies, angles, and modes to consider
  num_freq = 1 ; num_angle = 24 ; num_mode = 1
  call read_param(param_file, "INTERNAL_TIDE_FREQS", num_freq)
  call read_param(param_file, "INTERNAL_TIDE_ANGLES", num_angle)
  call read_param(param_file, "INTERNAL_TIDE_MODES", num_mode)
  if (.not.((num_freq > 0) .and. (num_angle > 0) .and. (num_mode > 0))) return
  CS%nFreq = num_freq ; CS%nAngle = num_angle ; CS%nMode = num_mode

  ! Allocate energy density array
  allocate(CS%En(isd:ied, jsd:jed, num_angle, num_freq, num_mode))
  CS%En(:,:,:,:,:) = 0.0

  ! Allocate phase speed array
  allocate(CS%cp(isd:ied, jsd:jed, num_freq, num_mode))
  CS%cp(:,:,:,:) = 0.0

  ! Allocate and populate frequency array (each a multiple of first for now)
  allocate(CS%frequency(num_freq))
  call get_param(param_file, mdl, "FIRST_MODE_PERIOD", period_1, units="s", scale=US%s_to_T)
  do fr=1,num_freq
    CS%frequency(fr) = (8.0*atan(1.0) * (real(fr)) / period_1) ! ADDED BDM
  enddo

  ! Read all relevant parameters and write them to the model log.

  CS%Time => Time ! direct a pointer to the current model time target

  call get_param(param_file, mdl, "INPUTDIR", CS%inputdir, default=".")
                 CS%inputdir = slasher(CS%inputdir)

  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "INTERNAL_TIDE_FREQS", num_freq, &
                 "The number of distinct internal tide frequency bands "//&
                 "that will be calculated.", default=1)
  call get_param(param_file, mdl, "INTERNAL_TIDE_MODES", num_mode, &
                 "The number of distinct internal tide modes "//&
                 "that will be calculated.", default=1)
  call get_param(param_file, mdl, "INTERNAL_TIDE_ANGLES", num_angle, &
                 "The number of angular resolution bands for the internal "//&
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
      "Inconsistent number of modes.")
  if (4*(num_angle/4) /= num_angle) call MOM_error(FATAL, &
    "Internal_tides_init: INTERNAL_TIDE_ANGLES must be a multiple of 4.")

  CS%diag => diag

  call get_param(param_file, mdl, "INTERNAL_TIDE_DECAY_RATE", CS%decay_rate, &
                 "The rate at which internal tide energy is lost to the "//&
                 "interior ocean internal wave field.", &
                 units="s-1", default=0.0, scale=US%T_to_s)
  call get_param(param_file, mdl, "INTERNAL_TIDE_VOLUME_BASED_CFL", CS%vol_CFL, &
                 "If true, use the ratio of the open face lengths to the "//&
                 "tracer cell areas when estimating CFL numbers in the "//&
                 "internal tide code.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_CORNER_ADVECT", CS%corner_adv, &
                 "If true, internal tide ray-tracing advection uses a "//&
                 "corner-advection scheme rather than PPM.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_SIMPLE_2ND_PPM", CS%simple_2nd, &
                 "If true, CONTINUITY_PPM uses a simple 2nd order "//&
                 "(arithmetic mean) interpolation of the edge values. "//&
                 "This may give better PV conservation properties. While "//&
                 "it formally reduces the accuracy of the continuity "//&
                 "solver itself in the strongly advective limit, it does "//&
                 "not reduce the overall order of accuracy of the dynamic "//&
                 "core.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_UPWIND_1ST", CS%upwind_1st, &
                 "If true, the internal tide ray-tracing advection uses "//&
                 "1st-order upwind advection.  This scheme is highly "//&
                 "continuity solver.  This scheme is highly "//&
                 "diffusive but may be useful for debugging.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_BACKGROUND_DRAG", &
                 CS%apply_background_drag, "If true, the internal tide "//&
                 "ray-tracing advection uses a background drag term as a sink.",&
                 default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_QUAD_DRAG", CS%apply_bottom_drag, &
                 "If true, the internal tide ray-tracing advection uses "//&
                 "a quadratic bottom drag term as a sink.", default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_WAVE_DRAG", CS%apply_wave_drag, &
                 "If true, apply scattering due to small-scale roughness as a sink.", &
                 default=.false.)
  call get_param(param_file, mdl, "INTERNAL_TIDE_FROUDE_DRAG", CS%apply_Froude_drag, &
                 "If true, apply wave breaking as a sink.", &
                 default=.false.)
  call get_param(param_file, mdl, "CDRAG", CS%cdrag, &
                 "CDRAG is the drag coefficient relating the magnitude of "//&
                 "the velocity field to the bottom stress.", units="nondim", &
                 default=0.003)
  call get_param(param_file, mdl, "INTERNAL_TIDE_ENERGIZED_ANGLE", CS%energized_angle, &
                 "If positive, only one angular band of the internal tides "//&
                 "gets all of the energy.  (This is for debugging.)", default=-1)
  call get_param(param_file, mdl, "USE_PPM_ANGULAR", CS%use_PPMang, &
                 "If true, use PPM for advection of energy in angular space.", &
                 default=.false.)
  call get_param(param_file, mdl, "GAMMA_ITIDES", CS%q_itides, &
                 "The fraction of the internal tidal energy that is "//&
                 "dissipated locally with INT_TIDE_DISSIPATION. "//&
                 "THIS NAME COULD BE BETTER.", &
                 units="nondim", default=0.3333)
  call get_param(param_file, mdl, "KAPPA_ITIDES", kappa_itides, &
               "A topographic wavenumber used with INT_TIDE_DISSIPATION. "//&
               "The default is 2pi/10 km, as in St.Laurent et al. 2002.", &
               units="m-1", default=8.e-4*atan(1.0), scale=US%L_to_m)
  call get_param(param_file, mdl, "KAPPA_H2_FACTOR", kappa_h2_factor, &
               "A scaling factor for the roughness amplitude with n"//&
               "INT_TIDE_DISSIPATION.",  units="nondim", default=1.0)

  ! Allocate various arrays needed for loss rates
  allocate(h2(isd:ied,jsd:jed)) ; h2(:,:) = 0.0
  allocate(CS%TKE_itidal_loss_fixed(isd:ied,jsd:jed))
    CS%TKE_itidal_loss_fixed = 0.0
  allocate(CS%TKE_leak_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode))
    CS%TKE_leak_loss(:,:,:,:,:) = 0.0
  allocate(CS%TKE_quad_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode))
    CS%TKE_quad_loss(:,:,:,:,:) = 0.0
  allocate(CS%TKE_itidal_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode))
    CS%TKE_itidal_loss(:,:,:,:,:) = 0.0
  allocate(CS%TKE_Froude_loss(isd:ied,jsd:jed,num_angle,num_freq,num_mode))
    CS%TKE_Froude_loss(:,:,:,:,:) = 0.0
  allocate(CS%tot_leak_loss(isd:ied,jsd:jed))   ; CS%tot_leak_loss(:,:) = 0.0
  allocate(CS%tot_quad_loss(isd:ied,jsd:jed) )  ; CS%tot_quad_loss(:,:) = 0.0
  allocate(CS%tot_itidal_loss(isd:ied,jsd:jed)) ; CS%tot_itidal_loss(:,:) = 0.0
  allocate(CS%tot_Froude_loss(isd:ied,jsd:jed)) ; CS%tot_Froude_loss(:,:) = 0.0

  ! Compute the fixed part of the bottom drag loss from baroclinic modes
  call get_param(param_file, mdl, "H2_FILE", h2_file, &
          "The path to the file containing the sub-grid-scale "//&
          "topographic roughness amplitude with INT_TIDE_DISSIPATION.", &
          fail_if_missing=.true.)
  filename = trim(CS%inputdir) // trim(h2_file)
  call log_param(param_file, mdl, "INPUTDIR/H2_FILE", filename)
  call MOM_read_data(filename, 'h2', h2, G%domain, timelevel=1, scale=US%m_to_Z)
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! Restrict rms topo to 10 percent of column depth.
    h2(i,j) = min(0.01*(G%bathyT(i,j))**2, h2(i,j))
    ! Compute the fixed part; units are [R L-2 Z3 ~> kg m-2] here
    ! will be multiplied by N and the squared near-bottom velocity to get into [R Z3 T-3 ~> W m-2]
    CS%TKE_itidal_loss_fixed(i,j) = 0.5*kappa_h2_factor*GV%Rho0 * US%L_to_Z*kappa_itides * h2(i,j)
  enddo ; enddo

  deallocate(h2)

  ! Read in prescribed coast/ridge/shelf angles from file
  call get_param(param_file, mdl, "REFL_ANGLE_FILE", refl_angle_file, &
               "The path to the file containing the local angle of "//&
               "the coastline/ridge/shelf with respect to the equator.", &
               fail_if_missing=.false., default='')
  filename = trim(CS%inputdir) // trim(refl_angle_file)
  allocate(CS%refl_angle(isd:ied,jsd:jed)) ; CS%refl_angle(:,:) = CS%nullangle
  if (file_exists(filename, G%domain)) then
    call log_param(param_file, mdl, "INPUTDIR/REFL_ANGLE_FILE", filename)
    call MOM_read_data(filename, 'refl_angle', CS%refl_angle, &
                       G%domain, timelevel=1)
  else
    if (trim(refl_angle_file) /= '' ) call MOM_error(FATAL, &
                                                     "REFL_ANGLE_FILE: "//trim(filename)//" not found")
  endif
  ! replace NANs with null value
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if (is_NaN(CS%refl_angle(i,j))) CS%refl_angle(i,j) = CS%nullangle
  enddo ; enddo
  call pass_var(CS%refl_angle,G%domain)

  ! Read in prescribed partial reflection coefficients from file
  call get_param(param_file, mdl, "REFL_PREF_FILE", refl_pref_file, &
               "The path to the file containing the reflection coefficients.", &
               fail_if_missing=.false., default='')
  filename = trim(CS%inputdir) // trim(refl_pref_file)
  allocate(CS%refl_pref(isd:ied,jsd:jed)) ; CS%refl_pref(:,:) = 1.0
  if (file_exists(filename, G%domain)) then
    call log_param(param_file, mdl, "INPUTDIR/REFL_PREF_FILE", filename)
    call MOM_read_data(filename, 'refl_pref', CS%refl_pref, G%domain, timelevel=1)
  else
    if (trim(refl_pref_file) /= '' ) call MOM_error(FATAL, &
                                                    "REFL_PREF_FILE: "//trim(filename)//" not found")
  endif
  !CS%refl_pref = CS%refl_pref*1 ! adjust partial reflection if desired
  call pass_var(CS%refl_pref,G%domain)

  ! Tag reflection cells with partial reflection (done here for speed)
  allocate(CS%refl_pref_logical(isd:ied,jsd:jed)) ; CS%refl_pref_logical(:,:) = .false.
  do j=jsd,jed
    do i=isd,ied
      ! flag cells with partial reflection
      if (CS%refl_angle(i,j) /= CS%nullangle .and. &
        CS%refl_pref(i,j) < 1.0 .and. CS%refl_pref(i,j) > 0.0) then
        CS%refl_pref_logical(i,j) = .true.
      endif
    enddo
  enddo

  ! Read in double-reflective (ridge) tags from file
  call get_param(param_file, mdl, "REFL_DBL_FILE", refl_dbl_file, &
               "The path to the file containing the double-reflective ridge tags.", &
               fail_if_missing=.false., default='')
  filename = trim(CS%inputdir) // trim(refl_dbl_file)
  allocate(ridge_temp(isd:ied,jsd:jed)) ; ridge_temp(:,:) = 0.0
  if (file_exists(filename, G%domain)) then
    call log_param(param_file, mdl, "INPUTDIR/REFL_DBL_FILE", filename)
    call MOM_read_data(filename, 'refl_dbl', ridge_temp, G%domain, timelevel=1)
  else
    if (trim(refl_dbl_file) /= '' ) call MOM_error(FATAL, &
                                                   "REFL_DBL_FILE: "//trim(filename)//" not found")
  endif
  call pass_var(ridge_temp,G%domain)
  allocate(CS%refl_dbl(isd:ied,jsd:jed)) ; CS%refl_dbl(:,:) = .false.
  do i=isd,ied; do j=jsd,jed
    if (ridge_temp(i,j) == 1) then; CS%refl_dbl(i,j) = .true.
    else ; CS%refl_dbl(i,j) = .false. ; endif
  enddo ; enddo

  ! Read in prescribed land mask from file (if overwriting -BDM).
  ! This should be done in MOM_initialize_topography subroutine
  ! defined in MOM_fixed_initialization.F90 (BDM)
  !call get_param(param_file, mdl, "LAND_MASK_FILE", land_mask_file, &
  !             "The path to the file containing the land mask.", &
  !             fail_if_missing=.false.)
  !filename = trim(CS%inputdir) // trim(land_mask_file)
  !call log_param(param_file, mdl, "INPUTDIR/LAND_MASK_FILE", filename)
  !G%mask2dCu(:,:) = 1 ; G%mask2dCv(:,:) = 1 ; G%mask2dT(:,:)  = 1
  !call MOM_read_data(filename, 'land_mask', G%mask2dCu, G%domain, timelevel=1)
  !call MOM_read_data(filename, 'land_mask', G%mask2dCv, G%domain, timelevel=1)
  !call MOM_read_data(filename, 'land_mask', G%mask2dT, G%domain, timelevel=1)
  !call pass_vector(G%mask2dCu, G%mask2dCv, G%domain, To_All+Scalar_Pair, CGRID_NE)
  !call pass_var(G%mask2dT,G%domain)

  ! Read in prescribed partial east face blockages from file (if overwriting -BDM)
  !call get_param(param_file, mdl, "dy_Cu_FILE", dy_Cu_file, &
  !             "The path to the file containing the east face blockages.", &
  !             fail_if_missing=.false.)
  !filename = trim(CS%inputdir) // trim(dy_Cu_file)
  !call log_param(param_file, mdl, "INPUTDIR/dy_Cu_FILE", filename)
  !G%dy_Cu(:,:) = 0.0
  !call MOM_read_data(filename, 'dy_Cu', G%dy_Cu, G%domain, timelevel=1, scale=US%m_to_L)

  ! Read in prescribed partial north face blockages from file (if overwriting -BDM)
  !call get_param(param_file, mdl, "dx_Cv_FILE", dx_Cv_file, &
  !             "The path to the file containing the north face blockages.", &
  !             fail_if_missing=.false.)
  !filename = trim(CS%inputdir) // trim(dx_Cv_file)
  !call log_param(param_file, mdl, "INPUTDIR/dx_Cv_FILE", filename)
  !G%dx_Cv(:,:) = 0.0
  !call MOM_read_data(filename, 'dx_Cv', G%dx_Cv, G%domain, timelevel=1, scale=US%m_to_L)
  !call pass_vector(G%dy_Cu, G%dx_Cv, G%domain, To_All+Scalar_Pair, CGRID_NE)

  ! Register maps of reflection parameters
  CS%id_refl_ang = register_diag_field('ocean_model', 'refl_angle', diag%axesT1, &
                 Time, 'Local angle of coastline/ridge/shelf with respect to equator', 'rad')
  CS%id_refl_pref = register_diag_field('ocean_model', 'refl_pref', diag%axesT1, &
                 Time, 'Partial reflection coefficients', '')
  CS%id_dx_Cv = register_diag_field('ocean_model', 'dx_Cv', diag%axesT1, &
                 Time, 'North face unblocked width', 'm', conversion=US%L_to_m)
  CS%id_dy_Cu = register_diag_field('ocean_model', 'dy_Cu', diag%axesT1, &
                 Time, 'East face unblocked width', 'm', conversion=US%L_to_m)
  CS%id_land_mask = register_diag_field('ocean_model', 'land_mask', diag%axesT1, &
                 Time, 'Land mask', 'logical')            ! used if overriding (BDM)
  ! Output reflection parameters as diags here (not needed every timestep)
  if (CS%id_refl_ang > 0)   call post_data(CS%id_refl_ang, CS%refl_angle, CS%diag)
  if (CS%id_refl_pref > 0)  call post_data(CS%id_refl_pref, CS%refl_pref, CS%diag)
  if (CS%id_dx_Cv > 0)      call post_data(CS%id_dx_Cv, G%dx_Cv, CS%diag)
  if (CS%id_dy_Cu > 0)      call post_data(CS%id_dy_Cu, G%dy_Cu, CS%diag)
  if (CS%id_land_mask > 0)  call post_data(CS%id_land_mask, G%mask2dT, CS%diag)

  ! Register 2-D energy density (summed over angles, freq, modes)
  CS%id_tot_En = register_diag_field('ocean_model', 'ITide_tot_En', diag%axesT1, &
                 Time, 'Internal tide total energy density', &
                 'J m-2', conversion=US%RZ3_T3_to_W_m2*US%T_to_s)
  ! Register 2-D drag scale used for quadratic bottom drag
  CS%id_itide_drag = register_diag_field('ocean_model', 'ITide_drag', diag%axesT1, &
                 Time, 'Interior and bottom drag internal tide decay timescale', 's-1', conversion=US%s_to_T)
  !Register 2-D energy input into internal tides
  CS%id_TKE_itidal_input = register_diag_field('ocean_model', 'TKE_itidal_input', diag%axesT1, &
                 Time, 'Conversion from barotropic to baroclinic tide, '//&
                 'a fraction of which goes into rays', &
                 'W m-2', conversion=US%RZ3_T3_to_W_m2)
  ! Register 2-D energy losses (summed over angles, freq, modes)
  CS%id_tot_leak_loss = register_diag_field('ocean_model', 'ITide_tot_leak_loss', diag%axesT1, &
                Time, 'Internal tide energy loss to background drag', &
                'W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_tot_quad_loss = register_diag_field('ocean_model', 'ITide_tot_quad_loss', diag%axesT1, &
                Time, 'Internal tide energy loss to bottom drag', &
                'W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_tot_itidal_loss = register_diag_field('ocean_model', 'ITide_tot_itidal_loss', diag%axesT1, &
                Time, 'Internal tide energy loss to wave drag', &
                'W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_tot_Froude_loss = register_diag_field('ocean_model', 'ITide_tot_Froude_loss', diag%axesT1, &
                Time, 'Internal tide energy loss to wave breaking', &
                'W m-2', conversion=US%RZ3_T3_to_W_m2)
  CS%id_tot_allprocesses_loss = register_diag_field('ocean_model', 'ITide_tot_allprocesses_loss', diag%axesT1, &
                Time, 'Internal tide energy loss summed over all processes', &
                'W m-2', conversion=US%RZ3_T3_to_W_m2)

  allocate(CS%id_En_mode(CS%nFreq,CS%nMode)) ; CS%id_En_mode(:,:) = -1
  allocate(CS%id_En_ang_mode(CS%nFreq,CS%nMode)) ; CS%id_En_ang_mode(:,:) = -1
  allocate(CS%id_itidal_loss_mode(CS%nFreq,CS%nMode)) ; CS%id_itidal_loss_mode(:,:) = -1
  allocate(CS%id_allprocesses_loss_mode(CS%nFreq,CS%nMode)) ; CS%id_allprocesses_loss_mode(:,:) = -1
  allocate(CS%id_itidal_loss_ang_mode(CS%nFreq,CS%nMode)) ; CS%id_itidal_loss_ang_mode(:,:) = -1
  allocate(CS%id_Ub_mode(CS%nFreq,CS%nMode)) ; CS%id_Ub_mode(:,:) = -1
  allocate(CS%id_cp_mode(CS%nFreq,CS%nMode)) ; CS%id_cp_mode(:,:) = -1

  allocate(angles(CS%NAngle)) ; angles(:) = 0.0
  Angle_size = (8.0*atan(1.0)) / (real(num_angle))
  do a=1,num_angle ; angles(a) = (real(a) - 1) * Angle_size ; enddo

  id_ang = diag_axis_init("angle", angles, "Radians", "N", "Angular Orienation of Fluxes")
  call define_axes_group(diag, (/ diag%axesT1%handles(1), diag%axesT1%handles(2), id_ang /), axes_ang)

  do fr=1,CS%nFreq ; write(freq_name(fr), '("freq",i1)') fr ; enddo
  do m=1,CS%nMode ; do fr=1,CS%nFreq
    ! Register 2-D energy density (summed over angles) for each freq and mode
    write(var_name, '("Itide_En_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy density in frequency ",i1," mode ",i1)') fr, m
    CS%id_En_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'J m-2', conversion=US%RZ3_T3_to_W_m2*US%T_to_s)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 3-D (i,j,a) energy density for each freq and mode
    write(var_name, '("Itide_En_ang_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide angular energy density in frequency ",i1," mode ",i1)') fr, m
    CS%id_En_ang_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 axes_ang, Time, var_descript, 'J m-2 band-1', conversion=US%RZ3_T3_to_W_m2*US%T_to_s)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 2-D energy loss (summed over angles) for each freq and mode
    ! wave-drag only
    write(var_name, '("Itide_wavedrag_loss_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy loss due to wave-drag from frequency ",i1," mode ",i1)') fr, m
    CS%id_itidal_loss_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'W m-2', conversion=US%RZ3_T3_to_W_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)
    ! all loss processes
    write(var_name, '("Itide_allprocesses_loss_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy loss due to all processes from frequency ",i1," mode ",i1)') fr, m
    CS%id_allprocesses_loss_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'W m-2', conversion=US%RZ3_T3_to_W_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 3-D (i,j,a) energy loss for each freq and mode
    ! wave-drag only
    write(var_name, '("Itide_wavedrag_loss_ang_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Internal tide energy loss due to wave-drag from frequency ",i1," mode ",i1)') fr, m
    CS%id_itidal_loss_ang_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 axes_ang, Time, var_descript, 'W m-2 band-1', conversion=US%RZ3_T3_to_W_m2)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 2-D period-averaged near-bottom horizonal velocity for each freq and mode
    write(var_name, '("Itide_Ub_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Near-bottom horizonal velocity for frequency ",i1," mode ",i1)') fr, m
    CS%id_Ub_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'm s-1', conversion=US%L_T_to_m_s)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

    ! Register 2-D horizonal phase velocity for each freq and mode
    write(var_name, '("Itide_cp_freq",i1,"_mode",i1)') fr, m
    write(var_descript, '("Horizonal phase velocity for frequency ",i1," mode ",i1)') fr, m
    CS%id_cp_mode(fr,m) = register_diag_field('ocean_model', var_name, &
                 diag%axesT1, Time, var_descript, 'm s-1', conversion=US%L_T_to_m_s)
    call MOM_mesg("Registering "//trim(var_name)//", Described as: "//var_descript, 5)

  enddo ; enddo

  ! Initialize wave_structure (not sure if this should be here - BDM)
  call wave_structure_init(Time, G, param_file, diag, CS%wave_structure_CSp)

end subroutine internal_tides_init

!> This subroutine deallocates the memory associated with the internal tides control structure
subroutine internal_tides_end(CS)
  type(int_tide_CS), pointer :: CS  !< A pointer to the control structure returned by a previous
                                    !! call to internal_tides_init, it will be deallocated here.

  if (associated(CS)) then
    if (associated(CS%En)) deallocate(CS%En)
    if (allocated(CS%frequency)) deallocate(CS%frequency)
    if (allocated(CS%id_En_mode)) deallocate(CS%id_En_mode)
    if (allocated(CS%id_Ub_mode)) deallocate(CS%id_Ub_mode)
    if (allocated(CS%id_cp_mode)) deallocate(CS%id_cp_mode)
    deallocate(CS)
  endif
  CS => NULL()
end subroutine internal_tides_end

end module MOM_internal_tides
