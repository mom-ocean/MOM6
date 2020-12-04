!> The routines here implement the offline tracer algorithm used in MOM6. These are called from step_offline
!! Some routines called here can be found in the MOM_offline_aux module.
module MOM_offline_main

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_ALE,                  only : ALE_CS, ALE_main_offline, ALE_offline_inputs
use MOM_checksums,            only : hchksum, uvchksum
use MOM_cpu_clock,            only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,            only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_cpu_clock,            only : CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diabatic_aux,         only : diabatic_aux_CS, set_pen_shortwave
use MOM_diabatic_driver,      only : diabatic_CS, extract_diabatic_member
use MOM_diabatic_aux,         only : tridiagTS
use MOM_diag_mediator,        only : diag_ctrl, post_data, register_diag_field
use MOM_domains,              only : sum_across_PEs, pass_var, pass_vector
use MOM_error_handler,        only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_error_handler,        only : callTree_enter, callTree_leave
use MOM_file_parser,          only : read_param, get_param, log_version, param_file_type
use MOM_forcing_type,         only : forcing
use MOM_grid,                 only : ocean_grid_type
use MOM_io,                   only : MOM_read_data, MOM_read_vector, CENTER
use MOM_offline_aux,          only : update_offline_from_arrays, update_offline_from_files
use MOM_offline_aux,          only : next_modulo_time, offline_add_diurnal_sw
use MOM_offline_aux,          only : update_h_horizontal_flux, update_h_vertical_flux, limit_mass_flux_3d
use MOM_offline_aux,          only : distribute_residual_uh_barotropic, distribute_residual_vh_barotropic
use MOM_offline_aux,          only : distribute_residual_uh_upwards, distribute_residual_vh_upwards
use MOM_opacity,              only : opacity_CS, optics_type
use MOM_open_boundary,        only : ocean_OBC_type
use MOM_time_manager,         only : time_type, real_to_time
use MOM_tracer_advect,        only : tracer_advect_CS, advect_tracer
use MOM_tracer_diabatic,      only : applyTracerBoundaryFluxesInOut
use MOM_tracer_flow_control,  only : tracer_flow_control_CS, call_tracer_column_fns, call_tracer_stocks
use MOM_tracer_registry,      only : tracer_registry_type, MOM_tracer_chksum, MOM_tracer_chkinv
use MOM_unit_scaling,         only : unit_scale_type
use MOM_variables,            only : thermo_var_ptrs
use MOM_verticalGrid,         only : verticalGrid_type

implicit none ; private

#include "MOM_memory.h"
#include "version_variable.h"

!> The control structure for the offline transport module
type, public :: offline_transport_CS ; private

  ! Pointers to relevant fields from the main MOM control structure
  type(ALE_CS),                  pointer :: ALE_CSp         => NULL()
          !< A pointer to the ALE control structure
  type(diabatic_CS),             pointer :: diabatic_CSp    => NULL()
          !< A pointer to the diabatic control structure
  type(diag_ctrl),               pointer :: diag            => NULL()
          !< Structure that regulates diagnostic output
  type(ocean_OBC_type),          pointer :: OBC             => NULL()
          !< A pointer to the open boundary condition control structure
  type(tracer_advect_CS),        pointer :: tracer_adv_CSp  => NULL()
          !< A pointer to the tracer advection control structure
  type(opacity_CS),              pointer :: opacity_CSp     => NULL()
          !< A pointer to the opacity control structure
  type(tracer_flow_control_CS),  pointer :: tracer_flow_CSp => NULL()
          !< A pointer to control structure that orchestrates the calling of tracer packages
  type(tracer_registry_type),    pointer :: tracer_Reg      => NULL()
          !< A pointer to the tracer registry
  type(thermo_var_ptrs),         pointer :: tv              => NULL()
          !< A structure pointing to various thermodynamic variables
  type(ocean_grid_type),         pointer :: G               => NULL()
          !< Pointer to a structure containing metrics and related information
  type(verticalGrid_type),       pointer :: GV              => NULL()
          !< Pointer to structure containing information about the vertical grid
  type(unit_scale_type),         pointer :: US              => NULL()
          !< structure containing various unit conversion factors
  type(optics_type),             pointer :: optics          => NULL()
          !< Pointer to the optical properties type
  type(diabatic_aux_CS),         pointer :: diabatic_aux_CSp => NULL()
          !< Pointer to the diabatic_aux control structure

  !> Variables related to reading in fields from online run
  integer :: start_index  !< Timelevel to start
  integer :: iter_no      !< Timelevel to start
  integer :: numtime      !< How many timelevels in the input fields
  type(time_type) :: accumulated_time !< Length of time accumulated in the current offline interval
  type(time_type) :: vertical_time !< The next value of accumulate_time at which to apply vertical processes
  ! Index of each of the variables to be read in with separate indices for each variable if they
  ! are set off from each other in time
  integer :: ridx_sum = -1 !< Read index offset of the summed variables
  integer :: ridx_snap = -1 !< Read index offset of the snapshot variables
  integer :: nk_input     !< Number of input levels in the input fields
  character(len=200) :: offlinedir  !< Directory where offline fields are stored
  character(len=200) :: & ! Names of input files
    surf_file,  &         !< Contains surface fields (2d arrays)
    snap_file,  &         !< Snapshotted fields (layer thicknesses)
    sum_file,   &         !< Fields which are accumulated over time
    mean_file             !< Fields averaged over time
  character(len=20)  :: redistribute_method !< 'barotropic' if evenly distributing extra flow
                                            !! throughout entire watercolumn, 'upwards',
                                            !! if trying to do it just in the layers above
                                            !! 'both' if both methods are used
  character(len=20) :: mld_var_name !< Name of the mixed layer depth variable to use
  logical :: fields_are_offset !< True if the time-averaged fields and snapshot fields are
                               !! offset by one time level
  logical :: x_before_y        !< Which horizontal direction is advected first
  logical :: print_adv_offline !< Prints out some updates each advection sub interation
  logical :: skip_diffusion    !< Skips horizontal diffusion of tracers
  logical :: read_sw           !< Read in averaged values for shortwave radiation
  logical :: read_mld          !< Check to see whether mixed layer depths should be read in
  logical :: diurnal_sw        !< Adds a synthetic diurnal cycle on shortwave radiation
  logical :: debug             !< If true, write verbose debugging messages
  logical :: redistribute_barotropic !< Redistributes column-summed residual transports throughout
                                     !! a column weighted by thickness
  logical :: redistribute_upwards    !< Redistributes remaining fluxes only in layers above
                                     !! the current one based as the max allowable transport
                                     !! in that cell
  logical :: read_all_ts_uvh  !< If true, then all timelevels of temperature, salinity, mass transports, and
                              !! Layer thicknesses are read during initialization
  !! Variables controlling some of the numerical considerations of offline transport
  integer :: num_off_iter   !< Number of advection iterations per offline step
  integer :: num_vert_iter  !< Number of vertical iterations per offline step
  integer :: off_ale_mod    !< Sets how frequently the ALE step is done during the advection
  real :: dt_offline        !< Timestep used for offline tracers [T ~> s]
  real :: dt_offline_vertical !< Timestep used for calls to tracer vertical physics [T ~> s]
  real :: evap_CFL_limit    !< Limit on the fraction of the water that can be fluxed out of the top
                            !! layer in a timestep [nondim].  This is Copied from diabatic_CS controlling
                            !! how tracers follow freshwater fluxes
  real :: minimum_forcing_depth !< The smallest depth over which fluxes can be applied [H ~> m or kg m-2].
                            !! This is copied from diabatic_CS controlling how tracers follow freshwater fluxes

  real :: Kd_max        !< Runtime parameter specifying the maximum value of vertical diffusivity
  real :: min_residual  !< The minimum amount of total mass flux before exiting the main advection routine
  !>@{ Diagnostic manager IDs for some fields that may be of interest when doing offline transport
  integer :: &
    id_uhr = -1, &
    id_vhr = -1, &
    id_ear = -1, &
    id_ebr = -1, &
    id_hr = -1,  &
    id_hdiff = -1, &
    id_uhr_redist = -1, &
    id_vhr_redist = -1, &
    id_uhr_end = -1, &
    id_vhr_end = -1, &
    id_eta_pre_distribute  = -1, &
    id_eta_post_distribute = -1, &
    id_h_redist = -1, &
    id_eta_diff_end = -1

  ! Diagnostic IDs for the regridded/remapped input fields
  integer :: &
    id_uhtr_regrid = -1, &
    id_vhtr_regrid = -1, &
    id_temp_regrid = -1, &
    id_salt_regrid = -1, &
    id_h_regrid = -1
  !>@}

  ! IDs for timings of various offline components
  integer :: id_clock_read_fields = -1   !< A CPU time clock
  integer :: id_clock_offline_diabatic = -1  !< A CPU time clock
  integer :: id_clock_offline_adv  = -1  !< A CPU time clock
  integer :: id_clock_redistribute = -1  !< A CPU time clock

  !> Zonal transport that may need to be stored between calls to step_MOM
  real, allocatable, dimension(:,:,:) :: uhtr
  !> Meridional transport that may need to be stored between calls to step_MOM
  real, allocatable, dimension(:,:,:) :: vhtr

  ! Fields at T-point
  real, allocatable, dimension(:,:,:) :: eatr
                   !< Amount of fluid entrained from the layer above within
                   !! one time step [H ~> m or kg m-2]
  real, allocatable, dimension(:,:,:) :: ebtr
                   !< Amount of fluid entrained from the layer below within
                   !! one time step [H ~> m or kg m-2]
  ! Fields at T-points on interfaces
  real, allocatable, dimension(:,:,:) :: Kd     !< Vertical diffusivity
  real, allocatable, dimension(:,:,:) :: h_end  !< Thicknesses at the end of offline timestep

  real, allocatable, dimension(:,:) :: netMassIn  !< Freshwater fluxes into the ocean
  real, allocatable, dimension(:,:) :: netMassOut !< Freshwater fluxes out of the ocean
  real, allocatable, dimension(:,:) :: mld        !< Mixed layer depths at thickness points [Z ~> m].

  ! Allocatable arrays to read in entire fields during initialization
  real, allocatable, dimension(:,:,:,:) :: uhtr_all !< Entire field of zonal transport
  real, allocatable, dimension(:,:,:,:) :: vhtr_all !< Entire field of mericional transport
  real, allocatable, dimension(:,:,:,:) :: hend_all !< Entire field of layer thicknesses
  real, allocatable, dimension(:,:,:,:) :: temp_all !< Entire field of temperatures
  real, allocatable, dimension(:,:,:,:) :: salt_all !< Entire field of salinities

end type offline_transport_CS

public offline_advection_ale
public offline_redistribute_residual
public offline_diabatic_ale
public offline_fw_fluxes_into_ocean
public offline_fw_fluxes_out_ocean
public offline_advection_layer
public register_diags_offline_transport
public update_offline_fields
public insert_offline_main
public extract_offline_main
public post_offline_convergence_diags
public offline_transport_init
public offline_transport_end

contains

!> 3D advection is done by doing flux-limited nonlinear horizontal advection interspersed with an ALE
!! regridding/remapping step. The loop in this routine is exited if remaining residual transports are below
!! a runtime-specified value or a maximum number of iterations is reached.
subroutine offline_advection_ale(fluxes, Time_start, time_interval, CS, id_clock_ale, h_pre, uhtr, vhtr, converged)
  type(forcing),    intent(inout)      :: fluxes        !< pointers to forcing fields
  type(time_type),  intent(in)         :: Time_start    !< starting time of a segment, as a time type
  real,             intent(in)         :: time_interval !< time interval
  type(offline_transport_CS), pointer  :: CS            !< control structure for offline module
  integer,          intent(in)         :: id_clock_ALE  !< Clock for ALE routines
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)), &
                    intent(inout)      :: h_pre         !< layer thicknesses before advection
                                                        !! [H ~> m or kg m-2]
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%GV)), &
                    intent(inout)      :: uhtr          !< Zonal mass transport [H m2 ~> m3 or kg]
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%GV)), &
                    intent(inout)      :: vhtr          !< Meridional mass transport [H m2 ~> m3 or kg]
  logical,          intent(  out)      :: converged     !< True if the iterations have converged

  ! Local pointers
  type(ocean_grid_type),      pointer :: G  => NULL() ! Pointer to a structure containing
                                                      ! metrics and related information
  type(verticalGrid_type),    pointer :: GV => NULL() ! Pointer to structure containing information
                                                      ! about the vertical grid
  ! Work arrays for mass transports
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%GV))   :: uhtr_sub
  ! Meridional mass transports
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%GV))   :: vhtr_sub

  real :: prev_tot_residual, tot_residual  ! Used to keep track of how close to convergence we are

  ! Variables used to keep track of layer thicknesses at various points in the code
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)) :: &
      h_new, &
      h_vol
  ! Fields for eta_diff diagnostic
  real, dimension(SZI_(CS%G),SZJ_(CS%G))         :: eta_pre, eta_end
  integer                                        :: niter, iter
  real                                           :: Inum_iter
  character(len=256) :: mesg  ! The text of an error message
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: isv, iev, jsv, jev ! The valid range of the indices.
  integer :: IsdB, IedB, JsdB, JedB
  logical :: z_first, x_before_y
  real :: evap_CFL_limit  ! Limit on the fraction of the water that can be fluxed out of the
                          ! top layer in a timestep [nondim]
  real :: minimum_forcing_depth ! The smallest depth over which fluxes can be applied [H ~> m or kg m-2]
  real :: dt_iter    ! The timestep to use for each iteration [T ~> s]

  integer :: nstocks
  real :: stock_values(MAX_FIELDS_)
  character(len=20) :: debug_msg
  call cpu_clock_begin(CS%id_clock_offline_adv)

  ! Grid-related pointer assignments
  G => CS%G
  GV => CS%GV

  x_before_y = CS%x_before_y

  ! Initialize some shorthand variables from other structures
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  evap_CFL_limit = CS%evap_CFL_limit
  minimum_forcing_depth = CS%minimum_forcing_depth

  niter = CS%num_off_iter
  Inum_iter = 1./real(niter)
  dt_iter = CS%dt_offline*Inum_iter

  ! Initialize working arrays
  h_new(:,:,:) = 0.0
  h_vol(:,:,:) = 0.0
  uhtr_sub(:,:,:) = 0.0
  vhtr_sub(:,:,:) = 0.0

  ! converged should only be true if there are no remaining mass fluxes
  converged = .false.

  ! Tracers are transported using the stored mass fluxes. Where possible, operators are Strang-split around
  ! the call to
  ! 1)  Using the layer thicknesses and tracer concentrations from the previous timestep,
  !     half of the accumulated vertical mixing (eatr and ebtr) is applied in the call to tracer_column_fns.
  !     For tracers whose source/sink terms need dt, this value is set to 1/2 dt_offline
  ! 2)  Half of the accumulated surface freshwater fluxes are applied
  !! START ITERATION
  ! 3)  Accumulated mass fluxes are used to do horizontal transport. The number of iterations used in
  !     advect_tracer is limited to 2 (e.g x->y->x->y). The remaining mass fluxes are stored for later use
  !     and resulting layer thicknesses fed into the next step
  ! 4)  Tracers and the h-grid are regridded and remapped in a call to ALE. This allows for layers which might
  !     'vanish' because of horizontal mass transport to be 'reinflated'
  ! 5)  Check that transport is done if the remaining mass fluxes equals 0 or if the max number of iterations
  !     has been reached
  !! END ITERATION
  ! 6)  Repeat steps 1 and 2
  ! 7)  Force a remapping to the stored layer thicknesses that correspond to the snapshot of the online model
  ! 8)  Reset T/S and h to their stored snapshotted values to prevent model drift

  ! Copy over the horizontal mass fluxes from the total mass fluxes
  do k=1,nz ; do j=jsd,jed ; do i=isdB,iedB
    uhtr_sub(I,j,k) = uhtr(I,j,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do j=jsdB,jedB ; do i=isd,ied
    vhtr_sub(i,J,k) = vhtr(i,J,k)
  enddo ; enddo ; enddo
  do k=1,nz ; do j=js,je ; do i=is,ie
    h_new(i,j,k) = h_pre(i,j,k)
  enddo ; enddo ; enddo

  if (CS%debug) then
    call hchksum(h_pre,"h_pre before transport",G%HI)
    call uvchksum("[uv]htr_sub before transport", uhtr_sub, vhtr_sub, G%HI)
  endif
  tot_residual = remaining_transport_sum(CS, uhtr, vhtr)
  if (CS%print_adv_offline) then
    write(mesg,'(A,ES24.16)') "Main advection starting transport: ", tot_residual
    call MOM_mesg(mesg)
  endif

  ! This loop does essentially a flux-limited, nonlinear advection scheme until all mass fluxes
  ! are used. ALE is done after the horizontal advection.
  do iter=1,CS%num_off_iter

    do k=1,nz ; do j=js,je ; do i=is,ie
      h_vol(i,j,k) = h_new(i,j,k) * G%US%L_to_m**2*G%areaT(i,j)
      h_pre(i,j,k) = h_new(i,j,k)
    enddo ; enddo ; enddo

    if (CS%debug) then
      call hchksum(h_vol,"h_vol before advect",G%HI)
      call uvchksum("[uv]htr_sub before advect", uhtr_sub, vhtr_sub, G%HI)
      write(debug_msg, '(A,I4.4)') 'Before advect ', iter
      call MOM_tracer_chkinv(debug_msg, G, GV, h_pre, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
    endif

    call advect_tracer(h_pre, uhtr_sub, vhtr_sub, CS%OBC, CS%dt_offline, G, GV, CS%US, &
        CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=1, &
        uhr_out=uhtr, vhr_out=vhtr, h_out=h_new, x_first_in=x_before_y)

    ! Switch the direction every iteration
    x_before_y = .not. x_before_y

    ! Update the new layer thicknesses after one round of advection has happened
    do k=1,nz ; do j=js,je ; do i=is,ie
      h_new(i,j,k) = h_new(i,j,k) / (G%US%L_to_m**2*G%areaT(i,j))
    enddo ; enddo ; enddo

    if (MODULO(iter,CS%off_ale_mod)==0) then
      ! Do ALE remapping/regridding to allow for more advection to occur in the next iteration
      call pass_var(h_new,G%Domain)
      if (CS%debug) then
        call hchksum(h_new,"h_new before ALE",G%HI)
        write(debug_msg, '(A,I4.4)') 'Before ALE ', iter
        call MOM_tracer_chkinv(debug_msg, G, GV, h_new, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
      endif
      call cpu_clock_begin(id_clock_ALE)
      call ALE_main_offline(G, GV, h_new, CS%tv, CS%tracer_Reg, CS%ALE_CSp, CS%OBC, CS%dt_offline)
      call cpu_clock_end(id_clock_ALE)

      if (CS%debug) then
        call hchksum(h_new,"h_new after ALE",G%HI)
        write(debug_msg, '(A,I4.4)') 'After ALE ', iter
        call MOM_tracer_chkinv(debug_msg, G, GV, h_new, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
      endif
    endif

    do k=1,nz ; do j=js,je ; do i=is,ie
      uhtr_sub(I,j,k) = uhtr(I,j,k)
      vhtr_sub(i,J,k) = vhtr(i,J,k)
    enddo ; enddo ; enddo
    call pass_var(h_new, G%Domain)
    call pass_vector(uhtr_sub,vhtr_sub,G%Domain)

    ! Check for whether we've used up all the advection, or if we need to move on because
    ! advection has stalled
    tot_residual = remaining_transport_sum(CS, uhtr, vhtr)
    if (CS%print_adv_offline) then
      write(mesg,'(A,ES24.16)') "Main advection remaining transport: ", tot_residual
      call MOM_mesg(mesg)
    endif
    ! If all the mass transports have been used u, then quit
    if (tot_residual == 0.0) then
      write(mesg,*) "Converged after iteration ", iter
      call MOM_mesg(mesg)
      converged = .true.
      exit
    endif
    ! If advection has stalled or the remaining residual is less than a specified amount, quit
    if ( (tot_residual == prev_tot_residual) .or. (tot_residual<CS%min_residual) ) then
      converged = .false.
      exit
    endif

    prev_tot_residual = tot_residual

  enddo

  ! Make sure that uhtr and vhtr halos are updated
  h_pre(:,:,:) = h_new(:,:,:)
  call pass_vector(uhtr,vhtr,G%Domain)

  if (CS%debug) then
    call hchksum(h_pre,"h after offline_advection_ale",G%HI)
    call uvchksum("[uv]htr after offline_advection_ale", uhtr, vhtr, G%HI)
    call MOM_tracer_chkinv("After offline_advection_ale", G, GV, h_pre, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
  endif

  call cpu_clock_end(CS%id_clock_offline_adv)

end subroutine offline_advection_ale

!> In the case where the main advection routine did not converge, something needs to be done with the remaining
!! transport. Two different ways are offered, 'barotropic' means that the residual is distributed equally
!! throughout the water column. 'upwards' attempts to redistribute the transport in the layers above and will
!! eventually work down the entire water column
subroutine offline_redistribute_residual(CS, h_pre, uhtr, vhtr, converged)
  type(offline_transport_CS), pointer       :: CS    !< control structure from initialize_MOM
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)), &
                              intent(inout) :: h_pre !< layer thicknesses before advection
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%GV)), &
                              intent(inout) :: uhtr  !< Zonal mass transport
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%GV)), &
                              intent(inout) :: vhtr  !< Meridional mass transport
  logical,                    intent(in   ) :: converged !< True if the iterations have converged

  type(ocean_grid_type),      pointer :: G  => NULL() ! Pointer to a structure containing
                                                      ! metrics and related information
  type(verticalGrid_type),    pointer :: GV => NULL() ! Pointer to structure containing information
                                                      ! about the vertical grid
  logical :: x_before_y
  ! Variables used to keep track of layer thicknesses at various points in the code
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)) :: &
      h_new, &
      h_vol

  ! Used to calculate the eta diagnostics
  real, dimension(SZI_(CS%G),SZJ_(CS%G)) :: eta_work
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%GV)) :: uhr  !< Zonal mass transport
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%GV)) :: vhr  !< Meridional mass transport

  character(len=256) :: mesg  ! The text of an error message
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz, iter
  real :: prev_tot_residual, tot_residual, stock_values(MAX_FIELDS_)
  integer :: nstocks

  ! Assign grid pointers
  G  => CS%G
  GV => CS%GV

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  x_before_y = CS%x_before_y

  if (CS%id_eta_pre_distribute>0) then
    eta_work(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      if (h_pre(i,j,k)>GV%Angstrom_H) then
        eta_work(i,j) = eta_work(i,j) + h_pre(i,j,k)
      endif
    enddo ; enddo ; enddo
    call post_data(CS%id_eta_pre_distribute,eta_work,CS%diag)
  endif

  ! These are used to find out how much will be redistributed in this routine
  if (CS%id_h_redist>0) call post_data(CS%id_h_redist, h_pre, CS%diag)
  if (CS%id_uhr_redist>0) call post_data(CS%id_uhr_redist, uhtr, CS%diag)
  if (CS%id_vhr_redist>0) call post_data(CS%id_vhr_redist, vhtr, CS%diag)

  if (converged) return

  if (CS%debug) then
    call MOM_tracer_chkinv("Before redistribute ", G, GV, h_pre, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
  endif

  call cpu_clock_begin(CS%id_clock_redistribute)

  if (CS%redistribute_upwards .or. CS%redistribute_barotropic) then
    do iter = 1, CS%num_off_iter

      ! Perform upwards redistribution
      if (CS%redistribute_upwards) then

        ! Calculate the layer volumes at beginning of redistribute
        do k=1,nz ; do j=js,je ; do i=is,ie
          h_vol(i,j,k) = h_pre(i,j,k)*G%US%L_to_m**2*G%areaT(i,j)
        enddo ; enddo ; enddo
        call pass_var(h_vol,G%Domain)
        call pass_vector(uhtr,vhtr,G%Domain)

        ! Store volumes for advect_tracer
        h_pre(:,:,:) = h_vol(:,:,:)

        if (CS%debug) then
          call MOM_tracer_chksum("Before upwards redistribute ", CS%tracer_Reg%Tr, CS%tracer_Reg%ntr, G)
          call uvchksum("[uv]tr before upwards redistribute", uhtr, vhtr, G%HI)
        endif

        if (x_before_y) then
          call distribute_residual_uh_upwards(G, GV, h_vol, uhtr)
          call distribute_residual_vh_upwards(G, GV, h_vol, vhtr)
        else
          call distribute_residual_vh_upwards(G, GV, h_vol, vhtr)
          call distribute_residual_uh_upwards(G, GV, h_vol, uhtr)
        endif

        call advect_tracer(h_pre, uhtr, vhtr, CS%OBC, CS%dt_offline, G, GV, CS%US, &
            CS%tracer_adv_CSp, CS%tracer_Reg, h_prev_opt = h_pre, max_iter_in=1, &
            h_out=h_new, uhr_out=uhr, vhr_out=vhr, x_first_in=x_before_y)

        if (CS%debug) then
          call MOM_tracer_chksum("After upwards redistribute ", CS%tracer_Reg%Tr, CS%tracer_Reg%ntr, G)
        endif

        ! Convert h_new back to layer thickness for ALE remapping
        do k=1,nz ; do j=js,je ; do i=is,ie
          uhtr(I,j,k) = uhr(I,j,k)
          vhtr(i,J,k) = vhr(i,J,k)
          h_vol(i,j,k) = h_new(i,j,k)
          h_new(i,j,k) = h_new(i,j,k) / (G%US%L_to_m**2*G%areaT(i,j))
          h_pre(i,j,k) = h_new(i,j,k)
        enddo ; enddo ; enddo

      endif ! redistribute upwards

      ! Perform barotropic redistribution
      if (CS%redistribute_barotropic) then

        ! Calculate the layer volumes at beginning of redistribute
        do k=1,nz ; do j=js,je ; do i=is,ie
          h_vol(i,j,k) = h_pre(i,j,k)*G%US%L_to_m**2*G%areaT(i,j)
        enddo ; enddo ; enddo
        call pass_var(h_vol,G%Domain)
        call pass_vector(uhtr,vhtr,G%Domain)

        ! Copy h_vol to h_pre for advect_tracer routine
        h_pre(:,:,:) = h_vol(:,:,:)

        if (CS%debug) then
          call MOM_tracer_chksum("Before barotropic redistribute ", CS%tracer_Reg%Tr, CS%tracer_Reg%ntr, G)
          call uvchksum("[uv]tr before upwards redistribute", uhtr, vhtr, G%HI)
        endif

        if (x_before_y) then
          call distribute_residual_uh_barotropic(G, GV, h_vol, uhtr)
          call distribute_residual_vh_barotropic(G, GV, h_vol, vhtr)
        else
          call distribute_residual_vh_barotropic(G, GV, h_vol, vhtr)
          call distribute_residual_uh_barotropic(G, GV, h_vol, uhtr)
        endif

        call advect_tracer(h_pre, uhtr, vhtr, CS%OBC, CS%dt_offline, G, GV, CS%US, &
            CS%tracer_adv_CSp, CS%tracer_Reg, h_prev_opt = h_pre, max_iter_in=1, &
            h_out=h_new, uhr_out=uhr, vhr_out=vhr, x_first_in=x_before_y)

        if (CS%debug) then
          call MOM_tracer_chksum("After barotropic redistribute ", CS%tracer_Reg%Tr, CS%tracer_Reg%ntr, G)
        endif

        ! Convert h_new back to layer thickness for ALE remapping
        do k=1,nz ; do j=js,je ; do i=is,ie
          uhtr(I,j,k) = uhr(I,j,k)
          vhtr(i,J,k) = vhr(i,J,k)
          h_vol(i,j,k) = h_new(i,j,k)
          h_new(i,j,k) = h_new(i,j,k) / (G%US%L_to_m**2*G%areaT(i,j))
          h_pre(i,j,k) = h_new(i,j,k)
        enddo ; enddo ; enddo

      endif ! redistribute barotropic

      ! Check to see if all transport has been exhausted
      tot_residual = remaining_transport_sum(CS, uhtr, vhtr)
      if (CS%print_adv_offline) then
        write(mesg,'(A,ES24.16)') "Residual advection remaining transport: ", tot_residual
        call MOM_mesg(mesg)
      endif
      ! If the remaining residual is 0, then this return is done
      if (tot_residual==0.0 ) then
        exit
      endif

      if ( (tot_residual == prev_tot_residual) .or. (tot_residual<CS%min_residual) ) exit
      prev_tot_residual = tot_residual

    enddo ! Redistribution iteration
  endif ! If one of the redistribution routines is requested

  if (CS%id_eta_post_distribute>0) then
    eta_work(:,:) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      if (h_pre(i,j,k)>GV%Angstrom_H) then
        eta_work(i,j) = eta_work(i,j) + h_pre(i,j,k)
      endif
    enddo ; enddo ; enddo
    call post_data(CS%id_eta_post_distribute,eta_work,CS%diag)
  endif

  if (CS%id_uhr>0) call post_data(CS%id_uhr,uhtr,CS%diag)
  if (CS%id_vhr>0) call post_data(CS%id_vhr,vhtr,CS%diag)

  if (CS%debug) then
    call hchksum(h_pre,"h_pre after redistribute",G%HI)
    call uvchksum("uhtr after redistribute", uhtr, vhtr, G%HI)
    call MOM_tracer_chkinv("after redistribute ", G, GV, h_new, CS%tracer_Reg%Tr, CS%tracer_Reg%ntr)
  endif

  call cpu_clock_end(CS%id_clock_redistribute)

end subroutine offline_redistribute_residual

!> Sums any non-negligible remaining transport to check for advection convergence
real function remaining_transport_sum(CS, uhtr, vhtr)
  type(offline_transport_CS), pointer  :: CS !< control structure for offline module
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%GV)), intent(in   )  :: uhtr  !< Zonal mass transport
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%GV)), intent(in   )  :: vhtr  !< Meridional mass transport

  ! Local variables
  integer :: i, j, k
  integer :: is, ie, js, je, nz
  real :: h_min !< A layer thickness below roundoff from GV type
  real :: uh_neglect !< A small value of zonal transport that effectively is below roundoff error
  real :: vh_neglect !< A small value of meridional transport that effectively is below roundoff error

  nz = CS%GV%ke
  is = CS%G%isc ; ie = CS%G%iec ; js = CS%G%jsc ; je = CS%G%jec

  h_min = CS%GV%H_subroundoff

  remaining_transport_sum = 0.
  do k=1,nz ; do j=js,je ; do i=is,ie
    uh_neglect = h_min*CS%G%US%L_to_m**2*MIN(CS%G%areaT(i,j),CS%G%areaT(i+1,j))
    vh_neglect = h_min*CS%G%US%L_to_m**2*MIN(CS%G%areaT(i,j),CS%G%areaT(i,j+1))
    if (ABS(uhtr(I,j,k))>uh_neglect) then
      remaining_transport_sum = remaining_transport_sum + ABS(uhtr(I,j,k))
    endif
    if (ABS(vhtr(i,J,k))>vh_neglect) then
      remaining_transport_sum = remaining_transport_sum + ABS(vhtr(i,J,k))
    endif
  enddo ; enddo ; enddo
  call sum_across_PEs(remaining_transport_sum)

end function remaining_transport_sum

!> The vertical/diabatic driver for offline tracers. First the eatr/ebtr associated with the interpolated
!! vertical diffusivities are calculated and then any tracer column functions are done which can include
!! vertical diffuvities and source/sink terms.
subroutine offline_diabatic_ale(fluxes, Time_start, Time_end, CS, h_pre, eatr, ebtr)

  type(forcing),    intent(inout)      :: fluxes     !< pointers to forcing fields
  type(time_type),  intent(in)         :: Time_start !< starting time of a segment, as a time type
  type(time_type),  intent(in)         :: Time_end   !< time interval
  type(offline_transport_CS), pointer  :: CS         !< control structure from initialize_MOM
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)), &
                    intent(inout)      :: h_pre      !< layer thicknesses before advection [H ~> m or kg m-2]
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)), &
                    intent(inout)      :: eatr       !< Entrainment from layer above [H ~> m or kg m-2]
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)), &
                    intent(inout)      :: ebtr       !< Entrainment from layer below [H ~> m or kg m-2]

  real, dimension(SZI_(CS%G),SZJ_(CS%G)) :: &
    sw, sw_vis, sw_nir !< Save old values of shortwave radiation [Q R Z T-1 ~> W m-2]
  real :: hval
  integer :: i,j,k
  integer :: is, ie, js, je, nz
  integer :: k_nonzero
  real :: stock_values(MAX_FIELDS_)
  real :: Kd_bot
  integer :: nstocks
  nz = CS%GV%ke
  is = CS%G%isc ; ie = CS%G%iec ; js = CS%G%jsc ; je = CS%G%jec

  call cpu_clock_begin(CS%id_clock_offline_diabatic)

  call MOM_mesg("Applying tracer source, sinks, and vertical mixing")

  if (CS%debug) then
    call hchksum(h_pre,"h_pre before offline_diabatic_ale",CS%G%HI)
    call hchksum(eatr,"eatr before offline_diabatic_ale",CS%G%HI)
    call hchksum(ebtr,"ebtr before offline_diabatic_ale",CS%G%HI)
    call MOM_tracer_chkinv("Before offline_diabatic_ale", CS%G, CS%GV, h_pre, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
  endif

  eatr(:,:,:) = 0.
  ebtr(:,:,:) = 0.
  ! Calculate eatr and ebtr if vertical diffusivity is read
  ! Because the saved remapped diagnostics from the online run assume a zero minimum thickness
  ! but ALE may have a minimum thickness. Flood the diffusivities for all layers with the value
  ! of Kd closest to the bottom which is non-zero
  do j=js,je ; do i=is,ie
    k_nonzero = nz+1
    ! Find the nonzero bottom Kd
    do k=nz+1,1,-1
      if (CS%Kd(i,j,k)>0.) then
        Kd_bot = CS%Kd(i,j,k)
        k_nonzero = k
        exit
      endif
    enddo
    ! Flood the bottom interfaces
    do k=k_nonzero,nz+1
      CS%Kd(i,j,k) = Kd_bot
    enddo
  enddo ; enddo

  do j=js,je ; do i=is,ie
    eatr(i,j,1) = 0.
  enddo ; enddo
  do k=2,nz ; do j=js,je ; do i=is,ie
    hval=1.0/(CS%GV%H_subroundoff + 0.5*(h_pre(i,j,k-1) + h_pre(i,j,k)))
    eatr(i,j,k) = (CS%GV%m_to_H**2*CS%US%T_to_s) * CS%dt_offline_vertical * hval * CS%Kd(i,j,k)
    ebtr(i,j,k-1) = eatr(i,j,k)
  enddo ; enddo ; enddo
  do j=js,je ; do i=is,ie
    ebtr(i,j,nz) = 0.
  enddo ; enddo

  ! Add diurnal cycle for shortwave radiation (only used if run in ocean-only mode)
  if (CS%diurnal_SW .and. CS%read_sw) then
    sw(:,:) = fluxes%sw(:,:)
    sw_vis(:,:) = fluxes%sw_vis_dir(:,:)
    sw_nir(:,:) = fluxes%sw_nir_dir(:,:)
    call offline_add_diurnal_SW(fluxes, CS%G, Time_start, Time_end)
  endif

  if (associated(CS%optics)) &
    call set_pen_shortwave(CS%optics, fluxes, CS%G, CS%GV, CS%US, CS%diabatic_aux_CSp, &
                           CS%opacity_CSp, CS%tracer_flow_CSp)

  ! Note that tracerBoundaryFluxesInOut within this subroutine should NOT be called
  ! as the freshwater fluxes have already been accounted for
  call call_tracer_column_fns(h_pre, h_pre, eatr, ebtr, fluxes, CS%MLD, CS%dt_offline_vertical, &
                              CS%G, CS%GV, CS%US, CS%tv, CS%optics, CS%tracer_flow_CSp, CS%debug)

  if (CS%diurnal_SW .and. CS%read_sw) then
    fluxes%sw(:,:) = sw(:,:)
    fluxes%sw_vis_dir(:,:) = sw_vis(:,:)
    fluxes%sw_nir_dir(:,:) = sw_nir(:,:)
  endif

  if (CS%debug) then
    call hchksum(h_pre,"h_pre after offline_diabatic_ale",CS%G%HI)
    call hchksum(eatr,"eatr after offline_diabatic_ale",CS%G%HI)
    call hchksum(ebtr,"ebtr after offline_diabatic_ale",CS%G%HI)
    call MOM_tracer_chkinv("After offline_diabatic_ale", CS%G, CS%GV, h_pre, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
  endif

  call cpu_clock_end(CS%id_clock_offline_diabatic)

end subroutine offline_diabatic_ale

!> Apply positive freshwater fluxes (into the ocean) and update netMassOut with only the negative
!! (out of the ocean) fluxes
subroutine offline_fw_fluxes_into_ocean(G, GV, CS, fluxes, h, in_flux_optional)
  type(offline_transport_CS), intent(inout) :: CS !< Offline control structure
  type(ocean_grid_type),      intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),    intent(in)    :: GV !< ocean vertical grid structure
  type(forcing),              intent(inout) :: fluxes !< Surface fluxes container
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                    optional, intent(in)    :: in_flux_optional !< The total time-integrated amount
                                                  !! of tracer that leaves with freshwater

  integer :: i, j, m
  real, dimension(SZI_(G),SZJ_(G)) :: negative_fw !< store all negative fluxes
  logical :: update_h !< Flag for whether h should be updated

  if ( present(in_flux_optional) ) &
    call MOM_error(WARNING, "Positive freshwater fluxes with non-zero tracer concentration not supported yet")

  ! Set all fluxes to 0
  negative_fw(:,:) = 0.

  ! Sort fluxes into positive and negative
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if (fluxes%netMassOut(i,j)<0.0) then
      negative_fw(i,j) = fluxes%netMassOut(i,j)
      fluxes%netMassOut(i,j) = 0.
    endif
  enddo ; enddo

  if (CS%debug) then
    call hchksum(h, "h before fluxes into ocean", G%HI)
    call MOM_tracer_chkinv("Before fluxes into ocean", G, GV, h, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
  endif
  do m = 1,CS%tracer_reg%ntr
    ! Layer thicknesses should only be updated after the last tracer is finished
    update_h = ( m == CS%tracer_reg%ntr )
    call applyTracerBoundaryFluxesInOut(G, GV, CS%tracer_reg%tr(m)%t, CS%dt_offline, fluxes, h, &
                                        CS%evap_CFL_limit, CS%minimum_forcing_depth, update_h_opt = update_h)
  enddo
  if (CS%debug) then
    call hchksum(h, "h after fluxes into ocean", G%HI)
    call MOM_tracer_chkinv("After fluxes into ocean", G, GV, h, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
  endif

  ! Now that fluxes into the ocean are done, save the negative fluxes for later
  fluxes%netMassOut(:,:) = negative_fw(:,:)

end subroutine offline_fw_fluxes_into_ocean

!> Apply negative freshwater fluxes (out of the ocean)
subroutine offline_fw_fluxes_out_ocean(G, GV, CS, fluxes, h, out_flux_optional)
  type(offline_transport_CS), intent(inout) :: CS !< Offline control structure
  type(ocean_grid_type),      intent(in)    :: G  !< Grid structure
  type(verticalGrid_type),    intent(in)    :: GV !< ocean vertical grid structure
  type(forcing),              intent(inout) :: fluxes !< Surface fluxes container
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h  !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), &
                    optional, intent(in)    :: out_flux_optional !< The total time-integrated amount
                                                  !! of tracer that leaves with freshwater

  integer :: m
  logical :: update_h !< Flag for whether h should be updated

  if ( present(out_flux_optional) ) &
    call MOM_error(WARNING, "Negative freshwater fluxes with non-zero tracer concentration not supported yet")

  if (CS%debug) then
    call hchksum(h,"h before fluxes out of ocean",G%HI)
    call MOM_tracer_chkinv("Before fluxes out of ocean", G, GV, h, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
  endif
  do m = 1, CS%tracer_reg%ntr
    ! Layer thicknesses should only be updated after the last tracer is finished
    update_h = ( m == CS%tracer_reg%ntr )
    call applyTracerBoundaryFluxesInOut(G, GV, CS%tracer_reg%tr(m)%t, CS%dt_offline, fluxes, h, &
                                        CS%evap_CFL_limit, CS%minimum_forcing_depth, update_h_opt = update_h)
  enddo
  if (CS%debug) then
    call hchksum(h,"h after fluxes out of ocean",G%HI)
    call MOM_tracer_chkinv("Before fluxes out of ocean", G, GV, h, CS%tracer_reg%Tr, CS%tracer_reg%ntr)
  endif

end subroutine offline_fw_fluxes_out_ocean

!> When in layer mode, 3D horizontal advection using stored mass fluxes must be used. Horizontal advection is
!! done via tracer_advect, whereas the vertical component is actually handled by vertdiff in tracer_column_fns
subroutine offline_advection_layer(fluxes, Time_start, time_interval, CS, h_pre, eatr, ebtr, uhtr, vhtr)
  type(forcing),              intent(inout)    :: fluxes        !< pointers to forcing fields
  type(time_type),            intent(in)       :: Time_start    !< starting time of a segment, as a time type
  real,                       intent(in)       :: time_interval !< Offline transport time interval
  type(offline_transport_CS), pointer          :: CS            !< Control structure for offline module
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)),  intent(inout) :: h_pre !< layer thicknesses before advection
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)),  intent(inout) :: eatr !< Entrainment from layer above
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)),  intent(inout) :: ebtr !< Entrainment from layer below
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%GV)), intent(inout) :: uhtr  !< Zonal mass transport
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%GV)), intent(inout) :: vhtr  !< Meridional mass transport
  ! Local pointers
  type(ocean_grid_type),      pointer :: G  => NULL() ! Pointer to a structure containing
                                                      ! metrics and related information
  type(verticalGrid_type),    pointer :: GV => NULL() ! Pointer to structure containing information
                                                      ! about the vertical grid
  ! Remaining zonal mass transports
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%GV))   :: uhtr_sub
  ! Remaining meridional mass transports
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%GV))   :: vhtr_sub

  real :: sum_abs_fluxes, sum_u, sum_v  ! Used to keep track of how close to convergence we are
  real :: dt_offline

  ! Local variables
  ! Vertical diffusion related variables
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)) :: &
      eatr_sub, &
      ebtr_sub
  ! Variables used to keep track of layer thicknesses at various points in the code
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)) :: &
      h_new, &
      h_vol
  ! Work arrays for temperature and salinity
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)) :: &
      temp_old, salt_old, &
      temp_mean, salt_mean, &
      zero_3dh     !
  integer :: niter, iter
  real    :: Inum_iter
  real    :: dt_iter  ! The timestep of each iteration [T ~> s]
  logical :: converged
  character(len=160) :: mesg  ! The text of an error message
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: isv, iev, jsv, jev ! The valid range of the indices.
  integer :: IsdB, IedB, JsdB, JedB
  logical :: z_first, x_before_y

  G => CS%G ; GV => CS%GV
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  dt_iter = CS%US%s_to_T * time_interval / real(max(1, CS%num_off_iter))

  do iter=1,CS%num_off_iter

    do k = 1, nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      eatr_sub(i,j,k) = eatr(i,j,k)
      ebtr_sub(i,j,k) = ebtr(i,j,k)
    enddo ; enddo ; enddo

    do k = 1, nz ; do j=js-1,je+1 ; do i=is-2,ie+1
      uhtr_sub(I,j,k) = uhtr(I,j,k)
    enddo ; enddo ; enddo

    do k = 1, nz ; do j=js-2,je+1 ; do i=is-1,ie+1
      vhtr_sub(i,J,k) = vhtr(i,J,k)
    enddo ; enddo ; enddo


    ! Calculate 3d mass transports to be used in this iteration
    call limit_mass_flux_3d(G, GV, uhtr_sub, vhtr_sub, eatr_sub, ebtr_sub, h_pre)

    if (z_first) then
      ! First do vertical advection
      call update_h_vertical_flux(G, GV, eatr_sub, ebtr_sub, h_pre, h_new)
      call call_tracer_column_fns(h_pre, h_new, eatr_sub, ebtr_sub, &
          fluxes, CS%mld, dt_iter, G, GV, CS%US, CS%tv, CS%optics, CS%tracer_flow_CSp, CS%debug)
      ! We are now done with the vertical mass transports, so now h_new is h_sub
      do k = 1, nz ; do j=js-1,je+1 ; do i=is-1,ie+1
        h_pre(i,j,k) = h_new(i,j,k)
      enddo ; enddo ; enddo
      call pass_var(h_pre,G%Domain)

      ! Second zonal and meridional advection
      call update_h_horizontal_flux(G, GV, uhtr_sub, vhtr_sub, h_pre, h_new)
      do k = 1, nz ; do i = is-1, ie+1 ; do j=js-1, je+1
        h_vol(i,j,k) = h_pre(i,j,k)*G%US%L_to_m**2*G%areaT(i,j)
      enddo ; enddo ; enddo
      call advect_tracer(h_pre, uhtr_sub, vhtr_sub, CS%OBC, dt_iter, G, GV, CS%US, &
          CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=30, x_first_in=x_before_y)

      ! Done with horizontal so now h_pre should be h_new
      do k = 1, nz ; do i=is-1,ie+1 ; do j=js-1,je+1
          h_pre(i,j,k) = h_new(i,j,k)
      enddo ; enddo ; enddo

    endif

    if (.not. z_first) then

      ! First zonal and meridional advection
      call update_h_horizontal_flux(G, GV, uhtr_sub, vhtr_sub, h_pre, h_new)
      do k = 1, nz ; do i = is-1, ie+1 ; do j=js-1, je+1
        h_vol(i,j,k) = h_pre(i,j,k)*G%US%L_to_m**2*G%areaT(i,j)
      enddo ; enddo ; enddo
      call advect_tracer(h_pre, uhtr_sub, vhtr_sub, CS%OBC, dt_iter, G, GV, CS%US, &
          CS%tracer_adv_CSp, CS%tracer_Reg, h_vol, max_iter_in=30, x_first_in=x_before_y)

      ! Done with horizontal so now h_pre should be h_new
      do k = 1, nz ; do i=is-1,ie+1 ; do j=js-1,je+1
          h_pre(i,j,k) = h_new(i,j,k)
      enddo ; enddo ; enddo

      ! Second vertical advection
      call update_h_vertical_flux(G, GV, eatr_sub, ebtr_sub, h_pre, h_new)
      call call_tracer_column_fns(h_pre, h_new, eatr_sub, ebtr_sub, &
          fluxes, CS%mld, dt_iter, G, GV, CS%US, CS%tv, CS%optics, CS%tracer_flow_CSp, CS%debug)
      ! We are now done with the vertical mass transports, so now h_new is h_sub
      do k = 1, nz ; do i=is-1,ie+1 ; do j=js-1,je+1
        h_pre(i,j,k) = h_new(i,j,k)
      enddo ; enddo ; enddo

    endif

    ! Update remaining transports
    do k = 1, nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      eatr(i,j,k) = eatr(i,j,k) - eatr_sub(i,j,k)
      ebtr(i,j,k) = ebtr(i,j,k) - ebtr_sub(i,j,k)
    enddo ; enddo ; enddo

    do k = 1, nz ; do j=js-1,je+1 ; do i=is-2,ie+1
      uhtr(I,j,k) = uhtr(I,j,k) - uhtr_sub(I,j,k)
    enddo ; enddo ; enddo

    do k = 1, nz ; do j=js-2,je+1 ; do i=is-1,ie+1
      vhtr(i,J,k) = vhtr(i,J,k) - vhtr_sub(i,J,k)
    enddo ; enddo ; enddo

    call pass_var(eatr,G%Domain)
    call pass_var(ebtr,G%Domain)
    call pass_var(h_pre,G%Domain)
    call pass_vector(uhtr,vhtr,G%Domain)
  !
    ! Calculate how close we are to converging by summing the remaining fluxes at each point
    sum_abs_fluxes = 0.0
    sum_u = 0.0
    sum_v = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      sum_u = sum_u + abs(uhtr(I-1,j,k))+abs(uhtr(I,j,k))
      sum_v = sum_v + abs(vhtr(i,J-1,k))+abs(vhtr(I,J,k))
      sum_abs_fluxes = sum_abs_fluxes + abs(eatr(i,j,k)) + abs(ebtr(i,j,k)) + abs(uhtr(I-1,j,k)) + &
          abs(uhtr(I,j,k)) + abs(vhtr(i,J-1,k)) + abs(vhtr(i,J,k))
    enddo ; enddo ; enddo
    call sum_across_PEs(sum_abs_fluxes)

    write(mesg,*) "offline_advection_layer: Remaining u-flux, v-flux:", sum_u, sum_v
    call MOM_mesg(mesg)
    if (sum_abs_fluxes==0) then
      write(mesg,*) 'offline_advection_layer: Converged after iteration', iter
      call MOM_mesg(mesg)
      exit
    endif

    ! Switch order of Strang split every iteration
    z_first = .not. z_first
    x_before_y = .not. x_before_y
  enddo

end subroutine offline_advection_layer

!> Update fields used in this round of offline transport. First fields are updated from files or from arrays
!! read during initialization. Then if in an ALE-dependent coordinate, regrid/remap fields.
subroutine update_offline_fields(CS, h, fluxes, do_ale)
  type(offline_transport_CS), pointer               :: CS !< Control structure for offline module
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)) :: h !< The regridded layer thicknesses
  type(forcing),        intent(inout) :: fluxes !< Pointers to forcing fields
  logical,              intent(in   ) :: do_ale !< True if using ALE
  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)) :: h_start
  is = CS%G%isc ; ie = CS%G%iec ; js = CS%G%jsc ; je = CS%G%jec ; nz = CS%GV%ke

  call cpu_clock_begin(CS%id_clock_read_fields)
  call callTree_enter("update_offline_fields, MOM_offline_main.F90")

  ! Store a copy of the layer thicknesses before ALE regrid/remap
  h_start(:,:,:) = h(:,:,:)

  ! Most fields will be read in from files
  call update_offline_from_files( CS%G, CS%GV, CS%nk_input, CS%mean_file, CS%sum_file, CS%snap_file, CS%surf_file,    &
                                  CS%h_end, CS%uhtr, CS%vhtr, CS%tv%T, CS%tv%S, CS%mld, CS%Kd, fluxes,                &
                                  CS%ridx_sum, CS%ridx_snap, CS%read_mld, CS%read_sw, .not. CS%read_all_ts_uvh, do_ale)
  ! If uh, vh, h_end, temp, salt were read in at the beginning, fields are copied from those arrays
  if (CS%read_all_ts_uvh) then
      call update_offline_from_arrays(CS%G, CS%GV, CS%nk_input, CS%ridx_sum, CS%mean_file, CS%sum_file, CS%snap_file, &
        CS%uhtr, CS%vhtr, CS%h_end, CS%uhtr_all, CS%vhtr_all, CS%hend_all, CS%tv%T, CS%tv%S, CS%temp_all, CS%salt_all)
    endif
  if (CS%debug) then
    call uvchksum("[uv]h after update offline from files and arrays", CS%uhtr, CS%vhtr, CS%G%HI)
  endif

  ! If using an ALE-dependent vertical coordinate, fields will need to be remapped
  if (do_ale) then
    ! These halo passes are necessary because u, v fields will need information 1 step into the halo
    call pass_var(h, CS%G%Domain)
    call pass_var(CS%tv%T, CS%G%Domain)
    call pass_var(CS%tv%S, CS%G%Domain)
    call ALE_offline_inputs(CS%ALE_CSp, CS%G, CS%GV, h, CS%tv, CS%tracer_Reg, CS%uhtr, CS%vhtr, CS%Kd, &
                            CS%debug, CS%OBC)
    if (CS%id_temp_regrid>0) call post_data(CS%id_temp_regrid, CS%tv%T, CS%diag)
    if (CS%id_salt_regrid>0) call post_data(CS%id_salt_regrid, CS%tv%S, CS%diag)
    if (CS%id_uhtr_regrid>0) call post_data(CS%id_uhtr_regrid, CS%uhtr, CS%diag)
    if (CS%id_vhtr_regrid>0) call post_data(CS%id_vhtr_regrid, CS%vhtr, CS%diag)
    if (CS%id_h_regrid>0) call post_data(CS%id_h_regrid, h, CS%diag)
    if (CS%debug) then
      call uvchksum("[uv]h after ALE regridding/remapping of inputs", CS%uhtr, CS%vhtr, CS%G%HI)
      call hchksum(h_start,"h_start after update offline from files and arrays", CS%G%HI)
    endif
  endif

  ! Update halos for some
  call pass_var(CS%h_end, CS%G%Domain)
  call pass_var(CS%tv%T, CS%G%Domain)
  call pass_var(CS%tv%S, CS%G%Domain)

  ! Update the read indices
  CS%ridx_snap = next_modulo_time(CS%ridx_snap,CS%numtime)
  CS%ridx_sum = next_modulo_time(CS%ridx_sum,CS%numtime)

  ! Apply masks/factors at T, U, and V points
  do k=1,nz ; do j=js,je ; do i=is,ie
    if (CS%G%mask2dT(i,j)<1.0) then
      CS%h_end(i,j,k) = CS%GV%Angstrom_H
    endif
  enddo ; enddo ; enddo

  do k=1,nz+1 ; do j=js,je ; do i=is,ie
    CS%Kd(i,j,k) = max(0.0, CS%Kd(i,j,k))
    if (CS%Kd_max>0.) then
      CS%Kd(i,j,k) = MIN(CS%Kd_max, CS%Kd(i,j,k))
    endif
  enddo ; enddo ; enddo

  do k=1,nz ; do J=js-1,je ; do i=is,ie
    if (CS%G%mask2dCv(i,J)<1.0) then
      CS%vhtr(i,J,k) = 0.0
    endif
  enddo ; enddo ; enddo

  do k=1,nz ; do j=js,je ; do I=is-1,ie
    if (CS%G%mask2dCu(I,j)<1.0) then
      CS%uhtr(I,j,k) = 0.0
    endif
  enddo ; enddo ; enddo

  if (CS%debug) then
    call uvchksum("[uv]htr_sub after update_offline_fields", CS%uhtr, CS%vhtr, CS%G%HI)
    call hchksum(CS%h_end, "h_end after update_offline_fields", CS%G%HI)
    call hchksum(CS%tv%T, "Temp after update_offline_fields", CS%G%HI)
    call hchksum(CS%tv%S, "Salt after update_offline_fields", CS%G%HI)
  endif

  call callTree_leave("update_offline_fields")
  call cpu_clock_end(CS%id_clock_read_fields)

end subroutine update_offline_fields

!> Initialize additional diagnostics required for offline tracer transport
subroutine register_diags_offline_transport(Time, diag, CS)

  type(offline_transport_CS), pointer :: CS   !< Control structure for offline module
  type(time_type),         intent(in) :: Time !< current model time
  type(diag_ctrl),         intent(in) :: diag !< Structure that regulates diagnostic output

  ! U-cell fields
  CS%id_uhr = register_diag_field('ocean_model', 'uhr', diag%axesCuL, Time, &
    'Zonal thickness fluxes remaining at end of advection', 'kg')
  CS%id_uhr_redist = register_diag_field('ocean_model', 'uhr_redist', diag%axesCuL, Time, &
    'Zonal thickness fluxes to be redistributed vertically', 'kg')
  CS%id_uhr_end = register_diag_field('ocean_model', 'uhr_end', diag%axesCuL, Time, &
    'Zonal thickness fluxes at end of offline step', 'kg')

  ! V-cell fields
  CS%id_vhr = register_diag_field('ocean_model', 'vhr', diag%axesCvL, Time, &
    'Meridional thickness fluxes remaining at end of advection', 'kg')
  CS%id_vhr_redist = register_diag_field('ocean_model', 'vhr_redist', diag%axesCvL, Time, &
    'Meridional thickness to be redistributed vertically', 'kg')
  CS%id_vhr_end = register_diag_field('ocean_model', 'vhr_end', diag%axesCvL, Time, &
    'Meridional thickness at end of offline step', 'kg')

  ! T-cell fields
  CS%id_hdiff  = register_diag_field('ocean_model', 'hdiff', diag%axesTL, Time, &
    'Difference between the stored and calculated layer thickness', 'm')
  CS%id_hr  = register_diag_field('ocean_model', 'hr', diag%axesTL, Time, &
    'Layer thickness at end of offline step', 'm')
  CS%id_ear  = register_diag_field('ocean_model', 'ear', diag%axesTL, Time, &
    'Remaining thickness entrained from above', 'm')
  CS%id_ebr  = register_diag_field('ocean_model', 'ebr', diag%axesTL, Time, &
    'Remaining thickness entrained from below', 'm')
  CS%id_eta_pre_distribute = register_diag_field('ocean_model','eta_pre_distribute', &
    diag%axesT1, Time, 'Total water column height before residual transport redistribution','m')
  CS%id_eta_post_distribute = register_diag_field('ocean_model','eta_post_distribute', &
    diag%axesT1, Time, 'Total water column height after residual transport redistribution','m')
  CS%id_eta_diff_end = register_diag_field('ocean_model','eta_diff_end', diag%axesT1, Time, &
    'Difference in total water column height from online and offline ' // &
    'at the end of the offline timestep','m')
  CS%id_h_redist = register_diag_field('ocean_model','h_redist', diag%axesTL, Time, &
    'Layer thicknesses before redistribution of mass fluxes','m')

  ! Regridded/remapped input fields
  CS%id_uhtr_regrid = register_diag_field('ocean_model', 'uhtr_regrid', diag%axesCuL, Time, &
                                          'Zonal mass transport regridded/remapped onto offline grid','kg')
  CS%id_vhtr_regrid = register_diag_field('ocean_model', 'vhtr_regrid', diag%axesCvL, Time, &
                                          'Meridional mass transport regridded/remapped onto offline grid','kg')
  CS%id_temp_regrid = register_diag_field('ocean_model', 'temp_regrid', diag%axesTL, Time, &
                                          'Temperature regridded/remapped onto offline grid','C')
  CS%id_salt_regrid = register_diag_field('ocean_model', 'salt_regrid', diag%axesTL, Time, &
                                          'Salinity regridded/remapped onto offline grid','g kg-1')
  CS%id_h_regrid = register_diag_field('ocean_model', 'h_regrid', diag%axesTL, Time, &
                                          'Layer thicknesses regridded/remapped onto offline grid','m')


end subroutine register_diags_offline_transport

!> Posts diagnostics related to offline convergence diagnostics
subroutine post_offline_convergence_diags(CS, h_off, h_end, uhtr, vhtr)
  type(offline_transport_CS), intent(in   ) :: CS     !< Offline control structure
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)),  intent(inout) :: h_off  !< Thicknesses at end of offline step
  real, dimension(SZI_(CS%G),SZJ_(CS%G),SZK_(CS%GV)),  intent(inout) :: h_end  !< Stored thicknesses
  real, dimension(SZIB_(CS%G),SZJ_(CS%G),SZK_(CS%GV)), intent(inout) :: uhtr   !< Remaining zonal mass transport
  real, dimension(SZI_(CS%G),SZJB_(CS%G),SZK_(CS%GV)), intent(inout) :: vhtr   !< Remaining meridional mass transport

  real, dimension(SZI_(CS%G),SZJ_(CS%G)) :: eta_diff
  integer :: i, j, k

  if (CS%id_eta_diff_end>0) then
    ! Calculate difference in column thickness
    eta_diff = 0.
    do k=1,CS%GV%ke ; do j=CS%G%jsc,CS%G%jec ; do i=CS%G%isc,CS%G%iec
      eta_diff(i,j) = eta_diff(i,j) + h_off(i,j,k)
    enddo ; enddo ; enddo
    do k=1,CS%GV%ke ; do j=CS%G%jsc,CS%G%jec ; do i=CS%G%isc,CS%G%iec
      eta_diff(i,j) = eta_diff(i,j) - h_end(i,j,k)
    enddo ; enddo ; enddo

    call post_data(CS%id_eta_diff_end, eta_diff, CS%diag)
  endif

  if (CS%id_hdiff>0) call post_data(CS%id_hdiff, h_off-h_end, CS%diag)
  if (CS%id_hr>0) call post_data(CS%id_hr, h_off, CS%diag)
  if (CS%id_uhr_end>0) call post_data(CS%id_uhr_end, uhtr, CS%diag)
  if (CS%id_vhr_end>0) call post_data(CS%id_vhr_end, vhtr, CS%diag)

end subroutine post_offline_convergence_diags

!> Extracts members of the offline main control structure. All arguments are optional except
!! the control structure itself
subroutine extract_offline_main(CS, uhtr, vhtr, eatr, ebtr, h_end, accumulated_time, vertical_time, &
                                dt_offline, dt_offline_vertical, skip_diffusion)
  type(offline_transport_CS), target, intent(in   ) :: CS !< Offline control structure
  ! Returned optional arguments
  real, dimension(:,:,:), optional, pointer       :: uhtr !< Remaining zonal mass transport [H m2 ~> m3 or kg]
  real, dimension(:,:,:), optional, pointer       :: vhtr !< Remaining meridional mass transport [H m2 ~> m3 or kg]
  real, dimension(:,:,:), optional, pointer       :: eatr !< Amount of fluid entrained from the layer above within
                                                          !! one time step [H ~> m or kg m-2]
  real, dimension(:,:,:), optional, pointer       :: ebtr !< Amount of fluid entrained from the layer below within
                                                          !! one time step [H ~> m or kg m-2]
  real, dimension(:,:,:), optional, pointer       :: h_end !< Thicknesses at the end of offline timestep
                                                          !! [H ~> m or kg m-2]
  type(time_type),        optional, pointer       :: accumulated_time !< Length of time accumulated in the
                                                          !! current offline interval
  type(time_type),        optional, pointer       :: vertical_time !< The next value of accumulate_time at which to
                                                          !! vertical processes
  real,                   optional, intent(  out) :: dt_offline !< Timestep used for offline tracers [T ~> s]
  real,                   optional, intent(  out) :: dt_offline_vertical !< Timestep used for calls to tracer
                                                          !! vertical physics [T ~> s]
  logical,                optional, intent(  out) :: skip_diffusion !< Skips horizontal diffusion of tracers

  ! Pointers to 3d members
  if (present(uhtr)) uhtr => CS%uhtr
  if (present(vhtr)) vhtr => CS%vhtr
  if (present(eatr)) eatr => CS%eatr
  if (present(ebtr)) ebtr => CS%ebtr
  if (present(h_end)) h_end => CS%h_end

  ! Pointers to integer members which need to be modified
  if (present(accumulated_time)) accumulated_time => CS%accumulated_time
  if (present(vertical_time)) vertical_time => CS%vertical_time

  ! Return value of non-modified integers
  if (present(dt_offline))  dt_offline = CS%dt_offline
  if (present(dt_offline_vertical)) dt_offline_vertical = CS%dt_offline_vertical
  if (present(skip_diffusion)) skip_diffusion = CS%skip_diffusion

end subroutine extract_offline_main

!> Inserts (assigns values to) members of the offline main control structure. All arguments
!! are optional except for the CS itself
subroutine insert_offline_main(CS, ALE_CSp, diabatic_CSp, diag, OBC, tracer_adv_CSp, &
                               tracer_flow_CSp, tracer_Reg, tv, G, GV, x_before_y, debug)
  type(offline_transport_CS), intent(inout) :: CS  !< Offline control structure
  ! Inserted optional arguments
  type(ALE_CS), &
            target, optional, intent(in   ) :: ALE_CSp  !< A pointer to the ALE control structure
  type(diabatic_CS), &
            target, optional, intent(in   ) :: diabatic_CSp !< A pointer to the diabatic control structure
  type(diag_ctrl), &
            target, optional, intent(in   ) :: diag     !< A pointer to the structure that regulates diagnostic output
  type(ocean_OBC_type), &
            target, optional, intent(in   ) :: OBC      !< A pointer to the open boundary condition control structure
  type(tracer_advect_CS), &
            target, optional, intent(in   ) :: tracer_adv_CSp !< A pointer to the tracer advection control structure
  type(tracer_flow_control_CS), &
            target, optional, intent(in   ) :: tracer_flow_CSp !< A pointer to the tracer flow control control structure
  type(tracer_registry_type), &
            target, optional, intent(in   ) :: tracer_Reg !< A pointer to the tracer registry
  type(thermo_var_ptrs), &
            target, optional, intent(in   ) :: tv       !< A structure pointing to various thermodynamic variables
  type(ocean_grid_type), &
            target, optional, intent(in   ) :: G        !< ocean grid structure
  type(verticalGrid_type), &
            target, optional, intent(in   ) :: GV       !< ocean vertical grid structure
  logical,          optional, intent(in   ) :: x_before_y !< Indicates which horizontal direction is advected first
  logical,          optional, intent(in   ) :: debug    !< If true, write verbose debugging messages


  if (present(ALE_CSp))         CS%ALE_CSp => ALE_CSp
  if (present(diabatic_CSp))    CS%diabatic_CSp => diabatic_CSp
  if (present(diag))            CS%diag => diag
  if (present(OBC))             CS%OBC => OBC
  if (present(tracer_adv_CSp))  CS%tracer_adv_CSp => tracer_adv_CSp
  if (present(tracer_flow_CSp)) CS%tracer_flow_CSp => tracer_flow_CSp
  if (present(tracer_Reg))      CS%tracer_Reg => tracer_Reg
  if (present(tv))              CS%tv => tv
  if (present(G))               CS%G => G
  if (present(GV))              CS%GV => GV
  if (present(x_before_y))      CS%x_before_y = x_before_y
  if (present(debug))           CS%debug = debug

end subroutine insert_offline_main

!> Initializes the control structure for offline transport and reads in some of the
! run time parameters from MOM_input
subroutine offline_transport_init(param_file, CS, diabatic_CSp, G, GV, US)

  type(param_file_type),           intent(in) :: param_file !< A structure to parse for run-time parameters
  type(offline_transport_CS),      pointer    :: CS !< Offline control structure
  type(diabatic_CS),               intent(in) :: diabatic_CSp !< The diabatic control structure
  type(ocean_grid_type),   target, intent(in) :: G  !< ocean grid structure
  type(verticalGrid_type), target, intent(in) :: GV !< ocean vertical grid structure
  type(unit_scale_type),   target, intent(in) :: US !< A dimensional unit scaling type

  character(len=40)  :: mdl = "offline_transport"
  character(len=20)  :: redistribute_method

  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB

  is   = G%isc   ; ie   = G%iec  ; js   = G%jsc  ; je   = G%jec ; nz = GV%ke
  isd  = G%isd   ; ied  = G%ied  ; jsd  = G%jsd  ; jed  = G%jed
  IsdB = G%IsdB  ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call callTree_enter("offline_transport_init, MOM_offline_control.F90")

  if (associated(CS)) then
    call MOM_error(WARNING, "offline_transport_init called with an associated "// &
      "control structure.")
    return
  endif
  allocate(CS)
  call log_version(param_file, mdl,version, "This module allows for tracers to be run offline")

  ! Determining the internal unit scaling factors for this run.
  CS%US => US

  ! Parse MOM_input for offline control
  call get_param(param_file, mdl, "OFFLINEDIR", CS%offlinedir, &
    "Input directory where the offline fields can be found",  fail_if_missing = .true.)
  call get_param(param_file, mdl, "OFF_SUM_FILE", CS%sum_file, &
    "Filename where the accumulated fields can be found",     fail_if_missing = .true.)
  call get_param(param_file, mdl, "OFF_SNAP_FILE", CS%snap_file, &
    "Filename where snapshot fields can be found",            fail_if_missing = .true.)
  call get_param(param_file, mdl, "OFF_MEAN_FILE", CS%mean_file, &
    "Filename where averaged fields can be found",            fail_if_missing = .true.)
  call get_param(param_file, mdl, "OFF_SURF_FILE", CS%surf_file, &
    "Filename where averaged fields can be found",            fail_if_missing = .true.)
  call get_param(param_file, mdl, "NUMTIME", CS%numtime, &
    "Number of timelevels in offline input files",            fail_if_missing = .true.)
  call get_param(param_file, mdl, "NK_INPUT", CS%nk_input, &
    "Number of vertical levels in offline input files", default = nz)
  call get_param(param_file, mdl, "DT_OFFLINE", CS%dt_offline, &
    "Length of time between reading in of input fields", units='s', scale=US%s_to_T, fail_if_missing = .true.)
  call get_param(param_file, mdl, "DT_OFFLINE_VERTICAL", CS%dt_offline_vertical, &
    "Length of the offline timestep for tracer column sources/sinks " //&
    "This should be set to the length of the coupling timestep for " //&
    "tracers which need shortwave fluxes", units="s", scale=US%s_to_T, fail_if_missing = .true.)
  call get_param(param_file, mdl, "START_INDEX", CS%start_index, &
    "Which time index to start from", default=1)
  call get_param(param_file, mdl, "FIELDS_ARE_OFFSET", CS%fields_are_offset, &
    "True if the time-averaged fields and snapshot fields "//&
    "are offset by one time level", default=.false.)
  call get_param(param_file, mdl, "REDISTRIBUTE_METHOD", redistribute_method, &
    "Redistributes any remaining horizontal fluxes throughout "    //&
    "the rest of water column. Options are 'barotropic' which "    //&
    "evenly distributes flux throughout the entire water column, " //&
    "'upwards' which adds the maximum of the remaining flux in "   //&
    "each layer above, both which first applies upwards and then " //&
    "barotropic, and 'none' which does no redistribution", &
    default='barotropic')
  call get_param(param_file, mdl, "NUM_OFF_ITER", CS%num_off_iter, &
    "Number of iterations to subdivide the offline tracer advection and diffusion", &
    default = 60)
  call get_param(param_file, mdl, "OFF_ALE_MOD", CS%off_ale_mod, &
    "Sets how many horizontal advection steps are taken before an ALE " //&
    "remapping step is done. 1 would be x->y->ALE, 2 would be"           //&
    "x->y->x->y->ALE", default = 1)
  call get_param(param_file, mdl, "PRINT_ADV_OFFLINE", CS%print_adv_offline, &
    "Print diagnostic output every advection subiteration",default=.false.)
  call get_param(param_file, mdl, "SKIP_DIFFUSION_OFFLINE", CS%skip_diffusion, &
    "Do not do horizontal diffusion",default=.false.)
  call get_param(param_file, mdl, "READ_SW", CS%read_sw, &
    "Read in shortwave radiation field instead of using values from the coupler"//&
    "when in offline tracer mode",default=.false.)
  call get_param(param_file, mdl, "READ_MLD", CS%read_mld, &
    "Read in mixed layer depths for tracers which exchange with the atmosphere"//&
    "when in offline tracer mode",default=.false.)
  call get_param(param_file, mdl, "MLD_VAR_NAME", CS%mld_var_name, &
    "Name of the variable containing the depth of active mixing",&
    default='ePBL_h_ML')
  call get_param(param_file, mdl, "OFFLINE_ADD_DIURNAL_SW", CS%diurnal_sw, &
    "Adds a synthetic diurnal cycle in the same way that the ice " // &
    "model would have when time-averaged fields of shortwave "    // &
    "radiation are read in", default=.false.)
  call get_param(param_file, mdl, "KD_MAX", CS%Kd_max, &
    "The maximum permitted increment for the diapycnal "//&
    "diffusivity from TKE-based parameterizations, or a "//&
    "negative value for no limit.", units="m2 s-1", default=-1.0)
  call get_param(param_file, mdl, "MIN_RESIDUAL_TRANSPORT", CS%min_residual, &
    "How much remaining transport before the main offline advection "// &
    "is exited. The default value corresponds to about 1 meter of "  // &
    "difference in a grid cell", default = 1.e9)
  call get_param(param_file, mdl, "READ_ALL_TS_UVH", CS%read_all_ts_uvh,  &
    "Reads all time levels of a subset of the fields necessary to run "    //      &
    "the model offline. This can require a large amount of memory "//      &
    "and will make initialization very slow. However, for offline "//      &
    "runs spanning more than a year this can reduce total I/O overhead",    &
    default = .false.)

  ! Concatenate offline directory and file names
  CS%snap_file = trim(CS%offlinedir)//trim(CS%snap_file)
  CS%mean_file = trim(CS%offlinedir)//trim(CS%mean_file)
  CS%sum_file = trim(CS%offlinedir)//trim(CS%sum_file)
  CS%surf_file = trim(CS%offlinedir)//trim(CS%surf_file)

  CS%num_vert_iter = CS%dt_offline/CS%dt_offline_vertical

  ! Map redistribute_method onto logicals in CS
  select case (redistribute_method)
    case ('barotropic')
      CS%redistribute_barotropic = .true.
      CS%redistribute_upwards    = .false.
    case ('upwards')
      CS%redistribute_barotropic = .false.
      CS%redistribute_upwards    = .true.
    case ('both')
      CS%redistribute_barotropic = .true.
      CS%redistribute_upwards    = .true.
    case ('none')
      CS%redistribute_barotropic = .false.
      CS%redistribute_upwards    = .false.
  end select

  ! Set the accumulated time to zero
  CS%accumulated_time = real_to_time(0.0)
  CS%vertical_time = CS%accumulated_time
  ! Set the starting read index for time-averaged and snapshotted fields
  CS%ridx_sum = CS%start_index
  if (CS%fields_are_offset) CS%ridx_snap = next_modulo_time(CS%start_index,CS%numtime)
  if (.not. CS%fields_are_offset) CS%ridx_snap = CS%start_index

  ! Copy members from other modules
  call extract_diabatic_member(diabatic_CSp, opacity_CSp=CS%opacity_CSp, optics_CSp=CS%optics, &
                               diabatic_aux_CSp=CS%diabatic_aux_CSp, &
                               evap_CFL_limit=CS%evap_CFL_limit, &
                               minimum_forcing_depth=CS%minimum_forcing_depth)

  ! Grid pointer assignments
  CS%G  => G
  CS%GV => GV

  ! Allocate arrays
  allocate(CS%uhtr(IsdB:IedB,jsd:jed,nz))   ; CS%uhtr(:,:,:) = 0.0
  allocate(CS%vhtr(isd:ied,JsdB:JedB,nz))   ; CS%vhtr(:,:,:) = 0.0
  allocate(CS%eatr(isd:ied,jsd:jed,nz))          ; CS%eatr(:,:,:) = 0.0
  allocate(CS%ebtr(isd:ied,jsd:jed,nz))          ; CS%ebtr(:,:,:) = 0.0
  allocate(CS%h_end(isd:ied,jsd:jed,nz))         ; CS%h_end(:,:,:) = 0.0
  allocate(CS%netMassOut(G%isd:G%ied,G%jsd:G%jed)) ; CS%netMassOut(:,:) = 0.0
  allocate(CS%netMassIn(G%isd:G%ied,G%jsd:G%jed))  ; CS%netMassIn(:,:) = 0.0
  allocate(CS%Kd(isd:ied,jsd:jed,nz+1)) ; CS%Kd = 0.
  if (CS%read_mld) then
    allocate(CS%mld(G%isd:G%ied,G%jsd:G%jed)) ; CS%mld(:,:) = 0.0
  endif

  if (CS%read_all_ts_uvh) then
    call read_all_input(CS)
  endif

  ! Initialize ids for clocks used in offline routines
  CS%id_clock_read_fields =      cpu_clock_id('(Offline read fields)',grain=CLOCK_MODULE)
  CS%id_clock_offline_diabatic = cpu_clock_id('(Offline diabatic)',grain=CLOCK_MODULE)
  CS%id_clock_offline_adv =      cpu_clock_id('(Offline transport)',grain=CLOCK_MODULE)
  CS%id_clock_redistribute =     cpu_clock_id('(Offline redistribute)',grain=CLOCK_MODULE)

  call callTree_leave("offline_transport_init")

end subroutine offline_transport_init

!> Coordinates the allocation and reading in all time levels of uh, vh, hend, temp, and salt from files. Used
!! when read_all_ts_uvh
subroutine read_all_input(CS)
  type(offline_transport_CS), intent(inout)  :: CS !< Control structure for offline module

  integer :: is, ie, js, je, isd, ied, jsd, jed, nz, t, ntime
  integer :: IsdB, IedB, JsdB, JedB

  nz = CS%GV%ke ; ntime = CS%numtime
  isd  = CS%G%isd   ; ied  = CS%G%ied  ; jsd  = CS%G%jsd  ; jed  = CS%G%jed
  IsdB = CS%G%IsdB  ; IedB = CS%G%IedB ; JsdB = CS%G%JsdB ; JedB = CS%G%JedB

  ! Extra safety check that we're not going to overallocate any arrays
  if (CS%read_all_ts_uvh) then
    if (allocated(CS%uhtr_all)) call MOM_error(FATAL, "uhtr_all is already allocated")
    if (allocated(CS%vhtr_all)) call MOM_error(FATAL, "vhtr_all is already allocated")
    if (allocated(CS%hend_all)) call MOM_error(FATAL, "hend_all is already allocated")
    if (allocated(CS%temp_all)) call MOM_error(FATAL, "temp_all is already allocated")
    if (allocated(CS%salt_all)) call MOM_error(FATAL, "salt_all is already allocated")

    allocate(CS%uhtr_all(IsdB:IedB,jsd:jed,nz,ntime))     ; CS%uhtr_all(:,:,:,:) = 0.0
    allocate(CS%vhtr_all(isd:ied,JsdB:JedB,nz,ntime))     ; CS%vhtr_all(:,:,:,:) = 0.0
    allocate(CS%hend_all(isd:ied,jsd:jed,nz,ntime))       ; CS%hend_all(:,:,:,:) = 0.0
    allocate(CS%temp_all(isd:ied,jsd:jed,nz,1:ntime))     ; CS%temp_all(:,:,:,:) = 0.0
    allocate(CS%salt_all(isd:ied,jsd:jed,nz,1:ntime))     ; CS%salt_all(:,:,:,:) = 0.0

    call MOM_mesg("Reading in uhtr, vhtr, h_start, h_end, temp, salt")
    do t = 1,ntime
      call MOM_read_vector(CS%snap_file, 'uhtr_sum', 'vhtr_sum', CS%uhtr_all(:,:,1:CS%nk_input,t), &
                       CS%vhtr_all(:,:,1:CS%nk_input,t), CS%G%Domain, timelevel=t)
      call MOM_read_data(CS%snap_file,'h_end', CS%hend_all(:,:,1:CS%nk_input,t), CS%G%Domain, &
        timelevel=t, position=CENTER)
      call MOM_read_data(CS%mean_file,'temp', CS%temp_all(:,:,1:CS%nk_input,t), CS%G%Domain, &
        timelevel=t, position=CENTER)
      call MOM_read_data(CS%mean_file,'salt', CS%salt_all(:,:,1:CS%nk_input,t), CS%G%Domain, &
        timelevel=t, position=CENTER)
    enddo
  endif

end subroutine read_all_input

!> Deallocates (if necessary) arrays within the offline control structure
subroutine offline_transport_end(CS)
  type(offline_transport_CS), pointer :: CS !< Control structure for offline module

  ! Explicitly allocate all allocatable arrays
  deallocate(CS%uhtr)
  deallocate(CS%vhtr)
  deallocate(CS%eatr)
  deallocate(CS%ebtr)
  deallocate(CS%h_end)
  deallocate(CS%netMassOut)
  deallocate(CS%netMassIn)
  deallocate(CS%Kd)
  if (CS%read_mld) deallocate(CS%mld)
  if (CS%read_all_ts_uvh) then
    deallocate(CS%uhtr_all)
    deallocate(CS%vhtr_all)
    deallocate(CS%hend_all)
    deallocate(CS%temp_all)
    deallocate(CS%salt_all)
  endif

  deallocate(CS)

end subroutine offline_transport_end

!> \namespace mom_offline_main
!! \section offline_overview Offline Tracer Transport in MOM6
!!  'Offline tracer modeling' uses physical fields (e.g. mass transports and layer thicknesses) saved
!!  from a previous integration of the physical model to transport passive tracers. These fields are
!!  accumulated or averaged over a period of time (in this test case, 1 day) and used to integrate
!!  portions of the MOM6 code base that handle the 3d advection and diffusion of passive tracers.
!!
!!  The distribution of tracers in the ocean modeled offline should not be expected to match an online
!!  simulation. Accumulating transports over more than one online model timestep implicitly assumes
!!  homogeneity over that time period and essentially aliases over processes that occur with higher
!!  frequency. For example, consider the case of a surface boundary layer with a strong diurnal cycle.
!!  An offline simulation with a 1 day timestep, captures the net transport into or out of that layer,
!!  but not the exact cycling. This effective aliasing may also complicate online model configurations
!!  which strongly-eddying regions. In this case, the offline model timestep must be limited to some
!!  fraction of the eddy correlation timescale. Lastly, the nonlinear advection scheme which applies
!!  limited mass-transports over a sequence of iterations means that tracers are not transported along
!!  exactly the same path as they are in the online model.
!!
!!  This capability has currently targeted the Baltic_ALE_z test case, though some work has also been
!!  done with the OM4 1/2 degree configuration. Work is ongoing to develop recommendations and best
!!  practices for investigators seeking to use MOM6 for offline tracer modeling.
!!
!!  \section offline_technical Implementation of offline routine in MOM6
!!
!!  The subroutine step_tracers that coordinates this can be found in MOM.F90 and is only called
!!  using the solo ocean driver. This is to avoid issues with coupling to other climate components
!!  that may be relying on fluxes from the ocean to be coupled more often than the offline time step.
!!  Other routines related to offline tracer modeling can be found in tracers/MOM_offline_control.F90
!!
!!  As can also be seen in the comments for the step_tracers subroutine, an offline time step
!!  comprises the following steps:
!!        -#  Using the layer thicknesses and tracer concentrations from the previous timestep,
!!            half of the accumulated vertical mixing (eatr and ebtr) is applied in the call to
!!            tracer_column_fns.
!!            For tracers whose source/sink terms need dt, this value is set to 1/2 dt_offline
!!        -#  Half of the accumulated surface freshwater fluxes are applied
!!        START ITERATION
!!        -#  Accumulated mass fluxes are used to do horizontal transport. The number of iterations
!!            used in advect_tracer is limited to 2 (e.g x->y->x->y). The remaining mass fluxes are
!!            stored for later use and resulting layer thicknesses fed into the next step
!!        -#  Tracers and the h-grid are regridded and remapped in a call to ALE. This allows for
!!            layers which might 'vanish' because of horizontal mass transport to be 'reinflated'
!!            and essentially allows for the vertical transport of tracers
!!        -#  Check that transport is done if the remaining mass fluxes equals 0 or if the max
!!            number of iterations has been reached
!!        END ITERATION
!!        -#  Repeat steps 1 and 2
!!        -#  Redistribute any residual mass fluxes that remain after the advection iterations
!!            in a barotropic manner, progressively upward through the water column.
!!        -#  Force a remapping to the stored layer thicknesses that correspond to the snapshot of
!!            the online model at the end of an accumulation interval
!!        -#  Reset T/S and h to their stored snapshotted values to prevent model drift
!!
!!  \section  offline_evaluation Evaluating the utility of an offline tracer model
!!  How well an offline tracer model can be used as an alternative to integrating tracers online
!!  with the prognostic model must be evaluated for each application. This efficacy may be related
!!  to the native coordinate of the online model, to the length of the offline timestep, and to the
!!  behavior of the tracer itself.
!!
!!  A framework for formally regression testing the offline capability still needs to be developed.
!!  However, as a simple way of testing whether the offline model is nominally behaving as expected,
!!  the total inventory of the advection test tracers (tr1, tr2, etc.) should be conserved between
!!  time steps except for the last 4 decimal places. As a general guideline, an offline timestep of
!!  5 days or less.
!!
!!  \section offline_parameters Runtime parameters for offline tracers
!!    - OFFLINEDIR:    Input directory where the offline fields can be found
!!    - OFF_SUM_FILE:  Filename where the accumulated fields can be found (e.g. horizontal mass transports)
!!    - OFF_SNAP_FILE: Filename where snapshot fields can be found (e.g. end of timestep layer thickness)
!!    - START_INDEX:   Which timelevel of the input files to read first
!!    - NUMTIME:       How many timelevels to read before 'looping' back to 1
!!    - FIELDS_ARE_OFFSET: True if the time-averaged fields and snapshot fields are offset by one
!!                        time level, probably not needed
!!    -NUM_OFF_ITER:  Maximum number of iterations to do for the nonlinear advection scheme
!!    -REDISTRIBUTE_METHOD: Redistributes any remaining horizontal fluxes throughout the rest of water column.
!!                          Options are 'barotropic' which "evenly distributes flux throughout the entire water
!!                          column,'upwards' which adds the maximum of the remaining flux in each layer above,
!!                          and 'none' which does no redistribution"

end module MOM_offline_main

