!> Initialization for the dyed_channel configuration
module dyed_channel_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_dyn_horgrid,     only : dyn_horgrid_type
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,     only : get_param, log_version, param_file_type
use MOM_get_input,       only : directories
use MOM_grid,            only : ocean_grid_type
use MOM_open_boundary,   only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary,   only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S, OBC_DIRECTION_E
use MOM_open_boundary,   only : OBC_segment_type, register_segment_tracer
use MOM_open_boundary,   only : OBC_registry_type, register_OBC
use MOM_time_manager,    only : time_type, time_type_to_real
use MOM_tracer_registry, only : tracer_registry_type, tracer_name_lookup
use MOM_tracer_registry, only : tracer_type
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public dyed_channel_set_OBC_tracer_data, dyed_channel_OBC_end
public register_dyed_channel_OBC, dyed_channel_update_flow

!> Control structure for dyed-channel open boundaries.
type, public :: dyed_channel_OBC_CS ; private
  real :: zonal_flow = 8.57         !< Mean inflow [L T-1 ~> m s-1]
  real :: tidal_amp = 0.0           !< Sloshing amplitude [L T-1 ~> m s-1]
  real :: frequency  = 0.0          !< Sloshing frequency [T-1 ~> s-1]
  logical :: OBC_transport_bug      !< If true and specified open boundary conditions are being
                                    !! used, use a 1 m (if Boussienesq) or 1 kg m-2 layer thickness
                                    !! instead of the actual thickness.
end type dyed_channel_OBC_CS

integer :: ntr = 0 !< Number of dye tracers
                   !! \todo This is a module variable. Move this variable into the control structure.

contains

!> Add dyed channel to OBC registry.
function register_dyed_channel_OBC(param_file, CS, US, OBC_Reg)
  type(param_file_type),     intent(in) :: param_file !< parameter file.
  type(dyed_channel_OBC_CS), pointer    :: CS         !< Dyed channel control structure.
  type(unit_scale_type),     intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),   pointer    :: OBC_Reg    !< OBC registry.
  ! Local variables
  logical                               :: register_dyed_channel_OBC
  character(len=32)  :: casename = "dyed channel"     ! This case's name.
  character(len=40)  :: mdl = "register_dyed_channel_OBC" ! This subroutine's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "register_dyed_channel_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  call get_param(param_file, mdl, "CHANNEL_MEAN_FLOW", CS%zonal_flow, &
                 "Mean zonal flow imposed at upstream open boundary.", &
                 units="m/s", default=8.57, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "CHANNEL_TIDAL_AMP", CS%tidal_amp, &
                 "Sloshing amplitude imposed at upstream open boundary.", &
                 units="m/s", default=0.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "CHANNEL_FLOW_FREQUENCY", CS%frequency, &
                 "Frequency of oscillating zonal flow.", &
                 units="s-1", default=0.0, scale=US%T_to_s)
  call get_param(param_file, mdl, "CHANNEL_FLOW_OBC_TRANSPORT_BUG", CS%OBC_transport_bug, &
                 "If true and specified open boundary conditions are being used, use a 1 m "//&
                 "(if Boussienesq) or 1 kg m-2 layer thickness instead of the actual thickness.", &
                 default=.true.)  !### Change the default to False.

  ! Register the open boundaries.
  call register_OBC(casename, param_file, OBC_Reg)
  register_dyed_channel_OBC = .true.

end function register_dyed_channel_OBC

!> Clean up the dyed_channel OBC from registry.
subroutine dyed_channel_OBC_end(CS)
  type(dyed_channel_OBC_CS), pointer :: CS    !< Dyed channel control structure.

  if (associated(CS)) then
    deallocate(CS)
  endif
end subroutine dyed_channel_OBC_end

!> This subroutine sets the dye and flow properties at open boundary conditions.
subroutine dyed_channel_set_OBC_tracer_data(OBC, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                                                !! to parse for model parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.
  ! Local variables
  character(len=40)  :: mdl = "dyed_channel_set_OBC_tracer_data" ! This subroutine's name.
  character(len=80)  :: name, longname
  integer :: m, n, ntr_id
  real :: dye  ! Inflow dye concentrations [arbitrary]
  type(tracer_type), pointer      :: tr_ptr => NULL()

  if (.not.associated(OBC)) call MOM_error(FATAL, 'dyed_channel_initialization.F90: '// &
        'dyed_channel_set_OBC_data() was called but OBC type was not initialized!')

  call get_param(param_file, mdl, "NUM_DYE_TRACERS", ntr, &
                 "The number of dye tracers in this run. Each tracer "//&
                 "should have a separate boundary segment.", default=0,   &
                 do_not_log=.true.)

  if (OBC%number_of_segments < ntr) then
    call MOM_error(WARNING, "Error in dyed_obc segment setup")
    return   !!! Need a better error message here
  endif

! ! Set the inflow values of the dyes, one per segment.
! ! We know the order: north, south, east, west
  do m=1,ntr
    write(name,'("dye_",I2.2)') m
    write(longname,'("Concentration of dyed_obc Tracer ",I2.2, " on segment ",I2.2)') m, m
    call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, name)

    do n=1,OBC%number_of_segments
      if (n == m) then
        dye = 1.0
      else
        dye = 0.0
      endif
      call register_segment_tracer(tr_ptr, ntr_id, param_file, GV, &
                                   OBC%segment(n), OBC_scalar=dye)
    enddo
  enddo

end subroutine dyed_channel_set_OBC_tracer_data

!> This subroutine updates the long-channel flow
subroutine dyed_channel_update_flow(OBC, CS, G, GV, US, h, Time)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(dyed_channel_OBC_CS),  pointer    :: CS  !< Dyed channel control structure.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness [H ~> m or kg m-2]
  type(time_type),            intent(in) :: Time !< model time.

  ! Local variables
  real :: flow      ! The OBC velocity [L T-1 ~> m s-1]
  real :: PI        ! 3.1415926535... [nondim]
  real :: time_sec  ! The elapsed time since the start of the calendar [T ~> s]
  real :: fixed_thickness ! A fixed layer thickness, hard-coded to 1 mks unit, that is used to
                    ! reproduce a bug with the older versions of this code [H ~> m or kg m-2]
  logical :: cross_channel  ! True if the segment runs across the channel
  integer :: turns    ! Number of index quarter turns
  integer :: i, j, k, l_seg, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB, is, ie, js, je
  type(OBC_segment_type), pointer :: segment => NULL()

  if (.not.associated(OBC)) call MOM_error(FATAL, 'dyed_channel_initialization.F90: '// &
        'dyed_channel_update_flow() was called but OBC type was not initialized!')

  time_sec = US%s_to_T * time_type_to_real(Time)
  PI = 4.0*atan(1.0)

  turns = modulo(G%HI%turns, 4)

  do l_seg=1, OBC%number_of_segments
    segment => OBC%segment(l_seg)
    if (.not. segment%on_pe) cycle
    if (segment%gradient) cycle
    if (segment%oblique .and. (.not. segment%nudged) .and. (.not. segment%Flather)) cycle

    if (CS%frequency == 0.0) then
      flow = CS%zonal_flow
    else
      flow = CS%zonal_flow + CS%tidal_amp * cos(2 * PI * CS%frequency * time_sec)
    endif
    if ((turns==2) .or. (turns==3)) flow = -1.0 * flow

    isd = segment%HI%isd ; ied = segment%HI%ied
    jsd = segment%HI%jsd ; jed = segment%HI%jed
    IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
    JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
    if (segment%is_E_or_W) then
      is = IsdB ; ie = IedB ; js = jsd ; je = jed
    else
      is = isd ; ie = ied ; js = JsdB ; je = JedB
    endif
    cross_channel = ((segment%is_E_or_W .and. ((turns==0) .or. (turns==2))) .or. &
                     (segment%is_N_or_S .and. ((turns==1) .or. (turns==3))))

    if ((segment%specified .or. segment%nudged) .and. cross_channel) then
      do k=1,GV%ke ; do j=js,je ; do I=is,ie
        segment%normal_vel(I,j,k) = flow
      enddo ; enddo ; enddo
    endif

    if (segment%specified .and. cross_channel) then
      if (CS%OBC_transport_bug) then
        fixed_thickness = 1.0 / GV%H_to_mks  ! This replicates the prevoius answers without rescaling.
        if ((segment%direction == OBC_DIRECTION_W) .or. (segment%direction == OBC_DIRECTION_E)) then
          do k=1,GV%ke ; do j=jsd,jed ; do I=IsdB,IedB
            segment%normal_trans(I,j,k) = flow * G%dyCu(I,j) * fixed_thickness
          enddo ; enddo ; enddo
        elseif ((segment%direction == OBC_DIRECTION_S) .or. (segment%direction == OBC_DIRECTION_N)) then
          do k=1,GV%ke ; do J=JsdB,JedB ; do i=isd,ied
            segment%normal_trans(i,J,k) = flow * G%dxCv(i,J) * fixed_thickness
          enddo ; enddo ; enddo
        endif
      else
        if (segment%direction == OBC_DIRECTION_W) then
          do k=1,GV%ke ; do j=jsd,jed ; do I=IsdB,IedB
            segment%normal_trans(I,j,k) = flow * G%dyCu(I,j) * h(i+1,j,k)
          enddo ; enddo ; enddo
        elseif (segment%direction == OBC_DIRECTION_E) then
          do k=1,GV%ke ; do j=jsd,jed ; do I=IsdB,IedB
            segment%normal_trans(I,j,k) = flow * G%dyCu(I,j) * h(i,j,k)
          enddo ; enddo ; enddo
        elseif (segment%direction == OBC_DIRECTION_S) then
          do k=1,GV%ke ; do J=JsdB,JedB ; do i=isd,ied
            segment%normal_trans(i,J,k) = flow * G%dxCv(i,J) * h(i,j+1,k)
          enddo ; enddo ; enddo
        elseif (segment%direction == OBC_DIRECTION_N) then
          do k=1,GV%ke ; do J=JsdB,JedB ; do i=isd,ied
            segment%normal_trans(i,J,k) = flow * G%dxCv(i,J) * h(i,j,k)
          enddo ; enddo ; enddo
        endif
      endif
    endif

    if (cross_channel) then
      do j=js,je ; do I=is,ie
        segment%normal_vel_bt(I,j) = flow
      enddo ; enddo
    else
      do J=js,je ; do i=is,ie
        segment%normal_vel_bt(i,J) = 0.0
      enddo ; enddo
    endif

  enddo

end subroutine dyed_channel_update_flow

!> \namespace dyed_channel_initialization
!!
!! Setting dyes, one for painting the inflow on each side.
end module dyed_channel_initialization
