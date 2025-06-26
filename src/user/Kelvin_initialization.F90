!> Configures the model for the Kelvin wave experiment.
!!
!! Kelvin = coastally-trapped Kelvin waves from the ROMS examples.
!! Initialize with level surfaces and drive the wave in at the west,
!! radiate out at the east.
module Kelvin_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE
use MOM_open_boundary,  only : OBC_segment_type, register_OBC, rotate_OBC_segment_direction
use MOM_open_boundary,  only : OBC_DIRECTION_N, OBC_DIRECTION_E, OBC_DIRECTION_S, OBC_DIRECTION_W
use MOM_open_boundary,  only : OBC_registry_type
use MOM_unit_scaling,   only : unit_scale_type
use MOM_verticalGrid,   only : verticalGrid_type
use MOM_time_manager,   only : time_type, time_type_to_real

implicit none ; private

#include <MOM_memory.h>

public Kelvin_set_OBC_data, Kelvin_initialize_topography
public register_Kelvin_OBC, Kelvin_OBC_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Control structure for Kelvin wave open boundaries.
type, public :: Kelvin_OBC_CS ; private
  integer :: mode = 0          !< Vertical mode
  real :: coast_angle = 0   !< Angle of coastline [rad]
  real :: coast_offset1 = 0 !< Longshore distance to coastal angle [L ~> m]
  real :: coast_offset2 = 0 !< Offshore distance to coastal angle [L ~> m]
  real :: H0 = 0            !< Bottom depth [Z ~> m]
  real :: F_0               !< Coriolis parameter [T-1 ~> s-1]
  real :: rho_range         !< Density range [R ~> kg m-3]
  real :: rho_0             !< Mean density [R ~> kg m-3]
  real :: wave_period       !< Period of the mode-0 waves [T ~> s]
  real :: ssh_amp           !< Amplitude of the sea surface height forcing for mode-0 waves [Z ~> m]
  real :: inflow_amp        !< Amplitude of the boundary velocity forcing for internal waves [L T-1 ~> m s-1]
  real :: OBC_nudging_time  !< The timescale with which the inflowing open boundary velocities are nudged toward
                            !! their intended values with the Kelvin wave test case [T ~> s], or a negative
                            !! value to retain the value that is set when the OBC segments are initialized.
  logical :: indexing_bugs  !< If true, retain several horizontal indexing bugs that were in the
                            !! original version of Kelvin_set_OBC_data.
end type Kelvin_OBC_CS

! This include declares and sets the variable "version".
#include "version_variable.h"

contains

!> Add Kelvin wave to OBC registry.
logical function register_Kelvin_OBC(param_file, CS, US, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(Kelvin_OBC_CS),      pointer    :: CS         !< Kelvin wave control structure.
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.

  ! Local variables
  logical :: enable_bugs  ! If true, the defaults for recently added bug-fix flags are set to
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  character(len=40)  :: mdl = "register_Kelvin_OBC"  !< This subroutine's name.
  character(len=32)  :: casename = "Kelvin wave"     !< This case's name.
  character(len=200) :: config

  if (associated(CS)) then
    call MOM_error(WARNING, "register_Kelvin_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "KELVIN_WAVE_MODE", CS%mode, &
                 "Vertical Kelvin wave mode imposed at upstream open boundary.", &
                 default=0)
  call get_param(param_file, mdl, "F_0", CS%F_0, &
                 default=0.0, units="s-1", scale=US%T_to_s, do_not_log=.true.)
  call get_param(param_file, mdl, "TOPO_CONFIG", config, fail_if_missing=.true., do_not_log=.true.)
  if (trim(config) == "Kelvin") then
    call get_param(param_file, mdl, "ROTATED_COAST_OFFSET_1", CS%coast_offset1, &
                   "The distance along the southern and northern boundaries "//&
                   "at which the coasts angle in.", &
                   units="km", default=100.0, scale=1.0e3*US%m_to_L)
    call get_param(param_file, mdl, "ROTATED_COAST_OFFSET_2", CS%coast_offset2, &
                   "The distance from the southern and northern boundaries "//&
                   "at which the coasts angle in.", &
                   units="km", default=10.0, scale=1.0e3*US%m_to_L)
    call get_param(param_file, mdl, "ROTATED_COAST_ANGLE", CS%coast_angle, &
                   "The angle of the southern bondary beyond X=ROTATED_COAST_OFFSET.", &
                   units="degrees", default=11.3, scale=atan(1.0)/45.) ! Convert to radians
  else
    CS%coast_offset1 = 0.0 ; CS%coast_offset2 = 0.0 ; CS%coast_angle = 0.0
  endif
  if (CS%mode == 0) then
    call get_param(param_file, mdl, "KELVIN_WAVE_PERIOD", CS%wave_period, &
                   "The period of the Kelvin wave forcing at the open boundaries.  "//&
                   "The default value is the M2 tide period.", &
                   units="s", default=12.42*3600.0, scale=US%s_to_T)
    call get_param(param_file, mdl, "KELVIN_WAVE_SSH_AMP", CS%ssh_amp, &
                   "The amplitude of the Kelvin wave sea surface height anomaly forcing "//&
                   "at the open boundaries.", units="m", default=1.0, scale=US%m_to_Z)
  else
    call get_param(param_file, mdl, "DENSITY_RANGE", CS%rho_range, &
                   units="kg m-3", default=2.0, scale=US%kg_m3_to_R, do_not_log=.true.)
    call get_param(param_file, mdl, "RHO_0", CS%rho_0, &
                   units="kg m-3", default=1035.0, scale=US%kg_m3_to_R, do_not_log=.true.)
    call get_param(param_file, mdl, "MAXIMUM_DEPTH", CS%H0, &
                   units="m", default=1000.0, scale=US%m_to_Z, do_not_log=.true.)
    call get_param(param_file, mdl, "KELVIN_WAVE_INFLOW_AMP", CS%inflow_amp, &
                   "The amplitude of the Kelvin wave sea surface inflow velocity forcing "//&
                   "at the open boundaries.", units="m s-1", default=1.0, scale=US%m_s_to_L_T)
  endif

  call get_param(param_file, mdl, "KELVIN_WAVE_VEL_NUDGING_TIMESCALE", CS%OBC_nudging_time, &
                 "The timescale with which the inflowing open boundary velocities are nudged toward "//&
                 "their intended values with the Kelvin wave test case, or a negative value to keep "//&
                 "the value that is set when the OBC segments are initialized.", &
                 units="s", default=1.0/(0.3*86400.), scale=US%s_to_T)
                 !### Change the default nudging timescale to -1. or another value?
  call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
  call get_param(param_file, mdl, "KELVIN_SET_OBC_INDEXING_BUGS", CS%indexing_bugs, &
                 "If true, retain several horizontal indexing bugs that were in the original "//&
                 "version of Kelvin_set_OBC_data.", default=enable_bugs)

  ! Register the Kelvin open boundary.
  call register_OBC(casename, param_file, OBC_Reg)
  register_Kelvin_OBC = .true.

  ! TODO: Revisit and correct the internal Kelvin wave test case.
  ! Specifically, using wave_speed() and investigating adding eta_anom
  ! noted in the comments below.

end function register_Kelvin_OBC

!> Clean up the Kelvin wave OBC from registry.
subroutine Kelvin_OBC_end(CS)
  type(Kelvin_OBC_CS), pointer    :: CS         !< Kelvin wave control structure.

  if (associated(CS)) then
    deallocate(CS)
  endif
end subroutine Kelvin_OBC_end

! -----------------------------------------------------------------------------
!> This subroutine sets up the Kelvin topography and land mask
subroutine Kelvin_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
  character(len=40)  :: mdl = "Kelvin_initialize_topography" ! This subroutine's name.
  real :: min_depth     ! The minimum and maximum depths [Z ~> m].
  real :: coast_angle   ! Angle of coastline [rad]
  real :: coast_offset1 ! Longshore distance to coastal angle [L ~> m]
  real :: coast_offset2 ! Offshore distance to coastal angle [L ~> m]
  integer :: i, j

  call MOM_mesg("  Kelvin_initialization.F90, Kelvin_initialize_topography: setting topography", 5)

  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "ROTATED_COAST_OFFSET_1", coast_offset1, &
                 units="km", default=100.0, do_not_log=.true.)
  call get_param(param_file, mdl, "ROTATED_COAST_OFFSET_2", coast_offset2, &
                 units="km", default=10.0, do_not_log=.true.)
  call get_param(param_file, mdl, "ROTATED_COAST_ANGLE", coast_angle, &
                 units="degrees", default=11.3, scale=(atan(1.0)/45.), do_not_log=.true.) ! Convert to radians

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    D(i,j) = max_depth
    ! Southern side
    if ((G%geoLonT(i,j) - G%west_lon > coast_offset1) .AND. &
        (atan2(G%geoLatT(i,j) - G%south_lat + coast_offset2, &
               G%geoLonT(i,j) - G%west_lon - coast_offset1) < coast_angle)) &
      D(i,j) = 0.5*min_depth
    ! Northern side
    if ((G%geoLonT(i,j) - G%west_lon < G%len_lon - coast_offset1) .AND. &
        (atan2(G%len_lat + G%south_lat + coast_offset2 - G%geoLatT(i,j), &
               G%len_lon + G%west_lon - coast_offset1 - G%geoLonT(i,j)) < coast_angle)) &
      D(i,j) = 0.5*min_depth

    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
  enddo ; enddo

end subroutine Kelvin_initialize_topography

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine Kelvin_set_OBC_data(OBC, CS, G, GV, US, h, Time)
  type(ocean_OBC_type),    pointer    :: OBC  !< This open boundary condition type specifies
                                              !! whether, where, and what open boundary
                                              !! conditions are used.
  type(Kelvin_OBC_CS),     pointer    :: CS   !< Kelvin wave control structure.
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness [H ~> m or kg m-2].
  type(time_type),         intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the Kelvin example.
  real :: time_sec     ! The time in the run [T ~> s]
  real :: cff          ! The wave speed [L T-1 ~> m s-1]
  real :: N0           ! Brunt-Vaisala frequency times a rescaling of slopes [L Z-1 T-1 ~> s-1]
  real :: lambda       ! Offshore decay scale, i.e. the inverse of the deformation radius of a mode [L-1 ~> m-1]
  real :: omega        ! Wave frequency [T-1 ~> s-1]
  real :: PI           ! The ratio of the circumference of a circle to its diameter [nondim]
  real :: depth_tot(SZI_(G),SZJ_(G))  ! The total depth of the ocean [Z ~> m]
  real :: depth_tot_vel ! The total depth of the ocean at a velocity point [Z ~> m]
  real :: depth_tot_corner ! The total depth of the ocean at a vorticity point [Z ~> m]
  real :: Cor_vel      ! The Coriolis parameter interpolated to a velocity point [T-1 ~> s-1]
  real    :: mag_SSH ! An overall magnitude of the external wave sea surface height at the coastline [Z ~> m]
  real    :: mag_int ! An overall magnitude of the internal wave at the coastline [L T-1 ~> m s-1]
  real    :: x1, y1  ! Various positions [L ~> m]
  real    :: x, y    ! Various positions [L ~> m]
  real    :: sin_wt  ! The sine-based periodicity factor [nondim]
  real    :: cos_wt  ! The cosine-based periodicity factor [nondim]
  real    :: val2    ! The local wave amplitude [Z ~> m]
  real    :: km_to_L_scale  ! A scaling factor from longitudes in km to L [L km-1 ~> 1e3]
  real    :: sina, cosa  ! The sine and cosine of the coast angle [nondim]
  real    :: normal_sign ! A variable that corrects the sign of normal velocities for rotation [nondim]
  real    :: trans_sign  ! A variable that corrects the sign of transverse velocities for rotation [nondim]
  type(OBC_segment_type), pointer :: segment => NULL()
  integer :: unrot_dir ! The unrotated direction of the segment
  integer :: turns    ! Number of index quarter turns
  integer :: i, j, k, n, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB, isq, ieq, jsq, jeq, is_vel, ie_vel, js_vel, je_vel

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) call MOM_error(FATAL, 'Kelvin_initialization.F90: '// &
        'Kelvin_set_OBC_data() was called but OBC type was not initialized!')
  if (G%grid_unit_to_L <= 0.) call MOM_error(FATAL, 'Kelvin_initialization.F90: '// &
          "Kelvin_set_OBC_data() is only set to work with Cartesian axis units.")

  time_sec = US%s_to_T*time_type_to_real(Time)
  PI = 4.0*atan(1.0)

  turns = modulo(G%HI%turns, 4)

  if (CS%indexing_bugs .and. (turns /= 0)) call MOM_error(FATAL, &
    "Kelvin_set_OBC_data does not support grid rotation when KELVIN_SET_OBC_INDEXING_BUGS is true.")

  do j=jsd,jed ; do i=isd,ied
    depth_tot(i,j) = 0.0
  enddo ; enddo
  do k=1,nz ; do j=jsd,jed ; do i=isd,ied
    depth_tot(i,j) = depth_tot(i,j) + GV%H_to_Z * h(i,j,k)
  enddo ; enddo ; enddo

  if (CS%mode == 0) then
    mag_SSH = CS%ssh_amp
    omega = 2.0 * PI / CS%wave_period
    sin_wt = sin(omega * time_sec)
  else
    mag_int = CS%inflow_amp
    N0 = sqrt((CS%rho_range / CS%rho_0) * (GV%g_Earth / CS%H0))
    lambda = PI * CS%mode * CS%F_0 / (CS%H0 * N0)
    ! Two wavelengths in domain
    omega = (4.0 * CS%H0 * N0)  / (CS%mode * (G%grid_unit_to_L*G%len_lon))
    ! If the modal wave speed were calculated via wave_speeds(), we should have
    !   lambda = CS%F_0 / CS%cg_mode
    !   omega = (4.0 * PI / (G%grid_unit_to_L*G%len_lon)) * CS%cg_mode
  endif
  cos_wt = cos(omega * time_sec)

  sina = sin(CS%coast_angle)
  cosa = cos(CS%coast_angle)
  do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle

   unrot_dir = segment%direction
   if (turns /= 0) unrot_dir = rotate_OBC_segment_direction(segment%direction, -turns)

    ! Apply values to the inflow end only.
    if ((unrot_dir == OBC_DIRECTION_E) .or. (unrot_dir == OBC_DIRECTION_N)) cycle

    ! Set variables that correct for sign changes during rotation.
    normal_sign = 1.0
    if ( (segment%is_E_or_W .and. ((turns == 1) .or. (turns == 2))) .or. &
         (segment%is_N_or_S .and. ((turns == 2) .or. (turns == 3))) ) normal_sign = -1.0

    ! If OBC_nudging_time is negative, the value of Velocity_nudging_timescale_in that was set
    ! when the segments are initialized is retained.
    if (CS%OBC_nudging_time >= 0.0) segment%Velocity_nudging_timescale_in = CS%OBC_nudging_time

    isd = segment%HI%isd ; ied = segment%HI%ied ; IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
    jsd = segment%HI%jsd ; jed = segment%HI%jed ; JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB

    if (unrot_dir == OBC_DIRECTION_W) then
      if (segment%is_E_or_W) then
        is_vel = IsdB ; ie_vel = IedB ; js_vel = jsd ; je_vel = jed
      else
        is_vel = isd ; ie_vel = ied ; js_vel = JsdB ; je_vel = JedB
      endif
      do j=js_vel,je_vel ; do I=is_vel,ie_vel
        if (segment%is_E_or_W) then
          x1 = G%grid_unit_to_L * G%geoLonCu(I,j)
          y1 = G%grid_unit_to_L * G%geoLatCu(I,j)
        else
          x1 = G%grid_unit_to_L * G%geoLonCv(i,J)
          y1 = G%grid_unit_to_L * G%geoLatCv(i,J)
        endif
        x = (x1 - CS%coast_offset1) * cosa + y1 * sina
        y = -(x1 - CS%coast_offset1) * sina + y1 * cosa
        if (CS%mode == 0) then
          ! Use inside bathymetry
          if (segment%direction == OBC_DIRECTION_W) then
            depth_tot_vel = depth_tot(i+1,j)
            Cor_vel = 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1))
          elseif (segment%direction == OBC_DIRECTION_S) then
            depth_tot_vel = depth_tot(i,j+1)
            Cor_vel = 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J))
          elseif (segment%direction == OBC_DIRECTION_E) then
            depth_tot_vel = depth_tot(i,j)
            Cor_vel = 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1))
          elseif (segment%direction == OBC_DIRECTION_N) then
            depth_tot_vel = depth_tot(i,j)
            Cor_vel = 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J))
          endif
          cff = sqrt(GV%g_Earth * depth_tot_vel )
          val2 = mag_SSH * exp(- Cor_vel * y / cff)
          segment%SSH(I,j) = val2 * cos_wt
          segment%normal_vel_bt(I,j) = (normal_sign*val2) * (sin_wt * cff * cosa / depth_tot_vel )
          if (segment%nudged) then
            do k=1,nz
              segment%nudged_normal_vel(I,j,k) = (normal_sign*val2) * (sin_wt * cff * cosa / depth_tot_vel )
            enddo
          elseif (segment%specified) then
            do k=1,nz
              segment%normal_vel(I,j,k) = (normal_sign*val2) * (sin_wt * cff * cosa / depth_tot_vel )
            enddo
          endif
        else
          ! Baroclinic, not rotated yet (and apparently not working as intended yet).
          segment%SSH(I,j) = 0.0
          segment%normal_vel_bt(I,j) = 0.0
          ! I suspect that the velocities in both of the following loops should instead be
          !   normal_vel(I,j,k) = CS%inflow_amp * CS%u_struct(k) * exp(-lambda * y) * cos_wt
          ! In addition, there should be a specification of the interface-height anomalies at the
          ! open boundaries that are specified as something like
          !   eta_anom(I,j,K) = (CS%inflow_amp*depth_tot/CS%cg_mode) * CS%w_struct(K) * &
          !                     exp(-lambda * y) * cos_wt
          ! In these expressions CS%u_struct and CS%w_struct could be returned from the subroutine wave_speeds
          ! in MOM_wave_speed() based on the horizontally uniform initial state.
          if (segment%nudged) then
            do k=1,nz
              segment%nudged_normal_vel(I,j,k) = (normal_sign*mag_int) * &
                   exp(-lambda * y) * cos(PI * CS%mode * (k - 0.5) / nz) * cos_wt
            enddo
          elseif (segment%specified) then
            do k=1,nz
              segment%normal_vel(I,j,k) = (normal_sign*mag_int) * &
                   exp(-lambda * y) * cos(PI * CS%mode * (k - 0.5) / nz) * cos_wt
            enddo
          endif
        endif
      enddo ; enddo
    endif

    if (unrot_dir == OBC_DIRECTION_S) then
      if (segment%is_E_or_W) then
        is_vel = IsdB ; ie_vel = IedB ; js_vel = jsd ; je_vel = jed
      else
        is_vel = isd ; ie_vel = ied ; js_vel = JsdB ; je_vel = JedB
      endif
      do J=js_vel,je_vel ; do i=is_vel,ie_vel
        if (segment%is_E_or_W) then
          x1 = G%grid_unit_to_L * G%geoLonCu(I,j)
          y1 = G%grid_unit_to_L * G%geoLatCu(I,j)
        else
          x1 = G%grid_unit_to_L * G%geoLonCv(i,J)
          y1 = G%grid_unit_to_L * G%geoLatCv(i,J)
        endif
        x = (x1 - CS%coast_offset1) * cosa + y1 * sina
        y = - (x1 - CS%coast_offset1) * sina + y1 * cosa
        if (CS%mode == 0) then
          if (segment%direction == OBC_DIRECTION_W) then
            depth_tot_vel = depth_tot(i+1,j)
            Cor_vel = 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1))
          elseif (segment%direction == OBC_DIRECTION_S) then
            depth_tot_vel = depth_tot(i,j+1)
            Cor_vel = 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J))
          elseif (segment%direction == OBC_DIRECTION_E) then
            depth_tot_vel = depth_tot(i,j)
            Cor_vel = 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1))
          elseif (segment%direction == OBC_DIRECTION_N) then
            depth_tot_vel = depth_tot(i,j)
            Cor_vel = 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J))
          endif
          cff = sqrt(GV%g_Earth * depth_tot_vel )
          val2 = mag_SSH * exp(- Cor_vel * y / cff)
          segment%SSH(I,j) = val2 * cos_wt
          segment%normal_vel_bt(I,j) = (sin_wt * cff * sina / depth_tot_vel ) * (normal_sign*val2)
          if (segment%nudged) then
            do k=1,nz
              segment%nudged_normal_vel(I,j,k) = (sin_wt * cff * sina / depth_tot_vel) * (normal_sign*val2)
            enddo
          elseif (segment%specified) then
            do k=1,nz
              segment%normal_vel(I,j,k) = (sin_wt * cff * sina / depth_tot_vel ) * (normal_sign*val2)
            enddo
          endif
        else
          ! Not rotated yet (also see the notes above on how this case might be improved)
          segment%SSH(i,J) = 0.0
          segment%normal_vel_bt(i,J) = 0.0
          if (segment%nudged) then
            do k=1,nz
              segment%nudged_normal_vel(i,J,k) = (normal_sign*mag_int) * &
                   exp(- lambda * y) * cos(PI * CS%mode * (k - 0.5) / nz) * cosa
              ! This is missing cos_wt
            enddo
          elseif (segment%specified) then
            do k=1,nz
              segment%normal_vel(i,J,k) = (normal_sign*mag_int) * &
                   exp(- lambda * y) * cos(PI * CS%mode * (k - 0.5) / nz) * cosa
              ! This is missing cos_wt
            enddo
          endif
        endif
      enddo ; enddo
    endif

    if (allocated(segment%tangential_vel)) then
      trans_sign = 1.0
      if (segment%is_E_or_W) then
        Isq = IsdB ; Ieq = IedB ; Jsq = JsdB+1 ; Jeq = JedB-1
        if ((turns == 2) .or. (turns == 3)) trans_sign = -1.0
      else
        Isq = IsdB+1 ; Ieq = IedB-1 ; Jsq = JsdB ; Jeq = JedB
        if ((turns == 1) .or. (turns == 2)) trans_sign = -1.0
      endif

      if ((unrot_dir == OBC_DIRECTION_W) .or. (unrot_dir == OBC_DIRECTION_S)) then
        do J=Jsq,Jeq ; do I=Isq,Ieq
          if (segment%direction == OBC_DIRECTION_W) then
            depth_tot_corner = 0.5*(depth_tot(i+1,j+1) + depth_tot(i+1,j))
          elseif (segment%direction == OBC_DIRECTION_E) then
            depth_tot_corner = 0.5*(depth_tot(i,j+1) + depth_tot(i,j))
          elseif (segment%direction == OBC_DIRECTION_S) then
            depth_tot_corner = 0.5*(depth_tot(i+1,j+1) + depth_tot(i,j+1))
          elseif (segment%direction == OBC_DIRECTION_N) then
            depth_tot_corner = 0.5*(depth_tot(i+1,j) + depth_tot(i,j))
          endif
          x1 = G%grid_unit_to_L * G%geoLonBu(I,J)
          y1 = G%grid_unit_to_L * G%geoLatBu(I,J)
          x = (x1 - CS%coast_offset1) * cosa + y1 * sina
          y = - (x1 - CS%coast_offset1) * sina + y1 * cosa
          cff = sqrt(GV%g_Earth * depth_tot_corner )
          val2 = (trans_sign*mag_SSH) * exp(- G%CoriolisBu(I,J) * y / cff)
          if (CS%indexing_bugs) then
            if (unrot_dir == OBC_DIRECTION_W) then
              cff = sqrt(GV%g_Earth * depth_tot(i+1,j) )
              val2 = (trans_sign*mag_SSH) * exp(- G%CoriolisBu(I,J) * y / cff)
            endif
            if (unrot_dir == OBC_DIRECTION_S) then
              cff = sqrt(GV%g_Earth * depth_tot(i,j+1) )
              val2 = (trans_sign*mag_SSH) * exp(- 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J)) * y / cff)
            endif
          endif
          if (CS%mode == 0) then ; do k=1,nz
            segment%tangential_vel(I,J,k) = (sin_wt * val2 * cff * sina) / depth_tot_corner
          enddo ; endif
        enddo ; enddo
      endif
    endif

    if (segment%specified .and. (.not.segment%nudged) .and. &
        ((unrot_dir == OBC_DIRECTION_S) .or. (unrot_dir == OBC_DIRECTION_W))) then
      if (segment%direction == OBC_DIRECTION_W) then
        do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
          segment%normal_trans(I,j,k) = segment%normal_vel(I,j,k) * h(i+1,j,k) * G%dyCu(I,j)
        enddo ; enddo ; enddo
      elseif (segment%direction == OBC_DIRECTION_E) then
        do k=1,nz ; do j=jsd,jed ; do I=IsdB,IedB
          segment%normal_trans(I,j,k) = segment%normal_vel(I,j,k) * h(i,j,k) * G%dyCu(I,j)
        enddo ; enddo ; enddo
      elseif (segment%direction == OBC_DIRECTION_S) then
        do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied
          segment%normal_trans(i,J,k) = segment%normal_vel(i,J,k) * h(i,j+1,k) * G%dxCv(i,J)
        enddo ; enddo ; enddo
      elseif (segment%direction == OBC_DIRECTION_N) then
        do k=1,nz ; do J=JsdB,JedB ; do i=isd,ied
          segment%normal_trans(i,J,k) = segment%normal_vel(i,J,k) * h(i,j,k) * G%dxCv(i,J)
        enddo ; enddo ; enddo
      endif
    endif

  enddo

end subroutine Kelvin_set_OBC_data

end module Kelvin_initialization
