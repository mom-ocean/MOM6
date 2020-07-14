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
use MOM_open_boundary,  only : OBC_segment_type, register_OBC
use MOM_open_boundary,  only : OBC_DIRECTION_N, OBC_DIRECTION_E
use MOM_open_boundary,  only : OBC_DIRECTION_S, OBC_DIRECTION_W
use MOM_open_boundary,  only : OBC_registry_type
use MOM_unit_scaling, only : unit_scale_type
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
  real    :: coast_angle = 0   !< Angle of coastline
  real    :: coast_offset1 = 0 !< Longshore distance to coastal angle
  real    :: coast_offset2 = 0 !< Longshore distance to coastal angle
  real    :: H0 = 0            !< Bottom depth
  real    :: F_0               !< Coriolis parameter
  real    :: rho_range         !< Density range
  real    :: rho_0             !< Mean density
  logical :: answers_2018    !< If true, use the order of arithmetic and expressions that recover the
                             !! answers from the end of 2018.  Otherwise, use expressions that give
                             !! rotational symmetry and eliminate apparent bugs.
end type Kelvin_OBC_CS

! This include declares and sets the variable "version".
#include "version_variable.h"

contains

!> Add Kelvin wave to OBC registry.
function register_Kelvin_OBC(param_file, CS, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(Kelvin_OBC_CS),      pointer    :: CS         !< Kelvin wave control structure.
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.

  ! Local variables
  logical :: register_Kelvin_OBC
  logical :: default_2018_answers
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
                 default=0.0, do_not_log=.true.)
  call get_param(param_file, mdl, "TOPO_CONFIG", config, do_not_log=.true.)
  if (trim(config) == "Kelvin") then
    call get_param(param_file, mdl, "ROTATED_COAST_OFFSET_1", CS%coast_offset1, &
                   "The distance along the southern and northern boundaries "//&
                   "at which the coasts angle in.", &
                   units="km", default=100.0)
    call get_param(param_file, mdl, "ROTATED_COAST_OFFSET_2", CS%coast_offset2, &
                   "The distance from the southern and northern boundaries "//&
                   "at which the coasts angle in.", &
                   units="km", default=10.0)
    call get_param(param_file, mdl, "ROTATED_COAST_ANGLE", CS%coast_angle, &
                   "The angle of the southern bondary beyond X=ROTATED_COAST_OFFSET.", &
                   units="degrees", default=11.3)
    CS%coast_angle = CS%coast_angle * (atan(1.0)/45.) ! Convert to radians
    CS%coast_offset1 = CS%coast_offset1 * 1.e3          ! Convert to m
    CS%coast_offset2 = CS%coast_offset2 * 1.e3          ! Convert to m
  endif
  call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.true.)
  call get_param(param_file, mdl, "KELVIN_WAVE_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use expressions that give rotational "//&
                 "symmetry and eliminate apparent bugs.", default=default_2018_answers)
  if (CS%mode /= 0) then
    call get_param(param_file, mdl, "DENSITY_RANGE", CS%rho_range, &
                   default=2.0, do_not_log=.true.)
    call get_param(param_file, mdl, "RHO_0", CS%rho_0, &
                   default=1035.0, do_not_log=.true.)
    call get_param(param_file, mdl, "MAXIMUM_DEPTH", CS%H0, &
                   default=1000.0, do_not_log=.true.)
  endif

  ! Register the Kelvin open boundary.
  call register_OBC(casename, param_file, OBC_Reg)
  register_Kelvin_OBC = .true.

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
                                   intent(out) :: D !< Ocean bottom depth in m or Z if US is present
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth in the units of D
  type(unit_scale_type), optional, intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
  character(len=40)  :: mdl = "Kelvin_initialize_topography" ! This subroutine's name.
  real :: m_to_Z  ! A dimensional rescaling factor.
  real :: min_depth ! The minimum and maximum depths [Z ~> m].
  real :: PI ! 3.1415...
  real :: coast_offset1, coast_offset2, coast_angle, right_angle
  integer :: i, j

  call MOM_mesg("  Kelvin_initialization.F90, Kelvin_initialize_topography: setting topography", 5)

  m_to_Z = 1.0 ; if (present(US)) m_to_Z = US%m_to_Z

  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0, scale=m_to_Z)
  call get_param(param_file, mdl, "ROTATED_COAST_OFFSET_1", coast_offset1, &
                 default=100.0, do_not_log=.true.)
  call get_param(param_file, mdl, "ROTATED_COAST_OFFSET_2", coast_offset2, &
                 default=10.0, do_not_log=.true.)
  call get_param(param_file, mdl, "ROTATED_COAST_ANGLE", coast_angle, &
                 default=11.3, do_not_log=.true.)

  coast_angle = coast_angle * (atan(1.0)/45.) ! Convert to radians
  right_angle = 2 * atan(1.0)

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
  type(unit_scale_type),   intent(in) :: US    !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h !< layer thickness [H ~> m or kg m-2].
  type(time_type),         intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the Kelvin example.
  real :: time_sec, cff
  real :: N0           ! Brunt-Vaisala frequency [s-1]
  real :: plx          !< Longshore wave parameter
  real :: pmz          !< Vertical wave parameter
  real :: lambda       !< Offshore decay scale
  real :: omega        !< Wave frequency [s-1]
  real :: PI
  integer :: i, j, k, n, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB
  real    :: fac, x, y, x1, y1
  real    :: val1, val2, sina, cosa
  type(OBC_segment_type), pointer :: segment => NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) call MOM_error(FATAL, 'Kelvin_initialization.F90: '// &
        'Kelvin_set_OBC_data() was called but OBC type was not initialized!')

  time_sec = time_type_to_real(Time)
  PI = 4.0*atan(1.0)
  fac = 1.0

  if (CS%mode == 0) then
    omega = 2.0 * PI / (12.42 * 3600.0)      ! M2 Tide period
    val1 = US%m_to_Z * sin(omega * time_sec)
  else
    N0 = US%L_to_m*US%s_to_T * sqrt((CS%rho_range / CS%rho_0) * GV%g_Earth * (US%m_to_Z * CS%H0))
    ! Two wavelengths in domain
    plx = 4.0 * PI / G%len_lon
    pmz = PI * CS%mode / CS%H0
    lambda = pmz * CS%F_0 / N0
    omega = CS%F_0 * plx / lambda

    ! lambda = PI * CS%mode * CS%F_0 / (CS%H0 * N0)
    ! omega = (4.0 * CS%H0 * N0)  / (CS%mode * G%len_lon)
  endif

  sina = sin(CS%coast_angle)
  cosa = cos(CS%coast_angle)
  do n=1,OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle
    ! Apply values to the inflow end only.
    if (segment%direction == OBC_DIRECTION_E) cycle
    if (segment%direction == OBC_DIRECTION_N) cycle

    ! This should be somewhere else...
    segment%Velocity_nudging_timescale_in = 1.0/(0.3*86400)

    if (segment%direction == OBC_DIRECTION_W) then
      IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
      jsd = segment%HI%jsd ; jed = segment%HI%jed
      JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
      do j=jsd,jed ; do I=IsdB,IedB
        x1 = 1000. * G%geoLonCu(I,j)
        y1 = 1000. * G%geoLatCu(I,j)
        x = (x1 - CS%coast_offset1) * cosa + y1 * sina
        y = - (x1 - CS%coast_offset1) * sina + y1 * cosa
        if (CS%mode == 0) then
          cff = sqrt(GV%g_Earth * 0.5 * (G%bathyT(i+1,j) + G%bathyT(i,j)))
          val2 = fac * exp(- US%T_to_s*CS%F_0 * US%m_to_L*y / cff)
          segment%eta(I,j) = val2 * cos(omega * time_sec)
          segment%normal_vel_bt(I,j) = (val2 * (val1 * cff * cosa / &
                 (0.5 * (G%bathyT(i+1,j) + G%bathyT(i,j)))) )
          if (segment%nudged) then
            do k=1,nz
              segment%nudged_normal_vel(I,j,k) = (val2 * (val1 * cff * cosa / &
                     (0.5 * (G%bathyT(i+1,j) + G%bathyT(i,j)))) )
            enddo
          elseif (segment%specified) then
            do k=1,nz
              segment%normal_vel(I,j,k) = (val2 * (val1 * cff * cosa / &
                     (0.5 * (G%bathyT(i+1,j) + G%bathyT(i,j)))) )
              segment%normal_trans(I,j,k) = segment%normal_vel(I,j,k) * h(i+1,j,k) * G%dyCu(I,j)
            enddo
          endif
        else
          ! Not rotated yet
          segment%eta(I,j) = 0.0
          segment%normal_vel_bt(I,j) = 0.0
          if (segment%nudged) then
            do k=1,nz
              segment%nudged_normal_vel(I,j,k) = US%m_s_to_L_T * fac * lambda / CS%F_0 * &
                   exp(- lambda * y) * cos(PI * CS%mode * (k - 0.5) / nz) * &
                   cos(omega * time_sec)
            enddo
          elseif (segment%specified) then
            do k=1,nz
              segment%normal_vel(I,j,k) = US%m_s_to_L_T * fac * lambda / CS%F_0 * &
                   exp(- lambda * y) * cos(PI * CS%mode * (k - 0.5) / nz) * &
                   cos(omega * time_sec)
              segment%normal_trans(I,j,k) = segment%normal_vel(I,j,k) * h(i+1,j,k) * G%dyCu(I,j)
            enddo
          endif
        endif
      enddo ; enddo
      if (associated(segment%tangential_vel)) then
        do J=JsdB+1,JedB-1 ; do I=IsdB,IedB
          x1 = 1000. * G%geoLonBu(I,J)
          y1 = 1000. * G%geoLatBu(I,J)
          x = (x1 - CS%coast_offset1) * cosa + y1 * sina
          y = - (x1 - CS%coast_offset1) * sina + y1 * cosa
          cff =sqrt(GV%g_Earth * 0.5 * (G%bathyT(i+1,j) + G%bathyT(i,j)))
          val2 = fac * exp(- US%T_to_s*CS%F_0 * US%m_to_L*y / cff)
          if (CS%mode == 0) then ; do k=1,nz
            segment%tangential_vel(I,J,k) = (val1 * val2 * cff * sina) / &
               ( 0.25*((G%bathyT(i,j) + G%bathyT(i+1,j+1)) +  (G%bathyT(i+1,j) + G%bathyT(i,j+1))) )

          enddo ; endif
        enddo ; enddo
      endif
    else
      isd = segment%HI%isd ; ied = segment%HI%ied
      JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
      do J=JsdB,JedB ; do i=isd,ied
        x1 = 1000. * G%geoLonCv(i,J)
        y1 = 1000. * G%geoLatCv(i,J)
        x = (x1 - CS%coast_offset1) * cosa + y1 * sina
        y = - (x1 - CS%coast_offset1) * sina + y1 * cosa
        if (CS%mode == 0) then
          cff = sqrt(GV%g_Earth * 0.5 * (G%bathyT(i,j+1) + G%bathyT(i,j)))
          val2 = fac * exp(- 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J)) * US%m_to_L*y / cff)
          segment%eta(I,j) = val2 * cos(omega * time_sec)
          segment%normal_vel_bt(I,j) = US%L_T_to_m_s * (val1 * cff * sina / &
                 (0.5*(G%bathyT(i+1,j) + G%bathyT(i,j)))) * val2
          if (segment%nudged) then
            do k=1,nz
              segment%nudged_normal_vel(I,j,k) = US%L_T_to_m_s * (val1 * cff * sina / &
                     (0.5*(G%bathyT(i+1,j) + G%bathyT(i,j)))) * val2
            enddo
          elseif (segment%specified) then
            do k=1,nz
              segment%normal_vel(I,j,k) = US%L_T_to_m_s * (val1 * cff * sina / &
                     (0.5*(G%bathyT(i+1,j) + G%bathyT(i,j)))) * val2
              segment%normal_trans(i,J,k) = segment%normal_vel(i,J,k) * h(i,j+1,k) * G%dxCv(i,J)
            enddo
          endif
        else
          ! Not rotated yet
          segment%eta(i,J) = 0.0
          segment%normal_vel_bt(i,J) = 0.0
          if (segment%nudged) then
            do k=1,nz
              segment%nudged_normal_vel(i,J,k) = US%m_s_to_L_T*fac * lambda / CS%F_0 * &
                   exp(- lambda * y) * cos(PI * CS%mode * (k - 0.5) / nz) * cosa
            enddo
          elseif (segment%specified) then
            do k=1,nz
              segment%normal_vel(i,J,k) = US%m_s_to_L_T*fac * lambda / CS%F_0 * &
                   exp(- lambda * y) * cos(PI * CS%mode * (k - 0.5) / nz) * cosa
              segment%normal_trans(i,J,k) = segment%normal_vel(i,J,k) * h(i,j+1,k) * G%dxCv(i,J)
            enddo
          endif
        endif
      enddo ; enddo
      if (associated(segment%tangential_vel)) then
        do J=JsdB,JedB ; do I=IsdB+1,IedB-1
          x1 = 1000. * G%geoLonBu(I,J)
          y1 = 1000. * G%geoLatBu(I,J)
          x = (x1 - CS%coast_offset1) * cosa + y1 * sina
          y = - (x1 - CS%coast_offset1) * sina + y1 * cosa
          cff = sqrt(GV%g_Earth * 0.5 * (G%bathyT(i,j+1) + G%bathyT(i,j)))
          val2 = fac * exp(- 0.5 * (G%CoriolisBu(I,J) + G%CoriolisBu(I-1,J)) * US%m_to_L*y / cff)
          if (CS%mode == 0) then ; do k=1,nz
            segment%tangential_vel(I,J,k) = ((val1 * val2 * cff * sina) / &
                ( 0.25*((G%bathyT(i,j) + G%bathyT(i+1,j+1)) + (G%bathyT(i+1,j) +  G%bathyT(i,j+1))) ))
          enddo ; endif
        enddo ; enddo
      endif
    endif
  enddo

end subroutine Kelvin_set_OBC_data

end module Kelvin_initialization
