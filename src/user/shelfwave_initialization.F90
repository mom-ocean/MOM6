!> Configures the model for the idealized shelfwave test case.
module shelfwave_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains,        only : sum_across_PEs
use MOM_dyn_horgrid,    only : dyn_horgrid_type
use MOM_error_handler,  only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_open_boundary,  only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_W
use MOM_open_boundary,  only : OBC_segment_type, register_OBC
use MOM_open_boundary,  only : OBC_registry_type, rotate_OBC_segment_direction
use MOM_time_manager,   only : time_type, time_type_to_real
use MOM_unit_scaling,   only : unit_scale_type
use MOM_verticalGrid,   only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

character(len=40) :: mdl = "shelfwave_initialization" !< This module's name.

! The following routines are visible to the outside world
public shelfwave_initialize_topography
public shelfwave_set_OBC_data
public register_shelfwave_OBC, shelfwave_OBC_end

!> Control structure for shelfwave open boundaries.
type, public :: shelfwave_OBC_CS ; private
  real :: my_amp        !< Amplitude of the open boundary current inflows [L T-1 ~> m s-1]
  real :: kk            !< Cross-shore wavenumber [km-1] or [m-1]
  real :: ll            !< Longshore wavenumber [km-1] or [m-1]
  real :: alpha         !< Exponential decay rate in the y-direction [km-1] or [m-1]
  real :: omega         !< Frequency of the shelf wave [T-1 ~> s-1]
  logical :: shelfwave_correct_amplitude !< If true, SHELFWAVE_AMPLITUDE gives the actual inflow
                        !! velocity, rather than giving an overall scaling factor for the flow.
end type shelfwave_OBC_CS

contains

!> Add shelfwave to OBC registry.
function register_shelfwave_OBC(param_file, CS, G, US, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(shelfwave_OBC_CS),   pointer    :: CS         !< shelfwave control structure.
  type(ocean_grid_type),    intent(in) :: G          !< The ocean's grid structure.
  type(unit_scale_type),    intent(in) :: US         !< A dimensional unit scaling type
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< Open boundary condition registry.
  logical                              :: register_shelfwave_OBC

  ! Local variables
  real :: PI      ! The ratio of the circumference of a circle to its diameter [nondim]
  character(len=32)  :: casename = "shelfwave"       !< This case's name.
  real :: jj      ! Cross-shore wave mode [nondim]
  real :: f0      ! Coriolis parameter [T-1 ~> s-1]
  real :: Lx      ! Long-shore length scale of bathymetry [km] or [m]
  real :: Ly      ! Cross-shore length scale [km] or [m]
  real :: default_amp  ! The default velocity amplitude [m s-1] or amplitude scaling factor [nondim]

  PI = 4.0*atan(1.0)

  if (associated(CS)) then
    call MOM_error(WARNING, "register_shelfwave_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Register the tracer for horizontal advection & diffusion.
  call register_OBC(casename, param_file, OBC_Reg)
  call get_param(param_file, mdl, "F_0", f0, &
                 default=0.0, units="s-1", scale=US%T_to_s, do_not_log=.true.)
  call get_param(param_file, mdl,"SHELFWAVE_X_WAVELENGTH", Lx, &
                 "Length scale of shelfwave in x-direction.",&
                 units=G%x_ax_unit_short, default=100.)
  call get_param(param_file, mdl, "SHELFWAVE_Y_LENGTH_SCALE", Ly, &
                 "Length scale of exponential dropoff of topography in the y-direction.", &
                 units=G%y_ax_unit_short, default=50.)
  call get_param(param_file, mdl, "SHELFWAVE_Y_MODE", jj, &
                 "Cross-shore wave mode.",               &
                 units="nondim", default=1.)
  call get_param(param_file, mdl, "SHELFWAVE_CORRECT_AMPLITUDE", CS%shelfwave_correct_amplitude, &
                 "If true, SHELFWAVE_AMPLITUDE gives the actual inflow velocity, rather than giving "//&
                 "an overall scaling factor for the flow.", default=.false.)  !### Make the default .true.?
  default_amp = 1.0 ; if (CS%shelfwave_correct_amplitude) default_amp = 0.1
  call get_param(param_file, mdl, "SHELFWAVE_AMPLITUDE", CS%my_amp, &
                 "Amplitude of the open boundary current inflows in the shelfwave configuration.", &
                 units="m s-1", default=default_amp, scale=US%m_s_to_L_T)

  CS%alpha = 1. / Ly
  CS%ll = 2. * PI / Lx
  CS%kk = jj * PI / G%len_lat
  CS%omega = 2 * CS%alpha * f0 * CS%ll / &
             (CS%kk*CS%kk + CS%alpha*CS%alpha + CS%ll*CS%ll)
  register_shelfwave_OBC = .true.

end function register_shelfwave_OBC

!> Clean up the shelfwave OBC from registry.
subroutine shelfwave_OBC_end(CS)
  type(shelfwave_OBC_CS), pointer    :: CS         !< shelfwave control structure.

  if (associated(CS)) then
    deallocate(CS)
  endif
end subroutine shelfwave_OBC_end

!> Initialization of topography.
subroutine shelfwave_initialize_topography( D, G, param_file, max_depth, US )
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
  real      :: y    ! Position relative to the southern boundary [km] or [m] or [degrees_N]
  real      :: rLy  ! Exponential decay rate of the topography [km-1] or [m-1] or [degrees_N-1]
  real      :: Ly   ! Exponential decay lengthscale of the topography [km] or [m] or [degrees_N]
  real      :: H0   ! The minimum depth of the ocean [Z ~> m]
  integer   :: i, j

  call get_param(param_file, mdl,"SHELFWAVE_Y_LENGTH_SCALE", Ly, &
                 units=G%y_ax_unit_short, default=50., do_not_log=.true.)
  call get_param(param_file, mdl,"MINIMUM_DEPTH", H0, &
                 units="m", default=10., scale=US%m_to_Z, do_not_log=.true.)

  rLy = 0. ; if (Ly>0.) rLy = 1. / Ly

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! Compute normalized zonal coordinates (x,y=0 at center of domain)
    y = ( G%geoLatT(i,j) - G%south_lat )
    D(i,j) = H0 * exp(2 * rLy * y)
  enddo ; enddo

end subroutine shelfwave_initialize_topography

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine shelfwave_set_OBC_data(OBC, CS, G, GV, US, h, Time)
  type(ocean_OBC_type),    pointer    :: OBC  !< This open boundary condition type specifies
                                              !! whether, where, and what open boundary
                                              !! conditions are used.
  type(shelfwave_OBC_CS),  pointer    :: CS   !< tidal bay control structure.
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness [H ~> m or kg m-2]
  type(time_type),         intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the shelfwave example.
  real :: time_sec ! The time in the run [T ~> s]
  real :: cos_wt, sin_wt ! Cosine and sine associated with the propagating x-direction structure [nondim]
  real :: cos_ky, sin_ky ! Cosine and sine associated with the y-direction structure [nondim]
  real :: x   ! Position relative to the western boundary [km] or [m] or [degrees_E]
  real :: y   ! Position relative to the southern boundary [km] or [m] or [degrees_N]
  real :: I_yscale  ! A factor to give the correct inflow velocity [km-1] or [m-1] or [degrees_N-1] or
                    ! to compensate for the variable units of the y-coordinate [km axis_unit-1], usually 1 [nondim]
  real :: my_amp    ! Amplitude of the open boundary current inflows, including sign changes
                    ! to account for grid rotation [L T-1 ~> m s-1]
  integer :: i, j, is, ie, js, je, n
  integer :: turns    ! Number of index quarter turns
  type(OBC_segment_type), pointer :: segment => NULL()

  if (.not.associated(OBC)) return

  turns = modulo(G%HI%turns, 4)
  my_amp = CS%my_amp ; if ((turns==2) .or. (turns==3)) my_amp = -CS%my_amp

  time_sec = US%s_to_T*time_type_to_real(Time)
  if (CS%shelfwave_correct_amplitude) then
    ! This makes the units and edge value of normal_vel_bt the same as my_amp.
    I_yscale = 1.0 / CS%kk
  else ! This preserves the previous answers.
    if (G%grid_unit_to_L == 0.0) call MOM_error(FATAL, &
          "shelfwave_set_OBC_data requires the use of Cartesian coordinates.")
    I_yscale = (1.0e3 * US%m_to_L) / G%grid_unit_to_L
  endif
  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (rotate_OBC_segment_direction(segment%direction, -turns) /= OBC_DIRECTION_W) cycle

    if (segment%is_E_or_W) then
      ! segment thicknesses are defined at cell face centers.
      is = segment%HI%isdB ; ie = segment%HI%iedB
      js = segment%HI%jsd ; je = segment%HI%jed
    else
      is = segment%HI%isd ; ie = segment%HI%ied
      js = segment%HI%jsdB ; je = segment%HI%jedB
    endif

    do j=js,je ; do I=is,ie
      if (segment%is_E_or_W) then
        x = G%geoLonCu(I,j) - G%west_lon
        y = G%geoLatCu(I,j) - G%south_lat
      else
        x = G%geoLonCv(i,J) - G%west_lon
        y = G%geoLatCv(i,J) - G%south_lat
      endif
      sin_wt = sin(CS%ll*x - CS%omega*time_sec)
      cos_wt = cos(CS%ll*x - CS%omega*time_sec)
      sin_ky = sin(CS%kk * y)
      cos_ky = cos(CS%kk * y)
      segment%normal_vel_bt(I,j) = my_amp * exp(- CS%alpha * y) * cos_wt * &
           ((CS%alpha * sin_ky + CS%kk * cos_ky) * I_yscale)
!     segment%tangential_vel_bt(I,j) = my_amp * (CS%ll * I_yscale) * exp(- CS%alpha * y) * sin_wt * sin_ky
!     segment%vorticity_bt(I,j) = my_amp * exp(- CS%alpha * y) * cos_wt * sin_ky * &
!           ((CS%ll**2 + CS%kk**2 + CS%alpha**2) * (I_yscale / G%grid_unit_to_L))
    enddo ; enddo
  enddo

end subroutine shelfwave_set_OBC_data

end module shelfwave_initialization
