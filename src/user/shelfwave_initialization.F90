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
use MOM_open_boundary,  only : OBC_registry_type
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
  real :: Lx = 100.0        !< Long-shore length scale of bathymetry.
  real :: Ly = 50.0         !< Cross-shore length scale.
  real :: f0 = 1.e-4        !< Coriolis parameter.
  real :: jj = 1            !< Cross-shore wave mode.
  real :: kk                !< Parameter.
  real :: ll                !< Longshore wavenumber.
  real :: alpha             !< 1/Ly.
  real :: omega             !< Frequency.
end type shelfwave_OBC_CS

contains

!> Add shelfwave to OBC registry.
function register_shelfwave_OBC(param_file, CS, OBC_Reg)
  type(param_file_type),    intent(in) :: param_file !< parameter file.
  type(shelfwave_OBC_CS),   pointer    :: CS         !< shelfwave control structure.
  type(OBC_registry_type),  pointer    :: OBC_Reg    !< OBC registry.
  logical                              :: register_shelfwave_OBC
  ! Local variables
  real :: kk, ll, PI, len_lat

  character(len=32)  :: casename = "shelfwave"       !< This case's name.

  PI = 4.0*atan(1.0)

  if (associated(CS)) then
    call MOM_error(WARNING, "register_shelfwave_OBC called with an "// &
                            "associated control structure.")
    return
  endif
  allocate(CS)

  ! Register the tracer for horizontal advection & diffusion.
  call register_OBC(casename, param_file, OBC_Reg)
  call get_param(param_file, mdl,"F_0",CS%f0, &
                 do_not_log=.true.)
  call get_param(param_file, mdl,"LENLAT",len_lat, &
                 do_not_log=.true.)
  call get_param(param_file, mdl,"SHELFWAVE_X_WAVELENGTH",CS%Lx, &
                 "Length scale of shelfwave in x-direction.",&
                 units="Same as x,y", default=100.)
  call get_param(param_file, mdl,"SHELFWAVE_Y_LENGTH_SCALE",CS%Ly, &
                 "Length scale of exponential dropoff of topography "//&
                 "in the y-direction.", &
                 units="Same as x,y", default=50.)
  call get_param(param_file, mdl,"SHELFWAVE_Y_MODE",CS%jj, &
                 "Cross-shore wave mode.",               &
                 units="nondim", default=1.)
  CS%alpha = 1. / CS%Ly
  CS%ll = 2. * PI / CS%Lx
  CS%kk = CS%jj * PI / len_lat
  CS%omega = 2 * CS%alpha * CS%f0 * CS%ll / &
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
                                   intent(out) :: D !< Ocean bottom depth in m or Z if US is present
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth in the units of D
  type(unit_scale_type), optional, intent(in)  :: US !< A dimensional unit scaling type

  ! Local variables
  real :: m_to_Z  ! A dimensional rescaling factor.
  integer   :: i, j
  real      :: y, rLy, Ly, H0

  m_to_Z = 1.0 ; if (present(US)) m_to_Z = US%m_to_Z

  call get_param(param_file, mdl,"SHELFWAVE_Y_LENGTH_SCALE",Ly, &
                 default=50., do_not_log=.true.)
  call get_param(param_file, mdl,"MINIMUM_DEPTH", H0, &
                 default=10., units="m", scale=m_to_Z, do_not_log=.true.)

  rLy = 0. ; if (Ly>0.) rLy = 1. / Ly

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    ! Compute normalized zonal coordinates (x,y=0 at center of domain)
    y = ( G%geoLatT(i,j) - G%south_lat )
    D(i,j) = H0 * exp(2 * rLy * y)
  enddo ; enddo

end subroutine shelfwave_initialize_topography

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine shelfwave_set_OBC_data(OBC, CS, G, GV, h, Time)
  type(ocean_OBC_type),    pointer    :: OBC !< This open boundary condition type specifies
                                             !! whether, where, and what open boundary
                                             !! conditions are used.
  type(shelfwave_OBC_CS),  pointer    :: CS  !< tidal bay control structure.
  type(ocean_grid_type),   intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV  !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< layer thickness.
  type(time_type),         intent(in) :: Time !< model time.

  ! The following variables are used to set up the transport in the shelfwave example.
  real :: my_amp, time_sec
  real :: cos_wt, cos_ky, sin_wt, sin_ky, omega, alpha
  real :: x, y, jj, kk, ll
  character(len=40)  :: mdl = "shelfwave_set_OBC_data" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, n
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment => NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return

  time_sec = time_type_to_real(Time)
  omega = CS%omega
  alpha = CS%alpha
  my_amp = 1.0
  jj = CS%jj
  kk = CS%kk
  ll = CS%ll
  do n = 1, OBC%number_of_segments
    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (segment%direction /= OBC_DIRECTION_W) cycle

    IsdB = segment%HI%IsdB ; IedB = segment%HI%IedB
    jsd = segment%HI%jsd ; jed = segment%HI%jed
    do j=jsd,jed ; do I=IsdB,IedB
      x = G%geoLonCu(I,j) - G%west_lon
      y = G%geoLatCu(I,j) - G%south_lat
      sin_wt = sin(ll*x - omega*time_sec)
      cos_wt = cos(ll*x - omega*time_sec)
      sin_ky = sin(kk * y)
      cos_ky = cos(kk * y)
      segment%normal_vel_bt(I,j) = G%US%m_s_to_L_T*my_amp * exp(- alpha * y) * cos_wt * &
           (alpha * sin_ky + kk * cos_ky)
!     segment%tangential_vel_bt(I,j) = G%US%m_s_to_L_T*my_amp * ll * exp(- alpha * y) * sin_wt * sin_ky
!     segment%vorticity_bt(I,j) = my_amp * exp(- alpha * y) * cos_wt * sin_ky&
!           (ll*ll + kk*kk + alpha*alpha)
    enddo ; enddo
  enddo

end subroutine shelfwave_set_OBC_data

end module shelfwave_initialization
