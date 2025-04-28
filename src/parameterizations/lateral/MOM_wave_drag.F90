!> Frequency-dependent linear wave drag

module MOM_wave_drag

use MOM_domains,       only : pass_vector, To_All, Scalar_Pair
use MOM_error_handler, only : MOM_error, NOTE
use MOM_file_parser,   only : get_param, log_param, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_io,            only : MOM_read_data, slasher, EAST_FACE, NORTH_FACE
use MOM_unit_scaling,  only : unit_scale_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

public wave_drag_init, wave_drag_calc

#include <MOM_memory.h>

!> Control structure for the MOM_wave_drag module
type, public :: wave_drag_CS ; private
  integer :: nf                                 !< Number of filters to be used in the simulation
  real, allocatable, dimension(:,:,:) :: coef_u !< frequency-dependent drag coefficients [H T-1 ~> m s-1]
  real, allocatable, dimension(:,:,:) :: coef_v !< frequency-dependent drag coefficients [H T-1 ~> m s-1]
end type wave_drag_CS

contains

!> This subroutine reads drag coefficients from file.
subroutine wave_drag_init(param_file, wave_drag_file, G, GV, US, CS)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  character(len=*),        intent(in)    :: wave_drag_file !< The file from which to read drag coefficients
  type(ocean_grid_type),   intent(inout) :: G          !< The ocean's grid structure
  type(verticalGrid_type), intent(in)    :: GV         !< The ocean's vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(wave_drag_CS),      intent(out)   :: CS         !< Control structure of MOM_wave_drag

  ! Local variables
  character(len=40)  :: mdl = "MOM_wave_drag"          !< This module's name
  character(len=50)  :: filter_name_str                !< List of drag coefficients to be used
  character(len=2),  allocatable, dimension(:) :: filter_names !< Names of drag coefficients
  character(len=80)  :: var_names(2)                   !< Names of variables in wave_drag_file
  character(len=200) :: mesg
  real               :: var_scale                      !< Scaling factors of drag coefficients [nondim]
  integer            :: c

  ! The number and names of drag coefficients should match those of the streaming filters.
  call get_param(param_file, mdl, "N_FILTERS", CS%nf, &
                 "Number of streaming band-pass filters to be used in the simulation.", &
                 default=0, do_not_log=.true.)
  call get_param(param_file, mdl, "FILTER_NAMES", filter_name_str, &
                 "Names of streaming band-pass filters to be used in the simulation.", &
                 do_not_log=.true.)

  allocate(CS%coef_u(G%IsdB:G%IedB,G%jsd:G%jed,CS%nf)) ; CS%coef_u(:,:,:) = 0.0
  allocate(CS%coef_v(G%isd:G%ied,G%JsdB:G%JedB,CS%nf)) ; CS%coef_v(:,:,:) = 0.0
  allocate(filter_names(CS%nf)) ; read(filter_name_str, *) filter_names

  if (len_trim(wave_drag_file) > 0) then
    do c=1,CS%nf
      call get_param(param_file, mdl, "BT_"//trim(filter_names(c))//"_DRAG_U", &
                     var_names(1), "The name of the variable in BT_WAVE_DRAG_FILE "//&
                     "for the drag coefficient of the "//trim(filter_names(c))//&
                     " frequency at u points.", default="")
      call get_param(param_file, mdl, "BT_"//trim(filter_names(c))//"_DRAG_V", &
                     var_names(2), "The name of the variable in BT_WAVE_DRAG_FILE "//&
                     "for the drag coefficient of the "//trim(filter_names(c))//&
                     " frequency at v points.", default="")
      call get_param(param_file, mdl, "BT_"//trim(filter_names(c))//"_DRAG_SCALE", &
                     var_scale, "A scaling factor for the drag coefficient of the "//&
                     trim(filter_names(c))//" frequency.", default=1.0, units="nondim")

      if (len_trim(var_names(1))+len_trim(var_names(2))>0 .and. var_scale>0.0) then
        call MOM_read_data(wave_drag_file, trim(var_names(1)), CS%coef_u(:,:,c), G%Domain, &
                           position=EAST_FACE, scale=var_scale*GV%m_to_H*US%T_to_s)
        call MOM_read_data(wave_drag_file, trim(var_names(2)), CS%coef_v(:,:,c), G%Domain, &
                           position=NORTH_FACE, scale=var_scale*GV%m_to_H*US%T_to_s)
        call pass_vector(CS%coef_u(:,:,c), CS%coef_v(:,:,c), G%domain, &
                         direction=To_All+SCALAR_PAIR)

        write(mesg, *) "MOM_wave_drag: ", trim(filter_names(c)), &
                       " coefficients read from file, scaling factor = ", var_scale
        call MOM_error(NOTE, trim(mesg))
      endif ! (len_trim(var_names(1))+len_trim(var_names(2))>0 .and. var_scale>0.0)
    enddo ! k=1,CS%nf
  endif ! (len_trim(wave_drag_file) > 0)

end subroutine wave_drag_init

!> This subroutine calculates the sum of the products of the tidal velocities and the scaled
!! frequency-dependent drag for each tidal constituent specified in MOM_input.
subroutine wave_drag_calc(u, v, drag_u, drag_v, G, CS)
  type(ocean_grid_type),           intent(in) :: G     !< The ocean's grid structure
  type(wave_drag_CS),              intent(in) :: CS    !< Control structure of MOM_wave_drag
  real, dimension(:,:,:), pointer, intent(in) :: u     !< Zonal velocity from the output of
                                                       !! streaming band-pass filters [L T-1 ~> m s-1]
  real, dimension(:,:,:), pointer, intent(in) :: v     !< Meridional velocity from the output of
                                                       !! streaming band-pass filters [L T-1 ~> m s-1]
  real, dimension(G%IsdB:G%IedB,G%jsd:G%jed), intent(out) :: drag_u !< Sum of products of filtered velocities
                                                       !! and scaled frequency-dependent drag [L2 T-2 ~> m2 s-2]
  real, dimension(G%isd:G%ied,G%JsdB:G%JedB), intent(out) :: drag_v !< Sum of products of filtered velocities
                                                       !! and scaled frequency-dependent drag [L2 T-2 ~> m2 s-2]

  ! Local variables
  integer :: is, ie, js, je, i, j, k

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  Drag_u(:,:) = 0.0 ; Drag_v(:,:) = 0.0

  !$OMP do
  do k=1,CS%nf ; do j=js,je ; do I=is-1,ie
    Drag_u(I,j) = Drag_u(I,j) + u(I,j,k) * CS%coef_u(I,j,k)
  enddo ; enddo ; enddo

  !$OMP do
  do k=1,CS%nf ; do J=js-1,je ; do i=is,ie
    Drag_v(i,J) = Drag_v(i,J) + v(i,J,k) * CS%coef_v(i,J,k)
  enddo ; enddo ; enddo

end subroutine wave_drag_calc

!> \namespace mom_wave_drag
!!
!! By Chengzhu Xu (chengzhu.xu@oregonstate.edu) and Edward D. Zaron, December 2024
!!
!! This module calculates the net effects of the frequency-dependent internal wave drag applied to
!! the tidal velocities, and returns the sum of products of frequency-dependent drag coefficients
!! and tidal velocities for each constituent to the MOM_barotropic module for further calculations.
!! It relies on the use of MOM_streaming_filter for determining the tidal velocities. Furthermore,
!! the number of drag coefficients cannot exceed that of the streaming filters, and the names of
!! drag coefficients should match those of the streaming filters. The frequency-dependent drag
!! coefficients are read from the same file for the linear drag coefficients in MOM_barotropic.

end module MOM_wave_drag

