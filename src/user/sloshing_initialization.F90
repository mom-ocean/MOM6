!> Initialization for the "sloshing" internal waves configuration.
module sloshing_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_domains, only : sum_across_PEs
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type

implicit none ; private

#include <MOM_memory.h>

! The following routines are visible to the outside world
public sloshing_initialize_topography
public sloshing_initialize_thickness
public sloshing_initialize_temperature_salinity

contains

!> Initialization of topography.
subroutine sloshing_initialize_topography( D, G, param_file, max_depth )
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth in the units of depth_max
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real,                    intent(in)  :: max_depth !< Maximum ocean depth in arbitrary units

  ! Local variables
  integer   :: i, j

  do i=G%isc,G%iec ; do j=G%jsc,G%jec
    D(i,j) = max_depth
  enddo ; enddo

end subroutine sloshing_initialize_topography


!> Initialization of thicknesses
!! This routine is called when THICKNESS_CONFIG is set to 'sloshing'
!!
!! This routine initializes layer positions to set off a sloshing motion in
!! the zonal direction in a rectangular basin. All layers have initially the
!! same thickness but all interfaces (except bottom and sea surface) are
!! displaced according to a half-period cosine, with maximum value on the
!! left and minimum value on the right. This sets off a regular sloshing motion.
subroutine sloshing_initialize_thickness ( h, G, GV, US, param_file, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2].
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.

  real    :: displ(SZK_(GV)+1)  ! The interface displacement [Z ~> m].
  real    :: z_unif(SZK_(GV)+1) ! Fractional uniform interface heights [nondim].
  real    :: z_inter(SZK_(GV)+1) ! Interface heights [Z ~> m]
  real    :: a0                 ! The displacement amplitude [Z ~> m].
  real    :: weight_z           ! A (misused?) depth-space weighting, in inconsistent units.
  real    :: x1, y1, x2, y2     ! Dimensonless parameters.
  real    :: x, t               ! Dimensionless depth coordinates?
  logical :: use_IC_bug         ! If true, set the initial conditions retaining an old bug.
  logical :: just_read          ! If true, just read parameters but set nothing.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "sloshing_initialization" !< This module's name.

  integer :: i, j, k, is, ie, js, je, nx, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params
  if (.not.just_read) call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "SLOSHING_IC_AMPLITUDE", a0, &
                 "Initial amplitude of sloshing internal interface height "//&
                 "displacements it the sloshing test case.", &
                 units='m', default=75.0, scale=US%m_to_Z, do_not_log=just_read)
  call get_param(param_file, mdl, "SLOSHING_IC_BUG", use_IC_bug, &
                 "If true, use code with a bug to set the sloshing initial conditions.", &
                 default=.false., do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  ! Define thicknesses
  do j=G%jsc,G%jec ; do i=G%isc,G%iec

    ! Define uniform interfaces
    do k = 0,nz
      z_unif(k+1) = -real(k)/real(nz)
    enddo

    ! 1. Define stratification
    do k = 1,nz+1

      ! Thin pycnocline in the middle
      !z_inter(k) = (2.0**(n-1)) * (z_unif(k) + 0.5)**n - 0.5

      ! Thin pycnocline in the middle (piecewise linear profile)
      x1 = 0.30; y1 = 0.48; x2 = 0.70; y2 = 0.52

      x = -z_unif(k)

      if ( x <= x1 ) then
        t = y1*x/x1
      elseif ( (x > x1 ) .and. ( x < x2 )) then
        t = y1 + (y2-y1) * (x-x1) / (x2-x1)
      else
        t = y2 + (1.0-y2) * (x-x2) / (1.0-x2)
      endif

      t = - z_unif(k)

      z_inter(k) = -t * G%max_depth

    enddo

    ! 2. Define displacement
    ! a0 is set via get_param; by default a0 is a 75m Displacement amplitude in depth units.
    do k = 1,nz+1

      weight_z = - 4.0 * ( z_unif(k) + 0.5 )**2 + 1.0

      x = G%geoLonT(i,j) / G%len_lon
      if (use_IC_bug) then
        displ(k) = a0 * cos(acos(-1.0)*x) + weight_z * US%m_to_Z
      else
        displ(k) = a0 * cos(acos(-1.0)*x) * weight_z
      endif

      if ( k == 1 ) then
        displ(k) = 0.0
      endif

      if ( k == nz+1 ) then
        displ(k) = 0.0
      endif

      z_inter(k) = z_inter(k) + displ(k)

    enddo

    ! 3. The last interface must coincide with the seabed
    z_inter(nz+1) = -G%bathyT(i,j)
    ! Modify interface heights to make sure all thicknesses are strictly positive
    do k = nz,1,-1
      if ( z_inter(k) < (z_inter(k+1) + GV%Angstrom_Z) ) then
        z_inter(k) = z_inter(k+1) + GV%Angstrom_Z
      endif
    enddo

    ! 4. Define layers
    do k = 1,nz
      h(i,j,k) = GV%Z_to_H * (z_inter(k) - z_inter(k+1))
    enddo

  enddo ; enddo

end subroutine sloshing_initialize_thickness


!> Initialization of temperature and salinity
!!
!! This subroutine initializes linear profiles for T and S according to
!! reference surface layer salinity and temperature and a specified range.
!! Note that the linear distribution is set up with respect to the layer
!! number, not the physical position).
subroutine sloshing_initialize_temperature_salinity ( T, S, h, G, GV, param_file, &
                                                      eqn_of_state, just_read_params)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure.
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T !< Potential temperature [degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S !< Salinity [ppt].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h !< Layer thickness [H ~> m or kg m-2].
  type(param_file_type),                     intent(in)  :: param_file !< A structure indicating the
                                                            !! open file to parse for model
                                                            !! parameter values.
  type(EOS_type),                            pointer     :: eqn_of_state !< Equation of state structure.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.

  integer :: i, j, k, is, ie, js, je, nz
  real    :: delta_S, delta_T
  real    :: S_ref, T_ref;      ! Reference salinity and temerature within
                                ! surface layer
  real    :: S_range, T_range;  ! Range of salinities and temperatures over the
                                ! vertical
  integer :: kdelta
  real    :: deltah
  real    :: xi0, xi1
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "initialize_temp_salt_linear" ! This subroutine's
                                                            ! name.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(param_file, mdl, "S_REF", S_ref, 'Reference value for salinity', &
                 default=35.0, units='1e-3', do_not_log=just_read)
  call get_param(param_file, mdl, "T_REF", T_ref, 'Reference value for temperature', &
                 units='degC', fail_if_missing=.not.just_read, do_not_log=just_read)

  ! The default is to assume an increase by 2 ppt for the salinity and a uniform temperature.
  call get_param(param_file, mdl,"S_RANGE",S_range,'Initial salinity range.', &
                 units='1e-3', default=2.0, do_not_log=just_read)
  call get_param(param_file, mdl,"T_RANGE",T_range,'Initial temperature range', &
                 units='degC', default=0.0, do_not_log=just_read)

  if (just_read) return ! All run-time parameters have been read, so return.

  ! Prescribe salinity
  !delta_S = S_range / ( GV%ke - 1.0 )

  !S(:,:,1) = S_ref
  !do k = 2,GV%ke
  !  S(:,:,k) = S(:,:,k-1) + delta_S
  !enddo

  deltah = G%max_depth / nz
  do j=js,je ; do i=is,ie
    xi0 = 0.0
    do k = 1,nz
      xi1 = xi0 + deltah / G%max_depth ! =  xi0 + 1.0 / real(nz)
      S(i,j,k) = 34.0 + 0.5 * S_range * (xi0 + xi1)
      xi0 = xi1
    enddo
  enddo ; enddo

  ! Prescribe temperature
  delta_T = T_range / ( GV%ke - 1.0 )

  T(:,:,1) = T_ref
  do k = 2,GV%ke
    T(:,:,k) = T(:,:,k-1) + delta_T
  enddo
  kdelta = 2
  T(:,:,GV%ke/2 - (kdelta-1):GV%ke/2 + kdelta) = 1.0

end subroutine sloshing_initialize_temperature_salinity

!> \namespace sloshing_initialization
!!
!! The module configures the model for the non-rotating sloshing test case.
end module sloshing_initialization
