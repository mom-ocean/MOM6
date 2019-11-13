!> Interface to CVMix convection scheme.
module MOM_CVMix_conv

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,      only : hchksum
use MOM_diag_mediator,  only : diag_ctrl, time_type, register_diag_field
use MOM_diag_mediator,  only : post_data
use MOM_EOS,            only : calculate_density
use MOM_error_handler,  only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser,    only : openParameterBlock, closeParameterBlock
use MOM_file_parser,    only : get_param, log_version, param_file_type
use MOM_grid,           only : ocean_grid_type
use MOM_unit_scaling,   only : unit_scale_type
use MOM_variables,      only : thermo_var_ptrs
use MOM_verticalGrid,   only : verticalGrid_type
use CVMix_convection,   only : CVMix_init_conv, CVMix_coeffs_conv
use CVMix_kpp,          only : CVMix_kpp_compute_kOBL_depth

implicit none ; private

#include <MOM_memory.h>

public CVMix_conv_init, calculate_CVMix_conv, CVMix_conv_end, CVMix_conv_is_used

!> Control structure including parameters for CVMix convection.
type, public :: CVMix_conv_cs

  ! Parameters
  real    :: kd_conv_const !< diffusivity constant used in convective regime [m2 s-1]
  real    :: kv_conv_const !< viscosity constant used in convective regime [m2 s-1]
  real    :: bv_sqr_conv   !< Threshold for squared buoyancy frequency
                           !! needed to trigger Brunt-Vaisala parameterization [s-2]
  real    :: min_thickness !< Minimum thickness allowed [m]
  logical :: debug         !< If true, turn on debugging

  ! Daignostic handles and pointers
  type(diag_ctrl), pointer :: diag => NULL() !< Pointer to diagnostics control structure
  !>@{ Diagnostics handles
  integer :: id_N2 = -1, id_kd_conv = -1, id_kv_conv = -1
  !!@}

  ! Diagnostics arrays
  real, allocatable, dimension(:,:,:) :: N2      !< Squared Brunt-Vaisala frequency [s-2]
  real, allocatable, dimension(:,:,:) :: kd_conv !< Diffusivity added by convection [Z2 T-1 ~> m2 s-1]
  real, allocatable, dimension(:,:,:) :: kv_conv !< Viscosity added by convection [Z2 T-1 ~> m2 s-1]

end type CVMix_conv_cs

character(len=40)  :: mdl = "MOM_CVMix_conv"     !< This module's name.

contains

!> Initialized the CVMix convection mixing routine.
logical function CVMix_conv_init(Time, G, GV, US, param_file, diag, CS)

  type(time_type),         intent(in)    :: Time       !< The current time.
  type(ocean_grid_type),   intent(in)    :: G          !< Grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Run-time parameter file handle
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics control structure.
  type(CVMix_conv_cs),     pointer       :: CS         !< This module's control structure.
  ! Local variables
  real    :: prandtl_conv !< Turbulent Prandtl number used in convective instabilities.
  logical :: useEPBL      !< If True, use the ePBL boundary layer scheme.

! This include declares and sets the variable "version".
#include "version_variable.h"

  if (associated(CS)) then
    call MOM_error(WARNING, "CVMix_conv_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  ! Read parameters
  call log_version(param_file, mdl, version, &
    "Parameterization of enhanced mixing due to convection via CVMix")
  call get_param(param_file, mdl, "USE_CVMix_CONVECTION", CVMix_conv_init, &
                 "If true, turns on the enhanced mixing due to convection "//&
                 "via CVMix. This scheme increases diapycnal diffs./viscs. "//&
                 "at statically unstable interfaces. Relevant parameters are "//&
                 "contained in the CVMix_CONVECTION% parameter block.", &
                 default=.false.)

  if (.not. CVMix_conv_init) return

  call get_param(param_file, mdl, "ENERGETICS_SFC_PBL", useEPBL, default=.false., &
                do_not_log=.true.)

  ! Warn user if EPBL is being used, since in this case mixing due to convection will
  ! be aplied in the boundary layer
  if (useEPBL) then
     call MOM_error(WARNING, 'MOM_CVMix_conv_init: '// &
           'CVMix convection may not be properly applied when ENERGETICS_SFC_PBL = True'//&
           'as convective mixing might occur in the boundary layer.')
  endif

  call get_param(param_file, mdl, 'DEBUG', CS%debug, default=.False., do_not_log=.True.)

  call get_param(param_file, mdl, 'MIN_THICKNESS', CS%min_thickness, default=0.001, do_not_log=.True.)

  call openParameterBlock(param_file,'CVMix_CONVECTION')

  call get_param(param_file, mdl, "PRANDTL_CONV", prandtl_conv, &
                 "The turbulent Prandtl number applied to convective "//&
                 "instabilities (i.e., used to convert KD_CONV into KV_CONV)", &
                 units="nondim", default=1.0)

  call get_param(param_file, mdl, 'KD_CONV', CS%kd_conv_const, &
                 "Diffusivity used in convective regime. Corresponding viscosity "//&
                 "(KV_CONV) will be set to KD_CONV * PRANDTL_TURB.", &
                 units='m2/s', default=1.00)

  call get_param(param_file, mdl, 'BV_SQR_CONV', CS%bv_sqr_conv, &
                 "Threshold for squared buoyancy frequency needed to trigger "//&
                 "Brunt-Vaisala parameterization.", &
                 units='1/s^2', default=0.0)

  call closeParameterBlock(param_file)

  ! set kv_conv_const based on kd_conv_const and prandtl_conv
  CS%kv_conv_const = CS%kd_conv_const * prandtl_conv

  ! allocate arrays and set them to zero
  allocate(CS%N2(SZI_(G), SZJ_(G), SZK_(G)+1)); CS%N2(:,:,:) = 0.
  allocate(CS%kd_conv(SZI_(G), SZJ_(G), SZK_(G)+1)); CS%kd_conv(:,:,:) = 0.
  allocate(CS%kv_conv(SZI_(G), SZJ_(G), SZK_(G)+1)); CS%kv_conv(:,:,:) = 0.

  ! Register diagnostics
  CS%diag => diag
  CS%id_N2 = register_diag_field('ocean_model', 'N2_conv', diag%axesTi, Time, &
      'Square of Brunt-Vaisala frequency used by MOM_CVMix_conv module', '1/s2')
  CS%id_kd_conv = register_diag_field('ocean_model', 'kd_conv', diag%axesTi, Time, &
      'Additional diffusivity added by MOM_CVMix_conv module', 'm2/s', conversion=US%Z2_T_to_m2_s)
  CS%id_kv_conv = register_diag_field('ocean_model', 'kv_conv', diag%axesTi, Time, &
      'Additional viscosity added by MOM_CVMix_conv module', 'm2/s', conversion=US%Z2_T_to_m2_s)

  call CVMix_init_conv(convect_diff=CS%kd_conv_const, &
                       convect_visc=CS%kv_conv_const, &
                       lBruntVaisala=.true.,    &
                       BVsqr_convect=CS%bv_sqr_conv)

end function CVMix_conv_init

!> Subroutine for calculating enhanced diffusivity/viscosity
!! due to convection via CVMix
subroutine calculate_CVMix_conv(h, tv, G, GV, US, CS, hbl)

  type(ocean_grid_type),                      intent(in)  :: G  !< Grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV !< Vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: h  !< Layer thickness [H ~> m or kg m-2].
  type(thermo_var_ptrs),                      intent(in)  :: tv !< Thermodynamics structure.
  type(CVMix_conv_cs),                            pointer :: CS !< The control structure returned
                                                                !! by a previous call to CVMix_conv_init.
  real, dimension(:,:),                 optional, pointer :: hbl!< Depth of ocean boundary layer [m]
  ! local variables
  real, dimension(SZK_(G)) :: rho_lwr !< Adiabatic Water Density, this is a dummy
                                      !! variable since here convection is always
                                      !! computed based on Brunt Vaisala.
  real, dimension(SZK_(G)) :: rho_1d  !< water density in a column, this is also
                                      !! a dummy variable, same reason as above.
  real, dimension(SZK_(G)+1) :: kv_col !< Viscosities at interfaces in the column [m2 s-1]
  real, dimension(SZK_(G)+1) :: kd_col !< Diffusivities at interfaces in the column [m2 s-1]
  real, dimension(SZK_(G)+1) :: iFaceHeight !< Height of interfaces [m]
  real, dimension(SZK_(G))   :: cellHeight  !< Height of cell centers [m]
  integer :: kOBL                        !< level of OBL extent
  real :: g_o_rho0  ! Gravitational acceleration divided by density in MKS units [m4 s-2]
  real :: pref, rhok, rhokm1, dz, dh, hcorr
  integer :: i, j, k

  g_o_rho0 = GV%mks_g_Earth / (US%R_to_kg_m3*GV%Rho0)

  ! initialize dummy variables
  rho_lwr(:) = 0.0; rho_1d(:) = 0.0

  if (.not. associated(hbl)) then
    allocate(hbl(SZI_(G), SZJ_(G)))
    hbl(:,:) = 0.0
  endif

  do j = G%jsc, G%jec
    do i = G%isc, G%iec

      ! set N2 to zero at the top- and bottom-most interfaces
      CS%N2(i,j,1) = 0.
      CS%N2(i,j,G%ke+1) = 0.

      ! skip calling at land points
      !if (G%mask2dT(i,j) == 0.) cycle

      pRef = 0.
      ! Compute Brunt-Vaisala frequency (static stability) on interfaces
      do k=2,G%ke

        ! pRef is pressure at interface between k and km1.
        pRef = pRef + GV%H_to_Pa * h(i,j,k)
        call calculate_density(tv%t(i,j,k), tv%s(i,j,k), pref, rhok, tv%eqn_of_state)
        call calculate_density(tv%t(i,j,k-1), tv%s(i,j,k-1), pref, rhokm1, tv%eqn_of_state)

        dz = ((0.5*(h(i,j,k-1) + h(i,j,k))+GV%H_subroundoff)*GV%H_to_m)
        CS%N2(i,j,k) = g_o_rho0 * (rhok - rhokm1) / dz ! Can be negative

      enddo

      iFaceHeight(1) = 0.0 ! BBL is all relative to the surface
      hcorr = 0.0
      ! compute heights at cell center and interfaces
      do k=1,G%ke
        dh = h(i,j,k) * GV%H_to_m ! Nominal thickness to use for increment
        dh = dh + hcorr ! Take away the accumulated error (could temporarily make dh<0)
        hcorr = min( dh - CS%min_thickness, 0. ) ! If inflating then hcorr<0
        dh = max( dh, CS%min_thickness ) ! Limit increment dh>=min_thickness
        cellHeight(k)    = iFaceHeight(k) - 0.5 * dh
        iFaceHeight(k+1) = iFaceHeight(k) - dh
      enddo

      ! gets index of the level and interface above hbl
      kOBL = CVMix_kpp_compute_kOBL_depth(iFaceHeight, cellHeight,hbl(i,j))

      kv_col(:) = 0.0 ; kd_col(:) = 0.0
      call CVMix_coeffs_conv(Mdiff_out=kv_col(:), &
                               Tdiff_out=kd_col(:), &
                               Nsqr=CS%N2(i,j,:), &
                               dens=rho_1d(:), &
                               dens_lwr=rho_lwr(:), &
                               nlev=G%ke,    &
                               max_nlev=G%ke, &
                               OBL_ind=kOBL)

      do K=1,G%ke+1
        CS%kv_conv(i,j,K) = US%m2_s_to_Z2_T * kv_col(K)
        CS%Kd_conv(i,j,K) = US%m2_s_to_Z2_T * kd_col(K)
      enddo
      ! Do not apply mixing due to convection within the boundary layer
      do k=1,kOBL
        CS%kv_conv(i,j,k) = 0.0
        CS%kd_conv(i,j,k) = 0.0
      enddo

    enddo
  enddo

  if (CS%debug) then
    call hchksum(CS%N2, "MOM_CVMix_conv: N2",G%HI,haloshift=0)
    call hchksum(CS%kd_conv, "MOM_CVMix_conv: kd_conv",G%HI,haloshift=0,scale=US%Z2_T_to_m2_s)
    call hchksum(CS%kv_conv, "MOM_CVMix_conv: kv_conv",G%HI,haloshift=0,scale=US%m2_s_to_Z2_T)
  endif

  ! send diagnostics to post_data
  if (CS%id_N2 > 0) call post_data(CS%id_N2, CS%N2, CS%diag)
  if (CS%id_kd_conv > 0) call post_data(CS%id_kd_conv, CS%kd_conv, CS%diag)
  if (CS%id_kv_conv > 0) call post_data(CS%id_kv_conv, CS%kv_conv, CS%diag)

end subroutine calculate_CVMix_conv

!> Reads the parameter "USE_CVMix_CONVECTION" and returns state.
!! This function allows other modules to know whether this parameterization will
!! be used without needing to duplicate the log entry.
logical function CVMix_conv_is_used(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  call get_param(param_file, mdl, "USE_CVMix_CONVECTION", CVMix_conv_is_used, &
                 default=.false., do_not_log = .true.)

end function CVMix_conv_is_used

!> Clear pointers and dealocate memory
subroutine CVMix_conv_end(CS)
  type(CVMix_conv_cs), pointer :: CS !< Control structure for this module that
                                     !! will be deallocated in this subroutine

  if (.not. associated(CS)) return

  deallocate(CS%N2)
  deallocate(CS%kd_conv)
  deallocate(CS%kv_conv)
  deallocate(CS)

end subroutine CVMix_conv_end

end module MOM_CVMix_conv
