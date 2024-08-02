!> Provides the ocean stochastic equation of state
module MOM_stoch_eos

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_diag_mediator,    only : register_diag_field, post_data, diag_ctrl
use MOM_error_handler,    only : MOM_error, FATAL
use MOM_file_parser,      only : get_param, param_file_type
use MOM_grid,             only : ocean_grid_type
use MOM_hor_index,        only : hor_index_type
use MOM_isopycnal_slopes, only : vert_fill_TS
use MOM_random,           only : PRNG, random_2d_constructor, random_2d_norm
use MOM_restart,          only : MOM_restart_CS, register_restart_field, is_new_run, query_initialized
use MOM_time_manager,     only : time_type
use MOM_unit_scaling,     only : unit_scale_type
use MOM_variables,        only : thermo_var_ptrs
use MOM_verticalGrid,     only : verticalGrid_type
!use random_numbers_mod,  only : getRandomNumbers, initializeRandomNumberStream, randomNumberStream

implicit none; private
#include <MOM_memory.h>

public MOM_stoch_eos_init
public MOM_stoch_eos_run
public stoch_EOS_register_restarts
public post_stoch_EOS_diags
public MOM_calc_varT

!> Describes parameters of the stochastic component of the EOS
!! correction, described in Stanley et al. JAMES 2020.
type, public :: MOM_stoch_eos_CS ; private
  real, allocatable :: l2_inv(:,:)  !< One over sum of the T cell side side lengths squared [L-2 ~> m-2]
  real, allocatable :: rgauss(:,:)  !< nondimensional random Gaussian [nondim]
  real        :: tfac=0.27          !< Nondimensional decorrelation time factor, ~1/3.7 [nondim]
  real        :: amplitude=0.624499 !< Nondimensional standard deviation of Gaussian [nondim]
  integer     :: seed               !< PRNG seed
  type(PRNG)  ::  rn_CS             !< PRNG control structure
  real, allocatable :: pattern(:,:) !< Random pattern for stochastic EOS [nondim]
  real, allocatable :: phi(:,:)     !< temporal correlation stochastic EOS [nondim]
  logical :: use_stoch_eos!< If true, use the stochastic equation of state (Stanley et al. 2020)
  real :: stanley_coeff   !< Coefficient correlating the temperature gradient
                          !! and SGS T variance [nondim]; if <0, turn off scheme in all codes
  real :: stanley_a       !< a in exp(aX) in stochastic coefficient [nondim]
  real :: kappa_smooth    !< A diffusivity for smoothing T/S in vanished layers [H Z T-1 ~> m2 s-1 or kg m-1 s-1]

  !>@{ Diagnostic IDs
  integer :: id_stoch_eos  = -1, id_stoch_phi  = -1, id_tvar_sgs = -1
  !>@}

end type MOM_stoch_eos_CS

contains

!> Initializes MOM_stoch_eos module, returning a logical indicating whether this module will be used.
logical function MOM_stoch_eos_init(Time, G, GV, US, param_file, diag, CS, restart_CS)
  type(time_type),         intent(in)    :: Time       !< Time for stochastic process
  type(ocean_grid_type),   intent(in)    :: G          !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV         !< Vertical grid structure
  type(unit_scale_type),   intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< structure indicating parameter file to parse
  type(diag_ctrl), target, intent(inout) :: diag       !< Structure used to control diagnostics
  type(MOM_stoch_eos_CS),  intent(inout) :: CS         !< Stochastic control structure
  type(MOM_restart_CS),    pointer       :: restart_CS !< A pointer to the restart control structure.

  ! local variables
  integer :: i,j

  MOM_stoch_eos_init = .false.

  CS%seed = 0

  call get_param(param_file, "MOM_stoch_eos", "STOCH_EOS", CS%use_stoch_eos, &
                 "If true, stochastic perturbations are applied "//&
                 "to the EOS in the PGF.", default=.false.)
  call get_param(param_file, "MOM_stoch_eos", "STANLEY_COEFF", CS%stanley_coeff, &
                 "Coefficient correlating the temperature gradient "//&
                 "and SGS T variance.", units="nondim", default=-1.0)
  call get_param(param_file, "MOM_stoch_eos", "STANLEY_A", CS%stanley_a, &
                 "Coefficient a which scales chi in stochastic perturbation of the "//&
                 "SGS T variance.", units="nondim", default=1.0, &
                 do_not_log=((CS%stanley_coeff<0.0) .or. .not.CS%use_stoch_eos))
  call get_param(param_file, "MOM_stoch_eos", "KD_SMOOTH", CS%kappa_smooth, &
                 "A diapycnal diffusivity that is used to interpolate "//&
                 "more sensible values of T & S into thin layers.", &
                 units="m2 s-1", default=1.0e-6, scale=GV%m2_s_to_HZ_T, &
                 do_not_log=(CS%stanley_coeff<0.0))

  ! Don't run anything if STANLEY_COEFF < 0
  if (CS%stanley_coeff >= 0.0) then
    if (.not.allocated(CS%pattern)) call MOM_error(FATAL, &
        "MOM_stoch_eos_CS%pattern is not allocated when it should be, suggesting that "//&
        "stoch_EOS_register_restarts() has not been called before MOM_stoch_eos_init().")

    allocate(CS%phi(G%isd:G%ied,G%jsd:G%jed), source=0.0)
    allocate(CS%l2_inv(G%isd:G%ied,G%jsd:G%jed), source=0.0)
    allocate(CS%rgauss(G%isd:G%ied,G%jsd:G%jed), source=0.0)
    call get_param(param_file, "MOM_stoch_eos", "SEED_STOCH_EOS", CS%seed, &
                 "Specfied seed for random number sequence ", default=0)
    call random_2d_constructor(CS%rn_CS, G%HI, Time, CS%seed)
    call random_2d_norm(CS%rn_CS, G%HI, CS%rgauss)
    ! fill array with approximation of grid area needed for decorrelation time-scale calculation
    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        CS%l2_inv(i,j) = 1.0 / ( (G%dxT(i,j)**2) + (G%dyT(i,j)**2) )
      enddo
    enddo

    if (.not.query_initialized(CS%pattern, "stoch_eos_pattern", restart_CS) .or. &
        is_new_run(restart_CS)) then
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        CS%pattern(i,j) = CS%amplitude*CS%rgauss(i,j)
      enddo ; enddo
    endif

    !register diagnostics
    CS%id_tvar_sgs = register_diag_field('ocean_model', 'tvar_sgs', diag%axesTL, Time, &
      'Parameterized SGS Temperature Variance ', 'None')
    if (CS%use_stoch_eos) then
      CS%id_stoch_eos = register_diag_field('ocean_model', 'stoch_eos', diag%axesT1, Time, &
        'random pattern for EOS', 'None')
      CS%id_stoch_phi = register_diag_field('ocean_model', 'stoch_phi', diag%axesT1, Time, &
        'phi for EOS', 'None')
    endif
  endif

  ! This module is only used if explicitly enabled or a positive correlation coefficient is set.
  MOM_stoch_eos_init = CS%use_stoch_eos .or. (CS%stanley_coeff >= 0.0)

end function MOM_stoch_eos_init

!> Register fields related to the stoch_EOS module for resarts
subroutine stoch_EOS_register_restarts(HI, param_file, CS, restart_CS)
  type(hor_index_type),    intent(in)    :: HI         !< Horizontal index structure
  type(param_file_type),   intent(in)    :: param_file !< structure indicating parameter file to parse
  type(MOM_stoch_eos_CS),  intent(inout) :: CS         !< Stochastic control structure
  type(MOM_restart_CS),    pointer       :: restart_CS !< A pointer to the restart control structure.

  call get_param(param_file, "MOM_stoch_eos", "STANLEY_COEFF", CS%stanley_coeff, &
                 "Coefficient correlating the temperature gradient "//&
                 "and SGS T variance.", units="nondim", default=-1.0, do_not_log=.true.)

  if (CS%stanley_coeff >= 0.0) then
    allocate(CS%pattern(HI%isd:HI%ied,HI%jsd:HI%jed), source=0.0)
    call register_restart_field(CS%pattern, "stoch_eos_pattern", .false., restart_CS, &
                                "Random pattern for stoch EOS", "nondim")
  endif

end subroutine stoch_EOS_register_restarts

!> Generates a pattern in space and time for the ocean stochastic equation of state
subroutine MOM_stoch_eos_run(G, u, v, delt, Time, CS)
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), &
                           intent(in)    :: u    !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), &
                           intent(in)    :: v    !< The meridional velocity [L T-1 ~> m s-1].
  real,                    intent(in)    :: delt !< Time step size for AR1 process [T ~> s].
  type(time_type),         intent(in)    :: Time !< Time for stochastic process
  type(MOM_stoch_eos_CS),  intent(inout) :: CS   !< Stochastic control structure

  ! local variables
  real    :: ubar, vbar ! Averaged velocities [L T-1 ~> m s-1]
  real    :: phi        ! A temporal correlation factor [nondim]
  integer :: i, j

  ! Return without doing anything if this capability is not enabled.
  if (.not.CS%use_stoch_eos) return

  call random_2d_constructor(CS%rn_CS, G%HI, Time, CS%seed)
  call random_2d_norm(CS%rn_CS, G%HI, CS%rgauss)

  ! advance AR(1)
  do j=G%jsc,G%jec
    do i=G%isc,G%iec
      ubar = 0.5*(u(I,j,1)*G%mask2dCu(I,j)+u(I-1,j,1)*G%mask2dCu(I-1,j))
      vbar = 0.5*(v(i,J,1)*G%mask2dCv(i,J)+v(i,J-1,1)*G%mask2dCv(i,J-1))
      phi = exp(-delt*CS%tfac * sqrt(((ubar**2) + (vbar**2))*CS%l2_inv(i,j)))
      CS%pattern(i,j) = phi*CS%pattern(i,j) + CS%amplitude*sqrt(1-phi**2)*CS%rgauss(i,j)
      CS%phi(i,j) = phi
    enddo
  enddo

end subroutine MOM_stoch_eos_run

!> Write out any diagnostics related to this module.
subroutine post_stoch_EOS_diags(CS, tv, diag)
  type(MOM_stoch_eos_CS), intent(in) :: CS  !< Stochastic control structure
  type(thermo_var_ptrs),  intent(in) :: tv  !< Thermodynamics structure
  type(diag_ctrl),        intent(inout) :: diag !< Structure to control diagnostics

  if (CS%id_stoch_eos > 0) call post_data(CS%id_stoch_eos, CS%pattern, diag)
  if (CS%id_stoch_phi > 0) call post_data(CS%id_stoch_phi, CS%phi, diag)
  if (CS%id_tvar_sgs > 0) call post_data(CS%id_tvar_sgs, tv%varT, diag)

end subroutine post_stoch_EOS_diags

!> Computes a parameterization of the SGS temperature variance
subroutine MOM_calc_varT(G, GV, US, h, tv, CS, dt)
  type(ocean_grid_type),   intent(in)   :: G   !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)   :: GV  !< Vertical grid structure
  type(unit_scale_type),   intent(in)   :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  &
                          intent(in)    :: h   !< Layer thickness [H ~> m]
  type(thermo_var_ptrs),  intent(inout) :: tv  !< Thermodynamics structure
  type(MOM_stoch_eos_CS), intent(inout) :: CS  !< Stochastic control structure
  real,                   intent(in)    :: dt  !< Time increment [T ~> s]

  ! local variables
  real, dimension(SZI_(G), SZJ_(G), SZK_(GV)) :: &
    T, &          !> The temperature (or density) [C ~> degC], with the values in
                  !! in massless layers filled vertically by diffusion.
    S             !> The filled salinity [S ~> ppt], with the values in
                  !! in massless layers filled vertically by diffusion.
  real :: hl(5)              !> Copy of local stencil of H [H ~> m]
  real :: dTdi2, dTdj2       !> Differences in T variance [C2 ~> degC2]
  integer :: i, j, k

  ! Nothing happens if a negative correlation coefficient is set.
  if (CS%stanley_coeff < 0.0) return

  ! This block does a thickness weighted variance calculation and helps control for
  ! extreme gradients along layers which are vanished against topography. It is
  ! still a poor approximation in the interior when coordinates are strongly tilted.
  if (.not. associated(tv%varT)) allocate(tv%varT(G%isd:G%ied, G%jsd:G%jed, GV%ke), source=0.0)
  call vert_fill_TS(h, tv%T, tv%S, CS%kappa_smooth*dt, T, S, G, GV, US, halo_here=1, larger_h_denom=.true.)

  do k=1,G%ke
    do j=G%jsc,G%jec
      do i=G%isc,G%iec
        hl(1) = h(i,j,k) * G%mask2dT(i,j)
        hl(2) = h(i-1,j,k) * G%mask2dCu(I-1,j)
        hl(3) = h(i+1,j,k) * G%mask2dCu(I,j)
        hl(4) = h(i,j-1,k) * G%mask2dCv(i,J-1)
        hl(5) = h(i,j+1,k) * G%mask2dCv(i,J)

        ! SGS variance in i-direction [C2 ~> degC2]
        dTdi2 = ( ( G%mask2dCu(I  ,j) * (G%IdxCu(I  ,j) * ( T(i+1,j,k) - T(i,j,k) )) &
                  + G%mask2dCu(I-1,j) * (G%IdxCu(I-1,j) * ( T(i,j,k) - T(i-1,j,k) )) &
                ) * G%dxT(i,j) * 0.5 )**2
        ! SGS variance in j-direction [C2 ~> degC2]
        dTdj2 = ( ( G%mask2dCv(i,J  ) * (G%IdyCv(i,J  ) * ( T(i,j+1,k) - T(i,j,k) )) &
                  + G%mask2dCv(i,J-1) * (G%IdyCv(i,J-1) * ( T(i,j,k) - T(i,j-1,k) )) &
                ) * G%dyT(i,j) * 0.5 )**2
        tv%varT(i,j,k) = CS%stanley_coeff * ( dTdi2 + dTdj2 )
        ! Turn off scheme near land
        tv%varT(i,j,k) = tv%varT(i,j,k) * (minval(hl) / (maxval(hl) + GV%H_subroundoff))
      enddo
    enddo
  enddo
  ! if stochastic, perturb
  if (CS%use_stoch_eos) then
    do k=1,G%ke
      do j=G%jsc,G%jec
        do i=G%isc,G%iec
          tv%varT(i,j,k) = exp(CS%stanley_a * CS%pattern(i,j)) * tv%varT(i,j,k)
        enddo
      enddo
    enddo
  endif
end subroutine MOM_calc_varT

end module MOM_stoch_eos
