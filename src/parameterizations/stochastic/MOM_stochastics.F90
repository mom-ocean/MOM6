!> Top-level module for the MOM6 ocean model in coupled mode.
module MOM_stochastics

! This file is part of MOM6. See LICENSE.md for the license.

! This is the top level module for the MOM6 ocean model.  It contains routines
! for initialization, update, and writing restart of stochastic physics. This
! particular version wraps all of the calls for MOM6 in the calls that had
! been used for MOM4.
!
use MOM_debugging,           only : hchksum, uvchksum, qchksum
use MOM_diag_mediator,       only : register_diag_field, diag_ctrl, time_type, post_data
use MOM_diag_mediator,       only : register_static_field, enable_averages, disable_averaging
use MOM_grid,                only : ocean_grid_type
use MOM_variables,           only : thermo_var_ptrs
use MOM_domains,             only : pass_var, pass_vector, CORNER, SCALAR_PAIR
use MOM_verticalGrid,        only : verticalGrid_type
use MOM_error_handler,       only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_error_handler,       only : callTree_enter, callTree_leave
use MOM_file_parser,         only : get_param, log_version, close_param_file, param_file_type
use mpp_domains_mod,         only : domain2d, mpp_get_layout, mpp_get_global_domain
use mpp_domains_mod,         only : mpp_define_domains, mpp_get_compute_domain, mpp_get_data_domain
use MOM_domains,             only : root_PE, num_PEs
use MOM_coms,                only : Get_PElist
use MOM_EOS,                 only : calculate_density, EOS_domain
use stochastic_physics,      only : init_stochastic_physics_ocn, run_stochastic_physics_ocn

#include <MOM_memory.h>

implicit none ; private

public stochastics_init, update_stochastics, apply_skeb

!> This control structure holds parameters for the MOM_stochastics module
type, public:: stochastic_CS
  logical :: do_sppt         !< If true, stochastically perturb the diabatic
  logical :: do_skeb         !< If true, stochastically perturb the horizontal velocity
  logical :: skeb_use_gm     !< If true, adds GM work to the amplitude of SKEBS
  logical :: skeb_use_frict  !< If true, adds viscous dissipation rate to the amplitude of SKEBS
  logical :: pert_epbl       !< If true, then randomly perturb the KE dissipation and genration terms
  integer :: id_sppt_wts    = -1 !< Diagnostic id for SPPT
  integer :: id_skeb_wts    = -1 !< Diagnostic id for SKEB
  integer :: id_skebu       = -1 !< Diagnostic id for SKEB
  integer :: id_skebv       = -1 !< Diagnostic id for SKEB
  integer :: id_diss        = -1 !< Diagnostic id for SKEB
  integer :: skeb_npass     = -1 !< number of passes of the 9-point smoother for the dissipation estimate
  integer :: id_psi         = -1 !< Diagnostic id for SPPT
  integer :: id_epbl1_wts   = -1 !< Diagnostic id for epbl generation perturbation
  integer :: id_epbl2_wts   = -1 !< Diagnostic id for epbl dissipation perturbation
  integer :: id_skeb_taperu = -1 !< Diagnostic id for u taper of SKEB velocity increment
  integer :: id_skeb_taperv = -1 !< Diagnostic id for v taper of SKEB velocity increment
  real    :: skeb_gm_coef     !< If skeb_use_gm is true, then skeb_gm_coef * GM_work is added to the
                              !! dissipation rate used to set the amplitude of SKEBS [nondim]
  real    :: skeb_frict_coef  !< If skeb_use_frict is true, then skeb_gm_coef * GM_work is added to the
                              !! dissipation rate used to set the amplitude of SKEBS [nondim]
  real, allocatable :: skeb_diss(:,:,:) !< Dissipation rate used to set amplitude of SKEBS [L2 T-3 ~> m2 s-2]
                                        !! Index into this at h points.
  ! stochastic patterns
  real, allocatable :: sppt_wts(:,:)  !< Random pattern for ocean SPPT
                                      !! tendencies with a number between 0 and 2
  real, allocatable :: skeb_wts(:,:)  !< Random pattern for ocean SKEB
  real, allocatable :: epbl1_wts(:,:) !< Random pattern for K.E. generation
  real, allocatable :: epbl2_wts(:,:) !< Random pattern for K.E. dissipation
  type(time_type), pointer :: Time !< Pointer to model time (needed for sponges)
  type(diag_ctrl), pointer :: diag=>NULL() !< A structure that is used to regulate the

  ! Taper array to smoothly zero out the SKEBS velocity increment near land
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: taperCu !< Taper applied to u component of
                                                            !! stochastic velocity increment
                                                            !! range [0,1], [nondim]
  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: taperCv !< Taper applied to v component of
                                                            !! stochastic velocity increment
                                                            !! range [0,1], [nondim]

end type stochastic_CS

contains

!!   This subroutine initializes the stochastics physics control structure.
subroutine stochastics_init(dt, grid, GV, CS, param_file, diag, Time)
  real, intent(in)                       :: dt      !< time step [T ~> s]
  type(ocean_grid_type),   intent(in)    :: grid    !< horizontal grid information
  type(verticalGrid_type), intent(in)    :: GV      !< vertical grid structure
  type(stochastic_CS), pointer, intent(inout) :: CS !< stochastic control structure
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(diag_ctrl), target, intent(inout) :: diag    !< structure to regulate diagnostic output
  type(time_type), target                :: Time    !< model time

  ! Local variables
  integer, allocatable :: pelist(:) ! list of pes for this instance of the ocean
  integer :: mom_comm          ! list of pes for this instance of the ocean
  integer :: num_procs         ! number of processors to pass to stochastic physics
  integer :: iret              ! return code from stochastic physics
  integer :: pe_zero           !  root pe
  integer :: nxT, nxB          ! number of x-points including halo
  integer :: nyT, nyB          ! number of y-points including halo
  integer :: i, j, k           ! loop indices
  real    :: tmp(grid%isdB:grid%iedB,grid%jsdB:grid%jedB) ! Used to construct tapers
  integer :: taper_width       ! Width (in cells) of the taper that brings the stochastic velocity
                               ! increments to 0 at the boundary.

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "ocean_stochastics_init"  ! This module's name.

  call callTree_enter("stochastic_init(), MOM_stochastics.F90")
  if (associated(CS)) then
    call MOM_error(WARNING, "MOM_stochastics_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%Time => Time
  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  ! get number of processors and PE list for stochastic physics initialization
  call get_param(param_file, mdl, "DO_SPPT", CS%do_sppt, &
                 "If true, then stochastically perturb the thermodynamic "//&
                 "tendencies of T,S, amd h.  Amplitude and correlations are "//&
                 "controlled by the nam_stoch namelist in the UFS model only.", &
                 default=.false.)
  call get_param(param_file, mdl, "DO_SKEB", CS%do_skeb, &
                 "If true, then stochastically perturb the currents "//&
                 "using the stochastic kinetic energy backscatter scheme.",&
                 default=.false.)
  call get_param(param_file, mdl, "SKEB_NPASS", CS%skeb_npass, &
                 "number of passes of a 9-point smoother of the "//&
                 "dissipation estimate.", default=3, do_not_log=.not.CS%do_skeb)
  call get_param(param_file, mdl, "SKEB_TAPER_WIDTH", taper_width, &
                 "number of cells over which the stochastic velocity increment "//&
                 "is tapered to zero.", default=4, do_not_log=.not.CS%do_skeb)
  call get_param(param_file, mdl, "SKEB_USE_GM", CS%skeb_use_gm, &
                 "If true, adds GM work rate to the SKEBS amplitude.", &
                 default=.false., do_not_log=.not.CS%do_skeb)
  if ((.not. CS%do_skeb) .and. (CS%skeb_use_gm)) call MOM_error(FATAL, "If SKEB_USE_GM is True "//&
                 "then DO_SKEB must also be True.")
  call get_param(param_file, mdl, "SKEB_GM_COEF", CS%skeb_gm_coef, &
               "Fraction of GM work that is added to backscatter rate.", &
               units="nondim", default=0.0, do_not_log=.not.CS%skeb_use_gm)
  call get_param(param_file, mdl, "SKEB_USE_FRICT", CS%skeb_use_frict, &
                 "If true, adds horizontal friction dissipation rate "//&
                 "to the SKEBS amplitude.", default=.false., do_not_log=.not.CS%do_skeb)
  if ((.not. CS%do_skeb) .and. (CS%skeb_use_frict)) call MOM_error(FATAL, "If SKEB_USE_FRICT is "//&
                 "True then DO_SKEB must also be True.")
  call get_param(param_file, mdl, "SKEB_FRICT_COEF", CS%skeb_frict_coef, &
               "Fraction of horizontal friction work that is added to backscatter rate.", &
               units="nondim", default=0.0, do_not_log=.not.CS%skeb_use_frict)
  call get_param(param_file, mdl, "PERT_EPBL", CS%pert_epbl, &
                 "If true, then stochastically perturb the kinetic energy "//&
                 "production and dissipation terms.  Amplitude and correlations are "//&
                 "controlled by the nam_stoch namelist in the UFS model only.", &
                 default=.false.)

  if (CS%do_sppt .OR. CS%pert_epbl .OR. CS%do_skeb) then
    num_procs = num_PEs()
    allocate(pelist(num_procs))
    call Get_PElist(pelist,commID = mom_comm)
    pe_zero = root_PE()
    nxT = grid%ied - grid%isd + 1
    nyT = grid%jed - grid%jsd + 1
    nxB = grid%iedB - grid%isdB + 1
    nyB = grid%jedB - grid%jsdB + 1
    call init_stochastic_physics_ocn(dt, grid%geoLonT, grid%geoLatT, nxT, nyT, GV%ke, &
                                     grid%geoLonBu, grid%geoLatBu, nxB, nyB, &
                                     CS%pert_epbl, CS%do_sppt, CS%do_skeb, pe_zero, mom_comm, iret)
    if (iret/=0)  then
      call MOM_error(FATAL, "call to init_stochastic_physics_ocn failed")
      return
    endif

    if (CS%do_sppt) allocate(CS%sppt_wts(grid%isd:grid%ied,grid%jsd:grid%jed))
    if (CS%do_skeb) allocate(CS%skeb_wts(grid%isdB:grid%iedB,grid%jsdB:grid%jedB))
    if (CS%do_skeb) allocate(CS%skeb_diss(grid%isd:grid%ied,grid%jsd:grid%jed,GV%ke), source=0.)
    if (CS%pert_epbl) then
      allocate(CS%epbl1_wts(grid%isd:grid%ied,grid%jsd:grid%jed))
      allocate(CS%epbl2_wts(grid%isd:grid%ied,grid%jsd:grid%jed))
    endif
  endif

  CS%id_sppt_wts = register_diag_field('ocean_model', 'sppt_pattern', CS%diag%axesT1, Time, &
       'random pattern for sppt', 'None')
  CS%id_skeb_wts = register_diag_field('ocean_model', 'skeb_pattern', CS%diag%axesB1, Time, &
       'random pattern for skeb', 'None')
  CS%id_epbl1_wts = register_diag_field('ocean_model', 'epbl1_wts', CS%diag%axesT1, Time, &
      'random pattern for KE generation', 'None')
  CS%id_epbl2_wts = register_diag_field('ocean_model', 'epbl2_wts', CS%diag%axesT1, Time, &
      'random pattern for KE dissipation', 'None')
  CS%id_skebu = register_diag_field('ocean_model', 'skebu', CS%diag%axesCuL, Time, &
       'zonal current perts', 'None')
  CS%id_skebv = register_diag_field('ocean_model', 'skebv', CS%diag%axesCvL, Time, &
       'zonal current perts', 'None')
  CS%id_diss = register_diag_field('ocean_model', 'skeb_amp', CS%diag%axesTL, Time, &
       'SKEB amplitude', 'm s-1')
  CS%id_psi  = register_diag_field('ocean_model', 'psi', CS%diag%axesBL, Time, &
       'stream function', 'None')
  CS%id_skeb_taperu = register_static_field('ocean_model', 'skeb_taper_u', CS%diag%axesCu1, &
       'SKEB taper u', 'None', interp_method='none')
  CS%id_skeb_taperv = register_static_field('ocean_model', 'skeb_taper_v', CS%diag%axesCv1, &
       'SKEB taper v', 'None', interp_method='none')

  ! Initialize the "taper" fields. These fields multiply the components of the stochastic
  ! velocity increment in such a way as to smoothly taper them to zero at land boundaries.
  if ((CS%do_skeb) .or. (CS%id_skeb_taperu > 0) .or. (CS%id_skeb_taperv > 0)) then
    ALLOC_(CS%taperCu(grid%IsdB:grid%IedB,grid%jsd:grid%jed))
    ALLOC_(CS%taperCv(grid%isd:grid%ied,grid%JsdB:grid%JedB))
    ! Initialize taper from land mask
    do j=grid%jsd,grid%jed ; do I=grid%isdB,grid%iedB
      CS%taperCu(I,j) = grid%mask2dCu(I,j)
    enddo ; enddo
    do J=grid%jsdB,grid%jedB ; do i=grid%isd,grid%ied
      CS%taperCv(i,J) = grid%mask2dCv(i,J)
    enddo ; enddo
    ! Extend taper land
    do k=1,(taper_width / 2)
      do j=grid%jsc-1,grid%jec+1 ; do I=grid%iscB-1,grid%iecB+1
        tmp(I,j) = minval(CS%taperCu(I-1:I+1,j-1:j+1))
      enddo ; enddo
      do j=grid%jsc,grid%jec ; do I=grid%iscB,grid%iecB
        CS%taperCu(I,j) = minval(tmp(I-1:I+1,j-1:j+1))
      enddo ; enddo
      do J=grid%jscB-1,grid%jecB+1 ; do i=grid%isc-1,grid%iec+1
        tmp(i,J) = minval(CS%taperCv(i-1:i+1,J-1:J+1))
      enddo ; enddo
      do J=grid%jscB,grid%jecB ; do i=grid%isc,grid%iec
        CS%taperCv(i,J) = minval(tmp(i-1:i+1,J-1:J+1))
      enddo ; enddo
      ! Update halo
      call pass_vector(CS%taperCu, CS%taperCv, grid%Domain, SCALAR_PAIR)
    enddo
    ! Smooth tapers. Each call smooths twice.
    do k=1,(taper_width - (taper_width/2))
      call smooth_x9_uv(grid, CS%taperCu, CS%taperCv, zero_land=.true.)
      call pass_vector(CS%taperCu, CS%taperCv, grid%Domain, SCALAR_PAIR)
    enddo
  endif

  !call uvchksum("SKEB taper [uv]", CS%taperCu, CS%taperCv, grid%HI)

  if (CS%id_skeb_taperu > 0) call post_data(CS%id_skeb_taperu, CS%taperCu, CS%diag, .true.)
  if (CS%id_skeb_taperv > 0) call post_data(CS%id_skeb_taperv, CS%taperCv, CS%diag, .true.)

  if (CS%do_sppt .OR. CS%pert_epbl .OR. CS%do_skeb) &
    call MOM_mesg('            === COMPLETED MOM STOCHASTIC INITIALIZATION =====')

  call callTree_leave("stochastic_init(), MOM_stochastics.F90")

end subroutine stochastics_init

!> update_ocean_model uses the forcing in Ice_ocean_boundary to advance the
!! ocean model's state from the input value of Ocean_state (which must be for
!! time time_start_update) for a time interval of Ocean_coupling_time_step,
!! returning the publicly visible ocean surface properties in Ocean_sfc and
!! storing the new ocean properties in Ocean_state.
subroutine update_stochastics(CS)
  type(stochastic_CS),      intent(inout) :: CS        !< diabatic control structure
  call callTree_enter("update_stochastics(), MOM_stochastics.F90")

! update stochastic physics patterns before running next time-step
  call run_stochastic_physics_ocn(CS%sppt_wts,CS%skeb_wts,CS%epbl1_wts,CS%epbl2_wts)

  call callTree_leave("update_stochastics(), MOM_stochastics.F90")

end subroutine update_stochastics

subroutine apply_skeb(grid,GV,CS,uc,vc,thickness,tv,dt,Time_end)

  type(ocean_grid_type),   intent(in)    :: grid   !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV     !< ocean vertical grid
  type(stochastic_CS),     intent(inout) :: CS     !< stochastic control structure

  real, dimension(SZIB_(grid),SZJ_(grid),SZK_(GV)), intent(inout) :: uc        !< zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(grid),SZJB_(grid),SZK_(GV)), intent(inout) :: vc        !< meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(grid),SZJ_(grid),SZK_(GV)),  intent(in)    :: thickness !< thickness [H ~> m or kg m-2]
  type(thermo_var_ptrs),                            intent(in)    :: tv       !< points to thermodynamic fields
  real,                                       intent(in)    :: dt       !< time increment [T ~> s]
  type(time_type),                            intent(in)    :: Time_end !< Time at the end of the interval
! locals

  real, dimension(SZIB_(grid),SZJB_(grid),SZK_(GV)) :: psi         !< Streamfunction for stochastic velocity increments
                                                                   !! [L2 T-1 ~> m2 s-1]
  real, dimension(SZIB_(grid),SZJ_(grid) ,SZK_(GV)) :: ustar       !< Stochastic u velocity increment [L T-1 ~> m s-1]
  real, dimension(SZI_(grid) ,SZJB_(grid),SZK_(GV)) :: vstar       !< Stochastic v velocity increment [L T-1 ~> m s-1]
  real, dimension(SZI_(grid),SZJ_(grid))            :: diss_tmp    !< Temporary array used in smoothing skeb_diss
                                                                   !! [L2 T-3 ~> m2 s-2]
  real, dimension(3,3) :: local_weights                            !< 3x3 stencil weights used in smoothing skeb_diss
                                                                   !! [L2 ~> m2]

  real    :: shr,ten,tot,kh
  integer :: i,j,k,iter
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state

  call callTree_enter("apply_skeb(), MOM_stochastics.F90")

  if ((.not. CS%skeb_use_gm) .and. (.not. CS%skeb_use_frict)) then
    ! fill in halos with zeros
    do k=1,GV%ke
      do j=grid%jsd,grid%jed ; do i=grid%isd,grid%ied
        CS%skeb_diss(i,j,k) = 0.0
      enddo ; enddo
    enddo

    !kh needs to be scaled

    kh=1!(120*111)**2
    do k=1,GV%ke
      do j=grid%jsc,grid%jec ; do i=grid%isc,grid%iec
        ! Shear
        shr = (vc(i,J,k)-vc(i-1,J,k))*grid%mask2dCv(i,J)*grid%mask2dCv(i-1,J)*grid%IdxCv(i,J)+&
              (uc(I,j,k)-uc(I,j-1,k))*grid%mask2dCu(I,j)*grid%mask2dCu(I,j-1)*grid%IdyCu(I,j)
        ! Tension
        ten = (vc(i,J,k)-vc(i-1,J,k))*grid%mask2dCv(i,J)*grid%mask2dCv(i-1,J)*grid%IdyCv(i,J)+&
              (uc(I,j,k)-uc(I,j-1,k))*grid%mask2dCu(I,j)*grid%mask2dCu(I,j-1)*grid%IdxCu(I,j)

        tot = sqrt( shr**2 + ten**2 ) * grid%mask2dT(i,j)
        CS%skeb_diss(i,j,k) = tot**3 * kh * grid%areaT(i,j)!!**2
      enddo ; enddo
    enddo
  endif ! Sets CS%skeb_diss without GM or FrictWork

  ! smooth dissipation skeb_npass times
  do iter=1,CS%skeb_npass
    if (mod(iter,2) == 1) call pass_var(CS%skeb_diss, grid%domain)
    do k=1,GV%ke
      do j=grid%jsc-1,grid%jec+1 ; do i=grid%isc-1,grid%iec+1
        ! This does not preserve rotational symmetry
        local_weights = grid%mask2dT(i-1:i+1,j-1:j+1)*grid%areaT(i-1:i+1,j-1:j+1)
        diss_tmp(i,j) = sum(local_weights*CS%skeb_diss(i-1:i+1,j-1:j+1,k)) / &
                       (sum(local_weights) + 1.E-16)
      enddo ; enddo
      do j=grid%jsc-1,grid%jec+1 ; do i=grid%isc-1,grid%iec+1
        if (grid%mask2dT(i,j)==0.) cycle
        CS%skeb_diss(i,j,k) = diss_tmp(i,j)
      enddo ; enddo
    enddo
  enddo
  call pass_var(CS%skeb_diss, grid%domain)

  ! call hchksum(CS%skeb_diss, "SKEB DISS", grid%HI, haloshift=2)
  ! call qchksum(CS%skeb_wts, "SKEB WTS", grid%HI, haloshift=1)

  do k=1,GV%ke
    do J=grid%jscB-1,grid%jecB ; do I=grid%iscB-1,grid%iecB
      psi(I,J,k) = sqrt(0.25 * dt * max((CS%skeb_diss(i  ,j  ,k) + CS%skeb_diss(i+1,j+1,k)) + &
                                        (CS%skeb_diss(i  ,j+1,k) + CS%skeb_diss(i+1,j  ,k)), 0.) ) &
                                  * CS%skeb_wts(I,J)
    enddo ; enddo
  enddo
  !call qchksum(psi,"SKEB PSI", grid%HI, haloshift=1)
  !call pass_var(psi, grid%domain, position=CORNER)
  do k=1,GV%ke
    do j=grid%jsc,grid%jec ; do I=grid%iscB,grid%iecB
      ustar(I,j,k) = - (psi(I,J,k) - psi(I,J-1,k)) * CS%taperCu(I,j) * grid%IdyCu(I,j)
      uc(I,j,k) = uc(I,j,k) + ustar(I,j,k)
    enddo ; enddo
    do J=grid%jscB,grid%jecB ; do i=grid%isc,grid%iec
      vstar(i,J,k) =   (psi(I,J,k) - psi(I-1,J,k)) * CS%taperCv(i,J) * grid%IdxCv(i,J)
      vc(i,J,k) = vc(i,J,k) + vstar(i,J,k)
    enddo ; enddo
  enddo

  !call uvchksum("SKEB increment [uv]", ustar, vstar, grid%HI)

  call enable_averages(dt, Time_end, CS%diag)
  if (CS%id_diss > 0) then
     call post_data(CS%id_diss, sqrt(dt * max(CS%skeb_diss(:,:,:), 0.)), CS%diag)
  endif
  if (CS%id_skeb_wts > 0) then
     call post_data(CS%id_skeb_wts, CS%skeb_wts, CS%diag)
  endif
  if (CS%id_skebu > 0) then
     call post_data(CS%id_skebu, ustar(:,:,:), CS%diag)
  endif
  if (CS%id_skebv > 0) then
     call post_data(CS%id_skebv, vstar(:,:,:), CS%diag)
  endif
  if (CS%id_psi > 0) then
     call post_data(CS%id_psi, psi(:,:,:), CS%diag)
  endif
  call disable_averaging(CS%diag)
  CS%skeb_diss(:,:,:) = 0.0 ! Must zero before next time step.

  call callTree_leave("apply_skeb(), MOM_stochastics.F90")

end subroutine apply_skeb

!> Apply a 9-point smoothing filter twice to a pair of velocity components to reduce
!! horizontal two-grid-point noise.
!! Note that this subroutine does not conserve angular momentum, so don't use it
!! in situations where you need conservation.  Also note that it assumes that the
!! input fields have valid values in the first two halo points upon entry.
subroutine smooth_x9_uv(G, field_u, field_v, zero_land)
  type(ocean_grid_type),             intent(in)    :: G         !< Ocean grid
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: field_u   !< u-point field to be smoothed[arbitrary]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: field_v   !< v-point field to be smoothed [arbitrary]
  logical,                 optional, intent(in)    :: zero_land !< If present and false, return the average
                                                                !! of the surrounding ocean points when
                                                                !! smoothing, otherwise use a value of 0 for
                                                                !! land points and include them in the averages.

  ! Local variables.
  real :: fu_prev(SZIB_(G),SZJ_(G))  ! The value of the u-point field at the previous iteration [arbitrary]
  real :: fv_prev(SZI_(G),SZJB_(G))  ! The value of the v-point field at the previous iteration [arbitrary]
  real :: Iwts             ! The inverse of the sum of the weights [nondim]
  logical :: zero_land_val ! The value of the zero_land optional argument or .true. if it is absent.
  integer :: i, j, s, is, ie, js, je, Isq, Ieq, Jsq, Jeq

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  zero_land_val = .true. ; if (present(zero_land)) zero_land_val = zero_land

  do s=1,0,-1
    fu_prev(:,:) = field_u(:,:)
    ! apply smoothing on field_u using rotationally symmetric expressions.
    do j=js-s,je+s ; do I=Isq-s,Ieq+s ; if (G%mask2dCu(I,j) > 0.0) then
      Iwts = 0.0625
      if (.not. zero_land_val) &
        Iwts = 1.0 / ( (4.0*G%mask2dCu(I,j) + &
                        ( 2.0*((G%mask2dCu(I-1,j) + G%mask2dCu(I+1,j)) + &
                               (G%mask2dCu(I,j-1) + G%mask2dCu(I,j+1))) + &
                         ((G%mask2dCu(I-1,j-1) + G%mask2dCu(I+1,j+1)) + &
                          (G%mask2dCu(I-1,j+1) + G%mask2dCu(I+1,j-1))) ) ) + 1.0e-16 )
      field_u(I,j) = Iwts * ( 4.0*G%mask2dCu(I,j) * fu_prev(I,j) &
                            + (2.0*((G%mask2dCu(I-1,j) * fu_prev(I-1,j) + G%mask2dCu(I+1,j) * fu_prev(I+1,j)) + &
                                    (G%mask2dCu(I,j-1) * fu_prev(I,j-1) + G%mask2dCu(I,j+1) * fu_prev(I,j+1))) &
                              + ((G%mask2dCu(I-1,j-1) * fu_prev(I-1,j-1) + G%mask2dCu(I+1,j+1) * fu_prev(I+1,j+1)) + &
                                 (G%mask2dCu(I-1,j+1) * fu_prev(I-1,j+1) + G%mask2dCu(I+1,j-1) * fu_prev(I-1,j-1))) ))
    endif ; enddo ; enddo

    fv_prev(:,:) = field_v(:,:)
    ! apply smoothing on field_v using rotationally symmetric expressions.
    do J=Jsq-s,Jeq+s ; do i=is-s,ie+s ; if (G%mask2dCv(i,J) > 0.0) then
      Iwts = 0.0625
      if (.not. zero_land_val) &
        Iwts = 1.0 / ( (4.0*G%mask2dCv(i,J) + &
                        ( 2.0*((G%mask2dCv(i-1,J) + G%mask2dCv(i+1,J)) + &
                               (G%mask2dCv(i,J-1) + G%mask2dCv(i,J+1))) + &
                         ((G%mask2dCv(i-1,J-1) + G%mask2dCv(i+1,J+1)) + &
                          (G%mask2dCv(i-1,J+1) + G%mask2dCv(i+1,J-1))) ) ) + 1.0e-16 )
      field_v(i,J) = Iwts * ( 4.0*G%mask2dCv(i,J) * fv_prev(i,J) &
                            + (2.0*((G%mask2dCv(i-1,J) * fv_prev(i-1,J) + G%mask2dCv(i+1,J) * fv_prev(i+1,J)) + &
                                    (G%mask2dCv(i,J-1) * fv_prev(i,J-1) + G%mask2dCv(i,J+1) * fv_prev(i,J+1))) &
                              + ((G%mask2dCv(i-1,J-1) * fv_prev(i-1,J-1) + G%mask2dCv(i+1,J+1) * fv_prev(i+1,J+1)) + &
                                 (G%mask2dCv(i-1,J+1) * fv_prev(i-1,J+1) + G%mask2dCv(i+1,J-1) * fv_prev(i-1,J-1))) ))
    endif ; enddo ; enddo
  enddo

end subroutine smooth_x9_uv

end module MOM_stochastics

