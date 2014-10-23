!> Implements the Mesoscale Eddy Kinetic Energy framework
module MOM_MEKE

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_domains,       only : create_group_pass, do_group_pass
use MOM_domains,       only : group_pass_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE, MOM_mesg
use MOM_file_parser, only : read_param, get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_io, only : vardesc
use MOM_restart, only : MOM_restart_CS, register_restart_field
use MOM_variables, only : vertvisc_type
use MOM_MEKE_types, only : MEKE_type

implicit none ; private

#include <MOM_memory.h>

public step_forward_MEKE, MEKE_init, MEKE_alloc_register_restart, MEKE_end

!> Control structure that contains MEKE parameters and diagnostics handles
type, public :: MEKE_CS ; private
  ! Parameters
  real :: MEKE_FrCoeff  !< Efficiency of conversion of ME into MEKE (non-dim)
  real :: MEKE_GMcoeff  !< Efficiency of conversion of PE into MEKE (non-dim)
  real :: MEKE_damping  !< Local depth-independent MEKE dissipation rate in s-1.
  real :: MEKE_Cd_scale !< The ratio of the bottom eddy velocity to the column mean
                        !! eddy velocity, i.e. sqrt(2*MEKE). This should be less than 1
                        !! to account for the surface intensification of MEKE.
  logical :: visc_drag  !< If true use the vertvisc_type to calculate bottom drag.
  logical :: Rd_as_max_scale !< If true the length scale can not exceed the
                        !! first baroclinic deformation radius.
  real :: cdrag         !< The bottom drag coefficient for MEKE (non-dim).
  real :: MEKE_BGsrc    !< Background energy source for MEKE in W/kg (= m2 s-3).
  real :: MEKE_dtScale  !< Scale factor to accelerate time-stepping (non-dim.)
  real :: MEKE_KhCoeff  !< Scaling factor to convert MEKE into Kh (non-dim.)
  real :: MEKE_Uscale   !< MEKE velocity scale for bottom drag (m/s)
  real :: MEKE_KH       !< Background lateral diffusion of MEKE (m^2/s)
  real :: KhMEKE_Fac    !< A factor relating MEKE%Kh to the diffusivity used for
                        !! MEKE itself (nondimensional).
  real :: viscosity_coeff !< The scaling coefficient in the expression for
                        !! viscosity used to parameterize lateral momentum mixing
                        !! by unresolved eddies represented by MEKE.

  ! Diagnostic handles
  type(diag_ctrl), pointer :: diag !< A pointer to shared diagnostics data
  integer :: id_MEKE = -1, id_Ue = -1, id_Kh = -1, id_src = -1
  integer :: id_GM_src = -1, id_mom_src = -1, id_decay = -1
  integer :: id_KhMEKE_u = -1, id_KhMEKE_v = -1, id_Ku = -1

  ! Infrastructure
  integer :: id_clock_pass !< Clock for group pass calls
  type(group_pass_type) :: pass_MEKE, pass_Kh, pass_Ku !< Type for group-halo pass calls
end type MEKE_CS

contains

!> Integrates forward-in-time the MEKE eddy energy equation.
!! See \ref section_MEKE_equations
subroutine step_forward_MEKE(MEKE, h, visc, dt, G, CS)
  type(MEKE_type),                       pointer       :: MEKE !< MEKE data.
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h    !< Layer thickness (m or kg m-2).
  type(vertvisc_type),                   intent(in)    :: visc !< The vertical viscosity type.
  real,                                  intent(in)    :: dt   !< Model(baroclinic) time-step (s).
  type(ocean_grid_type),                 intent(inout) :: G    !< Ocean grid.
  type(MEKE_CS),                         pointer       :: CS   !< MEKE control structure.

! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    mass, &         ! The total mass of the water column, in kg m-2.
    I_mass, &       ! The inverse of mass, in m2 kg-1.
    src, &          ! The sum of all MEKE sources, in m2 s-3.
    MEKE_decay, &   ! The MEKE decay timescale, in s-1.
    MEKE_GM_src, &  ! The MEKE source from thickness mixing, in m2 s-3.
    MEKE_mom_src, & ! The MEKE source from momentum, in m2 s-3.
    drag_rate_visc, &
    drag_rate       ! The MEKE spindown timescale due to bottom drag, in s-1.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    MEKE_uflux, &   ! The zonal diffusive flux of MEKE, in kg m2 s-3.
    Kh_u, &         ! The zonal diffusivity that is actually used, in m2 s-1.
    drag_vel_u      ! A (vertical) viscosity associated with bottom drag at
                    ! u-points, in m s-1.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    MEKE_vflux, &   ! The meridional diffusive flux of MEKE, in kg m2 s-3.
    Kh_v, &         ! The meridional diffusivity that is actually used, in m2 s-1.
    drag_vel_v      ! A (vertical) viscosity associated with bottom drag at
                    ! v-points, in m s-1.
  real :: Kh_here, Inv_Kh_max
  real :: cdrag2, bottomFac2
  real :: mass_neglect ! A negligible mass, in kg m-2.
  real :: ldamping  ! The MEKE damping rate in s-1.
  real :: Rho0      ! A density used to convert mass to distance, in kg m-3.
  real :: sdt  ! dt to use locally (could be scaled to accelerate)
  real :: sdt_damp  ! dt for damping (sdt could be split).
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  sdt = dt*CS%MEKE_dtScale ! Scaled dt to use for time-stepping
  Rho0 = G%H_to_kg_m2 * G%m_to_H
  mass_neglect = G%H_to_kg_m2 * G%H_subroundoff

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "MOM_MEKE: Module must be initialized before it is used.")
  if (.not.associated(MEKE)) call MOM_error(FATAL, &
         "MOM_MEKE: MEKE must be initialized before it is used.")

  ! Only integerate the MEKE equations if MEKE is required.
  if (associated(MEKE%MEKE)) then

    if (CS%MEKE_Cd_scale == 0.0 .and. .not. CS%visc_drag) then
!$OMP do
      do j=js,je ; do i=is,ie
        drag_rate(i,j) = 0.
      enddo ; enddo
    endif

    sdt = dt*CS%MEKE_dtScale ! Scaled dt to use for time-stepping
    Rho0 = G%H_to_kg_m2 * G%m_to_H
    mass_neglect = G%H_to_kg_m2 * G%H_subroundoff
    cdrag2 = CS%cdrag**2
    bottomFac2 = CS%MEKE_CD_SCALE**2
    ! With a depth-dependent (and possibly strong) damping, it seems
    ! advisable to use Strang splitting between the damping and diffusion.
    sdt_damp = sdt ; if (CS%MEKE_KH >= 0.0) sdt_damp = 0.5*sdt

!$OMP parallel default(none) shared(MEKE,CS,is,ie,js,je,nz,src,mass,G,h,I_mass, &
!$OMP                               sdt,drag_vel_u,visc,drag_vel_v,drag_rate_visc, &
!$OMP                               drag_rate,Rho0,MEKE_decay,sdt_damp,cdrag2, &
!$OMP                               bottomFac2) &
!$OMP                       private(ldamping)

!$OMP do
    do j=js,je ; do i=is,ie
      src(i,j) = CS%MEKE_BGsrc
    enddo ; enddo 

!$OMP do
    do j=js-1,je+1 
      do i=is-1,ie+1 ; mass(i,j) = 0.0 ; enddo 
      do k=1,nz ; do i=is-1,ie+1
        mass(i,j) = mass(i,j) + G%H_to_kg_m2 * h(i,j,k)
      enddo ; enddo
      do i=is-1,ie+1
        I_mass(i,j) = 0.0
        if (mass(i,j) > 0.0) I_mass(i,j) = 1.0 / mass(i,j)
        if (mass(i,j) < 0.0) mass(i,j) = 0.0  ! Should be unneccesary?
      enddo 
    enddo

    if (associated(MEKE%mom_src)) then
   !!$OMP do
   !  do j=js,je ; do i=is,ie
   !    MEKE%MEKE(i,j) = MEKE%MEKE(i,j) - I_mass(i,j) * &
   !        (sdt*CS%MEKE_FrCoeff)*MEKE%mom_src(i,j)
   !  enddo ; enddo
!$OMP do
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_FrCoeff*I_mass(i,j)*MEKE%mom_src(i,j)
      enddo ; enddo
    endif

    if (associated(MEKE%GM_src)) then
   !!$OMP do
   !  do j=js,je ; do i=is,ie
   !    MEKE%MEKE(i,j) = MEKE%MEKE(i,j) - I_mass(i,j) * &
   !        (sdt*CS%MEKE_GMcoeff)*MEKE%GM_src(i,j)
   !  enddo ; enddo
!$OMP do
      do j=js,je ; do i=is,ie
        src(i,j) = src(i,j) - CS%MEKE_GMcoeff*I_mass(i,j)*MEKE%GM_src(i,j)
      enddo ; enddo
    endif

   !!$OMP do
   !do j=js,je ; do i=is,ie
   !  MEKE%MEKE(i,j) = max(0.0,MEKE%MEKE(i,j) + (sdt*CS%MEKE_BGsrc)*G%mask2dT(i,j))
   !enddo ; enddo
!$OMP do
    do j=js,je ; do i=is,ie
      MEKE%MEKE(i,j) = max(0.0, MEKE%MEKE(i,j) + sdt*src(i,j) )*G%mask2dT(i,j)
    enddo ; enddo

! Following 3 lines seem identical to those above -AJA #######################################################
    ! With a depth-dependent (and possibly strong) damping, it seems
    ! advisable to use Strang splitting between the damping and diffusion.
    sdt_damp = sdt ; if (CS%MEKE_KH >= 0.0) sdt_damp = 0.5*sdt

    ! Calculate a viscous drag rate
    ! - drag_rate_visc(i,j) accounts for the model bottom mean flow
    if (CS%visc_drag) then
!$OMP do
      do j=js,je ; do I=is-1,ie
        drag_vel_u(I,j) = 0.0
        if ((G%mask2dCu(I,j) > 0.0) .and. (visc%bbl_thick_u(I,j) > 0.0)) &
          drag_vel_u(I,j) = visc%kv_bbl_u(I,j) / visc%bbl_thick_u(I,j)
      enddo ; enddo
!$OMP do
      do J=js-1,je ; do i=is,ie
        drag_vel_v(i,J) = 0.0
        if ((G%mask2dCu(i,J) > 0.0) .and. (visc%bbl_thick_v(i,J) > 0.0)) &
          drag_vel_v(i,J) = visc%kv_bbl_v(i,J) / visc%bbl_thick_v(i,J)
      enddo ; enddo

!$OMP do
      do j=js,je ; do i=is,ie
        drag_rate_visc(i,j) = (0.25*G%IareaT(i,j) * &
                ((G%areaCu(I-1,j)*drag_vel_u(I-1,j) + &
                  G%areaCu(I,j)*drag_vel_u(I,j)) + &
                 (G%areaCv(i,J-1)*drag_vel_v(i,J-1) + &
                  G%areaCv(i,J)*drag_vel_v(i,J)) ) )
        drag_rate(i,j) = (Rho0 * I_mass(i,j)) * sqrt( &
                 drag_rate_visc(i,j)**2               &
               + cdrag2 * ( max(0.0, 2.0*bottomFac2*MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2 ) )
      enddo ; enddo
    elseif (CS%MEKE_Cd_scale >= 0.0) then
!$OMP do
      do j=js,je ; do i=is,ie
        drag_rate(i,j) = (Rho0 * I_mass(i,j)) * sqrt( &
                 cdrag2 * ( max(0.0, 2.0*bottomFac2*MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2 ) )
      enddo ; enddo
    endif

    if (CS%MEKE_damping + CS%MEKE_Cd_scale > 0.0) then
      ! First stage of Strang splitting
!$OMP do
      do j=js,je ; do i=is,ie
        ldamping = CS%MEKE_damping + drag_rate(i,j) * bottomFac2
        MEKE%MEKE(i,j) = MEKE%MEKE(i,j) / (1.0 + sdt_damp*ldamping)
        MEKE_decay(i,j) = ldamping*G%mask2dT(i,j)
      enddo ; enddo
    endif
!$OMP end parallel
    if (CS%MEKE_KH >= 0.0) then
      call cpu_clock_begin(CS%id_clock_pass)
      call do_group_pass(CS%pass_MEKE, G%Domain)
      call cpu_clock_end(CS%id_clock_pass)
      
      ! Lateral diffusion of MEKE
      Kh_here = CS%MEKE_Kh
!$OMP parallel default(none) shared(is,ie,js,je,MEKE,CS,sdt,G,Kh_u,MEKE_uflux, &
!$OMP                               mass,mass_neglect,Kh_v,MEKE_vflux,I_mass, &
!$OMP                               sdt_damp,drag_rate,Rho0,drag_rate_visc,   &
!$OMP                               cdrag2,bottomFac2,MEKE_decay) &
!$OMP                       private(Kh_here,Inv_Kh_max,ldamping)
!$OMP do
      do j=js,je ; do I=is-1,ie
        ! Limit Kh to avoid CFL violations.
        if (associated(MEKE%Kh)) &
          Kh_here = CS%MEKE_Kh + CS%KhMEKE_Fac*0.5*(MEKE%Kh(i,j)+MEKE%Kh(i+1,j))
        Inv_Kh_max = 2.0*sdt * ((G%dy_Cu(I,j)*G%IdxCu(I,j)) * &
                     max(G%IareaT(i,j),G%IareaT(i+1,j)))
        if (Kh_here*Inv_Kh_max > 0.25) Kh_here = 0.25 / Inv_Kh_max
        Kh_u(I,j) = Kh_here

        MEKE_uflux(I,j) = ((Kh_here * (G%dy_Cu(I,j)*G%IdxCu(I,j))) * &
            ((2.0*mass(i,j)*mass(i+1,j)) / ((mass(i,j)+mass(i+1,j)) + mass_neglect)) ) * &
            (MEKE%MEKE(i,j) - MEKE%MEKE(i+1,j))
      enddo ; enddo
!$OMP do
      do J=js-1,je ; do i=is,ie
        if (associated(MEKE%Kh)) &
          Kh_here = CS%MEKE_Kh + CS%KhMEKE_Fac*0.5*(MEKE%Kh(i,j)+MEKE%Kh(i,j+1))
        Inv_Kh_max = 2.0*sdt * ((G%dx_Cv(i,J)*G%IdyCv(i,J)) * &
                     max(G%IareaT(i,j),G%IareaT(i,j+1)))
        if (Kh_here*Inv_Kh_max > 0.25) Kh_here = 0.25 / Inv_Kh_max
        Kh_v(i,J) = Kh_here

        MEKE_vflux(i,J) = ((Kh_here * (G%dx_Cv(i,J)*G%IdyCv(i,J))) * &
            ((2.0*mass(i,j)*mass(i,j+1)) / ((mass(i,j)+mass(i,j+1)) + mass_neglect)) ) * &
            (MEKE%MEKE(i,j) - MEKE%MEKE(i,j+1))
      enddo ; enddo
!$OMP do
      do j=js,je ; do i=is,ie
        MEKE%MEKE(i,j) = MEKE%MEKE(i,j) + (sdt*(G%IareaT(i,j)*I_mass(i,j))) * &
            ((MEKE_uflux(I-1,j) - MEKE_uflux(I,j)) + &
             (MEKE_vflux(i,J-1) - MEKE_vflux(i,J)))
      enddo ; enddo

      if ((CS%MEKE_damping + CS%MEKE_Cd_scale > 0.0) .and. (sdt>sdt_damp)) then
        ! Recalculate the drag rate, since MEKE has changed.
        if (CS%visc_drag) then
!$OMP do
          do j=js,je ; do i=is,ie
            drag_rate(i,j) = (Rho0 * I_mass(i,j)) * sqrt( &
                     drag_rate_visc(i,j)**2               &
                   + cdrag2 * ( max(0.0, 2.0*bottomFac2*MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2 ) )
          enddo ; enddo
        elseif (CS%MEKE_Cd_scale >= 0.0) then
!$OMP do
          do j=js,je ; do i=is,ie
            drag_rate(i,j) = (Rho0 * I_mass(i,j)) * sqrt( &
                     cdrag2 * ( max(0.0, 2.0*bottomFac2*MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2 ) )
          enddo ; enddo
        endif
        ! Second stage of Strang splitting
!$OMP do
        do j=js,je ; do i=is,ie
          ldamping = CS%MEKE_damping + drag_rate(i,j) * bottomFac2
          MEKE%MEKE(i,j) = MEKE%MEKE(i,j) / (1.0 + sdt_damp*ldamping)
          MEKE_decay(i,j) = 0.5 * G%mask2dT(i,j) * (MEKE_decay(i,j) + ldamping)
        enddo ; enddo
      endif
!$OMP end parallel
    endif

    call cpu_clock_begin(CS%id_clock_pass)
    call do_group_pass(CS%pass_MEKE, G%Domain)
    call cpu_clock_end(CS%id_clock_pass)

    ! Calculate diffusivity for main model to use
    if (CS%MEKE_KhCoeff>0.) then
      if (CS%Rd_as_max_scale) then
!$OMP parallel do default(none) shared(is,ie,js,je,MEKE,CS,G)
        do j=js-1,je+1 ; do i=is-1,ie+1
          MEKE%Kh(i,j) = (CS%MEKE_KhCoeff*sqrt(2.*max(0.,MEKE%MEKE(i,j))*G%areaT(i,j))) * &
                         min(MEKE%Rd_dx_h(i,j), 1.0)
        enddo ; enddo
      else
!$OMP parallel do default(none) shared(is,ie,js,je,MEKE,CS,G)
        do j=js-1,je+1 ; do i=is-1,ie+1
          MEKE%Kh(i,j) = CS%MEKE_KhCoeff*sqrt(2.*max(0.,MEKE%MEKE(i,j))*G%areaT(i,j))
        enddo ; enddo
      endif
      call cpu_clock_begin(CS%id_clock_pass)
      call do_group_pass(CS%pass_Kh, G%Domain)
      call cpu_clock_end(CS%id_clock_pass)
    endif

    ! Calculate viscosity for the main model to use
    if (CS%viscosity_coeff/=0.) then
!aja: should make range js:jeq, is:ieq
      do j=js-1,je+1 ; do i=is-1,ie+1
        MEKE%Ku(i,j) = CS%viscosity_coeff*sqrt(2.*max(0.,MEKE%MEKE(i,j))*G%areaT(i,j))
      enddo ; enddo
      call cpu_clock_begin(CS%id_clock_pass)
      call do_group_pass(CS%pass_Ku, G%Domain)
      call cpu_clock_end(CS%id_clock_pass)
    endif

! Offer fields for averaging.
    if (CS%id_MEKE>0) call post_data(CS%id_MEKE, MEKE%MEKE, CS%diag)
    if (CS%id_Ue>0) call post_data(CS%id_Ue, sqrt(2.0*MEKE%MEKE), CS%diag)
    if (CS%id_Kh>0) call post_data(CS%id_Kh, MEKE%Kh, CS%diag)
    if (CS%id_Ku>0) call post_data(CS%id_Ku, MEKE%Ku, CS%diag)
    if (CS%id_KhMEKE_u>0) call post_data(CS%id_KhMEKE_u, Kh_u, CS%diag)
    if (CS%id_KhMEKE_v>0) call post_data(CS%id_KhMEKE_v, Kh_v, CS%diag)
    if (CS%id_src>0) call post_data(CS%id_src, src, CS%diag)
    if (CS%id_decay>0) call post_data(CS%id_decay, MEKE_decay, CS%diag)
    if (CS%id_GM_src>0) call post_data(CS%id_GM_src, MEKE%GM_src, CS%diag)
    if (CS%id_mom_src>0) call post_data(CS%id_mom_src, MEKE%mom_src, CS%diag)

! else ! if MEKE%MEKE
!   call MOM_error(FATAL, "MOM_MEKE: MEKE%MEKE is not associated!")
  endif

end subroutine step_forward_MEKE

!> Initializes the MOM_MEKE module and reads parameters.
!! Returns True if module is to be used, otherwise returns False.
logical function MEKE_init(Time, G, param_file, diag, CS, MEKE)
  type(time_type),         intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),   intent(inout) :: G          !< The ocean's grid structure.
  type(param_file_type),   intent(in)    :: param_file !< Parameter file parser structure.
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(MEKE_CS),           pointer       :: CS         !< MEKE control structure.
  type(MEKE_type),         pointer       :: MEKE       !< MEKE-related fields.
! Local variables
  integer :: is, ie, js, je, isd, ied, jsd, jed, nz
  logical :: laplacian
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_MEKE" ! This module's name.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! Determine whether this module will be used
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "USE_MEKE", MEKE_init, &
                 "If true, turns on the MEKE scheme which calculates\n"// &
                 "a sub-grid mesoscale eddy kinetic energy budget.", &
                 default=.false.)
  if (.not. MEKE_init) return

  if (.not. associated(MEKE)) then
    ! The MEKE structure should have been allocated in MEKE_alloc_register_restart()
    call MOM_error(WARNING, "MEKE_init called with NO associated "// &
                            "MEKE-type structure.")
    return
  endif
  if (associated(CS)) then
    call MOM_error(WARNING, &
      "MEKE_init called with an associated control structure.")
    return
  else ; allocate(CS) ; endif

  call MOM_mesg("MEKE_init: reading parameters ", 5)

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mod, "MEKE_DAMPING", CS%MEKE_damping, &
                 "The local depth-indepented MEKE dissipation rate.", &
                 units="s-1", default=0.0)
  call get_param(param_file, mod, "MEKE_CD_SCALE", CS%MEKE_Cd_scale, &
                 "The ratio of the bottom eddy velocity to the column mean\n"//&
                 "eddy velocity, i.e. sqrt(2*MEKE). This should be less than 1\n"//&
                 "to account for the surface intensification of MEKE.", &
                 units="nondim", default=0.3)
  call get_param(param_file, mod, "MEKE_GMCOEFF", CS%MEKE_GMcoeff, &
                 "The efficiency of the conversion of potential energy \n"//&
                 "into MEKE by the thickness mixing parameterization. \n"//&
                 "If MEKE_GMCOEFF is negative, this conversion is not \n"//&
                 "used or calculated.", units="nondim", default=-1.0)
  call get_param(param_file, mod, "MEKE_FRCOEFF", CS%MEKE_FrCoeff, &
                 "The efficiency of the conversion of mean energy into \n"//&
                 "MEKE.  If MEKE_FRCOEFF is negative, this conversion \n"//&
                 "is not used or calculated.", units="nondim", default=-1.0)
  call get_param(param_file, mod, "MEKE_BGSRC", CS%MEKE_BGsrc, &
                 "A background energy source for MEKE.", units="W kg-1", &
                 default=0.0)
  call get_param(param_file, mod, "MEKE_KH", CS%MEKE_Kh, &
                 "A background lateral diffusivity of MEKE, or a \n"//&
                 "Use a negative value to not apply lateral diffusion to MEKE.", &
                 units="m2 s-1", default=-1.0)
  call get_param(param_file, mod, "MEKE_DTSCALE", CS%MEKE_dtScale, &
                 "A scaling factor to accelerate the time evolution of MEKE.", &
                 units="nondim", default=1.0)
  call get_param(param_file, mod, "MEKE_KHCOEFF", CS%MEKE_KhCoeff, &
                 "A scaling factor which is combined with the square root \n"//&
                 "of MEKE times the grid-cell area to give MEKE%Kh, or a \n"//&
                 "negative value not to calculate MEKE%Kh. \n"//&
                 "This factor must be >0 for MEKE to contribute to the \n"//&
                 "thickness/tracer mixing in the rest of the model. \n", &
                 units="nondim", default=-1.0)
  call get_param(param_file, mod, "MEKE_USCALE", CS%MEKE_Uscale, &
                 "The background velocity that is combined with MEKE to \n"//&
                 "calculate the bottom drag.", units="m s-1", default=0.0)
  call get_param(param_file, mod, "MEKE_VISC_DRAG", CS%visc_drag, &
                 "If true, use the vertvisc_type to calculate the bottom \n"//&
                 "drag acting on MEKE.", default=.true.)
  call get_param(param_file, mod, "MEKE_KHTH_FAC", MEKE%KhTh_fac, &
                 "A factor that maps MEKE%Kh to KhTh.", units="nondim", &
                 default=1.0)
  call get_param(param_file, mod, "MEKE_KHTR_FAC", MEKE%KhTr_fac, &
                 "A factor that maps MEKE%Kh to KhTr.", units="nondim", &
                 default=1.0)
  call get_param(param_file, mod, "MEKE_KHMEKE_FAC", CS%KhMEKE_Fac, &
                 "A factor that maps MEKE%Kh to Kh for MEKE itself.", &
                 units="nondim", default=0.0)
  call get_param(param_file, mod, "MEKE_RD_MAX_SCALE", CS%Rd_as_max_scale, &
                 "If true, the maximum length scale used by MEKE is \n"//&
                 "the deformation radius.", units="nondim", default=.true.)
  call get_param(param_file, mod, "MEKE_VISCOSITY_COEFF", CS%viscosity_coeff, &
                 "If non-zero, is the scaling coefficient in the expression for\n"//&
                 "viscosity used to parameterize lateral momentum mixing by\n"//&
                 "unresolved eddies represented by MEKE. Can be negative to\n"//&
                 "represent backscatter from the unresolved eddies.", &
                 units="nondim", default=0.0)

  ! Nonlocal module parameters
  call get_param(param_file, mod, "CDRAG", CS%cdrag, &
                 "CDRAG is the drag coefficient relating the magnitude of \n"//&
                 "the velocity field to the bottom stress.", units="nondim", &
                 default=0.003)
  call get_param(param_file, mod, "LAPLACIAN", laplacian, default=.false., do_not_log=.true.)
  if (CS%viscosity_coeff/=0. .and. .not. laplacian) call MOM_error(FATAL, &
                 "LAPLACIAN must be true if MEKE_VISCOSITY_COEFF is true.")

! In the case of a restart, these fields need a halo update
  if (associated(MEKE%MEKE)) then
    call create_group_pass(CS%pass_MEKE, MEKE%MEKE, G%Domain)
    call do_group_pass(CS%pass_MEKE, G%Domain)
  endif
  if (associated(MEKE%Kh)) then
    call create_group_pass(CS%pass_Kh, MEKE%Kh, G%Domain)
    call do_group_pass(CS%pass_Kh, G%Domain)
  endif
  if (associated(MEKE%Ku)) then
    call create_group_pass(CS%pass_Ku, MEKE%Ku, G%Domain)
    call do_group_pass(CS%pass_Ku, G%Domain)
  endif

! Register fields for output from this module.
  CS%diag => diag
  CS%id_MEKE = register_diag_field('ocean_model', 'MEKE', diag%axesT1, Time, &
     'Mesoscale Eddy Kinetic Energy', 'meter2 second-2')
  if (.not. associated(MEKE%MEKE)) CS%id_MEKE = -1
  CS%id_Kh = register_diag_field('ocean_model', 'MEKE_KH', diag%axesT1, Time, &
     'MEKE derived diffusivity', 'meter2 second-1')
  if (.not. associated(MEKE%Kh)) CS%id_Kh = -1
  CS%id_Ku = register_diag_field('ocean_model', 'MEKE_KU', diag%axesT1, Time, &
     'MEKE derived lateral viscosity', 'meter2 second-1')
  if (.not. associated(MEKE%Ku)) CS%id_Ku = -1
  CS%id_Ue = register_diag_field('ocean_model', 'MEKE_Ue', diag%axesT1, Time, &
     'MEKE derived eddy-velocity scale', 'meter second-1')
  if (.not. associated(MEKE%MEKE)) CS%id_Ue = -1
  CS%id_src = register_diag_field('ocean_model', 'MEKE_src', diag%axesT1, Time, &
     'MEKE energy source', 'meter2 second-3')
  CS%id_decay = register_diag_field('ocean_model', 'MEKE_decay', diag%axesT1, Time, &
     'MEKE decay rate', 'second-1')
  CS%id_KhMEKE_u = register_diag_field('ocean_model', 'KHMEKE_u', diag%axesCu1, Time, &
     'Zonal diffusivity of MEKE', 'meter2 second-1')
  CS%id_KhMEKE_v = register_diag_field('ocean_model', 'KHMEKE_v', diag%axesCv1, Time, &
     'Meridional diffusivity of MEKE', 'meter2 second-1')
  CS%id_GM_src = register_diag_field('ocean_model', 'MEKE_GM_src', diag%axesT1, Time, &
     'MEKE energy available from thickness mixing', 'Watt meter-2')
  if (.not. associated(MEKE%GM_src)) CS%id_GM_src = -1
  CS%id_mom_src = register_diag_field('ocean_model', 'MEKE_mom_src',diag%axesT1, Time, &
     'MEKE energy available from momentum', 'Watt meter-2')
  if (.not. associated(MEKE%mom_src)) CS%id_mom_src = -1

  CS%id_clock_pass = cpu_clock_id('(Ocean continuity halo updates)', grain=CLOCK_ROUTINE)

end function MEKE_init

!> Allocates memory and register restart fields for the MOM_MEKE module.
subroutine MEKE_alloc_register_restart(G, param_file, MEKE, restart_CS)
! Arguments
  type(ocean_grid_type), intent(inout) :: G          !< Grid structure
  type(param_file_type), intent(in)    :: param_file !< Parameter file parser structure.
  type(MEKE_type),       pointer       :: MEKE       !< A structure with MEKE-related fields.
  type(MOM_restart_CS),  pointer       :: restart_CS !< Restart control structure for MOM_MEKE.
! Local variables
  type(vardesc) :: vd
  real :: MEKE_GMcoeff, MEKE_FrCoeff, MEKE_KHCoeff, MEKE_viscCoeff
  logical :: useMEKE
  integer :: isd, ied, jsd, jed

! Determine whether this module will be used
  useMEKE = .false.; call read_param(param_file,"USE_MEKE",useMEKE)

! Read these parameters to determine what should be in the restarts
  MEKE_GMcoeff =-1.; call read_param(param_file,"MEKE_GMCOEFF",MEKE_GMcoeff)
  MEKE_FrCoeff =-1.; call read_param(param_file,"MEKE_FRCOEFF",MEKE_FrCoeff)
  MEKE_KhCoeff =-1.; call read_param(param_file,"MEKE_KHCOEFF",MEKE_KhCoeff)
  MEKE_viscCoeff =0.; call read_param(param_file,"MEKE_VISCOSITY_COEFF",MEKE_viscCoeff)

! Allocate control structure
  if (associated(MEKE)) then
    call MOM_error(WARNING, "MEKE_alloc_register_restart called with an associated "// &
                             "MEKE type.")
    return
  else; allocate(MEKE); endif

  if (.not. useMEKE) return
 
! Allocate memory
  call MOM_mesg("MEKE_alloc_register_restart: allocating and registering", 5)
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  allocate(MEKE%MEKE(isd:ied,jsd:jed)) ; MEKE%MEKE(:,:) = 0.0
  vd = vardesc("MEKE","Mesoscale Eddy Kinetic Energy",'h','1','s',"m2 s-2")
  call register_restart_field(MEKE%MEKE, vd, .false., restart_CS)
  if (MEKE_GMcoeff>=0.) then
    allocate(MEKE%GM_src(isd:ied,jsd:jed)) ; MEKE%GM_src(:,:) = 0.0
  endif
  if (MEKE_FrCoeff>=0.) then
    allocate(MEKE%mom_src(isd:ied,jsd:jed)) ; MEKE%mom_src(:,:) = 0.0
  endif
  if (MEKE_KhCoeff>=0.) then
    allocate(MEKE%Kh(isd:ied,jsd:jed)) ; MEKE%Kh(:,:) = 0.0
    vd = vardesc("MEKE_Kh","Lateral diffusivity from Mesoscale Eddy Kinetic Energy",'h','1','s',"m2 s-1")
    call register_restart_field(MEKE%Kh, vd, .false., restart_CS)
  endif
  allocate(MEKE%Rd_dx_h(isd:ied,jsd:jed)) ; MEKE%Rd_dx_h(:,:) = 0.0
  if (MEKE_viscCoeff/=0.) then
    allocate(MEKE%Ku(isd:ied,jsd:jed)) ; MEKE%Ku(:,:) = 0.0
    vd = vardesc("MEKE_Ah","Lateral viscosity from Mesoscale Eddy Kinetic Energy",'h','1','s',"m2 s-1")
    call register_restart_field(MEKE%Ku, vd, .false., restart_CS)
  endif

end subroutine MEKE_alloc_register_restart

!> Deallocates any variables allocated in MEKE_init or
!! MEKE_alloc_register_restart.
subroutine MEKE_end(MEKE, CS)
  type(MEKE_type), pointer :: MEKE !< A structure with MEKE-related fields.
  type(MEKE_CS),   pointer :: CS   !< The control structure for MOM_MEKE.

  if (associated(CS)) deallocate(CS)

  if (.not.associated(MEKE)) return

  if (associated(MEKE%MEKE)) deallocate(MEKE%MEKE)
  if (associated(MEKE%GM_src)) deallocate(MEKE%GM_src)
  if (associated(MEKE%mom_src)) deallocate(MEKE%mom_src)
  if (associated(MEKE%Kh)) deallocate(MEKE%Kh)
  if (associated(MEKE%Ku)) deallocate(MEKE%Ku)
  deallocate(MEKE)

end subroutine MEKE_end

!> \class mom_meke
!!
!! \section section_MEKE The Mesoscale Eddy Kinetic Energy (MEKE) framework
!!
!! The MEKE framework accounts for the mean potential energy removed by
!! the first order closures used to parameterize mesoscale eddies.
!! It requires closure at the second order, namely dissipation and transport
!! of eddy energy.
!!
!! Monitoring the sub-grid scale eddy energy budget provides a means to predict
!! a sub-grid eddy-velocity scale which can be used in the lower order closures.
!!
!! \subsection section_MEKE_equations MEKE equations
!!
!! The eddy kinetic energy equation is:
!! \f[ \partial_\tilde{t} E = 
!!   \overbrace{ \dot{E}_b + \gamma_\eta \dot{E}_\eta + \gamma_v \dot{E}_v
!!             }^\text{sources}
!! - \overbrace{ ( \lambda + C_d | U_d | \gamma_d^2 ) E 
!!             }^\text{local dissipation}
!! + \overbrace{ \nabla \cdot ( \kappa_E + \gamma_M \kappa_M ) \nabla E
!!             }^\text{smoothing}
!! \f]
!! where \f$ E \f$ is the eddy kinetic energy (variable <code>MEKE</code>) with units of
!! m<sup>2</sup>s<sup>-2</sup>,
!! and \f$\tilde{t} = a t\f$ is a scaled time. The non-dimensional factor
!! \f$ a\geq 1 \f$ is used to accelerate towards equilibrium.
!! 
!! The MEKE equation is two-dimensional and obtained by depth averaging the
!! the three-dimensional eddy energy equation. In the following expressions
!! \f$ \left< \phi \right> = \frac{1}{H} \int^\eta_{-D} \phi \, dz \f$ maps
!! three dimensional terms into the two-dimensional quantities needed.
!!
!! \subsubsection section_MEKE_source_terms MEKE source terms 
!!
!! The source term \f$ \dot{E}_b \f$ is a constant background source
!! of energy intended to avoid the limit \f$E\rightarrow 0\f$.
!!
!! The "GM" source term
!! \f[ \dot{E}_\eta = - \left< \overline{w^\prime b^\prime} \right> 
!! = \left< \kappa_h N^2S^2 \right>
!! \approx \left< \kappa_h g\prime |\nabla_\sigma \eta|^2 \right>\f]
!! equals the mean potential energy removed by the Gent-McWilliams closure,
!! and is excluded/included in the MEKE budget by the efficiency parameter
!! \f$ \gamma_\eta \in [0,1] \f$.
!!
!! The "frictional" source term
!! \f[ \dot{E}_{v} = \left<  u \cdot \tau_h \right> \f]
!! equals the mean kinetic energy removed by lateral viscous fluxes, and
!! is excluded/included in the MEKE budget by the efficiency parameter
!! \f$ \gamma_v \in [0,1] \f$.
!!
!! \subsubsection section_MEKE_dissipation_terms MEKE dissipation terms 
!!
!! The local dissipation of \f$ E \f$ is parameterized through a linear
!! damping, \f$\lambda\f$, and bottom drag, \f$ C_d | U_d | \gamma_d^2 \f$.
!! The \f$ \gamma_d \f$ accounts for the weak projection of the column-mean
!! eddy velocty to the bottom. In other words, the bottom velocity is
!! estimated as \f$ \gamma_d U_e \f$.
!! The bottom drag coefficient, \f$ C_d \f$ is the same as that used in the bottom
!! friction in the mean model equations.
!!
!! The bottom drag velocity scale, \f$ U_d \f$, has contributions from the
!! resolved state and \f$ E \f$:
!! \f[ U_d = \sqrt{ U_b^2 + |u|^2_{z=-D} + |\gamma_d U_e|^2 } .\f]
!! where the eddy velocity scale, \f$ U_e \f$, is given by:
!! \f[ U_e = \sqrt{ 2 E } .\f]
!! \f$ U_b \f$ is a constant background bottom velocity scale and is
!! typically not used (i.e. set to zero).
!!
!! \subsection section_MEKE_smoothing MEKE smoothing terms
!!
!! \f$ E \f$ is laterally diffused by a diffusivity \f$ \kappa_E + \gamma_M
!! \kappa_M \f$ where \f$ \kappa_E \f$ is a constant diffusivity and the term
!! \f$ \gamma_M \kappa_M \f$ is a "self diffusion" using a the diffusivity
!! calculated in the section \ref section_MEKE_diffucivity.
!!
!! \subsection section_MEKE_diffusivity Diffusivity derived from MEKE
!!
!! The predicted eddy velocity scale, \f$ U_e \f$, can be combined with a
!! mixing length scale to form a diffusivity:
!!
!! \f[  \kappa_M = \gamma_\kappa \sqrt{ U_e^2 A_\Delta } \f]
!! where \f$ A_\Delta \f$ is the area of the grid cell
!! and \f$ \gamma_\kappa \in [0,1] \f$ s a non-dimensional factor.
!!
!! \subsection section_MEKE_viscosity Viscosity derived from MEKE
!!
!! As for \f$ \kappa_M \f$, the predicted eddy velocity scale can be
!! used to form an eddy viscosity:
!!
!! \f[  \kappa_u = \gamma_u \sqrt{ U_e^2 A_\Delta } . \f]
!!
!! \subsection section_MEKE_limit_case Limit cases for local source-dissipative balance
!!
!! Note that in steady-state (or when \f$ a>>1 \f$) and there is no
!! diffusion of \f$ E \f$ then
!! \f[ \overline{E} \approx \frac{ \dot{E}_b + \gamma_\eta \dot{E}_\eta + 
!!               \gamma_v \dot{E}_v }{ \lambda + C_d|U_d|\gamma_d^2 } . \f]
!!
!! In the linear drag limit, where
!! \f$ U_e << \min(U_b, |u|_{z=-D}, C_d^{-1}\lambda) \f$, the equilibrium becomes 
!! \f$ \overline{E} \approx \frac{ \dot{E}_b + \gamma_\eta \dot{E}_\eta + 
!!               \gamma_v \dot{E}_v }{ \lambda + C_d \sqrt{ U_b^2 + |u|^2_{z=-D} } } \f$.
!!
!! In the nonlinear drag limit, where \f$ U_e >> \max(U_b, |u|_{z=-D}, C_d^{-1}\lambda) \f$,
!! the equilibrium becomes 
!! \f$ \overline{E} \approx ( \frac{ \dot{E}_b + \gamma_\eta \dot{E}_\eta +
!!               \gamma_v \dot{E}_v }{ \sqrt{2} C_d \gamma_d^3 } )^\frac{2}{3} \f$.
!!
!! \subsubsection section_MEKE_module_parameters MEKE module parameters
!!
!! | Symbol                | Module parameter |
!! | ------                | --------------- |
!! | -                     | <code>USE_MEKE</code> |
!! | \f$ a \f$             | <code>MEKE_DTSCALE</code> |
!! | \f$ \dot{E}_b \f$     | <code>MEKE_BGSRC</code> |
!! | \f$ \gamma_\eta \f$   | <code>MEKE_GMCOEFF</code> |
!! | \f$ \gamma_v \f$      | <code>MEKE_FrCOEFF</code> |
!! | \f$ \lambda \f$       | <code>MEKE_DAMPING</code> |
!! | \f$ U_b \f$           | <code>MEKE_USCALE</code> |
!! | \f$ \gamma_d \f$      | <code>MEKE_CD_SCALE</code> |
!! | \f$ \kappa_E \f$      | <code>MEKE_KH</code> |
!! | \f$ \gamma_\kappa \f$ | <code>MEKE_KHCOEFF</code> |
!! | \f$ \gamma_M \f$      | <code>MEKE_KHMEKE_FAC</code> |
!! | \f$ \gamma_u \f$      | <code>MEKE_VISCOSITY_COEFF</code> |
!! | -                     | <code>MEKE_KHTH_FAC</code> |
!! | -                     | <code>MEKE_KHTR_FAC</code> |
!!
!! | Symbol                | Model parameter |
!! | ------                | --------------- |
!! | \f$ C_d \f$           | <code>CDRAG</code> |

end module MOM_MEKE

