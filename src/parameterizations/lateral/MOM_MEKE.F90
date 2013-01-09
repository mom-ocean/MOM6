module MOM_MEKE
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*    This program contains the subroutine that calculates the         *
!*  effects of horizontal viscosity, including parameterizations of    *
!*  the value of the viscosity itself. mesosclae_EKE calculates        *
!*  the evolution of sub-grid scale mesoscale EKE.                     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ptrs, time_type
use MOM_domains, only : pass_var !, pass_vector
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
! public MEKE_type

type, public :: MEKE_CS ; private
  ! Parameters
  real :: MEKE_FrCoeff ! Efficiency of conversion of ME into MEKE (non-dim)
  real :: MEKE_GMcoeff ! Efficiency of conversion of PE into MEKE (non-dim)
  real :: MEKE_damping ! Local depth-independent MEKE dissipation rate in s-1.
  real :: MEKE_Cd_scale ! Scaling for bottom drag applied to MEKE, less than
                       ! 1 to account for surface intensification of MEKE or
                       ! potential (rather than kinetic) energy (non-dim).
  logical :: visc_drag ! If true use the vertvisc_type to calculate bottom drag.
  logical :: Rd_as_max_scale ! If true the length scale can not exceed the
                       ! first baroclinic deformation radius.
  real :: cdrag        ! The bottom drag coefficient for MEKE (non-dim).
  real :: MEKE_BGsrc   ! Background energy source for MEKE in W/kg (= m2 s-3).
  real :: MEKE_dtScale ! Scale factor to accelerate time-stepping (non-dim.)
  real :: MEKE_KhCoeff ! Scaling factor to convert MEKE into Kh (non-dim.)
  real :: MEKE_Uscale  ! MEKE velocity scale for bottom drag (m/s)
  real :: MEKE_KH      ! Background lateral diffusion of MEKE (m^2/s)
  real :: KhMEKE_Fac   ! A factor relating MEKE%Kh to the diffusivity used for
                       ! MEKE itself (nondimensional).
  ! Diagnostics
  type(diag_ptrs), pointer :: diag ! A pointer to a structure of shareable
  integer :: id_MEKE = -1, id_Kh = -1, id_src = -1
  integer :: id_GM_src = -1, id_mom_src = -1, id_decay = -1
  integer :: id_KhMEKE_u = -1, id_KhMEKE_v = -1
end type MEKE_CS

integer :: id_clock_pass

contains

subroutine step_forward_MEKE(MEKE, h, visc, dt, G, CS)
  type(MEKE_type),                       pointer       :: MEKE 
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h
  type(vertvisc_type),                   intent(in)    :: visc
  real,                                  intent(in)    :: dt
  type(ocean_grid_type),                 intent(inout) :: G
  type(MEKE_CS),                         pointer       :: CS

! Arguments: MEKE - A structure with MEKE-related fields (intent in/out).
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      visc - A structure containing vertical viscosities and related
!                   fields.
!  (in)      dt - The time step in s.
!  (in)      G -  The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 MEKE_init.

!  By A. Adcroft and R. Hallberg, 2006-2009.
!    This subroutine determines updates the Mesoscale Eddy Kinitic Energy
! (MEKE).

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
  real :: cdrag2
  real :: mass_neglect ! A negligible mass, in kg m-2.
  real :: ldamping  ! The MEKE damping rate in s-1.
  real :: Rho0      ! A density used to convert mass to distance, in kg m-3.
  real :: sdt  ! dt to use locally (could be scaled to accelerate)
  real :: sdt_damp  ! dt for damping (sdt could be split).
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not.associated(CS)) call MOM_error(FATAL, &
         "MOM_MEKE: Module must be initialized before it is used.")
  if (.not.associated(MEKE)) call MOM_error(FATAL, &
         "MOM_MEKE: MEKE must be initialized before it is used.")

  ! Only integerate the MEKE equations if MEKE is required.
  if (associated(MEKE%MEKE)) then

    sdt = dt*CS%MEKE_dtScale ! Scaled dt to use for time-stepping
    if (CS%id_src > 0) then ; do j=js,je ; do i=is,ie
      src(i,j) = CS%MEKE_BGsrc
    enddo ; enddo ; endif

    Rho0 = G%H_to_kg_m2 * G%m_to_H
    mass_neglect = G%H_to_kg_m2 * G%H_subroundoff

    do j=js-1,je+1 ; do i=is-1,ie+1 ; mass(i,j) = 0.0 ; enddo ; enddo
    do k=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      mass(i,j) = mass(i,j) + G%H_to_kg_m2 * h(i,j,k)
    enddo ; enddo ; enddo
    do j=js-1,je+1 ; do i=is-1,ie+1
      I_mass(i,j) = 0.0
      if (mass(i,j) > 0.0) I_mass(i,j) = 1.0 / mass(i,j)
      if (mass(i,j) < 0.0) mass(i,j) = 0.0  ! Should be unneccesary?
    enddo ; enddo

    if (associated(MEKE%mom_src)) then
      do j=js,je ; do i=is,ie
        MEKE%MEKE(i,j) = MEKE%MEKE(i,j) - I_mass(i,j) * &
            (sdt*CS%MEKE_FrCoeff)*MEKE%mom_src(i,j)
      enddo ; enddo
      if (CS%id_src > 0) then
        do j=js,je ; do i=is,ie
          src(i,j) = src(i,j) - CS%MEKE_FrCoeff*I_mass(i,j)*MEKE%mom_src(i,j)
        enddo ; enddo
      endif
    endif

    if (associated(MEKE%GM_src)) then
      do j=js,je ; do i=is,ie
        MEKE%MEKE(i,j) = MEKE%MEKE(i,j) - I_mass(i,j) * &
            (sdt*CS%MEKE_GMcoeff)*MEKE%GM_src(i,j)
      enddo ; enddo
      if (CS%id_src > 0) then
        do j=js,je ; do i=is,ie
          src(i,j) = src(i,j) - CS%MEKE_GMcoeff*I_mass(i,j)*MEKE%GM_src(i,j)
        enddo ; enddo
      endif
    endif

    do j=js,je ; do i=is,ie
      MEKE%MEKE(i,j) = max(0.0,MEKE%MEKE(i,j) + (sdt*CS%MEKE_BGsrc)*G%mask2dT(i,j))
    enddo ; enddo

    ! With a depth-dependent (and possibly strong) damping, it seems
    ! advisable to use Strang-splitting between the damping and diffusion.
    sdt_damp = sdt ; if (CS%MEKE_KH >= 0.0) sdt_damp = 0.5*sdt

    if (CS%visc_drag) then
      do j=js,je ; do I=is-1,ie
        drag_vel_u(I,j) = 0.0
        if ((G%mask2dCu(I,j) > 0.0) .and. (visc%bbl_thick_u(I,j) > 0.0)) &
          drag_vel_u(I,j) = visc%kv_bbl_u(I,j) / visc%bbl_thick_u(I,j)
      enddo ; enddo
      do J=js-1,je ; do i=is,ie
        drag_vel_v(i,J) = 0.0
        if ((G%mask2dCu(i,J) > 0.0) .and. (visc%bbl_thick_v(i,J) > 0.0)) &
          drag_vel_v(i,J) = visc%kv_bbl_v(i,J) / visc%bbl_thick_v(i,J)
      enddo ; enddo

      cdrag2 = CS%MEKE_Cd_scale * CS%cdrag**2
      do j=js,je ; do i=is,ie
        drag_rate_visc(i,j) = (0.25*G%IareaT(i,j) * &
                ((G%areaCu(I-1,j)*drag_vel_u(I-1,j) + &
                  G%areaCu(I,j)*drag_vel_u(I,j)) + &
                 (G%areaCv(i,J-1)*drag_vel_v(i,J-1) + &
                  G%areaCv(i,J)*drag_vel_v(i,J)) ) )
        drag_rate(i,j) = (Rho0 * I_mass(i,j)) * &
            sqrt(drag_rate_visc(i,j)**2 + 2.0*cdrag2 * max(0.0,MEKE%MEKE(i,j)))
      enddo ; enddo
    elseif (CS%MEKE_Cd_scale >= 0.0) then
      do j=js,je ; do i=is,ie
        drag_rate(i,j) = (Rho0 * I_mass(i,j)) * &
            (CS%cdrag*(2.0*sqrt(max(0.0,MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2)))
      enddo ; enddo
    endif

    if (CS%MEKE_damping + CS%MEKE_Cd_scale > 0.0) then
      do j=js,je ; do i=is,ie
        ldamping = CS%MEKE_damping
        if (CS%MEKE_Cd_scale > 0.0) &
          ldamping = CS%MEKE_damping + CS%MEKE_Cd_scale * drag_rate(i,j)

        MEKE%MEKE(i,j) = MEKE%MEKE(i,j) / (1.0 + sdt_damp*ldamping)
        MEKE_decay(i,j) = ldamping*G%mask2dT(i,j)
      enddo ; enddo
    endif

    if (CS%MEKE_KH >= 0.0) then
      call cpu_clock_begin(id_clock_pass)
      call pass_var(MEKE%MEKE, G%Domain)
      call cpu_clock_end(id_clock_pass)
      
      ! ### More elaborate prescriptions for Kh could be used here.
      Kh_here = CS%MEKE_Kh
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
      do j=js,je ; do i=is,ie
        MEKE%MEKE(i,j) = MEKE%MEKE(i,j) + (sdt*(G%IareaT(i,j)*I_mass(i,j))) * &
            ((MEKE_uflux(I-1,j) - MEKE_uflux(I,j)) + &
             (MEKE_vflux(i,J-1) - MEKE_vflux(i,J)))
      enddo ; enddo

      if ((CS%MEKE_damping + CS%MEKE_Cd_scale > 0.0) .and. (sdt>sdt_damp)) then
        ! Recalculate the drag rate, since MEKE has changed.
        if (CS%visc_drag) then
          do j=js,je ; do i=is,ie
            drag_rate(i,j) = (Rho0 * I_mass(i,j)) * &
                sqrt(drag_rate_visc(i,j)**2 + cdrag2 * 2.0*max(0.0,MEKE%MEKE(i,j)))
          enddo ; enddo
        elseif (CS%MEKE_Cd_scale >= 0.0) then
          do j=js,je ; do i=is,ie
            drag_rate(i,j) = (Rho0 * I_mass(i,j)) * &
                (CS%cdrag*2.0*(sqrt(max(0.0,MEKE%MEKE(i,j)) + CS%MEKE_Uscale**2)))
          enddo ; enddo
        endif

        do j=js,je ; do i=is,ie
          ldamping = CS%MEKE_damping
          if (CS%MEKE_Cd_scale > 0.0) &
            ldamping = CS%MEKE_damping + CS%MEKE_Cd_scale * drag_rate(i,j)

          MEKE%MEKE(i,j) = MEKE%MEKE(i,j) / (1.0 + sdt_damp*ldamping)
          MEKE_decay(i,j) = 0.5 * G%mask2dT(i,j) * (MEKE_decay(i,j) + ldamping)
        enddo ; enddo
      endif
    endif

    call cpu_clock_begin(id_clock_pass)
    call pass_var(MEKE%MEKE, G%Domain)
    call cpu_clock_end(id_clock_pass)

    if (CS%MEKE_KhCoeff>0.) then ; if (CS%Rd_as_max_scale) then
      do j=js-1,je+1 ; do i=is-1,ie+1
        MEKE%Kh(i,j) = (CS%MEKE_KhCoeff*sqrt(2.*max(0.,MEKE%MEKE(i,j))*G%areaT(i,j))) * &
                       min(MEKE%Rd_dx_h(i,j), 1.0)
      enddo ; enddo
    else
      do j=js-1,je+1 ; do i=is-1,ie+1
        MEKE%Kh(i,j) = CS%MEKE_KhCoeff*sqrt(2.*max(0.,MEKE%MEKE(i,j))*G%areaT(i,j))
      enddo ; enddo
    endif ; endif
    call cpu_clock_begin(id_clock_pass)
    call pass_var(MEKE%Kh, G%Domain)
    call cpu_clock_end(id_clock_pass)

! Offer fields for averaging.
    if (CS%id_MEKE>0) call post_data(CS%id_MEKE, MEKE%MEKE, CS%diag)
    if (CS%id_Kh>0) call post_data(CS%id_Kh, MEKE%Kh, CS%diag)
    if (CS%id_KhMEKE_u>0) call post_data(CS%id_KhMEKE_u, Kh_u, CS%diag)
    if (CS%id_KhMEKE_v>0) call post_data(CS%id_KhMEKE_v, Kh_v, CS%diag)
    if (CS%id_src>0) call post_data(CS%id_src, src, CS%diag)
    if (CS%id_decay>0) call post_data(CS%id_decay, MEKE_decay, CS%diag)
    if (CS%id_GM_src>0) then
      call post_data(CS%id_GM_src, MEKE%GM_src, CS%diag)
    endif
    if (CS%id_mom_src>0) then
      call post_data(CS%id_mom_src, MEKE%mom_src, CS%diag)
    endif

! else ! if MEKE%MEKE
!   call MOM_error(FATAL, "MOM_MEKE: MEKE%MEKE is not associated!")
  endif

end subroutine step_forward_MEKE

! ------------------------------------------------------------------------------

subroutine MEKE_init(Time, G, param_file, diag, CS, MEKE)
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(inout) :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ptrs), target, intent(inout) :: diag
  type(MEKE_CS),           pointer       :: CS
  type(MEKE_type),         pointer       :: MEKE
!   This subroutine allocates space for and claculates the static variables used
! by this module.  The metrics may be effectively 0, 1, or 2-D arrays,
! while fields like the background viscosities are 2-D arrays.
! ALLOC is a macro defined in MOM_memory.h for allocate or nothing with
! static memory.
!
! Arguments:
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure containing pointers to common diagnostic fields.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in/out)  MEKE - A structure with MEKE-related fields (intent in/out due to
!                   halo updates).

! Local variables
  integer :: is, ie, js, je, isd, ied, jsd, jed, nz
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'
  character(len=40)  :: mod = "MOM_MEKE" ! This module's name.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not. associated(MEKE)) then
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

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, tagname, "")
  call get_param(param_file, mod, "MEKE_DAMPING", CS%MEKE_damping, &
                 "The local depth-indepented MEKE dissipation rate.", &
                 units="s-1", default=0.0)
  call get_param(param_file, mod, "MEKE_CD_SCALE", CS%MEKE_Cd_scale, &
                 "A scaling for the bottom drag applied to MEKE.  This \n"//&
                 "should be less than 1 to account for the surface \n"//&
                 "intensification of MEKE and the fraction of MEKE that \n"//&
                 "may be temporarily stored as potential energy.", &
                 units="nondim", default=0.0)
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
                 "calculate the bottom drag.", units="m s-1", default=1.0)
  call get_param(param_file, mod, "CDRAG", CS%cdrag, &
                 "CDRAG is the drag coefficient relating the magnitude of \n"//&
                 "the velocity field to the bottom stress.", units="nondim", &
                 default=0.003)
  call get_param(param_file, mod, "MEKE_VISC_DRAG", CS%visc_drag, &
                 "If true, use the vertvisc_type to calculate the bottom \n"//&
                 "drag acting on MEKE.", default=.false.)
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

! In the case of a restart, these fields need a halo update
  if (associated(MEKE%MEKE)) call pass_var(MEKE%MEKE, G%Domain)
  if (associated(MEKE%Kh)) call pass_var(MEKE%Kh, G%Domain)

! Register fields for output from this module.
  CS%id_MEKE = register_diag_field('ocean_model', 'MEKE', G%axesT1, Time, &
     'Mesoscale Eddy Kinetic Energy', 'meter2 second-2')
  CS%id_Kh = register_diag_field('ocean_model', 'MEKE_KH', G%axesT1, Time, &
     'MEKE derived diffusivity', 'meter2 second-1')
  CS%id_src = register_diag_field('ocean_model', 'MEKE_src', G%axesT1, Time, &
     'MEKE energy source', 'meter2 second-3')
  CS%id_decay = register_diag_field('ocean_model', 'MEKE_decay', G%axesT1, Time, &
     'MEKE decay rate', 'second-1')
  CS%id_KhMEKE_u = register_diag_field('ocean_model', 'KHMEKE_u', G%axesCu1, Time, &
     'Zonal diffusivity of MEKE', 'meter2 second-1')
  CS%id_KhMEKE_v = register_diag_field('ocean_model', 'KHMEKE_v', G%axesCv1, Time, &
     'Meridional diffusivity of MEKE', 'meter2 second-1')
  if (associated(MEKE%GM_src)) &
    CS%id_GM_src = register_diag_field('ocean_model', 'MEKE_GM_src', G%axesT1, &
        Time, 'MEKE energy available from thickness mixing', 'Watt meter-2')
  if (associated(MEKE%mom_src)) &
    CS%id_mom_src = register_diag_field('ocean_model', 'MEKE_mom_src',G%axesT1,&
        Time, 'MEKE energy available from momentum', 'Watt meter-2')

  id_clock_pass = cpu_clock_id('(Ocean continuity halo updates)', grain=CLOCK_ROUTINE)

end subroutine MEKE_init

subroutine MEKE_alloc_register_restart(G, param_file, MEKE, restart_CS)
! Arguments
  type(ocean_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  type(MEKE_type), pointer             :: MEKE
  type(MOM_restart_CS), pointer       :: restart_CS
! Local variables
  type(vardesc) :: vd
  real :: MEKE_damping, MEKE_GMcoeff, MEKE_FrCoeff, MEKE_KHCoeff
  integer :: isd, ied, jsd, jed

! Read these parameters to determine what should be in the restarts
  MEKE_damping =-1.; call read_param(param_file,"MEKE_DAMPING",MEKE_damping)
  MEKE_GMcoeff =-1.; call read_param(param_file,"MEKE_GMCOEFF",MEKE_GMcoeff)
  MEKE_FrCoeff =-1.; call read_param(param_file,"MEKE_FRCOEFF",MEKE_FrCoeff)
  MEKE_KhCoeff =-1.; call read_param(param_file,"MEKE_KHCOEFF",MEKE_KhCoeff)

! Allocate control structure
  if (associated(MEKE)) then
    call MOM_error(WARNING, "MEKE_alloc_register_restart called with an associated "// &
                             "MEKE type.")
    return
  else; allocate(MEKE); endif

! Allocate memory
  call MOM_mesg("MEKE_alloc_register_restart: allocating and registering", 5)
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  if (MEKE_damping>=0.) then
    allocate(MEKE%MEKE(isd:ied,jsd:jed)) ; MEKE%MEKE(:,:) = 0.0
    vd = vardesc("MEKE","Mesoscale Eddy Kinetic Energy",'h','1','s',"m2 s-2")
    call register_restart_field(MEKE%MEKE, MEKE%MEKE, vd, .false., restart_CS)
  endif
  if (MEKE_GMcoeff>=0.) then
    allocate(MEKE%GM_src(isd:ied,jsd:jed)) ; MEKE%GM_src(:,:) = 0.0
  endif
  if (MEKE_FrCoeff>=0.) then
    allocate(MEKE%mom_src(isd:ied,jsd:jed)) ; MEKE%mom_src(:,:) = 0.0
  endif
  if (MEKE_KhCoeff>=0.) then
    allocate(MEKE%Kh(isd:ied,jsd:jed)) ; MEKE%Kh(:,:) = 0.0
    vd = vardesc("MEKE_Kh","Lateral diffusivity from Mesoscale Eddy Kinetic Energy",'h','1','s',"m2 s-1")
    call register_restart_field(MEKE%Kh, MEKE%Kh, vd, .false., restart_CS)
  endif
  allocate(MEKE%Rd_dx_h(isd:ied,jsd:jed)) ; MEKE%Rd_dx_h(:,:) = 0.0

end subroutine MEKE_alloc_register_restart

subroutine MEKE_end(MEKE, CS)
! This subroutine deallocates any variables allocated in MEKE_init.
! Arguments: MEKE - A structure with MEKE-related fields (intent in/out).
!  (in/out)  CS - The control structure returned by a previous call to
!                 MEKE_init.
  type(MEKE_type), pointer :: MEKE
  type(MEKE_CS),   pointer :: CS

  if (associated(CS)) deallocate(CS)

  if (.not.associated(MEKE)) return

  if (associated(MEKE%MEKE)) deallocate(MEKE%MEKE)
  if (associated(MEKE%GM_src)) deallocate(MEKE%GM_src)
  if (associated(MEKE%mom_src)) deallocate(MEKE%mom_src)
  if (associated(MEKE%KH)) deallocate(MEKE%Kh)
  deallocate(MEKE)

end subroutine MEKE_end

end module MOM_MEKE
