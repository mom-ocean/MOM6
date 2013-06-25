module MOM_mixed_layer_restrat
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
!
!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, June 2002, 2006                                *
!*                                                                     *
!*    The subroutine in this file implements a parameterization of     *
!*  unresolved viscous mixed layer restratification of the mixed layer *
!*  as described in Fox-Kemper, Ferrari and Hallberg (JPO, 2008), and  *
!*  whose impacts are described in Fox-Kemper et al. (Ocean Modelling, *
!*  2011).  This is derived in part from the older parameterizaton     *
!*  that is described in Hallberg (Aha Hulikoa, 2003), which this new  *
!*  parameterization surpasses, which in turn is based on the          *
!*  subinertial mixed layer theory of Young (JPO, 1994).  There is no  *
!*  net horizontal volume transport due to this parameterization, and  *
!*  no direct effect below the mixed layer.                            *
!*                                                                     *
!*    This parameterization, developed in detail by Fox-Kemper, et al.,*
!*  sets the restratification  timescale to agree with his high-       *
!*  resolution studies of mixed layer restratification.  The run-time  *
!*  parameter FOX_KEMPER_ML_RESTRAT_COEF is a nondimensional number    *
!*  of order a few tens, proportional to the ratio of the deformation  *
!*  radius or the gridscale (whichever is smaller to the dominant      *
!*  horizontal lengthscale of the submesoscale mixed layer             *
!*  instabilities.                                                     *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, vh, vav                               *
!*    j    x ^ x ^ x   At >:  u, uh, uav                               *
!*    j    > o > o >   At o:  h                                        *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator, only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator, only : register_diag_field, safe_alloc_ptr, time_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_EOS, only : calculate_density

implicit none ; private

#include <MOM_memory.h>

public mixedlayer_restrat, mixedlayer_restrat_init

type, public :: mixedlayer_restrat_CS ; private
  real    :: ml_restrat_coef !   A nondimensional factor by which the 
                             ! instability is enhanced over what would be
                             ! predicted based on the resolved  gradients.  This
                             ! increases with grid spacing^2, up to something
                             ! of order 500.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer :: id_urestrat_time , id_vrestrat_time 
  integer :: id_uhml = -1, id_vhml = -1
end type mixedlayer_restrat_CS

contains

subroutine mixedlayer_restrat(h, uhtr, vhtr, tv, fluxes, dt, G, CS)
  real, dimension(NIMEM_,NJMEM_,NKMEM_),  intent(inout) :: h
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(inout) :: uhtr
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(inout) :: vhtr
  type(thermo_var_ptrs),                  intent(in)    :: tv
  type(forcing),                          intent(in)    :: fluxes
  real,                                   intent(in)    :: dt
  type(ocean_grid_type),                  intent(in)    :: G
  type(mixedlayer_restrat_CS),            pointer       :: CS
!    This subroutine does interface depth diffusion.  The fluxes are
!  limited to give positive definiteness, and the diffusivities are
!  limited to guarantee stability.

! Arguments: h - Layer thickness, in m or kg m-2.  (Intent in/out.)
!                The units of h are referred to as H below.
!  (in/out)  uhtr - Accumulated zonal mass fluxes in m3 or kg.
!  (in/out)  vhtr - Accumulated meridional mass fluxes in m3 or kg.
!  (in)      tv - A structure containing the thermobaric variables.
!  (in)      fluxes - A structure containing pointers to any possible
!                     forcing fields.  Unused fields have NULL ptrs.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 mixedlayer_restrat_init.

  real :: uhml(SZIB_(G),SZJ_(G),SZK_(G)) ! The zonal and meridional mixed layer
  real :: vhml(SZI_(G),SZJB_(G),SZK_(G)) ! fluxes, in m3 s-1 or kg s-1.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    h_avail       ! The volume available for diffusion out of each face of each
                  ! sublayer of the mixed layer, divided by dt, in units
                  ! of H m2 s-1 (i.e., m3 s-1 or kg s-1).
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &       ! The sum of the thicknesses of layers in the mixed layer, H.
    Rml_av        ! g_Rho0 times the average mixed layer density, in m s-2.
  real :: g_Rho0  ! G_Earth/Rho0 in m4 s-2 kg-1.
  real :: Rho0(SZI_(G)) ! Potential density relative to the surface, in kg m-3.
  real :: p0(SZI_(G))   ! A pressure of 0, in Pa.

  real :: h_vel         ! htot interpolated onto velocity points, in m (not H).
  real :: absf          ! The absolute value of f, interpolated to velocity
                        ! points, in s-1.
  real :: u_star        ! The surface friction velocity, interpolated to velocity
                        ! points, in m s-1.
  real :: mom_mixrate   ! The rate at which momentum is homogenized within the
                        ! mixed layer in s-1.
  real :: timescale     ! The mixing growth timescale in s.
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected, in H.
  real :: dz_neglect    ! A thickness in m that is so small it is usually lost
                        ! in roundoff and can be neglected, in m.
  real :: I4dt          ! 1 / 4 dt
  real :: I2htot        ! Twice the total mixed layer thickness,
  real :: z_topx2       ! depth of the top of a layer, and
  real :: hx2           ! layer thickness, all at velocity points and in H.
  real :: a(SZK_(G))    ! A nondimensional value relating the overall flux
                        ! magnitudes (uDml & vDml) to the realized flux in a
                        ! layer.  The vertical sum of a() through the pieces of
                        ! the mixed layer must be 0.
  real :: uDml(SZIB_(G))  ! The zonal and meridional volume fluxes in the upper
  real :: vDml(SZI_(G))   ! half of the mixed layer in H m2 s-1 (m3 s-1 or kg s-1).
  real :: utimescale_diag(SZIB_(G),SZJ_(G)) ! The restratification timescales
  real :: vtimescale_diag(SZI_(G),SZJB_(G)) ! in the zonal and meridional
                                            ! directions, in s, stored in 2-D
                                            ! arrays for diagnostic purposes.
  logical :: use_EOS    ! If true, density is calculated from T & S using an
                        ! equation of state.
  integer :: i, j, k, is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "Module must be initialized before it is used.")
  if ((G%nkml<2) .or. (CS%ml_restrat_coef<=0.0)) return

  uDml(:) = 0.0 ; vDml(:) = 0.0
  I4dt = 0.25 / dt
  g_Rho0 = G%g_Earth/G%Rho0
  use_EOS = associated(tv%eqn_of_state)
  h_neglect = G%H_subroundoff
  dz_neglect = G%H_subroundoff*G%H_to_m

  if (.not.use_EOS) call MOM_error(FATAL, "MOM_mixedlayer_restrat: "// &
         "An equation of state must be used with this module.")

  ! Fix this later for G%nkml >= 3.

  p0(:) = 0.0

  do j=js-1,je+1
    do i=is-1,ie+1
      htot(i,j) = 0.0 ; Rml_av(i,j) = 0.0
    enddo
    do k=1,G%nkml
      call calculate_density(tv%T(:,j,k),tv%S(:,j,k),p0,Rho0(:),is-1,ie-is+3,tv%eqn_of_state)
      do i=is-1,ie+1
        Rml_av(i,j) = Rml_av(i,j) + h(i,j,k)*Rho0(i)
        htot(i,j) = htot(i,j) + h(i,j,k)
        h_avail(i,j,k) = max(I4dt*G%areaT(i,j)*(h(i,j,k)-G%Angstrom),0.0)
      enddo
    enddo

    do i=is-1,ie+1
      Rml_av(i,j) = (g_Rho0*Rml_av(i,j)) / (htot(i,j) + h_neglect)
    enddo
  enddo

! TO DO:
!   1. Mixing extends below the mixing layer to the mixed layer.  Find it!
!   2. Add exponential tail to streamfunction?

  do j=js,je ; do i=is,ie ; utimescale_diag(i,j) = 0.0 ; enddo ; enddo
  do j=js,je ; do i=is,ie ; vtimescale_diag(i,j) = 0.0 ; enddo ; enddo

!   U - Component
  do j=js,je ; do I=is-1,ie
    h_vel = 0.5*(htot(i,j) + htot(i+1,j)) * G%H_to_m

    u_star = 0.5*(fluxes%ustar(i,j) + fluxes%ustar(i+1,j))
    absf = 0.5*(abs(G%CoriolisBu(I,J-1)) + abs(G%CoriolisBu(I,J)))
    ! peak ML visc: u_star * 0.41 * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2 
    ! 0.41 is the von Karmen constant, 9.8696 = pi^2.
    mom_mixrate = (0.41*9.8696)*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)

    timescale = timescale * CS%ml_restrat_coef
!        timescale = timescale*(2?)*(L_def/L_MLI)*min(EKE/MKE,1.0 + G%dyCv(i,j)**2/L_def**2))

    utimescale_diag(I,j) = timescale

    uDml(I) = timescale * G%mask2dCu(I,j)*G%dyCu(I,j)* &
        G%IdxCu(I,j)*(Rml_av(i+1,j)-Rml_av(i,j)) * (h_vel**2 * G%m_to_H)

    if (uDml(i) == 0) then
      do k=1,G%nkml ; uhml(I,j,k) = 0.0 ; enddo
    else
      I2htot = 1.0 / (htot(i,j) + htot(i+1,j) + h_neglect)
      z_topx2 = 0.0
      ! a(k) relates the sublayer transport to uDml with a linear profile.
      ! The sum of a through the mixed layers must be 0.
      do k=1,G%nkml
        hx2 = (h(i,j,k) + h(i+1,j,k) + h_neglect)
        a(k) = (hx2 * I2htot) * (2.0 - 4.0*(z_topx2+0.5*hx2)*I2htot)
        z_topx2 = z_topx2 + hx2
        if (a(k)*uDml(I) > 0.0) then
          if (a(k)*uDml(I) > h_avail(i,j,k)) uDml(I) = h_avail(i,j,k) / a(k)
        else
          if (-a(k)*uDml(I) > h_avail(i+1,j,k)) uDml(I) = -h_avail(i+1,j,k)/a(k)
        endif
      enddo
      do k=1,G%nkml
        uhml(I,j,k) = a(k)*uDml(I)
        uhtr(I,j,k) = uhtr(I,j,k) + uhml(I,j,k)*dt
      enddo
    endif
  enddo ; enddo

!  V- component
  do J=js-1,je ; do i=is,ie
    h_vel = 0.5*(htot(i,j) + htot(i,j+1)) * G%H_to_m

    u_star = 0.5*(fluxes%ustar(i,j) + fluxes%ustar(i,j+1))
    absf = 0.5*(abs(G%CoriolisBu(I-1,J)) + abs(G%CoriolisBu(I,J)))
    ! peak ML visc: u_star * 0.41 * (h_ml*u_star)/(absf*h_ml + 4.0*u_star)
    ! momentum mixing rate: pi^2*visc/h_ml^2 
    ! 0.41 is the von Karmen constant, 9.8696 = pi^2.
    mom_mixrate = (0.41*9.8696)*u_star**2 / &
                  (absf*h_vel**2 + 4.0*(h_vel+dz_neglect)*u_star)
    timescale = 0.0625 * (absf + 2.0*mom_mixrate) / (absf**2 + mom_mixrate**2)

    timescale = timescale * CS%ml_restrat_coef
!        timescale = timescale*(2?)*(L_def/L_MLI)*min(EKE/MKE,1.0 + G%dyCv(i,j)**2/L_def**2))

    vtimescale_diag(i,J) = timescale

    vDml(i) = timescale * G%mask2dCv(i,J)*G%dxCv(i,J)* &
        G%IdyCv(i,J)*(Rml_av(i,j+1)-Rml_av(i,j)) * (h_vel**2 * G%m_to_H)
    if (vDml(i) == 0) then
      do k=1,G%nkml ; vhml(i,J,k) = 0.0 ; enddo
    else
      I2htot = 1.0 / (htot(i,j) + htot(i,j+1) + h_neglect)
      z_topx2 = 0.0
      ! a relates the sublayer transport to uDml with a linear profile.
      ! The sum of a through the mixed layers must be 0.
      do k=1,G%nkml
        hx2 = (h(i,j,k) + h(i,j+1,k) + h_neglect)
        a(k) = (hx2 * I2htot) * (2.0 - 4.0*(z_topx2+0.5*hx2)*I2htot)
        z_topx2 = z_topx2 + hx2
        if (a(k)*vDml(i) > 0.0) then
          if (a(k)*vDml(i) > h_avail(i,j,k)) vDml(i) = h_avail(i,j,k) / a(k)
        else
          if (-a(k)*vDml(i) > h_avail(i,j+1,k)) vDml(i) = -h_avail(i,j+1,k)/a(k)
        endif
      enddo
      do k=1,G%nkml
        vhml(i,J,k) = a(k)*vDml(i)
        vhtr(i,J,k) = vhtr(i,J,k) + vhml(i,J,k)*dt
      enddo
    endif
  enddo ; enddo

  do k=1,G%nkml ; do j=js,je ; do i=is,ie
    h(i,j,k) = h(i,j,k) - dt*G%IareaT(i,j) * &
        ((uhml(I,j,k) - uhml(I-1,j,k)) + (vhml(i,J,k) - vhml(i,J-1,k)))
  enddo ; enddo ; enddo

! Offer fields for averaging.
  if (query_averaging_enabled(CS%diag) .and. &
      ((CS%id_urestrat_time > 0)  .or. (CS%id_vrestrat_time > 0))) then
    call post_data(CS%id_urestrat_time, utimescale_diag, CS%diag)
    call post_data(CS%id_vrestrat_time, vtimescale_diag, CS%diag)
  endif
  if (query_averaging_enabled(CS%diag) .and. &
      ((CS%id_uhml>0) .or. (CS%id_vhml>0))) then
    do k=G%nkml+1,nz
      do j=js,je ; do I=Isq,Ieq ; uhml(I,j,k) = 0.0 ; enddo ; enddo
      do J=Jsq,Jeq ; do i=is,ie ; vhml(i,J,k) = 0.0 ; enddo ; enddo
    enddo
    if (CS%id_uhml > 0) call post_data(CS%id_uhml, uhml, CS%diag)
    if (CS%id_vhml > 0) call post_data(CS%id_vhml, vhml, CS%diag)
  endif

end subroutine mixedlayer_restrat

subroutine mixedlayer_restrat_init(Time, G, param_file, diag, CS)
  type(time_type),             intent(in)    :: Time
  type(ocean_grid_type),       intent(in)    :: G
  type(param_file_type),       intent(in)    :: param_file
  type(diag_ctrl), target,     intent(inout) :: diag
  type(mixedlayer_restrat_CS), pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                  for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_mixed_layer_restrat"  ! This module's name.
  character(len=48)  :: flux_units

  if (associated(CS)) then
    call MOM_error(WARNING, "mixedlayer_restrat_init called with an "// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  if (G%Boussinesq) then ; flux_units = "meter3 second-1"
  else ; flux_units = "kilogram second-1" ; endif

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "FOX_KEMPER_ML_RESTRAT_COEF", &
                 CS%ml_restrat_coef, &
             "A nondimensional coefficient that is proportional to \n"//&
             "the ratio of the deformation radius to the dominant \n"//&
             "lengthscale of the submesoscale mixed layer \n"//&
             "instabilities, times the minimum of the ratio of the \n"//&
             "mesoscale eddy kinetic energy to the large-scale \n"//&
             "geostrophic kinetic energy or 1 plus the square of the \n"//&
             "grid spacing over the deformation radius, as detailed \n"//&
             "by Fox-Kemper et al. (2010)", units="nondim", default=0.0)

  CS%id_uhml = register_diag_field('ocean_model', 'uhml', G%axesCuL, Time, &
      'Zonal Thickness Flux to Restratify Mixed Layer', flux_units)
  CS%id_vhml = register_diag_field('ocean_model', 'vhml', G%axesCvL, Time, &
      'Meridional Thickness Flux to Restratify Mixed Layer', flux_units)
  CS%id_urestrat_time = register_diag_field('ocean_model', 'MLu_restrat_time', G%axesCu1, Time, &
      'Mixed Layer Zonal Restratification Timescale', 'second')
  CS%id_vrestrat_time = register_diag_field('ocean_model', 'MLv_restrat_time', G%axesCu1, Time, &
      'Mixed Layer Meridional Restratification Timescale', 'second')

end subroutine mixedlayer_restrat_init

end module MOM_mixed_layer_restrat
